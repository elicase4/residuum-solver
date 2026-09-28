# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Residuum is a modular finite element (FEM) PDE solver written in C++20. It is organized as a
layered set of static libraries (mesh, io, linalg, fem, solver, equation) consumed by small
CLI applications (`mesh`, `heateq`) that are driven entirely by YAML configuration files.

## Goals and scope

Residuum is intended to grow into a general-purpose FEM/IGA (isogeometric analysis) framework for
solving PDEs, with a **high-performance, backend-agnostic (CPU/CUDA) codebase as the primary
design driver** — that takes priority over breadth when the two are in tension.

Target physics, in order of ambition:

1. Heat equation (current focus — `equation/heateq`).
2. Stress-strain / elasticity.
3. Incompressible RANS — a deliberate stretch goal meant to demonstrate VMS (variational
   multiscale) stabilization and effective preconditioning, not just "another equation module."

Beyond single-physics solves, the architecture should support **multiphysics coupling** (e.g.
thermal stress: heat equation coupled with stress-strain). New abstractions in `fem`, `solver`,
and `equation` should be evaluated against "would this still make sense for elasticity/RANS, and
for two of these coupled together?", not just against the heat equation.

## Build

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug   # or Release (default)
make -j$(nproc)
```

Build options (CMake `option()`, default ON): `BUILD_TESTS`, `BUILD_APPLICATIONS`.

External dependencies (must be discoverable via `find_package`/`find_path`): `yaml-cpp`, `GTest`,
and `exprtk` (header-only, located via `find_path(... NAMES exprtk.hpp)`).

Binaries land in `build/bin` (apps) and `build/bin/tests` (test executables); libraries in `build/lib`.

`cmake --install . --prefix <dir>` installs the compiled libraries, `mesh`/`heateq` binaries, and
the full `include/` tree, plus a CMake package config so a downstream project can
`find_package(Residuum)` and link `Residuum::residuum_<module>` (see `cmake/ResiduumConfig.cmake.in`).
Every module publicly links `residuum_headers`, an `INTERFACE` library that owns the
`include/`-tree path via `BUILD_INTERFACE`/`INSTALL_INTERFACE` generator expressions — link against
it (transitively, via any `residuum_*` library, or directly) rather than adding `include/` as an
include directory by hand.

Run an application against one of the example configs. Relative paths inside a config (mesh
files, field data) resolve against the current directory, not the config's own location, so run
from inside the example's directory:

```bash
cd examples/mesh/square && /path/to/build/bin/mesh config.yaml
cd examples/heateq/steady/constant_conductivity/2d && /path/to/build/bin/heateq config.yaml
```

## Testing

Tests use GoogleTest and are registered via the `add_residuum_test(<name> <source> [libs...])` helper
in `tests/CMakeLists.txt`, which wires up `gtest_discover_tests` (30s timeout per test) and exposes
`TEST_DATA_PATH`/`TEST_OUTPUT_PATH` compile definitions for tests that touch files under `tests/data`.

```bash
cd build
ctest                          # run the full suite
ctest -R UnitTest_Fem_Basis_Lagrange # run tests matching a regex
ctest -R IntegrationTest_Full_Heateq_Cpu_HeatEquationMinimal --output-on-failure
./bin/tests/UnitTest_Linalg_Solver_Iterative_Cg_Solver   # run a single test binary directly
```

Unit tests under `tests/unit/` mirror `include/`'s path exactly, down to the header being tested
(e.g. `include/fem/form/FormRegistry.hpp` → `tests/unit/fem/form/FormRegistry.cpp`); a test that
exercises a family of sibling headers together (e.g. `fem/basis/Lagrange{1D,Quad,Hex}.hpp`) lives
at that family's directory rather than any one header. Integration tests under
`tests/integration/{fem,full}/<physics>/cpu` are nested by physics module first (mirroring
`equation/<physics>`), then by tier — `fem/` exercises the FEM layer directly (mesh → assembler,
no config/IO/driver), `full/` goes through the complete app pipeline (config, I/O, driver) — then
by backend. New tests must be added explicitly via `add_residuum_test(...)` in
`tests/CMakeLists.txt`; there is no automatic glob discovery. No CUDA compilation is wired into the
build yet.

## Architecture

### Library layering

Dependencies flow bottom-up; each is its own CMake static library target (`add_subdirectory` in
the root `CMakeLists.txt`):

- **`residuum_mesh`** (`src/mesh`) — mesh data structures, block-mesh generation
  (`mesh/generator/BlockMesh2D`), and Gmsh import/conversion (`mesh/exchange/gmsh`).
- **`residuum_io`** (`src/io`) — mesh/field I/O: `MeshIO`, `GmshReader`, `VTKWriter`, `YAMLReader`.
  Depends on `residuum_mesh` and `yaml-cpp`.
- **`residuum_expression`** (`src/utils`) — runtime-evaluated math expressions (via exprtk) used for
  YAML-specified BC/source expressions such as `"(1-x)*(1-y)"`. Scalar/Vector/Tensor variants.
- **`residuum_solver`** (`src/solver`) — YAML parsers for every solver-related config section
  (`DriverConfigParser`, `DiscretizationConfigParser`, `MeshConfigParser`,
  `LinearSolverConfigParser`, `NonlinearSolverConfigParser`, `TimeStepperConfigParser`,
  `OutputConfigParser`). Depends on `residuum_expression`.
- **`fem`, `linalg`, `equation`** — header-only (no dedicated `src/` library); see below.

Application libraries (`src/application/{mesh,heateq}`) each build a `residuum_*_app` static library
plus a thin `add_executable` wrapping `main()`.

### `linalg` — backend-templated numerics

Core containers (`linalg/types/{Vector,Matrix,CSRMatrix,DistributedVector,DistributedCSRMatrix}`)
are templated on `<T, Backend>`, where `Backend` is a tag type (`linalg/types/backend/{CPU,CUDA}.hpp`)
providing static primitives like `Backend::alloc<T>(n)`. This same CPU/CUDA backend-tag pattern
repeats through `fem/assembly/backend/`, `fem/boundary/backend/`, `linalg/operations/backend/`, and
`solver/{nonlinear,timestepper}/backend/`. **Only the CPU backend is currently implemented/used** —
CUDA directories exist as scaffolding but nothing enables CUDA in CMake yet.

`linalg/operator` defines the operator abstraction (`Operator`, `FEMOperator`, `CSROperator`) that
linear solvers act on. `linalg/solver` has direct (`direct/lu`) and iterative
(`iterative/{cg,bicgstab,gmres}`) solvers plus `preconditioner/`.

### `fem` — discretization building blocks

`fem/basis`, `fem/quadrature`, `fem/dof`, `fem/geometry` provide the reference-element machinery.
`fem/form` defines the generic weak-form abstractions (`BilinearForm`, `LinearForm`,
`NonlinearForm`, `NonlinearTangentForm`, `FormRegistry`) that physics modules implement.
`fem/assembly/Assembler` (templated, with a CPU specialization in `backend/cpu/Assembler.tpp`)
walks the mesh and calls into forms to build the global system. `fem/boundary` implements
essential/natural BC application (`BoundaryApplicator`, `EssentialBoundaryRegistry`,
`NaturalBoundaryRegistry`) with the same CPU-backend-`.tpp` split — this module was the subject of
the most recent commits on this branch (`boundary` refactor).

### `equation` — physics-specific forms

Each physics module (currently `equation/heateq`) implements the `fem/form` interfaces for its
PDE: `form/{DiffusionForm,MassForm,SourceForm,FluxBoundaryForm}`, plus `evaluator/` (element/field/
quadrature-point evaluation, `ConductivityModel`, `SourceFunction`) and `boundary/`
(`BoundaryValueFunction`, `BoundaryFluxFunction`). `equation/heateq/HeatEquation.hpp` ties these
together. Adding a new PDE means adding a sibling directory here plus a new `application/<eq>`.

### `solver` — config, drivers, time integration

`solver/config` holds one struct per YAML section (mirrors `examples/heateq/steady/constant_conductivity/2d/config.yaml`'s top-level keys:
mesh/discretization/physics/boundary_conditions/solver/output/logging); `solver/parser` (in
`src/solver`) turns YAML nodes into those structs. `solver/driver` provides the top-level solve
strategies (`Steady`, `Transient`, base concepts `Driver`/`TransientDriver`), `solver/linear` /
`solver/nonlinear` wrap the linalg solvers with problem-level config (`NewtonSolver`, etc.), and
`solver/timestepper` has explicit/implicit integrators (`ForwardEuler`, `BackwardEuler`, `RK`,
`GeneralizedAlpha`). `solver/stage/Stage.hpp` defines the generic per-timestep lifecycle
(`initialize` / `assemble` / `solve` / `finalize`) that application-level stages implement; a
driver only ever holds one such stage, so a multi-physics solve (e.g. a segregated
pressure-velocity stage followed by a turbulence stage) is expressed as a single composite stage —
`solver/stage/SegregatedStage.hpp` sweeps a heterogeneous sequence of `Stage`-conforming sub-stages
each outer iteration and itself satisfies `Stage`, so it composes with `Steady`/`Transient` like
any other stage. Time-integration (`Steady` vs `Transient`) and stage-composition (single-physics
vs `SegregatedStage`) are orthogonal axes.

### `application` — CLI entry points

Each app (`heateq`, `mesh`) follows the same shape:

1. `main()` (`src/application/<app>/<App>Application.cpp`) takes a single YAML path argument.
2. A `<App>ConfigParser` (in `parser/`, using `residuum_expression`/`residuum_io`/`residuum_solver`) reads it into
   a `config::<App>Config` struct.
3. `<App>Application::run()` hands the config to `<App>Dispatcher::run()`, whose job is to resolve
   *runtime* config choices (basis type, quadrature rule, backend) into a concrete instantiation of
   a *compile-time templated* `Stage<BackendType, BasisType, QuadratureVolumeType,
   QuadratureBoundaryType>` (see `application/heateq/stage/HeatStage.tpp`) and drive its
   `initialize/assemble/solve/finalize` lifecycle.

This dispatcher/stage layer is under active development on `feature/solver` — `HeatDispatcher` and
`HeatStage` are currently stub implementations (see recent commit history), so expect incomplete
behavior here rather than a bug in surrounding code.

### Extensibility priorities

When a design choice trades off between these, weigh them in this order (matches the goals above):

1. **CUDA backend completion.** The CPU/CUDA backend-tag pattern (`linalg/types/backend`,
   `fem/assembly/backend`, `fem/boundary/backend`, `solver/{nonlinear,timestepper}/backend`) is
   scaffolded throughout but only CPU is implemented. New abstractions should be shaped so a CUDA
   specialization is a natural drop-in, not a redesign.
2. **Distributed/parallel solves.** `linalg/types/{DistributedVector,DistributedCSRMatrix}` and
   `parallel/` should stay viable for MPI-style domain decomposition — avoid designs that quietly
   assume a single address space or single-rank ownership of the mesh/DOFs.
3. **New PDE/physics modules.** The `equation/<physics>` + `application/<physics>` pattern
   established by `heateq` should generalize to elasticity and RANS, and eventually to coupled
   multiphysics, without requiring a rewrite of `fem/form`/`fem/assembly`.
4. **New discretizations/solvers.** Basis functions, quadrature rules, linear/nonlinear solvers,
   and timesteppers should remain swappable via the existing config-driven dispatch
   (`solver/config` + `solver/parser` + `Dispatcher`/`Stage`).

### Conventions

**Namespaces and files**

- Namespaces mirror directory paths exactly, e.g. `include/fem/boundary/BoundaryApplicator.hpp` →
  `residuum::fem::boundary::BoundaryApplicator`. New files should follow the same mapping.
- One primary class/struct per header, filename matches it exactly (`ConductivityModel` lives in
  `ConductivityModel.hpp`).
- Header guards are `RESIDUUM_<PATH>_<NAME>_HPP`, built from the full path under `include/` plus
  the filename, e.g. `include/solver/stage/SteadyStage.hpp` → `RESIDUUM_SOLVER_STAGE_STEADYSTAGE_HPP`.
  No shortened or module-specific prefixes.
- Template implementations are split into `.tpp` files (included at the bottom of/alongside their
  `.hpp`), following the backend-specialization pattern under `backend/cpu/` (and eventually
  `backend/cuda/`).
- Introduce a new namespace only when it matches a new directory (namespace and directory always
  move together); don't nest an extra namespace inside a file purely for scoping.
- `ref/` is git-ignored local reference material (course PDFs, etc.) — not part of the shipped repo.

**Formatting**

- Indentation is tabs, not spaces — throughout, including inside a single function; don't mix.
- Includes: the file's own header first (for `.cpp`/`.tpp` files implementing a `.hpp`), then a
  blank line, then project headers grouped by module with a blank line between groups, then
  standard-library/external headers last. See any existing stage/parser file for a concrete example.

**Naming**

- Types (classes, structs, enums, enum values): `PascalCase` — `HeatProblem`, `GatherMode::Free`.
- Functions and methods: `camelCase` — `assembleMatrix`, `numFreeDOFs`.
- Local variables and parameters: `camelCase`.
- Private/internal data members: `camelCase` with a trailing underscore — `problem_`, `K_`,
  `jacobianForms_`.
- `static constexpr` values that are properties of a *type* rather than a runtime constant (used
  in concept checks, e.g. `T::Order`, `EvalQP::NumAuxStates`) are `PascalCase`, matching the
  register of a nested typedef rather than a value — `NumDOFs`, `SpatialDim`, `Order`.
- Template parameters take a trailing `T` — `ProblemT`, `FormsT`, `ModelT`, `QuadratureT` — rather
  than `Type`, no suffix, or a bare descriptive name.
- Acronyms stay fully capitalized inside an identifier: `DOF`, `BC`, `CSR`, `CG`, `SI` — never
  `Dof`/`Bc`. `numFreeDOFs`, `applyEssentialBCs`.

**Comments**

- At most one line per comment, and only when the *why* isn't obvious from the code itself (a
  hidden constraint, a workaround, a non-obvious invariant) — a comment should read as a plain
  sentence, not `like this -- with a dash aside`.
- Never narrate a change, a fix, or "this used to be X" — that belongs in the commit message, not
  the code, and rots as the codebase evolves. A comment describes the code as it stands today.
  Deeper design rationale belongs in this file (or a future `docs/` write-up), not a code comment.
- Never restate *what* the code already says through naming.

**Config parsers**

- `src/*/parser/*ConfigParser.cpp` go through the shared `io::YAMLReader::required<T>(node, key)` /
  `optional<T>(node, key, default)` helpers for every scalar field — never a hand-rolled
  `node["key"] ? ... : default` ternary or raw `.as<T>()` call.
- A required *sub-section* (as opposed to a scalar field) is checked directly (`const YAML::Node&
  x = node["x"]; if (!x) throw ...;`, see `SolverConfigParser`/`TimeStepperConfigParser`) since
  `YAMLReader` itself only handles leaf values.

## Working with Claude Code on this project

This is a solo-developer project. The developer is deliberately building deep familiarity with
FEM, numerical linear algebra, and performance engineering, and is using Claude Code as an
architectural collaborator rather than an autopilot. Default behavior in this repo:

- **Lead with architecture, not code.** For anything beyond a small, localized fix, propose
  candidate designs and walk through their tradeoffs (performance, backend-agnosticism, and the
  extensibility priorities above) before writing any implementation. Get alignment on the design
  first.
- **Headers, yes; hot-loop implementations, no (by default).** Once a design is settled, writing
  header/interface files (class shapes, template parameters, method signatures) is useful and
  expected — they'll always be reviewed before use. For `.cpp`/`.tpp` implementation, prefer
  **pseudocode over working code**, especially for performance-critical hot loops (assembly
  kernels, solver inner loops, threading/parallelism). The developer wants to write those
  themselves — don't fill them in unless explicitly asked to.
- **Stubs are fine.** It's OK for new modules to land as skeletons in the style of the current
  `HeatDispatcher`/`HeatStage` (see Architecture > `application`) rather than fully implemented in
  one pass — matching the existing WIP pattern is preferred over forcing completeness.
- **Keep edits narrowly scoped.** Prefer small, focused changes over sweeping refactors, matching
  the existing commit history style. Commit cadence/authorship is handled by the developer —
  don't take over `git commit` unless asked.
- **Build/test verification is the developer's to drive.** They typically compile and run tests
  themselves as part of the hands-on implementation work, so you don't need to always build and
  run the full suite after every change — but do flag anything you're unsure compiles or is
  unverified.
