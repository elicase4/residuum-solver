<p align="center">
  <img src="docs/assets/logo.png" width="500" alt="Residuum">
</p>

<h1 align="center"></h1>

# Residuum

A modular finite element method (FEM) solver written in modern C++20, built around a backend-agnostic execution, concept-constrained template architecture, and YAML configuration at the application level.

## Highlights

- **Concept Constrained Core**: Core numerical modules are concept constrained, eliminating virtual dispatch overhead, while enforcing self-documenting type enforcement at compile time.
- **Backend Agnostic Design**: Core containers and kernels are templated on a `Backend` tag to allow extensibility for application execution on different hardware.
- **Runtime Dispatching**: A YAML config resolves solver settings at runtime to easily run mutliple solver configurations for an engineering analysis.
- **Composable solver architecture.**: `Steady`/`Transient` drivers, single-physics/multi-physics composed into a solver `Stage`.

## Status

Actively developed, solo project. The heat equation (`equation/heateq`) is the current focus and
is functionally complete: steady and transient (Backward Euler) solves, linear and nonlinear
(temperature-dependent) material models, essential/natural boundary conditions, file- and
expression-driven fields, monitors/quantities, VTK output. 

CUDA backend and MPI program execution are next on the roadmap, followed by Elasticity and an incompressible RANS solver.

## Quick start

```bash
git clone https://github.com/elicase4/residuum-solver.git && cd residuum-solver

cmake --preset release
cmake --build --preset release --parallel
export PATH="$PWD/build-release/bin:$PATH"

cd examples/mesh/square && mesh config.yaml && cd -
cd examples/heateq/steady/constant_conductivity/2d && heateq config.yaml && cd -
```

Dependencies: a C++20 compiler, CMake ≥ 3.18 
[`yaml-cpp`](https://github.com/jbeder/yaml-cpp),
[`GoogleTest`](https://github.com/google/googletest) (for the test suite), and
[`exprtk`](https://github.com/ArashPartow/exprtk) (header-only).

Every physics run is a YAML file. Mesh, materials, boundary conditions, solver, and output are
all declared within a configuration, see below snippet as an example.

```yaml
materials:
  conductivity:
    type: constant
    value: 1.0
    unit: W/(m*K)

boundary_conditions:
  - boundary: 0
    type: value
    read:
      mode: expression
      expression: "(1-x)*(1-y)"
      unit: K
    forms: [value_bc]
```

See `examples/` for steady, transient, nonlinear-conductivity, anisotropic-conductivity, and
Gmsh-imported-geometry configs with the Heat Equation. For a full walkthrough of one example end
to end, see the [transient heating tutorial](docs/examples/heating-block-tutotial/heating-block-tutorial.md).

## Architecture

```
equation/<physics>    PDE-specific weak forms, models, boundary functions (currently: heateq)
application/<app>     CLI entry points; config parsing → Dispatcher → templated Stage
solver/               Drivers (Steady/Transient), stages, linear/nonlinear solvers, timesteppers
fem/                  Assembly, basis functions, quadrature, boundary application, quantities
linalg/               Backend-templated containers, operators, iterative/direct solvers
mesh/, io/            Mesh generation/import (Gmsh), field & VTK I/O
```

Dependencies flow strictly bottom-up. 
`fem`, `linalg`, and `equation` are header-only, with template implementations living in `.tpp` files alongside their headers. 
`mesh`, `io`, and `solver config parsers compile into static libraries. 
A CMake package config (`find_package(Residuum)`) is available for downstream consumption — see `cmake/ResiduumConfig.cmake.in`.

## Testing

```bash
cd build
ctest                                          # full suite
ctest -R UnitTest_Fem_Basis_Lagrange           # by name/regex
ctest -R IntegrationTest_Full_Heateq_Cpu --output-on-failure
```

Unit tests mirror `include/`'s directory structure exactly. Integration tests are nested by physics module, then by tier, then by backend.

## License

MIT — see [LICENSE](LICENSE).
