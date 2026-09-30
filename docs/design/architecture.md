# Architecture

Residuum is organized as a small set of layered libraries, each built as its own static library
and consumed bottom up. Mesh and IO form the base. Linear algebra and finite element machinery sit
above them. Physics equations build on the finite element layer. The solver layer drives equations
through time or to convergence. Application binaries sit on top, translating a YAML configuration
into a running solve.

## A solve, traced end to end

A typical solve starts from a YAML configuration file. A parser reads it into a plain
configuration struct. A dispatcher reads the runtime choices in that struct (basis type,
quadrature rule, backend) and resolves them into a concrete, fully templated `Stage`. The stage's
`initialize`, `assemble`, `solve`, and `finalize` lifecycle drives the actual solve: assembling
matrices and vectors through the `Assembler`, solving the resulting linear system through a
`LinearOperator` and linear solver, and producing output.

The split between the dispatcher and the stage is what lets the hot assembly and solve loop stay
branch free. Once the dispatcher has resolved a configuration into a concrete template
instantiation, every subsequent call inside that instantiation is monomorphic. The dispatcher
pays the runtime cost once, at startup, and the stage pays none of it afterward.

## Concepts over virtual dispatch

Every extension point in the assembly and solver pipeline (quadrature point evaluators, weak
forms, linear operators, solver drivers, preconditioners) is expressed as a C++20 concept, checked
at the actual call site rather than left as a comment or a naming convention. A type that does not
satisfy the required concept fails to compile with a message pointing at the specific requirement
it is missing.

## Backend-agnostic design

Core containers and kernels are templated on a backend tag rather than hardcoding a particular
hardware target. `linalg::types::backend::CPU` is the only backend implemented today.

## Runtime dispatch, compile time templates

Basis type, quadrature rule, and backend are all runtime choices from the user's point of view:
they are set in a YAML file, not chosen at compile time by whoever is writing code. Internally,
the assembly and solve loops want to be fully templated on all three, so that the compiler
can specialize and inline aggressively. `fem::dispatch` bridges this gap. It is the one place in
the codebase where a runtime value gets matched against a fixed set of compile time possibilities
and used to instantiate the corresponding template. Everything downstream of that dispatch is
compile time specialized.

## Multi-physics composition

A driver such as `Steady` or `Transient` owns exactly one stage. For multi-physics problems, such
as a segregated pressure-velocity solve followed by a turbulence closure, `SegregatedStage`
composes a heterogeneous sequence of stages and itself satisfies the same `Stage` concept its
members do. Physics composition happens entirely at the stage level, and the time integration
and stage composition axes stay orthogonal.

## Benchmarking

`benchmark/` sits alongside `tests/` as a second kind of verification, aimed at performance rather
than correctness. Micro benchmarks are small, argument-free binaries exercising a single linear
algebra kernel against a synthetic problem. Macro benchmarks reuse the real application pipeline
directly, running an actual equation dispatcher against a real configuration, to measure how a
change behaves in the context it will actually run in. See [benchmarking.md](benchmarking.md) for
the detailed design.
