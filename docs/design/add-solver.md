# Adding a new solver method

There are several distinct extension points in the solver layer, and a new solver method usually
only touches one of them.

## Extension point application

A `Stage` implements the per-iteration lifecycle (`initialize`, `assemble`, `solve`, `finalize`)
for a particular equation and driver combination. A `LinearOperator` implements how a matrix or
matrix-free operator is applied to a vector. A linear solver, such as CG or GMRES, implements an
iterative method against any type satisfying `LinearOperator`. A nonlinear solver, such as Newton,
drives a sequence of linear solves. A timestepper advances a transient problem by one step. A
driver, such as `Steady` or `Transient`, owns the outer loop around a stage. These are independent
axes, and a new solver method should be placed at the axis it actually belongs to rather than
folded into an existing one.

## Stage and its variants

`SteadyStage` is the simplest case: assemble once, solve once. `BackwardEulerStage` adds the mass
matrix contribution and the residual construction a first-order implicit timestepper needs.
`SegregatedStage` composes a sequence of other stages and itself satisfies the `Stage` concept,
which is what allows a multi-physics problem to be driven by the same `Steady`/`Transient` drivers
as a single-physics one. A stage that also supports nonlinear iteration or transient time stepping
additionally satisfies `NonlinearCapableStage` or `TransientCapableStage`, which extend the base
concept with the additional hooks a nonlinear solver or timestepper needs to call into.

## LinearOperator

`CSROperator` and `FEMOperator` are the two currently supported shapes: an explicitly assembled CSR matrix,
or a matrix-free operator that assembles its action on the fly. Both satisfy the same
`LinearOperator` concept, checked with a `static_assert` placed inside a constructor.

## Preconditioner

`Identity` and `Jacobi` are the two currently supported preconditioners. Every preconditioner is called
through the same `update(const OperatorT&)` hook once per solve, which lets it pull whatever it
needs from the current operator, such as a diagonal Jacobi preconditioner re-extracting the
diagonal after the operator has changed between Newton iterations. An operator that wants to
support a diagonal-based preconditioner satisfies `DiagonalExtractable`, a small concept requiring
a `diagonal(VectorT&)` member. A preconditioner should be default constructible where possible,
sizing itself lazily on first use, rather than requiring its size up front.

## Wiring into the factory

Linear solvers, nonlinear solvers, and timesteppers are each resolved from configuration through
their own factory: `LinearSolverFactory`, `NonlinearSolverFactory`, and `TimeStepperFactory`.
`LinearSolverFactory` in particular needs two nested runtime-to-compile-time switches, one over the
solver type and one over the preconditioner type, since neither varies through virtual dispatch. A
runtime configuration value has to be turned into a compile time type choice before the
corresponding class template can be instantiated.
