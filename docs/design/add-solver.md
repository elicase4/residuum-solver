# Adding a new solver method

<!-- Outline — topics to discuss/draft, not final content. -->

- Which extension point actually applies: `Stage`, `LinearOperator`, linear solver, nonlinear
  solver, timestepper, or driver — these are distinct concepts, not one path
- `Stage`/`NonlinearCapableStage`/`TransientCapableStage` — what each requires, worked examples
  (`SteadyStage`, `BackwardEulerStage`, `SegregatedStage` as the composite case)
- `LinearOperator` — `CSROperator`/`FEMOperator` as the two shapes, the `static_assert` pattern
- Preconditioner — `Identity`/`Jacobi` as the two shapes; the uniform `update(const OperatorT&)`
  hook every preconditioner gets called through once per solve; `DiagonalExtractable` as the
  concept an operator must satisfy for a diagonal-based preconditioner (`CSROperator`/`FEMOperator`
  both do); why `Jacobi` is lazily-sized (`VectorT` has no default constructor) rather than taking
  a size in its constructor
- Wiring into the relevant factory/config (`LinearSolverFactory`, `NonlinearSolverFactory`,
  `TimeStepperFactory`) — `LinearSolverFactory` specifically needs two nested runtime-to-compile-time
  switches (solver type, then preconditioner type) since neither varies via virtual dispatch
- When a new `Driver` is actually warranted vs. composing a `SegregatedStage`
