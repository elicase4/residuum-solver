# Contributing

Residuum is a solo-developer research project.

Coding conventions live in `CLAUDE.md`, useful to read first. This file points at the deeper guides:

- [`docs/design/architecture.md`](docs/design/architecture.md) — how a solve actually happens, end to end: YAML config → `Dispatcher` → templated `Stage` → `Assembler` → linear solver.
- [`docs/design/add-equation.md`](docs/design/add-equation.md) — adding a new PDE (`equation/<physics>`): which `fem/form`/`fem/evaluator` concepts to implement, and what `equation/heateq` looks like as a worked example.
- [`docs/design/add-solver.md`](docs/design/add-solver.md) — adding a new `Stage`, linear operator, preconditioner, nonlinear solver, or timestepper.
- [`docs/design/benchmarking.md`](docs/design/benchmarking.md) — the `benchmark/` module: micro and macro kernels, the `Result`/`CSVReporter` harness, and how a new one fits in.

Before opening a PR: CI runs automatically against it with Debug, Release, and Sanitizer job. `gh pr checks` shows the status from the terminal.
