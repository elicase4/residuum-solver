# Benchmarking

<!-- Outline — topics to discuss/draft, not final content. -->

- Micro (`benchmark/linalg/`) vs. macro (`benchmark/equation/<physics>/`) — argument-free
  synthetic-problem binaries vs. reusing the real application pipeline directly
  (`HeatDispatcher::run(config)`), and why they need different harnesses
- Directory shape mirrors library layers, not headers (a kernel is inherently cross-cutting —
  operator + solver + preconditioner at once — so there's no single header it belongs under)
- Core harness (`include/benchmark/core/`): `Timer`, `Result`, `Runner::run`, `CSVReporter`; why
  the CSV is long-format (one row per axis combination) rather than wide, and how that trades off
  against plotting convenience (splitting into series happens at consumption time, not here)
- The `backend`/`equation` axis columns exist on `Result` from day one so CUDA and
  equation-specific macro benchmarks are additive rows later, not a schema change
- Macro-kernel fixtures (mesh files at different resolutions, generated once via the `mesh` app
  into the benchmark's own `fixtures/` directory) vs. reusing `examples/` directly, and why
- `CMakePresets.json` (`dev`/`release`/`benchmark`) — why benchmarking needs a separate
  optimized build from the Debug dev loop, and why the `benchmark` preset still needs
  `BUILD_APPLICATIONS=ON` despite the name
- Open direction: h/p-convergence studies and solver convergence-rate verification (e.g.
  quadratic Newton) — why these don't fit today's timing-shaped `Result`/`CSVReporter` and likely
  want a parallel type instead of a repurposed timing row
