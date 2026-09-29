# Benchmarking

`benchmark/` sits alongside `tests/` and measures performance rather than correctness.

## Micro and macro benchmarks

A micro benchmark is a small, argument-free binary that exercises a single linear algebra kernel,
such as a sparse matrix-vector product or an iterative solver, against a synthetic problem built
directly in code. A macro benchmark instead reuses the real application pipeline, running an
actual equation dispatcher against a real configuration, so that a measurement reflects the cost
of the full pipeline rather than an isolated kernel. 

## Directory shape

Benchmark kernels are organized by library layer, under `benchmark/linalg/` and
`benchmark/equation/<physics>/`, rather than mirrored one to one against a single header the way
unit tests are. 

## Core harness

The shared harness lives under `include/benchmark/core/`. `Timer` measures wall clock time.
`Runner::run` performs a warmup phase followed by a measured phase and returns a `Result`, which
records the kernel name, a set of named axis values, timing statistics, and an optional floating
point operation estimate used to derive a throughput figure. `CSVReporter` appends a `Result` to a
CSV file, validating that the axis columns of a file stay consistent across appends. Every benchmark
result carries a `backend` and an `equation` axis, even when only one backend or no specific equation 
currently applies. 

## Fixtures

A macro benchmark that needs mesh files at multiple resolutions generates them once, through the
`mesh` application, into its own `fixtures/` directory rather than reusing the meshes under
`examples/`. This keeps a benchmark's inputs under its own control, independent of anything an
example might need to demonstrate.

## Build configuration

`CMakePresets.json` provides a `dev` preset for day to day debug builds, a `release` preset for an
optimized build of everything, and a `benchmark` preset for an optimized build scoped to
benchmarking. The `benchmark` preset still needs application libraries built, since macro kernels
link directly against them, even though it does not build the standalone application executables
themselves.

## Extending the harness

Some measurements do not fit the timing-shaped `Result` described above. A convergence study,
verifying that a mesh refinement sequence produces the expected order of accuracy or that a
nonlinear solver converges at its expected rate, measures an error norm or a residual ratio rather
than a wall clock time. These additional benchmarking axes are expected to be added in future versions
