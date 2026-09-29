# Adding a new equation

Adding a new PDE to Residuum means adding a new `equation/<physics>` module and a matching
`application/<physics>` module. Neither the finite element assembly layer nor the solver layer
needs to change. `equation/heateq` is the reference implementation to read alongside this guide.

## Concepts a new equation must satisfy

A new equation module provides implementations of the necessary `fem::form` interfaces. 
It also provides an `fem::evaluator::EvalModel` implementation describing how material properties are evaluated at a quadrature point, together with the matching `EvalQuadraturePointVolume` and `EvalQuadraturePointBoundary` types. Quantities meant for monitoring or output, such as a boundary flux integral, are expressed through `fem::quantity::QuantityForm`.

Each of these interfaces is a concept, not a base class. A form or evaluator that does not satisfy
the required concept fails to compile at the point it is used, with the compiler pointing at the
specific member function or type alias that is missing.

## Directory shape

An equation module mirrors the shape of `equation/heateq`: a `form/` directory for the weak form
implementations, an `evaluator/` directory for quadrature point evaluation and material models, a
`boundary/` directory for boundary value and flux functions, and a `quantity/` directory for
monitorable quantities. A single header at the module root, analogous to `HeatEquation.hpp`, ties
these pieces together into the bundle type for use throught the codebase.

## Configuration

Every equation module has its own configuration struct under `application/<physics>/config`,
populated by a corresponding `application/<physics>/parser/*ConfigParser`. These parsers read the
equation-specific sections of a YAML file (material properties, boundary condition definitions,
source terms) using the shared `io::YAMLReader::required`/`optional` helpers.

## Wiring the application

The final piece is a `<Physics>Dispatcher`, which resolves the runtime discretization choices in a
parsed configuration into a concrete, templated `Stage` and drives it. This is the only place
equation-specific code touches the dispatch machinery in `fem::dispatch`. Everything else in that
machinery is equation agnostic.

## Test coverage

Unit tests for an equation module mirror its directory structure exactly, down to the header being
tested. Integration tests are organized by tier: tests under `fem/` exercise the assembly pipeline
directly from a mesh, without going through configuration parsing or the application layer, while
tests under `full/` exercise the complete pipeline from a YAML configuration through to a driver.
