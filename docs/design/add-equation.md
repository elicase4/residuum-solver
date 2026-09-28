# Adding a new equation

<!-- Outline — topics to discuss/draft, not final content. -->

- Concepts a new equation must satisfy (`fem::form::*`, `fem::evaluator::EvalModel`,
  `EvalQuadraturePointVolume`/`Boundary`, `fem::quantity::QuantityForm`)
- Directory shape (`form/`, `evaluator/`, `boundary/`, `quantity/`), using `equation/heateq` as
  the reference
- Config → parser → struct chain, `application/<physics>/config` + `*ConfigParser`
- Wiring `application/<physics>` (`<Physics>Dispatcher`)
- Test coverage expectations (unit mirrors `equation/<physics>`, integration tiers)
