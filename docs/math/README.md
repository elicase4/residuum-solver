# Math documentation

Derivations per equation, in Markdown with embedded LaTeX (GitHub renders `$...$`/`$$...$$` inline
and display math directly). One subdirectory per physics module, mirroring `equation/<physics>`.

## Conventions

<!-- Fill in as they're settled — notation should stay consistent across equations so a reader
     moving from heateq to elasticity doesn't have to relearn symbols. -->

- **Notation:** *(e.g. bold for vectors/tensors, $\Omega$ for the domain, $\Gamma$ for its
  boundary, $\Gamma_D$/$\Gamma_N$ for Dirichlet/Neumann portions — pin down here once decided.)*
- **Units:** SI throughout, matching the `unit:` fields required in YAML configs
  (`materials.conductivity.unit`, etc.) — see `include/utils/logging` unit handling and
  `fem/quantity/QuantityUnits.hpp`.
- **Discretization symbols:** try to keep these consistent with the actual C++ identifiers where
  it doesn't fight standard PDE notation — e.g. $\mathbf{K}$ for the stiffness matrix matches
  `SteadyStage::K_`, $\mathbf{F}$ for the load vector matches `HeatProblem::F()`.

## Index

- [`heateq/`](heateq/) — heat equation: [strong form](heateq/strong-form.md),
  [weak form](heateq/weak-form.md), [discretization](heateq/discretization.md).
