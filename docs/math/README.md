# Math documentation

Derivations per equation, in Markdown with embedded LaTeX. One subdirectory per physics module, mirroring `equation/<physics>`.

## Conventions

- **Notation:** : bold for vectors/tensors, $\Omega$ for the domain, $\Gamma$ for its
  boundary, $\Gamma_D$/$\Gamma_N$ for Dirichlet/Neumann portions.
- **Units:**: SI throughout, matching the `unit:` fields required in YAML configs.
- **Discretization symbols:**: Consistent with C++ identifiers $\mathbf{K}$ for the stiffness matrix matches, `SteadyStage::K_`, $\mathbf{F}$ for the load vector matches `HeatProblem::F()`.

## Index

- [`heateq/`](heateq/) — heat equation: [strong form](heateq/strong-form.md), [weak form](heateq/weak-form.md), [discretization](heateq/discretization.md).
