# Heat equation — strong form

<!-- STUB. Fill in below. -->

## Governing equation

<!-- Transient form, with the steady case as $\partial T/\partial t = 0$:

$$
\rho c_p \frac{\partial T}{\partial t} - \nabla \cdot (\mathbf{K} \nabla T) = s \quad \text{in } \Omega
$$

State which conductivity model this covers (constant / anisotropic-tensor / temperature-dependent
— see `equation/heateq/evaluator/ConductivityModel.hpp`'s Dependence x Symmetry axes) and note
where each term maps to a config field (`materials.conductivity`, `materials.specific_heat`,
`materials.density`, `physics.models.source`). -->

## Boundary conditions

<!-- Essential (Dirichlet) and natural (Neumann/flux) — map directly to
`boundary_conditions[].type: value` / `type: flux` in the YAML schema, and to
`fem::boundary::{EssentialBoundaryRegistry, NaturalBoundaryRegistry}` in code. -->

## Initial condition

<!-- Transient only — maps to `initial_condition` in YAML. -->
