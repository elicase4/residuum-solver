# Heat equation — discretization

<!-- STUB. Fill in below. -->

## Spatial discretization

<!-- Galerkin FEM: $T \approx T_h = \sum_a N_a T_a$, basis functions $N_a$ — connect to
`fem::basis::Lagrange{1D,Quad,Hex}` and the isoparametric mapping in
`fem::geometry::JacobianTransform`. Element matrices/vectors:

$$
K_{ab} = \int_{\Omega_e} \mathbf{K} \nabla N_a \cdot \nabla N_b \, d\Omega, \qquad
M_{ab} = \int_{\Omega_e} \rho c_p N_a N_b \, d\Omega
$$

and how quadrature approximates these integrals (`fem::quadrature::GaussQuadrature*`,
`fem::assembly::Assembler`). -->

## Time discretization

<!-- Backward Euler as currently implemented (`solver::timestepper::BackwardEuler`):

$$
\mathbf{M} \frac{T^{n+1} - T^n}{\Delta t} + \mathbf{K} T^{n+1} = F^{n+1}
$$

Note this is where `Udot_` (the rate-of-change aux state threaded through `Assembler`/
`FormRegistry`) comes from, and what changes for a future higher-order scheme
(`solver/timestepper/{ForwardEuler,RK,GeneralizedAlpha}.hpp` are scaffolded but not yet
implemented for heateq). -->

## Linear system

<!-- Resulting $\mathbf{K}T = F$ (steady) or the per-timestep linear solve — note SPD-ness for the
isotropic/constant case (why CG works) and non-symmetry for the anisotropic case (why GMRES,
task #121, is needed for `examples/heateq/steady/anisotropic_conductivity`). -->
