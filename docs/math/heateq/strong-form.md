# Heat equation, strong form

## Governing equation

The heat equation in its transient form is

$$
\rho c_p \frac{\partial T}{\partial t} - \nabla \cdot (\boldsymbol{\kappa} \nabla T) = f \quad \text{in } \Omega
$$

where $T$ is temperature, $\rho$ is density, $c_p$ is specific heat, $\boldsymbol{\kappa}$ is the
conductivity tensor, and $f$ is a volumetric heat source. The steady case follows by setting
$\partial T/\partial t = 0$, leaving

$$
-\nabla \cdot (\boldsymbol{\kappa} \nabla T) = f \quad \text{in } \Omega
$$

The conductivity $\boldsymbol{\kappa}$ may be a scalar (isotropic), a general symmetric tensor
(anisotropic), or a function of temperature itself, in which case the equation is nonlinear in
$T$.

## Boundary conditions

The domain boundary $\Gamma = \partial\Omega$ splits into a Dirichlet portion $\Gamma_D$ and a
Neumann portion $\Gamma_N$, with $\Gamma = \Gamma_D \cup \Gamma_N$.

An essential (Dirichlet) boundary condition prescribes temperature directly:

$$
T = \bar{T} \quad \text{on } \Gamma_D
$$

A natural (Neumann) boundary condition prescribes heat flux:

$$
-\boldsymbol{\kappa}\nabla T \cdot \mathbf{n} = \bar{q} \quad \text{on } \Gamma_N
$$

where $\mathbf{n}$ is the outward unit normal. Both $\bar{T}$ and $\bar{q}$ may vary over space
and time.

## Initial condition

The transient problem additionally requires an initial temperature field:

$$
T(\mathbf{x}, 0) = T_0(\mathbf{x}) \quad \text{in } \Omega
$$

The steady problem has no initial condition. It is a boundary value problem in space only.
