# Heat equation, discretization

## Spatial discretization

The Galerkin approach approximates the temperature field as $T \approx T_h = \sum_a N_a(\mathbf{x})
T_a$, where the $N_a$ are basis functions and the physical element geometry is reached through an isoparametric mapping.
Substituting this approximation into the weak form and taking $v = N_a$ for each basis function in turn produces the element level stiffness and mass matrices:

$$
K_{ab} = \int_{\Omega_e} \boldsymbol{\kappa} \nabla N_a \cdot \nabla N_b \, d\Omega, \qquad
M_{ab} = \int_{\Omega_e} \rho c_p N_a N_b \, d\Omega
$$

Both integrals are approximated numerically and accumulated over every element into the global $\mathbf{K}$, $\mathbf{M}$, and load vector $F$ by
the assembler.

## Time discretization

Backward Euler replaces the time derivative in the semi-discrete transient weak form with a backward difference:

$$
\dot{T} \approx \frac{T^{n+1} - T^n}{\Delta t}
$$

which gives the fully discrete system for one step:

$$
\left(\frac{\mathbf{M}}{\Delta t} + \mathbf{K}\right) T^{n+1}
= F^{n+1} + \frac{\mathbf{M}}{\Delta t} T^n
$$

This is first order accurate in time. Higher order schemes such as an explicit Runge-Kutta method
or a generalized-alpha integrator occupy the same timestepper interface and can be added without
changing anything in the spatial discretization above.

## Assembling a solve by mode

Crossing the steady/transient axis with the linear/nonlinear axis gives four combinations. Newton
is the nonlinear solver currently supported, so the nonlinear cases below are described in terms
of Newton's residual and tangent.

### Steady, linear

The element matrices assemble once into a single global system:

$$
\mathbf{K} T = F
$$

which is solved once. This is the path `SteadyStage` takes when the conductivity model does not
depend on temperature.

### Steady, nonlinear

When $\mathbf{K}$ depends on $T$, the system above is no longer solved directly. A residual is
formed instead:

$$
R(T) = \mathbf{K}(T) T - F = 0
$$

and driven to zero by Newton's method. Each iteration solves

$$
J \, \Delta T = -R(T^k), \qquad T^{k+1} = T^k + \Delta T
$$

for the tangent $J = \partial R/\partial T$. The stiffness matrix itself, evaluated at the current
temperature, uses the same formula as $K_{ab}$ above with $\boldsymbol{\kappa}$ replaced by
$\boldsymbol{\kappa}(T)$. Differentiating the residual with respect to nodal temperature produces one
additional term on top of that, since $\boldsymbol{\kappa}$ itself now depends on $T$:

$$
(K_T)_{ab} = \int_{\Omega_e} \frac{\partial \boldsymbol{\kappa}}{\partial T} \, N_b \,
\left(\nabla N_a \cdot \nabla T\right) \, d\Omega
$$

so that the full tangent is

$$
J = \mathbf{K}(T) + \mathbf{K}_T(T)
$$

`TangentDiffusionForm` computes only $(K_T)_{ab}$, added on top of what `DiffusionForm` already
computes for $K_{ab}(T)$.

### Transient, linear

The Backward Euler system above is exactly a steady linear solve with a modified operator and
right hand side: the effective stiffness becomes $\mathbf{M}/\Delta t + \mathbf{K}$, and the right
hand side gains the $\mathbf{M}T^n/\Delta t$ contribution from the previous step.
`BackwardEulerStage` assembles this system once per step and solves it once, the same way
`SteadyStage` does for the steady linear case.

### Transient, nonlinear

Combining both axes, the per-step residual is

$$
R(T^{n+1}) = \mathbf{M}(T^{n+1}) \dot{T}^{n+1} + \mathbf{K}(T^{n+1}) T^{n+1} - F^{n+1} = 0
$$

with $\dot{T}^{n+1}$ evaluated through the same Backward Euler relation and driven to zero by
Newton within each step. The tangent now includes both the steady tangent stiffness contribution
above and an analogous transient mass contribution:

$$
J = \mathbf{K}(T^{n+1}) + \mathbf{K}_T(T^{n+1})
+ \frac{\mathbf{M}(T^{n+1})}{\Delta t} + \mathbf{M}_T(T^{n+1})
$$

The mass matrix itself, evaluated at the current temperature, uses the same formula as $M_{ab}$
above with $c_p$ replaced by $c_p(T)$. Specific heat can depend on temperature just as conductivity
can, and differentiating the mass term of the residual with respect to nodal temperature, using
$\partial \dot{T}/\partial T = 1/\Delta t$ from the Backward Euler relation, produces one
additional term:

$$
(M_T)_{ab} = \int_{\Omega_e} \rho \, \frac{\partial c_p}{\partial T} \, \dot{T} \, N_a N_b \, d\Omega
$$

This term is scaled by the current rate of change $\dot{T}$ rather than by a gradient, since the
mass term itself has no spatial derivative to differentiate. `TangentMassForm` computes only
$(M_T)_{ab}$, added on top of what `MassForm` already computes for $M_{ab}(T)/\Delta t$.

`BackwardEulerStage` handles this general case when paired with a Newton-capable nonlinear solver
runner rather than a direct linear solve.

## Linear system properties

For an isotropic or general symmetric conductivity, the resulting operator ($\mathbf{K}$ alone, or
$\mathbf{M}/\Delta t + \mathbf{K}$) is symmetric positive definite, and CG is an appropriate linear
solver. An anisotropic conductivity tensor can make the operator non-symmetric, in which case GMRES
is required instead.
