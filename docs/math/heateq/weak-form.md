# Heat equation, weak form

## Derivation

Let $v$ be a test function drawn from an appropriate test space, vanishing on $\Gamma_D$.
Multiplying the strong form by $v$ and integrating over the domain gives

$$
\int_\Omega \rho c_p \frac{\partial T}{\partial t} v \, d\Omega
- \int_\Omega \nabla \cdot (\mathbf{\kappa} \nabla T) \, v \, d\Omega
= \int_\Omega f v \, d\Omega
$$

Integrating the diffusion term by parts moves one derivative from $T$ onto $v$ and produces a
boundary term,

$$
- \int_\Omega \nabla \cdot (\mathbf{\kappa} \nabla T) \, v \, d\Omega
= \int_\Omega \mathbf{\kappa}\nabla T \cdot \nabla v \, d\Omega
- \int_\Gamma (\mathbf{\kappa}\nabla T \cdot \mathbf{n}) \, v \, d\Gamma
$$

On $\Gamma_D$ the test function vanishes, so only the $\Gamma_N$ portion of the boundary integral
survives. Substituting the natural boundary condition turns that surviving term into a known flux
contribution.

The trial temperature $T$ is drawn from a space satisfying the essential boundary condition, and
the test function $v$ is drawn from the corresponding space with a homogeneous condition in its
place. In practice this is enforced by constraining the algebraic degrees of freedom associated
with $\Gamma_D$ directly, rather than by an explicit lifting function.

## Resulting forms

The weak form of the steady problem is: find $T$ such that for all admissible $v$,

$$
a(T, v) = L(v)
$$

with

$$
a(T, v) = \int_\Omega \mathbf{\kappa}\nabla T \cdot \nabla v \, d\Omega, \qquad
L(v) = \int_\Omega f v \, d\Omega + \int_{\Gamma_N} \bar{q} \, v \, d\Gamma
$$

The transient problem adds the mass term to the left-hand side,

$$
\int_\Omega \rho c_p \dot{T} v \, d\Omega + a(T, v) = L(v)
$$

Each term corresponds to a form class: the diffusion term above to `DiffusionForm`, the mass term
to `MassForm`, the volumetric source term to `SourceForm`, and the boundary flux term to
`FluxBoundaryForm`.

## Nonlinear case

When $\mathbf{\kappa}$ and/or $c_p$ depend on $T$, the diffusion and mass terms are nonlinear and the discrete problem is
instead written as a residual to be driven to zero,

$$
R(T) = \int_\Omega \rho c_p(T) \dot{T} v \, d\Omega
+ \int_\Omega \mathbf{\kappa}(T)\nabla T \cdot \nabla v \, d\Omega - L(v) = 0
$$

A Newton solver linearizes this residual at each iteration, requiring the Jacobian
$J = \partial R/\partial T$. Differentiating the diffusion term with respect to $T$ produces a
tangent contribution from $\partial \mathbf{\kappa}/\partial T$ in addition to the term already present
when $\mathbf{\kappa}$ is constant, corresponding to the tangent diffusion and tangent mass form
classes.
