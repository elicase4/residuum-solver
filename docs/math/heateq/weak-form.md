# Heat equation — weak form

<!-- STUB. Fill in below. -->

## Derivation

<!-- Multiply by test function $v \in H^1_0(\Omega)$ (or the appropriate trial/test space), 
integrate over $\Omega$, integrate the diffusion term by parts to move a derivative onto $v$ and
produce the boundary flux term. State the function spaces for trial ($T$) and test ($v$)
explicitly, including how the essential BC is handled (lifting function vs. constrained space —
match whichever `fem::boundary::EssentialBoundaryRegistry` actually implements). -->

## Resulting bilinear/linear forms

<!-- Write out $a(T, v)$ and $L(v)$ (steady) or the semi-discrete transient residual, then map
each term to its C++ form class:

| Term | Form class |
|---|---|
| $\int_\Omega \mathbf{K}\nabla T \cdot \nabla v \, d\Omega$ | `DiffusionForm` |
| $\int_\Omega \rho c_p \, \dot{T} \, v \, d\Omega$ | `MassForm` |
| $\int_\Omega s \, v \, d\Omega$ | `SourceForm` |
| $\int_{\Gamma_N} q \, v \, d\Gamma$ | `FluxBoundaryForm` |

and note where the nonlinear tangent terms ($\partial \mathbf{K}/\partial T$, etc. —
`TangentDiffusionForm`/`TangentMassForm`) come from for the temperature-dependent case. -->

## Nonlinear (temperature-dependent) case

<!-- Newton linearization: residual $R(T) = 0$, Jacobian $J = \partial R/\partial T$ — connect to
`solver/nonlinear/Newton.hpp` and the `K_T`/`M_T` tangent forms. -->
