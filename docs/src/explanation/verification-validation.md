# Verification and validation

The tutorials and the test suite carry out four different activities, and their
results mean different things.

- **Code verification** asks whether the implementation solves the intended discrete
  equations. Examples are exact view factors of catalog geometries, power conservation
  of the ray trace, symmetry selection rules that must hold to round-off, and
  derivatives checked against finite differences.
- **Solution verification** asks how large the discretization error of a given
  calculation is. Examples are the closure error of a view-factor matrix and its
  convergence with mesh refinement, and the agreement of the fast deposition operator
  with the ray trace.
- **Validation** asks whether the model reproduces an independent reference, either
  physical or analytic, for the intended use. Examples are the integrating-sphere
  multiplier, Stokes-type analytic profiles, and, eventually, comparisons with VisRad
  calculations and shot data.
- **Calibration** fits model parameters to reference data. The cross-beam energy
  transfer model is calibrated, not validated: its multiplier, flow and saturation
  clamp are fitted to measured transfer.

A tutorial states which of these it does in its closing section, "What this checks".

## Where each check lives

| Check | Activity | Location |
|:--|:--|:--|
| BVH ray casting against brute force | code verification | `test/runtests.jl`, "ray casting" |
| Coaxial-disk view factor, occlusion, point-to-polygon kernel | code verification | "view factors" |
| Integrating sphere, energy balance, closure enforcement | validation, code verification | "radiosity"; [Integrating sphere](@ref) |
| Beam power at any defocus; defocus second moments | code verification | "beam profile"; [Laser spot on a foil](@ref) |
| Fast deposition against the ray trace, with shadowing | solution verification | "deposition vs ray trace"; [Laser spot on a foil](@ref) |
| Legendre and Y_lm recovery; facet integration | code verification | "mode analysis", "spherical harmonics and 3D symmetry" |
| NIF D2 selection rule; superposition; influence matrix | code verification | "spherical harmonics and 3D symmetry"; [Dropped outer-cone bundle](@ref) |
| Gradients against finite differences | code verification | "hohlraum + gradients", "CBET", [Dropped outer-cone bundle](@ref) |
| CBET coupling against an independent Im Γ computation | code verification | "CBET" |
| CBET transfer against published points | calibration | `calibrate_cbet` |
| Clearance envelopes, wings, worst-case formulas | code verification | "clearance" |

## Not yet established

There is no comparison with VisRad calculations yet. The planned benchmark kit exports
published configurations to run in VisRad and compares the outputs, which will validate
absolute drive levels and mode amplitudes rather than only their internal consistency.
