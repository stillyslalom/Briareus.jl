# Verification and validation

The tutorials and the test suite carry out four different activities, and their
results mean different things.

- **Code verification** checks whether the implementation solves the intended discrete
  equations. Examples are exact view factors of catalog geometries, power conservation
  of the ray trace, symmetry selection rules that must hold to round-off, and
  derivatives checked against finite differences.
- **Solution verification** measures how large the discretization error of a given
  calculation is. Examples are the closure error of a view-factor matrix and its
  convergence with mesh refinement, and the agreement of the fast deposition operator
  with the ray trace.
- **Validation** tests whether the model reproduces an independent reference, either
  physical or analytic, for the intended use. Examples are the integrating-sphere
  multiplier, the exact convolution of a laser spot with its f-cone, and, eventually,
  comparisons with VisRad calculations and shot data.
- **Calibration** fits model parameters to reference data. The cross-beam energy
  transfer model is calibrated, not validated: its multiplier, flow and saturation
  clamp are fitted to measured transfer.

The closing section of each tutorial, "What this checks", lists which of these it carries
out.

## Where each check lives

| Check | Activity | Location |
|:--|:--|:--|
| BVH ray casting against brute force | code verification | `test/runtests.jl`, "ray casting" |
| Coaxial-disk view factor, occlusion, point-to-polygon kernel | code verification | "view factors" |
| Integrating sphere, energy balance, closure enforcement | validation, code verification | "radiosity"; [Integrating sphere](@ref) |
| LEH escape and capsule fractions against a Monte Carlo photon walk over LEH size (N210808 hohlraum); peak-power saving 3.64 → 3.1 mm against N210307/N210808 | code verification, approximate validation | `test/leh_power_balance.jl`; [Entrance-hole size and the power balance](@ref) |
| Beam power at any defocus; defocus second moments; oblique-incidence correction | code verification | "beam profile"; [Laser spot on a foil](@ref) |
| Defocused footprint against the exact f-cone convolution | validation | [Laser spot on a foil](@ref) |
| Fast deposition against the ray trace, with shadowing and at grazing incidence far from focus | solution verification | "deposition vs ray trace"; [Laser spot on a foil](@ref); [Freeing a bundle for radiography](@ref) |
| Legendre and Y_lm recovery; facet integration | code verification | "mode analysis", "spherical harmonics and 3D symmetry" |
| NIF D2 selection rule; superposition; influence matrix | code verification | "spherical harmonics and 3D symmetry"; [Freeing a bundle for radiography](@ref) |
| Gradients against finite differences | code verification | "hohlraum + gradients", "CBET" |
| Optimizer against the closed-form least-squares solution | code verification | [Freeing a bundle for radiography](@ref) |
| CBET coupling against an independent Im Γ computation | code verification | "CBET" |
| CBET transfer against published points | calibration | `calibrate_cbet` |
| Two-zone CBET against the N221204 HYDRA transfer (PRE 109 Fig. 9, digitized): cone fraction over the pulse, per-cone outer ratios not fitted | calibration | `test/cbet_zones.jl`; `examples/n221204_cbet_calibration.jl` |
| Clearance envelopes, wings, worst-case formulas | code verification | "clearance" |
| Differentiable entrance-hole keep-out against the clearance check | code verification | "clearance"; [Freeing a bundle for radiography](@ref) |

## Not yet established

There is no comparison with VisRad calculations yet. The planned benchmark kit exports
published configurations to run in VisRad and compares the outputs, which will validate
absolute drive levels and mode amplitudes rather than only their internal consistency.
