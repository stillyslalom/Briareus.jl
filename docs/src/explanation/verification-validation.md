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
  multiplier, the exact convolution of a laser spot with its f-cone, the comparison of
  laser deposition with VisRad calculations, and the N221204 benchmark against published
  shot data.
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
| Hard-sphere direct-drive illumination against Skupsky–Lee separability, Temporal et al. 2014 (one ring; LMJ), Murakami and Nishi 2017 | validation, solution verification | `test/direct_drive.jl` |
| NIF polar-direct-drive pointing (Barlow et al. 2024's Optimized Config) against an independent plane-wave quadrature; σ₁ compared with Barlow's, not fitted | solution verification, approximate validation | [Polar direct drive: a published NIF pointing](@ref) |
| Clearance envelopes, wings, worst-case formulas | code verification | "clearance" |
| Differentiable entrance-hole keep-out against the clearance check | code verification | "clearance"; [Freeing a bundle for radiography](@ref) |
| Laser deposition against VisRad 21 on NIF hohlraum workspaces | calibration (beam profile), validation (deposition) | research notes, "VisRad deposition parity: hohlraum runs" |
| N221204 drive history against PRE 109 Fig. 10 and Dante (blind) | approximate validation | `examples/n221204.jl`; research notes, "N221204 time-dependent benchmark" and later sections |

## Comparison with VisRad

VisRad 21 computed the laser deposition on workspaces built from Briareus geometry: an
N210808-like hohlraum (gold wall of radius 3.2 mm and length 11.24 mm, open LEHs, 2388
wall triangles of about 0.6 mm), its capsule, two window patches, a foil tilted 20°
and LEH shields, with 192 NIF beams and VisRad's built-in phase plates. One workspace
is a reference pointing, the other a re-pointing of it. Briareus read each workspace,
built the same target and beams, and compared its own deposition per element with the
direct laser source VisRad wrote into the file.

The beam profile away from focus is the one piece calibrated to VisRad.
[`VisRadProfile`](@ref) with `measured = true` tabulates the growth of the spot's semi-axes
and exponent with distance from focus, measured on a separate VisRad workspace that
deposits one beam per NIF quad on flat plates up to 15 mm from best focus. On the
hohlraum workspaces, which played no part in those tables:

- wall deposition differs from VisRad's by 0.9 % rms per element on both workspaces,
  with peaks of 779.1 against 778.1 TW/cm² and 911.5 against 915.9 TW/cm², and a
  correlation of 1.0000 on the reference pointing;
- total wall power is 0.07 % below VisRad's, window power 0.17 % below, and the tilted
  foil's 0.37 % above;
- absolute capsule P2 and P4 from the two codes' depositions, through the same
  Briareus radiosity, agree to 0.014 and 0.004 percentage points, and the re-pointing's
  change of P2 is 0.104 % against VisRad's 0.092 %.

The workspaces and comparison scripts are not public; the numbers, the setup and what
remains open (the law behind the tabulated semi-axes for plates other than the four
NIF ones, the per-beam renormalization on tilted surfaces) are in
`reference/research-notes.md`. [`write_visrad_kit`](@ref) and [`compare_visrad`](@ref)
set up the same kind of comparison for any case.

## The N221204 benchmark

`examples/n221204.jl` runs the time-dependent drive of NIF shot N221204 (Kritcher et
al., Phys. Rev. E 109, 025204 (2024)) on a public mock-up of the HYBRID-E pointing, with
the two-zone CBET model calibrated to that paper's Fig. 9, and the wall motion, LEH
closure and inner-beam bubble calibrated to other published data. Neither the paper's
Fig. 10 (the capsule flux P2/P0 against time) nor the measured Dante history is a fit
target. On the default settings (coarse mesh, measured backscatter), the P2 swing reaches
−0.181 at 5.25 ns and turns positive at 6.19 ns, against −0.170 at 5.6 ns and 6.67 ns
in Fig. 10; the synthetic Dante 1 at 8 ns reads 290.9 and 314.2 eV over the initial and
open LEH areas, against the measured 310.1 eV. After the laser falls, the synthetic Dante
reads 37–54 eV low at 9 ns on three shots: the quasi-static radiosity has no plasma
energy store. The research notes record each run and the inputs that remain
unconstrained.
