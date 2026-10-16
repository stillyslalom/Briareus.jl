# Briareus.jl

View-factor radiation and laser-pointing design for high-energy-density experiments:
indirect-drive hohlraums, backlighter foils, and other laser-irradiated targets. Briareus
covers the core of what VisRad (MacFarlane, JQSRT 81, 287, 2003) is used for — laser
deposition, 3D view-factor radiosity, capsule drive symmetry, synthetic images — and adds
gradient-based optimization of pointing, defocus and power balance.

```julia
using Briareus

target = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
quads = nif_quads()                                       # 48 NIF quads, official port table
pointing = GroupPointing(quads; fields = (:z, :r, :power))   # per-cone parameters
scenario = Scenario(target, quads; observe = :capsule, parameterization = pointing,
                    viewfactor_kw = (; closed = true))

x0 = [3.09, -1.09, 8.0,  3.63, -0.60, 8.0,  5.20, 0.0, 8.0,  5.01, 0.0, 8.0]
report(evaluate(scenario, x0))                            # power accounting + P_l/P0

result = optimize_pointing(scenario, symmetry_objective(modes = (2 => 1.0, 4 => 1.0)); x0)
```

See `examples/nif_hohlraum.jl` and `examples/backlighter.jl`.

## What it computes

- **Geometry.** Targets assembled from meshed primitives (cylinder, annulus/disk, sphere,
  rectangle, frustum) with per-component materials. `hohlraum()` and `foil()` build common
  targets. Surfaces are one- or two-sided; laser entrance holes are transparent
  *apertures* that measure the escaping x-ray flux.
- **Laser deposition.** Beams enter from facility ports with elliptical super-Gaussian
  phase-plate spots, an f-number cone, pointing, and defocus. Two methods:
  - `DepositionOperator` + `deposit` — precomputed quadrature with per-direction
    visibility masks. It is smooth in pointing, defocus, spot size and power, so
    ForwardDiff gradients flow through it. It takes about 3 ms for 48 quads on about
    3,000 elements.
  - `raytrace_deposition` — an energy-conserving geometric ray trace, used as the
    reference.

  Power that hits the outside of a target (LEH clipping), misses it, or strikes the
  capsule directly is accounted for separately.
- **Radiation.** Occlusion-aware view factors `view_factors` (the exact
  point-to-polygon kernel, so elements sharing an edge stay accurate) and gray diffuse
  radiosity `solve`. Wall albedo can be constant, tabulated in time, or the Marshak
  wall-loss scaling `WallLossAlbedo`, iterated to self-consistency with the local
  radiation temperature. For closed enclosures, `closed = true` removes discretization
  leakage while keeping reciprocity. `flux_response` precomputes the linear map from
  laser deposition to capsule flux, so optimizers never re-solve.
- **Analysis.** Legendre moments (`legendre_moments`, with exact facet solid-angle
  integration) and real spherical harmonics (`ylm_moments`, `mode_power`).
  `energy_balance` covers source, wall/capsule absorption, LEH loss and leakage.
  `radiation_temperature` converts flux to temperature.
- **Imaging.** `render` produces orthographic images of any per-element field (x-ray
  radiance, laser intensity) from any direction or port. `line_of_sight_power` gives
  the x-ray power toward a distant detector.
- **Laser–plasma corrections.** Applied between the laser setup and deposition, and
  differentiable in power and Δλ:
  - per-cone, time-dependent `Backscatter`;
  - an empirical `CalibratedTransfer` (transfer as a function of Δλ);
  - a reduced kinetic `CBET` model at the LEHs: multi-species Maxwellian ion-acoustic
    response, flow Doppler shift, polarization smoothing, and a smooth δn/n saturation
    clamp. `calibrate_cbet` fits its multiplier, flow and clamp to transfer data.
- **Shot-review checks.**
  - `clearance` gives the worst-case beam-to-hardware clearance, power and intensity
    over a linear error stack (pointing, target offset, hardware tolerance, focus). It
    uses intensity-contour envelopes and clearance-only spot wings (`nif_wings`), with
    pass/review/fail tiers.
  - `unconverted_light_check` covers the 1ω/2ω footprints, the 32 mm zone, and the
    curvature rule.
- **Beam-table I/O.** `write_beam_table`/`read_beam_table` read and write a canonical
  CSV, with column and unit mapping for foreign tables. Element and moment exports are
  also available.
- **Time dependence.** `PulseShape` per cone and `drive_history` give quasi-static
  T_r(t) and P_l(t).
- **Optimization.** `GroupPointing` maps parameters to beams, per cone, quad or beam:
  axial/radial/tangential aim, defocus, power and spot scale, with top/bottom
  mirroring. `symmetry_objective` targets capsule drive symmetry and
  `irradiance_objective` targets foil uniformity. `optimize_pointing` runs bounded
  L-BFGS-B with forward-mode AD gradients.

## Conventions

Units are mm, ns, TW (TW·ns = kJ), TW/mm² (`to_Wcm2` multiplies by 1e14) and eV.
Port angles are in degrees (θ from +z, φ counterclockwise from +x).

A beam's direction is fixed by its port (`k = -port`), and pointing translates its axis
through `aim`. Best focus is at `aim + defocus·k`. Spot axes are given in the port's
(θ̂, φ̂) frame.

Away from focus, the spot is the focal profile convolved with the f-cone. The
deposition operator samples that convolution at the cone directions: each surface point
receives, from each cone direction that reaches it unobstructed, the focal profile where
that ray crosses the best-focus plane, weighted by the ray's own incidence angle. Power
is conserved at any distance from focus and incidence, and shadows near focus are
correlated with the part of the spot they block.

## Facility data and provenance

| Data | Source | Notes |
|---|---|---|
| NIF quad ports (θ, φ) | NIF User Guide §4.1, Table 4-1 | Port 63 is misprinted as 62 there. |
| NIF per-beam offsets within a quad | Derived from the published 63.2 × 55.7 cm beam spacing and 7.7 m focal length | Per-beam tables are not public. |
| NIF CPP spots | 30°/44.5°/50°: Kemp et al., PoP 23, 101204 (2016); 23.5°: User Guide FWHM | Minor axis along θ̂, per Kemp; this is not officially documented. |
| OMEGA beams | Truncated icosahedron with edge ratio 1.2 in the hohlraum (P6–P7) frame; reproduces the 21.415/42.020/58.852/81.249° cones | The azimuthal origin is not tied to chamber coordinates. Load the real port table with `beams_from_table`. |
| OMEGA phase plates | OMEGA Users' Guide (2014) | SG4, SG8, E-IDI-300 |

## Provenance of physics defaults

`reference/research-notes.md` records sources and confidence levels for:
- the CBET formulas and benchmark;
- NIF Δλ conventions (quoted at 1ω since about 2013; `Beam.dlambda` is at the beam's
  own wavelength);
- backscatter practice;
- unconverted light;
- error budgets and clearance defaults.

Clearance thresholds and wing models are engineering defaults, not facility criteria.
Override them with your facility's review requirements and measured far fields.

## Scope and limitations

Briareus is a geometric view-factor and ray-trace code with reduced laser–plasma
corrections. It has no hydrodynamics, refraction, wall motion or LEH closure.

The CBET model assumes uniform LEH plasma and full pairwise overlap. Its unsaturated
coupling overpredicts transfer, so calibrate it (`calibrate_cbet`) against shot data
or reference transfer points before trusting magnitudes.

Radiosity is gray, with no spectral resolution. Visibility between elements uses a
single centroid ray.

## Tests

```bash
julia --project=. -t auto -e 'using Pkg; Pkg.test()'
```

The suite checks the code against analytic results:
- coaxial-disk and sphere view factors, and occlusion;
- the integrating-sphere radiosity solution and energy conservation;
- beam power normalization at any defocus;
- analytic-versus-ray-trace deposition, including shadowing near focus and hohlraum
  capsule symmetry;
- AD-versus-finite-difference gradients;
- Legendre and Y_lm recovery, and facility geometry checks.
