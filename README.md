# Briareus.jl

[![CI](https://github.com/stillyslalom/Briareus.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/stillyslalom/Briareus.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://stillyslalom.github.io/Briareus.jl/dev/)

View-factor radiation and laser-pointing design for high-energy-density experiments:
indirect-drive hohlraums, backlighter foils, and other laser-irradiated targets. Briareus
covers the core of what VisRad (MacFarlane, JQSRT 81, 287, 2003) is used for — laser
deposition, 3D view-factor radiosity, capsule drive symmetry, synthetic images — and adds
gradient-based optimization of pointing, defocus and power balance, worst-case beam
clearance and target-configuration checks, and file exchange with VisRad.

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

See `examples/nif_hohlraum.jl` and `examples/backlighter.jl`, and the documentation,
whose tutorials double as validation against exact results
(`julia --project=docs docs/make.jl`). `BriareusGUI/` is a desktop window (GLMakie + QML)
for the same workflow: target, beams, drive, optimization and shot review.

## What it computes

- **Geometry.** Targets assembled from meshed primitives (cylinder, annulus/disk, sphere,
  rectangle, frustum) with per-component materials. `hohlraum()` and `foil()` build common
  targets. Surfaces are one- or two-sided; laser entrance holes are transparent
  *apertures* that measure the escaping x-ray flux. A component can carry the exact
  surface its mesh approximates (`SphereSurface`, `CylinderSurface`, `DiskSurface`), and
  the radiation model then integrates over the true curved surface; hohlraums do this by
  default. `hohlraum(windows = [WallWindow(...)])` cuts diagnostic windows into the wall:
  closed ones re-emit with their own material, open ones are apertures that laser light
  passes through and x rays escape from, and the mesh conforms to their edges.
- **Laser deposition.** Beams enter from facility ports with elliptical super-Gaussian
  phase-plate spots, an f-number cone, pointing, and defocus. Two methods:
  - `DepositionOperator` + `deposit` — precomputed quadrature with per-direction
    visibility masks, refined where each beam can land. It is smooth in pointing,
    defocus, spot size and power, so ForwardDiff gradients flow through it; 48 quads on
    a 3,000-element hohlraum take about 30 ms, and each beam's power is conserved to
    1e-4.
  - `raytrace_deposition` — an energy-conserving geometric ray trace, used as the
    reference.

  Power that hits the outside of a target (LEH clipping), misses it, or strikes the
  capsule directly is accounted for separately.
- **Radiation.** Occlusion-aware view factors `view_factors` (the exact point-to-polygon
  kernel near the emitter, adaptive toward shared edges; about 0.4 s for a
  3,000-element hohlraum) and gray diffuse radiosity `solve`. Wall albedo can be
  constant, tabulated in time, or the Marshak wall-loss scaling `WallLossAlbedo`,
  iterated to self-consistency with the local radiation temperature. For closed
  enclosures, `closed = true` removes discretization leakage while keeping reciprocity.
  `flux_response` precomputes the linear map from laser deposition to capsule flux, so
  optimizers never re-solve.
- **Analysis.** Legendre moments (`legendre_moments`) and real spherical harmonics
  (`ylm_moments`, `flux_ylm`), both with exact facet solid-angle integration;
  `mode_amplitude` gives the fractional amplitude of each mode (the mode-1 drive
  asymmetry, for example). `influence_matrix` gives each beam's mode response.
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
- **Shot-review checks.** All are judged on the worst case over pointing, target
  positioning, hardware tolerance and focus errors, never on the nominal beams.
  - `clearance` gives the worst-case beam-to-hardware clearance, power and intensity,
    with intensity-contour envelopes, clearance-only spot wings (`nif_wings`) and
    pass/review/fail tiers.
  - `talis_check` applies NIF's published target-configuration rules: 3ω capture and
    incidence, 3ω and 1ω clearance, unconverted-light management (curvature, dimpling
    or the specular 20° rule), and the TAS stay-in zone.
  - `unconverted_light_check` covers the 1ω/2ω footprints and the 32 mm zone.
  - `ApertureKeepOut` puts the worst-case clearance of every beam from its entrance
    hole into an optimization objective.
- **I/O.** `write_beam_table`/`read_beam_table` read and write a canonical CSV, with
  column and unit mapping for foreign tables. VisRad interchange: NIF SST pointing and
  Laser Beam Summary CSVs, power files, and a benchmark kit (`write_visrad_kit`,
  `compare_visrad`) for comparing the same case run in both codes.
- **Time dependence.** `PulseShape` per cone and `drive_history` give quasi-static
  T_r(t) and P_l(t).
- **Optimization.** `GroupPointing` maps parameters to beams, per cone, quad or beam:
  axial/radial/tangential aim, defocus, power, spot scale and Δλ, with top/bottom
  mirroring. `symmetry_objective` and `asymmetry_objective` target capsule drive
  symmetry (axisymmetric and three-dimensional), and `irradiance_objective` targets foil
  uniformity. `optimize_pointing` runs bounded L-BFGS-B with forward-mode AD gradients,
  with a callback for progress and cancellation.

## Conventions

Units are mm, ns, TW (TW·ns = kJ), TW/mm² (`to_Wcm2` multiplies by 1e14) and eV.
Port angles are in degrees (θ from +z, φ counterclockwise from +x).

A beam's direction is fixed by its port (`k = -port`), and pointing translates its axis
through `aim`. Best focus is at `aim + defocus·k`; positive defocus moves the focus
along the beam, as in VisRad and the NIF CMT. Spot axes are given in the port's (θ̂, φ̂)
frame.

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
| NIF quad names | LMF chamber-geometry figures, confirmed by the LMF beam table | `nif_bundle`, `nif_bundles` |
| NIF per-beam directions | LMF Shot RI Resources, beam port locations table (2020) | `nif_beams()`, `data/nif_beams.csv` |
| NIF CPP spots | LMF Shot RI Resources, measured phase-plate fits (2020) | Converted from HWHM to 1/e semi-axes. |
| NIF target-configuration rules | LMF "TCC Configuration Requirements" (2018) | How the 0.5° tolerance is applied is not stated; see `talis_check`. |
| OMEGA beams | Truncated icosahedron with edge ratio 1.2 in the hohlraum (P6–P7) frame; reproduces the 21.415/42.020/58.852/81.249° cones | The azimuthal origin is not tied to chamber coordinates. Load the real port table with `beams_from_table`. |
| OMEGA phase plates | OMEGA Users' Guide (2014) | SG4, SG8, E-IDI-300 |
| VisRad formats | VisRad manual (prism-cs.com, v21.1.1) | Items to verify against real files are listed in the research notes. |

## Provenance of physics defaults

`reference/research-notes.md` records sources and confidence levels for:
- the CBET formulas and benchmark;
- NIF Δλ conventions (quoted at 1ω since about 2013; `Beam.dlambda` is at the beam's
  own wavelength);
- backscatter practice;
- unconverted light;
- error budgets, clearance defaults and the TaLIS rules;
- the VisRad file formats.

The `clearance` thresholds and wing models are engineering defaults, not facility
criteria. Override them with your facility's review requirements and measured far
fields.

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
- coaxial-disk, sphere and closed-cylinder view factors, and occlusion;
- the integrating-sphere radiosity solution and energy conservation;
- beam power normalization at any defocus, and per-beam power conservation of the
  deposition operator on fine and coarse meshes, with and without clipping;
- analytic-versus-ray-trace deposition, including shadowing near focus and grazing
  incidence far from focus;
- AD-versus-finite-difference gradients;
- Legendre and Y_lm recovery, the NIF port symmetry selection rule, and facility
  geometry checks;
- clearance, keep-out and TaLIS rules against analytic geometry, and VisRad file
  round trips.

The documentation's tutorials add exact references: the f-cone convolution of a laser
spot, the oblique-incidence footprint, the integrating sphere, and the symmetry of the
NIF ports.
