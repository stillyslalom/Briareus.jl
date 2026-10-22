# CLAUDE.md

Briareus.jl — view-factor radiation and laser-pointing design (a VisRad replacement/
complement) for indirect-drive hohlraums and backlighter foils. See [README.md](README.md)
for scope and [ROADMAP.md](ROADMAP.md) for the backlog.

## Commands

```bash
julia --project=. -t auto -e 'using Pkg; Pkg.test()'    # full suite (~5 min on 4 threads)
julia --project=. -t auto examples/nif_hohlraum.jl       # NIF symmetry optimization
julia --project=. -t auto examples/backlighter.jl        # OMEGA backlighter uniformity
julia --project=docs -t auto docs/make.jl                # docs; runs every tutorial
julia --project=docs -t auto docs/literate/<name>.jl     # one tutorial as a script
```

The package precompiles a workload (end of `src/Briareus.jl`); `BRIAREUS_NO_WORKLOAD=true`
skips it for quick rebuilds. Precompile caches do not key on that variable: after a
build without the workload, touch a source file before measuring time to first result.

Docs follow CompactLES.jl: Literate tutorials in `docs/literate/`, registered in
`docs/make.jl`, each ending in "What this checks" with exact references. The docs
environment has no ForwardDiff/Optim, and `Vec3` clashes with Makie there.

## Conventions

- Units: mm, ns, TW, TW/mm² (×1e14 → W/cm²), eV; angles in degrees at the API.
- Beam direction is fixed by the port (`k = -port`); pointing translates the axis through
  `aim`; best focus at `aim + defocus·k`; spot `a` axis at `angle` from port θ̂ toward φ̂.
  `DepositionOperator` is valid only while ports are unchanged.
- Elements = radiating patch sides. One-sided surfaces radiate from the front; their back
  side absorbs laser light as a per-component loss (`backside`, reported `<name>_outside`).
  Apertures are transparent to rays/visibility and have α = η = 0.
- Everything on the optimization path (`apply`, `deposit`, `flux_response` product,
  `flux_moments`, objectives) must stay generic in the number type for ForwardDiff;
  threaded loops go through a function barrier (a runtime-computed eltype captured in a
  `Threads.@threads` closure boxes and allocates per point).
- Mode analysis integrates Pₗ / Yₗₘ over each facet's exact solid angle (field constant
  per element); don't revert to centroid × A cosθ/R² — it leaks P0 into P2 at the
  0.5 % level on default capsule meshes.
- View factors: near pairs use the exact point-to-polygon kernel over sub-triangle points
  of the emitter (adaptive toward shared edges), mid/far pairs fixed quadrature, averaged
  for reciprocity. Components carrying a `CurvedSurface` (Sphere/Cylinder/Cone/DiskSurface)
  are integrated on the true surface — hohlraums do this by default; ray casting always
  uses the facets. `closed = true` (Sinkhorn balancing) only for closed enclosures; it is
  the `RadiationModel` default for targets declared `closed` (`hohlraum` does), since row-sum
  errors of ~1e-3 grow ~13× at wall albedo 0.925 (0.3% drive, ~0.05% P2 between meshes).
  Elements of more than two facets (remeshed CAD) emit from lumped cells and receive
  through their outline (or plane by plane when cut by the point's tangent plane);
  1–2-facet elements (every primitive) keep the per-facet path bit for bit.
  Visibility ignores facets of a cylinder/sphere/cone when both ray ends are inside it
  (faceting artifacts), and the centroid facing test is relaxed by each curved
  element's normal spread (a general triangle on a cylinder is tilted ~h/2R).
- Patches must be planar: on a cylinder, mesh anything but chord rectangles as
  triangles (a folded patch's centroid sits behind its facets and sees nothing).
  Circular `WallWindow`s use equal-area edge polygons for exact areas. The exception is
  remeshed CAD (`remesh_mesh`: elements folded across facets within `angle`): `Target`
  moves such an element's centroid along its normal onto its facets (a planar patch
  keeps its centroid bit for bit), or on a convex surface (a capsule) it would see
  nothing.
- CAD curved walls: `fit_surfaces` finds planes, cylinders, cones and spheres in a mesh;
  `curved_components` (or `visrad_target(...; remesh, fit = true)`) splits a part into
  `<part>/cyl1`, …, `<part>/cone1`, …, `<part>/sph1`, …, `<part>/rest`, each carrying
  its fitted surface. `ConeSurface` projects radially at constant height (rim chords
  land on the rim circles, so faceted frusta tile exactly); elements touching a cone's
  apex stay flat (the surface is singular there).
  Planes get none (flat elements are exact; `DiskSurface` assumes polar sectors). A
  curved element of more than two facets gets a lumped quadrature (as many points as a
  two-triangle element), so view-factor cost does not grow with the CAD faceting;
  elements of one or two facets (every primitive) keep the per-facet rule.
- Thin-foil windows (`BurnthroughAlbedo`): `drive_history` tracks each element's heat
  front (`HeatFront`, Hammer & Rosen 2003) from its own T(t). A single `solve` treats
  the foil as intact unless given the model's `burnthrough` or a per-element
  `burnthrough` vector. Once burned through, a foil transmits the fraction
  `transmission` T (default `eddington_transmission(optical_depth(m))`, Lester et al.
  2026: 4/(3τ+2), τ = 5 at 2.5 µm Au, linear in areal mass) and re-emits the rest with
  `before`: α = (1 − sT)α_before. An explicit `after` (e.g. 0, the open-hole bound)
  replaces this and books only α_before − α as transmitted. `sol.albedo` is the albedo
  B was solved with, `sol.intact` its before-burnthrough part, and `sol.transmission`
  the transmitted fraction that `energy_balance` splits out as `transmitted`.
  `RosselandTransmission` makes `T` follow the radiation temperature (iterated in
  `solve`). `GapClosure` (an open gap that closes at a given time) also transmits;
  albedo models report (albedo, intact, transmitted) through `_albedo_split`, and
  `_transmits` marks the ones `energy_balance` and `drive_history` book.
- Deposition (operator and ray trace) sums the best-focus profile over `ncone` (default
  13) f-cone directions, each with its own incidence factor `|n·dⱼ|/(k·dⱼ)`. Don't use the
  moment-matched `intensity` for deposition: it has the right second moment but the
  wrong shape far from focus (2% excess power, 10% misplaced, on hohlraum walls). It
  remains for envelopes. `profile = VisRadProfile(...)` (from `visrad_profile(ws)`)
  swaps in VisRad's model so exported designs reproduce in VisRad: a super-Gaussian
  cut at each beam's envelope ratio ρ (default (ln 1000)^(1/n)) and normalized per
  transverse plane, along the beam axis only. For the four NIF plates a(z), b(z), n(z)
  are VisRad 21's, measured on plates normal to 48 beams (`_VR_MEASURED`, research
  notes "VisRad beam-profile probe"; 0.9% wall rms on the hohlraum runs); other spots
  keep hyperbolic semi-axes with z_A × γ/ρ (γ = `growth` = 0.67, fitted before the
  probe). Spot rotation in VisRad files is relative to the beam's NIF quad horizontal
  (`_vr_reference_angle`, measured to 0.02°), not the beam's own φ̂. Don't refit the
  tables to hohlraum runs; remeasure with thor `visrad/build_profile_probe.jl`.
- VisRad's "Fraction of Total Beam Power" red/green model makes the power column the
  total; the 3ω cone carries 1 − f_G − f_R of it even when green/red are off.
- `DepositionOperator` samples unmodeled back sides (clipped power) only for their total:
  one point per triangle away from the beam, graded Gauss rules toward open edges
  (LEH rims) and their vertices, `resolution` within a spot semi-axis of them,
  `backside_resolution` elsewhere. Clipped light at 1e-4 lies within ~1/50 of a spot of
  the rim; uniform sub-triangles underestimated it by 10–60%. Element (front) points
  never depend on the back-side keywords. Check clipped power at 1e-4 against ray
  traces averaged over nspot 241–361: one nspot alone aliases by ±2% there.
- Validate against independent references (exact convolutions, closed forms, ray
  traces converged at nspot ≥ 121), never against the method under test; the ray trace
  aliases against patch grids at coarse nspot.
- Beam clearance is judged on the worst case over pointing, positioning and focus
  errors (`clearance`); `ApertureKeepOut` (LEH rims) and `SphereKeepOut` (the capsule)
  put the same worst case into objectives. The `direct` deposition penalty is
  quadratic in a tiny fraction and does not keep beams off the capsule.
- `optimize_pointing` bounds both iterations and objective evaluations (a line search
  on a kinked objective can take dozens per iteration). It hands L-BFGS-B value and
  gradient from one dual pass, under one fixed dual tag (`_PointingTag`), so compiled
  and precompiled code serves every objective.
- Performance of the optimization path: a `Scenario` keeps the flux response contracted
  with the cached Legendre/Y_lm weights (`moment_response`, `ylm_response`); the full
  `scenario.response` and `ev.flux` are computed lazily on first read, so test for
  radiation with the scenario's `model`/`observed`, never `ev.flux === nothing`. Dual
  mat-vecs go through BLAS (`_mul`). `deposit` skips points and cone directions off the
  spot on Float64 values (exact zeros); keep such screens exact. Keep heterogeneous
  collections (keep-out tuples, wing lists) out of generators and abstract-typed loops
  on hot paths: a `sum` over a generator on a mixed tuple took inference 27 s.
- Mesh sizes: `hohlraum(; h, capsule_h)`; the GUI presets come from
  `examples/mesh_convergence.jl` (research notes): the capsule mesh sets P4 (error ∝
  capsule_h²), fine (0.17/0.05 mm) is converged to ~0.02 percentage points.
  `capsule_mesh = :compact` (geodesic `icosphere_mesh`) reaches that P4 with ~1k
  capsule elements instead of ~8k; lat-long stays the default until the user decides.
- Facility data provenance is recorded in README; derived numbers (NIF per-beam offsets,
  OMEGA hohlraum frame) are labeled as such in docstrings.

## Architecture (src/, include order)

`geometry.jl` (Mesh, primitives) → `meshops.jl` (orient/coarsen imported meshes) → `bvh.jl` (BVH, closest/any hit) → `target.jl`
(materials, albedo models, Component, Target, element quadrature) → `targets.jl`
(hohlraum, foil) → `cadfit.jl` (fit_surfaces, curved_components) → `beams.jl` (SuperGaussianSpot, Beam, BeamFrame, intensity, cone
directions) → `facilities.jl` (NIF ports/spots, OMEGA) → `deposition.jl`
(DepositionOperator, deposit, raytrace_deposition) → `viewfactors.jl` → `enclosure.jl`
(enclosure: the elements bounding the free space around a point) → `radiation.jl`
(RadiationModel, solve, flux_response, energy_balance, enforce_closure!) → `analysis.jl`
(solid angles, Legendre/Y_lm) → `imaging.jl` (OrthoCamera, render, line_of_sight_power)
→ `optimize.jl` (GroupPointing, Scenario, objectives, optimize_pointing, report) →
`history.jl` (PulseShape, drive_history) → `retune.jl` (pulse tables, retune_drive) →
`clearance.jl` (envelopes, wings, worst-case clearance) → `keepout.jl` (differentiable
entrance-hole keep-out) → `uncertainty.jl` (beam-error sensitivities, robust objectives)
→ `talis.jl` → `lpi.jl` (LaserPlasmaModel, Backscatter, CalibratedTransfer) → `io.jl`
(CSV beam tables, write_results) → `visrad.jl`, `visrad_workspace.jl` (VisRad files and
workspaces, read and write) → `cbet.jl` (reduced kinetic CBET, calibrate_cbet) →
`unconverted.jl`.

Physics defaults and their sources live in `reference/research-notes.md`. Keep each
number's provenance and confidence there, and say in docstrings when a value is
derived or an engineering default rather than a facility criterion.

Editing pitfall on this Windows box: Python `open()` defaults to cp1252 and, in text
mode, writes CRLF. Use `encoding='utf-8', newline='
'`, or the Edit tool. Julia triple-quoted
strings dedent, so multi-line `replace` patterns built from them miss indented code.
