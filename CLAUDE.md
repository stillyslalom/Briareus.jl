# CLAUDE.md

Briareus.jl — view-factor radiation and laser-pointing design (a VisRad replacement/
complement) for indirect-drive hohlraums and backlighter foils. See [README.md](README.md)
for scope and [ROADMAP.md](ROADMAP.md) for the backlog.

## Commands

```bash
julia --project=. -t auto -e 'using Pkg; Pkg.test()'    # full suite (~2 min)
julia --project=. -t auto examples/nif_hohlraum.jl       # NIF symmetry optimization
julia --project=. -t auto examples/backlighter.jl        # OMEGA backlighter uniformity
julia --project=docs -t auto docs/make.jl                # docs; runs every tutorial
julia --project=docs -t auto docs/literate/<name>.jl     # one tutorial as a script
```

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
  for reciprocity. Components carrying a `CurvedSurface` (Sphere/Cylinder/DiskSurface)
  are integrated on the true surface — hohlraums do this by default; ray casting always
  uses the facets. `closed = true` (Sinkhorn balancing) only for closed enclosures.
  Visibility ignores facets of a cylinder/sphere when both ray ends are inside it
  (faceting artifacts), and the centroid facing test is relaxed by each curved
  element's normal spread (a general triangle on a cylinder is tilted ~h/2R).
- Patches must be planar: on a cylinder, mesh anything but chord rectangles as
  triangles (a folded patch's centroid sits behind its facets and sees nothing).
  Circular `WallWindow`s use equal-area edge polygons for exact areas.
- Thin-foil windows (`BurnthroughAlbedo`): `drive_history` tracks each element's heat
  front (`HeatFront`, Hammer & Rosen 2003) from its own T(t). A single `solve` treats
  the foil as intact unless given the model's `burnthrough` or a per-element
  `burnthrough` vector. `sol.albedo` is the albedo B was solved with, and `sol.intact`
  its before-burnthrough part, which `energy_balance` uses to split out `transmitted`.
- Deposition (operator and ray trace) sums the best-focus profile over `ncone` (default
  13) f-cone directions, each with its own incidence factor `|n·dⱼ|/(k·dⱼ)`. Don't use the
  moment-matched `intensity` for deposition: it has the right second moment but the
  wrong shape far from focus (2% excess power, 10% misplaced, on hohlraum walls). It
  remains for envelopes.
- Validate against independent references (exact convolutions, closed forms, ray
  traces converged at nspot ≥ 121), never against the method under test; the ray trace
  aliases against patch grids at coarse nspot.
- Beam clearance is judged on the worst case over pointing, positioning and focus
  errors (`clearance`); `ApertureKeepOut` puts the same worst case into objectives.
- Facility data provenance is recorded in README; derived numbers (NIF per-beam offsets,
  OMEGA hohlraum frame) are labeled as such in docstrings.

## Architecture (src/, include order)

`geometry.jl` (Mesh, primitives) → `bvh.jl` (BVH, closest/any hit) → `target.jl`
(materials, albedo models, Component, Target, element quadrature) → `targets.jl`
(hohlraum, foil) → `beams.jl` (SuperGaussianSpot, Beam, BeamFrame, intensity, cone
directions) → `facilities.jl` (NIF ports/spots, OMEGA) → `deposition.jl`
(DepositionOperator, deposit, raytrace_deposition) → `viewfactors.jl` → `radiation.jl`
(RadiationModel, solve, flux_response, energy_balance, enforce_closure!) → `analysis.jl`
(solid angles, Legendre/Y_lm) → `imaging.jl` (OrthoCamera, render, line_of_sight_power)
→ `optimize.jl` (GroupPointing, Scenario, objectives, optimize_pointing, report) →
`history.jl` (PulseShape, drive_history) → `clearance.jl` (envelopes, wings, worst-case
clearance) → `keepout.jl` (differentiable entrance-hole keep-out) → `lpi.jl` (LaserPlasmaModel, Backscatter, CalibratedTransfer) → `io.jl`
(CSV beam tables) → `cbet.jl` (reduced kinetic CBET, calibrate_cbet) → `unconverted.jl`.

Physics defaults and their sources live in `reference/research-notes.md`. Keep each
number's provenance and confidence there, and say in docstrings when a value is
derived or an engineering default rather than a facility criterion.

Editing pitfall on this Windows box: Python `open()` defaults to cp1252. Use
`encoding='utf-8'`, or the Edit tool, for files with Unicode. Julia triple-quoted
strings dedent, so multi-line `replace` patterns built from them miss indented code.
