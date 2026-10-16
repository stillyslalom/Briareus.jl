# CLAUDE.md

Briareus.jl — view-factor radiation and laser-pointing design (a VisRad replacement/
complement) for indirect-drive hohlraums and backlighter foils. See [README.md](README.md)
for scope and [ROADMAP.md](ROADMAP.md) for the backlog.

## Commands

```bash
julia --project=. -t auto -e 'using Pkg; Pkg.test()'    # full suite (~40 s)
julia --project=. -t auto examples/nif_hohlraum.jl       # NIF symmetry optimization
julia --project=. -t auto examples/backlighter.jl        # OMEGA backlighter uniformity
```

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
  of the emitter, far pairs from both centroids, averaged for reciprocity. `closed = true`
  (Sinkhorn balancing) only for closed enclosures.
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
`history.jl` (PulseShape, drive_history).
