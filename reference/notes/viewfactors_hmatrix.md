# Research notes: View factors, meshes, visibility and the H-matrix

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

## Hohlraum mesh work-precision (2026-10-16) *derived*

`examples/mesh_convergence.jl`: the GUI's NIF hohlraum (10.1 × 2.875 mm, LEH 1.55,
capsule 0.9, wall albedo 0.85, conversion 0.8), NIF quads at the example cone aims
("nominal") and at an optimizer-like pointing with the inner cones pulled to the waist
("tight"). Elements are sized with `hohlraum(; h, capsule_h)`. Build = deposition
operator + view factors + flux response, 16 threads.

Wall/cap element size `h` (capsule 0.07 mm), Legendre moments of the capsule flux in
percentage points:

| h (mm) | N | build | P2 nominal | P2 tight | P4 tight |
|---|---|---|---|---|---|
| 0.80 | 3,900 | 3.6 s | −14.64 | −4.54 | −1.127 |
| 0.60 | 4,276 | 2.4 s | −14.56 | −4.29 | −1.168 |
| 0.45 | 5,028 | 3.3 s | −14.60 | −4.14 | −1.196 |
| 0.35 | 5,796 | 3.9 s | −14.57 | −4.06 | −1.203 |
| 0.28 | 7,260 | 5.4 s | −14.54 | −3.97 | −1.205 |
| 0.22 | 9,428 | 12 s | −14.54 | −3.94 | −1.212 |
| 0.17 | 13,076 | 24 s | −14.52 | −3.90 | −1.216 |
| 0.13 | 19,572 | 58 s | −14.53 | −3.90 | −1.216 |

Capsule element size (wall 0.2 mm):

| capsule h (mm) | N | P2 nominal | P4 nominal | P4 tight |
|---|---|---|---|---|
| 0.30 | 7,208 | −13.94 | +1.238 | −0.687 |
| 0.22 | 7,372 | −14.20 | +0.814 | −0.909 |
| 0.16 | 7,656 | −14.38 | +0.508 | −1.057 |
| 0.12 | 8,160 | −14.49 | +0.325 | −1.167 |
| 0.09 | 9,056 | −14.52 | +0.245 | −1.205 |
| 0.07 | 10,452 | −14.53 | +0.216 | −1.213 |
| 0.05 | 13,620 | −14.55 | +0.185 | −1.216 |
| 0.04 | 17,232 | −14.55 | +0.177 | −1.226 |

Findings:
- The wall converges as about h^2.5–3 in P2 (tight: 0.64, 0.25, 0.075, 0.005 pp at
  h = 0.8, 0.45, 0.28, 0.17 against 0.13). The nominal design wanders by ±0.03 pp
  without a trend below h = 0.45.
- The capsule mesh dominates P4: a piecewise-constant flux on facets converges as
  capsule_h² (P4 nominal extrapolates to ≈ +0.16%; error 0.08 pp at 0.09 mm, 0.6 pp
  at 0.2 mm). The old "standard" capsule (24 × 32, ≈ 0.12 × 0.18 mm) was 0.15 pp off.
- The mean drive and laser power fractions converge to 1e-4 already at h = 0.45.

GUI presets chosen from these (errors against the finest meshes / extrapolation):
`:coarse` h = 0.6, capsule 0.2 (≈ 1,300 elements, < 1 s; P2 ≈ 0.4 pp, P4 ≈ 0.5 pp);
`:standard` h = 0.28, capsule 0.09 (≈ 5,900, ≈ 4 s; ≲ 0.1 pp); `:fine` h = 0.17,
capsule 0.05 (≈ 16,000, ≈ 40 s; ≲ 0.02 pp, converged for design work). Build times
are from before the 2026-10-17 performance pass; a fine `Scenario` (16,244 elements,
closed) now takes 15 s on 24 threads (view factors 4 s, LU 7 s), and an optimizer
iteration about 60 ms.

## Compact capsule mesh (2026-10-22) *derived*

`examples/capsule_mesh.jl`: the hohlraum and the nominal and tight pointings of the
mesh work-precision study above, wall fixed at h = 0.2 mm, closed, 16 threads. The
capsule is either the latitude-longitude `sphere_mesh` (`hohlraum(capsule_mesh =
:latlong)`, the default, sized by `capsule_counts`) or the geodesic `icosphere_mesh`
(`:compact`: icosahedron faces divided into n × n triangles projected onto the sphere,
20n² elements, edges within −16/+19% of 1.107R/n, n = ⌈1.107R/capsule_h⌉). Both carry
the `SphereSurface`. N is the capsule's element count; moments in percentage points.

| mesh | capsule h | N | P2 nominal | P4 nominal | P4 tight | build |
|---|---|---|---|---|---|---|
| lat-long | 0.30 | 200 | −13.938 | +1.238 | −0.687 | 12 s |
| lat-long | 0.16 | 648 | −14.383 | +0.508 | −1.057 | 9 s |
| lat-long | 0.12 | 1,152 | −14.488 | +0.325 | −1.167 | 10 s |
| lat-long | 0.09 | 2,048 | −14.524 | +0.245 | −1.205 | 11 s |
| lat-long | 0.07 | 3,444 | −14.533 | +0.216 | −1.213 | 15 s |
| lat-long | 0.05 | 6,612 | −14.546 | +0.185 | −1.219 | 22 s |
| lat-long | 0.04 | 10,224 | −14.549 | +0.177 | −1.226 | 35 s |
| compact | 0.30 | 320 | −14.386 | +0.203 | −1.149 | 9 s |
| compact | 0.22 | 500 | −14.454 | +0.190 | −1.185 | 9 s |
| compact | 0.16 | 980 | −14.499 | +0.177 | −1.211 | 9 s |
| compact | 0.12 | 1,620 | −14.518 | +0.171 | −1.221 | 10 s |
| compact | 0.09 | 2,880 | −14.532 | +0.167 | −1.228 | 12 s |
| compact | 0.07 | 4,500 | −14.539 | +0.165 | −1.231 | 15 s |
| compact | 0.05 | 8,000 | −14.544 | +0.164 | −1.234 | 26 s |
| compact | 0.04 | 12,500 | −14.547 | +0.163 | −1.235 | 45 s |

Findings:
- Both converge as c/N (∝ capsule_h²) to the same limits, extrapolated from each mesh's
  two finest points: P4 nominal +0.161 (lat-long) / +0.162 (compact), P4 tight −1.237 /
  −1.237, P2 nominal −14.553 / −14.551. The mean drive is within 2e-4 of its limit at every size on both.
- P4 error ≈ 160/N pp (nominal) and 120/N (tight) on the lat-long capsule, 16/N and
  24/N on the compact one: within 0.02 pp of the limit takes ~8,000 lat-long elements
  (capsule_h ≈ 0.045) against ~800 (nominal) to ~1,200 (tight) compact ones
  (capsule_h ≈ 0.12–0.16). The compact mesh is 5–10× cheaper for P4, consistent with the
  remeshed CAD capsule's 408 elements (another cavity, 2026-10-18 entry above).
- P2 converges at the same rate on both (≈ 50/N pp nominal; 0.02 pp at ~2,500
  elements): the capsule's piecewise-constant sampling, not element shape, sets it.
  The tight design's P2 is within 0.01 pp at every compact size.
- Odd moments stay at round-off: the icosphere is centrally symmetric. (Table above: a
  vertex on +z. Since it became the `hohlraum()` default, 2026-10-22, the icosahedron
  sits in its standard frame, mirror-symmetric in x, y and z, because the 5-fold axis on
  z leaked odd-m modes of ~6e-5 of the drive into the D2-symmetric NIF hohlraum; see
  "Default capsule mesh and visibility".)
- Element quality, not count, sets P4: the lat-long cells near the poles are slivers
  whose centroids sample P4 poorly; at 0.09 mm the "standard" lat-long capsule is still
  0.08 pp off in P4, the compact one 0.005.
- Build time is set by the wall at these sizes (~9 s up to ~3k capsule elements).

## Partial visibility (2026-10-22) *derived*

`view_factors(...; visibility = :sampled, visibility_grid = g)` (also through
`RadiationModel(t; visibility = :sampled)` and a `Scenario`'s `viewfactor_kw`). Each
element's area is lumped into a g × g grid of its tangent plane (area-weighted points); up
to four of those points (the extremes along ±e₁, ±e₂) are probes. Per facing pair, after the
centroid ray: probe rays from each element's probes to the other's (lifted) centroid. When
all agree with the centroid ray, the factor is exactly the `:centroid` factor; otherwise the
g² × g² sample rays give a visible fraction weighted by cosθᵢ cosθⱼ / r² that scales the
unobstructed factor (computed even when the centroid ray is blocked). With one fraction per
unordered pair, reciprocity is exact. Rays use the facets and the curved-surface skip rule
of the centroid ray. Default `:centroid` checked equal (`==`) to main's `viewfactors.jl` on
`hohlraum(h = 1, capsule_h = 0.4)` and the bored cavity at 1 mm.

The geometry is the primitive bored cavity of the CAD tutorial (R 2.875, L 9.4, 0.1 mm LEH
bores of radius 1.55 with `nphi ÷ 2` facets, 20 µm waist gap, capsule 0.9, apertures on
the LEHs and the gap exit), `closed = false`, 8 threads. Largest |row sum − 1| (rms):

| component | 0.5 mm, centroid | 0.5 mm, g = 4 | 0.25 mm, centroid | 0.25 mm, g = 4 | 0.25 mm, g = 6 |
|---|---|---|---|---|---|
| walls | 9.7e-3 (4.1e-3) | 1.9e-3 (5.0e-4) | 9.1e-3 (3.3e-3) | 1.4e-3 (6.9e-4) | 1.2e-3 (6.0e-4) |
| caps | 1.7e-3 (1.2e-3) | 1.4e-3 (7.8e-4) | 2.2e-3 (1.2e-3) | 1.2e-3 (5.4e-4) | 9.7e-4 (4.5e-4) |
| bores | 5.7e-3 | 2.6e-3 | 5.2e-3 | 2.1e-3 | 6.3e-4 |
| capsule | 1.6e-3 (9.2e-4) | 5.0e-4 (2.4e-4) | 2.0e-3 (1.1e-3) | 2.0e-3 (6.6e-4) | 2.0e-3 (6.6e-4) |
| LEH apertures | 5.9e-2 (3.2e-2) | 6.7e-3 (2.8e-3) | 6.1e-2 (3.5e-2) | 2.0e-2 (9.3e-3) | 1.5e-2 (7.4e-3) |
| gap faces | 5.0e-2 | 5.5e-2 | 5.2e-2 | 5.5e-2 | 3.3e-2 |
| gap exit | 2.3e-2 | 2.9e-2 | 3.9e-3 | 6.5e-3 | 6.1e-3 |
| time | 1.9 s | 4.0 s (g = 2: 2.5 s) | 4.5 s | 18.5 s (g = 3: 13.4 s) | 50 s |

(2,244 and 6,304 elements.) The 50 µm gap faces behind a 20 µm slot are not helped:
the slot is narrower than a wall sample cell. `hohlraum()` (3,264 elements, flush LEHs):
1.3 s → 4.6 s; at h = 0.5 its walls' 0.30% goes to 0.28% (most changed pairs there are
capsule shadows).
At 1 mm with 12-facet bores the option makes LEH rows worse (3.9% → 7.7%): the faceted
bore (sagitta 0.05 mm on a 0.1 mm bore) over-shadows rays from the LEH rim, which the
fractions then reproduce faithfully. With 64-facet bores at 0.5 mm the residual on the
walls (4.5e-3, near the waist) is the 0.3 mm capsule's faceting.

The reference is the visible fractions of 56 changed pairs at 0.5 mm (the 30 of largest
exchange change and 26 at random), against kernel-weighted visibility over 256 × 256 points
per facet pair (same rays, no grid): centroid error mean 0.42, max 0.49; sampled (g = 4)
mean 0.017, max 0.08. At 1 mm, one LEH element against every wall element: 0.020 / 0.16.
Half-shadowed pair (unit squares 2 apart, one turned 45°, a plate over half the mid-plane;
exactly ½ by point symmetry): 0.493 at g = 4, 0.500 at g = 8 (g = 3: 0.430, its center
samples tie on the edge).

Effect after Sinkhorn balancing (the tutorial's two-ring wall source): capsule flux
+0.05%, a₂/a₀ −0.175 pp (0.5 mm) and −0.155 pp (0.25 mm), a₄/a₀ +0.11 pp at both, LEH
loss −0.4%. The shift barely shrinks with the mesh, so balancing does not absorb the
centroid rule's visibility error; which limit is right needs a finer mesh or a
ray-traced radiosity reference. That reference is in the next entry, and the sampled
rule is closer to it. **M**

## Monte Carlo radiosity reference for partial visibility (2026-10-22) *derived*

> **Superseded in part (2026-10-25).** The comparisons stand, but the last bullet's
> reading does not: the flat facets of both geometries are cracked (24-sided LEH disks
> in 20-sided bores and on caps of `nphi` sides, a 180-sided gap aperture on 40-sided
> gap faces; `open_fraction` 7.5e-4 and 2.4e-4), so rays escape that the curved geometry
> keeps, and balancing them adds the exchange the cracks lose. On watertight facets and
> on the curved path balancing does not bias a₂/a₀ upward ("Closure on bored, windowed
> and CAD cavities (2026-10-25)").

`montecarlo_incident` (src/montecarlo.jl): diffuse transport on the
facets with no view factors. Paths start at a uniform point on a capsule element's
facets and go backward (reciprocity of the Lambertian kernel): cosine-weighted
direction, first hit by `closest_hit` over all triangles (apertures included: a hit
there is an escape), the hit element's source `ηP/A` added with the path weight, weight
× albedo, Russian roulette below 1/20. `reemit = :hit` continues from the hit point
(the continuum); `:element` from a uniform point on the hit element (the element model
with exact view factors and visibility). Statistical errors from 64 batches. Checked
(test/montecarlo.jl): a closed cube with uniform source gives S/(1 − α) within 1%; an
open-top cube against the closed-form factors (Hottel 0.199825 opposite,
(1 − 0.199825)/4 adjacent) within 5σ and 1%; 8 × 8 elements per face against `solve`
within 1.5% on face means (one element per face; the continuum differs by up to 26%).

Since the Monte Carlo traces facets, the comparison uses flat elements (the same meshes
with `surface = nothing`); the curved models are a different geometry (true surfaces)
and are listed only for the size of the visibility shift. Source: the CAD tutorial's
two Gaussian wall rings, integrated over each element, η = 0.75, wall α = 0.8, capsule
α = 0.3. 200k paths per capsule element (105–310 s on 16 threads).
`examples/mc_radiosity.jl`. Capsule incident flux, a₀ (TW/mm²), a₂/a₀ and a₄/a₀ (%):

| case | a₀ | a₂/a₀ | a₄/a₀ |
|---|---|---|---|
| **bored 0.5 mm (2,244 el.)**: MC `:element` | 0.85004 ± 0.00004 | 10.641 ± 0.012 | −0.510 ± 0.015 |
| MC `:hit` (continuum) | 0.85009 ± 0.00005 | 10.667 ± 0.011 | −0.530 ± 0.016 |
| flat, centroid, raw | +0.45% | +0.025 | −0.063 |
| flat, centroid, balanced | +0.08% | +0.508 | −0.274 |
| flat, sampled, raw | +0.29% | +0.108 | −0.063 |
| flat, sampled, balanced | +0.12% | +0.312 | −0.130 |
| **bored 0.25 mm (6,304 el.)**: MC `:element` | 0.84811 ± 0.00003 | 10.752 ± 0.008 | −0.600 ± 0.011 |
| MC `:hit` | 0.84835 ± 0.00003 | 10.746 ± 0.007 | −0.587 ± 0.011 |
| flat, centroid, raw | +0.02% | +0.070 | +0.019 |
| flat, centroid, balanced | +0.16% | +0.303 | −0.107 |
| flat, sampled, raw | +0.11% | +0.033 | +0.016 |
| flat, sampled, balanced | +0.21% | +0.142 | +0.007 |
| **`hohlraum()` (3,264 el.)**: MC `:element` | 0.84677 ± 0.00005 | 10.228 ± 0.012 | −0.465 ± 0.014 |
| MC `:hit` | 0.84673 ± 0.00004 | 10.238 ± 0.010 | −0.458 ± 0.014 |
| flat, centroid, raw | +0.46% | −0.091 | +0.016 |
| flat, centroid, balanced | +0.22% | +0.086 | +0.013 |
| flat, sampled, raw | +0.18% | +0.019 | +0.011 |
| flat, sampled, balanced | +0.23% | +0.083 | +0.018 |

(Model rows: a₀ relative to MC `:element`, moments as pp differences.) Curved models,
balanced, sampled − centroid: a₂/a₀ −0.175 / −0.155 pp, a₄/a₀ +0.11 / +0.11 pp (bored
0.5 / 0.25 mm), −0.003 / +0.004 pp on `hohlraum()`; flat: −0.196 / −0.161 and
+0.144 / +0.114 pp, −0.003 / −0.005 pp. The shift is the same on both geometries.

- The two Monte Carlo modes agree within 0.03 pp (2σ): the element discretization of the
  radiosity is negligible here, and MC `:element` is the reference for the view factors.
- After balancing, sampled visibility is closer in every bored case: a₂/a₀ error
  0.51 → 0.31 pp (0.5 mm) and 0.30 → 0.14 pp (0.25 mm), a₄/a₀ 0.27 → 0.13 and
  0.11 → 0.01 pp. The −0.16 pp shift is a correction, not an artifact; it is about
  half the balanced centroid model's error. On the flush-LEH control both visibilities
  balance to the same answer (+0.085 pp), as expected.
- Balancing itself biases a₂/a₀ upward: +0.06 to +0.48 pp over the raw factors of the
  same visibility, on every case including the control. Raw sampled factors are the
  closest of the four in a₂/a₀ on the 0.25 mm cavity and the control (0.02–0.03 pp);
  their a₀ is off 0.1–0.3%, balanced 0.1–0.25%. Presumably Sinkhorn spreads the
  residual row error (on the cavity, the waist gap's faces at 5% and the LEH apertures
  at 1–2%) over each row's partners symmetrically, which is not where the missing
  exchange is (not checked element by element). This holds for flat
  elements; the curved (default) path's row errors differ and a curved-surface Monte
  Carlo would be needed to check balancing there (done: next entry; on the curved path
  balancing moves toward the continuum). **M**

## Balancing on the curved path against a continuum Monte Carlo (2026-10-22) *derived*

Does `closed = true` (the `hohlraum()` default) bias the capsule drive on the default
curved path? Reference: `montecarlo_incident(...; reemit = :hit)` (the continuum) on a
flat copy of the hohlraum faceted at wall/cap h = 0.1 mm, capsule_h = 0.05 mm (33,292
elements, 6,612 on the capsule; the chord sagitta of the wall is 4e-4 mm). Geometry
10.1 × 5.75 mm, 3.1 mm LEHs, 0.9 mm capsule; wall gold α = 0.925 and 0.85, η = 0.8.
Source: 48 NIF quads at 8 TW at the example cone aims (`examples/mesh_convergence.jl`
"nominal"), deposited on each mesh by a separate `DepositionOperator` (totals agree to
4e-5 between meshes) and scaled to the reference total. The departure turns one outer-44 quad
off; by linearity the Monte Carlo traces only the source of that quad (paired batches give
its error). 298M paths at 0.925 (8.6 min on 16 threads), 198M at 0.85.
`examples/mc_closure.jl`.

Reference convergence in faceting (α = 0.925, 32M paths each): wall h 0.2 / 0.1 / 0.06
gives a₀ 8.4148 / 8.4116 / 8.4118 (± 0.0011) and a₂/a₀ −12.547 / −12.546 / −12.554
(± 0.027)%; capsule_h 0.035 instead of 0.05 moves a₀ +0.03 ± 0.02% and a₂/a₀
+0.04 ± 0.04 pp (a₄/a₀ +0.08 ± 0.05 pp, mostly the piecewise-constant moment analysis,
which `capsule_h` sets for model and reference alike). Converged to ~0.03% in a₀ and
~0.03 pp in a₂/a₀.

Capsule incident flux: reference a₀ (TW/mm²), a₂/a₀, a₄/a₀ (%), and the departure (design −
quad off, with Δa₀ in % and Δa₂/a₀ and Δa₄/a₀ in pp); model rows are differences from the
reference (a₀ in %, the rest in pp). Default visibility (centroid). Curved = the
`hohlraum()` default; flat = the same meshes without surfaces, a coarser polygonal geometry
than the reference (so its rows mix geometry error in).

| α = 0.925 | a₀ | a₂/a₀ | a₄/a₀ | Δa₀ | Δa₂/a₀ | Δa₄/a₀ |
|---|---|---|---|---|---|---|
| **MC `:hit`** | 8.4118 ± 0.0003 | −12.513 ± 0.010 | +0.225 ± 0.010 | +2.045 ± 0.001 | +0.233 ± 0.001 | −0.035 ± 0.001 |
| curved 0.45/0.15, raw | +0.56% | +0.161 | +0.266 | +0.000 | −0.006 | −0.003 |
| curved 0.45/0.15, balanced | −0.07% | +0.091 | +0.249 | +0.000 | −0.004 | −0.004 |
| curved 0.28/0.09, raw | −0.03% | +0.088 | +0.149 | −0.001 | −0.001 | −0.002 |
| curved 0.28/0.09, balanced | −0.05% | +0.029 | +0.064 | −0.000 | −0.001 | −0.002 |
| curved 0.17/0.05, raw | −0.09% | +0.036 | +0.006 | +0.000 | +0.001 | −0.001 |
| curved 0.17/0.05, balanced | −0.03% | +0.020 | +0.007 | +0.000 | +0.000 | −0.001 |
| flat 0.45/0.15, raw / balanced | +1.12 / +0.54% | +0.030 / +0.038 | +0.075 / +0.063 | +0.001 | +0.000 / +0.001 | +0.003 / +0.002 |
| flat 0.28/0.09, raw / balanced | +0.24 / +0.22% | −0.017 / +0.030 | +0.007 / −0.000 | −0.001 / +0.000 | +0.001 | +0.001 |
| flat 0.17/0.05, raw / balanced | +0.03 / +0.07% | +0.039 / +0.025 | −0.014 / −0.015 | +0.000 | +0.001 | −0.000 |

| α = 0.85 | a₀ | a₂/a₀ | a₄/a₀ | Δa₀ | Δa₂/a₀ | Δa₄/a₀ |
|---|---|---|---|---|---|---|
| **MC `:hit`** | 6.1316 ± 0.0003 | −14.547 ± 0.011 | +0.150 ± 0.014 | +2.013 ± 0.001 | +0.315 ± 0.001 | −0.046 ± 0.002 |
| curved 0.45/0.15, raw | +0.35% | +0.148 | +0.322 | −0.001 | −0.009 | −0.007 |
| curved 0.45/0.15, balanced | −0.06% | +0.077 | +0.306 | −0.001 | −0.008 | −0.007 |
| curved 0.28/0.09, raw | −0.03% | +0.074 | +0.178 | −0.002 | −0.003 | −0.005 |
| curved 0.28/0.09, balanced | −0.04% | +0.016 | +0.095 | −0.001 | −0.003 | −0.005 |
| curved 0.17/0.05, raw | −0.07% | +0.022 | +0.034 | −0.001 | −0.001 | −0.004 |
| curved 0.17/0.05, balanced | −0.02% | +0.010 | +0.035 | −0.001 | −0.001 | −0.004 |

Mean |row sum − 1| before balancing (curved): 7.6e-4 / 4.6e-4 / 2.0e-4, largest
3.2e-3 / 3.1e-3 / 1.5e-3 at the three meshes.

- On the curved path balancing does not bias the drive: it moves a₂/a₀ toward the
  continuum at every mesh and both albedos (0.161 → 0.091, 0.088 → 0.029,
  0.036 → 0.020 pp at 0.925) and fixes the study mesh's a₀ (+0.56% → −0.07%; the raw
  excess is the ~8e-4 row-sum error amplified by the wall albedo). The flat-element
  upward bias of the previous entry does not carry over: raw curved factors already
  sit above the reference in a₂/a₀, and balancing pulls them down.
- Balanced curved models are within 0.07% in a₀ and 0.09 pp in a₂/a₀ at every
  preset, 0.03% and 0.02 pp at the fine one (inside the ~0.03
  convergence of the reference). a₄/a₀ at the coarser presets is dominated by the capsule
  mesh (the moment analysis on 0.15 and 0.09 mm facets, cf. "Mesh convergence"); at
  0.17/0.05, whose capsule mesh is the reference's, it agrees to 0.01 pp at 0.925.
- Departures (one quad off) agree with the reference within 0.01 pp in Δa₂/a₀ and
  Δa₄/a₀ and 0.002% in Δa₀ for every model, raw or balanced: balancing is immaterial
  for departures here.
- Bored LEHs, gaps and CAD cavities are not covered (larger row errors; the flat bored
  cavity of the previous entry is the only data there), sampled visibility, windows.
  **H** for `hohlraum()` at these meshes and albedos; **M** beyond.

> **Covered since (2026-10-25):** bored LEHs, a waist gap, windows with gaps and a
> remeshed CAD cavity, centroid and sampled visibility, in "Closure on bored, windowed
> and CAD cavities (2026-10-25)".

## Default capsule mesh and visibility (2026-10-22) *decision*

User decision (Open question 9): `hohlraum()` builds its capsule as the geodesic
`icosphere_mesh` (`capsule_mesh = :compact`); `:latlong` reproduces earlier baselines.
`Target` carries a `visibility` rule that `view_factors` uses by default: `:sampled` for
`cad_target` and `visrad_target` (rims), `:centroid` for `hohlraum()` and hand-assembled
targets.

- **Orientation.** With a vertex on +z (the table in "Compact capsule mesh") the
  icosahedron's 5-fold axis lies along the hohlraum axis, and the D2-symmetric NIF
  hohlraum picked up odd-m capsule modes of ~6e-5 of the drive (test/runtests.jl and
  test/windows.jl selection rules, threshold 1e-10). The icosahedron now sits in its
  standard frame (vertices at the cyclic permutations of (0, ±1, ±ϕ)): its 2-fold axes
  and mirror planes are x, y and z, so the selection rules hold within round-off. P2/P4
  convergence rates are a property of element shape and are unchanged in kind; the
  table's compact rows were not rerun in the new frame.
- **GUI presets.** capsule_h 0.3 / 0.16 / 0.09 mm (320 / 980 / 2,880 capsule elements):
  P2 needs ~2,500 capsule elements for 0.02 pp on any capsule mesh (≈ 50/N pp), so the
  fine preset is not 0.12 as first proposed.
- **Effect on a study.** A backlighter study (N210808 + two diagnostic windows, h 0.45 /
  capsule_h 0.15, closed, tuned NIF beams) moved by ≤ 0.01 pp in every departure and
  drive between the two capsule meshes; a design optimized on the lat-long capsule
  evaluates at δ2 = 0.015% on the compact one.

## Large radiation models (2026-10-24) *derived*

ROADMAP "Scale". The goal is to make 1e4–5e4-element radiation models (the fine hohlraum
preset, remeshed CAD) practical on a laptop (Ryzen 9 4900HS, 4 threads of 8 used, 16 GB,
another job sharing the machine). Workload, as a design session runs it: the GUI NIF
hohlraum (10.1 × 2.875 mm, LEH 1.55, capsule 0.9) with a `WallLossAlbedo` wall, NIF quads;
deposition operator, view factors, closure, a temperature-dependent `solve` at 2 ns, a
6-step `drive_history`, a `Scenario` (62-row moment response), an optimizer gradient and
`ev.flux`. Script: the session's scratchpad `prof.jl` (one process per mesh,
`BRIAREUS_NO_WORKLOAD=true`, so first calls include compilation: "solve again" is the warm
figure).

### Where the time and memory went (before)

| N (h / capsule_h) | 3,712 (0.6/0.09) | 10,484 (0.22/0.07) |
|---|---|---|
| deposition operator | 2.4 s | 3.5 s |
| view factors | 1.2 s | 14.3 s |
| closure | 0.3 s | 2.1 s |
| solve (first, LU) / again | 0.9 / 0.4 s | 9.1 / 3.0 s |
| drive_history, 6 steps | 4.3 s | 20.9 s |
| Scenario (incl. a copy of F) | 3.4 s | 13.2 s |
| gradient | 0.02 s | 0.02 s |
| ev.flux (full response) | 0.9 s | 8.9 s |
| max RSS | 1.4 GB | 6.6 GB |

Dense Float64 holds F and up to two LU factors, 24N² bytes: 7.5 GB at 17.6k, 20 GB at
29k, 38 GB at 40k; one LU is 2N³/3 flops (≈40 s at 17.6k on 4 threads, ≈7 min at 40k).
The fine preset did not fit next to the other job, and 40k does not fit at all.

### What changed

1. **Float32 storage** (`Float32Matrix`, `view_factors(...; storage = Float32)`): Float32
   data read as Float64; every product accumulates in Float64, so it is exactly the
   product with `Float64.(F.data)`. It is assembled directly in Float32 (no Float64 peak)
   and is identical to `Float32.(view_factors(t))`. Balancing (`enforce_closure!`) runs in
   Float64 on the stored values and rounds once: row sums to ~1.5e-8 (Float64 ~1e-12),
   negligible beside the 1e-3 discretization error it removes.
2. **Factorization-free solves** (`RadiationModel(...; solver = :iterative)`): restarted
   block GMRES (`_gmres_block`), no preconditioner, 1e-12 relative residual like the
   LU-preconditioned path; the columns of a block (dual value and partials, spectral
   groups, the 62 moment rows) iterate in lockstep, each with a separate basis, sharing one
   product with F per iteration; bases kept under 256 MB (columns batched). Iterations
   at 17.6k: 17 for a cold solve, fewer warm (a history step starts from the last
   radiosity), 28 for the transposed moment rows. `ev.flux` on such a model is one solve
   instead of the full observed × N response (still available, a solve per row).
3. **Defaults**: above `Briareus.LARGE_MODEL[]` = 10,000 elements `RadiationModel`
   computes Float32 view factors and solves iteratively; below it nothing changes (dense
   Float64 + LU, with identical results). `storage`/`solver` override; a given `F` keeps its
   storage. `Scenario(...; model)` reuses a model (no copy of F).
4. **Kernels**: vector products threaded over row/column blocks at 28–33 GB/s (memory
   bound); a fused 4-column kernel for blocks of ≤ 48 columns (3–7× the widen-to-Float64
   + BLAS panels at k = 2–16; the panels win past ~50 columns). 17.6k: F x 43 ms,
   Fᵀx 37 ms, k = 2: 64/50 ms, k = 16: 280/240 ms, k = 62: 880/890 ms (panels).
5. **View-factor assembly** (both storages, identical results): the near-pair kernel's result
   was inferred `Any` (the CAD cell list is untyped), so every pair boxed its factor and
   dispatched the store dynamically: 13.7M allocations at 5k, and GC pressure across
   threads. A type assertion removed them: 10.5k 14.5 → 12.1 s (allocation 2.9 → 0.7
   GB), 17.6k 56 → 31 s (9.1 → 1.5 GB).

### After (Float32 + iterative; same script)

| N | 10,484 | 17,632 | 29,124 |
|---|---|---|---|
| deposition operator | 6.5 s | 6.6 s | 9.8 s |
| view factors | 13.0 s | 31.4 s | 101.6 s (140 s before the fix in 5) |
| closure | 1.6 s | 3.0 s | 6.2 s |
| solve again (T-dep albedo, warm) | 1.0 s | 2.5 s | 7.4 s |
| drive_history, 6 steps | 8.9 s | 14.7 s | 31.4 s |
| Scenario (62-row response) | 11.2 s | 17.2 s | 39.0 s |
| gradient | 0.02 s | 0.03 s | 0.04 s |
| ev.flux (one solve) | 0.37 s | 0.64 s | 1.4 s |
| max RSS | 1.3 GB | 2.3 GB | 4.5–5.3 GB |

(The before runs loaded a package image with the precompile workload, the after runs one
without it, so first calls here include compilation: ~3 s of the deposition operator,
~8 s of the first solve.) 40k extrapolates
to 6.4 GB for F, ~190 s of view factors, ~14 s per warm solve.

Accuracy of the Float32 path against the dense Float64 path at 10,484 elements (same
mesh, same process otherwise): capsule a₀ 1.3e-10 relative, P2/P0 and P4/P0 2e-10 and
1e-10 absolute (solve and Scenario), drive history moments ≤ 2e-10, Tr ≤ 2e-9 eV, H
1.2e-8 relative at worst (an LEH element). The Float64 iterative path agrees with LU to
the GMRES residual (~1e-12). Tests: `test/large_models.jl` (forced on at ~1k elements:
products, block GMRES, solve/flux_response/duals/drive_history/Scenario gradient/
`spectral_solve`/keyframe interpolation against the dense path).

### Tried and dropped

- **Two-grid preconditioner** (clusters of 16 or 4 elements by coordinate bisection,
  Galerkin coarse view factors, coarse LU per albedo): 17 → 12 (11) iterations per solve,
  28 → 26 (24) for the moment rows, no gain in time. The slow part of the spectrum is not
  the smooth one the clusters capture.
- **Row-blocked stores in the assembly** (rows written 16 at a time): the profile's
  samples on the stores were skid from the divisions; no change in time.

### Not done: GPU and hierarchical radiosity

- **GPU** (RTX 2060 Max-Q, 6 GB, ~1/32-rate FP64): Float32 F fits to ~35k elements; a
  mat-vec with Float64 accumulation needs 2 flops per 4-byte entry, ~130 GFLOP/s at the
  card's ~260 GB/s, about its FP64 rate, so ~6–8× the CPU's 30 GB/s is the ceiling. The
  solves are ~30% of the workload at 17–29k (Scenario response, history, solve), so a
  CUDA extension would save ≲ 25% end to end; view-factor assembly (BVH rays and exact
  contour kernels, CPU) remains the largest cost. Not worth CUDA.jl yet.
- **Hierarchical (H-matrix) view factors** are the real answer past ~30k: assembly and
  storage are N², and both grow 4× per halving of h. The design is a cluster tree by
  coordinate bisection (as tried above). Admissible cluster pairs (distance > η × diameter,
  and visibility tested by a few rays per pair and found all-or-nothing) are stored as
  low-rank blocks (ACA on the far-pair kernel, rank 1–8), and the rest (near pairs, partly
  occluded pairs straddling the capsule shadow or LEH rims) are stored dense. Memory and
  assembly are O(N log N). It needs only `mul!` with F and Fᵀ and row sums: the iterative
  path here already uses nothing else, so an H-matrix type plugs in as another
  `AbstractMatrix{Float64}` (balancing is a diagonal scaling s F s). Reciprocity per block
  (G = AF symmetric) halves it again. The next section, "Hierarchical view factors",
  describes the build.

## Hierarchical view factors (2026-10-24) *derived*

ROADMAP "Scale", the step after "Large radiation models": view factors whose assembly and
storage do not grow as N², as another matrix type (`HMatrix`, `src/hmatrix.jl`) that the
iterative path uses through products only. `view_factors(t; storage = :hmatrix)`,
`RadiationModel(t; storage = :hmatrix)`, and `RadiationModel`'s default above
`Briareus.HMATRIX_MODEL[]` = 15,000 elements. Workload and machine as in "Large radiation
models" (GUI NIF hohlraum 10.1 × 2.875 mm, LEH 1.55, capsule 0.9, `WallLossAlbedo` wall,
NIF quads aimed at z = ±4; 4 threads, another job sharing the machine); script: the
session's scratchpad `hmat/bench.jl` (warmed by a small target through the same path first).

### Design

- **What is stored.** The exchange matrix G = A F (symmetric), once per unordered pair:
  a cluster tree by median bisection of the centroids along the longest box axis (leaves
  ≤ 48 elements), blocks of the upper triangle only (Gᵀ blocks implied), `F = A⁻¹ S G S`
  with a balancing diagonal S. Reciprocity is exact; `enforce_closure!` (Sinkhorn, one
  product per iteration) updates S only. Dense blocks (near pairs, leaf pairs that do
  not compress) are Float32 with the entries `view_factors` computes (`_vf_pair`, now
  split from the row loop; the dense path is unchanged). Low-rank blocks are
  U Vᵀ (Float32 factors), and "scattered" dense blocks on index lists hold the rows a
  facing test cuts.
- **Admissibility.** Cluster boxes (centroids padded by element size) at least η × the
  larger diameter apart (η = 1). Then, per block:
  1. *Facing* (box test, no rays): rows whose element faces no corner of the other
     cluster's centroid box (beyond the facing test's slack) are zero; a block with
     none is zero. Small blocks (≤ 16 leaf²) take the facing test pair by pair: a block
     facing itself in under a quarter of its rows × columns (a capsule cluster at the
     terminator seen by a wall cluster: 93–100% zeros, ACA finds zero rows and stops)
     is stored dense on those.
  2. *Shaft culling* (Haines & Wallace 1994): the convex hull of the two boxes (the union
     box's faces plus the planes through one axis-parallel edge of each box) is searched
     in the BVH for facets that could occlude a segment between the clusters. Not
     occluders: facets of a curved surface (cylinder, sphere, cone) both clusters lie
     inside (the rule of `_vf_visible`); facets of a convex closed component with
     outward normals (a capsule; checked vertex against facet, ≤ 2e9 pairs) when one
     cluster lies on it with one facet per element (a segment leaving a face of a convex
     polyhedron outward does not come back); facets whose plane has both clusters
     strictly in front (by 2ε) or in that plane with its normal (caps and LEH discs,
     corner clusters with the wall in front of the cap). An empty shaft makes every
     facing pair wholly visible, so its entries skip the ray (`_vf_pair_free`: pairs
     facing beyond the slack; slack-admitted pairs still take theirs).
  3. *Hidden by a convex part*: if a convex component's inscribed ball (about the vertex
     mean) blocks every segment between the boxes — each corner of one box inside the
     ball's shadow cone from each corner of the other, beyond that corner's polar plane;
     the shadow beyond the polar plane is convex, which carries the corner test to the
     boxes — and both boxes lie outside its circumscribed ball, the block is zero.
  4. An occupied shaft refines to leaf pairs, which cast every pair's ray: wholly hidden
     → zero, wholly visible → cross approximation, partly seen (a shadow edge,
     `:sampled` fractions) → dense on the known mask.
  Visibility is never decided by sampling: rays between 7 representatives per cluster
  plus 32 random pairs (the first design) passed cap-to-cap blocks 0.5% occluded by the
  capsule (6 of 1280 entries, each a full entry), a 3% block error, and the sampled
  check of the cross approximation passed them too. Sampled hidden/visible decisions
  on occupied blocks of up to 4 leaves a side (with exact leaves) were 1.25× faster at
  28.6k but left element H errors of 5e-4 (moments 1e-7): dropped.
- **Low rank.** ACA with partial pivoting on the exact entries, stopped at a relative
  Frobenius estimate `tol` (1e-5), given up early when the terms' decay extrapolated
  (over the last 4) misses `3 tol` at `maxrank` (32, and ≤ 3mn/4(m+n) so it pays), QR +
  SVD recompression to `tol`, then a sampled check ((m+n)/2 entries, error ≤ 10 tol of
  the block). An empty result never proves a zero block (it did: capsule clusters whose
  sampled rows were all behind; 5.8e-4 of ‖F‖ dropped). Rows facing the other box
  wholly get the approximation and the few cut by the facing test (≤ half) stay dense.
- **Products.** Blocks grouped in row order into one run per thread of about equal
  cost, each into a separate buffer, summed in group order (deterministic for a thread
  count); 16 columns at a time; fused kernels (one pass over a dense block for its
  block and its transpose). Dual vectors go through the values-and-partials matrix.
  Keyframes of hierarchical models interpolate as an operator (`_BlendedExchange`, the
  two keyframe operators with per-row weights and a separate balancing), since their blocks
  differ.

### Accuracy (defaults tol 1e-5, η 1, leaf 48)

Against the dense Float32 path on the same mesh (balanced, T-dependent albedo solve at
2 ns), and at 28.6k against an `HMatrix` at tol 1e-7 (dense Float32 there needs 3.3 GB
of F next to the other job):

| N | capsule a₀ | P2/P0 | P4/P0 | H per element, max (median) | row sums | moment response |
|---|---|---|---|---|---|---|
| 10,484 | 3.1e-8 | 4.2e-8 | 2.6e-8 | 8.6e-5 (3.6e-8) | 1.7e-8 | 6.6e-5 |
| 17,632 | 1.4e-9 | 2.7e-10 | 2.8e-8 | 8.6e-6 (4.7e-8) | 1.1e-8 | 4.6e-6 |
| 28,628 (vs tol 1e-7) | 5.0e-9 | 9e-9 | 4.5e-9 | 3.5e-6 (6.3e-8) | — | 7.1e-6 |

(a₀ relative, P-ratios absolute, the rest relative.) tol 1e-4 moves the moments by
2–9e-7 for 3–5% less time; 1e-6 by ≤ 1e-8 for ~30% more. On small targets with rims,
bores, a waist gap, corner fillets, remeshed CAD (48 facets per capsule element) and
`:sampled` visibility (`test/hmatrix.jl`, leaves of 16), every G entry matches
`view_factors` to < 1e-5 of the largest.

### Time and memory (warm view-factor assembly; dense rows from "Large radiation models")

| N | 10,484 | 17,632 | 28,628 | 49,800 |
|---|---|---|---|---|
| Float32 dense: view factors | 13.4 s | 31.4 s | 101.6 s (29.1k) | ~300 s (est.) |
| Float32 dense: F | 419 MB | 1,186 MB | 3,278 MB | 9,920 MB (does not fit) |
| HMatrix: view factors | 7.6 s | 17.3 s | 34.8 s | 82.2 s |
| HMatrix: stored | 122 MB | 240 MB | 432 MB | 1,004 MB |
| HMatrix: max RSS, whole run | 0.85 GB | 1.00 GB | 1.36 GB | 2.26 GB |
| F x (HMatrix) | 7 ms | 19 ms | 34 ms | 81 ms |
| closure (Sinkhorn) | 0.5 s | 1.2 s | 2.1 s | 5.5 s |
| T-dep. solve, warm | 0.5 s | 0.7 s | 0.8 s | 1.6 s |
| 5-row moment response | 0.4 s | 1.1 s | 2.0 s | 5.7 s |
| mean rank of low-rank blocks | 7.5 | 6.9 | 6.0 | 5.6 |
| dense entries (incl. scattered) | 27M | 50M | 81M | 199M |

Assembly grows as N^1.4–1.6 and storage as N^1.2–1.5, not N log N: what remains dense
is the capsule's shadow edges (wall-to-wall leaf pairs partly hidden: 1.4k → 24k mixed
leaf pairs from 10k to 50k, the N^1.5 of a 3-D boundary in the 4-D pair space), leaf
pairs that do not compress (corner clusters of cap and wall seen past the capsule's
terminator; 3k → 100k ACA give-ups), and the near field. Wholly hidden blocks are rare
(the ball certificate fires on 4 → 2.7k blocks): the capsule's umbra is small against
the cavity. Profile at 28.6k (before the ray-free entries): rays were half the time;
with them, ACA ~45% (half of it on blocks that give up), dense leaves ~30%, leaf masks
~10%.

### Defaults

> **Superseded in part (2026-10-25).** The default `tol` is now **1e-6** and the cluster
> tree splits by surface first ("The H-matrix's N^1.5 remainder (2026-10-25)"); the
> thresholds below (`HMATRIX_MODEL[]` 15,000, Float32 from 10,000) are unchanged.

Above 15,000 elements `RadiationModel` takes the HMatrix (2× faster and 5× smaller at
17.6k, moments within 1e-8); 10–15k stays Float32 (1.8× faster and 3.4× smaller at
10.5k too, but the element-wise H error there reaches 8.6e-5 at an LEH element). Small
models are unchanged (dense Float64 + LU), and so is the dense path itself.
The precompile workload runs a hierarchical model (leaves of 16) on its coarse hohlraum.

### Open (at the time; see "The H-matrix's N^1.5 remainder" below)

- The N^1.5 terms: low rank of the shadow-edge blocks as kernel × mask (the mask is the
  costly part: a ray per pair), and a better split of corner clusters (cap + wall in one
  leaf) so their blocks compress.
- ACA give-ups cost ~20% of assembly; a cheaper predictor (the facing fraction, the
  shaft) could skip hopeless attempts.
- GUI `radiation_size` should budget an HMatrix (~0.2 GB per 10k elements here).

## The H-matrix's N^1.5 remainder (2026-10-25) *derived*

ROADMAP "Scale", the three open items above. Same workload, machine and script as
"Hierarchical view factors" (the session's scratchpad `hmat2/`: `attr.jl` time
attribution, `bench2.jl`, `vf.jl` interleaved timings, `rowchk.jl`); several other Julia
jobs shared the machine, so before/after times below are from runs interleaved in one
sitting (the old code from a worktree at 9e87912).

### Where the time went (before, 28.6k, 155 CPU-s inside block assembly on 4 threads)

| Part | CPU-s | share |
|---|---|---|
| near-field dense leaves (no ray / rays) | 18.1 / 3.9 | 14% |
| ACA that succeeded (free blocks; occupied leaves, with rays again) | 34.5 / 9.0 | 28% |
| ACA that gave up (non-leaf / leaf) | 16.8 / 7.3 | 16% |
| dense leaves after a give-up | 10.0 | 6% |
| leaf visibility masks (20M rays) | 14.3 | 9% |
| mixed (shadow-edge) leaves, kernel on the mask | 7.0 | 5% |
| rows cut by the facing test (split, sparse) | 19.3 | 12% |
| facing tests, shafts | 14.7 | 9% |

By cluster pair, the give-ups were capsule–wall (29.6k of 46.5k), cap/LEH–capsule and
mixed clusters (capsule + wall, cap + wall). The facing fraction predicts them well:
blocks with fewer than half their rows or columns wholly facing gave up 45.6k times
(23.6 s) and compressed 3.0k times (5.3 s). The shadow-edge leaves (11.8k) were mostly
wall–wall across the capsule. Why the capsule's ball certificate (below) could not
apply at first: 7k shafts per 10k-element model held capsule facets *and* a cluster
whose box reached into the capsule, because the median bisection put capsule and wall
elements in one cluster.

### What changed

1. **Clusters by surface.** `_hcluster_tree` splits a cluster whose elements lie on
   several surfaces (a component's curved surface, else the component; caps and their
   coplanar LEH discs share one) by surface first, the largest from the rest, then
   bisects as before. No cluster mixes capsule and wall or cap and wall at a corner.
   This is the "corner cluster" fix and most of the gain: far candidates 177k → 57k,
   low-rank blocks 42k → 10k, give-ups 46.5k → 15k at 28.6k. It is also more accurate
   element by element: the Frobenius tolerance of a mixed block let its small rows
   (capsule elements next to wall elements) be wrong relative to themselves. A capsule
   row at 28.6k had 2 entries dropped (G = 0 against 1.7e-8, 1.7% of the row's
   largest), its row sum 3.7e-5 low, in the old tree; 4e-11 after.
2. **Shadow edges decided by the capsule's balls** (`_shaft_occluder`, `_hvis_ball`,
   `_hmask`). The shaft search now tells apart an empty shaft, one holding facets of a
   single convex part only (the capsule) with both clusters outside its circumscribed
   ball, and anything else. In the second case every pair of the block is decided
   exactly, without sampling: a segment (between the lifted ray ends of `_vf_visible`)
   that misses the circumscribed ball hits no facet of the part, and one through the
   inscribed ball crosses its closed surface (margins 1e-9 relative, far above the ray
   test's 1e-10 barycentric slack). Only pairs between the two balls, or admitted by
   the facing slack, take the ray. The mask (Int8, blocks up to 64 leaf²) makes the
   block zero, low rank with entries off the mask (no more rays), or refined with the
   mask passed down to dense leaves. Children of an empty shaft inherit it, so no shaft
   is searched twice, and a BVH node whose facets are all of the certified part (or of a
   convex part a cluster lies on) is skipped whole. Leaf masks by rays now return the
   same mask type, so a wholly visible occupied leaf no longer re-casts its rays in ACA.
3. **Default `tol` 1e-5 → 1e-6.** With the surface tree, capsule–wall blocks compress
   where they used to give up and fall back to exact dense entries. At 1e-5 that moved
   P4/P0 by 1.7e-7 at 28.6k, against 1.6e-8 for the old tree. At 1e-6 every moment is
   within ~1e-8 of the references and the time is still well below the old tree's at
   1e-5.

### After (warm assembly, 4 threads; old = 9e87912 at tol 1e-5, new = tol 1e-6)

| N | 10,484 | 17,632 | 28,628 | 49,800 |
|---|---|---|---|---|
| view factors, old | 7.6–8.4 s | 16.6–17.5 s | 33–42 s | 90–98 s |
| view factors, new | 5.3–6.3 s | 13.2 s | 22–25 s | 56.5 s |
| stored, old → new | 122 → 99 MB | 240 → 212 MB | 432 → 439 MiB | 1004 → 1014 MiB |
| max RSS (whole bench), new | 0.80 GB | 0.94 GB | 1.25 GB | 2.15 GB |
| dense entries, new | 19M | 38M | 82M | 196M |
| F x, new | 8 ms | 18 ms | 42 ms | 114 ms |

Assembly now grows as N^1.4 from 17.6k to 49.8k (N^1.6 before). Storage is unchanged
(N^1.5): what is still dense is the near field and the shadow-edge leaves. Profile of
the new assembly at 28.6k (inclusive samples): ACA 35%, dense leaves 31%, rows cut by
the facing test (scattered) 23%, shafts 2%, ball masks 1.6%, ray masks ~0. The rays
are gone from the profile; what remains is evaluating the kernel itself.

### Accuracy (balanced, T-dependent solve at 2 ns; a₀ relative, P-ratios absolute)

| N, reference | a₀ | P2/P0 | P4/P0 | H per element, max (median) | moment response |
|---|---|---|---|---|---|
| 10,484, dense Float32 | 1.4e-9 | 2.6e-9 | 2.7e-9 | 1.2e-7 (4.7e-9) | 7.8e-7 |
| (old tree, tol 1e-5) | 3.1e-8 | 4.2e-8 | 2.6e-8 | 8.6e-5 (3.6e-8) | 6.6e-5 |
| 17,632, dense Float32 | 1.7e-10 | 6.5e-9 | 5.9e-9 | 1.3e-7 (6.7e-9) | 1.1e-6 |
| (old tree, tol 1e-5) | 1.4e-9 | 2.7e-10 | 2.8e-8 | 8.6e-6 (4.7e-8) | 4.6e-6 |
| 28,628, new tree at tol 1e-7 | 2.0e-9 | 9.7e-10 | 9.6e-10 | 7.8e-7 (6.6e-9) | 2.4e-6 |

The old tree's tol-1e-7 matrix at 28.6k is off from the new tree's by 1.5e-5 per element
whatever the tolerance; `rowchk.jl` puts that on the old tree (the dropped capsule
entries of item 1, checked against `_vf_exchange` pair by pair). The small-target test
entries (`test/hmatrix.jl`, leaves of 16) still match `view_factors` to < 1e-5 of the
largest, and the 3.9k hohlraum at leaves of 32, 2.3e-5 before, now matches to 1.5e-6.

### Tried and dropped

- **Skipping ACA on blocks the facing test cuts in most rows** (min wholly-facing
  fraction < 0.5; the predictor above). With the surface tree it skipped 22k attempts
  and left 279 give-ups. Time was unchanged (22.0–23.2 s against 22.1–23.0 s at
  28.6k), and it stored 6% more (464 MiB) and 10% more dense entries: the skipped
  blocks that would have compressed cost as much as the failed attempts saved.
- **Kernel × mask for shadow-edge leaves** (ACA on the unoccluded kernel, expanded on
  the mask, sampled check): 21.4 s against 22.2–22.5 s. With the balls the mask is
  cheap, the leaf's kernel evaluations are only ~15% of the dense entries, and the
  product of a Hadamard-stored block costs rank × a dense one. Not worth the code.

### Still open

- The near field and the rows cut by the facing test (a quarter of the time):
  quadrature-level reuse between neighboring pairs, or a coarser near-field rule where
  the exact contour kernel is not needed.
- `:sampled` visibility and non-convex occluders (LEH rims, CAD) still take a ray per
  pair of an occupied leaf. A rim certificate (a disc or annulus the segment must
  cross) would extend the ball test. (Done for centroid visibility: "H-matrix visibility
  certificates".)
- BriareusGUI `radiation_size` now budgets an HMatrix at 250 N^1.5 bytes (2.5 × the
  measured storage) above `HMATRIX_MODEL[]` and 4N² between 10k and 15k elements.

## H-matrix near field: zero blocks and cut rows (2026-10-26) *derived*

ROADMAP "Scale", the first "still open" item above. Measured here (same hohlraum workload and
`bench3.jl`-style script, 4 threads, `tol` 1e-6, leaf 48; old = 1f20b2d). Two changes:

1. **All-zero blocks dropped.** A dense leaf pair or scattered block whose Float32
   entries are all exactly zero is not stored (a capsule's near field against itself
   across the sphere, rows cut at a terminator against a box's far side). Exact in
   every entry: the product skipped `muladd(0, z, y)`. Only the split of blocks into
   product groups (by cost) moves, so threaded products sum in another order (round-off).
   At 17.6k: 6.6k of 25.7k dense leaves were zero, and 4.3M zero entries in scattered
   blocks.
2. **Facing-cut rows refined over the column tree** (`lowrank_index` blocks). The rows
   of a low-rank block that face only part of the other cluster (stored dense before)
   are re-tested against the column cluster's children, where most face wholly
   (indexed low-rank block, own ACA and sampled check) or not at all; dense only against
   a leaf or when ≤ 4 rows remain. These are new approximations at the same `tol`.

| N | 10,484 | 17,632 |
|---|---|---|
| stored, old → new | 99.4 → 81.8 MiB (−18%) | 212.0 → 175.1 MiB (−17%) |
| dense entries, old → new | 18.9M → 14.1M | 37.8M → 27.6M |
| indexed low-rank blocks | 52 | 215 |
| view factors (warm, 2 runs), old / new | 5.9–6.0 / 6.3–6.5 s | 11.5–11.6 / 11.7–11.8 s |
| F x, old → new | 6.6 → 5.6–5.9 ms | 16.6–17.1 → 13.5–13.7 ms |
| F X (16 columns), old → new | 93 → 80 ms | 256 → 212–217 ms |

New against old, balanced T-dependent solve at 2 ns (a₀ relative, P-ratios absolute):
a₀ −7.5e-12 / 1.4e-10, P2/P0 2.5e-11 / 3.8e-11, P4/P0 −1.5e-12 / 3.7e-11, H per element
max 1.4e-8 / 1.4e-8, moment response 1.4e-8 / 1.6e-8 (10.5k / 17.6k); against dense
Float32 at 10.5k unchanged (a₀ 1.4e-9, P2/P0 2.7e-9, H max 1.2e-7). The previous
session measured −15 to −27% storage over 10.5k–28.6k. Assembly time is unchanged to
+8%: the cut rows' extra ACA costs about what the dense entries did.

### Vectorized sector contour (H-matrix only)

A profile of the 10.5k assembly put ~45% of its work in `_contour_vf` over polar-sector
loops (34 edges, one scalar `atan` each: wall and capsule points onto the end caps). The
H-matrix assembly now gives sector receivers their loop as x, y, z arrays (`_SoAReceiver`)
and takes the contour in one `@simd` loop with a branch-free `atan(y ≥ 0, x)` (`_atan_nb`:
fdlibm's polynomial after reduction to |u| ≤ tan(π/8); ≤ 2 ulp from `Base.atan` over 10⁶
random arguments, exact at 0, π/4, π/2, π), 2.4× faster per contour (714 → 306 ns). The
results are not identical, since the terms are summed in another order. `view_factors`
(dense) keeps `_Receiver` and is unchanged (`==` on two hohlraums, 2.0k and 4.0k elements).
Measured against 5173a3d (4 threads, tol 1e-6, leaf 48):

| N | 10,484 | 17,632 |
|---|---|---|
| H-matrix build, warm (3 runs), old → new | 6.07–6.20 → 4.38–4.45 s (−28%) | 11.53–11.59 → 8.78–8.89 s (−24%) |
| entries, sampled rows (400 / 300 rows, all columns) | 0.77M of 2.0M differ; max abs ΔF 3.9e-17, max rel 2.5e-12 | 0.85M of 2.6M; 4.2e-17, 3.6e-12 |
| blocks, ranks, storage | identical | identical |
| a₀ (rel), P2/P0, P4/P0 (abs) | 2.1e-11, 9.4e-12, 3.9e-13 | −5.6e-16, 3.1e-14, 9.8e-15 |
| H per element max (median); moment response | 7.4e-9 (2.0e-11); 1.2e-8 | 2.6e-12 (3.6e-15); 2.5e-12 |

(Balanced T-dependent solve at 2 ns, `bench3.jl`.) The largest relative entry
differences are small factors whose contour terms cancel (7e-7 with 1.1e-17 absolute).

**Contiguous quadrature layout: tried, dropped.** The per-element point rules (`qnear`,
`qmid`, `qfar`: three small arrays per element) were packed into three arrays with offsets,
read as views. The change is neither faster nor identical in results: H-matrix build +2%
(4.47–4.55 s at 10.5k, 9.03–9.05 s at 17.6k, on top of the contour change); dense Float32
11.79 against 11.82 s and 30.2 against 30.8 s (one run each). Dense F changed in 15% of
entries (max rel 3.9e-13), from `_point_kernel` and `_self_kernel` with the weights as
views: StaticArrays' `dot` carries LLVM `reassoc contract` flags, so where FMAs are fused
depends on the surrounding code, and a change of data layout in a kernel does not
in general keep results identical.

## H-matrix near-field quadrature: where the time goes (2026-10-27) *derived*

ROADMAP "Scale", "the near-field quadrature itself". Same workload (`bench3.jl`, 4 threads
pinned to one CCX, `tol` 1e-6, leaf 48); old = d238f5d, timed interleaved in one sitting
(other jobs on the machine's other CCX, so absolute times run ~15% above the previous
sections').

**Profile (10.5k, before).** Attribution by category of `_hassemble!` × class of the pair
kernel (tier × pair type; instrumented, ~1.7× the real time, so shares only):
near-field dense leaves 27% of block assembly, ACA 37%, rows cut by the facing test
(refined, split, sparse) 13%, dense after a give-up 8%, shadow-edge leaves 6%. Of the
kernel time, a curved point (wall, capsule) onto an end-cap sector's 34-edge contour
("pc") is the largest class everywhere: far tier 6.3 s (3.7 s of it inside ACA), mid
2.1 s, near (adaptive rule) 1.3 s, against far point-to-point 2.3 s and the visibility
rays of occupied leaves 3.6 s. In the near-field leaves themselves: mid pc 1.6 s (0.5M
pairs, 3.2 µs), near pc 1.3 s (6.9k pairs, 194 µs each), mid point-to-point 0.8 s, far pc
0.5 s, rays 0.4 s, near point-to-point 0.3 s. Single pair costs (one thread, 20k random
pairs per class): far pc 1.24 µs, mid pc 2.84 µs, far sector–sector 2.37 µs, near pc
200 µs; of a near pc pair, ~half was `_distance_to` (32 point–triangle distances at
every cell the adaptive rule visits), a quarter the contours.

**Coarser rule: not viable.** The point rules' error is O((h/r)²): at the far tier's 12
element sizes a point-to-point or one-point-per-cell sector rule is off by ~1e-4 to 1e-3
relative, 100–1000× the block tolerance. Such a rule changes the discretization instead of
approximating the same matrix. An exact-enough quadrature of the same polygon (Gauss on a
sector to < 1e-9 needs ~16 points plus the 32 chord segments) costs as much as the contour.
Quadrature reuse between neighboring pairs has nothing to share in this workload: the
dominant pairs are distinct points onto distinct receivers.

**What changed** (all on the receiver side; the cost is in evaluating the same kernel):

1. `_closer_than` on a sector receiver (`_Receiver` with a loop) takes the distance from
   the loop in its plane's coordinates (height inside the polygon, else the nearest of
   34 edges, one pass) and tests the 32 triangles only when the verdict is within 1e-9
   (relative) of the threshold; `_adapt_tri` and `_adapt_sector` ask it through
   `_within`. The refinement, and so every factor, is **exactly** unchanged: `==`
   on `view_factors` (dense) of `hohlraum(h = 0.6, capsule_h = 0.3)` and both bored
   cavities, and 0 wrong verdicts of 2.9M (points 1e-4 to 1 mm from sectors, in and off
   their plane, thresholds at D(1 ± 1e-6), D/2, 2D and D). Near pc pair 167 → 94 µs.
2. The vectorized contour counts the loop's vertices on or behind the point's tangent
   plane in the same `@simd` loop (an `Int` sum; a floating-point `min`/`ifelse`
   reduction did not vectorize and made the loop 2× slower) instead of two scalar
   `all` passes over the loop.
3. A loop the tangent plane cuts (1% of points, 4 µs each through 32 clipped
   triangles) is clipped as one polygon (`_clipped_loop_vf`, streamed
   Sutherland–Hodgman; the coincident edges along the cutting line cancel), 2.5 µs.
   Against the triangles over 44.8k random cutting points: max abs 1.7e-14, within the
   test's 1e-12 rel + 1e-15 abs.
4. `_atan_nb` forms (t − 1)/(t + 1) of t = lo/hi as (lo − hi)/(lo + hi), one division
   instead of two; still ≤ 2 ulp from `Base.atan` (10⁶ arguments over 24 decades).

Pair costs after (same pairs): far pc 1.04 µs (−16%), mid pc 2.43 µs (−14%), far
sector–sector 1.83 µs (−23%), near pc 114 µs (−43%).

| | 10,484 | 17,632 |
|---|---|---|
| H-matrix build, warm, old → new (interleaved) | 4.30–4.43 → 3.89–4.06 s (−9%) | 8.84–9.08 → 8.15–8.35 s (−8%) |
| entries, 300 sampled rows (all columns) | 63 of 3.1M differ; max abs 8.9e-17 | 88 of 5.3M; max abs 1.9e-17 |
| blocks, ranks, storage | identical | identical |
| a₀ (rel), P2/P0, P4/P0 (abs) vs old | 0, −8.3e-17, 4.0e-17 | 0, 0, 5.2e-17 |
| H per element max (median); moment response vs old | 1.3e-15 (0); 2.5e-15 | 4.9e-15 (0); 3.9e-15 |

Against dense Float32 at 10.5k unchanged: a₀ −1.37e-9, P2/P0 2.67e-9, P4/P0 2.67e-9, H
max 1.19e-7 (median 4.7e-9), moment response 7.8e-7. (The largest relative entry
differences, up to 0.7, are factors of ~1e-17 whose contour terms cancel.) Dense
`view_factors` gains from item 1 only, with identical results:
`bored_cavity(hw = 0.8, hc = 0.25)` 1.72–1.75 → 1.35–1.38 s, `bored_cavity()` 2.35–2.36 →
1.98–1.99 s (4 threads).

**What is left.** After this the H assembly is mostly the contour itself (≈ 9 ns per
edge: a square root and two divisions per edge, vectorized), point-to-point pairs whose
cost is the scattered loads of three small arrays per element (the contiguous layout
above did not help), and the rays of occupied leaves. None of the remaining near-field
ideas measured here is worth 10%: a further step needs fewer kernel evaluations (rim
certificates for the rays; fewer ACA give-ups), not cheaper ones.

## Near flat pairs by the double contour (2026-10-27) *derived*

ROADMAP "CAD", "Adjacent flat CAD pairs". Scripts in the session scratchpad `adj/`
(`workloads.jl`: mock-ups only — a closed box of primitive quads, the flat and fitted
remeshed CAD hohlraum of the "Radiation in an imported CAD hohlraum" tutorial at
0.5/0.15 mm, the 2,000-strip sliver tube remeshed at 0.6; `profile.jl`, `compare.jl`,
`t2.jl`–`t5.jl`, `calib.jl`, `scope.jl`). Ryzen 9 4900HS, pinned to one CCX, other jobs on
the other.

**Share first (one thread, before; pair kernels classified by tier and by shared
vertices).** Pairs sharing an edge took 71% of a box's view-factor time (4 × 4 quads
per face: 72 pairs at 3.5 ms each), 58% of the flat CAD hohlraum's (1,959 pairs at
8 ms; all of them involve an element folded across the wall's strips) and 87% of the
sliver tube's (1,482 pairs at 74 ms, with 130-facet folded elements). The fitted CAD model
(walls on their cylinders) has almost none. Not a negligible share wherever flat
elements meet.

**Where the old rule's error was.** The cell rule refines toward the receiver to
maxdepth 8, but its one-point cells away from the shared edge set the error: one
element per plate, against Hottel's perpendicular rectangles, it is −5.6e-4 (1 × 1),
−4.4e-4 (1 × 0.6), −3.4e-4 (20:1 slivers), −5.0e-3 (shared edge 1/20 of the width);
maxdepth 13 changes none of these by more than 2e-5. On folded CAD elements the cell
rule emits from the element's mean plane: neighbors on the wall, folded by a few
degrees, exchanged 10–70% off a facet-level converged quadrature (median 16%, max 24%
over 12 sampled edge-sharing pairs; 7% / 10% for vertex-sharing, 3% / 53% for
near pairs apart). Their factors are small (G ~ 4e-4 mm² against 0.25 mm² elements),
so rows moved by ≤ 2.2e-3.

**Adopted: the double contour integral** (Sparrow 1963; Mitalas & Stephenson 1966):
G = A_i F_ij = (1/2π) Σ_a Σ_b (â·b̂) ∫_a∫_b ln r ds dt over the directed outline edges,
exact for any two surfaces that see each other wholly from in front, folded or not
(`_contour_pair_vf`). Used for every near pair of flat elements (not curved, not
sectors) whose vertices all lie in front of the other's facet planes (to 1e-9 of the
size) and whose own facets lie in front of one another (no self-hiding fold); otherwise
the cell rule as before. Inner integral closed form; outer adaptive Gauss–Kronrod 7–15
with a round-off floor; edges meeting at a point have their s ln s subtracted
analytically, collinear edges are closed form (x² ln|x|/2 − 3x²/4), an edge is split
where the other ends on it, and edges far apart for their lengths take product
Gauss–Legendre (n = 6…2 for ρ = 2|Δmid|/(La + Lb) ≥ 16, 24, 40, 120, 2000: worst error
over 3,000 random segment pairs ≤ 1e-13 La Lb against 40-point product Gauss on 4 × 4
panels, the round-off of the reference). Coordinates are scaled to the pair's size; the
constant part of ln r cancels over closed outlines. The reverse direction is the same
integral, so reciprocity is exact without averaging.

Accuracy (relative to the independent reference):

| case | old (cell rule) | new |
|---|---|---|
| perpendicular rectangles, Hottel: 1 × 1 × 1, 1 × 1 × 0.6, 1 × 0.05 × 0.02, 0.05 × 1 × 1 | −5.6e-4 … −5.0e-3 | ≤ 2.2e-15 |
| unit squares at 45°, 135°, 170°, 1 × 0.5 × 2 at 60°: cell rule against the contour | +4.5e-4, −9.4e-4, −1.0e-3, −1.8e-3 (maxdepth 13 the same) | (the prism below checks non-right angles) |
| CAD flat: edge / vertex / apart near pairs vs facet-level quadrature (mindepth 6, maxdepth 12) | median 1.6e-1 / 6.8e-2 / 2.9e-2 | 1.5e-5 / 1.4e-5 / 7.4e-6 (the reference's error) |
| box: same classes | 5.5e-4 / 1.7e-3 / 1.5e-3 | 6e-6 / 1.1e-5 / 5.7e-6 |
| closed cube, triangular prism (one element per face): row sums | — | 1 to < 1e-12 |
| cube with bottom + front one folded element: its row sum | — | 1 − F(Hottel) to 1e-12 |
| closed box, 4 × 4 per face: closure error max / mean | 5.3e-3 / 2.9e-3 | 5.3e-3 / 2.8e-3 |
| perpendicular plates n × n, Hottel (H = 1; 0.6), n = 2 / 4 / 8 | +2.3e-4 / +2.4e-3 / +2.2e-3; −1.7e-4 / +1.2e-3 / +2.2e-3 | 2e-16 / +2.2e-3 / +2.1e-3; 4e-16 / +1.0e-3 / +2.1e-3 |

What is left in the plate sums is the mid tier (3–12 element sizes, one centroid point
onto the exact polygon, O((h/r)²)): taking the contour there too makes every case above
exact to 3e-15 and the closed box's rows one to 7.7e-15, at 3.4× the box's time (0.07 →
0.24 s at 8 per face), 2.3× the CAD hohlraum's (2.0 → 4.9 s) and 16× the sliver tube's
(132-edge outlines; 6.1 → 98 s). Not adopted (an open decision).

Timings, dense `view_factors`, 4 threads (old = the same build with the contour
disabled, interleaved; identical otherwise):

| model | N | old | new | entries changed | max ΔF / max Δ row sum |
|---|---|---|---|---|---|
| `hohlraum(h = 0.6, capsule_h = 0.3, sombrero = true)` | 1,888 | 0.45 s | 0.45 s | 0 (`==`) | — |
| `bored_cavity(hw = 0.8, hc = 0.25)` | 1,068 | 1.32–1.40 s | 1.32 s | 0 (`==`) | — |
| box, 8 × 8 quads per face | 832 | 0.24 s | 0.07 s | 7,392 | 1.4e-4 / 4.3e-4 |
| CAD hohlraum, flat 0.5/0.15 | 1,693 | 7.3–7.4 s | 2.1–2.3 s | 25,402 | 3.5e-3 / 2.2e-3 |
| CAD hohlraum, fitted 0.5/0.15 | 1,693 | 2.60 s | 2.59 s | 248 | 1.7e-3 / 1.6e-3 |
| sliver tube 0.6 (open) | 518 | 34.1 s | 7.2 s | 14,190 | 4.9e-4 / 6.4e-4 |

Per pair (one thread): an edge-sharing pair of quads 3.3–7.6 ms → 2–100 µs (slivers the
most), folded CAD wall pairs 8 ms → 0.4 ms, tube pairs 74 ms → 3.8 ms (132 × 132 outline
edges). With non-touching near pairs, the box is faster (4 × 4 per face, 0.013 → 0.006 s
with them against edge-sharing pairs only), the CAD hohlraum's 25% faster, the tube's ~15%
slower (its total still 4.7× faster). Primitive hohlraums and bored cavities have no flat
non-curved elements and are unchanged; the H-matrix assembly goes through the same pair
kernel. The flat CAD tutorial's radiation metrics move by little (drive of the tutorial,
balanced F): capsule flux +4.2e-5, a₂/a₀ 10.587 → 10.584%, a₄/a₀ −0.788 → −0.787%, LEH and
gap losses unchanged to 1e-5. The sliver tube's mean row sum moves 0.72210 → 0.72197 (the
open cylinder's exact 0.72301; the rest is the mid and far tiers). Confidence: the method
and the closed forms **H**; the ρ thresholds and the 1e-9 front tolerance engineering
choices, **M**.

## H-matrix visibility certificates (2026-10-27) *derived*

ROADMAP "Scale", "rim certificates for non-convex occluders". Workloads: `bench3.jl`'s hohlraum
(10,484 and 17,632 elements) and a target where rims occlude, the bored cavity of
`test/partial_visibility.jl` (LEH bores, cap washers, a 20 µm waist gap with its washers;
latitude–longitude capsule) at `hw = 0.18, hc = 0.07` (11,360 elements). 4 threads pinned
to one CCX, `tol` 1e-6, leaf 48; old = abb3129 from a separate worktree, interleaved. Other
jobs shared the machine: process-to-process spread ~15%, so ranges below are of 3 builds
per process, 2 processes each.

**Measured first.** Visibility rays timed by instrumenting `_vf_visible` (CPU-s inside
it over a warm build):

| | hohlraum 10.5k | bored cavity 11.4k |
|---|---|---|
| rays (old) | 2.9M, ~3.8 CPU-s (~23% of assembly CPU) | 35.4M, 19–23 CPU-s (43%) |
| of which occupied leaf pairs | ~0.4% of assembly | ~32% far leaves + 9% near leaves |
| of which pairs only the facing slack admits, in empty shafts | ~23% (wall–capsule 1.4M, wall–wall 1.0M, caps–capsule) | ~10% |

So the rim question was a negative on the standard hohlraum (its occupied leaves are
~0.4%) and a large positive where rims occlude; and on the standard hohlraum the rays
were nearly all of a different kind: pairs that face each other only within the slack of
the centroid facing test on curved elements, whose rays (lifted farther off the facet) are
cast even in an empty shaft.

**What changed** (`src/hmatrix.jl`; all exact, so `HMatrix` is **unchanged**):

1. Certificates per component, built once (`_hrim_data`): a *planar* part (all facets
   within 1e-10 relative of one plane: washers, LEH caps, discs, foils) keeps its facets
   as 2-D edge lines binned on a polar grid; a *cylinder or cone* part keeps its axial and
   radial bounds about the surface axis; a convex part with balls keeps them. The shaft
   search of a leaf pair (`_shaft_occluder(…, full = true)`) now collects up to 8 such
   parts instead of giving up at the second one or at a non-convex one.
2. Each pair of a certified leaf is decided as `_vf_visible`'s ray (`any_hit` to
   1 − 1e-9, barycentric slack 1e-10) would: hidden if one certificate finds the segment
   through its facets, visible if all find it clear, else the ray. Planar: one side of
   the plane (margin 1e-9 relative plus twice the flatness), or the crossing outside every
   facet in its bin, or inside one with both ends off the plane — margins grown by
   1/sin of the crossing angle, and undecided when they exceed half a bin. Hull: the
   segment's part within the axial range (widened) inside the radial lower bound or
   outside the upper one (visible only). Ball: as before, with the per-pair condition that
   both ends lie outside the circumscribed ball (so a pair whose cluster box reaches into
   it no longer forces every pair to the ray), plus a separating plane: an end element
   lying on the convex part (its facets in one plane through its centroid) clears it when
   both ends are in front by 1e-8 × the target scale.
3. Pairs only the slack admits are judged on their own ray ends (`_hslack_ends`), when
   both ends stay inside their cluster boxes (lift ≤ the padding), every plane exclusion
   of the shaft held with the clusters' largest extra lift (`liftx`, recorded by the
   shaft search), and each convex part a cluster lies on (excluded from the search) is
   cleared by an end element's plane. This applies to empty shafts (the `_HOcc` of the
   shaft replaces `true` as the block's entry mode) and to certified masks; centroid
   visibility only.
4. Blocks above leaf level are unchanged (only the old single-ball masks), so the block
   structure, ranks and every stored entry stay as they were. (Superseded by the second
   pass below: certificates above leaf level, with a quadrature-tier check.)

| | hohlraum 10,484 | hohlraum 17,632 | bored cavity 11,360 |
|---|---|---|---|
| H-matrix build, warm, old | 4.46–4.69 s | 9.42–9.68 s | 9.79–11.40 s |
| new | 3.44–4.54 s (typ. 4.0) | 7.87–8.37 s | 7.74–7.92 s |
| change | ≈ −10% | −16% | −20 to −30% |
| rays old → new | 2.9M → 0.98M | | 35.4M → 1.34M (0.9 CPU-s) |
| stored blocks, all entries | identical (`==` on every dense, low-rank, scattered block) | identical | identical |

The remaining hohlraum rays are capsule elements at their terminator facing slightly away (c
≤ 0 within the slack), whose own plane does not separate; all visible. The wall-clock gain
is smaller than the ray CPU saved (on the bored cavity, ~20 CPU-s of rays → ~10 CPU-s wall ×
4): the certificates cost ~140 ns per pair, and the shaft search no longer stops early on a
leaf.

**Exactness checks** (`test/hmatrix.jl`, "visibility certificates"): 57k random segments
across the parts of the bored cavity and a hohlraum, and segments crossing a washer's plane
within 1e-3 to 1e-12 (relative) of its inner and outer rims: every decided one (49k, 5,023
of them hidden) agrees with `any_hit` against the facets of that part only; masks of every
certified leaf pair equal the ray masks (6,193 leaf pairs over three targets, a
latitude–longitude capsule among them); the 833k pairs of empty-shaft leaves get
`_vf_pair`'s factor (`===`). At the benchmark sizes the certified masks equal ray masks on
all 9,889 certified leaf pairs of the bored cavity (12.4M pairs).

**Certificates above leaf level (adopted, second pass).** The certificates also decide
blocks above leaf level (`_shaft_occluder(…, full = true)` for every far block under
centroid visibility; a mask for blocks up to `mask_size`), so a block of uniform
visibility gets a cross approximation instead of being refined to leaves. A first try
also counted latitude–longitude capsule quads (two coplanar facets) as "on the convex
part" in the shaft search; that rule is **not** adopted (below). The first try lost row
accuracy against dense Float32 (2.4e-7 → 2.3e-6); the loss was diagnosed as follows.

- *One entry.* Capsule element 7631 is a polar-fan triangle of size 0.265 (3–4× its
  neighbors). Wall element 1380 lies at r = 3.1825 from it, just past the far boundary
  12 × 0.2649 = 3.1783. The dense matrix evaluates that pair with one point per triangle
  (far tier) and the other 42 pairs of the same column with the level-1 rule (mid tier).
  The rules differ by 0.53% there: F = 4.2547e-4 (far) and 4.2319e-4 (mid). The
  converged value (level 4) is 4.2246e-4, so the dense entry is 0.71% off and the mid
  value 0.17%. In the low-rank block (43 wall × 107 capsule), the cross approximation's
  pivots and the check's 75 random samples all miss the isolated jump, and the block
  stores the mid-tier value. Its error against dense is 5.1e-5 of the block's norm; the
  other 4,728 low-rank blocks are all ≤ 7.6e-6. The row sum of 7631 is off by 2.28e-6.
  The dense row sum of that sliver is itself 1.00188 against 0.999998 converged, so the
  H-matrix's departure is ~1000× below the error of the dense quadrature itself, and toward
  the converged value.
- *Not the certificates.* The masks of every block that takes one are exact: 11,134
  blocks (8,763 leaves), 33.5M pairs (3.4M hidden), all equal to the ray masks. The 20.5M
  pairs of the shafts that the capsule rule empties get `_vf_pair`'s factor (`===`).
- *Not only the capsule rule.* Above-leaf certificates by themselves reach the same block:
  its parent's shaft holds capsule facets, which the balls certify pair by pair. They give
  the same 2.28e-6. With both rules the blocks are identical to those of above-leaf
  certificates only, so the capsule rule adds nothing there. By itself it saves ~3% (7.76
  against 7.97 s), and it is dropped.
- *Fix.* The low-rank check (`_hcheck`) now adds the exact error of the entries where a
  row's or a column's quadrature tier switches (`_vf_tier`, `_htier_entries`). These are
  the pairs off their line's most common tier, taking the lines with the fewest such
  pairs first (an isolated jump is the one the pivots miss), at most m + n of them. There
  are none, and no pass over the pairs, when the centroid boxes are farther apart than
  `mid` × the largest element. The block above is rejected and refined. A block with no
  such entries is checked exactly as before.

Bored cavity (11,360; 4 threads pinned, interleaved in one process, 2–3 builds each);
balanced solve as `bench3.jl` (NIF quads aimed at z = ±4 mm deposited by
`DepositionOperator`, `closed = true`, 2 ns), against the dense Float32 matrix:

| | committed | capsule rule | above-leaf | above-leaf + tier check (adopted) |
|---|---|---|---|---|
| build, warm | 7.0–8.2 s (typ. 7.9) | 7.8 s | 6.7 s | 5.7–6.6 s (typ. 6.4; ≈ −18%) |
| stored | 117.8 MB | 115.7 MB | 104.9 MB | 104.9 MB (−11%) |
| Frobenius error, 300 columns | 6.49e-8 | 7.07e-8 | 7.76e-8 | 7.75e-8 |
| row sums, max \|Δ\| | 2.42e-7 | 2.28e-6 | 2.28e-6 | 2.63e-7 |
| H per element, max (median) | 1.01e-7 (3.5e-9) | 5.5e-7 | 5.5e-7 (5.7e-9) | 1.10e-7 (5.7e-9) |
| capsule a0; P2/P0, P4/P0 | −1.7e-11; 1.1e-9, 6.6e-11 | | | −7.2e-10; 3.2e-9, −1.5e-9 |

Hohlraum: above-leaf certificates by themselves leave the blocks identical to the committed
build (10,484). The tier check rejects one more block at each size (10 → 11 at 10,484; 16 →
17 at 17,632), so the hohlraum is **not identical**. At 10,484 the Frobenius error
(2.679e-7) and the row sums (3.214e-7) are unchanged to four digits, storage is 85.7 → 85.8
MB, and the build stays at 4.0–4.2 s. At 17,632 storage stays 183.6 MB and the build 7.8–8.0
s. The dense matrix is unchanged (`_vf_tier` is `_vf_kernel`'s tier test, factored out).

Tests (`test/hmatrix.jl`): the mask-against-ray and empty-shaft checks of "visibility
certificates" now cover every block up to `mask_size`, not only leaf pairs. "Quadrature-
tier entries" uses a latitude–longitude hohlraum (3,254 far blocks): every entry returned is
off its row's or its column's majority tier; every pair that is the only one off its line's
majority is returned when the cap is not reached; and a block with no entries has one tier.
Confidence: the diagnosis and the measurements **H**; the m + n cap and the fewest-first
order are engineering choices, **M**.

**`:sampled` visibility (third pass, adopted).** Same bored cavity with
`visibility = :sampled` (11,360 elements, 4 threads pinned, warm).

*Measured first* (old build 106 s with `_vf_visible` instrumented): 48.7M calls,
404 CPU-s, 93% of the assembly CPU. By caller: 26.7M at occupied leaf pairs (227 CPU-s),
21.8M entries of near leaves and blocks without a certificate (176 CPU-s), 0.24M slack
pairs of empty shafts (1.3 CPU-s). 1.5M pairs disagreed with their probes and took all
256 sample rays. Under `:sampled` the old shaft search took no certificate but the
whole-ball one, so every occupied leaf took rays.

*What changed* (`src/hmatrix.jl`, `src/viewfactors.jl`):

1. Bundle certificates. All rays of a pair's `_visible_fraction` (the centroid ray, the
   probes, the 16 × 16 sample segments) run between points within `rvis[i]` of oi and
   `rvis[j]` of oj. Here `rvis` is the largest sample distance from the centroid. Each
   point of such a segment lies within ρ = max(rᵢ, rⱼ) of the centroid segment at the
   same parameter. So each certificate is applied to the centroid segment with margins
   grown by its tube:
   - ball: the ends outside the outer ball by rᵢ, rⱼ; clear past rout + ρ; hidden
     inside rin − ρ;
   - planar: the end heights shifted by rᵢ, rⱼ. The crossings lie within
     ρ (1 + |w|/|dh|) of the centroid crossing, since a point at height ≤ ρ on the
     axis lies that far along it from the crossing. The 1/sin margin uses
     (|w| + rᵢ + rⱼ)/(|dh| − rᵢ − rⱼ);
   - hull: axial range widened by ρ, rlo − ρ, rhi + ρ;
   - an end element's plane clearing a convex part: the samples of that end within
     `hdev` (≈ 0) of its plane, the other end's within rⱼ.

   Zero radii reproduce the old arithmetic exactly. A bundle is 0 if any part hides
   it whole, and 1 if all clear it whole. Otherwise `_visible_fraction` runs with each
   segment decided by the parts that did not clear the bundle (`_HSeen`), or by the ray.
2. Exclusions valid for the samples. In `:sampled` mode the shaft's plane tests take the
   box of every element's samples (and centroid) lifted by ε (`opts.qbox`) instead of the
   centroid box. A cluster's in-plane offset (`cplane`) also comes from those points.
   Per pair: the ends stay inside the boxes (rvis + lift ≤ the padding); a slack pair
   needs `occ.slack`; each excluded convex part must be cleared by an end element's plane
   for the bundle, else per segment, else the ray. Folded elements take `_vf_visible`.
3. Assembly: `:sampled` blocks take the full shaft search. Blocks up to `mask_size`
   with only certified parts get a fraction mask (`_HFrac`): wholly visible blocks get
   the cross approximation at any level, others are refined; near leaves are masked too.
4. Zero-weight sample pairs (facing away at an end) add nothing to either sum of the
   fraction, so their segments are tested only when no pair faces (the area-weighted
   fallback). The output is unchanged, and it removed 41M segments into the capsule.
5. One rounding for every caller. StaticArrays' `dot` is an `@simd` reduction, so it
   rounds as its caller is compiled: `ci` differed in the last bit between `_vf_visible`
   and the bundle code (a slack lift, so a fraction off by 6e-17). The ray ends
   (`_vf_ends`) and the sample weights (`_vis_weight`) are now `@noinline` with concrete
   arguments, one compiled instance for all.

| bored cavity 11,360, `:sampled` | old (7ec13bb) | new |
|---|---|---|
| H-matrix build, warm (interleaved) | 93.3, 93.6 s | 29.9, 29.7 s (3.1×) |
| stored | 121.0 MB | 108.3 MB |
| low-rank / masked blocks | 15,437 / 0 | 4,971 / 4,086 |
| `_vf_visible` calls | 48.7M | 0.61M |
| pairs decided by bundle (hidden / visible / per segment) | — | 1.41M / 23.5M / 11.3M |
| segments per segment (certified clear / hidden / ray) | — | 187M / 116M / 67M |
| Frobenius error, 300 dense columns | 5.98e-8 | 7.82e-8 |
| row sums, new − old | | ≤ 2.5e-7 |

The error is at the block tolerance (1e-6) as before; the comparable centroid change was
6.5e-8 → 7.8e-8. Most of the 67M rays left are segments through the 20 µm waist gap that
end on the wall. The wall's axisymmetric bounds cannot decide them, since the end lies
between the inscribed and circumscribed radii. A certificate for the inside of a convex
faceted wall (prism) would be needed: ~40M rays, roughly a quarter of the build.

*Exactness checks.* On the benchmark all 13,316 certified blocks (leaf pairs and blocks
up to `mask_size`) were checked. Their 33.2M facing pairs (2.92M hidden, 1.02M partial)
have `_hfrac` `===` `_vf_visible`'s fraction. The 7,644 empty-shaft blocks have
`_vf_pair_free` `===` `_vf_pair` per pair. No block of this target needs rays any more.
`test/hmatrix.jl` "sampled visibility bundles" (1,068 elements) checks:
- bundle verdicts against the rays of random segments inside the balls (some ends on
  their surfaces), against the facets of that part only;
- every certified block's fractions and every empty shaft's factors (`===`).

Centroid-visibility H-matrices are identical block by block (bored cavity 11,360,
hohlraum 10,484; 5.9 s and 3.5 s, unchanged). Dense centroid factors are unchanged.
Dense `:sampled` factors are not: 3.8% of entries move, by ≤ 4e-14 relative (item 5).
They are 5% faster from item 4 (10.05 → 9.50 s at 3,156 elements). The dense path has no
shafts, so the bundle certificates do not apply to it. Confidence **H** for exactness and
timings; the gap-certificate estimate is **M**.

## Near sector pairs by the double contour with arcs (2026-10-27) *derived*

ROADMAP "Scale", "the bored cavity's gap pairs". Scripts in the session scratchpad `gap/`
(`share.jl`, `attr.jl`, `near.jl`, `ref.jl`, `cmp.jl`, `blocks.jl`, `nocap.jl`); 4 threads
pinned to one CCX, old = fb62e1b from an exported copy, run back to back (not
interleaved; other jobs on the other CCX).

**Share first (before).** The "gap pairs" are the near pairs with a gap component: the
slot's aperture (180 chord rectangles of the cylinder at R_out, 20 µm tall, facing in)
onto the annuli's sectors (50 µm wide, 9:1 at 0.5 mm), sharing their outer arc. Timed
one pair at a time on one thread:

| | `bored_cavity()` 2,244 | `(hw = 0.8, hc = 0.25)` 1,068 | `(0.18, 0.07)` 11,360 |
|---|---|---|---|
| aperture ↔ annulus, near | 1,520 pairs × 2.28 ms | 1,536 × 2.25 ms | 1,584 × 2.26 ms |
| annulus ↔ annulus | 200 × 0.48 ms | 120 × 1.79 ms | 520 × 0.52 ms |
| all gap near pairs | 3.58 CPU-s | 3.69 CPU-s | 3.86 CPU-s |
| build, 4 threads | dense 1.96 s, H 1.05 s | dense 1.33 s, H 0.78 s | H (centroid) 5.85 s |
| share (CPU-s / 4 over the build) | dense ~45%, H ~85% | dense ~69% | H ~16% |

The cost is fixed by the gap (its 180 aperture elements), so the share falls with N but
is still above 10% at 11k. Other near pairs of the same kind: cap ↔ wall at the corner
(1,520 pairs, 184 µs), bore ↔ LEH disk (704, ~300 µs). All were the cell rule refined to
maxdepth 8 toward the shared arc onto a 34-edge sector receiver.

**Accuracy of the old rule.** Against a converged quadrature (2-D adaptive Gauss–Kronrod
over the emitter's true patch, (φ, z) on the cylinder or (r, φ) on the sector, of the
exact point-to-polygon factor onto the sector drawn with 1,024 chords per arc; tolerance
1e-8): aperture–annulus pairs −7.1e-4 (sharing the arc), −6.2e-4, −4.4e-2, +2.3e-3 at
0.04–0.9 mm; annulus pairs across the gap +8.4e-4 (opposed), −2.1e-2 (neighbors);
wall–cap −1.0e-3 … −4.9e-3; bore–LEH up to −1.9e-2. The 16-chord arcs of the receiver
leave a slit along the shared arc (the −0.38% of "Deposition-operator rays").

**Adopted: the double contour with exact arcs** (`_arc_pair_vf`). The contour form of
"Near flat pairs by the double contour" depends only on the two outlines, not on the
surface spanning them, so it also holds for a cylinder patch seen from inside, as long
as each element lies wholly in front of the other (every point of a sector's partner on
its front side; every point of an inward cylinder patch's partner within its radius).
Outlines of sectors and chord rectangles are arcs about the axis and straight edges:

- arc × arc (coaxial): the double integral of ρ₁ρ₂ cos(φ − ψ) ln r depends on
  u = φ − ψ only, so it is one integral of it weighted by the overlap w(u) of the two
  ranges (linear between four breakpoints). ln r² = ln(c + 4ρ₁ρ₂ sin²(u/2)),
  c = Δρ² + Δz², is singular at u = 0 on coincident circles (the shared arc) and sharply
  peaked across the 20 µm gap; ½ ln(c + ρ₁ρ₂u²) times w is integrated in closed form on
  each piece and the smooth rest by adaptive Gauss–Kronrod (a first version without the
  subtraction took 500 µs per gap pair; an x⁴ grading toward the singular ends was
  slower still);
- arc × line: the line's integral in closed form (as for flat pairs), the arc's adaptive,
  split at the azimuths of the line's ends; axial lines are perpendicular to every arc
  and radial line, so their terms vanish;
- line × line: the flat pairs' `_lnr_edges`.

Used for near pairs of a sector and a sector or an inward chord rectangle of a cylinder,
coaxial (to 1e-12), each in front of the other to 1e-9 of their size; otherwise the cell
rule as before (the wall above the slot and the annulus below it, outward cylinders,
cones, spheres, flat facets). Reciprocity is exact (one value for both directions).

Accuracy. Summed over every pair of the bored cavity at 0.5 mm (thousands of pairs, all
tiers, all taking the new rule), against the coaxial-disk closed forms: aperture → annulus
7e-13 to 1.7e-12, annulus ↔ annulus −2.9e-13, wall → cap 1.7e-14, bore → LEH disk (center
fans) 5.8e-14, cap ↔ cap 1.8e-12, LEH disk → far cap −9.8e-13. Per pair against the
converged quadrature above: ≤ 1.6e-7 on every pair checked (16 pairs, 0.02–2.1 mm apart),
the error of the reference itself. **H** for the method and checks.

Cost per near pair, one thread: aperture–annulus 2.29 ms → 27 µs (85×), annuli 0.48 ms
→ 81 µs, cap–wall 184 → 18 µs, bore–LEH ~300 → 20 µs, LEH disk–wall 23 → 11 µs; all such
near pairs of `bored_cavity()` 4.6 → 0.19 CPU-s.

| 4 threads, warm | old | new | entries changed | max \|ΔF\| / max \|Δ row sum\| |
|---|---|---|---|---|
| `bored_cavity()` 2,244, dense | 1.96 s | 0.86 s (2.3×) | 14,000 | 5.5e-4 / 8.7e-4 |
| `bored_cavity(hw = 0.8, hc = 0.25)` 1,068, dense | 1.33 s | 0.27 s (5×) | 10,032 | 1.3e-3 / 2.7e-3 |
| same, H-matrix | 1.05 s / 0.78 s | 0.53 s / 0.18 s | | |
| bored 11,360, H-matrix, centroid | 5.85 s | 5.15 s (−12%) | | |
| bored 11,360, H-matrix, `:sampled` | 33.0–34.1 s | 28.6 s (−14%) | | |
| `hohlraum(h = 0.6, capsule_h = 0.3)` 1,120 | 0.43–0.45 s | 0.33–0.35 s | 6,400 | 7.1e-4 / 1.1e-3 |
| with `sombrero = true` 1,888 | 0.49 s | 0.38–0.39 s | 6,400 | 7.1e-4 / 1.1e-3 |
| `hohlraum()` 4,116 | 2.35–2.37 s | 2.18–2.21 s | 10,176 | 5.6e-4 / 7.6e-4 |

**Not identical** wherever a sector meets a coaxial sector or inward cylinder near
by: every primitive hohlraum (cap–wall corner pairs; LEH disk–wall pairs), bored cavities,
coaxial disks. That is the point of the change. Whole-matrix component exchanges (all
tiers) against closed forms, old → new:

| | 2,244 | 1,068 |
|---|---|---|
| aperture → annulus | −7.9e-4 → −7.5e-6 | −3.2e-3 → −7.5e-6 |
| annulus ↔ annulus | +1.1e-4 → −1.4e-7 | −1.4e-4 → −5.6e-8 |
| wall → cap | +2.2e-4 → +4.6e-5 | +2.7e-4 → +6.9e-5 |
| bore → LEH disk | −2.7e-4 → −5.9e-5 | −8.3e-4 → −2.0e-5 |

(The rest is the mid and far tiers.) Row sums: the hohlraum without its capsule (closed
by its LEH apertures) closes better, max |1 − Σ| 8.7e-4 → 6.3e-4 and mean 3.4e-4 →
1.6e-4 at h = 0.6 (8.5e-4 → 5.8e-4, 2.5e-4 → 1.7e-4 at h = 0.3); its wall → cap exchange
+1.2e-4 → −1.2e-4 (+2.5e-4 → +1.3e-4 at h = 0.3). With the capsule, the changed rows of
`hohlraum(h = 0.6, capsule_h = 0.3)` move mean |1 − Σ| 4.0e-4 → 4.3e-4 and max 6.3e-4 →
1.7e-3 (cap rows +2.3e-4 → +4.5e-4: the old rule's error had canceled part of the
other tiers' on those rows); `hohlraum()` 5.5e-4 → 4.0e-4 mean. Bored cavity gap rows
change by < 1e-3 (their 5% closure error is the slot narrower than a sample cell,
"Partial visibility"). Tests: `test/runtests.jl` "view factors" (the four closed forms
summed over every pair to 1e-11, the out-of-front case, the near pair in
`view_factors`); `hmatrix.jl`, `partial_visibility.jl`, `corners.jl`, `meshops.jl`,
`cadfit.jl` pass unchanged.

Effect on the calibrated benchmark (coordinator, same day): `examples/n221204.jl 0.8 0.4
full`, main (d19d84f) against this change, every time of the CSV: P2/P0 within 1.4e-5,
P4/P0 within 4e-6, capsule Tr and Dante (both LEH areas) within 1.4e-3 eV, the P2
crossing unchanged to 1e-3 ns (6.241 ns), capsule energy unchanged (248.7 kJ, 14.0%).
The calibrations fitted with the old rule (wall loss, CBET, the N221204 defaults) are not
moved by it; nothing is refitted.

## Prism certificates for the waist gap (2026-10-26) *derived*

ROADMAP "Scale", "a certificate for the wall rays through the waist gap". Workload: the
`:sampled` bored cavity of "H-matrix visibility certificates" (`bored_cavity(hw = 0.18,
hc = 0.07)`, 11,360 elements), `tol` 1e-6, leaf 48, 4 threads pinned to one CCX; old =
3fcf470 from an exported copy, interleaved (another job shared the CCX at times: builds
of one process agree to 1%, between processes up to 50% slower, so the comparison is
of the quiet runs, with process CPU time beside them).

**Measured first** (`_HSeen`'s fallback rays by component pair and by the certificate
left undecided, instrumented build): 66.7M rays. The wall's axisymmetric bounds left
37.7M visible ones to the ray, all between the gap (its aperture, 16.7M; the washers,
21.0M) and the opposite half's wall; bore–wall, bore–cap and bore–washer pairs 5.4M
more (the bore's bounds); the capsule's balls ~21M (segments between the inscribed and
circumscribed balls). The ray of a gap segment costs 0.34 µs (one thread, 0.9M
segments of gap–wall pairs), so the gap rays were ~13 CPU-s of 112 (~12%), not the
quarter estimated above (that came from `_vf_visible`'s cost per call, several rays
each).

**What changed** (`src/hmatrix.jl`; exact, so `HMatrix` is **unchanged**). A hull
part (`_HHull`, a cylinder or cone component) whose facets are all parallel to its axis and,
seen along it, lie within 1e-10 (relative) of a set of lines that keep the axis inside
(a wall or bore of chord rectangles) also keeps those side lines (`_hprism`: facets
grouped by the azimuth of their outward normals; each line through the mid-range of its
vertices; `pf` the largest half-range, 2e-16 here). The hull verdict then also clears a
segment whose part within the axial range (the same clipping, widened by the margin)
has both ends at least m + 2 pf (+ the tube radius ρ for a bundle) inside every line:
the set where that holds is convex, so the whole part stays that far from every line,
and so from every facet; outside the axial range the facets are m away. A segment
through the gap passes the wall's edge below its axial range, so its part within the
range runs from just inside the edge to the far wall (lifted ε = 1e-7 × scale off it).
The lines are binned by azimuth: a bin lists every line whose normal lies within
Δ + 1e-4 of it, cos Δ = 0.999 pmin / rmax, so a point at radius r with r cos Δ ≤ pmin − M
is at least M inside every line its bin leaves out (≤ 5 lines per bin here); the bin is
found by an arctangent within 2e-5 rad (Abramowitz & Stegun 4.4.49 on the octant: a
division and a polynomial instead of `atan`), and a point inside the lines' inscribed
circle by M needs no bin. The verdict is "visible" only, like the bounds' (no hidden
verdict: the gap–wall rays that are hidden, 0.6M, still take the ray). Margins are the
bounds' (m = 1e-9 of the size), and so is the rounding argument (an estimate: the
facet test's rounding, ~1e-16 of the size over sin θ, reaches the margin only for a
segment within ~4e-7 rad of a facet's plane that also runs nearly level with an end of
the axial range).

Cost per segment (one thread, 0.9M gap–wall segments): bounds 41 ns, with the side
lines 68 ns, against the ray's 338 ns.

| `:sampled` bored cavity 11,360 | old | new |
|---|---|---|
| fallback rays (`_HSeen`) | 66.7M | 23.2M |
| H-matrix build, warm, 4 runs each | 28.4–28.7 s | 25.4–25.5 s (−11%) |
| process CPU over the build | 111.8–112.6 s | 99.3–100.3 s (−11%) |
| gap–wall pairs, `_hvis_bundle` CPU (4 classes, 1.6M pairs) | 32.3 s | 21.7 s |
| stored | 108.3 MB | 108.3 MB |
| blocks, every stored entry, row sums | | identical (`==`) |

So the errors against dense are those of "H-matrix visibility certificates" (`:sampled`:
Frobenius 7.8e-8 on 300 columns). Centroid builds are identical in every block and not
faster: bored cavity 11,360 5.07–5.10 → 5.13–5.15 s, hohlraum 10,484 (`bench3.jl`'s, no
gap) 3.33–3.39 → 3.37–3.42 s (its wall pairs lie inside the wall, whose facets the shaft
excludes). Pairs whose segments the lines cannot clear pay for the test (cap–washer
pairs 20.6–20.7 → 21.3–21.6 µs per pair).

*Exactness checks.* At the benchmark size (the `:sampled` pass's check, in eight shards) every block
that takes certified fractions, 33.2M facing pairs (2.92M hidden, 1.02M partial), has
`_hfrac` `===` `_vf_visible`'s fraction, and the 7,644 empty-shaft blocks `_vf_pair_free`
`===` `_vf_pair`. 800k random segments over the 1,068-element cavity's walls and bores
(ends on facets lifted −1e-9 to 1e-6, within 1e-12 to 1e-1 of their axial ends at the
polygon's vertex and inscribed radii, in the gap's slot, anywhere about the part): all
304k decided agree with `any_hit` against that part's facets. `test/hmatrix.jl` "prism
certificates" keeps 16k of these (half of the decided ones decided by the lines only)
and the arctangent's bound; "visibility certificates" and "sampled visibility bundles"
pass unchanged and now run through the lines too.

What is left is the capsule's ~21M rays (segments between its balls; a per-facet
certificate for the convex capsule would save ≲ 7% of this build) and the per-segment
certificate work itself (a gap–wall pair is 256 segments through up to six parts:
10–19 µs per pair after this). Confidence **H** for exactness and timings.

## Deposition-operator rays and sliver view factors (2026-10-25) *derived*

ROADMAP "Scale". The test-suite profile found two costs: building a `DepositionOperator`
(~6 s per 48-quad and ~22–30 s per 192-beam operator on one thread, almost all in the
visibility rays of every triangle × cone direction × beam) and the view factors of
cavities with sliver elements (the bored cavity's 20 µm waist gap: 70–80% of its
time). Scripts in the session scratchpad `depcost/` (`cmp.jl`, `harness.jl`, `vfcmp.jl`,
`closed.jl`, `conv2.jl`, `plates.jl`); timings on the Ryzen 9 4900HS with other jobs
sharing it (ratios are from interleaved runs).

### Deposition operator: packet visibility (identical results)

- **Screening on the profile's cut-off does not apply.** The operator's points and
  masks are reused as aims, spots and defocus move (by `reach` with full resolution,
  beyond it at `level`), so a point that receives an exact zero at the construction aim
  is still needed after re-aiming; dropping it changes `deposit` there. Only points whose
  mask is empty (no cone direction reaches them) are dispensable, and they are already
  dropped. For the record, the rays on triangles outside every beam's `near` region
  (zero for every aim within `reach`) are 2% of a quad's (Q15T) and 22–26% of a
  single NIF beam's: restricting the operator to `reach` would save at most a quarter on
  192 beams and lose the level-resolved deposition beyond it. Not adopted.
- **Where the rays go** (default `hohlraum()`, 6228 opaque triangles): ~0.5 M rays per
  quad, ~0.43 M per beam; 55–65% blocked (wall points looking back at the far wall),
  the rest lit. Per triangle and direction the answer is uniform for 98% of the rays
  (all points lit 3–7%, all blocked 33–46%, back sides ~50%, mixed 1.5–2%). A
  single `any_hit` costs ~300 ns (≈ 55 node boxes and 4.5 triangle tests; 1.7 of those
  on triangles the ray moves away from).
- **Adopted: one packet per triangle** (`_unblocked!`, all its quadrature points × the
  cone directions facing its lit side). Interior nodes are culled by a slab test that is
  conservative over the origins' box and the directions' inverse components, leaves by
  one such test per direction, and then each ray takes a separate slab test and triangle
  tests. fp subtraction and products are monotone, so no ray that passes a box is
  culled, and the slab test is monotone in box containment, so a ray passes a leaf only
  if it passes all its ancestors: every ray tests exactly the triangles a
  per-ray `any_hit` reaches until one blocks it. Each direction's last blocking leaf is tried
  first (a cache across a work unit's triangles; 38% of blocked rays settle there).
  The operator is identical (points, weights, normals, masks, destinations
  and `deposit`), on the 48-quad and 192-beam NIF operators. The facing test d·n < 0
  is now summed in a fixed order (`_facing_mask`): StaticArrays' `dot` is an `@simd`
  reduction whose rounding depends on the calling code, which flipped a direction lying
  in a back side's plane (d·n = 0 or −3e-17, a zero-weight contribution) when the loop
  was restructured.
- Tried (same answers; replaying 37k packets of six quads, `any_hit` 0.95–1.0 s): the
  adopted packet 0.62–0.65 s; one packet per direction 0.68; direction masks carried down
  the tree 0.67–0.71; without the leaf cache 0.63 (the cache barely matters once the
  packet is shared); visiting the child farther along the rays first 0.63 (the current
  order), nearer first 0.92. The occluder cache with single rays alone gave 1.1×.
- Not exact, not tried: skipping the triangle tests a ray provably moves away from (the
  facet the ray leaves), safe by an error bound of nine orders but not by construction.

| NIF operator, `hohlraum()` | 1 thread | 4 threads |
|---|---|---|
| 48 quads, before → after | 6.05 → 4.26 s (1.42×) | 1.61 → 1.14 s (1.41×) |
| 192 beams, before → after | 22.2 → 13.3 s (1.66×) | 5.72 → 3.54 s (1.61×) |

### View factors: slivers cut across their length

`_adapt_tri` (and `_adapt_sector`) split cells at their midpoints, which keeps a cell's
shape: a sliver stays a sliver at every depth, and the test `distance < ratio × longest
edge` holds for every thin row of it near a receiver, so a 5:1 element took ~10 rows of
cells where a square one takes ~3, and an annular sector 15:1 refined radially 16× more
than needed. Now a triangle with h² > 4.5 × twice its area (a chord rectangle beyond
~4:1) is cut across its length instead: a needle parallel to its short edge (its tip and
the two halves of the base trapezoid), a flat one at its longest edge; a sector whose
arc is more than 4× its width (or the reverse) is halved one way only. Either halves the
longest edge, so `mindepth`/`maxdepth` keep their meaning. Primitive triangles are under
the bound (walls ≤ 2.5, geodesic capsules 1.5, corner bands 2.45; latitude–longitude
polar caps reach 5.2 but are never adaptive emitters), and their midpoint children are
similar to them: **`view_factors` is unchanged** on `hohlraum()` (4116
elements), h 0.8/0.4, h 0.6/0.3 with `sombrero = true`, with `corner_radius = 0.5`, and
with `capsule_mesh = :latlong`. Reciprocity is exact as before (the pair is still
averaged over both directions); the bored cavities change by ≤ 3.5e-4 in F and ≤ 4.3e-4
in row sums.

Accuracy against closed forms (relative error of the summed exchange):

| case | before | after | converged quadrature |
|---|---|---|---|
| 1 × 0.05 plates 0.02 apart, one 20:1 element each (parallel opposed) | 1.0e-4 | 5.8e-4 | — |
| same, 5:1 elements | −8.6e-5 | 4.0e-4 | — |
| same, square elements (unchanged) | 8.8e-4 | 8.8e-4 | — |
| 1 × 0.05 plate to a 1 × 0.02 plate on its long edge, 20:1 | 8.8e-5 | −3.5e-4 | — |
| same, 5:1 | −2.6e-4 | −5.5e-4 | — |
| same, square elements (unchanged) | −1.6e-3 | −1.6e-3 | — |
| bored cavity (0.8/0.25): gap cylinder → gap annulus | −3.65e-3 | −3.22e-3 | −3.76e-3 |
| bored cavity: annulus ↔ annulus across the gap | −2.8e-4 | −1.4e-4 | −3.9e-5 |

Closed forms: directly opposed and perpendicular rectangles; coaxial parallel annuli
and a cylinder's inner surface to its end annulus from the disk–disk formula. The
"converged" column is the same kernel at maxdepth 12, ratio 8 (before and after agree
there): the −0.38% of the gap exchange is the receiver's polygonized arcs (16 chords
per arc leave a ~0.1 µm slit along the shared corner), not the emitter quadrature.
Slivers now carry the error of square elements of their width or less (3–6e-4, where
they were over-resolved to ~1e-4); the gap–annulus exchange is 5e-4 from its converged
value (1e-4 before) and the annulus pair 1e-4 (2.4e-4 before).

| view factors | 1 thread | 4 threads |
|---|---|---|
| `bored_cavity(hw = 0.8, hc = 0.25)` (1068 el.), before → after | 12.2 → 7.0 s | 3.19 → 1.87 s |
| `bored_cavity()` (2244 el.) | 14.6 → 9.5 s | 3.68 → 2.43 s |
| `hohlraum(h = 0.6, capsule_h = 0.3, sombrero = true)` (1888 el., identical) | 2.07 → 1.98 s | 0.53 → 0.55 s |

The gap pairs are still ~60% of the bored cavity's time: refinement toward the shared
long edge (maxdepth 8) onto sector receivers of ~34 edges each. A maxdepth relative to
the cell's width, or coarser arcs on receivers far from the point, would be the next
step; not tried. (Superseded 2026-10-27: those pairs take the double contour with exact
arcs, "Near sector pairs by the double contour with arcs".)

## Closure on bored, windowed and CAD cavities (2026-10-25) *derived*

Does `closed = true` (Sinkhorn balancing, `enforce_closure!`) apply beyond the flush-LEH
`hohlraum()` of "Balancing on the curved path against a continuum Monte Carlo"? CAD
decision 11 (2026-10-18, closed). `examples/mc_cavities.jl` (4 threads; ~30 min, cached).

The geometry is the bored cavity of "Partial visibility" (R 2.875, 9.4 mm long, 0.1 mm LEH
bores of radius 1.55, a 20 µm waist gap, capsule 0.9; apertures on the LEHs and the gap
exit); the same cavity from CAD-like facets (each half a closed solid of 360 strips, test
`cad_half`) through `cad_target(...; inside, remesh, fit = true)` (`:sampled`); and
`hohlraum()` of the same size with an open 1.5 mm circular window at the waist and a closed
1 × 2 mm window (α 0.5) at z = 1, each in an open 50 µm gap. Wall α 0.9, η 0.75, capsule α
0.3. Every model's capsule is the same geodesic mesh (`icosphere_mesh`, n = 15, 4,500
elements); source: the CAD tutorial's two Gaussian rings on the inner wall, integrated over
each element.

References: `montecarlo_incident(...; reemit = :hit)` (the continuum) on flat copies
faceted at 0.1 mm (capsule n = 29, 16,820 elements), watertight: bores, LEH disks and gap
aperture on the walls' polygon (`open_fraction` 0 and 3.4e-5). Bored cavity: 72.7M
capsule paths, a₀ 0.016511 ± 0.000001 per unit wall power, a₂/a₀ +5.334 ± 0.018%,
a₄/a₀ −0.472 ± 0.021%; the CAD cavity is the same geometry (its 360-gon lies within
1.1e-4 mm of the circle). Windows: 109M paths, a₀ 0.015232 ± 0.000001, a₂/a₀
+8.595 ± 0.015%, a₄/a₀ −0.686 ± 0.019%. The loss through the apertures (LEHs, gaps,
open window) per unit x-ray source, from the same Monte Carlo with the apertures as
receivers: 0.36880 ± 0.00005 and 0.40237 ± 0.00006.

**Is the enclosure closed?** In the continuum, yes: the apertures cover the openings
exactly, and the curved model integrates over the true surfaces. On facets, not as the
tests and examples build it: `open_fraction` (new; diffuse rays from the elements that
meet no facet, apertures counted as walls) is 9.2e-4 on the bored cavity (12–25% on
single bore elements: 24-sided LEH disks in 20-sided bores, bores of `nphi ÷ 2` sides on
caps of `nphi`, the 180-sided gap aperture on 40-sided gap faces), 7.6e-4 on the CAD
cavity (24-sided disks in 360-sided bores), 1–5e-4 on `hohlraum()` with or without
windows (LEH disks of `nphi ÷ 2` sides in caps of `nphi`), 2.7e-3 on `hohlraum(h = 1)`.
Openings without apertures: the CAD cavity without its LEH apertures 7e-2, the
40 µm waist gap of test/cadfit.jl 4.2e-3 (6.0e-4 with an aperture). Cracks and small
openings overlap, so the fraction locates leaks (`element`) but cannot decide closure.

Model − reference (a₀ and loss in %, moments in pp), raw → balanced:

| case | a₀ | a₂/a₀ | a₄/a₀ | loss |
|---|---|---|---|---|
| bored 0.5, centroid | +0.31 → −0.02 | −0.288 → +0.010 | +0.051 → −0.107 | +1.93 → +0.24 |
| bored 0.5, sampled | −0.01 → +0.06 | −0.153 → −0.168 | +0.050 → +0.023 | −0.29 → −0.01 |
| bored 0.25, centroid | −0.21 → −0.13 | +0.079 → +0.176 | +0.053 → −0.076 | +2.26 → +0.42 |
| bored 0.25, sampled | −0.05 → −0.06 | +0.028 → +0.022 | +0.050 → +0.028 | +0.60 → +0.16 |
| CAD 0.5, centroid | +0.34 → −0.10 | −0.262 → −0.012 | +0.022 → −0.092 | +2.82 → +0.33 |
| CAD 0.5, sampled | +0.16 → +0.03 | −0.254 → −0.176 | +0.048 → +0.039 | +0.73 → +0.09 |
| CAD 0.25, sampled | −0.06 → −0.06 | −0.014 → +0.020 | +0.045 → +0.032 | +0.76 → +0.19 |
| windows 0.45, centroid | +0.35 → +0.07 | −0.310 → −0.251 | +0.010 → +0.009 | +0.26 → +0.01 |
| windows 0.45, sampled | +0.17 → +0.07 | −0.237 → −0.249 | +0.006 → +0.011 | +0.12 → +0.00 |
| windows 0.28, centroid | +0.13 → −0.02 | −0.072 → −0.088 | +0.017 → +0.003 | +0.23 → +0.07 |
| windows 0.28, sampled | +0.06 → −0.01 | −0.052 → −0.089 | +0.016 → +0.004 | +0.16 → +0.06 |
| windows 0.17, centroid | −0.04 → −0.03 | −0.018 → −0.045 | −0.005 → −0.001 | +0.08 → +0.08 |
| windows 0.17, sampled | +0.04 → −0.03 | −0.045 → −0.042 | −0.004 → −0.002 | +0.16 → +0.08 |
| bored 0.5 flat, watertight, sampled | +0.71 → +0.23 | −0.071 → −0.015 | +0.043 → +0.040 | +0.53 → −0.08 |
| bored 0.5 flat, cracked, sampled | +0.71 → +0.49 | −0.072 → +0.128 | +0.043 → −0.051 | −0.63 → −0.61 |
| bored 0.5 flat, watertight, centroid | +1.04 → +0.18 | −0.207 → +0.158 | +0.044 → −0.060 | +2.51 → +0.10 |
| bored 0.5 flat, cracked, centroid | +1.04 → +0.41 | −0.207 → +0.329 | +0.044 → −0.220 | +1.59 → −0.37 |

(5,984 / 9,064 / 5,752 / 8,252 / 6,247 / 8,445 / 14,177 elements for the curved
models.) Largest |row sum − 1| before balancing: apertures 2–6% (centroid) and 1–3%
(sampled) on the rims, the 50 µm gap faces behind the 20 µm slot ~5% in every model
("Partial visibility"), 0.04–0.5% with windows. The reference's a₄/a₀ error (±0.02 pp)
is the size of most a₄ entries.

- **Balancing applies wherever apertures close the cavity, on the curved path.** It
  brings the aperture loss within 0.2% of the reference in every sampled case (raw
  −0.3 to +0.8%) and within 0.4% with centroid visibility (raw +1.9 to +2.8% on rims),
  and a₀ within 0.07% (sampled). With sampled visibility it moves a₂/a₀ by at most
  0.08 pp (CAD 0.5 mm, toward the reference); with centroid by 0.1–0.3 pp on rims,
  toward the reference at 0.5 mm and away at 0.25 mm (+0.08 → +0.18): sampled
  visibility is what gets the rims' a₂/a₀ right, and balancing then fixes the
  energetics. **H** for these geometries and albedos.
- **Windows and gaps in `hohlraum()`:** balanced loss within 0.08% and a₀ within 0.07%
  at every mesh; a₂/a₀ is set by the wall mesh (−0.25 pp at h 0.45, −0.04 at 0.17) and
  balancing changes it by ≤ 0.06 pp. `closed = true` stays the `hohlraum()` default with
  windows, gaps and lips (apertures) and sombreros (their rows are zero; test/sombrero.jl
  checks the hohlraum's solution unchanged to 1e-12).
- **Flat facets must be watertight.** On the cracked flat bored cavity balancing pushes
  a₂/a₀ up (+0.13 to +0.33 pp) and the loss low (−0.4 to −0.6%); the same mesh made
  watertight (matched polygons) balances to the loss within 0.1% and a₂/a₀ within
  0.02 pp (sampled). The upward bias of "Monte Carlo radiosity reference for partial
  visibility" is this (both its geometries are cracked). The curved path is immune: its
  apertures and rims are exact surfaces.
- **CAD:** the remeshed, fitted cavity behaves as the primitives (raw loss +0.7%,
  balanced +0.1–0.2%). `cad_target` cannot know whether the user's apertures close the
  cavity (a missing small opening and polygon cracks give the same `open_fraction`), so
  closure is declared: `cad_target(...; closed = true)` and `visrad_target(...; closed =
  true)` (a new keyword, default `false`, so existing targets are unchanged).
  The GUI's CAD path still builds raw targets.
- Thin-film windows that transmit are not covered (a transmitting aperture is a partial
  loss, not a closure) and VisRad workspaces.

## Compact capsule rows at capsule_h = h/3 (2026-10-25) *derived*

CAD decision 12 (2026-10-18, closed): on the compact capsule the LEH power-balance
tutorial's raw error grew from h 0.4 to 0.25 mm. The tutorial's N210808 hohlraum (6.4 ×
11.24 mm, 3.1 mm LEHs, capsule 1.13 mm), capsule_h = h/3 unless stated, Float64 view
factors, 4 threads (scripts in the session scratchpad). The reference for the capsule rows
rests on closure: the enclosure is closed, so every capsule element's row sums to 1, and the
capsule's exchange with the wall, the caps and the LEHs is the solid-angle fraction from its
center (cos θ_c, cos θ_L − cos θ_c, 1 − cos θ_L), exact for a sphere.

| h / capsule_h | mesh | capsule rows − 1 (mean; largest) | capsule→wall / solid angle − 1 | wall rows − 1 (mean) |
|---|---|---|---|---|
| 0.40 / 0.133 | compact | −9.2e-5; 4.4e-4 | −1.0e-4 | −2.3e-4 |
| 0.40 / 0.133 | lat-long | +1.3e-5; 1.3e-3 | | −2.9e-4 |
| 0.25 / 0.083 | compact | −3.5e-4; 6.8e-4 | −4.1e-4 | −3.0e-4 |
| 0.25 / 0.083 | lat-long | −2.3e-4; 1.9e-3 | −2.3e-4 | −1.6e-4 |
| 0.40 / 0.083 | compact | −7.2e-5 | −8.1e-5 | −3.9e-4 |
| 0.25 / 0.133 | compact | −5.6e-4 | −6.5e-4 | −1.9e-4 |
| 0.25 / 0.083, `mid = 24` | compact | −6.2e-6; 8.7e-5 | −1.9e-6 | −3.2e-4 |
| 0.40 / 0.133, `mid = 24` | compact | +1.4e-5; 5.4e-5 | +1.8e-5 | −3.0e-4 |
| 0.25 / 0.083, `level = 2` | compact | −4.5e-4 | −5.2e-4 | −3.6e-4 |

The capsule's exchange with the caps and the LEHs is within 4e-5 of the solid angles
in every case.

- **The capsule rows' error is the wall mesh's, not the geodesic triangles'.** It is the
  capsule–wall exchange, and it grows when the *wall* is refined (0.4 → 0.25 mm at either
  capsule size), not the capsule (0.133 → 0.083 at h 0.25 shrinks it). Curved pairs
  farther than `mid` (12) sizes of the larger element use one point per triangle
  (`qfar`), whose O((h/r)²) bias does not cancel; at h 0.25 that is 3 mm, so most
  capsule–wall pairs (2–6.5 mm apart) are far, while at h 0.4 (4.8 mm) most were
  mid-range (four points per triangle). `mid = 24` removes it (−6e-6), `level` does not
  (near pairs only). Lat-long elements are two triangles, two far points each, hence
  their smaller error at the same size.
- **The source-weighted loss is the wall rows', mostly the capsule's shadow.** Removing
  the capsule changes the wall–cap–LEH exchange by its shadow, which in a closed
  enclosure equals the capsule's area exactly. The centroid rule's shadow is
  0.4–0.9% too large: +5.7e-3 (compact) and +3.7e-3 (lat-long) at h 0.25, +5.8e-3 and
  +7.2e-3 at 0.4, +8.7e-3 for both at capsule_h 0.067 (h 0.4); on flat facets
  +8.3e-3 / +1.0e-2, so it is the all-or-nothing sampling of the silhouette, not the
  surfaces. Sampled visibility overshoots the other way (−1.7e-3 / −2.3e-3 at 0.4).
  Compact and lat-long differ by its scatter (wall–wall exchange per unit area, compact −
  lat-long, −1.1e-4 at 0.25 with centroid rays, −1.2e-5 sampled).
- **Tutorial errors now** (fractions of the source; 2 × 10⁷ histories of the tutorial's
  exact-geometry walk, σ 7e-5 LEH and 6e-5 capsule): h 0.4 raw LEH −0.8e-4, capsule
  −1.7e-4 (lost 6.1e-4), balanced +0.9e-4 / −0.1e-4; h 0.25 raw −2.8e-4 / −1.8e-4
  (lost 1.1e-3), balanced −1.0e-4 / +1.3e-4 (1.4σ / 2.2σ), the same on lat-long
  (−1.1e-4 / +1.2e-4) and with `mid = 24` (−1.1e-4 / +1.3e-4); the tutorial's default
  Float32 storage at 10,928 elements agrees to 1e-6. Against the seed-2
  walk of the tutorial the balanced capsule share is +0.5e-4. The ROADMAP's −3.3e-4 /
  −5.0e-4 (lost 1.4e-3) were of an earlier build.
- **Decision: documented, no default changed.** Balanced results at capsule_h = h/3 are
  within 2.2σ (1.3e-4 of the source) of the exact walk on either capsule mesh; the raw
  rows' 2–6e-4 come from the far-pair rule and the centroid shadow, neither specific to
  the compact mesh. `mid = 24` fixes the capsule rows for +30–60% view-factor time but
  not the loss (at h 0.4 it raises the raw loss, 6.1e-4 → 9.1e-4: the shadow's excess
  then stands alone), so it does not earn a default change. **H**

## Per-cavity balancing (2026-11-03) *derived*

`RadiationModel(target; closed = cavities)` and `enforce_closure!(F, A; cavities)`:
the Sinkhorn balancing of "Balancing on the curved path against a continuum Monte Carlo"
applied separately to the block of each cavity in the exchange matrix `Aᵢ Fᵢⱼ`, so the
rows of each cavity sum to one over its elements; elements in no cavity keep a factor of
one, and the one diagonal scaling applies to the whole matrix (reciprocity everywhere).
One cavity of every element is exactly the whole-target balancing. A cavity is the
caller's declaration (element indices, or component or part names);
`open_fraction(target; cavity)` checks it by rays: diffuse rays from the cavity's
elements that meet nothing, or first meet an element side outside the cavity (or a
one-sided surface from behind), leak. No value is fitted. Tests: `test/cavities.jl`
("per-cavity balancing"), GUI `test_halfraum.jl` ("per-cavity balancing").

**The GUI's cavities.** One per hohlraum part in the calculations: its pieces by
component name (entrance holes, lips and windows included), and every other element
whose centroid lies in the part's closed volume (its profile about its axis, ends
included, to 1e-6 of its length) and whose radiating side faces into it (1e-3 of the
length in front of the centroid still inside): a capsule, a package face in a
halfraum's open end, each side of a two-sided disk shared by two cavities. Overlapping
volumes give an element to the first part. The drive balances a cavity when
`open_fraction(t; cavity)` is at most `CAVITY_MAX_LEAK` (1e-2, the engineering default
of "Closure on bored, windowed and CAD cavities") and no element of its boundary is
left out of the calculations (a piece of the part, or another element in it whose
centroid one element size behind lies outside the volume; a capsule left out keeps the
cavity closed). Elements in no cavity (a package's outer face, a shield at an entrance
hole, free parts in the open) keep their raw rows. Leaks measured (`_closed_cavities`):
the default design 1.2e-3 (coarse) and 2.7e-4 (standard), the frustraum 3.8e-4 and
5.7e-5, each shock-tube halfraum 7.5–7.8e-3 and 1.1e-3 (the cracks where the holes'
disks, of half the cap's azimuthal count, meet the caps; the curved view factors do not
leak there), the opacity halfraum 0.17 (its open end). A ball of 0.9 mm reaching
through the default hohlraum's wall leaks 1.2e-2 (the wall buried in it sees it from
behind) and is left unbalanced, as before. The rays cost 0.09 s on the default design
(4796 elements, 4 threads).

Checks against independent references:
- **Two spheres with holes, side by side** (R 1.5 mm, 9 of 10 equal-area bands, 24
  azimuths, the hole closed by a flat disk aperture; albedo 0.8; a black disk in the open
  above one hole, which that cavity's wall sees through it). Every point of a sphere sees
  an area A of it with the view factor A/4πR², so the wall carries
  H = S w/(1 − α w), w = (1 + cos θ₀)/2, everywhere (closed form, every bounce). Raw
  in-cavity row sums miss by up to 1.6e-2 (the polar band), −2.2e-3 area-weighted. The
  wall's mean H: raw +3.3e-3, per-cavity balanced +3.7e-4; whole-target balancing,
  which forces the rows that include the free disk to one, −36%. The source of the lit
  sphere equals what the elements of that sphere absorb and let out, Σ(1 − α)HA, to 2.9e-13;
  the other receives nothing.
- **Two halfraums joined by a package, black** (2.5 mm long, radius 1.5 mm, holes of
  1 mm at the outer ends, the package 2 mm long; the package face emitting). From the
  face, the coaxial-disk view factors (r₁ = 1.5 to r₂ = 1.0 and 1.5 mm, 2.5 mm apart)
  give the hole 0.108027, the cap 0.111348 and the wall 0.780625 of its emission.
  Relative errors, raw → balanced: h 0.4 (348 elements per cavity) hole −6.0e-4 →
  −7.2e-5, cap −3.0e-4 → +3.1e-4, wall −3.2e-5 → −3.5e-5; h 0.28 hole −2.6e-4 →
  −4.7e-4, cap −2.8e-4 → +1.2e-4, wall +2.1e-4 → +4.9e-5; h 0.2 hole −1.6e-4 →
  −2.8e-4, cap −1.3e-4 → +1.9e-4, wall +2.0e-4 → +1.1e-5. Single-bounce view factors
  stay within their own error either way (balancing corrects row sums, not single
  entries); the loss falls from ~1.2e-4 of the source to 1e-13, and the other cavity
  receives exactly nothing.
- **The same pair reflecting** (walls 0.8, faces 0.3, the upper wall emitting 1 per unit
  area and the lower 0.5; `examples/mc_two_cavities.jl`). Flat watertight facets (each
  hole's disk of the cap's full azimuthal count: `open_fraction` 0) against the Monte
  Carlo of the element model with exact view factors (`reemit = :element`, 200,000 paths
  per receiver, ±1.4e-4 and ±8e-5): the package faces' mean flux raw +4.1e-3 and +3.4e-3
  (h 0.4 and 0.28), balanced −2.6e-4 and +1.0e-5. With the holes of the hohlraum itself (half
  the cap's count, cracks: `open_fraction` 2.7e-3) the balanced flat model is +8.8e-3
  off: balancing spreads the cracks' real loss over the cavity, so on flat facets it
  needs watertight meshes, as "Closure on bored, windowed and CAD cavities" found. On the
  default curved path the balanced face flux at h 0.4, 0.28, 0.2 lies +2.4e-4, +1.5e-4,
  +4.4e-5 from the balanced value at h 0.14 (6408 elements), the raw one +5.4 to
  +7.7e-4 at every h: raw does not converge to the same limit at the 1e-4 level (its row
  excess, about −2e-4 area-weighted, is amplified by the bounces).

As for changes to results, the default design, the frustraum and the opacity halfraum
(standard size) have the same drive and drive history, every serialized field, before
and after. The shock tube (standard, 1840 elements): the observed face's mean flux
−4.40e-4, its RMS 4.7098 → 4.7020%, peak −1.03e-3, lowest −4.7e-4; each halfraum's wall
Tr 311.2875 → 311.2444 eV (−1.38e-4); in the history the face's peak Tr
323.8916 → 323.8560 eV and each cavity's peak wall Tr 310.5660 → 310.5230 eV; Dante 1's Tr
at most 7.2e-5. A template-era session with a shield outside the hohlraum (GUI test data)
was unbalanced and is now balanced: its capsule's mean flux +0.33%.

## Line-of-sight power on partly hidden elements (2026-11-03) *derived*

- **Before.** `line_of_sight_power` counted an element whole when the ray from its
  centroid toward the detector met nothing, and not at all otherwise.
- **Now** (decision taken: sample points on the element, not `DanteView`'s pixel rays).
  Each element's unobstructed fraction is the area-weighted fraction of the centroids of
  `nsub²` = 16 equal sub-triangles per triangle that see the detector, and of
  `(refine·nsub)²` = 1024 where those first points disagree. Pixel rays err on every
  element boundary across the image and need a field and a resolution chosen for each
  view; sub-element points are exact on whole elements, put the error only on the
  partly hidden ones, and keep the function's `elems` subset and any direction.
- **Validation** (`test/line_of_sight.jl`): a 2 × 2 mm plate (3 × 3 elements) under an
  opaque sheet 0.5 mm above it whose straight edge runs at 20° to the mesh, seen straight
  on and at 30°, against the visible area in closed form (the square clipped by the
  shifted edge). Relative error of the power: old (centroid) +12.9% and +3.3%; one
  point per triangle −1.2% and +3.3%; the default 6.8e-5 and −9.5e-5; 64 per side
  −1.0e-4 and −1.9e-5; 128 per side 7.6e-6 and 6.3e-6 (first order).
- **Changed results.** Power from elements partly hidden; unobstructed and wholly hidden
  elements unchanged. Cost: 16 rays per triangle, 1024 on partly hidden ones.

## Entrance holes at the caps' azimuthal count (2026-11-03) *derived*

- **Before.** `hohlraum` meshed each entrance hole (and its closing lip) with half the
  cap's azimuthal count, so the hole's rim was an `nphi/2`-gon inside the cap's inner
  `nphi`-gon: a crack at every rim that rays leak through (the curved view factors,
  integrating the `DiskSurface` sectors, do not see it). The VisRad workspace's
  primitive hohlraum copied the rule.
- **Now.** Holes, lips and caps share the azimuths (`phibreaks` included) and every
  vertex. Diffuse-ray leak (`open_fraction`, by cavity): the GUI shock-tube starter's
  halfraums 7.5e-3 / 7.8e-3 → 0 (coarse) and 1.1e-3 → 0 (standard); the GUI default
  hohlraum 1.2e-3 → 0 (coarse) and 2.7e-4 → 0 (standard); core `hohlraum(h = 1)` 2.7e-3
  → 0. Cost: the holes' elements double (default design at the standard size 4,796 →
  5,228 elements, its drive 6.7 → 7.7 s).
- **Changed results** (the GUI's drives, balanced per cavity). Default design, standard:
  capsule flux +1.7e-5 relative, P2/P0 −14.5022 → −14.5016%, P4/P0 0.16691 → 0.16687%,
  wall Tr 267.7390 → 267.7407 eV; coarse: capsule flux +1.6e-4, ε₂ 14.33 → 14.37%, wall
  Tr 267.7287 → 267.7388 eV. Shock tube, coarse: observed face −7.6e-6, wall Tr 311.3656 →
  311.3661 eV; standard: +5.1e-6, 311.2444 → 311.2443 eV. Flat-facet references (Monte
  Carlo, `open_fraction`) no longer lose the rim cracks' share.
