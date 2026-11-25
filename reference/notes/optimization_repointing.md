# Research notes: Error budgets, retuning and constrained re-pointing

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

## Beam-error budgets (2026-10-17) *derived*

`mode_uncertainty` and the beam-errors tutorial: NIF hohlraum (10.1 × 5.75 mm, 3.1 mm
LEHs, 0.9 mm capsule, h = 0.45, capsule_h = 0.12), 48 quads at 8 TW, independent errors
per quad of 50 µm rms pointing (displacement length, transverse to the beam) and 3% rms
power.
- Per-beam sensitivities from three dual partials per beam agree with central
  differences of full evaluations to 2e-9 (relative).
- Linear propagation agrees with a 300-shot Monte Carlo of full evaluations to 1–2% in
  every mode l = 1–4 (the sampling error).
- Example pointing: rms of the mode amplitude (√Σ_m (δc_lm/c00)²) 0.28, 0.21, 0.13,
  0.03% for l = 1–4; pointing and power contribute about equally to l = 1 (0.20% each).
  ε₁ = mode_amplitude: median 0.43%, 90th percentile 0.67% (≈ 67 km/s of hot-spot
  velocity at MacGowan's 100 km/s per percent).
- After tuning the cones (l = 1–4 nulled at fixed drive): 0.23, 0.16, 0.09, 0.02%.
  The robust objective (expected mode power) started from the tuned design does not
  move it: with cone-level pointing, independent beam errors leave a floor set by the
  error sizes, not by where the cones point.

## Retuning to a reference drive (2026-10-17) *derived*

`retune_drive` and the retune tutorial: the same NIF hohlraum and mesh with
`WallLossAlbedo` (wall albedo 0.65 → 0.94 through the pulse), tutorial pointing, a
reference pulse from a total power (60 TW picket, 15 TW foot, 50 TW step, 420 TW peak)
and a cone fraction (0.35 in the foot to 0.25 at the peak, an engineering example, not a
shot), drive histories at 0.4 ns steps to 9.6 ns. Bundle 31 (Q31T, Q31B) dropped.
- Bundle down on the reference pulses: ε₁ 2.70% averaged over the peak (7.6–9.4 ns),
  5.17% at most (in the foot, where the dropped pair carries a larger share of the
  inner-cone drive); Tᵣ −3.2 eV at the peak. ε₂ and ε₄ barely change (10.20 → 10.00%,
  0.42 → 0.42%: the tutorial pointing is not tuned for P2).
- Retuned per quad, each within 1.3× its cone's reference power per quad: ε₁ 0.031%
  at the peak and 0.18% at most, ε₂, ε₄ and Tᵣ back to the reference's (ΔTᵣ +0.00 eV).
  The fit holds the albedos at the reference's; a drive history of the new pulses,
  with albedos following the drive, reproduces the predicted c_lm/c₀₀ to 8.4e-5, so the
  linearization is not the limit.
- The unconstrained retune is not unique and uses its bounds: 25 of 46 quads at 1.3×
  at the peak, some near zero. Weighting the departure from the reference powers
  (`stay`, per quad relative to the mean reference power):

  | stay | ε₁ peak | ε₁ max | powers / reference at 8.8 ns | quads at the limit |
  |---|---|---|---|---|
  | 0 | 0.031% | 0.18% | 0.00–1.30 | 25 |
  | 1e-5 | 0.13% | 0.26% | 0.78–1.30 | 7 |
  | 1e-4 | 0.34% | 0.46% | 0.93–1.30 | 3 |
  | 1e-3 | 1.31% | 1.57% | 1.01–1.19 | 0 |

  A dropped bundle's mode 1 is fixed mostly by the quads around it; a stay of order
  1e-4 keeps most of the correction with powers within a few percent of the reference
  for all but a few quads.
- The per-quad caps here are relative to the reference powers, not NIF's
  per-quad energy or power capability (which depends on the pulse and the laser's
  state). The view factors and pointing are the reference's through the pulse.
- Remeasured 2026-11-02 with the level held (next section) and the active-set polish;
  the 2026-10-17 table (soft level penalty) had ε₁ 0.023/0.12/0.34/1.24% at the peak.
  The numbers of the reference also moved (bundle-down ε₁ 2.02% then) with changes to
  the tutorial's pointing and pulse since; the trend in `stay` is unchanged.

## Retuning faces and foils; the held level (2026-11-02) *derived*

`retune_drive(quantity, ...)` with `CapsuleModes`, `FaceFlux` or `FoilIrradiance`; a
method change, no new constants.
- **Held level.** The level row was a penalty of weight `drive_weight` (default 1). A
  penalty pulls the level below its target by about the square of the nonuniformity
  (the shape terms scale with the level squared): negligible for capsule modes of a few
  percent, a percent on a face at 10% RMS, a third on a foil spot at 70% RMS (measured
  on the test foil, SG4 spots). Its target is now shifted by the secant method until the
  level is held (`drive_weight` removed). On `test/retune.jl`'s dropped-bundle fixture:
  free retune |drive − 1| 1.4e-5 → 2.9e-7, worst predicted ε₁ 3.5e-5 → 1.8e-5; with
  `stay = 1e-2` 1.4% → 5e-13 (ε₁ 3.62 → 3.63%); with `smooth = 0.1` 74% → 1e-8 (the
  penalty had let smoothing hold the powers near the first step's; ε₁ 1.2 → 3.2%).
  Without `stay` or `smooth` the 3ω capsule retune before the change was
  identical to that of the earlier code.
- **Solver.** L-BFGS-B stopped at a projected gradient near 2e-6 of the gradient's scale
  on these rank-deficient problems (46 quads, 25 rows); active-set Gauss–Newton steps
  after it take that to 1e-16–1e-10, and the held level then takes 3 solves per step.
- **Validation** (`test/retune.jl`): two beams symmetric by a half turn, the second at 0.6
  of the first's power in the reference, on a one-sided foil (SG 0.8 mm spots, 0.6 mm
  spot radius, RMS 14%) and on the closed bottom cap of a halfraum (albedo 0.8). The
  reference target recovers the powers to 2e-14 (foil); the uniform target gives equal
  powers at the asked level (1e-10 in the test, rounding measured), with the RMS a
  minimum along that level against ratios 0.9–1.1; reported mean, RMS and peak-to-mean
  equal those recomputed from depositions (foil) and from the history's flux on the
  face's elements (face) to 1e-12.

## Aperture keep-out rims (2026-11-03) *derived*

- **Before.** `ApertureKeepOut` took every hole as the circle through its mesh's
  outermost vertices about the area centroid, sampled at `nrim` = 96 points, and the
  clearance was the minimum over the samples. A rectangular, slotted or elliptical hole
  became its circumscribed circle (clearance overstated by up to the difference of the
  circle and the opening); a circle's clearance was overstated by the sampling, up to
  1.2e-4 mm on the 48 quads of `test/basics.jl`'s keep-out case (cones pointed through
  1.55 mm holes) against a 2e5-point reference on the circle.
- **Now.** The rim is the hole's boundary: edges of its triangles (with a coplanar lip
  about the same center) that no other triangle continues across. Shared edges cancel
  by their end points, and an edge with another triangle just beyond it (1e-6 of its
  length) is inside (T-junctions). On a `DiskSurface` hole, an edge whose ends lie at
  one radius from the center is the arc between them (the elements' exact sectors), and
  arcs of one circle with elements on both sides are inside (rings of different counts,
  whose chords leave a crack). The rim is sampled at about `nrim` points and the minimum
  refined by golden section (40 steps, Float64 abscissae, so ForwardDiff carries only
  the distance) on the four sample intervals with the smallest ends.
- **Validation** (`test/keepout.jl`, round spots so the envelope distance is exact):
  a 2 × 1 mm rectangle, beam down its normal: the closed-form edge distance to 1e-12
  relative, the gradient in aim (0, −1) to 1e-8 (the refined point's position); an
  oblique (35°), defocused f/8 beam: within 1e-9 mm of the minimum over 2e5 points on the
  analytic edges. An ellipse (1.5 × 0.8 mm) meshed as a 64-gon: the closed-form distance
  to the polygon's edges to 1e-12, and within the polygon's sag (1.8e-3 mm) of the true
  ellipse (6.8e-4 mm here); oblique within 1e-9 mm of the dense polygon. A circle on its
  `DiskSurface`, a 24-gon or a 16-gon disk in a 32-gon ring: within 1e-9 mm of 4e5
  points on the circle. The basics case agrees within 4e-12 mm with the circle reference.
- **Changed results.** Circular holes: clearances lower by the old sampling error (up to
  1.2e-4 mm); the keep-out penalty moves accordingly. Noncircular holes use the true rim in
  place of the circumscribed circle. Cost: 1.1 ms for 48 beams (about 2.7 times the old
  96 samples).

## Re-pointing studies on the public API (2026-11-18)

- **VisRad's NIF plates by cone** (`VISRAD_NIF_PLATES`, `visrad_nif_spot`). Rev3 Inner
  on 23.5°, 1.0 Rev1a Inner on 30°, Rev3 Outer on 44.5°, Rev1 Outer on 50°: the plates
  behind the measured CPP profiles of `nif_spot` (LMF phase-plate table, 2020-10-15).
  VisRad names the four among its NIF CPPs
  (VR/power_sources/import_nif_cmt_data_file.html: "Rev3 Inner Cone 2W CPP", "1.0 Rev1a
  Inner Cone 2W CPP", "Rev3 Outer Cone 2W CPP", "1.0 Rev1 Outer Cone 2W CPP"), but that
  page does not tie them to cones; the cone assignment is the LMF table's, and VisRad
  workspaces built from its NIF template carry the same four (§ VisRad deposition
  parity: hohlraum runs, "Spots"). A campaign default for VisRad parity, not a record of
  any shot's plates. **H** names and radii, **M** the cone assignment as VisRad's
  default.
- **Quad frame** (`nif_beams(; spot_frame = :quad)`): the turn of each beam's spot to its
  quad port's φ̂ (`_vr_reference_angle`), the rotation VisRad applies (§ VisRad
  beam-profile probe, within 0.02° on 48 plates). Whether the facility's spots follow it is not
  established. **H** for VisRad, **L** for the facility.
- **Measured-table fallback.** The tabulated VisRad profiles (`_VR_MEASURED`) were
  replaced by VisRad's closed-form law (§ VisRad's beam profile law), which applies to every
  spot, so no spot falls back to another profile and no warning applies.

## Capture and clearance rules as optimization quantities (2026-11-18) *derived*

Roadmap item 35. Each rule is a quantity of an evaluation (or of its beams) with a
natural limit, generic in the number type and thread-safe; the penalties are wrappers.
Validation in `test/design_rules.jl` on a mock-up: a coarse NIF-scale hohlraum (10.1 ×
2.875 mm, 1.55 mm holes) with a one-sided 3 × 2 mm shield plate beside the top hole,
facing away from the beams that pass it, and a closed 0.3 mm radius, 4 mm long stalk
(12 facets) across the hole; the 48 NIF quads at the cone aims of the basics test.

- **Penalty terms.** `keepout_penalty`, `keepout_clearance` and `direct_fraction` are
  the terms of the objectives (`_penalties` equals the public fractions exactly).
- **Power on named parts.** `power_fraction(ev, names...)` (absorbed on the parts'
  elements) plus `clipped_fraction(ev, names...)` (their back sides) is all the light the
  parts receive. Against `component_power`: 1e-12 relative; with the missed light, the
  parts' sum closes to 9e-16. Against the reference ray trace (nspot = 121, both sides):
  −2.0e-3 on the shield (4.4% of the delivered power), −7.7e-4 on the stalk (4.7%),
  the operator's sampling of the spot edge. Gradient in (z, r) per cone against central
  differences (h = 1e-5 mm): 9.5e-10 and 4.1e-10 relative.
- **Keep-out for named hardware** (`HardwareKeepOut`, `hardware_clearance`). Each
  beam's worst-case clearance from named parts, as `clearance` judges it with
  `wings = []`: the exact point–ellipse distance to the best-focus contour (1% of
  peak by default, the `ApertureKeepOut` default) grown by the f-cone, minimized over
  the parts upstream of the nominal axis's landing, less `pointing + target_offset +
  hardware_tolerance + focus/(2F)`. On a triangle on one side of best focus the
  distance is convex and has no interior stationary point off the medial strip (the
  major axis's segment ±(A² − B²)/A swept along the beam, at least B²/A inside the
  envelope), so the minimum lies on segments: edges, crossings with the strip, cuts by
  the landing plane, split at best focus. Best-first search over clusters of 16
  triangles and their segments by a lower bound (circumscribed circle and homothetic
  ellipse, Lipschitz √(1 + (1/2F)²)), each segment by Illinois regula falsi on the
  derivative: the exact minimum, not a sample. Derivatives by Danskin's theorem at the
  closest point; on the medial strip the mixture of the two closest points' that
  leaves the minimum stationary (minimax saddle).
- **Validation** (`test/keepout.jl`, `test/design_rules.jl`). Closed forms (plates normal to
  and along a beam, a plate cut by the landing plane, the center of an elliptical spot on a
  plate through the axis with its gradient [0, 0, 1/2F, 0, −ρ]): 1e-12 relative. An oblique
  (35°), defocused f/8 beam beside a tilted plate: within 1e-8 mm of 1.6e5 points on its
  rim. Against `clearance` (`rtol = 0`, resolution 1e-7 mm) on the mock-up's 48 quads:
  within 1e-6 mm (7.6e-7 measured) wherever −0.5 < clearance < 1 mm; the keep-out is lower
  (conservative) on 5 of 24 shield beams by 3e-4 to 0.40 mm, all above 1.2 mm, where
  `clearance` skips shield points the cap hides from the beam (the keep-out, as
  `talis_check`'s 3ω clearance, counts them). In deep violations (clearance < −0.75 mm) the
  report gives values up to 6e-5 mm lower: its point–ellipse distance (`_ellipse_distance`,
  Eberly's bisection) loses digits within ~1e-9 mm of the major axis inside the medial
  segment (2.4e-5 mm at a point 3e-13 mm off the axis, against the keep-out's value, which a
  1e6-point grid about the point confirms to 1e-10). *Since 2026-11-18* (§ Optimization
  defaults decided) the distance is exact there too and the two agree within 1e-7 mm on
  every row below 1 mm. Derivatives in aim, defocus and spot axes against central
  differences: within 1e-8 on every beam below 3 mm, medial-strip minima included; not where
  the minimum lies on the landing plane's cut (the cut does not move in the derivatives: one
  bottom beam at 10 mm from the stalk hidden behind the hohlraum).
- **Cost** (48 quads, one thread): 0.43 ms for the 48-triangle shield, 0.78 ms for the
  96-triangle stalk, 3.3 ms for a 3840-triangle stalk (the same with 12-partial duals:
  the search runs on values, the duals enter at the closest point only); the aperture
  keep-out 0.5 ms and one deposition of the coarse mock-up 48 ms.
- **Discontinuities.** The clearance is a minimum over a fixed surface, continuous in
  the beams; it jumps only where the nominal axis's landing moves from one component to
  another while the parts reach past the landing plane.
- **Toward the TaLIS 3ω clearance** (measured, not implemented). On the mock-up,
  `talis_check(...; capture = (:wall, :capsule), hardware = (:shield, :stalk))`'s
  `:clearance_3w` worst values (resolution 1e-7, `rtol = 0`) equal the hardware
  keep-out's clearance less the rule tolerance r tan 0.5° (r from TCC) taken at the
  keep-out's closest point to within 9.6e-5 mm (shield) and 7.2e-4 mm (stalk), on all
  48 rows: the two judge the same points (no reachability). Taken a posteriori the
  tolerance jumps where the closest point does (by tan 0.5° times the change in r); the
  exact rule needs the tolerance term inside the minimization (a concave term, so the
  per-segment convexity and the medial-strip argument hold only approximately).
- **TaLIS 3ω clearance as a keep-out** (`HardwareKeepOut(...; tolerance, pivot, miss)`). The
  rule tolerance −r tan α (r from the pivot: TCC, the aim point, the larger of the two, or a
  point) is subtracted inside the minimization. It is concave, and off the medial strip the
  envelope distance's Hessian has rank one, so in every plane the sum curves down along some
  direction: the minimum stays on the segments (edges, strip crossings, landing cuts). Along
  a segment the sum is not convex; each segment is minimized with the term replaced by its
  chord (a lower bound, exact at the ends: the convex search of the keep-out applies), and
  split at the chord's minimizer while the term there lies more than 1e-14 mm above the
  chord; cluster and segment bounds subtract tan α (r + radius). `miss = :all` judges a beam
  whose axis meets no target along its whole length (TaLIS: "everything along the beam is
  checked"), and a beam landing on one of the parts of the keep-out gets `Inf`
  (`talis_check` skips the landed component). With `targets = capture`, `tolerance = 0.5`,
  `miss = :all` and the defaults otherwise, on the mock-up's 48 quads against
  `talis_check(...; capture = (:wall, :capsule), hardware = (:shield, :stalk), omega1 =
  false, resolution = 1e-7, rtol = 0)`, every `:clearance_3w` row (24 shield and 25 stalk
  rows finite, the rest `Inf` in both): the keep-out is 6e-8 to 1.1e-7 mm above the report's
  lower bound on the 43 rows whose worst point lies off the major axis. On the 6 rows whose
  worst point lies within 3e-10 mm of the major axis inside the medial segment (−0.82 to
  −1.11 mm, deep violations), it is 1.7e-7 to 6.4e-6 mm above: there `_ellipse_distance`
  returns values that are low (`clearance.jl`, a known loss of digits left to the user), and
  the closed form on the axis, −B √(1 − u²/(A² − B²)), at the worst point of the report lies
  7e-9 to 2.2e-7 mm above the keep-out's minimum (the keep-out does not miss a lower point).
  *Since 2026-11-18* (§ Optimization defaults decided) those 6 rows (now 4, the worst points
  having moved) agree within 1.0e-7 mm as well. Continuity (closed form, two plates beside a
  beam down −z, their closest points at r = √10 and √82 mm from TCC): the keep-out equals
  the closed form min over the plates of (far edge distance − cone · axial − tan α r) to
  1e-12 relative on an aim scan of ±0.1 mm through the switch of the closest point, with no
  step larger than the aim step (1e-3 mm); the tolerance taken a posteriori at the
  tolerance-free closest point steps by 0.050 mm there. Derivatives in aim, defocus and spot
  axes against central differences (h = 1e-6): within 9.4e-8 on the beams below 3 mm for
  pivots TCC, aim, max and a point (1e-8 without the tolerance; the chord split's 1e-14 mm
  stop sets the FD noise). Cost (48 quads, one thread): 0.64 ms for the shield and 0.83 ms
  for the stalk (0.42 and 0.74 ms without the tolerance), 3.2 ms for a 3840-triangle stalk;
  the same with 8-partial duals. `talis_check` on the same mock-up and settings: 27 s on one
  thread. With `tolerance = 0` (the default) every existing result is unchanged (keep-out
  values and gradients, `clearance` and `talis_check` reports on the mock-up).
- **TaLIS 3ω capture as a quantity** (`CaptureKeepIn`, `capture_margin`). The part
  landed on is the first of the capture parts the nominal axis meets; the margin is the
  envelope distance at the 10% contour minimized over that part's free edges (as
  `talis_check` takes them), less the error set and the rule tolerance, by the keep-out's
  search (edges only, no landing cut). Defaults: the TaLIS values (10%, 0.5 mm, 0.5°
  about TCC). Where the axis misses every part `talis_check` reports −Inf; the quantity
  continues the margin there as the margin from the outline nearest the axis less twice
  the axis's transverse distance d from it: continuous where the axis crosses the
  outline (d = 0), falling at least as fast as the axis moves away (the margin is
  1-Lipschitz, −2d falls at rate 2), and on a straight edge the linear continuation of
  the landed margin. Design choice: −Inf would make a penalty infinite (and
  `keepout_penalty` skips non-finite values), and no gradient would bring the beam back.
  `-Inf` remains only when no capture part has a free edge. Checks: closed form (square
  foil normal to a round beam, aim across its edge, landed and missed): 2e-15 relative.
  Mock-up (a 5 × 4 mm one-sided foil at TCC tilted 20°, the 24 upper NIF quads aimed on
  a grid across and past its edges) against `talis_check` (resolution 1e-7, `rtol = 0`):
  the 13 landed rows 4.7e-8 to 1.0e-7 mm above the report's lower bound; the 11 missed
  rows (report −Inf) finite, −0.79 to −1.44 mm. An aim scan through the foil's end: slope
  −0.40 to 0.00 inside, −1.83 to −0.96 outside, no step at the crossing. Derivatives
  against central differences, landed and missed: 4.2e-8. Cost: 0.37 ms for 24 beams on
  the 50-triangle foil (one thread).
- **Refactor.** The keep-out's evaluation at the closest point is shared with the
  capture margin (`_envelope_clearance`). It changes `hardware_clearance` by at most
  4.4e-16 mm (on 2 to 12 of 48 beams per mock-up keep-out) and the penalty gradients by
  at most 2.3e-13 absolute (1e-16 relative): rounding of the same arithmetic compiled in
  another function, measured against the commit before.
- **TaLIS 3ω incidence as a quantity** (`IncidenceRule`, `incidence_angle`). As
  `talis_check` judges it: the landing part cut into 4^level sub-triangles with their
  normals (exact on analytic surfaces), the largest angle to the beam axis among those
  visible from upstream that the worst-case 10% footprint reaches, plus the f-cone's
  half-angle. The set reached steps as the beam moves; the quantity takes each
  sub-triangle's angle less `sharpness × max(0, d)` (d its exact clearance from the
  footprint, by the keep-out's search on the one sub-triangle), a continuous maximum that
  counts a sub-triangle in full once reached and fades it over Δθ/sharpness mm outside.
  `sharpness = 100°/mm` is an engineering default (a sub-triangle 0.01 mm outside the
  footprint counts 1° less). Best first over the sub-triangles by an upper bound (their
  angle less sharpness × a lower bound of d). Checks: closed form (a part of a plate normal
  to a round beam and a 60° plate beside it, with the clearance of the second plate from
  −0.1 to 0.7 mm): 1e-12 relative, the derivative in the aim 100°/mm to 1e-9; it lies at or
  above `talis_check`'s step and equals it outside the 0.6 mm fading band. Against
  `talis_check` on the design-rules hohlraum (wall and capsule, exact cylinder and sphere
  normals, 48 quads, 8 of them above 70°) and on the capture foil (13 landed rows): equal on
  every row (on each the largest angle lies inside the footprint, where no fading applies:
  the quantity is locally constant in the pointing, as the rule is). Cost: 12 ms for the 48
  quads on the hohlraum's 6976 sub-triangles, 0.63 ms for 24 beams on the foil's 640 (one
  thread).
- **TaLIS 1ω clearance as a quantity** (`Omega1KeepOut`, `omega1_clearance`). The
  footprint of `talis_check` (per NIF beamlet a square, near edge 4.45 mm from the 3ω
  axis, tapering to the 1ω focus, the focus error as ±1 mm shifts), less the error
  radius and the rule tolerance. The square distance is convex in the point but the
  minimum over a triangle can lie inside it (the square's medial lines), so the minimum
  is not found by a segment search: `talis_check`'s Lipschitz branch and bound runs to
  `resolution` (1e-9 mm by default, `rtol = 0`) and the quantity is the value at the
  best point found (within `resolution` above the minimum, continuous to that accuracy),
  its derivatives the footprint distance's there (each beamlet's focus moves with the aim
  and defocus). Checks: closed form (one beamlet's prism footprint and a plate between it
  and the beam): 1e-9 relative. Against `talis_check` (resolution 1e-7, `rtol = 0`) on
  the mock-up's 48 quads, shield and stalk (worst −0.25 to 3.8 mm): 9.5e-8 to 1.0e-7 mm
  above the report's lower bound on all 96 rows. An aim scan finds no step beyond the
  Lipschitz bound. Derivatives against central differences (h = 1e-4, resolution 1e-11):
  3e-11. Cost (48 quads, one thread): 160 ms (shield) and 184 ms (stalk) at 1e-9, 91 and
  124 ms at 1e-7: three times one deposition of the coarse mock-up, the price of the
  branch and bound.
- **Not represented.** 1ω management (curvature or specular reflection, a review), the
  diagnostics-to-1ω rule, the TAS stay-in zone (not a function of the pointing) and the
  LEH-edge review tier. The loop that reruns `talis_check` and tightens minimums is out
  of scope.
- **The snap removed** (2026-11-18, later; *measured*). The keep-outs took points within
  1e-6 B of the major axis onto it (the snap of § Optimization defaults decided), so that
  their medial-strip minima were known as such; next to a center of curvature the
  closest point moves fast along the ellipse, and the snapped distance read up to 1e-8
  of a semi-axis low on a 0.40/0.39 spot and 5.1e-7 on 0.4/0.399999 (against the
  BigFloat minimum on 7,020 points, seven ellipses; the exact path: 6.6e-16). On a
  `HardwareKeepOut` whose nearest point is a triangle's vertex 9.99e-7 B off the axis at
  the center of curvature of a 0.4/0.3999 spot, the keep-out was 1.02e-7 mm low. Now
  `_ellipse_closest` has no snap (the exact path of `clearance` is unchanged) and a
  keep-out's minimum takes the medial derivative where the closest point's mirror image
  across the major axis is as near to within 1e-10 of the semi-major axis, inside the
  medial segment (`_medial_tie`): the minima the search finds on the strip lie within
  3.1e-14 mm of the axis (all 716 such evaluations of `test/design_rules.jl`, none
  between that and 1e-3 B). The mixture of the two points' derivatives, stationary along
  the segment, keeps its value (the value at the kink). Measured against the base commit
  on the mock-up (48 quads; `HardwareKeepOut` plain and with the TaLIS tolerance about
  TCC and the aim, `CaptureKeepIn`, `IncidenceRule`, `Omega1KeepOut`): values within
  2.2e-16 mm (233 of 243 finite values identical), derivatives in aim, defocus and spot
  axes within 5.5e-11; their agreement with central differences unchanged (7.5e-10 to
  9.4e-8), and with `clearance` and `talis_check` unchanged to the printed digits
  (keep-out − report 5.8e-8 to 1.06e-7 mm off the hidden shield points). The vertex case
  now agrees with `envelope_distance` exactly (0 mm; `test/design_rules.jl` checks 1e-12).

## Constrained re-pointing: an irradiance ceiling and mode gates (2026-11-18) *measured*

Roadmap item 36. `optimize_pointing(...; constraints)` takes `Constraint`s (ceilings and
floors on any function of an `Evaluation`) and solves them by an augmented Lagrangian
around the same L-BFGS-B; `irradiance_ceiling` and `mode_gate` are the ready-made ones,
`drive_objective` the objective they go with. The public case is the dropped-bundle tutorial
(bundle 11 freed, the other 46 quads re-pointed in `:phi` and `:z` within ±0.3 mm, the
tutorial's mesh, 1,728 wall elements, and its coarse deposition operator), the largest
drive with each degree's departure from the 48-quad shape within the gate
`√(2l+1) × rms_mode` of `mode_uncertainty` at 50 µm rms pointing per quad (0.346,
0.395, 0.313, 0.095% for l = 1–4; the start departs 1.09, 1.07, 0.77, 0.12%), under
the 48-quad pointing's peak wall irradiance (8.273 TW/mm²; the start's hottest element
is the same, so the ceiling is active from the start). Scratch scripts, not in the
repository; 8 threads.

**Method** (`src/constraints.jl`). Powell–Hestenes–Rockafellar augmented Lagrangian for
inequalities on the normalized violations `gᵢ = ±(value − limit)/scale`, with
ALGENCAN's updates (Andreani, Birgin, Martínez and Schuverdt, SIAM J. Optim. 18, 1286,
2008): λ ← max(0, λ + μg); μ × 10 when `V = max |max(gᵢ, −λᵢ/μ)|` has not halved; first
μ = 10 max(1, |f|)/max(1, ½Σ max(0, g)²), clamped to [1e-6, 1e6]; infeasible once μ has
grown 1e8-fold. The inner runs are `_lbfgsb` (the unconstrained path, unchanged) warm-started,
projected-gradient tolerance from 1e-3 tightened tenfold per outer iteration to
`g_tol`, and a relative value tolerance of 1e-10 (Optim's `f_reltol`; engineering
default). Without it each inner run on the drive objective crept for hundreds of
iterations at 1e-6 gains: the first version ran out of its 1,000 iterations in 2–3 outer
iterations on every start (constraints met to 2e-4, but never "converged").

**Smooth stand-in for the peak.** The area-weighted p-norm `(Σ A I^p / Σ A)^(1/p)`,
scaled to equal the peak at the reference, against the squared hinge over the
elements above the limit (`Σ (A/Ā) max(0, I/L − 1)²`, limit 0, scale 1e-4). Augmented
Lagrangian runs to 1,000 iterations from the start and from one random start in the box
(peak over the elements after the run, relative to the limit):

| ceiling | start 1 | start 2 | gates |
|---|---|---|---|
| none | 1.272 | 1.383 | met |
| p = 16 | 1.029 | 1.027 | met |
| p = 32 | 1.010 | 1.005 | met |
| p = 64 | 1.0007 | 0.9991 | met |
| hinge | 1.0001 | 1.0006 | start 2 unmet (l = 2 and 4 at 2.0 and 1.7 × gate) |

The hinge holds the peak exactly where it converges but fails on the second start: its
gradient vanishes at the limit (a degenerate constraint for the multiplier update), so
it acts as a pure penalty. Chosen: p = 64 (`irradiance_norm` default).

Mesh: 41 designs near the limit (the start and random shifts of ±0.01–0.08 mm), wall
elements of 0.2, 0.3 and 0.45 mm (4,896, 2,176, 1,104 on the wall; default deposition
resolution, `reach = 0.6`). The 48-quad peak is 8.84, 8.67 and 8.15 TW/mm²: element
averages lower the peak 8% from 0.2 to 0.45 mm; the 64-norm moves 5% (8.03 → 7.61), the
32-norm 3.5%, the 16-norm 2.5%. Ratios to the reference on the same mesh, rms (max)
difference between meshes:

| quantity | 0.2 vs 0.45 mm | 0.2 vs 0.3 mm | correlation with the peak ratio (0.2 mm) | verdict ≤ 1 agrees with the peak's |
|---|---|---|---|---|
| peak | 0.0128 (0.034) | 0.0113 (0.035) | 1 | 41/41 across meshes |
| p = 16 | 0.0066 (0.025) | 0.0007 (0.002) | 0.930 | 19/41 |
| p = 32 | 0.0108 (0.041) | 0.0030 (0.008) | 0.969 | 31/41 |
| p = 64 | 0.0102 (0.038) | 0.0067 (0.018) | 0.991 | 37/41 |

Small random shifts raise the peak 1–2% (a spot moves onto a neighbor's) while the
norms barely move, so lower p is more mesh-steady but blind to the peak; p = 64 tracks
it (0.99) and its ratio is as mesh-steady as the peak's. Far from the limit (37 designs
shifted ±0.15–0.3 mm, peak ratios to 1.66) every quantity differs between 0.2 and
0.45 mm by about 6% rms (the coarse mesh averages the stacked spots). Confidence:
measured on this hohlraum and these spots only.

**Method comparison** (p = 64 ceiling, three starts: the tutorial's and two random in
the box; budget 2,000 evaluations for the penalty methods; the augmented Lagrangian's
counts are from before the convergence rule of § Constrained re-pointing: convergence
on a feasible, active start, kept as dated history):

| method | drive (× 48-quad) | departure / gate, l = 1–4 | peak / limit | evaluations |
|---|---|---|---|---|
| augmented Lagrangian (one call) | 0.96259–0.96261 | 1.000 1.000 1.000 1.000 | 0.9996–1.0025 | 488–765 (converged, `:box`) |
| penalty, final gates, one call (the study's form: 100 × squared excess over each gate², 0.01 × each ratio², 100 × squared ceiling excess) | 0.95921 | 0.002 0.120 0.005 0.066 | 1.0088–1.0095 | 1,071–1,081 (1,000 iterations) |
| the same, staged by hand (gates × 4, 2, 1.4, 1; 250 iterations each) | 0.95921 | 0.002 0.120 0.005 0.066 | 1.0084–1.0095 | 1,234–1,437 |
| penalty without the 0.01 term, one call | 0.95747–0.96153 | gates met | 1.0012–1.0199 | 2,000 (budget: 67–672 iterations, line searches stall on the stiff hinge) |
| penalty-weight continuation (½ w Σ max(0, g)², w = 1, 10, …), inner runs as the study's | 0.96261 | 1.000 (violation 2.4–3.0e-4) | 1.0008–1.0026 | 2,000 (budget; the first stage crept) |
| the same continuation, inner runs with the augmented Lagrangian's value tolerance | 0.96258–0.96260 | ≤ 1.000 (violation ≤ 2.4e-4) | 1.0003–1.0072 | 725–1,262 (2–4 stages) |

All methods reach the gates here: the public case does not reproduce the study's stall
at the bounds. What it shows is the cost of penalties: with the symmetry reward the
optimum sits deep inside the gates at 0.35% less drive, staging changes nothing but
the evaluations, and without the reward the stiff hinges stall the line search. The
augmented Lagrangian reaches the constrained optimum (every gate active) from every
start and reports it with its multipliers; under the later convergence rule (which
declares convergence only at the final inner tolerance) it takes 997 evaluations on
the tutorial's start against the penalty form's 1,072, a comparable count (the half
measured before that rule no longer holds), for the same 0.35% more drive.
Continuation in the weights, given the same inner tolerance, reaches the same design
in 1.0–1.7 times the augmented Lagrangian's evaluations on the same start (838, 725,
1,262 against 765, 698, 488), with no multiplier estimates and a schedule to choose;
continuation in the gates is the hand staging above. Chosen: the augmented Lagrangian.
Confidence: one public case, three starts; the stall seen in the study was not reproduced.

**Checks.** An equality-active quadratic program (48 quad powers fitted to 1.1 × the
drive's coefficients through l = 2, regularized, under a total-power ceiling): the run
matches the KKT solution to 4e-6 relative and its multiplier to 1e-5 (tested to 1e-4 and
1e-3); a box that keeps the powers above the ceiling ends `:infeasible_box` after 10
outer iterations. Constraint values' gradients against central differences: 1e-8–1e-9
relative. Unconstrained runs are identical to those of the base commit (three tutorial
runs: x, value, iterations, evaluations, reason, at-bounds).

**Bubble drivers.** The ceiling at the reference's peak irradiance is an engineering
limit of the study, not a published bubble criterion (the published scalings use the
outer cones' picket fluence on the wall: § The swing and the bubble without new data,
part 4, and `reference/litreview/bubble_drivers.md`): `weight` (a function of the
beam multiplying its deposited power; 0 leaves it out) gives the outer cones' or a
picket's fluence map on the same machinery, at one more deposit per evaluation.

**Keep-outs and the light on parts as constraints** (roadmap items 35 and 36 joined).
`clearance_floor(keepout; minimum = keepout.minimum, scale = 1 mm)` holds each beam's
`keepout_clearance` at or above the minimum: one component per beam of the evaluation,
labeled by its name (a `Constraint`'s labels may now be a function of the evaluation,
so the components follow the beams and one constraint applies to a stage on quads and one
on their beams). Any keep-out with a `keepout_clearance` method and a `minimum` field
takes part: `ApertureKeepOut`, `SphereKeepOut`, `HardwareKeepOut`, `CaptureKeepIn`,
`IncidenceRule`, `Omega1KeepOut`. A beam the keep-out does not apply to (`Inf`) has the
violation −∞: its term and multiplier stay zero. The scale of 1 mm puts the default
`feasibility = 1e-3` at 1 µm, finer than the keep-outs' agreement with `clearance`
(engineering default). `light_ceiling` and `light_floor` hold `light_fraction =
power_fraction + clipped_fraction` (all the light on named parts, from either side, a
fraction of the delivered power) at most or at least a limit, given or taken from a
reference evaluation; the scale `max(limit, 1e-3)` is relative above 1e-3 of the
delivered power and absolute below it, so a ceiling at zero is met to 1e-6 of the
delivered power (engineering default; a reference's light on a part meant to stay dark
is often 1e-13 or 0, where a relative scale fails). On the tutorial's hohlraum the
nominal light on the caps is 3.8e-13 of the delivered power for the 48 quads and 8e-9
for the design within limits whose worst case reaches 58 µm into the caps: the nominal
light misses what the worst case sees, so there the keep-out is the binding limit.

Per beam or the smallest clearance? The tutorial's first stage (46 quads, the inner
ones held, a power scale per cone in [0.8, 1], the four gates, the wall ceiling, the
cap keep-out at 10 µm and the light on the caps at most 1e-5), with the keep-out as the
per-beam vector, as the minimum over the beams (one component, kinked where the closest
beam changes) and as a log-sum-exp soft minimum with τ = 5 µm (smooth, below the
minimum by up to τ ln 46 = 19 µm); start 1 is the tutorial's, start 2 a random point in
the box; budget 3,000 evaluations:

| keep-out as | start | evaluations (outer) | drive (× 48 quads) | least clearance (mm) | departure / gate, l = 1–4 |
|---|---|---|---|---|---|
| per beam | 1 | 885 (7) | 0.95771 | 0.0099 | 1.000 1.000 0.980 1.000 |
| per beam | 2 | 745 (7) | 0.95774 | 0.0107 | 0.999 1.000 1.000 1.000 |
| minimum | 1 | 807 (7) | 0.95776 | 0.0100 | 0.998 1.000 0.990 1.000 |
| minimum | 2 | 1,154 (7) | 0.95779 | 0.0101 | 1.000 1.000 1.000 1.000 |
| soft minimum, τ = 5 µm (earlier stopping rule) | 1 | 715 (4) | 0.95778 | 0.0186 | 0.998 1.000 0.991 1.000 |

Every run ends `:box` with no limit violated. Per beam and minimum reach the same design
(drive within 1e-4) at no consistent difference in cost (885 and 745 against 807 and
1,154 evaluations); the soft minimum keeps 8.6 µm more clearance than asked, its bias.
In the per-beam polish of the second stage (184 beams, 200 evaluations) the keep-out is
not binding (21 µm) and the two forms give the same design to five digits.

Chosen: per beam, on grounds other than cost, which does not separate them here. Each
component is as smooth as its beam's clearance and carries a separate multiplier (zero
for beams clear of the minimum, which then add nothing to the inner problem), so the
report names the beams that bind; the soft minimum's bias wastes clearance, and the
hard minimum is a nonsmooth constraint that the augmented Lagrangian's convergence
theory does not cover. Confidence: one public case, two starts. A constraint of more than eight components displays
as its worst row and a count by status (a per-beam keep-out on 184 beams printed 184
rows). Tests (`test/constraints.jl`): the values against the quantities (the light
against `component_power`, 1e-12), gradients against central differences (1e-6), and
with one free parameter (the 50° cone's aim, pulled into the caps by the objective) a
constrained run ends with the keep-out, and separately the light ceiling, active at the
aim where the constraint reaches its limit, found by bisection, to 1e-5 mm.

**Stages: two calls, not a schedule** (roadmap item 34). The tutorial's study written
both ways, quads re-pointed, then each beam polished: (a) two calls of
`optimize_pointing`, the second on `remake(scenario; parameterization =
GroupPointing(nif_beams(result.beams); by = :beam, relative = true), deposition_kw)`
with the same constraints; (b) a schedule of stages (a parameterization built from the
previous beams, a box, budgets) run by one function that keeps each scenario and
result. (b) needs every stage-dependent piece as a closure of the previous beams (the
parameterization, the held set re-derived from the beams), the deposition options
anyway (a scenario does not record them) and a pair (scenario, result) per stage,
since `evaluate`, `check_constraints` and `pointing_summary` take the scenario; it
does not save a line of the tutorial's case and hides the check between the stages that
matters most (below). Chosen: (a), no new API; the change it needed was constraints
whose components follow the beams. Not built, so not measured: multipliers carried from
one stage to the next.

The step from quads to beams changes the deposition model, not only the parameters:
each beam enters through its port at f/20.7 instead of the quad's f/7.9 cone. On
the tutorial's case the first stage's aims, deposited per beam, move the drive from
0.9577 to 0.9568 of the 48 quads' and mode 1 to 1.016 times its gate; the polish (240
free parameters, ±0.1 mm per beam, 120 evaluations, ended by the budget) brings it
back within the gate at a drive of 0.9584, every limit met, the clearance check's worst
case 21 µm from the caps. Under the earlier stopping rule the same polish ended after
one iteration and three evaluations, at the first inner tolerance (1e-3 on the
projected gradient, met at once where each of 240 parameters moves little of the
drive), with the drive unchanged (0.95685 against 0.95688 at its start): found here
independently of the next section, whose fix (320335b) the tutorial runs on; with a
stopping rule that also accepted an inner run ended by its value tolerance (the rule
tried here first) the polish converged in 97 evaluations to the same drive (0.95842)
and the earlier run of the section kept its 765 evaluations, unchanged. The tutorial as a
script on 8 threads: 249 s, of which the study adds about 100 s (141 s before the
convergence fix and the study).

### The exact peak, restarts and the box (2026-11-18, item 36 continued) *measured*

What a trial on a staged re-pointing study kept outside the repository still lacked
after the constrained run had met every gate and the ceiling in one call: the exact
peak (the ceiling bounds the 64-norm, so the peak ended a fraction of a percent over),
a restart that keeps the multipliers, a wall-clock budget with a separate outcome, the
cost of the box, the drive and departures in the summary line, and a check that the
ceiling's reference shares the scenario's laser–plasma model. Public measurements on
the dropped-bundle test case of `test/constraints.jl` (bundle 11 freed, 46 quads in
`:phi` and `:z` within ±0.3 mm, 0.45 mm wall elements, coarse deposition operator, 8
threads) and its closed-form quadratic program.

**Correction: the stand-in's margin.** "Up to 1% at p = 64" above holds on the
tutorial's mesh (0.9991–1.0025 of the limit). On the test case the same run from the
same start ends at different designs that meet the same four gates, depending on the
floating-point rounding: plain (`julia`), the 64-norm at 0.974 of its limit and the
peak at 0.9974 (drive 0.96259136, 1,044 evaluations, `:interior`); under `Pkg.test`'s
`--check-bounds=yes`, the norm at its limit and the peak at 1.030 (the full suite's
failure of `peak ≤ 1.02 × limit`, 8.3958 against 8.1523 TW/mm², 12 threads). A plain
driver with `--check-bounds=yes` on 8 threads reaches a third (peak 0.9909, norm 0.976,
drive 0.96251211, 742 evaluations, `:interior`; `exact = true` there changes nothing).
The drives differ by 8e-5 relative, more than rounding: the gates leave a flat family of
designs, and which one a run ends at depends on rounding and thread layout. The stand-in
by itself does not bound the peak to a fixed margin; the test now holds the peak with
`exact = true`.

**Exact peak** (`irradiance_ceiling(...; exact = true)`). After the run has converged
(feasible, complementary, inner run at `g_tol`), the peak over the elements is compared
with the limit; while it exceeds the limit by more than `feasibility`, the norm's
working limit is multiplied by limit/peak and the run continues from the same design,
multipliers and penalty (the penalty is not raised for the new limit's violation: the
stall measure restarts). The norm-to-peak ratio changes little near the limit, so one
or two rounds land the peak within `feasibility`. Only a peak above the limit lowers
the working limit; a peak below it is left (conservative; raising the working limit to
make the peak active would buy the drive that the margin holds, and risks cycling where
the peak moves between elements). On the test case with the ceiling active (peak over
its limit after the plain run):

| start, limit | plain: peak / limit, outer, evaluations, drive | `exact = true` |
|---|---|---|
| tutorial's, 0.98 × the 48-quad peak | 1.0139, 7, 771, 0.96224 | 0.9961, 8, 842, 0.96224 |
| random in the box, the 48-quad peak | 1.0288, 7, 1,306, 0.96260 | 0.9872, 8, 1,363, 0.96259 |

One lowering round each, 4–9% more evaluations, the drive the same to 1e-5. The peak
ends below the limit, not on it (the norm-to-peak ratio moved with the design),
and one or two gates leave their limits by up to 3% of the gate. Default kept at
`exact = false` (today's behavior), as a decision for the user; the recommendation is `true`,
where the literal criterion is the peak over the elements. *2026-11-18, later:* the user
decided `exact = true` as the default (§ Optimization defaults decided).

**Warm restarts** (`from = result`). The result's `state` carries the multipliers, the
penalty parameter and its ceiling, the inner tolerance, the last infeasibility–
complementarity measure and the working limits. On the quadratic program (58 evaluations, 8
outer iterations uninterrupted): split by `outer = 2`, the continued run is identical to the
uninterrupted one (each inner run restarts its L-BFGS-B anyway); split inside an inner run
by the evaluation budget, the continued run loses that inner run's L-BFGS-B memory and
reaches the same optimum to 0.6–1.6e-7 relative (multiplier to 2e-5–3.5e-4), in 55, 24 and
20 evaluations after splits at 10, 20 and 40, against 49, 28 and 28 from `x` only, with the
multipliers relearned. On the dropped-bundle test case (the uninterrupted run took 1,044
evaluations and 7 outer iterations, drive 0.96259), split at 400 evaluations (3 outer
iterations): the restart with `from` converges in 311 more (5 outer iterations), from `x`
only in 394 (7); both end at the same design (drive 0.96247, the ceiling's peak at 1.0215 of
its limit, every gate within 0.6%), 1.6% in `x` from the uninterrupted run's: the flat
family of designs above, so a split run reaches the same design as the uninterrupted one
only where the optimum is unique (the quadratic program).

**Time limit** (`time_limit`, seconds). Checked after each evaluation, after the
evaluation budget; `reason` and `outcome` `:time_limit`, the best design kept, the state
kept for `from`.

**Bound sensitivities** (`bound_gradient`, `bound_sensitivity`). For a parameter at a bound,
the optimum's derivative with respect to the bound is the Lagrangian's partial derivative
there (envelope theorem, active set held), `∂L/∂xᵢ` with the multipliers of the normalized
violations; `∂f/∂xᵢ` is that of the objective only. Check: on the quadratic program, one
power's upper bound 0.5 below its optimum, against central differences of the re-optimized
objective in the bound (δ = 0.05 and 0.01): 1.2e-5 relative without constraints (limited by
the convergence of the runs, since it is the same at both δ), 2e-6–4e-5 with the total-power
ceiling active (the derivative of the objective only is there 1.76 times the optimum's). One
gradient each after the run, only when a parameter is at a bound.

**Summary line.** `pointing_summary(...; drive_reference)` (an evaluation or its
`flux_ylm` coefficients) adds the drive ratio (c₀₀ over the reference's; equal to the
ratio of `mean_flux` to 1e-10) and the departure per degree.

**Laser–plasma model of the reference.** A ceiling whose reference was evaluated
without the scenario's CBET model bounds another quantity (its peak and the norm's scale
come from the reference); on the study it read far above its limit at the start, with no
warning. The ceiling now records the reference scenario's `lpi` and time and throws an
`ArgumentError` on an evaluation with another model (compared field by field, so a model
built alike passes) or the same model at another time. An error, not a warning: the
mismatch has no use, and a warning inside an optimization repeats per evaluation or is
lost after the first. `mode_gate` takes coefficient vectors, which record no model: the
same trap there is open.

**On the study** (scratch scripts, numbers not recorded here): one call with
`exact = true` from the as-shot start, under a time limit, converged with every gate
active and the exact peak at (just below) the reference's, two parameters held by the
box, within the wall time its earlier norm-only run had been cut at plus 40%. Its
bound sensitivities failed their check there: restarts with `from` in the same box and
in a box widened at the two bound parameters both left the bounds, and the widened
box gained over an order of magnitude less than the first-order prediction. A
degenerate case: five active constraints and two bounds on a flat family of designs
(the same restart in the same box moved to another design of nearly equal drive), where
the multipliers are not unique, so the Lagrangian's derivative is not the optimum's.
The sensitivities hold where the active set pins the design (the quadratic program);
on a drive objective under several gates, a restart in the wider box is the measure.

### A verdict on conflicting constraints (2026-11-18, item 36 continued) *measured*

On the variants of a staged re-pointing study kept outside the repository, constrained
runs whose constraints were incompatible ran to their time limit instead of
ending `:infeasible`: the augmented Lagrangian declares infeasibility only once its
penalty has grown 10⁸-fold without progress, and on that study its first inner run
alone did not finish within the limit. The designer then needs which constraints
conflict and by how much, in their own units; the study found those by hand, several
runs each. Three forms compared:

- (a) **phase 1** (`optimize_pointing(...; phase1 = true)`): before the main run,
  minimize the infeasibility ψ = ½ Σ max(0, gᵢ)² over every component's normalized
  violation (the augmented Lagrangian's limit as μ → ∞; zero exactly on the feasible
  set), ignoring the objective, by the same L-BFGS-B with the inner runs' value
  tolerance. It stops at the first design with ψ ≤ ½ feasibility² (every violation
  within `feasibility`), where the main run starts. Converged above it, or stalled (the
  least ψ fallen by less than 1% over the last 100 evaluations; below), with a
  violation over `feasibility`: the run ends `:infeasible` (`:infeasible_box` with
  parameters at a bound) at that design, no outer iteration;
- (b) **elastic**: the objective plus γ Σ h(gᵢ), h a hinge smoothed quadratically below
  10⁻³, in one unconstrained run; the hinge's arguments at the end are the slacks;
- (c) **a stall rule in the augmented Lagrangian**: infeasible once the penalty has
  grown in each of the last two outer iterations and ψ has fallen by less than 10%
  over them.

Cases: the closed-form problem of `test/constraints.jl` (48 quad powers; the drive
`mean_flux` is linear in the powers, so the largest drive at a total power and the least
power for a drive are fills of the quads in decreasing drive per power, a closed form) with
a 120 TW cap and a drive floor 5% above the cap's largest drive, 1 thread; the
dropped-bundle test case of `test/constraints.jl` (46 quads in `:phi` and `:z`, 0.45 mm wall
elements, drive objective, the four gates and the wall ceiling; 8 threads), made infeasible
by a box of ±0.02 mm (the gates against the box) and by the ceiling at 0.8 of the 48-quad
peak (the ceiling against the box; by itself the ceiling reaches 0.84 at best from this
start), budget 3,000 evaluations and 90 s. Gates tightened to 0.3 and 0.1 of their values,
and the ceiling to 0.9 and 0.86, all stayed feasible: on this case the gates and the ceiling
hardly conflict.

| case | augmented Lagrangian (today) | (c) stall rule | (a) phase 1 | (b) elastic, γ = 10 |
|---|---|---|---|---|
| closed form, cap against floor | `:infeasible_box`, 71 evaluations | — | verdict in 20 | — |
| closed form, cap 0.95 × start against floor 1.05 × its drive | 38 | 26 | 9 | 50 (γ = 1: 58, the floor given up by 0.49) |
| box ±0.02 mm | 48 (1.6 s) | 27 | 22 | 24, ends `:box` with the slacks (no verdict) |
| ceiling 0.8 | none in 90 s (2,195, 4 outer) | 1,837 (75 s) | 577 (about 25 s; 2,681 without the stall rule) | 1,979, iteration cap, slack 0.050 |

Feasible cases (no false verdict in any; drive as a fraction of the 48 quads'; phase 1
seeds the main run):

| case | augmented Lagrangian | (a) phase 1 (its evaluations) | (c) | (b) γ = 10 |
|---|---|---|---|---|
| closed form, start 2.5 × the powers (far outside the cap) | 37 | 56 (7) | 37 | 240–381 (γ = 1–100), 5e-6 from the KKT point |
| tutorial's start | 1,044, 0.962592 | 1,037 (11), 0.962284 | as today | — |
| random start in the box | 1,238, 0.962592 | 1,140 (19), 0.962522 | as today | — |
| ceiling 0.86 | 841, 0.958760, `:box` | 1,085 (323), 0.958871, `:interior` | as today | time limit (2,249), violation 9e-6 |
| gates × 0.1, 60 s | time limit (1,491), violation 2.7e-2 | time limit (1,477), 1.4e-3 | as today | 429, violation 0.021 |

Reading. The verdict of the augmented Lagrangian itself is cheap where the box stops its
inner runs at once (the box case) and slow where the inner runs creep (the ceiling case, and
the study, where the first inner run alone outlasted the limit); (c) only shortens the
outer iterations' share, so it cannot help the study, whose single inner run took the
time. The elastic form needs γ above the multipliers (γ = 1 gave up a floor that
γ = 10 kept), cannot tell "infeasible" from "γ too small", and on the feasible cases
ended with violations of 1e-5 to 2e-2 where the augmented Lagrangian meets them. Phase
1 reaches the verdict first on every infeasible case, from the infeasibility alone,
and costs a feasible run from 8% fewer to 29% more evaluations on the dropped-bundle
case and 51% more on the small closed form (the main run starts elsewhere, and on a
flat family of designs ends at another, within 3.2e-4 of the drive).

The stall rule. On the ceiling case phase 1 converges slowly toward its local minimum:
ψ within 3% of its final value after 400 evaluations, within 1e-5 after 2,000, the
L-BFGS-B convergence after 2,681 (with only the ceiling, within 1% after 140, done
after 766). Every feasible case reached the limits within 323 phase-1 evaluations, with
ψ falling by orders of magnitude per hundred. Ending phase 1 when its least ψ falls by
less than 1% over 100 evaluations (engineering defaults) gave the verdict at 577
evaluations with the violations within a few percent of the converged ones (the
ceiling 5.26% over against 5.10%). A feasible run that stalls so would need thousands
of windows to reach the limits, which the augmented Lagrangian does not do either.

**Relaxations** (`relaxation = true`, implies `phase1`). After a verdict, each
constraint violated or active at phase 1's design is relaxed alone: a constrained run
from that design minimizing a smooth largest positive violation of its components (a
64-norm, at most 2.2% above the largest for four equal components) with every other
constraint as a constraint (phase 1 on them first, which reports when they conflict
among themselves). The result is achievable (the design meets the others) and locally
least. Active constraints are included because the least-squares compromise of phase 1
can leave a constraint of the conflict exactly at its limit: on the closed form with
the cap at 0.95 of the start, the floor carried the whole violation and the cap sat at
its limit. Checks against the closed forms (120 TW cap, floor 1.05 × the cap's largest
drive): the floor relaxed alone 0.04761904762 (closed form 1 − 1/1.05 =
0.04761904762), the cap relaxed alone 0.05056360668 (closed form 0.05056360668, the
least power for the floor over 120 TW), in 133 and 82 evaluations. On the box case the
gates relaxed alone come out at 2.23 times their limits (304 evaluations), and the
ceiling relaxed alone is reported as not enough (the gates alone cannot be met in that
box: their phase 1 converged in 3 evaluations). On the ceiling case one relaxation ran
past the 90 s limit (status `:stopped`, the ceiling's 5.1% at its best design), so a
time limit bounds the relaxations too: they cost one constrained run each, which is why
they are an option.

Chosen: (a) with relaxations as an option, default off for both (a feasible run with
phase 1 starts its main run elsewhere and ends at another design of a flat family;
whether to switch it on by default is the user's decision; recommended, as it cost at
most 29% on the drive problems and turns a run to the time limit into a verdict;
**Superseded**: phase 1 on by default since "Defaults changed (2026-11-19)" §5). Today's run is
unchanged without the keywords: unconstrained runs and default constrained runs are
identical to those of the base commit. Confidence: two public infeasible cases and the
closed form; the stall rule's 1% per 100 evaluations is an engineering default
measured on one creeping case.

**On the study** (scratch runs, numbers not recorded here; one call each with a
time limit, `VisRadProfile(; taper = 1)`, 8 threads). On its case with two windows and
the entrance-hole keep-out as a `clearance_floor`, phase 1 gave the verdict in under a
minute, where the run without it had gone to its time limit; the gates relaxed alone
matched the largest gate multiple found by hand (runs minimizing it) within 2%, the
keep-out relaxed alone fell inside the bracket the hand-found feasibility runs had set
(that relaxation stopped by the time limit, the other constraints within a few percent),
and the ceiling relaxed alone was reported as not enough: the conflict the hand-found
runs had named. On its case with only the outer quads re-pointed (the keep-out as a
floor), the verdict came in seconds and the whole call, relaxations included, in under
two minutes: the gates relaxed alone matched the hand-found multiple within 1%, and
the ceiling and the keep-out relaxed alone were each reported as not enough. The last
was new: the hand-found runs had named the ceiling as the gates' conflict, while the
keep-out, held by a penalty there, conflicts with the gates as well, so that only the
gates' relaxation resolves the case.

**`mode_gate` from an evaluation.** `mode_gate(reference::Evaluation, gates)` takes
the shape from the reference's `flux_ylm` and records its scenario's `lpi` and time,
as `irradiance_ceiling` does (§ "The exact peak, restarts and the box"): an evaluation
with another model, or the same at another time, throws an `ArgumentError` (the check
is now shared by the ceiling, the gate and the drive floor). The reference may be on
another target (a hohlraum without the design's windows, as in the study), since the
gate compares shapes, not elements. The coefficient-vector method stays unchecked.

**Drive floor.** On the study, an absolute ceiling with a power scale per cone (≤ 1) was met
by cutting power, at a large loss of drive. `drive_floor(reference; fraction)` holds
`mean_flux` at or above `fraction` times the reference's (scaled by the limit; the model
check is as for the gate). Test: on the closed form, the least total power that keeps 0.9 of
the start's drive, the floor active with a positive multiplier, against the fill of the
quads in decreasing drive per power (1e-6). There is no default fraction: a reference with
more beams than the design delivers more drive, so 1 is not always reachable.

**`pointing_summary` and a clearance of −Inf.** The summary's least clearance per
keep-out filtered every non-finite value, so a beam the rule cannot place (−Inf, a
capture failure) vanished from it. It now leaves out +Inf only (a beam the rule does
not apply to); −Inf and NaN propagate and show (`clearance -Inf mm`). Excluding +Inf
never changes a minimum but for an empty set, where the result stays `Inf`.

## Constrained re-pointing: convergence on a feasible, active start (2026-11-18) *measured*

The augmented Lagrangian of `optimize_pointing(...; constraints)` (§ Constrained
re-pointing) was declared converged as soon as an inner run had converged at its
current tolerance with every constraint met and complementary. The first inner
tolerance is 1e-3 on the projected gradient and the first multipliers are zero, so a
start that already meets the limits (a restart from a previous result, or a quad
result handed to per-beam shifts by `parameters(p, beams)`) passed at once wherever the
objective's gradient is below 1e-3, as a drive objective's is per beam. Found on a
staged re-pointing study kept outside the repository, whose per-beam polish from a
converged quad result ended after a few evaluations. Fixed by declaring convergence
only once the inner tolerance has reached `g_tol`. Measured on the public cases (the
distance is relative, to the closed-form KKT solution):

| case | before: outcome, outer, evaluations, result | after |
|---|---|---|
| closed-form QP (`test/constraints.jl`), start violated | interior, 3, 35, 1.1e-4 from the KKT point | interior, 7, 54, 1.2e-8 |
| the same, start on the cap | interior, 3, 31, 6.9e-6 | interior, 7, 49, 3.5e-9 |
| objective × 1e-4, start violated | interior, 3, 12, 7.5e-2 | interior, 7, 83, 3.0e-4 |
| objective × 1e-4, start on the cap | interior, 1, 1, 7.0e-2 (the start) | interior, 7, 26, 7.0e-5 |
| dropped bundle of `test/constraints.jl` (0.45 mm wall) | interior, 2, 186; drive 0.96168, δ₁ at 0.41 of its gate | interior, 7, 1,044; drive 0.96259, every gate active |

On the dropped-bundle test case the old rule stopped short of the constrained optimum
(convergence reported with mode 1 far inside its gate and 0.09% of the drive left);
with the new rule every gate ends active, as on the tutorial's mesh before (§
Constrained re-pointing, 0.96259–0.96261). The cost is more outer iterations, about 1.5
to 6 times the evaluations here; where the old rule had already reached the optimum the
design does not change.

## Optimization defaults decided (2026-11-18) *measured*

The user's decisions of 2026-11-18 on roadmap items 34–36, each measured on the public
cases below.

**`irradiance_ceiling`'s `exact` defaults to `true`** (item 36). The peak over the
elements is what a ceiling holds by default; `exact = false` keeps only the smooth norm
(§ Constrained re-pointing, "Exact peak", for the method and its cost). The
dropped-bundle test case of `test/constraints.jl` already ran with `exact = true` and
does not move; it passes plain and under `--check-bounds=yes` (121 tests each), which
reach different designs of the flat family, so it asserts only what every admissible
design meets. The tutorial (8 threads; before → after the default):

| run | before (`exact = false`) | after |
|---|---|---|
| within limits (46 quads) | 997 evaluations, drive 0.9626, peak within 0.3% of the limit, either side | 1,079 evaluations (8 outer), drive 0.96258, peak 0.99999 × limit (norm 0.998), `:box` |
| stage 1 (quads, four limits) | 885 evaluations, drive 0.9577 | 1,092 evaluations, drive 0.95733, peak 0.992 × limit, gates l = 2 and 4 active, cap keep-out 17 µm (none at its 10 µm minimum) |
| stage 2 (per-beam polish, 120 evaluations) | drive 0.9584, every reported limit met (the norm) | drive 0.9593, peak 1.0082 × limit: reported as violated |
| stage 2 run to convergence (scratch) | — | 861 evaluations (9 outer), 326 s, drive 0.9592, every limit met, peak 0.997 × limit |

The budget-capped polish shows what the default changes: the exact rounds run only
after a constrained run has converged, so a run stopped by its budget reports the peak
as it stands (0.8% over here), where with `exact = false` the report shows the norm
(the earlier tutorial printed the peak beside a report of every limit met). The tutorial
keeps its 120-evaluation polish (a converged one costs five and a half more minutes of
the docs build), and its text states so; a design to hand off runs to convergence.

**`GroupPointing`'s `unit` defaults to `:quad`** (items 23 and 34). The user's
condition: "Quads-as-units is a default for pointing, but quad-splitting shouldn't
require jumping over hurdles." What pointing a quad's beams apart takes after the
change, as a user meets it:

- *Per-beam aims*: `by = :beam` (one keyword, as before); each beam's `:r` and `:phi`
  are then in its quad's port azimuth. `unit = :beam` (one keyword) puts them in the port
  azimuth of the beam.
- *A split about the quad's aim point*: there was no cylindrical split field
  (`hybride_mock_pointing` splits by a fixed table; `:split_theta`/`:split_phi` are on the
  capsule sphere). The new field is `:split` (mm, along the quad's φ̂, each beam `±split` by
  the side of its port's azimuth from the quad's mean port; positive opens the quad in the
  sense of its ports, as `:split_phi`; `hybride_mock_pointing`'s crossing `s` is
  `split = -s`), absolute (it sets each beam's tangential offset from its quad's mean
  template aim) or a shift with `relative`, needing `unit = :quad` and quad names. Two
  lines: `p = GroupPointing(nif_beams(); by = :quad, fields = (:z, :r, :phi, :split))` and
  `apply(p, parameters(p, "Q31T" => (z = 5.8, split = 0.5)))`. Checked against beams split
  by hand (each aim moved by ±split along the quad's φ̂ computed in the test from the port
  table) and against `hybride_mock_pointing` per cone, to 1e-12 mm; the template's split
  read back by `initial_parameters` and `parameters(p, beams)` to 1e-12; d aim/d split
  equals the side times the quad's φ̂ to 1e-14 (ForwardDiff); a quad held as one beam gets
  exactly none; inferred and JET-clean on Float64 and the optimizer's duals, absolute and
  relative.
- *A per-beam table* (`read_visrad_sst(path; template = nif_beams(quads))`, an as-shot
  table): `relative = true` keeps its beams' differences under one shift per quad, in
  the quad's frame now; absolute `:z`, `:r`, `:phi` collapse each quad to one aim point,
  and absolute `:split` keeps the table's mean split.

Every caller built on quads held as one beam (`nif_quads`, the tests, examples
and tutorials but one) or on beams without quad names (OMEGA) is identical either way;
the spherical fields do not use the frame. Results change only for beam
lists holding a quad's beams separately with an aim field in play:
`test/parameters.jl`'s per-beam case (pinned to `unit = :beam`, the case it checks),
`test/basics.jl`'s relative per-beam table (the new default intended; its checks hold
either way) and the dropped-bundle tutorial's per-beam polish (the new default
intended: its tangential shifts are along the quad's φ̂; measured with the next
change, below). `BriareusGUI` already
passed `unit = :quad` everywhere (left explicit). `Scenario(target, beams)`'s default
parameterization on a 192-beam list now has the quad frame; `examples/n221204.jl`
does not build a parameterization, and its coarse CSV (`0.8 0.4 full`) is byte-identical to
the base commit's (8 threads, 2 min 12 s each).

**The point–ellipse distance next to the major axis** (item 35). `_ellipse_distance`
(`clearance.jl`, behind `clearance` and `talis_check`) bisected on Eberly's Lagrange
multiplier s; inside the ellipse next to the major axis s + 1 ~ y₁ is held only to an
absolute 1e-16, so x₁ = y₁/(s + 1) lost its digits. Against the signed distance as the
minimum over the ellipse in BigFloat (256 bits; every local minimum of a 20,000-point
Float64 sampling refined by golden section), on 4,260 points (five ellipses, a/b from
0.3/0.5 to 1.0/0.2 and near-circular 0.40/0.39 and 0.62/0.55; along the major axis to
past the center of curvature and the vertex, 0 to 0.05 mm off it, both signs, and 300
random points inside and outside):

| method | largest error (mm) | points over 1e-12 |
|---|---|---|
| bisection on s (the old one) | 7.0e-5 (a point 3e-13 off the axis, inside) | 588 |
| `design_rules.jl`'s: Newton on s, points within 1e-6 e₁ of the axis snapped onto it, distance projected on the normal | 1.0e-8 (near-circular spot, next to the center of curvature) | 464 |
| Newton on σ = s + 1 (no snap), projected on the normal | 4.4e-16 | 0 |

The snap is not the fix the roadmap proposed it to be: next to a center of curvature the
closest point moves fast along the ellipse, and the snap's error there is first order
(1e-8 e₁), not 1e-12. Chosen: Newton on σ, which holds σ to full relative precision, in
`_ellipse_closest`, now shared (moved to `clearance.jl`): `_ellipse_distance` calls it
with `snap = 0` and projects on the normal; `HardwareKeepOut` and the TaLIS rules keep
the snap (`snap = 1e-6`, the default), which marks the medial-axis minima their
derivatives need. Tested (`test/basics.jl`) against the BigFloat minimum on 59 points
next to the axis, to 1e-13 (4.4e-16 measured). *Since 2026-11-18 (later)* the snap is
gone: the keep-outs mark medial minima by a mirror point (§ Capture and clearance rules
as optimization quantities, "The snap removed").

Report movement, new against old, on the test cases (`ellipse_reports`, scratch): the
mock-up of `test/design_rules.jl` (`clearance`, resolution 1e-7: off the deep
violations ≤ 4.4e-16 mm, 23 of 48 rows identical; below −0.75 mm up to +6.0e-5 mm;
nominal values, deep in the envelope, up to +6.0e-5 mm at −0.65 mm; `talis_check`'s
`:clearance_3w` worst values ≤ 4.4e-16 off the deep rows, +6.3e-6 mm at −0.99 mm);
`test/talis.jl`'s sombrero, LEH-lip and OMEGA-foil reports ≤ 4.4e-16 mm. The old
distance read low; every moved value moved up. The shared `_ellipse_closest` (Newton on
σ in place of s) moves the design-rule quantities by ≤ 2.2e-16 mm in value and ≤ 1.1e-9
in their derivatives (aim, defocus, spot axes; the old Newton on s lost the same digits
just outside the snap's 1e-6 e₁). The design-rule quantities against the reports on the
mock-up (resolution 1e-7 mm; keep-out minus report, before → after): `clearance` below
−0.75 mm, shield 9e-8 to 1.1e-6 → 9e-8 to 1.0e-7, stalk 6e-8 to 6.0e-5 → 6e-8 to
1.0e-7; `talis_check`'s `:clearance_3w` worst values on the rows whose worst point lies
on the medial segment, 1.7e-7 to 6.4e-6 → 1.0e-7 (the other rows 6e-8 to 1.1e-7, as
before). Every row now agrees within the reports' resolution, deep violations
included (`test/design_rules.jl` tests 2e-7 on every `:clearance_3w` row and 1e-6 on
every `clearance` row below 1 mm, where it had excluded those below −0.5 mm and
allowed 1e-5 on the medial rows); `clearance` still gives values up to 0.40 mm higher than the
keep-out on shield beams above 1 mm, where it skips points the cap hides.

The tutorial was rerun with this change and the quad frame (8 threads, 240 s). The runs on
quads without a hardware keep-out are the same (within limits, 1,079 evaluations). Stage 1,
whose cap `HardwareKeepOut` derivatives moved by ≤ 1e-9, takes another path through the flat
family of designs (1,092 → 1,526 evaluations, drive 0.95733 → 0.9574, the keep-out 18 µm,
gates l = 2 and 4 active, no aim at the box). The 120-evaluation polish (its tangential
shifts now along each quad's φ̂) ends at drive 0.9593, with the wall's peak 0.44% over the
ceiling, mode 2 0.12% over its gate and the keep-out active on six beams (the clearance
check gives 9.5 µm, with no beam failing). The tutorial asserts a normalized violation of at
most 0.01 for the capped polish.

**`keepout_penalty` and `clearance_floor` on a clearance of −Inf** (item 35). The penalty
skipped every non-finite clearance, so a −Inf (a `CaptureKeepIn` beam whose axis misses
every capture part when no part has a free edge to continue the margin from) was dropped
silently, and the constraint form turned it into a violation of +Inf, an infinite augmented
Lagrangian without a message. The options were (a) a large finite penalty (the objective
stays finite, but the cause is hidden and the gradient is zero); (b) an `ArgumentError`
naming the keep-out and the beam; (c) skip it (the old behavior). Chosen: (b), for −Inf and
NaN alike, in both `keepout_penalty` and `clearance_floor` (the check is shared); +Inf, a
beam the keep-out does not apply to, still adds nothing and holds as a constraint. A −Inf
arises only where the rule has no finite margin to offer, so no step of the optimizer can
recover from it; the run stops at the first evaluation with the beam named. Existing results
do not move (the sum runs over the same terms in the same order).

**The aperture keep-out's slope along the aim** (checked on a report from a
re-pointing study whose objective's derivative disagreed with central differences when
`ApertureKeepOut`'s penalty, weight 1e3, was active: −1.229 against −1.169, −1.193,
−1.207 at h = 1e-3, 1e-4, 1e-5 mm). The public mock-up is the dropped-bundle tutorial's
hohlraum and keep-out (1% contour, minimum 10 µm), the 48 quads each moved along a
random direction to within ±20 µm of the minimum. Against central differences along 240
random aim directions (5 per quad) the forward-mode slope agrees to second order in h
(largest difference 8e-3 at h = 1e-3, 3.4e-6 at 1e-5, 1.1e-7 at 1e-6 mm): no derivative
is dropped. Scans over ±20 µm (144 directions, 1e-5 mm steps, the largest slope change
then resolved to 1e-8 mm and both one-sided AD slopes taken) find (a) one genuine kink:
on one quad the slope jumps from −0.079 to −0.309 across a point (AD each side within
5e-9 mm), the rim minimum switching between two separate local minima of the distance
along the rim (a minimum of two smooth branches, concave); and (b) narrow regions of
high curvature where two such minima are about to merge (slope changes of 0.2–0.3 over
a few micrometers, ~300 mm⁻¹), which are smooth (the jump falls linearly with the step)
but put central differences at h ≥ 1e-5 mm far off. The squared hinge adds a jump in
curvature (2·weight·c′², 500 at weight 1e3 and slope 0.5) where a clearance
crosses the minimum, which central differences straddling it resolve only to O(h).
The reported mismatch is what a line search encounters across such points on several beams
at once; it is not an AD error. `clearance_floor` per beam keeps the kink in the
active beam's component (it does not remove it) but confines it to one multiplier;
a floor on a minimum of two branches is the intersection of two smooth floors, so the
constrained optimum is well posed, though the inner problems still see one nonsmooth
component when the beam sits on the switch. The penalty's fixed large weight makes the
kink and the curvature jump part of the objective itself. A smooth minimum over the rim's local minima (a soft-min
with a length scale below the keep-out's accuracy) would remove both (a design change,
not made: user decision). **Superseded**: made, "Defaults changed (2026-11-19)" §3.

## One pointing over several phases of a pulse (2026-11-18) *measured*

**What the workarounds were.** A multi-phase re-pointing study (kept outside the
repository) held the capsule's drive shape at several times of one pulse (picket,
trough, second step, rise, and the peak before and after the wall's windows burn
through) with one set of aims and a power scale per quad and phase. On the API of the
day it needed four constructions by hand:
1. a copy of every beam per scaled phase under the same name, so that a per-beam
   `GroupPointing` gives the copies one aim;
2. a `PowerScales` whose `by` counted its calls (it does not receive an index) to key a
   scale on (quad, phase);
3. a main `Scenario` without radiation that only carried the parameters, with the real
   per-phase scenarios (with separate wall albedos and windows, so separate deposition
   operators and radiation models) evaluated inside the objective; the coarse deposition
   of every copy in the main scenario was computed and discarded at each evaluation;
4. a `WeakKeyDict` keyed on the evaluation's internal flux slot behind a lock, so that
   the objective and every constraint shared one evaluation per phase while gradient
   chunks ran concurrently.

The box on the scales (each beam at or below the laser's limit in every phase), the
per-phase reports and the staging between parameterizations were by hand too.

**Design.** A `PhasedScenario` holds a `Scenario` per `Phase` (a separate target or a ready
scenario, with separate Scenario options) and the phases' beam powers (a number, one per
beam, a function of the beam, or pulses at the phase's time), all on one pointing
parameterization (any of `GroupPointing`, `FramePointing`, `HeldParameters`, relative or
not), whose parameters set the aims of every phase. Its parameterization (`PhasedPointing`)
appends a power scale per (group key, phase) for the phases not fixed; `by` keys the groups
as `PowerScales` does, so one key for several quads ties their scales (a pulse shape, with
one multiplier per phase). `evaluate` returns a `PerPhase` of evaluations (indexed by phase
name), computed once per parameter vector, which `optimize_pointing` passes to the objective
and to every constraint: the shared evaluation needs neither a cache nor a lock, and each
gradient chunk runs a separate evaluation. `in_phase(phase, c)` binds any constraint or
objective of one evaluation to a phase, its name prefixed by the phase in the reports
(`picket.mode_gate[l = 2]`), an exact ceiling keeping its exact value. `bounds(param;
power_max)` caps each scale at the least `power_max / P` over its group's beams in its phase
(`power_max` a number or a function of the beam) and refuses a fixed phase above the limit.
`remake(scenario; pointing, by, scaled)` stages (per quad, then per beam, then tied scales)
on the same operators, and `parameters(param, beams)` reads a result's per-phase beams into
the next parameterization (scales as the mean of each group's power ratios).

Options considered. (a) Only a wrapper parameterization (copies hidden inside an
`apply`) keeps `Scenario` unchanged but leaves the main scenario's discarded
deposition and the cache; (b) a lazily evaluated `PerPhase` (each phase deposited on
first use) saves a phase that no term reads, at the cost of a lock per phase and
mutable shared state under concurrent chunks; (c) the eager collection chosen here,
every phase evaluated once per parameter vector, a pure function of the parameters.
`PerPhase` is one container for beams, evaluations and summaries, so `r.beams[:picket]`
and `ev[:picket]` read alike. The phase's `time` is also the built scenario's time
(materials' albedos and the laser–plasma model at that time); a ready scenario keeps
the ones it holds. A phase can also be a `remake` of one shared scenario with the phase's
albedos and time (`history_state` gives them from a drive history), so one target,
view-factor matrix and deposition operator are used for every phase.

**Checks** (`test/phased.jl`, public mock-up: the 48 NIF quads in a 10.1 mm hohlraum
at h = 1.2 mm, a picket at albedo 0.6 (0.1/0.4 TW per inner/outer quad), a step from
pulses at 2.5 ns at 0.8, and a peak at 2 TW fixed at 0.85).
- One phase at the powers of the pointing, unscaled, is exactly the plain `Scenario` in
  the evaluation and in a constrained run of the same objective (identical `x`,
  value, evaluation count and constraint values).
- Against the duplicate-beam construction at the same parameters (aims per quad,
  scales per quad and phase): deposition and flux coefficients of every phase are
  identical (relative difference 0.0, here a sum of two equal ports for each quad's
  azimuth is exact; the test allows 1e-12), and the gradient of a mode departure equal
  to 1e-10 relative.
- The gradient of an objective over all three phases against central differences
  (h = 1e-6) at aims and both phases' scales: agreement to 1e-6 of the largest.
- A constrained run by cone (16 parameters: gates on every phase, drive floors on
  two, an LEH clearance) from a start that violates the gates: `:interior` in 474
  evaluations, every row met; the box at the laser's limit puts every scaled beam at
  2 TW at its upper corner and none above it.

**Cost per evaluation** against the duplicate-beam construction on the mock-up (192
parameters, one thread, the same per-phase scenarios): an evaluation 0.054 s against
0.067 s (median of 20; the main scenario's deposition of the copies, 0.011 s, is the
difference), a gradient 1.47 s against 1.82 s (−19%). On a scratch rebuild of the
study (about 550 parameters, six phases, its per-phase scenarios shared by both
constructions, at its saved per-beam stage): every phase's flux coefficients equal to
1.4e-16 and its deposition to 6.2e-16 relative (the copies' quad azimuths sum each port
four times, a different rounding), the objective equal to all printed digits, the gradient to
6e-14, and the box on the scales identical; an evaluation 0.055 s against 0.089 s and a
gradient 2.8 s against 6.2 s at 16 threads (−54%: the copies' coarse deposition in
dual numbers and the cache's lock are gone; another run may have shared the
machine). A constrained run on the phased API from that stage, the study's limits
written as `in_phase` constraints, ran under a 180 s limit: 51 evaluations, the
objective down by 0.12%, the limits' rows reported per phase.

Existing paths are unchanged (the optimizer's internal functions take the
abstract scenario type and an untyped evaluation; no arithmetic moved).

## Re-pointing tools: the inner budget, the smallest shift, a pulse-shape budget (2026-11-18) *measured*

From the multi-phase re-pointing study of 2026-11-18 (items 34 and 36).

**The inner budget** (item 36). In `optimize_pointing(...; constraints)` each inner
L-BFGS-B run is capped at `iterations` (default 200), and the multipliers change only
between inner runs. On the multi-phase study (276 parameters, six phases, about forty
constraint rows; scratch runs) a cap of 1,000 completed one outer iteration in an hour
and a cap of 100 six in 40 minutes. Since 320335b a run converges only after an inner
run at the final tolerance `g_tol` has converged, so a cap cannot simply stay small.
Candidates, for constrained runs only (measured with an internal switch, not kept):

- *fixed*: every inner run capped at `iterations` (the rule kept);
- *grow10*: the first inner run capped at 10 iterations, doubled after each outer
  iteration that ends with the infeasibility–complementarity measure V above
  `feasibility`, `iterations` once one ends within it (with the cap in the restart
  state, so that a run continued with `from` is the uninterrupted run);
- *grow10, μ held*: the same, with the penalty not raised after an outer iteration
  whose inner run stopped at a cap below `iterations` (ALGENCAN's tenfold rise assumes
  an inner problem solved to its tolerance);
- *np3*: a cap of max(20, n/3) for n parameters (the study's 100 at n = 276);
- *vs20*, *vs50*: an inner run also ends once V at its evaluations has fallen by less
  than 10% over the last 20 (50) evaluations while above `feasibility` (the inexact
  augmented-Lagrangian practice of stopping an inner run when the violation stalls; V
  recorded from the constraint values of the augmented Lagrangian, no extra evaluation);
- *ls100*: phase 1's stall rule on the inner value (a fall of less than 1% of
  max(1, |value|) over 100 evaluations).

Public cases (evaluations to convergence; drive relative to the 48-quad design; every
limit met within `feasibility`; all ended `:interior` or `:box`; the tutorial's stages
with the tutorial's 8-thread script, the others from a scratch harness on the same
problems):

| case (parameters, constraint rows) | fixed | grow10 | grow10, μ held | np3 | vs20 | vs50 | ls100 |
|---|---|---|---|---|---|---|---|
| dropped bundle of `test/constraints.jl` (92, 5) | 1,044; 0.962592 | 856; 0.962583 | — | 798; 0.961874 | 1,503; 0.962581 | 1,044 (identical) | 1,108; 0.962045 |
| the tutorial, within limits (92, 5) | 1,079; 0.962575 | 843; 0.962503 | — | 706; 0.961429 | 1,115; 0.962041 | — | — |
| the tutorial, stage 1 (64, 52: gates, ceiling, cap keep-out per quad, light) | 1,092; 0.95733 | 1,315; 0.9558 | 417; 0.95370 | — | — | — | — |
| the tutorial, stage 2 (per-beam polish, 120 evaluations) | peak 1.0082 × ceiling | peak 1.0153 × (the tutorial's 1% assertion fails) | peak 1.0091 × | — | — | — | — |
| mock-up: per-quad aims and power scales (138, 52: gates, exact ceiling, total power, an entrance-hole keep-out per quad) | 1,629; 0.960880 | 718; 0.961427 | — | — | 1,341; 0.962665 | — | — |

On the small cases of `test/constraints.jl` (the closed-form QP, the QP started on its
cap, the boxed and drive-floor cases, the one-parameter keep-out and light runs) every
inner run converges within 20 iterations, and grow10 is identical to fixed; the composed
run capped at 40 iterations ends with the same value in 218 evaluations instead of 204.
grow10 saved 18–56% of the evaluations on three cases, at drives within −0.007% to
+0.057% of fixed's (other designs of the flat family that meets the gates, §
Constrained re-pointing), but on the tutorial's stage 1 it took 20% more evaluations
and ended 0.16% lower in drive, and its budget-capped polish ends 1.5% over the
ceiling. Holding the penalty on capped inner runs ends stage 1 in 417 evaluations but
0.38% lower. np3 caps the final inner runs too and loses 0.07–0.12%; vs20 cuts the
first inner runs short but then creeps for hundreds of iterations once feasible (the
multipliers are then right and V does not stall); ls100 ends at the feasibility edge
with less drive.

The study's quad stage as posed in its script (276 parameters), three runs in parallel
at 8 threads under one 15-minute acquisition, each with the same remaining time after
its setup: about 550 s and 45–47 evaluations each (about 12 s per evaluation with three
runs sharing the machine). A cap of 1,000 takes the same path as a cap of 100 until the
hundredth inner iteration, which none reached, so the two are one row: fixed completed
no outer iteration (one under way), grow10 two (a third under way), vs20 one; the worst
normalized violation at the time limit was within 3% of the same value for all three,
every run far from its limits. The comparison the study needs (its hour-long run under
a rule) did not fit an acquisition.

*Decision:* the default stays (every inner run capped at `iterations`; constrained and
unconstrained results are unchanged). No candidate was better on every case:
grow10, the best on evaluations, cost drive on the stage with power scales and broke a
budget-capped polish. Documented instead (`optimize_pointing`): the cap sets how long a
run works with each set of multipliers; on problems with hundreds of parameters, slow
evaluations and a start far from the limits, a cap of about 100 works better (the
study's finding).
A rule that scales the first caps with progress and keeps the penalty schedule sound
on capped runs, measured on the study converged, is left open.

**The smallest re-pointing** (item 34). `shift_objective` is the weighted mean squared
aim shift from reference beams matched by name, `Σ w_b |Δaim_b|² / Σ w_b` (mm²), by
default from the scenario's start; with gates as constraints, `optimize_pointing` then
finds the least re-pointing that meets them (the per-entrance-hole study of a
multi-phase re-pointing study used `mean(sum(abs2, b.aim - ref.aim))`, the unit-weight
case). The weights are read once from the reference beams (a per-beam or per-cone
`weight` function), so the objective contains no callable on the hot path; it is generic
in the number type and inferred and JET-clean on Float64 and the optimizer's duals. Its
rms (the square root at unit weights) is `PointingSummary.shift_rms`. Checked
(`test/pointing_summary.jl`): the value against `Σ w x² / Σ w` for shifts `x` set by
hand, its ForwardDiff gradient against `2 w x / Σ w`, both to 1e-12; and the least
weighted re-pointing for one linear constraint, the mean axial aim of the 16 upper
outer quads raised by δ = 0.1 mm with the 44.5° quads three times dearer, whose
solution is `Δ_b = δ |S| (1/w_b) / Σ_S (1/w)` on the quads of S and zero elsewhere
(0.15 and 0.05 mm): the constrained run ends `:interior` within 1e-5 mm of it, the
objective within 1e-4 relative of the closed form.

**A pulse-shape budget** (item 34). Each distinct pulse shape needs facility
validation, so per-quad, per-phase power multipliers are grouped into a few shapes.
`pulse_shape_groups` does it as the multi-phase study did by hand: k-means on the
natural log of the multipliers (a key is a point with one coordinate per phase), so a
shape's multiplier per phase is the geometric mean of its members' and the fit is least
squares in the relative error; per class (the inner and the outer cones follow
different base pulses, so they never share a shape), with the number of shapes split
between the classes by the smallest total squared error over every admissible split.
The reported rms is `√(SSE / (n_keys · n_phases))` in log multiplier, the study's
measure. Method, chosen for determinism and exactness where the sizes allow:

- exhaustive enumeration of the partitions of a class into k groups when there are at
  most `exact_limit` = 10⁶ of them (the Stirling number S(n, k): S(10, 3) = 9,330,
  S(14, 3) = 788,970, S(15, 3) ≈ 2.4e6: a cone of up to 14 quads into three shapes, or 12 into
  four);
- with one phase, dynamic programming over the sorted values (an optimal
  one-dimensional grouping is contiguous; O(k n²)), exact at any size;
- otherwise the best of `restarts` = 100 k-means++ seedings (Arthur and Vassilvitskii,
  SODA 2007), each refined by Hartigan's single-key transfers (a key moves from group a
  to b when n_b/(n_b+1)|x−c_b|² < n_a/(n_a−1)|x−c_a|², until none does: a local optimum
  at least as good as Lloyd's). The seedings draw from SplitMix64 (Steele, Lea and
  Flood, OOPSLA 2014) started at `seed`, a stream fixed on every Julia version
  (`Random`'s default generators may change between versions). Exact k-means in several
  dimensions is NP-hard, and the study's sizes (about two dozen quads per cone, four
  phases, three shapes per cone) are beyond enumeration (S(24, 3) ≈ 4.7e10), so this is
  the path the study takes; `exact` in the result records which held.

Checked (`test/pulse_shapes.jl`) against an independent exhaustive search over all kⁿ
assignments (7 keys into 3 groups in 4 phases, 8 into 2 in 2, 9 into 3 in 3, 8 into 4 in
1): the enumeration and the restarts (forced by `exact_limit = 0`) reach the exhaustive
minimum to 1e-12; the split against the exhaustive search over the three splits of four
shapes between classes of 4 and 6 keys; the shapes' multipliers against the geometric
means in closed form. `shape_pulses` writes the shapes out: each class's base pulse
times the shape's multiplier, interpolated linearly between the phase times and held
outside them, tabulated at the base pulse's points and the phase times (the table is the
multiplier applied at those points; the exact product of two linear pieces is
quadratic). A peak run at the laser's limit has no multiplier: `peak = t` ramps the
multiplier linearly from the last phase's value to 1 between the last phase time and
`t` and leaves the pulse unscaled after it, an assumption of the tool (the phases
sample the pulse at their times only), not a fitted quantity. The result goes to
`write_pulse_table` (ten significant digits) and, through `assignment`, to
`PowerScales(p; by = b -> g.assignment[b.quad])` for the tied re-optimization.

On the multi-phase study's per-quad multipliers (46 quads, four phases; scratch runs,
2026-11-18), six shapes come out as the k-means of the study had them (200 Lloyd
restarts, written by hand): the same split between the cones, the same members and the
same multipliers to the digits it printed, for seeds 1 to 4 alike, in 4 s with
compilation.

- `picket_irradiance`'s `:mean` and `:average` read alike; renaming `:mean` to
  `:weighted` is a break of the public API (user decision, recommended).
