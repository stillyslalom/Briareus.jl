# Research notes: Target construction, frames, the GUI and compile latency

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

## Compile latency (2026-10-25) *measured*

Where users and the suite pay compilation, and what removed it. Julia 1.13.2, Windows,
4 threads, on a machine shared with other runs (wall times ±10%; compile times are
`@timed`'s `compile_time`). Base = 3950215, new = the same plus this change, each in a
separate worktree and package image. Scripts in the session scratchpad (`compile/ttfr.jl`,
`cand.jl`, `bench_history.jl`, `bench_opt.jl`, SnoopCompile 3 in a throwaway
environment, not a dependency).

**Where it went.** `@snoop_inference` over the test files: `_drive_history` was the
largest single cost, 52 s of inference in energetics/retune/bubble/dynamics/capsule_drive
and 4 s more in uncertainty/parallel_gradient. Every call site's `beams_at` closure (and
each `lpi` type, `record` closure and `times` type) is a new type, so each recompiled the
whole march: ~2 s per call site, and no precompiled copy could be used for it. Next:
`visrad_target(...; remesh)` (42 s inclusive in visrad.jl: `remesh_mesh` 10 s and
`fit_surfaces` 10 s of inference over ~3,600 concrete method instances; nothing
pathological, just never precompiled), `talis_check` 7.6 s, `write_visrad_kit` 7.4 s,
`dante` on a spectral solution 7.6 s, `spectral_solve` 3.7 s, and every new objective
closure re-inferring `optimize_pointing`'s body and Optim's L-BFGS-B (~1 s). Not
addressable in src: tests that differentiate through their own closures
(`ForwardDiff.gradient(fc, x)`) get a new dual tag each, so the deposition stack is
recompiled per test closure (uncertainty.jl: `_sum_points*` 15 s, `lpi_deposit` on
duals 13 s); `optimize_pointing`'s fixed `_PointingTag` already avoids this for users.

**Changes.** (1) `_drive_history` is `Base.@nospecializeinfer` with `@nospecialize` on
`model_at`, `operator_at`, `beams_at`, `times`, `lpi` and `record`: one inferred and
compiled instance for every closure. `@nospecialize` alone is not enough: the keyword
wrapper calls the body with concrete types, and inference specializes on them anyway.
The callables' values reach `solve`, `apply_lpi` and `lpi_deposit` by dynamic dispatch
once per step. (2) `optimize_pointing` likewise on `objective` and `callback`. (3) The
precompile workload adds a wall-loss drive history (with and without CBET), a spectral
solve and `dante`, `talis_check`, a CAD `cad_target` (cavity selection, remesh, cylinder
and sphere fits) with its view factors, and `write_visrad_kit`.

**Results (first call in a fresh process, wall s / compile s).**

| workflow | base | new |
|---|---|---|
| `drive_history`, first | 10.2 / 10.0 | 0.56 / 0.35 |
| `drive_history`, each new closure | 2.1 / 1.9 | 0.17 / 0.01 |
| `optimize_pointing` (3 it.), new objective | 1.65 / 1.07 | 0.45 / 0.10 |
| same, objective with keep-outs | 3.95 / 1.38 | 2.15 / 0.13 |
| `spectral_solve` | 4.6 / 4.6 | 0.04 / 0.04 |
| `dante(dante1(), sp)` | 8.9 / 8.4 | 0.21 / 0.03 |
| `cad_target(...; remesh)` (small can + capsule) | 34.2 / 34.3 | 0.06 / 0.02 |
| `talis_check` | 6.1 / 8.5 | 0.04 / 0.00 |
| `write_visrad_kit` | 13.5 / 11.4 | 1.7 / 0.03 |
| `DepositionOperator`, `view_factors`, `deposit` | ~0 compile | ~0 compile |

Test files visrad, talis, energetics, retune, dynamics, capsule_drive, spectral,
parallel_gradient (4 threads, run side by side): compile 197 → 113 s, wall 247 → 168 s.
The package image builds in 162 → 192 s on the loaded machine (95 → 143 s quiet):
the cost moves from every first call to every rebuild, which `BRIAREUS_NO_WORKLOAD=true`
still skips.

**Hot loops unchanged.** Drive history, 1,300 elements, 24 steps, wall loss and a
recorder, warmed: min of 5 runs, four alternating pairs, base 1.80/1.77/1.97/1.89 s, new
1.98/1.79/1.82/1.84 s (noise), allocations 196 → 182 MiB. Results are unchanged (every
field of the history, `isequal`). `optimize_pointing` (default NIF hohlraum, keep-outs,
8 iterations): 0.379 vs 0.389 s, `x`, value and evaluation count are identical.

**Not done.** The GUI's 0.75 s first-import tick and 0.16 s first field-coloring frame
(ROADMAP "GUI on CAD") were not measured (they need the Qt window); the GUI's time panel and drives
call `drive_history`/`optimize_pointing` with their own closures and now reuse the
package's precompiled copies. `view_factors` on a CAD target still compiles ~2 s at
first call when its parts differ from the workload's (sampled visibility on other
element shapes). A robust objective (`uncertainty = ...`) is not in the workload.

### Update (2026-10-25): the workload cut to the common session

The user decided that the workload should cover what most sessions hit, and the VisRad kit
is not that. The common session is a GUI-less script: build a `hohlraum`, NIF quads, a
`Scenario` (deposition operator, view factors, radiation model), evaluate, moments,
`optimize_pointing` with keep-outs, solve, report, clearance, and a `drive_history`. All
nine tutorials and eleven examples build a `hohlraum`; four tutorials optimize pointing;
`drive_history` is in six examples and four tutorials. Dropped from the workload:
`write_visrad_kit` (GUI export only), `talis_check` (GUI only), `cad_target(...; remesh)`
and its view factors (one example, one tutorial), `spectral_solve`/`dante` (two/three
examples, the GUI's Dante panel) and the CBET drive history (three calibration examples).
The `@nospecializeinfer` changes stay: they are what make each new closure cheap, workload
or not. The pre-existing large-model, H-matrix and ray-trace items were not revisited.

Package image (touch a src file, `julia --project=. -e 'using Briareus'`, same machine
and session, another agent's runs alongside): 138 s → 96 s (twice 96 s; a first,
colder build of the old workload took 338 s and is discarded). First call in a fresh
process, wall s / compile s (scratchpad `ttfr.jl`, small `hohlraum(nz = 10, nphi = 16)`
for all but the optimize row):

| workflow | old workload | new |
|---|---|---|
| import | 1.8 | 1.5 |
| kept: `Scenario` (default-resolution NIF hohlraum) | 4.58 / 0.02 | 4.57 / 0.02 |
| kept: `optimize_pointing` with keep-outs, 3 it. | 0.17 / 0.09 | 0.17 / 0.09 |
| kept: `clearance`, `solve`, `report` | 0.35 / 0 | 0.34 / 0 |
| kept: `drive_history` (wall loss) | 0.43 / 0.31 | 0.44 / 0.33 |
| dropped: `drive_history(...; lpi = CBET(...))` | 0.13 / 0.07 | 1.9 / 2.2 |
| dropped: `spectral_solve` | 0.05 / 0.04 | 3.2 / 3.2 |
| dropped: `dante(dante1(), sp)` | 0.17 / 0.01 | 8.1 / 7.9 |
| dropped: `talis_check` | 0.03 / 0 | 4.2 / 5.8 |
| dropped: `cad_target(...; remesh)` | 0.06 / 0.02 | 24.7 / 24.7 |
| dropped: `view_factors` on that CAD target | 0.51 / 1.1 | 0.56 / 1.3 |
| dropped: `write_visrad_kit` | 1.3 / 0.09 | 7.7 / 6.5 |

The kept paths are unchanged; each dropped path pays its compile once per session that
uses it (CAD import the most, 25 s), against 42 s off every rebuild for everyone.

## GUI on CAD: import, parts and first frames, remeasured (2026-10-26) *measured*

Measured locally (8 threads pinned to one CCX; the NIF CAD workspace of the 2026-10-18
notes: 15 CAD objects and a sphere, the 133k-facet mount and two far parts drawn only,
12 parts in the calculations), coarse mesh (h = 0.6 mm).

**`remesh_mesh` on the CAD capsule (H).** Since the sliver merge of 2026-10-20
(`_merge_small!`/`_dissolve!`), remeshing the capsule part (38,200 facets, sphere and
fill tube) did not finish in 3.5 min at h = 0.6 or 0.28 mm (it took ~1 s before that
commit), so the GUI's import of the workspace never completed. The capsule's small
elements hold up to 4,000 facets; dissolving one tried every remaining facet against
every neighbor element after each move, each try a pass over the neighbor (≈2,400
facets): 1.9 M passes in 160 s, almost all rejected because the new facet itself lay
beyond `rmax` of the new centroid, and 267 of 298 dissolutions ended stuck and were
undone. Now a facet is retried only after a neighbor of it moved or an element it tried
grew, the element sums are extended by the new facet (the same left folds), and the
only the flatness and reach of the new facet are tested before the pass over the element (the
answers are those of the full test). Capsule: 4.1 s at 0.6 mm (663 elements), 8.4 s at
0.28 mm (978); before the merge commit ~1 s, so the merge pass still costs 3–7 s on it.
Results are unchanged: the remeshes of the 13 other CAD parts and two synthetic meshes at 0.28
and 0.6 mm equal the old code's (`==` on vertices, triangles and patches);
test/meshops.jl and test/cadfit.jl pass. GUI import of the workspace (coarse, no window,
first call): 12 s.

**Large small elements are no longer dissolved (H; decision 16).** Most of what was
left was still `_dissolve!`: 10.8 of 13.9 s at 0.28 mm, 5.2 of 6.7 s at 0.6 (an
instrumented copy of `_merge_small!` logging each dissolution's triangle count n,
outcome and time). The capsule's are dominated by a few of n = 10³–4.6 × 10³ at ~n²
element tests each, mostly stuck (0.28 mm: 446 of 520 stuck, 233 of them over 100
triangles; 0.6 mm: 656 of 803). Elements of more than `_DISSOLVE_MAX` = 256 refined
triangles are now left whole. The threshold, from the largest dissolutions that succeed
elsewhere (every remesh with dissolutions logged, no skip): the test meshes of
test/meshops.jl and test/cadfit.jl 112 (two of 233 over 100), the remeshes of
test/hmatrix.jl and test/partial_visibility.jl 49, the cad_hohlraum tutorial (0.5/0.15
and 0.25/0.1 mm) 24, the 13 other CAD parts of the workspace 60, `sphere_mesh(0.9;
ntheta = 60, nphi = 120)` at 0.25–0.6 mm 66, and a 40k-facet sphere
(`sphere_mesh(1.1; ntheta = 100, nphi = 200)`) at 0.6 mm 237; the sliver tubes (test
meshops' 360 strips, the GUI's 2,000) and `cylinder_mesh(...; nphi = 400)` have none at
0.25–0.6 mm. A dissolution that ends stuck is undone, so the result changes only where
a dissolution over 256 would have succeeded: all of the above are `==` the old code's
(46 meshes compared directly, the test and tutorial meshes by that criterion), and the
test files pass. Capsule time against the threshold (one run each, s, 0.28 / 0.6 mm):
none 13.9 / 6.7, 512 2.28 / 1.21, 256 2.07 / 1.01, 128 2.06 / 0.97, 100 2.02 / 0.98; at
256 the dissolutions take 0.1 / 0.05 s and clustering is the rest. 100 would change two
test meshes (n 108, 112), 128 the 40k-facet sphere (237) and more of the capsule's
successful dissolutions (n 131–142 and 192–255); 512 costs 0.2 s more on the capsule.
256 keeps every public case at no measurable cost over 100. Timed
three times after a warm-up (8 threads, the machine shared): before 12.4–12.7 s at
0.28 mm and 4.1–4.2 s at 0.6 (the 8.4 s above was a less loaded run); after 2.0–3.3 s
and 1.0–1.5 s.

The capsule mesh changes where a dissolution over 256 succeeded in the old run (0.28 mm:
14, n 318–2,094; 0.6 mm: 8, n 260–1,767): 36 elements (0.17% of the capsule's area) become
48 at 0.28 mm, 978 → 990 elements; 21 (0.18%) become 28 at 0.6, 663 → 670; all in one
thin band of the sphere. Elements under 0.3 size²: 445 → 457 and 620 → 627 (at 0.6 mm
the angle, not the size, bounds the capsule's elements, so almost all are "small": the
merge cannot help there). Unchanged: max tilt 9.995° / 9.997°, max reach 0.9998 /
0.9996 size, every centroid on its element, max diam/√A 22.3 / 29.9 → 29.8, the area
to 3e-14, closed and conforming (no open or non-manifold edges). Output triangles
+7% at 0.28 mm (59,958 → 64,280), +1.6% at 0.6. The skip also applies to finer CAD
spheres: `sphere_mesh(1.1; ntheta = 140, nphi = 280)` at 0.6 mm goes 1.6 → 1.2 s and
238 → 256 elements, at 200 × 400 5.3 → 2.2 s and 244 → 274 (a facet-count threshold
skips more on finer faceting). Scripts: scratchpad `dsk/` (`instr.jl`,
`survey_cad.jl`, `survey_tests.jl`, `survey_more.jl`, `time_capsule.jl`,
`cap_stats.jl`, `mock.jl`).

**A curvature-scaled minimum area: measured, not adopted (H).** The merge's minimum area,
`minarea` size² = size²/3, counts almost every element of a curved part as small where the
angle, not the size, bounds the elements. One variant tried was to scale it per element by
min(1, R sin(angle)/size)², R sin(angle) being the reach of a spherical cap within `angle`
of its axis (as `rmax` = size is on flat parts), with 1/R the area-weighted curvature of the
element's facets (each facet's largest dihedral angle to an edge neighbor within `angle`,
over the distance between their centroids; 0 on flat parts, so those stay unchanged).
Replayed on every remesh of test/meshops.jl, test/cadfit.jl, test/partial_visibility.jl,
test/hmatrix.jl (54 meshes) and the cad_hohlraum tutorial (15), synthetic spheres, tubes and
cylinders, and the NIF CAD capsule (timing and counts only). 48 of 85 unchanged (flat parts,
tubes, the tutorial capsule at 0.15/0.1 mm, the finest spheres); 37 changed, all by **more**
elements: the merges it stops are unions of angle-bounded elements that pass every
acceptance test, not slivers. `sphere_mesh(1.1; ntheta = 100, nphi = 200)` at 0.6 mm 229 →
270 elements, `sphere_mesh(0.9; ntheta = 60, nphi = 120)` at 0.4 and 0.6 mm 221 / 218 → 264,
the test meshes up to +12% (cadfit's sphere at 0.3 mm 474 → 531), the tutorial's wall parts
at 0.5 mm 457 → 462, 465 → 468, 21 → 24, and the CAD capsule 670 → 865 at 0.6 mm, 990 →
1,016 at 0.28. Elements a little more compact (that sphere at 0.6: diam/√A median 1.42 →
1.36, 95th percentile 1.53 → 1.46; max tilt 9.999 → 9.317°). The capsule's slivers (min area
3e-6 size², diam/√A to 29.8) are untouched: they are under either minimum, and no neighbor
takes them in either. Downstream, the inside of the spheres as an integrating sphere (F_ij =
A_j/A exactly): flat elements, closure max 6.5e-3 → 5.0e-3, 1.1e-2 → 4.9e-3, 5.8e-3 → 4.7e-3
(the three spheres at 0.6 mm), mean |closure| 4.3e-3 → 3.6e-3, but the rms relative error of
F_ij worse (6.8e-3 → 7.5e-3, 2.3e-3 → 1.4e-2, 4.5e-3 → 6.3e-3); fitted spheres exact to
1e-14 either way. The time of the capsule's remesh is unchanged (1.0 s at 0.6 mm, 2.0 s at
0.28), and `_DISSOLVE_MAX` is still needed (without it 4.1 s and 8.3 s: the slow
dissolutions are of elements small under either criterion). Net: +18–29% elements on curved
parts at sizes the angle bounds (view-factor cost ~1.4–1.65×) for slightly more compact
elements, no speed-up and no sliver removed; not adopted. Scripts: scratchpad `amin/`
(`record.jl` logs every remesh input of the test files and the tutorial, `replay.jl`
remeshes and compares, `vf.jl` the integrating sphere, `meshops_scaled_amin.jl` the
variant).

**The window's GUI-thread stalls (H).** Scripted Qt window (scratchpad
`gui_timing.jl`: a coarse hohlraum's drive, coloring by a field, a selection, the
workspace import, the largest part left out and put back, the template and the import
again), the longest gap between ticks after each step, first and warm, two runs each;
"before" is the GUI as of 576587b with the remesh fix. Most of the first-call stalls
were compilation on the GUI thread, not GUI work. The compilation was of two kinds:
(1) `_run_job(job, spawn)` was specialized on each job's closure, and inferring its
inline branch inferred the whole job (the import's reading and meshing, a part's
remeshing, the drive's scenario) on the GUI thread before spawning it; it is now
`@nospecializeinfer` on `job`, so the job compiles on its worker. (2) A worker's result
is applied on the GUI thread by a closure compiled at its first call there; the worker
now compiles it first (`_deliver_compiled`, a `precompile` of the commit for the
result's type; a lookup once compiled). Compilation the GUI thread still needs waits
behind the workers' (the code generator's lock), so a few ms of it can stall for
hundreds while a worker compiles the CAD remesh.

| step (window) | before | after | warm |
|---|---|---|---|
| drive requested (call) | 0.39–0.40 s | 0.01 s | — |
| first coloring by a field | 0.36 s | 0.21–0.24 s | 0.03 s |
| first selection | 0.22–0.23 s | 0.20 s | 0.03 s |
| import requested (call) | 0.72–0.73 s | 0 | 0 |
| import's commit tick | 2.2–2.4 s (commit 0.95–1.0 s) | 1.5–1.6 s (commit 0.27–0.28 s) | 0.14–0.19 s |
| a part left out (call; tick) | 0.35–0.38 s; 0.36–0.38 s | 0.09 s; 0.17–0.21 s | 0.03–0.05 s; 0.07–0.16 s |
| a part put back (call; tick) | 0.86–0.87 s | 0; 0.13–0.17 s | 0; 0.12–0.15 s |

Without a window (scratchpad `gui_import_ticks.jl`, a session with a queue, first
call): import call 1.33 s → 8 ms, its commit 0.85 → 0.05 s; leave out 0.33 →
0.05–0.09 s, put back 1.16 s → 0. Also, `import_visrad_workspace!` is
`@nospecializeinfer` on `include` and `done` (each new `done` closure recompiled its
body, 0.85 s); the GUI workload leaves its disk part out and puts it back through a
window-style queue and starts an import of a missing file (the call's path), and
precompiles the import's `set_beams!` call (0.17 s of its commit). The field colorbar
warm-up at the window's first frame now sets the field's colormap, scale and label and
a second scale at the next frame (the colormap and gradient texture updates were
compiled at the first coloring). Left: the import's commit tick is Makie's first
drawing of the CAD scene (the ghost and highlight meshes, the beams' line segments,
camera fitting: compute-graph edges of 50–80 ms each) and the shell's first rows of a
VisRad part, ~1.3 s; the first field coloring and the first selection, ~0.2 s each, are
GL-side first drawings too. Drawing the highlight from a degenerate triangle at the
first frame, as the ghost is, did not remove the first selection's stall (0.57 s in
one run) and was dropped.

## Frames and placement (2026-10-28)

Frames, frame-relative pointing and placement sweeps (`src/frames.jl`,
`src/frame_pointing.jl`, `src/sweep.jl`; design in `reference/frames-design.md`). The
conventions below are assumed, not taken from a facility document.
- **Port frame x̂** **L**. `port_frame(θ, φ; roll)` puts `z` toward the port and `x`
  along the port's θ̂ (turned by `roll` toward φ̂). That equals
  `euler_frame(origin, θ, φ, roll)`, VisRad's angle order (Polar, Azimuthal, Rotation)
  read as for its objects. VisRad's port (pointing-method 1) frames have not been seen
  in a sample workspace; check the x̂ there before relying on a port frame's roll. The
  origin at TCC and `z` toward the port are **H**: VisRad's coordinate transform tool
  documents its diagnostic-port-axis frame so
  (prism-cs.com/Manuals/VisRad/additional_topics/coord_transform_tool.html, read
  2026-11-02); the page gives no `x`. `port_frame(facility, port)` takes NIF quad ports
  from `NIF_PORTS` (NIF User Guide Table 4-1) and OMEGA ports from `omega_ports`.
- **Euler order and sense** **L**. `euler_rotation(polar, azimuth, rotation) =
  Rz(azimuth) Ry(polar) Rz(rotation)`, right-handed, the reading already assumed for
  VisRad object orientation ("VisRad workspace files (2026-10-17)"); the VisRad
  workspace reader now calls it and gets the same matrices as before.
- **Cylindrical coordinates** *decision*. `coords_of`/`point_of` order cylindrical
  coordinates (r, φ, z), as the GUI's part positions do; VisRad's summary table is
  read as (r, z, φ) (`_vrw_point`, a separate assumption). The two are kept apart.
- **`component_frame` rules** *decision*. Origin at the area-weighted centroid; `z`
  from the exact surface's axis, else the area-weighted mean normal when its length is
  at least 0.1 of the area (an engineering threshold: a flat foil 1, a hemispherical
  cup 0.5, a closed body 0), else the principal axis of the area moments that differs
  most from the other two, signed toward chamber `+z`. A sphere or cube falls back to
  chamber `z`. Ambiguous by nature; scripts that care give the frame.
- **Sweep defaults** *decision*. A sweep's re-optimization box is `bounds` about the
  initial parameters of each step (not about the warm start), so the box does not drift
  along the sweep. The backlighter tutorial's clearance check drops the NIF clearance
  wings (`nif_wings`, NIF engineering defaults) for its OMEGA-style mock-up.
- **Placement steps** *decision* (`src/placement.jl`). Steps apply in order, each in
  closed form, not through a simultaneous solver. A turn is the smallest rotation taking
  the frame's `z` to the target direction; `z` exactly opposite turns 180° about the
  frame `x` itself (a branch chosen for determinism). A relation conflicts with an
  earlier one when it raises that relation's residual by more than 1e-9 mm or 1e-7°
  (engineering tolerances, far above Float64 round-off on mm-scale frames and far below
  any placement that matters); relative moves (`Offset`, `Rotate`, `Roll` without a
  reference) never conflict and leave residuals instead. Degrees of freedom are counted
  as fixed translation directions and rotation axes at the orientation after each step,
  a linearization: a later tilt of `z` can change a fixed roll, which the residual check
  then catches as a conflict.
- **Placement steps in the GUI** *decision*. A new step starts where it leaves the part in
  place when the step allows: a coincident point at the origin's current place, the
  current distance, the current roll from the reference frame's `x` (its `y` when `x`
  lies along the part's `z`), flips and the facing side chosen to the nearer turn; an
  on-axis step projects (no distance); relative moves start at zero. A step sweep
  re-resolves the frames at each value and moves each part whose frame changes by
  `F(v) ∘ F(v₀)⁻¹`; the tests check that against the design rebuilt with the value, to
  1e-12 mm on vertices and aims. Session version 4 stores the steps.

## Halfraum ends and starters (2026-11-02) *derived*

`hohlraum(; top_end, bottom_end)` and the GUI's halfraum starters (assemblies phase C).
- **Open end.** An end given `:open` has no cap and no entrance hole; the target is not
  declared `closed`, so its view factors are not balanced and the power that leaves
  through the end is `energy_balance`'s `lost`. Checked against the coaxial-disk closed
  form (equal disks of radius R a distance L apart, `F = (X − √(X² − 4))/2`,
  `X = 2 + (L/R)²`; 0.17157 at R 1.5, L 3 mm): a halfraum closed at the top and open at
  the bottom, h 0.3, 0.2 and 0.12 mm, gives the cap's escaping fraction within 0.75–0.95e-3
  of it and the wall's escaping area `A_d (1 − F)` within 0.7–1.3e-3; the power lost
  through the open end with the wall driven matches the power escaping through an
  aperture disk closing the end (with view factors for the disk) within 0.8–1.5e-3. The
  errors are the view factors' (they do not fall steadily with h). `test/halfraum.jl`.
- **Per-end radius.** An end given a radius other than `leh_radius` takes separate ring
  counts from `h`; one given `leh_radius` is exactly the default end.
- **Closed end.** The cap and the hole's mesh (and its lip's) merged, of the wall's
  material: the GUI hohlraum part's closed end since phase B, moved into core unchanged.
- **Starters: illustrative mock-ups, no experiment's dimensions.** *Halfraum pair with
  shock tube*: two halfraums 2.5 mm long, radius 1.5 mm, entrance holes of 1 mm radius at
  their outer ends, facing each other across a 2 mm package (a cylinder of the halfraums'
  radius, absorbing) whose end disks close the open ends; the top disk is observed.
  *Opacity halfraum*: a halfraum 3 mm long, radius 1.5 mm, an entrance hole of 1 mm
  radius at the top and the bottom open, a 1 mm square sample (albedo 0.3, no laser
  conversion) in the open end's plane facing in. The scale (millimeter halfraums, mm-size
  packages and samples) is that of the public descriptions of NIF's halfraum-driven
  hydrodynamics and opacity platforms
  (https://lmf.llnl.gov/about/keys-to-success/targets); the numbers are chosen, not taken
  from any shot. Package and sample albedos are the starter capsule's 0.3 (an engineering
  default).
- **Balanced per cavity (2026-11-03).** Until then the GUI balanced view factors only for
  one hohlraum part with no open end, and left the shock-tube pair (two cavities, each
  closed by a package face) unbalanced, its rows missing by the view-factor error
  (~1e-3). Each halfraum and its package face are now balanced on their own; see
  "Per-cavity balancing" below. The opacity halfraum stays unbalanced: 17% of its rays
  leave through the open end around the sample.
- **The seam (2026-11-02).** The package disks have the halfraums' radius and lie in their
  open ends' planes exactly (z ±1.0 mm); no overlap. The GUI meshed disks, washers,
  cylinders and cones with azimuthal counts in multiples of 4 and core `hohlraum` its wall
  in multiples of 8 (`hohlraum_counts`), so at the standard element size (0.28 mm) the
  disk's rim was a 36-gon inside the wall's 40-gon end ring: a gap of 6.8e-3 mm², 9.6e-4
  of the end's area (closed form: the difference of the inscribed polygons' areas,
  `n/2 R² sin(2π/n)`); at 0.6 and 0.17 mm the counts agreed (16, 56). Measured on the
  upper halfraum closed by its disk (the enclosure alone, unbalanced view factors, the
  wall driven; the energy balance exact to 5e-15): the gap raised the view-factor area
  leaving the enclosure, Σ A(1 − ΣF), by 3.8e-4 mm² and the power lost to open space by
  2.3e-5 of the source, against a closure error of the view factors themselves of
  6.4e-4 of the source. The GUI's circular parts now take core's rule (`_nphi`: multiples
  of 8), so a part of a hohlraum part's radius and element size meets the wall ring
  at every vertex; at 0.28 mm the closed configuration then equals the one with the
  disk remeshed at the wall's count exactly. Circular parts meshed at sizes where the two
  rules differ change their meshes (results change at the view-factor error).

## Hohlraum profiles: frustraums and multi-section walls (2026-11-02) *derived*

`hohlraum(; profile)` (`HohlraumProfile`, `frustraum_profile`) and the GUI's frustraum
starter (assemblies phase D). Sources read 2026-11-02; the paper's numbers from its text,
checked against the page image of Fig. 1.
- **Frustraum geometry (public).** Baker et al., Phys. Plasmas 30, 092708 (2023),
  doi:10.1063/5.0163396 (OSTI 2222631): two conical frustum-shaped hohlraums joined at the
  capsule's equator (Fig. 1); in all but one shot 9.76 mm in diameter at the equator,
  12.24 mm long, 23° wall angles, entrance holes 3.1 mm in diameter at both ends, wall
  area quoted as 270 mm²; across the campaign wall angles of 21°, 23° and 25° and
  equatorial diameters of 7, 9.26 and 9.76 mm; HDC capsules of 1.2 mm inner radius and
  ~81 µm shell (1.4 mm in some shots). Confidence: high (stated in the text). Design
  paper: Amendt et al., Phys. Plasmas 26, 082707 (2019), doi:10.1063/1.5099934 (not read;
  paywalled at OSTI 1557379). LLNL news, 2019-11-13
  (https://lmf.llnl.gov/news/frustraum-hohlraum-design-shaping): the name dates from
  2018; full-scale 9.2 mm wide, subscale 7 mm.
- **Derived end radius.** The paper does not say what its wall angle is measured from.
  Taken from the axis, the end radius is 4.88 − 6.12 tan 23° = 2.282 mm (*derived*,
  medium confidence). The cones' area is then 299.2 mm² and the caps' 17.6 mm² (holes of
  1.55 mm radius), 317 mm² against the quoted 270 mm²: the quoted area is **not
  reproduced** from the quoted dimensions (a wall angle from the end plane is too steep
  to fit at all). The drawing's ends are not resolved well enough to explain it. The
  starter uses the stated dimensions as an illustrative mock-up of the shape; its wall
  area is not the paper's.
- **No published comparison used as a check.** Baker et al.'s APS DPP 2023 abstract
  (LLNL-PROC-851230, OSTI 2204932) gives the frustraum's wall area as ~25% smaller than
  a 6.4 × 11.3 mm cylinder's and its efficiency ~14% higher, while the 2023 paper gives
  270 mm² against 257 mm² for the 6.4 × 11.24 mm HYBRID-E cylinder (slightly larger).
  The two disagree and neither specifies the geometry fully, so neither is a check here.
- **Four-section (tentless) hohlraum: not public beyond a caption.** The LMF targets page
  (https://lmf.llnl.gov/about/keys-to-success/targets) captions a photo of a prototype
  tentless target whose hohlraum is made in four sections, the top and bottom ones
  sandwiching 50 nm polyimide membranes (with an amorphous-carbon layer) that hold the
  capsule. No section lengths or radii are given. Science & Technology Review, March 2018
  (https://str.llnl.gov/2018-03/stadermann) describes a hohlraum made in four segments
  whose walls clamp the tents' edges, also without dimensions. The four sections are a
  fabrication split; whether their radii differ is not public. A cylindrical hohlraum part
  switched to *Sections* in the GUI therefore starts as the same cylinder in four equal
  sections (a radius at each of five breaks, all the cylinder's), with no radii invented;
  a radius per break is then the user's.
- **Mesh.** Each section between breaks is a component (`wall_1`, … from the bottom,
  `part = :wall`) with its exact surface: `CylinderSurface`, `ConeSurface` or, for a step
  in radius, `DiskSurface`. A cone section's patches are trapezoids between parallel
  chords at the same azimuths (planar); all sections share the azimuthal grid, so they
  share every vertex. A profile of one cylindrical section is the cylindrical
  hohlraum exactly (`test/profiles.jl`).
- **Validation, closed forms.** A closed frustraum (length 4, waist radius 2, end radius
  1.2 mm, both ends closed, no capsule): every exchange follows from the coaxial-disk
  factor alone, because the ends see each other unobstructed and the waist plane splits
  the cavity (end to end F(re → re, L); end to the far cone F(re → rw, L/2) − F_ee; end to
  the near cone 1 − F(re → rw, L/2); cone to cone through the waist disk by reciprocity;
  cone to itself the remainder). Exact values 0.07672, 0.37864, 0.54464, 0.40607 and
  0.40105. With the cone surfaces, h 0.4 / 0.25 / 0.16 mm (576 / 1512 / 3200 elements):
  end to end −1.8e-4 / −7.3e-5 / −4.0e-5; end to far cone −7.4e-5 / −1.3e-5 / −5.9e-6; end
  to own cone 9.8e-4 / 9.1e-4 / 7.3e-4; cone to cone 1.7e-3 / 1.5e-3 / 9.4e-4; cone to
  itself 5.7e-4 / 5.4e-4 / 3.4e-4; cone areas exact within round-off, reciprocity to 1e-16.
  Flat facets instead: end to end −9.8e-3 / −3.1e-3 / −1.5e-3 and cone areas short by
  2.3e-3 / 7.4e-4 / 3.6e-4.
- **Validation, transport.** The frustraum with entrance holes of 0.7 mm radius and a
  0.6 mm capsule, the wall driven (albedo 0.8), on flat facets against
  `montecarlo_incident(...; reemit = :element)` (16 × 250 paths per receiver): with
  balanced view factors at h 0.3 mm the LEH loss differs by 1.2e-3 and the capsule's
  incident power by 1.9e-3 (statistical error 6e-4); the cylinder of the same length and
  holes (radius 1.6 mm) by 1.8e-3 and 2.6e-3. Unbalanced and at h 0.4 the frustraum
  stays within 6e-3, the cylinder within 1.9e-2. The energy balance (absorbed + escaped
  + lost = source) closes to 5e-15.

## Frustraum cone pointing: derived (2026-11-03) *derived*

BriareusGUI's NIF presets in a frustraum part's entrance holes (`FRUSTRAUM_CONE_POINTING`,
`preset_pointing`); derivation `examples/frustraum_pointing.jl`, test
`BriareusGUI/test/test_profiles.jl`. Sources read 2026-11-03.
- **Public pointing: the 50° cone only.** Baker et al., Phys. Plasmas 30, 092708 (2023)
  (OSTI 2222631), Sec. III A 2 and Fig. 3: the outer cones' pointing was moved along the
  axis by 0.2 mm between N200108 and N200217 and by 0.1 mm between N200721 and N201018
  (text); Fig. 3(b)'s abscissa, "50 deg. z pointing (mm)", read from the page image at
  400 dpi (±0.01 mm): N200217 and N200721 5.99, N201018 6.09, N200108 6.19. Fig. 3(a)
  sketches the 50° beams' footprint on the wall for two shots, without numbers. The text
  states that the 44.5° cone's z pointing was farther from the capsule than the 50° cone's in
  the layered shots (no value), and the 25° wall moves the outer spots by an effective
  ~0.13 mm in z. No inner-cone pointing, radial offsets or defocus are given. The shots'
  frustraum dimensions are not tabulated per shot (the paper's main geometry is the
  starter's). Baker et al.'s Frustraum 1100 paper (LLNL conference record, OSTI 2315082;
  journal version High Energy Density Phys., OSTI 2447656) gives quad splitting, not
  aims. Amendt et al. 2019 (the design paper) was not read (paywalled). Confidence in the
  50° values: medium (figure-read; the reference point of "z pointing" taken as TCC on
  the axis, which the paper does not state).
- **Derivation (engineering default).** No full public pointing, so the four cones are
  derived from the example cylinder pointing (section "Example NIF cone pointing:
  unsourced"; itself an illustrative example) on its 10.1 × 5.75 mm hohlraum: each cone's
  spot centroid (power-weighted mean height of the wall deposition of the cone's quads
  in one hemisphere) over the half-length is held at the cylinder's, on the frustraum
  starter (12.24 mm, waist radius 4.88 mm, ends 2.282 mm, holes 1.55 mm, capsule
  1.281 mm). Each cone's line is set by where it crosses the hole's plane; its aim keeps
  the cylinder aim's axial distance from the hole (so the spot size at the hole is the
  cylinder's). Where the fraction would bring the beams closer to the hole's rim than the
  cylinder pointing's smallest keep-out clearance (0.0656 mm, its 23.5° cone, by
  `ApertureKeepOut` defaults: 1% contour, 0.075 mm transverse and 1 mm focus error), the
  line stops at that clearance instead. Aims rounded to 0.01 mm.
- **Result** (aim z, r toward the port azimuth, mm, frustraum at TCC; fraction of the
  half-length, cylinder → frustraum, ray trace nspot = 121, wall mesh h = 0.2 mm;
  keep-out clearance, mm): 23.5° (4.16, −1.09), −0.187 → −0.372, 0.066; 30° (4.70,
  −1.16), −0.042 → −0.154, 0.069; 44.5° (6.27, −0.25), 0.454 → 0.449, 0.067; 50° (6.08,
  −0.06), 0.518 → 0.518, 0.229. All cones deposit their power on the wall. The 50° cone
  is at the cylinder's fraction; the 44.5° cone 0.005 short of it (0.03 mm), at the
  clearance limit; the inner cones cannot reach it. On the frustraum the 23.5° beams run
  nearly parallel to the 23° wall of their own half, so they land on the far half past
  the waist; reaching the cylinder's fractions (−0.19, −0.04) would take a hole-plane
  crossing near the rim (clearances −1.2 and −0.4 mm). The 23.5° line crosses the hole
  where the cylinder's does, so its clearance is the cylinder's exactly.
- **Check, not fitted.** The derived 50° aim, z = 6.08 mm, lies within the published
  5.99–6.19 mm (Fig. 3(b)), and the 44.5° aim (6.27 mm) lies farther from the capsule than
  the 50° aim, as the text states for the layered shots. Neither was a target of the
  derivation. Confidence in the derived table: low as a model of the platform's pointing
  (an engineering default for illustrative designs), the 50° cone consistent with
  public data.

## Visibility default on the GUI presets (2026-11-03) *measured*

Whether primitive targets should default to `visibility = :sampled` (roadmap item 3).
The GUI's default design (the hohlraum starter, flush entrance holes, capsule) and the
shock-tube starter, built by `BriareusGUI.Controllers.build_target`, curved surfaces,
4 threads; a wall source of two Gaussian rings (the CAD tutorial's), balanced per cavity
as the drive balances it. Max |ΣF − 1| before balancing (rms in parentheses):

| preset | rule | build | wall | caps | LEHs | capsule |
|---|---|---|---|---|---|---|
| coarse (1,152) | centroid | 0.34 s | 6.6e-3 (2.7e-3) | 1.3e-3 | 5.4e-4 | 3.6e-4 |
| coarse | sampled | 1.16 s | 3.8e-3 (1.5e-3) | 5.8e-4 | 9.5e-4 | 3.7e-4 |
| standard (4,796) | centroid | 3.81 s | 3.1e-3 (8.5e-4) | 7.4e-4 | 7.7e-4 | 1.3e-3 |
| standard | sampled | 16.07 s | 1.2e-3 (4.4e-4) | 7.0e-4 | 4.6e-4 | 1.3e-3 |

Sampled minus centroid, balanced (what the GUI solves): wall Tr +2.4e-6 / +2.6e-5
relative, LEH loss −5.5e-5 / −1.6e-4, capsule flux +6.4e-5 / −2.4e-6, a₂/a₀ +0.025 /
+0.003 pp, a₄/a₀ −0.010 / +0.001 pp (coarse / standard). Unbalanced the rule moves the
capsule flux by 9.1e-3 / 4.3e-4 and a₂/a₀ by −0.33 / +0.03 pp. Against a Monte Carlo of
the element model on the coarse flat facets (`montecarlo_incident(...; reemit =
:element)`, 20,000 paths per capsule element: a₀ ± 2.5e-4, a₂/a₀ ± 0.058 pp) the two
balanced rules differ from each other by less than the statistics (a₂/a₀ +0.14 and +0.16
pp; a₀ +6.0e-3 and +6.1e-3, the coarse mesh's cracks at the entrance-hole rims, which
balancing fills and the rays leak through). The shock-tube starter (standard) has no
pair a rim shadows: identical factors, 0.73 → 1.36 s.

Decision: the default stays `:centroid` for primitive targets. `:sampled` costs 3.4–4.2× the
build on the hohlraum presets (over the 2× bound set for the change), and after balancing it
does not move any result by more than 3e-5 relative or 0.03 pp. The 1–6% row sums at LEH
rims belong to bored entrance holes and gaps (the bored cavity of "Partial visibility"),
which `cad_target` and `visrad_target` already sample; on flush holes the largest row-sum
error is the wall's 0.66% at the coarse size.
