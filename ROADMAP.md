# Roadmap

## Done (v0.1, 2026-10-15)
Core pipeline: geometry/BVH, NIF + OMEGA beam sets, differentiable deposition operator +
reference ray trace, occluded view factors with exact near-field kernel, gray radiosity
with constant/tabulated/wall-loss albedo, flux response operator, Legendre/Y_lm
analysis with exact facet integration, orthographic rendering, quasi-static drive
history, grouped pointing parameterization + L-BFGS-B optimization, examples, 676 tests.

## Done (2026-10-16)
Shot-review layer:
- worst-case clearance over a linear error stack, with intensity-contour envelopes, spot
  wings and pass/review/fail tiers;
- the unconverted-light check;
- beam-table I/O;
- backscatter, calibrated transfer, and a reduced kinetic CBET model with calibration.

## Done (2026-10-16, overnight)
- **Facility data:** the LMF per-beam table, NIF quad names and bundles, and the
  measured CPP fits.
- **Radiation on true curved surfaces:** sphere, cylinder and disk sectors, with
  adaptive near-pair quadrature. View factors are 19× faster and use 165× less memory,
  with a SAH ray tree for every target.
- **Real spherical-harmonic analysis:** `mode_amplitude`, the influence matrix and
  `asymmetry_objective`.
- **Deposition operator rebuilt:**
  - It sums the best-focus profile over the f-cone directions, with per-direction
    incidence factors; the moment-matched profile had the wrong shape far from focus.
  - It refines triangles where each beam can land.
  - Per-beam power is conserved to 1e-4 unclipped and 2e-3 when clipped, on fine and
    coarse meshes.
- **Shot review:**
  - `ApertureKeepOut`, a differentiable worst-case LEH clearance for optimization;
  - `talis_check`, NIF's published target-configuration rules.
- **VisRad interchange:** SST pointing, Laser Beam Summary and power files, plus the
  benchmark kit.
- **Docs:** a Documenter site with tutorials that double as validation (laser spot,
  integrating sphere, freeing a bundle for radiography).
- **`BriareusGUI`:** a desktop window built on GLMakie + QML.
- **Hohlraum diagnostic windows** (`WallWindow`), exact cylinder-patch areas, cached
  radiosity factorizations (repeated solves 45× faster), and progress/cancel hooks on
  `clearance` and `optimize_pointing`.

## Done (2026-10-16, windows)
- **Circular windows and gaps:** `WallWindow(; diameter, gap)`. Disks are drawn on
  the unrolled wall and re-meshed inside a box stitched to the conforming grid. An
  equal-area edge polygon gives exact window and gap areas.
- **View-factor visibility on curved surfaces:** facets of a cylinder or sphere no
  longer occlude rays between points inside it. The centroid facing test is relaxed
  by the spread of true normals across a facet. A tilted triangle on the wall had been
  culling real exchanges, and capsule rows now close 10× better.
- **Thin-foil burnthrough:** `BurnthroughAlbedo` and `HeatFront` use Hammer & Rosen's
  subsonic gold heat front, integrated over each element's own T_r(t) in
  `drive_history`. `energy_balance` reports `transmitted`, and `drive_history`
  reports `window_loss`.
- **Tutorial "Diagnostic windows in the hohlraum wall":** the front is checked against
  the exact ramp solutions (2.2% and 0.7%) and the Nova foil-burnthrough data. It runs
  1.7, 2.5 and 3.6 µm windows with and without a 50 µm gap.
- **Adaptive drive histories:** `drive_history(...; transition, flux_change)` inserts
  steps through burnthrough transitions (planned from predicted burnthrough times) or
  wherever the observed flux shape moves too far. Each step's albedo iteration starts
  from the previous step's temperatures (`solve(...; guess)`), about a third fewer
  factorizations.
- **Post-burnthrough transmission (default):** a burned-through foil transmits
  `T = 4/(3τ + 2)` (gray planar Eddington with a Marshak boundary; Lester et al., Phys.
  Rev. E 114, L023201 (2026)), τ = 5 at 2.5 µm of gold and linear in areal mass, and
  re-emits the rest with the wall's albedo: T = 0.33, 0.235, 0.17 at 1.7, 2.5, 3.6 µm.
  `eddington_transmission`, `optical_depth`, `sol.transmission`; `after = 0` keeps the
  open-hole bound. In the tutorial the windows' flux-weighted ΔP2 drops from +4.1–5.3%
  (open hole) to +0.7–1.7%, and the capsule-energy cost from 12–15% to 2.4–5.5%.
- **Threaded deposition:** `deposit` hands out work units (element blocks, back-side
  pieces) dynamically; each unit owns its outputs, so results are bitwise independent
  of the thread count. Operator construction and the ray trace are balanced the same
  way.
- **Back-side sampling:** the outside of a target (clipped power) is sampled for its
  total, with Gauss rules graded toward open edges (LEH rims). Clipped power at 1e-4 of
  a beam had been 10–60% low with spot-resolving sub-triangles; it now matches
  converged ray traces to 1%, with 44% fewer operator points on the NIF hohlraum.
  Element values are unchanged bit for bit.

## Done (2026-10-17, GUI as a design tool)
- **BriareusGUI as a CAD-style tool:** the step rail is gone. Model tree (components,
  beams by cone/quad/beam, analyses with their state), Properties for the selection,
  result tabs, menus and a toolbar; any operation at any time, stale results flagged;
  coloring by a drive field computes the drive on request (as VisRad does).
- **Viewer:** own camera (`CadCamera`): zoom toward the pointer, orbit about the
  selection or the picked point, pan the grabbed point, double-click fit, standard and
  along-beam views, perspective/orthographic, near/far planes from the bounding
  sphere (no clipping), ray-cast picking, an orientation triad. Beams as rays to where
  they land, or as envelope volumes cut at the target.
- **Free parts** (sphere, disk, washer, rectangle, cylinder, cone) placed with
  VisRad's polar/azimuth/rotation angles, with their own material, sidedness,
  aperture flag and element size. Round diagnostic windows by default.
- **Beam groups** (VisRad's clusters, kept apart from the cone `group` the LPI models
  use): any set of beams, Ctrl-selected in the tree or the view, moved to a named
  group and pointed together, relative to the ports or at a common aim point; the
  optimizer's per-group parameters follow them. NIF quads split into their four
  beams and merge back. The tree carries the component and group colors (no color
  key on the canvas).
- **Capsule keep-out** (`SphereKeepOut`): the optimizer had moved inner cones onto the
  capsule (0.75% of their power) because the `direct` penalty is quadratic in a tiny
  fraction; the keep-out holds the worst-case envelope off it (on by default in the
  GUI). `optimize_pointing` also caps objective evaluations (`evaluations`), which a
  GUI run had overrun with 100 iterations of long line searches.
- **Mesh work-precision study** (`examples/mesh_convergence.jl`): element-size meshing
  (`hohlraum(; h, capsule_h)`), GUI presets coarse/standard/fine at ~0.5/0.1/0.02
  percentage points in P2 and P4.

## Done (2026-10-17, VisRad workspaces)
- **`read_visrad_workspace`** reads VisRad 20/21 `.vrw` workspaces, which are keyed
  text, not the JSON the manual describes. It reads the target components, including
  bundled CAD as generalized meshes (cm to mm, placed in the chamber frame through
  reference frames), and the beams with their pointing, 3ω power, f-number and spots.
  NIF CPP radii are recognized as HWHM by their numbers. `visrad_target` builds a
  `Target`, meshing spheres, rectangles and cylindrical hohlraums on VisRad's own grid
  (element counts checked against VisRad's power-source table). Checked on a real NIF
  CAD workspace: ports match `nif_beams` to 0.006°, and a ray trace of its 184 beams
  through all 210k CAD elements lands them in the hohlraum. CAD normals there are
  inconsistent (the inner wall faces out), hence `two_sided`.

## Done (2026-10-17, performance)
Measured on the NIF hohlraum at the GUI's standard mesh (5,864 elements, 48 quads,
12 parameters, both keep-outs), 24 threads (i9-12900K):

| | before | after |
|---|---|---|
| first objective call (compile) | 14 s, 12 GiB | 4 ms |
| gradient | 31 ms | 8 ms |
| optimizer, 20 iterations | 4.7 s | 1.05 s |
| `clearance` (48 quads) | 6.0 s, 15 GiB | 1.8 s, 0.15 GiB |
| `Scenario` setup | 2.1 s | 1.3 s |
| first `optimize_pointing` in a session | +1.5 s compile | none |

- A sum over a generator of a mixed keep-out tuple took Julia 1.13's inference 27 s.
- Objectives read the observed flux's moments from the response contracted with the
  Legendre/Y_lm weights; the full response (`scenario.response`, `ev.flux`) is computed
  on first read (8 s and 2.5 GiB on the fine mesh), and dual mat-vecs go through BLAS.
- L-BFGS-B gets value and gradient from one dual pass, under one package dual tag, and
  the package precompiles a coarse design workload.
- (2026-10-23) Gradient chunks run concurrently, each nesting the threaded deposit and
  CBET pair loop; bitwise equal to the serial chunk loop. NIF `hohlraum()`, 46 quads,
  CBET, `asymmetry_objective` + `ApertureKeepOut`, 16 threads: gradient 50 → 25 ms (46
  parameters), 148 → 72 ms (138), CPU 45 → 85%; 5 L-BFGS-B iterations 2.1 → 1.3 s.
  Chunk 12 with nested threads beat chunk 6, 24 and 46 and un-nested concurrency.
- Deposition skips points and cone directions off the spot on values (78% of points,
  exact zeros); the focal profile uses `exp(h log q)`; smaller work units.
- BVH slab tests without NaN-aware min/max (ray traces 2×); `any_hit` specialized on
  its skip predicate (view factors: 13.8M allocations to 0.6M).
- Radiosity and closure without n² temporaries.

## Done (2026-10-17, windows and CAD)
- **Window design studies** (roadmap item 2, in the diagnostic-windows tutorial):
  - `RosselandTransmission`: the post-burnthrough optical depth follows the radiation
    temperature as gold's Rosseland opacity (`T^-1.5`, Hammer & Rosen), normalized to
    Lester's point at 325 eV; Hammer & Rosen's absolute fit gives τ = 5.9 there
    (Lester: ∼5). About a tenth less window cost than the constant in the tutorial.
  - `GapClosure`: an open gap that closes at a given time; closed before the peak a
    50 µm gap costs a fifth (at 7 ns) or nothing (at 3 ns) of an open one.
  - `Scenario(...; albedo)` linearizes about a solution (e.g. a drive history's peak,
    windows burned through), reproducing its P2 exactly.
  - Restoring the drive: per-quad powers (±20%) or quad-pair re-pointing (≤ 0.3 mm)
    remove the windows' P2 and m = 2 at the peak (not their l = 4); settings found at
    the peak overshoot before burnthrough, and power trims switched on with the
    burnthrough hold the no-window shape through the pulse.
- **CAD repair** (part of item 5): `orient_mesh` (consistent windings, closed solids
  facing out of their material, open sheets toward a point) and `coarsen_mesh`
  (near-coplanar facets merged into elements of a given size); `visrad_target(...;
  orient, coarsen)`. On the NIF CAD workspace the inner hohlraum wall now faces the
  cavity without `two_sided`.

## Done (2026-10-17, GUI sessions and CAD parts)
- **Session files:** the whole design (target, windows, parts, beams with groups and
  on/off state, analysis settings, view) as a plain TOML file (`save_session`,
  `load_session!`; File ▸ Save / Open, Ctrl + S / Ctrl + O). Ports are stored as
  vectors, so a reloaded design is bit for bit the saved one.
- **Undo and redo** of every design edit (Edit ▸ Undo / Redo, Ctrl + Z / Ctrl + Y):
  snapshots of the design taken after each action (`record!`); the view is not part
  of the history.
- **STL parts:** `read_stl` (binary and ASCII) and a `:mesh` part kind (Target ▸ Add
  part ▸ Mesh from STL file): placed like any part, windings repaired, facets merged
  into elements of the part's size.
- **Part frames** (VisRad's reference frames): a part is placed relative to the chamber,
  either entrance hole (z out of the hohlraum, so one placement mirrors) or a
  backlighter foil, and follows it when the target's dimensions change; changing the
  frame keeps the part in place (`set_part_anchor!`).
- **X-ray images in the window** (View ▸ X-ray image along the view, `xray_image!`):
  the radiance seen along the 3D view's direction, framed as the view, rendered again
  as the drive is re-evaluated; the map or the image can sit beside the 3D view in a
  second canvas.

## Done (2026-10-17, beam errors)
- `mode_sensitivity` (exact per-beam sensitivities of the drive shape to pointing and
  power errors, one forward-mode pass per beam), `mode_uncertainty` (linear
  propagation of independent errors: rms per coefficient and per mode) and
  `sample_modes` (Monte Carlo check); the objectives' `uncertainty` keyword optimizes
  the expected symmetry (robust design). Tutorial "Pointing and power errors": linear
  and Monte Carlo agree to 1–2%; NIF 50 µm / 3% errors leave 0.23–0.28% rms of mode 1,
  a floor cone pointing does not lower.

## Done (2026-10-17, design workflow)
- **Retune to a reference drive:** `retune_drive` finds per-cone, quad or beam power
  histories that reproduce a reference shot's capsule drive (Y_lm shape and level)
  under power limits, linear about the reference's converged albedos
  (`drive_history(...; keep_albedo)`), optionally staying near the reference powers.
  Pulses come from a total power and a cone fraction (`cone_fraction_pulses`) or a
  CSV pulse table (`read_pulse_table`/`write_pulse_table`), and are applied per quad
  or beam (`with_pulses(...; by)`). Tutorial "Retuning a drive to a reference shot":
  bundle 31 down takes ε₁ from 0 to 2% at the peak, and retuning the other 46 quads
  brings it to 0.02%, confirmed by a drive history to 3e-5. In the window:
  Compute ▸ Retune to a reference drive.
- **Exported calculations:** `write_results` (per-element laser, radiosity and flux,
  per-component and per-beam power, Legendre and Y_lm modes, a summary) and File ▸
  Export results, with the session that produced them.
- **VisRad workspaces in the window:** File ▸ Import (components become VisRad object
  parts, beams the table) and Export, which writes the beams into a copy of the
  imported workspace (`write_visrad_workspace`: only beam rows change; round trip
  exact).
- **Frames and coordinates:** a part placed in another part's frame (followed when it
  moves, loops refused); positions in Cartesian, cylindrical or spherical coordinates
  per part; beams aimed at a component's centroid (foil, LEH, part) and following it.
- **Window:** facility selector (NIF, OMEGA, OMEGA EP, custom: presets and checks per
  facility, TaLIS on NIF only); right-click Show and Include in calculations per
  component or beam (excluded components drawn as ghosts); Shift/Ctrl multi-select
  with multi-part edits; cutaway width, azimuth and cut components; View ▸ Look from
  port; the x-ray image under Compute, unchecking it restores the previous view; menus
  sized to their entries; a stuck pan after losing focus with Shift held fixed.

## Done (2026-10-17, VisRad deposition)
- **3ω power under VisRad's red/green models:** under "Fraction of Total Beam Power" the
  power column is the total and the blue cone carries 1 − f_G − f_R of it, whether or
  not those cones are on. `read_visrad_workspace` and `write_visrad_workspace` now
  apply it (both sample workspaces use it: 0.9 TW per 1 TW beam).
- **Checked against a VisRad deposition table** (four NIF beams on a tilted 4 mm foil,
  9 × 9 elements): totals agree to 1e-8, and elements above 5% of peak to 1.4%. The
  faint wings differ by up to ~20%; see the research notes, "The faint wings".
- **`VisRadProfile`**, VisRad's beam model, as a deposition option
  (`DepositionOperator(...; profile)`, `Scenario(...; deposition_kw = (; profile))`,
  `visrad_profile(workspace)`). The super-Gaussian keeps its shape along the beam,
  with hyperbolic semi-axes, an envelope cut-off and per-plane power normalization,
  instead of being convolved over the f-cone. It is differentiable like the default
  model and tested against the closed form.

## Done (2026-10-18, CAD and STL)
- **STL interchange** (3ffcefa): `write_stl` writes a mesh, component or target as
  binary STL or ASCII with one named solid per component (normals from the winding,
  coordinates divided by `scale`, two-sided parts written once, apertures as their own
  solids unless `apertures = false`); `read_stl_solids` returns the named solids of
  multi-body CAD exports. The ASCII reader takes any case, tabs and CRLF and parses
  200k facets in 0.08 s (`%e`) to 0.3 s (full precision), from 1.2 s; binary in 6 ms;
  truncated files get a clear error.
- **CAD remeshing** (208fb44): `remesh_mesh` cuts sliver-tessellated CAD (curved walls
  as full-length strips) into compact, near-planar elements of a given size without
  moving the surface: conforming refinement with T-junctions joined, then Lloyd-relaxed
  clusters halved until each element is within `angle` of planar, within `size` of its
  centroid and has its centroid on it; areas kept to round-off.
  `visrad_target(...; remesh)`, and the GUI's STL and VisRad CAD parts. On a NIF CAD
  workspace the inner wall of the two hohlraum halves gets 2,756 elements at 0.3 mm (at
  most 0.53 mm across, against 4.5 mm slivers), 0.4 s per half.
- **Cavity selection** (5bca2e8): `enclosure(target, point)` keeps the element sides
  bounding the free space around a point. Each element gathers, through its own diffuse
  rays, the flux an isotropic source there would put on it (geometric albedo 0.9; kept
  above 5% of the cavity's flux), so the choice ignores element size and exterior
  structure cannot carry itself in. `visrad_target(ws; inside)` selects on the facets
  and remeshes only what is kept. On a NIF CAD workspace (capsule tent left out): 15k
  of 151k elements, open area 14.5 mm² (the LEHs), 3 s; with `remesh = 0.3`, 6,950
  elements in 11 s.
- **Straight cuts on CAD strips** (aa85e31): `remesh_mesh` cuts parallel long edges
  (and their strips' diagonals) on shared plane families and grows elements in bands
  about `size` long, so cuts across sliver strips are single straight chords shared by
  both neighbours (equal per-edge spacing had staggered them: a strip's diagonal took
  one point more than its sides). On a 2,000-strip tube: outline edges 275 → 132 per
  element, flat view factors 30 → 11 s (8 threads); planar strips become rectangles.
  The NIF CAD workspace and the CAD tutorial change within their stated agreements.
- **GUI on CAD, third pass** (a89155f): targets of more than 30k triangles not built
  before, and CAD capsule sphere fits, run on workers (the edit applies at once, the
  previous target is drawn until the new one is ready, the drive and checks wait for
  it). On a NIF CAD workspace: leaving a part out 733 → 145 ms, element-size edits
  0.36 + 0.35 s → 27 ms, capsule role 772 → 43 ms. A CAD part that is one cylinder or
  sphere as a whole keeps its name and gets that surface; the drive's status flags flat
  parts of more than 50 facets per element. The ~1 s gap at the end of a drive was
  GUI-thread compilation plus full model-tree resets (precompiled drive path,
  incremental row inserts): 1.21 → 0.31 s.
- **Analytic surfaces in CAD** (534c6a6): `fit_surfaces` segments a mesh into plane,
  cylinder and sphere regions (region growing, least squares refined by Gauss–Newton,
  rms-gated); remeshed CAD is fitted on its original facets (`reference`).
  `curved_components` and `visrad_target(...; remesh, fit = true)` split parts into
  `<part>/cylN`, `/sphN` and `/rest` carrying the fitted surfaces. Curved elements of
  more than two facets get a lumped quadrature, and a folded patch's centroid is moved
  onto its facets (remeshed convex CAD, e.g. a capsule, saw nothing before); primitives
  are bit for bit unchanged. A remeshed CAD sphere gives exact F and a uniform
  integrating-sphere flux (flat: 10% off); a CAD hohlraum's view factors run 24× faster
  than flat, and closed, its capsule flux matches the primitive's to 7e-6.
- **Tutorial "Radiation in an imported CAD hohlraum"** (d877790): a synthetic CAD
  assembly (thick halves with a waist gap, LEH bores, capsule, tent, mount; sliver
  tessellation, mixed windings) through multi-solid STL, `orient_mesh`, `enclosure`,
  `remesh_mesh`, `curved_components` and radiosity. At equal element size it matches
  the same cavity built from primitives (capsule flux 0.06%, a₂/a₀ 0.12 points, LEH
  loss 0.15% at 0.5 mm; 0.012%, 0.07 and 0.04% at 0.25 mm), capsule view factors equal
  solid angles within 1.2e-3, and fitted view factors cost a tenth of flat;
  `hohlraum()` loses 3.3% more through flush LEHs than bored ones. 53 s in the docs.
- **Flat view factors on remeshed CAD** (f0e67f1): elements of more than two facets
  receive through their outline when the point sees all of it (exactly the facet sum),
  otherwise plane piece by plane piece closed across the tangent-plane cut; near pairs
  emit from cells that lump the facets. Tutorial flat view factors 12.4 → 2.7 s,
  synthetic CAD hohlraum 47 → 4.1 s (fitted 2.3 → 1.3 s); closed forms within 5e-5 of
  before; primitives bit for bit unchanged. A flat CAD pair still costs ~5× a
  primitive pair (exact square–polygon clipping on planar elements would cut it).
- **CAD cones** (ac8ba2b): `ConeSurface` (radial projection at constant height, exact
  patch areas: faceted frusta tile the true cone to 4e-15); `fit_surfaces` finds cones
  (ties: plane, cylinder, cone, sphere), `curved_components` gives `<part>/coneN`;
  elements touching an apex stay flat. A closed 45° frustum's cone→cone exchange is
  +2.5e-4 curved vs +4.1e-3 flat. On a NIF CAD workspace: an 84.5° end-wall cone at
  the LEH and a 61° waist chamfer on the halves, a 77° rim chamfer on the washers.
  `curved_components(...; maxradius, minfraction)` leave near-flat and small fits flat
  (off by default; `minfraction = 1e-3`: 138 → 25 components for two halves).
- **CAD fit cost on fan-tessellated parts** (8ad51af): `fit_surfaces` visits a patch's
  neighbours through its vertices, each vertex once per ring and once per growth, in
  the same order (results bit for bit unchanged). A NIF CAD capsule part (fill tube in
  800-facet fans) fits in 0.11 s (was 0.5 s) with default kinds, 20 s (was 121 s) with
  cylinders and spheres only.
- **One call for CAD** (70b9dc0): `cad_target(parts; orient, inside, remesh/coarsen,
  fit, material, two_sided)` takes CAD meshes (`read_stl_solids` pairs) and primitive
  components to a radiation-ready target: orient, cavity selection on the facets,
  remesh only the kept facets, select again on the elements, fit cylinders and spheres
  on the kept facets; settings per part as a value, function or dictionary of the name.
  `visrad_target` uses it for workspace CAD (bitwise unchanged on the fixture), and the
  CAD tutorial ends with it.
- **GUI on CAD workspaces** (ffef88c): VisRad import on worker threads with progress;
  far mounting hardware drawn only (raw facets), pickable, and put back by *Include in
  calculations*; the drive of a CAD design runs on the cavity around a point
  (`enclosure`), and a pocket sealed by a capsule tent is opened automatically (the
  tent left out as a thin film, said in the status); view factors past half the
  memory fall back to a laser-only drive with the reason. Part roles (Capsule, a sphere
  fitted to a CAD capsule; Backlighter foil; Thin film), File ▸ Export target as STL,
  one part per solid on STL import, renamed VisRad parts rename their components,
  cutaway on imported designs, off beams hidden when most are off, 1 ms Windows timer
  while the window is open (steady 17 ms ticks). On a NIF CAD workspace: 1.7 ms
  frames at 211k triangles, import off the UI thread (0.4 s call), 9.8k-element cavity
  drive in 92 s, capsule ε₂ 15.4% → 9.2% in 40 evaluations (14 s).
- **GUI on CAD, second pass** (53811fa): VisRad primitives take a part's element size
  (`visrad_component(...; primitive_h)`; VisRad's grid by default). CAD remeshing after
  element-size edits, undo/redo and session loads runs on workers, the part drawn from
  its facets until ready and a queued drive starting then (all ten CAD parts: 1.97 s
  UI stall → 0.33 s). Targets are reused for identical components (undo/redo 1.1 →
  0.05 s, the drive stays current). The 0.7–0.9 s first-use stalls (drive result,
  field coloring, import) were JIT compilation, now warmed before the window opens.

## Done (2026-10-21, VisRad deposition parity)
- **Checked against VisRad 21 runs on a NIF hohlraum** (thor `visrad/workspaces`:
  n210808, n210808_reference, bu34_repointed; 192 CPP beams, 2388 wall triangles).
  With the manual's z_A the wall elements differed by 10% rms and Briareus's peaks were
  10% high: VisRad's spots grow faster away from focus. `VisRadProfile` now scales z_A
  by γ/ρ (γ = `growth` = 0.67, ρ the beam's envelope ratio), inferred from n210808 and
  checked on the other two: wall 2.7–2.8% rms, correlation 0.9997, peaks within 0.5–2%;
  bu34_repointed's capsule departures δ1..4 0.007, 0.090, 0.019, 0.037% against
  VisRad's 0.006, 0.092, 0.018, 0.046% (were 0.006, 0.054, 0.012, 0.024).
- **Beam-profile probe** (thor `visrad/profile_probe.vrw`: 48 beams onto plates normal
  to them, |z| = 0–15 mm). VisRad's spot rotation is relative to the quad's horizontal,
  measured to 0.02° and now applied by the readers and writers. n(z) grows linearly
  (manual form, z₀ ≈ 18.5/12.7 mm inner/outer). The semi-axes are tabulated per NIF
  plate (`measured = true`). Results:
  - hohlraum wall rms 0.9% (from 2.8%);
  - bu34_repointed δ1..4 0.008, 0.104, 0.020, 0.045% against VisRad's 0.006, 0.092,
    0.018, 0.046%.
- **Envelope per beam.** VisRad's default ratio is (ln 1000)^(1/n) (`visrad_envelope`;
  1.27–1.54 for the NIF plates). `VisRadProfile(; envelope)` takes a number, a `Dict`
  by beam name, or `nothing` for the default; `visrad_profile(ws)` reads each beam's.

## Done (2026-10-22, correlated beam errors)
- **Shared errors in `mode_uncertainty` and `sample_modes`** (roadmap item 11):
  - `group_pointing`: one displacement shared by a group's beams, transverse to the
    group's mean direction;
  - `group_power`: one relative power error per group;
  - `groups`: a beam field (`:quad`, `:group` for the cone) or a function, e.g. the bundle;
  - `positioning`: the target offset from chamber center, isotropic or per axis.
  Shared displacements move beams along themselves too, so the per-beam sensitivities
  gain an axial column (computed only when shared errors are asked for). `rms_by` splits
  the budget by source. The robust objectives (`symmetry_objective`,
  `asymmetry_objective`) take the same keywords in `uncertainty`. Tests: a one-beam
  group reproduces the independent errors to
  1e-9, a target shift matches central differences, and every source together matches
  Monte Carlo within 15%. Tutorial "Beam errors": 20 µm per half-cone gives ~40% of
  the l = 1 of 50 µm per beam. A shared displacement needs a group with a common beam
  direction: a cone over both hemispheres (`groups = :group`) is refused.

## Done (2026-10-22, GUI beam model)
- **Beam-model choice in the GUI** (Drive ▸ Beam model): Briareus's f-cone or VisRad's
  `VisRadProfile`, passed to the drive's `Scenario` (`deposition_kw`), so the
  optimization and the VisRad kit use it too. Importing a VisRad workspace sets
  VisRad's with `visrad_profile(ws)`; switching back keeps the parameters. Saved in
  session files (older sessions load as f-cone), undoable, and a change makes the
  drive stale.

## Done (2026-10-22, VisRad results)
- **VisRad's deposition read back:** `visrad_power_sources` reads a computed workspace's
  `[Power Source Data]` tables (columns named from the file's header; an error, or
  `nothing`, when not up to date), and `visrad_deposition(ws, target)` maps the direct
  laser source onto a `visrad_target` as absorbed TW per element, the convention of
  `deposit(...).element`. CAD meshes map in file order and rectangles local-x fastest;
  spheres and hohlraums with laser power are refused (order unknown). On N210808 it
  equals the ad hoc study reader bit for bit (wall 428.4920, THOR windows 6.2526 TW each).
  `write_visrad_workspace` resets a computed template's results (power sources,
  radiosities, capture/clearance) to "not up to date" when it writes any row, so old
  results cannot be read back as the new beams'.

## Done (2026-10-22, inverse beam errors)
- **`infer_errors`** (roadmap item 11): the most probable beam, group and positioning
  errors behind a measured shape departure (components or whole degrees, with a
  measurement sigma), under `mode_uncertainty`'s priors, by the linear-Gaussian
  posterior in data-space (Woodbury) form. Returns the errors (per beam, per group, the
  target offset), posterior and prior rms, the predicted departure and residual, the
  sources ranked by the share they explain, the number of combinations resolved, and
  the beams with the errors applied. Tests: a target offset recovered from a full
  evaluation's l = 1, 3 to ~1% (within 3σ), a hemisphere cone's shared displacement
  against every cone's, MAP and covariance equal to the dense error-space solve to
  1e-8, and full evaluations of the inferred beams agree to second order. Tutorial
  "Beam errors": a measured P1 from a 30 µm target offset gives 29 ± 3 µm under a
  positioning-only prior, but 7 µm (l = 1) or 11 µm (l = 1, 3) when beam and cone
  errors may share the blame; the median error keeps 99% of its prior spread.

## Done (2026-10-22, compact capsule mesh)
- **`icosphere_mesh`** (geodesic: icosahedron faces divided n × n, 20n² near-equilateral
  triangles) and `hohlraum(; capsule_mesh = :compact)` sized by `capsule_h`
  (`icosphere_frequency`); lat-long stays the default. `examples/capsule_mesh.jl`: P4
  within 0.02 pp at ~1k elements against ~8k lat-long, both converging as 1/N to the
  same limits. Tests: closed, outward, exact curved area, convex faceting (F_ij = A_j/4πR²),
  uniform field has no l ≥ 1 moments, hohlraum capsule drive and P2 match lat-long.

## Done (2026-10-22, partial visibility)
- **`view_factors(...; visibility = :sampled, visibility_grid = 4)`**: probe rays flag
  pairs that straddle an occluder; those get a kernel-weighted visible fraction from
  g² × g² sample rays, the rest keep the centroid factor bit for bit. Bored-LEH cavity at
  0.5 mm, largest row-sum error before balancing: walls 0.97% → 0.19%, LEH apertures
  5.9% → 0.7%, capsule 0.16% → 0.05%; 2–4× the time. Fractions within 0.02 (mean) of a
  256²-point reference. Tests: default unchanged, reciprocity, half-shadowed pair = ½.

## Done (2026-10-22, Monte Carlo radiosity)
- **`Briareus.montecarlo_incident`** (unexported): backward diffuse path tracing on the
  facets, no view factors, re-emission from the hit point (continuum) or the hit element
  (element model with exact factors); batch errors. `examples/mc_radiosity.jl` settles
  partial visibility (Open question 9): sampled is closer after balancing on the bored
  cavity at 0.5 and 0.25 mm; balancing biases a₂/a₀ by +0.06 to +0.48 pp on flat
  elements. Tests: closed cube S/(1 − α), open cube against closed-form factors,
  continuum limit of an 8 × 8 cube.

## Done (2026-10-22, rounded hohlraum corners)
- **`hohlraum(; corner_radius, nr_corner)`**: quarter-circle wall–cap fillets as frustum
  bands with exact `ConeSurface`s (the arc's chord polygon; about one wall cell of arc
  per band), wall and caps shortened to meet them vertex for vertex; `corner_radius = 0`
  is the square hohlraum element for element. Tests: areas (wall and caps exact, the
  fillet just under the quarter torus), facing, raw row sums, D2 selection rules, and a
  steep beam on the fillet (operator vs ray trace). For N210808 (user: 1.053 mm).

## Done (2026-10-23, overnight: hohlraum dynamics, implosion metrics, N221204)
Details and numbers in the research notes, under the section named in each entry.
- **Implosion metrics** (item 11; "Implosion metrics"): `src/implosion.jl` turns mode 1
  into a hot-spot velocity (100 km/s per %, away from the bright side, or a user Green's
  function over a drive history), with Hurricane's piston-model yield factors (no-α and
  α), Nora's burn amplification (low branch, an upper estimate), P2 swing metrics against
  Milovich's 7.5 % picket limit, the velocity covariance and expected yield factor from
  `mode_uncertainty`, `infer_errors` from a measured nToF velocity with Casey's
  capsule-thickness nuisance, and a velocity penalty for the objectives. Young's N170601
  gives 44.0 km/s (44.7 VF + G) and Kritcher's N230729 70 km/s (70–75). Tutorial
  "Implosion metrics".
- **Moving wall** (item 9; "Moving wall on keyframes"): `WallMotion`, `moving_hohlraum`
  and `Keyframes`. The wall cylinder shrinks at 0.10 mm/ns by default; the exchange
  matrix AᵢFᵢⱼ is interpolated between keyframes (reciprocity kept) and deposition
  operators are built at each step time. The error falls as the square of the keyframe
  spacing (11 keyframes: 0.012 pp in P2/P0 against a rebuild at every step), and six
  fine keyframes add ~14 % to a fine history. A still wall reproduces the fixed target to
  1e-12. `drive_history` runs on keyframes.
- **Closing entrance holes** (item 9; "Closing entrance holes"): `hohlraum(; leh_closure,
  leh_open_radius = 0.84 r₀)` splits each LEH into an open hole and a lip carrying a
  `GapClosure` albedo; transparent to beams, booked as escaped until it closes, the
  keep-out keeping the full rim. In the Monte Carlo reference, `write_stl` and the GUI.
  On the retune tutorial's hohlraum, closing at 5 ns raises Tr 3.3 % and P2/P0 2.5 points.
- **Inner-beam bubble absorption** (item 9; "Inner-beam bubble absorption"):
  `BubbleAbsorption`, an LPI model (composes with CBET). Izumi's bubble tip, each inner
  beam's own clipped fraction at the outer footprint's plane, a transparency fitted to
  both of Izumi's anchors (exact), the absorbed power moved to that hemisphere's outer
  beams; a stagnation depth and a ridge stop (power deposited on a wall band after the
  spots separate). Differentiable. Applied by the new `lpi_deposit` in `drive_history`,
  `Scenario`, the beam-error sensitivities and `influence_matrix`.
- **Mock HYBRID-E pointing** (item 16; "Mock HYBRID-E pointing"):
  `hybride_mock_pointing()`, seven numbers fitted to the symmetry of a non-public export
  (only summary agreement in the repository): capsule P2/P0 within 0.14 pp (foot) and
  0.03 pp (peak), Y_lm (l ≤ 6) within 0.06 pp, CBET transfer within 0.03 %, LEH beam
  overlaps within 7.8 %.
- **Two-zone CBET** (item 10; "Two-zone CBET calibrated to PRE 109 Fig. 9"):
  `TwoZoneCBET` (outside the LEH, supersonic, inner → outer; inside, Mach 0, outer →
  inner), so Δλ = 0 is not zero transfer, and a multi-epoch `calibrate_cbet` fitted to
  the post-transfer cone fraction against time of PRE 109 Fig. 9 (digitized,
  `n221204_fig9()`): rms 0.016, outer-cone ratios matching Fig. 9a without being fitted,
  dCF/dΔλ 0.16 per Å (target 0.15–0.18). `n221204_cbet()` holds the fit.
- **N221204 time-dependent benchmark** (item 16; "N221204 time-dependent benchmark"):
  `examples/n221204.jl` adds wall, lips and bubble in turn on the mock pointing with the
  calibrated CBET, against PRE 109 Fig. 10. The swing has the right shape, mostly from
  the bubble, but the minimum is shallow (−0.09 at 5.5 ns against −0.17), the curve turns
  positive about 0.5 ns early, and it ends about +0.07 high (+0.17 against +0.09 at 8 ns).
  Capsule Tr peaks at 305–312 eV against 313. Nothing was tuned to Fig. 10.
- **Concurrent gradient chunks** (see the 2026-10-17 performance section): gradients ~2×
  faster, bitwise equal to serial.
- **Faster albedo iterations** (`solve`, `drive_history`): GMRES preconditioned with
  the cached LU of a nearby albedo replaces a fresh LU per iteration. A 50-step fine-mesh
  history went from 32 to 3.6 s/step (standard mesh: 3.1 to 1.3), matching the old path
  to ~1e-12.
- **VisRad kits for item 1** (thor `visrad/kits_2026-10-23`): six cases (tilted-plate
  probe, normal foil, focus sweep, envelope 2.0, steep foil, shadowing) with a comparison
  script. Waiting for runs in VisRad.

## Done (2026-10-24, LEH shield)
- **`hohlraum(; sombrero)`, `Sombrero`** (item 2, first half; research notes "LEH shield
  (sombrero)"): a flat annular shield around each end, `:sombrero_top` and
  `:sombrero_bottom`, one-sided (facing away from the waist), exact `DiskSurface`,
  `BLACK` by default. Defaults measured from a THOR-class symcap CAD: bore 4.000 mm, rim
  9.200 mm, a 39 µm sheet whose mid-plane lies 0.623 mm inside the LEH plane. It is flat
  on average, with a hexagonal array of 0.21 mm dimples (±0.12 mm about the mean plane)
  that the model leaves out. The shield sits beside the wall, so it sees nothing of the
  cavity: zero view-factor rows and columns, `closed` balancing and the solution
  unchanged. Off by default. It blocks beams and lines of sight and is `clearance`
  hardware. Cone-averaged NIF drive pointing clears it by ≥ 2 mm (core 1e-4 envelope);
  only the engineering wings' scatter floor reaches it (~1e-4 per beam, `:review`).
  Tests: areas, facing, zero exchange, beams on both faces (operator vs ray trace),
  clearance pass and fail.

## Done (2026-10-24, energetics models)
Research notes, "Wall material, conversion and capsule energetics (2026-10-24)". Item 16.
- **`du_wall_loss()` / `uranium()`**: `WallLossAlbedo(; factor, factor_exponent, liner)`
  scales gold's loss by f (T/300 eV)^-0.2 (Rosen 2005, Schein 2007), gold while the heat
  front is inside a 0.7 µm liner (Döppner 2015). f = 0.785 *calibrated* to Döppner's
  6.5 % peak-power equivalence on the high-foot hohlraum; checked on his T-1 Au/DU pair
  (+4.9 eV against +5 ± 6). Worth 8.4–9 % of power on HYBRID-E. Gold unchanged, bit for bit.
- **`Conversion`**: η(t, I) = η₀(I)((t − start)/1 ns)^s, η₀ constant or a log-I table
  (Dewald 2008); `nif_conversion()` = 0.7 t^0.12 ≤ 0.95 (0.85 / 0.90 at 5 / 8 ns, Olson
  2012). Scalar conversions bit for bit; generic in the number type; `flux_response(...;
  laser)` for intensity tables. Checked against the wall-dominated closed form
  T ∝ t^((0.38 + s)/3.3); Lindl's open-LEH Nova t^0.15 is not reproduced (0.086).
- **`CapsuleAlbedo`** (α_c, effective radius schedule) and `TemperatureAlbedo`;
  `energy_balance` books E_cap; `drive_history` records laser, source, absorbed and
  setup power, and **`energy_budget`** integrates them.
- **`examples/n221204.jl`** takes `wall=`, `conversion=`, `capsule=`, `capsule_albedo=`,
  `degrade=`, `backscatter=` and `energetics=new`, and prints the energy budget and the
  rms against Fig. 10. Before → after (standard mesh): split 57.6 / 22.2 / 20.1 → 50.8 /
  33.0 / 16.2 % (HYDRA 58.6 / 25.2 / 16.2), capsule Tr 301.8 → 334.2 eV; Fig. 10 (blind):
  minimum −0.168 → −0.134 (−0.170), crossing 6.26 → 6.31 ns (6.67), +0.123 → +0.046 at
  8 ns (+0.095), rms 0.036 → 0.044. The example's defaults are unchanged.

## Done (2026-10-24, hohlraum dynamics calibrated)
Research notes, "Hohlraum dynamics calibrated to published data (2026-10-24)". Items 9,
10 and 16; PRE 109 Fig. 10 used only as the blind check.
- **Anchors digitized as data**: Young 2020 Fig. 3a (the LLNL dynamic view-factor radii,
  `young2020_radii()`), PRE 109 Fig. 1a (HYBRID-E wall edge at 6 ns,
  `pre109_fig1a_wall()`), Izumi 2021 Figs. 8 and 10a (bubble trajectories, HYDRA absorbed
  fraction; vector graphics, `izumi2021_bubble()`, `izumi2021_absorption()`).
- **`WallMotion` saturates**: ballistic, then stopping smoothly. Defaults fitted to
  Young: 0.185 mm/ns from 0.8 ns after the turn-on, stopping 0.45 mm in (HYBRID-E: 0.38,
  the Fig. 1a waist edge). The motion happens in the foot; the old 0.10 mm/ns from the
  turn-on moved N221204's wall 0.86 mm by 9 ns.
- **`BubbleAbsorption` refitted**: speed 0.36 mm/ns and launch 0.25 ns after the turn-on
  (Izumi Fig. 8a, 28 points); density exponent 1.8 from Izumi's Eq. 7; opacity 0.111 to
  the HYDRA absorbed fraction (0.02 / 0.23 / 0.55 at 4 / 5 / 6 ns); the collision time
  scales with the capsule–wall gap and picket energy (`capsule_radius` replaces
  `capsule_scale`), reproducing all three of Izumi 2018's separation times. HYBRID-E's
  wider gap gives t_c 6.09 ns (the old capsule-only extrapolation 5.41).
- **N221204 benchmark rerun** (`examples/n221204.jl`, new `dynamics` and `shot`
  arguments): minimum −0.168 (Fig. 10 −0.173), +0.123 at 8 ns (+0.095), rms over
  2–8.5 ns 0.035 (was 0.051); the upward zero crossing stays ~0.4 ns early (6.26 against
  6.67 ns). Proxies of N220919 and N210808 (Δλ 2.5 / 1.8 Å).
- **The 7.4 keV CBET temperature is a degeneracy, not data**: a refit with Te ≤ 5.5 keV
  fits Fig. 9b as well (rms 0.016) with lower multipliers, but moves the Δλ null from
  1.2 to 0.5 Å (NC24: 1–1.55). `n221204_cbet()` is unchanged.

## Done (2026-10-24, spectral radiosity and a synthetic Dante)
Research notes, "Spectral radiosity, M-band and a synthetic Dante (2026-10-24)". Item 13
and item 16's first absolute check.
- **`spectral_solve`**: multigroup radiosity about a gray solution (Planckian source split
  at each element's brightness temperature, group albedos, the rest of the gray albedo
  lost or re-emitted thermally). One group is the gray solution exactly; equal albedos or
  `thermalize` keep the gray totals (tests to 1e-11, dense-solve checks).
- **Laser-produced Au M-band** (`MBandEmission`, `nif_mband`, `mband_power`): per-cone
  conversion 0.20 outer / 0.042 inner (Dewald 2008/2020). On Dewald 2020's setup the
  M-band albedo 0.5 (thermalized) gives his F_M ≈ 2, and then Dante 1 sees 12 % above
  1.8 keV (measured 12 %); the capsule M-band is pole-hot by +0.15 to +0.21 relative to
  the thermal drive (Dewald ≈ +0.21 at the ablation front), `examples/mband_dewald2020.jl`.
- **Synthetic Dante** (`dante1`, `dante2`, `dante`, `record =` in `drive_history`):
  visible projected area through the LEH by ray casting, lips at their open fraction,
  Tr over the initial and the open LEH area, Gaussian channels at Rubery 2023's Dante 1
  energies. N221204: Dante 1 peaks at 286 eV (initial area) / 306 eV (open area) at
  6.25 ns against 313 measured at 8.4 ns; capsule 302 at 5.5 ns. The drive falls through
  the peak where the data rise: item 16's late-time energetics.

## Open questions (from the 2026-10-18 CAD work)
Choices made provisionally overnight or left open; the lean is in brackets.
**Decided (user, 2026-10-22):** the bracketed leans stand. On 1: parts may be split
for computational efficiency, but GUI display and interaction must stay elegant (show and
pick the original parts, not `<part>/cylN` pieces).
1. `fit = true` by default when `remesh` is given (cad_target, visrad_target, GUI)?
   Up to 24× faster view factors and exact curved walls, but parts are split and
   renamed `<part>/cylN`, `/coneN`, `/sphN`, `/rest`. [Yes eventually; off now.]
   **Looked at 2026-10-24, still off** (research notes, "Fitting CAD by default"). On
   the NIF CAD cavity the gain is small: view factors 1.0–1.5× faster (the 24× is a
   sliver tube fitted whole, which the GUI already does), the halves' flux 0.8 → 0.2%
   from a Monte Carlo reference, fits 7 s against ~220 s of selection and remeshing.
   The cost is the names: 10 parts become 68 components, and every name-keyed surface
   sees the pieces: `Target`/`Enclosure` display, `energy_balance`, `component_power`
   and `<name>_outside`, `write_results` (components.csv rows, elements.csv's
   component column), `report`, clearance/TaLIS/unconverted rows; `elements(t, :part)`
   and `Scenario(...; observe = :part)` throw; the `direct` penalty silently switches
   off when a named part is split; `SphereKeepOut` takes every fitted sphere (9 here,
   R 0.94–6.4 mm, not only the capsule). `visrad_deposition` and
   `write_visrad_workspace` are unaffected (remeshed targets are refused; only beams
   are written). The GUI keys its tree, picking, colors, hidden set, exclusions, films,
   the cavity's sealing film and saved sessions on component names. Plan, in order:
   (a) a part name on `Component` (`part`, default `name`; `curved_components` and
   every rebuild keep it) and `parts(t)`; (b) `elements` and the name-taking APIs
   (observe, `power_fraction`/`direct`, keep-outs, clearance hardware, TaLIS) accept a
   part for all its pieces; `SphereKeepOut` defaults to whole-part spheres; (c)
   reports aggregate pieces to parts (`energy_balance`, `component_power`, CSVs,
   `report`, `show`), pieces reachable by name; (d) a part fitted as one region keeps
   its name (as the GUI's `_whole_surface` does); (e) the GUI keeps its part-keyed
   state and maps pieces back to parts for drawing and picking (or fits only the
   drive's computational target, mapping elements back as the cavity's `index` does).
   Then flip.
2. ~~Default `minfraction = 1e-3` in cad_target/visrad_target/GUI fits?~~ **Done
   2026-10-24:** `curved_components` (so every fit) defaults to `minfraction = 1e-3`.
   Two NIF CAD hohlraum halves: 138 → 25 components, 1.2% of the curved area left flat;
   the NIF cavity 158 → 68, capsule flux within 3e-4 of fitting everything.
   `maxradius` stays off.
3. Cones on by default in `fit_surfaces` (22 more components on the NIF halves)? [Keep.]
4. Capsule tents and LEH windows in CAD are opaque closed solids that seal the capsule
   in a pocket. The GUI opens it automatically (tent as a thin film); the library
   does not. [Library stays explicit; add a transmitting thin-film option.]
5. `enclosure` threshold 0.05 (6.9k remeshed elements on the NIF workspace) or 0.01
   (1.4×, return estimate 0.8 → 0.12 mm² against 14.5 mm² open)? [0.05.]
6. GUI defaults: far CAD parts drawn only on import [keep; reversible per part]; the
   foil objective with per-beam parameters [yes]; the cavity for non-CAD designs [no];
   the 30k-triangle background-build threshold [keep]; a one-frame colorbar flash at
   window open to warm it [accept]; partial curved fits of CAD parts in the GUI (needs
   a part-to-component name mapping or a per-element surface) [wait for a workspace
   that benefits]; the whole-part fit guard (`_on_one_surface`) [keep].
7. `coarsen = x` route to remeshing, or deprecate `coarsen_mesh` for CAD? [Keep both;
   recommend `remesh` in the docs.]
8. Remeshed planar strips are now 2–4-facet rectangles on the primitive view-factor
   path, as accurate as primitive quads (test/cadfit.jl Hottel check loosened 1e-3 →
   2e-3). The alternative, a finer near rule for ≤ 2-facet elements, moves every
   primitive baseline. [Accept.]
9. Baseline moves to offer as options first: a compact capsule mesh in `hohlraum()`
   (offered 2026-10-22 as `capsule_mesh = :compact`; make it the default?) and
   partial-visibility sampling for pairs whose centroid ray grazes a rim (offered
   2026-10-22 as `visibility = :sampled`; make it the default?). With bored LEHs it cuts
   row-sum errors before balancing 5–9× on walls and LEHs (not the 20 µm gap's faces)
   for 2–4× the view-factor time, and after balancing moves a₂/a₀ by −0.16 pp and a₄/a₀
   by +0.11 pp at both 0.5 and 0.25 mm (research notes, "Partial visibility"). Checked
   against a Monte Carlo radiosity on the same facets (flat elements; research notes,
   "Monte Carlo radiosity reference"): the shift is toward the reference in every bored
   case, the balanced a₂/a₀ error going 0.51 → 0.31 pp (0.5 mm) and 0.30 → 0.14 pp
   (0.25 mm), a₄/a₀ 0.27 → 0.13 and 0.11 → 0.01 pp; on flush-LEH `hohlraum()` both
   visibilities balance to the same answer. What remains is mostly balancing's own
   bias (+0.06 to +0.48 pp in a₂/a₀ over the raw factors, the control included), a
   separate question for the `closed` default that needs a curved-surface reference
   (see 11).
   The geodesic capsule reaches P4 within 0.02 pp with ~0.8–1.2k elements against
   ~8k lat-long, same limits, same P2 and drive (research notes, "Compact capsule
   mesh"). Switching moves every capsule baseline and the GUI presets' meaning (a
   compact preset would be capsule_h ≈ 0.12). **Decided (user, 2026-10-22):** compact
   is the `hohlraum()` default, GUI presets re-chosen (capsule_h 0.3/0.16/0.09: P2 also
   needs ~2,500 capsule elements for 0.02 pp); `Target` carries a `visibility` rule,
   `:sampled` for `cad_target`/`visrad_target`, `:centroid` for `hohlraum()` and
   hand-assembled targets (pass `Target(...; visibility = :sampled)` for primitive
   cavities with bores).
10. `orient_mesh` faces a hollow solid's inner surface into the material; containment
    parity would face it into the gas, but `enclosure` from TCC would then select a
    hollow capsule's interior (the GUI relies on the current behavior). [Keep;
    document.]
11. `closed = true` the default for `hohlraum()`? Against a continuum Monte Carlo
    (`examples/mc_closure.jl`; research notes, "Balancing on the curved path"), on the
    curved default path with NIF quads at α = 0.925 and 0.85, balancing moves a₀ and
    a₂/a₀ toward the reference at every preset (study mesh 0.45/0.15: a₀ +0.56% →
    −0.07%, a₂/a₀ +0.16 → +0.09 pp; fine: within 0.03% and 0.02 pp); departures agree
    within 0.01 pp raw or balanced. The flat-element bias does not carry over.
    Unchecked for bored LEHs, gaps and CAD cavities. [Keep the default.]
12. LEH power-balance tutorial on the compact capsule (2026-10-22): against its photon
    walk, the error grows from 0.4 to 0.25 mm (LEH −1.2e-4 → −3.3e-4, capsule −2.8e-4 →
    −5.0e-4, 3–5σ) as the view factors' lost fraction grows (7.5e-4 → 1.4e-3); on the
    lat-long capsule it fell to the sampling error. Look at the capsule's raw row sums at
    capsule_h = h/3 (0.083 mm) before trusting balanced results at that resolution.
13. **Bug, found 2026-10-24: `:sampled` visibility on remeshed CAD** (the
    `cad_target`/`visrad_target` default). On the NIF CAD cavity at 0.6 mm the outer
    capsule's row sums are 0.21 (0.0–0.99) against 0.955 with `:centroid`, and its flux
    7× too low (`montecarlo_incident` agrees with `:centroid` to 0.03%); the halves'
    rows 0.79 against 0.91, flux −15%. Fitting does not help. `_visibility_samples`
    averages each grid cell's sub-triangle points; on an element folded across a
    convex part's facets that mean lies behind them (research notes, "Fitting CAD by
    default"). Put each sample on the element's facets (as `Target` does the centroid)
    and add a CAD capsule to the partial-visibility tests. The GUI builds `:centroid`
    targets and is unaffected.
    **Resolved (2026-10-24).** Cause confirmed on a synthetic CAD cavity (the test
    suite's `cad_half` halves at 0.6 mm, a 56 × 112 lat-long `cad_sphere` capsule at
    0.3 mm, 48 facets per capsule element): 91% of the capsule's samples sat up to
    0.9 µm behind its facets (ray offset 0.55 nm) and were blocked along their own
    normal; no centroid was. `_visibility_samples` now moves each sample of an element
    of more than two facets along the normal onto the facet that line meets (or to the
    nearest sub-triangle point of its cell), and such an element's own facets no longer
    occlude its sample rays (with the samples on the facets, every remaining blocked
    ray of the reduced capsule rows was the emitter's own fold toward its horizon,
    which the element model, radiating from the mean plane, already treats). Primitive
    `:sampled` factors are unchanged bit for bit (bored cavity, `hohlraum()`).
    Synthetic, against `montecarlo_incident` (`:element`, 8,000 paths): capsule rows
    0.027 → 0.9991 (`:centroid` 0.9986), capsule flux −97.5% → −0.33% (`:centroid`
    −0.45%, σ 0.06%), driven wall −0.32% (`:centroid` −0.41%). NIF CAD cavity at 0.6 mm:
    outer capsule rows (its 259 outward elements) 0.16 → 0.961 (`:centroid` 0.956), flux −86% → −0.31%
    (`:centroid` +0.03%, σ 0.05%), halves rows 0.79 → 0.915 (0.910), flux −15% →
    −0.67% (`:centroid` −0.81%), driven inner wall −0.13 to −0.29% by band
    (`:centroid` +0.33 to +0.39%: its capsule figure is a cancellation of high walls
    and low capsule rows); view factors 56 s against 27 s. Test:
    test/partial_visibility.jl, "sampled visibility on remeshed CAD".

## Next
1. **VisRad deposition parity: needs runs in VisRad.** The goal: optimize in
   Briareus with `profile = visrad_profile(ws)`, export with `write_visrad_workspace`,
   and get the same deposition in VisRad. Near focus that holds today. What VisRad
   does away from focus, and at envelopes, shadows and tilted surfaces, needs VisRad
   output, all on the foil workspace unless noted. Kits for the tilted-plate probe, the
   normal foil, the focus sweep, the envelope, the steep foil and shadowing are built
   (thor `visrad/kits_2026-10-23`, 2026-10-23) and wait for runs in VisRad:
   - **Intersection-point export** (File ▸ Export, laser beam–target intersection
     points, with "write beam spatial profile parameters" checked in Preferences ▸
     Lasers; VR/files/laser_target_intpts_file.html). It lists a(z), b(z), n(z),
     f_SQ(z), ΔL_W(z), ΔL_H(z) and M(z) at each hit. On a hohlraum workspace it would
     confirm or replace the inferred growth γ/ρ (a(z) should come out ~6% above a₀ at
     4.5 mm for the 50° cone), say whether it follows the envelope ratio or the
     exponent (rerun one cone with a non-default envelope), and give z₀ and f_SQ.
     Remaining after the 2026-10-21 profile probe:
     - the law behind the tabulated a(z), b(z), so plates other than the four NIF ones
       work;
     - δ2 0.012 pp high on bu34_repointed, almost all from the wall;
     - totals: wall 0.07% low, tilted foil 0.37% high;
     - the foil wings.
     A rerun of the probe with a non-default envelope or lens ratio would show which
     of them set the growth.
   - **Foil normal to the beams** (tilted to face the quad's mean port). With the
     distance from focus nearly constant across the spot, this separates the focal
     profile from the defocus law.
   - **Focus sweep**: focus offsets of ±1, ±3 and ±10 mm on the same foil. Hohlraum
     walls sit this far from focus, and that is where VisRad's shape-keeping profile
     and Briareus's f-cone convolution differ most.
   - **Envelope**: the same run with the envelope ratio changed (e.g. 2.0), on a
     finer foil (e.g. 41 × 41). VisRad's effective cut-off lay near 1.64 1/e radii
     against the file's 1.53645 (= (ln 1000)^(1/4.5), the 10⁻³ radius for n = 4.5).
     The ratio's exact definition is open.
   - **Tilted-plate probe** (todo, 2026-10-22): plates tilted to cos i ≈ 0.6–0.8 about
     4–5 mm from focus, where the outer cones land on the wall; build it with thor
     `visrad/build_profile_probe.jl`. On N230729 at 0.1 mm (thor `n230729_fine`) the
     tuned model is 1.16% rms per wall element (peak 599 vs 609 TW/cm²), almost all
     on the 44.5°/50° shoulders (model −3 to −16% at r/ρ 0.6–0.9). It is not
     sampling (centroid vs integrated agree), and per-cone refits to that run
     degrade the held-out runs and the normal plates (thor `visrad/refit_profile.jl`).
   - **Power on tilted surfaces**: VisRad's foil total was exact to 1e-8, while
     per-plane normalization gives 1.1e-4 too much at ~40°. A steeper foil (~70°)
     would show whether VisRad renormalizes each beam's deposited total.
   - **Shadowing**: a washer or sphere in front of the foil, and a hohlraum case with
     LEH clipping. VisRad decides occultation on its 25 × 24 ray grid; Briareus uses
     the beam axis through each point.
   - **End-to-end**: a primitive hohlraum with a few quads, optimized here, exported,
     and run in VisRad, with the workspace saved after computing. The
     `[Power Source Data]` tables read back with `visrad_power_sources` and
     `visrad_deposition` (see Done), which makes it a routine comparison on CAD and
     rectangles; a primitive hohlraum needs its element order measured first (e.g.
     one beam on the wall), or the hohlraum given as CAD.
2. **Backlighter geometry.** The LEH shield is done (`hohlraum(; sombrero)`, see Done
   2026-10-24). Left: the backlighter foil's placement relative to the diagnostic line
   of sight. The line-of-sight and clearance checks in the radiography tutorial, with
   the shields in the target, then decide which bundles can drive the foil.
3. **Diagnostic windows: what remains.** The studies are in the tutorial (see Done).
   Open: the expanded foil's density (held at Lester's 1 g/cm³ through the
   normalization), a backing that absorbs and re-emits part of what the gold
   transmits, a computed gap-closure time (plasma filling the gap), and the windows'
   l = 4, which neither quad powers nor pointing removed.
4. **Confirm against facility practice:**
   - the TaLIS 0.5° tolerance (pivot, and whether it replaces or adds to the error set);
   - the 1ω footprint side and shape;
   - which quads count as open for the 20° rule;
   - the VisRad items listed in the research notes (summary-CSV columns after Power,
     spot-rotation sense, SST preamble, element-export layout; in workspaces, the
     Rotation/Spin angle mapping, the cylindrical position order, primitive grid
     spacing and the spot-size definition of custom profiles: one small workspace
     with a rotated cylinder and a rotated CAD part would settle most of them);
   - the NIF CMT XML schema.
5. **Benchmark against VisRad.** Run `write_visrad_kit` cases at work and bring back the
   outputs for `compare_visrad`. This validates absolute drive levels and mode
   amplitudes.
6. **VisRad workspace import: what `read_visrad_workspace` skips.** Each is warned
   about or refused today. Most need a sample workspace that uses the feature, since
   the layout of the extra data is not seen yet.
   - Unconverted light: beams with 2ω/1ω on (the green/red power columns) as
     `nif_unconverted` colors or extra beams.
   - Time-dependent beam powers (power mode ≠ 0, `laser power table values` with
     N > 0; VisRad's time unit is seconds), energy-normalized beams, and named pulse
     shapes, as `PulseShape`s.
   - DPR beams (split into two cones, as VisRad does for OMEGA).
   - Pointing method 1 (distance toward a diagnostic port, from `[Diagnostic Port
     Data]`), and target-mount or diagnostic-port reference frames for objects and beams.
   - Objects positioned by a key point other than the centroid (cylinder top/bottom,
     box and rectangle corners).
   - Time-dependent albedo, x-ray conversion, lengths and power sources (the ten
     `Time-Dependent …` tables), laser reflectivity, and power self-emission.
   - Wall thickness (inner and side walls), partial spheres, and the other primitive
     types: cylinder, cone, disk, box, polygon, halfraum, rugby and slotted hohlraums,
     surfaces of revolution, multi-element objects.
   - CAD after `remesh_mesh` (see Done): view-factor receivers could use an element's
     outline instead of its ~10–14 CAD triangles (far-pair cost scales with them); an
     adaptive `angle` for thin rods and fillets (many narrow elements now); closing
     cracks where two CAD faces tessellate one curved edge differently (moves geometry
     slightly).
   - CAD capsule tents and LEH windows are opaque closed solids: a thin-film option
     that transmits (or automatic thin-solid detection in the library; the GUI already
     opens a sealed pocket) would let the cavity include the tent region.
   - CAD fits (see Done): tori for fillets (fillet pieces now fit as small tilted cones
     and spheres; `minfraction` hides them); one surface shared across parts (the
     halves' inner walls); a quadrature at a cone apex; TALIS curvature kind for cones;
     fitting by default (Open question 1; `minfraction = 1e-3` is the default since
     2026-10-24).
   - Adjacent flat CAD pairs: near-kernel evaluations scale as shared-edge length over
     r/256 (~700 per pair on a sliver tube even with straight cuts); an exact or coarser
     treatment along straight shared edges would cut them. Banding misses most of the
     NIF halves' long facets (long edges ~90° apart, fans, cones): circular plane
     families would extend it. `_earclip` falls back to raw triangles on pinched
     outlines.
   - Laser deposition on sliver CAD is slow in BVH `any_hit` (~5 s per beam on the full
     157k-element NIF CAD target; moot on the cavity).
7. **GUI on CAD: what remains.** Partial curved fits of CAD parts (needs a
   part-to-component name mapping, or a per-element surface on `Component`; on the NIF
   workspace 1.5× faster view factors and 0.6% on the halves' flux, the plan in Open
   question 1); a 0.75 s first-import tick and a 0.16 s first field-coloring
   frame (compilation); `include_components!` costs 0.3 s on the UI thread; a "fit the
   hohlraum template to CAD parts" action; the foil objective with one shared parameter
   set for beams aimed at one point does not move (per-beam parameters).
   `fit_surfaces` repeats rejected growths: a finely faceted face of none of the
   requested kinds (a tapered tube without cones, a fine fillet or torus) costs
   quadratic time, every patch seeding a growth the rms test rejects (20 s on a NIF
   capsule's fill tube with cylinders and spheres only). Skipping seeds inside rejected
   regions is ~5× faster but lost a true small cone on a NIF half; an rms test during
   growth, or a torus kind, might be safe.
8. **GUI requests:** frames on a port axis (the chamber,
   the entrance holes, the foil and other parts are done); clipping volumes; STEP import (needs a CAD kernel; STL
   and VisRad workspaces cover CAD meanwhile); an x-ray image that follows the camera
   as it moves; time in the window (pulse tables, drive histories and the retune's
   power histories plotted, rather than written to CSV).
   Adaptive mesh refinement of the radiation mesh where deposition varies fastest
   (the uniform fine mesh is converged but costs ~16k elements).
Items 9–13 and 16–17 were rewritten after the 2026-10-22 literature review against
Kritcher et al. 2024 (research notes, "Literature review against Kritcher et al. 2024";
full reports in `reference/litreview/`). Retuning to a reference drive (see Done) assumes
the reference's view factors hold through the pulse; items 9 and 10 remove that assumption.

9. **Time-dependent view factor: the P2 swing.** Done (2026-10-23, see Done): the moving
   wall on keyframes (`WallMotion`, `Keyframes`; error ∝ spacing²), closing LEH lips,
   and `BubbleAbsorption` with stagnation and ridge stop (`lpi_deposit`), compared with
   PRE 109 Fig. 10 for N221204 (item 16). Calibrated against published data
   (2026-10-24, see Done): wall speed, start and saturation (Young 2020, PRE 109
   Fig. 1a), bubble speed, launch and opacity (Izumi 2021), collision time (Izumi 2018),
   lips (Young 2020's LEH history). What remains:
   - the swing still turns positive ~0.4 ns early (6.26 against 6.67 ns) and its
     minimum comes ~0.3 ns early: inner absorption rises from 0.22 to 0.57 over
     5.5–6.5 ns on HYBRID-E. The HYBRID-E bubble speed is the weakest link (`energy` 1.1
     from the picket energy and quad splitting; PRE 109 Fig. 1a hides N221204's tip under
     its rays, ≥ 1.55 mm at 6 ns against the model's 1.96); a HYBRID-E bubble or thin-wall
     measurement would settle it;
   - the lip closure history (only Young's end value 0.857 r₀ is measured; the ramp is
     assumed); a computed lip closure (plasma filling the hole);
   - the capsule's shrinking radius (Young Fig. 3a has one: 0.83 mm at 5.5 ns, 0.22 at
     7.75 ns; PRE 109 Fig. 1a ~0.92 mm at 6 ns), Izumi 2023's low-contrast stall and
     relaunch;
   - Callahan's and Izumi's figures of merit beside the computed swing;
   - validation (2026-10-24, research notes "Young 2020 validation"): Young's
     static-vs-dynamic test run on `Keyframes` with his digitized radii and albedos: the
     85 % radii reproduce the moving geometry to 2 % in velocity, but his dynamic wall
     albedo (0.85–0.88 at peak against 0.8) lowers mode 1 by ~20 % where he reports
     4.6 % (open: needs N170827's delivered per-beam powers, unpublished; a proxy
     imbalance was used). Izumi 2018's SymCap changes, through Callahan's 2.7 µm per
     1 % CF as a transfer: signs all right (1.5× picket −30 µm, 1.18× capsule −12,
     0.15 mg/cc +24 against −21/−24, −21/−24, +10/+15), not blind (the bubble was fitted
     to Izumi). Callahan's 39-shot curve cannot be recomputed (per-shot inputs not
     tabulated); its line is digitized (101.9 − 85.0 x µm). Still skipped: Ralph 2018's
     18.5 % cone fraction for zero flux P2 (needs the high-foot pointing and measured
     wall radius, unpublished; his 51–62 % peak-averaged absorption is matched: 0.575
     on N221204). Fig. 10's N220919 is close
     on the proxy (minimum −0.140 against −0.145); N210808's proxy is not (foot +0.03
     against −0.05, minimum −0.09 against −0.14): its pulse and foot cone fraction are
     only approximated.
10. **Shot-data calibration.** Fit effective parameters (cone multipliers, CBET
   multiplier/clamp, bubble transparency, wall-motion speed, albedo) to measured in-flight
   and hot-spot P2/P4 and keyhole shock symmetry, with uncertainty. Reproduce the
   P2-vs-Δλ sequence N210808 → N220919 → N221204 → N230729 (1.8 / 2.5 / 2.75 / 2.9 Å at
   1ω: hot-spot P2 −2.3 / −10.8 / +0.3 / +6 µm) and Ruof 2026's rise-cone-fraction
   response (0.21 → 0.24 raises shell P2 9.1 → 27.5 µm; 0.21 → 0.16 leaves it). HYDRA's
   own knobs for comparison: clamp 0.002 (foot) / 0.0075 (peak), 2DConA fits need
   0.005–0.007; laser degraded 9–13 % or peak multipliers 0.825–0.86.
   CBET (`reference/litreview/cbet_calibration.md`). Done (2026-10-23, see Done):
   `TwoZoneCBET` (Δλ = 0 is not zero transfer) and a multi-epoch `calibrate_cbet`, fitted
   to N221204's post-transfer cone fraction against time in PRE 109 Fig. 9 (rms 0.016,
   dCF/dΔλ 0.16 per Å, within the 0.15–0.18 target). What remains:
   - a joint fit to N210808 and N220919 (multi-shot), which would also test the model's
     Δλ sensitivity (~2.5× HYDRA's cross-shot difference);
   - the fitted Te(t) reaches 7.4 keV late in the peak, above every published LEH value.
     Investigated (2026-10-24, see Done): Fig. 9b fixes only Te's rise with time; Te ≤
     5.5 keV fits as well with smaller multipliers, but moves the Δλ null from 1.2 to
     0.5 Å while bringing the 1.8 Å response closer to HYDRA's N210808. The joint fit
     with the null as a constraint decides it;
   - a clamp schedule driven by the `PulseShape` (the fit uses the deck's fixed
     schedule), a hard-minimum saturation option so HYDRA clamps carry over, plasma-state
     tables as inputs, outer-cone depletion and an inner-quad SBS proxy (the N220919
     limit and the 2.9 Å overshoot), and the 23.5°/30° split (the model favours 30°,
     HYDRA splits evenly).
   Item 9's dynamics: calibrated against published data (2026-10-24, see Done). Left:
   the HYBRID-E bubble speed (the picket-fluence transfer from the HDC basecamp is the
   least constrained input) and the lip closure history.
11. **Implosion metrics and UQ.** Done (2026-10-23, see Done; research notes,
   "Implosion metrics"). Open:
   - ~~two comparisons with Young 2020 do not match~~ **Done (2026-10-24)**: a setup
     difference (the tutorial's hohlraum). On Young's static N170827 setup with his
     digitized pulses and Green's function (`young2020_kernel`), dropped-quad
     velocities 336/349/161/152 against 317/337/170/146 km/s (ratio 2.2 against 2.07),
     offsets 0.09/1.30 (target) and 0.10/1.07 (capsule) against 0.1/1.3 and 0.1/1.1
     km/s/µm. One genuine error fixed: modes of an offset capsule were taken about TCC
     (`mode_center`). Test suite;
   - ~~Nora's high-yield branch and shape term~~ **Done (2026-10-24)**: Eq. 4 prints
     rounded coefficients (Fig. 6's are used; 40–80 refused), Table IV's burn column is
     reproduced; Eq. 2 and Table I in `shape_amplification`;
   - time-dependent error budgets and multi-shot trending.
12. **Backlighter physics.** CE(Z, I, τ) models with uncertainty, and the photometric
   chain to the detector.
13. ~~**Spectral radiosity.**~~ **Done (2026-10-24, see Done)**: multigroup radiosity,
    the laser-produced Au M-band, a synthetic Dante. Still open: measured Dante channel
    responses (filters, mirrors, XRD cathodes) in place of Gaussians and an unfold;
    XRD lines of sight; angle dependence of the M-band conversion (sin⁴θ, Dewald 2020) and
    its emission from the outer bubble (inner-beam clipping); emission of the closing LEH
    plasma (Dante's initial- vs open-area brackets differ by 20 eV at N221204's peak); a
    time-resolved M-band P2 history against Dewald's Fig. 10 on a moving hohlraum; the
    M-band split in the `CapsuleDrive` handover (item 17).
14. **Scale.** Hierarchical or iterative radiosity, Float32 storage for more than 10⁴
    elements, and GPU kernels.
15. **Polar direct drive mode**, benchmarked against published GA optimizers.
16. **Absolute energetics and an N221204 benchmark** (`reference/litreview/energetics.md`,
    `n221204_benchmark.md`). Done (2026-10-23, see Done): the mock HYBRID-E pointing
    (`hybride_mock_pointing()`) and `examples/n221204.jl` against PRE 109 Fig. 10. The
    P2/P0 swing has the right shape. With the dynamics calibrated to published data
    (2026-10-24, see Done) its minimum is −0.168 (Fig. 10 −0.173) and it ends +0.03 high
    (+0.123 at 8 ns against +0.095), but it still turns positive ~0.4 ns early (item 9).
    At 4.5 ns it sits 0.06 above Fig. 10 (DN24's input cone fraction there is below
    Fig. 9b's); the foot sits ~0.03 below (mock pointing, foot cone fraction read by
    eye). The energetics models are in (2026-10-24, see Done: `du_wall_loss`,
    `Conversion`/`nif_conversion`, `CapsuleAlbedo`, `energy_budget`, backscatter and
    degradation in the example). With all of them on, the split is wall 50.8 / LEH 33.0 /
    capsule 16.2 % (HYDRA 58.6 with gas / 25.2 / 16.2, the capsule calibrated) and the
    capsule Tr peaks at 334 eV (HYDRA ~300 internal, Dante 313): the published pieces
    together leave the wall too reflective (x-ray wall loss ~1.6× short at HYDRA's
    temperature). Still to do, each with its data:
    - the missing sinks and wall loss behind that excess: He gas (4.5 %), LEH windows, a
      wall loss above Lindl's gold at late times (the DU factor calibrated to Döppner's
      6.5 % came out 0.785, below Rosen's cocktail 0.83, a sign of the same thing), and
      the M-band (item 13); then recalibrate the capsule albedo (16.2 % needs α_c ≈ 0 at
      a drive that matches the LEH share);
    - ~~a synthetic Dante (item 13) reported beside the capsule-seen temperature~~
      **Done (2026-10-24)**: Dante 1 286 / 306 eV (initial / open LEH area) at 6.25 ns
      against 313 measured at 8.4 ns, capsule 302 eV at 5.5 ns; the model's drive peaks
      ~2.5 ns early and falls ~30 eV by 8.5 ns (default energetics);
    - the Zylstra 2022 LEH pair: the model needs ×1.10–1.15 more power for the 3.64 mm
      LEH (lips reduce it), the data ×1.19 or more; the shots' different histories and
      Dante's view are the likely gap, so a time-dependent run of both pulses is next.
    Benchmark targets still to check: bang time 9.53 ns as the end of the drive window,
    and hot-spot P2 ≈ 0 and velocity 36 km/s at (151°, 54°) through item 11. The pointing
    stays the public mock-up; the non-public export's aim points never enter the
    repository. CPPs, window geometry, corners and the liner thickness are unpublished
    and stay at repository defaults.
17. **Briareus → CompactLES** (`reference/litreview/hydro_coupling.md`). LLNL's two-step
    method with Briareus as the integrated step and CompactLES as the capsule-only step.
    Coupling (user, 2026-10-22): both codes are Julia, so couple in memory through a
    package extension (`BriareusCompactLESExt`, or the reverse; pick the side with the
    lighter dependency) rather than a static export. HDF5 output is optional, as an
    extension, for provenance and non-Julia consumers. The drive to hand over:
    T_r(t), capsule-averaged incident flux, a_l/a_0 for l ≤ 8 (P2 > 0
    pole-hot) on a ≤ 50 ps grid, the radius the flux was evaluated at (later a table over
    radius, and a Planck + M-band split from item 13), with provenance. Stages on the
    CompactLES side: (0) an absorbed-flux energy source in a thin layer at the ablation
    front with carbon exponents P_abl ∝ T^3.3 t^−0.09 and the prefactor calibrated to
    HYBRID-E (first shock > 12 Mbar at ~150 eV, ~140–150 Mbar at peak), checking v_imp
    380–400 km/s, 5–6.6 % mass remaining, 19–33 kJ coupled, Hurricane's free-fall
    condition and Masse's P2 sign flip near 4.5 ns; (1) implicit gray flux-limited
    diffusion with a Marshak boundary carrying the angle-resolved flux, flux-limited
    conduction at 0.15, QEOS/FEOS tables and TOPS opacities, validated first on the
    Shestakov 1998 Be/D capsule at 160 eV, then Kritcher 2024 Table III; (2) multigroup
    with the M-band. Bang time (9.53 ns) is an input in LLNL's capsule-only runs and the
    first prediction here. Pitfalls: apply the asymmetry near the ablation front, not at
    a distant fixed radius; the M-band has its own larger P2.
    *Stage 0 started (2026-10-24):* Briareus `CapsuleDrive`/`capsule_drive` (785684d);
    the extension lives in CompactLES (`CompactLESBriareusExt`, branch `briareus-drive`,
    unmerged): the code is CompactLES physics and Briareus stays free of MPI/KA. The
    1-D N221204-class run (static Briareus drive, 297 eV peak) holds the shell pressure
    on `P_abl` (17 Mbar foot, 145 peak) but reaches 315 km/s with 19.6 % ablator left:
    Lindl's ṁ under-ablates. Next: calibrate ṁ to 5–6.6 % remaining, then v_imp; the
    calibrated-dynamics/CBET drive; r–θ with the a_l (Masse sign flip); r–z with AMR.
    Research notes, "Capsule drive hand-over and CompactLES stage 0".
