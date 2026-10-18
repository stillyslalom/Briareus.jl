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

## Next
1. **VisRad deposition parity: needs runs in VisRad.** The goal: optimize in
   Briareus with `profile = visrad_profile(ws)`, export with `write_visrad_workspace`,
   and get the same deposition in VisRad. Near focus that holds today. What VisRad
   does away from focus, and at envelopes, shadows and tilted surfaces, needs VisRad
   output, all on the foil workspace unless noted:
   - **Intersection-point export** (File ▸ Export, laser beam–target intersection
     points; VR/files/laser_target_intpts_file.html). It lists a(z), b(z), n(z),
     f_SQ(z), ΔL_W(z), ΔL_H(z) and M(z) at each hit, which pins down the parts of the
     defocus law the manual leaves out: z₀ in n(z), and Γ, ΔL(z), M(z) in the
     square-beam blend. Then add them to `VisRadProfile` and recheck the foil wings
     (now up to ~20% off at 0.2–0.6 mm from focus).
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
   - **Power on tilted surfaces**: VisRad's foil total was exact to 1e-8, while
     per-plane normalization gives 1.1e-4 too much at ~40°. A steeper foil (~70°)
     would show whether VisRad renormalizes each beam's deposited total.
   - **Shadowing**: a washer or sphere in front of the foil, and a hohlraum case with
     LEH clipping. VisRad decides occultation on its 25 × 24 ray grid; Briareus uses
     the beam axis through each point.
   - **End-to-end**: a primitive hohlraum with a few quads, optimized here, exported,
     and run in VisRad, with the workspace saved after computing. Reading the
     `[Power Source Data]` tables back (a reader still to write; the rectangle's
     elements run local-x fastest, the hohlraum's order is unknown) would turn this
     into a routine comparison.
2. **Sombrero and backlighter geometry.** Model the LEH shield as a component (a
   dimpled annulus), and the backlighter foil's placement relative to the diagnostic
   line of sight. The line-of-sight and clearance checks in the radiography tutorial
   then decide which bundles can drive the foil. This needs real dimensions.
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
     default filters (see the decisions in the overnight log).
   - `orient_mesh` on hollow solids (a shell with its inner surface modeled) faces the
     inner surface into the material; containment parity would face it into the gas,
     but `enclosure` from TCC would then select a hollow capsule's interior (the GUI's
     cavity from TCC relies on the current behavior). A user decision.
   - Adjacent flat CAD pairs: near-kernel evaluations scale as shared-edge length over
     r/256 (~700 per pair on a sliver tube even with straight cuts); an exact or coarser
     treatment along straight shared edges would cut them. Banding misses most of the
     NIF halves' long facets (long edges ~90° apart, fans, cones): circular plane
     families would extend it. `_earclip` falls back to raw triangles on pinched
     outlines.
   - Laser deposition on sliver CAD is slow in BVH `any_hit` (~5 s per beam on the full
     157k-element NIF CAD target; moot on the cavity).
7. **GUI on CAD: what remains.** Partial curved fits of CAD parts (needs a
   part-to-component name mapping, or a per-element surface on `Component`; neutral on
   the NIF workspace); a 0.75 s first-import tick and a 0.16 s first field-coloring
   frame (compilation); `include_components!` costs 0.3 s on the UI thread; a "fit the
   hohlraum template to CAD parts" action; the foil objective with one shared parameter
   set for beams aimed at one point does not move (per-beam parameters).
   `fit_surfaces` repeats rejected growths: a finely faceted face of none of the
   requested kinds (a tapered tube without cones, a fine fillet or torus) costs
   quadratic time, every patch seeding a growth the rms test rejects (20 s on a NIF
   capsule's fill tube with cylinders and spheres only). Skipping seeds inside rejected
   regions is ~5× faster but lost a true small cone on a NIF half; an rms test during
   growth, or a torus kind, might be safe.
8. **GUI requests:** a beam-model choice (Briareus's f-cone or `VisRadProfile`), set
   to VisRad's with the workspace's parameters when one is imported, and passed to the
   drive's `Scenario` (`controllers/drive.jl`); frames on a port axis (the chamber,
   the entrance holes, the foil and other parts are done); clipping volumes; STEP import (needs a CAD kernel; STL
   and VisRad workspaces cover CAD meanwhile); an x-ray image that follows the camera
   as it moves; time in the window (pulse tables, drive histories and the retune's
   power histories plotted, rather than written to CSV).
   Adaptive mesh refinement of the radiation mesh where deposition varies fastest
   (the uniform fine mesh is converged but costs ~16k elements).
9. **Time-dependent view factor.** Wall, LEH and spot motion from rad-hydro or a
   parametric model; P2(t) swing. Retuning to a reference drive (see Done) assumes
   the reference's view factors hold through the pulse.
10. **Shot-data calibration.** Fit effective parameters (cone multipliers, CBET
   multiplier/clamp, albedo) to measured P2/P4 and shock symmetry, with uncertainty.
   Reproduce P2-vs-Δλ curves.
11. **Robust optimization and UQ: what remains.** Correlated errors (per bundle, per
   cone; target positioning, which moves every beam at once), time-dependent budgets
   (errors through `drive_history`), and inverse solves for a measured mode (MacGowan
   2021) as a documented workflow.
12. **Backlighter physics.** CE(Z, I, τ) models with uncertainty, and the photometric
   chain to the detector.
13. **Spectral radiosity.** Multigroup emission, detector responses along Dante and XRD
    lines of sight, and brightness temperatures.
14. **Scale.** Multi-ray partial visibility for view factors, hierarchical or iterative
    radiosity, Float32 storage for more than 10⁴ elements, and GPU kernels.
15. **Polar direct drive mode**, benchmarked against published GA optimizers.
