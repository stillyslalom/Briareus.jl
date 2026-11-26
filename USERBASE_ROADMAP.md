# Briareus.jl Userbase Wedge Roadmap

This file is the working copy of the usability and adoption plan handed over on
2026-11-25, kept on the branch so that a session can resume from it. The status section
comes first and the plan follows as handed over. The cleanup plan that preceded it is
archived at [reference/archive/cleanup-roadmap.md](reference/archive/cleanup-roadmap.md).

## Status

Marks: **done**, **partial** (with what remains), **open**.

| Item | Status | Note |
|---|---|---|
| A.1 witness geometry, shared 1ω contour | done | `omega1_contour`, `Omega1Approach`, `omega1_closest_approach`; the application's `unconverted_footprint` takes the mode, the view toggle for it belongs to B.2 |
| A.2 unsigned curvature, specular accuracy | partial | cones and tori, the normal spread, the reflection record with port and areas, the mesh estimate for free-form CAD, and the reflected-flux estimate (`reflected_fraction` of the record, `omega1_reflected_power` with a reflectivity of 1 by default, an upper bound without shadowing) are in. Left: the estimate runs on every facet of an unknown component instead of the exposed ones only (0.25 s per 80,000 facets) |
| A.3 hardware disposition field | done | `omega1_disposition`, per mesh triangle, with the shadow flag; 85,000 triangles under 48 quads take about 4 s at `level = 0` and 26 s at `level = 2` on four threads |
| A.4 VisRad 1ω/2ω beam import | done | `visrad_unconverted`; the import sets the participation flags and the drive ratios. Only the fraction-of-total model is measured on a workspace |
| A.5 STEP B-rep import | open | optional |
| B.1 disposition coloring | partial | the two color modes, the tight threshold and the curved-safe filter are in the 3D view, and a click on a facet gives its verdict in the status line; there is no hover |
| B.2 closest-approach faces, calipers, cones | done | a click on a TaLIS row draws the nearest side, the clearance plane, the caliper and the reflection cone, with a ring about the threatened port; the margin terms and the port angles show on hover. The Pinch toggle adds every other beam within the tight threshold of the part. View ▸ Beams switches the drawn 1ω light between the nominal footprint, the 10% contour and the worst-case contour. The application has no chamber indicator, so the port is marked in the 3D view |
| B.3 daylight and grazing views, slice | done | the Daylight and Grazing buttons turn the camera to the row under inspection, and the TaLIS tab draws the 2-D slice of a selected 1ω row (`omega1_slice`; `write_omega1_slice_svg` for a file) |
| B.4 safe-corridor map | partial | `corridor_map` in core, on a sphere or a plane, over 3ω envelopes, 1ω footprints and lines of sight (`SightLine`, `sight_line`), with `corridor_margin` for one point and `write_corridor_map_svg`; 192 beams on 90 by 180 samples take about 0.2 s. The clearance page draws the sphere about TCC, or a plane through a point with a normal and a half-width, in the 3D view, and the Clearance tab draws either flat (polar disks or a square) with an Export SVG button. The status line names a ticked diagnostic that gives no line of sight. The 3D view draws only the side of the sphere that faces the eye (all of it from inside), since the far side drawn through the faint open cells made the lower half look open. Left: hardware shadowing |
| B.5 VisRad parity and design diff | partial | `visrad_residuals` gives the element residuals of the laser deposition and of VisRad's element tables, with `resample = true` for a target meshed differently from VisRad, and the coloring Laser, Briareus − VisRad paints the laser residual in the window. The Drive tab shows the laser rows of `compare_visrad` per object (absorbed powers, rms, largest difference and correlation of the flux, from the `objects` of the resampled residuals). `design_diff` compares parts, beams, TaLIS rows and the physics of each role, with `write_design_diff` and `briareus diff`. Left: the residuals of the x-ray flux in the window (it has no reader for VisRad's element tables), and the diff as an overlay in the 3D view |
| B.6 look from a diagnostic, CAD prep | partial | `sight_intrusions` gives the margin of each part to the volume of a line of sight. View ▸ Look from diagnostic looks along the line of sight of any definition of the facility that gives one, selects the hardware within the target positioning error of the view and names it and the observed parts in the status line. Left: the synthetic image laid over the CAD, the frame of the collimator or crystal, and the guided three-step CAD preparation |
| C.1 outage engine | done | `roe_scan` over cases (a scenario and beams per epoch or observed face) and `ROECriterion` values: single beams, single quads, pairs, companion drops, roles by absorbing part. Linear tables without a laser–plasma model, full evaluation with one (`method`) |
| C.2 criteria and tiers | partial | built-in criteria for capsule drive, modes, hot-spot velocity and yield, face flux and nonuniformity, temperature, flux balance, foil irradiance, nonuniformity, centroid and He-α signal, hard-sphere and corona nonuniformity; tiers with the binding criterion, pair escalation and a per-role budget. `roe_default_criteria` selects a starting set from the observed part and the foils of a design. Left: a crystal (Bragg) throughput for CBI, the ablation-pressure proxy as a built-in, triple outages. Every built-in scores the change that an outage makes to the full shot, and the velocity thresholds (30 and 65 km/s) are the speeds of a 2% and a 10% loss of the piston-model yield. Velocities past 0.75% of mode 1 (75 km/s) are linear upper bounds; the tiers hold since both thresholds lie below it (research notes) |
| C.3 port map and export | done | `roe_table`, `write_roe_csv`, `roe_markdown`, the chamber port map as SVG for any facility (`port_map`, `port_map_svg`), the HTML sheet (`roe_html`, `write_roe_packet`) and the `briareus roe` command. In the window, Scan outages colors the beams by tier and a selected beam shows its binding criterion; the Outages tab draws the port map by tier, and marks clicked on it form a hypothetical outage of any number of units, scored by `roe_outage` against every criterion. A NIF quad held as one beam stands for its four beams (`stands_for` of `roe_scan`, `roe_stands_for`): its quad unit and a beam unit of a quarter of its power, so the sheet of the default preset shows 48 red quads and 48 amber single beams where it showed 48 red "beams", and each mark of the map shows the tier of one beam in its center. The marks carry no role badges |
| D.3 handoff loop | done | the SVG slice of a 1ω approach is written from core without a graphics context, `briareus design.briar --tab <name>` opens the window on a results tab, and File ▸ Export Julia script (`write_session_script`) writes the design as a script for `briareus run`. The script rebuilds the target through the controllers from the session state, and CAD files are named by path |
| D.1 command line | partial | `briareus talis`, `review-1w`, `roe`, `diff`, `optimize` and `run script.jl` on saved sessions or VisRad workspaces, with `--out` and `--json`, each one core call with the settings of the design; `optimize` runs the saved optimization and writes the optimized session, SST pointing and, for a workspace, a copy of it. The commands load the whole application package first. Left: the startup time of a packaged app |
| D.2 HTML packets | partial | `write_omega1_review`, `write_roe_packet` and `write_design_diff` write one self-contained page each (inline SVG slices and port map, sortable tables, nothing fetched, the same text for the same input), with CSV, SVG and JSON beside it. `briareus review-1w --images` adds a rendered Daylight view of each pair, drawn offscreen. File ▸ Export 1ω review report writes the packet from the window without the pictures, since the offscreen drawing and the Qt window cannot share a process. Left: pictures from the window, and per-role summaries beyond the tables |
| D.4 LLM kit | partial | `docs/src/llms.md`, `llms.txt` and the Claude skill `skills/briareus/SKILL.md` (both written by `julia dev/llms.jl`), guarded by `test/docrefs.jl` against drift and names that stop existing. Left: the optional MCP server |
| E.1 prebuilt Windows and macOS releases | partial | `dev/release/build.jl` (PackageCompiler `create_app` with `julia_main`, a precompile run of every headless command and the window), `smoke.sh` (the launcher with no Julia on PATH), `sign_macos.sh` with its entitlements, and `.github/workflows/Release.yml` for Windows x64, macOS arm64 and Linux x64, with signing gated on secrets. Left: a build that finishes (the system image ran out of memory at about 14 GB on a 15 GB Linux machine in all three configurations tried), the smoke test and launch time on built applications, an installer, a `Briareus.app` bundle, and the signing secrets (dev/release/README.md) |
| E.2 web and Tauri frontend | open | not started; a frontend rewrite over `Controllers`, planned in months in the plan |

## The plan as handed over

A staged execution plan for turning Briareus.jl into an adoption-ready tool for VisRad/HED target designers and LLM-assisted shot-review workflows.

**Core thesis:**
1. **Visual trust is the primary adoption wedge.** Physicists trust a 3-D geometry/clearance check when they can see the exact margin with their own eyes—especially on an existing VisRad `.vrw` workspace without rebuilding the target.
2. **1ω clearance + curvature/specular visualization is the killer first feature.** Replacing "inflate 48–192 1ω volumes and manually orbit for daylight" with hardware disposition maps, closest-approach beam faces, exact calipers, and one-click orthographic Daylight Views solves an immediate, high-friction NIF shot-review pain point.
3. **LLM scriptability & shot-review automation are the compounding moat.** Unlike VisRad, Briareus can round-trip cleanly between human visual inspection (`.briar` / `.vrw`) and LLM-driven sweeps, optimizations, and automated review packets—even on machines with no separate Julia installation.
4. **Distribution must meet users where they are.** Eliminate the `juliaup` + 10-minute compile + `PATH` edit + 20-second startup barrier first via `PackageCompiler.jl` bundles, then eliminate Qt/GLMakie packaging fragility permanently via a Web/Tauri frontend over the existing headless `Controllers` submodule.

---

## Ground Rules & Post-Cleanup Baseline

- **Execution timing:** This roadmap begins **after all phases (Phases 1–6) of the cleanup roadmap are complete and merged.** Assume the cleaned codebase baseline:
  - Unified beam-envelope and 1ω setup helpers (`src/clearance.jl`, `src/talis.jl`, `src/design_rules.jl`, `src/unconverted.jl`), no duplicate/shadowed utilities or data races.
  - Trimmed top-level `Briareus` exports (paper/shot fixtures moved to `Briareus.Benchmarks`) and a formalized `Briareus` ↔ `BriareusGUI` interface (never add new `Briareus._*` reaches from `BriareusGUI`).
  - Split `BriareusGUI/src/controllers/` and `BriareusGUI/src/qt/` modules (`platform.jl`, `models.jl`, `refresh.jl`, `callbacks.jl`, `window.jl`) with granular dirty flags and cleaned session schema.
- **Design for the full menagerie of NIF/HED experiments (never hohlraum-only):** Every core engine, GUI visualization, diff, ROE table, and review packet must compose cleanly across compound targets with multiple beam roles, multiple observed surfaces/cavities/foils, and distinct timing epochs—covering indirect-drive hohlraums/frustraums, single and dual halfraums (shock tubes, opacity/EOS packages), single and double backlighter foils, Crystal Backlighter Imager (CBI) / narrow-band crystal radiography setups, polar/symmetric direct drive, and Discovery Science platforms. Tie into the multi-part/role model in `reference/assemblies-design.md` (multiple observed surfaces, multiple foils, per-cavity `Tr`, per-role beam groups).
- Follow `CLAUDE.md` rules throughout (`test/docrefs.jl`, `test/inference.jl`, `test/aqua.jl`, `dev/prose/lint.jl`, 0.x no-compat-shims policy).

---

## Workstream A: 1ω Clearance, Unsigned Curvature & VisRad Parity (Core `src/`)

Build the exact geometric and rule engine that powers both the 3-D inspection UI and automated shot-review packets.

### A.1 Witness geometry & unified 1ω contours (`src/talis.jl`, `src/unconverted.jl`)
- **Problem:** `talis_check` (`src/talis.jl:551-576`) returns only `row.point` (the worst-case hardware point) and discards which 1ω beamlet won, the axial station $s^*$, the closest point on the 1ω square boundary $x_{\text{beam}}$, which face(s) of the square pyramid were closest, and which focus extreme ($0, \pm\text{focus}$) drove the worst case. Meanwhile, `BriareusGUI/src/controllers/beam_geometry.jl:132-135` draws the nominal 1ω square (`edge = 4.8 mm`, `size = 22.0 mm`) instead of `talis_check`'s 10% contour (`4.45 mm` / `22.7 mm`) or worst-case inflated contour.
- **Action:**
  1. Define a single shared 1ω contour helper (`omega1_contour(; mode = :talis_10pct)` supporting `:nominal`, `:talis_10pct`, and worst-case inflation) used by `src/talis.jl`, `src/unconverted.jl`, `src/design_rules.jl`, and `BriareusGUI`.
  2. Extend `_omega1_square` / `_omega1_distance` (`src/talis.jl:551-576`) to compute the exact closest point $x_{\text{beam}}$ on the 1ω square at station $s^*$ (clamping transverse $(p, q)$ to $[-h(s^*), +h(s^*)]^2$ in the $(\hat{o}, \hat{q})$ plane) and face bitmask (which of the 4 planar sides/edges owns the closest point).
  3. Add a lightweight struct `Omega1Approach` and query `omega1_closest_approach(target, beams, row_or_pair; ...)` returning:
     - `beam`, `beamlet` `(k, f3w, ô)`, `component`
     - `x_hw`, `x_beam`, `s_star`, `face_mask`, `focus_shift`
     - `nominal`, `worst`, `required`, `margin`
     - `budget` NamedTuple (`pointing`, `target_offset`, `hardware_tolerance`, `rotation_tol`, `focus_delta`, `required`) for transparent UI/report display.
  4. Verify via tests that `norm(x_hw - x_beam)` in the transverse plane matches `row.nominal` (at `focus_shift = 0`) and `row.worst + budget...` to `1e-12` mm.

### A.2 Unsigned local curvature & specular-reflection accuracy on real CAD (`src/talis.jl`, `src/cadfit.jl`)
- **Domain rule:** The NIF 1ω management concern is far-field specular reflection into another beamline at the ~7.7 m chamber/optics distance, not local millimeter-scale focusing. Both convex and concave smooth surfaces with unsigned $|R_{1D}| < 2\text{ mm}$ (cylinder/cone) or $|R_{2D}| < 5\text{ mm}$ (sphere/torus/compliant dimple with `half_side / |R| ≥ 0.25`) diverge rapidly into the far field and pass.
- **Action:**
  1. **Extend `_talis_curvature` (`src/talis.jl:531-547`) to `ConeSurface` and `TorusSurface`:**
     - `ConeSurface`: evaluate local principal radius $|R(x)| = r(x) / \cos\alpha$ per 1ω-exposed triangle against the 2 mm 1-D limit (since a tapered stalk can cross 2 mm along its length).
     - `TorusSurface`: evaluate both unsigned principal radii ($r_{\text{minor}}$ and $|R_{\text{major}}/\cos\psi + r_{\text{minor}}|$) on exposed triangles.
  2. **Add per-exposed-triangle local curvature fallback for unsplit/free-form CAD (`<name>/rest`):**
     - First filter triangles with the fast spatial test `_talis_reaches` (`src/talis.jl:509`), narrowing 100k+ CAD triangles to only the exposed subset.
     - On exposed triangles lacking a `CurvedSurface`, fit a local 2-ring vertex quadric / Weingarten shape operator to estimate unsigned principal radii $(|R_1|, |R_2|)$ and fit residual:
       - Low residual & $\max(|R_1|, |R_2|) < 5\text{ mm}$ (both principal curvatures strongly nonzero) → exempt (2-D scattering).
       - Low residual & $\min(|R_1|, |R_2|) < 2\text{ mm}$ → exempt (1-D scattering).
       - Low residual & gentle/flat ($|R| \ge$ threshold) → run specular reflection check.
       - High residual or overly coarse faceting → flag `:review` explicitly.
     - This avoids per-component `:unknown` failures on mixed CAD solids (e.g., a 1.5 mm stalk attached to a flat mounting flange in one STL) and avoids paying `fit_surfaces`'s quadratic cost on non-exposed mounts (`reference/ROADMAP.md:863-864`).
  3. **Account for facet-normal spread in the 20° specular check (`src/talis.jl:507-518`):**
     - On gentle curved surfaces ($|R_{1D}| \ge 2\text{ mm}$ or $|R_{2D}| \ge 5\text{ mm}$), add the continuous normal sweep across each exposed facet (from vertex normals / `CurvedSurface` normals) into the angular margin against open ports (`αw - 2 * tolerance - facet_spread - reflection_angle`), preventing coarse meshes from falsely clearing a port.
     - Also return the threatened open port ID/direction and reflected vector `r` on `_talis_management` rows (or via a companion query) so the GUI can draw the reflection cone.
  4. **Flux/area-aware reporting (no single-sliver false alarms):** Because the physical hazard is reflected 1ω flux entering an opposing beamline rather than isolated microscopic facets or small-area corners/dihedrals, report the total **exposed projected area** (and optional reflected-flux estimate when 1ω power is known) of each offending flat/gentle patch alongside its angular margin—so tiny CAD bevels/slivers don't drown out real specular reflectors.

### A.3 Hardware 1ω disposition field API (`src/talis.jl` or `src/unconverted.jl`)
- Add `omega1_disposition(target, beams; ...)` returning, per hardware element (or triangle):
  - `status::Symbol`:
    - `:clear` — worst-case margin `> tight_threshold` (default `1.0 mm`)
    - `:tight` — `0 < margin ≤ tight_threshold`
    - `:curved_safe` — inside worst-case 1ω volume, exempt by unsigned $|R_{1D}| < 2\text{ mm}$ or $|R_{2D}| < 5\text{ mm}$
    - `:specular_clear` — exposed flat/gentle face, specular bounce clears all open ports by `> 20°` (with tolerance)
    - `:review` — exposed, curvature/dimple-pitch uncertain or requires TaLIS review
    - `:fail` — exposed flat/gentle face reflecting within `20°` of an open port, or clearance fail on no-1ω hardware (e.g., diagnostics)
  - `margin`, `beam_idx`, `beamlet_idx`, `R1`, `R2`, `nearest_port`, `port_angle`, `shadowed::Bool` (optional fast BVH visibility check so reviewers can distinguish directly lit vs upstream-shadowed hardware).

### A.4 Complete VisRad `.vrw` 1ω/2ω beam import (`src/visrad_workspace.jl`, `src/visrad_json.jl`)
- **Problem:** `reference/ROADMAP.md` item 8 notes that VisRad import currently reads 3ω powers only and skips VisRad's green (2ω) and red (1ω) beam powers/flags.
- **Action:** Parse VisRad's 1ω/2ω beam power columns / flags on `.vrw` import (both text and JSON formats) into `VisRadWorkspace` and wire them into `BriareusGUI`'s `import_visrad_workspace!` so dropping in a real NIF `.vrw` activates the exact right 1ω beams immediately with zero manual setup.

### A.5 Direct STEP B-rep Import (Optional Extension / Sidecar)
- **Motivation:** STL strips away CAD's exact analytic B-rep geometry (cylinders, cones, spheres, tori, planes), forcing `src/cadfit.jl` to reverse-engineer them from triangles and forcing local quadric approximations for 1ω curvature.
- **Action:** Add an optional STEP importer (via `OpenCASCADE_jll` / `Gmsh_jll` or a lightweight sidecar converter) that reads target-fab `.step` files directly into `(Mesh, exact CurvedSurface per B-rep face)`—delivering exact, instant 1ω curvature classification and curved-surface view factors without slow `fit_surfaces` runs on 100k+-facet mounts (`reference/ROADMAP.md:864`).

---

## Workstream B: Visual Trust & Shot-Review GUI (`BriareusGUI/`)

Turn 1ω shot review, VisRad verification, CAD setup, and diagnostic checks into fast, unambiguous visual workflows.

### B.1 4-State Hardware 1ω Disposition Coloring
- **Files:** `BriareusGUI/src/controllers/view.jl`, `BriareusGUI/src/controllers/talis.jl`, `BriareusGUI/src/canvas/scene_canvas.jl`, `BriareusGUI/src/qml/panels/TalisPanel.qml`
- Add `:omega1_disposition` (categorical: Gray = clear, Teal/Blue = curved-safe, Amber = tight / specular-clear / review, Red = fail) and `:omega1_margin` (continuous distance heat-map near 0 mm) to `COLOR_MODES`.
- Hovering/picking any colored hardware face reports: local $(|R_1|, |R_2|)$, worst-case 1ω margin, owning beamlet name, exposed projected area, and (if exposed) nearest open port + reflection angle.
- Include a compact filter toggle: *"Hide curved-safe (`|R| < 2/5 mm`)"* and *"Margin threshold (mm)"* slider so reviewers focus strictly on tight/failing hardware.

### B.2 Closest-Approach Beam Surfaces, Calipers & Specular Cones
- **Files:** `BriareusGUI/src/controllers/beam_geometry.jl`, `BriareusGUI/src/canvas/scene_canvas.jl`
- Add `closest_approach_surface(approach::Omega1Approach; span = 6.0)` that builds:
  1. **Clipped closest planar face(s)** of *only* the winning 1ω beamlet (from `face_mask`, axially clipped to `[s_star - span, s_star + span]`)—drawing the **nominal 10% contour** as a crisp outline and the **worst-case inflated clearance plane** (errors + `r·tan(0.5°)` + `0.2 mm` at worst focus) as a translucent tinted sheet.
  2. **Multi-beam pinch mode:** Optionally show the closest-approach face of *every* beamlet with `margin < threshold` (e.g., `< 1.0 mm`) against the selected hardware part—clean surgical keep-out walls instead of full 3-D box soup.
  3. **3-D caliper segment `(x_hw, x_beam)`** with endpoint markers and an on-screen/inspector margin budget breakdown (`nominal`, `-pointing`, `-pos`, `-0.5° tol`, `-0.2 req = margin`).
  4. **Specular reflection cone (for `:fail` / `:specular_clear` rows):** Draw the incident 1ω ray into `x_hw`, the reflected specular axis `r`, and the `±(20° + 2×0.5° + facet_spread)` exclusion cone pointing toward the threatened open port (highlighted on the chamber sphere/direction indicator).
- Fix `unconverted_footprint(1)` in `BriareusGUI/src/controllers/beam_geometry.jl:132-135` to support `:nominal`, `:talis_10pct`, and `:talis_worst` modes so full-beam 1ω toggles agree visually with `talis_check`.

### B.3 One-Click Orthographic "Daylight View", "Grazing View" & 2-D Slice Inset
- **Files:** `BriareusGUI/src/canvas/cad_camera.jl`, `BriareusGUI/src/canvas/scene_canvas.jl`, `BriareusGUI/src/canvas/map_canvas.jl`, `BriareusGUI/src/qml/panels/TalisPanel.qml`
- Clicking any 1ω row (or hitting a **Daylight view** button) computes canonical camera poses from `(x_hw, x_beam, k)`:
  1. **Daylight View (down-the-beamlet orthographic):** Camera looking along $\pm k$, `ortho = true`, centered on `(x_hw + x_beam) / 2`, zoomed to `~4×` the local caliper/feature extent. Projects the 1ω square edge-on in the exact transverse plane where TaLIS defines clearance—visible daylight equals exact transverse margin.
  2. **Grazing View (edge-on to closest face):** Orthographic camera parallel to the closest 1ω face and perpendicular to `x_hw - x_beam`, showing the hardware profile grazing the keep-out wall axially.
- **Paired 2-D transverse cross-section at $s^*$:** Render a clean 2-D engineering slice in the beamlet's $(\hat{o}, \hat{q})$ plane at $s^*$: nominal 10% 1ω square, worst-case inflated square, hardware cross-section/projected silhouette near $s^*$, and caliper dimension arrow.

### B.4 Safe-Corridor / Keep-Out Volume Mapping (Inverse Clearance)
- **Goal:** Answer *"where am I allowed to place a new stalk, shield, foil, CBI crystal/snout, or diagnostic?"* rather than only checking existing hardware.
- Compute and display a **2-D polar safe-zone map** at a user-chosen radius $R$ from TCC (or on an arbitrary plane slice through the mounting zone) coloring **3ω beam keep-outs**, **1ω beam keep-outs**, and **all active diagnostic lines of sight** (Dante, SXI, CBI / backlit radiography LOS, VISAR/SOP), leaving **open green corridors** with exact angular/linear margins.

### B.5 Live VisRad Parity Check & Design Diff Overlay
- **VisRad Parity Check:** When loading a `.vrw` file that includes VisRad's saved deposition/radiosity (`src/visrad_results.jl`), offer a one-click **Compare with VisRad** view that paints element-by-element residuals (`Briareus − VisRad`) on the 3-D mesh alongside a role-aware summary table (capsule modes/drive, per-cavity halfraum `Tr`, package face flux/uniformity, foil irradiance, per-component powers)—building immediate user trust on their own shot file.
- **Design Diff (`briareus diff a.{vrw,briar} b.{vrw,briar}` + GUI overlay):** Load a reference design beside the current design to show moved CAD parts as translucent ghosts, beam aim-shift vectors, clearance/TaLIS margin deltas (newly failing vs. improved), and role-aware physics deltas across every active cavity, capsule, package face, and backlighter foil (`ΔTr`, `ΔP2/P0`, `Δmode-1`, `Δface_flux`, `Δbacklighter_signal`, `Δpicket_fluence`). Also makes LLM-generated modifications immediately auditable.

### B.6 "Look Through the Diagnostic" View & Guided CAD Prep Flow
- **Look from Diagnostic:** One-click dropdown covering every configured diagnostic in the shot (Dante 1/2, SXI, G-LEH, CBI / pinhole radiography LOS, VISAR, or any TOML definition in `data/diagnostics/`) that snaps the 3-D viewport to that diagnostic's exact port, FOV, and collimator/crystal/pinhole frame, highlights any stalk/shield/hardware intruding into the FOV (with target-positioning tolerance), and overlays the synthetic x-ray image, radiograph projection, or visible LEH/package area directly on top of the CAD.
- **Guided CAD-to-Target Prep:** Wrap `connected_pieces`, `orient_mesh`, `enclosure` (supporting multi-cavity detection such as dual halfraums), `remesh_mesh`, and `curved_components` into a 3-step import flow (Orient & split → Detect cavities, foils & external hardware → Assign multiple observed/foil/wall/hardware roles & preview curvature fits) so raw target-fab assemblies become analysis-ready in seconds.

---

## Workstream C: Pre-Shot Beam Rules of Engagement (ROE) Planner (Core + GUI + CLI)

Because shot-day operations are planned in advance (no day-of repointing or power tweaks), the operational need is pre-computing **which beams and quads are critical (delay shot), which require Shot RI approval to drop, and which can be dropped automatically**—generalized across the full menagerie of NIF experiments and compound multi-role shots.

### C.1 Composable, Role-Agnostic Outage Engine (`src/roe.jl`)
- Build ROE around a **composable `ROECriterion` protocol** evaluated over any combination of target parts, cavities, foils, diagnostic LOS, beam subsets, and timing epochs (never hardcoding a single hohlraum + capsule):
  - Leverage `Scenario`'s precomputed `ylm_response` / `moment_response`, `FaceFlux` responses, `influence_matrix` (`src/optimize.jl`), and per-foil deposition linearity for millisecond linear evaluation, paired with optional full LPI (`apply_lpi`) re-coupling verification.
  - **Single-beam drops** (all active beams) and **single-quad drops** (all active quads), grouped by beam role (drive vs. backlighter vs. secondary halfraum vs. CBI).
  - **Pairwise outage scan** (all quad pairs and configurable beam pairs, scoped within or across beam roles) to catch constructive interference where two individually safe beams become unsafe when dropped together (e.g., same-side mode-1 kick in a hohlraum, or losing 2 of 4 beams on a single foil/halfraum).
  - **Optional pre-planned symmetric companion-drop policies:** Evaluate whether pairing a dropped beam/quad with a pre-approved opposing drop restores symmetry/balance inside limits without any power or pointing tweaks.

### C.2 Built-In Criteria Across the NIF Experiment Menagerie + Custom Hooks
Score every beam, quad, and pair simultaneously against user-configured `(green, amber)` thresholds across every active role in the shot:
- **Indirect-drive hohlraums / frustraums:** Capsule `ΔTr` / `Δflux`, `|ΔP2/P0|`, `|ΔP4/P0|`, `hotspot_velocity` (`km/s`), `yield_factor` (`src/implosion.jl`), and per-phase (picket vs. peak) drop impact.
- **Halfraums, dual-halfraum shock tubes & opacity/EOS packages:** Per-cavity `wall_Tr` (`cavities`), `FaceFlux` mean flux and `rms_nonuniformity` on *each* observed package/sample face, and upper-vs-lower / left-vs-right drive balance across paired halfraums.
- **Backlighter foils & multi-foil / double-backlighter setups:** Per-foil `backlighter_signal` (photon flux drop `%`), foil irradiance `rms_nonuniformity` and peak/mean intensity, evaluated independently per foil and per timing window.
- **CBI (Crystal Backlighter Imager) & narrow-band crystal radiography:** Line-emission photon budget (`kshell_ce` + Bragg/filter chain), source-size / illuminated-area coverage, and **intensity-weighted source-centroid shift** on the foil (critical for narrow-band crystal alignment and rocking-curve throughput).
- **Direct drive / Polar Direct Drive (PDD):** `CoronaOperator` or hard-sphere illumination `mode_rms`, total absorbed power, and `ablation_pressure_proxy`.
- **Custom user/LLM criterion hook:** Any `ROECriterion(name, f; green, amber)` taking an evaluation/state so Discovery Science, ramp-compression, or ARC experiments plug directly into the same ROE matrix, chamber map, and HTML report.

Assign every beam and quad an explicit operational tier plus its **binding criterion** and **multi-outage budget rule** (e.g., *"up to 2 non-adjacent green drive beams; 0 CBI beams"*):
- **Green — Drop without asking Shot RI:** All active criteria across every role stay within tight auto-drop margins.
- **Amber — Ask Shot RI before dropping:** Crosses an auto-drop threshold, stays within hard shot-abort limits (or escalates only when combined with another outage).
- **Red — Critical (delay shot until fixed):** Violates a hard criterion alone (e.g., dedicated CBI or high-share backlighter quad, pushes `v_HS` over `mode1_tolerance`, breaks halfraum balance or package `Tr` floor).

### C.3 Visual Chamber Port Map & Shot RI ROE Export
- **GUI:** Color the 3-D beams and a 2-D NIF/OMEGA chamber port map by ROE tier (Red / Amber / Green), with role badges (Drive, Halfraum A/B, Foil 1/2, CBI). Clicking any beam/quad (or multi-selecting 2–3 hypothetical outages) displays all impacted criteria live (mode-1 kick vector, per-cavity `ΔTr`, package nonuniformity, foil/CBI photon drop & centroid shift).
- **Export & CLI:** One-click export (and `briareus roe input.{vrw,briar} --out dir/`) producing a ready-to-attach Shot RI ROE table (CSV, Markdown, and self-contained HTML with chamber port map) grouped by beam role and *Critical / Ask RI / Auto-Drop* tier with exact numerical justifications and pairwise caveats.

---

## Workstream D: LLM Scriptability, Automated Review Packets & Shareable HTML

Make Briareus trivially driveable by LLMs and batch scripts—and always visually auditable by humans, including reviewers who do not have Briareus installed.

### D.1 Bundled CLI & Zero-Install Script Runner (`BriareusGUI/src/app/app.jl` or `src/cli.jl`)
Extend the `briareus` command so the packaged app bundle doubles as a fast, zero-Julia-install automation tool using the exact same precompiled sysimage:
- `briareus [design.briar]` — open GUI (existing behavior).
- `briareus run script.jl [args...]` — execute a Julia script with `Briareus` (and `Controllers`) preloaded in ~1–2 s, no separate `juliaup` install required.
- `briareus review-1w input.{vrw,briar} --out dir/ [--threshold 1.0] [--json]` — run `talis_check` + `omega1_disposition` + render Daylight View and 2-D vector SVG slice figures for every failing or tight pair.
- `briareus roe input.{vrw,briar} --out dir/ [--json]` — run the multi-role Beam ROE Planner headlessly.
- `briareus diff a.{vrw,briar} b.{vrw,briar} --out dir/ [--json]` — compare two designs headlessly.
- `briareus talis input.{vrw,briar} --out report.csv [--json]` — headless TaLIS check for CI/batch pipelines.
- `briareus optimize ...` — run a saved optimization config headlessly and emit both `output.briar` and `output.vrw` / SST pointing files.

### D.2 Self-Contained Shareable HTML Shot-Review Packets
- Have `briareus review-1w`, `briareus roe`, `briareus diff`, and GUI **Export Review Report** emit a **single self-contained `.html` file** driven by the same composable role/criterion model (sortable TaLIS/1ω/ROE/diff tables, embedded Daylight View PNGs + pure-Julia vector SVG 2-D slices, chamber port map, and per-cavity/capsule/package/foil/CBI summaries—plus optional embedded Three.js mini-viewer once Stage B lands).
- **Zero-install committee distribution:** Any reviewer on the shot-review team opens the `.html` packet directly in their browser to inspect tight 1ω clearances, specular bounces, and Beam ROE justifications without installing Briareus.

### D.3 GUI ↔ LLM Handoff Loop
- **GUI → Script/LLM:** Add **File ▸ Export Julia script…** in `BriareusGUI/src/controllers/workspace.jl` emitting a clean standalone script reproducing the current session so a user can hand it straight to an LLM for sweeps, ROE updates, or optimization.
- **LLM → GUI:** Support `briareus design.briar --tab talis` (or `--tab roe` / auto-focusing a saved review view) so an LLM writes `result.briar` and the physicist opens it in one click to verify visually.
- **Deterministic headless 2-D SVG slices:** Emit the 2-D transverse $s^*$ cross-section slices (`B.3` / `D.1`) as pure vector SVGs directly from core Julia (zero GL context required).

### D.4 LLM Ergonomics Kit (`docs/src/llms.md`, `llms.txt`, optional MCP server)
- Publish a concise `llms.txt` and Claude skill file built against the post-cleanup public namespace covering:
  - Units, coordinate/angle conventions, and canonical workflows across platforms (hohlraum, halfraum/package, backlighter/CBI, direct drive): `read_visrad_workspace`, `talis_check`, `omega1_closest_approach`, `omega1_disposition`, `beam_roe`, `ROECriterion`, `Scenario` / `PhasedScenario`, `optimize_pointing`, `drive_history`, `write_visrad_workspace`, `save_session`.
  - Explicit "do/don't" traps (`dlambda` at 3ω vs 1ω, `history_state` for `WallLossAlbedo` scenarios, `VisRadProfile(; taper)`).
- *(Optional after CLI proves itself)* Add a thin `briareus mcp` stdio JSON-RPC server (~250 lines using `JSON`, zero new heavy deps) exposing `inspect_workspace`, `check_1w_and_talis`, `compute_beam_roe`, `diff_designs`, `render_1w_slices`, `run_drive`, `optimize_pointing`, and `save_briar_session`.

---

## Workstream E: Two-Stage Distribution Strategy (Windows & macOS First)

Eliminate installation friction for Windows and macOS desktop users first, then remove the Qt/GLMakie distribution bottleneck permanently.

### E.1 Stage A (Immediate — Weeks): Prebuilt Windows & macOS Relocatable Releases
- **GitHub Actions `PackageCompiler.jl` workflow (`.github/workflows/Release.yml`):**
  - **Primary targets:** **Windows x64** (`.zip` + InnoSetup/MSI installer) and **macOS arm64** (signed/notarized `.dmg` containing `Briareus.app`, resolving Qt6 framework rpaths in `reference/ROADMAP.md` item 2).
  - **Secondary target (optional):** Linux x64 `.tar.gz` and an offline Docker/Apptainer container (`dev/docker/Dockerfile`) for HPC/enclave users as needed.
  - Bundle both the GUI (`briareus`) and headless CLI (`briareus run`, `briareus review-1w`, `briareus roe`, `briareus diff`) on one sysimage.
  - Because `PackageCompiler` bakes compiled native code into the sysimage, drop `BriareusGUI/src/qt/first_frame.jl` replay on compiled builds and verify cold start drops from ~21 s to **< 4 s** on Windows and macOS.
  - Add clean-runner smoke tests in CI on Windows and macOS: unpack the release artifact with *no* Julia installed, run `briareus review-1w` on a fixture, and confirm zero missing DLLs/framework rpaths.
  - Add Windows Authenticode and macOS Developer ID signing + notarization steps (or explicit IT allowlisting docs) so SmartScreen and Gatekeeper do not block first-time users.

### E.2 Stage B (Strategic — Months): Web UI (Three.js) + Tauri Desktop App over Headless `Controllers`
Once Stage A validates the 1ω and ROE visual workflows with real VisRad users, replace `src/canvas/` + `src/qt/` + `src/qml/` with a Web frontend talking to the existing framework-free `Controllers` submodule over a lightweight Julia HTTP/WebSocket server:
1. **Protocol layer (`BriareusWeb` or `BriareusGUI/src/web/`):**
   - Expose `Controllers.DesignSession` state diffs, binary mesh/buffer endpoints (`Float32` vertices/normals/colors over WebSocket/HTTP—fast for 150k-triangle CAD), and action RPCs (`edit_part!`, `run_task!`, `import_visrad_workspace!`, `omega1_closest_approach`, `beam_roe`).
   - Zero `GLMakie`, zero `QML`, zero `QMLMakie`, zero `Qt6_*_jll`.
2. **Frontend (Three.js / WebGL2 + lightweight TS UI):**
   - Native GPU clipping planes (replacing CPU triangle clipping in `cutaway.jl`), 60 fps CAD orbiting, instanced/transparent 1ω closest-approach faces, orthographic Daylight/Grazing camera transitions, and crisp 2-D SVG/Canvas cross-sections & port maps.
   - Permanently eliminates all Win32/GLFW/Qt/AMD-driver/glyph-atlas workarounds (`BriareusGUI/src/qt/platform.jl`).
3. **Desktop-first + zero-install web channels from one build:**
   - **Native Windows & macOS desktop installers via Tauri (primary):** Bundles the static web assets + a slim `PackageCompiler` sysimage of *only* `Briareus + Controllers + HTTP` (~100–150 MB instead of ~700 MB, ~1–2 s launch, native WebView2 on Windows / WKWebView on macOS).
   - **Public web sandbox:** Hosted demo on public mock-up targets—physicists click a link in `README.md` and orbit a hohlraum + 1ω clearance view in 3 seconds.
   - **Lab-enclave / HPC portal (bonus from the same build):** Single container behind internal HTTPS / Open OnDemand when non-public `.vrw`/CAD data must stay on an enclave.

---

## Suggested Execution Sequence & Milestones (Post-Cleanup)

| Milestone | Scope | Concrete Acceptance Gate |
|---|---|---|
| **M1 — 1ω Engine & VisRad Parity** | Workstream A (A.1–A.4) | Real NIF `.vrw` auto-imports 1ω powers; `omega1_closest_approach` returns exact `(x_hw, x_beam, budget)`; cones, tori, and unsplit CAD faces classify unsigned local `|R₁|, |R₂|` accurately; facet-normal spread and exposed projected area/flux reported on 20° specular checks. |
| **M2 — 1ω Visual Wedge, Parity & Diff in GUI** | Workstream B (B.1–B.5) | 4-state hardware disposition coloring + closest-approach beam faces + 3-D caliper + specular reflection cone + one-click orthographic Daylight/Grazing Views + 2-D $s^*$ slice + safe-corridor map + role-aware VisRad parity & design diff working across hohlraum, halfraum, backlighter, and CAD targets. |
| **M3 — Multi-Platform Pre-Shot Beam ROE Planner** | Workstream C (C.1–C.3) | Composable `ROECriterion` outage engine (single-beam, single-quad, pairwise; linear + LPI) classifies every beam/quad into *Critical / Ask Shot RI / Auto-Drop* across hohlraums, halfraums/packages, multi-foil backlighters, CBI, direct drive, and custom metrics; visual chamber port map + exportable ROE sheet & `briareus roe` CLI. |
| **M4 — Windows/macOS Releases, HTML Packets & LLM Kit** | Workstream D + E.1 (+ A.5 / B.6 as desired) | Clean Windows x64 and macOS arm64 machines with no Julia install unpack CI releases → GUI launches in < 4 s; `briareus review-1w` and `briareus roe` emit self-contained shareable `.html` packets with vector SVG slices; `briareus run script.jl` runs LLM scripts in < 2 s; `llms.txt` published. |
| **M5 — Web/Tauri Frontend over `Controllers`** | Workstream E.2 | Feature-parity Three.js + TS frontend over headless `Controllers`; Windows & macOS Tauri desktop bundles < 150 MB with < 2 s launch; live public demo URL. |
