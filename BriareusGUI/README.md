# BriareusGUI.jl

**BriareusGUI** provides a desktop window for [Briareus.jl](../README.md) laser-pointing
and view-factor radiation design: build a hohlraum or backlighter target, point the
beams, see the laser deposition and x-ray drive, optimize the pointing, and check beam
clearance from the hardware and the TaLIS target-configuration rules, all in one
window, with VisRad pointing, beam-summary and benchmark-kit files in and out.

## Install the `briareus` command

No Julia knowledge is needed past these steps; they take a terminal, git, and about
ten minutes (most of it compiling, once).

1. Install Julia 1.12 or later through juliaup: on Windows `winget install --id
   9NJNWW8PVKMN -e` (Julia from the Microsoft Store), on macOS or Linux
   `curl -fsSL https://install.julialang.org | sh`. Open a new terminal afterward.
2. Clone the repository and install the app (the `launcher/` package) from the clone,
   then compile it:

   ```
   git clone https://github.com/stillyslalom/Briareus.jl
   julia -e "using Pkg; Pkg.Apps.develop(path=\"Briareus.jl/BriareusGUI/launcher\")"
   julia --project=Briareus.jl/BriareusGUI/launcher -e "using Pkg; Pkg.precompile()"
   ```

3. Put Julia's app folder on the `PATH` (the install prints a warning when it is
   missing). Windows, in PowerShell:

   ```powershell
   $p = [Environment]::GetEnvironmentVariable("Path", "User")
   [Environment]::SetEnvironmentVariable("Path", "$p;$HOME\.julia\bin", "User")
   ```

   macOS (zsh): `echo 'export PATH="$HOME/.julia/bin:$PATH"' >> ~/.zprofile`; Linux
   (bash): the same line appended to `~/.bashrc`. Open a new terminal afterward.
4. Add Briareus to the desktop and the applications menu: `briareus --shortcut`. On
   Windows this makes Start-menu and desktop shortcuts, on macOS
   `~/Applications/Briareus.app` (Launchpad and Spotlight list it), on Linux a menu entry.

Then the shortcut, or `briareus` in a terminal, opens the window on a fresh design, and
`briareus design.briareus.toml` on a saved session. A splash shows within about two
seconds while the window loads (about 20 s; the first start after an update compiles
for a few minutes). The command runs Julia with all threads (`--threads=auto`) and
returns when the window closes. A shortcut has no terminal: its output goes to
`%LOCALAPPDATA%\Briareus\briareus.log` (macOS: `~/Library/Logs/Briareus/`), and a
failure to start shows in the splash with the log's location.

To update, pull the clone and compile again:

```
git -C Briareus.jl pull
julia --project=Briareus.jl/BriareusGUI/launcher -e "using Pkg; Pkg.resolve(); Pkg.precompile()"
```

To uninstall, run `julia -e "using Pkg; Pkg.Apps.rm(\"BriareusLauncher\")"`, delete the
shortcuts, and delete the clone.

The app is installed from a clone rather than with `pkg> app add <url>`: Briareus is
not registered, and the GUI finds it through the `[sources]` path to the parent
checkout, which a subdirectory installed from a URL does not have. The launcher is a
package of its own so that the splash shows before BriareusGUI loads: `julia -m` loads
the app's whole package, and BriareusGUI's dependencies take 15 s to load.

## Open the Briareus window from Julia

Install Briareus and the GUI in the same Julia environment. While Briareus is not
registered, develop both from the repository root:

```julia
pkg> dev . ./BriareusGUI
```

Start Julia with several threads (`julia -t auto`), return to the Julia prompt with
Backspace, and open the window:

```julia
using BriareusGUI
sess = briareus()                          # or briareus(resolution = :fine)
```

The window needs a graphical session. It is laid out like a CAD tool, in the manner
of VisRad, rather than as a sequence of steps: everything is reachable at any time,
and each analysis says when it is out of date.

- **Model tree** (left): the target's components (a hohlraum's wall, caps, LEHs and
  capsule; diagnostic windows; free parts) with a swatch of their color in the 3D view
  (the tree is the color key), the beams by group, quad and beam with on/off boxes, and
  the analyses (drive, optimization, clearance, TaLIS) with their state. Selecting a
  node selects it in the 3D view and shows it in Properties; Ctrl + click adds
  components or beams to the selection and Shift + click a range, and Properties then
  edits them together; a double click fits the view to it; right click for its
  actions, among them *Show* (draw it or not) and *Include in calculations*: a
  component left out is drawn as a gray ghost and taken out of the target (an open
  enclosure then), and beams left out are switched off.
- **Properties** (right): the selection's settings.
  - *Target:* template (hohlraum with capsule, backlighter foil, or custom parts
    only), dimensions and wall material, and the mesh: presets from a work-precision
    study on the NIF hohlraum (`examples/mesh_convergence.jl`: coarse, standard and
    fine give capsule P2/P0 and P4/P0 to about 0.5, 0.1 and 0.02 percentage points;
    fine is converged) or element sizes in mm for the wall and the capsule.
  - *Windows* (`WallWindow`): azimuth, height, round (the default; exact disk area) or
    rectangular, open or closed with their own albedo and conversion. Overlapping windows, or one past the end of the wall, are
    refused with the reason.
  - *Parts:* spheres, disks, washers, rectangles, cylinders and cones placed anywhere
    (center, axis polar and azimuth angles and a rotation, as VisRad orders them),
    with their own albedo, conversion, sidedness, aperture flag and element size,
    e.g. a shield, a stalk, a detector aperture or a second foil. A part is placed in
    the chamber's frame, at an entrance hole (z pointing out of the hohlraum: one
    placement mirrors between the two), at the backlighter foil, or in another part's
    frame, and moves with it when the target's dimensions or that part change. Its
    position is given in Cartesian, cylindrical or spherical coordinates of that
    frame, chosen per part. A *VisRad object* part is a component of an imported
    VisRad workspace: a VisRad rectangle, sphere or cylindrical hohlraum keeps VisRad's
    own grid (its element table) unless the part is given an element size. A *mesh* part reads an
    STL file (a CAD export; Target ▸ Add part ▸ Mesh from STL file), scaled to mm, with
    its windings repaired (`orient_mesh`) and its surface remeshed into elements of
    the part's size (`remesh_mesh`); a multi-body export (an ASCII STL of named solids)
    gives one part per solid. A part's *Role* names it for the analyses: *Capsule*
    (the observed surface; a CAD capsule is replaced in the calculations by a sphere
    fitted to it, `fit_sphere`, so a fill tube does not count) or *Backlighter foil*
    (the foil objective's spot), and *Thin film* leaves it out of the drive's cavity
    (below).
  - *Beams* (a group, quad, beam, or any Ctrl-selected set): aim relative to the ports
    (z, r, tangential offset) or a common aim point (x, y, z, as VisRad points a
    cluster at a reference point), or *Aim at* a component's centroid (the foil, an
    entrance hole, a part), which the beams follow when it moves; defocus, power,
    on/off, and where the axis lands;
    "Look along beam" views the target from its port. Beams belong to groups, as
    VisRad's clusters: the cone by default, or any named set ("backlighter"), moved
    with Move to group; a group is edited, switched and optimized together (its
    members keep their cone for the laser–plasma models). A NIF quad held as one beam
    splits into its four beams (own ports, a quarter of the power each) to point them
    apart, and merges back.
  - *Beams:* the facility (Beams ▸ Facility: NIF, OMEGA, OMEGA EP or custom), which
    decides the presets offered (NIF quads or beams at the published cone aims, the
    OMEGA hohlraum or backlighter sets; OMEGA EP's port geometry is not built in) and
    the checks (the TaLIS rules are NIF's), and how beams are drawn.
  - *Analyses:* the settings and Run of the drive, the optimization, the clearance
    check and the TaLIS rules (as described below).
- **3D view** (center), with CAD navigation: the wheel zooms toward the point under
  the pointer (it stays under the pointer), a left or middle drag orbits about the
  selection or the point under the pointer, a right drag, Shift + drag or Ctrl +
  middle drag pans the grabbed point with the pointer, a click selects a component
  or beam, a double click fits it, `f` fits everything and the arrow keys orbit. Near
  and far planes follow the eye, so nothing is clipped however close the view.
  Standard views, a perspective or orthographic projection, turntable (z up) or free
  rotation, a cutaway (right-click the Cutaway tool for its width, its azimuth or the
  side facing the view, and which components it cuts), show/hide per component,
  View ▸ Look from port (a beam's port or any θ, φ), and an orientation triad. Beams are drawn as rays (each axis from the lens side to the
  surface it lands on, the aim point marked) or as volumes (the 5 % contour of the
  spot carried along the f-cone through best focus, each edge cut where it meets the
  target, so a beam clipping an LEH lip or the capsule shows), colored by group or by
  their worst clearance or TaLIS status (beam groups take the tree's colors). The
  target is colored by component (keyed in the tree) or by
  absorbed laser intensity, x-ray radiosity or incident flux (W/cm²), which computes
  the drive when it is missing. The capsule θ–φ flux map (or the foil face) shares
  the canvas, and so does an x-ray image (Compute ▸ X-ray image along the view,
  `xray_image!`; unchecking it returns to the view before): the radiance of the surfaces seen along the 3D view's direction,
  framed as the view frames them, through transparent apertures (a hohlraum is seen
  inside through its entrance holes). It follows the drive as beams are edited. View ▸
  Map or image beside the 3D view shows either of them in a second canvas next to
  the 3D target instead of in its place.
- **CAD workspaces.** A NIF CAD workspace brings 15 parts and 2 × 10⁵ facets, most of
  them mounting hardware reaching 10² mm from the target. On import, a CAD part that
  reaches farther than twice the farthest aim point of the beams VisRad fires is
  drawn only (from its file's facets, not remeshed) and left out of the
  calculations; the status names these parts and *Include in calculations* puts one
  back (remeshed on workers). Beams VisRad leaves off are not drawn when most are
  (View ▸ Beams ▸ Show beams that are off). The drive of a workspace with CAD runs on
  its cavity around TCC (Drive ▸ Cavity, `set_cavity!`, `Briareus.enclosure`): the
  elements bounding the free space around the point, without the assembly's outer
  faces (≈ 7,000 of 20,000 elements on that workspace); a thin film modeled as a solid
  that seals the point into a pocket (a capsule tent) is found and left out, and the
  Drive properties list what the cavity kept of each part. The cavity's view factors
  are balanced once aperture parts (facing into the cavity) close its entrance holes,
  under 1 % of its rays leaking (`CAVITY_MAX_LEAK`); the status says which. A target
  whose x-ray view factors would not fit in half the memory (`radiation_size`:
  ≈ 160,000 elements on 32 GB) gets the laser deposition alone, with the reason and the fix in the status.
  A CAD part whose element size changes (its own, the target's, an undo or a loaded
  session bringing back a size not meshed before) is remeshed on worker threads: the
  edit takes effect at once with the part drawn from its file's facets ("remeshing …"
  in the status), and a drive asked for meanwhile starts when the parts are ready.
  Building a target of 10⁵ facets (its BVH) takes half a second, so a large target not
  built before is built on a worker as well: the edit's settings apply at once, the
  view shows the previous target until the new one is ready ("building the target"),
  and the drive and the checks wait for it; so does a CAD capsule's sphere fit
  ("fitting the capsule to …"). A CAD part that is one cylinder or sphere as a whole (a
  tube, a shell; `fit_surfaces` on its file's facets) is integrated on that surface in
  the view factors under its own name; parts of several faces stay flat, and the
  drive's status names flat parts carrying more than 50 facets per element (a curved
  face tessellated in slivers: slow view factors).
- **Results** (bottom tabs):
  - *Beam table:* every cone, quad or beam with its pointing, defocus, power and
    on/off, edited in place.
  - *Drive:* the view factors are computed once per target and beam ports (on a
    worker; seconds on the standard mesh), then beam edits re-evaluate at once: the
    capsule mode amplitudes ε₁…ε₄ and P_l/P₀, the foil spot statistics, the laser
    power accounting, and **Write VisRad benchmark kit…** (`write_visrad_kit`: SST
    pointing, beam summary, beam table, a README with the geometry and the VisRad
    steps, and the Briareus results), reusing the drive's view factors.
  - *Time:* the drive in time. **Pulse table…** fires the beams with a pulse table
    (`read_pulse_table`: power per beam against time per cone group, quad or beam;
    disabled beams stay off), and **N221204 shape** with N221204's delivered pulse shape
    (`n221204_pulse`), each beam reaching its table power at the peak. **Run history**
    solves the drive at each time of the pulses with the drive's view factors
    (`drive_history`, on a worker with progress and cancel; the drive is computed first
    when it is missing). Four charts, colored as the model tree keys the beam groups and
    components: *Pulses* (each group's power and the total, TW), *Drive history*
    (capsule Tᵣ and a synthetic Dante 1, `dante(dante1(), …)`, in eV; a₂/a₀ and a₄/a₀ in
    % on the right axis), *Energy* (the laser energy delivered by each time and where it
    went, kJ, `energy_budget`: heat in the walls and the capsule, x-rays out of the
    entrance holes and windows) and *Retune* (below). A click or drag on the chart, or
    the slider under it, moves the time cursor: the 3D view is colored by that step of
    the history (by the incident x-ray flux when it was colored by component), on one
    color scale for the whole history, labeled with the time; each step's fields are
    kept with the history, so scrubbing computes nothing. **Now** returns to the drive
    at the table's powers. A beam, pulse or target edit marks the history out of date
    and the 3D view leaves it until it is run again.
  - *Optimization:* free aim, defocus and power per cone, quad or beam within ±
    bounds; minimize capsule asymmetry (axisymmetric P_l or all-m 3D modes, optionally
    holding the drive) or foil-spot nonuniformity with `optimize_pointing` in the
    background, with progress and cancel (canceling keeps the best design so far). A
    capsule keep-out (`SphereKeepOut`, on by default) keeps every beam's worst-case
    envelope off the capsule; an LEH keep-out (`ApertureKeepOut`) keeps it off the LEH
    rims. The run stops at the iteration limit or the evaluation limit, whichever
    comes first. Compare the modes, clipping and capsule clearance before and after,
    read the parameter table, and apply the result to the beam table.
  - *Clearance:* worst-case beam-to-hardware clearance (`clearance`) over pointing,
    positioning, hardware and focus errors, on the components you choose, in the
    background with progress per batch of beams and cancel (keeping the beams already
    checked); rows colored pass/review/fail, the closest beam per component, and CSV
    export (`write_clearance_csv`).
  - *TaLIS:* the NIF target-configuration checks (`talis_check`): 3ω capture,
    incidence and clearance, 1ω clearance and unconverted-light management, and the
    TAS stay-in zone. Choose the capture components and the components to check, the
    rule tolerance and how it becomes a distance (`pivot`: r from TCC, from the beam's
    aim point, the larger of the two, or from a given point); the error stack and
    sampling are the clearance check's. The report shows a per-rule summary, every
    check colored by status with its margin, and the notes: the tolerance
    interpretation is marked ASSUMED and stays at the top, then the error set, the
    checks not evaluated, and the requirements out of scope (`TALIS_OUT_OF_SCOPE`).
    CSV export with `write_talis_csv`.
- **Menus and toolbar:** new targets, parts and windows; beam presets; beam tables
  (`read_beam_table`/`write_beam_table`) and VisRad files in and out (a NIF SST
  pointing CSV, `read_visrad_sst`: aim points and defocus; a Laser Beam Summary CSV,
  `read_visrad_summary`: pointing, defocus, power and spot; beams off or absent in the
  file are switched off at their power, and the readers' warnings are listed); VisRad
  workspaces (File ▸ Import VisRad workspace: its components as VisRad object parts
  and its beams, read and meshed on worker threads with progress; File ▸ Export
  VisRad workspace writes the beams back into a copy of the imported one,
  `write_visrad_workspace`); File ▸ Export target as STL (`write_stl`, mm); File ▸ Export results (the drive's
  per-element laser, radiosity and flux, power per component and beam, and capsule
  modes as CSV, `write_results`, with the session that produced them); views; and the
  analyses.
- **Retune to a reference drive** (Compute ▸): given the pulse table of a reference
  shot that fired every beam of the table (`read_pulse_table`: power per beam against
  time per cone group, quad or beam), the enabled beams' powers are retuned at each
  time of the table to restore the reference's capsule drive level and its modes
  through l = 4 (`retune_drive`), within a limit relative to the reference powers and
  optionally staying near them. The disabled beams are the ones the shot loses. The
  new pulses are written as a pulse table; the dialog reports the worst-time mode 1
  and Tᵣ change before and after, from drive histories of both, and the Time tab plots
  them: each group's power on the reference pulses (dashed) and on the retuned ones,
  and ε₁ before and after.
- **Sessions, undo and redo:** File ▸ Save session (Ctrl + S) writes the whole design
  (target, windows, parts, beams with their groups and on/off state, the analyses'
  settings and the view) to a plain TOML file (`save_session`), and File ▸ Open session
  (Ctrl + O) reads it back (`load_session!`); results are computed again. Edit ▸ Undo
  (Ctrl + Z) and Redo (Ctrl + Y) step through every design edit, whether typed, from a
  menu, from a file or an applied optimization; the view (coloring, selection, camera)
  is not part of the history.

`briareus` returns its session (a `DesignSession`) when the window closes; the
controllers in `BriareusGUI.Controllers` run the same steps from a script:

```julia
using BriareusGUI, BriareusGUI.Controllers
sess = DesignSession(resolution = :coarse)       # hohlraum + NIF quads
add_window!(sess.target; phi = 60.0, open = true) # a 1 × 1 mm hole in the wall
add_part!(sess.target, :disk; z = 7.5, r1 = 2.5, theta = 180)   # a shield above the top LEH
edit_mesh!(sess.target, :capsule_h, "0.05")      # element sizes in mm
edit_beam_node!(sess.beams, :group, "inner23", :z, "3.2")   # inner 23.5° cone aim z (mm)
import_visrad_sst!(sess.beams, "pointing.csv")   # or apply a VisRad pointing file
build_scenario!(sess)                            # view factors + evaluation
drive_summary(sess.drive.result[])
write_visrad_kit!(sess, "visrad_kit")            # benchmark kit for VisRad
start_optimize!(sess; spawn = false)             # blocks; spawn = true runs on a worker
accept_optimum!(sess)
start_clearance!(sess; spawn = false)
worst_rows(sess.review)
start_talis!(sess; spawn = false)                # TaLIS rules, same error stack
talis_rule_rows(sess.talis)
talis_notes(sess.talis)                          # the ASSUMED tolerance caveat first
export_results(sess, "results")                  # per-element and per-mode CSV
retune_to_reference!(sess, "reference.csv", "retuned.csv")  # disabled beams are the dropped ones
load_pulses!(sess, "reference.csv")              # or use_pulse_shape!(sess): N221204's shape
start_history!(sess)                             # drive_history on the drive's view factors
time_chart(sess, :drive)                         # what the Time tab plots (also :pulses, :energy, :retune)
set_time_cursor!(sess, 6.0)                      # the 3D view at the history step nearest 6 ns
save_session("design.briareus.toml", sess)       # the design as TOML; load_session! reads it
record!(sess); undo!(sess)                       # undo/redo steps (a window records them itself)
```

`edit_*` functions take text-field text and return `false` with the reason in the
controller's `error` when the geometry or a check rejects the value.

## Development

On Julia ≥ 1.11 the `[sources]` entry in `Project.toml` couples this package to the
parent Briareus checkout:

```
julia --project=BriareusGUI -t auto -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

The tests drive every controller without a window and render the canvases offscreen
(they need a GL context; on a headless Linux machine wrap them in `xvfb-run`). The Qt
window test runs in its own process when `BRIAREUSGUI_QT_TESTS=true` is set on a
machine with a display (`test/qt_window.jl`, also runnable alone; set
`BRIAREUSGUI_SHOTS` to keep its window images). `BRIAREUSGUI_NO_WORKLOAD=true` skips
the precompile workload for faster rebuilds while developing.

Time to first window: the workload builds the window's shell and item models (QML's Qt
application exists while precompiling), but it has no GL context, so the first frame's
rendering code on the window's screen type comes from `precompile` signatures recorded
by `julia dev/gui_first_frame.jl` into `src/qt/first_frame.jl` and replayed by the
workload. Rerun the recorder after updating GLMakie, QMLMakie or Makie, or after
changing what the canvases plot. On a 16-thread Windows machine the first frame came
at 50 s without either, 21 s with both; the rest is loading the package images (15 s,
GLMakie and QML about equally) and QML's `__init__`. The launcher (`launcher/`, the
`briareus` command) has its own tests: `julia --project=BriareusGUI/launcher -e 'using
Pkg; Pkg.test()'`.

The package follows HammerheadGUI's architecture: all state and logic lives in the
framework-free `Controllers` submodule (Briareus + Observables, testable without a GL
context), the Makie canvases (`src/canvas/`) create every plot up front and only update
their inputs (the 3D view sets its own camera matrices from a `CadCamera`, which takes
its gestures from the scene's mouse and keyboard events and picks by ray casting into
the target's BVH), and `src/qt/shell.jl` bridges controllers to the QML pages (`src/qml/`)
through one property map, item models and `br_*` callbacks. Background jobs run on
worker threads and hand their observable updates to a queue that the window's 16 ms
tick applies on the GUI thread. Long core calls report progress and cancel through
their own callbacks (`optimize_pointing`'s per objective evaluation, `clearance`'s per
batch of beams), which run on the job's task. `talis_check` has no callback, so its
Cancel drops the result instead. GLMakie windows (or offscreen renders) and the Qt
window cannot share a process.
