# BriareusGUI.jl

**BriareusGUI** provides a desktop window for [Briareus.jl](../README.md) laser-pointing
and view-factor radiation design: build a hohlraum or backlighter target, point the
beams, see the laser deposition and x-ray drive, optimize the pointing, and check beam
clearance from the hardware and the TaLIS target-configuration rules, all in one
window, with VisRad pointing, beam-summary and benchmark-kit files in and out.

## Open the Briareus window

Install Briareus and the GUI in the same Julia environment. While Briareus is not
registered, develop both from the repository root:

```julia
pkg> dev . ./BriareusGUI
```

Start Julia with several threads (`julia -t auto`), return to the Julia prompt with
Backspace, and open the window:

```julia
using BriareusGUI
sess = briareus()                          # or briareus(resolution = :coarse)
```

The window needs a graphical session. It walks through one design in five steps, with
a 3D viewer (drag to rotate, wheel to zoom, a cutaway of the hohlraum shell) that
follows the step:

- **Target:** a hohlraum (length, radius, LEH and capsule radii, wall albedo and
  laser → x-ray conversion) or a backlighter foil facing TCC with a black sphere for
  the main target; coarse, standard or fine mesh. A hohlraum takes diagnostic windows
  in its wall (`WallWindow`): a small table of windows (azimuth φ, height z, width,
  height; closed with their own albedo and conversion, or open holes) that rebuilds
  the target. A window that overlaps another or runs past the end of the wall is
  refused with the reason. The viewer colors the mesh by component (each window is
  its own component, `window_1`, …) and outlines the windows.
- **Beams:** NIF quads or beams (pointed at the published cone aims), the OMEGA
  hohlraum or backlighter sets, or a beam table imported from CSV
  (`read_beam_table`); edit aim (z, r, tangential offset), defocus, power and on/off
  per cone, quad or beam, and export the table (`write_beam_table`). VisRad files
  apply to the loaded table: a NIF SST pointing CSV (`read_visrad_sst`: aim points and
  defocus) or a Laser Beam Summary CSV (`read_visrad_summary`: pointing, defocus,
  power and spot). Beams that are off or absent in the file are switched off at their
  power, and the readers' warnings (unmatched beams, port-angle mismatches) are listed
  under the table. Both files export too (`write_visrad_sst`, `write_visrad_summary`).
  The viewer draws each beam's axis from its port to its aim point.
- **Drive:** compute the view factors once per target and beam set (on a worker;
  seconds on the coarse mesh, tens of seconds on the standard one), then see the
  target colored by absorbed laser intensity, x-ray radiosity or incident flux
  (W/cm²), the capsule flux map over θ–φ (or the foil's face), the mode amplitudes
  ε₁…ε₄ and P_l/P₀, and the laser power accounting. Beam edits re-evaluate at once.
  **Write VisRad benchmark kit…** writes the case into a folder (`write_visrad_kit`:
  SST pointing, beam summary, beam table, a README with the geometry and the VisRad
  steps, and Briareus's results to compare against), reusing the drive's view factors.
- **Optimize:** free aim, defocus and power per cone, quad or beam within ± bounds;
  minimize capsule asymmetry (axisymmetric P_l or all-m 3D modes, optionally holding
  the drive and keeping beams off the LEH rims) or foil-spot nonuniformity with
  `optimize_pointing` in the background, with progress and cancel through its callback
  (canceling keeps the best design so far). Compare the modes before and after, read
  the parameter table, and apply the result to the beam table.
- **Review**, two tabs:
  - *Clearance:* worst-case beam-to-hardware clearance (`clearance`) over pointing,
    positioning, hardware and focus errors, on the components you choose, in the
    background with progress per batch of beams and cancel (keeping the beams already
    checked); rows colored pass/review/fail, the closest beam per component, and CSV
    export (`write_clearance_csv`).
  - *TaLIS rules:* the NIF target-configuration checks (`talis_check`): 3ω capture,
    incidence and clearance, 1ω clearance and unconverted-light management, and the
    TAS stay-in zone. Choose the capture components (by default the hohlraum wall or
    the backlighter foil) and the components to check, the rule tolerance and how it
    becomes a distance (`pivot`: r from TCC, from the beam's aim point, the larger of
    the two, or from a given point); the error stack and sampling are the Clearance
    tab's. The check runs in the background (Cancel drops its result). The tab shows a
    per-rule summary, every check colored by status with its margin, and the report's
    notes: the tolerance interpretation is marked ASSUMED and stays at the top, then
    the error set, the checks not evaluated, and the requirements out of scope
    (`TALIS_OUT_OF_SCOPE`). CSV export with `write_talis_csv`.

  The viewer colors the beams by their worst status in the tab shown.

`briareus` returns its session (a `DesignSession`) when the window closes; the
controllers in `BriareusGUI.Controllers` run the same steps from a script:

```julia
using BriareusGUI, BriareusGUI.Controllers
sess = DesignSession(resolution = :coarse)       # hohlraum + NIF quads
add_window!(sess.target; phi = 60.0, open = true) # a 1 × 1 mm hole in the wall
edit_beams!(sess.beams, 1, :z, "3.2")            # inner 23.5° cone aim z (mm)
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
```

`edit_*` functions take text-field text and return `false` with the reason in the
step's `error` when the geometry or a check rejects the value.

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

The package follows HammerheadGUI's architecture: all state and logic lives in the
framework-free `Controllers` submodule (Briareus + Observables, testable without a GL
context), the Makie canvases (`src/canvas/`) create every plot up front and only update
their inputs, and `src/qt/shell.jl` bridges controllers to the QML pages (`src/qml/`)
through one property map, item models and `br_*` callbacks. Background jobs run on
worker threads and hand their observable updates to a queue that the window's 16 ms
tick applies on the GUI thread. Long core calls report progress and cancel through
their own callbacks (`optimize_pointing`'s per objective evaluation, `clearance`'s per
batch of beams), which run on the job's task. `talis_check` has no callback, so its
Cancel drops the result instead. GLMakie windows (or offscreen renders) and the Qt
window cannot share a process.
