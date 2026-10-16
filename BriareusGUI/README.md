# BriareusGUI.jl

**BriareusGUI** provides a desktop window for [Briareus.jl](../README.md) laser-pointing
and view-factor radiation design: build a hohlraum or backlighter target, point the
beams, see the laser deposition and x-ray drive, optimize the pointing, and check beam
clearance from the hardware, all in one window.

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
  the main target; coarse, standard or fine mesh. The viewer colors the mesh by
  component.
- **Beams:** NIF quads or beams (pointed at the published cone aims), the OMEGA
  hohlraum or backlighter sets, or a beam table imported from CSV
  (`read_beam_table`); edit aim (z, r, tangential offset), defocus, power and on/off
  per cone, quad or beam, and export the table (`write_beam_table`). The viewer draws
  each beam's axis from its port to its aim point.
- **Drive:** compute the view factors once per target and beam set (on a worker;
  seconds on the coarse mesh, tens of seconds on the standard one), then see the
  target colored by absorbed laser intensity, x-ray radiosity or incident flux
  (W/cm²), the capsule flux map over θ–φ (or the foil's face), the mode amplitudes
  ε₁…ε₄ and P_l/P₀, and the laser power accounting. Beam edits re-evaluate at once.
- **Optimize:** free aim, defocus and power per cone, quad or beam within ± bounds;
  minimize capsule asymmetry (axisymmetric P_l or all-m 3D modes, optionally holding
  the drive and keeping beams off the LEH rims) or foil-spot nonuniformity with
  `optimize_pointing` in the background, with progress and cancel (canceling keeps the
  best design so far). Compare the modes before and after, read the parameter table,
  and apply the result to the beam table.
- **Review:** worst-case beam-to-hardware clearance (`clearance`) over pointing,
  positioning, hardware and focus errors, on the components you choose; rows colored
  pass/review/fail, the closest beam per component, beams colored by status in the
  viewer, and CSV export (`write_clearance_csv`).

`briareus` returns its session (a `DesignSession`) when the window closes; the
controllers in `BriareusGUI.Controllers` run the same steps from a script:

```julia
using BriareusGUI, BriareusGUI.Controllers
sess = DesignSession(resolution = :coarse)       # hohlraum + NIF quads
edit_beams!(sess.beams, 1, :z, "3.2")            # inner 23.5° cone aim z (mm)
build_scenario!(sess)                            # view factors + evaluation
drive_summary(sess.drive.result[])
start_optimize!(sess; spawn = false)             # blocks; spawn = true runs on a worker
accept_optimum!(sess)
start_clearance!(sess; spawn = false)
worst_rows(sess.review)
```

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
tick applies on the GUI thread. GLMakie windows (or offscreen renders) and the Qt
window cannot share a process.
