# BriareusGUI.jl

**BriareusGUI** provides a desktop window for [Briareus.jl](../README.md) laser-pointing
and view-factor radiation design: build a hohlraum or backlighter target, point the
beams, see the laser deposition and x-ray drive, optimize the pointing, and check beam
clearance from the hardware and the TaLIS target-configuration rules, all in one
window, with VisRad pointing, beam-summary, workspace and benchmark-kit files in and out.

![The Briareus window](../docs/src/assets/gui/overview.png)

The [Desktop application](https://stillyslalom.github.io/Briareus.jl/dev/gui/) page of
the documentation describes the window, a first session and every panel.

## Install the `briareus` command

1. Install Julia 1.12 or later through juliaup: on Windows `winget install --id
   9NJNWW8PVKMN -e`, on macOS or Linux `curl -fsSL https://install.julialang.org | sh`.
   Open a new terminal afterward.
2. Install the app (`BriareusLauncher`, in `launcher/`) and compile it, about ten
   minutes once: `julia -e "using Pkg; Pkg.Apps.add(\"BriareusLauncher\")"`, or
   `pkg> app add BriareusLauncher` at the Julia prompt.
3. Put Julia's app folder (`~/.julia/bin`) on the `PATH`; the install prints a warning
   when it is missing. Windows, in PowerShell:

   ```powershell
   $p = [Environment]::GetEnvironmentVariable("Path", "User")
   [Environment]::SetEnvironmentVariable("Path", "$p;$HOME\.julia\bin", "User")
   ```

   macOS (zsh): `echo 'export PATH="$HOME/.julia/bin:$PATH"' >> ~/.zprofile`; Linux
   (bash): the same line appended to `~/.bashrc`.
4. `briareus --shortcut` adds Briareus to the desktop and the applications menu.

`briareus` opens the window on a fresh design, and `briareus design.briareus.toml` on a
saved session. A shortcut has no terminal: its output goes to
`%LOCALAPPDATA%\Briareus\briareus.log` (macOS: `~/Library/Logs/Briareus/`).

## Open the window from Julia

```julia
pkg> add BriareusGUI

julia> using BriareusGUI   # in a session started with `julia -t auto`

julia> sess = briareus()   # returns the session when the window closes
```

## Development

The `[sources]` entries in `Project.toml` and `launcher/Project.toml` couple each
package to the checkout above it (Pkg honors them for the active project only):

```
julia --project=BriareusGUI -t auto -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia -e 'using Pkg; Pkg.Apps.develop(path="BriareusGUI/launcher")'   # `briareus` from the checkout
```

The tests drive every controller without a window and render the canvases offscreen
(they need a GL context; on a headless Linux machine wrap them in `xvfb-run`). The Qt
window test runs in its own process when `BRIAREUSGUI_QT_TESTS=true` is set on a
machine with a display (`test/qt_window.jl`, also runnable alone; set
`BRIAREUSGUI_SHOTS` to keep its window images). `BRIAREUSGUI_NO_WORKLOAD=true` skips
the precompile workload for faster rebuilds while developing. The launcher has its own
tests: `julia --project=BriareusGUI/launcher -e 'using Pkg; Pkg.test()'`.
`julia dev/gui_screenshots.jl` regenerates the documentation's screenshots.

Time to first window: the workload builds the window's shell and item models (QML's Qt
application exists while precompiling), but it has no GL context, so the first frame's
rendering code on the window's screen type comes from `precompile` signatures recorded
by `julia dev/gui_first_frame.jl` into `src/qt/first_frame.jl` and replayed by the
workload. Rerun the recorder after updating GLMakie, QMLMakie or Makie, or after
changing what the canvases plot. On a 16-thread Windows machine the first frame came
at 50 s without either, 21 s with both; the rest is loading the package images (15 s,
GLMakie and QML about equally) and QML's `__init__`. The launcher is a package of its
own so that its splash shows before BriareusGUI loads: `julia -m` loads an app's whole
package, and BriareusGUI's dependencies take 15 s to load.

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
