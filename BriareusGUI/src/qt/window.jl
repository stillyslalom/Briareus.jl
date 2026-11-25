# This file defines `briareus`, which opens the window.

# ---------------------------------------------------------------- the window

# The session `briareus()` starts from when none is given (`preset = nothing`: the default).
_launch_session(; starter, resolution, preset, beams) =
    preset === nothing ? DesignSession(; starter, resolution, beams) :
    DesignSession(; starter, resolution, preset, beams)

"""
    briareus(sess = DesignSession(); starter = :hohlraum, resolution = :standard,
             preset = nothing, beams = nothing, build = false, shown = nothing) -> DesignSession

Open the Briareus window and return its session when the window closes. The window is
a design tool in the manner of VisRad; it is not a sequence of steps:

- **Model tree** (left): the target's components (a hohlraum's wall, caps, LEHs and
  capsule, its diagnostic windows, free parts; each with its display mode), the beams by
  cone, quad and beam (with VisRad's dots: 1ω, 2ω and 3ω light in the calculations and
  drawn), and the analyses with their state. Selecting a node selects it in
  the 3D view and shows it in Properties.
- **Properties** (right): the selection's settings: the target template, dimensions,
  materials and mesh element sizes; a window or part's size, placement and material;
  a cone, quad or beam's pointing, defocus and power; an analysis's settings and Run.
- **3D view** (center): CAD-style navigation (wheel zooms toward the pointer, left or
  middle drag orbits about the selection or the point under the pointer, right drag
  or Shift-drag pans, double click fits the clicked part, `f` fits all), standard
  views, a perspective or orthographic projection, beams as rays to where they land
  or as volumes, coloring by component or by a drive field (computed on request), a
  cutaway (a pie slice or a section plane, clipped exactly), an orientation triad, and
  a placement triad on the selected part (arrows move it, rings turn it, one edit per
  drag).
  The capsule θ–φ (or foil) map shares the canvas.
- **Results** (bottom): the beam table, the drive (modes, power accounting), time
  (pulses, the drive history and the retune plotted against time, with a time cursor
  that colors the 3D view by that step), the optimization, the clearance check and the
  TaLIS rules.
- **Menus and toolbar**: new targets, parts and windows, beam presets and imports
  (CSV, VisRad SST and Laser Beam Summary) and exports, views, and the analyses.

The session starts from `sess`, or a fresh `DesignSession(; starter, resolution,
preset, beams)`; `build = true` starts computing the drive at once, and `shown()`, when
given, runs once on the window's first tick (the `briareus` command closes its splash
then). The call blocks while the window is open (Qt runs its event loop on this thread); start Julia with
several threads (`julia -t auto`) so background jobs leave the window responsive.
Closing the window cancels an optimization, clearance or TaLIS check in progress.

The window throws an `ArgumentError` when this Julia session already has a GLMakie
screen (a GLMakie window or an offscreen `colorbuffer`): GLMakie's windows and the Qt
window cannot share a process.
"""
function briareus(sess::Union{Nothing,DesignSession} = nothing; starter::Symbol = :hohlraum,
                  resolution::Symbol = :standard, preset = nothing, beams = nothing,
                  build::Bool = false, shown = nothing)
    sess === nothing && (sess = _launch_session(; starter, resolution, preset, beams))
    _SHELL[] === nothing || throw(ArgumentError("a BriareusGUI window is already open"))
    _check_no_glfw_screens()
    # From here on background work goes to the window's queue.
    queue = Channel{Any}(Inf)
    previous_deliver, previous_spawn = sess.deliver[], sess.spawn[]
    sess.deliver[] = f -> put!(queue, f)
    sess.spawn[] = true
    # Windows ticks timers at 15.6 ms unless a process asks for finer: the window's 16 ms
    # tick then fires every 31 ms (30 frames per second while orbiting). Ask for 1 ms while
    # the window is open, as browsers do.
    finer = Sys.iswindows() && ccall((:timeBeginPeriod, "winmm"), UInt32, (UInt32,), 1) == 0
    try
        _qt_init!()
        _warm_glyph_atlas!()
        sh = BriareusShell(sess; queue)
        _warm_canvases!(sh)
        _SHELL[] = sh
        _CLOSE_REQUESTED[] = false
        _register_qml_functions()
        build && build_scenario!(sess)
        _ICON_TRIES[] = 120                    # about two seconds of ticks
        _ON_SHOWN[] = shown
        loadqml(joinpath(QML_DIR, "BriareusWindow.qml"); app = sh.app,
                (k => m for (k, m) in sh.models)...)
        exec()
    finally
        finer && ccall((:timeEndPeriod, "winmm"), UInt32, (UInt32,), 1)
        _SHELL[] = nothing
        _ON_SHOWN[] = nothing
        Controllers._abandon_jobs!(sess)
        sess.deliver[] = previous_deliver
        sess.spawn[] = previous_spawn
        # For a target left building, the returned session holds the target of the
        # design, and parts left remeshing are meshed when the drive requests them.
        ts = sess.target
        ts.building[] && isempty(ts.pending[]) && await_target!(ts)
        _release_qml_screens!()
    end
    return sess
end
