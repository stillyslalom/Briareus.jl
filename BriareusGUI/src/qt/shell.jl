# Qt Quick shell for the Briareus window (QML.jl + QMLMakie).
#
# Threading: Qt's event loop owns the main thread while the window is open. QML→Julia
# calls only change state and return. Background jobs hand their observable updates to
# `deliver`, which queues them; `br_tick` (a QML Timer, ~60 Hz) drains the queue on the
# GUI thread and refreshes what QML shows.
#
# Canvas: the QMLMakie canvas (`MakieArea`) is never destroyed while the window is
# open — destroying one leaves a dangling callback in jlqml that crashes at window
# teardown. The 3D viewer and the flux map are two figures shown in turn on it.

const QML_DIR = joinpath(@__DIR__, "..", "qml")

# ---------------------------------------------------------------- Qt setup

const _QT_READY = Ref(false)

# GLFW contexts (GLMakie windows and offscreen `colorbuffer` screens) and Qt's GL
# contexts must not share a process: once a GLFW context had existed, Qt's render
# thread crashed in the AMD driver (see HammerheadGUI). The window refuses to open next
# to a GLFW screen.
_glfw_screen_count() = count(s -> !(s isa GLMakie.Screen{QMLMakie.QMLWindow}), GLMakie.ALL_SCREENS)

function _check_no_glfw_screens()
    n = _glfw_screen_count()
    n == 0 && return
    throw(ArgumentError(
        "this Julia session has created $n GLMakie screen" * (n == 1 ? "" : "s") *
        " (a GLMakie window or an offscreen render), and GLMakie windows and the Qt " *
        "window cannot share a process (the graphics driver can crash). Restart Julia " *
        "to open the Briareus window."))
end

function _qt_init!()
    _QT_READY[] && return
    if Sys.iswindows()
        # The FluentWinUI3 style plugin depends on a library outside the DLL search
        # path; load it first so the native Windows 11 style is used.
        impl = joinpath(Qt6Declarative_jll.artifact_dir, "bin",
                        "Qt6QuickControls2FluentWinUI3StyleImpl.dll")
        isfile(impl) && Libdl.dlopen(impl)
        # Qt reads msvcrt's copy of the environment, which Julia's ENV (Win32
        # SetEnvironmentVariable) does not update, so set both.
        style = get(ENV, "QT_QUICK_CONTROLS_STYLE", "FluentWinUI3")
        ENV["QT_QUICK_CONTROLS_STYLE"] = style
        ccall((:_putenv_s, "msvcrt"), Cint, (Cstring, Cstring), "QT_QUICK_CONTROLS_STYLE", style)
    end
    _QT_READY[] = true
    return
end

# Glyphs a canvas may show (ASCII, the symbols of the labels, Makie's tick minus). A
# glyph first laid out while the window is open reaches the atlas texture through a
# callback that runs outside Qt's render, where no GL context is current, and draws
# as garbage; glyphs in the atlas before the first render upload with it.
const _CANVAS_GLYPHS = vcat(Char(32):Char(126), collect("°±²³·×θφεΔσ–—−…→≤≥≈"))

# Rendering these glyphs takes seconds, so the warmed atlas is cached on disk next to
# Makie's own atlas cache and loaded on the next start, provided no screen uses the
# session's atlas yet and the cache holds every glyph it has.
const _ATLAS_KEY = (2048, 64)                        # GLMakie's atlas
_atlas_cache_path() = Makie.get_cache_path(_ATLAS_KEY...) * ".briareusgui"

function _load_cached_atlas!(atlas)
    path = _atlas_cache_path()
    (isempty(atlas.font_render_callback) && isfile(path)) || return atlas
    cached = try
        Makie.load_texture_atlas(path)
    catch
        return atlas
    end
    all(k -> haskey(cached.mapping, k), keys(atlas.mapping)) || return atlas
    Makie.TEXTURE_ATLASES[_ATLAS_KEY] = cached
    return cached
end

function _warm_glyph_atlas!()
    atlas = _load_cached_atlas!(Makie.get_texture_atlas(_ATLAS_KEY...))
    n = length(atlas.mapping)
    fonts = Makie.theme(:fonts)
    for key in keys(fonts)
        font = Makie.to_font(fonts[key][])
        for c in _CANVAS_GLYPHS
            Makie.insert_glyph!(atlas, c, font)
        end
    end
    if length(atlas.mapping) > n
        path = _atlas_cache_path()
        try
            tmp = path * ".$(getpid()).tmp"
            Makie.store_texture_atlas(tmp, atlas)
            mv(tmp, path; force = true)
        catch
        end
    end
    return
end

# After the window closes Qt has destroyed its GL contexts: forget the canvas's
# screens so GLMakie neither renders to them nor frees them at exit.
function _release_qml_screens!()
    for s in GLMakie.ALL_SCREENS
        s isa GLMakie.Screen{QMLMakie.QMLWindow} && (s.glscreen.context.valid = false)
    end
    filter!(s -> !(s isa GLMakie.Screen{QMLMakie.QMLWindow}), GLMakie.ALL_SCREENS)
    for key in collect(keys(GLMakie.atlas_texture_cache))
        key[2] isa QMLMakie.QMLWindow || continue
        _, callback = pop!(GLMakie.atlas_texture_cache, key)
        Makie.remove_font_render_callback!(key[1], callback)
    end
    return
end

# QML FileDialog URLs → local paths.
function _url_to_path(url::AbstractString)
    s = String(url)
    startswith(s, "file://") || return s
    s = s[8:end]                                   # drop "file://"
    Sys.iswindows() && startswith(s, "/") && occursin(r"^/[A-Za-z]:", s) && (s = s[2:end])
    return _percent_decode(s)
end

function _percent_decode(s::String)
    io = IOBuffer()
    i = firstindex(s)
    while i <= lastindex(s)
        if s[i] == '%' && i + 2 <= lastindex(s)
            write(io, parse(UInt8, s[i+1:i+2]; base = 16))
            i += 3
        else
            write(io, s[i])
            i = nextind(s, i)
        end
    end
    return String(take!(io))
end

# ---------------------------------------------------------------- bridge rows

mutable struct StepRow
    key::String
    label::String
    status::String                                 # :todo/:ok/:attention/:busy ("state" clashes with QML)
    summary::String
end

# One beam-table row; numbers travel as text-field text. `on`: 0 off, 1 on, 2 mixed.
mutable struct BeamRow
    number::Int                                    # "index" is reserved in QML delegates
    key::String
    count::Int
    theta::String
    phi::String
    z::String
    r::String
    phiOffset::String
    defocus::String
    power::String
    on::Int
end

mutable struct ModeRow
    l::Int
    epsilon::String
    legendre::String
    before::String                                 # Optimize: ε_l before → after
    after::String
end

mutable struct PowerRow
    name::String
    power::String
    percent::String
end

mutable struct ParamRow
    name::String
    start::String
    optimized::String
    lower::String
    upper::String
    atBound::Bool
end

mutable struct CompRow
    name::String
    role::String                                   # "target", "hardware", "none"
end

mutable struct ClearRow
    beam::String
    component::String
    status::String
    worst::String
    nominal::String
    power::String
    intensity::String
end

_same(a, b) = all(f -> getfield(a, f) == getfield(b, f), fieldnames(typeof(a)))

# Replace the rows of an item model: in place when the count is unchanged (only
# changed rows are re-sent), otherwise with a full reset.
function _sync_rows!(rows::Vector{R}, model::JuliaItemModel, new::Vector{R}) where {R}
    if length(new) == length(rows)
        for (i, r) in enumerate(new)
            _same(rows[i], r) && continue
            rows[i] = r
            model[i] = r
        end
    else
        empty!(rows)
        append!(rows, new)
        QML.force_model_update(model)
    end
    return
end

# ---------------------------------------------------------------- shell

const _MODEL_ROWS = (stepModel = StepRow, beamModel = BeamRow, modeModel = ModeRow,
                     powerModel = PowerRow, optModeModel = ModeRow, paramModel = ParamRow,
                     compModel = CompRow, clearModel = ClearRow, worstModel = ClearRow)

"""
State shared between a `DesignSession` and its QML window: the property map `app`,
the item models (by their QML names) and their rows, the 3D and map canvases, and
the queue of background updates that `br_tick` drains.
"""
mutable struct BriareusShell
    sess::DesignSession
    scene::SceneCanvas
    map::MapCanvas
    app::JuliaPropertyMap
    rows::Dict{Symbol,Vector}
    models::Dict{Symbol,JuliaItemModel}
    queue::Channel{Any}
    dirty::Base.RefValue{Bool}
    shown::Dict{String,Any}
    figure::Any                                    # the figure on the canvas
end

function BriareusShell(sess::DesignSession; queue::Channel{Any} = Channel{Any}(Inf))
    rows = Dict{Symbol,Vector}(k => R[] for (k, R) in pairs(_MODEL_ROWS))
    models = Dict{Symbol,JuliaItemModel}(k => JuliaItemModel(rows[k]) for k in keys(_MODEL_ROWS))
    scene = scene_canvas(sess)
    sh = BriareusShell(sess, scene, map_canvas(sess), JuliaPropertyMap(), rows, models, queue,
                       Ref(true), Dict{String,Any}(), scene.fig)
    sh.app["canvas"] = scene.fig
    _set!(sh, "grabPath", "")
    _set!(sh, "cutaway", true)
    mark = (_...) -> (sh.dirty[] = true)
    ts, bt, dr, st, rv = sess.target, sess.beams, sess.drive, sess.optimize, sess.review
    for obs in (sess.step, sess.status, ts.kind, ts.params, ts.resolution, ts.target, ts.error,
                bt.beams, bt.enabled, bt.preset, bt.by, bt.path, bt.error,
                dr.scenario, dr.building, dr.result, dr.evaluating, dr.field, dr.view,
                dr.spot_radius, dr.status,
                st.by, st.fields, st.deltas, st.objective, st.weights, st.hold_drive, st.keepout,
                st.keepout_min, st.spill, st.iterations, st.running, st.evaluations, st.status,
                st.error, st.result,
                rv.targets, rv.hardware, rv.peak_fraction, rv.wings, rv.errors, rv.quality,
                rv.running, rv.progress, rv.status, rv.error, rv.report)
        on(mark, obs)
    end
    _refresh!(sh)
    return sh
end

# Assign a QML property only when its value changed.
function _set!(sh::BriareusShell, key::String, val)
    haskey(sh.shown, key) && isequal(sh.shown[key], val) && return
    sh.shown[key] = val
    sh.app[key] = val
    return
end

_num(x::Real) = isfinite(x) ? Controllers.display_number(x) : ""
_fix(x::Real, d::Integer = 3) = isfinite(x) ? @sprintf("%.*f", d, x) : "–"
_pct(x::Real, d::Integer = 2) = isfinite(x) ? @sprintf("%.*f %%", d, 100x) : "–"
_g(x::Real) = isfinite(x) ? @sprintf("%.3g", x) : "–"

function _refresh!(sh::BriareusShell)
    sess = sh.sess
    _set!(sh, "title", "Briareus | " * target_summary(sess.target))
    _set!(sh, "step", String(sess.step[]))
    _set!(sh, "status", sess.status[])
    _refresh_target!(sh)
    _refresh_beams!(sh)
    _refresh_drive!(sh)
    _refresh_optimize!(sh)
    _refresh_review!(sh)
    steps = [StepRow(String(s), step_label(sess, s), String(st), txt)
             for s in DESIGN_STEPS for (st, txt) in (step_status(sess, s),)]
    _sync_rows!(sh.rows[:stepModel], sh.models[:stepModel], steps)
    # the flux map on Drive/Optimize when chosen there; the 3D view otherwise
    fig = sess.step[] in (:drive, :optimize) && sess.drive.view[] === :map ? sh.map.fig : sh.scene.fig
    if fig !== sh.figure
        sh.figure = fig
        sh.app["canvas"] = fig
    end
    return
end

function _refresh_target!(sh::BriareusShell)
    ts = sh.sess.target
    _set!(sh, "targetKind", String(ts.kind[]))
    _set!(sh, "targetResolution", String(ts.resolution[]))
    _set!(sh, "targetSummary", target_summary(ts))
    _set!(sh, "targetError", ts.error[])
    for (k, v) in ts.params[]
        _set!(sh, "t_" * String(k), _num(v))
    end
    return
end

function _refresh_beams!(sh::BriareusShell)
    bt = sh.sess.beams
    _set!(sh, "preset", String(bt.preset[]))
    _set!(sh, "presetLabel", get(PRESET_LABELS, bt.preset[], String(bt.preset[])))
    _set!(sh, "beamsBy", String(bt.by[]))
    _set!(sh, "beamsSummary", beams_summary(bt))
    _set!(sh, "beamsError", bt.error[])
    _set!(sh, "beamsPath", bt.path[])
    rows = [BeamRow(i, r.key, r.count, @sprintf("%.2f", r.theta), isnan(r.phi) ? "–" : @sprintf("%.2f", r.phi),
                    _fix(r.z), _fix(r.r), _fix(r.phi_offset), _fix(r.defocus), _num(r.power),
                    r.enabled === missing ? 2 : r.enabled ? 1 : 0)
            for (i, r) in enumerate(beam_rows(bt))]
    _sync_rows!(sh.rows[:beamModel], sh.models[:beamModel], rows)
    return
end

function _refresh_drive!(sh::BriareusShell)
    sess = sh.sess
    dr = sess.drive
    res = dr.result[]
    names = component_names(sess.target.target[])
    _set!(sh, "hasCapsule", :capsule in names)
    _set!(sh, "hasFoil", FOIL_NAME in names)
    _set!(sh, "driveBuilding", dr.building[])
    _set!(sh, "driveElapsed", dr.building[] ? Controllers._seconds(time() - dr.started[]) : "")
    _set!(sh, "driveEvaluating", dr.evaluating[])
    _set!(sh, "driveBuilt", dr.scenario[] !== nothing)
    _set!(sh, "driveStale", drive_stale(sess))
    _set!(sh, "driveStatus", dr.status[])
    _set!(sh, "driveHasResult", res !== nothing)
    _set!(sh, "driveField", String(dr.field[]))
    _set!(sh, "driveView", String(dr.view[]))
    _set!(sh, "spotRadius", _num(dr.spot_radius[]))
    _set!(sh, "driveSummary", res === nothing ? "" : drive_summary(res))
    _set!(sh, "driveMean", res === nothing || isnan(res.mean_flux) ? "" :
                           @sprintf("%.3g W/cm² · Tr %.1f eV", to_Wcm2(res.mean_flux),
                                    radiation_temperature(res.mean_flux)))
    s = res === nothing ? nothing : res.spot
    _set!(sh, "spotText", s === nothing ? "" :
                          @sprintf("mean %.3g W/cm², peak %.3g W/cm², rms nonuniformity %.1f %%, %.1f %% of the delivered power",
                                   s.mean, s.peak, 100 * s.rms, 100 * s.fraction))
    modes = res === nothing ? ModeRow[] :
            [ModeRow(m.l, _pct(m.epsilon), _pct(m.legendre), "", "") for m in res.modes]
    _sync_rows!(sh.rows[:modeModel], sh.models[:modeModel], modes)
    power = res === nothing ? PowerRow[] :
            [PowerRow(p.name, @sprintf("%.4g", p.power), _pct(p.fraction)) for p in res.power]
    _sync_rows!(sh.rows[:powerModel], sh.models[:powerModel], power)
    _set!(sh, "delivered", res === nothing ? "" :
                           @sprintf("%.4g TW delivered", res.evaluation.deposition.delivered))
    return
end

function _refresh_optimize!(sh::BriareusShell)
    st = sh.sess.optimize
    _set!(sh, "optBy", String(st.by[]))
    for f in BEAM_FIELDS
        _set!(sh, "optFree_" * String(f), f in st.fields[])
        _set!(sh, "optDelta_" * String(f), _num(st.deltas[][f]))
    end
    _set!(sh, "optObjective", String(st.objective[]))
    for l in 1:4
        _set!(sh, "optW$l", _num(get(st.weights[], l, 0.0)))
    end
    _set!(sh, "optHold", st.hold_drive[])
    _set!(sh, "optKeepout", st.keepout[])
    _set!(sh, "optKeepoutMin", _num(st.keepout_min[]))
    _set!(sh, "optSpill", _num(st.spill[]))
    _set!(sh, "optIterations", st.iterations[])
    _set!(sh, "optRunning", st.running[])
    _set!(sh, "optProgress", optimize_progress(st))
    _set!(sh, "optStatus", st.status[])
    _set!(sh, "optError", st.error[])
    r = st.result[]
    _set!(sh, "optHasResult", r !== nothing)
    _set!(sh, "optApplied", r !== nothing && r.applied)
    _set!(sh, "optValues", r === nothing ? "" :
                           @sprintf("objective %.4g → %.4g · %d evaluations · %s%s", r.initial_value,
                                    r.value, r.evaluations, Controllers._seconds(r.seconds),
                                    r.completed ? "" : " · canceled (best so far)"))
    _set!(sh, "optClipped", r === nothing ? "" :
                            "clipped laser power " * _pct(r.before.clipped) * " → " * _pct(r.after.clipped))
    _set!(sh, "optSpot", r === nothing || r.after.spot === nothing ? "" :
                         @sprintf("spot rms nonuniformity %.1f %% → %.1f %%", 100 * r.before.spot.rms,
                                  100 * r.after.spot.rms))
    modes = r === nothing || isempty(r.after.modes) ? ModeRow[] :
            [ModeRow(a.l, _pct(a.epsilon), _pct(a.legendre), _pct(b.epsilon), _pct(a.epsilon))
             for (b, a) in zip(r.before.modes, r.after.modes)]
    _sync_rows!(sh.rows[:optModeModel], sh.models[:optModeModel], modes)
    params = r === nothing ? ParamRow[] :
             [ParamRow(p.name, _fix(p.start), _fix(p.optimized), _fix(p.lower), _fix(p.upper), p.at_bound)
              for p in parameter_rows(r)]
    _sync_rows!(sh.rows[:paramModel], sh.models[:paramModel], params)
    return
end

function _clear_row(x)
    return ClearRow(x.beam, String(x.component), String(x.status), _fix(x.clearance),
                    _fix(x.nominal), @sprintf("%.3g", x.power), @sprintf("%.3g", x.peak_intensity))
end

function _refresh_review!(sh::BriareusShell)
    sess = sh.sess
    rv = sess.review
    e = rv.errors[]
    _set!(sh, "rvPeak", @sprintf("%g", rv.peak_fraction[]))
    _set!(sh, "rvWings", rv.wings[])
    for k in (:pointing, :target_offset, :hardware_tolerance, :focus)
        _set!(sh, "rv_" * String(k), _num(e[k]))
    end
    _set!(sh, "rvQuality", String(rv.quality[]))
    _set!(sh, "rvRunning", rv.running[])
    _set!(sh, "rvDone", rv.progress[][1])
    _set!(sh, "rvTotal", rv.progress[][2])
    _set!(sh, "rvStatus", rv.status[])
    _set!(sh, "rvError", rv.error[])
    _set!(sh, "rvHasReport", rv.report[] !== nothing)
    _set!(sh, "rvBrief", review_brief(rv))
    _set!(sh, "rvStale", review_stale(sess))
    comps = [CompRow(String(n), n in rv.targets[] ? "target" : n in rv.hardware[] ? "hardware" : "none")
             for n in surface_names(sess.target.target[])]
    _sync_rows!(sh.rows[:compModel], sh.models[:compModel], comps)
    _sync_rows!(sh.rows[:clearModel], sh.models[:clearModel], ClearRow[_clear_row(x) for x in clearance_rows(rv)])
    _sync_rows!(sh.rows[:worstModel], sh.models[:worstModel], ClearRow[_clear_row(x) for x in worst_rows(rv)])
    return
end

# ---------------------------------------------------------------- QML callbacks

const _SHELL = Ref{Union{Nothing,BriareusShell}}(nothing)
# `_TICK_HOOK[](shell)` runs on every tick (scripted smoke tests drive the window
# through it); `request_close()` closes the window from Julia.
const _TICK_HOOK = Ref{Any}(nothing)
const _CLOSE_REQUESTED = Ref(false)

"""
    request_close()

Close the open Briareus window (from a callback or a tick hook).
"""
request_close() = (_CLOSE_REQUESTED[] = true; nothing)

"""
    request_grab(path)

Save an image of the open window's body (step page and viewer) to `path` (PNG) on the
window's next tick, through QML `grabToImage`. The file appears asynchronously.
"""
request_grab(path::AbstractString) = (sh = _SHELL[]; sh === nothing || _set!(sh, "grabPath", String(path)); nothing)
br_grab_started() = (sh = _SHELL[]; sh === nothing || _set!(sh, "grabPath", ""); nothing)

# Run a callback body on the current shell, reporting errors in the status line
# instead of letting them escape into Qt.
function _with_shell(f)
    sh = _SHELL[]
    sh === nothing && return nothing
    try
        f(sh)
    catch err
        sh.sess.status[] = "error: " * Controllers._errmsg(err)
        @error "BriareusGUI callback failed" exception = (err, catch_backtrace())
    end
    sh.dirty[] = true
    return nothing
end

# Returns 0 (nothing to draw), 1 (redraw the canvas), or 2 (close the window).
function br_tick()
    sh = _SHELL[]
    sh === nothing && return 0
    # Only thread 1 runs libuv's event loop, and it sits in Qt's exec(): pump it here.
    Base.process_events()
    if _TICK_HOOK[] !== nothing
        try
            _TICK_HOOK[](sh)
        catch err
            @error "tick hook failed" exception = (err, catch_backtrace())
            request_close()
        end
    end
    _CLOSE_REQUESTED[] && (_CLOSE_REQUESTED[] = false; return 2)
    while isready(sh.queue)
        f = take!(sh.queue)
        try
            f()
        catch err
            @error "BriareusGUI background update failed" exception = (err, catch_backtrace())
        end
    end
    sess = sh.sess
    # running jobs report elapsed time on the rail
    (sess.drive.building[] || sess.optimize.running[] || sess.review.running[]) && (sh.dirty[] = true)
    sh.dirty[] && (sh.dirty[] = false; _refresh!(sh))
    redraw = sh.scene.dirty[] || sh.map.dirty[]
    sh.scene.dirty[] = false
    sh.map.dirty[] = false
    return redraw ? 1 : 0
end

_sym(x) = Symbol(String(x))

br_set_step(name) = _with_shell(sh -> set_step!(sh.sess, _sym(name)))
# Target
br_set_target_kind(kind) = _with_shell(sh -> set_target_kind!(sh.sess.target, _sym(kind)))
br_set_resolution(r) = _with_shell(sh -> set_resolution!(sh.sess.target, _sym(r)))
br_edit_target(key, text) = _with_shell(sh -> edit_target!(sh.sess.target, _sym(key), string(text)))
# Beams
br_load_preset(name) = _with_shell(sh -> load_preset!(sh.sess, _sym(name)))
br_set_grouping(by) = _with_shell(sh -> set_grouping!(sh.sess.beams, _sym(by)))
function br_edit_beam(row, field, text)
    _with_shell() do sh
        k = round(Int, row)
        # a rejected edit leaves the beams unchanged: re-send the row so the field
        # shows the kept value again
        edit_beams!(sh.sess.beams, k, _sym(field), string(text)) ||
            (k <= length(sh.rows[:beamModel]) && (sh.models[:beamModel][k] = sh.rows[:beamModel][k]))
    end
end
br_set_beam_on(row, on) = _with_shell(sh -> set_enabled!(sh.sess.beams, round(Int, row), Bool(on)))
br_import_beams(url) = _with_shell(sh -> (import_beams!(sh.sess.beams, _url_to_path(String(url)));
                                          sh.sess.status[] = "imported " * basename(sh.sess.beams.path[])))
br_export_beams(url) = _with_shell(sh -> (export_beams(sh.sess.beams, _url_to_path(String(url)));
                                          sh.sess.status[] = "exported " * basename(sh.sess.beams.path[])))
# Drive
br_build() = _with_shell(sh -> build_scenario!(sh.sess))
br_set_field(f) = _with_shell(sh -> set_field!(sh.sess.drive, _sym(f)))
br_set_view(v) = _with_shell(sh -> set_view!(sh.sess.drive, _sym(v)))
br_set_spot_radius(text) = _with_shell(sh -> (sh.sess.drive.spot_radius[] =
                                                  Controllers._parse_number(string(text), "the spot radius"; positive = true)))
# Optimize
br_opt_set_by(by) = _with_shell(sh -> (by = _sym(by); by in BEAM_GROUPINGS &&
                                                      Controllers._assign!(sh.sess.optimize.by, by)))
br_opt_free(field, on) = _with_shell(sh -> set_free!(sh.sess.optimize, _sym(field), Bool(on)))
br_opt_edit(key, text) = _with_shell(sh -> edit_optimize!(sh.sess.optimize, _sym(key), string(text)))
br_opt_objective(k) = _with_shell(sh -> set_objective!(sh.sess.optimize, _sym(k)))
br_opt_hold(on) = _with_shell(sh -> (sh.sess.optimize.hold_drive[] = Bool(on)))
br_opt_keepout(on) = _with_shell(sh -> (sh.sess.optimize.keepout[] = Bool(on)))
br_opt_start() = _with_shell(sh -> start_optimize!(sh.sess))
br_opt_cancel() = _with_shell(sh -> cancel_optimize!(sh.sess))
br_opt_accept() = _with_shell(sh -> accept_optimum!(sh.sess))
# Review
br_rv_role(name, role) = _with_shell(sh -> set_role!(sh.sess.review, _sym(name), _sym(role)))
br_rv_edit(key, text) = _with_shell(sh -> edit_review!(sh.sess.review, _sym(key), string(text)))
br_rv_wings(on) = _with_shell(sh -> (sh.sess.review.wings[] = Bool(on)))
br_rv_quality(q) = _with_shell(sh -> (q = _sym(q); haskey(REVIEW_QUALITY, q) &&
                                                   Controllers._assign!(sh.sess.review.quality, q)))
br_rv_start() = _with_shell(sh -> start_clearance!(sh.sess))
br_rv_cancel() = _with_shell(sh -> cancel_clearance!(sh.sess))
br_rv_export(url) = _with_shell(sh -> export_clearance(sh.sess.review, _url_to_path(String(url))))
# Viewer
br_set_cutaway(on) = _with_shell(sh -> (set_cutaway!(sh.scene, sh.sess, Bool(on)); _set!(sh, "cutaway", Bool(on))))
br_reset_view() = _with_shell(sh -> sh.figure === sh.map.fig ? (reset_limits!(sh.map.ax); sh.map.dirty[] = true) :
                                    reset_view!(sh.scene))

function _register_qml_functions()
    @qmlfunction br_tick br_grab_started br_set_step br_set_target_kind br_set_resolution
    @qmlfunction br_edit_target br_load_preset br_set_grouping br_edit_beam br_set_beam_on
    @qmlfunction br_import_beams br_export_beams br_build br_set_field br_set_view
    @qmlfunction br_set_spot_radius br_opt_set_by br_opt_free br_opt_edit br_opt_objective
    @qmlfunction br_opt_hold br_opt_keepout br_opt_start br_opt_cancel br_opt_accept
    @qmlfunction br_rv_role br_rv_edit br_rv_wings br_rv_quality br_rv_start br_rv_cancel
    @qmlfunction br_rv_export br_set_cutaway br_reset_view
    return
end

# ---------------------------------------------------------------- launcher

"""
    briareus(sess = DesignSession(); kind = :hohlraum, resolution = :standard,
             preset = nothing, beams = nothing, build = false) -> DesignSession

Open the Briareus window and return its session when the window closes. Its steps:

- **Target** — a hohlraum (length, radius, LEH and capsule radii, wall albedo and
  conversion) or a backlighter foil, and the mesh resolution; the 3D view colors the
  mesh by component.
- **Beams** — a facility preset or an imported beam table (CSV), edited per cone, quad
  or beam: aim, defocus, power, on/off; the 3D view draws each beam's axis from its
  port to its aim point.
- **Drive** — the view factors (computed on a worker; seconds to tens of seconds),
  then the target colored by absorbed laser intensity, x-ray radiosity or incident
  flux (W/cm²), the capsule flux map over θ–φ, the mode amplitudes ε₁…ε₄ and Pₗ/P₀, and
  the laser power accounting. Later beam edits re-evaluate at once.
- **Optimize** — free fields and bounds per group, an objective (capsule symmetry,
  axisymmetric or 3D, or foil uniformity) with mode weights, run in the background
  with progress and cancel; before/after modes and the parameter table; accepting
  puts the optimized pointing into the beam table.
- **Review** — worst-case beam-to-hardware clearance on chosen hardware, in the
  background with progress and cancel; rows colored by status, the closest beam per
  component, and CSV export.

The session starts from `sess`, or a fresh `DesignSession(; kind, resolution,
preset, beams)`; `build = true` starts computing the drive at once. The call blocks
while the window is open (Qt runs its event loop on this thread); start Julia with
several threads (`julia -t auto`) so background jobs leave the window responsive.
Closing the window cancels an optimization or clearance check in progress.

The window throws an `ArgumentError` when this Julia session already has a GLMakie
screen (a GLMakie window or an offscreen `colorbuffer`): GLMakie's windows and the Qt
window cannot share a process.
"""
function briareus(sess::Union{Nothing,DesignSession} = nothing; kind::Symbol = :hohlraum,
                  resolution::Symbol = :standard, preset = nothing, beams = nothing,
                  build::Bool = false)
    if sess === nothing
        sess = preset === nothing ? DesignSession(; kind, resolution, beams) :
               DesignSession(; kind, resolution, preset, beams)
    end
    _SHELL[] === nothing || throw(ArgumentError("a BriareusGUI window is already open"))
    _check_no_glfw_screens()
    # From here on background work goes to the window's queue.
    queue = Channel{Any}(Inf)
    previous_deliver, previous_spawn = sess.deliver[], sess.spawn[]
    sess.deliver[] = f -> put!(queue, f)
    sess.spawn[] = true
    try
        _qt_init!()
        _warm_glyph_atlas!()
        sh = BriareusShell(sess; queue)
        _SHELL[] = sh
        _CLOSE_REQUESTED[] = false
        _register_qml_functions()
        build && build_scenario!(sess)
        loadqml(joinpath(QML_DIR, "BriareusWindow.qml"); app = sh.app,
                (k => m for (k, m) in sh.models)...)
        exec()
    finally
        _SHELL[] = nothing
        Controllers._abandon_jobs!(sess)
        sess.deliver[] = previous_deliver
        sess.spawn[] = previous_spawn
        _release_qml_screens!()
    end
    return sess
end
