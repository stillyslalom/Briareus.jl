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
# The 3D viewer's mouse and keyboard gestures reach Makie's events through it
# (canvas/cad_camera.jl).

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
        ccall((:SetCurrentProcessExplicitAppUserModelID, "shell32"), Clong, (Cwstring,),
              APP_USER_MODEL_ID)
        _ENUM_QT_WINDOWS[] = _enum_qt_windows_ptr()
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

# Compile what the canvases first do with a drive or a CAD import before the window
# opens, so the first drive result, the first coloring by a field and the first parts
# left out of the calculations do not stall the window for most of a second each (JIT
# compilation, measured on a NIF CAD workspace; the package's precompile workload does
# not keep these). CPU side: a filled flux map and a field colorbar on throwaway
# canvases (their own session). GL side, at the window's first frames: the transparent
# ghost plot is drawn from a degenerate triangle (invisible), and the colorbar is shown
# as a drive field shows it, given a second scale (as the next field coloring changes
# it: the colormap and gradient textures' updates) and hidden again.
function _warm_canvases!(sh)
    ws = DesignSession(resolution = :coarse)
    mc = map_canvas(ws)
    x, y = collect(1.0:4.0), collect(1.0:3.0)
    _update!(mc.heatmap, x, y, Float32[k + j for k in 1:4, j in 1:3]; colorrange = (1.0, 7.0))
    mc.ax.aspect[] = AxisAspect(2)
    reset_limits!(mc.ax)
    sc = scene_canvas(ws)
    _field_colorbar!(sc.colorbar, 1.3e14)
    _show_colorbar!(sc, true)
    # the window's own canvas
    c = sh.scene
    sh.sess.target.ghost[] === nothing && _update!(c.ghost, _NOTRI...; visible = true)
    if !c.colorbar.blockscene.visible[]
        _field_colorbar!(c.colorbar, 1.3e14)
        _show_colorbar!(c, true)
        ticks = Ref(0)
        once = Ref{Any}(nothing)
        once[] = on(events(c.fig.scene).tick) do _
            # a tick precedes each frame's drawing: a second scale at the second, hidden
            # at the third, unless a drive field colors the body by then
            k = (ticks[] += 1)
            k < 2 && return
            field = _body_coloring(sh.sess) !== nothing
            if k == 2
                field || (_field_colorbar!(c.colorbar, 2.6e14); c.dirty[] = true)
                return
            end
            off(once[])
            field || (_show_colorbar!(c, false); c.dirty[] = true)
        end
    end
    return nothing
end

# A colorbar set as a drive field sets it (`_draw_body!`).
function _field_colorbar!(cb, hi)
    cb.colormap[] = FIELD_COLORMAP
    cb.limits[] = (0.0, hi)
    cb.ticks[] = Makie.automatic
    cb.label[] = Controllers.DRIVE_FIELD_LABELS[:incident]
    return cb
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

# One visible row of the model tree (see `Controllers.model_tree`).
mutable struct TreeItem
    id::String
    depth::Int
    label::String
    detail::String
    kind::String
    expandable::Bool
    expanded::Bool
    on::Int                                        # -1 none, 0 off, 1 on, 2 mixed
    status::String
    selected::Bool
    color::String                                  # swatch (hex), "" for none
    drawmode::String                               # a component's display mode, "" for none (not
                                                   # `display`: Qt's DisplayRole name)
    matrix::String                                 # a beam row's 2 × 3 dots, "" for none
end

# An editable field of the inspector: key, label, unit, value as text.
mutable struct FieldRow
    key::String
    label::String
    unit::String
    value::String
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
    before::String                                 # optimization: ε_l before → after
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
    name::String                                   # the component's name, for the callbacks
    label::String                                  # as the model tree names it
    role::String                                   # "target"/"capture", "hardware", "diagnostic", "none"
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

# One diagnostic window; numbers as text-field text.
mutable struct WindowRow
    number::Int
    name::String
    phi::String
    z::String
    width::String
    height::String
    open::Bool
    albedo::String
    conversion::String
    round::Bool
end

# TaLIS: one rule's summary, one check, one note.
mutable struct RuleRow
    label::String
    requirement::String
    pass::Int
    review::Int
    fail::Int
    worst::String                                  # worst margin with unit
    where_::String                                 # its beam / component
    status::String                                 # worst status
end

mutable struct TalisRow
    beam::String
    rule::String
    component::String
    status::String
    margin::String
    worst::String
    nominal::String
    required::String
    note::String
end

mutable struct NoteRow
    kind::String                                   # tolerance, settings, skipped, scope, source
    text::String
end

# One step of the placement sweep; numbers as text.
mutable struct SweepRowItem
    number::Int
    value::String
    clearance::String
    worst::String
    fails::Int
    reviews::Int
    los::String
    blocked::String
    covered::String
    captured::String
    objective::String
end

# One entry of the time chart's key.
mutable struct LegendRow
    label::String
    color::String
    style::String                                  # solid, dash, dot
end

_same(a, b) = all(f -> getfield(a, f) == getfield(b, f), fieldnames(typeof(a)))

# How many leading rows of `old` and `new` are the same nodes (equal first fields): rows
# are inserted or removed after them.
function _common_lead(old::AbstractVector, new::AbstractVector)
    p = 0
    while p < min(length(old), length(new)) && getfield(old[p+1], 1) == getfield(new[p+1], 1)
        p += 1
    end
    return p
end

# Replace the rows of an item model: only changed rows are re-sent, and a change in the
# count inserts or removes rows where the old and new rows part (aligned on their first
# field: a tree node's id), so the view keeps the delegates of the others. A full reset
# rebuilt every delegate of the model tree: 0.2 s on a CAD workspace each time a node
# was expanded or a part added. `new` may be an empty comprehension whose eltype
# inference gave up on (`Vector{Union{}}`).
function _sync_rows!(rows::Vector{R}, model::JuliaItemModel, new::AbstractVector) where {R}
    new = convert(Vector{R}, new)
    n, m = length(rows), length(new)
    if n != m && n > 0 && m > 0
        p = _common_lead(rows, new)
        if m > n
            QML.begin_insert_rows(model, p + 1, p + m - n)
            for k in 1:m-n
                insert!(rows, p + k, new[p+k])
            end
            QML.@deferredcall QML.end_insert_rows(model)
        else
            QML.remove_rows!(model, p + 1, n - m)
        end
    end
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

const _MODEL_ROWS = (treeModel = TreeItem, targetFieldModel = FieldRow, partFieldModel = FieldRow,
                     windowModel = WindowRow, beamModel = BeamRow, modeModel = ModeRow,
                     powerModel = PowerRow, optModeModel = ModeRow, paramModel = ParamRow,
                     compModel = CompRow, clearModel = ClearRow, worstModel = ClearRow,
                     talisCompModel = CompRow, ruleModel = RuleRow, talisModel = TalisRow,
                     noteModel = NoteRow, timeLegendModel = LegendRow, sweepModel = SweepRowItem,
                     sweepLegendModel = LegendRow)

"""
State shared between a `DesignSession` and its QML window: the property map `app`,
the item models (by their QML names) and their rows, the 3D and map canvases, and
the queue of background updates that `br_tick` drains.
"""
mutable struct BriareusShell
    sess::DesignSession
    scene::SceneCanvas
    map::MapCanvas
    time::TimeCanvas
    sweep::TimeCanvas
    app::JuliaPropertyMap
    rows::Dict{Symbol,Vector}
    models::Dict{Symbol,JuliaItemModel}
    queue::Channel{Any}
    dirty::Base.RefValue{Bool}
    shown::Dict{String,Any}
    figure::Any                                    # the figure on the canvas
    session_path::String                           # the session file, once saved or opened
    side_figure::Any                               # the figure on the side canvas
    blank::Any                                     # what the side canvas shows when hidden
end

function BriareusShell(sess::DesignSession; queue::Channel{Any} = Channel{Any}(Inf))
    rows = Dict{Symbol,Vector}(k => R[] for (k, R) in pairs(_MODEL_ROWS))
    models = Dict{Symbol,JuliaItemModel}(k => JuliaItemModel(rows[k]) for k in keys(_MODEL_ROWS))
    scene = scene_canvas(sess)
    sh = BriareusShell(sess, scene, map_canvas(sess), time_canvas(sess), sweep_canvas(sess),
                       JuliaPropertyMap(), rows, models, queue, Ref(true), Dict{String,Any}(),
                       scene.fig, "", nothing, Figure(size = (40, 40)))
    sh.app["canvas"] = scene.fig
    sh.app["timeCanvas"] = sh.time.fig
    sh.app["sweepCanvas"] = sh.sweep.fig
    # a figure is shown on one canvas at a time: the side canvas holds a blank one
    # unless the map sits beside the 3D view
    sh.side_figure = sh.blank
    sh.app["sideCanvas"] = sh.blank
    _set!(sh, "sideView", false)
    _set!(sh, "grabPath", "")
    _set!(sh, "resultsRequest", "")
    _set!(sh, "sessionPath", "")
    _set!(sh, "canUndo", false)
    _set!(sh, "canRedo", false)
    _set!(sh, "visradSource", "")
    _set!(sh, "cadPrompt", "")
    for k in ("compInfo", "selTitle", "partKind", "partName", "partFile", "bnBy", "bn_z", "bn_r", "bn_phi",
              "bn_defocus", "bn_power", "bn_aim_x", "bn_aim_y", "bn_aim_z", "bnInfo",
              "bnCluster", "bnQuad")
        _set!(sh, k, "")
    end
    for k in ("bnMirrored", "bnCanSplit", "bnCanMerge")
        _set!(sh, k, false)
    end
    _set!(sh, "bnCount", 0)
    _set!(sh, "clusterNames", String[])
    _set!(sh, "aimChoices", ["Target chamber (TCC)"])
    _set!(sh, "bnAimFrame", 0)
    _set!(sh, "bnAimFrameLabel", "")
    _set!(sh, "beamNames", String[])
    _set!(sh, "cutSettingsRequest", false)
    _set!(sh, "cutGrabPath", "")
    _set!(sh, "facility", "nif")
    for k in ("partTwoSided", "partInward", "partAperture", "partCopy", "partSourceExists")
        _set!(sh, k, false)
    end
    _set!(sh, "bnOn", 0)
    _set!(sh, "bnMatrix", "")
    _set!(sh, "compDisplay", 0)
    _set!(sh, "displayLabels", [DISPLAY_LABELS[m] for m in DISPLAY_MODES])
    _set!(sh, "partNumber", 0)
    _set!(sh, "partFrame", 0)
    _set!(sh, "moveFrames", String[])
    _set!(sh, "moveFrame", 0)
    _set!(sh, "moveAxis", 2)
    _set!(sh, "moveStep", "0.1")
    _set!(sh, "partRole", 0)
    _set!(sh, "partFilm", false)
    _set!(sh, "cavitySummary", "")
    _set!(sh, "compIncluded", true)
    _set!(sh, "roleLabels", [PART_ROLE_LABELS[r] for r in PART_ROLES])
    _set!(sh, "partCoords", 0)
    _set!(sh, "frameLabels", [FRAME_LABELS[a] for a in TEMPLATE_FRAMES])
    _set!(sh, "windowNumber", 0)
    _set!(sh, "ortho", false)
    _set!(sh, "turntable", true)
    mark = (_...) -> (sh.dirty[] = true)
    ts, bt, dr, st, rv, tl, v, tm = sess.target, sess.beams, sess.drive, sess.optimize, sess.review,
                                    sess.talis, sess.view, sess.time
    for obs in (sess.status, sess.meshing, sess.cad_changed, ts.kind, ts.params, ts.mesh, ts.resolution,
                ts.windows, ts.parts, ts.target, ts.error, ts.pending, ts.building, ts.fitting,
                bt.beams, bt.enabled, bt.preset, bt.by, bt.path, bt.error, bt.warnings,
                dr.scenario, dr.building, dr.result, dr.evaluating, dr.field, dr.view,
                dr.cavity, dr.cavity_point, dr.films, dr.beam_model, dr.profile,
                dr.spot_radius, dr.status, dr.retune_running, dr.unconverted_ratio,
                dr.retune_status,
                st.by, st.fields, st.deltas, st.objective, st.weights, st.hold_drive, st.keepout,
                st.keepout_min, st.capsule_keepout, st.capsule_min, st.spill, st.iterations,
                st.max_evaluations, st.running, st.evaluations, st.iteration, st.status,
                st.error, st.result, bt.cluster, bt.aim_frame, ts.frames, v.extra,
                sess.sweep.subject, sess.sweep.frame, sess.sweep.axis, sess.sweep.start, sess.sweep.stop,
                sess.sweep.step, sess.sweep.clearance, sess.sweep.los, sess.sweep.los_theta,
                sess.sweep.los_phi, sess.sweep.los_object, sess.sweep.capture, sess.sweep.optimize,
                sess.sweep.quality, sess.sweep.running, sess.sweep.progress, sess.sweep.status,
                sess.sweep.error, sess.sweep.result, sess.sweep.selected,
                rv.targets, rv.hardware, rv.peak_fraction, rv.wings, rv.errors, rv.quality,
                rv.running, rv.progress, rv.status, rv.error, rv.report,
                tl.capture, tl.hardware, tl.diagnostics, tl.tolerance, tl.pivot, tl.pivot_point, tl.omega1,
                tl.running, tl.status, tl.error, tl.report,
                v.color, v.beams, v.beam_color, v.show_off, v.cut_components, v.beam_length, v.cutaway, v.display, v.selection,
                v.hidden_beams, v.shown_1w, v.shown_2w, bt.unconverted, v.cut_kind, v.cut_axis,
                v.cut_range, v.cut_point, v.cut_normal, v.cut_shift, v.cut_line, v.cut_rotation,
                v.cut_azimuth, v.cut_width,
                v.expanded, tm.pulses, tm.history, tm.retune, tm.chart, tm.cursor, tm.dante, tm.running,
                tm.progress, tm.status)
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
# A fraction as a percentage; one that rounds to zero prints without a sign ("-0.00 %").
function _pct(x::Real, d::Integer = 2)
    isfinite(x) || return "–"
    t = @sprintf("%.*f %%", d, 100x)
    return startswith(t, "-") && all(c -> c in "0.", t[2:end-2]) ? t[2:end] : t
end
# A laser power (TW) and its fraction of the delivered power: four significant digits
# and two decimals of a percent; a power below 1 GW as "< 0.001" and a fraction below
# 0.005 % as "< 0.01 %" (when the beams reach only the wall, the other parts sum to
# round-off, ~10⁻¹³ TW).
_power_text(p::Real) = abs(p) < 1e-3 ? "< 0.001" : @sprintf("%.4g", p)
_power_pct(x::Real) = isfinite(x) && abs(x) < 5e-5 ? "< 0.01 %" : _pct(x)
_g(x::Real) = isfinite(x) ? @sprintf("%.3g", x) : "–"

function _refresh!(sh::BriareusShell)
    sess = sh.sess
    file = isempty(sh.session_path) ? "" : basename(sh.session_path) * " | "
    _set!(sh, "title", "Briareus | " * file * target_summary(sess.target))
    _set!(sh, "sessionPath", sh.session_path)
    _set!(sh, "canUndo", can_undo(sess))
    _set!(sh, "canRedo", can_redo(sess))
    _set!(sh, "visradSource", basename(cad_source(sess, sess.visrad_source[])))
    changed = sess.cad_changed[]
    _set!(sh, "cadPrompt", isempty(changed) ? "" :
                           "These CAD files changed on disk since the session was saved:\n\n" *
                           join(cad_source.(Ref(sess), changed), "\n") *
                           "\n\nThe design reads the copies saved in the session. Read the changed files instead?")
    _set!(sh, "status", sess.status[])
    _refresh_view!(sh)
    _refresh_target!(sh)
    _refresh_beams!(sh)
    _refresh_selection!(sh)
    _refresh_drive!(sh)
    _refresh_time!(sh)
    _refresh_optimize!(sh)
    _refresh_review!(sh)
    _refresh_sweep!(sh)
    tree = [TreeItem(r...) for r in model_tree(sess)]
    _set!(sh, "selectionCount", length(selected_ids(sess.view)))
    _sync_rows!(sh.rows[:treeModel], sh.models[:treeModel], tree)
    # the flux map (or the x-ray image) when chosen, in place of the 3D view or beside it
    side = sess.drive.side[]
    fig = !side && sess.drive.view[] in (:map, :image) ? sh.map.fig : sh.scene.fig
    sidefig = side ? sh.map.fig : sh.blank
    if fig !== sh.figure
        sh.figure = fig
        sh.app["canvas"] = fig
    end
    if sidefig !== sh.side_figure
        sh.side_figure = sidefig
        sh.app["sideCanvas"] = sidefig
    end
    _set!(sh, "sideView", side)
    return
end

function _refresh_view!(sh::BriareusShell)
    v = sh.sess.view
    _set!(sh, "colorBy", String(v.color[]))
    _set!(sh, "beamStyle", String(v.beams[]))
    _set!(sh, "beamColor", String(v.beam_color[]))
    _set!(sh, "beamLength", v.beam_length[] > 0 ? _num(v.beam_length[]) : "")
    _set!(sh, "cutaway", v.cutaway[])
    _set!(sh, "showOff", v.show_off[])
    _set!(sh, "cutWidth", v.cut_width[])
    _set!(sh, "cutAzimuth", v.cut_azimuth[])
    _set!(sh, "cutKind", String(v.cut_kind[]))
    _set!(sh, "cutAxis", String(v.cut_axis[]))
    lo, hi = v.cut_range[]
    _set!(sh, "cutLo", isfinite(lo) ? _num(lo) : "")
    _set!(sh, "cutHi", isfinite(hi) ? _num(hi) : "")
    n = v.cut_normal[]
    # (components below the digits shown print as 0, not as -0 or 1e-17)
    _set!(sh, "cutNormal", join((@sprintf("%.3g", abs(x) < 5e-4 ? 0.0 : x) for x in n), ", "))
    _set!(sh, "cutOffset", cut_offset(v))
    θ, φ = cut_angles(v)
    _set!(sh, "cutTheta", @sprintf("%.2f", θ))
    _set!(sh, "cutPhi", @sprintf("%.2f", φ))
    _set!(sh, "cutPoint", [_num(round(x; digits = 4) + 0.0) for x in v.cut_point[]])
    line = v.cut_line[]
    _set!(sh, "cutLine", any(!=(0), line))
    lθ, lφ = acosd(clamp(line[3], -1.0, 1.0)), mod(atand(line[2], line[1]), 360)
    _set!(sh, "cutLineText", any(!=(0), line) ? @sprintf("θ = %.1f°, φ = %.1f°", lθ, lφ) : "")
    # the Direction fields start at the line's angles (or 90°, 0°), its origin at TCC or not
    _set!(sh, "cutLineTheta", any(!=(0), line) ? @sprintf("%.2f", lθ) : "90")
    _set!(sh, "cutLinePhi", any(!=(0), line) ? @sprintf("%.2f", lφ) : "0")
    _set!(sh, "cutAtTcc", norm(Briareus.Vec3(v.cut_point[]...)) < 1e-9)
    _set!(sh, "cutRotation", v.cut_rotation[])
    sel = selected_component(sh.sess)
    _set!(sh, "cutSelected", sel === nothing ? "" : component_label(sel))
    _set!(sh, "cutScale", 1.1 * sh.sess.target.target[].scale)
    t = sh.sess.target.target[]
    _set!(sh, "cutChoices", [String(c.name) for c in t.components if c.kind === :surface])
    _set!(sh, "cutLabels", [component_label(c.name) for c in t.components if c.kind === :surface])
    _set!(sh, "cutSet", sort([String(n) for n in v.cut_components[]]))
    busy = String[]
    dr, st, rv, tl = sh.sess.drive, sh.sess.optimize, sh.sess.review, sh.sess.talis
    isempty(sh.sess.meshing[]) || push!(busy, sh.sess.meshing[])
    waiting = Controllers.target_waiting(sh.sess.target)
    isempty(waiting) || push!(busy, waiting)
    dr.building[] && push!(busy, "drive " * Controllers._seconds(time() - dr.started[]))
    st.running[] && push!(busy, "optimizing")
    rv.running[] && push!(busy, "clearance $(rv.progress[][1])/$(rv.progress[][2])")
    sw = sh.sess.sweep
    sw.running[] && push!(busy, "sweep $(sw.progress[][1])/$(sw.progress[][2])")
    tl.running[] && push!(busy, "TaLIS")
    tm = sh.sess.time
    tm.running[] && push!(busy, "drive history $(tm.progress[][1])/$(tm.progress[][2])")
    _set!(sh, "busy", join(busy, " · "))
    return
end

function _refresh_target!(sh::BriareusShell)
    ts = sh.sess.target
    _set!(sh, "targetKind", String(ts.kind[]))
    _set!(sh, "targetResolution", String(ts.resolution[]))
    _set!(sh, "mesh_h", _num(ts.mesh[].h))
    _set!(sh, "mesh_capsule_h", _num(ts.mesh[].capsule_h))
    _set!(sh, "targetSummary", target_summary(ts))
    _set!(sh, "targetError", ts.error[])
    fields = [FieldRow(String(k), label, unit, _num(ts.params[][k]))
              for (k, label, unit, _) in getfield(TARGET_FIELDS, ts.kind[])]
    _sync_rows!(sh.rows[:targetFieldModel], sh.models[:targetFieldModel], fields)
    ws = ts.kind[] === :hohlraum ? ts.windows[] : WindowSpec[]
    rows = [WindowRow(i, String(window_name(i)), _num(w.phi), _num(w.z), _num(w.width),
                      _num(w.height), w.open, _num(w.albedo), _num(w.conversion), w.round)
            for (i, w) in enumerate(ws)]
    _sync_rows!(sh.rows[:windowModel], sh.models[:windowModel], rows)
    _set!(sh, "windowCount", length(ws))
    _set!(sh, "partCount", length(ts.parts[]))
    names = component_names(ts.target[])
    _set!(sh, "hasCapsule", :capsule in names)
    _set!(sh, "hasFoil", FOIL_NAME in names)
    _set!(sh, "lehClosure", get(ts.params[], :leh_closure, 0.0) != 0)
    for (k, _, _, _) in LEH_CLOSURE_FIELDS
        _set!(sh, String(k), _num(ts.params[][k]))
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
    _set!(sh, "beamsWarnings", join(bt.warnings[], "\n"))
    _set!(sh, "hasQuads", !isempty(bt.beams[]) && all(b -> !isempty(b.quad), bt.beams[]))
    _set!(sh, "beamNames", [b.name for b in bt.beams[]])
    _set!(sh, "facility", String(bt.facility[]))
    rows = [BeamRow(i, r.key, r.count, @sprintf("%.2f", r.theta), isnan(r.phi) ? "–" : @sprintf("%.2f", r.phi),
                    _fix(r.z), _fix(r.r), _fix(r.phi_offset), _fix(r.defocus), _num(r.power),
                    r.enabled === missing ? 2 : r.enabled ? 1 : 0)
            for (i, r) in enumerate(beam_rows(bt))]
    _sync_rows!(sh.rows[:beamModel], sh.models[:beamModel], rows)
    return
end

const _BY_OF_KIND = Dict(:group => :group, :quad => :quad, :beam => :beam)

# The inspector: what is selected and its editable values.
function _refresh_selection!(sh::BriareusShell)
    sess = sh.sess
    id = sess.view.selection[]
    kind, key = selection_kind(id)
    _set!(sh, "selection", id)
    _set!(sh, "selKind", String(kind))
    _set!(sh, "selKey", key)
    t = sess.target.target[]
    if kind in (:comp, :part, :window)
        name = selected_component(sess)
        # a component left out of the calculations is described from the ghost
        g = sess.target.ghost[]
        tc = any(c -> c.name === name, t.components) || g === nothing ? t : g
        ci = findfirst(c -> c.name === name, tc.components)
        _set!(sh, "compIncluded", !(name in sess.target.excluded[]))
        if ci !== nothing
            c = tc.components[ci]
            # an entrance hole is shown with its closing lip (a split for the calculations)
            li = findfirst(l -> is_lip(l) && shown_component(l.name) === name, tc.components)
            es = li === nothing ? elements(tc, name) : [elements(tc, name); elements(tc, tc.components[li].name)]
            A = sum(tc.elem_area[e] for e in es; init = 0.0)
            what = li !== nothing ?
                   @sprintf("closing entrance hole (rays pass; its outer lip lets x-rays escape until %.3g ns, then re-emits them)",
                            tc.components[li].material.albedo.closure) :
                   c.kind === :aperture ? "aperture (transparent: rays pass, x-rays escape)" :
                   c.two_sided ? "two-sided surface" : "one-sided surface"
            mat = c.kind === :aperture ? "" :
                  @sprintf(" · albedo %.3g · conversion %.3g", albedo(c.material), c.material.conversion)
            surf = c.surface === nothing ? "faceted" : "exact " * replace(string(nameof(typeof(c.surface))), "Surface" => "") * " surface"
            _set!(sh, "compInfo", @sprintf("%s · %d elements · %.4g mm² · %s%s%s", what, length(es), A,
                                           surf, mat, tc === t ? "" : " · not in the calculations"))
            _set!(sh, "compDisplay", findfirst(==(display_mode(sess.view, name)), DISPLAY_MODES) - 1)
            _set!(sh, "selTitle", String(name))
        end
    end
    if kind === :part
        i = part_index(sess.target, Symbol(key))
        if i !== nothing
            p = sess.target.parts[][i]
            _set!(sh, "partNumber", i)
            _set!(sh, "partKind", PART_KIND_LABELS[p.kind])
            _set!(sh, "partName", String(p.name))
            src = cad_source(sess, p.file)
            _set!(sh, "partFile", src)
            _set!(sh, "partCopy", src != p.file)            # reads the copy saved in the session
            _set!(sh, "partSourceExists", src != p.file && isfile(src))
            # its frame: the target's, or another part's
            frames = _part_frame_choices(sess, p)
            _set!(sh, "frameLabels", [_frame_label(f) for f in frames])
            _set!(sh, "partFrame", something(findfirst(==(p.frame), frames), 0) - 1)
            # the move row: along an axis of any frame (the part's own frame by default)
            mf = _move_frame_choices(sess)
            _set!(sh, "moveFrames", [_frame_label(f) for f in mf])
            # (the part's own frame when another part is selected; kept by name as frames come and go)
            get(sh.shown, "_movePart", nothing) === p.name || (sh.shown["_movePart"] = p.name;
                                                                sh.shown["_moveFrameName"] = p.name)
            _set!(sh, "moveFrame", something(findfirst(==(_move_frame(sh)), mf), 1) - 1)
            _set!(sh, "partCoords", findfirst(==(p.coords), PART_COORDINATES) - 1)
            _set!(sh, "partTwoSided", p.two_sided)
            _set!(sh, "partInward", p.inward)
            _set!(sh, "partAperture", p.aperture)
            _set!(sh, "partRole", findfirst(==(part_role(p)), PART_ROLES) - 1)
            _set!(sh, "partFilm", p.name in sess.drive.films[])
            pos = Dict(part_position(p))
            fields = [[FieldRow(String(k), label, unit, _num(p.values[k]))
                       for (k, label, unit, _) in PART_SIZE_FIELDS[p.kind]];
                      [FieldRow(String(k), label, unit, _num(pos[k]))
                       for (k, label, unit) in PART_POSITION_FIELDS[p.coords]];
                      [FieldRow(String(k), label, k === :h ? part_size_unit(p) : unit, _num(p.values[k]))
                       for (k, label, unit, _) in PART_FIELDS if !(k in (:x, :y, :z))]]
            _sync_rows!(sh.rows[:partFieldModel], sh.models[:partFieldModel], fields)
        end
    elseif kind === :window
        _set!(sh, "windowNumber", parse(Int, key))
    elseif kind in (:group, :quad, :beam)
        bt = sess.beams
        idx = selected_beams(sess)
        if !isempty(idx)
            n = beam_set_values(bt, idx)
            bs = bt.beams[]
            nsel = length(selected_ids(sess.view))
            gname = split(key, '@')
            _set!(sh, "selTitle", nsel > 1 ? "$(nsel) selections · $(n.count) beams" :
                                  kind === :group ? "Group " * key :
                                  kind === :quad ? "Quad " * gname[1] * (length(gname) > 1 ? " (" * gname[2] * ")" : "") :
                                  "Beam " * key)
            _set!(sh, "bnBy", nsel > 1 ? "set" : String(_BY_OF_KIND[kind]))
            for (f, v) in ((:z, n.z), (:r, n.r), (:phi, n.phi), (:defocus, n.defocus),
                           (:aim_x, n.aim_local[1]), (:aim_y, n.aim_local[2]), (:aim_z, n.aim_local[3]))
                _set!(sh, "bn_" * String(f), _fix(v))
            end
            _set!(sh, "bn_power", _num(n.power))
            _set!(sh, "bnMirrored", n.mirrored)
            _set!(sh, "bnOn", n.enabled === missing ? 2 : n.enabled ? 1 : 0)
            _set!(sh, "bnMatrix", beam_matrix(sess, idx))
            _set!(sh, "bnCount", n.count)
            cls = unique(bt.cluster[][idx])
            _set!(sh, "bnCluster", length(cls) == 1 ? String(cls[1]) : "")
            _set!(sh, "clusterNames", String.(clusters(bt)))
            # the frame the beams aim in (its axes and origin for the aim point's fields)
            choices = _aim_choices(sess)
            _set!(sh, "aimChoices", [_frame_label(c) for c in choices])
            f = n.aim_frame
            _set!(sh, "bnAimFrame", f === nothing ? -1 : something(findfirst(==(f), choices), 1) - 1)
            _set!(sh, "bnAimFrameLabel", f === nothing ? "the chamber (the beams aim in different frames)" :
                                         f === :tcc ? "the chamber" : _frame_label(f))
            # quads: one NIF quad beam can be split; a quad's four beams merged
            qs = unique(bs[i].quad for i in idx)
            q = length(qs) == 1 ? qs[1] : ""
            _set!(sh, "bnQuad", q)
            nq = isempty(q) ? 0 : count(b -> b.quad == q, bs)
            _set!(sh, "bnCanSplit", nq == 1 && startswith(q, "Q") && length(idx) == 1)
            _set!(sh, "bnCanMerge", nq == 4 && startswith(q, "Q"))
            b = bs[idx[1]]
            p = beam_path(t, b)
            land = n.count > 1 ? "" : p.component === :missed ? " · misses the target" :
                   " · axis lands on " * _in_text(p.component)
            spot = @sprintf("spot %.3g × %.3g mm (order %.3g), f/%.3g", b.spot.a, b.spot.b,
                            b.spot.order, b.fnumber)
            θ, φ = port_angles(b)
            angles = n.count == 1 ? @sprintf("port θ %.2f°, φ %.2f°", θ, φ) : ""
            _set!(sh, "bnInfo", join(filter(!isempty, [@sprintf("%d beam%s", n.count, n.count == 1 ? "" : "s"),
                                                        angles, spot]), " · ") * land)
        end
    end
    return
end

function _refresh_drive!(sh::BriareusShell)
    sess = sh.sess
    dr = sess.drive
    res = dr.result[]
    _set!(sh, "driveBuilding", dr.building[])
    _set!(sh, "driveElapsed", dr.building[] ? Controllers._seconds(time() - dr.started[]) : "")
    _set!(sh, "driveEvaluating", dr.evaluating[])
    _set!(sh, "driveBuilt", dr.scenario[] !== nothing)
    _set!(sh, "driveStale", drive_stale(sess))
    _set!(sh, "driveCurrent", Controllers._drive_current(sess))
    _set!(sh, "driveStatus", dr.status[])
    _set!(sh, "driveHasResult", res !== nothing)
    _set!(sh, "driveField", String(dr.field[]))
    _set!(sh, "driveView", String(dr.view[]))
    _set!(sh, "spotRadius", _num(dr.spot_radius[]))
    r1, r2 = dr.unconverted_ratio[]
    _set!(sh, "ucRatio1", r1 > 0 ? _num(r1) : "")
    _set!(sh, "ucRatio2", r2 > 0 ? _num(r2) : "")
    _set!(sh, "ucLive", [r1 > 0, r2 > 0])
    _set!(sh, "unconvertedText", res === nothing ? "" : unconverted_summary(res))
    _set!(sh, "driveSummary", res === nothing ? "" : drive_summary(res))
    _set!(sh, "driveMean", res === nothing || isnan(res.mean_flux) ? "" :
                           @sprintf("%.3g W/cm² · Tr %.1f eV", to_Wcm2(res.mean_flux),
                                    radiation_temperature(res.mean_flux)))
    s = res === nothing ? nothing : res.spot
    _set!(sh, "spotText", s === nothing ? "" :
                          @sprintf("in r ≤ %.3g mm: mean %.2e W/cm², peak %.2e W/cm², rms nonuniformity %.1f %%, %.1f %% of the delivered power",
                                   s.radius, s.mean, s.peak, 100 * s.rms, 100 * s.fraction))
    modes = res === nothing ? ModeRow[] :
            [ModeRow(m.l, _pct(m.epsilon), _pct(m.legendre), "", "") for m in res.modes]
    _sync_rows!(sh.rows[:modeModel], sh.models[:modeModel], modes)
    power = res === nothing ? PowerRow[] :
            [PowerRow(p.name, _power_text(p.power), _power_pct(p.fraction))
             for p in res.power]
    _sync_rows!(sh.rows[:powerModel], sh.models[:powerModel], power)
    _set!(sh, "delivered", res === nothing ? "" :
                           @sprintf("%.4g TW delivered", res.evaluation.deposition.delivered))
    _set!(sh, "driveCavity", dr.cavity[])
    _set!(sh, "beamModel", String(dr.beam_model[]))
    _set!(sh, "cavityPoint", join((_num(x) for x in dr.cavity_point[]), ", "))
    # the cavity's account changes with the drive's build only
    if get(sh.shown, "_enclosure", missing) !== dr.enclosure[]
        sh.shown["_enclosure"] = dr.enclosure[]
        _set!(sh, "cavitySummary", cavity_summary(dr))
    end
    _set!(sh, "retuneRunning", dr.retune_running[])
    _set!(sh, "retuneStatus", dr.retune_status[])
    return
end

function _refresh_time!(sh::BriareusShell)
    sess = sh.sess
    tm = sess.time
    ch = time_chart(sess)
    h = tm.history[]
    _set!(sh, "timeChart", String(ch.chart))
    _set!(sh, "timeTitle", ch.title)
    _set!(sh, "timeEmpty", ch.empty)
    _set!(sh, "timeRunning", tm.running[])
    k, n = tm.progress[]
    _set!(sh, "timeProgress", tm.running[] ? @sprintf("step %d of %d · %s", k, n,
                                                       Controllers._seconds(time() - tm.started[])) : "")
    _set!(sh, "timeFraction", n > 0 ? k / n : 0.0)
    _set!(sh, "timeStatus", tm.status[])
    _set!(sh, "timeHasPulses", tm.pulses[] !== nothing)
    _set!(sh, "timeHasHistory", h !== nothing)
    _set!(sh, "timeStale", history_stale(sess))
    _set!(sh, "timeSteps", h === nothing ? 0 : length(h.time))
    step = time_step(tm)
    _set!(sh, "timeStep", step === nothing ? 0 : step)
    _set!(sh, "timeLabel", time_label(sess))
    _set!(sh, "timeDante", tm.dante[])
    _sync_rows!(sh.rows[:timeLegendModel], sh.models[:timeLegendModel],
                LegendRow[LegendRow(r.label, r.color, String(r.style)) for r in ch.legend])
    return
end

function _refresh_optimize!(sh::BriareusShell)
    st = sh.sess.optimize
    _set!(sh, "optBy", String(st.by[]))
    for f in OPTIMIZE_FIELDS
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
    _set!(sh, "optCapsule", st.capsule_keepout[])
    _set!(sh, "optCapsuleMin", _num(st.capsule_min[]))
    _set!(sh, "optSpill", _num(st.spill[]))
    _set!(sh, "optIterations", st.iterations[])
    _set!(sh, "optMaxEvaluations", st.max_evaluations[])
    _set!(sh, "optRunning", st.running[])
    _set!(sh, "optProgress", optimize_progress(st))
    _set!(sh, "optStatus", st.status[])
    _set!(sh, "optError", st.error[])
    r = st.result[]
    _set!(sh, "optHasResult", r !== nothing)
    _set!(sh, "optApplied", r !== nothing && r.applied)
    _set!(sh, "optValues", r === nothing ? "" : st.status[])
    _set!(sh, "optClipped", r === nothing ? "" :
                            "clipped laser power " * _pct(r.before.clipped) * " → " * _pct(r.after.clipped) *
                            (isfinite(r.after.capsule_clearance) ?
                             " · on the capsule " * _pct(r.before.capsule, 3) * " → " * _pct(r.after.capsule, 3) *
                             @sprintf(" · capsule worst-case clearance %.3f → %.3f mm",
                                      r.before.capsule_clearance, r.after.capsule_clearance) : ""))
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

# A component named within a sentence: "the top cap"-style labels in lower case, a
# user's part name as it is.
_in_text(n::Symbol) = (l = component_label(n); l == String(n) ? l : lowercasefirst(l))

function _clear_row(x)
    return ClearRow(x.beam, component_label(x.component), String(x.status), _fix(x.clearance),
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
    _set!(sh, "rvStale", review_stale(rv, sess.target.target[], active_beams(sess.beams)))
    comps = [CompRow(String(n), component_label(n), n in rv.targets[] ? "target" : n in rv.hardware[] ? "hardware" : "none")
             for n in surface_names(sess.target.target[])]
    _sync_rows!(sh.rows[:compModel], sh.models[:compModel], comps)
    _sync_rows!(sh.rows[:clearModel], sh.models[:clearModel], ClearRow[_clear_row(x) for x in clearance_rows(rv)])
    _sync_rows!(sh.rows[:worstModel], sh.models[:worstModel], ClearRow[_clear_row(x) for x in worst_rows(rv)])
    _refresh_talis!(sh)
    return
end

function _refresh_sweep!(sh::BriareusShell)
    sess = sh.sess
    sw = sess.sweep
    subs = sweep_subjects(sess)
    _set!(sh, "swSubjects", [s === :backlighter && sess.target.kind[] === :foil ? "Backlighter foil (template)" :
                             String(s) for s in subs])
    _set!(sh, "swSubject", something(findfirst(==(sweep_subject(sess)), subs), 0) - 1)
    mf = _move_frame_choices(sess)
    _set!(sh, "swFrames", [_frame_label(f) for f in mf])
    _set!(sh, "swFrame", something(findfirst(==(sweep_frame(sess)), mf), 1) - 1)
    _set!(sh, "swAxis", something(findfirst(==(sw.axis[]), SWEEP_AXES), 3) - 1)
    for k in (:start, :stop, :step, :los_theta, :los_phi)
        _set!(sh, "sw_" * String(k), _num(getfield(sw, k)[]))
    end
    for k in (:clearance, :los, :capture, :optimize)
        _set!(sh, "sw_" * String(k), getfield(sw, k)[])
    end
    objs = [:none; surface_names(sess.target.target[])]
    _set!(sh, "swObjects", [o === :none ? "(none)" : component_label(o) for o in objs])
    _set!(sh, "swObject", something(findfirst(==(sweep_los_object(sess)), objs), 1) - 1)
    _set!(sh, "swQuality", String(sw.quality[]))
    _set!(sh, "swRunning", sw.running[])
    _set!(sh, "swDone", sw.progress[][1])
    _set!(sh, "swTotal", sw.progress[][2])
    _set!(sh, "swStatus", sw.status[])
    _set!(sh, "swError", sw.error[])
    _set!(sh, "swHasResult", sw.result[] !== nothing)
    _set!(sh, "swBrief", sweep_brief(sw))
    _set!(sh, "swStale", sweep_stale(sess))
    _set!(sh, "swSelected", sw.selected[])
    steps = try
        length(sweep_values(sw))
    catch
        0
    end
    _set!(sh, "swSteps", steps)
    if get(sh.shown, "_sweepResult", missing) !== sw.result[]
        sh.shown["_sweepResult"] = sw.result[]
        rows = [SweepRowItem(k, @sprintf("%+.4g", w.value), _fix(w.clearance), w.worst, w.fails, w.reviews, w.los,
                             _pct(w.blocked, 1), _pct(w.covered, 1), _pct(w.captured, 2),
                             isnan(w.objective_after) ? "–" : @sprintf("%.3g → %.3g", w.objective_before, w.objective_after))
                for (k, w) in enumerate(sweep_rows(sw))]
        _sync_rows!(sh.rows[:sweepModel], sh.models[:sweepModel], rows)
    end
    ch = sweep_chart(sess)
    _set!(sh, "swEmpty", ch.empty)
    _set!(sh, "swTitle", ch.title)
    _sync_rows!(sh.rows[:sweepLegendModel], sh.models[:sweepLegendModel],
                LegendRow[LegendRow(r.label, r.color, String(r.style)) for r in ch.legend])
    return
end

_unit(u::Symbol) = u === :deg ? "°" : String(u)
# A margin or distance; ±Inf means the rule found nothing to measure (no free edge, all
# hardware downstream of the landing point).
_val(x::Real, u::Symbol) = isfinite(x) ? @sprintf("%+.3f %s", x, _unit(u)) : x > 0 ? "n/a" : "–"
_val0(x::Real, u::Symbol) = isfinite(x) ? @sprintf("%.3f %s", x, _unit(u)) : x > 0 ? "n/a" : "–"

function _talis_row(x::TaLISRow)
    return TalisRow(x.beam, get(TALIS_RULE_LABELS, x.rule, String(x.rule)), component_label(x.component),
                    String(x.status), _val(x.margin, x.unit), _val0(x.worst, x.unit),
                    _val0(x.nominal, x.unit), _val0(x.required, x.unit), x.note)
end

function _refresh_talis!(sh::BriareusShell)
    sess = sh.sess
    tl = sess.talis
    _set!(sh, "tlTolerance", _num(tl.tolerance[]))
    _set!(sh, "tlPivot", String(tl.pivot[]))
    _set!(sh, "tlPivotPoint", join((_num(v) for v in tl.pivot_point[]), ", "))
    _set!(sh, "tlOmega1", tl.omega1[])
    _set!(sh, "tlRunning", tl.running[])
    _set!(sh, "tlElapsed", tl.running[] ? Controllers._seconds(time() - tl.started[]) : "")
    _set!(sh, "tlStatus", tl.status[])
    _set!(sh, "tlError", tl.error[])
    _set!(sh, "tlHasReport", tl.report[] !== nothing)
    _set!(sh, "tlBrief", talis_brief(tl))
    _set!(sh, "tlStale", talis_stale(tl, sess.review, sess.target.target[], active_beams(sess.beams)))
    comps = [CompRow(String(n), component_label(n), String(talis_role(tl, n)))
             for n in surface_names(sess.target.target[])]
    _sync_rows!(sh.rows[:talisCompModel], sh.models[:talisCompModel], comps)
    # the rows change only with the report (hundreds of them: not on every refresh)
    if get(sh.shown, "_talisReport", missing) !== tl.report[]
        sh.shown["_talisReport"] = tl.report[]
        rules = [RuleRow(r.label, r.requirement, r.pass, r.review, r.fail,
                         _val(r.worst.margin, r.worst.unit), r.worst.beam * " / " * component_label(r.worst.component),
                         r.fail > 0 ? "fail" : r.review > 0 ? "review" : "pass")
                 for r in talis_rule_rows(tl)]
        _sync_rows!(sh.rows[:ruleModel], sh.models[:ruleModel], rules)
        _sync_rows!(sh.rows[:talisModel], sh.models[:talisModel],
                    TalisRow[_talis_row(x) for x in talis_rows(tl)])
    end
    _sync_rows!(sh.rows[:noteModel], sh.models[:noteModel],
                NoteRow[NoteRow(String(n.kind), n.text) for n in talis_notes(tl)])
    return
end

# ---------------------------------------------------------------- QML callbacks

const _SHELL = Ref{Union{Nothing,BriareusShell}}(nothing)
# `_TICK_HOOK[](shell)` runs on every tick (scripted smoke tests drive the window
# through it); `request_close()` closes the window from Julia.
const _TICK_HOOK = Ref{Any}(nothing)
const _ON_SHOWN = Ref{Any}(nothing)            # `briareus(; shown)`, until the first tick
const _CLOSE_REQUESTED = Ref(false)
const _LAST_REFRESH = Ref(0.0)

"""
    request_close()

Close the open Briareus window (from a callback or a tick hook).
"""
request_close() = (_CLOSE_REQUESTED[] = true; nothing)

"""
    request_grab(path)

Save an image of the open window's body (panels and viewer) to `path` (PNG) on the
window's next tick, through QML `grabToImage`. The file appears asynchronously.
"""
request_grab(path::AbstractString) = (sh = _SHELL[]; sh === nothing || _set!(sh, "grabPath", String(path)); nothing)
br_grab_started() = (sh = _SHELL[]; sh === nothing || _set!(sh, "grabPath", ""); nothing)

"""
    request_tab(key)

Show a results tab of the open window (`"beams"`, `"drive"`, `"time"`, `"optimize"`,
`"clearance"`, `"sweep"`, `"talis"`), as its menus do.
"""
request_tab(key::AbstractString) = (sh = _SHELL[]; sh === nothing ||
                                    (_set!(sh, "resultsRequest", ""); _set!(sh, "resultsRequest", String(key))); nothing)

"""
    request_cut_settings(path = "")

Open the cutaway's Cut settings popup of the open window, as the Cutaway button's
context menu does, and with a `path` save an image of the popup there (PNG, written
asynchronously, as by [`request_grab`](@ref)).
"""
function request_cut_settings(path::AbstractString = "")
    sh = _SHELL[]
    sh === nothing && return nothing
    _set!(sh, "cutGrabPath", String(path))
    _set!(sh, "cutSettingsRequest", false)
    _set!(sh, "cutSettingsRequest", true)
    return nothing
end

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

# The Briareus icon (docs/make_icons.jl writes it).
const ICON_ICO = joinpath(QML_DIR, "icons", "briareus.ico")
const _ICON_TRIES = Ref(0)

# The taskbar identity of the process (Windows): set before the first window, it keeps
# the window out of julia.exe's taskbar group, whose button shows julia.exe's icon.
const APP_USER_MODEL_ID = "Briareus.DesignWindow"

# Windows: show the icon in the title bar and taskbar. Qt Quick windows have no icon
# property and QML.jl does not expose QGuiApplication's, so the icon is sent to this
# process's top-level Qt windows on the first ticks after the window opens, and set on
# their window class, so that a native window Qt creates later carries it too. Windows
# are found by process, not by title: the title changes while a design loads. Returns
# true once a visible window has the icon. Elsewhere the window keeps the platform
# default.
function _set_window_icon()
    Sys.iswindows() || return true
    empty!(_QT_HWNDS)
    ccall((:EnumWindows, "user32"), Cint, (Ptr{Cvoid}, Int), _ENUM_QT_WINDOWS[], 0)
    isempty(_QT_HWNDS) && return false
    shown = false
    for (which, index, metric) in ((0, -34, 49), (1, -14, 11))  # ICON_SMALL/GCLP_HICONSM/SM_CXSMICON, ICON_BIG/GCLP_HICON/SM_CXICON
        px = ccall((:GetSystemMetrics, "user32"), Cint, (Cint,), metric)
        icon = ccall((:LoadImageW, "user32"), Ptr{Cvoid},
                     (Ptr{Cvoid}, Cwstring, Cuint, Cint, Cint, Cuint),
                     C_NULL, ICON_ICO, 1, px, px, 0x10)      # IMAGE_ICON, LR_LOADFROMFILE
        icon == C_NULL && continue
        for hwnd in _QT_HWNDS
            ccall((:SetClassLongPtrW, "user32"), Csize_t, (Ptr{Cvoid}, Cint, Ptr{Cvoid}),
                  hwnd, index, icon)
            ccall((:SendMessageW, "user32"), Ptr{Cvoid}, (Ptr{Cvoid}, Cuint, Csize_t, Ptr{Cvoid}),
                  hwnd, 0x0080, which, icon)                 # WM_SETICON
            shown |= ccall((:IsWindowVisible, "user32"), Cint, (Ptr{Cvoid},), hwnd) != 0
        end
    end
    return shown
end

const _QT_HWNDS = Ptr{Cvoid}[]

# EnumWindows callback: collect this process's top-level Qt windows (class "Qt…QWindow…").
function _collect_qt_window(hwnd::Ptr{Cvoid}, ::Int)::Cint
    pid = Ref{UInt32}(0)
    ccall((:GetWindowThreadProcessId, "user32"), UInt32, (Ptr{Cvoid}, Ref{UInt32}), hwnd, pid)
    pid[] == getpid() || return 1
    buf = Vector{UInt16}(undef, 256)
    n = ccall((:GetClassNameW, "user32"), Cint, (Ptr{Cvoid}, Ptr{UInt16}, Cint), hwnd, buf, 256)
    cls = transcode(String, buf[1:max(n, 0)])
    startswith(cls, "Qt") && occursin("QWindow", cls) && push!(_QT_HWNDS, hwnd)
    return 1
end

const _ENUM_QT_WINDOWS = Ref(C_NULL)
_enum_qt_windows_ptr() = @cfunction(_collect_qt_window, Cint, (Ptr{Cvoid}, Int))

# Returns 0 (nothing to draw), 1 (redraw the canvas), or 2 (close the window).
function br_tick()
    sh = _SHELL[]
    sh === nothing && return 0
    # Only thread 1 runs libuv's event loop, and it sits in Qt's exec(): pump it here.
    Base.process_events()
    if _ICON_TRIES[] > 0
        _ICON_TRIES[] = _set_window_icon() ? 0 : _ICON_TRIES[] - 1
    end
    if _ON_SHOWN[] !== nothing
        f, _ON_SHOWN[] = _ON_SHOWN[], nothing
        try
            f()
        catch err
            @error "shown callback failed" exception = (err, catch_backtrace())
        end
    end
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
    # after any action: note the design for undo, then show the new state
    if sh.dirty[]
        sh.dirty[] = false
        record!(sess)
        _refresh!(sh)
        _LAST_REFRESH[] = time()
    elseif (sess.drive.building[] || sess.optimize.running[] || sess.review.running[] ||
            sess.talis.running[] || sess.time.running[] || sess.sweep.running[] ||
            !target_ready(sess.target)) &&
           time() - _LAST_REFRESH[] > 0.25
        # running jobs report elapsed time, four times a second (not a design edit: a
        # snapshot for undo every tick took a third of the GUI thread during a drive)
        _refresh!(sh)
        _LAST_REFRESH[] = time()
    end
    redraw = sh.scene.dirty[] || sh.map.dirty[]
    chart = sh.time.dirty[]
    swchart = sh.sweep.dirty[]
    sh.scene.dirty[] = false
    sh.map.dirty[] = false
    sh.time.dirty[] = false
    sh.sweep.dirty[] = false
    return (redraw ? 1 : 0) | (chart ? 4 : 0) | (swchart ? 8 : 0)
end

_sym(x) = Symbol(String(x))
_parse(text, what) = Controllers._parse_number(string(text), what)

# Model tree and selection
br_select(id) = _with_shell(sh -> select!(sh.sess.view, String(id)))
br_select_add(id) = _with_shell(sh -> select_add!(sh.sess.view, String(id)))
br_select_range(id) = _with_shell(sh -> select_range!(sh.sess.view, [r.id for r in model_tree(sh.sess)],
                                                      String(id)))

# The nodes a tree action on `id` applies to: the whole selection when `id` is in it.
function _acted_nodes(sh, id)
    ids = selected_ids(sh.sess.view)
    return String(id) in ids ? ids : [String(id)]
end
# Component name of a component node (`nothing` for other nodes).
function _node_component(sess, id)
    kind, key = selection_kind(id)
    kind in (:comp, :part) && return Symbol(key)
    kind === :window && return window_name(parse(Int, key))
    return nothing
end
_node_beam_names(sess, id) = [sess.beams.beams[][i].name for i in node_beams(sess.beams, id)]

# Whether a node is shown in the 3D view and takes part in the calculations, for the
# tree's right-click menu: [shown, included].
function br_node_flags(id)
    sh = _SHELL[]
    sh === nothing && return [true, true]
    sess = sh.sess
    c = _node_component(sess, String(id))
    if c !== nothing
        return [display_mode(sess.view, c) !== :invisible, !(c in sess.target.excluded[])]
    end
    idx = node_beams(sess.beams, String(id))
    isempty(idx) && return [true, true]
    bs = sess.beams.beams[]
    return [all(i -> beam_drawn(sess.view, bs[i].name, 3), idx), all(i -> sess.beams.enabled[][i], idx)]
end
# Show or hide nodes: a component becomes opaque or invisible, a beam's 3ω light drawn or not.
function br_node_show(id, on)
    _with_shell() do sh
        sess = sh.sess
        for n in _acted_nodes(sh, id)
            c = _node_component(sess, n)
            if c !== nothing
                _show_component!(sess.view, c, Bool(on))
            else
                set_beams_drawn!(sess.view, _node_beam_names(sess, n), 3, Bool(on))
            end
        end
    end
end
# Showing an invisible component makes it opaque; another mode is kept.
_show_component!(v, c, on) = on ? (display_mode(v, c) === :invisible && set_display!(v, c, :opaque)) :
                                  set_display!(v, c, :invisible)
# Display mode of the components among the nodes acted on (the tree's or Properties').
function br_set_display(id, mode)
    _with_shell() do sh
        m = mode isa Real ? DISPLAY_MODES[round(Int, mode) + 1] : _sym(mode)
        cs = filter(!isnothing, [_node_component(sh.sess, n) for n in _acted_nodes(sh, id)])
        isempty(cs) || set_display!(sh.sess.view, cs, m)
    end
end
# A dot of the beam matrix: `row` "use" (takes part in the calculations) or "draw", of
# `harmonic` 1, 2 or 3, for the beams under the nodes acted on (every beam for "beams").
function br_beam_dot(id, row, h, on)
    _with_shell() do sh
        sess = sh.sess
        idx = sort!(unique!(reduce(vcat, (node_beams(sess.beams, n) for n in _acted_nodes(sh, id));
                                   init = Int[])))
        isempty(idx) && return
        h, on = Int(h), Bool(on)
        if String(row) == "use"
            set_participation!(sess.beams, idx, h, on)
        else
            set_beams_drawn!(sess.view, [sess.beams.beams[][i].name for i in idx], h, on)
        end
    end
end
function br_node_include(id, on)
    _with_shell() do sh
        sess = sh.sess
        comps = Symbol[]
        for n in _acted_nodes(sh, id)
            c = _node_component(sess, n)
            if c !== nothing
                push!(comps, c)
            else
                idx = node_beams(sess.beams, n)
                isempty(idx) || set_set_enabled!(sess.beams, idx, Bool(on))
            end
        end
        # CAD parts put back are remeshed on workers (the target follows when ready)
        isempty(comps) || include_components!(sess, comps, Bool(on))
    end
end
br_toggle_expand(id) = _with_shell(sh -> toggle_expanded!(sh.sess.view, String(id)))
function br_tree_on(id, on)
    _with_shell() do sh
        sess = sh.sess
        kind, key = selection_kind(String(id))
        on = Bool(on)
        if kind in (:group, :quad, :beam)
            idx = String(id) == sess.view.selection[] ? selected_beams(sess) : node_beams(sess.beams, String(id))
            isempty(idx) || set_set_enabled!(sess.beams, idx, on)
        elseif kind === :beams
            Controllers._assign!(sess.beams.enabled, fill(on, length(sess.beams.beams[])))
        elseif kind in (:comp, :part, :window)
            # every selected component when the box is on a selected one
            for n in _acted_nodes(sh, id)
                c = _node_component(sess, n)
                c === nothing || _show_component!(sess.view, c, on)
            end
        end
    end
end
# Target
br_set_target_kind(kind) = _with_shell(sh -> set_target_kind!(sh.sess.target, _sym(kind)))
br_set_resolution(r) = _with_shell(sh -> set_resolution!(sh.sess.target, _sym(r)))
br_edit_mesh(key, text) = _with_shell(sh -> edit_mesh!(sh.sess.target, _sym(key), string(text)))
br_edit_target(key, text) = _with_shell(sh -> edit_target!(sh.sess.target, _sym(key), string(text)))
br_leh_closure(on) = _with_shell(sh -> set_leh_closure!(sh.sess.target, Bool(on)))
function br_add_window()
    _with_shell() do sh
        n = length(sh.sess.target.windows[])
        add_window!(sh.sess.target) && select!(sh.sess.view, "window:$(n + 1)")
    end
end
br_remove_window(row) = _with_shell(sh -> (remove_window!(sh.sess.target, round(Int, row));
                                           select!(sh.sess.view, "target")))
function br_edit_window(row, key, text)
    _with_shell() do sh
        k = round(Int, row)
        # a rejected edit keeps the window: re-send its row so the field shows it again
        edit_window!(sh.sess.target, k, _sym(key), string(text)) ||
            (k <= length(sh.rows[:windowModel]) && (sh.models[:windowModel][k] = sh.rows[:windowModel][k]))
    end
end
function br_window_round(row, on)
    _with_shell() do sh
        k = round(Int, row)
        set_window_round!(sh.sess.target, k, Bool(on)) ||
            (k <= length(sh.rows[:windowModel]) && (sh.models[:windowModel][k] = sh.rows[:windowModel][k]))
    end
end
function br_window_open(row, on)
    _with_shell() do sh
        k = round(Int, row)
        set_window_open!(sh.sess.target, k, Bool(on)) ||
            (k <= length(sh.rows[:windowModel]) && (sh.models[:windowModel][k] = sh.rows[:windowModel][k]))
    end
end
function br_add_part(kind)
    _with_shell() do sh
        ts = sh.sess.target
        if add_part!(ts, _sym(kind))
            select!(sh.sess.view, "part:" * String(ts.parts[][end].name))
            sh.sess.status[] = "added " * String(ts.parts[][end].name) * " at TCC: set its size and place in Properties"
        end
    end
end
function br_add_stl(url)
    _with_shell() do sh
        ts = sh.sess.target
        path = _url_to_path(String(url))
        n0 = length(ts.parts[])
        # one part, or one per named solid of a multi-body export
        if add_stl_parts!(ts, path)
            name = String(ts.parts[][end].name)
            n = length(ts.parts[]) - n0
            select!(sh.sess.view, "part:" * String(ts.parts[][n0+1].name))
            sh.sess.status[] = n == 1 ?
                "added $name from $(basename(path)): set its scale (mm per file unit) and place in Properties" :
                "added $n parts from $(basename(path)), one per solid: set their scale (mm per file unit) and place in Properties"
        end
    end
end
function _resend_part_fields!(sh)
    rows = sh.rows[:partFieldModel]
    for i in eachindex(rows)
        sh.models[:partFieldModel][i] = rows[i]
    end
end
# The selected parts (indices), the one shown in Properties first: edits apply to all.
function _selected_parts(sh)
    ts = sh.sess.target
    idx = Int[]
    for id in selected_ids(sh.sess.view)
        kind, key = selection_kind(id)
        kind === :part || continue
        i = part_index(ts, Symbol(key))
        i === nothing || push!(idx, i)
    end
    return idx
end
function br_edit_part(key, text)
    _with_shell() do sh
        ts = sh.sess.target
        idx = _selected_parts(sh)
        isempty(idx) && return
        k = _sym(key)
        if k === :name                                  # names are per part
            i = idx[1]
            edit_part!(ts, i, k, string(text)) ?
                select!(sh.sess.view, "part:" * String(ts.parts[][i].name)) : _resend_part_fields!(sh)
            return
        end
        ok = true
        for i in idx
            haskey(ts.parts[][i].values, k) || continue  # a field of other kinds
            ok &= edit_part!(ts, i, k, string(text))
        end
        ok || _resend_part_fields!(sh)
    end
end
function br_part_flag(flag, on)
    _with_shell() do sh
        ts = sh.sess.target
        for i in _selected_parts(sh)
            set_part_flag!(ts, i, _sym(flag), Bool(on))
        end
    end
end
# A frame's name in the menus: a template frame's label, or the part's name.
_frame_label(f::Symbol) = haskey(FRAME_LABELS, f) ? FRAME_LABELS[f] : "Part: " * String(f)
# Frames a part can be placed in: the template's, then the other parts.
function _part_frame_choices(sess, p)
    fr = part_frames(sess.target)
    return [[f for f in TEMPLATE_FRAMES if haskey(fr, f) && !(f in (q.name for q in sess.target.parts[]))];
            [q.name for q in sess.target.parts[] if q.name !== p.name]]
end
# Frames to move along or aim in: the template's, then every part's.
function _move_frame_choices(sess)
    fr = part_frames(sess.target)
    names = [q.name for q in sess.target.parts[]]
    return [[f for f in TEMPLATE_FRAMES if haskey(fr, f) && !(f in names)]; names]
end
# Move the selected part to another frame, keeping it where it is (index into
# `frameLabels`).
function br_part_frame(index)
    _with_shell() do sh
        ts = sh.sess.target
        idx = _selected_parts(sh)
        isempty(idx) && return
        frame = _part_frame_choices(sh.sess, ts.parts[][idx[1]])[Int(index)+1]
        ok = true
        for i in idx
            ok &= set_part_frame!(ts, i, frame)
        end
        ok || _set!(sh, "partFrame", -1)
    end
end
# The move row: the frame (index into `moveFrames`) and axis (0 x, 1 y, 2 z) it moves
# along, and the step as text.
br_move_frame(index) = _with_shell(sh -> (sh.shown["_moveFrameName"] = _move_frame_choices(sh.sess)[Int(index)+1]))
# The frame the move row moves along: the chosen one while the design has it, else the chamber's.
function _move_frame(sh)
    f = get(sh.shown, "_moveFrameName", :tcc)
    return f in _move_frame_choices(sh.sess) ? f : :tcc
end
br_move_axis(index) = _with_shell(sh -> _set!(sh, "moveAxis", Int(index)))
br_move_step(text) = _with_shell(sh -> (_parse(text, "the step"); _set!(sh, "moveStep", String(strip(string(text))))))
# Move the selected part by ± the step along the chosen axis (`sign` +1 or -1): one edit.
function br_nudge_part(sign)
    _with_shell() do sh
        ts = sh.sess.target
        idx = _selected_parts(sh)
        isempty(idx) && return
        Δ = Float64(sign) * Controllers._parse_number(sh.shown["moveStep"], "the step"; positive = true)
        frame = _move_frame(sh)
        axis = SWEEP_AXES[sh.shown["moveAxis"]+1]
        for i in idx
            nudge_part!(ts, i, Δ; frame, axis)
        end
        sh.sess.status[] = @sprintf("moved %s by %+g mm along %s of %s", ts.parts[][idx[1]].name, Δ, axis,
                                    _frame_label(frame))
    end
end
# Set up a sweep of the selected part along the move row's axis, and show it.
function br_sweep_part()
    _with_shell() do sh
        idx = _selected_parts(sh)
        isempty(idx) && return
        p = sh.sess.target.parts[][idx[1]]
        set_sweep!(sh.sess.sweep; subject = p.name,
                   frame = _move_frame(sh),
                   axis = SWEEP_AXES[sh.shown["moveAxis"]+1])
        select!(sh.sess.view, "task:sweep")
        request_tab("sweep")
    end
end
# The selected part's role in the analyses (index into `roleLabels`); a refused role
# shows the part's own again.
function br_part_role(index)
    _with_shell() do sh
        idx = _selected_parts(sh)
        isempty(idx) && return
        ts = sh.sess.target
        role = PART_ROLES[Int(index)+1]
        name = ts.parts[][idx[1]].name
        # once given (a CAD capsule's sphere may be fitted on a worker first): select the
        # part holding the role, the fitted sphere for a CAD capsule
        done = function ()
            i = role === :none ? nothing : part_index(ts, Controllers._ROLE_NAMES[role])
            p = i === nothing ? ts.parts[][something(part_index(ts, name), idx[1])] : ts.parts[][i]
            select!(sh.sess.view, "part:" * String(p.name))
            sh.sess.status[] = "$(p.name): " * lowercase(PART_ROLE_LABELS[part_role(p)])
            sh.dirty[] = true
        end
        if set_part_role!(ts, idx[1], role; done)
            name in ts.fitting[] && (sh.sess.status[] = "$name: fitting the capsule's sphere…")
        else
            delete!(sh.shown, "partRole")
        end
    end
end
function br_part_coords(index)
    _with_shell() do sh
        for i in _selected_parts(sh)
            set_part_coords!(sh.sess.target, i, PART_COORDINATES[Int(index)+1])
        end
    end
end
function br_remove_selected()
    _with_shell() do sh
        sess = sh.sess
        kind, key = selection_kind(sess.view.selection[])
        if kind === :part
            i = part_index(sess.target, Symbol(key))
            i === nothing || (remove_part!(sess.target, i); select!(sess.view, "target"))
        elseif kind === :window
            remove_window!(sess.target, parse(Int, key))
            select!(sess.view, "target")
        else
            sess.status[] = "only parts and windows can be removed (hide a component instead)"
        end
    end
end
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
# Edit every selected beam (a group, quad, beam or Ctrl-selected set) at once.
function br_edit_node(field, text)
    _with_shell() do sh
        idx = selected_beams(sh.sess)
        isempty(idx) && return
        edit_beam_set!(sh.sess.beams, idx, _sym(field), string(text)) ||
            (delete!(sh.shown, "bn_" * String(field)))       # re-send the kept value
    end
end
function br_set_cluster(name)
    _with_shell() do sh
        sess = sh.sess
        idx = selected_beams(sess)
        isempty(idx) && return
        if set_cluster!(sess.beams, idx, string(name))
            g = strip(string(name))
            "group:" * g in sess.view.expanded[] || toggle_expanded!(sess.view, "group:" * g)
            select!(sess.view, "group:" * g)
            sess.status[] = "moved $(length(idx)) beam$(length(idx) == 1 ? "" : "s") to group $g"
        end
    end
end
function br_rename_cluster(new)
    _with_shell() do sh
        kind, key = selection_kind(sh.sess.view.selection[])
        kind === :group && rename_cluster!(sh.sess.beams, key, string(new)) &&
            select!(sh.sess.view, "group:" * strip(string(new)))
    end
end
_aim_choices(sess) = _move_frame_choices(sess)
# Hold the selected beams' aims in a frame, where they are (index into `aimChoices`).
function br_aim_frame(index)
    _with_shell() do sh
        idx = selected_beams(sh.sess)
        isempty(idx) && return
        f = _aim_choices(sh.sess)[Int(index)+1]
        set_aim_frame!(sh.sess, idx, f)
        sh.sess.status[] = f === :tcc ? "the beams aim in the chamber frame" :
                           "the beams aim in the frame of " * _frame_label(f) * " and follow it"
    end
end
# Aim the selected beams at their frame's point (a part's centroid, a template frame's origin).
function br_aim_at_frame()
    _with_shell() do sh
        idx = selected_beams(sh.sess)
        isempty(idx) && return
        fs = unique(sh.sess.beams.aim_frame[][idx])
        length(fs) == 1 || throw(ArgumentError("the beams aim in different frames: choose one first"))
        set_aim_frame!(sh.sess, idx, fs[1]; keep = false)
        sh.sess.status[] = "the beams aim at the center of " * _frame_label(fs[1])
    end
end
function br_set_facility(f)
    _with_shell() do sh
        s = sh.sess
        set_facility!(s.beams, _sym(f); target_kind = s.target.kind[], params = s.target.params[])
        s.status[] = "facility: " * FACILITY_LABELS[_sym(f)]
    end
end
function br_split_quad()
    _with_shell() do sh
        q = get(sh.shown, "bnQuad", "")
        isempty(q) && return
        g = sh.sess.beams.cluster[][findfirst(b -> b.quad == q, sh.sess.beams.beams[])]
        if split_quad!(sh.sess.beams, q)
            id = "quad:" * q * "@" * String(g)
            "group:" * String(g) in sh.sess.view.expanded[] || toggle_expanded!(sh.sess.view, "group:" * String(g))
            id in sh.sess.view.expanded[] || toggle_expanded!(sh.sess.view, id)
            select!(sh.sess.view, id)
            sh.sess.status[] = "split $q into its four beams (the drive must be computed again)"
        end
    end
end
function br_merge_quad()
    _with_shell() do sh
        q = get(sh.shown, "bnQuad", "")
        isempty(q) && return
        merge_quad!(sh.sess.beams, q) && (select!(sh.sess.view, "beam:" * q);
                                          sh.sess.status[] = "merged $q into one quad beam")
    end
end
# Session files and undo
br_undo() = _with_shell(sh -> undo!(sh.sess) || (sh.sess.status[] = "nothing to undo"))
br_redo() = _with_shell(sh -> redo!(sh.sess) || (sh.sess.status[] = "nothing to redo"))
function br_open_session(url)
    _with_shell() do sh
        path = _url_to_path(String(url))
        load_session!(sh.sess, path)
        sh.session_path = path
        select!(sh.sess.view, "target")
    end
end
# An empty url saves to the current session file; `refs` saves the CAD files by path only.
function br_save_session(url, refs)
    _with_shell() do sh
        path = isempty(String(url)) ? sh.session_path : _url_to_path(String(url))
        isempty(path) && throw(ArgumentError("choose a file to save the session to"))
        save_session(path, sh.sess; embed = !Bool(refs) && !endswith(path, ".toml"))
        sh.session_path = path
    end
end
# CAD copies saved with a session: read the changed sources, or keep the copies.
# (the prompt is answered either way: a file the target rejects stays a saved copy)
br_reload_cad_sources() = _with_shell(sh -> (reload_cad_sources!(sh.sess); sh.sess.cad_changed[] = String[]))
br_keep_cad_copies() = _with_shell(sh -> (sh.sess.cad_changed[] = String[]))
# The selected part's CAD file read from `newfile` (every part reading that file moves).
function _relink_part(sh, newfile)
    i = part_index(sh.sess.target, Symbol(selection_kind(sh.sess.view.selection[])[2]))
    i === nothing && return
    file = sh.sess.target.parts[][i].file
    relink_cad!(sh.sess, file, newfile) &&
        (sh.sess.status[] = "read the parts of $(basename(cad_source(sh.sess, file))) from $newfile")
end
br_relink_part(url) = _with_shell(sh -> _relink_part(sh, _url_to_path(String(url))))
function br_part_from_source()
    _with_shell() do sh
        i = part_index(sh.sess.target, Symbol(selection_kind(sh.sess.view.selection[])[2]))
        i === nothing || _relink_part(sh, cad_source(sh.sess, sh.sess.target.parts[][i].file))
    end
end
br_import_beams(url) = _with_shell(sh -> (import_beams!(sh.sess.beams, _url_to_path(String(url)));
                                          sh.sess.status[] = "imported " * basename(sh.sess.beams.path[])))
br_export_beams(url) = _with_shell(sh -> (export_beams(sh.sess.beams, _url_to_path(String(url)));
                                          sh.sess.status[] = "exported " * basename(sh.sess.beams.path[])))
function _visrad_io(f, url, verb)
    _with_shell() do sh
        path = _url_to_path(String(url))
        f(sh.sess.beams, path)
        sh.sess.status[] = verb * " " * basename(path)
    end
end
br_import_sst(url) = _visrad_io(import_visrad_sst!, url, "applied VisRad pointing from")
br_import_summary(url) = _visrad_io(import_visrad_summary!, url, "applied the VisRad beam summary")
br_export_sst(url) = _visrad_io(export_visrad_sst, url, "wrote VisRad pointing to")
br_export_summary(url) = _visrad_io(export_visrad_summary, url, "wrote the VisRad beam summary to")
# Drive
br_build() = _with_shell(sh -> build_scenario!(sh.sess))
br_set_view(v) = _with_shell(sh -> set_view!(sh.sess.drive, _sym(v)))
# Leave the x-ray image for the view shown before it.
br_close_image() = _with_shell(sh -> set_view!(sh.sess.drive, get(sh.shown, "_viewBeforeImage", :scene)))
br_side_view(on) = _with_shell(sh -> set_side_view!(sh.sess.drive, Bool(on)))
# An x-ray image seen along the 3D view's direction, framing what the view frames.
function br_xray_image()
    _with_shell() do sh
        c = sh.scene.camera
        f, _, u = camera_frame(c)
        width = 2 * norm(c.lookat - c.eye) * tand(c.fov / 2)
        current = Controllers._drive_current(sh.sess)
        v = sh.sess.drive.view[]
        v === :image || (sh.shown["_viewBeforeImage"] = v)
        xray_image!(sh.sess; direction = Tuple(f), up = Tuple(u), center = Tuple(c.lookat), width)
        current || (sh.sess.status[] = "computing the drive for the x-ray image…")
        current && sh.sess.drive.image[] !== nothing && (sh.sess.status[] = "x-ray image along the view")
    end
end
br_cavity(on) = _with_shell(sh -> set_cavity!(sh.sess.drive, Bool(on)))
function br_cavity_point(text)
    _with_shell() do sh
        xs = split(strip(string(text)), r"[,\s]+")
        length(xs) == 3 || throw(ArgumentError("give the cavity point as x, y, z (mm)"))
        set_cavity!(sh.sess.drive, sh.sess.drive.cavity[];
                    point = [Controllers._parse_number(x, "a coordinate") for x in xs])
        delete!(sh.shown, "cavityPoint")
    end
end
# The selected parts as thin films (left out of the drive's cavity) or not.
function br_part_film(on)
    _with_shell() do sh
        dr = sh.sess.drive
        f = copy(dr.films[])
        for i in _selected_parts(sh)
            n = sh.sess.target.parts[][i].name
            Bool(on) ? push!(f, n) : delete!(f, n)
        end
        set_cavity!(dr, dr.cavity[]; films = f)
    end
end
br_beam_model(model) = _with_shell(sh -> set_beam_model!(sh.sess.drive, _sym(model)))
# The 1ω (`h` = 1) or 2ω (2) power as a fraction of the 3ω power, as text (empty: off).
function br_unconverted_ratio(h, text)
    _with_shell() do sh
        dr = sh.sess.drive
        x = isempty(strip(string(text))) ? 0.0 :
            Controllers._parse_number(string(text), "the power ratio"; nonnegative = true)
        r = collect(dr.unconverted_ratio[])
        r[round(Int, h)] = x
        set_unconverted_ratio!(dr, r...)
    end
end
br_set_spot_radius(text) = _with_shell(sh -> (sh.sess.drive.spot_radius[] =
                                                  Controllers._parse_number(string(text), "the spot radius"; positive = true)))
function br_import_visrad_workspace(url)
    _with_shell() do sh
        # read and meshed on a worker (seconds for CAD); the design changes when it is ready
        import_visrad_workspace!(sh.sess, _url_to_path(String(url));
                                 done = () -> (select!(sh.sess.view, "target"); fit_all!(sh.scene)))
    end
end
br_export_visrad_workspace(url) = _with_shell(sh -> export_visrad_workspace(sh.sess, _url_to_path(String(url))))
br_export_results(url) = _with_shell(sh -> export_results(sh.sess, _url_to_path(String(url))))
br_export_stl(url) = _with_shell(sh -> export_target_stl(sh.sess, _url_to_path(String(url))))
# The retune's errors show in its dialog, where it was asked for.
function br_retune(table, out, by, cap, stay)
    _with_shell() do sh
        try
            retune_to_reference!(sh.sess, _url_to_path(String(table)), _url_to_path(String(out));
                                 by = _sym(by),
                                 cap = Controllers._parse_number(cap, "the power limit"; positive = true),
                                 stay = Controllers._parse_number(stay, "the weight on the reference powers";
                                                                  nonnegative = true))
        catch err
            sh.sess.drive.retune_status[] = "retune failed: " * Controllers._errmsg(err)
        end
    end
end
# Time panel
br_time_chart(key) = _with_shell(sh -> set_chart!(sh.sess.time, _sym(key)))
br_time_pulses(url) = _with_shell(sh -> load_pulses!(sh.sess, _url_to_path(String(url))))
br_time_shape() = _with_shell(sh -> use_pulse_shape!(sh.sess))
br_time_run() = _with_shell(sh -> start_history!(sh.sess))
br_time_cancel() = _with_shell(sh -> cancel_history!(sh.sess))
br_time_step(k) = _with_shell(sh -> set_time_step!(sh.sess, round(Int, k)))
br_time_now() = _with_shell(sh -> set_time_cursor!(sh.sess, NaN))
br_time_dante(on) = _with_shell(sh -> Controllers._assign!(sh.sess.time.dante, Bool(on)))
# Optimize
br_opt_set_by(by) = _with_shell(sh -> (by = _sym(by); by in BEAM_GROUPINGS &&
                                                      Controllers._assign!(sh.sess.optimize.by, by)))
br_opt_free(field, on) = _with_shell(sh -> set_free!(sh.sess.optimize, _sym(field), Bool(on)))
br_opt_edit(key, text) = _with_shell(sh -> edit_optimize!(sh.sess.optimize, _sym(key), string(text)))
br_opt_objective(k) = _with_shell(sh -> set_objective!(sh.sess.optimize, _sym(k)))
br_opt_hold(on) = _with_shell(sh -> (sh.sess.optimize.hold_drive[] = Bool(on)))
br_opt_keepout(on) = _with_shell(sh -> (sh.sess.optimize.keepout[] = Bool(on)))
br_opt_capsule(on) = _with_shell(sh -> (sh.sess.optimize.capsule_keepout[] = Bool(on)))
br_opt_start() = _with_shell(sh -> start_optimize!(sh.sess))
br_opt_cancel() = _with_shell(sh -> cancel_optimize!(sh.sess))
br_opt_accept() = _with_shell(sh -> accept_optimum!(sh.sess))
# Clearance
br_rv_role(name, role) = _with_shell(sh -> set_role!(sh.sess.review, _sym(name), _sym(role)))
br_rv_edit(key, text) = _with_shell(sh -> edit_review!(sh.sess.review, _sym(key), string(text)))
br_rv_wings(on) = _with_shell(sh -> (sh.sess.review.wings[] = Bool(on)))
br_rv_quality(q) = _with_shell(sh -> (q = _sym(q); haskey(REVIEW_QUALITY, q) &&
                                                   Controllers._assign!(sh.sess.review.quality, q)))
br_rv_start() = _with_shell(sh -> (start_clearance!(sh.sess);
                                   Controllers._assign!(sh.sess.view.beam_color, :clearance)))
br_rv_cancel() = _with_shell(sh -> cancel_clearance!(sh.sess))
br_rv_export(url) = _with_shell(sh -> export_clearance(sh.sess.review, _url_to_path(String(url))))
# Sweep
br_sw_subject(i) = _with_shell(sh -> set_sweep!(sh.sess.sweep; subject = sweep_subjects(sh.sess)[Int(i)+1]))
br_sw_frame(i) = _with_shell(sh -> set_sweep!(sh.sess.sweep; frame = _move_frame_choices(sh.sess)[Int(i)+1]))
br_sw_axis(i) = _with_shell(sh -> set_sweep!(sh.sess.sweep; axis = SWEEP_AXES[Int(i)+1]))
br_sw_edit(key, text) = _with_shell(sh -> edit_sweep!(sh.sess.sweep, _sym(key), string(text)))
br_sw_check(key, on) = _with_shell(sh -> set_sweep!(sh.sess.sweep; (_sym(key) => Bool(on),)...))
br_sw_object(i) = _with_shell(sh -> set_sweep!(sh.sess.sweep;
                                              los_object = [:none; surface_names(sh.sess.target.target[])][Int(i)+1]))
br_sw_quality(q) = _with_shell(sh -> set_sweep!(sh.sess.sweep; quality = _sym(q)))
br_sw_start() = _with_shell(sh -> (start_sweep!(sh.sess); request_tab("sweep")))
br_sw_cancel() = _with_shell(sh -> cancel_sweep!(sh.sess))
br_sw_select(k) = _with_shell(sh -> select_sweep_step!(sh.sess.sweep, round(Int, k)))
br_sw_apply(k) = _with_shell(sh -> apply_sweep_step!(sh.sess, round(Int, k)))
# TaLIS
br_tl_role(name, role) = _with_shell(sh -> set_talis_role!(sh.sess.talis, _sym(name), _sym(role)))
br_tl_edit(key, text) = _with_shell(sh -> edit_talis!(sh.sess.talis, _sym(key), string(text)))
br_tl_pivot(p) = _with_shell(sh -> set_talis_pivot!(sh.sess.talis, _sym(p)))
br_tl_omega1(on) = _with_shell(sh -> Controllers._assign!(sh.sess.talis.omega1, Bool(on)))
br_tl_start() = _with_shell(sh -> (start_talis!(sh.sess);
                                   Controllers._assign!(sh.sess.view.beam_color, :talis)))
br_tl_cancel() = _with_shell(sh -> cancel_talis!(sh.sess))
br_tl_export(url) = _with_shell(sh -> export_talis(sh.sess.talis, _url_to_path(String(url))))
br_cancel_all() = _with_shell(sh -> (cancel_optimize!(sh.sess); cancel_clearance!(sh.sess);
                                     cancel_talis!(sh.sess); cancel_sweep!(sh.sess)))
# Viewer
br_set_color(mode) = _with_shell(sh -> set_color!(sh.sess, _sym(mode)))
br_beam_style(s) = _with_shell(sh -> (s = _sym(s); s in BEAM_STYLES && Controllers._assign!(sh.sess.view.beams, s)))
br_beam_color(s) = _with_shell(sh -> (s = _sym(s); s in BEAM_COLORINGS &&
                                                   Controllers._assign!(sh.sess.view.beam_color, s)))
br_beam_length(text) = _with_shell(sh -> Controllers._assign!(sh.sess.view.beam_length,
                                                              isempty(strip(string(text))) ? 0.0 :
                                                              Controllers._parse_number(string(text), "the beam length";
                                                                                        nonnegative = true)))
br_show_off(on) = _with_shell(sh -> Controllers._assign!(sh.sess.view.show_off, Bool(on)))
br_set_cutaway(on) = _with_shell(sh -> Controllers._assign!(sh.sess.view.cutaway, Bool(on)))
# Cut parameters (the Cutaway button's right-click menu); any change turns the cut on.
br_cut_width(w) = _with_shell(sh -> (set_cut!(sh.sess.view; width = Float64(w), kind = :sector);
                                     Controllers._assign!(sh.sess.view.cutaway, true)))
br_cut_azimuth(a) = _with_shell(sh -> (set_cut!(sh.sess.view; azimuth = Float64(a), kind = :sector);
                                       Controllers._assign!(sh.sess.view.cutaway, true)))
# Center the cut on the side the 3D view looks from.
function br_cut_toward_view()
    _with_shell() do sh
        c = sh.scene.camera
        d = c.eye - c.lookat
        _, e1, e2 = Controllers._cut_basis(sh.sess.view.cut_axis[])
        x, y = dot(d, e1), dot(d, e2)
        hypot(x, y) > 0 || return
        set_cut!(sh.sess.view; azimuth = atand(y, x), kind = :sector)
        Controllers._assign!(sh.sess.view.cutaway, true)
    end
end
# The cut's kind ("sector" or "plane") and the sector's axis ("x", "y", "z").
br_cut_kind(k) = _with_shell(sh -> (set_cut!(sh.sess.view; kind = _sym(k));
                                    Controllers._assign!(sh.sess.view.cutaway, true)))
br_cut_axis(a) = _with_shell(sh -> (set_cut!(sh.sess.view; axis = _sym(a), kind = :sector);
                                    Controllers._assign!(sh.sess.view.cutaway, true)))
# The sector's bounds along its axis, as text (empty: unbounded).
function br_cut_range(lo, hi)
    _with_shell() do sh
        val(x, d, what) = isempty(strip(string(x))) ? d : _parse(x, what)
        set_cut!(sh.sess.view; range = (val(lo, -Inf, "the lower bound"), val(hi, Inf, "the upper bound")))
        Controllers._assign!(sh.sess.view.cutaway, true)
    end
end
# The cut plane: normal along an axis ("x", "y", "z", "-x", …) or toward the view
# ("view": the side facing the eye is removed), through the plane's point at offset 0;
# or the other side kept ("flip": the same plane, its normal turned over).
function br_cut_plane(which)
    _with_shell() do sh
        v = sh.sess.view
        w = String(which)
        if w == "flip"
            flip_cut!(v)
        else
            n = if w == "view"
                c = sh.scene.camera
                Briareus.Vec3((c.eye - c.lookat)...)
            else
                (startswith(w, "-") ? -1 : 1) * Controllers._axis_vector(Symbol(lstrip(w, '-')))
            end
            set_cut_plane!(v; normal = n, offset = 0.0)
        end
        Controllers._assign!(v.cutaway, true)
    end
end
br_cut_offset(x) = _with_shell(sh -> (set_cut_plane!(sh.sess.view; offset = Float64(x));
                                      Controllers._assign!(sh.sess.view.cutaway, true)))
# The plane's normal typed as port angles θ, φ (degrees) or as a vector "x, y, z".
br_cut_angles(theta, phi) =
    _with_shell(sh -> (set_cut_plane!(sh.sess.view; angles = (_parse(theta, "the polar angle θ"),
                                                               _parse(phi, "the azimuth φ")));
                       Controllers._assign!(sh.sess.view.cutaway, true)))
br_cut_normal(text) = _with_shell(sh -> (set_cut_plane!(sh.sess.view; normal = _parse_vector(text, "the normal"));
                                         Controllers._assign!(sh.sess.view.cutaway, true)))
# The plane's point, typed (mm), or the target center (TCC).
br_cut_point(x, y, z) =
    _with_shell(sh -> (set_cut_plane!(sh.sess.view; point = (_parse(x, "x"), _parse(y, "y"), _parse(z, "z")));
                       Controllers._assign!(sh.sess.view.cutaway, true)))
function _parse_vector(text, what)
    xs = split(strip(string(text), ['(', ')', '[', ']', ' ']), r"[,\s]+"; keepempty = false)
    length(xs) == 3 || throw(ArgumentError("$what is three numbers, x, y, z"))
    return Tuple(_parse(x, what) for x in xs)
end
# The plane through the selected component's centroid, with its normal if it is flat.
function br_cut_through_selected()
    _with_shell() do sh
        name = selected_component(sh.sess)
        name === nothing && throw(ArgumentError("select a component, window or part to cut through"))
        Controllers.await_target!(sh.sess.target)
        v = sh.sess.view
        cut_through_component!(v, sh.sess.target.target[], name)
        Controllers._assign!(v.cutaway, true)
        sh.sess.status[] = "cut plane through $(component_label(name))" *
                           (v.cut_line[] == (0.0, 0.0, 0.0) ? "" : ", containing the line of sight")
    end
end
# The plane containing the line of sight at port angles θ, φ from the target center or
# (`from_selected`) the selected component's centroid, turned `rotation` degrees about it.
function br_cut_line(theta, phi, from_selected, rotation)
    _with_shell() do sh
        v = sh.sess.view
        origin = (0.0, 0.0, 0.0)
        if Bool(from_selected)
            name = selected_component(sh.sess)
            name === nothing && throw(ArgumentError("select a component, window or part for the line's origin"))
            Controllers.await_target!(sh.sess.target)
            origin = Tuple(component_plane(sh.sess.target.target[], name)[1])
        end
        set_cut_line!(v; angles = (_parse(theta, "the polar angle θ"), _parse(phi, "the azimuth φ")),
                      origin, rotation = Float64(rotation))
        Controllers._assign!(v.cutaway, true)
    end
end
br_cut_rotation(a) = _with_shell(sh -> (set_cut_line!(sh.sess.view; rotation = Float64(a));
                                        Controllers._assign!(sh.sess.view.cutaway, true)))
function br_cut_component(name, on)
    _with_shell() do sh
        v = sh.sess.view
        s = copy(v.cut_components[])
        Bool(on) ? push!(s, Symbol(String(name))) : delete!(s, Symbol(String(name)))
        set_cut!(v; components = s)
        Controllers._assign!(v.cutaway, true)
    end
end
# Camera actions bring the 3D view back to the canvas (unless the map sits beside it).
_show_scene!(sh) = sh.sess.drive.side[] || set_view!(sh.sess.drive, :scene)
br_standard_view(name) = _with_shell(sh -> (_show_scene!(sh); standard_view!(sh.scene, _sym(name))))
br_fit_all() = _with_shell(sh -> sh.figure === sh.map.fig ? (reset_limits!(sh.map.ax); sh.map.dirty[] = true) :
                                 fit_all!(sh.scene))
br_fit_selection() = _with_shell(sh -> fit_selection!(sh.scene, sh.sess))
function br_look_from_port(theta, phi)
    _with_shell() do sh
        θ, φ = _parse(theta, "the polar angle θ"), _parse(phi, "the azimuth φ")
        _show_scene!(sh)
        look_from_port!(sh.scene, sh.sess, θ, φ)
        sh.sess.status[] = @sprintf("looking from the port at θ = %g°, φ = %g°", θ, φ)
    end
end
# Port angles of beam `i` of the table (0-based, as the dialog lists them), as text.
function br_port_angles(i)
    sh = _SHELL[]
    sh === nothing && return ["", ""]
    bs = sh.sess.beams.beams[]
    0 <= i < length(bs) || return ["", ""]
    θ, φ = port_angles(bs[i+1])
    return [@sprintf("%.3f", θ), @sprintf("%.3f", mod(φ, 360))]
end
br_look_along_beam() = _with_shell(sh -> (_show_scene!(sh); look_along_beam!(sh.scene, sh.sess)))
br_projection(ortho) = _with_shell(sh -> (set_projection!(sh.scene, Bool(ortho)); _set!(sh, "ortho", Bool(ortho))))
br_turntable(on) = _with_shell(sh -> (set_turntable!(sh.scene, Bool(on)); _set!(sh, "turntable", Bool(on))))
# The canvas lost the keyboard (or the window lost focus): forget held keys and drags.
br_canvas_focus(on) = (sh = _SHELL[]; sh === nothing || Bool(on) || release_input!(sh.scene); nothing)

function _register_qml_functions()
    @qmlfunction br_tick br_grab_started br_select br_select_add br_toggle_expand br_tree_on
    @qmlfunction br_select_range br_node_flags br_node_show br_node_include
    @qmlfunction br_undo br_redo br_open_session br_save_session br_reload_cad_sources br_keep_cad_copies
    @qmlfunction br_relink_part br_part_from_source
    @qmlfunction br_set_cluster br_rename_cluster br_split_quad br_merge_quad br_window_round
    @qmlfunction br_set_facility br_aim_frame br_aim_at_frame
    @qmlfunction br_set_target_kind br_set_resolution br_edit_mesh br_edit_target br_leh_closure
    @qmlfunction br_add_window br_remove_window br_edit_window br_window_open
    @qmlfunction br_add_part br_add_stl br_edit_part br_part_flag br_part_frame br_part_coords
    @qmlfunction br_move_frame br_move_axis br_move_step br_nudge_part br_sweep_part
    @qmlfunction br_sw_subject br_sw_frame br_sw_axis br_sw_edit br_sw_check br_sw_object
    @qmlfunction br_sw_quality br_sw_start br_sw_cancel br_sw_select br_sw_apply
    @qmlfunction br_part_role br_export_stl br_cavity br_cavity_point br_part_film
    @qmlfunction br_remove_selected
    @qmlfunction br_load_preset br_set_grouping br_edit_beam br_set_beam_on br_edit_node
    @qmlfunction br_import_beams br_export_beams br_import_sst br_import_summary br_export_sst
    @qmlfunction br_export_summary br_build br_set_view br_side_view br_xray_image br_close_image br_set_spot_radius br_unconverted_ratio br_beam_model br_export_results
    @qmlfunction br_import_visrad_workspace br_export_visrad_workspace br_retune
    @qmlfunction br_time_chart br_time_pulses br_time_shape br_time_run br_time_cancel br_time_step
    @qmlfunction br_time_now br_time_dante
    @qmlfunction br_opt_set_by br_opt_free br_opt_edit br_opt_objective br_opt_hold
    @qmlfunction br_opt_keepout br_opt_capsule br_opt_start br_opt_cancel br_opt_accept
    @qmlfunction br_rv_role br_rv_edit br_rv_wings br_rv_quality br_rv_start br_rv_cancel
    @qmlfunction br_rv_export br_tl_role br_tl_edit br_tl_pivot br_tl_omega1 br_tl_start
    @qmlfunction br_tl_cancel br_tl_export br_cancel_all
    @qmlfunction br_set_color br_beam_style br_beam_color br_beam_length br_set_cutaway br_show_off
    @qmlfunction br_cut_width br_cut_azimuth br_cut_toward_view br_cut_component
    @qmlfunction br_cut_kind br_cut_axis br_cut_range br_cut_plane br_cut_offset
    @qmlfunction br_cut_angles br_cut_normal br_cut_point br_cut_through_selected br_cut_line
    @qmlfunction br_cut_rotation
    @qmlfunction br_set_display br_beam_dot
    @qmlfunction br_standard_view br_fit_all br_fit_selection br_look_along_beam
    @qmlfunction br_look_from_port br_port_angles
    @qmlfunction br_projection br_turntable br_canvas_focus
    return
end

# ---------------------------------------------------------------- the window

"""
    briareus(sess = DesignSession(); kind = :hohlraum, resolution = :standard,
             preset = nothing, beams = nothing, build = false, shown = nothing) -> DesignSession

Open the Briareus window and return its session when the window closes. The window is
a design tool in the manner of VisRad, not a sequence of steps:

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
  cutaway (a pie slice or a section plane, clipped exactly), and an orientation triad.
  The capsule θ–φ (or foil) map shares the canvas.
- **Results** (bottom): the beam table, the drive (modes, power accounting), time
  (pulses, the drive history and the retune plotted against time, with a time cursor
  that colors the 3D view by that step), the optimization, the clearance check and the
  TaLIS rules.
- **Menus and toolbar**: new targets, parts and windows, beam presets and imports
  (CSV, VisRad SST and Laser Beam Summary) and exports, views, and the analyses.

The session starts from `sess`, or a fresh `DesignSession(; kind, resolution,
preset, beams)`; `build = true` starts computing the drive at once, and `shown()`, when
given, runs once on the window's first tick (the `briareus` command closes its splash
then). The call blocks while the window is open (Qt runs its event loop on this thread); start Julia with
several threads (`julia -t auto`) so background jobs leave the window responsive.
Closing the window cancels an optimization, clearance or TaLIS check in progress.

The window throws an `ArgumentError` when this Julia session already has a GLMakie
screen (a GLMakie window or an offscreen `colorbuffer`): GLMakie's windows and the Qt
window cannot share a process.
"""
function briareus(sess::Union{Nothing,DesignSession} = nothing; kind::Symbol = :hohlraum,
                  resolution::Symbol = :standard, preset = nothing, beams = nothing,
                  build::Bool = false, shown = nothing)
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
        # a target left building: the session returned holds the design's (parts left
        # remeshing are meshed when the drive asks)
        ts = sess.target
        ts.building[] && isempty(ts.pending[]) && await_target!(ts)
        _release_qml_screens!()
    end
    return sess
end
