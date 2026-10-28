# The placement sweep: move a part (and what is placed in its frame, and the beams aimed
# in those frames) along an axis of a named frame through a range of shifts, and check
# each step on a worker (`placement_sweep`): beam clearance, the backlighter's line of
# sight, the power the foil captures, and optionally a re-optimized pointing. A step of
# the result can be applied to the design as one edit.

"""
Placement-sweep quality presets, on top of the clearance check's sampling of the same
name ([`REVIEW_QUALITY`](@ref)): `:quick` (line of sight from the foil's vertices
through 2 rings of object points; the captured fraction on a 0.6 mm deposition grid) and
`:standard` (4 sub-points per foil triangle and 4 rings; the deposition's default 0.3 mm
grid). Engineering defaults for a sweep's many steps, not converged settings.
"""
const SWEEP_QUALITY = Dict(:quick => (; level = 0, nrings = 2, resolution = 0.6),
                           :standard => (; level = 1, nrings = 4, resolution = 0.3))

"Axes a sweep moves along: those of the chosen frame."
const SWEEP_AXES = (:x, :y, :z)

"The most steps a sweep takes (each rebuilds the target and runs every chosen check)."
const SWEEP_MAX_STEPS = 201

"""
    SweepState()

State of the placement sweep. Settings: `subject` (the part moved, with the parts
placed in its frame and the beams aimed in those frames; or `:backlighter`, a backlighter
target's foil frame; `:none` until chosen, see [`sweep_subject`](@ref)), `frame` and
`axis` (the move is along `axis`, one of [`SWEEP_AXES`](@ref), of the named frame,
[`part_frames`](@ref)), `start`, `stop`, `step` (the shifts, mm; the sweep translates,
it does not turn), the checks `clearance` (the clearance check's components, envelope
and errors), `los` (the line of sight at polar angle `los_theta` and azimuth `los_phi`,
degrees, from the backlighter foil through `los_object`; `NaN`, the default, for the
direction from the foil's center through the object's, see [`sweep_los`](@ref)),
`capture` (the least fraction
of a beam's power the foil absorbs) and `optimize` (the pointing re-optimized at each
step with the optimization's settings, warm-started from the previous step when the
aims are in frames that move with the part), and `quality` (a [`SWEEP_QUALITY`](@ref)
key). Job: `running`, `progress` `(steps done, total)`, `started`, `status`, `error`,
`result` (a core `SweepResult`, partial with `canceled` after a cancel), `selected` (the
row shown on the chart, 0 for none) and `inputs` (what the result was computed for).
"""
struct SweepState
    subject::Observable{Symbol}
    frame::Observable{Symbol}
    axis::Observable{Symbol}
    start::Observable{Float64}
    stop::Observable{Float64}
    step::Observable{Float64}
    clearance::Observable{Bool}
    los::Observable{Bool}
    los_theta::Observable{Float64}
    los_phi::Observable{Float64}
    los_object::Observable{Symbol}
    capture::Observable{Bool}
    optimize::Observable{Bool}
    quality::Observable{Symbol}
    running::Observable{Bool}
    progress::Observable{Tuple{Int,Int}}
    started::Observable{Float64}
    status::Observable{String}
    error::Observable{String}
    result::Observable{Union{Nothing,SweepResult}}
    selected::Observable{Int}
    inputs::Base.RefValue{Any}
    cancel::Threads.Atomic{Bool}
end

SweepState() = SweepState(Observable(:none), Observable(:none), Observable(:z), Observable(-1.0),
                          Observable(1.0), Observable(0.5), Observable(true), Observable(false),
                          Observable(NaN), Observable(NaN), Observable(:none), Observable(false),
                          Observable(false), Observable(:quick), Observable(false), Observable((0, 0)),
                          Observable(0.0), Observable(""), Observable(""),
                          Observable{Union{Nothing,SweepResult}}(nothing), Observable(0), Ref{Any}(nothing),
                          Threads.Atomic{Bool}(false))

"""
    set_sweep!(sw::SweepState; subject, frame, axis, clearance, los, los_object, capture,
               optimize, quality) -> sw

Set any of the sweep's choices (see [`SweepState`](@ref)); the numbers are edited with
[`edit_sweep!`](@ref).
"""
function set_sweep!(sw::SweepState; subject = nothing, frame = nothing, axis = nothing,
                    clearance = nothing, los = nothing, los_object = nothing, capture = nothing,
                    optimize = nothing, quality = nothing)
    axis === nothing || axis in SWEEP_AXES || throw(ArgumentError("the axis is one of $SWEEP_AXES"))
    quality === nothing || haskey(SWEEP_QUALITY, quality) ||
        throw(ArgumentError("quality is one of $(Tuple(keys(SWEEP_QUALITY)))"))
    for (obs, v) in ((sw.subject, subject), (sw.frame, frame), (sw.axis, axis), (sw.clearance, clearance),
                     (sw.los, los), (sw.los_object, los_object), (sw.capture, capture),
                     (sw.optimize, optimize), (sw.quality, quality))
        v === nothing || _assign!(obs, v)
    end
    return sw
end

"""
    edit_sweep!(sw::SweepState, key, text) -> Bool

Edit `:start`, `:stop`, `:step` (mm) or the line of sight's `:los_theta`, `:los_phi`
(degrees; empty text: from the foil through the object, see [`sweep_los`](@ref)) from
text. A rejected value sets `error` and returns `false`.
"""
function edit_sweep!(sw::SweepState, key::Symbol, text)
    try
        key in (:start, :stop, :step, :los_theta, :los_phi) ||
            throw(ArgumentError("unknown sweep setting :$key"))
        v = key in (:los_theta, :los_phi) && isempty(strip(string(text))) ? NaN :
            _parse_number(text, key === :step ? "the step" : "the $(replace(String(key), '_' => ' '))";
                          positive = key === :step)
        _assign!(getfield(sw, key), v)
        sw.error[] = ""
        return true
    catch err
        err isa ArgumentError || rethrow()
        sw.error[] = _errmsg(err)
        return false
    end
end

"""
    sweep_values(sw::SweepState) -> Vector{Float64}

The shifts of the sweep: `start`, `start + step`, … up to `stop` (within 1e-9 of a step,
so a range that divides evenly ends on `stop`). At most [`SWEEP_MAX_STEPS`](@ref).
"""
function sweep_values(sw::SweepState)
    a, b, h = sw.start[], sw.stop[], sw.step[]
    h > 0 || throw(ArgumentError("the step must be positive"))
    b >= a || throw(ArgumentError("the sweep must stop at or after its start"))
    n = floor(Int, (b - a) / h + 1e-9) + 1
    n <= SWEEP_MAX_STEPS ||
        throw(ArgumentError("$n steps: at most $SWEEP_MAX_STEPS (a larger step, or a shorter range)"))
    return [a + (k - 1) * h for k in 1:n]
end

"""
    sweep_subjects(sess) -> Vector{Symbol}

What a sweep can move: the backlighter foil's frame (`:backlighter`, a foil target's)
and every part.
"""
sweep_subjects(sess) =
    [sess.target.kind[] === :foil ? [:backlighter] : Symbol[]; [p.name for p in sess.target.parts[]]]

"""
    sweep_subject(sess) -> Symbol

The subject the sweep moves: the chosen one when the design still has it, otherwise the
first of [`sweep_subjects`](@ref) (the foil of a backlighter target), or `:none`.
"""
function sweep_subject(sess)
    subs = sweep_subjects(sess)
    s = sess.sweep.subject[]
    s in subs && return s
    return isempty(subs) ? :none : subs[1]
end

"""
    sweep_frame(sess) -> Symbol

The frame whose axis the sweep moves along: the chosen one when the design still has it,
otherwise (`:none`, the default) the subject's own frame (a part's, the foil's), or
`:tcc`.
"""
function sweep_frame(sess)
    fr = part_frames(sess.target)
    f = sess.sweep.frame[]
    f !== :none && haskey(fr, f) && return f
    s = sweep_subject(sess)
    return haskey(fr, s) ? s : :tcc
end

"""
    sweep_los(sess) -> (θ, φ)

The direction of the sweep's line of sight (degrees, from the object toward the
detector): `los_theta`, `los_phi` when both are set, otherwise from the backlighter
foil's centroid through the object's ([`sweep_los_object`](@ref)), as the foil
backlights it, at the sweep's start.
"""
function sweep_los(sess)
    sw = sess.sweep
    θ, φ = sw.los_theta[], sw.los_phi[]
    isfinite(θ) && isfinite(φ) && return (θ, φ)
    t = sess.target.target[]
    o = sweep_los_object(sess)
    (o === :none || !(FOIL_NAME in component_names(t))) && return (90.0, 0.0)
    d = part_centroid(t, o) - part_centroid(t, FOIL_NAME)
    return (acosd(clamp(d[3] / norm(d), -1, 1)), mod(atand(d[2], d[1]), 360))
end

"""
    sweep_los_object(sess) -> Symbol

The object the sweep's line of sight passes: the chosen one when the target has it,
otherwise a foil target's main sphere (`:main`) or a capsule, or `:none`.
"""
function sweep_los_object(sess)
    names = component_names(sess.target.target[])
    o = sess.sweep.los_object[]
    o in names && return o
    for d in (:main, :capsule)
        d in names && return d
    end
    return :none
end

# The frames that move with the subject: its own and every part's placed in one of them.
function _moved_frames(parts, subject::Symbol)
    moved = Set{Symbol}((subject,))
    grew = true
    while grew
        grew = false
        for p in parts
            p.frame in moved && !(p.name in moved) && (push!(moved, p.name); grew = true)
        end
    end
    return moved
end

# The target components that move with the subject: the parts of the moved frames (a
# VisRad part's components carry its name as a prefix) and, for `:backlighter`, the foil.
function _moved_components(t::Target, parts, frames::Set{Symbol})
    moved = Set{Symbol}()
    :backlighter in frames && push!(moved, FOIL_NAME)
    for p in parts
        p.name in frames || continue
        for c in t.components
            (c.name === p.name || p.kind === :visrad && startswith(String(c.name), String(p.name))) &&
                push!(moved, c.name)
        end
    end
    return moved
end

# The design of a sweep step: the moved components translated by `value * dir`, the
# followed beams' aims with them, and (when re-optimizing) the step's parameterization:
# a `FramePointing` in the beams' frames (the followed ones moved), else a `GroupPointing`.
struct _SweepDesign
    target::Target
    moved::Set{Symbol}
    dir::Vec3
    beams::Vector{Beam{Float64}}
    follow::Vector{Bool}
    frames::Vector{Frame}
    pointing::Any                      # nothing, or (; key, fields, by)
end

function (d::_SweepDesign)(value)
    Δ = value * d.dir
    t = d.target
    comps = Component[c.name in d.moved ? translate(c, Δ) : c for c in t.components]
    target = Target(comps; closed = t.closed, visibility = t.visibility)
    beams = Beam{Float64}[m ? remake(b; aim = b.aim + Δ) : b for (b, m) in zip(d.beams, d.follow)]
    pointing = if d.pointing === nothing
        nothing
    else
        frames = Frame[m ? Frame(F.R, F.origin + Δ) : F for (F, m) in zip(d.frames, d.follow)]
        _parameterization(beams, frames, d.pointing.key, d.pointing.by, d.pointing.fields)
    end
    return (; target, beams, pointing)
end

# The optimization's objective, made for each step's target (the foil's center, the
# keep-outs and the held drive are the step's): made once per target, under a lock, as
# gradient chunks evaluate it concurrently.
struct _StepObjective
    cfg::NamedTuple
    spot_radius::Float64
    made::IdDict{Any,Any}
    lock::ReentrantLock
end
_StepObjective(cfg, spot_radius) = _StepObjective(cfg, Float64(spot_radius), IdDict{Any,Any}(), ReentrantLock())

function (o::_StepObjective)(ev)
    s = ev.scenario
    f = lock(o.lock) do
        get!(() -> _make_objective(o.cfg, s, nothing, o.spot_radius), o.made, s)
    end
    return f(ev)
end

# Everything a sweep reads, captured on the GUI thread (and compared for staleness).
function _sweep_inputs(sess)
    sw, ts, bt = sess.sweep, sess.target, sess.beams
    subject = sweep_subject(sess)
    subject === :none && throw(ArgumentError("add a part (or use a backlighter target) to sweep"))
    frame = sweep_frame(sess)
    any((sw.clearance[], sw.los[], sw.capture[], sw.optimize[])) ||
        throw(ArgumentError("choose at least one check"))
    values = sweep_values(sw)
    names = component_names(ts.target[])
    (sw.los[] || sw.capture[]) && !(FOIL_NAME in names) &&
        throw(ArgumentError("the line of sight and the captured power need the backlighter foil " *
                            "(a foil target's, or a part with the role Backlighter foil)"))
    obj = sweep_los_object(sess)
    sw.los[] && obj === :none && throw(ArgumentError("choose the object the line of sight passes"))
    idx = findall(bt.enabled[])
    isempty(idx) && throw(ArgumentError("no beams are enabled"))
    rv = sess.review
    if sw.clearance[]
        isempty(rv.targets[]) && throw(ArgumentError("choose the components the beams land on (Clearance)"))
        isempty(rv.hardware[]) && throw(ArgumentError("choose hardware to check (Clearance)"))
    end
    q = sw.quality[]
    st = sess.optimize
    opt = sw.optimize[] ? (; by = st.by[], fields = Tuple(st.fields[]), deltas = copy(st.deltas[]),
                           cfg = _objective_settings(st), iterations = st.iterations[],
                           evaluations = st.max_evaluations[], spot_radius = sess.drive.spot_radius[]) :
          nothing
    if opt !== nothing
        isempty(opt.fields) && throw(ArgumentError("free at least one field (Optimization)"))
        _check_grouping(bt.beams[][idx], opt.by)
    end
    return (; target = ts.target[], frames = part_frames(ts), parts = ts.parts[], subject, frame,
            axis = sw.axis[], values, idx, beams = bt.beams[][idx], aim_frame = bt.aim_frame[][idx],
            clusters = bt.cluster[][idx],
            clearance = sw.clearance[] ? merge(clearance_kwargs(rv), REVIEW_QUALITY[q]) : nothing,
            los = sw.los[] ? sweep_los(sess) : nothing, los_object = obj,
            capture = sw.capture[], quality = q, optimize = opt, profile = _deposition_profile(sess.drive))
end

"""
    start_sweep!(sess::DesignSession; spawn = true)

Run the placement sweep with the current settings ([`SweepState`](@ref)) on the enabled
beams: at each shift the subject ([`sweep_subject`](@ref)), the parts placed in its
frame and the beams aimed in those frames move by the shift along the frame's axis
([`sweep_frame`](@ref)), and `placement_sweep` checks the step, on a worker with
`spawn`. Progress goes to `progress` after each step; [`cancel_sweep!`](@ref) stops after
the step in flight and keeps the steps done. Settings the design cannot serve set
`status` and start nothing.
"""
function start_sweep!(sess; spawn::Bool = true)
    sw = sess.sweep
    sw.running[] && return sess
    await_target!(sess.target)
    inp = try
        _sweep_inputs(sess)
    catch err
        err isa ArgumentError || rethrow()
        sw.status[] = _errmsg(err)
        return sess
    end
    t = inp.target
    moved_frames = _moved_frames(inp.parts, inp.subject)
    moved = _moved_components(t, inp.parts, moved_frames)
    dir = Vec3(to_chamber(inp.frames[inp.frame], _axis_vector3(inp.axis); direction = true))
    follow = Bool[f in moved_frames for f in inp.aim_frame]
    bframes = Frame[inp.frames[f] for f in inp.aim_frame]
    opt = inp.optimize
    key = opt === nothing ? nothing :
          opt.by === :group ? (d = Dict(b.name => c for (b, c) in zip(inp.beams, inp.clusters)); b -> d[b.name]) :
          opt.by
    design = _SweepDesign(t, moved, dir, inp.beams, follow, bframes,
                          opt === nothing ? nothing : (; key, by = opt.by, fields = opt.fields))
    q = SWEEP_QUALITY[inp.quality]
    los_kw = inp.los === nothing ? (;) :
             (; foil = FOIL_NAME, object = inp.los_object, center = part_centroid(t, inp.los_object),
              level = q.level, nrings = q.nrings)
    optimize = if opt === nothing
        nothing
    else
        frame_fields = _frame_fields(opt.fields)
        d = opt.deltas
        scenario_kw = opt.cfg.kind === :irradiance ? (; deposition_kw = (; inp.profile)) :
                      (; observe = :capsule, radiation = true, deposition_kw = (; inp.profile))
        bounds_kw = frame_fields ?
                    (; u = d[:u], v = d[:v], w = d[:w], defocus = d[:defocus], power = d[:power]) :
                    (; z = d[:z], r = d[:r], phi = d[:phi], defocus = d[:defocus], power = d[:power])
        # local aims in frames that move with the part keep their meaning from step to step
        (; objective = _StepObjective(opt.cfg, opt.spot_radius), scenario_kw, bounds_kw,
         warm = frame_fields, iterations = opt.iterations, evaluations = opt.evaluations)
    end
    capture = inp.capture ? FOIL_NAME : nothing
    capture_kw = (; resolution = q.resolution, inp.profile)
    deliver = sess.deliver[]
    sw.cancel[] = false
    sw.inputs[] = (; inp..., dir, follow)
    sw.progress[] = (0, length(inp.values))
    sw.started[] = time()
    sw.status[] = "sweeping…"
    sw.error[] = ""
    sw.running[] = true
    job = function ()
        progress = function (info)
            done, total = info.done, info.total
            deliver(() -> (sw.progress[] = (done, total)))
            return sw.cancel[]
        end
        out = _try_job(() -> placement_sweep(design, inp.values; clearance_kw = inp.clearance,
                                             los = inp.los, los_kw, capture, capture_kw, optimize,
                                             callback = progress))
        deliver(() -> _finish_sweep!(sw, out))
    end
    _run_job(job, spawn)
    return sess
end

function _finish_sweep!(sw::SweepState, out)
    if out.err !== nothing
        sw.status[] = "failed: " * _errmsg(out.err)
    else
        r = out.value
        sw.result[] = isempty(r.rows) ? nothing : r
        _assign!(sw.selected, 0)
        done, total = sw.progress[]
        sw.status[] = r.canceled ? "canceled after $done of $total steps" : "done: " * sweep_brief(sw)
    end
    sw.running[] = false
    return sw
end

"""
    cancel_sweep!(sess::DesignSession)

Stop the sweep after the step in flight (`placement_sweep`'s callback returns `true`);
the steps done are kept, the result marked `canceled`.
"""
cancel_sweep!(sess) = (sess.sweep.cancel[] = true; sess)

"""
    sweep_stale(sess::DesignSession) -> Bool

Whether the design or a setting the last sweep read changed since: the target, the
enabled beams and their frames, the sweep's settings, or the clearance check's and the
optimization's settings it used. `false` before the first sweep.
"""
function sweep_stale(sess)
    sw = sess.sweep
    inp = sw.inputs[]
    (inp === nothing || sw.running[]) && return false
    now = try
        _sweep_inputs(sess)
    catch err
        err isa ArgumentError || rethrow()
        return true
    end
    return now.target !== inp.target || now.beams != inp.beams || now.aim_frame != inp.aim_frame ||
           now.frames != inp.frames || now.subject !== inp.subject || now.frame !== inp.frame ||
           now.axis !== inp.axis || now.values != inp.values || now.clearance != inp.clearance ||
           now.los != inp.los || now.los_object !== inp.los_object || now.capture != inp.capture ||
           now.quality !== inp.quality || !isequal(now.optimize, inp.optimize) ||
           _profile_state(now.profile) != _profile_state(inp.profile)
end

"""
    apply_sweep_step!(sess::DesignSession, k) -> Bool

Put step `k` of the sweep into the design, as one edit (undone as one): the subject moves
by the step's shift along the sweep's axis from where it was when the sweep ran (a part's
center; the backlighter foil's distance, which moves it along chamber `x` only), the
parts placed in its frame and the beams aimed in those frames with it, and with
re-optimization on, the enabled beams take the step's optimized pointing. Refused (with
`status` set, returning `false`) when the design changed since the sweep.
"""
function apply_sweep_step!(sess, k::Integer)
    sw = sess.sweep
    r = sw.result[]
    r === nothing && throw(ArgumentError("run the sweep first"))
    1 <= k <= length(r.rows) || throw(ArgumentError("no sweep step $k"))
    if sweep_stale(sess)
        sw.status[] = "the design changed since the sweep: run it again to apply a step"
        return false
    end
    inp = sw.inputs[]
    ts = sess.target
    row = r.rows[k]
    Δ = row.value * inp.dir
    ok = if inp.subject === :backlighter && ts.kind[] === :foil
        abs(abs(inp.dir[1]) - 1) < 1e-12 ||
            throw(ArgumentError("the backlighter foil moves along chamber x only (its distance from TCC)"))
        p = copy(ts.params[])
        p[:distance] += Δ[1]
        _try_build!(ts; params = p)
    else
        i = part_index(ts, inp.subject)
        i === nothing && throw(ArgumentError("no part $(inp.subject)"))
        q = ts.parts[][i]
        c = part_center(q) + to_frame(inp.frames[q.frame], Δ; direction = true)
        vals = copy(q.values)
        vals[:x], vals[:y], vals[:z] = c
        ps = copy(ts.parts[])
        ps[i] = _with(q; values = vals)
        _try_build!(ts; parts = ps)
    end
    ok || (sw.status[] = "the step was rejected: " * ts.error[]; return false)
    if inp.optimize !== nothing
        bt = sess.beams
        bs = copy(bt.beams[])
        bs[inp.idx] = row.beams
        _set_aims!(bt, bs)
    end
    sw.status[] = @sprintf("applied the step at %+.4g mm", row.value)
    return true
end

"""
    sweep_brief(sw::SweepState) -> String

The sweep's outcome in a line: steps, the least clearance and where, the steps with a
blocked line of sight.
"""
function sweep_brief(sw::SweepState)
    r = sw.result[]
    r === nothing && return "not run"
    rows = r.rows
    parts = String["$(length(rows)) of $(length(r.values)) steps"]
    cl = filter(w -> isfinite(w.min_clearance), rows)
    isempty(cl) || (w = argmin(w -> w.min_clearance, cl);
                    push!(parts, @sprintf("least clearance %.3f mm at %+.4g", w.min_clearance, w.value)))
    nb = count(w -> w.los_clear === false, rows)
    any(w -> w.los_clear !== missing, rows) && push!(parts, nb == 0 ? "line of sight clear" :
                                                            "line of sight blocked at $nb")
    return join(parts, " · ")
end

"""
    sweep_rows(sw::SweepState) -> Vector{NamedTuple}

The result's rows for the table: `value` (mm), `clearance` (mm), `worst` (beam and
part), `fails`, `reviews`, `los` (`"clear"`, `"blocked"` or `""`), `blocked` and
`covered` (fractions), `captured` (fraction), `objective` (before and after), each a
number or `NaN` when the check did not run.
"""
function sweep_rows(sw::SweepState)
    r = sw.result[]
    r === nothing && return NamedTuple[]
    return [(; value = w.value, clearance = w.min_clearance,
             worst = w.worst_beam == "" ? "" : string(w.worst_beam, " / ", w.worst_component),
             fails = w.fails, reviews = w.reviews,
             los = w.los_clear === missing ? "" : w.los_clear ? "clear" : "blocked",
             blocked = w.blocked_fraction, covered = w.covered, captured = w.min_captured,
             objective_before = w.objective_before, objective_after = w.objective_after)
            for w in r.rows]
end

"""
    sweep_chart(sess) -> NamedTuple

What the sweep panel plots, in the form of [`time_chart`](@ref) (`chart = :sweep`): the
least clearance (mm, left axis) and the line of sight's blocked fraction (%, right
axis) against the shift, with the selected step as the cursor.
"""
function sweep_chart(sess)
    sw = sess.sweep
    r = sw.result[]
    inp = sw.inputs[]
    xlabel = inp === nothing ? "shift (mm)" : "shift along $(inp.axis) of $(inp.frame) (mm)"
    r === nothing && return (; chart = :sweep, title = "Placement sweep", xlabel,
                             left = (; label = "", series = NamedTuple[]), right = nothing,
                             legend = NamedTuple[], cursor = NaN,
                             empty = "Set up the sweep in Properties (select it in the model tree) and run it")
    x = [w.value for w in r.rows]
    left = NamedTuple[]
    any(w -> isfinite(w.min_clearance), r.rows) &&
        push!(left, _series("least clearance", x, [w.min_clearance for w in r.rows], _INK))
    right = NamedTuple[]
    any(w -> isfinite(w.blocked_fraction), r.rows) &&
        push!(right, _series("line of sight blocked", x, [100 * w.blocked_fraction for w in r.rows],
                             "#D55E00"; style = :dash))
    k = sw.selected[]
    cursor = 1 <= k <= length(x) ? x[k] : NaN
    stale = sweep_stale(sess) ? " (out of date: run it again)" : ""
    return (; chart = :sweep, title = "Placement sweep: " * sweep_brief(sw) * stale, xlabel,
            left = (; label = "Least clearance (mm)", series = left),
            right = isempty(right) ? nothing : (; label = "Blocked (%)", series = right),
            legend = [_key(s) for s in vcat(left, right)], cursor,
            empty = isempty(left) && isempty(right) ? "No clearance or line-of-sight check in this sweep" : "")
end

"""
    select_sweep_step!(sw::SweepState, k) -> sw

Show step `k` (0 for none) as the sweep chart's cursor and the table's selected row.
"""
select_sweep_step!(sw::SweepState, k::Integer) = (_assign!(sw.selected, Int(k)); sw)

# The sweep's settings as plain data (sessions, undo), and back. The line of sight's
# angles are left out while they follow the foil and the object (`NaN`).
function _sweep_state(sw::SweepState)
    d = Dict{String,Any}("subject" => String(sw.subject[]), "frame" => String(sw.frame[]),
                         "axis" => String(sw.axis[]), "start" => sw.start[], "stop" => sw.stop[],
                         "step" => sw.step[], "clearance" => sw.clearance[], "los" => sw.los[],
                         "los_object" => String(sw.los_object[]), "capture" => sw.capture[],
                         "optimize" => sw.optimize[], "quality" => String(sw.quality[]))
    for k in (:los_theta, :los_phi)
        x = getfield(sw, k)[]
        isfinite(x) && (d[String(k)] = x)
    end
    return d
end

function _restore_sweep!(sw::SweepState, d::AbstractDict)
    for k in (:subject, :frame, :axis, :los_object, :quality)
        haskey(d, String(k)) && _assign!(getfield(sw, k), Symbol(d[String(k)]))
    end
    for k in (:start, :stop, :step)
        haskey(d, String(k)) && _assign!(getfield(sw, k), Float64(d[String(k)]))
    end
    for k in (:los_theta, :los_phi)
        _assign!(getfield(sw, k), Float64(get(d, String(k), NaN)))
    end
    for k in (:clearance, :los, :capture, :optimize)
        haskey(d, String(k)) && _assign!(getfield(sw, k), Bool(d[String(k)]))
    end
    return sw
end
