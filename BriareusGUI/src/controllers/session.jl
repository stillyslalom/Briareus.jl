# The Briareus window's state: the target, the beams, the analyses on them (drive,
# optimization, clearance, TaLIS) and what the window shows, the links between them (a
# target or port change makes the scenario stale; pointing edits re-evaluate it), and
# the background-job plumbing (`spawn`, `deliver`) shared by every job. Nothing here
# imposes an order: any analysis can be started at any time, and each reports whether
# it is current.

"""
Analyses of a design, in the model tree's order: `:drive` (view factors, laser
deposition and x-ray drive), `:optimize` (pointing optimization), `:clearance`
(worst-case beam clearance), `:sweep` (a placement sweep, [`SweepState`](@ref)) and
`:talis` (NIF target-configuration rules).
"""
const TASKS = (:drive, :optimize, :clearance, :sweep, :talis)

"""
    session_tasks(sess) -> Tuple

The analyses that apply to the session's facility: all of [`TASKS`](@ref) for NIF,
without the TaLIS rules (NIF's) elsewhere.
"""
session_tasks(sess) = sess.beams.facility[] === :nif ? TASKS : filter(!=(:talis), TASKS)

const TASK_LABELS = Dict(:drive => "Drive", :optimize => "Optimization", :clearance => "Clearance",
                         :sweep => "Placement sweep", :talis => "TaLIS")

"""
    EditHistory(; limit = 200)

Undo and redo stacks of design snapshots ([`session_state`](@ref) without the view:
coloring, selection and the camera are not edits). [`record!`](@ref) takes a snapshot
when the design has changed since the last one; [`undo!`](@ref) and [`redo!`](@ref)
step through them.

`changes` counts the changes of the session's observables that a snapshot reads (a
[`DesignSession`](@ref) watches them all), and `seen` is the count at which `current`
last matched the design: with no change since, the design is `current` without taking
a snapshot, so [`record!`](@ref) and [`can_undo`](@ref) cost nothing between edits.
"""
mutable struct EditHistory
    undo::Vector{Dict{String,Any}}
    redo::Vector{Dict{String,Any}}
    current::Union{Nothing,Dict{String,Any}}
    limit::Int
    changes::Threads.Atomic{Int}
    seen::Int
end
EditHistory(; limit::Integer = 200) =
    EditHistory(Dict{String,Any}[], Dict{String,Any}[], nothing, limit, Threads.Atomic{Int}(0), -1)

"""
    DesignSession(; starter = :hohlraum, resolution = :standard, preset = :nif_quads,
                  beams = nothing, deliver = f -> f())

State of the Briareus window: `target` ([`TargetSetup`](@ref)), `beams`
([`BeamTable`](@ref)), the analyses `drive` ([`DriveState`](@ref)), `optimize`
([`OptimizeState`](@ref)), `review` ([`ReviewState`](@ref), the clearance check) and
`talis` ([`TaLISState`](@ref)), `sweep` ([`SweepState`](@ref), the placement sweep),
`time` ([`TimeState`](@ref): pulses, the drive history
and the time cursor), `view` ([`ViewState`](@ref): coloring, beam drawing,
selection, model tree), `status` (the last message), `history` (the
[`EditHistory`](@ref) of [`undo!`](@ref) and [`redo!`](@ref)) and `meshing` (the
progress of a VisRad workspace import or of CAD parts put back into the calculations,
`""` when none runs; see [`import_visrad_workspace!`](@ref) and
[`include_components!`](@ref)). Designs are saved and
loaded with [`save_session`](@ref) and [`load_session!`](@ref); `cad_sources` maps each CAD
copy extracted from a saved session to the file it was saved from, and `cad_changed` lists
the copies whose source changed since the save.

The table starts from `beams` (a vector of `Beam`s) or a [`BEAM_PRESETS`](@ref) entry
pointed for the target ([`preset_pointing`](@ref)); `starter` is the parts the target
starts from ([`TargetSetup`](@ref)).

Background work: `spawn[]` (`false` by default; a window sets it) runs the scenario
build and its evaluations on worker tasks, handing their results to
`deliver[]` (a window queues them for its GUI thread). Without it they run inline.
The optimization, the clearance check and the TaLIS checks take an explicit `spawn`
keyword.
"""
struct DesignSession
    target::TargetSetup
    beams::BeamTable
    drive::DriveState
    optimize::OptimizeState
    review::ReviewState
    talis::TaLISState
    time::TimeState
    view::ViewState
    status::Observable{String}
    deliver::Base.RefValue{Any}
    spawn::Base.RefValue{Bool}
    history::EditHistory
    visrad_source::Observable{String}      # the imported VisRad workspace (export template)
    meshing::Observable{String}            # a workspace import or CAD meshing in progress ("" when none)
    mesh_generation::Base.RefValue{Int}    # import and meshing requests (latest wins)
    build_queued::Base.RefValue{Bool}      # a drive asked for while CAD parts remesh
    cad_sources::Dict{String,String}       # a saved CAD copy => the file it was saved from
    cad_changed::Observable{Vector{String}} # saved copies whose source changed since the save
    sweep::SweepState
end

function DesignSession(; starter::Symbol = :hohlraum, resolution::Symbol = :standard,
                       preset::Symbol = starter === :foil ? :omega_backlighter : :nif_quads,
                       beams = nothing, deliver = f -> f())
    ts = TargetSetup(; starter, resolution)
    bt = if beams === nothing
        p = preset_pointing(preset; target = ts)
        b = BeamTable(p.beams; preset, frames = part_frames(ts), aim_frame = p.aim_frame)
        b.enabled[] = p.enabled
        b
    else
        BeamTable(beams; frames = part_frames(ts))
    end
    # the target remeshes CAD parts on the session's workers (the same switches)
    ts.deliver[] = deliver
    sess = DesignSession(ts, bt, DriveState(), OptimizeState(), ReviewState(), TaLISState(),
                         TimeState(), ViewState(; cut_components = default_cut(ts)), Observable(""), ts.deliver,
                         ts.spawn, EditHistory(),
                         Observable(""), Observable(""), Ref(0), Ref(false), Dict{String,String}(),
                         Observable(String[]), SweepState())
    _default_selections!(sess, ts.target[])
    _connect!(sess)
    _watch_changes!(sess)
    record!(sess)
    return sess
end

# The observables an undo snapshot leaves out (`_history_state`): the window as it was
# left and the diagnostics a drive history records.
_unrecorded(sess::DesignSession) = (sess.view.camera, sess.drive.view, sess.drive.image_request,
                                    sess.time.chart, sess.time.diagnostics)

# Count every change of the session's observables but those (`EditHistory`'s `changes`):
# every observable field of the session and of its parts, so a field added later is
# counted too. A change that a snapshot does not read only costs `record!` a snapshot.
function _watch_changes!(sess::DesignSession)
    n = sess.history.changes
    skip = _unrecorded(sess)
    for x in (sess, sess.target, sess.beams, sess.drive, sess.optimize, sess.review, sess.talis,
              sess.time, sess.view, sess.sweep), f in fieldnames(typeof(x))
        o = getfield(x, f)
        o isa Observables.AbstractObservable && !any(s -> s === o, skip) || continue
        on(_ -> (Threads.atomic_add!(n, 1); nothing), o)
    end
    return sess
end

function Base.show(io::IO, s::DesignSession)
    print(io, "DesignSession(", target_summary(s.target), "; ", beams_summary(s.beams), ")")
end

# The checks' default roles for the components of `t`, from the parts'
# roles (walls, the observed surface and the foils are where the beams land).
function _default_selections!(sess::DesignSession, t::Target)
    l = _landing(sess.target, t)
    default_selection!(sess.review, t; l.lands_on)
    default_selection!(sess.talis, t; l.capture, l.observed)
    return sess
end

"""
    follow_renames!(sess::DesignSession, rn::AbstractDict) -> sess

Carry the session's state keyed by part, frame or component names through the renames
`rn` (old name → new, a [`rename_map`](@ref)): display modes, the cutaway's components,
drawn triads, the model tree's selection and expanded nodes, the drive's films and wall
cavity, the clearance and TaLIS roles, the sweep's subject, frame and line-of-sight
object, and the parts the ticked diagnostics' views name (`center_part`, `area_part`):
these follow in the session's `time.view_parts`, the definition (the user's file) left
as it is, and the time status notes it. The target's left-out components and the active
foil follow in the edit itself.
"""
function follow_renames!(sess::DesignSession, rn::AbstractDict)
    r(n::Symbol) = get(rn, n, n)
    rs(xs) = map(r, xs)
    function rename!(o::Observable, f)
        x = f(o[])
        x == o[] || (o[] = x)
    end
    v, dr, sw = sess.view, sess.drive, sess.sweep
    rename!(v.display, d -> Dict{Symbol,Symbol}(r(k) => m for (k, m) in d))
    rename!(v.cut_components, s -> Set{Symbol}(rs(collect(s))))
    rename!(v.triads, s -> Set{Symbol}(rs(collect(s))))
    ids = Dict{String,String}()
    for (a, b) in rn, kind in ("part:", "assembly:", "frame:", "comp:", "window:")
        ids[kind * String(a)] = kind * String(b)
    end
    _rename_nodes!(v, ids)
    rename!(dr.films, s -> Set{Symbol}(rs(collect(s))))
    dr.wall_cavity[] === nothing || rename!(dr.wall_cavity, r)
    for st in (sess.review, sess.talis), f in (st isa ReviewState ? (:targets, :hardware) :
                                              (:capture, :hardware, :diagnostics))
        rename!(getfield(st, f), rs)
    end
    sess.review.known[] = rs(sess.review.known[])
    sess.talis.known[] = rs(sess.talis.known[])
    rename!(sw.subject, r)
    rename!(sw.frame, r)
    rename!(sw.los_object, r)
    _follow_view_parts!(sess, rn)
    return sess
end

# The `[view]` keys of a diagnostic definition that name a part of the design.
const _VIEW_PART_KEYS = ("center_part", "area_part")

# The parts the views of the ticked diagnostics (and those already overridden) name
# follow the renames `rn`, in `time.view_parts`; the definitions stay as they are.
function _follow_view_parts!(sess::DesignSession, rn::AbstractDict)
    tm = sess.time
    vp = Dict{String,Dict{String,String}}(id => copy(d) for (id, d) in tm.view_parts[])
    moved = String[]
    for id in union(tm.diagnostics[], keys(vp))
        def = _definition(sess, id)
        def === nothing && continue
        for (k, name) in _view_parts(def, get(vp, id, nothing))
            new = get(rn, Symbol(name), nothing)
            new === nothing && continue
            get!(vp, id, Dict{String,String}())[k] = String(new)
            push!(moved, "$id's $k → $new")
        end
    end
    isempty(moved) && return sess
    _assign!(tm.view_parts, vp)
    tm.status[] = "the diagnostics follow the rename in this session, their files unchanged: " * join(moved, ", ")
    return sess
end

# The parts definition `def`'s view names, `key => name`, with the session's own `vp`
# (`time.view_parts[id]`, or `nothing`) in place of the file's.
function _view_parts(def, vp)
    v = get(def.data, "view", nothing)
    out = Pair{String,String}[]
    for k in _VIEW_PART_KEYS
        name = vp !== nothing && haskey(vp, k) ? vp[k] : v isa AbstractDict && haskey(v, k) ? String(v[k]) : nothing
        name === nothing || push!(out, k => name)
    end
    return out
end

function _connect!(sess::DesignSession)
    ts, bt, dr = sess.target, sess.beams, sess.drive
    # beams aimed in a frame move and turn with it (the frames change with the committed
    # inputs, before a target built on a worker is in place)
    on(fr -> follow_frames!(bt, fr; renames = ts.renames[]), ts.frames)
    # what is keyed by a renamed part's names follows the rename (before the target of the
    # new names prunes what it no longer has)
    on(ts.frames) do _
        rn = ts.renames[]
        isempty(rn) && return
        follow_renames!(sess, rn)
        # (and the role holders below: a renamed wall is not a new wall; they wait for the
        # target of the new names)
        r(s) = Set{Symbol}(get(rn, n, n) for n in s)
        foils[] = r(foils[])
        held[] = map(r, held[])
        observed[] === nothing || (observed[] = get(rn, observed[], observed[]))
        old = component_names(ts.target[])
        any(n -> haskey(rn, n), old) && (renamed[] = ts.target[])
    end
    # a cutaway left at its default cuts the hohlraum parts' shells as the parts change
    cut = Ref(default_cut(ts))
    on(ts.parts) do _
        now = default_cut(ts)
        now == cut[] && return
        sess.view.cut_components[] == cut[] && (sess.view.cut_components[] = now)
        cut[] = now
    end
    # an assembly removed or ungrouped is no longer selected
    on(ts.parts) do ps
        kind, key = selection_kind(sess.view.selection[])
        kind === :assembly && !any(p -> is_assembly(p) && String(p.name) == key, ps) &&
            select!(sess.view, "target")
    end
    on(ts.target) do t
        _default_selections!(sess, t)
        # forget display modes and a selection of components that no longer exist
        names = Set(all_component_names(ts))
        d = filter(p -> first(p) in names, sess.view.display[])
        d == sess.view.display[] || (sess.view.display[] = d)
        c = selected_component(sess)
        c === nothing || c in names || select!(sess.view, "target")
        _drive_inputs_changed!(sess)
    end
    # a beam selection whose beams went away (a preset, a split, a regrouping)
    onany(bt.beams, bt.cluster) do _...
        kind, _ = selection_kind(sess.view.selection[])
        kind in (:beam, :quad, :group) || return
        gone = filter(id -> isempty(node_beams(bt, id)), selected_ids(sess.view))
        isempty(gone) && return
        sess.view.selection[] in gone ? select!(sess.view, "beams") :
        (sess.view.extra[] = filter(id -> !(id in gone), sess.view.extra[]))
    end
    onany((_...) -> _drive_inputs_changed!(sess), bt.beams, bt.enabled, bt.unconverted)
    # a part given or relieved of a role (its components unchanged, so the target is not
    # rebuilt): the clearance and TaLIS checks take a wall or foil as a landing surface,
    # the observed surface as one the clearance check lands on and TaLIS leaves out, and a
    # part that lost its role as hardware; another observed surface makes the drive stale
    # (its scenario observes it), another wall the drive history (`_history_key`), and
    # another active foil re-evaluates the drive's foil statistics
    foils = Ref(Set(foil_components(ts)))
    held = Ref(_held_sets(ts))
    observed = Ref{Union{Nothing,Symbol}}(observed_name(ts))
    active = Ref{Union{Nothing,Symbol}}(active_foil(ts))
    renamed = Ref{Any}(nothing)                   # the target of a rename's old names
    onany(ts.parts, ts.target, ts.active_foil) do _...
        # (a renamed part is committed before the target of its new names: the roles
        # follow the rename itself, not a comparison with the target of the old names)
        if renamed[] !== nothing
            ts.target[] === renamed[] && return
            renamed[] = nothing
        end
        names = Set(component_names(ts.target[]))
        now = Set(foil_components(ts))
        if now != foils[]
            for n in setdiff(now, foils[])
                set_role!(sess.review, n, :target)
                set_talis_role!(sess.talis, n, :capture)
            end
            for n in setdiff(foils[], now)
                n in names || continue
                set_role!(sess.review, n, :hardware)
                set_talis_role!(sess.talis, n, :hardware)
            end
            foils[] = now
        end
        # (walls and the observed surface: a role the checks give a new component by
        # default, as after a rebuild, is left as it is)
        now = _held_sets(ts)
        if now != held[]
            was = held[]
            rv, tl = sess.review, sess.talis
            for n in setdiff(union(now...), union(was...))
                n in names && !(n in rv.targets[]) && set_role!(rv, n, :target)
            end
            for n in setdiff(now.walls, was.walls)
                n in names && !(n in tl.capture[]) && set_talis_role!(tl, n, :capture)
            end
            for n in setdiff(now.observed, was.observed)
                n in names && (n in tl.capture[] || n in tl.hardware[]) && set_talis_role!(tl, n, :none)
            end
            for n in setdiff(union(was...), union(now...))
                n in names || continue
                n in rv.hardware[] || set_role!(rv, n, :hardware)
                n in tl.hardware[] || set_talis_role!(tl, n, :hardware)
            end
            held[] = now
        end
        o = observed_name(ts)
        o === observed[] || (observed[] = o; _drive_inputs_changed!(sess))
        a = active_foil(ts)
        a === active[] || (active[] = a; _request_evaluation!(sess))
    end
    # a rejected edit's message belongs to what was being edited
    on(sess.view.selection) do _
        isempty(ts.error[]) || (ts.error[] = "")
        isempty(bt.error[]) || (bt.error[] = "")
    end
    on(_ -> _request_evaluation!(sess), dr.spot_radius)
    onany((_...) -> _drive_inputs_changed!(sess), dr.cavity, dr.cavity_point, dr.films,
          dr.beam_model, dr.profile)
    on(_ -> _drive_inputs_changed!(sess), dr.unconverted_ratio)
    # a drive asked for while CAD parts remeshed, the target was built or a capsule was
    # fitted runs when the target is ready
    onany(ts.pending, ts.building, ts.fitting) do _...
        target_ready(ts) && sess.build_queued[] && (sess.build_queued[] = false; build_scenario!(sess))
    end
    # inputs committed whose target is still being built: the drive is not current
    on(b -> b && _drive_inputs_changed!(sess), ts.building)
    # a drive history asked for before the drive was computed runs when it is
    on(dr.result) do res
        sess.time.queued[] && res !== nothing && _drive_current(sess) && start_history!(sess)
        sess.time.check_queued[] && res !== nothing && _drive_current(sess) && start_irradiance!(sess)
    end
    # the x-ray image follows the drive
    on(dr.result) do res
        dr.image_request[] === nothing && return
        dr.image[] = _render_image(dr)
        res !== nothing && res.solution === nothing &&
            (sess.status[] = "no x-ray image: the drive is laser only" * _radiation_note(sess))
    end
    return sess
end

"""
    set_aim_frame!(sess::DesignSession, idx, frame::Symbol; keep = true) -> sess

Hold the aims of beams `idx` in the frame named `frame` (`:tcc`, a part's or a hohlraum
part's piece frame, [`part_frames`](@ref)), so they move and turn with it. With `keep` the aims stay
where they are; with `keep = false` the beams aim at the frame's point
([`frame_point`](@ref): a part's centroid, a piece frame's origin).
"""
function set_aim_frame!(sess::DesignSession, idx, frame::Symbol; keep::Bool = true)
    haskey(part_frames(sess.target), frame) || throw(ArgumentError("no frame named $frame"))
    at = keep ? nothing : frame_point(sess, frame)
    set_aim_frame!(sess.beams, idx, frame; at)
    return sess
end

"""
    frame_point(sess::DesignSession, frame::Symbol) -> Vec3

The point a beam aims at in a frame when told to aim at the frame itself
([`set_aim_frame!`](@ref) with `keep = false`): a part's centroid (`part_centroid`, the
center of a disk or a foil; a CAD part's origin may lie anywhere), or the frame's origin
for TCC, a piece frame (a hohlraum part's entrance hole or window center) and a part
left out of the calculations.
"""
function frame_point(sess::DesignSession, frame::Symbol)
    F = part_frames(sess.target)[frame]
    i = part_index(sess.target, frame)
    i === nothing && return F.origin
    await_target!(sess.target)
    t = sess.target.target[]
    frame in Briareus.parts(t) || return F.origin
    return part_centroid(t, frame)
end

"""
    set_color!(sess::DesignSession, mode::Symbol)

Color the 3D view by `mode` (one of [`COLOR_MODES`](@ref)). A drive field needs the
drive: if it is missing or stale it is computed now (in the background in a window),
as VisRad computes what a display asks for.
"""
function set_color!(sess::DesignSession, mode::Symbol)
    mode in COLOR_MODES || throw(ArgumentError("color mode must be one of $COLOR_MODES"))
    _assign!(sess.view.color, mode)
    if mode !== :component
        set_field!(sess.drive, mode)
        _drive_current(sess) || sess.drive.building[] || build_scenario!(sess)
    end
    return sess
end

"""
    xray_image!(sess::DesignSession; direction, up = (0, 0, 1), center = (0, 0, 0),
                width, nx = 320) -> sess

Show an x-ray image of the target: the radiance (`B/π`, W/cm²/sr) of the surfaces seen
along `direction` (from the viewer into the scene), as an orthographic camera at
`center` imaging a `width` mm square (`Briareus.render`; apertures are transparent,
so a hohlraum is seen inside through its entrance holes). Like a coloring by a drive
field, it computes the drive if it is missing or stale, and it follows the drive: an
edit that re-evaluates the drive renders the image again.
"""
function xray_image!(sess::DesignSession; direction, up = (0.0, 0.0, 1.0), center = (0.0, 0.0, 0.0),
                     width::Real, nx::Integer = 320)
    width > 0 || throw(ArgumentError("the image width must be positive"))
    d = Vec3(direction...)
    norm(d) > 0 || throw(ArgumentError("the viewing direction must not be zero"))
    sess.drive.image_request[] = (; direction = d / norm(d), up = Vec3(up...),
                                  center = Vec3(center...), width = Float64(width), nx = Int(nx))
    _assign!(sess.drive.view, :image)
    if _drive_current(sess) && sess.drive.result[] !== nothing
        sess.drive.image[] = _render_image(sess.drive)
        sess.drive.result[].solution === nothing &&
            (sess.status[] = "no x-ray image: the drive is laser only" * _radiation_note(sess))
    else
        sess.drive.image[] = nothing
        sess.drive.building[] || build_scenario!(sess)
    end
    return sess
end

"""
    task_status(sess::DesignSession, task) -> (state, summary)

`state` is `:todo`, `:ok`, `:attention`, or `:busy`, with a summary short enough for
the model tree's default width (`task` one of [`TASKS`](@ref)); the analysis's
Properties and its results tab give it in full.
"""
function task_status(sess::DesignSession, task::Symbol)
    task in TASKS || throw(ArgumentError("unknown analysis :$task"))
    if task === :drive
        dr = sess.drive
        sess.build_queued[] && return (:busy, dr.status[])
        dr.building[] && return (:busy, "computing view factors… " * _seconds(time() - dr.started[]))
        dr.scenario[] === nothing && return (:todo, isempty(dr.status[]) ? "not computed" : dr.status[])
        drive_stale(sess) && return (:attention, "out of date: compute again")
        dr.evaluating[] && return (:busy, "evaluating…")
        dr.result[] === nothing && return (:attention, isempty(dr.status[]) ? "no evaluation" : dr.status[])
        return (:ok, drive_summary(dr.result[]; short = true))
    elseif task === :optimize
        st = sess.optimize
        st.running[] && return (:busy, optimize_progress(st))
        r = st.result[]
        r === nothing && return (:todo, isempty(st.status[]) ? "not run" : st.status[])
        txt = @sprintf("%.3g → %.3g", r.initial_value, r.value)
        optimize_stale(sess) && return (:attention, "out of date: " * txt)
        r.applied && return (:ok, txt * ", applied")
        return (:attention, txt * ", not applied")
    elseif task === :clearance
        rv = sess.review
        rv.running[] && return (:busy, "checking: $(rv.progress[][1]) of $(rv.progress[][2]) beams")
        rv.report[] === nothing && return (:todo, isempty(rv.status[]) ? "not run" : rv.status[])
        rows = rv.report[].rows
        txt = _counts_brief(rows) *
              (isempty(rows) ? "" : @sprintf(" · worst %.3f mm", minimum(x -> x.clearance, rows)))
        rv.report[].canceled && (txt = "partial: " * txt)
        review_stale(rv, sess.target.target[], active_beams(sess.beams)) &&
            return (:attention, "out of date: " * txt)
        bad = rv.report[].canceled || any(x -> x.status !== :pass, rv.report[].rows)
        return (bad ? :attention : :ok, txt)
    elseif task === :sweep
        sw = sess.sweep
        sw.running[] && return (:busy, "step $(sw.progress[][1]) of $(sw.progress[][2])")
        r = sw.result[]
        r === nothing && return (:todo, isempty(sw.status[]) ? "not run" : sw.status[])
        nf = count(w -> w.fails > 0 || w.los_clear === false, r.rows)
        txt = "$(length(r.rows)) steps" * (nf == 0 ? "" : ", $nf with a failure or blocked")
        r.canceled && (txt = "partial: " * txt)
        sweep_stale(sess) && return (:attention, "out of date: " * txt)
        return (nf > 0 || r.canceled ? :attention : :ok, txt)
    else
        tl = sess.talis
        tl.running[] && return (:busy, "checking… " * _seconds(time() - tl.started[]))
        tl.report[] === nothing && return (:todo, isempty(tl.status[]) ? "not run" : tl.status[])
        txt = _counts_brief(tl.report[].rows)
        talis_stale(tl, sess.review, sess.target.target[], active_beams(sess.beams)) &&
            return (:attention, "out of date: " * txt)
        bad = any(x -> x.status !== :pass, tl.report[].rows)
        return (bad ? :attention : :ok, txt)
    end
end

# A check's outcome in a few words, for the model tree: "8 fail, 48 review", "all pass".
function _counts_brief(rows)
    isempty(rows) && return "no checks"
    nf = count(x -> x.status === :fail, rows)
    nr = count(x -> x.status === :review, rows)
    parts = String[]
    nf > 0 && push!(parts, "$nf fail")
    nr > 0 && push!(parts, "$nr review")
    return isempty(parts) ? "all pass" : join(parts, ", ")
end

# An optimization result no longer describes the beam table: its beams left the table
# (a new preset or target) or moved to other ports, or, once applied, the
# optimized beams were edited since.
function _optimum_stale(bs::AbstractVector, r)
    all(i -> i <= length(bs), r.indices) || return true
    r.applied && return any(k -> bs[r.indices[k]] != r.beams[k], eachindex(r.indices))
    return any(k -> bs[r.indices[k]].port != r.beams[k].port, eachindex(r.indices))
end

# The design an optimization runs on: the drive's scenario key (target, ports, cavity,
# deposition profile, observed surface), the spot radius and the foil.
_optimize_key(sess::DesignSession) = (sess.drive.key[], sess.drive.spot_radius[], active_foil(sess.target))

"""
    optimize_stale(sess::DesignSession) -> Bool

Whether the optimization result is out of date: the drive it ran on is not the current
design's (another target, ports, cavity, deposition profile, observed surface, spot
radius or foil), or its beams left the beam table or moved to other ports, or, once
applied, the optimized beams were edited since. An undo or redo keeps the result: it is
current again when the restored design is the one it ran on.
"""
function optimize_stale(sess::DesignSession)
    r = sess.optimize.result[]
    r === nothing && return false
    r.key === nothing || (!drive_stale(sess) && r.key == _optimize_key(sess)) || return true
    return _optimum_stale(sess.beams.beams[], r)
end

# After an undo or redo: the optimization result is applied when the beam table holds its
# optimized beams.
function _restored_optimum!(sess::DesignSession)
    r = sess.optimize.result[]
    r === nothing && return sess
    bs = sess.beams.beams[]
    held = all(i -> i <= length(bs), r.indices) && all(k -> bs[r.indices[k]] == r.beams[k], eachindex(r.indices))
    held == r.applied || (r.applied = held; notify(sess.optimize.result))
    return sess
end

# ---------------------------------------------------------------- drive

"""
    drive_stale(sess::DesignSession) -> Bool

Whether the scenario was built for another target, other beam ports, cavity, beam
model or observed surface ([`observed_name`](@ref)) (it must be computed again), or the
target of the design is still being built. `false` before the first build.
"""
function drive_stale(sess::DesignSession)
    k = sess.drive.key[]
    k === nothing && return false
    return sess.target.building[] || k[1] !== sess.target.target[] || k[2] != _ports(sess.beams.beams[]) ||
           k[3] != _cavity_key(sess.drive) || k[4] != _profile_state(_deposition_profile(sess.drive)) ||
           k[5] !== observed_name(sess.target)
end

# The cavity settings a scenario was built with (`nothing`: the whole target).
_cavity_key(dr::DriveState) = dr.cavity[] ? (dr.cavity_point[], dr.films[]) : nothing

# Whether the scenario matches the current target and ports.
_drive_current(sess::DesignSession) = sess.drive.scenario[] !== nothing && !drive_stale(sess)

function _drive_inputs_changed!(sess::DesignSession)
    if _drive_current(sess)
        _request_evaluation!(sess)
    else
        sess.drive.eval_generation[] += 1                 # drop evaluations in flight
        sess.drive.evaluating[] = false
        sess.drive.result[] === nothing || (sess.drive.result[] = nothing)
    end
    return sess
end

"""
    build_scenario!(sess::DesignSession)

Build the core `Scenario` for the current target and every beam of the table
(deposition operator, view factors, capsule flux response; seconds to tens of
seconds), on a worker when `spawn[]` is set, then evaluate the current pointing. With
the drive's cavity on ([`set_cavity!`](@ref)) the scenario is built on the cavity's
elements, and its laser deposition uses the drive's beam model ([`set_beam_model!`](@ref)).
A target too large for the x-ray view factors ([`radiation_size`](@ref)) gets
the laser deposition alone, and `status` says why and what to do.
"""
function build_scenario!(sess::DesignSession)
    dr = sess.drive
    beams = sess.beams.beams[]
    isempty(beams) && (dr.status[] = "no beams: load a preset or import a beam table first"; return sess)
    ts = sess.target
    isempty(ts.parts[]) && (dr.status[] = "the target is empty: insert a part first"; return sess)
    # parts left remeshing, or a target left building, when a window closed: now
    sess.spawn[] || await_target!(ts)
    if !target_ready(ts)
        sess.build_queued[] = true
        dr.generation[] += 1                              # drop a build in flight
        dr.building[] = false
        pend = ts.pending[]
        dr.status[] = !isempty(pend) ?
                      "waiting for " * join(sort!(String.(collect(pend))), ", ") * " to be remeshed, then computing" :
                      "waiting for the target (" * target_waiting(ts) * "), then computing"
        return sess
    end
    t = ts.target[]
    sess.build_queued[] = false
    g = (dr.generation[] += 1)
    cavity = dr.cavity[] ? (dr.cavity_point[], copy(dr.films[])) : nothing
    profile = _deposition_profile(dr)
    dr.building[] = true
    dr.started[] = time()
    dr.status[] = cavity !== nothing ? "selecting the cavity, then its view factors…" :
                  radiation_size(t).feasible ? "computing view factors…" :
                  "computing the laser deposition (too many elements for the x-ray view factors)…"
    dr.status[] *= _facet_note(t)
    runner = _runner(sess.spawn, sess.deliver)
    t0 = time()
    cavities = _cavity_volumes(ts)
    mesh, parts = ts.mesh[], ts.parts[]
    foil = active_foil(ts)
    observed = observed_name(ts)
    runner(function ()
               # the CAD parts' file facets, for the drive's fit of their curved faces
               references = _fit_references(t, mesh, parts)
               _build_scenario(t, beams; cavities, cavity, profile, references, foil, observed)
           end, function (out)
        g == dr.generation[] || return
        dr.building[] = false
        if out.err !== nothing
            dr.status[] = "failed: " * _errmsg(out.err)
            return
        end
        (; s, lookup, enc, index, film, escape, leak) = out.value
        dr.build_seconds[] = time() - t0
        # a film found sealing the cavity joins the films (and the scenario's key)
        film === nothing || (_assign!(dr.films, union(dr.films[], (film,)));
                             cavity = (cavity[1], copy(dr.films[])))
        dr.key[] = (t, _ports(beams), cavity, _profile_state(profile), observed)
        dr.lookup[] = lookup
        dr.enclosure[] = enc
        dr.cavity_index[] = index
        dr.scenario[] = s
        dr.status[] = @sprintf("built in %s: %d elements, %d beams", _seconds(dr.build_seconds[]),
                               nelements(s.target), length(beams)) *
                      (enc === nothing ? "" : @sprintf(" (the cavity: %d of the target's %d elements)",
                                                       nelements(s.target), nelements(t))) *
                      (cavity !== nothing && enc === nothing ?
                       @sprintf(" · on the whole target: the cavity point is in the open (%.0f %% of its rays escape)",
                                100escape) : "") *
                      (film === nothing ? "" :
                       " · $film sealed the point in a pocket (a capsule tent or window modeled as a solid): " *
                       "left out of the cavity as a thin film (Properties)") *
                      (enc !== nothing && cavity_sealed(dr) ?
                       " · the cavity is a sealed pocket: mark the film that closes it (Properties ▸ Thin film)" : "") *
                      (enc === nothing || s.model === nothing ? "" : _closure_note(leak, out.value.closed)) *
                      (s.model === nothing ? "" : _cavities_note(out.value.cavities; enclosure = enc !== nothing)) *
                      (s.model === nothing ? _radiation_note(s.target) : _facet_note(s.target))
        _request_evaluation!(sess)
    end)
    return sess
end

# A cavity volume (`_CavityVolume`, drive.jl) of a hohlraum part.
function _CavityVolume(p::PartSpec, frame::Frame, ghost)
    prof = try
        hohlraum_profile(p)
    catch err
        err isa ArgumentError || rethrow()
        nothing                               # a profile the geometry refuses
    end
    pieces = Set(piece_name(p.name, s) for s in hohlraum_pieces(p))
    v = _CavityVolume(p.name, frame, _hohlraum_length(p), p.values[:radius], prof, pieces, false)
    ghost === nothing && return v
    # a piece left out, or another element left out that would bound the cavity (a
    # package face): a hole. One wholly inside it (a capsule) leaves it closed.
    holed = any(c -> c.name in pieces, ghost.components) ||
            any(e -> _in_cavity(v, ghost, e) && _on_boundary(v, ghost, e), 1:nelements(ghost))
    return _CavityVolume(v.name, v.frame, v.length, v.radius, v.profile, v.pieces, holed)
end

# The cavities of the design (`_CavityVolume`), one per hohlraum part in the calculations
# that has a frame, in the parts' order.
function _cavity_volumes(ts::TargetSetup)
    frames = part_frames(ts)
    g = ts.ghost[]
    return _CavityVolume[_CavityVolume(p, frames[p.name], g) for p in ts.parts[]
                         if p.kind === :hohlraum && !(p.name in ts.excluded[]) && haskey(frames, p.name)]
end

# Whether point `x` (chamber) lies in the cavity's closed volume, to `tol` (mm).
function _in_volume(v::_CavityVolume, x, tol)
    q = to_frame(v.frame, x)
    L = v.length
    abs(q[3]) <= L / 2 + tol || return false
    r = v.profile === nothing ? v.radius : Briareus._profile_radius(v.profile, clamp(q[3], -L / 2, L / 2))
    return hypot(q[1], q[2]) <= (isfinite(r) ? r : v.radius) + tol
end

# Whether element `e` of `t` belongs to the cavity by its place and side: its centroid in
# the closed volume (to 10⁻⁶ of the length) and its radiating side facing into it (10⁻³
# of the length in front of the centroid still inside). A package face in a halfraum's
# open end belongs to that halfraum; each side of a two-sided disk between two cavities
# to the cavity it faces.
function _in_cavity(v::_CavityVolume, t::Target, e::Integer)
    tol = 1e-6 * v.length
    x = t.elem_centroid[e]
    return _in_volume(v, x, tol) && _in_volume(v, x + 1e-3 * v.length * t.elem_normal[e], tol)
end

# Whether element `e` of the cavity lies on its boundary: one element size behind its
# centroid lies outside the volume (a wall, a cap, a package face), not inside (a capsule).
_on_boundary(v::_CavityVolume, t::Target, e::Integer) =
    !_in_volume(v, t.elem_centroid[e] - sqrt(t.elem_area[e]) * t.elem_normal[e], 1e-6 * v.length)

# The elements of `t` in each cavity of `vs`: the pieces of its hohlraum part by component
# name (entrance holes, lips and windows included), and every other element by
# `_in_cavity`, the first cavity that takes it holding it.
function _cavity_elements(t::Target, vs::AbstractVector{_CavityVolume})
    pieces = Dict(n => k for (k, v) in enumerate(vs) for n in v.pieces)
    out = [Int[] for _ in vs]
    for e in 1:nelements(t)
        k = get(pieces, t.components[t.elem_comp[e]].name, 0)
        if k == 0
            k = something(findfirst(v -> _in_cavity(v, t, e), vs), 0)
        end
        k > 0 && push!(out[k], e)
    end
    return out
end

# Each cavity of `vs` on the drive's target `t`: `(; name, elements, leak, holed,
# balanced)`, `leak` the fraction of its diffuse rays that leave it (core
# `open_fraction(t; cavity)`: they meet nothing, or an element outside it first), and
# `balanced` whether the drive balances it: no element of its boundary left out, and
# `leak` at most `CAVITY_MAX_LEAK` (an opening without an aperture, an open end that no
# part closes, or an element the selection missed shows as a leak). The rays check the
# geometric selection.
function _closed_cavities(t::Target, vs::AbstractVector{_CavityVolume})
    map(zip(vs, _cavity_elements(t, vs))) do (v, es)
        leak = isempty(es) ? 1.0 : open_fraction(t; cavity = es).fraction
        (; name = v.name, elements = es, leak, v.holed, balanced = !v.holed && !isempty(es) && leak <= CAVITY_MAX_LEAK)
    end
end
_closed_cavities(ts::TargetSetup) = _closed_cavities(ts.target[], _cavity_volumes(ts))

# The view factors' `closed` for the cavities balanced (`_closed_cavities`): `false` for
# none, `true` for one holding every element, else their element sets.
function _closure_kw(t::Target, cs)
    sets = [c.elements for c in cs if c.balanced]
    isempty(sets) && return false
    length(sets) == 1 && length(sets[1]) == nelements(t) && return true
    return sets
end

# The drive status's note on the cavities' balancing: nothing for one cavity balanced
# (a hohlraum and what lies inside it) of the whole target, else each cavity balanced,
# and why not. On a cavity point's `enclosure`, the cavities balanced in it are named.
function _cavities_note(cs; enclosure::Bool = false)
    isempty(cs) && return ""
    length(cs) == 1 && cs[1].balanced &&
        return enclosure ? " · view factors balanced on the cavity of $(cs[1].name)" : ""
    bal = [String(c.name) for c in cs if c.balanced]
    out = isempty(bal) ? "" : " · view factors balanced per cavity: " * join(bal, ", ")
    for c in cs
        c.balanced && continue
        out *= " · $(c.name) not balanced: " *
               (c.holed ? "a piece of its boundary is left out of the calculations" :
                @sprintf("%.2g %% of its rays leave it (an opening without an aperture, or an open end)", 100c.leak))
    end
    return out
end

"""
Elements the x-ray view factors may take (`0`: from the memory, see
[`radiation_size`](@ref)). Set it to change the limit for a session (tests do).
"""
const RADIATION_MAX_ELEMENTS = Ref(0)

"""
    radiation_size(t::Target) -> (; elements, bytes, limit, feasible)

Whether the drive can include the x-ray radiation of `t`: up to `Briareus.LARGE_MODEL[]`
elements (10⁴) the model keeps dense view factors with its radiosity factorizations,
three `N × N` Float64 matrices (24 N² bytes); up to `Briareus.HMATRIX_MODEL[]` (1.5 × 10⁴)
dense Float32 with factorization-free iterative solves (4 N² bytes, research notes,
"Large radiation models"); above it hierarchical view factors (`HMatrix`), budgeted at
250 N^1.5 bytes (2.5 × their storage on a hohlraum; research notes, "Hierarchical view
factors"). `bytes` must fit in half of this machine's memory; `limit` is the largest
element count that does (or [`RADIATION_MAX_ELEMENTS`](@ref), which then sets
`feasible`): about 160,000 elements on 32 GB. The mesh presets of the starters stay far
below (16,000 elements at `:fine`), but an imported CAD workspace can hold 10⁵ elements.
Above the limit [`build_scenario!`](@ref) computes the laser deposition alone.
"""
function radiation_size(t::Target)
    N = nelements(t)
    bytes = _radiation_bytes(N)
    RADIATION_MAX_ELEMENTS[] > 0 &&
        return (; elements = N, bytes, limit = RADIATION_MAX_ELEMENTS[], feasible = N <= RADIATION_MAX_ELEMENTS[])
    mem = 0.5 * Sys.total_memory()
    # the largest count that fits, on the cheapest path that reaches it
    large, hmat = Briareus.LARGE_MODEL[], Briareus.HMATRIX_MODEL[]
    nh = floor(Int, (mem / _HMATRIX_BYTES)^(2 / 3))
    limit = nh > hmat ? nh :
            (n32 = floor(Int, sqrt(mem / 4)); n32 > large ? min(n32, hmat) : min(large, floor(Int, sqrt(mem / 24))))
    return (; elements = N, bytes, limit, feasible = bytes <= mem)
end

# Bytes of a radiation model of N elements (see `radiation_size`): dense Float64 with two
# LU factors up to `LARGE_MODEL[]`, Float32 up to `HMATRIX_MODEL[]`, then a hierarchical
# matrix, budgeted as `_HMATRIX_BYTES` N^1.5 (its storage grows as N^1.5 on a hohlraum:
# 92 MB at 10.5k elements, 409 MB at 28.6k; the budget is 2.5 × that, for assembly and
# Krylov bases; research notes, "Hierarchical view factors").
_radiation_bytes(N::Integer) =
    N > Briareus.HMATRIX_MODEL[] ? _HMATRIX_BYTES * Float64(N)^1.5 :
    (N > Briareus.LARGE_MODEL[] ? 4.0 : 24.0) * Float64(N)^2
const _HMATRIX_BYTES = 250.0

_gb(x) = x >= 2^30 ? @sprintf("%.0f GB", x / 2^30) : @sprintf("%.0f MB", x / 2^20)

# Why a target's drive is laser only, and what to do ("" when it is not), naming the
# components that hold most elements.
function _radiation_note(t::Target)
    r = radiation_size(t)
    r.feasible && return ""
    ps = Briareus.parts(t)                  # the drive's fitted pieces of a part as the part
    counts = [length(elements(t, p)) for p in ps]
    big = sort(eachindex(counts); by = i -> -counts[i])[1:min(3, length(counts))]
    who = join(("$(ps[i]) ($(counts[i]))" for i in big if counts[i] > 0), ", ")
    return " · laser only: x-ray view factors for $(r.elements) elements would need " *
           "$(_gb(r.bytes)) (limit $(r.limit) elements); most are in $who. " *
           "Leave the mounting hardware out of the calculations (right-click ▸ Include in calculations) " *
           "or give the CAD parts a larger element size"
end
_radiation_note(sess::DesignSession) = _radiation_note(sess.target.target[])

"""
Flat elements (no fitted surface) carrying more facets than this on average make a
component's view factors slow: a curved CAD face cut into slivers. The drive's status
names such components ([`build_scenario!`](@ref)).
"""
const SLOW_FACETS_PER_ELEMENT = 50

# Why the view factors may be slow ("" when nothing stands out): flat components whose
# elements carry many facets each (≈ 20× the cost of the same tube fitted, research notes).
function _facet_note(t::Target)
    ne = zeros(Int, length(t.components))
    nt = zeros(Int, length(t.components))
    for e in 1:nelements(t)
        ne[t.elem_comp[e]] += 1
    end
    for k in eachindex(t.tri_comp)
        nt[t.tri_comp[k]] += 1
    end
    # (a part's flat piece left by the drive's fit is named by its part)
    slow = [(c.part, nt[i] / ne[i]) for (i, c) in enumerate(t.components)
            if c.surface === nothing && c.kind === :surface && ne[i] > 0 && nt[i] / ne[i] > SLOW_FACETS_PER_ELEMENT]
    isempty(slow) && return ""
    return " · slow view factors: " * join((@sprintf("%s (%.0f facets per element)", n, f) for (n, f) in slow), ", ") *
           ": a curved face tessellated in slivers on the flat path (the drive fits cylinders, cones and " *
           "spheres to CAD parts; this face is none of them, or the fit is off)"
end

"""
    unconverted_ratios(sess) -> Union{Nothing,NamedTuple}

The unconverted power of each beam of the table as a fraction of its 3ω power,
`(; omega1, omega2)` (vectors in table order): the drive's ratio
([`set_unconverted_ratio!`](@ref)) for the beams whose 1ω or 2ω light takes part
([`participates`](@ref)), zero for the others; `nothing` when no unconverted light takes
part (the drive is then the 3ω drive alone, bit for bit).
"""
function unconverted_ratios(sess)
    bt = sess.beams
    r1, r2 = sess.drive.unconverted_ratio[]
    o1 = [r1 * u[1] for u in bt.unconverted[]]
    o2 = [r2 * u[2] for u in bt.unconverted[]]
    any(>(0), o1) || any(>(0), o2) || return nothing
    return (; omega1 = o1, omega2 = o2)
end

# The scenario's flux-map lookup for the active foil `foil`: the one built with the
# scenario, made again when another foil became active since (a capsule or face map stays).
function _foil_lookup!(dr::DriveState, s::Scenario, foil)
    lk = dr.lookup[]
    lk !== nothing && lk.kind in (:capsule, :face) && return lk
    _lookup_foil(lk) === foil && return lk
    return dr.lookup[] = flux_lookup(s; foil)
end

# Evaluate the table's current pointing on the scenario (latest request wins).
function _request_evaluation!(sess::DesignSession)
    dr = sess.drive
    _drive_current(sess) || return sess
    s, radius = dr.scenario[], dr.spot_radius[]
    foil = active_foil(sess.target)
    lookup = _foil_lookup!(dr, s, foil)
    beams = display_beams(sess.beams)
    ratio = unconverted_ratios(sess)
    cavities = cavity_walls(sess.target)
    g = (dr.eval_generation[] += 1)
    dr.evaluating[] = true
    runner = _runner(sess.spawn, sess.deliver)
    runner(() -> _evaluate_drive(s, beams, lookup, radius; ratio, foil, cavities), function (out)
        g == dr.eval_generation[] || return
        dr.evaluating[] = false
        if out.err === nothing
            dr.result[] = out.value
        else
            dr.result[] = nothing
            dr.status[] = "evaluation failed: " * _errmsg(out.err)
        end
    end)
    return sess
end

"""
    drive_cavity(sess::DesignSession) -> Union{Symbol,Nothing}

The cavity (a hohlraum part, [`cavity_walls`](@ref)) whose wall temperature the Drive
panel shows: the one chosen with [`set_drive_cavity!`](@ref) while the latest drive has
it, else the one that holds the observed surface ([`observed_cavity`](@ref)), else the
first; `nothing` before a drive with a cavity. The choice selects what is shown only:
the drive and its history compute every cavity's temperature.
"""
function drive_cavity(sess::DesignSession)
    res = sess.drive.result[]
    (res === nothing || isempty(res.cavities)) && return nothing
    names = [r.name for r in res.cavities]
    c = sess.drive.wall_cavity[]
    c in names && return c
    o = observed_cavity(sess.target)
    return o in names ? o : first(names)
end

"""
    set_drive_cavity!(sess::DesignSession, name::Symbol) -> Bool

Show the wall temperature of cavity `name` (a hohlraum part in the calculations,
[`cavity_walls`](@ref)) in the Drive panel ([`drive_cavity`](@ref)). Another name sets
the drive's status and returns `false`. A display choice: not saved with the session and
not undone.
"""
function set_drive_cavity!(sess::DesignSession, name::Symbol)
    if !(name in first.(cavity_walls(sess.target)))
        sess.drive.status[] = "$name is not a hohlraum part in the calculations"
        return false
    end
    _assign!(sess.drive.wall_cavity, name)
    return true
end

"""
    include_components!(sess::DesignSession, names, include::Bool; spawn = sess.spawn[]) -> sess

Put the components `names` back into the calculations (`include = true`) or leave them
out (see [`set_excluded!`](@ref)). A name may be a part's (its components) or an
assembly's (the parts of its subtree: the model tree's include and exclude); a part put
back comes back whole. A CAD part left out is drawn from its file's facets;
putting it back repairs its windings and remeshes it, which takes tens of seconds for
10⁵ facets: with `spawn` that runs on worker threads, with progress in `meshing`, and
the target is replaced when it is ready. A rejected change sets the target's `error`.
"""
function include_components!(sess::DesignSession, names, include::Bool; spawn::Bool = sess.spawn[])
    ts = sess.target
    ex, names = _with_excluded(ts, ts.excluded[], names, !include)
    ex == ts.excluded[] && return sess
    if !spawn || !include
        _try_build!(ts; excluded = ex)
        return sess
    end
    mesh, parts = ts.mesh[], ts.parts[]
    g = (sess.mesh_generation[] += 1)
    deliver = sess.deliver[]
    label = join(string.(names), ", ")
    sess.meshing[] = "meshing $label for the calculations…"
    report(msg) = deliver(() -> g == sess.mesh_generation[] && _assign!(sess.meshing, msg))
    job = function ()
        back = [p for p in parts if p.name in names]
        _mesh_parts(back, mesh, ex, report, label; all = parts)
        return _split_components(_target_components(mesh, parts; raw = ex), ex)
    end
    commit = function (out)
        g == sess.mesh_generation[] || return
        _assign!(sess.meshing, "")
        if out.err !== nothing
            ts.error[] = _errmsg(out.err)
            sess.status[] = "could not put $label back: " * ts.error[]
            return
        end
        same = ts.mesh[] == mesh && ts.parts[] === parts
        _try_build!(ts; excluded = ex, prebuilt = same ? out.value : nothing) &&
            (sess.status[] = "$label back in the calculations" * _radiation_note(sess))
    end
    _run_job(() -> _deliver_compiled(deliver, commit, _try_job(job)), true)
    return sess
end

# Names of the CAD (mesh) objects reaching farther from TCC than twice the farthest aim
# point of the beams that fire (none without such beams): mounting hardware.
function _far_objects(objs, beams)
    aims = [norm(b.aim) for b in beams if b.power > 0]
    isempty(aims) && return String[]
    reach = 2 * max(maximum(aims), 1.0)
    return [o.name for o in objs if o.kind === :mesh && o.mesh !== nothing &&
                                    maximum(norm, o.mesh.vertices) > reach]
end

"""
    import_visrad_workspace!(sess::DesignSession, path; include = o -> o.radiosity || o.laser,
                             display_only = :far, spawn = sess.spawn[], done = () -> nothing)
        -> sess

Replace the design with a VisRad workspace's (`read_visrad_workspace`): a custom target
whose parts are the workspace's objects (`:visrad` parts: each keeps its place, albedo
and conversion; CAD windings are repaired and surfaces remeshed to the target's element
size), in one assembly named after the file (its frame at TCC, the objects placed in it,
so moving the assembly moves them all: [`add_assembly!`](@ref)), and the workspace's
beams, in their clusters as groups, switched on where VisRad delivers power. Objects for which `include` is false (by default those VisRad
leaves out of both its laser and radiosity calculations) are imported but left out of
the calculations, drawn as ghosts. The workspace becomes the template of
[`export_visrad_workspace`](@ref), and the drive deposits with VisRad's beam model and
the workspace's parameters (`visrad_profile`; [`set_beam_model!`](@ref) switches back).

CAD workspaces bring the mounting hardware along: parts of 10⁵ facets reaching tens of
mm from the target, which would take most of the elements (and the x-ray view factors
out of reach, see [`radiation_size`](@ref)) and tens of seconds of remeshing. With
`display_only = :far` (the default) a CAD part that reaches farther from TCC than twice
the farthest aim point of the beams VisRad fires is left out of the calculations too,
drawn from its file's facets; `status` names these parts, and
[`include_components!`](@ref) or *Include in calculations* puts one back (laser,
clearance and the x-ray drive then see it). `display_only = :none` follows VisRad's
flags alone.

Reading a CAD workspace and remeshing its parts takes seconds, so with `spawn` (a
window sets it) the file is read and the parts meshed on worker threads, with progress
in `meshing`, and the design is replaced on the GUI thread when it is ready (then
`done()` runs); a later import supersedes one in progress. Errors are reported in
`status` then. Inline (`spawn = false`) the import completes before returning and
throws `ArgumentError` for a workspace it cannot use.
"""
Base.@nospecializeinfer function import_visrad_workspace!(
        sess::DesignSession, path::AbstractString;
        @nospecialize(include = o -> o.radiosity || o.laser), display_only::Symbol = :far,
        spawn::Bool = sess.spawn[], @nospecialize(done = () -> nothing))
    # one compiled copy serves every `include` and `done` (a window's, a script's): the
    # body is long, and specialized on each closure it took 0.8 s to compile per caller
    display_only in (:far, :none) || throw(ArgumentError("display_only must be :far or :none"))
    path = abspath(String(path))
    ts = sess.target
    mesh = ts.mesh[]
    g = (sess.mesh_generation[] += 1)
    deliver = sess.deliver[]
    file = basename(path)
    sess.meshing[] = "reading $file…"
    report(msg) = spawn && deliver(() -> g == sess.mesh_generation[] && _assign!(sess.meshing, msg))
    job = function ()
        ws = _visrad_workspace(path)
        supported = [o for o in ws.objects if o.kind in (:mesh, :sphere, :rectangle, :cylindrical_hohlraum)]
        isempty(supported) && throw(ArgumentError("$file has no target components Briareus reads"))
        any(include, supported) ||
            throw(ArgumentError("VisRad includes none of the workspace's objects in its calculations"))
        # one assembly of the objects, named after the file, its frame at TCC (an object
        # named capsule is the observed surface, as VisRad workspaces name it)
        asm = _free_identifier(first(splitext(file)), Set{Symbol}(Symbol(o.name) for o in supported), "workspace")
        parts = [[PartSpec(:visrad; name = Symbol(o.name), file = path, object = o.name,
                           values = Dict(:albedo => o.albedo, :conversion => o.conversion),
                           observed = o.name == "capsule", frame = asm, assembly = asm)
                  for o in supported]; PartSpec(:assembly; name = asm)]
        # the components of objects VisRad leaves out are left out of the calculations,
        # and so (by default) is CAD hardware far from where the beams point
        names = Set{Symbol}(n for o in supported if !include(o) for n in visrad_component_names(o))
        far = display_only === :far ? _far_objects(filter(include, supported), ws.beams) : String[]
        length(far) < count(include, supported) || (far = String[])
        union!(names, Symbol.(far))
        _mesh_parts(parts, mesh, names, report, file)
        report("$file: building the target…")
        built = _split_components(_target_components(mesh, parts; raw = names),
                                  names)
        return (; ws, parts, names, far, built, profile = visrad_profile(ws))
    end
    commit = function (out)
        g == sess.mesh_generation[] || return
        _assign!(sess.meshing, "")
        if out.err !== nothing
            sess.status[] = "import of $file failed: " * _errmsg(out.err)
            spawn || throw(out.err)
            return
        end
        (; ws, parts, names, far, built, profile) = out.value
        # a target edit while the worker ran: build again from the current settings
        same = ts.mesh[] == mesh
        if !_try_build!(ts; parts, excluded = names,
                        prebuilt = same ? built : nothing)
            sess.status[] = "the workspace's target was rejected: $(ts.error[])"
            spawn || throw(ArgumentError(sess.status[]))
            return
        end
        # the workspace's assembly is shown open in the model tree
        _assign!(sess.view.expanded, union(sess.view.expanded[], ("assembly:" * String(parts[end].name),)))
        bs = ws.beams
        set_beams!(sess.beams, bs; preset = :custom, enabled = [b.power > 0 for b in bs],
                   cluster = Symbol[b.group for b in bs])
        sess.beams.path[] = path
        sess.visrad_source[] = path
        # a CAD assembly's drive runs on its cavity around TCC (outer faces and mounting
        # hardware out of the view factors)
        set_cavity!(sess.drive, any(o -> o.kind === :mesh, ws.objects); point = (0.0, 0.0, 0.0),
                    films = Symbol[])
        # deposit as VisRad does, with the workspace's envelopes and focal length
        set_beam_model!(sess.drive, :visrad; profile)
        # VisRad workspaces list every port: draw the beams that fire when most do not
        _assign!(sess.view.show_off, 2 * count(sess.beams.enabled[]) >= length(bs))
        # the cutaway cuts every component of an imported design (its names are VisRad's)
        obs = observed_name(ts)
        set_cut!(sess.view; components = [n for n in surface_names(ts.target[]) if n !== obs])
        sess.status[] = "imported $file: $(length(parts) - 1) objects, $(length(bs)) beams" *
                        (isempty(far) ? "" :
                         "; drawn only, far from the target (Include in calculations puts one back " *
                         "for the laser, clearance and drive): " * join(far, ", ")) *
                        _radiation_note(sess)
        done()
    end
    if spawn
        _run_job(() -> _deliver_compiled(deliver, commit, _try_job(job)), true)
    else
        commit(_try_job(job))
    end
    return sess
end

"""
    export_visrad_workspace(sess::DesignSession, path) -> path

Write the design's beams (pointing, focus, power, on/off) into a copy of the VisRad
workspace the design was imported from ([`import_visrad_workspace!`](@ref)), with
`write_visrad_workspace`: VisRad's target components, frames, spots and settings stay
as they were. Disabled beams are written off.
"""
function export_visrad_workspace(sess::DesignSession, path::AbstractString)
    src = sess.visrad_source[]
    isempty(src) && throw(ArgumentError("import a VisRad workspace first: the export writes into a copy of it"))
    bs = [on ? b : remake(b; power = 0.0) for (b, on) in zip(sess.beams.beams[], sess.beams.enabled[])]
    write_visrad_workspace(path, src, bs)
    sess.status[] = "wrote the beams into a copy of $(basename(src)): $(basename(path))"
    return path
end

"""
    export_target_stl(sess::DesignSession, path; binary = true, apertures = true) -> path

Write the design's target as an STL file in millimeters (`Briareus.write_stl`): the
components in the calculations (those left out are not written), each triangle wound
toward its radiating side. Binary STL by default; `binary = false` writes ASCII STL with
one named solid per component, which `read_stl_solids` (or a CAD tool) reads back as
parts. Apertures (entrance holes, open windows) are written as their own solids unless
`apertures = false`.
"""
function export_target_stl(sess::DesignSession, path::AbstractString; binary::Bool = true,
                           apertures::Bool = true)
    await_target!(sess.target)
    t = sess.target.target[]
    write_stl(path, t; binary, apertures)
    sess.status[] = "wrote the target ($(length(t.components)) components, $(length(t.tri_area)) " *
                    "triangles, mm) to $(basename(path))"
    return path
end

"""
    export_results(sess::DesignSession, dir) -> Vector{String}

Write the current drive's results to the folder `dir` for postprocessing
(`Briareus.write_results`: per-element laser and x-ray fields, per-component and per-beam
power, the capsule flux's Legendre and Y_lm modes, a summary) and the design that
produced them (`session.briar.toml`, see [`save_session`](@ref)). With a current drive
history its results follow ([`export_history`](@ref): the drive, the energy and each
recorded diagnostic against time, the fields per element and step, the imagers'
images), and the image in the image view when there is one. The drive must be computed
and current; disabled beams appear at zero power. Returns the files written.
"""
function export_results(sess::DesignSession, dir::AbstractString)
    dr = sess.drive
    _drive_current(sess) && dr.result[] !== nothing ||
        throw(ArgumentError("compute the drive first: the export writes its results"))
    res = dr.result[]
    ev = res.evaluation
    t = ev.scenario.target
    observe = _observed_parts(ev.scenario)
    files = write_results(dir, t, ev.deposition; sol = res.solution, beams = ev.effective,
                          observe = res.solution === nothing || isempty(observe) ? nothing : observe)
    h = sess.time.history[]
    with_history = h !== nothing && !history_stale(sess)
    with_history && append!(files, export_history(sess, dir))
    img = dr.image[]
    img === nothing || push!(files, _write_image(joinpath(dir, "image.csv"), img.x, img.y, img.values,
                                                 img.title * ", " * img.unit))
    push!(files, save_session(joinpath(dir, "session.briar.toml"), sess))
    sess.status[] = "exported the results to $(basename(normpath(dir)))" *
                    (with_history ? ", with the drive history" :
                     h === nothing ? "" : " (the drive history is out of date: not exported)")
    return files
end

# ---------------------------------------------------------------- retune

# The key a pulse table's columns are: the beams' groups, quads or names.
function _pulse_key(beams, pulses)
    for by in (:group, :quad, :beam)
        key = Briareus._beam_key(by)
        all(b -> haskey(pulses, key(b)), beams) && return by
    end
    throw(ArgumentError("the pulse table's columns are not the beams' groups, quads or names"))
end

"""
    retune_pulses(s::Scenario, beams, enabled, pulses; by = :quad, cap = 1.3, stay = 0.0,
                  target = :reference, level = nothing, wall = (), unconverted = nothing,
                  prepulse = nothing, foil = nothing, spot_radius = 0.5, on = :auto) -> NamedTuple

Retune the enabled beams' powers to the drive of a reference shot in which every beam of
`beams` (the scenario's, in its order) fired `pulses` (power per beam against time, keyed
by group, quad or beam name, e.g. `read_pulse_table`): at each step of a drive history
of the reference ([`start_history!`](@ref)'s steps: the pulses' times, the prepulse
window every 0.5 ns and one step after it, and the zero-power steps between, from the
first time anything fires to the last), where the reference's 3ω light fires (elsewhere
the retuned powers are zero), the drive is held (`Briareus.retune_drive`) by one power
per `by` key (`:quad`, `:group` or `:beam`) within `cap` times its reference power,
preferring the reference powers with weight `stay`. `target = :reference` holds the
reference's drive shape, `:uniform` a uniform drive (no modes; a flat flux), the least
nonuniformity the powers allow; `level` sets the drive's level at every step (TW/mm²:
the face's or the foil spot's mean, or the capsule's mean incident flux, whose `c_00` is
`√(4π)` times it), and `nothing` holds the reference's level at each step. Where the
power limits cannot reach the level the retune reports the level they allow
(`retune.drive`). Disabled beams are the dropped ones. A cone's or a quad's column of
`pulses` is the power per `unit` ([`PULSE_UNITS`](@ref): per quad by default, shared
among the quad's beams), and so are the new pulses' columns (the result's `unit`); a key
whose quads are held some as one beam and some as four is refused (its beams would not
share one power). What is held (`kind`):

- `:capsule`, an observed sphere: its modes through l = 4 (`Briareus.CapsuleModes`);
- `:face`, an observed flat face ([`observed_flat`](@ref)): its incident x-ray flux on
  each element, so its RMS nonuniformity, peak-to-mean and mean (`Briareus.FaceFlux`);
- `:foil`, the backlighter foil `foil` (a part or component of the scenario's target,
  e.g. [`active_foil`](@ref)): its absorbed laser irradiance (3ω, and the unconverted
  light below) within `spot_radius` (mm) of the foil's area centroid on its front, as the
  drive's spot statistics (`Briareus.FoilIrradiance`); laser only, no radiation model
  needed.

`on = :auto` holds the observed surface when the scenario has a radiation model and an
observed surface, else the foil; `:observed` or `:foil` chooses. Returns the new `pulses`
(by `by`, on the table's times and the steps, zero outside the steps and for the keys of
the dropped beams, so that every beam of `beams` has one), the `kind`, the `retune`, and the `reference`,
the enabled beams on the reference pulses (`before`) and on the new ones (`after`): for
a capsule, drive histories, whose wall temperature is that of the names `wall` in the
scenario's target (a session passes its [`wall_names`](@ref)), else the capsule's; for
a face or a foil, `(; time, mean, peak, low, rms)` per time (W/cm², the RMS relative), as
the drive's face and spot statistics, the face's from drive histories; and the `target`
and `level` asked for. With
`unconverted` (per color, the unconverted power as a fraction of each beam's 3ω power,
`(; omega1, omega2)`: a number for every beam, or a vector with one per beam of
`beams`), every history and the retune take the unconverted 1ω and 2ω light
(`drive_history`'s and `retune_drive`'s `unconverted`, NIF's colors; the enabled beams'
fractions in the retune); for a foil, each color's absorbed power adds to the spot's
irradiance in the reference, the retune and the statistics, as `drive_history` deposits
it. `prepulse` adds light not proportional to the 3ω pulse, `(; omega1, omega2,
window)`: per color a power per beam of `beams` (TW, a vector), flat over `window =
(start, stop)` (ns); the histories and the retune take it as `drive_history`'s and
`retune_drive`'s `prepulse`, keyed by `by` (the beams of one key must carry the same
prepulse), and the retune does not retune it. Without either they are the 3ω ones, bit
for bit.
"""
function retune_pulses(s::Scenario, beams::AbstractVector{<:Beam}, enabled::AbstractVector{Bool},
                       pulses::AbstractDict; by::Symbol = :quad, cap::Real = 1.3, stay::Real = 0.0,
                       unit::Symbol = :quad,
                       target::Symbol = :reference, level::Union{Nothing,Real} = nothing,
                       wall = (), unconverted::Union{Nothing,NamedTuple} = nothing,
                       prepulse::Union{Nothing,NamedTuple} = nothing, foil::Union{Nothing,Symbol} = nothing,
                       spot_radius::Real = 0.5, on::Symbol = :auto)
    kind = _retune_kind(s, foil, on)
    kind === nothing && throw(ArgumentError(_RETUNE_NONE))
    by in (:quad, :group, :beam) || throw(ArgumentError("retune by :quad, :group or :beam, not :$by"))
    target in (:reference, :uniform) || throw(ArgumentError("the retune's target is :reference or :uniform, not :$target"))
    level === nothing || (isfinite(level) && level > 0) || throw(ArgumentError("the level must be positive"))
    cap > 0 || throw(ArgumentError("the power limit must be positive"))
    any(enabled) || throw(ArgumentError("every beam is disabled"))
    unit in PULSE_UNITS || throw(ArgumentError("a pulse table's unit is one of $PULSE_UNITS"))
    ref_by = _pulse_key(beams, pulses)
    refkey, key = Briareus._beam_key(ref_by), Briareus._beam_key(by)
    # each beam's share of its column, in the reference table and in the retuned one
    # (`_pulse_shares`); a retuned key gives one power to its beams, so they must share
    # alike: a cone whose quads are held some as one beam and some as four cannot
    rsh, osh = _pulse_shares(beams, ref_by, unit), _pulse_shares(beams, by, unit)
    keyshare = Dict{Any,Float64}()
    for (b, sh) in zip(beams, osh)
        get!(keyshare, key(b), sh) == sh ||
            throw(ArgumentError("$(key(b)) holds quads as one beam and as four: split or merge them all, " *
                                "or retune with one power per beam"))
    end
    nodes = sort(unique(reduce(vcat, [p.t for p in values(pulses)])))
    # the steps of a drive history of the reference (`_history_times`): the table's times,
    # the prepulse window and one step after it, gaps where nothing fires; a step without
    # 3ω light in the reference has none in the retune (its limits are zero)
    times = _history_times((; kind = :table, table = pulses, by = ref_by, times = nodes, unit), prepulse, beams,
                           fill(true, length(beams)))
    any(τ -> sum(b -> pulses[refkey(b)](τ), beams) > 0, times) ||
        throw(ArgumentError("the reference pulses have no power"))
    op = s.deposition
    keep = findall(enabled)
    kept = beams[keep]
    opk = op[keep]
    # a key's reference power: its beams' mean (the table may be finer than the retune)
    members = Dict(k => [(refkey(beams[i]), rsh[i]) for i in keep if key(beams[i]) == k] for k in unique(key.(kept)))
    nominal(k, τ) = sum(sh * pulses[r](τ) for (r, sh) in members[k]) / length(members[k])
    limits(k, τ) = (0.0, cap * nominal(k, τ))
    # the level in the quantity's units: a capsule's c_00 = √(4π) × its mean flux
    lv = level === nothing ? nothing : kind === :capsule ? sqrt(4π) * Float64(level) : Float64(level)
    lkw = (; by, limits, stay, nominal = stay > 0 ? nominal : nothing, target, level = lv)
    refbeams(τ) = _table_beams(beams, pulses, τ; by = ref_by, shares = rsh)
    oldbeams(τ) = _table_beams(kept, pulses, τ; by = ref_by, shares = rsh[keep])
    # the unconverted light of every beam (the reference's) and of the kept beams, as
    # `drive_history`'s keywords, and the latter as the retune's (its `by` is its own)
    uc = _light_kw(unconverted, prepulse, beams; by)
    uck = _light_kw(_kept_unconverted(unconverted, keep), _kept_unconverted(prepulse, keep), kept; by)
    ucr = Base.structdiff(uck, NamedTuple{(:by,)})
    if kind === :foil
        t = s.target
        es = _spot_elements(t, foil, spot_radius)
        es === nothing && throw(ArgumentError("the foil $foil has no front elements within $spot_radius mm of its center"))
        ref = (; time = times, deposition = _absorbed_light(op, refbeams, times, uc, by))
        r = retune_drive(FoilIrradiance(foil; radius = spot_radius), opk, kept, ref; lkw..., ucr...)
        newbeams(τ) = with_pulses(kept, r.pulses, τ; by)
        series(o, at, kw) = _spot_series(t, es, times, [d.element for d in _absorbed_light(o, at, times, kw, by)])
        reference, before, after = series(op, refbeams, uc), series(opk, oldbeams, uck), series(opk, newbeams, uck)
    else
        observe = _observed_parts(s)
        # the histories' wall: the walls in the scenario's target (`wall_names`), else the
        # observed surface
        names = component_names(s.target)
        walls = Tuple(n for n in wall if n in names)
        wall = isempty(walls) ? observe : walls
        face = kind === :face
        # a face's statistics per step, as the drive's (`_face_stats`)
        hkw = face ? (; record = sol -> _face_stats(s, sol)) : (; ylm = true)
        hist(o, at, kw) = drive_history(s.model, o, at, times; observe, wall, lmax = 4, hkw..., kw...)
        ref = hist(op, refbeams, (; keep_albedo = true, observed_flux = face, uc...))
        q = face ? FaceFlux(s.model, observe) : CapsuleModes(s.model; observe)
        r = retune_drive(q, opk, kept, ref; lkw..., ucr...)
        newbeams2(τ) = with_pulses(kept, r.pulses, τ; by)
        reference, before, after = ref, hist(opk, oldbeams, uck), hist(opk, newbeams2, uck)
        if face
            reference, before, after = _face_series(reference), _face_series(before), _face_series(after)
        end
    end
    # the new pulses on the table's times and the steps, off outside the steps, and off
    # for the keys of the dropped beams (so that the table keys every beam); in the
    # table's unit: a key's power per beam over its beams' share of the column
    at(g, τ) = (k = findfirst(==(τ), times); k === nothing ? 0.0 : r.power[g, k])
    ts = sort(unique([nodes; times]))
    out = Dict{Symbol,PulseShape}(Symbol(k) => PulseShape(ts, [at(g, τ) / keyshare[k] for τ in ts])
                                  for (g, k) in enumerate(r.keys))
    for k in unique(key.(beams))
        haskey(out, k) || (out[k] = PulseShape(ts, zeros(length(ts))))
    end
    return (; pulses = out, unit = by === :beam ? :beam : unit, kind, target, level, retune = r, reference, before,
            after)
end

# What a retune of scenario `s` holds (`retune_pulses`'s `kind`; `on`: `:auto`,
# `:observed` or `:foil`), or `nothing` when it has neither an observed surface under a
# radiation model nor the foil `foil`.
function _retune_kind(s::Scenario, foil, on::Symbol)
    on in (:auto, :observed, :foil) || throw(ArgumentError("retune on :auto, :observed or :foil, not :$on"))
    observed = s.model !== nothing && !isempty(s.observed) ? (_observed_flat(s) ? :face : :capsule) : nothing
    hasfoil = _has_part(s.target, foil) ? :foil : nothing
    on === :observed && return observed
    on === :foil && return hasfoil
    return something(observed, hasfoil, Some(nothing))
end

# What a retune of the drive's scenario can hold: `(; on, kind, label)` per choice (`on`,
# `retune_pulses`'s: `:observed` or `:foil`), the one `on = :auto` takes first; none
# without a scenario.
function _retune_quantities(sess::DesignSession)
    s = sess.drive.scenario[]
    s === nothing && return NamedTuple[]
    foil = active_foil(sess.target)
    out = NamedTuple[]
    obs = _retune_kind(s, foil, :observed)
    name = String(something(observed_name(sess.target), :observed))
    obs === nothing ||
        push!(out, (; on = :observed, kind = obs, label = obs === :capsule ? "$name's modes" : "$name's flux"))
    _retune_kind(s, foil, :foil) === nothing || push!(out, (; on = :foil, kind = :foil, label = "$foil's irradiance"))
    return out
end

# A retune's `level` from the dialog's text in the units of the quantity `kind`: a
# capsule's Tᵣ (eV) as its mean flux σTᵣ⁴, a face's or a foil spot's mean (TW/mm²); an
# empty text is the reference's level (`nothing`).
function _retune_level(kind::Symbol, text::AbstractString)
    isempty(strip(text)) && return nothing
    x = _parse_number(text, kind === :capsule ? "the level's Tᵣ" : "the level"; positive = true)
    return kind === :capsule ? Briareus.SIGMA_SB * x^4 : x
end

# Why a drive cannot be retuned: nothing to hold.
const _RETUNE_NONE = "the retune holds an observed surface's x-ray drive or a backlighter foil's laser " *
                     "irradiance: observe a part (Properties ▸ Role), or make one a backlighter foil"

# A face's statistics per step of a drive history recorded with `_face_stats`.
_finite_max(x) = maximum(filter(isfinite, x); init = 0.0)

_face_series(h) = (; time = h.time, mean = [f.mean for f in h.record], peak = [f.peak for f in h.record],
                   low = [f.low for f in h.record], rms = [f.rms for f in h.record])

# A foil spot's statistics (W/cm²; the RMS relative) at each of `times` from each time's
# absorbed laser power per element `P`, over the spot elements `es`, as `_spot_stats`.
function _spot_series(t::Target, es, times, Ps)
    A = t.elem_area[es]
    I = [P[es] ./ A for P in Ps]
    mean = [sum(i .* A) / sum(A) for i in I]
    return (; time = collect(times), mean = to_Wcm2.(mean), peak = to_Wcm2.(maximum.(I)),
            low = to_Wcm2.(minimum.(I)), rms = [m > 0 ? rms_nonuniformity(i, A) : NaN for (i, m) in zip(I, mean)])
end

# The absorbed laser power on the target of `op` from the beams `at(τ)` (3ω setup powers)
# at each of `times`: the 3ω light's `Deposition`, and with the unconverted light `kw`
# (`_light_kw`'s; empty: none, the 3ω deposition itself) each color's absorbed power
# added per element, `(; element)`, as `drive_history` deposits it (one
# `_UnconvertedLight` per call: its unit-power patterns serve every time).
function _absorbed_light(op, at, times, kw, by)
    isempty(kw) && return [deposit(op, at(τ)) for τ in times]
    ul = Briareus._UnconvertedLight(kw.unconverted, get(kw, :prepulse, nothing), get(kw, :by, by), kw.colors, (;))
    return map(times) do τ
        bs = at(τ)
        E = copy(deposit(op, bs).element)
        u = Briareus._unconverted_at!(ul, op.target, bs, τ)
        for c in ul.colors
            E .+= u[c.name].element
        end
        (; element = E)
    end
end

# A retune's face or foil statistics `f` (`retune_pulses`'s `before`, `after`): the mean's
# relative departure at each step from the reference's mean `ref.mean`, or with a
# `level` (TW/mm²) from it where the retune solved (`solved`; NaN elsewhere).
_mean_departure(f, level, solved, ref) =
    level === nothing ? f.mean ./ ref.mean .- 1 : [s ? m / to_Wcm2(level) - 1 : NaN for (m, s) in zip(f.mean, solved)]

# The retune's level as text: the reference's, or the set one.
_level_text(v) = v.level === nothing ? "the reference's level" :
                 v.kind === :capsule ? @sprintf("Tᵣ %.1f eV", radiation_temperature(v.level)) :
                 @sprintf("%.3g TW/mm²", v.level)

# The unconverted light `u` (`retune_pulses`'s `unconverted`) of the beams `keep`: each
# color's fraction per beam, if it has one, at those beams.
_kept_unconverted(::Nothing, keep) = nothing
_kept_unconverted(u::NamedTuple, keep) = map(f -> f isa AbstractVector ? f[keep] : f, u)

"""
    retune_to_reference!(sess::DesignSession, table, out; by = :quad, cap = 1.3, stay = 0.0,
                         unit = :quad, target = :reference, level = nothing, on = :auto, load = true)

Retune the enabled beams to the reference shot whose pulses are the table `table`
(`read_pulse_table`: a `time_ns` column and one per group, quad or beam; a cone's or a
quad's column the power per `unit`, [`PULSE_UNITS`](@ref), as the table written), as
[`retune_pulses`](@ref) with the drive's scenario (which must be current), its spot
radius and the active foil ([`active_foil`](@ref)), to the reference's drive shape
(`target = :reference`) or a uniform drive (`:uniform`), at the reference's level or at
`level` (TW/mm², `retune_pulses`'s), and write the new pulses to `out`
(`write_pulse_table`), on a worker when `spawn[]` is set. With `load` the beams then fire
the new pulses: they become the Time panel's (as [`load_pulses!`](@ref) of `out`, before
the CSV's rounding), an edit that undo takes back. The disabled beams are the ones
the shot loses. The retune holds the observed surface's drive (a capsule's modes, a flat
face's flux) when the drive has a radiation model, else the active foil's irradiance;
`on = :observed` or `:foil` chooses (`retune_pulses`'s `on`) when the design has both.
The unconverted 1ω and 2ω light takes part as in the drive
([`set_unconverted_ratio!`](@ref); `retune_pulses`'s `unconverted`), with the prepulse
([`set_prepulse!`](@ref); `retune_pulses`'s `prepulse`), in the reference shot and in
the retune, on the beams whose 1ω or 2ω light takes part ([`participates`](@ref)), for a
foil as each color's absorbed irradiance on its spot. Progress and the outcome go to
`drive.retune_running`/`retune_status`: for a capsule the worst-time mode 1 before and
after, for a face or a foil the worst RMS nonuniformity and the worst departure of the
mean from the reference's (from `level` when it is set; for a capsule the Tᵣ's likewise).
A drive with neither is refused with a status message.
"""
function retune_to_reference!(sess::DesignSession, table::AbstractString, out::AbstractString;
                              by::Symbol = :quad, cap::Real = 1.3, stay::Real = 0.0, unit::Symbol = :quad,
                              target::Symbol = :reference, level::Union{Nothing,Real} = nothing,
                              on::Symbol = :auto, load::Bool = true)
    dr = sess.drive
    dr.retune_running[] && return sess
    if !_drive_current(sess)
        dr.retune_status[] = "compute the drive first: the retune reuses its operator and radiation model"
        return sess
    end
    s = dr.scenario[]
    foil = active_foil(sess.target)
    kind = _retune_kind(s, foil, on)
    kind === nothing && (dr.retune_status[] = on === :auto ? _RETUNE_NONE :
                                              "the design has no $(on === :foil ? "foil" : "observed surface") to hold"; return sess)
    # the unconverted light of every beam (the reference fires them all)
    unconverted = _history_unconverted(sess, eachindex(sess.beams.beams[]))
    prepulse = _history_prepulse(sess, eachindex(sess.beams.beams[]))
    pulses = read_pulse_table(table)
    beams = sess.beams.beams[]
    enabled = copy(sess.beams.enabled[])
    out = String(out)
    spot_radius = dr.spot_radius[]
    dr.retune_running[] = true
    dr.retune_status[] = "retuning to $(basename(table))…"
    runner = _runner(sess.spawn, sess.deliver)
    wall = Tuple(wall_names(sess.target))
    runner(() -> retune_pulses(s, beams, enabled, pulses; by, cap, stay, unit, target, level, wall, unconverted,
                               prepulse, foil, spot_radius, on), function (res)
        dr.retune_running[] = false
        if res.err !== nothing
            dr.retune_status[] = "retune failed: " * _errmsg(res.err)
            return
        end
        v = res.value
        write_pulse_table(out, v.pulses)
        # the time panel plots the retune
        bt = sess.beams
        nodes = sort(unique(reduce(vcat, [p.t for p in values(pulses)])))
        sess.time.retune[] = (; label = basename(table), v.kind, v.target, v.level,
                              solved = isfinite.(v.retune.residual), by, unit, pulses = v.pulses, reference_pulses = pulses,
                              reference_by = _pulse_key(beams, pulses), nodes, beams, enabled,
                              table = (; cluster = copy(bt.cluster[]), colors = cluster_colors(bt)),
                              before = v.before, after = v.after, reference = v.reference)
        # the beams fire the new pulses (the table's beams unchanged while the worker ran)
        loaded = load && length(bt.beams[]) == length(beams) &&
                 all(b -> haskey(v.pulses, Briareus._beam_key(by)(b)), bt.beams[])
        loaded && _set_pulses!(sess, _table_pulses(v.pulses, by, basename(out); unit))
        set_chart!(sess.time, :retune)
        n = length(v.retune.keys)
        aim = v.target === :uniform ? "uniform, " : ""
        solved = isfinite.(v.retune.residual)
        if v.kind === :capsule
            ε1(h) = 100 * maximum(mode_amplitude(h.ylm[:, k], 1) for k in eachindex(h.time))
            # the Tᵣ against the reference's, or against the level's where the retune solves
            ΔT(h) = v.level === nothing ? maximum(abs.(h.Tr .- v.reference.Tr)) :
                    _finite_max(abs.(h.Tr[solved] .- radiation_temperature(v.level)))
            dr.retune_status[] = @sprintf("%d %s powers (%s%s): worst ε₁ %.2f%% → %.2f%%, |ΔTᵣ| %.1f → %.1f eV; wrote %s",
                                          n, by, aim, _level_text(v), ε1(v.before), ε1(v.after), ΔT(v.before),
                                          ΔT(v.after), basename(out))
        else
            # (over the steps where the reference lights the face or the spot)
            rms(f) = 100 * _finite_max(f.rms)
            Δ(f) = 100 * _finite_max(abs.(_mean_departure(f, v.level, solved, v.reference)))
            dr.retune_status[] = @sprintf("%d %s powers (%s%s): %s worst RMS %.2f%% → %.2f%% (reference %.2f%%), |Δmean| %.2f%% → %.2f%%; wrote %s",
                                          n, by, aim, _level_text(v),
                                          v.kind === :face ? "face" : "foil spot", rms(v.before), rms(v.after),
                                          rms(v.reference), Δ(v.before), Δ(v.after), basename(out))
        end
        sess.status[] = "wrote the retuned pulses to " * basename(out) *
                        (loaded ? "; the beams fire them (Undo takes them back)" : "")
    end)
    return sess
end

# ---------------------------------------------------------------- optimize

"""
    start_optimize!(sess::DesignSession; spawn = true)

Optimize the enabled beams on the current scenario (see
[`start_optimize!(::OptimizeState, …)`](@ref start_optimize!)).
"""
function start_optimize!(sess::DesignSession; spawn::Bool = true)
    st = sess.optimize
    # a group's power is one power per beam of it (`GroupPointing`'s `:power`): one power
    # per quad only where the group's quads are held alike
    if :power in st.fields[] && st.by[] === :group
        mixed = _mixed_quads(sess.beams)
        mixed === nothing ||
            (st.status[] = "$mixed holds quads as one beam and as four: split or merge them all, " *
                           "or free the power per quad or per beam"; return sess)
    end
    if !_drive_current(sess)
        st.status[] = "compute the drive first (Compute ▸ Drive)"
        return sess
    end
    s = sess.drive.scenario[]
    # objectives the drive cannot serve: say why and what to do, instead of a failed job
    foil = active_foil(sess.target)
    if st.objective[] === :irradiance
        _has_part(s.target, foil) ||
            (st.status[] = "the foil objective needs a foil: a foil target's, or a part with the role " *
                           "Backlighter foil (Properties ▸ Role)"; return sess)
    elseif s.model === nothing
        st.status[] = (st.objective[] === :face ? "the face objective" : "capsule symmetry") *
                      " needs the x-ray drive, and this target's drive is laser only" * _radiation_note(s.target)
        return sess
    elseif st.objective[] === :face
        _observed_flat(s) ||
            (st.status[] = "the face objective needs a flat observed face: a disk or rectangle part with " *
                           "the role Observed (Properties ▸ Role)"; return sess)
    elseif isempty(s.observed) || _observed_flat(s)
        _observed_flat(s) &&
            (st.status[] = "the observed surface is a flat face: choose the objective Observed face uniformity";
             return sess)
        st.status[] = "capsule symmetry needs a hohlraum with a capsule, or a sphere part with the " *
                      "role Observed (Properties ▸ Role)" *
                      (isempty(s.observed) ? "" : ": the observed surface is a flat face")
        return sess
    end
    bt = sess.beams
    idx = findall(bt.enabled[])
    start_optimize!(st, s, bt.beams[][idx], idx;
                    spot_radius = sess.drive.spot_radius[], foil,
                    observed = observed_name(sess.target), clusters = bt.cluster[][idx],
                    frames = beam_frames(bt)[idx], key = _optimize_key(sess), deliver = sess.deliver[], spawn)
    return sess
end

cancel_optimize!(sess::DesignSession) = cancel_optimize!(sess.optimize)

"""
    accept_optimum!(sess::DesignSession)

Put the optimized pointing into the beam table (the enabled beams it was computed
for); the drive then re-evaluates.
"""
function accept_optimum!(sess::DesignSession)
    r = sess.optimize.result[]
    r === nothing && throw(ArgumentError("run the optimization first"))
    bt = sess.beams
    bs = copy(bt.beams[])
    all(i -> i <= length(bs), r.indices) && all(k -> bs[r.indices[k]].port == r.beams[k].port,
                                                eachindex(r.indices)) ||
        throw(ArgumentError("the beam table changed since the optimization; run it again"))
    bs[r.indices] = r.beams
    _set_aims!(bt, bs)
    r.applied = true
    notify(sess.optimize.result)
    sess.status[] = "applied the optimized pointing to $(length(r.indices)) beams"
    return sess
end

# ---------------------------------------------------------------- review

"""
    start_clearance!(sess::DesignSession; spawn = true, batch = review_batch())

Check the enabled beams' clearance from the selected hardware (see
[`start_clearance!(::ReviewState, …)`](@ref start_clearance!)).
"""
start_clearance!(sess::DesignSession; spawn::Bool = true, batch::Integer = review_batch()) =
    (await_target!(sess.target);
     start_clearance!(sess.review, sess.target.target[], active_beams(sess.beams);
                      deliver = sess.deliver[], spawn, batch); sess)

cancel_clearance!(sess::DesignSession) = cancel_clearance!(sess.review)

"""
    start_talis!(sess::DesignSession; spawn = true)

Run the TaLIS checks on the enabled beams with the clearance check's error stack (see
[`start_talis!(::TaLISState, …)`](@ref start_talis!)).
"""
start_talis!(sess::DesignSession; spawn::Bool = true) =
    (await_target!(sess.target);
     start_talis!(sess.talis, sess.review, sess.target.target[], active_beams(sess.beams);
                  deliver = sess.deliver[], spawn); sess)

cancel_talis!(sess::DesignSession) = (cancel_talis!(sess.talis); sess)

"""
    run_task!(sess::DesignSession, task::Symbol; spawn = true) -> sess

Start analysis `task`, one of [`TASKS`](@ref): `:drive` builds the scenario
([`build_scenario!`](@ref)), `:optimize`, `:clearance`, `:sweep` and `:talis` start their
jobs ([`start_optimize!`](@ref), [`start_clearance!`](@ref), [`start_sweep!`](@ref),
[`start_talis!`](@ref)); the clearance and TaLIS checks color the beams by their reports.
The model tree's *Run* and the Compute commands go through it.
"""
function run_task!(sess::DesignSession, task::Symbol; spawn::Bool = true)
    task in TASKS || throw(ArgumentError("unknown analysis :$task (one of $TASKS)"))
    if task === :drive
        build_scenario!(sess)
    elseif task === :optimize
        start_optimize!(sess; spawn)
    elseif task === :clearance
        start_clearance!(sess; spawn)
        _assign!(sess.view.beam_color, :clearance)
    elseif task === :sweep
        start_sweep!(sess; spawn)
    else
        start_talis!(sess; spawn)
        _assign!(sess.view.beam_color, :talis)
    end
    return sess
end

"""
    review_stale(sess::DesignSession) -> Bool

Whether the target, the enabled beams or a setting changed since the clearance or TaLIS
check.
"""
function review_stale(sess::DesignSession)
    t, bs = sess.target.target[], active_beams(sess.beams)
    return review_stale(sess.review, t, bs) || talis_stale(sess.talis, sess.review, t, bs)
end

"""
    load_preset!(sess::DesignSession, preset::Symbol)

Load a beam preset pointed for the current target (see [`preset_beams`](@ref)).
"""
load_preset!(sess::DesignSession, preset::Symbol) =
    (load_preset!(sess.beams, preset; target = sess.target); sess)

"""
The facility each starter is a mock-up for, and the beam preset a new target of it gets
when the beams are another facility's ([`new_target!`](@ref)): the hohlraum, frustraum
and halfraum starters are at NIF scale (the NIF quads), the backlighter foil is the
OMEGA example's (the OMEGA backlighter beams). No OMEGA-scale hohlraum starter yet: it
needs published dimensions.
"""
const STARTER_PRESETS = Dict(:hohlraum => :nif_quads, :frustraum => :nif_quads, :shock_tube => :nif_quads,
                             :opacity => :nif_quads, :foil => :omega_backlighter)

"""
    new_target!(sess::DesignSession, starter = :hohlraum) -> Bool

Start a new target from `starter` ([`new_target!`](@ref) on the session's target). An
optimization of the previous target is cleared. The starters belong to a facility
([`STARTER_PRESETS`](@ref)): when the beams are another facility's preset, or a preset
of the same facility that does not drive this kind of target (the OMEGA hohlraum beams
on a foil), the starter's own preset replaces them, pointed for the new target; a
hohlraum preset whose beams do not aim in the new target's entrance holes as
[`preset_pointing`](@ref) aims them (a halfraum pair after a hohlraum) is loaded again
for the new target. Imported beams are kept, as are the beams of an empty target
(`:none`).
"""
function new_target!(sess::DesignSession, starter::Symbol = :hohlraum)
    new_target!(sess.target, starter) || return false
    clear_optimize!(sess.optimize)
    p = sess.beams.preset[]
    want = get(STARTER_PRESETS, starter, nothing)
    (want === nothing || p === :custom) && return true
    drives = starter === :foil ? p === :omega_backlighter : p in HOHLRAUM_PRESETS
    if facility_of(p) !== facility_of(want) || !drives
        load_preset!(sess, want)
    elseif p in HOHLRAUM_PRESETS && sess.beams.aim_frame[] != preset_pointing(p; target = sess.target).aim_frame
        load_preset!(sess, p)
    end
    return true
end

"""
    remove_nodes!(sess::DesignSession, ids = selected_ids(sess.view)) -> NamedTuple

Remove the model-tree nodes `ids` that can be removed, in one edit: parts, assemblies
(with their contents) and user frames ([`remove_parts!`](@ref)) and diagnostic windows
([`remove_window!`](@ref)); other nodes (pieces, beams, analyses) are skipped. Returns
`(; removed, replaced, skipped, ok)`: the names removed, the parts that stayed and were
placed in another frame where they are ([`replaced_parts`](@ref)), the ids skipped, and
whether the design accepted the removal (`sess.target.error` says why not). The axes of
a removed frame are no longer drawn, and the selection goes to the target.
"""
function remove_nodes!(sess::DesignSession, ids = selected_ids(sess.view))
    ts = sess.target
    names, windows, skipped = Symbol[], Tuple{Symbol,Int}[], String[]
    for id in ids
        kind, key = selection_kind(id)
        if kind in (:part, :assembly, :frame) && part_index(ts, Symbol(key)) !== nothing
            push!(names, Symbol(key))
        elseif kind === :window && window_of(ts, Symbol(key)) !== nothing
            i, k = window_of(ts, Symbol(key))
            push!(windows, (ts.parts[][i].name, k))
        else
            push!(skipped, String(id))
        end
    end
    unique!(names)
    removed = String[]
    ok = true
    # the windows of the parts that stay, the highest number first (the others keep theirs)
    for (part, k) in sort!(unique!(windows); by = w -> (w[1], -w[2]))
        part in names && continue
        label = String(piece_name(part, window_name(k)))
        remove_window!(ts, k; part) ? push!(removed, label) : (ok = false)
    end
    replaced = Symbol[]
    if !isempty(names)
        replaced = replaced_parts(ts, names)
        if remove_parts!(ts, names)
            append!(removed, String.(names))
            foreach(n -> set_triad!(sess.view, n, false), names)
        else
            ok = false
            empty!(replaced)
        end
    end
    isempty(removed) || select!(sess.view, "target")
    return (; removed, replaced, skipped, ok)
end

"""
    insert_starter!(sess::DesignSession, kind; frame = :tcc) -> Bool

Insert starter `kind` ([`insert_starter!`](@ref) on the session's target). When a hohlraum
preset is active ([`HOHLRAUM_PRESETS`](@ref)) and no beam aimed in an entrance hole
before, the beams held in `:tcc` then aim in the new entrance holes
([`aim_in_lehs!`](@ref)); beams that already drive a hohlraum stay with it.
"""
function insert_starter!(sess::DesignSession, kind::Symbol; frame::Symbol = :tcc)
    bt = sess.beams
    driving = any(f -> occursin(r"\.leh_(top|bottom)$", String(f)), bt.aim_frame[])
    insert_starter!(sess.target, kind; frame) || return false
    driving || aim_in_lehs!(bt, sess.target)
    return true
end

# Abandon jobs in flight: their results are dropped when they arrive (a closing
# window restores inline delivery).
function _abandon_jobs!(sess::DesignSession)
    sess.mesh_generation[] += 1
    _assign!(sess.meshing, "")
    # parts still remeshing stay `pending`, and a target being built `building`: the
    # drive and the checks finish them inline (`await_target!`); capsule fits are dropped
    sess.target.mesh_job[] += 1
    isempty(sess.target.fitting[]) ||
        (sess.status[] = "the capsule fit of " * join(sort!(String.(collect(sess.target.fitting[]))), ", ") *
                         " was abandoned when the window closed: give the role again")
    _cancel_fits!(sess.target)
    sess.build_queued[] = false
    sess.drive.generation[] += 1
    sess.drive.eval_generation[] += 1
    sess.drive.building[] = false
    sess.drive.evaluating[] = false
    sess.sweep.cancel[] = true
    for (job, cancel!) in ((sess.optimize, cancel_optimize!), (sess.review, cancel_clearance!),
                           (sess.talis, cancel_talis!), (sess.sweep, _ -> nothing))
        job.running[] || continue
        cancel!(job)                       # the worker stops at its next check
        job.running[] = false
        job.status[] = "abandoned when the window closed"
    end
    tm = sess.time
    tm.queued[] = false
    tm.check_queued[] = false
    tm.running[] && (tm.generation[] += 1; tm.cancel[] = true; tm.running[] = false;
                     tm.status[] = "abandoned when the window closed")
    tm.checking[] && (tm.check_generation[] += 1; tm.checking[] = false;
                      tm.status[] = "abandoned when the window closed")
    return sess
end
