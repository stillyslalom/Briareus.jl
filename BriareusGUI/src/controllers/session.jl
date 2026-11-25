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

# The parts definition `def`'s view names, `key => name`, with the session `vp`
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
        # The role holders below wait for the target of the new names, since a renamed
        # wall is not a new wall.
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
    # a beam selection whose beams went away, after a preset, a split or a regrouping
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
        # A renamed part is committed before the target of its new names, so the roles
        # follow the rename itself; a comparison with the target of the old names would
        # not.
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
        # For walls and the observed surface, a role that the checks give a new component
        # by default, as after a rebuild, is left as it is.
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

# After an undo or redo, the optimization result is applied when the beam table holds its
# optimized beams.
function _restored_optimum!(sess::DesignSession)
    r = sess.optimize.result[]
    r === nothing && return sess
    bs = sess.beams.beams[]
    held = all(i -> i <= length(bs), r.indices) && all(k -> bs[r.indices[k]] == r.beams[k], eachindex(r.indices))
    held == r.applied || (r.applied = held; notify(sess.optimize.result))
    return sess
end
