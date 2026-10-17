# The Briareus window's state: the target, the beams, the analyses on them (drive,
# optimization, clearance, TaLIS) and what the window shows, the links between them (a
# target or port change makes the scenario stale; pointing edits re-evaluate it), and
# the background-job plumbing (`spawn`, `deliver`) shared by every job. Nothing here
# imposes an order: any analysis can be started at any time, and each reports whether
# it is current.

"""
Analyses of a design, in the model tree's order: `:drive` (view factors, laser
deposition and x-ray drive), `:optimize` (pointing optimization), `:clearance`
(worst-case beam clearance) and `:talis` (NIF target-configuration rules).
"""
const TASKS = (:drive, :optimize, :clearance, :talis)

const TASK_LABELS = Dict(:drive => "Drive (view factors)", :optimize => "Pointing optimization",
                         :clearance => "Beam clearance", :talis => "TaLIS rules")

"""
    DesignSession(; kind = :hohlraum, resolution = :standard, preset = :nif_quads,
                  beams = nothing, deliver = f -> f())

State of the Briareus window: `target` ([`TargetSetup`](@ref)), `beams`
([`BeamTable`](@ref)), the analyses `drive` ([`DriveState`](@ref)), `optimize`
([`OptimizeState`](@ref)), `review` ([`ReviewState`](@ref), the clearance check) and
`talis` ([`TaLISState`](@ref)), `view` ([`ViewState`](@ref): coloring, beam drawing,
selection, model tree) and `status` (the last message).

The table starts from `beams` (a vector of `Beam`s) or a [`BEAM_PRESETS`](@ref) entry
pointed for the target. Changing the target kind swaps a NIF preset for the OMEGA
backlighter set (and back), unless the beams were imported.

Background work: `spawn[]` (`false` by default; a window sets it) runs the scenario
build, its evaluations and the VisRad kit on worker tasks, handing their results to
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
    view::ViewState
    status::Observable{String}
    deliver::Base.RefValue{Any}
    spawn::Base.RefValue{Bool}
end

function DesignSession(; kind::Symbol = :hohlraum, resolution::Symbol = :standard,
                       preset::Symbol = kind === :foil ? :omega_backlighter : :nif_quads,
                       beams = nothing, deliver = f -> f())
    ts = TargetSetup(; kind, resolution)
    bt = beams === nothing ?
         BeamTable(preset_beams(preset; target_kind = kind, params = ts.params[]); preset) :
         BeamTable(beams)
    sess = DesignSession(ts, bt, DriveState(), OptimizeState(), ReviewState(), TaLISState(),
                         ViewState(), Observable(""), Ref{Any}(deliver), Ref(false))
    default_selection!(sess.review, ts.target[])
    default_selection!(sess.talis, ts.target[])
    _connect!(sess)
    return sess
end

function Base.show(io::IO, s::DesignSession)
    print(io, "DesignSession(", target_summary(s.target), "; ", beams_summary(s.beams), ")")
end

function _connect!(sess::DesignSession)
    ts, bt, dr = sess.target, sess.beams, sess.drive
    on(ts.kind) do kind
        p = bt.preset[]
        if kind === :foil && p in (:nif_quads, :nif_beams, :omega_hohlraum)
            load_preset!(bt, :omega_backlighter; target_kind = kind, params = ts.params[])
        elseif kind === :hohlraum && p === :omega_backlighter
            load_preset!(bt, :nif_quads; target_kind = kind, params = ts.params[])
        end
    end
    on(ts.target) do t
        default_selection!(sess.review, t)
        default_selection!(sess.talis, t)
        # forget hidden components and a selection that no longer exist
        names = Set(component_names(t))
        h = filter(in(names), sess.view.hidden[])
        h == sess.view.hidden[] || (sess.view.hidden[] = h)
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
    onany((_...) -> _drive_inputs_changed!(sess), bt.beams, bt.enabled)
    # a rejected edit's message belongs to what was being edited
    on(sess.view.selection) do _
        isempty(ts.error[]) || (ts.error[] = "")
        isempty(bt.error[]) || (bt.error[] = "")
    end
    on(_ -> _request_evaluation!(sess), dr.spot_radius)
    return sess
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
    task_status(sess::DesignSession, task) -> (state, summary)

`state` is `:todo`, `:ok`, `:attention`, or `:busy`, with a one-line summary for the
model tree (`task` one of [`TASKS`](@ref)).
"""
function task_status(sess::DesignSession, task::Symbol)
    task in TASKS || throw(ArgumentError("unknown analysis :$task"))
    if task === :drive
        dr = sess.drive
        dr.building[] && return (:busy, "computing view factors… " * _seconds(time() - dr.started[]))
        dr.scenario[] === nothing && return (:todo, isempty(dr.status[]) ? "not computed" : dr.status[])
        drive_stale(sess) && return (:attention, "target or beam ports changed: compute again")
        dr.evaluating[] && return (:busy, "evaluating…")
        dr.result[] === nothing && return (:attention, isempty(dr.status[]) ? "no evaluation" : dr.status[])
        return (:ok, drive_summary(dr.result[]))
    elseif task === :optimize
        st = sess.optimize
        st.running[] && return (:busy, optimize_progress(st))
        r = st.result[]
        r === nothing && return (:todo, isempty(st.status[]) ? "not run" : st.status[])
        txt = @sprintf("objective %.3g → %.3g", r.initial_value, r.value)
        r.applied && return (:ok, txt * " · applied")
        return (:attention, txt * " · not applied")
    elseif task === :clearance
        rv = sess.review
        rv.running[] && return (:busy, "checking: $(rv.progress[][1]) of $(rv.progress[][2]) beams")
        rv.report[] === nothing && return (:todo, isempty(rv.status[]) ? "not run" : rv.status[])
        txt = review_brief(rv)
        review_stale(rv, sess.target.target[], active_beams(sess.beams)) &&
            return (:attention, "beams or target changed since: " * txt)
        bad = rv.report[].canceled || any(x -> x.status !== :pass, rv.report[].rows)
        return (bad ? :attention : :ok, txt)
    else
        tl = sess.talis
        tl.running[] && return (:busy, "checking… " * _seconds(time() - tl.started[]))
        tl.report[] === nothing && return (:todo, isempty(tl.status[]) ? "not run" : tl.status[])
        txt = talis_brief(tl)
        talis_stale(tl, sess.target.target[], active_beams(sess.beams)) &&
            return (:attention, "beams or target changed since: " * txt)
        bad = any(x -> x.status !== :pass, tl.report[].rows)
        return (bad ? :attention : :ok, txt)
    end
end

# ---------------------------------------------------------------- drive

"""
    drive_stale(sess::DesignSession) -> Bool

Whether the scenario was built for another target or other beam ports (it must be
computed again). `false` before the first build.
"""
function drive_stale(sess::DesignSession)
    k = sess.drive.key[]
    k === nothing && return false
    return k[1] !== sess.target.target[] || k[2] != _ports(sess.beams.beams[])
end

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
seconds), on a worker when `spawn[]` is set, then evaluate the current pointing.
"""
function build_scenario!(sess::DesignSession)
    dr = sess.drive
    t = sess.target.target[]
    beams = sess.beams.beams[]
    isempty(beams) && (dr.status[] = "add beams first"; return sess)
    g = (dr.generation[] += 1)
    dr.building[] = true
    dr.started[] = time()
    dr.status[] = "computing view factors…"
    runner = _runner(sess.spawn, sess.deliver)
    t0 = time()
    # free parts lie outside the hohlraum's enclosure: its rows no longer sum to one
    closed = sess.target.kind[] === :hohlraum && isempty(sess.target.parts[])
    runner(() -> _build_scenario(t, beams; closed), function (out)
        g == dr.generation[] || return
        dr.building[] = false
        if out.err !== nothing
            dr.status[] = "failed: " * _errmsg(out.err)
            return
        end
        s, lookup = out.value
        dr.build_seconds[] = time() - t0
        dr.key[] = (t, _ports(beams))
        dr.lookup[] = lookup
        dr.scenario[] = s
        dr.status[] = @sprintf("built in %s: %d elements, %d beams", _seconds(dr.build_seconds[]),
                               nelements(t), length(beams))
        _request_evaluation!(sess)
    end)
    return sess
end

# Evaluate the table's current pointing on the scenario (latest request wins).
function _request_evaluation!(sess::DesignSession)
    dr = sess.drive
    _drive_current(sess) || return sess
    s, lookup, radius = dr.scenario[], dr.lookup[], dr.spot_radius[]
    beams = display_beams(sess.beams)
    g = (dr.eval_generation[] += 1)
    dr.evaluating[] = true
    runner = _runner(sess.spawn, sess.deliver)
    runner(() -> _evaluate_drive(s, beams, lookup, radius), function (out)
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
    write_visrad_kit!(sess::DesignSession, dir)

Write a VisRad benchmark kit for the current target and beam table into `dir` with
`write_visrad_kit` (SST pointing, beam summary, the canonical beam table, a README with
the geometry and VisRad steps, and Briareus's results to compare against), on a worker
when `spawn[]` is set. It reuses the drive's deposition operator and radiation model,
so the drive must be computed and current. Disabled beams are written at zero power.
The capsule (or a foil target's main sphere) is the observed surface. Progress and
the outcome go to `drive.kit_running`/`drive.kit_status`.
"""
function write_visrad_kit!(sess::DesignSession, dir::AbstractString)
    dr = sess.drive
    dr.kit_running[] && return sess
    if !_drive_current(sess)
        dr.kit_status[] = "compute the drive first: the kit reuses its view factors"
        return sess
    end
    s = dr.scenario[]
    t = sess.target.target[]
    beams = display_beams(sess.beams)
    names = component_names(t)
    observe = :capsule in names ? :capsule : :main in names ? :main : :none
    dir = String(dir)
    dr.kit_running[] = true
    dr.kit_status[] = "writing the VisRad kit…"
    notes = "Written by BriareusGUI: " * target_summary(sess.target) * "; " * beams_summary(sess.beams)
    runner = _runner(sess.spawn, sess.deliver)
    runner(() -> write_visrad_kit(dir, t, beams; op = s.deposition, model = s.model, observe, notes),
           function (out)
               dr.kit_running[] = false
               if out.err === nothing
                   n = length(out.value.files)
                   dr.kit_status[] = "wrote $n files to $dir"
                   sess.status[] = "wrote the VisRad kit to " * basename(rstrip(dir, ['/', '\\']))
               else
                   dr.kit_status[] = "kit failed: " * _errmsg(out.err)
               end
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
    if !_drive_current(sess)
        st.status[] = "compute the drive first (Compute ▸ Drive)"
        return sess
    end
    bt = sess.beams
    idx = findall(bt.enabled[])
    start_optimize!(st, sess.drive.scenario[], bt.beams[][idx], idx;
                    spot_radius = sess.drive.spot_radius[], clusters = bt.cluster[][idx],
                    deliver = sess.deliver[], spawn)
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
    bt.beams[] = bs
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
    (start_clearance!(sess.review, sess.target.target[], active_beams(sess.beams);
                      deliver = sess.deliver[], spawn, batch); sess)

cancel_clearance!(sess::DesignSession) = cancel_clearance!(sess.review)

"""
    start_talis!(sess::DesignSession; spawn = true)

Run the TaLIS checks on the enabled beams with the clearance check's error stack (see
[`start_talis!(::TaLISState, …)`](@ref start_talis!)).
"""
start_talis!(sess::DesignSession; spawn::Bool = true) =
    (start_talis!(sess.talis, sess.review, sess.target.target[], active_beams(sess.beams);
                  deliver = sess.deliver[], spawn); sess)

cancel_talis!(sess::DesignSession) = (cancel_talis!(sess.talis); sess)

"""
    review_stale(sess::DesignSession) -> Bool

Whether the target or the enabled beams changed since the clearance or TaLIS check.
"""
function review_stale(sess::DesignSession)
    t, bs = sess.target.target[], active_beams(sess.beams)
    return review_stale(sess.review, t, bs) || talis_stale(sess.talis, t, bs)
end

"""
    load_preset!(sess::DesignSession, preset::Symbol)

Load a beam preset pointed for the current target (see [`preset_beams`](@ref)).
"""
load_preset!(sess::DesignSession, preset::Symbol) =
    (load_preset!(sess.beams, preset; target_kind = sess.target.kind[],
                  params = sess.target.params[]); sess)

# Abandon jobs in flight: their results are dropped when they arrive (a closing
# window restores inline delivery).
function _abandon_jobs!(sess::DesignSession)
    sess.drive.generation[] += 1
    sess.drive.eval_generation[] += 1
    sess.drive.building[] = false
    sess.drive.evaluating[] = false
    for (job, cancel!) in ((sess.optimize, cancel_optimize!), (sess.review, cancel_clearance!),
                           (sess.talis, cancel_talis!))
        job.running[] || continue
        cancel!(job)                       # the worker stops at its next check
        job.running[] = false
        job.status[] = "abandoned when the window closed"
    end
    sess.drive.kit_running[] && (sess.drive.kit_running[] = false;
                                 sess.drive.kit_status[] = "abandoned when the window closed")
    return sess
end
