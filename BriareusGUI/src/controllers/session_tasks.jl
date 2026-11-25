# This file covers the optimization and the checks of a session, the tasks, and new targets,
# starters and removed nodes.

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
on a foil), the preset of the starter replaces them, pointed for the new target; a
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

# Abandon jobs in flight. Their results are dropped when they arrive, and a closing window
# restores inline delivery.
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
