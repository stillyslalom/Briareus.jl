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

"""
    session_tasks(sess) -> Tuple

The analyses that apply to the session's facility: all of [`TASKS`](@ref) for NIF,
without the TaLIS rules (NIF's) elsewhere.
"""
session_tasks(sess) = sess.beams.facility[] === :nif ? TASKS : filter(!=(:talis), TASKS)

const TASK_LABELS = Dict(:drive => "Drive (view factors)", :optimize => "Pointing optimization",
                         :clearance => "Beam clearance", :talis => "TaLIS rules")

"""
    EditHistory(; limit = 200)

Undo and redo stacks of design snapshots ([`session_state`](@ref) without the view:
coloring, selection and the camera are not edits). [`record!`](@ref) takes a snapshot
when the design has changed since the last one; [`undo!`](@ref) and [`redo!`](@ref)
step through them.
"""
mutable struct EditHistory
    undo::Vector{Dict{String,Any}}
    redo::Vector{Dict{String,Any}}
    current::Union{Nothing,Dict{String,Any}}
    limit::Int
end
EditHistory(; limit::Integer = 200) = EditHistory(Dict{String,Any}[], Dict{String,Any}[], nothing, limit)

"""
    DesignSession(; kind = :hohlraum, resolution = :standard, preset = :nif_quads,
                  beams = nothing, deliver = f -> f())

State of the Briareus window: `target` ([`TargetSetup`](@ref)), `beams`
([`BeamTable`](@ref)), the analyses `drive` ([`DriveState`](@ref)), `optimize`
([`OptimizeState`](@ref)), `review` ([`ReviewState`](@ref), the clearance check) and
`talis` ([`TaLISState`](@ref)), `view` ([`ViewState`](@ref): coloring, beam drawing,
selection, model tree), `status` (the last message) and `history` (the
[`EditHistory`](@ref) of [`undo!`](@ref) and [`redo!`](@ref)). Designs are saved and
loaded with [`save_session`](@ref) and [`load_session!`](@ref).

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
    history::EditHistory
    visrad_source::Base.RefValue{String}   # the imported VisRad workspace (export template)
end

function DesignSession(; kind::Symbol = :hohlraum, resolution::Symbol = :standard,
                       preset::Symbol = kind === :foil ? :omega_backlighter : :nif_quads,
                       beams = nothing, deliver = f -> f())
    ts = TargetSetup(; kind, resolution)
    bt = beams === nothing ?
         BeamTable(preset_beams(preset; target_kind = kind, params = ts.params[]); preset) :
         BeamTable(beams)
    sess = DesignSession(ts, bt, DriveState(), OptimizeState(), ReviewState(), TaLISState(),
                         ViewState(), Observable(""), Ref{Any}(deliver), Ref(false), EditHistory(),
                         Ref(""))
    default_selection!(sess.review, ts.target[])
    default_selection!(sess.talis, ts.target[])
    _connect!(sess)
    record!(sess)
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
        names = Set(all_component_names(ts))
        h = filter(in(names), sess.view.hidden[])
        h == sess.view.hidden[] || (sess.view.hidden[] = h)
        c = selected_component(sess)
        c === nothing || c in names || select!(sess.view, "target")
        # beams that follow a component move with it
        any(!=(:none), bt.aim_ref[]) && follow_references!(bt, component_points(t))
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
    # the x-ray image follows the drive
    on(dr.result) do _
        dr.image_request[] === nothing || (dr.image[] = _render_image(dr))
    end
    return sess
end

"""
    component_points(t::Target) -> Dict{Symbol,Vec3}

Reference point of each component, for aiming beams at it: the area-weighted centroid
of its surface (the center of an entrance hole, a foil, a disk).
"""
function component_points(t::Target)
    acc = Dict{Symbol,Tuple{Vec3,Float64}}()
    for k in eachindex(t.tri_area)
        n = t.components[t.tri_comp[k]].name
        c, a = get(acc, n, (zero(Vec3), 0.0))
        A = t.tri_area[k]
        acc[n] = (c + A * (t.tri_v0[k] + (t.tri_e1[k] + t.tri_e2[k]) / 3), a + A)
    end
    return Dict(n => (a > 0 ? c / a : c) for (n, (c, a)) in acc)
end

"""
    aim_at_component!(sess::DesignSession, idx, name::Symbol) -> sess

Aim beams `idx` at component `name`'s reference point ([`component_points`](@ref)) and
keep them on it as the target changes (see [`set_aim_reference!`](@ref)); `:none`
releases them.
"""
function aim_at_component!(sess::DesignSession, idx, name::Symbol)
    name === :none && return (set_aim_reference!(sess.beams, idx, :none); sess)
    pts = component_points(sess.target.target[])
    haskey(pts, name) || throw(ArgumentError("no component $name in the calculations"))
    set_aim_reference!(sess.beams, idx, name, pts[name])
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
    else
        sess.drive.image[] = nothing
        sess.drive.building[] || build_scenario!(sess)
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
    # free parts lie outside the hohlraum's enclosure, and an excluded component leaves a
    # hole in it: its rows no longer sum to one
    closed = sess.target.kind[] === :hohlraum && isempty(sess.target.parts[]) &&
             sess.target.ghost[] === nothing
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
    import_visrad_workspace!(sess::DesignSession, path; include = o -> o.radiosity || o.laser)
        -> sess

Replace the design with a VisRad workspace's (`read_visrad_workspace`): a custom target
whose parts are the workspace's objects (`:visrad` parts: each keeps its place, albedo
and conversion; CAD windings are repaired and facets merged to the target's element
size), and the workspace's beams, in their clusters as groups, switched on where
VisRad delivers power. Objects for which `include` is false (by default those VisRad
leaves out of both its laser and radiosity calculations) are imported but left out of
the calculations, drawn as ghosts. The workspace becomes the template of
[`export_visrad_workspace`](@ref).
"""
function import_visrad_workspace!(sess::DesignSession, path::AbstractString;
                                  include = o -> o.radiosity || o.laser)
    ws = read_visrad_workspace(path)
    supported = [o for o in ws.objects if o.kind in (:mesh, :sphere, :rectangle, :cylindrical_hohlraum)]
    isempty(supported) && throw(ArgumentError("$(basename(path)) has no target components Briareus reads"))
    parts = [PartSpec(:visrad; name = Symbol(o.name), file = abspath(path), object = o.name,
                      values = Dict(:albedo => o.albedo, :conversion => o.conversion))
             for o in supported]
    ts = sess.target
    # the components of objects VisRad leaves out are left out of the calculations
    names = Set{Symbol}(n for o in supported if !include(o) for n in visrad_component_names(o))
    any(include, supported) ||
        throw(ArgumentError("VisRad includes none of the workspace's objects in its calculations"))
    _try_build!(ts; kind = :custom, windows = WindowSpec[], parts, excluded = names) ||
        throw(ArgumentError("the workspace's target was rejected: $(ts.error[])"))
    bs = ws.beams
    set_beams!(sess.beams, bs; preset = :custom, enabled = [b.power > 0 for b in bs],
               cluster = Symbol[b.group for b in bs])
    sess.beams.path[] = String(path)
    sess.visrad_source[] = abspath(path)
    sess.status[] = "imported $(basename(path)): $(length(parts)) objects, $(length(bs)) beams"
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
    export_results(sess::DesignSession, dir) -> Vector{String}

Write the current drive's results to the folder `dir` for postprocessing
(`Briareus.write_results`: per-element laser and x-ray fields, per-component and per-beam
power, the capsule flux's Legendre and Y_lm modes, a summary) and the design that
produced them (`session.briareus.toml`, see [`save_session`](@ref)). The drive must be
computed and current; disabled beams appear at zero power. Returns the files written.
"""
function export_results(sess::DesignSession, dir::AbstractString)
    dr = sess.drive
    _drive_current(sess) && dr.result[] !== nothing ||
        throw(ArgumentError("compute the drive first: the export writes its results"))
    res = dr.result[]
    ev = res.evaluation
    t = ev.scenario.target
    names = component_names(t)
    observe = :capsule in names ? :capsule : nothing
    files = write_results(dir, t, ev.deposition; sol = res.solution, beams = ev.effective,
                          observe = res.solution === nothing ? nothing : observe)
    push!(files, save_session(joinpath(dir, "session.briareus.toml"), sess))
    sess.status[] = "exported the results to $(basename(normpath(dir)))"
    return files
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
    retune_pulses(s::Scenario, beams, enabled, pulses; by = :quad, cap = 1.3, stay = 0.0)
        -> NamedTuple

Retune the enabled beams' powers to the drive of a reference shot in which every beam of
`beams` (the scenario's, in its order) fired `pulses` (power per beam against time, keyed
by group, quad or beam name, e.g. `read_pulse_table`): at each time of the pulses where
the reference fires, the capsule flux's modes through l = 4 and its level are matched
(`Briareus.retune_drive`) by one power per `by` key (`:quad`, `:group` or `:beam`) within
`cap` times its reference power, preferring the reference powers with weight `stay`.
Disabled beams are the dropped ones. Returns the new `pulses` (by `by`, zero where the
reference is), the `retune`, and drive histories of the `reference`, of the enabled
beams on the reference pulses (`before`) and on the new ones (`after`).
"""
function retune_pulses(s::Scenario, beams::AbstractVector{<:Beam}, enabled::AbstractVector{Bool},
                       pulses::AbstractDict; by::Symbol = :quad, cap::Real = 1.3, stay::Real = 0.0)
    s.model !== nothing && :capsule in component_names(s.target) ||
        throw(ArgumentError("the retune matches the capsule's drive: compute the drive of a target with a capsule"))
    by in (:quad, :group, :beam) || throw(ArgumentError("retune by :quad, :group or :beam, not :$by"))
    cap > 0 || throw(ArgumentError("the power limit must be positive"))
    any(enabled) || throw(ArgumentError("every beam is disabled"))
    ref_by = _pulse_key(beams, pulses)
    refkey, key = Briareus._beam_key(ref_by), Briareus._beam_key(by)
    nodes = sort(unique(reduce(vcat, [p.t for p in values(pulses)])))
    times = filter(τ -> sum(b -> pulses[refkey(b)](τ), beams) > 0, nodes)
    isempty(times) && throw(ArgumentError("the reference pulses have no power"))
    op = s.deposition
    ref = drive_history(s.model, op, τ -> with_pulses(beams, pulses, τ; by = ref_by), times;
                        lmax = 4, ylm = true, keep_albedo = true)
    keep = findall(enabled)
    kept = beams[keep]
    opk = op[keep]
    # a key's reference power: its beams' mean (the table may be finer than the retune)
    members = Dict(k => [refkey(b) for b in kept if key(b) == k] for k in unique(key.(kept)))
    nominal(k, τ) = sum(pulses[r](τ) for r in members[k]) / length(members[k])
    limits(k, τ) = (0.0, cap * nominal(k, τ))
    r = retune_drive(s.model, opk, kept, ref; by, limits, stay,
                     nominal = stay > 0 ? nominal : nothing)
    before = drive_history(s.model, opk, τ -> with_pulses(kept, pulses, τ; by = ref_by), times;
                           lmax = 4, ylm = true)
    after = drive_history(s.model, opk, τ -> with_pulses(kept, r.pulses, τ; by), times;
                          lmax = 4, ylm = true)
    # the new pulses on the table's times, off where the reference is
    at(g, τ) = (k = findfirst(==(τ), times); k === nothing ? 0.0 : r.power[g, k])
    out = Dict{Symbol,PulseShape}(Symbol(k) => PulseShape(nodes, [at(g, τ) for τ in nodes])
                                  for (g, k) in enumerate(r.keys))
    return (; pulses = out, retune = r, reference = ref, before, after)
end

"""
    retune_to_reference!(sess::DesignSession, table, out; by = :quad, cap = 1.3, stay = 0.0)

Retune the enabled beams to the reference shot whose pulses are the table `table`
(`read_pulse_table`: a `time_ns` column and one per group, quad or beam), as
[`retune_pulses`](@ref) with the drive's scenario (which must be current), and write the
new pulses to `out` (`write_pulse_table`), on a worker when `spawn[]` is set. The
disabled beams are the ones the shot loses. Progress and the outcome (the capsule
drive's worst-time mode 1 before and after) go to `drive.retune_running`/`retune_status`.
"""
function retune_to_reference!(sess::DesignSession, table::AbstractString, out::AbstractString;
                              by::Symbol = :quad, cap::Real = 1.3, stay::Real = 0.0)
    dr = sess.drive
    dr.retune_running[] && return sess
    if !_drive_current(sess)
        dr.retune_status[] = "compute the drive first: the retune reuses its operator and radiation model"
        return sess
    end
    pulses = read_pulse_table(table)
    s = dr.scenario[]
    beams = sess.beams.beams[]
    enabled = copy(sess.beams.enabled[])
    out = String(out)
    dr.retune_running[] = true
    dr.retune_status[] = "retuning to $(basename(table))…"
    runner = _runner(sess.spawn, sess.deliver)
    runner(() -> retune_pulses(s, beams, enabled, pulses; by, cap, stay), function (res)
        dr.retune_running[] = false
        if res.err !== nothing
            dr.retune_status[] = "retune failed: " * _errmsg(res.err)
            return
        end
        v = res.value
        write_pulse_table(out, v.pulses)
        ε1(h) = 100 * maximum(mode_amplitude(h.ylm[:, k], 1) for k in eachindex(h.time))
        ΔT(h) = maximum(abs.(h.Tr .- v.reference.Tr))
        dr.retune_status[] = @sprintf("%d %s powers: worst ε₁ %.2f%% → %.2f%%, |ΔTᵣ| %.1f → %.1f eV; wrote %s",
                                      length(v.retune.keys), by, ε1(v.before), ε1(v.after),
                                      ΔT(v.before), ΔT(v.after), basename(out))
        sess.status[] = "wrote the retuned pulses to " * basename(out)
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
