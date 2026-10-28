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

const TASK_LABELS = Dict(:drive => "Drive", :optimize => "Optimization", :clearance => "Clearance",
                         :talis => "TaLIS")

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
`talis` ([`TaLISState`](@ref)), `time` ([`TimeState`](@ref): pulses, the drive history
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
pointed for the target. Changing the target kind swaps a NIF preset for the OMEGA
backlighter set (and back), unless the beams were imported.

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
    visrad_source::Base.RefValue{String}   # the imported VisRad workspace (export template)
    meshing::Observable{String}            # a workspace import or CAD meshing in progress ("" when none)
    mesh_generation::Base.RefValue{Int}    # import and meshing requests (latest wins)
    build_queued::Base.RefValue{Bool}      # a drive asked for while CAD parts remesh
    cad_sources::Dict{String,String}       # a saved CAD copy => the file it was saved from
    cad_changed::Observable{Vector{String}} # saved copies whose source changed since the save
end

function DesignSession(; kind::Symbol = :hohlraum, resolution::Symbol = :standard,
                       preset::Symbol = kind === :foil ? :omega_backlighter : :nif_quads,
                       beams = nothing, deliver = f -> f())
    ts = TargetSetup(; kind, resolution)
    bt = beams === nothing ?
         BeamTable(preset_beams(preset; target_kind = kind, params = ts.params[]); preset) :
         BeamTable(beams)
    # the target remeshes CAD parts on the session's workers (the same switches)
    ts.deliver[] = deliver
    sess = DesignSession(ts, bt, DriveState(), OptimizeState(), ReviewState(), TaLISState(),
                         TimeState(), ViewState(), Observable(""), ts.deliver, ts.spawn, EditHistory(),
                         Ref(""), Observable(""), Ref(0), Ref(false), Dict{String,String}(),
                         Observable(String[]))
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
    onany((_...) -> _drive_inputs_changed!(sess), dr.cavity, dr.cavity_point, dr.films,
          dr.beam_model, dr.profile)
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
    await_target!(sess.target)
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
        _optimum_stale(sess.beams.beams[], r) && return (:attention, "out of date: " * txt)
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
# (a new preset or target template) or moved to other ports, or, once applied, the
# optimized beams were edited since.
function _optimum_stale(bs::AbstractVector, r)
    all(i -> i <= length(bs), r.indices) || return true
    r.applied && return any(k -> bs[r.indices[k]] != r.beams[k], eachindex(r.indices))
    return any(k -> bs[r.indices[k]].port != r.beams[k].port, eachindex(r.indices))
end

# ---------------------------------------------------------------- drive

"""
    drive_stale(sess::DesignSession) -> Bool

Whether the scenario was built for another target, other beam ports, cavity or beam
model (it must be computed again), or the target of the design is still being built. `false` before the
first build.
"""
function drive_stale(sess::DesignSession)
    k = sess.drive.key[]
    k === nothing && return false
    return sess.target.building[] || k[1] !== sess.target.target[] || k[2] != _ports(sess.beams.beams[]) ||
           k[3] != _cavity_key(sess.drive) || k[4] != _profile_state(_deposition_profile(sess.drive))
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
    isempty(beams) && (dr.status[] = "add beams first"; return sess)
    ts = sess.target
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
    closed = _closed_enclosure(sess.target)
    kind, params, mesh, parts = ts.kind[], ts.params[], ts.mesh[], ts.parts[]
    runner(function ()
               # the CAD parts' file facets, for the drive's fit of their curved faces
               references = _fit_references(t, kind, params, mesh, parts)
               _build_scenario(t, beams; closed, cavity, profile, references)
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
        dr.key[] = (t, _ports(beams), cavity, _profile_state(profile))
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
                      (s.model === nothing ? _radiation_note(s.target) : _facet_note(s.target))
        _request_evaluation!(sess)
    end)
    return sess
end

# Whether the drive's target is a closed hohlraum enclosure, whose view factors are
# balanced so its rows sum to one: the template's components only (free parts lie
# outside it), none of them left out (that leaves a hole).
const _TEMPLATE_HOHLRAUM = (:wall, :cap_top, :cap_bottom, :leh_top, :leh_bottom, :leh_top_lip,
                            :leh_bottom_lip, :capsule)
_template_component(n::Symbol) = n in _TEMPLATE_HOHLRAUM || occursin(r"^window_\d+$", String(n))
function _closed_enclosure(ts::TargetSetup)
    ts.kind[] === :hohlraum || return false
    all(c -> _template_component(c.name), ts.target[].components) || return false
    g = ts.ghost[]
    return g === nothing || !any(c -> _template_component(c.name), g.components)
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
`feasible`): about 160,000 elements on 32 GB. The mesh presets of the templates stay far
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
    include_components!(sess::DesignSession, names, include::Bool; spawn = sess.spawn[]) -> sess

Put the components `names` back into the calculations (`include = true`) or leave them
out (see [`set_excluded!`](@ref)). A CAD part left out is drawn from its file's facets;
putting it back repairs its windings and remeshes it, which takes tens of seconds for
10⁵ facets: with `spawn` that runs on worker threads, with progress in `meshing`, and
the target is replaced when it is ready. A rejected change sets the target's `error`.
"""
function include_components!(sess::DesignSession, names, include::Bool; spawn::Bool = sess.spawn[])
    ts = sess.target
    ex = copy(ts.excluded[])
    include ? setdiff!(ex, names) : union!(ex, names)
    ex == ts.excluded[] && return sess
    if !spawn || !include
        _try_build!(ts; excluded = ex)
        return sess
    end
    kind, params, mesh, windows, parts = ts.kind[], ts.params[], ts.mesh[], ts.windows[], ts.parts[]
    g = (sess.mesh_generation[] += 1)
    deliver = sess.deliver[]
    label = join(string.(names), ", ")
    sess.meshing[] = "meshing $label for the calculations…"
    report(msg) = deliver(() -> g == sess.mesh_generation[] && _assign!(sess.meshing, msg))
    job = function ()
        back = [p for p in parts if p.name in names]
        _mesh_parts(back, kind, params, mesh, ex, report, label; all = parts)
        return _split_components(_target_components(kind, params, mesh, windows, parts; raw = ex), ex)
    end
    commit = function (out)
        g == sess.mesh_generation[] || return
        _assign!(sess.meshing, "")
        if out.err !== nothing
            ts.error[] = _errmsg(out.err)
            sess.status[] = "could not put $label back: " * ts.error[]
            return
        end
        same = ts.kind[] === kind && ts.params[] === params && ts.mesh[] == mesh &&
               ts.windows[] === windows && ts.parts[] === parts
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
size), and the workspace's beams, in their clusters as groups, switched on where
VisRad delivers power. Objects for which `include` is false (by default those VisRad
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
    params, mesh = ts.params[], ts.mesh[]
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
        parts = [PartSpec(:visrad; name = Symbol(o.name), file = path, object = o.name,
                          values = Dict(:albedo => o.albedo, :conversion => o.conversion))
                 for o in supported]
        # the components of objects VisRad leaves out are left out of the calculations,
        # and so (by default) is CAD hardware far from where the beams point
        names = Set{Symbol}(n for o in supported if !include(o) for n in visrad_component_names(o))
        far = display_only === :far ? _far_objects(filter(include, supported), ws.beams) : String[]
        length(far) < count(include, supported) || (far = String[])
        union!(names, Symbol.(far))
        _mesh_parts(parts, :custom, params, mesh, names, report, file)
        report("$file: building the target…")
        built = _split_components(_target_components(:custom, params, mesh, WindowSpec[], parts; raw = names),
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
        same = ts.params[] === params && ts.mesh[] == mesh
        if !_try_build!(ts; kind = :custom, windows = WindowSpec[], parts, excluded = names,
                        prebuilt = same ? built : nothing)
            sess.status[] = "the workspace's target was rejected: $(ts.error[])"
            spawn || throw(ArgumentError(sess.status[]))
            return
        end
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
        set_cut!(sess.view; components = [n for n in surface_names(ts.target[]) if n !== :capsule])
        sess.status[] = "imported $file: $(length(parts)) objects, $(length(bs)) beams" *
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
produced them (`session.briar.toml`, see [`save_session`](@ref)). The drive must be
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
    push!(files, save_session(joinpath(dir, "session.briar.toml"), sess))
    sess.status[] = "exported the results to $(basename(normpath(dir)))"
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
        # the time panel plots the retune (and fires the reference pulses when it had none)
        bt = sess.beams
        nodes = sort(unique(reduce(vcat, [p.t for p in values(pulses)])))
        sess.time.retune[] = (; label = basename(table), by, pulses = v.pulses, reference_pulses = pulses,
                              reference_by = _pulse_key(beams, pulses), nodes, beams, enabled,
                              table = (; cluster = copy(bt.cluster[]), colors = cluster_colors(bt)),
                              before = v.before, after = v.after, reference = v.reference)
        sess.time.pulses[] === nothing && length(bt.beams[]) == length(beams) &&
            (sess.time.pulses[] = (; kind = :table, table = pulses, by = _pulse_key(beams, pulses),
                                   shape = nothing, nodes, times = nodes, label = basename(table)))
        set_chart!(sess.time, :retune)
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
    s = sess.drive.scenario[]
    # objectives the drive cannot serve: say why and what to do, instead of a failed job
    if st.objective[] === :irradiance
        FOIL_NAME in component_names(s.target) ||
            (st.status[] = "the foil objective needs a foil: a foil target's, or a part with the role " *
                           "Backlighter foil (Properties ▸ Role)"; return sess)
    elseif s.model === nothing
        st.status[] = "capsule symmetry needs the x-ray drive, and this target's drive is laser only" *
                      _radiation_note(s.target)
        return sess
    elseif isempty(s.observed)
        st.status[] = "capsule symmetry needs a hohlraum with a capsule, or a sphere part with the " *
                      "role Capsule (Properties ▸ Role)"
        return sess
    end
    bt = sess.beams
    idx = findall(bt.enabled[])
    start_optimize!(st, s, bt.beams[][idx], idx;
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
    (load_preset!(sess.beams, preset; target_kind = sess.target.kind[],
                  params = sess.target.params[]); sess)

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
    for (job, cancel!) in ((sess.optimize, cancel_optimize!), (sess.review, cancel_clearance!),
                           (sess.talis, cancel_talis!))
        job.running[] || continue
        cancel!(job)                       # the worker stops at its next check
        job.running[] = false
        job.status[] = "abandoned when the window closed"
    end
    tm = sess.time
    tm.queued[] = false
    tm.running[] && (tm.generation[] += 1; tm.cancel[] = true; tm.running[] = false;
                     tm.status[] = "abandoned when the window closed")
    return sess
end
