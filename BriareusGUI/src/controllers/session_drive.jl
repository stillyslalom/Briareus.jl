# This file covers the drive of a session, with the scenario and its cavities, the evaluation,
# the parts left out of the calculations, and the VisRad workspace, STL and results files.

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
only the laser deposition, and `status` explains why and what to do.
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
                       @sprintf(" · on the whole target: the cavity point is in the open (%.0f%% of its rays escape)",
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

# Whether element `e` of the cavity lies on its boundary, as walls, caps and package faces
# do: the point one element size behind its centroid lies outside the volume. For a
# capsule that point lies inside.
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
                @sprintf("%.2g%% of its rays leave it (an opening without an aperture, or an open end)", 100c.leak))
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
Above the limit [`build_scenario!`](@ref) computes only the laser deposition.
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
part (the drive is then identical to the drive of only the 3ω light).
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
the target is replaced when it is ready. A rejected change sets the `error` of the target.
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
clearance and the x-ray drive then see it). `display_only = :none` follows only the
VisRad flags.

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
parts. Apertures (entrance holes, open windows) are written as separate solids unless
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
