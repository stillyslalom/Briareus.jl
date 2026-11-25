# This file defines the setup of a target, with its state, the build of the `Target` from the
# parts, and the summary.

# ---------------------------------------------------------------- the setup

"""
    TargetSetup(; starter = :hohlraum, resolution = :standard, mesh = nothing,
                parts = PartSpec[])

The target being designed: `parts` (its parts, [`PartSpec`](@ref)s; when none are given,
the parts of `starter`, one of [`STARTERS`](@ref) placed at TCC
([`starter_parts`](@ref)), or with `starter = :none` a sphere at TCC: a target is never
empty), `mesh` (the element sizes, a [`MeshSpec`](@ref)), `resolution` (the preset the
sizes match, one of [`MESH_RESOLUTIONS`](@ref), or `:custom`), `frames` (every named
frame of the design in the chamber, [`part_frames`](@ref): updated with the parts,
before the target is built), `renames` (while an edit updates `frames`, the parts it
renamed, old name → new, so that beams aimed in their frames follow; empty otherwise),
the built `target`, `error` (the last rejected edit), `excluded` (component names left
out of the calculations, see [`set_excluded!`](@ref)), `ghost` (those components, for
drawing, or `nothing`) and `active_foil` (the foil chosen with
[`set_active_foil!`](@ref), or `nothing`; it follows a rename, and
[`active_foil`](@ref) gives the foil the analyses take).

CAD parts (VisRad CAD, STL) are remeshed when their element size, the target's, or
their placement changes, and when an undo or a session brings back a size not meshed
before: seconds to minutes for 10⁵ facets. With `spawn[]` set (a window shares its
session's) that runs on worker threads: the edit is committed at once with those parts
drawn from their file's facets, `pending` names them, and the target is built again
from the current design when they are ready (through `deliver[]`, on the GUI thread).
Without it the edit waits for the meshing.

Building the `Target` itself (its BVH and element tables) takes a second for 10⁵
triangles, so with `spawn[]` a target of more than [`BACKGROUND_TRIANGLES`](@ref)
triangles that has not been built before is built on a worker too: the edit's inputs
(`parts`, `mesh`, …) are committed at once, `target` and `ghost` keep the previous
design's until the new ones are ready, and `building` is set meanwhile. A CAD part
given the role of the capsule has its sphere fitted on a worker as well (`fitting`
names it; the role is committed with the fit). The drive and the checks wait for both
([`target_ready`](@ref)).

Edit with [`insert_starter!`](@ref), [`new_target!`](@ref), [`set_resolution!`](@ref),
[`edit_mesh!`](@ref), the window functions ([`add_window!`](@ref), …) and the part
functions ([`add_part!`](@ref), …).
"""
struct TargetSetup
    mesh::Observable{MeshSpec}
    resolution::Observable{Symbol}
    parts::Observable{Vector{PartSpec}}
    target::Observable{Target}
    error::Observable{String}
    excluded::Observable{Set{Symbol}}
    ghost::Observable{Union{Nothing,Target}}
    pending::Observable{Set{Symbol}}       # CAD parts being remeshed (drawn from their facets)
    spawn::Base.RefValue{Bool}             # remesh on workers (a window's session sets it)
    deliver::Base.RefValue{Any}            # hands a worker's result to the GUI thread
    mesh_job::Base.RefValue{Int}           # the latest remeshing or build job (earlier ones are dropped)
    building::Observable{Bool}             # the target of the committed inputs is being built
    fitting::Observable{Set{Symbol}}       # CAD parts whose capsule sphere is being fitted
    fit_jobs::Dict{Symbol,Int}             # their jobs (a fit whose job is gone is dropped)
    after_build::Vector{Any}               # run once the target being built is committed
    frames::Observable{Dict{Symbol,Frame}} # the named frames of the committed design (`part_frames`)
    renames::Base.RefValue{Dict{Symbol,Symbol}} # old => new part names of the edit committing `frames`
    active_foil::Observable{Union{Nothing,Symbol}} # the chosen foil (`active_foil` resolves it)
end

function TargetSetup(; starter::Symbol = :hohlraum, resolution::Symbol = :standard,
                     mesh::Union{Nothing,MeshSpec} = nothing,
                     parts::AbstractVector{PartSpec} = PartSpec[])
    if mesh === nothing
        resolution in MESH_RESOLUTIONS ||
            throw(ArgumentError("resolution must be one of $MESH_RESOLUTIONS"))
        mesh = MESH_PRESETS[resolution]
    end
    ps = isempty(parts) ? _starter_or_none(starter, mesh) : collect(PartSpec, parts)
    t = build_target(mesh, ps)
    return TargetSetup(Observable(mesh), Observable(_preset_of(mesh)), Observable(ps),
                       Observable(t), Observable(""), Observable(Set{Symbol}()),
                       Observable{Union{Nothing,Target}}(nothing), Observable(Set{Symbol}()),
                       Ref(false), Ref{Any}(f -> f()), Ref(0), Observable(false),
                       Observable(Set{Symbol}()), Dict{Symbol,Int}(), Any[],
                       Observable(part_frames(ps)), Ref(Dict{Symbol,Symbol}()),
                       Observable{Union{Nothing,Symbol}}(nothing))
end

# The parts of a new target, those of a starter, or none for `:none`.
function _starter_or_none(starter::Symbol, mesh::MeshSpec)
    starter === :none && return PartSpec[]
    starter in STARTERS || throw(ArgumentError("a starter is one of $STARTERS, or :none"))
    return starter_parts(starter; capsule_h = mesh.capsule_h)
end

part_frames(ts::TargetSetup) = ts.frames[]

"""
Targets of more triangles than this (≈ 30 ms to build) are built on a worker in a window
(see [`TargetSetup`](@ref)); smaller ones, such as the starters' meshes, at once. Set it
to change the threshold for a session (tests set 0).
"""
const BACKGROUND_TRIANGLES = Ref(30_000)

"""
    target_ready(ts::TargetSetup) -> Bool

Whether `target` is the target of the committed design: no CAD part is being remeshed
(`pending`), no target is being built (`building`) and no capsule sphere is being
fitted (`fitting`). The drive and the checks wait for it in a window.
"""
target_ready(ts::TargetSetup) = isempty(ts.pending[]) && !ts.building[] && isempty(ts.fitting[])

# What the target is waiting for, for status lines ("" when ready).
function target_waiting(ts::TargetSetup)
    out = String[]
    isempty(ts.pending[]) || push!(out, "remeshing " * join(sort!(String.(collect(ts.pending[]))), ", "))
    isempty(ts.fitting[]) ||
        push!(out, "fitting the capsule to " * join(sort!(String.(collect(ts.fitting[]))), ", "))
    ts.building[] && isempty(ts.pending[]) && push!(out, "building the target")
    return join(out, " · ")
end

# Run `f()` once the target of the committed design is in `target` (at once if it is).
function _after_build!(f, ts::TargetSetup)
    ts.building[] ? push!(ts.after_build, f) : f()
    return nothing
end

"""
    await_target!(ts::TargetSetup) -> Bool

Make `target` the committed design's now if a worker is still building it, as the checks
that use it need. The build in flight is waited for instead of being repeated. CAD parts
still remeshing stay `pending`, with their facets drawn, unless the session no longer
spawns jobs (a closed window), in which case they are meshed in this call. Returns
whether the target is current.
"""
function await_target!(ts::TargetSetup)
    (ts.building[] || (!ts.spawn[] && !isempty(ts.pending[]))) || return true
    return _try_build!(ts; sync = true)
end

"""
    build_target(mesh, parts) -> Target
    build_target(ts::TargetSetup) -> Target

The `Target` of `parts` ([`PartSpec`](@ref)s) meshed with the element sizes `mesh` (a
[`MeshSpec`](@ref) or a preset name). Throws `ArgumentError` for parts the geometry
rejects (overlapping windows, a window past the end of a hohlraum's wall, a duplicate
component name, …). No parts at all give a target with no components: the empty design
that File ▸ New target ▸ Empty starts from.
"""
build_target(resolution::Symbol, parts) = build_target(MESH_PRESETS[resolution], parts)

build_target(mesh::MeshSpec, parts::AbstractVector{PartSpec}; raw = ()) =
    _target_of(_target_components(mesh, parts; raw))
build_target(ts::TargetSetup) = build_target(ts.mesh[], ts.parts[])

# The components of a design (`raw`: names of CAD parts drawn from their file's facets).
function _target_components(mesh::MeshSpec, parts::AbstractVector{PartSpec}; raw = ())
    mesh.h > 0 && mesh.capsule_h > 0 || throw(ArgumentError("element sizes must be positive"))
    comps = Component[]
    frames = isempty(parts) ? Dict{Symbol,Frame}() : part_frames(parts)
    for part in parts
        new = part_components(part, mesh.h; frame = frames[part.name], raw = part.name in raw,
                              capsule_h = mesh.capsule_h)
        for c in new
            any(d -> d.name === c.name, comps) && throw(ArgumentError("two components are named $(c.name)"))
        end
        append!(comps, new)
    end
    return comps
end

# The lips' `GapClosure`: closed, they re-emit with the wall's albedo.
function _leh_closure(p::AbstractDict, wall::Material)
    p[:leh_closure_width] >= 0 || throw(ArgumentError("the closure duration must not be negative"))
    0 < p[:leh_open_fraction] < 1 ||
        throw(ArgumentError("the open radius must be between 0 and 1 of the LEH radius"))
    return GapClosure(closure = p[:leh_closure_time], width = p[:leh_closure_width], after = wall.albedo)
end

# Rebuild with changed inputs, committing them only if the geometry accepts them.
# `prebuilt` is the `(target, ghost)` of exactly these inputs when a worker built it
# already ([`import_visrad_workspace!`](@ref)). With `spawn[]`, CAD parts not meshed yet
# are drawn from their file's facets and meshed on workers (`pending`), and the target
# is built again when they are ready; a large target not built before is built on a
# worker (`building`) unless `sync`, the inputs committed at once.
function _try_build!(ts::TargetSetup; mesh = ts.mesh[], parts = ts.parts[], excluded = ts.excluded[],
                     prebuilt = nothing, sync::Bool = !ts.spawn[], renames = Dict{Symbol,Symbol}())
    try
        frames = part_frames(parts)
        waiting = prebuilt === nothing && ts.spawn[] ? _unmeshed_parts(mesh, parts, excluded) : Symbol[]
        # components left out of the calculations (those still in the target)
        built, lists = if prebuilt === nothing
            keep, gh = _split_lists(_target_components(mesh, parts; raw = union(Set{Symbol}(excluded), waiting)),
                                    excluded)
            !sync && _build_later(keep, gh) ? (nothing, (keep, gh)) : (_targets_of(keep, gh), nothing)
        else
            prebuilt, nothing
        end
        _assign!(ts.excluded, Set{Symbol}(excluded))
        built === nothing || (ts.ghost[] === built[2] || (ts.ghost[] = built[2]))
        # This comes first, so that what follows a frame sees the frames of the design,
        # and beams aimed in the frame of a renamed part follow its new name.
        ts.renames[] = Dict{Symbol,Symbol}(renames)
        try
            _assign!(ts.frames, frames)
        finally
            ts.renames[] = Dict{Symbol,Symbol}()
        end
        if mesh != ts.mesh[]
            ts.mesh[] = mesh
            _assign!(ts.resolution, _preset_of(mesh))
        end
        # the chosen foil keeps its role through a rename
        a = ts.active_foil[]
        a !== nothing && haskey(renames, a) && (ts.active_foil[] = renames[a])
        parts === ts.parts[] || (ts.parts[] = parts)
        g = (ts.mesh_job[] += 1)                  # a meshing or build job for earlier inputs is dropped
        ts.error[] = ""
        if built === nothing
            # the target and ghost of the previous inputs stay until these are built
            _assign!(ts.building, true)
            _assign!(ts.pending, Set{Symbol}(waiting))
            deliver = ts.deliver[]
            job = () -> _targets_of(lists...)
            _run_job(() -> _deliver_compiled(deliver, out -> _built!(ts, g, out), _try_job(job)), true)
        else
            built[1] === ts.target[] || (ts.target[] = built[1])
            _assign!(ts.pending, Set{Symbol}(waiting))
            _built_now!(ts)
        end
        isempty(waiting) || _mesh_in_background!(ts, g, waiting, mesh, parts)
        return true
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

# A worker built the target of job `g`: commit it unless a later edit superseded it.
function _built!(ts::TargetSetup, g::Integer, out)
    g == ts.mesh_job[] || return false
    if out.err !== nothing                       # a bug, since the inputs were checked
        ts.error[] = "building the target failed: " * _errmsg(out.err)
        empty!(ts.after_build)
        _assign!(ts.building, false)
        return false
    end
    t, ghost = out.value
    ts.ghost[] === ghost || (ts.ghost[] = ghost)
    t === ts.target[] || (ts.target[] = t)
    _built_now!(ts)
    return true
end

# The target of the committed design is in place: run what waited for it.
function _built_now!(ts::TargetSetup)
    _assign!(ts.building, false)
    isempty(ts.after_build) && return
    fs = copy(ts.after_build)
    empty!(ts.after_build)
    foreach(f -> f(), fs)
    return
end

# Mesh the parts `names` on workers, then build the target again from the design as it
# is then (job `g`: a later edit's job supersedes it; a part both ask for is meshed once).
function _mesh_in_background!(ts::TargetSetup, g::Integer, names, mesh, parts)
    deliver = ts.deliver[]
    todo = [p for p in parts if p.name in names]
    job = () -> _mesh_parts(todo, mesh, Set{Symbol}(), _ -> nothing, ""; all = parts)
    _run_job(() -> _deliver_compiled(deliver, out -> _meshed!(ts, g, out), _try_job(job)), true)
    return nothing
end

function _meshed!(ts::TargetSetup, g::Integer, out)
    g == ts.mesh_job[] || return false
    if out.err !== nothing
        ts.error[] = "remeshing " * join(sort!(String.(collect(ts.pending[]))), ", ") * " failed: " *
                     _errmsg(out.err)
        return false
    end
    return _try_build!(ts)
end

"""
    set_excluded!(ts::TargetSetup, name::Symbol, excluded::Bool) -> Bool

Leave component `name` out of the calculations (`excluded = true`) or put it back.
An excluded component is not part of `target` (no laser, view factors, clearance or
checks see it) but stays listed and drawn as a ghost (`ghost`), so it can be put
back; a hohlraum with an excluded component is no longer treated as closed. `name` may
be a part (every component of it) or an assembly (every part of its subtree,
[`assembly_parts`](@ref)); a part put back comes back whole, including pieces that were
left out separately.
"""
function set_excluded!(ts::TargetSetup, name::Symbol, excluded::Bool)
    ex = first(_with_excluded(ts, ts.excluded[], (name,), excluded))
    ex == ts.excluded[] && return true
    return _try_build!(ts; excluded = ex)
end

# The set of names left out, `ex`, with `names` (components, parts or assemblies) left
# out (`out`) or put back, and the names put in or taken out of it: an assembly stands for
# the parts of its subtree; a part put back takes its components with it.
function _with_excluded(ts::TargetSetup, ex, names, out::Bool)
    ex = Set{Symbol}(ex)
    ps = ts.parts[]
    leaves = Symbol[]
    for n in names
        k = part_index(ts, Symbol(n))
        k !== nothing && ps[k].kind === :assembly ? append!(leaves, assembly_parts(ps, ps[k].name)) :
                                                    push!(leaves, Symbol(n))
    end
    if out
        union!(ex, leaves)
    else
        setdiff!(ex, leaves)
        back = Set(leaves)
        for t in (ts.target[], ts.ghost[])
            t === nothing && continue
            foreach(c -> c.part in back && delete!(ex, c.name), t.components)
        end
    end
    return ex, leaves
end

"Every component name of the design: in the calculations, then excluded ones."
all_component_names(ts::TargetSetup) =
    [component_names(ts.target[]); ts.ghost[] === nothing ? Symbol[] : component_names(ts.ghost[])]

"""
    set_leh_closure!(ts::TargetSetup, on::Bool; part = nothing) -> Bool

Close hohlraum part `part`'s entrance holes (by default the design's first hohlraum
part's) during the pulse (`on`), with the settings of [`LEH_CLOSURE_FIELDS`](@ref), or
keep them open (its value `leh_closure`). Each hole is then split into the hole and its
closing lip (`<part>_leh_top_lip`, …), shown, selected and left out with the hole.
"""
function set_leh_closure!(ts::TargetSetup, on::Bool; part = nothing)
    try
        i = _window_part(ts, part)
        p = ts.parts[][i]
        (p.values[:leh_closure] != 0) == on && return true
        return edit_part!(ts, i, :leh_closure, on ? "1" : "0")
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

"""
    new_target!(ts::TargetSetup, starter = :hohlraum) -> Bool

Start a new target (File ▸ New target): only the parts of `starter` (one of
[`STARTERS`](@ref), placed at TCC, [`starter_parts`](@ref)), or with `starter = :none`
no parts at all (an empty design, to build from Insert); nothing is left out of the
calculations, and the element sizes are kept.
"""
function new_target!(ts::TargetSetup, starter::Symbol = :hohlraum)
    parts = try
        _starter_or_none(starter, ts.mesh[])
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    _cancel_fits!(ts)
    _try_build!(ts; parts, excluded = Set{Symbol}()) || return false
    ts.active_foil[] = nothing
    return true
end

"""
    set_resolution!(ts::TargetSetup, resolution::Symbol) -> Bool

Mesh element sizes from a preset, one of [`MESH_RESOLUTIONS`](@ref) (see
[`MESH_PRESETS`](@ref)).
"""
function set_resolution!(ts::TargetSetup, resolution::Symbol)
    resolution in MESH_RESOLUTIONS ||
        throw(ArgumentError("resolution must be one of $MESH_RESOLUTIONS"))
    ts.mesh[] == MESH_PRESETS[resolution] && return true
    return _try_build!(ts; mesh = MESH_PRESETS[resolution])
end

"""
    edit_mesh!(ts::TargetSetup, key::Symbol, text) -> Bool

Set an element size (`:h` or `:capsule_h`, mm) from text; the resolution becomes the
preset the sizes match, or `:custom`.
"""
function edit_mesh!(ts::TargetSetup, key::Symbol, text)
    key in (:h, :capsule_h) || throw(ArgumentError("unknown mesh field :$key"))
    try
        v = _parse_number(text, key === :h ? "the element size" : "the capsule element size";
                          positive = true)
        m = ts.mesh[]
        m = key === :h ? MeshSpec(v, m.capsule_h) : MeshSpec(m.h, v)
        return _try_build!(ts; mesh = m)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end

# ---------------------------------------------------------------- summary

"""
    target_summary(ts::TargetSetup) -> String

The hohlraum parts and their windows, the other parts, and the element count.
"""
function target_summary(ts::TargetSetup)
    t = ts.target[]
    ps = ts.parts[]
    nh = count(p -> p.kind === :hohlraum, ps)
    nw = sum(p -> length(p.windows), ps; init = 0)
    no = count(p -> !(p.kind in (:hohlraum, :frame, :assembly)), ps)
    kind = nh == 0 ? "target" : nh == 1 ? "hohlraum" : "$nh hohlraums"
    nw > 0 && (kind *= nw == 1 ? " with 1 window" : " with $nw windows")
    no > 0 && (kind *= nh == 0 ? (no == 1 ? " of 1 part" : " of $no parts") : no == 1 ? " + 1 part" : " + $no parts")
    w = target_waiting(ts)
    return "$kind · $(nelements(t)) elements" * (isempty(w) ? "" : " · " * w)
end

"""
Names of a target's components, in order: its parts (the drive's fitted pieces of a CAD
part are the part), then the pieces of the parts made in pieces (a hohlraum part's
`<part>_wall`, …), which are also addressed separately.
"""
component_names(t::Target) = unique!([Briareus.parts(t); [c.name for c in t.components if _is_piece(c)]])

"""
Names of a target's opaque (non-aperture) components, as the checks list them: each part
(the drive's fitted pieces of a CAD part are the part), but each piece of a part made in
pieces (a hohlraum part's wall and caps).
"""
surface_names(t::Target) = unique!([_is_piece(c) ? c.name : c.part for c in t.components if c.kind === :surface])

# Whether `c` is a piece of a part made in pieces (named `<part>_<piece>`; a hohlraum
# part's), as opposed to a whole part or a fitted piece of a CAD part (`<part>/cyl1`, …).
_is_piece(c::Component) = c.part !== c.name && startswith(String(c.name), String(c.part) * "_")

# The component of `names` that a name in a diagnostic definition (written for a hohlraum's
# pieces: `leh_top`) refers to: itself when the target has it, else the first hohlraum
# part's piece of that name (`hohlraum_leh_top`), else itself.
function _piece_named(names, n::Symbol)
    n in names && return n
    k = findfirst(m -> endswith(String(m), "_" * String(n)), names)
    return k === nothing ? n : names[k]
end

# A diagnostic's view (`Briareus.diagnostic`) on a target with components `names`: the
# components it refers to (`through`, `center`, `area`) resolved with `_piece_named`, with
# a separate cache; the view itself when it refers to none the target lacks.
function _view_on(v, names)
    T = typeof(v)
    refs = (:through, :center, :area)
    fs = fieldnames(T)
    any(f -> f in refs && getfield(v, f) isa Symbol && !(getfield(v, f) in names), fs) || return v
    vals = map(fs) do f
        x = getfield(v, f)
        f in refs && x isa Symbol ? _piece_named(names, x) :
        x isa Union{IdDict,WeakKeyDict} ? typeof(x)() : x isa ReentrantLock ? ReentrantLock() : x
    end
    return T(vals...)
end
