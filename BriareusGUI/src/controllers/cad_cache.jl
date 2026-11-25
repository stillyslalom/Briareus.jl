# This file meshes the parts and keeps the results in caches. It covers the components of each
# part, the CAD read from VisRad workspaces and STL files and remeshed, the parts meshed on
# workers, and the targets built from the components.

# A cache that keeps its `kept` most recently used entries, as the cache of built targets
# did with a dictionary and a separate list of keys. Looking an entry up (`get`, `getindex`)
# makes it the most recent; an entry stored beyond `kept` drops the least recent. Every
# access holds `_STL_LOCK` (`_cached`, `_iscached`, `clear_caches!`).
struct _LRU{K,V}
    data::Dict{K,V}
    order::Vector{K}                 # least recently used first
    kept::Int
end
_LRU{K,V}(kept::Integer) where {K,V} = _LRU{K,V}(Dict{K,V}(), K[], Int(kept))

Base.haskey(c::_LRU, key) = haskey(c.data, key)
Base.length(c::_LRU) = length(c.data)
Base.isempty(c::_LRU) = isempty(c.data)
Base.empty!(c::_LRU) = (empty!(c.data); empty!(c.order); c)
Base.delete!(c::_LRU, key) = (delete!(c.data, key); filter!(!=(key), c.order); c)

function _touch!(c::_LRU, key)
    isempty(c.order) || c.order[end] == key || (filter!(!=(key), c.order); push!(c.order, key))
    return c
end

function Base.get(c::_LRU, key, default)
    haskey(c.data, key) || return default
    _touch!(c, key)
    return c.data[key]
end

function Base.getindex(c::_LRU, key)
    v = c.data[key]
    _touch!(c, key)
    return v
end

function Base.setindex!(c::_LRU, v, key)
    haskey(c.data, key) ? _touch!(c, key) : push!(c.order, key)
    c.data[key] = v
    while length(c.order) > c.kept
        delete!(c.data, popfirst!(c.order))
    end
    return c
end

# The entries kept per CAD part (`_VRW_PARTS`, `_STL_CACHE`: a design of many CAD parts
# holds two per part, the remeshed elements and the file facets drawn while they mesh) and
# per CAD file (`_VRW_CACHE`, `_STL_SOLIDS`). An entry dropped is read or meshed again when
# a design asks for it. These are engineering defaults.
const _CAD_PARTS_KEPT = 512
const _CAD_FILES_KEPT = 16

const _PART_MEMO = Dict{Any,Vector{Component}}()

# Whether a part is meshed as the capsule (`capsule_h`, geodesic).
_capsule_meshed(p::PartSpec) = p.kind === :sphere && p.values[:capsule_mesh] != 0

function _part_key(p::PartSpec, h::Real, frame, raw::Bool, hc::Real = 0.0)
    stamp = isfile(p.file) ? (mtime(p.file), filesize(p.file)) : nothing
    return (p.kind, p.name, sort!(collect(p.values)), p.two_sided, p.inward, p.aperture, p.file,
            p.object, stamp, Float64(h), frame.R, frame.origin, raw, Float64(hc), p.windows, p.sections)
end

function _part_components(p::PartSpec, h::Real, frame::Frame, raw::Bool, hc::Real = 0.0)
    _no_geometry(p) && return Component[]                 # a user frame, an assembly: no geometry
    R, center = frame.R, frame.origin
    v = p.values
    hp = v[:h] > 0 ? v[:h] : Float64(h)
    hp > 0 || throw(ArgumentError("$(p.name): the element size must be positive"))
    for (k, label, _, _) in PART_SIZE_FIELDS[p.kind]
        v[k] > 0 || throw(ArgumentError("$(p.name): the $(lowercase(label)) must be positive"))
    end
    0 <= v[:albedo] <= 1 || throw(ArgumentError("$(p.name): the albedo must be between 0 and 1"))
    0 <= v[:conversion] <= 1 ||
        throw(ArgumentError("$(p.name): the conversion must be between 0 and 1"))
    p.kind === :hohlraum && return _hohlraum_components(p, hp, frame)
    out = !p.inward
    k = p.kind
    if k === :sphere && hc > 0
        r = v[:r1]
        m = icosphere_mesh(r; n = icosphere_frequency(r, v[:h] > 0 ? v[:h] : Float64(hc)), outward = out)
        surf = SphereSurface(r)
    elseif k === :sphere
        r = v[:r1]
        (; capsule_ntheta, capsule_nphi) = capsule_counts(r, hp)
        m = sphere_mesh(r; ntheta = capsule_ntheta, nphi = capsule_nphi, outward = out)
        surf = SphereSurface(r)
    elseif k === :disk
        r = v[:r1]
        m = disk_mesh(r, 0.0; nr = _count(r, hp), nphi = _nphi(r, hp), up = out)
        surf = DiskSurface()
    elseif k === :annulus
        ri, ro = v[:r1], v[:r2]
        ro > ri || throw(ArgumentError("$(p.name): the outer radius must exceed the inner"))
        m = annulus_mesh(ri, ro, 0.0; nr = _count(ro - ri, hp), nphi = _nphi(ro, hp),
                         up = out)
        surf = DiskSurface()
    elseif k === :rectangle
        w, ht = v[:r1], v[:r2]
        m = rectangle_mesh(w, ht; nu = _count(w, hp), nv = _count(ht, hp))
        out || (m = transform(m, [1 0 0; 0 -1 0; 0 0 -1]))          # face −z
        surf = nothing
    elseif k === :visrad
        cs = _visrad_components(p, hp, transpose(R) * (zero(Vec3) - center); raw)
        return [transform(c, R, center) for c in cs]
    elseif k === :mesh
        v[:scale] > 0 || throw(ArgumentError("$(p.name): the scale must be positive"))
        # facing the target center, remeshed into elements of the part's size
        toward = transpose(R) * (zero(Vec3) - center)
        m = _stl_mesh(p, hp, toward; raw)
        # a part that is one cylinder or sphere (a tube, a shell) gets it for the view factors
        surf = raw ? nothing : _whole_surface(m, _stl_mesh(p, hp, toward; raw = true))
    elseif k === :cylinder
        r, L = v[:r1], v[:length]
        m = cylinder_mesh(r, -L / 2, L / 2; nz = _count(L, hp), nphi = _nphi(r, hp),
                          inward = !out)
        surf = CylinderSurface(r)
    else
        r1, r2, L = v[:r1], v[:r2], v[:length]
        m = frustum_mesh(r1, -L / 2, r2, L / 2; nz = _count(L, hp),
                         nphi = _nphi(max(r1, r2), hp), inward = !out)
        surf = nothing
    end
    c = p.aperture ? aperture(p.name, m; surface = surf) :
        Component(p.name, m, Material(p.name; albedo = v[:albedo], conversion = v[:conversion]);
                  two_sided = p.two_sided, surface = surf)
    return [transform(c, R, center)]
end

# The components of a `:hohlraum` part: core `hohlraum` (gold wall, no capsule) with
# elements of `h` and the part's ends (core `top_end`, `bottom_end`: a closed end's
# entrance hole and lip merged into its cap, of the wall's material; an open end with
# neither), each component renamed `<part>_<name>` with `part = <part>`. Placed by
# `frame`; at the chamber frame the meshes are not moved.
function _hohlraum_components(p::PartSpec, h::Real, frame::Frame)
    v = p.values
    wall = gold(albedo = v[:albedo], conversion = v[:conversion])
    closing = v[:leh_closure] != 0
    lip = closing ? _leh_closure(v, wall) : nothing
    fopen = closing ? v[:leh_open_fraction] : 0.84
    t = try
        rb = v[:bottom_leh_radius]
        rb >= 0 || throw(ArgumentError("the bottom LEH radius must be positive (0: the LEH radius)"))
        ends = map((:top_end, :bottom_end)) do e
            k = _end_kind(v, e)
            k === :leh ? (e === :bottom_end && rb > 0 && rb != v[:leh_radius] ? rb : nothing) : k
        end
        prof = hohlraum_profile(p)
        dims = prof === nothing ? (; length = v[:length], radius = v[:radius]) : (; profile = prof)
        hohlraum(; dims..., leh_radius = v[:leh_radius], capsule_radius = 0.0, wall,
                 windows = [_wall_window(w, k) for (k, w) in enumerate(p.windows)], h,
                 leh_closure = lip, leh_open_radius = fopen * v[:leh_radius],
                 top_end = ends[1], bottom_end = ends[2])
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("$(p.name): " * err.msg))
    end
    out = Component[Component(piece_name(p.name, c.name), c.mesh, c.material; kind = c.kind,
                              two_sided = c.two_sided, surface = c.surface, part = p.name)
                    for c in t.components]
    frame == CHAMBER && return out
    return Component[transform(c, frame.R, frame.origin) for c in out]
end

# A VisRad object's components in the chamber frame (windings repaired toward `toward`,
# CAD surfaces remeshed into elements of `h`; with `raw` the file's facets as they are;
# primitives on VisRad's grid unless the part sets a separate element size), kept per
# workspace version and settings: every edit rebuilds the target, and remeshing a CAD
# part of 10⁵ facets takes tens of seconds.
# The key of `_visrad_key`: the file and its modification time, the object, the element
# size, the point open sheets face, two-sidedness, albedo, conversion, the part and `raw`.
const _VrwPartKey = Tuple{String,Float64,String,Float64,Union{Nothing,NTuple{3,Float64}},Bool,Float64,Float64,
                          Symbol,Bool}
const _VRW_PARTS = _LRU{_VrwPartKey,Vector{Component}}(_CAD_PARTS_KEPT)

# `cache[key]`, computed by `f()` outside the lock when missing (parts mesh in parallel,
# and a long meshing never holds up the GUI thread's lookups of other parts). A key
# being computed is computed once: other callers wait for it (an edit made while a
# part remeshes in the background asks for the same part again).
const _INFLIGHT = Dict{Any,Base.Event}()

function _cached(f, cache, key)
    # (no cache holds `nothing`; the entry is taken under the lock that found it, since
    # a bounded cache may drop it afterwards)
    found, ev, mine = lock(_STL_LOCK) do
        v = get(cache, key, nothing)
        v === nothing || return (v, nothing, false)
        e = get(_INFLIGHT, key, nothing)
        e === nothing || return (nothing, e, false)
        e = _INFLIGHT[key] = Base.Event()
        return (nothing, e, true)
    end
    found === nothing || return found
    if !mine
        wait(ev)
        v = lock(() -> get(cache, key, nothing), _STL_LOCK)
        return v === nothing ? _cached(f, cache, key) : v      # the other task failed, so build it in this one
    end
    try
        v = f()
        lock(() -> (cache[key] = v), _STL_LOCK)
        return v
    finally
        lock(() -> delete!(_INFLIGHT, key), _STL_LOCK)
        notify(ev)
    end
end

_iscached(cache, key) = lock(() -> haskey(cache, key), _STL_LOCK)

"""
    clear_caches!()

Empty the caches the target builds keep between edits: the components of each part, the
VisRad workspaces and STL files read, the CAD parts remeshed from them, the built
targets, and the curved pieces fitted to CAD parts. They release their memory (a remeshed
NIF CAD workspace holds 10⁵ facets and more), and the next build reads and meshes its CAD
again; a test calls this for the state of a new process. A part being meshed on a worker
is stored when it is done.
"""
function clear_caches!()
    lock(_STL_LOCK) do
        foreach(empty!, (_PART_MEMO, _VRW_PARTS, _VRW_CACHE, _STL_CACHE, _STL_SOLIDS, _TARGETS, _FIT_MEMO))
    end
    return nothing
end

# The cache key of a VisRad part's components, the element size a primitive takes
# (`nothing`: VisRad's grid, which reproduces VisRad's element table) and whether the
# part is drawn from its file's facets (`raw`, CAD only).
function _visrad_key(p::PartSpec, o::VisRadObject, h, toward, raw::Bool)
    v = p.values
    cad = o.kind === :mesh
    ph = cad || v[:h] <= 0 ? nothing : v[:h]
    raw &= cad
    key = (abspath(p.file), mtime(p.file), p.object, raw ? 0.0 : cad ? Float64(h) : something(ph, 0.0),
           raw || !cad ? nothing : round.(Tuple(toward); digits = 9), p.two_sided, v[:albedo],
           v[:conversion], p.name, raw)
    return key, ph, raw
end

function _visrad_components(p::PartSpec, h, toward; raw::Bool = false)
    v = p.values
    o = _visrad_object(p.file, p.object)
    cad = o.kind === :mesh
    key, ph, raw = _visrad_key(p, o, h, toward, raw)
    _cached(_VRW_PARTS, key) do
        mat = Material(p.name; albedo = v[:albedo], conversion = v[:conversion])
        # one component per CAD part (`fit = false`): the window shows and keys the parts;
        # the drive fits their curved faces in a separate target (`_fitted_target`)
        cs = raw ? visrad_component(o; material = mat, two_sided = p.two_sided) :
             visrad_component(o; material = mat, two_sided = p.two_sided, orient = !p.two_sided,
                              toward, remesh = cad ? h : nothing, primitive_h = ph, fit = false)
        isempty(cs) && throw(ArgumentError("$(p.name): VisRad object type $(o.type) is not supported"))
        # a CAD part that is one cylinder or sphere (a tube, a shell) gets it for the view
        # factors, fitted on its file's facets
        if cad && !raw
            ref = visrad_component(o; material = mat, two_sided = p.two_sided)
            cs = [length(ref) == 1 ? _with_surface(c, _whole_surface(c.mesh, only(ref).mesh)) : c for c in cs]
        end
        # Components carry the name of the part, renamed or given a role, instead of the
        # name of the object, and the part as their `part`. The pieces of a primitive,
        # `<part>_wall`, …, are named as those of a hohlraum part are, so they follow the
        # part through a rename and a left-out part.
        Component[String(p.name) == o.name && c.part === p.name ? c :
                  Component(String(p.name) * String(c.name)[ncodeunits(o.name)+1:end], c.mesh, c.material;
                            kind = c.kind, two_sided = c.two_sided, surface = c.surface, part = p.name)
                  for c in cs]
    end
end

"""
    part_size_unit(p::PartSpec) -> String

The unit text of a part's element-size field, with what 0 means for it (the target's
size, or for a VisRad primitive the VisRad grid).
"""
function part_size_unit(p::PartSpec)
    p.kind === :visrad || return "mm (0: target's)"
    cad = try
        _visrad_object(p.file, p.object).kind === :mesh
    catch err
        err isa ArgumentError || rethrow()
        true
    end
    return cad ? "mm (0: target's)" : "mm (0: VisRad's grid)"
end

"Names of the components a VisRad object makes (without meshing its CAD)."
visrad_component_names(o::VisRadObject) =
    o.kind === :mesh ? [Symbol(o.name)] : [c.name for c in visrad_component(o)]

# A VisRad workspace's object, the workspace read once per file version.
const _FileKey = Tuple{String,Float64,Int64}          # path, modification time, size
const _VRW_CACHE = _LRU{_FileKey,VisRadWorkspace}(_CAD_FILES_KEPT)

function _visrad_workspace(file::AbstractString)
    isfile(file) || throw(ArgumentError("VisRad workspace not found: $file"))
    key = (abspath(file), mtime(file), filesize(file))
    return _cached(() -> read_visrad_workspace(file), _VRW_CACHE, key)
end

function _visrad_object(file::AbstractString, object::AbstractString)
    ws = _visrad_workspace(file)
    i = findfirst(o -> o.name == object, ws.objects)
    i === nothing && throw(ArgumentError("no object $object in $(basename(file))"))
    return ws.objects[i]
end

# An STL part's mesh: its file's facets (a whole file, or solid `object` of a
# multi-solid file) scaled to mm, oriented (`orient_mesh`: closed solids face out of
# their material, open sheets toward `toward`, flipped for `inward`) and remeshed into
# elements of `h`; with `raw` the facets as they are (drawing only). Every edit rebuilds
# the target, so the result is kept per file version and settings.
# (the key of `_stl_key`: the file's path, modification time and size, the solid, the
# scale, the element size, the point open sheets face, and whether the mesh faces inward)
const _StlKey = Tuple{String,Float64,Int64,String,Float64,Float64,Union{Nothing,NTuple{3,Float64}},Bool}
const _STL_CACHE = _LRU{_StlKey,Mesh}(_CAD_PARTS_KEPT)
const _STL_LOCK = ReentrantLock()

function _stl_key(p::PartSpec, h, toward, raw::Bool)
    isfile(p.file) || throw(ArgumentError("STL file not found: $(p.file)"))
    return (abspath(p.file), mtime(p.file), filesize(p.file), p.object, p.values[:scale],
            raw ? 0.0 : Float64(h), raw ? nothing : round.(Tuple(toward); digits = 9),
            !raw && p.inward)
end

function _stl_mesh(p::PartSpec, h, toward; raw::Bool = false)
    scale = p.values[:scale]
    _cached(_STL_CACHE, _stl_key(p, h, toward, raw)) do
        m0 = if isempty(p.object)
            read_stl(p.file; scale)
        else
            solids = _stl_solids(p.file)
            k = parse(Int, p.object)
            1 <= k <= length(solids) || throw(ArgumentError("$(basename(p.file)) has no solid $k"))
            m = last(solids[k])
            Mesh([scale * v for v in m.vertices], m.triangles, m.patches)
        end
        raw && return m0
        m = remesh_mesh(orient_mesh(m0; toward); size = h)
        p.inward ? Mesh(m.vertices, [(a, c, b) for (a, b, c) in m.triangles], m.patches) : m
    end
end

"""
    _whole_surface(mesh, reference) -> Union{Nothing,CylinderSurface,SphereSurface}

The cylinder or sphere a remeshed CAD part lies on when the whole part is one
(`fit_surfaces` on its file's facets, `reference`, whose vertices lie on the true surface;
every element of `mesh` in that one region), else `nothing`. The part then keeps its name
and elements and is integrated on the true surface: a tube or shell tessellated in slivers
otherwise carries hundreds of facets per element, and its view factors cost ~20× more on the
flat path (research notes, "Finely faceted curved CAD on the flat path"). A part of several
faces, such as a hohlraum half with its bore, walls and ends, stays one flat component in
the window, which shows and keys parts; the drive splits it into its fitted pieces in a
separate target (`_fitted_target`).
"""
function _whole_surface(m::Mesh, ref::Mesh)
    npatches(m) >= 6 && length(ref.triangles) >= 6 || return nothing
    # `fit_surfaces` grows every face of a part (minutes for a 4 × 10⁴-facet CAD capsule
    # with its fill tube): only a part whose vertices are all on one sphere or cylinder
    # is worth it
    _on_one_surface(ref) || return nothing
    regs = try
        fit_surfaces(m; reference = ref, kinds = (:cylinder, :sphere))
    catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
    length(regs) == 1 && length(regs[1].patches) == npatches(m) || return nothing
    s = regs[1].surface
    # a nearly flat part fits a sphere or cylinder of any large radius; visibility ignores
    # a curved part's facets between points inside it, so only a true tube or shell counts
    lo, hi = extrema(v -> v[1], ref.vertices), extrema(v -> v[2], ref.vertices)
    zl = extrema(v -> v[3], ref.vertices)
    extent = hypot(lo[2] - lo[1], hi[2] - hi[1], zl[2] - zl[1])
    return s !== nothing && s.radius <= extent ? s : nothing
end

# Whether every vertex of `m` lies on one sphere or one cylinder, within twice
# `fit_surfaces`' default tolerance (max(0.1 µm, 10⁻³ R)): algebraic least-squares fits to
# all vertices, the cylinder's axis the direction the facet normals are all
# perpendicular to. Milliseconds for 10⁵ vertices.
function _on_one_surface(m::Mesh)
    P = unique(m.vertices)
    length(P) >= 6 || return false
    tol(r) = 2 * max(1e-4, 1e-3 * r)
    # sphere: |p|² = 2 c·p + (r² − |c|²)
    A = [2p[k] for p in P, k in 1:3]
    x = [A ones(length(P))] \ [dot(p, p) for p in P]
    c = Vec3(x[1], x[2], x[3])
    r = sqrt(max(x[4] + dot(c, c), 0.0))
    r > 0 && maximum(p -> abs(norm(p - c) - r), P) <= tol(r) && return true
    # For a cylinder, the facet normals span the plane across the axis.
    M = zeros(3, 3)
    for (a, b, d) in m.triangles
        n = cross(m.vertices[b] - m.vertices[a], m.vertices[d] - m.vertices[a])   # area-weighted
        M .+= n * transpose(n) ./ max(norm(n), eps())
    end
    F = eigen(Symmetric(M))
    F.values[1] <= 1e-3 * F.values[3] || return false
    ax = Vec3(F.vectors[:, 1]...)
    e1 = Vec3(F.vectors[:, 2]...)
    e2 = Vec3(F.vectors[:, 3]...)
    U = [dot(p, e1) for p in P]
    V = [dot(p, e2) for p in P]
    y = [2U 2V ones(length(P))] \ (U .^ 2 .+ V .^ 2)
    rc = sqrt(max(y[3] + y[1]^2 + y[2]^2, 0.0))
    return rc > 0 && maximum(k -> abs(hypot(U[k] - y[1], V[k] - y[2]) - rc), eachindex(P)) <= tol(rc)
end

_with_surface(c::Component, surf) =
    surf === nothing ? c : Component(c.name, c.mesh, c.material; kind = c.kind, two_sided = c.two_sided,
                                     surface = surf)

# The solids of a multi-solid (ASCII) STL file, read once per version.
const _STL_SOLIDS = _LRU{_FileKey,Vector{Pair{String,Mesh}}}(_CAD_FILES_KEPT)

function _stl_solids(file::AbstractString)
    isfile(file) || throw(ArgumentError("STL file not found: $file"))
    key = (abspath(file), mtime(file), filesize(file))
    return _cached(() -> read_stl_solids(file), _STL_SOLIDS, key)
end

# Whether part `p` (element size `h`, placed by its chamber `frame`) is meshed
# already: CAD parts (VisRad CAD, STL) remesh for seconds to minutes and are cached; the
# other kinds mesh at once. A part the geometry rejects counts as ready (the build says why).
function _part_ready(p::PartSpec, h::Real, frame::Frame)
    p.kind in (:visrad, :mesh) || return true
    _iscached(_PART_MEMO, _part_key(p, h, frame, false)) && return true
    try
        R, center = frame.R, frame.origin
        hp = p.values[:h] > 0 ? p.values[:h] : Float64(h)
        toward = transpose(R) * (zero(Vec3) - center)
        if p.kind === :visrad
            o = _visrad_object(p.file, p.object)
            o.kind === :mesh || return true
            return _iscached(_VRW_PARTS, first(_visrad_key(p, o, hp, toward, false)))
        end
        return _iscached(_STL_CACHE, _stl_key(p, hp, toward, false))
    catch err
        err isa ArgumentError || rethrow()
        return true
    end
end

# Names of the CAD parts of a design that are in the calculations and not meshed yet.
function _unmeshed_parts(mesh::MeshSpec, parts, excluded)
    out = Symbol[]
    any(p -> p.kind in (:visrad, :mesh) && !(p.name in excluded), parts) || return out
    frames = try
        part_frames(parts)
    catch err
        err isa ArgumentError || rethrow()
        return out                                    # (the build says why)
    end
    for p in parts
        (p.kind in (:visrad, :mesh) && !(p.name in excluded)) || continue
        _part_ready(p, mesh.h, frames[p.name]) || push!(out, p.name)
    end
    return out
end

# Mesh `parts` (of the design's parts `all`) on worker threads, each once (the part
# caches then serve the target build, its commit and later edits); `raw` parts (names)
# keep their file's facets. Progress goes to `report`.
function _mesh_parts(parts, mesh, raw, report, label; all = parts)
    n = length(parts)
    finished = Threads.Atomic{Int}(0)
    report("$label: meshing $n parts…")
    frames = part_frames(all)
    tasks = map(parts) do p
        Threads.@spawn begin
            part_components(p, mesh.h; frame = frames[p.name], raw = p.name in raw,
                            capsule_h = mesh.capsule_h)
            k = Threads.atomic_add!(finished, 1) + 1
            k < n && report("$label: meshed $k of $n parts…")
        end
    end
    for t in tasks
        try
            fetch(t)
        catch err
            err isa TaskFailedException ? throw(err.task.exception) : rethrow()
        end
    end
    return nothing
end

"""
    add_stl_parts!(ts::TargetSetup, file; values...) -> Bool

Add the parts of an STL file: one `:mesh` part for a file of one solid (as
[`add_part!`](@ref)), or one per named solid of an ASCII file of several (a multi-body CAD
export, `read_stl_solids`), each named after its solid (letters, digits and `_`, made
unique) and placed and sized by `values` alike, as one assembly named after the file
(its frame at TCC, the parts placed in it: [`add_assembly!`](@ref)). The parts are added
together: one the geometry rejects (an empty solid, …) adds none, sets `error` and
returns `false`.
"""
function add_stl_parts!(ts, file::AbstractString; values...)      # ts::TargetSetup (defined below)
    solids = try
        _stl_solids(file)
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
    length(solids) <= 1 && return add_part!(ts, :mesh; file, values...)
    taken = union(Set(all_component_names(ts)), (p.name for p in ts.parts[]), keys(part_frames(ts)))
    asm = _free_identifier(first(splitext(basename(file))), taken, "assembly")
    parts = copy(ts.parts[])
    vals = Dict{Symbol,Float64}(k => Float64(v) for (k, v) in values)
    for (k, (sname, _)) in enumerate(solids)
        name = _free_identifier(sname, taken, "solid")
        push!(parts, PartSpec(:mesh; name, file, object = string(k), values = vals, frame = asm, assembly = asm))
    end
    push!(parts, PartSpec(:assembly; name = asm))
    return _try_build!(ts; parts)
end

# A name from text `s` (letters, digits and `_`, starting with a letter or `_`: `prefix`
# and `_` put in front otherwise) not in `taken` (`name_2`, `name_3`, … when it is),
# added to `taken`.
function _free_identifier(s::AbstractString, taken::AbstractSet, prefix::AbstractString)
    base = replace(strip(s), r"[^A-Za-z0-9_]" => "_")
    (isempty(base) || !occursin(r"^[A-Za-z_]", base)) && (base = prefix * "_" * base)
    name, j = Symbol(base), 1
    while name in taken
        j += 1
        name = Symbol(base, :_, j)
    end
    push!(taken, name)
    return name
end

# The `Target` of `comps`, built once for the same components (identical objects): the
# few latest are kept, so an edit to a part left out of the calculations, or an undo,
# finds its targets again. The targets of a NIF CAD workspace take about a second to
# build, with BVHs of 2 × 10⁵ facets. A target being built on a worker is waited for
# instead of being built twice.
# Keyed by the components' object ids: the entry holds the components, so no other
# object takes their ids while it is kept.
const _TARGETS_KEPT = 6                        # calculation targets and ghosts of 3 designs
const _TARGETS = _LRU{Vector{UInt},Pair{Vector{Component},Target}}(_TARGETS_KEPT)

_target_key(comps::Vector{Component}) = UInt[objectid(c) for c in comps]
_target_built(comps::Vector{Component}) = _iscached(_TARGETS, _target_key(comps))

function _target_of(comps::Vector{Component})
    key = _target_key(comps)
    return _cached(() -> comps => Target(comps), _TARGETS, key).second
end

# A design's components split into the target of those in the calculations and that of
# those left out (`excluded`, drawn as a ghost: `nothing` when none is).
function _split_components(comps::Vector{Component}, excluded)
    keep, gh = _split_lists(comps, excluded)
    return _targets_of(keep, gh)
end

# The components in the calculations and those left out (`nothing` when none is).
function _split_lists(comps::Vector{Component}, excluded)
    # a lip with its hole; a piece with its part
    out(c) = (is_lip(c) ? shown_component(c.name) : c.name) in excluded || c.part in excluded
    keep = [c for c in comps if !out(c)]
    # (a design without parts is the empty target; one whose every component is left out is not)
    isempty(keep) && !isempty(comps) && throw(ArgumentError("every component is left out of the calculations"))
    length(keep) == length(comps) && return keep, nothing
    return keep, [c for c in comps if out(c)]
end

_targets_of(keep, gh) = _target_of(keep), gh === nothing ? nothing : _target_of(gh)

# Whether building the targets of these lists is worth a worker: not built yet, and
# more than `BACKGROUND_TRIANGLES` triangles to build.
function _build_later(keep, gh)
    n = 0
    for cs in (keep, gh)
        (cs === nothing || _target_built(cs)) && continue
        n += sum(c -> length(c.mesh.triangles), cs; init = 0)
    end
    return n > 0 && n > BACKGROUND_TRIANGLES[]
end
