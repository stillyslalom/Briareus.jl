# The radiating cavity around a point: the elements of a target (an imported CAD
# assembly, typically) that bound the free space the point is in.

"""
    Enclosure

Result of `enclosure`: the reduced `target`, and how it was cut from `source`.

- `target`: the components of `source` cut to their kept patches (components left
  empty are dropped).
- `index[i]`: the element of `source` that element `i` of `target` is.
- `kept[e]`: whether element `e` of `source` passed the selection. A two-sided patch
  keeps both sides when either passes, so `target` can hold sides with `kept` false.
- `flux[e]`: the relative flux on element `e` of `source` that decided it (`NaN` for
  elements beyond `radius`, which cast no rays); see `enclosure`.
- `escape`: the fraction of the rays from the point that leave the bounding sphere
  (or reach infinity): the solid-angle fraction of the openings seen from the point.
- `open_area`: Σ over kept elements of area × fraction of the element's rays that
  escape (mm²). By reciprocity this estimates the open area of the cavity (its
  entrance holes, where no aperture or window closes them).
- `leak_area`: the same sum over rays that end on dropped elements (mm²): radiation
  the kept set sends to surfaces the reduced target no longer has.
- `return_area`: what those surfaces would send back, Σ over kept elements of area ×
  the mean over the element's rays of `albedo × flux` at the dropped element hit (mm²,
  in units of the reference flux): the error of the cut, to compare with `open_area`
  (what leaves through the openings) and the cavity's area.
- `back_area`: the sum over rays that end on one-sided surfaces seen from behind
  (mm²), which absorb without radiating in the full target too.
- `back`: the fraction of the rays from the point whose first hit was a one-sided
  surface seen from behind (1 for a point inside the capsule; otherwise a sign of
  wrongly oriented surfaces).
- `rays`: the rays cast per element.
"""
struct Enclosure
    target::Target
    source::Target
    point::Vec3
    index::Vector{Int}
    kept::BitVector
    flux::Vector{Float64}
    escape::Float64
    open_area::Float64
    leak_area::Float64
    return_area::Float64
    back_area::Float64
    back::Float64
    rays::Int
end

"""
    enclosure(target, point = (0, 0, 0); threshold = 0.05, albedo = 0.9, rays = 64,
              points = 8, point_rays = 4096, radius = Inf, seed = 1) -> Enclosure
    enclosure(components, point; kwargs...)

The elements of `target` that bound the free-space region containing `point` (mm), as a
new `Target` (`enc.target`) with each component cut to those patches, with
diagnostics (`Enclosure`; `show` prints elements and area kept per component). Imported
CAD assemblies model a target as closed solids, with outer faces, end faces and
mounting structure that cannot exchange radiation with the cavity, while view factors
cost O(n²) in the elements. Only element sides that face the free space are found
(fronts of one-sided surfaces, either side of two-sided ones), so CAD windings must be
consistent first (`orient_mesh`, `visrad_target(...; orient = true)`).

The selection is the flux an isotropic source at the point would put on each element,
directly and after diffuse re-emission from the surfaces with the uniform albedo
`albedo` (a geometric weight, not the materials'). Every element within `radius` of the
point casts `rays` rays from random points on it in cosine-weighted directions,
stratified in angle, against the target's facets (apertures are transparent, as
everywhere), and the flux is gathered along them:

    Eⱼ = E⁰ⱼ + albedo × (mean of E over the element sides j's rays hit first),

with rays that escape (leave the sphere of `radius` about the point, or hit nothing)
or meet a one-sided surface from behind contributing nothing. `E⁰ⱼ` is the direct flux,
cos θ/(4π d²) averaged over `points` random points on the element, zero where the line
of sight to the point is blocked; lines of sight (and the rays from the point) pass
through surfaces seen from behind, so the point may lie inside the capsule. The
iteration (Gauss–Seidel in element order, to 10⁻⁶) starts from `E⁰`. Elements with
`E ≥ threshold × E_ref` are kept, where `E_ref` is the mean flux on the element sides
that `point_rays` rays from the point (a Fibonacci lattice) first hit: the cavity's own
scale. `enc.flux` holds `E/E_ref`. On a NIF CAD workspace (155k elements) this takes
about 3 s on 8 threads.

Each element gathers through its own rays, which sample its whole view whatever its
size or shape, so a small patch or a sliver deep in a recess or behind a shield is
found through the surfaces around it, which see the cavity. Surfaces outside, such as
outer walls and the far sides of mounts, get no flux at all. A surface that sees the
cavity only through an opening, such as mounting structure beyond an entrance hole,
gets about `albedo` × its view factor to the cavity, plus what its neighbors outside
gather in turn, attenuated by `albedo` at each step, so structure outside does not
carry itself in by seeing itself. The default threshold 0.05 drops a surface facing a
hole of radius a from more than about 4a away on its axis; the rim of a hole in a thick
wall spends about half its view on the cavity and is kept. `threshold = 0` keeps every
element with any flux at all (whatever is connected to the cavity by lines of sight,
however weakly), which `radius` bounds. The random numbers come from `seed` and the
element index, so the result does not depend on the thread count; an element whose
flux is within the sampling noise of the threshold can fall either way.

Removing elements the cavity cannot see changes nothing in the radiosity of the rest;
removing a surface seen through an opening turns the radiation it would have returned
into an escape, which `return_area` estimates against `open_area` (lower `threshold`
when it matters). Thin films that CAD models as solids, such as a capsule tent or an
entrance-hole window, bound the free space like walls: with the point inside a tent the
result is the tent's pocket around the capsule (`open_area` near zero). Leave such
parts out of the target (they transmit x-rays at their thickness) or give a point on
the side wanted.
"""
function enclosure(t::Target, point = (0.0, 0.0, 0.0); threshold::Real = 0.05, albedo::Real = 0.9,
                   rays::Integer = 64, points::Integer = 8, point_rays::Integer = 4096,
                   radius::Real = Inf, seed::Integer = 1)
    rays > 0 && points > 0 && point_rays > 0 ||
        throw(ArgumentError("rays, points and point_rays must be positive"))
    threshold >= 0 || throw(ArgumentError("threshold must not be negative"))
    0 <= albedo < 1 || throw(ArgumentError("albedo must lie in [0, 1)"))
    p = Vec3(point...)
    R = Float64(radius)
    N = nelements(t)
    ε = _ray_eps(t)
    cand = Bool[norm(t.elem_centroid[e] - p) <= R for e in 1:N]
    hits = zeros(Int32, rays, N)
    E0 = zeros(N)
    _cast_elements!(hits, E0, t, cand, p, R, ε, UInt64(seed), points)
    E = _gather_flux(hits, E0, Float64(albedo), threshold == 0)
    firsts, escape, back = _point_rays(t, p, R, ε, point_rays)
    seen = filter(>(0), firsts)
    Eref = isempty(seen) ? 0.0 : sum(E[e] for e in seen) / length(seen)
    Eref > 0 || throw(ArgumentError("the point $(Tuple(p)) sees no element side facing it " *
                                    "($(round(100back; digits = 1))% of its rays first meet a " *
                                    "surface from behind): check the point and the surface " *
                                    "orientation (orient_mesh)"))
    flux = [cand[e] ? E[e] / Eref : NaN for e in 1:N]
    kept = BitVector(cand[e] && (threshold == 0 ? E[e] > 0 : flux[e] >= threshold) for e in 1:N)
    open_area = leak_area = return_area = back_area = 0.0
    for e in 1:N
        kept[e] || continue
        w = t.elem_area[e] / rays
        for h in view(hits, :, e)
            if h == 0
                open_area += w
            elseif h < 0
                back_area += w
            elseif !kept[h]
                leak_area += w
                return_area += w * albedo * E[h] / Eref
            end
        end
    end
    sub, index = _cut_target(t, kept)
    return Enclosure(sub, t, p, index, kept, flux, escape, open_area, leak_area, return_area,
                     back_area, back, Int(rays))
end

enclosure(cs::AbstractVector{<:Component}, point = (0.0, 0.0, 0.0); kwargs...) =
    enclosure(Target(collect(Component, cs)), point; kwargs...)

# splitmix64: a counter-based stream per element, so samples do not depend on threads.
@inline function _splitmix(s::UInt64)
    s += 0x9e3779b97f4a7c15
    z = s
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return s, z ⊻ (z >> 31)
end
@inline function _uniform(s::UInt64)
    s, z = _splitmix(s)
    return s, (z >> 11) * 0x1.0p-53
end

# Orthonormal vectors perpendicular to the unit vector n.
@inline function _tangents(n::Vec3)
    a = abs(n[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)
    b1 = normalize(cross(n, a))
    return b1, cross(n, b1)
end

# Parameter at which the ray o + s d (|d| = 1) leaves the sphere of radius R about p
# (0 when it starts outside and moves away; Inf for an unbounded region).
@inline function _sphere_exit(o, d, p, R)
    isinf(R) && return Inf
    q = o - p
    b = dot(q, d)
    disc = b^2 - (dot(q, q) - R^2)
    return max(-b + sqrt(max(disc, 0.0)), 0.0)
end

# The element side a ray along d hits on triangle k: the front's element when it meets
# the front, else the back's; -1 for an unmodeled (one-sided) back.
@inline function _side_hit(t::Target, k, d)
    e = dot(d, t.tri_normal[k]) < 0 ? t.tri_front[k] : t.tri_back[k]
    return e == 0 ? Int32(-1) : Int32(e)
end

# Diffuse rays and direct flux from every candidate element: column e of `hits`
# (rays × elements) and E0[e] are written only by element e, in blocks handed out
# dynamically (`_tforeach`).
function _cast_elements!(hits::Matrix{Int32}, E0::Vector{Float64}, t::Target, cand::Vector{Bool},
                         p::Vec3, R::Float64, ε::Float64, seed::UInt64, npts::Int)
    N = size(hits, 2)
    block = 64
    _tforeach(cld(N, block)) do b
        for e in (b-1)*block+1:min(b * block, N)
            cand[e] && _cast_element!(hits, E0, t, e, p, R, ε, seed, npts)
        end
    end
    return hits
end

# A uniform random point on element e (a triangle by area, a point on it), the normal
# of its radiating side there, and the advanced random state.
@inline function _element_point(t::Target, e, side, tris, Atot, s)
    s, u0 = _uniform(s); s, u1 = _uniform(s); s, u2 = _uniform(s)
    j, acc = 0, u0 * Atot
    for jj in tris
        t.tri_area[jj] > 0 || continue
        j = jj
        acc -= t.tri_area[jj]
        acc < 0 && break
    end
    a = sqrt(u1)
    x = t.tri_v0[j] + a * (1 - u2) * t.tri_e1[j] + a * u2 * t.tri_e2[j]
    return x, side * t.tri_normal[j], s
end

function _cast_element!(hits, E0, t::Target, e::Int, p::Vec3, R::Float64, ε::Float64, seed::UInt64,
                        npts::Int)
    k = size(hits, 1)
    tris = t.elem_tris[e]
    side = t.elem_side[e]
    Atot = sum(t.tri_area[j] for j in tris)
    s = _splitmix(seed ⊻ (UInt64(e) * 0xd1b54a32d192ed03))[2]
    for r in 1:k
        x, n, s = _element_point(t, e, side, tris, Atot, s)
        s, u3 = _uniform(s); s, u4 = _uniform(s)
        # cosine-weighted direction, stratified in sin²θ
        w = (r - 1 + u3) / k
        st, ct = sqrt(w), sqrt(1 - w)
        φ = 2π * u4
        b1, b2 = _tangents(n)
        d = ct * n + st * (cos(φ) * b1 + sin(φ) * b2)
        o = x + ε * n
        th, id = closest_hit(t.bvh, o, d, _sphere_exit(o, d, p, R))
        hits[r, e] = id == 0 ? Int32(0) : _side_hit(t, id, d)
    end
    # direct flux of a unit isotropic source at p
    E = 0.0
    for _ in 1:npts
        x, n, s = _element_point(t, e, side, tris, Atot, s)
        v = p - x
        L = norm(v)
        c = dot(n, v) / L
        c > 0 || continue
        _first_side(t, p, -v / L, L * (1 - 1e-9) - 2ε, ε)[1] == 0 || continue
        E += c / (4π * L^2)
    end
    E0[e] = E / npts
    return nothing
end

# The first element side a ray from p along d (unit) meets within `left`, passing
# through one-sided surfaces seen from behind (0: none), and whether the first surface
# it met was one of those.
function _first_side(t::Target, p::Vec3, d::Vec3, left::Float64, ε::Float64)
    o = p
    back = false
    for i in 1:256
        left > 0 || break
        th, id = closest_hit(t.bvh, o, d, left)
        id == 0 && break
        e = _side_hit(t, id, d)
        e > 0 && return e, back
        i == 1 && (back = true)
        o += (th + ε) * d
        left -= th + ε
    end
    return Int32(0), back
end

# The element side each of n rays from the point (a Fibonacci lattice) first meets
# (0: none within the sphere of radius R), the escaping fraction, and the fraction
# whose first surface was seen from behind.
function _point_rays(t::Target, p::Vec3, R::Float64, ε::Float64, n::Integer)
    first = zeros(Int32, n)
    nback = 0
    golden = π * (3 - sqrt(5))
    for i in 0:n-1
        z = 1 - (2i + 1) / n
        ρ = sqrt(max(1 - z^2, 0.0))
        d = Vec3(ρ * cos(i * golden), ρ * sin(i * golden), z)
        e, back = _first_side(t, p, d, R, ε)
        first[i+1] = e
        nback += back
    end
    return first, count(==(0), first) / n, nback / n
end

# Gather the flux E = E0 + albedo × (mean of E over each element's ray hits), by
# Gauss–Seidel sweeps in element order from E0 (monotone from below) to a relative
# change of 1e-6; with `positive`, also until no element turns positive.
function _gather_flux(hits::Matrix{Int32}, E0::Vector{Float64}, albedo::Float64, positive::Bool;
                      tol = 1e-6, maxiter::Integer = 100_000)
    k, N = size(hits)
    E = copy(E0)
    w = albedo / k
    for _ in 1:maxiter
        δ = 0.0
        turned = false
        for j in 1:N
            g = 0.0
            @inbounds for r in 1:k
                h = hits[r, j]
                h > 0 && (g += E[h])
            end
            v = E0[j] + w * g
            δ = max(δ, v - E[j])
            turned |= E[j] == 0 && v > 0
            E[j] = v
        end
        δ <= tol * maximum(E) && !(positive && turned) && break
    end
    return E
end

# The components of t cut to the patches with a kept element (both sides of a two-sided
# patch), and the source element of each element of the result.
function _cut_target(t::Target, kept)
    comps = Component[]
    off = 0
    keep_elem = falses(nelements(t))
    for c in t.components
        m = c.mesh
        keep = falses(npatches(m))
        for (pi, r) in enumerate(m.patches)
            isempty(r) && continue
            k = first(r) + off
            f, b = t.tri_front[k], t.tri_back[k]
            keep[pi] = (f > 0 && kept[f]) || (b > 0 && kept[b])
            if keep[pi]
                f > 0 && (keep_elem[f] = true)
                b > 0 && (keep_elem[b] = true)
            end
        end
        off += ntriangles(m)
        any(keep) || continue
        mesh = all(keep) ? m : _select_patches(m, keep)
        push!(comps, Component(c.name, mesh, c.material; kind = c.kind, two_sided = c.two_sided,
                               surface = c.surface))
    end
    return Target(comps; t.visibility), findall(keep_elem)
end

# The mesh made of the patches p with keep[p], in order, sharing vertices as before.
function _select_patches(m::Mesh, keep)
    out = Mesh()
    vmap = zeros(Int, length(m.vertices))
    for (pi, r) in enumerate(m.patches)
        keep[pi] || continue
        t0 = length(out.triangles) + 1
        for k in r
            push!(out.triangles, map(m.triangles[k]) do vi
                vmap[vi] == 0 && (push!(out.vertices, m.vertices[vi]); vmap[vi] = length(out.vertices))
                vmap[vi]
            end)
        end
        push!(out.patches, t0:length(out.triangles))
    end
    return out
end

function Base.show(io::IO, ::MIME"text/plain", enc::Enclosure)
    t = enc.source
    nk = count(enc.kept)
    println(io, "Enclosure of ", Tuple(enc.point), ": ", nelements(enc.target), " of ", nelements(t),
            " elements (", nk, " selected), ", length(enc.target.components), " of ",
            length(t.components), " components")
    @printf(io, "  escape %.4f of the rays from the point; open area %.4g mm², leak area %.4g mm² (returning %.3g mm²), back sides %.4g mm² (%d rays per element)\n",
            enc.escape, enc.open_area, enc.leak_area, enc.return_area, enc.back_area, enc.rays)
    inres = falses(nelements(t))
    inres[enc.index] .= true
    for (i, c) in enumerate(t.components)
        es = t.comp_elems[i]
        isempty(es) && continue
        n = count(e -> inres[e], es)
        A = sum(t.elem_area[e] for e in es)
        Ak = sum(t.elem_area[e] for e in es if inres[e]; init = 0.0)
        @printf(io, "  %-20s %8d / %-8d elements  area %10.4g / %-10.4g mm²\n", c.name, n, length(es), Ak, A)
    end
end
Base.show(io::IO, enc::Enclosure) =
    print(io, "Enclosure(", nelements(enc.target), " of ", nelements(enc.source), " elements)")

"""
    cad_target(parts; orient = true, toward = (0, 0, 0), inside = nothing,
               remesh = nothing, coarsen = nothing, fit = false, material = gold(),
               two_sided = false, visibility = :sampled, enclosure_kw = (;)) -> Target

A radiation-ready target from CAD meshes in one call: what the
[Radiation in an imported CAD hohlraum](@ref) tutorial does by hand, and what
`visrad_target` does with the CAD a VisRad workspace bundles. `parts` holds
`name => mesh` pairs, the CAD parts (as `read_stl_solids` returns them), and
`Component`s, taken as they are (primitives, apertures, parts prepared by hand): those
take part in the cavity selection but are not oriented, remeshed or fitted. The
target's components keep the order of `parts`.

Each CAD part becomes a component of `material` radiating from its front (from both
sides when `two_sided`), through these steps:

1. `orient = true` repairs its windings with [`orient_mesh`](@ref): closed solids face
   out of their material, open sheets face the point `toward` (mm).
2. With `inside`, a point (mm) in the cavity, [`enclosure`](@ref) keeps only the
   elements that bound the free space around it, among the facets of all the parts as
   oriented (`enclosure_kw` are its keywords). The selection needs consistent windings.
3. `remesh` (mm) divides the kept facets into compact, near-planar elements of about
   that size ([`remesh_mesh`](@ref)); `coarsen` (mm) instead merges whole facets into
   elements ([`coarsen_mesh`](@ref)). Give at most one. Only the kept facets are
   divided, which is much cheaper when most of an assembly lies outside the cavity.
   With `inside`, the selection is then repeated on the elements: an element straddling
   the edge of the lit region averages its flux over its own area.
4. `fit = true` (or a named tuple of keywords of [`curved_components`](@ref) and
   [`fit_surfaces`](@ref), such as `minfraction = 0` to fit even the curved faces of
   under 10⁻³ of a part's area, which stay flat by default) splits each part into its
   cylinders, cones, spheres and the rest, fitted on the kept CAD facets, whose vertices
   lie on the true surfaces, so the radiation model integrates over the curved walls.
   The pieces are components of their own, named `<part>/cyl1`, …, `<part>/rest`:
   reports, `elements` and the keep-outs see them, not the part (`SphereKeepOut` by
   default takes every fitted sphere). Fitting is off by default for that reason.

The target's view factors default to `visibility = :sampled` (see `view_factors`): CAD
cavities have rims (entrance-hole bores, gaps) that partly shadow element pairs, which
a single centroid ray decides all or nothing.

`material`, `two_sided`, `remesh` and `coarsen` take one value for every CAD part, a
function of the part's name (a `Symbol`), or a dictionary keyed by name (`Symbol` or
`String`). A part missing from a `remesh` or `coarsen` dictionary keeps its facets, one
missing from a `two_sided` dictionary is one-sided, and one missing from a `material`
dictionary is an error. For example, `remesh = n -> n === :capsule ? 0.15 : 0.5`
remeshes a capsule finer than the walls. For the account of what the cavity selection
drops, call `enclosure` on the target built without `inside`.
"""
function cad_target(parts; orient::Bool = true, toward = (0.0, 0.0, 0.0), inside = nothing,
                    remesh = nothing, coarsen = nothing, fit = false, material = gold(),
                    two_sided = false, visibility = :sampled, enclosure_kw = (;))
    coarsen === nothing || remesh === nothing ||
        throw(ArgumentError("give one of `coarsen` and `remesh`, not both"))
    fit isa Union{Bool,NamedTuple} ||
        throw(ArgumentError("`fit` must be true, false or a named tuple of curved_components keywords"))
    comps = Component[]
    cad = Bool[]                         # per component: a CAD part
    for p in parts
        if p isa Component
            push!(comps, p)
            push!(cad, false)
        elseif p isa Pair && last(p) isa Mesh
            name = Symbol(first(p))
            mat = _part_value(material, name, nothing)
            mat isa Material || throw(ArgumentError("no material for CAD part $name"))
            mesh = orient ? orient_mesh(last(p); toward) : last(p)
            push!(comps, Component(name, mesh, mat; two_sided = Bool(_part_value(two_sided, name, false))))
            push!(cad, true)
        else
            throw(ArgumentError("a part must be a `name => mesh` pair or a Component, not a $(typeof(p))"))
        end
    end
    isempty(comps) && throw(ArgumentError("no parts given"))
    select = inside !== nothing
    if select
        enc = enclosure(comps, inside; enclosure_kw...)
        comps, cad = enc.target.components, cad[_survivors(enc)]
    end
    facets = [c.mesh for c in comps]     # the CAD facets kept: the reference for the fits
    redone = falses(length(comps))       # the parts divided into new elements
    if remesh !== nothing || coarsen !== nothing
        comps = copy(comps)
        for (i, c) in enumerate(comps)
            cad[i] || continue
            sr = remesh === nothing ? nothing : _part_value(remesh, c.name, nothing)
            sc = coarsen === nothing ? nothing : _part_value(coarsen, c.name, nothing)
            sr === nothing && sc === nothing && continue
            m = sc === nothing ? remesh_mesh(c.mesh; size = sr) : coarsen_mesh(c.mesh; size = sc)
            comps[i] = Component(c.name, m, c.material; two_sided = c.two_sided)
            redone[i] = true
        end
        if select
            # select again at the elements' own granularity
            enc = enclosure(comps, inside; enclosure_kw...)
            s = _survivors(enc)
            comps, cad, facets, redone = enc.target.components, cad[s], facets[s], redone[s]
        end
    end
    fit === false && return Target(comps; visibility)
    fkw = fit === true ? (;) : fit
    return Target(reduce(vcat, (cad[i] ? curved_components(comps[i]; reference = redone[i] ? facets[i] : nothing, fkw...) :
                                [comps[i]] for i in eachindex(comps)); init = Component[]); visibility)
end

# A per-part setting: one value, a function of the part's name, or a dictionary by name.
_part_value(x::AbstractDict, name::Symbol, default) =
    haskey(x, name) ? x[name] : get(x, String(name), default)
_part_value(x::Function, name::Symbol, default) = x(name)
_part_value(x, name::Symbol, default) = x

# The components of `enc.source` that `enc.target` keeps (those with a kept element), in order.
_survivors(enc::Enclosure) = sort!(unique(enc.source.elem_comp[enc.index]))
