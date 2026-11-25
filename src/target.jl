# This file defines target components and the assembled `Target`, which holds the elements
# and a ray-casting BVH.

"""
    Component(name, mesh, material; kind = :surface, two_sided = false, surface = nothing,
              part = name)

A named piece of a target. `kind = :surface` components are opaque: they block laser
rays and lines of sight, absorb laser light, and radiate from their front side (and
back side when `two_sided`). `kind = :aperture` components (e.g. laser entrance holes)
are transparent diagnostic surfaces that record the x-ray flux crossing them and
re-emit nothing, or, with a `GapClosure` or `LipClosure` albedo, re-emit what a closing lip of plasma
would (see `aperture`). `surface` is the exact `CurvedSurface` the mesh approximates, if any;
the radiation model then uses the curved geometry (see `view_factors`).

`part` names the user's part the component belongs to: the component name, unless it is a
piece of a part split for computation ([`curved_components`](@ref) names the pieces of
a part `p` `p/cyl1`, …, `p/rest`, each with `part = p`). Name-taking functions
(`elements`, `Scenario(...; observe)`, keep-outs, clearance hardware) accept the name of a
part for all its pieces, and reports (`energy_balance`, `component_power`,
`write_results`, `report`) are keyed by part; see [`parts`](@ref).
"""
struct Component
    name::Symbol
    mesh::Mesh
    material::Material
    kind::Symbol
    two_sided::Bool
    surface::Union{Nothing,CurvedSurface}
    part::Symbol
    function Component(name, mesh, material; kind::Symbol = :surface, two_sided::Bool = false,
                       surface::Union{Nothing,CurvedSurface} = nothing, part = name)
        kind in (:surface, :aperture) || throw(ArgumentError("kind must be :surface or :aperture"))
        kind === :aperture && two_sided && throw(ArgumentError("apertures are one-sided"))
        new(Symbol(name), mesh, material, kind, two_sided, surface, Symbol(part))
    end
end

"""
    aperture(name, mesh; surface = nothing, closure = nothing)

Transparent diagnostic surface (e.g. a laser entrance hole): rays pass through it and
its incident x-ray flux is the flux escaping through it.

With `closure`, a `GapClosure`, the aperture is an opening that plasma closes to x-rays
during the pulse while laser light and lines of sight still pass. An example is the lip
of a closing laser entrance hole, `hohlraum(; leh_closure)`. The aperture re-emits the
fraction ``α(t) = s α_{after}`` of its incident flux, and the fraction `1 - s` escapes and
is booked as escaped by `energy_balance`. The aperture does not convert laser light. In
the Monte Carlo reference ([`montecarlo_incident`](@ref)), paths escape through open
apertures and are re-emitted diffusely from a lip with the albedo of the lip (given
`albedo = albedos(target, time)`). With a `LipClosure`, the closed
fraction of each element is `drive_history`'s state instead of a function of time.
"""
aperture(name, mesh; surface = nothing, closure::Union{Nothing,GapClosure,LipClosure} = nothing) =
    Component(name, mesh, closure === nothing ? BLACK : Material(:plasma; albedo = closure, conversion = 0.0);
              kind = :aperture, surface)

"Component with its mesh rigidly moved (`x ↦ R*x + t`)."
transform(c::Component, R::AbstractMatrix, t = zero(Vec3)) =
    Component(c.name, transform(c.mesh, R, Vec3(t...)), c.material; kind = c.kind,
              two_sided = c.two_sided, surface = _transform_surface(c.surface, R, Vec3(t...)), c.part)
translate(c::Component, t) = transform(c, SMatrix{3,3,Float64}(I), t)

"""
    Target(components; closed = false, visibility = :centroid)

Assembled target geometry. Each patch side that radiates becomes an *element* (the
unknowns of the radiation model): one element per patch for one-sided surfaces and
apertures, two for two-sided surfaces. The back side of a one-sided surface absorbs
laser light as a loss but is not modeled radiatively.

`closed = true` declares that the elements (apertures included) form a closed enclosure, so
`RadiationModel` balances its view factors by default (`enforce_closure!`). `hohlraum` sets
it; a target assembled from its components does not, since added parts such as foils and
sombreros may lie outside the enclosure. Declare it (here, or
`cad_target(...; closed = true)`) when apertures cover every opening
([`open_fraction`](@ref) shows where facets leak). On bored entrance holes, waist gaps,
windows and CAD cavities so closed, balancing brings the loss through the apertures within
0.2% of a continuum Monte Carlo with `:sampled` visibility (raw factors up to 0.8% off) and
the capsule drive within 0.07%, and moves the capsule a₂/a₀ by at most 0.08 pp (toward the
reference where it moves it most). With `:centroid`, rims leave the raw loss 2–3% high and
balancing brings it within 0.4%, but moves a₂/a₀ by 0.1–0.3 pp, toward the reference or away
from it. Flat facets must be watertight for it: balancing a cracked mesh adds the exchange
the cracks lose (research notes, "Closure on bored, windowed and CAD cavities").

`visibility` is the `view_factors` visibility rule `RadiationModel` uses by default:
`:centroid` (one ray per pair) or `:sampled` (the visible fraction of pairs partly
shadowed by a rim). Targets with rims set `:sampled`: `cad_target` and `visrad_target`
(bored entrance holes, gaps, CAD walls), where one centroid ray leaves row sums up to
1% off on walls and 6% on LEH apertures; `hohlraum` and assembled targets keep
`:centroid`: on flush entrance holes the largest row-sum error is on the wall, 0.3–0.7%,
and after balancing `:sampled` moves the drive by less than 3e-5 relative and a₂/a₀ by
less than 0.03 pp for 3–4 times the build (research notes, "Partial visibility" and
"Visibility default on the GUI presets").

A `Target` is immutable in all its fields but has identity (a `mutable struct` of
`const` fields), so that caches of derived data (the footprint of a Dante view, the pixels
of an imager) hold it weakly and release their entries once the target is unreachable.
"""
mutable struct Target
    const components::Vector{Component}
    # triangles
    const tri_v0::Vector{Vec3}
    const tri_e1::Vector{Vec3}
    const tri_e2::Vector{Vec3}
    const tri_normal::Vector{Vec3}
    const tri_area::Vector{Float64}
    const tri_comp::Vector{Int32}
    const tri_front::Vector{Int32}       # element on the front side (0 = none)
    const tri_back::Vector{Int32}        # element on the back side (0 = none)
    const tri_opaque::BitVector
    # elements
    const elem_centroid::Vector{Vec3}
    const elem_normal::Vector{Vec3}      # unit normal of the radiating side
    const elem_area::Vector{Float64}
    const elem_comp::Vector{Int32}
    const elem_side::Vector{Int8}        # +1 front, -1 back
    const elem_tris::Vector{UnitRange{Int}}
    const elem_curved::Vector{Int8}      # radiating side relative to the curved surface's outward normal (0 if flat)
    const comp_elems::Vector{Vector{Int}}
    const bvh::BVH
    const scale::Float64                 # bounding-sphere radius about the origin (mm)
    const closed::Bool                   # the elements form a closed enclosure
    const visibility::Symbol             # default `view_factors` visibility rule
end

function Target(components::AbstractVector{Component}; closed::Bool = false,
                visibility::Symbol = :centroid)
    visibility in (:centroid, :sampled) ||
        throw(ArgumentError("visibility must be :centroid or :sampled, got :$visibility"))
    names = [c.name for c in components]
    allunique(names) || throw(ArgumentError("component names must be unique: $names"))
    for c in components
        c.part === c.name && continue
        any(d -> d.name === c.part && d.part !== c.part, components) &&
            throw(ArgumentError("part $(c.part) is also the name of a component of another part"))
    end
    tv0, te1, te2, tn, ta = Vec3[], Vec3[], Vec3[], Vec3[], Float64[]
    tcomp, tfront, tback = Int32[], Int32[], Int32[]
    topaque = BitVector()
    ec, en, ea, ecomp, eside, ecurv = Vec3[], Vec3[], Float64[], Int32[], Int8[], Int8[]
    etris = UnitRange{Int}[]
    comp_elems = Vector{Int}[]
    for (ci, c) in enumerate(components)
        m = c.mesh
        t_offset = length(tv0)
        for (a, b, d) in m.triangles
            v0, v1, v2 = m.vertices[a], m.vertices[b], m.vertices[d]
            cr = cross(v1 - v0, v2 - v0)
            push!(tv0, v0); push!(te1, v1 - v0); push!(te2, v2 - v0)
            ar = norm(cr) / 2
            push!(ta, ar); push!(tn, ar > 0 ? cr / (2ar) : zero(Vec3))
            push!(tcomp, ci); push!(tfront, 0); push!(tback, 0)
            push!(topaque, c.kind === :surface)
        end
        mine = Int[]
        for side in (c.two_sided ? (1, -1) : (1,))
            for r in m.patches
                gr = r .+ t_offset
                A = sum(ta[gr])
                A > 0 || continue
                cen = sum(ta[k] * (tv0[k] + (te1[k] + te2[k]) / 3) for k in gr) / A
                nrm = normalize(sum(ta[k] * tn[k] for k in gr))
                cen = _onto_patch(tv0, te1, te2, gr, cen, nrm)
                curv = Int8(0)
                if c.surface !== nothing && !_at_apex(c.surface, tv0, te1, te2, gr)
                    # Radiation uses the area of the curved patch, with normals on the
                    # side the flat element radiates from.
                    A = c.surface isa DiskSurface ?
                        _sector_area(_sector(c.surface, _patch_vertices(tv0, te1, te2, gr))) :
                        sum(_curved_area(c.surface, tv0[k], tv0[k] + te1[k], tv0[k] + te2[k]) for k in gr)
                    curv = Int8(sign(dot(side * nrm, last(surface_point(c.surface, cen)))))
                end
                push!(ec, cen); push!(en, side * nrm); push!(ea, A); push!(ecurv, curv)
                push!(ecomp, ci); push!(eside, side); push!(etris, gr)
                e = length(ec)
                push!(mine, e)
                for k in gr
                    side == 1 ? (tfront[k] = e) : (tback[k] = e)
                end
            end
        end
        push!(comp_elems, mine)
    end
    opaque_idx = findall(topaque)
    bvh = _sah_bvh(BVH(tv0[opaque_idx], te1[opaque_idx], te2[opaque_idx], Int32.(opaque_idx)))
    scale = isempty(tv0) ? 1.0 :
            maximum(max(norm(v), norm(v + a), norm(v + b)) for (v, a, b) in zip(tv0, te1, te2))
    return Target(collect(components), tv0, te1, te2, tn, ta, tcomp, tfront, tback, topaque,
                  ec, en, ea, ecomp, eside, etris, ecurv, comp_elems, bvh, scale, closed, visibility)
end

Target(cs::Component...; kwargs...) = Target(collect(cs); kwargs...)

# Test whether an element touches the apex of the cone it lies on. The true surface is
# singular at the apex, where the kernel between points on either side grows as 1/s².
# The point rules of curved elements miss this singularity by a fixed fraction however
# fine the mesh; the ring around the apex of a pointed cone received 3.4% too little. Such
# an element therefore stays flat, and near pairs on it use the exact polygon kernel with
# adaptive emitters.
_at_apex(s, tv0, te1, te2, gr) = false
function _at_apex(s::ConeSurface, tv0, te1, te2, gr)
    dmin, dmax = Inf, 0.0
    for k in gr, v in (tv0[k], tv0[k] + te1[k], tv0[k] + te2[k])
        d = norm(v - s.apex)
        dmin, dmax = min(dmin, d), max(dmax, d)
    end
    return dmin <= 1e-3 * dmax
end

# A folded patch (a remeshed CAD element spans facets within a few degrees of one
# plane) has its area centroid off its facets: behind them on a convex surface, where
# rays from the centroid (visibility, imaging) start inside the solid and are blocked by
# the patch itself. Such a centroid is moved along the mean normal onto the facet the
# line through it meets (`remesh_mesh` makes sure it meets one). The centroid of a planar
# patch is left unchanged.
function _onto_patch(tv0, te1, te2, gr, c, n)
    length(gr) > 1 || return c
    h = 0.0
    for k in gr, v in (tv0[k], tv0[k] + te1[k], tv0[k] + te2[k])
        h = max(h, norm(v - c))
    end
    best = Inf
    for k in gr
        nt = cross(te1[k], te2[k])
        den = dot(n, nt)
        den == 0 && continue
        a = tv0[k]
        s = dot(a - c, nt) / den
        p = c + s * n
        b, d = a + te1[k], a + te2[k]
        tol = -1e-9 * dot(nt, nt)
        dot(_tri_cross(a, b, p), nt) >= tol && dot(_tri_cross(b, d, p), nt) >= tol &&
            dot(_tri_cross(d, a, p), nt) >= tol && abs(s) < abs(best) && (best = s)
    end
    return isfinite(best) && abs(best) > 1e-9 * h ? c + best * n : c
end

"Number of radiating elements (patch sides) of a target."
nelements(t::Target) = length(t.elem_area)
ntriangles(t::Target) = length(t.tri_area)

"""
    component_index(target, name) -> Int

Index of the component named `name` (the name of a piece; [`part_pieces`](@ref)
gives the pieces of a part).
"""
function component_index(t::Target, name::Symbol)
    i = findfirst(c -> c.name === name, t.components)
    i === nothing && throw(KeyError(name))
    return i
end

"""
    parts(target) -> Vector{Symbol}

The user's parts of a target, in the order of their first components: the components'
`part`s. The parts of a primitive target are its component names; a CAD part that
[`curved_components`](@ref) split into `p/cyl1`, …, `p/rest` is the one part `p`.
"""
parts(t::Target) = unique!([c.part for c in t.components])

"""
    part_pieces(target, name) -> Vector{Int}

Indices of the components that make up `name`: the pieces of the part `name`, or else
the one component of that name, so a piece stays reachable by its name. A
`KeyError` for a name that is neither.
"""
function part_pieces(t::Target, name::Symbol)
    is = findall(c -> c.part === name, t.components)
    isempty(is) || return is
    return [component_index(t, name)]
end

"""
    elements(target, names...) -> Vector{Int}

Element indices of the named parts or components (both sides of two-sided
components); the name of a part stands for all its pieces (see [`part_pieces`](@ref)).
"""
elements(t::Target, names::Symbol...) =
    reduce(vcat, (t.comp_elems[ci] for n in names for ci in part_pieces(t, n)); init = Int[])

"""
    element_material(target, e)
"""
element_material(t::Target, e::Integer) = t.components[t.elem_comp[e]].material
is_aperture(t::Target, e::Integer) = t.components[t.elem_comp[e]].kind === :aperture

function Base.show(io::IO, ::MIME"text/plain", t::Target)
    ps = parts(t)
    split = length(ps) < length(t.components)
    println(io, "Target with ", length(ps), split ? " parts (" : " components, ",
            split ? string(length(t.components), " components), ") : "", nelements(t),
            " elements, ", ntriangles(t), " triangles")
    for p in ps
        cis = part_pieces(t, p)
        c = t.components[first(cis)]
        ne = sum(i -> length(t.comp_elems[i]), cis)
        A = sum(t.elem_area[e] for i in cis for e in t.comp_elems[i] if t.elem_side[e] == 1; init = 0.0)
        @printf(io, "  %-14s %-8s %-9s %6d elements  area %9.4f mm²%s%s\n", p, c.kind,
                c.material.name, ne, A, c.two_sided ? "  (two-sided)" : "",
                length(cis) > 1 ? "  ($(length(cis)) pieces)" : "")
    end
end
Base.show(io::IO, t::Target) = print(io, "Target(", length(t.components), " components, ",
                                     nelements(t), " elements)")

# Small offset used to start rays just off a surface (mm).
_ray_eps(t::Target) = 1e-7 * max(t.scale, 1.0)

"""
    visible(target, a, b) -> Bool

Whether the open segment between points `a` and `b` is free of opaque surfaces
(apertures are transparent): `visible(t, Vec3(0, 0, 0), Vec3(0, 0, 8))` is `true` when
no opaque surface lies between the target center and a point 8 mm up the axis.
"""
visible(t::Target, a::Vec3, b::Vec3) = !any_hit(t.bvh, a, b - a, 1 - 1e-9)

"""
    element_quadrature(target, e, level) -> (points, weights, normals)

Quadrature over element `e` from `4^level` sub-triangles per triangle: sub-triangle
centroids weighted by area, with the element normal. On a component with a
`CurvedSurface`, the points are projected onto the surface and carry its normals
(on the radiating side of the element) and curved-patch areas; an element of more than two
facets (imported CAD) instead gets as many points as a two-triangle element, from the
pieces of its facets lumped into cells of the element tangent plane. Weights sum to the
element area.
"""
function element_quadrature(t::Target, e::Integer, level::Integer)
    pts = Vec3[]
    ws = Float64[]
    for k in t.elem_tris[e]
        _subtri!(pts, ws, t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k], level)
    end
    s = t.components[t.elem_comp[e]].surface
    if s === nothing || t.elem_curved[e] == 0
        return pts, ws, fill(t.elem_normal[e], length(pts))
    end
    if s isa DiskSurface
        r0, r1, a, b = element_sector(t, e)
        nr, nφ = 2^level, 2^(level + 1)
        rs = sqrt.(range(r0^2, r1^2; length = nr + 1))        # equal-area rings
        φs = range(a, b; length = nφ + 1)
        cpts, cws = Vec3[], Float64[]
        for i in 1:nr, k in 1:nφ
            push!(cpts, _sector_point(s, sqrt((rs[i]^2 + rs[i+1]^2) / 2), (φs[k] + φs[k+1]) / 2))
            push!(cws, (rs[i+1]^2 - rs[i]^2) / 2 * (φs[k+1] - φs[k]))
        end
        return cpts, cws, fill(t.elem_curved[e] * s.normal, length(cpts))
    end
    length(t.elem_tris[e]) > 2 && return _lumped_quadrature(t, e, level, s)
    cpts, cws, cns = Vec3[], Float64[], Vec3[]
    for k in t.elem_tris[e]
        _curved_subtri!(cpts, cws, cns, s, t.elem_curved[e], t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k],
                        t.tri_v0[k] + t.tri_e2[k], level)
    end
    return cpts, cws, cns
end

# Build the quadrature of a curved element made of more than two facets, such as a CAD
# element remeshed across many strips of a cylinder. The quadrature has as many points as
# a two-triangle element would have, wherever the facets lie. It is built in three steps:
# - Each facet is bisected along its longest edge until that edge is under half a cell.
# - The pieces, with exact curved areas on the surface, are lumped into square cells of
#   side size / (2^level √2) in the tangent plane of the element. Here `size` is √2 times
#   the largest centroid-to-vertex distance, as in `view_factors`, so a square element has
#   2·4^level cells.
# - Each cell becomes one point, the projection of the weighted mean point of its pieces,
#   with their total area.
# The weights sum to the exact curved area of the element.
function _lumped_quadrature(t::Target, e::Integer, level::Integer, s::CurvedSurface)
    C, n, sgn = t.elem_centroid[e], t.elem_normal[e], t.elem_curved[e]
    r = 0.0
    for k in t.elem_tris[e], v in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
        r = max(r, norm(v - C))
    end
    hq = r / 2^level                                   # √2 r / (2^level √2)
    e1, e2 = _ortho_basis(n)
    cells = Dict{Tuple{Int,Int},Tuple{Float64,Vec3}}()
    for k in t.elem_tris[e]
        a = t.tri_v0[k]
        _lump_pieces!(cells, s, C, e1, e2, hq, a, a + t.tri_e1[k], a + t.tri_e2[k])
    end
    pts, ws, ns = Vec3[], Float64[], Vec3[]
    for key in sort!(collect(keys(cells)))
        w, wp = cells[key]
        w > 0 || continue
        p, nn = surface_point(s, wp / w)
        push!(pts, p); push!(ws, w); push!(ns, sgn * nn)
    end
    return pts, ws, ns
end

function _lump_pieces!(cells, s, C, e1, e2, hq, a, b, c)
    lab, lbc, lca = sum(abs2, b - a), sum(abs2, c - b), sum(abs2, a - c)
    if max(lab, lbc, lca) > (hq / 2)^2
        if lab >= lbc && lab >= lca
            m = (a + b) / 2
            _lump_pieces!(cells, s, C, e1, e2, hq, a, m, c)
            _lump_pieces!(cells, s, C, e1, e2, hq, m, b, c)
        elseif lbc >= lca
            m = (b + c) / 2
            _lump_pieces!(cells, s, C, e1, e2, hq, a, b, m)
            _lump_pieces!(cells, s, C, e1, e2, hq, a, m, c)
        else
            m = (c + a) / 2
            _lump_pieces!(cells, s, C, e1, e2, hq, a, b, m)
            _lump_pieces!(cells, s, C, e1, e2, hq, m, b, c)
        end
        return
    end
    g = (a + b + c) / 3
    key = (floor(Int, dot(g - C, e1) / hq), floor(Int, dot(g - C, e2) / hq))
    w = _curved_area(s, a, b, c)
    p = first(surface_point(s, g))
    w0, p0 = get(cells, key, (0.0, zero(Vec3)))
    cells[key] = (w0 + w, p0 + w * p)
    return
end

_patch_vertices(tv0, te1, te2, gr) =
    unique(vcat([[tv0[k], tv0[k] + te1[k], tv0[k] + te2[k]] for k in gr]...))

"""
    element_sector(target, e) -> (r0, r1, φa, φb)

Polar sector of an element of a component with a `DiskSurface`.
"""
element_sector(t::Target, e::Integer) =
    _sector(t.components[t.elem_comp[e]].surface,
            _patch_vertices(t.tri_v0, t.tri_e1, t.tri_e2, t.elem_tris[e]))

function _curved_subtri!(pts, ws, ns, s, sgn, a, b, c, level)
    if level == 0
        p, n = surface_point(s, (a + b + c) / 3)
        push!(pts, p); push!(ns, sgn * n); push!(ws, _curved_area(s, a, b, c))
        return
    end
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    _curved_subtri!(pts, ws, ns, s, sgn, a, ab, ca, level - 1)
    _curved_subtri!(pts, ws, ns, s, sgn, ab, b, bc, level - 1)
    _curved_subtri!(pts, ws, ns, s, sgn, ca, bc, c, level - 1)
    _curved_subtri!(pts, ws, ns, s, sgn, ab, bc, ca, level - 1)
end

function _subtri!(pts, ws, a, b, c, level)
    if level == 0
        push!(pts, (a + b + c) / 3)
        push!(ws, norm(cross(b - a, c - a)) / 2)
        return
    end
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    _subtri!(pts, ws, a, ab, ca, level - 1)
    _subtri!(pts, ws, ab, b, bc, level - 1)
    _subtri!(pts, ws, ca, bc, c, level - 1)
    _subtri!(pts, ws, ab, bc, ca, level - 1)
end
