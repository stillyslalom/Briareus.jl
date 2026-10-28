# The cutaway's geometry: the region a sector or plane cut removes, and triangles clipped
# exactly to what remains. No GL: the canvas draws the pieces.
#
# A sector narrower than 180° removes the intersection of two half-spaces, so what
# remains is a union and one set of GPU clip planes cannot express it. The triangles
# are split instead: every triangle of a cut component is split by each plane that
# bounds the removed region (the sector's two half-planes extended through the axis,
# the bounds along the axis, or the cut plane). The removed region is a Boolean
# combination of those planes' sides, so each convex cell of the split lies wholly in
# it or wholly out of it, and the cell's centroid decides which. The kept cells carry
# their triangle's element (its color), the cut's edges are known exactly (for section
# lines), and picking tests the hit point against the same region.

"""
    CutRegion

What a cutaway removes: `kind` (`:none`, `:sector` or `:plane`), the `components` it
cuts, the sector's axis `u` and azimuth basis `e1`, `e2`, center `azimuth` and `width`
(degrees) and bounds `lo`, `hi` along `u`, or the plane's normal `n` and offset `c`
(the removed side is n · p > c), the `planes` (n, c) that bound the region, and the
length scale `tol` below which a point counts as on a plane. Built by
[`cut_region`](@ref).
"""
struct CutRegion
    kind::Symbol
    components::Set{Symbol}
    u::Vec3
    e1::Vec3
    e2::Vec3
    azimuth::Float64
    width::Float64
    lo::Float64
    hi::Float64
    n::Vec3
    c::Float64
    planes::Vector{Tuple{Vec3,Float64}}
    tol::Float64
end

# Azimuth bases about each axis: (axis, from, toward), cyclic.
_cut_basis(a::Symbol) = a === :x ? (Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)) :
                        a === :y ? (Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 0, 0)) :
                        (Vec3(0, 0, 1), Vec3(1, 0, 0), Vec3(0, 1, 0))

"""
    cut_region(v::ViewState; scale = 1.0, on = v.cutaway[]) -> CutRegion

The region the view's cutaway removes (`kind = :none` with the cutaway off, or `on =
false`). `scale` (mm, the target's size) sets the tolerance of the on-plane tests (10⁻⁹
of it).
"""
function cut_region(v::ViewState; scale::Real = 1.0, on::Bool = v.cutaway[])
    u, e1, e2 = _cut_basis(v.cut_axis[])
    lo, hi = v.cut_range[]
    az, w = v.cut_azimuth[], v.cut_width[]
    n = Vec3(v.cut_normal[]...)
    c = dot(n, Vec3(v.cut_point[]...)) + v.cut_shift[]
    kind = on ? v.cut_kind[] : :none
    planes = Tuple{Vec3,Float64}[]
    if kind === :sector
        # the half-planes at azimuths az ∓ w/2, extended through the axis
        if w < 360
            for β in (az - w / 2, az + w / 2)
                push!(planes, (-sind(β) * e1 + cosd(β) * e2, 0.0))
            end
        end
        isfinite(lo) && push!(planes, (u, lo))
        isfinite(hi) && push!(planes, (u, hi))
    elseif kind === :plane
        push!(planes, (n, c))
    end
    return CutRegion(kind, copy(v.cut_components[]), u, e1, e2, az, w, lo, hi, n, c, planes,
                     1e-9 * max(Float64(scale), 1e-3))
end

"""
    is_cut(r::CutRegion, name::Symbol, point) -> Bool

Whether the cut removes component `name` at `point`. A sector removes the points within
`width / 2` of azimuth of its center about its axis (not the axis itself) and within
its bounds along the axis; a plane, the points on the side its normal points to.
"""
function is_cut(r::CutRegion, name::Symbol, p)
    r.kind === :none && return false
    name in r.components || return false
    return _removed(r, Vec3(p...))
end

function _removed(r::CutRegion, p::Vec3)
    if r.kind === :sector
        z = dot(r.u, p)
        r.lo <= z <= r.hi || return false
        x, y = dot(r.e1, p), dot(r.e2, p)
        hypot(x, y) > 0 || return false
        δ = mod(atand(y, x) - r.azimuth + 180, 360) - 180
        return abs(δ) < r.width / 2
    elseif r.kind === :plane
        return dot(r.n, p) > r.c
    end
    return false
end

"""
    cut_pieces(r::CutRegion, a, b, c) -> Vector{NamedTuple}

The triangle (a, b, c) split by the planes of the cut `r` into convex cells, each wholly
removed or kept: `polygon` (its vertices, in the triangle's winding), `edge` (per edge
from vertex i to i + 1: whether it is part of the triangle's own edges) and `removed`.
A triangle that no plane crosses is one cell. The cells' areas sum to the triangle's.
"""
function cut_pieces(r::CutRegion, a, b, c)
    polys = [(Vec3[Vec3(a...), Vec3(b...), Vec3(c...)], Bool[true, true, true])]
    for (n, d) in r.planes
        out = similar(polys, 0)
        for (P, E) in polys
            f = [dot(n, p) - d for p in P]
            if all(>=(-r.tol), f) || all(<=(r.tol), f)
                push!(out, (P, E))
                continue
            end
            for s in (1, -1)
                Q = _clip_side(P, E, f, s, r.tol)
                length(Q[1]) >= 3 && push!(out, Q)
            end
        end
        polys = out
    end
    return [(polygon = P, edge = E, removed = _removed(r, sum(P) / length(P))) for (P, E) in polys]
end

# The part of polygon `P` (edge flags `E`, plane values `f`) on side `s` (1: f ≥ 0, −1:
# f ≤ 0) of a plane, Sutherland–Hodgman, with the edge flags carried: a piece of an
# input edge keeps its flag, the chord along the plane gets `false`.
function _clip_side(P, E, f, s, tol)
    Q, F = Vec3[], Bool[]
    m = length(P)
    side(x) = x > tol ? 1 : x < -tol ? -1 : 0
    for i in 1:m
        j = mod1(i + 1, m)
        si, sj = s * side(f[i]), s * side(f[j])
        if si >= 0
            push!(Q, P[i])
            push!(F, sj < 0 && si == 0 ? false : E[i])
        end
        if si * sj < 0
            push!(Q, P[i] + (f[i] / (f[i] - f[j])) * (P[j] - P[i]))
            push!(F, si > 0 ? false : E[i])
        end
    end
    return Q, F
end

"""
    cut_triangle(r::CutRegion, name, a, b, c) -> (triangles, wires, sections)

What the 3D view draws of the triangle (a, b, c) of component `name` under the cut `r`:
the kept `triangles` (vertex triples; the whole triangle when the cut does not touch
it), the `wires` (the kept parts of the triangle's edges and the cut's edges, for a
wireframe) and the `sections` (segments where the cut crosses the triangle, the edge
of the cut surface).
"""
function cut_triangle(r::CutRegion, name::Symbol, a, b, c)
    tris = NTuple{3,Vec3}[]
    wires = NTuple{2,Vec3}[]
    secs = NTuple{2,Vec3}[]
    A, B, C = Vec3(a...), Vec3(b...), Vec3(c...)
    if r.kind === :none || !(name in r.components)
        push!(tris, (A, B, C))
        push!(wires, (A, B), (B, C), (C, A))
        return tris, wires, secs
    end
    δ = 1e3 * r.tol                        # how far across an edge its other side is probed
    for (P, E, rm) in cut_pieces(r, A, B, C)
        rm && continue
        m = length(P)
        for i in 2:m-1
            push!(tris, (P[1], P[i], P[i+1]))
        end
        cen = sum(P) / m
        for i in 1:m
            p, q = P[i], P[mod1(i + 1, m)]
            mid = (p + q) / 2
            out = mid - cen
            # an edge on a cut plane with the removed region across it: the cut's edge
            onplane = any(((n, d),) -> abs(dot(n, p) - d) <= r.tol && abs(dot(n, q) - d) <= r.tol, r.planes)
            sec = onplane && norm(out) > 0 && _removed(r, mid + δ * out / norm(out))
            sec && push!(secs, (p, q))
            (E[i] || sec) && push!(wires, (p, q))
        end
    end
    return tris, wires, secs
end

"""
    CUT_FLATNESS

How flat a component must be for a cut plane through it to take its normal
([`component_plane`](@ref)): the least length of the area-weighted mean of its
triangles' unit normals (1 for a plane). An engineering default: 0.9 passes a window on
a hohlraum wall (a 1 mm window on a 3 mm radius gives about 0.995) and fails a closed
shell or a cylinder.
"""
const CUT_FLATNESS = 0.9

"""
    component_plane(t::Target, name::Symbol; flatness = CUT_FLATNESS) -> (point, normal)

The plane of part or component `name` (all the pieces of a part, see `part_pieces`):
`point` is its area-weighted centroid and `normal` the area-weighted mean of its
triangles' front normals, normalized, or `nothing` when that mean is shorter than
`flatness` (the component is not flat). An `ArgumentError` for a name not in `t`.
"""
function component_plane(t::Target, name::Symbol; flatness::Real = CUT_FLATNESS)
    any(c -> c.part === name || c.name === name, t.components) ||
        throw(ArgumentError("no component $name in the target"))
    pieces = Set(part_pieces(t, name))
    c, nA, A = zero(Vec3), zero(Vec3), 0.0
    for k in eachindex(t.tri_area)
        t.tri_comp[k] in pieces || continue
        a = t.tri_area[k]
        c += a * (t.tri_v0[k] + (t.tri_e1[k] + t.tri_e2[k]) / 3)
        nA += a * t.tri_normal[k]
        A += a
    end
    A > 0 || throw(ArgumentError("component $name has no surface"))
    return c / A, norm(nA) >= flatness * A ? normalize(nA) : nothing
end

"""
    cut_through_component!(v::ViewState, t::Target, name::Symbol) -> v

Cut with the plane through component `name`'s centroid ([`component_plane`](@ref)),
with the component's normal if it is flat (a window, a foil, a disk) or else the
current normal, at offset 0.
"""
function cut_through_component!(v::ViewState, t::Target, name::Symbol)
    p, n = component_plane(t, name)
    return n === nothing ? set_cut_plane!(v; point = p, offset = 0.0) :
           set_cut_plane!(v; normal = n, point = p, offset = 0.0)
end
