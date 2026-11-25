# This file defines patch-structured triangle meshes and the geometric primitives that
# targets are built from.
#
# A *patch* is a planar convex polygon stored as a fan of triangles; patches become the
# radiosity elements of a `Target`. The triangle winding defines the front side of the
# patch, with front normal (v2 - v1) × (v3 - v1).

"""
    Vec3

A point or direction in 3D, in mm for positions, stored as `SVector{3,Float64}`.
"""
const Vec3 = SVector{3,Float64}

"""
    Mesh(vertices, triangles, patches)

Triangle mesh whose triangles are grouped into planar patches. `patches[p]` is the
contiguous range of triangle indices forming patch `p`. Lengths are in mm.
"""
struct Mesh
    vertices::Vector{Vec3}
    triangles::Vector{NTuple{3,Int}}
    patches::Vector{UnitRange{Int}}
end

Mesh() = Mesh(Vec3[], NTuple{3,Int}[], UnitRange{Int}[])

"Number of patches of a mesh."
npatches(m::Mesh) = length(m.patches)
ntriangles(m::Mesh) = length(m.triangles)

function Base.show(io::IO, m::Mesh)
    print(io, "Mesh(", npatches(m), " patches, ", ntriangles(m), " triangles)")
end

triangle_vertices(m::Mesh, t::Integer) = (m.vertices[m.triangles[t][1]],
                                          m.vertices[m.triangles[t][2]],
                                          m.vertices[m.triangles[t][3]])

_tri_cross(a, b, c) = cross(b - a, c - a)

function _polygon_normal(poly)
    n = zero(Vec3)
    for k in 2:length(poly)-1
        n += _tri_cross(poly[1], poly[k], poly[k+1])
    end
    return n
end

"""
    add_polygon!(mesh, poly; facing = nothing)

Appends the convex planar polygon `poly` as one fan-triangulated patch. If `facing` is
a vector, the winding is reversed when needed so the front normal of the patch points along
`facing`.
"""
function add_polygon!(m::Mesh, poly::AbstractVector; facing = nothing)
    pts = Vec3.(poly)
    if facing !== nothing && dot(_polygon_normal(pts), facing) < 0
        reverse!(pts)
    end
    i0 = length(m.vertices)
    append!(m.vertices, pts)
    t0 = length(m.triangles) + 1
    for k in 2:length(pts)-1
        push!(m.triangles, (i0 + 1, i0 + k, i0 + k + 1))
    end
    push!(m.patches, t0:length(m.triangles))
    return m
end

"""
    merge_meshes(meshes...)

Concatenates meshes into one, preserving patch structure and order.
"""
function merge_meshes(ms::Mesh...)
    out = Mesh()
    for m in ms
        nv = length(out.vertices)
        nt = length(out.triangles)
        append!(out.vertices, m.vertices)
        append!(out.triangles, [(a + nv, b + nv, c + nv) for (a, b, c) in m.triangles])
        append!(out.patches, [r .+ nt for r in m.patches])
    end
    return out
end

"""
    CurvedSurface

The exact surface that a faceted component approximates, one of `SphereSurface`,
`CylinderSurface`, `ConeSurface`, `TorusSurface` and `DiskSurface`. When a component has
one, the radiation model integrates over the true surface, using the points, normals and
areas of the curved patches. Ray casting and visibility still use the flat facets. On a
curved concave surface the view-factor kernel stays bounded, and flat facets would
misrepresent the exchange between nearby elements.
"""
abstract type CurvedSurface end

"""
    SphereSurface(radius; center = (0, 0, 0))
"""
struct SphereSurface <: CurvedSurface
    center::Vec3
    radius::Float64
end
SphereSurface(radius::Real; center = (0.0, 0.0, 0.0)) = SphereSurface(Vec3(center...), Float64(radius))

"""
    CylinderSurface(radius; origin = (0, 0, 0), axis = (0, 0, 1))

Infinite circular cylinder through `origin` along `axis`.
"""
struct CylinderSurface <: CurvedSurface
    origin::Vec3
    axis::Vec3
    radius::Float64
end
CylinderSurface(radius::Real; origin = (0.0, 0.0, 0.0), axis = (0.0, 0.0, 1.0)) =
    CylinderSurface(Vec3(origin...), normalize(Vec3(axis...)), Float64(radius))

"""
    ConeSurface(half_angle; apex = (0, 0, 0), axis = (0, 0, 1))

One nappe of a right circular cone. Its points are seen from `apex` at `half_angle`
degrees (0 < θ < 90) from `axis`, at distance z tan θ from the axis a distance z along
it. Entrance-hole cones, conical shields and tapered bores have this shape. As on a
cylinder, points are projected onto the cone radially, keeping their height along the
axis. The chords of a circle about the axis, such as the rims of a faceted cone,
therefore land on the circle, and a faceted frustum covers the true one exactly. The
facet patches have exact areas.
"""
struct ConeSurface <: CurvedSurface
    apex::Vec3
    axis::Vec3
    half_angle::Float64             # degrees
    function ConeSurface(apex, axis, half_angle)
        0 < half_angle < 90 || throw(ArgumentError("a cone's half-angle must lie between 0 and 90 degrees"))
        return new(Vec3(apex...), normalize(Vec3(axis...)), Float64(half_angle))
    end
end
ConeSurface(half_angle::Real; apex = (0.0, 0.0, 0.0), axis = (0.0, 0.0, 1.0)) = ConeSurface(apex, axis, half_angle)

"""
    TorusSurface(major_radius, minor_radius; center = (0, 0, 0), axis = (0, 0, 1))

A torus, the points at `minor_radius` r from its spine. The spine is the circle of
`major_radius` R about `axis` through `center`. A CAD fillet is a quarter of the tube,
where a wall meets a cap. Points are projected in their meridian half-plane:

- A point inside the tube circle is projected to the nearest point of the circle. For a
  ring torus, with R > r, this is the nearest point of the torus.
- A point outside the tube circle is projected as on a cone, radially from the axis at
  its height, onto the side of the tube it is on. Beyond the top or bottom of the tube
  it goes to the pole.

Normals point away from the spine. The rims of a faceted fillet are chords of circles
about the axis, which lie at the rim height toward the axis. The chords lie outside the
tube for rims on the poles, where a quarter fillet meets a cap, and on the inner half.
They lie inside the tube, along the radius, on the outer equator, where a fillet meets a
wall. On all of these rims the chords land on the rim circle, and the elements of a
fillet cover the torus band exactly. Elsewhere on the outer half, where a fillet meets a
cone, the rim chords stay off the circle by their sag. The band is then covered to
second order in the facet size. A projection of points cannot both cover those rims
exactly and keep facets near the poles second order. The facet patches have exact
areas under the projection.
"""
struct TorusSurface <: CurvedSurface
    center::Vec3
    axis::Vec3
    major_radius::Float64
    minor_radius::Float64
    function TorusSurface(center, axis, R, r)
        (R > 0 && r > 0) || throw(ArgumentError("a torus needs positive major and minor radii"))
        return new(Vec3(center...), normalize(Vec3(axis...)), Float64(R), Float64(r))
    end
end
TorusSurface(R::Real, r::Real; center = (0.0, 0.0, 0.0), axis = (0.0, 0.0, 1.0)) = TorusSurface(center, axis, R, r)

"""
    DiskSurface(; center = (0, 0, 0), normal = (0, 0, 1))

Plane through `center` whose elements are polar sectors about `center`, as produced by
`annulus_mesh` and `disk_mesh`: the exact shape of each element is the annular sector
spanned by its vertices, with circular arcs where the mesh has chords. Giving the end caps
and laser entrance holes of a hohlraum this surface makes them meet a `CylinderSurface` wall
exactly at the rim.
"""
struct DiskSurface <: CurvedSurface
    center::Vec3
    normal::Vec3
end
DiskSurface(; center = (0.0, 0.0, 0.0), normal = (0.0, 0.0, 1.0)) =
    DiskSurface(Vec3(center...), normalize(Vec3(normal...)))

"""
    surface_point(surface, p) -> (point, normal)

Returns the projection of `p` onto the surface and the outward unit normal at that
point. The projection is radial from the center or the axis. On a torus it is radial
from the spine inside the tube and from the axis outside it. The normal points away from
the center, the axis or the spine.
"""
function surface_point(s::SphereSurface, p)
    d = normalize(p - s.center)
    return s.center + s.radius * d, d
end
function surface_point(s::CylinderSurface, p)
    d = p - s.origin
    a = dot(d, s.axis)
    n = normalize(d - a * s.axis)
    return s.origin + a * s.axis + s.radius * n, n
end

# Returns the height z of a point along the cone axis from the apex, its distance ρ from
# the axis, and the unit vector from the axis toward it. A point on the axis gets an
# arbitrary perpendicular.
function _cone_coordinates(s::ConeSurface, p)
    d = p - s.apex
    z = dot(d, s.axis)
    q = d - z * s.axis
    ρ = norm(q)
    if ρ > 0
        return z, ρ, q / ρ
    end
    a = s.axis
    return z, ρ, first(_ortho_basis(a))
end

function surface_point(s::ConeSurface, p)
    sα, cα = sincosd(s.half_angle)
    z, ρ, ρ̂ = _cone_coordinates(s, p)
    return s.apex + max(z, 0.0) * (s.axis + (sα / cα) * ρ̂), cα * ρ̂ - sα * s.axis
end

# Returns the offset w of a point from the spine point in its meridian half-plane, that
# spine point, and ρ̂, the unit vector from the axis toward the point.
function _torus_coordinates(s::TorusSurface, p)
    d = p - s.center
    z = dot(d, s.axis)
    q = d - z * s.axis
    ρ = norm(q)
    a = s.axis
    ρ̂ = ρ > 0 ? q / ρ : first(_ortho_basis(a))
    c = s.center + s.major_radius * ρ̂
    return p - c, c, ρ̂
end

# Returns (cos ψ, sin ψ) of the tube point onto which a point is projected. The point
# lies at (x, z) from the spine point in its meridian half-plane.
# - Inside the tube circle, x² + z² < r², the result is the nearest point, radially from
#   the spine point.
# - On or outside the circle, the result is the point at the same height z on the same
#   side of the spine, x' = ±√(r² − z²) with the sign of x. When |z| ≥ r it is the pole.
# The two maps agree on the circle, so the map is continuous. A rim is a circle ψ = ψ0
# about the axis. Its chords lie at the rim height, displaced from it toward the axis.
# For ψ0 in [π/2, 3π/2] that is outside the tube, and the chords map onto the rim. This
# range includes the poles ψ0 = ±π/2, where the rim of a quarter fillet meets a cap. On
# the outer equator, ψ0 = 0, the chords are inside the tube on the rim radius from the
# spine point. Facet points outside the tube lie horizontally within the chord sag of it,
# so the horizontal map moves them by at most that sag. The error is therefore second
# order, as for the nearest point.
function _torus_foot(r, x, z)
    l = hypot(x, z)
    l < r && return x / l, z / l
    abs(z) >= r && return 0.0, copysign(1.0, z)
    return copysign(sqrt((r - z) * (r + z)), x) / r, z / r
end

function surface_point(s::TorusSurface, p)
    w, c, ρ̂ = _torus_coordinates(s, p)
    x, z = dot(w, ρ̂), dot(w, s.axis)
    if x == 0 && z == 0
        return c + s.minor_radius * ρ̂, ρ̂
    end
    cψ, sψ = _torus_foot(s.minor_radius, x, z)
    n = cψ * ρ̂ + sψ * s.axis
    return c + s.minor_radius * n, n
end

_transform_surface(::Nothing, R, t) = nothing
_transform_surface(s::SphereSurface, R, t) = SphereSurface(Vec3(R * s.center + t), s.radius)
_transform_surface(s::CylinderSurface, R, t) =
    CylinderSurface(Vec3(R * s.origin + t), normalize(Vec3(R * s.axis)), s.radius)
_transform_surface(s::ConeSurface, R, t) = ConeSurface(Vec3(R * s.apex + t), Vec3(R * s.axis), s.half_angle)
_transform_surface(s::DiskSurface, R, t) = DiskSurface(Vec3(R * s.center + t), normalize(Vec3(R * s.normal)))
_transform_surface(s::TorusSurface, R, t) =
    TorusSurface(Vec3(R * s.center + t), Vec3(R * s.axis), s.major_radius, s.minor_radius)

surface_point(s::DiskSurface, p) = (p - dot(p - s.center, s.normal) * s.normal, s.normal)

_disk_basis(s::DiskSurface) = _ortho_basis(s.normal)

"""
    _sector(s::DiskSurface, vertices) -> (r0, r1, φa, φb)

Polar sector spanned by the vertices of an element: radii from the smallest to the largest
vertex radius (zero when a vertex sits at the center), angles from the first to the last
vertex angle, unwrapped about the mean direction of the element.
"""
function _sector(s::DiskSurface, verts)
    e1, e2 = _disk_basis(s)
    rs = Float64[]
    φs = Float64[]
    for v in verts
        d = v - s.center
        r = hypot(dot(d, e1), dot(d, e2))
        push!(rs, r)
        r > 1e-12 * maximum(norm(w - s.center) for w in verts) && push!(φs, atan(dot(d, e2), dot(d, e1)))
    end
    ref = atan(sum(sin, φs), sum(cos, φs))
    rel = [mod(φ - ref + π, 2π) - π for φ in φs]
    r0 = minimum(rs) < 1e-9 * maximum(rs) ? 0.0 : minimum(rs)
    return r0, maximum(rs), ref + minimum(rel), ref + maximum(rel)
end

_sector_area((r0, r1, a, b)) = (r1^2 - r0^2) / 2 * (b - a)
_sector_point(s::DiskSurface, r, φ) = (e = _disk_basis(s); s.center + r * (cos(φ) * e[1] + sin(φ) * e[2]))

# Returns the triangles of the sector, with each arc divided into `m` straight segments.
function _sector_triangles(s::DiskSurface, (r0, r1, a, b); m = 16)
    φs = range(a, b; length = m + 1)
    tris = NTuple{3,Vec3}[]
    for k in 1:m
        o1, o2 = _sector_point(s, r1, φs[k]), _sector_point(s, r1, φs[k+1])
        if r0 == 0
            push!(tris, (s.center, o1, o2))
        else
            i1, i2 = _sector_point(s, r0, φs[k]), _sector_point(s, r0, φs[k+1])
            push!(tris, (i1, o1, o2), (i1, o2, i2))
        end
    end
    return tris
end

# Returns the area of the patch obtained by projecting the flat triangle (a, b, c) onto
# the surface.
# Radial projection maps a flat triangle onto a spherical triangle exactly, whose area is
# its solid angle times R².
function _curved_area(s::SphereSurface, a, b, c)
    x, y, z = a - s.center, b - s.center, c - s.center
    lx, ly, lz = norm(x), norm(y), norm(z)
    num = abs(dot(x, cross(y, z)))
    den = lx * ly * lz + dot(x, y) * lz + dot(x, z) * ly + dot(y, z) * lx
    return 2 * atan(num, den) * s.radius^2
end

# On a cylinder the projection keeps the axial coordinate z and maps each point to its
# azimuth φ about the axis. By Green's theorem around the projected edges, the patch has
# area R ∬ dφ dz = R |∮ z dφ|. Along a straight edge from P to Q, let u(t) = u_P + tΔ be
# the part of the point perpendicular to the axis. Then dφ/dt = C/|u|² with
# C = (u_P × u_Q)·axis, and z is linear in t, so
#     ∫ z dφ = z_P Δφ + (z_Q - z_P) T,   T = ∫₀¹ t C dt / (E + Bt + At²),
# with A = |Δ|², B = 2u_P·Δ, E = |u_P|² and Δφ the angle from u_P to u_Q. In closed form,
# with G = |u_Q|² - |u_P|² = Δ·(u_Q + u_P) = A + B,
#     T = Δφ/2 + [C log(1 + G/E) - G Δφ] / 2A.
# For a chord of the circle, with both ends on the surface as in `cylinder_mesh`, G = 0
# and ∫ z dφ is the trapezoid (z_P + z_Q)Δφ/2. The curved patch of a triangle with its
# vertices on the cylinder therefore has R times the area of the straight triangle in
# (φ, z). Writing B through G avoids the cancellation in 2u_P·Δ, where terms of size R|Δ|
# sum to |Δ|². That cancellation would cost a short chord ε R/|Δ| of relative accuracy.
# The bracket still loses digits ~|B|/A for short edges seen obliquely from the axis. For
# those edges the poles of the integrand lie far from [0, 1], and Gauss–Legendre is exact
# within round-off.
const _GL16 = let n = 16
    F = eigen(SymTridiagonal(zeros(n), [k / sqrt(4k^2 - 1) for k in 1:n-1]))
    (; x = (F.values .+ 1) ./ 2, w = F.vectors[1, :] .^ 2)          # on [0, 1], Σw = 1
end

# Returns ∫ z dφ along the edge from P to Q. The arguments are the perpendicular parts
# u_P and u_Q, the perpendicular part Δ of the edge, z_P, and the axial rise dz. Δ and dz
# are taken from Q - P directly, which is exact for nearby points; differencing u_Q - u_P
# would carry the ε R rounding of each projection.
function _cylinder_edge_integral(uP, uQ, Δ, zP, dz, ax)
    C = dot(cross(uP, Δ), ax)          # = (u_P × u_Q)·axis, without its cancellation
    Δφ = atan(C, dot(uP, uQ))
    A = dot(Δ, Δ)
    A == 0 && return zP * Δφ
    E = dot(uP, uP)
    G = dot(Δ, uQ + uP)
    B = G - A
    T = if abs(B) <= 16A
        Δφ / 2 + (C * log1p(G / E) - G * Δφ) / 2A
    else
        s = 0.0
        for k in eachindex(_GL16.x)
            x = _GL16.x[k]
            s += _GL16.w[k] * x * C / (E + B * x + A * x^2)
        end
        s
    end
    return zP * Δφ + dz * T
end

function _curved_area(s::CylinderSurface, a, b, c)
    ax = s.axis
    perp(d) = d - dot(d, ax) * ax
    ua, ub, uc = perp(a - s.origin), perp(b - s.origin), perp(c - s.origin)
    ab, bc, ca = b - a, c - b, a - c
    # Σ Δφ = 0 around the patch, so z may be measured from its first vertex. This keeps a
    # small patch far from z = 0 from canceling z Δφ terms of size |z| Δφ.
    zb, zc = dot(ab, ax), dot(c - a, ax)
    acc = _cylinder_edge_integral(ua, ub, perp(ab), 0.0, zb, ax) +
          _cylinder_edge_integral(ub, uc, perp(bc), zb, dot(bc, ax), ax) +
          _cylinder_edge_integral(uc, ua, perp(ca), zc, dot(ca, ax), ax)
    return s.radius * abs(acc)
end

# On a cone the projection keeps the height z along the axis from the apex and the
# azimuth φ. The surface element is ρ dφ dz / cos θ = (sin θ / cos²θ) z dz dφ, so by
# Green's theorem the patch has area (sin θ / 2cos²θ) |∮ z² dφ|. Along a straight edge
# from P to Q, z is linear and dφ/dt = C/|u|² as on the cylinder. ∫ z^k dφ is computed by
# Gauss–Legendre. The edge is halved until the poles of 1/|u|², at t = (−B ± 2i|C|)/2A,
# lie outside the Bernstein ellipse of parameter 5 of each piece. On such pieces 16 points
# converge to rounding level, and 8 points do beyond parameter 20. A chord that is short
# compared with its distance from the axis does not need halving. With z = z₀ + ζ, z₀ the
# lowest vertex, and ∮ dφ = 0 around the patch,
#     ∮ z² dφ = ∮ ζ² dφ + 2 z₀ ∮ ζ dφ.
# This keeps a small patch far from the apex from canceling terms of size z₀² Δφ. Such a
# patch lies on a slightly tapered cone of radius ρ₀ = z₀ tan θ. As θ → 0 the area tends
# to the cylinder area ρ₀ ∮ ζ dφ.
const _GL8 = let n = 8
    F = eigen(SymTridiagonal(zeros(n), [k / sqrt(4k^2 - 1) for k in 1:n-1]))
    (; x = (F.values .+ 1) ./ 2, w = F.vectors[1, :] .^ 2)
end

function _curved_area(s::ConeSurface, a, b, c)
    sα, cα = sincosd(s.half_angle)
    ax = s.axis
    perp(d) = d - dot(d, ax) * ax
    vs = (a, b, c)
    zs = map(v -> dot(v - s.apex, ax), vs)
    i0 = argmin(zs)
    z0, v0 = zs[i0], vs[i0]
    acc1 = 0.0                          # ∮ ζ dφ
    acc2 = 0.0                          # ∮ ζ² dφ
    for k in 1:3
        P, Q = vs[k], vs[mod1(k + 1, 3)]
        uP = perp(P - s.apex)
        Δ = perp(Q - P)
        C = dot(cross(uP, Δ), ax)
        C == 0 && continue
        i1, i2 = _cone_edge_integral(uP, Δ, C, dot(P - v0, ax), dot(Q - P, ax), 0.0, 1.0, 0)
        acc1 += i1
        acc2 += i2
    end
    return abs(sα / (2cα^2) * acc2 + (z0 * sα / cα^2) * acc1)
end

# Returns (∫ ζ dφ, ∫ ζ² dφ) over t ∈ [t0, t1] of the edge u(t) = u_P + tΔ,
# ζ(t) = ζP + dz t.
function _cone_edge_integral(uP, Δ, C, ζP, dz, t0, t1, depth)
    A = dot(Δ, Δ)
    B = 2dot(uP, Δ)
    # Map a pole of 1/|u|² to [-1, 1] and find the Bernstein ellipse through it.
    x = (2 * complex(-B / 2A, abs(C) / A) - t0 - t1) / (t1 - t0)
    ρB = abs(x + sqrt(x - 1) * sqrt(x + 1))
    ρB = max(ρB, 1 / ρB)
    if ρB < 5 && depth < 40
        tm = (t0 + t1) / 2
        a1, a2 = _cone_edge_integral(uP, Δ, C, ζP, dz, t0, tm, depth + 1)
        b1, b2 = _cone_edge_integral(uP, Δ, C, ζP, dz, tm, t1, depth + 1)
        return a1 + b1, a2 + b2
    end
    GL = ρB < 20 ? _GL16 : _GL8
    i1, i2 = 0.0, 0.0
    for k in eachindex(GL.x)
        t = t0 + (t1 - t0) * GL.x[k]
        u = uP + t * Δ
        ζ = ζP + dz * t
        dφ = GL.w[k] * C / dot(u, u)
        i1 += ζ * dφ
        i2 += ζ^2 * dφ
    end
    return (t1 - t0) * i1, (t1 - t0) * i2
end

# On a torus the projection keeps the azimuth φ and maps the point to a tube angle ψ,
# measured from the equator outward (`_torus_foot`). The surface element is
# r (R + r cos ψ) dψ dφ, so by Green's theorem the patch has area
# r |∮ (R ψ + r sin ψ) dφ|. Along a straight edge, dφ/dt = C/|u|² as on the cylinder. ψ is
# smooth wherever the edge stays off the spine and the axis and on one side of the tube
# circle, and pieces are split where the edge crosses it. ∮ is computed by Gauss–Legendre
# on pieces no longer than half the smaller of r and the distance of the edge from the
# axis. On such a piece the nearest singularity lies outside the Bernstein ellipse of
# parameter ~4, so 16 points converge to rounding level. ψ is measured from its value at
# the first vertex, and ∮ dφ = 0 around the patch. This keeps a small patch from
# canceling terms of size R ψ Δφ.
function _curved_area(s::TorusSurface, a, b, c)
    ax, R, r = s.axis, s.major_radius, s.minor_radius
    perp(d) = d - dot(d, ax) * ax
    ca, sa = _torus_foot(r, norm(perp(a - s.center)) - R, dot(a - s.center, ax))
    ψa = atan(sa, ca)
    acc = 0.0
    vs = (a, b, c)
    for k in 1:3
        P, Q = vs[k], vs[mod1(k + 1, 3)]
        uP = perp(P - s.center)
        Δ = perp(Q - P)
        C = dot(cross(uP, Δ), ax)
        C == 0 && continue
        zP, dz = dot(P - s.center, ax), dot(Q - P, ax)
        L = norm(Q - P)
        h = min(r, norm(uP), norm(uP + Δ))
        n = h > 0 ? clamp(ceil(Int, 2L / h), 1, 256) : 256
        for j in 1:n
            acc += _torus_edge_piece(uP, Δ, C, zP, dz, R, r, ψa, sa, (j - 1) / n, j / n)
        end
    end
    return r * abs(acc)
end

# Returns the integrand (R ψ + r sin ψ) dφ/dt at t on the edge, with ψ relative to ψa and
# sin ψ relative to sa, and whether the point is inside the tube circle.
@inline function _torus_integrand(uP, Δ, C, zP, dz, R, r, ψa, sa, t)
    u = uP + t * Δ
    z = zP + t * dz
    x = norm(u) - R
    cψ, sψ = _torus_foot(r, x, z)
    ψ = mod(atan(sψ, cψ) - ψa + π, 2π) - π
    return (R * ψ + r * (sψ - sa)) * C / dot(u, u), hypot(x, z) < r
end

# Integrates over t ∈ [t0, t1] by Gauss–Legendre, split where the edge crosses the tube
# circle. At a crossing the projection switches from the nearest point to the horizontal
# one. The switch puts a kink in the integrand, which would cost Gauss–Legendre its
# convergence to rounding level. Crossings between neighboring nodes are located by
# bisection.
function _torus_edge_piece(uP, Δ, C, zP, dz, R, r, ψa, sa, t0, t1)
    f(t) = _torus_integrand(uP, Δ, C, zP, dz, R, r, ψa, sa, t)
    x, w = _GL16.x, _GL16.w
    acc = 0.0
    cls1 = last(f(t0 + (t1 - t0) * x[1]))
    mixed = false
    for i in eachindex(x)
        v, cls = f(t0 + (t1 - t0) * x[i])
        acc += w[i] * v
        mixed |= cls != cls1
    end
    mixed || return (t1 - t0) * acc
    breaks = [t0]
    for i in 1:length(x)-1
        lo, hi = t0 + (t1 - t0) * x[i], t0 + (t1 - t0) * x[i+1]
        clo = last(f(lo))
        clo == last(f(hi)) && continue
        for _ in 1:60
            mid = (lo + hi) / 2
            (mid == lo || mid == hi) && break
            last(f(mid)) == clo ? (lo = mid) : (hi = mid)
        end
        push!(breaks, (lo + hi) / 2)
    end
    push!(breaks, t1)
    acc = 0.0
    for k in 1:length(breaks)-1
        b0, b1 = breaks[k], breaks[k+1]
        sub = 0.0
        for i in eachindex(x)
            sub += w[i] * first(f(b0 + (b1 - b0) * x[i]))
        end
        acc += (b1 - b0) * sub
    end
    return acc
end

"""
    split_mesh(mesh, keep) -> (selected, rest)

Splits a mesh by patch. `keep(centroid)` is called with the area-weighted centroid of
each patch, and the patch goes to `selected` when it returns true and to `rest`
otherwise. Use it to cut a port or a window out of a surface and make it an `aperture`.
"""
function split_mesh(m::Mesh, keep)
    a, b = Mesh(), Mesh()
    for r in m.patches
        poly = Vec3[]
        for k in r
            for vi in m.triangles[k]
                v = m.vertices[vi]
                any(p -> p == v, poly) || push!(poly, v)
            end
        end
        areas = [norm(_tri_cross(triangle_vertices(m, k)...)) / 2 for k in r]
        cen = sum(areas[i] * sum(triangle_vertices(m, k)) / 3 for (i, k) in enumerate(r)) / sum(areas)
        dest = keep(cen) ? a : b
        i0 = length(dest.vertices)
        t0 = length(dest.triangles) + 1
        vmap = Dict{Int,Int}()
        for k in r, vi in m.triangles[k]
            haskey(vmap, vi) || (push!(dest.vertices, m.vertices[vi]); vmap[vi] = length(dest.vertices))
        end
        for k in r
            push!(dest.triangles, map(vi -> vmap[vi], m.triangles[k]))
        end
        push!(dest.patches, t0:length(dest.triangles))
    end
    return a, b
end

"""
    transform(mesh, R, t)

Applies the rigid map `x ↦ R*x + t` to every vertex. `R` must be a proper rotation for
the patch orientations to be preserved.
"""
transform(m::Mesh, R::AbstractMatrix, t = zero(Vec3)) =
    Mesh([Vec3(R * v + t) for v in m.vertices], copy(m.triangles), copy(m.patches))

"""
    translate(mesh_or_component, t)

Returns a copy shifted by the vector `t` (mm).
"""
translate(m::Mesh, t) = Mesh([v + Vec3(t) for v in m.vertices], copy(m.triangles), copy(m.patches))

"""
    rotation_to(n)

Returns the proper rotation matrix that takes `ẑ` to the unit vector along `n` by the
minimal rotation.
"""
function rotation_to(n)
    z = normalize(Vec3(n))
    c = z[3]
    if c > 1 - 1e-14
        return SMatrix{3,3,Float64}(I)
    elseif c < -1 + 1e-14
        return SMatrix{3,3,Float64}(1, 0, 0, 0, -1, 0, 0, 0, -1)   # rotation by π about x̂
    end
    a = cross(Vec3(0, 0, 1), z)          # axis * sin(angle)
    s = norm(a)
    k = a / s
    K = SMatrix{3,3,Float64}(0, k[3], -k[2], -k[3], 0, k[1], k[2], -k[1], 0)
    return SMatrix{3,3,Float64}(I) + s * K + (1 - c) * K * K
end

# ---------------------------------------------------------------------------------------
# Primitives are built in local coordinates; place them with `transform` or `translate`.
# ---------------------------------------------------------------------------------------

"""
    conforming_grid(a, b, n, breaks; periodic = false) -> Vector{Float64}

Returns grid lines from `a` to `b` with about `n` cells that include every break line in
`breaks`. Each interval between consecutive breaks is divided uniformly into
`max(1, round(len / h))` cells, with `h = (b - a)/n`. Breaks closer than `1e-9 (b - a)` are
merged. `periodic` is for an azimuth with `b - a` one period. The breaks are then wrapped
into the period, and the grid extends from the first break to the first break plus the
period, so the end points `a` and `b` are not grid lines themselves. Without breaks, use
`range(a, b; length = n + 1)`.
"""
function conforming_grid(a, b, n::Integer, breaks; periodic::Bool = false)
    P = b - a
    tol = 1e-9 * abs(P)
    xs = periodic ? sort!([a + mod(x - a, P) for x in breaks]) :
         [a; sort!([Float64(x) for x in breaks if a + tol < x < b - tol]); b]
    keep = Float64[]
    for x in xs
        (isempty(keep) || x - keep[end] > tol) && push!(keep, x)
    end
    if periodic
        isempty(keep) && throw(ArgumentError("periodic grid needs at least one break"))
        keep[end] - keep[1] > P - tol && length(keep) > 1 && pop!(keep)
        push!(keep, keep[1] + P)
    end
    h = P / n
    grid = [keep[1]]
    for i in 1:length(keep)-1
        k = max(1, round(Int, (keep[i+1] - keep[i]) / h))
        append!(grid, range(keep[i], keep[i+1]; length = k + 1)[2:end])
    end
    return grid
end

"""
    cylinder_mesh(radius, z0, z1; nz, nphi, inward = true, phi0 = 0.0, zbreaks = (),
                  phibreaks = ())

Builds an open cylinder of the given radius about the z axis spanning `z0 ≤ z ≤ z1`,
with `nz` axial by `nphi` azimuthal quadrilateral patches. Front faces point toward the
axis when `inward`, as for a hohlraum wall, and away from it otherwise.

`zbreaks` (mm) and `phibreaks` (radians) are lines that the mesh must contain, such as
the edges of a window cut into the wall. The grid is then uniform between consecutive
breaks with about the default spacing (`conforming_grid`), so whole patches tile a
region bounded by break lines exactly. Without breaks the grid is uniform from `phi0`.
"""
function cylinder_mesh(radius, z0, z1; nz::Integer, nphi::Integer, inward::Bool = true,
                       phi0 = 0.0, zbreaks = (), phibreaks = ())
    m = Mesh()
    zs = isempty(zbreaks) ? range(z0, z1; length = nz + 1) : conforming_grid(z0, z1, nz, zbreaks)
    φs = isempty(phibreaks) ? range(phi0, phi0 + 2π; length = nphi + 1) :
         conforming_grid(phi0, phi0 + 2π, nphi, phibreaks; periodic = true)
    nz, nphi = length(zs) - 1, length(φs) - 1
    for i in 1:nz, j in 1:nphi
        a, b = φs[j], φs[j+1]
        p = [Vec3(radius * cos(a), radius * sin(a), zs[i]),
             Vec3(radius * cos(b), radius * sin(b), zs[i]),
             Vec3(radius * cos(b), radius * sin(b), zs[i+1]),
             Vec3(radius * cos(a), radius * sin(a), zs[i+1])]
        φc = (a + b) / 2
        out = Vec3(cos(φc), sin(φc), 0)
        add_polygon!(m, p; facing = inward ? -out : out)
    end
    return m
end

"""
    annulus_mesh(r_inner, r_outer, z; nr, nphi, up = true, phi0 = 0.0, phibreaks = ())

Builds a flat annulus in the plane `z`, or a disk when `r_inner == 0`, with `nr` radial by
`nphi` azimuthal patches. Front faces point along `+ẑ` when `up`, `-ẑ` otherwise. Rings are
equal-area so patch areas are uniform. `phibreaks` (radians) are azimuths the mesh must
contain, as in `cylinder_mesh`; giving both the same `nphi` and `phibreaks` makes an end cap
meet a wall at shared vertices.
"""
function annulus_mesh(r_inner, r_outer, z; nr::Integer, nphi::Integer, up::Bool = true,
                      phi0 = 0.0, phibreaks = ())
    m = Mesh()
    rs = sqrt.(range(r_inner^2, r_outer^2; length = nr + 1))
    φs = isempty(phibreaks) ? range(phi0, phi0 + 2π; length = nphi + 1) :
         conforming_grid(phi0, phi0 + 2π, nphi, phibreaks; periodic = true)
    nphi = length(φs) - 1
    facing = up ? Vec3(0, 0, 1) : Vec3(0, 0, -1)
    for i in 1:nr, j in 1:nphi
        a, b = φs[j], φs[j+1]
        r0, r1 = rs[i], rs[i+1]
        if r0 == 0
            p = [Vec3(0, 0, z), Vec3(r1 * cos(a), r1 * sin(a), z), Vec3(r1 * cos(b), r1 * sin(b), z)]
        else
            p = [Vec3(r0 * cos(a), r0 * sin(a), z), Vec3(r1 * cos(a), r1 * sin(a), z),
                 Vec3(r1 * cos(b), r1 * sin(b), z), Vec3(r0 * cos(b), r0 * sin(b), z)]
        end
        add_polygon!(m, p; facing)
    end
    return m
end

"""
    disk_mesh(radius, z; nr, nphi, up = true)

Builds a flat disk in the plane `z`, an `annulus_mesh` with zero inner radius.
"""
disk_mesh(radius, z; nr::Integer, nphi::Integer, up::Bool = true) =
    annulus_mesh(0.0, radius, z; nr, nphi, up)

"""
    sphere_mesh(radius; ntheta, nphi, outward = true)

Builds a latitude–longitude sphere about the origin with `ntheta` bands uniform in
`cos θ`, so all `ntheta × nphi` patches have equal area. Polar bands are triangle fans.
"""
function sphere_mesh(radius; ntheta::Integer, nphi::Integer, outward::Bool = true)
    m = Mesh()
    μs = range(1.0, -1.0; length = ntheta + 1)
    φs = range(0.0, 2π; length = nphi + 1)
    pt(μ, φ) = (s = sqrt(max(0.0, 1 - μ^2)); radius * Vec3(s * cos(φ), s * sin(φ), μ))
    for i in 1:ntheta, j in 1:nphi
        μ0, μ1 = μs[i], μs[i+1]
        a, b = φs[j], φs[j+1]
        if i == 1
            p = [Vec3(0, 0, radius), pt(μ1, a), pt(μ1, b)]
        elseif i == ntheta
            p = [pt(μ0, a), Vec3(0, 0, -radius), pt(μ0, b)]
        else
            p = [pt(μ0, a), pt(μ1, a), pt(μ1, b), pt(μ0, b)]
        end
        c = sum(p) / length(p)
        add_polygon!(m, p; facing = outward ? c : -c)
    end
    return m
end

"""
    cubed_sphere_mesh(radius; n, outward = true)

Builds a sphere about the origin from the six faces of a cube, each divided into `n × n`
equiangular cells projected onto the sphere. The four corners of a projected cell are
not coplanar. Each cell is therefore split into two triangular patches along alternating
diagonals, chosen to keep the faceted sphere convex, for `12n²` planar patches. The mesh
has no poles and no slivers, which suits radiation problems where every part of the
sphere matters equally.
"""
function cubed_sphere_mesh(radius; n::Integer, outward::Bool = true)
    m = Mesh()
    axes = (Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1))
    ts = range(-π / 4, π / 4; length = n + 1)
    for (fi, ax) in enumerate(axes), sgn in (1, -1)
        u = axes[mod1(fi + 1, 3)]
        v = axes[mod1(fi + 2, 3)]
        pt(a, b) = radius * normalize(sgn * ax + tan(a) * u + tan(b) * v)
        for i in 1:n, j in 1:n
            p = [pt(ts[i], ts[j]), pt(ts[i+1], ts[j]), pt(ts[i+1], ts[j+1]), pt(ts[i], ts[j+1])]
            # Take the diagonal that keeps the faceted sphere convex. The fourth corner must
            # lie inside the plane of the other three, because a reflex fold would hide
            # neighboring elements from each other.
            n123 = cross(p[2] - p[1], p[3] - p[1])
            n123 = dot(n123, p[1]) < 0 ? -n123 : n123          # outward
            halves = dot(n123, p[4] - p[1]) <= 0 ? ((p[1], p[2], p[3]), (p[1], p[3], p[4])) :
                     ((p[1], p[2], p[4]), (p[2], p[3], p[4]))
            for tri in halves
                c = sum(tri) / 3
                add_polygon!(m, collect(tri); facing = outward ? c : -c)
            end
        end
    end
    return m
end

"""
    icosphere_mesh(radius; n, outward = true)

Builds a geodesic sphere about the origin. Each face of an icosahedron is divided into
`n × n` triangles whose vertices are projected onto the sphere, for `20n²` triangular
patches. The icosahedron is in its standard frame, with vertices at the cyclic
permutations of (0, ±1, ±ϕ), and is mirror-symmetric in x, y and z. The triangles are
nearly equilateral everywhere, with edges within about ±20% of `1.107 radius / n`, and
the mesh has no poles and no slivers. On a capsule the high Legendre moments (P4) of the
flux converge with far fewer elements than on `sphere_mesh`, whose polar cells are
slivers. `icosphere_frequency` gives `n` for an
element size.
"""
function icosphere_mesh(radius; n::Integer, outward::Bool = true)
    n >= 1 || throw(ArgumentError("the subdivision n must be at least 1"))
    m = Mesh()
    ϕ = (1 + sqrt(5)) / 2
    # Place the icosahedron in its standard frame, with vertices at the cyclic
    # permutations of (0, ±1, ±ϕ). Its 2-fold axes and mirror planes are x, y and z, the
    # symmetries of a hohlraum with a NIF-like (D2) beam set. The mesh therefore does not
    # add modes beyond those of the beams and wall. With a vertex on z, which is a 5-fold
    # axis, odd-m modes leaked at ~6e-5 of the drive.
    V = [normalize(v) for s1 in (-1, 1) for s2 in (-1, 1)
         for v in (Vec3(0, s1, s2 * ϕ), Vec3(s1, s2 * ϕ, 0), Vec3(s2 * ϕ, 0, s1))]
    e = minimum(norm(V[i] - V[j]) for i in 1:12 for j in i+1:12)
    faces = [(i, j, k) for i in 1:12 for j in i+1:12 for k in j+1:12
             if all(d -> isapprox(d, e; rtol = 1e-9), (norm(V[i] - V[j]), norm(V[j] - V[k]), norm(V[i] - V[k])))]
    # Build a point from integer barycentric weights, summed in vertex-index order so that
    # points on shared edges come out bitwise identical from either face.
    pt(i, j, k, a, b, c) = radius * normalize(a * V[i] + b * V[j] + c * V[k])
    for (i, j, k) in faces
        P(a, b) = pt(i, j, k, n - a - b, a, b)
        for a in 0:n-1, b in 0:n-1-a
            tris = [(P(a, b), P(a + 1, b), P(a, b + 1))]
            a + b < n - 1 && push!(tris, (P(a + 1, b), P(a + 1, b + 1), P(a, b + 1)))
            for t in tris
                c = sum(t) / 3
                add_polygon!(m, collect(t); facing = outward ? c : -c)
            end
        end
    end
    return m
end

"""
    icosphere_frequency(radius, h) -> n

Returns the subdivision of `icosphere_mesh` for triangles with sides of about `h` (mm).
An icosahedron edge spans 1.107 radians of arc.
"""
function icosphere_frequency(r, h)
    h > 0 || throw(ArgumentError("the element size must be positive"))
    return max(1, ceil(Int, atan(2.0) * r / h))
end

"""
    rectangle_mesh(width, height; nu, nv)

Builds the rectangle `|x| ≤ width/2, |y| ≤ height/2` in the plane `z = 0`, with its front
face along `+ẑ`, divided into `nu × nv` patches.
"""
function rectangle_mesh(width, height; nu::Integer, nv::Integer)
    m = Mesh()
    xs = range(-width / 2, width / 2; length = nu + 1)
    ys = range(-height / 2, height / 2; length = nv + 1)
    for i in 1:nu, j in 1:nv
        add_polygon!(m, [Vec3(xs[i], ys[j], 0), Vec3(xs[i+1], ys[j], 0),
                         Vec3(xs[i+1], ys[j+1], 0), Vec3(xs[i], ys[j+1], 0)];
                     facing = Vec3(0, 0, 1))
    end
    return m
end

"""
    frustum_mesh(r0, z0, r1, z1; nz, nphi, inward = true)

Builds a conical frustum about the z axis from radius `r0` at `z0` to `r1` at `z1`, or a
cone when one radius is zero. It is used for shields, cones and tapered hohlraum sections.
"""
function frustum_mesh(r0, z0, r1, z1; nz::Integer, nphi::Integer, inward::Bool = true)
    m = Mesh()
    ts = range(0.0, 1.0; length = nz + 1)
    φs = range(0.0, 2π; length = nphi + 1)
    rz(t) = (r0 + (r1 - r0) * t, z0 + (z1 - z0) * t)
    for i in 1:nz, j in 1:nphi
        (ra, za), (rb, zb) = rz(ts[i]), rz(ts[i+1])
        a, b = φs[j], φs[j+1]
        p = Vec3[]
        ra > 0 && append!(p, [Vec3(ra * cos(a), ra * sin(a), za), Vec3(ra * cos(b), ra * sin(b), za)])
        ra == 0 && push!(p, Vec3(0, 0, za))
        rb > 0 && append!(p, [Vec3(rb * cos(b), rb * sin(b), zb), Vec3(rb * cos(a), rb * sin(a), zb)])
        rb == 0 && push!(p, Vec3(0, 0, zb))
        φc = (a + b) / 2
        # Compute the outward surface normal of the cone wall.
        out = normalize(Vec3((z1 - z0) * cos(φc), (z1 - z0) * sin(φc), -(r1 - r0)))
        add_polygon!(m, p; facing = inward ? -out : out)
    end
    return m
end
