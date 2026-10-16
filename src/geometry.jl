# Patch-structured triangle meshes and the geometric primitives targets are built from.
#
# A *patch* is a planar convex polygon stored as a fan of triangles; patches become the
# radiosity elements of a `Target`. Triangle winding defines the patch's front side:
# the front normal is (v2 - v1) × (v3 - v1).

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

Append the convex planar polygon `poly` as one patch (fan-triangulated). If `facing` is
a vector, the winding is reversed when needed so the patch's front normal points along
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

Concatenate meshes into one, preserving patch structure and order.
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
    split_mesh(mesh, keep) -> (selected, rest)

Split a mesh by patch: `keep(centroid)` decides, from a patch's area-weighted
centroid, whether the patch goes to `selected` or `rest`. Use it to cut a port or a
window out of a surface and make it an `aperture`.
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

Apply the rigid map `x ↦ R*x + t` to every vertex. `R` must be a proper rotation for
the patch orientations to be preserved.
"""
transform(m::Mesh, R::AbstractMatrix, t = zero(Vec3)) =
    Mesh([Vec3(R * v + t) for v in m.vertices], copy(m.triangles), copy(m.patches))

"""
    translate(mesh_or_component, t)

Copy shifted by the vector `t` (mm).
"""
translate(m::Mesh, t) = Mesh([v + Vec3(t) for v in m.vertices], copy(m.triangles), copy(m.patches))

"""
    rotation_to(n)

Proper rotation matrix taking `ẑ` to the unit vector along `n` (minimal rotation).
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
# Primitives (all in local coordinates; place them with `transform`/`translate`)
# ---------------------------------------------------------------------------------------

"""
    cylinder_mesh(radius, z0, z1; nz, nphi, inward = true, phi0 = 0.0)

Open cylinder of the given radius about the z axis spanning `z0 ≤ z ≤ z1`, with `nz`
axial by `nphi` azimuthal quadrilateral patches. Front faces point toward the axis
when `inward` (a hohlraum wall), away from it otherwise.
"""
function cylinder_mesh(radius, z0, z1; nz::Integer, nphi::Integer, inward::Bool = true,
                       phi0 = 0.0)
    m = Mesh()
    zs = range(z0, z1; length = nz + 1)
    φs = range(phi0, phi0 + 2π; length = nphi + 1)
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
    annulus_mesh(r_inner, r_outer, z; nr, nphi, up = true)

Flat annulus (a disk when `r_inner == 0`) in the plane `z`, with `nr` radial by `nphi`
azimuthal patches. Front faces point along `+ẑ` when `up`, `-ẑ` otherwise. Rings are
equal-area so patch areas are uniform.
"""
function annulus_mesh(r_inner, r_outer, z; nr::Integer, nphi::Integer, up::Bool = true,
                      phi0 = 0.0)
    m = Mesh()
    rs = sqrt.(range(r_inner^2, r_outer^2; length = nr + 1))
    φs = range(phi0, phi0 + 2π; length = nphi + 1)
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

Flat disk in the plane `z`; an `annulus_mesh` with zero inner radius.
"""
disk_mesh(radius, z; nr::Integer, nphi::Integer, up::Bool = true) =
    annulus_mesh(0.0, radius, z; nr, nphi, up)

"""
    sphere_mesh(radius; ntheta, nphi, outward = true)

Latitude–longitude sphere about the origin with `ntheta` bands uniform in `cos θ`, so
all `ntheta × nphi` patches have equal area. Polar bands are triangle fans.
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
    rectangle_mesh(width, height; nu, nv)

Rectangle `|x| ≤ width/2, |y| ≤ height/2` in the plane `z = 0`, front face along `+ẑ`,
divided into `nu × nv` patches.
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

Conical frustum (a cone when one radius is zero) about the z axis from radius `r0` at
`z0` to `r1` at `z1`. Useful for shields, cones and tapered hohlraum sections.
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
        # outward surface normal of the cone wall
        out = normalize(Vec3((z1 - z0) * cos(φc), (z1 - z0) * sin(φc), -(r1 - r0)))
        add_polygon!(m, p; facing = inward ? -out : out)
    end
    return m
end
