# Synthetic images and line-of-sight radiance.

"""
    OrthoCamera(direction; up = ẑ, center = 0, width, height = width, nx = 256, ny = nx)

Orthographic camera looking along `direction` (from the camera toward the scene) at
`center`, imaging a `width × height` mm field.
"""
struct OrthoCamera
    dir::Vec3
    right::Vec3
    up::Vec3
    center::Vec3
    width::Float64
    height::Float64
    nx::Int
    ny::Int
end

function OrthoCamera(direction; up = Vec3(0, 0, 1), center = zero(Vec3), width,
                     height = width, nx::Integer = 256, ny::Integer = nx)
    d = normalize(Vec3(direction))
    u = Vec3(up)
    if abs(dot(u, d)) > 0.999
        u = abs(d[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)
    end
    r = normalize(cross(d, u))
    u = cross(r, d)
    return OrthoCamera(d, r, u, Vec3(center), Float64(width), Float64(height), nx, ny)
end

"""
    render(target, values, camera; background = NaN) -> Matrix

Image of a per-element field: each pixel shows `values[e]` for the element seen first
along its line of sight (apertures are transparent), or `background`. Row 1 is the top
of the image. Pass `radiance(sol)` for an x-ray image or `laser_flux(dep, target)` for a
deposition map.
"""
function render(t::Target, values::AbstractVector, cam::OrthoCamera; background = NaN)
    img = fill(Float64(background), cam.ny, cam.nx)
    far = 4 * (t.scale + norm(cam.center)) + 1
    Threads.@threads for j in 1:cam.nx
        x = (j - 0.5) / cam.nx * cam.width - cam.width / 2
        for i in 1:cam.ny
            y = cam.height / 2 - (i - 0.5) / cam.ny * cam.height
            o = cam.center + x * cam.right + y * cam.up - far * cam.dir
            _, tri = closest_hit(t.bvh, o, cam.dir)
            tri == 0 && continue
            e = dot(cam.dir, t.tri_normal[tri]) < 0 ? t.tri_front[tri] : t.tri_back[tri]
            e > 0 && (img[i, j] = ForwardDiff.value(values[e]))
        end
    end
    return img
end

"""
    radiance(sol) -> Vector

Lambertian x-ray radiance `B/π` per element (TW mm⁻² sr⁻¹).
"""
radiance(s::RadiationSolution) = s.B ./ π

"""
    line_of_sight_power(target, B, direction; elems = all, nsub = 4, refine = 8) -> Real

X-ray power per steradian (TW/sr) radiated toward a distant detector in `direction`
(unit vector from the target toward the detector) by the given elements with
radiosity `B` (TW/mm²), each counted with its unobstructed fraction `fᵢ`:
`Σ Bᵢ/π · fᵢ Aᵢ cosθᵢ`.

`fᵢ` is the area-weighted fraction of points on the element from which a ray toward the
detector meets no opaque surface (apertures are transparent): the centroids of the
`nsub²` equal sub-triangles of each of its triangles, and of `(refine·nsub)²` on an
element whose first points disagree (an element partly hidden behind an edge). The
error is first order in the sub-triangle size, on partly hidden elements only; an
occluder that misses every first point of an element is not seen.
"""
function line_of_sight_power(t::Target, B::AbstractVector, direction;
                             elems = 1:nelements(t), nsub::Integer = 4, refine::Integer = 8)
    nsub >= 1 && refine >= 1 || throw(ArgumentError("nsub and refine must be positive"))
    d = normalize(Vec3(direction))
    s = zero(eltype(B))
    for e in elems
        c = dot(t.elem_normal[e], d)
        c > 0 || continue
        f = _visible_fraction(t, e, d, nsub)
        0 < f < 1 && refine > 1 && (f = _visible_fraction(t, e, d, refine * nsub))
        f > 0 && (s += B[e] / π * t.elem_area[e] * c * f)
    end
    return s
end

# The area-weighted fraction of the centroids of the n² equal sub-triangles of element
# `e`'s triangles that see along `d` without meeting an opaque surface.
function _visible_fraction(t::Target, e::Integer, d::Vec3, n::Integer)
    off = _ray_eps(t) * t.elem_normal[e]
    seen = total = 0.0
    for k in t.elem_tris[e]
        v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
        m = 0
        for i in 0:n-1, j in 0:n-1-i
            # the upward sub-triangle at (i, j), and the downward one beside it
            m += !any_hit(t.bvh, v0 + ((i + 1 / 3) / n) * e1 + ((j + 1 / 3) / n) * e2 + off, d)
            i + j <= n - 2 &&
                (m += !any_hit(t.bvh, v0 + ((i + 2 / 3) / n) * e1 + ((j + 2 / 3) / n) * e2 + off, d))
        end
        seen += t.tri_area[k] * m / n^2
        total += t.tri_area[k]
    end
    return seen / total
end
