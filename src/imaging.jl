# This file renders synthetic images and computes line-of-sight radiance.

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
    render(target, values, camera; background = NaN, transmission = nothing) -> Matrix

Image of a per-element field. Each pixel shows `values[e]` for the element first met
along its line of sight, or `background` where there is none; apertures are
transparent. Row 1 is the top of the image. Pass `radiance(sol)` for an x-ray image or
`laser_flux(dep, target)` for a deposition map.

`transmission` is a per-element open fraction, the `transmission` of a
`RadiationSolution`. With it, a pixel seen through the closing lip of an entrance hole
(an aperture with a `GapClosure` or `LipClosure` albedo) is weighted by the open
fraction of the lip element that its line of sight crosses, as in [`DanteView`](@ref).
"""
function render(t::Target, values::AbstractVector, cam::OrthoCamera; background = NaN,
                transmission::Union{Nothing,AbstractVector} = nothing)
    elems, lips = _first_elements(t, cam, transmission !== nothing)
    img = fill(Float64(background), cam.ny, cam.nx)
    for j in 1:cam.nx, i in 1:cam.ny
        e = elems[i, j]
        e > 0 || continue
        v = Float64(ForwardDiff.value(values[e]))
        l = lips[i, j]
        img[i, j] = l > 0 ? v * Float64(transmission[l]) : v
    end
    return img
end

# Returns the first element met along the line of sight of each pixel of `cam` (0 for
# none). With `with_lips` it also returns the closing-lip element that the line of sight
# crosses before that element (0 for none); otherwise that matrix is all zero.
function _first_elements(t::Target, cam::OrthoCamera, with_lips::Bool)
    elems, lips = zeros(Int32, cam.ny, cam.nx), zeros(Int32, cam.ny, cam.nx)
    far = 4 * (t.scale + norm(cam.center)) + 1
    ltris = with_lips ? _lip_triangles(t) : Int[]
    Threads.@threads for j in 1:cam.nx
        x = (j - 0.5) / cam.nx * cam.width - cam.width / 2
        for i in 1:cam.ny
            y = cam.height / 2 - (i - 0.5) / cam.ny * cam.height
            o = cam.center + x * cam.right + y * cam.up - far * cam.dir
            th, tri = closest_hit(t.bvh, o, cam.dir)
            tri == 0 && continue
            e = dot(cam.dir, t.tri_normal[tri]) < 0 ? t.tri_front[tri] : t.tri_back[tri]
            e > 0 || continue
            elems[i, j] = e
            isempty(ltris) || (lips[i, j] = _first_lip(t, ltris, o, cam.dir, th))
        end
    end
    return elems, lips
end

# Returns the triangles of the transmitting apertures of a target, which are the closing
# lips of its entrance holes.
_lip_triangles(t::Target) =
    [k for (ci, comp) in enumerate(t.components)
     if comp.kind === :aperture && _transmits(comp.material.albedo)
     for e in t.comp_elems[ci] for k in t.elem_tris[e]]

# Returns the element of the lip triangle among `ltris` that the ray from `o` along `d`
# crosses first within the distance `tmax`, or 0 when it crosses none.
function _first_lip(t::Target, ltris::Vector{Int}, o::Vec3, d::Vec3, tmax)
    l, tl = 0, tmax
    for k in ltris
        s = _intersect(o, d, t.tri_v0[k], t.tri_e1[k], t.tri_e2[k], tl)
        s < tl && ((tl, l) = (s, Int(max(t.tri_front[k], t.tri_back[k]))))
    end
    return l
end

"""
    radiance(sol) -> Vector

Lambertian x-ray radiance `B/π` per element (TW mm⁻² sr⁻¹).
"""
radiance(s::RadiationSolution) = s.B ./ π

"""
    line_of_sight_power(target, B, direction; elems = all, nsub = 4, refine = 8,
                        transmission = nothing) -> Real

X-ray power per steradian (TW/sr) radiated toward a distant detector in `direction`
(unit vector from the target toward the detector) by the given elements with
radiosity `B` (TW/mm²), each counted with its unobstructed fraction `fᵢ`:
`Σ Bᵢ/π · fᵢ Aᵢ cosθᵢ`.

`fᵢ` is the area-weighted fraction of points on the element from which a ray toward the
detector does not meet an opaque surface; apertures are transparent. The points are the
centroids of the `nsub²` equal sub-triangles of each of its triangles. On an element
whose first points give different results, such as an element partly hidden behind an
edge, `(refine·nsub)²` points are used. The error is first order in the sub-triangle
size, on partly hidden elements only. An occluder that misses every first point of an
element is not detected.

`transmission` is a per-element open fraction, the `transmission` of a
`RadiationSolution`. With it, a point seen through the closing lip of an entrance hole
counts with the open fraction of the lip element that its ray crosses first, as in
[`DanteView`](@ref) and [`render`](@ref). Without it, a closing lip is transparent.
"""
function line_of_sight_power(t::Target, B::AbstractVector, direction;
                             elems = 1:nelements(t), nsub::Integer = 4, refine::Integer = 8,
                             transmission::Union{Nothing,AbstractVector} = nothing)
    nsub >= 1 && refine >= 1 || throw(ArgumentError("nsub and refine must be positive"))
    d = normalize(Vec3(direction))
    ltris = transmission === nothing ? Int[] : _lip_triangles(t)
    s = zero(eltype(B))
    for e in elems
        c = dot(t.elem_normal[e], d)
        c > 0 || continue
        f = _visible_fraction(t, e, d, nsub, ltris, transmission)
        0 < f < 1 && refine > 1 && (f = _visible_fraction(t, e, d, refine * nsub, ltris, transmission))
        f > 0 && (s += B[e] / π * t.elem_area[e] * c * f)
    end
    return s
end

# Returns the area-weighted fraction of the centroids of the n² equal sub-triangles of
# the triangles of element `e` from which a ray along `d` does not meet an opaque
# surface. With lip triangles `ltris`, a point whose ray crosses a closing lip counts with
# the open fraction `τ` of the lip element crossed first.
function _visible_fraction(t::Target, e::Integer, d::Vec3, n::Integer, ltris::Vector{Int} = Int[], τ = nothing)
    off = _ray_eps(t) * t.elem_normal[e]
    seen = total = 0.0
    function point(p)
        any_hit(t.bvh, p, d) && return 0.0
        isempty(ltris) && return 1.0
        l = _first_lip(t, ltris, p, d, Inf)
        return l > 0 ? Float64(τ[l]) : 1.0
    end
    for k in t.elem_tris[e]
        v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
        m = 0.0
        for i in 0:n-1, j in 0:n-1-i
            # The upward sub-triangle at (i, j) is tested, then the downward one beside it.
            m += point(v0 + ((i + 1 / 3) / n) * e1 + ((j + 1 / 3) / n) * e2 + off)
            i + j <= n - 2 && (m += point(v0 + ((i + 2 / 3) / n) * e1 + ((j + 2 / 3) / n) * e2 + off))
        end
        seen += t.tri_area[k] * m / n^2
        total += t.tri_area[k]
    end
    return seen / total
end
