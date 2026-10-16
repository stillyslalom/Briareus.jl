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
    camera_from_port(theta, phi; kwargs...)

Camera at the port (θ, φ) in degrees looking at TCC — e.g. an x-ray imager's line of
sight.
"""
camera_from_port(θ, φ; kwargs...) = OrthoCamera(-port_vector(θ, φ); kwargs...)

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
    line_of_sight_power(target, B, direction; elems = all) -> Real

X-ray power per steradian (TW/sr) radiated toward a distant detector in `direction`
(unit vector from the target toward the detector) by the given elements with
radiosity `B` (TW/mm²), counting only unobstructed elements: `Σ Bᵢ/π · Aᵢ cosθᵢ`.
"""
function line_of_sight_power(t::Target, B::AbstractVector, direction;
                             elems = 1:nelements(t))
    d = normalize(Vec3(direction))
    ε = _ray_eps(t)
    s = zero(eltype(B))
    for e in elems
        c = dot(t.elem_normal[e], d)
        c > 0 || continue
        any_hit(t.bvh, t.elem_centroid[e] + ε * t.elem_normal[e], d) && continue
        s += B[e] / π * t.elem_area[e] * c
    end
    return s
end
