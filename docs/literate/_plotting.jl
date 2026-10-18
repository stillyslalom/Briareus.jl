# Figure helpers shared by the tutorials. Not a tutorial page: the tutorials `include`
# it. Plotting stays out of `src/`, so these use a few of Briareus's internals (the BVH
# and the triangle arrays of a `Target`).

using Briareus
using CairoMakie
using LinearAlgebra

set_theme!()   # every page starts from Makie's default theme, whatever ran before it

"""
    VertexField(target, values, component)

Values of the per-element field `values` (on every element of the target, or on those
of `component` in order) on the elements of `component`, averaged with
area weights onto the mesh vertices so that an image can interpolate linearly across
each triangle (Gouraud shading). Vertices are not shared across components, so the
shading does not blur a field across the edge between, say, a wall and an aperture.
"""
struct VertexField
    target::Target
    elems::BitVector
    corner::Dict{Int,NTuple{3,Float64}}   # triangle => values at v0, v0 + e1, v0 + e2
end

function VertexField(t::Target, values::AbstractVector, component::Symbol)
    elems = falses(nelements(t))
    elems[elements(t, component)] .= true
    if length(values) != nelements(t)      # values given on the component only
        length(values) == count(elems) || throw(DimensionMismatch("values"))
        values = setindex!(fill(NaN, nelements(t)), values, elements(t, component))
    end
    key(v, e) = (round.(Tuple(v) ./ t.scale; digits = 8) .+ 0.0, t.elem_side[e])   # +0.0: no -0.0
    corners(k) = (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
    acc = Dict{Any,Tuple{Float64,Float64}}()
    for e in findall(elems), k in t.elem_tris[e], v in corners(k)
        s, w = get(acc, key(v, e), (0.0, 0.0))
        acc[key(v, e)] = (s + t.tri_area[k] * values[e], w + t.tri_area[k])
    end
    vertex(v, e) = (sw = acc[key(v, e)]; sw[1] / sw[2])
    corner = Dict{Int,NTuple{3,Float64}}()
    for e in findall(elems), k in t.elem_tris[e]
        corner[k] = map(v -> vertex(v, e), corners(k))
    end
    return VertexField(t, elems, corner)
end

"Value of the shaded field where the ray from `o` along `d` first meets the target."
function sample(f::VertexField, o, d)
    t = f.target
    o, d = Briareus.Vec3(o), Briareus.Vec3(d)
    s, k = Briareus.closest_hit(t.bvh, o, d)
    k == 0 && return NaN
    e = dot(d, t.tri_normal[k]) < 0 ? t.tri_front[k] : t.tri_back[k]
    (e > 0 && f.elems[e]) || return NaN
    r = o + s * d - t.tri_v0[k]
    e1, e2 = t.tri_e1[k], t.tri_e2[k]
    a, b, c = dot(e1, e1), dot(e1, e2), dot(e2, e2)
    p, q = dot(r, e1), dot(r, e2)
    u, v = (c * p - b * q) / (a * c - b^2), (a * q - b * p) / (a * c - b^2)
    f0, f1, f2 = f.corner[k]
    return (1 - u - v) * f0 + u * f1 + v * f2
end

"""
    render_smoothed(target, values, component, camera; σ) -> Matrix

Like [`render`](@ref), but each pixel showing `component` takes a Gaussian-weighted
(width `σ`, in mm) and area-weighted average of the values of that component's
elements around the point it sees, instead of the one value of the element it hits.
Other components are left blank. Row 1 is the top of the image.
"""
function render_smoothed(t::Target, values, component::Symbol, cam::OrthoCamera; σ)
    els = elements(t, component)
    inside = falses(nelements(t))
    inside[els] .= true
    c, v, a = t.elem_centroid[els], values[els], t.elem_area[els]
    img = fill(NaN, cam.ny, cam.nx)
    far = 4 * (t.scale + norm(cam.center)) + 1
    Threads.@threads for j in 1:cam.nx
        x = (j - 0.5) / cam.nx * cam.width - cam.width / 2
        for i in 1:cam.ny
            y = cam.height / 2 - (i - 0.5) / cam.ny * cam.height
            o = cam.center + x * cam.right + y * cam.up - far * cam.dir
            s, k = Briareus.closest_hit(t.bvh, o, cam.dir)
            k == 0 && continue
            e = dot(cam.dir, t.tri_normal[k]) < 0 ? t.tri_front[k] : t.tri_back[k]
            (e > 0 && inside[e]) || continue
            p = o + s * cam.dir
            num = den = 0.0
            for q in eachindex(c)
                r2 = sum(abs2, c[q] - p)
                r2 < 9σ^2 || continue
                w = a[q] * exp(-r2 / 2σ^2)
                num += w * v[q]
                den += w
            end
            img[i, j] = num / den
        end
    end
    return img
end

"Pixel centers of a camera's image, in mm across and up, for `heatmap!`."
camera_axes(cam::OrthoCamera) = (range(-cam.width / 2, cam.width / 2; length = cam.nx + 1)[1:end-1] .+ cam.width / 2cam.nx,
                                 range(-cam.height / 2, cam.height / 2; length = cam.ny + 1)[1:end-1] .+ cam.height / 2cam.ny)

"`heatmap!` an image from [`render`](@ref) (row 1 at the top) in camera coordinates."
camera_image!(ax, img, cam::OrthoCamera; kwargs...) =
    heatmap!(ax, camera_axes(cam)..., permutedims(reverse(img; dims = 1)); kwargs...)

# ## Maps of a sphere
#
# The Hammer projection is equal-area: a patch of the map covers a solid angle in
# proportion to its area, so the eye weighs hot and cold regions as the capsule does.
# The map is centered on azimuth `φ0`, with φ increasing to the right as seen from
# outside the sphere and θ = 0 at the top.

function hammer(θ, φ; φ0 = 0.0)
    lat, lon = deg2rad(90 - θ), deg2rad(mod(φ - φ0 + 180, 360) - 180)
    s = sqrt(1 + cos(lat) * cos(lon / 2))
    return (2sqrt(2) * cos(lat) * sin(lon / 2) / s, sqrt(2) * sin(lat) / s)
end

function inverse_hammer(x, y; φ0 = 0.0)
    (x / 2sqrt(2))^2 + (y / sqrt(2))^2 > 1 && return (NaN, NaN)
    z = sqrt(1 - (x / 4)^2 - (y / 2)^2)
    return (90 - rad2deg(asin(clamp(z * y, -1, 1))), φ0 + rad2deg(2atan(z * x, 2(2z^2 - 1))))
end

"""
    sphere_map(target, values, component; nx = 640, φ0 = 0) -> (xs, ys, image)

Equal-area (Hammer) map of a per-element field on a spherical `component` seen from
outside, shaded from its vertex-averaged values.
"""
function sphere_map(t::Target, values, component::Symbol; nx = 640, φ0 = 0.0)
    f = VertexField(t, values, component)
    els = elements(t, component)
    w = t.elem_area[els]
    center = sum(w .* t.elem_centroid[els]) / sum(w)
    R = maximum(norm(t.tri_v0[k] - center) for e in els for k in t.elem_tris[e])
    xs = range(-2sqrt(2), 2sqrt(2); length = nx)
    ys = range(-sqrt(2), sqrt(2); length = nx ÷ 2)
    img = fill(NaN, length(xs), length(ys))
    Threads.@threads for i in eachindex(xs)
        for j in eachindex(ys)
            θ, φ = inverse_hammer(xs[i], ys[j]; φ0)
            isnan(θ) && continue
            n = [sind(θ) * cosd(φ), sind(θ) * sind(φ), cosd(θ)]
            img[i, j] = sample(f, center + 1.5R * n, -n)
        end
    end
    return xs, ys, img
end

"""
    sphere_map(f; nx = 640, φ0 = 0) -> (xs, ys, image)

Equal-area (Hammer) map of a function `f(θ, φ)` of the direction, in degrees.
"""
function sphere_map(f::Function; nx = 640, φ0 = 0.0)
    xs = range(-2sqrt(2), 2sqrt(2); length = nx)
    ys = range(-sqrt(2), sqrt(2); length = nx ÷ 2)
    img = fill(NaN, length(xs), length(ys))
    for i in eachindex(xs), j in eachindex(ys)
        θ, φ = inverse_hammer(xs[i], ys[j]; φ0)
        isnan(θ) || (img[i, j] = f(θ, φ))
    end
    return xs, ys, img
end

"Radius (mm) of the largest circle about the center of a camera image with no blank pixel."
function inscribed_radius(img, cam::OrthoCamera)
    xs, ys = camera_axes(cam)
    return minimum(hypot(xs[j], ys[end+1-i]) for i in axes(img, 1), j in axes(img, 2)
                   if isnan(img[i, j]))
end

"Blank the pixels of a camera image outside the circle of `radius` (mm) about its center."
function crop_to_circle!(img, cam::OrthoCamera, radius)
    xs, ys = camera_axes(cam)
    for i in axes(img, 1), j in axes(img, 2)
        hypot(xs[j], ys[end+1-i]) >= radius && (img[i, j] = NaN)
    end
    return img
end

"An axis for a Hammer map, with a graticule every 30° in θ and 45° in φ."
function sphere_axis(pos; title = "", φ0 = 0.0)
    ax = Axis(pos; title, aspect = DataAspect(), limits = ((-3.3, 2.9), (-1.48, 1.48)))
    hidedecorations!(ax)
    hidespines!(ax)
    style = (color = (:black, 0.3), linewidth = 0.6)
    for θ in 30:30:150
        lines!(ax, [Point2f(hammer(θ, φ; φ0)) for φ in range(φ0 - 179.999, φ0 + 179.999; length = 181)]; style...)
    end
    for φ in φ0 .+ (-135:45:135)
        lines!(ax, [Point2f(hammer(θ, φ; φ0)) for θ in 0:1:180]; style...)
    end
    lines!(ax, [Point2f(hammer(θ, φ0 + s * 179.999; φ0)) for s in (-1, 1) for θ in (s < 0 ? (0:180) : (180:-1:0))];
           color = :black, linewidth = 0.8)
    for φ in φ0 .+ (-135:45:135)
        text!(ax, Point2f(hammer(90, φ; φ0)) .- Point2f(0, 0.1); text = "$(Int(mod(φ, 360)))°",
              align = (:center, :top), fontsize = 10, color = (:black, 0.8))
    end
    for θ in 30:30:150
        text!(ax, Point2f(hammer(θ, φ0 - 179.999; φ0)) .- Point2f(0.04, 0); text = "$(θ)°",
              align = (:right, :center), fontsize = 10, color = (:black, 0.75))
    end
    return ax
end

"Draw a Hammer map from [`sphere_map`](@ref) into `ax`, beneath its graticule."
function sphere_heatmap!(ax, (xs, ys, img); kwargs...)
    hm = heatmap!(ax, xs, ys, img; nan_color = :transparent, kwargs...)
    translate!(hm, 0, 0, -10)
    return hm
end

"Mark the directions (θ, φ) in degrees on a Hammer map."
mark_directions!(ax, θφ; φ0 = 0.0, kwargs...) =
    scatter!(ax, [Point2f(hammer(θ, φ; φ0)) for (θ, φ) in θφ]; kwargs...)

# ## Cutaway drawings of meshes
#
# CairoMakie draws polygons in the order given, so a mesh seen from one side is drawn
# by the painter's algorithm: farthest polygon first. Sorting by centroid depth is
# right for the compact elements of a target, and good enough for a picture of long
# CAD slivers.

"""
    patch_outline(mesh, patch) -> Vector{Vec3}

Boundary loop of a mesh patch (an element): the edges used by one of its triangles,
chained in the triangles' winding (the longest loop when there are several).
"""
function patch_outline(m::Briareus.Mesh, p::Integer)
    r = m.patches[p]
    length(r) == 1 && return [m.vertices[i] for i in m.triangles[first(r)]]
    key(i) = round.(Tuple(m.vertices[i]); digits = 9)
    edges(k) = zip(m.triangles[k], circshift(collect(m.triangles[k]), -1))
    uses = Dict{Any,Int}()
    for k in r, (a, b) in edges(k)
        e = minmax(key(a), key(b))
        uses[e] = get(uses, e, 0) + 1
    end
    next = Dict{Any,Any}()
    point = Dict{Any,Briareus.Vec3}()
    for k in r, (a, b) in edges(k)
        uses[minmax(key(a), key(b))] == 1 || continue
        next[key(a)] = key(b)
        point[key(a)] = m.vertices[a]
    end
    best = Briareus.Vec3[]
    seen = Set{Any}()
    for s in keys(next)
        s in seen && continue
        loop, v = Briareus.Vec3[], s
        while !(v in seen) && haskey(next, v)
            push!(seen, v)
            push!(loop, point[v])
            v = next[v]
        end
        length(loop) > length(best) && (best = loop)
    end
    return best
end

"""
    cutaway!(ax, polygons, colors; view = (0, 1, -0.3), strokewidth = 0.3, strokecolor)

Draw 3-D polygons (vectors of points, mm) in a 2-D axis as seen along `view`, farthest
first, in screen coordinates across (perpendicular to `view` and to z) and up.
"""
function cutaway!(ax, polys, colors; view = (0.0, 1.0, -0.3), strokewidth = 0.3,
                  strokecolor = (:black, 0.45))
    f = normalize(Briareus.Vec3(view...))
    r = normalize(cross(f, Briareus.Vec3(0, 0, 1)))
    u = cross(r, f)
    depth = [sum(dot(p, f) for p in P) / length(P) for P in polys]
    order = sortperm(depth; rev = true)
    poly!(ax, [Makie.GeometryBasics.Polygon([Point2f(dot(p, r), dot(p, u)) for p in polys[i]]) for i in order];
          color = colors[order], strokewidth, strokecolor)
    return ax
end
