# Monte Carlo radiosity: an independent reference for the view factors and the radiosity
# solve, by ray tracing diffuse transport on the facets.

"""
    montecarlo_incident(target, laser_power, receivers; paths = 10_000, batches = 32,
                        reemit = :hit, albedo = albedos(target), seed = 1)
        -> (H, batch)

Element-averaged incident x-ray flux (TW/mm²) on the `receivers` by Monte Carlo
transport, with no view factors: for checking `view_factors` and `solve` against
independent ray tracing. Sources are those of `solve`: `ηPᵢ/Aᵢ`, uniform over each
element, from the absorbed laser power `laser_power` (per element).

Paths are traced backward from the receivers (reciprocity of the Lambertian kernel): a
uniform point on the receiver's facets, a cosine-weighted direction, the first facet hit
(`closest_hit` over every triangle, apertures included), whose element contributes its
source; the path continues diffusely with its weight multiplied by the element's albedo
(Russian roulette below 1/20), and ends on an aperture (escape), a one-sided surface's
back, or a miss. A closing lip (`aperture(...; closure)`, a `GapClosure`) is the
exception among apertures: like a surface, it re-emits diffusely the fraction `albedo`
of what reaches it (its own ``α(t) = s α_{after}``, from `albedos(target, time)` as `solve`
uses) and lets the rest escape; it has no source. Apertures stay transparent to laser
rays and lines of sight elsewhere; here the lip plane is the enclosure's boundary.
`albedo` defaults to the albedos at time 0; pass `albedos(target, time)` (or a
solution's `albedo`) to compare with a solve at a later time. Each path's sum is an unbiased estimate of the incident flux, as a
forward photon tally would be, but every sample lands on a receiver.

`reemit = :hit` continues from the hit point: the transport of the continuum, with
radiosity varying within elements (the limit the element model converges to).
`reemit = :element` continues from a uniform point on the hit element: the element
model with exact view factors, visibility included, which `solve` should reproduce
up to its view-factor error.

Geometry is the facets throughout: a `CurvedSurface` is ignored, so compare against a
target with flat elements (the same meshes without surfaces). `H` is the mean over
`batches` independent batches of `paths ÷ batches` paths per receiver (`batch`, one
row per batch, gives statistical errors of any functional of `H`). Batch `b` draws from
its own `Xoshiro` stream (from `seed` and `b`) and batches run on threads: results do
not depend on the thread count.
"""
function montecarlo_incident(t::Target, P::AbstractVector, receivers::AbstractVector{<:Integer};
                             paths::Integer = 10_000, batches::Integer = 32, reemit::Symbol = :hit,
                             albedo = albedos(t), seed::Integer = 1)
    reemit in (:hit, :element) || throw(ArgumentError("reemit must be :hit or :element, got :$reemit"))
    length(P) == nelements(t) || throw(DimensionMismatch("one laser power per element"))
    S = conversions(t) .* P ./ t.elem_area
    α = collect(Float64, albedo)
    nt = ntriangles(t)
    bvh = _sah_bvh(BVH(t.tri_v0, t.tri_e1, t.tri_e2, Int32.(1:nt)))   # apertures included
    # open apertures end a path; a closing lip re-emits its albedo like a surface
    aper = [(c = t.components[t.tri_comp[k]]; c.kind === :aperture && !_transmits(c.material.albedo))
            for k in 1:nt]
    # facet choice by area within each element
    cum = [cumsum(t.tri_area[t.elem_tris[e]]) for e in 1:nelements(t)]
    n = max(1, paths ÷ batches)
    out = zeros(batches, length(receivers))
    Threads.@threads for b in 1:batches
        _mc_batch!(view(out, b, :), t, bvh, aper, cum, S, α, receivers, n, reemit,
                   Random.Xoshiro(1_000_003 * seed + b))
    end
    return vec(sum(out; dims = 1)) ./ batches, out
end

function _mc_batch!(out, t, bvh, aper, cum, S, α, receivers, n, reemit, rng)
    eps = _ray_eps(t)
    for (r, e) in enumerate(receivers)
        acc = 0.0
        for _ in 1:n
            x, nrm = _mc_point(t, cum, e, rng)
            w = 1.0
            while true
                d = _cosine_direction(nrm, rng)
                s, k = closest_hit(bvh, x + eps * nrm, d)
                (k == 0 || aper[k]) && break                  # lost, or out through an open aperture
                dn = dot(d, t.tri_normal[k])
                e2 = dn < 0 ? t.tri_front[k] : t.tri_back[k]
                e2 == 0 && break                              # back of a one-sided surface
                acc += w * S[e2]
                w *= α[e2]
                w == 0 && break                               # a black element, or an open lip
                if w < 0.05
                    rand(rng) < 20w || break
                    w = 0.05
                end
                if reemit === :hit
                    x = x + eps * nrm + s * d
                    nrm = dn < 0 ? t.tri_normal[k] : -t.tri_normal[k]
                else
                    x, nrm = _mc_point(t, cum, e2, rng)
                end
            end
        end
        out[r] = acc / n
    end
    return out
end

# A uniform point on element e's facets, with the normal of its radiating side.
function _mc_point(t, cum, e, rng)
    c = cum[e]
    i = min(searchsortedfirst(c, rand(rng) * c[end]), length(c))
    k = t.elem_tris[e][i]
    u, v = rand(rng), rand(rng)
    u + v > 1 && ((u, v) = (1 - u, 1 - v))
    return t.tri_v0[k] + u * t.tri_e1[k] + v * t.tri_e2[k], t.elem_side[e] * t.tri_normal[k]
end

function _cosine_direction(n::Vec3, rng)
    a = abs(n[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)
    t1 = normalize(cross(n, a))
    t2 = cross(n, t1)
    u, φ = rand(rng), 2π * rand(rng)
    ρ = sqrt(u)
    return ρ * cos(φ) * t1 + ρ * sin(φ) * t2 + sqrt(1 - u) * n
end
