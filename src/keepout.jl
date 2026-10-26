# Differentiable aperture keep-out: worst-case clearance of beam envelopes from the rims
# of entrance holes, smooth enough to optimize against.

"""
    ApertureKeepOut(target; apertures = planar aperture components, peak_fraction = 0.01,
                    intensity = nothing, pointing = 0.05, target_offset = 0.025,
                    hardware_tolerance = 0.0, focus = 1.0, margin = 0.0,
                    minimum = 0.0, weight = 100.0, nrim = 96)

Keep-out constraint for beams entering through circular apertures (laser entrance
holes). `aperture_clearance(keepout, beams)` gives each beam's worst-case clearance
(mm) from the rim of its entrance aperture — the aperture on the beam's port side —
judged as `clearance` judges it: the beam envelope is the best-focus contour at
`peak_fraction` of the peak (or at `intensity`, W/cm²) grown by the f-cone radius
`|s|/(2F)` and `margin`; the clearance is the smallest distance from `nrim` rim points
to the envelope, less the transverse error radius `pointing + target_offset +
hardware_tolerance` and `focus/(2F)` for the focus error. The rim is the circle
through the aperture mesh's outermost vertices.

The distance to the envelope ellipse is first order (Sampson): exact on the envelope
and for round spots, and smooth in aim, defocus and spot size. Confirm a design with
`clearance`, which also accounts for wings and power on hardware. Pass it as
`keepout` to `symmetry_objective` or `asymmetry_objective`, which then add
`weight × Σᵢ max(0, minimum − cᵢ)²`.
"""
struct ApertureKeepOut
    names::Vector{Symbol}
    centers::Vector{Vec3}
    rims::Vector{Vector{Vec3}}
    peak_fraction::Union{Nothing,Float64}
    intensity::Union{Nothing,Float64}
    transverse::Float64
    focus::Float64
    margin::Float64
    minimum::Float64
    weight::Float64
end

function ApertureKeepOut(t::Target; apertures = nothing, peak_fraction = nothing,
                         intensity = nothing, pointing = 0.05, target_offset = 0.025,
                         hardware_tolerance = 0.0, focus = 1.0, margin = 0.0,
                         minimum = 0.0, weight = 100.0, nrim::Integer = 96)
    intensity === nothing && peak_fraction === nothing && (peak_fraction = 0.01)
    # by default every planar aperture; open windows in a hohlraum wall (apertures on a
    # `CylinderSurface`) are not entrance holes, and a closing lip (an aperture with a
    # `GapClosure`) widens the hole whose plane and center it shares
    lip(c) = c.kind === :aperture && _transmits(c.material.albedo)
    names = apertures === nothing ?
            unique!([c.part for c in t.components
                     if c.kind === :aperture && !(c.surface isa CylinderSurface) && !lip(c)]) :
            Symbol[apertures...]
    lips = [findall(==(ci), t.tri_comp) for (ci, c) in enumerate(t.components) if lip(c)]
    isempty(names) && throw(ArgumentError("target has no apertures"))
    centers, rims = Vec3[], Vector{Vec3}[]
    for name in names
        cis = Set(part_pieces(t, name))           # a part: all its pieces
        tris = findall(in(cis), t.tri_comp)
        A = sum(t.tri_area[tris])
        c = sum(t.tri_area[k] * (t.tri_v0[k] + (t.tri_e1[k] + t.tri_e2[k]) / 3) for k in tris) / A
        n = t.tri_normal[tris[1]]
        reach(ks) = maximum(max(norm(t.tri_v0[k] - c), norm(t.tri_v0[k] + t.tri_e1[k] - c),
                                norm(t.tri_v0[k] + t.tri_e2[k] - c)) for k in ks)
        r = reach(tris)
        for ks in lips
            cl = sum(t.tri_area[k] * (t.tri_v0[k] + (t.tri_e1[k] + t.tri_e2[k]) / 3) for k in ks) /
                 sum(t.tri_area[ks])
            norm(cl - c) < 1e-6 * r && abs(dot(t.tri_normal[ks[1]], n)) > 1 - 1e-9 && (r = max(r, reach(ks)))
        end
        e1 = normalize(abs(n[1]) < 0.9 ? cross(n, Vec3(1, 0, 0)) : cross(n, Vec3(0, 1, 0)))
        e2 = cross(n, e1)
        push!(centers, c)
        push!(rims, [c + r * (cos(2π * j / nrim) * e1 + sin(2π * j / nrim) * e2) for j in 0:nrim-1])
    end
    return ApertureKeepOut(names, centers, rims, peak_fraction, intensity,
                           pointing + target_offset + hardware_tolerance, focus, margin,
                           minimum, weight)
end

"""
    aperture_clearance(keepout, beams) -> Vector

Worst-case clearance (mm) of each beam from the rim of its entrance aperture (see
`ApertureKeepOut`); `Inf` for a beam with no aperture on its port side. Differentiable
in the beams' aim, defocus and spot size.
"""
aperture_clearance(k::ApertureKeepOut, beams::AbstractVector{<:Beam}) =
    [_aperture_clearance(k, b) for b in beams]

# `Inf` in the beam's number type: a Float64 `Inf` beside dual clearances made the
# clearance a `Union` and the objectives' keep-out penalty infer as `Any`.
_no_clearance(::Beam{T}) where {T} = T(Inf)

function _aperture_clearance(k::ApertureKeepOut, b::Beam)
    side = [dot(c, b.port) for c in k.centers]
    j = argmax(side)
    side[j] > 0 || return _no_clearance(b)
    ρ = envelope_radius(b.spot, b.power; intensity = k.intensity, peak_fraction = k.peak_fraction)
    A, B = ρ * b.spot.a, ρ * b.spot.b
    kv = direction(b)
    e1, e2 = spot_axes(b)
    focus = b.aim + b.defocus * kv
    cone = 1 / (2b.fnumber)
    best = nothing
    for x in k.rims[j]
        d = x - focus
        s, u, v = dot(d, kv), dot(d, e1), dot(d, e2)
        dist = _sampson_distance(u, v, A, B) - cone * abs(s) - k.margin
        best = best === nothing ? dist : min(best, dist)
    end
    best === nothing && return _no_clearance(b)   # an aperture without a rim
    return best - k.transverse - k.focus * cone
end

# First-order (Sampson) signed distance from (u, v) to the ellipse with semi-axes A, B:
# exact on the ellipse and for circles, smooth everywhere but the center.
function _sampson_distance(u, v, A, B)
    q = sqrt((u / A)^2 + (v / B)^2)
    g = sqrt((u / A^2)^2 + (v / B^2)^2)
    return q > 1e-9 ? (q - 1) * q / g : -min(A, B) + zero(q)
end

"""
    SphereKeepOut(target; spheres = parts on a `SphereSurface`, peak_fraction = 0.01,
                  intensity = nothing, pointing = 0.05, target_offset = 0.025,
                  hardware_tolerance = 0.0, focus = 1.0, margin = 0.0, minimum = 0.0,
                  weight = 100.0)

Keep-out constraint for spherical obstacles inside a target, such as the capsule of a
hohlraum. `sphere_clearance(keepout, beams)` gives each beam's worst-case clearance (mm)
from the nearest of the `spheres` it passes before landing, with the envelope and error
set of `ApertureKeepOut`: the distance from the sphere's center to the beam envelope
(best-focus contour at `peak_fraction` of the peak, or at `intensity`, grown by the
f-cone at the center's axial position and by `margin`) less the sphere's radius, the
transverse error radius `pointing + target_offset + hardware_tolerance` and
`focus/(2F)`. The distance to the envelope is first order (Sampson), as for apertures;
the f-cone tilt of the envelope surface (cos of the cone half-angle) is neglected.

By default the spheres are the parts that are spheres: those whose largest spherical
component carries at least 90% of the part's area (a primitive capsule, or a CAD capsule
fitted with `curved_components`, its fill tube a piece of its own), not the small
spherical fillets fitted on hohlraum CAD. `spheres` names parts or components; a part
stands for its largest sphere, and its other pieces are not landing surfaces.

A sphere counts only when its center lies upstream of where the beam axis first meets
another opaque component (the wall a hohlraum beam lands on), found on the nominal
beam; a sphere behind that surface is shadowed and ignored. The penalty added by the
objectives is `weight × Σᵢ max(0, minimum − cᵢ)²`. The deposition penalty `direct`
(squared fraction of power on the capsule) vanishes as the spot edge starts to clip and
cannot keep a margin; use this keep-out when the capsule must stay clear. Confirm with
`clearance` on the capsule.
"""
struct SphereKeepOut
    names::Vector{Symbol}
    centers::Vector{Vec3}
    radii::Vector{Float64}
    target::Target
    others::Set{Int}                     # components a beam lands on
    peak_fraction::Union{Nothing,Float64}
    intensity::Union{Nothing,Float64}
    transverse::Float64
    focus::Float64
    margin::Float64
    minimum::Float64
    weight::Float64
end

function SphereKeepOut(t::Target; spheres = nothing, peak_fraction = nothing,
                       intensity = nothing, pointing = 0.05, target_offset = 0.025,
                       hardware_tolerance = 0.0, focus = 1.0, margin = 0.0, minimum = 0.0,
                       weight = 100.0)
    intensity === nothing && peak_fraction === nothing && (peak_fraction = 0.01)
    area(ci) = sum(t.elem_area[e] for e in t.comp_elems[ci] if t.elem_side[e] == 1; init = 0.0)
    # the largest sphere among the components `cis` (nothing if none)
    function sphere_of(cis)
        sph = [ci for ci in cis if t.components[ci].kind === :surface &&
                                   t.components[ci].surface isa SphereSurface]
        isempty(sph) ? nothing : sph[argmax(area.(sph))]
    end
    names, comps = Symbol[], Int[]
    if spheres === nothing
        # parts that are spheres: a fitted CAD capsule (with its fill tube) is, the small
        # spherical fillets fitted on hohlraum CAD are not
        for p in parts(t)
            cis = part_pieces(t, p)
            k = sphere_of(cis)
            k === nothing || area(k) < 0.9 * sum(area, cis) || (push!(names, p); push!(comps, k))
        end
    else
        for n in spheres
            k = sphere_of(part_pieces(t, Symbol(n)))
            k === nothing && throw(ArgumentError("component :$n does not carry a SphereSurface"))
            push!(names, Symbol(n)); push!(comps, k)
        end
    end
    isempty(names) && throw(ArgumentError("target has no spherical components"))
    centers, radii = Vec3[], Float64[]
    for k in comps
        s = t.components[k].surface
        push!(centers, s.center)
        push!(radii, s.radius)
    end
    idx = Set(ci for n in names for ci in part_pieces(t, n))
    others = Set(i for (i, c) in enumerate(t.components) if c.kind === :surface && !(i in idx))
    return SphereKeepOut(names, centers, radii, t, others, peak_fraction, intensity,
                         pointing + target_offset + hardware_tolerance, focus, margin,
                         minimum, weight)
end

"""
    sphere_clearance(keepout, beams) -> Vector

Worst-case clearance (mm) of each beam from the spheres of a `SphereKeepOut` it passes
before landing; `Inf` for a beam that passes none. Differentiable in the beams' aim,
defocus and spot size.
"""
sphere_clearance(k::SphereKeepOut, beams::AbstractVector{<:Beam}) =
    [_sphere_clearance(k, b) for b in beams]

function _sphere_clearance(k::SphereKeepOut, b::Beam)
    ρ = envelope_radius(b.spot, b.power; intensity = k.intensity, peak_fraction = k.peak_fraction)
    A, B = ρ * b.spot.a, ρ * b.spot.b
    kv = direction(b)
    e1, e2 = spot_axes(b)
    focus = b.aim + b.defocus * kv
    cone = 1 / (2b.fnumber)
    s_stop = _landing(k.target, BeamEnvelope(_value_beam(b)), k.others)
    best = nothing
    for (c, R) in zip(k.centers, k.radii)
        d = c - focus
        s, u, v = dot(d, kv), dot(d, e1), dot(d, e2)
        ForwardDiff.value(s) <= s_stop || continue          # behind the landing surface
        dist = _sampson_distance(u, v, A, B) - cone * abs(s) - k.margin - R
        best = best === nothing ? dist : min(best, dist)
    end
    best === nothing && return _no_clearance(b)
    return best - k.transverse - k.focus * cone
end

_keepout_clearance(k::ApertureKeepOut, beams) = aperture_clearance(k, beams)
_keepout_clearance(k::SphereKeepOut, beams) = sphere_clearance(k, beams)

_keepout_penalty(::Nothing, beams) = 0.0

# Several keep-outs (a tuple or vector) add their penalties. A tuple is unrolled with
# `map`: a generator over a heterogeneous tuple took Julia 1.13's inference 27 s.
_keepout_penalty(ks::Tuple, beams) = sum(map(k -> _keepout_penalty(k, beams), ks); init = 0.0)

function _keepout_penalty(ks::AbstractVector, beams)
    J = 0.0
    for k in ks
        J += _keepout_penalty(k, beams)
    end
    return J
end

function _keepout_penalty(k::Union{ApertureKeepOut,SphereKeepOut}, beams)
    c = _keepout_clearance(k, beams)
    # in the beams' number type also when no clearance is finite (a Float64 0.0 beside
    # dual penalties made the sum over keep-outs a runtime dispatch)
    J = k.weight * sum(max(zero(ci), k.minimum - ci)^2 for ci in c if isfinite(ci); init = 0.0)
    return convert(eltype(c), J)
end
