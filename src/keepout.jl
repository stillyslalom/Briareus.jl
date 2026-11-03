# Differentiable aperture keep-out: worst-case clearance of beam envelopes from the rims
# of entrance holes, smooth enough to optimize against.

# A piece of a hole's rim, sampled at `taus` in [0, 1]: the straight edge from `p` to `q`,
# or (`arc`) the circular arc about `o` of radius `ρ` from angle `θ0` through `Δθ` in the
# hole plane's basis `e1`, `e2`.
struct _RimPiece
    arc::Bool
    p::Vec3
    q::Vec3
    o::Vec3
    e1::Vec3
    e2::Vec3
    ρ::Float64
    θ0::Float64
    Δθ::Float64
    taus::Vector{Float64}
end

_rim_point(r::_RimPiece, τ) =
    r.arc ? r.o + r.ρ * (cos(r.θ0 + τ * r.Δθ) * r.e1 + sin(r.θ0 + τ * r.Δθ) * r.e2) :
            r.p + τ * (r.q - r.p)

"""
    ApertureKeepOut(target; apertures = planar aperture components, peak_fraction = 0.01,
                    intensity = nothing, pointing = 0.05, target_offset = 0.025,
                    hardware_tolerance = 0.0, focus = 1.0, margin = 0.0,
                    minimum = 0.0, weight = 100.0, nrim = 96)

Keep-out constraint for beams entering through planar apertures of any shape (laser
entrance holes, slots). `aperture_clearance(keepout, beams)` gives each beam's
worst-case clearance (mm) from the rim of its entrance aperture — the aperture on the
beam's port side — judged as `clearance` judges it: the beam envelope is the best-focus
contour at `peak_fraction` of the peak (or at `intensity`, W/cm²) grown by the f-cone
radius `|s|/(2F)` and `margin`; the clearance is the smallest distance from the rim to
the envelope, less the transverse error radius `pointing + target_offset +
hardware_tolerance` and `focus/(2F)` for the focus error.

The rim is the hole's boundary: the edges of the aperture's triangles (with a closing
lip in its plane about the same center, of the union) that no other of its triangles
continues across. On a `DiskSurface` aperture an edge whose ends lie at one radius from
the surface's center is the circular arc between them (the element's exact shape), so
a circular hole meshed as a polygon keeps its circle; every other edge is straight
(rectangular, slotted or elliptical openings, CAD apertures: the mesh is the geometry).
The rim is sampled at about `nrim` points, every straight edge at both ends, and the
smallest distance is refined by golden section on the sample intervals beside the
closest samples, to about 1e-9 mm along the rim.

The distance to the envelope ellipse is first order (Sampson): exact on the envelope
and for round spots, and smooth in aim, defocus and spot size. Confirm a design with
`clearance`, which also accounts for wings and power on hardware. Pass it as
`keepout` to `symmetry_objective` or `asymmetry_objective`, which then add
`weight × Σᵢ max(0, minimum − cᵢ)²`.
"""
struct ApertureKeepOut
    names::Vector{Symbol}
    centers::Vector{Vec3}
    rims::Vector{Vector{_RimPiece}}
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
    nrim >= 3 || throw(ArgumentError("nrim must be at least 3"))
    # by default every planar aperture; open windows in a hohlraum wall (apertures on a
    # `CylinderSurface`) are not entrance holes, and a closing lip (an aperture with a
    # `GapClosure`) widens the hole whose plane and center it shares
    lip(c) = c.kind === :aperture && _transmits(c.material.albedo)
    names, holes = apertures === nothing ? _default_apertures(t, lip) :
                   (ns = Symbol[apertures...]; (ns, [Set(part_pieces(t, n)) for n in ns]))
    lips = [ci for (ci, c) in enumerate(t.components) if lip(c)]
    isempty(names) && throw(ArgumentError("target has no apertures"))
    centroid(ks) = sum(t.tri_area[k] * (t.tri_v0[k] + (t.tri_e1[k] + t.tri_e2[k]) / 3) for k in ks) /
                   sum(t.tri_area[ks])
    centers, rims = Vec3[], Vector{_RimPiece}[]
    for (name, cis) in zip(names, holes)
        tris = findall(in(cis), t.tri_comp)
        isempty(tris) && throw(ArgumentError("aperture $name has no triangles"))
        c = centroid(tris)
        n = t.tri_normal[tris[1]]
        r = maximum(norm(x - c) for k in tris
                    for x in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k]))
        comps = collect(cis)
        for ci in lips
            ks = findall(==(ci), t.tri_comp)
            isempty(ks) && continue
            norm(centroid(ks) - c) < 1e-6 * r && abs(dot(t.tri_normal[ks[1]], n)) > 1 - 1e-9 &&
                push!(comps, ci)
        end
        push!(centers, c)
        push!(rims, _hole_rim(t, comps, n, nrim))
    end
    return ApertureKeepOut(names, centers, rims, peak_fraction, intensity,
                           pointing + target_offset + hardware_tolerance, focus, margin,
                           minimum, weight)
end

# The rim of the hole made of the coplanar components `comps` (normal `n`): the edges of
# their triangles that no other of them continues across (an edge two triangles share,
# or one with another triangle just beyond it: a T-junction, a crack between rings), each
# a straight piece or an arc (see `ApertureKeepOut`), sampled with about `nrim` points.
function _hole_rim(t::Target, comps, n, nrim)
    e1 = normalize(abs(n[1]) < 0.9 ? cross(n, Vec3(1, 0, 0)) : cross(n, Vec3(0, 1, 0)))
    e2 = cross(n, e1)
    tris = Tuple{Vec3,Vec3,Vec3,Int}[]
    for ci in comps, k in findall(==(ci), t.tri_comp)
        push!(tris, (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k], ci))
    end
    edgekey(a, b) = Tuple(a) < Tuple(b) ? (a, b) : (b, a)
    uses = Dict{Tuple{Vec3,Vec3},Int}()
    for (a, b, c, _) in tris, (x, y) in ((a, b), (b, c), (c, a))
        uses[edgekey(x, y)] = get(uses, edgekey(x, y), 0) + 1
    end
    function inside(x, tri)
        a, b, c = tri[1] - x, tri[2] - x, tri[3] - x
        s1, s2, s3 = dot(cross(a, b), n), dot(cross(b, c), n), dot(cross(c, a), n)
        return (s1 > 0 && s2 > 0 && s3 > 0) || (s1 < 0 && s2 < 0 && s3 < 0)
    end
    angle(x, o) = atan(dot(x - o, e2), dot(x - o, e1))
    sweep(x, y, o) = rem2pi(angle(y, o) - angle(x, o), RoundNearest)
    # (ends, arc: center and radius, the side of the arc its element lies on)
    cand = Tuple{Vec3,Vec3,Union{Nothing,Tuple{Vec3,Float64}},Int}[]
    for (a, b, c, ci) in tris, (x, y, z) in ((a, b, c), (b, c, a), (c, a, b))
        uses[edgekey(x, y)] == 1 || continue
        # a point just beyond the edge, away from its triangle: inside another one, the
        # edge lies inside the hole
        m = (x + y) / 2
        out = normalize(cross(y - x, n))
        dot(out, z - m) > 0 && (out = -out)
        probe = m + 1e-6 * norm(y - x) * out
        any(tri -> inside(probe, tri), tris) && continue
        s = t.components[ci].surface
        arc, side = nothing, 0
        if s isa DiskSurface
            o = s.center - dot(s.center - x, n) * n
            ρx, ρy = norm(x - o), norm(y - o)
            if ρx > 0 && abs(ρx - ρy) <= 1e-9 * ρx
                arc = (o, (ρx + ρy) / 2)
                side = norm(z - o) < arc[2] ? -1 : 1
            end
        end
        push!(cand, (x, y, arc, side))
    end
    # arcs of one circle with elements on both sides (rings of different counts: their
    # exact sectors meet though their chords leave a crack) lie inside the hole
    function covered((x, y, arc, side))
        arc === nothing && return false
        o, ρ = arc
        θm = angle(x, o) + sweep(x, y, o) / 2
        return any(cand) do (p, q, a2, s2)
            (a2 === nothing || s2 != -side || norm(a2[1] - o) > 1e-9 * ρ || abs(a2[2] - ρ) > 1e-9 * ρ) &&
                return false
            Δ = sweep(p, q, o)
            θs = Δ >= 0 ? angle(p, o) : angle(q, o)
            return rem2pi(θm - θs + 1e-9, RoundDown) < abs(Δ) + 2e-9     # a midpoint at a vertex
        end
    end
    edges = [(x, y, arc) for (x, y, arc, side) in filter(!covered, cand)]
    len((x, y, arc)) = arc === nothing ? norm(y - x) : abs(sweep(x, y, arc[1])) * arc[2]
    pieces = _RimPiece[]
    isempty(edges) && return pieces
    P = sum(len, edges)
    for e in edges
        x, y, arc = e
        taus = collect(range(0.0, 1.0; length = max(1, ceil(Int, nrim * len(e) / P)) + 1))
        push!(pieces, arc === nothing ? _RimPiece(false, x, y, x, e1, e2, 0.0, 0.0, 0.0, taus) :
                      _RimPiece(true, x, y, arc[1], e1, e2, arc[2], angle(x, arc[1]),
                                sweep(x, y, arc[1]), taus))
    end
    return pieces
end

# The default entrance holes: the planar aperture components (not lips), one hole per part
# and plane. A part's pieces in one plane are one hole (a CAD part's hole split into
# pieces); a part with holes in several planes (a hohlraum part: `Component.part` is the
# part on every piece) has one per plane, named by its first piece there. Returns the
# names and each hole's component indices.
function _default_apertures(t::Target, lip)
    cs = [ci for (ci, c) in enumerate(t.components)
          if c.kind === :aperture && !(c.surface isa CylinderSurface) && !lip(c)]
    plane(ci) = (k = findfirst(==(ci), t.tri_comp); n = t.tri_normal[k]; (n, dot(n, t.tri_v0[k])))
    names, holes, keys = Symbol[], Set{Int}[], Tuple{Symbol,Vec3,Float64}[]
    for ci in cs
        p = t.components[ci].part
        n, d = plane(ci)
        j = findfirst(((q, m, e),) -> q === p && abs(dot(m, n)) > 1 - 1e-9 &&
                                      abs(e - sign(dot(m, n)) * d) <= 1e-9 * max(1.0, abs(d)), keys)
        if j === nothing
            push!(keys, (p, n, d)); push!(holes, Set(ci)); push!(names, p)
        else
            push!(holes[j], ci)
        end
    end
    # a part with several holes: each by its first piece's name
    for j in eachindex(names)
        count(==(keys[j][1]), first.(keys)) > 1 && (names[j] = t.components[minimum(holes[j])].name)
    end
    return names, holes
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
    rim = k.rims[j]
    isempty(rim) && return _no_clearance(b)        # an aperture without a rim
    ρ = envelope_radius(b.spot, b.power; intensity = k.intensity, peak_fraction = k.peak_fraction)
    A, B = ρ * b.spot.a, ρ * b.spot.b
    kv = direction(b)
    e1, e2 = spot_axes(b)
    focus = b.aim + b.defocus * kv
    cone = 1 / (2b.fnumber)
    # the distance from the rim point at `τ` on piece `r` to the envelope
    function dist(r, τ)
        d = _rim_point(r, τ) - focus
        s, u, v = dot(d, kv), dot(d, e1), dot(d, e2)
        return _sampson_distance(u, v, A, B) - cone * abs(s) - k.margin
    end
    # the samples, and the intervals beside the closest ones (four: both sides of the
    # closest sample, on two pieces where it is a vertex), refined by golden section
    best = dist(rim[1], 0.0)
    far = (oftype(best, Inf), 1, 1)         # (interval's smaller end, piece, interval)
    near = (far, far, far, far)
    for (i, r) in enumerate(rim)
        prev = dist(r, r.taus[1])
        best = min(best, prev)
        for m in 1:length(r.taus)-1
            next = dist(r, r.taus[m+1])
            best = min(best, next)
            lo = min(prev, next)
            w = argmax(map(first, near))
            lo < near[w][1] && (near = Base.setindex(near, (lo, i, m), w))
            prev = next
        end
    end
    for (_, i, m) in near
        r = rim[i]
        best = min(best, _golden_min(τ -> dist(r, τ), r.taus[m], r.taus[m+1]))
    end
    return best - k.transverse - k.focus * cone
end

# The smallest value of `g` on [a, b] found by golden-section search (the abscissae stay
# Float64; only `g`'s values carry its number type).
function _golden_min(g, a::Float64, b::Float64; iterations = 40)
    φ = (sqrt(5.0) - 1) / 2
    c, d = b - φ * (b - a), a + φ * (b - a)
    gc, gd = g(c), g(d)
    for _ in 1:iterations
        if gc < gd
            b, d, gd = d, c, gc
            c = b - φ * (b - a)
            gc = g(c)
        else
            a, c, gc = c, d, gd
            d = a + φ * (b - a)
            gd = g(d)
        end
    end
    return min(gc, gd)
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
