# This file defines differentiable keep-outs. The aperture keep-out gives the worst-case
# clearance of beam envelopes from the rims of entrance holes, smooth enough to optimize
# against.

# A piece of the rim of a hole, sampled at `taus` in [0, 1]. The piece is the straight
# edge from `p` to `q`, or with `arc` the circular arc about `o` of radius `ρ` from angle
# `θ0` through `Δθ` in the basis `e1`, `e2` of the hole plane.
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

# d(rim point)/dτ
_rim_tangent(r::_RimPiece, τ) =
    r.arc ? r.ρ * r.Δθ * (-sin(r.θ0 + τ * r.Δθ) * r.e1 + cos(r.θ0 + τ * r.Δθ) * r.e2) : r.q - r.p

# The rim samples of a hole, piece after piece. Each sample holds the rim point, its
# derivative in τ, τ and the piece. For an arc sample after the first of its piece, it
# also holds the arc length from the sample before, which bounds the change of s over
# the interval, and the sagitta of that arc. On an edge these are -Inf and 0.
struct _RimSamples
    x::Vector{Vec3}
    dx::Vector{Vec3}
    τ::Vector{Float64}
    piece::Vector{Int}
    arc::Vector{Float64}
    sag::Vector{Float64}
end

function _RimSamples(rim::Vector{_RimPiece})
    s = _RimSamples(Vec3[], Vec3[], Float64[], Int[], Float64[], Float64[])
    for (p, r) in enumerate(rim), (i, τ) in enumerate(r.taus)
        push!(s.x, _rim_point(r, τ)); push!(s.dx, _rim_tangent(r, τ)); push!(s.τ, τ); push!(s.piece, p)
        θ = r.arc && i > 1 ? abs(r.Δθ) * (τ - r.taus[i-1]) : 0.0
        push!(s.arc, r.arc && i > 1 ? r.ρ * θ : -Inf)
        push!(s.sag, r.ρ * (1 - cos(θ / 2)))
    end
    return s
end

"""
    ApertureKeepOut(target; apertures = planar aperture components, peak_fraction = 0.01,
                    intensity = nothing, pointing = 0.05, target_offset = 0.025,
                    hardware_tolerance = 0.0, focus = 1.0, margin = 0.0,
                    minimum = 0.0, weight = 100.0, nrim = 96, softness = 1e-4)

Keep-out constraint for beams entering through planar apertures of any shape, such as
laser entrance holes and slots. `aperture_clearance(keepout, beams)` gives the
worst-case clearance (mm) of each beam from the rim of its entrance aperture, which is
the aperture on the port side of the beam. The clearance is judged as `clearance` judges
it. The beam envelope is the best-focus contour at `peak_fraction` of the peak (or at
`intensity`, W/cm²), grown by the f-cone radius `|s|/(2F)` and `margin`. The clearance
is the smallest distance from the rim to the envelope, less the transverse error radius
`pointing + target_offset + hardware_tolerance` and `focus/(2F)` for the focus error.

The rim is the boundary of the hole. It consists of the edges of the aperture triangles
that no other triangle of the aperture continues across; a closing lip in the plane of
the aperture about the same center is part of the union. On a `DiskSurface` aperture,
an edge whose ends lie at one radius from the center of the surface is the circular arc
between them, the exact shape of the element. The rim of a circular hole meshed as a
polygon is therefore its circle. Every other edge is straight, as for rectangular, slotted or
elliptical openings and CAD apertures, where the mesh is the geometry. The rim is
sampled at about `nrim` points, with every straight edge sampled at both ends. The least
distance on each sample interval is found exactly, at an end or where the slope of the
distance along the rim rises through zero. That point is located from the slopes at the
samples to 1e-12 of the interval. An interval is first split at these points:

- where the distance has a kink, which is where the rim crosses the plane of best focus
  and the f-cone term turns;
- where the rim crosses an axis of an elongated envelope, since the first-order
  distance below turns sharply there;
- inside the envelope, where the rim passes nearest the beam axis, since the first-order
  distance dips steeply there;
- where the cubic through the distances and slopes at the ends turns twice.

The rim of a dual beam is searched on its values. The partials of the result follow from
its sensitivities to the focus and the semi-axes of the envelope, so their cost does not
grow with the number of partials.

The clearance is a smooth minimum of these interval minima, `−softness · log Σ σₖ
exp(−dₖ/softness)` (`softness` in mm). The weight σₖ of each interval is 1 within
10 `softness` of the least and falls smoothly to 0 at 20 `softness`, where a term
weighs at most e⁻¹⁰ of the least. The exact minimum over the rim has kinks where the
worst rim point jumps between separate local minima of the distance, and the smooth
minimum rounds them off. It never exceeds the exact minimum, lies below it by at most
`softness` times the log of the number of intervals counted, and is the exact minimum
for `softness = 0`. The default, 1e-4 mm, is an engineering default well below the
agreement of the keep-out with `clearance`, which is a few micrometers. On NIF quads
through the entrance holes of a hohlraum, the smooth minimum lies 0 to 0.18 µm below the
exact one. Over 144 scans of the aim, the largest change of slope in a step falls from
0.21 (a kink) to 0.002, which remains where one interval contains two local minima
(research notes, "Defaults changed (2026-11-19)").

The distance to the envelope ellipse is first order (Sampson). It is exact on the
envelope and for round spots, and smooth in aim, defocus and spot size. A design should
be confirmed with `clearance`, which also accounts for wings and power on hardware. The
keep-out is passed as `keepout` to `symmetry_objective` or `asymmetry_objective`, which
then add `weight × Σᵢ max(0, minimum − cᵢ)²`.
"""
struct ApertureKeepOut <: AbstractKeepOut
    names::Vector{Symbol}
    centers::Vector{Vec3}
    rims::Vector{Vector{_RimPiece}}
    samples::Vector{_RimSamples}
    peak_fraction::Union{Nothing,Float64}
    intensity::Union{Nothing,Float64}
    transverse::Float64
    focus::Float64
    margin::Float64
    minimum::Float64
    weight::Float64
    softness::Float64      # length (mm) of the smooth minimum over the rim (0 for the exact one)
end

function ApertureKeepOut(t::Target; apertures = nothing, peak_fraction = nothing,
                         intensity = nothing, pointing = 0.05, target_offset = 0.025,
                         hardware_tolerance = 0.0, focus = 1.0, margin = 0.0,
                         minimum = 0.0, weight = 100.0, nrim::Integer = 96, softness::Real = 1e-4)
    intensity === nothing && peak_fraction === nothing && (peak_fraction = 0.01)
    nrim >= 3 || throw(ArgumentError("nrim must be at least 3"))
    softness >= 0 || throw(ArgumentError("softness must be non-negative"))
    # By default the holes are every planar aperture. Open windows in a hohlraum wall,
    # which are apertures on a `CylinderSurface`, are not entrance holes. A closing lip,
    # an aperture with a `GapClosure`, widens the hole whose plane and center it shares.
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
    return ApertureKeepOut(names, centers, rims, _RimSamples.(rims), peak_fraction, intensity,
                           pointing + target_offset + hardware_tolerance, focus, margin,
                           minimum, weight, Float64(softness))
end

# Returns the rim of the hole made of the coplanar components `comps` (normal `n`). The
# rim consists of the edges of their triangles that no other of their triangles continues
# across. An edge that two triangles share is excluded, and so is an edge with another
# triangle just beyond it, as at a T-junction or a crack between rings. Each rim edge is
# a straight piece or an arc (see `ApertureKeepOut`), sampled with about `nrim` points.
function _hole_rim(t::Target, comps, n, nrim)
    e1, e2 = _ortho_basis(n)
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
    # (ends, arc as center and radius, and the side of the arc on which its element lies)
    cand = Tuple{Vec3,Vec3,Union{Nothing,Tuple{Vec3,Float64}},Int}[]
    for (a, b, c, ci) in tris, (x, y, z) in ((a, b, c), (b, c, a), (c, a, b))
        uses[edgekey(x, y)] == 1 || continue
        # A point just beyond the edge, away from its triangle, is tested. When it lies
        # inside another triangle, the edge lies inside the hole.
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
    # Arcs of one circle with elements on both sides lie inside the hole. This happens for
    # rings of different counts, whose exact sectors meet though their chords leave a
    # crack.
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

# Returns the default entrance holes, which are the planar aperture components other than
# lips, with one hole per part and plane. The pieces of a part in one plane form one hole,
# as when the hole of a CAD part is split into pieces. A part with holes in several
# planes, such as a hohlraum part (`Component.part` is the part on every piece), has one
# hole per plane, named by its first piece there. Returns the names and the component
# indices of each hole.
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
    # each hole of a part with several holes is named by its first piece
    for j in eachindex(names)
        count(==(keys[j][1]), first.(keys)) > 1 && (names[j] = t.components[minimum(holes[j])].name)
    end
    return names, holes
end

"""
    aperture_clearance(keepout, beams) -> Vector

Worst-case clearance (mm) of each beam from the rim of its entrance aperture (see
`ApertureKeepOut`), and `Inf` for a beam with no aperture on its port side.
Differentiable in the aim, defocus and spot size of the beams.
"""
function aperture_clearance(k::ApertureKeepOut, beams::AbstractVector{<:Beam})
    # the envelope radius once per distinct spot and power (values), with a root search each
    seen = Tuple{SuperGaussianSpot{Float64},Float64,Float64}[]
    return [_aperture_clearance(k, b, _envelope_radius!(seen, k, b)) for b in beams]
end

function _envelope_radius!(seen, k::ApertureKeepOut, b::Beam)
    sp, P = _value_spot(b.spot), _fvalue(b.power)
    for (s, p, ρ) in seen
        s === sp && p === P && return ρ
    end
    ρ = envelope_radius(sp, P; intensity = k.intensity, peak_fraction = k.peak_fraction)
    push!(seen, (sp, P, ρ))
    return ρ
end

# `Inf` in the number type of the beam. A Float64 `Inf` beside dual clearances made the
# clearance a `Union`, and the keep-out penalty of the objectives then inferred as `Any`.
_no_clearance(::Beam{T}) where {T} = T(Inf)

function _aperture_clearance(k::ApertureKeepOut, b::Beam, ρ::Float64)
    j, side = 0, 0.0
    for (i, c) in enumerate(k.centers)
        s = dot(c, b.port)
        (j == 0 || s > side) && ((j, side) = (i, s))
    end
    side > 0 || return _no_clearance(b)
    rim = k.rims[j]
    isempty(rim) && return _no_clearance(b)        # an aperture without a rim
    fv = _rim_frame(k, b isa Beam{Float64} ? b : _value_beam(b), ρ)
    # The rim is searched on the values of the beam, and the terms found enter in its
    # number type.
    terms, imin = _rim_terms(k, rim, k.samples[j], fv)
    return _rim_result(k, terms, imin, fv, b, ρ) - k.transverse - k.focus * fv.cone
end

# Returns the envelope of the beam for the distance from a rim point (`_BeamConeFrame`),
# in the number type of the beam. `ρ` is the envelope radius, computed from the values of
# the beam.
function _rim_frame(k::ApertureKeepOut, b::Beam{T},
                    ρ = envelope_radius(b.spot, b.power; intensity = k.intensity, peak_fraction = k.peak_fraction)) where {T}
    kv = direction(b)
    e1, e2 = spot_axes(b)
    return _BeamConeFrame{T}(b.aim + b.defocus * kv, kv, e1, e2, ρ * b.spot.a, ρ * b.spot.b, 1 / (2b.fnumber), k.margin)
end

# The distance from the rim point `x` to the envelope.
function _rim_distance(f::_BeamConeFrame, x::Vec3)
    d = x - f.focus
    s, u, v = dot(d, f.k), dot(d, f.e1), dot(d, f.e2)
    return _sampson_distance(u, v, f.A, f.B) - f.cone * abs(s) - f.margin
end
_rim_distance(f::_BeamConeFrame, r::_RimPiece, τ::Float64) = _rim_distance(f, _rim_point(r, τ))

# The smooth minimum counts at full weight the sample intervals whose least distance lies
# within `_SOFT_FULL` softness of the least. Intervals up to `_SOFT_ZERO` count with a
# weight that falls (C¹) to 0 there. Where a weight falls, its term is at most
# e⁻¹⁰ = 4.5e-5 of the least term.
const _SOFT_FULL = 10.0
const _SOFT_ZERO = 20.0

# Returns the terms of the smooth minimum, found on the values of the beam (frame `f`).
# For each sample interval of the rim, the term is its least distance, at an end or at a
# minimum inside. The term is kept with its rim point if it lies within `_SOFT_ZERO`
# softness of the least; with `softness = 0` only the least is kept. Returns the terms
# (distance, point) and the index of the least.
function _rim_terms(k::ApertureKeepOut, rim::Vector{_RimPiece}, rs::_RimSamples, f::_BeamConeFrame{Float64})
    reach = _SOFT_ZERO * k.softness
    sf = _SlopeFrame(f)
    terms = sizehint!(Tuple{Float64,Vec3}[], 32)
    lo, cap = Inf, 32
    x0 = rs.x[1]
    d0, S0, s0, t0, w0 = _rim_slope(sf, x0, rs.dx[1])
    @inbounds for j in 2:length(rs.x)
        x1 = rs.x[j]
        d1, S1, s1, t1, w1 = _rim_slope(sf, x1, rs.dx[j])
        if rs.piece[j] == rs.piece[j-1]           # the interval from sample j − 1
            g, xg = d1 < d0 ? (d1, x1) : (d0, x0)
            g, xg = _interval_min(sf, rim, rs.piece[j], rs.τ[j-1], rs.τ[j], rs.arc[j], rs.sag[j],
                                  (S0, s0, t0, w0, d0), (S1, s1, t1, w1, d1), g, xg)
            if g <= lo + reach
                push!(terms, (g, xg))
                lo = min(lo, g)
                # (terms kept before the least was reached are dropped now and then)
                length(terms) == cap && (_drop_above!(terms, lo + reach); cap = max(32, 2length(terms)))
            end
        end
        x0, d0, S0, s0, t0, w0 = x1, d1, S1, s1, t1, w1
    end
    _drop_above!(terms, lo + reach)
    # a NaN distance, from an envelope of zero size, passes into the result
    isempty(terms) && push!(terms, (NaN, rs.x[1]))
    imin = 1
    for i in eachindex(terms)
        terms[i][1] < terms[imin][1] && (imin = i)
    end
    return terms, imin
end

_drop_above!(terms, cut) = filter!(t -> t[1] <= cut, terms)

# The Float64 frame of the search, with the inverse squares of the semi-axes.
struct _SlopeFrame
    focus::Vec3
    k::Vec3
    e1::Vec3
    e2::Vec3
    ia::Float64
    ib::Float64
    ia2::Float64
    ib2::Float64
    minab::Float64
    cone::Float64
    margin::Float64
    turn::Float64           # a quarter of the angle (envelope units) over which the direction factor turns at its peak
end
function _SlopeFrame(f::_BeamConeFrame{Float64})
    e = max(f.A, f.B) / min(f.A, f.B)
    return _SlopeFrame(f.focus, f.k, f.e1, f.e2, inv(f.A), inv(f.B), inv(f.A^2), inv(f.B^2), min(f.A, f.B), f.cone,
                       f.margin, e > 1 ? 0.25 / sqrt(e^2 - 1) : Inf)
end

# Returns, at the rim point `x` with tangent `t` (dx/dτ), the distance to the envelope,
# the slope of its Sampson part in τ, the axial position s and its slope, and the point in
# envelope units (u/A, v/B, inside the envelope within 1 of the axis). The slope of the
# distance is S′ − cone·sign(s)·s′, since |s| has a concave kink where s = 0.
@inline function _rim_slope(f::_SlopeFrame, x::Vec3, t::Vec3)
    d = x - f.focus
    s, u, v = dot(d, f.k), dot(d, f.e1), dot(d, f.e2)
    ts, tu, tv = dot(t, f.k), dot(t, f.e1), dot(t, f.e2)
    ua, vb = u * f.ia2, v * f.ib2
    q = sqrt(u * ua + v * vb)
    g = sqrt(ua * ua + vb * vb)
    if q > 1e-9
        D = (q - 1) * q / g
        dq = (ua * tu + vb * tv) / q
        dg = (ua * f.ia2 * tu + vb * f.ib2 * tv) / g
        dD = ((2q - 1) * dq - D * dg) / g
    else
        D, dD = -f.minab, 0.0
    end
    return D - f.cone * abs(s) - f.margin, dD, s, ts, SVector(u * f.ia, v * f.ib)
end

# Returns the point and tangent of piece `r` at τ, with one sincos on an arc.
function _rim_point_tangent(r::_RimPiece, τ::Float64)
    r.arc || return r.p + τ * (r.q - r.p), r.q - r.p
    sn, cs = sincos(r.θ0 + τ * r.Δθ)
    return r.o + r.ρ * (cs * r.e1 + sn * r.e2), r.ρ * r.Δθ * (-sn * r.e1 + cs * r.e2)
end

# Returns the slope of the distance at τ on piece `r`, on the branch `σ` = sign(s) of |s|.
function _slope_at(f::_SlopeFrame, r::_RimPiece, τ::Float64, σ)
    x, t = _rim_point_tangent(r, τ)
    _, S, _, ts = _rim_slope(f, x, t)
    return S - f.cone * σ * ts
end

# Returns the least distance on the sample interval [τa, τb] of piece `p` of `rim`. `arc`
# is the arc length of the interval (-Inf on an edge) and `sag` the sagitta of the arc.
# `g` at `xg` is the distance at the lesser end. `ea` and `eb` hold, for each end, the
# Sampson slope, s, s′, the point in envelope units and the distance. The interval is
# split (`_split_min`) at these points:
# - where |s| has a kink;
# - where the rim crosses an axis of the envelope, if the direction factor of the Sampson
#   distance can turn sharply there. For an elongated spot the factor peaks over about
#   1/√(e² − 1) of the direction in envelope units, with e the ratio of the semi-axes,
#   and the interval is split if the direction turns through a quarter of that;
# - if the interval reaches inside the envelope, where it passes nearest the beam axis.
# Each part on which the distance is smooth is searched by `_part_min`.
@inline function _interval_min(f::_SlopeFrame, rim, p, τa, τb, arc, sag, ea, eb, g, xg)
    Sa, sa, ta, wa, da = ea
    Sb, sb, tb, wb, db = eb
    c = f.cone
    # A kink lies inside when s changes sign, or when on an arc s passes an extremum near
    # enough to zero to cross it; |s′| is at most the arc length of the interval per unit τ.
    kink = c > 0 && ((sa > 0 && sb < 0) || (sa < 0 && sb > 0) || (ta * tb < 0 && min(abs(sa), abs(sb)) <= arc))
    # nearest approach of the chord to the beam axis, in envelope units, with the sag of
    # the arc
    Δ = wb - wa
    t = clamp(-dot(wa, Δ) / dot(Δ, Δ), 0.0, 1.0)
    inside = sum(abs2, wa + t * Δ) < (1 + sag * max(f.ia, f.ib))^2
    # an axis crossing (u = 0 or v = 0) through a turn of the direction wider than `turn`
    axes = (wa[1] * wb[1] < 0 || wa[2] * wb[2] < 0) &&
           abs(wa[1] * wb[2] - wa[2] * wb[1]) > f.turn * sqrt(sum(abs2, wa) * sum(abs2, wb))
    (kink || inside || axes) &&
        return _split_min(f, rim[p], τa, τb, ea, eb, g, xg, kink, inside || axes, inside ? τa + t * (τb - τa) : NaN)
    σ = sa != 0 ? sign(sa) : sign(sb)
    ml, mr = Sa - c * σ * ta, Sb - c * σ * tb
    if ml < 0 < mr
        r = rim[p]
        return _smaller(f, r, _slope_root(f, r, τa, τb, ml, mr, σ), g, xg)
    end
    _hermite_turns(τb - τa, da, db, ml, mr) || return g, xg
    return _part_min(f, rim, p, τa, τb, da, db, ml, mr, σ, g, xg)
end

# Whether the cubic Hermite interpolant on an interval of length `L` with end values `dl`,
# `dr` and slopes `ml`, `mr` of one sign turns twice inside, which is when its slope, a
# quadratic, changes sign at its vertex.
@inline function _hermite_turns(L, dl, dr, ml, mr)
    (ml < 0) == (mr < 0) || return false
    a0, a1, δ = L * ml, L * mr, dl - dr
    α, β = 6δ + 3a0 + 3a1, -6δ - 4a0 - 2a1
    α == 0 && return false
    tv = -β / (2α)
    return 0 < tv < 1 && (a0 > 0) != (a0 + tv * (β + α * tv) > 0)
end

function _smaller(f::_SlopeFrame, r::_RimPiece, τ, g, xg)
    x, t = _rim_point_tangent(r, τ)
    d = first(_rim_slope(f, x, t))
    return d < g ? (d, x) : (g, xg)
end

# Returns the least distance below `g` (at `xg`) on a part [l, rr] of piece `p` of `rim` on
# which the distance is smooth, with `σ` the branch of |s|. The search uses the distances
# `dl`, `dr` and slopes `ml`, `mr` at the ends, and a minimum inside lies where the slope
# rises through zero (`_slope_root`). With slopes of one sign at both ends, the cubic
# Hermite interpolant may turn twice inside, which gives a minimum and a maximum between
# the samples. The part is then split at the two turning points, each a candidate, and
# each piece is searched again, once.
function _part_min(f::_SlopeFrame, rim, p, l, rr, dl, dr, ml, mr, σ, g, xg, again::Bool = true)
    if ml < 0 < mr
        r = rim[p]
        return _smaller(f, r, _slope_root(f, r, l, rr, ml, mr, σ), g, xg)
    end
    (again && (ml < 0) == (mr < 0)) || return g, xg
    # the interpolant's slope α t² + β t + γ on t ∈ [0, 1]
    L = rr - l
    a0, a1, δ = L * ml, L * mr, dl - dr
    α, β, γ = 6δ + 3a0 + 3a1, -6δ - 4a0 - 2a1, a0
    disc = β^2 - 4α * γ
    (α != 0 && disc > 0) || return g, xg
    q = -(β + copysign(sqrt(disc), β)) / 2
    t1, t2 = minmax(q / α, γ / q)
    0 < t1 && t2 < 1 || return g, xg
    r = rim[p]
    τ1, τ2 = l + t1 * L, l + t2 * L
    x1, s1 = _rim_point_tangent(r, τ1)
    d1, S1, _, ts1 = _rim_slope(f, x1, s1)
    x2, s2 = _rim_point_tangent(r, τ2)
    d2, S2, _, ts2 = _rim_slope(f, x2, s2)
    d1 < g && ((g, xg) = (d1, x1))
    d2 < g && ((g, xg) = (d2, x2))
    m1, m2 = S1 - f.cone * σ * ts1, S2 - f.cone * σ * ts2
    g, xg = _part_min(f, rim, p, l, τ1, dl, d1, ml, m1, σ, g, xg, false)
    g, xg = _part_min(f, rim, p, τ1, τ2, d1, d2, m1, m2, σ, g, xg, false)
    return _part_min(f, rim, p, τ2, rr, d2, dr, m2, mr, σ, g, xg, false)
end

# `_interval_min` on an interval split at the kinks of |s|. With `axes`, the interval is
# also split where the rim crosses the axes of the envelope (u = 0, v = 0). There the
# direction factor of the Sampson distance is at its extremes, and between them it turns
# within the distance of the axis, sharply for an elongated spot. A further split is at
# `τq`, the nearest approach of the chord to the beam axis inside the envelope, where the
# Sampson distance (q − 1)·h falls to −h. Each split point is a candidate, and each part
# is searched by `_part_min` on its branch of |s|.
function _split_min(f::_SlopeFrame, r::_RimPiece, τa, τb, ea, eb, g, xg, kink, axes, τq)
    k1, k2 = kink ? _crossings(f, r, τa, τb, f.k) : (NaN, NaN)
    u1, u2 = axes ? _crossings(f, r, τa, τb, f.e1) : (NaN, NaN)
    v1, v2 = axes ? _crossings(f, r, τa, τb, f.e2) : (NaN, NaN)
    m = τa < τq < τb ? τq : NaN
    splits = sort(SVector(k1, k2, u1, u2, v1, v2, m, τb); by = τ -> isnan(τ) ? Inf : τ)
    c = f.cone
    rim = (r,)
    l, dl, Sl, tl = τa, ea[5], ea[1], ea[3]  # left end of the part: τ, distance, Sampson slope, s′
    for rr in splits
        (isnan(rr) || rr <= l) && continue
        if rr == τb
            dr, Sr, tr = eb[5], eb[1], eb[3]
        else
            x, t = _rim_point_tangent(r, rr)
            dr, Sr, _, tr = _rim_slope(f, x, t)
            dr < g && ((g, xg) = (dr, x))
        end
        # the branch of |s| on (l, rr)
        σ = kink ? sign(dot(_rim_point(r, (l + rr) / 2) - f.focus, f.k)) : ea[2] != 0 ? sign(ea[2]) : sign(eb[2])
        g, xg = _part_min(f, rim, 1, l, rr, dl, dr, Sl - c * σ * tl, Sr - c * σ * tr, σ, g, xg)
        rr == τb && break
        l, dl, Sl, tl = rr, dr, Sr, tr
    end
    return g, xg
end

# Returns the points inside (τa, τb) at which the rim piece `r` crosses the plane through
# the focus normal to `n` (dot(x − focus, n) = 0), in ascending order (NaN for none). The
# crossing condition is linear in τ on an edge and C + M cos(θ − φ) on an arc. With `n`
# the beam direction, the points are the kinks of |s|.
function _crossings(f::_SlopeFrame, r::_RimPiece, τa, τb, n)
    if !r.arc
        ds = dot(r.q - r.p, n)
        τ = ds == 0 ? NaN : dot(f.focus - r.p, n) / ds
        return τa < τ < τb ? (τ, NaN) : (NaN, NaN)
    end
    C = dot(r.o - f.focus, n)
    a1, a2 = r.ρ * dot(r.e1, n), r.ρ * dot(r.e2, n)
    M = hypot(a1, a2)
    M > abs(C) || return (NaN, NaN)
    φ, δ = atan(a2, a1), acos(-C / M)
    period = 2π / abs(r.Δθ)
    t1, t2 = _arc_tau(r, φ - δ, τa, τb, period), _arc_tau(r, φ + δ, τa, τb, period)
    isnan(t1) && return (t2, NaN)
    isnan(t2) && return (t1, NaN)
    return (min(t1, t2), max(t1, t2))
end

# Returns the τ in (τa, τb) at which arc `r` reaches the angle θ (modulo a turn, with
# `period` in τ), or NaN.
function _arc_tau(r::_RimPiece, θ, τa, τb, period)
    τ = τa + mod((θ - r.θ0) / r.Δθ - τa, period)
    return τa < τ < τb ? τ : NaN
end

# Returns the point in (a, b) where the slope `S′ − cone·σ·s′` of the distance along piece
# `r` rises through zero, with `la` < 0 at a and `lb` > 0 at b. The search is regula falsi
# with the Illinois step, to 1e-12 of the interval.
function _slope_root(f::_SlopeFrame, r::_RimPiece, a, b, la, lb, σ)
    tol = 1e-12 * (b - a)
    side = 0
    τ = a
    for it in 1:100
        τn = (a * lb - b * la) / (lb - la)
        a < τn < b || (τn = (a + b) / 2)
        step = abs(τn - τ)
        τ = τn
        l = _slope_at(f, r, τ, σ)
        l == 0 && break
        if l < 0
            a, la = τ, l
            side == -1 && (lb /= 2)
            side = -1
        else
            b, lb = τ, l
            side == 1 && (la /= 2)
            side = 1
        end
        (b - a <= tol || (it > 1 && step <= tol)) && break
    end
    return τ
end

# Returns the smooth minimum of the distances of the terms (`ApertureKeepOut`), evaluated
# in the frame `f`, as −λ log Σₖ σₖ exp(−dₖ/λ). The weight σₖ of each term is 1 within
# `_SOFT_FULL` λ of the least and falls (C¹) to 0 at `_SOFT_ZERO` λ.
function _rim_soft_min(k::ApertureKeepOut, terms, imin, f::_BeamConeFrame)
    λ = k.softness
    lo = _rim_distance(f, terms[imin][2])
    λ == 0 && return lo
    S = zero(lo)
    for (_, x) in terms
        z = (_rim_distance(f, x) - lo) / λ
        S += _soft_weight(z) * exp(-z)
    end
    return lo - λ * log(S)
end

function _soft_weight(z)
    z <= _SOFT_FULL && return one(z)
    z >= _SOFT_ZERO && return zero(z)
    t = (z - _SOFT_FULL) / (_SOFT_ZERO - _SOFT_FULL)
    return 1 - t^2 * (3 - 2t)
end

# Dual numbers of the keep-out. They hold the sensitivities of the smooth minimum to the
# focus (three) and to the semi-axes of the envelope (two), the only quantities of the
# beam that it uses.
struct _RimTag end
const _RimDual = ForwardDiff.Dual{_RimTag,Float64,5}

# Returns the smooth minimum in the number type of the beam. On Float64 it uses the terms
# on the values. On a dual beam it uses the sensitivities to the five dual quantities of
# the frame (`_RimDual`, with Float64 partials), carried to the partials of the beam by
# the chain rule, so the cost does not grow with the partials of the beam. On nested
# duals it uses the terms in the number type of the beam.
_rim_result(k::ApertureKeepOut, terms, imin, fv, ::Beam{Float64}, ρ) = _rim_soft_min(k, terms, imin, fv)

function _rim_result(k::ApertureKeepOut, terms, imin, fv, b::Beam{<:ForwardDiff.Dual{<:Any,Float64}}, ρ)
    e(i) = ForwardDiff.Partials{5,Float64}(ntuple(j -> Float64(j == i), 5))
    fs = _BeamConeFrame{_RimDual}(SVector(_RimDual(fv.focus[1], e(1)), _RimDual(fv.focus[2], e(2)), _RimDual(fv.focus[3], e(3))),
                                  fv.k, fv.e1, fv.e2, _RimDual(fv.A, e(4)), _RimDual(fv.B, e(5)), fv.cone, fv.margin)
    c = _rim_soft_min(k, terms, imin, fs)
    p = ForwardDiff.partials(c)
    fb = _rim_frame(k, b, ρ)
    Δ(x) = x - ForwardDiff.value(x)          # only the partials, at value 0
    return ForwardDiff.value(c) + p[1] * Δ(fb.focus[1]) + p[2] * Δ(fb.focus[2]) + p[3] * Δ(fb.focus[3]) +
           p[4] * Δ(fb.A) + p[5] * Δ(fb.B)
end

_rim_result(k::ApertureKeepOut, terms, imin, fv, b::Beam, ρ) = _rim_soft_min(k, terms, imin, _rim_frame(k, b, ρ))

# Returns the first-order (Sampson) signed distance from (u, v) to the ellipse with
# semi-axes A, B. It is exact on the ellipse and for circles, and smooth everywhere
# except at the center.
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
hohlraum. `sphere_clearance(keepout, beams)` gives the worst-case clearance (mm) of each
beam from the nearest of the `spheres` that it passes before landing, with the envelope
and error set of `ApertureKeepOut`. The clearance is the distance from the center of the
sphere to the beam envelope, less the radius of the sphere, the transverse error radius
`pointing + target_offset + hardware_tolerance` and `focus/(2F)`. The envelope is the
best-focus contour at `peak_fraction` of the peak, or at `intensity`, grown by the
f-cone at the axial position of the center and by `margin`. The distance to the
envelope is first order (Sampson), as for apertures, and the f-cone tilt of the envelope
surface (the cosine of the cone half-angle) is neglected.

By default the spheres are the parts that are spheres, which are those whose largest
spherical component covers at least 90% of the part area. Examples are a primitive
capsule, or a CAD capsule fitted with `curved_components` with its fill tube a separate
piece. The small spherical fillets fitted on hohlraum CAD are excluded. `spheres` names
parts or components. A part stands for its largest sphere, and its other pieces are not
landing surfaces.

A sphere counts only when its center lies upstream of the point where the beam axis
first meets another opaque component, such as the wall on which a hohlraum beam lands.
That point is found on the nominal beam, and a sphere behind that surface is shadowed
and ignored. The penalty added by the objectives is `weight × Σᵢ max(0, minimum − cᵢ)²`.
The deposition penalty `direct` (the squared fraction of power on the capsule) vanishes
as the spot edge starts to clip, so it cannot hold a margin. This keep-out serves when
the capsule must stay clear, and a design is confirmed with `clearance` on the capsule.
"""
struct SphereKeepOut <: AbstractKeepOut
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
    # the largest sphere among the components `cis` (nothing if there is none)
    function sphere_of(cis)
        sph = [ci for ci in cis if t.components[ci].kind === :surface &&
                                   t.components[ci].surface isa SphereSurface]
        isempty(sph) ? nothing : sph[argmax(area.(sph))]
    end
    names, comps = Symbol[], Int[]
    if spheres === nothing
        # A fitted CAD capsule with its fill tube counts as a sphere, and the small
        # spherical fillets fitted on hohlraum CAD do not.
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

Worst-case clearance (mm) of each beam from the spheres of a `SphereKeepOut` that it
passes before landing, and `Inf` for a beam that passes none. Differentiable in the aim,
defocus and spot size of the beams.
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

"""
    keepout_clearance(keepout, beams) -> Vector
    keepout_clearance(keepout, ev::Evaluation) -> Vector

The worst-case clearance (mm) of each beam from a keep-out, which is the quantity that
its penalty holds at or above the `minimum` of the keep-out. It is
[`aperture_clearance`](@ref) for an [`ApertureKeepOut`](@ref),
[`sphere_clearance`](@ref) for a [`SphereKeepOut`](@ref) and
[`hardware_clearance`](@ref) for a [`HardwareKeepOut`](@ref), and `Inf` for a beam to
which the keep-out does not apply. It is a pure function of the beams, or of the setup
`beams` of an `Evaluation`, in their number type. A constraint
`keepout_clearance(keepout, ev) .≥ minimum` can therefore be handed to an optimizer as
it stands.
"""
keepout_clearance(k::ApertureKeepOut, beams::AbstractVector{<:Beam}) = aperture_clearance(k, beams)
keepout_clearance(k::SphereKeepOut, beams::AbstractVector{<:Beam}) = sphere_clearance(k, beams)
keepout_clearance(k, ev::Evaluation) = keepout_clearance(k, ev.beams)

"""
    keepout_penalty(keepout, beams) -> Real
    keepout_penalty(keepout, ev::Evaluation) -> Real

The penalty that a keep-out adds to an objective, `weight × Σᵢ max(0, minimum − cᵢ)²`
over the clearances `cᵢ = keepout_clearance(keepout, beams)`, with the `weight` and
`minimum` of the keep-out. For a tuple or vector of keep-outs it is the sum of their
penalties, and for `nothing` it is 0. The `keepout` keyword of
[`symmetry_objective`](@ref), [`asymmetry_objective`](@ref) and
[`face_flux_objective`](@ref) adds this penalty, so a hand-written objective adds
`keepout_penalty(keepout, ev)` to the same effect. In the number type of the beams
(ForwardDiff duals) and thread-safe. A keep-out of another type takes part through a
`keepout_clearance` method and `weight` and `minimum` fields.

A clearance of `+Inf`, for a beam to which the keep-out does not apply, adds nothing. A
clearance of `-Inf` or `NaN` throws an `ArgumentError` naming the keep-out and the beam.
`-Inf` marks a beam that fails a rule with no finite margin to continue, such as a
[`CaptureKeepIn`](@ref) beam whose axis misses every part when no part has a free edge.
Either value would give the optimizer a non-finite objective, and a large finite
penalty in its place would hide the cause.
"""
keepout_penalty(::Nothing, beams::AbstractVector{<:Beam}) = 0.0
keepout_penalty(k, ev::Evaluation) = keepout_penalty(k, ev.beams)

# The penalties of several keep-outs (a tuple or vector) are added. A tuple is unrolled
# with `map`, since a generator over a heterogeneous tuple made inference in Julia 1.13
# take tens of seconds.
keepout_penalty(ks::Tuple, beams::AbstractVector{<:Beam}) =
    sum(map(k -> keepout_penalty(k, beams), ks); init = 0.0)

function keepout_penalty(ks::AbstractVector, beams::AbstractVector{<:Beam})
    J = 0.0
    for k in ks
        J += keepout_penalty(k, beams)
    end
    return J
end

function keepout_penalty(k, beams::AbstractVector{<:Beam})
    c = keepout_clearance(k, beams)
    _check_clearances(k, beams, c)
    # The penalty is in the number type of the beams also when no clearance is finite. A
    # Float64 0.0 beside dual penalties made the sum over keep-outs a runtime dispatch.
    J = k.weight * sum(max(zero(ci), k.minimum - ci)^2 for ci in c if ci != Inf; init = 0.0)
    return convert(eltype(c), J)
end

# Checks the clearances `c` of `beams` from a keep-out for values that neither a penalty
# nor a constraint can take. `-Inf` and `NaN` are an `ArgumentError` naming the keep-out
# and the first such beam. `+Inf`, for a beam to which the keep-out does not apply, passes.
function _check_clearances(k, beams, c)
    for i in eachindex(c, beams)
        ci = ForwardDiff.value(c[i])
        (isnan(ci) || ci == -Inf) && _clearance_error(k, beams[i], ci)
    end
    return nothing
end

@noinline function _clearance_error(k, b::Beam, ci)
    names = hasfield(typeof(k), :names) ? getfield(k, :names) : nothing
    what = names isa AbstractVector ? string(nameof(typeof(k)), " (", join(names, ", "), ")") :
           string(nameof(typeof(k)))
    why = isnan(ci) ? "NaN" :
          "-Inf (the beam fails the rule with no finite margin: for a CaptureKeepIn, its axis misses " *
          "every part and no part has a free edge to continue the margin from)"
    throw(ArgumentError("keep-out $what: the clearance of beam $(b.name) is $why"))
end
