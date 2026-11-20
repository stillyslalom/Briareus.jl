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

# d(rim point)/dτ
_rim_tangent(r::_RimPiece, τ) =
    r.arc ? r.ρ * r.Δθ * (-sin(r.θ0 + τ * r.Δθ) * r.e1 + cos(r.θ0 + τ * r.Δθ) * r.e2) : r.q - r.p

# A hole's rim samples, piece after piece: the rim point, its derivative in τ, τ, the
# piece, and for an arc sample after the first of its piece the arc length from the
# sample before (it bounds the change of s over the interval) and the sagitta of that
# arc (-Inf and 0 on an edge).
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
least distance on each sample interval is found exactly: at an end, or where the
distance's slope along the rim rises through zero, located from the slopes at the
samples to 1e-12 of the interval. An interval is first split where the distance has a
kink (where the rim crosses the plane of best focus, the f-cone term turning), where
the rim crosses an axis of an elongated envelope (the first-order distance below turns
sharply there), inside the envelope also where the rim passes nearest the beam axis
(there the first-order distance dips steeply), and where the cubic through the ends'
distances and slopes turns twice. A dual beam's rim is searched on its values; the
result's partials follow from its sensitivities to the focus and the envelope's
semi-axes, so their cost does not grow with the number of partials.

The clearance is a smooth minimum of these interval minima, `−softness · log Σ σₖ
exp(−dₖ/softness)` (`softness` in mm), each interval's weight σₖ 1 within 10 `softness`
of the least and falling smoothly to 0 at 20 `softness` (a term there weighs at most
e⁻¹⁰ of the least). The exact minimum over the rim has kinks where the worst rim point
jumps between separate local minima of the distance, which the smooth minimum rounds
off. It never exceeds the exact minimum, lies below it by at most `softness` times the
log of the number of intervals counted, and is the exact minimum for `softness = 0`.
The default, 1e-4 mm, is an engineering default well below the keep-out's agreement
with `clearance` (a few micrometers): on NIF quads through a hohlraum's entrance holes
the smooth minimum lies 0 to 0.18 µm below the exact one, and over 144 scans of the aim
the largest change of slope in a step falls from 0.21 (a kink) to 0.002 (what remains
where one interval holds two local minima; research notes, "Defaults changed
(2026-11-19)").

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
    samples::Vector{_RimSamples}
    peak_fraction::Union{Nothing,Float64}
    intensity::Union{Nothing,Float64}
    transverse::Float64
    focus::Float64
    margin::Float64
    minimum::Float64
    weight::Float64
    softness::Float64      # length (mm) of the smooth minimum over the rim (0: the exact one)
end

function ApertureKeepOut(t::Target; apertures = nothing, peak_fraction = nothing,
                         intensity = nothing, pointing = 0.05, target_offset = 0.025,
                         hardware_tolerance = 0.0, focus = 1.0, margin = 0.0,
                         minimum = 0.0, weight = 100.0, nrim::Integer = 96, softness::Real = 1e-4)
    intensity === nothing && peak_fraction === nothing && (peak_fraction = 0.01)
    nrim >= 3 || throw(ArgumentError("nrim must be at least 3"))
    softness >= 0 || throw(ArgumentError("softness must be non-negative"))
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
    return ApertureKeepOut(names, centers, rims, _RimSamples.(rims), peak_fraction, intensity,
                           pointing + target_offset + hardware_tolerance, focus, margin,
                           minimum, weight, Float64(softness))
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
function aperture_clearance(k::ApertureKeepOut, beams::AbstractVector{<:Beam})
    # the envelope radius once per distinct spot and power (values): a root search each
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

# `Inf` in the beam's number type: a Float64 `Inf` beside dual clearances made the
# clearance a `Union` and the objectives' keep-out penalty infer as `Any`.
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
    fv = _RimFrame(k, b isa Beam{Float64} ? b : _value_beam(b), ρ)
    # the rim is searched on the beam's values; the terms found enter in its number type
    terms, imin = _rim_terms(k, rim, k.samples[j], fv)
    return _rim_result(k, terms, imin, fv, b, ρ) - k.transverse - k.focus * fv.cone
end

# What the distance from a rim point to a beam's envelope reads of the beam (in its number
# type): the envelope's semi-axes at best focus, the beam frame, the focus and the f-cone.
struct _RimFrame{T}
    A::T
    B::T
    kv::Vec3
    e1::Vec3
    e2::Vec3
    focus::SVector{3,T}
    cone::Float64
    margin::Float64
end

# (`ρ`: the envelope radius, read from the beam's values)
function _RimFrame(k::ApertureKeepOut, b::Beam{T},
                   ρ = envelope_radius(b.spot, b.power; intensity = k.intensity, peak_fraction = k.peak_fraction)) where {T}
    kv = direction(b)
    e1, e2 = spot_axes(b)
    return _RimFrame{T}(ρ * b.spot.a, ρ * b.spot.b, kv, e1, e2, b.aim + b.defocus * kv, 1 / (2b.fnumber), k.margin)
end

# The distance from the rim point `x` to the envelope.
function _rim_distance(f::_RimFrame, x::Vec3)
    d = x - f.focus
    s, u, v = dot(d, f.kv), dot(d, f.e1), dot(d, f.e2)
    return _sampson_distance(u, v, f.A, f.B) - f.cone * abs(s) - f.margin
end
_rim_distance(f::_RimFrame, r::_RimPiece, τ::Float64) = _rim_distance(f, _rim_point(r, τ))

# The smooth minimum counts the sample intervals whose least distance lies within
# `_SOFT_FULL` softness of the least at full weight, and those up to `_SOFT_ZERO` at a
# weight falling (C¹) to 0 there: where a weight falls, its term is at most e⁻¹⁰ = 4.5e-5
# of the least's.
const _SOFT_FULL = 10.0
const _SOFT_ZERO = 20.0

# The smooth minimum's terms, found on the beam's values (frame `f`): for each sample
# interval of the rim its least distance (at an end, or at a minimum inside), kept with its
# rim point if within `_SOFT_ZERO` softness of the least (`softness = 0`: the least
# alone). Returns the terms (distance, point) and the index of the least.
function _rim_terms(k::ApertureKeepOut, rim::Vector{_RimPiece}, rs::_RimSamples, f::_RimFrame{Float64})
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
                # (the terms kept before the least was reached, dropped now and then)
                length(terms) == cap && (_drop_above!(terms, lo + reach); cap = max(32, 2length(terms)))
            end
        end
        x0, d0, S0, s0, t0, w0 = x1, d1, S1, s1, t1, w1
    end
    _drop_above!(terms, lo + reach)
    # a NaN distance (an envelope of zero size) passes into the result
    isempty(terms) && push!(terms, (NaN, rs.x[1]))
    imin = 1
    for i in eachindex(terms)
        terms[i][1] < terms[imin][1] && (imin = i)
    end
    return terms, imin
end

_drop_above!(terms, cut) = filter!(t -> t[1] <= cut, terms)

# The Float64 frame of the search, with the semi-axes' inverse squares.
struct _SlopeFrame
    focus::Vec3
    kv::Vec3
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
function _SlopeFrame(f::_RimFrame{Float64})
    e = max(f.A, f.B) / min(f.A, f.B)
    return _SlopeFrame(f.focus, f.kv, f.e1, f.e2, inv(f.A), inv(f.B), inv(f.A^2), inv(f.B^2), min(f.A, f.B), f.cone,
                       f.margin, e > 1 ? 0.25 / sqrt(e^2 - 1) : Inf)
end

# At the rim point `x` with tangent `t` (dx/dτ): the distance to the envelope, the slope of
# its Sampson part in τ, the axial position s and its slope, and the point in the
# envelope's units (u/A, v/B: inside it within 1 of the axis). The distance's slope is
# S′ − cone·sign(s)·s′: |s| has a concave kink where s = 0.
@inline function _rim_slope(f::_SlopeFrame, x::Vec3, t::Vec3)
    d = x - f.focus
    s, u, v = dot(d, f.kv), dot(d, f.e1), dot(d, f.e2)
    ts, tu, tv = dot(t, f.kv), dot(t, f.e1), dot(t, f.e2)
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

# The point and tangent of piece `r` at τ (one sincos on an arc).
function _rim_point_tangent(r::_RimPiece, τ::Float64)
    r.arc || return r.p + τ * (r.q - r.p), r.q - r.p
    sn, cs = sincos(r.θ0 + τ * r.Δθ)
    return r.o + r.ρ * (cs * r.e1 + sn * r.e2), r.ρ * r.Δθ * (-sn * r.e1 + cs * r.e2)
end

# The distance's slope at τ on piece `r`, on the branch `σ` = sign(s) of |s|.
function _slope_at(f::_SlopeFrame, r::_RimPiece, τ::Float64, σ)
    x, t = _rim_point_tangent(r, τ)
    _, S, _, ts = _rim_slope(f, x, t)
    return S - f.cone * σ * ts
end

# The least distance on the sample interval [τa, τb] of piece `p` of `rim` (`arc`: its arc
# length, -Inf on an edge; `sag`: the arc's sagitta), given the lesser end's (`g` at `xg`)
# and the ends' Sampson slope, s, s′, point in envelope units and distance (`ea`, `eb`).
# The interval is split (`_split_min`) where |s| has a kink; where the rim crosses an
# axis of the envelope, if the Sampson distance's direction factor can turn sharply there
# (an elongated spot: it peaks over about 1/√(e² − 1) of the direction in envelope units,
# e the ratio of the semi-axes; the interval split if it turns through a quarter of
# that); and, if it reaches inside the envelope, also where it passes nearest the beam
# axis. Each part on which the distance is smooth is searched by `_part_min`.
@inline function _interval_min(f::_SlopeFrame, rim, p, τa, τb, arc, sag, ea, eb, g, xg)
    Sa, sa, ta, wa, da = ea
    Sb, sb, tb, wb, db = eb
    c = f.cone
    # a kink inside: s changes sign, or (on an arc) passes an extremum near enough zero to
    # cross it (|s′| is at most the interval's arc length per unit τ)
    kink = c > 0 && ((sa > 0 && sb < 0) || (sa < 0 && sb > 0) || (ta * tb < 0 && min(abs(sa), abs(sb)) <= arc))
    # the chord's nearest approach to the beam axis, in envelope units (with the arc's sag)
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
# `dr` and slopes `ml`, `mr` of one sign turns twice inside (its slope, a quadratic, changes
# sign at its vertex).
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

# The least distance below `g` (at `xg`) on a part [l, rr] of piece `p` of `rim` on which
# the distance is smooth (`σ`: the branch of |s|), from its ends' distances `dl`, `dr` and
# slopes `ml`, `mr`: a minimum inside where the slope rises through zero (`_slope_root`).
# Slopes of one sign at both ends whose cubic Hermite interpolant turns twice inside (a
# minimum and a maximum between the samples) split the part at its two turning points,
# each a candidate, and each piece is searched again (once).
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

# `_interval_min` on an interval split at the kinks of |s| and, with `axes`, where the rim
# crosses the envelope's axes (u = 0, v = 0: the Sampson distance's direction factor at
# its extremes, between which it turns within the axis's distance, sharply for an
# elongated spot), and at `τq`, the chord's nearest approach to the beam axis inside the
# envelope (where the Sampson distance (q − 1)·h falls to −h); each split point is a
# candidate, and each part is searched by `_part_min` on its branch of |s|.
function _split_min(f::_SlopeFrame, r::_RimPiece, τa, τb, ea, eb, g, xg, kink, axes, τq)
    k1, k2 = kink ? _crossings(f, r, τa, τb, f.kv) : (NaN, NaN)
    u1, u2 = axes ? _crossings(f, r, τa, τb, f.e1) : (NaN, NaN)
    v1, v2 = axes ? _crossings(f, r, τa, τb, f.e2) : (NaN, NaN)
    m = τa < τq < τb ? τq : NaN
    splits = sort(SVector(k1, k2, u1, u2, v1, v2, m, τb); by = τ -> isnan(τ) ? Inf : τ)
    c = f.cone
    rim = (r,)
    l, dl, Sl, tl = τa, ea[5], ea[1], ea[3]  # the part's left end: τ, distance, Sampson slope, s′
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
        σ = kink ? sign(dot(_rim_point(r, (l + rr) / 2) - f.focus, f.kv)) : ea[2] != 0 ? sign(ea[2]) : sign(eb[2])
        g, xg = _part_min(f, rim, 1, l, rr, dl, dr, Sl - c * σ * tl, Sr - c * σ * tr, σ, g, xg)
        rr == τb && break
        l, dl, Sl, tl = rr, dr, Sr, tr
    end
    return g, xg
end

# The points inside (τa, τb) at which the rim piece `r` crosses the plane through the focus
# normal to `n` (dot(x − focus, n) = 0), ascending (NaN for none): linear in τ on an edge,
# C + M cos(θ − φ) on an arc. With `n` the beam direction, the kinks of |s|.
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

# The τ in (τa, τb) at which arc `r` reaches the angle θ (modulo a turn, `period` in τ),
# or NaN.
function _arc_tau(r::_RimPiece, θ, τa, τb, period)
    τ = τa + mod((θ - r.θ0) / r.Δθ - τa, period)
    return τa < τ < τb ? τ : NaN
end

# The point in (a, b) where the slope `S′ − cone·σ·s′` of the distance along piece `r`
# rises through zero (`la` < 0 at a, `lb` > 0 at b): regula falsi with the Illinois step,
# to 1e-12 of the interval.
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

# The smooth minimum of the terms' distances (`ApertureKeepOut`), evaluated in the frame
# `f`: −λ log Σₖ σₖ exp(−dₖ/λ), each term's weight σₖ 1 within `_SOFT_FULL` λ of the
# least and falling (C¹) to 0 at `_SOFT_ZERO` λ.
function _rim_soft_min(k::ApertureKeepOut, terms, imin, f::_RimFrame)
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

# Dual numbers of the keep-out: the smooth minimum's sensitivities to the focus (three)
# and the envelope's semi-axes (two), the only quantities of the beam it reads.
struct _RimTag end
const _RimDual = ForwardDiff.Dual{_RimTag,Float64,5}

# The smooth minimum in the beam's number type. On Float64 the terms on the values; on a
# dual beam its sensitivities to the frame's five dual quantities (`_RimDual`, Float64
# partials), carried to the beam's partials by the chain rule, so the cost does not grow
# with the beam's partials; on nested duals the terms in the number type of the beam.
_rim_result(k::ApertureKeepOut, terms, imin, fv, ::Beam{Float64}, ρ) = _rim_soft_min(k, terms, imin, fv)

function _rim_result(k::ApertureKeepOut, terms, imin, fv, b::Beam{<:ForwardDiff.Dual{<:Any,Float64}}, ρ)
    e(i) = ForwardDiff.Partials{5,Float64}(ntuple(j -> Float64(j == i), 5))
    fs = _RimFrame{_RimDual}(_RimDual(fv.A, e(4)), _RimDual(fv.B, e(5)), fv.kv, fv.e1, fv.e2,
                             SVector(_RimDual(fv.focus[1], e(1)), _RimDual(fv.focus[2], e(2)), _RimDual(fv.focus[3], e(3))),
                             fv.cone, fv.margin)
    c = _rim_soft_min(k, terms, imin, fs)
    p = ForwardDiff.partials(c)
    fb = _RimFrame(k, b, ρ)
    Δ(x) = x - ForwardDiff.value(x)          # the partials alone, at value 0
    return ForwardDiff.value(c) + p[1] * Δ(fb.focus[1]) + p[2] * Δ(fb.focus[2]) + p[3] * Δ(fb.focus[3]) +
           p[4] * Δ(fb.A) + p[5] * Δ(fb.B)
end

_rim_result(k::ApertureKeepOut, terms, imin, fv, b::Beam, ρ) = _rim_soft_min(k, terms, imin, _RimFrame(k, b, ρ))

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
fitted with `curved_components`, its fill tube a separate piece), not the small
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

"""
    keepout_clearance(keepout, beams) -> Vector
    keepout_clearance(keepout, ev::Evaluation) -> Vector

Each beam's worst-case clearance (mm) from a keep-out, the quantity its penalty holds
at or above the keep-out's `minimum`: [`aperture_clearance`](@ref) for an
[`ApertureKeepOut`](@ref), [`sphere_clearance`](@ref) for a [`SphereKeepOut`](@ref),
[`hardware_clearance`](@ref) for a [`HardwareKeepOut`](@ref); `Inf` for a beam the
keep-out does not apply to. A pure function of the beams (of the
evaluation's setup `beams` for an `Evaluation`), in their number type, so a constraint
`keepout_clearance(keepout, ev) .≥ minimum` can be handed to an optimizer as it stands.
"""
keepout_clearance(k::ApertureKeepOut, beams::AbstractVector{<:Beam}) = aperture_clearance(k, beams)
keepout_clearance(k::SphereKeepOut, beams::AbstractVector{<:Beam}) = sphere_clearance(k, beams)
keepout_clearance(k, ev::Evaluation) = keepout_clearance(k, ev.beams)

"""
    keepout_penalty(keepout, beams) -> Real
    keepout_penalty(keepout, ev::Evaluation) -> Real

The penalty a keep-out adds to an objective: `weight × Σᵢ max(0, minimum − cᵢ)²` over
the clearances `cᵢ = keepout_clearance(keepout, beams)`, with the keep-out
`weight` and `minimum`; for a tuple or vector of keep-outs the sum of theirs, and 0 for
`nothing`. It is what the `keepout` keyword of [`symmetry_objective`](@ref),
[`asymmetry_objective`](@ref) and [`face_flux_objective`](@ref) adds, so a hand-written
objective adds `keepout_penalty(keepout, ev)` to the same effect. In the beams' number
type (ForwardDiff duals) and thread-safe. A keep-out of another type takes part through
a `keepout_clearance` method and `weight` and `minimum` fields.

A clearance of `+Inf` (a beam the keep-out does not apply to) adds nothing. A clearance
of `-Inf` (a beam that fails a rule with no finite margin to continue, as a
[`CaptureKeepIn`](@ref) beam whose axis misses every part when no part has a free edge)
or `NaN` throws an `ArgumentError` naming the keep-out and the beam: either would hand
the optimizer a non-finite objective, and a large finite penalty in its place would
hide the cause.
"""
keepout_penalty(::Nothing, beams::AbstractVector{<:Beam}) = 0.0
keepout_penalty(k, ev::Evaluation) = keepout_penalty(k, ev.beams)

# Several keep-outs (a tuple or vector) add their penalties. A tuple is unrolled with
# `map`: a generator over a heterogeneous tuple took Julia 1.13's inference 27 s.
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
    # in the beams' number type also when no clearance is finite (a Float64 0.0 beside
    # dual penalties made the sum over keep-outs a runtime dispatch)
    J = k.weight * sum(max(zero(ci), k.minimum - ci)^2 for ci in c if ci != Inf; init = 0.0)
    return convert(eltype(c), J)
end

# A keep-out's clearances `c` of `beams`, checked for what neither a penalty nor a
# constraint can take: `-Inf` and `NaN` are an `ArgumentError` naming the keep-out and the
# first such beam (`+Inf`, a beam the keep-out does not apply to, passes).
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
