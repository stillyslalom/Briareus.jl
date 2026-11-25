# This file defines capture and clearance rules as quantities of an evaluation, for
# objectives and constraints written by hand. Each rule is a pure function of the
# evaluation or of its beams, generic in the number type and thread-safe. Each has a
# natural limit, such as a ceiling or a floor on a fraction, or a minimum clearance.

"""
    direct_fraction(ev) -> Real
    direct_fraction(ev, names) -> Real

Fraction of the delivered laser power absorbed by the observed elements of the scenario (the
capsule whose drive is measured), or with `names` (a collection of part or component
names) by those parts instead. The `direct` penalty of [`symmetry_objective`](@ref),
[`asymmetry_objective`](@ref) and [`face_flux_objective`](@ref) weighs the square of this
fraction, and the `direct_names` of those objectives corresponds to `names`. A name that
is not in the target is an `ArgumentError`. The fraction is relative to the effective
(post-LPI) power, as every deposition fraction is; see [`clipped_fraction`](@ref) and
[`missed_fraction`](@ref) for the other two penalized fractions.
"""
direct_fraction(ev::Evaluation, names = nothing) = _direct_fraction(ev, names)

"""
    clipped_fraction(ev, names...) -> Real

Fraction of the delivered laser power striking the outside (back side) of the named
parts or components (a part stands for all its pieces, see [`part_pieces`](@ref)); with
no names, of every component. Beside [`power_fraction`](@ref)`(ev, names...)`, the power
absorbed on their elements, the sum

    power_fraction(ev, names...) + clipped_fraction(ev, names...)

is all the laser light the parts receive, from either side: the quantity of a ceiling
on hardware meant to stay dark, such as closed windows, shields and stalks (at most the
value of an accepted design), or of a floor on capture by a part meant to be hit. Both are
totals over the beams, in the number type of the evaluation (differentiable), and relative
to the effective (post-LPI) power, as [`component_power`](@ref) reports the same powers
in TW (`name` and `Symbol(name, "_outside")`). An open window is an aperture: light
passes through it, so it receives none.
"""
clipped_fraction(ev::Evaluation, names::Symbol...) =
    sum(ev.deposition.backside[ci] for n in names for ci in part_pieces(ev.scenario.target, n);
        init = zero(ev.deposition.delivered)) / ev.deposition.delivered

# ---------------------------------------------------------------------------------------
# Keep-out for named hardware
# ---------------------------------------------------------------------------------------

# A cluster is a group of the triangles of a keep-out, with a bounding sphere, the
# triangles and their edges.
struct _TriangleCluster
    center::Vec3
    radius::Float64
    tris::Vector{NTuple{3,Vec3}}
    edges::Vector{NTuple{2,Vec3}}
end

"""
    HardwareKeepOut(target, parts; targets = every other opaque surface,
                    peak_fraction = 0.01, intensity = nothing, pointing = 0.05,
                    target_offset = 0.025, hardware_tolerance = 0.0, focus = 1.0,
                    margin = 0.0, tolerance = 0.0, pivot = :tcc, miss = :aim,
                    minimum = 0.0)

Keep-out constraint for named hardware such as shields, stalks, diagnostic snouts and closed
windows. [`hardware_clearance`](@ref)`(keepout, beams)` gives the worst-case clearance (mm)
of each beam from `parts`, judged as [`clearance`](@ref) judges it with `wings = []`.
`parts` is a name or a collection of part or component names, and a part stands for all its
pieces. The beam envelope is the best-focus contour at `peak_fraction` of the peak, or at
`intensity` (W/cm²), grown by the f-cone radius `|s|/(2F)` and `margin`. The clearance is
the smallest signed distance from the surface of the parts to the envelope, transverse to
the beam and negative inside it, less the transverse error radius
`pointing + target_offset + hardware_tolerance` and `focus/(2F)` for the focus error.

Only the surface upstream of the landing of the beam counts. As in `clearance`, the landing
is the plane normal to the beam where its nominal axis first meets one of `targets`, by
default every opaque surface not among `parts`. A beam whose axis meets none of them is
judged up to its aim plane, or along its whole length with `miss = :all`. A beam whose
landing is on one of `parts` itself (when `targets` names it) gets `Inf`. The defaults are
those of [`ApertureKeepOut`](@ref). They are engineering defaults and are not set by a
facility criterion; the 1% contour is the level of the TaLIS 3ω clearance rule.

`tolerance` (degrees) is the rule tolerance of [`talis_check`](@ref): the parts may sit
rotated off nominal by up to that angle about `pivot` (`:tcc`, `:aim` the beam aim
point, `:max` the farther of the two, or a point), which costs each point `r·tan α`, `r`
its distance from the pivot, subtracted at each point inside the minimization (so the
clearance stays continuous where the closest point jumps). With `targets` the capture
components, `miss = :all`, `tolerance = 0.5` and `minimum = 0.5` (all other defaults
equal) the keep-out is the TaLIS 3ω clearance rule of `talis_check` for each beam:

    HardwareKeepOut(t, hardware; targets = capture, tolerance = 0.5, miss = :all,
                    minimum = 0.5)

a surrogate checked against `talis_check` in this package (its `:clearance_3w` worst
values; see the research notes), not against the facility tool. `talis_check` reports
each hardware part on a separate row, where the keep-out takes the minimum over its parts
(one keep-out per part gives the rows), and it passes a clearance strictly above
0.5 mm; on an `leh` component its failure is a review, which a keep-out does not
distinguish.

The distance to the envelope ellipse is exact (the point–ellipse distance `clearance`
uses), and the minimum over the parts is found, not sampled. On a triangle not crossing
best focus the distance is convex (a convex function of the transverse position less a
linear one of the axial), and it has no stationary point inside the triangle except on
the medial strip, where it is not differentiable: the segment of the major axis between
`±(A² − B²)/A` (the axis itself for a round spot) swept along the beam, at least
`B²/A` inside the envelope. So the minimum lies on the triangle edges, their crossings
with the medial strip, or their cuts by the landing plane, each split at best focus:
line segments on which the distance is convex. The triangles are grouped in spatially
compact clusters of up to 16; clusters and segments are taken best first by a lower
bound (the circumscribed circle and the homothetic ellipse through the center of a cluster
or the midpoint of a segment, less the Lipschitz constant `√(1 + 1/(2F)²)` times its radius
or half length), and each segment is minimized by a safeguarded root search on the
derivative, until the bound reaches the smallest value found. The result is the exact
minimum to about 1e-12 mm.

The tolerance term `−r·tan α` is concave in the point. With it the minimum still lies off
the triangle interiors, because off the medial strip the Hessian of the distance has rank
one, so in every plane the sum curves down along some direction. The term does spoil the
convexity along a segment. Each segment is therefore minimized with the term replaced by
its chord, a lower bound
equal to it at the ends, and split at the minimizer of the chord until the two agree to
1e-14 mm; the bounds of clusters and segments subtract `tan α (r + radius)`.

The clearance is continuous in the aim, defocus and spot size of the beams (a minimum over a
fixed surface), and its derivatives are those of the distance at the closest point
(Danskin's theorem); on the medial strip, where two points of the ellipse are equally near,
they are the mixture of theirs that leaves the distance stationary at the minimum, as at the
saddle of a minimax. They are one-sided only where two separate points of the parts are
equally close (a kink of the minimum itself). The landing is found on the nominal beam: the
clearance jumps where the landing of the axis moves from one component to another
while the parts reach past the landing plane, and the cut by the landing plane does
not move with the beam in the derivatives.

The keep-out departs from `clearance` in two ways, both conservative. Every point of the
parts upstream of the landing counts, also where opaque surfaces hide it from the beam;
`clearance` skips points the beam cannot reach, and `talis_check` judges 3ω clearance as
this keep-out does. The clearance-only far-field wings are not modeled. Pass the keep-out
as `keepout` to the objectives, by itself or in a tuple with other keep-outs; the
objectives then add
`weight × Σᵢ max(0, minimum − cᵢ)²` ([`keepout_penalty`](@ref)).
"""
struct HardwareKeepOut <: AbstractKeepOut
    names::Vector{Symbol}
    target::Target
    parts::Set{Int}                      # the parts' components
    others::Set{Int}                     # components a beam can land on
    clusters::Vector{_TriangleCluster}   # the parts' triangles
    peak_fraction::Union{Nothing,Float64}
    intensity::Union{Nothing,Float64}
    transverse::Float64
    focus::Float64
    margin::Float64
    tolerance::_TaLISTol                 # the rule tolerance r·tan α
    through::Bool                        # an unlanded beam judged along its whole length
    minimum::Float64
end

function HardwareKeepOut(t::Target, parts; targets = nothing, peak_fraction = nothing,
                         intensity = nothing, pointing = 0.05, target_offset = 0.025,
                         hardware_tolerance = 0.0, focus = 1.0, margin = 0.0, tolerance = 0.0,
                         pivot = :tcc, miss::Symbol = :aim, minimum = 0.0)
    intensity === nothing && peak_fraction === nothing && (peak_fraction = 0.01)
    miss in (:aim, :all) || throw(ArgumentError("miss must be :aim or :all"))
    names = _part_names(parts)
    idx = Set(ci for n in names for ci in part_pieces(t, n))
    others = targets === nothing ?
             Set(i for (i, c) in enumerate(t.components) if c.kind === :surface && !(i in idx)) :
             Set(ci for n in targets for ci in part_pieces(t, Symbol(n)))
    tris = _part_triangles(t, idx)
    isempty(tris) && throw(ArgumentError("parts $(join(repr.(names), ", ")) have no triangles"))
    return HardwareKeepOut(names, t, idx, others, _clusters(tris), peak_fraction, intensity,
                           pointing + target_offset + hardware_tolerance, focus, margin,
                           _talis_tol(tolerance, pivot), miss === :all, minimum)
end

function _part_names(parts)
    names = parts isa Union{Symbol,AbstractString} ? [Symbol(parts)] : Symbol[Symbol(p) for p in parts]
    isempty(names) && throw(ArgumentError("name at least one part"))
    return names
end

# Return the triangles of the components `idx` as vertex triples.
function _part_triangles(t::Target, idx)
    tris = NTuple{3,Vec3}[]
    for k in eachindex(t.tri_area)
        t.tri_comp[k] in idx || continue
        v0 = t.tri_v0[k]
        push!(tris, (v0, v0 + t.tri_e1[k], v0 + t.tri_e2[k]))
    end
    return tris
end

# Group triangles, or segments, into spatially compact clusters of at most `leaf`. A group
# that is too large is halved at the median centroid along its longest extent. Each
# cluster has a bounding sphere and its edges, each edge listed once; for segments the
# edges are the segments themselves.
function _clusters(tris::Vector{NTuple{N,Vec3}}; leaf::Int = 16) where {N}
    out = _TriangleCluster[]
    isempty(tris) && return out
    key(a, b) = Tuple(a) < Tuple(b) ? (a, b) : (b, a)
    stack = [tris]
    while !isempty(stack)
        ts = pop!(stack)
        lo = reduce((x, y) -> min.(x, y), (v for tri in ts for v in tri))
        hi = reduce((x, y) -> max.(x, y), (v for tri in ts for v in tri))
        if length(ts) <= leaf
            c = (lo + hi) / 2
            r = maximum(norm(v - c) for tri in ts for v in tri)
            push!(out, N == 2 ? _TriangleCluster(c, r, NTuple{3,Vec3}[], ts) :
                       _TriangleCluster(c, r, ts, unique!([key(x, y) for (a, b, d) in ts
                                                           for (x, y) in ((a, b), (b, d), (d, a))])))
            continue
        end
        ax = argmax(hi - lo)
        ts = sort(ts; by = tri -> sum(v[ax] for v in tri))
        h = length(ts) ÷ 2
        push!(stack, ts[1:h], ts[h+1:end])
    end
    return out
end

"""
    hardware_clearance(keepout, beams) -> Vector

Worst-case clearance (mm) of each beam from the parts of a [`HardwareKeepOut`](@ref),
the rule tolerance included; `Inf` for a beam whose landing lies upstream of all of them
or on one of them. The clearance is differentiable in the aim, defocus and spot size of
the beams and has their number type.
"""
hardware_clearance(k::HardwareKeepOut, beams::AbstractVector{<:Beam}) =
    [_hardware_clearance(k, b) for b in beams]

keepout_clearance(k::HardwareKeepOut, beams::AbstractVector{<:Beam}) = hardware_clearance(k, beams)

# Return the envelope of a beam on values at the contour of rule `k`, which sets the
# `intensity` or `peak_fraction` and the `margin`. Also return the contour radius in units
# of the spot semi-axes and the aim point on values.
function _envelope_frame(k, b::Beam)
    bv = _value_beam(b)
    ρ = envelope_radius(bv.spot, bv.power; intensity = k.intensity, peak_fraction = k.peak_fraction)
    kv = direction(b)
    e1, e2 = spot_axes(b)
    f0 = Vec3(bv.aim) + bv.defocus * kv
    return _BeamConeFrame{Float64}(f0, kv, e1, e2, ρ * bv.spot.a, ρ * bv.spot.b, 1 / (2b.fnumber), k.margin), ρ,
           Vec3(bv.aim)
end

# Return the axis of an envelope frame as a `BeamEnvelope`, which is all that
# `_landing_hit` uses.
_axis_envelope(ef::_BeamConeFrame{Float64}, aim::Vec3) = BeamEnvelope(ef.focus, aim, ef.k, ef.e1, ef.e2, ef.A, ef.B, ef.cone, 0.0)

function _hardware_clearance(k::HardwareKeepOut, b::Beam)
    ef, ρ, aim = _envelope_frame(k, b)
    s_stop, hit = _landing_hit(k.target, _axis_envelope(ef, aim), k.others; miss = k.through ? Inf : nothing)
    hit in k.parts && return _no_clearance(b)
    return _envelope_clearance(k, b, ef, ρ, aim, k.clusters, s_stop)
end

# Return the worst-case clearance of beam `b` from the triangles and segments of the
# clusters `cls` upstream of `s_stop`, less the error set and rule tolerance of rule `k`.
# The envelope of the beam is `ef` at the contour radius `ρ`, and its aim on values is
# `aim`. The result is the distance at the point found on values, in the number type of
# the beam, and `Inf` if no point counts.
function _envelope_clearance(k, b::Beam, ef::_BeamConeFrame{Float64}, ρ, aim::Vec3, cls, s_stop)
    found, x, w = _hardware_closest(cls, k.tolerance, ef, s_stop, aim)
    found || return _no_clearance(b)
    τ = k.tolerance
    kv, e1, e2 = ef.k, ef.e1, ef.e2
    # The slope of the tolerance along `w` enters the balance on the medial strip.
    gw = iszero(τ.tanα) ? 0.0 : -dot(_tol_gradient(τ, x, aim), w)
    # Evaluate the distance at the closest point in the number type of the beams.
    δ = x - (b.aim + b.defocus * kv)
    δ0 = x - ef.focus
    u0, v0, σ = dot(δ0, e1), dot(δ0, e2), dot(δ0, kv) < 0 ? -1 : 1
    u, v, A, B = dot(δ, e1), dot(δ, e2), ρ * b.spot.a, ρ * b.spot.b
    X0, X1, onaxis = _ellipse_closest(u0, v0, ef.A, ef.B)
    n1, n2 = _ellipse_unit_normal(u0, v0, ef.A, ef.B, X0, X1)
    E = _ellipse_signed(u, v, A, B, ef.A, ef.B, X0, X1, n1, n2)
    # Mirror the closest point across the major axis.
    Y0, Y1 = ef.A >= ef.B ? (X0, -X1) : (-X0, X1)
    if onaxis || _medial_tie(u0, v0, ef.A, ef.B, X0, X1, Y0, Y1)
        # On the medial axis two closest points are equally near. The derivative is then
        # the mixture of their derivatives that leaves the distance stationary along `w`,
        # as at the saddle of a minimax. `w` is the direction of the segment, or the
        # direction across the medial strip in the plane of the triangle. Because the
        # mixture is stationary along `w`, its value at a point found next to the kink
        # equals the value at the kink.
        m1, m2 = _ellipse_unit_normal(u0, v0, ef.A, ef.B, Y0, Y1)
        h⁺ = dot(_envelope_gradient(ef, n1, n2, σ), w) + gw
        h⁻ = dot(_envelope_gradient(ef, m1, m2, σ), w) + gw
        λ = h⁺ == h⁻ ? 0.5 : clamp(h⁻ / (h⁻ - h⁺), 0.0, 1.0)
        E = λ * E + (1 - λ) * _ellipse_signed(u, v, A, B, ef.A, ef.B, Y0, Y1, m1, m2)
    end
    s = dot(δ, kv)
    c = E - ef.cone * σ * s - k.margin - k.transverse - k.focus * ef.cone
    return iszero(τ.tanα) ? c : c - _tol_value(τ, x, b.aim)
end

# Return the allowance r·tan α of the rule tolerance at `x`, in the number type of `aim`.
# The pivots `:aim` and `:max` move with the aim.
function _tol_value(τ::_TaLISTol, x::Vec3, aim::SVector{3,T}) where {T}
    p = τ.pivot === :tcc ? zero(aim) : τ.pivot === :aim ? aim :
        τ.pivot === :point ? SVector{3,T}(τ.point) :
        norm(x) >= norm(x - ForwardDiff.value.(aim)) ? zero(aim) : aim
    return τ.tanα * norm(x - p)
end

# `_tol_gradient` returns the gradient of the allowance in `x`, and `_tol_bound` an upper
# bound of the allowance on the ball of radius `ρ` about `x`; r is 1-Lipschitz in the point.
function _tol_gradient(τ::_TaLISTol, x::Vec3, aim::Vec3)
    d = τ.pivot === :tcc ? x : τ.pivot === :aim ? x - aim : τ.pivot === :point ? x - τ.point :
        norm(x) >= norm(x - aim) ? x : x - aim
    r = norm(d)
    return r > 0 ? (τ.tanα / r) * d : zero(Vec3)
end
_tol_bound(τ::_TaLISTol, x::Vec3, aim::Vec3, ρ) = τ(x, aim) + τ.tanα * ρ

# Find the point of the triangles and segments of the clusters upstream of `s_stop` that is
# closest to the envelope `ef`, less the rule tolerance `τ`, whose pivot `aim` is the aim of
# the beam. Return whether such a point exists, the point, and the direction along which
# the distance is stationary at the point. That direction is the segment direction, or the
# direction across the medial strip in the plane of a triangle (see `HardwareKeepOut`).
function _hardware_closest(clusters::Vector{_TriangleCluster}, τ::_TaLISTol, ef::_BeamConeFrame{Float64},
                           s_stop::Float64, aim::Vec3)
    L = sqrt(1 + ef.cone^2)
    # Search best first over the clusters, bounded by their spheres, and the segments of
    # the opened clusters. Each heap entry contains the bound, the ends, the balancing
    # direction, whether the entry is a medial-strip crossing, and the cluster.
    heap = Tuple{Float64,Vec3,Vec3,Vec3,Bool,Int}[]
    for (i, c) in enumerate(clusters)
        _axial(ef, c.center) - c.radius > s_stop && continue
        _heap_push!(heap, (_point_bound(ef, c.center) - L * c.radius - _tol_bound(τ, c.center, aim, c.radius),
                           c.center, c.center, c.center, false, i))
    end
    segs = Tuple{Float64,Vec3,Vec3,Vec3,Bool}[]
    best, xbest, wbest = Inf, zero(Vec3), zero(Vec3)
    while !isempty(heap)
        lb, p, q, w, strip, i = _heap_pop!(heap)
        lb >= best && break
        if i > 0
            empty!(segs)
            _cluster_segments!(segs, clusters[i], ef, s_stop)
            for (sl, sp, sq, sw, ss) in segs
                sl -= _tol_bound(τ, (sp + sq) / 2, aim, norm(sq - sp) / 2)
                sl < best && _heap_push!(heap, (sl, sp, sq, sw, ss, 0))
            end
        elseif iszero(τ.tanα)
            # The ends of a crossing lie on edges or cuts. The search on those edges and
            # cuts finds the ends, with the direction along which a minimum at an end is
            # stationary.
            v, x = _segment_min(ef, p, q, best, !strip)
            v < best && ((best, xbest, wbest) = (v, x, w))
        else
            # The tolerance term −r·tan α is concave. The distance less the chord of the
            # term is convex, bounds the true value from below and equals it at the ends.
            # Where the minimum of that bound lies inside the segment and below the best
            # by more than the gap at that point, the segment is split there, and both
            # halves go back on the heap.
            gp, gq = -τ(p, aim), -τ(q, aim)
            v, x = _segment_min(ef, p, q, best, !strip, gp, gq)
            isfinite(v) || continue
            Δ = q - p
            θ = clamp(dot(x - p, Δ) / dot(Δ, Δ), 0.0, 1.0)
            gap = -τ(x, aim) - (gp + θ * (gq - gp))          # the term less its chord, ≥ 0
            fx = v + gap
            fx < best && ((best, xbest, wbest) = (fx, x, w))
            if v < best - 1e-14 && gap > 1e-14 && norm(Δ) > 1e-9
                θ = clamp(θ, 0.1, 0.9)                       # each half at most 0.9 of it
                m = p + θ * Δ
                if θ != dot(x - p, Δ) / dot(Δ, Δ)
                    σ = _axial(ef, (p + q) / 2) < 0 ? -1 : 1
                    fm = first(_envelope_value_gradient(ef, m, σ)) - τ(m, aim)
                    fm < best && ((best, xbest, wbest) = (fm, m, w))
                end
                _heap_push!(heap, (v, p, m, w, strip, 0))
                _heap_push!(heap, (v, m, q, w, strip, 0))
            end
        end
    end
    return isfinite(best), xbest, wbest
end

# Push onto `segs` the segments of cluster `c` on which the minimum of the envelope
# distance over its triangles upstream of `s_stop` lies (see `HardwareKeepOut`). Each
# entry contains the bound, the ends, the balancing direction, and whether the segment is
# a medial-strip crossing.
function _cluster_segments!(segs, cl::_TriangleCluster, ef::_BeamConeFrame{Float64}, s_stop)
    for (p, q) in cl.edges
        _segment_pieces!(segs, ef, s_stop, p, q, q - p, false)
    end
    # The medial strip is the plane through the axis and the major axis, within ±tmax of
    # the axis.
    amaj, tmax = ef.A >= ef.B ? (ef.e1, (ef.A^2 - ef.B^2) / ef.A) : (ef.e2, (ef.B^2 - ef.A^2) / ef.B)
    nm = cross(ef.k, amaj)
    for (a, b, c) in cl.tris
        sa, sb, sc = _axial(ef, a), _axial(ef, b), _axial(ef, c)
        # Add the cut by the landing plane.
        if max(sa, sb, sc) > s_stop && min(sa, sb, sc) < s_stop
            p, q, n = _plane_cut((a, b, c), (sa - s_stop, sb - s_stop, sc - s_stop))
            n == 2 && _segment_pieces!(segs, ef, s_stop, p, q, q - p, false)
        end
        # Add the crossing with the medial strip, clipped to its width.
        h = (dot(nm, a - ef.focus), dot(nm, b - ef.focus), dot(nm, c - ef.focus))
        (minimum(h) > 0 || maximum(h) < 0) && continue
        p, q, n = _plane_cut((a, b, c), h)
        n == 0 && continue
        tp, tq = dot(amaj, p - ef.focus), dot(amaj, q - ef.focus)
        tp > tq && ((p, q, tp, tq) = (q, p, tq, tp))
        (tq < -tmax || tp > tmax) && continue
        if tq > tp
            p, q = p + (max(tp, -tmax) - tp) / (tq - tp) * (q - p), p + (min(tq, tmax) - tp) / (tq - tp) * (q - p)
        end
        N = cross(b - a, c - a)
        _segment_pieces!(segs, ef, s_stop, p, q, cross(N, cross(N, nm)), true)
    end
    return segs
end

_axial(ef::_BeamConeFrame{Float64}, x) = dot(x - ef.focus, ef.k)

# Return the crossing of triangle `tri` with the plane on which the signed heights `h` of
# its vertices vanish. The function returns the two end points and their number, 0, 1 or
# 2. A vertex on the plane counts once, and an edge in the plane counts by its ends.
function _plane_cut(tri, h)
    p = q = tri[1]
    n = 0
    for (i, j) in ((1, 2), (2, 3), (3, 1))
        x = h[i] == 0 ? tri[i] : h[i] * h[j] < 0 ? tri[i] + h[i] / (h[i] - h[j]) * (tri[j] - tri[i]) :
            continue
        n == 0 ? (p = q = x) : (q = x)
        n += 1
    end
    return p, q, min(n, 2)
end

# Push the pieces of segment `p`–`q` upstream of `s_stop`, split at best focus, each with
# a lower bound of the envelope distance on it, the balancing direction `w` and `strip`.
function _segment_pieces!(segs, ef::_BeamConeFrame{Float64}, s_stop, p, q, w, strip::Bool)
    sp, sq = _axial(ef, p), _axial(ef, q)
    (sp > s_stop && sq > s_stop) && return segs
    if sp > s_stop || sq > s_stop
        x = p + (s_stop - sp) / (sq - sp) * (q - p)
        sp > s_stop ? ((p, sp) = (x, s_stop)) : ((q, sq) = (x, s_stop))
    end
    if (sp < 0 && sq > 0) || (sp > 0 && sq < 0)
        x = p + sp / (sp - sq) * (q - p)
        push!(segs, (_segment_bound(ef, p, x), p, x, w, strip), (_segment_bound(ef, x, q), x, q, w, strip))
    else
        push!(segs, (_segment_bound(ef, p, q), p, q, w, strip))
    end
    return segs
end

# Return a lower bound of the envelope distance at `x`. The bound is the larger of the
# distance to the circumscribed circle and, outside the ellipse, the gap (g − 1) min(A, B)
# to the homothetic ellipse through `x`. Both are at most the distance to the ellipse, and
# both are 1-Lipschitz across the beam.
function _point_bound(ef::_BeamConeFrame{Float64}, x)
    d = x - ef.focus
    u, v, s = dot(d, ef.e1), dot(d, ef.e2), dot(d, ef.k)
    lo = hypot(u, v) - max(ef.A, ef.B)
    if ef.A > 0 && ef.B > 0
        g = hypot(u / ef.A, v / ef.B)
        lo = max(lo, g >= 1 ? (g - 1) * min(ef.A, ef.B) : -min(ef.A, ef.B))
    end
    return lo - ef.cone * abs(s) - ef.margin
end

# Return a lower bound of the envelope distance on segment `p`–`q`, the bound at the
# midpoint less the Lipschitz constant √(1 + cone²) times the half length.
_segment_bound(ef::_BeamConeFrame{Float64}, p, q) =
    _point_bound(ef, (p + q) / 2) - sqrt(1 + ef.cone^2) * norm(q - p) / 2

# Return the gradient of the envelope distance where the unit normal of the ellipse is
# (n1, n2), on the side `σ` (±1) of best focus.
_envelope_gradient(ef::_BeamConeFrame{Float64}, n1, n2, σ) = n1 * ef.e1 + n2 * ef.e2 - (σ * ef.cone) * ef.k

# Return the envelope distance at `x` and its gradient on the side `σ` of best focus.
function _envelope_value_gradient(ef::_BeamConeFrame{Float64}, x, σ)
    d = x - ef.focus
    u, v = dot(d, ef.e1), dot(d, ef.e2)
    X0, X1, _ = _ellipse_closest(u, v, ef.A, ef.B)
    n1, n2 = _ellipse_unit_normal(u, v, ef.A, ef.B, X0, X1)
    return n1 * (u - X0) + n2 * (v - X1) - ef.cone * abs(dot(d, ef.k)) - ef.margin,
           _envelope_gradient(ef, n1, n2, σ)
end

# Return the smallest envelope distance on segment `p`–`q` and its point. The segment lies
# on one side of best focus, so the distance is convex along it. With `ends = false` the
# function returns the smallest value inside the segment, or `Inf` if that value is at an
# end. A bound at or above `best` ends the search early. With `gp` and `gq` the function
# minimizes the distance plus the linear function that takes those values at `p` and `q`.
function _segment_min(ef::_BeamConeFrame{Float64}, p, q, best, ends::Bool = true, gp = 0.0, gq = 0.0)
    Δ = q - p
    σ = _axial(ef, (p + q) / 2) < 0 ? -1 : 1
    fa, ga = _envelope_value_gradient(ef, p, σ)
    fb, gb = _envelope_value_gradient(ef, q, σ)
    fa, fb = fa + gp, fb + gq
    va, xa = !ends ? (Inf, p) : fa <= fb ? (fa, p) : (fb, q)
    ha, hb = dot(ga, Δ) + (gq - gp), dot(gb, Δ) + (gq - gp)
    (ha >= 0 || hb <= 0) && return va, xa                  # the minimum at an end
    # The tangents at the ends bound the convex function from below.
    fa + ha * (fb - hb - fa) / (ha - hb) >= min(best, va) && return va, xa
    # Use Illinois regula falsi on the derivative, and bisection when it stalls.
    a, b = 0.0, 1.0
    side = 0
    for _ in 1:100
        τ = (a * hb - b * ha) / (hb - ha)
        (a < τ < b) || (τ = (a + b) / 2)
        x = p + τ * Δ
        f, g = _envelope_value_gradient(ef, x, σ)
        f += gp + τ * (gq - gp)
        f < va && ((va, xa) = (f, x))
        h = dot(g, Δ) + (gq - gp)
        if h < 0
            a, ha = τ, h
            side == -1 && (hb /= 2)
            side = -1
        elseif h > 0
            b, hb = τ, h
            side == 1 && (ha /= 2)
            side = 1
        else
            break
        end
        b - a <= 1e-13 && break
    end
    return va, xa
end

# `_ellipse_closest` and `_ellipse_unit_normal` give the closest point of an ellipse and
# the normal there. They are defined in `clearance.jl`, beside `_ellipse_distance`, which
# shares them.

# Whether (u, v) lies on the medial axis of the ellipse with semi-axes (A, B), to within
# a distance of 1e-10 of the semi-major axis. The point must lie inside the segment of the
# major axis between the centers of curvature. The mirror image (Y0, Y1) of its closest
# point (X0, X1) across the major axis must also be as near as that point, to within the
# tolerance. A minimum the
# search finds on the medial strip lies within ~1e-13 mm of the axis (3e-14 measured on
# the test mock-up), its mirror point within 1e-13 mm as near. Next to a center of
# curvature the two points nearly coincide, and so do their distances and derivatives.
function _medial_tie(u, v, A, B, X0, X1, Y0, Y1)
    e0, e1, t = A >= B ? (A, B, u) : (B, A, v)
    abs(t) * e0 < (e0 - e1) * (e0 + e1) || return false
    return hypot(u - Y0, v - Y1) - hypot(u - X0, v - X1) <= 1e-10 * e0
end

# The signed distance, negative inside, from (u, v) to the ellipse with semi-axes (A, B),
# in their number type. It is the projection on the unit normal (n1, n2) at the closest
# point (X0, X1). That point is found on the values, the semi-axes `A0` and `B0`, and is
# held by its parameter. The value is the distance, and the first derivatives are those of
# the distance, because the closest point is stationary (Danskin's theorem).
_ellipse_signed(u, v, A, B, A0::Float64, B0::Float64, X0, X1, n1, n2) =
    n1 * (u - A * (A0 > 0 ? X0 / A0 : 0.0)) + n2 * (v - B * (B0 > 0 ? X1 / B0 : 0.0))

# ---------------------------------------------------------------------------------------
# Capture by the struck part
# ---------------------------------------------------------------------------------------

"""
    CaptureKeepIn(target, parts; peak_fraction = 0.1, intensity = nothing, pointing = 0.05,
                  target_offset = 0.025, hardware_tolerance = 0.0, focus = 1.0,
                  margin = 0.0, tolerance = 0.5, pivot = :tcc, minimum = 0.5)

Capture constraint for the parts meant to be struck by beams, such as backlighter foils
and capture plates. [`capture_margin`](@ref)`(keepin, beams)` gives the worst-case margin
(mm) by which the envelope of each beam stays inside the substrate edge of the struck
part. The margin is the TaLIS 3ω capture rule of [`talis_check`](@ref) (`:capture_3w`) as
a quantity, with `parts` in the role of the `capture` of `talis_check`. The struck part is
the first of `parts` (a name or a collection of part or component names) that the nominal
beam axis meets; its substrate
edge is its free edges (edges of a single opaque triangle, as `talis_check` takes them).
The margin is the smallest signed distance from that edge to the beam envelope (the
best-focus contour at `peak_fraction` of the peak, or at `intensity` W/cm², grown by the
f-cone and `margin`; negative inside it), less the transverse error radius `pointing +
target_offset + hardware_tolerance`, `focus/(2F)` for the focus error and the rule
tolerance `r·tan(tolerance)` (`r` the distance from `pivot`, as in
[`HardwareKeepOut`](@ref)), point by point. It is minimized over the edge as
`HardwareKeepOut` minimizes over the segments of its parts: exactly, continuous in the
aim, defocus and spot size of the beams, with the derivatives of the distance at the
closest point.
A part without a free edge (a closed surface) gives `Inf`, as `talis_check` passes it.

A beam whose axis meets none of `parts` fails the rule (`talis_check` reports `-Inf`). The
quantity continues the margin there, to keep the penalty of a stray beam finite and its
gradient pointing back: the margin from the outline of the part nearest the axis, less twice
the transverse distance `d` of the axis from that outline. Where the axis crosses the
outline `d` vanishes, so the margin is continuous there (on a straight edge it continues
linearly), and it falls at least as fast as the axis moves away. The margin is `-Inf` when
no part has a free edge. The margin jumps where the landing moves from one part to another.

The defaults are the TaLIS values (the 10% contour at least 0.5 mm inside the edge
with a 0.5° tolerance about TCC; `talis_check` passes a margin at or above 0.5 mm). The
quantity is a surrogate checked against `talis_check` in this package (see the
research notes), not against the facility tool. Pass it as `keepout` to the objectives,
which add `weight × Σᵢ max(0, minimum − mᵢ)²` ([`keepout_penalty`](@ref)).
"""
struct CaptureKeepIn <: AbstractKeepOut
    names::Vector{Symbol}
    target::Target
    comps::Set{Int}                              # the parts' components
    part_of::Dict{Int,Int}                       # part of each component (in `names`)
    clusters::Vector{Vector{_TriangleCluster}}   # each part's free edges
    peak_fraction::Union{Nothing,Float64}
    intensity::Union{Nothing,Float64}
    transverse::Float64
    focus::Float64
    margin::Float64
    tolerance::_TaLISTol
    minimum::Float64
end

function CaptureKeepIn(t::Target, parts; peak_fraction = nothing, intensity = nothing,
                       pointing = 0.05, target_offset = 0.025, hardware_tolerance = 0.0,
                       focus = 1.0, margin = 0.0, tolerance = 0.5, pivot = :tcc, minimum = 0.5)
    intensity === nothing && peak_fraction === nothing && (peak_fraction = 0.1)
    names = _part_names(parts)
    pieces = [part_pieces(t, n) for n in names]
    part_of = Dict(ci => j for (j, ps) in enumerate(pieces) for ci in ps)
    fe = _free_edges(t, collect(keys(part_of)))
    clusters = [_clusters(reduce(vcat, (fe[ci] for ci in ps); init = NTuple{2,Vec3}[])) for ps in pieces]
    return CaptureKeepIn(names, t, Set(keys(part_of)), part_of, clusters, peak_fraction, intensity,
                         pointing + target_offset + hardware_tolerance, focus, margin,
                         _talis_tol(tolerance, pivot), minimum)
end

"""
    capture_margin(keepin, beams) -> Vector

Worst-case capture margin (mm) of each beam on the part of a [`CaptureKeepIn`](@ref) struck
by the beam, the rule tolerance included; continued past the outline of the part for a beam
whose axis misses every part, `Inf` on a part without a free edge. The margin is
differentiable in the aim, defocus and spot size of the beams and has their number type.
"""
capture_margin(k::CaptureKeepIn, beams::AbstractVector{<:Beam}) = [_capture_margin(k, b) for b in beams]

keepout_clearance(k::CaptureKeepIn, beams::AbstractVector{<:Beam}) = capture_margin(k, beams)

const _NO_TOLERANCE = _TaLISTol(0.0, :tcc, zero(Vec3))

function _capture_margin(k::CaptureKeepIn, b::Beam)
    ef, ρ, aim = _envelope_frame(k, b)
    _, hit = _landing_hit(k.target, _axis_envelope(ef, aim), k.comps; miss = Inf)
    hit > 0 && return _envelope_clearance(k, b, ef, ρ, aim, k.clusters[k.part_of[hit]], Inf)
    # When the axis misses, find the part whose outline lies nearest the axis. The
    # distance is taken to the envelope of a point-like spot along the axis, without the
    # f-cone.
    ax = _BeamConeFrame{Float64}(aim, ef.k, ef.e1, ef.e2, 0.0, 0.0, 0.0, 0.0)
    j, d, xa = 0, Inf, aim
    for (i, cl) in enumerate(k.clusters)
        isempty(cl) && continue
        _, x, _ = _hardware_closest(cl, _NO_TOLERANCE, ax, Inf, aim)
        dx = hypot(dot(x - aim, ef.e1), dot(x - aim, ef.e2))
        dx < d && ((j, d, xa) = (i, dx, x))
    end
    j == 0 && return -_no_clearance(b)
    δ = xa - b.aim
    return _envelope_clearance(k, b, ef, ρ, aim, k.clusters[j], Inf) - 2 * hypot(dot(δ, ef.e1), dot(δ, ef.e2))
end

# ---------------------------------------------------------------------------------------
# Incidence on the struck part
# ---------------------------------------------------------------------------------------

"""
    IncidenceRule(target, parts; peak_fraction = 0.1, intensity = nothing, pointing = 0.05,
                  target_offset = 0.025, hardware_tolerance = 0.0, focus = 1.0, margin = 0.0,
                  max_incidence = 70.0, sharpness = 100.0, level = 2, minimum = 0.0)

Incidence limit on the parts struck by beams. [`incidence_angle`](@ref)`(rule, beams)`
gives the worst-case incidence (degrees) of each beam on the struck part, the first of
`parts` that the nominal beam axis meets. The quantity is the TaLIS 3ω incidence rule of
[`talis_check`](@ref) (`:incidence_3w`, at most `max_incidence`). As in `talis_check`, the
part is cut into `4^level` sub-triangles per mesh triangle, each with its normal, which is
exact on sphere, cylinder, cone and torus surfaces. The incidence is the largest angle
between the beam axis and the normal of a sub-triangle that is visible from upstream and
within reach of the worst-case footprint, plus the f-cone half-angle `atan(1/(2F))`, at
most 90°. The worst-case footprint is the envelope at `peak_fraction` of the peak, grown
by the f-cone and `margin` and enlarged by the transverse error radius
`pointing + target_offset + hardware_tolerance` and `focus/(2F)`. The rule tolerance does
not enter.

The set of sub-triangles within reach of the footprint changes as a beam moves, so the
largest angle over it steps. The quantity smooths the step. Each sub-triangle `f` enters
with its angle less `sharpness × max(0, d_f)` in degrees, where `d_f` (mm) is the
clearance of the sub-triangle from the footprint. `d_f` is the exact minimum over the
sub-triangle, found as in [`HardwareKeepOut`](@ref). A sub-triangle thus counts in full
once the footprint reaches it, and its weight fades over `Δθ/sharpness` mm outside it. The
incidence is continuous in the aim, defocus and spot size of the beams, and it is never
below the incidence of `talis_check`, to its resolution. It equals that incidence where no
sub-triangle outside the footprint lies within `Δθ/sharpness` of it, with `Δθ` the amount
by which the angle of that sub-triangle exceeds those of the reached sub-triangles. On a
flat part the incidence is the tilt of the part to the beam plus the half-angle,
independent of the pointing. `sharpness` (100°/mm) is an engineering default and is not set
by a facility criterion. The other defaults are those of `talis_check`. A beam whose axis
meets none of `parts` has no incidence (`-Inf`); [`CaptureKeepIn`](@ref) gives the capture
rule that such a beam fails.

As a keep-out ([`keepout_clearance`](@ref), [`keepout_penalty`](@ref)) the quantity of the
rule is the margin `max_incidence − incidence` in degrees, held at or above `minimum`. The
quantity is a surrogate checked against `talis_check` in this package (see the research
notes), not against the facility tool.
"""
struct IncidenceRule <: AbstractKeepOut
    names::Vector{Symbol}
    target::Target
    comps::Set{Int}
    part_of::Dict{Int,Int}
    leaves::Vector{Vector{NTuple{3,Vec3}}}       # each part's sub-triangles
    normals::Vector{Vector{Vec3}}                # and their normals
    peak_fraction::Union{Nothing,Float64}
    intensity::Union{Nothing,Float64}
    transverse::Float64
    focus::Float64
    margin::Float64
    tolerance::_TaLISTol                         # zero: incidence has no tolerance
    max_incidence::Float64
    sharpness::Float64
    minimum::Float64
end

function IncidenceRule(t::Target, parts; peak_fraction = nothing, intensity = nothing,
                       pointing = 0.05, target_offset = 0.025, hardware_tolerance = 0.0,
                       focus = 1.0, margin = 0.0, max_incidence = 70.0, sharpness = 100.0,
                       level::Integer = 2, minimum = 0.0)
    intensity === nothing && peak_fraction === nothing && (peak_fraction = 0.1)
    sharpness > 0 || throw(ArgumentError("sharpness must be > 0"))
    names = _part_names(parts)
    pieces = [part_pieces(t, n) for n in names]
    part_of = Dict(ci => j for (j, ps) in enumerate(pieces) for ci in ps)
    leaves = [NTuple{3,Vec3}[] for _ in pieces]
    normals = [Vec3[] for _ in pieces]
    for (j, ps) in enumerate(pieces), ci in ps
        t.components[ci].kind === :surface || continue
        g = _talis_samples(t, ci, level)
        append!(leaves[j], ((a, b, c) for (a, b, c, _) in g.leaves))
        append!(normals[j], g.normals)
    end
    return IncidenceRule(names, t, Set(keys(part_of)), part_of, leaves, normals, peak_fraction,
                         intensity, pointing + target_offset + hardware_tolerance, focus, margin,
                         _NO_TOLERANCE, max_incidence, sharpness, minimum)
end

"""
    incidence_angle(rule, beams) -> Vector

Worst-case incidence (degrees) of each beam on the part of an [`IncidenceRule`](@ref)
struck by the beam, the f-cone half-angle included; `-Inf` when no part is struck.
The incidence is continuous in the aim, defocus and spot size of the beams and has their
number type.
"""
incidence_angle(k::IncidenceRule, beams::AbstractVector{<:Beam}) = [_incidence_angle(k, b) for b in beams]

keepout_clearance(k::IncidenceRule, beams::AbstractVector{<:Beam}) =
    [k.max_incidence - θ for θ in incidence_angle(k, beams)]

function _incidence_angle(k::IncidenceRule, b::Beam{T}) where {T}
    ef, ρ, aim = _envelope_frame(k, b)
    _, hit = _landing_hit(k.target, _axis_envelope(ef, aim), k.comps; miss = Inf)
    hit > 0 || return T(-Inf)
    j = k.part_of[hit]
    leaves, normals = k.leaves[j], k.normals[j]
    t = k.target
    shift = k.transverse + k.focus * ef.cone
    L = sqrt(1 + ef.cone^2)
    far = 4 * (t.scale + norm(aim)) + 1
    ε = 1e3 * _ray_eps(t)
    # Search best first by an upper bound of the smoothed angle of each sub-triangle. The
    # bound is its angle less sharpness × a lower bound of its clearance from the footprint.
    θs = [acosd(clamp(abs(dot(n, ef.k)), 0.0, 1.0)) for n in normals]
    heap = Tuple{Float64,Int}[]
    for (i, (a, bb, c)) in enumerate(leaves)
        m = (a + bb + c) / 3
        r = max(norm(a - m), norm(bb - m), norm(c - m))
        lb = _point_bound(ef, m) - L * r - shift
        _heap_push!(heap, (-(θs[i] - k.sharpness * max(0.0, lb)), i))
    end
    best, ibest, cl = -Inf, 0, _TriangleCluster[]
    while !isempty(heap)
        nub, i = _heap_pop!(heap)
        -nub <= best && break
        a, bb, c = leaves[i]
        m = (a + bb + c) / 3
        visible(t, m - far * ef.k, m - ε * ef.k) || continue
        cl = _clusters([leaves[i]])
        found, x, _ = _hardware_closest(cl, k.tolerance, ef, Inf, aim)
        d = first(_envelope_value_gradient(ef, x, 1)) - shift
        v = θs[i] - k.sharpness * max(0.0, d)
        v > best && ((best, ibest) = (v, i))
    end
    ibest == 0 && return T(-Inf)
    hc = atand(1 / (2b.fnumber))
    a, bb, c = leaves[ibest]
    cl = _clusters([leaves[ibest]])
    # Evaluate the clearance in the number type of the beams where the smoothing acts.
    d = _envelope_clearance(k, b, ef, ρ, aim, cl, Inf)
    θ = θs[ibest] - k.sharpness * max(zero(d), d)
    return min(θ + hc, T(90))
end

# ---------------------------------------------------------------------------------------
# Keep-out from the 1ω footprints
# ---------------------------------------------------------------------------------------

"""
    Omega1KeepOut(target, parts; pointing = 0.05, target_offset = 0.025,
                  hardware_tolerance = 0.0, focus = 1.0, tolerance = 0.5, pivot = :tcc,
                  omega1_edge = 4.45, omega1_size = nothing, omega1_focus = nothing,
                  lens_distance = NIF_LENS_DISTANCE, offsets = nothing, resolution = 1e-9,
                  minimum = 0.2)

Keep-out of named parts from the 1ω footprints of the beams.
[`omega1_clearance`](@ref)`(keepout, beams)` gives the worst-case clearance (mm) of `parts`
from the 10% contour of the unconverted 1ω light of each beam, where `parts` is a name or
a collection of part or component names. The quantity is the TaLIS 1ω clearance rule of
[`talis_check`](@ref) (`:clearance_1w`, more than 0.2 mm). The footprint and its keywords
are those of `talis_check`, described there under "1ω footprint". Per beamlet of a NIF
quad, the footprint is a square whose near edge lies `omega1_edge` from the 3ω axis. It
tapers to the 1ω focus `omega1_focus` beyond the 3ω focus and moves by `±focus` with the
3ω focal plane. The clearance is the smallest transverse distance from the parts to the
footprint, less the transverse error radius `pointing + target_offset + hardware_tolerance`
and the rule tolerance `r·tan(tolerance)` about `pivot` at each point. Beams without a NIF
quad name or an `offsets` entry have no 1ω footprint and get `Inf`.

The minimum over the parts is found by `talis_check`'s branch and bound on the Lipschitz
constant of the footprint distance, to `resolution` (mm), not by an exact search: the
clearance is the smallest value found, within `resolution` of the minimum and
continuous in the aim and defocus of the beams to that accuracy. Its derivatives are
those of the footprint distance at that point (Danskin's theorem; one-sided at the
corners of the footprint and where the footprints of two beamlets are equally near). Not
clearing calls for unconverted-light management (`:management_1w`), a review outside the
scope of this quantity. The defaults are `talis_check`'s; the quantity is a surrogate
checked against `talis_check` in this package (see the research notes), not against the
facility tool.
"""
struct Omega1KeepOut{O,B} <: AbstractKeepOut
    names::Vector{Symbol}
    target::Target
    tris::Vector{NTuple{3,Vec3}}
    transverse::Float64
    focus::Float64
    tolerance::_TaLISTol
    c0::Float64                      # the footprint center's distance from the 3ω axis
    h0::Float64                      # its half side
    f1::Float64                      # the 1ω focus beyond the 3ω focus
    L::Float64                       # the lens distance (chief ray)
    lip::Float64                     # Lipschitz constant of the footprint distance
    offsets::O                       # footprint directions by beam name, or `nothing`
    table::B                         # the NIF beam table, for the beamlets
    resolution::Float64
    minimum::Float64
end

function Omega1KeepOut(t::Target, parts; pointing = 0.05, target_offset = 0.025,
                       hardware_tolerance = 0.0, focus = 1.0, tolerance = 0.5, pivot = :tcc,
                       omega1_edge = 4.45, omega1_size = nothing, omega1_focus = nothing,
                       lens_distance = NIF_LENS_DISTANCE, offsets = nothing, resolution = 1e-9, minimum = 0.2)
    resolution > 0 || throw(ArgumentError("resolution must be > 0"))
    names = _part_names(parts)
    idx = Set(ci for n in names for ci in part_pieces(t, n))
    tris = [tri for (k, tri) in zip(findall(in(idx), t.tri_comp), _part_triangles(t, idx)) if t.tri_area[k] > 0]
    isempty(tris) && throw(ArgumentError("parts $(join(repr.(names), ", ")) have no triangles"))
    (; c0, h0, f1) = _omega1_params(omega1_edge, omega1_size, omega1_focus)
    L = Float64(lens_distance)
    return Omega1KeepOut(names, t, tris, pointing + target_offset + hardware_tolerance, focus,
                         _talis_tol(tolerance, pivot), c0, h0, f1, L, _omega1_lipschitz(c0, h0, f1, L),
                         offsets, _nif_beam_table(), resolution, minimum)
end

"""
    omega1_clearance(keepout, beams) -> Vector

Worst-case clearance (mm) of the parts of an [`Omega1KeepOut`](@ref) from the 1ω
footprint of each beam, the rule tolerance included; `Inf` for a beam without a 1ω
footprint. The clearance is differentiable in the aim and defocus of the beams and has
their number type.
"""
omega1_clearance(k::Omega1KeepOut, beams::AbstractVector{<:Beam}) = [_omega1_clearance(k, b) for b in beams]

keepout_clearance(k::Omega1KeepOut, beams::AbstractVector{<:Beam}) = omega1_clearance(k, beams)

function _omega1_clearance(k::Omega1KeepOut, b::Beam{T}) where {T}
    bv = _value_beam(b)
    bl = _omega1_beamlets(bv, k.offsets, k.table)
    bl === nothing && return T(Inf)
    aim = Vec3(bv.aim)
    cx = (; focus = k.focus, c0 = k.c0, h0 = k.h0, f1 = k.f1, L = k.L)
    τ = k.tolerance
    f = x -> _omega1_distance(x, bl, cx)[2] - k.transverse - τ(x, aim)
    lb, x = _lipschitz_minimum(f, k.tris, k.lip + τ.tanα, (k.resolution, 0.0, k.minimum), aim)
    isfinite(lb) || return T(Inf)
    # Evaluate the footprint distance at that point in the number type of the beams. The
    # 3ω focus of each beamlet moves with the aim and defocus.
    d = T(Inf)
    for (kb, _, ô) in bl
        fb = b.aim + b.defocus * kb
        d = min(d, _omega1_square(x - fb, kb, ô, k.c0, k.h0, k.f1, k.L))
        if k.focus > 0
            d = min(d, _omega1_square(x - fb + k.focus * kb, kb, ô, k.c0, k.h0, k.f1, k.L),
                    _omega1_square(x - fb - k.focus * kb, kb, ô, k.c0, k.h0, k.f1, k.L))
        end
    end
    return d - k.transverse - _tol_value(τ, x, b.aim)
end
