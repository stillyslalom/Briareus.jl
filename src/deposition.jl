# This file computes laser deposition with a differentiable quadrature operator and with
# an energy-conserving ray trace.

# Run `f(i)` for `i in 1:n` on all threads, handing out indices one at a time in order.
# This balances the load dynamically; with heterogeneous cores, a static split would wait
# on the slowest core.
# Each `f(i)` must write only outputs owned by `i` for results to be deterministic.
function _tforeach(f, n::Integer)
    nt = min(Threads.nthreads(), n)
    if nt <= 1
        for i in 1:n
            f(i)
        end
        return nothing
    end
    next = Threads.Atomic{Int}(1)
    @sync for _ in 1:nt
        Threads.@spawn while (i = Threads.atomic_add!(next, 1)) <= n
            f(i)
        end
    end
    return nothing
end

# The points of each beam are sorted by destination, with elements ascending and then
# back sides by component. `deposit` splits the work into units that own their outputs:
# - A block of consecutive elements takes the points of every beam in it and is the only
#   writer of those elements.
# - The back-side points of each beam are cut into pieces with separate partial sums.
# Units are handed to threads dynamically, largest first. Every output is summed in a
# fixed order, as a serial loop would sum it: each element over beams, then over points.
# The result is therefore identical for any thread count and scheduling.
#
# `tiles[bi][t]` bounds the points `_TILE (t - 1) + 1 : _TILE t` of beam `bi`, in their
# stored order by destination, by a ball with a center and radius. A screen can then skip
# a run of points that all lie beyond the cut-off of a beam (`_vr_far_ball`).
struct _DepositPlan
    boff::Matrix{Int}                     # (nblock+1) × nb: first point of each element block; last row: first back-side point
    back::Vector{Tuple{Int,UnitRange{Int}}}  # back-side pieces (beam, points)
    order::Vector{Int}                    # units by decreasing size; ≤ nblock: element block, else back-side piece
    tiles::Vector{Vector{NTuple{4,Float64}}}  # per beam, per run of `_TILE` points: ball center and radius
end

const _TILE = 8

# Returns the balls of `_DepositPlan.tiles` for the points of one beam. Each ball has the
# center of the bounding box and the largest distance from it, grown by 1e-12 relative.
# The ball is a bound and is not the least enclosing ball.
function _point_tiles(P::Vector{Vec3})
    tiles = Vector{NTuple{4,Float64}}(undef, cld(length(P), _TILE))
    for t in eachindex(tiles)
        r = _TILE*(t-1)+1:min(_TILE * t, length(P))
        lo, hi = P[first(r)], P[first(r)]
        for q in r
            lo, hi = min.(lo, P[q]), max.(hi, P[q])
        end
        c = (lo + hi) / 2
        R = maximum(q -> norm(P[q] - c), r)
        tiles[t] = (c[1], c[2], c[3], R * (1 + 1e-12) + 1e-12)
    end
    return tiles
end

"""
    DepositionOperator(target, beams; level = 1, ncone = 13, resolution = 0.3,
                       reach = 2.0, maxlevel = 5, shadow = 0,
                       backside_resolution = 1.0, edge_resolution = 0.01,
                       profile = nothing)

Precomputed laser-deposition quadrature. For each beam, every opaque triangle side that
faces the beam and belongs to an element is sampled at `4^level` points. More points are
used where the beam can land:

- A triangle within `reach` (mm) of the beam footprint is divided until its sub-triangles
  are no longer than `resolution` times the smaller semi-axis of the spot, up to
  `4^maxlevel` points. The footprint is 2.5 spot semi-axes, grown by the f-cone along the
  beam axis through the aim.
- With `shadow > 0`, such a triangle is divided `shadow` times more where a shadow edge
  crosses it, so that its points see different cone directions.

`reach` is the aim motion over which the operator stays accurate. Beyond it, deposition is
resolved only at `level`. The defaults conserve the power of each beam to 1e-4 when nothing
clips and to 2e-3 when an entrance-hole lip clips several percent of it, on fine and coarse
hohlraum meshes alike.

Unmodeled back sides, such as the outside of a hohlraum, carry the clipped `backside`
loss. Only their total per component is needed, and the open edges that the light passes
determine it. When an entrance-hole lip clips 1e-4 of a beam, the clipped light falls
within ~1/50 of a spot semi-axis of the rim. Back sides away from the beam are sampled at
one point per triangle. Within `reach` they are sampled as follows:

- A triangle on an open edge, or touching one at a vertex, gets a Gauss rule graded
  geometrically toward the edge, from `resolution` down to `edge_resolution` spot
  semi-axes. An open edge is an edge that no other opaque triangle shares, such as an
  entrance-hole rim.
- A triangle within one larger spot semi-axis of an open edge, where the shoulder of a
  clipped spot is deposited, is divided down to `resolution`.
- The rest are divided down to `backside_resolution`.

On the NIF hohlraum, the clipped power then matches ray traces, averaged over nspot
241–361, to 1% from 1e-4 to 5% of a beam, 2.5% at 1e-6, and 0.6% up to 50%. It is also
smooth in the aim. Element values do not depend on these two keywords.

Each point stores its area weight, the normal of the lit side and a bitmask of which of
the `ncone` f-cone directions reach it from the lens without hitting another surface.
These depend only on the target and the beam *directions*, so the operator stays valid
as pointing, defocus, spot size and power change. Deposition is then a smooth function of
those parameters (see `deposit`).

The beam is the superposition of `ncone` collimated sub-beams, one per f-cone
direction `dⱼ` (`cone_directions`), each carrying the best-focus profile. A point
receives from each visible direction the focal profile where the ray along `dⱼ`
through it crosses the best-focus plane, times the incidence factor of that ray
`|n·dⱼ|/(k·dⱼ)`. This is the geometric-optics convolution of the focal spot with the
f-cone, sampled at the cone directions as the reference ray trace samples it. The power
on any surface that intercepts the whole beam is therefore the beam power, at any
distance from focus and any incidence. Shadows cast by objects near focus, such as
capsules at TCC and LEH edges, are correlated with the part of the spot that they block.

`profile = VisRadProfile(...)` replaces this with VisRad's beam model, so that a design
can be checked in VisRad and its deposition reproduced. The model is a super-Gaussian
whose semi-axes and exponent grow away from focus and which blends toward a square. It is
cut off at 10⁻³ of its local peak, or tapered to the cut for optimization with
`VisRadProfile(; taper)`. Points are then sampled only along the beam axis, `ncone` is
ignored, and the beam axis through each point determines its shadowing. The two models
agree near focus. Away from focus the f-cone convolution softens the edges, while
VisRad's spot turns square and grows faster; on NIF hohlraum walls the peaks of the
f-cone model are 10% higher than VisRad's. Power on a surface that intercepts the whole
beam equals the beam power in planes normal to the beam and only approximately on tilted
ones. This follows from the per-plane normalization, which is 0.2% high at 40° incidence
and 2.2% at 78°, where VisRad's totals are exact.

The beam ports must not change after construction.
"""
struct DepositionOperator
    target::Target
    ports::Vector{Vec3}
    dirs::Vector{Vector{Vec3}}     # cone directions per beam
    pts::Vector{Vector{Vec3}}
    w::Vector{Vector{Float64}}     # area weight (mm²)
    normal::Vector{Vector{Vec3}}   # normal of the lit side
    mask::Vector{Vector{UInt16}}   # visible cone directions (bit j-1 set: direction j)
    dest::Vector{Vector{Int32}}    # > 0: element; < 0: -(component) for an unmodeled back side
    plan::_DepositPlan             # work units of the threaded `deposit`
    profile::Union{Nothing,VisRadProfile}   # nothing: the f-cone model
    envelope::Vector{Float64}      # per beam: VisRad envelope ratio (empty for the f-cone model)
end

const _DEPOSIT_UNIT = 4096   # points per work unit

function _deposit_plan(ne::Integer, dest::Vector{Vector{Int32}}, pts::Vector{Vector{Vec3}};
                       unit::Integer = _DEPOSIT_UNIT)
    nb = length(dest)
    cnt = zeros(Int, ne)
    for D in dest, d in D
        d > 0 && (cnt[d] += 1)
    end
    bounds = [1]
    acc = 0
    for e in 1:ne
        acc += cnt[e]
        if acc >= unit && e < ne
            push!(bounds, e + 1)
            acc = 0
        end
    end
    push!(bounds, ne + 1)
    nblock = length(bounds) - 1
    boff = Matrix{Int}(undef, nblock + 1, nb)
    back = Tuple{Int,UnitRange{Int}}[]
    for (bi, D) in enumerate(dest)
        q = 1
        for j in 1:nblock+1
            while q <= length(D) && 0 < D[q] < bounds[j]
                q += 1
            end
            boff[j, bi] = q
        end
        issorted(@view D[1:q-1]) && all(<(0), @view D[q:end]) ||
            throw(ArgumentError("operator points are not sorted by destination"))
        nq = length(D) - q + 1
        npiece = cld(nq, unit)
        for p in 1:npiece
            push!(back, (bi, q + (p - 1) * nq ÷ npiece:q + p * nq ÷ npiece - 1))
        end
    end
    cost = vcat([sum(boff[j+1, bi] - boff[j, bi] for bi in 1:nb; init = 0) for j in 1:nblock],
                [length(r) for (_, r) in back])
    return _DepositPlan(boff, back, sortperm(cost; rev = true), [_point_tiles(P) for P in pts])
end

DepositionOperator(t::Target, ports, dirs, pts, w, normal, mask, dest,
                   profile::Union{Nothing,VisRadProfile} = nothing, envelope = Float64[]) =
    DepositionOperator(t, ports, dirs, pts, w, normal, mask, dest, _deposit_plan(nelements(t), dest, pts),
                       profile, envelope)
DepositionOperator(t::Target, ports, dirs, pts, w, normal, mask, dest, plan::_DepositPlan) =
    DepositionOperator(t, ports, dirs, pts, w, normal, mask, dest, plan, nothing, Float64[])

function DepositionOperator(t::Target, beams::AbstractVector{<:Beam}; level::Integer = 1,
                            ncone::Integer = 13, resolution = 0.3, reach = 2.0,
                            maxlevel::Integer = 5, shadow::Integer = 0,
                            backside_resolution = 1.0, edge_resolution = 0.01,
                            profile::Union{Nothing,VisRadProfile} = nothing)
    1 <= ncone <= 16 || throw(ArgumentError("ncone must be in 1:16"))
    profile === nothing || (ncone = 1)               # VisRad's model: the beam axis only
    maxlv = max(maxlevel, level)
    nb = length(beams)
    # Build the sub-triangle quadratures in barycentric form, shared by all triangles.
    rules = [_bary_quadrature(l) for l in 0:maxlv]
    ε = _ray_eps(t)
    geo = [_operator_beam(b, ncone, profile) for b in beams]
    flags, rimdist = _open_boundary(t)
    back = (; resolution, coarse = backside_resolution, edge = edge_resolution, flags, rimdist)
    # Work units are pairs of a beam and a triangle range. Beams differ in cost, and the
    # triangles near the footprint of a beam are refined, so units are small and handed
    # out dynamically. The points of each unit are stored apart and joined in order, so
    # the result does not depend on scheduling.
    tris = findall(t.tri_opaque)
    nchunk = clamp(cld(length(tris), 256), 1, 16)
    ranges = [tris[(i-1)*length(tris)÷nchunk+1:i*length(tris)÷nchunk] for i in 1:nchunk]
    parts = Matrix{_OperatorPoints}(undef, nchunk, nb)
    _tforeach(nchunk * nb) do u
        ic, bi = mod1(u, nchunk), cld(u, nchunk)
        parts[ic, bi] = _operator_points(t, geo[bi], ranges[ic], rules, ε, level, maxlv,
                                         resolution, reach, shadow, back)
    end
    ne = nelements(t)
    key(d) = d > 0 ? Int(d) : ne - Int(d)
    pts = Vector{Vector{Vec3}}(undef, nb)
    w = Vector{Vector{Float64}}(undef, nb)
    normal = Vector{Vector{Vec3}}(undef, nb)
    mask = Vector{Vector{UInt16}}(undef, nb)
    dest = Vector{Vector{Int32}}(undef, nb)
    _tforeach(nb) do bi
        col = @view parts[:, bi]
        D = reduce(vcat, [p.D for p in col])
        perm = sortperm(key.(D))                      # stable: keeps point order per element
        pts[bi] = reduce(vcat, [p.P for p in col])[perm]
        w[bi] = reduce(vcat, [p.W for p in col])[perm]
        normal[bi] = reduce(vcat, [p.N for p in col])[perm]
        mask[bi] = reduce(vcat, [p.M for p in col])[perm]
        dest[bi] = D[perm]
    end
    envelope = profile === nothing ? Float64[] : [_vr_envelope(profile, b) for b in beams]
    return DepositionOperator(t, [b.port for b in beams], [g.ds for g in geo], pts, w, normal,
                              mask, dest, profile, envelope)
end

const _OperatorPoints = @NamedTuple{P::Vector{Vec3}, W::Vector{Float64}, N::Vector{Vec3},
                                    M::Vector{UInt16}, D::Vector{Int32}}

# Returns the beam quantities, as Float64 values, that the sampling of the operator
# depends on.
function _operator_beam(b::Beam, ncone, profile)
    k = direction(b)
    ds, _ = cone_directions(b, ncone)
    sp = b.spot
    amin = ForwardDiff.value(min(sp.a, sp.b))
    # The footprint is 2.5 semi-axes grown by the f-cone. VisRad's profile is nonzero only
    # where α < ln 1000. That region is inside its ellipse q < ρ_d², with
    # ρ_d = (ln 1000)^(1/n₀) since n ≥ n₀, or inside its square
    # |u|, |v| < ΔL (ln 1000)^(1/M) < 1.22 ΔL, since M ≥ 10. The profile is therefore
    # nonzero only within ρ_d a(z) ≤ ρ_d (a₀ + |z|/(√2 F)) or within
    # 1.22 √2 ΔL = 1.22 (ρ a₀ + √2 · 0.99605 |z|/(2F)) of the axis.
    amax = ForwardDiff.value(max(sp.a, sp.b))
    if profile === nothing
        r, cone = 2.5, 1 / (2b.fnumber)
    else
        ρd = visrad_envelope(sp.order)
        r = max(ρd, 1.22 * _vr_envelope(profile, b))
        cone = sqrt(2) * max(ρd, 1.22 * _VR_DL_SLOPE) / (2b.fnumber)
    end
    extent = r * amax * (sp.halo > 0 ? max(1, sp.halo_scale) : 1)
    aim = Vec3(ForwardDiff.value.(b.aim))
    focus = aim + ForwardDiff.value(b.defocus) * k
    return (; k, ds, amin, extent, focus, cone)
end

# Returns the lit quadrature points of one beam on the triangles `tris`.
function _operator_points(t::Target, g, tris, rules, ε, level, maxlevel, resolution, reach, shadow,
                          back)
    (; k, ds, amin, extent, focus, cone) = g
    P, W, N, M, D = Vec3[], Float64[], Vec3[], UInt16[], Int32[]
    # Cache the rays toward the lens and the last leaf that blocked each one.
    cache = (; dirs = -ds, invs = [_safe_inv(-d) for d in ds], leaf = zeros(Int32, length(ds)))
    for tri in tris
        n = t.tri_normal[tri]
        c = dot(k, n)
        c == 0 && continue
        side_n = c < 0 ? n : -n            # normal of the lit side
        e = c < 0 ? t.tri_front[tri] : t.tri_back[tri]
        target_id = e > 0 ? e : -t.tri_comp[tri]
        A = t.tri_area[tri]
        v0, e1, e2 = t.tri_v0[tri], t.tri_e1[tri], t.tri_e2[tri]
        lv = level
        L = max(norm(e1), norm(e2), norm(e2 - e1))
        cen = v0 + (e1 + e2) / 3
        rad = max(norm(v0 - cen), norm(v0 + e1 - cen), norm(v0 + e2 - cen))
        d = cen - focus
        s = dot(d, k)
        near = norm(d - s * k) - rad <= extent + cone * abs(s) + reach
        if e <= 0
            # On an unmodeled back side only the total matters (see `_backside_points!`).
            lv = 0
            if near
                _backside_points!(P, W, N, M, D, t, tri, back, amin, side_n, target_id, ds, ε, cache) &&
                    continue
                # The shoulder of a clipped spot is deposited within a spot semi-axis of an
                # open edge, so those triangles get the finer resolution.
                res = back.rimdist[tri] <= extent / 2.5 ? back.resolution : back.coarse
                lv = clamp(ceil(Int, log2(L / (res * amin))), 0, maxlevel)
            end
        elseif near && L > resolution * amin * 2^level
            lv = clamp(ceil(Int, log2(L / (resolution * amin))), level, maxlevel)
        end
        xs, ms = _visible_directions(t, rules[lv+1][1], v0, e1, e2, side_n, ds, ε, cache)
        # Resolve the triangle more finely when a shadow edge crosses it.
        if near && shadow > 0 && lv < maxlevel && any(!=(ms[1]), ms)
            lv = min(lv + shadow, maxlevel)
            xs, ms = _visible_directions(t, rules[lv+1][1], v0, e1, e2, side_n, ds, ε, cache)
        end
        for (x, m, wq) in zip(xs, ms, rules[lv+1][2])
            m == 0 && continue
            push!(P, x); push!(W, A * wq); push!(N, side_n); push!(M, m); push!(D, target_id)
        end
    end
    return (; P, W, N, M, D)
end

# Returns the quadrature points of a triangle and, for each point, the bitmask of cone
# directions that reach it unobstructed. These are the directions that face the lit side
# and whose ray toward the lens hits nothing. `any_hit` is called for all points and
# directions as one packet. `cache` holds the reversed directions and keeps the last
# blocking leaf of each direction across a work unit.
function _visible_directions(t::Target, bary, v0, e1, e2, side_n, ds, ε, cache)
    n = length(bary)
    xs = Vector{Vec3}(undef, n)
    os = Vector{Vec3}(undef, n)
    for (i, λ) in enumerate(bary)
        x = v0 + λ[1] * e1 + λ[2] * e2
        xs[i] = x
        os[i] = x + ε * side_n
    end
    ms = fill(_facing_mask(ds, side_n), n)
    _unblocked!(ms, t.bvh, os, cache.dirs, cache.invs, cache.leaf)
    return xs, ms
end

# Returns the cone directions that strike the side with normal `n`, where d·n < 0. The
# product is summed in a fixed order, so the result for a direction in the plane does not
# depend on how `dot` is compiled.
function _facing_mask(ds, n)
    m = UInt16(0)
    for (j, d) in enumerate(ds)
        (d[1] * n[1] + d[2] * n[2]) + d[3] * n[3] < 0 && (m |= UInt16(1) << (j - 1))
    end
    return m
end

# Finds the open boundary of the opaque surfaces. It consists of the triangle edges that
# no other opaque triangle shares, such as entrance-hole rims and foil edges, and the
# vertices on them. Per triangle, bit j-1 flags edge j, from vertex j to vertex j+1, and
# bit j+2 flags vertex j; the vertices are v0, v0+e1 and v0+e2.
function _open_boundary(t::Target)
    tol = 1e-9 * max(t.scale, 1.0)
    key(v) = (round(Int, v[1] / tol), round(Int, v[2] / tol), round(Int, v[3] / tol))
    corners(tri) = (key(t.tri_v0[tri]), key(t.tri_v0[tri] + t.tri_e1[tri]),
                    key(t.tri_v0[tri] + t.tri_e2[tri]))
    ekey(a, b) = a < b ? (a, b) : (b, a)
    opaque = findall(t.tri_opaque)
    count = Dict{NTuple{2,NTuple{3,Int}},Int}()
    for tri in opaque
        c = corners(tri)
        for j in 1:3
            kk = ekey(c[j], c[mod1(j + 1, 3)])
            count[kk] = get(count, kk, 0) + 1
        end
    end
    onedge = Set{NTuple{3,Int}}()
    for ((a, b), n) in count
        n == 1 && (push!(onedge, a); push!(onedge, b))
    end
    flags = zeros(UInt8, length(t.tri_area))
    for tri in opaque
        c = corners(tri)
        f = 0x00
        for j in 1:3
            count[ekey(c[j], c[mod1(j + 1, 3)])] == 1 && (f |= 0x01 << (j - 1))
            c[j] in onedge && (f |= 0x01 << (j + 2))
        end
        flags[tri] = f
    end
    # Bound the distance of each triangle from the open boundary from below.
    segs = [(a, b) for ((a, b), n) in count if n == 1]
    pos(kk) = Vec3(kk[1] * tol, kk[2] * tol, kk[3] * tol)
    S = [(pos(a), pos(b)) for (a, b) in segs]
    rimdist = fill(Inf, length(t.tri_area))
    for tri in opaque
        t.tri_front[tri] > 0 && t.tri_back[tri] > 0 && continue   # no unmodeled side
        v0, e1, e2 = t.tri_v0[tri], t.tri_e1[tri], t.tri_e2[tri]
        cen = v0 + (e1 + e2) / 3
        rad = max(norm(v0 - cen), norm(v0 + e1 - cen), norm(v0 + e2 - cen))
        dmin = Inf
        for (a, b) in S
            ab = b - a
            λ = clamp(dot(cen - a, ab) / dot(ab, ab), 0.0, 1.0)
            dmin = min(dmin, norm(cen - a - λ * ab))
        end
        rimdist[tri] = max(dmin - rad, 0.0)
    end
    return flags, rimdist
end

const _GL2 = (0.5 - sqrt(3) / 6, 0.5 + sqrt(3) / 6)

# Returns the panel ends on [0, 1] for `n` equal panels. With `lo` the first panel, and
# with `hi` the last, is divided geometrically with ratio `q` toward its end, down to
# width `wmin`.
function _graded_cuts(n, lo, hi, wmin, q)
    n = max(n, lo && hi ? 2 : 1)
    cuts = Float64[]
    for i in 0:n
        x = i / n
        if lo && i == 1
            m = max(0, ceil(Int, log(x / wmin) / log(1 / q)))
            append!(cuts, x * q^j for j in m:-1:1)
        end
        push!(cuts, x)
        if hi && i == n - 1
            m = max(0, ceil(Int, log((1 - x) / wmin) / log(1 / q)))
            append!(cuts, 1 - (1 - x) * q^j for j in 1:m)
        end
    end
    return cuts
end

# Adds the back-side points of a triangle on the open boundary. At an open edge, such as
# an entrance-hole lip that clips the tail of a beam, the light on a back side lies within
# a few e-folding lengths of the tail of the edge. At 1e-4 clipped that is ~1/50 of a spot
# semi-axis. This is far finer than spot-resolving sub-triangles, whose centroids then
# miss most of the light; they underestimated 1e-4 clipping by 10–60%.
#
# The triangle is mapped from the unit square, collapsed onto one vertex C, with u from C
# to the opposite edge and v across. The coordinates that reach open edges and vertices
# are graded geometrically toward them with ratio 1/2, down to `edge` spot semi-axes, and
# each panel has two Gauss points. A vertex on the rim needs this grading as much as an
# edge, because the triangles that fan out from it hold a share ~λ/L of the light at the
# rim, with λ the tail length and L the edge. Panels are `resolution` spot semi-axes
# across the triangle and `coarse` along an open edge. Returns false, with no points
# added, for a triangle off the open boundary.
function _backside_points!(P, W, N, M, D, t::Target, tri, back, amin, side_n, target_id, ds, ε,
                           cache)
    f = back.flags[tri]
    f == 0 && return false
    E(j) = (f >> (mod1(j, 3) - 1)) & 0x1 == 0x1      # edge j: vertex j to j+1
    V(j) = (f >> (mod1(j, 3) + 2)) & 0x1 == 0x1
    nE = E(1) + E(2) + E(3)
    C = nE >= 2 ? findfirst(j -> E(j) && E(j - 1), 1:3) :
        nE == 1 ? mod1(findfirst(E, 1:3) + 2, 3) : findfirst(V, 1:3)
    Pv, Qv = mod1(C + 1, 3), mod1(C + 2, 3)
    ulo = V(C)
    uhi = E(Pv) || (V(Pv) && !E(C)) || (V(Qv) && !E(Qv))
    vlo = E(C) || (V(Pv) && !E(Pv))
    vhi = E(Qv) || (V(Qv) && !E(Pv))
    v0, e1, e2 = t.tri_v0[tri], t.tri_e1[tri], t.tri_e2[tri]
    X = (v0, v0 + e1, v0 + e2)
    B = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    Lu = max(norm(X[Pv] - X[C]), norm(X[Qv] - X[C]))
    Lv = norm(X[Qv] - X[Pv])
    h, halong, wmin = back.resolution * amin, back.coarse * amin, back.edge * amin
    # Use two-point Gauss panels in u, and in v as many panels as the width of the u panel
    # needs.
    ucuts = _graded_cuts(ceil(Int, Lu / (vlo || vhi ? halong : h)), ulo, uhi, wmin / Lu, 0.5)
    hv = uhi ? halong : h
    bary = Tuple{Float64,Float64}[]
    wq = Float64[]
    for k in 1:length(ucuts)-1
        a, b = ucuts[k], ucuts[k+1]
        vcuts = _graded_cuts(ceil(Int, b * Lv / hv), vlo, vhi, wmin / (b * Lv), 0.5)
        for gu in _GL2, l in 1:length(vcuts)-1, gv in _GL2
            u = a + gu * (b - a)
            v = vcuts[l] + gv * (vcuts[l+1] - vcuts[l])
            push!(bary, ntuple(i -> B[C][i] + u * ((1 - v) * (B[Pv][i] - B[C][i]) +
                                                   v * (B[Qv][i] - B[C][i])), 2))
            push!(wq, 2u * (b - a) / 2 * (vcuts[l+1] - vcuts[l]) / 2)
        end
    end
    xs, ms = _visible_directions(t, bary, v0, e1, e2, side_n, ds, ε, cache)
    A = t.tri_area[tri]
    for (x, m, w) in zip(xs, ms, wq)
        m == 0 && continue
        push!(P, x); push!(W, A * w); push!(N, side_n); push!(M, m); push!(D, target_id)
    end
    return true
end

function _bary_quadrature(level)
    pts = Vec3[]
    ws = Float64[]
    _subtri!(pts, ws, Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), level)
    return [(p[1], p[2]) for p in pts], ws ./ sum(ws)
end

_fullmask(ds) = UInt16((1 << length(ds)) - 1)

"Operator restricted to the beams `idx` (same order as `idx`)."
Base.getindex(op::DepositionOperator, idx::AbstractVector{<:Integer}) =
    DepositionOperator(op.target, op.ports[idx], op.dirs[idx], op.pts[idx], op.w[idx],
                       op.normal[idx], op.mask[idx], op.dest[idx], op.profile,
                       op.profile === nothing ? Float64[] : op.envelope[idx])

function Base.show(io::IO, op::DepositionOperator)
    np = sum(length, op.w; init = 0)
    npart = sum(count(!=(_fullmask(ds)), M) for (M, ds) in zip(op.mask, op.dirs); init = 0)
    print(io, "DepositionOperator(", length(op.ports), " beams, ", np,
          " lit quadrature points, ", npart, " partially shadowed",
          op.profile === nothing ? "" : ", VisRad profile", ")")
end

"""
    Deposition

The absorbed laser power, with fields:

- `element` (TW), the power absorbed by each element;
- `backside` (TW), the power on the unmodeled back side of each component, such as the
  outside of a hohlraum, which receives the LEH clipping;
- `per_beam` (TW), the total of each beam deposited on any opaque surface.
"""
struct Deposition{T}
    element::Vector{T}
    backside::Vector{T}
    per_beam::Vector{T}
    delivered::T
end

"""
    deposit(op, beams; absorption = nothing) -> Deposition

Evaluates laser deposition for `beams`, with the same ports and order as at
construction. The result is differentiable in the `aim`, `power`, `defocus` and spot axes
of every beam. Power that does not strike an opaque surface, for instance power that
passes through a hohlraum and out the other LEH, is `delivered - sum(per_beam)`.
Threads share the work dynamically, and every sum is taken in a fixed order: the result
is identical on any number of threads.

By default all the power of every ray is deposited where the ray strikes. With
`absorption = η⊥` (0 < η⊥ ≤ 1), only the fraction

    η_a(ψ) = 1 - (1 - η⊥)^(cos³ψ)

of each ray is kept, with ψ its angle of incidence. That angle is the angle of each
f-cone direction, or of the beam axis with a VisRad profile.

This is inverse-bremsstrahlung absorption in a thin corona whose density falls
exponentially with radius. The absorption at normal incidence is
η⊥ = 1 - exp(-8``ν_{ei}`` L/3c), with ``ν_{ei}`` the collision frequency at critical
density and L the scale length (Murakami & Nishi, Matter Radiat. Extremes 2, 55 (2017),
Eq. 29, after Kruer's text). It is a model for hard-sphere direct-drive illumination,
derived for that profile, and it is neither measured nor calibrated. η⊥ therefore has no
default and is the user's estimate.

The model holds only for thin coronas, with scale length L ≪ the critical radius
``r_c``. Against ray tracing through a spherical exponential corona, the absorption of an
oblique ray approaches ``η_a`` only as ``L/r_c`` → 0. At sin ψ = 0.9 it is 4% off at
``L/r_c`` = 0.0025 and 13% off at 0.01, and at 0.1 it is 0.33 against ``η_a`` = 0.18.
Rays that miss the critical sphere also still absorb. On a 40-beam OMEGA PDD pointing
with ``L/r_c`` = 0.2, ``η_a`` on a hard sphere at ``r_c`` gives σ = 5.7% (l ≤ 30), where
the corona gives 10.4% (research notes, "Corona model prototype: impact-parameter kernels
vs 3-D ray tracing"). Polar direct drive gets its equatorial drive from exactly those
grazing rays.

The unabsorbed light is dropped, and no reflection is followed. `element`, `backside` and
`per_beam` are therefore absorbed powers, and `delivered - sum(per_beam)` includes the
unabsorbed light. `η⊥` may be a dual number, and the result is then differentiable in it
too. `η⊥ = 1` gives the model of the default (``η_a ≡ 1``) and equals it to rounding;
the default itself is unchanged by this option.
"""
deposit(op::DepositionOperator, beams::AbstractVector{<:Beam};
        absorption::Union{Nothing,Real} = nothing) = _deposit(op, beams, nothing, absorption)

# Computes `deposit`, and with `split = (key, EG)` also the share of each beam group per
# element in the same pass. The element powers of beam `bi`, as values, are added to
# column `key[bi]` of `EG`, an elements × groups matrix, in the order `deposit` sums them.
# Beams with `key[bi] == 0` are not added. Each column is therefore exactly `deposit` of
# only the beams of that group, with the same `absorption`, since the beams at zero power
# add exact zeros. With `split` a `_GroupSplit`, the same holds for sets of groups in the
# number type of the deposition. The deposition itself is identical to that of
# `deposit`.
function _deposit(op::DepositionOperator, beams::AbstractVector{<:Beam}, split,
                  absorption::Union{Nothing,Real} = nothing)
    length(beams) == length(op.ports) || throw(DimensionMismatch("beam count changed"))
    for (b, p) in zip(beams, op.ports)
        b.port == p || throw(ArgumentError("beam $(b.name) port differs from the operator's"))
    end
    absorption === nothing || 0 < absorption <= 1 ||
        throw(ArgumentError("absorption must be in (0, 1], got $absorption"))
    T = _with_eltype(mapreduce(b -> typeof(b.power), promote_type, beams), absorption)
    t = op.target
    ne, nc = nelements(t), length(t.components)
    nb = length(beams)
    plan = op.plan
    nblock, nback = size(plan.boff, 1) - 1, length(plan.back)
    frames = [BeamFrame(b) for b in beams]
    screens = [_screen_frame(f, b) for (f, b) in zip(frames, beams)]
    slopes = [_cone_slopes(op.dirs[bi], frames[bi]) for bi in 1:nb]
    E = zeros(T, ne)
    PB = zeros(T, nb, nblock)     # per-beam power by element block
    LB = zeros(T, nc, nback)      # back-side power by piece
    SB = zeros(T, nback)          # piece totals
    # A beam at exactly zero power, in value and partials, deposits exact zeros and is
    # skipped. The other cone of a picket measure is such a beam.
    live = Bool[b.power !== zero(b.power) for b in beams]
    _deposit_units!(E, PB, LB, SB, op, frames, screens, slopes, absorption, split, live)
    # Reduce the partial sums in a fixed order.
    per_beam = T[sum(@view PB[bi, :]; init = zero(T)) for bi in 1:nb]   # typed: 1.10 infers Any here
    L = zeros(T, nc)
    for (u, (bi, _)) in enumerate(plan.back)
        per_beam[bi] += SB[u]
        L .+= @view LB[:, u]
    end
    return Deposition(E, L, per_beam, convert(T, sum(b.power for b in beams)))
end

# Returns the element type of a deposition, which is the beam type `T` promoted with the
# type of a dual `absorption`.
_with_eltype(T, ::Nothing) = T
_with_eltype(T, η) = promote_type(T, typeof(η))

# Returns the power absorbed by a ray at incidence cos ψ = `c`, per unit of incident flux
# along it. It is `c` for full absorption, or `η_a(ψ) c` with `η_a` of MN17 Eq. 29 (see
# `deposit`). With `nothing`, the expressions that use it are exactly as they are without
# the option.
_absorbed(::Nothing, c) = c
_absorbed(η, c) = (1 - (1 - η)^(c^3)) * c

# This function is the barrier for the threaded loop over the work units of `op.plan`.
function _deposit_units!(E, PB, LB, SB, op::DepositionOperator, frames, screens, slopes,
                         η = nothing, split = nothing, live = nothing)
    plan = op.plan
    nblock = size(plan.boff, 1) - 1
    _tforeach(length(plan.order)) do k
        u = plan.order[k]
        if u <= nblock
            # In an element block every point is on an element, so E is never indexed
            # below zero.
            for bi in axes(plan.boff, 2)
                r = plan.boff[u, bi]:plan.boff[u+1, bi]-1
                (isempty(r) || live !== nothing && !live[bi]) && continue
                PB[bi, u] = split === nothing ?
                            _deposit_points!(E, E, op, bi, frames[bi], screens[bi], slopes[bi], r, η) :
                            _deposit_points_split!(E, split, op, bi, frames[bi], screens[bi], slopes[bi], r, η)
            end
        else
            v = u - nblock
            bi, r = plan.back[v]
            live === nothing || live[bi] || return
            SB[v] = _deposit_points!(E, view(LB, :, v), op, bi, frames[bi], screens[bi],
                                     slopes[bi], r, η)
        end
    end
end

# Holds the element powers of a deposition by sets of beam groups, in the number type of
# the deposition, for the split of a `Scenario` in `evaluate`. `element[j]` sums the beams
# whose group is in `sets[j]`, which is sorted and unique, and the sets in `cols[bi]` hold
# beam `bi`. The element powers of each set are summed in the order `deposit` sums them.
# `element[j]` therefore equals exactly `deposit` of the beams of that set with the others
# at zero power, since the zero-power beams are skipped. The deposition itself is
# unchanged.
struct _GroupSplit{T}
    sets::Vector{Vector{Symbol}}
    cols::Vector{Vector{Int}}
    element::Vector{Vector{T}}
end

function _GroupSplit(sets::Vector{Vector{Symbol}}, beams::AbstractVector{<:Beam}, ne::Integer, ::Type{T}) where {T}
    cols = [findall(s -> b.group in s, sets) for b in beams]
    return _GroupSplit{T}(sets, cols, [zeros(T, ne) for _ in sets])
end

# Returns a beam-group set in the form that `_GroupSplit` uses as a key.
_group_set(gs) = sort!(unique!(collect(Symbol, gs)))

# Computes `deposit` of several sets of beams on one operator. The sets differ only in
# the beam powers; beam by beam they share the aim, focus, spot and port. The phases of a
# `PhasedScenario` on one mesh are such sets (`_deposit_groups`). The group sets of set k
# are kept in `splits[k]`, and a `_GroupSplit` with no sets keeps none.
#
# The deposition of a beam is its power times that of the beam at unit power, because the
# shape, cut-off and screens of the profile do not depend on the power. The intensity of
# each point is therefore evaluated once, at unit power, and every set adds it times the
# power of its beam. The deposition of set k then equals `deposit(op, sets[k])` to
# rounding, since the power enters as a last factor and not through the profile peak. A
# beam at exactly zero power in a set, in value and partials, adds nothing to that set,
# as `deposit` skips it.
function _deposit_sets(op::DepositionOperator, sets::Vector{Vector{Beam{T}}},
                       splits::Vector{_GroupSplit{T}}) where {T}
    b1 = first(sets)
    for beams in sets
        length(beams) == length(op.ports) || throw(DimensionMismatch("beam count changed"))
    end
    for (b, p) in zip(b1, op.ports)
        b.port == p || throw(ArgumentError("beam $(b.name) port differs from the operator's"))
    end
    t = op.target
    ne, nc, nb, nk = nelements(t), length(t.components), length(b1), length(sets)
    plan = op.plan
    nblock, nback = size(plan.boff, 1) - 1, length(plan.back)
    P = [T[b.power for b in beams] for beams in sets]           # per set, per beam
    frames = [BeamFrame(remake(b; power = one(T))) for b in b1]
    screens = [_screen_frame(f, b) for (f, b) in zip(frames, b1)]
    slopes = [_cone_slopes(op.dirs[bi], frames[bi]) for bi in 1:nb]
    live = [Int[k for k in 1:nk if P[k][bi] !== zero(T)] for bi in 1:nb]   # per beam, its sets with power
    E = [zeros(T, ne) for _ in 1:nk]
    PB = zeros(T, nb, nblock)                 # per beam at unit power, by element block
    LB = [zeros(T, nc, nback) for _ in 1:nk]  # per set, back-side power by piece
    SB = zeros(T, nback)                      # piece totals at unit power
    _deposit_sets_units!(E, PB, LB, SB, op, frames, screens, slopes, P, splits, live)
    unit = T[sum(@view PB[bi, :]; init = zero(T)) for bi in 1:nb]
    for (u, (bi, _)) in enumerate(plan.back)
        unit[bi] += SB[u]
    end
    out = Vector{Deposition{T}}(undef, nk)
    for k in 1:nk
        L = zeros(T, nc)
        for u in 1:nback
            L .+= @view LB[k][:, u]
        end
        out[k] = Deposition(E[k], L, T[P[k][bi] * unit[bi] for bi in 1:nb],
                            convert(T, sum(b.power for b in sets[k])))
    end
    return out
end

# This function is the barrier for the threaded loop of `_deposit_sets` over the work
# units of `op.plan`. Each unit writes only its elements of every set and of their group
# sets, or only its back-side piece.
function _deposit_sets_units!(E::Vector{Vector{T}}, PB::Matrix{T}, LB::Vector{Matrix{T}}, SB::Vector{T},
                              op::DepositionOperator, frames, screens, slopes, P::Vector{Vector{T}},
                              splits::Vector{_GroupSplit{T}}, live::Vector{Vector{Int}}) where {T}
    plan = op.plan
    nblock = size(plan.boff, 1) - 1
    _tforeach(length(plan.order)) do k
        u = plan.order[k]
        if u <= nblock
            for bi in axes(plan.boff, 2)
                r = plan.boff[u, bi]:plan.boff[u+1, bi]-1
                (isempty(r) || isempty(live[bi])) && continue
                PB[bi, u] = _deposit_points_sets!(E, splits, P, live[bi], op, bi, frames[bi], screens[bi],
                                                  slopes[bi], r)
            end
        else
            v = u - nblock
            bi, r = plan.back[v]
            isempty(live[bi]) && return
            ks = live[bi]
            SB[v] = _sum_points(op, bi, frames[bi], screens[bi], slopes[bi], r, T) do d, val
                for j in ks
                    @inbounds LB[j][-d, v] += P[j][bi] * val
                end
            end
        end
    end
end

# Deposits the points of beam `bi` in an element block, all on elements, at unit power.
# Each point is added, times the beam power of the set, to the element powers of every
# set in `ks` that gives the beam power, and to those of the group sets of that set that
# hold the beam. Returns the unit-power sum.
function _deposit_points_sets!(E::Vector{Vector{T}}, splits::Vector{_GroupSplit{T}}, P::Vector{Vector{T}},
                               ks::Vector{Int}, op::DepositionOperator, bi, f::BeamFrame,
                               fv::BeamFrame{Float64}, slopes, r) where {T}
    return _sum_points(op, bi, f, fv, slopes, r, T) do d, val
        for k in ks
            @inbounds v = P[k][bi] * val
            @inbounds E[k][d] += v
            @inbounds s = splits[k]
            @inbounds for j in s.cols[bi]
                s.element[j][d] += v
            end
        end
    end
end

# Deposits the points of beam `bi` in an element block, all on elements, and also adds
# them to the element powers of the group sets that hold the beam.
function _deposit_points_split!(E::AbstractVector{T}, s::_GroupSplit{T}, op::DepositionOperator,
                                bi, f::BeamFrame, fv::BeamFrame{Float64}, slopes, r, η = nothing) where {T}
    cs = s.cols[bi]
    isempty(cs) && return _deposit_points!(E, E, op, bi, f, fv, slopes, r, η)
    G = s.element
    return _sum_points(op, bi, f, fv, slopes, r, T, η) do d, val
        @inbounds E[d] += val
        for j in cs
            @inbounds G[j][d] += val
        end
    end
end

# Deposits the points of beam `bi` in an element block, all on elements, and also adds
# them to the column of its group.
function _deposit_points_split!(E::AbstractVector{T}, (key, EG), op::DepositionOperator,
                                bi, f::BeamFrame, fv::BeamFrame{Float64}, slopes, r, η = nothing) where {T}
    g = key[bi]
    return _sum_points(op, bi, f, fv, slopes, r, T, η) do d, val
        @inbounds E[d] += val
        g > 0 && (@inbounds EG[d, g] += ForwardDiff.value(val))
    end
end

# Returns the transverse slopes of the cone directions in the beam frame, which the port
# fixes.
function _cone_slopes(ds, f::BeamFrame)
    κ = [dot(d, f.k) for d in ds]
    α = [dot(d, f.e1) / dot(d, f.k) for d in ds]
    β = [dot(d, f.e2) / dot(d, f.k) for d in ds]
    return (κ, α, β)
end

# Returns the Float64 frame used to screen points and cone directions (`_beyond`,
# `_live_directions`).
_screen_frame(f::BeamFrame{Float64}, b) = f
_screen_frame(::BeamFrame, b) = BeamFrame(_value_beam(b))

# Deposits the points `r` of beam `bi` into the elements `Ec` and back sides `Lc` and
# returns their sum. Most points of a beam lie off its spot, and a bound on all cone
# directions at once skips them. For dual beams, each remaining direction is screened on
# values. Skipped points and directions would contribute exact zeros.
function _deposit_points!(Ec::AbstractVector{T}, Lc::AbstractVector{T}, op::DepositionOperator,
                          bi, f::BeamFrame, fv::BeamFrame{Float64}, slopes, r, η = nothing) where {T}
    return _sum_points(op, bi, f, fv, slopes, r, T, η) do d, val
        @inbounds d > 0 ? (Ec[d] += val) : (Lc[-d] += val)
    end
end

# Runs the deposition loop. It calls `sink(destination, power)` for every point of beam
# `bi` in `r` that receives light (`destination` > 0 is an element, < 0 the back side of
# a component) and returns the sum of those powers (of type `T`); `η` is the `absorption`
# of `deposit`. The coordinates of a point in the beam frame are taken once, on values
# (for the screens). Those of a dual beam are built from them with the partials that
# every point shares (`_cone_coords`).
@inline function _sum_points(sink::F, op::DepositionOperator, bi, f::BeamFrame,
                             fv::BeamFrame{Float64}, (κ, α, β), r, ::Type{T}, η = nothing) where {F,T}
    prof = op.profile
    prof === nothing ||
        return _sum_points_visrad(sink, op, bi, _VRFrame(f, op.envelope[bi], prof.taper),
                                  _VRFrame(fv, op.envelope[bi], prof.taper), r, T, η)
    ds = op.dirs[bi]
    P, W, Nr, M, D = op.pts[bi], op.w[bi], op.normal[bi], op.mask[bi], op.dest[bi]
    αmax, βmax = maximum(abs, α), maximum(abs, β)
    o = _origin_coords(f)
    acc = zero(T)
    @inbounds for q in r
        x = P[q]
        suv = _cone_coords(fv, x)
        _beyond(fv, αmax, βmax, suv) && continue
        m = T === Float64 ? M[q] : _live_directions(fv, α, β, M[q], suv)
        m == 0 && continue
        val = W[q] * _cone_intensity(f, ds, α, β, κ, m, _cone_coords(f, o, x, suv), Nr[q], η)
        sink(D[q], val)
        acc += val
    end
    return acc
end

# Returns the coordinates (s, u₀, v₀) of a point in a beam frame, with s measured along
# the beam from best focus and u₀ and v₀ across it.
@inline function _cone_coords(f::BeamFrame, x)
    d = x - f.focus
    return dot(d, f.k), dot(d, f.e1), dot(d, f.e2)
end

# Returns the same coordinates in the frame `f` from those in its Float64 frame, `suv`,
# which is built from the beam values. On duals the result takes those values with the
# partials of the coordinates `o` of the origin (`_origin_coords`). Every point shares
# these partials, because the axes of the frame carry none, so three dual dot products per
# point are avoided.
@inline _origin_coords(f::BeamFrame) = _cone_coords(f, zero(Vec3))
@inline _cone_coords(f::BeamFrame{Float64}, o, x, suv) = suv
@inline function _cone_coords(f::BeamFrame{ForwardDiff.Dual{G,Float64,N}}, o, x, suv) where {G,N}
    D = ForwardDiff.Dual{G,Float64,N}
    return D(suv[1], ForwardDiff.partials(o[1])), D(suv[2], ForwardDiff.partials(o[2])),
           D(suv[3], ForwardDiff.partials(o[3]))
end
@inline _cone_coords(f::BeamFrame, o, x, suv) = _cone_coords(f, x)

# Checks whether the ray of every cone direction through `x` crosses the focal plane
# outside the profile cutoff. The point can also be given by its coordinates
# `(s, u₀, v₀)` (`_cone_coords`). The crossing (u₀ - sαⱼ, v₀ - sβⱼ) is at least
# |u₀| - |s| max|α| and |v₀| - |s| max|β| off axis, with a relative allowance for
# rounding.
@inline _beyond(f::BeamFrame{Float64}, αmax, βmax, x) = _beyond(f, αmax, βmax, _cone_coords(f, x))
@inline function _beyond(f::BeamFrame{Float64}, αmax, βmax, (s, u0, v0)::NTuple{3,Float64})
    s = abs(s)
    u = max(abs(u0) - s * αmax, 0.0)
    v = max(abs(v0) - s * βmax, 0.0)
    _outside(c) = u^2 / c.a2 + v^2 / c.b2 > c.qcut * (1 + 1e-9)
    return _outside(f.core) && (!f.has_halo || _outside(f.halo))
end

# Returns the visible directions `m` whose ray through a point, given by its coordinates
# `(s, u₀, v₀)`, crosses the focal plane inside the profile cutoff. The test uses the beam
# values. The other directions contribute exact zeros in value and partials, so dual
# arithmetic is spent only on the ~1/6 of point-direction pairs under the spot. A Float64
# beam is not screened, because `_cone_intensity` tests the cutoff itself.
@inline function _live_directions(f::BeamFrame, α, β, m::UInt16, (s, u0, v0))
    live = UInt16(0)
    for j in eachindex(α)
        (m >> (j - 1)) & 0x1 == 0x1 || continue
        u = u0 - s * α[j]
        v = v0 - s * β[j]
        (_sg_inside(f.core, u, v) || (f.has_halo && _sg_inside(f.halo, u, v))) &&
            (live |= UInt16(1) << (j - 1))
    end
    return live
end

# Deposits a beam with VisRad's profile. Each stored point sees the beam axis, the one
# direction of the operator, and absorbs the profile at its distance from focus times
# |n·k|. Points outside the envelope, tested on values, contribute exact zeros. The
# points are screened in three steps:
# - Each run of `_TILE` points of a beam is skipped when a bound without logarithms puts
#   its ball beyond the cut-off.
# - Each remaining point is skipped when the same bound puts it beyond the cut-off
#   (`_vr_far`). The figures for this bound on case B of the research notes are 99% of
#   the points beyond the cut and 90% of the points left lit.
# - The profile law is evaluated for the rest and is exactly zero, in value and partials,
#   at the points beyond its cut. Screening those points on values first would cost more
#   than the few evaluations it would skip.
# The coordinates of a point are taken once, on values. Those of a dual beam carry the
# partials that every point shares (`_vr_coords(f, x, uvs)`).
function _sum_points_visrad(sink::F, op::DepositionOperator, bi, f::_VRFrame, fv::_VRFrame{Float64},
                            r, ::Type{T}, η = nothing) where {F,T}
    P, W, Nr, D = op.pts[bi], op.w[bi], op.normal[bi], op.dest[bi]
    tiles = op.plan.tiles[bi]
    k = f.k
    acc = zero(T)
    isempty(r) && return acc
    @inbounds for t in fld(first(r) - 1, _TILE)+1:fld(last(r) - 1, _TILE)+1
        _vr_far_ball(fv, tiles[t]) && continue
        for q in max(first(r), _TILE * (t - 1) + 1):min(last(r), _TILE * t)
            x = P[q]
            uvs = _vr_coords(fv, x)
            _vr_far(fv, uvs...) && continue
            val = W[q] * _vr_intensity(f, _vr_coords(f, x, uvs)...) * _absorbed(η, abs(dot(Nr[q], k)))
            sink(D[q], val)
            acc += val
        end
    end
    return acc
end

_deposit_beam!(Ec, Lc, op::DepositionOperator, bi, f::BeamFrame) =
    _deposit_points!(Ec, Lc, op, bi, f, _screen_frame(f, nothing), _cone_slopes(op.dirs[bi], f),
                     eachindex(op.w[bi]))

# Returns the absorbed power per unit area at `x`, with lit-side normal `n`. It is the mean
# over the visible cone directions dⱼ of the best-focus profile where the ray along dⱼ
# through `x` crosses the focal plane, times the incidence factor |n·dⱼ|/(k·dⱼ) of the
# ray. The sub-beam of each direction is a bundle of parallel rays. The visibility mask
# that selects it is therefore exactly consistent with where its light goes, and power is
# conserved at shadow edges. Giving each sub-beam a share of the f-cone spread would smooth
# the footprint but break that consistency, and 0.6% of a beam clipped by an LEH edge
# would be lost. With `η`, the `absorption` of `deposit`, the incidence factor of each
# direction includes its absorbed fraction η_a. `(s, u₀, v₀)` are the coordinates of the
# point (`_cone_coords`).
@inline function _cone_intensity(f::BeamFrame, ds, α, β, κ, m::UInt16, (s, u0, v0), n, η = nothing)
    acc = zero(_absorbed(η, f.core.power_norm * s))
    for j in eachindex(ds)
        (m >> (j - 1)) & 0x1 == 0x1 || continue
        u = u0 - s * α[j]                     # crossing of the best-focus plane
        v = v0 - s * β[j]
        I = _sg_focal(f.core, u, v)
        f.has_halo && (I += _sg_focal(f.halo, u, v))
        acc += I * (_absorbed(η, abs(dot(n, ds[j]))) / κ[j])
    end
    return acc / length(ds)
end

"""
    laser_flux(dep, target) -> Vector

Returns the absorbed laser intensity per element (TW/mm²; multiply by 1e14 for W/cm²).
"""
laser_flux(d::Deposition, t::Target) = d.element ./ t.elem_area

"""
    component_power(dep, target; by = :part) -> Dict{Symbol,Float64}

Returns the absorbed laser power (TW) per part, plus back-side (outside) losses keyed as
`Symbol(name, "_outside")` and `:missed` for power that strikes nothing. The pieces of
a part split for computation (`curved_components`) are summed into the part; `by =
:component` keys by component instead. The two agree on targets of unsplit parts.
"""
function component_power(d::Deposition, t::Target; by::Symbol = :part)
    by in (:part, :component) || throw(ArgumentError("by must be :part or :component"))
    out = Dict{Symbol,Float64}()
    for (ci, c) in enumerate(t.components)
        c.kind === :aperture && continue
        n = by === :part ? c.part : c.name
        out[n] = get(out, n, 0.0) + sum(ForwardDiff.value(d.element[e]) for e in t.comp_elems[ci]; init = 0.0)
        bs = ForwardDiff.value(d.backside[ci])
        bs > 0 && (k = Symbol(n, :_outside); out[k] = get(out, k, 0.0) + bs)
    end
    out[:missed] = ForwardDiff.value(d.delivered - sum(d.per_beam))
    return out
end

# ---------------------------------------------------------------------------------------
# Ray trace (reference)
# ---------------------------------------------------------------------------------------

"""
    raytrace_deposition(target, beams; nspot = 41, ncone = 13, extent = nothing) -> Deposition

Computes the reference deposition by geometric ray tracing. Each spot component (core, and
halo if present) is sampled with an `nspot × nspot` grid of best-focus points spanning
`±extent` 1/e semi-axes (by default out to 1e-12 of the component peak), each launched along
`ncone` f-cone directions and weighted by the focal intensity. Rays travel straight (no
refraction) to the first opaque surface. Power is conserved exactly, and defocus is exact
within the sampling. The result is not differentiable.
"""
function raytrace_deposition(t::Target, beams::AbstractVector{<:Beam}; nspot::Integer = 41,
                             ncone::Integer = 13, extent = nothing)
    ne, nc, nb = nelements(t), length(t.components), length(beams)
    E = [zeros(ne) for _ in 1:nb]
    L = [zeros(nc) for _ in 1:nb]
    per_beam = zeros(nb)
    _tforeach(nb) do bi
        per_beam[bi] = _raytrace_beam!(E[bi], L[bi], t, beams[bi], nspot, ncone, extent)
    end
    return Deposition(sum(E), sum(L), per_beam, sum(ForwardDiff.value(b.power) for b in beams))
end

"""
    raytrace_component_power(target, beams; nspot = 41, ncone = 13, extent = nothing) -> Matrix

Returns the power (TW) that each beam (rows) deposits on each component (columns, both
sides), from the reference ray trace. The power accounting is exact and independent of
mesh resolution.
"""
function raytrace_component_power(t::Target, beams::AbstractVector{<:Beam}; nspot::Integer = 41,
                                  ncone::Integer = 13, extent = nothing)
    nc = length(t.components)
    M = zeros(length(beams), nc)
    _tforeach(length(beams)) do bi
        E = zeros(nelements(t))
        L = zeros(nc)
        _raytrace_beam!(E, L, t, beams[bi], nspot, ncone, extent)
        for e in eachindex(E)
            M[bi, t.elem_comp[e]] += E[e]
        end
        M[bi, :] .+= L
    end
    return M
end

# Each spot component (core, halo) is sampled on a separate nspot × nspot grid spanning
# `extent` 1/e radii, by default out to where the component falls to 1e-12 of its peak.
function _raytrace_beam!(E, L, t::Target, b::Beam, nspot, ncone, extent)
    val = ForwardDiff.value
    sp = b.spot
    a, bb = val(sp.a), val(sp.b)
    e1, e2 = spot_axes(b)
    k = direction(b)
    focus = Vec3(val.(b.aim)) + val(b.defocus) * k
    dirs, dw = cone_directions(b, ncone)
    far = 4 * (t.scale + norm(focus)) + 1
    P = val(b.power)
    comps = [(1 - sp.halo, 1.0, sp.order)]
    sp.halo > 0 && push!(comps, (sp.halo, sp.halo_scale, sp.halo_order))
    total = 0.0
    for (frac, sc, n) in comps
        ext = extent === nothing ? 1.02 * log(1e12)^(1 / n) : extent
        us = range(-ext, ext; length = nspot)
        wsum = sum(exp(-(u^2 + w^2)^(n / 2)) for u in us, w in us)
        for u in us, w in us
            wp = frac * exp(-(u^2 + w^2)^(n / 2)) / wsum
            wp < 1e-15 && continue
            pf = focus + (sc * u * a) * e1 + (sc * w * bb) * e2
            for (d, wd) in zip(dirs, dw)
                pw = P * wp * wd
                _, tri = closest_hit(t.bvh, pf - far * d, d)
                tri == 0 && continue
                total += pw
                e = dot(d, t.tri_normal[tri]) < 0 ? t.tri_front[tri] : t.tri_back[tri]
                e > 0 ? (E[e] += pw) : (L[t.tri_comp[tri]] += pw)
            end
        end
    end
    return total
end
