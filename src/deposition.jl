# Laser deposition: a differentiable quadrature operator and an energy-conserving ray trace.

# Run `f(i)` for `i in 1:n` on all threads, handing out indices one at a time in order
# (dynamic load balancing; with heterogeneous cores a static split waits on the slowest).
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

# Each beam's points are sorted by destination (elements ascending, then back sides by
# component). `deposit` splits the work into units that own their outputs: a block of
# consecutive elements takes every beam's points in it (and is the only writer of those
# elements), and each beam's back-side points are cut into pieces with their own partial
# sums. Units are handed to threads dynamically, largest first; every output is summed
# in a fixed order (an element over beams, then points, as a serial loop would), so the
# result is bitwise independent of thread count and scheduling.
struct _DepositPlan
    boff::Matrix{Int}                     # (nblock+1) × nb: first point of each element block; last row: first back-side point
    back::Vector{Tuple{Int,UnitRange{Int}}}  # back-side pieces (beam, points)
    order::Vector{Int}                    # units by decreasing size; ≤ nblock: element block, else back-side piece
end

"""
    DepositionOperator(target, beams; level = 1, ncone = 13, resolution = 0.3,
                       reach = 2.0, maxlevel = 5, shadow = 0,
                       backside_resolution = 1.0, edge_resolution = 0.01,
                       profile = nothing)

Precomputed laser-deposition quadrature. For each beam, every opaque triangle side that
faces the beam and belongs to an element is sampled at `4^level` points, or more where
the beam can land: a triangle within `reach` (mm) of the beam's footprint (2.5 spot
semi-axes grown by the f-cone along the beam axis through the aim) is divided until its
sub-triangles are no longer than `resolution` times the spot's smaller semi-axis, up to
`4^maxlevel` points; with `shadow > 0`, such a triangle crossed by a shadow edge (its
points see different cone directions) is divided `shadow` times more. `reach` is the aim
motion the operator stays accurate for; beyond it, deposition is resolved only at
`level`. The defaults conserve each beam's power to 1e-4 when nothing clips and to 2e-3
when an entrance-hole lip clips several percent of it, on fine and coarse hohlraum
meshes alike.

Unmodeled back sides (the outside of a hohlraum, whose power is the clipped `backside`
loss) need only their total per component, and that is set by the open edges the light
passes: an entrance-hole lip clipping 1e-4 of a beam puts the clipped light within
~1/50 of a spot semi-axis of the rim. Back sides are sampled at one point per triangle
away from the beam; within `reach`,
- a triangle on an open edge (an edge no other opaque triangle shares, such as an
  entrance-hole rim) or touching one at a vertex gets a Gauss rule graded geometrically
  toward it, from `resolution` down to `edge_resolution` spot semi-axes;
- a triangle within one (larger) spot semi-axis of an open edge, where the shoulder of a
  clipped spot lands, is divided down to `resolution`;
- the rest are divided down to `backside_resolution`.
On the NIF hohlraum, clipped power then matches ray traces (averaged over nspot
241–361) to 1% from 1e-4 to 5% of a beam, 2.5% at 1e-6, and 0.6% up to 50%, and is
smooth in the aim. Element values do not depend on these two keywords.

Each point stores its area weight, the normal of the lit side and a bitmask of which of
the `ncone` f-cone directions reach it from the lens without hitting another surface.
These depend only on the target and the beam *directions*, so the operator stays valid
as pointing, defocus, spot size, and power change — deposition is then a smooth
function of those parameters (see `deposit`).

The beam is the superposition of `ncone` collimated sub-beams, one per f-cone
direction `dⱼ` (`cone_directions`), each carrying the best-focus profile. A point
receives from each visible direction the focal profile where the ray along `dⱼ`
through it crosses the best-focus plane, times that ray's incidence factor
`|n·dⱼ|/(k·dⱼ)`. This is the geometric-optics convolution of the focal spot with the
f-cone, sampled at the cone directions as the reference ray trace samples it, so the
power on any surface that intercepts the whole beam is the beam power, at any
distance from focus and any incidence; and shadows cast by objects near focus — a
capsule at TCC, an LEH edge — are correlated with the part of the spot they block.

`profile = VisRadProfile(...)` replaces this with VisRad's beam model, so that a
design can be checked in VisRad and its deposition reproduced. The model keeps the
super-Gaussian's shape along the beam, with semi-axes that grow on a hyperbola and an
envelope cut-off. Points are then sampled along the beam axis only (`ncone` is
ignored), and the beam axis through each point decides its shadowing. The two models
agree near focus; away from it VisRad's spot keeps its shape while the f-cone
convolution softens the edges, and VisRad's spot grows faster (on NIF hohlraum walls
the f-cone model's peaks are 10% higher than VisRad's). Power on a surface that
intercepts the whole beam is the beam's power in planes normal to the beam, and only
approximately on tilted ones (per-plane normalization: 0.15% high on a foil tilted 20°
across ±2.5 mm of focus, where VisRad's totals are exact).

The beams' ports must not change after construction.
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

function _deposit_plan(ne::Integer, dest::Vector{Vector{Int32}}; unit::Integer = _DEPOSIT_UNIT)
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
    return _DepositPlan(boff, back, sortperm(cost; rev = true))
end

DepositionOperator(t::Target, ports, dirs, pts, w, normal, mask, dest,
                   profile::Union{Nothing,VisRadProfile} = nothing, envelope = Float64[]) =
    DepositionOperator(t, ports, dirs, pts, w, normal, mask, dest, _deposit_plan(nelements(t), dest),
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
    # sub-triangle quadratures in barycentric form, shared by all triangles
    rules = [_bary_quadrature(l) for l in 0:maxlv]
    ε = _ray_eps(t)
    geo = [_operator_beam(b, ncone, profile) for b in beams]
    flags, rimdist = _open_boundary(t)
    back = (; resolution, coarse = backside_resolution, edge = edge_resolution, flags, rimdist)
    # Work units are (beam, triangle range): beams differ in cost and the triangles near
    # a beam's footprint are refined, so units are small and handed out dynamically.
    # Each unit's points are stored apart and joined in order: the result does not
    # depend on scheduling.
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

# Beam quantities the operator's sampling depends on (Float64 values).
function _operator_beam(b::Beam, ncone, profile)
    k = direction(b)
    ds, _ = cone_directions(b, ncone)
    sp = b.spot
    amin = ForwardDiff.value(min(sp.a, sp.b))
    # footprint: 2.5 semi-axes grown by the f-cone; VisRad's envelope grows as a(z) does,
    # within ρ (a₀ + |z| a₀/z_A), a₀/z_A < 1/(2f#) × max(1, ρ/γ), or for a measured plate
    # within ρ times its steepest growth of either semi-axis
    env = profile === nothing ? 1.0 : _vr_envelope(profile, b)
    r = max(2.5, env)
    amax = ForwardDiff.value(max(sp.a, sp.b))
    extent = r * amax * (sp.halo > 0 ? max(1, sp.halo_scale) : 1)
    aim = Vec3(ForwardDiff.value.(b.aim))
    focus = aim + ForwardDiff.value(b.defocus) * k
    cone = (profile === nothing ? 1.0 : max(1.0, env) / sqrt(min(1.0, _vr_zfactor(profile, env)))) /
           (2b.fnumber)
    m = profile === nothing ? nothing : _vr_measured(profile, sp)
    if m !== nothing
        slope(y, r0) = maximum((y[i+1] - y[i]) / (m.z[i+1] - m.z[i]) for i in 1:length(m.z)-1) * r0
        cone = max(cone, 1.1 * env * max(slope(m.a_ratio, m.a), slope(m.b_ratio, m.b)))
    end
    return (; k, ds, amin, extent, focus, cone)
end

# Lit quadrature points of one beam on the triangles `tris`.
function _operator_points(t::Target, g, tris, rules, ε, level, maxlevel, resolution, reach, shadow,
                          back)
    (; k, ds, amin, extent, focus, cone) = g
    P, W, N, M, D = Vec3[], Float64[], Vec3[], UInt16[], Int32[]
    # rays toward the lens, and the last leaf that blocked each
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
            # unmodeled back side: only its total matters (see `_backside_points!`)
            lv = 0
            if near
                _backside_points!(P, W, N, M, D, t, tri, back, amin, side_n, target_id, ds, ε, cache) &&
                    continue
                # within a spot semi-axis of an open edge: the shoulder of a clipped spot lands here
                res = back.rimdist[tri] <= extent / 2.5 ? back.resolution : back.coarse
                lv = clamp(ceil(Int, log2(L / (res * amin))), 0, maxlevel)
            end
        elseif near && L > resolution * amin * 2^level
            lv = clamp(ceil(Int, log2(L / (resolution * amin))), level, maxlevel)
        end
        xs, ms = _visible_directions(t, rules[lv+1][1], v0, e1, e2, side_n, ds, ε, cache)
        # a shadow edge crosses the triangle: resolve it more finely
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

# Quadrature points of a triangle and, for each, the bitmask of cone directions that
# reach it unobstructed: those facing the lit side whose ray toward the lens hits
# nothing (`any_hit`, asked of all points and directions as one packet; `cache` holds
# the reversed directions and keeps each one's last blocking leaf across a work unit).
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

# Cone directions that strike the side with normal `n` (d·n < 0, summed in a fixed order
# so the answer for a direction in the plane does not depend on how `dot` is compiled).
function _facing_mask(ds, n)
    m = UInt16(0)
    for (j, d) in enumerate(ds)
        (d[1] * n[1] + d[2] * n[2]) + d[3] * n[3] < 0 && (m |= UInt16(1) << (j - 1))
    end
    return m
end

# Open boundary of the opaque surfaces: triangle edges that no other opaque triangle
# shares (an entrance hole's rim, a foil's edge), and the vertices on them. Per triangle,
# bit j-1 flags edge j (vertex j to j+1; vertices v0, v0+e1, v0+e2) and bit j+2 vertex j.
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
    # distance of each triangle from the open boundary (lower bound)
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

# Panel ends on [0, 1]: `n` equal panels, the first (`lo`) and/or last (`hi`) divided
# geometrically (ratio `q`) toward that end down to width `wmin`.
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

# Back-side points of a triangle on the open boundary. The light a back side intercepts
# at an open edge (an entrance-hole lip clipping a beam's tail) lies within a few
# e-folding lengths of the tail of the edge: ~1/50 of a spot semi-axis at 1e-4 clipped,
# far finer than spot-resolving sub-triangles, whose centroids then miss most of it
# (they underestimated 1e-4 clipping by 10–60%). The triangle is mapped from the unit
# square, collapsing onto one vertex C (u from C to the opposite edge, v across), with
# the coordinates that reach open edges and vertices graded geometrically (ratio 1/2)
# toward them, down to `edge` spot semi-axes, and two Gauss points per panel. A vertex
# on the rim needs this as much as an edge: the triangles fanning out from it hold a
# share ~λ/L of the rim's light (λ the tail length, L the edge). Panels are `resolution`
# spot semi-axes across the triangle and `coarse` along an open edge. Returns false (no
# points added) for a triangle off the open boundary.
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
    # two-point Gauss panels in u; in v, as many as the u panel's width across needs
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

Absorbed laser power by element (`element`, TW), power striking the unmodeled back side
of each component (`backside`, TW — e.g. the outside of a hohlraum, i.e. LEH clipping),
and per-beam totals (`per_beam`, TW, power that landed on any opaque surface).
"""
struct Deposition{T}
    element::Vector{T}
    backside::Vector{T}
    per_beam::Vector{T}
    delivered::T
end

"""
    deposit(op, beams; absorption = nothing) -> Deposition

Evaluate laser deposition for `beams` (same ports and order as at construction). The
result is differentiable in every beam's `aim`, `power`, `defocus` and spot axes.
Power that strikes no opaque surface (e.g. passes through a hohlraum and out the
other LEH) is `delivered - sum(per_beam)`.
Threads share the work dynamically, and every sum is taken in a fixed order: the result
is the same, bit for bit, on any number of threads.

By default every ray deposits all its power where it lands. `absorption = η⊥`
(0 < η⊥ ≤ 1) keeps only the fraction

    η_a(ψ) = 1 - (1 - η⊥)^(cos³ψ)

of each ray, with ψ its angle of incidence (each f-cone direction's, or the beam axis's
with a VisRad profile). This is inverse-bremsstrahlung absorption in a thin corona whose
density falls exponentially with radius, η⊥ = 1 - exp(-8ν_ei L/3c) being the absorption
at normal incidence (ν_ei the collision frequency at critical density, L the scale
length): Murakami & Nishi, Matter Radiat. Extremes 2, 55 (2017), Eq. 29, after Kruer's
text. It is a model for hard-sphere direct-drive illumination, derived for that profile,
not a measured or calibrated absorption: η⊥ has no default and is the user's estimate.
The unabsorbed light is dropped (no reflection is followed), so `element`, `backside`
and `per_beam` are absorbed powers and `delivered - sum(per_beam)` includes it. `η⊥` may
be a dual number (the result is differentiable in it too). `η⊥ = 1` is the default's
model (η_a ≡ 1), equal to it to rounding; the default itself is unchanged by this option.
"""
function deposit(op::DepositionOperator, beams::AbstractVector{<:Beam};
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
    _deposit_units!(E, PB, LB, SB, op, frames, screens, slopes, absorption)
    # partial sums reduced in a fixed order
    per_beam = [sum(@view PB[bi, :]; init = zero(T)) for bi in 1:nb]
    L = zeros(T, nc)
    for (u, (bi, _)) in enumerate(plan.back)
        per_beam[bi] += SB[u]
        L .+= @view LB[:, u]
    end
    return Deposition(E, L, per_beam, convert(T, sum(b.power for b in beams)))
end

# Element type of a deposition: the beams' (`T`), promoted with a dual `absorption`'s.
_with_eltype(T, ::Nothing) = T
_with_eltype(T, η) = promote_type(T, typeof(η))

# The power a ray at incidence cos ψ = `c` absorbs, per unit of incident flux along it:
# `c` (all of it), or `η_a(ψ) c` with `η_a` of MN17 Eq. 29 (see `deposit`). `nothing`
# leaves the expressions it enters as they were without the option (bit for bit).
_absorbed(::Nothing, c) = c
_absorbed(η, c) = (1 - (1 - η)^(c^3)) * c

# Function barrier for the threaded loop over `op.plan`'s work units.
function _deposit_units!(E, PB, LB, SB, op::DepositionOperator, frames, screens, slopes,
                         η = nothing)
    plan = op.plan
    nblock = size(plan.boff, 1) - 1
    _tforeach(length(plan.order)) do k
        u = plan.order[k]
        if u <= nblock
            # an element block: every point lands on an element (E is never indexed < 0)
            for bi in axes(plan.boff, 2)
                r = plan.boff[u, bi]:plan.boff[u+1, bi]-1
                isempty(r) || (PB[bi, u] = _deposit_points!(E, E, op, bi, frames[bi], screens[bi],
                                                            slopes[bi], r, η))
            end
        else
            v = u - nblock
            bi, r = plan.back[v]
            SB[v] = _deposit_points!(E, view(LB, :, v), op, bi, frames[bi], screens[bi],
                                     slopes[bi], r, η)
        end
    end
end

# Transverse slopes of the cone directions in the beam frame (fixed by the port).
function _cone_slopes(ds, f::BeamFrame)
    κ = [dot(d, f.k) for d in ds]
    α = [dot(d, f.e1) / dot(d, f.k) for d in ds]
    β = [dot(d, f.e2) / dot(d, f.k) for d in ds]
    return (κ, α, β)
end

# Float64 frame for screening points and cone directions (`_beyond`, `_live_directions`).
_screen_frame(f::BeamFrame{Float64}, b) = f
_screen_frame(::BeamFrame, b) = BeamFrame(_value_beam(b))

# Deposit beam `bi`'s points `r` into elements `Ec` and back sides `Lc`; returns their sum.
# Most points of a beam lie off its spot: a bound on all cone directions at once skips
# them, and for dual beams each remaining direction is screened on values. Skipped
# points and directions would contribute exact zeros.
function _deposit_points!(Ec::AbstractVector{T}, Lc::AbstractVector{T}, op::DepositionOperator,
                          bi, f::BeamFrame, fv::BeamFrame{Float64}, slopes, r, η = nothing) where {T}
    return _sum_points(op, bi, f, fv, slopes, r, T, η) do d, val
        @inbounds d > 0 ? (Ec[d] += val) : (Lc[-d] += val)
    end
end

# The deposition loop: `sink(destination, power)` for every point of beam `bi` in `r`
# that receives light (`destination` > 0 an element, < 0 a component's back side), and
# the sum of those powers (of type `T`); `η` the `absorption` of `deposit`.
@inline function _sum_points(sink::F, op::DepositionOperator, bi, f::BeamFrame,
                             fv::BeamFrame{Float64}, (κ, α, β), r, ::Type{T}, η = nothing) where {F,T}
    op.profile === nothing ||
        return _sum_points_visrad(sink, op, bi, _VRFrame(f, op.profile, op.envelope[bi]),
                                  _VRFrame(fv, op.profile, op.envelope[bi]), r, T, η)
    ds = op.dirs[bi]
    P, W, Nr, M, D = op.pts[bi], op.w[bi], op.normal[bi], op.mask[bi], op.dest[bi]
    αmax, βmax = maximum(abs, α), maximum(abs, β)
    acc = zero(T)
    @inbounds for q in r
        _beyond(fv, αmax, βmax, P[q]) && continue
        m = T === Float64 ? M[q] : _live_directions(fv, α, β, M[q], P[q])
        m == 0 && continue
        val = W[q] * _cone_intensity(f, ds, α, β, κ, m, P[q], Nr[q], η)
        sink(D[q], val)
        acc += val
    end
    return acc
end

# Whether every cone direction's ray through `x` crosses the focal plane outside the
# profile's cutoff: the crossing (u₀ - sαⱼ, v₀ - sβⱼ) is at least |u₀| - |s| max|α| and
# |v₀| - |s| max|β| off axis (with a relative allowance for rounding).
@inline function _beyond(f::BeamFrame{Float64}, αmax, βmax, x)
    d = x - f.focus
    s = abs(dot(d, f.k))
    u = max(abs(dot(d, f.e1)) - s * αmax, 0.0)
    v = max(abs(dot(d, f.e2)) - s * βmax, 0.0)
    _outside(c) = u^2 / c.a2 + v^2 / c.b2 > c.qcut * (1 + 1e-9)
    return _outside(f.core) && (!f.has_halo || _outside(f.halo))
end

# The visible directions `m` whose ray through `x` crosses the focal plane inside the
# profile's cutoff, from the beam's values: the others contribute exact zeros (value
# and partials), so dual arithmetic is spent only on the ~1/6 of point-direction pairs
# under the spot. (A Float64 beam is not screened: `_cone_intensity` tests the cutoff
# itself.)
@inline function _live_directions(f::BeamFrame, α, β, m::UInt16, x)
    d = x - f.focus
    s = dot(d, f.k)
    u0 = dot(d, f.e1)
    v0 = dot(d, f.e2)
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

# VisRad's profile: each stored point sees the beam axis (the operator's one direction)
# and absorbs the profile at its own distance from focus times |n·k|. Points outside
# the envelope (on values) contribute exact zeros: a dual beam skips them on values;
# a Float64 one evaluates them (the profile is exactly zero there, and screening first
# would compute the semi-axes twice).
function _sum_points_visrad(sink::F, op::DepositionOperator, bi, f::_VRFrame, fv::_VRFrame{Float64},
                            r, ::Type{T}, η = nothing) where {F,T}
    P, W, Nr, D = op.pts[bi], op.w[bi], op.normal[bi], op.dest[bi]
    k = f.k
    acc = zero(T)
    @inbounds for q in r
        T === Float64 || !_vr_beyond(fv, P[q]) || continue
        val = W[q] * _vr_intensity(f, P[q]) * _absorbed(η, abs(dot(Nr[q], k)))
        sink(D[q], val)
        acc += val
    end
    return acc
end

_deposit_beam!(Ec, Lc, op::DepositionOperator, bi, f::BeamFrame) =
    _deposit_points!(Ec, Lc, op, bi, f, _screen_frame(f, nothing), _cone_slopes(op.dirs[bi], f),
                     eachindex(op.w[bi]))

# Absorbed power per unit area at `x` (lit-side normal `n`): the mean over the visible
# cone directions dⱼ of the best-focus profile where the ray along dⱼ through `x`
# crosses the focal plane, times the ray's incidence factor |n·dⱼ|/(k·dⱼ). Each
# direction's sub-beam is a bundle of parallel rays, so the visibility mask that
# selects it is exactly consistent with where its light goes, and power is conserved
# at shadow edges. (Giving each sub-beam a share of the f-cone spread would smooth the
# footprint but break that consistency: 0.6% of a beam clipped by an LEH edge is lost.)
# With `η` (the `absorption` of `deposit`), each direction's incidence factor carries
# its absorbed fraction η_a.
@inline function _cone_intensity(f::BeamFrame, ds, α, β, κ, m::UInt16, x, n, η = nothing)
    d = x - f.focus
    s = dot(d, f.k)
    u0 = dot(d, f.e1)
    v0 = dot(d, f.e2)
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

Absorbed laser intensity per element (TW/mm²; multiply by 1e14 for W/cm²).
"""
laser_flux(d::Deposition, t::Target) = d.element ./ t.elem_area

"""
    component_power(dep, target; by = :part) -> Dict{Symbol,Float64}

Absorbed laser power (TW) per part, plus back-side (outside) losses keyed as
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

Reference deposition by geometric ray tracing. Each spot component (core, and halo if
present) is sampled with an `nspot × nspot` grid of best-focus points spanning
`±extent` 1/e semi-axes (default: out to 1e-12 of the component's peak), each
launched along `ncone` f-cone directions and weighted by the focal intensity. Rays
travel straight (no refraction) to the first opaque surface. Power is conserved
exactly; defocus is exact within the sampling. Not differentiable.
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

Power (TW) each beam (rows) deposits on each component (columns, both sides), by the
reference ray trace — exact power accounting independent of mesh resolution.
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

# Each spot component (core, halo) is sampled on its own nspot × nspot grid spanning
# `extent` 1/e radii (default: out to where the component falls to 1e-12 of its peak).
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
