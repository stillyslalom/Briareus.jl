# This file builds the cluster tree of the hierarchical view factors and the visibility
# tests between clusters: shaft culling and the certificates of whole visibility.

# ------------------------------------------------------------------ cluster tree

struct _HCluster
    lo::Int                 # tree positions lo:hi
    hi::Int
    bmin::Vec3              # bounding box of the elements (centroids padded by size)
    bmax::Vec3
    cmin::Vec3              # bounding box of the centroids
    cmax::Vec3
    kids::NTuple{2,Int}     # (0, 0) for a leaf
    center::Int             # the element nearest the box center (the first pivot)
end

_hdiam(c::_HCluster) = norm(c.bmax - c.bmin)
_hdist(a::_HCluster, b::_HCluster) =
    norm(max.(zero(Vec3), max.(a.bmin - b.bmax, b.bmin - a.bmax)))
_hleaf(c::_HCluster) = c.kids == (0, 0)

# Bisects the centroids along the longest axis of their box, at the median. A cluster
# whose elements lie on several surfaces (`group`) is first split by surface, the largest
# from the rest. No cluster then mixes, for example, capsule and wall elements, or a cap
# and the wall at a corner. Such blocks seldom compress, and a capsule element in a wall
# cluster would bar the cluster from the ball certificate of the capsule.
function _hcluster_tree(C::Vector{Vec3}, pad::Vector{Float64}, leaf::Int, group::Vector{Int} = zeros(Int, length(C)))
    n = length(C)
    perm = collect(1:n)
    nodes = _HCluster[]
    function build(lo, hi)
        idx = view(perm, lo:hi)
        cmin = reduce((a, b) -> min.(a, b), (C[e] for e in idx))
        cmax = reduce((a, b) -> max.(a, b), (C[e] for e in idx))
        p = maximum(pad[e] for e in idx)
        mid = (cmin + cmax) / 2
        center = idx[argmin([norm(C[e] - mid) for e in idx])]
        push!(nodes, _HCluster(lo, hi, cmin .- p, cmax .+ p, cmin, cmax, (0, 0), center))
        me = length(nodes)
        g1 = group[idx[1]]
        if hi > lo && any(e -> group[e] != g1, idx)
            counts = Dict{Int,Int}()
            for e in idx
                counts[group[e]] = get(counts, group[e], 0) + 1
            end
            gmax = minimum(g for (g, k) in counts if k == maximum(values(counts)))
            sort!(idx; by = e -> (group[e] != gmax, e))
            m = lo + counts[gmax] - 1
        else
            hi - lo + 1 <= leaf && return me
            ax = argmax(cmax - cmin)
            m = (lo + hi) ÷ 2
            sort!(idx; by = e -> (C[e][ax], e))
        end
        a = build(lo, m)
        b = build(m + 1, hi)
        c = nodes[me]
        nodes[me] = _HCluster(c.lo, c.hi, c.bmin, c.bmax, c.cmin, c.cmax, (a, b), c.center)
        return me
    end
    n > 0 && build(1, n)
    return nodes, perm
end

# Returns, per element, an index of the surface it lies on: the curved surface of its
# component, or the component if it has none. Coplanar caps and LEH discs therefore share
# one index.
function _hsurface_groups(t::Target)
    keys = Any[]
    gc = map(t.components) do c
        k = c.surface === nothing ? c : c.surface
        i = findfirst(x -> x == k, keys)
        i === nothing ? (push!(keys, k); length(keys)) : i
    end
    return [gc[c] for c in t.elem_comp]
end

# Returns the cluster pairs (s, t), s before t or s = t, that cover the upper triangle of
# G. They are geometric candidates for low rank (`far`, decided per block in assembly)
# and leaf pairs (`near`).
function _hblock_pairs(nodes, η)
    far, near = Tuple{Int,Int}[], Tuple{Int,Int}[]
    function rec(s, t)
        a, b = nodes[s], nodes[t]
        if s != t && _hdist(a, b) >= η * max(_hdiam(a), _hdiam(b))
            push!(far, (s, t))
        elseif _hleaf(a) && _hleaf(b)
            push!(near, (s, t))
        elseif s == t
            k1, k2 = a.kids
            rec(k1, k1); rec(k1, k2); rec(k2, k2)
        elseif _hleaf(b) || (!_hleaf(a) && _hdiam(a) >= _hdiam(b))
            rec(a.kids[1], t); rec(a.kids[2], t)
        else
            rec(s, b.kids[1]); rec(s, b.kids[2])
        end
    end
    isempty(nodes) || rec(1, 1)
    return far, near
end

# ------------------------------------------------------------------ assembly

# Visibility between two clusters. The shaft, the convex hull of their boxes, contains
# every segment between their elements. If no facet that could occlude one of those
# segments meets it, the pair is wholly visible wherever it faces. The test is exact and
# conservative (the shaft culling of Haines & Wallace). Facets of the elements of the two
# clusters are not occluders, since a ray leaves the facet of its element and ends on the
# receiving facet. The facets of a curved surface that both clusters lie inside are not
# occluders either (the rule of `_vf_visible`).

# Returns the planes (n, d), with n·x ≤ d inside, of the convex hull of two boxes. They
# are the faces of the union box and the planes through an edge of each box, parallel to
# an axis; every face of the hull is one of these (Haines & Wallace 1994). Each candidate
# is kept only if it supports both boxes, so the planes always bound a region that
# contains the hull.
function _shaft_planes(alo, ahi, blo, bhi)
    lo, hi = min.(alo, blo), max.(ahi, bhi)
    planes = Tuple{Vec3,Float64}[]
    for k in 1:3
        ek = Vec3(ntuple(i -> i == k ? 1.0 : 0.0, 3))
        push!(planes, (ek, hi[k]), (-ek, -lo[k]))
    end
    ca = _box_corners(alo, ahi)
    cb = _box_corners(blo, bhi)
    L = norm(hi - lo)
    tol = 1e-12 * L
    for k in 1:3
        ek = Vec3(ntuple(i -> i == k ? 1.0 : 0.0, 3))
        for pa in ca, pb in cb
            (pa[k] == alo[k] && pb[k] == blo[k]) || continue    # one edge per box, by its lower end
            nv = cross(ek, pb - pa)
            nn = norm(nv)
            nn <= 1e-9 * L && continue
            nv = nv / nn
            d = dot(nv, pa)
            hmax, hmin = -Inf, Inf
            for p in ca
                h = dot(nv, p) - d; hmax = max(hmax, h); hmin = min(hmin, h)
            end
            for p in cb
                h = dot(nv, p) - d; hmax = max(hmax, h); hmin = min(hmin, h)
            end
            if hmax <= tol
                push!(planes, (nv, d + tol))
            elseif hmin >= -tol
                push!(planes, (-nv, -d + tol))
            end
        end
    end
    return planes
end

_box_corners(lo, hi) = (Vec3(lo[1], lo[2], lo[3]), Vec3(hi[1], lo[2], lo[3]), Vec3(lo[1], hi[2], lo[3]),
                        Vec3(hi[1], hi[2], lo[3]), Vec3(lo[1], lo[2], hi[3]), Vec3(hi[1], lo[2], hi[3]),
                        Vec3(lo[1], hi[2], hi[3]), Vec3(hi[1], hi[2], hi[3]))

# Returns what occludes the shaft between clusters a and b, as an `_HOcc`. The result is
# one of:
# - nothing that could occlude a segment between them;
# - only facets of up to 8 components with a certificate. The certificates are the balls
#   of a convex part, a planar part, and the axisymmetric bounds of a cylinder or cone.
#   `_hcert_verdict` decides each segment or leaves it to the ray;
# - a facet of another component.
#
# The tests apply to the pairs that `_vf_pair_free` evaluates without a ray, which face
# beyond the slack and so leave their facets. With `:sampled` visibility they apply to
# every ray between the samples of the pair lifted by ε. The plane tests below then use
# the box of those points, `opts.qbox`, and the least height of the samples of a cluster
# off a plane, `cplane`. A facet is not an occluder if:
# - it lies on a curved surface (cylinder, sphere, cone) that both clusters lie inside;
# - it belongs to a convex component (every vertex behind every facet) on which one of
#   the clusters lies entirely, with the facets of each element in one plane through its
#   centroid with its normal. A segment that leaves a face of a convex polyhedron outward
#   stays outside it;
# - each cluster lies strictly in front of its plane, or is planar in that plane facing
#   the same way. The lifted ray ends are then in front, so the segment does not cross
#   the plane.
# For the pairs admitted by the slack (`_hslack_verdict`, whose rays start farther off
# their facets), the result also records whether every plane exclusion holds with that
# margin (`liftx`), and which convex components on which a cluster lies were met in the
# shaft. With `full`, up to 8 certified components are allowed. Otherwise only the balls
# of a convex part with both boxes outside its circumscribed ball are allowed, and the
# search stops at the first other facet, as for the tests of an empty shaft.
function _shaft_occluder(ctx, opts, sa::Int, sb::Int, full::Bool = false)
    (; t, bvh, tri_surf, ε) = ctx
    a, b = opts.nodes[sa], opts.nodes[sb]
    allin, cplane, ccomp, convex = opts.allin, opts.cplane, opts.ccomp, opts.convex
    planes = _shaft_planes(a.bmin, a.bmax, b.bmin, b.bmax)
    outside_box(lo, hi) = any(((n, d),) -> dot(n, ifelse.(n .> 0, lo, hi)) > d, planes)
    # (boxes of the ray ends, which are centroids, or with `:sampled` visibility the
    # lifted samples)
    ca, cb = _box_corners(opts.qbox[sa]...), _box_corners(opts.qbox[sb]...)
    L = norm(max.(a.cmax, b.cmax) - min.(a.cmin, b.cmin))
    # 0 for not clear of the plane, 1 for clear, 2 for clear also with the slack lifts
    function clear(c, corners, n, d)
        h = minimum(x -> dot(n, x) - d, corners)
        h > 2ε + 1e-12 * L && return h > 2ε + 1e-12 * L + opts.liftx[c] ? 2 : 1
        for (pn, pd, off) in cplane[c]
            dot(pn, n) > 1 - 1e-12 && abs(pd - d) <= 1e-12 * (L + abs(d)) && off > 2ε + 1e-12 * L &&
                return off > 2ε + 1e-12 * L + opts.liftx[c] ? 2 : 1
        end
        return 0
    end
    # Components whose facets need no test are the convex ones on which a cluster lies,
    # which are never occluders, and the certified ones once one of their facets is in
    # the shaft.
    excluded(c) = c > 0 && convex[c] && (ccomp[sa] == c || ccomp[sb] == c)
    isempty(bvh.nodes) && return _HOCC_EMPTY
    found = ntuple(_ -> 0, 8)
    nf = 0
    slackok = true
    xa = xb = 0
    certifies(c) = full ? opts.cert[c] != 0 || opts.ball[c] != 0 : _hball_whole(opts, c, a, b)
    stack = Int32[1]
    while !isempty(stack)
        ni = pop!(stack)
        node = bvh.nodes[ni]
        sig = opts.bvh_comp[ni]
        sig > 0 && sig in found && continue
        outside_box(node.lo, node.hi) && continue
        if excluded(sig)
            sig == ccomp[sa] ? (xa = sig) : (xb = sig)
            continue
        end
        if node.b < 0
            for k in node.a:(node.a - node.b - 1)
                it = bvh.items[k]
                v0 = bvh.v0[it]
                v1, v2 = v0 + bvh.e1[it], v0 + bvh.e2[it]
                any(((n, d),) -> dot(n, v0) > d && dot(n, v1) > d && dot(n, v2) > d, planes) && continue
                id = bvh.ids[it]
                sf = tri_surf[id]
                sf != 0 && allin[sf, sa] && allin[sf, sb] && continue
                c = t.tri_comp[id]
                if excluded(c)
                    c == ccomp[sa] ? (xa = Int(c)) : (xb = Int(c))
                    continue
                end
                n = t.tri_normal[id]
                d = dot(n, v0)
                ka = clear(sa, ca, n, d)
                kb = ka == 0 ? 0 : clear(sb, cb, n, d)
                if kb != 0
                    slackok &= ka == 2 && kb == 2
                    continue
                end
                c in found && continue
                nf < (full ? 8 : 1) && certifies(c) || return _HOCC_RAYS
                nf += 1
                found = Base.setindex(found, Int(c), nf)
            end
        else
            push!(stack, node.a, node.b)
        end
    end
    return _HOcc(nf, found, slackok, (xa, xb))
end

# Returns how the entries of a block are computed from its shaft: the empty shaft itself
# (no rays), or else `false` (a ray per pair).
_hfree(occ) = occ.n == 0 ? occ : false

# What occludes a shaft (`_shaft_occluder`). `n = 0` means nothing, and `n = -1` a facet
# that no certificate decides, so that every pair takes a ray. Otherwise `c[1:n]` are the
# components whose certificates decide most pairs (`_hvis_cert`): the balls of a convex
# part, the facets of a planar part (`_HPlanar`), or the axisymmetric bounds of a part
# (`_HHull`). For the pairs that the facing slack admits (`_hslack_verdict`), `slack`
# records whether the plane exclusions hold with the lifts of their rays, and `x` the
# convex components excluded because a cluster lies on them.
struct _HOcc
    n::Int
    c::NTuple{8,Int}
    slack::Bool
    x::NTuple{2,Int}
end
const _HOCC_EMPTY = _HOcc(0, ntuple(_ -> 0, 8), true, (0, 0))
const _HOCC_RAYS = _HOcc(-1, ntuple(_ -> 0, 8), false, (0, 0))

# Whether the balls of convex component c decide every pair between clusters a and b,
# except those between the balls. This holds when c has balls (`occluders`) and both
# boxes lie outside its circumscribed ball. It is the only certificate of
# `_shaft_occluder` without `full`.
function _hball_whole(opts, c, a::_HCluster, b::_HCluster)
    k = opts.ball[c]
    k == 0 && return false
    _, ctr, _, rout = opts.occluders[k]
    return _box_distance(ctr, a.bmin, a.bmax) > rout && _box_distance(ctr, b.bmin, b.bmax) > rout
end

# A planar component, such as a washer, a disc or a foil. Every facet lies within `tp` of
# the plane n·(x - c) = 0, with in-plane basis (u, v) and in-plane distances from c within
# [rlo, rhi]. Its facets as seen in the plane (inward edge lines, radial range) are binned
# on a polar grid (nθ × nρ from ρ0 by dρ, CSR `start`/`idx`). Each facet is placed in the
# bins that its angular span and radial range overlap and in their neighbors.
struct _HPlanar
    c::Vec3
    n::Vec3
    u::Vec3
    v::Vec3
    tp::Float64
    rlo::Float64
    rhi::Float64
    nθ::Int
    nρ::Int
    ρ0::Float64
    dρ::Float64
    start::Vector{Int32}
    idx::Vector{Int32}
    edges::Vector{NTuple{9,Float64}}    # (nx, ny, d) × 3, inside where nx X + ny Y > d
    rr::Vector{NTuple{2,Float64}}       # least and largest in-plane distance from c
end

# The axisymmetric bounds of a component about the axis (c, a) of its surface, with every
# facet within a·(x - c) ∈ [z0, z1] and at distance [rlo, rhi] from the axis. The facets
# may all be parallel to the axis and lie, seen along it, on one of `ns` > 0 lines, as in
# a prism such as a wall of chord rectangles. The bounds then also hold those side lines
# in the basis (u, v) of the plane normal to the axis, as (νx, νy, p), the outward unit
# normal and the distance from the axis. The axis lies inside every line (p ≥ pmin > 0),
# and each facet lies within `pf` of its line. Lines are binned by the azimuth of ν (`nb`
# bins, CSR `start`/`idx`). A bin lists every line whose ν is within Δ of the bin, with
# cos Δ = `cosΔ`. A point at distance r from the axis with r cos Δ ≤ pmin - M therefore
# lies at least M inside every line that its bin leaves out (`_hprism_inside`).
struct _HHull
    c::Vec3
    a::Vec3
    z0::Float64
    z1::Float64
    rlo::Float64
    rhi::Float64
    u::Vec3
    v::Vec3
    ns::Int
    pmin::Float64
    pf::Float64
    cosΔ::Float64
    nb::Int
    start::Vector{Int32}
    idx::Vector{Int32}
    sides::Vector{NTuple{3,Float64}}
end

# Returns the distance from the origin to the triangle (A, B, C) in the plane.
function _origin_tri_distance(A, B, C)
    function seg(P, Q)
        e = Q - P
        ee = dot(e, e)
        s = ee > 0 ? clamp(-dot(P, e) / ee, 0.0, 1.0) : 0.0
        return norm(P + s * e)
    end
    d = min(seg(A, B), seg(B, C), seg(C, A))
    cr(P, Q) = P[1] * Q[2] - P[2] * Q[1]
    s1, s2, s3 = cr(A, B), cr(B, C), cr(C, A)
    inside = (s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0)
    return inside ? 0.0 : d
end

# Returns the certificate of component facets `ks` if they are planar, or else nothing.
function _hplanar(t::Target, ks, surface)
    n = t.tri_normal[ks[1]]
    all(k -> abs(dot(t.tri_normal[k], n)) > 1 - 1e-12, ks) || return nothing
    tri(k) = (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
    c = surface isa DiskSurface ? surface.center : sum(v for k in ks for v in tri(k)) / (3 * length(ks))
    R = maximum(norm(v - c) for k in ks for v in tri(k))
    tp = maximum(abs(dot(n, v - c)) for k in ks for v in tri(k))
    tp <= 1e-10 * (R + norm(c)) || return nothing
    u, v = _ortho_basis(n)
    P(x) = SVector(dot(u, x - c), dot(v, x - c))
    nf = length(ks)
    edges = Vector{NTuple{9,Float64}}(undef, nf)
    rr = Vector{NTuple{2,Float64}}(undef, nf)
    span = Vector{NTuple{2,Float64}}(undef, nf)        # angular span (lo, hi), or (-Inf, Inf) for all
    for (f, k) in enumerate(ks)
        A, B, C = map(P, tri(k))
        rr[f] = (_origin_tri_distance(A, B, C), max(norm(A), norm(B), norm(C)))
        σ = sign((B - A)[1] * (C - A)[2] - (B - A)[2] * (C - A)[1])
        line(p, q) = (e = q - p; l = norm(e); nx = -σ * e[2] / l; ny = σ * e[1] / l; (nx, ny, nx * p[1] + ny * p[2]))
        # (A degenerate facet gets NaN lines. It is never decided, so its pairs take the ray.)
        edges[f] = σ == 0 ? ntuple(_ -> NaN, 9) : (line(A, B)..., line(B, C)..., line(C, A)...)
        if rr[f][1] > 0
            θ1 = atan(A[2], A[1])
            δ(p) = rem2pi(atan(p[2], p[1]) - θ1, RoundNearest)
            lo, hi = θ1 + min(0.0, δ(B), δ(C)), θ1 + max(0.0, δ(B), δ(C))
            span[f] = hi - lo < π ? (lo, hi) : (-Inf, Inf)
        else
            span[f] = (-Inf, Inf)
        end
    end
    rlo, rhi = minimum(first, rr), maximum(last, rr)
    wθ = [s[2] - s[1] for s in span if isfinite(s[1])]
    wρ = [r[2] - r[1] for r in rr]
    nθ = isempty(wθ) ? 8 : clamp(round(Int, 2π / max(sort!(wθ)[cld(end, 2)], 1e-6)), 8, 4096)
    dr = sort!(wρ)[cld(end, 2)]
    nρ = clamp(dr > 0 ? round(Int, (rhi - rlo) / dr) : 1, 1, max(1, min(256, 2^20 ÷ nθ)))
    dρ = max(rhi - rlo, 1e-12 * (R + 1)) / nρ
    bins = [Int32[] for _ in 1:(nθ * nρ)]
    for f in 1:nf
        lo, hi = span[f]
        bθs = isfinite(lo) ? ((floor(Int, (lo + π) / 2π * nθ) - 1):(floor(Int, (hi + π) / 2π * nθ) + 1)) : (0:(nθ - 1))
        length(bθs) >= nθ && (bθs = 0:(nθ - 1))
        bρs = max(0, floor(Int, (rr[f][1] - rlo) / dρ) - 1):min(nρ - 1, floor(Int, (rr[f][2] - rlo) / dρ) + 1)
        for bρ in bρs, bθ in bθs
            push!(bins[bρ * nθ + mod(bθ, nθ) + 1], Int32(f))
        end
    end
    start = Int32[1; 1 .+ cumsum(Int32.(length.(bins)))]
    return _HPlanar(c, n, u, v, tp, rlo, rhi, nθ, nρ, rlo, dρ, start, reduce(vcat, bins), edges, rr)
end

# Returns the axisymmetric bounds of component facets `ks` about the axis (c, a).
function _hhull(t::Target, ks, c, a)
    u, v = _ortho_basis(a)
    z0, z1, rlo, rhi = Inf, -Inf, Inf, 0.0
    for k in ks
        ps = (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
        for p in ps
            z = dot(a, p - c)
            z0, z1 = min(z0, z), max(z1, z)
        end
        A, B, C = map(p -> SVector(dot(u, p - c), dot(v, p - c)), ps)
        rlo = min(rlo, _origin_tri_distance(A, B, C))
        rhi = max(rhi, norm(A), norm(B), norm(C))
    end
    return _HHull(c, a, z0, z1, rlo, rhi, u, v, _hprism(t, ks, c, a, u, v)...)
end

# Returns the side lines of facets `ks` (the ns, pmin, pf, cosΔ, nb, start, idx and sides
# of `_HHull`). They are found if every facet is parallel to the axis (c, a) and, seen
# along it, lies within 1e-10 (relative) of a line that keeps the axis inside, with pmin
# > 1e-6 of the largest radius. Otherwise ns = 0. Facets are grouped into lines by the
# azimuth of their outward normals, in runs with one normal to 1e-12 and one offset to
# 1e-10. A line split into several runs is kept several times, which is harmless. Each
# line passes through the mid-range of the vertices of its facets, and `pf` is the
# largest half-range.
function _hprism(t::Target, ks, c, a, u, v)
    none = (0, 0.0, 0.0, -1.0, 0, Int32[1], Int32[], NTuple{3,Float64}[])
    P(x) = SVector(dot(u, x - c), dot(v, x - c))
    facets = Tuple{Float64,SVector{2,Float64},NTuple{3,SVector{2,Float64}}}[]
    rmax = 0.0
    for k in ks
        n = t.tri_normal[k]
        abs(dot(n, a)) <= 1e-9 || return none
        ν = SVector(dot(u, n), dot(v, n))
        ν = ν / norm(ν)
        vs = (P(t.tri_v0[k]), P(t.tri_v0[k] + t.tri_e1[k]), P(t.tri_v0[k] + t.tri_e2[k]))
        sum(x -> dot(ν, x), vs) < 0 && (ν = -ν)
        rmax = max(rmax, maximum(norm, vs))
        push!(facets, (atan(ν[2], ν[1]), ν, vs))
    end
    isempty(facets) && return none
    sort!(facets; by = first)
    tol = 1e-10 * (rmax + norm(c))
    sides = NTuple{3,Float64}[]
    pf = 0.0
    ν0, lo, hi = facets[1][2], Inf, -Inf
    for (k, (_, ν, vs)) in enumerate(facets)
        if k > 1 && !(dot(ν, ν0) > 1 - 1e-12 && abs(sum(x -> dot(ν0, x), vs) / 3 - (lo + hi) / 2) <= tol)
            push!(sides, (ν0[1], ν0[2], (lo + hi) / 2))
            pf = max(pf, (hi - lo) / 2)
            ν0, lo, hi = ν, Inf, -Inf
        end
        for x in vs
            h = dot(ν0, x)
            lo, hi = min(lo, h), max(hi, h)
        end
    end
    push!(sides, (ν0[1], ν0[2], (lo + hi) / 2))
    pf = max(pf, (hi - lo) / 2)
    pmin = minimum(last, sides)
    ns = length(sides)
    (pf <= tol && pmin > 1e-6 * rmax && ns <= 100_000) || return none
    cosΔ = 0.999 * pmin / rmax
    Δ = acos(cosΔ) + 1e-4
    nb = clamp(ns, 8, 4096)
    β = 2π / nb
    bins = [Int32[] for _ in 1:nb]
    for (k, (νx, νy, _)) in enumerate(sides)
        α = atan(νy, νx)
        bs = (floor(Int, (α - Δ + π) / β) - 1):(floor(Int, (α + Δ + π) / β) + 1)
        length(bs) >= nb && (bs = 0:(nb - 1))
        for b in bs
            push!(bins[mod(b, nb) + 1], Int32(k))
        end
    end
    start = Int32[1; 1 .+ cumsum(Int32.(length.(bins)))]
    return ns, pmin, pf, cosΔ, nb, start, reduce(vcat, bins), sides
end

# Whether the points xa and xb (from the axis, normal to it) lie at least M inside every
# side line of the prism H, which is M beyond the facets of each line seen along the axis.
function _hprism_inside(H::_HHull, xa, xb, M)
    M += 2 * H.pf
    pm = H.pmin - M
    pm > 0 || return false
    Xa, Ya, Xb, Yb = dot(H.u, xa), dot(H.v, xa), dot(H.u, xb), dot(H.v, xb)
    ra, rb = Xa^2 + Ya^2, Xb^2 + Yb^2
    c2 = H.cosΔ * abs(H.cosΔ)
    (ra * c2 <= pm^2 && rb * c2 <= pm^2) || return false
    # (a point inside the inscribed circle of the lines by M lies inside every line)
    return (ra <= pm^2 || _hprism_sides(H, Xa, Ya, M)) && (rb <= pm^2 || _hprism_sides(H, Xb, Yb, M))
end

# Whether (X, Y) lies at least M inside the lines of its bin.
@inline function _hprism_sides(H::_HHull, X, Y, M)
    b = clamp(floor(Int, (_atan_approx(Y, X) + π) / (2π) * H.nb), 0, H.nb - 1) + 1
    @inbounds for q in H.start[b]:(H.start[b+1] - 1)
        νx, νy, p = H.sides[H.idx[q]]
        p - (νx * X + νy * Y) >= M || return false
    end
    return true
end

# Returns atan(y, x) within 2e-5 rad (Abramowitz & Stegun 4.4.49 on the octant, 1e-5),
# for the bin of `_hprism_sides`, whose lists have 1e-4 of slack.
function _atan_approx(y, x)
    ax, ay = abs(x), abs(y)
    a = min(ax, ay) / max(ax, ay, floatmin(Float64))
    s = a * a
    r = a * (0.9998660 + s * (-0.3302995 + s * (0.1801410 + s * (-0.0851330 + s * 0.0208351))))
    r = ifelse(ay > ax, π / 2 - r, r)
    r = ifelse(x < 0, π - r, r)
    return ifelse(y < 0, -r, r)
end

# Returns the certificate of each component (> 0 into `planar`, < 0 into `hull`, 0 for
# none). Convex parts with balls keep their balls, and apertures are not occluders.
function _hrim_data(t::Target, ball)
    planar, hull = _HPlanar[], _HHull[]
    cert = zeros(Int, length(t.components))
    tris = [Int[] for _ in t.components]
    for (k, c) in enumerate(t.tri_comp)
        push!(tris[c], k)
    end
    for (c, comp) in enumerate(t.components)
        (comp.kind === :aperture || ball[c] != 0 || isempty(tris[c])) && continue
        s = comp.surface
        p = _hplanar(t, tris[c], s)
        if p !== nothing
            push!(planar, p)
            cert[c] = length(planar)
        elseif s isa CylinderSurface || s isa ConeSurface
            push!(hull, _hhull(t, tris[c], s isa CylinderSurface ? s.origin : s.apex, s.axis))
            cert[c] = -length(hull)
        end
    end
    return cert, planar, hull
end

# Returns the ball verdict of the segment from oi to oj for a convex part with inscribed
# ball (ctr, rin) and circumscribed ball rout. The verdict is 1 if the segment misses the
# circumscribed ball, and so every facet. It is 0 if the segment passes through the
# inscribed ball with both ends outside the circumscribed one, so that it crosses the
# closed surface, and -1 otherwise. Margins of 1e-9 relative are added to those in
# `rin`, `rout`.
#
# Bundles apply to every certificate. With `ri`, `rj` > 0, the verdict holds for every
# segment from a point within `ri` of oi to one within `rj` of oj. With `:sampled`
# visibility these are the segments between the ray samples of two elements. The point
# of such a segment at parameter s lies within ρ = max(ri, rj) of the point of the
# segment from oi to oj at s, in a tube about it. The heights of its ends over a plane
# lie within ri, rj of those of oi, oj. Zero radii give the verdict of the segment,
# identical to that without them.
function _ball_verdict(oi, oj, ctr, rin, rout, ri = 0.0, rj = 0.0)
    ro = rout * (1 + 1e-9)
    (dot(oi - ctr, oi - ctr) > (ro + ri)^2 && dot(oj - ctr, oj - ctr) > (ro + rj)^2) || return Int8(-1)
    ρ = max(ri, rj)
    w = oj - oi
    u = clamp(dot(ctr - oi, w) / dot(w, w), 0.0, 1.0)
    x = oi + u * w - ctr
    q2 = dot(x, x)
    q2 > (ro + ρ)^2 && return Int8(1)
    rr = rin * (1 - 1e-9) - ρ
    rr > 0 && q2 < rr^2 && return Int8(0)
    return Int8(-1)
end

# Returns the planar verdict of the segment from oi to oj, as the ray of `_vf_visible`
# finds the facets (`any_hit` to 1 - 1e-9, barycentric slack 1e-10). The verdict is 1 if
# the segment stays on one side of the plane or crosses it outside every facet, 0 if it
# crosses inside a facet with its end and start off the plane, and -1 otherwise. Margins
# of 1e-9 relative and the flatness `tp` of the plane are grown by 1/sin of the crossing
# angle in the plane. The crossings of a bundle (radii ri, rj) lie within
# ρ (1 + |w|/|dh|) of the crossing of the segment from oi to oj, with ρ = max(ri, rj),
# since a point at height ≤ ρ on the axis of the tube is that far along it from the
# crossing.
function _planar_verdict(P::_HPlanar, oi, oj, ri = 0.0, rj = 0.0)
    hi, hj = dot(P.n, oi - P.c), dot(P.n, oj - P.c)
    L = max(norm(oi - P.c) + ri, norm(oj - P.c) + rj, P.rhi)
    m = 1e-9 * L + 2 * P.tp
    ((hi - ri > m && hj - rj > m) || (hi + ri < -m && hj + rj < -m)) && return Int8(1)
    ((hi - ri > m && hj + rj < -m) || (hi + ri < -m && hj - rj > m)) || return Int8(-1)
    w = oj - oi
    dh = hi - hj
    x = oi + (hi / dh) * w - P.c
    X, Y = dot(P.u, x), dot(P.v, x)
    ρ = sqrt(X^2 + Y^2)
    nw, adh = norm(w), abs(dh)
    M = m * (1 + (nw + ri + rj) / (adh - ri - rj)) + max(ri, rj) * (1 + nw / adh)
    (ρ > P.rhi + M || ρ < P.rlo - M) && return Int8(1)
    dθ = 2π / P.nθ
    2M < min(P.dρ, ρ * dθ) || return Int8(-1)
    bθ = clamp(floor(Int, (atan(Y, X) + π) / dθ), 0, P.nθ - 1)
    bρ = clamp(floor(Int, (ρ - P.ρ0) / P.dρ), 0, P.nρ - 1)
    b = bρ * P.nθ + bθ + 1
    # the crossing at a ray parameter well inside (0, 1 - 1e-9)
    ends = min(abs(hi) - ri, abs(hj) - rj) > 4e-9 * (adh + ri + rj) + m
    out = true
    @inbounds for q in P.start[b]:(P.start[b+1] - 1)
        f = P.idx[q]
        r0, r1 = P.rr[f]
        (ρ > r1 + M || ρ < r0 - M) && continue
        e = P.edges[f]
        s = min(e[1] * X + e[2] * Y - e[3], e[4] * X + e[5] * Y - e[6], e[7] * X + e[8] * Y - e[9])
        s < -M && continue
        s > M && ends && return Int8(0)
        out = false
    end
    return out ? Int8(1) : Int8(-1)
end

# Returns the axisymmetric verdict of the segment from oi to oj. The verdict is 1 if the
# part of the segment within the axial range of the facets, widened by the margin, lies
# inside rlo or outside rhi. On a prism it is also 1 if that part lies inside every side
# line by the margin, as for a segment through a waist gap onto the wall, which passes
# the edge of the wall outside its axial range. Otherwise the verdict is -1. A bundle
# (tube radius ρ) is decided as its axis against bounds widened by ρ.
function _hull_verdict(H::_HHull, oi, oj, ρ = 0.0)
    pi, pj = oi - H.c, oj - H.c
    zi, zj = dot(H.a, pi), dot(H.a, pj)
    L = max(norm(pi) + ρ, norm(pj) + ρ, H.rhi, abs(H.z0), abs(H.z1))
    m = 1e-9 * L
    lo, hi = H.z0 - m - ρ, H.z1 + m + ρ
    dz = zj - zi
    if abs(dz) <= 1e-6 * L
        (max(zi, zj) < lo || min(zi, zj) > hi) && return Int8(1)
        sa, sb = 0.0, 1.0
    else
        s1, s2 = (lo - zi) / dz, (hi - zi) / dz
        sa, sb = max(0.0, min(s1, s2) - 1e-9), min(1.0, max(s1, s2) + 1e-9)
        sa > sb && return Int8(1)
    end
    Ri = pi - zi * H.a
    Rw = (pj - zj * H.a) - Ri
    xa, xb = Ri + sa * Rw, Ri + sb * Rw
    rlo = H.rlo - ρ
    rlo > m && max(dot(xa, xa), dot(xb, xb)) < (rlo - m)^2 && return Int8(1)
    # on a prism, that part lies inside every side line (it is convex), so beyond its facets
    H.ns > 0 && _hprism_inside(H, xa, xb, m + ρ) && return Int8(1)
    ww = dot(Rw, Rw)
    s = ww > 0 ? clamp(-dot(Ri, Rw) / ww, sa, sb) : sa
    x = Ri + s * Rw
    dot(x, x) > (H.rhi + ρ + m)^2 && return Int8(1)
    return Int8(-1)
end

# Returns the verdict of the certificate of component c on the segment from oi to oj
# between elements i and j. A convex part is also cleared by the plane of the element of
# an end that lies on it (`onconvex`, with its facets in that plane). The part lies behind
# that plane, within the convexity tolerance, so a segment with both ends in front misses
# it. With radii ri, rj, the verdict is that of the bundle of segments between points
# within ri of oi and rj of oj, which are the ray samples of elements i and j
# (`_hvis_bundle`). The samples of an end element lie in its plane within `opts.hdev`.
function _hcert_verdict(ctx, opts, c, i, j, oi, oj, ri = 0.0, rj = 0.0)
    k = opts.ball[c]
    if k != 0
        _, ctr, rin, rout = opts.occluders[k]
        (; C, Nrm, hsep) = ctx
        opts.onconvex[i] == c && dot(Nrm[i], oi - C[i]) - min(ri, opts.hdev[i]) > hsep &&
            dot(Nrm[i], oj - C[i]) - rj > hsep && return Int8(1)
        opts.onconvex[j] == c && dot(Nrm[j], oi - C[j]) - ri > hsep &&
            dot(Nrm[j], oj - C[j]) - min(rj, opts.hdev[j]) > hsep && return Int8(1)
        return _ball_verdict(oi, oj, ctr, rin, rout, ri, rj)
    end
    k = opts.cert[c]
    return k > 0 ? _planar_verdict(opts.planar[k], oi, oj, ri, rj) : _hull_verdict(opts.hull[-k], oi, oj, max(ri, rj))
end

# Returns, per BVH node, the component of every facet below it, or -1 if the facets
# belong to several components.
function _bvh_components(bvh, tri_comp)
    sig = zeros(Int, length(bvh.nodes))
    function rec(ni)
        node = bvh.nodes[ni]
        if node.b < 0
            cs = unique(Int(tri_comp[bvh.ids[bvh.items[k]]]) for k in node.a:(node.a - node.b - 1))
            sig[ni] = length(cs) == 1 ? cs[1] : -1
        else
            x, y = rec(node.a), rec(node.b)
            sig[ni] = x == y ? x : -1
        end
        return sig[ni]
    end
    isempty(bvh.nodes) || rec(1)
    return sig
end

# Whether every segment between the boxes a and b crosses the ball (c, r), exactly for
# boxes outside the ball. For a point q outside the ball, the points x past the ball as
# seen from q, where the segment qx crosses it, include the convex set K_q. K_q is the
# cone from q over the ball, beyond the polar plane of q, (x - c)·(q - c) < r². Suppose
# every corner of a lies in K_q for every corner q of b. Then each segment from a corner
# of b to any point of a crosses the ball. Each corner of b then lies past the ball from
# any point p of a, and so in K_p, since a point outside the ball past it from p lies
# beyond the polar plane of p. All of b then lies in K_p, since K_p is convex, so every
# segment from a to b crosses the ball.
function _ball_blocks(alo, ahi, blo, bhi, c, r)
    for q in _box_corners(blo, bhi)
        w = c - q
        D2 = dot(w, w)
        D2 > r^2 || return false
        cosmax = sqrt(1 - r^2 / D2)                 # cone half-angle
        for x in _box_corners(alo, ahi)
            v = x - q
            nv = norm(v)
            nv > 0 || return false
            dot(v, w) >= cosmax * nv * sqrt(D2) || return false
            dot(x - c, q - c) < r^2 || return false
        end
    end
    return true
end

# Whether a convex component (closed, outward normals) hides every pair of clusters a and
# b from each other. Its inscribed ball (about the vertex mean) must block every segment
# between their boxes, with both boxes outside its circumscribed ball, so that each
# segment enters the component through a facet. Neither cluster may lie on it.
function _convex_hides(opts, sa, sb)
    a, b = opts.nodes[sa], opts.nodes[sb]
    for (c, ctr, rin, rout) in opts.occluders
        (opts.ccomp[sa] == c || opts.ccomp[sb] == c) && continue
        _box_distance(ctr, a.bmin, a.bmax) > rout && _box_distance(ctr, b.bmin, b.bmax) > rout || continue
        _ball_blocks(a.bmin, a.bmax, b.bmin, b.bmax, ctr, rin) && return true
    end
    return false
end
# Whether every vertex `vs` lies behind (within tol) the plane of every facet `ks`.
function _convex_mesh(t, ks, vs, tol)
    ok = Threads.Atomic{Bool}(true)
    Threads.@threads for k in ks
        ok[] || continue
        n, v0 = t.tri_normal[k], t.tri_v0[k]
        h = dot(n, v0) + tol
        for v in vs
            if dot(n, v) > h
                ok[] = false
                break
            end
        end
    end
    return ok[]
end

_box_distance(p, lo, hi) = norm(max.(zero(Vec3), max.(lo - p, p - hi)))

# Returns, per cluster, its plane (n, d) if every element is a flat patch of facets in
# one plane with one normal (else a zero normal), and its component if every element is
# a one-facet element of one component (else 0). Per component, it returns whether the
# mesh is convex with outward normals, for sphere-surface components only.
function _hcluster_data(t::Target, nodes, perm, pts = nothing)
    cplane = Vector{Tuple{Vec3,Float64,Float64}}[]
    ccomp = Int[]
    # elements that are flat patches (flat facets or disk sectors) of facets in one plane
    flat = map(1:nelements(t)) do e
        s = t.components[t.elem_comp[e]].surface
        (t.elem_curved[e] == 0 || s isa DiskSurface) || return false
        n, C = t.elem_normal[e], t.elem_centroid[e]
        tol = 1e-12 * (norm(C) + 1)
        all(k -> dot(t.tri_normal[k], n) > 1 - 1e-12 && abs(dot(n, t.tri_v0[k] - C)) <= tol, t.elem_tris[e])
    end
    for nd in nodes
        es = perm[nd.lo:nd.hi]
        L = norm(nd.cmax - nd.cmin) + norm(nd.cmax) + norm(nd.cmin)
        # The planes of its flat elements, each with the least height above it of the
        # other elements of the cluster (Inf if none). The heights are those of their
        # centroids, or with `:sampled` visibility of the lifted ray samples `pts[e]`.
        planes = Tuple{Vec3,Float64,Float64}[]
        for e in es
            flat[e] || continue
            n = t.elem_normal[e]
            d = dot(n, t.elem_centroid[e])
            any(((pn, pd, _),) -> dot(pn, n) > 1 - 1e-12 && abs(pd - d) <= 1e-12 * L, planes) && continue
            off = Inf
            for f in es
                on = flat[f] && dot(t.elem_normal[f], n) > 1 - 1e-12 && abs(dot(n, t.elem_centroid[f]) - d) <= 1e-12 * L
                on && continue
                if pts === nothing
                    off = min(off, dot(n, t.elem_centroid[f]) - d)
                else
                    for x in pts[f]
                        off = min(off, dot(n, x) - d)
                    end
                end
            end
            push!(planes, (n, d, off))
            length(planes) > 4 && break
        end
        push!(cplane, planes)
        c = t.elem_comp[es[1]]
        push!(ccomp, all(e -> t.elem_comp[e] == c && length(t.elem_tris[e]) == 1, es) ? Int(c) : 0)
    end
    convex = falses(length(t.components))
    occluders = Tuple{Int,Vec3,Float64,Float64}[]
    for (c, comp) in enumerate(t.components)
        comp.surface isa SphereSurface || continue
        ks = findall(==(c), t.tri_comp)
        vs = unique!([v for k in ks for v in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])])
        length(ks) * length(vs) <= 2_000_000_000 || continue
        R = maximum(v -> norm(v - vs[1]), vs)
        convex[c] = _convex_mesh(t, ks, vs, 1e-9 * R)
        if convex[c] && length(ks) >= 4
            ctr = sum(vs) / length(vs)
            rin = minimum(k -> dot(t.tri_normal[k], t.tri_v0[k] - ctr), ks)
            rout = maximum(v -> norm(v - ctr), vs)
            rin > 0 && push!(occluders, (c, ctr, rin * (1 - 1e-9), rout * (1 + 1e-9)))
        end
    end
    return cplane, ccomp, convex, occluders
end
