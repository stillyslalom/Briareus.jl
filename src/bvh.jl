# Bounding-volume hierarchy over triangles for ray casting (closest hit and any hit).

struct BVHNode
    lo::Vec3
    hi::Vec3
    a::Int32   # leaf: first index into `items`; interior: left child
    b::Int32   # leaf: -(item count);           interior: right child
end

"""
    BVH(v0, e1, e2, ids)

Median-split BVH over triangles given by vertex `v0` and edges `e1 = v1 - v0`,
`e2 = v2 - v0`. `ids[k]` is the caller's identifier for triangle `k` (returned by the
hit queries).
"""
struct BVH
    nodes::Vector{BVHNode}
    items::Vector{Int32}       # triangle slots, leaf-contiguous
    v0::Vector{Vec3}
    e1::Vector{Vec3}
    e2::Vector{Vec3}
    ids::Vector{Int32}
end

const _LEAF_SIZE = 4

function BVH(v0::Vector{Vec3}, e1::Vector{Vec3}, e2::Vector{Vec3}, ids::Vector{Int32})
    n = length(v0)
    items = Int32.(collect(1:n))
    nodes = BVHNode[]
    if n == 0
        return BVH(nodes, items, v0, e1, e2, ids)
    end
    lo = [min.(v0[k], v0[k] + e1[k], v0[k] + e2[k]) for k in 1:n]
    hi = [max.(v0[k], v0[k] + e1[k], v0[k] + e2[k]) for k in 1:n]
    cen = [(lo[k] + hi[k]) / 2 for k in 1:n]
    _build!(nodes, items, lo, hi, cen, 1, n)
    return BVH(nodes, items, v0, e1, e2, ids)
end

function _build!(nodes, items, lo, hi, cen, first, last)
    blo = Vec3(Inf, Inf, Inf)
    bhi = Vec3(-Inf, -Inf, -Inf)
    clo = Vec3(Inf, Inf, Inf)
    chi = Vec3(-Inf, -Inf, -Inf)
    for k in first:last
        it = items[k]
        blo = min.(blo, lo[it]); bhi = max.(bhi, hi[it])
        clo = min.(clo, cen[it]); chi = max.(chi, cen[it])
    end
    idx = length(nodes) + 1
    count = last - first + 1
    if count <= _LEAF_SIZE || maximum(chi - clo) == 0
        push!(nodes, BVHNode(blo, bhi, Int32(first), Int32(-count)))
        return idx
    end
    push!(nodes, BVHNode(blo, bhi, 0, 0))     # placeholder
    axis = argmax(chi - clo)
    mid = (first + last) ÷ 2
    v = view(items, first:last)
    partialsort!(v, mid - first + 1; by = it -> cen[it][axis])
    l = _build!(nodes, items, lo, hi, cen, first, mid)
    r = _build!(nodes, items, lo, hi, cen, mid + 1, last)
    nodes[idx] = BVHNode(blo, bhi, Int32(l), Int32(r))
    return idx
end

@inline function _slab(node::BVHNode, o::Vec3, inv::Vec3, tmax)
    tn, tf = _slab_span(node, o, inv)
    return tn <= tf && tf >= 0 && tn <= tmax
end

# Parameter interval of the ray inside a node's box. Every operand is finite (`_safe_inv`),
# so plain comparisons select the min and max (Base's handle NaN and signed zeros, and
# cost twice as much in this, the innermost loop of every ray query).
@inline function _slab_span(node::BVHNode, o::Vec3, inv::Vec3)
    t1 = (node.lo - o) .* inv
    t2 = (node.hi - o) .* inv
    lo = map(_fmin, t1, t2)
    hi = map(_fmax, t1, t2)
    return _fmax(_fmax(lo[1], lo[2]), lo[3]), _fmin(_fmin(hi[1], hi[2]), hi[3])
end
@inline _fmin(a, b) = ifelse(a < b, a, b)
@inline _fmax(a, b) = ifelse(a > b, a, b)

# Zero direction components map to a huge finite slope so slab tests never form 0*Inf.
@inline _safe_inv(d::Vec3) = map(x -> iszero(x) ? 1e300 : 1 / x, d)

const _BARY_TOL = 1e-10

# Möller–Trumbore; returns (t, cosine sign) with t = Inf on a miss.
@inline function _intersect(o::Vec3, d::Vec3, v0::Vec3, e1::Vec3, e2::Vec3, tmax)
    p = cross(d, e2)
    det = dot(e1, p)
    det == 0 && return Inf
    inv = 1 / det
    s = o - v0
    u = dot(s, p) * inv
    (u < -_BARY_TOL || u > 1 + _BARY_TOL) && return Inf
    q = cross(s, e1)
    v = dot(d, q) * inv
    (v < -_BARY_TOL || u + v > 1 + _BARY_TOL) && return Inf
    t = dot(e2, q) * inv
    return (t > 0 && t < tmax) ? t : Inf
end

"""
    closest_hit(bvh, origin, dir, tmax = Inf) -> (t, id)

Nearest intersection along `origin + t*dir` for `0 < t < tmax`. Returns `id = 0` on a
miss. `t` is in units of `|dir|`.
"""
function closest_hit(b::BVH, o::Vec3, d::Vec3, tmax = Inf)
    isempty(b.nodes) && return (Inf, Int32(0))
    inv = _safe_inv(d)
    stack = MVector{64,Int32}(undef)
    sp = 1
    stack[1] = 1
    best = Float64(tmax)
    hit = Int32(0)
    @inbounds while sp > 0
        ni = stack[sp]; sp -= 1
        node = b.nodes[ni]
        _slab(node, o, inv, best) || continue
        if node.b < 0
            for k in node.a:(node.a - node.b - 1)
                it = b.items[k]
                t = _intersect(o, d, b.v0[it], b.e1[it], b.e2[it], best)
                if t < best
                    best = t
                    hit = b.ids[it]
                end
            end
        else
            stack[sp+1] = node.a
            stack[sp+2] = node.b
            sp += 2
        end
    end
    return (best, hit)
end

_skipped(::Nothing, id) = false
_skipped(f, id) = f(id)

"""
    any_hit(bvh, origin, dir, tmax = Inf, skip = nothing) -> Bool

Whether any triangle intersects `origin + t*dir` for `0 < t < tmax`, ignoring triangles
whose identifier `id` has `skip(id) == true`.
"""
function any_hit(b::BVH, o::Vec3, d::Vec3, tmax = Inf, skip::S = nothing) where {S}
    isempty(b.nodes) && return false
    inv = _safe_inv(d)
    stack = MVector{64,Int32}(undef)
    sp = 1
    stack[1] = 1
    @inbounds while sp > 0
        ni = stack[sp]; sp -= 1
        node = b.nodes[ni]
        _slab(node, o, inv, tmax) || continue
        if node.b < 0
            for k in node.a:(node.a - node.b - 1)
                it = b.items[k]
                _intersect(o, d, b.v0[it], b.e1[it], b.e2[it], tmax) < Inf &&
                    !_skipped(skip, b.ids[it]) && return true
            end
        else
            stack[sp+1] = node.a
            stack[sp+2] = node.b
            sp += 2
        end
    end
    return false
end

"""
    _unblocked!(open, bvh, origins, dirs, invs, cache) -> open

Packet form of `any_hit` with `tmax = Inf` for rays sharing a few directions: on entry
bit `j-1` of `open[i]` selects the ray `origins[i] + t*dirs[j]` (`length(dirs) ≤ 16`,
`invs[j] = _safe_inv(dirs[j])`); on exit only the bits of rays that hit nothing remain,
exactly as `!any_hit(bvh, origins[i], dirs[j])` answers. `cache[j]` is a leaf (node
index) that blocked direction `j` before, or 0; it is tested first and updated.

The packet is traversed once. Interior nodes are culled by a slab test conservative
over the origins' bounding box and the directions' inverse components, leaves by one
per direction; then each ray takes its own slab test and triangle tests. fp subtraction
and products are monotone, so the conservative tests never cull a node a ray passes,
and the slab test is monotone in box containment, so a ray passes a leaf's test only if
it passes every ancestor's: each ray tests only the triangles its own `any_hit`
traversal reaches, and each of them until one blocks it.
"""
function _unblocked!(open::AbstractVector{UInt16}, b::BVH, os::AbstractVector{Vec3},
                     dirs::AbstractVector{Vec3}, invs::AbstractVector{Vec3},
                     cache::AbstractVector{Int32})
    n, nd = length(os), length(dirs)
    nd <= 16 || throw(ArgumentError("at most 16 directions"))
    (isempty(b.nodes) || n == 0) && return open
    alive = 0
    used = UInt16(0)
    @inbounds for i in 1:n
        m = open[i]
        for j in 1:nd
            (m >> (j - 1)) & 0x1 == 0x1 || continue
            c = cache[j]
            c > 0 && _leaf_blocks(b, b.nodes[c], os[i], dirs[j], invs[j]) &&
                (m &= ~(UInt16(1) << (j - 1)))
        end
        open[i] = m
        used |= m
        alive += count_ones(m)
    end
    alive == 0 && return open
    ilo, ihi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
    @inbounds for j in 1:nd
        if (used >> (j - 1)) & 0x1 == 0x1
            ilo = map(_fmin, ilo, invs[j]); ihi = map(_fmax, ihi, invs[j])
        end
    end
    olo = ohi = os[1]
    @inbounds for i in 2:n
        olo = map(_fmin, olo, os[i]); ohi = map(_fmax, ohi, os[i])
    end
    stack = MVector{64,Int32}(undef)
    sp = 1
    stack[1] = 1
    @inbounds while sp > 0
        ni = stack[sp]; sp -= 1
        node = b.nodes[ni]
        _packet_slab(node, olo, ohi, ilo, ihi) || continue
        if node.b < 0
            live = UInt16(0)
            for j in 1:nd
                (used >> (j - 1)) & 0x1 == 0x1 || continue
                _packet_slab(node, olo, ohi, invs[j], invs[j]) && (live |= UInt16(1) << (j - 1))
            end
            live == 0 && continue
            for i in 1:n
                m = open[i]
                r = m & live
                r == 0 && continue
                for j in 1:nd
                    (r >> (j - 1)) & 0x1 == 0x1 || continue
                    if _leaf_blocks(b, node, os[i], dirs[j], invs[j])
                        m &= ~(UInt16(1) << (j - 1))
                        cache[j] = ni
                        alive -= 1
                    end
                end
                open[i] = m
            end
            alive == 0 && return open
        else
            stack[sp+1] = node.a
            stack[sp+2] = node.b
            sp += 2
        end
    end
    return open
end

# Whether the ray `o + t*d` (t > 0) passes leaf `node`'s slab test and hits one of its
# triangles: `any_hit`'s test at that leaf.
@inline function _leaf_blocks(b::BVH, node::BVHNode, o::Vec3, d::Vec3, inv::Vec3)
    _slab(node, o, inv, Inf) || return false
    @inbounds for k in node.a:(node.a - node.b - 1)
        it = b.items[k]
        _intersect(o, d, b.v0[it], b.e1[it], b.e2[it], Inf) < Inf && return true
    end
    return false
end

# Slab test passed whenever a ray with origin in [olo, ohi] and inverse direction in
# [ilo, ihi] passes `_slab(node, o, inv, Inf)`: per axis, the ray's plane crossings lie
# between the extreme corner products (an axis whose inverse components change sign
# does not cull).
@inline function _packet_slab(node::BVHNode, olo::Vec3, ohi::Vec3, ilo::Vec3, ihi::Vec3)
    tn, tf = -Inf, Inf
    @inbounds for a in 1:3
        p, q = ilo[a], ihi[a]
        (p > 0 || q < 0) || continue
        x1, x2 = node.lo[a] - ohi[a], node.lo[a] - olo[a]
        x3, x4 = node.hi[a] - ohi[a], node.hi[a] - olo[a]
        lo = _fmin(_fmin(_fmin(x1 * p, x1 * q), _fmin(x2 * p, x2 * q)),
                   _fmin(_fmin(x3 * p, x3 * q), _fmin(x4 * p, x4 * q)))
        hi = _fmax(_fmax(_fmax(x1 * p, x1 * q), _fmax(x2 * p, x2 * q)),
                   _fmax(_fmax(x3 * p, x3 * q), _fmax(x4 * p, x4 * q)))
        tn = _fmax(tn, lo)
        tf = _fmin(tf, hi)
    end
    return tn <= tf && tf >= 0
end

"""
    _sah_bvh(bvh) -> BVH

The same triangles as `bvh` regrouped by a binned surface-area heuristic, which
`Target` uses for every ray query. A median-split tree has fat interior boxes inside a
concave enclosure, so a chord across a hohlraum visits ~700 nodes and ~140 triangles;
this tree visits ~55 and ~10. Hit answers are unchanged: the triangle tests are the
same and only the culling boxes differ, which matters only for a ray grazing a
triangle within round-off.
"""
function _sah_bvh(b::BVH; leaf::Integer = 4, nbins::Integer = 16)
    n = length(b.v0)
    items = Int32.(collect(1:n))
    nodes = BVHNode[]
    n == 0 && return BVH(nodes, items, b.v0, b.e1, b.e2, b.ids)
    lo = [min.(b.v0[k], b.v0[k] + b.e1[k], b.v0[k] + b.e2[k]) for k in 1:n]
    hi = [max.(b.v0[k], b.v0[k] + b.e1[k], b.v0[k] + b.e2[k]) for k in 1:n]
    cen = [(lo[k] + hi[k]) / 2 for k in 1:n]
    bins = (zeros(Int, nbins), Vector{Vec3}(undef, nbins), Vector{Vec3}(undef, nbins),
            zeros(nbins))
    _sah_build!(nodes, items, lo, hi, cen, 1, n, leaf, bins)
    return BVH(nodes, items, b.v0, b.e1, b.e2, b.ids)
end

_box_area(lo, hi) = (e = hi - lo; 2 * (e[1] * e[2] + e[2] * e[3] + e[3] * e[1]))
_bin(c, ax, clo, chi, nb) = clamp(floor(Int, (c[ax] - clo[ax]) / (chi[ax] - clo[ax]) * nb) + 1, 1, nb)

function _sah_build!(nodes, items, lo, hi, cen, first, last, leaf, bins)
    blo, bhi = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf)
    clo, chi = blo, bhi
    for k in first:last
        it = items[k]
        blo = min.(blo, lo[it]); bhi = max.(bhi, hi[it])
        clo = min.(clo, cen[it]); chi = max.(chi, cen[it])
    end
    idx = length(nodes) + 1
    count = last - first + 1
    if count <= leaf || maximum(chi - clo) == 0
        push!(nodes, BVHNode(blo, bhi, Int32(first), Int32(-count)))
        return idx
    end
    push!(nodes, BVHNode(blo, bhi, 0, 0))     # placeholder
    # split between centroid bins at the least Σ count × box area over the two sides
    bc, bl, bh, rcost = bins
    nb = length(bc)
    best, bax, bs = Inf, 0, 0
    for ax in 1:3
        chi[ax] > clo[ax] || continue
        fill!(bc, 0); fill!(bl, Vec3(Inf, Inf, Inf)); fill!(bh, Vec3(-Inf, -Inf, -Inf))
        for k in first:last
            it = items[k]
            s = _bin(cen[it], ax, clo, chi, nb)
            bc[s] += 1; bl[s] = min.(bl[s], lo[it]); bh[s] = max.(bh[s], hi[it])
        end
        rl, rh, nr = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf), 0
        for s in nb:-1:2
            rl = min.(rl, bl[s]); rh = max.(rh, bh[s]); nr += bc[s]
            rcost[s] = nr == 0 ? Inf : nr * _box_area(rl, rh)
        end
        ll, lh, nl = Vec3(Inf, Inf, Inf), Vec3(-Inf, -Inf, -Inf), 0
        for s in 1:nb-1
            ll = min.(ll, bl[s]); lh = max.(lh, bh[s]); nl += bc[s]
            nl == 0 && continue
            c = nl * _box_area(ll, lh) + rcost[s+1]
            c < best && ((best, bax, bs) = (c, ax, s))
        end
    end
    v = view(items, first:last)
    if bax == 0
        mid = (first + last) ÷ 2
        axis = argmax(chi - clo)
        partialsort!(v, mid - first + 1; by = it -> cen[it][axis])
    else
        m = 0
        for k in eachindex(v)
            _bin(cen[v[k]], bax, clo, chi, nb) <= bs && (m += 1; (v[m], v[k]) = (v[k], v[m]))
        end
        mid = first + m - 1
    end
    l = _sah_build!(nodes, items, lo, hi, cen, first, mid, leaf, bins)
    r = _sah_build!(nodes, items, lo, hi, cen, mid + 1, last, leaf, bins)
    nodes[idx] = BVHNode(blo, bhi, Int32(l), Int32(r))
    return idx
end
