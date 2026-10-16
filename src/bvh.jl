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
    t1 = (node.lo - o) .* inv
    t2 = (node.hi - o) .* inv
    tn = maximum(min.(t1, t2))
    tf = minimum(max.(t1, t2))
    return tn <= tf && tf >= 0 && tn <= tmax
end

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
function any_hit(b::BVH, o::Vec3, d::Vec3, tmax = Inf, skip = nothing)
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
