# Hierarchical (H-matrix) view factors: near pairs dense, far cluster pairs of uniform
# visibility low rank (adaptive cross approximation on the exact pair kernel), so
# storage and assembly grow as N^1.3–1.5 rather than N² on a hohlraum (the capsule's
# shadow edges and the near field stay dense). Only products are offered (`mul!` with F and Fᵀ), which is
# all the iterative radiosity path uses.

"""
    HMatrix

View factors in hierarchical form (`view_factors(target; storage = :hmatrix)`, and
`RadiationModel` by default above `Briareus.HMATRIX_MODEL[]` elements): the symmetric
exchange matrix `G = diag(A) F` (`Gᵢⱼ = Aᵢ Fᵢⱼ`, stored once per unordered pair, so
reciprocity is exact) partitioned over a cluster tree of the elements, and `F = diag(A)⁻¹
S G S` with a diagonal balancing `S` (`enforce_closure!` updates `S` only). Reads as an
`AbstractMatrix{Float64}`; products with Float64 vectors and matrices (`*`, `mul!`, also
transposed, and with dual vectors) accumulate in Float64, which is all the iterative
radiosity path uses. Indexing one entry searches the blocks (for inspection only);
`Matrix(H)` expands it by products with the identity.

Clusters split the elements by the surface they lie on (a component's curved surface, or
the component), then bisect the centroids along the longest axis of their box down to
`leaf` elements, so that no cluster mixes capsule and wall, or cap and wall at a corner
(such blocks seldom compress, and their low-rank error falls on the small rows). A pair of clusters whose boxes are `eta` × the larger diameter apart is
a candidate for low rank, decided block by block:
- rows and columns whose elements face no centroid of the other cluster are zero (a box
  test), and so is a block that a convex closed part (a capsule) hides wholly (its
  inscribed ball blocks every segment between the boxes, a cone test from each corner);
- the *shaft* between the boxes (their convex hull) is searched for facets that could
  occlude a segment between them (Haines & Wallace's shaft culling, with the facets of
  a curved surface both lie inside, of a convex part one of them lies on, and of a plane
  both lie in front of or in excluded). An empty shaft makes every facing pair wholly
  visible: its entries need no ray, and the block is cross-approximated (ACA with partial
  pivoting on the exact pair kernel, recompressed by SVD to `tol` relative Frobenius
  error, checked on sampled entries); a block that does not compress within `maxrank`
  terms is refined;
- a shaft holding facets of one convex part only (the capsule, both clusters outside its
  circumscribed ball) resolves each pair exactly by the part's balls: a segment missing
  the circumscribed ball is visible, one through the inscribed ball hidden, and only the
  few between take a ray. The mask makes the block zero, low rank without further rays,
  or refined with it down to dense leaf pairs (a shadow edge);
- any other occupied shaft is refined down to leaf pairs. A leaf pair whose shaft holds
  facets of up to 8 parts with a certificate is decided the same way, pair by pair and
  exactly as the ray would: a planar part (a washer, an LEH cap, a disc) by where the
  segment crosses its plane, inside or outside its facets with margin; a cylinder or
  cone part by its axisymmetric bounds; a convex part by its balls or the plane of an end
  element lying on it. Other leaf pairs take every pair's visibility ray. Wholly hidden
  pairs give zero, wholly visible ones a cross approximation, and partly seen ones (a
  rim, `:sampled` fractions) stay dense. With `visibility = :sampled` the certificates
  resolve a pair's whole bundle of sample rays at once (margins grown by the samples'
  spread), else each of its rays, and the fractions are those of `view_factors`;
- pairs the facing slack admits (rays that start farther off curved facets) are judged
  on their own ray ends when the shaft's exclusions hold with that margin, and take the
  ray otherwise;
- rows cut by the facing test (a sphere's terminator) are refined against the column
  cluster's children beside the low-rank rest (low rank on the rows facing a child
  wholly, dense against a leaf or when few), and a small block with only a few facing
  pairs is stored on those.
Near pairs are dense (Float32), with the same entries as `view_factors` (`_vf_pair`);
blocks whose entries are all zero are not stored.

On the hohlraum presets (10–29k elements) the default `tol = 1e-6` changes the capsule
moments by ≤ 1e-8 and the incident flux by ≤ 2e-6 relative per element against the
Float32 matrix or a tighter `HMatrix` (research notes, "Hierarchical view factors" and,
for the surface-split cluster tree and `tol`, "The H-matrix's N^1.5 remainder"; for the
zero blocks and cut rows, "H-matrix near field"; for the visibility certificates,
"H-matrix visibility certificates").
"""
struct HMatrix <: AbstractMatrix{Float64}
    n::Int
    perm::Vector{Int}                   # tree position → element
    iperm::Vector{Int}                  # element → tree position
    area::Vector{Float64}               # element areas (element order)
    scale::Vector{Float64}              # balancing S (element order)
    dense::Vector{Tuple{Int,Int,Matrix{Float32}}}               # (row0, col0, D), tree order
    lowrank::Vector{Tuple{Int,Int,Matrix{Float32},Matrix{Float32}}}   # (row0, col0, U, V)
    scattered::Vector{Tuple{Vector{Int},Vector{Int},Matrix{Float32}}}  # (rows, cols, D), off the diagonal
    lowrank_index::Vector{Tuple{Vector{Int},Vector{Int},Matrix{Float32},Matrix{Float32}}}  # (rows, cols, U, V)
    groups::Vector{Vector{Tuple{Int,Int}}}    # product work per task: (kind, index)
    info::NamedTuple
end

Base.size(H::HMatrix) = (H.n, H.n)
Base.copy(H::HMatrix) = HMatrix(H.n, H.perm, H.iperm, H.area, copy(H.scale), H.dense, H.lowrank, H.scattered,
                                H.lowrank_index, H.groups, H.info)

# Stored G entry at tree positions (p, q).
function _hentry(H::HMatrix, p::Int, q::Int)
    for (r0, c0, D) in H.dense
        m, n = size(D)
        r0 <= p < r0 + m && c0 <= q < c0 + n && return Float64(D[p-r0+1, q-c0+1])
        r0 <= q < r0 + m && c0 <= p < c0 + n && return Float64(D[q-r0+1, p-c0+1])
    end
    for (rs, cs, D) in H.scattered
        a, b = findfirst(==(p), rs), findfirst(==(q), cs)
        a !== nothing && b !== nothing && return Float64(D[a, b])
        a, b = findfirst(==(q), rs), findfirst(==(p), cs)
        a !== nothing && b !== nothing && return Float64(D[a, b])
    end
    # the indexed blocks before the contiguous ones: a low-rank block's cut rows are zero
    # in it and stored in scattered or indexed blocks
    for (rs, cs, U, V) in H.lowrank_index
        a, b = findfirst(==(p), rs), findfirst(==(q), cs)
        a !== nothing && b !== nothing && return _hlr_entry(U, V, a, b)
        a, b = findfirst(==(q), rs), findfirst(==(p), cs)
        a !== nothing && b !== nothing && return _hlr_entry(U, V, a, b)
    end
    for (r0, c0, U, V) in H.lowrank
        m, n = size(U, 1), size(V, 1)
        a, b = r0 <= p < r0 + m && c0 <= q < c0 + n ? (p - r0 + 1, q - c0 + 1) :
               r0 <= q < r0 + m && c0 <= p < c0 + n ? (q - r0 + 1, p - c0 + 1) : (0, 0)
        a == 0 && continue
        return _hlr_entry(U, V, a, b)
    end
    return 0.0
end

function _hlr_entry(U, V, a, b)
    s = 0.0
    for l in axes(U, 2)
        s += Float64(U[a, l]) * Float64(V[b, l])
    end
    return s
end

function Base.getindex(H::HMatrix, i::Int, j::Int)
    @boundscheck checkbounds(H, i, j)
    ip = H.iperm
    return H.scale[i] * _hentry(H, ip[i], ip[j]) * H.scale[j] / H.area[i]
end

Base.Matrix(H::HMatrix) = _expand(H)
function _expand(H::AbstractMatrix{Float64})
    n = size(H, 1)
    M = Matrix{Float64}(undef, n, n)
    for c0 in 1:256:n
        cols = c0:min(c0 + 255, n)
        E = zeros(n, length(cols))
        for (k, j) in enumerate(cols)
            E[j, k] = 1.0
        end
        M[:, cols] .= H * E
    end
    return M
end
Base.Matrix{Float64}(H::HMatrix) = Matrix(H)
Float32Matrix(H::HMatrix) = Float32Matrix(Matrix{Float32}(Matrix(H)))

function Base.show(io::IO, ::MIME"text/plain", H::HMatrix)
    print(io, H.n, "×", H.n, " HMatrix: ", length(H.dense) + length(H.scattered), " dense and ",
          length(H.lowrank) + length(H.lowrank_index), " low-rank blocks, ", round(H.info.bytes / 2^20; digits = 1), " MiB (",
          round(100 * H.info.bytes / (4.0 * H.n^2); digits = 1), "% of Float32 dense)")
end

_rowsums(H::HMatrix) = H * ones(H.n)
_balance!(H::HMatrix, s::AbstractVector{Float64}) = (H.scale .*= s; H)

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

# Bisection of the centroids along the longest axis of their box, at the median. A
# cluster whose elements lie on several surfaces (`group`) is first split by surface,
# the largest one from the rest, so that no cluster mixes, say, capsule and wall
# elements, or a cap and the wall at a corner (such blocks seldom compress, and a
# capsule element in a wall cluster keeps it from the capsule's ball certificate).
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

# Per element: an index of the surface it lies on (its component's curved surface, or
# its component if it has none), so that coplanar caps and LEH discs share one.
function _hsurface_groups(t::Target)
    keys = Any[]
    gc = map(t.components) do c
        k = c.surface === nothing ? c : c.surface
        i = findfirst(x -> x == k, keys)
        i === nothing ? (push!(keys, k); length(keys)) : i
    end
    return [gc[c] for c in t.elem_comp]
end

# Cluster pairs (s, t), s before t or s = t, covering the upper triangle of G: geometric
# candidates for low rank (`far`, decided per block in assembly) and leaf pairs (`near`).
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

# Visibility between two clusters. The shaft (the convex hull of their boxes) holds every
# segment between their elements; if no facet that could occlude one of those segments
# meets it, the pair is wholly visible wherever it faces (exact, conservative: Haines &
# Wallace's shaft culling). Facets of the two clusters' own elements are not occluders
# (a ray leaves its element's facet, and the receiving facet is its end), nor are those
# of a curved surface both clusters lie inside (`_vf_visible`'s rule).

# Planes (n, d), n·x ≤ d inside, of the convex hull of two boxes: the faces of their
# union box and the planes through an edge of each, parallel to an axis (every face of
# the hull is one of these; Haines & Wallace 1994). Each candidate is kept only if it
# supports both boxes, so the planes bound a region containing the hull in any case.
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
            (pa[k] == alo[k] && pb[k] == blo[k]) || continue    # one edge each, by its lower end
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

# What occludes the shaft between clusters a and b (an `_HOcc`): nothing that could
# occlude a segment between them; only facets of up to 8 components with a certificate
# (a convex part's balls, a planar part, a cylinder's or cone's axisymmetric bounds:
# `_hcert_verdict` decides each segment, or leaves it to the ray); or a facet of another
# component. A facet is not an occluder (for the pairs `_vf_pair_free` evaluates
# without a ray: facing beyond the slack, so leaving their facets; with `:sampled`
# visibility, for every ray between their samples lifted by ε: the plane tests below
# take the box of those points, `opts.qbox`, and the least height of a cluster's samples
# off a plane, `cplane`) if
# - it lies on a curved surface (cylinder, sphere, cone) both clusters lie inside;
# - it belongs to a convex component (every vertex behind every facet) one of the
#   clusters lies on entirely, each element's facets in one plane through its centroid
#   with its normal: a segment leaving a face of a convex polyhedron outward stays
#   outside it;
# - each cluster lies strictly in front of its plane, or is planar in that plane facing
#   the same way: the lifted ray ends are in front, so the segment does not cross it.
# For the pairs admitted by the slack (`_hslack_verdict`, whose rays start farther off
# their facets) it also records whether every plane exclusion holds with that margin
# (`liftx`), and which convex components a cluster lies on were met in the shaft.
# `full`: up to 8 certified components, else only a convex part's balls with both boxes
# outside its circumscribed ball (the search stops at the first other facet: the tests
# for an empty shaft).
function _shaft_occluder(ctx, opts, sa::Int, sb::Int, full::Bool = false)
    (; t, bvh, tri_surf, ε) = ctx
    a, b = opts.nodes[sa], opts.nodes[sb]
    allin, cplane, ccomp, convex = opts.allin, opts.cplane, opts.ccomp, opts.convex
    planes = _shaft_planes(a.bmin, a.bmax, b.bmin, b.bmax)
    outside_box(lo, hi) = any(((n, d),) -> dot(n, ifelse.(n .> 0, lo, hi)) > d, planes)
    # (the ray ends' boxes: centroids, or with `:sampled` visibility the samples lifted)
    ca, cb = _box_corners(opts.qbox[sa]...), _box_corners(opts.qbox[sb]...)
    L = norm(max.(a.cmax, b.cmax) - min.(a.cmin, b.cmin))
    # 0: not clear of the plane; 1: clear; 2: clear also with the slack lifts
    function clear(c, corners, n, d)
        h = minimum(x -> dot(n, x) - d, corners)
        h > 2ε + 1e-12 * L && return h > 2ε + 1e-12 * L + opts.liftx[c] ? 2 : 1
        for (pn, pd, off) in cplane[c]
            dot(pn, n) > 1 - 1e-12 && abs(pd - d) <= 1e-12 * (L + abs(d)) && off > 2ε + 1e-12 * L &&
                return off > 2ε + 1e-12 * L + opts.liftx[c] ? 2 : 1
        end
        return 0
    end
    # components whose facets need no test: convex ones a cluster lies on (never
    # occluders), and the certified ones once a facet of them is in the shaft
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

# A block's entries by its shaft: the empty shaft itself (no rays), else `false` (a ray
# per pair).
_hfree(occ) = occ.n == 0 ? occ : false

# What occludes a shaft (`_shaft_occluder`): `n = 0` nothing, `n = -1` a facet no
# certificate decides (every pair takes a ray), else the components `c[1:n]` whose
# certificates decide most pairs (`_hvis_cert`): a convex part's balls, a planar part's
# facets (`_HPlanar`), a part's axisymmetric bounds (`_HHull`). For the pairs the facing
# slack admits (`_hslack_verdict`): `slack`, whether the plane exclusions hold with their
# rays' lifts, and `x`, the convex components excluded because a cluster lies on them.
struct _HOcc
    n::Int
    c::NTuple{8,Int}
    slack::Bool
    x::NTuple{2,Int}
end
const _HOCC_EMPTY = _HOcc(0, ntuple(_ -> 0, 8), true, (0, 0))
const _HOCC_RAYS = _HOcc(-1, ntuple(_ -> 0, 8), false, (0, 0))

# Whether convex component c's balls decide every pair between clusters a and b but those
# between the balls: c has balls (`occluders`) and both boxes lie outside its
# circumscribed ball. (The one certificate of `_shaft_occluder` without `full`.)
function _hball_whole(opts, c, a::_HCluster, b::_HCluster)
    k = opts.ball[c]
    k == 0 && return false
    _, ctr, _, rout = opts.occluders[k]
    return _box_distance(ctr, a.bmin, a.bmax) > rout && _box_distance(ctr, b.bmin, b.bmax) > rout
end

# A planar component (a washer, a disc, a foil): every facet within `tp` of the plane
# n·(x - c) = 0, in-plane basis (u, v), in-plane distances from c within [rlo, rhi]. Its
# facets as seen in the plane (inward edge lines, radial range) are binned on a polar
# grid (nθ × nρ from ρ0 by dρ; CSR `start`/`idx`), each in the bins its angular span and
# radial range overlap and their neighbors.
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

# A component's axisymmetric bounds about the axis (c, a) of its surface: every facet
# within a·(x - c) ∈ [z0, z1] and at distance [rlo, rhi] from the axis.
struct _HHull
    c::Vec3
    a::Vec3
    z0::Float64
    z1::Float64
    rlo::Float64
    rhi::Float64
end

# Distance from the origin to the triangle (A, B, C) in the plane.
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

_hbasis(n) = (u = normalize(cross(n, abs(n[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0))); (u, cross(n, u)))

# The certificate of component facets `ks` if they are planar (else nothing).
function _hplanar(t::Target, ks, surface)
    n = t.tri_normal[ks[1]]
    all(k -> abs(dot(t.tri_normal[k], n)) > 1 - 1e-12, ks) || return nothing
    tri(k) = (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
    c = surface isa DiskSurface ? surface.center : sum(v for k in ks for v in tri(k)) / (3 * length(ks))
    R = maximum(norm(v - c) for k in ks for v in tri(k))
    tp = maximum(abs(dot(n, v - c)) for k in ks for v in tri(k))
    tp <= 1e-10 * (R + norm(c)) || return nothing
    u, v = _hbasis(n)
    P(x) = SVector(dot(u, x - c), dot(v, x - c))
    nf = length(ks)
    edges = Vector{NTuple{9,Float64}}(undef, nf)
    rr = Vector{NTuple{2,Float64}}(undef, nf)
    span = Vector{NTuple{2,Float64}}(undef, nf)        # angular span (lo, hi), or (-Inf, Inf): all
    for (f, k) in enumerate(ks)
        A, B, C = map(P, tri(k))
        rr[f] = (_origin_tri_distance(A, B, C), max(norm(A), norm(B), norm(C)))
        σ = sign((B - A)[1] * (C - A)[2] - (B - A)[2] * (C - A)[1])
        line(p, q) = (e = q - p; l = norm(e); nx = -σ * e[2] / l; ny = σ * e[1] / l; (nx, ny, nx * p[1] + ny * p[2]))
        # (a degenerate facet gets NaN lines: never decided, so its pairs take the ray)
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

# The axisymmetric bounds of component facets `ks` about the axis (c, a).
function _hhull(t::Target, ks, c, a)
    u, v = _hbasis(a)
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
    return _HHull(c, a, z0, z1, rlo, rhi)
end

# Per component: its certificate (> 0 into `planar`, < 0 into `hull`, 0 none). Convex
# parts with balls keep those; apertures are no occluders.
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

# The ball verdict of the segment oi → oj for a convex part with inscribed ball (ctr,
# rin) and circumscribed ball rout: 1 if it misses the circumscribed ball (so every
# facet), 0 if it passes through the inscribed ball with both ends outside the
# circumscribed one (so it crosses the closed surface), else -1. Margins 1e-9 relative
# on top of those in `rin`, `rout`.
#
# Bundles (every certificate): with `ri`, `rj` > 0 the verdict holds for every segment
# from a point within `ri` of oi to one within `rj` of oj (`:sampled` visibility: the
# segments between two elements' ray samples). Such a segment's point at parameter s lies
# within ρ = max(ri, rj) of the point of oi → oj at s (a tube about it), and its ends'
# heights over a plane within ri, rj of those of oi, oj. Zero radii give the segment's
# verdict, bit for bit as without them.
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

# The planar verdict of the segment oi → oj, as `_vf_visible`'s ray (`any_hit` to
# 1 - 1e-9, barycentric slack 1e-10) finds the facets: 1 if the segment stays on one
# side of the plane or crosses it outside every facet, 0 if it crosses inside a facet
# with its end and start off the plane, else -1. Margins 1e-9 relative (and the plane's
# flatness `tp`), grown by 1/sin of the crossing angle in the plane. A bundle's crossings
# (radii ri, rj) lie within ρ (1 + |w|/|dh|) of oi → oj's (ρ = max(ri, rj)): a point at
# height ≤ ρ on the tube's axis is that far along it from the crossing.
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

# The axisymmetric verdict of the segment oi → oj: 1 if its part within the facets'
# axial range (widened by the margin) lies inside rlo or outside rhi, else -1. A bundle
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
    ww = dot(Rw, Rw)
    s = ww > 0 ? clamp(-dot(Ri, Rw) / ww, sa, sb) : sa
    x = Ri + s * Rw
    dot(x, x) > (H.rhi + ρ + m)^2 && return Int8(1)
    return Int8(-1)
end

# The verdict of component c's certificate on the segment oi → oj between elements i
# and j. A convex part is also cleared by the plane of an end's element lying on it
# (`onconvex`: its facets in that plane): the part lies behind it (within the convexity
# tolerance), so a segment with both ends in front misses it. With radii ri, rj: the
# verdict of the bundle of segments between points within ri of oi and rj of oj, the
# ray samples of elements i and j (`_hvis_bundle`); an end element's own samples lie in
# its plane within `opts.hdev`.
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

# Per BVH node: the component of every facet below it, or -1 if they are of several.
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

# Whether every segment between the boxes a and b crosses the ball (c, r) (exactly, for
# boxes outside the ball). For a point q outside the ball, the points x past the ball as
# seen from q (the segment qx crosses it) include the convex set K_q: the cone from q
# over the ball, beyond q's polar plane (x - c)·(q - c) < r². If every corner of a lies
# in K_q for every corner q of b, each segment from a corner of b to any point of a
# crosses the ball; each corner of b then lies past the ball from any point p of a, so
# in K_p (a point outside the ball past it from p is beyond p's polar plane), and so does
# all of b (K_p is convex): every segment from a to b crosses it.
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
# b from each other: its inscribed ball (about the vertex mean) blocks every segment
# between their boxes, both boxes lying outside its circumscribed ball (so each segment
# enters the component through a facet). Neither cluster may lie on it.
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

# Per cluster: its plane (n, d) if every element is a flat patch of facets in one plane
# with one normal (else a zero normal), and its component if every element is a
# one-facet element of one component (else 0). Per component: whether its mesh is convex
# with outward normals (sphere-surface components only).
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
        # the planes of its flat elements, each with the least height above it of the
        # cluster's other elements (Inf if none): of their centroids, or of the points
        # `pts[e]` (`:sampled` visibility: the ray samples lifted)
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

# The exact visibility of every pair of a block (`_vf_visible`'s verdict, by tree
# positions from (r0, c0)): -1 not facing, 0 hidden, 1 visible.
struct _HMask
    v::Matrix{Int8}
    r0::Int
    c0::Int
end

# The visible fractions (`_vf_visible`'s, -1 if not facing) of every pair of a block, by
# tree positions from (r0, c0) (`:sampled` visibility).
struct _HFrac
    v::Matrix{Float64}
    r0::Int
    c0::Int
end

# G at tree positions (p, q); `ev`: the pair's shaft holds no occluder (an empty `_HOcc`:
# no ray but for some pairs the facing slack admits), `false` (a ray per pair), or the
# block's visibility mask.
function _hentry_eval(ctx, perm, p, q, ::Bool)
    return _vf_exchange(ctx, perm[p], perm[q])
end
function _hentry_eval(ctx, perm, p, q, occ::_HOcc)
    i, j = perm[p], perm[q]
    i == j && return _vf_exchange(ctx, i, j)
    return i < j ? _vf_pair_free(ctx, i, j, occ) : _vf_pair_free(ctx, j, i, occ)
end
function _hentry_eval(ctx, perm, p, q, M::_HMask)
    i, j = minmax(perm[p], perm[q])
    i == j && return _vf_exchange(ctx, i, j)     # (a near-field diagonal block)
    M.v[p-M.r0+1, q-M.c0+1] > 0 || return 0.0
    d = ctx.C[j] - ctx.C[i]
    return _vf_kernel(ctx, i, j, 1.0, d, dot(d, d))
end
function _hentry_eval(ctx, perm, p, q, M::_HFrac)
    i, j = minmax(perm[p], perm[q])
    i == j && return _vf_exchange(ctx, i, j)
    f = M.v[p-M.r0+1, q-M.c0+1]
    f > 0 || return 0.0
    d = ctx.C[j] - ctx.C[i]
    return _vf_kernel(ctx, i, j, f, d, dot(d, d))
end

# The ray ends of `_vf_visible` for the pair i < j (centroids at distance r, facing
# cosines ci, cj) that only the facing slack admits, and whether its shaft's exclusions
# (`_shaft_occluder`) hold for that segment: the ends stay inside the clusters' boxes
# (lifted by at most the boxes' padding), the plane exclusions had the lifts' margin
# (`occ.slack`), and each convex component a cluster lies on is cleared by the plane of
# an end's element on it (`ctx.onplane`: its facets in that plane, which the convex part
# lies behind) with both ends in front by `ctx.hsep`.
function _hslack_ends(ctx, i, j, r, ci, cj, occ::_HOcc)
    (; t, C, Nrm, size, ε, onplane, hsep) = ctx
    li, lj = ε + size[i] * max(0.0, -ci), ε + size[j] * max(0.0, -cj)
    oi, oj = C[i] + li * Nrm[i], C[j] + lj * Nrm[j]
    ok = occ.slack && li <= size[i] / sqrt(2) && lj <= size[j] / sqrt(2)
    for x in occ.x
        (ok && x != 0) || continue
        ok = (t.elem_comp[i] == x && onplane[i] && dot(Nrm[i], oj - C[i]) > hsep && dot(Nrm[i], oi - C[i]) > hsep) ||
             (t.elem_comp[j] == x && onplane[j] && dot(Nrm[j], oi - C[j]) > hsep && dot(Nrm[j], oj - C[j]) > hsep)
    end
    return ok, oi, oj
end

# Visibility of the pair i < j whose shaft holds facets of the certified components
# `occ.c[1:occ.n]` only (`_shaft_occluder`), as `_vf_visible`'s ray finds it: hidden if
# one certificate finds the segment through its facets, visible if every one finds it
# clear (`_hcert_verdict`: a convex part's balls, a planar part's facets, a part's
# axisymmetric bounds). A pair the facing slack admits is judged on its own ray ends
# (`_hslack_ends`). Only pairs left undecided take the ray. Centroid visibility only.
function _hvis_cert(ctx, opts, i, j, occ::_HOcc)
    (; C, Nrm, slack, ε) = ctx
    d = C[j] - C[i]
    r2 = dot(d, d)
    r2 == 0 && return Int8(-1)
    r = sqrt(r2)
    ci, cj = dot(Nrm[i], d) / r, -dot(Nrm[j], d) / r
    ci > -slack[i] && cj > -slack[j] || return Int8(-1)
    if ci > slack[i] && cj > slack[j]
        # the ray ends of `_vf_visible` (lifted by ε)
        ok, oi, oj = true, C[i] + ε * Nrm[i], C[j] + ε * Nrm[j]
    else
        ok, oi, oj = _hslack_ends(ctx, i, j, r, ci, cj, occ)
    end
    if ok
        clear = true
        for k in 1:occ.n
            v = _hcert_verdict(ctx, opts, occ.c[k], i, j, oi, oj)
            v == 0 && return Int8(0)
            clear &= v > 0
        end
        clear && return Int8(1)
    end
    return first(_vf_visible(ctx, i, j)) > 0 ? Int8(1) : Int8(0)
end

# The visibility mask of the block rows × cols whose shaft `occ` certifies.
function _hmask(ctx, opts, perm, rows::UnitRange{Int}, cols::UnitRange{Int}, occ::_HOcc)
    V = Matrix{Int8}(undef, length(rows), length(cols))
    for (b, q) in enumerate(cols), (a, p) in enumerate(rows)
        V[a, b] = _hvis_cert(ctx, opts, minmax(perm[p], perm[q])..., occ)
    end
    return _HMask(V, first(rows), first(cols))
end

# The `:sampled` visible fraction of the pair i < j whose shaft holds facets of the
# certified components `occ.c[1:occ.n]` only (`_shaft_occluder`, searched with the
# clusters' sample boxes), `===` `_vf_visible`'s. Its rays (the centroid ray, the probes',
# those between the samples) all run between points within `opts.rvis` of the two ray
# ends, so each certificate decides the bundle (`_hcert_verdict` with radii): hidden if
# one finds every segment through its facets, visible if all find every segment clear.
# Otherwise `_visible_fraction` runs as `_vf_visible` does, each segment decided by the
# certificates that did not clear the bundle (`_HSeen`), or by its ray. The shaft's
# exclusions must hold for the samples: their ends inside the cluster boxes (rvis + the
# lift within the padding), the plane exclusions with the slack lifts' margin for a pair
# the slack admits (`occ.slack`), and each convex part a cluster lies on cleared by an
# end element's plane (per bundle, else per segment). Folded elements (whose own facets
# the samples' rays skip) take `_vf_visible`.
function _hvis_bundle(ctx, opts, i, j, occ::_HOcc)
    (; t, bvh, C, Nrm, size, slack, ε, tri_surf, inside, vis, multi, onplane, hsep) = ctx
    facing, _, _, ci, cj, oi, oj = _vf_ends(C, Nrm, size, slack, ε, i, j)    # `_vf_visible`'s
    facing || return -1.0
    (multi[i] || multi[j]) && return first(_vf_visible(ctx, i, j))
    li, lj = oi - C[i], oj - C[j]
    ri, rj = opts.rvis[i], opts.rvis[j]
    ok = ((ci > slack[i] && cj > slack[j]) || occ.slack) &&
         ri + norm(li) <= size[i] / sqrt(2) && rj + norm(lj) <= size[j] / sqrt(2)
    ok || return first(_vf_visible(ctx, i, j))
    # convex parts a cluster lies on: cleared for the bundle, else checked per segment
    xs = (0, 0)
    for (k, x) in enumerate(occ.x)
        x == 0 && continue
        whole = (t.elem_comp[i] == x && onplane[i] && dot(Nrm[i], li) - opts.hdev[i] > hsep &&
                 dot(Nrm[i], oj - C[i]) - rj > hsep) ||
                (t.elem_comp[j] == x && onplane[j] && dot(Nrm[j], lj) - opts.hdev[j] > hsep &&
                 dot(Nrm[j], oi - C[j]) - ri > hsep)
        whole || (xs = Base.setindex(xs, x, k))
    end
    vb = ntuple(_ -> Int8(1), 8)
    clear = true
    for k in 1:occ.n
        v = _hcert_verdict(ctx, opts, occ.c[k], i, j, oi, oj, ri, rj)
        v == 0 && return 0.0
        vb = Base.setindex(vb, v, k)
        clear &= v > 0
    end
    clear && xs == (0, 0) && return 1.0
    skip = k -> (s = tri_surf[k]; s != 0 && inside[s, i] && inside[s, j])
    seen = _HSeen(ctx, opts, occ, vb, xs, i, j, _RaySeen(bvh, skip))
    return _visible_fraction(seen, vis[i], vis[j], oi, oj, li, lj, Nrm[i], Nrm[j], !seen(oi, oj))
end

# Whether the segment p → q between the ray samples of elements i and j is unobstructed:
# by the certificates of `occ` that did not clear the bundle (`vb[k] != 1`), once the
# convex parts `xs` a cluster lies on are cleared by an end element's plane; else by the
# ray.
struct _HSeen{C,O,R}
    ctx::C
    opts::O
    occ::_HOcc
    vb::NTuple{8,Int8}
    xs::NTuple{2,Int}
    i::Int
    j::Int
    ray::R
end
function (s::_HSeen)(p, q)
    (; ctx, opts, occ, vb, xs, i, j) = s
    (; t, C, Nrm, onplane, hsep) = ctx
    for x in xs
        x == 0 && continue
        (t.elem_comp[i] == x && onplane[i] && dot(Nrm[i], p - C[i]) > hsep && dot(Nrm[i], q - C[i]) > hsep) ||
            (t.elem_comp[j] == x && onplane[j] && dot(Nrm[j], p - C[j]) > hsep && dot(Nrm[j], q - C[j]) > hsep) ||
            return s.ray(p, q)
    end
    clear = true
    for k in 1:occ.n
        vb[k] == 1 && continue
        v = _hcert_verdict(ctx, opts, occ.c[k], i, j, p, q)
        v == 0 && return false
        clear &= v > 0
    end
    return clear || s.ray(p, q)
end

# The visible fractions of the block rows × cols whose shaft `occ` certifies.
function _hfrac(ctx, opts, perm, rows::UnitRange{Int}, cols::UnitRange{Int}, occ::_HOcc)
    V = Matrix{Float64}(undef, length(rows), length(cols))
    for (b, q) in enumerate(cols), (a, p) in enumerate(rows)
        V[a, b] = _hvis_bundle(ctx, opts, minmax(perm[p], perm[q])..., occ)
    end
    return _HFrac(V, first(rows), first(cols))
end

# The certified visibility of a block: a mask (centroid visibility) or fractions.
_hcertified(ctx, opts, perm, rows, cols, occ) =
    ctx.vis === nothing ? _hmask(ctx, opts, perm, rows, cols, occ) : _hfrac(ctx, opts, perm, rows, cols, occ)

# `_vf_pair` of a pair i < j whose shaft `occ` holds no occluder: a pair facing beyond the
# slack of the facing test leaves its elements' facets, and its ray would be
# unobstructed; one admitted by the slack is visible if its own ray ends keep the
# shaft's exclusions (`_hslack_ends`; centroid visibility), else it takes the ray.
function _vf_pair_free(ctx, i, j, occ::_HOcc)
    (; C, Nrm, slack) = ctx
    d = C[j] - C[i]
    r2 = dot(d, d)
    r2 == 0 && return 0.0
    r = sqrt(r2)
    ci, cj = dot(Nrm[i], d) / r, -dot(Nrm[j], d) / r
    ci > -slack[i] && cj > -slack[j] || return 0.0
    if !(ci > slack[i] && cj > slack[j])
        (ctx.vis === nothing && first(_hslack_ends(ctx, i, j, r, ci, cj, occ))) || return _vf_pair(ctx, i, j)
    end
    return _vf_kernel(ctx, i, j, 1.0, d, r2)
end

# Dense block of G at tree positions rows × cols (Float32); a diagonal block is filled
# from its upper triangle.
function _hdense(ctx, perm, rows::UnitRange{Int}, cols::UnitRange{Int}, free)
    D = Matrix{Float32}(undef, length(rows), length(cols))
    if rows == cols
        for b in eachindex(cols), a in 1:b
            g = Float32(_hentry_eval(ctx, perm, rows[a], cols[b], free))
            D[a, b] = g
            D[b, a] = g
        end
    else
        for b in eachindex(cols), a in eachindex(rows)
            D[a, b] = Float32(_hentry_eval(ctx, perm, rows[a], cols[b], free))
        end
    end
    return D
end

# Dense block at tree positions rows × cols (index lists of two different clusters).
function _hdense_index(ctx, perm, rows::Vector{Int}, cols::Vector{Int}, free)
    D = Matrix{Float32}(undef, length(rows), length(cols))
    for b in eachindex(cols), a in eachindex(rows)
        D[a, b] = Float32(_hentry_eval(ctx, perm, rows[a], cols[b], free))
    end
    return D
end

# Dense block rows × cols (s ≠ t) whose visibility `vis` (`_vf_visible`'s fraction, by
# position in the block) is known.
function _hdense_masked(ctx, perm, rows, cols, vis)
    D = zeros(Float32, length(rows), length(cols))
    for b in eachindex(cols), a in eachindex(rows)
        f = vis[a, b]
        f > 0 || continue
        i, j = minmax(perm[rows[a]], perm[cols[b]])
        d = ctx.C[j] - ctx.C[i]
        D[a, b] = Float32(_vf_kernel(ctx, i, j, f, d, dot(d, d)))
    end
    return D
end

# The positions in `rows` whose element may face some centroid of cluster c: those
# with a corner of c's centroid box in front of the element's plane (by the facing test's
# slack); the others' rows of the block are zero. And those that face all of c's box
# beyond the slack: their rows have no facing cut-off (a kink in the kernel that the
# cross approximation resolves slowly, at a sphere's terminator).
_hactive(ctx, perm, rows, c::_HCluster) = _hactive(ctx, perm, rows, c.cmin, c.cmax)
function _hactive(ctx, perm, rows, cmin::Vec3, cmax::Vec3)
    corners = _box_corners(cmin, cmax)
    act, full = Int[], Int[]
    for p in rows
        e = perm[p]
        n, x0, sl = ctx.Nrm[e], ctx.C[e], ctx.slack[e]
        mx = maximum(dot(n, x - x0) for x in corners)
        mn = minimum(dot(n, x - x0) for x in corners)
        far = maximum(norm(x - x0) for x in corners)
        mx > -sl * far && push!(act, p)
        mn > sl * far && push!(full, p)
    end
    return act, full
end

# The rows and columns (of `rows` × `cols`) with some pair that passes the facing test.
function _hfacing(ctx, perm, rows, cols)
    (; C, Nrm, slack) = ctx
    rk, ck = falses(length(rows)), falses(length(cols))
    for (b, q) in enumerate(cols), (a, p) in enumerate(rows)
        rk[a] && ck[b] && continue
        i, j = perm[p], perm[q]
        d = C[j] - C[i]
        r = norm(d)
        r > 0 || continue
        dot(Nrm[i], d) > -slack[i] * r && -dot(Nrm[j], d) > -slack[j] * r || continue
        rk[a] = ck[b] = true
    end
    return rows[rk], cols[ck]
end

# Adaptive cross approximation with partial pivoting of the block rows × cols of G to
# relative Frobenius accuracy `tol`, at most `kmax` terms. Returns (U, V, converged).
function _haca(ctx, perm, rows, cols, tol, kmax, i0, free)
    m, n = length(rows), length(cols)
    Us, Vs = Vector{Float64}[], Vector{Float64}[]
    used = falses(m)
    row = zeros(n)
    i = i0
    nrm2 = 0.0
    zeros_seen = 0
    hist = Float64[]
    while length(Us) < kmax
        used[i] = true
        for b in 1:n
            row[b] = _hentry_eval(ctx, perm, rows[i], cols[b], free)
        end
        for l in eachindex(Us)
            u = Us[l][i]
            u == 0 || (row .-= u .* Vs[l])
        end
        j = argmax(abs.(row))
        if row[j] == 0
            # a zero residual row: the block may be zero (rank 0) or this row alone
            zeros_seen += 1
            left = findall(!, used)
            (zeros_seen > 3 || isempty(left)) && return Us, Vs, true
            i = left[cld(length(left), 2)]
            continue
        end
        v = row ./ row[j]
        u = Vector{Float64}(undef, m)
        for a in 1:m
            u[a] = _hentry_eval(ctx, perm, rows[a], cols[j], free)
        end
        for l in eachindex(Us)
            vl = Vs[l][j]
            vl == 0 || (u .-= vl .* Us[l])
        end
        uu, vv = dot(u, u), dot(v, v)
        cross = 0.0
        for l in eachindex(Us)
            cross += dot(u, Us[l]) * dot(v, Vs[l])
        end
        nrm2 += uu * vv + 2cross
        push!(Us, u)
        push!(Vs, v)
        uu * vv <= tol^2 * nrm2 && return Us, Vs, true
        # give up early when the terms' decay, extrapolated, misses tol at kmax: a block
        # that will not compress is refined, and its children pay for it again
        e = sqrt(uu * vv / nrm2)
        push!(hist, e)
        k = length(hist)
        if k >= 8
            rate = min(1.0, (hist[k] / hist[k-4])^(1 / 4))
            hist[k] * rate^(kmax - k) > 3tol && return Us, Vs, false
        end
        best, i = -1.0, 0
        for a in 1:m
            used[a] && continue
            abs(u[a]) > best && (best = abs(u[a]); i = a)
        end
        i == 0 && return Us, Vs, true
    end
    return Us, Vs, false
end

# Recompress U Vᵀ (QR of both factors, SVD of the small core) to relative accuracy tol.
function _hrecompress(Us, Vs, tol)
    k = length(Us)
    k == 0 && return zeros(Float32, 0, 0), zeros(Float32, 0, 0)
    U, V = reduce(hcat, Us), reduce(hcat, Vs)
    qu, qv = qr(U), qr(V)
    Ru, Rv = qu.R, qv.R
    F = svd(Ru * Rv')
    σ = F.S
    tot = sum(abs2, σ)
    r = k
    tail = 0.0
    while r > 1 && tail + σ[r]^2 <= (tol / 2)^2 * tot
        tail += σ[r]^2
        r -= 1
    end
    Qu = Matrix(qu.Q)[:, 1:k]
    Qv = Matrix(qv.Q)[:, 1:k]
    Uc = Qu * (F.U[:, 1:r] .* σ[1:r]')
    Vc = Qv * F.V[:, 1:r]
    return Matrix{Float32}(Uc), Matrix{Float32}(Vc)
end

# Estimated Frobenius error of U Vᵀ against the block, from `ns` sampled entries, plus the
# error of the entries where a row's or a column's quadrature tier switches
# (`_htier_entries`): G jumps there by the quadrature rules' difference, and an isolated
# jump (one pair of a capsule's polar sliver past the far boundary) is missed by the
# cross approximation's pivots and the random samples alike.
function _hcheck(ctx, perm, rows, cols, U, V, ns, seed, free)
    m, n = length(rows), length(cols)
    rng = Xoshiro(seed)
    e2 = 0.0
    g2 = 0.0
    for _ in 1:ns
        a, b = rand(rng, 1:m), rand(rng, 1:n)
        g = _hentry_eval(ctx, perm, rows[a], cols[b], free)
        s = 0.0
        for l in axes(U, 2)
            s += Float64(U[a, l]) * Float64(V[b, l])
        end
        e2 += (g - s)^2
        g2 += g^2
    end
    x2 = 0.0
    for (a, b) in _htier_entries(ctx, perm, rows, cols)
        g = _hentry_eval(ctx, perm, rows[a], cols[b], free)
        s = 0.0
        for l in axes(U, 2)
            s += Float64(U[a, l]) * Float64(V[b, l])
        end
        x2 += (g - s)^2
    end
    return sqrt(e2 * m * n / ns + x2), sqrt(g2 * m * n / ns)
end

# The positions (a, b) of the block rows × cols where G may jump: in each row (column)
# whose pairs fall in more than one quadrature tier (`_vf_tier`), the pairs outside its
# most common tier; at most m + n of them. None, without a pass over the pairs, if the
# centroid boxes are farther apart than `mid` times the largest element (all far).
function _htier_entries(ctx, perm, rows, cols)
    out = Tuple{Int,Int}[]
    (; C, size, mid2) = ctx
    m, n = length(rows), length(cols)
    alo, ahi = _hcentroid_box(C, perm, rows)
    blo, bhi = _hcentroid_box(C, perm, cols)
    gap = max.(zero(Vec3), max.(blo - ahi, alo - bhi))
    smax = max(maximum(p -> size[perm[p]], rows), maximum(q -> size[perm[q]], cols))
    dot(gap, gap) >= mid2 * smax^2 && return out
    code(s) = s === :near ? 0x01 : s === :mid ? 0x02 : 0x03
    T = Matrix{UInt8}(undef, m, n)
    for (b, q) in enumerate(cols), (a, p) in enumerate(rows)
        i, j = perm[p], perm[q]
        d = C[j] - C[i]
        T[a, b] = code(_vf_tier(ctx, i, j, dot(d, d)))
    end
    # lines (a row a > 0 or a column -a) with their pairs outside the most common tier,
    # the fewest first: an isolated jump is the one the cross approximation misses
    lines = Tuple{Int,Int,UInt8}[]
    function line!(k, line)
        c1, c2, c3 = count(==(0x01), line), count(==(0x02), line), count(==(0x03), line)
        (c1 > 0) + (c2 > 0) + (c3 > 0) > 1 || return
        maj = c1 >= c2 && c1 >= c3 ? 0x01 : c2 >= c3 ? 0x02 : 0x03
        push!(lines, (length(line) - max(c1, c2, c3), k, maj))
    end
    for a in 1:m
        line!(a, view(T, a, :))
    end
    for b in 1:n
        line!(-b, view(T, :, b))
    end
    sort!(lines)
    seen = falses(m, n)
    for (_, k, maj) in lines
        for x in 1:(k > 0 ? n : m)
            a, b = k > 0 ? (k, x) : (x, -k)
            (T[a, b] == maj || seen[a, b]) && continue
            seen[a, b] = true
            push!(out, (a, b))
            length(out) >= m + n && return out
        end
    end
    return out
end

# The exact visibility of the block rows × cols by a ray per pair (`_vf_visible`) on the
# facing positions ra × cb (the others do not face); centroid visibility only.
function _hmask_rays(ctx, perm, rows::UnitRange{Int}, cols::UnitRange{Int}, ra, cb)
    V = fill(Int8(-1), length(rows), length(cols))
    for q in cb, p in ra
        f = first(_vf_visible(ctx, minmax(perm[p], perm[q])...))
        V[p-first(rows)+1, q-first(cols)+1] = f < 0 ? Int8(-1) : f > 0 ? Int8(1) : Int8(0)
    end
    return _HMask(V, first(rows), first(cols))
end

# Visible and hidden facing pairs of a mask on positions ra × cb.
function _hmask_counts(M::_HMask, ra, cb)
    nv = nh = 0
    for q in cb, p in ra
        f = M.v[p-M.r0+1, q-M.c0+1]
        nv += f > 0
        nh += f == 0
    end
    return nv, nh
end
# (fractions: visible, and hidden or partly visible)
function _hmask_counts(M::_HFrac, ra, cb)
    nv = nh = 0
    for q in cb, p in ra
        f = M.v[p-M.r0+1, q-M.c0+1]
        nv += f > 0
        nh += 0 <= f < 1
    end
    return nv, nh
end

# The blocks of cluster pair (s, t) into `dense`/`lowrank`: low rank when the pair is
# admissible, of uniform visibility and its cross approximation converges and passes a
# sampled check; otherwise its children, down to dense leaf pairs. `known`: what a parent
# block found of the pair's visibility (`true`: an empty shaft; a mask), else `nothing`.
# `sub`: the rows of s that a facing test cut against a parent of t (`split` below),
# refined over t only: against t's smaller boxes most of them face wholly (low rank) or
# not at all. Blocks go into `out` = (dense, lowrank, scattered, lowrank_index).
function _hassemble!(out, ctx, nodes, perm, s, t, far, opts, known = nothing,
                     sub::Union{Nothing,Vector{Int}} = nothing)
    dense, lowrank, scattered, lowrank_index = out
    a, b = nodes[s], nodes[t]
    rows, cols = sub === nothing ? (a.lo:a.hi) : sub, b.lo:b.hi
    m, n = length(rows), length(cols)
    st = opts.stats
    leaves = _hleaf(a) && _hleaf(b)
    vis = nothing        # visible fractions of a leaf pair (`:sampled` visibility)
    ev = false           # entries: `true` (empty shaft, no rays), `false` (rays), a mask
    ra, cb = Int[], Int[]
    if far
        Threads.atomic_add!(st.candidates, 1)
        # rows (columns) whose element faces no centroid of the other cluster are zero
        amin, amax = sub === nothing ? (a.cmin, a.cmax) : _hcentroid_box(ctx.C, perm, sub)
        (ra, fa), (cb, fb) = _hactive(ctx, perm, rows, b), _hactive(ctx, perm, cols, amin, amax)
        if !isempty(ra) && !isempty(cb) && length(ra) * length(cb) <= opts.facing_size
            # small blocks: the facing test pair by pair (no rays), which leaves the rows
            # and columns of some facing pair; a block facing itself in a few pairs only
            # (a sphere's terminator) is kept dense on those
            ra, cb = _hfacing(ctx, perm, ra, cb)
            fa, fb = intersect(fa, ra), intersect(fb, cb)
            if !isempty(ra) && 4 * length(ra) * length(cb) <= m * n
                ev = known !== nothing ? known : _hfree(_shaft_occluder(ctx, opts, s, t))
                _hpush_nonzero!(scattered, (ra, cb, _hdense_index(ctx, perm, ra, cb, ev)))
                Threads.atomic_add!(st.sparse, 1)
                return
            end
        end
        if isempty(ra) || isempty(cb)
            Threads.atomic_add!(st.facing, 1)
            return
        end
        if known !== nothing
            ev = known
        else
            occ = _shaft_occluder(ctx, opts, s, t, true)
            if occ.n != 0 && _convex_hides(opts, s, t)
                Threads.atomic_add!(st.shadowed, 1)
                return                                   # wholly behind a convex part
            end
            if occ.n > 0 && m * n <= opts.mask_size
                # only certified parts' facets in the shaft: their certificates decide
                # most pairs (`:sampled`: whole sample bundles), at any level (a block of
                # uniform visibility then gets the cross approximation instead of being
                # refined)
                ev = _hcertified(ctx, opts, perm, rows, cols, occ)
                Threads.atomic_add!(st.masked, 1)
            elseif occ.n != 0 && leaves
                # the exact visibility of every pair: a ray each
                if ctx.vis === nothing
                    ev = _hmask_rays(ctx, perm, rows, cols, ra, cb)
                else
                    vis = fill(-1.0, m, n)
                    for q in cb, p in ra
                        vis[p-a.lo+1, q-b.lo+1] = first(_vf_visible(ctx, minmax(perm[p], perm[q])...))
                    end
                    if !any(>(0), vis)
                        Threads.atomic_add!(st.hidden, 1)
                        return
                    end
                    all(f -> f < 0 || f == 1, vis) || (Threads.atomic_add!(st.mixed, 1); @goto dense)
                end
            else
                ev = _hfree(occ)
            end
        end
        if ev isa _HMask || ev isa _HFrac
            # zero if wholly hidden, low rank if wholly visible (no more rays), else
            # refined with the mask down to dense leaf pairs
            nv, nh = _hmask_counts(ev, ra, cb)
            if nv == 0
                Threads.atomic_add!(st.hidden, 1)
                return
            end
            if nh > 0
                leaves && Threads.atomic_add!(st.mixed, 1)
                @goto dense
            end
        end
        if ev !== false || vis !== nothing
            # rows and columns facing the other cluster wholly get the cross approximation,
            # the few cut by the facing test (if they are few) are kept dense
            split = (length(fa) < length(ra) || length(fb) < length(cb)) &&
                    2length(fa) >= length(ra) && 2length(fb) >= length(cb) && !isempty(fa) && !isempty(fb)
            ar, ac = split ? (fa, fb) : (ra, cb)
            ma, na = length(ar), length(ac)
            kmax = min(opts.maxrank, (3 * ma * na) ÷ (4 * (ma + na)))
            if kmax >= 1
                i0 = something(findfirst(==(opts.iperm[a.center]), ar), cld(ma, 2))
                Us, Vs, ok = _haca(ctx, perm, ar, ac, opts.tol, kmax, i0, ev)
                # a cross approximation that found nothing does not prove a zero block
                ok &= !isempty(Us)
                if ok
                    Uc, Vc = _hrecompress(Us, Vs, opts.tol)
                    ns = min(ma * na, max(32, (ma + na) ÷ 2))
                    seed = sub === nothing ? hash((a.lo, b.lo)) : hash((ar[1], b.lo, ma))
                    err, nrm = _hcheck(ctx, perm, ar, ac, Uc, Vc, ns, seed, ev)
                    if err <= 10 * opts.tol * nrm
                        if size(Uc, 2) > 0 && sub !== nothing
                            push!(lowrank_index, (ar, ac, Uc, Vc))
                        elseif size(Uc, 2) > 0
                            U, V = zeros(Float32, m, size(Uc, 2)), zeros(Float32, n, size(Vc, 2))
                            U[ar .- (a.lo - 1), :] .= Uc
                            V[ac .- (b.lo - 1), :] .= Vc
                            push!(lowrank, (a.lo, b.lo, U, V))
                        end
                        if split
                            pa, pb = setdiff(ra, fa), setdiff(cb, fb)
                            if isempty(pa)
                            elseif _hleaf(b) || length(pa) <= 4 || !(ev isa _HOcc || ev isa _HMask || ev isa _HFrac)
                                _hpush_nonzero!(scattered, (pa, cb, _hdense_index(ctx, perm, pa, cb, ev)))
                            else
                                # the cut rows against t's children, where most face wholly
                                # or not at all
                                for k in b.kids
                                    _hassemble!(out, ctx, nodes, perm, s, k, true, opts, ev, pa)
                                end
                            end
                            isempty(pb) || _hpush_nonzero!(scattered, (fa, pb, _hdense_index(ctx, perm, fa, pb, ev)))
                        end
                        Threads.atomic_add!(st.lowrank, 1)
                        return
                    end
                    Threads.atomic_add!(st.rejected, 1)
                else
                    Threads.atomic_add!(st.unconverged, 1)
                end
            end
        end
    end
    @label dense
    if sub !== nothing
        # cut rows: dense against a leaf or if they are few (low rank would not pay), else
        # refined over t
        if _hleaf(b) || length(ra) <= 4
            _hpush_nonzero!(scattered, (ra, cb, _hdense_index(ctx, perm, ra, cb, ev)))
        else
            for k in b.kids
                _hassemble!(out, ctx, nodes, perm, s, k, true, opts, ev, ra)
            end
        end
    elseif leaves
        if vis !== nothing
            D = _hdense_masked(ctx, perm, rows, cols, vis)
        elseif far
            D = _hdense(ctx, perm, rows, cols, ev)
        else
            # a near leaf pair: no ray in an empty shaft, the certificates' mask in a
            # certified one, else a ray per pair
            occ = _shaft_occluder(ctx, opts, s, t, true)
            D = _hdense(ctx, perm, rows, cols, occ.n == 0 ? occ :
                        occ.n > 0 ? _hcertified(ctx, opts, perm, rows, cols, occ) : false)
        end
        any(!iszero, D) && push!(dense, (a.lo, b.lo, D))     # (a capsule's near field is zero)
    else
        # what this block found of the visibility holds for its children
        kn = far && ev !== false ? ev : nothing
        if _hleaf(b) || (!_hleaf(a) && m >= n)
            for k in a.kids
                _hassemble!(out, ctx, nodes, perm, k, t, far, opts, kn)
            end
        else
            for k in b.kids
                _hassemble!(out, ctx, nodes, perm, s, k, far, opts, kn)
            end
        end
    end
    return
end

# A scattered block, unless all its entries are zero (rows cut at a sphere's terminator
# against the far side of a box).
_hpush_nonzero!(v, blk) = (any(!iszero, blk[3]) && push!(v, blk); v)

# Bounding box of the centroids at tree positions ps.
function _hcentroid_box(C, perm, ps)
    lo = hi = C[perm[ps[1]]]
    for p in ps
        lo, hi = min.(lo, C[perm[p]]), max.(hi, C[perm[p]])
    end
    return lo, hi
end

"""
    Briareus.HMATRIX_DEFAULTS

Defaults of the hierarchical view factors (`view_factors(...; storage = :hmatrix)`):
`tol` (relative Frobenius accuracy of each low-rank block), `eta` (admissibility:
cluster boxes at least `eta` × the larger diameter apart), `leaf` (largest leaf
cluster), `maxrank` (most terms of a cross approximation before the block is refined).
Engineering defaults, from the work–precision sweeps in the research notes
("Hierarchical view factors"; `tol` 1e-5 → 1e-6 in "The H-matrix's N^1.5 remainder").
"""
const HMATRIX_DEFAULTS = (tol = 1e-6, eta = 1.0, leaf = 48, maxrank = 32)

function _hmatrix_view_factors(ctx; tol = HMATRIX_DEFAULTS.tol, eta = HMATRIX_DEFAULTS.eta,
                               leaf::Integer = HMATRIX_DEFAULTS.leaf,
                               maxrank::Integer = HMATRIX_DEFAULTS.maxrank)
    0 < tol < 1 || throw(ArgumentError("tol must be in (0, 1), got $tol"))
    eta > 0 || throw(ArgumentError("eta must be positive, got $eta"))
    leaf >= 2 || throw(ArgumentError("leaf must be at least 2, got $leaf"))
    t = ctx.t
    ctx, nodes, perm, far, near, opts = _hsetup(ctx, tol, eta, leaf, maxrank)
    tasks = [[(p, true) for p in far]; [(p, false) for p in near]]
    # largest first
    cost((p, _)) = (nodes[p[1]].hi - nodes[p[1]].lo + 1) * (nodes[p[2]].hi - nodes[p[2]].lo + 1)
    sort!(tasks; by = cost, rev = true)
    nt = Threads.nthreads()
    outs = [(Tuple{Int,Int,Matrix{Float32}}[], Tuple{Int,Int,Matrix{Float32},Matrix{Float32}}[],
             Tuple{Vector{Int},Vector{Int},Matrix{Float32}}[],
             Tuple{Vector{Int},Vector{Int},Matrix{Float32},Matrix{Float32}}[]) for _ in 1:nt]
    next = Threads.Atomic{Int}(1)
    @sync for w in 1:min(nt, max(length(tasks), 1))
        Threads.@spawn begin
            _hwork!(outs[w], ctx, nodes, perm, tasks, next, opts)
        end
    end
    return _hfinish(t, perm, opts, outs, tol, eta, leaf, maxrank)
end

# The assembly context (the view-factor context with sector receivers vectorized and the
# data of `_hslack_ends`), cluster tree, block pairs and assembly data (visibility
# certificates) of a view-factor context.
function _hsetup(ctx, tol, eta, leaf, maxrank)
    t = ctx.t
    # per element: whether its facets lie in one plane through its centroid with its normal
    onplane = map(1:nelements(t)) do e
        n, C = t.elem_normal[e], t.elem_centroid[e]
        δ = 1e-12 * (norm(C) + 1)
        all(k -> dot(t.tri_normal[k], n) > 1 - 1e-12 && abs(dot(n, t.tri_v0[k] - C)) <= δ, t.elem_tris[e])
    end
    # sectors receive through the vectorized contour (`_SoAReceiver`; round-off); `hsep`:
    # the margin of a plane clearing a convex part (10× its convexity tolerance, 1e-9 of
    # its size, and a tenth of the rays' lift ε)
    ctx = merge(ctx, (; recv = Union{_Receiver,_SoAReceiver,_PolyReceiver}[_soa_receiver(r) for r in ctx.recv],
                      onplane = BitVector(onplane), hsep = 1e-8 * max(t.scale, 1.0)))
    pad = ctx.size ./ sqrt(2)
    nodes, perm = _hcluster_tree(ctx.C, pad, Int(leaf), _hsurface_groups(t))
    iperm = invperm(perm)
    far, near = _hblock_pairs(nodes, eta)
    stats = NamedTuple{(:candidates, :lowrank, :facing, :sparse, :shadowed, :masked, :hidden, :mixed, :unconverged,
                        :rejected)}(ntuple(_ -> Threads.Atomic{Int}(0), 10))
    # allin[k, c]: every element of cluster c inside curved surface k
    ns = size(ctx.inside, 1)
    allin = BitMatrix(undef, ns, length(nodes))
    for (c, nd) in enumerate(nodes), k in 1:ns
        allin[k, c] = all(p -> ctx.inside[k, perm[p]], nd.lo:nd.hi)
    end
    # `:sampled` visibility: per element the ends of its sample rays lifted by ε (as
    # `_vf_visible` lifts them for a pair facing beyond the slack), their largest distance
    # from the lifted centroid (`rvis`) and the samples' from the element's plane (`hdev`);
    # per cluster the box of those ends, for the shaft's plane tests (`qbox`)
    N = nelements(t)
    if ctx.vis === nothing
        pts = nothing
        rvis, hdev = zeros(N), zeros(N)
        qbox = [(nd.cmin, nd.cmax) for nd in nodes]
    else
        pts = [push!([p + ctx.ε * ctx.Nrm[e] for p in ctx.vis[e].points], ctx.C[e] + ctx.ε * ctx.Nrm[e]) for e in 1:N]
        rvis = [maximum(p -> norm(p - ctx.C[e]), ctx.vis[e].points) for e in 1:N]
        hdev = [maximum(p -> abs(dot(ctx.Nrm[e], p - ctx.C[e])), ctx.vis[e].points) for e in 1:N]
        qbox = map(nodes) do nd
            ps = (x for p in nd.lo:nd.hi for x in pts[perm[p]])
            (reduce((a, b) -> min.(a, b), ps), reduce((a, b) -> max.(a, b), ps))
        end
    end
    cplane, ccomp, convex, occluders = _hcluster_data(t, nodes, perm, pts)
    ball = zeros(Int, length(t.components))
    for (k, o) in enumerate(occluders)
        ball[o[1]] = k
    end
    cert, planar, hull = _hrim_data(t, ball)
    # per element: the convex part with balls it lies on, its facets in one plane (else 0)
    onconvex = [ball[t.elem_comp[e]] != 0 && onplane[e] ? Int(t.elem_comp[e]) : 0 for e in 1:nelements(t)]
    # per cluster: the largest extra lift of a ray the facing slack admits (`_vf_visible`)
    liftx = [maximum(p -> ctx.size[perm[p]] * ctx.slack[perm[p]], nd.lo:nd.hi) for nd in nodes]
    opts = (; tol = Float64(tol), maxrank = Int(maxrank), iperm, stats, nodes, allin, cplane, ccomp, convex, occluders,
            ball, cert, planar, hull, onconvex, liftx, qbox, rvis, hdev, bvh_comp = _bvh_components(ctx.bvh, t.tri_comp),
            facing_size = 16 * Int(leaf)^2, mask_size = 64 * Int(leaf)^2)
    return ctx, nodes, perm, far, near, opts
end

# The HMatrix of the assembled blocks.
function _hfinish(t, perm, opts, outs, tol, eta, leaf, maxrank)
    N, iperm, stats, nt = nelements(t), opts.iperm, opts.stats, length(outs)
    dense = sort!(reduce(vcat, getindex.(outs, 1)); by = x -> (x[1], x[2]))
    lowrank = sort!(reduce(vcat, getindex.(outs, 2)); by = x -> (x[1], x[2]))
    scattered = sort!(reduce(vcat, getindex.(outs, 3)); by = x -> (x[1][1], x[2][1], length(x[1]), length(x[2])))
    lowrank_index = sort!(reduce(vcat, getindex.(outs, 4)); by = x -> (x[1][1], x[2][1], length(x[1]), length(x[2])))
    groups = _hgroups(dense, lowrank, scattered, lowrank_index, nt)
    bytes = 4 * (sum(length(D) for (_, _, D) in dense; init = 0) +
                 sum(length(U) + length(V) for (_, _, U, V) in lowrank; init = 0)) +
            sum(4 * length(D) + 8 * (length(r) + length(c)) for (r, c, D) in scattered; init = 0) +
            sum(4 * (length(U) + length(V)) + 8 * (length(r) + length(c)) for (r, c, U, V) in lowrank_index; init = 0)
    ranks = [[size(U, 2) for (_, _, U, _) in lowrank]; [size(U, 2) for (_, _, U, _) in lowrank_index]]
    info = (; bytes, tol = Float64(tol), eta = Float64(eta), leaf = Int(leaf),
            maxrank = Int(maxrank), meanrank = isempty(ranks) ? 0.0 : sum(ranks) / length(ranks),
            map(x -> x[], stats)...,
            dense_entries = sum(length(D) for (_, _, D) in dense; init = 0) +
                            sum(length(D) for (_, _, D) in scattered; init = 0))
    return HMatrix(N, perm, iperm, copy(t.elem_area), ones(N), dense, lowrank, scattered, lowrank_index, groups, info)
end

# Assembly tasks taken one at a time (a function barrier: no captured variables).
function _hwork!(out, ctx, nodes, perm, tasks, next, opts)
    while (k = Threads.atomic_add!(next, 1)) <= length(tasks)
        (s, u), isfar = tasks[k]
        _hassemble!(out, ctx, nodes, perm, s, u, isfar, opts)
    end
    return nothing
end

# Blocks split into contiguous runs (by row) of about equal product cost.
function _hgroups(dense, lowrank, scattered, lowrank_index, nt)
    items = [[(0, k, dense[k][1], 2 * length(dense[k][3])) for k in eachindex(dense)];
             [(1, k, lowrank[k][1], 2 * (length(lowrank[k][3]) + length(lowrank[k][4])))
              for k in eachindex(lowrank)];
             [(2, k, scattered[k][1][1], 2 * length(scattered[k][3])) for k in eachindex(scattered)];
             [(3, k, lowrank_index[k][1][1], 2 * (length(lowrank_index[k][3]) + length(lowrank_index[k][4])))
              for k in eachindex(lowrank_index)]]
    sort!(items; by = x -> (x[3], x[1], x[2]))
    total = sum(x[4] for x in items; init = 0)
    ng = total < 200_000 ? 1 : nt
    groups = [Tuple{Int,Int}[] for _ in 1:ng]
    acc, g = 0, 1
    for x in items
        push!(groups[g], (x[1], x[2]))
        acc += x[4]
        g < ng && acc >= g * total / ng && (g += 1)
    end
    return groups
end

# ------------------------------------------------------------------ products

# Y += G Z (tree order) over the blocks of one group; Z, Y n × k.
function _hgroup_mul!(Y, H::HMatrix, Z, group, tmp)
    k = size(Z, 2)
    for (kind, idx) in group
        if kind == 3
            rs, cs, U, V = H.lowrank_index[idx]
            r = size(U, 2)
            for c in 1:k
                @inbounds for l in 1:r
                    sv, su = 0.0, 0.0
                    for b in eachindex(cs)
                        sv = muladd(Float64(V[b, l]), Z[cs[b], c], sv)
                    end
                    for a in eachindex(rs)
                        su = muladd(Float64(U[a, l]), Z[rs[a], c], su)
                    end
                    tmp[l] = sv
                    tmp[r+l] = su
                end
                @inbounds for l in 1:r
                    sv, su = tmp[l], tmp[r+l]
                    for a in eachindex(rs)
                        Y[rs[a], c] = muladd(Float64(U[a, l]), sv, Y[rs[a], c])
                    end
                    for b in eachindex(cs)
                        Y[cs[b], c] = muladd(Float64(V[b, l]), su, Y[cs[b], c])
                    end
                end
            end
        elseif kind == 2
            rs, cs, D = H.scattered[idx]
            for c in 1:k
                @inbounds for b in eachindex(cs)
                    zb = Z[cs[b], c]
                    acc = 0.0
                    for a in eachindex(rs)
                        d = Float64(D[a, b])
                        Y[rs[a], c] = muladd(d, zb, Y[rs[a], c])
                        acc = muladd(d, Z[rs[a], c], acc)
                    end
                    Y[cs[b], c] += acc
                end
            end
        elseif kind == 1
            r0, c0, U, V = H.lowrank[idx]
            m, n, r = size(U, 1), size(V, 1), size(U, 2)
            for c in 1:k
                # Y[rows] += U (Vᵀ Z[cols]),  Y[cols] += V (Uᵀ Z[rows])
                @inbounds for l in 1:r
                    sv, su = 0.0, 0.0
                    @simd for b in 1:n
                        sv = muladd(Float64(V[b, l]), Z[c0+b-1, c], sv)
                    end
                    @simd for a in 1:m
                        su = muladd(Float64(U[a, l]), Z[r0+a-1, c], su)
                    end
                    tmp[l] = sv
                    tmp[r+l] = su
                end
                @inbounds for l in 1:r
                    sv, su = tmp[l], tmp[r+l]
                    @simd for a in 1:m
                        Y[r0+a-1, c] = muladd(Float64(U[a, l]), sv, Y[r0+a-1, c])
                    end
                    @simd for b in 1:n
                        Y[c0+b-1, c] = muladd(Float64(V[b, l]), su, Y[c0+b-1, c])
                    end
                end
            end
        else
            r0, c0, D = H.dense[idx]
            m, n = size(D)
            diag = r0 == c0
            for c in 1:k
                if diag
                    @inbounds for b in 1:n
                        zb = Z[c0+b-1, c]
                        @simd for a in 1:m
                            Y[r0+a-1, c] = muladd(Float64(D[a, b]), zb, Y[r0+a-1, c])
                        end
                    end
                else
                    @inbounds for b in 1:n
                        zb = Z[c0+b-1, c]
                        s = 0.0
                        @simd for a in 1:m
                            d = Float64(D[a, b])
                            Y[r0+a-1, c] = muladd(d, zb, Y[r0+a-1, c])
                            s = muladd(d, Z[r0+a-1, c], s)
                        end
                        Y[c0+b-1, c] += s
                    end
                end
            end
        end
    end
    return Y
end

# W = G Z in tree order (n × k), threaded over the groups, each into its own buffer,
# summed in group order (deterministic for a given thread count).
function _hmul_tree(H::HMatrix, Z::Matrix{Float64})
    n, k = size(Z)
    ng = length(H.groups)
    maxr = max(maximum((size(x[3], 2) for x in H.lowrank); init = 0),
               maximum((size(x[3], 2) for x in H.lowrank_index); init = 0))
    W = zeros(n, k)
    ng == 0 && return W
    if ng == 1
        return _hgroup_mul!(W, H, Z, H.groups[1], zeros(2maxr))
    end
    bufs = [g == 1 ? W : zeros(n, k) for g in 1:ng]
    @sync for g in 1:ng
        Threads.@spawn _hgroup_mul!(bufs[g], H, Z, H.groups[g], zeros(2maxr))
    end
    for g in 2:ng
        W .+= bufs[g]
    end
    return W
end

# Y = F X = A⁻¹ S G S X (tr = false) or Fᵀ X = S G S A⁻¹ X (tr = true), in column chunks.
function _hmul!(Y::AbstractMatrix{Float64}, H::HMatrix, X::AbstractMatrix{Float64}, tr::Bool)
    n, k = H.n, size(X, 2)
    p, s, A = H.perm, H.scale, H.area
    for c0 in 1:16:max(k, 1)
        cc = c0:min(c0 + 15, k)
        isempty(cc) && break
        Z = Matrix{Float64}(undef, n, length(cc))
        for (q, c) in enumerate(cc), i in 1:n
            e = p[i]
            Z[i, q] = tr ? s[e] * X[e, c] / A[e] : s[e] * X[e, c]
        end
        W = _hmul_tree(H, Z)
        for (q, c) in enumerate(cc), i in 1:n
            e = p[i]
            Y[e, c] = tr ? s[e] * W[i, q] : s[e] * W[i, q] / A[e]
        end
    end
    return Y
end

# ------------------------------------------------- interpolated keyframe exchange

# The exchange factors between two keyframes of hierarchical models, Fᵢⱼ = c₀ᵢ F₀ᵢⱼ +
# c₁ᵢ F₁ᵢⱼ with c₀ = (1 - w) A₀/A and c₁ = w A₁/A (`_interpolate_exchange`), kept as the two
# operators (their blocks differ) with a symmetric balancing `S F S`.
struct _BlendedExchange{M0<:AbstractMatrix{Float64},M1<:AbstractMatrix{Float64}} <: AbstractMatrix{Float64}
    F0::M0
    F1::M1
    c0::Vector{Float64}
    c1::Vector{Float64}
    scale::Vector{Float64}
end

const _OperatorMatrix = Union{HMatrix,_BlendedExchange}

Base.size(B::_BlendedExchange) = size(B.F0)
Base.Matrix(B::_BlendedExchange) = _expand(B)
Base.Matrix{Float64}(B::_BlendedExchange) = _expand(B)
Base.getindex(B::_BlendedExchange, i::Int, j::Int) =
    B.scale[i] * (B.c0[i] * B.F0[i, j] + B.c1[i] * B.F1[i, j]) * B.scale[j]
Base.copy(B::_BlendedExchange) = _BlendedExchange(B.F0, B.F1, B.c0, B.c1, copy(B.scale))
_rowsums(B::_BlendedExchange) = B * ones(size(B, 1))
_balance!(B::_BlendedExchange, s::AbstractVector{Float64}) = (B.scale .*= s; B)

_interpolate_exchange(F0::HMatrix, F1::AbstractMatrix, A0, A1, A, w) =
    _BlendedExchange(F0, F1, (1 - w) .* A0 ./ A, w .* A1 ./ A, ones(length(A)))

# Y = F X (tr = false) or Fᵀ X
function _hmul!(Y::AbstractMatrix{Float64}, B::_BlendedExchange, X::AbstractMatrix{Float64}, tr::Bool)
    s = B.scale
    if tr
        Y .= s .* (transpose(B.F0) * (B.c0 .* s .* X) .+ transpose(B.F1) * (B.c1 .* s .* X))
    else
        SX = s .* X
        Y .= s .* (B.c0 .* (B.F0 * SX) .+ B.c1 .* (B.F1 * SX))
    end
    return Y
end

const _HT = Union{Transpose{Float64,<:_OperatorMatrix},Adjoint{Float64,<:_OperatorMatrix}}

function LinearAlgebra.mul!(y::StridedVector{Float64}, H::_OperatorMatrix, x::StridedVector{Float64},
                            α::Number, β::Number)
    length(x) == size(H, 1) == length(y) || throw(DimensionMismatch("matrix–vector sizes"))
    z = zeros(size(H, 1))
    _hmul!(reshape(z, :, 1), H, reshape(collect(x), :, 1), false)
    return _axpby!(y, z, α, β)
end
function LinearAlgebra.mul!(y::StridedVector{Float64}, Ht::_HT, x::StridedVector{Float64},
                            α::Number, β::Number)
    H = parent(Ht)
    length(x) == size(H, 1) == length(y) || throw(DimensionMismatch("matrix–vector sizes"))
    z = zeros(size(H, 1))
    _hmul!(reshape(z, :, 1), H, reshape(collect(x), :, 1), true)
    return _axpby!(y, z, α, β)
end
function LinearAlgebra.mul!(Y::StridedMatrix{Float64}, H::_OperatorMatrix, X::StridedMatrix{Float64},
                            α::Number, β::Number)
    size(X, 1) == size(H, 1) && size(Y) == (size(H, 1), size(X, 2)) || throw(DimensionMismatch("matrix sizes"))
    Z = similar(Y)
    _hmul!(Z, H, X, false)
    return _axpby!(Y, Z, α, β)
end
function LinearAlgebra.mul!(Y::StridedMatrix{Float64}, Ht::_HT, X::StridedMatrix{Float64},
                            α::Number, β::Number)
    H = parent(Ht)
    size(X, 1) == size(H, 1) && size(Y) == (size(H, 1), size(X, 2)) || throw(DimensionMismatch("matrix sizes"))
    Z = similar(Y)
    _hmul!(Z, H, X, true)
    return _axpby!(Y, Z, α, β)
end
# X F = (Fᵀ Xᵀ)ᵀ
function LinearAlgebra.mul!(Y::StridedMatrix{Float64}, X::StridedMatrix{Float64}, H::_OperatorMatrix,
                            α::Number, β::Number)
    size(X, 2) == size(H, 1) && size(Y) == (size(X, 1), size(H, 1)) || throw(DimensionMismatch("matrix sizes"))
    Z = Matrix{Float64}(undef, size(H, 1), size(X, 1))
    _hmul!(Z, H, Matrix(transpose(X)), true)
    return _axpby!(Y, transpose(Z), α, β)
end

# Products with vectors of first-order Float64 duals (values and partials as columns).
for T in (:_OperatorMatrix, :_HT)
    @eval function Base.:*(A::$T, x::AbstractVector{ForwardDiff.Dual{Tag,Float64,N}}) where {Tag,N}
        X = Matrix{Float64}(undef, length(x), N + 1)
        for (i, d) in enumerate(x)
            X[i, 1] = ForwardDiff.value(d)
            p = ForwardDiff.partials(d)
            for k in 1:N
                X[i, k+1] = p[k]
            end
        end
        Y = A * X
        return [ForwardDiff.Dual{Tag}(Y[i, 1], ntuple(k -> Y[i, k+1], Val(N))...) for i in axes(Y, 1)]
    end
end

