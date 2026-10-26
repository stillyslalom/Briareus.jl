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
  circumscribed ball) decides each pair exactly by the part's balls: a segment missing
  the circumscribed ball is visible, one through the inscribed ball hidden, and only the
  few between take a ray. The mask makes the block zero, low rank without further rays,
  or refined with it down to dense leaf pairs (a shadow edge);
- any other occupied shaft is refined down to leaf pairs, which take every pair's
  visibility ray: wholly hidden pairs give zero, wholly visible ones a cross
  approximation, and partly seen ones (a rim, `:sampled` fractions) stay dense;
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
zero blocks and cut rows, "H-matrix near field").
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

# What occludes the shaft between clusters a and b: 0 if nothing can occlude a segment
# between them, c > 0 if only facets of the convex component c can and its balls certify
# each segment (`_hball_certifies`: neither cluster on it, both outside its circumscribed
# ball), else -1. A facet is not an occluder (for the pairs `_vf_pair_free` evaluates
# without a ray: facing beyond the slack, so leaving their facets) if
# - it lies on a curved surface (cylinder, sphere, cone) both clusters lie inside;
# - it belongs to a convex component (every vertex behind every facet) one of the
#   clusters lies on entirely, with one facet per element: a segment leaving a face of a
#   convex polyhedron outward stays outside it;
# - each cluster lies strictly in front of its plane, or is planar in that plane facing
#   the same way: the lifted ray ends are in front, so the segment does not cross it.
function _shaft_occluder(ctx, opts, sa::Int, sb::Int)
    (; t, bvh, tri_surf, ε) = ctx
    a, b = opts.nodes[sa], opts.nodes[sb]
    allin, cplane, ccomp, convex = opts.allin, opts.cplane, opts.ccomp, opts.convex
    planes = _shaft_planes(a.bmin, a.bmax, b.bmin, b.bmax)
    outside_box(lo, hi) = any(((n, d),) -> dot(n, ifelse.(n .> 0, lo, hi)) > d, planes)
    ca, cb = _box_corners(a.cmin, a.cmax), _box_corners(b.cmin, b.cmax)
    L = norm(max.(a.cmax, b.cmax) - min.(a.cmin, b.cmin))
    function clear(c, corners, n, d)
        all(x -> dot(n, x) - d > 2ε + 1e-12 * L, corners) && return true
        for (pn, pd, off) in cplane[c]
            dot(pn, n) > 1 - 1e-12 && abs(pd - d) <= 1e-12 * (L + abs(d)) && off > 2ε + 1e-12 * L && return true
        end
        return false
    end
    # components whose facets need no test: convex ones a cluster lies on (never
    # occluders), and the certified one once a facet of it is in the shaft
    excluded(c) = c > 0 && convex[c] && (ccomp[sa] == c || ccomp[sb] == c)
    isempty(bvh.nodes) && return 0
    found = 0
    stack = Int32[1]
    while !isempty(stack)
        ni = pop!(stack)
        node = bvh.nodes[ni]
        sig = opts.bvh_comp[ni]
        (excluded(sig) || (found > 0 && sig == found)) && continue
        outside_box(node.lo, node.hi) && continue
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
                excluded(c) && continue
                n = t.tri_normal[id]
                d = dot(n, v0)
                clear(sa, ca, n, d) && clear(sb, cb, n, d) && continue
                c == found && continue
                found == 0 && _hball_certifies(opts, c, a, b) || return -1
                found = Int(c)
            end
        else
            push!(stack, node.a, node.b)
        end
    end
    return found
end
_shaft_occupied(ctx, opts, sa::Int, sb::Int) = _shaft_occluder(ctx, opts, sa, sb) != 0

# Whether the balls of convex component c decide every pair between clusters a and b
# (`_hvis_ball`): c has balls (`occluders`), neither cluster lies on c, and both boxes
# lie outside its circumscribed ball.
function _hball_certifies(opts, c, a::_HCluster, b::_HCluster)
    k = opts.ball[c]
    k == 0 && return false
    _, ctr, _, rout = opts.occluders[k]
    return _box_distance(ctr, a.bmin, a.bmax) > rout && _box_distance(ctr, b.bmin, b.bmax) > rout
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
function _hcluster_data(t::Target, nodes, perm)
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
        # cluster's other elements (Inf if none)
        planes = Tuple{Vec3,Float64,Float64}[]
        for e in es
            flat[e] || continue
            n = t.elem_normal[e]
            d = dot(n, t.elem_centroid[e])
            any(((pn, pd, _),) -> dot(pn, n) > 1 - 1e-12 && abs(pd - d) <= 1e-12 * L, planes) && continue
            off = Inf
            for f in es
                on = flat[f] && dot(t.elem_normal[f], n) > 1 - 1e-12 && abs(dot(n, t.elem_centroid[f]) - d) <= 1e-12 * L
                on || (off = min(off, dot(n, t.elem_centroid[f]) - d))
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

# G at tree positions (p, q); `free`: the pair's shaft holds no occluder (no ray needed),
# or the block's visibility mask.
function _hentry_eval(ctx, perm, p, q, free::Bool)
    i, j = perm[p], perm[q]
    (!free || i == j) && return _vf_exchange(ctx, i, j)
    return i < j ? _vf_pair_free(ctx, i, j) : _vf_pair_free(ctx, j, i)
end
function _hentry_eval(ctx, perm, p, q, M::_HMask)
    M.v[p-M.r0+1, q-M.c0+1] > 0 || return 0.0
    i, j = minmax(perm[p], perm[q])
    d = ctx.C[j] - ctx.C[i]
    return _vf_kernel(ctx, i, j, 1.0, d, dot(d, d))
end

# Visibility of the pair i < j whose shaft holds facets of one convex component only, its
# inscribed ball (center `ctr`, radius `rin`) and circumscribed ball (`rout`) both
# clusters' boxes lie outside of (`_hball_certifies`): a segment missing the
# circumscribed ball hits no facet of it, and one through the inscribed ball crosses the
# closed surface. Only pairs between the two, or admitted by the facing slack (whose rays
# start farther off their facets), take `_vf_visible`'s ray. Centroid visibility only.
function _hvis_ball(ctx, i, j, ctr, rin, rout)
    (; C, Nrm, slack, ε) = ctx
    d = C[j] - C[i]
    r2 = dot(d, d)
    r2 == 0 && return Int8(-1)
    r = sqrt(r2)
    ci, cj = dot(Nrm[i], d) / r, -dot(Nrm[j], d) / r
    ci > -slack[i] && cj > -slack[j] || return Int8(-1)
    if ci > slack[i] && cj > slack[j]
        # the ray ends of `_vf_visible` (lifted by ε)
        oi, oj = C[i] + ε * Nrm[i], C[j] + ε * Nrm[j]
        w = oj - oi
        u = clamp(dot(ctr - oi, w) / dot(w, w), 0.0, 1.0)
        x = oi + u * w - ctr
        q2 = dot(x, x)
        q2 > (rout * (1 + 1e-9))^2 && return Int8(1)
        q2 < (rin * (1 - 1e-9))^2 && return Int8(0)
    end
    return first(_vf_visible(ctx, i, j)) > 0 ? Int8(1) : Int8(0)
end

# The visibility mask of the block rows × cols whose shaft component c certifies.
function _hmask(ctx, opts, perm, rows::UnitRange{Int}, cols::UnitRange{Int}, c::Int)
    _, ctr, rin, rout = opts.occluders[opts.ball[c]]
    V = Matrix{Int8}(undef, length(rows), length(cols))
    for (b, q) in enumerate(cols), (a, p) in enumerate(rows)
        V[a, b] = _hvis_ball(ctx, minmax(perm[p], perm[q])..., ctr, rin, rout)
    end
    return _HMask(V, first(rows), first(cols))
end

# `_vf_pair` of a pair no facet can occlude: a pair facing beyond the slack of the
# facing test leaves its elements' facets, and its ray would be unobstructed; the others
# (admitted by the slack) take the ray.
function _vf_pair_free(ctx, i, j)
    (; C, Nrm, slack) = ctx
    d = C[j] - C[i]
    r2 = dot(d, d)
    r2 == 0 && return 0.0
    r = sqrt(r2)
    ci, cj = dot(Nrm[i], d) / r, -dot(Nrm[j], d) / r
    ci > -slack[i] && cj > -slack[j] || return 0.0
    ci > slack[i] && cj > slack[j] || return _vf_pair(ctx, i, j)
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

# Estimated Frobenius error of U Vᵀ against the block, from `ns` sampled entries.
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
    return sqrt(e2 * m * n / ns), sqrt(g2 * m * n / ns)
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
                ev = known !== nothing ? known : !_shaft_occupied(ctx, opts, s, t)
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
            occ = _shaft_occluder(ctx, opts, s, t)
            if occ != 0 && _convex_hides(opts, s, t)
                Threads.atomic_add!(st.shadowed, 1)
                return                                   # wholly behind a convex part
            end
            if occ > 0 && ctx.vis === nothing && m * n <= opts.mask_size
                # only a convex part's facets in the shaft: its balls decide most pairs
                # (centroid visibility; `:sampled` fractions take the rays below)
                ev = _hmask(ctx, opts, perm, rows, cols, occ)
                Threads.atomic_add!(st.masked, 1)
            elseif occ != 0 && leaves
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
                ev = occ == 0
            end
        end
        if ev isa _HMask
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
                            elseif _hleaf(b) || length(pa) <= 4 || !(ev === true || ev isa _HMask)
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
        else
            D = _hdense(ctx, perm, rows, cols, far ? ev : !_shaft_occupied(ctx, opts, s, t))
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
    N = nelements(t)
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
    cplane, ccomp, convex, occluders = _hcluster_data(t, nodes, perm)
    ball = zeros(Int, length(t.components))
    for (k, o) in enumerate(occluders)
        ball[o[1]] = k
    end
    opts = (; tol = Float64(tol), maxrank = Int(maxrank), iperm, stats, nodes, allin, cplane, ccomp, convex, occluders,
            ball, bvh_comp = _bvh_components(ctx.bvh, t.tri_comp), facing_size = 16 * Int(leaf)^2,
            mask_size = 64 * Int(leaf)^2)
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

