# This file computes hierarchical (H-matrix) view factors. Near pairs are dense, and far
# cluster pairs of uniform visibility are low rank, by adaptive cross approximation on
# the exact pair kernel. Storage and assembly then grow as N^1.3–1.5 on a hohlraum,
# against N² for a dense matrix; the shadow edges of the capsule and the near field stay
# dense. Only products are offered (`mul!` with F and Fᵀ), which is all that the
# iterative radiosity path uses.

"""
    HMatrix

View factors in hierarchical form, from `view_factors(target; storage = :hmatrix)` and
by default from `RadiationModel` above `Briareus.HMATRIX_MODEL[]` elements. The form
stores the symmetric exchange matrix `G = diag(A) F` (`Gᵢⱼ = Aᵢ Fᵢⱼ`), once per
unordered pair so that reciprocity is exact, partitioned over a cluster tree of the
elements. Then `F = diag(A)⁻¹ S G S` with a diagonal balancing `S`, and
`enforce_closure!` updates only `S`. The type is an `AbstractMatrix{Float64}`.
Products with Float64 vectors and matrices (`*`, `mul!`, also transposed, and with dual
vectors) accumulate in Float64, which is all that the iterative radiosity path uses.
Indexing one entry searches the blocks, for inspection only, and `Matrix(H)` expands
the matrix by products with the identity.

Clusters split the elements by the surface they lie on (the curved surface of a
component, or the component). They then bisect the centroids along the longest axis of
their box down to `leaf` elements. No cluster therefore mixes capsule and wall, or cap
and wall at a corner; such blocks seldom compress, and their low-rank error falls on
the small rows. A pair of clusters whose boxes are `eta` × the larger diameter apart is
a candidate for low rank, decided block by block:

- Rows and columns whose elements face no centroid of the other cluster are zero, by a
  box test. A block that a convex closed part, such as a capsule, hides wholly is also
  zero. The test is that the inscribed ball blocks every segment between the boxes, by
  a cone test from each corner.
- The *shaft* between the boxes, their convex hull, is searched for facets that could
  occlude a segment between them (the shaft culling of Haines & Wallace). The search
  excludes the facets of a curved surface that both clusters lie inside, of a convex
  part that one of them lies on, and of a plane that both lie in front of or in. An
  empty shaft makes every facing pair wholly visible, so its entries need no ray, and
  the block is cross-approximated. The approximation is ACA with partial pivoting on the
  exact pair kernel, recompressed by SVD to `tol` relative Frobenius error and checked
  on sampled entries. A block that does not compress within `maxrank` terms is refined.
- A shaft with facets of only one convex part, such as the capsule with both clusters
  outside its circumscribed ball, resolves each pair exactly by the balls of the part. A
  segment that misses the circumscribed ball is visible, one through the inscribed ball
  is hidden, and only the few between take a ray. With the resulting mask the block is
  zero, low rank without further rays, or refined with the mask down to dense leaf pairs
  at a shadow edge.
- Any other occupied shaft is refined down to leaf pairs. A leaf pair whose shaft
  contains facets of up to 8 parts with a certificate is decided in the same way, pair
  by pair and exactly as the ray would decide it. The certificates are:
  - for planar parts such as washers, LEH caps and discs, the point where the segment
    crosses the plane of the part, inside or outside its facets with margin;
  - for a cylinder or cone part, its axisymmetric bounds. A prism of facets parallel to
    its axis, such as a wall of chord rectangles, is also decided by the side lines of
    its cross-section. A segment qualifies when its part within the axial range of the
    prism stays inside the half-plane of every side, as a segment through a waist gap
    onto the wall does;
  - for a convex part, its balls or the plane of an end element lying on it.

  Other leaf pairs take the visibility ray of every pair. Wholly hidden pairs give zero,
  wholly visible ones a cross approximation, and partly visible ones, such as at a rim
  or with `:sampled` fractions, stay dense. With `visibility = :sampled` the
  certificates resolve the whole bundle of sample rays of a pair at once, with margins
  grown by the spread of the samples, and otherwise each of its rays. The fractions are
  those of `view_factors`.
- Pairs that the facing slack admits, whose rays start farther off curved facets, are
  judged on their ray ends when the exclusions of the shaft hold with that margin, and
  take the ray otherwise.
- Rows cut by the facing test, at the terminator of a sphere, are refined against the
  children of the column cluster beside the low-rank rest. They are low rank on the rows
  that face a child wholly, and dense against a leaf or when few. A small block with
  only a few facing pairs is stored on those pairs.

Near pairs are dense (Float32), with the same entries as `view_factors` (`_vf_pair`).
Blocks whose entries are all zero are not stored.

On the hohlraum presets (10–29k elements) the default `tol = 1e-6` changes the capsule
moments by ≤ 1e-8 and the incident flux by ≤ 2e-6 relative per element against the
Float32 matrix or a tighter `HMatrix` (research notes, "Hierarchical view factors" and,
for the surface-split cluster tree and `tol`, "The H-matrix's N^1.5 remainder"; for the
zero blocks and cut rows, "H-matrix near field"; for the visibility certificates,
"H-matrix visibility certificates").
"""
struct HMatrix <: AbstractMatrix{Float64}
    n::Int
    perm::Vector{Int}                   # element at each tree position
    iperm::Vector{Int}                  # tree position of each element
    area::Vector{Float64}               # element areas (element order)
    scale::Vector{Float64}              # balancing S (element order)
    dense::Vector{Tuple{Int,Int,Matrix{Float32}}}               # (row0, col0, D), tree order
    lowrank::Vector{Tuple{Int,Int,Matrix{Float32},Matrix{Float32}}}   # (row0, col0, U, V)
    scattered::Vector{Tuple{Vector{Int},Vector{Int},Matrix{Float32}}}  # (rows, cols, D), off the diagonal
    lowrank_index::Vector{Tuple{Vector{Int},Vector{Int},Matrix{Float32},Matrix{Float32}}}  # (rows, cols, U, V)
    groups::Vector{Vector{Tuple{Int,Int}}}    # product work per task, as (kind, index)
    info::NamedTuple
end

Base.size(H::HMatrix) = (H.n, H.n)
Base.copy(H::HMatrix) = HMatrix(H.n, H.perm, H.iperm, H.area, copy(H.scale), H.dense, H.lowrank, H.scattered,
                                H.lowrank_index, H.groups, H.info)

# Returns the stored G entry at tree positions (p, q).
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
    # The indexed blocks come before the contiguous ones, since the cut rows of a low-rank
    # block are zero in it and stored in scattered or indexed blocks.
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

# The exact visibility of every pair of a block, the verdict of `_vf_visible` by tree
# positions from (r0, c0): -1 for not facing, 0 for hidden and 1 for visible.
struct _HMask
    v::Matrix{Int8}
    r0::Int
    c0::Int
end

# The visible fractions of every pair of a block, those of `_vf_visible` with -1 for not
# facing, by tree positions from (r0, c0) (`:sampled` visibility).
struct _HFrac
    v::Matrix{Float64}
    r0::Int
    c0::Int
end

# Returns G at tree positions (p, q). `ev` is an empty `_HOcc` when the shaft of the pair
# holds no occluder, with no ray except for some pairs that the facing slack admits. It is
# `false` for a ray per pair, or the visibility mask of the block.
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

# Returns the ray ends of `_vf_visible` for the pair i < j (centroids at distance r,
# facing cosines ci, cj) that only the facing slack admits. It also returns whether the
# exclusions of its shaft (`_shaft_occluder`) hold for that segment. They hold when:
# - the ends stay inside the boxes of the clusters, lifted by at most the padding of
#   the boxes;
# - the plane exclusions had the margin of the lifts (`occ.slack`);
# - each convex component on which a cluster lies is cleared by the plane of the element
#   of an end on it, with both ends in front by `ctx.hsep`. `ctx.onplane` marks elements
#   with their facets in that plane, behind which the convex part lies.
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

# Returns the visibility of the pair i < j whose shaft contains only facets of the
# certified components `occ.c[1:occ.n]` (`_shaft_occluder`), as the ray of `_vf_visible`
# finds it. The pair is hidden if one certificate finds the segment through its facets,
# and visible if every certificate finds it clear (`_hcert_verdict`). The certificates
# are the balls of a convex part, the facets of a planar part and the axisymmetric bounds
# of a part. A pair that the facing slack admits is judged on its ray ends
# (`_hslack_ends`). Only pairs left undecided take the ray. For centroid visibility only.
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

# Returns the visibility mask of the block rows × cols whose shaft `occ` certifies.
function _hmask(ctx, opts, perm, rows::UnitRange{Int}, cols::UnitRange{Int}, occ::_HOcc)
    V = Matrix{Int8}(undef, length(rows), length(cols))
    for (b, q) in enumerate(cols), (a, p) in enumerate(rows)
        V[a, b] = _hvis_cert(ctx, opts, minmax(perm[p], perm[q])..., occ)
    end
    return _HMask(V, first(rows), first(cols))
end

# Returns the `:sampled` visible fraction of the pair i < j whose shaft contains only
# facets of the certified components `occ.c[1:occ.n]` (`_shaft_occluder`, searched with
# the sample boxes of the clusters). The result is `===` that of `_vf_visible`. The rays
# of the pair, which are the centroid ray, the rays of the probes and those between the
# samples, all run between points within `opts.rvis` of the two ray ends. Each
# certificate therefore decides the bundle (`_hcert_verdict` with radii). The pair is
# hidden if one certificate finds every segment through its facets, and visible if all
# find every segment clear. Otherwise `_visible_fraction` runs as `_vf_visible` does,
# with each segment decided by the certificates that did not clear the bundle (`_HSeen`)
# or by its ray. The exclusions of the shaft must hold for the samples:
# - their ends lie inside the cluster boxes (rvis plus the lift within the padding);
# - for a pair that the slack admits, the plane exclusions hold with the margin of the
#   slack lifts (`occ.slack`);
# - each convex part on which a cluster lies is cleared by the plane of an end element,
#   per bundle or else per segment.
# Folded elements, whose facets the rays of the samples skip, take `_vf_visible`.
function _hvis_bundle(ctx, opts, i, j, occ::_HOcc)
    (; t, bvh, C, Nrm, size, slack, ε, tri_surf, inside, vis, multi, onplane, hsep) = ctx
    facing, _, _, ci, cj, oi, oj = _vf_ends(C, Nrm, size, slack, ε, i, j)    # as in `_vf_visible`
    facing || return -1.0
    (multi[i] || multi[j]) && return first(_vf_visible(ctx, i, j))
    li, lj = oi - C[i], oj - C[j]
    ri, rj = opts.rvis[i], opts.rvis[j]
    ok = ((ci > slack[i] && cj > slack[j]) || occ.slack) &&
         ri + norm(li) <= size[i] / sqrt(2) && rj + norm(lj) <= size[j] / sqrt(2)
    ok || return first(_vf_visible(ctx, i, j))
    # convex parts on which a cluster lies are cleared for the bundle, or else checked per
    # segment
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

# Whether the segment from p to q between the ray samples of elements i and j is
# unobstructed. It is decided by the certificates of `occ` that did not clear the bundle
# (`vb[k] != 1`), once the convex parts `xs` on which a cluster lies are cleared by the
# plane of an end element, and otherwise by the ray.
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

# Returns the visible fractions of the block rows × cols whose shaft `occ` certifies.
function _hfrac(ctx, opts, perm, rows::UnitRange{Int}, cols::UnitRange{Int}, occ::_HOcc)
    V = Matrix{Float64}(undef, length(rows), length(cols))
    for (b, q) in enumerate(cols), (a, p) in enumerate(rows)
        V[a, b] = _hvis_bundle(ctx, opts, minmax(perm[p], perm[q])..., occ)
    end
    return _HFrac(V, first(rows), first(cols))
end

# Returns the certified visibility of a block, as a mask (centroid visibility) or as
# fractions.
_hcertified(ctx, opts, perm, rows, cols, occ) =
    ctx.vis === nothing ? _hmask(ctx, opts, perm, rows, cols, occ) : _hfrac(ctx, opts, perm, rows, cols, occ)

# `_vf_pair` of a pair i < j whose shaft `occ` holds no occluder. A pair facing beyond the
# slack of the facing test leaves the facets of its elements, and its ray would be
# unobstructed. A pair admitted by the slack is visible if its ray ends satisfy the
# exclusions of the shaft (`_hslack_ends`, centroid visibility), and otherwise it takes
# the ray.
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

# Returns the dense block of G at tree positions rows × cols (Float32). A diagonal block
# is filled from its upper triangle.
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

# Returns the dense block at tree positions rows × cols (index lists of two different
# clusters).
function _hdense_index(ctx, perm, rows::Vector{Int}, cols::Vector{Int}, free)
    D = Matrix{Float32}(undef, length(rows), length(cols))
    for b in eachindex(cols), a in eachindex(rows)
        D[a, b] = Float32(_hentry_eval(ctx, perm, rows[a], cols[b], free))
    end
    return D
end

# Returns the dense block rows × cols (s ≠ t) whose visibility `vis` is known, as the
# fraction of `_vf_visible` by position in the block.
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

# Returns the positions in `rows` whose element may face some centroid of cluster c,
# which are those with a corner of the centroid box of c in front of the plane of the
# element, by the slack of the facing test. The rows of the other positions in the block
# are zero. It also returns the positions that face all of the box of c beyond the slack.
# Their rows have no facing cut-off, a kink in the kernel that the cross approximation
# resolves slowly at the terminator of a sphere.
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

# Returns the rows and columns (of `rows` × `cols`) with some pair that passes the facing
# test.
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
            # after a zero residual row, the block may be zero (rank 0), or only this row
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
        # The approximation stops early when the extrapolated decay of the terms misses tol
        # at kmax. A block that will not compress is refined, and its children then pay
        # for it again.
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

# Recompresses U Vᵀ (QR of both factors, SVD of the small core) to relative accuracy tol.
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

# Returns the estimated Frobenius error of U Vᵀ against the block, from `ns` sampled
# entries, plus the error of the entries where the quadrature tier of a row or a column
# switches (`_htier_entries`). G jumps there by the difference of the quadrature rules.
# An isolated jump, such as one pair of a polar sliver of a capsule past the far
# boundary, is missed by the pivots of the cross approximation and the random samples
# alike.
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

# Returns the positions (a, b) of the block rows × cols where G may jump. In each row
# (column) whose pairs fall in more than one quadrature tier (`_vf_tier`), these are the
# pairs outside its most common tier, at most m + n of them. There are none, without a
# pass over the pairs, if the centroid boxes are farther apart than `mid` times the
# largest element, so that all pairs are far.
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
    # Lines (a row a > 0 or a column -a) with their pairs outside the most common tier,
    # the fewest first, since the cross approximation misses an isolated jump.
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

# Returns the exact visibility of the block rows × cols by a ray per pair (`_vf_visible`)
# on the facing positions ra × cb, since the others do not face. For centroid visibility
# only.
function _hmask_rays(ctx, perm, rows::UnitRange{Int}, cols::UnitRange{Int}, ra, cb)
    V = fill(Int8(-1), length(rows), length(cols))
    for q in cb, p in ra
        f = first(_vf_visible(ctx, minmax(perm[p], perm[q])...))
        V[p-first(rows)+1, q-first(cols)+1] = f < 0 ? Int8(-1) : f > 0 ? Int8(1) : Int8(0)
    end
    return _HMask(V, first(rows), first(cols))
end

# Returns the visible and hidden facing pairs of a mask on positions ra × cb.
function _hmask_counts(M::_HMask, ra, cb)
    nv = nh = 0
    for q in cb, p in ra
        f = M.v[p-M.r0+1, q-M.c0+1]
        nv += f > 0
        nh += f == 0
    end
    return nv, nh
end
# (for fractions, the visible pairs, and the hidden or partly visible pairs)
function _hmask_counts(M::_HFrac, ra, cb)
    nv = nh = 0
    for q in cb, p in ra
        f = M.v[p-M.r0+1, q-M.c0+1]
        nv += f > 0
        nh += 0 <= f < 1
    end
    return nv, nh
end

# Assembles the blocks of cluster pair (s, t) into `dense`/`lowrank`. The block is low
# rank when the pair is admissible and of uniform visibility, and its cross approximation
# converges and passes a sampled check. Otherwise its children are assembled, down to
# dense leaf pairs. `known` is what a parent block found of the visibility of the pair
# (`true` for an empty shaft, or a mask), or else `nothing`. `sub` holds the rows of s
# that a facing test cut against a parent of t (`split` below), refined only over t.
# Against the smaller boxes of t, most of them face wholly (low rank) or not at all.
# Blocks go into `out` = (dense, lowrank, scattered, lowrank_index).
function _hassemble!(out, ctx, nodes, perm, s, t, far, opts, known = nothing,
                     sub::Union{Nothing,Vector{Int}} = nothing)
    dense, lowrank, scattered, lowrank_index = out
    a, b = nodes[s], nodes[t]
    rows, cols = sub === nothing ? (a.lo:a.hi) : sub, b.lo:b.hi
    m, n = length(rows), length(cols)
    st = opts.stats
    leaves = _hleaf(a) && _hleaf(b)
    vis = nothing        # visible fractions of a leaf pair (`:sampled` visibility)
    ev = false           # entries from `true` (empty shaft, no rays), `false` (rays) or a mask
    ra, cb = Int[], Int[]
    if far
        Threads.atomic_add!(st.candidates, 1)
        # rows (columns) whose element faces no centroid of the other cluster are zero
        amin, amax = sub === nothing ? (a.cmin, a.cmax) : _hcentroid_box(ctx.C, perm, sub)
        (ra, fa), (cb, fb) = _hactive(ctx, perm, rows, b), _hactive(ctx, perm, cols, amin, amax)
        if !isempty(ra) && !isempty(cb) && length(ra) * length(cb) <= opts.facing_size
            # Small blocks get the facing test pair by pair (no rays), which leaves the rows
            # and columns of some facing pair. A block facing itself in only a few pairs,
            # at the terminator of a sphere, is kept dense on those pairs.
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
                # With only facets of certified parts in the shaft, their certificates
                # decide most pairs (with `:sampled`, whole sample bundles) at any level. A
                # block of uniform visibility then gets the cross approximation and is not
                # refined.
                ev = _hcertified(ctx, opts, perm, rows, cols, occ)
                Threads.atomic_add!(st.masked, 1)
            elseif occ.n != 0 && leaves
                # the exact visibility of every pair, with a ray each
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
            # The block is zero if wholly hidden, low rank if wholly visible (no more rays),
            # and otherwise refined with the mask down to dense leaf pairs.
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
            # Rows and columns that face the other cluster wholly get the cross
            # approximation, and those cut by the facing test, if few, are kept dense.
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
                                # the cut rows against the children of t, where most face
                                # wholly or not at all
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
        # Cut rows are dense against a leaf or if they are few, since low rank would not
        # pay, and are otherwise refined over t.
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
            # A near leaf pair needs no ray in an empty shaft, uses the mask of the
            # certificates in a certified one, and otherwise takes a ray per pair.
            occ = _shaft_occluder(ctx, opts, s, t, true)
            D = _hdense(ctx, perm, rows, cols, occ.n == 0 ? occ :
                        occ.n > 0 ? _hcertified(ctx, opts, perm, rows, cols, occ) : false)
        end
        any(!iszero, D) && push!(dense, (a.lo, b.lo, D))     # (zero in the near field of a capsule)
    else
        # the visibility found for this block also holds for its children
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

# Pushes a scattered block, unless all its entries are zero, as for rows cut at the
# terminator of a sphere against the far side of a box.
_hpush_nonzero!(v, blk) = (any(!iszero, blk[3]) && push!(v, blk); v)

# Returns the bounding box of the centroids at tree positions ps.
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

- `tol`: the relative Frobenius accuracy of each low-rank block;
- `eta`: the admissibility, with cluster boxes at least `eta` × the larger diameter
  apart;
- `leaf`: the largest leaf cluster;
- `maxrank`: the most terms of a cross approximation before the block is refined.

These are engineering defaults, from the work–precision sweeps in the research notes
("Hierarchical view factors"). The change of `tol` from 1e-5 to 1e-6 is in "The
H-matrix's N^1.5 remainder".
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

# Returns the assembly context, cluster tree, block pairs and assembly data (visibility
# certificates) of a view-factor context. The assembly context is the view-factor
# context with sector receivers vectorized and the data of `_hslack_ends`.
function _hsetup(ctx, tol, eta, leaf, maxrank)
    t = ctx.t
    # per element, whether its facets lie in one plane through its centroid with its normal
    onplane = map(1:nelements(t)) do e
        n, C = t.elem_normal[e], t.elem_centroid[e]
        δ = 1e-12 * (norm(C) + 1)
        all(k -> dot(t.tri_normal[k], n) > 1 - 1e-12 && abs(dot(n, t.tri_v0[k] - C)) <= δ, t.elem_tris[e])
    end
    # Sectors receive through the vectorized contour (`_SoAReceiver`, within round-off).
    # `hsep` is the margin of a plane that clears a convex part: 10× its convexity
    # tolerance, 1e-9 of its size, and a tenth of the lift ε of the rays.
    ctx = merge(ctx, (; recv = Union{_Receiver,_SoAReceiver,_PolyReceiver}[_soa_receiver(r) for r in ctx.recv],
                      onplane = BitVector(onplane), hsep = 1e-8 * max(t.scale, 1.0)))
    pad = ctx.size ./ sqrt(2)
    nodes, perm = _hcluster_tree(ctx.C, pad, Int(leaf), _hsurface_groups(t))
    iperm = invperm(perm)
    far, near = _hblock_pairs(nodes, eta)
    stats = NamedTuple{(:candidates, :lowrank, :facing, :sparse, :shadowed, :masked, :hidden, :mixed, :unconverged,
                        :rejected)}(ntuple(_ -> Threads.Atomic{Int}(0), 10))
    # allin[k, c] marks every element of cluster c inside curved surface k
    ns = size(ctx.inside, 1)
    allin = BitMatrix(undef, ns, length(nodes))
    for (c, nd) in enumerate(nodes), k in 1:ns
        allin[k, c] = all(p -> ctx.inside[k, perm[p]], nd.lo:nd.hi)
    end
    # For `:sampled` visibility, per element: the ends of its sample rays lifted by ε, as
    # `_vf_visible` lifts them for a pair facing beyond the slack, their largest distance
    # from the lifted centroid (`rvis`), and the largest distance of the samples from the
    # plane of the element (`hdev`). Per cluster: the box of those ends, for the plane
    # tests of the shaft (`qbox`).
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
    # per element, the convex part with balls on which it lies, with its facets in one
    # plane (else 0)
    onconvex = [ball[t.elem_comp[e]] != 0 && onplane[e] ? Int(t.elem_comp[e]) : 0 for e in 1:nelements(t)]
    # per cluster, the largest extra lift of a ray that the facing slack admits (`_vf_visible`)
    liftx = [maximum(p -> ctx.size[perm[p]] * ctx.slack[perm[p]], nd.lo:nd.hi) for nd in nodes]
    opts = (; tol = Float64(tol), maxrank = Int(maxrank), iperm, stats, nodes, allin, cplane, ccomp, convex, occluders,
            ball, cert, planar, hull, onconvex, liftx, qbox, rvis, hdev, bvh_comp = _bvh_components(ctx.bvh, t.tri_comp),
            facing_size = 16 * Int(leaf)^2, mask_size = 64 * Int(leaf)^2)
    return ctx, nodes, perm, far, near, opts
end

# Returns the HMatrix of the assembled blocks.
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

# Takes assembly tasks one at a time, behind a function barrier without captured
# variables.
function _hwork!(out, ctx, nodes, perm, tasks, next, opts)
    while (k = Threads.atomic_add!(next, 1)) <= length(tasks)
        (s, u), isfar = tasks[k]
        _hassemble!(out, ctx, nodes, perm, s, u, isfar, opts)
    end
    return nothing
end

# Splits the blocks into contiguous runs (by row) of about equal product cost.
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

# Computes Y += G Z (tree order) over the blocks of one group, with Z, Y n × k.
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

# Computes W = G Z in tree order (n × k), threaded over the groups. Each group writes a
# separate buffer, and the buffers are summed in group order, which is deterministic for
# a given thread count.
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

# Computes Y = F X = A⁻¹ S G S X (tr = false) or Fᵀ X = S G S A⁻¹ X (tr = true), in
# column chunks.
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

# The exchange factors between two keyframes of hierarchical models,
# Fᵢⱼ = c₀ᵢ F₀ᵢⱼ + c₁ᵢ F₁ᵢⱼ with c₀ = (1 - w) A₀/A and c₁ = w A₁/A (`_interpolate_exchange`).
# They are kept as the two operators, since their blocks differ, with a symmetric
# balancing `S F S`.
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

# Products with strided vectors of first-order Float64 duals, with values and partials as
# columns. The vectors are strided, so that FillArrays zero vectors keep their separate
# methods.
for T in (:_OperatorMatrix, :_HT)
    @eval function Base.:*(A::$T, x::StridedVector{ForwardDiff.Dual{Tag,Float64,N}}) where {Tag,N}
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

