# Mesh repair and coarsening for imported (CAD) meshes: consistent orientation of the
# patches, and merging of near-coplanar patches into elements of a chosen size.

# Vertex ids welded by position (to `tol` times the mesh's extent), so patches that
# share an edge in space share it in the mesh even when the file duplicates vertices.
function _welded_ids(m::Mesh; tol = 1e-9)
    isempty(m.vertices) && return Int[]
    lo = reduce((a, b) -> min.(a, b), m.vertices)
    hi = reduce((a, b) -> max.(a, b), m.vertices)
    h = tol * max(norm(hi - lo), eps())
    # cells of size h; a vertex joins a representative within h in its own or a
    # neighboring cell, so points straddling a cell boundary still weld
    cells = Dict{NTuple{3,Int},Vector{Int}}()
    rep = Vec3[]
    ids = Vector{Int}(undef, length(m.vertices))
    for (i, v) in enumerate(m.vertices)
        c = Tuple(floor.(Int, (v - lo) ./ h))
        found = 0
        for dx in -1:1, dy in -1:1, dz in -1:1
            for r in get(cells, (c[1] + dx, c[2] + dy, c[3] + dz), Int[])
                norm(rep[r] - v) <= h && (found = r; break)
            end
            found > 0 && break
        end
        if found == 0
            push!(rep, v)
            found = length(rep)
            push!(get!(cells, c, Int[]), found)
        end
        ids[i] = found
    end
    return ids
end

# Patch of each triangle.
function _triangle_patch(m::Mesh)
    owner = zeros(Int, length(m.triangles))
    for (p, r) in enumerate(m.patches), k in r
        owner[k] = p
    end
    return owner
end

# Edges between patches: undirected welded edge => [(patch, +1 if traversed from the
# lower to the higher vertex id, else -1)], over every triangle edge that is not interior
# to its patch.
function _patch_edges(m::Mesh, id)
    owner = _triangle_patch(m)
    edges = Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
    for (k, (a, b, c)) in enumerate(m.triangles)
        owner[k] == 0 && continue
        for (u, v) in ((id[a], id[b]), (id[b], id[c]), (id[c], id[a]))
            u == v && continue
            push!(get!(edges, minmax(u, v), Tuple{Int,Int}[]), (owner[k], u < v ? 1 : -1))
        end
    end
    # an edge traversed in both directions within one patch is interior to it
    for (e, list) in edges
        ps = first.(list)
        keep = [x for x in list if !(count(==(x[1]), ps) == 2 && any(y -> y[1] == x[1] && y[2] == -x[2], list))]
        length(keep) == length(list) || (edges[e] = keep)
    end
    filter!(kv -> !isempty(kv[2]), edges)
    return edges
end

_patch_area_vector(m::Mesh, r) =
    sum(_tri_cross(triangle_vertices(m, k)...) for k in r; init = zero(Vec3)) / 2
_patch_centroid(m::Mesh, r) =
    (A = sum(norm(_tri_cross(triangle_vertices(m, k)...)) for k in r);
     A > 0 ? sum(norm(_tri_cross(triangle_vertices(m, k)...)) * sum(triangle_vertices(m, k)) / 3
                 for k in r) / A : sum(sum(triangle_vertices(m, k)) / 3 for k in r) / length(r))

"""
    orient_mesh(mesh; toward = nothing, tol = 1e-9) -> Mesh

The mesh with its patches' windings made consistent: across every edge two patches
share (vertices welded by position to `tol` times the mesh's extent), the patches
traverse it in opposite directions, so their front sides lie on the same side of the
surface. Imported CAD is often stitched from faces of either orientation (a radiating
element's front is the side its winding defines).

Each connected surface then gets one orientation as a whole:
- a closed surface (every edge shared by exactly two patches) faces out of the volume
  it encloses (positive signed volume), so a thick hohlraum wall, meshed as a closed
  solid, faces into the cavity on its inner side;
- an open surface faces the point `toward` (e.g. the target center for a thin hohlraum
  wall) by area-weighted majority, or, without `toward`, keeps the orientation of the
  larger part of its area.

Edges shared by more than two patches (non-manifold) do not constrain the orientation,
and a surface that cannot be oriented consistently (a Möbius-like stitching) keeps the
first orientation found for each patch; both are reported in a warning.
"""
function orient_mesh(m::Mesh; toward = nothing, tol = 1e-9)
    np = npatches(m)
    np == 0 && return m
    id = _welded_ids(m; tol)
    edges = _patch_edges(m, id)
    # neighbors: (patch, same) with `same` true when the two must keep relative winding
    nbr = [Tuple{Int,Bool}[] for _ in 1:np]
    boundary = falses(np)
    nonmanifold = 0
    for list in values(edges)
        if length(list) == 1
            boundary[list[1][1]] = true
        elseif length(list) == 2
            (p, dp), (q, dq) = list
            p == q && continue
            push!(nbr[p], (q, dp != dq))       # opposite traversal: consistent as they are
            push!(nbr[q], (p, dp != dq))
        else
            nonmanifold += 1
            for (p, _) in list
                boundary[p] = true
            end
        end
    end
    flip = falses(np)
    shell = zeros(Int, np)
    conflicts = 0
    nshell = 0
    for s in 1:np
        shell[s] == 0 || continue
        nshell += 1
        shell[s] = nshell
        queue = [s]
        while !isempty(queue)
            p = popfirst!(queue)
            for (q, same) in nbr[p]
                want = same ? flip[p] : !flip[p]
                if shell[q] == 0
                    shell[q] = nshell
                    flip[q] = want
                    push!(queue, q)
                elseif flip[q] != want
                    conflicts += 1
                end
            end
        end
    end
    conflicts ÷= 2                                  # each seen from both sides
    # one orientation per connected surface
    areavec = [_patch_area_vector(m, r) for r in m.patches]
    sgn(p) = flip[p] ? -1 : 1
    for s in 1:nshell
        ps = findall(==(s), shell)
        closed = !any(boundary[p] for p in ps)
        score = if closed
            sum(sgn(p) * sum(dot(triangle_vertices(m, k)[1], cross(triangle_vertices(m, k)[2],
                                                                    triangle_vertices(m, k)[3]))
                             for k in m.patches[p]) for p in ps)       # 6 × signed volume
        elseif toward !== nothing
            c = Vec3(toward...)
            sum(norm(areavec[p]) * sign(sgn(p) * dot(areavec[p], c - _patch_centroid(m, m.patches[p])))
                for p in ps)
        else
            sum((flip[p] ? -1 : 1) * norm(areavec[p]) for p in ps)   # keep the larger part
        end
        if score < 0
            for p in ps
                flip[p] = !flip[p]
            end
        end
    end
    (nonmanifold > 0 || conflicts > 0) &&
        @warn "orient_mesh: $nonmanifold non-manifold edges, $conflicts inconsistent patch pairs"
    owner = _triangle_patch(m)
    tris = [owner[k] > 0 && flip[owner[k]] ? (a, c, b) : (a, b, c)
            for (k, (a, b, c)) in enumerate(m.triangles)]
    return Mesh(m.vertices, tris, m.patches)
end

"""
    coarsen_mesh(mesh; size, angle = 10.0, tol = 1e-9) -> Mesh

The mesh with neighboring patches merged into larger ones (the radiosity elements of a
target): each new patch grows from a seed patch over shared edges, taking patches whose
normals are within `angle` degrees of the seed's and whose vertices lie within
`size / 2` mm of the seed's centroid. Triangles are unchanged; only their grouping into
patches (and their order) is. Imported CAD often has hundreds of thousands of facets
(fine tessellation of curved faces, slivers along edges), far more than view factors
need: coarsening to a few hundred micrometres makes such parts usable for radiation.

Patches are merged only across edges where both are oriented consistently (see
`orient_mesh`), so the front side of a merged patch is well defined, and an element
stays within `angle` of planar: the flat-patch error of a merged element on a curved
surface is that of a facet of the same size. Seeds are taken by decreasing area, so
large flat faces are kept whole up to `size`.
"""
function coarsen_mesh(m::Mesh; size::Real, angle::Real = 10.0, tol = 1e-9)
    size > 0 || throw(ArgumentError("size must be positive"))
    np = npatches(m)
    np == 0 && return m
    id = _welded_ids(m; tol)
    nbr = [Int[] for _ in 1:np]
    for list in values(_patch_edges(m, id))
        length(list) == 2 || continue
        (p, dp), (q, dq) = list
        p != q && dp != dq || continue                 # consistently oriented neighbors
        push!(nbr[p], q)
        push!(nbr[q], p)
    end
    av = [_patch_area_vector(m, r) for r in m.patches]
    area = norm.(av)
    nrm = [a > 0 ? n / a : zero(Vec3) for (n, a) in zip(av, area)]
    cen = [_patch_centroid(m, r) for r in m.patches]
    pverts = [unique(v for k in r for v in m.triangles[k]) for r in m.patches]
    cosmax = cosd(angle)
    rmax = size / 2
    group = zeros(Int, np)
    order = sortperm(area; rev = true)
    ngroup = 0
    for s in order
        group[s] == 0 || continue
        ngroup += 1
        group[s] = ngroup
        queue = [s]
        while !isempty(queue)
            p = popfirst!(queue)
            for q in nbr[p]
                group[q] == 0 || continue
                dot(nrm[q], nrm[s]) >= cosmax || continue
                all(v -> norm(m.vertices[v] - cen[s]) <= rmax, pverts[q]) || continue
                group[q] = ngroup
                push!(queue, q)
            end
        end
    end
    tris = NTuple{3,Int}[]
    patches = UnitRange{Int}[]
    members = [Int[] for _ in 1:ngroup]
    for p in 1:np
        push!(members[group[p]], p)
    end
    for ps in members
        t0 = length(tris) + 1
        for p in ps, k in m.patches[p]
            push!(tris, m.triangles[k])
        end
        push!(patches, t0:length(tris))
    end
    return Mesh(m.vertices, tris, patches)
end

# ---------------------------------------------------------------- remeshing

"""
    remesh_mesh(mesh; size, angle = 10.0, tol = 1e-9) -> Mesh

The mesh's surface divided into compact, near-planar radiosity elements of about `size`
mm, however its facets are shaped: on a smooth surface the elements average `size`² in
area (a square element of side `size` has the same) and every vertex lies within `size`
of its element's centroid. Imported CAD tessellates curved walls as long slivers (a
hohlraum wall as axial strips spanning its length), which `coarsen_mesh` cannot
shorten: it only groups whole facets. Here:

1. vertices are welded by position (to `tol` times the mesh's extent), and one-sided
   edges with other such edges' vertices on them (T-junctions between CAD faces) are
   split there, so faces that meet along an edge are connected;
2. every facet is split, in its own plane, into triangles with edges of at most
   `size / 3`; an edge gets the same points from both facets sharing it, so the refined
   triangulation is conforming;
3. the refined triangles are clustered: seeds `0.85 size` apart (a Poisson-disk sample)
   claim triangles nearest first across shared, consistently wound edges (see
   `orient_mesh`), taking those whose normals are within `angle` degrees of the seed's;
   Lloyd iterations move each seed to its element's centroid and mean normal. An
   element with a facet more than `angle` from its mean normal, a vertex farther than
   `size` from its centroid, or a centroid that does not project (along the normal)
   onto it, is halved until it has none;
4. within each element, the pieces of facets in one plane are joined and
   re-triangulated from their outline (a point is dropped only where it lies on a
   straight outline in every element using it, so the mesh stays conforming).

The surface is unchanged: every output triangle lies in the plane of an input facet,
and the elements tile the input facets exactly (areas agree to round-off). Facets of
zero area (height below 1e-12 of their longest edge) and triangles outside every patch
are dropped; slivers thinner than 1e-6 of their length have no orientation to keep and
join any neighboring element. The input's patches do not constrain the elements; its
windings do (orient it first). Faces narrower than `size` (a wall's end, a fillet) give
elements as narrow, and elements on a faceted curved wall keep its faceting (a facet of
width w on radius R lies up to w²/8R off the true surface).
"""
function remesh_mesh(m::Mesh; size::Real, angle::Real = 10.0, tol = 1e-9)
    size > 0 || throw(ArgumentError("size must be positive"))
    0 < angle < 90 || throw(ArgumentError("angle must be between 0 and 90 degrees"))
    npatches(m) == 0 && return m
    P, tris = _welded_triangles(m, tol)
    isempty(tris) && return Mesh()
    lo = reduce((a, b) -> min.(a, b), P)
    hi = reduce((a, b) -> max.(a, b), P)
    chains = _tjunction_chains(P, tris, tol * max(norm(hi - lo), eps()))
    rtris, parent = _refine_facets!(P, tris, chains, size / 3)
    region = _cluster_elements(P, rtris, Float64(size), Float64(angle))
    return _assemble_elements(P, rtris, _coplanar_classes(P, tris, parent, region, size), region)
end

# Welded positions and the patch triangles in welded ids, without degenerate facets
# (repeated vertices, or height below 1e-12 of the longest edge).
function _welded_triangles(m::Mesh, tol)
    id = _welded_ids(m; tol)
    P = Vector{Vec3}(undef, isempty(id) ? 0 : maximum(id))
    for i in reverse(eachindex(id))
        P[id[i]] = m.vertices[i]                  # the first vertex of each weld
    end
    owner = _triangle_patch(m)
    tris = NTuple{3,Int}[]
    for (k, (a, b, c)) in enumerate(m.triangles)
        owner[k] == 0 && continue
        t = (id[a], id[b], id[c])
        (t[1] == t[2] || t[2] == t[3] || t[1] == t[3]) && continue
        pa, pb, pc = P[t[1]], P[t[2]], P[t[3]]
        l2 = max(sum(abs2, pb - pa), sum(abs2, pc - pb), sum(abs2, pa - pc))
        norm(_tri_cross(pa, pb, pc)) > 1e-12 * l2 || continue
        push!(tris, t)
    end
    return P, tris
end

# One-sided edges with other one-sided edges' vertices on them (within `δ`):
# edge (lo, hi) => the vertices along it from lo to hi.
function _tjunction_chains(P, tris, δ)
    cnt = Dict{Tuple{Int,Int},Int}()
    for (a, b, c) in tris, (u, v) in ((a, b), (b, c), (c, a))
        e = minmax(u, v)
        cnt[e] = get(cnt, e, 0) + 1
    end
    open = [e for (e, n) in cnt if n == 1]
    chains = Dict{Tuple{Int,Int},Vector{Int}}()
    isempty(open) && return chains
    vs = unique!(vcat(first.(open), last.(open)))
    lens = sort!([norm(P[v] - P[u]) for (u, v) in open])
    cell = max(lens[cld(length(lens), 2)], 4δ, eps())
    lo = reduce((a, b) -> min.(a, b), P[vs])
    key(p) = Tuple(floor.(Int, (p - lo) ./ cell))
    grid = Dict{NTuple{3,Int},Vector{Int}}()
    for v in vs
        push!(get!(grid, key(P[v]), Int[]), v)
    end
    cand = Int[]
    hits = Tuple{Float64,Int}[]
    for (u, v) in open
        a, b = P[u], P[v]
        d = b - a
        L = norm(d)
        L > 2δ || continue
        c0, c1 = key(min.(a, b) .- δ), key(max.(a, b) .+ δ)
        empty!(cand)
        if prod(c1 .- c0 .+ 1) > length(vs)
            append!(cand, vs)
        else
            for i in c0[1]:c1[1], j in c0[2]:c1[2], k in c0[3]:c1[3]
                append!(cand, get(grid, (i, j, k), Int[]))
            end
        end
        empty!(hits)
        for w in cand
            (w == u || w == v) && continue
            s = dot(P[w] - a, d) / L
            δ < s < L - δ || continue
            norm(a + (s / L) * d - P[w]) <= δ || continue
            push!(hits, (s, w))
        end
        isempty(hits) && continue
        sort!(hits)
        chains[(u, v)] = [u; last.(hits); v]
    end
    return chains
end

# Points of segment p → q (with p, without q), equally spaced at most `h` apart; each
# segment's points are made once, so both facets sharing it use the same ones.
function _segment_points!(loop, P, segs, p, q, h)
    k = minmax(p, q)
    ids = get(segs, k, nothing)
    if ids === nothing
        a, b = P[k[1]], P[k[2]]
        n = max(1, ceil(Int, norm(b - a) / h - 1e-9))
        ids = Vector{Int}(undef, n + 1)
        ids[1], ids[end] = k[1], k[2]
        for i in 1:n-1
            push!(P, a + (b - a) * (i / n))
            ids[i+1] = length(P)
        end
        segs[k] = ids
    end
    if p == k[1]
        append!(loop, @view ids[1:end-1])
    else
        append!(loop, @view ids[end:-1:2])
    end
    return loop
end

function _edge_points!(loop, P, segs, chains, a, b, h)
    k = minmax(a, b)
    ch = get(chains, k, nothing)
    ch === nothing && return _segment_points!(loop, P, segs, a, b, h)
    seq = a == k[1] ? ch : reverse(ch)
    for i in 1:length(seq)-1
        _segment_points!(loop, P, segs, seq[i], seq[i+1], h)
    end
    return loop
end

# Every facet split in its plane into triangles with edges of at most `h`, conforming
# across shared edges: the triangles, and the facet each came from.
function _refine_facets!(P, tris, chains, h)
    segs = Dict{Tuple{Int,Int},Vector{Int}}()
    rtris = NTuple{3,Int}[]
    parent = Int[]
    loop = Int[]
    for (k, (a, b, c)) in enumerate(tris)
        empty!(loop)
        _edge_points!(loop, P, segs, chains, a, b, h)
        _edge_points!(loop, P, segs, chains, b, c, h)
        _edge_points!(loop, P, segs, chains, c, a, h)
        if length(loop) == 3
            push!(rtris, (a, b, c))
            push!(parent, k)
        else
            _split_polygon!(rtris, parent, k, P, copy(loop), h)
        end
    end
    return rtris, parent
end

_cross2(u, v) = u[1] * v[2] - u[2] * v[1]

# A convex planar polygon (a facet with points on its edges, counterclockwise about the
# facet normal) split by chords between its points into triangles with edges of at
# most `h`. While a piece is wider than `h` it is cut by the shortest chord leaving at
# least a third of its area on each side (with points added along the chord every
# `h`), so slivers are cut across; smaller pieces by the chord of least length² / area.
function _split_polygon!(rtris, parent, k, P, ids, h)
    a = P[ids[1]]
    nv = zero(Vec3)
    for i in 2:length(ids)-1
        nv += _tri_cross(a, P[ids[i]], P[ids[i+1]])
    end
    n = normalize(nv)
    far = argmax(i -> sum(abs2, P[ids[i]] - a), eachindex(ids))
    e1 = normalize(P[ids[far]] - a)
    e2 = cross(n, e1)
    xy = [SVector(dot(P[i] - a, e1), dot(P[i] - a, e2)) for i in ids]
    D = sqrt(maximum(sum(abs2, q) for q in xy))
    h2 = h^2
    stack = [collect(eachindex(ids))]
    S = Float64[]
    while !isempty(stack)
        poly = pop!(stack)
        np = length(poly)
        if np == 3
            push!(rtris, (ids[poly[1]], ids[poly[2]], ids[poly[3]]))
            push!(parent, k)
            continue
        end
        resize!(S, np + 1)
        S[1] = 0.0
        for l in 1:np
            S[l+1] = S[l] + _cross2(xy[poly[l]], xy[poly[mod1(l + 1, np)]])
        end
        A = S[np+1] / 2
        bi = bj = qi = qj = 0
        bscore = qlen = Inf
        diam2 = 0.0
        for i in 1:np-1, j in i+1:np
            d2 = sum(abs2, xy[poly[j]] - xy[poly[i]])
            diam2 = max(diam2, d2)
            (j == i + 1 || (i == 1 && j == np)) && continue
            A1 = (S[j] - S[i] + _cross2(xy[poly[j]], xy[poly[i]])) / 2
            amin = min(A1, A - A1)
            amin > 1e-10 * sqrt(d2) * D || continue          # not along a straight side
            if amin >= A / 3 && d2 < qlen
                qi, qj, qlen = i, j, d2
            end
            if d2 / amin < bscore
                bi, bj, bscore = i, j, d2 / amin
            end
        end
        if diam2 > h2 * (1 + 1e-9) && qi > 0
            bi, bj = qi, qj
        end
        if bi == 0               # too thin for any chord to tell its sides apart: a fan
            for l in 2:np-1
                push!(rtris, (ids[poly[1]], ids[poly[l]], ids[poly[l+1]]))
                push!(parent, k)
            end
            continue
        end
        u, v = poly[bi], poly[bj]
        nc = max(1, ceil(Int, sqrt(sum(abs2, xy[v] - xy[u])) / h - 1e-9))
        chord = Int[]                                  # new local points from u to v
        for t in 1:nc-1
            push!(P, P[ids[u]] + (P[ids[v]] - P[ids[u]]) * (t / nc))
            push!(ids, length(P))
            push!(xy, xy[u] + (xy[v] - xy[u]) * (t / nc))
            push!(chord, length(ids))
        end
        push!(stack, [poly[bi:bj]; reverse(chord)])
        push!(stack, [poly[bj:end]; poly[1:bi]; chord])
    end
    return rtris
end

# Binary min-heap of (key, triangle, seed) for growing elements nearest first.
struct _GrowHeap
    key::Vector{Float64}
    tri::Vector{Int}
    seed::Vector{Int}
end
_GrowHeap() = _GrowHeap(Float64[], Int[], Int[])

function _heap_swap!(H::_GrowHeap, i, j)
    H.key[i], H.key[j] = H.key[j], H.key[i]
    H.tri[i], H.tri[j] = H.tri[j], H.tri[i]
    H.seed[i], H.seed[j] = H.seed[j], H.seed[i]
end

function _heap_push!(H::_GrowHeap, k, t, s)
    push!(H.key, k); push!(H.tri, t); push!(H.seed, s)
    i = length(H.key)
    while i > 1 && H.key[i>>1] > H.key[i]
        _heap_swap!(H, i, i >> 1)
        i >>= 1
    end
end

function _heap_pop!(H::_GrowHeap)
    k, t, s = H.key[1], H.tri[1], H.seed[1]
    n = length(H.key)
    _heap_swap!(H, 1, n)
    pop!(H.key); pop!(H.tri); pop!(H.seed)
    n -= 1
    i = 1
    while 2i <= n
        c = 2i + 1 <= n && H.key[2i+1] < H.key[2i] ? 2i + 1 : 2i
        H.key[i] <= H.key[c] && break
        _heap_swap!(H, i, c)
        i = c
    end
    return k, t, s
end

# Neighbors of each triangle across edges shared by exactly two triangles traversing
# them in opposite directions (0 where there is none).
function _triangle_neighbors(tris)
    N = length(tris)
    E = Vector{NTuple{4,Int}}(undef, 3N)              # (lo, hi, triangle, forward)
    for (t, (a, b, c)) in enumerate(tris), (j, (u, v)) in enumerate(((a, b), (b, c), (c, a)))
        E[3(t-1)+j] = (min(u, v), max(u, v), t, Int(u < v))
    end
    sort!(E)
    nbr = zeros(Int, 3, N)
    i = 1
    while i <= 3N
        j = i
        while j < 3N && E[j+1][1] == E[i][1] && E[j+1][2] == E[i][2]
            j += 1
        end
        if j == i + 1 && E[i][4] != E[j][4] && E[i][3] != E[j][3]
            p, q = E[i][3], E[j][3]
            nbr[findfirst(==(0), view(nbr, :, p)), p] = q
            nbr[findfirst(==(0), view(nbr, :, q)), q] = p
        end
        i = j + 1
    end
    return nbr
end

# Grow the given seeds over unowned triangles, nearest (to the seed point) first.
function _grow_seeds!(owner, best, H, seeds, P, rtris, cen, nrm, thin, nbr, spos, snrm, stri, cosmax, rmax)
    for s in seeds
        t = stri[s]
        owner[t] == 0 || continue
        k = norm(cen[t] - spos[s])
        k < best[t] || continue
        best[t] = k
        _heap_push!(H, k, t, s)
    end
    while !isempty(H.key)
        _, t, s = _heap_pop!(H)
        owner[t] == 0 || continue
        owner[t] = s
        for q in view(nbr, :, t)
            q == 0 && break
            owner[q] == 0 || continue
            thin[q] || dot(nrm[q], snrm[s]) >= cosmax || continue
            a, b, c = rtris[q]
            max(norm(P[a] - spos[s]), norm(P[b] - spos[s]), norm(P[c] - spos[s])) <= rmax || continue
            k = norm(cen[q] - spos[s])
            k < best[q] || continue
            best[q] = k
            _heap_push!(H, k, q, s)
        end
    end
    return owner
end

# Element of each refined triangle. Poisson-disk seeds `spacing` apart (a seed per face
# direction) grow nearest first across shared edges, within `angle` of the seed normal
# (triangles no seed reaches seed elements of their own), and Lloyd iterations move
# each seed to its element's centroid and mean normal; teeth of the outlines are then
# smoothed. Elements more than `angle` from planar, reaching beyond `rmax` of their
# centroid, or with a centroid off the element are halved until they are not. The
# spacing 0.85 `size` makes the mean element area ≈ `size`² on smooth surfaces.
function _cluster_elements(P, rtris, size::Float64, angle::Float64; iterations::Int = 5,
                           spacing::Float64 = 0.85size, rmax::Float64 = size,
                           smooth::Int = 3)
    N = length(rtris)
    cen = Vector{Vec3}(undef, N)
    nrm = Vector{Vec3}(undef, N)
    area = Vector{Float64}(undef, N)
    thin = falses(N)              # height below 1e-6 of the length: no orientation to keep
    for (i, (a, b, c)) in enumerate(rtris)
        cr = _tri_cross(P[a], P[b], P[c])
        A = norm(cr)
        area[i] = A / 2
        nrm[i] = A > 0 ? cr / A : zero(Vec3)
        cen[i] = (P[a] + P[b] + P[c]) / 3
        thin[i] = A <= 1e-6 * max(sum(abs2, P[b] - P[a]), sum(abs2, P[c] - P[b]), sum(abs2, P[a] - P[c]))
    end
    nbr = _triangle_neighbors(rtris)
    cosmax = cosd(angle)
    order = randperm(MersenneTwister(1), N)
    # Poisson-disk seeds: a triangle seeds unless a seed facing the same way is closer
    # than the spacing
    spos, snrm, stri = Vec3[], Vec3[], Int[]
    lo = reduce((a, b) -> min.(a, b), cen)
    key(p) = Tuple(floor.(Int, (p - lo) ./ spacing))
    grid = Dict{NTuple{3,Int},Vector{Int}}()
    for t in order
        thin[t] && continue
        c = key(cen[t])
        taken = false
        for dx in -1:1, dy in -1:1, dz in -1:1
            for s in get(grid, (c[1] + dx, c[2] + dy, c[3] + dz), Int[])
                if norm(spos[s] - cen[t]) < spacing && dot(snrm[s], nrm[t]) >= cosmax
                    taken = true
                    break
                end
            end
            taken && break
        end
        taken && continue
        push!(spos, cen[t]); push!(snrm, nrm[t]); push!(stri, t)
        push!(get!(grid, c, Int[]), length(spos))
    end
    owner = zeros(Int, N)
    best = fill(Inf, N)
    H = _GrowHeap()
    reach = 2size                 # bounds a seed's reach (halving enforces `rmax`)
    for it in 0:iterations
        fill!(owner, 0)
        fill!(best, Inf)
        _grow_seeds!(owner, best, H, eachindex(spos), P, rtris, cen, nrm, thin, nbr, spos, snrm, stri, cosmax, reach)
        for pass in 1:2, t in order                      # triangles no seed reached
            owner[t] == 0 && (pass == 2 || !thin[t]) || continue
            push!(spos, cen[t]); push!(snrm, nrm[t]); push!(stri, t)
            _grow_seeds!(owner, best, H, length(spos):length(spos), P, rtris, cen, nrm, thin, nbr,
                         spos, snrm, stri, cosmax, reach)
        end
        it == iterations && break
        # Lloyd: each seed to its element's centroid and mean normal
        ns = length(spos)
        wa, wc, wn = zeros(ns), fill(zero(Vec3), ns), fill(zero(Vec3), ns)
        for t in 1:N
            s = owner[t]
            wa[s] += area[t]
            wc[s] += area[t] * cen[t]
            wn[s] += area[t] * nrm[t]
        end
        keep = [s for s in 1:ns if wa[s] > 0 && norm(wn[s]) > 0]
        dist = fill(Inf, ns)
        newtri = zeros(Int, ns)
        for t in 1:N
            s = owner[t]
            d = norm(cen[t] - wc[s] / wa[s])
            d < dist[s] && ((dist[s], newtri[s]) = (d, t))
        end
        spos = [wc[s] / wa[s] for s in keep]
        snrm = [normalize(wn[s]) for s in keep]
        stri = [newtri[s] for s in keep]
    end
    # smooth the outlines: a triangle with two of its edges on one other element (a
    # tooth of a staircase outline) joins that element when its normal is within
    # `angle` of the element's seed normal
    for _ in 1:smooth
        moved = 0
        for t in 1:N
            q1, q2, q3 = nbr[1, t], nbr[2, t], nbr[3, t]
            q3 == 0 && continue
            r1, r2, r3 = owner[q1], owner[q2], owner[q3]
            o = owner[t]
            r = r1 == r2 || r1 == r3 ? r1 : r2 == r3 ? r2 : 0
            (r == 0 || r == o || stri[o] == t) && continue
            thin[t] || dot(nrm[t], snrm[r]) >= cosmax || continue
            owner[t] = r
            moved += 1
        end
        moved == 0 && break
    end
    members = [Int[] for _ in eachindex(spos)]
    for t in 1:N
        push!(members[owner[t]], t)
    end
    fill!(owner, 0)
    nel = 0
    stack = [ms for ms in members if !isempty(ms)]
    lab = zeros(Int, N)
    while !isempty(stack)
        ms = pop!(stack)
        A = sum(area[t] for t in ms)
        c = A > 0 ? sum(area[t] * cen[t] for t in ms) / A : sum(cen[ms]) / length(ms)
        nv = sum(area[t] * nrm[t] for t in ms)
        nbar = norm(nv) > 0 ? nv / norm(nv) : nv
        flat = all(t -> thin[t] || dot(nrm[t], nbar) >= cosmax - 1e-12, ms)
        r = maximum(t -> maximum(v -> norm(P[v] - c), rtris[t]), ms)
        if length(ms) == 1 || (flat && r <= rmax && _centroid_on(P, rtris, ms, c, nbar))
            nel += 1
            owner[ms] .= nel
            continue
        end
        # halve: from the two triangles farthest apart in normal (if not flat) or in
        # position, each half grown nearest first in that measure
        if flat
            a = argmax(t -> norm(cen[t] - c), ms)
            b = argmax(t -> norm(cen[t] - cen[a]), ms)
        else
            a = argmin(t -> thin[t] ? 1.0 : dot(nrm[t], nbar), ms)
            b = argmin(t -> thin[t] ? 1.0 : dot(nrm[t], nrm[a]), ms)
        end
        a == b && (b = ms[ms[1] == a ? 2 : 1])
        lab[ms] .= -1
        empty!(H.key); empty!(H.tri); empty!(H.seed)
        _heap_push!(H, 0.0, a, 1)
        _heap_push!(H, 0.0, b, 2)
        while !isempty(H.key)
            _, t, s = _heap_pop!(H)
            lab[t] == -1 || continue
            lab[t] = s
            o = s == 1 ? a : b
            for q in view(nbr, :, t)
                q == 0 && break
                lab[q] == -1 || continue
                _heap_push!(H, flat ? norm(cen[q] - cen[o]) : 1 - dot(nrm[q], nrm[o]), q, s)
            end
        end
        for s in (1, 2, -1)                  # -1: unreached (not edge-connected)
            part = [t for t in ms if lab[t] == s]
            isempty(part) || push!(stack, part)
        end
        lab[ms] .= 0
    end
    return owner
end

# Does the line through `c` along `n` meet one of the triangles `ts`? (An element's
# centroid, moved along its normal, should land on the element: radiation is exchanged
# from there.)
function _centroid_on(P, rtris, ts, c, n)
    for t in ts
        a, b, d = P[rtris[t][1]], P[rtris[t][2]], P[rtris[t][3]]
        nt = _tri_cross(a, b, d)
        den = dot(n, nt)
        den > 0 || continue
        p = c + (dot(a - c, nt) / den) * n
        tol = -1e-9 * dot(nt, nt)
        dot(_tri_cross(a, b, p), nt) >= tol && dot(_tri_cross(b, d, p), nt) >= tol &&
            dot(_tri_cross(d, a, p), nt) >= tol && return true
    end
    return false
end

# Is `v` on the straight segment from `p` to `q` (to round-off)?
function _on_segment(p, v, q)
    d = q - p
    return norm(cross(d, v - p)) <= 1e-10 * sum(abs2, d) && dot(v - p, q - v) > 0
end

# The outline of a set of triangles as one counterclockwise vertex loop, or `nothing`
# when it is not a single simple loop.
function _outline(rtris, ts)
    edges = Set{Tuple{Int,Int}}()
    for t in ts
        a, b, c = rtris[t]
        push!(edges, (a, b), (b, c), (c, a))
    end
    next = Dict{Int,Int}()
    for (u, v) in edges
        (v, u) in edges && continue
        haskey(next, u) && return nothing                # pinched
        next[u] = v
    end
    isempty(next) && return nothing
    start = first(keys(next))
    loop = [start]
    v = next[start]
    while v != start
        push!(loop, v)
        length(loop) > length(next) && return nothing
        v = get(next, v, 0)
        v == 0 && return nothing
    end
    return length(loop) == length(next) ? loop : nothing
end

# Ear-clipping triangulation of a simple planar polygon (counterclockwise about `n`),
# taking the best-shaped ear each time; `nothing` if it gets stuck.
function _earclip(P, poly, n)
    m = length(poly)
    m < 3 && return nothing
    e1 = normalize(cross(n, abs(n[1]) < 0.9 ? Vec3(1, 0, 0) : Vec3(0, 1, 0)))
    e2 = cross(n, e1)
    o = P[poly[1]]
    xy = [SVector(dot(P[v] - o, e1), dot(P[v] - o, e2)) for v in poly]
    scale = maximum(sum(abs2, q) for q in xy)
    tol = 1e-12 * scale
    idx = collect(1:m)
    out = NTuple{3,Int}[]
    while length(idx) > 3
        best, bq = 0, -Inf
        r = length(idx)
        for t in 1:r
            i0, i1, i2 = idx[mod1(t - 1, r)], idx[t], idx[mod1(t + 1, r)]
            a, b, c = xy[i0], xy[i1], xy[i2]
            cr = _cross2(b - a, c - a)
            cr > tol || continue
            inside = false
            for s in idx
                (s == i0 || s == i1 || s == i2) && continue
                p = xy[s]
                if _cross2(b - a, p - a) >= -tol && _cross2(c - b, p - b) >= -tol &&
                   _cross2(a - c, p - c) >= -tol
                    inside = true
                    break
                end
            end
            inside && continue
            q = cr / (sum(abs2, b - a) + sum(abs2, c - b) + sum(abs2, a - c))
            q > bq && ((best, bq) = (t, q))
        end
        best == 0 && return nothing
        r = length(idx)
        push!(out, (poly[idx[mod1(best - 1, r)]], poly[idx[best]], poly[idx[mod1(best + 1, r)]]))
        deleteat!(idx, best)
    end
    _cross2(xy[idx[2]] - xy[idx[1]], xy[idx[3]] - xy[idx[1]]) > tol || return nothing
    push!(out, (poly[idx[1]], poly[idx[2]], poly[idx[3]]))
    return out
end

# Plane of each refined triangle within its element: facets of one element that lie in
# one plane (to round-off) share a class, so their pieces are re-triangulated together.
# Returns the classes and their unit normals.
function _coplanar_classes(P, tris, parent, region, size)
    plane = zeros(Int, length(parent))
    normals = Vec3[]
    points = Vec3[]
    byregion = sortperm(region)
    cls = Dict{Int,Int}()                       # facet => class, within one element
    for (i, t) in enumerate(byregion)
        i > 1 && region[byregion[i-1]] != region[t] && empty!(cls)
        p = parent[t]
        c = get(cls, p, 0)
        if c == 0
            a, b, d = tris[p]
            n = normalize(_tri_cross(P[a], P[b], P[d]))
            for q in values(cls)
                if dot(normals[q], n) > 0 && all(v -> abs(dot(normals[q], P[v] - points[q])) <= 1e-12 * size, (a, b, d))
                    c = q
                    break
                end
            end
            if c == 0
                push!(normals, n); push!(points, P[a])
                c = length(normals)
            end
            cls[p] = c
        end
        plane[t] = c
    end
    return plane, normals
end

# The output mesh: one patch per element, each plane's pieces in it re-triangulated
# from their outline without the points that are straight in every outline using them.
function _assemble_elements(P, rtris, (plane, pnormal), region)
    ord = sortperm(collect(zip(region, plane)))
    groups = UnitRange{Int}[]
    i = 1
    while i <= length(ord)
        j = i
        while j < length(ord) && region[ord[j+1]] == region[ord[i]] && plane[ord[j+1]] == plane[ord[i]]
            j += 1
        end
        push!(groups, i:j)
        i = j + 1
    end
    kept = falses(length(P))
    loops = Vector{Union{Nothing,Vector{Int}}}(nothing, length(groups))
    keepall(g) = for k in groups[g], v in rtris[ord[k]]
        kept[v] = true
    end
    for (g, r) in enumerate(groups)
        lp = length(r) == 1 ? nothing : _outline(rtris, view(ord, r))
        if lp === nothing
            keepall(g)
            continue
        end
        L = length(lp)
        for k in 1:L
            _on_segment(P[lp[mod1(k - 1, L)]], P[lp[k]], P[lp[mod1(k + 1, L)]]) || (kept[lp[k]] = true)
        end
        loops[g] = lp
    end
    clipped = Vector{Union{Nothing,Vector{NTuple{3,Int}}}}(nothing, length(groups))
    changed = true
    while changed                       # a failed clip keeps its points; redo the rest
        changed = false
        for g in eachindex(groups)
            loops[g] === nothing && continue
            tri = _earclip(P, filter(v -> kept[v], loops[g]), pnormal[plane[ord[first(groups[g])]]])
            if tri === nothing
                loops[g] = nothing
                clipped[g] = nothing
                keepall(g)
                changed = true
            else
                clipped[g] = tri
            end
        end
    end
    tris = NTuple{3,Int}[]
    patches = UnitRange{Int}[]
    t0 = 1
    for (g, r) in enumerate(groups)
        if clipped[g] === nothing
            append!(tris, (rtris[ord[k]] for k in r))
        else
            append!(tris, clipped[g])
        end
        if g == length(groups) || region[ord[first(groups[g+1])]] != region[ord[first(r)]]
            push!(patches, t0:length(tris))
            t0 = length(tris) + 1
        end
    end
    used = zeros(Int, length(P))
    verts = Vec3[]
    for (a, b, c) in tris, v in (a, b, c)
        used[v] == 0 && (push!(verts, P[v]); used[v] = length(verts))
    end
    return Mesh(verts, [(used[a], used[b], used[c]) for (a, b, c) in tris], patches)
end
