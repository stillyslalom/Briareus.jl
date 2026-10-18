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
