# Mesh repair for imported CAD: consistent orientation and coarsening into elements.

# a random subset of patches with reversed winding (as stitched CAD often has)
function scramble(m::Mesh, rng)
    owner = Briareus._triangle_patch(m)
    flip = rand(rng, Bool, npatches(m))
    return Mesh(m.vertices, [flip[owner[k]] ? (a, c, b) : (a, b, c)
                             for (k, (a, b, c)) in enumerate(m.triangles)], m.patches)
end
patch_normal(m, p) = normalize(Briareus._patch_area_vector(m, m.patches[p]))
patch_center(m, p) = Briareus._patch_centroid(m, m.patches[p])

@testset "orient_mesh" begin
    rng = Random.MersenneTwister(7)
    # a closed cube, every face scrambled: all faces out, positive volume
    cube = Mesh()
    for (n, u, v) in ((Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)), (Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 0, 0)),
                      (Vec3(0, 0, 1), Vec3(1, 0, 0), Vec3(0, 1, 0))), s in (1, -1)
        for i in 0:1, j in 0:1                     # four patches per face (shared edges)
            a, b = -1 + i, -1 + j
            add_polygon!(cube, [s * n + a * u + b * v, s * n + (a + 1) * u + b * v,
                                s * n + (a + 1) * u + (b + 1) * v, s * n + a * u + (b + 1) * v])
        end
    end
    for trial in 1:3
        o = orient_mesh(scramble(cube, rng))
        @test all(dot(patch_normal(o, p), patch_center(o, p)) > 0 for p in 1:npatches(o))
        @test o.patches == cube.patches && o.vertices == cube.vertices
    end
    # a thick hohlraum wall meshed as a closed solid (outer and inner cylinders joined by
    # end annuli, each a separate mesh with its own vertices): out of the material, i.e.
    # the inner wall faces the axis and the outer wall away from it
    R, r, L, nphi = 3.0, 2.875, 10.0, 24
    shell = merge_meshes(cylinder_mesh(R, -L / 2, L / 2; nz = 8, nphi, inward = false),
                         cylinder_mesh(r, -L / 2, L / 2; nz = 8, nphi, inward = true),
                         annulus_mesh(r, R, L / 2; nr = 1, nphi, up = true),
                         annulus_mesh(r, R, -L / 2; nr = 1, nphi, up = false))
    o = orient_mesh(scramble(shell, rng))
    radial(p) = (c = patch_center(o, p); Vec3(c[1], c[2], 0))
    inner = [p for p in 1:npatches(o) if abs(norm(radial(p)) - r * cos(π / nphi)) < 0.01]
    outer = [p for p in 1:npatches(o) if abs(norm(radial(p)) - R * cos(π / nphi)) < 0.01]
    @test length(inner) == length(outer) == 8nphi
    @test all(dot(patch_normal(o, p), radial(p)) < 0 for p in inner)
    @test all(dot(patch_normal(o, p), radial(p)) > 0 for p in outer)
    # the thin (open) wall faces the point it is given; without one, the larger part of
    # its area keeps its winding
    wall = cylinder_mesh(r, -L / 2, L / 2; nz = 8, nphi, inward = true)
    for toward in ((0, 0, 0), (0.0, 0.5, 2.0))
        o = orient_mesh(scramble(wall, rng); toward)
        @test all(dot(patch_normal(o, p), patch_center(o, p) .* Vec3(1, 1, 0)) < 0 for p in 1:npatches(o))
    end
    o = orient_mesh(scramble(wall, rng); toward = (0, 0, 10r))    # far along the axis: ambiguous
    signs = [sign(dot(patch_normal(o, p), patch_center(o, p) .* Vec3(1, 1, 0))) for p in 1:npatches(o)]
    @test all(==(signs[1]), signs)                                 # consistent, either way
    w2 = scramble(wall, rng)
    o = orient_mesh(w2)
    flipped = count(p -> dot(patch_normal(o, p), patch_normal(w2, p)) < 0, 1:npatches(o))
    @test flipped <= npatches(o) ÷ 2
    # a Target built from the repaired shell radiates into the cavity from the inner wall
    t = Target([Component(:shell, orient_mesh(scramble(shell, rng)), gold())])
    e = [i for i in 1:nelements(t) if abs(norm(t.elem_centroid[i] .* Vec3(1, 1, 0)) - r * cos(π / nphi)) < 0.01]
    @test all(dot(t.elem_normal[i], t.elem_centroid[i] .* Vec3(1, 1, 0)) < 0 for i in e)
    # a non-manifold edge (three sheets on one edge) warns but orients the rest
    fan = Mesh()
    for ϕ in (0.0, 2π / 3, 4π / 3)
        add_polygon!(fan, [Vec3(0, 0, 0), Vec3(0, 0, 1), Vec3(cos(ϕ), sin(ϕ), 1), Vec3(cos(ϕ), sin(ϕ), 0)])
    end
    @test_logs (:warn, r"non-manifold") orient_mesh(fan)
end

@testset "coarsen_mesh" begin
    # a finely divided plate: patches merged up to the size, triangles and area unchanged
    plate = rectangle_mesh(4.0, 2.0; nu = 40, nv = 20)
    c = coarsen_mesh(plate; size = 0.8)
    @test sort(c.triangles) == sort(plate.triangles) && c.vertices == plate.vertices
    @test 2 * 4 * 2 / 0.8^2 <= npatches(c) < npatches(plate) / 10
    span(m, p) = (vs = unique(v for k in m.patches[p] for v in m.triangles[k]);
                  maximum(norm(m.vertices[a] - m.vertices[b]) for a in vs, b in vs))
    @test all(span(c, p) <= 0.8 + 1e-12 for p in 1:npatches(c))
    t, tc = Target([Component(:p, plate, gold())]), Target([Component(:p, c, gold())])
    @test sum(tc.elem_area) ≈ sum(t.elem_area) rtol = 1e-12
    @test all(n -> n ≈ Vec3(0, 0, 1), tc.elem_normal)
    # a curved surface: merged elements stay within the angle of planar
    s = sphere_mesh(1.0; ntheta = 48, nphi = 96)
    cs = coarsen_mesh(s; size = 0.3, angle = 10.0)
    @test npatches(cs) < npatches(s) / 3
    for p in 1:npatches(cs)
        n = patch_normal(cs, p)
        @test all(dot(normalize(Briareus._patch_area_vector(cs, k:k)), n) >= cosd(10.0) - 1e-12
                  for k in cs.patches[p] if norm(Briareus._patch_area_vector(cs, k:k)) > 0)
    end
    # never across a fold: the cube's faces stay apart however large the size
    cube = Mesh()
    for (n, u, v) in ((Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)), (Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 0, 0)),
                      (Vec3(0, 0, 1), Vec3(1, 0, 0), Vec3(0, 1, 0))), sg in (1, -1)
        for i in 0:1, j in 0:1
            a, b = -1 + i, -1 + j
            add_polygon!(cube, [sg * n + a * u + b * v, sg * n + (a + 1) * u + b * v,
                                sg * n + (a + 1) * u + (b + 1) * v, sg * n + a * u + (b + 1) * v];
                         facing = sg * n)
        end
    end
    cc = coarsen_mesh(cube; size = 100.0)
    @test npatches(cc) == 6
    # nor across an inconsistent winding
    rng = Random.MersenneTwister(3)
    @test npatches(coarsen_mesh(scramble(plate, rng); size = 0.8)) > npatches(c)
    @test npatches(coarsen_mesh(orient_mesh(scramble(plate, rng)); size = 0.8)) == npatches(c)
    @test_throws ArgumentError coarsen_mesh(plate; size = 0)
end

# one patch per triangle with its own vertices, as an STL file has them
function facet_soup(tris)
    m = Mesh()
    foreach(t -> add_polygon!(m, collect(t)), tris)
    return m
end
mesh_area(m) = sum(norm(Briareus._patch_area_vector(m, k:k)) for k in eachindex(m.triangles))
# welded edges used by one triangle (count, total length) and by more than two
function open_edges(m)
    id = Briareus._welded_ids(m)
    pos = Dict(id[i] => m.vertices[i] for i in eachindex(id))
    uses = Dict{Tuple{Int,Int},Int}()
    for (a, b, c) in m.triangles, (u, v) in ((id[a], id[b]), (id[b], id[c]), (id[c], id[a]))
        uses[minmax(u, v)] = get(uses, minmax(u, v), 0) + 1
    end
    one = [e for (e, n) in uses if n == 1]
    return length(one), sum((norm(pos[u] - pos[v]) for (u, v) in one); init = 0.0), count(>(2), values(uses))
end
# an element's area, largest facet tilt from its mean normal (degrees; facets thinner
# than 1e-6 of their length have no orientation), largest vertex distance from its
# centroid, diameter, and whether the centroid projects along the normal onto it
function element_quality(m, p)
    r = m.patches[p]
    tv(k) = Briareus.triangle_vertices(m, k)
    A = sum(norm(Briareus._patch_area_vector(m, k:k)) for k in r)
    n = normalize(Briareus._patch_area_vector(m, r))
    c = Briareus._patch_centroid(m, r)
    tilt, on = 0.0, false
    for k in r
        a, b, d = tv(k)
        nk = cross(b - a, d - a)
        norm(nk) > 1e-6 * maximum(sum(abs2, e) for e in (b - a, d - b, a - d)) &&
            (tilt = max(tilt, acosd(clamp(dot(normalize(nk), n), -1, 1))))
        q = c + dot(a - c, nk) / dot(n, nk) * n
        on |= all(dot(cross(y - x, q - x), nk) >= -1e-9 * dot(nk, nk) for (x, y) in ((a, b), (b, d), (d, a)))
    end
    vs = unique(v for k in r for v in m.triangles[k])
    reach = maximum(norm(m.vertices[v] - c) for v in vs)
    diam = maximum(norm(m.vertices[u] - m.vertices[v]) for u in vs, v in vs)
    return (; A, tilt, reach, diam, on)
end
# edges of an element's outline, straight runs merged (as its view-factor receiver has them)
function outline_edges(m, p)
    T = m.triangles[m.patches[p]]
    return length(Briareus._merge_straight(m.vertices, Briareus._boundary_edges(T, eachindex(T))))
end

@testset "remesh_mesh" begin
    rng = Random.MersenneTwister(11)
    # a 4 × 2 mm plate tessellated as 40 strips of two slivers spanning its length
    # (l²/A = 320), each facet with its own vertices, windings scrambled
    strips = vcat([[(Vec3(0, y, 0), Vec3(4, y, 0), Vec3(4, y + 0.05, 0)),
                    (Vec3(0, y, 0), Vec3(4, y + 0.05, 0), Vec3(0, y + 0.05, 0))] for y in 0:0.05:1.95]...)
    plate = orient_mesh(scramble(facet_soup(strips), rng); toward = (2, 1, 1))
    @test npatches(coarsen_mesh(plate; size = 0.4)) == 80          # merging cannot shorten them
    size = 0.4
    r = remesh_mesh(plate; size)
    @test mesh_area(r) ≈ 8 rtol = 1e-12
    q = [element_quality(r, p) for p in 1:npatches(r)]
    @test sum(x -> x.A, q) ≈ 8 rtol = 1e-12
    # about one element per size² (a square element of side `size` has that area), each
    # within `size` of its centroid, compact, and with the centroid on it
    @test 0.9 < npatches(r) * size^2 / 8 < 1.15
    @test all(x -> x.reach <= size * (1 + 1e-12) && x.diam / sqrt(x.A) < 2.2 && x.on, q)
    @test all(x -> x.A > 0.6size^2, q)
    # cut across the strips on common planes, along them at strip edges: rectangles
    @test all(p -> outline_edges(r, p) == 4, 1:npatches(r))
    t = Target([Component(:p, r, gold())])
    @test all(n -> n ≈ Vec3(0, 0, 1), t.elem_normal)
    # the same surface: its only open edges are the plate's outline
    @test open_edges(r)[2] ≈ 12 && open_edges(r)[3] == 0

    # a closed can (radius 1, length 2) as one sheet: wall strips spanning its length,
    # caps as fans of slivers from the axis, every facet with its own vertices, windings
    # scrambled; every other wall strip is split at the midpoint of its bottom chord,
    # which the cap's facet does not have (a T-junction, as between CAD faces)
    nphi = 48
    facets = NTuple{3,Vec3}[]
    for j in 1:nphi
        a, b = 2π * (j - 1) / nphi, 2π * j / nphi
        p1, p2 = Vec3(cos(a), sin(a), -1), Vec3(cos(b), sin(b), -1)
        p3, p4 = Vec3(cos(b), sin(b), 1), Vec3(cos(a), sin(a), 1)
        if iseven(j)
            m = (p1 + p2) / 2
            push!(facets, (p1, p3, m), (m, p3, p2))
        else
            push!(facets, (p1, p3, p2))
        end
        push!(facets, (p1, p4, p3), (Vec3(0, 0, -1), p2, p1), (Vec3(0, 0, 1), p4, p3))
    end
    can = orient_mesh(scramble(facet_soup(facets), rng))
    @test open_edges(can)[1] == 3nphi ÷ 2                           # the T-junctions
    size = 0.25
    r = remesh_mesh(can; size)
    A = mesh_area(can)
    @test mesh_area(r) ≈ A rtol = 1e-12
    @test open_edges(r)[1] == 0 && open_edges(r)[3] == 0             # closed and conforming
    q = [element_quality(r, p) for p in 1:npatches(r)]
    @test 0.9 < npatches(r) * size^2 / A < 1.2
    @test all(x -> x.tilt <= 10.0 + 1e-9 && x.reach <= size * (1 + 1e-12) && x.on, q)
    @test all(x -> x.diam / sqrt(x.A) < 2.2, q)
    # radiation inside it (turned inward): every element sees the whole enclosure (rows
    # sum to 1 within 0.6%, mean 0.14%; the prism below: 0.33%, mean 0.16%), and the wall
    # sends the caps what the same prism meshed by `cylinder_mesh` sends
    inside = Target([Component(:can, Mesh(r.vertices, [(a, c, b) for (a, b, c) in r.triangles], r.patches),
                               gold())])
    F = view_factors(inside)
    @test maximum(abs, sum(F; dims = 2) .- 1) < 1e-2
    @test sum(abs, sum(F; dims = 2) .- 1) / nelements(inside) < 2e-3
    prism = Target([Component(:wall, cylinder_mesh(1.0, -1, 1; nz = 8, nphi), gold()),
                    Component(:caps, merge_meshes(annulus_mesh(0.0, 1.0, -1; nr = 4, nphi, up = true),
                                                  annulus_mesh(0.0, 1.0, 1; nr = 4, nphi, up = false)), gold())])
    Fp = view_factors(prism)
    wall(t) = [e for e in 1:nelements(t) if abs(t.elem_normal[e][3]) < 0.5]
    cap(t) = [e for e in 1:nelements(t) if abs(t.elem_normal[e][3]) > 0.5]
    exchange(t, F) = sum(t.elem_area[i] * F[i, j] for i in wall(t), j in cap(t)) / sum(t.elem_area[wall(t)])
    @test exchange(inside, F) ≈ exchange(prism, Fp) rtol = 2e-3

    # a tube of 360 full-length sliver strips (radius 2, length 3: 35 µm by 3 mm), its
    # length a multiple of size/3, so a strip's diagonal, a hair longer than its sides,
    # would take one point more under equal spacing and stagger every cut. The strips are
    # cut on common planes and elements fill bands between them: every cut across the
    # strips is one straight chord per strip, shared by the elements on both sides
    nphi = 360
    ring(z) = [Vec3(2cos(2π * j / nphi), 2sin(2π * j / nphi), z) for j in 0:nphi-1]
    tube = Mesh([ring(-1.5); ring(1.5)], NTuple{3,Int}[], UnitRange{Int}[])
    for j in 1:nphi
        k = mod1(j + 1, nphi)
        push!(tube.triangles, (j, nphi + j, k), (nphi + j, nphi + k, k))
        push!(tube.patches, length(tube.triangles)-1:length(tube.triangles))
    end
    size = 0.3
    r = remesh_mesh(tube; size)
    A = mesh_area(tube)
    @test mesh_area(r) ≈ A rtol = 1e-12
    @test open_edges(r)[2] ≈ open_edges(tube)[2] && open_edges(r)[3] == 0
    q = [element_quality(r, p) for p in 1:npatches(r)]
    @test 0.9 < npatches(r) * size^2 / A < 1.1
    @test all(x -> x.tilt <= 10.0 + 1e-9 && x.reach <= size * (1 + 1e-12) && x.on && x.diam / sqrt(x.A) < 1.6, q)
    zext = [extrema(r.vertices[v][3] for k in r.patches[p] for v in r.triangles[k]) for p in 1:npatches(r)]
    @test length(unique(round(z; digits = 6) for e in zext for z in e)) == round(Int, 3 / size) + 1
    # outlines: a chord across each strip at the two cuts, two runs along strips
    strip_of(c) = floor(Int, mod(atan(c[2], c[1]), 2π) * nphi / 2π)
    nstrips(p) = length(unique(strip_of(sum(r.vertices[v] for v in r.triangles[k])) for k in r.patches[p]))
    @test all(p -> outline_edges(r, p) <= 2nstrips(p) + 2, 1:npatches(r))

    @test_throws ArgumentError remesh_mesh(plate; size = 0)
    @test_throws ArgumentError remesh_mesh(plate; size = 0.4, angle = 0)
end

@testset "read_stl" begin
    # a unit cube's 12 facets, every other one wound backward (as exported CAD can be)
    corners = [Vec3(x, y, z) for x in 0:1, y in 0:1, z in 0:1]
    quads = (((1, 3, 4, 2)), ((5, 6, 8, 7)), ((1, 2, 6, 5)), ((3, 7, 8, 4)), ((1, 5, 7, 3)), ((2, 4, 8, 6)))
    facets = NTuple{3,Vec3}[]
    for (i, q) in enumerate(quads)
        a, b, c, d = corners[collect(q)]
        push!(facets, (a, b, c), isodd(i) ? (a, c, d) : (a, d, c))
    end
    dir = mktempdir()
    ascii = joinpath(dir, "cube_ascii.stl")
    open(ascii, "w") do io
        println(io, "solid cube")
        for (a, b, c) in facets
            println(io, "  facet normal 0 0 0\n    outer loop")
            for v in (a, b, c)
                println(io, "      vertex ", join(v, " "))
            end
            println(io, "    endloop\n  endfacet")
        end
        println(io, "endsolid cube")
    end
    binary = joinpath(dir, "cube_binary.stl")
    open(binary, "w") do io
        write(io, rpad("solid but binary", 80))           # binary files may start with "solid"
        write(io, UInt32(length(facets)))
        for (a, b, c) in facets
            write(io, zeros(Float32, 3), Float32.(a), Float32.(b), Float32.(c), UInt16(0))
        end
    end
    for path in (ascii, binary)
        m = read_stl(path; scale = 2.0)
        @test npatches(m) == 12 && length(m.vertices) == 36
        @test extrema(v[1] for v in m.vertices) == (0.0, 2.0)
        o = orient_mesh(m)
        vol(mm) = sum(dot(Briareus.triangle_vertices(mm, k)[1], cross(Briareus.triangle_vertices(mm, k)[2],
                                                                     Briareus.triangle_vertices(mm, k)[3]))
                      for k in eachindex(mm.triangles)) / 6
        @test vol(o) ≈ 8.0 rtol = 1e-12
        t = Target([Component(:cube, o, gold())])
        @test all(dot(t.elem_normal[e], t.elem_centroid[e] - Vec3(1, 1, 1)) > 0 for e in 1:nelements(t))
    end
    bad = joinpath(dir, "bad.stl")
    write(bad, "not a mesh")
    @test_throws ArgumentError read_stl(bad)
end

