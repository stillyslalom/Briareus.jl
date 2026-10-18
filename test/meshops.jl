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
    @test sum(tc.elem_area) ≈ sum(t.elem_area) rtol = 1e-14
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

