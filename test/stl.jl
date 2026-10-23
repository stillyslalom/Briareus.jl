using Printf: @sprintf

@testset "STL interchange" begin
    dir = mktempdir()
    ntri(m) = length(m.triangles)
    tri_area(m, k) = norm(cross(Briareus.triangle_vertices(m, k)[2] - Briareus.triangle_vertices(m, k)[1],
                                Briareus.triangle_vertices(m, k)[3] - Briareus.triangle_vertices(m, k)[1])) / 2
    tri_normal(m, k) = normalize(cross(Briareus.triangle_vertices(m, k)[2] - Briareus.triangle_vertices(m, k)[1],
                                       Briareus.triangle_vertices(m, k)[3] - Briareus.triangle_vertices(m, k)[1]))
    total_area(m) = sum(k -> tri_area(m, k), eachindex(m.triangles); init = 0.0)

    @testset "round trip" begin
        # a small hohlraum: wall, end caps, two LEH apertures, a capsule
        t = hohlraum(nz = 6, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8)
        names = [String(c.name) for c in t.components]
        @test any(c -> c.kind === :aperture, t.components)
        for binary in (true, false), scale in (1.0, 10.0)
            path = write_stl(joinpath(dir, "h_$(binary)_$scale.stl"), t; binary, scale)
            raw = read_stl(path)                        # in file units
            m = read_stl(path; scale)
            all_tris = sum(c -> ntri(c.mesh), t.components)
            @test ntri(m) == npatches(m) == all_tris
            @test maximum(norm, raw.vertices) ≈ maximum(norm, m.vertices) / scale rtol = 1e-6
            # areas and windings of every triangle, in component order; Float32 coordinates
            # in binary move a vertex by ≤ 6e-8 of its size (here ≲ 5 mm): 1e-5 relative
            # on areas and 1e-5 on unit normals of the smallest facets is conservative
            tol = binary ? 1e-5 : 1e-12
            src = [(c.mesh, j) for c in t.components for j in eachindex(c.mesh.triangles)]
            @test all(isapprox(tri_area(m, k), tri_area(s...); rtol = tol) for (k, s) in enumerate(src))
            @test maximum(norm(tri_normal(m, k) - tri_normal(s...)) for (k, s) in enumerate(src)) < tol
            if !binary
                solids = read_stl_solids(path; scale)
                @test first.(solids) == names                  # names and order
                for ((n, sm), c) in zip(solids, t.components)
                    @test total_area(sm) ≈ total_area(c.mesh) rtol = 1e-12
                end
            else
                solids = read_stl_solids(path; scale)
                @test first.(solids) == ["h_true_$scale"]      # binary: one solid, the file name
                @test startswith(String(read(path)[1:80]), "Briareus.jl")   # not "solid"
            end
        end
        # facet normals in the file come from the winding (front side)
        c = Component(:plate, rectangle_mesh(2.0, 2.0; nu = 1, nv = 1), gold())
        path = write_stl(joinpath(dir, "plate.stl"), c; binary = false)
        @test occursin("facet normal 0.0 0.0 1.0", read(path, String))
        @test first.(read_stl_solids(path)) == ["plate"]
        bpath = write_stl(joinpath(dir, "plate.stl"), c)
        @test reinterpret(Float32, read(bpath)[85:96]) == Float32[0, 0, 1]
        # apertures: written as their own solids by default, left out on request
        noleh = read_stl_solids(write_stl(joinpath(dir, "noleh.stl"), t; binary = false, apertures = false))
        @test first.(noleh) == [String(c.name) for c in t.components if c.kind !== :aperture]
        # closing lips are written with their entrance holes: each hole whole
        tl = hohlraum(nz = 6, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8,
                      leh_closure = GapClosure(closure = 5.0))
        sl = read_stl_solids(write_stl(joinpath(dir, "lips.stl"), tl; binary = false))
        @test first.(sl) == names
        leh(ss) = only(m for (n, m) in ss if n == "leh_top")
        @test total_area(leh(sl)) ≈ total_area(leh(read_stl_solids(joinpath(dir, "h_false_1.0.stl")))) rtol = 1e-9
        @test first.(read_stl_solids(write_stl(joinpath(dir, "lips.stl"), tl; binary = false, apertures = false))) ==
              first.(noleh)
        lone = read_stl_solids(write_stl(joinpath(dir, "lone.stl"), [tl.components[component_index(tl, :leh_top_lip)]];
                                         binary = false))
        @test first.(lone) == ["leh_top_lip"]
        # a two-sided component is written once; a mesh alone, pairs, a component vector
        f = foil(:f, (1.0, 0.0, 0.0), (1, 0, 0); width = 1.0, n = 2)
        @test f.two_sided
        @test ntri(read_stl(write_stl(joinpath(dir, "f.stl"), f))) == ntri(f.mesh)
        mesh = t.components[1].mesh
        @test first.(read_stl_solids(write_stl(joinpath(dir, "m.stl"), mesh; binary = false))) == ["briareus"]
        @test first.(read_stl_solids(write_stl(joinpath(dir, "m.stl"), mesh; binary = false, name = "wall"))) == ["wall"]
        p = write_stl(joinpath(dir, "pairs.stl"), ["a" => mesh, :b => f.mesh]; binary = false)
        @test first.(read_stl_solids(p)) == ["a", "b"]
        @test ntri(read_stl(p)) == ntri(mesh) + ntri(f.mesh)
        p = write_stl(joinpath(dir, "comps.stl"), [f, t.components[1]]; binary = false)
        @test first.(read_stl_solids(p)) == ["f", String(t.components[1].name)]
        @test_throws ArgumentError write_stl(joinpath(dir, "x.stl"), mesh; scale = 0)
    end

    @testset "foreign files" begin
        facet(io, a, b, c; kw = ("facet normal", "outer loop", "vertex", "endloop", "endfacet"), nl = "\n") =
            print(io, "  ", kw[1], " 1.0E+00 -0.0e0 0.5e-3", nl, "    ", kw[2], nl,
                  ("      " * kw[3] * " " * join(v, " ") * nl for v in (a, b, c))..., "    ", kw[4], nl,
                  "  ", kw[5], nl)
        a, b, c, d = Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(1, 1, 0)
        # several solids from one CAD exporter: CRLF, tabs, upper case, exponents,
        # a solid name with spaces, an empty solid, an unnamed solid
        path = joinpath(dir, "assembly.stl")
        open(path, "w") do io
            print(io, "solid outer shell\r\n")
            facet(io, a, b, c; nl = "\r\n")
            facet(io, b, d, c; nl = "\r\n")
            print(io, "endsolid outer shell\r\n")
            print(io, "SOLID\tPART_2\n")
            facet(io, 2a, 2c, 2b; kw = ("FACET\tNORMAL", "OUTER LOOP", "VERTEX\t", "ENDLOOP", "ENDFACET"))
            print(io, "ENDSOLID PART_2\n")
            print(io, "solid empty\nendsolid empty\n")
            print(io, "solid\n")
            facet(io, Vec3(1e-3, 0, 0), Vec3(0, 2.5e2, 0), Vec3(0, 0, -3.0e-1))
            print(io, "endsolid")                       # no final newline
        end
        solids = read_stl_solids(path)
        @test first.(solids) == ["outer shell", "PART_2", "empty", "solid4"]
        @test ntri.(last.(solids)) == [2, 1, 0, 1]
        @test tri_normal(solids[1][2], 2) ≈ Vec3(0, 0, 1)
        @test tri_normal(solids[2][2], 1) ≈ Vec3(0, 0, -1)
        @test solids[4][2].vertices == [Vec3(1e-3, 0, 0), Vec3(0, 250, 0), Vec3(0, 0, -0.3)]
        m = read_stl(path; scale = 25.4)
        @test ntri(m) == npatches(m) == 4
        @test m.vertices[7:9] == 25.4 .* [2a, 2c, 2b]
        # an all-empty file is an empty mesh
        e = joinpath(dir, "empty.stl")
        write(e, "solid nothing\nendsolid nothing\n")
        @test ntri(read_stl(e)) == 0 && first.(read_stl_solids(e)) == ["nothing"]

        # a binary file whose header begins with "solid" (as SolidWorks writes them)
        bin = joinpath(dir, "solidworks.stl")
        tris = [(a, b, c), (b, d, c)]
        open(bin, "w") do io
            write(io, rpad("solid PART exported by a CAD tool", 80))
            write(io, UInt32(length(tris)))
            for (p, q, r) in tris
                write(io, Float32[0, 0, 1], Float32.(p), Float32.(q), Float32.(r), UInt16(0))
            end
        end
        @test read_stl(bin).vertices == [a, b, c, b, d, c]
        @test first.(read_stl_solids(bin)) == ["solidworks"]
        # ... truncated by a byte: neither binary (size) nor ASCII (NUL bytes): clear error
        trunc = joinpath(dir, "truncated.stl")
        write(trunc, read(bin)[1:end-1])
        err = try read_stl(trunc); nothing catch e; e end
        @test err isa ArgumentError && occursin("truncated", err.msg)
        nohdr = joinpath(dir, "odd.stl")
        write(nohdr, vcat(fill(0x01, 80), reinterpret(UInt8, [UInt32(3)]), zeros(UInt8, 7)))
        @test_throws ArgumentError read_stl(nohdr)

        # the fast number parser agrees with Base's parse bit for bit, in every format
        rng = Random.MersenneTwister(7)
        xs = [randn(rng) * 10.0^rand(rng, -30:30) for _ in 1:2000]
        strs = vcat(string.(xs), [@sprintf("%e", x) for x in xs], [@sprintf("%.9E", x) for x in xs],
                    [@sprintf("%.3f", x) for x in xs[1:200]], [@sprintf("%+.6g", x) for x in xs],
                    ["0", "-0", "+1", "1.", ".5", "-.5e+3", "1e22", "1e23", "9007199254740993",
                     "123456789012345678901", "1e-400", "4.9e-324", "inf", "-Inf", "NaN"])
        @test all(s -> isequal(Briareus._stl_float(s, 1, ncodeunits(s)), tryparse(Float64, s)), strs)
        @test all(s -> Briareus._stl_float(s, 1, ncodeunits(s)) === nothing, ["e5", "1e", "1e+", "1.2.3", "-", "1x"])

        # malformed ASCII
        bad = joinpath(dir, "bad.stl")
        write(bad, "solid x\n facet normal 0 0 1\n outer loop\n vertex 0 0 0\n vertex 1 0 0\n endloop\n endfacet\nendsolid x\n")
        err = try read_stl(bad); nothing catch e; e end
        @test err isa ArgumentError && occursin("3 vertices, not 2", err.msg) && occursin("line 7", err.msg)
        write(bad, "solid x\n facet normal 0 0 1\n outer loop\n vertex 0 zero 0\n")
        err = try read_stl(bad); nothing catch e; e end
        @test err isa ArgumentError && occursin("line 4", err.msg)
        write(bad, "solid x\n vertex 0 0")
        @test_throws ArgumentError read_stl(bad)
        write(bad, "")
        @test_throws ArgumentError read_stl(bad)
        write(bad, "not a mesh")
        @test_throws ArgumentError read_stl(bad)
    end
end
