# The guided preparation of a CAD target: an STL file of one solid that holds a hohlraum
# (its wall and caps with their entrance holes), a capsule and a foil, as a target
# fabrication export without bodies would give them.

@testset "guided CAD preparation (no GL)" begin
    dir = mktempdir()
    wall = merge_meshes(cylinder_mesh(2.5, -4.5, 4.5; nz = 12, nphi = 24, inward = true),
                        annulus_mesh(1.2, 2.5, 4.5; nr = 2, nphi = 24, up = false),
                        annulus_mesh(1.2, 2.5, -4.5; nr = 2, nphi = 24, up = true))
    capsule = sphere_mesh(0.9; ntheta = 10, nphi = 16)
    foil = translate(rectangle_mesh(3.0, 3.0; nu = 3, nv = 3), Briareus.Vec3(0.0, 0.0, -9.0))
    path = joinpath(dir, "assembly.stl")
    write_stl(path, merge_meshes(wall, capsule, foil))
    sess = DesignSession(starter = :none, resolution = :coarse)
    ts = sess.target
    @test add_stl_parts!(ts, path; h = 0.6)
    # step 1: one part of three connected pieces
    rows = cad_prep_rows(ts)
    p0 = rows[1].name
    @test length(rows) == 1 && rows[1].pieces == 3 &&
          rows[1].facets == sum(length(m.triangles) for m in (wall, capsule, foil))
    @test split_cad_part!(ts, rows[1].index)
    rows = cad_prep_rows(ts)
    @test [r.name for r in rows] == [Symbol(p0, :_, k) for k in 1:3] && all(r -> r.pieces == 1, rows)
    @test is_assembly(ts.parts[][end]) && ts.parts[][end].name === p0 &&
          all(p -> p.assembly === p0, ts.parts[][1:end-1])
    @test sum(r.facets for r in rows) == sum(length(m.triangles) for m in (wall, capsule, foil))
    @test !split_cad_part!(ts, rows[1].index) && occursin("one piece already", ts.error[])
    # step 2: the wall bounds the cavity at the center, the capsule is round and closed,
    # the foil is a sheet; which piece is which is told by its reach from the axis
    reach(n) = maximum(v -> hypot(v[1], v[2]), (v for c in ts.target[].components if c.part === n for v in c.mesh.vertices))
    names = Dict(r.name => (reach(r.name) > 2.3 ? :wall : reach(r.name) < 1 ? :capsule : :foil) for r in rows)
    s = suggest_cad_roles(sess)
    @test length(s.cavities) == 1 && s.cavities[1].walls == [k for (k, v) in names if v === :wall]
    @test all(r -> r.role === Dict(:wall => :wall, :capsule => :observed, :foil => :foil)[names[r.name]], s.rows)
    @test all(r -> r.current === :none, s.rows) && occursin("cavity 1 at (0, 0, 0) mm", s.text)
    # step 3: the roles are given, the capsule fitted with a sphere, and the wall shows
    # its cylinder in the preview
    @test isempty(apply_cad_roles!(sess, s.rows))
    @test all(r -> part_role(ts.parts[][part_index(ts, r.name)]) === r.role || names[r.name] === :capsule, s.rows)
    @test observed_name(ts) === :capsule && :foil in [names[n] for n in foil_names(ts)]
    pv = curvature_preview(sess)
    w = only(r for r in pv if get(names, r.part, :none) === :wall)
    # (the caps left flat: their share of the area of the wall)
    tri(x, y) = sqrt((x[2]y[3] - x[3]y[2])^2 + (x[3]y[1] - x[1]y[3])^2 + (x[1]y[2] - x[2]y[1])^2) / 2
    area(m) = sum(tri(m.vertices[b] - m.vertices[a], m.vertices[c] - m.vertices[a]) for (a, b, c) in m.triangles)
    caps = area(annulus_mesh(1.2, 2.5, 4.5; nr = 2, nphi = 24))
    @test w.cylinders == 1 && w.spheres == 0 && abs(w.flat - 2caps / area(wall)) < 0.005
    @test only(r for r in pv if r.part === :capsule).spheres == 1
end
