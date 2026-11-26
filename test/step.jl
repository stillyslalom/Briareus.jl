# STEP import through the Gmsh extension. The shapes are made with the OpenCASCADE
# kernel of Gmsh and written as STEP, so the exact surfaces are known.

function _step_tests()
    gmsh = Gmsh.gmsh
    dir = mktempdir()
    path = joinpath(dir, "shapes.step")
    # cylinder r 1.5 on z ∈ [0, 2]; cone r 1 at z = 0 to 0.4 at z = 3 about x = 5; sphere
    # r 1.2 at (0, 5, 0); torus R 2, r 0.5 at (5, 5, 1); box; half of a cylinder (one
    # face per half: a cylinder of angle π twice); a lofted free-form solid; a sheet
    begin
        gmsh.initialize(["-v", "0"])
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("shapes")
        occ = gmsh.model.occ
        occ.addCylinder(0, 0, 0, 0, 0, 2, 1.5)
        occ.addCone(5, 0, 0, 0, 0, 3, 1.0, 0.4)
        occ.addSphere(0, 5, 0, 1.2)
        occ.addTorus(5, 5, 1, 2.0, 0.5)
        occ.addBox(-5, -5, 0, 1, 2, 3)
        a = occ.addCylinder(-5, 5, 0, 0, 0, 1, 0.8, -1, π)
        b = occ.addCylinder(-5, 5, 0, 0, 0, 1, 0.8, -1, π)
        occ.rotate([(3, b)], -5, 5, 0, 0, 0, 1, π)
        gmsh.option.setNumber("Geometry.OCCUnionUnify", 0)     # keep the two faces of each half
        occ.fuse([(3, a)], [(3, b)])
        c1 = occ.addCircle(10, 0, 0, 1.0)
        e2 = occ.addEllipse(10, 0, 2, 1.0, 0.6)
        occ.addThruSections([occ.addCurveLoop([c1]), occ.addCurveLoop([e2])])
        occ.addRectangle(10, 5, 0, 2, 1)
        occ.synchronize()
        gmsh.write(path)
        gmsh.finalize()
    end
    ss = read_step(path; size = 0.4, angle = 10)
    @test length(ss) == 8
    names = [s.name for s in ss]
    @test allunique(names)
    kinds(s) = sort([f.kind for f in s.faces])
    cyl, cone, sph, tor, box, half, loft, sheet = ss[1], ss[2], ss[3], ss[4], ss[5], ss[6], ss[7], ss[end]
    @test kinds(cyl) == [:cylinder, :plane, :plane]
    @test kinds(box) == fill(:plane, 6)
    @test any(==(:other), kinds(loft))
    @test kinds(sheet) == [:plane]

    # the exact surfaces against the parameters of the shapes
    surf(s, kind) = only(f.surface for f in s.faces if f.kind === kind)
    sc = surf(cyl, :cylinder)
    @test sc isa CylinderSurface && isapprox(sc.radius, 1.5; rtol = 1e-12)
    @test abs(abs(sc.axis[3]) - 1) < 1e-12 && hypot(sc.origin[1], sc.origin[2]) < 1e-12
    sk = surf(cone, :cone)
    @test sk isa ConeSurface && isapprox(sk.half_angle, atand(0.2); rtol = 1e-12)
    # the file stores the half-angle to 11 digits, so the apex 5 mm away is good to 1e-10
    @test norm(sk.apex - Vec3(5, 0, 5)) < 1e-9 && sk.axis ≈ Vec3(0, 0, -1)
    ss_ = surf(sph, :sphere)
    @test ss_ isa SphereSurface && isapprox(ss_.radius, 1.2; rtol = 1e-12) && norm(ss_.center - Vec3(0, 5, 0)) < 1e-12
    st = surf(tor, :torus)
    @test st isa TorusSurface && isapprox(st.major_radius, 2.0; rtol = 1e-12) && isapprox(st.minor_radius, 0.5; rtol = 1e-12)
    @test norm(st.center - Vec3(5, 5, 1)) < 1e-12 && abs(abs(st.axis[3]) - 1) < 1e-12

    # the mesh vertices lie on the surfaces, and every solid faces out of its material
    @test all(r.maxdev < 1e-10 for s in ss for r in s.regions)
    vol(m) = sum(dot(Briareus.triangle_vertices(m, k)[1], cross(Briareus.triangle_vertices(m, k)[2],
                                                               Briareus.triangle_vertices(m, k)[3]))
                 for k in eachindex(m.triangles)) / 6
    # the facets cut inside the curved faces by up to 1.5% of the volume at this size
    @test vol(cyl.mesh) ≈ π * 1.5^2 * 2 rtol = 0.02
    @test vol(sph.mesh) ≈ 4π / 3 * 1.2^3 rtol = 0.02
    @test vol(tor.mesh) ≈ 2π^2 * 2 * 0.5^2 rtol = 0.02
    @test vol(box.mesh) ≈ 6 rtol = 1e-12
    @test all(r.facing == 1 for s in ss[1:end-1] for r in s.regions)
    # each face on a surface is in one region; the two halves of the split cylinder merge
    @test count(==(:cylinder), kinds(half)) == 2
    @test count(r -> r.kind === :cylinder, half.regions) == 1
    @test all(f.kind === :other || count(r -> issubset(f.patches, r.patches), s.regions) == 1
              for s in ss for f in s.faces)
    @test sum(npatches(s.mesh) for s in ss) == sum(length(f.patches) for s in ss for f in s.faces)

    # the same surfaces as fit_surfaces finds on the same triangles written as STL (a
    # whole torus is not compared: fit_surfaces finds fillet quarters, not whole tubes)
    for (s, kind, exact) in ((cyl, :cylinder, sc), (cone, :cone, sk), (sph, :sphere, ss_))
        stl = write_stl(joinpath(dir, "$kind.stl"), s.mesh; binary = false)
        fit = [r for r in fit_surfaces(read_stl(stl)) if r.kind === kind]
        @test length(fit) == 1
        R(x) = x isa ConeSurface ? (x.half_angle, x.apex...) : x isa SphereSurface ? (x.radius, x.center...) :
               (x.radius,)
        @test all(isapprox.(R(fit[1].surface), R(exact); atol = 1e-6))
    end

    # components and targets carry the exact surfaces
    cs = step_components(cyl, gold())
    @test [c.surface isa CylinderSurface for c in cs] == [true, false]
    @test all(c.part === Symbol(cyl.name) for c in cs)
    t = cad_target([ss[1:4]; Component(:extra, sphere_mesh(0.5; ntheta = 4, nphi = 6), gold())];
                   material = Dict(Symbol(s.name) => gold() for s in ss))
    @test count(c -> c.surface !== nothing, t.components) == 4
    @test Set(parts(t)) == Set([Symbol.(names[1:4]); :extra])
    @test_throws "no material" cad_target(ss[1:1]; material = Dict{Symbol,Briareus.Material}())
    @test_throws "minsize" read_step(path; size = 0.4, minsize = 1.0)
    @test occursin("3 faces: 2 plane, 1 cylinder", sprint(show, cyl))
end

@testset "STEP import" begin
    @test_throws "needs the package Gmsh" read_step(joinpath(@__DIR__, "step.jl"))
    if Base.find_package("Gmsh") === nothing
        @info "Gmsh is not installed: the STEP import tests are skipped"
        @test_skip false
    else
        @eval using Gmsh
        Base.invokelatest(_step_tests)
    end
end
