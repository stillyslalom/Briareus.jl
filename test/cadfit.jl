# Analytic surfaces in CAD-like meshes (fit_surfaces, curved_components). Expects
# `facet_soup` from meshops.jl.

# CAD-like facets: curved faces as full-length strips, flat annuli as radial slivers,
# every facet its own vertices (as STL has them), wound to face `facing`.
function cad_quad!(F, a, b, c, d, facing)
    dot(cross(b - a, c - a), facing) < 0 && ((a, b, c, d) = (d, c, b, a))
    push!(F, (a, b, c), (a, c, d))
end
cad_pt(r, φ, z) = Vec3(r * cos(φ), r * sin(φ), z)

# upper half of a hohlraum with a thick wall, as a closed solid: inner wall r = R
# (0 ≤ z ≤ h), outer wall Rout (0 ≤ z ≤ h + t), end wall h ≤ z ≤ h + t pierced by an
# entrance hole of radius Rleh (its bore a short cylinder), ring at the waist z = 0
function cad_half(; R, Rout, Rleh, h, t, nphi)
    F = NTuple{3,Vec3}[]
    φs = range(0, 2π; length = nphi + 1)
    for j in 1:nphi
        a, b = φs[j], φs[j+1]
        out = cad_pt(1, (a + b) / 2, 0)
        cad_quad!(F, cad_pt(R, a, 0), cad_pt(R, b, 0), cad_pt(R, b, h), cad_pt(R, a, h), -out)
        cad_quad!(F, cad_pt(Rout, a, 0), cad_pt(Rout, b, 0), cad_pt(Rout, b, h + t), cad_pt(Rout, a, h + t), out)
        cad_quad!(F, cad_pt(Rleh, a, h), cad_pt(Rleh, b, h), cad_pt(Rleh, b, h + t), cad_pt(Rleh, a, h + t), -out)
        cad_quad!(F, cad_pt(Rleh, a, h), cad_pt(R, a, h), cad_pt(R, b, h), cad_pt(Rleh, b, h), Vec3(0, 0, -1))
        cad_quad!(F, cad_pt(Rleh, a, h + t), cad_pt(Rout, a, h + t), cad_pt(Rout, b, h + t), cad_pt(Rleh, b, h + t),
                  Vec3(0, 0, 1))
        cad_quad!(F, cad_pt(R, a, 0), cad_pt(Rout, a, 0), cad_pt(Rout, b, 0), cad_pt(R, b, 0), Vec3(0, 0, -1))
    end
    return F
end
cad_mirror(F) = [(a .* Vec3(1, 1, -1), c .* Vec3(1, 1, -1), b .* Vec3(1, 1, -1)) for (a, b, c) in F]

# latitude-longitude sphere: quads of two facets, fans of slivers at the poles
function cad_sphere(r; ntheta, nphi, outward = true)
    F = NTuple{3,Vec3}[]
    θs = range(0, π; length = ntheta + 1)
    φs = range(0, 2π; length = nphi + 1)
    p(θ, φ) = r * Vec3(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
    for i in 1:ntheta, j in 1:nphi
        a, b, c, d = p(θs[i], φs[j]), p(θs[i+1], φs[j]), p(θs[i+1], φs[j+1]), p(θs[i], φs[j+1])
        f = (outward ? 1 : -1) * (a + b + c + d)
        tri(x, y, z) = dot(cross(y - x, z - x), f) >= 0 ? (x, y, z) : (x, z, y)
        i == 1 ? push!(F, tri(a, b, c)) : i == ntheta ? push!(F, tri(a, b, d)) : push!(F, tri(a, b, c), tri(a, c, d))
    end
    return F
end

# through binary STL and back: Float32 coordinates, one patch per facet
function stl_roundtrip(F)
    path = joinpath(mktempdir(), "part.stl")
    write_stl(path, facet_soup(F))
    return read_stl(path)
end

@testset "fit_surfaces" begin
    R, Rout, Rleh, h, t = 2.875, 2.975, 1.55, 4.7, 0.05
    raw = orient_mesh(stl_roundtrip(cad_half(; R, Rout, Rleh, h, t, nphi = 360)))
    A = sum(norm(Briareus._patch_area_vector(raw, r)) for r in raw.patches)
    for (mesh, rtol) in ((raw, 1e-6), (remesh_mesh(raw; size = 0.3), 3e-5))
        regs = fit_surfaces(mesh)
        @test sum(length(r.patches) for r in regs) == npatches(mesh)      # every patch on a surface
        @test sum(r.area for r in regs) ≈ A rtol = 1e-9
        cyl = sort([r for r in regs if r.kind === :cylinder]; by = r -> r.surface.radius)
        @test length(cyl) == 3 && count(r -> r.kind === :plane, regs) == 3 && length(regs) == 6
        # the bore, the inner wall (both facing the axis) and the outer wall (facing out)
        @test [r.surface.radius for r in cyl] ≈ [Rleh, R, Rout] rtol = rtol
        @test [r.facing for r in cyl] == [-1, -1, 1]
        @test all(r -> abs(r.surface.axis[3]) > 1 - 1e-10 && norm(r.surface.origin[1:2]) < 1e-5, cyl)
        @test all(r -> r.maxdev < 1.2 * R^2 * (2π / 360)^2 / 8, cyl)        # within the facets' sagitta
        planes = sort([r for r in regs if r.kind === :plane]; by = r -> r.centroid[3])
        @test [r.normal[3] for r in planes] ≈ [-1, -1, 1]
        @test [r.centroid[3] for r in planes] ≈ [0, h, h + t] atol = 1e-6
        @test all(r -> r.maxdev < 1e-6, planes)
    end
    # coarse strips (0.25 mm on the 2.875 mm wall, 2.7 µm sagitta): remeshed points on
    # the facets lie too far inside the wall to fit it, so the regions are fitted on the
    # CAD facets and carried over to the elements
    coarse = orient_mesh(stl_roundtrip(cad_half(; R, Rout, Rleh, h, t, nphi = 72)))
    rc = remesh_mesh(coarse; size = 0.3)
    regs = fit_surfaces(rc; reference = coarse)
    @test sum(length(r.patches) for r in regs) == npatches(rc)
    cyl = sort([r for r in regs if r.kind === :cylinder]; by = r -> r.surface.radius)
    @test [r.surface.radius for r in cyl] ≈ [Rleh, R, Rout] rtol = 1e-6
    @test all(r -> r.maxdev <= R^2 * (2π / 72)^2 / 8 * 1.01, cyl)       # the sagitta
    @test count(r -> r.kind === :cylinder && r.surface.radius ≈ R, fit_surfaces(rc)) == 0
    # a sphere, raw and remeshed
    sph = orient_mesh(stl_roundtrip(cad_sphere(0.9; ntheta = 56, nphi = 112)))
    for mesh in (sph, remesh_mesh(sph; size = 0.3))
        regs = fit_surfaces(mesh)
        @test length(regs) == 1 && regs[1].kind === :sphere && length(regs[1].patches) == npatches(mesh)
        @test regs[1].surface.radius ≈ 0.9 rtol = 1e-5
        @test norm(regs[1].surface.center) < 1e-6
        @test regs[1].facing == 1
    end
    # a free-form surface (a paraboloid) holds no plane, cylinder or sphere
    para = Mesh()
    for i in 0:19, j in 0:19
        z(x, y) = 0.3 * x^2 + 0.1 * y^2 + 0.05 * x * y
        p(x, y) = Vec3(x, y, z(x, y))
        x0, y0, d = -1 + 0.1i, -1 + 0.1j, 0.1
        add_polygon!(para, [p(x0, y0), p(x0 + d, y0), p(x0 + d, y0 + d)])
        add_polygon!(para, [p(x0, y0), p(x0 + d, y0 + d), p(x0, y0 + d)])
    end
    @test isempty(fit_surfaces(para))
    @test_throws ArgumentError fit_surfaces(para; kinds = (:cone,))
end

@testset "curved_components" begin
    R, Rout, Rleh, h, t = 2.875, 2.975, 1.55, 4.7, 0.05
    raw = orient_mesh(stl_roundtrip(cad_half(; R, Rout, Rleh, h, t, nphi = 360)))
    el = remesh_mesh(raw; size = 0.5)
    c = Component(:half, el, gold())
    cs = curved_components(c; reference = raw)
    @test [x.name for x in cs] == [Symbol("half/cyl1"), Symbol("half/cyl2"), Symbol("half/cyl3"), Symbol("half/rest")]
    @test [x.surface.radius for x in cs[1:3]] ≈ [Rout, R, Rleh] rtol = 1e-6
    @test cs[4].surface === nothing && sum(npatches(x.mesh) for x in cs) == npatches(el)
    @test all(x -> x.material === c.material && x.kind === :surface && !x.two_sided, cs)
    flat = Target([c])
    fitted = Target(cs)
    # the same facets; the curved elements get the exact areas of the wall
    @test sum(fitted.tri_area) ≈ sum(flat.tri_area) rtol = 1e-12
    @test sum(fitted.elem_area[elements(fitted, Symbol("half/cyl2"))]) ≈ 2π * R * h rtol = 1e-7
    # components with a surface, or without curved faces, are kept
    tp = hohlraum(; h = 1.0, capsule_h = 0.5)
    @test curved_components(tp).components == tp.components
    plate = Component(:plate, rectangle_mesh(2.0, 1.0; nu = 4, nv = 2), gold())
    @test curved_components(plate) == [plate]
    @test_throws ArgumentError curved_components(c; kinds = (:plane,))
end

@testset "radiation on fitted CAD" begin
    # interior of a CAD sphere (latitude-longitude strips, remeshed): on the fitted
    # sphere F_ij = A_j / 4πR² for every pair, self included, and an integrating sphere
    # has the same incident flux everywhere wherever the source is; flat elements miss
    # it by percents next to the source
    Rs = 2.0
    sraw = stl_roundtrip(cad_sphere(Rs; ntheta = 30, nphi = 60, outward = false))
    sel = remesh_mesh(sraw; size = 0.5)
    w = Material(:w; albedo = 0.5, conversion = 0.7)
    flat = Target([Component(:s, sel, w)])
    fitted = Target(curved_components(Component(:s, sel, w); reference = sraw))
    @test length(fitted.components) == 1 && fitted.components[1].surface isa SphereSurface
    @test sum(fitted.elem_area) ≈ 4π * Rs^2 rtol = 1e-7          # the radius fitted to Float32 vertices
    P = zeros(nelements(fitted)); P[1:3] .= 1.0
    spread = Float64[]
    for t in (flat, fitted)
        F = view_factors(t)
        A = t.elem_area
        sol = solve(RadiationModel(t; F), P)
        exact = 0.7 * sum(P) / (0.5 * sum(A))
        push!(spread, maximum(abs, sol.H ./ exact .- 1))
        t === fitted || continue
        @test maximum(abs.(F .- (A ./ sum(A))')) < 1e-12
        @test maximum(abs, closure_error(F)) < 1e-12
    end
    @test spread[2] < 1e-12 && spread[1] > 0.02

    # closed CAD can (R = 1, L = 2; wall strips, end disks as fans of slivers) against
    # the exact disk-to-disk and wall-to-wall exchange
    Rc, L, nphi = 1.0, 2.0, 90
    facets = NTuple{3,Vec3}[]
    φs = range(0, 2π; length = nphi + 1)
    for j in 1:nphi
        a, b = φs[j], φs[j+1]
        cad_quad!(facets, cad_pt(Rc, a, 0), cad_pt(Rc, b, 0), cad_pt(Rc, b, L), cad_pt(Rc, a, L),
                  -cad_pt(1, (a + b) / 2, 0))
        for (z, f) in ((0.0, 1), (L, -1))
            p, q, c = Vec3(0, 0, z), cad_pt(Rc, a, z), cad_pt(Rc, b, z)
            push!(facets, f * cross(q - p, c - p)[3] > 0 ? (p, q, c) : (p, c, q))
        end
    end
    craw = stl_roundtrip(facets)
    cel = remesh_mesh(craw; size = 0.25)
    X = 1 + (1 + (Rc / L)^2) / (Rc / L)^2
    dd = (X - sqrt(X^2 - 4)) / 2
    ww = 1 - π * Rc^2 * (1 - dd) / (π * Rc * L)
    errs = map((Target([Component(:can, cel, gold())]),
                Target(curved_components(Component(:can, cel, gold()); reference = craw)))) do t
        F = view_factors(t)
        A = t.elem_area
        wl = [e for e in 1:nelements(t) if abs(t.elem_normal[e][3]) < 0.5]
        bt = [e for e in 1:nelements(t) if t.elem_normal[e][3] > 0.5]
        tp = [e for e in 1:nelements(t) if t.elem_normal[e][3] < -0.5]
        Fdd = sum(A[i] * F[i, j] for i in bt, j in tp) / sum(A[bt])
        Fww = sum(A[i] * F[i, j] for i in wl, j in wl) / sum(A[wl])
        (; dd = Fdd / dd - 1, ww = Fww / ww - 1, closure = sum(abs, closure_error(F)) / nelements(t))
    end
    @test abs(errs[2].dd) < 3e-3 && abs(errs[2].ww) < 1e-3
    @test abs(errs[2].ww) < abs(errs[1].ww) && errs[2].closure < errs[1].closure

    # a CAD hohlraum (the cavity side of two thick halves, LEH apertures, a CAD
    # capsule), fitted, against the exact capsule rows (a sphere radiates by solid
    # angle: L/√(4R²+L²) of its power onto the wall) and against the primitive hohlraum
    # of the same dimensions
    R, Rout, Rleh, h, t, rc = 2.875, 2.975, 1.55, 4.7, 0.05, 0.9
    cavity(m) = split_mesh(m, c -> norm(c[1:2]) < R + 0.01 && abs(c[3]) < h + 1e-6)[1]
    comps = Component[]
    for (name, F) in ((:top, cad_half(; R, Rout, Rleh, h, t, nphi = 360)),
                      (:bottom, cad_mirror(cad_half(; R, Rout, Rleh, h, t, nphi = 360))),
                      (:capsule, cad_sphere(rc; ntheta = 56, nphi = 112)))
        m = orient_mesh(stl_roundtrip(F))
        append!(comps, curved_components(Component(name, cavity(remesh_mesh(m; size = 0.5)),
                                                   name === :capsule ? ablator() : gold()); reference = m))
    end
    push!(comps, aperture(:leh_top, disk_mesh(Rleh, h; nr = 3, nphi = 24, up = false); surface = DiskSurface(center = (0, 0, h))),
          aperture(:leh_bottom, disk_mesh(Rleh, -h; nr = 3, nphi = 24, up = true); surface = DiskSurface(center = (0, 0, -h))))
    tc = Target(comps)
    @test [c.name for c in tc.components] == Symbol.(["top/cyl1", "top/rest", "bottom/cyl1", "bottom/rest",
                                                      "capsule/sph1", "leh_top", "leh_bottom"])
    F = view_factors(tc)
    A = tc.elem_area
    cap = elements(tc, Symbol("capsule/sph1"))
    @test sum(A[cap]) ≈ 4π * rc^2 rtol = 1e-7
    wall = elements(tc, Symbol("top/cyl1"), Symbol("bottom/cyl1"))
    @test sum(A[wall]) ≈ 4π * R * h rtol = 1e-7
    caps = elements(tc, Symbol("top/rest"), Symbol("bottom/rest"))
    leh = elements(tc, :leh_top, :leh_bottom)
    row(J) = sum(A[i] * F[i, j] for i in cap, j in J) / sum(A[cap])
    @test row(wall) ≈ 2h / sqrt(4R^2 + 4h^2) rtol = 1e-3
    @test row(caps) ≈ h / sqrt(h^2 + Rleh^2) - h / sqrt(h^2 + R^2) rtol = 2e-3
    @test row(leh) ≈ 1 - h / sqrt(h^2 + Rleh^2) rtol = 2e-3
    # with the wall driven uniformly, the capsule's flux and the LEH loss match the
    # primitive hohlraum's (both closed by Sinkhorn balancing)
    prim = hohlraum(; length = 2h, radius = R, leh_radius = Rleh, capsule_radius = rc, h = 0.5, capsule_h = 0.15)
    drive(t, w) = (P = zeros(nelements(t)); P[w] .= t.elem_area[w]; P)
    sp = solve(RadiationModel(prim; closed = true), drive(prim, elements(prim, :wall)))
    Fc = copy(F)
    enforce_closure!(Fc, A)
    sc = solve(RadiationModel(tc; F = Fc), drive(tc, wall))
    flux(sol, t, es) = sum(sol.H[es] .* t.elem_area[es]) / sum(t.elem_area[es])
    @test flux(sc, tc, cap) ≈ flux(sp, prim, elements(prim, :capsule)) rtol = 3e-4
    @test flux(sc, tc, leh) ≈ flux(sp, prim, elements(prim, :leh_top, :leh_bottom)) rtol = 3e-4
end
