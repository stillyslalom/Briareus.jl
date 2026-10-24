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

# conical frustum from radius r0 at z0 to r1 at z1 as CAD tessellates it: nz bands of
# nphi strips (one band: vertices on the two rims only), facing the axis or away
function cad_frustum(r0, z0, r1, z1; nz = 1, nphi, inward = true)
    F = NTuple{3,Vec3}[]
    φs = range(0, 2π; length = nphi + 1)
    for i in 1:nz, j in 1:nphi
        (ra, za), (rb, zb) = ((r0 + (r1 - r0) * x, z0 + (z1 - z0) * x) for x in ((i - 1) / nz, i / nz))
        a, b = φs[j], φs[j+1]
        out = Vec3((z1 - z0) * cos((a + b) / 2), (z1 - z0) * sin((a + b) / 2), -(r1 - r0))
        cad_quad!(F, cad_pt(ra, a, za), cad_pt(ra, b, za), cad_pt(rb, b, zb), cad_pt(rb, a, zb), inward ? -out : out)
    end
    return F
end
# the fan of slivers CAD makes of a disk of radius r at height z, facing ±z
function cad_fan(r, z, up; nphi)
    F = NTuple{3,Vec3}[]
    φs = range(0, 2π; length = nphi + 1)
    for j in 1:nphi
        p, q, c = Vec3(0, 0, z), cad_pt(r, φs[j], z), cad_pt(r, φs[j+1], z)
        push!(F, (cross(q - p, c - p)[3] > 0) == up ? (p, q, c) : (p, c, q))
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
    @test_throws ArgumentError fit_surfaces(para; kinds = (:torus,))

    # cones: an entrance-hole cone (39°) and a 45° one as single bands (vertices on two
    # rims only, which a sphere fits as well: the cone is taken), a slightly tapered bore
    # (2.9°, apex 40 mm away) and a near-flat 83° cone in several bands, facing in or out
    for (r0, z0, r1, z1, nz, nphi, inward) in ((1.55, 4.7, 2.2, 5.5, 1, 360, true), (1.0, 0.0, 2.0, 1.0, 1, 120, false),
                                               (2.0, 0.0, 1.9, 2.0, 8, 180, true), (0.5, 0.0, 3.0, 0.3, 3, 120, true))
        k = (r1 - r0) / (z1 - z0)
        raw = orient_mesh(stl_roundtrip(cad_frustum(r0, z0, r1, z1; nz, nphi, inward)))
        for mesh in (raw, remesh_mesh(raw; size = 0.3))
            regs = fit_surfaces(mesh; reference = mesh === raw ? nothing : raw)
            @test length(regs) == 1 && regs[1].kind === :cone && length(regs[1].patches) == npatches(mesh)
            c = regs[1].surface
            @test c.half_angle ≈ atand(abs(k)) rtol = 1e-6
            @test c.apex ≈ Vec3(0, 0, z0 - r0 / k) atol = 1e-5 * max(1, abs(z0 - r0 / k))
            @test c.axis ≈ Vec3(0, 0, sign(k)) atol = 1e-7
            @test regs[1].facing == (inward ? -1 : 1)
        end
    end
end

# a torus section (tube radius r about a circle of radius R at height zc, tube angles
# ψ0..ψ1 in nψ bands of nphi strips), as CAD tessellates a fillet
function cad_torus(R, r, zc, ψ0, ψ1; nψ, nphi)
    F = NTuple{3,Vec3}[]
    φs = range(0, 2π; length = nphi + 1)
    ψs = range(ψ0, ψ1; length = nψ + 1)
    p(ψ, φ) = cad_pt(R + r * cos(ψ), φ, zc + r * sin(ψ))
    for i in 1:nψ, j in 1:nphi
        ψ, φ = (ψs[i] + ψs[i+1]) / 2, (φs[j] + φs[j+1]) / 2
        cad_quad!(F, p(ψs[i], φs[j]), p(ψs[i], φs[j+1]), p(ψs[i+1], φs[j+1]), p(ψs[i+1], φs[j]), cad_pt(cos(ψ), φ, sin(ψ)))
    end
    return F
end

@testset "fit_surfaces on faces of none of the kinds" begin
    opts = (; rtol = 1e-3, atol = 1e-4, cosmax = cosd(30.0))
    fit(m, kinds; kw...) = Briareus._fit_surfaces(m, kinds, opts, 6, 0.0, 1e-9; kw...)
    # a 0.15° tapered tube (r = 50 µm, 2 mm long, like a capsule's fill tube) without
    # cones: every patch seeds a cylinder growth of ~1000 patches that the rms test
    # rejects, which made the cost quadratic (4× the work per doubling of the facets,
    # 171 s at 18k facets); covered seeds now give up their growths early
    taper(n) = facet_soup(cad_frustum(0.05, 0.0, 0.05 - 2tand(0.15), 2.0; nz = 4n, nphi = n))
    m1, m2 = taper(34), taper(48)
    @test 1.9 < length(m2.triangles) / length(m1.triangles) < 2.1
    (r1, s1), (r2, s2) = fit(m1, (:cylinder, :sphere)), fit(m2, (:cylinder, :sphere))
    @test isempty(r1) && isempty(r2)
    @test s2.work < 3s1.work
    @test s1.givenup > 0.9npatches(m1) && s2.givenup > 0.9npatches(m2)
    @test isempty(fit_surfaces(m1; kinds = (:cylinder, :sphere)))
    # a true cone (half-angle 70°, two bands) 10° off a fine torus fillet that ends on a
    # flat cap, all four kinds: the fillet's growths are covered and given up (here with
    # covers forced small) and the cone is still found whole
    n = 56
    F = cad_torus(2.0, 0.5, 0.0, deg2rad(-60), π / 2; nψ = n ÷ 2, nphi = 4n)
    append!(F, cad_frustum(2.0, 0.5, 0.3, 0.5; nphi = 4n))
    ρ0, z0 = 2.0 + 0.5cosd(-60), 0.5sind(-60)
    append!(F, cad_frustum(ρ0, z0, ρ0 - 0.1sind(70), z0 - 0.1cosd(70); nz = 2, nphi = 4n))
    m = facet_soup(F)
    regs, st = fit(m, (:plane, :cylinder, :cone, :sphere); covermin = 16, check = 4)
    @test st.givenup > 1000
    cone = only(r for r in regs if r.kind === :cone && r.area > 0.5)
    @test cone.patches == collect(npatches(m)-16n+1:npatches(m))
    @test cone.surface.half_angle ≈ 70 rtol = 1e-6
    @test cone.surface.apex ≈ Vec3(0, 0, z0 - ρ0 / tand(70)) atol = 1e-6
    @test any(r -> r.kind === :plane && length(r.patches) == 8n, regs)     # the cap
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
    # dropping gently curved and small regions into the rest: the part's extent is
    # 9.7 mm, its area 215 mm², the bore's 0.49 mm²
    @test [x.name for x in curved_components(c; reference = raw, maxradius = 0.2)] == Symbol.(["half/cyl1", "half/rest"])
    @test [x.surface.radius for x in curved_components(c; reference = raw, minfraction = 0.01)[1:2]] ≈ [Rout, R] rtol = 1e-6
    @test length(curved_components(c; reference = raw, minfraction = 0.01)) == 3
    # by default a curved region of under 1e-3 of the part's area stays flat: a sphere
    # (50.3 mm²) with a small tube stub (0.031 mm², 6e-4 of it)
    stub = orient_mesh(stl_roundtrip([cad_sphere(2.0; ntheta = 24, nphi = 48);
                                      cad_frustum(0.05, 5.0, 0.05, 5.1; nphi = 24, inward = false)]))
    cs = Component(:s, stub, gold())
    @test [x.name for x in curved_components(cs)] == Symbol.(["s/sph1", "s/rest"])
    @test [x.name for x in curved_components(cs; minfraction = 0)] == Symbol.(["s/sph1", "s/cyl1"])
    # a cone, with the cylinders of a part
    cone = orient_mesh(stl_roundtrip([cad_frustum(1.0, 0.0, 2.0, 1.0; nphi = 120); cad_frustum(2.0, 1.0, 2.0, 2.5; nphi = 120)]))
    cs = curved_components(Component(:shield, remesh_mesh(cone; size = 0.3), gold()); reference = cone)
    @test [x.name for x in cs] == Symbol.(["shield/cyl1", "shield/cone1"])
    @test cs[2].surface isa ConeSurface
    @test cs[2].surface.half_angle ≈ 45 rtol = 1e-7
    @test sum(Target(cs).elem_area) ≈ π * 3 * sqrt(2) + 6π rtol = 1e-7
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
    # (the capsule is one sphere: it keeps its name)
    @test [c.name for c in tc.components] == Symbol.(["top/cyl1", "top/rest", "bottom/cyl1", "bottom/rest",
                                                      "capsule", "leh_top", "leh_bottom"])
    @test parts(tc) == [:top, :bottom, :capsule, :leh_top, :leh_bottom]
    F = view_factors(tc)
    A = tc.elem_area
    cap = elements(tc, :capsule)
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

    # the halves are split for computation only: the name-taking functions and the
    # reports see the parts, and the pieces stay reachable by their own names
    top, toprest = Symbol("top/cyl1"), Symbol("top/rest")
    @test elements(tc, :top) == elements(tc, top, toprest)
    @test part_pieces(tc, :top) == [1, 2] && part_pieces(tc, toprest) == [2]
    @test all(c -> c.part === :top, tc.components[1:2]) && tc.components[1].name === top
    @test_throws KeyError elements(tc, :nothing)
    @test occursin("(2 pieces)", sprint(show, MIME"text/plain"(), tc))
    @test_throws ArgumentError Target([Component(:a, tc.components[1].mesh, gold(); part = :b),
                                       Component(:b, tc.components[2].mesh, gold(); part = :c)])
    bal, balc = energy_balance(sc), energy_balance(sc; by = :component)
    @test sort(collect(keys(bal.absorbed))) == [:bottom, :capsule, :top]
    @test bal.absorbed[:top] ≈ balc.absorbed[top] + balc.absorbed[toprest] rtol = 1e-12
    @test sum(values(bal.absorbed)) ≈ sum(values(balc.absorbed)) rtol = 1e-12
    beams = nif_quads()
    dep = deposit(DepositionOperator(tc, beams), beams)
    cp, cpc = component_power(dep, tc), component_power(dep, tc; by = :component)
    @test !any(k -> occursin('/', String(k)), keys(cp))
    @test cp[:top] ≈ cpc[top] + cpc[toprest] rtol = 1e-12
    @test sum(values(cp)) ≈ sum(values(cpc)) rtol = 1e-12
    dir = mktempdir()
    write_results(dir, tc, dep; sol = sc, observe = :capsule)
    @test [first(split(l, ',')) for l in readlines(joinpath(dir, "components.csv"))[2:end]] ==
          ["top", "bottom", "capsule", "leh_top", "leh_bottom"]
    @test split(readline(joinpath(dir, "elements.csv")), ',')[1:4] == ["element", "component", "part", "side"]
    @test first.(read_stl_solids(write_stl(joinpath(dir, "t.stl"), tc; binary = false))) == String.(parts(tc))
    # observe, the direct penalty and the keep-outs take the parts
    sn = Scenario(tc, beams; observe = :top, radiation = false)
    @test sn.observed == elements(tc, :top)
    ev = evaluate(sn, initial_parameters(sn.param))
    @test power_fraction(ev, :top) > 0
    @test Briareus._penalties(ev; clip = 0, missed = 0, direct = 1.0, direct_names = (:top,)) ≈
          power_fraction(ev, :top)^2
    @test SphereKeepOut(tc).names == [:capsule]
    @test ApertureKeepOut(tc).names == [:leh_top, :leh_bottom]
    tr = talis_check(tc, beams[1:4]; capture = (:top, :bottom), omega1 = false)
    @test all(r -> r.component in (:top, :bottom, :capsule, :none), tr.rows)
    @test any(r -> r.rule === :clearance_3w && r.component === :capsule, tr.rows)
end

@testset "radiation in a fitted CAD cone" begin
    # a CAD frustum (45°, r 1 → 2 over z 0 → 1, one band of 180 strips) closed by two
    # fans, remeshed at 0.3 mm: the fitted cone against the exact exchanges (the disks'
    # coaxial factor; the cone's by closure and reciprocity), and against flat elements
    r1, r2, L = 1.0, 2.0, 1.0
    Ac, A1, A2 = π * (r1 + r2) * hypot(r2 - r1, L), π * r1^2, π * r2^2
    X = 1 + (1 + (r2 / L)^2) / (r1 / L)^2
    f12 = (X - sqrt(X^2 - 4(r2 / r1)^2)) / 2
    cd1, cd2 = A1 * (1 - f12) / Ac, (A2 - A1 * f12) / Ac
    raw = [stl_roundtrip(f) for f in (cad_frustum(r1, 0.0, r2, L; nphi = 180), cad_fan(r1, 0.0, true; nphi = 180),
                                      cad_fan(r2, L, false; nphi = 180))]
    el = [remesh_mesh(m; size = 0.3) for m in raw]
    disks = [Component(:d1, el[2], gold()), Component(:d2, el[3], gold())]
    errs = map((Target([Component(:cone, el[1], gold()); disks]),
                Target([curved_components(Component(:cone, el[1], gold()); reference = raw[1]); disks]))) do t
        F = view_factors(t)
        A = t.elem_area
        c = [e for e in 1:nelements(t) if abs(t.elem_normal[e][3]) < 0.9]
        d1, d2 = elements(t, :d1), elements(t, :d2)
        flux(I, J) = sum(A[i] * F[i, j] for i in I, j in J) / sum(A[I])
        (; names = [x.name for x in t.components], area = sum(A[c]) / Ac - 1, d2 = flux(c, d2) / cd2 - 1,
         cc = flux(c, c) / (1 - cd1 - cd2) - 1, back = flux(d2, c) / (1 - A1 * f12 / A2) - 1)
    end
    flat, fitted = errs
    @test fitted.names == [:cone, :d1, :d2]             # one cone: the part keeps its name
    @test abs(fitted.area) < 1e-8                         # (Float32 vertices)
    @test abs(fitted.d2) < 1e-3 && abs(fitted.cc) < 1e-3 && abs(fitted.back) < 1e-3
    @test abs(fitted.d2) < abs(flat.d2) / 3 && abs(fitted.back) < abs(flat.back) / 3
end

@testset "view factors on remeshed CAD elements" begin
    # a rectangle o + [0,1]u + [0,1]v as CAD strips (long along v), facing `facing`
    function strips(o, u, v, facing; n = 40)
        F = NTuple{3,Vec3}[]
        for i in 0:n-1
            a, b = o + u * i / n, o + u * (i + 1) / n
            cad_quad!(F, a, b, b + v, a + v, facing)
        end
        return facet_soup(F)
    end
    X, Y, Z = Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)
    flux(t, F, I, J) = sum(t.elem_area[i] * F[i, j] for i in I, j in J) / sum(t.elem_area[I])
    # perpendicular rectangles with a common edge of length 1, widths 1 and 0.6
    function hottel(W, H)
        Q = W^2 + H^2
        return (W * atan(1 / W) + H * atan(1 / H) - sqrt(Q) * atan(1 / sqrt(Q)) +
                log((1 + W^2) * (1 + H^2) / (1 + Q) * (W^2 * (1 + Q) / ((1 + W^2) * Q))^(W^2) *
                    (H^2 * (1 + Q) / ((1 + H^2) * Q))^(H^2)) / 4) / (π * W)
    end
    plates = [(strips(0X, Y, X, Z), strips(0X, 0.6Z, X, Y)),        # strips along the common edge
              (strips(0X, X, Y, Z), strips(0X, X, 0.6Z, Y))]        # and across it
    # remeshed planar strips are rectangles cut on common planes (two to four facets, the
    # extra ones at neighbors' corners): they take the two-facet path, and match Hottel
    # as primitive quads of the same size do (1.2e-3 for 4 × 4 quads, 1.9e-3 for 5 × 5;
    # the many-facet path, at a finer level, gave 4e-4 to 7e-4 on the staggered elements
    # remeshing used to make)
    for (a, b) in plates
        t = Target(Component(:a, remesh_mesh(a; size = 0.25), gold()), Component(:b, remesh_mesh(b; size = 0.25), gold()))
        @test maximum(length, t.elem_tris) <= 4
        F = view_factors(t)
        @test flux(t, F, elements(t, :a), elements(t, :b)) ≈ hottel(1.0, 0.6) rtol = 2e-3
        @test flux(t, F, elements(t, :b), elements(t, :a)) ≈ hottel(0.6, 1.0) rtol = 2e-3
    end
    # coaxial disks (radii 1 and 0.8, 0.4 apart) as fans of 360 slivers
    function fan(r, z, facing)
        F = NTuple{3,Vec3}[]
        φs = range(0, 2π; length = 361)
        for j in 1:360
            p, q, c = Vec3(0, 0, z), cad_pt(r, φs[j], z), cad_pt(r, φs[j+1], z)
            push!(F, dot(cross(q - p, c - p), facing) > 0 ? (p, q, c) : (p, c, q))
        end
        return facet_soup(F)
    end
    t = Target(Component(:d1, remesh_mesh(fan(1.0, 0.0, Z); size = 0.25), gold()),
               Component(:d2, remesh_mesh(fan(0.8, 0.4, -Z); size = 0.25), gold()))
    R1, R2 = 1.0 / 0.4, 0.8 / 0.4
    S = 1 + (1 + R2^2) / R1^2
    @test flux(t, view_factors(t), elements(t, :d1), elements(t, :d2)) ≈ (S - sqrt(S^2 - 4(R2 / R1)^2)) / 2 rtol = 4e-3
    # a closed box of strips: rows sum to one
    faces = [(0X, X, Y), (Z, Y, X), (0X, Z, X), (Y, X, Z), (0X, Y, Z), (X, Z, Y)]
    box = Target([Component(Symbol(:f, i), remesh_mesh(strips(o, u, v, cross(u, v)); size = 0.25), gold())
                  for (i, (o, u, v)) in enumerate(faces)])
    @test all(e -> dot(box.elem_normal[e], Vec3(0.5, 0.5, 0.5) - box.elem_centroid[e]) > 0, 1:nelements(box))
    F = view_factors(box)
    @test maximum(abs, closure_error(F)) < 6e-3
    @test sum(abs, closure_error(F)) / nelements(box) < 4e-3

    # The outline receiver is the sum over the facets to round-off, at points anywhere
    # (seen whole, cut by the point's tangent plane, behind), for planar elements and
    # elements folded across a faceted cylinder
    R = 2.875
    wall = orient_mesh(stl_roundtrip(cad_half(; R, Rout = 2.975, Rleh = 1.55, h = 4.7, t = 0.05, nphi = 120)))
    t = Target(Component(:p, remesh_mesh(plates[2][1]; size = 0.25), gold()),
               Component(:w, remesh_mesh(wall; size = 0.5), gold()))
    offset = cumsum([0; [length(c.mesh.triangles) for c in t.components]])
    rng = Xoshiro(7)
    worst, npoly, nfolded, nedges, nfacet_edges = 0.0, 0, 0, 0, 0
    for e in 1:nelements(t)
        length(t.elem_tris[e]) > 2 || continue
        r = Briareus._poly_receiver(t, e, offset)
        npoly += 1
        nfolded += length(r.pieces) > 1
        nedges += length(r.edges)
        nfacet_edges += 3 * length(t.elem_tris[e])
        tris = [(t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k]) for k in t.elem_tris[e]]
        for j in 1:nelements(t)
            # points on nearby elements with their normals, and points around the element
            j != e && (norm(t.elem_centroid[j] - r.center) < 3r.radius || rand(rng) < 0.002) || continue
            k = rand(rng, t.elem_tris[j])
            u, v = rand(rng), rand(rng)
            u + v > 1 && ((u, v) = (1 - u, 1 - v))
            for (p, n) in ((t.tri_v0[k] + u * t.tri_e1[k] + v * t.tri_e2[k], t.elem_normal[j]),
                           (r.center + 2r.radius * randn(rng, Vec3), normalize(randn(rng, Vec3))))
                worst = max(worst, abs(Briareus._receiver_vf(p, n, r) - sum(_point_polygon_vf(p, n, x) for x in tris)))
            end
        end
    end
    @test worst < 1e-12
    @test nfolded > 100 && npoly - nfolded > 100
    # cost: the outlines have under 40% of the facets' edges, and the emitter's root cells
    # are about as many as a square two-facet element has (2·4^depth, plus the cells its
    # edges cut), however many facets an element has
    @test nedges < 0.4 * nfacet_edges
    square = Target(Component(:q, rectangle_mesh(0.3, 0.3; nu = 1, nv = 1), gold()))
    ncells(t, e) = length(Briareus._cell_emitter(t, e, 3).cells)
    multi = filter(e -> length(t.elem_tris[e]) > 2, 1:nelements(t))
    @test all(e -> sum(c[1] for c in Briareus._cell_emitter(t, e, 3).cells) ≈ t.elem_area[e], multi)
    @test maximum(e -> length(t.elem_tris[e]), multi) > 20
    @test sum(e -> ncells(t, e), multi) / length(multi) < ncells(square, 1)
    @test maximum(e -> ncells(t, e), multi) < 1.5 * ncells(square, 1)
end

@testset "cad_target" begin
    # a multi-solid STL assembly with mixed windings, as CAD exports it: two thick-walled
    # hohlraum halves 40 µm apart at the waist, a capsule and a mounting block outside,
    # and primitive entrance-hole apertures passed with the parts as components
    R, Rout, Rleh, h, t, g, rc = 2.875, 2.975, 1.55, 4.7, 0.05, 0.02, 0.9
    top = [map(v -> v + Vec3(0, 0, g), f) for f in cad_half(; R, Rout, Rleh, h, t, nphi = 72)]
    F = NTuple{3,Vec3}[]
    lo, hi = Vec3(Rout + 0.3, -1, -1.5), Vec3(Rout + 2.3, 1, 1.5)
    c = [Vec3(x, y, z) for x in (lo[1], hi[1]), y in (lo[2], hi[2]), z in (lo[3], hi[3])]
    for q in ((1, 3, 7, 5), (2, 6, 8, 4), (1, 5, 6, 2), (3, 4, 8, 7), (1, 2, 4, 3), (5, 7, 8, 6))
        cad_quad!(F, c[q[1]], c[q[2]], c[q[3]], c[q[4]], sum(c[collect(q)]) / 4 - (lo + hi) / 2)
    end
    rng = Xoshiro(7)
    soup(F) = facet_soup([rand(rng) < 0.3 ? (a, c, b) : (a, b, c) for (a, b, c) in shuffle(rng, F)])
    names = ["hohlraum_top", "hohlraum_bottom", "capsule", "mount"]
    path = write_stl(joinpath(mktempdir(), "assembly.stl"),
                     [n => soup(f) for (n, f) in zip(names, (top, cad_mirror(top), cad_sphere(rc; ntheta = 16, nphi = 24), F))];
                     binary = false)
    solids = read_stl_solids(path)
    @test first.(solids) == names
    zl = h + t + g
    leh = [aperture(:leh_top, disk_mesh(Rleh, zl; nr = 2, nphi = 24, up = false)),
           aperture(:leh_bottom, disk_mesh(Rleh, -zl; nr = 2, nphi = 24, up = true))]
    mat = n -> n === :capsule ? ablator() : gold()
    size_of = n -> n === :capsule ? 0.2 : 0.5
    tc = cad_target([solids; leh]; inside = (0, 0, 0), remesh = size_of, fit = true, material = mat)

    # the cavity: inner walls, bores, caps and waist faces of the halves, the capsule and
    # the apertures as given; the outer walls (no cylinder of radius Rout) and the mount
    # are dropped
    @test !any(c -> startswith(String(c.name), "mount"), tc.components)
    @test tc.components[end-1:end] == leh
    cyl = [c for c in tc.components if c.surface isa CylinderSurface]
    @test sort([c.surface.radius for c in cyl]) ≈ [Rleh, Rleh, R, R] rtol = 1e-6
    sph = only(c for c in tc.components if c.surface isa SphereSurface)
    @test sph.name === :capsule && sph.part === :capsule
    @test sph.surface.radius ≈ rc rtol = 1e-6
    @test sph.material == ablator() && all(c -> c.material == gold(), cyl)
    # the curved elements have the true surfaces' areas
    area(c) = sum(tc.elem_area[elements(tc, c.name)])
    @test all(isapprox(area(c), 2π * R * h; rtol = 1e-6) for c in cyl if c.surface.radius > 2)
    @test area(sph) ≈ 4π * rc^2 rtol = 1e-6
    # compact elements of about the size asked
    @test 0.4 < sqrt(area(cyl[1]) / npatches(cyl[1].mesh)) < 0.6

    # the same steps by hand: orient, select on the facets, remesh the kept ones, select
    # again on the elements, fit on the kept facets
    kept = enclosure([[Component(n, orient_mesh(m; toward = (0, 0, 0)), mat(Symbol(n))) for (n, m) in solids]; leh]).target.components
    el = enclosure([c.kind === :aperture ? c : Component(c.name, remesh_mesh(c.mesh; size = size_of(c.name)), c.material)
                    for c in kept]).target.components
    ref = Dict(c.name => c.mesh for c in kept)
    hand = Target(reduce(vcat, [c.kind === :aperture ? [c] : curved_components(c; reference = ref[c.name]) for c in el]))
    @test [c.name for c in tc.components] == [c.name for c in hand.components]
    @test tc.elem_area == hand.elem_area && tc.elem_centroid == hand.elem_centroid

    # without `inside`, every part, oriented, in order; per-part settings by dictionary
    t0 = cad_target(solids; two_sided = Dict(:mount => true))
    @test [c.name for c in t0.components] == Symbol.(names)
    @test all(npatches(c.mesh) == npatches(m) for (c, (_, m)) in zip(t0.components, solids))
    @test [c.two_sided for c in t0.components] == [false, false, false, true]
    @test all(c -> c.material == gold(), t0.components)
    @test t0.components[3].mesh.triangles == orient_mesh(solids[3][2]; toward = (0, 0, 0)).triangles
    @test cad_target(solids; orient = false).components[1].mesh === solids[1][2]
    # coarsening instead of remeshing, the capsule only (a dictionary keyed by String)
    tk = cad_target(solids; inside = (0, 0, 0), coarsen = Dict("capsule" => 0.3), material = mat)
    tf = cad_target(solids; inside = (0, 0, 0), material = mat)
    cap(t) = t.components[findfirst(c -> c.name === :capsule, t.components)]
    @test npatches(cap(tk).mesh) < npatches(cap(tf).mesh)
    @test sum(tk.elem_area) ≈ sum(tf.elem_area) rtol = 1e-10
    @test [npatches(c.mesh) for c in tk.components if c.name !== :capsule] ==
          [npatches(c.mesh) for c in tf.components if c.name !== :capsule]
    @test_throws ArgumentError cad_target(solids; remesh = 0.5, coarsen = 0.5)
    @test_throws ArgumentError cad_target(solids; material = Dict("capsule" => ablator()))
    @test_throws ArgumentError cad_target([solids[1][2]])
    @test_throws ArgumentError cad_target(Pair{String,Mesh}[])
end
