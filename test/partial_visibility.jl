# Partial visibility (`view_factors(...; visibility = :sampled)`): pairs that see each
# other past a rim get the sampled visible fraction instead of the centroid ray's 0 or 1.

using Briareus, LinearAlgebra, Test

# A hohlraum with 0.1 mm LEH bores and a 20 µm waist gap, closed by apertures: the
# primitive cavity of the "Radiation in an imported CAD hohlraum" tutorial at 0.5 mm
# (2,244 elements).
function bored_cavity(; hw = 0.5, hc = 0.15)
    R, Rout, Rleh, h, g2, rc = 2.875, 2.925, 1.55, 4.7, 0.01, 0.9
    zL = h + 0.1
    nphi = 8 * ceil(Int, 2π * R / hw / 8)
    nz, nr = ceil(Int, (h - g2) / hw), ceil(Int, (R - Rleh) / hw)
    (; capsule_ntheta, capsule_nphi) = capsule_counts(rc, hc)
    au = gold()
    wall(z0, z1) = cylinder_mesh(R, z0, z1; nz, nphi, inward = true)
    disk(z) = DiskSurface(center = (0, 0, z))
    return Target([Component(:wall_top, wall(g2, h), au; surface = CylinderSurface(R)),
                   Component(:wall_bottom, wall(-h, -g2), au; surface = CylinderSurface(R)),
                   Component(:cap_top, annulus_mesh(Rleh, R, h; nr, nphi, up = false), au; surface = disk(h)),
                   Component(:cap_bottom, annulus_mesh(Rleh, R, -h; nr, nphi, up = true), au; surface = disk(-h)),
                   Component(:bore_top, cylinder_mesh(Rleh, h, zL; nz = 1, nphi = nphi ÷ 2), au; surface = CylinderSurface(Rleh)),
                   Component(:bore_bottom, cylinder_mesh(Rleh, -zL, -h; nz = 1, nphi = nphi ÷ 2), au; surface = CylinderSurface(Rleh)),
                   Component(:gap_top, annulus_mesh(R, Rout, g2; nr = 1, nphi, up = false), au; surface = disk(g2)),
                   Component(:gap_bottom, annulus_mesh(R, Rout, -g2; nr = 1, nphi, up = true), au; surface = disk(-g2)),
                   Component(:capsule, sphere_mesh(rc; ntheta = capsule_ntheta, nphi = capsule_nphi), ablator();
                             surface = SphereSurface(rc)),
                   aperture(:leh_top, disk_mesh(Rleh, zL; nr = 3, nphi = 24, up = false); surface = disk(zL)),
                   aperture(:leh_bottom, disk_mesh(Rleh, -zL; nr = 3, nphi = 24, up = true); surface = disk(-zL)),
                   aperture(:gap, cylinder_mesh(Rout, -g2, g2; nz = 1, nphi = 180, inward = true);
                            surface = CylinderSurface(Rout))])
end

@testset "partial visibility" begin
    @test_throws ArgumentError view_factors(hohlraum(h = 2.0, capsule_h = 0.8); visibility = :rays)

    # nothing in the way: identical to the centroid rule, bit for bit
    s = Target(Component(:s, sphere_mesh(2.0; ntheta = 8, nphi = 12, outward = false), gold();
                         surface = SphereSurface(2.0)))
    @test view_factors(s; visibility = :sampled) == view_factors(s)

    # bored LEHs and a waist gap (closed by apertures): row sums before balancing
    t = bored_cavity()
    F0 = view_factors(t)
    F1 = view_factors(t; visibility = :sampled)
    A = t.elem_area
    @test maximum(abs, A .* F1 .- (A .* F1)') < 1e-15
    # pairs the probes find wholly visible or hidden keep the centroid factor exactly;
    # the others only shrink from the unobstructed factor (or appear, if the centroid
    # ray was blocked)
    changed = count(F1 .!= F0)
    @test 0 < changed < length(F0) ÷ 4
    rows(F, names) = (es = reduce(vcat, [elements(t, n) for n in names]); maximum(abs, closure_error(F)[es]))
    walls = (:wall_top, :wall_bottom)
    lehs = (:leh_top, :leh_bottom)
    # largest |row sum − 1|, measured: walls 9.7e-3 → 1.9e-3, LEH apertures 5.9e-2 →
    # 6.7e-3, capsule 1.6e-3 → 5.0e-4 (the gap's faces, 50 µm wide behind a 20 µm slot,
    # stay at 5%: the slot is narrower than a sample cell of the walls)
    @test rows(F1, walls) < 0.25rows(F0, walls)
    @test rows(F1, lehs) < 0.2rows(F0, lehs)
    @test rows(F1, (:capsule,)) < 0.5rows(F0, (:capsule,))
    # its facets crack where the 24-sided LEH disks meet the 20-sided bores and the
    # 180-sided gap aperture the 40-sided gap faces (open fraction measured 7.5e-4, 12%
    # on the worst bore element); without the LEH apertures the holes take ~7% of the view
    o = open_fraction(t)
    @test 0 < o.fraction < 2e-3 && maximum(o.element[elements(t, :bore_top, :bore_bottom)]) > 0.02
    @test open_fraction(Target(filter(c -> !startswith(String(c.name), "leh"), t.components))).fraction > 0.03

    # A pair half shadowed by a thin plate: two unit squares facing each other 2 apart,
    # and an opaque plate in the mid-plane covering x < 0. A segment from p to q crosses
    # the mid-plane at (p + q)/2, so (p, q) → (−p, −q) swaps shadowed and seen pairs of
    # equal kernel: half the unshadowed exchange (both squares are symmetric about their
    # centers; the upper one is turned 45° so that no two samples mirror each other and
    # tie on the edge). The edge is moved 1e-3 past the centroid ray, which is blocked.
    # Measured: 0.493 (g = 4), 0.500 (g = 8).
    RotZ45 = [cosd(45) -sind(45) 0; sind(45) cosd(45) 0; 0 0 1]
    lo = Component(:lo, rectangle_mesh(1.0, 1.0; nu = 1, nv = 1), gold())
    hi = Component(:hi, transform(rectangle_mesh(1.0, 1.0; nu = 1, nv = 1), rotation_to((0, 0, -1)) * RotZ45, Vec3(0, 0, 2)), gold())
    plate = translate(Component(:plate, rectangle_mesh(4.0, 4.0; nu = 1, nv = 1), gold(); two_sided = true),
                      (-1.999, 0, 1))
    open_ = view_factors(Target(lo, hi))
    shut = Target(lo, hi, plate)
    i, j = elements(shut, :lo)[1], elements(shut, :hi)[1]
    @test view_factors(shut)[i, j] == 0
    for g in (4, 8)
        F = view_factors(shut; visibility = :sampled, visibility_grid = g)
        @test isapprox(F[i, j] / open_[i, j], 0.5; atol = 0.01)
        @test F[i, j] ≈ F[j, i]                              # equal areas
    end
end

# Remeshed CAD (`cad_half`, `cad_mirror`, `cad_sphere` from cadfit.jl, `facet_soup` from
# meshops.jl): elements folded across many facets (48 per capsule element here). Their
# sample cells' mean points lay behind a convex part's facets, which blocked 91% of the
# capsule's sample rays: capsule row sums 0.03 (0.0–0.999) and flux 40× too low (ROADMAP
# open question 13). Samples now lie on the facets, and an element's own facets do not
# occlude its samples.
@testset "sampled visibility on remeshed CAD" begin
    R, Rout, Rleh, h, th, rc = 2.875, 2.975, 1.55, 4.7, 0.05, 0.9
    cavity(m) = split_mesh(m, c -> norm(c[1:2]) < R + 0.01 && abs(c[3]) < h + 1e-6)[1]
    part(name, F, size) = Component(name, cavity(remesh_mesh(orient_mesh(facet_soup(F)); size)),
                                    name === :capsule ? ablator() : gold())
    t = Target([part(:top, cad_half(; R, Rout, Rleh, h, t = th, nphi = 360), 0.6),
                part(:bottom, cad_mirror(cad_half(; R, Rout, Rleh, h, t = th, nphi = 360)), 0.6),
                part(:capsule, cad_sphere(rc; ntheta = 56, nphi = 112), 0.3),
                aperture(:leh_top, disk_mesh(Rleh, h; nr = 3, nphi = 24, up = false)),
                aperture(:leh_bottom, disk_mesh(Rleh, -h; nr = 3, nphi = 24, up = true))]; visibility = :sampled)
    cap = elements(t, :capsule)
    A = t.elem_area
    capmean(x) = sum(x[cap] .* A[cap]) / sum(A[cap])
    Fs, Fc = view_factors(t), view_factors(t; visibility = :centroid)
    rs, rsc = vec(sum(Fs; dims = 2)), vec(sum(Fc; dims = 2))
    # closed cavity (LEH apertures): measured capsule rows 0.9991 (0.9974–1.0005),
    # centroid rule 0.9986
    @test minimum(rs[cap]) > 0.995
    @test abs(capmean(rs) - capmean(rsc)) < 1e-3
    # a two-ring source on the walls: the capsule's mean flux against the Monte Carlo
    # radiosity on the same facets (the element model with exact factors), measured
    # −0.33% sampled, −0.45% centroid (σ 0.06%)
    P = zeros(nelements(t))
    for e in elements(t, :top, :bottom)
        p, n = t.elem_centroid[e], t.elem_normal[e]
        dot(n, Vec3(p[1], p[2], 0)) < 0 && (P[e] = exp(-((abs(p[3]) - 3.5) / 0.8)^2) * A[e])
    end
    Hs = capmean(solve(RadiationModel(t; F = Fs, closed = false), P).H)
    Hc = capmean(solve(RadiationModel(t; F = Fc, closed = false), P).H)
    if LONG_TESTS
        _, batch = Briareus.montecarlo_incident(t, P, cap; paths = 8_000, reemit = :element)
        m = batch * (A[cap] ./ sum(A[cap]))
        nb = length(m)
        Hm = sum(m) / nb
        se = sqrt(sum(abs2, m .- Hm) / (nb - 1) / nb) / Hm
        @test abs(Hs / Hm - 1) < 0.006
        @test abs(Hs / Hm - 1) < abs(Hc / Hm - 1) + 3se
    end
end
