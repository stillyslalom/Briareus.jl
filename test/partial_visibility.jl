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
