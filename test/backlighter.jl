# Backlighter foils on a line of sight (`backlighter_foil`) and the line-of-sight
# obstruction check (`line_of_sight`), against closed-form geometry.

using Briareus, LinearAlgebra, Test

@testset "backlighter line of sight" begin
    # placement: on the far side of TCC, normal along the line of sight, width along θ̂
    for (θ, φ) in ((90.0, 78.0), (90.0, 348.0), (37.0, 211.0), (0.0, 0.0))
        d = port_vector(θ, φ)
        eθ, eφ = Briareus.port_basis(d)
        c = backlighter_foil((θ, φ); standoff = 8.0, width = 2.0, height = 1.5, n = 4)
        t = Target([c])
        fr = findall(==(1), t.elem_side)                   # front sides (the foil is two-sided)
        @test sum(t.elem_area[fr]) ≈ 3.0 rtol = 1e-12
        @test sum(t.elem_area[fr] .* t.elem_centroid[fr]) / 3.0 ≈ -8.0 * d atol = 1e-12
        @test all(e -> t.elem_normal[e] ≈ d, fr)
        vs = c.mesh.vertices
        @test maximum(v -> abs(dot(v, eθ)), vs) ≈ 1.0 && maximum(v -> abs(dot(v, eφ)), vs) ≈ 0.75
    end
    # rotated off the line of sight by β toward φ̂ (tilt_toward = 90°): the normal turns
    # by β, the center stays, and the foil presents height · cos β across the line of sight
    β = 35.0
    d = port_vector(90.0, 0.0)
    eθ, eφ = Briareus.port_basis(d)
    c = backlighter_foil((90.0, 0.0); standoff = 5.0, width = 2.0, height = 1.0, tilt = β,
                         tilt_toward = 90.0, n = 4)
    t = Target([c])
    fr = findall(==(1), t.elem_side)
    @test all(e -> t.elem_normal[e] ≈ cosd(β) * d + sind(β) * eφ, fr)
    @test sum(t.elem_area[fr] .* t.elem_centroid[fr]) / 2.0 ≈ -5.0 * d atol = 1e-12
    @test maximum(v -> abs(dot(v, eφ)), c.mesh.vertices) ≈ 0.5 * cosd(β)
    @test maximum(v -> abs(dot(v, eθ)), c.mesh.vertices) ≈ 1.0
    @test_throws ArgumentError backlighter_foil((90, 0); standoff = 0.0)

    # coverage of a foil rotated off the line of sight: a foil of width w turned by β
    # about its height axis covers a disk of radius a iff (w/2) cos β ≥ a (height ample).
    # The disk samples are the center and rings k = 1…4 of 8k points from θ̂, so the
    # uncovered samples are those with |a (k/4) cos ψ| > (w/2) cos β.
    β = 60.0
    c = backlighter_foil((90.0, 0.0); standoff = 8.0, width = 2.0, height = 4.0, tilt = β, n = 8)
    t = Target([c])
    edge = cosd(β)                                        # (w/2) cos β = 0.5 mm
    for (a, full) in ((edge * (1 - 1e-6), true), (edge * (1 + 1e-6), false), (0.9, false))
        r = line_of_sight(t, (90.0, 0.0); object = nothing, radius = a)
        n_out = count(abs(a * k / 4 * cos(2π * j / 8k)) > edge for k in 1:4 for j in 0:8k-1)
        @test r.covered ≈ 1 - n_out / 81
        @test (r.covered == 1) == full && r.clear == full
        @test r.samples == (9 * 17, 81) && r.rays == 9 * 17 * 81 && r.blocked == 0      # 8 × 16 patches
    end
    @test count(abs(edge * (1 + 1e-6) * cos(2π * j / 32)) > edge for j in 0:31) == 2   # ψ = 0, π
    # turned edge-on past the line of sight, the bundle would leave through the back
    c = backlighter_foil((90.0, 0.0); standoff = 8.0, width = 2.0, tilt = 100.0, n = 2)
    r = line_of_sight(Target([c]), (90.0, 0.0); object = nothing, radius = 0.5)
    @test !r.clear && r.covered == 0 && (:backlighter => r.rays) in r.blocked_by

    # A bundle exactly grazing a shield's rim. Line of sight +x, foil (half-width w) at
    # x = −D, disk of radius a at TCC, an annular shield in the plane z = Z facing +z
    # (bore 4 mm, rim radius ρ). The steepest rising ray, from the foil's lower edge
    # (−D, 0, −w) through the disk's top (0, 0, a), crosses z = Z at
    # x* = D (Z + w)/(a + w) − D; every other ray crosses it farther out. The shield
    # clears the bundle iff ρ < x*; just past that, it blocks that one ray (the rim
    # polygon has a vertex at φ = 0, and the ray crosses there).
    D, w, a, Z = 8.0, 1.0, 0.9, 3.0
    xs = D * (Z + w) / (a + w) - D
    fl = backlighter_foil((90.0, 0.0); standoff = D, width = 2w, n = 8)
    for (ρ, nblocked) in ((xs * (1 - 1e-7), 0), (xs * (1 + 1e-7), 1))
        shield = Component(:shield, annulus_mesh(4.0, ρ, Z; nr = 2, nphi = 96), BLACK)
        capsule = Component(:capsule, icosphere_mesh(a; n = 8), ablator())
        t = Target([fl, shield, capsule])
        r = line_of_sight(t, (90.0, 0.0); radius = a)
        @test r.blocked == nblocked && r.clear == (nblocked == 0)
        @test r.blocked_by == [:shield => nblocked]
        @test r.covered == 1
        # the imaging rays of a distant detector run parallel to the line of sight, at
        # most a above the axis: they clear the shield either way
        ri = line_of_sight(t, (90.0, 0.0); radius = a, detector = Inf)
        @test ri.clear && ri.rays == 81
    end
    # the default disk radius is the object's projected radius (an icosphere's vertices
    # lie on the sphere, one of them on the great circle normal to x)
    t = Target([fl, Component(:capsule, icosphere_mesh(a; n = 8), ablator())])
    @test line_of_sight(t, (90.0, 0.0)).radius ≈ a rtol = 1e-12
    @test_throws ArgumentError line_of_sight(t, (90.0, 0.0); object = nothing)

    # a point detector: the imaging rays converge on it, so a plate beside the line of
    # sight between the object and the detector at distance L blocks them iff it reaches
    # within the rays' envelope there, a (1 − s/L) at distance s. The rays spread to
    # a (1 + D/L) = 1.26 mm at the foil, so a 2 mm foil leaves the disk's rim uncovered
    L, s = 20.0, 10.0
    r = line_of_sight(Target([fl]), (90.0, 0.0); object = nothing, radius = a, detector = L)
    @test !r.clear && r.blocked == 0
    @test r.covered ≈ 1 - count(a * (k / 4) * (1 + D / L) * max(abs(cos(2π * j / 8k)), abs(sin(2π * j / 8k))) > w
                                for k in 1:4 for j in 0:8k-1) / 81
    fl = backlighter_foil((90.0, 0.0); standoff = D, width = 3.0, n = 6)
    reach = a * (1 - s / L)
    for (y0, ok) in ((reach * (1 + 1e-6), true), (reach * (1 - 1e-6), false))
        plate = Component(:plate, transform(rectangle_mesh(2.0, 4.0; nu = 1, nv = 1),
                                            Briareus.rotation_to(Vec3(1, 0, 0)), Vec3(s, y0 + 2.0, 0)),
                          BLACK; two_sided = true)
        t = Target([fl, plate])
        r = line_of_sight(t, (90.0, 0.0); object = nothing, radius = a, detector = L)
        @test r.clear == ok && r.covered == 1
        @test r.blocked == (ok ? 0 : 1)            # the ray through the disk's φ̂ extreme
    end
end
