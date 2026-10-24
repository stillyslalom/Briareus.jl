# Entrance-hole shields ("sombreros"): `hohlraum(; sombrero)`.

using Briareus, LinearAlgebra, Test

@testset "sombrero" begin
    kw = (; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9, h = 0.6, capsule_h = 0.3)
    # off by default: the same target, element for element
    t0 = hohlraum(; kw...)
    for s in (nothing, false)
        ts = hohlraum(; kw..., sombrero = s)
        @test [c.name for c in ts.components] == [c.name for c in t0.components]
        @test ts.elem_centroid == t0.elem_centroid && ts.elem_area == t0.elem_area
    end

    s = Sombrero()
    @test (s.inner_radius, s.outer_radius, s.offset) == (4.0, 9.2, 0.623)
    t = hohlraum(; kw..., sombrero = true)
    @test [c.name for c in t.components] == [[c.name for c in t0.components]; :sombrero_top; :sombrero_bottom]
    # the hohlraum's elements come first and are unchanged
    n0 = nelements(t0)
    @test t.elem_centroid[1:n0] == t0.elem_centroid && t.elem_area[1:n0] == t0.elem_area
    z = 10.1 / 2 - 0.623
    for (name, sgn) in ((:sombrero_top, 1), (:sombrero_bottom, -1))
        es = elements(t, name)
        @test length(es) == 4 * 96
        # exact annulus area (DiskSurface sectors), in the plane ±z, facing away from the waist
        @test sum(t.elem_area[es]) ≈ π * (9.2^2 - 4.0^2) rtol = 1e-12
        @test all(e -> t.elem_centroid[e][3] ≈ sgn * z && t.elem_normal[e] ≈ Vec3(0, 0, sgn), es)
        @test t.components[component_index(t, name)].kind === :surface       # opaque
        ρ = [hypot(v[1], v[2]) for v in t.components[component_index(t, name)].mesh.vertices]
        @test all(isapprox.(extrema(ρ), (4.0, 9.2); rtol = 1e-14))
    end

    # outside the cavity: the shields and the hohlraum do not see each other, so the
    # closed enclosure, its balancing and its solution are those without them
    m0, m = RadiationModel(t0), RadiationModel(t)
    sh = n0+1:nelements(t)
    @test all(iszero, m.F[sh, :]) && all(iszero, m.F[:, sh])
    @test m.F[1:n0, 1:n0] ≈ m0.F rtol = 1e-12
    P0 = zeros(n0)
    P0[elements(t0, :wall)] .= 1.0
    sol0, sol = solve(m0, P0), solve(m, [P0; zeros(length(sh))])
    @test sol.H[1:n0] ≈ sol0.H rtol = 1e-12
    @test all(iszero, sol.H[sh])

    # geometry checks
    @test_throws ArgumentError Sombrero(inner_radius = 9.2)
    @test_throws ArgumentError hohlraum(; kw..., sombrero = Sombrero(inner_radius = 2.8))   # cuts the wall
    @test hohlraum(; kw..., sombrero = Sombrero(inner_radius = 2.0, offset = -0.5)) isa Target  # beyond the hole
    @test_throws ArgumentError hohlraum(; kw..., sombrero = Sombrero(inner_radius = 1.5, offset = -0.5))
    @test_throws ArgumentError hohlraum(; kw..., sombrero = Sombrero(offset = 6.0))

    # it blocks beams: beams from ports on the shield's side, aimed at the top shield, land
    # on its front from 30° and on its waist side (the back-side loss) from 120°; operator
    # and ray trace agree
    spot = SuperGaussianSpot(0.3; order = 4.0)
    top = component_index(t, :sombrero_top)
    for (θ, front) in ((30, true), (120, false))
        b = Beam("s$θ"; theta = θ, phi = 180, power = 1.0, spot, fnumber = 20.0, aim = (-6.5, 0.0, z))
        d = deposit(DepositionOperator(t, [b]), [b])
        r = raytrace_deposition(t, [b]; nspot = 61)
        on = front ? sum(d.element[elements(t, :sombrero_top)]) : d.backside[top]
        # back sides are sampled only for their total (`backside_resolution`)
        @test on ≈ 1 rtol = (front ? 1e-4 : 2e-3)
        @test (front ? sum(r.element[elements(t, :sombrero_top)]) : r.backside[top]) ≈ 1 rtol = 1e-6
    end

    # clearance: shields are hardware like any other surface. Cone-averaged drive
    # pointing clears them (core 1e-4 envelope, no wings: > 1.5 mm worst case); a beam
    # aimed past the hole's edge onto the shield fails
    qs = nif_quads()
    rz = Dict(:inner23 => (1.09, 3.09), :inner30 => (0.6, 3.63), :outer44 => (0.0, 5.2), :outer50 => (0.0, 5.01))
    qs = [(b = qs[findfirst(b -> b.group === g, qs)]; φ = port_angles(b)[2]; (r0, z0) = rz[g];
           remake(b; aim = (-r0 * cosd(φ), -r0 * sind(φ), sign(b.port[3]) * z0)))
          for g in (:inner23, :inner30, :outer44, :outer50)]
    hw = (:sombrero_top, :sombrero_bottom)
    rep = clearance(t, qs; targets = (:wall, :capsule), hardware = hw, peak_fraction = 1e-4, wings = [],
                    nangles = 8)
    @test all(x -> x.clearance > 1.5 && x.status === :pass, rep.rows)
    b = remake(qs[end]; aim = (-6.0 * cosd(port_angles(qs[end])[2]), -6.0 * sind(port_angles(qs[end])[2]), z))
    bad = clearance(t, [b]; targets = (:wall, :capsule), hardware = hw, peak_fraction = 1e-4, wings = [], nangles = 8)
    @test minimum(x.clearance for x in bad.rows) < 0 && any(x -> x.status === :fail, bad.rows)
end
