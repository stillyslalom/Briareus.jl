# Unconverted-light deposition: exact geometric checks and a Monte Carlo ray trace.

@testset "unconverted-light deposition" begin
    bl = first(nif_beams())
    k = direction(bl)
    ô = Briareus._own_side(bl.port, Briareus._nif_quad_port(bl.quad))
    ŵ = cross(k, ô)
    fb = Vec3(bl.aim...) - bl.defocus * bl.port
    # the NIF colors absorbing all their light (the geometry alone)
    uc = [UnconvertedColor(c.name; footprint = c.footprint, edge = c.edge, focal_shift = c.focal_shift)
          for c in nif_unconverted()]

    # A plate in the 3ω focal plane facing the lens: the footprint squares are exact
    # (22 mm with its near edge 4.8 mm off the 3ω axis; 12 mm at 2.9 mm), uniform, and
    # take the whole power
    plate = Component(:plate, transform(rectangle_mesh(80.0, 80.0; nu = 40, nv = 40), hcat(ŵ, ô, -k),
                                        Vec3(0, 0, 0)), gold())
    t = Target(plate)
    d = unconverted_deposition(t, [bl]; ratio = (omega1 = 0.1, omega2 = 0.02), colors = uc)
    @test isempty(d.skipped)
    uv(x) = (dot(x - fb, ô), dot(x - fb, ŵ))
    tris = [findall(==(e), t.tri_front) for e in 1:nelements(t)]
    corners(e) = [uv(v) for i in tris[e] for v in (t.tri_v0[i], t.tri_v0[i] + t.tri_e1[i], t.tri_v0[i] + t.tri_e2[i])]
    for (col, r) in zip(uc, (0.1, 0.02))
        dc = getfield(d, col.name)
        P = r * bl.power
        @test dc.delivered == P && isapprox(sum(dc.element), P; rtol = 1e-12) && isapprox(dc.per_beam[1], P; rtol = 1e-12)
        h = col.footprint / 2
        c = col.edge + h
        inside(p) = abs(p[1] - c) <= h && abs(p[2]) <= h
        full = [e for e in 1:nelements(t) if all(inside, corners(e))]
        dark = [e for e in 1:nelements(t) if any(side -> all(side, corners(e)),
                                                  (p -> p[1] <= c - h, p -> p[1] >= c + h, p -> p[2] <= -h, p -> p[2] >= h))]
        @test !isempty(full) && all(e -> dc.element[e] ≈ P * t.elem_area[e] / col.footprint^2, full)
        @test all(e -> dc.element[e] <= 1e-12 * P, dark)                # edges on element edges: rounding
    end

    # a color's absorption scales what is absorbed, not what is delivered
    weak = UnconvertedColor(:omega1; footprint = 22.0, edge = 4.8, focal_shift = 460.0, absorption = 0.4)
    dw = unconverted_deposition(t, [bl]; ratio = (omega1 = 0.1,), colors = [weak]).omega1
    @test dw.delivered == 0.1 * bl.power && isapprox(sum(dw.element), 0.04 * bl.power; rtol = 1e-12)
    @test_throws ArgumentError UnconvertedColor(:omega1; footprint = 22.0, edge = 4.8, focal_shift = 460.0,
                                                absorption = 1.5)

    # A tilted two-sided plate 10 mm upstream, partly shadowed by a sphere 60 mm upstream:
    # against a Monte Carlo trace of the same beam (uniform rays through the footprint
    # square toward the 1ω focus, each to its first hit)
    col = uc[1]
    h0, f = col.footprint / 2, col.focal_shift
    c0 = col.edge + h0
    A = fb + f * k + c0 * (1 + f / 7700) * ô
    tilt = [1 0 0; 0 cosd(40) -sind(40); 0 sind(40) cosd(40)]
    plateB = Component(:plate, transform(rectangle_mesh(60.0, 60.0; nu = 30, nv = 30), hcat(ô, ŵ, k) * tilt,
                                         -10.0 * k), gold(); two_sided = true)
    ball = Component(:ball, translate(sphere_mesh(3.0; ntheta = 16, nphi = 32), fb + c0 * ô - 60.0 * k), gold())
    tB = Target(plateB, ball)
    dB = unconverted_deposition(tB, [bl]; ratio = Dict(:omega1 => 1.0), colors = uc[1:1]).omega1
    n = 1_000_000
    rng = Xoshiro(7)
    E = zeros(nelements(tB))
    for _ in 1:n
        q = fb + (c0 - h0 + 2h0 * rand(rng)) * ô + (-h0 + 2h0 * rand(rng)) * ŵ
        dir = normalize(A - q)
        _, id = Briareus.closest_hit(tB.bvh, q - 3000 * dir, dir)
        id == 0 && continue
        e = dot(tB.tri_normal[id], dir) < 0 ? tB.tri_front[id] : tB.tri_back[id]
        E[e] += bl.power / n
    end
    # (the shadow edge is resolved to `resolution`: conserved to ~1e-9)
    @test isapprox(sum(dB.element), bl.power; rtol = 1e-6) && isapprox(sum(E), bl.power; rtol = 1e-9)
    pb = sum(dB.element[elements(tB, :ball)])
    σb = bl.power * sqrt(pb / bl.power * (1 - pb / bl.power) / n)
    @test 0.04 < pb / bl.power < 0.05 && abs(pb - sum(E[elements(tB, :ball)])) < 4σb
    σ = sqrt.(dB.element * bl.power / n)
    lit = findall(>(1e-6), dB.element)
    @test all(i -> abs(dB.element[i] - E[i]) < 5σ[i], lit)
    @test all(i -> dB.element[i] > 1e-6 || E[i] == 0, eachindex(E))

    # A NIF quad is its four beamlets, each footprint on its own side of the quad axis
    q = only(b for b in nif_quads() if b.quad == bl.quad)
    wide = Component(:wide, transform(rectangle_mesh(120.0, 120.0; nu = 60, nv = 60), hcat(ŵ, ô, -k),
                                      Vec3(0, 0, 0)), gold())
    # (with the NIF colors: 1ω absorbs 0.88 of what reaches the plate)
    dq = unconverted_deposition(Target(wide), [q]; ratio = (omega1 = 0.05, omega2 = 0.0)).omega1
    @test nif_unconverted()[1].absorption == 0.88 && isapprox(sum(dq.element), 0.88 * 0.05 * q.power; rtol = 1e-12)
    tw = Target(wide)
    side = [dot(tw.elem_centroid[e] - fb, ô) for e in 1:nelements(tw)]
    lhs, rhs = sum(dq.element[side .< 0]), sum(dq.element[side .> 0])
    @test isapprox(lhs, rhs; rtol = 0.05) && isapprox(lhs + rhs, 0.88 * 0.05 * q.power; rtol = 1e-12)

    # beams that cannot be placed are listed; an offset places them
    b = Beam("b"; theta = 0, phi = 0, power = 1.0)
    ds = unconverted_deposition(tB, [b]; ratio = (omega1 = 0.1, omega2 = 0.1), colors = uc)
    @test ds.skipped == ["b"] && ds.omega1.delivered == 0 && sum(ds.omega1.element) == 0
    @test unconverted_deposition(tB, [b]; ratio = (omega1 = 0.1, omega2 = 0.1), colors = uc,
                                 offsets = Dict("b" => 0.0)).omega1.delivered ≈ 0.1
    @test_throws ArgumentError unconverted_deposition(tB, [bl]; ratio = (omega1 = 0.1,))
    @test_throws ArgumentError unconverted_deposition(tB, [bl]; ratio = (omega1 = -0.1, omega2 = 0.0))
    @test_throws DimensionMismatch unconverted_deposition(tB, [bl]; ratio = (omega1 = [0.1, 0.2], omega2 = 0.0))
end
