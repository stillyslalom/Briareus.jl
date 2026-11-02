# Unconverted light over time: irradiance histories against closed forms, against
# `unconverted_deposition` at each time, and threshold crossings against hand computations.

@testset "unconverted-light histories" begin
    bl = first(nif_beams())                       # one beamlet: its whole power in one square
    k = direction(bl)
    ô = Briareus._own_side(bl.port, Briareus._nif_quad_port(bl.quad))
    ŵ = cross(k, ô)
    fb = Vec3(bl.aim...) - bl.defocus * bl.port
    # A plate in the 3ω focal plane facing the lens (2 mm elements)
    plate = Component(:plate, transform(rectangle_mesh(80.0, 80.0; nu = 40, nv = 40), hcat(ŵ, ô, -k),
                                        Vec3(0, 0, 0)), gold())
    t = Target(plate)
    uv(x) = (dot(x - fb, ô), dot(x - fb, ŵ))
    tris = [findall(==(e), t.tri_front) for e in 1:nelements(t)]
    corners(e) = [uv(v) for i in tris[e] for v in (t.tri_v0[i], t.tri_v0[i] + t.tri_e1[i], t.tri_v0[i] + t.tri_e2[i])]
    inside(col) = p -> abs(p[1] - col.edge - col.footprint / 2) <= col.footprint / 2 && abs(p[2]) <= col.footprint / 2
    uc = nif_unconverted()
    full1 = [e for e in 1:nelements(t) if all(inside(uc[1]), corners(e))]
    full2 = [e for e in 1:nelements(t) if all(inside(uc[2]), corners(e))]

    # A shaped pulse (a foot and a peak), a 1ω fraction that falls as the power rises
    # and a constant 2ω fraction
    pulse = PulseShape([0.0, 0.5, 2.0, 2.5, 4.0, 4.3], [0.0, 0.2, 0.2, 2.0, 2.0, 0.0])
    f1 = PulseShape([0.0, 5.0], [0.12, 0.04])
    times = collect(-0.5:0.1:5.0)
    h = unconverted_history(t, [bl], pulse, times; ratio = (omega1 = f1, omega2 = 0.01))
    @test isempty(h.skipped) && h.colors == [:omega1, :omega2]

    # closed form: inside each square the incident irradiance is the color's power over
    # the square's area, 1e14 P/footprint² W/cm² (mm² → cm²), times the pulse at each time
    I1 = [1e14 * f1(τ) * pulse(τ) / 22.0^2 for τ in times]
    I2 = [1e14 * 0.01 * pulse(τ) / 12.0^2 for τ in times]
    @test !isempty(full1) && !isempty(full2)
    @test isapprox(h.element[:omega1][full1, :], repeat(I1', length(full1)); rtol = 1e-12)
    @test isapprox(h.element[:omega2][full2, :], repeat(I2', length(full2)); rtol = 1e-12)
    # the incident power in total: the whole color, at every time
    P1 = vec(sum(h.element[:omega1] .* t.elem_area; dims = 1)) ./ 1e14
    @test isapprox(P1, [f1(τ) * pulse(τ) for τ in times]; rtol = 1e-12)
    @test all(iszero, h.element[:omega1][:, times .< 0]) && all(iszero, h.backside[:omega1])

    # a constant fraction: absorption × irradiance × area is `unconverted_deposition` of
    # the beams at each time's power, for both colors, from a quad (four beamlets) too
    q = only(b for b in nif_quads() if b.quad == bl.quad)
    ratio = (omega1 = 0.05, omega2 = 0.02)
    hq = unconverted_history(t, [q], pulse, times; ratio)
    ref = [unconverted_deposition(t, with_pulses([q], pulse, τ); ratio) for τ in times]
    for c in (:omega1, :omega2)
        A = hq.absorption[c] .* hq.element[c] .* t.elem_area ./ 1e14
        R = reduce(hcat, [getfield(d, c).element for d in ref])
        @test isapprox(A, R; rtol = 1e-12) && sum(R) > 0
    end

    # Beams on their own pulses: one deposition per schedule, the sum of the separate
    # histories, and the same as `unconverted_deposition` at each time
    b2 = nif_beams()[2]
    late = PulseShape([1.0, 1.5, 3.0], [0.0, 1.0, 1.0])
    pulses = Dict(Symbol(bl.name) => pulse, Symbol(b2.name) => late)
    hb = unconverted_history(t, [bl, b2], pulses, times; ratio, by = :beam)
    ha = unconverted_history(t, [bl], pulse, times; ratio)
    hl = unconverted_history(t, [b2], late, times; ratio)
    @test isapprox(hb.element[:omega1], ha.element[:omega1] .+ hl.element[:omega1]; rtol = 1e-12)
    refb = [unconverted_deposition(t, with_pulses([bl, b2], pulses, τ; by = :beam); ratio).omega2.element
            for τ in times]
    @test isapprox(hb.element[:omega2] .* t.elem_area ./ 1e14, reduce(hcat, refb); rtol = 1e-12)
    # a fraction per beam (keyed by name): 1ω on the first beam only, on a falling
    # fraction, and 2ω on the second only; against `unconverted_deposition` at each time
    # with those fractions per beam, and the same as a vector per beam
    n1, n2 = Symbol(bl.name), Symbol(b2.name)
    rb = (omega1 = Dict(n1 => f1, n2 => 0.0), omega2 = Dict(n1 => 0.0, n2 => 0.02))
    hk = unconverted_history(t, [bl, b2], pulses, times; ratio = rb, by = :beam)
    refk = [unconverted_deposition(t, with_pulses([bl, b2], pulses, τ; by = :beam);
                                   ratio = (omega1 = [f1(τ), 0.0], omega2 = [0.0, 0.02])) for τ in times]
    for c in (:omega1, :omega2)
        A = hk.absorption[c] .* hk.element[c] .* t.elem_area ./ 1e14
        R = reduce(hcat, [getfield(d, c).element for d in refk])
        @test isapprox(A, R; rtol = 1e-12) && sum(R) > 0
    end
    hkv = unconverted_history(t, [bl, b2], pulses, times; ratio = (omega1 = [f1, 0.0], omega2 = [0.0, 0.02]), by = :beam)
    @test isequal(hk.element, hkv.element) && isequal(hk.backside, hkv.backside)
    @test_throws KeyError unconverted_history(t, [bl, b2], pulses, times; ratio = (omega1 = Dict(n1 => f1), omega2 = 0.0),
                                              by = :beam)

    # Threshold crossing on a linear ramp: P = P0 (τ - τ0)/Δ, so the 1ω irradiance inside
    # the square reaches I* at τ0 + Δ I* 22² / (1e14 r P0), between the sampled times
    P0, τ0, Δ, r = 4.0, 1.0, 2.0, 0.05
    ramp = PulseShape([τ0, τ0 + Δ], [0.0, P0])
    hr = unconverted_history(t, [bl], ramp, 0.0:0.25:3.5; ratio = (omega1 = r, omega2 = 0.0))
    Istar = 3e9
    τstar = τ0 + Δ * Istar * 22.0^2 / (1e14 * r * P0)
    x = threshold_crossing(hr, Istar; colors = (:omega1,))
    @test isapprox(x.time, τstar; rtol = 1e-12) && x.part === :plate && x.element in full1
    @test threshold_crossing(hr, 1e20).time == Inf
    s = only(specular_check(hr; threshold = Istar))
    @test !s.specular && isapprox(s.crossing, τstar; rtol = 1e-12) && s.peak_time == 3.0 &&
          isapprox(s.peak, 1e14 * r * P0 / 22.0^2; rtol = 1e-12)
    @test only(specular_check(hr; threshold = 5e10)).specular

    # A 1ω prepulse ramping to 1 GW per beamlet over the 5 ns before the pulse: it crosses
    # 1e8 W/cm² where 1e14 P/22² = 1e8, P = 4.84e-4 TW, at τ = -5 + 5 · 0.484
    pre = PulseShape([-5.0, 0.0, 0.0], [0.0, 1e-3, 0.0])
    tp = collect(-6.0:0.25:5.0)
    hp = unconverted_history(t, [bl], pulse, tp; ratio = (omega1 = 0.0, omega2 = 0.0),
                             prepulse = (omega1 = pre,))
    onset = -5 + 5 * 4.84e-4 / 1e-3
    @test isapprox(threshold_crossing(hp, 1e8).time, onset; rtol = 1e-12)
    # the beam's own 3ω reaches 1e8 W/cm² on the plate within the first foot sample: on a
    # 0.25 ns grid, at the ramp's first sampled point past zero (τ = 0.25 ns, the peak
    # element's 3ω irradiance at 0.1 TW ≫ 1e8), linearly from τ = 0
    op = DepositionOperator(t, [bl])
    h3 = irradiance_history(op, [bl], pulse, tp)
    # (`deposit` at each time's power, and the whole pulse on the plate: the spot's
    # quadrature conserves its power to ~1e-8)
    R3 = reduce(hcat, [deposit(op, with_pulses([bl], pulse, τ)).element for τ in tp])
    @test isapprox(h3.element[:omega3] .* t.elem_area ./ 1e14, R3; rtol = 1e-12)
    @test isapprox(vec(sum(R3; dims = 1)), pulse.(tp); rtol = 1e-7)
    I3 = maximum(h3.element[:omega3][:, tp .== 0.25])
    @test isapprox(threshold_crossing(h3, 1e8).time, 0.25 * 1e8 / I3; rtol = 1e-12)
    pc = preplasma_check(hp, h3; foil = :plate)
    @test pc.preplasma && isapprox(pc.onset, onset; rtol = 1e-12) &&
          isapprox(pc.lead, pc.arrival - onset; rtol = 1e-12) && 0 < pc.arrival < 0.01
    nope = preplasma_check(hp, -4.0; foil = :plate)
    @test !nope.preplasma && nope.arrival == -4.0
    @test preplasma_check(h, 1.0; foil = :plate, threshold = 1e20).lead == -Inf

    # inputs
    @test_throws ArgumentError unconverted_history(t, [bl], pulse, times; ratio = (omega1 = 0.1,))
    @test_throws ArgumentError unconverted_history(t, [bl], pulse, times; ratio,
                                                   prepulse = (omega3 = pre,))
    # keys are the colors' Symbols: a String key is an error, never a dropped entry
    @test_throws ArgumentError unconverted_history(t, [bl], pulse, times; ratio,
                                                   prepulse = Dict("omega1" => pre))
    @test_throws ArgumentError unconverted_history(t, [bl], pulse, times;
                                                   ratio = Dict(:omega1 => 0.1, :omega2 => 0.0, "omega1" => 0.2))
    @test_throws ArgumentError unconverted_deposition(t, [bl]; ratio = Dict("omega1" => 0.1, "omega2" => 0.0))
    @test_throws ArgumentError unconverted_deposition(t, [bl]; ratio = (omega1 = 0.1, omega2 = 0.0, omega3 = 0.1))
    @test_throws ArgumentError unconverted_history(t, [bl], pulse, [1.0, 0.0]; ratio)
    @test_throws ArgumentError unconverted_history(t, [bl], pulse, times; ratio = (omega1 = -0.1, omega2 = 0.0))
    @test_throws ArgumentError threshold_crossing(h, 1e8; colors = (:omega3,))
end
