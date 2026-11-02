# Per-color x-ray sources in `solve`: a two-plate enclosure solved by hand, the 3ω path
# bit for bit, and the equivalent-3ω weighting for conversions independent of intensity.

@testset "per-color x-ray sources" begin
    uc = nif_unconverted()
    r1, r2 = uc[1].conversion, uc[2].conversion

    # Two directly opposed unit squares 1 mm apart, one element each. View factor between
    # them: the closed form for directly opposed parallel rectangles at X = Y =
    # side/spacing = 1 (Howell's catalog of configuration factors).
    X = 1.0
    fx = 2 / (π * X^2) * (log(sqrt((1 + X^2)^2 / (1 + 2X^2))) + 2X * sqrt(1 + X^2) * atan(X / sqrt(1 + X^2)) -
                          2X * atan(X))
    @test isapprox(fx, 0.19982; atol = 1e-5)
    F = [0.0 fx; fx 0.0]
    # plate 1: Dewald 2008's intensity table; plate 2: a time-dependent conversion
    tab = Conversion(intensities = (1e14, 1e15), efficiencies = (0.95, 0.80))
    m1 = Material(:Au; albedo = 0.6, conversion = tab)
    m2 = Material(:Au; albedo = 0.3, conversion = Conversion(efficiency = 0.5, time_exponent = 0.2))
    sq = rectangle_mesh(1.0, 1.0; nu = 1, nv = 1)
    t2 = Target(Component(:a, sq, m1),
                Component(:b, transform(sq, [1 0 0; 0 -1 0; 0 0 -1], Vec3(0, 0, 1)), m2))
    @test nelements(t2) == 2 && all(≈(1.0), t2.elem_area)
    m = RadiationModel(t2, F)
    P = [3.0, 0.4]                                 # TW on 1 mm²: 3e14 and 4e13 W/cm²
    u = (omega1 = [0.5, 0.1], omega2 = [2.0, 0.0], colors = uc)
    τ = 2.0
    sol = solve(m, P; time = τ, unconverted = u)
    # by hand: η(I) = 0.95 - 0.15 log₁₀(I/1e14), held at 0.95 below 1e14 W/cm²
    η(I) = 0.95 - 0.15 * log10(clamp(I, 1e14, 1e15) / 1e14)
    S1 = η(3e14) * 3.0 + r1 * 0.95 * 0.5 + r2 * η(2e14) * 2.0
    S2 = 0.5 * τ^0.2 * (0.4 + r1 * 0.1)
    α1, α2 = 0.6, 0.3
    B1 = (S1 + α1 * fx * S2) / (1 - α1 * α2 * fx^2)
    B2 = S2 + α2 * fx * B1
    @test isapprox(sol.source, [S1, S2]; rtol = 1e-12)
    @test isapprox(sol.B, [B1, B2]; rtol = 1e-12) && isapprox(sol.H, [fx * B2, fx * B1]; rtol = 1e-12)
    @test sol.laser == P .+ (u.omega1 .+ u.omega2)
    # the equivalent 3ω power, at its own intensity, is another source on plate 1
    Seq = solve(m, P .+ r1 .* u.omega1 .+ r2 .* u.omega2; time = τ).source
    @test abs(Seq[1] - S1) > 1e-3 * S1 && isapprox(Seq[2], S2; rtol = 1e-12)
    # a `Deposition` per color and a `Dict` give the same solution
    dz = zeros(2)
    ud = Dict(:omega1 => Deposition(u.omega1, [0.0, 0.0], [0.6], 0.6), :omega2 => Deposition(u.omega2, dz, [2.0], 2.0),
              :colors => uc)
    @test isequal(solve(m, P; time = τ, unconverted = ud).B, sol.B)
    # derivative in the 3ω power on plate 1: d(ηP)/dP = η + dη/dlnI = η(I) - 0.15/ln 10
    dB = ForwardDiff.derivative(x -> solve(m, [x, 0.4]; time = τ, unconverted = u).B, 3.0)
    dS1 = η(3e14) - 0.15 / log(10)
    @test isapprox(dB, [dS1, α2 * fx * dS1] ./ (1 - α1 * α2 * fx^2); rtol = 1e-12)
    @test_throws ArgumentError solve(m, P; unconverted = (omega1 = u.omega1, colors = uc))
    @test_throws DimensionMismatch solve(m, P; unconverted = (omega1 = [1.0], omega2 = [1.0], colors = uc))
    # the conversions are those of the colors the light comes with: none, an error;
    # an entry for no color, an error; other conversions (0.5) change the source by hand
    @test_throws ArgumentError solve(m, P; unconverted = (omega1 = u.omega1, omega2 = u.omega2))
    @test_throws ArgumentError solve(m, P; unconverted = (u..., omega3 = u.omega1))
    half = [UnconvertedColor(c.name; footprint = c.footprint, edge = c.edge, focal_shift = c.focal_shift,
                             absorption = c.absorption, conversion = 0.5) for c in uc]
    sh = solve(m, P; time = τ, unconverted = (u..., colors = half))
    S1h = η(3e14) * 3.0 + 0.5 * 0.95 * 0.5 + 0.5 * η(2e14) * 2.0
    @test isapprox(sh.source, [S1h, 0.5 * τ^0.2 * (0.4 + 0.5 * 0.1)]; rtol = 1e-12)

    # A hohlraum: the 3ω path bit for bit, and the equivalent 3ω power for conversions
    # independent of intensity (constant gold; nif_gold's time-dependent conversion with a
    # temperature-dependent albedo, iterated)
    th = hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8)
    mh = RadiationModel(th)
    rng = Xoshiro(11)
    Ph = 0.1 .* rand(rng, nelements(th))
    uh = (omega1 = 0.01 .* rand(rng, nelements(th)), omega2 = 0.002 .* rand(rng, nelements(th)), colors = uc)
    zero_u = (omega1 = zeros(nelements(th)), omega2 = zeros(nelements(th)), colors = uc)
    s0, sz = solve(mh, Ph), solve(mh, Ph; unconverted = zero_u)
    @test all(isequal(getfield(s0, f), getfield(sz, f)) for f in (:laser, :source, :B, :H, :albedo))
    Peq = Ph .+ r1 .* uh.omega1 .+ r2 .* uh.omega2
    se, sc = solve(mh, Peq), solve(mh, Ph; unconverted = uh)
    @test isapprox(sc.B, se.B; rtol = 1e-12) && isapprox(sc.H, se.H; rtol = 1e-12)
    @test sc.laser == Ph .+ (uh.omega1 .+ uh.omega2)
    mw = RadiationModel(hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8, wall = nif_gold()); F = mh.F)
    sew, scw = solve(mw, Peq; time = 3.0), solve(mw, Ph; time = 3.0, unconverted = uh)
    @test isapprox(scw.H, sew.H; rtol = 1e-12) && isapprox(scw.albedo, sew.albedo; rtol = 1e-12)
    # linear at fixed albedo and constant conversion: d H / d(scale of P) = H(P), with
    # the unconverted source fixed (the dual path)
    obs = th.comp_elems[findfirst(c -> c.name === :capsule, th.components)]
    dH = ForwardDiff.derivative(x -> sum(solve(mh, x .* Ph; unconverted = uh).H[obs]), 1.0)
    @test isapprox(dH, sum(s0.H[obs]); rtol = 1e-12)
end

# Two directly opposed unit squares one side length apart: Howell's closed form at X = 1.
const _UC_FX = let X = 1.0
    2 / (π * X^2) * (log(sqrt((1 + X^2)^2 / (1 + 2X^2))) + 2X * sqrt(1 + X^2) * atan(X / sqrt(1 + X^2)) - 2X * atan(X))
end

@testset "per-color flux response and converted share" begin
    uc = nif_unconverted()
    r1, r2 = uc[1].conversion, uc[2].conversion
    fx = _UC_FX
    F = [0.0 fx; fx 0.0]
    tab = Conversion(intensities = (1e14, 1e15), efficiencies = (0.95, 0.80))
    m1 = Material(:Au; albedo = 0.6, conversion = tab)
    m2 = Material(:Au; albedo = 0.3, conversion = Conversion(efficiency = 0.5, time_exponent = 0.2))
    sq = rectangle_mesh(1.0, 1.0; nu = 1, nv = 1)
    t2 = Target(Component(:a, sq, m1),
                Component(:b, transform(sq, [1 0 0; 0 -1 0; 0 0 -1], Vec3(0, 0, 1)), m2))
    m = RadiationModel(t2, F)
    P = [3.0, 0.4]
    u = (omega1 = [0.5, 0.1], omega2 = [2.0, 0.0], colors = uc)
    τ = 2.0
    sol = solve(m, P; time = τ, unconverted = u)
    R = flux_response(m, [1, 2]; time = τ, laser = P, unconverted = u)
    @test R isa Dict{Symbol,Matrix{Float64}} && sort(collect(keys(R))) == [:omega1, :omega2, :omega3]
    @test isequal(R[:omega3], flux_response(m, [1, 2]; time = τ, laser = P))
    Hr(u) = R[:omega3] * P .+ R[:omega1] * u.omega1 .+ R[:omega2] * u.omega2
    @test isapprox(Hr(u), sol.H; rtol = 1e-12)
    # away from the secant: the 2ω flux on plate 1 rises from 2e14 to 2.2e14 W/cm², on
    # the table; the 1ω flux (5.5e13) stays below it and plate 2's conversion depends on
    # time only. The response's error is the source change r₂ [η(2e14) − η(2.2e14)] 2.2
    # on plate 1 through the 2×2 radiosity solve.
    u2 = (omega1 = 1.1 .* u.omega1, omega2 = 1.1 .* u.omega2, colors = uc)
    η(I) = 0.95 - 0.15 * log10(clamp(I, 1e14, 1e15) / 1e14)
    ΔS = r2 * (η(2e14) - η(2.2e14)) * 2.2
    ΔB1 = ΔS / (1 - 0.6 * 0.3 * fx^2)
    ΔH = [fx * 0.3 * fx * ΔB1, fx * ΔB1]
    @test isapprox(Hr(u2) .- solve(m, P; time = τ, unconverted = u2).H, ΔH; rtol = 1e-9)
    # the converted share is the per-color source (no re-emission here), where the 3ω
    # conversion of the total power differs
    @test isequal(sol.converted, sol.source)
    @test spectral_solve(sol).converted == sol.source
    @test abs(Briareus._xray_source(t2, sol.laser, τ)[1] - sol.source[1]) > 1e-3 * sol.source[1]

    # a conversion independent of intensity: the response is the per-color solve at any
    # powers
    th = hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8)
    mh = RadiationModel(th)
    rng = Xoshiro(5)
    Ph = 0.1 .* rand(rng, nelements(th))
    uh = (omega1 = 0.01 .* rand(rng, nelements(th)), omega2 = 0.002 .* rand(rng, nelements(th)), colors = uc)
    obs = elements(th, :capsule)
    Rh = flux_response(mh, obs; unconverted = uh)
    Hh(s) = Rh[:omega3] * Ph .+ s .* (Rh[:omega1] * uh.omega1 .+ Rh[:omega2] * uh.omega2)
    Hs(s) = solve(mh, Ph; unconverted = (omega1 = s .* uh.omega1, omega2 = s .* uh.omega2, colors = uc)).H[obs]
    @test isapprox(Hh(1.0), Hs(1.0); rtol = 1e-12) && isapprox(Hh(1.7), Hs(1.7); rtol = 1e-12)
    @test_throws ArgumentError flux_response(mh, obs; unconverted = (omega1 = uh.omega1, colors = uc))
    @test_throws ArgumentError flux_response(mh, obs; unconverted = (omega1 = uh.omega1, omega2 = uh.omega2))
end

@testset "unconverted light in drive histories" begin
    uc = nif_unconverted()
    r1, r2 = uc[1].conversion, uc[2].conversion
    a1, a2 = uc[1].absorption, uc[2].absorption
    fx = _UC_FX

    # Two one-element plates in one beamlet's 3ω focal plane, facing the lens, under its
    # 1ω and 2ω squares (1ω: 4.8–26.8 mm along ô from the 3ω axis; 2ω: 2.9–14.9 mm):
    # plate a at 10 mm in both, plate b at 20 mm in the 1ω square only. Inside a square
    # a plate in the focal plane absorbs absorption × fraction × P × area / footprint².
    # The view factors are given (the opposed-squares closed form), so the 2×2
    # radiosity is solved by hand at each of two steps.
    bl = first(nif_beams())
    k = direction(bl)
    ô = Briareus._own_side(bl.port, Briareus._nif_quad_port(bl.quad))
    ŵ = cross(k, ô)
    fb = Vec3(bl.aim...) - bl.defocus * bl.port
    sq = rectangle_mesh(1.0, 1.0; nu = 1, nv = 1)
    tab = Conversion(intensities = (1e9, 1e11), efficiencies = (0.9, 0.6))
    ηa(I) = 0.9 - 0.15 * log10(clamp(I, 1e9, 1e11) / 1e9)
    ηb(τ) = 0.5 * τ^0.2
    αa, αb = 0.6, 0.3
    Rot = hcat(ŵ, ô, -k)
    tp = Target(Component(:a, transform(sq, Rot, fb + 10ô), Material(:Au; albedo = αa, conversion = tab)),
                Component(:b, transform(sq, Rot, fb + 20ô),
                          Material(:Au; albedo = αb, conversion = Conversion(efficiency = 0.5, time_exponent = 0.2))))
    @test all(≈(1.0), tp.elem_area)
    mp = RadiationModel(tp, [0.0 fx; fx 0.0])
    op = DepositionOperator(tp, [bl])
    pulse = PulseShape([0.0, 3.0], [1.0, 3.0])
    f1 = PulseShape([0.0, 3.0], [0.06, 0.03])
    beams_at(τ) = with_pulses([bl], pulse, τ)
    hp = drive_history(mp, op, beams_at, [1.0, 2.0]; observe = :b, wall = (:a,), lmax = 0, mband = nothing,
                       record = sol -> sol, unconverted = (omega1 = f1, omega2 = 0.01))
    @test hp.time == [1.0, 2.0] && isempty(hp.unconverted.skipped)
    dev = map(enumerate(hp.time)) do (j, τ)
        P = pulse(τ)
        P3 = deposit(op, beams_at(τ)).element
        P1, P2 = a1 * f1(τ) * P / 22.0^2, a2 * 0.01 * P / 12.0^2
        Sa = ηa(1e14 * P3[1]) * P3[1] + r1 * ηa(1e14 * P1) * P1 + r2 * ηa(1e14 * P2) * P2
        Sb = ηb(τ) * (P3[2] + r1 * P1)
        Ba = (Sa + αa * fx * Sb) / (1 - αa * αb * fx^2)
        Bb = Sb + αb * fx * Ba
        s = hp.record[j]
        rel(x, y) = maximum(abs.(x .- y) ./ abs.(y))
        # the hand values sit on the intensity table (1e9–1e11 W/cm²)
        (1e9 < 1e14 * P1 < 1e11 && 1e9 < 1e14 * P2 < 1e11) || return Inf
        max(rel(s.source, [Sa, Sb]), rel(s.B, [Ba, Bb]), rel(s.H, [fx * Bb, fx * Ba]), rel(hp.flux[j], fx * Ba),
            rel(s.laser, P3 .+ [P1 + P2, P1]), rel(hp.laser[:a][j], P3[1] + P1 + P2),
            rel(hp.unconverted.delivered[:omega1][j], f1(τ) * P), rel(hp.unconverted.delivered[:omega2][j], 0.01 * P))
    end
    @test maximum(dev) < 1e-12

    # A 1ω prepulse alone (zero fractions) on the same plates, before and during a 3ω pulse
    # that starts at 1 ns: the 1ω light is absorption × prepulse × area / footprint² on
    # both plates, whatever the 3ω power, by hand as above
    main = PulseShape([1.0, 3.0], [1.0, 3.0])
    pre = PulseShape([0.0, 3.0], [0.05, 0.2])
    hq = drive_history(mp, op, τ -> with_pulses([bl], main, τ), [0.5, 2.0]; observe = :b, wall = (:a,), lmax = 0,
                       mband = nothing, record = sol -> sol, unconverted = (omega1 = 0.0, omega2 = 0.0),
                       prepulse = (omega1 = pre,))
    devq = map(enumerate(hq.time)) do (j, τ)
        P3 = deposit(op, with_pulses([bl], main, τ)).element
        P1 = a1 * pre(τ) / 22.0^2
        Sa = ηa(1e14 * P3[1]) * P3[1] + r1 * ηa(1e14 * P1) * P1
        Sb = ηb(τ) * (P3[2] + r1 * P1)
        Ba = (Sa + αa * fx * Sb) / (1 - αa * αb * fx^2)
        Bb = Sb + αb * fx * Ba
        s = hq.record[j]
        rel(x, y) = maximum(abs.(x .- y) ./ abs.(y))
        1e9 < 1e14 * P1 < 1e11 || return Inf
        max(rel(s.source, [Sa, Sb]), rel(s.B, [Ba, Bb]), rel(s.laser, P3 .+ P1),
            rel(hq.unconverted.delivered[:omega1][j], pre(τ))) + abs(hq.unconverted.delivered[:omega2][j])
    end
    @test hq.delivered[1] == hq.unconverted.delivered[:omega1][1] && maximum(devq) < 1e-12

    # A hohlraum with an intensity-dependent conversion and constant albedos (each step a
    # direct solve): the history against `solve` with `unconverted_deposition` of that
    # step's beams, for a constant fraction and for a time-varying 1ω fraction on per-cone
    # pulses (beams on two schedules)
    tab2 = Conversion(intensities = (1e11, 1e15), efficiencies = (0.9, 0.5))
    th = hohlraum(nz = 8, nphi = 12, capsule_ntheta = 6, capsule_nphi = 8,
                  wall = Material(:Au; albedo = 0.8, conversion = tab2))
    mh = RadiationModel(th)
    quads = nif_quads(power = 8.0)
    oph = DepositionOperator(th, quads)
    late = PulseShape([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 4.0, 4.0])
    pulses = Dict(:inner23 => pulse, :inner30 => pulse, :outer44 => late, :outer50 => late)
    fω1 = PulseShape([0.0, 4.0], [0.08, 0.02])
    times = [0.5, 1.0, 2.0, 3.0]
    cases = ((τ -> with_pulses(quads, pulse, τ), (omega1 = 0.05, omega2 = 0.01), τ -> (omega1 = 0.05, omega2 = 0.01)),
             (τ -> with_pulses(quads, pulses, τ), (omega1 = fω1, omega2 = 0.01), τ -> (omega1 = fω1(τ), omega2 = 0.01)))
    for (bat, ratio, ratio_at) in cases
        h = drive_history(mh, oph, bat, times; record = sol -> sol, unconverted = ratio)
        dev = map(enumerate(times)) do (j, τ)
            b = bat(τ)
            d = deposit(oph, b)
            u = unconverted_deposition(th, b; ratio = ratio_at(τ))
            ref = solve(mh, d; time = τ, unconverted = u)
            s = h.record[j]
            lz = d.delivered + u.omega1.delivered + u.omega2.delivered
            tot = sum(d.element) + sum(d.backside) + sum(u.omega1.element) + sum(u.omega1.backside) +
                  sum(u.omega2.element) + sum(u.omega2.backside)
            maximum(f -> maximum(abs.(getfield(s, f) .- getfield(ref, f))) / maximum(abs.(getfield(ref, f))),
                    (:source, :B, :H, :laser, :converted)) +
            abs(h.delivered[j] - lz) / lz + abs(sum(h.laser[n][j] for n in keys(h.laser)) - tot) / tot
        end
        @test maximum(dev) < 1e-12
        @test h.unconverted.delivered[:omega1][end] > 0 && isempty(h.unconverted.skipped)
    end

    # a prepulse per cone on per-cone pulses (1ω on the outer cones, 2ω on all): each step
    # against `solve` with each color's light from `unconverted_deposition` of the beams at
    # that color's power, r_c P₃ω + prepulse_c. NIF's 1ω and 2ω squares miss this target's
    # inside (the wedged lenses keep them ≥ 2.9 mm off the 3ω axis), so mock squares near
    # the axis put the light on the walls and capsule
    mock = [UnconvertedColor(:omega1; footprint = 2.0, edge = 0.0, focal_shift = 100.0, absorption = 0.88,
                             conversion = 0.88),
            UnconvertedColor(:omega2; footprint = 1.5, edge = 0.3, focal_shift = 60.0, absorption = 1.0, conversion = 0.9)]
    # the deposition carries its colors' conversions into `solve`: the mock squares with
    # conversion 0.5 (NIF's are 0.88 and 0.9), against each color's conversion by hand
    m05 = [UnconvertedColor(c.name; footprint = c.footprint, edge = c.edge, focal_shift = c.focal_shift,
                            absorption = c.absorption, conversion = 0.5) for c in mock]
    u05 = unconverted_deposition(th, quads; ratio = (omega1 = 0.05, omega2 = 0.01), colors = m05)
    P05 = deposit(oph, quads).element
    ηs(p) = conversions(th, 1.0, p) .* p ./ th.elem_area
    S05 = ηs(P05) .+ 0.5 .* (ηs(u05.omega1.element) .+ ηs(u05.omega2.element))
    s05 = solve(mh, P05; time = 1.0, unconverted = u05)
    @test u05.colors == m05 && sum(u05.omega1.element) > 1e-3 * sum(P05)
    @test isapprox(s05.source, S05; rtol = 1e-12)
    # the M-band is the 3ω light's on both paths: a direct `spectral_solve` of a step's
    # solution (with unconverted light) and the history's own, per component; the
    # M-band of all the absorbed light differs
    em = MBandEmission()
    hm = drive_history(mh, oph, τ -> with_pulses(quads, pulse, τ), times[2:3]; record = sol -> sol,
                       unconverted = (omega1 = 0.05, omega2 = 0.01), colors = m05, mband = MBandSplit(em))
    mdev = map(enumerate(hm.time)) do (j, τ)
        s = hm.record[j]
        sm = spectral_solve(s; edges = [0.0, em.band[1], Inf], mband = em)
        mc(x) = Dict(c.name => sum(x[e] * th.elem_area[e] for e in th.comp_elems[ci]; init = 0.0)
                     for (ci, c) in enumerate(th.components))
        d, w = mc(sm.mband), mc(min.(mband_power(em, th, s.laser) ./ th.elem_area, s.source))
        (maximum(abs(d[n] - hm.mband[n][j]) for n in keys(d)) / hm.mband[:wall][j],
         abs(w[:wall] - hm.mband[:wall][j]) / hm.mband[:wall][j])
    end
    @test maximum(first, mdev) < 1e-12 && minimum(last, mdev) > 1e-4
    pre1 = Dict(:outer44 => PulseShape([0.0, 4.0], [0.02, 0.01]), :outer50 => PulseShape([0.0, 4.0], [0.03, 0.0]),
                :inner23 => τ -> 0.0, :inner30 => τ -> 0.0)
    pre2 = PulseShape([0.0, 4.0], [0.004, 0.004])
    bat = τ -> with_pulses(quads, pulses, τ)
    hpp = drive_history(mh, oph, bat, times; record = sol -> sol, unconverted = (omega1 = fω1, omega2 = 0.01),
                        prepulse = (omega1 = pre1, omega2 = pre2), colors = mock)
    inside = Float64[]
    devp = map(enumerate(times)) do (j, τ)
        b = bat(τ)
        d = deposit(oph, b)
        pc = (omega1 = [fω1(τ) * x.power + pre1[x.group](τ) for x in b], omega2 = [0.01 * x.power + pre2(τ) for x in b])
        u = Dict{Symbol,Any}(n => unconverted_deposition(th, [remake(x; power = p) for (x, p) in zip(b, pc[n])];
                                             ratio = (omega1 = 1.0, omega2 = 1.0), colors = mock)[n]
                 for n in (:omega1, :omega2))
        u[:colors] = mock
        push!(inside, sum(u[:omega1].element) / sum(d.element))
        ref = solve(mh, d; time = τ, unconverted = u)
        s = hpp.record[j]
        lz = d.delivered + u[:omega1].delivered + u[:omega2].delivered
        maximum(f -> maximum(abs.(getfield(s, f) .- getfield(ref, f))) / maximum(abs.(getfield(ref, f))),
                (:source, :B, :H, :laser, :converted)) + abs(hpp.delivered[j] - lz) / lz +
        abs(hpp.unconverted.delivered[:omega1][j] - u[:omega1].delivered) / u[:omega1].delivered
    end
    @test maximum(devp) < 1e-12 && minimum(inside) > 1e-3
    # a zero prepulse: the history without one, bit for bit
    hz0 = drive_history(mh, oph, bat, times; record = sol -> sol.B, unconverted = (omega1 = fω1, omega2 = 0.01),
                        colors = mock)
    hzp = drive_history(mh, oph, bat, times; record = sol -> sol.B, unconverted = (omega1 = fω1, omega2 = 0.01),
                        colors = mock, prepulse = (omega1 = τ -> 0.0,))
    @test all(isequal(hz0[f], hzp[f]) for f in keys(hz0))
    @test_throws ArgumentError drive_history(mh, oph, bat, times; prepulse = (omega1 = pre2,))
    @test_throws ArgumentError drive_history(mh, oph, bat, times; unconverted = (omega1 = 0.0, omega2 = 0.0),
                                             prepulse = (omega3 = pre2,))
    # keys are the colors' Symbols: a String key is an error, never a dropped prepulse
    @test_throws ArgumentError drive_history(mh, oph, bat, times; unconverted = (omega1 = 0.0, omega2 = 0.0),
                                             prepulse = Dict("omega1" => pre2))
    @test_throws ArgumentError drive_history(mh, oph, bat, times; unconverted = Dict("omega1" => 0.0, "omega2" => 0.0))

    # zero fractions: the 3ω history, bit for bit
    h0 = drive_history(mh, oph, τ -> with_pulses(quads, pulse, τ), times; record = sol -> sol.B)
    hz = drive_history(mh, oph, τ -> with_pulses(quads, pulse, τ), times; record = sol -> sol.B,
                       unconverted = (omega1 = 0.0, omega2 = 0.0))
    @test h0.unconverted === nothing && all(isequal(h0[f], hz[f]) for f in keys(h0) if f !== :unconverted)
    @test_throws ArgumentError drive_history(mh, oph, τ -> with_pulses(quads, pulse, τ), times;
                                             unconverted = (omega1 = 0.05,))
    @test_throws ArgumentError drive_history(mh, oph, τ -> with_pulses(quads, pulse, τ), times;
                                             unconverted_kw = (; coarse = 0.2))
    moved(τ) = [remake(b; power = pulse(τ), aim = b.aim + Vec3(0.0, 0.0, 0.1τ)) for b in quads]
    @test_throws ArgumentError drive_history(mh, oph, moved, times; unconverted = (omega1 = 0.05, omega2 = 0.0))
end
