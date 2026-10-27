# Implosion metrics: mode 1 → hot-spot velocity against published transfer numbers,
# the piston model's tolerance table, kernel and swing integrals against independent
# quadrature, and the velocity uncertainty, inverse and penalty against full
# evaluations.

@testset "implosion metrics: published numbers" begin
    # a dipole drive 1 + ε n·r̂ on a sphere, through the exact facet Y_lm weights
    ts = Target(Component(:s, icosphere_mesh(1.0; n = 8), gold()))
    els = collect(1:nelements(ts))
    W = ylm_weights(ts, els; lmax = 1)
    n = normalize(Vec3(0.3, -0.5, 0.8))
    dipole(ε) = W * [1 + ε * dot(n, normalize(ts.elem_centroid[e])) for e in els]
    c = dipole(0.0044)
    @test mode1_vector(c) ≈ 0.0044 * n rtol = 2e-3
    v = hotspot_velocity(c)
    @test normalize(v) ≈ -n rtol = 1e-6                # away from the high-flux side
    # Young et al. 2020, N170601: 0.44 % → 37.6 km/s (3D HYDRA), 44.7 km/s (VF + G)
    @test 37.6 <= norm(v) <= 44.7
    @test abs(norm(v) / 44.7 - 1) < 0.02
    # Kritcher et al. 2024, N230729: ~0.7 % → 70–75 km/s burn-off
    @test 70 * (1 - 2e-3) <= norm(hotspot_velocity(dipole(0.007))) <= 75
    # Spears et al. 2014: 2 % at peak → 120 km/s, saturated; the linear default
    # (valid below ~0.75 %, Young) over-predicts it
    @test norm(hotspot_velocity(dipole(0.02))) > 1.5 * 120

    # Hurricane et al. 2020, Table I (v_imp = 380 km/s, near ignition, χ ≈ 1–2): the
    # allowed f to two digits and v_HS (rounded to 10 km/s in the table) at χ = 2
    for (loss, f, vmax) in ((0.05, 0.08, 30), (0.10, 0.11, 40), (0.15, 0.14, 50), (0.20, 0.16, 60))
        tol = mode1_tolerance(loss; vimp = 380, a = 3.3, chi = 2.0)
        @test round(tol.f; digits = 2) == f
        @test abs(tol.velocity - vmax) <= 5
        m = implosion_metrics(Vec3(0.0, 0.0, tol.velocity); vimp = 380, chi = 2.0)
        @test m.total ≈ 1 - loss rtol = 1e-9
    end
    # Casey et al. 2023: a = 6 − 1.5 ln T gives Hurricane's 3.3 at 6 keV
    @test round(reactivity_exponent(6.0); digits = 1) == 3.3
    # Hurricane Eq. 20, 23: closed forms
    m = implosion_metrics([30.0, -40.0, 0.0]; vimp = 400, ke = 20.0)
    @test m.f ≈ 0.125 && m.speed ≈ 50
    @test m.yield_factor ≈ (1 - 0.125^2)^3.3 && m.rke ≈ 0.66 * 0.125^2 * 20 && m.total == m.yield_factor
    @test yield_factor(1.2) == 0
    # Nora et al. 2023 Eq. 4, low branch: its Table IV lists 1/Vhs_amp^burn per shot
    # (N130530, N140520, N210207, N211107, N220919, N210808)
    for (Y, inv) in ((1.13, 0.956), (2.13, 0.885), (4.42, 0.789), (13.2, 0.649), (25.0, 0.580), (27.6, 0.570))
        @test abs(1 / burn_amplification(Y) - inv) <= 1e-3     # Y_amp itself is rounded
    end
    # high branch with the coefficients of the fit drawn in Fig. 6 (Eq. 4 rounds them)
    @test burn_amplification(100.0) ≈ 13.105 * log(1685.225) - 95.539
    @test 1.6 < burn_amplification(80.0) < 1.7 && 4.5 < burn_amplification(500.0) < 4.7
    @test_throws DomainError burn_amplification(60.0)      # between the fits
    @test_throws DomainError burn_amplification(600.0)
    # shape term: Eq. 2 (2D fit) is continuous at 0.35 to the printed digits
    @test shape_amplification(0.0) ≈ 0.71
    @test abs(shape_amplification(0.35 - 1e-9) - shape_amplification(0.35 + 1e-9)) < 5e-3
    # Table I (3D): its nodes, the equatorial fold, and N130530 of Table IV (P2 −2 %,
    # flow at 90°: 1/Vhs_amp^shape = 0.965)
    @test shape_amplification(-0.5656, 0.0) ≈ 0.046
    @test shape_amplification(0.4028, 30.0) ≈ 1.027
    @test shape_amplification(0.1607, 120.0) ≈ shape_amplification(0.1607, 60.0)
    @test round(1 / shape_amplification(-0.02, 90.0); digits = 3) == 0.965
    # generic in the number type
    @test ForwardDiff.derivative(x -> yield_factor(x / 380), 50.0) ≈
          -3.3 * (1 - (50 / 380)^2)^2.3 * 2 * 50 / 380^2
end

@testset "implosion metrics: kernels and swing" begin
    # a mode-1 ramp along z, 0 → 1 % over 10 ns, on steps that miss the kernel's edges
    times = collect(0.0:1.5:10.5)
    ε = 0.001 .* times
    ylm = zeros(4, length(times))
    ylm[1, :] .= 1.0
    ylm[mode_index(1, 0), :] .= ε ./ sqrt(3)
    moments = vcat(ones(1, length(times)), ε')
    hist = (; time = times, ylm, moments)
    # step kernel: −5 km/s/%/ns to 4 ns, 20 after; ∫ G · 100 ε dt = −4 + 84 = 80 km/s
    G = ResponseKernel([0.0, 4.0, 4.0, 10.0], [-5.0, -5.0, 20.0, 20.0])
    @test G(3.9) == -5 && G(4.0) == 20 && G(10.1) == 0
    @test hotspot_velocity(hist; kernel = G) ≈ Vec3(0.0, 0.0, -80.0) atol = 1e-10
    @test hotspot_velocity((; hist..., ylm = nothing); kernel = G) ≈ Vec3(0.0, 0.0, -80.0) atol = 1e-10
    # a function kernel against an independent midpoint rule
    g(t) = sin(t)
    ref = -100 * sum(g(τ) * 0.001τ for τ in range(0, 10.5; length = 200_001)[1:end-1] .+ 10.5 / 400_000) * 10.5 / 200_000
    @test hotspot_velocity(hist; kernel = g)[3] ≈ ref rtol = 1e-4
    # the flat kernel: a constant 1 % over the window gives per_percent
    flat = (; time = times, ylm = vcat(ones(1, length(times)), zeros(1, length(times)),
                                      fill(0.01 / sqrt(3), 1, length(times)), zeros(1, length(times))),
            moments = nothing)
    @test hotspot_velocity(flat; window = (2.0, 7.0))[3] ≈ -100
    @test hotspot_velocity(flat; window = (2.0, 7.0), per_percent = 110, amplification = 1.3)[3] ≈ -143

    # P2 swing against dense sampling of the piecewise-linear history
    tp = [0.0, 2.0, 4.0, 6.0, 8.0]
    p2 = [0.05, -0.02, 0.0, 0.04, -0.01]
    h2 = (; time = tp, moments = vcat(fill(2.0, 1, 5), zeros(1, 5), 2 .* p2'), ylm = nothing)
    sw = swing_metrics(h2; window = (1.0, 7.0))
    dense = range(1.0, 7.0; length = 600_001)
    mid = (dense[1:end-1] .+ dense[2:end]) ./ 2
    yv = [Briareus._interp(tp, p2, τ) for τ in mid]
    μ = sum(yv) / length(yv)
    @test sw.mean ≈ μ rtol = 1e-6
    @test sw.deviation ≈ sum(abs.(yv .- μ)) / length(yv) rtol = 1e-6
    @test sw.max ≈ 0.04 && sw.min ≈ -0.02 && sw.swing ≈ 0.06
    @test sw.picket_mean ≈ 0.015 && sw.picket_max ≈ 0.05 && sw.picket_ok
    @test !swing_metrics(h2; picket = (0.0, 1.0), limit = 0.03).picket_ok
    @test sw.response == (;)
    # Masse et al. 2019's sign structure: an early pole-hot P2 alone gives opposite
    # hot-spot and ρR responses when the hot-spot kernel flips sign before 4.5 ns
    early = (; time = tp, moments = vcat(ones(1, 5), zeros(1, 5), [0.05 0.05 0.0 0.0 0.0]), ylm = nothing)
    r = swing_metrics(early; kernels = (; hotspot = ResponseKernel([0, 4.5, 4.5, 8], [-1, -1, 1, 1]),
                                         rhoR = ResponseKernel([0, 8], [1, 1]))).response
    @test r.hotspot < 0 < r.rhoR
    @test r.rhoR ≈ 100 * (0.05 * 2 + 0.05 * 2 / 2)
end

@testset "implosion metrics: uncertainty, inverse, penalty" begin
    t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 wall = gold(albedo = 0.85, conversion = 0.8), h = 0.6, capsule_h = 0.2)
    pub = Dict(:outer50 => (5.01, 0.0), :outer44 => (5.20, 0.0),
               :inner30 => (3.63, -0.60), :inner23 => (3.09, -1.09))
    p0 = GroupPointing(nif_quads(power = 8.0); fields = (:z, :r))
    x0 = reduce(vcat, [[pub[k]...] for k in p0.keys])
    quads = apply(p0, x0)
    s = Scenario(t, GroupPointing(quads; fields = (:z, :r)); observe = :capsule,
                 viewfactor_kw = (; closed = true))

    # a constant drive history reproduces the static velocity (a dropped inner quad)
    dropped = [b.name == quads[3].name ? remake(b; power = 0.0) : b for b in quads]
    vstatic = hotspot_velocity(evaluate(s, dropped))
    @test norm(vstatic) > 50
    h = drive_history(s.model, s.deposition, τ -> dropped, [0.0, 1.0, 2.5]; ylm = true, lmax = 2)
    @test hotspot_velocity(h) ≈ vstatic rtol = 1e-6
    @test implosion_metrics(h; vimp = 380).speed ≈ norm(vstatic) rtol = 1e-6

    # velocity covariance against a Monte Carlo of full evaluations
    vimp = 60.0                                 # small, so the yield factor varies
    u = velocity_uncertainty(s, quads; pointing = 0.03, power = 0.02, vimp, chi = 1.0)
    @test u.rms_by.beam ≈ u.rms && u.rms_by.group == 0
    if LONG_TESTS
        samples = sample_modes(s, quads; pointing = 0.03, power = 0.02, n = 400, lmax = 1,
                               rng = MersenneTwister(5))
        vs = [hotspot_velocity(samples[:, k]) for k in axes(samples, 2)]
        rms_mc = sqrt(sum(sum(abs2, v - u.velocity) for v in vs) / length(vs))
        @test abs(rms_mc / u.rms - 1) < 0.08
        ey = sum(exp(-3.3 * sum(abs2, v) / vimp^2) for v in vs) / length(vs)
        @test abs(ey / u.expected_yield - 1) < 0.03
        et = sum(exp(-(3.3 + 2.3) * sum(abs2, v) / vimp^2) for v in vs) / length(vs)
        @test abs(et / u.expected_total - 1) < 0.05
    end
    @test u.expected_total < u.expected_yield < 1

    # the inverse: a 20 µm vertical target offset, measured as a full evaluation's velocity
    D = Vec3(0, 0, 0.02)
    shifted = [remake(b; aim = Vec3(b.aim) + D) for b in quads]
    vm = hotspot_velocity(evaluate(s, shifted))
    ie = infer_errors(s, quads, vm; sigma = 0.5, pointing = 0.0, positioning = 0.05, thickness = 0.0)
    @test ie.positioning[3] ≈ 0.02 rtol = 0.05
    @test norm(ie.velocity.residual) < 2
    @test ie.velocity.chi2 ≈ sum(abs2, vm) / 0.25 / 3 rtol = 1e-6       # nominal is round
    # thickness alone (Casey 2021): +x thicker moves the hot spot toward +x, opposite
    # to a +x flux mode 1
    tt = Vec3(0.003, 0.0, 0.0)
    vt = 100 * 110 * tt
    it = infer_errors(s, quads, vt; sigma = 1.0, pointing = 0.0, thickness = 0.01)
    @test it.thickness ≈ tt rtol = 0.01 atol = 1e-6
    @test it.sources[1].source == "capsule thickness"
    @test hotspot_velocity([1.0, 0.0, 0.0, 0.01])[1] < 0 < vt[1]
    @test 0 < it.velocity.chi2_prior < it.velocity.chi2

    # the penalty: |v|²/scale², plus the errors' variance; gradients through duals
    P = VelocityPenalty(; scale = 20.0)
    ev = evaluate(s, dropped)
    @test P(ev) ≈ sum(abs2, vstatic) / 400
    Pu = VelocityPenalty(; scale = 20.0, uncertainty = (; pointing = 0.03, power = 0.02))
    @test Pu(evaluate(s, quads)) - P(evaluate(s, quads)) ≈ u.rms^2 / 400 rtol = 1e-8
    pd = GroupPointing(dropped; fields = (:z, :r))   # mode 1 from the dropped quad
    xg = initial_parameters(pd)
    x1 = xg .+ 0.05 .* sin.(1:length(xg))
    J(x) = P(evaluate(s, apply(pd, x)))
    gr = ForwardDiff.gradient(J, x1)
    k = argmax(abs.(gr))
    e = zeros(length(xg)); e[k] = 1e-5
    @test gr[k] ≈ (J(x1 + e) - J(x1 - e)) / 2e-5 rtol = 1e-5
end

@testset "mode analysis about an offset capsule" begin
    # a pure P2 field about a displaced sphere has no mode 1 about the sphere's center;
    # analyzed about the origin it would leak into one
    c = Vec3(0.3, 0.0, 0.2)
    m = translate(icosphere_mesh(0.8; n = 8), c)
    ts = Target(Component(:capsule, m, ablator(); surface = SphereSurface(0.8; center = c)))
    els = collect(1:nelements(ts))
    @test mode_center(ts, els) == c
    f = [1 + 0.1 * (3 * normalize(ts.elem_centroid[e] - c)[3]^2 - 1) / 2 for e in els]
    a = ylm_weights(ts, els; lmax = 2) * f
    @test mode_amplitude(a, 1) < 1e-4
    @test mode_amplitude(ylm_weights(ts, els; lmax = 2, center = zero(Vec3)) * f, 1) > 5e-3
    @test legendre_moments(f, ts, els; lmax = 2)[3] ≈ 0.1 rtol = 5e-3
    # a hohlraum's capsule carries its center; the wall alone has none
    th = hohlraum(; length = 6.0, radius = 2.0, leh_radius = 1.0, capsule_radius = 0.5,
                  capsule_center = (0.1, 0.0, -0.05), h = 1.0, capsule_h = 0.4)
    @test mode_center(th, elements(th, :capsule)) == Vec3(0.1, 0.0, -0.05)
    @test mode_center(th, elements(th, :wall)) == zero(Vec3)
end

@testset "Young et al. 2020: static view factors and their Green's function" begin
    # N170827's static model (radii at 85 % of Table II, albedos 0.8 / 0.3, conversion
    # 0.9, no windows), the requested per-beam pulses (Fig. 2e) and the mode-1 Green's
    # function (Fig. 5), all digitized; the pointing is the tutorials' scaled to the
    # 11.3 mm length (N170827's is not published). Research notes, "Young 2020 validation".
    K = young2020_kernel()
    @test 112 < sum((K.t[k+1] - K.t[k]) * (K.g[k] + K.g[k+1]) / 2 for k in 1:length(K.t)-1) < 122
    mk(cc = (0.0, 0.0, 0.0)) = hohlraum(; length = 11.3, radius = 0.85 * 3.1, leh_radius = 0.85 * 1.82,
        capsule_radius = 0.85 * 0.98, corner_radius = 0.85 * 1.05, capsule_center = cc,
        wall = gold(albedo = 0.8, conversion = 0.9), h = 0.6, capsule_h = 0.2)
    pub = Dict(:outer50 => (5.01, 0.0), :outer44 => (5.20, 0.0), :inner30 => (3.63, -0.60), :inner23 => (3.09, -1.09))
    p = GroupPointing(nif_quads(power = 8.0); fields = (:z, :r))
    quads = apply(p, reduce(vcat, [[1.1188 * pub[k][1], pub[k][2]] for k in p.keys]))
    P = young2020_pulses()
    pul = (; inner = PulseShape(P.inner.t, 4 .* P.inner.p), outer = PulseShape(P.outer.t, 4 .* P.outer.p))
    side(b) = startswith(String(b.group), "inner") ? :inner : :outer
    times = collect(0.25:0.25:7.75)
    beams_at(bs, drop = 0) = t -> [j == drop ? remake(b; power = 0.0) : b
                                   for (j, b) in enumerate(with_pulses(bs, pul, t; by = side))]
    hist(m, bs, drop = 0; op = DepositionOperator(m.target, bs)) =
        drive_history(m, op, beams_at(bs, drop), times; ylm = true, lmax = 1)
    model = RadiationModel(mk(); closed = true)
    op = DepositionOperator(model.target, quads)        # (aims move below, ports do not)
    Y = young2020_dropped_quads()
    v = Dict{Symbol,Float64}()
    for (g, vy) in zip(Y.cones, Y.velocity)
        g in (:inner23, :outer50) || continue
        h = hist(model, quads, findfirst(q -> q.group === g, quads); op)
        m1 = [100norm(mode1_vector(view(h.ylm, :, k))) for k in axes(h.ylm, 2)]
        y = [Briareus._interp(Y.t, getproperty(Y, g), t) for t in times]
        ok = times .>= 0.5
        @test sqrt(sum(abs2, (m1 - y)[ok]) / count(ok)) < 0.3          # Fig. 10a, % of mode 1
        v[g] = norm(hotspot_velocity(h; kernel = K))
        @test abs(v[g] / vy - 1) < 0.1                                # Fig. 10c
    end
    @test 1.8 < v[:inner23] / v[:outer50] < 2.6                       # Young 317/146 = 2.17
    # whole-target offsets (every aim moves against the target): Young 1.3 and 0.1 km/s/µm
    δ = 0.05
    slope(f) = norm(f(1) - f(-1)) / (2δ * 1000)
    shift(d) = s -> hotspot_velocity(hist(model, [remake(b; aim = b.aim .- s * δ .* d) for b in quads]; op); kernel = K)
    @test abs(slope(shift([0.0, 0.0, 1.0])) - 1.3) < 0.25
    @test slope(shift([1.0, 0.0, 0.0])) < 0.25
    # capsule offset in the hohlraum: Young 1.1 (vertical) and 0.1 (horizontal) km/s/µm;
    # about TCC instead of the capsule's center the horizontal one was ~0.7
    cap(d) = s -> hotspot_velocity(hist(RadiationModel(mk(s * δ .* d); closed = true), quads); kernel = K)
    @test abs(slope(cap((0.0, 0.0, 1.0))) - 1.1) < 0.25
    @test slope(cap((1.0, 0.0, 0.0))) < 0.25
end
