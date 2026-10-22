# Drive uncertainty from beam errors: exact per-beam sensitivities, linear propagation
# against Monte Carlo, and robust objectives.

std_rows(A) = sqrt.(sum(abs2, A .- sum(A; dims = 2) ./ size(A, 2); dims = 2) ./ (size(A, 2) - 1))

@testset "mode uncertainty" begin
    t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 wall = gold(albedo = 0.85, conversion = 0.8), nz = 16, nphi = 24,
                 capsule_ntheta = 10, capsule_nphi = 20)
    pub = Dict(:outer50 => (5.01, 0.0), :outer44 => (5.20, 0.0),
               :inner30 => (3.63, -0.60), :inner23 => (3.09, -1.09))
    p0 = GroupPointing(nif_quads(power = 8.0); fields = (:z, :r))
    x0 = reduce(vcat, [[pub[k]...] for k in p0.keys])
    quads = apply(p0, x0)
    s = Scenario(t, quads; observe = :capsule, viewfactor_kw = (; closed = true),
                 parameterization = GroupPointing(quads; fields = (:z, :r)))
    lmax = 4
    sens = mode_sensitivity(s, quads; lmax)
    nm = (lmax + 1)^2
    @test size(sens.J) == (nm, 3 * length(quads)) && length(sens.names) == 3 * length(quads)
    @test sens.c ≈ flux_ylm(evaluate(s, quads); lmax) ./ flux_ylm(evaluate(s, quads); lmax)[1]
    @test all(iszero, sens.J[1, :])                                # c00/c00 does not move
    # against central differences of full evaluations
    shape(bs) = (c = flux_ylm(evaluate(s, bs); lmax); c ./ c[1])
    for (bi, k) in ((1, 1), (7, 2), (20, 3), (41, 1))
        b = quads[bi]
        e1, e2 = spot_axes(b)
        move(δ) = k == 3 ? remake(b; power = b.power * (1 + δ)) :
                  remake(b; aim = Vec3(b.aim) + δ * (k == 1 ? e1 : e2))
        h = 1e-4
        fd = (shape([j == bi ? move(h) : quads[j] for j in eachindex(quads)]) .-
              shape([j == bi ? move(-h) : quads[j] for j in eachindex(quads)])) ./ 2h
        col = sens.J[:, 3(bi-1)+k]
        @test maximum(abs, fd .- col) <= 1e-6 * max(maximum(abs, col), 1e-6)
    end
    # a power error moves no mode when all beams' powers scale together
    @test maximum(abs, sum(sens.J[:, 3:3:end]; dims = 2)) < 1e-12

    # linear propagation against Monte Carlo (small errors: linear regime)
    u = mode_uncertainty(s, quads; pointing = 0.02, power = 0.02, lmax)
    @test u.sigma ≈ sqrt.(diag(u.covariance)) && length(u.rms_mode) == lmax
    @test u.sigma[1] == 0
    u2 = mode_uncertainty(s, quads; pointing = 0.04, power = 0.04, lmax)
    @test u2.sigma ≈ 2u.sigma rtol = 1e-12                          # linear in the errors
    samples = sample_modes(s, quads; pointing = 0.02, power = 0.02, n = 400, lmax,
                           rng = MersenneTwister(11))
    sd = vec(std_rows(samples))
    big = findall(>(0.3 * maximum(u.sigma)), u.sigma)
    @test length(big) >= 3
    @test all(abs.(sd[big] ./ u.sigma[big] .- 1) .< 0.15)          # MC: ±7% at 2σ for n = 400
    rms_mc = [sqrt(sum(abs2, samples[l^2+1:(l+1)^2, :] .- sens.c[l^2+1:(l+1)^2]) / size(samples, 2))
              for l in 1:lmax]
    @test all(abs.(rms_mc ./ u.rms_mode .- 1) .< 0.12)

    # the robust objectives add exactly the propagated variance
    errs = (; pointing = 0.03, power = 0.01)
    ua = mode_uncertainty(s, quads; errs..., lmax)
    for (modes, comps) in (((1 => 1.0, 2 => 0.5, 4 => 2.0), ()), ((), ((2, 2) => 1.0, (3, -1) => 3.0)))
        nominal = asymmetry_objective(; modes, components = comps)
        robust = asymmetry_objective(; modes, components = comps, uncertainty = errs)
        ev = evaluate(s, quads)
        extra = sum(w * ua.rms_mode[l]^2 for (l, w) in modes; init = 0.0) +
                sum(w * ua.sigma[mode_index(l, m)]^2 for ((l, m), w) in comps; init = 0.0)
        @test robust(ev) - nominal(ev) ≈ extra rtol = 1e-9
    end
    # correlated errors. One beam per group: a shared transverse displacement is the
    # beam's own pointing error
    ug = mode_uncertainty(s, quads; pointing = 0.0, group_pointing = 0.02, groups = :name, lmax)
    up = mode_uncertainty(s, quads; pointing = 0.02, lmax)
    @test ug.covariance ≈ up.covariance rtol = 1e-9
    @test ug.rms_by.beam == zeros(lmax) && ug.rms_by.group ≈ ug.rms_mode
    # one shared power scale for every beam moves nothing
    @test maximum(mode_uncertainty(s, quads; pointing = 0.0, group_power = 0.05, groups = b -> 1, lmax).sigma) < 1e-12
    # target positioning along x alone: σ = |∂c/∂x| against central differences of
    # every aim moved together
    ux = mode_uncertainty(s, quads; pointing = 0.0, positioning = (1.0, 0.0, 0.0), lmax)
    h = 1e-4
    fdx = (shape([remake(b; aim = Vec3(b.aim) + Vec3(h, 0, 0)) for b in quads]) .-
           shape([remake(b; aim = Vec3(b.aim) - Vec3(h, 0, 0)) for b in quads])) ./ 2h
    @test maximum(abs, ux.sigma .- abs.(fdx)) <= 1e-6 * maximum(abs, fdx)
    @test ux.rms_by.positioning ≈ ux.rms_mode
    # all sources together against Monte Carlo (per cone: groups = :group)
    errc = (; pointing = 0.02, power = 0.01, group_pointing = 0.02, group_power = 0.01,
            groups = :group, positioning = 0.01)
    uc = mode_uncertainty(s, quads; errc..., lmax)
    @test uc.rms_mode ≈ sqrt.(uc.rms_by.beam .^ 2 .+ uc.rms_by.group .^ 2 .+ uc.rms_by.positioning .^ 2)
    @test all(uc.rms_by.group .> 0) && all(uc.rms_by.positioning .> 0)
    sc = sample_modes(s, quads; errc..., n = 400, lmax, rng = MersenneTwister(5))
    rms_c = [sqrt(sum(abs2, sc[l^2+1:(l+1)^2, :] .- sens.c[l^2+1:(l+1)^2]) / size(sc, 2)) for l in 1:lmax]
    @test all(abs.(rms_c ./ uc.rms_mode .- 1) .< 0.15)
    @test_throws ArgumentError mode_uncertainty(s, quads; positioning = (0.1, 0.1), lmax)
    # robust objectives with shared errors add exactly the propagated variance
    nominal = asymmetry_objective(; modes = (1 => 1.0, 3 => 2.0))
    robustc = asymmetry_objective(; modes = (1 => 1.0, 3 => 2.0), uncertainty = errc)
    ev = evaluate(s, quads)
    @test robustc(ev) - nominal(ev) ≈ uc.rms_mode[1]^2 + 2uc.rms_mode[3]^2 rtol = 1e-9
    xg = initial_parameters(s.param)
    fc(x) = robustc(evaluate(s, x))
    gc = ForwardDiff.gradient(fc, xg)
    for i in (2, 5)
        h = 1e-5
        e = [j == i ? h : 0.0 for j in eachindex(xg)]
        @test gc[i] ≈ (fc(xg .+ e) - fc(xg .- e)) / 2h rtol = 1e-5
    end

    # and the Legendre one, with the moments' own sensitivities
    sym = symmetry_objective(modes = (2 => 1.0, 4 => 1.0), odd = 0.0)
    symr = symmetry_objective(modes = (2 => 1.0, 4 => 1.0), odd = 0.0, uncertainty = errs)
    ev = evaluate(s, quads)
    @test symr(ev) > sym(ev)
    @test_throws ArgumentError symmetry_objective(uncertainty = errs, center = (0, 0, 0.1))
    @test_throws ArgumentError asymmetry_objective(uncertainty = (; pointing = -1.0))
    # differentiable (nested duals): the gradient against central differences
    x = initial_parameters(s.param)
    f(x) = symr(evaluate(s, x))
    g = ForwardDiff.gradient(f, x)
    for i in (1, 4, 7)
        h = 1e-5
        e = [j == i ? h : 0.0 for j in eachindex(x)]
        @test g[i] ≈ (f(x .+ e) - f(x .- e)) / 2h rtol = 1e-5
    end
end

