# Corona absorption operator (polar direct drive, Stage 1; ROADMAP "Polar direct drive",
# reference/litreview/polar_direct_drive.md decision 3 (b); research notes "Corona
# absorption operator"). References, all independent of the operator's kernels, lattice
# and rotations:
# - closed forms: the radial ray (b = 0) of the truncated exponential corona absorbs
#   1 - exp(-2(ν_c L/c)(2√w - (2/3)w^{3/2})), w = 1 - e^{-(R - r_c)/L} (Kruer's planar
#   cos³ law at normal incidence, (8/3)ν_c L/c, as R → ∞); in n_e/n_c = a/r² a ray turns
#   at r_t = √(a + b²), sweeps 2(b/r_t) arccos(r_t/R) and has optical depth
#   2(ν_c/c)a²[√(R² - r_t²)/(2r_t²R²) + arcsec(R/r_t)/(2r_t³)] inside radius R;
# - the thin-corona limit L/r_c → 0: a hard sphere at r_c absorbing η_a(ψ) I cos ψ
#   (Murakami & Nishi 2017, Eq. 29), by 1-D quadrature about the beam axis;
# - 3-D ray tracing by RK4 through the same corona (test/corona_raytrace.jl, the
#   prototype's model (c));
# - ForwardDiff against central differences.

include("corona_raytrace.jl")

# ∫₀¹ by Gauss–Legendre, and Legendre polynomials (local to this file).
function _co_gauss01(n)
    x, w = Briareus._gauss_legendre(n)
    return (x .+ 1) ./ 2, w ./ 2
end
function _co_P(l, μ)
    l == 0 && return 1.0
    p0, p1 = 1.0, μ
    for n in 1:l-1
        p0, p1 = p1, ((2n + 1) * μ * p1 - n * p0) / (n + 1)
    end
    return p1
end

@testset "corona absorption" begin
    rc, L = 0.40, 0.08
    η⊥ = 1 - exp(-2.4)                     # (8/3) ν_c L/c = 2.4
    νc = 2.4 / (8 / 3 * L)
    prof = CoronaProfile(; critical_radius = rc, scale_length = L, absorption = η⊥)
    R = rc + 14L

    @testset "kernels: closed forms and energy" begin
        # radial ray of the truncated exponential corona
        w = 1 - exp(-(R - rc) / L)
        α, f, esc, θ = Briareus._ray_deposits(prof, 0.0, 200)
        @test isapprox(1 - esc, 1 - exp(-2νc * L * (2sqrt(w) - 2 / 3 * w^1.5)); rtol = 1e-12)
        @test θ == 0 && all(==(0), α)
        # n_e/n_c = a/r² (functions, not a table): swept angle and depth at every b
        a, ν, Ra = 0.04, 2.0, 2.0
        ne(r) = r > Ra ? 0.0 : a / r^2
        rate(r) = ν * ne(r)^2
        errs = map((0.0, 0.1, 0.3, 0.6, 1.0, 1.5)) do b
            rt = sqrt(a + b^2)
            α, f, esc, θ = Briareus._ray_deposits(ne, rate, Ra, 0.1, b, 200)
            τ = 2ν * a^2 * (sqrt(Ra^2 - rt^2) / (2rt^2 * Ra^2) + acos(rt / Ra) / (2rt^3))
            (abs(θ - asin(b / Ra) - 2(b / rt) * acos(rt / Ra)), abs(-log(esc) / τ - 1),
             abs(sum(f) + esc - 1))
        end
        @test maximum(e -> e[1], errs) < 1e-10      # 2.5e-12 rad measured
        @test maximum(e -> e[2], errs) < 1e-10      # 1.8e-12
        @test maximum(e -> e[3], errs) < 1e-14      # absorbed + escaped = 1
        # a tabulated isothermal snapshot is the analytic corona: ln n_e is linear in r,
        # which the Hermite table reproduces on any grid
        Te, Z, lnΛ = 3600.0, 3.5, 8.0
        nc = 1.115e21 / 0.351^2                       # NRL: n_c = 1.115e21/λ² cm⁻³
        νt = 2.91e-6 * Z * nc * lnΛ / Te^1.5 / 2.998e11   # NRL ν_e × Z, per mm of light path
        r = sort!(vcat(rc - L, R, rc - L .+ (R - rc + L) .* rand(MersenneTwister(3), 30)))
        tab = CoronaProfile(r, exp.(-(r .- rc) ./ L), fill(Te, length(r)); Z, lnΛ)
        ana = CoronaProfile(; critical_radius = rc, scale_length = L,
                            absorption = 1 - exp(-8νt * L / 3))
        @test isapprox(Briareus._corona_rate(tab, 0.5), νt * exp(-2(0.5 - rc) / L); rtol = 2e-3)
        fa = [1 - Briareus._ray_deposits(ana, b, 100)[3] for b in 0:0.1:1.4]
        ft = [1 - Briareus._ray_deposits(tab, b, 100)[3] for b in 0:0.1:1.4]
        @test isapprox(ft, fa; rtol = 2e-3)                 # the constants' rounding
        @test_throws ArgumentError CoronaProfile(r, fill(0.5, length(r)), fill(Te, length(r)); Z, lnΛ)
    end

    # three OMEGA beams, one per PDD ring, shifted normal to the beam toward the equator
    # (Mn14's ring shifts; research notes "Corona model prototype")
    spot = SuperGaussianSpot(0.308; order = 3.66)
    all40 = omega_beams(; cones = (:cone1, :cone2, :cone3), spot, power = 1.0, fnumber = 6.7)
    base = [all40[findfirst(b -> b.group == g, all40)] for g in (:cone1, :cone2, :cone3)]
    U = [sign(b.port[3]) * Briareus.port_basis(b.port)[1] for b in base]
    pdd(Δ) = [remake(b; aim = Δ[i] * U[i]) for (i, b) in enumerate(base)]
    Δ0 = [0.09, 0.15, 0.15]
    beams = pdd(Δ0)
    op = CoronaOperator(prof, beams; ncone = 3)
    dep = deposit(op, beams)

    @testset "against 3-D ray tracing" begin
        c, absorbed, _ = _rt_corona(beams, (rc, L, νc, R); ncone = 3, hlat = 0.04, hτ = 0.008,
                                    nθ = 720)
        # measured: absorbed 5e-7 relative, max |Δc_lm|/c₀₀ 2.9e-5 (l ≤ 30; the ray
        # trace's angular grid; 1.3e-5 at half its lattice spacing)
        @test isapprox(sum(dep.per_beam), absorbed; rtol = 5e-6)
        @test maximum(abs, c ./ c[1] .- dep.ylm ./ dep.ylm[1]) < 1e-4
        # energy: the harmonics carry the absorbed power, which is below the delivered
        @test isapprox(sqrt(4π) * dep.ylm[1], sum(dep.per_beam); rtol = 1e-12)
        @test all(0 .< dep.per_beam .< 1) && dep.delivered == 3
    end

    @testset "thin-corona limit: the η_a weight" begin
        # one collimated beam aimed at the center; reference: a hard sphere at r_c
        # absorbing η_a(μ) I cos ψ, Legendre spectrum about the beam axis
        sp = SuperGaussianSpot(0.25; order = 4.0)
        b1 = [Beam("t"; theta = 37.0, phi = 20.0, spot = sp, power = 1.0, fnumber = 6.7)]
        η = 0.6
        μ, wμ = _co_gauss01(400)
        F = [rc^2 * m * focal_profile(sp, 1.0, rc * sqrt(1 - m^2) / 0.25) * (1 - (1 - η)^(m^3)) for m in μ]
        I0 = sum(wμ .* F)
        σref = [sqrt(2l + 1) * abs(sum(wμ .* F .* _co_P.(l, μ))) / I0 for l in 0:30]
        err = map((0.01, 0.005)) do Lr
            pr = CoronaProfile(; critical_radius = rc, scale_length = Lr * rc, absorption = η)
            r = deposit(CoronaOperator(pr, b1; ncone = 1), b1)
            (abs(sum(r.per_beam) - 2π * I0), maximum(abs.(mode_power(ylm_matrix(r.ylm, 30)) .- σref)))
        end
        # first order in L/r_c (measured 2.5e-4, 1.3e-4 absorbed; 2.7e-2, 1.4e-2 in σ_l)
        @test err[2][1] < 2e-4 && err[2][2] < 2e-2
        @test 1.8 < err[1][1] / err[2][1] < 2.2 && 1.8 < err[1][2] / err[2][2] < 2.2
    end

    @testset "derivatives and checks" begin
        fσ(x) = mode_rms(deposit(op, pdd(x)).ylm, 30)
        g = ForwardDiff.gradient(fσ, Δ0)
        h = 1e-6
        fd = [(fσ(Δ0 .+ h .* (1:3 .== i)) - fσ(Δ0 .- h .* (1:3 .== i))) / 2h for i in 1:3]
        @test isapprox(g, fd; rtol = 1e-6)                  # 1.1e-8 measured
        # dual beams give the Float64 values (to round-off: CI CPUs fuse the dual
        # arithmetic differently, ~1e-16 of the norm)
        dd = deposit(op, pdd(ForwardDiff.Dual.(Δ0, 1.0)))
        @test isapprox(ForwardDiff.value.(dd.ylm), dep.ylm; rtol = 1e-12)
        @test_throws ArgumentError deposit(op, [remake(beams[1]; port = beams[2].port), beams[2:3]...])
        @test_throws DimensionMismatch deposit(op, beams[1:2])
    end

    @testset "deposition per radius" begin
        # radial ray of the truncated exponential corona: per shell, from the one-pass
        # depth outside radius r, τ(r) = (ν_c L/c)(H(x(r)) - H(x(R))), x = n_e/n_c,
        # H(x) = -2√(1-x) + (2/3)(1-x)^{3/2} (∫ x/√(1-x) dx, with dr = -L dx/x)
        H(x) = -2sqrt(1 - x) + 2 / 3 * (1 - x)^1.5
        τin(r) = νc * L * (H(exp(-(clamp(r, rc, R) - rc) / L)) - H(exp(-(R - rc) / L)))
        T = τin(rc)
        edges = [rc - L; rc .+ L .* [0.25, 0.5, 1, 2, 4, 8]; R]
        ref = [exp(-τin(b)) - exp(-τin(a)) + exp(-(2T - τin(a))) - exp(-(2T - τin(b)))
               for (a, b) in zip(edges[1:end-1], edges[2:end])]
        _, f, _, _, sh = Briareus._ray_deposits(prof, 0.0, 200, edges)
        @test maximum(abs, [sum(f[sh .== s]) for s in eachindex(ref)] .- ref) < 1e-12   # 5.7e-14 measured
        # oblique rays in n_e/n_c = a/r²: the one-pass depth inside radius ρ,
        # ν a² [√(ρ² - r_t²)/(2r_t²ρ²) + arcsec(ρ/r_t)/(2r_t³)], r_t = √(a + b²)
        a, ν, Ra = 0.04, 2.0, 2.0
        ne(r) = r > Ra ? 0.0 : a / r^2
        rate(r) = ν * ne(r)^2
        e2 = collect(range(0.1, Ra; length = 12))
        errs = map((0.0, 0.3, 0.6, 1.0, 1.5)) do b
            rt = sqrt(a + b^2)
            τ1(r) = (ρ = clamp(r, rt, Ra); ν * a^2 * (sqrt(ρ^2 - rt^2) / (2rt^2 * ρ^2) + acos(rt / ρ) / (2rt^3)))
            T1 = τ1(Ra)
            r1 = [exp(-(T1 - τ1(hi))) - exp(-(T1 - τ1(lo))) + exp(-(T1 + τ1(lo))) - exp(-(T1 + τ1(hi)))
                  for (lo, hi) in zip(e2[1:end-1], e2[2:end])]
            _, f, _, _, sh = Briareus._ray_deposits(ne, rate, Ra, 0.1, b, 200, e2)
            maximum(abs, [sum(f[sh .== s]) for s in eachindex(r1)] .- r1)
        end
        @test maximum(errs) < 1e-12                          # 3.9e-14 measured
        # the operator: shells partition the totals
        ops = CoronaOperator(prof, beams; ncone = 3, shells = [rc - L; rc .+ L .* (-0.5:0.5:4); R])
        ds = deposit(ops, beams)
        @test isapprox(sum(ds.radial), sum(ds.per_beam); rtol = 1e-12)
        @test maximum(abs, sum(eachcol(ds.shell_ylm)) .- ds.ylm) < 1e-12 * ds.ylm[1]
        @test isapprox(sqrt(4π) .* ds.shell_ylm[1, :], ds.radial; rtol = 1e-12)
        # against the operator without shells: absorbed power to round-off, harmonics to
        # the refined panels' angular placement (2.2e-7 of c₀₀ measured)
        @test isapprox(sum(ds.per_beam), sum(dep.per_beam); rtol = 1e-12)
        @test maximum(abs, ds.ylm .- dep.ylm) < 1e-6 * dep.ylm[1]
        @test isempty(dep.radial) && size(dep.shell_ylm) == (31^2, 0)
        @test_throws ArgumentError CoronaOperator(prof, beams; shells = [rc, R])
    end

    @testset "drive proxies" begin
        # the product quadrature reproduces band-limited coefficients
        c = randn(MersenneTwister(5), 49)
        c[1] = 10.0
        @test maximum(abs, Briareus._sphere_transform((v, μ) -> v, c, 6) .- c) < 1e-12   # 8.4e-14
        # Ba24 Eqs. 2a–2b by hand: an axisymmetric (l ≤ 1) deposition in three shells,
        # weights (n_e/n_c)^{2/3}/r_c² at the shells' middles, Legendre projection of
        # 24.7 F^{2/3} by 1-D quadrature
        e = [rc - L, rc, rc + L, R]
        o3 = CoronaOperator(prof, beams; ncone = 1, lmax = 8, shells = e)
        S(t) = (M = zeros(typeof(t), 81, 3); M[1, :] .= (1.0, 2.0, 0.5); M[3, :] .= (t, -0.3, 0.1); M)
        P = ablation_pressure_proxy(o3, (; shell_ylm = S(0.2)))
        w = exp.(-2 / 3 .* ([rc - L / 2, rc + L / 2, (rc + L + R) / 2] .- rc) ./ L) ./ rc^2
        F(μ) = sum(w .* (S(0.2)[1, :] ./ sqrt(4π) .+ S(0.2)[3, :] .* sqrt(3 / 4π) .* μ))
        μ, wμ = _co_gauss01(200)
        μ, wμ = 2 .* μ .- 1, 2 .* wμ
        ref = [2π * sum(wμ .* 24.7 .* F.(μ) .^ (2 / 3) .* sqrt((2l + 1) / 4π) .* _co_P.(l, μ)) for l in 0:8]
        @test maximum(abs, P[[mode_index(l, 0) for l in 0:8]] .- ref) < 1e-12 * ref[1]   # 4.1e-15
        @test maximum(abs, P[[mode_index(l, m) for l in 1:8 for m in -l:l if m != 0]]) < 1e-12 * ref[1]
        # Ye21 Eq. 3 by hand: E = A Y₀₀ + B Y₁₀ + C Y₁₁, Ē(θ) = A Y₀₀ + B Y₁₀, and a mean
        # displacement d̄(θ); hand-written harmonics on a finer product grid
        A, B, C = 1.0, 0.15, 0.2
        E(t) = (v = zeros(typeof(t), 81); v[1] = A; v[3] = B; v[mode_index(1, 1)] = t; v)
        dbar(θ) = 1 + 0.05cosd(2θ)
        d = shell_displacement_proxy(E(C), 8; mean_displacement = dbar)
        x, wx = Briareus._gauss_legendre(100)
        φs = 2π .* (0:199) ./ 200
        Yh = ((0, 0) => (μ, φ) -> 1 / sqrt(4π), (1, 0) => (μ, φ) -> sqrt(3 / 4π) * μ,
              (1, 1) => (μ, φ) -> sqrt(3 / 4π) * sqrt(1 - μ^2) * cos(φ),
              (2, 0) => (μ, φ) -> sqrt(5 / 16π) * (3μ^2 - 1),
              (2, 1) => (μ, φ) -> sqrt(15 / 4π) * μ * sqrt(1 - μ^2) * cos(φ),
              (2, 2) => (μ, φ) -> sqrt(15 / 16π) * (1 - μ^2) * cos(2φ))
        Ē(μ) = A / sqrt(4π) + B * sqrt(3 / 4π) * μ
        dh(μ, φ) = dbar(acosd(μ)) * ((Ē(μ) + C * sqrt(3 / 4π) * sqrt(1 - μ^2) * cos(φ)) / Ē(μ))^0.72
        errs = [abs(sum(wx[i] * 2π / 200 * dh(x[i], φ) * Y(x[i], φ) for i in eachindex(x), φ in φs) -
                    d[mode_index(l, m)]) for ((l, m), Y) in Yh]
        @test maximum(errs) < 1e-12                          # 1.4e-14 measured
        # an axisymmetric E gives d = d̄ (here 1)
        d0 = shell_displacement_proxy(E(0.0), 8)
        @test isapprox(d0[1], sqrt(4π); rtol = 1e-13) && maximum(abs, d0[2:end]) < 1e-13
        # derivatives through both, against central differences
        fP(t) = mode_rms(ablation_pressure_proxy(o3, (; shell_ylm = S(t))), 8)
        fd(t) = mode_rms(shell_displacement_proxy(E(t), 8; mean_displacement = dbar), 8; m = :nonzero)
        @test all(isapprox(ForwardDiff.derivative(g, 0.2), (g(0.2 + 1e-6) - g(0.2 - 1e-6)) / 2e-6;
                           rtol = 1e-7) for g in (fP, fd))      # 5e-11, 3e-11 measured
        # errors: no shells; a deposition that leaves part of the sphere dark
        @test_throws ArgumentError ablation_pressure_proxy(op, dep)
        @test_throws ArgumentError ablation_pressure_proxy(o3, deposit(o3, beams))
        @test_throws ArgumentError shell_displacement_proxy(dep.ylm, 30)
    end
end
