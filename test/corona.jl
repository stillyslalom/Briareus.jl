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
end
