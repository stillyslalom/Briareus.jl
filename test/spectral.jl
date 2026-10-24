# Multigroup (spectral) radiosity, the laser-produced M-band and the synthetic Dante:
# Planck fractions against quadrature, the gray solution reproduced by one group and by
# the group sums, in-group and thermalized group solves against dense solves of their
# defining systems, the M-band conversion table, and a blackbody cavity seen through
# its entrance hole. Research notes, "Spectral radiosity, M-band and a synthetic Dante
# (2026-10-24)".

using Briareus, Test, LinearAlgebra
using Briareus: RadiationSolution, _planck_tail, _channel_matrix, _color_temperature

@testset "spectral radiosity" begin
    @testset "Planck fractions" begin
        @test planck_fraction(300.0, 0, Inf) === 1.0
        @test planck_fraction(0.0, 0, 100) == 1.0 && planck_fraction(0.0, 100, Inf) == 0.0
        E = spectral_edges()
        @test E[1] == 0 && E[end] == Inf && 1800.0 in E && issorted(E)
        for T in (30.0, 150.0, 300.0, 1000.0)
            @test sum(planck_fractions(T, E)) ≈ 1 atol = 1e-14
        end
        # against midpoint quadrature of x³/(eˣ - 1), on both sides of the series switch
        f(x) = x^3 / expm1(x)
        cdf(x) = (n = 200_000; h = x / n; sum(f((k - 0.5) * h) for k in 1:n) * h * 15 / π^4)
        for x in (0.1, 0.49, 0.51, 1.0, 3.0, 6.0)
            @test 1 - _planck_tail(x) ≈ cdf(x) rtol = 1e-8
        end
        @test abs(_planck_tail(0.5 - 1e-12) - _planck_tail(0.5 + 1e-12)) < 1e-11
        # Wien tail: (15/π⁴) e⁻ˣ (x³ + 3x² + 6x + 6) for x ≫ 1
        x = 40.0
        @test _planck_tail(x) ≈ 15 / π^4 * exp(-x) * (x^3 + 3x^2 + 6x + 6) rtol = 1e-12
        # 1.8 keV is ~14 % of a 300 eV Planckian
        @test 0.13 < planck_fraction(300.0, 1800, Inf) < 0.15
    end

    @testset "M-band conversion" begin
        m = MBandEmission()
        @test mband_conversion(m, 1e14) ≈ 0.01 && mband_conversion(m, 1e15) ≈ 0.20
        @test mband_conversion(m, 1e13) == 0.01 && mband_conversion(m, 5e15) == 0.20
        @test mband_conversion(m, 3e14) ≈ 0.01 * 3^log10(20) rtol = 1e-12     # log–log
        n = nif_mband()
        @test n.groups[:outer44] == n.groups[:outer50] == 0.20
        @test n.groups[:inner23] == n.groups[:inner30]
        @test n.groups[:inner23] ≈ 0.0418 atol = 1e-4
        @test_throws ArgumentError MBandEmission(fraction = [0.5, 1.5])
    end

    # a coarse hohlraum: constant albedos (cached LU path) and the Lindl wall loss
    # (iterated albedos: the GMRES path)
    for wall in (gold(albedo = 0.85, conversion = 0.8), gold(albedo = WallLossAlbedo(), conversion = 0.8))
        t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.6, capsule_radius = 0.9,
                     h = 0.9, capsule_h = 0.45, wall)
        beams = hybride_mock_pointing(nif_beams(power = 2.0))
        op = DepositionOperator(t, beams)
        m = RadiationModel(t)
        sol = solve(m, deposit(op, beams); time = 5.0)
        n = nelements(t)
        cap = elements(t, :capsule)
        @testset "gray reproduction ($(wall.albedo isa WallLossAlbedo ? "wall loss" : "constant"))" begin
            s1 = spectral_solve(sol)
            @test ngroups(s1) == 1
            @test s1.B == reshape(sol.B, n, 1) && s1.H == reshape(sol.H, n, 1)    # exactly
            E = spectral_edges(n = 16)
            G = length(E) - 1
            absorbed(s) = vec(sum((1 .- s.albedo) .* s.H; dims = 2)) .* t.elem_area
            gray_abs = (1 .- sol.albedo) .* sol.H .* t.elem_area
            for kw in ((;), (; thermalize = true), (; mband = MBandEmission()),
                       (; mband = MBandEmission(), thermalize = true))
                s = spectral_solve(sol; edges = E, kw...)
                @test vec(sum(s.source; dims = 2)) ≈ sol.source rtol = 1e-13
                @test maximum(abs.(vec(sum(s.B; dims = 2)) .- sol.B)) < 1e-11 * maximum(sol.B)
                @test maximum(abs.(vec(sum(s.H; dims = 2)) .- sol.H)) < 1e-11 * maximum(sol.H)
                @test sum(absorbed(s)) ≈ sum(gray_abs) rtol = 1e-11
                @test all(s.B .>= -1e-12 * maximum(sol.B))
            end
            # thermalized, the totals are the gray ones whatever the group albedos
            k = findfirst(==(1800.0), E)
            low = Dict(c => [g >= k ? 0.2 : NaN for g in 1:G] for c in (:wall, :cap_top, :cap_bottom))
            s = spectral_solve(sol; edges = E, albedo = low, thermalize = true, mband = MBandEmission())
            @test maximum(abs.(vec(sum(s.B; dims = 2)) .- sol.B)) < 1e-9 * maximum(sol.B)
            @test all(s.albedo[elements(t, :wall), k:end] .== 0.2)
            # lost instead, the hard groups lose energy
            s = spectral_solve(sol; edges = E, albedo = low)
            @test sum(s.B[cap, :]) < sum(sol.B[cap])
        end

        @testset "group solves against dense solves" begin
            E = [0.0, 500.0, 1800.0, 4000.0, Inf]
            G = 4
            T = radiation_temperature.(sol.B) .* (1 .+ 0.1 .* sin.(1:n))   # any fixed split
            a = Dict(:wall => [NaN, 0.7, 0.3, 0.1], :capsule => [0.2, 0.25, 0.3, 0.35])
            M(αg) = I - Diagonal(αg) * m.F
            # in group: B_g = (I - diag(α_g) F)⁻¹ S_g
            s = spectral_solve(sol; edges = E, albedo = a, temperature = T, mband = MBandEmission())
            for g in 1:G
                @test s.B[:, g] ≈ M(s.albedo[:, g]) \ s.source[:, g] rtol = 1e-10
            end
            @test s.H ≈ m.F * s.B
            # thermalized: the coupled system of all groups, dense
            s = spectral_solve(sol; edges = E, albedo = a, temperature = T, thermalize = true,
                               mband = MBandEmission())
            b = reduce(hcat, [planck_fractions(T[e], E) for e in 1:n])'
            α = sol.albedo
            K = zeros(n * G, n * G)
            for g in 1:G, h in 1:G
                blk = (g == h ? Diagonal(s.albedo[:, g]) : zero(Diagonal(α))) +
                      Diagonal(b[:, g] .* (α .- s.albedo[:, h]))
                K[(g-1)*n+1:g*n, (h-1)*n+1:h*n] .= blk * m.F
            end
            x = (I - K) \ vec(s.source)
            @test vec(s.B) ≈ x rtol = 1e-9
            # all re-emission thermal
            z = Dict(c.name => zeros(G) for c in t.components if c.kind === :surface)
            s = spectral_solve(sol; edges = E, albedo = z, temperature = T, thermalize = true)
            @test s.B ≈ s.source .+ b .* (α .* sol.H)
        end

        @testset "M-band source" begin
            mb = nif_mband()
            PM = mband_power(mb, op, beams)
            wall_e = elements(t, :wall)
            @test all(PM[cap] .== 0) && sum(PM) > 0
            # the per-group sum: each group deposited alone
            P = zeros(n)
            for (g, f) in mb.groups
                d = deposit(op, [b.group === g ? b : remake(b; power = 0.0) for b in beams])
                for e in vcat(wall_e, elements(t, :cap_top), elements(t, :cap_bottom))
                    P[e] += f * d.element[e]
                end
            end
            @test PM ≈ P
            @test_throws ArgumentError mband_power(mb, t, sol.laser)
            s = spectral_solve(sol; edges = spectral_edges(n = 16), mband = mb, mband_power = PM)
            @test s.mband ≈ min.(PM ./ t.elem_area, sol.source)
            # laser-produced M-band is in its band, and F_M is finite
            @test sum(band_flux(s, 1800.0; field = :source)) > sum(s.mband) * 0.999
            @test 0 < mband_recirculation(s, wall_e) < Inf
        end
    end

    @testset "synthetic Dante" begin
        t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.6, capsule_radius = 0.9,
                     h = 0.9, capsule_h = 0.45)
        m = RadiationModel(t)
        n = nelements(t)
        # a blackbody cavity at T0: the hole radiates σT0⁴ over its area
        T0 = 250.0
        B = [c.kind === :aperture ? 0.0 : SIGMA_SB * T0^4 for c in t.components[t.elem_comp]]
        z = zeros(n)
        sol = RadiationSolution(m, z, z, B, z, z, z, z, 0.0)
        # the opening is the caps' inscribed inner polygon (the equal-area hole is π r²)
        rim = unique(round.((p[1], p[2]); digits = 9) for e in elements(t, :cap_top) for k in t.elem_tris[e]
                     for p in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
                     if abs(hypot(p[1], p[2]) - 1.6) < 1e-9)
        Apoly = length(rim) / 2 * 1.6^2 * sin(2π / length(rim))
        for v in (dante1(), dante2(), DanteView(143, 274; pixels = 200))
            r = dante(v, sol)
            hole = v.through === :leh_bottom ? elements(t, :leh_bottom) : elements(t, :leh_top)
            @test r.area ≈ sum(t.elem_area[hole]) ≈ π * 1.6^2
            @test r.Tr ≈ T0 * (Apoly / r.area)^0.25 rtol = 5e-4
            @test r.Tr_open == r.Tr
            @test sum(values(visible_area(v, t))) ≈ Apoly * r.cos rtol = 2e-3
        end
        @test dante1().direction ⋅ [0, 0, -1] ≈ cosd(37) atol = 1e-3   # 37° to the axis
        # a closing lip transmits its open fraction of what lies behind it
        tl = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.6, capsule_radius = 0.9,
                      h = 0.9, capsule_h = 0.45, leh_closure = GapClosure(closure = 5.0),
                      leh_open_radius = 1.2)
        nl = nelements(tl)
        Bl = [c.kind === :aperture ? 0.0 : SIGMA_SB * T0^4 for c in tl.components[tl.elem_comp]]
        zl = zeros(nl)
        for f in (0.0, 0.5, 1.0)
            τ = zeros(nl)
            τ[elements(tl, :leh_bottom_lip)] .= f
            r = dante(dante1(), RadiationSolution(RadiationModel(tl; F = zeros(nl, nl)), zl, zl, Bl, zl, zl, zl, τ, 0.0))
            Ao, Al = (sum(tl.elem_area[elements(tl, c)]) for c in (:leh_bottom, :leh_bottom_lip))
            open = (Ao + f * Al) / (Ao + Al)
            @test r.Tr ≈ T0 * open^0.25 rtol = 4e-3
            f > 0 && @test r.Tr_open ≈ r.Tr / open^0.25 rtol = 1e-12
        end
        # spectrum: a Planckian's channel signals unfold to its temperature
        E = spectral_edges()
        ch = dante1_channels()
        W = _channel_matrix(ch, E)
        for T in (100.0, 200.0, 300.0)
            sig = W * planck_fractions(T, E)
            @test _color_temperature(ch, E, W, sig) ≈ T rtol = 1e-6
        end
        @test all(W[:, end] .== 0)
    end

    @testset "drive history record" begin
        t = hohlraum(nz = 10, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12,
                     wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
        q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
        pulse = PulseShape([0.0, 1.0, 5.0, 6.0], [0.0, 2.0, 8.0, 0.0])
        v = dante1()
        h = drive_history(RadiationModel(t), DepositionOperator(t, q), τ -> with_pulses(q, pulse, τ),
                          [0.5, 2.0, 4.0]; record = sol -> dante(v, sol).Tr)
        @test length(h.record) == 3 && issorted(h.record)
        @test drive_history(RadiationModel(t), DepositionOperator(t, q), τ -> with_pulses(q, pulse, τ),
                            [1.0]).record === nothing
    end
end
