# Dante channel responses from their hardware and the Dante unfold: Henke attenuation
# against NIST, Fresnel mirrors, the XRD model against the digitized Rubery 2023 curves,
# channel 2's published calibration, and Planckians unfolded in closed form. Research
# notes, "Dante channel responses and the unfold (2026-10-24)".

using Briareus, Test, LinearAlgebra
using Briareus: _rubery_dante, _response, _DANTE_GRID, _planck_shape, _channel_matrix,
                RadiationSolution, SIGMA_SB, SpectralSolution

@testset "Dante channel responses and unfold" begin
    @testset "Henke attenuation and mirrors" begin
        # photoabsorption (Henke) against NIST (with scattering) at 3 keV
        for m in (:Be, :Al, :Ti, :Fe, :Cu, :polyethylene)
            @test henke_mass_attenuation(m, 3000.0) ≈ mass_attenuation(m, 3.0; table = :nist) rtol = 0.03
        end
        @test mass_attenuation(:Al, 3.0) == mass_attenuation(:Al, 3.0; table = :nist)   # :auto keeps NIST
        @test mass_attenuation(:Al, 1.0) == henke_mass_attenuation(:Al, 1000.0)          # and falls back
        @test henke_mass_attenuation(:Al, 1560.0) / henke_mass_attenuation(:Al, 1559.0) > 9   # sharp K edge
        # a compound is the mass-weighted sum of its elements
        μ = (2 * 12.011 * henke_mass_attenuation(:C, 800.0) + 4 * 1.008 * henke_mass_attenuation(:H, 800.0)) /
            (2 * 12.011 + 4 * 1.008)
        @test henke_mass_attenuation(:polyethylene, 800.0) ≈ μ rtol = 1e-12
        @test filter_transmission([:Mg => 0.022], 1.2; table = :henke) ≈
              exp(-henke_mass_attenuation(:Mg, 1200.0) * 1.738 * 0.0022)
        # total external reflection below the critical angle, a cutoff above it, roughness
        @test mirror_reflectivity(:SiO2, 0.3, 1000.0) > 0.9
        @test mirror_reflectivity(:SiO2, 3.5, 2000.0) < 1e-3
        @test mirror_reflectivity(:SiO2, 2.0, 300.0; roughness = 4) < mirror_reflectivity(:SiO2, 2.0, 300.0)
        @test_throws ArgumentError henke_mass_attenuation(:Xx, 1000.0)
    end

    @testset "XRDs and channel 2 against Rubery 2023" begin
        d = _rubery_dante()
        for (n, c, tol) in (("cr_xrd", :Cr, 0.13), ("al_xrd", :Al, 0.08))
            E, s = d[n]
            k = [i for i in eachindex(E) if E[i] > 75 && !(1540 < E[i] < 1750)]
            @test sqrt(sum(abs2, log.(xrd_sensitivity.(c, E[k]) ./ 1e6 ./ s[k])) / length(k)) < tol
        end
        ch = dante1_responses()
        @test [c.number for c in ch] == 1:18
        @test ch[3].mirror == (:SiO2, 3.5)
        @test ch[3].grid == 0.088
        @test ch[16].solid_angle == 1.77e-6
        # duplicated hardware (Rubery: 1/8 differ by the grid, 10/14 identical)
        @test ch[8].R ≈ 0.088 .* ch[1].R
        @test ch[10].R == ch[14].R
        # channel 2: 0.97 TW/sr per volt for power flat over 2–4 keV (Rubery 2023 §II.D)
        E = _DANTE_GRID
        V = ch[2].solid_angle * sum(_response(ch[2], E[k]) * 0.5 * (E[k+1] - E[k]) for k in 1:length(E)-1
                                    if 2000 <= E[k] < 4000)
        @test 1 / V ≈ 0.97 rtol = 0.03
        # and its response against Fig. 8b's (mirror × 47 µm polyimide × Al XRD)
        Eg, g = d["mlm_system"]
        sel = 2000 .< Eg .< 4000
        r = [_response(ch[2], e) / 1e3 for e in Eg[sel]]
        @test sum(r) / sum(g[sel]) ≈ 1 atol = 0.08
        # each band peaks within 25 % of its nominal energy (16 and 18, broad, excepted)
        for c in ch
            c.number in (11, 16, 18) && continue
            @test c.E[argmax(c.R)] / c.energy ≈ 1 atol = 0.25
        end
        @test dante1_channels(model = :gaussian)[3] == DanteChannel(250.0, 75.0)
        @test dante1_channels() === dante1_responses()
    end

    @testset "unfold: Planckians in closed form" begin
        ch = dante1_responses()
        A, c = π * 1.55^2, cosd(37)
        E = _DANTE_GRID
        w = [k == 1 ? (E[2] - E[1]) / 2 : k == length(E) ? (E[end] - E[end-1]) / 2 : (E[k+1] - E[k-1]) / 2
             for k in eachindex(E)]
        for T in (100.0, 200.0, 300.0)
            S = 1e3 * A * c / π * SIGMA_SB * T^4 .* _planck_shape.(E, T)
            V = [h.solid_angle * sum(w .* h.R .* S) for h in ch]
            # the unfold integrates 50 eV – 20 keV: the Planckian's share of it
            Tx = T * planck_fraction(T, 50.0, 20_000.0)^0.25
            u = dante_unfold(ch, V; area = A, cos = c, method = :smooth)
            @test u.start_temperature ≈ T rtol = 1e-3
            @test u.Tr ≈ Tx rtol = 1e-6
            @test u.misfit < 1e-9
            @test all(>=(0), u.spectrum)
            u0 = dante_unfold(ch, V; area = A, cos = c, mband = (0.0, 2500.0, 600.0))
            @test u0.Tr ≈ Tx rtol = 1e-6
            T >= 200 || continue
            # Harris's start (a 1 GW/sr/eV M-band) is unfolded away to < 0.05 eV
            ut = dante_unfold(ch, V; area = A, cos = c)
            @test ut.Tr ≈ Tx atol = 0.2
            T == 300 && @test abs(ut.Tr - Tx) < 0.05
            T == 300 && @test ut.misfit < 1e-6
        end
        @test_throws ArgumentError dante_unfold(ch, zeros(18); area = A)
    end

    @testset "synthetic Dante through the unfold" begin
        t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.6, capsule_radius = 0.9,
                     h = 0.9, capsule_h = 0.45)
        m = RadiationModel(t)
        n = nelements(t)
        T0 = 250.0
        B = [c.kind === :aperture ? 0.0 : SIGMA_SB * T0^4 for c in t.components[t.elem_comp]]
        z = zeros(n)
        sol = RadiationSolution(m, z, z, B, z, z, z, z, 0.0)
        # every element at T0: each group carries its Planck fraction of B
        Eg = spectral_edges(n = 60)
        Bg = B .* planck_fractions(T0, Eg)'
        Z = zeros(size(Bg))
        s = SpectralSolution(sol, Eg, Z, Bg, Z, Z, fill(T0, n), z, zeros(length(Eg) - 1), false)
        r = dante(dante1(), s)
        @test length(r.channels) == 18
        @test all(>(0), r.channels)
        @test r.Tr_unfold ≈ r.Tr * planck_fraction(T0, 50.0, 20_000.0)^0.25 rtol = 2e-3
        rs = dante(dante1(), s; unfold = :smooth)
        @test rs.Tr_unfold ≈ r.Tr * planck_fraction(T0, 50.0, 20_000.0)^0.25 rtol = 2e-3
        @test !haskey(dante(dante1(), s; unfold = nothing), :unfold)
        # Sentinel on the polar DIM: down the axis through the upper hole, onto the capsule
        # and out of the lower hole (same radius)
        sp = dante(sentinel(), sol)
        @test sentinel().through === :leh_top
        @test sentinel().direction[3] ≈ 1
        @test sp.cos ≈ 1
        @test sp.Tr ≈ T0 * (0.9^2 / 1.6^2)^0.25 rtol = 0.03     # the capsule; the rest is the lower hole
        @test_throws ArgumentError sentinel(:equatorial_78)
        @test sentinel(:equatorial_315; through = :leh_top).phi == 315
        g = dante(dante1(channels = dante1_channels(model = :gaussian)), s)
        @test !haskey(g, :unfold)
        @test g.color_temperature ≈ T0 rtol = 0.02
    end
end
