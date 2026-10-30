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
        @test m.materials == [:Au, :DU]   # gold, and the gold liner of a NIF uranium wall
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
    # (iterated albedos: the GMRES path). The walls differ in material only, so they
    # share the 192-beam deposition operator and the view factors.
    coarse(; kw...) = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.6, capsule_radius = 0.9,
                               h = 0.9, capsule_h = 0.45, kw...)
    beams = hybride_mock_pointing(nif_beams(power = 2.0))
    op = DepositionOperator(coarse(), beams)
    F0 = view_factors(coarse())
    for wall in (gold(albedo = 0.85, conversion = 0.8), gold(albedo = WallLossAlbedo(), conversion = 0.8))
        t = coarse(; wall)
        m = RadiationModel(t; F = F0, closed = true)               # ≡ RadiationModel(t)
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
                       (; mband = MBandEmission(), thermalize = true),
                       (; mband = MBandEmission(), spot = :brightness),
                       (; mband = MBandEmission(), spot = :wall, thermalize = true))
                s = spectral_solve(sol; edges = E, kw...)
                @test vec(sum(s.source; dims = 2)) ≈ sol.source rtol = 1e-13
                @test maximum(abs.(vec(sum(s.B; dims = 2)) .- sol.B)) < 1e-11 * maximum(sol.B)
                @test maximum(abs.(vec(sum(s.H; dims = 2)) .- sol.H)) < 1e-11 * maximum(sol.H)
                @test sum(absorbed(s)) ≈ sum(gray_abs) rtol = 1e-11
                @test all(s.B .>= -1e-12 * maximum(sol.B))
            end
            # thermalized, the totals are the gray ones whatever the group albedos
            k = findfirst(==(1800.0), E)
            low = Dict(c => [g >= k ? 0.2 : nothing for g in 1:G] for c in (:wall, :cap_top, :cap_bottom))
            s = spectral_solve(sol; edges = E, albedo = low, thermalize = true, mband = MBandEmission())
            @test maximum(abs.(vec(sum(s.B; dims = 2)) .- sol.B)) < 1e-9 * maximum(sol.B)
            @test all(s.albedo[elements(t, :wall), k:end] .== 0.2)
            # lost instead, the hard groups lose energy
            s = spectral_solve(sol; edges = E, albedo = low)
            @test sum(s.B[cap, :]) < sum(sol.B[cap])
        end

        @testset "laser-spot spectrum ($(wall.albedo isa WallLossAlbedo ? "wall loss" : "constant"))" begin
            E = spectral_edges(n = 16)
            mb = MBandEmission()
            sh = Briareus._band_shares(E, mb.band)
            Pf(T) = reduce(hcat, [planck_fractions(x, E) for x in T])'
            lit = findall(>(0), sol.laser)
            dark = findall(iszero, sol.laser)
            @test !isempty(lit) && !isempty(dark)
            # the earlier rule: everything at the brightness temperature of B
            sb = spectral_solve(sol; edges = E, mband = mb, spot = :brightness)
            TB = radiation_temperature.(sol.B)
            @test sb.temperature == sb.spot_temperature == TB
            @test sb.source ≈ (sol.source .- sb.mband) .* Pf(TB) .+ sb.mband .* sh' rtol = 1e-13
            # converted x-rays less the M-band at the radiation temperature, the rest at the
            # wall temperature (B − ηP/A); here the source is all converted laser
            s = spectral_solve(sol; edges = E, mband = mb)
            @test s.spot === :radiation
            @test s.converted ≈ sol.source rtol = 1e-14
            @test s.spot_temperature ≈ radiation_temperature.(sol.H) rtol = 1e-14
            @test s.temperature ≈ radiation_temperature.(sol.B .- sol.source) rtol = 1e-12
            @test s.source ≈ (sol.source .- s.mband) .* Pf(s.spot_temperature) .+ s.mband .* sh' rtol = 1e-12
            # spots (B > H) are cooler than their brightness temperature, so softer
            k = findfirst(==(1800.0), E)
            spots = findall(e -> sol.laser[e] > 0 && sol.B[e] > sol.H[e], 1:n)
            @test !isempty(spots)
            hot = vec(sum(sb.source[spots, k:end]; dims = 2))
            @test all(vec(sum(s.source[spots, k:end]; dims = 2)) .<= hot .* (1 + 1e-12))
            # without laser light the rules agree
            sw = spectral_solve(sol; edges = E, mband = mb, spot = :wall)
            @test s.B[dark, :] != sb.B[dark, :]          # (the spots' light reaches them)
            @test s.source[dark, :] == sb.source[dark, :] == sw.source[dark, :]
            @test sw.spot_temperature == sw.temperature
            # a stored-energy emission δ on the wall (a re-emitting wall of drive_history)
            # is split at the wall temperature, the converted x-rays as before
            δ = zeros(n)
            δ[elements(t, :wall)] .= 0.1 * maximum(sol.source)
            sol2 = RadiationSolution(m, sol.laser, sol.source .+ δ, sol.B .+ δ, sol.H, sol.albedo,
                                     sol.intact, sol.transmission, sol.time)
            s2 = spectral_solve(sol2; edges = E, mband = mb)
            @test s2.converted ≈ sol.source rtol = 1e-14
            @test s2.temperature ≈ radiation_temperature.(sol.B .- sol.source .+ δ) rtol = 1e-12
            @test s2.source ≈ (sol.source .- s2.mband) .* Pf(s2.spot_temperature) .+ δ .* Pf(s2.temperature) .+
                              s2.mband .* sh' rtol = 1e-12
            @test vec(sum(s2.source; dims = 2)) ≈ sol2.source rtol = 1e-13
            @test_throws ArgumentError spectral_solve(sol; edges = E, spot = :plasma)
        end

        @testset "group solves against dense solves" begin
            E = [0.0, 500.0, 1800.0, 4000.0, Inf]
            G = 4
            T = radiation_temperature.(sol.B) .* (1 .+ 0.1 .* sin.(1:n))   # any fixed split
            a = Dict(:wall => [nothing, 0.7, 0.3, 0.1], :capsule => [0.2, 0.25, 0.3, 0.35])
            M(αg) = I - Diagonal(αg) * m.F
            # in group: B_g = (I - diag(α_g) F)⁻¹ S_g
            s = spectral_solve(sol; edges = E, albedo = a, temperature = T, mband = MBandEmission())
            for g in 1:G
                @test s.B[:, g] ≈ M(s.albedo[:, g]) \ s.source[:, g] rtol = 1e-10
            end
            @test s.H ≈ m.F * s.B
            # a NaN group albedo (computed, say) throws: `nothing` is the gray albedo
            @test_throws ArgumentError spectral_solve(sol; edges = E, albedo = Dict(:wall => [NaN, 0.7, 0.3, 0.1]))
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
            # one group at the gray albedo (taken as the gray radiosity less the other), dense
            E2 = [0.0, 1800.0, Inf]
            s = spectral_solve(sol; edges = E2, albedo = Dict(:wall => [nothing, 0.5]), temperature = T,
                               thermalize = true, mband = MBandEmission())
            b2 = reduce(hcat, [planck_fractions(T[e], E2) for e in 1:n])'
            K2 = zeros(2n, 2n)
            for g in 1:2, h in 1:2
                blk = (g == h ? Diagonal(s.albedo[:, g]) : zero(Diagonal(α))) +
                      Diagonal(b2[:, g] .* (α .- s.albedo[:, h]))
                K2[(g-1)*n+1:g*n, (h-1)*n+1:h*n] .= blk * m.F
            end
            @test vec(s.B) ≈ (I - K2) \ vec(s.source) rtol = 1e-9
            @test vec(sum(s.B; dims = 2)) ≈ sol.B rtol = 1e-12
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

        wall.albedo isa WallLossAlbedo && continue
        @testset "drive history M-band" begin
            # outer → inner transfer growing in time: the cone balance changes
            cbet = CalibratedTransfer((:outer44, :outer50), (:inner23, :inner30), (Δλ, τ) -> 0.15τ)
            pulse = PulseShape([0.0, 1.0, 5.0], [0.0, 2.0, 2.0])
            beams_at = τ -> with_pulses(beams, pulse, τ)
            mb = nif_mband()
            times = [1.0, 3.0]
            h0 = drive_history(m, op, beams_at, times; lpi = cbet, mband = nothing)
            h = drive_history(m, op, beams_at, times; lpi = cbet)          # MBandSplit() by default
            @test h0.mband === nothing && h0.mband_moments === nothing
            # the gray history is untouched, bit for bit
            @test all(getfield(h, k) == getfield(h0, k) for k in keys(h0) if !(k in (:mband, :mband_moments)))
            # the default split is the calibrated pairing, spelled out (α_M as a number and
            # as the group albedos it stands for)
            E2 = [0.0, 1800.0, Inf]
            αM = Dict(c => [nothing, 0.55] for c in (:wall, :cap_top, :cap_bottom))
            @test mband_albedo(t, mb, E2) == αM
            sd = MBandSplit()
            @test sd.emission.groups == mb.groups && sd.edges == E2 && sd.thermalize && sd.albedo == 0.55
            hs = drive_history(m, op, beams_at, times; lpi = cbet,
                               mband = MBandSplit(mb; edges = E2, thermalize = true, albedo = αM))
            @test hs.mband_moments == h.mband_moments && hs.mband == h.mband
            @test_throws ArgumentError mband_albedo(t, mb, [0.0, 1000.0, Inf])
            @test_throws ArgumentError MBandSplit(mb; edges = [0.0, 1000.0, Inf])
            @test_throws ArgumentError MBandSplit(mb; albedo = Dict(:wall => [nothing, 0.5, 0.5]))
            # off where there is nothing to convert: no beam of the model's groups, or no emitter
            @test drive_history(m, op, beams_at, times[1:1];
                                mband = MBandSplit(MBandEmission(groups = Dict(:none => 0.2)))).mband === nothing
            @test drive_history(m, op, beams_at, times[1:1];
                                mband = MBandSplit(nif_mband(materials = (:Pt,)))).mband_moments === nothing
            # reference: a static solve at each time and the M-band summed cone by cone by
            # hand from the effective (post-transfer) beams
            au = vcat(elements(t, :wall), elements(t, :cap_top), elements(t, :cap_bottom))
            ref = map(times) do τ
                eff = apply_lpi(cbet, beams_at(τ), τ)
                P = zeros(n)
                for (g, f) in mb.groups
                    d = deposit(op, [b.group === g ? b : remake(b; power = 0.0) for b in eff])
                    P[au] .+= f .* d.element[au]
                end
                s = spectral_solve(solve(m, deposit(op, eff); time = τ); edges = E2, thermalize = true,
                                   albedo = αM, mband = mb, mband_power = P)
                (; a = legendre_moments(band_flux(s, 1800.0), t, cap; lmax = 6),
                   PM = sum(s.mband .* t.elem_area))
            end
            @test h.mband_moments ≈ reduce(hcat, [r.a for r in ref]) rtol = 1e-10
            @test sum(h.mband[k] for k in keys(h.mband)) ≈ [r.PM for r in ref] rtol = 1e-12
            # the transfer to the low-conversion inner cones lowers the M-band per laser power
            PL(k) = sum(h.laser[c][k] for c in keys(h.laser))
            @test h.mband[:wall][2] / PL(2) < h.mband[:wall][1] / PL(1)
            # with every cone at fraction 1 the per-cone sum is the deposition
            eff = apply_lpi(cbet, beams_at(3.0), 3.0)
            one_f = MBandEmission(groups = Dict(g => 1.0 for g in keys(mb.groups)))
            @test mband_power(one_f, op, eff; lpi = cbet, time = 3.0)[au] ≈ deposit(op, eff).element[au] rtol = 1e-12
            # each group through the operator restricted to its beams: the sums of every
            # beam deposited with the others at zero power, exactly
            Pz = zeros(n)
            for (g, f) in mb.groups
                d = lpi_deposit(cbet, op, [b.group === g ? b : remake(b; power = 0.0) for b in eff], 3.0)
                Pz[au] .+= f .* d.element[au]
            end
            @test mband_power(mb, op, eff; lpi = cbet, time = 3.0) == Pz
            # the split pass, with and without an absorption weight: the deposition is
            # `deposit`'s and each column the group deposited alone, exactly
            for η in (nothing, 0.6)
                key, EG = Briareus._mband_split(mb, t, eff)
                d = Briareus._deposit(op, eff, (key, EG), η)
                d0 = deposit(op, eff; absorption = η)
                alone = [deposit(op, [b.group === g ? b : remake(b; power = 0.0) for b in eff]; absorption = η).element
                         for g in keys(mb.groups)]
                @test d.element == d0.element && d.backside == d0.backside && d.per_beam == d0.per_beam
                @test all(EG[:, k] == alone[k] for k in eachindex(alone))
            end
            # the split's keywords belong to the MBandSplit; `spectral` alone has no effect
            @test_throws ArgumentError drive_history(m, op, beams_at, times; spectral = (; edges = [0.0, Inf]))
            @test_throws ArgumentError drive_history(m, op, beams_at, times; mband = nothing, spectral = (; spot = :wall))
            @test_throws ArgumentError drive_history(m, op, beams_at, times; mband = mb)
        end
    end

    @testset "synthetic Dante" begin
        t = coarse()
        m = RadiationModel(t; F = F0, closed = true)
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

    @testset "Dante field of view" begin
        # the blackbody cavity of "synthetic Dante": a 60 mm field (a NIF Dante's) sees
        # the cavity only through the hole (the outside of every one-sided surface is
        # dark), so it reads what the view through the hole reads, to pixel resolution
        t = coarse()
        n = nelements(t)
        T0 = 250.0
        B = [c.kind === :aperture ? 0.0 : SIGMA_SB * T0^4 for c in t.components[t.elem_comp]]
        z = zeros(n)
        sol = RadiationSolution(RadiationModel(t; F = zeros(n, n)), z, z, B, z, z, z, z, 0.0)
        hole, wide = dante(dante1(), sol), dante(dante1(field_of_view = 60.0, pixels = 800), sol)
        @test wide.area == hole.area && wide.cos == hole.cos        # the same reference
        @test wide.intensity ≈ hole.intensity rtol = 2e-3
        @test Set(keys(visible_area(dante1(field_of_view = 60.0, pixels = 800), t))) ==
              Set(keys(visible_area(dante1(), t)))
        @test_throws ArgumentError DanteView(143, 274; through = nothing)
        @test isnan(dante(dante1(field_of_view = 60.0, through = nothing), sol).Tr)
        @test_throws ArgumentError dante(DanteView(143, 274; through = :window), sol)

        # a flat Lambertian disk facing the detector: I = B/π × the area seen
        R, Bd = 2.0, 1.0
        disk = Component(:emitter, disk_mesh(R, 0.0; nr = 4, nphi = 32), gold())
        lone = Target(disk)
        # half of it hidden by a dark two-sided sheet over y > 0 (the pixel grid is
        # symmetric about the disk's center, and so is the 32-gon)
        shield = foil(:shield, (0.0, 3R / 2, 1.0), (0, 0, 1); width = 3R, n = 6)
        shielded = Target(disk, shield)
        sol_of(tt) = (m = nelements(tt);
                      Bt = [tt.components[tt.elem_comp[e]].name === :emitter ? Bd : 0.0 for e in 1:m];
                      RadiationSolution(RadiationModel(tt; F = zeros(m, m)), zeros(m), zeros(m), Bt,
                                        zeros(m), zeros(m), zeros(m), zeros(m), 0.0))
        A32 = 16 * R^2 * sin(2π / 32)                          # the 32-gon's area
        top = DanteView(0, 0; through = nothing, field_of_view = 100.0, pixels = 1000)
        full = dante(top, sol_of(lone))
        @test full.intensity ≈ Bd / π * A32 rtol = 2e-3
        # against the same sheet behind the disk (the same extent, so the same rays)
        behind = Target(disk, foil(:shield, (0.0, 3R / 2, -1.0), (0, 0, 1); width = 3R, n = 6))
        unblocked = dante(DanteView(0, 0; through = nothing, field_of_view = 100.0, pixels = 1000), sol_of(behind))
        half = dante(DanteView(0, 0; through = nothing, field_of_view = 100.0, pixels = 1000), sol_of(shielded))
        @test unblocked.intensity ≈ full.intensity rtol = 2e-3
        @test half.intensity ≈ unblocked.intensity / 2 rtol = 1e-12
        # Tr over the part's visible projected area: the emitter's own temperature, exactly
        Tb = (π * Bd / π / SIGMA_SB)^0.25
        for tt in (lone, shielded)
            r = dante(DanteView(0, 0; through = nothing, field_of_view = 100.0, pixels = 200, area = :emitter),
                      sol_of(tt))
            @test r.Tr ≈ Tb rtol = 1e-12
        end
        # a field narrower than the disk: the field's disk is what is seen
        D = 2.0
        clip = dante(DanteView(0, 0; through = nothing, field_of_view = D, pixels = 1000), sol_of(lone))
        @test clip.intensity ≈ Bd / π * π * (D / 2)^2 rtol = 2e-3
        off = DanteView(0, 0; through = nothing, field_of_view = D, pixels = 400, center = (0.0, 0.0, 5.0))
        @test dante(off, sol_of(lone)).intensity ≈ clip.intensity rtol = 1e-2   # along the line of sight
        side = DanteView(0, 0; through = nothing, field_of_view = D, pixels = 400, center = (10.0, 0.0, 0.0))
        @test dante(side, sol_of(lone)).intensity == 0
        @test dante(DanteView(0, 0; through = nothing, field_of_view = D, pixels = 400, center = :emitter),
                    sol_of(lone)).intensity ≈ clip.intensity rtol = 1e-2
    end

    @testset "drive history record" begin
        t = hohlraum(nz = 10, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12,
                     wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
        q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
        pulse = PulseShape([0.0, 1.0, 5.0, 6.0], [0.0, 2.0, 8.0, 0.0])
        v = dante1()
        m, op = RadiationModel(t), DepositionOperator(t, q)
        h = drive_history(m, op, τ -> with_pulses(q, pulse, τ), [0.5, 2.0, 4.0]; record = sol -> dante(v, sol).Tr)
        @test length(h.record) == 3 && issorted(h.record)
        @test drive_history(m, op, τ -> with_pulses(q, pulse, τ), [1.0]).record === nothing
    end
end
