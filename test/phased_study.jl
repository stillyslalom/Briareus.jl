# What a multi-phase re-pointing study writes on the phased scenario: a drive band, gates
# at a fraction of the noise floor (or physical gates above it), a shifted target shape,
# keep-outs over every phase, the power scales inside an objective and in a table, a
# phase from a drive history, the summary against gates and floors, and a picket ceiling
# whose reference sits on another target.
#
# References: the quantities by hand (mean flux, keep-out clearances per phase, the
# pulses' powers, the history's state and the scenario remade from it); the Legendre
# moments against the spherical-harmonic coefficients on the same sub-triangles for the
# shift's units; the closed-form gradient of a quadratic in the scales.

@testset "phased study API" begin
    geom = (; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9, h = 1.2, capsule_h = 0.9)
    wall(a) = Material(:Au; albedo = a, conversion = 0.8)
    tP, tK = hohlraum(; geom..., wall = wall(0.6)), hohlraum(; geom..., wall = wall(0.85))   # picket, later phases
    aims = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20), :inner30 => (0.60, 3.63), :inner23 => (1.09, 3.09))
    quads = map(nif_quads(power = 1.0)) do q
        r, z = aims[q.group]
        φ = port_angles(q)[2]
        remake(q; aim = (-r * cosd(φ), -r * sind(φ), sign(q.port[3]) * z))
    end
    inner(b) = b.group in (:inner23, :inner30)
    gq = GroupPointing(quads; by = :quad, fields = (:z, :phi), relative = true, mirror = false)
    dk = (; resolution = 0.6, reach = 1.0, level = 0)
    sP = Scenario(tP, gq; observe = :capsule, viewfactor_kw = (; closed = true), deposition_kw = dk)
    vf = (; F = sP.model.F, closed = true)
    pulses = Dict(:inner23 => PulseShape([0.0, 2.0, 4.0], [0.1, 0.3, 2.0]),
                  :inner30 => PulseShape([0.0, 2.0, 4.0], [0.1, 0.35, 2.0]),
                  :outer44 => PulseShape([0.0, 2.0, 4.0], [0.4, 0.8, 2.0]),
                  :outer50 => PulseShape([0.0, 2.0, 4.0], [0.4, 0.7, 2.0]))
    picket_power(b) = inner(b) ? 0.1 : 0.4
    ps = PhasedScenario(gq, [Phase(:picket, sP; power = picket_power),
                             Phase(:step, tK; pulses, time = 2.5, viewfactor_kw = vf),
                             Phase(:peak, tK; power = 2.0, scaled = false, viewfactor_kw = vf)];
                        by = :quad, observe = :capsule, deposition_kw = dk)
    pp = ps.param
    nq = nparams(gq)
    y = vcat(0.15 .* sin.(1:nq), 1 .+ 0.1 .* cos.(1:length(pp.keys)))
    e0, ev = evaluate(ps, initial_parameters(pp)), evaluate(ps, y)
    backscatter(s) = remake(s; lpi = Backscatter(:outer50 => 0.5))

    @testset "a drive band" begin
        lo, hi = drive_floor(e0[:step]; fraction = 1.0), drive_ceiling(e0[:step]; fraction = 1.02)
        @test hi.upper && !lo.upper && hi.name === :drive_ceiling
        @test only(hi.limit) == 1.02 * mean_flux(e0[:step]) && only(hi.scale) == only(hi.limit)
        @test hi.f(ev[:step]) == mean_flux(ev[:step])
        rows = check_constraints(ev, [in_phase(:step, lo), in_phase(:step, hi)])
        @test [r.name for r in rows] == ["step.drive", "step.drive_ceiling"]
        L = 1.02 * mean_flux(e0[:step])
        @test rows[2].violation ≈ (mean_flux(ev[:step]) - L) / L rtol = 1e-12
        # the reference's laser–plasma model, checked as for the floor
        @test_throws ArgumentError hi.f(evaluate(backscatter(ps[:step]), ev[:step].beams))
        @test_throws ArgumentError drive_ceiling(e0[:step]; fraction = 0.0)
    end

    @testset "gates: the noise floor and a shifted target" begin
        u = (; rms_mode = [0.002, 0.0012, 0.0015, 0.0004])     # a floor as mode_uncertainty reports it
        fl = [sqrt(2l + 1) * r for (l, r) in enumerate(u.rms_mode)]
        @test mode_gate(e0[:picket], u).limit == fl
        @test mode_gate(e0[:picket], u; fraction = 0.8).limit == 0.8 .* fl
        phys = [0.004, 0.004, 0.004, 0.0005]
        @test mode_gate(e0[:picket], phys; floor = u, fraction = 0.8).limit == max.(phys, 0.8 .* fl)
        @test mode_gate(e0[:picket], phys; floor = fl).limit == max.(phys, fl)
        @test mode_gate(flux_ylm(e0[:picket]), u; fraction = 0.5).limit == 0.5 .* fl
        @test_throws ArgumentError mode_gate(e0[:picket], phys; fraction = 0.8)     # no floor to scale
        @test_throws ArgumentError mode_gate(e0[:picket], u; floor = u)             # the floor twice
        @test_throws ArgumentError mode_gate(e0[:picket], u; fraction = -1.0)
        @test_throws DimensionMismatch mode_gate(e0[:picket], phys[1:3]; floor = u)
        # a shifted target: the gate's value is the departure from the shifted shape, and
        # the reference's model is still checked
        c0 = flux_ylm(e0[:peak]; lmax = 4)
        g0, gs = mode_gate(e0[:peak], phys), mode_gate(e0[:peak], phys; shift = 2 => 0.024)
        @test g0.f(ev[:peak]) == mode_departure(ev[:peak], c0)
        @test gs.f(ev[:peak]) == mode_departure(ev[:peak], shift_modes(c0, 2 => 0.024))
        @test mode_gate(c0, phys; shift = 2 => 0.024).f(ev[:peak]) == gs.f(ev[:peak])
        @test_throws ArgumentError gs.f(evaluate(backscatter(ps[:peak]), ev[:peak].beams))
        # the shift in its degree alone, for m = 0 and m ≠ 0; the drive unchanged
        @test mode_departure(shift_modes(c0, 2 => 0.024), c0) ≈ [0, 0.024, 0, 0] atol = 1e-15
        @test mode_departure(shift_modes(c0, [(3, -2) => 0.01, 1 => -0.005]), c0) ≈ [0.005, 0, 0.01, 0] atol = 1e-15
        @test shift_modes(c0, 2 => 0.024)[1] == c0[1]
        # Legendre units: aₗ/a₀ of the capsule flux's Legendre moments is √(2l+1) c_l0/c_00
        # (both integrated on the same sub-triangles), so the shift moves P2/P0 by itself
        p2(c) = sqrt(5) * c[mode_index(2, 0)] / c[1]
        e = e0[:peak]
        obs, te = e.scenario.observed, e.scenario.target
        a = Briareus.legendre_weights(te, obs; lmax = 4, level = 3) * e.flux
        cy = Briareus.ylm_weights(te, obs; lmax = 4, level = 3) * e.flux
        @test p2(cy) ≈ a[3] / a[1] rtol = 1e-10
        @test p2(shift_modes(c0, 2 => 0.024)) - p2(c0) ≈ 0.024 rtol = 1e-12
        @test_throws ArgumentError shift_modes(c0, (2, 3) => 0.01)
        @test_throws ArgumentError shift_modes(c0, 5 => 0.01)                     # past the coefficients
    end

    @testset "keep-outs over every phase" begin
        ko = ApertureKeepOut(tP; minimum = -0.8)
        clf = clearance_floor(ko; name = :leh)
        rows = check_constraints(ev, clf)
        @test length(rows) == length(quads) && all(r -> startswith(r.name, "leh["), rows)
        # a keep-out at a fraction of the peak: every phase's clearance, bit for bit
        @test clf.f(ev) == in_phase(:picket, clf).f(ev) == keepout_clearance(ko, ev[:peak].beams)
        # at an absolute intensity the envelope widens with the power: the brightest phase
        # sets each beam's clearance
        ki = ApertureKeepOut(tP; intensity = 0.05, minimum = -0.8)
        per = [keepout_clearance(ki, e.beams) for e in ev]
        @test clearance_floor(ki).f(ev) == min.(per...) && per[3] != per[1]
        # under the optimizer's duals (partials in an aim and a power scale, so the phases'
        # powers differ in value and partials): the phases the keep-out sees as the first
        # are skipped, and the value is the least over every phase, values and partials
        # bit for bit; at an intensity the third phase is not skipped
        yd = [ForwardDiff.Dual{:ko}(v, Float64(i == 1), Float64(i == length(y))) for (i, v) in enumerate(y)]
        evd = evaluate(ps, yd)
        perd = [keepout_clearance(ko, e.beams) for e in evd]
        got, least = clf.f(evd), reduce((a, b) -> min.(a, b), perd)
        @test all(i -> got[i] === least[i], eachindex(got, least))
        same(k, j) = all(i -> Briareus._same_clearance(k, evd[1].beams[i], evd[j].beams[i]), eachindex(evd[1].beams))
        @test same(ko, 2) && same(ko, 3) && !same(ki, 3)
        gi, leasti = clearance_floor(ki).f(evd), reduce((a, b) -> min.(a, b), [keepout_clearance(ki, e.beams) for e in evd])
        @test all(i -> gi[i] === leasti[i], eachindex(gi, leasti))
    end

    @testset "the power scales in an objective and a table" begin
        s = power_scales(ps, y)
        @test s == y[nq+1:end]
        @test power_scales(ps, apply(pp, y)) ≈ s rtol = 1e-14
        @test power_scales(ps, ev) ≈ s rtol = 1e-14
        # differentiable: the gradient of a quadratic in the scales, in closed form
        G = Briareus._scenario_gradient(ps, e -> sum(abs2, power_scales(ps, e) .- 1), y)
        @test all(iszero, G[1:nq]) && isapprox(G[nq+1:end], 2 .* (s .- 1); rtol = 1e-12)
        # the phases' powers before the scales
        @test phase_powers(ps, :picket) == [picket_power(q) for q in quads]
        @test phase_powers(ps, :step) == [b.power for b in with_pulses(quads, pulses, 2.5)]
        @test phase_powers(ps, :peak) == fill(2.0, length(quads))
        # a key dark in its phase reads 1; a phase on the pointing's own powers has no
        # powers to read the scales against
        dark = PhasedScenario(gq, Phase(:dark, sP; power = b -> inner(b) ? 0.0 : 0.4); by = :quad)
        xd = vcat(zeros(nq), fill(1.3, length(quads)))
        @test power_scales(dark, apply(dark.param, xd)) ≈ [inner(q) ? 1.0 : 1.3 for q in quads] rtol = 1e-14
        own = PhasedScenario(gq, Phase(:own, sP); by = :quad)
        @test_throws ArgumentError power_scales(own, apply(own.param, initial_parameters(own.param)))
        @test power_scales(own, initial_parameters(own.param)) == ones(length(quads))
        # per key and phase
        tb = power_scale_table(ps, y)
        @test tb.keys == [q.quad for q in quads] && tb.phases == [:picket, :step]
        @test tb.scales[:, 1] == s[1:48] && tb.scales[:, 2] == s[49:96]
        @test power_scale_table(ps, ev).scales ≈ tb.scales rtol = 1e-14
        txt = sprint(show, MIME("text/plain"), tb)
        @test occursin("step", txt) && occursin("mean", txt) && occursin(quads[1].quad, txt)
    end

    @testset "the summary against gates and floors" begin
        gates = [in_phase(k, mode_gate(e0[k], fill(0.01, 4); shift = k === :peak ? (2 => 0.01) : nothing)) for k in keys(ps)]
        cs = [gates; in_phase(:step, drive_floor(e0[:step]; fraction = 1.0))]
        floors = PerPhase(:picket => fill(0.002, 4), :step => fill(0.003, 4), :peak => fill(0.004, 4))
        bs = apply(pp, y)
        sm = pointing_summary(ps, bs; drive_reference = e0, gates = cs, floor = floors)
        @test all(k -> sm[k].gate == fill(0.01, 4) && sm[k].floor == floors[k], keys(ps))
        # a gate's own departure (the peak's from its shifted shape), the others' the summary's
        @test sm[:picket].gate_departure ≈ sm[:picket].departure rtol = 1e-12
        @test sm[:peak].gate_departure ≈ mode_departure(ev[:peak], shift_modes(flux_ylm(e0[:peak]; lmax = 4), 2 => 0.01)) rtol = 1e-12
        @test sm[:peak].departure ≈ mode_departure(ev[:peak], flux_ylm(e0[:peak])) rtol = 1e-12
        txt = sprint(show, MIME("text/plain"), sm)
        @test occursin("over the gates", txt) && occursin("over the floor", txt)
        # one gate vector for every phase; a phase without an entry
        s1 = pointing_summary(ps, bs; drive_reference = e0, gates = fill(0.01, 4))
        @test all(k -> s1[k].gate_departure == s1[k].departure, keys(ps))
        s2 = pointing_summary(ps, bs; drive_reference = e0, gates = (; picket = fill(0.01, 4)))
        @test s2[:picket].gate == fill(0.01, 4) && isempty(s2[:step].gate)
        @test_throws ArgumentError pointing_summary(ps, bs; gates = fill(0.01, 4))      # no drive reference
        @test isempty(pointing_summary(ps, bs; gates = cs)[:step].floor)
    end

    @testset "a picket ceiling on another target" begin
        # the design's hohlraum has a window in its wall; the reference's has none
        tW = hohlraum(; geom..., wall = wall(0.6), windows = [WallWindow(phi = 0, diameter = 2.0, material = wall(0.3))])
        sW = Scenario(tW, gq; deposition_kw = dk)
        bP = apply(pp, y)[:picket]
        eW, ref = evaluate(sW, bP), e0[:picket]
        ca = picket_ceiling(ref, :wall, :window_1; reference_names = (:wall,), measure = :average)
        @test only(ca.limit) == picket_irradiance(ref, :wall; measure = :average)
        @test ca.f(eW) == picket_irradiance(eW, :wall, :window_1; measure = :average) == check_constraints(eW, ca)[1].value
        @test_throws ArgumentError ca.f(ref)                      # the reference's target has no window
        # the other integral measures, the same names on both targets
        cw, cc = picket_ceiling(ref, :wall; measure = :weighted), picket_ceiling(ref, :wall; measure = :corridor)
        @test cw.f(eW) == picket_irradiance(eW, :wall; measure = :weighted)
        @test cc.f(eW) == picket_irradiance(eW, :wall; measure = :corridor)
        @test cw.f(ev[:picket]) == picket_irradiance(ev[:picket], :wall; measure = :weighted)
        # the peak holds on the reference's mesh only
        @test_throws ArgumentError picket_ceiling(ref, :wall).f(eW)
        # with lpi = true the reference's laser–plasma model and time are checked
        cl = picket_ceiling(ref, :wall; measure = :average, lpi = true)
        @test cl.f(eW) == picket_irradiance(eW, :wall; measure = :average, lpi = true)
        @test_throws ArgumentError cl.f(evaluate(backscatter(sW), bP))
    end

    @testset "a phase from a drive history" begin
        kw = (; nz = 12, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12)
        th = hohlraum(; kw..., wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
        m = RadiationModel(th; closed = true)
        qh = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
        pulse = PulseShape([0.0, 0.2, 0.6, 0.8, 1.6, 2.0, 3.0], [0.0, 3.0, 3.0, 0.5, 0.5, 8.0, 8.0])
        h = drive_history(m, DepositionOperator(th, qh), τ -> with_pulses(qh, pulse, τ), collect(0.25:0.25:2.5);
                          keep_albedo = true, lmax = 4)
        gh = GroupPointing(qh; by = :group, fields = (:z,), mirror = false)
        lpi = Backscatter(:outer50 => 0.1)        # a model, so that the reference checks compare times
        shared = Scenario(th, gh; observe = :capsule, model = m, lpi, albedo = :reference)   # the phases' from the history
        τ = 1.3
        st = history_state(th, h, τ)
        ph = Phase(:trough, shared; history = h, time = τ, pulses = pulse)
        @test ph.time == τ && ph.source.time == τ && ph.source.albedo == st.albedo
        @test ph.source.moment_response == remake(shared; st.time, st.albedo).moment_response
        psh = PhasedScenario(gh, [ph, Phase(:peak, shared; history = h, time = 2.5, power = 8.0, scaled = false)])
        eh = evaluate(psh, initial_parameters(psh.param))
        # the reference remade at the phase's time passes the gate's time check
        ref = remake(shared; history = h, time = τ)
        rb = phase_beams(psh, qh)[:trough]
        gate = in_phase(:trough, mode_gate(evaluate(ref, rb), fill(0.01, 4)))
        @test gate.f(eh) == mode_departure(eh[:trough], flux_ylm(evaluate(ref, rb); lmax = 4))
        @test_throws ArgumentError in_phase(:trough, mode_gate(evaluate(shared, rb), fill(0.01, 4))).f(eh)
        # a phase built from the target takes the history's albedos the same way
        pt = PhasedScenario(gh, [Phase(:trough, th; history = h, time = τ, pulses = pulse, model = m)];
                            observe = :capsule, lpi)
        @test pt[:trough].albedo == st.albedo && pt[:trough].time == τ
        @test flux_ylm(evaluate(pt, initial_parameters(pt.param))[:trough]) ≈ flux_ylm(eh[:trough]) rtol = 1e-12
        @test_throws ArgumentError Phase(:x, shared; history = h, pulses = pulse)             # no time
        @test_throws ArgumentError remake(shared; history = h)
        @test_throws ArgumentError remake(shared; history = h, time = τ, st.albedo)
        @test_throws ArgumentError Phase(:x, th; history = h, time = τ, pulses = pulse, st.albedo)
    end
end
