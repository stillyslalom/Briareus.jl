# One pointing over several phases of a pulse (`PhasedScenario`).
#
# References: the plain `Scenario` (one phase at the pointing's own powers, bit for bit);
# the duplicate-beam construction it replaces at the same parameters (a copy of every
# beam per scaled phase under one quad name, so that a per-quad pointing gives the copies
# one aim, a `PowerScales` keyed by (quad, phase) through a counting `by`, and each
# phase's scenario evaluated on its copies); central finite differences for the
# gradient; the beams at the box's upper corner for the laser limit.

@testset "phased scenarios" begin
    geom = (; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9, h = 1.2, capsule_h = 0.9)
    wall(a) = Material(:Au; albedo = a, conversion = 0.8)
    tP, tS, tK = (hohlraum(; geom..., wall = wall(a)) for a in (0.6, 0.8, 0.85))   # picket, step, peak
    aims = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20), :inner30 => (0.60, 3.63), :inner23 => (1.09, 3.09))
    quads = map(nif_quads(power = 1.0)) do q
        r, z = aims[q.group]
        φ = port_angles(q)[2]
        remake(q; aim = (-r * cosd(φ), -r * sind(φ), sign(q.port[3]) * z))
    end
    inner(b) = b.group in (:inner23, :inner30)
    gq = GroupPointing(quads; by = :quad, fields = (:z, :phi), relative = true, mirror = false)
    dk = (; resolution = 0.6, reach = 1.0, level = 0)   # a coarse deposition, for the runs' cost
    sP = Scenario(tP, gq; observe = :capsule, viewfactor_kw = (; closed = true), deposition_kw = dk)
    vf = (; F = sP.model.F, closed = true)        # the step's and the peak's targets differ in albedo only
    pulses = Dict(:inner23 => PulseShape([0.0, 2.0, 4.0], [0.1, 0.3, 2.0]),
                  :inner30 => PulseShape([0.0, 2.0, 4.0], [0.1, 0.35, 2.0]),
                  :outer44 => PulseShape([0.0, 2.0, 4.0], [0.4, 0.8, 2.0]),
                  :outer50 => PulseShape([0.0, 2.0, 4.0], [0.4, 0.7, 2.0]))
    picket_power(b) = inner(b) ? 0.1 : 0.4
    phases() = [Phase(:picket, sP; power = picket_power),
                 Phase(:step, tS; pulses, time = 2.5, viewfactor_kw = vf),
                 Phase(:peak, tK; power = 2.0, scaled = false, viewfactor_kw = vf)]
    ps = PhasedScenario(gq, phases(); by = :quad, observe = :capsule, deposition_kw = dk)
    gg = GroupPointing(quads; by = :group, fields = (:z, :phi), relative = true, mirror = false)   # by cone
    pp = ps.param
    nq = nparams(gq)
    y = vcat(0.15 .* sin.(1:nq), 1 .+ 0.1 .* cos.(1:length(pp.keys)))   # aims and scales
    rel(a, b) = maximum(abs, a .- b) / maximum(abs, b)

    @testset "parameterization" begin
        @test keys(ps) == [:picket, :step, :peak] && ps[:step].target === tS
        @test nparams(pp) == nq + 2 * 48 && pp.n == nq
        @test length(pp.keys) == 96 && pp.keys[1] == (quads[1].quad, :picket) && pp.keys[49] == (quads[1].quad, :step)
        names = parameter_names(pp)
        @test names[1:nq] == parameter_names(gq) && names[nq+1] == "$(quads[1].quad).power_scale.picket"
        bs = apply(pp, y)
        @test bs isa PerPhase{Vector{Beam{Float64}}} && keys(bs) == [:picket, :step, :peak]
        # one set of aims; the phases' powers, scaled where they have scales
        @test all(k -> all(i -> bs[k][i].aim == bs[1][i].aim, eachindex(quads)), 2:3)
        @test [b.power for b in bs[:picket]] ≈ [picket_power(q) * y[nq+i] for (i, q) in enumerate(quads)] rtol = 1e-14
        st = with_pulses(quads, pulses, 2.5)
        @test [b.power for b in bs[:step]] ≈ [st[i].power * y[nq+48+i] for i in eachindex(quads)] rtol = 1e-14
        @test all(b -> b.power == 2.0, bs[:peak])
        @test phase_beams(ps, apply(gq, initial_parameters(gq))) == apply(pp, initial_parameters(pp))
        # named values: a key's scales by phase, its aim through the pointing
        q = quads[3].quad
        x = parameters(pp, q => (z = 0.2, power_scale = (step = 1.3,)))
        @test x[nq+48+3] == 1.3 && x[nq+3] == 1.0 && count(!=(0.0), x[1:nq]) == 1
        @test parameters(pp, q => (power_scale = 0.9,))[[nq + 3, nq + 51]] == [0.9, 0.9]
        @test_throws ArgumentError parameters(pp, q => (power_scale = (peak = 1.1,),))   # a fixed phase
        @test_throws ArgumentError Phase(:a, tP; power = 1.0, pulses, time = 1.0)
        @test_throws ArgumentError Phase(:a, tP; pulses)                                  # no time
        @test_throws ArgumentError Phase(:a, sP; observe = :capsule)                      # a ready scenario
        @test_throws ArgumentError PhasedScenario(gq, Phase(:a, sP), Phase(:a, sP))
        @test_throws ArgumentError evaluate(ps, PerPhase(:picket => bs[1], :peak => bs[3], :step => bs[2]))
        @test occursin("picket", sprint(show, MIME("text/plain"), bs))
    end

    @testset "one phase is the plain scenario" begin
        # bit for bit: the evaluation and a constrained run of the same objective
        sG = remake(sP; parameterization = gg)
        ps1 = PhasedScenario(gg, Phase(:a, sG; scaled = false))
        x = [0.25, -0.1, -0.2, 0.1, 0.15, 0.0, -0.1, 0.05]
        e0, e1 = evaluate(sG, x), evaluate(ps1, x)[:a]
        @test e1.deposition.element == e0.deposition.element && flux_ylm(e1) == flux_ylm(e0) && e1.flux == e0.flux
        ref = evaluate(sG, initial_parameters(gg))
        gate = mode_gate(ref, fill(0.02, 4))
        obj = drive_objective(; drive = mean_flux(ref))
        lo, hi = bounds(gg; z = 0.5, phi = 0.5)
        kw = (; x0 = x, lower = lo, upper = hi, evaluations = 25)
        r0 = optimize_pointing(sG, obj; constraints = (gate,), kw...)
        r1 = optimize_pointing(ps1, in_phase(:a, obj); constraints = (in_phase(:a, gate),), kw...)
        @test r1.x == r0.x && r1.value == r0.value && r1.evaluations == r0.evaluations
        @test r1.beams[:a] == r0.beams
        @test [c.value for c in r1.constraints] == [c.value for c in r0.constraints]
        @test [c.name for c in r1.constraints] == ["a." * c.name for c in r0.constraints]
    end

    @testset "against the duplicate-beam construction" begin
        # the scaled phases' beams once each, under the same quad names (so one aim per
        # quad), scaled by (quad, phase) through a `by` that counts its calls
        early = (:picket, :step)
        at(k) = k === :picket ? [remake(b; power = picket_power(b)) for b in quads] : with_pulses(quads, pulses, 2.5)
        copies = vcat(at(:picket), at(:step))
        gc = GroupPointing(copies; by = :quad, fields = (:z, :phi), relative = true, mirror = false)
        i = Ref(0)
        pc = PowerScales(gc; by = b -> (i[] += 1; (b.quad, early[(i[] - 1) ÷ 48 + 1])))
        @test pc.keys == pp.keys
        wb = apply(pc, y)
        ev = evaluate(ps, y)
        slot(k) = (findfirst(==(k), early) - 1) * 48 .+ (1:48)
        work = Dict(k => evaluate(ps[k], wb[slot(k)]) for k in early)
        work[:peak] = evaluate(ps[:peak], [remake(b; power = 2.0) for b in wb[slot(:picket)]])
        @test all(k -> rel(ev[k].deposition.element, work[k].deposition.element) < 1e-12, keys(ps))
        @test all(k -> rel(flux_ylm(ev[k]), flux_ylm(work[k])) < 1e-12, keys(ps))
        # the gradient of a mode departure through both, at the optimizer's tag and chunks
        cref = flux_ylm(evaluate(ps, initial_parameters(pp))[:step])
        f(e) = sum(abs2, mode_departure(e, cref))
        g1 = Briareus._scenario_gradient(ps, in_phase(:step, f), y)
        g0 = ForwardDiff.gradient(z -> f(evaluate(ps[:step], apply(pc, z)[slot(:step)])), y)
        @test rel(g1, g0) < 1e-10
    end

    @testset "gradient against finite differences" begin
        e0 = evaluate(ps, initial_parameters(pp))
        cref, D = flux_ylm(e0[:picket]), mean_flux(e0[:step])
        drive = in_phase(:peak, drive_objective(; drive = mean_flux(e0[:peak])))
        J(ev) = sum(abs2, mode_departure(ev[:picket], cref)) + (mean_flux(ev[:step]) / D - 1)^2 + drive(ev)
        G = Briareus._scenario_gradient(ps, J, y)
        idx = [1, 2, 17, nq - 1, nq + 1, nq + 30, nq + 49, nq + 90]   # aims and both phases' scales
        h = 1e-6
        fd = map(idx) do i
            e = zeros(length(y)); e[i] = h
            (J(evaluate(ps, y .+ e)) - J(evaluate(ps, y .- e))) / 2h
        end
        @test maximum(abs, G[idx] .- fd) < 1e-6 * maximum(abs, fd)
    end

    @testset "the laser-limit box" begin
        lo, hi = bounds(pp, y; scale = Inf, scale_min = 0.2, power_max = 2.0, z = 0.5, phi = 0.5)
        @test lo[nq+1:end] == fill(0.2, 96) && lo[1:nq] ≈ y[1:nq] .- 0.5
        bs = apply(pp, hi)
        # every beam of every phase at or below the limit, each scale's beam at it
        @test all(k -> all(b -> b.power <= 2.0 * (1 + 1e-12), bs[k]), keys(ps))
        @test all(b -> isapprox(b.power, 2.0; rtol = 1e-12), bs[:picket]) && all(b -> isapprox(b.power, 2.0; rtol = 1e-12), bs[:step])
        # a limit per beam; a relative box inside it
        lim(b) = inner(b) ? 3.0 : 2.0
        lo2, hi2 = bounds(pp, y; scale = 1.0, power_max = lim)
        bs2 = apply(pp, hi2)
        @test all(k -> all(b -> b.power <= lim(b) * (1 + 1e-12), bs2[k]), keys(ps))
        @test hi2[nq+1] ≈ 2 * y[nq+1] rtol = 1e-14                          # an inner picket quad: below its cap
        @test hi2[nq+48+9] ≈ 2.0 / with_pulses(quads, pulses, 2.5)[9].power rtol = 1e-14   # an outer step quad: capped
        # a fixed phase above the limit, and a cap below the floor
        @test_throws ArgumentError bounds(pp; power_max = 1.5)
        @test_throws ArgumentError bounds(pp; scale = Inf, scale_min = 6.0, power_max = 2.0)
    end

    # by cone: a few parameters, for the runs
    psg = PhasedScenario(gg, phases(); by = :group, observe = :capsule, deposition_kw = dk)
    @testset "constrained run meets every phase's gates" begin
        e0 = evaluate(psg, initial_parameters(psg.param))
        x0 = vcat([0.25, -0.1, -0.2, 0.1, 0.15, 0.0, -0.1, 0.05], ones(8))
        ex = evaluate(psg, x0)
        dep(e, k) = mode_departure(e[k], flux_ylm(e0[k]))
        # gates the start violates and the reference pointing meets
        gates = [in_phase(k, mode_gate(e0[k], max.(0.4 .* dep(ex, k), 1e-3))) for k in keys(psg)]
        floors = [in_phase(k, drive_floor(e0[k]; fraction = 1.0)) for k in (:picket, :step)]
        clear = in_phase(:peak, clearance_floor(ApertureKeepOut(tK); minimum = -0.8))   # below the start's
        @test any(r -> r.status === :violated, check_constraints(ex, gates))
        lo, hi = bounds(psg.param, x0; scale = Inf, scale_min = 0.5, power_max = 2.0, z = 0.5, phi = 0.5)
        r = optimize_pointing(psg, in_phase(:peak, drive_objective(; drive = mean_flux(e0[:peak])));
                              x0, lower = lo, upper = hi, constraints = [gates; floors; clear], evaluations = 1500)
        @test r.outcome in (:interior, :box)
        @test all(c -> c.status !== :violated, r.constraints)
        @test any(c -> startswith(c.name, "step.mode_gate[l = "), r.constraints)
        @test r.beams isa PerPhase && all(b -> b.power <= 2.0 * (1 + 1e-12), r.beams[:step])
        # the result's report per phase
        sm = pointing_summary(psg, r; drive_reference = e0, keepout = ApertureKeepOut(tK))
        @test keys(sm) == keys(psg) && all(s -> s isa PointingSummary, sm)
        @test all(k -> sm[k].departure ≈ mode_departure(evaluate(psg, r.x)[k], flux_ylm(e0[k])), (:picket, :step))
        @test occursin("step", sprint(show, MIME("text/plain"), sm))
        # staging: by cone → by quad (scales per quad), the same beams; then tied scales
        psq = remake(psg; pointing = gq, by = :quad)
        xq = parameters(psq.param, r.beams)
        bq = apply(psq.param, xq)
        @test all(k -> rel([b.aim[j] for b in bq[k] for j in 1:3], [b.aim[j] for b in r.beams[k] for j in 1:3]) < 1e-12 &&
                       rel([b.power for b in bq[k]], [b.power for b in r.beams[k]]) < 1e-12, keys(psg))
        # the finer parameterization starts where the coarse run ended
        cs = [gates; floors; clear]
        @test all(i -> isapprox(check_constraints(evaluate(psq, xq), cs)[i].value, r.constraints[i].value; rtol = 1e-9),
                  eachindex(r.constraints))
        # a restart under a time limit
        r2 = optimize_pointing(psg, in_phase(:peak, drive_objective(; drive = mean_flux(e0[:peak]))); from = r,
                               constraints = cs, evaluations = 10, time_limit = 600)
        @test r2.evaluations <= 10 && r2.beams isa PerPhase
        cone(b) = inner(b) ? :inner : :outer
        pst = remake(psq; by = cone)
        @test length(pst.param.keys) == 4 && pst.param.keys[1] == (:inner, :picket)
        xt = parameters(pst.param, r.beams)
        sc(k, c) = [r.beams[k][i].power / psg.param.power[k === :picket ? 1 : 2][i] for i in eachindex(quads) if cone(quads[i]) === c]
        @test xt[end-3:end] ≈ [sum(sc(k, c)) / length(sc(k, c)) for k in (:picket, :step) for c in (:inner, :outer)] rtol = 1e-12
        @test_throws ArgumentError remake(psg; scaled = (:late,))
        @test length(remake(psg; scaled = (:step,)).param.keys) == 4
        # an unreachable drive floor: phase 1's verdict at the laser's limit, and the floor
        # relaxed alone
        high = in_phase(:step, drive_floor(e0[:step]; fraction = 3.0))
        rv = optimize_pointing(psg, in_phase(:peak, drive_objective(; drive = mean_flux(e0[:peak]))); x0,
                               lower = lo, upper = hi, constraints = (high,), relaxation = true, evaluations = 300)
        @test rv.outcome === :infeasible_box && rv.phase1.verdict
        @test only(rv.relaxations).name == "step.drive" && only(rv.relaxations).status in (:relaxed, :stopped)
    end

    @testset "per-phase quantities and subsets" begin
        ev = evaluate(ps, y)
        e0 = evaluate(ps, initial_parameters(pp))
        ceil = in_phase(:peak, irradiance_ceiling(e0[:peak], :wall))
        light = in_phase(:picket, light_ceiling(e0[:picket], :cap_top))
        rows = check_constraints(ev, (ceil, light))
        @test rows[1].name == "peak.irradiance" && rows[1].value == peak_irradiance(ev[:peak], :wall)
        @test rows[2].name == "picket.light[cap_top]" && rows[2].value ≈ light_fraction(ev[:picket], :cap_top)
        @test map(mean_flux, ev) == PerPhase(:picket => mean_flux(ev[:picket]), :step => mean_flux(ev[:step]),
                                             :peak => mean_flux(ev[:peak]))
        # a subset: every phase's scenario restricted, the powers and the groups kept
        keep = findall(b -> b.group !== :outer50, quads)
        sub = subset(ps, keep)
        @test sub.param.n == 2 * length(keep) && length(sub.param.keys) == 2 * length(keep)
        yk = parameters(sub.param, PerPhase(k => apply(pp, y)[k][keep] for k in keys(ps)))
        es = evaluate(sub, yk)
        @test all(k -> rel(es[k].deposition.element, evaluate(subset(ps[k], keep; parameterization = sub.param.pointing),
                                                              apply(pp, y)[k][keep]).deposition.element) < 1e-12, keys(ps))
        # phases as remakes of one shared scenario, each with its own albedos and time
        # (as `history_state` supplies them): the operator, view factors and target shared
        shared = [Phase(:picket, sP; power = picket_power),
                  Phase(:step, remake(sP; albedo = albedos(tS, 2.5), time = 2.5); pulses, time = 2.5),
                  Phase(:peak, remake(sP; albedo = albedos(tK, 0.0)); power = 2.0, scaled = false)]
        pss = PhasedScenario(gq, shared; by = :quad)
        @test pss[:step].target === tP && pss[:step].deposition === sP.deposition
        esh = evaluate(pss, y)
        @test all(k -> rel(flux_ylm(esh[k]), flux_ylm(ev[k])) < 1e-10, keys(ps))
    end
end
