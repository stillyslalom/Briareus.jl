# Constraints in `optimize_pointing` (an augmented Lagrangian around L-BFGS-B), the
# wall-irradiance ceiling and the per-degree mode gate.
#
# References: the closed-form solution and multiplier of an equality-active quadratic
# program (the drive's coefficients are linear in the quad powers, `influence_matrix`);
# for the reported multipliers and bound sensitivities, the Karush–Kuhn–Tucker solution
# of a drive-like objective linear in the powers under two active linear ceilings, and
# restarts with a bound moved;
# for phase 1's relaxations and the drive floor, fills of the quads in decreasing drive
# per power (the largest drive at a total power, the least power for a drive);
# central finite differences for the gradients; the power-mean inequality and its limit
# for the irradiance norm; synthetic coefficient vectors for the mode departure.

@testset "constraints" begin
    geom = (; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)

    @testset "quantities" begin
        @test_throws ArgumentError Constraint(:a, identity)
        @test_throws ArgumentError Constraint(:a, identity; upper = 1.0, lower = 0.0)
        @test_throws ArgumentError Constraint(:a, identity; upper = 0.0)        # a zero limit needs a scale
        @test_throws ArgumentError Constraint(:a, identity; upper = Inf)
        @test_throws DimensionMismatch Constraint(:a, identity; upper = [1.0, 2.0], scale = [1.0, 2.0, 3.0])
        c = Constraint(:a, identity; lower = [1.0, 2.0], labels = ("x", "y"))
        @test !c.upper && c.scale == [1.0, 2.0] && c.labels == ["x", "y"]

        # One coefficient off by δ c₀₀ gives a departure of √(2l+1) |δ| at that degree and
        # nowhere else, whatever the overall levels of the two shapes.
        r = [1.0; 0.01 .* sin.(1:24)]
        x = 3.0 .* r
        x[mode_index(2, 1)] += 3.0 * 0.004
        d = mode_departure(x, r; lmax = 4)
        @test d[2] ≈ sqrt(5) * 0.004 rtol = 1e-12
        @test all(abs.(d[[1, 3, 4]]) .< 1e-15)
        @test mode_departure(r, r; lmax = 4) == zeros(4)
        # An axisymmetric difference is |aₗ/a₀ - aₗ'/a₀'|, with c_l0 = aₗ √(4π/(2l+1)).
        y = copy(r)
        y[mode_index(3, 0)] += 0.02 * sqrt(4π / 7) / sqrt(4π) * r[1]
        @test mode_departure(y, r; lmax = 3)[3] ≈ 0.02 rtol = 1e-12
        @test_throws ArgumentError mode_departure(r[1:9], r; lmax = 4)
        @test_throws ArgumentError mode_gate(r, [0.01, 0.0])
        @test mode_gate(r, (; rms_mode = [0.001, 0.002])).limit ≈ [sqrt(3) * 0.001, sqrt(5) * 0.002]
    end

    # A quadratic program with a closed form fits the quad powers to 1.1 times the
    # coefficients of the 48-quad drive through l = 2. The coefficients are linear in the
    # powers through the influence matrix M. The fit is regularized toward the start, under
    # a ceiling on the total power that the fit would exceed. The KKT conditions are H P + ν
    # 1 = r and 1ᵀP = cap.
    @testset "augmented Lagrangian against a closed form" begin
        t = hohlraum(; geom..., h = 0.8, capsule_h = 0.9)
        quads = nif_quads(power = 8.0)
        s = Scenario(t, GroupPointing(quads; by = :quad, fields = (:power,)); observe = :capsule,
                     viewfactor_kw = (; closed = true))
        P0 = initial_parameters(s.param)
        M = influence_matrix(s, quads; lmax = 2)
        b = 1.1 .* (M * P0)
        ε = 1e-3
        f(ev) = sum(abs2, flux_ylm(ev; lmax = 2) .- b) / b[1]^2 +
                ε * sum(abs2, (q.power - 8.0) / 8.0 for q in ev.beams)
        cap = 0.95 * sum(P0)
        # The ceiling is on a quantity that depends only on the beams, and the floor stays
        # inactive.
        total = Constraint(:total_power, ev -> sum(q.power for q in ev.beams); upper = cap)
        drive = Constraint(:drive, mean_flux; lower = 0.5 * mean_flux(evaluate(s, P0)))
        H = 2 .* (M' * M) ./ b[1]^2 + (2ε / 64) * I
        rhs = 2 .* (M' * b) ./ b[1]^2 + (2ε / 64) .* P0
        o = ones(length(P0))
        ν = (o' * (H \ rhs) - cap) / (o' * (H \ o))
        Pstar = H \ (rhs - ν .* o)
        lo, hi = zeros(48), fill(20.0, 48)
        res = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = (total, drive),
                                feasibility = 1e-8)
        @test res.outcome === :interior && res.reason === :converged && res.outer > 1
        @test norm(res.x - Pstar) / norm(Pstar) < 1e-4
        # the multiplier of the normalized constraint (ΣP - cap)/cap is ν cap
        @test res.constraints[1].multiplier ≈ ν * cap rtol = 1e-3
        @test [r.status for r in res.constraints] == [:active, :satisfied]
        @test res.constraints[2].multiplier == 0 && res.constraints[2].violation < -1
        @test res.value ≈ f(evaluate(s, res.x)) && res.initial_value ≈ f(evaluate(s, P0))
        @test occursin("constraints met inside the box", sprint(show, MIME"text/plain"(), res))
        @test check_constraints(evaluate(s, P0), (total, drive))[1].status === :violated
        # The start is already feasible and on the cap, and the objective is scaled so that
        # its gradient lies below the first inner tolerance (1e-3). The run still reaches
        # the solution; it does not end at the start with zero multipliers on the first
        # inner run.
        small(ev) = 1e-4 * f(ev)
        Pa = P0 .* (cap / sum(P0))
        on_cap = optimize_pointing(s, small; x0 = Pa, lower = lo, upper = hi, constraints = (total,))
        @test on_cap.outcome === :interior && norm(on_cap.x - Pstar) / norm(Pstar) < 1e-3
        @test norm(Pa - Pstar) / norm(Pstar) > 0.05

        # Warm restarts. Each inner run starts its L-BFGS-B afresh from the state stored in
        # the result, so a run split at an outer iteration and continued is identical to the
        # uninterrupted run. A run split inside an inner run by the evaluation budget loses
        # the L-BFGS-B memory of that inner run and takes another path to the same strictly
        # convex optimum. The designs agree within the convergence tolerance (feasibility
        # 1e-8, with the run 1e-8 from the KKT point), tested at 1e-6.
        cs2 = (total, drive)
        a = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = cs2,
                              feasibility = 1e-8, outer = 2)
        @test a.outcome === :iterations && a.outer == 2 && length(a.state.λ) == length(a.constraints)
        ab = optimize_pointing(s, f; from = a, constraints = cs2, feasibility = 1e-8)
        @test ab.x == res.x && ab.outer == res.outer - 2 && a.evaluations + ab.evaluations == res.evaluations
        @test ab.lower == lo && ab.upper == hi
        a2 = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = cs2,
                               feasibility = 1e-8, evaluations = 20)
        @test a2.outcome === :evaluations && a2.state.μ > 0
        b2 = optimize_pointing(s, f; from = a2, constraints = cs2, feasibility = 1e-8)
        @test b2.outcome === :interior && norm(b2.x - res.x) / norm(res.x) < 1e-6
        @test b2.constraints[1].multiplier ≈ res.constraints[1].multiplier rtol = 1e-3
        @test_throws DimensionMismatch optimize_pointing(s, f; from = a, constraints = (total,))
        # a time limit ends a run at its first evaluation past it, with or without
        # constraints, keeping the best design (and the state to continue from)
        for c in ((), cs2)
            tl = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = c, time_limit = 1e-9)
            @test tl.reason === tl.outcome === :time_limit && tl.stopped && tl.evaluations == 1
            @test occursin("time limit", sprint(show, MIME"text/plain"(), tl))
        end
        @test_throws ArgumentError optimize_pointing(s, f; x0 = P0, time_limit = 0)
        # without constraints, `from` sets the start and the box only
        @test optimize_pointing(s, f; from = a, iterations = 1).lower == lo

        # the box keeps every power above the cap's share: infeasible at the bounds
        boxed = optimize_pointing(s, f; x0 = P0, lower = 0.99 .* P0, upper = hi, constraints = (total,))
        @test boxed.outcome === :infeasible_box && boxed.reason === :converged
        @test boxed.constraints[1].status === :violated && length(boxed.at_bounds) == 48
        @test occursin("ended with constraints violated", sprint(show, MIME"text/plain"(), boxed))
        # With phase 1 the verdict is the same. Only the cap is relaxed, since there is no
        # other constraint, to the least total power the box allows, with every power at its
        # lower bound.
        bv = optimize_pointing(s, f; x0 = P0, lower = 0.99 .* P0, upper = hi, constraints = (total,),
                               relaxation = true)
        @test bv.outcome === :infeasible_box && only(bv.relaxations).status === :relaxed
        @test only(bv.relaxations).violation ≈ 0.99 * sum(P0) / (0.95 * sum(P0)) - 1 rtol = 1e-9

        # Phase 1 and the relaxations are checked on a conflict with closed forms. The drive
        # (`mean_flux`) is linear in the powers, with `m` per quad at unit power. The
        # largest drive at a given total power and the least total power for a given drive
        # are therefore fills of the quads in decreasing `m` up to the 20 TW of the box. The
        # cap is 120 TW and the drive floor is 5% above the largest drive under the cap.
        # When only the floor is relaxed, it falls to that drive, by 1 - 1/1.05 of its
        # value. When only the cap is relaxed, it rises to the least power for the floor.
        m = [mean_flux(evaluate(s, Float64.(1:48 .== i))) for i in 1:48]
        @test dot(m, P0) ≈ mean_flux(evaluate(s, P0)) rtol = 1e-12
        order = sortperm(m; rev = true)
        function drive_at(P)
            D = 0.0
            for i in order
                p = min(20.0, P); D += m[i] * p; P -= p
            end
            return D
        end
        function power_for(D)
            P = 0.0
            for i in order
                m[i] > 0 || continue
                p = clamp(D / m[i], 0.0, 20.0); P += p; D -= m[i] * p
            end
            return P
        end
        cap120 = Constraint(:total_power, ev -> sum(q.power for q in ev.beams); upper = 120.0)
        Dfl = 1.05 * drive_at(120.0)
        above = Constraint(:drive, mean_flux; lower = Dfl)
        @test power_for(Dfl) < 20 * count(>(0), m)       # the floor by itself is feasible
        # For comparison, the augmented Lagrangian runs without phase 1.
        today = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = (cap120, above),
                                  feasibility = 1e-8, phase1 = false)
        v = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = (cap120, above),
                              feasibility = 1e-8, relaxation = true)
        @test today.outcome in (:infeasible, :infeasible_box) && today.phase1 === nothing && isempty(today.relaxations)
        @test v.outcome in (:infeasible, :infeasible_box) && v.reason === :converged && v.outer == 0
        @test v.phase1.verdict && !v.phase1.feasible && v.phase1.evaluations < today.evaluations
        @test all(r -> r.status === :violated, v.constraints)
        r1, r2 = v.relaxations
        @test (r1.name, r2.name) == ("total_power", "drive") && r1.status === r2.status === :relaxed
        @test r1.violation ≈ power_for(Dfl) / 120 - 1 rtol = 1e-5
        @test r1.relaxed_limit ≈ [power_for(Dfl)] rtol = 1e-7
        @test r2.violation ≈ 1 - 1 / 1.05 rtol = 1e-5
        @test r2.relaxed_limit ≈ [drive_at(120.0)] rtol = 1e-7
        # each relaxation's design meets the other constraint
        @test r1.constraints[2].status !== :violated && r2.constraints[1].status !== :violated
        @test v.evaluations == v.phase1.evaluations + r1.evaluations + r2.evaluations
        shown = sprint(show, MIME"text/plain"(), v)
        @test occursin("no design meets every constraint", shown) && occursin("relaxed alone", shown)
        # Phase 1 has a budget separate from the inner runs' `iterations`: with inner runs
        # capped at 3 iterations it still reaches the verdict, at the same violations
        # (before, it stopped at 3 iterations without one and the main run began at the
        # start). Stopped by the run's evaluation budget first, it says so.
        v3 = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = (cap120, above),
                               feasibility = 1e-8, iterations = 3, evaluations = 3000)
        @test v3.phase1.verdict && v3.outer == 0 && v3.phase1.reason in (:converged, :stalled)
        @test [r.violation for r in v3.constraints] ≈ [r.violation for r in v.constraints] rtol = 1e-4
        v5 = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = (cap120, above),
                               feasibility = 1e-8, evaluations = 5)
        @test !v5.phase1.verdict && v5.phase1.reason === :evaluations && v5.outcome === :evaluations
        @test occursin("the evaluation budget ended it first", sprint(show, MIME"text/plain"(), v5))
        # no verdict on a feasible run from far outside the cap: the KKT point as without
        # phase 1, its main run started from phase 1's first design within the limits
        far = optimize_pointing(s, f; x0 = clamp.(2.5 .* P0, lo, hi), lower = lo, upper = hi,
                                constraints = cs2, feasibility = 1e-8)
        @test far.outcome === :interior && far.phase1.feasible && !far.phase1.verdict
        @test norm(far.x - Pstar) / norm(Pstar) < 1e-4
        # phase 1 is the default (`res` started over the cap); `phase1 = false` skips it
        @test res.phase1.feasible && !res.phase1.verdict && isempty(res.relaxations)
        @test today.phase1 === nothing

        # The iterate of least violation. Without phase 1 the conflicting run ends at the
        # penalty's ceiling: it returns its iterate of least violation (the largest
        # normalized violation of its rows, against `check_constraints` at its last
        # iterate, `state.x`), and no iterate met both limits.
        lastrows = check_constraints(evaluate(s, today.state.x), (cap120, above); feasibility = 1e-8)
        @test today.least_violation.x == today.x && today.last_feasible === nothing
        @test today.least_violation.violation == maximum(r.violation for r in today.constraints)
        @test today.least_violation.violation <= maximum(r.violation for r in lastrows)
        @test occursin("no iterate met every limit", sprint(show, MIME"text/plain"(), today))
        # A run stopped by its budget above a limit that its start met: the start on the
        # cap, the objective pulling the powers above it. The last iterate (`state.x`)
        # exceeds the cap; the result is the iterate of least violation, here one within
        # the cap, and a restart continues from the last iterate.
        Pa = P0 .* (cap / sum(P0))
        cut = optimize_pointing(s, f; x0 = Pa, lower = lo, upper = hi, constraints = (total,),
                                evaluations = 8, phase1 = false)
        @test cut.outcome === :evaluations && cut.least_violation.x == cut.x
        @test only(check_constraints(evaluate(s, cut.state.x), total)).status === :violated
        @test only(cut.constraints).status !== :violated && cut.last_feasible !== nothing
        @test occursin("returned the iterate of least violation", sprint(show, MIME"text/plain"(), cut))
        @test Briareus._restart_x(cut) == cut.state.x

        # Inner runs ended by their cap above the limits hold the penalty (the result's
        # `capped` counts the capped runs): with two iterations per inner run the feasible
        # quadratic program still converges to the KKT point, where escalating the penalty
        # after every capped run ended it at the penalty's ceiling with a false verdict.
        two = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = cs2, feasibility = 1e-8,
                                iterations = 2, evaluations = 3000, outer = 100, phase1 = false)
        @test two.outcome === :interior && two.capped > 0 && norm(two.x - Pstar) / norm(Pstar) < 1e-3
        @test res.capped == 0 && res.last_feasible.x == res.x == res.state.x

        # a drive floor active at the solution: the least total power that keeps 0.9 of the
        # start's drive, a fill of the quads in decreasing `m`
        ref0 = evaluate(s, P0)
        df = drive_floor(ref0; fraction = 0.9)
        @test !df.upper && only(df.limit) == 0.9 * mean_flux(ref0) && df.name === :drive
        lp = optimize_pointing(s, ev -> sum(q.power for q in ev.beams) / 384; x0 = P0, lower = lo, upper = hi,
                               constraints = (df,), feasibility = 1e-8)
        @test lp.outcome === :box && only(lp.constraints).status === :active && only(lp.constraints).multiplier > 0
        @test sum(lp.x) ≈ power_for(0.9 * mean_flux(ref0)) rtol = 1e-6
        @test_throws ArgumentError drive_floor(ref0; fraction = 0)

        # the budget ends a run first
        short = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, evaluations = 10, constraints = [total])
        @test short.outcome === :evaluations && short.stopped && short.evaluations == 10
        # Without constraints the run is the plain run, and only the box sets its outcome.
        plain = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi)
        @test plain.outcome === :interior && isempty(plain.constraints) && plain.outer == 0
        @test plain.least_violation === plain.last_feasible === nothing && plain.capped == 0
        @test sum(plain.x) > cap
        @test isempty(plain.bound_sensitivity) && isempty(plain.bound_gradient) && isnan(plain.stationarity)

        # bound sensitivities against central differences of the optimum in the bound (the
        # optimal value of a quadratic program with a fixed active set is quadratic in the
        # bound, so the difference is exact up to the runs' convergence): one power's upper
        # bound below its optimum, without and with the constraints
        i = argmax(plain.x)
        δ = 0.05
        for (c, xref, fz) in (((), plain.x, 1e-12), (cs2, res.x, 1e-8))
            u = xref[i] - 0.5
            opt(u) = optimize_pointing(s, f; x0 = P0, lower = lo, upper = [k == i ? u : 20.0 for k in 1:48],
                                       constraints = c, feasibility = fz)
            r0 = opt(u)
            j = findfirst(==(r0.names[i]), r0.at_bounds)
            @test r0.outcome === :box && r0.at_bounds == [r0.names[i]]
            fd = (opt(u + δ).value - opt(u - δ).value) / 2δ
            @test r0.bound_sensitivity[j] < 0
            @test r0.bound_sensitivity[j] ≈ fd rtol = 1e-4
            # the derivative of the objective, against its central difference at the result
            fx(h) = f(evaluate(s, r0.x .+ h .* (1:48 .== i)))
            @test r0.bound_gradient[j] ≈ (fx(1e-4) - fx(-1e-4)) / 2e-4 rtol = 1e-6
            @test isempty(c) ? r0.bound_gradient == r0.bound_sensitivity :
                  abs(r0.bound_gradient[j] - r0.bound_sensitivity[j]) > 1e-3 * abs(fd)
            @test occursin("d objective/d bound", sprint(show, MIME"text/plain"(), r0))
        end

        # a ceiling's reference and the optimized scenario share the laser–plasma model
        # and its time: another model, or none against one, is an error; a model built
        # alike is the same model
        bs(f) = Backscatter(:inner23 => f, :inner30 => ([0.0, 5.0], [0.0, 0.1]))
        ceil0 = irradiance_ceiling(evaluate(s, P0), :wall)
        sb = remake(s; lpi = bs(0.1))
        @test_throws ArgumentError ceil0.f(evaluate(sb, P0))
        @test_throws ArgumentError check_constraints(evaluate(sb, P0), ceil0)
        ceilb = irradiance_ceiling(evaluate(sb, P0), :wall)
        @test_throws ArgumentError ceilb.f(evaluate(s, P0))
        @test_throws ArgumentError ceilb.f(evaluate(remake(s; lpi = bs(0.2)), P0))
        @test ceilb.f(evaluate(remake(s; lpi = bs(0.1)), P0)) == ceilb.f(evaluate(sb, P0))
        @test_throws ArgumentError ceilb.f(evaluate(remake(sb; time = 2.0), P0))
        # so do a mode gate and a drive floor from a reference evaluation; the gate's values
        # are those from the reference's coefficients, which record no model (unchecked)
        gate_ev = mode_gate(evaluate(sb, P0), [0.01, 0.01])
        gate_c = mode_gate(flux_ylm(evaluate(sb, P0); lmax = 2), [0.01, 0.01])
        xb = P0 .* (1 .+ 0.1 .* sin.(1:48))
        @test gate_ev.f(evaluate(sb, xb)) == gate_c.f(evaluate(sb, xb)) && gate_ev.limit == gate_c.limit
        @test gate_ev.f(evaluate(remake(s; lpi = bs(0.1)), xb)) == gate_ev.f(evaluate(sb, xb))
        @test_throws ArgumentError gate_ev.f(evaluate(s, xb))
        @test_throws ArgumentError gate_ev.f(evaluate(remake(sb; time = 2.0), xb))
        @test length(gate_c.f(evaluate(s, xb))) == 2
        @test mode_gate(evaluate(sb, P0), (; rms_mode = [0.001, 0.002])).limit ≈ [sqrt(3) * 0.001, sqrt(5) * 0.002]
        dfb = drive_floor(evaluate(sb, P0); fraction = 0.9)
        @test dfb.f(evaluate(sb, xb)) == mean_flux(evaluate(sb, xb))
        @test_throws ArgumentError dfb.f(evaluate(s, xb))
    end

    # The multipliers from the converged active set (`PointingResult.constraints`) are
    # checked against the Karush–Kuhn–Tucker multipliers of a closed form. The objective is
    # drive-like, linear in the quad powers with a small regularization toward the start. It
    # is held under ceilings on the total power and on the power of the inner cones, both
    # active at the optimum, and above an inactive floor on the power of the 50° cone. The
    # KKT conditions are ε (P - P₀) = c - Aᵀν with A P = b on the active rows, and the
    # normalized constraint (aᵀP - b)/b has the multiplier ν b. The stationarity residual is
    # computed from the closed-form gradients; the optimizer gradients are not used.
    @testset "multipliers from the converged active set" begin
        t = hohlraum(; geom..., nz = 8, nphi = 16, capsule_ntheta = 6, capsule_nphi = 8)
        quads = nif_quads(power = 8.0)
        s = Scenario(t, GroupPointing(quads; by = :quad, fields = (:power,)))
        P0 = initial_parameters(s.param)
        n = length(P0)
        inner = [q.group in (:inner23, :inner30) for q in quads]
        o50 = [q.group === :outer50 for q in quads]
        c = (1 .+ 0.3 .* sin.(1:n)) .* (1 .+ 0.5 .* inner) ./ 384
        ε = 2e-4
        f(ev) = -sum(c[i] * q.power for (i, q) in enumerate(ev.beams)) +
                ε / 2 * sum(abs2, q.power - 8.0 for q in ev.beams)
        cap, capI = 1.01 * 384, 130.0
        cs = (Constraint(:total, ev -> sum(q.power for q in ev.beams); upper = cap),
              Constraint(:inner, ev -> sum(q.power for (i, q) in enumerate(ev.beams) if inner[i]); upper = capI),
              Constraint(:outer50, ev -> sum(q.power for (i, q) in enumerate(ev.beams) if o50[i]); lower = 50.0))
        function kkt(A, b)
            ν = (A * A') \ (A * c .+ ε .* (A * P0 .- b))
            return P0 .+ (c .- A' * ν) ./ ε, ν
        end
        A = [ones(n)'; Float64.(inner)']
        Pk, ν = kkt(A, [cap, capI])
        λk = [ν .* [cap, capI]; 0.0]
        lo, hi = zeros(n), fill(100.0, n)
        @test all(lo .< Pk .< hi) && sum(Pk[o50]) > 60 && all(>(0), ν)
        # |∇f + Σ λⱼ ∇gⱼ|∞ at P, in closed form
        resid(P, λ) = norm(-c .+ ε .* (P .- 8.0) .+ λ[1] / cap .+ λ[2] .* inner ./ capI .- λ[3] .* o50 ./ 50, Inf)
        mults(r) = [x.multiplier for x in r.constraints]
        r = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = cs)
        @test r.outcome === :interior && norm(r.x - Pk) / norm(Pk) < 1e-3
        @test all(isapprox.(mults(r)[1:2], λk[1:2]; rtol = 1e-3))
        # By complementarity, the multiplier of the inactive floor is zero.
        @test r.constraints[3].status === :satisfied && r.constraints[3].multiplier == 0
        @test resid(r.x, mults(r)) <= resid(r.x, r.state.λ)
        # the result's stationarity is that residual (every parameter inside the box)
        @test r.stationarity ≈ resid(r.x, mults(r)) rtol = 1e-6
        @test occursin("stationarity", sprint(show, MIME"text/plain"(), r))
        # the display notes a row `:satisfied` with a positive multiplier (an exact
        # ceiling's, whose multiplier is its stand-in's at the working limit), here
        # written into a copy of the result; the rows themselves are unchanged
        @test !occursin("stand-in", sprint(show, MIME"text/plain"(), r))
        x1 = r.constraints[1]
        odd = ConstraintRow(x1.name, x1.value, x1.limit, x1.upper, -0.01, 0.5, :satisfied)
        ro = PointingResult((k === :constraints ? [odd; r.constraints[2:end]] : getfield(r, k)
                             for k in fieldnames(PointingResult))...)
        shown = sprint(show, MIME"text/plain"(), ro)
        note = filter(l -> occursin("stand-in", l), split(shown, '\n'))
        @test length(note) == 1 && occursin("total", only(note)) && ro.constraints[1] === odd

        # The old failure started from that design at a large penalty, μ = 10⁶, the order a
        # drive run under gates reaches. The inner run ended on its value tolerance a few
        # 1e-6 inside the limit of the inner cones, and the first-order update λ + μ g
        # zeroed the multiplier of that limit. That multiplier, `state.λ`, is the value a
        # further restart continues from. Its stationarity residual was thousands of times
        # that of the reported multipliers, which still match the closed form. The penalty
        # is set by copying the result with another `state`.
        big = PointingResult((k === :state ? merge(r.state, (; μ = 1e6, μmax = 1e14)) : getfield(r, k)
                              for k in fieldnames(PointingResult))...)
        rb = optimize_pointing(s, f; from = big, constraints = cs)
        @test rb.outcome === :interior && norm(rb.x - Pk) / norm(Pk) < 1e-3
        @test all(isapprox.(mults(rb)[1:2], λk[1:2]; rtol = 1e-3)) && mults(rb)[3] == 0
        @test abs(rb.state.λ[2] - λk[2]) > 0.5 * λk[2]
        @test resid(rb.x, mults(rb)) < 1e-2 * resid(rb.x, rb.state.λ)

        # Bound sensitivities, from the reported multipliers: one quad's upper bound 1 TW
        # below its optimum. Closed form: the bound as a third active row, ∂f*/∂u = -ν_u.
        # Direct: restarts from the result with the bound moved by ±δ (`from`), a central
        # difference, exact for a quadratic program whose active set stays.
        i = argmax(Pk)
        u = Pk[i] - 1.0
        hiu = [k == i ? u : 100.0 for k in 1:n]
        Pu, νu = kkt([A; (1:n .== i)'], [cap, capI, u])
        rbnd = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hiu, constraints = cs)
        @test rbnd.outcome === :box && rbnd.at_bounds == [rbnd.names[i]] && norm(rbnd.x - Pu) / norm(Pu) < 1e-3
        @test all(isapprox.(mults(rbnd)[1:2], νu[1:2] .* [cap, capI]; rtol = 1e-3))
        @test only(rbnd.bound_sensitivity) ≈ -νu[3] rtol = 1e-3
        δ = 0.2
        moved(d) = optimize_pointing(s, f; from = rbnd, upper = [k == i ? u + d : 100.0 for k in 1:n],
                                     constraints = cs)
        rp, rm = moved(δ), moved(-δ)
        @test rp.at_bounds == rm.at_bounds == [rbnd.names[i]]
        @test only(rbnd.bound_sensitivity) ≈ (rp.value - rm.value) / 2δ rtol = 1e-2
        @test only(rbnd.bound_gradient) ≈ -c[i] + ε * (rbnd.x[i] - 8.0) rtol = 1e-6
        @test only(rbnd.bound_sensitivity) < 0

        # The inner value tolerance is `inner_reltol`. Its default given explicitly gives a
        # run identical to the default run. With the plain `c` the cap of 130 TW on the
        # inner cones binds with a small multiplier. At the default the run ends 0.8% from
        # the KKT point (research notes). Continued through `from` with the value tolerance
        # off, so that the inner runs go to `g_tol`, it reaches the KKT point and
        # multipliers. Any converged run of this strictly convex problem does, so no
        # rounding-dependent path is tested.
        @test optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = cs,
                                inner_reltol = 1e-10).x == r.x
        c0 = (1 .+ 0.3 .* sin.(1:n)) ./ 384
        f0(ev) = -sum(c0[i] * q.power for (i, q) in enumerate(ev.beams)) +
                 ε / 2 * sum(abs2, q.power - 8.0 for q in ev.beams)
        ν0 = (A * A') \ (A * c0 .+ ε .* (A * P0 .- [cap, capI]))
        P0k, λ0 = P0 .+ (c0 .- A' * ν0) ./ ε, ν0 .* [cap, capI]
        @test all(lo .< P0k .< hi) && all(>(0), ν0)
        r0 = optimize_pointing(s, f0; x0 = P0, lower = lo, upper = hi, constraints = cs)
        rt = optimize_pointing(s, f0; from = r0, constraints = cs, inner_reltol = 0)
        @test rt.outcome === :interior && norm(rt.x - P0k) / norm(P0k) < 1e-5
        @test all(isapprox.(mults(rt)[1:2], λ0; rtol = 1e-3)) && mults(rt)[3] == 0
        @test_throws ArgumentError optimize_pointing(s, f0; from = r0, constraints = cs, inner_reltol = -1.0)

        # An inner evaluation cap is set with `inner_evaluations`. A cap above every inner
        # run, larger than the evaluation count of the whole run, leaves the run unchanged.
        # At 20, below the longest inner runs of the default run (about 40), every inner run
        # stays within the cap, as the callback evaluations per outer iteration show, and
        # the run still converges to the KKT point.
        big_cap = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = cs,
                                    inner_evaluations = r.evaluations + 1)
        @test big_cap.x == r.x && big_cap.evaluations == r.evaluations && big_cap.outer == r.outer
        per = Dict{Int,Int}()
        count_inner = st -> (per[st.outer] = get(per, st.outer, 0) + 1; false)
        rcap = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = cs, inner_evaluations = 20,
                                 outer = 100, callback = count_inner)
        @test !isempty(per) && all(<=(20), values(per)) && sum(values(per)) == rcap.evaluations
        @test rcap.outcome === :interior && norm(rcap.x - Pk) / norm(Pk) < 1e-3
        @test_throws ArgumentError optimize_pointing(s, f; x0 = P0, constraints = cs, inner_evaluations = 0)
        @test_throws ArgumentError optimize_pointing(s, f; x0 = P0, inner_evaluations = 20)   # no constraints
    end

    # The dropped-bundle tutorial's design on a coarser wall mesh: bundle 11 freed, the
    # other 46 quads re-pointed within ±0.3 mm for the largest drive with each mode's
    # departure from the 48-quad shape within the spread of 50 µm pointing errors, under
    # the 48-quad pointing's wall irradiance.
    @testset "dropped bundle: gates under a ceiling in one call" begin
        t = hohlraum(; geom..., wall = gold(albedo = 0.85, conversion = 0.8), h = 0.45)
        aims = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20), :inner30 => (0.60, 3.63),
                    :inner23 => (1.09, 3.09))
        quads = map(nif_quads(power = 8.0)) do q
            r, z = aims[q.group]
            φ = port_angles(q)[2]
            remake(q; aim = (-r * cosd(φ), -r * sind(φ), sign(q.port[3]) * z))
        end
        s = Scenario(t, GroupPointing(quads; by = :quad, fields = (:power,)); observe = :capsule,
                     deposition_kw = (; resolution = 0.6, reach = 0.5, level = 0))
        ref = evaluate(s, quads)
        c48 = flux_ylm(ref; lmax = 4)
        keep = setdiff(eachindex(quads), nif_bundles(quads)["11"])
        rp = subset(s, keep; parameterization = GroupPointing(quads[keep]; by = :quad, fields = (:phi, :z)))
        x0 = initial_parameters(rp.param)
        lower, upper = bounds(rp.param, x0; phi = 0.3, z = 0.3)
        ceiling = irradiance_ceiling(ref, :wall)
        gates = mode_gate(c48, mode_uncertainty(s, quads; pointing = 0.05, lmax = 4))

        # The power means of the irradiance lie below the peak and rise with p, and p = 1
        # gives the area mean. They are homogeneous in the beam powers. A weight of one
        # reproduces the deposition, and a weight of zero on the inner cones reproduces the
        # deposition with those cones at zero power.
        pk = peak_irradiance(ref, :wall)
        wall = elements(t, :wall)
        norms = [irradiance_norm(ref, :wall; p) for p in (8, 32, 64)]
        @test issorted(norms) && norms[end] < pk
        @test irradiance_norm(ref, :wall; p = 1) ≈ sum(ref.deposition.element[wall]) / sum(t.elem_area[wall]) rtol = 1e-12
        @test peak_irradiance(ref, :wall; weight = q -> 2.0) ≈ 2pk rtol = 1e-12
        @test irradiance_norm(ref, :wall; weight = q -> 1.0) == norms[end]
        outer = q -> q.group in (:outer44, :outer50) ? 1.0 : 0.0
        outer_only = [remake(q; power = outer(q) * q.power) for q in quads]
        @test peak_irradiance(ref, :wall; weight = outer) ≈ peak_irradiance(evaluate(s, outer_only), :wall) rtol = 1e-12
        @test ceiling.f(ref) ≈ pk rtol = 1e-12
        @test only(ceiling.limit) == pk
        @test_throws ArgumentError peak_irradiance(ref, :no_such_part)

        # gradients of the constraint values against central differences
        d = sin.(1:length(x0))
        along(f) = τ -> f(evaluate(rp, x0 .+ τ .* d))
        for f in (ceiling.f, ev -> sum(gates.f(ev) .* (1:4)))
            ad = ForwardDiff.derivative(along(f), 0.0)
            fd = (along(f)(1e-6) - along(f)(-1e-6)) / 2e-6
            @test ad ≈ fd rtol = 1e-5
        end

        # At the start every gate is violated and the ceiling is at its limit, because the
        # freed quads do not reach the hottest elements.
        start = check_constraints(evaluate(rp, x0), (ceiling, gates))
        @test count(r -> r.status === :violated, start) == 4 && start[1].status !== :violated
        @test peak_irradiance(evaluate(rp, x0), :wall) ≈ pk rtol = 1e-12

        # With the exact peak held (`exact = true`, the default), the working limit of the
        # stand-in is lowered until the peak over the elements is within `feasibility` of
        # the limit. The stand-in by itself (`exact = false`) does not bound the peak to a
        # fixed margin. Designs with the same gates differ in the peak by several percent.
        # From this start, the peaks are 0.997 and 1.030 of the limit, with the norm at
        # 0.974 and at the limit, under two roundings (research notes). Under
        # `--check-bounds=yes` the run reaches another design of the same flat family, so
        # the test checks only what every admissible design meets.
        plain = irradiance_ceiling(ref, :wall; exact = false)
        @test ceiling.f.exact && !plain.f.exact
        res = optimize_pointing(rp, drive_objective(; drive = mean_flux(ref)); x0, lower, upper,
                                constraints = (ceiling, gates), iterations = 1000, evaluations = 3000)
        ev = evaluate(rp, res.x)
        @test res.outcome in (:interior, :box)
        @test all(r -> r.status !== :violated, res.constraints)
        @test all(mode_departure(ev, c48) .<= gates.limit .* (1 + 1e-3))
        @test peak_irradiance(ev, :wall) <= pk * (1 + 1e-3)
        @test res.value <= res.initial_value + 1e-3
        # complementarity of the reported multipliers: zero for every gate inside its
        # limit by more than `feasibility`. The exact ceiling's multiplier is its
        # stand-in's at the working limit, which can be active while the peak the row
        # reports ends inside the limit (1.2% inside under one CI rounding): only its
        # sign is tested.
        @test all(r -> r.multiplier >= 0, res.constraints)
        @test all(r -> r.status !== :satisfied || r.multiplier == 0, res.constraints[2:end])
        @test res.stationarity >= 0
        # In the rows, the exact ceiling reports the peak against the limit, and the plain
        # ceiling reports its norm.
        @test res.constraints[1].value == peak_irradiance(ev, :wall) && res.constraints[1].limit == pk
        @test check_constraints(ev, ceiling)[1].value == peak_irradiance(ev, :wall)
        @test check_constraints(ev, plain)[1].value == plain.f(ev) == ceiling.f(ev)
    end

    # Keep-outs and the light on named parts as constraints, on a coarse hohlraum. The
    # references: the quantities that define them; for the constrained runs, the root of
    # the one free parameter's constraint value found by bisection (the objective pulls
    # the 50° cone's aim toward the center, monotonically, so the solution is where the
    # constraint becomes active); central differences for the gradients.
    @testset "keep-outs and the light on parts as constraints" begin
        t = hohlraum(; geom..., nz = 8, nphi = 16, capsule_ntheta = 6, capsule_nphi = 8)
        aims = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20), :inner30 => (0.60, 3.63),
                    :inner23 => (1.09, 3.09))
        quads = map(nif_quads(power = 8.0)) do q
            r, z = aims[q.group]
            φ = port_angles(q)[2]
            remake(q; aim = (-r * cosd(φ), -r * sind(φ), sign(q.port[3]) * z))
        end
        p = GroupPointing(quads; fields = (:z, :phi))
        s = Scenario(t, p; observe = :capsule)
        xp = initial_parameters(p)
        ref = evaluate(s, xp)
        caps = (:cap_top, :cap_bottom)
        hw = HardwareKeepOut(t, caps; targets = (:wall, :capsule), minimum = 0.05)
        ak = ApertureKeepOut(t; minimum = 0.05)
        sk = ApertureKeepOut(t; apertures = (:leh_top,))       # Inf for the lower hemisphere

        # values against the quantities; one component per beam, labeled by its name
        x = xp .+ 0.05 .* sin.(1:length(xp))
        ev = evaluate(s, x)
        cf, ca, cs = clearance_floor(hw), clearance_floor(ak; name = :aperture, minimum = 0.02),
                     clearance_floor(sk; name = :top)
        @test cf.f(ev) == hardware_clearance(hw, ev.beams) && ca.f(ev) == aperture_clearance(ak, ev.beams)
        @test !cf.upper && cf.limit == [0.05] && cf.scale == [1.0] && ca.limit == [0.02]
        rows = check_constraints(ev, (cf, ca, cs))
        @test length(rows) == 3 * 48
        @test [r.name for r in rows[1:48]] == ["clearance[$(b.name)]" for b in quads]
        @test [r.value for r in rows[49:96]] == aperture_clearance(ak, ev.beams)
        @test all(r -> r.violation ≈ (0.05 - r.value) / 1.0, filter(r -> isfinite(r.value), rows[1:48]))
        # A beam that the keep-out does not apply to (Inf, with no aperture on its side) is
        # satisfied.
        top = aperture_clearance(sk, ev.beams)
        @test count(isinf, top) == 24 && all(r -> r.status === :satisfied, rows[97:end][isinf.(top)])
        @test_throws MethodError clearance_floor(42)
        # any keep-out with a `keepout_clearance` method and a `minimum`: a capture rule
        kc = CaptureKeepIn(t, :wall; minimum = 0.0)
        cc = clearance_floor(kc; name = :capture)
        @test cc.f(ev) == capture_margin(kc, ev.beams) && cc.limit == [0.0]
        # The display shows a constraint of more than eight components as its worst row.
        shown = sprint(show, MIME"text/plain"(), rows)
        worst = rows[argmax([r.violation for r in rows[1:48]])]
        @test count(==('\n'), shown) == 3 && occursin("the worst of 48", shown) && occursin(worst.name, shown)
        @test count(==('\n'), sprint(show, MIME"text/plain"(), rows[1:5])) == 5
        lc = light_ceiling(ref, caps...)
        cp = component_power(ev.deposition, t)
        onCaps = sum(get(cp, k, 0.0) for k in (:cap_top, :cap_bottom, :cap_top_outside, :cap_bottom_outside))
        @test lc.f(ev) ≈ onCaps / ev.deposition.delivered rtol = 1e-12
        @test light_fraction(ev, caps...) == power_fraction(ev, caps...) + clipped_fraction(ev, caps...)
        @test only(lc.limit) == light_fraction(ref, caps...) && lc.upper
        @test only(lc.scale) == max(only(lc.limit), 1e-3)
        lf = light_floor(0.25, :wall; name = :capture)
        @test !lf.upper && lf.limit == [0.25] && lf.scale == [0.25]
        @test check_constraints(ev, lf)[1].name == "capture[wall]"
        @test_throws ArgumentError light_ceiling(ref, :no_such_part)
        @test_throws ArgumentError light_ceiling(1.5, :wall)

        # gradients of the values through the optimizer's duals against central differences
        d = cos.(1:length(xp))
        along(f) = τ -> f(evaluate(s, x .+ τ .* d))
        for f in (ev -> sum(cf.f(ev)[1:8] .* (1:8)), lc.f, ev -> light_fraction(ev, :capsule))
            ad = ForwardDiff.derivative(along(f), 0.0)
            fd = (along(f)(1e-6) - along(f)(-1e-6)) / 2e-6
            @test ad ≈ fd rtol = 1e-6
        end

        # constrained runs with one free parameter (the 50° cone's z; the rest held),
        # whose objective pulls the cone's aim toward the center: the constraint ends
        # active, at the parameter where its value reaches the limit
        h = HeldParameters(p, n -> n != "outer50.z")
        sh = remake(s; parameterization = h)
        z0 = initial_parameters(h)
        @test only(z0) ≈ 5.01 rtol = 1e-12
        pull(ev) = sum(b.aim[3]^2 for b in ev.beams if b.group === :outer50) / 25
        at(z) = evaluate(sh, [z])
        function root(g, a, b)              # bisection for g = 0, g(a) < 0 < g(b)
            @assert g(a) < 0 < g(b)
            for _ in 1:60
                m = (a + b) / 2
                g(m) < 0 ? (a = m) : (b = m)
            end
            return (a + b) / 2
        end
        # The keep-out floor is the smallest clearance at the start less 0.1 mm, and the 50°
        # cone reaches it first.
        least(z) = minimum(filter(isfinite, hardware_clearance(hw, at(z).beams)))
        floor_ = clearance_floor(hw; minimum = least(5.01) - 0.1)
        zk = root(z -> least(z) - only(floor_.limit), 4.0, 5.01)
        rk = optimize_pointing(sh, pull; lower = [4.0], upper = [5.3], constraints = (floor_,),
                               feasibility = 1e-6)
        @test rk.outcome === :interior && abs(only(rk.x) - zk) < 1e-5
        act = filter(r -> r.status === :active, rk.constraints)
        @test !isempty(act) && all(r -> startswith(r.name, "clearance[Q") && r.multiplier > 0, act)
        @test all(r -> abs(r.value - r.limit) <= 1e-6, act)
        @test all(r -> r.status !== :violated, rk.constraints)
        # The light on the caps has a ceiling at 0.1% of the delivered power.
        ceil_ = light_ceiling(1e-3, caps...)
        zl = root(z -> 1e-3 - light_fraction(at(z), caps...), 4.0, 5.01)
        rl = optimize_pointing(sh, pull; lower = [4.0], upper = [5.3], constraints = [ceil_],
                               feasibility = 1e-6)
        @test rl.outcome === :interior && abs(only(rl.x) - zl) < 1e-5
        @test only(rl.constraints).status === :active && only(rl.constraints).multiplier > 0
        @test abs(only(rl.constraints).value - 1e-3) <= 1e-6 * 1e-3

        # The pieces compose. A constrained run uses held parameters of power scales per
        # cone, with the keep-out and the light as constraints, and the result is
        # summarized.
        ps = HeldParameters(PowerScales(p), n -> startswith(n, "inner") && !endswith(n, "power_scale"))
        sc = remake(s; parameterization = ps)
        y0 = initial_parameters(sc.param)
        @test length(y0) == 8 && parameter_names(sc.param)[end] == "outer50.power_scale"
        lo, hi = bounds(sc.param, y0; z = 0.3, phi = 0.3, scale = 0.2)
        total = Constraint(:power, ev -> sum(b.power for b in ev.beams); upper = 384.0)
        rc = optimize_pointing(sc, drive_objective(; drive = mean_flux(ref)); lower = lo, upper = hi,
                               constraints = (total, clearance_floor(hw), light_ceiling(ref, :capsule)),
                               iterations = 40, evaluations = 300)
        again = apply(sc.param, rc.x)
        @test rc.names == parameter_names(sc.param)
        @test [(b.aim, b.power) for b in rc.beams] == [(b.aim, b.power) for b in again]
        @test length(rc.constraints) == 2 + 48 && rc.constraints[1].name == "power"
        @test rc.constraints[end].name == "light[capsule]"
        sm = pointing_summary(sc, rc; keepout = hw)
        c = hardware_clearance(hw, rc.beams)
        @test sm.clearance == [minimum(filter(isfinite, c))] && sm.margin == sm.clearance .- 0.05
        @test sm.power ≈ sum(b.power for b in rc.beams) rtol = 1e-14
        start = apply(sc.param, y0)
        @test sm.reference_clearance == [minimum(filter(isfinite, hardware_clearance(hw, start)))]
        @test sm.shift_max ≈ maximum(norm(b.aim - a.aim) for (b, a) in zip(rc.beams, start)) rtol = 1e-14
    end

    # Under the VisRad hard cut (`smooth_cut`, the default), the run uses the profile
    # tapered below the cut on the points of the same operator and is reported under the
    # scenario cut. The reference is the same run on a scenario built with the tapered
    # profile and a newly built operator, and the two are identical.
    @testset "VisRad's hard cut optimized tapered" begin
        t = hohlraum(; geom..., nz = 8, nphi = 16, capsule_ntheta = 6, capsule_nphi = 8)
        quads = nif_quads(power = 8.0)
        p = GroupPointing(quads; fields = (:z,))
        dk(taper) = (; profile = VisRadProfile(; taper), resolution = 0.6, reach = 0.5, level = 0)
        hard = Scenario(t, p; observe = :capsule, deposition_kw = dk(0.0))
        soft = remake(hard; deposition_kw = dk(1.0))
        @test soft.deposition.pts == hard.deposition.pts          # the taper keeps the support
        ref = evaluate(hard, quads)
        obj = drive_objective(; drive = mean_flux(ref))
        cs = (irradiance_ceiling(ref, :wall),)
        x0 = initial_parameters(p) .+ 0.1
        kw = (; x0, constraints = cs, iterations = 30, evaluations = 300)
        r = optimize_pointing(hard, obj; kw...)
        rs = optimize_pointing(soft, obj; kw...)
        r0 = optimize_pointing(hard, obj; kw..., smooth_cut = false)
        @test r.x == rs.x && r.beams == rs.beams && r.evaluations == rs.evaluations
        @test r.tapered_value == rs.value && isnan(rs.tapered_value) && isnan(r0.tapered_value)
        ev = evaluate(hard, r.x)
        @test r.value == obj(ev) && r.initial_value == obj(evaluate(hard, x0))
        @test [c.value for c in r.constraints] == [peak_irradiance(ev, :wall)]
        @test [c.multiplier for c in r.constraints] == [c.multiplier for c in rs.constraints]
        @test r.outcome === Briareus._outcome(r.reason, r.constraints, r.at_bounds)
        @test occursin("tapered", sprint(show, MIME"text/plain"(), r)) &&
              !occursin("tapered", sprint(show, MIME"text/plain"(), r0))
        # In a phased scenario each phase is tapered, and the report uses the hard cut.
        ps = PhasedScenario(p, Phase(:a, hard; power = 8.0), Phase(:b, hard; power = 4.0, scaled = false))
        @test all(sc -> sc.deposition.profile.taper == 1, Briareus._tapered(ps).scenarios)
        pobj = ev -> obj(ev[:a])
        rp = optimize_pointing(ps, pobj; iterations = 3, evaluations = 30)
        @test rp.value == pobj(evaluate(ps, rp.x)) && isfinite(rp.tapered_value)
    end
end
