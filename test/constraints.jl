# Constraints in `optimize_pointing` (an augmented Lagrangian around L-BFGS-B), the
# wall-irradiance ceiling and the per-degree mode gate.
#
# References: the closed-form solution and multiplier of an equality-active quadratic
# program (the drive's coefficients are linear in the quad powers, `influence_matrix`);
# for phase 1's relaxations, fills of the quads in decreasing drive per power (the
# largest drive at a total power, the least power for a drive);
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

        # the departure: one coefficient off by δ c₀₀ departs by √(2l+1) |δ| at that
        # degree and nowhere else, whatever the two shapes' overall levels
        r = [1.0; 0.01 .* sin.(1:24)]
        x = 3.0 .* r
        x[mode_index(2, 1)] += 3.0 * 0.004
        d = mode_departure(x, r; lmax = 4)
        @test d[2] ≈ sqrt(5) * 0.004 rtol = 1e-12
        @test all(abs.(d[[1, 3, 4]]) .< 1e-15)
        @test mode_departure(r, r; lmax = 4) == zeros(4)
        # an axisymmetric difference: |aₗ/a₀ - aₗ'/a₀'| (c_l0 = aₗ √(4π/(2l+1)))
        y = copy(r)
        y[mode_index(3, 0)] += 0.02 * sqrt(4π / 7) / sqrt(4π) * r[1]
        @test mode_departure(y, r; lmax = 3)[3] ≈ 0.02 rtol = 1e-12
        @test_throws ArgumentError mode_departure(r[1:9], r; lmax = 4)
        @test_throws ArgumentError mode_gate(r, [0.01, 0.0])
        @test mode_gate(r, (; rms_mode = [0.001, 0.002])).limit ≈ [sqrt(3) * 0.001, sqrt(5) * 0.002]
    end

    # A quadratic program with a closed form: quad powers fitted to 1.1 times the
    # 48-quad drive's coefficients through l = 2 (linear in the powers: the influence
    # matrix M), regularized toward the start, under a ceiling on the total power that
    # the fit would exceed. KKT: H P + ν 1 = r, 1ᵀP = cap.
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
        # a ceiling on a quantity of the beams alone, and a floor that stays inactive
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
        # a start already feasible and on the cap, the objective scaled so that its gradient
        # lies below the first inner tolerance (1e-3): the run still reaches the solution
        # rather than ending at the start on the first inner run, with zero multipliers
        small(ev) = 1e-4 * f(ev)
        Pa = P0 .* (cap / sum(P0))
        on_cap = optimize_pointing(s, small; x0 = Pa, lower = lo, upper = hi, constraints = (total,))
        @test on_cap.outcome === :interior && norm(on_cap.x - Pstar) / norm(Pstar) < 1e-3
        @test norm(Pa - Pstar) / norm(Pstar) > 0.05

        # warm restarts. Split at an outer iteration (each inner run starts its L-BFGS-B
        # afresh, from the state the result carries), the continued run is the
        # uninterrupted one bit for bit. Split inside an inner run by the evaluation budget,
        # the restart loses that run's L-BFGS-B memory and takes another path to the same
        # strictly convex optimum: the designs agree to the convergence tolerance
        # (feasibility 1e-8; the run lies 1e-8 from the KKT point), tested at 1e-6.
        cs2 = (total, drive)
        a = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = cs2,
                              feasibility = 1e-8, outer = 2)
        @test a.outcome === :iterations && a.outer == 2 && a.state.λ == [r.multiplier for r in a.constraints]
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
        # with phase 1 the same verdict, and the cap relaxed alone (no other constraint) to
        # the least total power the box allows, every power at its lower bound
        bv = optimize_pointing(s, f; x0 = P0, lower = 0.99 .* P0, upper = hi, constraints = (total,),
                               relaxation = true)
        @test bv.outcome === :infeasible_box && only(bv.relaxations).status === :relaxed
        @test only(bv.relaxations).violation ≈ 0.99 * sum(P0) / (0.95 * sum(P0)) - 1 rtol = 1e-9

        # Phase 1 and the relaxations on a conflict with closed forms. The drive
        # (`mean_flux`) is linear in the powers, `m` per quad by unit powers, so the largest
        # drive at a total power and the least total power for a drive are fills of the
        # quads in decreasing `m` to the box's 20 TW. A cap of 120 TW and a drive floor 5 %
        # above the cap's largest drive: the floor relaxed alone falls to that drive (by
        # 1 - 1/1.05 of itself), the cap relaxed alone rises to the floor's least power.
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
        @test power_for(Dfl) < 20 * count(>(0), m)       # the floor alone is feasible
        today = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, constraints = (cap120, above),
                                  feasibility = 1e-8)
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
        # no verdict on a feasible run from far outside the cap: the KKT point as without
        # phase 1, its main run started from phase 1's first design within the limits
        far = optimize_pointing(s, f; x0 = clamp.(2.5 .* P0, lo, hi), lower = lo, upper = hi,
                                constraints = cs2, feasibility = 1e-8, phase1 = true)
        @test far.outcome === :interior && far.phase1.feasible && !far.phase1.verdict
        @test norm(far.x - Pstar) / norm(Pstar) < 1e-4
        @test res.phase1 === nothing && isempty(res.relaxations)

        # the budget ends a run first
        short = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, evaluations = 10, constraints = [total])
        @test short.outcome === :evaluations && short.stopped && short.evaluations == 10
        # without constraints: the plain run, its outcome from the box alone
        plain = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi)
        @test plain.outcome === :interior && isempty(plain.constraints) && plain.outer == 0
        @test sum(plain.x) > cap
        @test isempty(plain.bound_sensitivity) && isempty(plain.bound_gradient)

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
            # the objective's own derivative, against its central difference at the result
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
        # so does a mode gate from a reference evaluation; its values are those from the
        # reference's coefficients, which record no model (unchecked)
        gate_ev = mode_gate(evaluate(sb, P0), [0.01, 0.01])
        gate_c = mode_gate(flux_ylm(evaluate(sb, P0); lmax = 2), [0.01, 0.01])
        xb = P0 .* (1 .+ 0.1 .* sin.(1:48))
        @test gate_ev.f(evaluate(sb, xb)) == gate_c.f(evaluate(sb, xb)) && gate_ev.limit == gate_c.limit
        @test gate_ev.f(evaluate(remake(s; lpi = bs(0.1)), xb)) == gate_ev.f(evaluate(sb, xb))
        @test_throws ArgumentError gate_ev.f(evaluate(s, xb))
        @test_throws ArgumentError gate_ev.f(evaluate(remake(sb; time = 2.0), xb))
        @test length(gate_c.f(evaluate(s, xb))) == 2
        @test mode_gate(evaluate(sb, P0), (; rms_mode = [0.001, 0.002])).limit ≈ [sqrt(3) * 0.001, sqrt(5) * 0.002]
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

        # the irradiance: power means below the peak, rising with p; p = 1 the area mean;
        # homogeneous in the beams' powers; a weight of one reproduces the deposition and
        # a weight of zero on the inner cones their deposition at zero power
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

        # the start: every gate violated, the ceiling at its limit (the freed quads
        # do not reach the hottest elements)
        start = check_constraints(evaluate(rp, x0), (ceiling, gates))
        @test count(r -> r.status === :violated, start) == 4 && start[1].status !== :violated
        @test peak_irradiance(evaluate(rp, x0), :wall) ≈ pk rtol = 1e-12

        # the exact peak held (`exact = true`, the default): the stand-in's working limit
        # lowered until the peak over the elements is within `feasibility` of the limit. The
        # stand-in alone (`exact = false`) does not bound the peak to a fixed margin: designs
        # with the same gates differ in the peak by several percent (from this start, 0.997
        # and 1.030 of the limit with the norm at 0.974 and at the limit, under two
        # roundings: research notes). Under `--check-bounds=yes` the run reaches another
        # design of the same flat family: only what every admissible design meets is tested.
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
        # the rows: the exact ceiling reports the peak against the limit, the plain one its norm
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
        # a beam the keep-out does not apply to (Inf: no aperture on its side) holds
        top = aperture_clearance(sk, ev.beams)
        @test count(isinf, top) == 24 && all(r -> r.status === :satisfied, rows[97:end][isinf.(top)])
        @test_throws ArgumentError clearance_floor(42)
        # any keep-out with a `keepout_clearance` method and a `minimum`: a capture rule
        kc = CaptureKeepIn(t, :wall; minimum = 0.0)
        cc = clearance_floor(kc; name = :capture)
        @test cc.f(ev) == capture_margin(kc, ev.beams) && cc.limit == [0.0]
        # the display: a constraint of more than eight components as its worst row
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
        # the keep-out: the smallest clearance at the start less 0.1 mm, reached by the 50°
        # cone first
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
        # the light on the caps: a ceiling at 0.1 % of the delivered power
        ceil_ = light_ceiling(1e-3, caps...)
        zl = root(z -> 1e-3 - light_fraction(at(z), caps...), 4.0, 5.01)
        rl = optimize_pointing(sh, pull; lower = [4.0], upper = [5.3], constraints = [ceil_],
                               feasibility = 1e-6)
        @test rl.outcome === :interior && abs(only(rl.x) - zl) < 1e-5
        @test only(rl.constraints).status === :active && only(rl.constraints).multiplier > 0
        @test abs(only(rl.constraints).value - 1e-3) <= 1e-6 * 1e-3

        # the pieces compose: held parameters of power scales per cone in a constrained
        # run, with the keep-out and the light as constraints, and the summary of the result
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
end
