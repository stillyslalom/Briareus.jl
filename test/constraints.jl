# Constraints in `optimize_pointing` (an augmented Lagrangian around L-BFGS-B), the
# wall-irradiance ceiling and the per-degree mode gate.
#
# References: the closed-form solution and multiplier of an equality-active quadratic
# program (the drive's coefficients are linear in the quad powers, `influence_matrix`);
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

        # the box keeps every power above the cap's share: infeasible at the bounds
        boxed = optimize_pointing(s, f; x0 = P0, lower = 0.99 .* P0, upper = hi, constraints = (total,))
        @test boxed.outcome === :infeasible_box && boxed.reason === :converged
        @test boxed.constraints[1].status === :violated && length(boxed.at_bounds) == 48
        @test occursin("ended with constraints violated", sprint(show, MIME"text/plain"(), boxed))
        # the budget ends a run first
        short = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi, evaluations = 10, constraints = [total])
        @test short.outcome === :evaluations && short.stopped && short.evaluations == 10
        # without constraints: the plain run, its outcome from the box alone
        plain = optimize_pointing(s, f; x0 = P0, lower = lo, upper = hi)
        @test plain.outcome === :interior && isempty(plain.constraints) && plain.outer == 0
        @test sum(plain.x) > cap
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

        # the exact peak held (`exact = true`): the stand-in's working limit lowered until the
        # peak over the elements is within `feasibility` of the limit. The stand-in alone
        # does not bound the peak to a fixed margin: designs with the same gates differ in
        # the peak by several percent (from this start, 0.997 and 1.030 of the limit with
        # the norm at 0.974 and at the limit, under two roundings: research notes).
        exact = irradiance_ceiling(ref, :wall; exact = true)
        res = optimize_pointing(rp, drive_objective(; drive = mean_flux(ref)); x0, lower, upper,
                                constraints = (exact, gates), iterations = 1000, evaluations = 3000)
        ev = evaluate(rp, res.x)
        @test res.outcome in (:interior, :box)
        @test all(r -> r.status !== :violated, res.constraints)
        @test all(mode_departure(ev, c48) .<= gates.limit .* (1 + 1e-3))
        @test peak_irradiance(ev, :wall) <= pk * (1 + 1e-3)
        @test res.value <= res.initial_value + 1e-3
        # the rows: the exact ceiling reports the peak against the limit, the plain one its norm
        @test res.constraints[1].value == peak_irradiance(ev, :wall) && res.constraints[1].limit == pk
        @test check_constraints(ev, exact)[1].value == peak_irradiance(ev, :wall)
        @test check_constraints(ev, ceiling)[1].value == ceiling.f(ev)
    end
end
