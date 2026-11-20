# `pointing_summary`: aim shifts against closed forms (rigid shifts per quad), power and
# cone fraction, keep-out clearances and the clipped fraction as the functions that
# define them give them, and the summary of an optimizer's result and of a subset; the
# smallest-shift objective against closed forms (below).

using Briareus, LinearAlgebra, ForwardDiff, Test

# A keep-out whose rule cannot place the first beam (-Inf), does not apply to the second
# (+Inf) and clears the others by 0.2 mm, or with `fail = false` the first by 0.2 mm too.
struct _SummaryKeepOut
    minimum::Float64
    fail::Bool
end
Briareus.keepout_clearance(k::_SummaryKeepOut, beams::AbstractVector{<:Beam}) =
    [i == 1 && k.fail ? -Inf : i == 2 ? Inf : 0.2 for i in eachindex(beams)]

@testset "pointing summary" begin
    t = hohlraum(nz = 8, nphi = 16, capsule_ntheta = 6, capsule_nphi = 8)
    quads = hybride_mock_pointing(nif_quads(power = 8.0))
    p = GroupPointing(quads; by = :quad, fields = (:z, :phi), mirror = false, relative = true)
    s = Scenario(t, p)
    names = parameter_names(p)
    # every quad up 0.3 mm, one also 0.4 mm along its φ̂: shifts 0.3 and 0.5 mm
    x = [endswith(n, ".z") ? 0.3 : n == "Q31B.phi" ? 0.4 : 0.0 for n in names]
    beams = apply(p, x)
    ko = ApertureKeepOut(t; minimum = 0.1)
    sm = pointing_summary(s, beams; keepout = ko)
    @test sm.nbeams == 48
    @test sm.shift_max ≈ 0.5 rtol = 1e-12
    @test sm.shift_mean ≈ (47 * 0.3 + 0.5) / 48 rtol = 1e-12
    @test sm.power == 384.0
    @test sm.cone_fraction ≈ 16 / 48 rtol = 1e-14
    c, c0 = aperture_clearance(ko, beams), aperture_clearance(ko, quads)
    @test sm.clearance == [minimum(filter(isfinite, c))]
    @test sm.reference_clearance == [minimum(filter(isfinite, c0))]
    @test sm.margin == sm.clearance .- 0.1
    ev = evaluate(s, beams)
    @test (sm.clipped, sm.missed) == (clipped_fraction(ev), missed_fraction(ev))
    # several keep-outs, none
    two = pointing_summary(s, beams; keepout = (ko, SphereKeepOut(t)))
    @test length(two.clearance) == 2 && two.clearance[1] == sm.clearance[1]
    @test isempty(pointing_summary(s, beams).clearance)
    @test pointing_summary(s, quads).shift_max < 1e-12      # the start rebuilt by `apply`
    @test pointing_summary(s, quads; reference = quads).shift_max == 0.0
    # The display is one line.
    shown = sprint(show, sm)
    @test !occursin('\n', shown) && occursin("aim shifts mean", shown) && occursin("cone fraction 0.333", shown)
    # a beam the rule cannot place (-Inf) is the least clearance and shows; a beam it does
    # not apply to (+Inf) is left out
    fails = pointing_summary(s, beams; keepout = _SummaryKeepOut(0.1, true))
    @test fails.clearance == [-Inf] && fails.margin == [-Inf] && fails.reference_clearance == [-Inf]
    @test occursin("clearance -Inf mm", sprint(show, fails))
    @test pointing_summary(s, beams; keepout = _SummaryKeepOut(0.1, false)).clearance == [0.2]

    # a subset's result, summarized against the full set's start (beams matched by name)
    keep = findall(b -> nif_bundle(b) != "31", quads)
    sk = subset(s, keep)
    r = optimize_pointing(sk, ev -> clipped_fraction(ev); iterations = 1)
    rs = pointing_summary(sk, r; reference = quads, keepout = ko)
    bs = pointing_summary(sk, r.beams; reference = quads, keepout = ko)
    @test rs.nbeams == length(keep) == bs.nbeams
    @test (rs.shift_mean, rs.shift_max, rs.power, rs.clearance, rs.clipped) ==
          (bs.shift_mean, bs.shift_max, bs.power, bs.clearance, bs.clipped)
    @test_throws ArgumentError pointing_summary(s, [remake(quads[1]; name = "X1")])
    @test isnan(sm.drive_ratio) && isempty(sm.departure) && !occursin("drive", shown)

    # the capsule drive against a reference design: the ratio of the solid-angle means
    # (`mean_flux`, by the Legendre weights) and the departure per degree, by an
    # evaluation or by its coefficients alike; a design against itself departs by zero
    so = Scenario(t, p; observe = :capsule)
    ref = evaluate(so, quads)
    d = pointing_summary(so, beams; drive_reference = ref)
    @test d.drive_ratio ≈ mean_flux(evaluate(so, beams)) / mean_flux(ref) rtol = 1e-10
    @test d.departure == mode_departure(evaluate(so, beams), flux_ylm(ref; lmax = 4)) && length(d.departure) == 4
    dc = pointing_summary(so, beams; drive_reference = flux_ylm(ref; lmax = 2), lmax = 2)
    @test dc.drive_ratio ≈ d.drive_ratio rtol = 1e-12
    @test dc.departure ≈ d.departure[1:2] rtol = 1e-12
    same = pointing_summary(so, quads; drive_reference = ref)
    @test same.drive_ratio ≈ 1 rtol = 1e-12
    @test all(same.departure .< 1e-12)
    # The split by order is checked by computing the m = 0 part by hand from the
    # coefficients, and the squares of the parts add to the whole.
    cd, cr = flux_ylm(evaluate(so, beams); lmax = 4), flux_ylm(ref; lmax = 4)
    zonal = [sqrt(2l + 1) * abs(cd[mode_index(l, 0)] / cd[1] - cr[mode_index(l, 0)] / cr[1]) for l in 1:4]
    @test d.departure_zero ≈ zonal rtol = 1e-12
    @test d.departure_zero .^ 2 .+ d.departure_nonzero .^ 2 ≈ d.departure .^ 2 rtol = 1e-12
    @test d.departure_zero == mode_departure(cd, cr; m = :zero) && d.departure_nonzero == mode_departure(cd, cr; m = :nonzero)
    @test isempty(sm.departure_zero) && isempty(sm.departure_nonzero)
    # On coefficient vectors, a change of only c_11 is m ≠ 0 (√3 Δ at l = 1), and a
    # change of only c_20 is m = 0.
    c1 = [1.0; zeros(24)]
    c2 = copy(c1); c2[mode_index(1, 1)] = 0.004; c2[mode_index(2, 0)] = -0.003
    @test mode_departure(c2, c1; m = :zero) ≈ [0, sqrt(5) * 0.003, 0, 0] atol = 1e-16
    @test mode_departure(c2, c1; m = :nonzero) ≈ [sqrt(3) * 0.004, 0, 0, 0] atol = 1e-16
    @test mode_departure(c2, c1; m = :all) == mode_departure(c2, c1)
    @test_throws ArgumentError mode_departure(c2, c1; m = :odd)
    ds = sprint(show, d)
    @test !occursin('\n', ds) && occursin("of the reference, departure (l = 1–4)", ds)
    @test isapprox(sm.shift_rms, sqrt((47 * 0.3^2 + 0.5^2) / 48); rtol = 1e-12) && occursin("rms 0.30", shown)
    # groups moved past thresholds, by hand: 47 quads by 0.3 mm and Q31B by 0.5 mm; per cone,
    # every cone by 0.3 mm and Q31B's by 0.5 mm
    mq = pointing_summary(s, beams; moved = (0.25, 0.4, 0.6))
    @test mq.moved == [48, 1, 0] && mq.ngroups == 48 && mq.moved_thresholds == [0.25, 0.4, 0.6]
    mc = pointing_summary(s, beams; moved = [0.25, 0.4], moved_by = :group)
    @test mc.moved == [4, 1] && mc.ngroups == 4
    @test pointing_summary(s, beams; moved = (0.4,), moved_by = b -> b.port[3] > 0).moved == [1]
    @test isempty(sm.moved) && sm.ngroups == 0
    @test occursin("quads moved > 0.25 / 0.4 / 0.6 mm: 48 / 1 / 0 of 48", sprint(show, mq))
    @test_throws ArgumentError pointing_summary(s, beams; moved = (-0.1,))
end

# `shift_objective`: the weighted mean squared aim shift against closed forms, its
# gradient against the closed form 2 w Δ / Σw, and the least re-pointing that meets one
# linear constraint on the aims, whose solution has a closed form: with the mean axial aim
# of a set S of quads raised by δ, each quad of S shifts by δ |S| (1/w) / Σ_S (1/w) and
# every other quad stays.
@testset "shift objective" begin
    t = hohlraum(nz = 8, nphi = 16, capsule_ntheta = 6, capsule_nphi = 8)
    quads = hybride_mock_pointing(nif_quads(power = 8.0))
    p = GroupPointing(quads; by = :quad, fields = (:z,), mirror = false, relative = true)
    s = Scenario(t, p)
    names = parameter_names(p)
    x = [0.01 * i * (-1)^i for i in eachindex(names)]
    ev = evaluate(s, x)
    o = shift_objective(s)
    @test o(ev) ≈ sum(abs2, x) / 48 rtol = 1e-12
    @test o(ev) ≈ pointing_summary(s, ev.beams).shift_rms^2 rtol = 1e-12
    dearer = b -> b.group === :outer44 ? 3.0 : 1.0
    w = [dearer(q) for q in quads]
    ow = shift_objective(quads; weight = dearer)
    @test ow(ev) ≈ sum(w .* x .^ 2) / sum(w) rtol = 1e-12
    gw = ForwardDiff.gradient(y -> ow(evaluate(s, y)), x)
    @test gw ≈ 2 .* w .* x ./ sum(w) rtol = 1e-12
    # Beams are matched by name, both for a subset against the full set and for a beam
    # without a reference beam.
    keep = findall(b -> b.group !== :inner23, quads)
    @test o(evaluate(subset(s, keep), ev.beams[keep])) ≈ sum(abs2, x[keep]) / length(keep) rtol = 1e-12
    @test_throws ArgumentError shift_objective(quads[1:4])(ev)
    @test_throws ArgumentError shift_objective(quads; weight = b -> -1.0)
    @test_throws ArgumentError shift_objective([quads; quads[1]])
    @test_throws ArgumentError shift_objective(quads; weight = b -> 0.0)(ev)

    # the least re-pointing for a mean axial lift δ of the upper outer quads
    top(b) = b.port[3] > 0 && b.group in (:outer44, :outer50)
    S = findall(top, quads)
    δ = 0.1
    zbar(ev) = sum(b.aim[3] for b in ev.beams if top(b)) / length(S)
    lift = Constraint(:lift, zbar; lower = zbar(evaluate(s, zeros(48))) + δ)
    r = optimize_pointing(s, ow; constraints = (lift,), feasibility = 1e-8)
    expect = [i in S ? δ * length(S) / w[i] / sum(1 ./ w[S]) : 0.0 for i in 1:48]
    @test r.outcome === :interior && maximum(abs, r.x - expect) < 1e-5
    @test r.value ≈ sum(w .* expect .^ 2) / sum(w) rtol = 1e-4
    @test only(r.constraints).status === :active
end
