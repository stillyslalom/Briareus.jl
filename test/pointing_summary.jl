# `pointing_summary`: aim shifts against closed forms (rigid shifts per quad), power and
# cone fraction, keep-out clearances and the clipped fraction as the functions that
# define them give them, and the summary of an optimizer's result and of a subset.

using Briareus, LinearAlgebra, Test

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
    # display: one line
    shown = sprint(show, sm)
    @test !occursin('\n', shown) && occursin("aim shifts mean", shown) && occursin("cone fraction 0.333", shown)

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
    ds = sprint(show, d)
    @test !occursin('\n', ds) && occursin("of the reference, departure (l = 1–4)", ds)
end
