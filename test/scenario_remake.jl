# `remake(::Scenario)`: each replaced part gives the scenario built from scratch with it,
# what does not depend on the change is shared, and changes that need a new scenario are
# refused.

using Briareus, Test

@testset "remake(::Scenario)" begin
    t = hohlraum(nz = 8, nphi = 16, capsule_ntheta = 6, capsule_nphi = 8)
    quads = nif_quads(power = 8.0)
    p = GroupPointing(quads; fields = (:z, :r, :power))
    x0 = parameters(p; inner23 = (z = 3.09, r = -1.09), inner30 = (z = 3.63, r = -0.60),
                    outer44 = (z = 5.20,), outer50 = (z = 5.01,))
    s = Scenario(t, p; observe = :capsule)
    moments(sc, x) = flux_moments(evaluate(sc, x))

    # an evaluation computes its moments and Y_lm coefficients once: every read is the
    # scenario's response product with the deposition, a copy of the stored one
    ev = evaluate(s, x0)
    E = ev.deposition.element
    m1, y1 = flux_moments(ev; lmax = 4), flux_ylm(ev; lmax = 3)
    @test m1 == (s.moment_response*E)[1:5] && y1 == (s.ylm_response*E)[1:16]
    m1[1] = y1[1] = 0.0
    @test mean_flux(ev) == (s.moment_response*E)[1] && flux_ylm(ev) == s.ylm_response * E

    # the laser–plasma model: nothing rebuilt
    lpi = CalibratedTransfer((:outer44, :outer50), (:inner23, :inner30), linear_transfer(0.0, 0.05))
    sl = remake(s; lpi)
    @test sl.lpi === lpi && sl.deposition === s.deposition && sl.model === s.model
    @test sl.moment_response === s.moment_response
    @test moments(sl, x0) == moments(Scenario(t, p; observe = :capsule, model = s.model, lpi), x0)
    @test remake(sl; lpi = nothing).lpi === nothing

    # a parameterization on the same ports keeps the operator (built where the original
    # beams aimed), so it evaluates as the original scenario on the same beams
    pp = GroupPointing(apply(p, x0); by = :quad, fields = (:power,))
    sp = remake(s; parameterization = pp)
    @test sp.param === pp && sp.deposition === s.deposition
    xp = initial_parameters(pp) .* 1.1
    @test moments(sp, xp) == flux_moments(evaluate(s, apply(pp, xp)))
    # other ports (or another order) need the operator rebuilt: refused without deposition_kw
    @test_throws ArgumentError remake(s; parameterization = GroupPointing(reverse(quads)))
    @test_throws ArgumentError remake(s; parameterization = GroupPointing(omega_beams()))
    @test_throws ArgumentError remake(s; parameterization = GroupPointing(nif_quads(fnumber = 10.0)))

    # deposition options: the operator rebuilt for the parameterization's beams
    kw = (; ncone = 7, resolution = 0.4)
    sd = remake(s; deposition_kw = kw)
    @test sd.deposition !== s.deposition && sd.moment_response === s.moment_response
    @test evaluate(sd, x0).deposition.element ==
          evaluate(Scenario(t, p; deposition_kw = kw), x0).deposition.element
    so = remake(s; parameterization = GroupPointing(omega_beams()), deposition_kw = (;))
    @test length(so.deposition.ports) == length(omega_beams())

    # time, albedos, observed elements and the model: the responses recomputed
    α = fill(0.8, nelements(t))
    st = remake(s; time = 1.0, albedo = α)
    ref = Scenario(t, p; observe = :capsule, model = s.model, time = 1.0, albedo = α)
    @test st.time == 1.0 && st.albedo == α && st.deposition === s.deposition
    @test st.moment_response ≈ ref.moment_response rtol = 1e-12
    @test st.ylm_response ≈ ref.ylm_response rtol = 1e-12
    @test st.response ≈ ref.response rtol = 1e-12
    sw = remake(s; observe = :wall)
    @test sw.observed == elements(t, :wall)
    @test sw.moment_weights == Scenario(t, p; observe = :wall, model = s.model).moment_weights
    sn = remake(s; model = nothing)
    @test sn.model === nothing && sn.response === nothing && sn.observed == s.observed

    # a subset keeps its restricted operator
    keep = findall(b -> nif_bundle(b) != "31", quads)
    ss = remake(subset(s, keep); lpi)
    @test length(ss.deposition.ports) == length(keep) && ss.lpi === lpi

    # what needs a new scenario
    @test_throws ArgumentError remake(s; target = t)
    @test_throws ArgumentError remake(s; viewfactor_kw = (; closed = true))
    @test_throws ArgumentError remake(s; radiation = false)
    other = Target([foil(:foil, (0, 0, 0), (0, 0, 1); width = 4.0, n = 4)])
    @test_throws ArgumentError remake(s; model = RadiationModel(other))
end
