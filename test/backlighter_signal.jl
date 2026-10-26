# He-α conversion efficiency (kshell_ce), filter attenuation and the photometric chain
# (backlighter_signal), against closed forms and the tabulated data.

using Briareus, LinearAlgebra, Test

@testset "backlighter signal" begin
    # --- the CE fit -----------------------------------------------------------------------
    d = kshell_ce_data()
    m = kshell_ce_model()
    k = findall(==(1.0), d.fit)
    @test m.n == length(k) == 68
    # weights 1/(s² + τ²), s the stated relative error (W01's 10% ⊕ 20% where blank), τ
    # such that χ² per degree of freedom is 1; a constant per class zeroes the weighted
    # residual sum in each class; the scatter is τ ⊕ the rms stated error
    res = [log(kshell_ce(d.hnu[i]; irradiance = 1e15 * d.irradiance[i],
                         prepulse = d.class[i] == "prepulsed").ce / (d.ce[i] / 100)) for i in k]
    s = [isnan(d.u[i]) ? sqrt(0.1^2 + 0.2^2) : d.u[i] / d.ce[i] for i in k]
    w = 1 ./ (s .^ 2 .+ m.tau^2)
    @test sum(w .* res .^ 2) / (68 - 7) ≈ 1 rtol = 1e-9
    pre = [d.class[i] == "prepulsed" for i in k]
    @test abs(sum(w .* res .* pre)) < 1e-8 && abs(sum(w .* res .* .!pre)) < 1e-8
    @test m.scatter ≈ sqrt(m.tau^2 + sum(abs2, s) / 68) rtol = 1e-12
    @test 0.15 < m.tau < 0.35 && 0.3 < m.scatter < 0.45
    # held-out check on independent data: the form fitted to the 2026-10-24 compilation
    # alone (Workman & Kyrala, Barrios et al.) predicts Coppari et al.'s 31 OMEGA points
    # within their stated error and the fit's predictive scatter
    c19 = d.source .== "C19"
    d0 = merge(d, (; fit = [c ? 0.0 : f for (c, f) in zip(c19, d.fit)]))
    m0 = Briareus._fit_kshell_ce(d0)
    @test m0.n == 37
    z = [log(kshell_ce(d.hnu[i]; irradiance = 1e15 * d.irradiance[i], model = m0,
                       prepulse = d.class[i] == "prepulsed").ce / (d.ce[i] / 100)) for i in findall(c19)]
    zs = [kshell_ce(d.hnu[i]; irradiance = 1e15 * d.irradiance[i], model = m0,
                    prepulse = d.class[i] == "prepulsed").sigma for i in findall(c19)]
    @test sqrt(sum(abs2, z ./ hypot.(zs, 0.2)) / length(z)) < 1.2
    # the held-out Ti foil (below the fitted lines) is flagged: the fit overpredicts it
    iTi = findfirst(==("Ti"), d.element)
    @test d.fit[iTi] == 0 && kshell_ce(:Ti; irradiance = 1.32e15).extrapolated
    # the line energies of the data are the tabulated ones
    for (el, hν) in zip(d.element, d.hnu)
        @test kshell_line_energy(Symbol(el)) == hν
    end
    # the optimum is a maximum, with the peak CE there; and σ ≥ the scatter
    for line in (:Fe, :Zn, :Ge, :Zr), prepulse in (false, true)
        r = kshell_ce(line; irradiance = 2e15, prepulse)
        @test r.peak ≥ r.ce
        @test kshell_ce(line; irradiance = r.optimum, prepulse).ce ≈ r.peak rtol = 1e-12
        for f in (0.7, 1.4)
            @test kshell_ce(line; irradiance = f * r.optimum, prepulse).ce < r.peak
        end
        @test r.sigma ≥ m.scatter && r.lo ≈ r.ce * exp(-r.sigma) && r.hi ≈ r.ce * exp(r.sigma)
    end
    # the peak falls with the line as the published Z scalings, hν^-(2.5–4) (Workman &
    # Kyrala 2001) and hν^-(3–4) (Barrios et al. 2013), across Fe–Ge (foil) and Zn–Mo
    # (prepulsed)
    slope(a, b, p) = -log(kshell_ce(b; irradiance = 1e15, prepulse = p).peak /
                          kshell_ce(a; irradiance = 1e15, prepulse = p).peak) /
                     log(kshell_line_energy(b) / kshell_line_energy(a))
    @test 2.5 < slope(:Fe, :Ge, false) < 4.5
    @test 2.5 < slope(:Zn, :Mo, true) < 4.5
    # extrapolation flags and growing σ outside the data
    @test !kshell_ce(:Ge; irradiance = 1.5e15, prepulse = true, pulse = 2).extrapolated
    @test kshell_ce(:Ge; irradiance = 1.5e15, prepulse = true, pulse = 0.5).extrapolated
    @test kshell_ce(:Ti; irradiance = 1e15, prepulse = true).extrapolated
    @test kshell_ce(:Ag; irradiance = 1e16).extrapolated
    @test kshell_ce(:Ag; irradiance = 1e16).sigma > kshell_ce(:Zn; irradiance = 1e16).sigma
    @test_throws ArgumentError kshell_line_energy(:Xx)
    @test kshell_ce(10.279; irradiance = 1e15).ce == kshell_ce(:Ge; irradiance = 1e15).ce

    # --- attenuation ------------------------------------------------------------------------
    @test mass_attenuation(:Be, 10.0) ≈ 0.6466                     # NIST nodes
    @test mass_attenuation(:Cu, 8.979) ≈ 278.4                     # on an edge: above it
    @test mass_attenuation(:Cu, 8.97) < 40                         # just below the K edge
    @test mass_attenuation(:Mo, 20.0 - 1e-9) ≈ 13.08 rtol = 1e-6
    # log-log interpolation: a power law between nodes
    E = 12.0
    @test mass_attenuation(:Al, E) ≈ 26.23 * (7.955 / 26.23)^(log(E / 10) / log(1.5)) rtol = 1e-12
    @test mass_attenuation(:kapton, 7.0) == mass_attenuation(:mylar, 7.0)
    @test_throws ArgumentError mass_attenuation(:Al, 50.0)
    # transmission: Beer–Lambert, multiplicative over layers, density override
    T1 = filter_transmission((:Be => 0.2,), 9.0)
    @test T1 ≈ exp(-mass_attenuation(:Be, 9.0) * 1.848 * 0.02)
    @test filter_transmission((:Be => 0.2, :Al => 0.01), 9.0) ≈ T1 * filter_transmission((:Al => 0.01,), 9.0)
    @test filter_transmission((:Be => 0.1, :Be => 0.1), 9.0) ≈ T1
    @test filter_transmission((:Al => (0.01, 2 * 2.699),), 9.0) ≈ filter_transmission((:Al => 0.02,), 9.0)
    @test filter_transmission((), 9.0) == 1

    # --- the photometric chain --------------------------------------------------------------
    # one-sided square foil (w = 1 mm) facing +x at x = −10, uniformly irradiated
    w, D = 1.0, 10.0
    t = Target([backlighter_foil((90.0, 0.0); standoff = D, width = w, n = 4, two_sided = false)])
    fe = elements(t, :backlighter)
    A = sum(t.elem_area[fe])
    @test A ≈ w^2
    Ptot, τ = 10.0, 2.0                                            # TW, ns
    P = Ptot .* t.elem_area ./ A
    dep = Briareus.Deposition(P, Float64[], [Ptot], Ptot)
    I = 1e14 * Ptot / A                                            # W/cm²
    ce = kshell_ce(:Ge; irradiance = I, prepulse = true)
    hνJ = 10.279 * 1.602176634e-16
    EL = 1e3 * Ptot * τ                                            # J
    Ex = ce.ce * EL
    # point projection: Σ_e E_e/(4π r_e²) (isotropic), against the element sum by hand
    L = 500.0
    s = backlighter_signal(t, dep, (90.0, 0.0); line = :Ge, pulse = τ, prepulse = true,
                           distance = L, object = nothing, pixel = 0.05)
    @test s.laser_energy ≈ EL && s.irradiance ≈ I && s.ce ≈ ce.ce && s.sigma ≈ ce.sigma
    @test s.xray_energy ≈ Ex && s.facing ≈ 1
    Pd = [L, 0, 0]
    byhand = sum(Ex * t.elem_area[e] / A / (4π * sum(abs2, Pd - t.elem_centroid[e])) for e in fe)
    @test s.fluence ≈ byhand rtol = 1e-12
    @test s.fluence ≈ Ex / (4π * (L + D)^2) rtol = 1e-5           # a point source at this distance
    @test s.photons ≈ s.fluence / hνJ && s.photons_pixel ≈ s.photons * 0.05^2
    @test s.photons_lo ≈ s.photons * exp(-ce.sigma) && s.photons_hi ≈ s.photons * exp(ce.sigma)
    @test s.per_sr ≈ Ex / 4π                                       # isotropic, along the normal
    # Lambertian, along the normal: f = 4
    sl = backlighter_signal(t, dep, (90.0, 0.0); line = :Ge, pulse = τ, prepulse = true,
                            distance = L, object = nothing, angular = :lambertian)
    @test sl.per_sr ≈ Ex / π
    # energy conservation over the front hemisphere: ∫ per_sr dΩ = Ex/2 (isotropic) and Ex
    # (Lambertian). per_sr(θ) on a ring of detectors at polar angle θ from the normal,
    # integrated in μ = cos θ by the midpoint rule, exact here (f is constant or linear in μ,
    # and the same for every element of a flat foil)
    for (ang, total) in ((:isotropic, Ex / 2), (:lambertian, Ex))
        nμ = 200
        acc = 0.0
        for j in 1:nμ
            μ = (j - 0.5) / nμ
            dir = (μ, sqrt(1 - μ^2), 0.0)
            r = backlighter_signal(t, dep, dir; line = :Ge, pulse = τ, prepulse = true, distance = L,
                                   object = nothing, angular = ang)
            acc += r.per_sr * 2π / nμ
        end
        @test acc ≈ total rtol = 1e-10
    end
    # with a CE given, the chain is CE·E_laser
    s5 = backlighter_signal(t, dep, (90.0, 0.0); line = :Ge, pulse = τ, distance = L,
                            object = nothing, ce = 0.01)
    @test s5.xray_energy ≈ 0.01 * EL && s5.sigma == 0
    # filters scale the fluence by their transmission
    sf = backlighter_signal(t, dep, (90.0, 0.0); line = :Ge, pulse = τ, prepulse = true,
                            distance = L, object = nothing, filters = (:Be => 0.25, :Al => 0.05))
    @test sf.fluence ≈ s.fluence * filter_transmission((:Be => 0.25, :Al => 0.05), 10.279)

    # pinhole imaging of a Lambertian foil: uniform radiance R = Ex/(π A), each image
    # point receives R·A_p/q² at any obliquity (the cos θ of the emitter cancels)
    p, dpin = 100.0, 0.02
    q = L - p
    sp = backlighter_signal(t, dep, (90.0, 0.0); line = :Ge, pulse = τ, prepulse = true, distance = L,
                            object = nothing, radius = 0.4, pinhole = dpin, pinhole_distance = p,
                            angular = :lambertian)
    Ap = π * (dpin / 2)^2
    @test sp.fluence ≈ Ex / (π * A) * Ap / q^2 rtol = 1e-12
    @test all(x -> x ≈ sp.photons, sp.photons_range)
    @test sp.solid_angle ≈ Ap / (p + D)^2
    # energy through the pinhole = fluence × the foil's image area (magnification q/(p+D))
    @test sp.fluence * A * (q / (p + D))^2 ≈ sl.per_sr * Ap / (p + D)^2 rtol = 1e-12
    # a disk larger than the foil's image: the uncovered points see nothing
    so = backlighter_signal(t, dep, (90.0, 0.0); line = :Ge, pulse = τ, prepulse = true, distance = L,
                            object = nothing, radius = 2.0, pinhole = dpin, pinhole_distance = p,
                            angular = :lambertian)
    @test so.photons_range[1] == 0 && so.fluence < sp.fluence

    # the object (a sphere at TCC) does not attenuate; the same sphere declared hardware
    # (object = nothing) shadows the foil's center
    tb = Target([backlighter_foil((90.0, 0.0); standoff = D, width = w, n = 4, two_sided = false),
                 Component(:capsule, icosphere_mesh(0.9; n = 4), ablator())])
    Pb = zeros(nelements(tb)); fb = elements(tb, :backlighter)
    Pb[fb] .= Ptot .* tb.elem_area[fb] ./ A
    db = Briareus.Deposition(Pb, Float64[], [Ptot], Ptot)
    sb = backlighter_signal(tb, db, (90.0, 0.0); line = :Ge, pulse = τ, prepulse = true, distance = L)
    @test sb.fluence ≈ s.fluence rtol = 1e-12
    sh = backlighter_signal(tb, db, (90.0, 0.0); line = :Ge, pulse = τ, prepulse = true, distance = L,
                            object = nothing)
    @test sh.fluence < s.fluence                              # the sphere now shadows the center
    @test_throws ArgumentError backlighter_signal(t, dep, (90.0, 0.0); line = :Ge, pulse = τ, distance = L,
                                                  object = nothing, pinhole = 0.01)
end
