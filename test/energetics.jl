# Absolute energetics: the wall-material factor of `WallLossAlbedo` (Au-lined DU against
# Döppner 2015's 6.5 % peak-power equivalence), time- and intensity-dependent laser →
# x-ray conversion (closed-form Marshak balance, Dewald 2008's table, ForwardDiff), the
# capsule's effective-radius albedo and the time-integrated `energy_budget` of a
# `drive_history`. Research notes, "Wall material, conversion and capsule energetics
# (2026-10-24)".

using Briareus, Test, ForwardDiff, LinearAlgebra

# Unit laser source on a hohlraum's wall in two pairs of bands (fraction, |z| range).
function _band_source(t, R, bands)
    P = zeros(nelements(t))
    for e in elements(t, :wall)
        pts = [p for k in t.elem_tris[e] for p in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k],
                                                    t.tri_v0[k] + t.tri_e2[k])]
        z0, z1 = extrema(p[3] for p in pts)
        for (f, a, b) in bands, (lo, hi) in ((a, b), (-b, -a))
            P[e] += f / 2 * max(0.0, min(z1, hi) - max(z0, lo)) / (b - a) *
                    t.elem_area[e] / (2π * R * (z1 - z0))
        end
    end
    return P ./ sum(P)
end
_capsule_Tr(sol, t) = (c = elements(t, :capsule);
                       radiation_temperature(sum(sol.H[c] .* t.elem_area[c]) / sum(t.elem_area[c])))

@testset "wall material factor" begin
    g = WallLossAlbedo()
    # the defaults are gold's, bit for bit, and the factor scales gold's clamped loss
    for (τ, T) in ((0.05, 90.0), (1.0, 150.0), (8.0, 313.0), (14.0, 300.0))
        loss = g.coefficient * (T / 100)^g.temperature_exponent * max(τ, g.t0)^g.time_exponent
        @test albedo(g, τ, T) === clamp(1 - loss, g.range...)
        @test wall_loss_factor(g, τ, T) ≈ 1
        d = WallLossAlbedo(factor = 0.8, factor_exponent = -0.2)
        @test 1 - albedo(d, τ, T) ≈ 0.8 * (T / 300)^-0.2 * (1 - albedo(g, τ, T))
    end
    # the cap applies to gold's loss: a DU wall may exceed it
    @test albedo(g, 14.0, 290.0) < 0.95 < albedo(du_wall_loss(), 14.0, 290.0)
    # liner: gold while the front is inside it, the substrate's ratio once far past it
    d = du_wall_loss()
    f = d.front
    mF(τ, T) = f.coefficient * (T / 100)^f.temperature_exponent * τ^f.time_exponent
    @test wall_loss_factor(d, 0.2, 150.0) == 1                 # m_F < liner
    w = d.liner / mF(8.0, 300.0)
    @test 0 < w < 0.1
    @test wall_loss_factor(d, 8.0, 300.0) ≈ w + (1 - w) * Briareus.DU_LOSS_FACTOR
    nl = du_wall_loss(liner = 0.0)
    @test wall_loss_factor(nl, 8.0, 250.0) ≈ Briareus.DU_LOSS_FACTOR * (250 / 300)^-0.2
    @test uranium().name === :DU && uranium().albedo isa WallLossAlbedo

    # Döppner 2015: an Au-lined DU wall is worth 6.5 % of peak laser power on the NIF
    # high-foot hohlraum (5.75 × 9.43 mm, 3.1 mm LEHs) at peak, 13 ns; the calibration
    # mesh was h = 0.4, here a coarser one (the ratio moves by < 0.3 %)
    R, h = 2.875, 0.5
    mk(alb) = hohlraum(; length = 9.43, radius = R, leh_radius = 1.55, capsule_radius = 1.125, h,
                       capsule_h = h / 3, wall = Material(:w; albedo = alb, conversion = 0.9),
                       capsule = ablator(albedo = 0.3, conversion = 0.0))
    ta, td = mk(WallLossAlbedo()), mk(du_wall_loss())
    ma = RadiationModel(ta)
    md = RadiationModel(td; F = ma.F)
    s = _band_source(ta, R, ((0.30, 0.5, 1.7), (0.70, 3.0, 4.2)))
    P = 0.86 * 390
    Tr(m, p) = _capsule_Tr(solve(m, p .* s; time = 13.0), ta)
    T0 = Tr(md, P)
    lo, hi = P, 1.2P
    for _ in 1:30
        mid = (lo + hi) / 2
        Tr(ma, mid) < T0 ? (lo = mid) : (hi = mid)
    end
    @test lo / P ≈ 1.065 atol = 0.006
    # Döppner's T-1 pair at equal power: +5 ± 6 eV measured
    @test 3.5 < T0 - Tr(ma, P) < 6.5
end

@testset "conversion models" begin
    c = Conversion(efficiency = 0.7, time_exponent = 0.12)
    @test conversion_efficiency(c, 1.0) ≈ 0.7
    @test conversion_efficiency(c, 8.0) ≈ 0.7 * 8^0.12
    @test conversion_efficiency(c, 0.0) ≈ 0.7 * 0.1^0.12                   # t clamped at t0
    @test conversion_efficiency(Conversion(efficiency = 0.7, time_exponent = 0.12, start = 1.0), 9.0) ≈
          conversion_efficiency(c, 8.0)
    n = nif_conversion()
    @test conversion_efficiency(n, 5.0) ≈ 0.85 atol = 0.005
    @test conversion_efficiency(n, 8.0) ≈ 0.90 atol = 0.005
    @test conversion_efficiency(n, 30.0) == 0.95
    d = Conversion(intensities = (1e14, 1e15), efficiencies = (0.95, 0.80))      # Dewald 2008
    @test conversion_efficiency(d, 1.0, 1e14) ≈ 0.95
    @test conversion_efficiency(d, 1.0, 1e15) ≈ 0.80
    @test conversion_efficiency(d, 1.0, sqrt(10) * 1e14) ≈ 0.875
    @test conversion_efficiency(d, 1.0, 1e13) ≈ 0.95 && conversion_efficiency(d, 1.0, 1e17) ≈ 0.80
    @test_throws ArgumentError conversion_efficiency(d, 1.0)
    @test_throws ArgumentError Conversion(intensities = (1e14, 1e15))
    @test Material(:x; conversion = 1).conversion === 1.0

    t0 = hohlraum(; h = 0.8, capsule_h = 0.4)
    mats(η) = hohlraum(; h = 0.8, capsule_h = 0.4, wall = gold(albedo = WallLossAlbedo(), conversion = η))
    P = zeros(nelements(t0))
    w = elements(t0, :wall)
    P[w] .= range(0.0, 2.0, length(w))
    # a constant Conversion reproduces the number, bit for bit
    s1 = solve(RadiationModel(mats(0.75)), P; time = 3.0)
    m2 = RadiationModel(mats(Conversion(efficiency = 0.75)))
    s2 = solve(m2, P; time = 3.0)
    @test s1.source == s2.source && s1.H == s2.H
    @test conversions(m2.target, 3.0) == conversions(s1.model.target)

    # time dependence in solve and flux_response
    mt = RadiationModel(mats(c); F = s1.model.F)
    st = solve(mt, P; time = 4.0)
    @test st.source[w] ≈ 0.7 * 4^0.12 .* P[w] ./ t0.elem_area[w]
    obs = elements(t0, :capsule)
    R = flux_response(mt, obs; time = 4.0, albedo = st.albedo)
    @test R * P ≈ st.H[obs]

    # intensity dependence: per-element source, the secant response, and duals
    mi = RadiationModel(mats(d); F = s1.model.F)
    si = solve(mi, P; time = 4.0)
    I = to_Wcm2.(P[w] ./ t0.elem_area[w])
    @test si.source[w] ≈ [conversion_efficiency(d, 4.0, x) for x in I] .* P[w] ./ t0.elem_area[w]
    @test_throws ArgumentError flux_response(mi, obs; time = 4.0)
    Ri = flux_response(mi, obs; time = 4.0, albedo = si.albedo, laser = P)
    @test Ri * P ≈ si.H[obs]
    α = si.albedo
    f(x) = sum(solve(mi, x .* P; time = 4.0, albedo = α).H[obs])
    g = ForwardDiff.derivative(f, 1.0)
    @test g ≈ (f(1 + 1e-6) - f(1 - 1e-6)) / 2e-6 rtol = 1e-6
    @test !(g ≈ f(1.0))                             # the intensity dependence is in the derivative

    # closed form: with wall loss ∝ T^-0.7 t^-0.38 dominating, η ∝ t^s gives T ∝
    # t^((0.38 + s)/3.3) (Lindl 2004, Fig. 4-4: 0.7 t^0.12 makes T ∝ t^0.15); a scale-1
    # Nova hohlraum with pinhole LEHs and a reflecting capsule, 20 TW on the wall
    for (s, k) in ((0.0, 0.38 / 3.3), (0.12, 0.50 / 3.3))
        t = hohlraum(; length = 2.55, radius = 0.8, leh_radius = 0.05, capsule_radius = 0.25,
                     h = 0.15, capsule_h = 0.08, capsule = ablator(albedo = 0.999),
                     wall = Material(:Au; albedo = WallLossAlbedo(),
                                     conversion = Conversion(efficiency = 0.7, time_exponent = s)))
        m = RadiationModel(t)
        Q = zeros(nelements(t))
        ww = elements(t, :wall)
        Q[ww] .= 20.0 .* t.elem_area[ww] ./ sum(t.elem_area[ww])
        T3, T10 = (_capsule_Tr(solve(m, Q; time = τ), t) for τ in (0.3, 1.0))
        @test log(T10 / T3) / log(1 / 0.3) ≈ k atol = 0.008
    end
end

@testset "capsule albedo and energy budget" begin
    ca = CapsuleAlbedo(albedo = 0.2, radius = ([0.0, 5.0, 9.0], [1.0, 1.0, 0.5]))
    @test albedo(ca, 2.0) ≈ 0.2
    @test albedo(ca, 7.0) ≈ 1 - 0.8 * 0.75^2
    @test capsule_radius_fraction(ca, 20.0) == 0.5
    @test albedo(CapsuleAlbedo(radius = τ -> 0.5), 1.0) ≈ 1 - 0.8 / 4
    ta = TemperatureAlbedo([100.0, 300.0], [0.4, 0.2])
    @test albedo(ta, 0.0, 200.0) ≈ 0.3 && albedo(ta, 0.0) ≈ 0.25
    @test albedo(CapsuleAlbedo(albedo = ta), 0.0, 300.0) ≈ 0.2

    t = hohlraum(; h = 0.8, capsule_h = 0.4, wall = gold(albedo = WallLossAlbedo()),
                 capsule = Material(:ablator; albedo = CapsuleAlbedo(albedo = ta, radius = ca.radius)))
    m = RadiationModel(t)
    P = zeros(nelements(t))
    w = elements(t, :wall)
    P[w] .= 300.0 / length(w)
    sol = solve(m, P; time = 7.0)
    c = elements(t, :capsule)
    # E_cap: the smaller sphere's absorption at its own (iterated) albedo
    αc = [albedo(ta, 7.0, radiation_temperature(sol.H[e])) for e in c]
    @test energy_balance(sol).absorbed[:capsule] ≈
          sum((1 .- αc) .* 0.75^2 .* sol.H[c] .* t.elem_area[c]) rtol = 1e-3

    # drive_history records the budget; energy_budget closes it
    beams = nif_beams()
    op = DepositionOperator(t, beams)
    pulse = PulseShape([0.0, 1.0, 6.0, 7.0], [0.0, 1.5, 1.5, 0.0])
    lpi = Backscatter(:inner23 => 0.1, :inner30 => ([0.0, 7.0], [0.0, 0.2]))
    h = drive_history(m, op, τ -> with_pulses(beams, pulse, τ), 0.0:0.5:7.0; lpi)
    b = energy_budget(h)
    @test b.setup ≈ 192 * pulse_energy(pulse)
    @test 0 < b.backscattered < 0.1 * b.setup
    @test b.coupled ≈ b.setup - b.backscattered
    total = b.backscattered + b.unabsorbed + sum(values(b.heating)) + sum(values(b.escaped)) +
            sum(values(b.transmitted); init = 0.0) + b.lost
    @test total ≈ b.setup rtol = 1e-6
    @test b.absorbed[:capsule] > 0 && haskey(b.escaped, :leh_top)
    @test b.heating[:wall] ≈ b.laser[:wall] - b.source[:wall] + b.absorbed[:wall]
    @test energy_budget(h; from = 1.0, to = 6.0).setup ≈ 192 * 1.5 * 5
end

@testset "wall-loss history" begin
    t = hohlraum(; h = 0.8, capsule_h = 0.4, wall = gold(albedo = WallLossAlbedo()))
    m = RadiationModel(t)
    beams = nif_beams()
    op = DepositionOperator(t, beams)
    # a foot at 1/8 of the peak power, then a step: the wall heated in the foot holds
    # little energy at the peak temperature
    pulse = PulseShape([0.0, 0.2, 4.0, 4.2, 9.0], [0.0, 0.2, 0.2, 1.6, 1.6])
    times = [0.5:0.5:4.0; 4.25; 4.5:0.5:9.0]
    rec = sol -> [radiation_temperature(x) for x in sol.H]
    on = drive_history(m, op, τ -> with_pulses(beams, pulse, τ), times; keep_albedo = true, record = rec)
    off = drive_history(m, op, τ -> with_pulses(beams, pulse, τ), times; keep_albedo = true, record = rec,
                        wall_history = false)
    wm = WallLossAlbedo()
    w = elements(t, :wall)
    s = 3.3 / 0.62
    # each wall element's albedo is the constant-T loss at its effective time
    # t_eff = ∫ (T/1 heV)^s dt / (T/1 heV)^s (trapezoidal over the steps; constant T
    # before the first step)
    err_on = err_off = 0.0
    for e in w
        T = [r[e] for r in on.record]
        I = times[1] * (T[1] / 100)^s
        for k in eachindex(times)
            k > 1 && (I += (times[k] - times[k-1]) * ((T[k-1] / 100)^s + (T[k] / 100)^s) / 2)
            err_on = max(err_on, abs(on.albedo[k][e] - albedo(wm, I / (T[k] / 100)^s, T[k])))
            err_off = max(err_off, abs(off.albedo[k][e] - albedo(wm, times[k], off.record[k][e])))
        end
    end
    @test err_on < 2e-3 && err_off < 2e-3
    # identical in the foot's first step; after the step the wall loses more, and the
    # drive then rises toward the constant-T history's
    @test on.Tr[1] ≈ off.Tr[1] rtol = 1e-8
    k1, k2 = findfirst(==(4.5), times), findfirst(==(9.0), times)
    @test on.wall_albedo[k1] < off.wall_albedo[k1] - 0.03
    @test on.Tr[k1] < off.Tr[k1] - 5
    @test on.Tr[k2] - on.Tr[k1] > off.Tr[k2] - off.Tr[k1] + 5
    @test abs(on.Tr[k2] - off.Tr[k2]) < abs(on.Tr[k1] - off.Tr[k1])
end

@testset "wall re-emission" begin
    t = hohlraum(; h = 0.8, capsule_h = 0.4, wall = gold(albedo = WallLossAlbedo()))
    m = RadiationModel(t)
    beams = nif_beams()
    op = DepositionOperator(t, beams)
    w = elements(t, :wall)
    η = conversions(t)
    # per step: each element's incident-flux temperature and net flux into it,
    # (1 − α) H − E with E the emission the history adds to the source
    rec = sol -> (T = radiation_temperature.(sol.H),
                  F = (1 .- sol.albedo) .* sol.H .- (sol.source .- η .* sol.laser ./ t.elem_area))
    # a rising drive: no element falls below its running maximum, bit for bit
    rise = PulseShape([0.0, 0.2, 4.0, 4.2, 9.0], [0.0, 0.2, 0.2, 1.6, 1.6])
    times = [0.5:0.5:4.0; 4.25; 4.5:0.5:9.0]
    on = drive_history(m, op, τ -> with_pulses(beams, rise, τ), times; keep_albedo = true)
    off = drive_history(m, op, τ -> with_pulses(beams, rise, τ), times; keep_albedo = true,
                        wall_release = false)
    @test on.albedo == off.albedo && on.Tr == off.Tr && on.moments == off.moments
    @test all(iszero, on.released[:wall]) && energy_budget(on).released[:wall] == 0

    # a 40 % droop after a flat peak
    drop = PulseShape([0.0, 0.2, 4.0, 4.2, 7.0, 7.5, 9.0], [0.0, 0.2, 0.2, 1.6, 1.6, 0.96, 0.96])
    times = [0.5:0.5:7.0; 7.1:0.1:9.0]
    on = drive_history(m, op, τ -> with_pulses(beams, drop, τ), times; record = rec, tol = 1e-10)
    off = drive_history(m, op, τ -> with_pulses(beams, drop, τ), times; wall_release = false, tol = 1e-10)
    k7, k75, k9 = (findfirst(≈(x), times) for x in (7.0, 7.5, 9.0))
    @test on.Tr[1:k7] == off.Tr[1:k7]
    # the wall gives energy back: the drive falls less, and never above the peak's
    @test on.Tr[k75] > off.Tr[k75] + 2
    # once the drive holds again the cooled wall takes energy back in: the two converge
    @test abs(on.Tr[k9] - off.Tr[k9]) < 0.25 * (on.Tr[k75] - off.Tr[k75])
    @test maximum(on.Tr[k7+1:end]) < on.Tr[k7] + 0.5
    @test on.wall_albedo[k75] > off.wall_albedo[k75]
    b = energy_budget(on)
    total = b.backscattered + b.unabsorbed + sum(values(b.heating)) + sum(values(b.escaped)) +
            sum(values(b.transmitted); init = 0.0) + b.lost
    @test total ≈ b.setup rtol = 1e-9
    @test 0 < b.released[:wall] < b.absorbed[:wall]
    @test b.absorbed[:wall] < energy_budget(off).absorbed[:wall]
    # each element's net flux is the heating loss at its running maximum less the
    # Duhamel response to its deficit, evaluated here from the recorded temperatures
    # (wall loss Lindl's: c = 0.32, p = −0.7, q = −0.38; a = 3.3, b = 0.62; release 0.6)
    wm = WallLossAlbedo()
    a_, b_ = 3.3, 0.62
    φ(T) = (T / 100)^a_
    r(T) = (T / 100)^(a_ / b_)
    err = 0.0
    for e in w
        T = [x.T[e] for x in on.record]
        Tm = accumulate(max, T)
        D = φ.(Tm) .- φ.(T)
        I = 0.0
        for k in eachindex(times)
            I = k == 1 ? times[1] * r(T[1]) : I + (times[k] - times[k-1]) * (r(Tm[k-1]) + r(Tm[k])) / 2
            te = I / r(Tm[k])
            L = (1 - albedo(wm, te, Tm[k])) * SIGMA_SB * Tm[k]^4
            C = wm.release * L * te^(1 - b_) / (b_ * φ(Tm[k]))
            G = sum((D[j] - D[j-1]) / (times[j] - times[j-1]) *
                    ((times[k] - times[j-1])^b_ - (times[k] - times[j])^b_) for j in 2:k; init = 0.0)
            F = L - C * G
            err = max(err, abs(on.record[k].F[e] - F) / (SIGMA_SB * T[k]^4))
        end
    end
    @test err < 1e-6
    # the release converges with the step: 0.1 against 0.05 ns steps after the droop
    fine = [0.5:0.5:7.0; 7.05:0.05:9.0]
    hf = drive_history(m, op, τ -> with_pulses(beams, drop, τ), fine)
    for x in (7.5, 8.0, 9.0)
        @test hf.Tr[findfirst(≈(x), fine)] ≈ on.Tr[findfirst(≈(x), times)] atol = 0.1
    end
end

@testset "N221204 delivered pulse and Dante 1" begin
    p = n221204_pulse()
    @test keys(p) == (:t, :power, :dante1)
    @test p.t[1] == 0 && p.t[end] ≈ 9.6 && all(diff(p.t) .≈ 0.05)
    @test pulse_energy(PulseShape(p.t, p.power)) ≈ 2050 rtol = 0.005      # 2.05 MJ
    @test maximum(p.power) ≈ 448 atol = 3
    @test all(>(425), p.power[5.6 .<= p.t .<= 7.9])
    ok = isfinite.(p.dante1)
    @test p.t[ok][1] ≈ 3.3 && p.t[ok][end] ≈ 9.4
    k = argmax(replace(p.dante1, NaN => -Inf))
    @test p.dante1[k] ≈ 311.9 atol = 0.1
    @test p.t[k] ≈ 8.3
end
