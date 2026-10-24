using Briareus, ForwardDiff, Test

@testset "two-zone CBET" begin
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.0)))
         for b in nif_quads()]
    inner = (:inner23, :inner30)
    out = CBET(leh_z = 5.55, mach = 2.0, multiplier = 0.1, clamp = nothing)
    ins = CBET(leh_z = 4.55, mach = 0.0, clamp = nothing)
    two = TwoZoneCBET(out, ins)
    r(m, dl) = receiver_ratio(m, q, dl; receivers = inner)
    # no flow and no detuning: no transfer
    @test r(ins, 0.0) ≈ 1 atol = 1e-12
    # supersonic outflow alone moves power from the inner to the outer cones; inflow
    # (negative Mach) the other way
    @test r(out, 0.0) < 0.98
    @test r(CBET(leh_z = 5.55, mach = -2.0, multiplier = 0.1, clamp = nothing), 0.0) > 1.02
    # inside, the red-detuned inner cones gain
    @test r(ins, 0.5) > 1.02
    # two zones: Δλ = 0 is not zero transfer (net inner → outer), detuning turns it around
    @test r(two, 0.0) ≈ r(out, 0.0) rtol = 1e-12
    @test r(two, 0.0) < 0.98
    @test r(two, 1.0) > 1.02
    # the zones apply in the order the beams meet them; power is conserved
    shifted = [remake(b; dlambda = b.group in inner ? 0.6 : 0.0) for b in q]
    eff = apply_lpi(two, shifted, 5.0)
    @test [b.power for b in eff] == [b.power for b in apply_lpi(ins, apply_lpi(out, shifted, 5.0), 5.0)]
    @test sum(b.power for b in eff) ≈ sum(b.power for b in shifted) rtol = 1e-12
    kz = TwoZoneCBET(leh_z = 5.05, offset = 0.5, inside = (multiplier = 0.5,))
    @test kz.outside.leh_z ≈ 5.55 && kz.inside.leh_z ≈ 4.55
    @test kz.outside.mach == 2.0 && kz.inside.mach == 0.0 && kz.inside.multiplier == 0.5
    @test TwoZoneCBET(leh_z = t -> 5.0 + t, offset = 0.5).inside.leh_z(1.0) ≈ 5.5
    # forward-mode gradients through both zones, including a plasma-state table
    function obj(x)
        m = TwoZoneCBET(CBET(leh_z = 5.55, mach = x[1], multiplier = x[2],
                             Te = ([0.0, 10.0], [x[3], x[3] + 2000])),
                        CBET(leh_z = 4.55, mach = x[4], multiplier = x[5]))
        return sum(b.power^2 for b in apply_lpi(m, shifted, 5.0) if b.group in inner)
    end
    x = [1.8, 0.1, 4000.0, 0.2, 0.9]
    g = ForwardDiff.gradient(obj, x)
    hs = [1e-6, 1e-6, 1e-2, 1e-6, 1e-6]
    fd = [(obj(x .+ hs[i] .* (1:5 .== i)) - obj(x .- hs[i] .* (1:5 .== i))) / 2hs[i] for i in 1:5]
    @test g ≈ fd rtol = 1e-5
    @test all(abs.(g) .> 0)
    # multi-epoch calibration recovers a known outside multiplier and a shared Te history
    # from synthetic cone fractions
    truth = TwoZoneCBET(leh_z = 5.05, outside = (multiplier = 0.1, clamp = nothing,
                                                 Te = ([0.0, 10.0], [3000.0, 6000.0])),
                        inside = (clamp = nothing, Te = ([0.0, 10.0], [3000.0, 6000.0])))
    beams_at(t) = [remake(b; power = b.power * (1 + 0.1t)) for b in shifted]
    frac(m, t) = (e = apply_lpi(m, beams_at(t), t);
                  sum(b.power for b in e if b.group in inner) / sum(b.power for b in e))
    data = [(t, frac(truth, t)) for t in 0.0:2.0:10.0]
    start = TwoZoneCBET(leh_z = 5.05, outside = (multiplier = 0.2, clamp = nothing,
                                                 Te = ([0.0, 10.0], [4000.0, 4000.0])),
                        inside = (clamp = nothing, Te = ([0.0, 10.0], [4000.0, 4000.0])))
    fitted, res = calibrate_cbet(start, beams_at, data; receivers = inner,
                                 fit = ((:outside, :multiplier), (:both, :Te)))
    @test fitted.outside.multiplier ≈ 0.1 rtol = 1e-3
    @test fitted.inside.Te[2] ≈ [3000.0, 6000.0] rtol = 1e-3
    @test fitted.outside.Te[2] == fitted.inside.Te[2]
    @test fitted.inside.multiplier == 1.0                     # not fitted
    @test_throws ArgumentError calibrate_cbet(start, beams_at, data; receivers = inner,
                                              fit = ((:outside, :length),))
    # joint multi-shot calibration: two shots with different powers, detunings and
    # timing (the second 0.5 ns early) recover one model; the penalty enters the loss
    isin2(b) = b.group in inner
    shot2(t) = [remake(b; power = b.power * (1 + 0.15t), dlambda = isin2(b) ? 0.4 : 0.0) for b in shifted]
    frac2(m, t, d) = (e = apply_lpi(m, shot2(t), t - d);
                      sum(b.power for b in e if b.group in inner) / sum(b.power for b in e))
    shots = [(; beams_at, data), (; beams_at = shot2, data = [(t, frac2(truth, t, -0.5)) for t in 0.0:2.0:10.0], delay = -0.5)]
    jfit, jres = calibrate_cbet(start, shots; receivers = inner, fit = ((:outside, :multiplier), (:both, :Te)))
    @test jfit.outside.multiplier ≈ 0.1 rtol = 1e-3
    @test jfit.inside.Te[2] ≈ [3000.0, 6000.0] rtol = 1e-3
    loss = Briareus._joint_loss(shots, inner, m -> 0.25)
    @test loss(truth) ≈ 0.25 atol = 1e-12
    @test Briareus._joint_loss(shots, inner, nothing)(start) > 1e-6
end

@testset "CBET saturation, clamp schedule and Δλ null" begin
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.0)),
                dlambda = b.group in (:inner23, :inner30) ? 0.5 : 0.0) for b in nif_quads()]
    inner = (:inner23, :inner30)
    P(m) = sum(b.power for b in apply_lpi(m, q, 0.0) if b.group in inner)
    lin = P(CBET(leh_z = 5.05, clamp = nothing))
    # hard minimum: no effect while δn stays below the clamp (bit for bit linear theory up
    # to summation order), and never less transfer than the smooth factor
    @test P(CBET(leh_z = 5.05, clamp = 1.0, saturation = :min)) ≈ lin rtol = 1e-12
    @test P(CBET(leh_z = 5.05, clamp = 1.0)) != lin
    P0 = sum(b.power for b in q if b.group in inner)
    for c in (1e-4, 3e-4, 1e-3)
        hard, smooth = P(CBET(leh_z = 5.05, clamp = c, saturation = :min)), P(CBET(leh_z = 5.05, clamp = c))
        @test 0 < abs(smooth - P0) < abs(hard - P0) < abs(lin - P0)
    end
    @test_throws ArgumentError P(CBET(leh_z = 5.05, saturation = :hard))
    # forward-mode gradients through the hard clamp
    g = ForwardDiff.derivative(c -> sum(b.power for b in apply_lpi(CBET(leh_z = 5.05, clamp = c, saturation = :min), q, 0.0)
                                        if b.group in inner), 3e-4)
    h = 1e-9
    fd = (P(CBET(leh_z = 5.05, clamp = 3e-4 + h, saturation = :min)) -
          P(CBET(leh_z = 5.05, clamp = 3e-4 - h, saturation = :min))) / 2h
    @test g ≈ fd rtol = 1e-4
    # the pulse-driven clamp schedule: the rise's knee is at 4.1 ns on N221204's
    # delivered pulse and at 4.0 ns on PRE 109 Fig. 9's (degraded) inputs, where the HYDRA
    # deck switches 0.002 → 0.0075 (3.9–4.1 ns)
    p = n221204_pulse()
    @test pulse_rise((p.t, p.power)) ≈ 4.11 atol = 0.03
    f9 = n221204_fig9()
    @test pulse_rise((f9.t, (f9.p44_in .+ f9.p50_in) ./ (1 .- f9.cf_in))) ≈ 4.0 atol = 0.06
    @test pulse_rise(PulseShape([0.0, 1.0, 2.0], [0.0, 0.0, 10.0])) ≈ 1.0      # the knee
    @test pulse_rise(PulseShape([0.0, 1.0, 2.0], [0.0, 0.0, 10.0])) ≈
          pulse_rise(PulseShape([0.0, 1.0, 2.0], [0.0, 0.0, 3.0])) rtol = 1e-12    # scale-free
    sch = cbet_clamp_schedule((p.t, p.power))
    @test sch[2] == [0.002, 0.0075]
    @test sch[1] ≈ [4.01, 4.21] atol = 0.03
    m = CBET(leh_z = 5.05, clamp = sch)
    @test Briareus._tvalue(m.clamp, 1.0) == 0.002 && Briareus._tvalue(m.clamp, 8.0) == 0.0075
    # the Δλ null: where the two-zone model leaves the inner power unchanged
    two = TwoZoneCBET(leh_z = 5.05, outside = (multiplier = 0.1, clamp = nothing), inside = (clamp = nothing,))
    d0 = cbet_null(two, q, 0.0; receivers = inner)
    @test 0 < d0 < 1.5
    @test receiver_ratio(two, q, d0; receivers = inner) ≈ 1 atol = 1e-4
    @test isnan(cbet_null(two, q, 0.0; receivers = inner, lo = 0.0, hi = d0 / 2))
end

@testset "N221204 CBET calibration (PRE 109 Fig. 9)" begin
    fig9 = n221204_fig9()
    @test fig9.t[1] ≈ 0.7 && fig9.t[end] ≈ 9.0 && length(fig9.t) == 84
    k = argmin(abs.(fig9.t .- 5.6))
    @test 0.33 < fig9.cf_in[k] < 0.34 && 0.63 < fig9.cf_post[k] < 0.67
    # Fig. 9a and 9b agree: per-cone powers give the plotted cone fraction at peak
    cfa = (fig9.p23_post[k] + fig9.p30_post[k]) /
          (fig9.p23_post[k] + fig9.p30_post[k] + fig9.p44_post[k] + fig9.p50_post[k])
    @test cfa ≈ fig9.cf_post[k] atol = 0.02
    base = [remake(b; dlambda = b.group in (:inner23, :inner30) ? 2.75 / 3 : 0.0)
            for b in hybride_mock_pointing(nif_beams())]
    m = n221204_cbet()
    isin(b) = b.group in (:inner23, :inner30)
    cfr(x) = sum(b.power for b in x if isin(b)) / sum(b.power for b in x)
    # the calibrated model reproduces Fig. 9b over 0.8–8.4 ns within the fit's stated
    # error (rms 0.016, worst 0.041 at 8.4 ns; digitization ±0.01)
    ts = 0.8:0.2:8.4
    res = map(ts) do t
        bs = n221204_fig9_beams(base, fig9, t)
        i = argmin(abs.(fig9.t .- t))
        @test cfr(bs) ≈ fig9.cf_in[i] rtol = 1e-12
        return cfr(apply_lpi(m, bs, t)) - fig9.cf_post[i]
    end
    @test sqrt(sum(abs2, res) / length(res)) < 0.02
    @test maximum(abs, res) < 0.045
    @test all(abs(res[findfirst(≈(t), ts)]) < 0.03 for t in (2.0, 3.0, 5.6, 7.0, 8.0))
    # not fitted: the 44° cone gives up more than the 50° cone at peak, as in Fig. 9a
    bs = n221204_fig9_beams(base, fig9, fig9.t[k])
    e = apply_lpi(m, bs, fig9.t[k])
    ratio(g) = sum(b.power for b in e if b.group === g) / sum(b.power for b in bs if b.group === g)
    @test ratio(:outer44) ≈ fig9.p44_post[k] / fig9.p44_in[k] atol = 0.05
    @test ratio(:outer50) ≈ fig9.p50_post[k] / fig9.p50_in[k] atol = 0.05
    @test ratio(:outer44) < ratio(:outer50)
    # Δλ = 0 moves power the other way (the outside zone's flow)
    cf0 = cfr(apply_lpi(m, [remake(b; dlambda = 0.0) for b in bs], fig9.t[k]))
    @test cf0 < fig9.cf_in[k] - 0.1
end

@testset "HYBRID-E multi-shot CBET data" begin
    a, b, f1, k22 = n221204_fig9(), n220919_fig9(), pre109_fig1(), n210808_fig1()
    at(tab, c, t) = getfield(tab, c)[argmin(abs.(tab.t .- t))]
    # N220919 (Fig. 9, black): same peak input cone fraction as N221204, slightly less
    # transfer at 2.5 Å; its curves end at 7.4 ns
    @test maximum(b.t) ≈ 7.4 && all(isnan, b.cf_post[b.t .> 7.45])
    @test at(b, :cf_in, 5.6) ≈ at(a, :cf_in, 5.6) atol = 0.01
    @test 0 < at(a, :cf_post, 5.6) - at(b, :cf_post, 5.6) < 0.03
    @test at(b, :cf_in, 3.4) < at(a, :cf_in, 3.4) - 0.02          # lower rise cone fraction
    @test all(x -> x in ("v", "p", "m", "i", "c", "h", ""), b.f_cf_post)
    # PRE 109 Fig. 1e reproduces Fig. 9b's N221204 curves (two rasters) to 0.01
    for t in (2.0, 3.4, 5.6, 7.0)
        @test at(f1, :cf_post_n221204, t) ≈ at(a, :cf_post, t) atol = 0.01
        @test at(f1, :cf_in_n210808, t) ≈ at(k22, :cf_in, t) atol = 0.01   # same incident pulse
    end
    # N210808 at 1.8 Å: less transfer than N221204 at 2.75 Å, but only ~0.05 at peak start
    d = at(a, :cf_post, 5.6) - at(f1, :cf_post_n210808, 5.6)
    @test 0.02 < d < 0.08
    # the two papers' N210808 post-CBET curves differ (different HYDRA runs)
    @test at(k22, :cf_post, 3.4) - at(f1, :cf_post_n210808, 3.4) > 0.05
    # PRE 106 Fig. 1d per-cone powers reproduce its Fig. 1e cone fraction at peak
    k = argmin(abs.(k22.t .- 5.6))
    cfa = (k22.p23_post[k] + k22.p30_post[k]) / (k22.p23_post[k] + k22.p30_post[k] + k22.p44_post[k] + k22.p50_post[k])
    @test cfa ≈ k22.cf_post[k] atol = 0.02
    s = hybride_shots()
    @test s.shot == ["N210808", "N220919", "N221204", "N230729"]
    @test s.dlambda == [1.8, 2.5, 2.75, 2.9]
    @test s.hs_p2[2] ≈ -10.8 && isnan(s.u_hs_p2[4]) && s.clamp_peak[3] == 0.0075
end

@testset "PRE 109 Fig. 10 capsule flux P2/P0" begin
    f = n221204_fig10()
    @test keys(f) == (:t, :n210808, :u210808, :n220919, :u220919, :n221204, :u221204)
    @test all(length(c) == length(f.t) for c in f) && length(f.t) > 80
    @test issorted(f.t; lt = <=) && f.t[1] ≈ 0.8 && f.t[end] ≈ 8.9
    @test all(0 < u <= 0.015 for u in (f.u210808..., f.u220919..., f.u221204...))
    # N221204's minimum: −0.17 near 5.5 ns, deeper than the other two shots'
    k = argmin(f.n221204)
    @test f.n221204[k] ≈ -0.173 atol = 0.005
    @test 5.45 < f.t[k] < 5.7
    @test minimum(f.n220919) ≈ -0.149 atol = 0.005
    @test minimum(f.n210808) > f.n221204[k] + 0.02
    # foot near zero, late positive swing
    foot = 1.0 .<= f.t .<= 3.0
    @test all(abs.(f.n221204[foot]) .< 0.02)
    @test f.n221204[findfirst(≈(8.0), f.t)] ≈ 0.095 atol = 0.01
    @test maximum(f.n220919) ≈ 0.153 atol = 0.005
end

@testset "joint HYBRID-E CBET calibration (hybride_cbet)" begin
    base = hybride_mock_pointing(nif_beams())
    isin(b) = b.group in (:inner23, :inner30)
    cfr(x) = sum(b.power for b in x if isin(b)) / sum(b.power for b in x)
    k22, f1 = n210808_fig1(), pre109_fig1()
    at(c, t) = (i = argmin(abs.(k22.t .- t)); abs(k22.t[i] - t) < 0.051 ? c[i] : NaN)
    t808 = (t = f1.t, cf_in = f1.cf_in_n210808, cf_post = f1.cf_post_n210808,
            p44_in = [at(k22.p44_in, t) for t in f1.t], p50_in = [at(k22.p50_in, t) for t in f1.t])
    total(tab) = (ok = isfinite.(tab.cf_in) .& isfinite.(tab.p44_in) .& isfinite.(tab.p50_in);
                  (tab.t[ok], ((tab.p44_in .+ tab.p50_in) ./ (1 .- tab.cf_in))[ok]))
    ref = n221204_fig9()
    m = hybride_cbet()
    function rmsres(tab, dλ)
        delay = pulse_rise(total(tab)) - pulse_rise(total(ref))
        r = Float64[]
        for t in 0.8:0.2:8.4
            i = argmin(abs.(tab.t .- t))
            (abs(tab.t[i] - t) < 0.051 && all(isfinite, (tab.cf_post[i], tab.cf_in[i], tab.p44_in[i], tab.p50_in[i]))) || continue
            bs = [remake(b; dlambda = isin(b) ? dλ / 3 : 0.0) for b in n221204_fig9_beams(base, tab, t)]
            push!(r, cfr(apply_lpi(m, bs, t - delay)) - tab.cf_post[i])
        end
        return sqrt(sum(abs2, r) / length(r))
    end
    # the joint fit's stated cone-fraction rms per shot (0.027 / 0.016 / 0.047)
    @test rmsres(ref, 2.75) < 0.03
    @test rmsres(n220919_fig9(), 2.5) < 0.02
    @test rmsres(t808, 1.8) < 0.05
    # the Δλ null at peak inside (or at the edge of) NC24's 1–1.55 Å
    bs = [remake(b; dlambda = 0.0) for b in n221204_fig9_beams(base, ref, 5.6)]
    @test 0.85 < 3cbet_null(m, bs, 5.6; receivers = (:inner23, :inner30)) < 1.1
    # hard clamp, binding at peak; delay shifts the time tables
    @test m.inside.saturation === :min && m.inside.clamp[2][2] < 1e-3
    d = hybride_cbet(delay = -0.2)
    @test d.inside.Te[1] ≈ m.inside.Te[1] .- 0.2 && d.outside.clamp[1] ≈ m.outside.clamp[1] .- 0.2
end
