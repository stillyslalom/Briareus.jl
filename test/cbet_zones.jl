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
