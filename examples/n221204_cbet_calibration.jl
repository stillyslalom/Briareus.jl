# Two-zone CBET calibrated to the N221204 HYDRA transfer of Kritcher et al., PRE 109,
# 025204 (2024) Fig. 9 (digitized: data/pre109_fig9_n221204.csv). Fits the outside
# zone's multiplier and outflow Mach number, the inside zone's multiplier and one LEH
# electron-temperature history shared by both zones to the post-transfer inner cone
# fraction (Fig. 9b) over 0.8–8.4 ns, then checks what was not fitted: the per-cone
# ratios of Fig. 9a and the response to Δλ. The result is `n221204_cbet()`; research
# notes, "Two-zone CBET calibrated to PRE 109 Fig. 9 (2026-10-23)". ~6 min on 8 threads.
#
#   julia --project=. -t auto examples/n221204_cbet_calibration.jl

using Briareus, Printf

fig9 = n221204_fig9()
# Mock HYBRID-E pointing, Δλ 2.75 Å (1ω) on the inner cones = 0.917 Å at 3ω
base = [remake(b; dlambda = b.group in (:inner23, :inner30) ? 2.75 / 3 : 0.0)
        for b in hybride_mock_pointing(nif_beams())]
beams_at(t; dλ = 2.75 / 3) = [remake(b; dlambda = b.group in (:inner23, :inner30) ? dλ : 0.0)
                              for b in n221204_fig9_beams(base, fig9, t)]
isin(b) = b.group in (:inner23, :inner30)
cf(x) = sum(b.power for b in x if isin(b)) / sum(b.power for b in x)
post(t) = fig9.cf_post[argmin(abs.(fig9.t .- t))]

# Start: the deck's clamp schedule in both zones, a guessed Te history
clamp = ([0.0, 3.9, 4.1, 20.0], [0.002, 0.002, 0.0075, 0.0075])
Te0 = ([1.0, 3.0, 4.0, 4.75, 5.5, 6.5, 8.0], [3000.0, 3000.0, 3500.0, 4200.0, 4800.0, 5200.0, 6300.0])
m0 = TwoZoneCBET(leh_z = 5.62, outside = (mach = 1.4, multiplier = 0.28, clamp, Te = Te0),
                 inside = (multiplier = 1.1, clamp, Te = Te0))
ts = collect(0.8:0.2:8.4)
fit = ((:outside, :multiplier), (:inside, :multiplier), (:outside, :mach), (:both, :Te))
t0 = time()
m, res = calibrate_cbet(m0, beams_at, [(t, post(t)) for t in ts]; receivers = (:inner23, :inner30), fit)
@printf("fit: %d iterations, %.0f s\n", Briareus.Optim.iterations(res), time() - t0)
@printf("outside: multiplier %.3f, Mach %.3f; inside: multiplier %.3f\n",
        m.outside.multiplier, m.outside.mach, m.inside.multiplier)
println("Te (eV) at ", m.inside.Te[1], " ns: ", round.(m.inside.Te[2]))

r = [cf(apply_lpi(m, beams_at(t), t)) - post(t) for t in ts]
@printf("cone fraction residual: rms %.3f, worst %+.3f at %.1f ns\n",
        sqrt(sum(abs2, r) / length(r)), r[argmax(abs.(r))], ts[argmax(abs.(r))])
for m′ in (m, n221204_cbet())
    q = [cf(apply_lpi(m′, beams_at(t), t)) - post(t) for t in ts]
    @printf("  %s: rms %.4f\n", m′ === m ? "this fit   " : "n221204_cbet", sqrt(sum(abs2, q) / length(q)))
end

println("\nper-cone post/input (model | PRE 109 Fig. 9a; not fitted)")
for t in (4.5, 5.6, 7.0, 8.0)
    bs = beams_at(t)
    e = apply_lpi(m, bs, t)
    i = argmin(abs.(fig9.t .- t))
    ratio(g) = sum(b.power for b in e if b.group === g) / sum(b.power for b in bs if b.group === g)
    @printf("%.1f ns  23° %.2f 30° %.2f 44° %.2f 50° %.2f | %.2f %.2f %.2f %.2f\n", t,
            ratio(:inner23), ratio(:inner30), ratio(:outer44), ratio(:outer50),
            fig9.p23_post[i] / fig9.p23_in[i], fig9.p30_post[i] / fig9.p30_in[i],
            fig9.p44_post[i] / fig9.p44_in[i], fig9.p50_post[i] / fig9.p50_in[i])
end

println("\npost-transfer cone fraction vs Δλ (1ω), Fig. 9 powers")
println("Δλ (Å)  5.6 ns  7.0 ns  8.0 ns  2.0 ns")
for dλ in (0.0, 0.5, 1.0, 1.4, 1.8, 2.25, 2.75)
    @printf("%.2f   %s\n", dλ, join([@sprintf("%.3f  ", cf(apply_lpi(m, beams_at(t; dλ = dλ / 3), t)))
                                     for t in (5.6, 7.0, 8.0, 2.0)]))
end
