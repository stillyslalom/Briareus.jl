# Joint multi-shot calibration of the two-zone CBET model to HYDRA's post-transfer cone
# fraction for three HYBRID-E shots at three detunings: N210808 (Δλ 1.8 Å, PRE 109,
# 025204 Fig. 1e with Kritcher et al. PRE 106, 025201 Fig. 1d's incident powers),
# N220919 (2.5 Å) and N221204 (2.75 Å) (PRE 109 Fig. 9). One plasma history (Te(t)), one pair of zone multipliers and one outflow Mach
# number serve all shots; per shot only what physically differs: its input powers, its
# Δλ, and its pulse timing (the model's time tables follow each shot's rise to peak).
# Physical constraints: Te ≤ `temax` (published LEH values), and the Δλ null inside NC24's
# 1–1.55 Å for HYBRID-E (a penalty on N221204's peak beams). PRE 109 Fig. 10 is not used.
# Research notes, "Joint multi-shot CBET calibration (2026-10-24)". 15–90 min per fit on
# 4 threads (`clamp`: 98 iterations; the result is `hybride_cbet()`).
#
#   julia --project=. -t 4 examples/hybride_cbet_calibration.jl [variant ...] [key=value ...]
#
# Variants: `free` (Te ≤ 20 keV, no null constraint, the deck's smooth clamp schedule),
# `tebound` (free with Te ≤ temax), `constrained` (Te ≤ temax and the null penalty, the
# deck's clamps, which do not bind), `clamp` (constrained, plus the foot and peak clamp
# values fitted with HYDRA's hard-minimum saturation on each shot's pulse-driven
# schedule, started where the clamp binds). Options: `temax=5500`, `iterations=150`, `n210808=pre109|pre106` (below),
# `laplace=on|off` (Gauss–Newton covariance at the optimum), `profile=on|off` (best loss
# with the null pinned at 0.5–2.0 Å). Writes `hybride_cbet_<variant>.txt`.

using Briareus, Printf, ForwardDiff, LinearAlgebra

opt = Dict(String.(split(a, '=', limit = 2)) for a in ARGS if occursin('=', a))
variants = [a for a in ARGS if !occursin('=', a)]
isempty(variants) && (variants = ["constrained", "clamp"])
temax = parse(Float64, get(opt, "temax", "5500"))
iterations = parse(Int, get(opt, "iterations", "150"))

const INNER = (:inner23, :inner30)
isin(b) = b.group in INNER
cf(x) = sum(b.power for b in x if isin(b)) / sum(b.power for b in x)
base = hybride_mock_pointing(nif_beams())

# Shots: HYDRA's input powers per cone (outer cones and the input cone fraction), Δλ at
# 1ω, and the pulse's rise (`pulse_rise` on the input total power, inner + outer).
function shot(name, tab, dλ)
    ok = [isfinite(tab.cf_in[i]) && isfinite(tab.p44_in[i]) && isfinite(tab.p50_in[i]) for i in eachindex(tab.t)]
    total = [(tab.p44_in[i] + tab.p50_in[i]) / (1 - tab.cf_in[i]) for i in eachindex(tab.t) if ok[i]]
    rise = pulse_rise((tab.t[ok], total))
    beams_at(t; d = dλ) = [remake(b; dlambda = isin(b) ? d / 3 : 0.0) for b in n221204_fig9_beams(base, tab, t)]
    post(t) = (i = argmin(abs.(tab.t .- t)); abs(tab.t[i] - t) < 0.051 ? tab.cf_post[i] : NaN)
    ts = [t for t in 0.8:0.2:8.4 if isfinite(post(t)) && ok[argmin(abs.(tab.t .- t))]]
    return (; name, tab, dλ, rise, beams_at, post, ts)
end
# N210808 (`n210808=pre109`, default): the cone fractions of PRE 109 Fig. 1e, the same
# 2024 HYDRA model as Fig. 9 (and Fig. 10), with the outer cones' incident powers of
# PRE 106 Fig. 1d (Fig. 1e's incident cone fraction agrees with PRE 106's to 0.004);
# `n210808=pre106`: PRE 106 Fig. 1(d, e) alone (the 2022 model; post-CBET cone
# fraction +0.03 at peak, +0.09 at 3.4 ns).
function n210808_table(source)
    k22 = n210808_fig1()
    source == "pre106" && return k22
    f1 = pre109_fig1()
    at(c, t) = (i = argmin(abs.(k22.t .- t)); abs(k22.t[i] - t) < 0.051 ? c[i] : NaN)
    return (t = f1.t, cf_in = f1.cf_in_n210808, cf_post = f1.cf_post_n210808,
            p23_in = [at(k22.p23_in, t) for t in f1.t], p23_post = [at(k22.p23_post, t) for t in f1.t],
            p30_in = [at(k22.p30_in, t) for t in f1.t], p30_post = [at(k22.p30_post, t) for t in f1.t],
            p44_in = [at(k22.p44_in, t) for t in f1.t], p44_post = [at(k22.p44_post, t) for t in f1.t],
            p50_in = [at(k22.p50_in, t) for t in f1.t], p50_post = [at(k22.p50_post, t) for t in f1.t])
end
shots = [shot("N221204", n221204_fig9(), 2.75), shot("N220919", n220919_fig9(), 2.5),
         shot("N210808", n210808_table(get(opt, "n210808", "pre109")), 1.8)]
ref = shots[1]
for s in shots
    @printf("%s: Δλ %.2f Å, rise %.2f ns (delay %+.2f), %d rows %.1f–%.1f ns\n", s.name, s.dλ,
            s.rise, s.rise - ref.rise, length(s.ts), first(s.ts), last(s.ts))
end
cal = [(; beams_at = s.beams_at, data = [(t, s.post(t)) for t in s.ts], delay = s.rise - ref.rise)
       for s in shots]

# The Δλ null constraint (NC24: 1–1.55 Å for HYBRID-E): with N221204's peak beams the
# inner cones must not gain at 1.0 Å and must gain at 1.55 Å. Squared hinge on the
# inner-power ratio, weight 10 (a 1 % violation costs as much as ~3 cone-fraction
# residuals of 0.01).
peak_beams = [(t, ref.beams_at(t)) for t in (5.6, 7.0)]
function null_penalty(m; lo = 1.0, hi = 1.55, w = 10.0)
    acc = 0.0
    for (t, bs) in peak_beams
        acc += w * max(0.0, receiver_ratio(m, bs, lo / 3; receivers = INNER, time = t) - 1)^2
        acc += w * max(0.0, 1 - receiver_ratio(m, bs, hi / 3; receivers = INNER, time = t))^2
    end
    return acc
end

Te0 = ([1.0, 3.0, 4.0, 4.75, 5.5, 6.5, 8.0], [2790.0, 3550.0, 3720.0, 5190.0, 5500.0, 6030.0, 7440.0])
deck_clamp = ([3.9, 4.1], [0.002, 0.0075])
function start(variant)
    Te = variant == "free" ? Te0 : (Te0[1], min.(Te0[2], 0.95temax))
    if variant == "clamp"
        # start where the hard clamp binds at peak (δn/n of the strongest pairs ~1e-3): a
        # clamp above every pair's δn has zero gradient and the fit leaves it there
        clamp = cbet_clamp_schedule((ref.tab.t, (ref.tab.p44_in .+ ref.tab.p50_in) ./ (1 .- ref.tab.cf_in));
                                    foot = 1e-3, peak = 3e-4)
        Tc = ([1.0, 3.0, 4.0, 4.75, 5.5, 6.5, 8.0], min.([1975.0, 2370.0, 2739.0, 3946.0, 4158.0, 4215.0, 5500.0], 0.95temax))
        return TwoZoneCBET(leh_z = 5.62, outside = (mach = 1.5, multiplier = 0.1, Te = Tc, clamp, saturation = :min),
                           inside = (multiplier = 1.5, Te = Tc, clamp, saturation = :min))
    end
    return TwoZoneCBET(leh_z = 5.62, outside = (mach = 1.23, multiplier = 0.30, Te, clamp = deck_clamp),
                       inside = (multiplier = 1.23, Te, clamp = deck_clamp))
end
fitspec(variant) = variant == "clamp" ?
    ((:outside, :multiplier), (:inside, :multiplier), (:outside, :mach), (:both, :Te), (:both, :clamp)) :
    ((:outside, :multiplier), (:inside, :multiplier), (:outside, :mach), (:both, :Te))
function bounds_for(variant)
    b = copy(Briareus._CBET_BOUNDS)
    variant == "free" || (b[:Te] = (300.0, temax))
    b[:clamp] = (1e-4, 0.05)
    return b
end

residuals(m, s) = [cf(apply_lpi(m, s.beams_at(t), t - (s.rise - ref.rise))) - s.post(t) for t in s.ts]
rms(r) = sqrt(sum(abs2, r) / length(r))
function report(io, name, m)
    println(io, "== ", name)
    z = m.outside
    @printf(io, "outside multiplier %.4f, Mach %.4f; inside multiplier %.4f\n", z.multiplier, z.mach, m.inside.multiplier)
    println(io, "Te (eV) at ", m.inside.Te[1], " ns: ", round.(m.inside.Te[2]))
    println(io, "clamp: ", m.inside.clamp, " saturation ", m.inside.saturation)
    allr = Float64[]
    for s in shots
        r = residuals(m, s)
        append!(allr, r)
        k = argmax(abs.(r))
        @printf(io, "  %s rms %.4f, worst %+.3f at %.1f ns; peak-mean residual (5.4–7.4 ns) %+.4f\n", s.name, rms(r),
                r[k], s.ts[k], sum(r[i] for i in eachindex(r) if 5.4 <= s.ts[i] <= 7.4) /
                               count(t -> 5.4 <= t <= 7.4, s.ts))
    end
    @printf(io, "  all rms %.4f (%d rows); null penalty %.2e\n", rms(allr), length(allr), null_penalty(m))
    for t in (5.6, 7.0, 8.0)
        bs = ref.beams_at(t)
        nl = 3cbet_null(m, bs, t; receivers = INNER)
        c(d) = cf(apply_lpi(m, ref.beams_at(t; d), t))
        @printf(io, "  N221204 beams %.1f ns: null %.2f Å; CF at Δλ 0 / 1 / 1.8 / 2.5 / 2.75 Å: %.3f / %.3f / %.3f / %.3f / %.3f; dCF/dΔλ 1–1.8 %.3f, 2.5–2.75 %.3f per Å\n",
                t, nl, c(0.0), c(1.0), c(1.8), c(2.5), c(2.75), (c(1.8) - c(1.0)) / 0.8, (c(2.75) - c(2.5)) / 0.25)
    end
    # cross-shot check against HYDRA: each shot's own beams at peak
    for t in (5.6, 7.0)
        @printf(io, "  %.1f ns, each shot's beams: ", t)
        for s in shots
            e = apply_lpi(m, s.beams_at(t), t - (s.rise - ref.rise))
            @printf(io, "%s %.3f (HYDRA %.3f)  ", s.name, cf(e), s.post(t))
        end
        println(io)
    end
    # per-cone ratios (not fitted) against each shot's Fig. 9a / K22 Fig. 1d
    for s in shots, t in (5.6, 8.0)
        i = argmin(abs.(s.tab.t .- t))
        bs = s.beams_at(t)
        e = apply_lpi(m, bs, t - (s.rise - ref.rise))
        ratio(g) = sum(b.power for b in e if b.group === g) / sum(b.power for b in bs if b.group === g)
        @printf(io, "  %s %.1f ns post/input 23/30/44/50: %.2f %.2f %.2f %.2f | HYDRA %.2f %.2f %.2f %.2f\n", s.name, t,
                ratio(:inner23), ratio(:inner30), ratio(:outer44), ratio(:outer50),
                s.tab.p23_post[i] / s.tab.p23_in[i], s.tab.p30_post[i] / s.tab.p30_in[i],
                s.tab.p44_post[i] / s.tab.p44_in[i], s.tab.p50_post[i] / s.tab.p50_in[i])
    end
    flush(io)
end

# Gauss–Newton (Laplace) covariance in the fit's parameter space: residual Jacobian J
# by forward-mode AD, σ² from the residuals (model error included), parameters at a
# bound held fixed. Returns the standard deviations of the fitted values (multiplicative
# for log-space parameters: exp(±σ_u)).
function laplace(io, m, variant)
    u0, ulo, uhi, with = Briareus._calibration_space(m, fitspec(variant), bounds_for(variant))
    free = findall(i -> u0[i] > ulo[i] + 1e-6 && u0[i] < uhi[i] - 1e-6, eachindex(u0))
    rvec(v) = (u = collect(promote_type(eltype(v), Float64), u0); u[free] = v;
               mm = with(u); vcat([[cf(apply_lpi(mm, s.beams_at(t), t - (s.rise - ref.rise))) - s.post(t)
                                     for t in s.ts] for s in shots]...))
    t0 = time()
    J = ForwardDiff.jacobian(rvec, u0[free])
    r = rvec(u0[free])
    σ2 = sum(abs2, r) / (length(r) - length(free))
    C = σ2 * inv(Symmetric(J' * J))
    @printf(io, "Laplace (%.0f s): σ_resid %.4f, %d free of %d parameters (log-space σ, i.e. ×exp(±σ)):\n",
            time() - t0, sqrt(σ2), length(free), length(u0))
    println(io, "  free index: ", free)
    println(io, "  u: ", round.(u0[free], digits = 4))
    println(io, "  σ_u: ", round.(sqrt.(diag(C)), digits = 4))
    corr = C ./ sqrt.(diag(C) * diag(C)')
    println(io, "  correlations:")
    for i in axes(corr, 1)
        println(io, "   ", join([@sprintf("%+.2f", c) for c in corr[i, :]], " "))
    end
    # derived: the null and the peak slope, by finite differences along the covariance
    q(mm) = (bs = ref.beams_at(5.6); c(d) = cf(apply_lpi(mm, ref.beams_at(5.6; d), 5.6));
             [3cbet_null(mm, bs, 5.6; receivers = INNER), (c(1.8) - c(1.0)) / 0.8, (c(2.75) - c(2.5)) / 0.25])
    q0 = q(m)
    G = zeros(3, length(free))
    for (j, i) in enumerate(free)
        h = 1e-3
        up = copy(u0); up[i] += h
        dn = copy(u0); dn[i] -= h
        G[:, j] = (q(with(up)) - q(with(dn))) / 2h
    end
    Cq = G * C * G'
    @printf(io, "  null at 5.6 ns %.3f ± %.3f Å; dCF/dΔλ 1–1.8 Å %.3f ± %.3f; 2.5–2.75 Å %.3f ± %.3f per Å\n",
            q0[1], sqrt(Cq[1, 1]), q0[2], sqrt(Cq[2, 2]), q0[3], sqrt(Cq[3, 3]))
    flush(io)
end

for variant in variants
    tag = variant * (get(opt, "n210808", "pre109") == "pre106" ? "_pre106" : "")
    open("hybride_cbet_$(tag).txt", "w") do io
        tee = IOContext(io)
        m0 = start(variant)
        t0 = time()
        pen = variant in ("free", "tebound") ? nothing : null_penalty
        m, res = calibrate_cbet(m0, cal; receivers = INNER, fit = fitspec(variant), bounds = bounds_for(variant),
                                penalty = pen, iterations)
        @printf(io, "%s: %d iterations, %d evaluations, %.0f s, loss %.5f, converged %s\n", variant,
                Briareus.Optim.iterations(res), Briareus.Optim.f_calls(res), time() - t0,
                Briareus.Optim.minimum(res), Briareus.Optim.converged(res))
        report(io, variant, m)
        report(stdout, variant, m)
        println(io, "model = ", repr(m))
        get(opt, "laplace", "on") == "on" && laplace(io, m, variant)
        if get(opt, "profile", "off") == "on" && !(variant in ("free", "tebound"))
            # profile: best loss with the null pinned near each Δλ (± 0.05 Å), warm start
            for nl in (0.5, 0.8, 1.2, 1.55, 2.0)
                pp(mm) = null_penalty(mm; lo = nl - 0.05, hi = nl + 0.05, w = 100.0)
                mp, rp = calibrate_cbet(m, cal; receivers = INNER, fit = fitspec(variant),
                                        bounds = bounds_for(variant), penalty = pp, iterations = 60)
                r = vcat([residuals(mp, s) for s in shots]...)
                @printf(io, "profile null %.2f Å: loss %.5f, rms %.4f (%s)\n", nl, Briareus.Optim.minimum(rp), rms(r),
                        join([@sprintf("%s %.4f", s.name, rms(residuals(mp, s))) for s in shots], ", "))
                flush(io)
            end
        end
    end
    print(read("hybride_cbet_$(tag).txt", String))
end
