# N221204-class time-dependent benchmark: capsule flux P2/P0(t) of a HYBRID-E hohlraum
# with CBET, a moving wall, closing entrance holes and inner-beam bubble absorption, each
# piece added in turn, against the HYDRA curve of Kritcher et al., PRE 109, 025204 (2024)
# Fig. 10. Inputs: reference/litreview/n221204_benchmark.md §5; results: research notes,
# "N221204 time-dependent benchmark (2026-10-23)" and "Wall material, conversion and
# capsule energetics (2026-10-24)".
#
#   julia --project=. -t auto examples/n221204.jl [h capsule_h [runs [transfer [table [shot]]]]]
#
# `h`, `capsule_h`: mesh (default 0.28 / 0.16, the standard mesh); `runs`: `all`
# (default) or `full` (only the full model, e.g. for a fine-mesh confirmation);
# `transfer`: `calibrated` (the default, `joint` for Zylstra's pair below; the two-zone
# reduced CBET calibrated to PRE 109 Fig. 9, `n221204_cbet()`), `fig9b` (the
# outer-to-inner transfer that gives PRE 109 Fig. 9b's
# post-transfer cone fraction, HYDRA's own CBET result for this shot), `joint` (the
# joint three-shot calibration, `hybride_cbet()`) or `iraum` (the same with the Δλ null
# at I-Raum's ≤ 0.5 Å, `hybride_cbet_iraum()`); `table` as 5th
# argument solves only the comparison times (any other word, e.g. `x`, the default
# times); `shot`: `n221204` (default),
# `n220919` or `n210808` (proxies, see below), or `n191118` / `n210601`, Zylstra 2022's
# LEH-size pair of symcaps (3.64 and 3.1 mm LEHs, their measured pulses and Dante;
# research notes, "Zylstra 2022 LEH pair, time-dependent (2026-10-25)"). Writes
# `<shot>_<transfer>_<h>_<capsule_h>[_<options>].csv` (P2/P0, P4/P0 and
# capsule Tr per run) to the current directory.
#
# Energetics (research notes, "Wall material, conversion and capsule energetics
# (2026-10-24)", "Energetics defaults and the capsule share (2026-10-25)" and "Capsule
# absorption, conversion and the ablating-wall release (2026-10-25)"): the Au-lined DU
# wall of a NIF drive history (`nif_uranium()`: `du_wall_loss()` and the rising
# `nif_conversion()`), each wall element's loss from its own temperature history with
# re-emission under a falling drive (`drive_history`'s defaults), the capsule absorbing
# at the effective radius of its ablation front with HYDRA's ablator albedo
# (`CapsuleAlbedo`, below), HYDRA's laser degradation and, on N221204 and its proxies,
# the measured backscatter. Options, as `key=value` arguments anywhere on the line:
# - `capsule=model|mesh`: `CapsuleAlbedo` with albedo `capsule_albedo` absorbing at the
#   effective radius (default), or the capsule mesh itself shrinking along the schedule
#   on the keyframes (the flux on the shrinking sphere; the static runs keep R₀);
#   `capsule_albedo=young|<number>`: α_c, Young 2020's quasi-1D HYDRA ablator albedo on
#   a time map to the shot's bang (default; published, not fitted) or a constant;
#   `capsule_outside=<a>`: `CapsuleAlbedo`'s `outside`, the absorptance of the blow-off
#   outside the front (default 0: unpublished);
# - `degrade=<fraction>`: HYDRA's laser degradation (PRE 109 App. 3: 9–13 %), default
#   0.11, the middle of the range.
# Sensitivities of the P2 swing (research notes, "The P2 swing: capsule radius,
# backscatter and the collision (2026-10-25)"):
# - `backscatter=measured|history|uniform|off`: the measured inner-cone SBS/SRS energies
#   on shapes of our own (default), the same energies on N210808's measured per-cone SBS
#   time histories (`history`, Kemp 2026 Fig. 2), or a uniform 48 kJ power cut on every
#   beam, the treatment PRE 109's HYDRA used (with `measured` or `history`, the full model
#   is also run with the uniform cut, for Fig. 10); `sbs_energy=<kJ>` (default 42);
# - `refraction=<f>`, `refraction_shift=<ns>`: a sensitivity, the inner power lost at the
#   late-time on-axis density spike (Kemp 2026 §III C; research notes, "Inner-cone SBS time
#   histories and the late-time refraction (2026-10-25)");
# - `bubble_energy=<E>` (default 1.1), `collision=<ns>` (default 6.09): the bubble's
#   picket fluence and the ablator–wall collision time.
# Further options (research notes, "N221204 late-time drive (2026-10-24)"):
# - `lips=prescribed|computed`: the lips close along the prescribed ramp (default) or by
#   the computed `LipClosure`; `lag=<ns>`: the power through a first-order lag, a
#   sensitivity estimate of plasma-stored energy (research notes, "Computed LEH closure
#   (2026-10-24)");
# - `unfold=on`: Dante 1's spectrum, channel voltages and their unfold at 5–8.5 ns in the
#   full model (research notes, "Dante channel responses and the unfold (2026-10-24)");
# - `metrics=on`: the steps run on to bang time and the full model prints ROADMAP 11's
#   implosion metrics (hot-spot velocity, flux P2 over the drive window, the drive at
#   bang) over the window from the rise to bang (research notes, "Zylstra 2022 LEH pair,
#   time-dependent (2026-10-25)");
# - `mband=on`: the full model's capsule M-band (> 1.8 keV) P2/P0 at every step from each
#   cone's own post-CBET, post-bubble deposition (`drive_history`'s `mband`), beside the
#   earlier proxy (the cones weighted by their pre-CBET deposition), against Dewald 2020's
#   +0.15 at 6 ns on the Mo sphere (a check; research notes, "Per-cone M-band in the drive
#   history (2026-10-27)").
# Every run prints its time-integrated energy budget (`energy_budget`) against PRE 109
# App. 3's N221204 split (wall 54.1 %, LEH 25.2 %, gas 4.5 %, capsule 16.2 %), and its
# synthetic Dante 1 Tr peak (`dante1()`) beside the capsule-seen Tr, against the measured
# 313 eV (research notes, "Spectral radiosity, M-band and a synthetic Dante").

using Briareus, Printf, DelimitedFiles

pos = filter(a -> !occursin('=', a), ARGS)
opt = Dict(String.(split(a, '=', limit = 2)) for a in ARGS if occursin('=', a))
given = sort(collect(keys(opt)))          # the options on the command line (the file tag)
shot = get(pos, 6, "n221204")
# Zylstra 2022's LEH pair (`n191118`, `n210601`): no published backscatter, so none for
# them (N221204's measured SBS/SRS would not apply).
pair = shot in ("n191118", "n210601")
for k in keys(opt)
    k in ("capsule", "capsule_albedo", "capsule_outside", "degrade", "lips", "lag", "unfold", "metrics", "dt", "alpha_m", "mband",
          "backscatter", "collision", "bubble_energy", "sbs_energy", "refraction", "refraction_shift") ||
        error("unknown option $k=$(opt[k])")
end
h = parse(Float64, get(pos, 1, "0.28"))
capsule_h = parse(Float64, get(pos, 2, "0.16"))
which = get(pos, 3, "all")
transfer = get(pos, 4, pair ? "joint" : "calibrated")   # the pair: the joint multi-shot fit (below)
transfer in ("calibrated", "fig9b", "joint", "iraum") ||
    error("transfer calibrated, fig9b, joint or iraum, not $transfer")

# Geometry (mm; n221204_benchmark.md §1, PRE 109 Fig. 1a): 11.24 × 6.4 mm Au-lined DU
# hohlraum, 3.1 mm LEHs, 1.135 mm capsule. Unpublished, so repository defaults: square
# corners, no diagnostic windows. The wall is `nif_uranium()`: the Au-lined DU loss with
# the rising conversion of its gold liner.
# The capsule: a sphere absorbing at its ablation front's radius R(t) inside the fixed
# R₀ mesh (`CapsuleAlbedo`'s `radius`). HYDRA's front for N221204 is unpublished; at 6 ns
# it lies between 0.84 R₀ (no ablator left outside the ablator–fuel interface,
# 0.906 × 1.050/1.135) and 0.906 R₀ (the shell keeping its relative thickness). The
# schedule (research notes, "Capsule absorption, conversion and the ablating-wall release
# (2026-10-25)"): N221204's own HYDRA radius history (Kritcher et al., PoP 31, 070502
# (2024) Fig. 1, `kritcher2024_radius()`, the ablator–fuel interface or fuel centre of
# mass scaled to its start) applied to the whole capsule, on the shot's own time axis: the
# upper bound on the front, with the published timing (held until the first shock breaks
# out, ~3.6 ns; minimum at bang). Zylstra's pair: K24's thinner N210808 capsule (their
# design) on the same drive, its time axis scaled so that its minimum (9.18 ns) falls at
# the shot's bang time (`ScaledRadius`, derived, *L*). (`schedule=young`, Young 2020's
# N170827 front, removed 2026-10-26: it lies below the lower bound.)
zy = zylstra2022_fig2()
bang = pair ? getfield(zy.bang, Symbol(shot)) : 9.53
capsule_schedule = if pair
    ScaledRadius(kritcher2024_radius(:n210808); time_scale = bang / 9.18)
else
    kritcher2024_radius(:n221204)
end
# The albedo α_c (`capsule_albedo`): Young 2020's quasi-1D HYDRA ablator albedo of the
# LLNL dynamic view-factor model for N170827 (Fig. 3a top, `young2020_albedos()`: 0.01–0.16
# to 6 ns, 0.30 at 7, 0.40 from 7.5 ns), its time axis mapped piecewise linearly so that
# the model's 5.75 ns falls at 6.0 ns and its 7.75 ns at bang time (default; published,
# not fitted), or a constant. `capsule_outside`: the absorptance of the blow-off outside the
# front (`CapsuleAlbedo`'s `outside`), default 0. `capsule=mesh`: the schedule shrinks the
# capsule mesh itself on the keyframes (`moving_hohlraum`'s `capsule_shrink`; the flux is
# evaluated on the shrinking sphere) and the albedo is α_c alone; `capsule=model` keeps
# the mesh at R₀ and absorbs at the effective radius (`CapsuleAlbedo`'s `radius`).
ymap(τ) = τ <= 5.75 ? τ * 6.0 / 5.75 : 6.0 + (τ - 5.75) * (bang - 6.0) / (7.75 - 5.75)
ya = young2020_albedos()
capsule_alpha = get(opt, "capsule_albedo", "young") == "young" ? TabulatedAlbedo(ymap.(ya.t), ya.capsule) :
                parse(Float64, opt["capsule_albedo"])
capsule_kind = get(opt, "capsule", "model")
capsule_kind in ("model", "mesh") || error("capsule=model or capsule=mesh, not $capsule_kind")
capsule = Material(:ablator; conversion = 0.3,
                   albedo = CapsuleAlbedo(albedo = capsule_alpha,
                                          radius = capsule_kind == "model" ? capsule_schedule : nothing,
                                          outside = parse(Float64, get(opt, "capsule_outside", "0"))))
shrink = capsule_kind == "mesh" ? (; capsule_shrink = capsule_schedule) : (;)
# Zylstra's pair: N191118's 3.64 mm LEH (N210601's 3.1 mm) and the N210808-design capsule,
# 1.1285 mm (K22 Fig. 4; N191118's, of the N210207 design, taken the same, *L*).
leh_radius = shot == "n191118" ? 1.82 : 1.55
capsule_radius = pair ? 1.1285 : 1.135
geometry = (length = 11.24, leh_radius, capsule_radius, h, capsule_h,
            wall = nif_uranium(), capsule)
R0 = 3.2

# Laser: 192 beams at the public mock HYBRID-E pointing; Δλ 2.75 Å (1ω) on the inner
# cones, 0.917 Å at 3ω (the 23.5°/30° and 44.5°/50° splits are unpublished: none).
# `shot` (7th argument) `n220919` or `n210808` runs a proxy of the two other shots of
# Fig. 10 (*L*): Δλ 2.5 or 1.8 Å on the same mock pointing, CBET calibration and
# hohlraum; N220919 on N221204's pulse (the same pulse with re-tuned foot powers,
# K24), N210808 with the rise 0.2 ns earlier and the peak ending 0.45 ns earlier (P24
# Fig. 1, read by eye). Zylstra's pair: Δλ 1.25 Å (N191118) and 1.8 Å (N210601) (Zylstra
# 2022), on the same mock pointing (N191118's outer beams sat farther from the midplane
# for its larger LEH; unpublished, *L*).
dλ = Dict("n221204" => 2.75, "n220919" => 2.5, "n210808" => 1.8, "n191118" => 1.25, "n210601" => 1.8)[shot]
beams = [remake(b; dlambda = b.group in (:inner23, :inner30) ? dλ / 3 : 0.0)
         for b in hybride_mock_pointing(nif_beams())]

# Delivered outer (128 beams) and inner (64 beams) power, DN24 Fig. 1 read by eye (TW):
# N210808's pulse (shifted, above) and the cone fraction of every shot's split.
pt = [0.0, 0.4, 0.8, 1.3, 2.0, 2.5, 3.3, 3.9, 4.5, 5.0, 5.5, 8.0, 8.3, 8.6, 8.9, 9.1, 9.25, 9.35]
shot == "n210808" && (pt = [t >= 8.0 ? t - 0.45 : t >= 5.0 ? t - 0.2 : t for t in pt])
outer = [0, 2, 50, 50, 25, 26, 60, 65, 150, 270, 296, 294, 250, 110, 40, 0, 0, 0]
inner = [0, 2, 4, 4, 4.5, 5, 7, 9, 25, 90, 146, 146, 146, 146, 146, 146, 140, 0]
total = Float64.(outer .+ inner)
# N221204 and N220919: the total power as delivered, read from the vector paths of Pak et
# al., PRE 109, 025203 (2024) Fig. 1 (`n221204_pulse()`, 2.051 MJ), split by DN24's cone
# fraction above (the split is published only there).
measured = n221204_pulse()
if shot in ("n221204", "n220919")
    oi, ii = PulseShape(pt, Float64.(outer)), PulseShape(pt, Float64.(inner))
    frac(τ) = (s = oi(τ) + ii(τ); s > 0 ? ii(τ) / s : τ > 5 ? 1.0 : 0.5)
    pt, total = measured.t, measured.power
    inner = total .* frac.(pt)
elseif pair
    # Zylstra 2022 Fig. 2's total power (`zylstra2022_fig2()`), split by DN24's cone
    # fraction with its rise and droop moved to the shot's own (the split of these shots
    # is unpublished, *L*): shot time maps piecewise linearly onto N221204's through the
    # starts of the rise (`pulse_rise`) and the half-power points of the fall.
    oi, ii = PulseShape(pt, Float64.(outer)), PulseShape(pt, Float64.(inner))
    cfrac(τ) = (s = oi(τ) + ii(τ); s > 0 ? ii(τ) / s : τ > 5 ? 1.0 : 0.5)
    fall(t, P) = t[findlast(>=(maximum(P) / 2), P)]
    ps = Float64.(getfield(zy, Symbol("power_", shot)))
    knots_s = [0.0, pulse_rise((zy.t, ps)), fall(zy.t, ps)]
    knots_r = [0.0, pulse_rise((measured.t, measured.power)), fall(measured.t, measured.power)]
    to_ref(τ) = τ >= knots_s[3] ? τ - knots_s[3] + knots_r[3] :
                (k = τ >= knots_s[2] ? 2 : 1;
                 knots_r[k] + (τ - knots_s[k]) * (knots_r[k+1] - knots_r[k]) / (knots_s[k+1] - knots_s[k]))
    pt, total = collect(zy.t), ps
    inner = total .* cfrac.(to_ref.(pt))
    measured = (; t = zy.t, power = ps, dante1 = getfield(zy, Symbol("dante_", shot)))
    @printf("%s: rise %.2f ns, half-power fall %.2f ns (N221204 %.2f, %.2f); bang %.2f ns\n",
            shot, knots_s[2], knots_s[3], knots_r[2], knots_r[3], bang)
end
# `lag=<ns>`: a sensitivity estimate of energy stored in the hohlraum plasma (gas and
# gold blow-off) and given back as x-rays: the delivered power (total and inner) is
# replaced by its first-order lag of that time constant, so the stored energy is
# P·lag (energy conserved, though the tail past the last step is not solved; CBET,
# backscatter and the deposition see the lagged power too). *L*, a sensitivity, not a
# calibrated model; research notes, "Computed LEH closure (2026-10-24)".
lag = parse(Float64, get(opt, "lag", "0"))
if lag > 0
    tl = collect(0.0:0.005:pt[end] + 5lag)
    Pt, Pi = PulseShape(pt, Float64.(total)), PulseShape(pt, Float64.(inner))
    function lagged(f)
        y = zeros(length(tl))
        for k in 2:length(tl)
            δ = tl[k] - tl[k-1]
            y[k] = y[k-1] + (f(tl[k]) - y[k-1]) * (1 - exp(-δ / lag))
        end
        return y
    end
    pt, total, inner = tl, lagged(Pt), lagged(Pi)
end
pulses = cone_fraction_pulses(beams, pt, total, inner ./ max.(total, eps()))
beams_at = τ -> with_pulses(beams, pulses, τ)
@printf("laser energy %.3f MJ, peak %.0f TW\n",
        sum(pulse_energy(pulses[b.group]) for b in beams) / 1e3, maximum(total))

# Backscatter (N221204 and its proxies) and degradation, applied before CBET and the bubble.
# SBS: 42 ± 11 kJ "early in peak power and mostly on the inner cones", SRS 6 ± 1 kJ
# "throughout peak" (AS24): on the inner cones, SBS as a fraction rising over 5.0–5.5 ns,
# flat to 6.5 ns and gone by 7.5 ns, SRS flat over 5.5–9.0 ns (the shapes are *L*), each
# scaled to its measured energy on the delivered inner power. Degradation: every beam's
# power × (1 − degrade), the HYDRA multiplier of PRE 109 App. 3 (9–13 %).
inner_power = PulseShape(pt, Float64.(inner))
fine = range(0.0, 9.4, 9401)
∫(f) = sum(f(τ) for τ in fine) * step(fine)
sbs_shape = PulseShape([0.0, 5.0, 5.5, 6.5, 7.5, 20.0], [0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
srs_shape = PulseShape([0.0, 5.4, 5.5, 9.0, 9.1, 20.0], [0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
# `sbs_energy=<kJ>` (default 42, AS24's measured SBS): the SBS energy of `measured`,
# `history` and (with the 6 kJ SRS) `uniform`, below.
sbs_energy = parse(Float64, get(opt, "sbs_energy", "42"))
f_sbs = sbs_energy / ∫(τ -> sbs_shape(τ) * inner_power(τ))
f_srs = 6.0 / ∫(τ -> srs_shape(τ) * inner_power(τ))
bs_fraction = PulseShape(collect(fine), [f_sbs * sbs_shape(τ) + f_srs * srs_shape(τ) for τ in fine])
# `backscatter=history`: the measured per-cone SBS time histories of N210808 (Kemp et al.,
# PoP 33, 082705 (2026) Fig. 2, `kemp2026_sbs()`; research notes, "Inner-cone SBS time
# histories and the late-time refraction (2026-10-25)"): all backscattered SBS power of a
# 23° and of a 30° quad, 2.8 and 0.9 kJ per quad, starting as the inner power reaches
# peak, peaking at 5.3 and 5.75 ns and gone by 6.6–6.9 ns. On the shot's time axis: shifted
# by the time the shot's inner power reaches half its peak, less N210808's 4.74 ns on the
# paper's axis (*L*: the shut-off follows the hydrodynamics, timed here by the rise);
# scaled to `sbs_energy` with N210808's 23°:30° split (3.1:1, *L* for N221204, whose
# per-cone split is unpublished). Each cone's fraction is its SBS power over its share of
# the delivered inner power (eight quads per cone, equal quad powers). SRS as above.
k26 = kemp2026_sbs()
function half_rise(t, P)
    k = findfirst(>=(maximum(P) / 2), P)
    return t[k-1] + (maximum(P) / 2 - P[k-1]) * (t[k] - t[k-1]) / (P[k] - P[k-1])
end
sbs_shift = half_rise(pt, Float64.(inner)) - 4.74
E_quad = (sum(k26.inner23) + sum(k26.inner30)) * step(range(k26.t[1], k26.t[end], length(k26.t)))   # kJ, a 23° and a 30° quad
sbs_scale = sbs_energy / (8E_quad)
sbs_cone = Dict(g => PulseShape(k26.t .+ sbs_shift, getfield(k26, g)) for g in (:inner23, :inner30))
history_fraction(g) = PulseShape(collect(fine), [inner_power(τ) > 0 ?
    min(1.0, sbs_scale * sbs_cone[g](τ) / (inner_power(τ) / 16)) + f_srs * srs_shape(τ) : 0.0 for τ in fine])
degrade = parse(Float64, get(opt, "degrade", "0.11"))
# `backscatter=measured|history|uniform|off`: the measured SBS and SRS on the inner cones
# with the shapes above (default, the physics: modeled backscatter is per cone and time
# dependent; none for Zylstra's pair), the same energies on the published per-cone SBS
# histories (`history`, above), or a uniform 48 kJ power cut, the same fraction of every
# beam's power: one simple treatment, the one PRE 109's HYDRA used (backscatter "as part
# of an effective total laser power multiplier", §II), so the like-for-like setting for
# Fig. 10. With `measured` or `history` the full model is run a second time with the
# uniform cut, and Fig. 10 is compared with that run.
bs_mode = get(opt, "backscatter", pair ? "off" : "measured")
bs_mode in ("measured", "history", "uniform", "off") ||
    error("backscatter=measured, history, uniform or off, not $bs_mode")
backscatter = bs_mode in ("measured", "history")
measured_cut = Backscatter(:inner23 => bs_fraction, :inner30 => bs_fraction)
history_cut = Backscatter(:inner23 => history_fraction(:inner23), :inner30 => history_fraction(:inner30))
uniform_cut = Backscatter((g => (sbs_energy + 6.0) / sum(pulse_energy(pulses[b.group]) for b in beams)
                           for g in (:inner23, :inner30, :outer44, :outer50))...)
degradation = degrade > 0 ? Any[Backscatter((g => degrade for g in (:inner23, :inner30, :outer44, :outer50))...)] : Any[]
# `refraction=<f>` (default 0), a sensitivity: the inner power lost at the on-axis density
# spike between capsule and LEH that refracts the inner beams into the gold bubble and
# absorbs them (Kemp 2026 §III C): a fraction f of the 30° power at the spike's peak,
# following the spike's density on axis in the paper's MSH simulation of N210808 (Fig. 6a
# inset: n_e/n_c 0.01 to 6.0 ns, 0.02 at 6.5, 0.08 at 6.75, 0.15 at 7.0, 0.14 at 7.75; held
# after that, *L*),
# the 23° cone the same 0.25 ns later (its SBS ends 0.25 ns after the 30°'s), on the
# shot's axis as `history`; `refraction_shift=<ns>` moves it (−1: the paper's hypothesis
# that the spike forms ~1 ns before the hydrodynamics predicts). The paper's only number:
# without the 7 ns spike, twice the 30° power reaches the ablator (AAC 2025, OSTI
# 3388898), so f = 0.5 there. The light is booked as lost (an upper bound on its effect
# on the wall drive: part of it is absorbed in the bubble and converts there). *L*.
refraction = parse(Float64, get(opt, "refraction", "0"))
refraction_shift = parse(Float64, get(opt, "refraction_shift", "0"))
spike = PulseShape([0.0, 6.0, 6.5, 6.75, 7.0, 7.75, 20.0], [0.0, 0.0, 0.01, 0.07, 0.14, 0.13, 0.13] ./ 0.14)
spike_at(τ, lag) = spike(τ - sbs_shift - refraction_shift - lag)
refraction_loss = Backscatter(:inner30 => PulseShape(collect(fine), [refraction * spike_at(τ, 0.0) for τ in fine]),
                              :inner23 => PulseShape(collect(fine), [refraction * spike_at(τ, 0.25) for τ in fine]))
lpi_pre(mode) = Any[(mode == "measured" ? [measured_cut] : mode == "history" ? [history_cut] :
                     mode == "uniform" ? [uniform_cut] : [])...;
                    (refraction > 0 ? [refraction_loss] : [])...; degradation...]
pre = lpi_pre(bs_mode)
# the laser–plasma models of a run: `p` (`pre`) first; exactly the given model(s) without it
with_pre(ms...; p = pre) = isempty(p) ? (length(ms) == 1 ? ms[1] : ms) : (p..., ms...)
bs_mode == "measured" &&
    @printf("backscatter: SBS %.1f %% of the inner power over 5.5–6.5 ns, SRS %.2f %% over the peak\n", 100f_sbs, 100f_srs)
if bs_mode == "history"
    e(g) = ∫(τ -> sbs_scale * sbs_cone[g](τ)) * 8
    pk(g) = maximum(history_fraction(g).(fine)) - f_srs
    @printf("backscatter: per-cone SBS histories (Kemp 2026, N210808) shifted %+.2f ns, ×%.2f: 23° %.1f kJ (peak %.1f %%), 30° %.1f kJ (peak %.1f %%); SRS %.2f %% over the peak\n",
            sbs_shift, sbs_scale, e(:inner23), 100pk(:inner23), e(:inner30), 100pk(:inner30), 100f_srs)
end
refraction > 0 && @printf("refraction: %.0f %% of the 30° power lost at the spike's peak (%.2f ns on this axis)\n",
                          100refraction, 7.0 + sbs_shift + refraction_shift)
degrade > 0 && @printf("degradation %.0f %%\n", 100degrade)

# CBET (`transfer=reduced`, the single-zone model at the deck's settings, which reverses
# the transfer, removed 2026-10-26). `fig9b`: move the fraction of the outer power that turns the input cone
# fraction into PRE 109 Fig. 9b's post-transfer one (read by eye, n221204_benchmark.md
# §2 and the Fig. 9b notes), shared by the inner beams in proportion to their power.
cf9b = PulseShape([0.0, 1.0, 2.0, 2.5, 3.3, 3.9, 4.5, 5.0, 5.6, 7.0, 8.0, 8.5, 9.0, 9.4],
                  [0.17, 0.17, 0.22, 0.20, 0.20, 0.30, 0.40, 0.55, 0.65, 0.55, 0.52, 0.52, 0.90, 1.0])
inner_in, total_in = PulseShape(pt, Float64.(inner)), PulseShape(pt, total)
function fig9b_fraction(Δλ, t)
    c = inner_in(t) / max(total_in(t), eps())
    return clamp((cf9b(t) - c) / max(1 - c, eps()), 0.0, 1.0)
end
fig9b = CalibratedTransfer((:outer44, :outer50), (:inner23, :inner30), fig9b_fraction)
# The two-zone model (supersonic outflow outside the LEH, slow plasma inside) with its
# multipliers, outside Mach number and LEH electron-temperature history calibrated to
# PRE 109 Fig. 9b (research notes, "Two-zone CBET calibrated to PRE 109 Fig. 9"). The
# calibration used HYDRA's degraded input powers (~370 TW at peak); here it acts on the
# delivered 442 TW.
calibrated = n221204_cbet(leh_z = geometry.length / 2)
# `joint`: the two-zone model fitted jointly to N210808, N220919 and N221204 with the Te
# and Δλ-null constraints (`hybride_cbet()`; research notes, "Joint multi-shot CBET
# calibration"); its time tables follow each shot's rise (N210808's 0.2 ns early).
# Zylstra's pair: the tables follow the shot's own rise (`pulse_rise`, against N221204's).
cbet_delay = shot == "n210808" ? -0.2 :
             pair ? pulse_rise((pt, total)) - pulse_rise((n221204_pulse().t, n221204_pulse().power)) : 0.0
joint = hybride_cbet(leh_z = geometry.length / 2, delay = cbet_delay)
# `iraum`: the same joint fit with the null held at I-Raum's ≤ 0.5 Å, NC24's 1–1.55 Å
# being the Δλ of round hot spots rather than the null (`hybride_cbet_iraum()`; research
# notes, "Joint multi-shot CBET calibration: the null")
iraum = hybride_cbet_iraum(leh_z = geometry.length / 2, delay = cbet_delay)
cbet = transfer == "fig9b" ? fig9b : transfer == "joint" ? joint : transfer == "iraum" ? iraum : calibrated

# Hohlraum dynamics, calibrated against published data other
# than Fig. 10 (research notes, "Hohlraum dynamics calibrated to published data
# (2026-10-24)"):
# - wall: `WallMotion`'s defaults (fitted to the LLNL dynamic view-factor radius of Young
#   2020) from 0.8 ns after the outer turn-on (0.4 ns), stopping 0.38 mm in, HYBRID-E's
#   own HYDRA waist edge at 6 ns (PRE 109 Fig. 1a, `pre109_fig1a_wall()`);
# - LEH lips: closing to 0.857 r₀ along the smoothstep fitted to Young 2020's LEH radius
#   (closure −0.90 ns, width 7.85 ns), shifted 0.8 ns later for N221204's later peak
#   (peak reached at ~5.4 ns against ~4.6);
# - bubble: `BubbleAbsorption`'s defaults (Izumi 2021 Figs. 8a and 10a) with an outer
#   picket fluence 1.1× the HDC basecamp's: 1.5× the energy (~40 against ~26 kJ in the
#   first ns after the turn-on, DN24 Fig. 1 and Izumi 2021 Fig. 3a read by eye) spread
#   over split outer quads (−27 % fluence on Izumi's split shots); the tip, 1.97 mm in at
#   6 ns, against 1.85 ± 0.1 for N210808 in PRE 109 Fig. 1a (N221204's is hidden under
#   its rays, ≥ 1.55). The collision time scales with HYBRID-E's capsule–wall gap
#   (3.2 − 1.135 mm): t_c 6.09 ns.
# (`late_lips`, the 2026-10-23 lip closure at 5.0 ns, and with it the `dynamics`
# argument, removed 2026-10-26.)
turn_on = 0.4
motion = WallMotion(start = turn_on + 0.8, saturation = 0.38)
lip, r_open = GapClosure(closure = -0.1, width = 7.85, after = WallLossAlbedo()), 0.857
# `bubble_energy=<E>` (default 1.1) and `collision=<ns>` (default the gap scaling,
# 6.09 ns) for sensitivity runs (research notes, "The P2 swing: capsule radius,
# backscatter and the collision (2026-10-25)")
bkw = haskey(opt, "collision") ? (; collision = parse(Float64, opt["collision"])) : (;)
bubble = BubbleAbsorption(; radius = R0, turn_on, wall = motion, capsule_radius,
                          energy = parse(Float64, get(opt, "bubble_energy", "1.1")), bkw...)
# `lips=computed`: the lips' closure computed instead of prescribed (`LipClosure`, its
# speed, depth and hold fitted to Chen 2020's measured LEH radii, `chen2020_leh()`): free
# expansion at the lip's own gold sound speed in the foot, held 0.24 mm inside the rim
# while the beams pass, free again as the power falls; the lip reaches down to 0.6 r₀
# (research notes, "Computed LEH closure (2026-10-24)").
get(opt, "lips", "prescribed") == "computed" && (lip = :computed; r_open = 0.6)
geometry = (geometry..., leh_open_radius = r_open * geometry.leh_radius)
@printf("dynamics: wall %.2f mm in at 6 ns, %.2f at 9 ns; bubble t_c %.2f ns, tip %.2f mm at 6 ns\n",
        wall_displacement(motion, 6.0), wall_displacement(motion, 9.0), bubble.collision,
        bubble_depth(bubble, 6.0))
# Every step is quasi-static (no burnthrough fronts here), so `table` (5th argument)
# solves only the times of the comparison table, for a fine-mesh check.
# The default steps run to the end of the pulse (9.35 ns), for the energy budget.
# `dt=<ns>` (default 0.25) sets the step of the default times (a step-size check)
tstep = parse(Float64, get(opt, "dt", "0.25"))
times = get(pos, 5, "") == "table" ? [2.0, 4.0, 4.5, 5.5, 6.5, 7.0, 8.0, 8.5] :
        sort(unique([collect(0.5:tstep:9.25); 9.35; 0.5:0.25:9.25]))
# `metrics=on`: the drive window runs on to bang time (N221204 9.53 ns, K24E Table I; the
# pair's from Zylstra Fig. 2), the full model records every Y_lm, and ROADMAP 11's
# implosion metrics are printed over the window from the rise to bang (below).
metrics = get(opt, "metrics", "off") == "on"
metrics && (times = sort(unique([times; filter(>(times[end]), [9.45, bang]); bang])))
keyframes = range(turn_on, 9.0, 11)

p2(r) = r.moments[3, :] ./ r.moments[1, :]
p4(r) = r.moments[5, :] ./ r.moments[1, :]
# Synthetic Dante 1 (θ 143°, φ 274°, 37° through the lower LEH; `dante1()`): the
# brightness temperature of the power per steradian it sees, over the LEH's initial
# area (Dante's convention) and over its open area (the closing lips' open fraction: the
# temperature behind the hole), recorded at every step of every run.
d1 = dante1()
dante_tr(sol) = (d = dante(d1, sol); (d.Tr, d.Tr_open))
# `unfold=on`: in the full model, at the times below, also the spectrum Dante 1 sees
# (`spectral_solve`, 42 groups; without and with the intensity-form laser M-band
# `MBandEmission()`), its 18 channel voltages through the hardware responses
# (`dante1_responses()`) and their unfold (`dante_unfold`, traditional and smooth)
# against the input spectrum's Tr (research notes, "Dante channel responses and the
# unfold (2026-10-24)"). Laser spots split by `spectral_solve`'s default `spot =
# :radiation` (converted x-rays thermal at the local radiation temperature, plus the
# M-band) and, for comparison, by the earlier `:brightness` rule (each element a
# Planckian at its brightness temperature; research notes, "Laser-spot spectrum").
# Also, as calibrated on Dewald's setup (`examples/mband_dewald2020.jl`): the per-cone
# M-band (`nif_mband()`, each element's conversion weighted by the cones' shares of its
# pre-CBET deposition, a proxy: CBET and the bubble change the cone powers, not which
# cone lights an element) with the M-band group albedo `alpha_m` of the gold
# (thermalized), 0.55 with the thermal spots and 0.5 with the brightness rule.
unfold_times = get(opt, "unfold", "off") == "on" ? [5.0, 6.0, 7.0, 8.0, 8.5] : Float64[]
unfolds = Dict{Float64,Any}()
d1g = dante1(channels = dante1_channels(model = :gaussian))   # the Gaussian fallback, for comparison
alpha_m = parse(Float64, get(opt, "alpha_m", "0.55"))
function cone_mband(mb, sol)
    t = sol.model.target
    bs = beams_at(sol.time)
    op = DepositionOperator(t, bs)
    D = [deposit(op, [b.group === g ? b : remake(b; power = 0.0) for b in bs]).element for g in keys(mb.groups)]
    f = collect(values(mb.groups))
    tot = sum(D)
    return [tot[e] > 0 && !Briareus.is_aperture(t, e) && Briareus.element_material(t, e).name in (:Au, :DU) ?
            sum(f[i] * D[i][e] for i in eachindex(D)) / tot[e] * Float64(sol.laser[e]) : 0.0 for e in 1:nelements(t)]
end
function dante_spectral(sol)
    out = Any[]
    E = spectral_edges()
    k = findfirst(==(1800.0), E)
    t = sol.model.target
    mband_albedo(a) = Dict(c.name => [g >= k ? a : NaN for g in 1:length(E)-1]
                           for c in t.components if c.kind === :surface && c.material.name in (:Au, :DU))
    mb = nif_mband()
    PM = cone_mband(mb, sol)
    cases = [(; spot = :brightness, mband = nothing), (; spot = :brightness, mband = MBandEmission()),
             (; spot = :radiation, mband = nothing), (; spot = :radiation, mband = MBandEmission()),
             (; spot = :brightness, mband = mb, mband_power = PM, albedo = mband_albedo(0.5), thermalize = true),
             (; spot = :radiation, mband = mb, mband_power = PM, albedo = mband_albedo(alpha_m), thermalize = true)]
    for kw in cases
        s = spectral_solve(sol; edges = E, kw...)
        d = dante(d1, s)
        gauss = dante(d1g, s).color_temperature
        push!(out, (; d.Tr, d.Tr_open, d.mband_fraction, d.color_temperature, gauss, traditional = d.Tr_unfold,
                    smooth = dante_unfold(d1.channels, d.channels; area = d.area, cos = d.cos, method = :smooth).Tr,
                    misfit = d.unfold.misfit, solver = d.unfold.solver,
                    laser_mband = sum(s.mband .* t.elem_area) / sum(Float64.(sol.laser))))
    end
    return out
end
dante_rec(sol) = (any(≈(sol.time), unfold_times) && (unfolds[sol.time] = dante_spectral(sol)); dante_tr(sol))
runs = Pair{String,Any}[]
t0 = time()
function run!(name, f)
    push!(runs, name => f())
    @printf("%-22s %6.0f s\n", name, time() - t0)
    flush(stdout)
end
if which == "all"
    static = hohlraum(; radius = R0, geometry...)
    model, op = RadiationModel(static), DepositionOperator(static, beams)
    run!("static, no CBET", () -> drive_history(model, op, beams_at, times; lpi = isempty(pre) ? nothing : Tuple(pre),
                                                record = dante_tr))
    run!("static", () -> drive_history(model, op, beams_at, times; lpi = with_pre(cbet), record = dante_tr))
    kw = Keyframes(moving_hohlraum(motion; radius = R0, geometry..., shrink...), keyframes, beams)
    run!("+wall", () -> drive_history(kw, beams_at, times; lpi = with_pre(cbet), record = dante_tr))
end
kl = Keyframes(moving_hohlraum(motion; radius = R0, geometry..., leh_closure = lip, shrink...), keyframes, beams)
which == "all" && run!("+lips", () -> drive_history(kl, beams_at, times; lpi = with_pre(cbet), record = dante_tr))
# `mband=on`: the full model's M-band from each cone's own post-CBET, post-bubble
# deposition at every step (`drive_history`'s `mband`; Dewald's per-cone conversions
# `nif_mband()` on the gold liner of the DU wall), two groups split at 1.8 keV with the
# M-band group albedo `alpha_m` thermalized (the > 1.8 keV flux of a two-group solve is
# that of the 42-group one: the soft groups keep the gray albedo); and, at the same steps,
# the earlier proxy `cone_mband` (the cones weighted by their pre-CBET deposition).
mband_on = get(opt, "mband", "off") == "on"
mb_cones = nif_mband(materials = (:Au, :DU))
mband_spectral = (; edges = [0.0, 1800.0, Inf], thermalize = true,
                  albedo = Dict(c.name => [NaN, alpha_m] for c in kl.reference.components
                                if c.kind === :surface && c.material.name in (:Au, :DU)))
mband_proxy = Dict{Float64,Vector{Float64}}()
function bubble_rec(sol)
    if mband_on
        s = spectral_solve(sol; mband_spectral..., mband = mb_cones, mband_power = cone_mband(mb_cones, sol))
        tg = sol.model.target
        mband_proxy[sol.time] = legendre_moments(band_flux(s, 1800.0), tg, elements(tg, :capsule); lmax = 6)
    end
    return dante_rec(sol)
end
mband_kw = mband_on ? (; mband = mb_cones, spectral = mband_spectral) : (;)
run!("+bubble", () -> drive_history(kl, beams_at, times; lpi = with_pre(cbet, bubble), record = bubble_rec,
                                     ylm = metrics, mband_kw...))
# the full model with the uniform backscatter cut, for Fig. 10 (`backscatter`, above)
fig10_run = backscatter && hasproperty(n221204_fig10(), Symbol(shot))
fig10_run && run!("+bubble, uniform cut",
                  () -> drive_history(kl, beams_at, times; lpi = with_pre(cbet, bubble; p = lpi_pre("uniform")),
                                      record = dante_tr, ylm = metrics))

# PRE 109 Fig. 10, N221204 (HYDRA, post-CBET), digitized (`n221204_fig10()`: ±0.002
# where the curve is visible, up to ±0.015 where N220919's curve covers it).
fig10 = n221204_fig10()
if hasproperty(fig10, Symbol(shot))     # no HYDRA curve for Zylstra's pair
    ref, uref = getfield(fig10, Symbol(shot)), getfield(fig10, Symbol("u", shot[2:end]))
    println("\ncapsule flux P2/P0 ($shot)", fig10_run ? "; Fig. 10's like-for-like run: +bubble, uniform cut" : "")
    println(rpad("t (ns)", 7), join([rpad(n, 22) for (n, _) in runs]), "PRE 109 Fig. 10")
    for t in [2.0, 4.0, 4.5, 5.5, 6.5, 7.0, 8.0, 8.5]
        k = findfirst(≈(t), times)
        i = findfirst(≈(t), fig10.t)
        k === nothing && continue
        println(rpad(@sprintf("%.1f", t), 7), join([rpad(@sprintf("%+.3f", p2(r)[k]), 22) for (_, r) in runs]),
                @sprintf("%+.3f ± %.3f", ref[i], uref[i]))
    end
    # swing shape: minimum after 4.5 ns, the upward zero crossing after it, the value at 8 ns
    function swing(t, p)
        s = findall(k -> 4.5 <= t[k] <= 8.95 && isfinite(p[k]), eachindex(t))
        km = s[argmin(p[s])]
        kz = findfirst(k -> k > km && p[k - 1] < 0 <= p[k], eachindex(p))
        tz = kz === nothing ? NaN : t[kz - 1] - p[kz - 1] * (t[kz] - t[kz - 1]) / (p[kz] - p[kz - 1])
        return (min = p[km], tmin = t[km], zero = tz)
    end
    grid = findall(t -> 0.75 < t < 8.95 && isinteger(round(t * 10, digits = 6)), fig10.t)
    sh = swing(fig10.t[grid], ref[grid])
    @printf("swing (min, at, upward zero): Fig. 10 %+.3f at %.2f ns, %.2f ns\n", sh.min, sh.tmin, sh.zero)
    for (n, r) in runs
        s = swing(times, p2(r))
        ks = [k for k in eachindex(times) if 2.0 <= times[k] <= 8.5 && findfirst(≈(times[k]), fig10.t) !== nothing]
        d = [p2(r)[k] - ref[findfirst(≈(times[k]), fig10.t)] for k in ks]
        @printf("  %-22s %+.3f at %.2f ns, %.2f ns; %+.3f at 8 ns; rms against Fig. 10 (2–8.5 ns) %.3f\n",
                n, s.min, s.tmin, s.zero, p2(r)[findfirst(≈(8.0), times)], sqrt(sum(abs2, d) / length(d)))
    end
end
println("\nradiation temperature (eV): capsule Tr / Dante 1 over the initial LEH area / over the open area; measured Dante ",
        pair ? "(Zylstra PRE 106 Fig. 2)" : "1 (Pak PRE 109 Fig. 1)")
println(rpad("t (ns)", 7), join([rpad(n, 24) for (n, _) in runs]), "measured")
for t in [4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]
    k = findfirst(≈(t), times)
    k === nothing && continue
    i = findfirst(≈(t), measured.t)
    println(rpad(@sprintf("%.1f", t), 7),
            join([rpad(@sprintf("%.1f / %.1f / %.1f", r.Tr[k], r.record[k]...), 24) for (_, r) in runs]),
            @sprintf("%.1f", measured.dante1[i]))
end
if mband_on
    r = last(runs[findfirst(p -> first(p) == "+bubble", runs)])
    mp2 = r.mband_moments[3, :] ./ r.mband_moments[1, :]
    xp2 = [mband_proxy[τ][3] / mband_proxy[τ][1] for τ in times]
    println("\ncapsule M-band (> 1.8 keV) P2/P0, full model: per cone after CBET and the bubble / the earlier",
            " pre-CBET cone weighting / difference; M-band share of the capsule flux; laser M-band over absorbed laser",
            " (Dewald 2020 Fig. 10, Mo sphere, its own shot: mostly > 0, +0.15 at 6 ns)")
    for τ in (2.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0)
        k = findfirst(≈(τ), times)
        k === nothing && continue
        PL = sum(r.laser[c][k] for c in keys(r.laser))
        @printf("  %4.1f ns  %+.3f / %+.3f / %+.4f;  share %.3f;  laser M-band %.3f\n", τ, mp2[k], xp2[k],
                mp2[k] - xp2[k], r.mband_moments[1, k] / r.moments[1, k], sum(r.mband[c][k] for c in keys(r.mband)) / PL)
    end
    w = findall(τ -> 2.0 <= τ <= 9.0, times)
    @printf("  max |difference| over 2–9 ns %.4f (at %.2f ns); mean difference over 4.5–8.5 ns %+.4f\n",
            maximum(abs.(mp2[w] .- xp2[w])), times[w[argmax(abs.(mp2[w] .- xp2[w]))]],
            (ww = findall(τ -> 4.5 <= τ <= 8.5, times); sum(mp2[ww] .- xp2[ww]) / length(ww)))
end
for (n, r) in runs
    isempty(r.lip_radius) && continue
    rl = r.lip_radius[:leh_bottom_lip]
    println("computed lip edge radius, lower LEH (mm), ", n, ": ",
            join([@sprintf("%.2f ns %.3f", t, rl[k]) for t in (2.0, 4.5, 5.5, 8.0, 8.5, 9.0)
                  for k in (findfirst(≈(t), times),) if k !== nothing], ", "))
end
# Implosion metrics (`metrics=on`; ROADMAP 11, research notes "Implosion metrics"), over
# the drive window from the start of the rise (`pulse_rise`) to bang time, against the
# measured N221204 hot spot (P24 Fig. 2, K24 Table II): velocity 36 ± 6 km/s toward
# (θ, φ) = (151°, 54°), P2 0.3 ± 0.1 µm. The velocity is the flat-kernel default
# (100 km/s per % of mode 1, burn-off); the mock pointing is up-down and azimuthally
# symmetric and the hohlraum has no windows, so its mode 1 is the mesh's and the beams'
# residual. The flux P2 is compared as a time average with Fig. 10's (HYDRA), and as a
# hot-spot P2 shift from it at the −9.2 µm per % of flux P2 over 4.5–6.5 ns transfer of
# the research notes ("Callahan 2018 and Izumi 2018", a sensitivity from another
# platform, *L*); Briareus ships no hot-spot P2 kernel.
if metrics
    rise = pulse_rise((pt, total))
    println("\nimplosion metrics, drive window $(round(rise, digits = 2))–$(bang) ns (rise to bang)")
    for (n, r) in runs
        r.ylm === nothing && continue
        v = hotspot_velocity(r; window = (rise, bang))
        m = implosion_metrics(r; window = (rise, bang))
        sp = sqrt(sum(abs2, v))
        θ, φ = sp > 0 ? (acosd(v[3] / sp), mod(atand(v[2], v[1]), 360)) : (NaN, NaN)
        @printf("  %-22s hot-spot velocity %.2g km/s toward (%.0f°, %.0f°) (measured 36 ± 6 at (151°, 54°)); no-α yield factor %.4f\n",
                n, sp, θ, φ, m.yield_factor)
        sw = swing_metrics(r; window = (rise, bang))
        @printf("  %-22s flux P2/P0 over the window: mean %+.4f, min %+.3f, max %+.3f; P4/P0 mean %+.4f\n",
                "", sw.mean, sw.min, sw.max, swing_metrics(r; l = 4, window = (rise, bang)).mean)
        k = findfirst(≈(bang), times)
        @printf("  %-22s at bang: capsule Tr %.1f eV, Dante %.1f (open area %.1f)", "", r.Tr[k], r.record[k]...)
        i = findfirst(τ -> abs(τ - bang) < 0.026, measured.t)
        i === nothing || isnan(measured.dante1[i]) || @printf(", measured Dante %.1f", measured.dante1[i])
        println()
        if hasproperty(n221204_fig10(), Symbol(shot))
            f10 = n221204_fig10()
            hyd = getfield(f10, Symbol(shot))
            for w in ((rise, 8.9), (4.5, 6.5))
                kk = findall(τ -> w[1] - 1e-9 <= τ <= w[2] + 1e-9, f10.t)
                th, yh = f10.t[kk], hyd[kk]
                mh = sum((th[j+1] - th[j]) * (yh[j] + yh[j+1]) / 2 for j in 1:length(th)-1) / (th[end] - th[1])
                mb = swing_metrics(r; window = w).mean
                @printf("  %-22s %.2f–%.2f ns: mean flux P2/P0 %+.4f, Fig. 10 %+.4f; difference %+.4f → hot-spot P2 %+.1f µm from HYDRA's\n",
                        "", w..., mb, mh, mb - mh, -9.2 * 100 * (mb - mh))
            end
        end
    end
end
println("\ncapsule Tr peak (eV; drive Trad 313 eV, AS24)")
for (n, r) in runs
    k = argmax(r.Tr)
    @printf("  %-22s %.1f at %.2f ns\n", n, r.Tr[k], times[k])
end
println("\nDante 1 Tr peak (eV; measured 313 ± 2, Rubery PRL 132, 065104; Pak PRE 109 Fig. 1: 311.9 at 8.30 ns)")
for (n, r) in runs
    dt, dto = first.(r.record), last.(r.record)
    k = argmax(dt)
    @printf("  %-22s %.1f at %.2f ns (over the open area %.1f); capsule Tr there %.1f (Dante − capsule %+.1f eV)\n",
            n, dt[k], times[k], dto[k], r.Tr[k], dt[k] - r.Tr[k])
    k = findfirst(≈(8.5), times)
    k === nothing || @printf("  %-22s at 8.5 ns: Dante %.1f (open area %.1f), capsule %.1f\n", "", dt[k], dto[k], r.Tr[k])
end
if !isempty(unfolds)
    println("\nDante 1 unfolded (full model; eV over the initial LEH area): input spectrum Tr / traditional unfold /",
            " smooth unfold; color temperature (hardware, Gaussian channels); M-band share > 1.8 keV;",
            " worst channel misfit (traditional)")
    labels = ("brightness", "brightness + M-band", "spots thermal", "spots thermal + M-band",
              "brightness, cones α_M 0.5", @sprintf("thermal, cones α_M %.2g", alpha_m))
    for t in sort(collect(keys(unfolds))), (u, lab) in zip(unfolds[t], labels)
        @printf("  %.1f ns %-25s %.1f (open %.1f) / %.1f / %.1f (%+.2f / %+.2f); Tc %.1f, %.1f; M-band %.3f (laser → M-band %.3f); misfit %.1e (%s)\n",
                t, lab, u.Tr, u.Tr_open, u.traditional, u.smooth, u.traditional - u.Tr, u.smooth - u.Tr,
                u.color_temperature, u.gauss, u.mband_fraction, u.laser_mband, u.misfit, u.solver)
    end
end

# Energy budget over the solved times (0.5–9.35 ns), against PRE 109 App. 3's HYDRA
# split of N221204's incident laser energy: wall 54.1 % (Au lining 4.4, DU 49.7), LEH
# 25.2 %, He gas 4.5 % (not modeled here: the bubble's absorbed light is conducted to
# the wall, so compare the wall with wall + gas, 58.6 %), capsule 16.2 %. Fractions of
# the coupled energy (after backscatter and degradation) and of the delivered.
println("\nenergy budget (kJ; % of coupled / of delivered)")
surfaces(b) = [n for n in keys(b.heating) if n !== :capsule]
for (n, r) in runs
    b = energy_budget(r)
    wallh = sum(b.heating[k] for k in surfaces(b))
    leh = sum(values(b.escaped)) + sum(values(b.transmitted); init = 0.0)
    pc(x) = @sprintf("%7.1f (%4.1f / %4.1f)", x, 100x / b.coupled, 100x / b.setup)
    @printf("  %-22s delivered %.0f, backscatter+degradation %.0f, coupled %.0f\n", n, b.setup, b.backscattered, b.coupled)
    println("    wall heating      ", pc(wallh), "   [laser not converted ", pc(sum(b.laser[k] - b.source[k] for k in surfaces(b))), "]")
    println("    LEH x-rays        ", pc(leh))
    println("    capsule (E_cap)   ", pc(b.absorbed[:capsule]))
    println("    laser out the LEH ", pc(b.unabsorbed), "   x-rays lost ", pc(b.lost))
    println("    wall re-emitted   ", pc(sum(values(b.released); init = 0.0)), "   (stored energy given back; wall heating is net of it)")
    @printf("    mean conversion %.3f\n", sum(values(b.source)) / sum(values(b.laser)))
end

# CBET and bubble diagnostics: post-transfer cone fraction (PRE 109 Fig. 9b: 0.65 at
# 5.6 ns, 0.55 at 7 ns, 0.52 at 8 ns) and the inner power the bubble takes.
println("\nt (ns)  CF in  CF post-CBET  bubble x (mm)  inner23/30 absorbed  ridge g")
isin(b) = b.group in (:inner23, :inner30)
for t in (2.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 8.5)
    bs = beams_at(t)
    ec = apply_lpi(with_pre(cbet), bs, t)
    cf(x) = sum(b.power for b in x if isin(b)) / sum(b.power for b in x)
    f = bubble_fractions(bubble, ec, t)
    a(g) = (i = findall(b -> b.group == g, ec); sum(f.absorbed[i] .* [ec[j].power for j in i]) /
                                                 sum(ec[j].power for j in i))
    @printf("%4.1f    %.3f  %.3f         %.3f          %.3f / %.3f        %.2f\n",
            t, cf(bs), cf(ec), f.depth, a(:inner23), a(:inner30), f.ridge)
end

tag = join(("$k-$(opt[k])" for k in given), "_")
out = "$(shot)_$(transfer)_$(h)_$(capsule_h)$(isempty(tag) ? "" : "_" * tag).csv"
header = ["time"; ["P2/P0 $n" for (n, _) in runs]; ["P4/P0 $n" for (n, _) in runs]; ["Tr $n" for (n, _) in runs];
          ["Dante $n" for (n, _) in runs]; ["Dante open $n" for (n, _) in runs]]
writedlm(out, vcat(permutedims(header),
                   hcat(times, [p2(r) for (_, r) in runs]..., [p4(r) for (_, r) in runs]...,
                        [r.Tr for (_, r) in runs]..., [first.(r.record) for (_, r) in runs]...,
                        [last.(r.record) for (_, r) in runs]...)), ',')
println("\nwrote ", out)
