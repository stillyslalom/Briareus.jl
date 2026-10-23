# N221204-class time-dependent benchmark: capsule flux P2/P0(t) of a HYBRID-E hohlraum
# with CBET, a moving wall, closing entrance holes and inner-beam bubble absorption, each
# piece added in turn, against the HYDRA curve of Kritcher et al., PRE 109, 025204 (2024)
# Fig. 10. Inputs: reference/litreview/n221204_benchmark.md §5; results: research notes,
# "N221204 time-dependent benchmark (2026-10-23)".
#
#   julia --project=. -t auto examples/n221204.jl [h capsule_h [runs [transfer [table]]]]
#
# `h`, `capsule_h`: mesh (default 0.28 / 0.16, the standard mesh); `runs`: `all`
# (default) or `full` (only the full model, e.g. for a fine-mesh confirmation);
# `transfer`: `reduced` (default: the reduced CBET model at the deck's settings) or
# `fig9b` (the outer-to-inner transfer that gives PRE 109 Fig. 9b's post-transfer cone
# fraction, HYDRA's own CBET result for this shot). Writes
# `n221204_<transfer>_<h>_<capsule_h>.csv` (P2/P0, P4/P0 and capsule Tr per run) to the
# current directory.

using Briareus, Printf, DelimitedFiles

h = parse(Float64, get(ARGS, 1, "0.28"))
capsule_h = parse(Float64, get(ARGS, 2, "0.16"))
which = get(ARGS, 3, "all")
transfer = get(ARGS, 4, "reduced")

# Geometry (mm; n221204_benchmark.md §1, PRE 109 Fig. 1a): 11.24 × 6.4 mm Au-lined DU
# hohlraum, 3.1 mm LEHs, 1.135 mm capsule. Unpublished, so repository defaults: square
# corners, no diagnostic windows, a gold wall with the Au wall-loss albedo, the default
# ablator. The capsule stays at its initial radius (HYDRA's ablation front is near
# 0.92 mm at 6 ns).
geometry = (length = 11.24, leh_radius = 1.55, capsule_radius = 1.135, h, capsule_h,
            wall = gold(albedo = WallLossAlbedo()))
R0 = 3.2

# Laser: 192 beams at the public mock HYBRID-E pointing; Δλ 2.75 Å (1ω) on the inner
# cones, 0.917 Å at 3ω (the 23.5°/30° and 44.5°/50° splits are unpublished: none).
beams = [remake(b; dlambda = b.group in (:inner23, :inner30) ? 2.75 / 3 : 0.0)
         for b in hybride_mock_pointing(nif_beams())]

# Delivered outer (128 beams) and inner (64 beams) power, DN24 Fig. 1 read by eye (TW).
pt = [0.0, 0.4, 0.8, 1.3, 2.0, 2.5, 3.3, 3.9, 4.5, 5.0, 5.5, 8.0, 8.3, 8.6, 8.9, 9.1, 9.25, 9.35]
outer = [0, 2, 50, 50, 25, 26, 60, 65, 150, 270, 296, 294, 250, 110, 40, 0, 0, 0]
inner = [0, 2, 4, 4, 4.5, 5, 7, 9, 25, 90, 146, 146, 146, 146, 146, 146, 140, 0]
total = Float64.(outer .+ inner)
pulses = cone_fraction_pulses(beams, pt, total, inner ./ max.(total, eps()))
beams_at = τ -> with_pulses(beams, pulses, τ)
@printf("laser energy %.3f MJ, peak %.0f TW\n",
        sum(pulse_energy(pulses[b.group]) for b in beams) / 1e3, maximum(total))

# CBET: the reduced model at the deck's settings, δn/n clamp 0.002 until the start of
# the rise (≈4.0 ns) and 0.0075 after (PRE 109 App. 1); the other plasma parameters are
# the model's defaults (repository defaults, not N221204's).
reduced = CBET(leh_z = geometry.length / 2, clamp = ([0.0, 3.9, 4.1, 20.0], [0.002, 0.002, 0.0075, 0.0075]))
# The alternative: move the fraction of the outer power that turns the input cone
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
cbet = transfer == "fig9b" ? fig9b : reduced

# Hohlraum dynamics, all at repository defaults: the wall moves in at 0.10 mm/ns from
# the outer turn-on (0.4 ns); the LEH lips (0.84 r₀) close at 5.0 ns over 0.5 ns (the
# example closure of the research notes, not measured for N221204); the bubble model
# with its defaults, launched by the same turn-on, footprint following the wall.
turn_on = 0.4
motion = WallMotion(start = turn_on)
lip = GapClosure(closure = 5.0, width = 0.5, after = WallLossAlbedo())
bubble = BubbleAbsorption(; radius = R0, turn_on, wall = motion)
# Every step is quasi-static (no burnthrough fronts here), so `table` (5th argument)
# solves only the times of the comparison table, for a fine-mesh check.
times = get(ARGS, 5, "") == "table" ? [2.0, 4.0, 4.5, 5.5, 6.5, 7.0, 8.0, 8.5] : collect(0.5:0.25:9.0)
keyframes = range(turn_on, 9.0, 11)

p2(r) = r.moments[3, :] ./ r.moments[1, :]
p4(r) = r.moments[5, :] ./ r.moments[1, :]
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
    run!("static, no CBET", () -> drive_history(model, op, beams_at, times))
    run!("static", () -> drive_history(model, op, beams_at, times; lpi = cbet))
    kw = Keyframes(moving_hohlraum(motion; radius = R0, geometry...), keyframes, beams)
    run!("+wall", () -> drive_history(kw, beams_at, times; lpi = cbet))
end
kl = Keyframes(moving_hohlraum(motion; radius = R0, geometry..., leh_closure = lip), keyframes, beams)
which == "all" && run!("+lips", () -> drive_history(kl, beams_at, times; lpi = cbet))
run!("+bubble", () -> drive_history(kl, beams_at, times; lpi = (cbet, bubble)))

# PRE 109 Fig. 10, N221204 (HYDRA, post-CBET), digitized off the figure (±0.001 where
# the curve is visible, ±0.006–0.024 where covered, 7.3–8.9 ns).
fig10 = Dict(2.0 => -0.002, 4.0 => -0.013, 4.5 => -0.091, 5.5 => -0.168, 6.5 => -0.057,
             7.0 => 0.039, 8.0 => 0.090, 8.5 => 0.117)
println("\ncapsule flux P2/P0")
println(rpad("t (ns)", 7), join([rpad(n, 17) for (n, _) in runs]), "PRE 109 Fig. 10")
for t in sort(collect(keys(fig10)))
    k = findfirst(≈(t), times)
    println(rpad(@sprintf("%.1f", t), 7), join([rpad(@sprintf("%+.3f", p2(r)[k]), 17) for (_, r) in runs]),
            @sprintf("%+.3f", fig10[t]))
end
println("\ncapsule Tr peak (eV; drive Trad 313 eV, AS24)")
for (n, r) in runs
    k = argmax(r.Tr)
    @printf("  %-22s %.1f at %.2f ns\n", n, r.Tr[k], times[k])
end

# CBET and bubble diagnostics: post-transfer cone fraction (PRE 109 Fig. 9b: 0.65 at
# 5.6 ns, 0.55 at 7 ns, 0.52 at 8 ns) and the inner power the bubble takes.
println("\nt (ns)  CF in  CF post-CBET  bubble x (mm)  inner23/30 absorbed  ridge g")
isin(b) = b.group in (:inner23, :inner30)
for t in (2.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 8.5)
    bs = beams_at(t)
    ec = apply_lpi(cbet, bs, t)
    cf(x) = sum(b.power for b in x if isin(b)) / sum(b.power for b in x)
    f = bubble_fractions(bubble, ec, t)
    a(g) = (i = findall(b -> b.group == g, ec); sum(f.absorbed[i] .* [ec[j].power for j in i]) /
                                                 sum(ec[j].power for j in i))
    @printf("%4.1f    %.3f  %.3f         %.3f          %.3f / %.3f        %.2f\n",
            t, cf(bs), cf(ec), f.depth, a(:inner23), a(:inner30), f.ridge)
end

out = "n221204_$(transfer)_$(h)_$(capsule_h).csv"
header = ["time"; ["P2/P0 $n" for (n, _) in runs]; ["P4/P0 $n" for (n, _) in runs]; ["Tr $n" for (n, _) in runs]]
writedlm(out, vcat(permutedims(header),
                   hcat(times, [p2(r) for (_, r) in runs]..., [p4(r) for (_, r) in runs]...,
                        [r.Tr for (_, r) in runs]...)), ',')
println("\nwrote ", out)
