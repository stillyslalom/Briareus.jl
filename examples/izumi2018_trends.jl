# Izumi 2018 (LLNL-PROC-751038) single-parameter SymCap changes against Briareus's flux
# P2 with `BubbleAbsorption`, turned into hot-spot P2 with Callahan et al. 2018's 2.7 µm
# per 1 % cone fraction over Briareus's own d(flux P2/P0)/dCF. A trend check, not blind
# (the bubble model was fitted to Izumi's data). Research notes, "Young 2020 validation".
#
#     julia --project=. -t auto examples/izumi2018_trends.jl     # ~5 min on 4 threads
using Briareus, LinearAlgebra, Printf
example = Dict(:outer50 => (5.01, 0.0), :outer44 => (5.20, 0.0), :inner30 => (3.63, -0.60), :inner23 => (3.09, -1.09))
p = GroupPointing(nif_quads(power = 8.0); fields = (:z, :r))
quads = apply(p, reduce(vcat, [[example[k]...] for k in p.keys]))
P = young2020_pulses(); side(b) = startswith(String(b.group), "inner") ? :inner : :outer
times = collect(0.5:0.25:7.5)
pk = (times .>= 4.5) .& (times .<= 6.5)
function run(; rcap = 0.908, fill = 0.3, energy = 1.0, finner = 1.0, bubble = true, cache = Dict())
    m, op = get!(cache, rcap) do
        t = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = rcap,
                     wall = gold(albedo = 0.8, conversion = 0.9), h = 0.45, capsule_h = 0.15)
        mm = RadiationModel(t; closed = true); (mm, DepositionOperator(t, quads))
    end
    pul = (; inner = PulseShape(P.inner.t, 4finner .* P.inner.p), outer = PulseShape(P.outer.t, 4 .* P.outer.p))
    lpi = bubble ? BubbleAbsorption(; radius = 2.875, turn_on = 0.25, fill, energy, capsule_radius = rcap) : nothing
    h = drive_history(m, op, t -> with_pulses(quads, pul, t; by = side), times; lpi, lmax = 2)
    p2 = h.moments[3, :] ./ h.moments[1, :]
    return 100 * sum(p2[pk]) / count(pk)
end
cache = Dict()
base = run(; cache)
cf0 = let i = sum(P.inner.(times[pk])), o = sum(P.outer.(times[pk])); i / (i + 2o) end
cfp = let i = 1.05 * sum(P.inner.(times[pk])), o = sum(P.outer.(times[pk])); i / (i + 2o) end
s = (run(; finner = 1.05, cache) - base) / (100 * (cfp - cf0))     # % flux P2 per % CF
k = 2.7 / s                                                      # hot-spot µm per % flux P2
@printf("base peak flux P2/P0 %.2f %%; CF %.3f; d(P2 %%)/d(CF %%) = %.3f; transfer 2.7 µm/%%CF -> %.2f µm per %% flux P2\n", base, cf0, s, k)
for (name, kw, meas) in (("1.5x picket", (; energy = 1.5), "-21 to -24"), ("1.18x capsule", (; rcap = 1.064), "-21 to -24"),
                         ("1.18x capsule, no bubble change", (; rcap = 1.064, bubble = :novf), ""), ("0.15 mg/cc fill", (; fill = 0.15), "+10 to +15"))
    if get(kw, :bubble, true) === :novf
        d = run(; rcap = 1.064, bubble = false, cache) - run(; bubble = false, cache)
        @printf("%-32s flux P2 %+6.2f pp (view factors only) -> %+6.1f µm\n", name, d, k * d)
        continue
    end
    d = run(; kw..., cache) - base
    @printf("%-32s flux P2 %+6.2f pp -> hot-spot %+6.1f µm  [Izumi 2018 core: %s µm]\n", name, d, k * d, meas)
end
