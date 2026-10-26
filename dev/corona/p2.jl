# What a spherical corona drops: the OMEGA 40-beam PDD case (public inputs; review §4.1:
# rings 21.4°/42.0°/58.9° of `omega_beams` shifted 90/150/150 µm normal to the beam toward
# the equator, SG4 r₀ = 308 µm, n = 3.66, f/6.7) in an isothermal exponential corona
# (r_c = 0.40 mm, L = 0.08 mm, η⊥ = 1 - e^{-2.4}, truncated at r_c + 14L; an engineering
# snapshot) whose critical radius carries a P2, n_e/n_c = exp(-(r - r_c(1 + ε P₂))/L):
# 3-D ray tracing (Corona.jl) against the package's spherical `CoronaOperator`, which
# returns the ε = 0 result for every ε. Numbers: research notes, "Corona model prototype".
#   julia --project=. -t 8 dev/corona/p2.jl [ε] [hlat] [hτ]
include(joinpath(@__DIR__, "Corona.jl"))
using .Corona, Briareus, Printf
const C = Corona

ε = parse(Float64, get(ARGS, 1, "0.04"))
hlat = parse(Float64, get(ARGS, 2, "0.04"))
hτ = parse(Float64, get(ARGS, 3, "0.004"))
rc, L = 0.40, 0.08
νc = 2.4 / (8 / 3 * L)
ring = Dict(:cone1 => 0.090, :cone2 => 0.150, :cone3 => 0.150)
base = omega_beams(; cones = (:cone1, :cone2, :cone3), spot = SuperGaussianSpot(0.308; order = 3.66),
                   power = 1.0, fnumber = 6.7)
beams = [remake(b; aim = ring[b.group] * sign(b.port[3]) * Briareus.port_basis(b.port)[1]) for b in base]

op = CoronaOperator(CoronaProfile(; critical_radius = rc, scale_length = L, absorption = 1 - exp(-2.4)),
                    beams; ncone = 7)
d = deposit(op, beams)
G = C.AngGrid(720, 30)
E, inc, nr, _ = C.trace_beams(C.P2ExpProfile(rc, L, νc, rc + 14L, ε), beams, G; ncone = 7, hlat, hτ)
c = C.ylm_coeffs(E, G)
σb, σc = mode_power(ylm_matrix(d.ylm, 30)), mode_power(ylm_matrix(c, 30))
@printf("spherical operator: absorbed %.5f, σ (l ≤ 30) %.4f %%, σ₂ %.4f %%, c₂₀/c₀₀ %+.5f\n",
        sum(d.per_beam) / d.delivered, 100mode_rms(d.ylm, 30), 100σb[3], d.ylm[mode_index(2, 0)] / d.ylm[1])
@printf("3-D trace, ε = %+.3f (%d rays): absorbed %.5f, σ %.4f %%, σ₂ %.4f %%, c₂₀/c₀₀ %+.5f\n", ε, nr,
        sum(E) / inc, 100mode_rms(c, 30), 100σc[3], c[mode_index(2, 0)] / c[1])
