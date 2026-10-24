# NIF Dante 1 channel responses from their hardware: refits the XRD and mirror
# parameters of `src/dante_response.jl` to the curves digitized from Rubery et al., RSI
# 2023 (`data/rubery2023_dante.csv`), checks channel 2 against its published calibration,
# lists the channels' bands and unfolds Planckians. Research notes, "Dante channel
# responses and the unfold (2026-10-24)".
#
#   julia --project=. examples/dante_channels.jl

using Briareus, Optim, Printf
using Briareus: _rubery_dante, _xrd_model, _response, _DANTE_GRID, _planck_shape, _XRD_SURFACE,
                _XRD_FIT, _SIO2_ANGLE_SCALE, _SIO2_ROUGHNESS, SIGMA_SB

d = _rubery_dante()
rms(x) = sqrt(sum(abs2, x) / length(x))

# --- SiO2 mirror: effective grazing angle and roughness (Fig. 3a, "3.5°") ----------------
Em, Rm = d["sio2_mirror"]
k = Em .> 60
fm(p) = rms([mirror_reflectivity(:SiO2, p[1], e; roughness = abs(p[2])) for e in Em[k]] .- Rm[k])
r = optimize(fm, [2.2, 1.0], NelderMead())
θ, σ = Optim.minimizer(r)[1], abs(Optim.minimizer(r)[2])
@printf("SiO2 \"3.5°\" mirror: effective angle %.3f° (scale %.4f), roughness %.2f nm, rms %.3f (src: %.4f, %.2f nm)\n",
        θ, θ / 3.5, σ, Optim.minimum(r), _SIO2_ANGLE_SCALE, _SIO2_ROUGHNESS)
@printf("  at the nominal 3.5°, no roughness: rms %.3f\n", fm([3.5, 0.0]))

# --- Cr XRD (Fig. 3b): K, L, surface C and O, γ ---------------------------------------------
Ec, Sc = d["cr_xrd"]
k = Ec .> 75
bad(p) = !(0 <= p[3] <= 0.9 && 0 <= p[4] <= 0.9 && p[3] + p[4] <= 0.95)
fc(p) = bad(p) ? 1e9 : rms(log.(_xrd_model.(:Cr, Ec[k], 1e6exp(p[1]), p[5]; L = exp(p[2]), wC = p[3], wO = p[4]) ./ 1e6) .- log.(Sc[k]))
best = nothing
for L0 in (3.0, 10.0, 30.0), wc in (0.05, 0.3, 0.6), g in (-0.5, 0.0, 0.5)
    o = optimize(fc, [log(5.0), log(L0), wc, 0.1, g], NelderMead(), Optim.Options(iterations = 6000))
    global best = best === nothing || Optim.minimum(o) < Optim.minimum(best) ? o : best
end
p = Optim.minimizer(best)
@printf("Cr XRD: K %.4g MV/GW, L %.4g nm, w_C %.4f, w_O %.4f, γ %.4f; rms log %.3f\n",
        exp(p[1]), exp(p[2]), p[3], p[4], p[5], Optim.minimum(best))
@printf("  src: K %.2f MV/GW, L %.1f nm, w_C %.4f, w_O %.4f, γ %.4f\n", _XRD_FIT[:Cr].K / 1e6,
        _XRD_SURFACE.L, _XRD_SURFACE.wC, _XRD_SURFACE.wO, _XRD_FIT[:Cr].γ)

# --- Al XRD (Fig. 8b): K, γ on Cr's surface layer (the Al K edge excluded) ---------------
Ea, Sa = d["al_xrd"]
k = .!(1540 .< Ea .< 1750)
fa(q) = rms(log.(_xrd_model.(:Al, Ea[k], 1e6exp(q[1]), q[2]; L = exp(p[2]), wC = p[3], wO = p[4]) ./ 1e6) .- log.(Sa[k]))
o = optimize(fa, [log(40.0), 0.0], NelderMead())
q = Optim.minimizer(o)
@printf("Al XRD: K %.4g MV/GW, γ %.4f; rms log %.3f (src: %.2f, %.4f)\n", exp(q[1]), q[2], Optim.minimum(o),
        _XRD_FIT[:Al].K / 1e6, _XRD_FIT[:Al].γ)

# --- channel 2 (2–4 keV multilayer) against its calibration and Fig. 8b -----------------
ch = dante1_responses()
E = _DANTE_GRID
V = ch[2].solid_angle * sum(_response(ch[2], E[i]) * 0.5 * (E[i+1] - E[i]) for i in 1:length(E)-1 if 2000 <= E[i] < 4000)
@printf("channel 2: %.3f TW/sr per volt for power flat over 2–4 keV (Rubery 2023: 0.97)\n", 1 / V)
Eg, g = d["mlm_system"]
k = 2000 .< Eg .< 4000
rr = [_response(ch[2], e) / 1e3 for e in Eg[k]]
@printf("  response / Fig. 8b over 2–4 keV: mean %.3f, rms %.3f\n", sum(rr) / sum(g[k]), rms(rr ./ g[k] .- 1))

# --- the channels ---------------------------------------------------------------------------
println("\nchannel  nominal  peak  half-maximum band (eV)  peak ΩR (V per GW/sr/eV)")
for c in ch
    rc = c.R
    m = argmax(rc)
    i = m; while i > 1 && rc[i] > rc[m] / 2; i -= 1; end
    j = m; while j < length(rc) && rc[j] > rc[m] / 2; j += 1; end
    @printf("%5d  %7.0f  %5.0f  %6.0f – %6.0f        %.3g\n", c.number, c.energy, c.E[m], c.E[i], c.E[j],
            c.solid_angle * rc[m])
end

# --- Planckians through the unfold ------------------------------------------------------------
A, cs = π * 1.55^2, cosd(37)
w = [i == 1 ? (E[2] - E[1]) / 2 : i == length(E) ? (E[end] - E[end-1]) / 2 : (E[i+1] - E[i-1]) / 2 for i in eachindex(E)]
println("\nPlanckian T (eV): exact Tr over 50 eV – 20 keV / traditional / smooth unfold; channel 1 and 3 volts")
for T in (100.0, 150.0, 200.0, 250.0, 300.0, 350.0)
    S = 1e3 * A * cs / π * SIGMA_SB * T^4 .* _planck_shape.(E, T)
    Vp = [c.solid_angle * sum(w .* c.R .* S) for c in ch]
    ut = dante_unfold(ch, Vp; area = A, cos = cs)
    us = dante_unfold(ch, Vp; area = A, cos = cs, method = :smooth)
    @printf("%5.0f   %.3f / %.3f (%s, misfit %.1e) / %.3f    %.3g V, %.3g V\n", T,
            T * planck_fraction(T, 50.0, 20_000.0)^0.25, ut.Tr, ut.solver, ut.misfit, us.Tr, Vp[1], Vp[3])
end
