# Gold M-band preheat asymmetry at the capsule, against Dewald et al., "Direct
# measurement of hohlraum M-band preheat asymmetry at the ICF capsule", PoP 27 (2020)
# (LLNL-JRNL-813556). Multigroup radiosity (`spectral_solve`) with the per-cone
# laser-produced M-band (`nif_mband`) on their setup, a static solve at peak power.
# Research notes, "Spectral radiosity, M-band and a synthetic Dante (2026-10-24)".
#
#   julia --project=. -t auto examples/mband_dewald2020.jl [h capsule_h] [spot=radiation|wall|brightness]
#
# Setup (Dewald 2020 §II–V): 10.1 × 5.75 mm Au hohlraum, 3.37 mm LEHs, 360 TW peak with
# inner/outer power 1/2 (cone fraction 1/3), the wall ~0.2 mm in, the HDC ablation front
# at 0.815 mm (4.35 ns). Their pointing is not published: the mock HYBRID-E pointing
# (`hybride_mock_pointing`) scaled to this hohlraum, its outer aims moved along the axis
# until the outer-beam ring sits at the measured z = ±3.3 mm (time-integrated hard x-ray
# image, their Fig. 2; θ_o ≈ 40°, 4.4 mm from the centre). Square corners, no window
# or Mo sphere. At 5.5 ns, mid-plateau (4.9–6.0 ns): `nif_conversion()` there (0.86) and
# gold's `WallLossAlbedo` at the wall's effective age, 0.8 ns (*L*: `drive_history`'s
# wall history after the rise at ~4.4 ns; N221204's measured Dante gives 0.76 ns at the
# start of its plateau), not at the 5.5 ns since the pulse start (research notes,
# "Energetics defaults and the capsule share (2026-10-25)").
#
# Their numbers: Dante sees 12 % of the flux above 1.8 keV at the 290 eV peak (15 % in
# their Fig. 3); the wall M-band P2/P0 is +0.53 (Mo-sphere P2 0.15 / smoothing 0.28),
# which with outer spots at 40° requires F_M ≈ 2 (recirculating / laser-produced M-band);
# at the HDC ablation front the smoothing is 0.4 (their Fig. 6), so ≈ +0.21 there, and
# +0.28 at its inner surface (volumetric absorber, their Fig. 14).
#
# The M-band group albedo α_M (> 1.8 keV, Au wall and caps; the rest of the gray albedo
# re-emitted thermally, so totals are the gray solution's) is the one free number: it is
# scanned, and the α_M with Dewald's F_M ≈ 2 is the calibrated default of the notes.
#
# `spot`: the spectrum of the laser spots' converted x-rays (`spectral_solve`'s `spot`):
# thermal at the local radiation temperature plus the M-band (`radiation`, default), at
# the wall temperature (`wall`), or the whole radiosity as a Planckian at its brightness
# temperature (`brightness`, the rule of the first calibration). Research notes,
# "Laser-spot spectrum (2026-10-24)".

using Briareus, Printf

pos = filter(a -> !occursin('=', a), ARGS)
opt = Dict(String.(split(a, '=', limit = 2)) for a in ARGS if occursin('=', a))
h = parse(Float64, get(pos, 1, "0.3"))
capsule_h = parse(Float64, get(pos, 2, "0.15"))
spot = Symbol(get(opt, "spot", "radiation"))
L, R = 10.1, 2.875 - 0.2
t = hohlraum(; length = L, radius = R, leh_radius = 3.37 / 2, capsule_radius = 0.815, h, capsule_h,
             wall = gold(albedo = WallLossAlbedo(), conversion = conversion_efficiency(nif_conversion(), 5.5)))
t_eff = 0.8
wall, cap = elements(t, :wall), elements(t, :capsule)
isinner(b) = b.group in (:inner23, :inner30)
zw = [t.elem_centroid[e][3] for e in wall]

# pointing: the mock HYBRID-E table scaled to 10.1 × 5.75 mm, outer aims shifted by dz
function pointing(dz)
    table = Dict(k => (z = v.z * L / 11.24 + (k in (:outer44, :outer50) ? dz : 0.0),
                       r = v.r * 2.875 / 3.2, s = v.s) for (k, v) in Briareus.HYBRIDE_MOCK_POINTING)
    return hybride_mock_pointing(nif_beams(); table)
end
function ring(beams)
    op = DepositionOperator(t, beams)
    d = deposit(op, [remake(b; power = isinner(b) ? 0.0 : b.power) for b in beams])
    w, up = d.element[wall], zw .> 0
    return sum(w[up] .* zw[up]) / sum(w[up]), op
end
dz, beams, op, zo = 0.0, nothing, nothing, 0.0
for _ in 1:4
    global beams = pointing(dz)
    global zo, op = ring(beams)
    abs(zo - 3.3) < 0.02 && break
    global dz += 3.3 - zo
end
@printf("outer aims %+.2f mm along the axis: ring at z = %.2f mm (3.3), θ_o %.1f° (40), %.2f mm from TCC (4.4)\n",
        dz, zo, atand(R, zo), hypot(R, zo))

m = RadiationModel(t)
E = spectral_edges()
k = findfirst(==(1800.0), E)
au = (:wall, :cap_top, :cap_bottom)
mband_albedo(a) = Dict(c => [g >= k ? a : NaN for g in 1:length(E)-1] for c in au)
boundary = vcat((elements(t, c) for c in (:wall, :cap_top, :cap_bottom, :leh_top, :leh_bottom))...)
# P2/P0 of a radiosity over the boundary as seen from the centre (LEHs dark), and of the
# capsule's incident flux
wallP2(x) = (a = legendre_moments(x, t, boundary; lmax = 2, center = zero(Vec3)); a[3] / a[1])
capP(x) = (a = legendre_moments(x, t, cap; lmax = 4); (a[3] / a[1], a[5] / a[1]))
mean_cap(x) = sum(x[cap] .* t.elem_area[cap]) / sum(t.elem_area[cap])
d1 = DanteView(143, 274)
mb = nif_mband()
for cf in (1 / 3, 0.25)
    bs = [remake(b; power = isinner(b) ? cf * 360 / 64 : (1 - cf) * 360 / 128) for b in beams]
    dep = deposit(op, bs)
    PM = mband_power(mb, op, bs)
    sol = solve(m, dep; time = t_eff)
    g2, g4 = capP(sol.H)
    @printf("\ncone fraction %.3f: laser M-band %.1f TW (%.3f of the absorbed laser); capsule Tr %.1f eV, Dante-1 Tr %.1f eV\n",
            cf, sum(PM), sum(PM) / sum(dep.element), radiation_temperature(mean_cap(sol.H)), dante(d1, sol).Tr)
    @printf("  gray: capsule P2/P0 %+.3f, P4/P0 %+.3f; boundary P2/P0 %+.3f\n", g2, g4, wallP2(sol.B))
    @printf("  spots at :%s\n", spot)
    println("  α_M      F_M   capsule M-band: fraction  P2/P0   P4/P0  (first flight)  boundary P2/P0   Dante-1 > 1.8 keV, Tc")
    for (name, kw) in (("gray", (;)), ("0.7", (; albedo = mband_albedo(0.7), thermalize = true)),
                       ("0.6", (; albedo = mband_albedo(0.6), thermalize = true)),
                       ("0.55", (; albedo = mband_albedo(0.55), thermalize = true)),
                       ("0.5", (; albedo = mband_albedo(0.5), thermalize = true)),
                       ("0.3", (; albedo = mband_albedo(0.3), thermalize = true)),
                       ("0", (; albedo = mband_albedo(0.0), thermalize = true)))
        s = spectral_solve(sol; edges = E, mband = mb, mband_power = PM, spot, kw...)
        HM = band_flux(s, 1800.0)
        p2, p4 = capP(HM)
        d2, _ = capP(m.F * s.mband)
        ds = dante(d1, s)
        @printf("  %-6s %5.2f   %25.3f  %+.3f  %+.3f  (%+.3f)  %+15.3f   %10.3f, %.1f eV\n", name,
                mband_recirculation(s, wall), mean_cap(HM) / mean_cap(sol.H), p2, p4, d2,
                wallP2(band_flux(s, 1800.0; field = :B)), ds.mband_fraction, ds.color_temperature)
    end
end
println("\nDewald 2020: F_M ≈ 2; wall M-band P2/P0 +0.53; ≈ +0.21 at the ablation front (×0.4), +0.28 at the inner HDC surface; Dante 12 % above 1.8 keV at 290 eV")
