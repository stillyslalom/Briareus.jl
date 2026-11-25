"""
    Briareus.Benchmarks

Published shot data and calibration fixtures: the digitized figures and tables of the
papers that the models were fitted to or are checked against, the fitted CBET models of
those shots, and the mock-up pointing of the examples. The names are exported from this
module only, `using Briareus.Benchmarks`, and each docstring cites its source.
"""
module Benchmarks

using ..Briareus
using ..Briareus: _BASECAMP, _data_table
using LinearAlgebra

export n221204_fig9, n220919_fig9, n210808_fig1, pre109_fig1, hybride_shots, n221204_fig10, n221204_pulse, zylstra2022_fig2, n221204_fig9_beams, n221204_cbet, hybride_cbet, hybride_cbet_iraum, young2020_radii, young2020_albedos, chen2020_leh, pre109_fig1a_wall, kritcher2024_capsule, kritcher2024_radius, izumi2021_bubble, izumi2021_absorption, young2020_green, young2020_kernel, young2020_pulses, young2020_dropped_quads, kemp2026_sbs, basecamp_picket_fluence, hybride_mock_pointing

# ---- from cbet.jl

"""
    n221204_fig9() -> NamedTuple

PRE 109, 025204 (2024) Fig. 9 for N221204 (HYDRA), digitized
(`data/pre109_fig9_n221204.csv`, sources in its header). The columns are `t` (ns, 0.7–9.0
in 0.1 ns steps), inner cone fraction `cf_in`/`cf_post` before and after CBET (Fig. 9b), and
per-cone powers `p23_in`, `p23_post`, …, `p50_post` (TW, Fig. 9a). The reading error is
about ±0.01 in cone fraction and ±2 TW per cone; the foot inner-cone powers (< 8 TW) are
below the resolution of the figure.
"""
n221204_fig9() = _data_table("pre109_fig9_n221204.csv")

"""
    n220919_fig9() -> NamedTuple

PRE 109, 025204 (2024) Fig. 9 for N220919 (HYDRA; the black curves, Δλ 2.5 Å), digitized
(`data/pre109_fig9_n220919.csv`, sources in its header). The columns are those of
[`n221204_fig9`](@ref) (`t`, `cf_in`, `cf_post`, `p23_in`, …, `p50_post`; 0.7–7.4 ns,
where the black curves end), each with a 1σ reading error `u_*` and a reading flag
`f_*` (`"v"` seen, `"p"` partly covered, `"h"` hidden under N221204's red curve, value
N221204's; see the header); `NaN` where not readable.
"""
n220919_fig9() = _data_table("pre109_fig9_n220919.csv")

"""
    n210808_fig1() -> NamedTuple

Kritcher et al., PRE 106, 025201 (2022) Fig. 1(d, e) for N210808 (HYDRA, Δλ 1.8 Å):
incident and post-CBET power per cone and inner cone fraction, digitized
(`data/pre106_fig1_n210808.csv`, sources in its header). The columns are those of
[`n221204_fig9`](@ref). The incident powers carry HYDRA's multipliers (Fig. 1 caption),
as Fig. 9's do.
"""
n210808_fig1() = _data_table("pre106_fig1_n210808.csv")

"""
    pre109_fig1() -> NamedTuple

PRE 109, 025204 (2024) Fig. 1(d, e), digitized (`data/pre109_fig1_n210808.csv`, sources in
its header). The columns are the inner cone fraction before and after CBET
(`cf_in_n210808`, `cf_post_n210808`, `cf_in_n221204`, `cf_post_n221204`) and the power of
the 23° cone (`p23_in_*`, `p23_post_*`, TW) against `t` (ns), each with reading errors
`u_*` and flags `f_*`. N210808's post-CBET curve in this table is that of the 2024 HYDRA
model (the one of Fig. 9 and Fig. 10); it differs from PRE 106's ([`n210808_fig1`](@ref))
by +0.04 on average (PRE 106 higher; 0.27 against 0.17 at 3.4 ns, +0.03 in the peak) while
the incident curves agree to 0.004.
"""
pre109_fig1() = _data_table("pre109_fig1_n210808.csv")

"""
    hybride_shots() -> NamedTuple

Scalar per-shot observables of N210808, N220919, N221204 and N230729
(`data/hybride_shots.csv`, sources in its header). The columns are `shot` (strings),
`dlambda` (Å at 1ω, with values 1.8, 2.5, 2.75 and 2.9), `laser_energy` (MJ), `peak_power`
(TW), hot-spot P2 `hs_p2` and its error `u_hs_p2` (µm, `NaN` = not quoted), `yield` (MJ)
and HYDRA's published clamps `clamp_foot`, `clamp_peak` (`NaN` where not published).
"""
hybride_shots() = _data_table("hybride_shots.csv")

"""
    n221204_fig10() -> NamedTuple

Capsule radiation-flux P2/P0 against time from HYDRA (post-CBET) in PRE 109, 025204 (2024)
Fig. 10, for N210808, N220919 and N221204, digitized (`data/pre109_fig10.csv`, sources in
its header). The columns are `t` (ns, 0.8–8.9 in 0.1 ns steps plus the extrema of each
curve), `n210808`, `n220919`, `n221204` (P2/P0; negative is waist-hot) and their reading
errors `u210808`, `u220919`, `u221204`. The error is ±0.002 where a curve is visible and
0.004–0.010 where another curve covers it. N221204 lies under the other two for much of
6.6–8.9 ns. Along the late spikes the error is 0.006, and 0.015 for N221204's spike, which
is hidden under N220919's except its left flank. N221204's minimum is −0.173 at 5.58 ns.
The curves are HYDRA calculations; none is measured.
"""
n221204_fig10() = _data_table("pre109_fig10.csv")

"""
    n221204_pulse() -> NamedTuple

N221204's delivered laser power and measured Dante 1 radiation temperature, Pak et al.,
PRE 109, 025203 (2024) Fig. 1, read from the vector paths of the figure
(`data/pak2024_fig1_n221204.csv`, sources in its header). The columns are `t` (ns, 0–9.6
in 0.05 ns steps), `power` (TW, total of the 192 beams; 2.051 MJ) and `dante1` (eV; `NaN`
outside 3.3–9.4 ns, the plotted range). Dante 1 peaks at 311.9 eV at 8.30 ns. Only the
total power is published in the paper; the inner/outer split is not.
"""
n221204_pulse() = _data_table("pak2024_fig1_n221204.csv")

"""
    zylstra2022_fig2() -> NamedTuple

The HYBRID-E LEH-size pair of Zylstra et al., PRE 106, 025202 (2022) Fig. 2, read from
the vector paths of the figure (`data/zylstra2022_fig2.csv`, sources in its header): the
symcaps N191118 (3.64 mm LEH) and N210601 (3.1 mm LEH). The fields are `t` (ns,
0–10.85 in 0.05 ns steps), `power_n191118`, `power_n210601` (TW, total delivered; 1.623 and
1.663 MJ, peaks 470.6 and 398.9 TW), `dante_n191118`, `dante_n210601` (eV, measured
with Dante; `NaN` before 3.50 and 1.60 ns, where the curves start) and `bang`, the
SPIDER bang times of 9.37 and 9.48 ns in the figure. Only the total power is published; the
inner/outer split is not.
"""
zylstra2022_fig2() =
    (; _data_table("zylstra2022_fig2.csv")..., bang = (n191118 = 9.37, n210601 = 9.48))

"""
    n221204_fig9_beams(beams, fig9, t) -> Vector{Beam}

`beams` (NIF cones `:inner23`, `:inner30`, `:outer44`, `:outer50`) with the powers HYDRA
fed to CBET at time `t` in PRE 109 Fig. 9 (`fig9 = n221204_fig9()`, nearest sample). The
outer-cone powers come from Fig. 9a and the inner-cone total from the input cone fraction
of Fig. 9b, because Fig. 9a's foot inner powers are below its resolution. The inner total
is split evenly between 23.5° and 30°, as Fig. 9a shows at peak, and the power of each
cone is shared equally by its beams. These are HYDRA's degraded powers (~370 TW at peak,
against ~440 TW delivered).
"""
function n221204_fig9_beams(beams::AbstractVector{<:Beam}, fig9, t)
    i = argmin(abs.(fig9.t .- t))
    o44, o50, c = fig9.p44_in[i], fig9.p50_in[i], fig9.cf_in[i]
    inner = c / (1 - c) * (o44 + o50) / 2
    P = Dict(:outer44 => o44, :outer50 => o50, :inner23 => inner, :inner30 => inner)
    n = Dict(g => count(b -> b.group === g, beams) for g in keys(P))
    return [remake(b; power = P[b.group] / n[b.group]) for b in beams]
end

# The calibrated two-zone models below share one form. Both zones have one Te(t) table and
# one clamp schedule (shot times, ns, shifted by `delay`). The outside zone has the fitted
# outflow Mach number and multiplier, and the inside zone is at Mach 0 with a separate
# multiplier.
# Their LEH planes default to ±5.62 mm, half the 11.24 mm length of the HYBRID-E
# hohlraum of the shots they were fitted to (`hybride_mock_pointing`'s target); a call
# site with a target passes its half-length.
const _HYBRIDE_LEH_Z = 5.62

function _calibrated_two_zone(fit; leh_z, offset, delay)
    Te = (fit.Te[1] .+ delay, fit.Te[2])
    clamp = (fit.clamp[1] .+ delay, fit.clamp[2])
    return TwoZoneCBET(; leh_z, offset,
                       outside = (mach = fit.mach, multiplier = fit.outside, Te, clamp,
                                  saturation = fit.saturation),
                       inside = (multiplier = fit.inside, Te, clamp, saturation = fit.saturation))
end

# The N221204 calibration contains the LEH electron temperature (eV) at shot times (ns),
# the clamp schedule of the deck and smooth saturation.
const _N221204_FIT = (Te = ([1.0, 3.0, 4.0, 4.75, 5.5, 6.5, 8.0],
                            [2790.0, 3550.0, 3720.0, 5190.0, 5500.0, 6030.0, 7440.0]),
                      clamp = ([0.0, 3.9, 4.1, 20.0], [0.002, 0.002, 0.0075, 0.0075]),
                      mach = 1.23, outside = 0.30, inside = 1.23, saturation = :smooth)

"""
    n221204_cbet(; leh_z = 5.62, offset = 0.5) -> TwoZoneCBET

Two-zone CBET calibrated to the N221204 HYDRA transfer (PRE 109, 025204 Fig. 9b, as
digitized in [`n221204_fig9`](@ref)), on [`hybride_mock_pointing`](@ref) beams with
Fig. 9's input powers ([`n221204_fig9_beams`](@ref)) and Δλ 2.75 Å (1ω, 0.917 Å at 3ω)
on the inner cones. Both zones have the clamp schedule of the deck (δn/n 0.002 until
3.9 ns, 0.0075 from 4.1 ns), CH at ne/nc 0.06 and Ti 1.5 keV (repository defaults), and
one fitted electron-temperature history Te(t) (2.8 keV at 1 ns, 3.7 at 4, 5.5 at 5.5,
7.4 at 8 ns; constant outside). The outside zone has a fitted outflow Mach number of 1.23
and a fitted multiplier of 0.30; the inside zone is at Mach 0 with a fitted multiplier of
1.23. Over 0.8–8.4 ns the cone-fraction rms of the fit is 0.016 and its worst residual is
0.041 (8.4 ns); the outer-cone post/input ratios (not fitted) match Fig. 9a to 0.03.
`leh_z` (mm) sets the LEH planes (`TwoZoneCBET`). Its default is half the 11.24 mm length
of the HYBRID-E hohlraum fitted, and a caller with a target passes its half-length, as
`examples/n221204.jl` does with `geometry.length / 2`.

*Derived* (research notes, "Two-zone CBET calibrated to PRE 109 Fig. 9"): the model is an
effective calibration to one HYDRA simulation of one shot and does not describe a
measured plasma state. Times are shot times, so the model applies to pulses timed like
N221204's.
"""
n221204_cbet(; leh_z = _HYBRIDE_LEH_Z, offset::Real = 0.5) =
    _calibrated_two_zone(_N221204_FIT; leh_z, offset, delay = 0.0)

# The joint fit contains the LEH electron temperature (eV), the clamp knots (ns, in
# N221204's timing) and HYDRA's hard-minimum saturation.
const _HYBRIDE_FIT = (Te = ([1.0, 3.0, 4.0, 4.75, 5.5, 6.5, 8.0],
                            [2054.0, 2502.0, 2799.0, 3949.0, 3357.0, 3622.0, 5500.0]),
                      clamp = ([3.88, 4.08], [0.00143, 0.000476]),
                      mach = 1.54, outside = 0.142, inside = 0.729, saturation = :min)

"""
    hybride_cbet(; leh_z = 5.62, offset = 0.5, delay = 0.0) -> TwoZoneCBET

Two-zone CBET fitted jointly to HYDRA's post-transfer cone fraction of three HYBRID-E
shots at three detunings: N210808 (1.8 Å, PRE 109 Fig. 1e with PRE 106 Fig. 1d's incident
powers, [`pre109_fig1`](@ref), [`n210808_fig1`](@ref)), N220919 (2.5 Å,
[`n220919_fig9`](@ref)) and N221204 (2.75 Å, [`n221204_fig9`](@ref)), with one plasma
history for all, under two published constraints: Te ≤ 5.5 keV (the LEH values of Michel
2009 and Barrios 2016) and the Δλ null inside NC24's 1–1.55 Å for HYBRID-E. That band is
the Δλ of round hot spots in the 2020–21 shots and lies above the transfer null; in
[`hybride_cbet_iraum`](@ref) the null is held at I-Raum's ≤ 0.5 Å instead. The model has
HYDRA's hard-minimum saturation (`saturation = :min`), with the clamp schedule switched at
the rise ([`cbet_clamp_schedule`](@ref)) and both clamp values fitted: 1.4×10⁻³ in the
foot (not binding, unconstrained) and 4.8×10⁻⁴ at peak (binding; ±27%). The outside zone
has Mach 1.54 and multiplier 0.14 (×/÷ 2), and the inside zone Mach 0 and multiplier 0.73
(±18%). The fitted Te(t) is 2.1 keV at 1 ns, 2.8 at 4, 3.4–3.9 in the peak and 5.5 (the
bound) at 8 ns (±6–14%).

The fit (107 rows, 0.8–8.4 ns) has a cone-fraction rms of 0.027 (N221204), 0.016 (N220919)
and 0.047 (N210808), against 0.016 / 0.030 / 0.068 for [`n221204_cbet`](@ref); the null is
0.97 ± 0.10 Å at 5.6 ns. The model does not reproduce HYDRA's flat cross-shot response.
Its dCF/dΔλ at peak is 0.17–0.20 per Å, whereas HYDRA's N220919 and N221204 differ by
0.03–0.06 per Å, so N210808's peak is 0.07–0.10 low and N221204's 0.03 high. `delay` (ns)
shifts the time tables for a pulse whose rise comes `delay` later than N221204's (a
negative `delay` for an earlier rise); `leh_z` is as in [`n221204_cbet`](@ref).

*Derived* (research notes, "Joint multi-shot CBET calibration (2026-10-24)"): the
parameters are effective values fitted to HYDRA and do not describe a measured plasma
state.
"""
hybride_cbet(; leh_z = _HYBRIDE_LEH_Z, offset::Real = 0.5, delay::Real = 0.0) =
    _calibrated_two_zone(_HYBRIDE_FIT; leh_z, offset, delay)

# The joint fit with the null at I-Raum's ≤ 0.5 Å
const _HYBRIDE_IRAUM_FIT = (Te = ([1.0, 3.0, 4.0, 4.75, 5.5, 6.5, 8.0],
                                  [3223.0, 3899.0, 4522.0, 3806.0, 2544.0, 3031.0, 5500.0]),
                            clamp = ([3.88, 4.08], [0.000912, 0.0001416]),
                            mach = 1.489, outside = 0.0723, inside = 1.601, saturation = :min)

"""
    hybride_cbet_iraum(; leh_z = 5.62, offset = 0.5, delay = 0.0) -> TwoZoneCBET

The joint three-shot fit of [`hybride_cbet`](@ref) (same data, model, Te ≤ 5.5 keV bound
and fitted hard-minimum clamps) with the Δλ null constraint corrected. NC24's
1–1.55 Å for HYBRID-E is the detuning that minimized hot-spot P2 in the 2020–21 shots,
which also compensates the bubble absorption of the inner beams, so it lies above the
transfer null. The published statement on the null itself is I-Raum's (0.5 Å "mostly to
negate the flow-induced CBET", NC24). Here the null at N221204's peak (5.6 and 7.0 ns) is
held in 0.25–0.5 Å. The lower edge was chosen in Briareus so that a flow-driven transfer
remains at Δλ = 0.

The outside zone has Mach 1.49 (± 0.8) and multiplier 0.072 (×/÷ 5), and the inside zone
Mach 0 and multiplier 1.60 (± 31%). The fitted clamps are 9.1×10⁻⁴ (foot, not binding,
unconstrained) and 1.4×10⁻⁴ (peak, binding, ± 29%). The peak clamp lies below the
5×10⁻⁴–7.5×10⁻³ of HYDRA decks and is correlated with the inside multiplier and Te
(correlations −0.92 to −0.97). The fitted Te(t) is 3.2 keV at 1 ns, 4.5 at 4, 2.5–3.8 in
the peak and 5.5 (the bound) at 8 ns (± 7–15%).

The fit (107 rows, 0.8–8.4 ns) has a cone-fraction rms of 0.024 (N221204), 0.016 (N220919)
and 0.028 (N210808), and 0.024 over all shots, against 0.027 / 0.016 / 0.047 and 0.033
over all shots for `hybride_cbet`. Its loss is 0.059, against 0.120 with the null in
1–1.55 Å. The null is 0.24 ± 0.20 Å at 5.6 ns (0.29 at 7 ns), and dCF/dΔλ is 0.137 ± 0.011
per Å at peak (2.5–2.75 Å), linear from 1 Å. HYDRA's cross-shot response stays 2–3×
flatter, which is not supported by the measured hot-spot P2 sensitivity (research notes).
`delay` and `leh_z` are as in [`hybride_cbet`](@ref).

*Derived* (research notes, "Joint multi-shot CBET calibration: the null"): the parameters
are effective values fitted to HYDRA and do not describe a measured plasma state.
"""
hybride_cbet_iraum(; leh_z = _HYBRIDE_LEH_Z, offset::Real = 0.5, delay::Real = 0.0) =
    _calibrated_two_zone(_HYBRIDE_IRAUM_FIT; leh_z, offset, delay)

# ---- from dynamics.jl

"""
    young2020_radii() -> NamedTuple

The radii of the LLNL dynamic view-factor model for N170827 against time, from Young et
al., PoP 27, 082702 (2020), Fig. 3a, digitized (`data/young2020_fig3a.csv`, sources in
its header). The columns are `t` (ns, 0–8 every 0.25), `wall` (the hohlraum radius, the
average of the HYDRA laser-absorption and x-ray-emission surfaces), `leh` (from SXI, of
which only the end value is measured) and `capsule` (ablation front), all in mm, and
their reading errors `u_wall`, `u_leh`, `u_capsule` (0.03 mm, 0.05 where the static line
hides the curve). The plotted initial wall radius is 3.0 mm, so displacements should be
used. The table is the calibration target of the defaults of [`WallMotion`](@ref).
"""
young2020_radii() = _data_table("young2020_fig3a.csv")

"""
    young2020_albedos() -> NamedTuple

The albedos of the LLNL dynamic view-factor model for N170827 against time, from Young
et al., PoP 27, 082702 (2020), Fig. 3a top panel, digitized
(`data/young2020_fig3a_albedo.csv`, sources in its header). The columns are `t` (ns,
0–8 every 0.25, plus 0.55–0.65 on the drop of the wall), `wall` and `capsule`
(quasi-1D HYDRA albedos of gold and the HDC ablator) and their reading errors `u_wall`,
`u_capsule` (0.03; 0.05 where the static line hides the curve or on the drop). The
static model uses 0.8 and 0.3. With [`young2020_radii`](@ref) the table defines the
dynamic model of Young, which `TabulatedAlbedo` builds.
"""
young2020_albedos() = _data_table("young2020_fig3a_albedo.csv")

"""
    chen2020_leh() -> NamedTuple

Measured effective LEH radius against time on five HDC basecamp shots (5.75 × 10.13 mm Au
hohlraum, 3.37 mm LEHs, r₀ = 1.685 mm, about 400 TW from about 4.2 to 6 ns), from
H. Chen et al., Phys. Plasmas 27, 022702 (2020), Fig. 27. The values were taken from the
vector graphics of the accepted manuscript (`data/chen2020_fig27.csv`, sources in its
header). The columns are:
- `shot`: N171010, N171212, N180102, N180103, N180117;
- `kind`: `"gleh"` for the edge from the gated LEH imager, with an analysis error of
  10–20 µm; `"lasnex_f003_regional"`, `"lasnex_f003_global"` and `"lasnex_f015_global"`
  for LASNEX images analyzed the same way; `"laser"` for the total delivered power;
- `t` (ns) and `value` (radius in mm, power in TW).

The table is the calibration data of [`LipClosure`](@ref).
"""
chen2020_leh() = _data_table("chen2020_fig27.csv")

"""
    pre109_fig1a_wall() -> NamedTuple

The inner edge of the wall material at 6 ns in the HYDRA simulations of N210808 and
N221204 drawn in PRE 109, 025204 (2024) Fig. 1a, digitized (`data/pre109_fig1a_wall.csv`,
sources in its header). Columns `z` (mm from the waist, every 0.1 mm), `n210808`,
`n221204` (edge radius, mm; initial wall 3.2 mm) and `region` (`"waist"`, `"inner"`,
`"collision"`, `"outer"`, `"bubble"`, `"lip"`). Laser rays are drawn over the plasma, so
the edge is where they end (the laser-absorption surface) or, under crossing rays, an
upper bound on the radius. The waist rows (z < 0.9 mm) give the usable number, which
puts the wall 0.38 ± 0.04 mm in at 6 ns.
"""
pre109_fig1a_wall() = _data_table("pre109_fig1a_wall.csv")

"""
    kritcher2024_capsule() -> NamedTuple

The scaled capsule radius R/R(0) from HYDRA against time on the 2.05 MJ drive of
N221204, from Kritcher et al., PoP 31, 070502 (2024) Fig. 1, digitized
(`data/kritcher2024_fig1.csv`, sources in its header). The columns are `t` (ns, 0–9.5
every 0.25), `n221204` (the as-shot, thicker capsule), `n210808` (the thinner N210808
capsule on the same drive, NaN after its minimum at 9.18 ns) and `u` (reading error of
`n221204`). The paper does not state which radius is plotted. The radius stays constant
until about 3.6 ns, near the breakout of the first shock from the ablator into the fuel
(3.86–3.91 ns in keyhole N230110, the N221204 design). That is about 2 ns before the
ice–gas surface can move. The curve therefore follows the ablator–fuel interface
(1.050 mm initially) or the center of mass of the fuel, and the ablation front is ruled
out. It reaches its minimum at bang time (9.53 ns). It is the only published radius
history of the N221204 capsule. See the research notes, "Capsule radius schedules:
published histories and the API (2026-10-26)"; [`kritcher2024_radius`](@ref) gives it
as a schedule.
"""
kritcher2024_capsule() = _data_table("kritcher2024_fig1.csv")

"""
    kritcher2024_radius(capsule = :n221204) -> TabulatedRadius

[`kritcher2024_capsule`](@ref) as a [`RadiusSchedule`](@ref). `:n221204` gives the
R/R(0) of the as-shot capsule on the drive of N221204 (0–9.5 ns). `:n210808` gives the
thinner N210808 capsule on the same drive, digitized to 9.0 ns and closed at its
minimum, 0.055 at 9.18 ns; the curve after the minimum is not digitized. These are
published HYDRA histories (Kritcher et al., PoP 31, 070502 (2024) Fig. 1, digitized).
Applied to the outer radius of the capsule, they are the upper bound on the ablation
front, for a shell of constant relative thickness. For another shot,
[`ScaledRadius`](@ref) moves the minimum to its bang time (derived).
"""
function kritcher2024_radius(capsule::Symbol = :n221204)
    k = kritcher2024_capsule()
    capsule === :n221204 &&
        return TabulatedRadius(k.t, k.n221204;
                               source = "Kritcher et al., PoP 31, 070502 (2024) Fig. 1, N221204 (HYDRA, published)")
    capsule === :n210808 || throw(ArgumentError("capsule must be :n221204 or :n210808, not :$capsule"))
    ok = findall(isfinite, k.n210808)
    return TabulatedRadius([k.t[ok]; 9.18], [k.n210808[ok]; 0.055];
                           source = "Kritcher et al., PoP 31, 070502 (2024) Fig. 1, N210808 capsule on N221204's drive (HYDRA, published)")
end

# ---- from bubble.jl

"""
    izumi2021_bubble() -> NamedTuple

Au bubble tip displacement from the initial wall against time on the HDC SymCap series of
Izumi et al., PoP 28, 022706 (2021), Fig. 8 (5.75 mm Au hohlraum, outer turn-on 0.3 ns),
extracted from the vector graphics of the figure (`data/izumi2021_fig8.csv`, sources in
its header). The columns are `shot`, `source` (`"data"` for measured values, Fig. 8a,
error bars ±40 µm; `"hydra"` for 2D HYDRA, Fig. 8b), `energy` (outer picket fluence over
that of the basecamp), `fill` (mg/cc), `capsule` (radius over that of the basecamp),
`t` (ns), `x_bub` (mm) and `u` (mm). The table is the calibration target of the `speed`
and `delay` of [`BubbleAbsorption`](@ref).
"""
izumi2021_bubble() = _data_table("izumi2021_fig8.csv")

"""
    izumi2021_absorption() -> NamedTuple

The fraction of the inner-cone power absorbed in the bubble against time for the
0.3 mg/cc basecamp of Izumi et al., PoP 28, 022706 (2021), Fig. 10a, from the vector
graphics of the figure (`data/izumi2021_fig10a.csv`). The columns are `source`, which is
`"hydra"` for the 2D simulation (0, 0.02, 0.23, 0.55 at 3, 4, 5, 6 ns) or `"model"` for
Eq. 6 of the paper (every 0.1 ns), `t` (ns), `absorbed` and `u`. The table is the
calibration target of the `opacity` of [`BubbleAbsorption`](@ref).
"""
izumi2021_absorption() = _data_table("izumi2021_fig10a.csv")

# ---- from implosion.jl

"""
    young2020_green() -> NamedTuple
    young2020_kernel() -> ResponseKernel

The mode-1 Green's function of NIF shot N170827 (km/s of final hot-spot velocity per % of
flux mode 1 per ns) from 59 2D capsule-only HYDRA impulse runs in Young et al., Phys.
Plasmas 27, 082702 (2020), Fig. 5, digitized twice (`data/young2020_fig5.csv`, sources and
reading errors in its header). `young2020_green()` has columns `t` (ns, 0.25–7.875 every
0.125), `panel_b` and `panel_d` (`NaN` where the curve is hidden under another stroke);
`young2020_kernel()` is their mean, or either one where the other is missing, as a
`ResponseKernel` for `hotspot_velocity(history; kernel)`. The kernel is negative in the
foot and dominated by the rise and peak. Its integral is ~117 km/s per %, the
constant-mode-1 sensitivity on this platform, against 100 for
`HOTSPOT_VELOCITY_PER_PERCENT`. The kernel is platform specific: that of N170601 is almost
identical, and that of N180128 puts more weight on the rise.
"""
young2020_green() = _data_table("young2020_fig5.csv")

@doc (@doc young2020_green)
function young2020_kernel()
    d = young2020_green()
    g = [isnan(b) ? c : isnan(c) ? b : (b + c) / 2 for (b, c) in zip(d.panel_b, d.panel_d)]
    keep = .!isnan.(g)
    return ResponseKernel(d.t[keep], g[keep])
end

"""
    young2020_pulses() -> NamedTuple

The requested laser power per beam of NIF shot N170827 (TW), `(; inner, outer)` as
`PulseShape`s (inner for 23.5° and 30°, outer for 44.5° and 50°), from Young et al., Phys.
Plasmas 27, 082702 (2020), Fig. 2e, digitized (`data/young2020_fig2e.csv`, 0.25–7.75 ns
every 0.125, reading error ~0.03 TW). The pulses are passed to
`with_pulses(beams, pulses, t; by)` with
`by = b -> startswith(String(b.group), "inner") ? :inner : :outer`.
"""
function young2020_pulses()
    d = _data_table("young2020_fig2e.csv")
    pulse(v) = (k = .!isnan.(v); PulseShape(d.t[k], v[k]))
    return (; inner = pulse(d.inner), outer = pulse(d.outer))
end

"""
    young2020_dropped_quads() -> NamedTuple

Flux mode 1 (%, the norm of the l = 1 vector) on the capsule of N170827 with one quad
dropped from each cone and the requested pulses, in the LLNL static view-factor model
(radii at 85%, albedos 0.8/0.3), from Young et al., Phys. Plasmas 27, 082702 (2020),
Fig. 10a, digitized (`data/young2020_fig10a.csv`, sources in its header). The columns are
`t` (ns, 0.25–7.75 every 0.125), `inner23`, `inner30`, `outer44`, `outer50` (`NaN` where
hidden). Fig. 10c's hot-spot velocities (km/s, read off the markers ±2) are in `velocity`
(requested pulses) and `velocity_delivered`, in the same cone order.
"""
function young2020_dropped_quads()
    d = _data_table("young2020_fig10a.csv")
    return (; d..., cones = (:inner23, :inner30, :outer44, :outer50),
            velocity = (316.9, 336.9, 169.7, 145.9), velocity_delivered = (347.5, 358.0, 203.2, 181.2))
end

# ---- from lpi.jl

"""
    kemp2026_sbs() -> NamedTuple

Measured inner-cone SBS backscatter of N210808 per quad against time, digitized from
Kemp et al., PoP 33, 082705 (2026) Fig. 2 into `data/kemp2026_fig2.csv`; the sources
are listed in its header. The columns are `t` (ns, 4.4–7.1 every 0.05), `inner23` and
`inner30`. Each column contains all backscattered SBS power of one 23° or 30° quad in TW,
cone-averaged over both hemispheres; it is the in-lens signal × 3.4. The energies in
the paper are 2.8 and 0.9 kJ per quad, with eight quads per cone. Single-shot errors
reach 100%, and N210808 repeats vary by 30–50%.

SBS starts as the inner power reaches its peak. On this axis the pre-CBET quad power is
at half its 9 TW peak at 4.74 ns. SBS peaks at 5.3 ns (23°) and 5.75 ns (30°) and ends
by 6.6–6.9 ns. The 30° cone ends first, which the paper attributes to refraction and
absorption by an on-axis density spike (§III C). `examples/n221204.jl` uses this table
for `backscatter=history`. Research notes, "Inner-cone SBS time histories and the
late-time refraction (2026-10-25)".
"""
kemp2026_sbs() = _data_table("kemp2026_fig2.csv")

# ---- from picket_fluence.jl

# The HDC SymCap basecamp of Izumi 2021, from the header of `data/izumi2021_fig8.csv`.
# The fields are:
# - the length and LEH radius (mm), with the wall and capsule radii of `_BASECAMP`;
# - the outer picket energy (kJ, all outer beams, first ns after the outer turn-on), read
#   by eye from Izumi 2021 Fig. 3a (research notes, "The P2 swing: capsule radius,
#   backscatter and the collision");
# - the outer aims (mm on the axis, upper hemisphere), the proxy pointing of the package
#   for a 5.75 mm hohlraum.
const _BASECAMP_PICKET = (length = 10.13, leh_radius = 1.685, energy = 26.0,
                          aim = (outer44 = 5.20, outer50 = 5.01))

"""
    basecamp_picket_fluence(; energy = 26.0, h = 0.2, capsule_h = 0.4) -> PicketFluence

The outer picket fluence of the calibration platform of [`BubbleAbsorption`](@ref),
where the `picket_ratio` is E = 1. The platform is the HDC SymCap basecamp of Izumi et
al., PoP 28, 022706 (2021), a 5.75 × 10.13 mm Au hohlraum with 3.37 mm entrance holes and
a 0.908 mm capsule. The mesh size is `h` (mm), and `capsule_h` is that of the capsule.
The 128 outer beams ([`nif_beams`](@ref), with the NIF outer phase-plate spots) carry
equal shares of `energy`, the outer picket energy in kJ over all outer beams. The value
is about 26 kJ in the first ns after the outer turn-on, read by eye from Izumi 2021
Fig. 3a (*L*). The quads are unsplit and aimed on the axis at z = ±5.20 mm (44.5°) and
±5.01 mm (50°). This is the proxy pointing of the package for this hohlraum, the one
with which the bubble opacity was calibrated; the basecamp pointing of the facility is
unpublished (*L*). The result is the reference of [`picket_ratio`](@ref) for an absolute
E, and the map of the design should be computed on the same `h`.
"""
function basecamp_picket_fluence(; energy::Real = _BASECAMP_PICKET.energy, h::Real = 0.2,
                                 capsule_h::Real = 0.4)
    t = hohlraum(; length = _BASECAMP_PICKET.length, radius = _BASECAMP.radius,
                 leh_radius = _BASECAMP_PICKET.leh_radius, capsule_radius = _BASECAMP.capsule,
                 h, capsule_h)
    outer = [b for b in nif_beams() if b.group in (:outer44, :outer50)]
    beams = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * getfield(_BASECAMP_PICKET.aim, b.group)),
                    power = energy / length(outer)) for b in outer]
    return picket_fluence(DepositionOperator(t, beams), beams)
end

# ---- from facilities.jl

# Mock HYBRID-E pointing per cone, with aim z (mm, upper hemisphere), aim radius r (mm,
# along the port azimuth of the quad, negative across the axis) and tangential split s
# (mm). It is an illustrative stand-in and does not reproduce a facility pointing.
const HYBRIDE_MOCK_POINTING = Dict(
    :inner23 => (z = 3.3, r = -1.3, s = 0.0),
    :inner30 => (z = 4.2, r = -0.8, s = 0.0),
    :outer44 => (z = 5.7, r = 0.0, s = 0.5),
    :outer50 => (z = 5.8, r = 0.0, s = 0.5),
)

"""
    hybride_mock_pointing(beams = nif_beams(); table = HYBRIDE_MOCK_POINTING) -> Vector{Beam}

A mock of the HYBRID-E (N221204/N230729-class) NIF pointing. `beams` (NIF beams or
quads, with `group` the cone) are re-aimed by cone with defocus 0, and power, spot and
`dlambda` are kept. The quads of each cone share one rule in the frame of each quad,
mirrored top to bottom. The aim lies at `±z` (mm, the bottom hemisphere at `-z`), at
`r` along the port azimuth of the quad (negative across the axis, toward the opposite
wall) and, for the outer cones, at `±s` along φ̂ of the quad. Each beam aims to the side
opposite its port (`-sign(φ_beam − φ_quad) s`), so the beams of the quad cross. The quad
azimuth is that of the mean port of its beams in `beams`, so a quad given as one beam
(`nif_quads()`) is not split.

| cone | z | r | s |
|---|---|---|---|
| 23.5° | 3.3 | −1.3 | 0 |
| 30° | 4.2 | −0.8 | 0 |
| 44.5° | 5.7 | 0 | 0.5 |
| 50° | 5.8 | 0 | 0.5 |

*Derived*; the rule is an illustrative HYBRID-E-like rule on a 0.1 mm grid, for the
11.24 × 6.4 mm hohlraum with 3.1 mm LEHs, and does not reproduce a facility pointing.
"""
function hybride_mock_pointing(beams::AbstractVector{<:Beam} = nif_beams();
                               table = HYBRIDE_MOCK_POINTING)
    sums = Dict{String,Vec3}()
    for b in beams
        sums[b.quad] = get(sums, b.quad, zero(Vec3)) + Vec3(b.port...)
    end
    return map(beams) do b
        p = table[b.group]
        q = sums[b.quad]
        φq = atan(q[2], q[1])
        δ = mod(atan(b.port[2], b.port[1]) - φq + π, 2π) - π
        side = abs(δ) < 1e-9 ? 0.0 : sign(δ)
        er, eφ = Vec3(cos(φq), sin(φq), 0), Vec3(-sin(φq), cos(φq), 0)
        aim = p.r * er - side * p.s * eφ + Vec3(0, 0, sign(b.port[3]) * p.z)
        remake(b; aim, defocus = zero(b.defocus))
    end
end

end # module Benchmarks
