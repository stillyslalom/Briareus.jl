# Research notes: Hohlraum geometry in time: wall motion, bubbles, LEH closure and windows

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

## HYBRID-E (N210808) hohlraum and LEH size (2026-10-20)

Source: Kritcher et al., Phys. Rev. E 106, 025201 (2022), Table I, Secs. II A and II C
(full text, MIT PSFC copy). **H**
- Hohlraum: Au-lined depleted uranium, 6.4 mm diameter, 11.24 mm long, 0.3 mg/cm³ He.
- LEH diameter: 3.64 mm on N210207 and N210307, 3.1 mm on N210808. Callahan's HEDS
  slides (2022) quote 3.65 → 3.1 mm, 27% less area. **H**
- Capsule: HDC, 1.05 mm inner radius, ablator 76–80 µm (density 3.5 MCD / 3.32 NCD);
  the tutorial's 1.13 mm outer radius is *derived*.
- Peak laser power / energy: N210307 487 TW / 1.909 MJ; N210808 441 TW / 1.917 MJ.
  The smaller LEH "was used to operate at lower peak power ... while still achieving
  similar radiation temperature" (qualitative). Δλ 1.55 → 1.8 Å (1ω); the 50° beams
  moved toward the midplane to clear the smaller LEH; peak cone fraction 33%.

`docs/literate/leh_power_balance.jl` and `test/leh_power_balance.jl` (*derived*):
- Inputs that are engineering estimates, not shot data: wall albedo 0.925
  (`WallLossAlbedo` at 300 eV, 6 ns; the Au-lined DU wall is not modeled), capsule
  albedo 0.3, the x-ray source as uniform bands at 0.6 ≤ |z| ≤ 1.8 mm (inner, 33%)
  and 3.4 ≤ |z| ≤ 4.6 mm (outer). Moving the outer bands ±0.5 mm moves the fractions
  by < 0.005.
- Holes are open and static; LEH closure and the gold bubble shrink the effective hole
  during the pulse, so the computed LEH loss is the open-hole bound.
- At 3.1 mm: 34.4% of the x-ray source escapes through the LEHs, 23.7% is absorbed by
  the capsule; per 0.1 mm of diameter over 3.0–3.6 mm, +1.6% escape, −0.5% capsule.
- Briareus (`closed = true`, h = 0.4) against a Monte Carlo photon walk in the exact
  geometry: within 1e-4 (10⁸ histories, σ 3.5e-5). The uniform-flux power balance
  gives the LEH loss 6.5–7% low and the capsule's share ~1% high.
- Required-power ratio 3.1 vs 3.64 mm at fixed capsule drive: 0.890 (α_w 0.925),
  0.905 (α_w 0.90), against 441/487 = 0.906. This is consistent in size only, since the
  shots also changed Δλ, pointing and pulse length. **L**

## Literature review against Kritcher et al. 2024 (2026-10-22)

Kritcher et al., Phys. Plasmas 31, 070502 (2024), doi:10.1063/5.0210904 (the N221204
design review; CC-BY, full text **H**) was read for the modeling gaps it exposes in a
view-factor code, and six focused reviews were run against it. The full reports, with
every number's source and confidence, are in `reference/litreview/`:
`kritcher2024_notes.md` (the paper), `hohlraum_dynamics.md`, `cbet_calibration.md`,
`energetics.md`, `asymmetry_transfer.md`, `n221204_benchmark.md`, `hydro_coupling.md`.
What follows is the digest that sets roadmap items 9–13 and 16–17.

### The modeling chain of the paper **H**
- Symmetry is held to ~1% time-varying flux asymmetry; the wall "ingresses into the path
  of laser beams"; late-time symmetry comes from a semi-empirical model (Callahan 2018)
  calibrated on years of data, not first-principles LPI.
- The P2 swing (time-varying P2) is the largest unknown, calibrated by in-flight
  radiography at ~200 µm; for N221204 that calibration was missing, and the N230729
  post-shot was matched by reducing the swing by hand.
- HYDRA settings: inline CBET with δn/n clamp 0.002 until the rise, 0.0075 after (PRE 109
  rounds to 0.008 in the PoP); electron flux limiter 0.15 (conduction, not radiation);
  laser multipliers fit to in-flight symmetry, velocity and keyhole shock velocities;
  the laser is degraded 9–13% to match data (PRE 109, 025204).
- Mode-1 sources named: laser delivery imbalance, capsule and ice thickness, non-axisymmetric
  diagnostic windows, "manufacturing gaps on the waist" (Hurricane PRL 132, 065103).
- An applied drive mode 1 of ~0.7% reproduces >100 km/s burn-on (≈75 km/s burn-off).

### Hohlraum dynamics (`hohlraum_dynamics.md`)
Two surfaces move at different speeds, and a view-factor model should treat them apart:
- **Au/He bubble tip** (absorbs the inner beams): ballistic at 300–460 µm/ns. 182 ± 11
  µm/ns on a 3× picket CH shot (Ralph PoP 25, 082701, 2018, **H**); 350/400 µm/ns at
  0.3/0.15 mg/cc (Izumi RSI 89, 10K111, 2018, **H**); 321–461 km/s measured against
  LASNEX f = 0.03/0.15 (Meezan PoP 27, 2020, **H**); x_bub ∝ E_picket^0.2 ρ₀^−0.2
  (Izumi PoP 28, 022706, 2021, **H**). Izumi 2021 gives an "opaque iris" clipping
  fraction of the inner cones and an absorbed fraction with bubble transparency 29%
  (0.3 mg/cc) and 54% (0.15) at 5.5 ns; Ralph puts inner-beam absorption in the bubble
  at 51–62% and finds a view factor at the measured peak wall position needs cone
  fraction 18.5% for zero flux P2.
- **X-ray-emitting wall**: 50–100 µm/ns; 400–600 µm inward by peak, and VisRad models of
  LEH images needed the wall radius reduced by 300–700 µm (Chen PoP 27, 022702, 2020, **H**).
- **LEH**: closes in the low-power foot, reopens at the rise, varies ≤ 4% through peak,
  about 0.84 r₀ at peak (Schneider PoP 22, 122705, 2015, **H**); closure runs ~100 µm in
  diameter faster than simulated (Swadling "Build-A-Hohlraum", PoP 32, 052707, 2025, **H**).
- **Callahan's metric** (PoP 25, 056305, 2018, **H**): P2 = 0 at cone fraction 33% when
  √(E_pic,out/(A_out ρ_fill))·(τ/R)(r_cap/R) ≈ 1.2 (kJ, mm², mg/cc, ns, mm); 39 shots,
  σ = 8 µm; 2.7 µm hot-spot P2 per 1% cone fraction. Izumi restates it as
  P2 − 270(CF − 0.333) = 85(1.2 − FoM) µm.
- **Prior art**: Young et al., PoP 27, 082702 (2020, **H**) is LLNL's "dynamic view factor":
  VisRad at 100 times with as-delivered powers, the cylinder shrinking to the mean of the
  absorption and emission surfaces from HYDRA, LEH closure from SXI; a static model with
  all radii at 85%, albedo 0.8/0.3 and conversion 0.9 agrees to 0.14% in mode 1 at peak
  and < 5% in hot-spot velocity. The research-notes entry crediting IRAD3D with spot
  motion is not supported by its abstracts (**L**).
- A simulated flux P2/P0(t) for N210808/N220919/N221204 exists only as PRE 109, 025204
  Fig. 10 (**L**, by eye): ≈ 0 in the foot, −0.17 at 5.5 ns, 0 at 6.5 ns, +0.12–0.15 at
  8 ns. No tabulated swing is public.

### CBET calibration (`cbet_calibration.md`)
- **Clamps are per-phase effective gains, not limits** (**L**): linear δn/n at NIF
  conditions is 1e-4 to 1e-3, below every ignition-era clamp. Values: Kritcher PRE 98,
  053206 (2018, 0.6 mg/cc CH, **H**): 0.01 in the foot; 5e-4 to 2e-3 in rise/peak gives
  1.36–1.6× inner power, 2–3e-3 gives 1.67–1.9× (all HYDRA). Pickworth PoP 27, 102702
  (2020, **H**): 0.01 throughout plus peak cone-fraction multipliers 0.86 (0 Å) and 0.80
  (1 Å). N221204: 0.002 foot, 0.0075 from the rise (PRE 109, 025204, **H**). Strozzi PRL
  118, 025002 (2017) caps δn with a hard minimum; Briareus's smooth clamp is √2 lower at
  the clamp.
- **N221204 transfer** (PRE 109 Figs. 9–10, by eye, **L**): cone fraction input → post-CBET
  0.13–0.15 → 0.18 (foot), 0.34 → 0.64 (5.5 ns), 0.34 → 0.52 (8 ns); per cone at early
  peak 23° ×2.0, 30° ×1.9, 44° ×0.45, 50° ×0.76; late peak 44° ×0.62, 50° ×0.84. Symmetry
  by power alone would need the outers cut to 150–200 TW (**H**).
- **Hot-spot P2 per Å** (**H**): 40–50 µm (HYBRID-E, 0.3 mg/cc), 24 ± 3 µm (CH, 0.6 mg/cc,
  Pickworth), 14 µm (CH, 0.45 mg/cc, Izumi, at 30–40% lower peak power: transfer scales
  with intensity). Converted at 2.7 µm per % cone fraction: dCF_eff/dΔλ ≈ 0.15–0.18 per Å
  for HYBRID-E, 0.085 (0.6 mg/cc CH), 0.05 (0.45 mg/cc) (**L**).
- **Where it happens**: outside the LEH the supersonic outflow drives transfer from inner to
  outer (into the 50° cone), and inside it the transfer runs from outer to inner, so Δλ = 0
  is not zero transfer; I-Raum needed 0.5 Å to null the flow-driven part (Pickworth HYDRA;
  Ralph Nat. Commun. 2024, **H**). The 44° cone gives up more than the 50°.
- **Plasma state**: no measured low-fill LEH ne, Te, Ti or flow history is public; Michel
  PoP 20, 056308 (2013, high fill, **H**) gives ne 0.03 nc, Te 2.8 keV, Ti 0.8 keV at the
  start of the main pulse with ion heating at 6–7 keV/ns, after which the gain is linear in
  Δλ; Barrios PoP 23, 056307 (2016) measured peak Te 4.2 keV, 0.5 keV above HYDRA.
- **Backscatter** (PRE 106, 025201, **H**): N210808 total 13.4 kJ of 1.92 MJ (0.7%); 50°
  SBS 4.8 ± 2.4 kJ (9.5 ± 3.4 on N210207); the only damage number found is ~1 J/cm² SBS
  fluence (Pickworth); the criterion that limited N220919 is not public.

### Energetics (`energetics.md`)
- **N221204 energy split** (HYDRA calibrated to data; PRE 109, 025204 App. 3, **H**):
  wall 54.1% (Au lining 4.4, DU 49.7), LEH 25.2%, He 4.5%, capsule 16.2% of incident
  laser energy; 246 kJ "coupled to the implosion". The 245–258 kJ "capsule absorbed
  energy" quoted elsewhere (PoP 31 Table III; Hurricane PRL 132 Table I, Y/G_cap) is the
  energy left in the unablated shell, not total x-ray absorption. A 0-D estimate with
  capsule albedo 0.2–0.3 on the outer area and Dante's T(t) gives 316–361 kJ, bracketing
  the 16.2% (≈ 290–330 kJ): calibrate a capsule model to the 16.2%, not to 246 kJ (**L**).
- **Wall material**: U-rich and cocktail walls lose 0.83–0.88 × gold at 270–310 eV, ratio
  ∝ T^−0.22 t^−0.05 (Rosen 2005; Schein PRL 98, 175003, 2007, **H**); DU is worth 6.5% of
  laser power at peak (Döppner PRL 115, 055001, 2015, **H**, 0.7 µm Au liner). Implied
  multiplier on gold's 1 − α at HYBRID-E peak ≈ 0.80–0.85 (**L**).
- **Conversion efficiency**: 0.7 t^0.12–0.2 (Lindl 2004; Rosen 2005), 0.85 t^0.2 (Schein),
  plateau 95% at 10¹⁴ and 80% at 10¹⁵ W/cm² (Dewald 2008, **M**), 85–90% in NIF vacuum
  hohlraums (Olson 2012, **M**). The Briareus default 0.75 is a Nova 1 ns number.
- **The only published VisRad energetics calibration**: η = 0.95, wall albedo 0.77 match
  Dante 1/2 peak intensity on a 5.75 mm Au hohlraum at 400 TW (Swadling 2025, **H**); its
  0-D model uses capsule albedo ≈ 0.2 and gives α_wall ≈ 0.84 at 288 eV, with
  (1/P_cap) dP_cap/dα_w ≈ 4.5.
- **Trad**: drive Trad 307 eV (N210808 average) and 313 eV (N220919, N221204)
  (Abu-Shawareb PRL 132, 065102 Table I, **H**); Dante 1 313 ± 2 eV at θ 143°, φ 274°
  (37° through the LEH), Dante 2 308 ± 3 eV (Rubery PRL 132, 065104, **H**); the simulated
  internal Trad is ~13 eV below Dante 1 at peak and 15–20 eV below in the trough (PRE 109
  Fig. 1c, **L**). LEH pair at fixed Trad: 3.64 mm at 465 TW → 295 eV, 3.1 mm at 390 TW →
  300 eV (Zylstra PRE 106, 025202 Fig. 2, **H-img**).

### Flux asymmetry → implosion observables (`asymmetry_transfer.md`)
- **Mode 1 and the hot-spot velocity**: ~100 km/s per 1% of a₁/a₀ (burn-off), linear to about
  0.75% then saturating (Young 2020: 0.44% → 37.6 km/s 3D HYDRA vs 44.7 km/s model;
  Kritcher 2024: 0.7% → 75 km/s burn-off, **H**; ratio **L**). Young's recipe is a
  Green's function G(t) in km/s per % per ns from 59 capsule-only impulse runs: negative in
  the foot, dominated by rise and peak. Source sensitivities (N170827): 5.4 km/s per 10%
  wall–window albedo mismatch; 0.1 / 1.1–1.3 km/s per µm horizontal / vertical offset;
  1.6 km/s per degree of tilt; a dropped inner quad ≈ 2× an outer one.
- **Capsule thickness** mode 1: 110 (sim) to 140 ± 30 (data) km/s per %, opposite in
  sign to a flux mode 1 (+1% thickness at a pole ≈ −1% flux there) (Casey PRL 126,
  025002, 2021, **H**).
- **Piston model** (Hurricane PoP 27, 062704, 2020; PoP 29, 012703, 2022, **H**):
  f = v_HS/v_imp = √(RKE/(c·KE)), c ≈ 0.66; no-α yield ∝ (1 − f²)^3.3; α factor
  exp(−(a − 1)χ^1.2 f²) (Casey PRE 108, L053203, 2023 has an extra 1.2 in the exponent);
  at v_imp 380 km/s, 30/40/50/60 km/s cost 5/10/15/20% of yield. For general asymmetry
  RKE/(c·KE) = 1 − ρR_WHM/ρR_ave.
- **Burn-on amplification**: Nora PoP 30, 092705 (2023, **H**) fits
  0.324 ln(Y_amp + 2.239) + 0.653; Kritcher 2024 Table II implies 1.27–1.39 (**L**).
- **P2/P4**: Masse PoP 26 (2019, **H**) Green's functions: hot-spot P2 response flips sign
  for drive before ~4.5 ns (an early waist-hot drive makes it oblate), ρR response keeps one
  sign; sensitivity peaks at the rise. Milovich PoP 23, 032701 (2016, **H**): picket P2 <
  7.5%; +6% P2 per 200 µm outward outer pointing; 2% P2 per 10% albedo change. Scott PRL
  110, 075001 (2013, **H**): 10% flux P4 → ~20 µm a₄, up to 15× yield loss, sign set by
  timing. Ralph Nat. Commun. 15, 2975 (2024, **H**): yield vs hot-spot P2/P0 is Gaussian, W
  = 0.306, optimum +11.8% (HYBRID-E).

### N221204 benchmark (`n221204_benchmark.md`)
- Input deck (**H** unless marked): hohlraum 11.24 mm × 3.2 mm radius, LEH radius 1.55 mm,
  0.3 mg/cc He; capsule outer radius 1.135 mm (N210808: 1.1285), inner 1.050 mm, ablator
  85.4 µm (7.7 / 15.1 doped / 62.6; Pak PRE 109, 025203 Table I); W 0.54–0.6 at.%; DT ice
  64.5 µm; Δλ 2.75 Å at 1ω (0.917 Å at 3ω); 2.05 MJ, 440 TW; SBS 42 ± 11 kJ (inner, early
  peak), SRS 6 ± 1 kJ, 98% coupled. Per-cone delivered powers are digitizable from Di
  Nicola et al., HEDP (2024), OSTI 2440791, Fig. 1 (outer 296 TW and inner 146 TW at
  peak, cone fraction 1/3; tabulated by eye in the report, **L**). PRE 109 Fig. 9 gives
  the design cone fraction vs time, input and post-CBET (~0.65 at 5.5 ns, **L**).
- Validation targets (**H**): bang time 9.53 ± 0.07 ns (N210808 9.28 ± 0.07); drive Trad
  313 eV; capsule 16.2% / 258 kJ; hot-spot P2 0.3 ± 0.1 µm, velocity 36 ± 6 km/s at
  (151°, 54°); keyhole N230110 shock-1 breakout 3.864 ns (equator) / 3.906 ns (pole);
  2DConA N240304 (2.9 Å) shell velocity 355 ± 5 km/s, shell P2 21.7 ± 1.3 µm, stagnation
  9.876 ns; rise cone fraction 0.21 → 0.24 raises shell P2 9.1 → 27.5 µm, 0.21 → 0.16
  leaves it (Ruof et al., PoP 2026, OSTI 3367476).
- These remain unpublished: absolute pointing, defocus and CPPs per cone; any 23.5°/30° or
  44.5°/50° wavelength split; window sizes and azimuths; corner/lip geometry; Au liner
  thickness; tabulated powers (data "available on request"). The benchmark uses
  an illustrative pointing, `hybride_mock_pointing()` ("Mock HYBRID-E pointing").

### Briareus → CompactLES (`hydro_coupling.md`)
- LLNL's capsule-only step applies "a prescribed inward x-ray flux" on a spherical
  boundary, tuned to shock timing, in-flight velocity and shape, the M-band fraction and
  the DT bang time (Clark PoP 26, 050601, 2019, **H**), with asymmetries as
  F = Σ a_l(t) P_l at 100 µm from the ablation front (Masse 2019, **H**). Briareus's
  T_r(t) + a_l(t) export already has that shape; bang time, an input there, becomes the
  first prediction here.
- Stage-0 closures: P_abl ≈ 3 T_heV^3.5 Mbar, ṁ ≈ 3×10⁵ T³ g cm⁻² s⁻¹ (Lindl 2004 eqs.
  6-2/6-3, **H**, exponents partly reconstructed); carbon P_abl ∝ T^3.3 t^−0.09 to 360 eV
  (Trickey HEDP 44, 100995, 2022, **H**, no prefactor); NIC point design P_abl 140 Mbar,
  HYBRID-E hydro efficiency ~14%, imploding KE 19–33 kJ (Lindl PoP 33, 050503, 2026, **H**).
- The one fully specified public imposed-drive capsule test: Shestakov et al.,
  UCRL-JC-131429 (1998), Be/D capsule, gray diffusion, Marshak BC E − (2/c)F·n = aT_r⁴ at
  160 eV (**H**). M-band is ~15% of the spectrum above 1.8 keV at peak and pole-hot with
  a distinct P2 (+0.53 at the wall) (Dewald PoP 27, 2020, **H**).
- CompactLES (2026-10-22) has no radiation, laser, plasma EOS or ablation model; its
  SESAME/IONMIX readers exist but no solver uses them.

## Moving wall on keyframes (2026-10-23) *derived*

`WallMotion` / `moving_hohlraum` / `Keyframes` (ROADMAP "The P2 swing", first piece). The radiating
wall cylinder shrinks at `speed` from `start`; caps keep their planes, LEHs their radius,
the capsule stays; mesh counts are frozen at R₀ so every time has the same elements.
Between keyframes the exchange matrix Aᵢ Fᵢⱼ is interpolated linearly in time and divided
by the areas at that time (reciprocity kept; rebalanced for closed targets); deposition
operators are built at each time.

**Default speed 0.10 mm/ns** (superseded 2026-10-24 by a saturating wall fitted to Young 2020: "Hohlraum dynamics calibrated to published data" below; engineering default, **M**): Chen 2020 (500 µm in by 5.5 ns
at 0.3 mg/cc, 400 at 0.45: ~0.07–0.09 mm/ns); PRE 109, 025204 Fig. 1a (HYBRID-E,
measured off the figure ±0.15 mm: wall-material edge at the waist 2.80–2.85 mm at 6 ns from
3.2 mm, ~0.1 mm/ns; under the outer spots 2.2–2.5 mm).

**Cost** (NIF hohlraum 10.1 × 5.75 mm, WallLossAlbedo, 48 quads, 50-step history, 16
threads): coarse (0.6/0.3, 1,152 elements) view factors 0.3 s, deposition operator 0.6 s,
history 4.3 s; standard (0.28/0.16, 4,796) 2.1 / 0.7 / 63 s; fine (0.17/0.09, 12,512)
12.4 / 1.2 / 945 s. The history dominates (fresh LUs in the albedo iteration, ~19 s a step
at fine), so keyframes on the solve mesh are cheap: six fine keyframes plus per-step
deposition operators add ~14% to a fine history; a rebuild at every step ~70%.

**Accuracy** (standard mesh, 0.1 mm/ns from 0.3 ns, 39 steps 0.25–9.75 ns; max error
against a rebuild at every step):

| keyframes (spacing) | interpolated deposition: P2/P0, P4/P0 (pp) | per-time deposition: P2/P0, P4/P0 (pp), Tr (eV) |
|---|---|---|
| 3 (4.9 ns) | 1.56, 0.40 | 0.086, 0.110, 0.15 |
| 6 (1.9 ns) | 0.24, 0.060 | 0.035, 0.022, 0.03 |
| 11 (0.97 ns) | 0.068, 0.020 | 0.012, 0.007, 0.01 |
| 21 (0.49 ns) | 0.022, 0.008 | 0.009, 0.003, 0.005 |

Interpolating depositions blurs the outer footprints, which slide toward the LEH by
ΔR/tanθ; per-time operators (adopted) remove most of the error. A still wall reproduces
the fixed-target history to 1e-12.

**Effect** (tutorial pointing of `docs/literate/retune.jl`, not tuned): the moving wall
makes capsule P2/P0 more positive through the pulse, −0.091 → +0.044 at 9 ns against the
static wall's −0.091 (shift +0.13; 0.005 at 1 ns, 0.04 at 5 ns, 0.08 at 7 ns), and Tr
+1 to +7 eV (smaller cavity). For scale, PRE 109 Fig. 10's HYDRA swing (wall, bubble and
CBET together) is −0.17 at 5.5 ns to +0.08–0.10 at 8 ns.

## Closing entrance holes (2026-10-23)

- **Model** (`hohlraum(; leh_closure, leh_open_radius)`, `aperture(...; closure)`): each
  entrance hole is split into the hole that stays open (radius `leh_open_radius`) and a
  lip, the annulus out to `leh_radius`. The lip is an aperture with a `GapClosure`
  albedo: transparent to laser light, lines of sight and view-factor rays at all times,
  it lets all of its incident x-ray flux escape until `closure`, then across `width`
  re-emits α = s·α_after (s the smooth step) and lets 1 − s escape. `energy_balance`
  books τP as escaped and (1 − α − τ)P as absorbed; `drive_history` resolves the
  closure with `transition` steps and counts the lip with the entrance holes (not in
  `window_loss`). `ApertureKeepOut` keeps the full `leh_radius` rim (beams must still
  clear the hardware). It does not convert laser light. *Engineering model.* The plasma that
  closes the hole is not computed; closure time and radius are inputs. The Monte Carlo
  reference (`montecarlo_incident`) re-emits from a lip with its albedo at the time
  given (`albedo = albedos(t, time)`); test/montecarlo.jl checks it against `solve`.
- **Default `leh_open_radius` = 0.84 r₀.** Schneider et al., PoP 22, 122705 (2015): the
  LEH closes in the low-power part of the pulse and reopens at the rise. Its radius at
  peak is ~0.84 r₀ and varies ≤ 4% through peak. The measured open area is
  1.25–1.4× the modeled (HYDRA) area. **H** (via the 2026-10-22 lit review).
  Young et al., PoP 27, 082702 (2020): a static view-factor model with all radii at 85%
  reproduces the dynamic model at peak. **H**. 0.84 lies between them. **M**
- **Timing.** Chen et al., PoP 27, 022702 (2020): at 0.15 mg/cc a ~1550 µm-radius LEH stays
  open until ~4.5 ns, then its x-ray-emitting size drops rapidly over 4.5–5 ns. Where the
  laser passes, the beams keep it open. **M** (secondary transcription). Swadling et al.,
  PoP 32, 052707 (2025): closure runs faster than simulated (~100 µm in diameter). **H**
  Briareus has no default `closure`. The example below closes at 5.0 ns over 0.5 ns (the
  rise from the foot) with the wall's `WallLossAlbedo`. **L**: an engineering default
  to be set per design from SXI data.
- **What it changes** (*derived*, 2026-10-23). Hohlraum and pointing of
  `docs/literate/retune.jl` (10.1 × 5.75 mm, 3.1 mm LEH, `WallLossAlbedo` wall,
  h = 0.45, capsule_h = 0.15, its reference pulses, steps of 0.4 ns,
  `closed = true`). Compare open holes with `leh_closure = GapClosure(closure = 5.0,
  width = 0.5, after = WallLossAlbedo())` at 0.84 r₀, averaged over 7.6–9.4 ns:
  - capsule Tr goes from 310.2 to 320.4 eV (+10.2 eV, +3.3%);
  - P2/P0 goes from −9.19% to −6.67% (+2.5 points; closing the holes adds pole drive);
  - P4/P0 goes from −0.09% to −0.15%.
  Before 4.8 ns the two agree to < 0.01 eV and 0.003 points. A 16% smaller LEH radius
  is therefore worth ~2.5 points of P2, which is comparable to the cone-fraction
  retunes of the tutorials. Set the closure from data before trusting absolute P2.

## Inner-beam bubble absorption (2026-10-23)

`BubbleAbsorption` (`src/bubble.jl`, ROADMAP "The P2 swing", third piece; `test/bubble.jl`). A
laser–plasma model (`apply_lpi`), so it enters `drive_history(...; lpi)`,
`Scenario(...; lpi)` and composes with CBET (`lpi = (cbet, bubble)`). Sources are in
`reference/litreview/hohlraum_dynamics.md` §1 and "What Briareus should implement" §3.

- **Model.** Tip x_bub(t) = v_b·max(0, t − t_L), t_L = outer turn-on + 0.2 ns,
  v_b = v_ref (E/E_ref)^0.2 (ρ₀/0.3)^−0.2. In the plane z_b where the `footprint`
  (default 50°) beams' axes reach the wall (the moved wall when given a `WallMotion`),
  the bubble is the annulus r > R₀ − x_bub. Each inner beam loses f·(1 − T_b), with f its
  power crossing z_b inside the annulus, from the beam itself: best-focus spot (core and
  halo, 12 equal-power radii × 24 angles) carried along 13 f-cone directions, as the
  deposition model does. T_b = exp(−τ_b (ρ₀/0.3)^q (x/1 mm)^0.5 (R₀/(R₀ − x))²).
  *Engineering model.*
- **Departures from the lit-review spec, and why.**
  - A density factor (ρ₀/0.3)^q in T_b. Without it the 0.15 mg/cc bubble, being deeper,
    is *more* opaque, so no τ_b can meet both of Izumi's anchors. Two anchors fix τ_b and
    q: **τ_b = 0.1441, q = 1.862** (*derived*, `_bubble_anchor_fit`: R = 2.88 mm, launch
    0.5 ns, x = 1.75 / 2.01 mm at 5.5 ns). q ≈ 2 is the ρ² of inverse bremsstrahlung,
    which supports the form (**L**).
  - The absorbed power is not deposited on a separate wall band: it is moved to the
    footprint beams of the inner beam's hemisphere, in equal shares, so the existing
    deposition operator puts it where those beams land. Beam power is conserved
    exactly (per hemisphere), the band follows the outer pointing and the moving wall,
    and no new operator or keyframe data is needed. Equal shares (not power-weighted)
    keep it defined when the outer power is zero. The LEH loss of the footprint beams
    (a few 1e-4) applies to the moved power.
  - The bubble edge is a logistic step 0.05 mm wide (`edge`), so f is smooth in time,
    depth and pointing (ForwardDiff through `speed`, `opacity` and the beams).
  - The tip is measured from the initial wall; only the footprint plane follows a
    moving wall.
- **Defaults** (superseded 2026-10-24: speed, launch, opacity, q and collision refitted, "Hohlraum dynamics calibrated to published data" below; engineering, **M**): v_ref 0.35 mm/ns at 0.3 mg/cc (Izumi PoP 28,
  022706, 2021, 345 µm/ns; Meezan 2020 measured 321–461 µm/ns); scalings E^0.2, ρ₀^−0.2
  (Izumi 2021, ±0.03 each); delay 0.2 ns (Izumi: intercept 0.5 ± 0.1 ns for a 0.3 ns
  turn-on); p = 0.5 (Izumi Eq. 6); (R/(R − x))² stands for Izumi's ρ_bub ∝ (R/R_bub)^k
  with k = 1 (the paper's k is illegible in the extraction, **L**).
- **Checks** (*derived*):
  - Anchors: T_b = 0.290 (0.3 mg/cc) and 0.540 (0.15) at 5.5 ns, exact by construction.
  - f against a closed form: a Gaussian spot on a collimated vertical beam gives
    exp(−ρ²/a²) to 3e-3 (nr = 400).
  - Ralph et al. 2018 (6.72 mm CH high-foot, 0.6 mg/cc). Geometry **L**: R = 3.36,
    L = 11.3 mm, the retune tutorial's pointing scaled in z by 11.3/10.1. On Ralph's
    measured trajectory (0.182 mm/ns, at the wall at 2.47 ns), inner absorption reaches
    Ralph's 51% and 62% at x = 1.93 and 2.05 mm, i.e. 13.1 and 13.7 ns. That is
    consistent with his anchor if mid-peak falls there (the high-foot pulse timing is
    not in our notes, **L**). The 30° cone first meets the bubble (1% clipped) at
    r = 2.11 mm, against Ralph's 2.4 mm; the scaled pointing is a guess.
  - Swing direction. Retune-tutorial hohlraum (10.1 × 5.75 mm, `WallLossAlbedo`,
    h 0.45 / capsule_h 0.15, 0.25 ns steps), wall 0.1 mm/ns from 0.3 ns (11 keyframes),
    lips `GapClosure(closure = 5.0, width = 0.5, after = WallLossAlbedo())`, bubble
    turn-on 0 (launch 0.2 ns), footprint following the wall. Capsule flux P2/P0:

    | t (ns) | static | wall | wall + lips | wall + lips + bubble | PRE 109 Fig. 10 (N221204) |
    |---|---|---|---|---|---|
    | 2.0 | −0.173 | −0.167 | −0.167 | −0.167 | ≈0 |
    | 4.0 | −0.151 | −0.129 | −0.129 | −0.107 | −0.013 |
    | 5.5 | −0.122 | −0.075 | −0.048 | +0.051 | −0.168 |
    | 6.5 | −0.099 | −0.029 | −0.004 | +0.131 | −0.057 |
    | 8.0 | −0.092 | +0.012 | +0.035 | +0.172 | +0.090 |

    The bubble moves P2/P0 positive, as HYDRA's swing does, by +0.10 at 5.5 ns and
    +0.14 from 6.5 ns. Inner absorption is 54% (23°) and 44% (30°) at 5.5 ns, 88/79%
    at 6.5 ns, and > 99% by 8 ns. Capsule Tr falls 2–8 eV against wall + lips (the moved
    power is deposited near the LEH). The tutorial pointing is not HYBRID-E's and its foot
    is not tuned (−0.17 against ≈0), so only the direction and the late rise compare.
- **Stagnation and ridge stop (added 2026-10-23, same day).** The ballistic tip crossed a
  2.9 mm hohlraum by ~8 ns and clipped every inner beam (> 99% absorbed), so late-peak
  P2 plateaued at full inner absorption. Two options, on by default (*engineering*, **M**):
  - *Stagnation depth.* `bubble_depth` is ballistic up to the knee x_max − w, then x_max −
    w·exp(−(x_bal − x_max + w)/w) (C¹, stops at x_max). Defaults x_max = 0.80 R, w = 0.08 R:
    exactly ballistic to 0.72 R (2.07 mm at R = 2.88), past both anchor depths (1.75 / 2.01
    mm at 5.5 ns), so **the anchors stay exact** (T_b 29% / 54%, tested for exact equality
    with the ballistic tip), and near the ~2.1 mm the ballistic fit to Izumi's 0.3 mg/cc
    trajectories reaches by 6.5 ns; it stops between that and the ~2.4 mm (0.83 R) of 0.15
    mg/cc (end of the measured range, 1.2–2.4 mm). In Izumi 2023, the tip slows on
    back-pressure near the axis and is not ballistic in low-contrast pulses (stall in the
    trough, relaunch ≤ 300 µm/ns at ~6 ns; not modeled, set `speed`/`turn_on` for such
    pulses). *Not* x_max = R − 0.73 R: the ridge radius is measured at the waist, and 0.78
    mm would stop the tip before it reaches the inner beams at z_b (1.3 mm from the wall,
    Izumi 2021) and contradict the measured ballistic trajectories (1.2–2.4 mm over 4.5–6.5
    ns). Izumi 2018's point that attenuation follows the line-integrated density, not the
    displacement, is what the (R/(R − x))² density factor stands for; it is not modeled
    further.
  - *Ridge stop.* After t_c a fraction g(t) (C¹ smoothstep over t_c ± 0.2 ns) of each
    inner beam's power left by the bubble is deposited on an axisymmetric band of the
    wall component centered at the height where the beam's axis crosses r = 0.73 R
    (Izumi 2023: ridge at 2.1 mm in 2.88 mm), Gaussian in z with rms 0.5 mm, area-weighted.
    The heated ridge is taken to radiate at the wall's conversion from the wall just
    outside it (Briareus has no volume emitter; a lower ridge conversion would lower
    the waist drive further). Power is conserved and booked on the beam's `per_beam`.
    t_c = turn_on + 5.5 − 2.0 (r_cap/r_basecamp − 1) − 0.7 (E − 1) ns (*derived*, Izumi
    2018: spots separate at 5.8 ns for a 0.3 ns turn-on, ~0.4 ns earlier for 1.2× r_cap,
    5.45 ns with a 1.5× picket). The 0.5 mm band width and 0.2 ns rise are engineering
    choices (**L**). It is applied by `lpi_deposit` (new, `lpi.jl`), which
    `drive_history` (fixed and keyframed) and `Scenario` evaluations now call instead
    of `deposit`. Since the same day the beam-error sensitivities (`uncertainty.jl`:
    per-beam dual partials take g of an inner beam off the operator and put it on its
    ridge band) and `influence_matrix` see it too (checked against central
    differences of `lpi_deposit` in test/uncertainty.jl).
  - `stagnation` and `collision` are generic (`remake(m; speed, opacity, stagnation,
    collision)`), for the Next 10 fits; `stagnation = Inf`, `collision = Inf` restore
    the ballistic, ridge-free model.
  - Swing rerun (same retune hohlraum, wall + lips, bubble turn-on 0 so t_c = 5.5 ns):

    | t (ns) | wall + lips | + ballistic bubble (before) | + stagnation only | + stagnation + ridge (defaults) | PRE 109 Fig. 10 |
    |---|---|---|---|---|---|
    | 2.0 | −0.167 | −0.167 | −0.167 | −0.167 | ≈0 |
    | 4.0 | −0.129 | −0.107 | −0.107 | −0.107 | −0.013 |
    | 5.5 | −0.048 | +0.051 | +0.051 | +0.057 | −0.168 |
    | 6.5 | −0.004 | +0.131 | +0.126 | +0.128 | −0.057 |
    | 8.0 | +0.035 | +0.172 | +0.152 | +0.150 | +0.090 |
    | 9.0 | +0.065 | +0.190 | +0.169 | +0.163 | |

    Inner absorption now plateaus at ~90% (23°) and ~81% (30°) from ~7.5 ns (tip
    2.29 mm at 8 ns) instead of > 99%; Ralph 2018 measured 51–62% at peak. The late
    rise from 6.5 to 9 ns falls from +0.059 to +0.035 (the wall and lips alone give
    +0.069). The ridge stop moves P2/P0 by < 0.01: its band sits ~1–2 mm from where the
    inner beams would hit the wall, both near the waist, and converts at the wall's
    efficiency. The early offset (+0.10 above Fig. 10 at 5.5 ns, where HYDRA is at its
    minimum) is unchanged: the bubble absorbs 54 / 44% of the inner power by 5.5 ns.
    Defaults were not tuned to Fig. 10, and the pointing is not HYBRID-E's.
- **Open** (addressed 2026-10-24, "Hohlraum dynamics calibrated to published data": Ralph's 51–62% is a peak average; N221204's is 0.575). The plateau absorption (~80–90%) sits above Ralph's 51–62%: with Izumi's
  anchors fixed, any stagnation past the anchor depths leaves T_b ≲ 1% on clipped
  light. Fitting v_ref, τ_b, x_max (and t_c) to shot data (Next 10) is the next step.
  Low-contrast relaunch (Izumi 2023) is not modeled.

## Mock HYBRID-E pointing (2026-10-23) *derived*

`hybride_mock_pointing()` (`src/facilities.jl`, table `HYBRIDE_MOCK_POINTING`) is an
illustrative HYBRID-E-like pointing for the N221204 benchmark (ROADMAP "Absolute energetics on N221204"), not a facility
or shot pointing. The published papers give only relative statements about the HYBRID-E
pointing (`n221204_benchmark.md`), so the rule is an engineering choice of the usual
indirect-drive form.
- Form: one rule per cone in the frame of each quad, mirrored top to bottom, defocus 0:
  aim z, aim radius r along the quad's port azimuth (inner cones across the axis) and,
  for the outer cones, a tangential split ±s with each beam aiming to the far side of
  its quad (the 450 µm outer-quad split of Kritcher et al. 2021 is this kind). No
  top/bottom difference beyond the port pattern, and each cone's quads are rotated copies
  of each other.
- Values: the seven numbers (four z, two inner r, one outer s) on a 0.1 mm grid;
  the beams clip ~1e-5 of the laser power on the LEHs.
- Geometry (`n221204_benchmark.md` §5): 11.24 × 6.4 mm hohlraum, 3.1 mm LEHs, 1.135 mm
  capsule. Corner radius, windows and the Au/DU liner are unpublished: repository
  defaults (square corners, no windows, `gold()` wall, `ablator()` capsule). All 192
  beams (`nif_beams()`, f/20.7, CPP spots), 440 TW, inner power split evenly over the
  23.5°/30° cones and outer over 44.5°/50°, one static solve per cone fraction.
- Use: a stand-in for pointing-level studies of the benchmark (moments, CBET, LEH
  clearance), not a facility pointing. **L** for any claim that it resembles the fielded
  pointing.

## LEH shield (sombrero) (2026-10-24)

`Sombrero()` and `hohlraum(; sombrero = true)`. The defaults are
representative engineering values for a NIF symcap LEH shield, measured from a CAD model
of one (design intent, not metrology). **M**
- **Identification.** Each shield is a thin sheet around one end of the hohlraum, 18.4 mm
  across. The two shields are mirror images to the tessellation's precision.
- **Measurement.** The sheet's facets were cut by meridional half-planes every 5°, and
  height maps of both faces were ray cast on a 0.02 mm grid. The r–z sections and maps
  were checked by eye.
- **Bore and rim.** Inner radius 4.000 mm, outer 9.200 mm (vertex radii 3.9997–4.0002
  and 9.2000–9.2003 over all azimuths). Over one azimuthal sector, where the shield is
  mounted, the sheet runs on past 9.2 mm; that is not modeled.
- **Thickness.** 38.7 µm median (1–99%: 37.7–40.7 µm).
- **Axial position.** The area-mean mid-plane lies 0.623 mm (±0.005 mm, grid and
  dimple sampling) inside the inner face of the entrance-hole washer. It is toward the
  waist, beside the hohlraum wall, and not beyond the hole. That face is the LEH plane
  (`length/2`) of the matching `hohlraum` model.
- **Profile.** Flat on average: no cone or axisymmetric dimple. The sheet carries a
  hexagonal array of shallow dimples, 2.48 mm pitch, their centers sunk toward the
  waist. Dimple floor to the ridges between dimples is 0.21 mm (top face: floor −0.094,
  ridge points +0.119 mm from the mean of that face). The faces stay
  within ±0.12 mm of the mid-plane (0.5–99.5%). Two mounting recesses go ~0.3 mm deeper.
  The rise from a dimple's floor is close to linear in the distance from its center:
  ~4° near the floor, ~12° toward the ridges.
- **Model.** `Sombrero` keeps the mean plane: a flat annulus (`annulus_mesh`, exact
  `DiskSurface`). The dimples are not axisymmetric; they move the true faces by up to
  ±0.12 mm, about the size of a hardware tolerance in `clearance`.
- **Material.** The CAD carries none. The default is `BLACK` (engineering default).
- **Radiation.** It has no effect on the radiation. The shield sits beside the wall, 0.623
  mm inside the LEH plane, so the hohlraum wall hides the hole from it. Its view-factor rows
  and columns with the hohlraum are exactly zero (`test/sombrero.jl`), `closed = true`
  balancing leaves those rows alone, and the hohlraum's solution is unchanged to 1e-12.
  One-sided, facing away from the waist: light from the port side is booked on its elements,
  light from the waist side (e.g. backlighter beams) as `sombrero_*_outside`.
- **Drive-beam clearance** *derived*. The test's cone-averaged NIF pointing (Build-A-
  Hohlraum example, 10.1 mm hohlraum, example cone aims) was run as 192 beams. The
  core 1e-4 envelope clears both shields by ≥ 2.09 mm in the worst case (outer 50°;
  44° 2.50, 30° 2.41, 23° 3.25), and the ray trace deposits no core power on them. With
  the engineering `nif_wings()` the envelope clearance is ≥ 1.24 mm. The wings'
  10⁻³ scatter floor puts ~1e-4 of each beam on a shield, so every beam is `:review`
  (above 1e-5) at 1.2e8 W/cm². Drive beams
  reach the shield only past the hohlraum's outside, so it matters for beams aimed
  outside the hohlraum (backlighters), not for the drive.

## Hohlraum dynamics calibrated to published data (2026-10-24) *derived*

ROADMAP "The P2 swing" and "Shot-data calibration". Each piece of the hohlraum dynamics (wall, LEH lips, bubble, ridge
stop, two-zone CBET) is refitted to a separate published anchor. PRE 109 Fig. 10 is not used
in any fit; it is the blind check at the end. Full texts are the accepted manuscripts.

### Anchors extracted

| Anchor | Source | Form | Value | Stored |
|---|---|---|---|---|
| Dynamic VF hohlraum radius vs t (N170827) | Young et al., PoP 27, 082702 (2020) Fig. 3a | raster, axes from frame and gridlines, overlay checked | flat to 0.9 ns, in at 0.185 mm/ns, plateau 0.43–0.48 mm in through the peak, 0.55 at 7.5 ns | `young2020_radii()` |
| LEH radius vs t (same model) | Young Fig. 3a | same | 1.83 → 1.57 mm (0.857 r₀) linearly by ~6 ns; only the end value is SXI data | same |
| HYBRID-E wall edge at 6 ns | PRE 109 Fig. 1a (HYDRA) | raster, axis/waist/outline calibration, 1 px = 0.018 mm | waist 0.38 ± 0.04 mm in (laser-absorption edge); bubble tip ≈1.85 ± 0.1 mm (N210808), ≥ 1.55 (N221204, under its rays) at z ≈ 2.9–3.5 mm | `pre109_fig1a_wall()` |
| Bubble tip x_bub(t), 7 shots | Izumi et al., PoP 28, 022706 (2021) Fig. 8a/b | vector graphics (exact) | 1.38–2.28 mm over 4.65–6.19 ns; ±40 µm | `izumi2021_bubble()` |
| Inner-power fraction absorbed in the bubble | Izumi 2021 Fig. 10a (HYDRA circles, Eq. 6 line) | vector | 0, 0.020, 0.230, 0.550 at 3, 4, 5, 6 ns (basecamp, 0.3 mg/cc) | `izumi2021_absorption()` |
| Spot separation (ablator–wall collision) | Izumi et al., RSI 89, 10K111 (2018) | text | 5.8 ns basecamp (turn-on 0.3), 5.45 ns at 1.5× picket, ~0.4 ns earlier with the 1.18× capsule | docstring, tests |
| Peak-averaged inner absorption | Ralph et al., PoP 25, 082701 (2018) | text | 51 ± 2% to 62 ± 2% (CH high foot, 6.72 mm, 0.6 mg/cc, over the peak FWHM) | notes |
| Cone fraction for zero flux P2 | Ralph 2018 (VisRad at the measured peak wall position) | text | 18.5%; peak F_in 0.34–0.38 → F' 0.14–0.18 (Fig. 5) | notes |
| Bubble density | Izumi 2021 Eq. 7 | text | ρ_bub = S ρ₀ R/(R − x), S ≈ 3 ∝ ρ₀^-0.1 (k = 1, settling the "illegible" exponent) | `density_exponent` 1.8 |

Digitization: Izumi 2021's figures are vector graphics in the accepted manuscript, so
markers and lines are read from the PDF paths (axes from the tick marks; exact to the
drawing). Young Fig. 3a and PRE 109 Fig. 1a are rasters: axes from frames, gridlines and
the drawn outline, stroke centers per column, overlays checked by eye.

### Fits (each against a separate anchor)

| Parameter | Before | Now | Uncertainty | Anchor |
|---|---|---|---|---|
| Wall speed | 0.10 mm/ns | 0.185 mm/ns | 0.175–0.195 | Young Fig. 3a |
| Wall start | outer turn-on | turn-on + 0.8 ns (fit 0.78) | 0.68–0.85 ns | Young Fig. 3a |
| Wall saturation | none | 0.45 mm (default); 0.38 for HYBRID-E | 0.445–0.455 (±0.03 reading); HYBRID-E ±0.04, upper bound | Young; PRE 109 Fig. 1a waist |
| Bubble speed | 0.35 mm/ns | 0.36 mm/ns | 0.347–0.372 | Izumi Fig. 8a, 28 points, χ²/dof 0.89 |
| Bubble launch delay | 0.2 ns | 0.25 ns | 0.07–0.42 (correlated with speed) | Izumi Fig. 8a |
| Density exponent q | 1.862 (two-anchor fit) | 1.8 (physics) | — | Izumi Eq. 7: ρ_bub² ∝ (Sρ₀)², S ∝ ρ₀^-0.1 |
| Opacity τ_b | 0.144 (fit to T_b 29%/54%) | 0.111 | 0.100–0.124 stat.; 0.085–0.158 with plane ±0.2 mm and launch | Izumi Fig. 10a HYDRA absorbed fraction |
| Collision t_c | turn_on + 5.5 − 2.0(s − 1) − 0.7(E − 1) | turn_on + 5.5 (R − r_cap)/1.967 mm · E^-0.16 | ±0.15 ns (frame spacing) | Izumi 2018: 5.8, 5.45 (1.5× picket), −0.4 ns (1.18× capsule; the gap gives −0.44) |
| HYBRID-E picket fluence `picket_ratio` | 1 | 1.1 (example) | 0.8–1.5 | 1.5× basecamp energy (DN24 vs Izumi Fig. 3a) × 0.73 for split outer quads; Fig. 1a N210808 tip 1.85 ± 0.1 mm at 6 ns |
| LEH lips (example) | 5.0 ns step over 0.5 ns to 0.84 r₀ | smoothstep from −0.1 ns over 7.85 ns to 0.857 r₀ | half-closed 3.0–3.8 ns (timing mapping) | Young Fig. 3a LEH radius (end value SXI; ramp assumed) shifted 0.8 ns |

Checks of the fits:
- *Wall.* Young's radius is flat to 0.9 ns, falls at 0.185 mm/ns to 0.45 mm in by 3.3 ns
  and then stays at 0.43–0.48 mm through the peak: the motion happens in the foot, not at
  a constant speed through the pulse. The fit's residuals are < 0.03 mm over 0.5–6.75 ns
  (`test/dynamics.jl`). On N221204 the old wall (0.10 mm/ns from the turn-on) was
  0.56 mm in at 6 ns and 0.86 at 9 ns; the new one is 0.38 from 3.3 ns on.
- *Bubble tip.* The seven shots' residuals to the global fit (E^0.2 ρ₀^-0.2 scalings,
  one speed, one launch) are −0.09 to +0.07 mm, rms 0.036 mm, inside the ±0.04 error bars.
  Per shot, the slopes are 0.33–0.43 mm/ns and the intercepts 0.14–1.02 ns.
- *Absorbed fraction.* The opacity is degenerate with the bubble's plane and launch: with
  the plane where the proxy pointing's 50° axes reach the wall (z = 2.59 mm instead of
  Izumi's 3 mm) τ_b = 0.086. The absorbed fraction against time is the same in every
  variant (0.03 / 0.20 / 0.36 / 0.56 at 4 / 5 / 5.5 / 6 ns, against HYDRA's 0.02 / 0.23 /
  – / 0.55 and Izumi's Eq. 6 0.04 / 0.20 / 0.32 / 0.56). The old defaults gave 0.09 / 0.32
  / 0.50 / 0.69 on the proxy geometry (plane 2.59 mm): too much, too early.
- *Collision.* Izumi 2018's 1.18× capsule (outer radius 1.064 mm against 0.908) narrows
  the capsule–wall gap by 8%, which at constant closing speed (0.36 mm/ns) brings the
  collision 0.44 ns earlier: the measured "roughly 400 ps". The 1.5× picket's 5.45 ns
  gives the energy exponent −0.16. HYBRID-E's gap (3.2 − 1.135 = 2.065 mm) is 5%
  *wider* than the basecamp's, so t_c = 0.4 + 5.77 · 1.1^-0.16 = 6.09 ns. The old
  extrapolation in capsule radius alone gave 5.41 ns.
- *Ralph 2018.* Energy-weighted over N221204's peak FWHM (4.62–8.76 ns), after CBET, the
  calibrated bubble absorbs 0.575 of the inner power (old defaults 0.50, E = 1.5 0.65),
  inside Ralph's 51–62%. The platforms differ (CH high foot, 0.6 mg/cc, 14 ns pulse), so
  this is a plausibility check, not a fit. The earlier "85% against 51–62%" compared an
  instantaneous value at 8 ns with a peak average. Ralph's 18.5% cone fraction for zero
  flux P2 at the measured peak wall position is not reproduced: it needs the high-foot
  pointing and the measured wall radius, neither of which is published.
- *HYBRID-E bubble (PRE 109 Fig. 1a, 6 ns).* N210808's tip is 1.85 ± 0.1 mm from the
  wall at z ≈ 2.9–3.5 mm; N221204's is ≥ 1.55 mm (its denser inner-beam rays are drawn
  over the plasma). The calibrated model with `picket_ratio` 1.1 puts it at 1.96 mm, in the
  plane z_b = 3.49 mm.

### Two-zone CBET: the 7.4 keV LEH temperature

`n221204_cbet()`'s Te(t) reaches 7.44 keV at 8 ns, above every published LEH value. Refit
to the same Fig. 9b data with Te ≤ 5.5 keV (`calibrate_cbet` with an upper bound,
16 min on 4 threads: `examples/n221204_cbet_calibration.jl` with `bounds[:Te] = (300, 5500)`):

| Fit | Fig. 9b rms | Te (keV, 1 → 8 ns) | multipliers out / in, Mach | 44° / 50° at 5.6 ns (Fig. 9a 0.43 / 0.64) | CF at Δλ 0 | Δλ null | CF at 1.8 Å, 5.6 ns |
|---|---|---|---|---|---|---|---|
| `n221204_cbet()` | 0.0163 | 2.79 → 7.44 | 0.30 / 1.23, 1.23 | 0.41 / 0.68 | 0.090 | 1.20 Å | 0.461 |
| A: Te ≤ 5.5 keV | 0.0160 | 1.75 → 5.01 | 0.07 / 0.41, 0.91 | 0.47 / 0.61 | 0.281 | 0.50 Å | 0.499 |
| B: A + ne/nc(t) | 0.0115 | 4.87, 1.87 … 5.50, 4.88 | 0.04 / 1.09, 0.33 | 0.48 / 0.58 | 0.333 | 0.05 Å | 0.513 |

- Fig. 9b does not set the absolute Te: A fits it as well with every Te ~0.67× lower,
  the outside multiplier 4× smaller and the inside one 3× smaller. The 7.4 keV is the
  overall gain scale and the balance between the zones, not a measured temperature.
  The rise of Te with time (×2.7–2.9 from 1 to 8 ns) is what the data constrain.
- What separates the fits is the Δλ response, which is not fitted. The reference puts
  the null at 1.2 Å, inside NC24's 1–1.55 Å for HYBRID-E. A moves it to 0.5 Å. A
  matches HYDRA's N210808 (1.8 Å) post-CBET cone fraction better (0.50 against ~0.55–0.6
  in Fig. 1e; the reference 0.46), i.e. it halves the excess Δλ sensitivity.
- B fits Fig. 9b better only with an unphysical density history (ne/nc 0.15 at 1 ns,
  0.01 at 4.75–5.5 ns) and moves the null to 0.05 Å: rejected.
- Not adopted: `n221204_cbet()` is unchanged. Constraining Te trades one unfitted check
  (the null) for another (N210808). A joint fit to N210808 / N220919 with the null as a
  constraint (ROADMAP "Shot-data calibration") is what can fix the temperature scale. Done: "Joint multi-shot
  CBET calibration (2026-10-24)" below. Te ≤ 5.5 keV costs the three-shot fit little; the
  null in 1–1.55 Å is what the cross-shot data resist.

### N221204 benchmark rerun: the blind check against Fig. 10

`examples/n221204.jl` (default `dynamics = calibrated`; `old` reruns the 2026-10-23
settings), calibrated two-zone CBET, mock HYBRID-E pointing, standard mesh (h 0.28 /
capsule_h 0.16, 31 min for the five runs on 4 threads). Settings: wall from 1.2 ns at
0.185 mm/ns stopping 0.38 mm in; lips to 0.857 r₀ along Young's ramp shifted 0.8 ns
(closure −0.1 ns, width 7.85 ns); bubble defaults with `picket_ratio` 1.1 and
`capsule_radius` 1.135 (t_c 6.09 ns). `picket_ratio` was set from the picket energy and quad
splitting before the comparison was run.

| t (ns) | static | +wall | +lips | +bubble (now) | +bubble (2026-10-23) | Fig. 10 |
|---|---|---|---|---|---|---|
| 2.0 | −0.057 | −0.034 | −0.031 | −0.031 | −0.032 | −0.001 |
| 4.0 | −0.072 | −0.026 | −0.017 | −0.016 | −0.029 | −0.018 |
| 4.5 | −0.088 | −0.053 | −0.040 | −0.034 | −0.043 | −0.092 |
| 5.5 | −0.258 | −0.266 | −0.239 | −0.166 | −0.086 | −0.169 |
| 6.5 | −0.227 | −0.233 | −0.202 | +0.029 | +0.031 | −0.045 |
| 7.0 | −0.207 | −0.210 | −0.179 | +0.085 | +0.093 | +0.046 |
| 8.0 | −0.177 | −0.175 | −0.145 | +0.123 | +0.168 | +0.095 |
| 8.5 | −0.220 | −0.229 | −0.196 | +0.126 | +0.176 | +0.109 |

| Swing metric (full model) | 2026-10-23 (coarse mesh) | now (standard) | Fig. 10 |
|---|---|---|---|
| minimum | −0.164 at 5.00 ns | −0.168 at 5.25 ns | −0.173 at 5.58 ns |
| upward zero crossing | 6.25 ns | 6.26 ns | 6.67 ns |
| at 8 ns | +0.162 | +0.123 | +0.095 |
| mean over 5.25–8.5 ns | +0.054 | +0.024 | −0.001 |
| rms against Fig. 10, 2–8.5 ns | 0.051 | 0.035 | — |

The coarse mesh (0.8 / 0.4) reproduces these within 0.006. Capsule Tr peaks at 301.8 eV
(5.5 ns) with everything on, 307.2 with wall and lips (drive Trad 313 eV).

**Why the swing changed.**
1. *Wall.* Saturating at 0.38 mm takes +0.09 off the late P2 (the +wall column now
   stays at the static value after 5.5 ns; before, it rose to −0.083 at 8 ns). This is
   most of the late improvement (+0.168 → +0.123 at 8 ns).
2. *Collision.* The ridge stop now starts at 6.09 ns instead of 5.41, so at 5.5 ns the
   inner light still reaches the wall where it is not absorbed: the 5.5 ns value drops
   from −0.086 to −0.166 and the minimum deepens to Fig. 10's level.
3. *Bubble.* On HYBRID-E's geometry the refitted bubble absorbs about as much as before
   by 5.5 ns (0.22 against 0.17) but more later (0.57 at 6.5 ns against 0.47, 0.92 at
   8 ns against 0.84). With the ridge on from ~6.2 ns this sets the rise from −0.17 to
   +0.03 between 5.5 and 6.5 ns, 0.4 ns ahead of HYDRA's.
4. *Lips.* Closing along Young's ramp instead of a step at 5.0 ns moves P2 by < 0.01
   from 4 ns on (sensitivity run `late_lips`: minimum −0.169 at 5.0 ns, crossing 6.25,
   +0.118 at 8 ns, coarse mesh).

**Sensitivity to the least constrained input** (coarse mesh, full model). The HYBRID-E
picket fluence `picket_ratio` is 1.1 from the picket energy (1.5×, by eye) and split outer
quads (−27%); PRE 109 Fig. 1a's N210808 bubble tip (1.85 ± 0.1 mm at 6 ns) alone gives 0.8
(0.6–1.05). E = 1.5: minimum −0.138 at 5.25 ns, crossing 5.99 ns, +0.121 at
8 ns. E = 0.8: −0.184 at 5.5 ns, 6.62 ns, +0.109 (rms 0.024). The crossing time follows
the bubble speed: 0.6 ns per factor 1.9 in E, i.e. per ±7% in speed. The default is
kept at 1.1, the value set before the comparison; 0.8 matches Fig. 10 better and is
consistent with Fig. 1a, but choosing it now would be tuning to the blind check.

**N220919 and N210808 proxies** (coarse, `shot` argument; same mock pointing, CBET
calibration and hohlraum; N220919 on N221204's pulse with Δλ 2.5 Å; N210808 with the rise
0.2 ns and the peak end 0.45 ns earlier and Δλ 1.8 Å):

| | minimum | crossing | at 8 ns | rms 2–8.5 ns |
|---|---|---|---|---|
| N220919 now (old) | −0.140 at 5.25 (−0.139 at 5.0) | 6.23 (6.14) | +0.118 (+0.165) | 0.033 (0.053) |
| Fig. 10 N220919 | −0.145 at 5.6 | 6.66 | +0.094 | |
| N210808 now (old) | −0.092 at 5.25 (−0.095 at 5.0) | 6.09 (5.67) | +0.119 (+0.166) | 0.061 (0.077) |
| Fig. 10 N210808 | −0.141 at 5.5 | 6.68 | +0.084 | |

N220919 behaves like N221204. N210808's proxy has a positive foot (+0.03 at 4 ns against
−0.05) and a shallow minimum: its pulse and foot cone fraction are only approximated
(N221204's DN24 powers, shifted), and the CBET calibration is N221204's (its Δλ
sensitivity is ~2.5× HYDRA's across shots, ROADMAP "Shot-data calibration").

**Remaining discrepancies and likely causes.**
- Crossing ~0.4 ns early: the HYBRID-E bubble speed (above), and the bubble plane:
  z_b follows the mock pointing's 50° axes (3.49 mm), which sets how much of the inner
  beams the bubble clips; the absorbed-fraction calibration was done on a 5.75 mm
  hohlraum with Izumi's plane.
- 4.5 ns 0.06 above Fig. 10: DN24's delivered input cone fraction (0.13–0.14) is below
  Fig. 9b's (0.17–0.18) in the rise, so the CBET gives less transfer there.
- Foot ~0.03 below: mock pointing and the foot cone fraction (read by eye).
- Late +0.02–0.03 high: the capsule does not shrink (Young's ablation front is at 0.83 mm
  at 5.5 ns, falling fast after 6 ns), no laser degradation or backscatter.
- Every dynamic piece is still one effective parameter set from other platforms (HDC
  SymCap in 5.75 mm Au, N170827); only the wall saturation and the bubble-tip check come
  from HYBRID-E itself, and both from the HYDRA picture of Fig. 1a.

## Computed LEH closure (2026-10-24) *calibrated*

ROADMAP "The P2 swing" ("a computed lip closure") and the late-drive suspects of ROADMAP "Absolute energetics on N221204" ("Wall
re-emission and the capsule radius": the LEH closing faster after the laser turns off,
the energy stored in the hohlraum plasma). The prescribed lips (`GapClosure`, Young's
ramp shifted 0.8 ns) stay available and remain the example's default; the computed one
is `hohlraum(; leh_closure = :computed)`. Dante 1 and Fig. 10 are checks, never fit.
Confidence: **H** primary text or vector data, **L** our construction.

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| Effective LEH radius vs t, 5 basecamp shots (5.75 × 10.13 mm Au, 3.37 mm LEHs, r₀ 1.685 mm, ~400 TW from 4.2 to 6.05 ns) | 1.54–1.61 mm at 3.6 ns (0.92–0.95 r₀); 1.51–1.58 at 4.0–4.2; 1.40–1.50 through the peak (mean 1.44, 0.855 r₀; shot means ±0.035); after the laser turns off (6.05 ns) N171010 1.46 → 1.44 / 1.42 / 1.43 at 7.07 / 7.44 / 7.75 ns, N171212 1.40 at 6.6–6.9 | Chen et al., PoP 27, 022702 (2020) Fig. 27, vector paths (`chen2020_leh()`, `data/chen2020_fig27.csv`); analysis error 10–20 µm (Sec. II.2b) | H |
| LASNEX on the same shots (images analyzed as the data) | closes at 0.07–0.11 mm/ns after laser-off; rms against the 31 points 0.045 (f = 0.15) to 0.056 mm (f = 0.03), after laser-off 0.09–0.12 mm, all too closed | same | H (rms: L) |
| Mechanism | a dense ridge forms on the rise where the outer beams cross the lip blow-off; "once the ridge is formed the effective LEH radius is fairly stable, closing only very slowly, until the laser turns off"; then (in LASNEX) "the ridge diffuses inward"; "heating of the plasma inside the LEH region by the lasers keeps the LEH open"; the edge is set by hydrodynamics, not opacity | Chen 2020 Sec. II.2a, Fig. 5b–6 | H |
| Closure at peak, other platforms | 1.83 → 1.57 mm (SXI end value, depth 0.26 mm; 3.64 mm LEH) | Young et al., PoP 27, 082702 (2020) Fig. 3a | H |
| | ~0.84 r₀ of 1.55 mm (depth 0.25), ≤ 4% variation through peak; closes in the low-power part, reopens at the rise | Schneider et al., PoP 22, 122705 (2015) ("Closing entrance holes" above) | H |
| SXI | time-integrated (3–5 keV clear aperture, 870 eV channel); Dante's correction f (R_LEH/R_CA)² ≈ 1.15 (+10 eV at 300 eV) | Schneider et al., HTPD 2012 (OSTI 1043636); Schneider et al., RSI 81, 10E538 (2010) (OSTI 1009651) | H |
| Hohlraum gas | He 0.3 mg/cc in N221204; 4.5% of the incident laser energy heats the gas | n221204_benchmark.md; PRE 109, 025204 App. 3 | H |
| Gas and window heat capacity | 25–30% of the picket energy in low-fill hohlraums, a small fraction of the peak's | Callahan et al., PoP 27, 072704 (2020) Sec. II | H |
| Cooling phase | the 0-D model misses the Dante data "when the laser power drops and the hohlraum starts to cool ... the radiation lost out of the LEH becomes more important" | Callahan 2020 Sec. II | H |

OSTI searches ("laser entrance hole", "hole closure", "LEH closure", hohlraum filling)
found no other time-resolved LEH-radius data; Hall et al. 2017 (OSTI 1361586) and the
SXI papers are time-integrated. No published lip temperature or blow-off speed.

### Model (`LipClosure`, `lip_depth`; `drive_history` per lip) *L*

The lip's plasma edge sits a depth d inside r₀ and moves as

    dd/dt = (1 − h) v(T) + h (depth − d)/hold_time,    v(T) = speed (T/100 eV)^½,

T the radiation temperature on the lip (the sound speed of x-ray-ablated gold, c_s ∝
T^½ at fixed Z̄), h the smooth step of the laser intensity through the hole (beams whose
port is on that side, after `lpi`, over π r₀²) between `hold` = (6, 17) TW/mm². Below
`hold[1]` (the foot, and as the power falls) the blow-off expands freely; above it the
edge is held at `depth`, from either side (a lip that closed further in a long foot is
pushed back, Schneider's reopening). The lip is meshed in `nr_lip` = 6 equal-area rings
down to 0.6 r₀; each ring's closed fraction is the share of its area outside r₀ − d,
so the open area is exact (to 1e-10 against the chord polygon, `test/leh_closure.jl`).
`drive_history` integrates d per lip with the trapezoidal mean of v and h over each step
(predictor: last step's T; corrector: the T of the current step, re-solved until d moves <
1e-5 mm), exactly in time for the hold's relaxation, and records `lip_radius`.

### Calibration (`examples`: none; scripts in the session scratchpad) *calibrated*

The lip temperature comes from the model itself: a coarse basecamp history (h 0.8 / 0.4;
10.13 × 5.75 mm, 3.37 mm LEHs, 0.908 mm capsule, `WallLossAlbedo` wall, `nif_conversion()`,
the Izumi-trends pointing, each shot's Fig. 27 power split by N170827's cone fraction from
`young2020_pulses()`, prescribed lips): lip T 120–125 eV in the foot, 170–180 at 3–3.5 ns,
275–288 through the peak, 226 / 170 / 138 eV at 6.25 / 7.0 / 7.75 ns (capsule Tr peak
292–294 eV); the five shots' lip temperatures differ by < 3 eV.

Grid fit (speed 0–0.06 mm/ns, depth 0.18–0.32 mm, hold_time 0.05–0.8 ns) of the 31 GLEH
points at hold = (6, 17) TW/mm² (107–303 TW on two 3.37 mm LEHs, the rise):

| Parameter | Fit | Range | Note |
|---|---|---|---|
| speed (at 100 eV) | 0.022 mm/ns | 0.010–0.030 (rms within +10%); 0.020–0.024 over the hold variants | ~0.55 c_s of gold at 100 eV (Z̄ ≈ 30, 39 km/s; *L*) |
| depth | 0.240 mm | 0.235–0.245 | Young 0.26 (3.64 mm LEH), Schneider 0.25 (3.1 mm) |
| hold_time | 0.40 ns | 0.2–0.6 (degenerate with `hold`) | the rise's 0.1 mm drop over ~0.5 ns |
| hold | (6, 17) TW/mm² | (3–10, 12–22) give rms 0.032–0.034 | free at 3.6 ns (110 TW), held at peak |

rms 0.033 mm over the 31 points (rise 0.034, peak 0.034, after laser-off 0.026), mean
−0.001; LASNEX 0.045–0.056 (after laser-off 0.09–0.12). The floor is the shot scatter:
at peak N171212 sits at 1.40 and N180117 at 1.49 mm with the same pulse and
temperatures. After the laser turns off the model's edge moves in at 0.02–0.03 mm/ns
(T falling from 280 to 140 eV), measured ≤ 0.02 (N171010: −0.03 mm in 1.7 ns), LASNEX
0.07–0.11: **the measured LEH does not close faster after laser-off**, and the model,
fitted to the foot and the peak, agrees with that within the analysis error. The
integrated `drive_history` on the same basecamp with `leh_closure = :computed`
reproduces the offline integration to 2e-5 mm and the data to rms 0.033 mm.

A depth of 0.24–0.26 mm on 3.1–3.64 mm LEHs is also a fraction of 0.84–0.86 r₀; the three
platforms cannot tell an absolute depth from a fraction. The depth is the model's form
(the ridge stands where the beams meet the blow-off, which does not scale with the
hole). On the Zylstra 2022 pair (§5 of "Wall material, conversion and capsule
energetics") a common depth gives ×1.143 against ×1.122 for a common fraction, the
data ≥ ×1.19: the right direction, not enough.

The thermal-band radius is smaller than the GLEH (3–5 keV, > 5 keV) edge (longer mean
free path at higher photon energy; Chen Sec. II.2); calibrating on GLEH, like Young's
SXI, makes the hole somewhat too open for the x-ray loss (unquantified, *L*).

### N221204 (check, not fit)

`examples/n221204.jl 0.8 0.4 full energetics=new` (coarse; the full model with wall
re-emission), prescribed lips (Young's ramp shifted 0.8 ns, to 0.857 r₀) against
`lips=computed`. Edge radius of the lower LEH (mm; prescribed: the radius of equal open
area):

| t (ns) | 2.0 | 4.5 | 5.5 | 8.0 | 8.5 | 9.0 |
|---|---|---|---|---|---|---|
| prescribed | 1.513 | 1.415 | 1.375 | 1.328 | 1.328 | 1.328 |
| computed | 1.513 | 1.420 | 1.320 | 1.310 | 1.310 | 1.306 |

The computed lip follows the prescribed ramp through the foot (free expansion at
120–180 eV), is held 0.24 mm in from the rise on (prescribed, still closing to 8 ns), and
stays held through the droop: the intensity through each hole is 17.8 TW/mm² at 8.6 ns
and 8.8 at 9.2 ns against the hold's 6–17, so the edge is free again only after ~9 ns.

Capsule Tr / Dante 1 over the initial / open LEH area (eV):

| t (ns) | prescribed | computed | measured Dante 1 |
|---|---|---|---|
| 5.0 | 260.9 / 248.4 / 261.8 | 261.5 / 244.0 / 261.9 | 237.9 |
| 5.5 | 279.9 / 264.4 / 280.7 | 280.9 / 259.5 / 281.1 | 263.1 |
| 6.0 | 285.8 / 271.0 / 289.4 | 286.9 / 267.0 / 290.1 | 275.5 |
| 7.0 | 297.8 / 283.5 / 305.5 | 298.7 / 281.5 / 306.2 | 291.3 |
| 8.0 | 310.6 / 294.6 / 318.2 | 311.5 / 293.2 / 319.0 | 310.1 |
| 8.5 | 299.5 / 283.1 / 305.8 | 300.4 / 281.9 / 306.6 | 311.3 |
| 9.0 | 277.0 / 260.4 / 281.3 | 278.0 / 259.0 / 282.2 | 301.5 |
| LEH / capsule share (%) | 25.5 / 11.7 | 24.8 / 11.8 | HYDRA 25.2 / 16.2 |
| Fig. 10: min (ns), zero, 8 ns, rms | −0.181 (5.25), 6.24, +0.073, 0.041 | −0.176 (5.25), 6.23, +0.073, 0.041 | −0.170 (5.6), 6.67, +0.095 |

+1 eV through the peak (the hole is 0.02 mm smaller from 5.5 ns), the LEH share 0.7
points lower (HYDRA's 25.2% lies between), Fig. 10 within 0.005. **The late gap is not
the LEH**: 8.0 → 8.5 ns the model still falls 11 eV (capsule), Dante rises 1. As a bound, to
hold the drive the losses must fall by ~4 × 3.5% ≈ 14% (T^~3.5 at fixed source); with
the LEH 25% of the loss that is −56% of the LEH area, an edge 0.44 mm further in within
0.45 ns (~1 mm/ns): 25× the fitted free speed at 300 eV, 10× LASNEX's post-laser-off
closure, and against Chen's data, which show no faster closure after the laser turns
off. Even LASNEX's 0.1 mm/ns from 8.05 ns would add ~1.5 eV by 8.5 ns.

Standard mesh (0.28 / 0.16), computed lips: capsule Tr / Dante 1 initial / open area
261.6 / 244.4 / 262.2 at 5.0 ns, 311.4 / 293.4 / 319.2 at 8.0, 300.3 / 282.0 / 306.7 at
8.5, 278.0 / 258.9 / 282.2 at 9.0 (prescribed, "Wall re-emission and the capsule radius":
261.0 / 248.8 / 262.3, 310.5 / 295.0 / 318.6, 299.4 / 283.4 / 306.1, 277.0 / 260.5 /
281.4); LEH / capsule share 24.9 / 11.8% (25.6 / 11.7); Fig. 10 −0.183 at 5.25 ns,
crossing 6.22, +0.077 at 8 ns, rms 0.042 (−0.189, 6.23, +0.077, 0.042). Coarse agrees
within 0.3 eV and 0.007 in P2.

### Plasma-stored energy (estimate; not implemented) *L*

The quasi-static radiosity converts each step's laser power at once; energy held in the
hohlraum plasma (He gas, gold blow-off) and given back after the droop has no store.

- **He gas** (the one published inventory). 0.3 mg/cc in 0.356 cm³ (6.4 × 11.24 mm less
  the capsule) is 0.107 mg, 1.6e19 ions; fully ionized at T_e 2.5 / 3 / 5 keV (the LEH
  T_e ≤ 5.5 keV of the two-zone CBET calibration; Chen's bubble interior 3–6 keV) it
  contains 4.5 n T_e = 29 / 35 / 58 kJ. PRE 109 App. 3 books 4.5% of the incident energy
  (~92 kJ) into the gas over the pulse, consistent with a store of that size drained into
  the wall as it goes. He radiates nothing to speak of (bremsstrahlung ~0.5 TW at
  n_e 9e19 cm⁻³); it gives its energy to the gold by conduction, at most the
  flux-limited f n_e T_e v_e = 5e12 W/cm² (f 0.05, 3 keV) over ~2.75 cm², 14 TW, so
  ~35 kJ drains over ~2.5 ns: ≲ 12 TW of x-rays after conversion against the ~100 TW
  the x-ray source loses from 8.05 to 8.5 ns (laser 448 → ~320 TW × 0.89 × 0.87), worth
  ~1–2 eV at 8.5 ns. Not the gap.
- **A fast reservoir** (the laser-heated gold corona, which converts the laser at
  Te 3–5 keV) would act as a lag τ of the source with E = P τ. Sensitivity
  (`examples/n221204.jl ... lag=<ns>`, the delivered power replaced by its first-order
  lag; coarse, prescribed lips): capsule Tr / Dante 1 open area (eV), Fig. 10.

  | lag (ns) | stored at 395 TW coupled (kJ) | 4.5 ns | 5.0 ns | 8.0 ns | 8.5 ns | 9.0 ns | Fig. 10 min (5.25 ns), zero, 8 ns, rms |
  |---|---|---|---|---|---|---|---|
  | 0 | — | 216.5 / 222.1 | 260.9 / 261.8 | 310.6 / 318.2 | 299.5 / 305.8 | 277.0 / 281.3 | −0.181, 6.24, +0.073, 0.041 |
  | 0.08 | 32 | 209.9 / 216.0 | 254.9 / 256.6 | 310.4 / 318.1 | 302.7 / 309.4 | 280.8 / 285.4 | −0.169, 6.24, +0.073, 0.041 |
  | 0.15 | 59 | 204.8 / 211.1 | 249.5 / 251.9 | 310.2 / 317.9 | 304.8 / 311.7 | 284.5 / 289.4 | −0.156, 6.24, +0.074, 0.042 |
  | 0.25 | 99 | 198.9 / 205.5 | 242.3 / 245.4 | 309.7 / 317.4 | 306.6 / 313.7 | 289.3 / 294.6 | −0.138, 6.23, +0.075, 0.044 |
  | 0.15, computed lips | 59 | 204.5 / 210.9 | 249.9 / 251.8 | 311.0 / 318.6 | 305.7 / 312.5 | 285.5 / 290.3 | −0.152, 6.23, +0.074, 0.042 |
  | measured Dante 1 / Fig. 10 | | 203.4 | 237.9 | 310.1 | 311.3 | 301.5 | −0.170 (5.6), 6.67, +0.095 |

  A lag of ~0.15 ns (~60 kJ at 395 TW coupled) brings Dante 1 at 8.5 ns to the measured
  311 eV, lowers the early-peak excess at 5 ns from +24 to +14 eV (open area), puts the
  4.5 ns readings around the data (201.7 / 211.1 eV initial / open area against 203.4),
  and moves Fig. 10's minimum to −0.156 (Fig. 10 −0.170)
  without moving the crossing. One time constant explaining both the late rise and the
  early lag (the model's rise was ~0.15 ns ahead) points at a fast plasma store of
  ~60 kJ, the size of the gas inventory but needing a ~0.15 ns return, which only the
  radiating gold corona has. Its energy content and the conversion's response time are
  not published (PRE 109 gives integrals only; Callahan 2020's 0-D model notes the
  cooling phase as where it fails), so the lag stays a sensitivity in the example, not
  a model. The energy is conserved: the run's delivered energy drops by the 19 kJ the
  lag moves past the last step. 9.0 ns stays 7–12 eV low for lags of 0.15–0.25 ns (a
  slower store, the gas, would add there).

### What remains

- A measured HYBRID-E LEH history (3.1 mm, 440 TW, a longer peak) to test the hold
  depth's independence of the LEH size and power; the thermal-band radius against the
  GLEH edge.
- The hohlraum plasma's energy content and the conversion's response time (a HYDRA
  inventory, or Dante's delay behind a step in the laser power) to turn the lag into a
  model.
- Shot scatter at peak (±0.035 mm across the five shots) is not modeled: N180117's split
  quads sit 0.05 mm more open, N171212's 1.5× picket 0.05 mm more closed.

## The swing and the bubble without new data (2026-11-18) *derived*

ROADMAP items 10 (the P2 swing and the gold bubble), 36 (constrained re-pointing) and 12
(the notes fix, done in "N221204 late-time drive ▸ What remains"). What can move on the
swing without bubble-tip histories or HYDRA's multipliers: the current state remeasured,
the computed lips against the band they were fitted in, the published bubble inputs not
yet used, and what the literature ties bubble growth to. PRE 109 Fig. 10 is the blind
check throughout; nothing below is fitted to it and no default or calibrated constant
changed. Literature for §4: `reference/litreview/bubble_drivers.md`. Runs:
`examples/n221204.jl 0.8 0.4 full calibrated x n221204 [options]` (coarse mesh, full
model, defaults of commit 85866ae; ~85 s per option set on 8 threads, both the measured
backscatter and the uniform cut), and for options the example does not expose a
scratch copy that passes `BubbleAbsorption`'s `delay` and `LipClosure`'s `depth`
(nothing else changed). The coarse mesh reproduces the standard one within ~0.01 ns in
the crossing and 0.006 in P2 ("The P2 swing: capsule radius, backscatter and the
collision"); no standard-mesh run was made. Confidence: **H** primary text, **L** ours.

### 1. The current state, remeasured

| Metric (coarse, full model) | measured backscatter (default) | uniform 48 kJ cut (the setting matching Fig. 10) | PRE 109 Fig. 10 |
|---|---|---|---|
| minimum | −0.181 at 5.25 ns | −0.195 at 5.25 ns | −0.170 at 5.60 ns |
| upward zero crossing | 6.19 ns | 6.24 ns | 6.67 ns |
| P2/P0 at 8 ns | +0.080 | +0.080 | +0.095 |
| rms against Fig. 10, 2–8.5 ns | 0.047 | 0.044 | — |
| capsule / LEH share | 13.9 / 24.2% | 14.0 / 24.1% | HYDRA 16.2 / 25.2% |
| Dante 1 at 8 ns, initial / open area | 289.7 / 312.9 eV | 288.6 / 311.7 eV | measured 310.1 |

- **The swing metrics have not moved** since "Capsule absorption, conversion and the
  ablating-wall release (2026-10-25)" §2: the four swing numbers of both runs agree to
  the printed digit. Later commits (the M-band split, the entrance holes at the caps'
  azimuthal count, Dante 1 as the whole-field view) changed the energetics read-outs
  only: capsule share 13.8 → 13.9%, LEH share 24.7 → 24.2%, Dante 1 at 8 ns 290.9 /
  314.2 → 289.7 / 312.9 eV.
- **Rise and foot** (P2/P0 minus Fig. 10, measured backscatter; the uniform cut is within
  0.004):

  | t (ns) | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 | 3.5 | 4.0 | 4.5 | 5.0 |
  |---|---|---|---|---|---|---|---|---|---|
  | model | +0.008 | +0.015 | −0.017 | −0.003 | +0.052 | +0.046 | +0.009 | −0.030 | −0.174 |
  | model − Fig. 10 | −0.007 | +0.003 | −0.016 | −0.005 | +0.054 | +0.015 | +0.027 | +0.063 | −0.075 |

  The roadmap's "0.06 above and ~0.03 below Fig. 10" (the 2026-10-24 baseline: 4.5 ns
  +0.058, 2.0 ns −0.030) is now: the rise +0.063 above at 4.5 ns and −0.075 below at
  5.0 ns (the model's drop into the minimum comes later and steeper than HYDRA's: Fig. 10
  falls −0.018 → −0.092 → −0.099 over 4.0–5.0 ns, the model +0.009 → −0.030 → −0.174),
  mean −0.006 over 4.25–5.0 ns; the foot within ±0.016 over 1.0–2.5 ns and *above* Fig. 10
  over 3.0–4.0 ns (+0.054, +0.015, +0.027), where the outer power steps from 30 to 60–80
  TW. The change from the 2026-10-24 baseline most likely came with the energetics of
  2026-10-24/25 (wall history, DU wall, conversion, measured pulse), not with later
  commits: the 2026-10-25 runs match today's swing metrics to the digit, but they did not
  record the foot (*L*). The foot's causes stay those named then (mock pointing, DN24's foot
  cone fraction read by eye).

### 2. The computed lips against their own calibration band

**What was fitted to what.** `LipClosure`'s speed (0.022 mm/ns at 100 eV), depth
(0.240 mm) and hold time (0.40 ns) are fitted to Chen et al., PoP 27, 022702 (2020)
Fig. 27: 31 G-LEH effective radii on five basecamp shots, where the edge is the peak
gradient of the image in the G-LEH band, 3–5 keV (Ti filter on N171010 and N171212) or > 5 keV
(Al on N180102, N180103 and N180117) (Chen Table III, **H**); rms 0.033 mm ("Computed LEH
closure"). The example's default prescribed lips (0.857 r₀ along Young's ramp) end at
Young 2020's SXI value, the SXI's clear aperture being measured in its 3–5 keV channel
(Schneider, RSI 81, 10E538 (2010), and HTPD 2012, **H**; *L* that Young's value is that
channel's). The x-ray loss through the hole is set by the thermal band (~0.9 keV at 300 eV).

**What the literature gives for the thermal band.**
- Chen 2020 §II.2 (**H**): the apparent LEH radius grows with photon energy (the 3ω
  critical-density edge, then the ~0.9 keV thermal peak, then the few-keV tail G-LEH
  sees) because the mean free path grows. No thermal-band radius is published.
- Chen 2020 Fig. 7 (**H**, a simulated image of the 0.15 mg/cc basecamp): lowering the
  opacity of the warm (< 0.5 keV), sub-solid (< 1 g/cm³) lip blow-off a hundredfold
  moves the inferred 3–5 keV edge out to the solid LEH radius, i.e. by the whole closure
  depth (~0.24 mm on that platform at peak, *L*). The same paragraph calls the edge only
  weakly sensitive to the opacity and set by the ridge's hydrodynamics. Read as an e-folding length of the edge in opacity, Fig. 7 gives
  ≥ 0.24/ln 100 = 0.052 mm per e-fold (*L*, a lower bound since the edge cannot move past
  the solid rim).
- Opacity ratio between the bands (*L*, cold gold as a proxy for the warm blow-off):
  the NIST X-ray mass attenuation table for Z = 79 (Hubbell and Seltzer) gives μ/ρ
  4652 cm²/g at 1 keV against 2049 / 1144 / 666 at 3 / 4 / 5 keV
  and 425 at 6 keV, a ratio of 2.3–7 for the Ti band and ~7–11 above 5 keV.
- Schneider et al., HTPD 2012 (OSTI 1043636, **H**): the 870 eV SXI channel sees part of
  the thermal flux from outside the 3–5 keV clear aperture; the Dante correction
  f (R_LEH/R_CA)² is ~1.15 whatever the LEH size. With R_CA 0.84–0.86 r₀, f ≈ 0.81–0.85
  (*L*): the lip blow-off emits in the thermal band, and Dante's open-area value
  brackets that emission; it does not bound a smaller hole.

**The thermal-band radius the model implies** (*L*): the edge at the G-LEH value minus
Δ = 0.052 mm × ln(ratio), i.e. Δ ≈ 0.04–0.10 mm for the Ti-band shots and 0.10–0.13 mm
for the Al-band ones, between 0 (Chen's hydrodynamic edge) and ~0.13 mm. On N221204 the
computed lower-LEH edge is 1.513 mm at 2.0 ns, 1.420 at 4.5, 1.320 at 5.5 and stays at
1.310 mm (0.845 r₀) from 8.0 to 9.0 ns (1.305); the thermal-band edge would be
~0.06–0.12 mm further in, ~1.19–1.25 mm through the peak. Since the fit mixes the two
bands (2 Ti, 3 Al shots), the implied depth for the x-ray loss is ~0.30–0.36 mm instead
of 0.24.

**What it does** (`lips=computed`, `LipClosure(depth = …)`, the other lip parameters
fitted; coarse; measured backscatter, the uniform cut in parentheses):

| depth (mm) | edge at 5.5 / 8.0 ns (mm) | minimum (at) | upward zero | 8 ns | rms vs Fig. 10 | capsule Tr 8.0 / 8.5 ns (eV) | Dante 1 at 8 ns, initial / open area | LEH / capsule share |
|---|---|---|---|---|---|---|---|---|
| 0.24 (fitted, `lips=computed`) | 1.320 / 1.310 | −0.175 (−0.189) at 5.25 | 6.18 (6.23) | +0.081 (+0.080) | 0.047 (0.043) | 307.1 / 299.3 | 288.2 / 313.5 | 23.5 / 14.0% |
| 0.30 (thermal, Δ 0.06) | 1.265 / 1.250 | −0.170 (−0.184) at 5.25 | 6.17 (6.22) | +0.082 (+0.082) | 0.047 (0.044) | 309.4 / 301.7 | 282.8 / 314.9 | 22.1 / 14.2% |
| 0.36 (thermal, Δ 0.12) | 1.210 / 1.190 | −0.165 (−0.178) at 5.25 | 6.16 (6.21) | +0.084 (+0.083) | 0.048 (0.044) | 311.7 / 304.1 | 277.2 / 316.4 | 20.7 / 14.5% |
| prescribed lips (default) | — | −0.181 (−0.195) at 5.25 | 6.19 (6.24) | +0.080 (+0.080) | 0.047 (0.044) | 306.3 / 298.4 | 289.7 / 312.9 | 24.2 / 13.9% |
| HYDRA / measured | | −0.170 at 5.60 | 6.67 | +0.095 | — | | measured 310.1 | 25.2 / 16.2% |

- **The swing barely moves.** Per 0.06 mm of thermal offset the crossing moves −0.01 ns,
  the minimum +0.005 shallower and the 8 ns value +0.001 to +0.002: a smaller hole adds
  pole drive a little through the whole peak. The thermal band cannot carry any of the
  0.45 ns crossing gap; it moves the crossing the wrong way.
- **The energetics move.** The LEH share falls 1.4 points per 0.06 mm (23.5 → 22.1 →
  20.7%, HYDRA 25.2%), the capsule share rises 0.2–0.3 points (14.0 → 14.5%, HYDRA
  16.2%), capsule Tr rises 2.3 eV per 0.06 mm at 8 ns. Dante 1 over the initial area
  falls 5–6 eV per 0.06 mm (measured 310.1 eV at 8 ns), over the open area it rises
  1.4 eV. Only the hold depth is offset here: the edge in the foot (free expansion,
  1.513 mm at 2 ns) is unchanged, so the foot's swing is unchanged by construction.
- So the thermal band trades the LEH share against the capsule share, in opposite
  directions relative to HYDRA, and it widens the gap between Dante's initial-area
  reading and the measurement. Neither a thermal offset nor the fitted G-LEH edge is
  preferred by these checks; the band difference stays an *L* bracket, and the fitted
  depth is unchanged.

### 3. Published bubble inputs not yet used

**HYBRID-E bubble speed and collision time: still no measurement.** Searches (2026-11-18:
HYBRID-E, 6.4 mm Au-lined DU, N210808, N221204 with bubble, G-LEH, thin-wall, inner-spot
motion) found no HYBRID-E bubble trajectory or spot separation; the inputs are those of
"The P2 swing" (PRE 109 Fig. 1a's tips; the picket fluence; the collision near 6.1 ns in
Fig. 1a). One more basecamp datum: Ralph et al., PoP 28, 072704 (2021) (OSTI 1860800,
**H**) measured the Au bubble on the hohlraum-science platform with the 0.908 mm HDC
capsule (0.3 mg/cc, N190325 and repeats; the basecamp's size, *L*) at 364 km/s averaged from t = 0, and 399 km/s with a Ta₂O₅ liner (~10%
faster). The calibrated model's average from t = 0 on that geometry is 0.32–0.33 mm/ns at
5–6 ns (0.36 mm/ns from 0.55 ns), so Ralph's number points to a ~10% faster bubble
(*L*: the frame times are not in the text); it is the calibration platform, not
HYBRID-E, and the wrong direction for Fig. 10.

**Pak's picket** (`bubble_energy=1.23`: N221204's delivered pulse split by DN24's cone
fraction gives 43.0 kJ of outer power in the first ns after the turn-on, against the
basecamp's ~26 kJ, × 0.73 for split quads; "The P2 swing", not run then): t_c 5.99 ns,
tip 2.01 mm at 6 ns. Crossing 6.10 ns (−0.09; uniform cut 6.17, −0.07), minimum −0.171
(uniform −0.184) at 5.25 ns, rms 0.053 / 0.049: earlier, as expected.

**Izumi 2023's stall and relaunch** (Izumi et al., PoP 30 (2023), doi:10.1063/5.0136145,
OSTI 1971298, **H**): in a low-trough 4-shock CH design (1 TW trough after a high picket,
280 TW peak) the picket-launched tip stalls in the trough and is relaunched at ~6 ns, at
most 300 µm/ns, slowing near the axis; extrapolated ballistically it appears launched
~7 ns late. The trajectory is ballistic for foot-to-peak intensity contrasts below ~15×,
and the 3-shock HDC basecamp (the calibration platform) is in that class. N221204's
delivered pulse (`n221204_pulse()`, split by DN24's cone fraction): trough 28.3 TW total
(24.5 TW outer) at 1.90 ns, below 40 TW for ~1 ns (1.75–2.75 ns), against 433–448 TW
(outer 290–300 TW) in the peak: contrast 16× total, 12× outer, at the criterion's edge,
and 2× below the picket of that shot (55 TW). By the published criterion the bubble is
ballistic; a stall is not excluded.

*No new keyword is needed.* A stall from t_s to t_r with relaunch at the same speed
leaves the tip at v_b (t − t_L − (t_r − t_s)) from t_r on, which is `delay` + (t_r − t_s);
before t_r the two differ only while the tip is shallower than its stall depth (0.40 mm
at t_s = 1.75 ns), where the inner beams lose ≤ 3.3e-9 of their power to the bubble
(`bubble_fractions` on N221204's geometry: the inner beams start to be clipped at
1e-3 at 1.05 mm, 3.5 ns). A different relaunch speed is `speed` with the matching
`delay`. So the existing keywords carry a stall exactly where it matters, and
`src/bubble.jl` is unchanged. The collision time is kept (the waist's ablator–wall
collision does not depend on the outer bubble).

| Run (coarse; measured backscatter / uniform cut) | tip at 6 ns | minimum (at) | upward zero | 8 ns | rms vs Fig. 10 | capsule share |
|---|---|---|---|---|---|---|
| default (delay 0.25 ns) | 1.96 mm | −0.181 (5.25) / −0.195 (5.25) | 6.19 / 6.24 | +0.080 / +0.080 | 0.047 / 0.044 | 13.9 / 14.0% |
| half speed through the trough (delay +0.5 ns) | 1.78 mm | −0.217 (5.25) / −0.233 (5.50) | 6.58 / 6.73 | +0.071 / +0.070 | 0.043 / 0.046 | 14.1 / 14.3% |
| full stall through the trough (delay +1.0 ns, a bound) | 1.60 mm | −0.235 (5.25) / −0.260 (5.50) | 7.22 / 7.27 | +0.048 / +0.048 | 0.055 / 0.061 | 14.3 / 14.5% |
| Pak's picket (`bubble_energy=1.23`) | 2.01 mm | −0.171 (5.25) / −0.184 (5.25) | 6.10 / 6.17 | +0.082 / +0.081 | 0.053 / 0.049 | 13.8 / 13.9% |
| PRE 109 Fig. 10 / Fig. 1a | N210808 1.85 ± 0.1; N221204 ≥ 1.55 | −0.170 (5.60) | 6.67 | +0.095 | — | HYDRA 16.2% |

- **The crossing moves by about the delay of the tip** (+0.39 / +0.49 ns for +0.5 ns,
  +1.03 ns for +1.0 ns) with t_c held at 6.09 ns: the ridge stop alone does not turn the
  swing positive; the bubble's absorption does. A delay of ~0.5 ns, i.e. the tip at half
  speed through N221204's ~1 ns trough, puts the uniform-cut crossing at 6.73 ns against
  Fig. 10's 6.67, with the tip at 1.78 mm at 6 ns (inside N210808's 1.85 ± 0.1 of
  Fig. 1a and N221204's ≥ 1.55).
- **It deepens the minimum** by 0.036–0.038 (−0.233 against Fig. 10's −0.170): more inner
  light reaches the waist through 5–6 ns. The minimum's time moves to 5.50 ns (uniform
  cut; Fig. 10 5.6). Timing and depth do not close together; the depth's other levers
  (the backscatter history, −0.029; the lag, "Computed LEH closure") are separate.
- Capsule Tr rises by ≤ 2.0 eV (+0.5 ns) and ≤ 4.0 eV (+1.0 ns) over 5.5–7.5 ns, Dante 1
  by ≤ 1.3 eV at 8 ns: a later bubble moves less inner power to the outer footprint near
  the LEH.
- **Not adopted** (no default changed): a 0.5 ns delay is the gap Fig. 10 leaves, and
  choosing it would be fitting to the blind check; Izumi 2023's criterion puts N221204
  in the ballistic class, so no published datum supports a stall. Recorded as the size
  of what an unmodeled slow-down in the trough would have to be.

### 4. What drives bubble growth (for item 36)

From `reference/litreview/bubble_drivers.md` (Callahan 2018, Ralph 2018, Izumi 2021,
Izumi 2023, Ralph 2021; **H** unless noted):
- **Measure.** Every published scaling uses the outer cones' *picket fluence* (picket
  energy over the FWHM, net of outer backscatter, per phase-plate spot area, i.e. the
  outer irradiance integrated over the picket, Ralph 2018), not a peak irradiance.
  Izumi 2021's split-quad shot (27% lower wall fluence at the same beam energy, tip
  −120 µm) shows the fluence *on the wall* (spot overlap) counts (*L* for re-pointing:
  no published re-pointing study).
- **Epoch.** The picket, for pulses below ~15× foot-to-peak contrast (the tip then moves
  at constant speed through the rise and the peak, extrapolating back to the outer
  turn-on); the main drive as well in low-trough pulses (Izumi 2023's relaunch).
- **Scaling.** Speed ∝ fluence^0.2 ± 0.03 (Izumi 2021, seven shots); ∝ I^0.36 (Ralph
  2018's model; data fit I^1/3 to I^2/3 with I measured as the picket fluence); bubble
  size ∝ (E/A)^0.5 in Callahan 2018's 39-shot metric (R² 84%; 1/3 and 2/3 give 76 and
  79%). Laser mass ablation and exhaust speed ∝ I^1/3 (Izumi 2023); 10% less ablated
  mass extends inner-beam propagation by ~0.4 ns (Izumi 2023).
- **Against `src/bubble.jl`.** The model's only intensity dependence is the scalar
  `picket_ratio` E (speed ∝ E^0.2, Izumi 2021; collision ∝ E^-0.16, Izumi 2018), set by
  hand; it does not read the computed wall irradiance. It does follow the pointing: the
  bubble plane sits where the 50° axes reach the (moving) wall and the clipping is
  computed from the actual inner beams. On N221204, E 1.1 → 1.23 moves the crossing by
  −0.09 ns (§3).
- **For item 36** (*L*): a ceiling on the peak irradiance over wall elements during the
  peak is not the quantity the bubble data scale with. The supported measure is the
  outer picket fluence on the wall under the outer footprint, with a weak exponent
  (0.2–0.5): a 20% rise in local picket fluence is a ~4% faster tip, ~−0.15 ns in the
  crossing with the collision following E as in the model (scaled from the E = 1.23
  run, −0.09 ns for +12%). Two options, the user's: put
  `BubbleAbsorption` in the study's laser–plasma chain (it already moves the bubble's
  effect with the pointing) and add a picket-fluence ceiling on the outer footprint, or
  keep a peak-irradiance ceiling as an engineering limit labeled as such (LPI and
  damage thresholds scale with peak intensity; the bubble does not on the published
  data).

### 5. The bubble's speed from the pointing's picket fluence (item 10)

`picket_fluence`, `picket_ratio`, `basecamp_picket_fluence` (`src/picket_fluence.jl`,
`test/picket_fluence.jl`) and `remake(bubble; picket_ratio)` close the link of §4 as an
opt-in: E, `BubbleAbsorption`'s `picket_ratio`, from the outer picket fluence of a
pointing on the wall. *Derived*: a link built on the published definitions and exponents,
not a facility criterion; no default, exponent or calibrated constant changed, and no
new constant was fitted.

**Definition.** F is the deposition of the outer beams (44.5° and 50°) at their picket
energies (pulse integrated over a window, here the first ns after the outer turn-on)
divided by each wall element's area (kJ/mm²; ×1e5 → J/cm²). The summary the ratio reads
is the energy-weighted mean F̄ = ∫F² dA / ∫F dA = E/A_eff, with A_eff = (∫F dA)²/∫F² dA
the footprint's effective area. Which papers' definition (**H** for the definitions,
`bubble_drivers.md` §2; **L** that F̄ is their computed counterpart):

| Source | Fluence as published | F̄ reproduces |
|---|---|---|
| Ralph 2018 (speed ∝ I^0.36) | outer picket energy (FWHM, net of backscatter) over the outer phase-plate spot area | one isolated spot: its energy over its effective area (flat top: E/A exactly) |
| Callahan 2018 (size ∝ (E/A)^0.5) | outer picket energy over the one-side outer spot area, 53 mm² (all of one side's outer beams) | a footprint of many spots: energy over the footprint's effective area |
| Izumi 2021 (speed ∝ E^0.2) | picket fluence relative to the basecamp's; −27% from quad and cone splitting at the same energy | the overlap of spots on the wall: k coincident spots k× one, separated spots 1× |

The peak over elements was not chosen (it is one element and mesh-dependent; the dz −0.5 mm
variant below has a peak ratio 1.126 at h 0.4 mm and 1.143 at h 0.2 mm, where F̄'s ratio
moves 1.1629 → 1.1499), and the mean over the area above half the peak (inherits the peak's
threshold, and is not smooth in the pointing). F̄ is an integral: smooth in the aims
(ForwardDiff through the pointing, checked against central differences to 1e-5), and its
mesh error is O(h²).

**Checks** (`test/picket_fluence.jl`; flat plate, normal incidence, collimated beams,
elements 0.04 mm, 0.06 mm for the pair):
- one spot of energy 3 kJ: F̄ = 2^(−2/n) F₀ with F₀ = E/(π a b Γ(1 + 2/n)): Gaussian
  (a = 0.5 mm) −8.9e-4 relative, the element-average estimate h²/(6a²) being −1.1e-3;
  super-Gaussian (0.6 × 0.4 mm, n 4.5) −1.7e-3; both converge as h² (n = 50/100/200 cells:
  −3.6e-3/−8.9e-4/−2.2e-4 for the Gaussian);
- two Gaussian spots a distance d apart: F̄ = (1 + e^(−d²/2a²)) × one spot's: 2 exactly
  when they coincide, 1.6071 against 1.6065 at d = a, 1.1359 against 1.1353 at d = 2a,
  1.000004 at d = 5a;
- a 50° quad (four beams, `nif_beams`) on the basecamp's cylinder (h 0.2 mm), split
  tangentially by ±s per beam pair (`hybride_mock_pointing`'s rule): F̄ falls to 0.728 of
  the unsplit quad's at s = 0.25 mm, 0.549 at 0.5, 0.521 at 0.75 (pairs 1.5 mm apart,
  toward half), with the energy on the wall unchanged to 3e-4. Izumi 2021's split shot
  lost 27% of its wall fluence; a 0.5 mm separation of the pairs (Kritcher 2021: 450 µm
  on HYBRID-E) gives the same −27% on this geometry (**L**: the split geometry of Izumi's
  shot, N180117, is not published, and it split the cones as well).

**Mesh** (F̄ on the hohlraum only, separate from the radiation mesh; kJ/mm²):

| h (mm) | basecamp, 26 kJ | N221204 mock pointing, 43.04 kJ | E = ratio |
|---|---|---|---|
| 0.4 | 0.5040 | 0.5861 | 1.1629 |
| 0.2 | 0.5204 | 0.5985 | 1.1502 |
| 0.1 | 0.5250 | 0.6023 | 1.1472 |
| h² extrapolation | 0.5265 | 0.6035 | 1.1462 |

F̄ is 1.2% (basecamp) and 0.8% (N221204) low at h 0.2 mm, the ratio 0.35% high; the
ratios between pointings on the same mesh (table below) agree between h 0.2 and 0.1 to
0.1–0.7%. `basecamp_picket_fluence` defaults to h 0.2 mm (2–3 s on 8 threads; the
N221204 operator for 192 beams takes 3 s at h 0.2, 10 s at h 0.1). The coarse radiation
mesh of the benchmark (0.8 mm) is too coarse for the map.

**The reference** (`basecamp_picket_fluence`). E = 1 is where the tip speed was fitted
(Izumi 2021's basecamp, "Hohlraum dynamics calibrated to published data"): 5.75 × 10.13 mm
Au hohlraum, 3.37 mm LEHs, 0.908 mm capsule, 128 outer beams unsplit at the package's
proxy pointing for this hohlraum (z 5.20 / 5.01 mm on the axis, the pointing the opacity
calibration used; the facility's is unpublished, **L**), 26 kJ of outer picket in the
first ns after the turn-on (Izumi 2021 Fig. 3a read by eye, "The P2 swing", **L**). The
picket energy is the least certain input of an absolute E: Young 2020's requested pulse
for N170827, another 5.75 mm HDC shot (`young2020_pulses()`), carries 42.2 kJ of outer
power in the first ns after its 0.25 ns turn-on; if the basecamp's picket were that, every
absolute E would be 0.62 of the values below. Ratios between pointings do not depend on
it.

**The N221204 pointing against the hand-set value.** The delivered pulse split by
DN24's cone fraction (`examples/n221204.jl`) gives 43.04 kJ of outer picket in
0.4–1.4 ns (§3's 43.0), 1.655× the basecamp's 26 kJ. On the mock HYBRID-E pointing
(quads split ±0.5 mm) F̄ gives **E = 1.147** (h 0.1 mm) against the hand-set 1.1 (1.5×
energy by eye × 0.73 for split quads) and Pak's 1.23 (1.69 × 0.73): energy 1.655 × the
larger hohlraum and its unsplit footprint 0.826 (E 1.368 with s = 0) × the mock's split
0.839. The hand-set value is 4% below the computed one, i.e. +0.8% in tip speed.

**Re-pointing the outer cones** (h 0.1 mm; F̄ relative to the mock pointing; energy on
the wall; then the full model, coarse mesh 0.8 / 0.4, measured backscatter, E fixed at
1.1 or following the pointing, `picket_ratio(design, mock; at = 1.1)`, through
`remake(bubble; picket_ratio)` so the collision follows E^-0.16):

| Outer cones (both, mirrored) | wall picket energy (kJ) | F̄ / mock | E linked | E 1.1: min (at) / upward zero / 8 ns | E linked: min / upward zero / 8 ns | crossing, linked − 1.1 |
|---|---|---|---|---|---|---|
| mock (z 5.7 / 5.8, s 0.5) | 43.04 | 1 | 1.1 | −0.181 (5.25) / 6.19 / +0.080 | the same (E = 1.1 exactly) | 0 |
| z −0.5 mm | 43.03 | 0.9998 | 1.0998 | −0.152 (5.00) / 6.14 / +0.063 | identical to 0.0001 | +0.0002 ns |
| z −0.3 mm | 43.04 | 1.0000 | 1.1000 | −0.161 (5.25) / 6.15 / +0.074 | identical | 0.0000 ns |
| z +0.3 mm | 42.91 | 1.0015 | 1.1017 | −0.205 (5.25) / 6.25 / +0.076 | −0.205 / 6.25 / +0.076 | −0.001 ns |
| z +0.5 mm (clips 2.8% of the picket on the LEH) | 41.85 | 1.0097 | 1.1107 | −0.227 (5.25) / 6.42 / +0.067 | −0.226 / 6.41 / +0.067 | −0.013 ns |
| unsplit quads (s 0) | 43.04 | 1.1923 | 1.3115 | −0.177 (5.25) / 6.19 / +0.080 | −0.160 / 6.00 / +0.082 | −0.189 ns |
| s 0.225 (pairs 0.45 mm apart, Kritcher 2021) | 43.04 | 1.0587 | 1.1646 | −0.178 (5.25) / 6.19 / +0.080 | −0.173 / 6.15 / +0.081 | −0.041 ns |
| s 0.75 | 42.78 | 1.0869 | 1.1955 | −0.187 (5.25) / 6.21 / +0.080 | −0.180 / 6.15 / +0.082 | −0.057 ns |
| s 1.0 (not run: clips 8% on the LEH) | 39.58 | 1.1784 | (1.296) | | | |
| mock, absolute E (`basecamp_picket_fluence`, h 0.1) | 43.04 | | 1.1472 | | −0.177 (5.25) / 6.17 / +0.081 | −0.028 ns |
| PRE 109 Fig. 10 | | | | −0.170 (5.60) / 6.67 / +0.095 | | |

- **Moving both outer cones along z barely changes the picket fluence** (F̄ within
  −0.02% to +1% over ±0.5 mm: the two cones move together, so their overlap and the
  quads' splits are unchanged and the incidence changes little), and the bubble link
  moves the crossing by ≤ 0.013 ns. The same re-pointing moves the swing through the
  geometry by −0.05 to +0.23 ns (crossing 6.14–6.42 ns, minimum −0.15 to −0.23): what
  a z re-pointing does to the swing is the deposition's, not the bubble speed's.
- **Overlap is what moves it.** Unsplitting the outer quads raises F̄ 19% and moves the
  crossing −0.19 ns through the bubble only (−0.01 ns through the geometry only); a 0.45 mm
  split +6% (−0.04 ns); s 0.75 +9% (−0.06 ns: the split pairs begin to land on the
  neighboring quads' spots, so F̄ rises again past s ≈ 0.5, and at s 1.0 the beams clip
  8% on the entrance hole). Per unit of ln E the crossing moves −0.7 to −1.1 ns
  (−0.12 to −0.20 ns for +20% in F̄, against §4's −0.15 ns estimate); the minimum rises
  by 0.005–0.017 and the 8 ns value by ≤ 0.002.
- **The peak irradiance would have told a different story.** Over the same variants the
  peak element fluence moves ≤ 0.02% with z (h 0.1 mm) but +39% without the split
  (1.30 against 0.94 kJ/mm²), +19% at s 0.225 and +20% at s 0.75, i.e. 2–3× the
  mean's relative change: a ceiling on the peak follows the overlap in direction, not in
  size.
- **The computed absolute E** (1.147) moves the crossing −0.03 ns from the hand-set
  1.1; Fig. 10's crossing stays 0.5 ns later. The link cannot close the gap (it moves
  E up, not down, on this pointing), and nothing here was fitted to Fig. 10.

Confidence **L**: the exponents (E^0.2, E^-0.16) are published for pickets and splits on
one platform; applying them to a fluence changed by re-pointing at fixed energy has no
published test (Izumi 2021's split shot is the nearest, and it also split the cones);
the pointing is the mock HYBRID-E rule; the absolute E rests on a proxy basecamp pointing
and a picket energy read by eye (above). The relative numbers (F̄ between pointings on
one mesh) are **derived**, converged to < 1%. Not adopted: the example keeps the
hand-set 1.1 (`bubble_energy`); adopting the link, absolute or relative, is the user's
decision.

**The picket ceiling, the linked bubble and the prediction** (later on 2026-11-18;
`test/picket_fluence.jl`). The user's decision on item 36's quantity: the re-pointing
study's ceiling holds the *maximum laser irradiance at the picket* (not at peak power),
and the model chain should predict what a re-pointing does to the bubble. Three opt-in
pieces, *derived* (no default, exponent or calibrated constant changed; nothing fitted):
- `picket_irradiance(ev, names...; pulses, window | power, groups, measure)` and
  `picket_ceiling(reference, names...; ...)`: the irradiance during the picket, each
  beam at its mean power over the picket window (or a given picket power), the outer
  cones by default, deposited by the scenario's operator with the setup beams and
  without the scenario's laser–plasma model (set for the time of that model). Its ratios are
  the picket fluence's. `measure = :peak` (the user's quantity; the 64-norm stand-in with
  the exact peak held, `exact = true` by default) or `:weighted` (F̄, what the bubble reads; named `:mean` until 2026-11-18).
- `PicketLinkedBubble(bubble, target, reference; pulses, window)`: `BubbleAbsorption`
  with E = `bubble.picket_ratio` × F̄(beams)/F̄(reference) recomputed from the beams of
  every evaluation (two picket depositions per evaluation), so that `evaluate`, mode
  gates, objectives and `drive_history` see the bubble a re-pointing grows; at the
  reference it is exactly `bubble` (`picket_ratio` now forms `at × (F̄/F̄_ref)`, `at`
  exactly at the reference: a change of rounding only). The per-beam sensitivities
  (`mode_sensitivity`) deposit beams one at a time and refuse it (`linked_bubble`
  resolves it for given beams).
- `bubble_prediction(model, design, reference; history, fixed)`: the peak and mean
  picket ratios, E, the tip-speed ratio, the collision times and, with a drive-history
  function, P2/P0 of both pointings with the swing's minimum and upward zero crossing,
  split with `fixed = true` into the change through the deposition (E held) and through
  the bubble.

Checks: the values against `picket_fluence` to 1e-12; gradients of both ceilings and of
an evaluation through the linked bubble against central differences to 1e-5; one-
parameter constrained runs (the 50° cone's z pulled toward the waist on a 0.45 mm
hohlraum) end at the bisection root of the held quantity to 1e-5 mm with `:weighted` and with
the peak's stand-in (`exact = false`); with the exact peak the run ends 0.4% below the
limit, the stand-in running above the peak near that limit and the working limit only
being lowered (the same machinery as `irradiance_ceiling(...; exact = true)`). The linked
bubble equals `remake(bubble; picket_ratio)` by hand to 1e-12, and a drive history with
it equals that of the resolved bubble exactly.

*Inference (2026-11-18, later).* `BubbleAbsorption`'s `apply_lpi`, `lpi_deposit` and
`bubble_fractions` formed their number type by splatting over the beams in `promote_type`
and inferred abstract containers; they now dispatch on the beams' type
(`AbstractVector{Beam{S}}`; a vector of mixed types is converted to one type first, which
also fixes that conversion for beams whose aim and power differ in type). They infer
concretely on Float64 and on the optimizer's duals, for the bubble and the linked bubble
(`test/inference.jl` runs `@inferred`, and JET from Julia 1.11; the clipping is also checked
on 1.10). Results are unchanged for an evaluation with the bubble at 6 ns (Float64, and
12-partial duals' values and partials) and the N221204 coarse CSV, against the base commit.
What it bought is the guard, not time: one `evaluate` at 6 ns on the inference fixture's 48
quads (one thread, best of 10) 23.6 → 23.1 ms (22.7 ms without the bubble), with duals 44.4
→ 44.2 ms; `apply_lpi` alone 0.57 → 0.56 ms, its allocations 52 → 46 kB (Float64) and 123 →
101 kB (duals). The per-beam iris integrals were already behind a function barrier.

**The benchmark's table, linked column, through `bubble_prediction`** (coarse mesh, full
model, measured backscatter; F̄ and the peak on the h 0.1 map; the reference is the mock
pointing at E = 1.1, minimum −0.181 at 5.25 ns, upward zero 6.194 ns):

| Outer cones | F̄ / mock | peak / mock | E | tip speed | min (at) / upward zero / 8 ns | shift: min / zero |
|---|---|---|---|---|---|---|
| z −0.5 mm | 0.9998 | 0.9999 | 1.0998 | ×1.0000 | −0.152 (5.00) / 6.143 / +0.063 | +0.029 / −0.051 ns |
| z −0.3 mm | 1.0000 | 0.9999 | 1.1000 | ×1.0000 | −0.161 (5.25) / 6.154 / +0.074 | +0.020 / −0.040 ns |
| z +0.3 mm | 1.0015 | 1.0000 | 1.1017 | ×1.0003 | −0.205 (5.25) / 6.248 / +0.076 | −0.024 / +0.054 ns |
| z +0.5 mm | 1.0097 | 1.0000 | 1.1107 | ×1.0019 | −0.226 (5.25) / 6.406 / +0.067 | −0.045 / +0.212 ns |
| unsplit (s 0) | 1.1922 | 1.3913 | 1.3115 | ×1.0358 | −0.160 (5.25) / 5.996 / +0.082 | +0.021 / −0.198 ns |
| s 0.225 | 1.0587 | 1.1879 | 1.1646 | ×1.0115 | −0.173 (5.25) / 6.147 / +0.081 | +0.008 / −0.047 ns |
| s 0.75 | 1.0869 | 1.1952 | 1.1955 | ×1.0168 | −0.180 (5.25) / 6.149 / +0.081 | +0.001 / −0.045 ns |

- **Reproduced**: every E, minimum, crossing and 8 ns value of the linked column above
  (the earlier sweep set E by hand from the same map), to the digits printed there and
  to 1e-5 in the saved values; the s 0.75 row's 8 ns value is +0.0815, i.e. +0.081 (the
  table above prints +0.082). The shift columns include the deposition's change, which
  for a z re-pointing is nearly all of it (§5 above).
- **The peak does not bound the mean.** Moving the outer cones up by 0.5 mm leaves the
  picket's peak element unchanged (1.0000) and raises F̄ by 1% (with 2.8% of the
  picket clipped on the entrance hole; the causes were not separated); splitting changes the peak 2–3 times as much as
  F̄. A ceiling at the reference's peak holds the overlap of the spots, the main lever
  on F̄ here, but a design can grow its bubble by ~1% in F̄ (0.2% in tip speed) under
  it. The bubble's prediction reads F̄; the peak is the user's limit.

**The dropped-bundle tutorial's design, given a pulse** (scratch; the tutorial's `limited`
run reproduced on the mesh and options of the tutorial: 997 evaluations, ending at the box,
under its peak-power ceiling with `exact = false`). Its scenario has no pulse; Young et al.
2020's requested N170827 pulse (`young2020_pulses()`, a 5.75 mm HDC hohlraum like the
tutorial's) per quad, the picket the first ns after the outer turn-on at 0.25 ns,
`BubbleAbsorption(radius = 2.875, turn_on = 0.25, capsule_radius = 0.9)` at E = 1 for the 48
quads. On the scenario's mesh, freeing bundle 11 leaves the picket's peak (×1.0000) and
lowers F̄ to ×0.9865 (its outer quads are gone); the re-pointed design raises the picket's
peak to ×1.0110 (its peak-power irradiance ×1.0017) and F̄ to ×0.9939. On a 0.2 mm map the
design's F̄ is ×0.9934 and its peak ×0.9983 of the 48 quads': the peak's ratio depends on
the mesh by more than the change in the ratio itself. The predicted tip speed is ×0.9987,
and the collision +0.006 ns. A drive history on the tutorial's static model (constant albedo
0.85, no wall motion, 1–8 ns): the upward zero of P2/P0 moves +0.153 ns (6.176 → 6.328 ns),
of which +0.144 ns through the deposition and +0.009 ns through the bubble; the minimum
(−0.136 at 4.25 ns) moves by −0.001. The static model is a proxy for the swing's timing (no
moving wall, no closing holes, no CBET); the split between the two causes is the result.

Confidence **L** for the link itself (above); the picket powers of a design without a
separate pulse are an assumption to state with each use.

### What remains

- A HYBRID-E bubble-tip or inner-spot history (unchanged: the one datum that separates
  the picket, a trough slow-down and the collision).
- A published thermal-band (~0.9 keV) LEH radius, or the warm lip's opacity, to replace
  the 0–0.13 mm bracket of §2.
- The foot's +0.03–0.05 over 3–4 ns (mock pointing, DN24's foot cone fraction).
- §5's link: whether the example reads E from the pointing (a decision), and a published
  bubble speed at fixed picket energy with a different spot overlap or phase plate (the
  datum that would test it).

## Backed windows: what delays the emergent flux, and eq. 3 (2026-11-18) *derived*

ROADMAP item 28. `HeatFront` puts a 2.5 µm gold window through at 4.6–4.9 ns on an
ignition-scale drive; the window-emergent flux that Lester et al. (Phys. Rev. E 114,
L023201 (2026)) extracted from Dante 2 has its leading edge at 7.5 ns. Sources, what was
read and at what confidence: `reference/litreview/window_burnthrough.md`. Calculations:
`examples/tamped_foil.jl` (`nova`, `table [scale]`) and the closed-form testset "window
modes against the closed form" in `test/windows.jl`. Nothing below is fitted to the
7.5 ns edge except the two calibrations labeled as such; the backing's burn rates are
fitted to Olson et al.'s burnthrough measurements (SAND2007-7414C), an independent data
set. Confidence: **H** primary text, **H-img** read from a figure, **L** derived here.

### What the 7.5 ns measures **H**
- The windows are 2.5–3.6 µm of gold (2.5 µm on the two-window ignition shot in the
  Letter; Haines et al., Phys. Plasmas 33, 112712 (2026), Sec. VII, give "≈ 3 µm" for
  the same shot) backed by 35 µm of EPON at 1.185 g/cm³ (4.15 mg/cm², about the gold's
  4.83 mg/cm² at 2.5 µm). The backing was chosen to hold the membrane in place.
- The leading edge is the first emergence, seen from outside through gold and backing,
  of a background-subtracted difference signal; it is near zero to about 7 ns.

### The account
1. **Holding the gold dense does not delay its front.** In mass coordinates the diffusion
   flux depends on density only through κ(ρ)e(ρ) ∝ ρ^0.06 (Hammer & Rosen's fits), so a
   static front heats m ∝ ρ^−0.03 (*L*). The 1-D foil of `examples/tamped_foil.jl` (Hammer &
   Rosen gold, flux-limited diffusion, the backing as a radiatively inert ideal gas of its
   areal mass) on Nova's drive: transmission onset (1%) of the backed foil 0.03–0.07 ns
   *earlier* than the bare one, its rear gold at 0.8–1.1 g/cm³ against 0.30–0.35 bare (the
   backing holds it near the 1 g/cm³ the Letter assumes). Against Nova's measured
   half-signal times the 1-D onset is 0.04–0.23 ns late (bare), `HeatFront` 0.09–0.23 ns
   early. The result is converged: doubling the zones or halving the step moves the onset
   and the emergence by ≤ 0.01 ns; the backing's γ (1.4–3) by ≤ 0.005 ns.
2. **The leading edge is the backing's burnthrough.** Cold epoxy of 4.15 mg/cm² is opaque
   below ~2 keV (τ ≈ 10 at 1 keV). The flux the gold transmits burns it through at the
   measured low-Z rate ṁ = A (F/σ(200 eV)⁴)^(3/4); A = 4.8 mg/cm²/ns for an epoxy
   (Olson's polyimide, the closest CHNO polymer measured), 5.9 for Be, ×0.7–1.3 for the
   scatter and the material (**H-img**, fits **L**). After the gold's 1-D onset the burn
   takes 1.8–2.4 ns for 35 µm of epoxy and 2.7–3.3 ns for 50 µm of Be (1.5–3 µm of gold,
   table below): it is most of the delay.
3. **Prediction for the published case** (1-D gold plus the backing at Olson's rate, no
   parameter set from the 7.5 ns): on the public drive below, 2.5 µm with 35 µm of epoxy
   emerges at 7.00 ns (6.68–7.51 for A ×1.3–×0.7); with the drive ×0.95 (the waist
   wall's incident temperature a few percent below Dante's, the scale at which the bare
   front lands at the 4.8 ns the roadmap quotes) at 7.41 ns (7.06–7.97). For ≈ 3 µm
   (Haines 2026) 7.41 and 7.86 ns. Measured 7.5: the account reproduces the edge within
   its uncertainty, which is dominated by A (±0.4 ns), the drive (±0.4 ns) and the
   thickness (2.5 or 3 µm, +0.4 ns).
4. **The front-plus-backing estimate from `BurnthroughAlbedo` is about 1 ns early**
   (6.24 and 6.60 ns): its transition gives the gold its full gray transmission 0.4 ns
   after the front. In 1-D the transmission starts 0.6 ns after `HeatFront`'s arrival on
   this drive (0.2–0.45 ns on Nova's) and rises slowly: for 2.5 µm 0.04 at +0.65 ns,
   0.08 at +0.9, 0.11 at +1.4, 0.12–0.16 to the end of the peak, against 0.17–0.22 for
   `RosselandTransmission`; the energy transmitted 4–9 ns is 0.62 of the transition's
   with `RosselandTransmission` (2.0 and 2.5 µm; 0.49–0.54 of the constant default's,
   1.5–3 µm: "A diffusion transmission for burned-through windows (2026-11-18)"). The Letter's extracted window flux, about 1.3 TW/sr near 8.5–9 ns
   (**H-img**), is 0.18 of a 2 mm Lambertian disk at Dante's ~305 eV seen from Dante 2's
   angle (*L*; the Letter notes its solid-angle calibration is pending), between the two.
5. **For the hohlraum the backing changes nothing to first order.** Until the backing
   burns through, it absorbs what the gold transmits; that flux is lost to the hohlraum
   either way (what the hot backing returns through the gold is of order f² α_b ≲ 1% of
   the incident flux, *L*). The window's loss starts at the *gold's* burnthrough, not at
   the emergence. Both calibrations of the heat front to the 7.5 ns edge therefore start
   the loss too late: for 2 µm at 6.6–6.9 ns, 2.4–2.9 ns after `HeatFront`'s arrival
   (3.9–4.2 ns) and 1.6–2.1 ns after the 1-D onset (4.7–5.0 ns).
6. A tamped-foil term in `BurnthroughAlbedo` is not supported: the gold's front is not
   delayed by the backing, and the albedo does not depend on the backing. What the
   backing delays is the outside signal. The open question for the albedo is the
   transition's onset and rise under a NIF drive (item 4), whose `width` (0.4 ns) is
   calibrated on bare Nova foils at 260 eV; not changed here.

### Table: public ignition-scale drive
N221204's Dante 1 from 3.3 ns (Pak 2024 Fig. 1), before that N210601's (Zylstra 2022
Fig. 2) ×0.958 to join it, a linear rise over 0.2–0.8 ns (143, 203, 263, 283, 302,
311 eV at 2.5, 4.5, 5.5, 6.5, 7.5, 8.5 ns). Times in ns. Calibrations fitted on this
drive to 2.5 µm through at 7.5 ns: coefficient ×0.393, offset +2.89 ns (×0.434, +2.62 ns
on the drive ×0.95). "Front + backing": `HeatFront`, `BurnthroughAlbedo`'s transition and
`RosselandTransmission` drive the backing. "> 9": not through by the end of the run.

| Au (µm) | backing | `HeatFront` | coefficient cal. | offset cal. | 1-D onset | front + backing | 1-D emergence (A ×1.3–×0.7) |
|---|---|---|---|---|---|---|---|
| 1.5 | epoxy 35 µm | 2.84 | 5.67 | 5.73 | 3.71 | 5.52 | 6.14 (5.87–6.61) |
| 1.5 | Be 50 µm | 2.84 | 5.67 | 5.73 | 3.69 | 6.22 | 7.00 (6.58–7.67) |
| 2.0 | epoxy 35 µm | 3.87 | 6.57 | 6.76 | 4.67 | 5.87 | 6.57 (6.28–7.07) |
| 2.0 | Be 50 µm | 3.87 | 6.57 | 6.76 | 4.64 | 6.64 | 7.48 (7.05–8.18) |
| 2.5 | epoxy 35 µm | 4.61 | 7.50 | 7.50 | 5.21 | 6.24 | 7.00 (6.68–7.51) |
| 2.5 | Be 50 µm | 4.61 | 7.50 | 7.50 | 5.19 | 7.05 | 7.93 (7.48–8.66) |
| 3.0 | epoxy 35 µm | 5.08 | 8.42 | 7.97 | 5.60 | 6.62 | 7.41 (7.09–7.93) |
| 3.0 | Be 50 µm | 5.08 | 8.42 | 7.97 | 5.58 | 7.46 | 8.36 (7.90–> 9) |

Drive ×0.95: `HeatFront` 3.19, 4.22, 4.88, 5.32 ns; 1-D emergence with epoxy 6.46, 6.94,
7.41, 7.86 ns and with Be 7.41, 7.94, 8.43, 8.95 ns (1.5–3.0 µm); the calibrations' times
are those above to 0.1 ns (the coefficient calibration's exactly: a scaled drive scales
the drive integral, which the fitted coefficient absorbs).

**What separates the calibrations.** Neither depends on the backing; the 1-D emergence
does (50 µm of Be through 0.9–1.0 ns later than 35 µm of epoxy). Across thickness the
coefficient calibration spreads 2–3 µm by 1.85 ns, the offset by 1.2 ns, the 1-D
emergence by 0.85 ns. The decisive measurement is the emergence through windows of the
same gold with another backing (or none), then another gold thickness with the same
backing; a window shot's Dante with the backing's areal mass varied separates the
backing's burn from the gold's front directly.

### Eq. 3 in the convention of the Letter
Eq. 3 of the Letter, ΔP₂/P₀ ≈ (2.5 n A S₂/4πr²)(Δα/α)(F/(F + 1)), is the axial Legendre
moment of n small windows at 90° from the hohlraum axis: a_l/a₀ = −(2l + 1) P_l(cos θ_w)
(nΩ/4π) S_l f w, with Ω = A/r² the window's solid angle from the capsule, Δα/α = −f, S_l
the hohlraum's smoothing of mode l at the capsule and w = F/(F + 1) the re-emitted share
of the wall's radiance (*L*, derivation in the litreview §5). The 2.5 is 5 |P₂(0)|: as
printed the formula carries the sign of Δα/α, so its "−0.9%" is a magnitude. The
deficit on the waist is a **positive** P₂ of the drive about the hohlraum axis (more flux
at the poles), which pushes a prolate core rounder, as the Letter observed. In Briareus
terms the quantity is a₂/a₀ = √5 c₂₀/c₀₀ (real orthonormal Y_lm, `ylm_weights`). Two
windows at 90° are not axisymmetric about the hohlraum axis: c₂,±₂ carries three
quarters of the l = 2 power, and about the window axis the same pattern is a₂/a₀ =
−2 × (eq. 3's), the drive deficit along the window axis that the neutron image's
elongation along it reflects.

Closed-form check (`test/windows.jl`, "window modes against the closed form"): a sphere
of radius 3 mm with a uniform source, two windows of transmission 0.3 and 20° half-angle
at (90°, 150°) and (90°, 330°), a 0.6 mm absorbing capsule. Briareus's c₂₀/c₀₀ is within
1% and c₂,±₂/c₀₀ within 3–4.5% of the first-order closed form (window facets' exact
solid angles, S₂(x) from the capsule's cosine kernel, w from the sphere's energy
balance), converged over meshes 30 × 60 to 60 × 120; the residual in m = ±2 is the
capsule's shadow in the wall-to-wall exchange (order x² = 4%), which the closed form
leaves out. The m = ±2 share of the l = 2 power matches the closed form's 0.75 to 3%.

Public mock-up (scratch run): the HYBRID-E-scale hohlraum of `examples/n221204.jl`
(11.24 × 6.4 mm, 3.1 mm LEHs, 1.135 mm capsule, h = 0.25, capsule_h = 0.12),
`hybride_mock_pointing(nif_beams())`, a constant wall albedo α and conversion 0.8, two
2 mm windows at (90°, 150°) and (90°, 330°) burned through with eq. 1–2's transmission
for 2.5 µm (f = 0.235, albedo α(1 − f)), against the same hohlraum without windows:

| α | drive | c₂₀/c₀₀ | c₂,₂/c₀₀, c₂,−₂/c₀₀ | a₂/a₀ hohlraum axis | a₂/a₀ window axis | δ₂, δ₄ | w | eq. 3 with S₂(sphere), w | eq. 3's S₂ w |
|---|---|---|---|---|---|---|---|---|---|
| 0.88 | 0.976 | +0.467% | −0.397, +0.683% | +1.04% | −2.05% | 2.05, 0.56% | 0.754 | +1.20% | 0.363 |
| 0.80 | 0.983 | +0.413% | −0.361, +0.622% | +0.92% | −1.86% | 1.86, 0.51% | 0.674 | +1.07% | 0.322 |

(w: one minus the capsule flux with every albedo zero, i.e. the spots' direct share, over
the full flux. S₂ = 0.555 for a sphere at the capsule-to-wall ratio 0.355. "eq. 3's
S₂ w": the product that reproduces Briareus's a₂.) The Letter's −0.9% implies S₂ w =
0.313. **Once the sign and the axis are matched, Briareus and eq. 3 agree: +0.92 to
+1.04% against 0.9% (2–16% apart, α 0.80–0.88).** A sphere's S₂ overstates the
cylinder's by 13–14% (eq. 3 with the sphere's S₂ and Briareus's w is 1.07–1.20%); an
open hole (f = 1) gives the same ratios (the response is linear in f; a₂ 3.99–4.52%).
Scaled to 2 µm (f = 0.286) the mock-up's δ₂ is 2.5%, the roadmap's figure on the
study's pointing.

The core-shape comparison is no check of this size. Haines et al. (Phys. Plasmas 33,
112712 (2026), Fig. 9; section "Drive asymmetry inferred from post-shot capsule
simulations: Haines 2026" above) put the applied flux P₂/P₀ that
reproduces each shot's neutron image at −0.024 for the reference shot, +0.003 for the
window shot and +0.022 for a nominal no-window repeat (negative: prolate in their
convention; assumed a Legendre ratio, √5 c₂₀/c₀₀): the scatter between nominal repeats
(0.046) exceeds the windows' +0.01, though the window shot's offset from the reference
has the sign above. The Letter's neutron M2/M0 along the window axis is likewise a check
only; nothing here is fitted to it.

With the 1-D foil's slower transmission (item 4 of the account above: 0.62 of the
transition's transmitted energy over 4–9 ns with `RosselandTransmission`, about half the
constant default's) the time-averaged window signal would be
about half of the static numbers here (which use the constant 0.235); that is a statement about the transition, not
about eq. 3.

## The published picket fluence, the trough and the inner-beam corridor (2026-11-18) *derived*

ROADMAP items 36 (constrained re-pointing) and 10 (the bubble). A multi-phase re-pointing
study (scratch runs) held the outer cones' picket fluence as the published bubble
scalings define it, an energy over a fixed spot area, after finding that holding the
energy-weighted F̄ = ∫F² dA/∫F dA of "The swing and the bubble without new data" §5
penalized uneven per-quad picket powers and cost picket drive. This
section records the definitions from the papers, the quantities added for them, what
the two summaries do on public cases, and what the bubble link would give if it read the
published one. No default, exponent or calibrated constant changed; nothing was fitted;
PRE 109 Fig. 10 is a blind check only. Confidence: **H** primary text read (Izumi 2021
and 2023 reread today from OSTI 1773252 and 1971298), **L** ours.

### 1. The published quantity

| Source | What is divided by what | Where |
|---|---|---|
| Callahan et al., PoP 25, 056305 (2018) | outer picket energy over the one-side outer spot area, 53 mm² (Table I's 60 kJ design: √(60/(53 × 0.45)) = 1.586, the outer energy of both sides over one side's area, **L** on "both sides") | Eqs. 2–4, Table I (`hohlraum_dynamics.md` §1), **H** |
| Ralph et al., PoP 25, 082701 (2018) | outer picket energy over its FWHM, net of backscatter, over the outer phase-plate spot area | absorption metric, **H** |
| Izumi et al., PoP 28, 022706 (2021) | the outer picket fluence averaged over the wall, over the first ns (Sec. III, defining the FoM); the picket as the initial ~1 ns (Sec. II) | **H** |

- **Fixed area, with one exception.** Callahan's and Ralph's areas are the phase plates'
  nominal spots, fixed for a pointing family; their ratio between two designs is the ratio
  of the outer picket energies. Izumi 2021 counts its cone- and quad-split shots as a 27%
  lower overlapped outer-cone fluence on the wall at the same energy (Sec. III B and Fig.
  8a, with the tip at −120 ± 30 µm), i.e. a larger footprint, and reduces it to 20 ± 3% for
  the best P2 fit (Sec. III A, Fig. 5b). The same paper argues for an average over the
  azimuthal structure: the 50° spots are 24% more intense (smaller spots), but by the time
  the bubble clips the inner beams it has converged 2–3× and the inner quads overlap around
  the azimuth, so an average intensity is the appropriate measure (Sec. III B, **H**). So
  the published quantity holds the energy over a fixed area and treats a split of the quads
  as a change of area; it does not resolve spikes.
- **In Briareus terms** (`picket_irradiance(...; measure = :average)`, **L** for each
  choice): the outer beams' power *absorbed on the named parts* (the wall), over a fixed
  `area`. Light clipped on the entrance holes' outside, landing on the capsule or the LEH
  hardware, or leaving the target does not ablate the wall under the outer spots, so it
  is not counted; the papers count delivered energy, and on their shots clipping in the
  picket is small, so they do not separate the two. *Net of backscatter*: the powers as
  given are taken as net (`pulses`/`power` net of the measured picket backscatter), or
  `lpi = true` deposits the evaluation's effective beams, after the scenario's
  laser–plasma model at the time of that scenario (a `Backscatter`, a transfer). The default
  `area` is Callahan's 53 mm², so that the irradiance × the picket's length (ns) is his E/A
  in kJ/mm²; a ratio to a reference does not depend on it. On the N221204 mock pointing (h
  0.1 map) the 43.04 kJ of outer picket absorbed on the wall give 0.812 kJ/mm² over 53 mm²
  (60 kJ and 1.13 kJ/mm² in Callahan's Table I design), against F̄ 0.602 kJ/mm².
- `picket_ceiling(...; measure = :average)` holds it against a reference
  (`:picket_average`), in any phase whose powers the evaluation or the keywords give
  (picket, trough, second step). It is linear in the powers and moves with the pointing only
  through clipping and the power deposited off the wall: at fixed powers it is flat until
  the outer spots reach the entrance holes (the 1-parameter test: flat over the 50° cone's z
  4.6–5.4 mm, −1% at 5.8, −7% at 6.0, −18% at 6.2 on the 0.45 mm test hohlraum).

### 2. The exponents and the size of a 10% change

- **Core P2 against fluence: P2 shift ∝ F^0.5.** Callahan 2018: R_bub ∝ √(E/(A ρ_fill)) τ/R
  (Eq. 2), P2 ∝ −R_bub r_cap/R (Eq. 3), zero at the metric 1.2 for CF 33% (Eq. 4); 39
  shots, σ 8 µm, the stated goal ±5 µm (**H**). Izumi 2021 restates it as
  P2 − 270(CF − 0.333) = 85(1.2 − FoM) µm (Eq. 1) and keeps the square root by building
  its new FoM as x_bub^2.5 with x_bub ∝ E^0.2, FoM ∝ E^0.5 (Eq. 5, Sec. III B), slope
  −115 µm per unit FoM over its seven shots (Fig. 5b; −100 in Fig. 5a, **H**).
- **Bubble position ∝ F^0.2.** Izumi 2021 Fig. 8a (Sec. III B): 1.5× picket +160 µm
  (+9 ± 1.5% at x̄ ≈ 1.8 mm) gives E^0.2±0.03; the split shots' −27% fluence −120 ± 30 µm
  (−6 ± 1.5%) E^0.2±0.05; fill ρ₀^−0.2±0.03 (**H**). Ralph 2018's speed ∝ I^0.36 (data
  fit I^1/3–I^2/3) and Izumi 2023's ablation rate ∝ I^1/3 are the alternatives (§4 of the
  bubble section).
- **Per +10% of fluence** (**L**, from the published slopes at the P2 = 0 point,
  FoM ≈ 1.2): 85 × 1.2 × (√1.1 − 1) = −5.0 µm (Eq. 1); −6.7 µm with Fig. 5b's −115 (its
  redefined FoM's value at P2 = 0 is not given); Izumi 2018's measured 1.5× picket (−21 to
  −24 µm core P2, "Callahan 2018 and Izumi 2018" above) scaled as √: −4.6 to −5.2 µm. The
  study read −4 to −7 µm per +10% from the papers' figures (user, 2026-11-18; read from
  figures, not tabulated): consistent. Against tolerances: Callahan's ±5 µm goal
  (attributed to Izumi in the study; Izumi 2021's text states none, **H** for the
  search), and Kritcher 2024's −5 µm of P2 ≈ 1.9× yield on N230729 (`asymmetry_transfer.md`).
  So +10% of outer picket fluence is about one tolerance of core P2.

### 3. The trough and the corridor (Izumi et al., PoP 30, 032706 (2023), **H**)

- **The trough counts.** A sufficiently low trough lessens bubble growth by starving the
  bubble of ablated mass (Sec. II); on the 4-shock CH design the picket-launched bubble
  stalls in a 1 TW trough and is relaunched near 6 ns (Sec. III B, Fig. 9b); ballistic
  growth from the picket holds for < 15× foot-to-peak contrast. So the drive that grows
  the bubble includes the trough and the next step when they differ between designs;
  the study holds the published quantity in the picket, the trough and the second step
  (user's decision). `picket_irradiance`/`picket_ceiling` take any phase's powers.
- **Spikes in line with the 23.5° beams.** With azimuthal beam phasing (half the 50°
  quads on early, at 4× power) the 50° bubble's 8-fold pattern became 4-fold, its tips
  within 5.6° of the 23.5° beams (Sec. III C, Fig. 10c); the spike regions advanced 300 µm
  and the troughs fell 500 µm behind (Fig. 9b); the 23.5° beams' thin-wall signal fell
  and the core went from a2 −17 ± 3 µm to −28 ± 3 µm (Sec. III D, Fig. 11), which the paper
  attributes to the 23.5° beams absorbed at the bubble tips. The same section's
  ablation argument (mass ∝ I^1/3, so concentrating the 50° power on half the quads
  lowers the foot's ablated mass by 2^(−2/3) at fixed energy, Sec. II) shows that
  concentration does not grow the bubble as a whole; it grows it *where* the power is.
- **Consequence for the two summaries** (**L**). F̄ rises with concentration (overlap or
  uneven powers) everywhere on the wall; the fixed-area energy ignores it; neither
  locates it. A spike matters where an inner beam passes. Hence the corridor: the outer
  irradiance near the inner centerlines, `measure = :corridor`, and the hottest
  elements with their distance from the nearest inner centerline,
  `picket_hot_elements`.
- **The corridor's form.** Mean outer irradiance ∫w I dA/∫w dA over the named parts, w the
  union over inner beams (with power) of a quintic smoothstep in the 3-D distance from the
  element centroid to the beam's centerline (aim, direction): 1 within `corridor − edge/2`,
  0 beyond `corridor + edge/2`. Chosen over membership fixed at the reference because the
  inner beams move with the design (by a sizable fraction of the radius in a re-pointing)
  and the quantity is meant to follow them; the smooth edge keeps it C² in their aims (a
  hard 1.3 mm membership jumps as centroids cross it). Elements beyond every edge are
  skipped exactly (the step is flat there). Defaults: `corridor = 1.3` mm, the study's
  choice (an engineering default, not published: Izumi 2021's 1.3 mm is the distance from
  the wall to the inner beams' outer edge at the bubble's z, Sec. III B, a different
  quantity), `edge = 0.3` mm (wider than an element of a 0.2 mm map; on the 0.45 mm test
  hohlraum the 0.3 mm edge gives a mean 0.7% below the hard 1.3 mm mean, and 0.1 mm a mean
  0.13% below). Checks (`test/picket_fluence.jl`): against the hard mean by hand where no
  centroid lies on the edge (1e-12), the weight against the smoothstep union by hand,
  gradients against central differences, a one-parameter run ending at the bisection root.

### 4. Published against energy-weighted on public cases

Ratios to the reference pointing on the same map; outer cones only; `measure_a.jl`
(scratch). N221204: the mock HYBRID-E pointing, 43.04 kJ of outer picket equally per
beam, h 0.1 map. Tutorial: the dropped-bundle tutorial's 48 quads in the 5.75 mm
hohlraum, Young 2020's requested picket (0.25–1.25 ns, 1.318 TW per outer quad), h 0.2
map. "Uneven": per-quad log-normal picket scales of the outer quads (seed 1), mean 1.

| Case | published (`:average`) | energy-weighted F̄ | peak | corridor (1.3 mm) |
|---|---|---|---|---|
| N221204, outer cones z −0.5 / −0.3 mm | 0.9998 / 1.0000 | 0.9998 / 1.0000 | 0.9999 / 0.9999 | 1.2045 / 1.1340 |
| N221204, z +0.3 / +0.5 mm | 0.9972 / 0.9724 | 1.0015 / 1.0097 | 1.0000 / 1.0000 | 0.8598 / 0.7612 |
| N221204, unsplit quads (s 0) | 1.0000 | 1.1922 | 1.3913 | 1.0213 |
| N221204, s 0.225 / 0.75 | 1.0000 / 0.9941 | 1.0587 / 1.0869 | 1.1879 / 1.1952 | 1.0140 / 0.9745 |
| N221204, uneven rms 0.1 / 0.2 / 0.3 / 0.5 | 1.0000 | 1.0115 / 1.0499 / 1.1186 / 1.3634 | 1.27 / 1.58 / 1.93 / 2.76 | 1.004 / 1.007 / 1.010 / 1.014 |
| N221204, bundle 11's outer quads off, the rest ×1.067 | 1.0000 | 1.0368 | 1.0667 | 1.0143 |
| tutorial, bundle 11 freed, as is | 0.9375 | 0.9873 | 1.0000 | 0.9623 |
| tutorial, freed, the rest ×1.0667 (outer energy kept) | 1.0000 | 1.0532 | 1.0667 | 1.0265 |
| tutorial, its `limited` re-pointed design | 0.9375 | 0.9933 | 1.0026 | 0.9816 |
| tutorial, outer cones z −0.3 / +0.3 mm | 1.0000 / 1.0000 | 1.0000 / 1.0000 | 0.979 / 0.982 | 1.1893 / 0.7982 |
| tutorial, 50° quads at the 44° quads' z (5.20 mm) | 1.0000 | 0.9454 | 0.9841 | 0.9235 |
| tutorial, uneven rms 0.2 / 0.5 | 1.0000 | 1.0660 / 1.4610 | 1.69 / 3.29 | 0.997 / 0.992 |

- **Uneven powers part the two most.** At the same outer energy, per-quad scales of rms
  0.1/0.2/0.3/0.5 (log) raise F̄ by 1.2/5.0/11.9/36% on N221204 and 6.6/46% (rms
  0.2/0.5) on the tutorial: holding F̄ at the reference with uniform outer powers would
  cost 1.1/4.8/10.6/27% of the outer picket energy (6.2/32% in the tutorial). The scratch
  study's finding (picket drive lost under an F̄ ceiling with per-bundle picket scales)
  is of this kind. Dropping a bundle and restoring its outer energy on the
  others raises F̄ 3.7% (N221204) and 5.3% (tutorial) at the same published value.
- **Spot crowding** (splits, unsplit quads, one cone on the other's ring) moves F̄ by
  −5 to +19% and the published quantity by < 0.6%: only F̄ responds to overlap, which is the
  one effect Izumi 2021 measured on the bubble (the split shot).
- **z shifts** move neither (≤ 0.2%) until the spots clip on the holes (z +0.5: the
  published quantity −2.8%, the energy that leaves the wall; F̄ +1%). They move the
  corridor most (±20% per ∓0.5 mm on N221204, ±19% per ∓0.3 mm on the tutorial): the
  outer spots slide along the wall toward or away from the inner paths.
- The peak follows neither: it moves 2–3× F̄ with overlap and 5–23× it with uneven powers.

### 5. The bubble link reading the published quantity (item 10; a measurement)

`PicketLinkedBubble` reads F̄. Read instead E = 1.1 × (published ratio), on §5 of the
bubble section's benchmark variants (`bench_pub.jl`, scratch: the coarse full model,
measured backscatter, the h 0.1 map; the reference is the mock pointing at E = 1.1,
minimum −0.181 at 5.25 ns, upward zero 6.194 ns). Variants whose published ratio is 1
to < 1e-3 (z −0.5, −0.3, unsplit, s 0.225) give E = 1.1, i.e. the "E 1.1" column of that
table: unsplit 6.19 ns against the linked 6.00, s 0.225 6.19 against 6.15.

| Outer cones | published / mock | E published | upward zero, E published (shift) | E 1.1 | E linked to F̄ | minimum, E published (at) | 8 ns |
|---|---|---|---|---|---|---|---|
| z +0.3 mm | 0.9972 | 1.0969 (speed ×0.9994, t_c 6.089 ns) | 6.250 (+0.057) | 6.25 | 6.248 | −0.205 (5.25) | +0.076 |
| z +0.5 mm | 0.9724 | 1.0697 (×0.9944, 6.112) | 6.455 (+0.262) | 6.42 | 6.406 | −0.228 (5.25) | +0.067 |
| s 0.75 | 0.9941 | 1.0935 (×0.9988, 6.092) | 6.210 (+0.016) | 6.21 | 6.149 | −0.187 (5.25) | +0.080 |
| unsplit, s 0.225, z −0.5, z −0.3 | 1.0000 | 1.1 | as E 1.1 | 6.19, 6.19, 6.14, 6.15 | 6.00, 6.15, 6.14, 6.15 | | |
| PRE 109 Fig. 10 | | | 6.67 | | | −0.170 (5.60) | +0.095 |

- **The swing moves by less than the deposition itself changes.** Through the published
  quantity the bubble responds only where light leaves the wall: z +0.5 (2.8% clipped)
  puts the crossing 0.035 ns later than E 1.1 (6.455 against 6.42) and 0.049 ns later
  than the F̄ link (6.406), the minimum 0.001–0.002 deeper; z +0.3 and s 0.75 are within
  0.01 ns of E 1.1. The F̄ link instead moves the crossing −0.19 ns for unsplit quads
  and −0.04 to −0.06 ns for the splits, which the published link does not see.

- Read through the published quantity, the link does not respond to the one measured
  re-pointing-like change (Izumi 2021's split quads: −27% in the overlapped fluence,
  −120 µm of tip), and it moves the bubble the other way from F̄ when the spots clip
  (E down where F̄ rises). Fig. 10's crossing (6.67 ns) stays 0.2–0.7 ns later under either link.
- **Recommendation (user decision; no default changed):** keep `PicketLinkedBubble` on
  F̄ (its only published re-pointing-like test is the split shot, which F̄ reproduces
  and the fixed area cannot), and use the published `:average` as the study's ceiling
  in each early phase, as decided, with the `:corridor` reported or held where the
  outer light approaches the inner paths. A link on a different summary would be a new
  model choice to test against a bubble-tip datum at fixed energy, which no paper gives
  ("The swing and the bubble without new data" §5, What remains).

### What remains

- A published bubble trajectory for a re-pointing at fixed energy (unchanged), and one
  with uneven quad powers, to choose between F̄ and the fixed-area energy for the link.
- Whether the corridor's mean, or its peak, is the better proxy for Izumi 2023's spikes;
  its 1.3 mm and 0.3 mm are engineering defaults.
- `picket_irradiance`'s `:mean` and `:average` read alike: the energy-weighted measure
  is `:weighted` since 2026-11-18 (user), and its constraint `:picket_weighted`.

## A diffusion transmission for burned-through windows (2026-11-18) *derived*

ROADMAP item 28, the user's decision on the previous section's item 4: the default
transition of `BurnthroughAlbedo` stays (Nova's `width`, the constant gray `f`), its
docstring carries a caution, and a 1-D-informed transition is opt-in. Implemented:
`DiffusionTransmission` and `diffusion_transmission` (`src/target.jl`), tests in
`test/burnthrough.jl` ("diffusion foil transmission") and `test/history_state.jl`.
Existing results are unchanged (a new transmission type; the default is untouched).

### The law (*derived*, no fitted parameter)
- Steady planar diffusion through a foil whose Rosseland opacity varies as T^(−α) across
  it, the front at the drive temperature T_D and the rear radiating into vacuum
  (T_rear⁴ = t T_D⁴/2, the boundaries that give Lester et al.'s eq. 2): with x = T⁴/T_D⁴
  and s the optical depth at T_D, dx/ds = (3/4) t x^(−α/4), so
  (3τ/4) t = (4/(4 + α)) [1 − (t/2)^((4 + α)/4)]. α = 0 is eq. 2; for τ ≫ 1,
  t → 4/(4 + α) of eq. 2, 8/11 for Hammer & Rosen's α = 1.5 (their eq. 31). The rear,
  colder than the front, is more opaque. τ at T_D scales as `RosselandTransmission`'s
  (τ = 5 at 325 eV for 2.5 µm, Lester's point). The root is never above the gray value
  (h(u) = 1 − k − u + k u^p has h′ ≤ 0 on [0, 1] because k p = 1), so Newton's method
  from the gray value converges monotonically.
- Tests: α = 0 against `eddington_transmission`; α = 1.5 and 3 against an RK4 shooting
  solution of the ODE (to 1e-9); the 8/11 limit; the implicit derivative through dual
  numbers.

### What it replaces, and what it does not
- **Magnitude and its time dependence: replaced.** The fraction follows the flux
  temperature on each element and the foil's profile.
- **Onset: not replaced.** The 1-D foil's onset lag behind `HeatFront` grows with
  thickness on Nova's drive (1% transmission 0.20–0.45 ns after the front for
  0.77–3 µm), while Nova's measured half signal is 0.09–0.23 ns after it with no trend.
  A law that followed the 1-D onset would be 0.1–0.3 ns late on the data the default was
  calibrated on; an onset at a fixed ratio of the front's mass to the foil's
  (m_F/M ≈ 1.24, the measured mean) fitted to Nova scatters 0.09 ns rms against the
  default's 0.06. Nova's measurement outranks the 1-D model, so the option keeps the
  front's arrival and the 0.4 ns `width`: check (a) below is therefore met with the
  residuals of the default. (The Nova times are measured from 10% of the drive flux to
  50% of the rear flux's peak in two narrow bands, 210–240 and 430–570 eV; BUCKY
  report, fti.neep.wisc.edu fdm1017, **H**.)

### Checks
(a) Nova's half-signal times (0.77–3 µm): option and default identical, residuals
−0.11 to +0.03 ns (rms 0.06).

(b) The 1-D foil (`examples/tamped_foil.jl`, backed by 35 µm of epoxy), transmitted
energy over 4–9 ns, 1-D over each law (onset and `width` the default's):

| drive | Au (µm) | / constant default | / `RosselandTransmission` | / `DiffusionTransmission` |
|---|---|---|---|---|
| public ×1 | 1.5 / 2.0 / 2.5 / 3.0 | 0.533 / 0.536 / 0.538 / 0.532 | 0.617 / 0.624 / 0.619 / 0.606 | 0.760 / 0.784 / 0.789 / 0.779 |
| public ×0.95 | 1.5 / 2.0 / 2.5 / 3.0 | 0.501 / 0.501 / 0.499 / 0.487 | 0.620 / 0.621 / 0.611 / 0.591 | 0.768 / 0.784 / 0.781 / 0.763 |

Transmission history, 2.5 µm on the public drive ×1 (front 4.61 ns, 1-D 1% at 5.21 ns),
1-D / default / Rosseland / diffusion: 5.01 ns 0.000 / 0.235 / 0.155 / 0.119; 5.41 ns
0.069 / 0.235 / 0.174 / 0.135; 6.0 ns 0.108 / 0.235 / 0.189 / 0.147; 7.0 ns 0.129 /
0.235 / 0.203 / 0.159; 8.0 ns 0.147 / 0.235 / 0.221 / 0.174; 8.5 ns 0.157 / 0.235 /
0.222 / 0.175. From 0.8 ns after the 1-D onset the 1-D foil is 0.73–0.90 of the option
(0.90 late in the peak); the rest of the energy difference is the onset. On Nova's drive the 1-D foil
transmits less than all three laws (to 2 ns, 0.21–0.57 of the option's energy for
0.77–3 µm; at 2 ns 0.27 against 0.46 for 0.77 µm, 0.07 against 0.13 for 3 µm): no
measurement constrains the transmitted magnitude there.

(c) Blind: the Letter's extracted window flux, about 0.18 of the flux on the window near
8.5–9 ns (previous section, item 4; its solid-angle calibration pending). The option
gives 0.17–0.175 for 2.5 µm at the public drive's 305–311 eV, the 1-D foil 0.16, the
default 0.235. Nothing is fitted to it.

### What it changes for a hohlraum
Scratch run (`drive_history`, 8 threads, 35–38 s per run): the HYBRID-E-scale mock-up of
the previous section (h = 0.3, capsule_h = 0.15, `hybride_mock_pointing`, N221204's
delivered pulse shared equally by the 192 beams, gold `WallLossAlbedo` wall, conversion
0.8, `transition = 0.05`, `wall_history = false`), two 2 mm windows of 2.5 µm gold at
(90°, 150°) and (90°, 330°), against the same hohlraum with wall in place of the windows.
The windows burn through at 4.49–4.50 ns in all three. Window loss (TW) and the
capsule's change δ₂, δ₄ and a₂/a₀ about the hohlraum axis (%), default / Rosseland /
diffusion:

| t (ns) | P (TW) | Tr (eV) | window loss | δ₂ | δ₄ | a₂/a₀ |
|---|---|---|---|---|---|---|
| 4.6 | 227 | 243 | 0.76 / 0.50 / 0.38 | 0.34 / 0.22 / 0.17 | 0.09 / 0.06 / 0.05 | +0.17 / +0.11 / +0.09 |
| 5.0 | 353 | 276 | 7.63 / 5.88 / 4.60 | 2.12 / 1.62 / 1.27 | 0.58 / 0.45 / 0.35 | +1.08 / +0.83 / +0.65 |
| 5.8 | 443 | 295 | 9.99 / 8.41 / 6.63 | 2.15 / 1.80 / 1.41 | 0.59 / 0.49 / 0.39 | +1.10 / +0.92 / +0.72 |
| 7.0 | 444 | 297 | 10.31 / 8.78 / 6.93 | 2.16 / 1.83 / 1.44 | 0.59 / 0.50 / 0.39 | +1.11 / +0.94 / +0.73 |
| 7.6 | 444 | 298 | 10.44 / 8.92 / 7.04 | 2.17 / 1.84 / 1.45 | 0.60 / 0.51 / 0.40 | +1.11 / +0.94 / +0.74 |
| 8.6 | 269 | 262 | 6.23 / 4.49 / 3.50 | 2.16 / 1.54 / 1.20 | 0.59 / 0.42 / 0.33 | +1.10 / +0.79 / +0.61 |
| 9.0 | 173 | 233 | 3.90 / 2.40 / 1.86 | 2.15 / 1.30 / 1.00 | 0.59 / 0.36 / 0.28 | +1.10 / +0.67 / +0.51 |

The window loss over the pulse is 39.0 / 32.1 / 25.3 kJ. The option takes the windows' peak
δ₂ from 2.2 to 1.4% and δ₄ from 0.6 to 0.4%; the eq. 3 estimate (0.9%) lies between
the default's a₂ (+1.1%) and the option's (+0.73%). Before the windows burn through,
the three runs are the hohlraum without windows to 0.00%. With `wall_history = true`
(the default) the intact windows differ from the wall before burnthrough (δ₂ 0.44% at
4.6 ns in a coarse check): `drive_history` follows the temperature history of the wall itself
only for `WallLossAlbedo` elements, so a foil's `before` loss is evaluated at the time
since the pulse start (a brighter foil early in the peak). Not changed here (history
code); noted for the roadmap.

### history_state
`history_state` works unchanged: it reads each step's `transmission`, which `solve` computes
from the model. Tested: a drive history with a `DiffusionTransmission` window,
`keep_albedo`, and `history_state` at the last step gives each window element the option's
fraction at the flux temperature of that step (to the solve's tolerance), below the gray
value there.

## Intact windows follow the wall's history (2026-11-19) *derived*

ROADMAP item 28, the intact foil (user decision: "Can't think of a physical reason why an
intact window would be brighter than the wall"). With `wall_history` on (the default),
`drive_history` now routes a thin-foil element whose `BurnthroughAlbedo` `before` is a
`WallLossAlbedo` through the same per-element wall history as a wall element of that
loss (`src/history.jl`, `_FoilHistory` and `_foil_history_split` in `src/target.jl`, a
branch in `_albedos`). Intact it is that wall element; across the transition and after
it `α_before` is the history's, so the documented blend α = (1 − s f) α_before holds with
it; the stored energy it re-emits under a falling drive (`wall_release`) is the wall
element's scaled by α/α_before (1 − s f; with an explicit `after`, its blend). The heat
front, the transmission, a standalone `solve` (no history: `before` at the time since
the pulse start, as documented) and `wall_history = false` are unchanged; `history_state`
reads the history's recorded albedos and transmissions and does not need a change.

Decision (the release scaling): (a) the wall's emission × α/α_before, chosen: the foil's
radiosity is then (1 − s f) of the wall element's, the same factor as its albedo, and an
intact foil is the wall exactly; (b) unscaled (a burned-through foil would re-emit all of
its stored energy inward while transmitting s f of its flux); (c) no release once the
transition starts (a step in the foil's emission at burnthrough). The choice matters only
where a burned-through window is exposed to a falling drive (after the laser turns off): on
the mock-up below, the rows at 8.6–9.0 ns, where it acts together with the change of
`before` itself (the two were not separated).

Tests (`test/burnthrough.jl`, "intact foil follows the wall's history"): a small
hohlraum with a 2 µm window against the same mesh with wall in its place, through a
picket, a trough in which the wall re-emits, and a peak: capsule modes and the window
elements' albedos equal to 1e-12 at every step before the front arrives; with
`wall_release = false`, at every step, the intact albedo equals Lindl's loss at the
effective time rebuilt from the recorded temperatures (1e-6), the transmitted fraction
s f and the albedo (1 − s f) α_before; the release scaling of a cooling state; the
standalone `solve` at `time`.

Measured on the HYBRID-E-scale mock-up of "A diffusion transmission for burned-through
windows (2026-11-18)" (h = 0.3, N221204 pulse, two 2.5 µm windows at the waist,
`transition = 0.05`), now with `wall_history = true`, base (main before the change) and
new. δ₂, δ₄ and a₂/a₀ (hohlraum axis) are the capsule's change against the same
hohlraum with wall in place of the windows (%), default transmission / `DiffusionTransmission`:

| t (ns) | P (TW) | δ₂ base | δ₂ new | δ₄ base | δ₄ new | a₂/a₀ base | a₂/a₀ new |
|---|---|---|---|---|---|---|---|
| 3.0 | 58 | 0.871 / 0.871 | 0 / 0 | 0.239 / 0.239 | 0 / 0 | −0.437 / −0.437 | 0 / 0 |
| 4.0 | 92 | 0.628 / 0.628 | 0 / 0 | 0.173 / 0.173 | 0 / 0 | −0.311 / −0.311 | 0 / 0 |
| 4.6 | 227 | 0.961 / 0.961 | 0 / 0 | 0.264 / 0.264 | 0 / 0 | −0.476 / −0.476 | 0 / 0 |
| 5.0 | 353 | 0.330 / 0.307 | 1.132 / 0.581 | 0.090 / 0.085 | 0.310 / 0.159 | +0.168 / −0.151 | +0.568 / +0.292 |
| 5.4 | 434 | 1.262 / 0.341 | 1.920 / 1.084 | 0.345 / 0.092 | 0.526 / 0.297 | +0.637 / +0.174 | +0.964 / +0.544 |
| 6.2 | 433 | 1.564 / 0.693 | 2.002 / 1.183 | 0.429 / 0.189 | 0.550 / 0.325 | +0.794 / +0.353 | +1.013 / +0.598 |
| 7.0 | 444 | 1.708 / 0.877 | 2.046 / 1.252 | 0.469 / 0.240 | 0.562 / 0.344 | +0.870 / +0.448 | +1.039 / +0.636 |
| 7.6 | 444 | 1.785 / 0.971 | 2.070 / 1.287 | 0.490 / 0.266 | 0.569 / 0.353 | +0.910 / +0.497 | +1.054 / +0.657 |
| 8.2 | 405 | 1.956 / 1.128 | 2.116 / 1.307 | 0.536 / 0.309 | 0.581 / 0.359 | +1.004 / +0.580 | +1.085 / +0.670 |
| 8.6 | 269 | 2.433 / 1.504 | 2.237 / 1.288 | 0.663 / 0.409 | 0.612 / 0.353 | +1.276 / +0.790 | +1.171 / +0.675 |
| 9.0 | 173 | 2.797 / 1.752 | 2.330 / 1.231 | 0.758 / 0.474 | 0.634 / 0.336 | +1.494 / +0.937 | +1.242 / +0.656 |

- Before burnthrough the windows' signature vanishes: base δ₂ 0.25–0.96% (a₂/a₀ to
  −0.48%: the intact foils, evaluated at the time since the pulse start, re-emitted more
  than the wall around them); new 0.000 at every step.
- After it the effect of the windows alone remains, flat through the peak: δ₂ 2.0–2.1%
  (default) and 1.2–1.3% (`DiffusionTransmission`), δ₄ 0.55–0.58 and 0.32–0.36%; the base
  rows were that effect less the intact foils' opposite one, still decaying through the
  peak (it is why the base rose from 1.3 to 1.8% over 5.4–7.6 ns).
- Window loss over the pulse 30.45 → 30.31 kJ (default), 18.32 → 18.21 kJ (diffusion);
  the windows' incident flux history changes a little, and the burnthrough times move by
  +0.005 ns (mean 4.7670 → 4.7724 ns, range 4.762–4.774 → 4.767–4.779 ns).
- With `wall_history = false` (the earlier section's table, and every history in the
  diagnostic-windows tutorial) nothing changes.

Nothing changes elsewhere: histories without thin-foil windows are unchanged
(`examples/n221204.jl`, which has no `BurnthroughAlbedo`, at 0.8/0.4, `full`, `table`: CSV
identical to main's, `cmp`); the diagnostic-windows tutorial runs its histories with
`wall_history = false`, so its numbers and assertions are unchanged (run as a script against
this change).
