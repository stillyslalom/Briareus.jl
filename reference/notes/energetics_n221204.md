# Research notes: Wall loss, energetics and the N221204 benchmark

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

## Gold wall loss and thin-foil burnthrough (2026-10-16)

Implemented in `src/target.jl` (`WallLossAlbedo`, `HeatFront`, `BurnthroughAlbedo`).
Units: T in heV (100 eV), t in ns, areal mass in g/cm², 1 hJ/mm² = 10⁴ J/cm².
Confidence: **H** read in the primary text; **H-img** read from a rendered page image of
a scanned or garbled PDF; **M** secondary; **L** inference.

### Sources
- Hammer & Rosen, Phys. Plasmas 10, 1829 (2003), via the preprint
  (https://www.osti.gov/servlets/purl/15002150). **H-img**
- Lindl et al., Phys. Plasmas 11, 339 (2004), Sec. IV.A. **H-img**
- Shussman & Heizler, Phys. Plasmas 22, 082109 (2015), arXiv:1505.05524. **H**
- Rosen, "Fundamentals of ICF hohlraums" (2005), https://www.osti.gov/servlets/purl/885380.
  **H**
- Rosen, 33rd EPS Conf. Plasma Phys., P2.006 (2006). **H-img**
- Rosen, "Marshak waves: constant flux vs constant T" (1994),
  https://www.osti.gov/servlets/purl/104716. **H-img**
- Lindl, Phys. Plasmas 2, 3933 (1995): not read.
- R. S. Lester et al., Phys. Rev. E 114, L023201 (2026), doi:10.1103/k8dz-pjdx
  (https://journals.aps.org/pre/abstract/10.1103/k8dz-pjdx), p. L023201-2, eqs. 1–2
  and the text between them. **H** (text, checked against the rendered page)

### Gold power-law fits (Hammer & Rosen eq. 31) **H-img**
- Specific energy e = f T^β ρ^−μ with f = 3.4 MJ/g, β = 1.6, μ = 0.14.
- Opacity 1/K = g T^α ρ^−λ with g = 1/7200 g/cm², α = 1.5, λ = 0.2.
- The fits hold for 1–2 heV. The subsonic solution adds r = P/(ρe) = 0.25.

### Subsonic heat front, T_s = T₀ tᵏ (Hammer & Rosen eq. 112, q₀ = 5.5k) **H-img**
- m = m₀ T_s(t)^1.914 t^0.5156 and F = F₀ T_s(t)^3.346 t^−0.4115, with T_s the
  instantaneous surface temperature.
- q₀ = 0 (constant T): m₀ = 9.51e-4 g/cm² and F₀ = 3.44e-3 MJ/ns/cm². These are the
  `HeatFront` defaults.
- q₀ = 1 (k = 0.182): m₀ = 7.13e-4. Their Fig. 4b gives an ablated mass of about
  2.2–2.4e-3 g/cm² at 4 ns for T_s = t^0.182, from both HYDRA and the analytic curve.
- Shussman & Heizler's exact solution, same fits (eqs. 35–36) **H**:
  - constant T: m_F = 10.17e-4 T^1.91 t^0.52 and E = 0.59 T^3.35 t^0.59 hJ/mm²
  - constant flux (T ∝ t^0.123): m_F = 8.32e-4 T₀^1.91 t^0.75
- Rosen 2005 quotes m_F0 = 9.9e-4 and 7.4e-4 for k = 0 and 0.18. **H**
- Lindl 2004 eq. 4-24 has m₀ = 1.20e-3 T^1.86 τ^(0.54+1.86P) (1+3.44P)^−0.46, with T
  at 1 ns and T ∝ t^P. **H-img** At 300 eV and 1 ns this is about 20% more mass than
  Hammer & Rosen.

### Time-varying drive (*derived*) **L**
- `HeatFront` integrates m_F^(1/b) = c^(1/b) ∫ (T/1 heV)^(a/b) dt. This is exact at
  constant T. For T ∝ tᵏ it gives the self-similar exponent a k + b and the
  coefficient c (1 + k a/b)^(−b).
- Against the exact solutions, m(k)/m(0) is +2.2% for Hammer & Rosen at k = 0.182
  and +0.7% for Shussman & Heizler at k = 0.123. A single drive integral cannot carry
  the profile-shape change with drive history.
- Lindl's (1+3.44P)^−0.46 factor is worse: +6.7% at k = 0.182.
- Not validated for a drive falling after its peak (Hammer & Rosen Sec. IIc).
- No published subsonic formula covers an arbitrary T(t).

### Wall loss and consistency with `WallLossAlbedo`
- Lindl 2004 eq. 4-21: E_w = 5.2e-3 T^3.3 τ^0.62 MJ/cm² at constant T. Its time
  derivative over σT⁴ is eq. 4-31, the instantaneous albedo
  α = 1 − 0.32 T^−0.7 τ^−0.38 (the `WallLossAlbedo` defaults). Eq. 4-32 gives the
  time-averaged 1 − 0.52 T^−0.7 τ^−0.38. **H-img**
- Hammer & Rosen's flux gives 1 − α = 0.334 T^−0.654 t^−0.4115, and Shussman &
  Heizler 0.338 T^−0.65 t^−0.41 (*derived*). Over 1–3 heV and 1–10 ns these lie within
  about 10% of Lindl's loss, so the burnthrough front and the wall albedo come from the
  same physics.
- The 20% difference in heated mass between Lindl and Hammer & Rosen is a difference
  in the mean specific energy E_w/m_F: Lindl 4.3 T^1.44 MJ/g against Hammer & Rosen
  6.1 T^1.43 MJ/g (*derived*). `HeatFront` takes Hammer & Rosen. The Lindl values are
  in its docstring for comparison.
- Lindl's T is the incident (drive) temperature, and Hammer & Rosen's the surface
  temperature, with T_D⁴ ≈ T_s⁴ (1 + (1−α)/2). Briareus drives with the incident-flux
  temperature, a few percent difference for gold late in the pulse.

### Measured burnthrough (Nova) **H-img**
- Gold patches of 0.77–3 µm on the wall of a hohlraum near 260 eV, with T rising
  about as t^0.15 under a flat pulse (Lindl 2004 Fig. 4-6; Sigel et al., PRL 1990;
  Porter & Thiessen 1992).
- Half-signal burnthrough times, read by eye to ±0.03 ns:

  | Au (µm) | Time (ns) |
  |---|---|
  | 0.77 | 0.33 |
  | 1.05–1.1 | 0.55–0.58 |
  | 1.55 | 0.78 |
  | 1.93 | 0.86 |
  | 2.05 | 1.03 |
  | 2.27 | 1.14 |
  | 3.0 | 1.47–1.58 |

- `HeatFront` at 260 eV (t/1 ns)^0.15 puts the front 0.09–0.23 ns (mean 0.19)
  earlier (*derived*). This sets the default 0.4 ns transition width of
  `BurnthroughAlbedo`, which reaches half transmission at t_b + 0.2 ns.
- Lindl's analytic curve, x = 2.45 t^0.82 µm, is also about 0.2 ns early ("about a
  10 eV temperature difference").
- A constant 260 eV estimate is 3× too early.

### Modeling choices
- Before burnthrough the foil uses the thick-wall albedo. Ahead of a nonlinear
  diffusion front the cold layer does not change the flux into the wall, and the Nova
  data show no early transmission.
- After burnthrough (default since 2026-10-16) a fraction T of the incident flux is
  transmitted out of the hohlraum and the rest re-emits with the `before` albedo:
  α = (1 − sT) α_before, transmitted sT, with s the smooth step. See the next section.
- `after = 0` (a burned-through foil is an open hole) is an upper bound on the loss,
  kept as an option. With an explicit `after`, only the loss beyond the intact foil's,
  (α_before − α)HA, is booked as transmitted (the earlier default's accounting).
  `transmission = 1` gives the same albedo and books the whole loss as transmitted.

### Post-burnthrough transmission (Lester et al. 2026)
- Eq. 1: τ(L) = L / λ_Au(T_h), "where λ_Au is the Rosseland mean free path of gold at
  radiation temperature T_h (typically 300–350 eV)", in a "simplified gray-opacity
  model". **H**
- "At L = 2.5 µm, and assuming the Au density is ≈1 g/cm³, the resulting τ ∼ 5". **H**
- Eq. 2: F_trans = (4/3) F_inc / (τ + 2/3), "the planar Eddington approximation with a
  Marshak boundary condition", i.e. T = 4/(3τ + 2) = 4/17 = 0.235 at τ = 5
  (*derived*; the paper does not print T). **H**
- The paper uses eq. 2 to benchmark 1D radiation-hydrodynamics transmission estimates
  after burnthrough. It is a diffusion result: T > 1 below τ = 2/3, where
  `eddington_transmission` caps it at 1.
- Density and thickness: the paper pairs the as-built thickness L = 2.5 µm with λ_Au
  at an assumed (expanded) density of ≈1 g/cm³. Literally, ρκL at 1 g/cm³ over
  2.5 µm is 0.25 mg/cm², not the as-built foil's 4.8 mg/cm². Briareus keeps the
  paper's pairing as the reference point: τ = 5 for the areal mass of 2.5 µm of
  solid gold (19.3 g/cm³). **H** (the pairing) / **L** (the interpretation)
- Thickness scaling (*derived*, **L**): τ is linear in areal mass at fixed λ_Au and
  density (eq. 1), so τ = 5 × ρ·thickness / (19.3 g/cm³ × 2.5 µm). T = 0.328, 0.235
  and 0.169 at 1.7, 2.5 and 3.6 µm (τ = 3.4, 5, 7.2).
- Not modeled by the default: λ_Au's dependence on T_h and density (held at the
  paper's 300–350 eV evaluation; the gray mean free path grows with temperature), time
  dependence of T after burnthrough, absorption or re-emission by a low-Z backing,
  non-gold foils (pass `transmission`). **L**

### Temperature-dependent transmission (`RosselandTransmission`, 2026-10-17)
- τ(T) = τ_ref (T/325 eV)^−1.5: the exponent of Hammer & Rosen's gold Rosseland
  opacity fit, 1/K = g T^1.5 ρ^−0.2 with g = 1/7200 g/cm² (T in heV; eq. 31, fitted
  over 1–2 heV), and Lester's τ ∼ 5 at 2.5 µm as the reference point at 325 eV, the
  middle of their 300–350 eV. At the reference it is the constant default. *derived*
- Cross-check (*derived*, extrapolating the fit to 3.25 heV): the as-built areal mass
  of 2.5 µm of solid gold, 4.83 mg/cm², with the opacity at Lester's expanded
  1 g/cm³, gives τ = 4.83e-3 × 7200 × 3.25^−1.5 = 5.9, within 20% of their τ ∼ 5. This
  supports the reading of their pairing in the previous section: the as-built areal
  mass with the opacity of the expanded foil. At solid density K would be 1.8× larger
  (ρ^0.2), τ ≈ 11.
- The density of the expanded foil is held fixed (1 g/cm³ implicitly, through the
  normalization). In the windows tutorial (NIF-like drive, 2 mm windows) the foils
  burn through at wall temperatures of 150–285 eV and the wall peaks near 310 eV, so
  the temperature-following transmission is about a tenth smaller over the pulse than
  the constant one (flux-weighted ΔP2 1.48 vs 1.69% for 1.7 µm, 1.00 vs 1.12% for
  2.5 µm).

### Closing gaps (`GapClosure`, 2026-10-17)
- A gap is open (all flux on it transmitted) until a given closure time, then wall-like
  across a smooth step (0.4 ns default). The closure time is an input; no gap-closure
  speed was found in the sources read. **L**
- Tutorial result (*derived*): 50 µm gaps around 2.5 µm windows closing at 3 ns cost
  nothing measurable beyond the windows; closing at 7 ns (after burnthrough, before the
  peak) they add 0.4% capsule energy and 0.09% ΔP2, a fifth of a gap open throughout
  (2.0%, 0.55%).

### Restoring the drive around a window pair (2026-10-17) *derived*
- Two 2 mm windows of 2.5 µm gold at the waist (φ = 78°, 258°), NIF-like pulse, the
  windows tutorial's pointing: at the peak they raise P2 by 1.2% and add m = 2 of 0.9%.
- Per-quad powers within ±20% (asymmetry objective against the no-window Y_lm shape,
  drive held), linearized at the peak with the converged albedos: l = 2 departure 1.09
  → 0.07%, l = 1, 3 below 0.02%, l = 4 unchanged (0.15 → 0.19%); the 30° cone +17%,
  32 of 48 quads at a limit, +3.9% peak power.
- Per-quad-pair re-pointing (z and along the wall, ≤ 0.3 mm, fixed power, keep-outs):
  l = 1–3 to 0.02%, l = 4 unchanged, 0.14/0.09 mm rms; cannot restore the 3.4% drive.
  This replaces the earlier look's "re-pointing only trims m = 2" (two open 2 × 2 mm
  windows, a much larger loss).
- Over time, settings found at the peak overshoot before burnthrough (ΔP2 −1.3 to −1.4%,
  m = 2 1.1% of the opposite phase, 5 ns to burnthrough). Per-quad power trims switched
  on with the burnthrough (NIF quads take their own pulse shapes) hold P2 and m = 2 at
  the no-window values throughout (m = 2 ≤ 0.04%), for +3.3% laser energy.

## N221204 time-dependent benchmark (2026-10-23) *derived*

`examples/n221204.jl` (ROADMAP "The P2 swing" and "Absolute energetics on N221204"): capsule flux P2/P0(t) of an N221204-class
HYBRID-E hohlraum, each hohlraum-dynamics piece added in turn, against the HYDRA curve
of PRE 109, 025204 Fig. 10 (N221204, post-CBET, digitized ±0.001 where visible, up to
±0.024 where covered; notes in the 2026-10-22 PRE 109 reading). These by-eye values are
superseded by `n221204_fig10()` ("PRE 109 Fig. 10 digitization" below). The table keeps
them as they were run.

> The example's `dynamics` argument (6th) was removed 2026-10-26: the shot is now the 6th
> argument and the CSV names drop the dynamics token; run commands here and below are as run.

- **Inputs** (`n221204_benchmark.md` §5). 11.24 × 6.4 mm hohlraum, 3.1 mm LEHs, 1.135 mm
  capsule (**H**). Unpublished, so repository defaults: square corners, no diagnostic
  windows, `gold(albedo = WallLossAlbedo())` wall (conversion 0.75), `ablator()` capsule,
  the capsule held at its initial radius (HYDRA's is near 0.92 mm at 6 ns, Fig. 1a).
  192 beams at `hybride_mock_pointing()` (an illustrative stand-in, not the shot's
  pointing). DN24 Fig. 1 delivered inner/outer powers (**L**; 2.03 MJ, 442 TW),
  `cone_fraction_pulses`; Δλ 2.75 Å (1ω) on the inner cones, 0.917 Å at 3ω, no
  23.5°/30° or 44.5°/50° split. Wall `WallMotion()` (0.10 mm/ns) from the outer turn-on
  0.4 ns, 11 keyframes; lips `GapClosure(closure = 5.0, width = 0.5, after =
  WallLossAlbedo())` at 0.84 r₀ (the research-notes example closure, not N221204 data);
  `BubbleAbsorption(; radius = 3.2, turn_on = 0.4, wall)` at its defaults (t_c = 5.9 ns,
  `capsule_scale` 1). Steps of 0.25 ns, 0.5–9 ns, standard mesh h 0.28 / capsule_h 0.16.
- **Two transfers.** (a) The reduced `CBET` at the deck's settings (leh_z 5.62 mm,
  clamp 0.002 → 0.0075 across 3.9–4.1 ns, other plasma parameters at the model's
  defaults). (b) The transfer that turns the input cone fraction into Fig. 9b's
  post-transfer one (0.65 at 5.6 ns, 0.55 at 7, 0.52 at 8; the HYDRA CBET result for
  this shot, a `CalibratedTransfer` from the outer to the inner cones). (b) takes an
  input from the same paper as Fig. 10, but not from Fig. 10, and no model default was
  changed.
- **Capsule flux P2/P0** (standard mesh; each column adds one piece to the one before):

  | t (ns) | static, no CBET | (a) static | (a) +wall | (a) +lips | (a) +bubble | (b) static | (b) +wall | (b) +lips | (b) +bubble | Fig. 10 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | 2.0 | −0.017 | −0.006 | +0.028 | +0.028 | +0.028 | −0.067 | −0.043 | −0.043 | −0.043 | −0.002 |
  | 4.0 | −0.010 | +0.014 | +0.081 | +0.081 | +0.081 | −0.114 | −0.083 | −0.083 | −0.082 | −0.013 |
  | 4.5 | −0.020 | +0.024 | +0.092 | +0.092 | +0.094 | −0.148 | −0.127 | −0.127 | −0.120 | −0.091 |
  | 5.5 | −0.107 | +0.010 | +0.080 | +0.094 | +0.108 | −0.233 | −0.230 | −0.193 | −0.138 | −0.168 |
  | 6.5 | −0.105 | +0.008 | +0.097 | +0.110 | +0.158 | −0.209 | −0.178 | −0.143 | +0.030 | −0.057 |
  | 7.0 | −0.105 | +0.007 | +0.106 | +0.119 | +0.176 | −0.193 | −0.145 | −0.111 | +0.092 | +0.039 |
  | 8.0 | −0.104 | +0.006 | +0.128 | +0.140 | +0.203 | −0.179 | −0.093 | −0.061 | +0.165 | +0.090 |
  | 8.5 | −0.166 | −0.086 | +0.047 | +0.069 | +0.203 | −0.181 | −0.073 | −0.042 | +0.182 | +0.117 |

  A coarse mesh (0.8 / 0.4) gives the (a) columns to within 0.006; a finer one (0.2 / 0.1,
  table times only, 408 s on 16 threads) the (b) full model to within 0.001. The
  standard run takes ~35 min per transfer on 8 threads (keyframed runs ~10 min each).
- **Swing shape, (b) full model against Fig. 10:** minimum −0.175 at 5.0 ns (−0.172 at
  5.57 ns); upward zero crossing 6.26 ns (6.69); +0.165 at 8 ns (+0.090); mean over
  5.25–8.5 ns +0.048 (−0.001); peak to peak ≈0.36 (0.31). Without the bubble (wall +
  lips) the curve never crosses zero before 8.5 ns.
- **Capsule Tr peak** (313 eV drive Trad, AS24; Dante ≈311 eV at 8.2 ns): (b) static
  295.1, +wall 304.5, +lips 312.7 eV (all at 8.0 ns), +bubble 306.1 eV at 5.75 ns and
  302.1 at 8 ns (the moved inner power is deposited near the LEH); (a) 288.7 → 302.1 (lips) →
  299.3 eV. Without CBET 293.0 eV. The +lips value matching 313 eV is partly the
  undegraded laser (HYDRA degrades it 9–13%) against the missing windows and lip
  hardware, so it is not a calibration.
- **Diagnosis.**
  1. *The reduced CBET reverses the transfer.* At the deck's settings it moves power
     from the inner to the outer cones: post-transfer cone fraction 0.14 at peak (input
     0.33) against Fig. 9b's 0.65, with the 50° cone gaining ×2.1 and the 23° cone
     keeping 29%. A static probe at peak (6 ns): Δλ = 0 gives 0.013, Δλ = 0.917 Å 0.14,
     clamp 0.002 0.15, Δλ = 2.0 Å (3ω) 0.66, and `mach = 0` with 0.917 Å 0.80. The
     default Mach-1 axial outflow Doppler-shifts the outer beams to the receiving side
     and overwhelms the 0.917 Å detuning. With (a) the static P2 is near zero and every
     dynamic piece pushes it positive, so (a) cannot show the swing; the comparison is
     informative only with (b). This is the CBET item of ROADMAP "Shot-data calibration": two crossing
     zones per LEH (from outer to inner inside), plasma-state (flow, Te, ne) inputs and a
     multi-epoch `calibrate_cbet` against the published post-CBET per-cone powers
     (Fig. 9a/b, Fig. 1d) would fix it (done: "Two-zone CBET calibrated to PRE 109
     Fig. 9" below; the example now defaults to it).
  2. *With Fig. 9b's transfer the swing is reproduced in shape and timing of the
     minimum; the late side rises too early and too far.* Contributions at 8 ns: wall
     +0.086, lips +0.032, bubble +0.226; at 6.5 ns: +0.031, +0.035, +0.173. The bubble
     dominates the swing, as the caption of Fig. 10 states ("as the wall enters and blocks the
     inner beams"), but its inner absorption (48 / 46% at 6.5 ns, 85 / 83% at 8 ns for
     23° / 30°) sits above Ralph 2018's 51–62% at peak, and the wall has moved 0.56 mm
     at 6 ns against 0.35–0.40 mm at the waist in PRE 109 Fig. 1a (constant 0.10 mm/ns
     from the turn-on overstates the motion in the low-power foot). Both put the zero
     crossing 0.4 ns early and the late plateau +0.07 high. Data that would fix it
     (Next 10): the HYBRID-E inner-beam propagation/absorption history (Ralph-type
     inner-beam imaging or HYDRA bubble transparency vs time) for v_ref, τ_b, x_max; the
     emitting-wall radius vs time (SXI/wall-motion imaging) for the wall speed and its
     start; LEH SXI for the lip closure time.
  3. *The foot is offset* by −0.04 to −0.07 with (b) (static −0.07 at 2 ns, −0.11 at
     4 ns) where HYDRA stays at ≈0: the foot P2 is set by pointing and foot cone fraction,
     and the mock pointing and Fig. 9b's foot CF (read by eye, ±0.02 ≈ ±0.02 in P2) are
     not the shot's tuned foot. The moving wall alone adds +0.03 by 2 ns (0.16 mm).
- **Caveats.** Mock pointing; unpublished geometry (corners, windows, liner, LEH
  washer) at repository defaults; the capsule does not shrink (Fig. 10's evaluation
  radius and flux definition are not stated); no backscatter and no laser degradation;
  Fig. 10 is a HYDRA calculation, not data. Nothing here was tuned to Fig. 10.

## PRE 109 Fig. 10 digitization (2026-10-24) *derived*

PRE 109, 025204 Fig. 10 shows HYDRA capsule flux P2/P0(t), post-CBET, for N210808
(black), N220919 (red) and N221204 (blue). It is now stored as data in
`data/pre109_fig10.csv` (`n221204_fig10()`), replacing the eight by-eye N221204 points
that were inline in `examples/n221204.jl`.

- **Method.** The figure is a raster (874 × 868 px) with no vector paths.
  - Axes: fitted linearly to 18 time ticks (0–9 ns, residual < 1 px) and 21 P2/P0 ticks
    (residual < 1.5 px). The scale is 77.7 px per ns and 1760 px per unit P2/P0, so one
    pixel is 0.013 ns or 0.0006.
  - Colors: each pixel's optical density is unmixed into blue, red and black. A 5 px
    opening removes the thin black laser-power and cone-fraction lines.
  - Centerline: in each pixel column, a 12 px stroke (the measured width) is fitted to
    cover only the curve color and no background. Curves drawn over it count as
    unknown: red is drawn over black, and black over blue.
  - Each sample averages ±0.027 ns of columns, every 0.1 ns over 0.8–8.9 ns.
  - The two narrow spikes (N220919 at 8.40 ns, N210808 at 7.66 ns) are traced row by
    row along their flanks.
  - Each curve's extrema are added as extra rows. Their values are the stroke's outer
    edge minus half the stroke width; a column center is biased by up to 0.003 at a
    sharp tip.
  - The overlay of the samples on the figure was checked by eye.
- **Accuracy** (per-point errors are in the `u*` columns).
  - ±0.002 where a curve is visible.
  - 0.004–0.012 where another curve covers it. N221204 lies under N210808 or N220919 at
    1.5–1.6, 2.4–2.6, 3.5–3.7 and 5.0–5.2 ns, and over much of 6.6–8.9 ns. N210808 is
    uncertain to 0.004–0.006 in the foot (under the thin lines) and near its 5.5 ns
    minimum (under N220919).
  - N221204's 8.2–8.55 ns spike is hidden under N220919's, apart from its left flank. It
    is N220919's spike moved 0.03 ns earlier with its peak capped at 0.147 ± 0.008
    (error 0.015).
  - The start transients before 0.75 ns go off scale and are not digitized.
- **N221204 values.** Minimum −0.173 at 5.58 ns. Hump +0.031 at 3.50 ns, shoulder −0.098
  at 4.77 ns. Upward zero crossing 6.67 ns. +0.095 at 8.0 ns. Peak ≈ +0.147 at 8.36 ns.
- **Comparison with the old by-eye points** (new − old):

  | t (ns) | 2.0 | 4.0 | 4.5 | 5.5 | 6.5 | 7.0 | 8.0 | 8.5 |
  |---|---|---|---|---|---|---|---|---|
  | old | −0.002 | −0.013 | −0.091 | −0.168 | −0.057 | +0.039 | +0.090 | +0.117 |
  | new | −0.001 | −0.018 | −0.092 | −0.169 | −0.045 | +0.046 | +0.095 | +0.109 |

  - The largest differences fall on the steep flanks: +0.012 at 6.5 ns and −0.005 at
    4.0 ns, where 0.03 ns of timing is worth 0.01. At 8.5 ns the curve is covered, and
    the difference is within its ±0.015.
  - None of the benchmark's conclusions above change.

## Young 2020 validation (2026-10-24) *derived*

ROADMAP "The P2 swing" (validation) and "Implosion metrics and UQ". Validation, not tuning: no
calibrated default changed. Full texts were read; data digitized from them are
in `data/` with source headers.

**New data.** Young et al., PoP 27, 082702 (2020), N170827: requested per-beam pulses
(Fig. 2e, `young2020_pulses()`), the mode-1 Green's function (Fig. 5b and 5d read
independently, agreeing to ~1 km/s/%/ns; `young2020_green()`, `young2020_kernel()`,
integral ~117 km/s per %), the dropped-quad mode-1 histories (Fig. 10a) and velocities
(Fig. 10c markers; `young2020_dropped_quads()`), the dynamic model's albedos (Fig. 3a top,
`young2020_albedos()`). Nora et al., PoP 30, 092705 (2023) Table I (`nora2023_table1.csv`,
transcribed). All but Table I are rasters; axes from the panel frames with at least two
tick or gridline checks per axis, overlays checked.

**Young's setup (what the comparisons must match).** Section III (dropped quads, offsets,
tilt, windows) uses the *static* VisRad model of N170827: Table II's 11.3 mm × 3.1 mm
hohlraum, 1.82 mm LEH, 1.05 mm corners and 0.98 mm capsule with "hohlraum, LEH, and capsule
radii all reduced to 85%", albedo 0.8 (Au) / 0.3 (ablator), conversion 0.9, the *requested*
(azimuthally symmetric) pulses, no windows, and the velocity from the N170827 Green's
function integrated over the whole pulse (negative in the foot, peaked at the rise and at
5.8 ns). Mode 1 is "the vector norm of the three P1m components", c_lm = 100 a_lm/a00
(Briareus `mode1_vector`; Young's x/y normalization of the associated Legendre functions was
not checked; the z one is the same). The offsets are of two kinds: the whole target moved
against fixed beams, or the capsule moved inside a fixed hohlraum. N170827's pointing is not
published; Briareus uses the tutorials' pointing with z scaled by 11.3/10.1, and unscaled as
a check.

**Comparisons** (Briareus static model, h 0.5 / capsule_h 0.15; the test suite runs
h 0.6 with the same tolerances):

| | Briareus (z scaled; unscaled) | Young 2020 |
|---|---|---|
| dropped quad, peak mode 1 (%): 23.5 / 30 / 44.5 / 50 | 2.85 / 2.96 / 1.39 / 1.31 (2.64 / 2.90 / 1.52 / 1.41) | 2.72 / 2.89 / 1.50 / 1.30 (Fig. 10a) |
| mode-1 history rms difference, 0.5–7.75 ns (%) | 0.21 / 0.18 / 0.15 / 0.05 (0.14 / 0.16 / 0.06 / 0.11) | — |
| dropped-quad velocity (km/s), G of Fig. 5 | 336 / 349 / 161 / 152 (312 / 342 / 176 / 163) | 317 / 337 / 170 / 146 (Fig. 10c) |
| inner/outer velocity ratio | 2.20 (1.93) | 2.07 ("almost doubles") |
| whole-target offset, horizontal / vertical (km/s/µm) | 0.09 / 1.30 | 0.1 / 1.3 |
| capsule offset in the hohlraum, horizontal / vertical | 0.10 / 1.07 | 0.1 / 1.1 |

**Findings.**
- *Setup difference:* the ROADMAP "Implosion metrics and UQ" mismatches (inner/outer
  1.3–1.6, horizontal offset 0.60 km/s/µm) came from the implosion-metrics tutorial's
  hohlraum (10.1 × 5.75 mm, 0.9 mm capsule, albedo 0.85, a static velocity at 100 km/s/%). A
  horizontal whole-target offset moves the inner and outer footprints so that their mode 1
  opposes; in Young's geometry they nearly cancel at peak (the inner beams by themselves
  0.0056%/µm, the outer by themselves 0.0067%/µm, together 0.0012), on the tutorial's they
  do not (0.030, 0.009, together 0.0059). The foot, with another cone balance, does not
  cancel (~0.5% per 100 µm in the foot and ~0.1% at peak in Young's Fig. 14a).
- *Model error, fixed:* the mode analysis (`legendre_weights`, `ylm_weights` and the
  `Scenario` caches, hence `drive_history`) measured directions from TCC even when the
  capsule was offset. That mixes the capsule's P2 (−0.1 here) into a spurious mode 1 at
  first order in the offset: 0.0074%/µm of fake horizontal mode 1 (0.74 km/s/µm)
  against 0.00005%/µm about the capsule's center, and the vertical one had the wrong
  sign. They now default to `mode_center` (the center of the observed sphere; TCC
  otherwise), so centered capsules are unchanged.
- The 100 km/s/% default against Young's kernel: a flat 100 km/s/% over 4–7.75 ns gives
  11–17% less than G on these traces (∫G ≈ 117 km/s/%). The default is unchanged.

**Static vs dynamic (Young Sec. II A, Fig. 3).** Young's dynamic model (radii and albedos
of Fig. 3a against time) is built with `Keyframes` at every 0.25 ns step: wall 3.1 mm minus
Fig. 3a's displacement, LEH and capsule scaled to Table II, corners scaled with the wall,
`TabulatedAlbedo` for wall and capsule. N170827's delivered per-beam powers are not
published, so the mode-1 source is a proxy: a fixed 5% rms per-quad imbalance on the
requested pulses (and, separately, a dropped 23.5° quad). Results (h 0.5):

| | imbalance proxy | dropped 23.5° quad | Young (delivered N170827) |
|---|---|---|---|
| mean abs. difference of the mode-1 norm, all times | 0.18% (mean norm 0.48) | 0.73% (2.29) | 0.70% |
| same, peak 4.5–7 ns | 0.10% (0.34) | 0.72% (2.15) | 0.14% |
| velocity static / dynamic (km/s) | 52.2 / 40.7 (+28%) | 336 / 260 (+30%) | 40.5 / 38.7 (+4.6%) |
| direction difference | 2.5° | 0.8° | (75°, 265°) vs (67°, 270°), ~8° |

Decomposition (imbalance proxy, static 52.2 km/s): dynamic radii with static albedos
51.1; with only the capsule shrinking, 50.4; static radii with Young's dynamic albedos 41.7;
static radii with a wall albedo of 0.87, 41.2. So **the 85% radii stand in for the
moving wall, closing LEH and ablating capsule to 2% in velocity**, as Young found; **the
albedo does not**: his dynamic wall albedo is 0.85–0.88 through the peak against the
static 0.8, and the extra re-emission smooths the mode 1 by ~20% (the same physics as
in his foot, where the lower albedo raises it). Young reports the static velocity only
4.6% higher. One point is not resolved: his delivered imbalance has a time structure (larger
in the foot and rise, where the albedos are closer), which the proxy lacks; the
delivered per-beam powers would settle it. The mean-difference numbers are on Young's
scale (0.10 against 0.14% at peak for a similar peak mode 1) but depend on the source's
amplitude.

**Nora 2023 (ROADMAP "Implosion metrics and UQ"): the transcription, checked against the paper.**
- Eq. 4's high branch is printed with rounded coefficients (13 ln(Y + 1585) − 95). The
  fit drawn in Fig. 6 has a = 13.105, b = 1585.225, c = −95.539; since ln(Y + b) ≈ 7.4,
  the rounding moves the result by −0.24 (1.43 against 1.67 at Y_amp = 80). The two
  branches still do not meet (the low one is 2.08 at 80), and the paper fits nothing
  between 40 and 80 (burn fractions ~6%, where the burn moves into the dense fuel).
  Implemented with the Fig. 6 coefficients on 80–500; 40–80 is refused.
- Table IV lists the inverses (factors from the measured to the unamplified speed):
  1/Eq. 4 reproduces its burn column to 1e-3 (N130530 0.956, N210808 0.570, N220919
  0.580; test suite), and total = shape × burn (N210808: 0.963 × 0.570 = 0.549; the
  text's "184%" is 1/0.549 = 1.82).
- The shape term: Eq. 2 is the 2D fit with P2 as a fraction (1.31 P2 + 0.71 below 0.35,
  −0.65 P2 + 1.4 above, continuous to 4e-3; 0.71 for a round hot spot because the 2D
  points are shifted from the 3D ones, Fig. 5). The "80–120% for −40% < P2 < 40%" of
  the conclusion describes Table I (3D), not Eq. 2; that is where the old transcription
  went wrong. Both are implemented (`shape_amplification(P2)`,
  `shape_amplification(P2, θ)`). Table IV's shape column follows Table I for equatorial
  flows (N130530: 0.965 exactly) but not for polar ones from the listed angles (N210808
  would need θ ≈ 77° where its flow is at 166°; no angle matches N210531): the angle
  convention of Table IV is not recoverable from the paper.

**Callahan 2018 and Izumi 2018 (trends; the repository had no hot-spot transfer).**
Callahan's Fig. 6 line (digitized from its dotted fit, 29 pixels at the two ends): P2_hs −
2.7(CF − 33) = 101.9 − 85.0 x µm, x = √(E_picket,outer/(A_outer ρ_fill)) τ/R · r_cap/R, zero
at x = 1.199 (his Eq. 4: 1.2), σ 8 µm over 39 shots. The per-shot inputs are not tabulated,
so the curve cannot be recomputed shot by shot; Briareus compares sensitivities. Transfer:
Callahan's 2.7 µm per 1% cone fraction divided by the Briareus d(flux P2/P0)/dCF (−0.293%
per %, peak 4.5–6.5 ns) gives −9.2 µm of hot-spot P2 per % of flux P2 (a pole-hot flux makes
it oblate). The setup is a 10.1 × 5.75 mm Au hohlraum (albedo 0.8), 0.908 mm capsule, the
tutorials' pointing, N170827's pulses, `BubbleAbsorption` with its calibrated defaults
(turn-on 0.25 ns), h 0.45. Izumi 2018 (LLNL-PROC-751038) single-parameter SymCaps, core P2
against the basecamp's +5.6 µm (the assignment of Fig. 6b's points to shots follows the
image order of Fig. 6a, **L**):

| change | Briareus flux P2 (pp) → hot spot (µm) | Callahan's line at x ≈ 1.2 | Izumi 2018 measured change |
|---|---|---|---|
| 1.5× picket energy | +3.25 → −30 | −23 | −21 to −24 |
| 1.18× capsule | +1.31 → −12 (view factors alone −0.59 → +5.5) | −18 | −21 to −24 |
| 0.15 mg/cc fill | −2.56 → +24 | −42 | +10 to +15 |
| quad splitting | not modeled (the outer fluence change is not given) | — | +10 to +16 |

All three Briareus signs match Izumi, magnitudes within a factor of 2. Callahan's fill term
gets the low-fill shot's sign wrong (a point Izumi makes: the 0.15 mg/cc shot is an outlier
of the bubble-size picture). The check is not blind: the bubble opacity's density exponent
was fitted to Izumi 2021's 0.15/0.3 mg/cc anchors and the collision time to Izumi 2018's
separation times (basecamp, 1.5× picket, 1.18× capsule), so this checks the composition
(bubble, view factors and transfer), not the inputs. Ralph 2018's 18.5% anchor stays skipped
(unpublished pointing).

## Wall material, conversion and capsule energetics (2026-10-24) *derived*

ROADMAP "Absolute energetics on N221204": model changes, each with its data check. Full texts read (accepted
manuscripts): Döppner et al., PRL 115,
055001 (2015) (OSTI 1209256); Schein et al., PRL 98, 175003 (2007) (OSTI 908388); Rosen,
"Fundamentals of ICF hohlraums" (2005, OSTI 885380), §VI; Callahan et al., PoP 27,
072704 (2020); Zylstra et al., PRE 106, 025202 (2022). Confidence: **H** primary text,
**L** our arithmetic or construction, *calibrated* / *derived* as marked.

> **Superseded in part (2026-10-25).** The models and inputs stand; these values do not:
> - the DU factor f = 0.785 (§1) is **0.836**, recalibrated at the wall's effective age
>   ("Energetics defaults and the capsule share (2026-10-25)");
> - the "0.92 mm ablation front at 6 ns" (§3) is the ice–gas interface; the example's
>   capsule schedule is now the K24 radius history of N221204 ("The P2 swing: capsule
>   radius, backscatter and the collision (2026-10-25)");
> - α_c = 0.22 (§3) became 0, then Young 2020's published ablator albedo in the example;
>   `CapsuleAlbedo`'s src default stays 0.2 ("Capsule absorption, conversion and the
>   ablating-wall release (2026-10-25)");
> - `nif_conversion()` (§2) is the conversion of the drive-history presets `nif_gold()`
>   and `nif_uranium()`; `gold()` and `uranium()` keep 0.75 (same section);
> - §4's `energetics=new` and the example's `wall=`/`conversion=` switches are gone; the
>   table's numbers are of their day (current: "Capsule absorption, conversion and the
>   ablating-wall release" §2).

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| DU wall equivalence | "equivalent to an additional 6.5% of laser power (≈25 TW) at the peak" | Döppner 2015, citing Kline 2013 and Callahan 2012 | H |
| Same, as a drive factor | f_DU = 1.065 in T^3.3 ∝ E f_cap f_LEH f_DU, "in both simulations and experiments" | Callahan 2020 Eq. 5 | H |
| DU hohlraum build | 7.0 µm DU lined with 0.7 µm Au, 22.3 µm Au backing; the liner keeps the picket as in Au | Döppner 2015 | H |
| Au/DU pair, high foot (T-1) | N140311 Au 391 TW, η_L 0.857, 317 ± 4 eV; N140520 DU 393 TW, 0.863, 322 ± 4 eV | Döppner 2015 Table I | H |
| Cocktail/Au loss ratio | E/A: Au (0.598, 0.398) T^3.3 t^0.6, U₀.₆Dy₀.₂Au₀.₂ (0.604, 0.407) T^3.1 t^0.57; ratio ∝ T^-0.22 t^-0.05; "about 0.83" for NIF | Rosen 2005 §VI | H |
| Same, LASNEX | 97.6 / 90.8 / 88.6% at 185 / 275 / 305 eV (1 ns); 17.5% less wall loss → 10% less laser energy at 300 eV | Schein 2007 | H |
| Conversion in time | 0.7 t^0.12 makes Nova's 1 ns T ∝ t^0.15 (Lindl 2004 Fig. 4-4); 0.7 t^0.2 (Rosen); 0.85 t^0.2 (Schein) | energetics.md §2.2 | H |
| Conversion in intensity | 0.95 at 10¹⁴, 0.80 at 10¹⁵ W/cm² (spheres, backscatter corrected) | Dewald 2008 abstract | M |
| NIF hohlraum η | 85–90% at 235–345 eV | Olson 2012 abstract | M |
| N221204 split | wall 54.1% (Au lining 4.4), LEH 25.2%, gas 4.5%, capsule 16.2%; laser degraded 9–13% | PRE 109 App. 3 | H |
| N221204 backscatter | SBS 42 ± 11 kJ early in the peak, mostly inner cones; SRS 6 ± 1 kJ | AS24 | H |
| LEH pair | N191118 3.64 mm LEH, ≈465 TW; N210601 3.1 mm, ≈390 TW; Dante Trad "comparable until ~7.5 ns" (≈295 / ≈300 eV) | Zylstra 2022 Fig. 2 and text | H, H-img |

### 1. Au-lined DU wall (`du_wall_loss()`, `uranium()`)

`WallLossAlbedo` gains `factor` f, `factor_exponent` δ (at `factor_temperature` 300 eV),
applied to gold's clamped loss, and an optional gold `liner` (areal mass) during which
the loss stays gold's (an energy-content blend w = min(1, liner/m_F) with m_F the
constant-T `HeatFront` mass; *engineering*). The defaults (gold) are unchanged.

- δ = −0.2 (Rosen −0.22; Schein's three LASNEX points give −0.19).
- **f = 0.785, calibrated** (*superseded*: 0.836 at the wall's effective age, see the
  note above) to the 6.5% equivalence: Döppner's high-foot hohlraum
  (5.75 × 9.43 mm, 3.1 mm LEHs, 1.125 mm capsule at albedo 0.3), 390 TW × 0.86 coupling,
  conversion 0.9, a wall source in two bands (30% at |z| 0.5–1.7 mm, 70% at 3.0–4.2),
  at 13 ns, h 0.4: the Au wall needs 6.5% more power for the DU wall's capsule Tr.
  Range 0.778–0.789 over 12–14 ns and conversion 0.75–0.95 (*L*).
- **Check, Döppner's T-1 pair** (not fitted): at equal power the DU wall raises the
  capsule Tr by 4.9 eV (294.2 → 299.1); measured 317 → 322 ± 4 eV (+5 ± 6) at 391 / 393 TW.
- **Rosen's 0.83 / Schein's 0.825 give only 5.0–5.1%** here (3.9 eV): in this model
  gold's Lindl loss leaves the wall ~30% of the x-ray source at a NIF peak (α ≈ 0.94),
  far from Schein's wall-dominated accounting (17.5% → 10% implies ~57%). Pure DU
  should be no better than the cocktail, so the calibrated f below Rosen's cocktail value
  indicates that the Lindl gold loss is low at late times, not that DU beats the cocktail.
- On HYBRID-E (11.24 × 6.4 mm, same bands, 350 TW) the same f is worth 8.4–9.0% of
  power (6.0–6.5 eV) at 6–8 ns: the larger case-to-LEH ratio gives the wall a larger
  share. In the N221204 example (full dynamics, nif conversion, ablator) DU raises the
  peak capsule Tr from 303.3 to 310.6 eV.
- Liner: 0.7 µm Au (Döppner; HYBRID-E's is unpublished) is crossed in ~0.4 ns at 150 eV;
  at 300 eV, 8 ns it is 6% of the heated mass (PRE 109 puts 4.4 of 54.1% of the wall
  heating, 8%, in the lining; *L*, consistent in size).

### 2. Conversion as a model (`Conversion`, `nif_conversion()`)

η(t, I) = clamp(η₀(I) ((t − start)/1 ns)^s, range), η₀ constant or a table in log₁₀ I of
the element's absorbed flux. For materials with a number, the scalar path is unchanged
(tested); intensity tables are generic in the number type (ForwardDiff through `solve`
tested); `flux_response(...; laser)` linearizes an intensity-dependent conversion (secant,
exact at `laser`); `Scenario`'s flux response evaluates a time-dependent one at its time.

- **Closed-form check**: with the wall loss dominating (pinhole LEHs, reflecting capsule,
  Nova scale-1 cylinder 1.6 × 2.55 mm, 20 TW), T ∝ t^((0.38 + s)/3.3): model slopes
  0.111 / 0.147 / 0.171 for s = 0 / 0.12 / 0.2 over 0.3–1 ns (closed form 0.115 / 0.152 /
  0.176).
- **Lindl's Nova statement** (0.7 t^0.12 → T ∝ t^0.15) is *not* reproduced with open
  LEHs: on the scale-1 hohlraum with 1.2 mm LEHs the slopes are 0.052 / 0.086 / 0.108 for
  s = 0 / 0.12 / 0.2. The LEHs take ~36% of the source and flatten T(t); Lindl's
  analysis is wall-dominated. Matching T ∝ t^0.15 with open holes would need s ≈ 0.3
  (*L*), so s stays an input.
- **NIF preset** `nif_conversion()` = 0.7 t^0.12 ≤ 0.95 (*derived*): 0.80 / 0.85 / 0.90
  at 3 / 5 / 8 ns, Olson's 85–90% at NIF peak; the cap is Dewald's 10¹⁴ plateau and
  Swadling 2025's VisRad value. Energy-weighted over N221204 it is 0.868. The materials'
  default stays 0.75 (changing it would move every baseline, and see §4).
- Dewald's intensity table is available as `Conversion(intensities = (1e14, 1e15),
  efficiencies = (0.95, 0.80))` but is not used in the example: it is a disk/sphere
  plateau (hohlraum η exceeds it, Lindl Fig. 4-10), and the element-mean absorbed flux
  depends on the mesh.

### 3. Capsule albedo with an effective radius (`CapsuleAlbedo`)

α_eff = 1 − (1 − α_c)(R(t)/R₀)²: the fixed capsule mesh absorbs what a sphere of the
effective radius would at the same incident flux; `energy_balance` books it,
`energy_budget` integrates E_cap. α_c a number or a model (`TabulatedAlbedo` in time,
the new `TemperatureAlbedo` in Tr); default 0.2 (Swadling 2025's f_cap 0.8).

- **Schedule (example, *L*; *superseded*: the anchor is the ice–gas interface, and the
  default schedule is K24's, see the note above)**: HYDRA's ablation front, 0.92 mm =
  0.81 R₀ at 6 ns (PRE 109 Fig. 1a), along the shape of Young 2020 Fig. 3a's N170827 front, its time axis mapped
  piecewise linearly (Young's 0.81 R₀ at 5.75 ns → 6.0 ns; its end, 0.23 R₀ at 7.75 ns →
  N221204's bang time 9.53 ns): R/R₀ = 0.94 at 4.2, 0.86 at 5.0, 0.81 at 6.0, 0.66 at 7.0, 0.51 at
  8.0, 0.31 at 9.0 ns.
- **Calibration** (*superseded*, see the note above): α_c = 0.22 gives E_cap = 16.2% of the coupled energy (285 kJ of
  1762) in the full energetics configuration below (α_c 0.20: 16.5%). It is degenerate
  with the schedule and with the hohlraum temperature: in the configuration whose LEH
  share matches HYDRA's (v4 below, Tr 307.7 eV) α_c = 0.2 gives 12.8%, and 16.2% would
  need α_c ≈ 0 (or a larger effective radius than the ablation front). So "α_c ≈ 0.2"
  is consistent with Swadling only together with the too-hot drive of §4.
- Kr24b's 16.2% is read as a fraction of the laser energy in the simulation (after the
  9–13% degradation): the four fractions sum to 100% with no laser loss. As a fraction
  of the delivered energy the same run gives 14.1%.

### 4. N221204 benchmark: energetics and the Fig. 10 blind effect

`examples/n221204.jl ... full` with the new `key=value` options; full model (calibrated
CBET, wall, lips, bubble); laser 2.03 MJ by eye. Backscatter: SBS 16.7% of the inner
power over 5.5–6.5 ns (ramps 5.0–5.5 and 6.5–7.5) and SRS 1.14% over 5.5–9.0 ns, scaled to
42 and 6 kJ (shapes *L*); degradation 11% on every beam (the middle of 9–13%); both
before CBET (whose calibration used HYDRA's degraded powers). Coarse mesh (0.8 / 0.4)
unless "std" (standard 0.28 / 0.16, 10 min per run on 4 threads; coarse within 0.006 in
P2 and 0.1 eV).

Energy split, % of the coupled energy (after backscatter and degradation); wall = all
surfaces' heating (unconverted laser + absorbed x-rays; the bubble's light is conducted
to the wall, so HYDRA's wall + gas, 58.6%, is the comparison):

| Run | wall | LEH | capsule | peak capsule Tr (eV) | P2/P0 min (ns) | zero | 8 ns | rms vs Fig. 10 |
|---|---|---|---|---|---|---|---|---|
| HYDRA (PRE 109 App. 3; Fig. 10) | 58.6 | 25.2 | 16.2 | ~300 internal, Dante 313 | −0.170 (5.6) | 6.67 | +0.095 | — |
| before, std (Au, η 0.75, ablator 0.3 at R₀, full laser) | 57.6 | 22.2 | 20.1 | 301.8 at 5.5 ns | −0.168 (5.25) | 6.26 | +0.123 | 0.036 |
| v1 + backscatter + 11% degradation | 58.2 | 22.1 | 19.8 | 288.9 at 8.0 | −0.138 (5.25) | 6.19 | +0.119 | 0.042 |
| v2 = v1 + nif conversion | 50.8 | 25.9 | 23.3 | 303.3 | −0.136 | 6.19 | +0.117 | 0.042 |
| v3 = v2 + DU wall | 46.4 | 28.2 | 25.4 | 310.6 | −0.127 | 6.20 | +0.106 | 0.041 |
| v4 = v1 + capsule model (α_c 0.2) | 62.1 | 25.1 | 12.8 | 307.7 | −0.140 | 6.28 | +0.058 | 0.041 |
| v5 = v4 + nif conversion | 55.4 | 29.7 | 15.0 | 323.4 | −0.139 | 6.28 | +0.056 | 0.042 |
| after, std (`energetics=new`: v5 + DU, α_c 0.22) | 50.8 | 33.0 | 16.2 | 334.2 at 8.0 | −0.134 (5.25) | 6.31 | +0.046 | 0.044 |

Findings:
- **The split has two independent numbers** (wall = 1 − LEH − capsule with no other
  sink): the LEH share fixes the drive temperature, the capsule share then fixes
  α_c(R/R₀)². With everything published switched on (`after`), the LEH share is 33%
  against 25% and the capsule Tr peaks at 334 eV against HYDRA's ~300 (Dante 313):
  (25.2/33.0)^¼ × 334 = 312 eV, i.e. matching the LEH share and matching Dante are the
  same correction. The x-ray wall loss is 37.6% of the coupled energy; HYDRA's split at
  this conversion needs ~45%, a wall loss ~1.6× this model's at the lower temperature.
- **What overshoots**: Lindl's gold loss (α ≈ 0.94 at the HYBRID-E peak) times the
  calibrated DU factor, with η 0.87, and the shrinking capsule that removes the late sink
  (v1 → v4: +19 eV). Each piece is from data; together they are too reflective. The
  configuration that matches the LEH share (v4: Au loss, η 0.75) does so by the old
  constant conversion standing in for the missing wall loss. The sinks missing here are the He
  gas (4.5%), the LEH windows, the M-band (Dewald 2020, not gray), and HYDRA's NLTE wall
  emission. The gap is not closed by tuning: no default was changed to close it.
- **Backscatter and degradation** (v1) cool the peak by 13 eV and move the swing minimum
  from −0.168 to −0.138 (less inner power after CBET); they barely change the late P2.
- **Fig. 10 (blind)**: with the full energetics the minimum is −0.134 at 5.25 ns
  (Fig. 10 −0.170 at 5.6), the upward crossing 6.31 ns (6.67), +0.046 at 8 ns (+0.095);
  rms 0.044 against 0.036 before. The late P2 drops by 0.065 with the capsule model: the
  shrinking absorber turns the capsule mesh into a reflector (α_eff 0.80 at 8 ns, 0.92 at 9), and a
  hohlraum with less sink smooths more. Part of that is physical (the real capsule
  absorbs less late), part is the fixed-mesh approximation (the flux is still evaluated
  at R₀; Fig. 10's radius is not stated). Nothing was tuned to Fig. 10.

### 5. LEH lips on the Zylstra 2022 pair

Static solves at 7 ns, 6.4 × 11.24 mm, h 0.4, wall bands as above; the power the 3.64 mm
LEH needs to give the 3.1 mm LEH's capsule Tr at 390 TW (data: 465/390 = ×1.192 for
"comparable" Dante Trad, the 3.1 mm one 5 eV hotter, so the equal-Tr ratio is ≥ 1.19):

| Wall, conversion | no lips, static wall | no lips, moving wall | lips (0.857 r₀), static | lips, moving |
|---|---|---|---|---|
| Au loss, 0.75 | ×1.124 | ×1.131 | ×1.099 | ×1.107 |
| DU, nif | ×1.140 | ×1.147 | ×1.115 | ×1.122 |

Lips closing to Young's 0.857 r₀ for both sizes *reduce* the model's LEH sensitivity
(both holes lose the same fraction of their area, so the difference shrinks); closing both by Young's absolute 0.26 mm gives ×1.143
(DU, nif, moving). All variants sit below the data's ×1.19. Reaching it at equal
temperature would need the LEHs to carry ~55% of the 3.64 mm hohlraum's loss at 7 ns,
implausible against PRE 109's 25% for the whole pulse; so the two shots' different
histories (N191118's pulse ends at ~7.5 ns, N210601's continues; different feet) and
Dante's LEH view (it sees the outer-beam bubble region, Ch24), not the lips, are the
likely difference. The pair checks the trend only (*L*).

## N221204 late-time drive (2026-10-24) *derived*

ROADMAP "Absolute energetics on N221204". Two results of the N221204 benchmark (`examples/n221204.jl`) disagreed about
the late drive. (a) Default energetics: capsule Tr peaked at 301.8 eV at 5.5 ns and was
flat to 8 ns; synthetic Dante 1 (open LEH area) peaked at 306 eV at 6.25 ns; measured
Dante 1 rises to 312 eV at 8.3 ns. (b) `energetics=new`: capsule Tr 334 eV at 8 ns and an
LEH share of 33% against HYDRA's 25.2%. Neither Dante nor Fig. 10 was used to fit
anything below; both stay checks. Confidence: **H** primary text or figure vector data,
**L** our construction.

> **Current state.** The fix (`wall_history = true`) stands. The numbers below were run
> on that day's energetics (α_c 0.22, DU factor 0.785, no re-emission, the misread
> capsule anchor); `energetics=new`/`old` no longer exist (the new energetics are the
> example's only settings). Current N221204 numbers: "Capsule absorption, conversion
> and the ablating-wall release (2026-10-25)" §2–3. The "What remains" items were taken
> up in "Wall re-emission and the capsule radius", "Late-time drive: the wall, the
> plasma and Dante" and "The P2 swing: capsule radius, backscatter and the collision".

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| N221204 delivered power | 2.051 MJ; 440–448 TW from 5.45 to 8.05 ns, 404 TW at 8.2, 269 at 8.6, 133 at 9.2, off at 9.45 ns | Pak et al., PRE 109, 025203 (2024) Fig. 1 top, read from the accepted manuscript's vector paths (`n221204_pulse()`, `data/pak2024_fig1_n221204.csv`) | H |
| N221204 Dante 1 Tr | 164 eV at 3.3 ns, 238 at 5.0, 263 at 5.5, 275.5 at 6.0, 291 at 7.0, 310 at 8.0, peak 311.9 at 8.30, 311.3 at 8.5, 301.5 at 9.0 | same, Fig. 1 bottom (its inset agrees to 0.05 eV) | H |
| Lindl wall loss | E_w = 5.2e-3 T^3.3 τ^0.62 MJ/cm² at constant T; α = 1 − 0.32 T^-0.7 τ^-0.38 is its derivative | Lindl 2004 eqs. 4-21, 4-31 ("Gold wall loss and thin-foil burnthrough" above) | H-img |

The example's by-eye DN24 pulse (2.028 MJ) is within a few % of the measured one: the
rise is ~0.1 ns late (175 against 213 TW at 4.5 ns) and the droop matches. HYDRA's
Fig. 9a input powers have the same shape (degraded ~0.84× in the peak only).

### Cause

`WallLossAlbedo` is Lindl's constant-temperature loss, and `drive_history` evaluated it
at the time since the pulse start: at 5.5 ns and 300 eV it used the loss of a wall held
at 300 eV for 5.5 ns (α ≈ 0.92). A NIF wall spends 4 ns of foot at 100–200 eV, which
stores little energy on the peak's scale (the loss integrand goes as T^5.3), so at the
start of the peak the wall loses like a wall ~1 ns old (α ≈ 0.85), and its loss then
falls as t_eff^-0.38 through the peak. The constant-T evaluation made the model hot
early in the peak (Dante 296 eV at 5.5 ns against 263 measured) and flat after it; with
the published energetics it also left the wall too reflective (LEH share 33%).
Disagreements (a) and (b) are the same error. The pulse was not the cause (below).

**Fix** (`drive_history(...; wall_history = true)`, the default): Lindl's wall energy is
carried to a varying T as E_w ∝ [∫ (T/100 eV)^(a/b) dt]^b (a = 3.3, b = 0.62; the form
of `HeatFront`'s drive integral, exact at constant T). The loss is then the constant-T
loss at each element's effective time t_eff = ∫ (T(t′)/T(t))^(a/b) dt′, integrated
trapezoidally from the incident-flux temperature of the element (constant T before the
first step). For T ∝ t^k it multiplies E_w by (1 + 5.32k)^-0.62 (*L*; the same form for
the heat front is within 2.2% of Hammer & Rosen's exact solution at k = 0.18). It is
not validated for a falling drive (a self-similar wall never gives energy back).
`wall_history = false` restores the old evaluation; the tests that compare a march with
standalone `solve`s at τ use it. The example now uses the measured pulse (`pulse=pak`,
the default), split inner/outer by DN24's cone fraction.

### Breakdown

Coarse mesh (0.8 / 0.4), full model (+bubble). Each cell is capsule Tr / Dante 1 over the
initial LEH area / Dante 1 over the open area (eV).

| t (ns) | A: old (DN24 pulse, constant-T loss) | B: A + measured pulse | C: B + wall history (new default) | E: `energetics=new`, old | D: `energetics=new`, pulse and history | measured Dante 1 |
|---|---|---|---|---|---|---|
| 5.0 | 283.3 / 269.2 / 283.8 | 281.8 / 267.9 / 282.4 | 255.0 / 241.8 / 255.0 | 293.5 / 279.6 / 294.8 | 260.7 / 248.2 / 261.7 | 237.9 |
| 6.5 | 298.8 / 284.7 / 305.7 | 297.9 / 283.9 / 304.8 | 278.0 / 266.2 / 285.8 | 315.0 / 299.4 / 321.4 | 289.8 / 276.9 / 297.3 | 282.9 |
| 8.0 | 298.4 / 283.6 / 306.3 | 298.4 / 283.6 / 306.3 | 283.2 / 270.3 / 292.0 | 334.3 / 315.7 / 341.0 | 310.6 / 294.6 / 318.2 | 310.1 |
| 8.5 | 270.4 / 257.2 / 277.8 | 271.4 / 258.1 / 278.8 | 263.8 / 251.4 / 271.6 | 306.1 / 288.9 / 312.1 | 292.6 / 277.0 / 299.2 | 311.3 |
| LEH / capsule share (%) | 22.2 / 20.1 | 22.2 / 20.1 | 17.4 / 15.6 | 33.0 / 16.2 | 24.8 / 11.6 | HYDRA 25.2 / 16.2 |

The dynamics pieces, capsule Tr (each column adds one piece; C, then D):

| t (ns) | static, no CBET | static | +wall | +lips | +bubble |
|---|---|---|---|---|---|
| 5.0 | 242.6; 248.2 | 247.0; 251.9 | 254.5; 259.7 | 256.5; 262.1 | 255.0; 260.7 |
| 6.5 | 271.6; 281.1 | 275.4; 283.9 | 283.5; 292.1 | 287.3; 297.3 | 278.0; 289.8 |
| 8.0 | 279.7; 302.3 | 281.8; 304.1 | 288.6; 312.6 | 293.7; 320.6 | 283.2; 310.6 |
| 8.5 | 261.8; 285.4 | 263.1; 286.6 | 269.8; 295.2 | 274.7; 303.3 | 263.8; 292.6 |

With the constant-T loss (A) the static column was flat from 5.5 to 8 ns (288.8 → 292.9).

- *Pulse*: ≤ 1.5 eV from 5 ns on (8 eV at 4.5 ns, from the earlier rise). Not the cause.
- *Wall history*: −27 eV at 5 ns, −20 at 6.5, −15 at 8 (full model, default
  energetics); −33 / −25 / −24 with the new energetics. It turns the static drive from
  flat into one that rises through the peak, 4% (default) to 13% (new energetics)
  from 5.5 to 8 ns, and moves ~10% of the energy from the LEH and capsule into the wall.
- *Moving wall* +7–8 eV, *lips* +2–8 eV, *bubble* −1.5 to −10 eV (it moves inner power to
  the outer footprint near the LEH), CBET +2–4 eV: none of them changes the shape.
- *`energetics=new` against the default, both with history* (D − C): +6 eV at 5 ns, +12 at
  6.5, +27 at 8, +29 at 8.5. The rising conversion (0.7 t^0.12) and the shrinking capsule
  (`CapsuleAlbedo` along HYDRA's ablation front) keep the drive rising through the peak,
  as Dante does. The default energetics (constant η 0.75, the capsule at R₀) have
  neither.
- *Dante convention*: from 6 ns on, the measured history lies between the initial-area
  reading (lips dark) and the open-area one (lips as bright as the wall behind them):
  in D at 6 / 7 / 8 ns, 271 / 283 / 295 and 290 / 306 / 318 eV against 275.5 / 291.3 /
  310.1.

### After the fix (standard mesh 0.28 / 0.16)

`energetics=new` (coarse D agrees within 0.1 eV and 0.007 in P2):

| t (ns) | 4.5 | 5.0 | 5.5 | 6.0 | 6.5 | 7.0 | 7.5 | 8.0 | 8.5 | 9.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| capsule Tr | 216.3 | 260.9 | 280.0 | 285.8 | 289.7 | 297.7 | 305.2 | 310.5 | 292.5 | 261.6 |
| Dante 1, initial area | 212.7 | 248.7 | 264.7 | 271.3 | 277.3 | 283.8 | 290.1 | 295.0 | 277.2 | 246.7 |
| Dante 1, open area | 222.6 | 262.2 | 281.0 | 289.8 | 297.8 | 305.9 | 313.3 | 318.6 | 299.5 | 266.4 |
| measured Dante 1 | 203.4 | 237.9 | 263.1 | 275.5 | 282.9 | 291.3 | 301.7 | 310.1 | 311.3 | 301.5 |

Capsule Tr peaks at 310.5 eV at 8.0 ns (HYDRA ~300 internal, AS24 drive Trad 313). The
split is wall 63.6 / LEH 24.8 / capsule 11.6% of the coupled energy (HYDRA 58.6 with
gas / 25.2 / 16.2). With the default energetics (standard mesh): capsule Tr 255.2 / 278.4 / 279.9 / 283.1 / 263.7 eV at 5 / 6 / 7 / 8 / 8.5 ns (peak 283.1 at 8.0 ns), Dante 1 over the open area 255.4 / 281.1 / 288.8 / 292.3 / 271.9; split 67.0 / 17.5 / 15.6%.

Fig. 10 (blind) metrics, full model:

| Run | minimum | upward zero | at 8 ns | rms 2–8.5 ns |
|---|---|---|---|---|
| PRE 109 Fig. 10 | −0.170 at 5.6 | 6.67 | +0.095 | — |
| default, before (A, coarse) | −0.162 at 5.25 | 6.28 | +0.117 | 0.034 |
| default, after (C, coarse) | −0.226 at 5.25 | 6.26 | +0.143 | 0.046 |
| default, after (standard) | −0.235 at 5.25 | 6.25 | +0.150 | 0.050 |
| `energetics=new`, before (E, coarse) | −0.129 at 5.25 | 6.34 | +0.043 | 0.045 |
| `energetics=new`, after (standard) | −0.189 at 5.25 | 6.23 | +0.077 | 0.041 |

The swing grows with the history: the wall is less reflective early in the peak, so the
hohlraum smooths the beams' pattern less (the static P2/P0 at 5.5 ns goes from −0.26 to
−0.33 with the default energetics). With the published energetics the minimum and the
8 ns value now sit on either side of Fig. 10's, closer than before (−0.189 / +0.077
against −0.170 / +0.095); the crossing is still ~0.4 ns early (the bubble, see
"Hohlraum dynamics calibrated to published data"). Nothing was tuned to Fig. 10.

The energetics pieces, each on top of the default with history and the measured pulse
(coarse; capsule Tr / Dante 1 initial / open area in eV; LEH / capsule share of the
coupled energy; Fig. 10 minimum, crossing, 8 ns value, rms):

| Run | 5.0 ns | 6.5 ns | 8.0 ns | 8.5 ns | LEH / capsule (%) | Fig. 10 |
|---|---|---|---|---|---|---|
| C: default | 255.0 / 241.8 / 255.0 | 278.0 / 266.2 / 285.8 | 283.2 / 270.3 / 292.0 | 263.8 / 251.4 / 271.6 | 17.4 / 15.6 | −0.226, 6.26, +0.143, 0.046 |
| F1: + backscatter, 11% degradation | 246.3 / 234.4 / 247.1 | 264.5 / 253.9 / 272.6 | 273.5 / 261.4 / 282.3 | 254.8 / 243.1 / 262.6 | 17.3 / 15.2 | −0.192, 6.19, +0.147, 0.049 |
| F2: F1 + `nif_conversion` | 254.1 / 241.9 / 255.0 | 275.5 / 264.4 / 283.9 | 286.8 / 274.0 / 296.0 | 267.6 / 255.2 / 275.7 | 20.2 / 17.8 | −0.192, 6.19, +0.145, 0.049 |
| F3: F2 + DU wall | 259.1 / 246.7 / 260.1 | 283.0 / 271.0 / 290.9 | 295.1 / 281.2 / 303.8 | 274.6 / 261.3 / 282.3 | 22.1 / 19.7 | −0.179, 6.19, +0.130, 0.044 |
| F4: F1 + capsule model (α_c 0.2) | 247.4 / 235.4 / 248.2 | 269.6 / 258.2 / 277.2 | 285.2 / 271.4 / 293.2 | 268.7 / 255.0 / 275.5 | 18.9 / 9.3 | −0.195, 6.23, +0.093, 0.041 |
| D: all (α_c 0.22) | 260.7 / 248.2 / 261.7 | 289.8 / 276.9 / 297.3 | 310.6 / 294.6 / 318.2 | 292.6 / 277.0 / 299.2 | 24.8 / 11.6 | −0.182, 6.24, +0.073, 0.041 |
| all, α_c 0.2 | 260.5 / 247.9 / 261.4 | 289.5 / 276.6 / 297.0 | 310.3 / 294.4 / 318.0 | 292.4 / 276.8 / 299.0 | 24.7 / 11.8 | −0.181, 6.24, +0.074, 0.041 |

Backscatter and degradation cool the peak by ~10 eV; the conversion (+3 eV at 5 ns,
+13 at 8) and the shrinking capsule (+12 at 8 ns, +14 at 8.5, nothing at 5) add the late
rise; the DU wall adds 5–8 eV throughout. The capsule model takes the late P2 down by
0.05 (less sink, more smoothing), as before the fix.

### What remains

> **Superseded (2026-11-18).** This list predates the later sections and its numbers are
> those of that day's energetics. The capsule share is 13.8% on the final defaults
> ("Capsule absorption, conversion and the ablating-wall release (2026-10-25)" §2, `k24`
> schedule with Young's albedo; 13.9% on the coarse mesh at commit 85866ae, "The swing
> and the bubble without new data (2026-11-18)"), not 11.6% (nor the 12.3% of the next
> section's "What remains"); the energetics
> named below became the example's only settings; the late drive after the droop is
> taken up in "Wall re-emission and the capsule radius" and "Late-time drive: the wall,
> the plasma and Dante"; the open items are kept in ROADMAP item 12.

- **After the laser droop (from 8.05 ns)** the model falls with the power (−18 eV by
  8.5 ns, −49 by 9 ns) while Dante stays at 311 eV to 8.6 ns and measures 301.5 at 9 ns. A
  quasi-static radiosity with a self-similar wall loss cannot give stored wall energy
  back; re-emission from the hot surface layer of a wall under a falling drive is the
  missing physics. Not addressed.
- **Early peak**: over the open area the model is 15–25 eV above Dante at 5.0–5.5 ns
  (initial area +11 / +2): the model's rise is ~0.15 ns ahead of the measured one.
- **Capsule share** 11.6% against 16.2% with α_c 0.22, which was calibrated on the
  old, too hot drive. Even α_c = 0 would give only ~15%, so 16.2% needs a larger
  effective radius than the ablation-front schedule, not a lower albedo. The
  recalibration is still open (ROADMAP "Absolute energetics on N221204") and is not done here.
- **Default energetics** (Au, η 0.75, the capsule at R₀) now under-drive the peak
  (capsule Tr 283 eV at 8 ns, Dante 270 / 292 eV). `energetics=new` is the configuration
  that matches Dante's level and shape through 8 ns and the LEH share. The example's
  default is unchanged (the 2026-10-24 benchmark); switching it is a decision for the
  next round.

## Wall re-emission and the capsule radius (2026-10-24) *derived*

ROADMAP "Absolute energetics on N221204", the two open items of "N221204 late-time drive": (1) after the laser droop
(8.05 ns) the model's drive fell with the power while Dante 1 stays at 311 eV to 8.6 ns;
(2) the capsule share, 11.6% against HYDRA's 16.2%. Dante 1 and Fig. 10 are checks
here, never fit. All runs `examples/n221204.jl ... full energetics=new` (the full model:
calibrated CBET, wall, lips, bubble; DU wall, `nif_conversion`, backscatter, 11%
degradation, `CapsuleAlbedo` on HYDRA's ablation-front schedule), coarse mesh 0.8 / 0.4
unless "std" (0.28 / 0.16). Confidence: **H** primary text, **L** our construction,
*calibrated* as marked.

> **Superseded in part (2026-10-25).** The re-emission closure (§1) stands, but its
> `release` is **1.0**, calibrated on an ablating wall; the static-wall 0.77 below no
> longer applies ("Capsule absorption, conversion and the ablating-wall release
> (2026-10-25)" §1). The capsule "ablation front" of §2 is the ice–gas interface and the
> example's schedule is now K24's ("The P2 swing: capsule radius, backscatter and the
> collision"). Numbers here are on that day's energetics (α_c 0.22, DU factor 0.785).

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| Re-emission and the Milne boundary | re-emitted flux φ_B(1 − F/2φ_B), incident φ_B(1 + F/2φ_B), absorbed F; α = 1 − F/φ; emission from ~2/3 mfp into the Marshak profile | Rosen, "Fundamentals of ICF hohlraums" (2005, OSTI 885380), Lecture 1, eqs. 1.34–1.38 | H |
| Gold power laws | e = 3.4 T^1.6 ρ^-0.14 MJ/g, κ_R = 7200 T^-1.5 ρ^0.2 cm²/g (T in heV) | Hammer & Rosen, PoP 10, 1829 (2003), as quoted by Rosen 2005 | H |
| Capsule ablation front | 0.92 mm (0.81 R₀) at 6 ns, N221204 | PRE 109, 025204 Fig. 1a | H-img |
| Dynamic-model capsule radius and albedo, N170827 | front 0.991 → 0.733 mm at 6 ns → 0.224 at 7.75; quasi-1D HYDRA ablator albedo 0.01–0.16 to 6 ns, 0.30 at 7, 0.40 from 7.5 ns | Young et al., PoP 27, 082702 (2020) Fig. 3a (`young2020_radii()`, `young2020_albedos()`) | H-img |
| Capsule energy in the rocket model | E_cap ∝ ∫ T_r⁴ r_cap(t)² dt with r_cap the rocket model's shell radius | Callahan et al., PoP 27, 072704 (2020) eq. 13 | H |

Rosen's lectures give the re-emission of a driven Marshak wave, not a cooling law: no
published closure for a wall under a falling drive was found. The model below is our
construction, checked against our own 1-D diffusion solve (an independent reference:
it shares only the material power laws with the closure, and no other code or assumption).

### 1. Wall re-emission (`drive_history(...; wall_release = true)`, default)

The self-similar loss with history (`wall_history`) never returns energy: when T falls,
t_eff grows and the loss only shrinks. The closure (*L*, in `_wall_split`):

- the history loss is the heating response, evaluated at each element's running maximum
  T_h = max T (exactly the old loss while T never falls below it);
- a drop is the deficit D = φ(T_h) − φ(T) in the wall-energy potential φ = (T/100 eV)^a
  (a = 4 + p = 3.3), and the wall answers it with the self-similar kernel of the loss,
  superposed (Duhamel, linear in φ): net flux in F = L(T_h, t_eff) −
  f_r b C ∫ (t − t′)^(b−1) dD(t′), C = L t_eff^(1−b)/(b φ(T_h)) (constant, c σ(100 eV)⁴/b,
  for gold's unclamped Lindl loss); a drop D at t_d gives back f_r C D (t − t_d)^b by t:
  fast at first (the hot surface layer), then as the t^b of the heat wave;
- D is piecewise linear between steps and the kernel integral exact for it. The term of
  the step itself is stiff (Δ^(b−1)); as an emission it made `solve`'s albedo iteration
  diverge (the first try gave albedos of 10⁴), so below T_h it goes into the albedo, 1 − f_r C
  Δ^(b−1) φ(T)/H ∝ T^-0.7 (iterates like Lindl's), and the T-independent rest into a source
  term (`solve` gained a per-element emission from the albedo state; `sol.source`
  includes it, so B = source + αH stays exact). The effective wall albedo exceeds 1
  where the wall gives back more than it takes.
- For bookkeeping, `drive_history` books `source` as ηP and `absorbed` net of the emission;
  `released` (new) is the re-emitted power per component and `energy_budget` integrates
  it. The budget still closes to 1e-15 on a falling drive. `keep_albedo` returns the
  effective albedo α + E/H, so a flux response about it reproduces the step exactly (a
  secant linearization; `retune_drive`'s test failed on the raw albedo, whose pulse
  ends in a ramp down).

**Calibration of f_r (*calibrated*, `WallLossAlbedo(; release = 0.77)`), against 1-D
nonlinear radiation diffusion** (`examples/wall_reemission.jl`): static gold
(Hammer–Rosen power laws, ρ 19.3, 300 zones geometric from 1e-8 g/cm², implicit, Picard
in T⁴), Marshak boundary on the incident flux σT_D⁴. Its constant-T energy goes as
T^3.57 t^0.5 (similarity 3.55, 0.5), which the closure is given (coefficient from the
same run). Drive: 100 eV foot to 4 ns, 300 eV by 5 ns, flat to 8 ns, then down by 5 / 10
/ 30% over 0.5 / 0.5 / 1.0 ns. Energy into the wall over 8–9.5 ns (MJ/cm²; 0.149 stored
at 8 ns):

| Drop | reference | no release | f_r = 1 | f_r = 0.77 | f_r the reference needs |
|---|---|---|---|---|---|
| 5% | +0.0193 | +0.0234 | +0.0156 | +0.0191 | 0.759 |
| 10% | +0.0086 | +0.0177 | +0.0023 | +0.0089 | 0.779 |
| 30% | −0.0148 | +0.0094 | −0.0291 | −0.0153 | 0.762 |

The closure is linear in f_r for a prescribed T_D, so f_r is read off exactly; the same
0.76–0.78 fits all three drops, and with 0.77 the net flux follows the reference point by
point (10% drop, F/σT_D⁴ at 8.1 / 8.5 / 9.0 ns: +0.011 / −0.010 / +0.015 against +0.013 /
−0.009 / +0.014; 30%: +0.003 / −0.035 / −0.095 against +0.006 / −0.032 / −0.093). The linear
response (f_r = 1) overshoots by ~1.6×: a cooled surface is more opaque and throttles the
release. Without the release the wall keeps absorbing after the drop. The caveats are a
static wall (the ablating wall's b = 0.62 is not tested) and the exponents of one material;
the history loss itself is 8% low against the reference early in the peak (0.112 against
0.122 of σT_D⁴ at 5 ns) and within 2% from 6 ns.

**Checks.** Rising drives are unchanged (`test/energetics.jl`, "wall re-emission"). Step
size: coarse N221204 at 0.125 against 0.25 ns steps (f_r 1): capsule Tr 302.1 / 283.9
against 301.7 / 284.0 eV at 8.5 / 9.0 ns. The per-element net flux matches the formula
evaluated from the recorded temperatures to 1e-6 of σT⁴.

**N221204 after the droop (check, not a fit).** Capsule Tr / Dante 1 over the initial /
open LEH area (eV), coarse:

| t (ns) | 8.0 | 8.5 | 9.0 |
|---|---|---|---|
| no release (before) | 310.6 / 294.6 / 318.2 | 292.6 / 277.0 / 299.2 | 261.6 / 246.5 / 266.3 |
| release, f_r 0.77 (after) | 310.6 / 294.6 / 318.2 | 299.5 / 283.1 / 305.8 | 277.0 / 260.4 / 281.3 |
| release, f_r 1 | 310.6 / 294.6 / 318.2 | 301.7 / 285.1 / 308.0 | 284.0 / 266.8 / 288.1 |
| measured Dante 1 | 310.1 | 311.3 | 301.5 |

**It does not hold 311 eV.** The release takes the fall from 8 to 8.5 ns from −19 to
−12 eV (open area) and from 8 to 9 ns from −52 to −37 eV; Dante rises 1 eV and then
falls 9.8. Over the open area the model is 5.5 eV low at 8.5 ns (12.1 before) and 20 eV
at 9.0 ns (35 before). Through 8 ns nothing moves by more than 0.2 eV (elements whose T
dips while the drive rises, near the lips and the bubble), so the early-peak excess at
5–5.5 ns is unchanged (+24 / +18 eV over the open area). Split (% of coupled): wall 63.6
→ 62.8, LEH 24.8 → 25.5, capsule 11.6 → 11.7; 103 kJ re-emitted in all (gross). Fig. 10:
unchanged to 8 ns (−0.181 at 5.25, crossing 6.24, +0.073 at 8 ns, rms 0.041); 8.5 ns
+0.056 → +0.047 (Fig. 10 +0.109). The remaining late gap needs something else: the
energy in the hohlraum plasma (gas and blow-off, which the radiosity has no store for),
the capsule's shrinking sink (it is in), and the LEH closing faster after the laser
turns off (Chen 2020); none is modeled here.

### 2. Capsule radius and albedo from published data

The schedule is unchanged: HYDRA's ablation front for N221204 (0.81 R₀ at 6 ns) along
Young's N170827 front, time-mapped so its end falls at bang time (9.53 ns). That end
mapping is the slowest plausible shrink (anchoring only at 6 ns, a pure shift, ends at
8.0 ns and absorbs less), and it agrees with N221204's kinematics to ~0.05 R₀ (*L*: from
0.92 mm at 6 ns to ~0.3 mm at peak velocity, ~9.0 ns, needs a mean ~190 km/s rising to
~380). What was not from data was α_c = 0.22, calibrated on the old, too hot drive.

| Run (coarse, release 0.77) | capsule α_c | capsule Tr 5.0 / 8.0 / 8.5 | Dante open 5.0 / 8.5 | LEH / capsule (%) | Fig. 10 min (ns), zero, 8 ns, rms |
|---|---|---|---|---|---|
| R₀ mesh, CapsuleAlbedo radius (`capsule=model`) | 0.22 (calibrated) | 260.9 / 310.6 / 299.5 | 261.8 / 305.8 | 25.5 / 11.7 | −0.181 (5.25), 6.24, +0.073, 0.041 |
| same | Young's HYDRA α(t) (`capsule_albedo=young`) | 259.1 / 310.8 / 299.8 | 260.1 / 306.1 | 25.3 / 12.3 | −0.180 (5.25), 6.23, +0.072, 0.043 |
| shrinking mesh (`capsule=mesh`) | 0.22 | 260.9 / 310.7 / 299.7 | 261.9 / 306.1 | 25.5 / 11.7 | −0.170 (5.25), 6.28, +0.047, 0.045 |
| shrinking mesh | Young's α(t) | 259.2 / 310.9 / 300.0 | 260.2 / 306.4 | 25.2 / 12.3 | −0.169 (5.25), 6.25, +0.046, 0.047 |
| HYDRA / PRE 109 Fig. 10 / Dante | — | ~300 internal | 237.9 / 311.3 | 25.2 / 16.2 | −0.170 (5.6), 6.67, +0.095 |

- **The published albedo raises the share to 12.3%**, not 16.2%: Young's ablator
  albedo is below 0.16 to 6 ns (more absorption early, −1.8 eV at 5 ns, which moves the
  early-peak excess by that much) and 0.3–0.4 late, where the small radius makes it
  matter little. With α_c = 0 the share is ~15% (before). **16.2% is not reachable
  with absorption at the ablation front** under this drive; HYDRA's "absorbed by the
  capsule" must include absorption outside the front (the blow-off, W-doped, which
  fills the low-fill hohlraum's center) or a larger effective radius than the front.
  No published radius of that absorbing surface was found, so nothing is fitted.
- **The fixed mesh is an energetics-exact stand-in**: shrinking the mesh itself on the
  keyframes (`moving_hohlraum(...; capsule_shrink)`, new) gives the same E_cap (11.7 /
  12.3%) and capsule Tr within 0.3 eV, and costs nothing (the keyframes already rebuild
  the target for the moving wall: 264 against 271 s per coarse run). The approximation
  `CapsuleAlbedo`'s radius makes is in the flux *pattern*: it evaluates the capsule
  flux on R₀ and re-emits the unabsorbed share from the R₀ sphere, as if light crossed
  the shell's surroundings undisturbed.
- **P2: the shrinking sphere sees less of the late pattern**: +0.073 → +0.047 at 8 ns,
  +0.047 → +0.026 at 8.5 ns, and the minimum −0.181 → −0.170 (Fig. 10 −0.170), the
  crossing 6.24 → 6.28 ns; rms 0.041 → 0.045 (the late values move away from Fig. 10's
  +0.095 / +0.109). Early (≤ 5 ns, R/R₀ ≥ 0.86) the two agree within 0.002. Fig. 10's
  flux radius is not stated (HYDRA's is at the ablation front), so the shrinking-mesh
  P2 is the matching comparison if Fig. 10 is the flux at the front.

### Standard mesh (0.28 / 0.16)

Capsule Tr / Dante 1 initial / open area (eV); "before" is the late-time-drive notes'
standard run (no release); coarse agrees within 0.6 eV and 0.008 in P2.

| t (ns) | 5.0 | 8.0 | 8.5 | 9.0 | LEH / capsule (%) | Fig. 10 min (ns), zero, 8 ns, rms |
|---|---|---|---|---|---|---|
| before | 260.9 / 248.7 / 262.2 | 310.5 / 295.0 / 318.6 | 292.5 / 277.2 / 299.5 | 261.6 / 246.7 / 266.4 | 24.8 / 11.6 | −0.189 (5.25), 6.23, +0.077, 0.041 |
| release (`energetics=new`, now) | 261.0 / 248.8 / 262.3 | 310.5 / 295.0 / 318.6 | 299.4 / 283.4 / 306.1 | 277.0 / 260.5 / 281.4 | 25.6 / 11.7 | −0.189 (5.25), 6.23, +0.077, 0.042 |
| + shrinking mesh, Young's α(t) | 259.3 / 247.2 / 260.7 | 310.9 / 295.4 / 319.0 | 300.0 / 283.9 / 306.7 | 277.5 / 261.1 / 282.1 | 25.2 / 12.3 | −0.176 (5.25), 6.24, +0.048, 0.047 |
| measured Dante 1 / HYDRA / Fig. 10 | 237.9 | 310.1 | 311.3 | 301.5 | 25.2 / 16.2 | −0.170 (5.6), 6.67, +0.095 |

Split with the release: wall 62.7 / LEH 25.6 / capsule 11.7% of the coupled energy,
103 kJ re-emitted (gross).

### What remains

- In the late drive, the model with the release falls 12 eV from 8 to 8.5 ns where Dante
  rises 1; the hohlraum plasma's stored energy and the LEH after the laser turns off are
  the candidates (neither modeled).
- The capsule share is 12.3% with published radius and albedo, against 16.2%; the
  absorbing radius outside the ablation front is the open input. (Superseded: 13.8% on
  the final defaults, "Capsule absorption, conversion and the ablating-wall release
  (2026-10-25)" §2; the release item was taken up there, §1.)
- The early peak (+18 to +24 eV over the open area at 5–5.5 ns, +1 to +11 over the
  initial area) is not moved by either change (the published capsule albedo takes
  1.7 eV off).
- The release closure is calibrated on a static wall; an ablating-wall (hydro) reference
  would test b = 0.62.

## Zylstra 2022 LEH pair, time-dependent (2026-10-25) *derived*

ROADMAP "Absolute energetics on N221204": the LEH-size pair run with each shot's measured pulse, the remaining
N221204 benchmark targets (bang time, hot-spot P2 and velocity), and the example's
default energetics. Dante, Fig. 10, the hot-spot data and HYDRA's split are checks
here; nothing was fitted to them. Confidence: **H** for the primary text or the figure's
vector data, **L** for our construction.

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| N191118 (3.64 mm LEH) delivered power | 1.623 MJ; peak 470.6 TW at 5.6 ns; 455–466 TW to 6.6 ns, falling from 7.0, half power at 7.65, off at 8.35 ns | Zylstra et al., PRE 106, 025202 (2022) Fig. 2 top, read from the accepted manuscript's vector paths (LLNL-JRNL-835348 p. 5; `zylstra2022_fig2()`, `data/zylstra2022_fig2.csv`) | H |
| N210601 (3.1 mm LEH) delivered power | 1.663 MJ (text: 1.7 MJ); peak 398.9 TW; 388–397 TW from 5.4 to 7.55 ns, half power at 8.25, off at 9.05 ns | same | H |
| Dante Tr, N191118 | 237.7 at 5.0 ns, 279.2 at 6.0, 293.0 at 7.0, peak 294.8 at 7.25, 281.4 at 8.0, 238.2 at 9.0, 226.0 at bang | same, Fig. 2 bottom | H |
| Dante Tr, N210601 | 233.8 at 5.0, 273.3 at 6.0, 288.5 at 7.0, peak 301.2 at 7.9, 292.3 at 8.5, 277.7 at 9.0, 253.8 at bang | same | H |
| Bang times (SPIDER) | N191118 9.37 ns, N210601 9.48 ns | same, the dashed lines | H |
| Δλ (1ω) | 1.25 Å (N191118), 1.8 Å (N210601) | Zylstra 2022 text | H |
| N221204 bang time | 9.53 ± 0.07 ns | K24E Table I | H |
| N221204 hot spot | P2 0.3 ± 0.1 µm; velocity 36 ± 6 km/s toward (θ, φ) = (151°, 54°) | P24 Fig. 2, K24 Table II | H |

Not published for the pair, so taken from N221204 (*L*): the inner/outer split (DN24's
cone fraction, with its rise and droop moved onto those of each shot by a piecewise-linear
time map through `pulse_rise` and the half-power fall), the pointing (the mock HYBRID-E
pointing; N191118's outer beams sat farther from the midplane), the capsule (the
N210808 design, 1.1285 mm, for both), the CBET (the joint fit `hybride_cbet`, its tables
shifted by the rise difference, −0.03 and +0.03 ns), and the capsule schedule (its end
moved to each shot's bang time). Backscatter is unpublished for the pair and is left
out; `energetics=new` otherwise (DU, `nif_conversion`, `CapsuleAlbedo` α_c 0.22, 11%
degradation).

### 1. The pair: did the histories and Dante's view close the gap?

The static result (§5 of "Wall material, conversion and capsule energetics"): the
3.64 mm LEH needs ×1.10–1.15 the power for the 3.1 mm LEH's capsule Tr, against "×1.19 or
more" in the data. That data number compared the peaks (≈295 eV at 7.3 ns, ≈300 at
8.0 ns) at different times, after N191118's laser had started to fall. **At equal times
the data show the 3.64 mm shot the hotter one** through the rise and peak: Dante
N191118 − N210601 = +3.9 / +7.8 / +5.9 / +6.1 / +4.5 eV at 5.0 / 5.5 / 6.0 / 6.5 / 7.0 ns
at 1.12–1.18× the power, −4.4 eV at 7.5 ns once its laser falls. The equal-Tr power
ratio this implies at 6.0–6.5 ns is 1.07–1.09 (P ∝ T^3.3–4, *L*); at 5.5 ns 1.05–1.07.
So the gap was mostly the reading of the data. The static model's ×1.10–1.15 sits at or
just above it.

Time-dependent, full model (+bubble), standard mesh (0.28 / 0.16), `energetics=new`.
Cells are capsule Tr / Dante 1 over the initial LEH area / over the open area, then
measured Dante (eV):

| t (ns) | N191118 | measured | N210601 | measured |
|---|---|---|---|---|
| 5.0 | 258.6 / 254.1 / 267.8 | 237.7 | 252.0 / 244.2 / 257.5 | 233.8 |
| 5.5 | 280.6 / 272.9 / 289.7 | 267.0 | 275.1 / 262.9 / 279.1 | 259.2 |
| 6.0 | 286.5 / 278.1 / 297.0 | 279.2 | 281.2 / 268.6 / 286.9 | 273.3 |
| 6.5 | 290.4 / 282.4 / 303.2 | 286.4 | 285.4 / 274.0 / 294.2 | 280.4 |
| 7.0 | 293.0 / 284.2 / 306.3 | 293.0 | 290.9 / 278.3 / 299.9 | 288.5 |
| 7.5 | 279.4 / 268.8 / 290.2 | 292.8 | 296.0 / 282.4 / 305.0 | 297.1 |
| 8.0 | 252.8 / 240.9 / 260.2 | 281.4 | 288.2 / 273.7 / 295.6 | 301.1 |
| 8.5 | 202.1 / 188.6 / 203.7 | 261.9 | 263.4 / 247.9 / 267.8 | 292.3 |
| 9.0 | 178.9 / 166.9 / 180.3 | 238.2 | 225.6 / 210.2 / 227.0 | 277.7 |
| LEH / capsule (% of coupled) | 31.4 / 12.1 | | 25.3 / 12.2 | |

The pair difference, N191118 − N210601 (eV; coarse mesh 0.8 / 0.4, within 0.5 eV of the
standard mesh), data Dante against model capsule Tr / Dante initial / Dante open:

| run | 5.0 | 5.5 | 6.0 | 6.5 | 7.0 | 7.5 |
|---|---|---|---|---|---|---|
| data (Dante) | +3.9 | +7.8 | +5.9 | +6.1 | +4.5 | −4.4 |
| `energetics=new` | +6.6 / +9.6 / +10.2 | +5.5 / +9.6 / +10.2 | +5.3 / +9.0 / +9.6 | +5.1 / +7.9 / +8.5 | +2.1 / +5.3 / +5.7 | −16.7 / −14.2 / −15.3 |
| `energetics=old` | +6.4 / +9.7 / +10.2 | +5.5 / +10.0 / +10.6 | +5.5 / +9.5 / +10.2 | +5.6 / +8.7 / +9.3 | +3.2 / +6.6 / +7.1 | −12.9 / −10.6 / −11.5 |
| new, `lag=0.15` | +7.2 / +10.0 / +10.5 | +5.8 / +9.8 / +10.4 | +5.7 / +9.4 / +10.1 | +5.3 / +8.2 / +8.8 | +3.6 / +6.7 / +7.3 | −10.3 / −7.4 / −7.9 |
| new, `lips=computed` | +6.0 / +12.9 / +10.1 | +4.9 / +12.9 / +10.0 | +4.5 / +12.1 / +9.2 | +4.0 / +10.9 / +7.8 | +0.8 / +8.2 / +4.8 | −18.1 / −11.5 / −16.3 |

Findings:
- **The histories close most of the gap.** On the reading of the instrument itself (Dante)
  the model puts the 3.64 mm shot 8–10 eV above the 3.1 mm one at 5–6.5 ns, against 4–8 eV
  measured: 2–4 eV too much, ~1% in T, ~3–4% in power (*L*), against ~6% before.
  The sign is the same as the static result's (the model's larger LEH loses slightly too
  little). The difference does not depend on the energetics (old and new agree within
  1 eV): it is geometry, not wall or conversion physics.
- **The view of Dante widens the model's difference**: the capsule
  sees +5 eV, Dante +9. Through the larger hole Dante 1 sees deeper into the hohlraum,
  so the 3.64 mm shot measures hotter relative to its capsule. Comparing the model's capsule
  Tr with the measured Dante would hide the residual.
- **Computed lips widen it further** (+11–13 eV on the initial area): `LipClosure`'s
  absolute depth (0.24 mm) closes the 3.1 mm hole by a larger fraction, so Dante's
  initial-area reading darkens more on N210601. The prescribed relative closure (0.857
  r₀) does better here; Young's 0.26 mm on 3.64 mm and Schneider's 0.25 on 3.1 mm cannot
  separate the two forms (see "Computed LEH closure").
- **When the laser falls the difference flips sign too hard** (7.5 ns: −14 against
  −4.4): N191118's model Tr drops with the power while its Dante stays up for ~0.3 ns more.
  A 0.15 ns lag halves that (−7.4). This is the late-time drive problem below.
- Absolute level, rise and peak: the measured Dante lies between the initial-area and
  open-area readings on both shots (N191118: +0.8 / +19 eV mean, initial / open,
  5–7 ns; N210601: −4.1 / +15 eV, 5–7.5 ns, coarse), as on N221204.

### 2. After the laser falls, and the drive at bang

On all three shots the model's drive falls far faster than Dante's after the power
drops. Model minus measured (initial / open area, `energetics=new`, standard mesh):
N191118 −41 / −21 eV at 8.0 ns, −73 / −58 at 8.5, −71 / −58 at 9.0; N210601 −27 / −5 at
8.0, −44 / −24 at 8.5, −68 / −51 at 9.0; N221204 −28 / −5 at 8.5, −41 / −20 at 9.0. At
bang time Dante measures 226.0 eV (N191118) and 253.8 eV (N210601), the model 154.7 / 167.1
and 185.5 / 200.4. The ordering Zylstra points to ("by bang time Trad is significantly
higher on N210601") is reproduced: model +31 eV (initial area), data +28. The level is
not.

Size of the gap (*L*): the LEH x-ray energy that Dante's Tr implies over the fall (σT⁴ over
both initial LEH areas, Lambertian) exceeds the model's by 125 kJ (N191118, 7.25–9.37 ns),
89 kJ (N210601, 7.75–9.48) and 71 kJ (N221204, 8.0–9.4). The LEHs carry 25–31% of the
model's losses, so holding the hohlraum at Dante's temperature would take ~250–450 kJ more
x-ray source after the laser falls, 15–25% of the laser energy. That is several times the He
gas's 30–60 kJ ("Computed LEH closure") and more than the wall already gives back in the
model (`wall_release`: 137–272 kJ). A 0.15 ns lag of the source (~60 kJ stored) recovers
only 10–15 eV of it (N191118 at 8.5 ns, coarse mesh, −59 instead of −73). So the tail is not
a reservoir of the published size. What remains is either physics a gray, quasi-static
radiosity lacks over a 1–2 ns coast (a hot, optically thick LEH plasma or gold blow-off in
Dante's view; the wall's re-emission when the drive is fully cut, which `release` = 0.77 was
calibrated for only up to 30% drops), or a Dante time base or analysis that differs after
laser-off. The pair, with N191118's laser off 1 ns before bang, is the cleanest data for it.
Open.

### 3. Bang time and the hot spot (ROADMAP "Implosion metrics and UQ" metrics)

`examples/n221204.jl metrics=on`: the steps run on to bang time, the full model records
every Y_lm, and the metrics use the drive window from the start of the rise
(`pulse_rise`, 4.11 ns on N221204) to bang (9.53 ns). Standard mesh:

| | `energetics=new` | `energetics=old` | data / HYDRA |
|---|---|---|---|
| hot-spot velocity (flat kernel, 100 km/s per %) | 2e-12 km/s | 1e-12 km/s | 36 ± 6 km/s at (151°, 54°) |
| mean flux P2/P0, rise to bang | −0.0149 | +0.0057 | — |
| mean flux P2/P0, rise to 8.9 ns | −0.0163 | −0.0016 | Fig. 10 −0.0147 |
| mean flux P2/P0, 4.5–6.5 ns | −0.1089 | −0.1402 | Fig. 10 −0.1160 |
| hot-spot P2 shift from HYDRA's at −9.2 µm per % (4.5–6.5 ns) | −6.5 µm | +22.3 µm | measured 0.3 ± 0.1 µm |
| Fig. 10 minimum / 8 ns / rms (blind) | −0.189 / +0.077 / 0.042 | −0.235 / +0.150 / 0.050 | −0.170 / +0.095 |
| capsule Tr / Dante 1 (initial, open) at bang | 233.6 / 217.5, 235.0 eV | 208.6 / 194.5, 210.1 | no Dante after 9.4 ns |

- **Velocity**: zero within round-off. The mock pointing and the hohlraum are mirror
  symmetric and there are no windows, so the model has no mode 1. N221204's 36 km/s
  comes from sources that are not in it: per-quad power imbalance (DN24 Fig. 1 left,
  not attributable per quad), capsule thickness mode 1 (0.13 µm), the three equatorial
  diagnostic windows and the waist manufacturing gaps (H24), none of them published in a
  usable form. This check cannot pass with public inputs. It can rule out the mode 1 of the
  drive itself as the cause only if those inputs are known.
- **Hot-spot P2**: with no shipped hot-spot kernel, the comparison is with HYDRA's flux
  history, whose design gave a near-round hot spot (0.3 µm measured). Over the whole
  window the published energetics match Fig. 10's mean to 0.002 (old: 0.013 more
  pole-hot). Over the peak window that the −9.2 µm per % transfer was derived on
  (Callahan 2018's line on another platform, good to a factor of ~2), `energetics=new`
  gives a −6.5 µm (oblate) shift and the old energetics +22 µm. Neither is zero within
  the transfer's uncertainty band; the published energetics are 3× closer.
- **Bang time as the end of the window**: extending the window from 8.9 ns to bang moves the
  mean flux P2 by 0.0014 (new): the coast's flux is low and its P2 small, so the window end
  hardly matters for P2. It matters for the drive. At bang the model's Tr is 70 eV below
  Dante on the pair (§2), and on N221204 the capsule Tr is 234 eV at bang against ~276 eV
  that HYDRA gives N210808 there (Kr22 Table I; N221204 not tabulated). This is the same
  late-time deficit. As for the coast time, the measured pulse's half-power fall is at 8.75
  ns, and 9.53 − 8.75 = 0.78 ns is AS24's coast time exactly (a consistency check of the
  pulse data, not of the model).

### 4. `energetics=new` as the example's default: decided, flipped

Checks independent of anything fitted (standard mesh; the pair's `energetics=old`
numbers are coarse):

| check | `energetics=new` | `energetics=old` |
|---|---|---|
| Dante rise through the peak (convention-free): N221204 5→8 ns, measured +72.2 eV | +46.2 initial / +56.3 open | +28.3 / +36.8 |
| N191118 5→7 ns, measured +55.3 | +30.1 / +38.5 | +19.4 / — |
| N210601 5→7.5 ns, measured +63.3 | +38.2 / +47.5 | +23.5 / — |
| Dante level, initial area (the instrument's convention), mean / rms over rise and peak: N221204, N191118, N210601 | −4.8 / 9.3, +0.8 / 8.9, −4.1 / 9.3 | −18.0 / 23.0, −10.0 / 16.0, −17.0 / 21.6 |
| same, open area | +14.8 / 15.5, +19.2 / 20.1, +14.6 / 15.5 | +0.7 / 11.5, +7.6 / 12.5, +0.7 / 10.6 |
| Dante after the fall, initial area mean | −35, −53, −47 | −60, −68, −68 |
| Fig. 10 (HYDRA, blind) rms, minimum, 8 ns | 0.042, −0.189, +0.077 | 0.050, −0.235, +0.150 |
| hot-spot P2 via the transfer | −6.5 µm | +22.3 µm |
| HYDRA LEH / capsule share (25.2 / 16.2%) | 25.9 / 11.7 | 18.2 / 16.3 |
| pair difference | same | same |

(The level rows use coarse runs; N221204's standard mesh is within 0.5 eV.)

The published energetics win on Dante's shape, the one check that does not depend on
how Dante's area is counted (they capture 64–78% of N221204's rise against 39–51%; the
pair alike). They also win on the level in the Dante convention on all three shots, on
the late drive, on Fig. 10, on the hot-spot transfer and on HYDRA's LEH share. The old
energetics win on the level only if the lips are taken to glow like the wall (open
area), and on the capsule share. That share is HYDRA's, a model, and `energetics=new`
calibrated α_c on a drive that was then too hot, so it is a recalibration issue
(ROADMAP "Absolute energetics on N221204"), not evidence for the old physics. **Decision**: `examples/n221204.jl` now
defaults to `energetics=new`; `energetics=old` restores the 2026-10-24 settings. Untagged
runs and CSVs before 2026-10-25 used the old energetics. No src default changed; the
pieces stay opt-in (`du_wall_loss`, `nif_conversion`, `CapsuleAlbedo` are
HYBRID-E-specific or calibrated, and the 0.75 materials default anchors every other
baseline). *Superseded the same day* ("Energetics defaults and the capsule share
(2026-10-25)"): `energetics=old` is removed, drive-following `gold()`/`uranium()` walls
default to `nif_conversion()`, the DU factor is 0.836 and the example's α_c is 0. *And
again* ("Capsule absorption, conversion and the ablating-wall release (2026-10-25)"):
`gold()`/`uranium()` convert at 0.75 and the drive-history walls are `nif_gold()`/
`nif_uranium()`; α_c is Young's albedo on the K24 schedule; `release` is 1.0.

Open, with their data: the late-time drive (§2: the pair is the test); the pair
residual of 2–4 eV (pointing and cone fraction of N191118 are unpublished); the capsule
share (α_c recalibration on the current drive, or a larger absorbing radius); the
velocity check (needs per-quad delivered powers or window geometry).

## Energetics defaults and the capsule share (2026-10-25) *calibrated*

ROADMAP "Absolute energetics on N221204". The user decided that the published ("new")
energetics become the defaults wherever they are the physics Briareus uses, the old paths go
without deprecation, and the constants calibrated on the old, too reflective wall
(constant-T loss at the time since the pulse start, before "N221204 late-time drive") are
recalibrated. Dante 1, PRE 109 Fig. 10 and HYDRA's LEH share are checks here, never fit.
Confidence: **H** primary text, **L** our construction.

### What became a default, and where

| Piece | Where now | Before |
|---|---|---|
| `nif_conversion()` (0.7 t^0.12 ≤ 0.95) | src: the default `conversion` of `gold()` and `uranium()` whose albedo follows the drive (any `AlbedoModel` but a constant) | 0.75 for every `gold()` / `uranium()` |
| static walls (`gold()`, `gold(albedo = 0.85)`) | unchanged, 0.75 (VisRad-style static solves; at `solve`'s default `time = 0` the rising conversion would read 0.53) | — |
| Au-lined DU wall | the material, `uranium()` (`du_wall_loss()`, now with the rising conversion); `examples/n221204.jl` builds HYBRID-E's wall with it | `gold(albedo = du_wall_loss())` behind `wall=du` |
| DU loss factor | `DU_LOSS_FACTOR` 0.785 → **0.836** (below) | calibrated at 13 ns, constant T |
| wall history, re-emission | already `drive_history`'s defaults; the example's `history=` / `release=` switches are gone. `wall_history = false` stays in src for comparing a march with standalone `solve`s (tests, the `diagnostic_windows` tutorial) | — |
| `CapsuleAlbedo` on HYDRA's ablation-front schedule, 11% degradation, measured backscatter | `examples/n221204.jl` (shot-specific inputs); its `energetics=`, `wall=`, `conversion=`, `backscatter=` options and `capsule=ablator` are removed | `energetics=new` |
| `CapsuleAlbedo` default α_c 0.2 | unchanged (Swadling 2025's published f_cap, not calibrated on the old drive) | — |

*Superseded in part (2026-10-25, "Capsule absorption, conversion and the ablating-wall
release")*: `gold()`/`uranium()` now convert at a fixed 0.75 and the drive-history walls
are the presets `nif_gold()`/`nif_uranium()`; the example's α_c is Young's published
albedo on N221204's K24 schedule, not the fitted 0; `release` is 1.0.

### DU loss factor at the wall's effective time (*calibrated*)

The 6.5% peak-power equivalence (Döppner 2015) was fitted on the high-foot hohlraum
with gold's loss at 13 ns, the time since the pulse start (α_Au 0.944): the
constant-T evaluation that "N221204 late-time drive" showed to leave the wall too
reflective. `drive_history` now uses the wall's effective age, t_eff = ∫(T/T_now)^5.32
dt. From N221204's measured Dante 1 (Pak 2024 Fig. 1; the foot stores almost nothing)
it is 0.76 / 1.04 / 1.64 / 2.04 ns at 5.5 / 6 / 7 / 8 ns (*L*). The high-foot pulse's
peak (~3–4 ns at 390 TW, ~1.7 MJ; the 90 eV foot negligible) gives t_eff ≈ 2–4 ns at
Dante's peak; 2.5 ns is taken (*L*). Same setup as before (h 0.4, 390 TW × 0.86,
conversion 0.9, two-band wall source, capsule albedo 0.3; session scratch `du_cal.jl`):

| t (ns) | f = 0.785: equivalence, ΔT at equal power | α_Au | f for 6.5% | ΔT at that f |
|---|---|---|---|---|
| 13 (old) | 6.51%, +4.9 eV | 0.944 | 0.785 | +4.9 |
| 2.0 | 9.16%, +6.6 | 0.880 | 0.838 | +4.7 |
| 2.5 | 8.95%, +6.5 | 0.891 | **0.836** | +4.8 |
| 3.0 | 8.74%, +6.4 | 0.899 | 0.833 | +4.8 |
| 4.0 | 8.36%, +6.2 | 0.910 | 0.827 | +4.8 |

- **f = 0.836, range 0.827–0.838** over t_eff 2–4 ns (the t_eff uncertainty dominates;
  the mesh moves the ratio < 0.3%). Check, not fitted: Döppner's T-1 pair, +4.8 eV at
  equal power against the measured +5 ± 6.
- With a less reflective wall the wall loss is a larger share of the source, so a given
  loss reduction is worth more power and the fit needs less of it. The new factor
  agrees with the independent cocktail estimates, Rosen's "about 0.83" and Schein's
  0.825. The earlier reading ("the calibrated f below Rosen's cocktail value says the
  Lindl gold loss is low at late times") was an artifact of the constant-T evaluation.
- `test/energetics.jl` "wall material factor" checks the 6.5% at 2.5 ns.

### Capsule albedo recalibrated on the current drive (*calibrated*, at its bound)

The target is the original fit's: E_cap = 16.2% of the coupled energy (PRE 109 App. 3,
read as a fraction of the degraded laser energy), in the full model of
`examples/n221204.jl` (calibrated CBET, wall motion, lips, bubble; `uranium()` with
f = 0.836, wall history and re-emission, 11% degradation, measured backscatter,
`CapsuleAlbedo` on HYDRA's ablation-front schedule). Coarse mesh 0.8 / 0.4:

| α_c | capsule share | LEH share | capsule Tr 8.0 ns | Dante 1 initial / open, 8.0 ns | Fig. 10 min, zero, 8 ns, rms |
|---|---|---|---|---|---|
| 0.22 (previous default; with f 0.785: 11.7%) | 11.4% | 24.8% | 307.8 | 292.1 / 315.6 | −0.185, 6.24, +0.077, 0.042 |
| 0 | 13.9% | 23.9% | 305.2 | 289.9 / 313.1 | −0.183, 6.21, +0.085, 0.045 |

- The share is linear in α_c to the accuracy needed (−11.4 points per unit α_c), so
  16.2% would need **α_c = −0.20 ± 0.03** (the range from the 9–13% degradation and
  the mesh, *L*): unphysical. **The bounded fit is α_c = 0**, with 2.3 points of the
  share unexplained; the example's default is now 0. (*Superseded*: the 2.3 points came
  from the misread capsule anchor, "The P2 swing" side result; the example now uses
  Young's published albedo on the K24 schedule, 13.8%.)
- So the missing capsule absorption is not an albedo. As found before ("Wall
  re-emission and the capsule radius"), HYDRA's share needs absorption outside the
  ablation front (the W-doped blow-off) or a larger effective radius than the front's.
  No published radius of that absorbing surface was found, so nothing else is fitted.
  α_c = 0 contradicts Young 2020's HYDRA ablator albedo (0.16 at 6 ns, 0.3–0.4 later;
  `capsule_albedo=young`, 12.3% before this change): the calibrated value stands in
  for absorption the schedule leaves out; it is not a material property.
- The two recalibrations barely interact: at α_c 0.22 the share moved 11.7 → 11.4%
  with f 0.785 → 0.836.

### N221204 with the new defaults (standard mesh 0.28 / 0.16, full model; checks)

| t (ns) | 4.5 | 5.0 | 5.5 | 6.0 | 6.5 | 7.0 | 7.5 | 8.0 | 8.5 | 9.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| capsule Tr | 212.8 | 256.4 | 274.7 | 280.0 | 283.9 | 291.8 | 299.4 | 305.1 | 294.5 | 272.8 |
| Dante 1, initial area | 209.6 | 244.5 | 259.7 | 266.1 | 272.3 | 278.7 | 285.2 | 290.2 | 279.0 | 256.8 |
| Dante 1, open area | 219.4 | 257.8 | 275.7 | 284.2 | 292.3 | 300.4 | 308.0 | 313.5 | 301.4 | 277.4 |
| measured Dante 1 | 203.4 | 237.9 | 263.1 | 275.5 | 282.9 | 291.3 | 301.7 | 310.1 | 311.3 | 301.5 |

| check | before (α_c 0.22, f 0.785; "Wall re-emission", std) | now | HYDRA / data |
|---|---|---|---|
| wall / LEH / capsule (% of coupled) | 62.7 / 25.6 / 11.7 | 62.2 / 23.9 / 13.9 | 58.6 with gas / 25.2 / 16.2 |
| capsule Tr peak | 310.5 at 8.0 ns | 305.1 at 8.0 ns | ~300 internal; drive Trad 313 |
| Dante 1 at 8.0 ns, initial / open | 295.0 / 318.6 | 290.2 / 313.5 | 310.1 |
| Dante 1 level 5–8 ns, initial / open area, mean | −4.8 / +14.8 eV (coarse) | −9.4 / +9.9 eV | — |
| Dante rise 5 → 8 ns, initial / open | +46.2 / +56.3 | +45.7 / +55.7 | +72.2 |
| Fig. 10 min, upward zero, 8 ns, rms (blind) | −0.189 (5.25), 6.23, +0.077, 0.042 | −0.190 (5.25), 6.20, +0.089, 0.046 | −0.170 (5.6), 6.67, +0.095 |
| wall re-emitted | 103 kJ | 101 kJ | — |

- The weaker DU factor takes ~3 eV off the drive and the black capsule ~2.5 eV more:
  the capsule now peaks at 305 eV (HYDRA ~300), and Dante 1 at 8 ns brackets the
  measured 310.1 (290 / 314). Over the rise and peak the initial-area reading is 9 eV
  low on average and the open-area one 10 eV high; the measured history still lies
  between the two conventions. The shape (the rise through the peak) is unchanged, and
  so is the deficit after the droop.
- The LEH share falls to 23.9% (HYDRA 25.2): the extra capsule sink takes from the
  LEH as well as from the wall, so an albedo alone cannot match both shares.
- Fig. 10: the late P2 rises by 0.012 at 8 ns (more capsule sink, less smoothing), the
  minimum is unchanged; rms 0.042 → 0.046.
- Not rerun: Zylstra's pair (its difference was independent of the energetics, within
  1 eV between old and new; its level should drop by about N221204's ~5 eV) and the
  implosion metrics.

### Other constants checked

- **Wall re-emission `release` 0.77**: calibrated against 1-D diffusion with the
  constant-T scaling of the same medium, independent of the drive. Unchanged. (Later
  2026-10-25: 1.0, the ablating wall's value; "Capsule absorption, conversion and the
  ablating-wall release".)
- **`LipClosure`** (speed, depth, hold): fitted on a gold basecamp history with
  `nif_conversion()` (gold walls are unaffected by the DU factor); the speed scales as
  the lip temperature's square root, so the few eV these changes move a drive would move
  it ~1%, inside its fitted range (0.010–0.030 mm/ns). Not refitted.
- **CBET (`n221204_cbet`, `hybride_cbet`), `BubbleAbsorption`, `WallMotion`**: fitted to
  cone fractions, HYDRA absorbed fractions and radii, not to the drive temperature.
  Unchanged.
- **M-band group albedo α_M** (0.55 with thermal spots, `examples/mband_dewald2020.jl`):
  calibrated on a static solve with gold's loss at 5.5 ns since the pulse start and
  η 0.75. With the wall at its effective age (0.8 ns, mid-plateau after a rise at
  ~4.4 ns, *L*) and `nif_conversion` at 5.5 ns (0.86), cone fraction 1/3, F_M at α_M
  0.55 / 0.6 is 1.90 / 2.16 (was 2.06 / 2.33): F_M ≈ 2 needs α_M ≈ 0.57. The move is
  below the scan's 0.05 step and F_M ≈ 2 is itself approximate, so α_M stays 0.55
  (negative result). The example now uses the consistent wall; its Dante 1 gives
  277 eV against Dewald's 290 eV peak (292 before, where the too reflective wall and the
  low conversion compensated), and the share above 1.8 keV at α_M 0.55 is 0.133
  (Dewald 0.12; 0.114 before).
- **`docs/literate/leh_power_balance.jl`**: its static wall albedo was gold's loss at
  6 ns since the start (0.925); it is now `du_wall_loss()` at the effective age from
  N221204's Dante (1.07 ns), 0.875. The predicted N210307 → N210808 power saving is
  8.3% at that albedo (measured 9.4%; 9.5% at 0.90, 11% at the old 0.925); the
  Monte Carlo agreement is unchanged (2×10⁻⁴ of the source).

### What remains

- The capsule share is 13.9% against 16.2% at the albedo's bound; the absorbing surface
  outside the ablation front is the open input (a published HYDRA absorption radius, or
  the blow-off's opacity). (Now 13.8% with Young's albedo on the K24 schedule; the
  remaining 2.4–4.2 points are `CapsuleAlbedo`'s `outside_absorptance`, unset: "Capsule
  absorption, conversion and the ablating-wall release" §2.)
- The late-time drive after the droop (unchanged by these defaults).

## Late-time drive: the wall, the plasma and Dante (2026-10-25) *derived*

ROADMAP "Absolute energetics on N221204", the late-time drive ("Zylstra 2022 LEH pair, time-dependent" §2): after the
laser falls the model's drive falls far faster than Dante on all three shots. Literature
in `reference/litreview/late_time_drive.md`. Dante, Fig. 10 and HYDRA's Tr are checks
here; nothing was fitted to them. Confidence: **H** primary text, **H-img** read from a
rendered figure, **L** our construction.

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| HYDRA internal Tr (1.1–1.5 capsule radii): max / at peak velocity / at bang | N210207 295 / 279 / 260 eV (coast 1.05 ns); N210307 294 / 270 / 254 (1.2 ns); N210808 305 / 291 / 276 (0.86 ns, 3.1 mm LEH) | Kritcher et al., PRE 106, 025201 (2022) Table I | H-img |
| Dante impulse response | "~400 ps wide Gaussian", cable and XRD; cross-timing ±50 ps | Rubery et al., PRL 132, 065104 (2024) | H |
| | cable removed by Wiener deconvolution; XRD response ~140 ps FWHM | Rubery et al., RSI 94, 031101 (2023) (OSTI 1975922) | H |
| Dante Tr convention | initial LEH dimensions, "a lower bound" | Rubery 2024; Swadling et al. arXiv:2501.10350 eq. 12 | H |
| LASNEX after laser-off | "falls more quickly after the lasers turn off" than Dante; simulated LEH closes too fast after the drive | Swadling et al. (BAH), arXiv:2501.10350 | H |
| Corona energy | P_corona = (1 − η)(P_L − P_bs) "retained in the hot, under-dense coronal plasma"; 15–30% of the energy on NIF (5–10% Nova); η 85–90% | Kline et al., PRL 106, 085003 (2011) | H |
| | "up to 40%" stored in the corona at NIF scale (XSN); after N091204's laser peak the XSN, f = 0.05 model follows Dante's fall, DCA f = 0.15 falls ~20–25% lower in flux by 1.5–2 ns | Town et al., PoP 18, 056302 (2011) (OSTI 1018446) Fig. 2 | H, H-img |
| Cooling term in 0-D models | T^3.3 ∝ exp[−3.3 (t − t_laser)/τ_decay] (τ_decay not printed); 0-D fails in the cooling phase | Callahan et al., PoP 27, 072704 (2020) eq. 9, Sec. II; after Moody et al., PoP 24, 042709 (2017) (closed) | H |
| Cooling rate | set by the laser fall, LEH size, wall area and opacity; drive-pressure cooling "typical[ly] 1.5 ns" | Lindl, Haan & Landen, PoP 30, 012705 (2023) (OSTI 1961756) | H |
| Ablating gold wall | P = r e/V, r = 0.25; E/A = 0.58 T^3.34 t^0.59 hJ/mm² (k = 0) | Rosen 2005 (OSTI 885380) eqs. 1.26–1.30 | H |
| Reheating after bang | W-doped detritus heated by the exploding core, Dante peak ~500 ps after bang; 220 ± 30 kJ on N221204 | Rubery 2024 | H |

### 1. Dante's time response: ruled out

The model's Dante 1 history (standard mesh, `metrics=on` runs of "Zylstra 2022 LEH pair")
convolved in flux with a Gaussian: 140 ps FWHM moves it ≤ 1 eV anywhere in 5–9 ns, 400 ps
≤ 5 eV (N210601 +5.0 eV at 9.0 ns, N191118 +3.5 at 8.5 ns, N221204 ≤ 1). Against 30–70 eV
this is not the cause.

### 2. The wall: an independent 0-D hohlraum (`examples/late_time_wall.jl`) *L*

A 0-D hohlraum on each shot's measured pulse, ηP = A_w F_wall(T) + A_open σT⁴ solved
every step (η 0.75; HYBRID-E 6.4 × 11.24 mm; capsule a 0.7-absorbing sphere at R₀, or
none with the LEHs at 0.857 r₀), with three walls: the static 1-D diffusion reference of
`wall_reemission.jl`; a new ablating reference (1-D Lagrangian radiation hydrodynamics,
Hammer–Rosen power laws, P = 0.25 e/V, flux-limited diffusion, Marshak boundary; energy
conserved to 1e-3; at constant 300 eV it stores E = 0.0039 T^3.48 t^0.67 MJ/cm², 9%
below Lindl's 5.2e-3 T^3.3 t^0.62 at 1–2 ns and 2.4× the static wall); and
`drive_history`'s closure (Lindl's gold loss, `release` 0.77). Each column falls from the
peak of that column; Dante is shown for the shape only.

T (eV) after the laser falls, capsule at R₀ (static / ablating / closure / Dante):

| t (ns) | N191118 | N210601 | N221204 |
|---|---|---|---|
| 7.0 | 298 / 273 / 266 / 293 | 299 / 269 / 263 / 289 | 307 / 278 / 271 / 291 |
| 8.0 | 251 / 239 / 227 / 281 | 290 / 265 / 257 / 301 | 311 / 282 / 275 / 310 |
| 8.5 | 191 / 198 / 180 / 262 | 261 / 243 / 232 / 292 | 297 / 272 / 264 / 311 |
| 9.0 | 171 / 177 / 160 / 238 | 214 / 213 / 195 / 278 | 270 / 253 / 242 / 302 |
| 9.5 | 158 / 163 / 145 / 221 | 187 / 190 / 172 / 254 | 222 / 221 / 204 / — |

- **No gold wall holds the drive.** From 7.0 to 9.0 ns on N191118 the ablating wall falls
  35%, the static 43%, Dante 19%; on N210601 from 7.5 to 9.5 ns 30% against 15%. With no
  capsule and the LEHs at 0.857 r₀ (the most favorable sinks) N191118's ablating wall still
  falls 30% (298 → 210 eV) by 9.0 ns. The ablating wall stores 620–820 kJ at the fall but
  gives back only 30–45 kJ in the next ns: most of its energy is kinetic or deep. So the
  slow fall needs a source outside the wall, as HYDRA's internal Tr (−29 to −40 eV from
  maximum to bang, against −94 to −127 eV for Briareus's capsule Tr on the pair) indicates
  it has.
- **The closure is good to ~10 eV.** It sits 2–3% below the ablating reference at the
  peak and 8–10% below it 1–1.5 ns after the fall (N191118 160 against 177 eV at 9 ns).

**The release factor on an ablating wall** (`late_time_wall.jl release`): the closure
given the constant-T scaling of the ablating wall needs `release` 1.07 / 0.99 / 0.91 for
10 / 30 / 50% drops (static wall 0.76–0.78): the cooled surface throttles the release
less when the heated layer is expanded. In the 0-D hohlraum itself (closure with the
ablating scaling against the ablating reference, T normalized at 7 ns), the fall to
8.5–9.5 ns is matched by 0.92–1.0 on N191118, 0.97–1.05 on N210601 and ~1.06 on N221204;
0.77 gives 6–15% too low a T 1.5–2 ns after the fall (*L*). *Calibrated* candidate:
`release` ≈ 1.0 (0.9–1.07). The present 0.77 was calibrated on a static wall while
`WallLossAlbedo`'s loss is Lindl's ablative one, so 1.0 is the consistent value.
**Not adopted** (a default change for the user to decide; see §3). *Adopted later the
same day: `release` = 1.0, "Capsule absorption, conversion and the ablating-wall
release" §1.*

### 3. The full model with `release` 1.0 (sensitivity)

`examples/n221204.jl 0.8 0.4 full <transfer> x calibrated <shot> metrics=on` on the new
defaults (commit 2e20e7a), `WallLossAlbedo`'s `release` 0.77 (default) against 1.0 (set
through a scratch copy of the source). Capsule Tr / Dante 1 initial / open area (eV),
coarse mesh; measured Dante in the last column:

| shot, t (ns) | release 0.77 | release 1.0 | measured |
|---|---|---|---|
| N191118 8.0 | 249.0 / 236.8 / 255.8 | 256.4 / 243.5 / 263.0 | 281.4 |
| N191118 8.5 | 199.4 / 186.0 / 201.0 | 220.6 / 205.9 / 222.4 | 261.9 |
| N191118 9.0 | 176.7 / 164.9 / 178.1 | 205.6 / 191.8 / 207.2 | 238.2 |
| N191118 bang 9.37 | 161.5 / 150.6 / 162.7 | 196.4 / 183.2 / 197.9 | 226.0 |
| N210601 8.5 | 259.1 / 243.9 / 263.4 | 265.7 / 249.8 / 269.9 | 292.3 |
| N210601 9.0 | 222.3 / 207.1 / 223.7 | 238.9 / 222.5 / 240.4 | 277.7 |
| N210601 bang 9.48 | 196.5 / 183.0 / 197.6 | 221.2 / 206.0 / 222.5 | 253.8 |
| N221204 8.5 | 294.6 / 278.8 / 301.1 | 296.7 / 280.7 / 303.2 | 311.3 |
| N221204 9.0 | 272.8 / 256.7 / 277.2 | 279.6 / 262.8 / 283.9 | 301.5 |
| N221204 bang 9.53 | 230.6 / 214.7 / 232.0 | 248.9 / 231.8 / 250.4 | — |

- Results are unchanged until the power falls (rising drives are unchanged by construction).
- The Dante deficit (initial area) at bang falls from −75 to −43 eV (N191118) and −71 to
  −48 eV (N210601); at 8.5 ns from −76 to −56 and −48 to −43; N221204 at 9.0 ns from −45
  to −39. So ~30–45% of the late gap at bang, 6–26% at 8.5 ns. The rest is the
  plasma (§4).
- Fig. 10 (blind) unchanged to 8.5 ns (−0.183 at 5.25, 6.21 ns, +0.085, rms 0.045); mean
  flux P2 from rise to bang −0.0081 against −0.0076 (−6.1 against −6.5 µm). Split (N221204,
  % of coupled): LEH 24.2 → 24.7, capsule 13.9 → 14.0; wall re-emitted 136 → 145 kJ
  (N191118: 270 → 295 kJ).

### 4. What remains, with its data

- **A plasma source after the fall.** The published evidence (Kline 2011's definition of
  1 − η as energy retained in the corona, 15–30% of the energy on NIF; Town 2011 Fig. 2:
  the model with more coronal energy follows Dante's fall) points to the laser-heated
  gold corona/bubble, with the capsule's W-doped exhaust (~100 kJ thermal; it is also
  where HYDRA's extra 2 points of capsule absorption likely go) and the gas (30–60 kJ) as
  smaller reservoirs. Size of order the gap (1 − η = 10–25% of 1.6–2 MJ); its return rate
  after laser-off is unpublished. Only the thermal part (~60 kJ, ~0.3 ns) is estimable,
  and that is the 0.15 ns lag already tested. **Not modeled**: no published number
  constrains the return without fitting Dante.
- **Data that would settle it**: HYDRA's energy inventory through the coast (radiation,
  wall internal and kinetic, corona/bubble, capsule exhaust, gas) at laser-off and at
  bang; a post-shot synthetic Dante 1 for N191118 or N210601; Dante spectra (color
  against brightness temperature, the >1.8 keV share) between laser-off and bang;
  Moody 2017's cooling constant (closed access); gated LEH/bubble images after
  laser-off.
- **Decision needed** (ROADMAP "Absolute energetics on N221204"): whether `release` 1.0 (the ablating-wall
  calibration) becomes the default. Not changed here. *Decided: 1.0 is the default
  ("Capsule absorption, conversion and the ablating-wall release").*

## The P2 swing: capsule radius, backscatter and the collision (2026-10-25) *derived*

> **Superseded in part (2026-10-26).** K24 Fig. 1's radius is not the ice–gas surface:
> it leaves R(0) at ~3.6 ns, ~2 ns before that surface can move, so it follows the
> ablator–fuel interface or the fuel's center of mass, and the lower bound on the front
> at 6 ns is 0.84 R₀, not 0.79 ("Capsule radius schedules: published histories and the
> API (2026-10-26)"). The runs and their numbers stand.

ROADMAP "The P2 swing" ("what remains"). The N221204 flux P2/P0 swing turns positive ~0.45 ns before
PRE 109 Fig. 10 and its minimum comes ~0.3 ns early. Each published input that bears
on it is checked here for how far it moves the upward zero crossing and the minimum.
Fig. 10 is the blind check throughout: nothing below is fitted to it, and no default
changed. Runs: `examples/n221204.jl 0.8 0.4 full <options>` (coarse mesh, full model,
current energetics defaults; 2–4 min each; coarse reproduces the standard mesh within
~0.01 ns in the crossing and 0.006 in P2). Confidence: **H** primary text, **L** our
construction.

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| N221204 capsule radius history (HYDRA) | R/R(0) flat to ~3.6 ns, 0.943 at 5.5, 0.906 at 6.0, 0.767 at 7.0, 0.537 at 8.0, 0.220 at 9.0, minimum 0.055 at 9.54 ns (bang 9.53); the N210808 capsule on the same drive reaches its minimum at 9.18 ns | Kritcher et al., PoP 31, 070502 (2024), doi:10.1063/5.0210904 (CC BY), Fig. 1 bottom right, digitized (`kritcher2024_capsule()`, `data/kritcher2024_fig1.csv`) | H-img |
| Which radius | the fuel's inner (ice–gas) surface: it holds until the first shock breaks out (keyhole N230110, the N221204 design: 3.86–3.91 ns, Ruof 2026 Table II), and its 0.906 at 6 ns matches the ice–gas interface drawn at 6 ns, 0.92 mm against 0.985 initially (63.0 px at 68.75 px/mm in K24 Fig. 1; PRE 109 Fig. 1a, the same picture, gave 0.92 mm before) | K24 Fig. 1 caption ("the position of the ice-gas interface (red)"); PRE 109, 025204 Fig. 1a caption (ablator light gray, DT ice blue) | H (caption), L (identification) |
| HYDRA's treatment of backscatter | not modeled explicitly, "but as part of an effective total laser power multiplier and δn/n" | PRE 109, 025204 §II | H |
| HYDRA's late-time symmetry | set with Callahan's data-based model and measured sensitivity curves; "artificial multipliers on the input laser power ... to match the in-flight symmetry"; late P2 "not currently predictive" (LEH closure, plasma near the bubble, CBET), HYBRID-E needing a 1.065 peak-CF multiplier and, with 4.0 mm LEHs, symmetrized drive from 7.3 ns | PRE 109 §II; Kritcher et al., PoP 28, 072706 (2021) §III | H |
| N221204 swing calibration | "a calibration experiment for N221204 had not yet been completed"; the swing level "is currently unknown" and was reduced by hand for N230729 | K24 §II | H |
| Inner-beam glint | < 8 TW (< 2% of peak) escapes by specular reflection, at 0.3–0.6 mg/cc | Lemos et al., PoP 29, 092705 (2022), doi:10.1063/5.0099937 (abstract) | H |

**The earlier capsule anchor was misread.** The research notes since "Wall material,
conversion and capsule energetics" call 0.92 mm (0.81 R₀) at 6 ns "HYDRA's ablation
front". It is the ice–gas interface; the ablation front lies outside it by the shell's
thickness (DT and remaining HDC, unpublished at 6 ns). Bounds on the front at 6 ns:
0.79 R₀ (zero thickness, 0.868 × 0.906) to 0.906 R₀ (the shell keeping its relative
thickness, i.e. K24's scaled curve applied to R₀). The example's Young-shaped schedule
(0.81 R₀ at 6 ns) sits near the lower bound; `schedule=k24` (new) is the upper.

**No published HYBRID-E bubble measurement was found** (searches 2026-10-25: thin-wall
windows, gated LEH bubble images, x-ray imaging of the hohlraum interior, for HYBRID-E,
N210808, N221204, Au-lined DU 6.4 mm). The thin-wall and G-LEH bubble data are all HDC
SymCap basecamp (Izumi 2018/2021, Meezan 2020, Chen 2020) or foam-lined / I-Raum
platforms. The HYBRID-E bubble speed stays inferred from the picket energy (`picket_ratio`
1.1) and PRE 109 Fig. 1a's N210808 tip. Reading the outer picket from Pak's delivered
pulse (`n221204_pulse()` split by DN24's cone fraction) gives 43.9 kJ in the first ns
after the turn-on (the example's "~40" was by eye), i.e. `picket_ratio` ≈ 1.69 × 0.73 = 1.23
against the basecamp's ~26 kJ (by eye): a faster bubble, the wrong direction for Fig. 10.

### Effect of each mechanism (coarse; crossing, minimum)

| Run | upward zero (ns) | minimum (at ns) | at 8 ns | capsule share | rms vs Fig. 10 |
|---|---|---|---|---|---|
| baseline (`capsule=model`, Young-shaped schedule) | 6.21 | −0.183 (5.25) | +0.085 | 13.9% | 0.045 |
| `capsule=mesh` (shrinking mesh, same schedule) | 6.22 (+0.01) | −0.171 (5.25) | +0.056 | 14.0% | 0.047 |
| `capsule=mesh schedule=k24` | 6.18 (−0.03) | −0.176 (5.25) | +0.062 | 16.0% | 0.050 |
| `schedule=k24` (model) | 6.17 (−0.04) | −0.179 (5.25) | +0.091 | 16.0% | 0.049 |
| `backscatter=uniform` (then `hydra`) (48 kJ as a uniform multiplier) | 6.26 (+0.05) | −0.196 (5.25) | +0.084 | 14.0% | 0.042 |
| `backscatter=off` | 6.27 (+0.06) | −0.199 (5.25) | +0.084 | 14.1% | 0.042 |
| transfer `fig9b` (HYDRA's post-CBET cone fraction) | 6.32 (+0.11) | −0.203 (5.00) | +0.080 | — | 0.044 |
| `fig9b` + `backscatter=uniform` (then `hydra`) (uniform cut) | 6.37 (+0.16) | −0.203 (5.00) | +0.079 | — | 0.043 |
| `collision=6.39` (t_c + 0.3 ns) | 6.38 (+0.17) | −0.183 (5.25) | +0.085 | 13.9% | 0.043 |
| `collision=5.79` (t_c − 0.3 ns) | 6.10 (−0.11) | −0.183 (5.25) | +0.085 | 13.9% | 0.052 |
| `bubble_energy=0.8` (t_c 6.38) | 6.54 (+0.33) | −0.205 (5.25) | +0.077 | 14.1% | 0.040 |
| `backscatter=uniform capsule=mesh schedule=k24` | 6.23 (+0.02) | −0.189 (5.25) | +0.062 | 16.2% | 0.047 |
| PRE 109 Fig. 10 (HYDRA) | 6.67 | −0.170 (5.60) | +0.095 | 16.2% (App. 3) | — |

1. **Capsule radius** (published in K24 Fig. 1): ≤ 0.04 ns either way, and the K24
   schedule moves the crossing *earlier* (a larger capsule sees more of the late
   pole-hot pattern). Minimum within 0.012. It does not close the gap. The old
   schedule's time map is confirmed: K24's minimum falls at bang time (9.54 ns), as
   the map assumed.
2. **Backscatter as HYDRA has it**: the measured inner-cone SBS (16.8% of the inner
   power over 5.5–6.5 ns, the shape *L*) pulls the post-CBET cone fraction 0.06 below
   Fig. 9b at 5.5–6.5 ns (0.59 / 0.56 / 0.53 at 5.5 / 6.0 / 6.5 ns against 0.65 / 0.62 /
   0.59). Fig. 10's HYDRA books backscatter as a total-power multiplier, so
   `backscatter=uniform` (then `hydra`) is the matching setting: +0.05 ns, minimum 0.013
   deeper. The measured SBS stays the default (the physics); PRE 109 itself notes that
   "including the measured backscatter time history could impact implosion shape".
3. **CBET's time dependence**: forcing the HYDRA post-transfer cone fraction (`fig9b`,
   read by eye) adds +0.06 more (6.37 with the uniform backscatter cut) but moves the
   minimum *earlier* (5.00 ns) and makes the rise too waist-hot (−0.128 at 4.5 ns against
   −0.092). The calibrated two-zone CBET follows Fig. 9b to rms 0.016 when fed HYDRA's
   inputs; the residual here is the input cone fraction (DN24 delivered against HYDRA
   requested) and the backscatter. Kritcher 2021: CBET shifts the swing, not its slope.
4. **Ablator–wall collision**: the crossing follows t_c at ~0.5 ns per ns (the ridge
   stop turns on at 6.09 ± 0.2 ns, just before the 6.2 ns crossing). Closing the gap
   by t_c alone would need t_c ≈ 6.9 ns; PRE 109 Fig. 1a already shows wall and ablator
   plasma meeting at z 2–3 mm at 6 ns (by eye, *L*), so a collision near 6.1 ns is the
   published picture. Not changed.
5. **Bubble speed**: `picket_ratio` 0.8 (from only the N210808 tip of PRE 109 Fig. 1a) gives
   +0.33 ns, the largest single lever, as found before; Pak's pulse points the other way
   (1.23). No HYBRID-E datum separates them.
6. **Glint** (Lemos 2022, estimated, not run): < 8 TW out of the inner beams, ~1.8 points of
   cone fraction at peak, −0.005 in P2 at the crossing (d(P2/P0)/dCF −0.29%/%), i.e. ≤ 0.03
   ns *earlier*. LEH lips (computed closure): ≤ 0.005 in P2 (before). The gas fill enters
   through the bubble speed and opacity scalings (0.3 mg/cc on both platforms).

So the published inputs move the crossing by −0.04 to +0.16 ns (HYDRA's CBET result and the
uniform backscatter cut together); ~0.3 ns remains, and only the two least constrained bubble inputs
(picket fluence, collision) can carry it. Fig. 10 is itself HYDRA with power and
cone-fraction multipliers set to in-flight symmetry data of other shots (PRE 109 §II),
on a design whose swing calibration had not been done (K24): the 0.4 ns is a difference
between two models whose late-time inner-beam propagation is calibrated on different
data, not a measured discrepancy. The minimum's timing (5.25 against 5.6 ns) moves with
none of these except the Fig. 9b transfer, which moves it the wrong way.

**Side result: the capsule share.** With K24's radius as the absorbing surface (the
upper bound), the capsule takes 16.0% of the coupled energy at α_c = 0 (16.2% with
the uniform backscatter cut): HYDRA's 16.2% (App. 3). The 2.3-point "missing capsule
absorption" of "Energetics defaults and the capsule share" comes from the misread
anchor: the Young-shaped schedule shrinks the absorbing radius too early (0.81 against
≤ 0.906 R₀ at 6 ns, and from 2 ns instead of the shock breakout at ~3.6 ns). Capsule Tr
drops 2–3.5 eV through the peak (276.5 against 280.0 at 6 ns, 303.3 against 305.2 at
8 ns) and the LEH share falls 23.9 → 23.1%. The front lies between the bounds, and with
the published HDC albedo (Young 2020, 0.16 at 6 ns) the k24 schedule would give ~14%.
Not adopted: the schedule and α_c = 0 belong to item 16 (decision below).

### Figures of merit beside the swing (derived, not fits)

| Metric | N221204 | Meaning |
|---|---|---|
| Callahan 2018 x = √(E_pk,out/(A_out ρ)) τ/R · r_cap/R (A 53 mm², ρ 0.3 mg/cc, R 3.2, r_cap 1.135) | 1.61–2.11 (picket 43.9 kJ in the first ns or 65.4 kJ to the trough at 2.1 ns; τ 8.75 ns half-power fall or 9.40 ns end) | his fitted line (P2_hs − 2.7(CF − 33) = 101.9 − 85.0x µm) gives −35 to −78 µm at 33% CF without CBET: +13 to +29 points of cone fraction needed from transfer for a round hot spot |
| Izumi 2021 FoM ∝ (E^0.2 t′/R)^2.5 (ρ/0.3)^0.13 (r_cap/R)^1.1, relative to the HDC basecamp (t′ 5 ns, R 2.875, r_cap 0.908) | 1.62 (E 1.1, t′ 6.33 ns: outer turn-on 0.4 to mid-peak 6.73 ns); 1.39–1.90 over E 0.8–1.5; 1.72 at Pak's E 1.23 | P2 = 90(3CF − FoM) (his Eq. 2, basecamp FoM 1): round at CF_eff ≈ FoM/3 = 0.54 (0.46–0.63) |
| post-CBET cone fraction, peak | 0.53–0.59 (Briareus, calibrated two-zone); 0.59–0.65 (Fig. 9b) | +20–26 (+26–32) points over 33% |
| measured hot-spot P2 | 0.3 ± 0.1 µm (P24) | round |

Both figures of merit, extrapolated beyond their data (Izumi's seven shots span
0.8–1.4 of the basecamp; Callahan's line was fitted without CBET), put N221204 where
the transfer HYDRA computes (and Briareus is calibrated to) makes the hot spot round:
consistent, and blind to the swing's timing, since both are time-integrated. They do
not discriminate between the crossing at 6.2 and at 6.7 ns.

### What would close it

- A HYBRID-E (6.4 mm Au-lined DU, 0.3 mg/cc) thin-wall or G-LEH measurement of the
  inner-beam spots or the bubble tip against time (Izumi 2018's platform on this
  hohlraum): it sets the bubble speed and the collision time directly.
- The N221204-design 2DConA in-flight shell P2 against radius (the calibration K24 reports as
  pending; Ruof 2026's N240304 surrogate gives only 270–160 µm, at 2.9 Å):
  integrated constraints on the swing, through a hydro transfer.
- The HYDRA cone-fraction and power multipliers behind Fig. 10 (unpublished), and
  whether HYDRA's flux is evaluated at the ablation front or a fixed radius.
- The shell thickness (or ablation-front radius) at 5–7 ns for N221204: it fixes the
  capsule schedule between the bounds, and the capsule share with it.

### Decisions needed

*Decided 2026-10-25 ("Capsule absorption, conversion and the ablating-wall release",
"What changed"): `schedule=k24` is the example's default (with Young's albedo); Fig. 10
is compared with a uniform backscatter cut (`backscatter=uniform`) beside the measured
default; `capsule=mesh` stays off.*

- Whether the example's capsule schedule becomes K24's (`schedule=k24`, the published
  N221204 history, an upper bound on the front) or a front between the bounds; with it,
  α_c = 0 is no longer at its bound (item 16).
- Whether `backscatter=uniform` (then `hydra`) (the HYDRA bookkeeping) is the setting for
  the Fig. 10 comparison (the measured SBS shape stays *L* either way).
- Whether to shrink the capsule mesh by default (`capsule=mesh`): crossing +0.01 ns,
  minimum 0.012 shallower, late P2 −0.03; the right choice depends on where Fig. 10's
  flux is evaluated, which is not stated.

## Capsule absorption, conversion and the ablating-wall release (2026-10-25) *derived*

> **Superseded in part (2026-10-26).** The bounds of §2 changed. The 0.81 R₀ of the
> `young` schedule at 6 ns lies below the lower bound, now 0.84 R₀ (K24's radius is the
> ablator–fuel interface or fuel center of mass, not the ice–gas surface). The schedules
> are now `RadiusSchedule`s with exactly the same values ("Capsule radius schedules:
> published histories and the API (2026-10-26)"). The numbers stand.

ROADMAP "Absolute energetics on N221204" (and 9): the user's four decisions of 2026-10-25 on "Energetics defaults and
the capsule share", "Late-time drive" and "The P2 swing", implemented together. Dante 1,
PRE 109 Fig. 10 and HYDRA's 16.2% capsule share are checks here; nothing is fitted to
them. Runs: `examples/n221204.jl 0.8 0.4 full <transfer> x calibrated <shot> metrics=on`
(coarse mesh, full model) on the final defaults below. Confidence: **H** primary text,
**L** our construction.

### What changed

| Piece | Now | Before |
|---|---|---|
| capsule albedo α_c (example) | Young 2020's quasi-1D HYDRA ablator albedo on the Young time map (`capsule_albedo=young`, the default; a number still accepted) | 0, fitted to HYDRA's share and stopped at its bound |
| absorption outside the front | `CapsuleAlbedo(; outside_absorptance = a_out)`: α_eff = 1 − (1 − α_c)s² − a_out(1 − s²), s = R/R₀; default 0, so the default is the front's alone (`capsule_outside=` in the example) | — |
| capsule radius schedule (example) | `schedule=k24`, the N221204 HYDRA history (Kritcher 2024 Fig. 1), the default; Zylstra's pair: K24's N210808 capsule, time × bang/9.18 (*L*) | Young's N170827 front through the misread 6 ns anchor (`schedule=young` kept it; removed 2026-10-26) |
| conversion of `gold()` / `uranium()` | 0.75 whatever the albedo | `nif_conversion()` when the albedo was an `AlbedoModel` |
| NIF drive-history walls | `nif_gold()` (`WallLossAlbedo()`, `nif_conversion()`), `nif_uranium()` (`du_wall_loss()`, `nif_conversion()`); the example uses `nif_uranium()`, `test/energetics.jl` `nif_gold()` | `gold(albedo = WallLossAlbedo())`, `uranium()` |
| `WallLossAlbedo` `release` | **1.0** (the ablating wall; §1) | 0.77 (static wall) |
| example options | `pulse=dn24` (N221204/N220919) and `dynamics=old` removed; the DN24 table stays as N210808's pulse and as every shot's cone-fraction split; `backscatter=hydra` renamed `backscatter=uniform`; `capsule=mesh` stays off | — |
| Fig. 10 comparison | with the measured backscatter (default) the full model is run again with a uniform 48 kJ power cut ("+bubble, uniform cut"), and Fig. 10 is compared with that run | one run |

**Backscatter wording.** A uniform power cut is one simple treatment, the one PRE 109's
HYDRA used ("as part of an effective total laser power multiplier", §II); modeled
backscatter is per cone and time dependent. Kemp et al., PoP 33, 082705 (2026) (LLNL
news, lmf.llnl.gov, "LLNL researchers improve predictive modeling of laser backscatter")
find that earlier single-species LASNEX/HYDRA inner-cone SBS models overpredicted it in
magnitude and duration, that a multispecies (Au/He/D/T) treatment fixes this, and that a
late-time density spike between capsule and LEH refracts the inner beams and ends the
backscatter early; see also Higginson et al., PRL 134, 165101 (2025) and the outer-cone
SBS paper, PoP 31, 042707 (2024). These are cited only; a per-cone, time-dependent backscatter
model from them is a separate task. The measured per-cone SBS/SRS (shapes *L*) stays the
physics default.

### 1. Wall re-emission: the ablating wall against Rosen's self-similar solution

`examples/late_time_wall.jl rosen` drives the 1-D Lagrangian ablating gold wall with its
surface temperature held (the solution's boundary condition) and pure diffusion,
against Hammer & Rosen's ablative self-similar solution as quoted by Rosen (OSTI 885380,
Lecture 1 after eq. 1.30; **H**): m_F = (9.9, 7.4)×10⁻⁴ T^1.91 t^0.52 g/cm² and E/A =
(0.58, 0.39) T^3.34 t^0.59 hJ/mm² (T in heV, t in ns) for k = 0 and the ramp
T = T₀ t^0.18. E is internal plus kinetic; m_F is read where T = T_B/2, divided by the
zero-order profile's 15/16 there. The solution is first order in ε = 0.29 (energy
conservation verified to ε², ~8%).

| case | t (ns) | E/A, 1-D over H&R | m_F, 1-D over H&R |
|---|---|---|---|
| T = 2 heV | 0.5 / 1 / 3 | 0.996 / 1.001 / 1.006 | 1.047 / 1.043 / 1.028 |
| T = 3 heV | 0.5 / 1 / 3 | 0.990 / 0.996 / 1.004 | 1.063 / 1.057 / 1.040 |
| T = 2.5 t^0.18 (to 3.05 heV) | 0.5 / 1 / 3 | 0.967 / 0.968 / 0.968 | 1.065 / 1.059 / 1.038 |

- **It agrees within the accuracy of the solution**: stored energy to 1–3%, heated
  mass 3–7% high, against ε² ≈ 8%. Halving the zone masses moves E ≤ 0.3% and
  m_F ≤ 2.5% (`fine=on`).
- With the flux limiter the reference otherwise uses (`limiter=on`) E is 9–12% lower at
  2–3 ns (0.71–0.77 at 0.25 ns): the limiter, absent from H&R, throttles the low-density
  blow-off. The fit E = 0.0039 T^3.48 t^0.67 MJ/cm² of "Late-time drive" refers to the
  incident T_D (Marshak boundary), not the surface T, hence its difference from H&R's.
- **The release** (`late_time_wall.jl release [limiter=off]`, the closure given the
  constant-T scaling of that wall, drops of 10 / 30 / 50%): 1.07 / 0.99 / 0.91 flux
  limited (as before), 0.98 / 0.90 / 0.84 in pure diffusion; the 0-D hohlraums gave
  0.92–1.06 (flux limited). **Adopted: `release` = 1.0** (*calibrated*), range
  0.84–1.07 (drop size and limiter). The static wall's 0.76–0.78 (`wall_reemission.jl`)
  no longer applies: Lindl's loss, which `WallLossAlbedo` uses, is the ablative one.

### 2. The capsule share with Young's albedo on both radius bounds

N221204, full model, final defaults (release 1.0, `nif_uranium()`), coarse; Dante 1 over
the initial / open LEH area:

| schedule (front at 6 ns) | capsule share | LEH share | capsule Tr 6.0 / 8.0 / 8.5 | Dante 8.0 | Dante 8.5 | Dante 9.0 | Dante at bang 9.53 | Fig. 10 min, zero, 8 ns, rms: measured backscatter / uniform cut |
|---|---|---|---|---|---|---|---|---|
| `k24` (0.906 R₀, upper bound; **default**) | **13.8%** (14.0 with the uniform cut) | 24.7% | 278.8 / 306.3 / 298.4 | 290.9 / 314.2 | 282.2 / 304.8 | 264.6 / 285.8 | 233.6 / 252.4 | −0.181 (5.25), 6.19, +0.080, 0.047 / −0.195 (5.25), 6.24, +0.080, 0.044 |
| `young` (0.81 R₀, near the lower bound 0.79) | 12.0% (12.1) | 25.4% | 282.0 / 307.9 / 299.4 | 292.3 / 315.7 | 283.0 / 305.7 | 265.0 / 286.3 | 234.0 / 252.8 | −0.184 (5.25), 6.23, +0.076, 0.044 / −0.197 (5.25), 6.32, +0.075, 0.041 |
| before (α_c 0, `young`, release 0.77, `uranium()`) | 13.9% | 23.9% | — / 305.2 / — | 289.9 / 313.1 | — | — | — | −0.183 (5.25), 6.21, +0.085, 0.045 |
| HYDRA / measured / Fig. 10 | 16.2% | 25.2% | ~300 internal | 310.1 | 311.3 | 301.5 | — | −0.170 (5.6), 6.67, +0.095 |

- With the published ablator albedo, absorption at the front gives 12.0–13.8% between
  the radius bounds, against HYDRA's 16.2%: 2.4–4.2 points remain. That is the share
  `outside_absorptance` would carry (a_out ≈ 0.1–0.2 by the share's near-linear response, *L*); it is
  not set, since no blow-off opacity or absorbing radius is published. Alternatively
  Young's N170827 albedo is too high for N221204's thicker W-doped ablator. Neither is
  fixed here.
- **Default schedule `k24`, on physical grounds**: it is the published N221204 HYDRA
  history with the right timing (the radius holds until the first shock breaks out,
  ~3.6–3.9 ns, and reaches its minimum at bang); `young` is another shot's front shape,
  anchored on a misidentified surface (the ice–gas interface) and shrinking from 2 ns.
  Its magnitude is the upper bound (a shell that keeps its relative thickness, whereas
  ablation and compression thin it), so its capsule share is an upper bound for
  absorption at the front. HYDRA's 16.2% played no part in the choice (both bounds miss
  it).
- The schedule moves the drive by 1–3 eV through the peak (the larger capsule is the
  larger sink) and the swing by ≤ 0.04 ns. The uniform cut moves the crossing +0.05 to
  +0.09 ns and deepens the minimum by 0.013–0.014, as before. Mean flux P2/P0 from the
  rise to bang with the uniform cut: −0.0144 against Fig. 10's −0.0147 (−0.3 µm in
  hot-spot P2 through the −9.2 µm/% transfer); with the measured backscatter −0.0082.

### 3. Late-time Dante on the three shots (final defaults)

Dante 1 over the initial LEH area (eV), coarse. The "before" columns are "Late-time
drive" §3 (α_c 0, `young` schedule, the same wall otherwise):

| shot, t (ns) | now | release 1.0 before | release 0.77 before | measured |
|---|---|---|---|---|
| N191118 8.0 | 244.1 | 243.5 | 236.8 | 281.4 |
| N191118 8.5 | 206.8 | 205.9 | 186.0 | 261.9 |
| N191118 9.0 | 192.8 | 191.8 | 164.9 | 238.2 |
| N191118 bang 9.37 | 185.5 | 183.2 | 150.6 | 226.0 |
| N210601 8.5 | 251.0 | 249.8 | 243.9 | 292.3 |
| N210601 9.0 | 223.9 | 222.5 | 207.1 | 277.7 |
| N210601 bang 9.48 | 207.3 | 206.0 | 183.0 | 253.8 |
| N221204 8.5 | 282.2 | 280.7 | 278.8 | 311.3 |
| N221204 9.0 | 264.6 | 262.8 | 256.7 | 301.5 |
| N221204 bang 9.53 | 233.6 | 231.8 | 214.7 | — |

- The deficit at bang is now −40.5 eV (N191118) and −46.5 eV (N210601); at 9.0 ns −45 /
  −54 / −37 eV. The release takes ~30–45% of the late gap, as in "Late-time drive" §3;
  the capsule changes (Young's albedo, the K24-shaped schedules) add 1–2 eV. The rest is
  the plasma source after the fall, still not modeled.
- The capsule shares on the pair are 14.9% (N191118) and 14.7% (N210601); the wall
  re-emitted 290 and 192 kJ (N221204 143 kJ).
- Tests: `test/energetics.jl` (presets, `outside`), `dynamics.jl`, `retune.jl`,
  `capsule_drive.jl`, `docrefs.jl` pass on the new defaults.

### Decisions needed

*Decided 2026-10-26 (user): `outside` stays 0 until a blow-off opacity or absorbing
radius is published (a candidate for the 13.8 vs 16.2% capsule share); the pair's
K24-shaped schedule is accepted, labeled derived, with room in the API for other
schedules (research pending, ROADMAP "Absolute energetics on N221204").*

- Whether to set `outside_absorptance` (a_out ≈ 0.1–0.2 would close HYDRA's share on the two
  bounds) without a published blow-off opacity or absorbing radius; it is 0.
- The pair's K24-shaped schedule (N210808 capsule, time scaled to bang) is *L*; the radius
  histories of the pair are unpublished.

## Capsule radius schedules: published histories and the API (2026-10-26) *derived*

ROADMAP "Absolute energetics on N221204", the user's decision of 2026-10-26 on "Capsule absorption, conversion and the
ablating-wall release": the pair's K24-shaped schedule is accepted, labeled derived,
with room in the API for other schedules. What published radius histories exist for the
HYBRID-E shots in use (N221204, N210808, the Zylstra pair N191118 / N210601), which
surface each follows, and what the new `RadiusSchedule` types support. No default or
result changed. Confidence: **H** primary text or table, **H-img** digitized figure,
**M** known from citing papers only, **L** our construction.

### Published radius histories and anchors

| Source | Shot(s) | Which radius | Times | Uncertainty | Use | Conf. |
|---|---|---|---|---|---|---|
| Kritcher et al., PoP 31, 070502 (2024), doi:10.1063/5.0210904, Fig. 1 (`kritcher2024_capsule()`) | N221204 capsule, and the thinner N210808 capsule, on N221204's 2.05 MJ drive (HYDRA) | not stated; the ablator–fuel interface or the fuel's center of mass (below), as R/R(0) | 0–9.5 ns (N210808 capsule to its minimum, 0.055 at 9.18 ns) | ±0.01 in R/R(0) (0.02 at 9.25 ns), ±0.04 ns | **schedule**: `kritcher2024_radius()`, the default `schedule=k24`; `:n210808` time-scaled for the pair | H-img (curve), L (surface) |
| Young et al., PoP 27, 082702 (2020) Fig. 3a (`young2020_radii()`) | N170827 (LLNL dynamic view-factor model) | the model's ablation front | 0–7.75 ns | ±0.03 mm | schedule on a piecewise time map (`schedule=young`, removed 2026-10-26); another shot's shape | H-img |
| Kritcher et al., PRE 106, 025201 (2022) Table I (HYDRA post-shot) | N210207 / N210307 / N210808 | limb-minimum radius at peak velocity (RPV) | point anchors | — | `TrajectoryRadius` anchors, checks: v_imp 393 / 387 / 391 km/s at 8.63 / 8.9 / 8.8 ns, RPV 195 / 211 / 185 µm, bang 9.05 / 9.35 / 9.26 ns (measured 9.09 / 9.31 / 9.26 ± 0.03), so peak velocity comes 0.42–0.46 ns before bang | H |
| Kritcher et al., PRE 109, 025204 (2024) Table I (as read in `litreview/n221204_benchmark.md`; OSTI 2474067 has no full text) | N210808 / N221204 | — | point anchors | — | anchors: v_imp 401 / 383 km/s (HYDRA, recalibrated on the surrogates; PRE 106 had 391 for N210808), ablator remaining 3.8 / 5.7%, bang 9.28 / 9.53 ± 0.07 ns measured, 9.27 / 9.61 simulated; no time or radius of peak velocity | H |
| Ruof et al., PoP 33, 052706 (2026) (OSTI 3367476, LLNL-JRNL-2019015) Table II, keyholes | N221106 (N210808 design), N230110 (N221204 design), N231016 | first shock's breakout from the ablator into the liquid D2 (VISAR) | 3.6–6.1 ns | ±0.006–0.06 ns | the time the ablator–fuel interface starts to move: equator / pole 3.714 / 3.650 ns (N221106), 3.864 ± 0.012 / 3.906 ± 0.018 (N230110; HYDRA 3.899 / 3.886), 3.575 / 3.588 (N231016); shock 1 at 32.7–34.0 km/s, merged 1–2 at ~50 km/s from 4.85 ns (N230110) | H |
| Ruof et al. 2026 Table IV and Figs. 8–9, 2DConA | N220618, N221010 (N210808 design), N240304 (N221204 design, 85.9 µm, D2 gas, lower peak power) | limb minimum and outer slope of the backlit shell (inner and outer boundaries), Legendre P0 | from R < 300 µm to stagnation (the last ~0.6–0.8 ns) | figure; HYDRA trajectory offsets −50 / −35 / +25 µm | **check only**: shell velocity 392 ± 24 / 381 ± 15 / 355 ± 5 km/s (HYDRA 401 / 408 / 344), stagnation 9.170 ± 0.061 / 9.230 ± 0.029 / 9.876 ± 0.090 ns, mass remaining 8.8 ± 1.0% (N240304; HYDRA 10.2). A surrogate without ice on the pulse of each shot, and R/R₀ < 0.27 (R² < 0.07, where the capsule's absorption hardly matters) | H (the extracted table's rows realigned to their labels by us) |
| Olson et al., SAND2011-8633C (OSTI 1141159) | NIC CH ConA | center of mass of the remaining ablator (ConA's radius) | ~19–22 ns, 900 → ~200 µm | — | context: baseline HYDRA over-predicted the ConA velocity by ~18%, matched with 85–92% peak-power multipliers; another platform | H (slides) |
| Hicks et al., PoP 19, 122702 (2012) (1-D streaked ConA) | NIC CH | remaining-mass center and shell thickness | — | — | not read; another platform | M |
| Callahan et al., PoP 27, 072704 (2020) §II.B (LLNL-JRNL-809797) | scaling model | shell radius of a rocket model, m dv/dt = −4πr² P_abl, dm/dt = −4πr² ṁ, driven by their T_r(t) fit | — | — | an analytic r(t), but it needs an ablation rate and pressure fitted per design and the drive; not adopted (below) | H |
| Zylstra et al., PRE 106, 025202 (2022) | N191118, N210601 | none published | — | — | bang times 9.37 / 9.48 ns, which place the pair's scaled schedule; no v_imp or RPV | H |

**Which radius K24 Fig. 1 follows.** The curve stays at 1.000 to 3.50 ns and is 0.995 at
3.75 (the N210808 capsule 0.999 at 3.50): it starts to move at ~3.6 ns, 0.2–0.3 ns before
the N221204 design's keyhole breakout (3.86–3.91 ns), as the N210808 capsule (~3.4 ns)
does before its keyhole breakout (3.65–3.71 ns), a common offset between HYDRA on N221204's
drive and the keyholes' pulses plus reading. The ice–gas surface cannot move until the first
shock has crossed the 65 µm of ice: at the keyhole's speeds (32.7 km/s to the 1–2 merge at
4.85 ns, ~32 µm, then ~50 km/s) that is ~5.5 ns (*L*). So the curve follows the
ablator–fuel interface (1.050 mm initially) or the fuel's center of mass, which start
moving when the shock enters the fuel. The 6 ns picture agrees (*L*): the interface at
0.906 × 1.050 = 0.951 mm leaves 31 µm of shocked ice (65 µm compressed 2.1×) to the
ice–gas surface drawn at 0.92 mm. Before today the curve was read as the ice–gas surface
("The P2 swing"), from that 0.92 mm alone. Consequences:
- The ablation front at 6 ns lies between **0.84 R₀** (no ablator left outside the
  interface, 0.906 × 1.050/1.135) and 0.906 R₀ (K24's scaled curve on the outer radius, the
  shell keeping its relative thickness): narrower than 0.79–0.906. `schedule=k24` stays
  the upper bound, and the `young` schedule's 0.81 R₀ now lies below the lower bound.
- The data file's header, `kritcher2024_capsule`'s docstring and the example's comments
  say so; the digitized numbers are unchanged.

**What a schedule can be built from.** Only K24 Fig. 1 is a radius history of the shots
in use, and only for N221204's drive (the N210808 capsule on it, not the N210808
shot itself). Young 2020's front belongs to N170827. The 2DConA trajectories (Ruof 2026)
cover the last ~0.7 ns of surrogates, where R² < 0.07, and are figures. The rest are point
anchors (breakout, peak velocity and its time and radius, bang).

### The API (`src/target.jl`, `src/dynamics.jl`)

- `RadiusSchedule`: abstract; `radius_fraction(s, t)` gives R/R₀. `CapsuleAlbedo(;
  radius)` and `moving_hohlraum(; capsule_shrink)` take a schedule, a `(times, R/R₀)`
  table (now stored as a `TabulatedRadius`) or a function, as before.
- `TabulatedRadius(t, r; source)`: linear, held at its ends, times strictly increasing;
  `source` records where the table comes from. `kritcher2024_radius(:n221204 | :n210808)` returns K24's
  two curves (the N210808 one closed at 0.055, 9.18 ns), labeled published.
- `ScaledRadius(base; time_scale, radius_scale, source)`: R = radius_scale ·
  base(t/time_scale); for a tabulated base the times are scaled once and interpolated as
  a table, which is what the example did, so the pair's schedule
  (`ScaledRadius(kritcher2024_radius(:n210808); time_scale = bang/9.18)`) is the same
  table, exactly. Its `source` records it as derived.
- `TrajectoryRadius(; start, peak_time, peak_radius, velocity, bang, minimum, radius)`:
  a kinematic interpolant through published anchors, for a shot with anchors but no
  history: R(0) to `start`, a cubic Hermite segment to `peak_radius` at `peak_time` with
  speed `velocity` (at rest at `start`), a second to `minimum` at `bang` (at rest), held
  after. Both segments must be monotone (speed ≤ 3× the mean, Fritsch–Carlson), else an
  `ArgumentError`. *Derived*. Checked against K24 Fig. 1, which it is not built from:

  | anchors | max \|ΔR/R₀\| (3 ns to peak velocity) | max \|Δ(R/R₀)²\| | at 5 / 6 / 7 / 8 ns (K24) |
  |---|---|---|---|
  | N221204: start 3.864 (keyhole), bang 9.53, 383 km/s (PRE 109), peak velocity 0.44 ns before bang at 197 µm (the PRE 106 rows' means), R/1.050 mm | 0.019 (7.0 ns) | 0.029 | 0.972 / 0.892 / 0.748 / 0.528 (0.965 / 0.906 / 0.767 / 0.537) |
  | N210808 capsule: start 3.714, 391 km/s at 8.8 ns, 185 µm (on its drive, PRE 106) | 0.028 (6.75 ns) | 0.043 | 0.959 / 0.861 / 0.693 / 0.442 (0.959 / 0.881 / 0.720 / 0.452; on N221204's drive) |

  N221204's time and radius of peak velocity are borrowed from the N210808 family, so the
  first row is a check of the form, not of N221204's anchors. The trajectory runs inside
  K24 by up to 0.02 from 5.5 to 8 ns. Not used by the example; for the pair it would need
  v_imp and RPV, which are unpublished, so the pair keeps `ScaledRadius`.
- **Left out: a rocket model.** Callahan 2020's rocket model needs an ablation rate and
  pressure fitted per design and the drive (a coupled march with `drive_history`, since
  the capsule's absorption depends on R). CompactLES stage 0 on our full-model drive, at
  the effective rate fitted to the remaining mass (ṁ = 2.85e5 T³), stagnates at 10.2 ns
  against bang 9.53 ("Capsule drive hand-over and CompactLES stage 0"), and the LLNL
  HYDRA times the surrogates with peak-power multipliers of 0.825–0.86, varying 5% shot
  to shot (Ruof 2026; Olson 2011 needed 0.85–0.92 on NIC): a rocket trajectory from our drive would be a 0.5–0.7 ns
  timing error, not a schedule. A `RadiusSchedule` subtype can take one later, e.g.
  from CompactLES's trajectory as a `TabulatedRadius`.

### Check of identical results

`examples/n221204.jl` now builds its schedules from the types (`kritcher2024_radius()`,
`ScaledRadius(...)` for the pair, `TabulatedRadius(...)` for `young`). The capsule
albedo and R/R₀ of all five configurations (k24; young at bang 9.53 and 9.37; pair at
9.37 and 9.48), each with `outside` 0 and 0.15, at 992 times from −1 to 12 ns (a 13.7 ps
grid plus every knot and the bang times), and the shrinking capsule mesh (area and a
centroid, k24 and the N191118 pair schedule, five times) were written from the old code
(0c1bd2a) and the new with `repr` (shortest round-trip): the 9930 lines are identical
(`cmp`). `test/energetics.jl` ("capsule radius schedules") checks the k24 and pair cases
with `==` against the old tables. The example was not rerun (2–4 min coarse per run on
an idle machine; the schedule enters only through these values).

### Open

- The front's position inside 0.84–0.906 R₀ (the shell thickness at 5–7 ns) is
  unpublished; it sets the capsule share between the bounds.
- `schedule=young` is now below the lower bound; removed 2026-10-26.
- The radius histories, v_imp and RPV of the pair are unpublished.

## Material defaults and their sources (2026-11-18)

ROADMAP 17, third bullet: a library of material defaults, each quantity with its data
source and uncertainty. `material_library()`, `material_default(name)`,
`MaterialDefault`, `SourcedValue`, `MATERIAL_BASES` (`src/material_library.jl`,
`test/material_library.jl`, reference page "Material defaults"). Nothing in the package
or the GUI changes: no default, constant or constructor is touched. Confidence as at the
top of this file (**H** primary text, **M** abstract or secondary, **L** inference);
sections cited by title hold the earlier sources and are not repeated. Literature
searched this session: the Z dependence of laser-to-x-ray conversion (Nakano & Kuroda,
PRA 27, 2168 (1983), abstract: periodic in Z, no numbers; the 0.53 µm and 0.26 µm
Z scans of PRA 33, 525 (1986) and PRA 34, 4184 (1986), abstracts only: undulatory in Z,
no numbers read) and low-Z (CH) conversion at 351 nm (no number found).

### Survey

"Today" is the value the package uses and where it is set; "inside" states whether it lies
within the published spread.

| Material | Quantity | Today (where) | Published (source, conf.) | Uncertainty / spread | Inside? |
|---|---|---|---|---|---|
| Au wall, static | x-ray albedo | 0.8 (`gold()`, `src/target.jl`; `hohlraum()` wall) | 0.8, the static VisRad model of N170827 (Young et al., PoP 27, 082702 (2020); "Young 2020 validation (2026-10-24)") **H** | 0.77 (Swadling et al., PoP 32, 052707 (2025), VisRad fitted to Dante at 400 TW, **H**) to 0.84 (same paper's 0-D model at 288 eV, **H-img**); Lindl's loss at a NIF peak's effective wall age 0.88–0.94 ("Energetics defaults and the capsule share (2026-10-25)", *derived*) | yes (equal to Young's) |
| Au wall, static | conversion | 0.75 (`gold()`, `uranium()`) | ≈0.7 at 1 ns on Nova, up to 0.85 for shaped pulses (Lindl et al., PoP 11, 339 (2004) §IV); 0.68 at 1 ns, Omega (Rosen, OSTI 885380 (2005)); "Energetics (`energetics.md`)" **H** | 1 ns hohlraums 0.68–0.85; NIF: 0.85–0.90 (Olson et al., PoP 19, 053301 (2012), **M**), 0.9 (Young 2020 static model, **H**), 0.95 (Swadling 2025, **H**) | yes for 1 ns pulses; below every NIF value |
| Au wall, drive history | albedo model | `WallLossAlbedo()` (`nif_gold()`) | Lindl 2004 eq. 4-31 ("Gold wall loss and thin-foil burnthrough (2026-10-16)") **H-img** | Hammer & Rosen and Shussman & Heizler losses within ~10% over 1–3 heV, 1–10 ns (*derived*); `release` 1.0, range 0.84–1.07 (*calibrated*) | yes |
| Au wall, drive history | conversion | `nif_conversion()`, 0.7 t^0.12 ≤ 0.95 | "Wall material, conversion and capsule energetics (2026-10-24)" §2 (*derived*) | 0.80 / 0.85 / 0.90 at 3 / 5 / 8 ns against Olson's 0.85–0.90 (**M**); exponents 0.12–0.2 published | yes at NIF peak |
| Au-lined DU wall | albedo model | `du_wall_loss()` (`uranium()`, `nif_uranium()`) | factor 0.836 at 300 eV *calibrated* to Döppner et al., PRL 115, 055001 (2015) ("Energetics defaults and the capsule share (2026-10-25)") | 0.827–0.838 over effective ages 2–4 ns; Rosen's 0.83 and Schein's 0.825 (**H**) | yes |
| Au-lined DU wall | conversion | 0.75 / `nif_conversion()` | as the Au wall (the liner converts) | as the Au wall | as the Au wall |
| Au, GUI hohlraum part | albedo / conversion | 0.85 / 0.8 (`BriareusGUI` `HOHLRAUM_VALUES`; "as the former hohlraum template had") | no separate source | albedo inside 0.77–0.94, conversion inside 0.68–0.95 | yes, but unsourced |
| Au, directly irradiated (foil, disk) | conversion | 0.05 (`foil`, `backlighter_foil`: `gold(conversion = 0.05)`); GUI foil starter 0.02 | 0.95 at 10¹⁴, 0.80 at 10¹⁵ W/cm², backscatter-corrected plateau after ~0.5 ns on Au, U and cocktail spheres (Dewald et al., PoP 15, 072706 (2008), abstract) **M** | open 3ω disks 0.57 / 0.59 / 0.51 of the absorbed energy at 2e13 / 1e14 / 4e14 W/cm² (Mead et al., LA-UR-86-300 (1985) Fig. 16, **H**; "Unconverted light: absorption and conversion vs. wavelength (2026-10-28)") | **no**: 10–19× below the published 0.51–0.95 (25–48× for 0.02) |
| Au | 3ω laser absorption | 1 (every surface; backscatter booked separately) | ≈0.93 / 0.95 / 0.91 at 3e13 / 1.5e14 / 3e14 W/cm², disks (Ze et al., UCRL-86414 (1981) Fig. 4) **H** | read from a figure (two of three marked approximate) | not a material input (CE data are per absorbed or per incident energy as their sources state) |
| Au | 1ω and 2ω, relative to 3ω | absorption 0.88 / 1.0, conversion 0.88 / 0.9 (`nif_unconverted`) | the rule of "Unconverted light: absorption and conversion vs. wavelength (2026-10-28)" (*engineering*, **L**) | 1ω absorption 0.6–1.0, conversion 0.65–1.0; 2ω 0.8–1.0, 0.72–1.16 | yes, by construction |
| Au | 1ω cold reflectance | not used (`specular_check` uses only a 10⁸ W/cm² threshold) | ≈0.98 at normal incidence (Johnson & Christy / Palik) **M** | approximate | — |
| Capsule ablator (HDC) | x-ray albedo | 0.3 (`ablator()`, GUI capsule); `CapsuleAlbedo()` 0.2 | 0.3, Young 2020's static model (**H**); f_cap ≈ 0.8, α ≈ 0.2 (Swadling 2025 0-D, **H**); quasi-1D HYDRA α_c(t) for HDC, 0.16 at 6 ns, 0.3–0.4 later (Young 2020 Fig. 3a, `young2020_albedos()`, **H-img**) | reading error 0.03–0.05; 0 with the fuel-layer area (Lindl 2004) | yes |
| Capsule ablator (HDC, CH, Be) | conversion | 0.3 (`ablator()`, GUI capsule) | **no public source found** | — | unsourced; the benchmark capsule receives no laser light (below) |
| CH, Be | x-ray albedo | as the ablator | no CH- or Be-specific value found beyond the generic ablator values above | — | — |
| Be, Ti | 1ω laser absorption | not used | 0.68 (Be), ≈0.76 (Ti) against 0.82 (Au) at 2.7e13 W/cm² (Ze 1981 Fig. 7) **H** | 1ω only | — |
| Al, Cu | 1ω cold reflectance | not used | ≈0.94 (Al), ≈0.97 (Cu), normal incidence (Johnson & Christy / Palik) **M** | approximate | — |
| Al, Cu, mid-Z foils (Ti, Fe, Zn, Ge, Ag) | broadband conversion (the `Material` input) | 0.05 (gold foil default; no material of their own) | **no public source found** at 3ω (the Z scans above are abstracts without numbers) | — | unsourced |
| Ti, Fe, Zn, Ge, Ag foils | He-α conversion (signal, not radiosity) | `kshell_ce` ("He-α CE: saturating peak (2026-10-26)") | fit to 68 published points (Workman 2001, Barrios 2013/2014, Coppari 2019) *engineering* | ×/÷ 1.43–1.50 inside the data; Ti held out 1.58 ± 0.66% against the fit's 2.97% (−1.1σ); Ag outside the fitted lines (flagged `extrapolated`) | yes inside the data |
| Si | 1ω | not used | nearly transparent before plasma forms **M** | — | — |
| SiO₂ | x-ray mirror reflectivity | `mirror_reflectivity` (Dante channels) | Henke optical constants ("Dante channel responses and the unfold (2026-10-24)") | — | not a laser material |
| black (`BLACK`, apertures, GUI free parts) | albedo / conversion | 0 / 0 (`BLACK`; `Material(name)` keywords; GUI part fields) | a definition (a perfect absorber that emits nothing) | — | by construction |

### The library

`material_library()` holds only entries whose every quantity has a source or a
derivation above. Each quantity is a `SourcedValue` (value, basis, source, uncertainty);
the bases are `:published`, `:range` (a round value inside a published range, not
itself printed anywhere: `gold()`'s conversion 0.75), `:calibrated`, `:derived` and
`:definition`. Entries: `:gold`, `:uranium`, `:nif_gold`, `:nif_uranium`, `:black`
(each equal in every field to the constructor's material today, tested), and the
published alternatives `:gold_young2020` (0.8 / 0.9), `:gold_swadling2025` (0.77 /
0.95) and `:gold_disk` (Lindl's loss; Dewald's 0.95 / 0.80 intensity table; with the 3ω
absorption and 1ω cold reflectance recorded but not used).

*Decision*: unsourced values are left out; they are not included as flagged engineering
defaults. A library entry is a claim that the numbers have a source; a flagged entry
would put unsourced numbers behind the same call and invite their use. The cost is that
`ablator()` has no entry: its albedo 0.3 is Young 2020's, but its conversion 0.3 has no
source, and a `Material` needs both.

### What switching would change (measured)

Static solve on the default-size hohlraum (`hohlraum(length = 10.1, radius = 2.875,
leh_radius = 1.55, capsule_radius = 0.9)`, h 0.8 / capsule_h 0.3, 824 elements; 48 NIF
quads at 8 TW on the GUI's NIF cone pointing, an engineering default), capsule Tr from
a₀ and P2/P0, P4/P0 by `legendre_moments`. Same view factors and deposition for every
row; the capsule receives 0 of the 384 TW. *derived*

| wall / capsule | capsule Tr (eV) | P2/P0 | ΔTr vs core | ΔP2/P0 |
|---|---|---|---|---|
| core `gold()` 0.8 / 0.75, `ablator()` 0.3 / 0.3 | 262.60 | −0.1556 | — | — |
| GUI hohlraum 0.85 / 0.8, capsule 0.3 / 0.3 | 277.77 | −0.1437 | +15.17 (+5.78%) | +0.0119 |
| `:gold_young2020` 0.8 / 0.9, capsule 0.3 / 0.3 | 274.84 | −0.1556 | +12.25 (+4.66%) | 0.0000 |
| `:gold_swadling2025` 0.77 / 0.95, capsule 0.3 / 0.3 | 272.85 | −0.1621 | +10.25 (+3.90%) | −0.0065 |
| core wall, capsule `CapsuleAlbedo()` 0.2 | 261.70 | −0.1542 | −0.89 | +0.0014 |
| core wall, capsule conversion 0 | 262.60 | −0.1556 | 0.00 | 0.0000 |

- Checks: a conversion change alone scales every source by the same factor, so Tr
  scales as its fourth root and P2/P0 is unchanged: (0.9/0.75)^¼ × 262.60 = 274.84
  (row 3, exact). On h 0.5 / 0.2 (1900 elements) every ΔTr agrees within 0.01 eV and
  every ΔP2/P0 within 2e-4.
- The albedo sets P2: 0.85 smooths the drive (+0.012 in P2/P0), 0.77 sharpens it.
- The foil's x-ray source is its conversion times its absorbed laser power, its only
  source, so a foil conversion change scales its emission exactly: 0.05 → 0.80
  (Dewald at 10¹⁵) is ×16 in power and ×2.0 in the Tr it drives on any surface; 0.02 →
  0.80 is ×40 and ×2.5.

Not changed; for the user (recommendations):
1. **GUI hohlraum 0.85 / 0.8**: switch to `:gold_young2020` (the published NIF static
   setting: Tr −2.9 eV, P2/P0 −0.012 from today's GUI) or to `:gold` to match the core
   (−15.2 eV). Recommended: `:gold_young2020`.
2. **`gold()` conversion 0.75** (sourced as a 1 ns range value, low for NIF): moving it
   to Young's 0.9 raises every static capsule Tr by 4.66% and leaves P2 unchanged; it
   moves every static baseline in the tests and tutorials. Recommended: adopt 0.9 for
   static NIF solves in a deliberate baseline change, or keep 0.75 and point static NIF
   users to `:gold_young2020`.
3. **Capsule conversion 0.3** (`ablator()`, GUI capsule): no source and no effect while
   the capsule is unlit. Recommended: keep, and say in `ablator`'s docstring that it is
   an engineering default (or 0 by convention); keep the albedo 0.3 (Young's static
   value; `CapsuleAlbedo()`'s 0.2 moves Tr −0.9 eV).
4. **Foil conversion 0.05 (core) and 0.02 (GUI)**: both 10–48× below every published
   value for gold. Recommended: `:gold_disk` for gold foils; for mid-Z foils, require
   the conversion from the user (no public total-conversion data were found). The `foil`
   docstring states `material = gold()` while the code uses `gold(conversion = 0.05)`.
5. **GUI free parts 0 / 0**: a black body by definition (`:black`). Recommended: keep,
   with the material picked from the library named in the part's properties.

*Decided the same day* (the user's decisions 1, 2 and 4 above; "Decisions of
2026-11-18" below): `gold()` and `uranium()` convert at 0.9, the GUI hohlraum is
`gold()`, and foils take `foil_material()`.

### Decisions of 2026-11-18: the moved defaults

The user's decisions on recommendations 1, 2 and 4 above; 3 and 5 stay open. Results
change from this commit on (version 0.2.1 unchanged).

1. **`gold()` and `uranium()` convert at 0.9** (was 0.75): Young 2020's static NIF value,
   inside the published NIF 0.85–0.95. The library's `:gold` is now Young's entry (basis
   `:published`); `:gold_young2020` is merged into it (two entries with the same numbers
   and source would drift apart). `uranium()` moves with it: its conversion is defined as
   the gold liner's (0.7 µm Au, Döppner 2015), and the DU factor 0.836 was itself
   calibrated with conversion 0.9 ("Energetics defaults and the capsule share",
   *calibrated*), so the same source applies. The default `efficiency` of `Conversion()`
   0.75 is not changed (a separate default; for the user; **Superseded**: 0.9 since
   "Defaults changed (2026-11-19)" §1). 0.75 remains available as
   `gold(conversion = 0.75)` for 1 ns Nova/Omega hohlraums.
2. **GUI hohlraum part 0.8 / 0.9** (was 0.85 / 0.8): `HOHLRAUM_VALUES` reads `gold()`,
   so the GUI follows the core. Closed diagnostic windows (`WindowSpec`) default to the
   same wall material (they had the old 0.85 / 0.8). Saved sessions keep their values: a
   part's `values` are saved in full since version 8, a window's albedo and conversion
   always, and sessions before version 8 take the template-era 0.85 / 0.8 from
   `_TEMPLATE_DEFAULTS`, which is unchanged.
3. **Foils: `foil_material(element = :Au)`**, the default of `foil` and
   `backlighter_foil` and of the GUI foil starter (was `gold(conversion = 0.05)`; the GUI
   starter 0 / 0.02). A `Material` serving every path needs one number: `flux_response`
   without the laser power, the scenarios' observed flux and `montecarlo_incident` take
   a conversion independent of intensity, and a constant albedo keeps the solve linear,
   so gold takes the static equivalent of `:gold_disk`:
   - conversion 0.8, Dewald's plateau at 10¹⁵ W/cm², the table's end, held beyond it,
     where backlighter foils are driven (He-α optima 1–5 × 10¹⁵ W/cm², "He-α CE:
     saturating peak"); basis `:published`, **M** (abstract);
   - albedo 0.8, Lindl's loss at 1 ns: 1 − 0.32 T^−0.7 = 0.76 / 0.80 / 0.85 at 150 /
     200 / 300 eV (*derived*); basis `:range`.
   The library entry `:gold_foil` equals `foil_material()` in every field (tested).
4. **Mid-Z foils, suggested ranges** (`:titanium_foil` … `:silver_foil`, basis
   `:suggested`, a new member of `MATERIAL_BASES`). The elements are those of
   `kshell_ce` that are solid as pure elements (Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ge,
   Zr, Nb, Mo, Ag; Ga, Br and Kr left out). No 3ω broadband conversion of flat mid-Z
   foils was found (searched again this session: the Z scans at 0.26 µm, PRA 34, 4184
   (1986), at 0.25 µm, Popil et al., PRA 35, 3874 (1987), and at 0.53 µm, Schwanda &
   Eidmann, LPB 9, 551 (1991), give no numbers in their abstracts). What is published,
   and the range it brackets:

   | Source | Conditions | Total conversion | Conf. |
   |---|---|---|---|
   | Eidmann & Kishimoto, APL 49, 377 (1986) | 0.53 µm, 3 ns, 3e13 W/cm², planar, Z 4–82 | 2% (Be) rising to nearly 50% (heaviest), of the incident energy | **M** |
   | Pérez et al., PoP 19, 083101 (2012), Table I | 351 nm, 1 ns, OMEGA, Dante/DMX 0–20 keV | Fe-oxide aerogels 0.35–0.60, steel-lined cavities 0.60–0.75 (Fe K-shell 1.2–3%) | **H** |
   | Mead et al., LA-UR-86-300 (1985) Fig. 16 | 3ω, open gold disks | 0.51–0.59 of the absorbed energy | **H** |
   | Tanabe et al., APL 93, 051505 (2008), for scale | planar Ti foil, K-shell 4–6 keV | (1.4 ± 0.9)% of the incident energy (TiO₂ nanofiber cotton 3.7 ± 0.5%): line output only, under the broadband | **M** |

   - Conversion: suggested range **0.1–0.6**, default its midpoint **0.35**. The upper end
     is gold's open-disk value and the volumetric Fe targets' (a flat foil neither traps
     its radiation as a cavity does nor heats in volume as an aerogel does); the lower
     end lies under 0.13–0.26, the linear-in-Z reading of Eidmann & Kishimoto over
     Z 22–47 at 2ω (3ω converts better), with room for the measured oscillation in Z
     (Nakano & Kuroda 1983; PRA 33, 525 (1986)) (*L*).
   - Albedo: suggested range **0.3–0.8**, default **0.55**: between the static albedos of
     a low-Z ablator and gold in Young 2020 (*L*); it enters only through x-rays incident
     from other surfaces.
   - Midpoints, not an end: a foil's x-ray emission is the wanted source in some designs
     and the unwanted background (preheat of a nearby target, detector background) in
     others, so neither end is conservative; the midpoint bounds the error at half the
     range either way.
   - The broadband conversion is the radiosity source; the He-α line signal is
     `kshell_ce`'s, fitted to separate data and not touched.

**Measured baseline moves** (*derived*; the benchmark above: default-size hohlraum, h
0.8 / capsule_h 0.3, 824 elements, 48 NIF quads at 8 TW on the GUI's cone pointing,
`ablator()` capsule):

| Case | Before | After | Check |
|---|---|---|---|
| capsule Tr, `gold()` | 262.596 eV | 274.842 eV (+12.25, +4.66%) | (0.9/0.75)^¼ × 262.596 = 274.842, difference 6e-14 eV |
| capsule P2/P0, P4/P0, `gold()` | −0.15559, +0.00010 | same (Δ < 2e-16) | a uniform source factor leaves the shape |
| capsule Tr, GUI hohlraum part | 277.768 eV (0.85 / 0.8) | 274.842 eV (−2.93) | the GUI part is `gold()` now |
| capsule P2/P0, GUI hohlraum part | −0.14370 | −0.15559 (−0.0119) | the albedo 0.85 → 0.8 sharpens the drive |
| foil case: mean flux on a 1 mm black sphere 6 mm from a 2 mm `backlighter_foil` driven by four 0.5 TW beams (2.0 TW absorbed) | 2.1961e-4 TW/mm² (21.50 eV) | 3.5138e-3 TW/mm² (42.99 eV) | ×16.000000 = 0.8/0.05 (the foil's only source), Tr ×2.00000 |
| same, the GUI's old foil starter (0 / 0.02) | 8.784e-5 TW/mm² (17.10 eV) | as above | ×40.0000 |
| same, `foil_material(:Ti)` | — | 1.5373e-3 TW/mm² (34.97 eV) | 0.35/0.8 of gold's (the albedo does not enter: the sphere is black) |

In the test suite and tutorials (every number moved was checked by the closed form, not
pasted from a run):

| File | Quantity | Before | After | Check |
|---|---|---|---|---|
| `test/montecarlo.jl` | closed cube, uniform source: H = ηP/A/(1 − α) | 0.75/0.2 | 0.9/0.2 | the run gave 0.2006 above the old closed form, 0.9/0.75 − 1 = 0.2 within its Monte Carlo error |
| `test/windows.jl` | the integrating-cavity η | 0.75 | 0.9 | the open/closed ratio is independent of η (it passed unchanged), but its tolerance `dev`, the solved wall flux's largest deviation from ηs/D, had grown by about 0.2 (the flux is 1.2 times the old closed form's): stale, not failing |
| `test/energetics.jl` | `gold().conversion`, `uranium().conversion` | 0.75 | 0.9 | the decision |
| `BriareusGUI/test/test_hohlraum.jl`, `test_halfraum.jl`, `test_profiles.jl` | the starter's wall against core `hohlraum`; a closed end's cap albedo | `gold(albedo = 0.85, conversion = 0.8)`; 0.85 | `gold()`; the value of the part | the decision (the GUI part is `gold()`) |
| `docs/literate/cad_hohlraum.jl` | "the gold converts 75%"; printed capsule H, six models | 0.75; 0.84761 … 0.84548 TW/mm² | 0.9; 1.01714 … 1.01457 TW/mm² | ×1.2 within the 5-digit printing (linear: constant albedos); a₂/a₀, a₄/a₀ and the loss fractions unchanged in every digit |

Every other test passed unchanged: the hohlraum tests compare with closed forms,
ratios or explicit materials, and the regression numbers of the energetics, N221204,
CBET and M-band calibrations were made with explicit walls (`nif_uranium()`,
`Material(:w; conversion = 0.9)`, `gold(albedo = WallLossAlbedo(), conversion = …)`).
No calibrated constant in these notes was fitted with `gold()`'s 0.75: the DU factor
used conversion 0.9, the capsule albedo and the M-band split the N221204 model's
`nif_uranium()`, the lip closure `nif_conversion()`, the dynamics and CBET fits no
radiation. **No call site is pinned.** The tutorials that set 0.85 / 0.8 explicitly
(`beam_errors`, `dropped_bundle`, `implosion_metrics`) and the examples with explicit
walls are unchanged; `examples/backlighter.jl` now uses `foil_material(:Ti)` in place of
an unsourced `Material(:Ti; conversion = 0.02)`.

## Absolute energetics without new data: the early peak and the missing sinks (2026-11-18) *derived*

ROADMAP item 12, the part that needs no HYDRA energy inventory, post-shot synthetic
Dante or blow-off opacity: the current state with the whole-field Dante 1 (§1), what the
model's existing inputs do to the early peak (§2), and the sinks the model lacks, sized
from published numbers (§3). Measured Dante 1 (Pak et al., PRE 109, 025203 (2024)
Fig. 1, `n221204_pulse()`) and PRE 109 Fig. 10 are checks; nothing below is fitted to
them, and no default or calibrated constant changed. Runs: a scratch copy of
`examples/n221204.jl 0.8 0.4 full calibrated x n221204` (coarse mesh, full model,
measured backscatter, defaults of commit 865f161; 45–55 s per variant on 8 threads;
the uniform-cut rerun skipped). The copy adds only switches the example lacks (the wall
conversion as another `Conversion`, the DU loss factor) and a per-step record: Dante 1
over its whole field and through the hole only, the share of its intensity that is
converted laser light seen directly (`sol.converted`), and one diffuse bounce of each
step's absorbed laser through the step's view factors. The coarse baseline reproduces
"The swing and the bubble without new data" §1 to the printed digit (capsule / LEH share
13.9 / 24.2%, Dante 1 at 8 ns 289.7 / 312.9 eV; Fig. 10 minimum −0.181 at 5.25 ns,
crossing 6.19 ns). Not run: the standard mesh (the coarse drive is within ~0.1 eV of it,
"N221204 late-time drive"), combinations of variants (additivity not measured),
Zylstra's pair. Confidence: **H** primary text, **M** abstract or secondary, **L** ours.

### 1. The current state with the whole-field Dante 1

Coarse, full model; Dante 1 (`dante1()`, the 60 mm field) over the initial LEH area (Dante's
convention, which Rubery et al., PRL 132, 065104 (2024), call "a lower bound") and over the
open area; eV. "Direct share": the fraction of Dante 1's intensity that is converted laser
light, the spots seen directly.

| t (ns) | capsule Tr | Dante 1, initial | Dante 1, open | measured | initial − measured | open − measured | open fraction A_open/A₀ | direct share |
|---|---|---|---|---|---|---|---|---|
| 4.5 | 212.8 | 208.4 | 218.1 | 203.4 | +5.0 | +14.7 | 0.833 | 0.35 |
| 5.0 | 255.1 | 242.5 | 255.7 | 237.9 | +4.6 | +17.8 | 0.809 | 0.31 |
| 5.25 | 267.5 | 252.7 | 267.3 | 252.1 | +0.6 | +15.3 | 0.798 | 0.29 |
| 5.5 | 273.2 | 257.6 | 273.5 | 263.1 | −5.5 | +10.4 | 0.787 | 0.27 |
| 6.0 | 278.8 | 263.9 | 281.9 | 275.5 | −11.6 | +6.4 | 0.768 | 0.25 |
| 6.5 | 282.5 | 269.5 | 289.4 | 282.9 | −13.4 | +6.5 | 0.753 | 0.24 |
| 7.0 | 291.0 | 276.4 | 297.9 | 291.3 | −14.9 | +6.6 | 0.741 | 0.20 |
| 7.5 | 299.5 | 283.6 | 306.3 | 301.7 | −18.1 | +4.6 | 0.735 | 0.18 |
| 8.0 | 306.3 | 289.7 | 312.9 | 310.1 | −20.4 | +2.9 | 0.734 | 0.16 |
| 8.5 | 298.4 | 281.2 | 303.7 | 311.3 | −30.1 | −7.6 | 0.734 | 0.13 |
| 9.0 | 281.6 | 263.7 | 284.9 | 301.5 | −37.8 | −16.6 | 0.734 | 0.09 |
| 9.35 | 265.0 | 247.1 | 266.9 | 293.8 | −46.7 | −26.9 | 0.734 | 0.04 |

Split (% of the coupled 1778 kJ; delivered 2046, backscatter and degradation 268): wall
61.9, LEH 24.2, capsule 13.9 (247 kJ); wall re-emitted 106 kJ; laser out of the holes
−0.1 kJ (none). Dante 1 peaks at 289.7 eV at 8.0 ns (open 312.9) against 311.9 at
8.30 ns measured.

- **Whole field and through the hole agree within 0.1 eV at every step** (4.5–9.35 ns).
  Through the lower LEH at 36.75° to the axis Dante 1 sees only the wall (the capsule
  lies outside the line of sight), and the target's outside is one-sided and dark, so
  the 60 mm field adds nothing on this target. The roadmap's "not rerun on N221204" is
  closed: the numbers are those of the 2026-11-18 swing section, and the move from
  290.9 / 314.2 came from the other commits named there (the M-band split, the LEH
  azimuthal count), not from the view.
- **The area convention is worth 4.7–8%** in Tr: the prescribed lips leave 0.81 of the
  initial area open at 5.0 ns and 0.734 (0.857²) from 7.5 ns, and the open reading is
  the initial one times (A₀/A_open)^¼.
- **The rise through the peak is too small on both readings**: from 5.0 to 8.0 ns +47.1 eV
  (initial area) and +57.2 eV (open area) against +72.2 eV measured. The early excess on
  the open area and the late deficit on the initial area are the same shape error.
- **After the laser falls** (initial area): −30 eV at 8.5 ns, −38 at 9.0, −47 at 9.35 ns,
  as in "Capsule absorption, conversion and the ablating-wall release" §3; the plasma
  source after the fall stays unmodeled ("Late-time drive", §4).

### 2. The early peak: what each existing input does

"13–20 eV above Dante at 5.0–5.5 ns" is an open-area statement. On the Dante
convention the model is +4.6 / +0.6 / −5.5 eV at 5.0 / 5.25 / 5.5 ns, and its rise
over 4.5–5.5 ns lines up with the measured one after a delay of 0.07 ns (rms 2.7 eV);
over the open area the delay is 0.245 ns (rms 3.1 eV) (*L*: the Δt minimizing the rms of
model(t − Δt) − measured on a 0.05 ns grid). Each input below is varied alone within its
documented range; change against the baseline, eV (Dante 1 over the initial area unless
marked):

| Variant (documented range) | 5.0 | 5.25 | 5.5 | 6.5 | 8.0 | 8.5 | 9.0 | open 5.25 / 8.0 | rise 5 → 8 (initial / open) | capsule / LEH (%) | delay Δt initial / open (ns) | Fig. 10 min (5.25), zero, 8 ns, rms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline (absolute values) | 242.5 | 252.7 | 257.6 | 269.5 | 289.7 | 281.2 | 263.7 | 267.3 / 312.9 | 47.1 / 57.2 | 13.9 / 24.2 | 0.070 / 0.245 | −0.181, 6.19, +0.080, 0.047 |
| conversion 0.90 (t/8 ns)^0.2 ≤ 0.95 (Rosen's and Schein's exponent at `nif_conversion()`'s 0.90 at 8 ns) | −3.0 | −2.8 | −2.5 | −1.5 | −0.4 | −0.2 | −0.1 | −2.9 / −0.4 | 49.8 / 60.0 | 13.4 / 23.7 | 0.020 / 0.195 | −0.183, 6.19, +0.081, 0.047 |
| conversion 0.875, constant (Olson's 0.85–0.90) | +2.9 | +2.5 | +2.1 | +0.5 | −1.4 | −1.6 | −1.6 | +2.7 / −1.5 | 42.8 / 52.8 | 14.3 / 24.4 | 0.125 / 0.295 | −0.179, 6.19, +0.079, 0.047 |
| conversion in intensity, 0.95 at 10¹⁴, 0.80 at 10¹⁵ W/cm² (Dewald 2008), no time dependence | +5.0 | +4.0 | +3.3 | +1.0 | −3.9 | −3.3 | −2.2 | +4.3 / −4.2 | 38.2 / 47.7 | 14.6 / 24.4 | 0.165 / 0.340 | −0.176, 6.22, +0.076, 0.044 |
| computed lips (`lips=computed`, `LipClosure`) | −4.5 | −5.2 | −5.2 | −3.0 | −1.5 | −1.3 | −1.6 | −0.1 / +0.6 | 50.1 / 58.0 | 14.0 / 23.5 | 0.020 / 0.240 | −0.175, 6.18, +0.081, 0.047 |
| capsule albedo 0.3, constant (Young 2020's static value; HYDRA's 0.01–0.16 before 6 ns is the default) | +3.0 | +3.2 | +3.3 | +2.6 | +0.9 | +0.7 | +0.7 | +3.4 / +0.9 | 45.0 / 54.9 | 12.1 / 24.9 | 0.115 / 0.285 | −0.184, 6.22, +0.078, 0.043 |
| SBS on N210808's per-cone histories (`backscatter=history`, Kemp 2026) | −3.8 | −3.9 | −2.1 | +1.5 | 0.0 | 0.0 | 0.0 | −4.1 / 0.0 | 51.0 / 61.2 | 13.8 / 24.2 | 0.035 / 0.220 | −0.152, 6.20, +0.080, 0.043 |
| SBS 53 kJ (42 + 11, AS24) | 0.0 | −0.4 | −0.9 | −1.1 | −0.1 | −0.1 | −0.1 | −0.4 / −0.2 | 47.0 / 57.1 | 13.8 / 24.2 | 0.070 / 0.240 | −0.176, 6.18, +0.081, 0.048 |
| DU loss factor 0.92 (×1.1; the loss models agree within ~10%) | −2.3 | −2.5 | −2.6 | −3.0 | −3.9 | −3.7 | −3.2 | −2.6 / −4.2 | 45.5 / 55.4 | 13.3 / 23.2 | 0.040 / 0.215 | −0.187, 6.19, +0.088, 0.047 |
| DU loss factor 0.752 (×0.9) | +2.4 | +2.7 | +2.8 | +3.2 | +4.1 | +3.9 | +3.4 | +2.8 / +4.5 | 48.9 / 59.2 | 14.5 / 25.4 | 0.105 / 0.275 | −0.175, 6.20, +0.073, 0.046 |
| source lag 0.15 ns (`lag=0.15`, a ~60 kJ store; a sensitivity, not a calibrated input) | −9.1 | −6.9 | −4.1 | −1.1 | −0.4 | +4.3 | +5.9 | −7.3 / −0.4 | 55.9 / 66.4 | 13.2 / 24.1 | −0.080 / 0.105 | −0.154, 6.19, +0.082, 0.049 |
| measured (absolute values) | 237.9 | 252.1 | 263.1 | 282.9 | 310.1 | 311.3 | 301.5 | — | 72.2 | HYDRA 16.2 / 25.2 | — | −0.170 (5.6), 6.67, +0.095 |

- **No input within its documented range carries the excess.** The largest single
  reductions at 5.25 ns are the published SBS timing (−3.9 eV), the conversion's time
  exponent at 0.2 (−2.8) and the DU factor at +10% (−2.5, a level shift that deepens
  the 8 ns deficit by more); together, if additive, ~9 eV of the open area's 15. None
  steepens the 5 → 8 ns rise by more than +4 eV against the 15–25 eV missing.
- **The conversion's intensity dependence has the wrong sign.** With Dewald's table the
  element-mean absorbed flux rises through the peak, so η falls (0.88 at 5.25 ns, 0.84
  at 8 ns, against `nif_conversion()`'s 0.85 / 0.90): +4 eV early, −4 eV late, the rise
  9 eV smaller. The coarse elements' mean flux sets η here, so the table is also mesh
  dependent ("Wall material, conversion and capsule energetics" §2). A constant η shows
  the converse: the rising `nif_conversion()` already supplies +4.3 eV of the rise.
- **The lips move the convention, not the drive.** Computed lips close faster in the rise
  (open fraction 0.754 against 0.809 at 5.0 ns): the initial-area reading drops 4.5–5.2 eV
  at 5–5.5 ns; the open-area reading and the capsule Tr move ≤ 0.8 eV.
- **The capsule albedo history is already the favorable one**: HYDRA's dark early ablator
  (0.01–0.16 to 6 ns) keeps the early drive 3 eV below a constant 0.3.
- **Dante's view of the spots weighs the early peak most.** Converted laser light seen
  directly is 31% of Dante 1's intensity at 5.0 ns, 29% at 5.25, 16% at 8.0 and 13%
  at 8.5 ns (the lips cover the spots near the hole, and the wall's re-emission grows).
  A term that dims the spots therefore lowers Dante early about twice as much as late,
  but removing the 15 eV open-area excess at 5.25 ns this way would need the spots'
  emission 72% lower (*L*: 1 − (252.1/267.3)⁴ = 0.21 of the intensity, over a direct
  share of 0.29), which no conversion datum allows.
- **Only the source lag acts on the shape**: −9 eV at 5.0 ns and +6 eV at 9.0 ns, the
  rise +9 eV, and the open-area delay 0.245 → 0.105 ns. As in "Computed LEH closure"
  (measured there on the 2026-10-24 energetics), a ~60 kJ store with a ~0.15 ns return
  helps both ends, but its content and return time are unpublished, so it stays a
  sensitivity. It makes the swing minimum shallower (−0.154 against −0.181; the
  backscatter history gives −0.152; Fig. 10 −0.170) and leaves the crossing.

**Reading.** The model's rise into the peak is in time with Dante on the Dante
convention (0.07 ns) and early by ~0.25 ns only on the open area. What the existing
inputs cannot supply is the rise through the peak (+47 to +57 eV against +72 eV from 5
to 8 ns). The measured power is flat over the plateau (440–448 TW over 5.45–8.05 ns), so
that rise needs losses that fall, or a source that grows, faster in time than the
model's (the wall's effective age and the conversion's t^0.12 are its two such terms),
or an LEH whose effective area in the thermal band closes more than the lips' 0.857 r₀
(*L*). No published number fixes either. Recommendation: no change; quote the early
peak on both conventions, and track the 5 → 8 ns rise as the quantity.

### 3. Sinks the model lacks, sized from published numbers

Sensitivity used for the drive (*L*, from the gas run below): a uniform 4.5% cut of the
coupled power lowers Dante 1 and the capsule Tr by 3.1–3.7 eV (1.2–1.3%, T ∝ P^0.28)
at every time.

| Sink | Published size | Model accounting | kJ | TW at peak | eV at 5.5 / 8.0 ns | Run |
|---|---|---|---|---|---|---|
| He gas fill, 0.3 mg/cc | 4.5% of HYDRA's laser energy heats the gas (Kritcher et al., PRE 109, 025204 (2024) App. 3, **H**); "∼5%" in gas and membranes together (Hurricane et al., PRL 132, 065103 (2024), **H**) | none: no gas; the bubble's absorbed light goes to the wall and converts | 80 (4.5% of the coupled 1778; 92 of the delivered) | ~17 (4.5% of ~380 TW coupled, if proportional to the power, *L*) | −3.3 / −3.7 (initial area); open −3.4 / −4.0 | `degrade=0.15` (0.89 × 0.955 of the delivered): capsule 13.9 → 13.7%, LEH 24.2 → 24.1%, Dante 1 at 8 ns 289.7 → 286.0 (open 312.9 → 308.9); Fig. 10 −0.175 (5.25), 6.18, +0.081, rms 0.048 |
| LEH windows (and the capsule tent) | ≲ 0.5%: Hurricane's ∼5% for gas and membranes less Kritcher's 4.5% for the gas (**H** inputs, **L** difference); 0.5 µm polyimide windows that burn through in the first few hundred ps of the foot (an LLNL target-fabrication update; Bradley et al., APS DPP 2008, abstract; **M**) | none (open holes from t = 0) | ≲ 9; two 0.5 µm, 7.5 mm² polyimide films (5.4 µg each), fully ionized at 2 keV, hold ~1 kJ each as heat (*L*) | 0 (gone in the foot) | 0 / 0 | not run (outside the peak; below the coarse mesh's resolution) |
| Laser scattered or reflected inside the hohlraum, re-deposited | Au disks absorb 0.93 / 0.95 / 0.91 of 3ω light at 3e13 / 1.5e14 / 3e14 W/cm² (Ze et al., UCRL-86414 (1981) Fig. 4, **H**): ≤ 5–9% leaves the first hit (backscatter into the lens included) | every surface absorbs all laser light at its first hit (`deposit`'s default), so nothing is reflected, lost or re-deposited | 5–11 out of the holes; 4–14 onto the capsule (below) | 1–2 lost; 1–3 onto the capsule | −0.4 to −0.9 (*L*: x-ray source −0.45 to −1.1% through the escape and the capsule's 0.3 conversion) | one diffuse bounce through the step's view factors (below) |
| Specular glint escaping the far LEH | inner-beam glint < 8 TW (< 2% of peak) at 0.3–0.6 mg/cc, none at 1.2 mg/cc (Lemos et al., PoP 29 (2022), doi:10.1063/5.0099937, abstract, **M**); outer-cone glint in the picket lands on the capsule poles and is "insignificant" for symmetry (Chen et al., PoP 31 (2024), doi:10.1063/5.0201962, abstract, **M**) | none | ≤ 24 if all of it escapes over the 5–8 ns peak (*L*) | ≤ 8 | at most −1.8 (≤ 2.1% of the coupled power, *L*) | not run (an upper bound) |

**Where diffusely scattered laser light would land** (each step's absorbed laser
re-emitted once, diffusely, through the view factors of that step; fractions of the
absorbed laser; *L*):

| t (ns) | wall | caps | capsule | out of the holes | lips |
|---|---|---|---|---|---|
| 4.5 | 0.724 | 0.126 | 0.071 | 0.067 | 0.013 |
| 5.0 | 0.724 | 0.116 | 0.087 | 0.059 | 0.013 |
| 5.5 | 0.724 | 0.117 | 0.086 | 0.058 | 0.015 |
| 6.5 | 0.722 | 0.133 | 0.061 | 0.064 | 0.020 |
| 8.0 | 0.718 | 0.144 | 0.046 | 0.068 | 0.023 |

With a scattered fraction r = 0.05–0.09 (Ze's disks), 0.3–0.6% of the coupled energy
leaves through the holes and 0.2–0.8% is deposited on the capsule as 3ω light (+0.2–0.8
points of capsule share if the ablator keeps it; the capsule share counts only x-rays
today). Both are below the 2.4–4.2 points of the capsule share and below the coarse
mesh's resolution of the drive. *Decision*: no opt-in diffuse laser reflectance (a 3ω
radiosity on the same view factors) was added; it would touch the deposition and the
history, and what it would move is ≤ 1 eV and ≤ 1 point.

**The capsule share and the LEH share move against each other** under every sink and
input above. `capsule_outside=0.15` (the blow-off absorptance that would close the
capsule share, "Capsule absorption, conversion and the ablating-wall release" §2) gives
15.8% (HYDRA 16.2) but takes the LEH share from 24.2 to 23.6% (HYDRA 25.2) and
Dante 1 at 8 ns down 2.9 eV (−4.2 at 9.0 ns); Fig. 10 +0.092 at 8 ns (+0.080 before). The
gas moves neither share (−0.2 / −0.1 points). Only a lower wall loss (DU factor ×0.9)
raises both (14.5 / 25.4%) and Dante 1 (+2.4 to +4.1 eV); Döppner's 6.5% equivalence
fixes the factor to 0.827–0.838 ("Energetics defaults and the capsule share"), so that
is not a lever without new data.

### What remains

- The rise through the peak: Dante 1 rises +72 eV from 5 to 8 ns, the model +47 (initial
  area) to +57 (open area); no existing input within its range adds more than +4 eV. The
  "early peak" is the open-area face of it; on Dante's convention the rise is in time
  (0.07 ns). A thermal-band LEH radius or a HYDRA loss history through the peak would
  separate the hole from the wall.
- After the laser falls: −30 / −38 / −47 eV at 8.5 / 9.0 / 9.35 ns (initial area);
  the plasma's return rate is unpublished ("Late-time drive" §4).
- The capsule share (13.9%, HYDRA 16.2) and the LEH share (24.2%, 25.2) together: no
  sink or input moves both up except a lower wall loss, which the DU calibration rules
  out.
- Of the sinks, the gas is 80 kJ and 3–4 eV at peak and at 8 ns (sized, not modeled); windows,
  diffuse scattering and glint are each ≤ 2 eV and ≤ 1 point.

## Drive histories: sliver steps and the albedo of a phase (2026-11-18) *measured*

**Sliver steps through a burnthrough.** A quasi-static `drive_history` of a hohlraum with
two thin gold windows (`BurnthroughAlbedo`, `WallLossAlbedo` wall, `wall_history` and
`wall_release` on, `transition = 0.05`, 11k elements: a Float32, iterative model) stopped in
a scratch study with "iterative radiosity solve (radiosity) did not converge in 400
iterations" at 0.1 ns steps; the same run at 0.2 ns steps, or with two more beams,
converged. The cause was that `transition` planned a step ending 7.8e-6 ns short of one of
`times`, and that remainder was taken as a separate step. The wall cooled slightly as the
windows began to burn through, so the wall elements were in `wall_release`'s re-emitting
branch, whose implicit part puts `C Δ^(b−1) φ(T)/H` (b − 1 = −0.38) into the albedo: at Δ =
7.8e-6 ns the factor `Δ^(b−1)` is 87 and the wall albedo was −27 on 7,310 elements. The
unpreconditioned restarted GMRES of a large model does not solve that system: on the failing
step's albedos and source it did not converge in 3,000 iterations either, while the same
albedo pattern scaled to Δ = 1e-4, 1e-3, 1e-2 and 0.1 ns (albedo minima −9.6, −3.4, −0.85,
0) converged in 80, 22, 11 and 9 iterations. Raising the cap would not help. The fix is in
the planner: a planned step that would leave less than `min_step` (1e-3 ns) before the next
of `times` runs to it instead (the transition then advances by at most a `min_step`'s worth
beyond the tolerance, inside the 1.5× check). With it the study's run completes (108 steps,
shortest 4.7e-3 ns, lowest albedo of any converged step −0.057). Public mock-up
(`test/burnthrough.jl`, "burnthrough steps leave no sliver"): a 760-element hohlraum with a
1 µm gold window 2 mm across, NIF quads, a pulse falling from 6 to 2 TW per quad between 0.6
and 0.8 ns while the window burns through (0.44 ns onward), an iterative model, and a goal
5e-6 ns past an inserted step: before the fix every such goal on a cooling step failed (13
of 13), after it all complete, with the lowest albedo −1.4. Given as times, the sliver's
albedo is −33, and on a direct (LU) model the gray solve then succeeds but the M-band's
thermalized group solve (an unpreconditioned GMRES as well) does not converge. Histories
without `transition` are unchanged (the N221204 example's CSV byte-identical), and with it,
only where a planned step ended within `min_step` of a goal.

**The wall albedo of a phase.** A linear scenario at one time of a pulse (a "phase":
the picket, the trough, the peak) needs albedos. Three sources disagree: a standalone
`solve(model, deposition; time)` iterates `WallLossAlbedo` to the element's present
temperature at the time since the pulse start, i.e. as if the wall had been at that
temperature since t = 0, without re-emission; a `Scenario` without `albedo` evaluates it
at the material's `reference_temperature` (250 eV) at `time`; `drive_history` follows
the temperature history of each element (`wall_history`) and its re-emission
(`wall_release`), and is the answer of the model itself. The public mock-up is the
diagnostic-windows tutorial's NIF-like hohlraum without windows (10.1 × 5.75 mm, 3.1 mm
entrance holes, 1,908 elements), 48 quads pointed as there, its pulse (picket 1.2 TW per
quad to 1.0 ns, trough 0.3 TW to 4.5 ns, second step 1.0 TW 5.0–6.0 ns, peak 9.5 TW 7.5–9.5
ns), history at 0.1 ns steps. Capsule drive relative to the history, and P2/P0, P4/P0 (%):

| time | wall albedo: history / standalone / scenario default | drive: standalone, default, wall mean | P2/P0: history, standalone, default, wall mean |
|---|---|---|---|
| picket 0.8 ns | 0.671 / 0.741 / 0.817 | 1.161, 1.425, 0.993 | −18.70, −17.07, −15.27, −18.18 |
| trough 3.0 ns | 0.858 / 0.807 / 0.889 | 0.846, 1.129, 0.998 | −14.54, −15.61, −13.45, −14.31 |
| second step 5.5 ns | 0.769 / 0.882 / 0.912 | 1.440, 1.642, 0.999 | −16.45, −13.69, −12.80, −16.30 |
| peak 8.5 ns | 0.872 / 0.936 / 0.925 | 1.338, 1.267, 0.999 | −14.02, −12.11, −12.41, −13.92 |

P4/P0 (%) at the four times: history 0.044, 0.155, 0.103, 0.164; standalone 0.084,
0.124, 0.170, 0.206; wall mean 0.047, 0.156, 0.103, 0.165. The history's albedos per
element, frozen, reproduce the history at the beams of the history (drive 1.0000, P2/P0 and
P4/P0 to the digits shown). In the trough the standalone albedo is the lower one (the
history's wall re-emits what it stored in the picket); everywhere else it is higher
(a wall heated in a low foot stores less energy than one at today's temperature since
t = 0, so it loses more). For the change a pointing perturbation makes (the inner cones
0.2 mm farther from the waist, the perturbed history the truth), the frozen per-element
albedos give ΔP2/P0 −1.086, −0.693, −0.891, −0.661 against the history's −1.099,
−0.699, −0.896, −0.664 percentage points (within 1.2%: the response of the albedos
themselves to the new drive is what a frozen state leaves out), the wall mean within 0.6%,
the standalone solve 13–29% off (low but in the trough, 16% high there) and the scenario
default 11–37% low. The choice made is that a phase
scenario uses the history's per-element albedos at the phase time (`history_state`),
interpolated between steps; the wall mean as one constant (`mean = :wall`) is within
0.5 percentage points of P2/P0 (largest at the picket, where the wall albedo varies
most, 0.647–0.682 over the wall) and 0.7% of the drive. Thin-foil windows need no
separate treatment: their albedo in the history carries the transition, and
`history_state` reports each window's transition and transmitted fraction. Confidence:
the comparison is the model against itself (which albedo the same physics gives a
fixed-time scenario); the accuracy of the wall-loss model is in "N221204 late-time drive
(2026-10-24)".
