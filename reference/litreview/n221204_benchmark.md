# N221204 (and N210808) benchmark inputs for Briareus: literature digest

Compiled 2026-10-22. Confidence: **H** = primary text read (numbers quoted from it), **M** = abstract/secondary,
**L** = inference or read by eye off a figure. Δλ is quoted at 1ω (inner − outer) unless marked 3ω.

## Addendum (PRE 109 025204, PRL 132 065103), what changed

Both papers were read in full: Kritcher et al. PRE 109, 025204 (2024) (key **K24E**) and Hurricane et al. PRL 132, 065103 (2024) (key **H24**).
Figs. 1, 9 and 10 of K24E were rendered and read by eye.

**Filled (H):**
- **Bang time** (K24E Table I, data): N210808 9.28 ± 0.07 ns, N221204 9.53 ± 0.07 ns. Simulated: 9.27 and 9.61 ns.
- **Capsule outer radius:** 1135 µm (K24E Fig. 1a label).
- **W dopant:** 0.44 at.% (N210808) and 0.54–0.6 at.% (N221204, batch KC952).
- **Hohlraum energy balance, N221204** (K24E App. 3, HYDRA calibrated to data).
  - Fractions of incident laser energy: 54.1 % heats the wall (4.4 % Au lining, 49.7 % DU), 25.2 % leaves as radiation through the LEHs, 4.5 % heats the He, 16.2 % is absorbed by the capsule.
  - The laser energy had to be degraded by 9–13 % to match data.
  - Downstream: 246 kJ (~12 %) is coupled to the implosion, fuel kinetic energy 15.33 kJ, 20 kJ of work on the DT, 8.2 kJ to the hot spot.
- **Capsule absorbed energy:** 258 kJ for N221204 (H24 Table I caption: Y/G_cap), and 19.4 kJ into the DT fuel.
  - H24's generic low-fill accounting: ~50 % wall, 30–35 % out the LEHs, 10–15 % capsule, ~5 % gas plus LEH windows.
  - H24: f₂ ≈ 0.07 asymmetry roll-up for N221204.
- **Modeling settings:** CBET clamp δn/n = 0.0075 after the start of the rise (K24E). The PoP rounds it to 0.008.
- **Outer-power alternative:** symmetry by power alone would need the outers cut to 150–200 TW.

**Filled (L, digitized by eye):**
- **N221204 per-cone input and post-CBET powers** (K24E Fig. 9a).
- **Input and post-transfer cone fraction vs time**, N221204 vs N220919 (Fig. 9b) and vs N210808 (Fig. 1e).
- **Simulated flux P2/P0(t)** for N210808, N220919 and N221204 (Fig. 10). This is a direct comparison target for Briareus moments.
- **Simulated Trad(t)** (Fig. 1c).

**Still not public:**
- Neither paper gives absolute pointing, defocus, CPPs, a per-cone (23.5°/30° or 44.5°/50°) wavelength split, window sizes or azimuths, corner/lip geometry, or the Au liner thickness.
- **Pointing, non-public (user, 2026-10-22):** an SST export of the N230729 (HYBRID-E, same
  platform) beam pointing is at `C:\Users\alexa\Dev\thor\I_Stag_DT_HyE_S08-BeamPointing.csv`
  (192 beams, aim points in chamber mm, quad and beam defocus 0). It reads with
  `read_visrad_sst(path)` and is the pointing to use for the benchmark. Do not copy its
  numbers into public documents.
- K24E treats Δλ as one inner−outer number.
- H24 lists "manufacturing gaps on the waist of the hohlraum" as an additional mode-1 source.

## Sources (short keys)

| Key | Reference | Read |
|---|---|---|
| K24 | Kritcher et al., Phys. Plasmas 31, 070502 (2024), doi:10.1063/5.0210904 | full text, H |
| AS24 | Abu-Shawareb et al., PRL 132, 065102 (2024), doi:10.1103/PhysRevLett.132.065102 (copy: api.creol.ucf.edu/Publications/18162.pdf) | full text incl. App. A–C, H |
| P24 | Pak et al., PRE 109, 025203 (2024), doi:10.1103/PhysRevE.109.025203; accepted MS, OSTI 2377242 (osti.gov/servlets/purl/2377242) | accepted MS, H |
| K22 | Kritcher et al., PRE 106, 025201 (2022), doi:10.1103/PhysRevE.106.025201 (www2.psfc.mit.edu/research/hedp/Papers/Kritcher_PhysRevE-2022.pdf) | full text + figures, H |
| Z22 | Zylstra et al., PRE 106, 025202 (2022), OSTI 1959535 | accepted MS, H |
| AS22 | Abu-Shawareb et al., PRL 129, 075001 (2022), OSTI 1880954 | accepted MS, H |
| K21 | Kritcher et al., PoP 28, 072706 (2021) (HYBRID-E), OSTI 1818404 | accepted MS, H |
| R26 | Ruof et al., "Constraining hydrodynamic models of ICF implosions using capsule surrogate experiments", PoP (2026), LLNL-JRNL-2019015, OSTI 3367476 (osti.gov/pages/servlets/purl/3367476) | MS dated 25 Apr 2026, H |
| DN24 | Di Nicola et al., "Delivering laser performance conditions to enable fusion ignition…", HEDP (2024), LLNL-JRNL-862033, OSTI 2440791 | accepted MS, H |
| KE26 | Kemp et al., "Modeling SBS backscatter from the inner laser cones…", PoP (2026), LLNL-JRNL-2014353, OSTI 3441608 | MS, H |
| SW25 | Swadling et al., "Build-A-Hohlraum", arXiv:2501.10350 | full text, H |
| B23 | Baker et al., PoP 30, 092708 (2023), doi:10.1063/5.0163396 (frustum, not HYBRID-E) | full text, H for frustum |
| K24E | Kritcher et al., PRE 109, 025204 (2024), doi:10.1103/PhysRevE.109.025204 (design paper) | full text + Figs. 1, 9, 10, H |
| H24 | Hurricane et al., PRL 132, 065103 (2024), doi:10.1103/PhysRevLett.132.065103 | full text, H |

**Not read:** Rubery PRL 132, 065104 (hohlraum reheating); MacGowan HEDP 40, 100944 (2021) and Milovich PPCF 63 (2021) on
diagnostic windows.

## 1. Target geometry

| Item | N210808 | N221204 | Src, conf. |
|---|---|---|---|
| Hohlraum | Au-lined DU, 6.4 mm ID × 11.24 mm | same | AS24 Fig. 1 & p. 3; K22 Fig. 4; K24 §II, H |
| LEH diameter | 3.1 mm (3.64 mm on N210207/N210307) | 3.1 mm ("3.11 mm" in R26) | AS24 Fig. 1; K22 Fig. 1a (radius 1.55 vs 1.82 mm); R26 §II, H |
| He fill | 0.3 mg/cc | 0.3 mg/cc | AS24 p. 3; K22 §II C, H |
| Au liner thickness, DU thickness | not given | not given | gap |
| Corner radius / LEH lip | Schematics show rounded corners and an "LEH washer" (Au insert); no dimensions | same | AS22 Fig. 1a; SW25 Fig. 3, L |
| LEH window | 500 nm polyimide + 25 nm Al flash ("standard on high-performance ID-ICF targets"); glued to an Al washer | same | SW25 §II.4, H |
| Storm window | 100 nm polyimide + 32 nm C, ~700 µm outside the pressure membrane; its Al retainer blocks ~27 % of Dante-2's view | same | SW25, H |
| Diagnostic windows | "three equatorial hohlraum diagnostic windows were repositioned" (new design first on N210808); "low Z" windows, sizes and location changed | same design | AS24 App. A; AS22; Z22 §I C; K22 §II C, H (no sizes or azimuths published in what I read) |
| Tent | two 45-nm Formvar membranes | same | AS24 p. 3; AS22, H |
| Fill tube | 2 µm diameter (was 5 µm) | 2 µm | P24 Table I; AS22, H |
| Capsule inner radius | 1050 µm (Fig. 4 label 1049 µm) | 1050 µm | K22 Fig. 4; P24 App., H |
| Ablator layers (in→out), µm | 6 / 18.7 doped / 54.9 (= 79.6) | 7.7 / 15.1 doped / 62.6 (= 85.4) | P24 Table I, H |
| Ablator radii | 1055, 1073.5, outer 1128.5 µm (OD 2.257 mm) | outer 1135 µm (K24E Fig. 1a, H); doped layer at 1057.7–1072.8 µm (L, from P24) | K22 Fig. 4; K24E Fig. 1a |
| Thickness change | — | +5.75 µm (N221204), +6.9 µm (N230729), +~5 µm more (N231029/N240210) | K24 §II, III; R26 abstract, H |
| W dopant | 0.42 % (K22) / 0.44 at.% (K24E), KC789, NCD; 7.9 %W·µm | 0.54–0.6 at.% (K24E), 0.6 at.% (AS24), 0.62 % (K24); NCD batch KC952; 9.2 %W·µm; ~1.19× optical depth | K22; K24E App. 2; P24 Table I & App.; AS24, H |
| Capsule mass | 3927 µg | 4253.1 µg ("4.25 mg") | P24 Table I; AS24, H |
| DT ice | 65 µm; inner radius 983 µm | 64.5 µm, 220 µg | K22 Fig. 4; AS24, H |
| DT gas / ice density | 0.45 mg/cc / 0.25 g/cc | same | AS24 Fig. 1, H |
| Inclusions (>200 µm³) | 5 | 105; >4000 small W–C inclusions | P24 Table I, H |
| Capsule mode 1 (wall) | 0.2 µm | 0.13 µm (<0.25 % of wall) | P24 Table I, H |

## 2. Laser

**Energy and power (H unless noted).**

| Shot | E_laser | P_peak | Δλ (1ω) | Src |
|---|---|---|---|---|
| N210307 | 1.91 MJ | 487 TW | 1.55 Å | AS24 T. I; K22 |
| N210808 | 1.89 MJ (AS24, P24, corrected) / 1.917 (AS22) / 1.92 (K24) | 441 TW | 1.8 Å | AS24 T. I; AS22 T. I |
| N220919 | 2.05 MJ | 440 TW | 2.5 Å | AS24; K24; P24 |
| N221204 | 2.05 MJ ± 2 % | 440 TW | 2.75 Å | AS24; DN24 |
| N230729 | 2.05 MJ | ~440 TW | 2.9 Å | K24 |
| N231029, N240210 | 2.2 MJ | 440 TW | 4.00 Å for the 2.2-MJ keyhole N231016 | K24 §II; R26 Table I |

AS24 App. B: the reported energies were "corrected a few percent downward due to neutron background on the laser calorimeters".
DN24: for N221204 the calorimeter read 3.7 % high before correction. N210808 is therefore 1.89 MJ in the 2024 papers and 1.917 MJ in 2022.

**Pulse shape, total power** (P24 Fig. 1 top, delivered N210808/N220919/N221204, read by eye, **L**).
- Picket: ≈55 TW at 0.7–1.3 ns, the same for all three shots.
- Trough: ≈28 TW. 1.8–2.5 ns on N210808, to ≈2.6 ns on N221204.
- Second step: ≈80 TW. Plateau ≈3.2–3.9 ns on N210808, ≈3.4–3.9 ns on N221204.
- Rise: from ≈4.0 ns to peak at ≈5.2 ns (N210808) or ≈5.4 ns (N221204).
- Peak: N210808 ≈430 TW, falling at ≈7.6 ns, half power ≈8.3 ns, off by ≈8.85 ns. N221204 ≈440–445 TW to ≈8.0 ns, half power ≈8.7 ns, off by ≈9.3 ns.

K24 (H): the peak was "extended by about 300 ps" and the trough lengthened. On N230729 the early pulse was a further ~30 ps longer.
K22: the N210808 trough was 100 ps shorter than N210307's.

**Per-cone delivered power, N221204** (DN24 Fig. 1 right, "Outer cone – LEH reported (128 beams)" and "Inner cone (64 beams)"
against "Specified", right axis in TW, read by eye, **L**):

| t (ns) | 0.4 | 0.8–1.3 | 2.0 | 2.5 | 3.3 | 3.9 | 4.5 | 5.0 | 5.5 | 8.0 | 8.3 | 8.6 | 8.9 | 9.1 | 9.25 | 9.35 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Outer (TW) | ~2 | 50 | 25 | 26 | 60 | 65 | 150 | 270 | 296 | 294 | 250 | ~110 | ~40 | 0 | 0 | 0 |
| Inner (TW) | ~2 | 4 | 4.5 | 5 | 7 | 9 | 25 | 90 | 146 | 146 | 146 | 146 | 146 | 146 | ~140 | 0 |
| CF (L) | — | 0.07 | 0.15 | 0.16 | 0.10 | 0.12 | 0.14 | 0.25 | 0.33 | 0.33 | 0.37 | 0.57 | 0.78 | 1 | 1 | — |

- Cross-check: DN24 Fig. 1 left gives per-quad power of ≈9.3 TW for both cones at peak (32 outer and 16 inner quads), so 297 + 149 TW = 446 TW and CF = 1/3.
- R26 (H): "the cone fraction during the peak power is limited to 0.33, where the peak laser power is distributed evenly to all the beams".
- Requested vs delivered energy, N221204 (DN24 Fig. 1, H).
  - Inner: +1.5 %, −2.6 %, −7.4 %, −4.8 % over the foot epochs; +0.6 % at peak.
  - Outer: +1.1 %, +0.9 %, +0.6 %, +1.5 %; −0.5 % at peak.
- Quad-to-quad power balance rms, N221204 (DN24, H).
  - Inner: 12.1 → 13.7 → 10.0 → 5.6 % in the foot, 2.2 % at peak.
  - Outer: 6.4 → 6.8 → 3.7 → 3.4 %, 2.1 % at peak.
  - Requirements: foot 22 % inner / 14 % outer, peak 2.5 % / 3.5 %. Pointing accuracy <50 µm.

**Cone fraction, other sources.**
- R26 Fig. 5 (N210808 keyhole and 2DConA pulses, read by eye, L): CF ≈0.10 in the picket, ≈0.16 at ≈2.3 ns, ≈0.10 at 3 ns.
- "Rise CF" at ≈4 ns (H, R26 Tables I, III):
  - 0.21 for N210808 and for the N221204 keyhole N230110.
  - 0.17 on the N221204 2DConA surrogate N240304.
- Peak CF ≈0.33 (L), then →1 as the outers droop.
- K22 Figs. 1(e) and 3(a) show the same shape for N210808, incident and post-CBET.
  - Incident CF (L): ≈1 at t≈0, ≈0.1 at 1–3 ns, a bump to ≈0.17 at ≈2.3 ns, ≈0.2 at 5 ns, 0.33 over 6–8 ns, back to 1 at ≈9 ns.
  - Post-CBET CF (L): ≈0.2 in the foot, ≈0.5–0.6 at peak.
- K22 Fig. 1(d), per-cone incident peak powers (L): ≈60 TW (23°), ≈60 TW (30°), ≈110 TW (44°), ≈115 TW (50°).
  - These are simulation powers after multipliers. They sum to ≈0.8 × 441 TW, consistent with R26's peak-power multipliers of 0.825–0.86.

**Design per-cone powers and cone fraction, N221204** (K24E Fig. 9; all values read by eye, **L**).

Fig. 9a shows "input laser powers (solid) and calculated … post energy transfer (dashed)" for each cone, N221204 (red) vs N220919 (black).
Fig. 9b shows the requested input and post-transfer cone fraction. Fig. 1(d, e) show the same for N210808, and Fig. 1(b) the total, outer and inner powers.
Times are in ns, powers in TW per whole cone.

| t (ns) | 0.5 | 0.7–1.3 | 2.0 | 2.5 | 3.3 | 3.9 | 4.5 | 5.0 | 5.5 | 7.0 | 8.0 | 8.5 | 9.0 | 9.3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 23° input | ~1 | ~2 | ~2 | ~2 | ~4 | ~5 | ~12 | ~40 | 62 | 64 | 64 | 64 | 64 | 0 |
| 30° input | same as 23° | | | | | | | | | | | | | |
| 44° input | ~10 | 25 | 12 | 12 | 33 | 38 | 75 | 120 | 125 | 125 | 125 | ~70 | ~10 | 0 |
| 50° input | same as 44° (≈2–3 TW higher at peak) | | | | | | | | | | | | | |
| CF input | 0.10 | 0.07–0.08 | 0.15 | 0.13 | 0.12 | 0.17 | 0.20 | 0.25 | 0.33 | 0.33 | 0.33 | 0.45 | 0.9 | 1.0 |
| CF post-CBET | — | 0.17 | 0.22 | 0.20 | 0.20 | 0.30 | 0.40 | 0.55 | 0.65 | 0.57 | 0.53 | 0.52 | 0.9 | 1.0 |
| 23° post-CBET | | | | | | | | | 125 (peak ≈5.6 ns) | 105 | 100 | ~85 | | |
| 44° post-CBET | | | | | | | | | 80 (5.3 ns), dip 55 (6 ns) | 70 | 75 | | | |
| 50° post-CBET | | | | | | | | | 95, dip 85 | 95 | 100 | | | |

- Peak input per cone sums to ≈64 + 64 + 125 + 128 ≈ 380 TW (L), below the 440 TW delivered. HYDRA plots appear to carry the 9–13 % laser degradation that K24E App. 3 says is needed.
  - DN24's measured split at peak is 146 TW inner and 296 TW outer, i.e. 73 TW per inner cone and 148 TW per outer cone.
  - Use DN24 for absolute powers and K24E for the design CF shape.
- K24E text (H): the cone fraction "in the peak (33%)". The outers "droop" late while the inners are held, so CF → 1 at ≈9.2 ns.
- Foot CF for N221204 was *raised* relative to N220919 (≈0.13 vs 0.09 at ≈3.4 ns in Fig. 9b, L). This followed keyhole N221106, which showed less foot asymmetry than modeled; N230110 verified it.
- Fig. 1e (L): N221204 input CF lies within ≈0.02 of N210808 through the foot and peak. The N221204 droop (CF rise) is ≈0.3 ns later. Its post-transfer CF in the peak is ≈0.05 higher (≈0.6 vs ≈0.55).

**Pointing, spots, wavelengths.**
- K24E and H24 give no absolute or relative pointing, defocus or CPP values, and no per-cone color split. K24E Fig. 1a shows only the cone geometry (23°, 30°, 44°, 50°), 3.2 mm radius and 1.55 mm LEH radius.
- 50° beams moved toward the midplane to clear the 3.1 mm LEH. K22 Fig. 1a labels the move "250 µm vs N210307" (H).
- Outer-quad splitting (K21, H): the focal points of the four beams in each outer quad are separated in azimuth by 450 µm. The 44° and 50° focal points are separated by >430 µm along the axis.
- Absolute aim points, defocus and CPP types per cone: **not published** in any source read (gap).
- NIF tuning ranges (DN24, H): outer 1052.32–1052.95 nm, inner 1052.85–1053.7 nm (1ω). Since 2020 a 4th oscillator gives the 50° cone its own wavelength (the 44.5°/50° split, used against 50° SBS).
- Whether N221204 split 23.5°/30° or 44.5°/50°: **not stated** in any source read. K24, AS24 and P24 give a single Δλ (gap; likely in PRE 109 025204).
- Beam smoothing: no per-shot data found.

**Backscatter.**
- N221204, SBS: 42 ± 11 kJ (1σ), "early in peak power and mostly on the inner cones" (AS24, H).
- N221204, SRS: 6 ± 1 kJ, throughout peak (AS24, H).
- N221204, coupled fraction: 98 % ± 0.6 % (AS24, H).
- N210808, inner-cone SBS: 2.8 kJ/quad (23°) and 0.9 kJ/quad (30°), i.e. ≈30 kJ total over 8 + 8 inner quads (KE26 §II, H; total L).

## 3. Measured observables

- **Dante-1 Trad(t)** (P24 Fig. 1 bottom; read by eye, L).
  - N221204: ≈150 eV at 3 ns, ≈200 eV at 4.6 ns, ≈275 eV at 6 ns, peak ≈311 eV at ≈8.2 ns.
  - N210808-pulse repeats: peak ≈310 eV, ≈0.3 ns earlier.
  - "Dante 1 data are not available for N210808"; traces stop 100 ps before peak nuclear production.
  - AS24: N221204 reached "a 2 % higher peak x-ray drive temperature" (H).
- **"Drive Trad"** (AS24 Table I, H): 304 (N210207), 306 (N210307), 307 (N210808, average of similar shots), 313 (N220919), 313 eV (N221204).
  - AS22 gives the N210808 value as computed "from the measured laser pulse", not measured directly.
- **Coast time** (AS24, H): 0.93 ns (N210808), 0.89 ns (N220919), 0.78 ns (N221204). Defined from mid-fall of the laser to peak compression (AS22).
- **Bang time** (K24E Table I, measured, H): N210808 9.28 ± 0.07 ns, N221204 9.53 ± 0.07 ns. Post-shot HYDRA: 9.27 and 9.61 ns.
  - The earlier estimate from t_coast (≈9.2 and ≈9.5 ns) agrees.
- **Simulated Trad(t)** (K24E Fig. 1c, L): ≈130 eV at 1 ns, ≈125 eV at 2.5 ns, ≈200 eV at 4.6 ns, ≈260 eV at 6 ns, ≈300 eV at 8.5 ns. N221204 is slightly higher and longer than N210808 after 8 ns. Reheating from fusion shows after ≈9.3 ns.
- **Simulated capsule flux P2/P0(t)** (K24E Fig. 10, post-transfer, L): this is the most direct comparison for Briareus's flux moments. Negative means waist-hot.
  - ≈0 (±0.01) over 0.5–3 ns.
  - ≈+0.03–0.04 at ≈3.6 ns (N220919 +0.09).
  - Minimum near 5.5 ns: N221204 ≈−0.17, N220919 ≈−0.15, N210808 ≈−0.10.
  - Back through 0 at ≈6.5 ns.
  - ≈+0.10 at 7.5 ns and ≈+0.12–0.15 at 8.0–8.3 ns (pole-hot).
  - K24E: P2 is held in spec in the foot and "allowed to swing" in the peak. The swing comes from wall motion, which a static Briareus geometry will not reproduce.
- **Energy partition** (K24E App. 3, HYDRA calibrated, H): see the addendum (54.1 / 25.2 / 4.5 / 16.2 %).
  - R26 2DConA stagnation (H): 9.170 ± 0.061 ns for N220618 (N210808 surrogate) and 9.876 ± 0.090 ns for N240304 (N221204 surrogate, Δλ 2.9 Å).
- **Hot-spot P2 (PNI)** (K24 Table II, H): N220919 ~−10.5 µm (oblate), N221204 ~0, N230618 ~0, N230729 ~+6, N230904 ~+6 µm.
  - P24 (H): N210808 −2.3 ± 0.5 µm, N211024 +6.5 ± 1.1 µm, N220919 |P2| 10.8 ± 1.6 µm, N221204 0.3 ± 0.1 µm (17 % contour).
- **Hot-spot velocity, magnitude and direction (θ, φ)** (P24 Fig. 2, H): N210808 68 km/s (166°, 167°); N220919 51 (74°, 71°); N221204 36 (151°, 54°); N211024 133 (154°, 313°).
  - K24 Table II: N221204 36 ± 6 km/s (~27 without α heating); N230729 104 ± 9 km/s.
  - The Fig. 2 skymaps also mark the "mode-1 drive" direction, which gives a mode-1 drive vector to compare with Briareus's flux P1 including windows.
- **In-flight shell (2DConA, R26 Tables III–IV, H).** Values are data (simulation).

| Shot (DT design) | Δλ | rise CF | v_shell (km/s) | stag. (ns) | ⟨shell P2⟩ (µm) | ⟨shell P4⟩ | HS P0 / P2 / P4 (µm) |
|---|---|---|---|---|---|---|---|
| N220618 (N210808) | 1.80 | 0.21 | 392 ± 24 (401) | 9.170 ± 0.061 | 9.07 ± 0.30 (11.1) | 1.41 ± 1.89 | 80.2 / 5.7 ± 0.49 / 0.56 |
| N221010 (N210808) | 1.80 | 0.16 | 381 ± 15 (408) | 9.230 ± 0.029 | 10.92 ± 0.87 (11.9) | 3.12 ± 1.15 | 67.7 / 8.67 ± 0.34 / −0.4 |
| N230827 (N210808) | 1.80 | 0.24 | — | 9.240 ± 0.190 | 27.5 ± 2.37 (25.8) | 3.76 ± 3.19 | 65.3 / 25.4 ± 0.93 / 6.14 |
| N240304 (N221204) | 2.90 | 0.17 | 355 ± 5 (344) | 9.876 ± 0.090 | 21.7 ± 1.31 (19.5) | 4.16 ± 1.04 | 71.0 / 7.39 ± 1.85 / 0.07 |

- R26 Fig. 9 (N240304 P2 vs time, L): shell P2 ≈22–24 µm at 9.15 ns falling to ≈16 µm at 9.45 ns, at shell P0 270→160 µm.
- R26: matching these data needs HYDRA peak-power multipliers of 0.825–0.86 and peak CBET clamps of 0.005–0.007, varying "5 % and 30 %" shot to shot.
- **Shock timing, keyhole** (R26 Table II, H), N230110 = the N221204 design: 85.6 µm, Δλ 2.75 Å, rise CF 0.21.
  - Equator: shock-1 breakout 3.864 ± 0.012 ns, 1–2 merge 4.847 ± 0.003 ns, 1–2–3 merge 5.972 ± 0.008 ns; velocities 32.74 and 50.30 km/s.
  - Pole: 3.906 ± 0.018, 4.857 ± 0.01 and 6.050 ± 0.046 ns; velocities 34.00 and 49.90 km/s.
  - N221106 (N210808 design): equator 3.714 / 4.725 / 5.700 ns, pole 3.650 / 4.655 / 5.809 ns.
- **Implosion velocity:** 380–400 km/s (K24, AS24, H). K24 Table III (simulated N230729): 386 km/s, 245 kJ absorbed by the capsule.
- **Capsule absorbed energy** (H24, H): 258 kJ for N221204. For N210808, 1.37 MJ / G_cap 5.8 ≈ 236 kJ (L). G_cap is simulation-based.
  - K24E App. 3: the capsule absorbs 16.2 % of incident laser energy, i.e. ≈330 kJ of 2.05 MJ before the 9–13 % degradation (L). 246 kJ of that is coupled to the implosion.
  - K24E Table I simulated implosion velocity: 401 km/s (N210808) and 383 km/s (N221204). Remaining ablator mass: 3.8 % and 5.7 %.
- **LEH closure:** no N221204 GLEH data found. SW25 (H, BAH targets): simulations overpredict LEH size early and match it at peak.
  - SW25: measured bang times are 400–700 ps later than simulated across ID-ICF.
- **Sensitivities.**
  - 40–50 µm hot-spot P2 per Å (K21, H; already in research notes).
  - K22 Fig. 3b (N210808-class): ΔP2 vs model at Δλ = 0 is ≈27 µm at 1.25 Å, 43–48 µm at 1.55 Å and ≈56 µm at 1.75–1.8 Å (L).
  - N220919→N221204: +0.25 Å with re-tuned powers moved hot-spot P2 from −10.8 to +0.3 µm (P24, H). That is ≈44 µm/Å, confounded by the power changes (L).
  - Rise CF (R26 Fig. 10, H): 0.21→0.24 raised the shell P2 from 9.1 to 27.5 µm (≈6 µm per 0.01 CF). 0.21→0.16 left it unchanged (10.9 µm); HYDRA had predicted −10 µm.
  - Outer pointing (B23, frustum, H): 7–12 µm hot-spot P2 per 100 µm of outer z-pointing; 44 µm/Å.

## 4. Shot-to-shot deltas

| | N210808 | N220919 | N221204 | N230729 | N231029 / N240210 |
|---|---|---|---|---|---|
| E_laser (MJ) | 1.89–1.92 | 2.05 | 2.05 | 2.05 | 2.2 |
| P_peak (TW) | 441 | 440 | 440 | ~440 | 440 |
| Pulse | ref. | +~300 ps peak, longer trough | same as N220919, re-tuned powers | early pulse +~30 ps | longer |
| Δλ 1ω (Å) | 1.8 | 2.5 | 2.75 | 2.9 | (4.00 on keyhole N231016) |
| Ablator (µm) | 79.6 | ~85.4 (+6) | 85.4 (+5.75) | +6.9 vs N210808 | +~5 vs N230729 |
| W | 0.42 % | 0.62 % | 0.62 % | ~0.4 % | ? |
| LEH (mm) | 3.1 | 3.1 | 3.1 | 3.1 | 3.1 (L) |
| Yield (MJ) | 1.37 (AS22) / 1.3 ± 0.07 | 1.2 ± 0.06 | 3.1 ± 0.16 (3.15) | 3.88 | 5.2 (record, N240210 per K24 abstract) |
| HS P2 (µm) | −2.3 ± 0.5 | −10.5 / −10.8 | ~0 / 0.3 | ~+6 | ? |
| HS velocity (km/s) | 68 ± 7 | 51 ± 9 | 36 ± 6 | 104 ± 9 | ? |
| v_imp | 380–400 km/s design; −15 km/s for +6 µm (K24 Table I, simulated) | | | 386 (sim.) | |
| t_coast (ns) | 0.93 | 0.89 | 0.78 | — | — |
| Bang time (ns) | 9.28 ± 0.07 | — | 9.53 ± 0.07 | — | — |
| Drive Trad (eV) | 307 | 313 | 313 | — | — |

## 5. Briareus input deck (N221204; N210808 differences in brackets)

```julia
# Geometry (mm). Briareus takes radii, not diameters.
hohlraum(; length = 11.24, radius = 3.2, leh_radius = 1.55,        # [same]   H
          capsule_radius = 1.135,                                   # [1.1285] H (K24E Fig. 1a; K22 Fig. 4)
          corner_radius = ?,                                         # GAP: not published
          windows = (?),                                             # GAP: 3 equatorial low-Z windows, sizes/azimuths unpublished
          wall = <Au-lined DU>, capsule = ablator())                 # Au thickness: GAP
# The LEH is open in Briareus. Real LEH: 0.5 µm polyimide + 25 nm Al; storm window 0.1 µm polyimide, +0.7 mm.
# Laser: 192 NIF beams, NIF ports; groups :inner23 (32 beams), :inner30 (32), :outer44 (64), :outer50 (64).
times = [0.0, 0.4, 0.8, 1.3, 2.0, 2.5, 3.3, 3.9, 4.5, 5.0, 5.5, 8.0, 8.3, 8.6, 8.9, 9.1, 9.25, 9.35]   # ns, L (DN24 Fig. 1)
outer = [0, 2, 50, 50, 25, 26, 60, 65, 150, 270, 296, 294, 250, 110, 40, 0, 0, 0]                  # TW, 128 beams
inner = [0, 2, 4, 4, 4.5, 5, 7, 9, 25, 90, 146, 146, 146, 146, 146, 146, 140, 0]                     # TW, 64 beams
total = outer .+ inner;  fraction = inner ./ max.(total, eps())
pulses = cone_fraction_pulses(beams, times, total, fraction)    # per-beam TW
# Absolute powers follow DN24 (measured). The CF shape agrees with K24E Fig. 9b, the design input (L):
#   0.07–0.08 picket, 0.15 at 2 ns, 0.12 at 3.3 ns, 0.17–0.20 at 4–4.5 ns, 0.33 for 5.5–8.0 ns, then → 1 by 9.2 ns.
# Split inner power equally over 23°/30° and outer equally over 44°/50° (K24E Fig. 9a: within a few TW).
# [N210808: same picket/trough, peak ≈430–441 TW ending ≈0.45 ns earlier (P24 Fig. 1); CF from K22 Fig. 1e / K24E Fig. 1e,
#  within ≈0.02 of N221204 until the droop, which starts ≈0.3 ns earlier.]
# Validation clock: bang 9.53 ± 0.07 ns [9.28 ± 0.07]; compare flux moments before ≈9.3 ns (reheating after).
# Detuning, Beam.dlambda at 3ω = Δλ(1ω)/3, inner cones relative to outer:
dlambda_inner = 2.75 / 3     # = 0.917 Å   [N210808: 0.600; N220919: 0.833; N230729: 0.967; 2.2 MJ: 1.333]
dlambda_outer = 0.0          # GAP: any 44.5°/50° or 23.5°/30° split on N221204 unknown
# Pointing: 50° cone ≈0.25 mm closer to the midplane than N210307 (K22). Outer-quad beams split 0.45 mm in azimuth;
# 44°/50° foci >0.43 mm apart in z (K21). Absolute aims, defocus and CPPs: GAP (use the repo's current HYBRID-E defaults).
# Backscatter (subtract from incident): inner SBS 42 ± 11 kJ, early in peak; SRS 6 ± 1 kJ over peak.
# CBET (HYDRA practice): δn/n = 0.002 until the rise, 0.0075 after (K24E App. 1, H); 2DConA fits need 0.005–0.007 (R26).
# Energetics: HYDRA needed the laser degraded 9–13 % (K24E) or peak multipliers 0.825–0.86 (R26) to match data.
# Use the same knob, not the wall albedo, if Briareus's absolute drive comes out high.
```

## 6. Validation targets

| Observable | Value | Time | Source |
|---|---|---|---|
| Peak Trad (Dante-1, through LEH) | ≈311 eV (N221204); ≈2 % above N210808 | ≈8.2 ns | P24 Fig. 1 (L); AS24 (H) |
| Drive Trad (inferred) | 313 eV (N221204), 307 eV (N210808) | peak | AS24 Table I, H |
| Trad(t) foot/rise | ≈150 eV @ 3 ns, 200 @ 4.6, 275 @ 6 ns | — | P24 Fig. 1, L |
| Capsule absorbed energy | 258 kJ (N221204); ≈236 kJ (N210808) | integrated | H24 (H); N210808 from G_cap (L) |
| Laser energy fractions (N221204) | wall 54.1 % (Au 4.4, DU 49.7), LEH 25.2 %, He 4.5 %, capsule 16.2 % | integrated | K24E App. 3, H (HYDRA calibrated) |
| Flux P2/P0 (simulated) | ≈0 foot; ≈−0.17 at 5.5 ns; 0 at 6.5 ns; ≈+0.12–0.15 at 8 ns | 0.5–8.5 ns | K24E Fig. 10, L (model, not data) |
| Bang time | 9.53 ± 0.07 ns (N221204); 9.28 ± 0.07 ns (N210808) | — | K24E Table I, H |
| Laser coupled (1 − SBS) | 98 % ± 0.6 % | — | AS24, H |
| In-flight shell P2 (N221204 surrogate) | 21.7 ± 1.3 µm | 9.1–9.5 ns, R ≈ 270→160 µm | R26 T. IV, H |
| In-flight shell P2 (N210808 surrogate) | 9.07 ± 0.30 µm; 27.5 ± 2.4 µm at rise CF 0.24 | ~8.6–9.1 ns | R26, H |
| ΔP2 per 0.01 rise CF | ≈+6 µm (0.21→0.24); ~0 (0.21→0.16) | in flight | R26 Fig. 10, H |
| Hot-spot P2 | N210808 −2.3 ± 0.5; N220919 −10.8 ± 1.6; N221204 0.3 ± 0.1; N230729 ≈+6 µm | stagnation | P24; K24, H |
| ΔP2 / ΔΔλ | ≈40–50 µm/Å (hot spot) | stagnation | K21 (H); P24 pair (L) |
| Mode-1 (hot-spot velocity) | 36 ± 6 km/s at (151°, 54°) (N221204); 68 km/s at (166°, 167°) (N210808) | stagnation | P24 Fig. 2; K24, H |
| Mode-1 drive → velocity | 0.7 % applied mode 1 ↔ >100 km/s burn-on | — | K24 Fig. 5, H |
| Foot shock symmetry | pole − equator shock-1 breakout +42 ± 22 ps (N230110) | 3.9 ns | R26 Table II, H (difference L) |
| Coast time / 2DConA stagnation | 0.78 ns; 9.876 ± 0.090 ns | — | AS24; R26, H |

A view-factor code predicts flux moments, not shell or hot-spot P2 in µm. Use the ΔP2 rows (Δλ, rise CF, cone-fraction changes)
through a transfer function, and the Trad and capsule-energy rows for absolute energetics.

## 7. Open questions and figure-only data

1. **Tabulated cone powers.** None of the public sources tabulates inner or outer powers vs time, and no data release turned up (no OSTI dataset or Zenodo record).
   - K24 data availability: "available from the corresponding author upon reasonable request".
   - Digitizable figures:
     - DN24 Fig. 1 right (OSTI 2440791, p. 2): delivered and specified inner (64-beam) and outer (128-beam) power, N221204, 0–10 ns, 0–300 TW axis. The cleanest source.
     - DN24 Fig. 1 left: per-quad power, 0–12 TW axis, all 48 quads overlaid.
     - P24 Fig. 1: total power for N210808/N220919/N221204, and Dante-1 Trad.
     - K22 Fig. 1b, d, e: N210808 total/inner/outer power, per-cone incident and post-CBET, cone fraction.
     - R26 Figs. 3, 5: cone fraction (right axis, 0–0.4) for the keyhole and 2DConA surrogates.
     - AS24 Fig. 1: N210808 and N221204 power plus simulated Tr.
     - K24 Fig. 1 top right: inner/outer powers, N221204 vs N210808.
     - K24E Fig. 1(b, d, e): total, outer and inner power; 23° input vs post-transfer; CF; N221204 vs N210808. Digitized in §2.
     - K24E Fig. 9(a, b): all four cones, input vs post-transfer, and CF, N221204 vs N220919.
     - K24E Fig. 10: simulated flux P2/P0(t) with scaled power and CF.
2. **PRE 109, 025204** is now read. It gives no per-cone Δλ split, pointing, defocus or CPP. Its per-cone "input" powers sum to ≈380 TW at peak against 440 TW delivered. They are probably HYDRA-degraded, but the paper does not say so explicitly.
3. **Windows.** No public sizes, azimuths or materials for the three equatorial diagnostic windows; K24E and H24 do not add any, and the "THOR" naming did not appear in any source. H24 adds "manufacturing gaps on the waist of the hohlraum" as a mode-1 source. Look in MacGowan HEDP 40, 100944 and Milovich PPCF 63 (2021).
4. **Hohlraum details.** Corner and lip geometry, Au liner and DU thicknesses, LEH washer dimensions: unpublished in what I read. The SW25 Fig. 3 CAD (BAH hohlraum, 5.75 mm) shows the LEH hardware layout only.
5. **Per-cone Δλ splits** (23.5°/30°, 44.5°/50°) for N221204: unknown. The 4-color capability existed (DN24).
6. **Pointing, defocus, CPPs** for HYBRID-E: only relative statements (50° moved 250 µm; quad split 450 µm; cones >430 µm apart).
7. **Energy bookkeeping.** N210808 at 1.89 (corrected) vs 1.92 MJ. Use the 2024 values for both shots.
8. **2DConA surrogate.** N240304 used Δλ 2.9 Å (not 2.75) and rise CF 0.17, so it is a surrogate of the N221204 *design*, not of the shot.
9. **Bang times:** resolved (K24E Table I).
10. **Liner thickness:** still unpublished. K24E only splits the wall energy between the Au lining (4.4 % of laser energy) and the DU (49.7 %).
