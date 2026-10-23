# HYBRID-E hohlraum energetics for Briareus: literature review (2026-10-22)

> **Addendum (PRL 132 065103, PRE 109 025204).** Both papers were read in full text: Hurricane
> et al. ("Hu24") and Kritcher et al. ("Kr24b"). Changes:
> - New: the HYDRA energy balance for N221204, Kr24b Appendix §3, **H**. Of the incident laser energy:
>   - wall heating 54.1% (4.4% in the Au lining, 49.7% in DU)
>   - LEH radiation loss 25.2%
>   - He gas 4.5%
>   - **capsule absorbed 16.2%**
>
>   In addition, the simulated laser is degraded by 9–13%. The "useful energy coupled to the implosion
>   (and not in the ablated capsule material)" is ≈246 kJ (≈12%).
> - New: Hu24's generic low-fill accounting, citing Callahan et al., PoP 27, 072704 (2020), **H**:
>   ~50% wall, 30–35% LEH, 10–15% capsule, ~5% gas plus windows. Also E_abs = Y/G_cap = 258 kJ for
>   N221204.
> - §2.4 is revised. The 245–258 kJ "capsule absorbed energy" is the energy left in the unablated
>   capsule. The total x-ray energy absorbed by the capsule is 16.2% ≈ 300–330 kJ. My 316–361 kJ
>   estimate (f_cap 0.7–0.8, outer area, Dante Trad) therefore matches the total, and the "30–45%
>   overshoot" is resolved as a definition mismatch.
> - Kr24b Fig. 1(c) adds the simulated internal Trad(t) for N221204.
> - The LEH-loss bound is now compared with HYDRA's 25.2%. The validation table, §4 item 3 and the
>   open questions are updated.

Confidence follows `reference/research-notes.md`: **H** read in the primary text; **H-img** read from a
rendered page or figure (numbers digitized by eye are ±5 eV or ±10 TW); **M** abstract or secondary;
**L** inference or my own arithmetic. Rosen's "radiation hohlraum units" (rhu): T in heV, area in mm²,
time in ns, energy in hJ, σ ≈ 1 (1.03). BAH's units: σ = 10.26 TW cm⁻² heV⁻⁴.

Downloaded texts and rendered pages were kept only in the session scratchpad (not in the repository).

## 1. Sources read

| Key | Reference | Access | Conf. |
|---|---|---|---|
| Kr22 | A. L. Kritcher et al., PRE 106, 025201 (2022), doi:10.1103/PhysRevE.106.025201 | published PDF (earlier session's copy); Table I rendered | H, H-img |
| Kr24 | A. L. Kritcher et al., Phys. Plasmas 31, 070502 (2024), doi:10.1063/5.0210904 | given text | H |
| Kr24b | A. L. Kritcher et al., PRE 109, 025204 (2024), doi:10.1103/PhysRevE.109.025204 | APS featured-article PDF (full text); Fig. 1 rendered | H, H-img |
| Hu24 | O. A. Hurricane et al., PRL 132, 065103 (2024), doi:10.1103/PhysRevLett.132.065103 (CC-BY) | APS featured-article PDF (full text) | H |
| AS24 | H. Abu-Shawareb et al., PRL 132, 065102 (2024), doi:10.1103/PhysRevLett.132.065102 | https://www2.psfc.mit.edu/research/hedp/Papers/Abu-Shawareb_PRL-2024.pdf; Table I rendered | H, H-img |
| Pak24 | A. Pak et al. (as "Allen ... Zylstra"), PRE 109, 025203 (2024), doi:10.1103/PhysRevE.109.025203 | preprint LLNL-JRNL-856035 (OSTI 2377242); Fig. 1 rendered | H, H-img |
| Zy22 | A. B. Zylstra et al., PRE 106, 025202 (2022) | https://www.osti.gov/servlets/purl/1959535; Fig. 2 rendered | H, H-img |
| Ru24 | M. S. Rubery et al., PRL 132, 065104 (2024), "Hohlraum reheating from burning NIF implosions" | https://www.osti.gov/servlets/purl/2427473 | H |
| BAH | G. F. Swadling et al., Phys. Plasmas 32, 052707 (2025), "Build-A-Hohlraum" | arXiv:2501.10350; p. 3 rendered | H, H-img |
| Ch24 | H. Chen et al., PRE 110, L013201 (2024), "ViewFactor" | https://www.osti.gov/servlets/purl/2429360 | H |
| Ch20 | H. Chen et al., Phys. Plasmas 27, 022702 (2020), gated LEH images | https://www.osti.gov/servlets/purl/1632379 | H |
| Do15 | T. Döppner et al., PRL 115, 055001 (2015), DU hohlraums | https://www.osti.gov/pages/servlets/purl/1209256 | H |
| Sc07 | J. Schein et al., PRL 98, 175003 (2007), cocktail hohlraums | https://www.osti.gov/servlets/purl/908388 (UCRL-JRNL-229156) | H |
| Ro05 | M. D. Rosen, "Fundamentals of ICF hohlraums" (2005) | https://www.osti.gov/servlets/purl/885380 | H |
| Li04 | J. D. Lindl et al., Phys. Plasmas 11, 339 (2004), Sec. IV | Warwick copy (physplasmas_11_339.pdf) | H |
| To11 | R. P. J. Town et al., Phys. Plasmas 18, 056302 (2011), NIF energetics analysis | LLNL-JRNL-463439, UNT copy | H (citation: M) |
| De08 | E. L. Dewald et al., Phys. Plasmas 15, 072706 (2008), doi:10.1063/1.2943700 | OSTI abstract only | M |
| Ol12 | R. E. Olson et al., Phys. Plasmas 19, 053301 (2012), doi:10.1063/1.4704795 | OSTI abstract only | M |
| Kl13 | J. L. Kline et al., Phys. Plasmas 20, 056314 (2013), doi:10.1063/1.4803907 | OSTI abstract only | M |

Hu24 and Kr24b were added in the addendum. Not read: Callahan et al., PoP 27, 072704 (2020) (Hu24's
source for the generic accounting), Jones et al.,
PoP 24, 056312 (2017) (the LASNEX hohlraum template and drive multipliers), Rosen et al., HEDP 7,
180 (2011) (DCA energetics), and Suter et al., PoP 3, 2057 (1996) (only the abstract, M).

## 2. Findings

### 2.1 Wall material: Au-lined DU and cocktails versus Au

- **Hammer–Rosen fits (Ro05 §VI, H).** For U₀.₆Dy₀.₂Au₀.₂, κ = 5670 ρ^0.10/T^0.90 cm²/g and e = 0.95 e_Au.
  For Au, κ = 6544 ρ^0.18/T^1.43 and e = 3.33 T^1.54/ρ^0.15 MJ/g.
  - Wall energy per area (hJ/mm²):
    - Au: (0.598, 0.398) T^3.3 t^0.6
    - cocktail: (0.604, 0.407) T^3.1 t^0.57
    - The pairs are for k = 0 (constant T) and k = 0.18 (T ∝ t^0.18).
  - At 270 eV and 1 ns the cocktail/Au loss ratio is (0.84, 0.85). Multigroup simulations give (0.85, 0.87).
  - The ratio scales as **T^−0.22 t^−0.05**, so a cocktail is worse than Au at low T. For a NIF
    ignition hohlraum, "the ratio is about 0.83".
  - U₀.₈₆Nb₀.₁₄ (the Omega cocktail): E/A = 0.416 T^3.1 against Au 0.39 T^3.3 (k = 0.18, 1 ns).
  - Oxygen contamination raises e: e_c/e_Au = 1.06 with no O, 1.08 at 4% O, 1.22 at 40% O.
    This erased the early Omega gains.
- **Cocktail at NIF scale (Sc07, H).** A 300 eV NIF drive with U:Dy:Au 60:20:20 (atom %) gives
  **17.5% less wall loss** (LASNEX/STA). That means **10% less laser energy** once LEH and capsule
  losses are included.
  - Omega data: about +6 eV at 275 eV.
  - LASNEX cocktail/Au loss ratio at 1 ns: 97.6% at 185 eV, 90.8% at 275 eV, 88.6% at 305 eV,
    confirming the T^−0.2 scaling.
  - Their 0-D model: (E_L − E_scat)·η_CE = E_wall + E_LEH, with
    - E_wall,Au = 0.39 A_w T₀^3.3 t^1.18
    - E_LEH = 0.58 A_LEH T₀⁴ t^1.72
    - **η_CE = 0.85 t^0.2**
    - T₀ is the drive at 1 ns (heV), for T ∝ t^0.18.
- **DU on NIF (Do15, H).**
  - Construction: 7.0 µm DU, lined inside with **0.7 µm Au** and backed by 22.3 µm Au.
  - "The first DU hohlraum experiments during the NIC ... quantified the drive enhancement from DU to
    be equivalent to an additional **6.5% of laser power (≈25 TW) at the peak**" (citing Kl13 and
    Callahan et al., PoP 19, 056305 (2012)).
  - Their power balance (eq. 1): P_w = σT⁴[(1−α_w)A_w + A_LEH + (1−α_cap)A_cap] − P_LEH − P_cap.
  - The Au-equivalent laser power is 1.065 × P_DU.
- **Kl13 (M).** DU hohlraums at 520 TW reached Trad ≈ 330 ± 4 eV. Laser coupling was 84 ± 3% over
  350–520 TW.
- **HYBRID-E (Kr22 §II, H).** Au-lined DU "require[s] less energy in the ... picket to achieve the same
  radiation temperature compared to unlined DU ... while providing a higher radiation drive in the
  peak ... compared to solid Gold ... due to the higher albedo." The liner thickness for HYBRID-E is
  not given.
- **Implied factor (L).**
  - The DU or cocktail wall loss is about 0.83–0.88 × Au at 270–310 eV. Taking the Ro05 ratio
    0.857 at 300 eV, 1 ns, and t^−0.03 to 8 ns gives **≈0.80–0.85 at HYBRID-E peak**.
  - In `WallLossAlbedo` terms, the coefficient would be 0.32 × ~0.83 ≈ 0.27, with the temperature
    exponent shifted by about −0.2. At 313 eV and 8 ns, 1 − α would go from 0.065 to about 0.054.
  - Do15's 6.5% peak-power equivalence is consistent in size with Sc07's "17% wall loss → 10% laser
    energy".

### 2.2 Laser-to-x-ray conversion efficiency η

- **Li04 §IV (H).**
  - η ≈ 70% in 1 ns Nova hohlraums, and possibly **as high as 85%** for longer shaped pulses (PS22, PS26).
  - The 1 ns flat-pulse Trad curve T ∝ t^0.15 "requires a conversion efficiency **0.7 t^0.12**" (Fig. 4-4 caption).
  - Hohlraum η exceeds disk η at the same intensity, because the re-emission region grows beyond the
    laser spot (Fig. 4-10). Disk η is more sensitive to intensity.
  - Power balance, eqs. 4-27 and 4-28: η_ce E_L = E_w + E_ca + E_h, with E_h = 10² T_r⁴ A_h τ and
    E_ca = 10² T_r⁴ A_c τ. The capsule is assumed to absorb its entire incident flux.
  - **A_c is "the initial area of the inside of the fuel layer"**, chosen as "representative of the
    effective area during the bulk of the energy absorption".
  - Coupling efficiency E_ca/E_L ≈ **10–20%** for case/capsule area 15–30, hole/capsule area 1–2,
    and η 70–80% (eqs. 4-29 to 4-34). Eq. 4-33 adds capsule re-emission through N_c, the ratio of
    incident to absorbed capsule flux.
- **Ro05 (H).**
  - Nova LASNEX: "CE = 0.7 t_ns^0.2. The efficiency increases with time in part because the albedo
    behind the conversion layer builds up with time." Balancing this source against E_w ∝ T^3.3 t^0.6
    gives T ∝ t^0.18.
  - Omega example: about 8% reflectivity and 68% CE at 1 ns.
- **Sc07 (H).** η = 0.85 t^0.2 (Omega, 1 ns).
- **De08 (M).** Measured on direct-drive Au, U and cocktail spheres. After about 0.5 ns the
  backscatter-corrected η reaches "a slowly increasing plateau of **95% at 10¹⁴ W/cm²** and of **80% at
  10¹⁵ W/cm²**". At 10¹⁵ the M-band (2–5 keV) is 10% of the flux for U and 27% for Au. LASNEX with
  f = 0.15 matched.
- **Ol12 (M).** NIF vacuum hohlraums show **η ≈ 85–90%** at 235 < Tr < 345 eV. This is explained by
  DCA emissivity plus f = 0.15, which reduce the coronal energy.
- **To11 (H).**
  - The coronal plasma stores "up to 40% of the energy" at NIF scale against 20% at Omega in the old
    XSN model. This energy is what the "high flux" (DCA + f = 0.15) model returns to x-rays.
  - The model matched peak Dante flux within the 5% error. Dante errors are ±5% at peak and ±10% in
    the foot.
- **BAH (H-img, p. 3).** 0-D model: η_ce f_laser P_laser = P_cap + P_wall + P_LEH, with
  - P_cap = σT⁴A_cap f_cap, f_cap = 1 − α_cap
  - P_wall = σT⁴A_wall(1 − α_wall)
  - Typical values for the 5.75 mm × 10.1 mm Au hohlraum:
    - **η_ce ≈ 0.95**
    - **f_cap ≈ 0.8**
    - A_wall ≈ 2.2 cm², A_LEH ≈ 0.15 cm², A_cap ≈ 0.04 cm² "at 2/3 radial convergence"
  - At 400 TW the model gives Trad ≈ 288 eV, α_wall ≈ 0.84, P_LEH ≈ 108 TW, P_wall ≈ 250 TW and
    P_cap ≈ 23 TW. The sensitivity is (1/P_cap) dP_cap/dα_w ≈ 4.5: **an α_w error of 0.02 is a 10%
    error in capsule drive**.
  - Average albedo, eq. 7, printed as ᾱ ≈ 1 − 0.44 T₀ τ^−0.46. The T₀ exponent is unclear in the
    render, so do not use it.
- **BAH VisRad calibration (H, Fig. 7 caption).** In a VisRad model of these Au hohlraums at a
  constant 400 TW, "X-ray conversion efficiency was set to **0.95**. The hohlraum wall albedo
  parameter was adjusted to **0.77** to match ... peak radiant intensity along the Dante-1 & Dante-2
  lines of sight."
  - This is the only published VisRad energetics calibration I found. It is the natural first
    benchmark for Briareus, which reproduces VisRad.
- **Picket versus peak (L).** η rises with time (t^0.12–0.2) and falls with intensity (95% → 80% from
  10¹⁴ to 10¹⁵ W/cm²).
  - In the picket (~50 TW, low intensity, < 1 ns) η is lower, because it is in its rise and the wall
    behind the conversion layer is cold.
  - At peak (~10¹⁵ W/cm² outer spots), the measured hohlraum-integrated η is 0.85–0.95.
  - The Briareus default of 0.75 is a Nova-era 1 ns number and is low for NIF peak.
  - The "hot spot" re-emission (the laser spot's emission region growing beyond the beam) is part of
    η in all of these 0-D treatments. None separates it from the Marshak wall.

### 2.3 Measured and simulated Trad for HYBRID-E

- **AS24 Table I (H-img).** "Drive Trad" (eV): N210207 **304**, N210307 **306**, N210808 **307**
  (averaged over similar shots), N220919 **313**, N221204 **313**.
  - Peak power (TW) and energy (MJ): 470/1.93, 487/1.91, 441/1.89, 440/2.05, 440/2.05.
  - SBS on N221204 was 42 ± 11 kJ and SRS 6 ± 1 kJ, so **98 ± 0.6% coupled**.
  - The longer pulse reached a 2% higher peak Trad "as expected given the increasing hohlraum wall
    albedo in time".
- **Kr22 Table I (H-img).**
  - N210808: Tr,max 307 (data) against 305 (simulated, "internal"). Simulated Tr at bang time 276 eV
    and at peak velocity 291 eV.
  - N210307 simulated: 294/254/270 eV.
  - The internal Trad at bang time was 20 eV higher for N210808 than N210307.
- **Pak24 Fig. 1 (H-img), Dante 1 Trad(t).** Measured through the LEH. Hardware failure: no N210808
  data, so "N210808-equivalent" repeat shots are shown.
  - Trad: ≈145 eV at 2–2.5 ns, ≈175 at 4 ns, ≈250 at 5 ns, ≈275 at 6 ns, ≈295 at 7 ns.
  - Peak ≈ **310–313 eV at ≈8.2 ns** for the N210808 pulse and ≈8.4 ns for N220919/N221204, falling
    to ≈285 eV at 9.2 ns.
  - "During the foot at low Trad the data has a large uncertainty". The LEH reduction "increased the
    radiation drive on the capsule by 7%".
- **Zy22 Fig. 2 (H-img), the cleanest LEH energetics pair** (symcaps, Dante):
  - N191118 (3.64 mm LEH, peak ≈ **465 TW**) peaks at ≈ **295 eV** near 7.3 ns.
  - N210601 (3.1 mm LEH, peak ≈ **390 TW**, 1.7 MJ) peaks at ≈ **300 eV** at 8.0 ns.
  - Foot: ≈150 eV at 2–3 ns. The two are "comparable until ~7.5 ns".
- **Ru24 (H).**
  - Across nine Hybrid-E shots: **TR = 313 ± 2 eV (Dante 1)** and **308 ± 3 eV (Dante 2**, flux
    corrected +28% for an LEH obstruction).
  - Dante 1 sits at θ = 143°, φ = 274°, with a **37° view of the LEH**. Dante 2 is at θ = 64°,
    φ = 350° with a shallower view. Uncertainties at peak are 4% in flux and 1.5% in TR.
  - Hohlraum energy from the Dante flux F (their eq. 1, Lambertian): E_H ∝ ∫F dt · [(1−α)A_wall + 2A_LEH]/A_LEH.
    The equation is garbled in the text, and it carries a cos θ factor.
  - The authors take a constant α = 0.8.
  - Post-bang reheating is not laser-driven: 313 → ~350 eV, E_H = 220 ± 30 kJ on N221204.
- **Kr24b Fig. 1(c) (H-img), simulated internal Trad.**
  - Trace values: ≈130 eV through the picket/trough (1–3.5 ns), ≈200 eV at 5 ns, ≈270 eV at 6–7 ns,
    peak ≈ 300 eV near 8.3 ns (N221204). The N210808 curve is nearly identical but ends ~0.3 ns earlier.
  - The internal trough is ~15–20 eV below Dante 1's ≈145–150 eV (L), and the peak is close to it.
  - Hu24 (H): "Trad ∼ 310 eV at peak levels". Kr24b (H): coupling "approximately 96%–99%".
  - Capsule (Kr24b Fig. 1a): outer radius 1135 µm, inner radius 1050 µm, 65 µm DT ice. Hohlraum
    radius 3.2 mm; LEH radius 1.55 mm.
- **What Dante sees (H).** Ch24 and BAH: through the LEH, Dante 1 sees mostly the outer-beam
  gold-bubble region and the inner wall. It sees a "smaller region dominated by outer beam-driven gold
  bubble", so the M-band fraction is higher than in an open-end view.
  - Ch24: LASNEX (f = 0.15) overpredicts the total flux by >10% for the LEH view and by <10% for the
    open-end view, and the M-band by 30–50%. A gold M-band opacity multiplier of 0.80 fixes it (0.87
    in BAH).
  - Kr22 lists the 307 eV data next to a 305 eV "internal" simulation. No paper I read gives a
    separate "Trad seen by the capsule" for HYBRID-E.

### 2.4 Capsule absorbed energy and albedo

- **Kr24 Table III (H).** Simulated capsule absorbed energy: **245 kJ** (N230729), **257** (+4 µm HDC),
  **274** (+10 µm), all at 2.05 MJ. The ablator mass remaining is 5, 5.9 and 6.6%.
- **Implied by G_capsule = Y_total/E_abs** (Kr22 and AS24 define E_abs as the "capsule x-ray absorbed
  energy"). These divisions are mine (**L** in the arithmetic, H in the inputs):
  - N210207: 0.17/0.75 ≈ **227 kJ**
  - N210307: 0.14/0.69 ≈ **203 kJ**
  - N210808: 1.3/5.8 ≈ **224 kJ** (AS24), or 1.37/5.8 ≈ 236 kJ (Kr22 data column)
  - N210808 simulated: 1.353/6.3 ≈ **215 kJ**
  - N221204: 3.1/12 ≈ **258 kJ**
  - That is **11–13% of laser energy**.
- **Capsule albedo (H).** BAH uses f_cap = 1 − α_cap ≈ 0.8, i.e. **α_cap ≈ 0.2**. Li04 assumes full
  absorption (α_cap = 0) but uses the area of the inner fuel surface, which is an area correction in
  place of an albedo.
  - I found no published time history of α_cap for W-doped HDC.
  - Ru24 notes that W-doped ablator debris reradiates after bang time. That is outside the drive window.
- **HYDRA energy balance for N221204 (Kr24b Appendix §3, H).** "Estimated using HYDRA 2D
  simulations calibrated to experimental measurements." Fractions of the incident laser energy:
  - **54.1%** heats the hohlraum wall (≈4.4% in the gold lining, ≈49.7% in the DU wall)
  - **25.2%** is radiation loss out of the LEHs
  - **4.5%** heats the He gas fill
  - **≈16.2% is absorbed by the capsule**
  - "In addition ... the laser energy needs to be degraded in simulations by 9%–13%."
  - "The useful energy coupled to the implosion (and not in the ablated capsule material) is
    approximately **246 kJ** or approximately 12%."
  - Fuel kinetic energy ≈15.3 kJ, total work on the DT ≈20 kJ (1%), and ≈8.2 kJ to hot-spot internal
    energy.
  - Notes: the four fractions sum to 100%, so wall heating includes the laser energy that is not
    converted to x-rays (no separate η is given). Whether "incident" means before or after the
    9–13% degradation is not stated. 16.2% is **332 kJ** of 2.05 MJ, or **289–302 kJ** of the
    degraded energy (L).
- **Hu24 (H).**
  - "The typical low gas-fill hohlraum energy accounting [Callahan et al., PoP 27, 072704 (2020)]
    is ∼50% of the unscattered laser energy ends up in the hohlraum wall, ∼30%–35% ... radiated back
    out the LEHs, ∼10%–15% is absorbed by the capsule (depends upon hohlraum to capsule surface
    area), and ∼5% ... absorbed into both the helium gas ... and the thin plastic membranes".
  - Table I caption: "for N221204 the capsule absorbed energy is Y_total/G_cap = 258 kJ", with
    G_cap 12 (N221204) and 5.8 (N210808).
- **What "capsule absorbed energy" means (H inputs, L reading).** Kr24b separates the total x-ray
  energy absorbed by the capsule (16.2%, ≈290–330 kJ) from the energy left in the unablated
  implosion (≈246 kJ). The quoted 245 kJ (Kr24, N230729), 246 kJ (Kr24b) and Y/G_cap = 258 kJ
  (Hu24, AS24) all match the second.
  - Neither paper states a time window. Kr24b's accounting is for the whole drive.
- **Consistency check, revised (L).**
  - Inputs: Pak24's Dante Trad(t) for N221204 (digitized), capsule R₀ = 1.135 mm held to 6 ns then
    shrinking linearly to 2/3 R₀ at 9 ns, and f_cap = 0.7–0.8.
  - Result: E_cap ≈ 316–361 kJ, which brackets HYDRA's **total** capsule absorption of 290–330 kJ.
  - The earlier "30–45% overshoot" against 245–258 kJ was a definition mismatch: about 25% of the
    absorbed energy leaves with the ablated material.
  - A Briareus capsule model with α_cap ≈ 0.2–0.3 on the outer area, driven by a Trad close to
    Dante's, is therefore consistent with HYDRA to about ±10%.
  - Using the internal Trad of Kr24b Fig. 1(c), which is ~15–20 eV lower in the foot, lowers the
    estimate slightly. The foot carries little energy.

### 2.5 LEH loss and closure

- **AS24 Appendix A (H).** Reducing the LEH from 3.65 to 3.1 mm gave "a 7% increase in radiation flux";
  P_LEH = σT_r⁴A_LEH.
- **Kr22 (H).** The smaller LEH "was used to operate at lower peak power ... while still achieving
  similar Trad" (Zy22 gives the numbers above).
- **Ch20 (H).**
  - The effective LEH radius from gated LEH images (3–5 keV) is ≈1550 µm before the high-density
    ridge forms. It drops quickly at the rise to peak, holds nearly constant at peak, and closes faster
    after the laser turns off.
  - At peak, the measurement is within about 50 µm of f = 0.15 LASNEX (0.6 mg/cc). At 0.3 mg/cc
    "data and simulations agree".
  - **VisRad models matched the LEH images only with the wall radius reduced by 300–700 µm**, because
    the x-ray conversion surface moves in.
- **BAH (H).** Simulations overpredict the LEH size early, by up to 100 µm at the highest power, which
  is ≈2% in capsule coupling. They match at peak. Closure modeling cannot explain the drive deficit.
- **Open-hole bound (L).**
  - For HYBRID-E, A_LEH = 2π(1.55 mm)² = 0.151 cm².
  - At 313 eV, σT⁴ = 985 TW/cm², so **P_LEH ≈ 149 TW, about 34% of 440 TW**.
  - Integrated over the digitized N221204 Trad(t) to 9.3 ns: **E_LEH ≈ 560 kJ (27% of 2.05 MJ)**.
  - This matches the 34.4% open-hole fraction already in research-notes.md.
- **Against HYDRA (Kr24b, H; comparison L).** HYDRA's LEH loss for N221204 is **25.2%** of incident
  laser energy, about 517 kJ of 2.05 MJ. My open-hole bound is ≈560 kJ, so the net effect of closure
  plus the slightly lower internal T is only ≈ −8%. Hu24's generic range is 30–35%.

### 2.6 Effective albedo and the drive deficit

- **BAH and Ch24 (H).**
  - Bang times are 400–700 ps late against simulation. Matching them takes 15–20% laser-power
    multipliers, set against ~2% laser-energy error.
  - BAH found no evidence for errors in LTE wall albedo, mix or LEH closure. The deficit is
    attributed to NLTE M-band overemission: 2–4 keV is 30% low in data, and the bubble Te is 31% higher
    than simulated.
  - The total x-ray drive is ~6% low (Dante 1).
- No published HYDRA effective-albedo table for low-fill DU hohlraums was found. The adoptable
  numbers are BAH's VisRad (η 0.95, α 0.77, Au, 400 TW) and BAH's 0-D α ≈ 0.84 at 288 eV.
- **Briareus 0-D cross-check (L).** HYBRID-E at peak: A_w = 2.75 cm², A_LEH = 0.151 cm², 440 TW × 0.98,
  313 eV.
  - η = 0.95, A_cap at 2/3 convergence (0.072 cm²), f_cap 0.8: 1 − α_w = **0.075** (α = 0.925).
  - η = 0.75, same capsule: 1 − α_w = 0.044 (α = 0.956).
  - η = 0.95 with the full initial A_cap (0.162 cm²): 0.049.
  - Lindl's gold formula at 313 eV and 8 ns: 0.065. With a DU factor of ~0.83: 0.054.
  - So **η ≈ 0.9–0.95 with U-wall albedo ≈ 0.94–0.95** closes the peak power balance only if the
    capsule's effective absorbing area is ≤ ~½ of its initial area.
  - Addendum: this used peak power alone. Over the whole pulse, Kr24b puts the capsule at 16.2%
    against 25.2% for the LEH, a ratio of 0.64. The open-hole areas give f_cap·A_cap/A_LEH = 0.8 ×
    0.162/0.151 = 0.86 at the initial capsule area. Averaged over the pulse, the effective f_cap·A
    is therefore ≈ 0.75 × the initial outer value (L), not ½.
- **Energy-level target (Kr24b, L).** Briareus's integrated output split should reproduce:
  - wall (laser deposition plus x-ray loss) ≈ 54%
  - LEH ≈ 25%
  - capsule ≈ 16%
  - gas ≈ 4.5% (Briareus has no gas absorber; book it as a laser loss)
  - ...with the laser reduced by 9–13% as HYDRA does. That reduction is in effect a drive multiplier
    (BAH).

## 3. Validation targets (N210808 / N221204 benchmark)

| Quantity | N210808 | N221204 | Source | Conf. |
|---|---|---|---|---|
| Laser energy / peak power | 1.89 (1.917) MJ / 441 TW | 2.05 MJ / 440 TW | AS24 T.I (Kr22 T.I) | H |
| Coupled fraction (1 − backscatter) | "energetically low" | 98 ± 0.6% (SBS 42 kJ, SRS 6 kJ) | Kr24 App.; AS24 | H |
| Peak Dante 1 Trad | 307 eV (avg of repeats) | 313 eV (±2 eV over 9 shots; 1.5% 1σ) | AS24 T.I; Ru24 | H |
| Peak Dante 2 Trad | — | 308 ± 3 eV (after 28% flux correction) | Ru24 | H |
| Simulated internal Trad | max 305, 291 at peak velocity, 276 at BT | — | Kr22 T.I | H-img |
| Trad, end of picket / trough | ≈145–150 eV at 2–3 ns (large uncertainty) | same | Pak24 Fig. 1; Zy22 Fig. 2 | H-img |
| Trad at rise | ≈175 eV at 4 ns, ≈250 at 5 ns, ≈275 at 6 ns | similar, ~0.1–0.2 ns later | Pak24 Fig. 1 | H-img |
| Time of peak Trad | ≈8.2 ns | ≈8.4 ns | Pak24 Fig. 1 | H-img |
| LEH-size pair at fixed Trad | 3.64 mm, 465 TW → 295 eV; 3.1 mm, 390 TW → 300 eV (symcaps N191118, N210601) | — | Zy22 Fig. 2 | H-img |
| Simulated internal Trad (N221204) | — | ≈130 eV trough, ≈270 eV at 6–7 ns, peak ≈300 eV at ≈8.3 ns | Kr24b Fig. 1(c) | H-img |
| HYDRA energy split (fraction of incident laser) | generic: ~50% wall, 30–35% LEH, 10–15% capsule, ~5% gas+windows (Hu24) | **wall 54.1% (Au 4.4%, DU 49.7%), LEH 25.2%, gas 4.5%, capsule 16.2%**; laser degraded 9–13% | Kr24b App. §3; Hu24 | H |
| Total capsule x-ray absorption | — | 16.2% ≈ 332 kJ (≈289–302 kJ of the degraded laser) | Kr24b | H / L |
| "Capsule absorbed energy" (unablated capsule) | ≈215 kJ sim (1.353 MJ/6.3); Y/G_cap ≈224–236 kJ | ≈246 kJ (Kr24b), Y/G_cap = 258 kJ (Hu24); 245 kJ N230729 (Kr24) | Kr24b; Hu24; Kr22; Kr24 | H / L |
| Capsule share of laser | ≈12% (unablated) | 12% unablated, 16.2% total | Kr24b; derived | H / L |
| Fuel kinetic energy / work on DT / hot-spot internal energy | — | ≈15.3 kJ / ≈20 kJ / ≈8.2 kJ | Kr24b | H |
| LEH loss | open-hole bound ≈34% of peak power | HYDRA 25.2% (≈517 kJ); open-hole bound ≈27% (≈560 kJ) | Kr24b; 0-D | H / L |
| Wall share | — | 54.1% (laser deposition plus x-ray loss) | Kr24b | H |
| VisRad calibration (Au, 5.75 mm, 400 TW) | η 0.95, α_w 0.77 match Dante 1 and 2 peak intensity | — | BAH Fig. 7 | H |

Dante geometry for synthetic diagnostics:
- Dante 1: θ = 143°, φ = 274°, 37° to the hohlraum axis through the lower LEH.
- Dante 2: θ = 64°, φ = 350°.
- Ch24 gives 37.4°. The 0.8 keV to few-keV LEH size differs, so use the open-hole geometry with Ch20's
  radii.

## 4. What Briareus should implement (minimal)

1. **Wall-material loss factor.** Add `loss_factor` (default 1) and `loss_temperature_exponent`
   (default 0) to `WallLossAlbedo`, giving 1 − α = f c (T/100)^(p+δp) t^q, and a `uranium_lined()`
   preset.
   - Preset values: f ≈ 0.83 and δp ≈ −0.2, so the ratio to Au follows T^−0.2. Source: Ro05, Sc07.
   - Optionally keep the Au value until `HeatFront` mass exceeds the liner areal mass,
     ρ_Au·d = 19.3 × 0.7 µm = 1.35e-3 g/cm² (Do15 liner). By `HeatFront` this happens within
     ~0.5 ns at 150 eV, so it affects only the picket, consistent with Kr22 (L).
   - Validate with the Do15 Au/DU pair (6.5% Au-equivalent peak power) and the Kl13 330 eV
     (520 TW, DU) point.
2. **Time- and intensity-dependent conversion.** Make `conversion` accept a model:
   - η(t) = η₀ (t/1 ns)^s, clamped. Use s = 0.12–0.2 (Li04, Ro05, Sc07).
   - η(I) per element, interpolated in log I between (10¹⁴, 0.95) and (10¹⁵ W/cm², 0.80) (De08).
   - Raise the NIF gold default toward 0.85–0.9 at peak (Ol12, BAH).
   - Keep it generic in the number type; the absorbed flux per element is already known.
   - Validate with the BAH VisRad pair (η 0.95, α 0.77 → Dante 1/2 peak) and the Li04 1 ns Nova
     T ∝ t^0.15 curve.
3. **Capsule model.** Expose `ablator(albedo)` as a `TabulatedAlbedo`, default constant 0.2 per BAH,
   plus an effective-radius schedule R_cap(t) (initial R, then R₀(1 − (1/3)·ramp) toward 2/3 near peak
   velocity).
   - Report E_cap = ∫(1−α_cap)·(incident on capsule) dt in `energy_balance`.
   - Calibrate the **total** absorbed energy against Kr24b's 16.2% (≈290–330 kJ for N221204).
     Do not calibrate against the 245–258 kJ "capsule absorbed" (unablated) values: they exclude
     ≈25% that leaves with the ablated material.
   - Validate the whole split as well: wall 54%, LEH 25%, capsule 16%, gas 4.5%, with a 9–13%
     laser degradation.
4. **LEH effective-radius schedule.** Add a `TabulatedAperture` or a time-dependent scale on aperture
   area, built from Ch20/BAH-type GLEH data: open size in the foot, about 100–300 µm smaller at peak,
   closing after laser-off.
   - Also consider an inward wall-radius offset of 300–700 µm for the emitting surface (Ch20, which
     VisRad needed).
   - Validate with the Zy22 3.64 mm and 3.1 mm symcap pair: same Trad at 465 and 390 TW.
5. **Synthetic Dante.** Use `line_of_sight_power` or `OrthoCamera` from θ = 143°, φ = 274° through
   the lower LEH, report T_Dante = (F π / (σ A_LEH,proj))^¼, and compare with 313 ± 2 eV and the
   Pak24 Trad(t).
   - Report alongside it the capsule-incident T. This separates "Trad" definitions, which today
     differ silently.

## 5. Open questions

- What is the HYBRID-E Au liner thickness, and is the wall fully DU behind it? (Do15's 0.7/7.0/22.3 µm
  is for 2014 high-foot.)
- Partly settled by Kr24b: "capsule absorbed" (≈246 kJ) excludes the ablated material, and total
  capsule absorption is 16.2%. Still open:
  - Are Kr24b's fractions of the delivered laser or of the 9–13% degraded laser?
  - Is the 16.2% net of capsule re-emission? (It reads as net absorbed.)
  - What is the integration window?
- No W-doped HDC α_cap(t) is published. The total-energy check supports α_cap ≈ 0.2–0.3 on the outer
  area to ±10%. The time history is unconstrained.
- What Trad does the capsule see, against Dante 1's 37° LEH view? At peak, simulated internal and
  measured agree to within a few eV (Kr22: 305 against 307; Kr24b: ≈300 against 313 for N221204).
  In the trough, the internal ≈130 eV sits below Dante's ≈145–150 eV.
- How does Kr24b's 54.1% wall heating split between unconverted laser energy and x-ray Marshak loss?
  This would fix η directly. Callahan et al., PoP 27, 072704 (2020), the source of Hu24's generic
  accounting, may give it; it was not read.
- How to fold the 15–20% drive multipliers and NLTE M-band overemission (BAH, Ch24) into a gray
  radiosity model? A gray Briareus calibrated to Dante (as BAH's VisRad was) would absorb them into
  η and α. That can be checked against bang time only through a capsule model.
