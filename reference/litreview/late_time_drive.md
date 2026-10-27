# Late-time hohlraum drive and Dante after the laser falls: literature review (2026-10-25)

ROADMAP, "Absolute energetics on N221204". After the laser power falls, Briareus's drive falls much faster than Dante 1
on all three benchmark shots. Model minus measured over the initial LEH area, at 8.5 ns:
N191118 −73 eV, N210601 −44 eV, N221204 −28 eV. At bang time Dante reads 226 / 254 eV on
the Zylstra pair, the model 155 / 186 eV (research notes, "Zylstra 2022 LEH pair,
time-dependent"). Closing this takes ~250–450 kJ more x-ray source after the fall. This
review asks which published mechanism could supply it, how large each is, and which of
them Briareus's data can test.

Confidence: **H** primary text read; **H-img** read from a rendered figure (±5 eV,
±10 %); **M** abstract or secondary; **L** our arithmetic or inference. The papers are
cite-only. Downloaded texts and rendered pages stayed in the session scratchpad.

## 1. Sources read

| Key | Reference | Access | Used for |
|---|---|---|---|
| Ru24 | M. S. Rubery et al., "Hohlraum reheating from burning NIF implosions", PRL 132, 065104 (2024), doi:10.1103/PhysRevLett.132.065104 (CHORUS accepted manuscript; LLNL-JRNL-845705) | full text | Dante conventions, impulse response, cooling-hohlraum fits, reheating (W-doped detritus) |
| Ru23 | M. S. Rubery et al., "Soft x-ray power diagnostics for fusion experiments at NIF, Omega and Z", RSI 94, 031101 (2023) (OSTI 1975922, LLNL-JRNL-841387) | full text | cable Wiener deconvolution, XRD impulse response 140 ps |
| BAH | G. F. Swadling et al., "Resolving discrepancies in bang-time predictions ... Build-A-Hohlraum campaign", arXiv:2501.10350 (2025) | full text (preprint) | Dante 1 view filled by the gold bubble; LASNEX falls faster than Dante after laser-off; initial-LEH-area convention |
| Ch24 | H. Chen et al., "Understanding the deficiency in ICF hohlraum x-ray flux predictions ...", PRE 110, L013201 (2024), doi:10.1103/PhysRevE.110.L013201 (also OSTI 2429360) | full text | NLTE M-band over-emission of the bubble; GLEH after laser-off within 2 % |
| Ca20 | D. A. Callahan et al., "A simple model to scope out parameter space for indirect drive designs on NIF", PoP 27, 072704 (2020), doi:10.1063/5.0006217 (OSTI 1785922) | full text | exponential cooling term after laser-off; 0-D fails in the cooling phase; Fig. 4 Dante falls |
| LHL23 | J. D. Lindl, S. W. Haan, O. L. Landen, "Impact of hohlraum cooling on ignition metrics for inertial fusion implosions", PoP 30, 012705 (2023), doi:10.1063/5.0113138 (OSTI 1961756, LLNL-JRNL-833944) | full text | cooling rate set by laser fall, LEH, wall; "typical 1.5 ns" drive-pressure cooling time; HYBRID-E "slow cool" via smaller LEH and slower turn-off |
| Kr22 | A. L. Kritcher et al., PRE 106, 025201 (2022), doi:10.1103/PhysRevE.106.025201 | full text | HYDRA internal Tr at peak, peak velocity and bang (Table I) |
| Zy22 | A. B. Zylstra et al., PRE 106, 025202 (2022), doi:10.1103/PhysRevE.106.025202 | full text | the LEH pair's pulses and Dante (already in `zylstra2022_fig2()`) |
| Pa24 | A. Pak et al., PRE 109, 025203 (2024), doi:10.1103/PhysRevE.109.025203 | full text | N221204 pulse and Dante 1 (already in `n221204_pulse()`) |
| To11 | R. P. J. Town et al., "Analysis of the NIF ignition hohlraum energetics experiments", PoP 18, 056302 (2011) (OSTI 1018446, LLNL-JRNL-463439) | full text | corona stores "up to 40 %" at NIF scale (XSN); Fig. 2 Dante after laser-off vs the two models |
| Kl11 | J. L. Kline et al., "Observation of high soft x-ray drive in large-scale hohlraums at NIF", PRL 106, 085003 (2011), doi:10.1103/PhysRevLett.106.085003 | full text | P_corona = (1 − η)(P_L − P_bs); corona 15–30 % of the energy on NIF; η 85–90 % |
| Ro05 | M. D. Rosen, "Fundamentals of ICF hohlraums" (2005), OSTI 885380 | full text | ablative-wall equations and EOS (P = r e/V, r = 0.25), Hammer–Rosen coefficients |
| Ch20 | H. Chen et al., PoP 27, 022702 (2020) (OSTI 1632379) | full text (read before) | GLEH radius after laser-off |
| Mo17 | J. D. Moody et al., "Semi-empirical 'leaky-bucket' model of laser-driven x-ray cavities", PoP 24, 042709 (2017), doi:10.1063/1.4981221 | **abstract only (closed)** | the cooling term Ca20 adopts |
| RR22 | M. D. Rosen, M. S. Rubery, "Understanding the enhanced hohlraum drive that follows high yield", APS DPP 2022 BP11.127 | abstract | N210808 Dante rise after bang |

Not found in open form: Moody 2017 (the leaky-bucket cooling constant), Jones et al.
PoP 24, 056312 (2017) (OSTI 1375645 has no full text), the Aybar et al. DPP 2022 abstract
on the rate of rise (BP11.126, abstract text not served).

## 2. Mechanisms

### 2.1 What the published codes and data say about the fall itself

- **HYDRA keeps the hohlraum hot after the laser falls (Kr22 Table I, H-img).** Simulated
  internal Tr (the volume average over 1.1–1.5 capsule radii, Ca20's definition) at its
  maximum / at peak velocity / at bang: N210207 295 / 279 / 260 eV (coast 1.05 ns, 3.64 mm
  LEH); N210307 294 / 270 / 254 (coast 1.2 ns, 3.64 mm); N210808 305 / 291 / 276 (coast
  0.86 ns, 3.1 mm). Peak to bang: −35, −40, −29 eV. Dante's maxima (302, 294, 307 eV) match
  the internal maxima. On the pair the data fall −69 eV (N191118, coast 1.72 ns) and −47 eV
  (N210601, coast 1.23 ns) from peak to bang; the model's capsule Tr falls −127 and −94 eV
  (coarse, current defaults).
  So the slow fall is in the hohlraum, as HYDRA computes it, not only in Dante's view (*L*).
  Caveat: late in the implosion the 1.1–1.5 r_cap shell lies inside the capsule's own
  hot blow-off, so this "internal" Tr includes the exhaust's emission near the capsule.
- **LASNEX falls faster than Dante after laser-off (BAH, H).** "There is evidence that
  simulations ... overestimate closure speed once the laser drive turns off. This may help
  to explain why x-ray drive rises earlier in simulations and then falls more quickly after
  the lasers turn off." The LEH measured after the drive ends is larger than simulated.
  Ch24 (H): GLEH after laser-off agrees with LASNEX within 2 % in diameter. So the LASNEX
  shortfall is a view effect of a few %, of the opposite sign to ours (Briareus's lips do not
  close after laser-off, Ch20 and "Computed LEH closure").
- **The older LASNEX model, with more energy in the corona, fits the fall better (To11 Fig.
  2, H-img).** On N091204 (1.0 MJ) the Dante flux after the laser peak (19 ns) follows the
  data with the XSN, f = 0.05 ("Nova") model, while the DCA, f = 0.15 "high-flux" model,
  which matches the peak, falls ~20–25 % below the data in flux by 20.5–21 ns (~6 % in T,
  ~15 eV, *L*). To11's text: the corona stores "up to 40 %" of the energy at NIF scale in the
  XSN model; the more emissive model stores less. Kl11 (H): "absorbed laser power not
  converted into x-rays is retained in the hot, under-dense coronal plasma",
  P_corona = (1 − η)(P_L − P_bs), 15–30 % of the energy on NIF against 5–10 % on Nova;
  DCA with f = 0.15 halves the coronal specific energy. Kl11 Fig. 3 (H-img): the 2009
  vacuum-hohlraum data also sit above the DCA runs after laser-off at 150–635 kJ.
  This is the only published evidence on the mechanism, and it is qualitative: energy held
  in the corona during the pulse comes back after it.
- **0-D models fail exactly here (Ca20, H).** "The main disagreement between the model and
  the Dante data happens when the laser power drops and the hohlraum starts to cool ...
  During the cooling phase, the radiation lost out of the laser-entrance-hole becomes more
  important and is not accounted for in this simple model." Ca20's eq. 9 replaces the
  physics after laser-off by T_r^3.3(t) = T_r^3.3(t_laser) exp[−3.3 (t − t_laser)/τ_decay]
  (following Mo17); the τ_decay value is not printed. Its Fig. 4 Dante falls ~30–35 eV/ns
  after the peak (290 → 195 eV over ~3 ns on an HDC/DU shot: T e-folds in ~7.5 ns, *L*).
- **Cooling rate (LHL23, H).** "The hohlraum drive temperature drops at a rate that depends
  on the rate at which the laser power drops, the size of the laser entrance holes, and the
  hohlraum wall area and opacity"; the drive-pressure cooling time is "typical[ly] 1.5 ns";
  HYBRID-E's "slower cooling hohlraum was developed by reducing the LEH size and by turning
  off the laser more slowly" (Fig. 1f, HYDRA, 305 eV). No mechanism beyond the wall and LEH
  is named.

### 2.2 Dante conventions and instrument

- **Initial LEH area (Ru24, BAH, H).** Tr is computed from the flux "using the
  Stefan-Boltzmann relation combined with the initial LEH dimensions ... gold expansion
  throughout the drive will reduce the effective LEH area; TR values presented herein
  therefore represent a lower bound" (Ru24). BAH eq. 12 uses the initial LEH diameter for
  data and simulation alike. Briareus's `dante(...).Tr` follows this convention; `Tr_open`
  is the upper reading.
- **Time response.** Ru24: channel impulse response "well-approximated by a ~400 ps wide
  Gaussian, dominated by the bandwidth of the signal cable and XRD"; cross-timing ±50 ps.
  Ru23: the cable response is removed by Wiener deconvolution, and the XRD impulse response
  is ~140 ps FWHM. Applied to the model's Dante 1 history (standard mesh), a 140 ps Gaussian
  moves it ≤ 1 eV and a 400 ps one ≤ 5 eV anywhere in 5–9 ns (N210601 +5.0 eV at 9.0 ns,
  N191118 +3.5 at 8.5) (*L*). **Ruled out** as the cause.
- **Cooling hohlraum fits (Ru24, H).** To isolate reheating Ru24 fits "a power law to the
  flux before and after the onset of reheating" from 8 ns: the flux of a cooling hohlraum is
  smooth, and the reheating peak is ~500 ps after bang. So reheating does not enter Dante
  before bang (below).
- **Unfold after laser-off.** No published statement that UNSPEC is treated differently
  after laser-off. Ru24's reheating spectrum at 10 ns is "clearly thermal". No published
  Dante spectra between laser-off and bang were found, so whether the late signal carries a
  hot (M-band or bremsstrahlung) component is unknown.

### 2.3 Gold bubble and LEH plasma in Dante's view

- **The bubble fills Dante 1's view (BAH, H).** "As the experiment proceeds the gold bubble
  expands out to almost fill the Dante-1 line of sight." The bubble is ~3 keV (BAH), Te
  3.8–4.7 keV at peak in LASNEX with MHD (Ch24); LASNEX over-emits its M-band by 30–50 %
  (Ch24, BAH). Its emission after laser-off is not published.
- **Magnitude (L).** Its thermal energy is small: n_e ≈ 0.1 n_c (9e20 cm⁻³), 3.5 keV,
  ~0.08 cm³ (two annuli near the outer spots) gives ~60 kJ, which conduction (f v_e n_e T at
  f = 0.15, ~3 mm/ns) hands to the wall in ~0.3 ns. That is the 0.15 ns source lag already
  tested ("N221204 late-time drive": ~60 kJ, 10–15 eV at 8.5 ns). The bubble's kinetic
  energy is larger and comes back only where it stagnates (on the axis, against the capsule
  blow-off and the gas); no published power or timing.
- **LEH lips.** Chen 2020's GLEH data hold ~0.86 r₀ after laser-off; the model's open-area
  reading (lips as bright as the wall behind them) already bounds a glowing lip and is
  still 58 eV low at 8.5 ns on N191118. **Not sufficient.**

### 2.4 The capsule's ablated corona (W-doped HDC)

- **Reheating after bang (Ru24, H).** On shots with yield the remnant W-doped ablator,
  heated by the exploding core ("detritus drive"), reheats the walls; the Dante signal peaks
  ~500 ps after bang; E_H = 220 ± 30 kJ on N221204, 52 ± 12 kJ on N210808 (LASNEX 42 kJ);
  E_H ∝ Y^1.21. Prompt core x-rays < 5 keV are absorbed by the W-doped ablator; > 5 keV
  deposit too deep in the gold to show on Dante. The core's later stagnation on the gold
  is a several-ns effect (13 ns). Before bang none of this operates: it is not the cause of
  the 7.5–9.4 ns deficit, and N191118 / N210601 are symcaps with no yield.
- **Before bang (L).** The ablated HDC carries most of the ~300 kJ the capsule absorbs
  (16.2 % of 2.05 MJ, Kr24b; ~15 kJ ends as fuel kinetic energy): roughly 3.8 mg of carbon
  ablated, ~100 kJ thermal at 300 eV (Z = 6) and ~150–250 kJ kinetic. With W at 0.2–0.4 at.%
  it radiates efficiently, but it is heated by the drive and cooled by expansion, so it is
  a net emitter only where it is hotter than the radiation (stagnation against the gold or
  gas, or after the drive falls). Its emission after laser-off is not published. HYDRA's
  ablator albedo rises to 0.3–0.4 late in the drive (Young 2020, "Wall re-emission and the
  capsule radius"); after the fall it should exceed 1. Briareus's capsule share is
  13.9 % at α_c = 0 against HYDRA's 16.2 %: the missing ~2 points (~40 kJ on N221204) are
  absorbed outside the ablation front, plausibly in this blow-off (research notes,
  "Energetics defaults and the capsule share"), which is the same material whose stored
  energy would come back after the fall. **Unquantified.**

### 2.5 Hohlraum gas

He 0.3 mg/cc, 4.5 % of the laser energy (Kr24b): 30–60 kJ at 2.5–5 keV, returned by
conduction over ~2–3 ns ("Computed LEH closure"). ≲ 12 TW, ~1–2 eV (*L*).

### 2.6 The wall itself

No published closure for a gold wall under a cut-off drive was found (Ro05 and Hammer &
Rosen treat a driven Marshak wave). Briareus's closure (`wall_release`, `release` 0.77)
was calibrated on a static wall for 5–30 % drops. §3 tests the wall against independent
1-D references, including an ablating wall with Ro05's equations, inside a 0-D hohlraum
on each shot's measured pulse.

## 3. Tests with Briareus's data

Details and tables: research notes, "Late-time drive: the wall, the plasma and Dante
(2026-10-25)"; script `examples/late_time_wall.jl`.

| Mechanism | Published magnitude | Bears on Dante over the initial / open area | Test here | Result | Conf. |
|---|---|---|---|---|---|
| Dante time response | XRD 140 ps FWHM after cable deconvolution (Ru23); "~400 ps" Gaussian (Ru24); timing ±50 ps | same on both readings | convolve the model's Dante 1 history | ≤ 1 eV (140 ps), ≤ 5 eV (400 ps) | H / L |
| Initial vs open LEH area | initial area by convention; a lower bound (Ru24) | open reading is the upper bound | already reported (`Tr`, `Tr_open`) | open reading still 58 eV low at 8.5 ns (N191118) | H |
| LEH closing after laser-off | GLEH holds ~0.86 r₀ (Ch20); LASNEX closes too fast (BAH); GLEH within 2 % (Ch24) | initial-area reading darkens with closure | "Computed LEH closure" | ruled out (+0.8 eV) | H |
| Gold wall re-emission (static) | none published for a cut-off | both | 0-D hohlraum, 1-D static diffusion | falls 43 % (N191118, 7→9 ns) against Dante's 19 % | L |
| Gold wall re-emission (ablating) | E/A = 0.58 T^3.34 t^0.59 hJ/mm² (Ro05) | both | 0-D hohlraum, 1-D Lagrangian rad-hydro | falls 35 %; gives back 30–45 kJ in the ns after the fall of 620–820 kJ stored | L |
| Release factor (closure) | — | both | prescribed drops and the 0-D hohlraum on the ablating wall; full model with 1.0 | 0.9–1.07 needed against 0.77 (static); with 1.0 the full model's Dante deficit at bang shrinks −75 → −43 eV (N191118), −71 → −48 (N210601) | L, calibrated |
| Coronal (bubble) thermal energy | P_corona = (1 − η)P_abs retained (Kl11); 15–30 % (Kl11), up to 40 % in XSN (To11) of the energy | Dante 1's view is filled by the bubble (BAH) | 0.15 ns source lag ("Computed LEH closure", ~60 kJ) | 10–15 eV at 8.5 ns; thermal part ~60 kJ, conducted out in ~0.3 ns | H / L |
| Coronal kinetic energy, stagnation | part of the 15–40 %; no return rate published | Dante 1 sees the axis region through the LEH | not testable | — | — |
| Capsule exhaust (W-doped HDC) | ~300 kJ absorbed (Kr24b), ~100 kJ thermal (L); W at 0.2–0.4 at.%; emission after laser-off unpublished | capsule-region emission seen through the LEH | not testable | — | L |
| Detritus reheating | 220 ± 30 kJ (N221204), peak ~0.5 ns after bang (Ru24) | both | timing | after bang only; not in the 7.5–9.4 ns window, absent on the symcaps | H |
| He gas | 4.5 % of the laser energy (Kr24b), 30–60 kJ, conduction ~2–3 ns | both | estimate | ~1–2 eV | L |

The shot data in hand (N221204, N191118, N210601: delivered power and Dante 1 Tr, plus
HYDRA's internal Tr for the sister shots N210207/N210307/N210808) cannot separate the
plasma reservoirs: they only constrain their sum and its time scale, and fitting either
to Dante would make Dante a fit target.

## 4. Conclusions and what would settle it

1. **The gold wall cannot hold the drive.** In a 0-D hohlraum on each shot's measured
   pulse, neither a static nor an ablating 1-D gold wall (Ro05's equations) keeps T up
   after the laser falls: the fall is ~2× Dante's in fraction. HYDRA does keep it up
   (Kr22: −29 to −40 eV from maximum to bang). So the missing ~250–450 kJ is in the
   hohlraum plasma, which a radiosity model has no store for.
2. **The closure is not the main error, but its release factor is low for an ablating
   wall**: ~1.0 rather than 0.77 (0.9–1.07 over 10–50 % drops and in the 0-D hohlraum).
   In the full model it recovers 2–27 eV at 8.5–9 ns and 17–33 eV at bang: 30–45 % of the
   gap at bang, less at 8.5 ns.
3. **The published evidence points to energy held in the laser-heated corona** (Kl11's
   definition of 1 − η; To11 Fig. 2: the model storing more energy in the corona follows
   Dante's fall, the high-flux model falls ~20–25 % low in flux). Its size, 15–40 % of the
   laser energy (Kl11, To11) or 1 − η = 10–25 % with the measured η 0.85–0.90 (Kl11), is the
   right order (160–400 kJ on 1.6–2 MJ). Its return rate after laser-off is not published;
   only its thermal part (~60 kJ, ~0.3 ns) is estimable, and that is the lag already
   tested.
4. **The capsule's W-doped exhaust** (the user's suspect) holds ~100 kJ thermal; its
   emission before bang is unpublished. Its reheating role after bang (Ru24) is
   established but lies outside the window.

What would settle it, in order of value:
- HYDRA's (or LASNEX's) **energy inventory vs time** through the coast: radiation,
  wall internal + kinetic, Au corona/bubble internal + kinetic, capsule exhaust, gas, at
  laser-off and at bang. One table for N210808 or N221204 would fix both the reservoir
  size and its return rate.
- A **post-shot synthetic Dante 1** from HYDRA for N191118 or N210601 (to compare with
  HYDRA's own internal Tr and the data), or the internal Tr histories of Kr22 Fig. 6a
  read through the coast.
- **Dante spectra between laser-off and bang** (UNSPEC color temperature against
  brightness temperature, the >1.8 keV fraction): a hot plasma source in the view shows
  up as excess M-band or a hard tail; a cooling wall stays Planckian at Tr.
- **Mo17's τ_decay values** (closed access) and Callahan 2020's fitted τ for the
  HYBRID-E shots: a published empirical cooling rate to check against, not fit to.
- **Laser-off LEH and bubble images** (GLEH, SXI gated) to bound the bubble's emitting
  area and temperature in Dante 1's view after the laser falls.
