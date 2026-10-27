# Inner-cone SBS time histories and the late-time refraction: literature (2026-10-25)

ROADMAP "The P2 swing" and "Absolute energetics on N221204". `examples/n221204.jl` removes the measured N221204 backscatter
(AS24: SBS 42 ± 11 kJ "early in peak power and mostly on the inner cones", SRS 6 ± 1 kJ
"throughout peak") from the inner cones on time shapes of our own (*L*). This review asks
what is published about the per-cone, time-resolved backscatter of HYBRID-E shots and
about the late-time refraction of the inner beams, and whether either can be applied.
Results of applying them: research notes, "Inner-cone SBS time histories and the
late-time refraction (2026-10-25)".

Confidence: **H** primary text read; **H-img** read from a figure raster (digitized,
errors stated); **M** abstract only; **L** our inference. Cite-only papers: numbers with
citations, no text or figures reproduced.

## 1. Sources

| Key | Reference | Access | Used for |
|---|---|---|---|
| K26 | A. J. Kemp, T. Chapman, L. Divol, D. P. Higginson, E. Kur, N. Lemos, S. MacLaren, P. Michel, D. J. Strozzi, G. B. Zimmerman, "Modeling stimulated Brillouin backscatter from the inner laser cones during indirect-drive inertial confinement fusion experiments at the National Ignition Facility", PoP 33, 082705 (2026), doi:10.1063/5.0320730 | full text: accepted manuscript LLNL-JRNL-2014353, OSTI 3441608 | N210808 per-cone SBS histories (Fig. 2), SSH vs MSH (Fig. 4), the on-axis density spike (Fig. 6), axis/equator balance (Fig. 8) |
| K25s | Same authors, AAC 2025 talk (Sedona, 2025-05-11), LLNL-CFPRES-2006194, OSTI 3388898 (doi:10.2172/3388898) | full slides | per-time reflectivities, "2× incoming power" without the spike, N210404 |
| H25 | D. P. Higginson et al., "Direct evidence of multispecies hydrodynamics in ignition-scale hohlraums", PRL 134, 165101 (2025), doi:10.1103/PhysRevLett.134.165101 | full text: LLNL-JRNL-871162, OSTI 2571728 | the bubble–axis stagnation measured with DT gas (N230313-002) |
| I24 | N. Izumi et al., "Neutron imaging of the deuterium–tritium tamping gas volume in an inertial confinement fusion hohlraum", RSI 95, 103510 (2024), doi:10.1063/5.0219564 | full text: OSTI 2476138 | the imager; the shot's pulse (power halved after 6.4 ns) |
| K24o | A. J. Kemp, M. Belyaev, N. Lemos, T. Chapman, L. Divol, E. Kur, P. Michel, "Modeling stimulated Brillouin backscatter from outer-cone quads across multiple inertial confinement fusion hohlraum designs", PoP 31, 042707 (2024), doi:10.1063/5.0179693 (LLNL-JRNL-853952; open access, CC BY-NC-ND) | full text (publisher PDF, supplied by the user) | outer 50° SBS per platform (Figs. 1, 3) and N210808's Q36B history (Fig. 6) |
| LLNL26 | LLNL news, "LLNL researchers improve predictive modeling of laser backscatter" (lmf.llnl.gov, 2026) | web | summary of K26 |

## 2. Inner-cone SBS: what is measured (K26 §II, Fig. 2; **H**, **H-img**)

- Shot: N210808 (first ignition, HYBRID-E, 0.3 mg/cc He). The measured curves are
  cone-averaged (both hemispheres, all azimuths) DrD and FABS signals combined by a
  Gaussian-process / hierarchical statistical model (the paper's appendix), in-lens power
  × 3.4 (lens fraction f_lens = 1/3.4 from past NBI calibrations): all backscattered SBS
  light per quad.
- **Energies**: 2.8 kJ per quad (23°) and 0.9 kJ per quad (30°); eight quads per cone, so
  ~22 + 7 ≈ 30 kJ of inner SBS on N210808 (**L**, arithmetic). Single-shot error "up to
  100 %" (lens fraction and calibration); N210808 repeats vary by 30–50 %.
- **Time history** (digitized, `kemp2026_sbs()`, ±0.03 TW, ±0.02 ns): both cones start
  as the inner power reaches its peak (23° at 4.6 ns with a sharp rise the authors call
  a possible processing artifact; 30° at 4.8 ns). The pre-CBET inner power per quad on
  the same axis is ~9 TW, half of it at 4.74 ns. 23°: peak 2.18 TW at 5.3 ns, 1.0 TW at
  6.2, 0.15 TW at 6.55, zero by 6.95 ns. 30°: a shoulder of 0.41 TW at 5.1, peak 0.90 TW
  at 5.75, 0.46 at 6.2, zero by 6.7–7.0 ns. FWHM ~1.4 ns (23°, 4.8–6.15 ns) and ~0.8 ns (30°, 5.4–6.2 ns, above its 0.4 TW shoulder); the
  inner power stays at peak to ~8.8 ns, so inner SBS covers only the first ~2 ns of a
  ~4 ns peak.
- Peak instantaneous fraction of the pre-CBET quad power: 23° 24 %, 30° 10 % (**L**).
- Axis- vs equator-side beams of a quad: peak in-lens SBS ratio ~10 for both cones on
  HYBRID-E (Fig. 8a: the 30° in-lens total 0.26 TW at 5.75 ns, nearly all axis side).
- FABS (Q31B, 30°) spectra were not recorded on N210808; N211031 is used instead.
- N221204 itself: not in K26. Its only published backscatter is AS24's 42 ± 11 kJ SBS
  and 6 ± 1 kJ SRS (already in the example).

## 3. Inner-cone SBS: what is modeled (K26 §III; K25s; **H**, **H-img**)

- pF3D on frozen 2D LASNEX (LHT) plasma maps, 4–6 snapshots per shot, post-CBET quad
  power (inline CBET, |δn/n| clamp 3.6 × 10⁻³) as the boundary condition; the post-shot
  hydro matches drive with wall-opacity multipliers rather than the 85 % laser power of
  earlier work (K24o).
- **Single-species hydro (SSH)** over-predicts: 30° quad 2.42 TW at 6.0 ns against the
  measured 0.90 at 5.75 (Fig. 4); reflectivities R = 0.08 / 0.21 / 0.003 at 5 / 6 / 7 ns
  (K25s). "Seemingly no plausible adjustments to CBET" fix it (K25s); the CBET power
  differs by < 5 % between SSH and MSH.
- **Multispecies hydro (MSH**, 13-moment, Pearson-IV closure): He interpenetrates the Au
  bubble over tens of µm, Landau damping of the ion-acoustic waves rises, SBS falls
  exponentially with the He fraction. 30° quad: R = 0.06 / 0.09 / 10⁻⁶ at 5 / 6 / 7 ns
  (K25s); powers (Fig. 2) 23°: 0.55 (5.0), 1.18 (5.25), 1.70 (5.5), 2.00 (5.75), 1.62
  (6.0), 1.04 (6.5), 0 (7.0) TW; 30°: 0.48 (5.0), 1.20 (5.5), 1.06 (6.0), 0.79 (6.5),
  0.22 (6.75), 0 (7.0) TW. Energies 2.8 (23°, as measured) and 1.7 kJ per quad (30°,
  1.9× measured). A "hybrid" run (SSH densities and temperatures, MSH ion fractions)
  behaves like MSH: mixing, not the background plasma, sets the level.
- **Duration.** The fast drop shortly after peak power is "common across all
  ignition-relevant designs with a 0.3 mg/cc fill" and is reproduced by SSH and MSH
  alike: it is not mixing but the plasma filling the inner-beam path. The model's 23°
  onset trails the measurement by ~0.3 ns (DrD timing up to 0.2 ns, plus the bubble's
  expansion into the inner beams).
- **Spectra**: pF3D puts SBS from the Au wall near the capsule at 5 ns (spectral center
  +2.5 Å), Au bubble plus ablator at 5.5–6 ns (double-humped), ablator dominated by 6 ns
  (+7 Å); FABS does not show the 7–8 Å ablator peak at 6–7 ns. Hypothesis: the on-axis
  density feature forms ~1 ns earlier than the hydro predicts; adding it at 6 ns removes
  the 7–8 Å peak and reduces the SBS.

## 4. The late-time on-axis density spike and the refraction (K26 §III C, Fig. 6; K25s; H25)

- **Where and when** (MSH map of N210808, Fig. 6a): at ~7 ns the plasmas from the
  capsule ablator, the LEH window and the wall meet near the axis at z ≈ ±4.5 mm (LEH
  at ±5.62 mm). The electron density on axis there (inset) is n_e/n_c ≈ 0.01 to 6.0 ns,
  0.02 at 6.5, 0.08 at 6.75, 0.15 at 7.0 and 0.14 at 7.75 ns; similar in SSH and MSH.
- **What it does**: it "enhances inverse-bremsstrahlung absorption and causes refraction
  of inner-cone beams, reducing the incident laser power in the region where SBS
  backscatter originates". The 30° path passes closer to the spike, so its SBS drops
  ~0.25 ns before the 23°'s. Halving the spike's density at 7 ns raises the 30°
  reflectivity by four orders of magnitude. The talk (K25s, slide 11) titles it "inner
  cone SBS turns off when incident light is refracted into the Au bubble / absorbed by
  on-axis plasma" and gives the only number on the incident light: "smoothing over" the
  7 ns feature gives 2× the incoming power into the ablator.
- **Magnitude for the drive** (**L**): so at 7 ns about half of the 30° power that would
  reach the ablator region (between capsule and wall, where the beam meets the bubble) is
  refracted into the gold bubble or absorbed on axis. The transverse spectra (Fig. 6c, e)
  broaden from the quad's two lobes (k_y/k₀ ≈ ±0.03) to roughly twice that by 7 ns:
  deflections of a few degrees, i.e. a few tenths of a mm over the remaining few mm of
  path, smaller than the spot. So the published picture is a reduction of the inner power
  reaching the ablator/waist region (absorbed near the axis, or deposited in the bubble,
  closer to the LEH), not a displaced spot. No fraction for the 23° cone, no energy, and
  no effect on drive symmetry is published.
- **Timing caveat**: the spike rises over 6.5–7.0 ns in the hydro, while the spectra
  suggest it forms ~1 ns earlier in the experiment (K26 §III D, conclusions: "a more
  rapid Au bubble expansion in the underlying hydro model" would bring spectra and SBS
  history into agreement).
- **Measurement** (H25, I24): N230313-002 (N210404 design, DT at 0.3 mg/cc in the
  hohlraum, 2 mg/cc ⁴He in the capsule; laser 1.345 MJ, power cut to 1 TW per beam after
  6.15 ns): neutrons come mostly from the bubble–axis region at z ≈ ±3.4 mm (FWHM
  1.2–1.35 mm), burn peaking at 8.15 ± 0.05 ns (FWHM 1.86 ns); MSH reproduces yield
  (1.01 against 1.36 ± 0.1 × 10¹³) and shape, SSH over-predicts 3.5×. At 8 ns the Au
  compresses the DT on axis to 4.3–5.6 mg/cc (14–19×), i.e. n_e/n_c ≈ 0.12–0.15 at 3ω
  (**L**), the spike's density. Its timing on a full-power N221204 pulse is not measured.

## 5. Outer cones (K24o; K26; **H**, **H-img**)

- **What is measured.** 50° outer quads only (FABS on Q36B, DrDs elsewhere, the same
  cone-averaging statistics as K26; the 44° quads are not studied). SBS energy *into the
  lenses* of one 50° quad (Fig. 3, read by eye, ±10 J): HYBRID-E N210808 ~50 J
  (pF3D 40 J); N190218 ~10, N190320 ~55, N200125 ~105, N191110 ~230 (2019–2021 HYBRID-E
  shots; pF3D 140, 180, 145, 300 J); BigFoot/Hybrid-B N170702 ~125, N170601 ~265,
  N180429 ~175, N180319 ~540, N180128 ~715 J (pF3D 25, 85, 215, 310, 410). Errors ~30 %
  (simulation sampling) and the measured spread around the cone; model and measurement
  agree to within a factor ~2.
- **Peak reflectivity** (Fig. 1, R_peak = max over time of P_SBS / P_net, P_net the
  post-CBET quad power from HYDRA's CBET model): HYBRID-E 0.012–0.08 measured (pF3D
  0.02–0.05), Hybrid-B 0.05–0.16, BigFoot 0.04–0.35 (pF3D 0.1–0.4); the historical fall
  came from design changes and from Δλ moving power off the outer cones.
- **N210808's time history** (Fig. 6, one quad, Q36B, FABS; DrDs were noise-dominated):
  nothing above the noise (~0.02–0.05 TW) before ~6.5 ns; a single burst rising from
  ~6.8 ns to 0.20 TW at ~7.4 ns and gone by ~7.9 ns, i.e. late in the peak, after the
  inner SBS has ended. Integrated 6–8 ns: 170 J per quad, 50 J into the lens (×1/3.4).
  pF3D peaks ~0.5 ns earlier (R < 10⁻⁴ at 6 ns, 0.024 at 7 ns, from the Au LEH liner's
  expansion at 7 ns) and gives 40 J into the lens. On a ~9 TW pre-CBET quad (Fig. 4a:
  the outer cones' 32 quads carry ~300 TW), 0.2 TW is ~2 % at peak (**L**).
- **Scale** (**L**): 16 quads at 50° × 170 J ≈ 2.7 kJ of outer SBS on N210808 (unknown at
  44°; ≤ 5.4 kJ if the same), against ~30 kJ on the inner cones: ~10–20 %, consistent
  with AS24's "mostly on the inner cones" for N221204.
- K26: outer-cone SBS originates deeper in the Au bubble and MSH changes it less.
- **Use in Briareus** (**L**): the outer history is one quad on one shot, at one of the two
  outer cones, with the noise floor at a quarter of its peak; it is not cone-averaged and
  N221204's is unpublished. `backscatter=history` therefore stays inner-only. Its
  effect would be small: ~2 % of the 50° power for ~0.7 ns near 7.4 ns raises the
  post-CBET cone fraction by < 1 point there, ≲ 0.003 in flux P2/P0 by the
  −0.29 %/% sensitivity, after the swing's zero crossing.

## 6. Bearing on Briareus

1. The published inner-cone history is **early and short**: N210808's inner SBS starts
   with the peak, is largest in its first ~0.7 ns and is gone ~2 ns into a ~4 ns peak.
   AS24's "early in peak power" for N221204 says the same. The example's `measured`
   shape (ramp 5.0–5.5, flat to 6.5, gone by 7.5 ns, equal on both cones) puts it later,
   flatter and longer, and splits it equally between the cones where the 23° cone
   carries 3.1× the 30°'s. `backscatter=history` applies the published curves.
2. Kemp's conclusion that single-species models over-predict the magnitude and duration
   concerns pF3D on hydro maps. Briareus does not model SBS; it removes measured
   backscatter. What matters here is the measured history, now published per cone.
3. The refraction is a reduction of the inner power reaching the ablator/waist region at
   ~7 ns (hydro), possibly ~6 ns (spectra). Briareus's `BubbleAbsorption` already
   absorbs a rising share of the inner light in the bubble (0.22 → 0.57 over 5.5–6.5 ns
   on HYBRID-E) before the ridge stop: the spike's "refraction into the Au bubble" is
   partly the same physics, so a separate loss would double count. Tried as a bound
   (`refraction`, research notes): booked as lost with the spike ~1 ns early, it moves
   the crossing earlier and lowers the late drive, away from Fig. 10.
4. Neither paper gives N221204's per-cone SBS. Applying N210808's history to N221204
   needs a time map (rise-aligned, *L*) and an energy (AS24's 42 kJ).
