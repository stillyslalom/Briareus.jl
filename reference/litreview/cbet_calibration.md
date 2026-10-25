# CBET calibration data for low-fill NIF hohlraums: what a reduced model needs

Scope: the data and physics `CBET` / `CalibratedTransfer` need to be predictive for
0.3–0.6 mg/cc He hohlraums. Confidence: **H** = primary text read (incl. LLNL preprints on
OSTI), **M** = abstract/press/secondary, **L** = inference. "Digitized" = read off a
figure by eye (±10%), so **L** for the number even when the figure is H.
Δλ is at 1ω (inner − outer) unless stated. Items already in `reference/research-notes.md`
(Michel 2009 conditions, Nguyen gain, Pickworth 24 µm/Å, HYBRID-E 40–50 µm/Å) are only
cited where the table needs them.

## Sources read

| Key | Citation | Access | Conf. |
|---|---|---|---|
| K18 | A. L. Kritcher et al., PRE 98, 053206 (2018), doi:10.1103/PhysRevE.98.053206 | OSTI 1597608 (LLNL-JRNL-753089) | H |
| P20 | L. A. Pickworth et al., PoP 27, 102702 (2020), doi:10.1063/5.0004866 | OSTI 1769081 (LLNL-JRNL-804567) | H |
| HE21 | A. L. Kritcher et al., "HYBRID-E: Record hot spot energies…", PoP 28 (2021) | OSTI 1818404 (LLNL-JRNL-818898; preprint, some numbers garbled) | H |
| K22 | A. L. Kritcher et al., PRE 106, 025201 (2022), doi:10.1103/PhysRevE.106.025201 | journals.aps.org (open) | H |
| K24E | A. L. Kritcher et al., PRE 109, 025204 (2024), doi:10.1103/PhysRevE.109.025204 | CHORUS accepted manuscript | H |
| K24P | A. L. Kritcher et al., PoP 31, 070502 (2024), doi:10.1063/5.0210904 | given | H |
| I23 | N. Izumi et al., "Control of low-mode drive asymmetry in an efficient long pulse low-density filled hohlraum", PoP (2023) | OSTI 1971298 (LLNL-JRNL-842205) | H |
| S17 | D. J. Strozzi et al., PRL 118, 025002 (2017), doi:10.1103/PhysRevLett.118.025002 | arXiv:1607.06523 | H |
| M13 | P. Michel et al., PoP 20, 056308 (2013) | OSTI 1212855 | H |
| N21 | K. L. Nguyen et al., "CBET saturation by ion trapping-induced detuning" | arXiv:2104.07725 | H |
| B16 | M. A. Barrios et al., PoP 23, 056307 (2016) | OSTI 1842431 | H |
| NC24 | "Impact of low-mode symmetry on inertial fusion energy output…", Nat. Commun. 15 (2024), doi:10.1038/s41467-024-47302-8 | OSTI 2405189 | H |
| M26 | P. Michel et al., "Laser polarization effects on CBET in ICF", PoP (2026), doi:10.1063/5.0335309 | LLNL/ANS press only | M |
| R21 | J. E. Ralph, A. Kemp et al., PoP (2021), doi:10.1063/5.0044404 (Ta₂O₅ liner, SBS) | LLNL press only | M |
| Kw14 | R. Kirkwood et al., "Laser backscatter damage risk assessments…", APS DPP 2014 NP8.118 | abstract via search snippet | M |

Not obtained: Zylstra PRE 106, 025202; Pak PRE 109, 025203; the PRL 132, 065102
supplement; any NIF Thomson-scattering measurement in a low-fill LEH (none found).

## 1. Measured transfer vs Δλ in low fill

**K18 (H)**: CH, 0.6 mg/cc He, 6.72 × 11.24 mm, 3.64 mm LEH. 1ω wavelengths
(Table I): N161006 23°/30° 10530.5, outer 10528.5 Å (Δλ = 2); N161103 10529.8/10524.3
(5.5); N170706 23° 10530.5, 30° 10529.8, outer 10524.3 (6.2/5.5, "three-color").
Outers delayed 1.2 ns. Transfer is *inferred*: HYDRA (flux limiter 0.15, high-flux
model) computes CBET with full power, which is then imposed with measured backscatter
removed; foot multipliers came from no-CBET shots.
- Foot: clamp δn/n = 1×10⁻² ("effectively unsaturated", since CBET waves are
  10⁻⁴–10⁻³) matched VISAR shock symmetry at 2 and 5.5 Å with **no** symmetry multiplier.
- Rise + peak, N161103: clamps 5×10⁻⁴ to 2×10⁻³ match the final shock mergers, giving
  **1.36–1.6×** the average inner-beam power.
- N170706 in-flight P2 (radius ≈ 200 µm) needs clamp (2–3)×10⁻³, giving **1.67 (1.9)×** the
  inner power. HiFoot N140520 (1.6 mg/cc, 6.2 Å) gave 1.26×.
- N150826 (Δλ = 0, same platform, peak CF 34% vs 35%) is oblate in flight and at the hot spot.
- Predicted for a 1100 µm CH capsule at 1.8 MJ: hot-spot P2 −3 to +4.5 µm at 1 Å and
  +19 to +22 µm at 2 Å (clamp (2–2.5)×10⁻³).

**P20 (H)**: CH, 0.6 mg/cc, 4.0 mm LEH, 470 TW, CF 33% at peak (11–13% picket).
- ΔP2 between 0 and 1 Å: 27 µm at ~1.25 MJ, 23 µm at 1.56 MJ, 24 ± 3 µm/Å overall,
  i.e. CF 33 → 41.5% (+26% waist drive) via 2.7 µm per %CF.
- Picket VISAR pair (N180522 at 0 Å, N180412 at 1 Å): equator shock 5% faster, so the
  outers give **6% more** of their picket energy to the inners (outer share
  0.87 → 0.82).
- HYDRA used clamp 0.01, with peak cone-fraction multipliers 0.86 (0 Å) and 0.80 (1 Å).
- SBS fluence 0.1–0.19 J/cm² vs a damage threshold of "~1 J/cm²".
- 0.5 Å was predicted round.

**I23 (H)**: CH, 0.45 mg/cc, 90% scale, 14 ns, ~280 TW, peak CF 0.44.
- Core a2 goes from −17 ± 3 µm (N211129, 0 Å) to ~0 at 1.5 Å (N220301); **14 µm/Å**,
  "1.5–2× less" than other low-fill work.
- The authors attribute the difference to 30–40% lower peak power ("recall CBET transfer
  fraction ~ I").
- Simulations give ΔCF = +0.08 (+18%) for +20 µm between 0 and 1.5 Å.
- Foot CF lowered by 0.02 per Å to cancel picket CBET.

**HE21 (H, preprint)**: HDC, 0.3 mg/cc ("first … with Δλ > 0"), CCR ≈ 2.7.
- Sensitivity ~40–50 µm P2/Å; 1–2 Å suffices.
- Shrinking the LEH 4.0 → 3.64 mm made shapes more oblate; +0.5 Å recovered them,
  ending at 1.75 Å (1100 µm).
- "No sign of 'saturation' for CBET."
- A small P2 correction per 20 TW of peak power accounts for more transfer at lower
  power and cooler, denser LEH plasma. The number is illegible in the OSTI text.
- Peak CF multiplier 1.065 in simulations.
- The 44° cone gives more than the 50° "due to positioning"; there is some
  inner–inner and outer–outer transfer.

**K22 (H)**: N210207/N210307 at 1.55 Å (3.64 mm LEH, 470/487 TW); N210808 at 1.8 Å
(3.1 mm LEH, 441 TW). Simulated mode-2 flux asymmetry was multiplied by **~0.7 in the
peak** to match the data.

**K24E/K24P (H)**: N220919 2.5 Å, hot-spot P2 −12 µm (K24E; K24P Table II: −10.5);
N221204 2.75 Å, ~0; N230729 2.9 Å, +6 µm (thicker ablator, foot +30 ps, so confounded).
Implied local slope ≈ 40–48 µm/Å (*derived*, L).

## 2. Plasma state at the crossing region (Briareus holds it constant)

No published *measured* time history of ne, Te, Ti or flow in a low-fill LEH was found.
What exists:
- **M13 (H)**, LEH at the start of the main pulse:
  - Fully ionized He; ne = 0.03 nc, Te = 2.8 keV, Ti = 0.8 keV.
  - Quad intensities 5×10¹⁴ (inner) and 10¹⁵ W/cm² (outer), Δλ = 6 Å.
  - The 276 quad-pair beat waves reach |δn|/n ≈ 10⁻³.
  - Stochastic ion heating runs at 6–7 keV/ns in a ~1 mm³ LEH region and slows after
    ~0.5 ns. In C₅H₁₂, C ions pass 10 keV and H stays below 4 keV after 1 ns.
  - Lasnex at peak: Ti ≈ 4 keV in the Au plasma near the LEH, ~1 keV in the He where
    the beams overlap.
  - Flow is near sonic at the LEH and supersonic outside; outer→inner transfer
    accelerates it.
  - Once heated, the gain is **linear in Δλ** (phase velocity < v_Ti).
  - High-fill era.
- **S17 (H)**: off resonance, g ∝ (λᵢ−λ₀) Z ne Tᵢ^{1/2}/(Tᵢ + Z Te)², so hotter
  electrons cut CBET. Inline SRS Langmuir heating raises LEH Te and reduces CBET; LW
  power "is effectively outer-beam power" (deposited near the LEH).
- **B16 (H)**: Mn/Co dot spectroscopy, 0.8-scale target, C₅H₁₂ fill equivalent to
  1.6 mg/cc He, ~2 mm from TCC. Peak Te is 4.2 ± 0.2 keV, 0.5–0.6 keV above HYDRA; Te
  follows the drive and peaks late.
- **HE21 (H)**: lower power implies a cooler, denser LEH and more transfer.
  *Inference (L)*: in 0.3 mg/cc the crossing plasma is ablated Au/C, not He, early in
  the pulse. Briareus' CH default is an assumption, not a measurement.

## 3. Where transfer happens

- **P20 (H, HYDRA)**, at 9.5 ns:
  - *Outside* the LEH the transfer is inner → outer, entering the 50° cone; the 44°
    cone barely takes part. *Inside* it is outer → inner.
  - The net depends on the balance, "throughout the pulse".
- **NC24 (H)**: I-Raum needed 0.5 Å "mostly to negate the flow-induced CBET", and
  HYBRID-E 1–1.55 Å. So Δλ = 0 is **not** zero transfer in a low-fill hohlraum.
  Both numbers are the Δλ that *minimized hot-spot P2* ("In experiments, 1 ≤ Δλ ≤ 1.55 Å
  was required to minimize P2 in HYBRID-E designs"), not the transfer null: HYBRID-E's
  Δλ also compensates bubble absorption of the inner beams, so its transfer null lies
  below 1 Å; I-Raum's 0.5 Å (little bubble absorption) bounds the null from above
  (*L* for HYBRID-E). The same paper's pairs are steep and unsaturated: N201101 → N210307,
  1.75 → 1.55 Å, P2 15.7 ± 0.7 → 0.68 ± 2.8 µm (≈ 75 µm/Å).
- **K24E Fig. 10 (H, digitized L)**, N221204 per cone, input → post-CBET (TW):
  - 23°: 62 → ~125 at the start of the peak, ~105 at 7.5 ns.
  - 30°: 62 → ~118 / ~100, then ~80 at 8.5 ns.
  - 44°: 125 → ~55–60 at 6 ns, ~78 at 8 ns.
  - 50°: 125 → ~95 at 6 ns, ~105 at 8 ns.
  - The 44° cone loses ≈ 2× the fraction the 50° does. Absolute TW are as plotted; read
    the ratios.
- **I23 (H)**: at 1 Å the 30° spots get 4× brighter in x-rays while the 23° spots do
  not, so the transfer goes preferentially to 30°. Within-cone pairing and
  pointing/bubble shadowing matter.
- **M26 (M)**: CBET *within* a cone gives beam-to-beam power spread under linear
  polarization; circular polarization keeps the mean transfer and cuts the spread,
  which is relevant to SBS on the amplified beams.
- **Δλ as a power proxy (K24E, H)**: without CBET, reaching the same symmetry by
  power alone would need the outers cut to 150–200 TW (from 440 TW).
- Beam deflection and spot displacement from CBET: nothing quantitative found; none of
  the papers invoke it for P2.

## 4. Saturation: clamp values and physics

| Era / design | Fill | Clamp δn/n | Phase | Source |
|---|---|---|---|---|
| CH 2016–17 | 0.6 | 1×10⁻² | foot | K18 H |
| CH 2016–17 | 0.6 | 5×10⁻⁴–2×10⁻³ (shocks); 2–3×10⁻³ (in-flight) | rise+peak | K18 H |
| CH 2018 | 0.6 | 1×10⁻² all pulse, plus CFM 0.80–0.86 | all | P20 H |
| N121130 Lasnex inline SRS+CBET | 1.45 | 1×10⁻² (results same for 0.004–0.1) | all | S17 H |
| N221204 | 0.3 | 0.002 foot; 0.0075 (K24E) / 0.008 (K24P) rise+peak | split at rise | H |

- **Physical basis.** Ion trapping and two-ion-wave decay set in near δn/n ~ 10⁻² (S17,
  H). Below that, the operative saturation is detuning by stochastic ion heating (M13)
  and trapping-induced frequency shifts (N21, H: the gain peaks at only ~80% of linear
  in the small-v_p case, then detunes over > 20 ps).
- **How the clamp enters.** S17 uses a hard limit δn = min(δn_sat, (p/α)^{1/2}), so the
  coupling is min(I₀Iᵢ, β√(I₀Iᵢ)). Briareus uses 1/√(1+(δn/c)²), which is √2 lower at
  δn = c (L).
- **Ignition-era clamps are not limiters.** Linear δn/n is 10⁻⁴–10⁻³ (K18, M13), so
  0.0075–0.01 is nearly inactive and 0.002 binds only on the strongest pairs (L). Their
  role is to set a *phase-dependent effective gain*. K24E retuned the foot clamp from
  N221106 shock data ("modifying the CBET saturation") and verified it on N230110.
- **Clamp or multiplier?** The low-fill evidence (linear P2 vs Δλ, no roll-off, a slope
  that scales with power: HE21, I23, P20) favors a gain multiplier plus a correct
  plasma state over a clamp (L).
- **The cliff.** K24P: 2.75 → 2.9 Å overshot by ~6 µm. Candidates (L):
  - Pump depletion: the 44° cone is already at ~45–50% of input at the start of the
    peak, so transfer saturates in power, not δn.
  - Nonlinear inner SBS: I23 saw ×5.2 SBS for +36% power.
  - Confounding by the 30 ps foot change and the thicker ablator.

## 5. Backscatter and the SBS limit

- **K22 (H)**: on N210808, inner SBS was ~0–1 kJ per cone (23° up to 16 kJ error bar);
  SRS 2.8 ± 1.1 kJ (N210207: 0.74 ± 0.28).
  - 50° SBS: 4.8 ± 2.4 kJ (N210808) vs 9.5 ± 3.4 kJ (N210207), reduced because more
    outer power was transferred. 44° SBS was similar to 50°; outer SRS < 1 kJ.
  - Totals: 13.4 kJ of 1.92 MJ (0.7%) for N210808 and 15.5 kJ of 1.93 MJ for N210207,
    vs 216 kJ of 1.76 MJ for N140520 (1.6 mg/cc).
- **K24E (H)**: backscatter on N210808/N221204 is "energetically small", mostly inner
  SBS early in the peak.
  - It is not modeled explicitly; it is lumped into a total power multiplier and δn/n.
  - Simulations degrade the laser by 9–13%.
  - Energy balance: 54.1% to the wall (4.4% in the Au liner, 49.7% in the DU), 25.2%
    LEH loss, 4.5% gas, 16.2% capsule.
  - N220919's Δλ had "a maximum limit … due to concerns with laser damage from …
    SBS".
- **Damage criterion.** The only number found is "~1 J/cm²" SBS fluence (P20, H). HE21
  puts 50° SBS fluence ~13× below threshold (H). Risk review compares a proposed shot
  with a database of past shots plus extrapolation rules and asks for walk-up shots
  (Kw14, R21, M). The N220919 numbers are not public.
- **I23 (H)**: SBS + SRS ≈ 3% (0.45 mg/cc). 30° SBS rises and 50° SBS falls with Δλ.

## 6. Calibration dataset

CF = inner/total power. "Inferred" = post-processed by HYDRA, not measured.

| Design / shot | Fill, LEH | Δλ (Å) | Phase | Transfer measure | Hot-spot / in-flight P2 | Src |
|---|---|---|---|---|---|---|
| N150826 CH | 0.6, 3.64 | 0 | peak | — | oblate (in flight + HS) | K18 H |
| N161006 CH | 0.6, 3.64 | 2.0 | foot | shocks matched, clamp 0.01 | — | K18 H |
| N161103 CH | 0.6, 3.64 | 5.5 | rise+peak | inner ×1.36–1.6 (inferred) | shock symmetry | K18 H |
| N170706 CH | 0.6, 3.64 | 6.2/5.5 | rise+peak | inner ×1.67 (1.9) (inferred) | in-flight P2 matched by clamp 2–3e-3 | K18 H |
| N180522/N180412 CH | 0.6, 4.0 | 0 / 1 | picket | outer picket share 0.87 → 0.82 (+6%) | shock Δv 5% | P20 H |
| CH symmetry pairs | 0.6, 4.0 | 0 / 1 | integrated | CF_eff 33 → 41.5% | ΔP2 = 27 (1.25 MJ), 23 (1.56 MJ); 24±3 µm/Å | P20 H |
| N211129 → N220301 CH | 0.45, 3.64 | 0 → 1.5 | peak | sim ΔCF +0.08 (0.44 → 0.52) | a2 −17 → ~0; 14 µm/Å | I23 H |
| HYBRID-E HDC | 0.3, 4.0/3.64 | 1–1.75 | peak | — | 40–50 µm/Å; LEH 4.0 → 3.64 needs +0.5 Å | HE21 H |
| N210207, N210307 | 0.3, 3.64 | 1.55 | — | — | burning plasma, ~round | K22 H |
| N210808 | 0.3, 3.1 | 1.8 | — | sim P2 flux ×0.7 at peak | ~0 µm | K22, K24E H |
| N220919 | 0.3, 3.1 | 2.5 | — | SBS-limited | −12 (K24E) / −10.5 (K24P) µm | H |
| N221204 | 0.3, 3.1 | 2.75 | foot ~2–4 ns | CF 0.13–0.15 → 0.17–0.20 (sim, digitized) | ~0 µm | K24E H/L |
| N221204 | " | 2.75 | early peak ~5.5 ns | CF 0.34 → ~0.64; 23° ×2.0, 30° ×1.9, 44° ×0.45, 50° ×0.76 | " | K24E L |
| N221204 | " | 2.75 | late peak ~8 ns | CF 0.34 → ~0.52; 44° ×0.62, 50° ×0.84 | " | K24E L |
| N230729 | 0.3, 3.1 | 2.9 | — | — | +6 µm (prolate) | K24P H |

*Derived (L)*: N221204 moves ≈ 45% of outer power to the inners at the start of the peak,
≈ 27% late in the peak, and ≈ 5–7% in the foot. Via 2.7 µm/%CF the sensitivity targets
are dCF_eff/dΔλ ≈ 0.15–0.18 per Å (0.3 mg/cc, ~440–480 TW), 0.085 (0.6 mg/cc, 470 TW)
and 0.05 (0.45 mg/cc, 280 TW).

## What Briareus should implement

`CBET` already accepts time functions or tables for every plasma field, `multiplier` and
`clamp`. The gaps are defaults, structure and tests:

1. **A phase schedule preset**, e.g. `cbet_schedule(pulse; rise)`: clamp 0.002 before
   the start of the rise and 0.0075 after (K24E), driven by the same `PulseShape` used in
   `drive_history`. Document it as an HYDRA-equivalent setting, not physics.
2. **A hard-min saturation option** (`saturation = :min | :smooth`) matching S17's
   min(δn_sat, linear). Then HYDRA clamp values transfer one to one.
3. **Two crossing zones per LEH**: an outside slab (supersonic outflow, Mach > 1, where
   inner → outer occurs at small Δλ) and an inside slab (Mach ≲ 1), each with its own
   length and state. Regression test: Δλ = 0 gives net inner → outer, and an I-Raum-like
   case needs Δλ ≈ 0.5 Å to null it (NC24). Magnitude L.
4. **Plasma state tables as calibration inputs, not defaults**: ne/nc(t), Te(t), Ti(t),
   u(t), with a He/Au/C composition option. Ship M13's (0.03, 2.8, 0.8 keV, He) as a
   documented "main-pulse start, high fill" point and mark all low-fill values L. Add
   an optional ion-heating relaxation (Ti rising at ≤ 6–7 keV/ns until c_s exceeds the
   beat-wave phase velocities, M13). That makes the gain linear in Δλ, as the data show.
5. **Multi-epoch, multi-shot `calibrate_cbet`**: fit `multiplier` (optionally per phase)
   jointly to (Δλ, time, per-cone post/input ratio) rows. Reproduction tests with
   target numbers:
   - N221204 (2.75 Å): post-CBET CF ≈ 0.18 ± 0.03 (foot), 0.64 ± 0.06 (5.5 ns),
     0.52 ± 0.05 (8 ns); 44° ratio below the 50° ratio at the peak.
   - Slope dCF_eff/dΔλ in 0.15–0.18 per Å for HYBRID-E at 1–1.8 Å, linear within 10%.
   - Intensity scaling: the slope falls by 1.5–2× when peak power drops 30–40% (I23).
   - Then convert to hot-spot P2 at 2.7 µm/%CF and check 40–50 µm/Å.
6. **Report outer-cone depletion** (post/input per cone) and a post-CBET inner-quad SBS
   proxy (intensity × time against a user threshold). This flags the N220919-type limit
   and the 2.9 Å cliff, where linear extrapolation fails.
7. **`Backscatter` presets**: low-fill inner SBS during the early peak at ≲ 1% of
   energy; N210808 total 0.7%. Note that ignition-era HYDRA folds it into power
   multipliers, so do not double count with a calibrated `multiplier`.

## Open questions

- No measured low-fill LEH ne, Te, Ti or flow history exists in the open literature
  (NIF 5ω OTS reports, LLNL-TR-858245, not checked). Briareus' constant state cannot be
  validated, only calibrated.
- Is the post-CBET CF of K24E Fig. 9/10 the HYDRA *inline* value at the LEH or an
  effective one? Fig. 10's absolute cone powers sum to ~374 TW against 440 TW peak.
- HE21's power-correction coefficient (µm per 20 TW) is illegible in the preprint;
  check the published PoP 28 version.
- How was N220919's "maximum Δλ" set (SBS fluence model, database rule)? Kirkwood 2014
  is abstract-only.
- Does the 2.75 → 2.9 Å overshoot come from CBET nonlinearity, outer pump depletion,
  or the confounded foot/ablator changes? A Δλ scan at fixed design is not public.
- Do CBET-driven spot displacement or beam deflection matter for P2? Nothing found.
- How should polarization handle within-cone spread (M26), given that Briareus
  averages with (1+cos²θ)/4?
