# Drive asymmetry → implosion observables: published transfer functions

Scope: the quantities Briareus computes (flux Pₗ/Yₗₘ moments, their histories, their
error sensitivities) mapped to what NIF measures: hot-spot velocity v_HS (nToF Doppler
shift), hot-spot P2/P4 (x-ray/neutron images), ρR asymmetry (activation detectors),
RKE and yield. Confidence follows research-notes: **H** = primary text read (OSTI/arXiv
accepted manuscript or publisher full text), **M** = abstract or secondary, **L** = my
inference or arithmetic. Things already in research-notes (2.7 µm per 1% cone fraction,
40–50 µm/Å, MacGowan's ~100 km/s per 1%, the Kritcher 2024 swing remarks) are not repeated
except where a primary source changes them.

**Normalization.** Every mode-1 number below is the Legendre ratio a₁/a₀ of flux on the
capsule: 2% means +2% at one pole and −2% at the other (4% peak to valley). Spears 2014
states this explicitly. Young 2020 uses c₁ₘ = 100·a₁ₘ/a₀₀, the vector norm of the three
components. In Briareus that is `mode_amplitude(c, 1)`, not c₁₀/c₀₀, which is a₁/a₀ ÷ √3
in the orthonormal Yₗₘ.

## 1. Mode 1

### 1a. Asymmetric-piston model
Hurricane et al., *Phys. Plasmas* **27**, 062704 (2020), doi:10.1063/5.0001335 (OSTI
1637584, LLNL-JRNL-799257). **H**
- The shell asymmetry fraction is f = (ρR_max − ρR_min)/(ρR_max + ρR_min). Eq. 23:
  **f = v_HS/v_imp = √(RKE/(c·KE))**. c = 1 for incompressible pistons and c ≈ 0.66 for
  NIF-like simulations. The 1D hot-spot energy is about 66% of the fuel kinetic energy.
  RKE = ½ m_shell v_imp² f². The hot spot moves toward the low-ρR piston.
- Degradation factors for low α-heating, each relative to 1D:
  - P_stag ∝ (1−f²)^{5/2}
  - E_hs, Pτ and T ∝ (1−f²)
  - no-α yield ∝ **(1−f²)^a** (Eq. 20), where ⟨σv⟩ ∝ Tᵃ and a = 3.3 for 4–7 keV
- α-heating (Eq. 21): **Y_amp,3D/Y_amp,1D ≈ exp(χ₁D^{1.2}[(1−f²)^{a−1}−1]) ≈
  exp(−(a−1)χ₁D^{1.2} f²)**. Here χ is the generalized Lawson parameter (Y_amp ≈ exp(χ^{1.2})).
- Tolerance (Eq. 26): ΔY/Y₁D = 1 − (1−f²)ᵃ e^{−(a−1)χ₁D^{1.2}f²}. Table I gives the
  allowed f and v_HS for v_imp = 380 km/s, near ignition (χ ≈ 1–2):

  | ΔY/Y | f max | v_HS max |
  |---|---|---|
  | 0.05 | 0.08 | 30 km/s |
  | 0.10 | 0.11 | 40 km/s |
  | 0.15 | 0.14 | 50 km/s |
  | 0.20 | 0.16 | 60 km/s |
- The model over-predicts v_HS for some shots near 80 km/s, which the authors attribute
  to saturation.
- Hurricane et al., *Phys. Plasmas* **29**, 012703 (2022), OSTI 1840112. **H**
  - General asymmetry: RKE/(c·KE) = 1 − ρR_WHM/ρR_ave (Eq. 14).
  - ρR_WHM = ∫dA/∫(dA/ρR) is the area-weighted harmonic mean of shell ρR over the
    hot-spot surface.
  - Table I swaps (1−f²) for ρR_WHM/ρR_ave: P_stag ∝ (·)^{5/2}, yield ∝ (·)ᵃ.
  - "Swing" (time-dependent f): the shell asymmetry *at stagnation* dominates. Swing
    matters only for low-pressure implosions (b > 0.1, P_stag/P_pv ~ 10s), not when
    P_stag/P_pv ~ 100s.
  - Applied to NIF data, the typical yield is 0.66 ± 0.14 of the symmetric yield.

### 1b. Flux P1 → v_HS (linear transfer and its time weighting)
Young et al., *Phys. Plasmas* **27**, 082702 (2020), doi:10.1063/5.0009746 (OSTI
1650438). **H** This is the closest published analogue of a Briareus `drive_history` → v_HS layer.
- **Method.** VisRad view factors are run at 100 time points with as-delivered powers
  for all 192 beams. Windows are patches with albedo 0.3 against 0.8 for the wall, and
  the gaps are ignored. **V_x = Σᵢ G(t_f,tᵢ) c₁₁(tᵢ)Δtᵢ, V_y = −Σ G c₁,₋₁ Δt,
  V_z = −Σ G c₁₀ Δt** (Eq. 9). G is in km/s per % per ns.
- **How G was built.** 59 capsule-only HYDRA runs, each with a 2.5%, 120 ps c₁₀ impulse,
  stepped every 140 ps from 0.2 to 8.32 ns. v_HS comes from synthetic neutron spectra.
- **Shape of G.**
  - The sign inverts early: foot mode-1 makes a thin spot and pushes the hot spot
    *toward* the high-drive side.
  - The rise and peak dominate, and the first half of peak is weighted more for N170601
    and N170827.
  - The static model (all radii at 85%, albedo 0.8/0.3) differs from the dynamic one by
    0.70% (abs) in mode 1 over all time but 0.14% in peak, and by <5% in v_HS
    (38.7 against 40.5 km/s).
- **Linearity.** It holds to about 0.75% constant c₁₀ in 2D HYDRA, after which v_HS
  saturates. The model was accurate to about 100 km/s on N180128-family shots.
- **3D check (N170601).** The net rise-and-peak mode 1 was 0.44%. 3D HYDRA gives
  **37.6 km/s** and the VF+G model **44.7 km/s**, with directions within 1°. That is about
  85–100 km/s per % of peak-time P1. **H** for the numbers, **L** for the ratio.
- **Source sensitivities (N170827, static model).**
  - Capsule offset in the hohlraum: 0.1 km/s/µm horizontal, 1.1 km/s/µm vertical.
  - Whole-target offset from TCC: 0.1 and 1.3 km/s/µm.
  - Hohlraum tilt: 1.6 km/s per degree.
  - Windows: **5.4 km/s per 10% wall–window albedo mismatch**, 27 km/s total for N170827.
  - A dropped inner quad gives about twice the v_HS of a dropped outer quad.
  - A dropped single beam is 4–5× smaller than a dropped quad.
  - As-delivered imbalance adds about 30 km/s.
  - v_HS scales linearly with the beam-to-beam deviation multiplier m in
    P_adj = mPᵢ + (1−m)P_avg (Eq. 10).
- **nToF.** About ±14 km/s (1σ) in magnitude; direction ±13° in θ and ±29° in φ.

Spears et al., *Phys. Plasmas* **21**, 042702 (2014), OSTI 1229835. **H**
- 2% P1/P0, applied in the 4th pulse (peak), gives a neutron-weighted v_HS of 120 km/s
  and a jet of up to 500 km/s.
- Yield falls almost 2× (4.5e15 → 2.4e15 in the text; Table I lists 4.1e15 → 1.5e15).
- Recommended **x-ray P1/P0 < 0.5%**.

MacGowan et al., *HEDP* **40**, 100944 (2021), doi:10.1016/j.hedp.2021.100944
(LLNL-JRNL-804082). **M**: full text not accessible (Elsevier 403); abstract via search.
- Static view factors with measured foot and peak powers and pointing, on 39 HDC shots.
- v_HS of 20–100 µm/ns correlates with the metric.
- Each source is about 0.5%. Summed, they give up to 100 km/s and **35% yield
  reduction**.
- Casey 2023 (**H**) describes the combination: seeds add as **vectors**, and agreement is
  judged by χ²_ν = (v_HS − v_pred)ᵀC⁻¹(v_HS − v_pred).

Rinderknecht et al., *PRL* **124**, 145002 (2020), OSTI 1616030. **H**
- v_HS clusters in azimuth −20° < φ < 160° (centered near 70°), toward the windows.
- **The ρR mode-1 amplitude is 39 ± 7% per 100 km/s.**
- In 2D simulations, up to 2% peak-to-mean flux asymmetry spans the observed ≤120 km/s.

### 1c. Capsule-thickness equivalence
Casey et al., *PRL* **126**, 025002 (2021), OSTI 1769078 (LLNL-JRNL-810637). **H**
- HDC thickness mode-1 of about 0.5% of the total thickness gives ρR asymmetry of about
  25% and v_HS of about 100 km/s. N181104 had 0.78 ± 0.2% of 78 µm.
- Simulated sensitivity is about **110 km/s per % thickness mode-1**; the data give
  **140 ± 30 km/s/%**.
- **Flux mode 1 (in % of peak flux, whole pulse) has about the same sensitivity, with
  the opposite phase.**
  - The hot spot moves *toward the thicker side* and *away from the high-flux side*.
  - Rocket estimates: Δv/v ≈ −Δm₀/m₀ and Δv/v ≈ (7/8)ΔF/F.
  - **To first order, +1% thickness at a pole ≈ −1% flux there.** **H**/**L**
- About 75% of cases are explained by radiation (MacGowan); the remaining ~25% motivated
  the capsule seed.

### 1d. High-yield amplification of v_HS (burn on/off)
Casey et al., *Phys. Rev. E* **108**, L053203 (2023), OSTI 2283276. **H**
- Spherical-convergence extension (Eq. 4): **v_HS ≈ ½(p₁/p₀) Rᵢ cosh(π/2) ·
  1.12√(P[Gbar]·R[µm]/M[µg])** (km/s).
  - Rᵢ is the initial radius in µm; the last factor is 1/τ in ns⁻¹.
  - M = c·M_shell, with c ≈ 0.67.
  - p₁/p₀ is the seed: flux or mass mode-1.
- RKE: **f² = nRKE = (p₁/p₀)² (Rᵢ²/M)(PR/vᵢ²)** (Eq. 6). The yield factor is
  (1−f²)ᵃ e^{−1.2(a−1)χ^{1.2}f²} (Eq. 5), with a = 6 − 1.5 ln T[keV].
- v_HS sensitivity is up to 2× higher at yields above 1e17.
  - For Y₁D = 4 MJ, a 0.5% seed costs about **50%** of the yield; for 0.1 MJ it costs about 14%.
  - Fig. 4a reads about 120 km/s at 0.5% for 4 MJ and about 55 km/s for 0.1 MJ
    (read off the figure: **M**).
- Seeds, combined as vectors with these sensitivities, explain 4 of 6 N210808-family
  shots (χ²≈1), against 2 of 6 with static sensitivities. N210808 itself was degraded
  about 30%.

Nora et al., *Phys. Plasmas* **30**, 092705 (2023), doi:10.1063/5.0153421 (OSTI
2439462). **H**
- v_HS = v_HS,no-a₂,no-burn × Vhs_amp^shape × Vhs_amp^burn, decoupled.
- Shape term (Eq. 2, a 2D fit): **1.31 P₂^NIS + 0.71** for P₂^NIS < 0.35, and
  **−0.65 P₂^NIS + 1.4** for P₂^NIS > 0.35. P₂^NIS is the fractional hot-spot P2/P0 from
  the neutron image. The "0.8–1.2 for |P2| < 40%" of the conclusion is Table I (3D,
  against P2 and the polar angle of the mode 1), not Eq. 2 (checked against the paper
  2026-10-24; research notes, "Young 2020 validation").
- Burn term (Eq. 4): **0.324 ln(Y_amp + 2.239) + 0.653** for Y_amp < 40, and
  **13.105 ln(Y_amp + 1585.225) − 95.539** for 80 < Y_amp < 500 (Fig. 6's coefficients;
  Eq. 4 prints them rounded to 13, 1585, −95, which shifts the value by −0.24). Table IV
  lists the inverses per shot.
- Prolate or oblate stagnation deflects the flow (their Fig. 4, Tables I and II).
- In the N210808 family, measured v_HS was inflated up to 184%. Read naively it implies
  a₁/a₀ ≈ 1.8% with YOC 81%; unamplified it is about 0.78% with YOC 94%.
- Database mean total amplification is 1.07 ± 0.42.

Kritcher 2024 (paper.txt), Table II. **H** The burn-on/burn-off ratios are:

| Shot | Ratio |
|---|---|
| N220919 | 51/40 = 1.28 |
| N221204 | 36/27 = 1.33 |
| N230618 | 107/84 = 1.27 |
| N230729 | 104/75 = 1.39 |
| N230904 | 82/63 = 1.30 |

My arithmetic (**L**). At Y_amp ≈ 30 Nora's Eq. 4 gives ≈ 1.78, so the LLNL burn-on/off
simulations and Nora's fit disagree by about 30% (**L**; Y_amp for these shots not checked).
- Kritcher's "~0.7% applied mode 1 → >100 km/s burn-on, ~70–75 km/s burn-off" means
  **~105 km/s per % (burn-off)** and about 145 km/s/% burn-on. That is consistent with
  Young and MacGowan (**L**).

### 1e. Mode 1 in an empirical yield model
Ralph et al. (HYBRID-E/I-Raum burning plasma), *Nat. Commun.* **15**, 2975 (2024),
doi:10.1038/s41467-024-47302-8. **H**
- Θ_m1 = (1−f²)^{3.3}·Y_amp,3D/Y_amp,1D with f = v_m1/v_imp (Eq. 3; the extracted text is
  garbled; it cites the piston model).
- Θ_m1 = 0.99–0.65 for HYBRID-E (11–103 km/s) and 0.99–0.46 for I-Raum (17–123 km/s).
- N210220 deliberately retuned 3D laser power to cancel the known mode-1 seeds.

## 2. Modes 2 and 4: flux history → shape

Masse et al., *Phys. Plasmas* **26** (6) (2019; article number not checked), doi:10.1063/1.5092827 (OSTI
1885111). **H**
- Linear Green's functions G(t,a) map flux P2(a) to ablation-front, hot-spot, ρR,
  radiograph and self-emission P2(t): ΔR(t) = Σᵢ f(aᵢ)G(t,aᵢ) (Eq. 5).
- Built from 10%, 200-ps impulses every 175 ps in the 2-shock CH platform.
- **Maximum sensitivity is near the rise to peak, around 6 ns of 8 ns.**
- **The hot-spot P2 response at bang time changes sign before about 4.5 ns.** Early
  waist-hot drive gives an *oblate* hot spot, while the same perturbation at 6 ns gives a
  prolate one.
- **The ablation front and ρR keep one sign.** That is why hot spot and shell are not
  conformal, and why tuning the self-emission round "doesn't usually correct the swing
  nor the fuel ρR symmetry".
- Inverse problem: χ² against radiography and self-emission data, by simulated annealing
  or GA, with |P2| ≤ 0.25 and slope limits. It is ill-conditioned (they took 2000
  solutions).
- Proposed correction: early pole-hot P2 to cancel late pole-hot, giving a "round and
  swing-less" shape.
- Magnitudes appear only in figures. **H** for the method and sign structure.

Kritcher et al., *Phys. Plasmas* **21**, 042708 (2014), doi:10.1063/1.4871718. **M**
(abstract only).
- More than 600 2D runs with P2/P4.
- The metric combines hot-spot shape, ρR shape and RKE. Counter-asymmetries that round
  the core do not recover ρR or RKE.
- Degradation correlates with time-dependent shape: shock, in-flight and core.

Scott et al., *PRL* **110**, 075001 (2013), arXiv:1212.3123. **H**
- P4 flux applied 100 µm from the ablation front.
- **The sign of the stagnated P4 reverses with timing.** High flux during shock
  compression protrudes outward; during peak drive it protrudes inward.
- **10% flux P4 → hot-spot a₄ ≈ 20 µm and up to 15× yield loss.** The loss comes
  through RKE and lateral mass flow.
- Positive P4 aliases as a false oblate P2 in GXD images, and "correcting" it worsens
  the implosion.
- Expected capsule flux asymmetry was <3% except in the first 2 ns (up to 10%).

Milovich et al., *Phys. Plasmas* **23**, 032701 (2016), OSTI 1763929. **H**
- The picket (first 2 ns) and the peak set the core shape; the 2nd and 3rd shocks
  barely matter (±10% cone fraction).
- Early-time P2 must stay **< 7.5%** (P2/P0), keeping margin loss < 10% and hot-spot
  deviation < 10% rms.
- That corresponds to a cone-fraction window of ±0.035.
- View factors: **+6% P2 per 200 µm outward outer-cone pointing**, and **2% P2/P0 per
  10% albedo change**.

Kritcher 2024 (paper.txt). **H**
- Thicker capsules sample more of the late pole-hot drive and become more oblate;
  larger Δλ compensates with earlier waist-hot drive.
- The P2 swing "is currently unknown". It is calibrated with in-flight radiography at
  about 200 µm and was missing for N221204.
- −10.5 µm (oblate) to +6 µm (prolate) mattered.
- **Reducing hot-spot P2 by 5 µm → ~90% more yield (N230729).**
- Yield vs hot-spot P2 spreads with the assumed swing (Fig. 4 bottom).

Ralph et al. 2024 (above). **H**
- **Θ_P2 = A·exp(−[(P2/P0 − μ)/W]²) + B**. The form is my reading of the garbled Eq. 5;
  its parameters are given in Methods.
  - HYBRID-E: A = 0.9225, μ = 0.1180, W = 0.3058, B = 0.0775.
  - I-Raum: A = 0.9502, W = 0.1820, B = 0.0498.
- Hot-spot x-ray P2/P0 is a fraction. **50% yield at ±25% (HYBRID-E) and ±17% (I-Raum)
  from the optimum.**
- Simulated sensitivity is about 0.7 without α and 0.1–0.2 with α, for ±50% P2/P0.
- The simulated optimum sits at −7.3% (HYBRID-E) against +11% in experiment, an
  unexplained offset.

Lindl et al., *Phys. Plasmas* **33**, 050503 (2026), doi:10.1063/5.0316678. **M** (AIP
full text read through a summarizing fetch).
- m = 4 is "typically no larger than 3 µm" on cores and falls off as sin⁴θ toward the
  equator.
- About 1.5× yield came from reducing m = 4 between a pair of shots.
- The quantitative effects of P4 and m = 4 "remain inconclusive".
- N231119's large mode 1 lowered the average DSR by about 5%.

## 3. RKE → yield
- Piston model (**H**): f² = RKE/(c·KE), c ≈ 0.66. Y/Y₁D = (1−f²)ᵃ, times the α factor
  above. Mode-2/4 RKE enters through ρR_WHM (Hurricane 2022). Hurricane 2020 Fig. 3 shows
  2D runs with modes 1, 2 and 4 falling on the same (1−f²)^{3.3} curve when f² = RKE/E_sym.
- Kritcher 2024 Fig. 4 (top) gives yield against RKE (burn-off sims). Thicker,
  high-dopant capsules have margin at the RKE where N230729 sits; the "cliff" is
  further out. No closed form is published. **H**
- Casey PRL 2021 (**H**): ~25% ρR asymmetry (v_HS ≈ 100 km/s at 400 km/s) costs ~7% of
  E_hs and ~22% of no-α yield ((1−f²)^{10/3}). With α, at Y₀ = 2e16 and DSR 3.3%, it
  costs ~38%.

## 4. Error budgets
- MacGowan 2021: each source is about 0.5%; combined, up to 100 km/s and −35% yield. **M**
- Young 2020: N170827 windows 27 km/s, plus about 30 km/s of delivery imbalance on top
  of drop scenarios. **H**
- Casey PRL 2021: 75% of cases explained by radiation. Casey PRE 2023: seeds (laser
  power-balance view factor + measured capsule mass) explain about 67% (4/6) of
  N210808-family shots; about 25% historically stay anomalous. **H**
- Kritcher 2024: N230729's 104 km/s was "largely" laser delivery imbalance, aligned with
  capsule, ice and window seeds. N230618's came from capsule thickness plus delivery.
  **No per-cone/quad kJ or % is published there. H**

## 5. Transfer-coefficient table

| Drive quantity (Briareus) | Observable | Coefficient | Validity | Time weighting | Source |
|---|---|---|---|---|---|
| a₁/a₀ flux, whole pulse | v_HS (burn-off) | ~100 km/s per % (Young 3D: 37.6–44.7 km/s at 0.44%) | linear ≲0.75%, saturates ≳100 km/s | G(t): foot negative, rise and peak dominant | Young 2020 H; MacGowan 2021 M |
| a₁/a₀ peak | v_HS | 120 km/s at 2% (sat.) | CH NIC, 2D | 4th pulse | Spears 2014 H |
| a₁/a₀ | v_HS (burn-on), high yield | ½(p₁/p₀)Rᵢcosh(π/2)·1.12√(PR/M) | 1e15–1e18 | seed at peak velocity | Casey 2023 H |
| HDC thickness mode-1 (%) | v_HS | 110 (sim), 140±30 (data) km/s/%; sign opposite to flux | ≲100 km/s | — | Casey 2021 H |
| Window albedo mismatch | v_HS | 5.4 km/s per 10% | N170827 | static | Young 2020 H |
| Capsule/target offset | v_HS | 0.1 (horizontal), 1.1–1.3 (vertical) km/s/µm | <200 µm | static | Young 2020 H |
| v_HS | ρR mode-1 | 39±7% per 100 km/s | HDC/BF/CH | — | Rinderknecht 2020 H |
| f = v_HS/v_imp | no-α yield | (1−f²)^a, a≈3.3 | Y_amp<2 | stagnation shell | Hurricane 2020 H |
| f, χ | α yield factor | exp(−(a−1)χ^{1.2}f²) [Hurricane] or exp(−1.2(a−1)χ^{1.2}f²) [Casey] | moderate burn | — | H |
| Y_amp | v_HS burn amplification | 0.324 ln(Y_amp+2.239)+0.653 | Y_amp<40 | — | Nora 2023 H |
| hot-spot P2/P0 | v_HS shape multiplier | 1.31P₂+0.71 / −0.65P₂+1.4 | ±40% | — | Nora 2023 H |
| flux P2(t) | hot-spot, ρR, radiograph P2 | Green's function G(t,a) (figures only) | linear to stagnation | peak-dominated; hot spot sign flips before ~4.5 ns | Masse 2019 H |
| early flux P2 | margin | < 7.5% → <10% margin loss | picket | first 2 ns | Milovich 2016 H |
| flux P4 | hot-spot a₄, yield | 10% → ~20 µm, ≤15× loss | 2D CH | sign flips shock vs peak | Scott 2013 H |
| hot-spot P2/P0 | yield factor | Gaussian, W=0.306, μ=+11.8% | HYBRID-E | — | Ralph 2024 H |
| hot-spot P2 (µm) | yield | −5 µm → ×1.9 | N230729 | — | Kritcher 2024 H |

## 6. What Briareus should implement (minimal)

1. **`implosion_metrics(c; vimp, ...)`**, a small unexported-first layer (`implosion.jl`
   after `uncertainty.jl`). Inputs are mode coefficients or a mode history.
   - Mode-1 vector: `p1 = √3·(c₁ₘ/c₀₀)` mapped to Cartesian. Its norm is
     `mode_amplitude(c,1)`. The pushing direction is opposite the high-flux pole.
   - `v_HS = k_v · p1` (km/s), with `k_v` a documented engineering default of 100 km/s
     per % (Young/MacGowan; Kritcher's 0.7% → 75 km/s burn-off is the check). Say that
     it is platform-specific and linear below about 0.75%.
   - `f = v_HS/vimp`. No-α yield factor `(1−f²)^a` (a = 3.3, or 6 − 1.5 ln T).
   - Optional α factor `exp(−(a−1)χ^{1.2}f²)`, with a keyword for Casey's 1.2.
   - RKE fraction `f²·c`, c = 0.66.
   - Optional burn amplification (Nora Eq. 4) and shape multiplier (Nora Eq. 2) for
     comparison with measured v_HS.
   - Validate on:
     - Young's N170601 test (0.44% → 37.6/44.7 km/s);
     - Kritcher N230729 (0.7% → 75 km/s burn-off);
     - Hurricane Table I (f ↔ ΔY), as an exact closed-form unit test;
     - Spears (2% → 120 km/s, saturated: an expected over-prediction).
2. **A time-weighted mode-1 metric over `drive_history`.**
   - `v_HS = Σᵢ G₁(tᵢ)·p1(tᵢ)Δtᵢ` with a user-supplied `G₁` (km/s/%/ns) table.
   - Ship a placeholder: zero in the foot, flat over rise and peak, normalized so a
     constant 1% gives `k_v`. Label it as derived.
   - The real G is platform-specific; Young Fig. 5 holds three tables, read off the
     figure only.
   - Validate on Young's static and dynamic agreement (<5%) and on the quad/beam-drop
     ratio (inner about 2× outer).
3. **P2/P4 swing metrics.** These are generic Green's-function convolutions, the same code
   as item 2.
   - Report the time-integrated P2 over a user window (default: rise and peak).
   - Report the swing: max − min, and the integral of |P2(t) − its mean|.
   - Report the picket-only P2 against Milovich's 7.5%.
   - Accept response kernels for hot spot and ρR. The hot-spot kernel changes sign before
     the rise, while ρR keeps one sign, so the two can disagree.
   - With no kernel, warn that hot-spot P2 cannot be inferred; never output µm by default.
   - Validate qualitatively on Masse's sign structure. Quantitatively it needs kernels.
4. **Error budgets in `mode_uncertainty`.**
   - Propagate the per-source mode-1 covariance through `k_v` (or G₁) into a v_HS
     covariance (3×3).
   - Report the per-source v_HS rms and the expected yield factor
     E[(1−f²)ᵃ] ≈ 1 − a·E[f²], where E[f²] = (|v̄|² + tr Σ)/v_imp² (**L**).
   - Add window albedo and target offset as sources. Young's 5.4 km/s per 10% and
     ~1.2 km/s/µm vertical are checks on the existing derivatives.
5. **`infer_errors` from v_HS.**
   - Accept a measured v_HS vector with nToF covariance (±14 km/s; 13° in θ, 29° in φ)
     as the measurement. Map it to p1 by dividing by `k_v × amp_shape × amp_burn`.
   - Add a capsule-thickness mode-1 nuisance term with *opposite* sign (Casey 2021),
     with metrology covariance of 0.1–0.15 µm per axis.
   - Report MacGowan's χ²_ν.
   - Validate on Casey PRE Fig. 3 (N210808 family). Shots and seeds are listed but the
     per-shot seed values are not; this needs LLNL data.
6. **Objectives.** Add an optional `VelocityPenalty` (k_v·p1, or the G₁-weighted p1) so
   mode 1 is weighted by its time history, not equally at all times. It needs no new
   physics.

## 7. Open questions
- MacGowan 2021 full text (Elsevier) could not be read. Its exact metric (time window,
  foot/peak weights, how window losses are modeled, fit coefficients) is still **M**.
- The Green's-function values are only in figures (Young Fig. 4–5, Masse Fig. 5). No
  published HYBRID-E/N221204 G₁ or P2 kernel was found.
- The 1.2 factor in the α-amplification exponent differs between Hurricane 2020 Eq. 21
  and Casey 2023 Eq. 5.
- Nora's burn amplification (about 1.8 at Y_amp ≈ 30) disagrees with the ~1.3 burn-on/off
  ratios in Kritcher Table II. They may differ in what is held fixed (P2, Y_amp
  definition). Y_amp for the 2022–23 shots was not checked.
- No published per-cone or per-quad imbalance (kJ or %) for N230729.
- No closed-form RKE → yield for the N221204 design (only Kritcher Fig. 4). The P2-swing
  magnitude is "unknown" (Kritcher 2024).
- No transfer for hot-spot P2 in µm per % flux P2 for HYBRID-E beyond 2.7 µm per % cone
  fraction and 40–50 µm/Å (research notes).
- The Kritcher 2014 metrics paper (time-integrated P2 requirement, swing tolerance) was
  only seen as an abstract. The "~1% time-integrated P2" ignition spec was not found in a
  primary source this round.
- For P4/m = 4, only Scott 2013 (2D, CH) and the Lindl 2026 review give numbers. A
  Michel DPP 2024 m = 4 result was not found.
- The direction convention when mapping Briareus's real Yₗₘ (c₁,₋₁, c₁₀, c₁₁) to
  Young's Cartesian (V_x = G c₁₁, V_y = −G c₁,₋₁, V_z = −G c₁₀) must be checked against
  Young's associated-Legendre convention before comparing vectors.
