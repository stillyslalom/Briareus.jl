# Hohlraum dynamics and the P2 swing: what a time-dependent view factor needs

Literature review, 2026-10-22, for ROADMAP "Next 9" (time-dependent view factor) and
"Next 10" (shot calibration). Confidence follows `reference/research-notes.md`:
**H** = primary text read (the full LLNL preprint or the journal text), **M** = abstract or
secondary source only, **L** = our inference. Numbers marked *derived* are our own
computations. Full texts were read from the OSTI/arXiv copies cited below (not kept in the repository).

A distinction runs through everything below. Two different surfaces move:

- **the Au/He interface ("bubble tip")**, the low-density high-Z plasma that absorbs inner
  beams. It moves at about **300–460 µm/ns**.
- **the x-ray conversion surface** (about 1 g/cm³, the base of the bubble), which radiates.
  It moves at about **50–100 µm/ns**.

A view-factor code should move its radiating wall with the second and put an inner-beam
absorber at the first.

## 1. Bubble expansion: measured speeds and scaling laws

**Callahan et al., Phys. Plasmas 25, 056305 (2018), doi:10.1063/1.5020057** (preprint
LLNL-CONF-742209, OSTI 1788342, full text) — **H**
- Hypothesis: in low-fill (0.3–0.6 mg/cc) hohlraums, P2 is set by the Au/DU bubble size
  relative to R_hohl and by r_cap/R_hohl. The bubble is launched by the outer-cone picket
  and "grows approximately linearly throughout the duration of the laser pulse".
- Eq. 1 (vacuum): R_bub/R_hohl ∝ c_s τ/R_hohl ∝ √(E_picket,outer/(A_outer X_d ρ_wall)) · τ/R_hohl.
- Gas tamping, taken from simulations: the bubble speed scales as √(ρ_wall/ρ_fill).
- Eq. 2: R_bub/R_hohl ∝ √(E_pic,out/(A_out ρ_fill)) · τ/R_hohl.
- Eq. 3: P2 ∝ −(R_bub/R_hohl)(r_cap/R_hohl).
- Eq. 4: P2 = 0 at a 33% cone fraction (CF) when the metric √(E/(Aρ_fill))·(τ/R_hohl)(r_cap/R_hohl) ≈ 1.2.
- Units and inputs:
  - E in kJ, A in mm², using the one-side outer spot area of 53 mm² (Table I).
  - ρ_fill in mg/cc, τ = total pulse length in ns, R in mm.
- Data and fit:
  - 39 shots: CH/HDC/Be ablators, τ 5.5–14.1 ns, CCR 3.0–4.2, hohlraum diameter
    5.4–6.72 mm, Au, DU or Au-lined DU walls.
  - R² = 84%, standard deviation **8 µm** (the goal was ±5 µm).
  - Alternative fluence exponents give R² 76% (1/3) and 79% (2/3). The (r_cap/R)^1.5 form fits about as well.
- Cone-fraction sensitivity: 2.7 µm of hot-spot P2 per 1% CF (range 2.0–3.7, seven campaigns).
- Worked design (Table I): 6.72 mm hohlraum, r_cap 1076 µm, 60 kJ picket, 7.7 ns,
  0.45 mg/cc. Predicted P2 ≈ +2.4 µm.
  - *Derived* check: the metric is 1.586 × 2.292 × 0.320 = 1.16. With Izumi's slope
    (below) that gives +3 µm, so the units above are consistent.

**Ralph et al., Phys. Plasmas 25, 082701 (2018), doi:10.1063/1.5023008** (LLNL-JRNL-738440, OSTI 1469454) — **H**
- Platform: CH high-foot, 6.72 mm Au-lined/DU hohlraum, 0.6 mg/cc. Shots are 2DConA,
  keyhole and DT.
- Bubble trajectory (N160504, polar framing camera):
  - Constant inward speed of **182 ± 11 µm/ns**, extrapolating to the initial wall
    (R = 3357 µm) at 2.47 ± 0.24 ns.
  - That shot's picket had about 3× the intended energy. Nominal pickets give bubbles that
    stay outside the LEH view.
  - 30° beams meet the bubble when it reaches r = 2.4 mm.
- Speed model: U_bub = U_wall + U_las, with U_las = c_s ln(Q_las) (isothermal expansion)
  and U_wall = c_s ln(Q_wall), Q_wall ≈ 10.
  - The two terms are about equal at 0.6 mg/cc (T_r 100 eV, T_e 400 eV).
  - Result: **U_bub ∝ I^0.36**. Fits hold over I^1/3 to I^2/3.
- Absorption metric: α = β·τ·(E_out,pic/a)^0.36. Here a is the outer spot area,
  E_out,pic is the outer-cone picket energy over its FWHM net of backscatter, and τ runs
  from the pulse start to mid-peak.
- Effective cone fraction: CF_eff = CF(1 − α). This fixes **inner-beam absorption at
  51 ± 2% to 62 ± 2%**. The anchor is a view-factor calculation with the *measured wall
  position at peak*, which needs **CF = 18.5%** for zero flux P2.
- Shape: 2DConA in-flight P2 (at about 200 µm) and hot-spot P2 become more oblate together.
  - At 0.3 mg/cc, N160202: hot spot −29 ± 3 µm, in flight −35 ± 3 µm.
  - High picket, N151220: hot spot −22 ± 3 µm.
- Inverse-bremsstrahlung absorption in the bubble against the He: κ_wall/κ_He ≈ Z_wall/Z_He ≈ 20.

**Izumi et al., Rev. Sci. Instrum. 89, 10K111 (2018), doi:10.1063/1.5039364** (OSTI 1835034/1545514) — **H**
- Platform: HDC SymCap in a 5.75 mm Au hohlraum, 10.13 mm long, 3.37 mm LEH, 0.3 mg/cc
  basecamp, with 8 µm "thin-wall" windows at the waist.
- Au/He boundary speeds: **350 µm/ns at 0.3 mg/cc, 400 µm/ns at 0.15 mg/cc**. At low fill
  the boundary runs 240 µm ahead.
- Thin-wall inner-beam emission:
  - Brightest at 5.36 ns, then dims.
  - Upper and lower inner spots separate at 5.8 ns, when the ablator collides with the wall
    plasma (5.45 ns with 1.5× picket).
- Core P2 over the six single-parameter shots ranges from **−18 to +21 µm**, with +5.6 µm
  at basecamp. It tracks the time-integrated thin-wall signal.
- Inner-beam attenuation depends on the line-integrated density along the beam, not on
  bubble displacement alone.
  - 0.15 mg/cc moves the bubble further yet gives a 2× brighter thin-wall signal.
  - A 1.2× larger capsule collides about 400 ps earlier.

**Izumi et al., Phys. Plasmas 28, 022706 (2021), doi:10.1063/5.0030826** (LLNL-JRNL-812268, OSTI 1773252) — **H**
(the most complete parametric model found)
- Restates Callahan: Core P2 − 270(CF − 0.333) = 85(1.2 − FoM) µm, equivalently
  P2 = 90(3CF − FoM).
- Bubble tip x_bub (displacement from the wall):
  - ballistic, launched at the outer turn-on (+0.3 ns); fits intercept at 0.5 ± 0.1 ns;
  - baseline about **345 µm/ns**;
  - **x_bub ∝ ρ₀^−0.2±0.03**: 0.15 vs 0.3 mg/cc gives +250 ± 30 µm (14 ± 2%) at x̄ = 1.8 mm;
  - **x_bub ∝ E^0.2±0.03**: 1.5× picket gives +160 µm (+9%); a 27% fluence drop from
    cone/quad split gives −120 µm;
  - ablation scaling cited: dm/dt ∝ ((1 − η_x)E)^1/3.
- Inner-beam geometry at z_bub = ±3 mm (R_hohl = 2.88 mm): the inner cones form an annulus
  **1.3–2.4 mm from the wall** (R_inner = 1.1 mm).
- Clipping fraction (Eq. 3, "opaque iris"): f = 0 until x_bub > R_hohl − R_o, rising to 1 at
  x_bub = R_hohl − R_i.
  - Our reading of the garbled text: f = [R_o² − (R_hohl − x_bub)²]/(R_o² − R_i²) (**L** on
    the exact form).
  - Power-law fit (x_bub/(R_hohl − R_i))^2.5, within 15% for 0.3 < f < 1.
- Absorbed fraction (Eq. 6): CF_bub/CF = f[1 − exp(−G Z_Au ρ_bub² Δz (x_bub/1 mm)^0.5/T_e^1.5)].
  - Chord length 1.1 → 1.5 mm.
  - Bubble density ρ = S ρ₀ × (R_hohl/R_bub)^k, with S ≈ 3 from Rankine–Hugoniot at
    v_bub/c_s ≈ 345/140. The exponent k is illegible in the extracted text (**L**).
  - Computed bubble transparency at 5.5 ns: **29% (0.3 mg/cc) and 54% (0.15 mg/cc)**.
- Final FoM (HDC, 0.15–0.3 mg/cc): about (E^0.2 t′/R_hohl)^2.5 ρ₀^(∓0.13) (R_cap/R_hohl)^1.1,
  ordinate P2/M^1.5. Slope −115 µm per unit FoM, R² = 0.96 over 7 shots.
  - t′ is the time from the outer picket to mid-peak, 5 ns here.
  - The exponents are legible; the ρ₀ sign is ambiguous in the extraction (**M** on the form).

**Meezan et al., Phys. Plasmas 27 (10) (2020), doi:10.1063/5.0018733** (OSTI 1738889) — **H**
- Measured bubble speeds (µm/ns = km/s), against LASNEX with flux limiter f = 0.03 and f = 0.15:

  | Shot | Measured | f = 0.03 | f = 0.15 |
  |---|---|---|---|
  | N171010 baseline | 362 ± 40 | 334 | 215 |
  | N171212 1.5× picket | 356 | 369 | 356 |
  | N180102 1.2× r_cap | 367 | 308 | 211 |
  | N180103 0.15 mg/cc | 461 | 462 | 374 |
  | N180117 cone/quad split | 321 | 261 | 184 |

- Restricted conduction (f = 0.03 in Au and He) is needed. Bubble T_e is hotter than
  f = 0.15 predicts.

**Izumi et al., Phys. Plasmas 30 (3) (2023), doi:10.1063/5.0136145** (LLNL-JRNL-842205, OSTI 1971298) — **H**
- Bubble speed in a 3-shock HDC design: v_bub = 375 µm/ns (N171010). Scaling cited:
  v ∝ (E/ρ_He)^0.2.
- A low-trough 4-shock CH design is *not* ballistic.
  - The picket-launched bubble stagnates in the trough and relaunches at about 6 ns, at
    most 300 µm/ns.
  - It slows on back-pressure near the axis.
  - Constant speed holds only for pulses with < 15× foot-to-peak contrast.
- Ablator–wall stagnation ridge: r ≈ 2.1 mm in a 2.88 mm hohlraum. The transmitted inner
  power stops there late in the pulse.
- Outer-spot motion: the deposition region "move[s] back along [its] incident path due to
  the x-ray induced wall ablation". The angular motion scales as (ablated mass ∝ T_r^1.8)/R
  ∝ t/R. The exact expression is lost in extraction (**M**).
- Mitigation, temporal and azimuthal beam phasing: bubble spikes +300 µm, troughs −500 µm,
  net +100 µm for inner transport.

**Robey et al., I-Raum, Phys. Plasmas 25 (2018), doi:10.1063/1.5010922** — **M**
- A recessed pocket at the outer spots buys about 1 ns of inner-beam propagation.

## 2. Effective wall radius and LEH closure

**Chen et al., Phys. Plasmas 27 (2) (2020), doi:10.1063/1.5128501** (G-LEH, OSTI 1632379) — **H**
- X-ray-emitting wall radius:
  - **600 µm** inward on N170706 (0.6 mg/cc), "a speed of ~50 km/s".
  - **500 µm (0.3) vs 400 µm (0.45 mg/cc)** through peak at 5.5 ns (5.75 mm).
- VisRad models of LEH images needed the hohlraum radius **reduced by 300–700 µm**.
- Effective LEH at 0.15 mg/cc (simulation): about 1550 µm (open) until 4.5 ns, then a
  rapid drop at 4.5–5 ns as a high-density ridge forms where the outer beams cross the lip
  blow-off.
  - Nearly constant through peak; laser heating "keeps the LEH open".
  - Faster closure after laser off.
  - Data agree with LASNEX at 0.15–0.3 mg/cc; at 0.6 mg/cc the measured radius is larger
    (within about 50 µm of f = 0.15).

**Schneider et al., Phys. Plasmas 22, 122705 (2015), doi:10.1063/1.4937369** — **H** (journal full text via fetch)
- Platform: 5.75 mm × 9.425 mm hohlraum, 3.1 mm LEH (r₀ = 1.55 mm), 0.96/1.6 mg/cc.
- Radius vs time:
  - The LEH closes in the low-power part and reopens at peak.
  - Low Foot: about 10% smaller at ¼ peak power than at peak.
  - At peak it varies ≤ 4%.
- At 20 ns, 500 eV: **1.30 mm measured vs 1.13 mm modeled**. The measured area is
  1.25–1.4× the model's.
- *Derived*: (1.30/1.55)² = 0.70 of the initial area.

**Swadling et al., "Build-A-Hohlraum", Phys. Plasmas 32, 052707 (2025), arXiv:2501.10350** — **H**
- LEH closure is faster than simulated: simulations overpredict LEH size early, match at
  peak and underpredict after. The difference is about 100 µm in diameter at peak intensity.
- Au-bubble T_e is about 31% above simulation; 2–4 keV M-band is 30% below.

**Young et al., Phys. Plasmas 27 (8) (2020), doi:10.1063/5.0009746 (article number not verified)** (LLNL-JRNL-807170, OSTI 1650438) — **H**
- This is the published form of the LLNL "dynamic view factor": VisRad plus a mode-1 Green's
  function.
- Dynamic model:
  - 100 time slices.
  - The cylinder shrinks radially; its radius is the **average of the laser-absorption and
    x-ray-emission surfaces** from 2D HYDRA.
  - The capsule ablates uniformly.
  - Albedos come from a quasi-1D wall calculation.
  - LEH closure comes from SXI images.
- Static model: **hohlraum, LEH and capsule radii at 85% of initial**, albedo 0.8 (Au)
  and 0.3 (ablator), 3ω→x-ray conversion 0.9.
- Static against dynamic: mode 1 differs by 0.70% (mean absolute, all times) but **0.14% at
  peak**. The hot-spot velocity error is < 5% (38.7 vs 40.5 km/s).
- Windows are modeled as α = 0.3 patches.
- Bubble growth "tends to effectively redistribute power from the inner to outer cones".

**MacGowan et al., HEDP 40, 100944 (2021), doi:10.1016/j.hedp.2021.100944** — **M**
- Abstract only; the full text was 403. A *static* view factor using delivered foot and peak
  powers, 39 HDC shots.
- Same group and tool as Young 2020; that it uses the 85% static geometry is our inference (**L**).

**Dewald et al., PRL 95, 215004 (2005)** (filling limit; LLNL copy) — **H**
- Mean ablated-wall ion density n_i ≈ 0.0007 T_R^1.86 t^0.54/R, with n_i in n_c, T_R in
  heV, t in ns, R in cm.
- The laser deposition migrates to the LEH when the inverse-bremsstrahlung absorption length
  ℓ_IB < 0.7 r_LEH.

**Landen et al., Phys. Plasmas 6, 2137 (1999), doi:10.1063/1.873465** — **M**
- Analytic *time-dependent* P2 view-factor model with wall motion, time-varying albedo and
  hohlraum/capsule ratio, validated on Nova and Omega. Paywalled; equations not read.
- Spot-motion prior art: Suter et al., PRL 73, 2328 (1994), Nova symmetry database (cited by
  Izumi 2023, not read).

## 3. Magnitude and measurement of the swing

- **Kritcher et al., Phys. Plasmas 28, 072706 (2021)**, HYBRID-E (OSTI 1818404) — **H**
  - Foot P2/P0 is held within ±1%. At constant 33% CF in the peak, P2/P0 swings from
    negative (inners reach the waist) to positive (they no longer propagate).
  - CBET shifts the swing negative "and does not significantly change the slope".
  - Swing specification: ±4 µm in-flight→hot spot. Matching hot-spot P2 needed a **peak-CF
    multiplier of 1.065**.
  - "Hohlraum modeling is not currently predictive for the late-time P2", attributed to LEH
    closure, bubble plasma conditions or CBET.
  - In-flight P2 at P0 = 200 µm was within error bars; the numerical value is garbled in
    the extraction.
- **Kritcher et al., Phys. Plasmas 25, 056309 (2018)** (OSTI 1463835) — **H**
  - The in-flight-to-hot-spot swing in scaled HDC designs comes from a late +P2/P0 "during
    the last several hundred pico-seconds … due to bubble expansion and capsule ablation".
  - A peak drive multiplier of 0.86–0.92 for low fill (0.78 for high fill) correlates with
    LEH size, possibly from LEH closure.
- **Kritcher 2024** (motivating paper) — **H**
  - Swing is calibrated by in-flight radiography at about 200 µm. For N230729 the swing was
    reduced ad hoc to match yield.
- **Milovich et al.** (OSTI 1763929) — **M**
  - Picket-time P2 must stay < 7.5% to keep the margin loss < 10%; measured with re-emit spheres.
- **Kritcher et al., Phys. Plasmas 21, 042708 (2014)** — **M**
  - Time-dependent shape metrics; the abstract only.

## 4. Gaps in the prior art (relevant to Briareus)

- **IRAD3D** (arXiv 1606.00934, 1703.01579): the abstracts read show no spot motion or wall
  motion. The research-notes entry claiming "spot motion" should be downgraded to **L**
  unless a full text says otherwise.
- **OSTI 789033**, "Time-dependent view factor simulations of radiation symmetry in ICF
  hohlraums" (LANL subcontract LA-SUB-99-65, 1999), has no abstract and no full text.
- **Hilsabeck et al., arXiv 2608.11664 (2026)** — **M**
  - Cu blow-off into 0.3/0.6 mg/cc He on OMEGA EP.
  - HYDRA/Gorgon expand 20–50% faster than measured at 1–3 ns. (The NIF data above find
    LASNEX f = 0.15 *slower* than measured.)

## What Briareus should implement

A minimal parametric model with four pieces, all driven by `drive_history` time.

**1. A radiating wall that shrinks: `WallMotion`.**
- Equation: R_w(t) = R₀ − v_w·max(0, t − t₀), cylinder only, ends fixed.
- Defaults (engineering, from Chen 2020): v_w = 0.10 mm/ns at 0.3 mg/cc and 0.07 at 0.45
  (*derived* from 500/400 µm over about 5 ns); 0.05 mm/ns at 0.6 mg/cc. t₀ = start of the
  outer pulse.
- Cross-check: the time-average should land near Young's 0.85 R₀.
- Implementation: rebuild the target at K keyframes (K ≈ 6–10 through the peak) and
  interpolate the moment and Y_lm responses linearly in t.
  - A uniform scaling of all lengths leaves view factors unchanged. Only the *relative*
    shrink of R_w against r_cap, r_LEH and L matters, so K can be small (**L**).
- What it gives automatically: Briareus keeps k = −port and the aim fixed, so the outer-cone
  footprint on the moved cylinder slides **toward the LEH** by ΔR/tanθ (*derived*): 0.84ΔR
  at 50°, 1.02ΔR at 44.5°, i.e. 0.35–0.6 mm for ΔR = 0.4–0.6 mm. This is the "spot motion"
  P2 swing, with no new beam physics.
- Constraint: `DepositionOperator` must be rebuilt per keyframe (the target changes, not
  the ports).

**2. LEH radius(t) via the gap machinery.**
- Model the closing annulus r_eff(t) < r < r₀ as a component with a time-dependent albedo,
  reusing `GapClosure`-style `_albedo_split` with transmitted = open hole, rather than
  remeshing.
- Default: r_eff = r₀ until the rise to peak, ramp over about 0.5 ns to **0.84 r₀** (the
  Schneider peak value, matching Young's 0.85), hold through peak.
  - Label it an engineering default; Build-A-Hohlraum suggests a further ~50 µm of radius
    closing early.
- Laser clipping on the closed annulus is out of scope: the laser keeps it open (Chen).

**3. Bubble absorption of the inner cones: `BubbleAbsorption` (Izumi 2021, Eqs. 3 and 6).**
- Tip: x_bub(t) = v_b (t − t_L), with t_L = outer turn-on + 0.2 ns.
- Speed: v_b = v_ref (E/E_ref)^0.2 (ρ₀/ρ_ref)^−0.2, with v_ref = 0.35 mm/ns at
  E_ref = baseline HDC outer picket fluence, ρ_ref = 0.3 mg/cc.
- Placement: the bubble sits at the axial position of each outer (50°) footprint.
- Per inner beam:
  - compute the clipped fraction f(x_bub) geometrically from Briareus's own inner-beam
    envelope at z_bub, which is better than the annulus approximation;
  - absorbed = f·(1 − T_b(t)), with T_b a single calibration curve:
    exp(−τ_b (x_bub/x_ref)^p (R_hohl/(R_hohl − x_bub))^2), p ≈ 0.5, τ_b fitted;
  - anchors: transmission 29% at 5.5 ns (0.3 mg/cc) and 54% (0.15).
- Redeposit the absorbed power on the wall elements under the outer footprint at z_bub,
  so it becomes x-rays near the LEH. Ralph treats it as conducted to the wall, not lost.
- This is differentiable and stays generic in the number type (power fractions), so it can
  enter `optimize_pointing` and the Next 10 fits as parameters (v_ref, τ_b).
- Optional stage: a late "ridge stop", where inner power not absorbed is deposited at
  r ≈ 0.73 R_hohl (2.1/2.88) on an axisymmetric band after the ablator–wall collision time
  t_c. t_c is about 5.8 ns on basecamp and about 0.4 ns earlier per +20% r_cap.

**4. A time-integrated check against Callahan/Izumi.**
- Report the Callahan FoM and Izumi's FoM for each design next to Briareus's swing.

**Validation, in order of availability:**
- **(a)** Ralph's anchor: Briareus at the peak wall position for the 6.72 mm CH high-foot
  geometry should need CF ≈ 18.5% for zero flux P2.
- **(b)** Izumi 2018/2021 HDC SymCap set (N170427/N171010 basecamp, N171212, N180102,
  N180103, N180117, N180705, N180917): x_bub(t) trajectories 1.2–2.4 mm over 4.5–6.5 ns,
  thin-wall inner-signal time histories, and core P2 −18 to +21 µm. Seven shots with
  single-parameter changes test v_b scalings, clipping and the CF→P2 conversion (2.7 µm per
  1%, or 270 µm per unit CF).
- **(c)** Young 2020's static-vs-dynamic test reproduced in Briareus (85% radii against
  keyframes): mode 1 within about 0.14% at peak.
- **(d)** Callahan's 39-shot curve: P2 = 0 at FoM 1.2, σ = 8 µm.
- No published, numerically tabulated **flux P2/P0(t) swing curve** was found. All swing
  curves are HYDRA outputs in figures.

## Open questions

1. A digitizable P2/P0(t) flux swing for a HYBRID-E or N221204-class shot. Kritcher 2021
   Fig. 7 and 2024 Fig. 4 are figures only. 2DConA in-flight P2 vs radius for N210808 and
   N221204 was not found in open text.
2. The Landen 1999 analytic time-dependent P2 model equations (paywalled) and the Suter 1994
   spot-motion rates. Is there an open preprint (UCRL-JC)?
3. The Izumi 2021 Eq. 7 bubble-density exponent, and the sign of the ρ₀^0.13 FoM term
   (garbled extraction). The journal PDF would settle both.
4. Measured effective LEH radius vs time for low fill and HYBRID-E's 3.1 mm LEH (Chen 2020
   and Build-A-Hohlraum give it as figures). Whether DU vs Au walls change v_w.
5. MacGowan 2021 full text: which geometry and times it uses, and how well it predicts P2
   (not just mode 1).
6. Whether "MacLaren's empirical model" exists as such. Only the ViewFactor drive-deficit
   work (PRL 112, 105003, 2014) was found.
7. CBET coupling: does Δλ change the bubble's T_e and absorption, i.e. the swing slope?
   Kritcher 2021 says CBET shifts the swing but not its slope.
