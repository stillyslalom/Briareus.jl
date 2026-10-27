# Briareus → CompactLES: capsule-only implosions driven by a view-factor drive

Literature review, 2026-10-22. Confidence: **H** = primary text read (full text or accepted
manuscript), **M** = abstract, secondary summary, or an equation reconstructed from a garbled
PDF extraction, **L** = my inference. Full texts were read from the OSTI/arXiv copies cited
below (not kept in the repository).

## 0. What LLNL's "step two" actually is

- Kritcher et al., Phys. Plasmas 31, 070502 (2024), doi:10.1063/5.0210904, Appendix and
  Fig. 10 caption (**H**). Integrated HYDRA hohlraum+capsule runs, calibrated on tuning
  shots, give "the spatially, temporally, and frequency resolved radiation drive surrounding
  the capsule". That drive is "extracted and applied" to capsule-only runs. Those runs add the
  tent, fill tube, roughness, ablator and DT thickness non-uniformity, and drive modes 1, 2
  and 4. The physics list is detailed EOS, radiation and neutron transport, opacity models,
  electron–ion coupling, and "flux-limited electron heat transport with a limiter of 0.15".
  Laser multipliers are tuned to in-flight symmetry, velocity and shock velocities.
- The 0.15 is the **electron-conduction** flux limiter, not a radiation limiter.
  HYBRID-E: "Spitzer-Härm electron thermal conduction with a flux-limiter of 0.15 n_e T_e v_Te"
  (Kritcher et al., Phys. Plasmas 28 (2021; page not verified), doi:10.1063/5.0047841, OSTI 1818404,
  accepted manuscript, **H**). The same paper uses a 'high-flux model' with DCA NLTE for
  the hohlraum (**H**).
- Clark et al., Phys. Plasmas 26, 050601 (2019), doi:10.1063/1.5091449, OSTI 1756731 (**H**)
  states the method most plainly. "In these capsule-only simulations, the hohlraum is treated
  as a spherical boundary surface where a prescribed inward x-ray flux is applied." The flux
  is adjusted to four data sets:
  - VISAR shock timing (keyholes);
  - 1-D/2-D ConA in-flight velocity and shape;
  - the hard x-ray ("M-band") fraction from the View Factor platform;
  - the DT shot's bang time.

  Long-wavelength asymmetries come from separate 3-D HYDRA hohlraum runs with the
  as-delivered beam powers. The physics is multigroup diffusion, Monte Carlo burn products,
  and tabular EOS and opacity. Convergence was checked in zoning, energy groups and MC
  particles, but no group count is given. Clark calls the x-ray flux on the capsule "the
  largest source of uncertainty … the x-ray flux itself remains unmeasured" (**H**).
- Masse et al., Phys. Plasmas 26, no. 6 (2019; page not verified), doi:10.1063/1.5092827, OSTI 1885111 (**H**)
  gives the angular recipe: "A frequency dependent radiation drive that has been extracted
  from integrated hohlraum simulations is applied to the capsules only simulations. The
  asymmetries are applied to the radiation flux on the capsule at 100 µm from the ablation
  front", as F_X(θ,t) = Σ_l a_l(t) P_l(cos θ), with even l ≤ 6. The sign convention is that
  P2 > 0 means pole-hot, which gives an oblate shape.

**The angular part of Briareus's export already matches this interface:** a_l(t) Legendre
moments of the capsule flux, plus T_r(t). What Briareus lacks is the spectral content (see §1.4).

## 1. Stage 0 physics: ablation-pressure / rocket closures

### 1.1 Lindl scalings (Lindl et al., Phys. Plasmas 11, 339 (2004), Sec. VI.B, eqs. 6-1…6-13; full text, **H**, but equations garbled in extraction)
T_r is in heV and I₁₅ in 10¹⁵ W/cm². The equations hold for ablator albedo near zero. "When the
albedo of the capsule blowoff starts becoming significant, Eqs. 6-2 through 6-5 become time
dependent." They were validated on Nova at 200–300 eV.
- Eq. 6-1, rocket: V_imp = (P/ṁ) ln(m₀/m) = V_ex ln(m₀/m).
- Eq. 6-2: **P_abl ≈ 3 T_r^3.5 Mbar ≈ 170 I₁₅^(7/8) Mbar**.
- Eq. 6-3: **ṁ ≈ 3×10⁵ T_r³ g cm⁻² s⁻¹ ≈ 1.07×10⁶ I₁₅^(3/4)**.
- Eq. 6-4: V_imp ≈ 10⁷ T_r^(1/2) ln(m₀/m) cm/s. The exponent is reconstructed from
  P/ṁ = 3×10¹²/3×10⁵ T_r^0.5 (**M**).
- Eq. 6-5: V_abl = ṁ/ρ ≈ 2.3×10⁵ α^(3/5) T_r^(9/10) cm/s.
- Eq. 6-6: DT pressure at given density = α × Fermi pressure. The coefficient is garbled in
  my extraction; the widely quoted value is P ≈ 2.17 α ρ^(5/3) Mbar with ρ in g/cm³ (**M**).
- Eq. 6-13: V_imp ≈ 5.1×10⁵ α^(3/5) (R/ΔR) T_r^0.9 cm/s.
- Context from Sec. II.B (**H**). LASNEX 1-D capsule runs use P_N transport, in-line QEOS
  (More et al., 1988) and XSN average-atom opacities; "other opacity models such as OPAL …
  and STA … predict the same capsule performance". A non-Planckian, frequency-dependent
  source from 2-D hohlraum runs is used for RT studies. "Other than this effect, which can
  change short wavelength perturbation amplitudes by about a factor of 2 … the spectrum has
  little effect on target characteristics."

### 1.2 Carbon-specific, time-dependent (Trickey, Walsh, Pasley, High Energy Density Phys. 44, 100995 (2022), doi:10.1016/j.hedp.2022.100995, CC-BY, **H**)
These are 310 HYADES runs on 1 mm slabs under a constant T_r drive. Table 1 gives subsonic
scalings P_abl ∝ T_r^a t^b:
- CH: a = 3.6, b = −0.07, valid to 280 eV.
- Be: a = 3.3, b = −0.07, to 260 eV.
- **Carbon: a = 3.3, b = −0.09, to 360 eV.**
- Al and Ti (¾ solid density): a = 3.2, b = −0.17 and −0.21.

No prefactors are tabulated, so the absolute level must be calibrated (**H**). Above the
quoted temperature the ablation is trans- or supersonic and the scaling weakens.

### 1.3 Modern implosion metrics that a pressure model must hit (Lindl et al., "Key metrics of progress in the NIF ignition implosions…", Phys. Plasmas 33, 050503 (2026), doi:10.1063/5.0316678, OSTI 3367337, accepted ms, **H**; equations reconstructed from text → **M**)
- NIC point design: V = 370 km/s, adiabat α = 1.4, **P_abl = 140 Mbar**, M_imp = 390 µg
  (**H**).
- ITF ≈ 3.25 (M_imp/390 µg)^1.15 (α/1.4)^−1.9 (V/370 km/s)^8.5 (P_abl/140 Mbar)^0.88 (**M**).
- Rocket: V ~ [((Z+1)/A) T_r,eff]^0.5 ln(M₀/M_r). "T_r,eff … is reduced from the source
  temperature by the ablator albedo" (**H** text, **M** form). Max useful V ≈ 425 km/s at
  300 eV.
- Drive must stay on until ~0.5 ns before peak velocity (~0.7–0.8 ns before peak
  compression); this interval is the "coast time" (**H**).
- T_i(no-α) ≈ 4.0 keV (ρR_tot/1.4)^0.182 (V/370)^1.063 and ρR_hs ≈ 0.31 (…)^0.97 (V/370)^1.24
  g/cm² (**M**).
- **Energy accounting for HyE (**H**).** Energy in the imploding mass grew from 19 to 33 kJ.
  Coupling efficiency (laser → imploding KE) rose from 1.1 % to 1.6 %. "Hydrodynamic
  efficiency of about 14 % so the capsule absorbed energy varies from ~210–235 kJ",
  i.e. ~11 % of laser energy. HyE capsules absorbed ~210 kJ against ~140 kJ for the
  HDC (910 µm) capsules. Error bars are 3 % on V and 5 % on remaining mass. The inferences
  use a "power law approximation of the Rocket model".
- Hurricane, "How ignition and target gain > 1 were achieved", High Energy Density Phys.
  53, 101157 (2024), doi:10.1016/j.hedp.2024.101157, OSTI 2472991 (**H**):
  - Eq. 1: C m_shell v_imp²(1 − f²) ∼ (ρR)_hs T_hs ∼ p_hs R_hs ∼ R_hs², with C ≈ 0.66 the
    fraction of fuel KE into the hot spot.
  - f² = 1 − 1/[⟨(ρR)⁻¹⟩⟨ρR⟩] is the residual KE.
  - Coast/free-fall condition: p_hs(R_pv + Δ) ≈ p_abl(t_coast − R_pv/v_imp).
  - Optimum p_hs ∼ v_imp³ T_r^(7/5)(t_pv).

  These are directly checkable against a stage-0 run with a decaying P_abl(t).
- Hurricane et al., PRL 132, 065103 (2024): **not read** (paywalled; OSTI 2283993 has no
  full text). A search-engine excerpt says the ablation front drives "roughly 150 Mbar" (**M**).

### 1.4 Spectral issue: M-band
- Dewald et al., Phys. Plasmas 27, no. 12 (2020; page not verified), doi:10.1063/5.0027467, OSTI 1764324, accepted
  ms (**H**):
  - "15 % of the Au hohlraum spectrum represents M-band radiation > 1.8 keV" at peak.
  - M-band is pole-hot (P2 > 0) because the outer beams' wall irradiance is 5× the inner
    beams'.
  - At the Mo sphere, P2/P0 swings by d(P2/P0)/dt ~ 0.2/ns and reaches +0.15 at end of drive.
  - The inferred wall M-band P2/P0 is ≈ +0.53 (smoothing factor 0.28 at τ = 2).
  - "until recently high-resolution capsule only hydrodynamic simulations assume for
    simplicity that the M-band has similar drive symmetry to the softer x-ray drive that is
    typically only 1 % asymmetric".
- HDC capsules carry a ~15–20 µm buried layer of 0.3–0.5 % W "for shielding the inner undoped
  C from multi-keV X-rays" (Lindl 2026, **H**).
- Lindl 2004, Nova HEP4 (**H**): the 2–4 keV Au M-line fraction is the source of capsule
  preheat. A 3× underestimate of the hard fraction cut yield by only 30 % for 1.3 % Ge-doped
  CH.
- When it matters (**L**):
  - Stage 0: not at all, since a pressure BC carries no spectrum.
  - Stages 1–2: the M-band sets the inner-ablator preheat, the ablator–ice Atwood number, and
    a separate pole-hot P2 that differs from the gray drive's P2.
  - A gray drive with the right a_l(t) should still give v_imp and bang time to the few-%
    level (consistent with Lindl's "little effect" remark), but not interface stability.

## 2. Stage 1–2 physics: radiation diffusion, opacity, boundary condition

- **Gray or multigroup flux-limited diffusion.**
  - Equation: ∂E_g/∂t − ∇·(D_g∇E_g) = cκ_P,g(B_g − E_g), with D_g = cλ(R)/κ_R,g and
    R = |∇E|/(κ_R E).
  - Levermore–Pomraning limiter: λ = (coth R − 1/R)/R, or the rational (2+R)/(6+3R+R²).
    Original reference: Levermore & Pomraning, ApJ 248, 321 (1981) (form **M** via the
    survey arXiv:2401.05138 and others; citation **M**).
  - Larsen-type limiters: Olson, Auer & Hall, JQSRT 64, 619 (2000) (**M**, not read).
  - Implicit (backward Euler) treatment of the stiff matter–radiation coupling is standard.
    The time step otherwise falls to ~1/(cκ_P) in cold dense ablator (**L**).
- **Worked gray-diffusion capsule test:** Shestakov, Milovich, Kaiser & Prasad, UCRL-JC-131429
  (1998, NECDC), UNT digital library ark:/67531/metadc676941 (**H**). The setup:
  - D gas, ρ = 10⁻³ g/cc, out to 1.00 mm;
  - D fuel, 0.21 g/cc, 1.00–1.10 mm;
  - Be ablator, 1.85 g/cc, 1.10–1.21 mm;
  - SESAME EOS, initial T = 7 eV, Lee–More conductivity;
  - Rosseland opacity in D, Planck opacity in the coupling term;
  - hydro BC p = 58.22 GPa;
  - radiation BC **E_r − (2/c) F_r·n̂ = E_s** with E_s = a T_r⁴ and T_r = 160 eV held
    constant. This is the Marshak condition: the inward half-range current equals the
    incident σT_r⁴.

  Reported at t = 8 ns: T ≈ 0.716 keV near the origin and T_r ≈ 158 eV outside the
  ablation front. Inside the front T_r ≈ 18–30 eV, and ≈ 75 eV in the center. This is a
  public, fully specified problem with an imposed drive. It is the only one I found.
- **Angle-resolved BC (**L**).** Set E_s(θ,t) = 4F_inc(θ,t)/c from Briareus's
  F_inc = Σ a_l(t)P_l. For multigroup, split by a Planckian at T_r(t) plus an M-band
  component carrying its own a_l^M(t).
  - Pitfall 1: applied at a fixed outer radius through optically thin blowoff, diffusion
    smears the angular structure, and geometry smooths it by the case-to-capsule ratio
    (Dewald Fig. 6).
  - That is why Masse applies the asymmetry 100 µm outside the moving ablation front.
  - Pitfall 2: Briareus computes the flux on the initial capsule radius. Its
    flux-at-radius-R(t) must match where the BC sits (ROADMAP, "The P2 swing").
- **Opacity data access** (**M**):
  - LANL TOPS web (aphysics2.lanl.gov) returns gray Rosseland/Planck means by default, with
    multigroup optional; frequency-dependent runs are limited to 6 (T,ρ) points; data dated
    2015-08-19. Usable for C, W-doped C, DT and CH mixtures.
  - OPAL is LLNL, astrophysical (**L**).
  - CompactLES already parses SESAME and IONMIX tables (`SesameTable`, `IonmixTable`), but
    "a solver does not evaluate these tables" (CompactLES docs, Physics models, **H**).
- **EOS** (**M** unless noted):
  - Gaffney et al., "A review of equation-of-state models for ICF materials", HEDP 28, 7–24
    (2018), doi:10.1016/j.hedp.2018.08.001. Models differ by 5–10 %, more where ionization
    and dissociation occur; D is especially uncertain.
  - Clark 2019 (**H**): DT EOS choice shows "remarkably little sensitivity", but HDC EOS
    models "do show some sensitivity".
  - QEOS: More, Warren, Young & Zimmerman, Phys. Fluids 31, 3059 (1988), analytic TF +
    Cowan ions + bonding correction (**M**).
  - FEOS (Faik et al.; QEOS-based C++ table generator, Mendeley Data 6vjsv6v48p) handles
    mixtures (**M**).
  - Carbon multiphase EOS: Benedict et al., PRB 89, 224109 (2014), arXiv:1311.4577,
    3–20 g/cc (**M**).
  - SESAME needs a signed LANL agreement (**M**).
  - HYBRID-E chose a first shock > 12 Mbar "to avoid refreeze of the carbon", so HDC melt
    matters to the foot (Kritcher 2021, **H**).
  - Minimum viable EOS (**L**): ideal gas cannot hold a DT adiabat of 1.4–2.5. The floor is
    a cold Fermi term P_c = 2.17 ρ^(5/3) plus thermal ideal gas with a Thomas–Fermi or Saha
    Z̄(ρ,T). QEOS/FEOS tables are the first real target.

## 3. Validation numbers available in public text

| Quantity | Value | Source |
|---|---|---|
| v_imp, ablator mass remaining, capsule absorbed energy (sim, N230729 base 2.05 MJ, +0/+4/+10 µm HDC) | 386/383/374 km/s; 5/5.9/6.6 %; 245/257/274 kJ | Kritcher PoP 31 Table III (**H**) |
| HYBRID-E Table I (480 TW) | v 360–>400 km/s, mass remaining <5–6.5 %; first shock ~50 TW ≈ 150 eV ≈ >12 Mbar; adiabat ≈ 2.5 at peak velocity, ~5 at bang | Kritcher PoP 28 (2021) (**H**) |
| HyE energetics | KE_imp 19–33 kJ; η_hydro ≈ 14 %; E_abs 210–235 kJ ≈ 11 % of laser | Lindl PoP 33 050503 (**H**) |
| 2DConA surrogates, data (sim) | N220618: v 392±24 (401) km/s, t_stag 9.170±0.061 (9.214) ns. N221010: 381±15 (408), 9.230 (9.260), mass rem. 7.2±1.0 (7.52) %. N240304: 355±5 (344), 9.876±0.090 (9.876), 8.8±1.0 (10.20) %. Shell P2 9–27 µm, HS P2 5.7–25.4 µm (see table) | Ruof et al., PoP 33, 052706 (2026), OSTI 3367476, Table IV (**H**; row labels offset by one in extraction, re-aligned) |
| Capsule-only (drive-tuned) vs N170601 | bang 8.19 vs 8.16±0.03 ns; T_ion 4.97 vs 4.93±0.2 keV; DSR 3.11 vs 3.10±0.2 %; P0 30.9 vs 30±1 µm; Y 1.31 vs 1.46e16 (3-D) | Clark PoP 26 Table I (**H**) |
| Mass remaining N210808 → N221204/N230729 | 3.8 % → 5.7 % (capsule-only, surrogate-tuned) | Ruof 2026 (**H**) |
| Coast time (pre-N210808 HyE) | ~1.4 ns from 90 % peak T_r to peak compression, ~1.6 ns from laser-off | Zylstra et al., PRE 106, 025202 (2022), OSTI 1959535 (**H**) |
| Yield scaling | Y ∝ η^3.83 p_abl^0.64 v_imp^4.47 α_if^−1.44 S^4.67 (exponent pairing reconstructed) | Ruof 2026 Eq. 1 (**M**) |
| Kritcher PRE 106, 025201 and PRE 109, 025204 tables | not read (paywalled; OSTI 1880163 / 2474067 no full text) | — |

The bang time of capsule-only runs is an **input** (the drive multipliers are tuned to it,
Clark 2019, **H**). For a Briareus-driven run it is the first **prediction**.

**Open-source cross-checks.**
- MULTI-IFE (Ramis & Meyer-ter-Vehn, CPC 203, 226 (2016), doi:10.1016/j.cpc.2016.02.014) is
  open 1-D Lagrangian code with multigroup transport, DT burn and α diffusion. I found no
  published indirect-drive capsule test with its drive (**M**).
- I found no FLASH, Castro, xRAGE or Pyranda capsule test with a tabulated T_r(t) (**M**,
  negative search).
- The 1998 ICF3D problem above is the one fully specified imposed-drive test.

## 4. Applying an angle-resolved, time-dependent drive in 2-D

- **Green's-function response** (Masse 2019, **H**):
  - Shape response R_l(t) = Σ_i f_l(τ_i) G(t, τ_i), built from 45 HYDRA impulse runs over an
    8 ns drive. It is linear; Clark 2019 says "this linearity assumption is very accurate".
  - Peak sensitivity is near 6 ns, roughly peak velocity.
  - **Sign inversion:** a pole-hot P2 before ~4.5 ns gives an oblate hot spot but a prolate
    one later. This explains why hot spot and shell are non-conformal, the "swing".
  - Proposed fix: an early positive P2 to cancel a late positive P2, giving a round,
    swing-less hot spot and ρR. This was confirmed in direct HYDRA runs.
  - The inverse fit tracks integrated runs until ~4 ns, then deviates (N150105). It matches
    them for a 0.3 mg/cc fill (N150607).
- **Time-dependent mode-1 swing piston model** (Hurricane et al., Phys. Plasmas 29 (2022; page not verified),
  doi:10.1063/5.0067699, OSTI 1840112, **H** text; page number **M**):
  - Shell asymmetry at stagnation dominates, and a weighted harmonic mean of shell ρR is
    what matters.
  - Swing matters when P_stag/P_pv ~ 10 but is negligible at ~100.
- Implication for Briareus (**L**). P2(t) from `drive_history` convolved with G gives hot-spot
  P2 without a hydro code, a "stage −1". A CompactLES 2-D run's job is to compute G
  for the Briareus capsule and check the linear prediction.
- **Pitfalls of capsule-only with an imposed P2(t)** (from the above, **H**/**L**):
  - The drive is calibrated, not predicted; LLNL tunes multipliers per design (Ruof 2026).
  - The M-band has its own, larger P2 (Dewald).
  - The asymmetry must be applied near the ablation front, not at a distant wall.
  - The capsule–wall stagnation plasma at the equator (Dewald) and bubble or LEH geometry
    are absent.
  - A gray view-factor drive gives no spectral angular dependence.

## (b) Staged coupling plan

**Stage 0: pressure-driven surrogate (1-D spherical, then 2-D r–θ).**
- Briareus exports T_r(t) and a_l(t)/a_0(t), l ≤ 6. CompactLES applies
  P_abl(θ,t) = P₀ T_r,eff(θ,t)^3.3 (t/t₀)^−0.09, with carbon exponents from Trickey and
  T_eff⁴ = T_r⁴(1 + Σ a_l P_l/a_0).
- Mass loss: ṁ = 3×10⁵ T_r³ removed at the outer boundary, or more simply a fixed-mass
  "rocket" shell.
- Calibrate P₀ to HYBRID-E: first shock >12 Mbar at ~150 eV, P_abl ~140–150 Mbar at peak.
- Validate on:
  - v_imp 380–400 km/s with 5–6.6 % mass remaining;
  - coupled KE 19–33 kJ;
  - Hurricane's free-fall condition p_hs ≈ p_abl(t_coast − R_pv/v);
  - Masse-type G(t,τ) sign inversion near 4.5 ns.
- It **cannot** validate: shock timing (sensitive to the foot ablation physics), adiabat,
  absorbed energy, M-band effects, or absolute bang time to better than ~100 ps (**L**).
- CompactLES has no moving-boundary or Lagrangian ablation front. The practical Eulerian form
  is an `add_source!` energy deposition of the absorbed flux (1 − α)σT_r⁴ in a thin layer
  near the density gradient's maximum. The rocket then emerges from the hydro (**L**).

**Stage 1: gray flux-limited diffusion + table EOS (1-D, then 2-D).**
- Implicit gray radiation energy equation coupled to T_e (one-T to start), with κ_R and κ_P
  from TOPS for C (+0.4 % W) and DT.
- Marshak BC E − (2/c)F·n̂ = 4F_inc(θ,t)/c.
- Spitzer–Härm conduction flux-limited at f = 0.15.
- QEOS/FEOS tables for C, DT and He.
- Validate: first the 1998 Be/D 160 eV problem, then Lindl's 3 T^3.5 Mbar and 3×10⁵ T³ in
  planar C slabs, then Kritcher Table III (v_imp, mass remaining, E_abs ≈ 245 kJ) with a
  Briareus drive scaled to Dante T_r.

**Stage 2: multigroup + M-band.**
- Briareus spectral radiosity (ROADMAP "Spectral radiosity") supplies a Planck component plus an Au M-band
  component (~15 % > 1.8 keV, pole-hot, its own a_l).
- Multigroup diffusion; the group count is not public, so pick it by convergence (**L**).
- Two-temperature model with electron–ion coupling.
- Validate: W-layer preheat, Dewald's capsule-inner-surface P2 history, ablation-front growth
  vs HYBRID-E growth-factor curves, Ruof 2DConA trajectories (v, t_stag, shell P2), and the
  Clark N170601 observables as a full-chain check.

## (c) Interface sketch

**Briareus exports** (HDF5, or an Arrow/CSV pair; one file per run):
- Time grid t_k [ns], uniform Δt ≤ 50 ps (Masse used 45 impulses over 8 ns ≈ 180 ps; finer
  for the rise).
- T_r(t) [eV] (Dante-equivalent), and F_inc,0(t) [TW/mm²] = capsule-averaged incident flux.
- a_l(t)/a_0 for l = 1…8 (Legendre, axisymmetric), and real Y_lm for 3-D or m ≠ 0 diagnostics.
  Convention: P2 > 0 = pole-hot, θ from +z (hohlraum axis).
- The radius R_cap at which the flux was evaluated. Optionally a table over R ∈ [R_abl,min,
  R₀] so CompactLES can apply the flux at the moving front.
- Later: spectral split into groups ν_g, F_g(t), and a_l,g(t), at minimum Planck + M-band.
- Provenance: hohlraum, pulse, albedo model, mesh, and git hash.

**CompactLES must add:**
1. `AblationPressureSource` / `AbsorbedFluxSource <: Source`, via `add_source!(src, solver,
   dQ, Q, t)`, which receives t. It interpolates the Briareus table in t and θ.
2. Radiation module: gray, then multigroup, implicit FLD (Levermore–Pomraning limiter), with
   a Marshak face BC. This is a `correct_flux!`/`correct_rhs!` boundary. The docs list no
   time argument for BC hooks, so time must come through the solver state (**M**).
3. Electron conduction: implicit, nonlinear Spitzer–Härm with flux limiter f.
4. Table EOS: SESAME/QEOS-backed `EOS` subtype implementing `recover_primitives!`,
   `eos_phi`, etc. The readers exist but are not solver EOS.
5. Opacity tables (κ_R, κ_P vs ρ, T, group).
6. Spherical 2-D r–θ (PoleBC at both θ ends exists; AMR excludes spherical). Otherwise
   axisymmetric r–z with AMR, where the drive must be mapped to the sphere.
7. Diagnostics: Legendre moments of ablation-front and hot-spot contours, shell ρR(θ),
   v_imp, mass remaining, bang time (burn rate).

## (d) Open questions

1. How large is the P_abl prefactor error from using the incident T_r rather than
   T_r,eff (reduced by ablator albedo)? Briareus's capsule albedo model sets this.
2. Does a pressure or absorbed-flux source in an Eulerian compact-FD code survive the
   ~10³ density jump at the ablation front without artificial-property smearing that
   changes ṁ?
3. What group count and spectral binning does HYDRA use for HDC capsules? It is not in any
   text I read. Clark says only "converged".
4. Is P_abl ∝ T^3.3 t^−0.09 (constant-T slabs) valid through a shaped NIF pulse? Trickey
   warns real drives are neither as fast nor isothermal.
5. Where should the angle-resolved BC sit: a fixed radius with geometric smoothing, or
   tracking the ablation front as in Masse? This determines whether Briareus must supply
   F(θ; R) for several R.
6. What is the M-band a_l(t)? It needs a spectral/non-LTE source in Briareus, which emits
   gray only today.
7. Kritcher PRE 106/109 and Hurricane PRL 132 tables were not accessible. The numbers in §3
   are from PoP 31, PoP 28, PoP 33 (Lindl, Ruof), PoP 26 and HEDP 53 only.
8. Lindl 1995 (PoP 2, 3933), Betti/Zhou scaling papers, Saillard (Nucl. Fusion 46, 1017
   (2006), cited by Masse for hot-spot shape) and the Springer/Patel hot-spot models were
   not read.
