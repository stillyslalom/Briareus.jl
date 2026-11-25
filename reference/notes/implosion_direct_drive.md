# Research notes: Implosion metrics, hydrodynamics hand-over, direct drive and the corona

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

## Implosion metrics (2026-10-22)

`src/implosion.jl`, the implosion-metrics tutorial and `test/implosion.jl`. The sources
are summarized in `reference/litreview/asymmetry_transfer.md`. Confidence: **H** =
primary text read, **M** = abstract or secondary, **L** = inference or arithmetic.

- **`HOTSPOT_VELOCITY_PER_PERCENT` = 100 km/s per % of mode 1** (`a₁/a₀`, the norm of
  the l = 1 vector, burn-off, rise and peak). *Engineering default.*
  - Young et al., Phys. Plasmas 27, 082702 (2020) **H**: N170601 0.44% → 37.6 km/s
    (3D HYDRA), 44.7 km/s (VF + G); linear below ~0.75%.
  - MacGowan et al., HEDP 40, 100944 (2021) **M**: ~100 km/s per %.
  - Kritcher et al. 2024 **H**: N230729 ~0.7% → 70–75 km/s burn-off (~105 per %, **L**).
  - Spears et al. 2014 **H**: 2% → 120 km/s (saturated; the linear rule gives 200).
  - Platform specific; Briareus gets 44.0 and 70.0 km/s for the first two (test suite).
- **Direction**: the hot spot moves away from the high-flux side (Casey et al., PRL 126,
  025002, 2021, **H**). Briareus's real Y_lm have no Condon–Shortley phase, so
  `v = −k (c₁₁, c₁,₋₁, c₁₀) √3/c₀₀` in chamber x, y, z. Young's Eq. 9 signs
  (V_x = +G c₁₁) follow their own associated-Legendre convention, which was not
  checked; the physical rule (away from the bright side) is what is implemented. **L**
- **Piston model** (Hurricane et al., Phys. Plasmas 27, 062704, 2020, **H**):
  f = v_HS/v_imp; RKE/KE = c f², c = 0.66; no-α yield (1 − f²)ᵃ, a = 3.3; α factor
  exp(−q χ^1.2 f²), q = a − 1 (Eq. 21, 26). Casey et al. PRE 108, L053203 (2023) **H**
  write q = 1.2(a − 1) and a = 6 − 1.5 ln T; both are offered (`q`).
  - `vimp` = 380 km/s default: Hurricane's Table I value. *Engineering default.*
  - Table I is reproduced at χ = 2 (f to two digits, v_HS 29.4/42.0/52.2/61.1 against
    30/40/50/60, rounded); at χ = 1.5 the 5% row would give f = 0.09. So the table's
    "near ignition" χ is ≈ 2. **L**
  - Casey 2021's "~22% no-α yield at 100 km/s / 400 km/s" gives 19% with
    (1 − f²)^{10/3} at f = 0.25. Small gap, not resolved. **L**
- **Expected yield over errors**: `E[exp(−β f²)]` for Gaussian v in closed form,
  `det(I + 2bΣ)^{-1/2} exp(−b μᵀ(I + 2bΣ)⁻¹μ)`, b = β/v_imp². This uses
  (1 − f²)ᵃ ≈ exp(−a f²) (0.7% high at f = 0.25). *Derived.* Checked against a Monte
  Carlo of full evaluations (3%; 5% with α).
- **Burn amplification** (Nora et al., Phys. Plasmas 30, 092705, 2023, Eq. 4, **H** via
  the lit review): 0.324 ln(Y_amp + 2.239) + 0.653 for Y_amp ≤ 40.
  - It is ~1.8 at Y_amp = 30, against Kritcher 2024's burn-on/off ratios of 1.27–1.39.
    Use it as an upper estimate.
  - The Y_amp 80–500 branch as transcribed (13 ln(Y_amp + 1585) − 95) is 1.43 at 80,
    below the low branch's 1.87 at 40. *Resolved 2026-10-24* ("Young 2020 validation"
    below): Eq. 4 gives rounded coefficients; Fig. 6's (13.105, 1585.225, −95.539) are
    implemented on 80–500, with 40–80 refused.
  - The shape term (Eq. 2) as transcribed does not give the stated 0.8–1.2 range for
    |P2| < 0.4 (0.19 at −0.4). *Resolved 2026-10-24*: that range is Table I's (3D);
    Eq. 2 is the 2D fit. Both are in `shape_amplification`.
- **Thickness nuisance in `infer_velocity_errors`** (Casey 2021, **H**):
  - 110 km/s per % of thickness mode 1 (simulated; 140 ± 30 measured), toward the
    thicker side.
  - The prior is 0.2% per axis by default. *Derived* from ~0.15 µm metrology on a
    78 µm shell (**L**).
- **nToF covariance**: ±14 km/s, 13° θ, 29° φ (Young 2020 **H**), in the measured
  vector's radial/polar/azimuthal frame.
- **Picket P2 limit** 7.5% (Milovich et al., Phys. Plasmas 23, 032701, 2016, **H**).
  The paper does not say whether this limits the picket's average or its peak.
  `swing_metrics` judges the average and reports both. **L**
- **Comparisons with Young 2020 that do not match** (tutorial hohlraum: 10.1 × 5.75 mm,
  h 0.6 / capsule_h 0.2, wall albedo 0.85, no LPI; *derived*):
  - A dropped inner quad gives 190–217 km/s and an outer quad 139–150 km/s, a ratio of
    1.3–1.6 against Young's ~2 (N170827).
  - The whole-target offset gives 0.99 km/s/µm vertical (Young 1.3) and 0.60 km/s/µm
    along x (Young 0.1 horizontal).
  - The geometry and LPI differ from N170827. A horizontal offset may couple more
    strongly in a static single-albedo model; open. *Resolved 2026-10-24* ("Young 2020
    validation" below): a setup difference (on Young's static N170827 setup Briareus
    gives a ratio of 2.2 against 2.07 and offsets of 0.09/1.30 against 0.1/1.3), plus a
    fixed error in the mode analysis of offset capsules.
- **Not implemented** (no published values):
  - Young's G₁(t) and Masse's P2 kernels, which are printed only as figures. Briareus
    ships `flat_kernel` (flat over a window, the static sensitivity; *derived*) and
    takes user kernels.
  - Time-dependent error budgets and multi-shot trending remain (ROADMAP "Implosion metrics and UQ").

## Capsule drive hand-over and CompactLES stage 0 (2026-10-24) *derived*

ROADMAP "Briareus → CompactLES", stage 0. Sources: `reference/litreview/hydro_coupling.md` (all
literature numbers below are from there, with its confidence marks). Confidence here:
**L** our construction, *calibrated* as marked.

**Hand-over (`capsule_drive`).** The drive is the `moments` of the history (Legendre
about +z, about the capsule center, P2 > 0 pole-hot), interpolated linearly in aₗ and
divided by a₀; `flux = a₀`, `Tr = radiation_temperature(a₀)`, exactly the history's at
its times. Default grid 50 ps (hydro_coupling.md §(c)). The radius is the observed
`SphereSurface`'s (the initial capsule; the flux is not re-evaluated at the moving
front, pitfall 2). Tested: (1 − α) a₀ A equals `energy_balance`'s capsule absorption
to 2e-3 on the coarse test hohlraum (solid-angle vs area weights of the facets); beams
landing on the caps give P2 > 0, toward the waist P2 < 0.

**Extension side.** `CompactLESBriareusExt` in CompactLES (weak dependency on Briareus):
the glue is CompactLES physics (a source term and its calibration, later a Marshak
boundary) and evolves with that solver; Briareus's side is a plain struct. Load cost
is symmetric (an extension loads only with both), but this way Briareus never names
MPI/KernelAbstractions even as a weak dependency, and CompactLES's `AblationDrive`
takes plain vectors, so it does not need Briareus either (**L**).

**Stage-0 model (CompactLES `AblationDrive`).** P_abl = P₀ T^3.3 (t/1 ns)^-0.09 Mbar
(carbon, Trickey 2022), ṁ = 3e5 T³ g cm⁻² s⁻¹ (Lindl 2004 Eq. 6-3), T in heV,
T_eff⁴ = Tr⁴(1 + Σ aₗ/a₀ Pₗ) on r–θ grids. The front is the ρₐ(−∂ρₐ/∂r)⁺-weighted
radius per radial line; mass is converted from ablator to a driver gas in a 5 µm
Gaussian layer, and the driver gas's pressure is relaxed to P_abl (τ = 10 ps): a
piston at P_abl. A body-force form (force density on the layer) lagged the drive
(18 Mbar in the shell against 117 applied at 5 ns) and failed; kept for experiments.

**First run (*calibrated*: P₀ only).** Briareus: N221204-class hohlraum (11.24 × 6.4 mm,
3.1 mm LEH, 1.135 mm capsule), delivered pulse by eye, `hybride_mock_pointing(nif_quads())`,
static gold wall, no CBET, h 0.4/0.2: capsule Tr peaks at 297 eV (300 with the
calibrated two-zone CBET; P2/P0 −0.09 vs −0.25 at the peak, neither with wall motion).
P₀ = 4.81 Mbar puts the peak at 145 Mbar; 17 Mbar where Tr first reaches 150 eV
(> 12). CompactLES 1-D, 1400 points to 1.4 mm (1 µm), ideal-gas fully ionized C
(3.5 g/cc, 1.050–1.135 mm), DT ice (0.255, 0.985–1.050), DT gas 0.3 mg/cc, He driver
1 mg/cc, initial 0.1 Mbar:

| | stage 0 | target |
|---|---|---|
| shell pressure at the front | P_abl within 1% from 0.75 ns | — (by construction) |
| peak DT velocity | 315 km/s at 8.7 ns | 380–400 km/s |
| ablator remaining at peak velocity | 19.6% | 5–6.6% |
| ablator + DT KE at peak velocity | 36 kJ | 19–33 kJ |
| stagnation | ≈ 9.3 ns | (bang 9.53, not a prediction) |

Reading (**L**): with the pressure right, the rocket is too heavy: Lindl's ṁ (a Nova
fit for low-albedo ablators) removes 80% of the HDC ablator where HYBRID-E removes
~94%. Calibrating ṁ to the mass remaining fixes v through the rocket equation; it is
one more knob and the obvious next step, and the KE comparison is loose (our KE counts
the remaining ablator). Numerics: the run takes ~30 min on 4 threads (timestep set by
the hot driver gas, CFL 0.1 after the automatic retries in the foot); a 2 µm grid
(layer < 3 cells) fails in the foot.

**Update: robustness and the ṁ calibration (2026-10-24, CompactLES `briareus-drive`;
drive from Briareus 785684d, cached).** Two defects of the first run, both fixed in
CompactLES. (1) The foot failure was the driver gas, not the layer. Held at P_abl, the
1 mg/cc He is a 2–7 keV ideal gas against a 0.05 eV ablator. The contact cells had
negative internal energy from 0.1 ns and collapsed the timestep near 0.5 ns, at any
CFL. A 10 mg/cc driver (≈ 1% of the ablator's mass) runs the foot at CFL 0.4 with
no retry. The whole implosion is now one `run!` call of 5–8 min on one thread, with
1–2 automatic retries when the first shock reaches the origin (7.3 ns). (2) The
19.6% was mostly ablator stranded outside the front. A partly converted layer gets no
push in the uniform-pressure driver, so the accelerating shell left it behind (22% of
the ablator at 6 ns). It then pulled the front marker off the shell, and ablation of
the shell stopped. Now ablator beyond the front is converted first, and the piston
acts only at and beyond the front. Results (1 µm grid; ablator remaining = ablator
inside the front at peak DT velocity), *calibrated*:

| ṁ prefactor A (ṁ = A T³) | v_DT peak | ablator remaining | KE (ablator + DT inside front) |
|---|---|---|---|
| 2.6e5 | 320 km/s | 7.4% | 26.5 kJ |
| **2.85e5 (adopted)** | **336 km/s** | **5.7%** | **25.8 kJ** |
| 3.0e5 (Lindl Eq. 6-3) | 346 km/s | 4.8% | 25.1 kJ |
| 3.8e5 | 406 km/s | 0.9% | 19.2 kJ |

Reading (**L**): 5–6.6% remaining needs A = 2.7–3.0e5, so Lindl's 3e5 T³ was right
and the earlier "too little ablation" came from the stranding. The T³ exponent is
kept: Lindl's I^(3/4) is T³, and P/ṁ ∝ T^0.5 with P ∝ T^3.3 would give T^2.8, a
difference the 290–297 eV-dominated ablation does not resolve. The velocity check
fails without tuning: 328–344 km/s across the band against 380–400. Reaching 380–400
needs A ≈ 3.5–3.75e5 and leaves 1.5–2.5%. Lindl's ideal rocket gives 395 km/s at this
mass ratio. The gap is in P/ṁ, so the likely causes are the ideal-gas carbon adiabat,
the 145 Mbar peak and the gray drive, not ṁ. KE 25–26 kJ is inside 19–33 kJ. The grid
is not converged: at 2 µm the same A gives 298 km/s and 7.3%, and the 0.5 µm run fails
at 5.5 ns. Read v to ±10%. The calibration belongs to the drive. Briareus 3f25f29
(late-time wall-loss history) peaks at 285 eV, and with P₀ recalibrated to 5.53 Mbar
the same A gives 307 km/s, 8.2% and 25.6 kJ (the CompactLES example, end to end).
r–θ run (aₗ to l = 8, 2 µm × 64 cells, annulus from 0.25 mm, stopped at 7.6 ns).
The pole-hot foot (a₂/a₀ +0.03…+0.15 to 4.7 ns) makes the shell oblate, to −0.45%
P2/P0 at 5.5–6 ns. The pole-cold peak (−0.09) reverses it, and the shell is prolate
from ≈ 7.1 ns (+0.5% at 7.6). This is a swing, but not yet Masse's G(t,τ) inversion,
because the P2 of the drive itself changes sign at 4.7 ns. Telling them apart needs impulse
runs to stagnation (**L**).

**Update: every m, and the r–z and 3-D metrics (2026-10-24; CompactLES
`briareus-drive`).** `capsule_drive` now carries `ylm`, the history's real Y_lm
(`drive_history(...; ylm = true)`) as Schmidt semi-normalized ratios
bₗₘ = √((2l+1)/4π) cₗₘ/a₀, so bₗ₀ = aₗ/a₀ (the m = 0 rows are set to `modes`; the
history's Y_lm quadrature is one sub-triangle level coarser, < 1e-3 of a₀ apart) and
√Σₘ bₗₘ² = `mode_amplitude`. Tested: a known field's ratios recovered to its
element-constant quadrature, round trip, and a 45°-rotated three-quad pattern on a
45°-symmetric mesh rotates (cₗₘ, cₗ,₋ₘ) by mα to 1e-7. CompactLES's front search and
layer normalization are now sums over equal-solid-angle bins about the capsule center
with dV/R² weights (per steradian, so the per-line form is reproduced to 1e-14 in
1-D), which runs on r–z (axis) and 3-D Cartesian grids, and on refined levels through
a level-wide hook before each level's right-hand side. Validation (**L**), flat
250 eV drive, P₀ and ṁ not retuned, 2-cell layer: r–z 10 µm against 1-D 10 µm to
1.5 ns, front within 0.45 µm, DT shell within 0.15 µm, ablator left equal to 4 digits,
momentum within 1% (a tiled 10 µm band likewise to 0.6 ns, then failed); a b₂₀ = 0.2
drive gives the same front P2 on r–z as on r–θ to 10–15% (−0.55 vs −0.50% at 1.5
ns); a 3-D octant at 20 µm reproduces 1-D (front +0.5…+1.7 µm, ablator to 1e-4), and
its response to b₂₂ = 0.2 equals its response to b₂₀ = 0.2 to 1–4% (−0.123 vs
−0.122% of R_f at 0.6 ns), with ≤ 0.014% in Y₄₂ (cubic-grid coupling) and ≤ 0.003%
elsewhere. Two limits were met: a band's coarse level takes the timestep of its covered nodes
(8× the uniform grid's steps), and tiled levels failed on this problem; a 1 µm 3-D
shell (~10⁹ cells uniform) waits on those.

**Update: the full-model drive and an HDC ablation-rate source (2026-10-25; Briareus
c54e366, CompactLES `briareus-drive-2`).** The drive is now `examples/n221204.jl`'s full
model (+bubble: `nif_uranium()`, Young's albedo on the K24 schedule, release 1.0,
calibrated two-zone CBET, measured backscatter, 11% degradation, wall, lips, bubble),
coarse mesh 0.8/0.4, steps of 0.1 ns from 0.2 ns to bang (9.53 ns), l ≤ 6 with every m.
The capsule Tr peaks at 306 eV at 8.0 ns (297 before) and is still 282 eV at 9.0 ns; P2/P0
+0.05 near 3.2 ns, −0.18 at 5.2 ns, +0.09 over 7.2–7.7 ns. Stage 0, 1 µm, ṁ = 2.85e5 T³,
P₀ recalibrated to 145 Mbar (*calibrated* as before; the checks are PRE 109 Table I,
HYDRA, N221204: v_imp 383 km/s, ablator remaining 5.7%, bang 9.61 ns simulated,
9.53 ± 0.07 measured):

| | stage 0 (new drive) | check |
|---|---|---|
| peak DT velocity | 339.5 km/s at 9.26 ns | 383 km/s |
| ablator remaining at peak velocity | 4.97% | 5.7% |
| KE (ablator + DT inside the front) | 24.0 kJ | 19–33 kJ |
| stagnation (min. DT mean radius, peak p(0)) | 10.22 / 10.20 ns | bang 9.53–9.61 ns |

Reading (**L**): the stronger, later peak moves peak velocity to 9.26 ns and stagnation
0.6–0.7 ns past bang; no burn or radiative loss, so stagnation is not bang, but a
0.6 ns lag is a drive-or-closure timing error, not a definition. Velocity is unchanged
within the ±10% grid band (336 → 340). Ablation-rate source independent of v_imp
(cite-only): Olson et al., PoP 18, 032706 (2011) (OSTI 1114987, abstract: HDC 2–9
mg cm⁻² ns⁻¹ at 160–260 eV, ≈ T³, planar Omega halfraum burnthrough, **H**); the HDC
points of SAND2007-5659C (OSTI 1147367, IFSA07; 10 points read from the slide's vector
plot, 175–255 eV, e.g. 8.3 at 245 ± 5 eV) give ṁ = A T³ with A = 4.9–5.2×10⁵
g cm⁻² s⁻¹ (free exponent 3.8), ±11% from Dante's ±3.5% (**M**); SAND2011-8633C (OSTI
1141159) uses 0.415 T³ mg cm⁻² ns⁻¹ (4.15×10⁵) for HDC and Ge-CH together. Lindl's
3×10⁵ is 1.4–1.7× below the HDC data. At those rates this shell over-ablates (stage 0,
1 µm, same drive): 4.15e5 leaves 0% and 5.0e5 fails at 7.33 ns with 0.11% left. The
grid study at 2.85e5 fails too: nx 700 and 2000 end on `dt_collapse`, only nx 1400
completes. A = 2.85e5 is kept as an *effective* rate, *calibrated* to the ablator remaining
(5–6.6%, PRE 109 Table I / HYDRA), so that number is a fit target, not a check; v_imp,
KE and bang time remain checks (decided 2026-10-26). Replacing it with a measured rate
is CompactLES ROADMAP A22 (rocket-relation pressure tied to ṁ, real carbon EOS, grid
robustness).

**Where the stagnation lag comes from (2026-10-27, L).** Eight variants of the run above
(same drive, 1 µm, nx 1400, one closure changed each), timed against K24's HYDRA radius
history (Fig. 1, the ablator–fuel interface; lag = stage-0 time minus K24's at the same
R/R₀, DT mean radius; keyhole N230110 shock timing, Ruof 2026 Table II, as a check only):

| case (change) | shock into ice | v_peak | ablator left | stagnation | lag at R/R₀ 0.95 / 0.85 / 0.50 / 0.20 |
|---|---|---|---|---|---|
| baseline | 4.80 ns | 339.5 km/s | 4.97% | 10.22 ns | +0.98 / +0.68 / +0.57 / +0.55 |
| p30 (foot P ∝ T^3.0, peak kept) | 4.48 | 324.1 | 5.58 | 10.22 | +0.79 / +0.54 / +0.52 / +0.54 |
| foot (picket/foot Tr × 1.15 to 3 ns) | 4.16 | 332.7 | 4.64 | 10.12 | +0.65 / +0.46 / +0.44 / +0.44 |
| g3p30 (p30 + γ_C = 3, stiff carbon) | 3.70 | 301.4 | 7.33 | 10.52 | +0.67 / +0.47 / +0.67 / +0.76 |
| footg3 (foot + γ_C = 3) | 3.42 | 300.9 | 9.47 | 10.44 | +0.48 / +0.38 / +0.62 / +0.72 |
| pk160 (P × 1.10) | 4.72 | 344.0 | 5.72 | 10.08 | +0.87 / +0.56 / +0.45 / +0.42 |
| sh03 (whole drive 0.3 ns earlier) | 4.56 | 339.4 | 5.06 | 9.92 | +0.67 / +0.37 / +0.27 / +0.25 |
| r33 (ṁ = 3.3e5 T³) | 4.86 | 371.1 | 2.41 | 9.84 | +0.94 / +0.60 / +0.37 / +0.28 |
| r35 (ṁ = 3.5e5 T³) | 4.86 | 386.1 | 1.45 | 9.70 | +0.92 / +0.56 / +0.29 / +0.16 |

In the keyhole check, the first shock enters the ice at 3.86–3.91 ns. The early lag (R/R₀ ≥
0.85) is therefore the first shock: every change that brings the shock to the keyhole timing
(stiff carbon, stronger foot) cuts it. But the stiff-carbon shells then accelerate more
slowly (300 km/s, 7–9% left) and the lag grows again to +0.7–0.8 ns by R/R₀ 0.2: fixing
shock timing alone moves stagnation later, not earlier. Only burning more mass closes the
late lag (r35: 386 km/s, lag +0.16, but 1.45% left against 5–6.6%), or a drive 0.3 ns
earlier (sh03), which Dante rules out (±6 eV over 4.5–7.5 ns). No single closure change
matches shock timing, v_imp, remaining mass and stagnation together: in a 1-D rocket model
velocity and remaining mass are tied by the rocket relation through ṁ and the exhaust
velocity (P/ṁ), and with ideal-gas carbon no ṁ gives HYDRA's 383 km/s *and* 5.7% left. This
agrees with the user's reading that the lag is rooted in a 1-D rocket model tuned to 3-D
results (ṁ fitted to remaining mass), and points at A22's items 1–2 (P tied to ṁ, real
carbon EOS); the drive is not the cause. Side note: K24's late slope (335–345 km/s) and Ruof
2026's 2DConA 355 ± 5 km/s (N240304) are near stage 0's 340, so the 383 km/s target itself
is model-dependent.

## Hard-sphere direct-drive illumination (2026-10-27) *derived*

Polar direct drive, Stage 0 (ROADMAP "Polar direct drive"; `reference/litreview/polar_direct_drive.md`
§5). Tests: `test/direct_drive.jl`. The parameterization-independent part came first;
the spherical pointing fields of `GroupPointing` (decision 5(a)) and the LMJ check
(Te14b) followed on 2026-10-26, below.

**What exists and is used.** A bare sphere is a one-component target,
`Target(Component(:capsule, icosphere_mesh(R; n), material; surface = SphereSurface(R)))`.
`DepositionOperator`/`deposit` on it give hard-sphere illumination: straight rays, the
incident `I cos ψ` absorbed where it is first incident, the far side shadowed by the sphere
itself; `ncone = 1` makes the beams collimated (the parallel rays of MN17 and Te14).
Direct-drive pointings are aims (documented in the `Beam` docstring, no new names): a
surface point `s` is `aim = s` (first intersection when `port·(s − c) > 0`), a shift Δ
normal to the beam is `aim = c + Δ u` with `u ⟂ port`, an axis crossing the polar axis
at height z is `aim = (0, 0, z)`.

**New API.** `mode_rms(c, lmax; m = :all | :zero | :nonzero)` (flat real-harmonic
coefficients; σ = √(Σ_{1≤l≤lmax} Σ_m c_lm²)/|c₀₀|, generic in the number type, in
`test/inference.jl` with `@inferred` and JET) and `mode_power(C; m)` (the per-l spectrum
σ_l, total or its m = 0 / m ≠ 0 part; m = :zero is the Legendre spectrum
|a_l/a₀|/√(2l+1)). Both work on the exact-solid-angle facet weights of `ylm_weights`,
which now reuse their work array in `_real_ylm!` and evaluate cos mφ, sin mφ once per m
(the same weights exactly; n = 30 icosphere, lmax = 30: 5.3 s → 2.0 s).

**Checks** (R₀ = 0.5 mm; icospheres n = 20 and 30, edges 0.055 R and 0.037 R; σ_l
absolute in units of the mean):
- *Skupsky–Lee separability* (J. Appl. Phys. 54, 3662, 1983), OMEGA's 60 ports
  (`omega_ports`), SG4 spots (`omega_spot(:SG4)`), aimed at the center: the reference is σ_l
  = |a_l/a₀| G_l/√(2l+1) with a_l from a 1-D Gauss quadrature of one beam and G_l from the
  ports. σ (l ≤ 30): reference 0.0821%; Briareus 0.0792% (n = 20), 0.0807% (n = 30),
  Richardson-extrapolated 0.0819%. max_l |Δσ_l| = 3.5e-5 and 1.6e-5 (ratio 2.2, the h²
  rate). Absorbed fraction equals the enclosed fraction of the spot inside the sphere's disk
  to 1e-4. (An exploration at R = 0.43 mm gave 0.1874% for the reference and 0.1869% at n =
  30.)
- *Facet quadrature at l ≤ 30*: the zonal spectrum of `ylm_weights` (level 2) and the
  Legendre spectrum of `legendre_weights` (level 3) agree to 3.8e-6.
- *Temporal et al. 2014* (arXiv:1406.7762), one ring at 49°/131°, m = 4: the direct
  quadrature reference reproduces the published values with incident I cos ψ fully
  absorbed, the shift normal to the beam axis toward the equator and the a axis
  azimuthal: σ = 11.93% (circular, a = b = R₀) and 5.23% (elliptical, a = 2b, b = R₀)
  centered (published 11.9, 5.2); minima 1.088% at Δ/R₀ = 0.1374 and 1.341% at 0.0743
  (published 1.1% at 13.7%, 1.3% at 7.4%). The conventions of the review's §4.2 are
  therefore pinned. Briareus (8 beams per ring, mirrored, m = 0 part, n = 30) is within
  6e-5 of the reference at all four points. dσ/dΔ by ForwardDiff through `deposit`
  matches a central difference to 1e-8 relative.
- *Murakami & Nishi 2017* §4 (Matter Radiat. Extremes 2, 55), Eqs. 28–37 integrated
  directly (a = 0.6, b = 6, η⊥ = 0.9, each ring's b_l normalized to its mean as in
  Eq. 37): the printed solution j = 46°, 53°, z = −0.65, 0.80, E₂ = 0.40 leaves
  B₂…B₁₀ = 0.048, −0.188, 0.046, 0.009, 0.013, far beyond the rounding of the printed
  digits. Newton from there converges (residual 5e-15) to j = 48.926°, 51.685°,
  z = −0.5756, 0.7687, E₂ = 0.4208; weighting the rings by absorbed instead of
  normalized amplitude gives the same geometry with E₂ = 0.542. So either the printed
  digits are from a different model detail or the published root is imprecise (*L*).
  Briareus is checked per ring at the re-solved geometry against the same quadrature,
  without absorption (η⊥ = 1): max_l≤30 |Δ(b_l/b₀)| = 2.5e-3 / 5.1e-3 (n = 20),
  1.1e-3 / 2.3e-3 (n = 30), 4.8e-5 / 7.4e-5 extrapolated; and with `absorption = 0.9`
  (below): 3.2e-3 / 7.3e-3, 1.4e-3 / 3.3e-3, 1.6e-5 / 3.5e-5. The two rings combined
  as Eq. 37 (E₂ = 0.4208) give Briareus B₂…B₁₀ ≤ 1.8e-4 at n = 30 and ≤ 7.4e-7
  extrapolated: the absorption-weighted root is reproduced. MN17 does not give another
  weighted number (Figs. 5, 10, 12, 13 are figures only).

**Absorption weight η_a(ψ) (2026-10-27) *derived*.** `deposit(op, beams; absorption =
η⊥)` keeps the fraction η_a(ψ) = 1 − (1 − η⊥)^(cos³ψ) of each ray at incidence ψ (each
f-cone direction's; the beam axis's with a VisRad profile), MN17 Eq. 29 (they cite
Kruer's *The Physics of Laser Plasma Interactions* for it): inverse-bremsstrahlung
absorption in a thin corona with an exponential density profile, η⊥ = 1 − exp(−8ν_ei
L/3c) at normal incidence. Confidence: the formula **H** (read in MN17, p. 60, and
again for centered beams in Eq. 40); its use for real coronas **L** (thin corona,
parallel rays, no refraction of the deposition pattern). η⊥ has no default and is the
user's estimate, not a calibrated or facility number (MN17 use 0.9 in §4 and 0.95 in
Figs. 5, 12, 13). Unabsorbed light is dropped. Generic in the number type (beams and
η⊥ may be duals; `test/inference.jl` has `@inferred` and JET cases).
- *Off is unchanged.* With `absorption = nothing` (the default) the expressions are
  the previous ones: f-cone and VisRad-profile depositions on the NIF test hohlraum, their
  ForwardDiff Jacobians in (z, r, power) and OMEGA 60 on a sphere were `==` before and
  after the change. `absorption = 1` is the same model but a differently compiled
  expression: it agrees to ~1 ulp but is not identical.
- *Independent check* (`test/direct_drive.jl`): Skupsky–Lee separability still holds
  for centered beams (ψ equals the polar angle from the beam axis), with aₗ of
  I_L(R₀ sin θ) η_a(cos θ) cos θ by 1-D Gauss quadrature. OMEGA 60, SG4, R₀ = 0.5 mm,
  η⊥ = 0.6: max_l |Δσ_l| = 1.25e-4 (n = 20), 5.6e-5 (n = 30), 4.4e-6 extrapolated
  (O(h²); ~3× the unweighted error, because the facets' flat normals now enter η_a as
  well as cos ψ); absorbed fraction 0.42008 / 0.42018 / 0.4202646 against 0.4202646
  (enclosed fraction × the η_a-weighted share). The weight moves σ_l by up to 7.3e-3.
  dσ/dΔ and dσ/dη⊥ by ForwardDiff match central differences to 1e-6 relative (Te14 ring).

**Accuracy of element metrics.** The error is the facet averaging of the element
representation, second order in the element size; the deposition quadrature does not
set it (operator `level` 1, 2, 3 at n = 20: OMEGA max |Δσ_l| 3.5e-5, 3.4e-5, 3.4e-5). Sharp spots (order 6) need finer meshes
than SG4. For σ (l ≤ 30) to ~2% use n ≥ 30 (18 000 elements; the 961 × 18 000 weight
matrix is 140 MB), or extrapolate n = 20/30. Ba24 evaluates σ on HEALPix nside 64 (49 152
pixels, 0.9°); n = 30 elements are about 2.1°.

**Spherical pointing fields (2026-10-26) *derived*.** `GroupPointing(beams; radius = R,
center, fields = (:theta_s, :dphi_s, :split_theta, :split_phi, …))`, decision 5(a): the
aim is a point on the sphere, `center + R port_vector(θ, φ_port + Δφ)`, with θ mirrored
to 180° − θ_s in the lower hemisphere and Δφ measured from the port of each beam (Ba24
Table I's "polar angle, azimuthal offset" form, read as its caption states; the other
reading is ruled out below, "Ba24 σ₁ comparison and optimization"). The splits add ±split_theta / ±split_phi to each beam of
a quad with the sign of its port's offset from the quad's mean port (Co08 split quads,
Ye21 left/right Δφ, as angles; mapping Co08's µm shifts to angles is *L*). Defocus is
measured from the surface aim, as VisRad and NIF measure the focus offset from the
reference point (section "Pointing and focus conventions", **H**); the PDD papers do not
all state their defocus origin (*L*). Template angles are where each template axis first
meets the sphere. Checks (`test/direct_drive.jl`): hand-built aim sets (NIF quads at the
Ba24 Table I angles, relative shifts, a shifted center, the 192 beams with splits, signs
taken from the beam table) to 1e-13 mm; the Te14 ring pointed by θ_s = 49° + asin(Δ/R₀)
(the same axis as the shift Δ) gives the shifted-axis σ to 1e-10 relative and the
quadrature reference to 1e-4; dσ/dθ_s and dσ/dΔφ_s by ForwardDiff match central
differences to 1e-6 relative. The existing (non-spherical) paths are unchanged:
`apply`, `bounds`, `initial_parameters`, their Jacobians and the NIF-hohlraum depositions
and deposition Jacobians of five parameterizations were `isequal` before and after.

**LMJ, Temporal et al. 2014 (HPLSE 2, e8) (2026-10-26) *derived*.** Inputs (paper §2, §4,
**H**): option A, the 20 quads of the 49°/131° rings (10 each, the hemispheres rotated
18°); quads as single beams with I = exp(−[(x/Δa)² + (y/Δb)²]^(m/2)), Δb = 320 µm in the
meridian of the polar and beam axes, Δa = 2Δb, m = 4; static sphere of radius r₀; σ₀ the
rms/mean of I(θ, φ) over the sphere (their equation in §4), no beam errors, no PDD
shift. The quad azimuths are not in the text (Fig. 1 only); uniform 36° spacing is
assumed (*L*). Published: σ₀ = 4.6% is the minimum of the intrinsic curve (Fig. 5; its
minimum lies near r₀ = 350 µm), while r₀ = 320 µm is the optimum of the error-averaged σ̄
= 6.2%, read from Fig. 5; the review's §4.2 item 13 said "4.6% at r₀ = 320 µm", which
is wrong. With Te14's pinned convention (I cos ψ fully absorbed, collimated beams), a
direct 2-D quadrature of all 20 beams gives σ₀ = 5.231% at r₀ = 320 µm and a minimum of
4.589% at r₀ = 350.4 µm (rounds to the published 4.6%; converged to 1e-6 absolute);
with no hemisphere rotation the minimum is 4.64% (also 4.6). Circular Δ = 450 µm (same
1/e area) at r₀ = 320 µm: 10.19%, 1.95× the elliptical value ("almost double" in the
text). Briareus (n = 30 sphere, σ to l = 30, ncone = 1): 4.587% at the minimum and
5.232% at r₀ = 320 µm (within 1e-4 of the reference; the m ≠ 0 part is 0.15%). A test
(`test/direct_drive.jl`). Confidence: inputs **H**, the minimum location from a figure
reading (**M**), azimuths *L* (the value is insensitive to them at the published digit).
The PDD numbers (σ̄ = 3.6% at r₀ = 270 µm, δ/r₀ = 15%) average Monte Carlo beam errors
and are not deterministic; not used. LMJ port geometry was not added to the package: the
paper gives ring angles (33.2°, 49°, 59.5° and supplements, 10/10/2 quads per ring) but
not the azimuths.

**Ba24 σ₁ comparison and optimization (2026-10-26) *derived*.** Decisions 1(b), 2 and 6
of the review: a documented comparison, never a fit target. Tutorial
`docs/literate/polar_direct_drive.jl` (63 s as a script on 8 threads); optimization
`examples/ba24_pdd_optimization.jl`; scratch scripts in the session scratchpad `ba24/`.
- *Inputs* (Barlow et al., PRL 133, 175101, 2024, accepted manuscript OSTI 2475934 and
  arXiv:2311.18342v2, **H**): Table I, the Optimized Config (OC), per cone 23.5°/30°/
  44.5°/50°: polar angle 6.36/14.95/52.48/85.62°, azimuthal offset 52.02/160.75/
  −12.14/6.90°, defocus 8.22/9.07/4.67/6.83 mm, power 66.26/66.26/95.05/82.88% of
  4 TW/quad (the post-optimization rescaling to the requested pulse is taken as common to
  all cones, which leaves σ₁ unchanged, *L*); σ₁ = 6.03% (Table II), the rms/mean of the
  vacuum illumination of the 1.1 mm sphere over l ≤ 30; appendix: the beams of a quad
  overlap exactly at the aim point on the sphere, each beam a plane wave whose
  super-Gaussian 1/e radii grow as w₀(1 + κd), κ = 0.30/0.40/0.80/0.95 per 50 mm (fitted
  to pF3D). Neither w₀ nor the profile shapes (Ifriit's internal "NIF-inner-23" …
  "NIF-outer-50", named in Polaris's deck generator, **H**) nor the intensity convention
  is published.
- *Briareus setup.* `GroupPointing(nif_quads(); fields = (:theta_s, :dphi_s, :defocus,
  :power), radius = 1.1)` at Table I gives the 48 quad aims; each of the 192 `nif_beams`
  takes its quad's aim, defocus and a quarter of its power (Ba24's overlap; the quad's
  port is the mean of its beams', as in Polaris). Two beam models: *own* (LMF CPP fits of
  `nif_spot`, f/20.7 f-cone defocus from the surface aim, ncone 13) and *emulated* (the
  same fits scaled by 1 + κd at zero defocus, collimated, ncone = 1; decision 2's
  emulation, no new API). Icospheres n = 20, 30, Richardson-extrapolated.
- *Independent check* (emulated beams): a plane-wave quadrature written for the purpose
  (Gauss–Legendre 256 in μ × 512 in φ, normalized associated-Legendre recurrence; the
  projected σ₁ converged to 1e-7 relative from 128 to 512 nodes, the incident one to
  1e-4). Briareus's σ_l against it: max |Δσ_l| 7.3e-4 (n = 20), 3.3e-4 (n = 30), 4.0e-6
  extrapolated (units of the mean); σ₁ 7.8669% against 7.8674%.

| case (l ≤ 30, % of the mean) | σ₁ | m = 0 | m ≠ 0 | absorbed | largest σ_l |
|---|---|---|---|---|---|
| Ba24 (Ifriit, their spots) | 6.03 | — | — | — | Fig. 2a only |
| own beams (extrapolated) | 6.569 | 5.266 | 3.926 | 0.936 | σ₈ 3.91, σ₆ 3.32, σ₄ 1.98 |
| emulated (extrapolated) | 7.867 | 6.936 | 3.712 | 0.924 | σ₄ 3.95, σ₁₀ 3.72, σ₆ 3.55 |
| emulated, incident I (quadrature) | 14.25 | 13.66 | 4.06 | — | σ₄ 8.71, σ₂ 7.91 |
| own beams, η_a at η⊥ = 0.9 (n = 30) | 14.36 | — | — | 0.620 | — |

- *Which reading of Table I* (review §4 item 3). The caption gives the polar angle from the
  pole and the offset from the port's azimuth (reading A, `:theta_s`/`:dphi_s`). Method text
  and Polaris's `create_run_files_pdd` (**H**, code read): the aim is the port's surface
  point rotated by θ_o away from the port, in the direction φ_o about it, aim = Rz(φ_p)
  Ry(θ_p) Rz(φ_o) Ry(θ_o) R ẑ with (θ_p, φ_p) the quad's mean port, defined in the lower
  hemisphere's frame and carried to the upper by φ_o → 180° − φ_o ("symmetric") or φ_o +
  180° ("anti-symmetric", Polaris's default). Emulated beams, quadrature, projected σ₁: A
  7.87%; about the port, symmetric 144.7%, anti-symmetric 144.7%, and with φ_o counted from
  the toward-equator direction 30.0% (150% with the native Briareus beams, symmetric). About
  the port, the 50° cone's beams aim 85.6° from their ports and graze the sphere. *Settled*:
  Table I lists the published aims in the caption's reading (the about-the-port parameters
  are the search variables of the method; the table gives the resulting aim angles);
  confidence **H** for the factor 4–20 between the readings, *L* for the interpretation of
  the text. Reversing the offset's sense in the lower hemisphere (the ports' rotational
  symmetry instead of Ba24's "same direction in each hemisphere") gives 8.16% (m = 0
  unchanged, m ≠ 0 3.71 → 4.29%), so σ₁ does not decide it; the text does.
- *Spot sensitivity* (emulated, quadrature): all spots × 0.8/0.9/1.0/1.1/1.2/1.3/1.4 →
  σ₁ 15.98/8.78/7.87/9.26/10.50/11.25/11.65% (the LMF size sits near the minimum, so
  Ba24's w₀ are probably of that size, *L*); a common order 2/3/4/5/6 with the same 1/e
  radii → 6.15/6.01/7.20/8.64/9.93%. The 0.5- and 1.8-point differences from 6.03% are
  within this spread. Causes, in order (*L*): the best-focus spot shapes (LMF fits of
  orders 4.5–8 against Ifriit's internal models); the defocus convention (the f-cone blur
  smooths the discrete outer spots more than the linear growth, σ₁ 6.57 against 7.87%
  with the same w₀); CBET and plasma are absent from σ₁ in both. Ba24's mode spectrum is
  a figure (Fig. 2a), not compared.
- *Optimization of σ₁ only* (decision 6; own beams, n = 20 sphere, Ba24's 16 parameters
  and bounds θ_s 0–90°, Δφ_s ±180°, defocus 0–10 mm, power 2–4 TW/quad; 8 threads pinned
  to one CCX, shared with another job, so times are indicative). Evaluation 0.33 s, a
  ForwardDiff gradient in 16 parameters 1.26 s (≈ 4 evaluations; against central
  differences 8e-9 relative). The bounds do not enforce Ba24's "visible half of the
  target" rule; an aim behind the limb is still a well-defined axis.

| run | σ₁ n = 20 | σ₁ n = 30 | evaluations | wall | best so far |
|---|---|---|---|---|---|
| OC | 6.492% | 6.534% | — | — | — |
| L-BFGS-B from the OC (200-iteration limit) | 2.112% | 2.131% | 229 | 220 s | 4.95% at 10, 2.94 at 25, 2.61 at 50, 2.19 at 100, 2.12 at 150 |
| L-BFGS-B, 4 uniform random starts × 75 | 2.416% | 2.434% | 300 | 279 s | per start 5.25, 3.36, 2.42, 3.74% |
| Nelder–Mead from the OC (unit box) | 2.074% | 2.092% | 1 201 | 290 s | 5.19% at 200, 2.67 at 400, 2.11 at 700, 2.09 at 1 000 |
| particle swarm (24 particles, unit box) | 2.485% | 2.504% | 1 208 | 340 s | 4.81% at 200, 4.17 at 400, 2.82 at 700, 2.56 at 1 000 |

  The gradient needs ~5× fewer evaluations than Nelder–Mead for the same σ₁ but, at four
  evaluations per gradient in 16 parameters, only ~1.7× less wall time; it gains with
  more parameters (per-quad or split pointings). Different starts end at different σ₁
  (2.4–5.3% in 75 evaluations), so σ₁ has several optima; the L-BFGS-B optimum puts
  both inner cones near the pole (θ_s ≈ 0–1°) and the outer ones at 54° and 83°, with
  defocus and power at bounds. Ba24 report > 10 000 evaluations to settle σ₁ by itself with
  their GA and coordinate descent (*H* for the number; their evaluations are Ifriit runs
  with a plasma). These optima are CBET-free vacuum illuminations, not designs.

**Not done (Stage 0 remainder).** `absorption` through `Scenario`/`evaluate` (it is a
`deposit` keyword only).

## Corona model prototype: impact-parameter kernels vs 3-D ray tracing (2026-10-27) *measured*

Decision 3 of `reference/litreview/polar_direct_drive.md` (user, 2026-10-27): build (b),
prototype (c), compare on evidence. Research code in `dev/corona/` (not package API):
`Corona.jl` (both methods), `validate.jl` (closed forms), `compare.jl` (the OMEGA case,
by stages), `table.jl` (the tables below). Both use cold-plasma geometric optics,
n² = 1 − n_e/n_c, with the ray parameter τ (ds = n dτ): dx/dτ = v, dv/dτ = −∇(n_e/n_c)/2,
and inverse bremsstrahlung d ln P/dτ = −(ν_ei/c)(n_e/n_c) (Kruer, *The Physics of Laser
Plasma Interactions*, ch. 5: κ = ν_ei n_e/(n_c c n)), ν_ei the local collision frequency.
Absorbed power is projected radially onto the sphere (absorbed power per steradian,
decision 4).
- **(b)** For a spherical profile each ray stays in its plane and is fixed by its impact
  parameter b. Per b: the turning radius (r_t n(r_t) = b, bisection), and on N panels of
  r = r_t + (R_out − r_t)t² (Gauss, 4 points each; the t² substitution removes the
  1/√ singularity) the swept angle ∫ b dr/(r√(r²n² − b²)) and the depth
  ∫ g r dr/√(r²n² − b²), g = (ν_ei/c)(n_e/n_c). No ray integration. The deposit nodes
  (angle from the incoming direction, absorbed fraction) are tabulated on a uniform b
  grid, and their real-harmonic coefficients κ_lm(b) (m ≥ 0) in the frame z′ = −d. The
  operator takes its rays from a Briareus `DepositionOperator` built on an icosphere at
  R_out: each (point, live f-cone direction) is a ray of power W·I·|n·d|/κ/ncone (the
  expression of the operator, recomputed per pointing, ForwardDiff-generic); a ray at
  impact azimuth ψ about its sub-beam's axis adds κ_lm(b)·(cos mψ, sin mψ); each
  sub-beam's sum is rotated to the lab frame by fixed blocks (exact to l = 30 by a
  product quadrature). The result is exact harmonics, no angular bins.
- **(c)** RK4 in τ on (x, v, ∫g dτ) through any n_e(x) (radial, planar, P2-perturbed),
  rays on a square lattice of spacing h_lat in the plane normal to each beam through
  the target center, along Briareus's `cone_directions` (the same f-cone quadrature as
  the operator), weighted by the best-focus profile at their focal-plane crossing times
  h_lat²/ncone. The lattice is fixed by the port, so only the weights depend on the
  aims and ForwardDiff goes through them (paths in Float64). Absorption per step is
  deposited at the step midpoint, cloud-in-cell onto a 720 × 1440 (θ, φ) grid, then
  projected on Y_lm by rows (FFT-free DFT in φ, exact P̄_lm at bin centers). Nearest-bin
  assignment instead of cloud-in-cell left a first-order bias of up to 4e-4 of the mean
  in σ_l (l ≤ 12) on a coarse run, larger than σ₁₀ itself; with cloud-in-cell the grid
  against exact per-deposit Y_lm is 6e-3 relative at l = 10, 1.5e-3 at l = 30.

**Validation** (`validate.jl`; every reference independent of both methods):
- *V1, planar linear ramp* n_e/n_c = x/L (L = 0.1 mm, ν_c/c = 5 mm⁻¹), (c): absorption
  against 1 − exp(−(32/15)(ν_c L/c)cos⁵θ) (Kruer) within 4.5e-8 relative at hτ = 1e-3
  for θ = 0…75° (2.8e-7 at hτ = 4e-3); re-exit at 2L sin 2θ to 3e-6 mm. Starting a ray
  exactly on the ramp's kink (∇n_e = 0 in the first RK stage) makes the scheme first
  order (2e-3): rays must start just inside.
- *V2, planar exponential* n_e/n_c = e^{x/L}: absorption against
  1 − exp(−(8/3)(ν_c L/c)cos³θ) (Kruer; MN17 Eq. 29, the Stage-0 η_a weight) within
  6.7e-10 for θ = 0…80°; turning depth L ln cos²θ to 1e-6 mm.
- *V3, spherical n_e/n_c = a/r²* (derived closed forms: r_t = √(a + b²); swept angle
  inside radius R, 2(b/r_t)arccos(r_t/R); depth 2(ν_c/c)a²[√(R² − r_t²)/(2r_t²R²) +
  arcsec(R/r_t)/(2r_t³)]): (b) r_t to 2e-16 mm, angle to 4e-12 rad, depth to 3e-12;
  (c) angle to 5e-9 rad, depth to 1.5e-5 (the step that crosses the truncation radius).
- *V4, spherical exponential*, (b) against (c) ray by ray (r_c = 0.4 mm, L/r_c = 0.2,
  b/r_c = 0…2): absorbed fraction within 6.7e-7, exit angle within 4.2e-7 rad.
- *V5*: kernel-table interpolation in b, max error of the absorbed fraction 6.7e-5,
  4.1e-6, 2.6e-7 for Nb = 256, 1024, 4096 (second order); panel count N = 50…400 all
  converged to 3e-11; a tabulated profile (cubic Hermite in ln n_e and ln g, a profile
  with an outward-growing scale length and varying T_e) at Δr = 20, 10, 5, 2.5 µm moves
  the absorbed fraction by 7e-5, 1.5e-5, 2e-6, 4e-7.
- *Operator checks* (`compare.jl check`): the rotation blocks reproduce Y_lm(Qx) to
  1.2e-13 for l ≤ 30; the harmonic operator against exact per-ray harmonics on a hard
  sphere, σ (l ≤ 30) 19.195913% vs 19.195938% (max |Δc_lm|/c₀₀ 3.7e-4, from linear
  interpolation in b across the sphere's limb). dσ/dΔr through (b) by ForwardDiff
  matches central differences to 2.7e-9 relative.

**The Stage-0 η_a weight against a spherical corona** (V4 table: absorbed fraction of a
ray of impact parameter b in an isothermal exponential corona with (8/3)(ν_c L/c) = 2.4,
η⊥ = 0.909, against η_a(ψ) with sin ψ = b/r_c). Normal incidence agrees exactly at every
L/r_c; oblique rays converge to η_a as L/r_c → 0 (b/r_c = 0.9: 0.1868 at L/r_c =
0.0025 vs η_a 0.1803) but at L/r_c = 0.1 and 0.2 absorb 0.332 and 0.406 (η_a 0.180), at
b/r_c = 0.95 0.261 and 0.349 (η_a 0.070), and rays that miss the critical sphere
(b > r_c) still absorb (0.12 at b = 1.2 r_c, L/r_c = 0.2), which a hard sphere cannot
do. So η_a is consistent with the corona only for thin coronas (at b/r_c = 0.9: 4%
at L/r_c = 0.0025, 13% at 0.01); for PDD, where the equatorial drive comes from exactly these grazing
rays, it is not (σ 5.7% against 10.4% on the OMEGA case below; *derived*).

**OMEGA 40-beam PDD comparison** (`compare.jl`; public inputs, review §4.1): rings
21.4°/42.0°/58.9° of `omega_beams` shifted Δr = 90/150/150 µm normal to the beam toward
the PDD equator (Mn14 Table 1, Ma16's spherical-shell pointing), SG4 r₀ = 308 µm,
n = 3.66 (Ma16), f/6.7, ncone = 7, equal powers. Corona (an engineering snapshot, not
from a simulation): isothermal exponential, r_c = 0.40 mm, L = 0.08 mm, ν_c/c =
11.25 mm⁻¹ (η⊥ = 0.909 planar; about T_e 3.6 keV, Z 3.5, lnΛ 8 at 351 nm), truncated
at r_c + 14L = 1.52 mm. σ_l in % of the mean; the last column is max over 1 ≤ l ≤ 30 of
|σ_l − σ_l(ref)| with ref = (c) at h_lat = 0.01 mm. Odd l vanish by the up–down symmetry
((b) exactly; (c) ≤ 1e-6).

| model | absorbed | σ (l≤30) | σ₂ | σ₄ | σ₆ | σ₈ | σ₁₀ | σ₁₂ | σ₂₀ | σ₃₀ | max Δσ_l |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hard sphere at r_c, all absorbed | 0.8914 | 16.17 | 15.88 | 0.323 | 2.342 | 0.984 | 0.212 | 0.288 | 0.291 | 0.800 | 5.8 |
| hard sphere at r_c, η_a (Stage 0) | 0.5815 | 5.66 | 2.096 | 4.389 | 1.964 | 2.104 | 0.151 | 0.098 | 0.054 | 0.155 | 8.0 |
| (b), icosphere n = 20 (4.4 M rays) | 0.64592 | 10.4072 | 10.0753 | 1.5143 | 1.7733 | 1.1606 | 0.0860 | 0.0597 | 0.0026 | 0.0239 | 2.3e-2 |
| (b), n = 30, Nb = 4096 (9.9 M) | 0.64591 | 10.4087 | 10.0771 | 1.5136 | 1.7726 | 1.1614 | 0.0861 | 0.0591 | 0.0007 | 0.0033 | 2.2e-3 |
| (b), n = 40 (17.6 M) | 0.64591 | 10.4092 | 10.0777 | 1.5133 | 1.7724 | 1.1616 | 0.0862 | 0.0588 | 0.0006 | 0.0003 | 9.9e-4 |
| (c), h_lat = 0.08 mm (98 k rays) | 0.64591 | 29.13 | 10.0783 | 1.5125 | 1.7726 | 1.1623 | 0.0862 | 0.0566 | 1.948 | 21.8 | 22 |
| (c), h_lat = 0.04 (391 k) | 0.64591 | 10.4097 | 10.0783 | 1.5124 | 1.7723 | 1.1623 | 0.0857 | 0.0585 | 0.0010 | 0.0018 | 6.6e-4 |
| (c), h_lat = 0.02 (1.57 M) | 0.64591 | 10.4098 | 10.0783 | 1.5124 | 1.7723 | 1.1623 | 0.0858 | 0.0585 | 0.0009 | 0.0012 | 1.2e-4 |
| (c), h_lat = 0.01 (6.27 M) | 0.64591 | 10.4098 | 10.0783 | 1.5124 | 1.7723 | 1.1623 | 0.0858 | 0.0585 | 0.0009 | 0.0011 | ref |

(c) at h_lat = 0.04 with hτ = 2e-3, 4e-3, 8e-3 (L/hτ = 40, 20, 10): max Δσ_l 5.9e-4,
6.6e-4, 1.8e-3. The two methods agree on the absorbed fraction to 1e-5 and on every
σ_l to 1e-5 of the mean once each is converged; the residual is (b)'s sampling, which
is the `DepositionOperator`'s second-order quadrature over the R_out sphere: σ₂ moves
10.0753 → 10.0771 → 10.0777 for n = 20, 30, 40 and Richardson-extrapolates (h²) to
10.0785 against (c)'s 10.0783. (c)'s square lattice converges faster (the midpoint rule
on a smooth, compactly supported profile): 391 k rays (~1 400 per beam and f-cone
direction, spacing 0.13 r₀) give σ_l (l ≤ 30) to 7e-6 of the mean; 98 k rays alias at
high l (σ₃₀ 22% of the mean).

| work (8 threads, Ryzen CCX) | (b) n = 20 | (b) n = 40 | (c) h_lat = 0.04 |
|---|---|---|---|
| build | operator 0.6 s + kernels 0.6 s + map 8.6 s | 3.2 + 0.6 + 8.4 s | none |
| memory | 43 MB operator + 234 MB map | 171 + 638 MB | 8 MB grid per thread |
| per pointing change | 0.16 s | 0.28 s | 18 s |
| ForwardDiff, 3 ring shifts | 0.56 s | 1.5 s | 22.6 s |
| ForwardDiff, 80 aim coordinates | 11 s | 32 s | not run; ~7 chunks × ~22 s ≈ 150 s (*L*, estimate) |

About 3 s of the map is the 280 sub-beam rotation blocks (0.09 s each, fixed per port;
a quadrature, where a Wigner recursion would be cheaper), the rest the per-ray loop. (c)'s gradient costs little more than its
value because its paths do not depend on the aims; ForwardDiff's chunking retraces
once per chunk, hence the estimate for 80 parameters. The gradients agree:
dσ/dΔr = (1.04755, 0.76448, 0.58639) for (b) n = 40 and (1.04751, 0.76453, 0.58645)
for (c), 1e-4 relative.

**What (c) gives that (b) cannot: a P2-perturbed corona.** n_e/n_c =
exp(−(r − r_c(1 + ε P₂(cos θ)))/L) (critical surface +ε r_c at the poles, −ε r_c/2 at
the equator), (c) at h_lat = 0.04. (b) returns the ε = 0 result for every ε.

| ε | absorbed | Δ(c₂₀/c₀₀) | Δ(c₄₀/c₀₀) | σ₂ | σ₁₀ |
|---|---|---|---|---|---|
| 0 | 0.64591 | — | — | 10.078% | 0.086% |
| +0.02 (pole +8 µm) | 0.64117 | +0.00137 | +0.00050 | 9.941% | 0.068% |
| +0.04 (pole +16 µm) | 0.63627 | +0.00285 | +0.00099 | 9.793% | 0.062% |
| −0.04 | 0.65482 | −0.00249 | −0.00092 | 10.328% | 0.140% |

Linear in ε to 4% at ε = 0.02 (±0.04 differ by 13% in Δc₂₀): a P2 of the critical radius of ε moves σ₂ by about
−7 ε (percentage points) and the absorption by −0.24 ε, i.e. a 16 µm polar bulge on a
400 µm critical radius changes σ₂ by 0.29 points (3% of its value) and the absorption
by 1 point. (b) misses that; but it is a first-order, smooth response that a
perturbative correction measured with (c) would capture.

**Verdict.**
1. *Stage 1 operator: (b).* On a spherical corona both methods give the same answer
   (absorption to 1e-5, σ_l to 1e-5 of the mean), and both are validated against the
   closed forms; (b) evaluates a pointing change 60–110× faster than converged (c),
   gives exact harmonics, and its gradient in 80 aims costs 11 s against an estimated
   150 s. It keeps the operator structure, so it slots into the optimizer like
   `deposit`. One change from the prototype is to take (b)'s rays from a transverse lattice
   per sub-beam, as (c) does, in place of the points of the `DepositionOperator` on the
   R_out sphere. The lattice converged at 391 k rays where the sphere sampling needs
   17.6 M (n = 40) for the same l ≤ 30 accuracy, so the map and the evaluation would be
   ~10–40× smaller again (built 2026-10-27: 7× faster per evaluation, 250× less memory;
   next section).
2. *Hybrid: yes, (b) as the operator and (c) as the check.* (c) is ~150 lines of
   independent code (shared only with the beam model), already validated, and is the
   only way to measure what spherical symmetry drops (the P2 corona above: 3% in σ₂
   at a 4% critical-surface P2). Keep it as a test reference (converged (b) against
   (c) at h_lat = 0.04 on a spherical corona, a few seconds on a reduced beam set) and
   as an offline study tool for non-spherical snapshots.
3. *The Stage-0 η_a weight* agrees with the corona at normal incidence and in the
   thin-corona limit only; at L/r_c ~ 0.1–0.2 it underestimates grazing absorption by
   2–5× and ignores rays outside r_c. On the OMEGA case it gives σ 5.7% where the
   corona gives 10.4%. Its docstring's "thin corona" caveat should name this limit
   (L/r_c of a few 1e-3 for a few % at b/r_c = 0.9).
4. *Not settled here*: kernels per snapshot for a time-dependent corona (stage 2; a
   table costs ~1 s, so rebuilding per snapshot is cheap); non-isothermal ν_ei(r) is
   already supported by the table form (g(r) given separately).

*Superseded (2026-10-27)*: (b) is now package code (next section); `dev/corona/` keeps
only the 3-D tracer (c) (`Corona.jl`, `validate.jl` V1–V2) and `p2.jl`, the P2 study
above against the package operator. `compare.jl` and `table.jl` are removed; the
numbers above are their record.

## Corona absorption operator (2026-10-27) *measured*

Polar direct drive, Stage 1 (ROADMAP "Polar direct drive"; review decision 3 (b), with the
prototype's verdict above): `CoronaProfile`, `CoronaOperator` and
`deposit(op::CoronaOperator, beams) -> (; ylm, per_beam, delivered)` in `src/corona.jl`;
tests `test/corona.jl` with the 3-D tracer as a test helper (`test/corona_raytrace.jl`,
~150 lines, model (c) restricted to a spherical exponential corona). The output is the
absorbed power per steradian about the center as flat real-harmonic coefficients
(decision 4), for `mode_rms` / `mode_power`. CBET-free (decision 8).

**Model and design.** As (b) above: cold-plasma geometric optics, inverse
bremsstrahlung, radial projection; per impact parameter the turning radius by
bisection and N = `panels` (200) panels per pass of r = r_t + (R − r_t)t², Gauss 4 points
each, a deposit at each panel's mid angle. Changes from the prototype:
- *Rays from a transverse lattice per sub-beam* (the verdict's item 1): for each beam and
  f-cone direction d a square lattice of spacing h (default 0.1 × the smallest spot
  semi-axis) in the plane normal to d through the center, rays with b < R (the profile's
  outer radius), weight = best-focus profile at the focal-plane crossing × h²/(k·d)/ncone
  (the collimated sub-beam's power per lattice cell; the `DepositionOperator`'s
  expression with n = d). In the frame of each sub-beam (pole −d, x along the lattice's
  first axis) every sub-beam has the same lattice, so the impact parameters are the
  distinct values h√(i² + j²): 745 for 7 653 rays at h = 0.0308 mm, R = 1.52 mm. The
  kernels are computed at exactly those b, so the prototype's linear interpolation in b
  (Nb = 2048) is gone.
- *Rotation by a fixed matrix*: the frame is R_z(α)R_y(β) (Euler angles of −d) and its
  harmonic rotation Z(α) J Z(β) Jᵀ with J the blocks of R_x(−90°) (one 39 711-entry set
  for l ≤ 30, by product quadrature) and Z the cos/sin mixing per m. The prototype stored
  a block set per sub-beam (0.09 s and 318 kB each, so 790 MB for NIF's 192 × 13). The
  factorization reproduces the quadrature blocks of R_z(α)R_y(β) to 8.8e-14 (l ≤ 30).
- *Exact screen*: the lattice rows and columns that can reach the spot are bounded by
  Cauchy–Schwarz (|h i − c₁| ≤ √qcut ‖(a e₁·e₁′, b e₂·e₁′)‖, plus one index); rays
  outside add exact zeros. Dual beams are screened on values, as in `deposit`; the dual
  path's values equal the Float64 path's exactly (tested).
- *Profiles*: `CoronaProfile(r, ne_over_nc, Te; Z, lnΛ, wavelength)` tabulates a snapshot (the
  hook for the CompactLES example), cubic Hermite in ln(n_e/n_c) and ln of the rate with
  three-point slopes on a non-uniform grid, so exponential profiles are exact. Rate
  (ν_ei/c)(n_e/n_c) with Kruer's ν_ei = 4√(2π)Z e⁴ n_e lnΛ/(3 m_e^{1/2} T_e^{3/2}) from
  CODATA 2018 constants: 2.906e-6 Z n_e lnΛ/T_e^{3/2} s⁻¹ (n_e cm⁻³, T_e eV; the NRL
  Formulary's 2.91e-6), n_c = π m_e c²/(e²λ²) = 1.1149e21/λ² cm⁻³ (**H**, textbook
  formulas). `ne_over_nc[1] ≥ 1` is required (every ray turns inside the table; no hard
  surface). `CoronaProfile(; critical_radius, scale_length, absorption)` is the analytic
  isothermal exponential, its rate set by the planar normal-incidence absorption η⊥ of
  `deposit(...; absorption)`, truncated at r_c + 14L by default.

**Checks** (`test/corona.jl`; references independent of the kernels, lattice and
rotations):
- *Closed forms*: the radial ray of the truncated exponential corona,
  1 − exp(−2(ν_c L/c)(2√w − (2/3)w^{3/2})), w = 1 − e^{−(R−r_c)/L}, to 4e-14 at N =
  50–400; in n_e/n_c = a/r² (as functions) the swept angle to 2.5e-12 rad and the depth to
  1.8e-12 relative at b = 0…1.5 (the prototype's V3); absorbed + escaped = 1 to 1e-16.
  A tabulated isothermal snapshot on a random non-uniform grid (T_e 3.6 keV, Z 3.5, lnΛ 8,
  351 nm) reproduces the analytic corona of the same ν_c within 2e-3 (the NRL constants'
  rounding in the reference).
- *Against 3-D ray tracing* (three OMEGA beams, one per PDD ring with the shifts below,
  ncone = 3, corona as below): absorbed fraction 0.6502373 against 0.6502378 (5e-7);
  max_lm |Δc_lm|/c₀₀ (l ≤ 30) 2.9e-5 with the tracer at h_lat = 0.04, hτ = 0.008 and a
  720 × 1440 grid, 1.3e-5 at h_lat = 0.02, 1.6e-4 with a 360 × 720 grid: the residual is
  the tracer's angular binning. Test tolerances 5e-6 and 1e-4.
- *Truncation*: the profile's outer radius is a density step; (b) refracts at it by
  Bouguer's invariant, the tracer starts rays at |v| = 1 there. At R = r_c + 8L
  (n_e = 3.4e-4 n_c) the absorbed fractions then differ by 3e-4 relative; at 14L by
  5e-7. Compare the two only where n_e(R) is negligible.
- *Thin-corona limit* (one centered collimated SG4 beam, a = 0.25 mm, r_c = 0.4 mm,
  η⊥ = 0.6, against a hard sphere at r_c absorbing η_a(ψ) I cos ψ by 1-D quadrature):
  absorbed-fraction error 9.2e-4, 4.8e-4, 2.5e-4, 1.3e-4, 6.4e-5 and max_l |Δσ_l| 9.5e-2,
  5.1e-2, 2.7e-2, 1.4e-2, 6.9e-3 at L/r_c = 0.04, 0.02, 0.01, 0.005, 0.0025: first order
  in L/r_c (ratio 1.9–2.0), so η_a is the operator's thin limit (*derived*). The σ_l error
  at small L/r_c is the grazing rays (V4 above).
- *Derivatives*: dσ/dΔr (3 ring shifts) by ForwardDiff against central differences
  1.1e-8 relative; the absorbed power's 1.9e-10.

**OMEGA 40-beam PDD case** (the prototype's case and corona, r_c = 0.40, L = 0.08 mm,
R = r_c + 14L; ncone = 7): absorbed 0.645914; σ (l ≤ 30) 10.4102%, σ₂ 10.0788,
σ₄ 1.5127, σ₆ 1.7720, σ₈ 1.1621, σ₁₀ 0.0862, σ₁₂ 0.0584, σ₂₀ 0.0006% (identical to four
decimals at h = 0.04 and 0.02 mm). Against (c) at h_lat = 0.01 (the table above) the
listed σ_l differ by at most 1.0e-3% of the mean, as the prototype's (b) at n = 40
(9.9e-4); the convergence of the operator (max_l |Δσ_l| against h = 0.015 mm and
400 panels): 5.8e-5% at the defaults, 7.6e-5% at h = 0.04, 2.9e-4% at 100 panels,
1.2e-3% at 50. So the residual against (c) is (c)'s. dσ/dΔr = (1.04754, 0.764528,
0.586444) against (c)'s (1.04751, 0.76453, 0.58645).

| work (8 threads, CCX shared with another job) | operator, defaults | prototype (b), n = 40 | (c), h_lat = 0.04 |
|---|---|---|---|
| build | 0.8 s | 12.2 s | none |
| memory | 3.2 MB (+ 3.2 MB scratch per evaluation) | 809 MB | 8 MB per thread |
| per pointing change, ncone = 7 / 13 | 0.040 s / 0.075 s | 0.28 s | 18 s |
| ForwardDiff, 3 ring shifts | 0.065 s | 1.5 s | 22.6 s |
| ForwardDiff, 80 aim coordinates | 2.1 s | 32 s | ~150 s (estimate) |

NIF scale (cost only: 192 beams aimed at the center, SG4 a = 0.6 mm, f/22, r_c = 1.1,
L = 0.15 mm): 8 945 rays and 859 impact parameters per direction, build 0.6–0.8 s,
4 MB, an evaluation 0.14 s at ncone = 7 and 0.23 s at 13. `test/corona.jl` takes 18 s on
one thread, compilation included.

**Limits** (in the docstrings): spherical coronas only (a 4% P2 of r_c moves σ₂ by 3%,
above); bare spheres only (no other parts, which would need their own shadowing);
CBET-free. **Not done**: absorbed power per radius; the P_abl and Ye21 drive proxies
(decision 4's optional transforms); the CompactLES snapshot example; `Scenario` /
`evaluate` integration (objectives call `deposit` directly); a tutorial. (The first two
done 2026-11-02, next section.)

## Corona deposition per radius and drive proxies (2026-11-02) *measured*

Polar direct drive, the rest of Stage 1 (review §5; decision 4): `CoronaOperator(...;
shells)`, `deposit(op, beams).shell_ylm` / `.radial`, `ablation_pressure_proxy`,
`shell_displacement_proxy` (`src/corona.jl`, tests in `test/corona.jl`).

**Deposition per radius.** `shells` (edges covering the profile) become panel boundaries
of every ray (the t-grid of r = r_t + (R − r_t)t² takes t = √((e − r_t)/(R − r_t)) for
each edge inside the ray), so each panel lies in one shell and its deposit and mid angle
are as before; per shell the operator keeps a zonal kernel set and the absorbed fraction
per impact-parameter class, and `deposit` contracts and rotates each shell as it does
the totals. Without `shells` the operator and `deposit` are identical to the previous
code (kernels, `ylm`, `per_beam`, Float64 and duals; checked against the previous
`corona.jl` loaded beside the new one). With shells the totals come from the refined
panels: absorbed power equal within round-off, harmonics within 0.7–12 × 10⁻⁷ of c₀₀ for
14–448 shells on the OMEGA case below (the angular placement at panel mid angles).
Checks (references independent of the operator):
- radial ray of the truncated exponential corona: one-pass depth outside r,
  τ(r) = (ν_c L/c)[H(x(r)) − H(x(R))], x = n_e/n_c, H(x) = −2√(1−x) + (2/3)(1−x)^{3/2};
  shell s absorbs e^{−τ(b)} − e^{−τ(a)} + e^{−(2T−τ(a))} − e^{−(2T−τ(b))} (T = τ(r_c)):
  agreement 1e-14–6e-14 at N = 50–400 panels for shells from L/4 to 10L wide (*derived*);
- oblique rays in n_e/n_c = a/r² (a = 0.04, rate 2 n_e², R = 2, b = 0–1.5): one-pass
  depth inside ρ, νa²[√(ρ² − r_t²)/(2r_t²ρ²) + arcsec(ρ/r_t)/(2r_t³)]: 3.9e-14;
- partition: Σ radial = Σ per_beam, Σ_s shell_ylm = ylm (4e-16), √(4π) shell_ylm[1, s] =
  radial[s] (1e-15).

**Ba24 ablation-pressure proxy.** Barlow et al., PRL 133, 175101 (2024), appendix "on
approximating ablation pressure", Eqs. (2a)–(2b) of the accepted manuscript (OSTI
2475934); Eqs. (1a)–(1b) in arXiv:2311.18342v2; read from page images, both versions,
same formulas (**H**): P_abl = 24.7 Mbar (F/10¹⁴ W cm⁻²)^{2/3}, F = Σ_r r² n_r^{2/3} I_r /
(R₃₅₁² n₃₅₁^{2/3}), with I_r the intensity absorbed in a cell of radial extent Δr and solid
angle ΔΩ such that the power is Σ ΔΩ r² I_r (the manuscript's added definition), R₃₅₁ and
n₃₅₁ the 351-nm critical radius and density. The 2/3 power follows Manheimer, Colombant &
Gardner, Phys. Fluids 25, 1644 (1982) (their Ref. 33); the n^{2/3} weight is empirical,
"chosen to match pressure and asymmetries" of several of their 3-D runs (**H** that it
is what they did; *L* as physics elsewhere). Fitness (their Eq. 1, arXiv Eq. 2) uses σ₂,
the rms of P_abl over l ≤ 30, and P̄_abl. Implementation: F(Ω) = r_c⁻² Σ_s
(n_s/n_c)^{2/3} dP_s/dΩ in TW/mm² (= 10¹⁴ W/cm²), n_s at the middle of the shell's part
inside the profile, r_c the outermost n_e = n_c of the profile (n_c of the profile's
wavelength), the 2/3 power applied on a Gauss–Legendre product grid (2L + 2 × 4L + 4,
exact for results to degree 3L + 3) and projected to l ≤ lmax. Both exponents are
keywords (decision 4). Check: an axisymmetric l ≤ 1 deposition in three shells against a
1-D Legendre projection of 24.7 F^{2/3} with hand weights: 4e-15 (*derived*).

**Ye21 shell-displacement proxy.** Yeamans et al., Nucl. Fusion 61, 046031 (2021)
(LLNL-JRNL-814327, OSTI 1772693), Eq. (3), read from the page image (**H**):
d(θ, φ) = d̄(θ) · 0.72 E(θ, φ)/Ē(θ) as typeset, with the text: shell velocity scales as
the 0.72 power of laser energy (their Ref. 47, Tucker, LLE report 2011, from implosion
simulations), close to an experimental 0.71 they cite; Ē(θ) is the φ-average of the
deposited energy per steradian, d̄(θ) the 2-D run's mean distance moved. The typeset
form is ambiguous (the factor sits beside E/Ē, not as an exponent); the function uses
the power d̄ (E/Ē)^{0.72} the text describes (*M* for the reading). The linear reading
d̄ (1 + 0.72(E/Ē − 1)) agrees to first order and makes their "ΔV_m reduced by a factor of
0.72 from the ΔE_m" exact; on the OMEGA case below the power gives ΔV_{m≠0}/ΔE_{m≠0} =
1.06495/1.47559 = 0.7217. Their reported ΔV = 2.22% (m = 0: 1.23% from their 2-D
hydro, m ≠ 0: 1.85%) is not reproduced here (no comparable profile; the review's
benchmark). Check: E = A Y₀₀ + B Y₁₀ + C Y₁₁ with d̄(θ) = 1 + 0.05 cos 2θ against
hand-written harmonics on a 100 × 200 grid: 1.4e-14 (*derived*).

**OMEGA 40-beam PDD case** (the corona-operator case: r_c = 0.40, L = 0.08 mm, η⊥ =
1 − e^{−2.4}, Mn14 ring shifts 0.09/0.15/0.15 mm, ncone = 7; an engineering corona, not a
snapshot): absorbed 0.645914, σ₁ (absorbed, l ≤ 30) 10.410%; P̄_abl and σ₂ against the
shells:

| shells | P̄_abl (Mbar) | σ₂ (%) | `deposit` (8 threads) |
|---|---|---|---|
| none | — | — | 0.039 s |
| 14 uniform (1.0 L) | 100.70 | 5.046 | 0.31 s |
| 56 uniform (0.27 L) | 104.67 | 4.012 | 0.84 s |
| 224 uniform | 104.61 | 3.996 | 3.1 s |
| 448 uniform | 104.60 | 3.994 | 6.4 s |
| 30: L/4 to r_c + 4L, then L | 104.56 | 4.025 | 0.42 s |
| 60: L/8 to r_c + 4L, then L/2 | 104.60 | 3.993 | 0.86 s |

So shells of L/8 near the critical radius give σ₂ within 0.1% and P̄ within 0.01%, L/4
within 1%; the cost is about 14 ms per shell here (each shell's contraction and
rotation), the proxy 2 ms. σ₂ is 0.38 of σ₁ (the 2/3 power and the density weight).
Absolute P̄_abl depends on the corona (an engineering snapshot), not a prediction.

**Use in an optimization.** No `Scenario` field: an objective for `optimize_pointing`
captures the operator, `ev -> mode_rms(ablation_pressure_proxy(op, deposit(op,
ev.beams)), 30)^2`, on a scenario over a bare sphere (whose surface deposition is then
computed and unused). Checked: 40 beams, `:theta_s` per ring, 30 shells, three L-BFGS-B
iterations take σ₂² 1.68e-3 → 3.49e-4 (15 s with compilation). Not adopted as a test.

## Drive asymmetry inferred from post-shot capsule simulations: Haines 2026 (2026-11-18) *published; figure values read by eye*

Haines et al., "Relating metrology to performance in recent high-yield implosions on the
National Ignition Facility," Phys. Plasmas 33, 112712 (2026) (CC BY), post-shot
simulations of 24 layered implosions with xRAGE in 2-D, driven by frequency-dependent
sources from HYDRA runs calibrated to bang time. The paper concerns capsule metrology;
what bears on a view-factor drive model:

- **A transfer curve for one design** (their Fig. 9; a 2.05 MJ HybridE capsule, 85 µm
  ablator, the design of N230729). A flux asymmetry P2/P0 fixed in time is applied to
  the source; the synthetic neutron-image P2/P0 changes little for applied values of
  −0.03 to +0.035 (the hot spot is exploding at bang time on the ignition cliff) and
  steeply outside; the yield is unaffected for −0.01 to +0.015 and, read from the
  figure, falls to half near −0.021 and +0.028 and to under 5% beyond about −0.035 and
  +0.045. A negative applied P2 gives a prolate image. The neutron-image P2 of an
  igniting shot is therefore not a linear measure of the drive's P2: a linear
  flux-to-hot-spot transfer (roadmap item 11) holds for in-flight shapes or below the
  cliff only. *Confidence:* the ranges in the text are the paper's; the half-yield
  points are read from the figure (±0.003).
- **Three nominal repeats** of that design and drive (same figure, vertical lines; read
  ±0.003): the fixed applied P2/P0 that reproduces each shot's observed shape is about
  −0.024 (N230729, 9% prolate, 52% of yield lost to shape in their Table II), +0.003
  (N250622, the shot with two wall windows; no shape loss) and +0.022 (N250930; 30%).
  The spread between nominal repeats, 0.046, is about ten times the capsule P2 that
  50 µm of pointing error per quad produces in a view-factor model (0.4% in the same
  Legendre units on the dropped-bundle tutorial's case), and exceeds the m = 0 signal of the
  windows themselves: shot-to-shot shape at this scale is set by what a view-factor model
  does not carry (energy transfer and plasma conditions), and a difference between a
  windowed and an unwindowed shot cannot be read as the windows' effect. The pre-shot
  uncertainty that the paper states for the applied P2/P0 is 0.005 for an exact repeat plus
  0.005 for a change of pulse shape or wavelength detuning, calibrated on earlier shots, and
  its authors note that these repeats suggest more; an earlier ±3% assumption on the flux
  asymmetry "proved inadequate".
- **Scales for gates.** Mode gates set at a fraction of the pointing-noise floor
  (tenths of a percent of the drive per degree) are several times inside both the band
  with no yield effect and the repeatability above. They state "no change from the
  reference" in place of a physical tolerance; where such a gate costs drive, a gate in
  physical units (0.5% in P2/P0, the modelers' repeat uncertainty) is defensible. A
  reference shot's drive is not thereby a good target: N230729's inferred drive is 2.4%
  waist-heavy in the fixed-in-time equivalent.
- **The windowed shot** (their Sec. VII): two 2-mm-diameter windows of "≈ 3 µm" of gold
  (the Letter of Lester et al. gives 2.5 µm); a gap between window and wall, not
  measurable before the shot, estimated at no more than 400 µm as built and 200 µm
  after closure, treated as a lost-drive asymmetry in the pre-shot uncertainty; the shot
  matched the nominal prediction, which the authors read as the gap having closed.
- **Published histories** (their Fig. 1(a)): total laser power and Dante radiation
  temperature for N230729 through the reheating peak after the laser turns off (about
  55 TW in the picket, about 465 TW at peak; Dante near 150 eV at 3 ns, 275 eV at 6 ns,
  310 eV at 8.3 ns, 360 eV at the reheat peak near 10 ns; read by eye, ±10 eV): a second
  public drive history for the rise-through-the-peak error of item 12, as a check.
- Mode 1 from a diagnostic window is compensated by the same window opposite (their
  reference to Young et al. 2020), as the default hohlraum models it.
- The capsule-quality fit and the fill-tube results are not used (outside the scope). The
  paper states that integrated simulations predict neither the coupling nor the shape,
  and that the flux asymmetry is inferred from in-flight and stagnation images, not
  measured: none of its drive numbers is a measurement of the drive.
