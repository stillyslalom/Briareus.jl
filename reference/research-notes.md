# Research notes

These are literature findings gathered from public sources (October 2026) that set
Briareus models and defaults. Confidence levels:
- **H**: primary text read.
- **M**: abstract or secondary source only.
- **L**: inference.

Numbers marked *derived* are our own computations.

Sections are appended as work is done and are not rewritten afterwards: they are the
provenance of every number they hold. When a later section changes a value, the earlier
one carries a **Superseded** callout naming the replacement section; numbers in a
superseded section are of their day. Cross-references quote section titles (a `##`/`###`
heading or a bold bullet title).

## Contents

Status: *current* unless marked; "in part" means some values were replaced later (the
callout in the section says which).

**Laser–plasma, facility data and file formats**
- Cross-beam energy transfer (CBET); Wavelength detuning (NIF); Backscatter; Unconverted
  light (NIF); NIF quad names and bundles; Pointing and alignment; VisRad; State of the
  art (pointing design and tools); Beam clearance methodology (clearance defaults).
- Later search (2026-10-16): LMF Shot RI Resources and the VisRad manual (NIF beam
  table, CPP spots, TaLIS rules, OMEGA).
- VisRad interchange formats (2026-10-16): SST, summary, power and workspace formats;
  "VisRad workspace files (2026-10-17)"; "VisRad deposition parity: hohlraum runs
  (2026-10-21)" (`VisRadProfile`: growth 0.67, measured plates, the faint wings); "CAD
  import, remeshing and surface fits" (bullets "Flat CAD elements", "Tori", "Cost of
  fitting faces of none of the kinds", "Fitting CAD by default", "Sampled visibility on
  folded elements", "Parts of fitted CAD"); result tables; items to verify at work.

**Radiation, meshes and scale**
- Gold wall loss and thin-foil burnthrough (2026-10-16); Hohlraum mesh work-precision
  (2026-10-16); Compact capsule mesh (2026-10-22); Partial visibility (2026-10-22); Monte
  Carlo radiosity reference for partial visibility (2026-10-22) (*in part*: balancing on
  flat facets); Balancing on the curved path against a continuum Monte Carlo
  (2026-10-22); Default capsule mesh and visibility
  (2026-10-22) *decision*.
- Large radiation models (2026-10-24); Hierarchical view factors (2026-10-24) (*in
  part*: `tol`); The H-matrix's N^1.5 remainder (2026-10-25); H-matrix near field: zero
  blocks and cut rows (2026-10-26; vectorized sector contour; contiguous quadrature, dropped); H-matrix visibility certificates (2026-10-27); Deposition-operator rays
  and sliver view factors (2026-10-25); Compile latency (2026-10-25); Closure on bored,
  windowed and CAD cavities (2026-10-25); Compact capsule rows at capsule_h = h/3
  (2026-10-25).

**Design, errors and shot review**
- Beam-error budgets (2026-10-17); Retuning to a reference drive (2026-10-17); Implosion
  metrics (2026-10-22); Young 2020 validation (2026-10-24); LEH shield (sombrero)
  (2026-10-24); Backlighter foil on the line of sight (2026-10-24); Backlighter
  conversion efficiency and the photometric chain (2026-10-24); Weighted He-α CE
  compilation (2026-10-26); He-α CE: saturating peak (2026-10-26).

**Hohlraum dynamics and CBET calibration** (HYBRID-E)
- HYBRID-E (N210808) hohlraum and LEH size (2026-10-20); Literature review against
  Kritcher et al. 2024 (2026-10-22); Mock HYBRID-E pointing (2026-10-23).
- Moving wall on keyframes (2026-10-23) (*in part*: speed); Closing entrance holes
  (2026-10-23); Inner-beam bubble absorption (2026-10-23) (*in part*: defaults); both
  refitted in Hohlraum dynamics calibrated to published data (2026-10-24); Computed LEH
  closure (2026-10-24).
- Two-zone CBET calibrated to PRE 109 Fig. 9 (2026-10-23) (`n221204_cbet`); Joint
  multi-shot CBET calibration (2026-10-24) (`hybride_cbet`); Joint multi-shot CBET
  calibration: the null (2026-10-25) (`hybride_cbet_iraum`); Inner-cone SBS time
  histories and the late-time refraction (2026-10-25).

**Energetics and the N221204 benchmark** (in order; read the last for the current state)
- N221204 time-dependent benchmark (2026-10-23) (numbers of their day); PRE 109 Fig. 10
  digitization (2026-10-24).
- Wall material, conversion and capsule energetics (2026-10-24) (*in part*); N221204
  late-time drive (2026-10-24) (`wall_history`; numbers of their day); Wall re-emission
  and the capsule radius (2026-10-24) (*in part*: `release`); Zylstra 2022 LEH pair,
  time-dependent (2026-10-25); Energetics defaults and the capsule share (2026-10-25)
  (*in part*; the DU factor is current); Late-time drive: the wall, the plasma and Dante
  (2026-10-25); The P2 swing: capsule radius, backscatter and the collision
  (2026-10-25); **Capsule absorption, conversion and the ablating-wall release
  (2026-10-25)** (the current energetics defaults and N221204 numbers); Capsule radius
  schedules: published histories and the API (2026-10-26).

**Spectra and diagnostics**
- Spectral radiosity, M-band and a synthetic Dante (2026-10-24); Dante channel responses
  and the unfold (2026-10-24); Laser-spot spectrum (2026-10-24) (α_M 0.55 with thermal
  spots); Per-cone M-band in the drive history (2026-10-27).

**Hydrodynamics coupling**
- Capsule drive hand-over and CompactLES stage 0 (2026-10-24).

### Where the current defaults are set

| Default (src) | Value | Section |
|---|---|---|
| `WallLossAlbedo` `release` | 1.0 | Capsule absorption, conversion and the ablating-wall release §1 (closure: Wall re-emission and the capsule radius §1) |
| `drive_history` `wall_history`, `wall_release` | on | N221204 late-time drive; Wall re-emission and the capsule radius |
| `DU_LOSS_FACTOR` (`du_wall_loss`, `uranium`, `nif_uranium`) | 0.836 at 300 eV, T^-0.2, 0.7 µm Au liner | Energetics defaults and the capsule share (factor); Wall material, conversion and capsule energetics §1 (model) |
| `gold()`, `uranium()` conversion | 0.75 | Capsule absorption, conversion and the ablating-wall release ("What changed") |
| `nif_conversion()` (`nif_gold`, `nif_uranium`) | 0.7 t^0.12 ≤ 0.95 | Wall material, conversion and capsule energetics §2 |
| `CapsuleAlbedo` | α_c 0.2, `outside` 0 | Wall material, conversion and capsule energetics §3 (0.2); Capsule absorption, conversion and the ablating-wall release (`outside`) |
| `WallMotion` | 0.185 mm/ns from 0.8 ns, stopping 0.45 mm in | Hohlraum dynamics calibrated to published data |
| `BubbleAbsorption` | speed 0.36 mm/ns, delay 0.25 ns, q 1.8, opacity 0.111 | Hohlraum dynamics calibrated to published data |
| `LipClosure` | 0.022 mm/ns at 100 eV, depth 0.24 mm, hold 0.4 ns | Computed LEH closure |
| `VisRadProfile` | growth 0.67, measured plates | VisRad deposition parity: hohlraum runs (2026-10-21) |
| `hohlraum()` capsule mesh, `Target` visibility | `:compact`; `:centroid` (`:sampled` for CAD targets) | Default capsule mesh and visibility |
| `Target`/`hohlraum()` `closed`, `cad_target`/`visrad_target` `closed` | `hohlraum()` true; others false (declared when apertures close the cavity) | Closure on bored, windowed and CAD cavities |
| `curved_components` / CAD `fit` | fit on with `remesh`, `minfraction` 1e-3 | "Fitting CAD by default" (`minfraction`) and "Parts of fitted CAD" (fit on), under CAD import |
| `LARGE_MODEL[]`, `HMATRIX_MODEL[]`, `HMATRIX_DEFAULTS` | 10,000; 15,000; `tol` 1e-6, η 1, leaf 48, maxrank 32 | Large radiation models; Hierarchical view factors; The H-matrix's N^1.5 remainder (`tol`) |
| `nif_mband()` | outer 0.20, inner 0.042 | Spectral radiosity, M-band and a synthetic Dante §2 |
| `spectral_solve` `spot` | `:radiation` | Laser-spot spectrum |

## Cross-beam energy transfer (CBET)

### Kinetic ray-pair gain

Source: Nguyen et al., arXiv:2304.13268 (2023), Eqs. 1–5. **H**

- `(v_a·∇) I_a = g_ab I_b I_a`
- `g_ab = −e² k² η_ab Im Γ(Ω, k) / (2 ε₀ mₑ² v_b ω_b² ω_a)`
- `Γ = χe(1+χi) / (1+χe+χi)`
- `Ω = ω₀ − ω₁`, `k = k₀ − k₁`, `|k| = 2k₀ sin(θ/2)`
- Small-gain form: `G₀ = k₀ L sin²(θ/2) Im Γ |a₀|²`
- Practical units (*derived*, linear polarization, vacuum intensity; the √ε factors
  cancel):
  `g [µm⁻¹] = 4.59×10⁻¹⁸ λ_µm sin²(θ/2) Im Γ I [W/cm²]`

Energy flows from the higher-frequency beam to the lower-frequency one, as seen in
the flow frame.

### Susceptibilities (Maxwellian)

Source: Turnbull et al. 2018, OSTI 1434633. **H**

- `χ_s = (kλ_Ds)⁻² [1 + ξ_s Z(ξ_s)]`
- `ξ_s = (Ω − k·u) / (√2 k v_Ts)`
- `v_Ts = √(T_s/m_s)`

Flow enters only through the Doppler shift `k·u`.

### Fluid / damping-rate form

Source: Follett et al., PRE 98, 043202 (2018), Eq. 1. **H** It agrees with the kinetic
gain to about 25% for a NIF LEH case (*derived*).

### Michel 2009 paraxial model

Source: OSTI 957176. **H** Its `γ = (Δk²/8k₀)Γ` is the same as `G₀` above; the
preprint misprints `Δk²` as `Δk`.

### Polarization

- With polarization smoothing and uncorrelated speckle, multiply the gain by
  `(1+cos²θ)/4`. Source: Michel et al., PoP 20, 056308 (2013). **M** At small θ this
  halves the transfer, consistent with Michel 2009.
- SSD bandwidth (2.2 Å at 1ω) increases off-resonance transfer.

### Saturation and calibration

The saturation clamps (δn/n) used to match NIF data vary widely:

| Era / design | Clamp δn/n |
|---|---|
| 2011 | 3–4×10⁻⁴ |
| 2018 | 5×10⁻⁴ to 3×10⁻³ |
| Pickworth 2020 | 0.01 |
| N221204 | 0.002 (foot), 0.0075 (peak) |

Sources: Michel 2011, Kritcher 2018, Pickworth 2020, Kritcher 2024 (PRE 109, 025204).
**H**

Designs also apply ad hoc multipliers on laser power, fit to in-flight P2, velocity
and shock symmetry, and degrade laser energy by 9–13%. Reduced CBET models must
therefore expose a gain multiplier as a calibration knob.

### LEH conditions for reduced models

Source: Michel 2009. **H**
- `n_e/n_c ≈ 0.06`
- `T_e ≈ 5–6 keV`
- CH plasma
- Flow about Mach 1 at the LEH throat
- Quad intensities: 2.5–5×10¹⁴ W/cm² (inner), 0.6–1.1×10¹⁵ W/cm² (outer)

*Derived*: with no flow, resonance sits at Δλ(3ω) ≈ 2 Å for a 20° crossing at
T_e 5 keV, T_i 1.5 keV.

## Wavelength detuning (NIF)

- **Sign convention.** `Δλ = λ_inner − λ_outer`. A positive Δλ transfers power from
  outer to inner cones.
- **Quoted wavelength.** Since about 2013, Δλ is quoted at **1ω** (Δλ₃ω = Δλ₁ω/3).
  Michel 2009–2011 quoted it at 3ω.
- **Briareus convention.** `Beam.dlambda` is at the beam's own (3ω) wavelength.
- **N170706 example (1ω).** 23.5° = 10530.5 Å, 30° = 10529.8 Å, 44.5°/50° = 10524.3 Å.

Published Δλ values and sensitivities:

| Design | Δλ | Sensitivity |
|---|---|---|
| NIC 2009 | 1.7 Å (3ω) for round | Inner +60%, outer −30% |
| HiFoot | 6.2 Å (1ω) | ~9 µm P2/Å |
| Low-fill CH | — | 1 Å ≈ +24 µm P2 ≈ cone fraction 33→41.5% |
| HYBRID-E (0.3 mg/cc) | 1–2 Å (1ω) | 40–50 µm P2/Å |
| N210808 / N220919 / N221204 | 1.8 / 2.5 / 2.75 Å (1ω) | Peak cone fraction 33% |

Rule of thumb: about 2.7 µm of hot-spot P2 per 1% of cone fraction (Pickworth 2020).
**H**

**Split Δλ.** One Be shot ran the 30° cone 0.7 Å shorter than 23.5° to correct
azimuthal asymmetry (arXiv:1701.09179). **M** The m=4 "cloverleaf" is attributed to
the four-fold geometry of the inner cone (Michel DPP 2024). **M**

## Backscatter

Practice: designers subtract the measured, time-dependent FABS/NBI backscatter from
each cone's input power after CBET (Michel 2011). Recent designs fold it into the
power multipliers instead (Kritcher 2024).

Typical magnitudes:
- **High-fill NIC.** Inner SRS reached tens of percent at peak. **M**
- **Low-fill and recent ignition designs.** Small, mostly early inner-beam SBS; 96–99%
  coupling. **H**
- **HYBRID-E per cone, time resolved** (N210808, Kemp et al., PoP 33, 082705 (2026)
  Fig. 2, `kemp2026_sbs()`): 2.8 kJ per 23° quad and 0.9 kJ per 30° quad, from the start
  of peak power to ~2 ns into it. See "Inner-cone SBS time histories and the late-time
  refraction (2026-10-25)". **H**

## Unconverted light (NIF)

From NIF User Guide §4.5 and §5.6.1. **H**
- A wedged focus lens separates the harmonics.
- The nearest 2ω edge is ≥ 2.9 mm and the nearest 1ω edge ≥ 4.8 mm from the 3ω aim
  point.
- At the 3ω focal plane the footprints are about 12 mm (2ω) and 22 mm (1ω) square.
- Unconverted light reaches out to about 32 mm from TCC.

The 1ω focus is 46 cm beyond the 3ω focus (MacFarlane 2003). **H**

*Derived* from fused-silica dispersion:
- 2ω focus about +27 cm beyond 3ω.
- Centroid offsets about 16 mm (1ω) and 9 mm (2ω), i.e. about 2.0 and 1.2 mrad.
- **Offset direction** (later search, below): the footprints fall to the left and
  right of the 3ω spot along ±φ, each column of a quad to its own side. Sources: NIF
  User Guide Fig. 4-12 (top/side views) and VisRad's "NIF Unconverted Light in 3ω
  Focus Plane" sample. **H/M**

## NIF quad names and bundles

Source: LMF chamber-geometry page,
https://lmf.llnl.gov/for-users/experimental-capabilities/chamber-geometry.

Naming:
- Four beamlines form a quad and two quads form a bundle.
- Quads are named by cluster digit, bundle digit and hemisphere, e.g. Q31T and Q31B.
- φ is measured counterclockwise from the axis through TCC and Switchyard 1.

How the quad-to-port mapping was extracted:
- Figures 3 and 4 (LLNL 2012-040636 and 2012-040642) are top views labeling the top
  and bottom quads.
- Each label's projected radius (∝ sin θ, which separates the cones) and azimuth was
  measured and assigned to the nearest port of its cone.
- All 48 assignments are distinct and within 9° of the tabulated azimuths.
- The figure's "FABS (130-186)_36B" label fixes the bottom view's orientation: Q36B is
  port 53.

The mapping is `NIF_QUAD_NAMES` in `src/facilities.jl`. **H** (*derived* from public
figures)

Structure of the mapping:
- 12 bundles are outer–outer and pair the same cone top and bottom, 4 are
  inner–inner, and 8 are mixed (30° in one hemisphere, 50° in the other).
- The port set is symmetric under 180° rotations about x, y and z (D2), and under
  nothing larger.
- The outer bundles form three D2 orbits:
  - {11, 23, 34, 46} at 50°;
  - {12, 25, 32, 43} at 44.5°;
  - {14, 22, 35, 45} at 44.5°.

## Pointing and alignment

- NIF pointing: 50 µm rms per beam, specified and met (physics/0111104; Appl. Opt. 50,
  1136, 2011). **H**
- Target alignment: about 10 µm. **H**
- Pointing range: ±30 mm up/down, ±5 mm left/right, ±30 mm along the beam. **H**
- Power balance: required ≤ 8% rms over the pulse; achieved 1.6% rms in quad peak
  power (Dunne 2012). **H**
- SG-III: about 30–80 µm rms. **M**

## VisRad

*Corrected:* the full manual is online at https://prism-cs.com/Manuals/VisRad/
(v21.1.1). An earlier search missed it. Formats documented there are listed in
"Later search" below.

Confirmed features:
- Workspace files (`.vrw`). The manual calls them JSON from v21.0.0, but real 20.0.0
  and 21.0.0 workspaces are keyed text; see "VisRad workspace files" below.
- STEP CAD import.
- A "Laser-Target Capture Clearance Report" (capture fraction and maximum incidence
  angle) exportable as ASCII.
- Multi-column ASCII time tables.
- Per-element flux exports for rad-hydro input (MacFarlane 2003).

Not documented in the manual: the `.vrw` layout. Two sample workspaces (2026-10-17)
settled the text format; see "VisRad workspace files".

## State of the art (pointing design and tools)

- **LLNL view-factor tools.**
  - A "dynamic view factor" with wall and spot motion taken from rad-hydro (Young, DPP
    2018–19). **M**
  - MacGowan's 3D view-factor metric (HEDP 40, 100944, 2021): 1% mode-1 flux gives
    about 100 km/s hot-spot velocity, and deliberate ±1.2% mode-1 drive was used to
    cancel capsule asymmetry. **H/M**
- **China.** IRAD3D adds time-dependent albedo (arXiv:1606.00934, 1703.01579). **M**
  The spot-motion claim recorded here earlier is not supported by the abstracts
  (2026-10-22 review). **L**
- **Ignition-era tuning (Kritcher 2024).** **H**
  - Late-time symmetry is set from data-driven models (Callahan) and measured
    P2-vs-Δλ curves, not first-principles CBET.
  - P2 swings from waist-hot to pole-hot during peak power.
- **Optimization.**
  - The polar-direct-drive repointing optimizer of Barlow/Colaïtis (arXiv:2311.18342)
    uses a genetic algorithm with 550 IFRIIT evaluations and does no robustness
    analysis. **H**
  - LLNL uses Bayesian and transfer learning (Humbird). **M**
  - No published gradient or adjoint illumination optimizer was found. **M**
- **Practitioner gaps.**
  - Static view factors.
  - Ad hoc CBET and backscatter calibration.
  - Manual pointing iteration.
  - No pointing uncertainty quantification for hohlraums.
  - Separate bespoke metrics for mode-1 and m=4.

## Beam clearance methodology

*Corrected:* NIF's TaLIS criteria **are** public, in "TCC Configuration Requirements"
(NIF User Forum, 2018), linked from the LMF Shot RI Resources page; see "Later search"
below. LLE's campaign criteria were not found. The material in this section predates
that finding, and the engineering defaults derived from it are superseded for NIF by
the TaLIS rules. What else is public:

- **Qualitative LLNL rule.** "Spaced from the laser beams appropriately to account for
  beam intensity profile as well as possible misalignment of laser beams and target
  components." Shield assemblies were held to ≤ 250 µm, preferably < 100 µm, of
  stacked tolerance (Boehm et al., OSTI 1342996). **H**
- **Design-era encircled-energy specs.** NIF and LMJ: 98% of 3ω energy within a 600 µm
  circle, i.e. about 2% outside the nominal spot (Caird et al., UCRL-JC-124872).
- **Why a 99%-enclosed envelope is unsafe** (*derived*). For a super-Gaussian, the 99%
  contour sits at 4–9% of peak intensity: about 4×10¹³ W/cm² for a 2 TW beam with
  n = 4.7 and r₀ = 300 µm.

### Error budgets

| Item | Value | Source |
|---|---|---|
| NIF pointing | 50 µm rms spec, met | Haynam AO 2011; Di Nicola 2024 |
| NIF beam-to-beam pointing | 5.2 µrad (≈ 40 µm) | Di Nicola 2024 |
| NIF mean beams-to-target | ≤ 25 µm | SPIE 8505 |
| NIF target vertical positioning | ±22 µm | Di Nicola 2024 |
| NIF focus error (1994 spec) | ≤ 1 mm | Fusion Tech. 1994 |
| OMEGA pointing | 9–20 µm rms | LLE Review v100 |
| OMEGA EP pointing | < 75 µm | LLE |

How reviews combine these errors (box vs. 3σ, linear vs. RSS) is not public.

### Far field

- **OMEGA measured** (Theobald, OSTI 1846973):
  - SG5-650: R1/e = 285.6 µm, n = 4.19.
  - SG5-850: R1/e = 362.5 µm, n = 4.74.
  - Less than 1% of the energy lies in the non-super-Gaussian wing.
- **NIF:** no public wing measurement.
- **Upstream scatter:** limited by the spatial-filter pinholes, about ±200 µrad
  (≈ 1.5 mm at TCC).
- **Final-optics scatter:** unfiltered, up to about 3 mrad (≈ 23 mm).
- **NIF SSD smear** (*derived*, low confidence): about ±30 µm along the dispersion
  direction at 45 GHz.

### Unconverted light

- NIF uses a wedged focus lens.
- Footprints are 1ω ≈ 22 mm square (nearest edge ≈ 4.8 mm from the aim point) and
  2ω ≈ 12 mm square (nearest edge ≈ 2.9 mm), reaching about 32 mm from TCC.
- Curvature rule inside the 1ω footprint: 1D radius ≤ 2 mm, 2D radius ≤ 5 mm, or
  dimpled (§5.6.2).
- Off-TCC pointing can send 1ω into near-opposite final optics assemblies (OSTI 13900).

### Briareus defaults

These are engineering defaults, not facility criteria:

- **Envelope:** 10⁻⁴-of-peak contour, plus the f-cone.
- **Wings** (`nif_wings`):
  - power-law halo [1+(r/r95)²]^(-3/2) with 1% of the energy, truncated at 1.5 mm;
  - scatter floor of 10⁻³ of the energy, out to 23 mm.
- **Error set:** pointing 50 µm + target 25 µm + hardware tolerance (user-supplied),
  added linearly, plus focus ±1 mm.
- **Status tiers:**

| Tier | Intensity | Intercepted fraction |
|---|---|---|
| fail | > 10¹² W/cm² | > 10⁻³ |
| review | > 10⁹ W/cm² | > 10⁻⁵ |

## Later search (2026-10-16): LMF Shot RI Resources and the VisRad manual

The sources are the LMF Shot RI Resources page
(https://lmf.llnl.gov/for-users/shot-ri-resources, including xlsx and pdf files) and
the VisRad manual (https://prism-cs.com/Manuals/VisRad/). Everything below was checked
by reading the originals.

### NIF beam table

- Source: `beam-port-locations-updated-10-15-2020.xlsx`, sheet "Quads Port Locations".
  It gives θ and φ for B111–B468, with bundle and TOP/BOTTOM.
- The table is transcribed to `data/nif_beams.csv`.
- Each quad's mean beam direction matches the NIF User Guide port table to within
  0.004°. That independently confirms all 48 quad names read from the
  chamber-geometry figures. **H**
- Chamber map: +x toward Switchyard 1 (φ = 0), +y toward the control room (φ = 90°),
  z up.

### NIF CPP spots

- Source: `phase-plate-table-updated-10-15-2020.xlsx`, sheet "Measured profiles".
- Fit function: `c + m exp[-½ ρ^n]`, with radii = σ (2 ln 2)^(1/n), i.e. half widths
  at half maximum.
- Data: TAS images of rod shot N151022-002-006 (2ω), averaged per cone.

| Cone | n | Major HWHM (mm) | Minor HWHM (mm) |
|---|---|---|---|
| 23.5° | 7.37 | 0.886 | 0.656 |
| 30° | 8.05 | 0.858 | 0.606 |
| 44.5° | 4.79 | 0.637 | 0.372 |
| 50° | 4.52 | 0.619 | 0.362 |

- These are now `NIF_SPOTS`, converted to 1/e semi-axes. **H**
- The "Specs" sheet gives design 80/90/95% radii and ellipticities.
- Major axis horizontal: inferred from User Guide plots and VisRad's convention. **M**

### TaLIS rules

Source: "TCC Configuration Requirements", NIF User Forum, 27 Oct 2018. **H** All
distances are measured orthogonal to the beam.

- **3ω capture.** The 10%-of-peak contour must be ≥ 0.5 mm from the substrate edge,
  with incidence < 70° to the target normal.
- **3ω clearance.** The 1%-of-peak contour must clear a component by > 0.5 mm, plus up
  to 0.5° of metrology/alignment/positioning tolerance. A failure is assessed for
  blow-by by the backscatter working group.
- **1ω clearance.** The 10% contour of the 1ω footprint must clear by > 0.2 mm, plus the
  same 0.5°. That contour is offset 4.45 mm from the 3ω aim point. ARC must not
  intersect a NIF beam footprint on a target feature.
- **Uncleared 1ω surfaces.** Radius of curvature < 2 mm (cylinders) or < 5 mm
  (spherical dimpling). Otherwise the specular 1ω reflection must not come within 20°
  of any open quad (FOA).
- **Hardware.**
  - Positioner clearance: 25 mm inside the CIVS NFOV (NFOV ≤ 120 mm for TANDM-124,
    < 150 mm otherwise), 35 mm outside.
  - Target to TAS: 5 mm above/below, 10 mm radially.
  - Diagnostics to 1ω light: a cone of 10 mm at TCC extruded to 25 mm at 1 m.
  - Diagnostics to the shroud swept volume: 10 mm.
- **How the 0.5° enters.** VisRad's clearance report implements fixed and z-varying
  enlargements of a cone test volume; the documentation does not state how the 0.5° is
  applied. Briareus's interpretation, a lateral allowance of r tan 0.5° with r the
  distance from TCC, is an assumption to confirm.

#### Re-read of the source, and `talis_check` (2026-10-16)

The document was re-read in full:
https://lmf.llnl.gov/sites/lmf/files/2023-05/TCC-configuration-requirements.pdf. It is
16 slides (PDF title "TaLIS Configuration requirements_22Oct2018"), with input from
S. Ayers, D. Kalantar and D. LaTray. The rule slides (3–6) are text only, with no
figures. **H** Details beyond the summary above:

- **Origin of the rules.** Slides 3–6 come from the "TaLIS Configuration Assessments
  Guidance Document". That document is not public; it was not found on LMF, OSTI or by
  web search. The hardware rules are "Alignment Team requirements (Sept 2018)".
- **0.5° wording.** "Clearance allowing for up to 0.5 degrees
  metrology/alignment/positioning tolerance from the target component". No pivot or
  formula is given. **H** (wording); how it is applied is still not found.
- **Consequences.**
  - Capture or 3ω clearance violated, or incidence > 70°: "must be assessed for
    blow-by by the backscatter working group".
  - 1ω not cleared: "evaluated for unconverted light management".
  - Undimpled surface whose specular 1ω reflection comes within 20° of an open quad:
    "requires additional review and may not be allowed".
  - 1ω light that clears is "no concern". **H**
- **TAS stay-in zone** (slide 14). During beam alignment the maximum target height in
  TAS is about 38 mm.
  - No component may extend more than 19 mm above the lowest aim point of an upper
    beam, or more than 19 mm below the highest aim point of a lower beam.
  - Violations "may require a different alignment sequence or ... beam selection".
  - TAS withdrawal envelope: a cylinder of radius 60 mm (TARPOS 90-239, TANDM 90-348) or
    29 mm (C-TARPOS 90-015), target height 38 mm. The cylinder axis is not stated.
  - **H**
- **Not on the summary slides.**
  - The 4.45 mm 1ω 10% offset sits 0.35 mm inside the User Guide's 4.8 mm nominal edge.
  - The User Guide (§5.6.2, 2024 rev.) words the curvature rule as "not larger than"
    2 mm / 5 mm. TaLIS says "< 2 mm / < 5 mm"; Briareus uses the strict form.
  - The User Guide adds: dimple half-side to radius "cannot be smaller than 0.25"
    (printed "0.25 mm"); horizontal planar targets have conjugate beams in the same
    hemisphere. **H**

Other LMF material checked:
- Brown Bag charts (LLNL-PRES-673980, "Other constraints"): unconverted light must be
  managed against reflections from flat surfaces. No numbers.
- Laser Performance Shot RI training deck (NP1314-C): describes TaLIS as a NIF Expert
  Group and shows counter-propagating 1ω via refraction (N150917-005). No clearance
  numbers.
- NIF User Guide §3.3: TaLIS reviews configuration, unconverted light and debris.

Context for the 0.5°, not the rule itself:
- Alger et al., "NIF Target Assembly Metrology Methodology and Results" (LLNL-JRNL-431151,
  Fusion Sci. Technol. 2010; https://www.osti.gov/servlets/purl/1227013). Ignition
  targets were held to x, y, z ± 0.5 mm, roll ± 1.0° and pitch ± 0.25° relative to the
  base, with about 0.02° metrology uncertainty. **H**
  - These are rotations about the target base or mount, which supports offering a mount
    pivot.
- A NIF shroud flexure mechanism (https://www.osti.gov/servlets/purl/1343009) reportedly
  has 0.5° rotational travel. **L**: from a search snippet only; the PDF did not
  download cleanly.

Implementation (`src/talis.jl`, `talis_check`):

- **Tolerance.** Default: lateral allowance r tan(0.5°), r = |x − TCC|. Alternatives:
  `pivot = :aim` (about the beam's aim point), `:max` (larger of the two), or a point
  such as the component's mount. It is added linearly to the `clearance` error set:
  pointing 50 µm + target 25 µm + hardware tolerance, and focus/(2F).
  - That double counts the positioning part, conservatively.
  - Angle rules get +0.5° on incidence, and −1° on specular-reflection angles (a tilt by
    α turns a reflection by 2α).
- **1ω footprint model.** *Derived*, **M**.
  - Each beamlet's 10% 1ω region at its 3ω focal plane is a 22.7 mm square: the User
    Guide's 22 mm, grown 0.35 mm per side so that its near edge is the TaLIS 4.45 mm.
  - It is displaced along φ̂ to the beamlet's own side of the quad axis (User Guide
    Fig. 4-12 top view).
  - Away from that plane it is the cone converging to the 1ω focus 460 mm beyond the
    3ω focus, with its centre on the chief ray from the lens centre 7.7 m upstream. The
    centroid offset is then 15.8 mm at 2.05 mrad, consistent with the earlier estimate
    of about 16 mm and 2.0 mrad.
  - Near edge = 4.45 + 0.027 s mm: closer to the 3ω axis upstream of focus.
  - Quad-level beams are expanded into their four LMF-table beamlets, so both sides
    are lit.
- **Rule interpretations to confirm.**
  - Substrate edge = free mesh edge.
  - "Diagnostics to 1ω": clearance 10 + 0.015 r mm from the 1ω footprint region, r from
    TCC.
  - Open quads default to all 48 FOAs.
  - The "sombrero" dimpled LEH shield is modeled as a flat annulus flagged
    `(:dimpled, r)`.

All from the manual. **H**

- **Workspaces:** `.vrw`; the manual says JSON since v21.0.0, real 20.0.0/21.0.0 files
  are keyed text ("VisRad workspace files" below).
- **Laser Beam Summary CSV columns:** Cluster, Beam ID, On/Off, Polar Ang, Azm Ang,
  Pt Ref-CS, Pt Geom, Pt Param-1/2/3, Focus offset, TC Pt Pos-x/y/z (cm), Power, DPR,
  profile model, R0 major/minor, minor-axis rotation.
- **NIF SST CSV columns:** "Beam", "Pointing: Beam Aim Point X (mm)", "Y (mm)",
  "Z (mm)", "Pointing: Beam Aim Point Beam Defocus (mm)". Chamber Cartesian about TCC,
  written without the defocus applied.
- **NIF CMT XML:** chamber Cartesian aim (cm, mm or µm), defocus, CPP and on/off.
- **Defocus sign:** positive moves the focus along the propagation direction, the same
  as NIF CMT's positive (lens toward the target) and as Briareus's `defocus`.
- **Power files:** a time column (ns) then one TW column per beam (headed B115, or
  Beam_10 on OMEGA).
- **Other exports:** `.bti` beam–target intersection files, and surface-element flux
  files.
- **Capture/clearance report defaults:** equal to the TaLIS rules.

### OMEGA (partial)

- Beams 10–69, numbered by cluster digit and beam digit (Users' Guide Fig. 2.5).
- P1 is at (0°, 0°); P2–P6 at θ = 63.44°, φ = 54/126/198/270/342°.
- P6–P7 axis: P6 = (63.44°, 342°). **M**
- About 18 beam numbers are tied to angles via VisRad screenshots.
- The full 60-beam assignment in the chamber frame was not found. Places searched:
  OMEGA Users' Guide mirror, VisRad manual, LLE Review, GitHub; the PI Portal requires
  a login.

### Still not found

Places searched: Precision Diagnostic System and HDR focal-spot papers, Dixit CPP
design papers, OSTI full texts (library redirect), GitHub code search, User Guide
§4.5–4.6.

- The NIF far-field halo or wing fraction.
- The unconverted-light energy fraction as a function of pulse shape.
- OMEGA EP port angles.

## VisRad interchange formats (2026-10-16)

Implemented in `src/visrad.jl`. Source: the VISRAD 21.1.1 manual,
https://prism-cs.com/Manuals/VisRad/ (abbreviated `VR/`). All 311 linked pages were
crawled and searched; the screenshots cited below were downloaded and read, since
several format details exist only in them. Confidence: **H** stated in the text or
legible in a screenshot; **M** inferred from a worked example; **L** assumption.

### NIF SST pointing CSV

- Exported by File | Export Laser Beam Data | NIF SST Beam Data. One row per beam:
  aim point x, y, z about TCC and the defocus. "The position corresponds to that
  without the defocus applied." (VR/files/laser_beam_data_files.html) **H**
- Exact headers: `Beam`, `Pointing: Beam Aim Point X (mm)`, `… Y (mm)`, `… Z (mm)`,
  `Pointing: Beam Aim Point Beam Defocus (mm)`. These are legible in the v21.0.0
  import-mapping screenshot (VR/intro/images/updated_import_sst_file.jpg, also
  VR/power_sources/images/import_laser_date_csv_custom.jpg), with Units = mm. **H**
- On import, VisRad sets the reference frame to "Target Chamber" and turns a beam on or
  off by whether the file has pointing data for it
  (VR/power_sources/import_laser_beam_data_csv.html). **H**
- NIF beam names are `B###`, e.g. B334, B421 in the same screenshot. **H**
- Not documented: delimiter quoting, preamble lines, and whether the SST adds columns.
  The reader finds the header by "Beam Aim Point" and ignores the rest.
- No example data rows are published. The tests use synthetic rows.

### Pointing and focus conventions

- The Reference Point is the aim. The focus offset is measured along the beam
  direction; a positive value moves the focus along the propagation
  (VR/power_sources/set_laser_beam_pointing.html). This equals Briareus's `defocus`.
  **H**
- The worked example confirms the sign: offset −0.03 cm puts the focus 300 µm back
  toward the port (VR/intro/nif/sample_nif_half_lasers.html). **H**
- NIF CMT defocus: "a positive value moves the lens towards the target"
  (VR/power_sources/import_nif_cmt_data_file.html). That is the same sign. **H**
- VisRad's beam axis runs from the port position on the chamber wall through the
  Reference Point. Briareus keeps the axis parallel to the port. The difference is
  ≲ 1 µm at focus for mm-scale aims. **H** (the formula is stated; the size is
  estimated)
- Pointing refers to the 3ω (blue) cone. The 2ω and 1ω foci follow from the NIF optics
  (VR/power_sources/nif_beam_params.html). **H**
- VisRad's NIF chamber frame matches Briareus's: the sample problem views Q15T from
  θ = 23.5°, φ = 78.75° (VR/intro/nif/sample_nif_3wrect_lasers.html), as in
  `NIF_QUAD_NAMES`. **H**

### Laser Beam Summary CSV

- The export writes "the table data, as they are displayed in the table". Only shown
  beams are included (unused beams can be hidden), in display order
  (VR/files/laser_beam_data_files.html). **H**
- The file can be re-imported as the "Standard" format, selectively: General, Position,
  Energetics, Spatial profile (VR/power_sources/import_laser_beam_data_csv.html). **H**
- v21 headers legible in VR/power_sources/images/laser_summary.jpg: `Cluster`,
  `Beam ID`, `On/Off`, `Polar Ang`, `Azm Ang`, `Pt Ref-CS`, `Pt Geom`, `Pt Param-1/2/3`,
  `Focus offset`, `TC Pt Pos-x/y/z`, `Powe…`. **H**
- The later columns are cut off in the screenshot. The page lists them as quantities:
  power, randomized power multiplier, DPR, profile model, R0 major/minor, minor-axis
  rotation (VR/power_sources/laser_summary_table.html). Their exact headers are
  unknown. **L**
- An older layout (VR/intro/images/laser_summary.jpg) uses `Beam #` (numeric), `Pt Ref
  CS`, `Pt CS`, `Pt X-1/2/3`, `Sp Prof` and `Ro`. The reader accepts both layouts.
- **Units.**
  - Distances are in cm (VR/intro/conventions.html; the Pointing tab shows "(cm)"). **H**
  - Focus offset is in cm: the dialog shows `Offset: -0.148 (cm)`. **H**
  - R0 is in µm. The spatial-profile dialog shows "(microns)", and the old table shows
    `Ro = 50` with "default spot size 50 microns" in the v6.5.0 release notes. **H**
  - The user can switch input widgets to mm in Preferences. Whether that changes the
    CSV is unknown. **L**
- **TC Pt Pos is the focus position, not the Reference Point.** The manual's example
  rows show it. Beam_10 (58.85°, 5.94°), Beam_28 (58.85°, 318.06°) and Beam_31
  (42.02°, 342°) all point at (0, 0, 0) in `Pb_Bi_foil` with focus offset −0.215 cm.
  Their TC positions differ by 0.215 cm × (port₁ − port₂) to 4×10⁻⁵ cm, and all three
  decode to one aim. Beams at the LEH origin with zero offset share one TC position, and
  Beam_40 at +0.01 cm along the LEH axis sits 0.0100 cm from it. **M-H** (derived, but
  checked to the table's 4-digit precision; this is a test in `test/visrad.jl`)
- Pointing angles in cylindrical or spherical form are shown in degrees since v9.0.0, and
  the chamber-frame TC positions were added to the table in v8.0.0
  (VR/intro/revision_versions_800to1020.html). **H**
- A time-dependent power shows as "T-<peak>". **H**
- **Spot profile.**
  - At focus, `I = I₀ exp(−[(x'/a)² + (y'/b)²]^{n/2})` with a₀, b₀ the 1/e major and
    minor radii (VR/power_sources/set_laser_profile_props.html, eq. images
    profile_eq_A/B). This is the same function as `SuperGaussianSpot`. **H**
  - Spot size can be defined at 36.8 % (1/e), 50 % or 5 % of peak; the conversion is
    r_f = a (ln 1/f)^{1/n}. **H**
  - NIF CPP radii have been HWHM since v10.0.0 (VR/intro/revision_versions_800to1020.html). **H**
  - Away from focus VisRad uses a(z) = a₀√((z/z_A)² + 1) and n(z) = n₀[1 + 0.2 z/z₀]. This
    differs from Briareus's f-cone convolution. **H**
- **Rotation.**
  - The rotation is that of the minor axis from vertical, positive clockwise; zero has
    the major axis "nearly horizontal", "coincident with the upper edge of the square
    beam as it enters the target chamber" (same page; revision notes). **H**
  - Which side "clockwise" is viewed from is not stated. Briareus assumes the beam's-eye
    view (looking along propagation): ψ_major = 90° − ρ. `rotation_sense = :target`
    gives 90° + ρ. **L — verify at work**
  - NIF spots (major axis along φ̂) have ρ = 0 either way.
- DPR ("On" splits an OMEGA beam into two cones 85 µm apart) is in the CSV. **H**
  Briareus ignores it with a warning.

### Phase-plate table

- Table at VR/power_sources/phase_plate_params.html. Radii are in µm: the "(mm)"
  headers are "µm" in the Symbol font. **H**
  - OMEGA, OMEGA EP: 1/e radii.
  - NIF CPPs: HWHM.
  - LMJ: radii at 3 %.
  - GEKKO: definition not stated; 1/e assumed. **L**
- The OMEGA table has a "# radii at 5 %" footnote that no row references.
- VisRad's NIF CPPs agree with Briareus's measured fits to ≤ 2 % in 1/e radius:
  Rev3-Inner 890/660 µm HWHM with n = 7.6, against 886/656 µm with n = 7.37. **H**
- Encoded as `VISRAD_PHASE_PLATES` and `visrad_spot`.

### Power files

- CSV layout: a `time (ns)` header, then one column per beam named `B115` (NIF; `BL115`
  is accepted) or `Beam_10` (OMEGA). Rows hold times in ns and powers in TW
  (VR/files/laser_beam_data_files.html). **H**
- The import screenshot shows a NIF file `N121207-001-999_LEH_beam_powers.csv` with
  `BL115…` columns (VR/power_sources/images/import_tdep_laser_powers.jpg). **H**
- Negative powers are set to zero on import. For NIF the powers go to the blue cones
  (VR/power_sources/import_tdep_beam_powers.html). **H**
- Pulse-shape library files are two columns, W versus s
  (VR/files/laser_pulse_shape_files.html). The column order is not stated; time first
  is assumed. **L**
- Beam-energy files: Beam, then kJ columns; only one column is read. **H** Not
  implemented.

### NIF CMT XML (not implemented)

- Documented content (VR/power_sources/import_nif_cmt_data_file.html,
  import_nif_cmt_pointing_data.html): **H**
  - target STP file per mount
  - beam on/off
  - chamber Cartesian aim and defocus, in "cm", "mm" or "um"
  - CPP names: "Rev3 Outer Cone 2W CPP", "Rev3 Inner Cone 2W CPP", "1.0 Rev1 Outer Cone
    2W CPP", "1.0 Rev1a Inner Cone 2W CPP", "400um CPP", "No CPP"
- Beams in the file are turned on in all three cones; others are turned off. **H**
- The XML element names, nesting and unit encoding are not documented anywhere in the
  manual (the TaLIS controller screenshot shows only the import list). No public sample
  was found (web search 2026-10-16). Without the schema a reader would be guesswork, so
  none was written; a sample CMT file from work would settle it.

### Workspace `.vrw` and `.bti`

- The manual says `.vrw` is JSON since v21.0.0 and "can be readily modified using
  scripting tools" (VR/files/workspace_files.html). Real 20.0.0 and 21.0.0 workspaces
  are keyed text instead; `read_visrad_workspace` reads that format (next section). A
  JSON workspace (perhaps 21.1+) is rejected with an error. **H** for the text format
- The only key the manual names is "Beam cone F number" (v21.1.0 notes).
- `.bti` (File | Export Laser Beam Data → Beam-Target Intersection Points) gives each
  beam's intersection points in the component's and the chamber's frames, in cm and
  degrees (VR/files/laser_target_intpts_file.html). **H**
  - Format ID 2 adds the beam name and port angles.
  - An optional block adds a(z), b(z), n(z), fSQ, DLW, DLH, M.
  - The column layout is undocumented. Not implemented.

### VisRad workspace files (2026-10-17)

Read off two real workspaces (not committed): a NIF hohlraum with bundled STEP CAD
(VisRad 20.0.0, 17 objects, 192 beams, ~20 MB) and a primitive hohlraum + foil
(VisRad 21.0.0). Implemented in `src/visrad_workspace.jl`, `read_visrad_workspace`. The
test fixture rebuilds the layout from rows copied from these files.

- **Layout.** CRLF text. Line 1 `Workspace Format ID  = 11`; sections `[Name]:` …
  `[End Name]`: Header Data (`VisRad Version_21.0.0`), Target Component Data
  (`# File format ID:  29`), Laser Beam Data (`#VF LaserParmFile: Format ID = 22` in
  20.0.0, 23 in 21.0.0, same columns), Diagnostic Port, Chamber Component, Dialog
  Settings, Snapshots, Movie Maker, Memorized Views, Power Source Data (per-element
  source table), Radiosity Data, Beam Capture-Clearance Data, TVS Viewer, SRF Beam,
  View Factor Data. **H**
- **Units.** Lengths in cm and angles in degrees regardless of the GUI's "Ruler units =
  mm" preference (aims of 0.5905 = 5.905 mm at the LEH plane). **H**
- **Objects.** `Begin Parameters for Object = <name>` … `End Parameters for Object:
  <name>`, `key = value` lines: `Object type`, `Length 1..6 value`, `Number of mesh pts
  1..4`, `Normals in/out`, RGB, `Object visibility`, `Reference coord. system name`,
  `Set position of ID`, `Coord. sys. geometry type` (0 Cartesian, 2 spherical (r, θ, φ)
  — a marker at (0.32, 90, 150); 1 cylindrical assumed (r, z, φ)), `Position 1..3`,
  `Tilt angle`, `Rotation angle`, `Spin angle`, albedo / x-ray conversion / power
  source / laser reflectivity values with model IDs, and the three `Include in …
  calc.` flags. Then ten `Time-Dependent … for <name>:  # table times = N` lines
  (N = 0 in both files; the row layout of N > 0 is not seen), `Number of Bounding
  Volumes`. **H**
- **Generalized meshes (type 13, CAD).** `Gen. mesh object origin: x/y/z` (= −Position
  in every object), `Number of vertices`, `Number of surfaces`, `Node positions:`
  (`index x y z`, cm, centered on the object's bounding-box centroid, 0-based) and
  `Surface element connectivity:` (`index v1 v2 v3 -1`, all triangles). The manual: a
  GM object's position defaults to its bounding-box centroid
  (VR/grid/import_export_cad_files.html). Placement x = R·node + Position is consistent
  with every object (all angles are zero in the file, so R is untested on meshes). **H**
  for zero angles
- **CAD normals are not reliable.** In the NIF workspace the inner hohlraum wall
  (r = 3.21 mm) faces away from the cavity on every triangle, and the next sheet
  (r = 3.24 mm) is mixed (556 out, 259 in). Edges are not shared between CAD faces
  (tens of thousands of open edges per part). VisRad: "Surface normal directions are
  not editable" and radiosity results for such objects "may be inaccurate" (same
  page). `visrad_target(...; two_sided = true)` is the workaround. **H**
- **Type codes.** 13 = generalized mesh; 11 = the "reference cone" VisRad makes for a
  multi-shape STEP import (all CAD shapes use it as their reference frame; manual,
  same page). 1 = sphere (Length 1 = radius, Lengths 2–5 = 0, 360, 0, 180 angular
  limits), 4 = rectangle (Length 1 × 2), 5 = cylindrical hohlraum (Length 1 = radius,
  2 = length, 3 = LEH radius) are inferred from parameters. **M**
- **Primitive grids** from VisRad's element counts (Power Source table): a rectangle with
  10 × 10 mesh points has 81 elements; a hohlraum with (11, 12, 4) has 216 = 10 × 12
  wall + 2 × 3 × 12 end cap (albedo 0.7) + 2 × 12 LEH (albedo 0, indices 192–215).
  Periodic azimuth counts elements, other directions count nodes. Radial spacing on the
  caps (uniform in r or area) and the azimuthal phase are not fixed; Briareus uses
  equal-area rings from φ = 0. **H** counts, **L** spacing
- **Orientation.** The manual names Polar, Azimuthal and Rotation angles applied as
  Rz(Rotation), then Ry(Polar), then Rz(Azimuthal) about the original axes, then a
  translation of the centroid (VR/grid/specifying_orient_angles.html); 2D objects lie
  in x-y with normals +z, cylinders along z. The file's Tilt/Rotation/Spin are read as
  Polar/Azimuthal/Rotation with right-handed rotations. Only one nonzero angle was seen:
  the foil's tilt 90° turns its normal to +x, facing the four beams aimed at it (ports
  at x ≈ +0.78). **M** for tilt, **L** for the other two
- **Directory tree.** `Begin Target Component Directory Structure` … `End …`:
  `Number of target components = n`, n names, `Number of child folders = m`, then m ×
  (folder name, subtree). **H**
- **Beams.** After `Cluster Names:` (one per line, blank-terminated) each beam is
  `Laser Beam Name = B113`, `Cluster Name = Ring_050`, one whitespace row of 53 columns
  named by the file's own `!` header, `Ref. Coord. Sys.:  Target Chamber`, `D toward Port
  ID:  P 0 - 0`, `Pulse Shape Name:`, GUI text items, and three `… laser power table
  values:  # table times = N` lines (blue, green, red). Columns: port index, θ, φ, port
  radius (cm), custom flag, cone half-angle min/max (deg), f-number, RGB, transparency,
  pointing method, distance to port, coordinate geometry, position 1–3, focus offset
  (cm), red/green power model, power on/off ×3 (blue, green, red), visibility ×3, power
  mode ×3, power (TW) ×3, random multiplier, energy-normalized flag, energy (kJ), use
  multiplier ×3, multiplier ×3, spatial profile model, spot size model, spot size (µm),
  elliptical flag, minor radius (µm), rotation (deg), rotation model, rotation angle ID
  (blank in both files, leaving 52 tokens), super-Gaussian exponent, envelope and lens
  radius ratios, DPR on/off. The rows do not line up with the header's dash runs, so
  they are split on whitespace. **H**
- **Pointing methods** are 0-based in the order of the Pointing tab
  (VR/power_sources/set_laser_beam_pointing.html): 0 TCC, 1 distance toward a
  diagnostic port (not supported), 2 position in a coordinate system. Both files use
  2, Cartesian, chamber frame. **H** for 2
- **Red/green power model** (the column after the focus offset). The manual's four
  models (VR/power_sources/set_laser_beam_energetics.html): None (P_TOT = P_B),
  Absolute Power (P_G, P_R in TW), Fraction of Total Beam Power (f_G, f_R; the user gives
  P_TOT and f_B = 1 − f_G − f_R), Relative to Blue Cone Power (R_i = P_i/P_B). Only
  under the fraction model is the power column the total. Both workspaces use model 2
  with columns 1.000 / 0.050 / 0.050, and the foil workspace's deposition table holds
  4 × 0.900 TW for four beams (green off, red on). So the codes are numbered 0–3 in the
  manual's order, and the fractions are subtracted whether or not those cones are on.
  The reader takes the blue cone power as P × (1 − f_G − f_R) under model 2. This is the
  right design: the laser's output does not depend on which cones a calculation
  includes or displays, so a drive tuned in total power keeps its 3ω power when the
  unconverted cones are hidden or switched off.
  **H** for model 2 (measured), **M** for the numbering of the others
- **Deposition check (foil workspace, 2026-10-17).** The `[Power Source Data]` section
  lists, per object, one row per element: index, direct laser power source (TW/cm²),
  uniform source, albedo, x-ray conversion efficiency, laser reflectivity, pre-reflection
  power, numbers of red and blue cones hitting, and the minimum overlap angles at 1ω and
  3ω. The rectangle's elements run with its local x fastest (the transpose of
  `rectangle_mesh`'s order). Values are element averages, not centroid values (centroid
  values differ from it by up to 13% next to the peak). Against Briareus, four beams with a 620/360 µm HWHM
  n = 4.5 spot aimed at the 4 mm foil, tilted ~40° (9 × 9 elements of 0.44 mm):
  - Totals: 3.6000 TW in both.
  - Elements above 20 TW/cm² (peak 371): within 1.4%, centre within 0.1%.
  - The nine brightest elements' power fractions agree within 0.001.
  - Faint wings (2–30 TW/cm²) differ by 8–15%, and VisRad's wings are wider along the
    major axis and narrower across it. Elements near 1 TW/cm² differ by 20–27%, and
    those below 0.01 TW/cm² by up to 80%. Elements under 30 TW/cm² hold 4.38% of the
    power in VisRad and 4.29% in Briareus.
  - An independent analytic projection (super-Gaussian evaluated in the beam plane,
    × cos ι) reproduces Briareus's single-cone deposition to 0.05% in the core.
  **H** (measured locally)
- **Power Source Data layout** (`visrad_power_sources`; NIF CAD and primitive
  workspaces, VisRad 20–21). After `[Power Source Data]:` a line "Power sources are
  [not ]up to date for time = t"; a workspace not computed has no tables. Then per
  object `#  --- Power Source Summary for <name>`, `Format ID = 6`, `Number of
  surfaces = n`, a three-line comment header (names on two lines, units on the third,
  columns 14 characters wide) and n rows indexed from 0. The Beam Capture-Clearance
  section repeats the same summaries in the primitive workspace, so only the Power
  Source section is read. CAD row order is the connectivity order: Briareus's
  deposition matches VisRad's on the N210808 wall element by element (2388 elements,
  correlation 1.0000, 0.9% rms with `VisRadProfile`). Hohlraum and sphere element
  order: not established. **H** layout and CAD order
- **The faint wings (follow-up, 2026-10-17).** Rays and envelope:
  - *Rays.* "Incident angle" in the capture data is the largest ray-to-normal angle over
    a grid of rays from the final optic to the focal spot, out to the I/Imax boundary
    (0.1 for capture) (VR/power_sources/laser_capture_report.html, beam_capture.jpg).
    With the beam's 1/e radius at the lens on the f-number cone (half-angle 1.38°), the
    10% boundary lies (ln 10)^(1/4.5) × 1.38° = 1.66° off axis. That gives 40.04,
    39.38, 44.65 and 44.06° against the file's 40.02, 39.26, 44.61 and 43.93°, so the
    beam geometry agrees. **M**
  - *Envelope ratio.* 1.53645 = (ln 1000)^(1/4.5): the radius, in 1/e units, where an
    n = 4.5 profile falls to 10⁻³ of peak. **H** (arithmetic)
  - *Cone counts.* The per-element "# Blue Cones Hitting Surf" counts fit a 25 × 24 polar
    ray grid (nodes 0…ρe, azimuths offset by half a step) with ρe = 1.65–1.66 in 1/e
    units, against an elliptical footprint. Only two diagonal corner elements miss, with
    VisRad 3 and the model 0. The fitted ρe is ≈ 1.53645 × (ln 2)^(−1/4.5) = 1.667, as
    if the envelope ratio scaled the HWHM radius converted to 1/e once more. Elements
    just past that envelope along the major axis (1.11–1.56 mm) get 5×10⁻⁵ TW/cm²
    against Briareus's 1×10⁻³. **M**
- **Wing residuals.**
  - The residual follows each element's distance from best focus along the beam, not
    its position in the beam. The foil is tilted ~40°, so the distance reaches ±0.6 mm
    across the spot. VisRad/Briareus is 1.01 at 0.10 mm (major-axis edge), 1.08 and
    1.14 at 0.17 and 0.37 mm (the diagonals), and 0.79 at 0.54 mm (minor-axis edge).
    The table is exactly point-symmetric about the aim, so the dependence is even in
    that distance.
  - VisRad's profile changes with distance from focus
    (VR/power_sources/set_laser_profile_props.html, eqs. profile_eq_A–C):
    - I = I₀(z) exp(−α), with α = f_sq(z)[(x/ΔL)^M + (y/ΔL)^M] + (1 − f_sq)[(x′/a)² + (y′/b)²]^{n(z)/2}.
    - f_sq = [z/(z² + a(z)b(z))^{1/2}]^Γ.
    - n(z) = n₀[1 + 0.2(z/z₀)].
    - a(z) = a₀√((z/z_A)² + 1), z_A ≈ 2a₀f# = 28 mm.
    - I₀(z) is renormalized to the envelope's power at each z.
  - ΔL(z), M(z) and z₀ are not given; the manual's text gives Γ = 1.0 ("currently set
    to 1.0") and f_sq ≡ 0 for a beam circular at the final optics (re-read 2026-10-21;
    the equation images, profile_eq_A–C.jpg, carry no alt text). a(z) and b(z) are
    negligible within 0.6 mm, but f_sq is not unless Γ is large: z/√(z² + ab) = 0.5 at
    0.3 mm and 0.76 at 0.6 mm. (The hohlraum runs below show f_sq = 0 for NIF beams.)
  - Briareus instead convolves the best-focus profile over the f-cone. That moves these
    wings by ~1% (`ncone = 1` against 13).
  - A blend with box half-widths c·a and c·b, a constant M, and I₀(z) renormalized
    (fit: Γ ≈ 10, c ≈ 0.94, M ≈ 7.7) takes the minor-axis edge from −21% to −2% and
    keeps the core within 1.4%. It leaves the diagonals at +8% and +14%. Three
    parameters fitted to ~20 informative elements: suggestive, not identified. **L**
- **Ruled out:**
  - Centroid or few-point element quadrature (vertices, 2×2/3×3 Gauss, Simpson,
    p × p grids).
  - Averaging the profile at the ray hit points.
  - Fixed 600-cell beamlets.
  - An elliptical envelope cut alone (1.5–2.0).
  - Global or per-beam spot rotations, including quad-parallel square beams (±1.6–1.9°).
  - Axis or exponent changes. These fit the wings only by moving the core by 3–6%.
- **Power on the tilted foil.** VisRad's foil total is 3.5999999 TW for 3.6 TW
  delivered, exact to the table's rounding (~1e-6). Normalizing per transverse plane
  as the manual describes, with the hyperbolic semi-axes and the envelope
  (`VisRadProfile`), gives 3.60041 TW on this ~40° foil. VisRad therefore also
  renormalizes each beam's deposited total, perhaps by its ray-grid capture fraction
  (1.00000 for all four beams here). **M**
- **What would settle it.** The VisRad runs listed in ROADMAP, "VisRad deposition
  parity". The intersection-point export (VR/files/laser_target_intpts_file.html:
  a(z), b(z), n(z), f_SQ(z), ΔL_W(z), ΔL_H(z), M(z) at each hit) gives the defocus law
  directly. Rerunning with the foil normal to the beams, where the distance from focus
  barely varies, separates the focal profile from that law.

#### VisRad deposition parity: hohlraum runs (2026-10-21)

VisRad 21 ran laser deposition on three workspaces built from Briareus geometry in the
thor repo (`visrad/workspaces/{n210808,n210808_reference,bu34_repointed}.vrw`): an
N210808-like hohlraum (Au wall r = 3.2 mm, 11.24 mm long, open LEHs; 2388 wall
triangles of ~0.6 mm, one VisRad element each), capsule, two window patches, foam, a
4 × 4 mm foil tilted 20° (3200 elements), sombreros; 192 NIF beams with the built-in
CPP spots, f/20.755, port radius 769.5 cm. Each file's `[Power Source Data]` holds the
per-element direct laser source. Comparison scripts: thor `visrad/compare_*.jl`; the
fits below used Briareus's own operator points with a free-form evaluator.

- **Manual (re-read).** VR/power_sources/set_laser_profile_props.html: hybrid model,
  "the analytic model provides the spatial profile for the intensity at a surface
  element, while the ray-trace component is used to test occultation ... and energy
  loss ... through holes". a(z) = a₀[(z/z_A)² + 1]^½, z_A = FL[(FL/(2a₀f#))² − 1]^−½;
  n(z) = n₀[1 + k(z/z₀)^a], a = 1, k = 0.2; f_sq with Γ = 1.0, zero for a beam circular
  at the final optics; I₀(z) normalized so the envelope holds the power at each z;
  envelope "typically ~2–3"; "the ratio of the lens radius to the value of d(z) at the
  final optics ... changes z₀" (d(z) undefined).
  VR/power_sources/set_direct_pwr_calc_props.html: the ray trace uses a grid of
  N_r × N_θ points within each cone (default 25 × 24). Preferences, Lasers tab
  (VR/additional_topics/set_preferences_data.html): AMR thresholds segment length/d(z)
  0.5, distance/d(z) 1.5. **H** (manual text)
- **Envelope ratio = (ln 1000)^(1/n).** The four NIF plates' file ratios 1.53645,
  1.49576, 1.27326, 1.28955 (n = 4.5, 4.8, 8.0, 7.6) equal (ln 1000)^(1/n) to 6
  digits: VisRad's default puts the cut-off at 10⁻³ of peak (`visrad_envelope`).
  **H** (arithmetic on four values)
- **Lens radius ratio** is 0.4 for every beam (the row column after the envelope). Its
  role cannot be identified from these runs (all beams share it). **H** value
- **Before (manual z_A, one envelope):** wall elements rms difference 10.2% of mean
  (n210808), correlation 0.9969, Briareus peak 858.8 vs 778.1 TW/cm² (+10%); totals
  within 1e-4. Elements dominated (> 80%) by one quad, binned in focal coordinates:
  VisRad/Briareus 0.93–0.97 at the spot centre, 1.1–1.7 at q = 1–1.3 (q = elliptical
  radius in 1/e units), i.e. VisRad's spots are wider on the wall (3–8 mm past focus).
- **Free-form profile fit.** With intensity P/(a(z)b(z)) g(q) and g a hat-basis
  function per cone (17 nodes, q ≤ 2.4), least squares over wall and window elements
  reproduces VisRad to 1.09% rms: VisRad's profile is a function of q alone per cone,
  cut sharply at the envelope. Against the manual's a(z), the fitted g follows
  exp(−(q/s)^m) with s ≈ 1.09, m ≈ 5.3 (n₀ = 4.5) and s ≈ 1.03, m ≈ 8.6 (n₀ = 8):
  mostly a larger spot. With the semi-axes grown as below, the fitted g matches the n₀
  super-Gaussian within the fit's noise (±5–10%) at every q < 1.2.
- **What it is not** (n210808 wall rms, each at its best value):
  - the square blend with Γ = 1 (a square of half-width c·a(z), c·√(ab) or a rectangle
    c·a(z) × c·b(z), M = 4–20): 6–35%; 6.0% only for the rectangle with c = 1, M = 4,
    another way of lowering the edge order. The NIF beams here behave as circular at
    the final optics (f_sq = 0);
  - n(z) increasing as the manual writes it (z₀ = 3–10 mm): worse (11–12%);
  - a lower exponent: n × 0.6 helps the outer cones and not the inner ones (inner best
    0.8–1.0); n decreasing with |z| (z₀ ≈ 2.5 mm) gives 5.7–6% but is not the manual's
    form;
  - element sampling (with the growth below): centroid 4.7%, vertices 13.0%, edge
    midpoints 3.4%, against 2.8% for the element average and 2.7% for a 4-subtriangle
    rule: VisRad's values are element averages, as on the foil;
  - spot orientation (axes swapped): 18–35%; the envelope cut (ρ, 0.95ρ, 1.05ρ, grown
    with a(z) or not): ≤ 0.1%.
- **Faster growth, z_A → (γ/ρ) z_A.** Scaling z_A (and z_B) by a per-cone factor, best
  factor 0.45, 0.45, 0.52, 0.48 for n = 4.5, 4.8, 8.0, 7.6 on n210808 and 0.44, 0.45,
  0.52, 0.49 on bu34_repointed. Times the envelope ratio ρ: 0.69, 0.67, 0.66, 0.62
  (n210808) and 0.68, 0.67, 0.66, 0.63 (bu34), one constant γ ≈ 0.67. The form
  γ(ln 2)^(1/n)/ρ spreads more (0.65–0.75); z_A ∝ a₀² fails (≥ 5.8%). One γ for all
  beams on n210808: 0.60 4.2%, 0.64 3.2%, 0.68 2.80%, 0.72 3.0%, 0.76 3.5%. Adopted
  γ = 0.67 (`VisRadProfile(; growth)`). Every beam has the default ρ = (ln 1000)^(1/n),
  so a dependence on ρ cannot be told from one on n. What a(FL) ≈ ρ R_lens/γ ≈ 1.5 ρ
  R_lens means in VisRad's terms (d(z), the lens ratio 0.4) is unknown.
  **M** (inferred from VisRad output: one constant, fitted on one workspace)
- **Validation** (`visrad_profile(ws)`: per-beam envelopes, γ = 0.67, through
  `DepositionOperator`; before → after):

  | case | component | rms diff | peak, Briareus vs VisRad (TW/cm²) | correlation |
  |---|---|---|---|---|
  | n210808 (fit) | wall | 10.21% → 2.81% | 858.8 → 774.4 vs 778.1 | 0.9969 → 0.9997 |
  | n210808 | windows | 2.46% → 0.94% | 315.2 → 308.6 vs 308.5 | 0.9974 → 0.9995 |
  | n210808_reference | wall | 10.21% → 2.81% | as n210808 | |
  | bu34_repointed | wall | 10.66% → 2.68% | 942.1 → 897.4 vs 915.9 | 0.9969 → 0.9997 |
  | bu34_repointed | windows | 3.80% → 1.41% | 386.5 → 384.6 vs 374.6 | 0.9899 → 0.9989 |
  | bu34_repointed | foil (±2.5 mm of focus) | 2.81% → 2.60% | 1039.1 → 1039.2 vs 1032.0 | 0.99991 → 0.99993 |
  | basic.vrw foil, 4 beams | elements > 20 TW/cm² | max 1.43% → 1.31% | | |

  Capsule-drive departures of bu34_repointed from n210808_reference (Briareus
  radiosity on each deposition; drive ratio, then δ1..4 in %): VisRad 0.9261, 0.006,
  0.092, 0.018, 0.046; before 0.9261, 0.006, 0.054, 0.012, 0.024; after 0.9261, 0.007,
  0.090, 0.019, 0.037. n210808 δ2: VisRad 2.576, before 2.574, after 2.573. δ2 moves by
  ∓0.008 for γ ± 0.05. **M**
- **What remains.**
  - Wall rms 2.7–2.8%. Outer-dominated elements carry 82% of it (4.5% rms; the outer
    peak ring at |z| = 2.7 mm is 3–4% low), inner-dominated ones 1.2%. Per-cone γ
    would give 1.1–2.0% per cone. δ4 stays 0.009 pp below VisRad's.
  - Totals on tilted surfaces: per-plane normalization puts 0.15% too much on the 20°
    foil (VisRad: 20.80008 TW for 8 × 2.6, exact) and 0.06% on the 40° basic.vrw foil
    (0.01% with the manual z_A). VisRad renormalizes each beam's total; not modeled.
  - Window patches in bu34_repointed get 1.0% more than in VisRad (n210808: 0.02%).
  - Near-focus wings on the foil (10⁻³–10⁻² of peak): ±10–40%, as before.
  - Occultation by the ray grid: VisRad puts ≤ 7e-4 TW/cm² on the LEH washers'
    elements, Briareus nothing.
- **VisRad beam-profile probe (2026-10-21).** Workspace thor
  `visrad/workspaces/profile_probe.vrw` (`build_profile_probe.jl`):
  - **Setup.** One beam per NIF quad (48, 1 TW each, the template's CPP rows), each
    focused 60 mm out along its port and deposited on its own flat plate normal to the
    beam, at |z| = 0–12 mm (inner cones) or 0–15 mm (outer) from best focus, before and
    after it. The plates have 2 × 44² triangles, laser deposition only.
  - **Fit.** `fit_profiles.jl` fits each plate with a rotated, shifted super-Gaussian
    averaged over each element (16 sub-triangles), cut at ρ. It leaves 0.0–2.5% rms;
    the fitted power is the beam's within 0.13% and the center is on axis within 1 µm.
    VisRad catches 1.00009 of each beam's power on its plate.
  - **Rotation.** VisRad's spot is turned from Briareus's axes by exactly the angle
    between the beam's φ̂ and its quad port's φ̂ (NIF User Guide angles): within
    0.02° on all 48 plates. That angle is 4.3–5.3° (23° cone), 3.2–3.9° (30°),
    1.9–2.3° (44°), 1.6–1.9° (50°), with its sign varying by beam. VisRad's rotation 0
    is the quad's horizontal ("nearly horizontal" in the manual). **H** (48 beams)
  - **Symmetry.** −z and +z plates give the same a, b and n to the fit's precision.
  - **Exponent.** n(z)/n₀ − 1 is linear in |z|: 0.010–0.012/mm for both inner plates
    and 0.015–0.017/mm for both outer plates. That is the manual's n₀(1 + 0.2|z|/z₀)
    with z₀ ≈ 18.5 mm (inner) and 12.7 mm (outer), about 0.45–0.50 of the manual's z_A
    of the major axis. The subagent's earlier n(z) test used z₀ = 3–10 mm, too
    strong. **H**
  - **Semi-axes.** They follow no hyperbola, neither the 1/e radius nor any radius at
    a fixed intensity level, with or without n(z). a/a₀ = 1.50–1.54 and
    b/b₀ = 2.04–2.11 at 15 mm (outer); 1.16 and 1.25–1.28 at 12 mm (inner). The minor
    axis grows faster in relative terms but slower than γ/ρ makes it (2.26–2.35 at
    15 mm). **H** (tabulated; the law is unknown)
  - **Not an f-cone average.** Briareus's f-cone deposition through the same fit
    lowers n with |z| (to 0.5 n₀ at 15 mm); VisRad's rises.
  - **Model** (`VisRadProfile`, `measured = true`). The four plates' tables
    (`_VR_MEASURED`, normalized at z = 0), interpolated linearly in |z| and extended
    linearly beyond the last node, with the normalization tabulated per node. The
    `.vrw`/CSV readers and writers apply the quad-frame rotation
    (`_vr_reference_angle`).
  - **Validation through `DepositionOperator`** (thor `compare_deposition.jl`):
    - **Probe plates:** element rms 0.1–1.6% (inner), 0.2–4.5% (outer, largest at
      15 mm), against 2.6–20% with γ/ρ.
    - **n210808 wall:** 0.9% rms (from 2.8%), peak 779.1 against 778.1 TW/cm²,
      correlation 1.0000.
    - **bu34_repointed wall:** 0.9% rms, peak 911.5 against 915.9.
    - **Totals:** wall 0.07% below VisRad (γ/ρ: 1e-4); windows 0.17% below; the 20°
      foil 0.37% above.
    - **Capsule modes:** absolute P2/P4 agree with VisRad's deposition to
      0.014/0.004 pp (bu34_repointed).
    - **Departures of bu34_repointed from the reference:** δ1..4 0.008, 0.104, 0.020,
      0.045%, against VisRad's 0.006, 0.092, 0.018, 0.046% (γ/ρ: 0.007, 0.090, 0.019,
      0.037). Swapping one component's deposition at a time puts 0.010 of the δ2
      difference on the wall, 0.001 on each window and none on the foil.
    - **n210808 δ2:** 2.575 against 2.576.
    **M** (measured tables, checked on two hohlraum runs)
  - **Still unknown.** The law behind a(z) and b(z) (the tables cover only the four NIF
    plates at the default envelope and lens ratio 0.4); whether the envelope follows
    n(z); what moves δ2 by 0.01 pp on the wall; the per-beam renormalization on tilted
    surfaces.
- **Spots.** NIF profile models 49–52 carry the built-in plate's numbers in the row:
  49 = Rev3-Inner-2w (890/660 µm, n 7.6), 50 = 1.0 Rev1a-Inner-2w (860/610, 8.0),
  51 = Rev3-Outer-2w (640/370, 4.8), 52 = Rev1-Outer-2w (620/360, 4.5), i.e. HWHM radii
  (VR/power_sources/phase_plate_params.html). The reader identifies plates by those
  numbers, not the model ID. Nothing in the file states a spot-size definition for
  custom profiles (`spot size model` = 0 throughout); 1/e is assumed. **H** plates,
  **L** custom
- **Checks on the NIF workspace.** Port angles match `nif_beams()` to 0.0056° (the
  file's 3 decimals). Outer cones aim at z = ±5.905 (50°) and ±5.723 mm (44.5°), 0.45 mm
  off axis; inner cones at ±4.078 (30°) and ±3.462 mm (23.5°), 0.67 and 1.22 mm off
  axis; 8 of the 44.5° beams are off; 1 TW total each, 0.9 TW in 3ω (red/green model 2,
  below). A ray trace of the 184 beams through all 15 CAD parts (210k elements)
  conserves power: at 0.9 TW per beam, 130 TW on the two hohlraum halves, 28.4 TW on the
  waist part, 6.7 TW on two side parts, < 1 mW on the LEH washers (first run at 1 TW per
  beam, scaled).

#### CAD import, remeshing and surface fits (2026-10-17 to 2026-10-24)

- **CAD windings (2026-10-17).** In the NIF workspace all 15 CAD parts are closed
  solids (no boundary edges once vertices are welded by position; 247 non-manifold
  edges in the largest), with 1–25% of the shared edges joining patches of opposite
  winding (signed volumes up to 50× the bounding box before repair). After
  `orient_mesh` (consistent, closed solids facing out of their material) the volumes
  are plausible (an LEH washer 6.2 mm across and 50 µm thick: 84 → 1.1 mm³), and the
  hohlraum halves' inner wall (r = 3.20–3.21 mm) faces the axis and the outer
  (3.24 mm and beyond) away from it; as read, both faced out. `coarsen_mesh` at
  0.3 mm reduces the hohlraum halves only from 7,300 to 6,300 elements each: their
  tessellation is long slivers that no merging shortens. **H** (measured locally)
- **CAD remeshing (2026-10-17).** `remesh_mesh` (facets split to edges ≤ size/3,
  conforming, then clustered into elements and re-triangulated) on the same workspace,
  one thread: the hohlraum halves (7,300–7,600 facets, 245 mm² each, inner and outer
  surfaces) give 4,240–4,310 elements at 0.3 mm in 0.4 s and 7,970–8,010 at 0.2 mm in
  1.0 s, about 10 triangles per element; areas agree to 5e-15, every element is within
  10° of planar, reaches at most `size` from its centroid, and has its centroid on it;
  flat and gently curved faces give elements of area ≈ size² (diameter/√area 1.4–1.5),
  while faces narrower than the size (wall ends, fillets, thin rims) give narrow
  elements: 20% of the elements but 0.5% of the area. The largest part (133k facets,
  7,200 mm²) takes 20 s (128k elements) at 0.3 mm and 47 s (236k) at 0.2 mm; a thin rod
  part, curved faster than 10° per element across its width, gives 14 elements per
  size² (pass a larger `angle` there). The design default `spacing = 0.85 size`
  (Poisson-disk seeds) was set on synthetic plates, tubes and cans, where it gives
  0.9–1.1 elements per size²; `rmax = size` (rather than 0.75 size) avoids splitting
  Lloyd cells into fragments (0.75 doubled the count). **H** (measured locally; the
  constants are engineering choices, **M**)
- **Cavity selection (`enclosure`, 2026-10-17).** On the same workspace (oriented,
  `coarsen = 0.3`: 155k elements, 212k facets), 64 diffuse rays per element and the
  gathered flux take 3.0–3.5 s on 8 threads (7.8 s on 3; kept set and fluxes bitwise
  equal). From the capsule center the kept set is only the capsule and the inside of
  its tent (4.4k elements, open area 5e-4 mm²): the two tent films are closed solids in
  the CAD and seal the capsule into a pocket. With the tent left out, 15.3k elements
  (486 mm²) are kept: the inner walls and LEH-side end faces of the hohlraum halves
  (≈137 mm² each), the LEH washers' inner faces, the capsule and fill tube, two inner
  parts, part of the waist band through the 20 µm gap between the halves, and 62 mm²
  of the mounting structure (of 7,200) near the LEHs. The open-area estimate is
  14.5 mm², the two LEHs' area (r ≈ 1.5 mm); the dropped surfaces would return
  0.8 mm² of the reference flux (0.12 mm² at `threshold = 0.01`, which keeps 21.9k
  elements; 140k at `threshold = 0`, through 5e-4 mm² of gaps). `visrad_target(ws;
  orient = true, remesh = 0.3, inside = (0, 0, 0))` (tent excluded) selects on the
  205k oriented facets (6.5 s), remeshes the kept facets and selects again: 6,950
  elements (358 mm²) in 11 s, against 24 s to remesh everything (154k elements) and
  7,300 elements (388 mm²) selected after. A single selection on facets keeps 10.7k
  remeshed elements (492 mm²): a long facet straddling the edge of the lit region is
  kept whole. The defaults (geometric albedo 0.9, threshold 0.05 of the reference
  flux, 64 rays) are engineering choices, **M**; the measurements **H** (locally).
- **The window's CAD path (BriareusGUI, 2026-10-18).** Engineering defaults, **M**: a
  CAD part reaching farther from TCC than twice the farthest aim point of the beams
  VisRad fires is drawn only (its facets as read, not remeshed) and left out of the
  calculations on import; on the NIF workspace that is the 133k-facet mounting part
  (to 53 mm) and two far parts (to 97 and 14 mm), leaving 20k elements at the
  standard 0.28 mm in 12 parts (import 4.6 s on 8 threads, against ~25 s to remesh
  the mounting part too). The drive then runs on the cavity around TCC (`enclosure`
  on the remeshed elements, one selection): 7,155 elements with the tent left out; a
  pocket (nothing escapes, at most two components bound it) is opened by leaving out
  the one of its three largest components that gives the largest open cavity, which
  finds the tent. A capsule role on the CAD capsule fits a sphere to its vertices
  weighted by area (the fill tube is ~7 % of the area): center at TCC to 1e-4 mm,
  r = 1.1333 mm, the outer face of the shell (the median vertex radius within 1.5 mm
  is 1.1333). Cavity + 3,120-element capsule (h = 0.09): 9,726 elements, scenario
  built in 93 s on 8 threads, 3.1 GB live; an `optimize_pointing` evaluation ≈ 0.4 s.
  The x-ray view factors are refused above the element count whose three dense N×N
  Float64 matrices (view factors and two cached factorizations) fill half the memory
  (`radiation_size`: 26.8k elements on 32 GB, *derived*); the drive is then laser
  only.
- **The window's cavity and stalls (BriareusGUI, 2026-10-18, measured locally, H).**
  Cavity selection on the window's remeshed elements against the library's faceted
  two-pass (`visrad_target(...; remesh = 0.28, inside)`: select on the oriented CAD
  facets, remesh the kept facets, select again), same objects, tent left out: one
  selection 0.8 s (2.2 s with the pocket search) and 7,155 elements, 360 mm², against
  5.3 s and 6,908 elements, 339 mm². The window keeps each part remeshed whole (it
  draws it, deposits laser on it and checks clearance against it), so the faceted path
  would add the remeshing of the kept facets to every cavity build; the 6 % extra area
  is whole elements straddling the cavity's edge (mostly a cap washer, 413 vs 295).
  Element-size edits: building the `Target` of the 12 parts and of the 200k-facet
  ghost took 1.2–1.5 s per edit (BVHs), now reused for identical components; a CAD
  remesh runs on workers with the part drawn from its facets meanwhile. The ~1 s stalls
  at the first drive result, the first coloring by a field and the first ghost were
  JIT compilation (Makie's heatmap colormapping, the colorbar's first drawing, the
  transparent mesh's first upload), not GC (≤ 0.26 s pauses at 3.3 GB live); the
  window now compiles them before it opens.
- **The window's target builds, capsule fit and drive end (BriareusGUI, 2026-10-18,
  measured locally in the scripted Qt window on a NIF CAD workspace, 8 threads, H).**
  `Target` of the 12 calculation parts (275k facets, 20k elements) 0.39–0.40 s, of the
  4 drawn-only parts (141k facets) 0.18–0.19 s; the components themselves 1–5 ms (cached
  parts). The CAD capsule's sphere fit 0.13 s (0.79 s at the first call: compilation).
  `fit_surfaces` on that capsule (38k facets, sphere and fill tube) took 133 s
  single-threaded, so the window's whole-part curved fit first checks that every
  vertex lies on one algebraic sphere or cylinder (ms): none of the workspace's 14 CAD
  parts is one surface as a whole (hohlraum halves: bore, walls, ends), so the fit is
  neutral there; a sliver tube (the synthetic workspace) gets its cylinder. The ~1 s
  GUI-thread gap at a drive's end was inside the window's refresh, not Qt: compilation
  of the drive's result path (selection, drive and tree refreshes: 1.0 s and 0.9 s with
  1.3–1.5 s compiled, the evaluation compiling on a worker meanwhile) and a full reset
  of the model tree's rows (0.2–0.26 s, every delegate rebuilt) at every row-count
  change. A precompiled drive (the package's workload) and row inserts/removals where
  the tree changes bring the largest gap at the drive's end from 1.2 s to 0.31 s, of
  which 0.23 s is a GC pause.
- **CAD surface fitting (2026-10-18).** `fit_surfaces` (region growing over patches
  sharing a vertex, least-squares plane/cylinder/sphere fits refined by Gauss–Newton,
  tolerance `max(1e-4 mm, 1e-3 min(R, extent))`, regions kept when their rms is under a
  quarter of it) on one hohlraum half of the NIF workspace (7,342 facets, 245 mm²), one
  thread:
  - Raw facets, 0.7 s: 64 regions, 77% of the area (75% curved). The inner wall is a
    cylinder of r = 3.2100 mm (480 facets, 86.7 mm², facing the axis, rms 0.28 µm, axis
    within 1.4e-8 rad of z), the next sheet r = 3.2400 mm (65 mm², facing out, rms
    0.25 µm, max 1.3 µm: not quite one cylinder), a step at r = 3.348 mm (22 mm²), the
    entrance-hole bore r = 2.15 mm about z, holes of r = 1.0 and 0.25 mm on horizontal
    axes, and a few large-radius spheres and cylinders (R = 23–640 mm, rms 0.3–0.9 µm)
    on faces that are neither planar nor (within 0.1 µm) any of the three; those are
    harmless (within tolerance) but split off as small components.
  - Remeshed at 0.3 mm (4,237 elements), fitted directly: 1.6% of the area. Points
    `remesh_mesh` adds inside the wall's facets lie up to 5.8 µm inside the cylinder,
    beyond the 3.2 µm tolerance. With the raw facets as `reference` (fit there, then
    each element assigned to the region whose surface its vertices lie within tolerance
    plus that region's own facet deviation of), 1.3 s: 74% of the area, the inner wall
    881 elements and 84.6 of its 86.7 mm² (elements straddling the wall's edges stay
    flat). **H** (measured locally)
- **Curved quadrature on remeshed elements (2026-10-18).** The cylinder/sphere
  quadrature (sub-triangle points projected radially, exact curved areas via
  `_curved_area`) is exact in area for any element whose facets are a radial graph
  over the surface, not only chord rectangles. Its cost was not: per-facet points made
  a 12-facet element 6× (36× per pair) dearer than a primitive one. Elements of more than
  two facets now lump the facets' pieces into cells of size/(2^level √2), as many as a
  two-facet element has. Remeshed elements also need their centroid on their facets:
  on a convex surface (a CAD capsule) the area centroid of an element folded across
  facets lies inside the solid, and its visibility rays were blocked by the element
  itself (every capsule element saw nothing). **H** (synthetic, see below)
- **Fitted vs flat CAD (2026-10-18), synthetic, Float32 STL round trip.** Interior of a
  latitude-longitude sphere (R = 2 mm, 60 × 120 bands, remeshed at 0.5 mm, 257
  elements): fitted, F_ij = A_j/4πR² to 5e-17 and an integrating sphere's incident flux
  uniform to 1e-15 with the source on five elements; flat, entries off by up to 1.4×
  the exact value (5e-3 absolute), rows by 0.37% on average, and the flux 10.5% low on
  the worst element; view factors 0.4 s against 28 s (8 threads). Closed can (R = 1, L = 2, 180 strips,
  remeshed at 0.25): wall–wall exchange −9.0e-4 flat, +4.6e-4 fitted (exact closed form);
  disk–disk +2.6e-3 for both (the planar disks are flat elements in both). Hohlraum
  (two thick halves with entrance-hole bores, 360 strips, a 56 × 112 capsule, remeshed
  at 0.3, 2,902 cavity elements): capsule rows to wall, caps and LEHs within 6e-4 of the
  exact solid-angle fractions for both; with Sinkhorn closure the capsule flux under a
  uniform wall drive is 2.097957 fitted, 2.098185 flat, 2.097943 for the primitive
  `hohlraum(h = 0.3)` (fitted 7e-6 off, flat 1.1e-4); view factors 3.5 s fitted, 83 s
  flat, 0.7 s primitive (8 threads). Integrated hohlraum quantities are thus already good
  on flat remeshed CAD once closed; fitting fixes the local exchange (concave self-view,
  the integrating-sphere uniformity) and the cost. **H**
- **CAD tutorial assembly (2026-10-18), synthetic, ASCII STL round trip** (docs tutorial
  "Radiation in an imported CAD hohlraum"). Thick-walled halves (R = 2.875, 50 µm walls,
  20 µm waist gap, 120 full-length strips), 50 µm end walls and washers (0.1 mm LEH
  bores), a 40 × 80 capsule, a two-cup tent, a mount; 30% of the facets reversed.
  `enclosure` from TCC with the tent: pocket only (open area 0); without: exactly the
  cavity surfaces; open area 14.99 mm² at 64 rays/element (−2.8%) and 15.33 at 512
  (−0.6%) against 15.42 (LEH 120-gons 15.09 + gap exit 0.37 − 0.04 intercepted by the
  mount); its spread over seeds is ~2% at 64 rays, 0.45% at 256. Under a prescribed
  wall source (closed by Sinkhorn), against the same cavity built from primitives with
  exact surfaces: at 0.5 mm (capsule 0.15) capsule flux within 0.06%, a₂/a₀ within
  0.12 points, LEH loss within 0.15%, flat and fitted alike; at 0.25 mm (capsule 0.1 vs
  0.05) within 0.012%, 0.07 points, 0.04%; a₂/a₀ itself moves 0.35 points from 0.5 to
  0.25 mm (10.92% at 0.2). Flat view factors 13.3 s vs fitted 1.4 s (1,684 elements,
  8 threads); with 180 strips per wall the flat ones take 22 s: the flat path integrated
  over every facet of a remeshed element (until "Flat CAD elements" below). The thin-walled `hohlraum()` loses 3.3% more
  through flush LEHs and has a₂/a₀ 0.77 points lower: the bores and gap matter.
  **Capsule mesh and P4:** the remeshed CAD capsule (compact elements) gives
  a₄/a₀ = −0.78 to −0.80% at 0.15 and 0.1 mm (408 and ~900 elements; walls 0.25–0.5);
  the latitude-longitude `sphere_mesh` capsule of the same cavity gives −1.01% at
  0.15 mm, −0.90% at 0.1 (walls 0.35), −0.82% at 0.05 (6.6k elements; −0.815% with
  0.2 mm walls); see "Compact capsule mesh" (2026-10-22). **Visibility:** with bored LEHs and a
  gap, row sums before balancing miss one by up to 1.3% on walls/caps/capsule and 8% on
  LEH apertures and narrow faces, in the primitives as in the CAD (`hohlraum()`: 0.14%):
  pairs partly shadowed by a bore rim, decided by one centroid ray. **H** (synthetic)
- **Flat CAD elements (2026-10-18).** A remeshed CAD element carries 6–46 facets (mean
  15 on the tutorial's walls, 23 on its capsule), and the flat path summed the
  point-to-polygon factor over each receiving facet and integrated each emitting facet
  separately (16 sub-triangles each, then adaptive). Now: (1) the receiver is the
  element's outline when the point sees all of it from in front with every facet facing
  the same way (the contour integrand is divergence-free, so the facets' contours sum to
  the boundary's; exact for folded elements too), otherwise piece by piece (facets in one
  plane to 1e-12; a piece cut by the point's tangent plane is closed by one segment
  between its two crossings, more crossings sum its facets); straight runs merged. Over
  110k points on neighbors and around elements of the tutorial's flat model (37k of them
  cut) it equals the facet sum to 8e-13; outlines have 16 edges against 52 facet edges.
  (2) Near-pair emitters lump their facets into cells of the tangent plane and refine
  cells toward the receiver. At the two-facet resolution (level + 1) the cells carry the
  bias that primitive quads have: perpendicular plates with a common edge (Hottel's
  closed form), all pairs on the near path, +2.0e-4 for strips remeshed at 0.25 and for
  quads alike, against ≤ 3e-5 for the per-facet rule; rooted one level finer (level + 2)
  +6e-7 to +3e-5 (pieces bisected under a child's side when refining: under half
  changed nothing at 1e-6 and cost 70% more). At the default `near` (mid/far pairs by
  centroid, unchanged) the closed forms move by under 5e-5 from main: plates +6.45e-4 and
  +4.34e-4 (main +6.28e-4, +4.13e-4), coaxial disks −3.05e-3 (−2.98e-3, 90-gon fans),
  a closed box of strips max |row − 1| 4.60e-3 (4.64e-3). Curved (fitted) elements of
  many facets use the same cells, rooted at level + 1 (the resolution of their point
  rule). Times (8 threads, main → now): tutorial flat model 12.4 → 2.7 s, fitted
  1.2 → 1.1 s (primitives at 0.5 mm, 2,244 elements: 1.0 s); the synthetic
  thick-walled hohlraum of test/cadfit.jl at 0.5 mm (1,217 elements) flat 47.1 → 4.1 s,
  fitted 2.3 → 1.3 s; integrated results unchanged (capsule flux, a₂, a₄, LEH loss to
  ≤ 2e-6 relative or 1e-4 points; capsule rows identical to 4 digits). A flat CAD pair still costs ~5× a primitive pair: folded
  outlines have ~3× a quad's edges, and the cells' bookkeeping. Primitive targets
  (1–2 facets per element) are bit-for-bit unchanged. **H** (synthetic)
- **Finely faceted curved CAD on the flat path (2026-10-18).** A slow drive reported on
  the GUI's synthetic sliver workspace (`sliver_workspace(4000)` in
  BriareusGUI/test/test_workspace.jl: a tube of 4,000 full-length strips, r = 3 mm,
  2,000 around, 9.4 µm wide, across the hohlraum) is the cost of the flat many-facet
  path, not a hang or runaway refinement: refinement stops at depth 8, no loop is
  unbounded. Remeshed at the coarse preset the tube gives 492 of the 802 elements, with
  274 facets, 72 pieces (strips are never coplanar) and 275 outline edges per element
  (the cuts across strips are staggered: a boundary between axial neighbors spans
  0.15–0.2 mm along the strips, about size/3). View factors: 54 s on 4
  threads (210 CPU-s), 99% of it in the 9,177 tube–tube near pairs: 1,800 cell
  evaluations per pair (max 34k: the staggered boundaries interleave neighbors over a
  band, not a line) at 13 µs each (the receiver's pieces and outline). The NIF CAD
  cavity at the coarse preset (3,740 elements, ~7 pieces per wall element, 35 on the
  capsule): 19 s on 4 threads, 443 evaluations per near pair, 1.7 µs each on the walls.
  With the tube given its cylinder (`curved_components`, 3 s) its view factors take
  2.4 s (22× faster; mean tube row sum 0.6185 vs 0.6168 flat); fitting the NIF cavity's
  parts changes its time little (19.3 vs 19.7 s). Deposition on the tube target: 8.5 s
  for 3 beams. **H** (synthetic; the NIF numbers measured locally)
- **Straight cuts across sliver strips (2026-10-18).** The stagger above came from
  dividing each edge equally: a strip's diagonal (10.0000044 mm on the tube) took 51
  points where its sides (exactly 50 × size/3) took 50, so the chords across every strip
  were slanted (267 of an element's 275 outline edges). `remesh_mesh` now cuts parallel
  edges longer than size/3 (within 0.25°; a strip's diagonal joins its sides' family
  within 30°, and a family whose strip facets pair it mostly with a heavier one joins
  that) where they cross common planes ⟂ the family, spaced ≤ (size/3)·cos 0.5°/1.1, and
  grows elements within bands between every m-th plane, about `size` long (bands of
  under 2 size² connected area are dropped); seeds in a band are spaced across it only
  (0.75 size²/length: random sequential seeding in a row leaves gaps of 1.34 times the
  spacing), a sliver strip's piece in a band (mean width under size/5) goes whole to the
  element holding most of it, and an element left in pieces is split. Standalone tube
  (2,000 strips, r = 3, 10 mm, the GUI test's), size 0.6, main → now: 492 → 523
  elements (0.94 → 1.00 per size²), 274 → 130 facets and 275.5 → 132.0 outline edges per
  element (2.1×; the floor is two per strip crossed, ~65), reach ≤ 0.99 → 0.82 size,
  diam/√A ≤ 1.80 → 1.44, smallest element 0.47 → 0.80 size²; view factors 30.3 → 10.6 s
  (8 threads), cell evaluations per near kernel 999 → 693 (1.44×) and 17.2M → 9.9M in
  all, mean row sum 0.7195 → 0.7221. Evaluations do not halve: an adjacent pair refines
  along the shared boundary down to cells of r/256, so its count is the boundary's
  length over that, which straight cuts shorten only by the stagger (the same tube cut
  into eight faces at random heights, where equal spacing happens to align, gave 603).
  Plate of 40 strips (test/meshops.jl): 28.5 → 3.2 facets and 30.1 → 4.0 outline edges
  per element (rectangles); can of 48 strips: 13.3 → 8.3 facets. Remeshed planar strips
  are thus 2–4-facet rectangles on the two-facet view-factor path: perpendicular plates
  remeshed at 0.25 match Hottel to +7.8e-4 and +1.9e-3, as 4 × 4 and 5 × 5 primitive
  quads do (+1.2e-3, +1.9e-3); the staggered many-facet elements, on the finer cell
  path, gave +3.9e-4 and +6.6e-4. NIF CAD workspace (local): the hohlraum halves'
  facets are mostly not long parallel strips at these sizes (34 of 245 mm² banded; most
  long facets have long edges ~90° apart), and change little: 4,637/4,728 → 4,684/4,716
  elements at 0.28 mm, 2,147/2,209 → 2,115/2,225 at 0.6; the GUI's coarse cavity 3,740 →
  3,790 elements, outline edges 16.0 → 16.1 per many-facet element, view factors 10.9 →
  11.9 s (8 threads), 266 → 269 evaluations per near kernel. CAD tutorial: halves 449/454
  → 470/470 elements at 0.5 mm; flat model 1,684 → 1,723 elements, outline edges 16.2 →
  13.8 (walls 16.3 → 11.8), view factors 3.0 → 2.7 s; row sums before balancing within
  1.3e-2 → 8.4e-3 (narrow faces 7.7e-2 → 6.4e-2); capsule flux against the primitives
  −0.057 → −0.050% (flat, 0.5 mm) and −0.012 → −0.019% (fitted, 0.25 mm), a₂/a₀ 10.550 →
  10.570% (flat 0.5), LEH loss 0.25356 → 0.25321. The constants (0.25°, 30°, the 0.05 g
  end margin, band length `size`, 2 size², size/5) are engineering choices, **M**; the
  measurements **H** (synthetic; the NIF numbers measured locally)
- **Cones (2026-10-18).** `ConeSurface(half_angle; apex, axis)` projects points radially
  at constant height along the axis, as on a cylinder; the patch area is
  (sin θ / 2cos²θ)|∮ z² dφ| (Green's theorem in height and azimuth), each edge by
  Gauss–Legendre (16 or 8 points) halved until the poles of 1/ρ² lie outside the
  Bernstein ellipse of parameter 5, with z measured from the lowest vertex so a slightly
  tapered cone does not cancel. The nearest-point projection was tried first: rim chords
  then dip toward the apex and a 12-strip 45° frustum covered only 97.7% of the true
  one (−2.3e-2, −9e-4 at 60 strips); at constant height rim chords land on the rims and
  faceted frusta (12–360 strips, 1–8 bands; 45°, 14.0°, 0.57°, 89.3°, and a pointed
  26.6° cone) give π(r₀ + r₁)s to 4e-15, a 0.01° cone to 3e-14, and single triangles of
  cones of 5°, 30°, 75° match a 200/400 subdivision of their projection (Richardson)
  to 8e-12. Cost 450 ns per triangle against 310 ns for the cylinder's closed form.
  The apex is a singular point: point rules on the ring of elements around it missed
  3.4% of its exchange (ring closure −3.4e-2 at 4, 8 and 16 bands alike), so elements
  touching the apex stay flat (+3.4e-3, as the flat cone). Closed 45° frustum
  (r 1 → 2, z 0 → 1, sector disks, 8 bands × 48), against the coaxial-disk factor and
  closure/reciprocity: curved cone→d1 +1.1e-3, cone→d2 −5.7e-4, cone→cone +2.5e-4
  (flat +1.0e-2, −2.2e-3, +4.1e-3; at 16 × 96 curved +8.6e-4, −3.6e-4, +1.5e-4, flat
  +4.0e-3, −1.3e-3, +1.0e-3). Pointed cone on its base: F(base → cone) = 1 − 5e-5.
  **H** (synthetic)
- **Cone fits (2026-10-18).** `fit_surfaces` grows cones as ρ = R + kz about an axis
  point (k → 0 is a cylinder, so the fit stays well conditioned for slight tapers):
  axis from the centered patch normals (on a cone they trace a circle about the axis;
  right for a band around it) or the cylinder's axis (right for a small patch), an
  algebraic fit of ρ² (with and without the k²z² term: on a band with only two rims z²
  is a combination of 1 and z), then Gauss–Newton with step halving; the lower residual
  wins. A cone whose apex lies over 1000 extents away is left to the cylinder, one whose
  curvature radius R/cos θ does to the plane. Tie-break order plane, cylinder, cone,
  sphere: a band with vertices on two rims fits a sphere and a cone equally, and is taken
  for a cone (a tessellator would put rings inside a spherical zone that bulged beyond
  its chord tolerance; a 45° one-band frustum fitted as a sphere has 5% more area).
  Synthetic Float32 STL frusta (one band at 39° and 45°, 8 bands at 2.86°, 3 bands at
  83.2°, raw and remeshed at 0.3 mm with the facets as reference): half-angle to 1e-6
  relative, apex to 1e-5 mm, one region holding every element. CAD frustum (45°, 180
  strips, fan disks, remeshed at 0.3): fitted cone→d2 −5.9e-4, cone→cone +3.5e-4,
  d2→cone −3.8e-4 against flat −2.7e-3, −1.5e-3, −2.7e-3; view factors 0.6 s vs 1.3 s
  (8 threads). Cylinder, sphere and plane fits, and every result on targets without a
  cone, are bit-for-bit unchanged (seven primitive targets; the tutorial's half, washer,
  capsule and tent cup fitted raw and remeshed). **H** (synthetic)
- **Cones in the NIF CAD workspace (2026-10-18, local).** Raw facets, one thread; time
  with cones vs without: hohlraum halves 1.1/0.8 s vs 0.7/0.6 s. Each half: the inner
  end wall beside the LEH is a 84.47° cone facing the cavity (apex 5.789 mm from the
  waist plane, radius 2.153–2.363 mm, 3.04 mm², rms 0.25 µm; it meets the r = 2.15 mm
  LEH bore), and a 61.1° chamfer at the waist edge of the inner wall (r 3.21–3.25 mm,
  0.19 mm², axis within 1.7e-3 rad of z). The LEH washers: a 77.14° chamfer on the outer
  rim (r 3.065–3.082 mm, 0.35 mm², facing out). A waist sleeve: a 10.0° bore chamfer
  facing in (2.35 mm²) and two 45.0° chamfers facing out (1.72 mm² each), rms 2e-10 mm.
  Cones also absorb pieces of free-form faces (fillets) within tolerance: on the halves
  8–10 more cones of 0.001–0.12 mm², most with axes tilted 0.03–1.5 rad, rms 0.5–1 µm
  (previously unfitted or spheres). Two halves through `visrad_target(remesh = 0.3,
  fit = true)`: 116 components (102 cylinders, 12 spheres, 2 rest) before, 138 with 22
  cones; curved area 353.6 → 360.8 of 490 mm²; view factors no slower (10.8 s before,
  7.9 s with cones, 8 threads) and identical on the inner wall (self-view, wall-to-wall,
  6 digits). **H** (measured locally)
- **Near-flat and small fitted regions (2026-10-18, local).** `curved_components` now
  takes `maxradius` (relative to the part's bounding-box diagonal, 11 mm for a half)
  and `minfraction` (relative to the part's area, 245 mm²); both off by default. Two
  halves, remeshed at 0.3 mm, with cones: maxradius 1, 2, 5 alone → 127, 128, 131
  components (the R = 23–670 mm spheres and cylinders are only ~10 of them; most of the
  138 are small holes and edge faces); minfraction 1e-3, 3e-3, 1e-2 → 25, 16, 10
  components with 356.6, 353.2, 344.1 mm² curved (of 360.8); minfraction 1e-2 keeps the
  inner walls, outer sheets, steps and the two 84.5° cones. **H** (measured locally)
- **Tori (2026-10-24).** `fit_surfaces` fits tori (`:torus`, in the default kinds) and
  `TorusSurface` carries them into the radiation model. Fit: candidate axes from the
  patch normals (the line their normal lines all meet, a linear complex: the two least
  eigenvectors of the form left after minimizing over the axis moment; the cone's two
  axes), each with a circle fit in the meridian half-plane, and starts from the patch's
  shape operator (spine through x₀ + n₀/κ₁, axis in the plane of n₀ and the tube
  direction through x₀ + n₀/κ₂, its tilt scanned: on a small patch the axis and R show
  only in third-order terms); the three best distinct starts are refined by
  Levenberg–Marquardt on √((ρ − R)² + (z − z₀)²) − r. Plain Gauss–Newton with step
  halving stalled on 37-patch starts, and of three equally cheap scan starts on such a
  patch only one reached the exact torus (the others stop at rms 1e-4 of r), hence the
  distinct-start rule. A 7-parameter fit needs stricter acceptance than the others:
  - free-form faces: the 20 × 20 paraboloid patch of the tests fits tori at rms 0.23 of
    the relative tolerance (accepted at the others' 1/4);
  - neighbouring faces: within the relative tolerance (rtol min(r, L), 5 µm for r = 10
    mm) a large torus (R 13.5, r 10 mm) took a NIF part's cone and 40% of a cylinder
    beside it at rms 16 nm, and a fillet torus 9 boundary strips of the cylinder it
    blends into;
  - band tessellations: strips with vertices only on their rims lie on any circle
    through three rims, and four rims placed symmetrically (two cone bands at either end
    of a cylinder) lie on one circle too: with tolerance `atol` the large torus above
    came back through such rims, at round-off, over the cone and the cylinder's strips.
  So tori are held to CAD round-off: tolerance `atol`, rms under `atol/16` (6 nm; the
  NIF CAD tori fit at 0.1–5 nm), at least four evenly spaced points along the tube
  circle (largest step under 3× the median past the empty arc), no plane, cylinder,
  cone or sphere fitting the region, and merges that keep that rms. A torus within
  2 tol of a sphere (R), cylinder (L²/8R) or cone/plane (L²/8r) is left to it, and a
  torus growth that becomes one, or that a simpler surface fits to round-off, is given
  up at its checkpoint; the patches of a rejected torus region (half its radius) grow no
  torus again, and a torus starts from at most 64 patches (rings around CAD fan apexes
  hold hundreds). Without these the cylinder–sphere–torus fit of the NIF capsule took
  103–123 s (torus growths retracing the sphere and the fill tube's fans); with them
  21 s against 18 s. A boxed closure in the LM step (`ntuple` over a reassigned matrix)
  had made each torus fit 10× slower.
  Synthetic (`test/cadfit.jl`): tilted and moved quarter-to-150° sections of R = 2,
  r = 0.5 mm, fine (12.5k facets) and coarse (19° bands), through Float32 STL, facing out
  or in, raw and remeshed with `reference`: one torus, R, r, centre to 1e-6, axis to
  1e-10. Fine fillet (all kinds): 7.9k / 32.8k / 65.5k facets in 0.13 / 0.51 / 1.6 s, one
  region (without tori 4 s at 8k and ~80 s at 33k, quadratic; "Cost of fitting faces of
  none of the kinds", below).
  `TorusSurface` projection (rim-preserving since 2026-10-24): in the point's meridian
  half-plane, inside the tube circle the nearest point of the circle (radially from the
  spine point), on or outside it the circle's point at the same height on the point's
  side of the spine (radially from the axis, as on a cone; the pole when above or below
  the tube). The two agree on the circle. A rim (circle ψ = ψ0 about the axis) has its
  chords at its height, displaced toward the axis by their sag δ: outside the tube for
  ψ0 on the poles and the inner half (distance² r² + 2δr|cos ψ0| + δ²), so they map onto
  the rim; inside it on the outer half, where only the outer equator (ψ0 = 0, the
  horizontal is the radius) maps exactly. Facet points outside the tube lie within the
  azimuthal sag of it horizontally (the meridional sag points inward), so the horizontal
  map moves them by at most that sag: second order like the nearest point (max distance
  from the nearest point on a quarter fillet 3.8e-2 / 9.6e-3 / 2.4e-3 / 6.0e-4 mm at
  16 / 32 / 64 / 128 strips, ratio 4.00). Exact rims everywhere are not reachable by a
  projection of points alone: exactness on every outer-half rim forces the horizontal
  map throughout the inside of the outer half, whose tangential error ε tan ψ (depth ε)
  is first order next to the poles (on a quarter fillet's cap strip, ~0.2 r Δψ). A
  facet-aware map (barycentric in (φ, ψ)) would be exact everywhere but changes every
  `surface_point` caller. Patch areas r |∮ (Rψ + r sin ψ) dφ| by Gauss–Legendre on edge
  pieces, split by bisection where an edge crosses the tube circle (the integrand kinks
  there); a fine subdivision of the projected facet agrees to its own discretization.
  Covered to round-off (|error| ≤ 7e-16): closed torus, ψ ∈ [0, π], [−π/2, π/2],
  [π/2, 3π/2], and the quarter fillet at 16 / 32 / 64 strips (before: +7.0e-3 / 1.8e-3
  at 32 / 64, the cap-rim chords projecting past the rim); ψ ∈ [0, π/3] (rim at 60°)
  +4.6e-2 / 1.1e-2 / 2.8e-3 (second order). The rounded hohlraum (rc = 1 mm) with its
  corners as tori: fillet area exact (the frusta's chord polygon −15% / −3.8% / −1.0%
  with 1 / 2 / 4 bands per corner), raw row sums max |1 − ΣF| 5.6e-3 / 3.0e-3 with 1 / 4
  bands (frusta 7.1e-3 / 3.0e-3), on the fillets 6.0e-4 / 1.2e-3 (frusta 4.7e-3 /
  2.2e-3; with the earlier projection 4.5e-5 / 1.8e-3, the one-band value a
  cancellation of the cap rims' +1.8e-3 area against the quadrature). Exchange against
  the frusta (Σ A_i F_ij): capsule–fillets −2.9e-3 / +1.5e-3 / +7.1e-4, fillets–LEHs
  −9.6e-3 / −2.9e-3 / −1.5e-3 with 1 / 2 / 4 bands; `hohlraum` keeps its frusta (the
  default 4 bands at h = 0.5; baselines). NIF CAD (14 parts other than the mount):
  fits bit for bit unchanged (fitting uses its own model, not the projection); the 16
  tori's 11.2k elements 177.64 → 177.74 mm² (+6.0e-4): the full quarter fillets
  (r 0.9 mm) from +2.1e-3 of their band to −2e-7, inner-half fillets (r 1.5 mm) from
  −2.4e-3 to ≤ 1e-4, thin tori (r 25 µm, inner-half rims, azimuthal facets long against
  r) +8–13% (a full one from −11% of its band to −1.3e-5), the convex round between
  cones (rims −29°, +23°) unchanged to 1e-15, and a fillet (r 1.05 mm) whose region
  takes a 0.4° sliver of cap past the pole −0.30%: facets lying nearly in the pole's
  tangent plane are squeezed by the horizontal map (one 1e-3 mm² facet −91%), though
  the band they tile is kept. Tori are not convex: their facets still occlude pairs
  inside them.
  NIF CAD workspace, its 14 parts other than the mount, raw facets, default kinds with
  and without `:torus` (and cylinder–cone–sphere, cylinder–sphere): 9 parts bit for bit
  unchanged; on the other 5 every region that changes is a fillet patchwork replaced
  whole, 29 small cones, cylinders and spheres (0.11–5.8 mm² per part) by 16 tori
  (R 1.6–3.4 mm, r 0.025–1.5 mm, 0.31–61 mm² per part, which also take the fillet facets
  that stayed unfitted), and on one part two spheres refitted (same patches, centres
  1e-16 mm apart). The capsule (fill tube) is bit for bit unchanged with every kind set,
  no tori: with the default kinds its tapered tube is two cones (0.147°, 0.110°) at
  round-off, with cylinders and spheres the 25 regions of "Cost of fitting faces of none
  of the kinds" (tolerance-stopped cylinder pieces). Fit time for the 14 parts 9.0 → 4.3 s
  (fillets no longer seed rejected growths); cylinder–cone–sphere 8.3 → 5.0 s;
  cylinder–sphere 21.6 → 25 s (capsule 18 → 21 s). The 133k-facet mount (default
  kinds, one thread): 1911 → 1948 regions, 1278 → 1100 s; 85 tori (347 mm², R 0.02–5.9
  mm, r 0.02–1 mm, rms ≤ 6 nm) take 54 of the 61 regions that change (39 cones, 19
  cylinders, 3 spheres; 7.4 mm²) whole, and every other region that changes or appears
  shares vertices with a torus: a 106 mm² cylinder (R 3.58 mm) that
  stayed unfitted before, and 0.3 mm² of small cylinders, cones and a sphere re-cut or
  newly fitted beside fillets. **H** (measured locally and synthetic)
- **Cost of `fit_surfaces` on finely faceted CAD (2026-10-18).** Not the sphere: a
  finely faceted sphere is one region grown from its first seed (lat-long spheres of
  10k and 40k facets, an icosphere of 20k: 0.01–0.07 s). The NIF CAD capsule part
  (38,200 facets, raw, one thread, local) has a fill tube tessellated in fans: vertex
  degrees up to 801, so the patch-neighbor lists (patches sharing a vertex) held 4.4M
  entries (116 per patch), and each seed's two-ring start walked 1,600 lists of 1,600:
  83% of the time. Neighbors are now visited through each vertex once per pass, with
  the same order of first occurrence, so results are bit-for-bit unchanged (every region
  of the workspace's 14 CAD parts, raw; the tutorial's parts raw and remeshed; synthetic
  frusta, capsules with tubes, tori). Capsule 0.50 → 0.11 s with the default kinds and
  121 → 20 s with `kinds = (:cylinder, :sphere)`; hohlraum halves 1.39/0.84 →
  1.31/0.79 s; the tutorial's fit step 47 → 30 ms; a synthetic fan-tessellated tube
  (19k facets) 0.37 → 0.05 s. Refits already come at every ×1.25 in size, so a growth
  costs linear time in its region; an incremental least-squares fit would not change the
  order. What remains super-linear is a face of none of the requested kinds: every one
  of its patches seeds a growth that runs until the tolerance stops it and is then
  rejected by the rms test (the capsule's 0.15° tapered tube without cones: 13,836
  rejected cylinder growths of ~1,000 patches; a synthetic tapered tube 12 s at 16k
  facets, 45 s at 25k; a 40k-facet torus 14 s with all kinds). Not seeding the interior
  patches of such a rejected region again cut the capsule to 4.4 s and the halves to
  0.5 s, but dropped regions on 5 of the 14 parts, among them a true 6-patch 56.9° cone
  (rms 0.1 µm) on a half, so it was not adopted. **H** (measured locally and synthetic)
- **Cost of fitting faces of none of the kinds (2026-10-24).** `fit_surfaces` now
  keeps rejected growths as covers: one of ≥ 512 patches covers, for its kind, its
  patches within half its radius of its seed, while all its patches stay available
  (an accepted region taking any of them lifts it). A growth of that kind from a covered
  seed is still run, with its rms recorded at 128, 256, … members, and given up at the
  first such checkpoint where its rms exceeds half the covering growth's there: a seed
  on the same free-form face retraces it, a true surface fits its patches to the
  faceting and grows on. The growth's own steps are
  unchanged, so a growth that is not given up ends exactly as before. Synthetic
  (`test/cadfit.jl` shapes): the 0.15° tapered tube (r = 50 µm, 2 mm) with cylinders and
  spheres, 2k / 8k / 18k facets 1.8 / 27 / 171 s before, after 1.5 s (8k) and
  3.3 / 6.5 / 30 s at 18k / 37k / 131k; patches tested 0.99M / 1.93M / 3.80M at 9k / 18k /
  37k facets (×1.96 per doubling; 41.8M at 9k without covers). A quarter-torus fillet
  (R = 2, r = 0.5 mm, all kinds) is unchanged, 4 s at 8k and ~80 s at 33k facets: its
  rejected growths stay under 512 patches (~200–500), so it is quadratic until they
  pass it (a torus kind would end that). NIF CAD workspace, its 14 parts other than the
  133k-facet mount, raw facets, kinds default / cylinder–cone–sphere / cylinder–sphere
  (42 fits, shared machine): 41 bit for bit unchanged (kind, surface, centroid, normal,
  patches, area, rms, max deviation, facing of every region); 8.3 → 8.6, 8.1 → 8.3,
  59 → 22 s (capsule alone 56 → 18 s); the 133k-facet mount with the default kinds
  also unchanged (1911 regions), 1348 → 1278 s on one thread (that cost is still open).
  The capsule with cylinders and spheres changes:
  its tapered fill tube's tolerance-stopped cylinder pieces (rms ≈ tol/4) come out cut
  differently (17 pieces of 0.151 mm² → 16 of 0.173 mm², 28 → 25 regions), and in their
  wake a 400-patch R = 7 µm cylinder (max 7 nm) of the tube's end is split between a
  piece and a near-flat R = 388 mm one (no planes in those kinds). Why these thresholds:
  with covers from 128 or 256 patches, or a fixed budget of 128 patches for covered
  seeds, 2–4 parts changed with the default kinds: regions that only just pass the rms test
  (rms 0.17–0.24 of the tolerance; fillet pieces fit as 30–300-patch spheres and cones)
  are found in the full search from seeds inside standing rejected growths, and the
  rms trajectory cannot tell them from the face around them; skipping covered seeds
  outright changed 6–12 of the 42 fits. Letting covered growths continue while their rms
  is at the faceting level (≤ 0.03 of the tolerance) is safe but brought the tube back
  to quadratic (16 s at 33k, 266 s at 131k facets): 128 patches there are one ring of
  the tube, too short to show the taper. The thresholds are patch counts, so a part
  faceted several times finer than the NIF parts would bring such marginal pieces
  under covers. **H** (measured locally and synthetic)
- **Cost of fitting large CAD parts (2026-10-25).** Profiled `fit_surfaces` (default
  kinds, one thread) on the 133k-facet NIF mount and on samples of it (the N patches
  nearest its centroid). Before: 40k facets 26–30 s, 73% of it the merge of regions of
  one kind (every pair rescanned after each merge, each pair recomputing both regions'
  extents: quadratic in the regions times the merges). The merge now tries each pair of
  region versions once (a pair's outcome depends on its two regions only) and keeps the
  tolerances, and the final fit of a growth reuses the last fit when the members are
  unchanged. Bit for bit unchanged (kind, surface, centroid, normal, patches, area, rms,
  max deviation, facing) on 90 fits: the 14 other NIF parts with four kind sets, the
  test shapes and a synthetic boss plate, and the mount samples; the full mount
  identical to the 2026-10-24 run (1948 regions). Mount 1100 s → 29 s; samples 20k /
  40k / 80k facets 7.4 / 26 → 2.9 / 6.0 / 17 s (now linear). Of the remaining time 53%
  is the start fits on the seed rings (tori 29%, cones 13%), 37% refits while growing,
  7% the merge. Comparisons must use one compile mode: the same code gives tori differing
  by ~1e-13 compiled in the package image and compiled in the session (likely FMA).
  Variable-radius fillet (spine ρ = 2 + r, tube r = 0.5(1 + 0.3 cos φ), all kinds):
  16k / 26k / 37k / 50k facets 16 / 38 / 71 / 107 s; cylinder and cone growths end
  tolerance-stopped (rms 0.41 of tol) at 200–370 patches, under the 512-patch covers,
  so the cost is quadratic in that window and drops to ~N^1.3 once growths pass 512
  (50k: 8,994 given up). No exact bound gives such a growth up early (the tolerance
  grows with the region's extent, so its rms in tolerance units can fall); smaller
  covers changed NIF parts before ("Cost of fitting faces of none of the kinds"), so
  none was adopted. **H** (measured locally and synthetic)
- **Fitting CAD by default (2026-10-24, local).** Should `fit = true` follow `remesh`
  in `cad_target`/`visrad_target`/the GUI (ROADMAP open question 1)? The NIF CAD
  workspace, its 11 parts near the cavity (no mount, far parts or tent),
  `visrad_target(...; orient = true, remesh = 0.6 or 0.3 (capsule half that),
  inside = (0, 0, 0))`, 4 threads with other jobs running. Fit cost: `fit_surfaces` on
  the parts' raw facets 7 s in all (halves 1.9–2.7 s; the 38k-facet capsule with its
  fill tube 0.29 s with the default kinds: the quadratic case of the previous item needs
  a face of none of the kinds), against ~220 s (0.6) and ~420 s (0.3) for the cavity
  selection and remeshing. Components: 10 parts → 158 (0.6) / 159 (0.3) fitted, 68 with
  `minfraction = 1e-3` (40 cylinders, 9 cones, 9 spheres, 9 rests, one part unsplit);
  the capsule's mean flux and moments with 1e-3 are within 3e-4 and 0.04 points of
  fitting everything (`:sampled` runs, below), so 1e-3 became `curved_components`' default. Left in with
  `maxradius` off: cylinders of R 54–670 mm and spheres of R 0.94–6.4 mm besides the
  capsule's (all of which `SphereKeepOut` would take by default). View factors, default
  `:sampled` visibility: 47 → 49 s at 0.6, 177 → 161 s at 0.3; with `:centroid`
  visibility at 0.6, 33 → 22 s. The 22× of the sliver tube (above) is a whole-part fit,
  which the GUI already makes (`_whole_surface`). Accuracy against
  `montecarlo_incident` on the flat facets (α = 0.7, a source on the inner walls,
  4,000 paths per receiver), `:centroid` visibility: outer capsule flux flat +0.03%,
  fitted +0.26% (MC `:element` 0.5096 ± 0.0003; the fitted capsule is the true sphere,
  and `:hit` gives +1.3% over `:element`); halves (400 elements) flat −0.81%, fitted
  −0.22%. **With the default `:sampled` visibility the
  outer capsule's row sums are 0.21 (0.0–0.99) against 0.955 with `:centroid`, and its
  flux 7× too low, flat or fitted; the halves' rows 0.79 against 0.91 (flux −15%).**
  `_visibility_samples` averages each grid cell's sub-triangle points, and on an
  element folded across the facets of a convex part (13.5 facets per capsule element)
  that mean lies behind the facets, where every ray is blocked (the failure `Target`
  avoids for the centroid by moving it onto the facets). Fix before relying on
  `:sampled` CAD targets. **H** (measured locally) Fixed the same day: next entry.
- **Sampled visibility on folded elements (2026-10-24, local and synthetic).** Cause
  confirmed on a public synthetic cavity (test/partial_visibility.jl: two `cad_half`
  halves remeshed at 0.6 mm, a 56 × 112 lat-long `cad_sphere` capsule of 0.9 mm at
  0.3 mm, 48 facets per capsule element, 978 elements, LEH apertures): 3,608 of the
  capsule's 3,951 samples (91%) sat behind its facets (0–0.9 µm, mean 0.24 µm; the ray
  offset ε is 0.55 nm) and were blocked along their own normal; no centroid was.
  Fix in `_visibility_samples`: on an element of more than two facets that is not
  planar, each cell mean moves along the normal onto the facet the line meets (the
  rule of `_onto_patch`), or to the nearest sub-triangle point of its cell if it meets
  none. That alone gave capsule rows 0.9974 against 0.9986 for `:centroid`: every
  remaining blocked sample ray of the reduced capsule rows (561k of 561k) hit the
  emitter's own facets, rays toward the horizon of a sample's own facet that the
  element model (radiating from the mean plane) already counts. So an element of
  more than two facets no longer occludes its own sample rays (nor the receiver its
  own). Primitives keep `:sampled` bit for bit (`==` on the bored cavity at 0.5 mm and
  `hohlraum(h = 0.8, capsule_h = 0.3)`). Raw factors (`closed = false`), a two-ring
  source on the driven walls, against `montecarlo_incident` (`:element`, 8,000 paths
  per receiver):

  | synthetic | capsule rows | capsule flux vs MC | driven wall flux vs MC |
  |---|---|---|---|
  | `:centroid` | 0.9986 (0.9963–0.9999) | −0.45% | −0.41% |
  | `:sampled` before | 0.027 (0.000–0.999) | −97.5% | −4.7% |
  | `:sampled` fixed | 0.9991 (0.9974–1.0005) | −0.33% | −0.32% |

  (MC 0.51716 ± 0.00029 and 0.55423 ± 0.00025.) NIF CAD cavity at 0.6 mm (as in the
  previous entry, flat, α = 0.7, 4,000 paths): outer capsule rows 0.16 → 0.961
  (`:centroid` 0.956), flux −86% → −0.31% (`:centroid` +0.03%; MC σ 0.05%); halves
  (400 elements) rows 0.79 → 0.915 (0.910), flux −15% → −0.67% (`:centroid` −0.81%).
  View factors 56 s against 27 s for `:centroid` (4 threads, shared machine). On the
  NIF capsule `:sampled` is now 0.3% further from the reference than `:centroid`;
  its capsule rows' own change adds +0.04% at fixed radiosity, so the shift comes from
  the walls' rows. The driven inner wall (ρ 3.0–3.4, 3,000 paths, σ 0.07–0.08%) shows
  `:centroid`'s capsule figure to be a cancellation: waist / mid / end bands
  `:centroid` +0.33 / +0.39 / +0.33%, `:sampled` −0.18 / −0.13 / −0.29%, so `:sampled`
  is closer on every wall band and uniformly slightly low (the raw factors' loss), while
  `:centroid` is high on the walls and low in the capsule's rows.
  **H** (synthetic, measured; NIF measured locally)
- **Parts of fitted CAD (2026-10-24, local).** With the split made computational
  (`Component.part`; name-taking functions and reports work on parts, ROADMAP open
  question 1), fitting became the default when `remesh` is given. The NIF CAD cavity of
  "Fitting CAD by default" (11 parts selected, 10 in the cavity, `remesh = 0.6`, capsule
  0.3, `inside = (0, 0, 0)`, default `:sampled` visibility with the folded-element fix,
  4 threads, shared machine): build 253 s flat, 243 s fitted (the 8 s of fits are inside
  the noise of the selection and remeshing); 10 parts, 68 components fitted (none of the
  parts is one surface: the capsule carries its fill tube), 3,086 elements either way.
  View factors 48.2 s flat, 48.7 s fitted. Against `montecarlo_incident` on the flat
  facets (α = 0.7, the two-ring source on the inner walls, 4,000 paths per receiver,
  `:element` re-emission): outer capsule (259 elements) flat −0.31%, fitted −0.06% (MC
  0.50960 ± 0.00027); halves (400 elements) flat −0.67%, fitted −0.07% (MC 0.45706 ±
  0.00031). The GUI's route (its flat window target, the drive fitting each part on its
  whole file facets): 8.4 s, 65 components, view factors 45.8 s, capsule −0.02%,
  halves −0.06%. Reports on the fitted target: `energy_balance` keyed by the 10 parts
  (68 with `by = :component`), components.csv 10 rows and no piece names,
  `SphereKeepOut` default the capsule part alone (9 spherical components),
  `show` "10 parts (68 components)". So fitting costs nothing measurable here and cuts
  the flux error 5–10×. **H** (measured locally)

### Result tables and fits

- **Surface Element Properties.**
  - Saved via File | Save Data for Object / All Objects / Selected Rows as
    "multi-column ascii". `*.sed` files from a run have the same contents
    (VR/results/display_surf_elem_table.html, VR/simulation_times/config_sim_params.html).
    **H**
  - Positions are in the component's own coordinates, (r, θ, φ) for spheres. **H**
  - Headers in the screenshot: `Z (cm)`, `Azim. ang (deg)`, `Radius (cm)`,
    `Surf Index`, `Emis T (eV)`, `Rad T (eV)`, `Flux In (TW/cm2)`,
    `Flux Out (TW/cm2)`, `Pwr Src …` (VR/results/images/surf_elem_table.jpg). **H**
  - The file's delimiter and header layout are undocumented, and the polar-angle
    header for spheres is not shown. **L**
- **Power Summary** (`*.psd`): object rows with `Surface Area (cm2)`,
  `Incident Rad Pwr`, `Emitted Rad Pwr`, `Pwr Src Term`, `Laser Pwr Dep`,
  `Direct Laser Pwr` (VR/results/images/power_summary_table.jpg). **H** The file
  layout is undocumented. **L**
- **Fits** (VR/results/display_surf_elem_table.html). **H**
  - Legendre: least-squares fit `F(θ) = Σ Aₙ Pₙ(cos θ)` up to A4, each term selectable.
  - Real spherical harmonics up to L = 4, with `N_LM = √((2L+1)(L−M)!/(L+M)!)`, √2 for
    M ≠ 0, cos for M > 0 and sin for M < 0.
  - This normalization gives Y₀₀ = 1, so A₀₀ is the mean; the example shows mean 9.736
    and A₀₀ = 9.739.
  - It also gives A_L0 = A_L/√(2L+1); the example shows 1.143·√5 = 2.556 against
    A₂ = 2.563.
  - Not stated: the fit weights and the Condon–Shortley phase. `visrad_legendre_fit`
    and `visrad_ylm_fit` are unweighted and have no phase. **L**
- In VisRad, the grid size of Sphere, Cylinder, Cylindrical Hohlraum and similar
  objects is the outer wall; the inner wall lies inward by the Wall Thickness
  (VR/intro/revision_history.html, v21.0.0). **H** The benchmark kit's README says so.

### To verify with real files at work

1. Export a Laser Beam Summary CSV from a NIF workspace with elliptical CPPs:
   - the exact headers after "Power", and whether exponent and spot definition appear
   - whether "TC Pt Pos" is the focus (as derived above)
   - cm units for the focus offset with mm preferences set
2. Set one elliptical beam to rotation +30° and check the spot overlay. Is the minor
   axis clockwise from vertical when looking along the beam (Briareus's `:beam`) or
   from TCC (`:target`)?
3. Export a NIF SST file. Check the header strings, any preamble or extra columns, and
   whether off beams are listed.
4. Save Data for Object on a sphere. Check the delimiter, the header spellings
   (especially the polar-angle column) and whether rows carry an index.
5. Obtain a sample CMT XML to implement the CMT reader.


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
- Held fixed: the density of the expanded foil (1 g/cm³ implicitly, through the
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
- In time: settings found at the peak overshoot before burnthrough (ΔP2 −1.3 to −1.4%,
  m = 2 1.1% of the opposite phase, 5 ns to burnthrough). Per-quad power trims switched
  on with the burnthrough (NIF quads take their own pulse shapes) hold P2 and m = 2 at
  the no-window values throughout (m = 2 ≤ 0.04%), for +3.3% laser energy.

## Hohlraum mesh work-precision (2026-10-16) *derived*

`examples/mesh_convergence.jl`: the GUI's NIF hohlraum (10.1 × 2.875 mm, LEH 1.55,
capsule 0.9, wall albedo 0.85, conversion 0.8), NIF quads at the published cone aims
("nominal") and at an optimizer-like pointing with the inner cones pulled to the waist
("tight"). Elements are sized with `hohlraum(; h, capsule_h)`. Build = deposition
operator + view factors + flux response, 16 threads.

Wall/cap element size `h` (capsule 0.07 mm), Legendre moments of the capsule flux in
percentage points:

| h (mm) | N | build | P2 nominal | P2 tight | P4 tight |
|---|---|---|---|---|---|
| 0.80 | 3,900 | 3.6 s | −14.64 | −4.54 | −1.127 |
| 0.60 | 4,276 | 2.4 s | −14.56 | −4.29 | −1.168 |
| 0.45 | 5,028 | 3.3 s | −14.60 | −4.14 | −1.196 |
| 0.35 | 5,796 | 3.9 s | −14.57 | −4.06 | −1.203 |
| 0.28 | 7,260 | 5.4 s | −14.54 | −3.97 | −1.205 |
| 0.22 | 9,428 | 12 s | −14.54 | −3.94 | −1.212 |
| 0.17 | 13,076 | 24 s | −14.52 | −3.90 | −1.216 |
| 0.13 | 19,572 | 58 s | −14.53 | −3.90 | −1.216 |

Capsule element size (wall 0.2 mm):

| capsule h (mm) | N | P2 nominal | P4 nominal | P4 tight |
|---|---|---|---|---|
| 0.30 | 7,208 | −13.94 | +1.238 | −0.687 |
| 0.22 | 7,372 | −14.20 | +0.814 | −0.909 |
| 0.16 | 7,656 | −14.38 | +0.508 | −1.057 |
| 0.12 | 8,160 | −14.49 | +0.325 | −1.167 |
| 0.09 | 9,056 | −14.52 | +0.245 | −1.205 |
| 0.07 | 10,452 | −14.53 | +0.216 | −1.213 |
| 0.05 | 13,620 | −14.55 | +0.185 | −1.216 |
| 0.04 | 17,232 | −14.55 | +0.177 | −1.226 |

Findings:
- The wall converges as about h^2.5–3 in P2 (tight: 0.64, 0.25, 0.075, 0.005 pp at
  h = 0.8, 0.45, 0.28, 0.17 against 0.13). The nominal design wanders by ±0.03 pp
  without a trend below h = 0.45.
- The capsule mesh dominates P4: a piecewise-constant flux on facets converges as
  capsule_h² (P4 nominal extrapolates to ≈ +0.16 %; error 0.08 pp at 0.09 mm, 0.6 pp
  at 0.2 mm). The old "standard" capsule (24 × 32, ≈ 0.12 × 0.18 mm) was 0.15 pp off.
- The mean drive and laser power fractions converge to 1e-4 already at h = 0.45.

GUI presets chosen from these (errors against the finest meshes / extrapolation):
`:coarse` h = 0.6, capsule 0.2 (≈ 1,300 elements, < 1 s; P2 ≈ 0.4 pp, P4 ≈ 0.5 pp);
`:standard` h = 0.28, capsule 0.09 (≈ 5,900, ≈ 4 s; ≲ 0.1 pp); `:fine` h = 0.17,
capsule 0.05 (≈ 16,000, ≈ 40 s; ≲ 0.02 pp, converged for design work). Build times
are from before the 2026-10-17 performance pass; a fine `Scenario` (16,244 elements,
closed) now takes 15 s on 24 threads (view factors 4 s, LU 7 s), and an optimizer
iteration about 60 ms.

## Beam-error budgets (2026-10-17) *derived*

`mode_uncertainty` and the beam-errors tutorial: NIF hohlraum (10.1 × 5.75 mm, 3.1 mm
LEHs, 0.9 mm capsule, h = 0.45, capsule_h = 0.12), 48 quads at 8 TW, independent errors
per quad of 50 µm rms pointing (displacement length, transverse to the beam) and 3% rms
power.
- Per-beam sensitivities from three dual partials per beam agree with central
  differences of full evaluations to 2e-9 (relative).
- Linear propagation agrees with a 300-shot Monte Carlo of full evaluations to 1–2% in
  every mode l = 1–4 (the sampling error).
- Published pointing: rms of the mode amplitude (√Σ_m (δc_lm/c00)²) 0.28, 0.21, 0.13,
  0.03% for l = 1–4; pointing and power contribute about equally to l = 1 (0.20% each).
  ε₁ = mode_amplitude: median 0.43%, 90th percentile 0.67% (≈ 67 km/s of hot-spot
  velocity at MacGowan's 100 km/s per percent).
- After tuning the cones (l = 1–4 nulled at fixed drive): 0.23, 0.16, 0.09, 0.02%.
  The robust objective (expected mode power) started from the tuned design does not
  move it: with cone-level pointing, independent beam errors leave a floor set by the
  error sizes, not by where the cones point.

## Retuning to a reference drive (2026-10-17) *derived*

`retune_drive` and the retune tutorial: the same NIF hohlraum and mesh with
`WallLossAlbedo` (wall albedo 0.65 → 0.94 through the pulse), tutorial pointing, a
reference pulse from a total power (60 TW picket, 15 TW foot, 50 TW step, 420 TW peak)
and a cone fraction (0.35 in the foot to 0.25 at the peak, an engineering example, not a
shot), drive histories at 0.4 ns steps to 9.6 ns. Bundle 31 (Q31T, Q31B) dropped.
- Bundle down on the reference pulses: ε₁ 2.02% averaged over the peak (7.6–9.4 ns),
  5.2% at most (in the foot, where the dropped pair carries a larger share of the
  inner-cone drive); Tᵣ −3.3 eV at the peak. ε₂ and ε₄ barely change (9.17 → 8.99%,
  0.25 → 0.20%: the tutorial pointing is not tuned for P2).
- Retuned per quad, each within 1.3× its cone's reference power per quad: ε₁ 0.023%
  at the peak and 0.18% at most, ε₂, ε₄ and Tᵣ back to the reference's (ΔTᵣ +0.01 eV).
  The fit holds the albedos at the reference's; a drive history of the new pulses,
  with albedos following the drive, reproduces the predicted c_lm/c₀₀ to 3e-5, so the
  linearization is not the limit.
- The unconstrained retune is not unique and uses its bounds: 25 of 46 quads at 1.3×
  at the peak, some near zero. Weighting the departure from the reference powers
  (`stay`, per quad relative to the mean reference power):

  | stay | ε₁ peak | ε₁ max | powers / reference at 8.8 ns | quads at the limit |
  |---|---|---|---|---|
  | 0 | 0.023% | 0.18% | 0.00–1.30 | 25 |
  | 1e-5 | 0.12% | 0.25% | 0.82–1.30 | 6 |
  | 1e-4 | 0.34% | 0.45% | 0.95–1.30 | 2 |
  | 1e-3 | 1.24% | 1.55% | 1.01–1.14 | 0 |

  A dropped bundle's mode 1 is fixed mostly by the quads around it; a stay of order
  1e-4 keeps most of the correction with powers within a few percent of the reference
  for all but a few quads.
- Limits: the per-quad caps here are relative to the reference powers, not NIF's
  per-quad energy or power capability (which depends on the pulse and the laser's
  state). The view factors and pointing are the reference's through the pulse.

## HYBRID-E (N210808) hohlraum and LEH size (2026-10-20)

Source: Kritcher et al., Phys. Rev. E 106, 025201 (2022), Table I, Secs. II A and II C
(full text, MIT PSFC copy). **H**
- Hohlraum: Au-lined depleted uranium, 6.4 mm diameter, 11.24 mm long, 0.3 mg/cm³ He.
- LEH diameter: 3.64 mm on N210207 and N210307, 3.1 mm on N210808. Callahan's HEDS
  slides (2022) quote 3.65 → 3.1 mm, 27% less area. **H**
- Capsule: HDC, 1.05 mm inner radius, ablator 76–80 µm (density 3.5 MCD / 3.32 NCD);
  the tutorial's 1.13 mm outer radius is *derived*.
- Peak laser power / energy: N210307 487 TW / 1.909 MJ; N210808 441 TW / 1.917 MJ.
  The smaller LEH "was used to operate at lower peak power ... while still achieving
  similar radiation temperature" (qualitative). Δλ 1.55 → 1.8 Å (1ω); the 50° beams
  moved toward the midplane to clear the smaller LEH; peak cone fraction 33%.

`docs/literate/leh_power_balance.jl` and `test/leh_power_balance.jl` (*derived*):
- Inputs that are engineering estimates, not shot data: wall albedo 0.925
  (`WallLossAlbedo` at 300 eV, 6 ns; the Au-lined DU wall is not modeled), capsule
  albedo 0.3, the x-ray source as uniform bands at 0.6 ≤ |z| ≤ 1.8 mm (inner, 33%)
  and 3.4 ≤ |z| ≤ 4.6 mm (outer). Moving the outer bands ±0.5 mm moves the fractions
  by < 0.005.
- Holes are open and static; LEH closure and the gold bubble shrink the effective hole
  during the pulse, so the computed LEH loss is the open-hole bound.
- At 3.1 mm: 34.4% of the x-ray source escapes through the LEHs, 23.7% is absorbed by
  the capsule; per 0.1 mm of diameter over 3.0–3.6 mm, +1.6% escape, −0.5% capsule.
- Briareus (`closed = true`, h = 0.4) against a Monte Carlo photon walk in the exact
  geometry: within 1e-4 (10⁸ histories, σ 3.5e-5). The uniform-flux power balance
  reads the LEH loss 6.5–7% low and the capsule's share ~1% high.
- Required-power ratio 3.1 vs 3.64 mm at fixed capsule drive: 0.890 (α_w 0.925),
  0.905 (α_w 0.90), against 441/487 = 0.906. This is consistent in size only, since the
  shots also changed Δλ, pointing and pulse length. **L**

## Compact capsule mesh (2026-10-22) *derived*

`examples/capsule_mesh.jl`: the hohlraum and the nominal and tight pointings of the
mesh work-precision study above, wall fixed at h = 0.2 mm, closed, 16 threads. The
capsule is either the latitude-longitude `sphere_mesh` (`hohlraum(capsule_mesh =
:latlong)`, the default, sized by `capsule_counts`) or the geodesic `icosphere_mesh`
(`:compact`: icosahedron faces divided into n × n triangles projected onto the sphere,
20n² elements, edges within −16/+19% of 1.107R/n, n = ⌈1.107R/capsule_h⌉). Both carry
the `SphereSurface`. N is the capsule's element count; moments in percentage points.

| mesh | capsule h | N | P2 nominal | P4 nominal | P4 tight | build |
|---|---|---|---|---|---|---|
| lat-long | 0.30 | 200 | −13.938 | +1.238 | −0.687 | 12 s |
| lat-long | 0.16 | 648 | −14.383 | +0.508 | −1.057 | 9 s |
| lat-long | 0.12 | 1,152 | −14.488 | +0.325 | −1.167 | 10 s |
| lat-long | 0.09 | 2,048 | −14.524 | +0.245 | −1.205 | 11 s |
| lat-long | 0.07 | 3,444 | −14.533 | +0.216 | −1.213 | 15 s |
| lat-long | 0.05 | 6,612 | −14.546 | +0.185 | −1.219 | 22 s |
| lat-long | 0.04 | 10,224 | −14.549 | +0.177 | −1.226 | 35 s |
| compact | 0.30 | 320 | −14.386 | +0.203 | −1.149 | 9 s |
| compact | 0.22 | 500 | −14.454 | +0.190 | −1.185 | 9 s |
| compact | 0.16 | 980 | −14.499 | +0.177 | −1.211 | 9 s |
| compact | 0.12 | 1,620 | −14.518 | +0.171 | −1.221 | 10 s |
| compact | 0.09 | 2,880 | −14.532 | +0.167 | −1.228 | 12 s |
| compact | 0.07 | 4,500 | −14.539 | +0.165 | −1.231 | 15 s |
| compact | 0.05 | 8,000 | −14.544 | +0.164 | −1.234 | 26 s |
| compact | 0.04 | 12,500 | −14.547 | +0.163 | −1.235 | 45 s |

Findings:
- Both converge as c/N (∝ capsule_h²) to the same limits, extrapolated from each mesh's
  two finest points: P4 nominal +0.161 (lat-long) / +0.162 (compact), P4 tight −1.237 /
  −1.237, P2 nominal −14.553 / −14.551. The mean drive is within 2e-4 of its limit at every size on both.
- P4 error ≈ 160/N pp (nominal) and 120/N (tight) on the lat-long capsule, 16/N and
  24/N on the compact one: within 0.02 pp of the limit takes ~8,000 lat-long elements
  (capsule_h ≈ 0.045) against ~800 (nominal) to ~1,200 (tight) compact ones
  (capsule_h ≈ 0.12–0.16). The compact mesh is 5–10× cheaper for P4, consistent with the
  remeshed CAD capsule's 408 elements (another cavity, 2026-10-18 entry above).
- P2 converges at the same rate on both (≈ 50/N pp nominal; 0.02 pp at ~2,500
  elements): the capsule's piecewise-constant sampling, not element shape, sets it.
  The tight design's P2 is within 0.01 pp at every compact size.
- Odd moments stay at round-off: the icosphere is centrally symmetric. (Table above: a
  vertex on +z. Since it became the `hohlraum()` default, 2026-10-22, the icosahedron
  sits in its standard frame, mirror-symmetric in x, y and z, because the 5-fold axis on
  z leaked odd-m modes of ~6e-5 of the drive into the D2-symmetric NIF hohlraum; see
  "Default capsule mesh and visibility".)
- Element quality, not count, sets P4: the lat-long cells near the poles are slivers
  whose centroids sample P4 poorly; at 0.09 mm the "standard" lat-long capsule is still
  0.08 pp off in P4, the compact one 0.005.
- Build time is set by the wall at these sizes (~9 s up to ~3k capsule elements).

## Partial visibility (2026-10-22) *derived*

`view_factors(...; visibility = :sampled, visibility_grid = g)` (also through
`RadiationModel(t; visibility = :sampled)` and a `Scenario`'s `viewfactor_kw`). Each
element's area is lumped into a g × g grid of its tangent plane (area-weighted points);
up to four of those points (the extremes along ±e₁, ±e₂) are probes. Per facing pair,
after the centroid ray: probe rays from each element's probes to the other's (lifted)
centroid. All agreeing with the centroid ray → the `:centroid` factor, bit for bit;
otherwise the g² × g² sample rays give a visible fraction weighted by
cosθᵢ cosθⱼ / r² that scales the unobstructed factor (computed even when the centroid
ray is blocked). One fraction per unordered pair: reciprocity exact. Rays use the facets
and the curved-surface skip rule of the centroid ray. Default `:centroid` checked equal
(`==`) to main's `viewfactors.jl` on `hohlraum(h = 1, capsule_h = 0.4)` and the bored
cavity at 1 mm.

Geometry: the primitive bored cavity of the CAD tutorial (R 2.875, L 9.4, 0.1 mm LEH
bores of radius 1.55 with `nphi ÷ 2` facets, 20 µm waist gap, capsule 0.9, apertures on
the LEHs and the gap exit), `closed = false`, 8 threads. Largest |row sum − 1| (rms):

| component | 0.5 mm, centroid | 0.5 mm, g = 4 | 0.25 mm, centroid | 0.25 mm, g = 4 | 0.25 mm, g = 6 |
|---|---|---|---|---|---|
| walls | 9.7e-3 (4.1e-3) | 1.9e-3 (5.0e-4) | 9.1e-3 (3.3e-3) | 1.4e-3 (6.9e-4) | 1.2e-3 (6.0e-4) |
| caps | 1.7e-3 (1.2e-3) | 1.4e-3 (7.8e-4) | 2.2e-3 (1.2e-3) | 1.2e-3 (5.4e-4) | 9.7e-4 (4.5e-4) |
| bores | 5.7e-3 | 2.6e-3 | 5.2e-3 | 2.1e-3 | 6.3e-4 |
| capsule | 1.6e-3 (9.2e-4) | 5.0e-4 (2.4e-4) | 2.0e-3 (1.1e-3) | 2.0e-3 (6.6e-4) | 2.0e-3 (6.6e-4) |
| LEH apertures | 5.9e-2 (3.2e-2) | 6.7e-3 (2.8e-3) | 6.1e-2 (3.5e-2) | 2.0e-2 (9.3e-3) | 1.5e-2 (7.4e-3) |
| gap faces | 5.0e-2 | 5.5e-2 | 5.2e-2 | 5.5e-2 | 3.3e-2 |
| gap exit | 2.3e-2 | 2.9e-2 | 3.9e-3 | 6.5e-3 | 6.1e-3 |
| time | 1.9 s | 4.0 s (g = 2: 2.5 s) | 4.5 s | 18.5 s (g = 3: 13.4 s) | 50 s |

(2,244 and 6,304 elements.) The 50 µm gap faces behind a 20 µm slot are not helped:
the slot is narrower than a wall sample cell. `hohlraum()` (3,264 elements, flush LEHs):
1.3 s → 4.6 s; at h = 0.5 its walls' 0.30% goes to 0.28% (most changed pairs there are
capsule shadows).
At 1 mm with 12-facet bores the option makes LEH rows worse (3.9% → 7.7%): the faceted
bore (sagitta 0.05 mm on a 0.1 mm bore) over-shadows rays from the LEH rim, which the
fractions then reproduce faithfully. With 64-facet bores at 0.5 mm the residual on the
walls (4.5e-3, near the waist) is the 0.3 mm capsule's faceting.

Reference: visible fractions of 56 changed pairs at 0.5 mm (the 30 of largest exchange
change and 26 at random), against kernel-weighted visibility over 256 × 256 points per
facet pair (same rays, no grid): centroid error mean 0.42, max 0.49; sampled (g = 4)
mean 0.017, max 0.08. At 1 mm, one LEH element against every wall element: 0.020 / 0.16.
Half-shadowed pair (unit squares 2 apart, one turned 45°, a plate over half the
mid-plane; exactly ½ by point symmetry): 0.493 at g = 4, 0.500 at g = 8 (g = 3: 0.430,
its center samples tie on the edge).

Effect after Sinkhorn balancing (the tutorial's two-ring wall source): capsule flux
+0.05%, a₂/a₀ −0.175 pp (0.5 mm) and −0.155 pp (0.25 mm), a₄/a₀ +0.11 pp at both, LEH
loss −0.4%. The shift barely shrinks with the mesh, so balancing does not absorb the
centroid rule's visibility error; which limit is right needs a finer mesh or a
ray-traced radiosity reference (done: next entry, the sampled one is closer). **M**

## Monte Carlo radiosity reference for partial visibility (2026-10-22) *derived*

> **Superseded in part (2026-10-25).** The comparisons stand, but the last bullet's
> reading does not: the flat facets of both geometries are cracked (24-sided LEH disks
> in 20-sided bores and on caps of `nphi` sides, a 180-sided gap aperture on 40-sided
> gap faces; `open_fraction` 7.5e-4 and 2.4e-4), so rays escape that the curved geometry
> keeps, and balancing them adds the exchange the cracks lose. On watertight facets and
> on the curved path balancing does not bias a₂/a₀ upward ("Closure on bored, windowed
> and CAD cavities (2026-10-25)").

`Briareus.montecarlo_incident` (src/montecarlo.jl, unexported): diffuse transport on the
facets with no view factors. Paths start at a uniform point on a capsule element's
facets and go backward (reciprocity of the Lambertian kernel): cosine-weighted
direction, first hit by `closest_hit` over all triangles (apertures included: a hit
there is an escape), the hit element's source `ηP/A` added with the path weight, weight
× albedo, Russian roulette below 1/20. `reemit = :hit` continues from the hit point
(the continuum); `:element` from a uniform point on the hit element (the element model
with exact view factors and visibility). Statistical errors from 64 batches. Checked
(test/montecarlo.jl): a closed cube with uniform source gives S/(1 − α) within 1%; an
open-top cube against the closed-form factors (Hottel 0.199825 opposite,
(1 − 0.199825)/4 adjacent) within 5σ and 1%; 8 × 8 elements per face against `solve`
within 1.5% on face means (one element per face: the continuum differs by up to 26%).

Since the Monte Carlo traces facets, the comparison uses flat elements (the same meshes
with `surface = nothing`); the curved models are a different geometry (true surfaces)
and are listed only for the size of the visibility shift. Source: the CAD tutorial's
two Gaussian wall rings, integrated over each element, η = 0.75, wall α = 0.8, capsule
α = 0.3. 200k paths per capsule element (105–310 s on 16 threads).
`examples/mc_radiosity.jl`. Capsule incident flux, a₀ (TW/mm²), a₂/a₀ and a₄/a₀ (%):

| case | a₀ | a₂/a₀ | a₄/a₀ |
|---|---|---|---|
| **bored 0.5 mm (2,244 el.)**: MC `:element` | 0.85004 ± 0.00004 | 10.641 ± 0.012 | −0.510 ± 0.015 |
| MC `:hit` (continuum) | 0.85009 ± 0.00005 | 10.667 ± 0.011 | −0.530 ± 0.016 |
| flat, centroid, raw | +0.45% | +0.025 | −0.063 |
| flat, centroid, balanced | +0.08% | +0.508 | −0.274 |
| flat, sampled, raw | +0.29% | +0.108 | −0.063 |
| flat, sampled, balanced | +0.12% | +0.312 | −0.130 |
| **bored 0.25 mm (6,304 el.)**: MC `:element` | 0.84811 ± 0.00003 | 10.752 ± 0.008 | −0.600 ± 0.011 |
| MC `:hit` | 0.84835 ± 0.00003 | 10.746 ± 0.007 | −0.587 ± 0.011 |
| flat, centroid, raw | +0.02% | +0.070 | +0.019 |
| flat, centroid, balanced | +0.16% | +0.303 | −0.107 |
| flat, sampled, raw | +0.11% | +0.033 | +0.016 |
| flat, sampled, balanced | +0.21% | +0.142 | +0.007 |
| **`hohlraum()` (3,264 el.)**: MC `:element` | 0.84677 ± 0.00005 | 10.228 ± 0.012 | −0.465 ± 0.014 |
| MC `:hit` | 0.84673 ± 0.00004 | 10.238 ± 0.010 | −0.458 ± 0.014 |
| flat, centroid, raw | +0.46% | −0.091 | +0.016 |
| flat, centroid, balanced | +0.22% | +0.086 | +0.013 |
| flat, sampled, raw | +0.18% | +0.019 | +0.011 |
| flat, sampled, balanced | +0.23% | +0.083 | +0.018 |

(Model rows: a₀ relative to MC `:element`, moments as pp differences.) Curved models,
balanced, sampled − centroid: a₂/a₀ −0.175 / −0.155 pp, a₄/a₀ +0.11 / +0.11 pp (bored
0.5 / 0.25 mm), −0.003 / +0.004 pp on `hohlraum()`; flat: −0.196 / −0.161 and
+0.144 / +0.114 pp, −0.003 / −0.005 pp. The shift is the same on both geometries.

- The two Monte Carlo modes agree within 0.03 pp (2σ): the element discretization of the
  radiosity is negligible here, and MC `:element` is the reference for the view factors.
- After balancing, sampled visibility is closer in every bored case: a₂/a₀ error
  0.51 → 0.31 pp (0.5 mm) and 0.30 → 0.14 pp (0.25 mm), a₄/a₀ 0.27 → 0.13 and
  0.11 → 0.01 pp. The −0.16 pp shift is a correction, not an artifact; it is about
  half the balanced centroid model's error. On the flush-LEH control both visibilities
  balance to the same answer (+0.085 pp), as expected.
- Balancing itself biases a₂/a₀ upward: +0.06 to +0.48 pp over the raw factors of the
  same visibility, on every case including the control. Raw sampled factors are the
  closest of the four in a₂/a₀ on the 0.25 mm cavity and the control (0.02–0.03 pp);
  their a₀ is off 0.1–0.3%, balanced 0.1–0.25%. Presumably Sinkhorn spreads the
  residual row error (on the cavity, the waist gap's faces at 5% and the LEH apertures
  at 1–2%) over each row's partners symmetrically, which is not where the missing
  exchange is (not checked element by element). This holds for flat
  elements; the curved (default) path's row errors differ and a curved-surface Monte
  Carlo would be needed to check balancing there (done: next entry; on the curved path
  balancing moves toward the continuum). **M**

## Balancing on the curved path against a continuum Monte Carlo (2026-10-22) *derived*

Does `closed = true` (the `hohlraum()` default) bias the capsule drive on the default
curved path? Reference: `montecarlo_incident(...; reemit = :hit)` (the continuum) on a
flat copy of the hohlraum faceted at wall/cap h = 0.1 mm, capsule_h = 0.05 mm (33,292
elements, 6,612 on the capsule; the chord sagitta of the wall is 4e-4 mm). Geometry
10.1 × 5.75 mm, 3.1 mm LEHs, 0.9 mm capsule; wall gold α = 0.925 and 0.85, η = 0.8.
Source: 48 NIF quads at 8 TW at the published cone aims (`examples/mesh_convergence.jl`
"nominal"), deposited on each mesh by its own `DepositionOperator` (totals agree to
4e-5 between meshes) and scaled to the reference total. Departure: one outer-44 quad
off; by linearity the Monte Carlo traces that quad's source alone (paired batches give
its error). 298M paths at 0.925 (8.6 min on 16 threads), 198M at 0.85.
`examples/mc_closure.jl`.

Reference convergence in faceting (α = 0.925, 32M paths each): wall h 0.2 / 0.1 / 0.06
gives a₀ 8.4148 / 8.4116 / 8.4118 (± 0.0011) and a₂/a₀ −12.547 / −12.546 / −12.554
(± 0.027) %; capsule_h 0.035 instead of 0.05 moves a₀ +0.03 ± 0.02% and a₂/a₀
+0.04 ± 0.04 pp (a₄/a₀ +0.08 ± 0.05 pp, mostly the piecewise-constant moment analysis,
which `capsule_h` sets for model and reference alike). Converged to ~0.03% in a₀ and
~0.03 pp in a₂/a₀.

Capsule incident flux: reference a₀ (TW/mm²), a₂/a₀, a₄/a₀ (%), and the departure
(design − quad off: Δa₀ in %, Δa₂/a₀ and Δa₄/a₀ in pp); model rows are differences from
the reference (a₀ in %, the rest in pp). Default visibility (centroid). Curved = the
`hohlraum()` default; flat = the same meshes without surfaces, a coarser polygonal
geometry than the reference (so its rows mix geometry error in).

| α = 0.925 | a₀ | a₂/a₀ | a₄/a₀ | Δa₀ | Δa₂/a₀ | Δa₄/a₀ |
|---|---|---|---|---|---|---|
| **MC `:hit`** | 8.4118 ± 0.0003 | −12.513 ± 0.010 | +0.225 ± 0.010 | +2.045 ± 0.001 | +0.233 ± 0.001 | −0.035 ± 0.001 |
| curved 0.45/0.15, raw | +0.56% | +0.161 | +0.266 | +0.000 | −0.006 | −0.003 |
| curved 0.45/0.15, balanced | −0.07% | +0.091 | +0.249 | +0.000 | −0.004 | −0.004 |
| curved 0.28/0.09, raw | −0.03% | +0.088 | +0.149 | −0.001 | −0.001 | −0.002 |
| curved 0.28/0.09, balanced | −0.05% | +0.029 | +0.064 | −0.000 | −0.001 | −0.002 |
| curved 0.17/0.05, raw | −0.09% | +0.036 | +0.006 | +0.000 | +0.001 | −0.001 |
| curved 0.17/0.05, balanced | −0.03% | +0.020 | +0.007 | +0.000 | +0.000 | −0.001 |
| flat 0.45/0.15, raw / balanced | +1.12 / +0.54% | +0.030 / +0.038 | +0.075 / +0.063 | +0.001 | +0.000 / +0.001 | +0.003 / +0.002 |
| flat 0.28/0.09, raw / balanced | +0.24 / +0.22% | −0.017 / +0.030 | +0.007 / −0.000 | −0.001 / +0.000 | +0.001 | +0.001 |
| flat 0.17/0.05, raw / balanced | +0.03 / +0.07% | +0.039 / +0.025 | −0.014 / −0.015 | +0.000 | +0.001 | −0.000 |

| α = 0.85 | a₀ | a₂/a₀ | a₄/a₀ | Δa₀ | Δa₂/a₀ | Δa₄/a₀ |
|---|---|---|---|---|---|---|
| **MC `:hit`** | 6.1316 ± 0.0003 | −14.547 ± 0.011 | +0.150 ± 0.014 | +2.013 ± 0.001 | +0.315 ± 0.001 | −0.046 ± 0.002 |
| curved 0.45/0.15, raw | +0.35% | +0.148 | +0.322 | −0.001 | −0.009 | −0.007 |
| curved 0.45/0.15, balanced | −0.06% | +0.077 | +0.306 | −0.001 | −0.008 | −0.007 |
| curved 0.28/0.09, raw | −0.03% | +0.074 | +0.178 | −0.002 | −0.003 | −0.005 |
| curved 0.28/0.09, balanced | −0.04% | +0.016 | +0.095 | −0.001 | −0.003 | −0.005 |
| curved 0.17/0.05, raw | −0.07% | +0.022 | +0.034 | −0.001 | −0.001 | −0.004 |
| curved 0.17/0.05, balanced | −0.02% | +0.010 | +0.035 | −0.001 | −0.001 | −0.004 |

Mean |row sum − 1| before balancing (curved): 7.6e-4 / 4.6e-4 / 2.0e-4, largest
3.2e-3 / 3.1e-3 / 1.5e-3 at the three meshes.

- On the curved path balancing does not bias the drive: it moves a₂/a₀ toward the
  continuum at every mesh and both albedos (0.161 → 0.091, 0.088 → 0.029,
  0.036 → 0.020 pp at 0.925) and fixes the study mesh's a₀ (+0.56% → −0.07%; the raw
  excess is the ~8e-4 row-sum error amplified by the wall albedo). The flat-element
  upward bias of the previous entry does not carry over: raw curved factors already
  sit above the reference in a₂/a₀, and balancing pulls them down.
- Balanced curved models are within 0.07% in a₀ and 0.09 pp in a₂/a₀ at every
  preset, 0.03% and 0.02 pp at the fine one (inside the reference's own ~0.03
  convergence). a₄/a₀ at the coarser presets is dominated by the capsule mesh (the
  moment analysis on 0.15 and 0.09 mm facets, cf. "Mesh convergence"); at 0.17/0.05,
  whose capsule mesh is the reference's, it agrees to 0.01 pp at 0.925.
- Departures (one quad off) agree with the reference within 0.01 pp in Δa₂/a₀ and
  Δa₄/a₀ and 0.002% in Δa₀ for every model, raw or balanced: balancing is immaterial
  for departures here.
- Not covered: bored LEHs, gaps and CAD cavities (larger row errors; the flat bored
  cavity of the previous entry is the only data there), sampled visibility, windows.
  **H** for `hohlraum()` at these meshes and albedos; **M** beyond.

> **Covered since (2026-10-25):** bored LEHs, a waist gap, windows with gaps and a
> remeshed CAD cavity, centroid and sampled visibility, in "Closure on bored, windowed
> and CAD cavities (2026-10-25)".

## Default capsule mesh and visibility (2026-10-22) *decision*

User decision (Open question 9): `hohlraum()` builds its capsule as the geodesic
`icosphere_mesh` (`capsule_mesh = :compact`); `:latlong` reproduces earlier baselines.
`Target` carries a `visibility` rule that `view_factors` uses by default: `:sampled` for
`cad_target` and `visrad_target` (rims), `:centroid` for `hohlraum()` and hand-assembled
targets.

- **Orientation.** With a vertex on +z (the table in "Compact capsule mesh") the
  icosahedron's 5-fold axis lies along the hohlraum axis, and the D2-symmetric NIF
  hohlraum picked up odd-m capsule modes of ~6e-5 of the drive (test/runtests.jl and
  test/windows.jl selection rules, threshold 1e-10). The icosahedron now sits in its
  standard frame (vertices at the cyclic permutations of (0, ±1, ±ϕ)): its 2-fold axes
  and mirror planes are x, y and z, so the selection rules hold to round-off. P2/P4
  convergence rates are a property of element shape and are unchanged in kind; the
  table's compact rows were not rerun in the new frame.
- **GUI presets.** capsule_h 0.3 / 0.16 / 0.09 mm (320 / 980 / 2,880 capsule elements):
  P2 needs ~2,500 capsule elements for 0.02 pp on any capsule mesh (≈ 50/N pp), so the
  fine preset is not 0.12 as first proposed.
- **Effect on a study.** thor aps-dpp-2026 (N210808 + two THOR windows, h 0.45 /
  capsule_h 0.15, closed, tuned NIF beams) moved by ≤ 0.01 pp in every departure and
  drive (Round 8 vs Round 7); a design optimized on the lat-long capsule evaluates at
  δ2 = 0.015% on the compact one.

## Literature review against Kritcher et al. 2024 (2026-10-22)

Kritcher et al., Phys. Plasmas 31, 070502 (2024), doi:10.1063/5.0210904 (the N221204
design review; CC-BY, full text **H**) was read for the modeling gaps it exposes in a
view-factor code, and six focused reviews were run against it. The full reports, with
every number's source and confidence, are in `reference/litreview/`:
`kritcher2024_notes.md` (the paper), `hohlraum_dynamics.md`, `cbet_calibration.md`,
`energetics.md`, `asymmetry_transfer.md`, `n221204_benchmark.md`, `hydro_coupling.md`.
What follows is the digest that sets roadmap items 9–13 and 16–17.

### The paper's own modeling chain **H**
- Symmetry is held to ~1 % time-varying flux asymmetry; the wall "ingresses into the path
  of laser beams"; late-time symmetry comes from a semi-empirical model (Callahan 2018)
  calibrated on years of data, not first-principles LPI.
- The P2 swing (time-varying P2) is the largest unknown, calibrated by in-flight
  radiography at ~200 µm; for N221204 that calibration was missing, and the N230729
  post-shot was matched by reducing the swing by hand.
- HYDRA settings: inline CBET with δn/n clamp 0.002 until the rise, 0.0075 after (PRE 109
  rounds to 0.008 in the PoP); electron flux limiter 0.15 (conduction, not radiation);
  laser multipliers fit to in-flight symmetry, velocity and keyhole shock velocities;
  the laser is degraded 9–13 % to match data (PRE 109, 025204).
- Mode-1 sources named: laser delivery imbalance, capsule and ice thickness, non-axisymmetric
  diagnostic windows, "manufacturing gaps on the waist" (Hurricane PRL 132, 065103).
- An applied drive mode 1 of ~0.7 % reproduces >100 km/s burn-on (≈75 km/s burn-off).

### Hohlraum dynamics (`hohlraum_dynamics.md`)
Two surfaces move at different speeds, and a view-factor model should treat them apart:
- **Au/He bubble tip** (absorbs the inner beams): ballistic at 300–460 µm/ns. 182 ± 11
  µm/ns on a 3× picket CH shot (Ralph PoP 25, 082701, 2018, **H**); 350/400 µm/ns at
  0.3/0.15 mg/cc (Izumi RSI 89, 10K111, 2018, **H**); 321–461 km/s measured against
  LASNEX f = 0.03/0.15 (Meezan PoP 27, 2020, **H**); x_bub ∝ E_picket^0.2 ρ₀^−0.2
  (Izumi PoP 28, 022706, 2021, **H**). Izumi 2021 gives an "opaque iris" clipping
  fraction of the inner cones and an absorbed fraction with bubble transparency 29 %
  (0.3 mg/cc) and 54 % (0.15) at 5.5 ns; Ralph puts inner-beam absorption in the bubble
  at 51–62 % and finds a view factor at the measured peak wall position needs cone
  fraction 18.5 % for zero flux P2.
- **X-ray-emitting wall**: 50–100 µm/ns; 400–600 µm inward by peak, and VisRad models of
  LEH images needed the wall radius reduced by 300–700 µm (Chen PoP 27, 022702, 2020, **H**).
- **LEH**: closes in the low-power foot, reopens at the rise, varies ≤ 4 % through peak,
  about 0.84 r₀ at peak (Schneider PoP 22, 122705, 2015, **H**); closure runs ~100 µm in
  diameter faster than simulated (Swadling "Build-A-Hohlraum", PoP 32, 052707, 2025, **H**).
- **Callahan's metric** (PoP 25, 056305, 2018, **H**): P2 = 0 at cone fraction 33 % when
  √(E_pic,out/(A_out ρ_fill))·(τ/R)(r_cap/R) ≈ 1.2 (kJ, mm², mg/cc, ns, mm); 39 shots,
  σ = 8 µm; 2.7 µm hot-spot P2 per 1 % cone fraction. Izumi restates it as
  P2 − 270(CF − 0.333) = 85(1.2 − FoM) µm.
- **Prior art**: Young et al., PoP 27, 082702 (2020, **H**) is LLNL's "dynamic view factor":
  VisRad at 100 times with as-delivered powers, the cylinder shrinking to the mean of the
  absorption and emission surfaces from HYDRA, LEH closure from SXI; a static model with
  all radii at 85 %, albedo 0.8/0.3 and conversion 0.9 agrees to 0.14 % in mode 1 at peak
  and < 5 % in hot-spot velocity. The research-notes entry crediting IRAD3D with spot
  motion is not supported by its abstracts (**L**).
- A simulated flux P2/P0(t) for N210808/N220919/N221204 exists only as PRE 109, 025204
  Fig. 10 (**L**, by eye): ≈ 0 in the foot, −0.17 at 5.5 ns, 0 at 6.5 ns, +0.12–0.15 at
  8 ns. No tabulated swing is public.

### CBET calibration (`cbet_calibration.md`)
- **Clamps are per-phase effective gains, not limits** (**L**): linear δn/n at NIF
  conditions is 1e-4 to 1e-3, below every ignition-era clamp. Values: Kritcher PRE 98,
  053206 (2018, 0.6 mg/cc CH, **H**): 0.01 in the foot; 5e-4 to 2e-3 in rise/peak gives
  1.36–1.6× inner power, 2–3e-3 gives 1.67–1.9× (all HYDRA). Pickworth PoP 27, 102702
  (2020, **H**): 0.01 throughout plus peak cone-fraction multipliers 0.86 (0 Å) and 0.80
  (1 Å). N221204: 0.002 foot, 0.0075 from the rise (PRE 109, 025204, **H**). Strozzi PRL
  118, 025002 (2017) caps δn with a hard minimum; Briareus's smooth clamp is √2 lower at
  the clamp.
- **N221204 transfer** (PRE 109 Figs. 9–10, by eye, **L**): cone fraction input → post-CBET
  0.13–0.15 → 0.18 (foot), 0.34 → 0.64 (5.5 ns), 0.34 → 0.52 (8 ns); per cone at early
  peak 23° ×2.0, 30° ×1.9, 44° ×0.45, 50° ×0.76; late peak 44° ×0.62, 50° ×0.84. Symmetry
  by power alone would need the outers cut to 150–200 TW (**H**).
- **Hot-spot P2 per Å** (**H**): 40–50 µm (HYBRID-E, 0.3 mg/cc), 24 ± 3 µm (CH, 0.6 mg/cc,
  Pickworth), 14 µm (CH, 0.45 mg/cc, Izumi, at 30–40 % lower peak power: transfer scales
  with intensity). Converted at 2.7 µm per % cone fraction: dCF_eff/dΔλ ≈ 0.15–0.18 per Å
  for HYBRID-E, 0.085 (0.6 mg/cc CH), 0.05 (0.45 mg/cc) (**L**).
- **Where it happens**: outside the LEH the supersonic outflow drives inner → outer
  transfer (into the 50° cone), inside it runs outer → inner, so Δλ = 0 is not zero
  transfer; I-Raum needed 0.5 Å to null the flow-driven part (Pickworth HYDRA; Ralph
  Nat. Commun. 2024, **H**). The 44° cone gives up more than the 50°.
- **Plasma state**: no measured low-fill LEH ne, Te, Ti or flow history is public; Michel
  PoP 20, 056308 (2013, high fill, **H**) gives ne 0.03 nc, Te 2.8 keV, Ti 0.8 keV at the
  start of the main pulse with ion heating at 6–7 keV/ns, after which the gain is linear in
  Δλ; Barrios PoP 23, 056307 (2016) measured peak Te 4.2 keV, 0.5 keV above HYDRA.
- **Backscatter** (PRE 106, 025201, **H**): N210808 total 13.4 kJ of 1.92 MJ (0.7 %); 50°
  SBS 4.8 ± 2.4 kJ (9.5 ± 3.4 on N210207); the only damage number found is ~1 J/cm² SBS
  fluence (Pickworth); the criterion that limited N220919 is not public.

### Energetics (`energetics.md`)
- **N221204 energy split** (HYDRA calibrated to data; PRE 109, 025204 App. 3, **H**):
  wall 54.1 % (Au lining 4.4, DU 49.7), LEH 25.2 %, He 4.5 %, capsule 16.2 % of incident
  laser energy; 246 kJ "coupled to the implosion". The 245–258 kJ "capsule absorbed
  energy" quoted elsewhere (PoP 31 Table III; Hurricane PRL 132 Table I, Y/G_cap) is the
  energy left in the unablated shell, not total x-ray absorption. A 0-D estimate with
  capsule albedo 0.2–0.3 on the outer area and Dante's T(t) gives 316–361 kJ, bracketing
  the 16.2 % (≈ 290–330 kJ): calibrate a capsule model to the 16.2 %, not to 246 kJ (**L**).
- **Wall material**: U-rich and cocktail walls lose 0.83–0.88 × gold at 270–310 eV, ratio
  ∝ T^−0.22 t^−0.05 (Rosen 2005; Schein PRL 98, 175003, 2007, **H**); DU is worth 6.5 % of
  laser power at peak (Döppner PRL 115, 055001, 2015, **H**, 0.7 µm Au liner). Implied
  multiplier on gold's 1 − α at HYBRID-E peak ≈ 0.80–0.85 (**L**).
- **Conversion efficiency**: 0.7 t^0.12–0.2 (Lindl 2004; Rosen 2005), 0.85 t^0.2 (Schein),
  plateau 95 % at 10¹⁴ and 80 % at 10¹⁵ W/cm² (Dewald 2008, **M**), 85–90 % in NIF vacuum
  hohlraums (Olson 2012, **M**). The Briareus default 0.75 is a Nova 1 ns number.
- **The only published VisRad energetics calibration**: η = 0.95, wall albedo 0.77 match
  Dante 1/2 peak intensity on a 5.75 mm Au hohlraum at 400 TW (Swadling 2025, **H**); its
  0-D model uses capsule albedo ≈ 0.2 and gives α_wall ≈ 0.84 at 288 eV, with
  (1/P_cap) dP_cap/dα_w ≈ 4.5.
- **Trad**: drive Trad 307 eV (N210808 average) and 313 eV (N220919, N221204)
  (Abu-Shawareb PRL 132, 065102 Table I, **H**); Dante 1 313 ± 2 eV at θ 143°, φ 274°
  (37° through the LEH), Dante 2 308 ± 3 eV (Rubery PRL 132, 065104, **H**); the simulated
  internal Trad is ~13 eV below Dante 1 at peak and 15–20 eV below in the trough (PRE 109
  Fig. 1c, **L**). LEH pair at fixed Trad: 3.64 mm at 465 TW → 295 eV, 3.1 mm at 390 TW →
  300 eV (Zylstra PRE 106, 025202 Fig. 2, **H-img**).

### Flux asymmetry → implosion observables (`asymmetry_transfer.md`)
- **Mode 1 → hot-spot velocity**: ~100 km/s per 1 % of a₁/a₀ (burn-off), linear to about
  0.75 % then saturating (Young 2020: 0.44 % → 37.6 km/s 3D HYDRA vs 44.7 km/s model;
  Kritcher 2024: 0.7 % → 75 km/s burn-off, **H**; ratio **L**). Young's recipe is a
  Green's function G(t) in km/s per % per ns from 59 capsule-only impulse runs: negative in
  the foot, dominated by rise and peak. Source sensitivities (N170827): 5.4 km/s per 10 %
  wall–window albedo mismatch; 0.1 / 1.1–1.3 km/s per µm horizontal / vertical offset;
  1.6 km/s per degree of tilt; a dropped inner quad ≈ 2× an outer one.
- **Capsule thickness** mode 1: 110 (sim) to 140 ± 30 (data) km/s per %, opposite in
  sign to a flux mode 1 (+1 % thickness at a pole ≈ −1 % flux there) (Casey PRL 126,
  025002, 2021, **H**).
- **Piston model** (Hurricane PoP 27, 062704, 2020; PoP 29, 012703, 2022, **H**):
  f = v_HS/v_imp = √(RKE/(c·KE)), c ≈ 0.66; no-α yield ∝ (1 − f²)^3.3; α factor
  exp(−(a − 1)χ^1.2 f²) (Casey PRE 108, L053203, 2023 has an extra 1.2 in the exponent);
  at v_imp 380 km/s, 30/40/50/60 km/s cost 5/10/15/20 % of yield. For general asymmetry
  RKE/(c·KE) = 1 − ρR_WHM/ρR_ave.
- **Burn-on amplification**: Nora PoP 30, 092705 (2023, **H**) fits
  0.324 ln(Y_amp + 2.239) + 0.653; Kritcher 2024 Table II implies 1.27–1.39 (**L**).
- **P2/P4**: Masse PoP 26 (2019, **H**) Green's functions: hot-spot P2 response flips
  sign for drive before ~4.5 ns (early waist-hot → oblate), ρR response keeps one sign;
  sensitivity peaks at the rise. Milovich PoP 23, 032701 (2016, **H**): picket P2 < 7.5 %;
  +6 % P2 per 200 µm outward outer pointing; 2 % P2 per 10 % albedo change. Scott PRL
  110, 075001 (2013, **H**): 10 % flux P4 → ~20 µm a₄, up to 15× yield loss, sign set by
  timing. Ralph Nat. Commun. 15, 2975 (2024, **H**): yield vs hot-spot P2/P0 is Gaussian,
  W = 0.306, optimum +11.8 % (HYBRID-E).

### N221204 benchmark (`n221204_benchmark.md`)
- Input deck (**H** unless marked): hohlraum 11.24 mm × 3.2 mm radius, LEH radius 1.55 mm,
  0.3 mg/cc He; capsule outer radius 1.135 mm (N210808: 1.1285), inner 1.050 mm, ablator
  85.4 µm (7.7 / 15.1 doped / 62.6; Pak PRE 109, 025203 Table I); W 0.54–0.6 at.%; DT ice
  64.5 µm; Δλ 2.75 Å at 1ω (0.917 Å at 3ω); 2.05 MJ, 440 TW; SBS 42 ± 11 kJ (inner, early
  peak), SRS 6 ± 1 kJ, 98 % coupled. Per-cone delivered powers are digitizable from Di
  Nicola et al., HEDP (2024), OSTI 2440791, Fig. 1 (outer 296 TW and inner 146 TW at
  peak, cone fraction 1/3; tabulated by eye in the report, **L**). PRE 109 Fig. 9 gives
  the design cone fraction vs time, input and post-CBET (~0.65 at 5.5 ns, **L**).
- Validation targets (**H**): bang time 9.53 ± 0.07 ns (N210808 9.28 ± 0.07); drive Trad
  313 eV; capsule 16.2 % / 258 kJ; hot-spot P2 0.3 ± 0.1 µm, velocity 36 ± 6 km/s at
  (151°, 54°); keyhole N230110 shock-1 breakout 3.864 ns (equator) / 3.906 ns (pole);
  2DConA N240304 (2.9 Å) shell velocity 355 ± 5 km/s, shell P2 21.7 ± 1.3 µm, stagnation
  9.876 ns; rise cone fraction 0.21 → 0.24 raises shell P2 9.1 → 27.5 µm, 0.21 → 0.16
  leaves it (Ruof et al., PoP 2026, OSTI 3367476).
- Still unpublished: absolute pointing, defocus and CPPs per cone; any 23.5°/30° or
  44.5°/50° wavelength split; window sizes and azimuths; corner/lip geometry; Au liner
  thickness; tabulated powers (data "available on request"). Pointing is covered
  non-publicly by the N230729 SST export in thor (`I_Stag_DT_HyE_S08-BeamPointing.csv`,
  reads with `read_visrad_sst`; see `n221204_benchmark.md`).

### Briareus → CompactLES (`hydro_coupling.md`)
- LLNL's capsule-only step applies "a prescribed inward x-ray flux" on a spherical
  boundary, tuned to shock timing, in-flight velocity and shape, the M-band fraction and
  the DT bang time (Clark PoP 26, 050601, 2019, **H**), with asymmetries as
  F = Σ a_l(t) P_l at 100 µm from the ablation front (Masse 2019, **H**). Briareus's
  T_r(t) + a_l(t) export already has that shape; bang time, an input there, becomes the
  first prediction here.
- Stage-0 closures: P_abl ≈ 3 T_heV^3.5 Mbar, ṁ ≈ 3×10⁵ T³ g cm⁻² s⁻¹ (Lindl 2004 eqs.
  6-2/6-3, **H**, exponents partly reconstructed); carbon P_abl ∝ T^3.3 t^−0.09 to 360 eV
  (Trickey HEDP 44, 100995, 2022, **H**, no prefactor); NIC point design P_abl 140 Mbar,
  HYBRID-E hydro efficiency ~14 %, imploding KE 19–33 kJ (Lindl PoP 33, 050503, 2026, **H**).
- The one fully specified public imposed-drive capsule test: Shestakov et al.,
  UCRL-JC-131429 (1998), Be/D capsule, gray diffusion, Marshak BC E − (2/c)F·n = aT_r⁴ at
  160 eV (**H**). M-band is ~15 % of the spectrum above 1.8 keV at peak and pole-hot with
  its own P2 (+0.53 at the wall) (Dewald PoP 27, 2020, **H**).
- CompactLES (2026-10-22) has no radiation, laser, plasma EOS or ablation model; its
  SESAME/IONMIX readers exist but no solver uses them.

## Moving wall on keyframes (2026-10-23) *derived*

`WallMotion` / `moving_hohlraum` / `Keyframes` (ROADMAP Next 9, first piece). The radiating
wall cylinder shrinks at `speed` from `start`; caps keep their planes, LEHs their radius,
the capsule stays; mesh counts are frozen at R₀ so every time has the same elements.
Between keyframes the exchange matrix Aᵢ Fᵢⱼ is interpolated linearly in time and divided
by the areas at that time (reciprocity kept; rebalanced for closed targets); deposition
operators are built at each time.

**Default speed 0.10 mm/ns** (superseded 2026-10-24 by a saturating wall fitted to Young 2020: "Hohlraum dynamics calibrated to published data" below; engineering default, **M**): Chen 2020 (500 µm in by 5.5 ns
at 0.3 mg/cc, 400 at 0.45: ~0.07–0.09 mm/ns); PRE 109, 025204 Fig. 1a (HYBRID-E,
measured off the figure ±0.15 mm: wall-material edge at the waist 2.80–2.85 mm at 6 ns from
3.2 mm, ~0.1 mm/ns; under the outer spots 2.2–2.5 mm).

**Cost** (NIF hohlraum 10.1 × 5.75 mm, WallLossAlbedo, 48 quads, 50-step history, 16
threads): coarse (0.6/0.3, 1,152 elements) view factors 0.3 s, deposition operator 0.6 s,
history 4.3 s; standard (0.28/0.16, 4,796) 2.1 / 0.7 / 63 s; fine (0.17/0.09, 12,512)
12.4 / 1.2 / 945 s. The history dominates (fresh LUs in the albedo iteration, ~19 s a step
at fine), so keyframes on the solve mesh are cheap: six fine keyframes plus per-step
deposition operators add ~14 % to a fine history; a rebuild at every step ~70 %.

**Accuracy** (standard mesh, 0.1 mm/ns from 0.3 ns, 39 steps 0.25–9.75 ns; max error
against a rebuild at every step):

| keyframes (spacing) | interpolated deposition: P2/P0, P4/P0 (pp) | per-time deposition: P2/P0, P4/P0 (pp), Tr (eV) |
|---|---|---|
| 3 (4.9 ns) | 1.56, 0.40 | 0.086, 0.110, 0.15 |
| 6 (1.9 ns) | 0.24, 0.060 | 0.035, 0.022, 0.03 |
| 11 (0.97 ns) | 0.068, 0.020 | 0.012, 0.007, 0.01 |
| 21 (0.49 ns) | 0.022, 0.008 | 0.009, 0.003, 0.005 |

Interpolating depositions blurs the outer footprints, which slide toward the LEH by
ΔR/tanθ; per-time operators (adopted) remove most of the error. A still wall reproduces
the fixed-target history to 1e-12.

**Effect** (tutorial pointing of `docs/literate/retune.jl`, not tuned): the moving wall
makes capsule P2/P0 more positive through the pulse, −0.091 → +0.044 at 9 ns against the
static wall's −0.091 (shift +0.13; 0.005 at 1 ns, 0.04 at 5 ns, 0.08 at 7 ns), and Tr
+1 to +7 eV (smaller cavity). For scale, PRE 109 Fig. 10's HYDRA swing (wall, bubble and
CBET together) is −0.17 at 5.5 ns to +0.08–0.10 at 8 ns.

## Implosion metrics (2026-10-22)

`src/implosion.jl`, the implosion-metrics tutorial and `test/implosion.jl`. The sources
are summarized in `reference/litreview/asymmetry_transfer.md`. Confidence: **H** =
primary text read, **M** = abstract or secondary, **L** = inference or arithmetic.

- **`HOTSPOT_VELOCITY_PER_PERCENT` = 100 km/s per % of mode 1** (`a₁/a₀`, the norm of
  the l = 1 vector, burn-off, rise and peak). *Engineering default.*
  - Young et al., Phys. Plasmas 27, 082702 (2020) **H**: N170601 0.44 % → 37.6 km/s
    (3D HYDRA), 44.7 km/s (VF + G); linear below ~0.75 %.
  - MacGowan et al., HEDP 40, 100944 (2021) **M**: ~100 km/s per %.
  - Kritcher et al. 2024 **H**: N230729 ~0.7 % → 70–75 km/s burn-off (~105 per %, **L**).
  - Spears et al. 2014 **H**: 2 % → 120 km/s (saturated; the linear rule gives 200).
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
    30/40/50/60, rounded); at χ = 1.5 the 5 % row would give f = 0.09. So the table's
    "near ignition" χ is ≈ 2. **L**
  - Casey 2021's "~22 % no-α yield at 100 km/s / 400 km/s" gives 19 % with
    (1 − f²)^{10/3} at f = 0.25. Small gap, not resolved. **L**
- **Expected yield over errors**: `E[exp(−β f²)]` for Gaussian v in closed form,
  `det(I + 2bΣ)^{-1/2} exp(−b μᵀ(I + 2bΣ)⁻¹μ)`, b = β/v_imp². This uses
  (1 − f²)ᵃ ≈ exp(−a f²) (0.7 % high at f = 0.25). *Derived.* Checked against a Monte
  Carlo of full evaluations (3 %; 5 % with α).
- **Burn amplification** (Nora et al., Phys. Plasmas 30, 092705, 2023, Eq. 4, **H** via
  the lit review): 0.324 ln(Y_amp + 2.239) + 0.653 for Y_amp ≤ 40.
  - It is ~1.8 at Y_amp = 30, against Kritcher 2024's burn-on/off ratios of 1.27–1.39.
    Use it as an upper estimate.
  - The Y_amp 80–500 branch as transcribed (13 ln(Y_amp + 1585) − 95) is 1.43 at 80,
    below the low branch's 1.87 at 40. *Resolved 2026-10-24* ("Young 2020 validation"
    below): Eq. 4 prints rounded coefficients; Fig. 6's (13.105, 1585.225, −95.539) are
    implemented on 80–500, with 40–80 refused.
  - The shape term (Eq. 2) as transcribed does not give the stated 0.8–1.2 range for
    |P2| < 0.4 (0.19 at −0.4). *Resolved 2026-10-24*: that range is Table I's (3D);
    Eq. 2 is the 2D fit. Both are in `shape_amplification`.
- **Thickness nuisance in `infer_errors`** (Casey 2021, **H**):
  - 110 km/s per % of thickness mode 1 (simulated; 140 ± 30 measured), toward the
    thicker side.
  - The prior is 0.2 % per axis by default. *Derived* from ~0.15 µm metrology on a
    78 µm shell (**L**).
- **nToF covariance**: ±14 km/s, 13° θ, 29° φ (Young 2020 **H**), in the measured
  vector's radial/polar/azimuthal frame.
- **Picket P2 limit** 7.5 % (Milovich et al., Phys. Plasmas 23, 032701, 2016, **H**).
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
  - Time-dependent error budgets and multi-shot trending remain (ROADMAP 11).

## Closing entrance holes (2026-10-23)

- **Model** (`hohlraum(; leh_closure, leh_open_radius)`, `aperture(...; closure)`): each
  entrance hole is split into the hole that stays open (radius `leh_open_radius`) and a
  lip, the annulus out to `leh_radius`. The lip is an aperture with a `GapClosure`
  albedo: transparent to laser light, lines of sight and view-factor rays at all times,
  it lets all of its incident x-ray flux escape until `closure`, then across `width`
  re-emits α = s·α_after (s the smooth step) and lets 1 − s escape. `energy_balance`
  books τP as escaped and (1 − α − τ)P as absorbed; `drive_history` resolves the
  closure with `transition` steps and counts the lip with the entrance holes (not in
  `window_loss`). `ApertureKeepOut` keeps the full `leh_radius` rim (beams must still
  clear the hardware). It converts no laser light. *Engineering model.* The plasma that
  closes the hole is not computed; closure time and radius are inputs. The Monte Carlo
  reference (`montecarlo_incident`) re-emits from a lip with its albedo at the time
  given (`albedo = albedos(t, time)`); test/montecarlo.jl checks it against `solve`.
- **Default `leh_open_radius` = 0.84 r₀.** Schneider et al., PoP 22, 122705 (2015): the
  LEH closes in the low-power part of the pulse and reopens at the rise. Its radius at
  peak is ~0.84 r₀ and varies ≤ 4 % through peak. The measured open area is
  1.25–1.4× the modeled (HYDRA) area. **H** (via the 2026-10-22 lit review).
  Young et al., PoP 27, 082702 (2020): a static view-factor model with all radii at 85 %
  reproduces the dynamic model at peak. **H**. 0.84 sits between them. **M**
- **Timing.** Chen et al., PoP 27, 022702 (2020): at 0.15 mg/cc a ~1550 µm-radius LEH stays
  open until ~4.5 ns, then its x-ray-emitting size drops rapidly over 4.5–5 ns. Where the
  laser passes, the beams keep it open. **M** (secondary transcription). Swadling et al.,
  PoP 32, 052707 (2025): closure runs faster than simulated (~100 µm in diameter). **H**
  Briareus has no default `closure`. The example below closes at 5.0 ns over 0.5 ns (the
  rise from the foot) with the wall's `WallLossAlbedo`. **L**: an engineering default
  to be set per design from SXI data.
- **What it changes** (*derived*, 2026-10-23). Hohlraum and pointing of
  `docs/literate/retune.jl` (10.1 × 5.75 mm, 3.1 mm LEH, `WallLossAlbedo` wall,
  h = 0.45, capsule_h = 0.15, its reference pulses, steps of 0.4 ns,
  `closed = true`). Compare open holes with `leh_closure = GapClosure(closure = 5.0,
  width = 0.5, after = WallLossAlbedo())` at 0.84 r₀, averaged over 7.6–9.4 ns:
  - capsule Tr goes from 310.2 to 320.4 eV (+10.2 eV, +3.3 %);
  - P2/P0 goes from −9.19 % to −6.67 % (+2.5 points; closing the holes adds pole drive);
  - P4/P0 goes from −0.09 % to −0.15 %.
  Before 4.8 ns the two agree to < 0.01 eV and 0.003 points. A 16 % smaller LEH radius
  is therefore worth ~2.5 points of P2, which is comparable to the cone-fraction
  retunes of the tutorials. Set the closure from data before trusting absolute P2.

## Inner-beam bubble absorption (2026-10-23)

`BubbleAbsorption` (`src/bubble.jl`, ROADMAP Next 9, third piece; `test/bubble.jl`). A
laser–plasma model (`apply_lpi`), so it enters `drive_history(...; lpi)`,
`Scenario(...; lpi)` and composes with CBET (`lpi = (cbet, bubble)`). Sources are in
`reference/litreview/hohlraum_dynamics.md` §1 and "What Briareus should implement" §3.

- **Model.** Tip x_bub(t) = v_b·max(0, t − t_L), t_L = outer turn-on + 0.2 ns,
  v_b = v_ref (E/E_ref)^0.2 (ρ₀/0.3)^−0.2. In the plane z_b where the `footprint`
  (default 50°) beams' axes reach the wall (the moved wall when given a `WallMotion`),
  the bubble is the annulus r > R₀ − x_bub. Each inner beam loses f·(1 − T_b), with f its
  power crossing z_b inside the annulus, from the beam itself: best-focus spot (core and
  halo, 12 equal-power radii × 24 angles) carried along 13 f-cone directions, as the
  deposition model does. T_b = exp(−τ_b (ρ₀/0.3)^q (x/1 mm)^0.5 (R₀/(R₀ − x))²).
  *Engineering model.*
- **Departures from the lit-review spec, and why.**
  - A density factor (ρ₀/0.3)^q in T_b. Without it the 0.15 mg/cc bubble, being deeper,
    is *more* opaque, so no τ_b can meet both of Izumi's anchors. Two anchors fix τ_b and
    q: **τ_b = 0.1441, q = 1.862** (*derived*, `_bubble_anchor_fit`: R = 2.88 mm, launch
    0.5 ns, x = 1.75 / 2.01 mm at 5.5 ns). q ≈ 2 is the ρ² of inverse bremsstrahlung,
    which supports the form (**L**).
  - The absorbed power is not deposited on a separate wall band: it is moved to the
    footprint beams of the inner beam's hemisphere, in equal shares, so the existing
    deposition operator puts it where those beams land. Beam power is conserved
    exactly (per hemisphere), the band follows the outer pointing and the moving wall,
    and no new operator or keyframe data is needed. Equal shares (not power-weighted)
    keep it defined when the outer power is zero. The footprint beams' own LEH loss
    (a few 1e-4) applies to the moved power.
  - The bubble edge is a logistic step 0.05 mm wide (`edge`), so f is smooth in time,
    depth and pointing (ForwardDiff through `speed`, `opacity` and the beams).
  - The tip is measured from the initial wall; only the footprint plane follows a
    moving wall.
- **Defaults** (superseded 2026-10-24: speed, launch, opacity, q and collision refitted, "Hohlraum dynamics calibrated to published data" below; engineering, **M**): v_ref 0.35 mm/ns at 0.3 mg/cc (Izumi PoP 28,
  022706, 2021, 345 µm/ns; Meezan 2020 measured 321–461 µm/ns); scalings E^0.2, ρ₀^−0.2
  (Izumi 2021, ±0.03 each); delay 0.2 ns (Izumi: intercept 0.5 ± 0.1 ns for a 0.3 ns
  turn-on); p = 0.5 (Izumi Eq. 6); (R/(R − x))² stands for Izumi's ρ_bub ∝ (R/R_bub)^k
  with k = 1 (the paper's k is illegible in the extraction, **L**).
- **Checks** (*derived*):
  - Anchors: T_b = 0.290 (0.3 mg/cc) and 0.540 (0.15) at 5.5 ns, exact by construction.
  - f against a closed form: a Gaussian spot on a collimated vertical beam gives
    exp(−ρ²/a²) to 3e-3 (nr = 400).
  - Ralph et al. 2018 (6.72 mm CH high-foot, 0.6 mg/cc). Geometry **L**: R = 3.36,
    L = 11.3 mm, the retune tutorial's pointing scaled in z by 11.3/10.1. On Ralph's
    measured trajectory (0.182 mm/ns, at the wall at 2.47 ns), inner absorption reaches
    Ralph's 51 % and 62 % at x = 1.93 and 2.05 mm, i.e. 13.1 and 13.7 ns. That is
    consistent with his anchor if mid-peak falls there (the high-foot pulse timing is
    not in our notes, **L**). The 30° cone first meets the bubble (1 % clipped) at
    r = 2.11 mm, against Ralph's 2.4 mm; the scaled pointing is a guess.
  - Swing direction. Retune-tutorial hohlraum (10.1 × 5.75 mm, `WallLossAlbedo`,
    h 0.45 / capsule_h 0.15, 0.25 ns steps), wall 0.1 mm/ns from 0.3 ns (11 keyframes),
    lips `GapClosure(closure = 5.0, width = 0.5, after = WallLossAlbedo())`, bubble
    turn-on 0 (launch 0.2 ns), footprint following the wall. Capsule flux P2/P0:

    | t (ns) | static | wall | wall + lips | wall + lips + bubble | PRE 109 Fig. 10 (N221204) |
    |---|---|---|---|---|---|
    | 2.0 | −0.173 | −0.167 | −0.167 | −0.167 | ≈0 |
    | 4.0 | −0.151 | −0.129 | −0.129 | −0.107 | −0.013 |
    | 5.5 | −0.122 | −0.075 | −0.048 | +0.051 | −0.168 |
    | 6.5 | −0.099 | −0.029 | −0.004 | +0.131 | −0.057 |
    | 8.0 | −0.092 | +0.012 | +0.035 | +0.172 | +0.090 |

    The bubble moves P2/P0 positive, as HYDRA's swing does, by +0.10 at 5.5 ns and
    +0.14 from 6.5 ns. Inner absorption is 54 % (23°) and 44 % (30°) at 5.5 ns, 88/79 %
    at 6.5 ns, and > 99 % by 8 ns. Capsule Tr falls 2–8 eV against wall + lips (the moved
    power lands near the LEH). The tutorial pointing is not HYBRID-E's and its foot
    is not tuned (−0.17 against ≈0), so only the direction and the late rise compare.
- **Stagnation and ridge stop (added 2026-10-23, same day).** The ballistic tip crossed a
  2.9 mm hohlraum by ~8 ns and clipped every inner beam (> 99 % absorbed), so late-peak
  P2 plateaued at full inner absorption. Two options, on by default (*engineering*, **M**):
  - *Stagnation depth.* `bubble_depth` is ballistic up to the knee x_max − w, then
    x_max − w·exp(−(x_bal − x_max + w)/w) (C¹, stops at x_max). Defaults x_max = 0.80 R,
    w = 0.08 R: exactly ballistic to 0.72 R (2.07 mm at R = 2.88), past both anchor
    depths (1.75 / 2.01 mm at 5.5 ns), so **the anchors stay exact** (T_b 29 % / 54 %,
    tested bit for bit against the ballistic tip), and near the ~2.1 mm the ballistic fit to
    Izumi's 0.3 mg/cc trajectories reaches by 6.5 ns; it stops between that and the
    ~2.4 mm (0.83 R) of 0.15 mg/cc (end of the measured range, 1.2–2.4 mm).
    Izumi 2023: the tip slows on back-pressure near the axis and is not ballistic in
    low-contrast pulses (stall in the trough, relaunch ≤ 300 µm/ns at ~6 ns; not
    modeled, set `speed`/`turn_on` for such pulses). *Not* x_max = R − 0.73 R: the ridge
    radius is measured at the waist, and 0.78 mm would stop the tip before it reaches
    the inner beams at z_b (1.3 mm from the wall, Izumi 2021) and contradict the
    measured ballistic trajectories (1.2–2.4 mm over 4.5–6.5 ns). Izumi 2018's point that
    attenuation follows the line-integrated density, not the displacement, is what the
    (R/(R − x))² density factor stands for; it is not modeled further.
  - *Ridge stop.* After t_c a fraction g(t) (C¹ smoothstep over t_c ± 0.2 ns) of each
    inner beam's power left by the bubble is deposited on an axisymmetric band of the
    wall component centred at the height where the beam's axis crosses r = 0.73 R
    (Izumi 2023: ridge at 2.1 mm in 2.88 mm), Gaussian in z with rms 0.5 mm, area-weighted.
    The heated ridge is taken to radiate at the wall's conversion from the wall just
    outside it (Briareus has no volume emitter; a lower ridge conversion would lower
    the waist drive further). Power is conserved and booked on the beam's `per_beam`.
    t_c = turn_on + 5.5 − 2.0 (r_cap/r_basecamp − 1) − 0.7 (E − 1) ns (*derived*, Izumi
    2018: spots separate at 5.8 ns for a 0.3 ns turn-on, ~0.4 ns earlier for 1.2× r_cap,
    5.45 ns with a 1.5× picket). The 0.5 mm band width and 0.2 ns rise are engineering
    choices (**L**). It is applied by `lpi_deposit` (new, `lpi.jl`), which
    `drive_history` (fixed and keyframed) and `Scenario` evaluations now call instead
    of `deposit`. Since the same day the beam-error sensitivities (`uncertainty.jl`:
    per-beam dual partials take g of an inner beam off the operator and put it on its
    ridge band) and `influence_matrix` see it too (checked against central
    differences of `lpi_deposit` in test/uncertainty.jl).
  - `stagnation` and `collision` are generic (`remake(m; speed, opacity, stagnation,
    collision)`), for the Next 10 fits; `stagnation = Inf`, `collision = Inf` restore
    the ballistic, ridge-free model.
  - Swing rerun (same retune hohlraum, wall + lips, bubble turn-on 0 so t_c = 5.5 ns):

    | t (ns) | wall + lips | + ballistic bubble (before) | + stagnation only | + stagnation + ridge (defaults) | PRE 109 Fig. 10 |
    |---|---|---|---|---|---|
    | 2.0 | −0.167 | −0.167 | −0.167 | −0.167 | ≈0 |
    | 4.0 | −0.129 | −0.107 | −0.107 | −0.107 | −0.013 |
    | 5.5 | −0.048 | +0.051 | +0.051 | +0.057 | −0.168 |
    | 6.5 | −0.004 | +0.131 | +0.126 | +0.128 | −0.057 |
    | 8.0 | +0.035 | +0.172 | +0.152 | +0.150 | +0.090 |
    | 9.0 | +0.065 | +0.190 | +0.169 | +0.163 | |

    Inner absorption now plateaus at ~90 % (23°) and ~81 % (30°) from ~7.5 ns (tip
    2.29 mm at 8 ns) instead of > 99 %; Ralph 2018 measured 51–62 % at peak. The late
    rise from 6.5 to 9 ns falls from +0.059 to +0.035 (the wall and lips alone give
    +0.069). The ridge stop moves P2/P0 by < 0.01: its band sits ~1–2 mm from where the
    inner beams would hit the wall, both near the waist, and converts at the wall's
    efficiency. The early offset (+0.10 above Fig. 10 at 5.5 ns, where HYDRA is at its
    minimum) is unchanged: the bubble absorbs 54 / 44 % of the inner power by 5.5 ns.
    Defaults were not tuned to Fig. 10, and the pointing is not HYBRID-E's.
- **Open** (addressed 2026-10-24, "Hohlraum dynamics calibrated to published data": Ralph's 51–62 % is a peak average; N221204's is 0.575). The plateau absorption (~80–90 %) sits above Ralph's 51–62 %: with Izumi's
  anchors fixed, any stagnation past the anchor depths leaves T_b ≲ 1 % on clipped
  light. Fitting v_ref, τ_b, x_max (and t_c) to shot data (Next 10) is the next step.
  Low-contrast relaunch (Izumi 2023) is not modeled.

## Mock HYBRID-E pointing (2026-10-23) *derived*

`hybride_mock_pointing()` (`src/facilities.jl`, table `HYBRIDE_MOCK_POINTING`) is a
public stand-in for the HYBRID-E pointing of the N221204 benchmark (ROADMAP 16). The
user's brief (2026-10-22): the same symmetry, the same static capsule radiosity and
similar CBET interactions as the shot's pointing, with the non-public N230729 SST export
(thor, `read_visrad_sst`) used only to fit and check it. Its aim points are not in this
repository; only the aggregate agreement below is.
- Form: one rule per cone in each quad's own frame, mirrored top to bottom, defocus 0:
  aim z, aim radius r along the quad's port azimuth (inner cones across the axis) and,
  for the outer cones, a tangential split ±s with each beam aiming to the far side of
  its quad (the 450 µm outer-quad split of Kritcher et al. 2021 is this kind). The
  export has the same form, so the mock keeps its symmetry class: no top/bottom
  difference beyond the port pattern, and each cone's quads rotated copies of each
  other.
- Fit: the seven numbers (four z, two inner r, one outer s) on a 0.1 mm grid, exhaustive
  search of the per-cone grid combinations against the export's capsule-flux Legendre
  P1–P6 and real Y_lm (l ≤ 6) over P0 at a foot (cone fraction 0.15) and a peak (1/3)
  cone fraction, keeping only candidates whose LEH-plane beam overlaps agree within
  10 % and ranking ties by clipped power. Several grid points fit about equally well;
  the one chosen is not the export rounded to the grid.
- Geometry (`n221204_benchmark.md` §5): 11.24 × 6.4 mm hohlraum, 3.1 mm LEHs, 1.135 mm
  capsule. Corner radius, windows and the Au/DU liner are unpublished: repository
  defaults (square corners, no windows, `gold()` wall, `ablator()` capsule). All 192
  beams (`nif_beams()`, f/20.7, CPP spots), 440 TW, inner power split evenly over the
  23.5°/30° cones and outer over 44.5°/50°, one static solve per cone fraction.
- Agreement, mock − export, fine mesh (h 0.17, capsule_h 0.09; the standard mesh 0.3/0.16
  agrees with it to 0.01 pp): P2/P0 +0.14 pp (foot) and −0.03 pp (peak); P4/P0 and
  P6/P0 within 0.03 pp; odd P_l zero for both; every Y_lm/Y00 with l ≤ 6 within 0.06 pp
  (foot) and 0.02 pp (peak), and the largest azimuthal (m ≠ 0) term the same to
  0.001 pp; capsule-mean flux within 0.1 %; clipped power ~1e-5 of the laser for both.
- CBET: the reduced model (`CBET`, leh_z 5.62 mm, Δλ 0.917 Å inner at 3ω, clamps 0.0075,
  0.002 and none) transfers the same power to within 0.03 % of the transfer for every
  cone. Its pointing dependence is only through each beam's width at the LEH plane,
  which an aim near the LEH and f/20 optics leave insensitive to 0.1 mm moves; so the
  check also compares a geometric one, the mean overlap of Gaussian beam footprints
  (width from `_beam_at_plane`, axes where they cross the LEH plane) per cone pair:
  within 7.8 % (worst pair, outer44–outer44), 3.5 % rms over the ten pairs.
- Use: a stand-in for pointing-level studies of the benchmark (moments, CBET, LEH
  clearance), not a facility pointing. Fit and check scripts: thor `hybride_mock/`.

## N221204 time-dependent benchmark (2026-10-23) *derived*

`examples/n221204.jl` (ROADMAP Next 9 and 16): capsule flux P2/P0(t) of an N221204-class
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
  192 beams at `hybride_mock_pointing()` (a fitted public stand-in, not the shot's
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
  post-transfer one (0.65 at 5.6 ns, 0.55 at 7, 0.52 at 8; HYDRA's own CBET result for
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
  302.1 at 8 ns (the moved inner power lands near the LEH); (a) 288.7 → 302.1 (lips) →
  299.3 eV. Without CBET 293.0 eV. The +lips value matching 313 eV is partly the
  undegraded laser (HYDRA degrades it 9–13 %) against the missing windows and lip
  hardware, so it is not a calibration.
- **Diagnosis.**
  1. *The reduced CBET reverses the transfer.* At the deck's settings it moves power
     from the inner to the outer cones: post-transfer cone fraction 0.14 at peak (input
     0.33) against Fig. 9b's 0.65, with the 50° cone gaining ×2.1 and the 23° cone
     keeping 29 %. A static probe at peak (6 ns): Δλ = 0 gives 0.013, Δλ = 0.917 Å 0.14,
     clamp 0.002 0.15, Δλ = 2.0 Å (3ω) 0.66, and `mach = 0` with 0.917 Å 0.80. The
     default Mach-1 axial outflow Doppler-shifts the outer beams to the receiving side
     and overwhelms the 0.917 Å detuning. With (a) the static P2 is near zero and every
     dynamic piece pushes it positive, so (a) cannot show the swing; the comparison is
     informative only with (b). This is the CBET item of ROADMAP Next 10: two crossing
     zones per LEH (outer → inner inside), plasma-state (flow, Te, ne) inputs and a
     multi-epoch `calibrate_cbet` against the published post-CBET per-cone powers
     (Fig. 9a/b, Fig. 1d) would fix it (done: "Two-zone CBET calibrated to PRE 109
     Fig. 9" below; the example now defaults to it).
  2. *With Fig. 9b's transfer the swing is reproduced in shape and timing of the
     minimum; the late side rises too early and too far.* Contributions at 8 ns: wall
     +0.086, lips +0.032, bubble +0.226; at 6.5 ns: +0.031, +0.035, +0.173. The bubble
     dominates the swing, as Fig. 10's caption says ("as the wall enters and blocks the
     inner beams"), but its inner absorption (48 / 46 % at 6.5 ns, 85 / 83 % at 8 ns for
     23° / 30°) sits above Ralph 2018's 51–62 % at peak, and the wall has moved 0.56 mm
     at 6 ns against 0.35–0.40 mm at the waist in PRE 109 Fig. 1a (constant 0.10 mm/ns
     from the turn-on overstates the motion in the low-power foot). Both put the zero
     crossing 0.4 ns early and the late plateau +0.07 high. Data that would fix it
     (Next 10): the HYBRID-E inner-beam propagation/absorption history (Ralph-type
     inner-beam imaging or HYDRA bubble transparency vs time) for v_ref, τ_b, x_max; the
     emitting-wall radius vs time (SXI/wall-motion imaging) for the wall speed and its
     start; LEH SXI for the lip closure time.
  3. *The foot is offset* by −0.04 to −0.07 with (b) (static −0.07 at 2 ns, −0.11 at
     4 ns) where HYDRA holds ≈0: the foot P2 is set by pointing and foot cone fraction,
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
  - Colours: each pixel's optical density is unmixed into blue, red and black. A 5 px
    opening removes the thin black laser-power and cone-fraction lines.
  - Centreline: in each pixel column, a 12 px stroke (the measured width) is fitted to
    cover the curve's own colour and no background. Curves drawn over it count as
    unknown: red is drawn over black, and black over blue.
  - Each sample averages ±0.027 ns of columns, every 0.1 ns over 0.8–8.9 ns.
  - The two narrow spikes (N220919 at 8.40 ns, N210808 at 7.66 ns) are traced row by
    row along their flanks.
  - Each curve's extrema are added as extra rows. Their values are the stroke's outer
    edge minus half the stroke width; a column centre is biased by up to 0.003 at a
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

## Two-zone CBET calibrated to PRE 109 Fig. 9 (2026-10-23) *derived*

This fixes diagnosis 1 of the N221204 benchmark above. At the deck's settings the
single-zone reduced `CBET` sent the transfer the wrong way: post-CBET cone fraction 0.14
at peak, against HYDRA's 0.65. Its Mach-1 outflow outweighs the 0.917 Å (3ω) detuning.

- **Data.** PRE 109, 025204 Fig. 9 (N221204, HYDRA), digitized by pixel analysis of the
  figure's raster into `data/pre109_fig9_n221204.csv` (`n221204_fig9()`).
  - Coverage: 0.7–9.0 ns in 0.1 ns steps. Columns are the input and post-transfer inner
    cone fraction (Fig. 9b) and the per-cone power (Fig. 9a).
  - Reading error ≈ ±0.01 in cone fraction and ±2 TW per cone. Foot inner-cone powers
    (< 8 TW) are below the figure's resolution.
  - Checks: Fig. 9a's per-cone post powers reproduce Fig. 9b's cone fraction to 0.01 at
    peak. Total power is conserved post/input to < 1 %.
  - Cone fraction, input → post: 0.15 → 0.22 (2.0 ns), 0.12 → 0.20 (3.4 ns), 0.18 → 0.40
    (4.4 ns), 0.34 → 0.65 (5.6 ns), 0.34 → 0.57 (7.0 ns), 0.35 → 0.52 (8.0 ns).
  - Fig. 9a post/input per cone (23° / 30° / 44° / 50°): ×1.92 / ×1.88 / ×0.43 / ×0.64 at
    5.6 ns; ×1.53 / ×1.46 / ×0.64 / ×0.84 at 8 ns.
  - The plotted inputs are HYDRA's degraded powers: ≈370 TW at peak, against 442 TW
    delivered.
- **Model** (`TwoZoneCBET`, `src/cbet.jl`). Two full `CBET` zones per LEH, applied in
  the order the beams meet them.
  - *Outside*: crossing plane 0.5 mm outside the LEH, supersonic outflow. The outflow's
    Doppler shift makes the steeper inner beams bluer in the flow frame, so this zone
    moves power inner → outer even at Δλ = 0.
  - *Inside*: crossing plane 0.5 mm inside the LEH, Mach 0. Here the red-detuned inner
    beams gain.
  - This is P20's HYDRA picture. It is also why Δλ = 0 is not zero transfer at low fill
    (NC24).
  - Flow is signed: negative means into the hohlraum.
  - The single-zone `CBET` and its defaults are unchanged. One HYDRA shot is not enough
    evidence to change a general default, and the two-zone fit is tied to that shot's
    timing.
  - Also unchanged, at repository defaults (L): CH composition, ne/nc 0.06, Ti 1.5 keV.
- **Calibration** (`calibrate_cbet(model, beams_at, data; fit)`, the multi-epoch method).
  Script `examples/n221204_cbet_calibration.jl`: 67 L-BFGS-B iterations, 9 min on
  4 threads.
  - Beams: `hybride_mock_pointing(nif_beams())` with Fig. 9's input powers
    (`n221204_fig9_beams`). Outer cones come from Fig. 9a. The inner total comes from
    Fig. 9b's input cone fraction, split evenly between 23.5° and 30°. Δλ is 0.917 Å
    (3ω) on the inner cones.
  - Both zones use the deck's clamp schedule: 0.002 → 0.0075 across 3.9–4.1 ns.
  - Target: Fig. 9b's post-transfer cone fraction every 0.2 ns over 0.8–8.4 ns
    (39 points), with 10 fitted parameters:

  | Parameter | Value | Plausibility |
  |---|---|---|
  | outside multiplier | 0.30 | Published NIF practice treats gains as calibration factors. A value < 1 also stands in for the outside zone's lower density and shorter overlap. |
  | outside outflow Mach | 1.23 | Supersonic outside the LEH (M13: near-sonic at the throat, supersonic outside). |
  | inside multiplier | 1.23 | Order one: the linear kinetic gain needs little adjustment. |
  | Te(t), both zones (keV) | 2.79 (1 ns), 3.55 (3), 3.72 (4), 5.19 (4.75), 5.50 (5.5), 6.03 (6.5), 7.44 (8) | M13: 2.8 keV at the start of the main pulse (high fill). Michel 2009: 5–6 keV LEH in the reduced model. B16: 4.2 ± 0.2 keV at peak (0.8 scale, 1.6 mg/cc He equivalent); Te "follows the drive, peaks late". The fitted rise with the drive matches B16. The 7.4 keV late in the peak is above every published LEH value (L: low fill, no measurement exists). |

  - Fit quality: cone-fraction rms 0.016. The largest residuals are +0.041 at 8.4 ns
    (input cone fraction rising steeply), −0.025 at 2.0–2.6 ns (foot) and +0.035 at
    3.6–3.8 ns (before the clamp switch).
  - Outer-cone ratios against Fig. 9a (not fitted), at 5.6 / 7 / 8 ns:

    | Cone | Model | Fig. 9a |
    |---|---|---|
    | 44° | 0.41 / 0.54 / 0.65 | 0.43 / 0.54 / 0.64 |
    | 50° | 0.67 / 0.79 / 0.86 | 0.64 / 0.80 / 0.84 |

    The 44° cone loses more than the 50° cone, as in HYDRA.
  - Inner-cone ratios at 5.6 ns (not fitted): 23° 1.62 and 30° 2.20, against Fig. 9a's
    1.92 and 1.88. The model favours 30°, as I23 saw in x-ray spots, but HYDRA splits
    the transfer evenly.
  - Ion-wave amplitude inside at peak: 4×10⁻⁴ median, 1.2×10⁻³ max. The 0.0075 peak
    clamp therefore does not bind, as the lit review expected. The 0.002 foot clamp does
    bind, and the clamp switch produces Fig. 9b's step at 3.9 ns.
  - `n221204_cbet()` holds the fit rounded (rms 0.016, unchanged).
- **Δλ response.** Post-transfer cone fraction with Fig. 9's powers, Δλ at 1ω:

  | Δλ (Å) | 0 | 0.5 | 1.0 | 1.4 | 1.8 | 2.25 | 2.75 |
  |---|---|---|---|---|---|---|---|
  | 5.6 ns | 0.090 | 0.188 | 0.298 | 0.381 | 0.461 | 0.548 | 0.641 |
  | 7.0 ns | 0.121 | 0.207 | 0.298 | 0.364 | 0.426 | 0.491 | 0.561 |
  | 8.0 ns | 0.155 | 0.228 | 0.302 | 0.356 | 0.404 | 0.454 | 0.507 |
  | 2.0 ns | 0.128 | 0.142 | 0.153 | 0.160 | 0.169 | 0.179 | 0.194 |

  - dCF/dΔλ over 1–1.8 Å is 0.20 per Å at 5.6 ns, 0.16 at 7 ns and 0.13 at 8 ns (peak
    mean 0.16). The ROADMAP target is 0.15–0.18 per Å, derived from 40–50 µm/Å at
    2.7 µm/%CF.
  - The response is linear within 10 % over 0.5–2.75 Å at 5.6 ns.
  - At Δλ = 0 the net transfer is inner → outer (cone fraction 0.336 → 0.090 at peak).
    The null is at ≈1.2 Å, inside NC24's 1–1.55 Å for HYBRID-E.
  - Caveat: HYDRA's peak post-CBET cone fraction differs by only ≈0.07 between N210808
    (1.8 Å) and N221204 (2.75 Å): 0.55–0.6 against 0.65 in Fig. 1e. The calibrated model
    is therefore ~2.5× more sensitive to Δλ than HYDRA's cross-shot difference
    (different pulses, L). Beyond the clamp, the model has no saturation of transfer
    with Δλ.
- **Benchmark rerun.** `examples/n221204.jl` with transfer `calibrated`, now the
  example's default. Standard mesh (h 0.28 / capsule_h 0.16), 20 min on 8 threads.
  - The calibrated model acts on the example's delivered powers (DN24 Fig. 1, 442 TW),
    ~19 % above Fig. 9's.
  - Resulting post-transfer cone fraction: 0.21 (2 ns), 0.25 (4), 0.29 (4.5), 0.56 (5),
    0.69 (5.5), 0.66 (6), 0.63 (6.5), 0.59 (7), 0.52 (8).
  - At 4–4.5 ns this is below Fig. 9b's 0.35–0.40, because DN24's input cone fraction
    there (0.13–0.14) is below Fig. 9b's (0.17–0.18). At 5.5 ns it is 0.04 above,
    from the higher intensity.
  - Also changed in the example: the bubble's `capsule_scale` is now 1.135/0.91 = 1.25
    (*derived*). That is HYBRID-E's capsule against the ≈0.91 mm outer radius of the HDC
    basecamp SymCap (the "910 µm" HDC capsule class of the hydro-coupling review, M). It
    moves t_c to 5.41 ns from 5.9.

  | t (ns) | static, no CBET | static | +wall | +lips | +bubble | (b) Fig. 9b, +bubble, scale 1.25 | Fig. 10 |
  |---|---|---|---|---|---|---|---|
  | 2.0 | −0.017 | −0.057 | −0.032 | −0.032 | −0.032 | −0.043 | −0.002 |
  | 4.0 | −0.010 | −0.072 | −0.029 | −0.029 | −0.029 | −0.082 | −0.013 |
  | 4.5 | −0.020 | −0.088 | −0.048 | −0.048 | −0.043 | −0.120 | −0.091 |
  | 5.5 | −0.107 | −0.258 | −0.245 | −0.207 | −0.086 | −0.086 | −0.168 |
  | 6.5 | −0.105 | −0.227 | −0.187 | −0.151 | +0.031 | +0.030 | −0.057 |
  | 7.0 | −0.105 | −0.207 | −0.150 | −0.116 | +0.093 | +0.092 | +0.039 |
  | 8.0 | −0.104 | −0.177 | −0.083 | −0.052 | +0.168 | +0.165 | +0.090 |
  | 8.5 | −0.166 | −0.220 | −0.109 | −0.076 | +0.176 | +0.182 | +0.117 |

  Capsule Tr peak: 294.7 eV static, 303.7 eV +wall and 312.1 eV +lips (all at 8 ns),
  305.0 eV +bubble (at 5.5 ns).
- **Reading.**
  1. The reduced model now carries HYDRA's transfer.
     - Static P2/P0 at 5.5 ns is −0.26, against +0.01 with the single-zone deck
       settings and −0.23 with Fig. 9b's prescribed transfer.
     - With all dynamics, the calibrated and prescribed transfers agree within 0.007
       from 5.5 ns on.
     - The foot is closer to Fig. 10: −0.03, against −0.04 to −0.08 with the
       prescription (Fig. 10 ≈ 0). DN24's lower input cone fraction at 4–4.5 ns gives
       less transfer there than Fig. 9b's own.
  2. `capsule_scale` 1.25 moves the ridge stop 0.5 ns earlier (g = 0.83 already at
     5.5 ns).
     - The 5.5 ns minimum rises from −0.138 (scale 1) to −0.086, against Fig. 10's
       −0.168. The values at 6.5–8.5 ns are unchanged, since the ridge was already fully
       on there.
     - So the swing turns up ≈0.5 ns too early and ends 0.06–0.08 high.
     - The open items are the same as before: the bubble's inner absorption (48 % at
       6.5 ns and 85 % at 8 ns, above Ralph's 51–62 % at peak) and the wall motion
       (constant 0.10 mm/ns from the turn-on).
     - New: the collision time. Izumi's 0.4 ns per 20 % of capsule radius is
       extrapolated here to 25 % and to a 3.2 mm hohlraum. Whether r_cap or r_cap/R_hohl
       should set it is open; r_cap/R_hohl would give a scale of 1.12 and t_c ≈ 5.66 ns.
- **Open.**
  - This is one shot and one HYDRA run. Only the time history and the Δλ = 0 null
    constrain the split between the zones.
  - Te(t) absorbs everything time-dependent: ne, Ti, flow and the HYDRA power
    multipliers.
  - The 23°/30° split disagrees with Fig. 9a.
  - The model's Δλ sensitivity across shots is above HYDRA's (N210808 against N221204).
  - Next test of the Δλ response: a joint fit to N210808 (Fig. 1e) and N220919 (Fig. 9,
    black curves). Both are in the same figures.

## LEH shield (sombrero) (2026-10-24)

`Sombrero()` and `hohlraum(; sombrero = true)`. Dimensions measured from the CAD of a
THOR-class NIF symcap target (a non-public VisRad workspace; only the shield's own
dimensions are recorded here). The CAD is design intent, not metrology. **M**
- **Identification.** Each shield is a thin sheet around one end of the hohlraum, 18.4 mm
  across. The two shields are mirror images to the tessellation's precision.
- **Measurement.** The sheet's facets were cut by meridional half-planes every 5°, and
  height maps of both faces were ray cast on a 0.02 mm grid. The r–z sections and maps
  were checked by eye.
- **Bore and rim.** Inner radius 4.000 mm, outer 9.200 mm (vertex radii 3.9997–4.0002
  and 9.2000–9.2003 over all azimuths). Over one azimuthal sector, where the shield is
  mounted, the sheet runs on past 9.2 mm; that is not modeled.
- **Thickness.** 38.7 µm median (1–99 %: 37.7–40.7 µm).
- **Axial position.** The area-mean mid-plane lies 0.623 mm (±0.005 mm, grid and
  dimple sampling) inside the entrance-hole washer's inner face, i.e. toward the waist,
  beside the hohlraum wall rather than beyond the hole. That face is the LEH plane
  (`length/2`) of the matching `hohlraum` model.
- **Profile.** Flat on average: no cone or axisymmetric dimple. The sheet carries a
  hexagonal array of shallow dimples, 2.48 mm pitch, their centers sunk toward the
  waist. Dimple floor to the ridges between dimples is 0.21 mm (top face: floor −0.094,
  ridge points +0.119 mm from its own mean). The faces stay
  within ±0.12 mm of the mid-plane (0.5–99.5 %). Two mounting recesses go ~0.3 mm deeper.
  The rise from a dimple's floor is close to linear in the distance from its center:
  ~4° near the floor, ~12° toward the ridges.
- **Model.** `Sombrero` keeps the mean plane: a flat annulus (`annulus_mesh`, exact
  `DiskSurface`). The dimples are not axisymmetric; they move the true faces by up to
  ±0.12 mm, about the size of a hardware tolerance in `clearance`.
- **Material.** The CAD carries none. The default is `BLACK` (engineering default).
- **Radiation.** It plays no part. The shield sits beside the wall, 0.623 mm inside
  the LEH plane, so the hohlraum wall hides the hole from it. Its
  view-factor rows and columns with the hohlraum are exactly zero (`test/sombrero.jl`),
  `closed = true` balancing leaves those rows alone, and the hohlraum's solution is
  unchanged to 1e-12. One-sided, facing away from the waist: light from the port side
  is booked on its elements, light from the waist side (e.g. backlighter beams) as
  `sombrero_*_outside`.
- **Drive-beam clearance** *derived*. The test's cone-averaged NIF pointing (Build-A-
  Hohlraum example, 10.1 mm hohlraum, published cone aims) was run as 192 beams. The
  core 1e-4 envelope clears both shields by ≥ 2.09 mm in the worst case (outer 50°;
  44° 2.50, 30° 2.41, 23° 3.25), and the ray trace puts no core power on them. With
  the engineering `nif_wings()` the envelope clearance is ≥ 1.24 mm. The wings'
  10⁻³ scatter floor puts ~1e-4 of each beam on a shield, so every beam is `:review`
  (above 1e-5) at 1.2e8 W/cm². Drive beams
  reach the shield only past the hohlraum's outside, so it matters for beams aimed
  outside the hohlraum (backlighters), not for the drive.

## Hohlraum dynamics calibrated to published data (2026-10-24) *derived*

ROADMAP Next 9 and 10. Each piece of the hohlraum dynamics (wall, LEH lips, bubble, ridge
stop, two-zone CBET) is refitted to its own published anchor. PRE 109 Fig. 10 is not used
in any fit; it is the blind check at the end. Full texts are the accepted manuscripts in
the user's paper folder (not in the repository).

### Anchors extracted

| Anchor | Source | Form | Value | Stored |
|---|---|---|---|---|
| Dynamic VF hohlraum radius vs t (N170827) | Young et al., PoP 27, 082702 (2020) Fig. 3a | raster, axes from frame and gridlines, overlay checked | flat to 0.9 ns, in at 0.185 mm/ns, plateau 0.43–0.48 mm in through the peak, 0.55 at 7.5 ns | `young2020_radii()` |
| LEH radius vs t (same model) | Young Fig. 3a | same | 1.83 → 1.57 mm (0.857 r₀) linearly by ~6 ns; only the end value is SXI data | same |
| HYBRID-E wall edge at 6 ns | PRE 109 Fig. 1a (HYDRA) | raster, axis/waist/outline calibration, 1 px = 0.018 mm | waist 0.38 ± 0.04 mm in (laser-absorption edge); bubble tip ≈1.85 ± 0.1 mm (N210808), ≥ 1.55 (N221204, under its rays) at z ≈ 2.9–3.5 mm | `pre109_fig1a_wall()` |
| Bubble tip x_bub(t), 7 shots | Izumi et al., PoP 28, 022706 (2021) Fig. 8a/b | vector graphics (exact) | 1.38–2.28 mm over 4.65–6.19 ns; ±40 µm | `izumi2021_bubble()` |
| Inner-power fraction absorbed in the bubble | Izumi 2021 Fig. 10a (HYDRA circles, Eq. 6 line) | vector | 0, 0.020, 0.230, 0.550 at 3, 4, 5, 6 ns (basecamp, 0.3 mg/cc) | `izumi2021_absorption()` |
| Spot separation (ablator–wall collision) | Izumi et al., RSI 89, 10K111 (2018) | text | 5.8 ns basecamp (turn-on 0.3), 5.45 ns at 1.5× picket, ~0.4 ns earlier with the 1.18× capsule | docstring, tests |
| Peak-averaged inner absorption | Ralph et al., PoP 25, 082701 (2018) | text | 51 ± 2 % to 62 ± 2 % (CH high foot, 6.72 mm, 0.6 mg/cc, over the peak FWHM) | notes |
| Cone fraction for zero flux P2 | Ralph 2018 (VisRad at the measured peak wall position) | text | 18.5 %; peak F_in 0.34–0.38 → F' 0.14–0.18 (Fig. 5) | notes |
| Bubble density | Izumi 2021 Eq. 7 | text | ρ_bub = S ρ₀ R/(R − x), S ≈ 3 ∝ ρ₀^-0.1 (k = 1, settling the "illegible" exponent) | `density_exponent` 1.8 |

Digitization: Izumi 2021's figures are vector graphics in the accepted manuscript, so
markers and lines are read from the PDF paths (axes from the tick marks; exact to the
drawing). Young Fig. 3a and PRE 109 Fig. 1a are rasters: axes from frames, gridlines and
the drawn outline, stroke centres per column, overlays checked by eye.

### Fits (each against its own anchor)

| Parameter | Before | Now | Uncertainty | Anchor |
|---|---|---|---|---|
| Wall speed | 0.10 mm/ns | 0.185 mm/ns | 0.175–0.195 | Young Fig. 3a |
| Wall start | outer turn-on | turn-on + 0.8 ns (fit 0.78) | 0.68–0.85 ns | Young Fig. 3a |
| Wall saturation | none | 0.45 mm (default); 0.38 for HYBRID-E | 0.445–0.455 (±0.03 reading); HYBRID-E ±0.04, upper bound | Young; PRE 109 Fig. 1a waist |
| Bubble speed | 0.35 mm/ns | 0.36 mm/ns | 0.347–0.372 | Izumi Fig. 8a, 28 points, χ²/dof 0.89 |
| Bubble launch delay | 0.2 ns | 0.25 ns | 0.07–0.42 (correlated with speed) | Izumi Fig. 8a |
| Density exponent q | 1.862 (two-anchor fit) | 1.8 (physics) | — | Izumi Eq. 7: ρ_bub² ∝ (Sρ₀)², S ∝ ρ₀^-0.1 |
| Opacity τ_b | 0.144 (fit to T_b 29 %/54 %) | 0.111 | 0.100–0.124 stat.; 0.085–0.158 with plane ±0.2 mm and launch | Izumi Fig. 10a HYDRA absorbed fraction |
| Collision t_c | turn_on + 5.5 − 2.0(s − 1) − 0.7(E − 1) | turn_on + 5.5 (R − r_cap)/1.967 mm · E^-0.16 | ±0.15 ns (frame spacing) | Izumi 2018: 5.8, 5.45 (1.5× picket), −0.4 ns (1.18× capsule; the gap gives −0.44) |
| HYBRID-E picket fluence `energy` | 1 | 1.1 (example) | 0.8–1.5 | 1.5× basecamp energy (DN24 vs Izumi Fig. 3a) × 0.73 for split outer quads; Fig. 1a N210808 tip 1.85 ± 0.1 mm at 6 ns |
| LEH lips (example) | 5.0 ns step over 0.5 ns to 0.84 r₀ | smoothstep from −0.1 ns over 7.85 ns to 0.857 r₀ | half-closed 3.0–3.8 ns (timing mapping) | Young Fig. 3a LEH radius (end value SXI; ramp assumed) shifted 0.8 ns |

Checks of the fits:
- *Wall.* Young's radius is flat to 0.9 ns, falls at 0.185 mm/ns to 0.45 mm in by 3.3 ns
  and then holds (0.43–0.48 mm) through the peak: the motion happens in the foot, not at
  a constant speed through the pulse. The fit's residuals are < 0.03 mm over 0.5–6.75 ns
  (`test/dynamics.jl`). On N221204 the old wall (0.10 mm/ns from the turn-on) was
  0.56 mm in at 6 ns and 0.86 at 9 ns; the new one is 0.38 from 3.3 ns on.
- *Bubble tip.* The seven shots' residuals to the global fit (E^0.2 ρ₀^-0.2 scalings,
  one speed, one launch) are −0.09 to +0.07 mm, rms 0.036 mm, inside the ±0.04 error bars.
  Per shot, the slopes are 0.33–0.43 mm/ns and the intercepts 0.14–1.02 ns.
- *Absorbed fraction.* The opacity is degenerate with the bubble's plane and launch: with
  the plane where the proxy pointing's 50° axes reach the wall (z = 2.59 mm instead of
  Izumi's 3 mm) τ_b = 0.086. The absorbed fraction against time is the same in every
  variant (0.03 / 0.20 / 0.36 / 0.56 at 4 / 5 / 5.5 / 6 ns, against HYDRA's 0.02 / 0.23 /
  – / 0.55 and Izumi's Eq. 6 0.04 / 0.20 / 0.32 / 0.56). The old defaults gave 0.09 / 0.32
  / 0.50 / 0.69 on the proxy geometry (plane 2.59 mm): too much, too early.
- *Collision.* Izumi 2018's 1.18× capsule (outer radius 1.064 mm against 0.908) narrows
  the capsule–wall gap by 8 %, which at constant closing speed (0.36 mm/ns) brings the
  collision 0.44 ns earlier: the measured "roughly 400 ps". The 1.5× picket's 5.45 ns
  gives the energy exponent −0.16. HYBRID-E's gap (3.2 − 1.135 = 2.065 mm) is 5 %
  *wider* than the basecamp's, so t_c = 0.4 + 5.77 · 1.1^-0.16 = 6.09 ns. The old
  extrapolation in capsule radius alone gave 5.41 ns.
- *Ralph 2018.* Energy-weighted over N221204's peak FWHM (4.62–8.76 ns), after CBET, the
  calibrated bubble absorbs 0.575 of the inner power (old defaults 0.50, E = 1.5 0.65),
  inside Ralph's 51–62 %. The platforms differ (CH high foot, 0.6 mg/cc, 14 ns pulse), so
  this is a plausibility check, not a fit. The earlier "85 % against 51–62 %" compared an
  instantaneous value at 8 ns with a peak average. Ralph's 18.5 % cone fraction for zero
  flux P2 at the measured peak wall position is not reproduced: it needs the high-foot
  pointing and the measured wall radius, neither of which is published.
- *HYBRID-E bubble (PRE 109 Fig. 1a, 6 ns).* N210808's tip is 1.85 ± 0.1 mm from the
  wall at z ≈ 2.9–3.5 mm; N221204's is ≥ 1.55 mm (its denser inner-beam rays are drawn
  over the plasma). The calibrated model with `energy` 1.1 puts it at 1.96 mm, in the
  plane z_b = 3.49 mm.

### Two-zone CBET: the 7.4 keV LEH temperature

`n221204_cbet()`'s Te(t) reaches 7.44 keV at 8 ns, above every published LEH value. Refit
to the same Fig. 9b data with Te ≤ 5.5 keV (`calibrate_cbet` with an upper bound,
16 min on 4 threads: `examples/n221204_cbet_calibration.jl` with `bounds[:Te] = (300, 5500)`):

| Fit | Fig. 9b rms | Te (keV, 1 → 8 ns) | multipliers out / in, Mach | 44° / 50° at 5.6 ns (Fig. 9a 0.43 / 0.64) | CF at Δλ 0 | Δλ null | CF at 1.8 Å, 5.6 ns |
|---|---|---|---|---|---|---|---|
| `n221204_cbet()` | 0.0163 | 2.79 → 7.44 | 0.30 / 1.23, 1.23 | 0.41 / 0.68 | 0.090 | 1.20 Å | 0.461 |
| A: Te ≤ 5.5 keV | 0.0160 | 1.75 → 5.01 | 0.07 / 0.41, 0.91 | 0.47 / 0.61 | 0.281 | 0.50 Å | 0.499 |
| B: A + ne/nc(t) | 0.0115 | 4.87, 1.87 … 5.50, 4.88 | 0.04 / 1.09, 0.33 | 0.48 / 0.58 | 0.333 | 0.05 Å | 0.513 |

- Fig. 9b does not set the absolute Te: A fits it as well with every Te ~0.67× lower,
  the outside multiplier 4× smaller and the inside one 3× smaller. The 7.4 keV is the
  overall gain scale and the balance between the zones, not a measured temperature.
  The rise of Te with time (×2.7–2.9 from 1 to 8 ns) is what the data constrain.
- What separates the fits is the Δλ response, which is not fitted. The reference puts
  the null at 1.2 Å, inside NC24's 1–1.55 Å for HYBRID-E. A moves it to 0.5 Å. A
  matches HYDRA's N210808 (1.8 Å) post-CBET cone fraction better (0.50 against ~0.55–0.6
  in Fig. 1e; the reference 0.46), i.e. it halves the excess Δλ sensitivity.
- B fits Fig. 9b better only with an unphysical density history (ne/nc 0.15 at 1 ns,
  0.01 at 4.75–5.5 ns) and moves the null to 0.05 Å: rejected.
- Not adopted: `n221204_cbet()` is unchanged. Constraining Te trades one unfitted check
  (the null) for another (N210808). A joint fit to N210808 / N220919 with the null as a
  constraint (ROADMAP 10) is what can fix the temperature scale. Done: "Joint multi-shot
  CBET calibration (2026-10-24)" below. Te ≤ 5.5 keV costs the three-shot fit little; the
  null in 1–1.55 Å is what the cross-shot data resist.

### N221204 benchmark rerun: the blind check against Fig. 10

`examples/n221204.jl` (default `dynamics = calibrated`; `old` reruns the 2026-10-23
settings), calibrated two-zone CBET, mock HYBRID-E pointing, standard mesh (h 0.28 /
capsule_h 0.16, 31 min for the five runs on 4 threads). Settings: wall from 1.2 ns at
0.185 mm/ns stopping 0.38 mm in; lips to 0.857 r₀ along Young's ramp shifted 0.8 ns
(closure −0.1 ns, width 7.85 ns); bubble defaults with `energy` 1.1 and
`capsule_radius` 1.135 (t_c 6.09 ns). `energy` was set from the picket energy and quad
splitting before the comparison was run.

| t (ns) | static | +wall | +lips | +bubble (now) | +bubble (2026-10-23) | Fig. 10 |
|---|---|---|---|---|---|---|
| 2.0 | −0.057 | −0.034 | −0.031 | −0.031 | −0.032 | −0.001 |
| 4.0 | −0.072 | −0.026 | −0.017 | −0.016 | −0.029 | −0.018 |
| 4.5 | −0.088 | −0.053 | −0.040 | −0.034 | −0.043 | −0.092 |
| 5.5 | −0.258 | −0.266 | −0.239 | −0.166 | −0.086 | −0.169 |
| 6.5 | −0.227 | −0.233 | −0.202 | +0.029 | +0.031 | −0.045 |
| 7.0 | −0.207 | −0.210 | −0.179 | +0.085 | +0.093 | +0.046 |
| 8.0 | −0.177 | −0.175 | −0.145 | +0.123 | +0.168 | +0.095 |
| 8.5 | −0.220 | −0.229 | −0.196 | +0.126 | +0.176 | +0.109 |

| Swing metric (full model) | 2026-10-23 (coarse mesh) | now (standard) | Fig. 10 |
|---|---|---|---|
| minimum | −0.164 at 5.00 ns | −0.168 at 5.25 ns | −0.173 at 5.58 ns |
| upward zero crossing | 6.25 ns | 6.26 ns | 6.67 ns |
| at 8 ns | +0.162 | +0.123 | +0.095 |
| mean over 5.25–8.5 ns | +0.054 | +0.024 | −0.001 |
| rms against Fig. 10, 2–8.5 ns | 0.051 | 0.035 | — |

The coarse mesh (0.8 / 0.4) reproduces these within 0.006. Capsule Tr peaks at 301.8 eV
(5.5 ns) with everything on, 307.2 with wall and lips (drive Trad 313 eV).

**Why the swing changed.**
1. *Wall.* Saturating at 0.38 mm takes +0.09 off the late P2 (the +wall column now
   stays at the static value after 5.5 ns; before, it rose to −0.083 at 8 ns). This is
   most of the late improvement (+0.168 → +0.123 at 8 ns).
2. *Collision.* The ridge stop now starts at 6.09 ns instead of 5.41, so at 5.5 ns the
   inner light still reaches the wall where it is not absorbed: the 5.5 ns value drops
   from −0.086 to −0.166 and the minimum deepens to Fig. 10's level.
3. *Bubble.* On HYBRID-E's geometry the refitted bubble absorbs about as much as before
   by 5.5 ns (0.22 against 0.17) but more later (0.57 at 6.5 ns against 0.47, 0.92 at
   8 ns against 0.84). With the ridge on from ~6.2 ns this sets the rise from −0.17 to
   +0.03 between 5.5 and 6.5 ns, 0.4 ns ahead of HYDRA's.
4. *Lips.* Closing along Young's ramp instead of a step at 5.0 ns moves P2 by < 0.01
   from 4 ns on (sensitivity run `late_lips`: minimum −0.169 at 5.0 ns, crossing 6.25,
   +0.118 at 8 ns, coarse mesh).

**Sensitivity to the least constrained input** (coarse mesh, full model). The HYBRID-E
picket fluence `energy` is 1.1 from the picket energy (1.5×, by eye) and split outer
quads (−27 %); PRE 109 Fig. 1a's N210808 bubble tip (1.85 ± 0.1 mm at 6 ns) gives 0.8
(0.6–1.05) on its own. E = 1.5: minimum −0.138 at 5.25 ns, crossing 5.99 ns, +0.121 at
8 ns. E = 0.8: −0.184 at 5.5 ns, 6.62 ns, +0.109 (rms 0.024). The crossing time follows
the bubble speed: 0.6 ns per factor 1.9 in E, i.e. per ±7 % in speed. The default is
kept at 1.1, the value set before the comparison; 0.8 matches Fig. 10 better and is
consistent with Fig. 1a, but choosing it now would be tuning to the blind check.

**N220919 and N210808 proxies** (coarse, `shot` argument; same mock pointing, CBET
calibration and hohlraum; N220919 on N221204's pulse with Δλ 2.5 Å; N210808 with the rise
0.2 ns and the peak end 0.45 ns earlier and Δλ 1.8 Å):

| | minimum | crossing | at 8 ns | rms 2–8.5 ns |
|---|---|---|---|---|
| N220919 now (old) | −0.140 at 5.25 (−0.139 at 5.0) | 6.23 (6.14) | +0.118 (+0.165) | 0.033 (0.053) |
| Fig. 10 N220919 | −0.145 at 5.6 | 6.66 | +0.094 | |
| N210808 now (old) | −0.092 at 5.25 (−0.095 at 5.0) | 6.09 (5.67) | +0.119 (+0.166) | 0.061 (0.077) |
| Fig. 10 N210808 | −0.141 at 5.5 | 6.68 | +0.084 | |

N220919 behaves like N221204. N210808's proxy has a positive foot (+0.03 at 4 ns against
−0.05) and a shallow minimum: its pulse and foot cone fraction are only approximated
(N221204's DN24 powers, shifted), and the CBET calibration is N221204's (its Δλ
sensitivity is ~2.5× HYDRA's across shots, ROADMAP 10).

**Remaining discrepancies and likely causes.**
- Crossing ~0.4 ns early: the HYBRID-E bubble speed (above), and the bubble plane:
  z_b follows the mock pointing's 50° axes (3.49 mm), which sets how much of the inner
  beams the bubble clips; the absorbed-fraction calibration was done on a 5.75 mm
  hohlraum with Izumi's plane.
- 4.5 ns 0.06 above Fig. 10: DN24's delivered input cone fraction (0.13–0.14) is below
  Fig. 9b's (0.17–0.18) in the rise, so the CBET gives less transfer there.
- Foot ~0.03 below: mock pointing and the foot cone fraction (read by eye).
- Late +0.02–0.03 high: the capsule does not shrink (Young's ablation front is at 0.83 mm
  at 5.5 ns, falling fast after 6 ns), no laser degradation or backscatter.
- Every dynamic piece is still one effective parameter set from other platforms (HDC
  SymCap in 5.75 mm Au, N170827); only the wall saturation and the bubble-tip check come
  from HYBRID-E itself, and both from the HYDRA picture of Fig. 1a.

## Young 2020 validation (2026-10-24) *derived*

ROADMAP Next 9 ("validation still to run") and 11 (open items). Validation, not tuning: no
calibrated default changed. Full texts are in the user's paper folder (not in the
repository); data digitized from them are in `data/` with provenance headers.

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
hohlraum, 1.82 mm LEH, 1.05 mm corners and 0.98 mm capsule with "hohlraum, LEH, and
capsule radii all reduced to 85 %", albedo 0.8 (Au) / 0.3 (ablator), conversion 0.9, the
*requested* (azimuthally symmetric) pulses, no windows, and the velocity from the N170827
Green's function integrated over the whole pulse (negative in the foot, peaked at the
rise and at 5.8 ns). Mode 1 is "the vector norm of the three P1m components",
c_lm = 100 a_lm/a00 (Briareus `mode1_vector`; Young's x/y normalization of the associated
Legendre functions was not checked, the z one is the same). Offsets: the whole target
moved against fixed beams, or the capsule moved inside a fixed hohlraum. N170827's
pointing is not published; Briareus uses the tutorials' pointing with z scaled by
11.3/10.1, and unscaled as a check.

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
- *Setup difference:* the ROADMAP 11 mismatches (inner/outer 1.3–1.6, horizontal offset
  0.60 km/s/µm) came from the implosion-metrics tutorial's hohlraum (10.1 × 5.75 mm,
  0.9 mm capsule, albedo 0.85, a static velocity at 100 km/s/%). A horizontal
  whole-target offset moves the inner and outer footprints so that their mode 1
  opposes; in Young's geometry they nearly cancel at peak (inner beams alone
  0.0056 %/µm, outer alone 0.0067 %/µm, together 0.0012), on the tutorial's they do not
  (0.030, 0.009, together 0.0059). The foot, with another cone balance, does not cancel
  (Young's Fig. 14a: ~0.5 % per 100 µm in the foot, ~0.1 % at peak).
- *Model error, fixed:* the mode analysis (`legendre_weights`, `ylm_weights` and the
  `Scenario` caches, hence `drive_history`) measured directions from TCC even when the
  capsule was offset. That mixes the capsule's P2 (−0.1 here) into a spurious mode 1 at
  first order in the offset: 0.0074 %/µm of fake horizontal mode 1 (0.74 km/s/µm)
  against 0.00005 %/µm about the capsule's centre, and the vertical one had the wrong
  sign. They now default to `mode_center` (the centre of the observed sphere; TCC
  otherwise), so centred capsules are unchanged bit for bit.
- The 100 km/s/% default against Young's kernel: a flat 100 km/s/% over 4–7.75 ns gives
  11–17 % less than G on these traces (∫G ≈ 117 km/s/%). The default is unchanged.

**Static vs dynamic (Young Sec. II A, Fig. 3).** Young's dynamic model (radii and albedos
of Fig. 3a against time) is built with `Keyframes` at every 0.25 ns step: wall 3.1 mm minus
Fig. 3a's displacement, LEH and capsule scaled to Table II, corners scaled with the wall,
`TabulatedAlbedo` for wall and capsule. N170827's delivered per-beam powers are not
published, so the mode-1 source is a proxy: a fixed 5 % rms per-quad imbalance on the
requested pulses (and, separately, a dropped 23.5° quad). Results (h 0.5):

| | imbalance proxy | dropped 23.5° quad | Young (delivered N170827) |
|---|---|---|---|
| mean abs. difference of the mode-1 norm, all times | 0.18 % (mean norm 0.48) | 0.73 % (2.29) | 0.70 % |
| same, peak 4.5–7 ns | 0.10 % (0.34) | 0.72 % (2.15) | 0.14 % |
| velocity static / dynamic (km/s) | 52.2 / 40.7 (+28 %) | 336 / 260 (+30 %) | 40.5 / 38.7 (+4.6 %) |
| direction difference | 2.5° | 0.8° | (75°, 265°) vs (67°, 270°), ~8° |

Decomposition (imbalance proxy, static 52.2 km/s): dynamic radii with static albedos
51.1; the capsule shrinking alone 50.4; static radii with Young's dynamic albedos 41.7;
static radii with a wall albedo of 0.87, 41.2. So **the 85 % radii stand in for the
moving wall, closing LEH and ablating capsule to 2 % in velocity**, as Young found; **the
albedo does not**: his dynamic wall albedo is 0.85–0.88 through the peak against the
static 0.8, and the extra re-emission smooths the mode 1 by ~20 % (the same physics as
in his foot, where the lower albedo raises it). Young reports the static velocity only
4.6 % higher. Not resolved: his delivered imbalance has its own time structure (larger
in the foot and rise, where the albedos are closer), which the proxy lacks; the
delivered per-beam powers would settle it. The mean-difference numbers are on Young's
scale (0.10 against 0.14 % at peak for a similar peak mode 1) but depend on the source's
amplitude.

**Nora 2023 (ROADMAP 11): the transcription, checked against the paper.**
- Eq. 4's high branch is printed with rounded coefficients (13 ln(Y + 1585) − 95). The
  fit drawn in Fig. 6 has a = 13.105, b = 1585.225, c = −95.539; since ln(Y + b) ≈ 7.4,
  the rounding moves the result by −0.24 (1.43 against 1.67 at Y_amp = 80). The two
  branches still do not meet (the low one is 2.08 at 80), and the paper fits nothing
  between 40 and 80 (burn fractions ~6 %, where the burn moves into the dense fuel).
  Implemented with the Fig. 6 coefficients on 80–500; 40–80 is refused.
- Table IV lists the inverses (factors from the measured to the unamplified speed):
  1/Eq. 4 reproduces its burn column to 1e-3 (N130530 0.956, N210808 0.570, N220919
  0.580; test suite), and total = shape × burn (N210808: 0.963 × 0.570 = 0.549; the
  text's "184 %" is 1/0.549 = 1.82).
- The shape term: Eq. 2 is the 2D fit with P2 as a fraction (1.31 P2 + 0.71 below 0.35,
  −0.65 P2 + 1.4 above, continuous to 4e-3; 0.71 for a round hot spot because the 2D
  points are shifted from the 3D ones, Fig. 5). The "80–120 % for −40 % < P2 < 40 %" of
  the conclusion describes Table I (3D), not Eq. 2; that is where the old transcription
  went wrong. Both are implemented (`shape_amplification(P2)`,
  `shape_amplification(P2, θ)`). Table IV's shape column follows Table I for equatorial
  flows (N130530: 0.965 exactly) but not for polar ones from the listed angles (N210808
  would need θ ≈ 77° where its flow is at 166°; N210531 matches no angle): the angle
  convention of Table IV is not recoverable from the paper.

**Callahan 2018 and Izumi 2018 (trends; the repository had no hot-spot transfer).**
Callahan's Fig. 6 line (digitized from its dotted fit, 29 pixels at the two ends):
P2_hs − 2.7(CF − 33) = 101.9 − 85.0 x µm, x = √(E_picket,outer/(A_outer ρ_fill)) τ/R ·
r_cap/R, zero at x = 1.199 (his Eq. 4: 1.2), σ 8 µm over 39 shots. The per-shot inputs are
not tabulated, so the curve cannot be recomputed shot by shot; Briareus compares
sensitivities. Transfer: Callahan's 2.7 µm per 1 % cone fraction divided by Briareus's
own d(flux P2/P0)/dCF (−0.293 % per %, peak 4.5–6.5 ns) gives −9.2 µm of hot-spot P2 per %
of flux P2 (pole-hot flux → oblate). Setup: 10.1 × 5.75 mm Au hohlraum (albedo 0.8),
0.908 mm capsule, the tutorials' pointing, N170827's pulses, `BubbleAbsorption` with its
calibrated defaults (turn-on 0.25 ns), h 0.45. Izumi 2018 (LLNL-PROC-751038)
single-parameter SymCaps, core P2 against the basecamp's +5.6 µm (the assignment of
Fig. 6b's points to shots follows the image order of Fig. 6a, **L**):

| change | Briareus flux P2 (pp) → hot spot (µm) | Callahan's line at x ≈ 1.2 | Izumi 2018 measured change |
|---|---|---|---|
| 1.5× picket energy | +3.25 → −30 | −23 | −21 to −24 |
| 1.18× capsule | +1.31 → −12 (view factors alone −0.59 → +5.5) | −18 | −21 to −24 |
| 0.15 mg/cc fill | −2.56 → +24 | −42 | +10 to +15 |
| quad splitting | not modeled (the outer fluence change is not given) | — | +10 to +16 |

All three Briareus signs match Izumi, magnitudes within a factor of 2. Callahan's fill
term gets the low-fill shot's sign wrong (Izumi's own point: the 0.15 mg/cc shot is an
outlier of the bubble-size picture). Not blind: the bubble opacity's density exponent was
fitted to Izumi 2021's 0.15/0.3 mg/cc anchors and the collision time to Izumi 2018's
separation times (basecamp, 1.5× picket, 1.18× capsule), so this checks the composition
(bubble, view factors and transfer), not the inputs. Ralph 2018's 18.5 % anchor stays
skipped (unpublished pointing).

## Wall material, conversion and capsule energetics (2026-10-24) *derived*

ROADMAP Next 16's model changes, each with its data check. Full texts read (accepted
manuscripts, in the user's paper folder, not in the repository): Döppner et al., PRL 115,
055001 (2015) (OSTI 1209256); Schein et al., PRL 98, 175003 (2007) (OSTI 908388); Rosen,
"Fundamentals of ICF hohlraums" (2005, OSTI 885380), §VI; Callahan et al., PoP 27,
072704 (2020); Zylstra et al., PRE 106, 025202 (2022). Confidence: **H** primary text,
**L** our arithmetic or construction, *calibrated* / *derived* as marked.

> **Superseded in part (2026-10-25).** The models and inputs stand; these values do not:
> - the DU factor f = 0.785 (§1) is **0.836**, recalibrated at the wall's effective age
>   ("Energetics defaults and the capsule share (2026-10-25)");
> - the "0.92 mm ablation front at 6 ns" (§3) is the ice–gas interface; the example's
>   capsule schedule is now N221204's own K24 radius history ("The P2 swing: capsule
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
| DU wall equivalence | "equivalent to an additional 6.5 % of laser power (≈25 TW) at the peak" | Döppner 2015, citing Kline 2013 and Callahan 2012 | H |
| Same, as a drive factor | f_DU = 1.065 in T^3.3 ∝ E f_cap f_LEH f_DU, "in both simulations and experiments" | Callahan 2020 Eq. 5 | H |
| DU hohlraum build | 7.0 µm DU lined with 0.7 µm Au, 22.3 µm Au backing; the liner keeps the picket as in Au | Döppner 2015 | H |
| Au/DU pair, high foot (T-1) | N140311 Au 391 TW, η_L 0.857, 317 ± 4 eV; N140520 DU 393 TW, 0.863, 322 ± 4 eV | Döppner 2015 Table I | H |
| Cocktail/Au loss ratio | E/A: Au (0.598, 0.398) T^3.3 t^0.6, U₀.₆Dy₀.₂Au₀.₂ (0.604, 0.407) T^3.1 t^0.57; ratio ∝ T^-0.22 t^-0.05; "about 0.83" for NIF | Rosen 2005 §VI | H |
| Same, LASNEX | 97.6 / 90.8 / 88.6 % at 185 / 275 / 305 eV (1 ns); 17.5 % less wall loss → 10 % less laser energy at 300 eV | Schein 2007 | H |
| Conversion in time | 0.7 t^0.12 makes Nova's 1 ns T ∝ t^0.15 (Lindl 2004 Fig. 4-4); 0.7 t^0.2 (Rosen); 0.85 t^0.2 (Schein) | energetics.md §2.2 | H |
| Conversion in intensity | 0.95 at 10¹⁴, 0.80 at 10¹⁵ W/cm² (spheres, backscatter corrected) | Dewald 2008 abstract | M |
| NIF hohlraum η | 85–90 % at 235–345 eV | Olson 2012 abstract | M |
| N221204 split | wall 54.1 % (Au lining 4.4), LEH 25.2 %, gas 4.5 %, capsule 16.2 %; laser degraded 9–13 % | PRE 109 App. 3 | H |
| N221204 backscatter | SBS 42 ± 11 kJ early in the peak, mostly inner cones; SRS 6 ± 1 kJ | AS24 | H |
| LEH pair | N191118 3.64 mm LEH, ≈465 TW; N210601 3.1 mm, ≈390 TW; Dante Trad "comparable until ~7.5 ns" (≈295 / ≈300 eV) | Zylstra 2022 Fig. 2 and text | H, H-img |

### 1. Au-lined DU wall (`du_wall_loss()`, `uranium()`)

`WallLossAlbedo` gains `factor` f, `factor_exponent` δ (at `factor_temperature` 300 eV),
applied to gold's clamped loss, and an optional gold `liner` (areal mass) during which
the loss stays gold's (an energy-content blend w = min(1, liner/m_F) with m_F the
constant-T `HeatFront` mass; *engineering*). Defaults unchanged (gold, bit for bit).

- δ = −0.2 (Rosen −0.22; Schein's three LASNEX points give −0.19).
- **f = 0.785, calibrated** (*superseded*: 0.836 at the wall's effective age, see the
  note above) to the 6.5 % equivalence: Döppner's high-foot hohlraum
  (5.75 × 9.43 mm, 3.1 mm LEHs, 1.125 mm capsule at albedo 0.3), 390 TW × 0.86 coupling,
  conversion 0.9, a wall source in two bands (30 % at |z| 0.5–1.7 mm, 70 % at 3.0–4.2),
  at 13 ns, h 0.4: the Au wall needs 6.5 % more power for the DU wall's capsule Tr.
  Range 0.778–0.789 over 12–14 ns and conversion 0.75–0.95 (*L*).
- **Check, Döppner's T-1 pair** (not fitted): at equal power the DU wall raises the
  capsule Tr by 4.9 eV (294.2 → 299.1); measured 317 → 322 ± 4 eV (+5 ± 6) at 391 / 393 TW.
- **Rosen's 0.83 / Schein's 0.825 give only 5.0–5.1 %** here (3.9 eV): in this model
  gold's Lindl loss leaves the wall ~30 % of the x-ray source at a NIF peak (α ≈ 0.94),
  far from Schein's wall-dominated accounting (17.5 % → 10 % implies ~57 %). Pure DU
  should be no better than the cocktail, so the calibrated f below Rosen's cocktail value
  says the Lindl gold loss is low at late times, not that DU beats the cocktail.
- On HYBRID-E (11.24 × 6.4 mm, same bands, 350 TW) the same f is worth 8.4–9.0 % of
  power (6.0–6.5 eV) at 6–8 ns: the larger case-to-LEH ratio gives the wall a larger
  share. In the N221204 example (full dynamics, nif conversion, ablator) DU raises the
  peak capsule Tr from 303.3 to 310.6 eV.
- Liner: 0.7 µm Au (Döppner; HYBRID-E's is unpublished) is crossed in ~0.4 ns at 150 eV;
  at 300 eV, 8 ns it is 6 % of the heated mass (PRE 109 puts 4.4 of 54.1 % of the wall
  heating, 8 %, in the lining; *L*, consistent in size).

### 2. Conversion as a model (`Conversion`, `nif_conversion()`)

η(t, I) = clamp(η₀(I) ((t − start)/1 ns)^s, range), η₀ constant or a table in log₁₀ I of
the element's absorbed flux. Materials with a number keep the scalar path bit for bit
(tested); intensity tables are generic in the number type (ForwardDiff through `solve`
tested); `flux_response(...; laser)` linearizes an intensity-dependent conversion (secant,
exact at `laser`); `Scenario`'s flux response evaluates a time-dependent one at its time.

- **Closed-form check**: with the wall loss dominating (pinhole LEHs, reflecting capsule,
  Nova scale-1 cylinder 1.6 × 2.55 mm, 20 TW), T ∝ t^((0.38 + s)/3.3): model slopes
  0.111 / 0.147 / 0.171 for s = 0 / 0.12 / 0.2 over 0.3–1 ns (closed form 0.115 / 0.152 /
  0.176).
- **Lindl's Nova statement** (0.7 t^0.12 → T ∝ t^0.15) is *not* reproduced with open
  LEHs: on the scale-1 hohlraum with 1.2 mm LEHs the slopes are 0.052 / 0.086 / 0.108 for
  s = 0 / 0.12 / 0.2. The LEHs take ~36 % of the source and flatten T(t); Lindl's
  analysis is wall-dominated. Matching T ∝ t^0.15 with open holes would need s ≈ 0.3
  (*L*), so s stays an input.
- **NIF preset** `nif_conversion()` = 0.7 t^0.12 ≤ 0.95 (*derived*): 0.80 / 0.85 / 0.90
  at 3 / 5 / 8 ns, Olson's 85–90 % at NIF peak; the cap is Dewald's 10¹⁴ plateau and
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
- **Calibration** (*superseded*, see the note above): α_c = 0.22 gives E_cap = 16.2 % of the coupled energy (285 kJ of
  1762) in the full energetics configuration below (α_c 0.20: 16.5 %). It is degenerate
  with the schedule and with the hohlraum temperature: in the configuration whose LEH
  share matches HYDRA's (v4 below, Tr 307.7 eV) α_c = 0.2 gives 12.8 %, and 16.2 % would
  need α_c ≈ 0 (or a larger effective radius than the ablation front). So "α_c ≈ 0.2"
  is consistent with Swadling only together with the too-hot drive of §4.
- Kr24b's 16.2 % is read as a fraction of the laser energy in the simulation (after the
  9–13 % degradation): the four fractions sum to 100 % with no laser loss. As a fraction
  of the delivered energy the same run gives 14.1 %.

### 4. N221204 benchmark: energetics and the Fig. 10 blind effect

`examples/n221204.jl ... full` with the new `key=value` options; full model (calibrated
CBET, wall, lips, bubble); laser 2.03 MJ by eye. Backscatter: SBS 16.7 % of the inner
power over 5.5–6.5 ns (ramps 5.0–5.5 and 6.5–7.5) and SRS 1.14 % over 5.5–9.0 ns, scaled to
42 and 6 kJ (shapes *L*); degradation 11 % on every beam (the middle of 9–13 %); both
before CBET (whose calibration used HYDRA's degraded powers). Coarse mesh (0.8 / 0.4)
unless "std" (standard 0.28 / 0.16, 10 min per run on 4 threads; coarse within 0.006 in
P2 and 0.1 eV).

Energy split, % of the coupled energy (after backscatter and degradation); wall = all
surfaces' heating (unconverted laser + absorbed x-rays; the bubble's light is conducted
to the wall, so HYDRA's wall + gas, 58.6 %, is the comparison):

| Run | wall | LEH | capsule | peak capsule Tr (eV) | P2/P0 min (ns) | zero | 8 ns | rms vs Fig. 10 |
|---|---|---|---|---|---|---|---|---|
| HYDRA (PRE 109 App. 3; Fig. 10) | 58.6 | 25.2 | 16.2 | ~300 internal, Dante 313 | −0.170 (5.6) | 6.67 | +0.095 | — |
| before, std (Au, η 0.75, ablator 0.3 at R₀, full laser) | 57.6 | 22.2 | 20.1 | 301.8 at 5.5 ns | −0.168 (5.25) | 6.26 | +0.123 | 0.036 |
| v1 + backscatter + 11 % degradation | 58.2 | 22.1 | 19.8 | 288.9 at 8.0 | −0.138 (5.25) | 6.19 | +0.119 | 0.042 |
| v2 = v1 + nif conversion | 50.8 | 25.9 | 23.3 | 303.3 | −0.136 | 6.19 | +0.117 | 0.042 |
| v3 = v2 + DU wall | 46.4 | 28.2 | 25.4 | 310.6 | −0.127 | 6.20 | +0.106 | 0.041 |
| v4 = v1 + capsule model (α_c 0.2) | 62.1 | 25.1 | 12.8 | 307.7 | −0.140 | 6.28 | +0.058 | 0.041 |
| v5 = v4 + nif conversion | 55.4 | 29.7 | 15.0 | 323.4 | −0.139 | 6.28 | +0.056 | 0.042 |
| after, std (`energetics=new`: v5 + DU, α_c 0.22) | 50.8 | 33.0 | 16.2 | 334.2 at 8.0 | −0.134 (5.25) | 6.31 | +0.046 | 0.044 |

Findings:
- **The split has two independent numbers** (wall = 1 − LEH − capsule with no other
  sink): the LEH share fixes the drive temperature, the capsule share then fixes
  α_c(R/R₀)². With everything published switched on (`after`), the LEH share is 33 %
  against 25 % and the capsule Tr peaks at 334 eV against HYDRA's ~300 (Dante 313):
  (25.2/33.0)^¼ × 334 = 312 eV, i.e. matching the LEH share and matching Dante are the
  same correction. The x-ray wall loss is 37.6 % of the coupled energy; HYDRA's split at
  this conversion needs ~45 %, a wall loss ~1.6× this model's at the lower temperature.
- **What overshoots**: Lindl's gold loss (α ≈ 0.94 at the HYBRID-E peak) times the
  calibrated DU factor, with η 0.87, and the shrinking capsule that removes the late sink
  (v1 → v4: +19 eV). Each piece is from data; together they are too reflective. The
  configuration that matches the LEH share (v4: Au loss, η 0.75) does so by the old
  constant conversion standing in for the missing wall loss. Missing sinks here: the He
  gas (4.5 %), the LEH windows, the M-band (Dewald 2020, not gray), and HYDRA's NLTE wall
  emission. Not fixed by tuning: no default was changed to close it.
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
temperature would need the LEHs to carry ~55 % of the 3.64 mm hohlraum's loss at 7 ns,
implausible against PRE 109's 25 % for the whole pulse; so the two shots' different
histories (N191118's pulse ends at ~7.5 ns, N210601's continues; different feet) and
Dante's LEH view (it sees the outer-beam bubble region, Ch24), not the lips, are the
likely difference. The pair checks the trend only (*L*).

## Spectral radiosity, M-band and a synthetic Dante (2026-10-24) *derived*

ROADMAP Next 13 and item 16's first absolute check. Read (accepted manuscripts in the
user's paper folder, not in the repository): Dewald et al., "Direct measurement of
hohlraum M-band preheat asymmetry at the ICF capsule", PoP 27 (2020), LLNL-JRNL-813556;
Rubery et al., "Soft x-ray power diagnostics for fusion experiments at NIF, Omega and Z
facilities", RSI (2023), OSTI 1975922. Confidence: **H** primary text, **M** quoted
secondhand, **L** our construction, *calibrated* as marked.

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| M-band definition | > 1.8 keV (Au); Mo fluorescence sensitive > 2.5 keV | Dewald 2020 §II | H |
| Laser → M-band, Au, normal incidence | ~1 % at 10¹⁴, 20 % at 10¹⁵ W/cm² | Dewald 2020 §I citing Dewald PoP 15, 072706 (2008) | M |
| Wall irradiance | inner 3 × 10¹⁴, outer 1.5 × 10¹⁵ W/cm²; inner turning point 2.5× lower density (M-band ∝ n_e² ∝ sin⁴θ) | Dewald 2020 §I | H |
| Dante M-band share | 12 % of the flux at the 290 eV peak (Fig. 1b); 15 % of the spectrum, 1/3 of it > 2.5 keV (Fig. 3) | Dewald 2020 §II–III | H |
| Setup | 10.1 × 5.75 mm Au, 3.37 mm LEHs, 360 TW / 1 MJ, inner/outer power 1/2, outer ring at z = ±3.3 mm (EHXI), θ_o ≈ 40°, R_o 4.4 mm, wall ~0.2 mm in, HDC front 0.815 mm at 4.35 ns | Dewald 2020 §II, V | H |
| Mo-sphere P2/P0 | swings, mostly > 0, +0.15 at 6 ns | Dewald 2020 Fig. 10 | H-img |
| Wall M-band P2/P0 | +0.53 (= 0.15/0.28); a pure ring at 40° gives 1.9, so F_M ≈ 2 (recirculating / laser-produced) with the LEH deficit (−25 %) | Dewald 2020 §V.2 | H |
| Smoothing | P2: 0.4 at the HDC ablation front (CCR 5.4), 0.28 at the Mo sphere; inferred M-band P2/P0 at the inner HDC surface ≤ +0.28 | Dewald 2020 Figs. 6, 12–14 | H |
| Dante 1 channels | 18 channels, 0.25–3.66 keV and 14.4 keV (table of energies, filters, mirrors, XRD cathodes); the unfold fits a Planckian to channels < 1.5 keV then Gaussian responses | Rubery 2023 Fig. 4, §II.C | H |
| Dante lines of sight | Dante 1 θ 143°, φ 274°, 37° view of the LEH; Dante 2 θ 64°, φ 350°, shallow | Rubery 2023 §II | H |

### 1. Multigroup radiosity (`spectral_solve`)

View factors are spectrally independent, so each group g solves the gray system with
its own source and albedo, `B_g = S_g + α_g H_g`. The source ηP/A minus the laser M-band
is split by Planck fractions `b_g(T)` at each element's emission brightness temperature
`(B/σ)^¼` (iterated when the totals change); the M-band is flat over its band. With
`thermalize = true` the rest of the gray albedo, `α − α_g`, is re-emitted with `b_g(T)`,
so the totals are the gray solution's whatever the group albedos (*energy as calibrated
gray*); without it, it is lost to the wall. Groups at the gray albedo reuse the model's
cached LU (or GMRES preconditioned with it; the GMRES became generic in the operator,
the gray path bit for bit); hard groups are solved jointly by unpreconditioned GMRES.

Checks (`test/spectral.jl`): Planck tails against quadrature to 1e-8 (series/tail switch
continuous to 6e-12); one group returns the gray `B`, `H` exactly (`==`); 16 groups, with
and without the M-band, in group and thermalized: Σ_g B, Σ_g H to < 1e-11 of the gray,
absorbed power to 1e-11; thermalized with α_M = 0.2 to < 1e-9; in-group solves against
dense `(I − diag α_g F) \ S_g` to 1e-10 and the thermalized coupled system against a
dense (nG)² solve to 1e-9.

### 2. Laser-produced M-band (`MBandEmission`, `nif_mband`) *engineering*

Gold elements emit `f P` of their absorbed laser power `P` in a flat 1.8–2.85 keV band
(the upper edge puts 1/3 above 2.5 keV, **L**), taken out of their Planckian source.
`f(I)` is the Dewald 2008 table in log–log; element-averaged intensity is mesh dependent
(on Dewald's hohlraum at h 0.3 mm the outer spots average 4.4e14 W/cm² against his
1.5e15), so `nif_mband()` fixes f per cone at Dewald 2020's mean wall irradiances:
outer 0.20, inner 0.042 (the 5× n_e² factor of the inners is not applied: an upper
bound). `mband_power(m, op, beams)` deposits each cone separately.

Dewald's setup (`examples/mband_dewald2020.jl`, h 0.3/0.15, static at 360 TW, 5.5 ns):
the scaled mock HYBRID-E pointing put the outer ring at z 2.58 mm; its outer aims moved
+0.53 mm along the axis to the measured 3.3 mm (θ_o 39.1°, R_o 4.24 mm against 40°,
4.4). Laser M-band 52.6 TW, 0.147 of the absorbed laser at cone fraction 1/3.

| α_M (thermalized) | F_M | capsule M-band share | capsule M-band P2/P0 | boundary P2/P0 | Dante-1 > 1.8 keV |
|---|---|---|---|---|---|
| gray | 7.63 | 0.293 | −0.036 | −0.140 | 0.300 |
| 0.7 | 3.52 | 0.151 | +0.003 | −0.074 | 0.164 |
| **0.5** | **2.15** | 0.104 | **+0.039** | −0.012 | **0.119** |
| 0.3 | 1.41 | 0.078 | +0.076 | +0.051 | 0.094 |
| 0 | 0.79 | 0.057 | +0.132 | +0.145 | 0.073 |

(cone fraction 1/3; gray capsule P2/P0 −0.114, boundary −0.298; first flight of the laser
M-band +0.505 at the capsule.) At cone fraction 0.25 (gray capsule −0.070): α_M 0.5 gives
F_M 2.05, capsule M-band P2/P0 +0.145, boundary +0.200, Dante 0.125.

- **α_M ≈ 0.5 (thermalized) is calibrated to Dewald's F_M ≈ 2**, and then Dante 1 sees
  11.9–12.5 % of its flux above 1.8 keV against Dewald's measured 12 % at 290 eV: an
  independent agreement (Dante's M-band share was not fitted). (With the laser spots
  thermal at the local Tr, the 2026-10-24 default, F_M ≈ 2 takes α_M ≈ 0.55 and Dante
  sees 11.4–12.0 %: "Laser-spot spectrum" below.)
- The M-band is pole-hot relative to the thermal drive by +0.15 to +0.21 in P2/P0 at the
  capsule (M-band minus gray), against Dewald's ≈ +0.21 at the ablation front (0.53 ×
  0.4, **L**). The absolute M-band P2 here (+0.04 / +0.15) sits lower because this mock
  pointing's thermal drive is waist-hot (−0.11 / −0.07) and the recirculating M-band
  (2/3 of it) carries that shape, where Dewald takes it as uniform. Likewise the
  boundary P2/P0 (−0.01 / +0.20) against his +0.53, which assumes a uniform background.
- Not modeled: the outer-spot motion and inner-beam clipping in the outer bubble that
  Dewald invokes for the late swing; angle dependence (sin⁴θ) of the conversion.

### 3. Synthetic Dante (`DanteView`, `dante1`, `dante`)

Parallel rays (400 × 400) across the LEH give each element's visible projected area,
partial views included; a ray through a closing lip is weighted by its open fraction.
`Tr = (π I / (σ A cos θ_v))^¼` over the full initial LEH area (Dante's convention),
`Tr_open` over the open area. Channels: Gaussians at Rubery's 18 Dante 1 energies, FWHM
0.3 E (**L**, widths unpublished; since replaced by the hardware responses, "Dante channel
responses and the unfold" below); `color_temperature` fits a Planckian to the channels
below 1.5 keV. Checks: a blackbody cavity reads T₀ to 5e-4 (after the caps' inscribed
24-gon opening, 1.1 % below the equal-area hole on the coarse test mesh: on coarse
meshes the opening is smaller than π r²); lips at open fractions 0, ½, 1; a Planckian's
channel signals unfold to its T to 1e-6.

N221204 (`examples/n221204.jl` full model, standard mesh, default energetics):
| | peak | at 8.5 ns |
|---|---|---|
| capsule Tr (mean incident flux) | 301.8 eV at 5.50 ns | 270.3 |
| Dante 1, initial LEH area (Dante's convention) | 285.9 at 6.25 ns | 257.4 |
| Dante 1, open area (lips at their open fraction) | 306.2 at 6.25 ns | 278.1 |
| measured Dante 1 (Rubery 2024; Pak 2024 Fig. 1) | 313 ± 2 at ≈8.4 ns | ≈310 |

- Dante 1 sees only the wall (no capsule, no cap) at 37°. Over the open area it reads
  7 eV above the capsule at the Dante peak, the sign of PRE 109's "internal ~13 eV below
  Dante" (**L**: the simulated internal temperature is not exactly our capsule mean).
  Over the initial area it reads 13 eV *below* the capsule: the lips (closing to
  0.857 r₀) are dark to Dante here, while the real closing plasma radiates. Which of the
  two the measurement corresponds to depends on the LEH plasma's emission, which this
  model does not have; the two bracket the measured 313 only at the open-area end, and
  7 eV low there.
- The larger miss is the time history: the model's drive peaks at 5.5–6.25 ns and falls
  ~30 eV by 8.5 ns, where Dante's keeps rising to its 8.4 ns peak. The default
  energetics (Au wall loss, constant conversion 0.75, ablator at its initial radius, no
  backscatter) and the bubble's growing inner-beam absorption lose drive through the
  peak; this is the first absolute check of item 16, and it points at the late-time
  energetics (wall loss vs time, conversion rising in time, the `energetics=new` options)
  rather than at the line of sight. (Resolved in "N221204 late-time drive" below: the
  wall loss was evaluated at the time since the pulse start; these numbers are the old
  evaluation.)

## Capsule drive hand-over and CompactLES stage 0 (2026-10-24) *derived*

ROADMAP Next 17, stage 0. Sources: `reference/litreview/hydro_coupling.md` (all
literature numbers below are from there, with its confidence marks). Confidence here:
**L** our construction, *calibrated* as marked.

**Hand-over (`capsule_drive`).** The drive is the history's own `moments` (Legendre
about +z, about the capsule centre, P2 > 0 pole-hot), interpolated linearly in aₗ and
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
| shell pressure at the front | P_abl within 1 % from 0.75 ns | — (by construction) |
| peak DT velocity | 315 km/s at 8.7 ns | 380–400 km/s |
| ablator remaining at peak velocity | 19.6 % | 5–6.6 % |
| ablator + DT KE at peak velocity | 36 kJ | 19–33 kJ |
| stagnation | ≈ 9.3 ns | (bang 9.53, not a prediction) |

Reading (**L**): with the pressure right, the rocket is too heavy: Lindl's ṁ (a Nova
fit for low-albedo ablators) removes 80 % of the HDC ablator where HYBRID-E removes
~94 %. Calibrating ṁ to the mass remaining fixes v through the rocket equation; it is
one more knob and the obvious next step, and the KE comparison is loose (our KE counts
the remaining ablator). Numerics: the run takes ~30 min on 4 threads (timestep set by
the hot driver gas, CFL 0.1 after the automatic retries in the foot); a 2 µm grid
(layer < 3 cells) fails in the foot.

**Update: robustness and the ṁ calibration (2026-10-24, CompactLES `briareus-drive`;
drive from Briareus 785684d, cached).** Two defects of the first run, both fixed in
CompactLES. (1) The foot failure was the driver gas, not the layer. Held at P_abl, the
1 mg/cc He is a 2–7 keV ideal gas against a 0.05 eV ablator. The contact cells had
negative internal energy from 0.1 ns and collapsed the timestep near 0.5 ns, at any
CFL. A 10 mg/cc driver (≈ 1 % of the ablator's mass) runs the foot at CFL 0.4 with
no retry. The whole implosion is now one `run!` call of 5–8 min on one thread, with
1–2 automatic retries when the first shock reaches the origin (7.3 ns). (2) The
19.6 % was mostly ablator stranded outside the front. A partly converted layer gets no
push in the uniform-pressure driver, so the accelerating shell left it behind (22 % of
the ablator at 6 ns). It then pulled the front marker off the shell, and ablation of
the shell stopped. Now ablator beyond the front is converted first, and the piston
acts only at and beyond the front. Results (1 µm grid; ablator remaining = ablator
inside the front at peak DT velocity), *calibrated*:

| ṁ prefactor A (ṁ = A T³) | v_DT peak | ablator remaining | KE (ablator + DT inside front) |
|---|---|---|---|
| 2.6e5 | 320 km/s | 7.4 % | 26.5 kJ |
| **2.85e5 (adopted)** | **336 km/s** | **5.7 %** | **25.8 kJ** |
| 3.0e5 (Lindl Eq. 6-3) | 346 km/s | 4.8 % | 25.1 kJ |
| 3.8e5 | 406 km/s | 0.9 % | 19.2 kJ |

Reading (**L**): 5–6.6 % remaining needs A = 2.7–3.0e5, so Lindl's 3e5 T³ was right
and the earlier "too little ablation" came from the stranding. The T³ exponent is
kept: Lindl's I^(3/4) is T³, and P/ṁ ∝ T^0.5 with P ∝ T^3.3 would give T^2.8, a
difference the 290–297 eV-dominated ablation does not resolve. The velocity check
fails without tuning: 328–344 km/s across the band against 380–400. Reaching 380–400
needs A ≈ 3.5–3.75e5 and leaves 1.5–2.5 %. Lindl's ideal rocket gives 395 km/s at this
mass ratio. The gap is in P/ṁ, so the likely causes are the ideal-gas carbon adiabat,
the 145 Mbar peak and the gray drive, not ṁ. KE 25–26 kJ is inside 19–33 kJ. The grid
is not converged: at 2 µm the same A gives 298 km/s and 7.3 %, and the 0.5 µm run fails
at 5.5 ns. Read v to ±10 %. The calibration belongs to the drive. Briareus 3f25f29
(late-time wall-loss history) peaks at 285 eV, and with P₀ recalibrated to 5.53 Mbar
the same A gives 307 km/s, 8.2 % and 25.6 kJ (the CompactLES example, end to end).
r–θ run (aₗ to l = 8, 2 µm × 64 cells, annulus from 0.25 mm, stopped at 7.6 ns).
The pole-hot foot (a₂/a₀ +0.03…+0.15 to 4.7 ns) makes the shell oblate, to −0.45 %
P2/P0 at 5.5–6 ns. The pole-cold peak (−0.09) reverses it, and the shell is prolate
from ≈ 7.1 ns (+0.5 % at 7.6). This is a swing, but not yet Masse's G(t,τ) inversion,
because the drive's own P2 changes sign at 4.7 ns. Telling them apart needs impulse
runs to stagnation (**L**).

**Update: every m, and the r–z and 3-D metrics (2026-10-24; CompactLES
`briareus-drive`).** `capsule_drive` now carries `ylm`, the history's real Y_lm
(`drive_history(...; ylm = true)`) as Schmidt semi-normalized ratios
bₗₘ = √((2l+1)/4π) cₗₘ/a₀, so bₗ₀ = aₗ/a₀ (the m = 0 rows are set to `modes`; the
history's Y_lm quadrature is one sub-triangle level coarser, < 1e-3 of a₀ apart) and
√Σₘ bₗₘ² = `mode_amplitude`. Tested: a known field's ratios recovered to its
element-constant quadrature, round trip, and a 45°-rotated three-quad pattern on a
45°-symmetric mesh rotates (cₗₘ, cₗ,₋ₘ) by mα to 1e-7. CompactLES's front search and
layer normalization are now sums over equal-solid-angle bins about the capsule centre
with dV/R² weights (per steradian, so the per-line form is reproduced to 1e-14 in
1-D), which runs on r–z (axis) and 3-D Cartesian grids, and on refined levels through
a level-wide hook before each level's right-hand side. Validation (**L**), flat
250 eV drive, P₀ and ṁ not retuned, 2-cell layer: r–z 10 µm against 1-D 10 µm to
1.5 ns, front within 0.45 µm, DT shell within 0.15 µm, ablator left equal to 4 digits,
momentum within 1 % (a tiled 10 µm band likewise to 0.6 ns, then failed); a b₂₀ = 0.2
drive gives the same front P2 on r–z as on r–θ to 10–15 % (−0.55 vs −0.50 % at 1.5
ns); a 3-D octant at 20 µm reproduces 1-D (front +0.5…+1.7 µm, ablator to 1e-4), and
its response to b₂₂ = 0.2 equals its response to b₂₀ = 0.2 to 1–4 % (−0.123 vs
−0.122 % of R_f at 0.6 ns), with ≤ 0.014 % in Y₄₂ (cubic-grid coupling) and ≤ 0.003 %
elsewhere. Limits met: a band's coarse level takes the timestep of its covered nodes
(8× the uniform grid's steps), and tiled levels failed on this problem; a 1 µm 3-D
shell (~10⁹ cells uniform) waits on those.

**Update: the full-model drive and an HDC ablation-rate source (2026-10-25; Briareus
c54e366, CompactLES `briareus-drive-2`).** The drive is now `examples/n221204.jl`'s full
model (+bubble: `nif_uranium()`, Young's albedo on the K24 schedule, release 1.0,
calibrated two-zone CBET, measured backscatter, 11 % degradation, wall, lips, bubble),
coarse mesh 0.8/0.4, steps of 0.1 ns from 0.2 ns to bang (9.53 ns), l ≤ 6 with every m.
The capsule Tr peaks at 306 eV at 8.0 ns (297 before) and holds 282 eV at 9.0 ns; P2/P0
+0.05 near 3.2 ns, −0.18 at 5.2 ns, +0.09 over 7.2–7.7 ns. Stage 0, 1 µm, ṁ = 2.85e5 T³,
P₀ recalibrated to 145 Mbar (*calibrated* as before; the checks are PRE 109 Table I,
HYDRA, N221204: v_imp 383 km/s, ablator remaining 5.7 %, bang 9.61 ns simulated,
9.53 ± 0.07 measured):

| | stage 0 (new drive) | check |
|---|---|---|
| peak DT velocity | 339.5 km/s at 9.26 ns | 383 km/s |
| ablator remaining at peak velocity | 4.97 % | 5.7 % |
| KE (ablator + DT inside the front) | 24.0 kJ | 19–33 kJ |
| stagnation (min. DT mean radius, peak p(0)) | 10.22 / 10.20 ns | bang 9.53–9.61 ns |

Reading (**L**): the stronger, later peak moves peak velocity to 9.26 ns and stagnation
0.6–0.7 ns past bang; no burn or radiative loss, so stagnation is not bang, but a
0.6 ns lag is a drive-or-closure timing error, not a definition. Velocity is unchanged
within the ±10 % grid band (336 → 340). Ablation-rate source independent of v_imp
(cite-only): Olson et al., PoP 18, 032706 (2011) (OSTI 1114987, abstract: HDC 2–9
mg cm⁻² ns⁻¹ at 160–260 eV, ≈ T³, planar Omega halfraum burnthrough, **H**); the HDC
points of SAND2007-5659C (OSTI 1147367, IFSA07; 10 points read from the slide's vector
plot, 175–255 eV, e.g. 8.3 at 245 ± 5 eV) give ṁ = A T³ with A = 4.9–5.2×10⁵
g cm⁻² s⁻¹ (free exponent 3.8), ±11 % from Dante's ±3.5 % (**M**); SAND2011-8633C (OSTI
1141159) uses 0.415 T³ mg cm⁻² ns⁻¹ (4.15×10⁵) for HDC and Ge-CH together. Lindl's
3×10⁵ is 1.4–1.7× below the HDC data. At those rates this shell over-ablates (stage 0,
1 µm, same drive): 4.15e5 leaves 0 % and 5.0e5 fails at 7.33 ns with 0.11 % left. The
grid study at 2.85e5 fails too: nx 700 and 2000 end on `dt_collapse`, only nx 1400
completes. A = 2.85e5 is kept as an *effective* rate, *calibrated* to the ablator remaining
(5–6.6 %, PRE 109 Table I / HYDRA), so that number is a fit target, not a check; v_imp,
KE and bang time remain checks (decided 2026-10-26). Replacing it with a measured rate
is CompactLES ROADMAP A22 (rocket-relation pressure tied to ṁ, real carbon EOS, grid
robustness).

**Where the stagnation lag comes from (2026-10-27, L).** Eight variants of the run above
(same drive, 1 µm, nx 1400, one closure changed each), timed against K24's HYDRA radius
history (Fig. 1, the ablator–fuel interface; lag = stage-0 time minus K24's at the same
R/R₀, DT mean radius; keyhole N230110 shock timing, Ruof 2026 Table II, as a check only):

| case (change) | shock into ice | v_peak | ablator left | stagnation | lag at R/R₀ 0.95 / 0.85 / 0.50 / 0.20 |
|---|---|---|---|---|---|
| baseline | 4.80 ns | 339.5 km/s | 4.97 % | 10.22 ns | +0.98 / +0.68 / +0.57 / +0.55 |
| p30 (foot P ∝ T^3.0, peak kept) | 4.48 | 324.1 | 5.58 | 10.22 | +0.79 / +0.54 / +0.52 / +0.54 |
| foot (picket/foot Tr × 1.15 to 3 ns) | 4.16 | 332.7 | 4.64 | 10.12 | +0.65 / +0.46 / +0.44 / +0.44 |
| g3p30 (p30 + γ_C = 3, stiff carbon) | 3.70 | 301.4 | 7.33 | 10.52 | +0.67 / +0.47 / +0.67 / +0.76 |
| footg3 (foot + γ_C = 3) | 3.42 | 300.9 | 9.47 | 10.44 | +0.48 / +0.38 / +0.62 / +0.72 |
| pk160 (P × 1.10) | 4.72 | 344.0 | 5.72 | 10.08 | +0.87 / +0.56 / +0.45 / +0.42 |
| sh03 (whole drive 0.3 ns earlier) | 4.56 | 339.4 | 5.06 | 9.92 | +0.67 / +0.37 / +0.27 / +0.25 |
| r33 (ṁ = 3.3e5 T³) | 4.86 | 371.1 | 2.41 | 9.84 | +0.94 / +0.60 / +0.37 / +0.28 |
| r35 (ṁ = 3.5e5 T³) | 4.86 | 386.1 | 1.45 | 9.70 | +0.92 / +0.56 / +0.29 / +0.16 |

Keyhole check: first shock into the ice 3.86–3.91 ns. Reading: the early lag (R/R₀ ≥ 0.85)
is the first shock: every change that brings the shock to the keyhole timing (stiff
carbon, stronger foot) cuts it. But the stiff-carbon shells then accelerate more slowly
(300 km/s, 7–9 % left) and the lag grows again to +0.7–0.8 ns by R/R₀ 0.2: fixing shock
timing alone moves stagnation later, not earlier. Only burning more mass closes the late
lag (r35: 386 km/s, lag +0.16, but 1.45 % left against 5–6.6 %), or a drive 0.3 ns earlier
(sh03), which Dante rules out (±6 eV over 4.5–7.5 ns). No single closure change matches
shock timing, v_imp, remaining mass and stagnation together: in a 1-D rocket model
velocity and remaining mass are tied by the rocket relation through ṁ and the exhaust
velocity (P/ṁ), and with ideal-gas carbon no ṁ gives HYDRA's 383 km/s *and* 5.7 % left.
This agrees with the user's reading that the lag is rooted in a 1-D rocket model tuned to
3-D results (ṁ fitted to remaining mass), and points at A22's items 1–2 (P tied to ṁ, real
carbon EOS) rather than at the drive. Side note: K24's late slope (335–345 km/s) and Ruof
2026's 2DConA 355 ± 5 km/s (N240304) are near stage 0's 340, so the 383 km/s target itself
is model-dependent.

## Large radiation models (2026-10-24) *derived*

ROADMAP 14. Goal: make 1e4–5e4-element radiation models (the fine hohlraum preset,
remeshed CAD) practical on a laptop (Ryzen 9 4900HS, 4 threads of 8 used, 16 GB, another
job sharing the machine). Workload, as a design session runs it: the GUI NIF hohlraum
(10.1 × 2.875 mm, LEH 1.55, capsule 0.9) with a `WallLossAlbedo` wall, NIF quads;
deposition operator, view factors, closure, a temperature-dependent `solve` at 2 ns, a
6-step `drive_history`, a `Scenario` (62-row moment response), an optimizer gradient
and `ev.flux`. Script: the session's scratchpad `prof.jl` (one process per mesh,
`BRIAREUS_NO_WORKLOAD=true`, so first calls include compilation: "solve again" is the
warm figure).

### Where the time and memory went (before)

| N (h / capsule_h) | 3,712 (0.6/0.09) | 10,484 (0.22/0.07) |
|---|---|---|
| deposition operator | 2.4 s | 3.5 s |
| view factors | 1.2 s | 14.3 s |
| closure | 0.3 s | 2.1 s |
| solve (first, LU) / again | 0.9 / 0.4 s | 9.1 / 3.0 s |
| drive_history, 6 steps | 4.3 s | 20.9 s |
| Scenario (incl. a copy of F) | 3.4 s | 13.2 s |
| gradient | 0.02 s | 0.02 s |
| ev.flux (full response) | 0.9 s | 8.9 s |
| max RSS | 1.4 GB | 6.6 GB |

Dense Float64 holds F and up to two LU factors, 24N² bytes: 7.5 GB at 17.6k, 20 GB at
29k, 38 GB at 40k; one LU is 2N³/3 flops (≈40 s at 17.6k on 4 threads, ≈7 min at 40k).
The fine preset did not fit next to the other job, and 40k does not fit at all.

### What changed

1. **Float32 storage** (`Float32Matrix`, `view_factors(...; storage = Float32)`): Float32
   data read as Float64; every product accumulates in Float64, so it is exactly the
   product with `Float64.(F.data)`. Assembled directly in Float32 (no Float64 peak),
   bit-identical to `Float32.(view_factors(t))`. Balancing (`enforce_closure!`) runs in
   Float64 on the stored values and rounds once: row sums to ~1.5e-8 (Float64 ~1e-12),
   negligible beside the 1e-3 discretization error it removes.
2. **Factorization-free solves** (`RadiationModel(...; solver = :iterative)`): restarted
   block GMRES (`_gmres_block`), no preconditioner, 1e-12 relative residual like the
   LU-preconditioned path; the columns of a block (dual value and partials, spectral
   groups, the 62 moment rows) iterate in lockstep, each with its own basis, sharing one
   product with F per iteration; bases kept under 256 MB (columns batched). Iterations
   at 17.6k: 17 for a cold solve, fewer warm (a history step starts from the last
   radiosity), 28 for the transposed moment rows. `ev.flux` on such a model is one solve
   instead of the full observed × N response (still available, a solve per row).
3. **Defaults**: above `Briareus.LARGE_MODEL[]` = 10,000 elements `RadiationModel`
   computes Float32 view factors and solves iteratively; below it nothing changes (dense
   Float64 + LU, bit for bit). `storage`/`solver` override; a given `F` keeps its storage.
   `Scenario(...; model)` reuses a model (no copy of F).
4. **Kernels**: vector products threaded over row/column blocks at 28–33 GB/s (memory
   bound); a fused 4-column kernel for blocks of ≤ 48 columns (3–7× the widen-to-Float64
   + BLAS panels at k = 2–16; the panels win past ~50 columns). 17.6k: F x 43 ms,
   Fᵀx 37 ms, k = 2: 64/50 ms, k = 16: 280/240 ms, k = 62: 880/890 ms (panels).
5. **View-factor assembly** (both storages, bit for bit): the near-pair kernel's result
   was inferred `Any` (the CAD cell list is untyped), so every pair boxed its factor and
   dispatched the store dynamically: 13.7M allocations at 5k, and GC pressure across
   threads. A type assertion removed them: 10.5k 14.5 → 12.1 s (allocation 2.9 → 0.7
   GB), 17.6k 56 → 31 s (9.1 → 1.5 GB).

### After (Float32 + iterative; same script)

| N | 10,484 | 17,632 | 29,124 |
|---|---|---|---|
| deposition operator | 6.5 s | 6.6 s | 9.8 s |
| view factors | 13.0 s | 31.4 s | 101.6 s (140 s before the fix in 5) |
| closure | 1.6 s | 3.0 s | 6.2 s |
| solve again (T-dep albedo, warm) | 1.0 s | 2.5 s | 7.4 s |
| drive_history, 6 steps | 8.9 s | 14.7 s | 31.4 s |
| Scenario (62-row response) | 11.2 s | 17.2 s | 39.0 s |
| gradient | 0.02 s | 0.03 s | 0.04 s |
| ev.flux (one solve) | 0.37 s | 0.64 s | 1.4 s |
| max RSS | 1.3 GB | 2.3 GB | 4.5–5.3 GB |

(The before runs loaded a package image with the precompile workload, the after runs one
without it, so first calls here include compilation: ~3 s of the deposition operator,
~8 s of the first solve.) 40k extrapolates
to 6.4 GB for F, ~190 s of view factors, ~14 s per warm solve.

Accuracy of the Float32 path against the dense Float64 path at 10,484 elements (same
mesh, same process otherwise): capsule a₀ 1.3e-10 relative, P2/P0 and P4/P0 2e-10 and
1e-10 absolute (solve and Scenario), drive history moments ≤ 2e-10, Tr ≤ 2e-9 eV, H
1.2e-8 relative at worst (an LEH element). The Float64 iterative path agrees with LU to
the GMRES residual (~1e-12). Tests: `test/large_models.jl` (forced on at ~1k elements:
products, block GMRES, solve/flux_response/duals/drive_history/Scenario gradient/
`spectral_solve`/keyframe interpolation against the dense path).

### Tried and dropped

- **Two-grid preconditioner** (clusters of 16 or 4 elements by coordinate bisection,
  Galerkin coarse view factors, coarse LU per albedo): 17 → 12 (11) iterations per solve,
  28 → 26 (24) for the moment rows, no gain in time. The slow part of the spectrum is not
  the smooth one the clusters capture.
- **Row-blocked stores in the assembly** (rows written 16 at a time): the profile's
  samples on the stores were skid from the divisions; no change in time.

### Not done: GPU and hierarchical radiosity

- **GPU** (RTX 2060 Max-Q, 6 GB, ~1/32-rate FP64): Float32 F fits to ~35k elements; a
  mat-vec with Float64 accumulation needs 2 flops per 4-byte entry, ~130 GFLOP/s at the
  card's ~260 GB/s, about its FP64 rate, so ~6–8× the CPU's 30 GB/s is the ceiling. The
  solves are ~30 % of the workload at 17–29k (Scenario response, history, solve), so a
  CUDA extension would save ≲ 25 % end to end; view-factor assembly (BVH rays and exact
  contour kernels, CPU) remains the largest cost. Not worth CUDA.jl yet.
- **Hierarchical (H-matrix) view factors** are the real answer past ~30k: assembly and
  storage are N², and both grow 4× per halving of h. Design: cluster tree by coordinate
  bisection (as tried above); admissible cluster pairs (distance > η × diameter, and
  visibility tested by a few rays per pair and found all-or-nothing) stored as low-rank
  blocks (ACA on the far-pair kernel, rank 1–8), the rest (near pairs, partly occluded
  pairs straddling the capsule shadow or LEH rims) dense; O(N log N) memory and
  assembly. It needs only `mul!` with F and Fᵀ and row sums: the iterative path here
  already uses nothing else, so an H-matrix type plugs in as another
  `AbstractMatrix{Float64}` (balancing is a diagonal scaling s F s). Reciprocity per
  block (G = AF symmetric) halves it again. Built: next section, "Hierarchical view
  factors".

## Hierarchical view factors (2026-10-24) *derived*

ROADMAP 14, the step after "Large radiation models": view factors whose assembly and
storage do not grow as N², as another matrix type (`HMatrix`, `src/hmatrix.jl`) that the
iterative path uses through products only. `view_factors(t; storage = :hmatrix)`,
`RadiationModel(t; storage = :hmatrix)`, and `RadiationModel`'s default above
`Briareus.HMATRIX_MODEL[]` = 15,000 elements. Workload and machine as in "Large radiation
models" (GUI NIF hohlraum 10.1 × 2.875 mm, LEH 1.55, capsule 0.9, `WallLossAlbedo` wall,
NIF quads aimed at z = ±4; 4 threads, another job sharing the machine); script: the
session's scratchpad `hmat/bench.jl` (warm: a small target through the same path first).

### Design

- **What is stored.** The exchange matrix G = A F (symmetric), once per unordered pair:
  a cluster tree by median bisection of the centroids along the longest box axis (leaves
  ≤ 48 elements), blocks of the upper triangle only (Gᵀ blocks implied), `F = A⁻¹ S G S`
  with a balancing diagonal S. Reciprocity is exact; `enforce_closure!` (Sinkhorn, one
  product per iteration) updates S only. Dense blocks (near pairs, leaf pairs that do
  not compress) are Float32 with the entries `view_factors` computes (`_vf_pair`, now
  split from the row loop; the dense path is bit for bit unchanged). Low-rank blocks are
  U Vᵀ (Float32 factors), and "scattered" dense blocks on index lists hold the rows a
  facing test cuts.
- **Admissibility.** Cluster boxes (centroids padded by element size) at least η × the
  larger diameter apart (η = 1). Then, per block:
  1. *Facing* (box test, no rays): rows whose element faces no corner of the other
     cluster's centroid box (beyond the facing test's slack) are zero; a block with
     none is zero. Small blocks (≤ 16 leaf²) take the facing test pair by pair: a block
     facing itself in under a quarter of its rows × columns (a capsule cluster at the
     terminator seen by a wall cluster: 93–100% zeros, ACA finds zero rows and stops)
     is stored dense on those.
  2. *Shaft culling* (Haines & Wallace 1994): the convex hull of the two boxes (the union
     box's faces plus the planes through one axis-parallel edge of each box) is searched
     in the BVH for facets that could occlude a segment between the clusters. Not
     occluders: facets of a curved surface (cylinder, sphere, cone) both clusters lie
     inside (`_vf_visible`'s own rule); facets of a convex closed component with
     outward normals (a capsule; checked vertex against facet, ≤ 2e9 pairs) when one
     cluster lies on it with one facet per element (a segment leaving a face of a convex
     polyhedron outward does not come back); facets whose plane has both clusters
     strictly in front (by 2ε) or in that plane with its normal (caps and LEH discs,
     corner clusters with the wall in front of the cap). An empty shaft makes every
     facing pair wholly visible, so its entries skip the ray (`_vf_pair_free`: pairs
     facing beyond the slack; slack-admitted pairs still take theirs).
  3. *Hidden by a convex part*: if a convex component's inscribed ball (about the vertex
     mean) blocks every segment between the boxes — each corner of one box inside the
     ball's shadow cone from each corner of the other, beyond that corner's polar plane;
     the shadow beyond the polar plane is convex, which carries the corner test to the
     boxes — and both boxes lie outside its circumscribed ball, the block is zero.
  4. An occupied shaft refines to leaf pairs, which cast every pair's ray: wholly hidden
     → zero, wholly visible → cross approximation, partly seen (a shadow edge,
     `:sampled` fractions) → dense on the known mask.
  Visibility is never decided by sampling: rays between 7 representatives per cluster
  plus 32 random pairs (the first design) passed cap-to-cap blocks 0.5% occluded by the
  capsule (6 of 1280 entries, each a full entry), a 3% block error, and the sampled
  check of the cross approximation passed them too. Sampled hidden/visible decisions
  on occupied blocks of up to 4 leaves a side (with exact leaves) were 1.25× faster at
  28.6k but left element H errors of 5e-4 (moments 1e-7): dropped.
- **Low rank.** ACA with partial pivoting on the exact entries, stopped at a relative
  Frobenius estimate `tol` (1e-5), given up early when the terms' decay extrapolated
  (over the last 4) misses `3 tol` at `maxrank` (32, and ≤ 3mn/4(m+n) so it pays), QR +
  SVD recompression to `tol`, then a sampled check ((m+n)/2 entries, error ≤ 10 tol of
  the block). An empty result never proves a zero block (it did: capsule clusters whose
  sampled rows were all behind; 5.8e-4 of ‖F‖ dropped). Rows facing the other box
  wholly get the approximation and the few cut by the facing test (≤ half) stay dense.
- **Products.** Blocks grouped in row order into one run per thread of about equal
  cost, each into its own buffer, summed in group order (deterministic for a thread
  count); 16 columns at a time; fused kernels (one pass over a dense block for its
  block and its transpose). Dual vectors go through the values-and-partials matrix.
  Keyframes of hierarchical models interpolate as an operator (`_BlendedExchange`, the
  two keyframe operators with per-row weights and its own balancing), since their blocks
  differ.

### Accuracy (defaults tol 1e-5, η 1, leaf 48)

Against the dense Float32 path on the same mesh (balanced, T-dependent albedo solve at
2 ns), and at 28.6k against an `HMatrix` at tol 1e-7 (dense Float32 there needs 3.3 GB
of F next to the other job):

| N | capsule a₀ | P2/P0 | P4/P0 | H per element, max (median) | row sums | moment response |
|---|---|---|---|---|---|---|
| 10,484 | 3.1e-8 | 4.2e-8 | 2.6e-8 | 8.6e-5 (3.6e-8) | 1.7e-8 | 6.6e-5 |
| 17,632 | 1.4e-9 | 2.7e-10 | 2.8e-8 | 8.6e-6 (4.7e-8) | 1.1e-8 | 4.6e-6 |
| 28,628 (vs tol 1e-7) | 5.0e-9 | 9e-9 | 4.5e-9 | 3.5e-6 (6.3e-8) | — | 7.1e-6 |

(a₀ relative, P-ratios absolute, the rest relative.) tol 1e-4 moves the moments by
2–9e-7 for 3–5% less time; 1e-6 by ≤ 1e-8 for ~30% more. On small targets with rims,
bores, a waist gap, corner fillets, remeshed CAD (48 facets per capsule element) and
`:sampled` visibility (`test/hmatrix.jl`, leaves of 16), every G entry matches
`view_factors` to < 1e-5 of the largest. Non-public CAD was not run.

### Time and memory (warm view-factor assembly; dense rows from "Large radiation models")

| N | 10,484 | 17,632 | 28,628 | 49,800 |
|---|---|---|---|---|
| Float32 dense: view factors | 13.4 s | 31.4 s | 101.6 s (29.1k) | ~300 s (est.) |
| Float32 dense: F | 419 MB | 1,186 MB | 3,278 MB | 9,920 MB (does not fit) |
| HMatrix: view factors | 7.6 s | 17.3 s | 34.8 s | 82.2 s |
| HMatrix: stored | 122 MB | 240 MB | 432 MB | 1,004 MB |
| HMatrix: max RSS, whole run | 0.85 GB | 1.00 GB | 1.36 GB | 2.26 GB |
| F x (HMatrix) | 7 ms | 19 ms | 34 ms | 81 ms |
| closure (Sinkhorn) | 0.5 s | 1.2 s | 2.1 s | 5.5 s |
| T-dep. solve, warm | 0.5 s | 0.7 s | 0.8 s | 1.6 s |
| 5-row moment response | 0.4 s | 1.1 s | 2.0 s | 5.7 s |
| mean rank of low-rank blocks | 7.5 | 6.9 | 6.0 | 5.6 |
| dense entries (incl. scattered) | 27M | 50M | 81M | 199M |

Assembly grows as N^1.4–1.6 and storage as N^1.2–1.5, not N log N: what remains dense
is the capsule's shadow edges (wall-to-wall leaf pairs partly hidden: 1.4k → 24k mixed
leaf pairs from 10k to 50k, the N^1.5 of a 3-D boundary in the 4-D pair space), leaf
pairs that do not compress (corner clusters of cap and wall seen past the capsule's
terminator; 3k → 100k ACA give-ups), and the near field. Wholly hidden blocks are rare
(the ball certificate fires on 4 → 2.7k blocks): the capsule's umbra is small against
the cavity. Profile at 28.6k (before the ray-free entries): rays were half the time;
with them, ACA ~45% (half of it on blocks that give up), dense leaves ~30%, leaf masks
~10%.

### Defaults

> **Superseded in part (2026-10-25).** The default `tol` is now **1e-6** and the cluster
> tree splits by surface first ("The H-matrix's N^1.5 remainder (2026-10-25)"); the
> thresholds below (`HMATRIX_MODEL[]` 15,000, Float32 from 10,000) are unchanged.

Above 15,000 elements `RadiationModel` takes the HMatrix (2× faster and 5× smaller at
17.6k, moments within 1e-8); 10–15k stays Float32 (1.8× faster and 3.4× smaller at
10.5k too, but the element-wise H error there reaches 8.6e-5 at an LEH element). Small
models are unchanged bit for bit (dense Float64 + LU), and so is the dense path itself.
The precompile workload runs a hierarchical model (leaves of 16) on its coarse hohlraum.

### Open (at the time; see "The H-matrix's N^1.5 remainder" below)

- The N^1.5 terms: low rank of the shadow-edge blocks as kernel × mask (the mask is the
  costly part: a ray per pair), and a better split of corner clusters (cap + wall in one
  leaf) so their blocks compress.
- ACA give-ups cost ~20% of assembly; a cheaper predictor (the facing fraction, the
  shaft) could skip hopeless attempts.
- GUI `radiation_size` should budget an HMatrix (~0.2 GB per 10k elements here).

## The H-matrix's N^1.5 remainder (2026-10-25) *derived*

ROADMAP 14, the three open items above. Same workload, machine and script as
"Hierarchical view factors" (the session's scratchpad `hmat2/`: `attr.jl` time
attribution, `bench2.jl`, `vf.jl` interleaved timings, `rowchk.jl`); several other Julia
jobs shared the machine, so before/after times below are from runs interleaved in one
sitting (the old code from a worktree at 9e87912).

### Where the time went (before, 28.6k, 155 CPU-s inside block assembly on 4 threads)

| Part | CPU-s | share |
|---|---|---|
| near-field dense leaves (no ray / rays) | 18.1 / 3.9 | 14 % |
| ACA that succeeded (free blocks; occupied leaves, with rays again) | 34.5 / 9.0 | 28 % |
| ACA that gave up (non-leaf / leaf) | 16.8 / 7.3 | 16 % |
| dense leaves after a give-up | 10.0 | 6 % |
| leaf visibility masks (20M rays) | 14.3 | 9 % |
| mixed (shadow-edge) leaves, kernel on the mask | 7.0 | 5 % |
| rows cut by the facing test (split, sparse) | 19.3 | 12 % |
| facing tests, shafts | 14.7 | 9 % |

By cluster pair, the give-ups were capsule–wall (29.6k of 46.5k), cap/LEH–capsule and
mixed clusters (capsule + wall, cap + wall). The facing fraction predicts them well:
blocks with fewer than half their rows or columns wholly facing gave up 45.6k times
(23.6 s) and compressed 3.0k times (5.3 s). The shadow-edge leaves (11.8k) were mostly
wall–wall across the capsule. Why the capsule's ball certificate (below) could not
apply at first: 7k shafts per 10k-element model held capsule facets *and* a cluster
whose box reached into the capsule, because the median bisection put capsule and wall
elements in one cluster.

### What changed

1. **Clusters by surface.** `_hcluster_tree` splits a cluster whose elements lie on
   several surfaces (a component's curved surface, else the component; caps and their
   coplanar LEH discs share one) by surface first, the largest from the rest, then
   bisects as before. No cluster mixes capsule and wall or cap and wall at a corner.
   This is the "corner cluster" fix and most of the gain: far candidates 177k → 57k,
   low-rank blocks 42k → 10k, give-ups 46.5k → 15k at 28.6k. It is also more accurate
   element by element: the Frobenius tolerance of a mixed block let its small rows
   (capsule elements next to wall elements) be wrong relative to themselves. A capsule
   row at 28.6k had 2 entries dropped (G = 0 against 1.7e-8, 1.7 % of the row's
   largest), its row sum 3.7e-5 low, in the old tree; 4e-11 after.
2. **Shadow edges decided by the capsule's balls** (`_shaft_occluder`, `_hvis_ball`,
   `_hmask`). The shaft search now tells apart an empty shaft, one holding facets of a
   single convex part only (the capsule) with both clusters outside its circumscribed
   ball, and anything else. In the second case every pair of the block is decided
   exactly, without sampling: a segment (between the lifted ray ends of `_vf_visible`)
   that misses the circumscribed ball hits no facet of the part, and one through the
   inscribed ball crosses its closed surface (margins 1e-9 relative, far above the ray
   test's 1e-10 barycentric slack). Only pairs between the two balls, or admitted by
   the facing slack, take the ray. The mask (Int8, blocks up to 64 leaf²) makes the
   block zero, low rank with entries off the mask (no more rays), or refined with the
   mask passed down to dense leaves. Children of an empty shaft inherit it, so no shaft
   is searched twice, and a BVH node whose facets are all of the certified part (or of a
   convex part a cluster lies on) is skipped whole. Leaf masks by rays now return the
   same mask type, so a wholly visible occupied leaf no longer re-casts its rays in ACA.
3. **Default `tol` 1e-5 → 1e-6.** With the surface tree, capsule–wall blocks compress
   where they used to give up and fall back to exact dense entries. At 1e-5 that moved
   P4/P0 by 1.7e-7 at 28.6k, against 1.6e-8 for the old tree. At 1e-6 every moment is
   within ~1e-8 of the references and the time is still well below the old tree's at
   1e-5.

### After (warm assembly, 4 threads; old = 9e87912 at tol 1e-5, new = tol 1e-6)

| N | 10,484 | 17,632 | 28,628 | 49,800 |
|---|---|---|---|---|
| view factors, old | 7.6–8.4 s | 16.6–17.5 s | 33–42 s | 90–98 s |
| view factors, new | 5.3–6.3 s | 13.2 s | 22–25 s | 56.5 s |
| stored, old → new | 122 → 99 MB | 240 → 212 MB | 432 → 439 MiB | 1004 → 1014 MiB |
| max RSS (whole bench), new | 0.80 GB | 0.94 GB | 1.25 GB | 2.15 GB |
| dense entries, new | 19M | 38M | 82M | 196M |
| F x, new | 8 ms | 18 ms | 42 ms | 114 ms |

Assembly now grows as N^1.4 from 17.6k to 49.8k (N^1.6 before). Storage is unchanged
(N^1.5): what is still dense is the near field and the shadow-edge leaves. Profile of
the new assembly at 28.6k (inclusive samples): ACA 35 %, dense leaves 31 %, rows cut by
the facing test (scattered) 23 %, shafts 2 %, ball masks 1.6 %, ray masks ~0. The rays
are gone from the profile; what remains is evaluating the kernel itself.

### Accuracy (balanced, T-dependent solve at 2 ns; a₀ relative, P-ratios absolute)

| N, reference | a₀ | P2/P0 | P4/P0 | H per element, max (median) | moment response |
|---|---|---|---|---|---|
| 10,484, dense Float32 | 1.4e-9 | 2.6e-9 | 2.7e-9 | 1.2e-7 (4.7e-9) | 7.8e-7 |
| (old tree, tol 1e-5) | 3.1e-8 | 4.2e-8 | 2.6e-8 | 8.6e-5 (3.6e-8) | 6.6e-5 |
| 17,632, dense Float32 | 1.7e-10 | 6.5e-9 | 5.9e-9 | 1.3e-7 (6.7e-9) | 1.1e-6 |
| (old tree, tol 1e-5) | 1.4e-9 | 2.7e-10 | 2.8e-8 | 8.6e-6 (4.7e-8) | 4.6e-6 |
| 28,628, new tree at tol 1e-7 | 2.0e-9 | 9.7e-10 | 9.6e-10 | 7.8e-7 (6.6e-9) | 2.4e-6 |

The old tree's tol-1e-7 matrix at 28.6k is off from the new tree's by 1.5e-5 per element
whatever the tolerance; `rowchk.jl` puts that on the old tree (the dropped capsule
entries of item 1, checked against `_vf_exchange` pair by pair). The small-target test
entries (`test/hmatrix.jl`, leaves of 16) still match `view_factors` to < 1e-5 of the
largest, and the 3.9k hohlraum at leaves of 32, 2.3e-5 before, now matches to 1.5e-6.

### Tried and dropped

- **Skipping ACA on blocks the facing test cuts in most rows** (min wholly-facing
  fraction < 0.5; the predictor above). With the surface tree it skipped 22k attempts
  and left 279 give-ups. Time was unchanged (22.0–23.2 s against 22.1–23.0 s at
  28.6k), and it stored 6 % more (464 MiB) and 10 % more dense entries: the skipped
  blocks that would have compressed cost as much as the failed attempts saved.
- **Kernel × mask for shadow-edge leaves** (ACA on the unoccluded kernel, expanded on
  the mask, sampled check): 21.4 s against 22.2–22.5 s. With the balls the mask is
  cheap, the leaf's kernel evaluations are only ~15 % of the dense entries, and the
  product of a Hadamard-stored block costs rank × a dense one. Not worth the code.

### Still open

- The near field and the rows cut by the facing test (a quarter of the time):
  quadrature-level reuse between neighbouring pairs, or a coarser near-field rule where
  the exact contour kernel is not needed.
- `:sampled` visibility and non-convex occluders (LEH rims, CAD) still take a ray per
  pair of an occupied leaf. A rim certificate (a disc or annulus the segment must
  cross) would extend the ball test. (Done for centroid visibility: "H-matrix visibility
  certificates".)
- BriareusGUI `radiation_size` now budgets an HMatrix at 250 N^1.5 bytes (2.5 × the
  measured storage) above `HMATRIX_MODEL[]` and 4N² between 10k and 15k elements.

## H-matrix near field: zero blocks and cut rows (2026-10-26) *derived*

ROADMAP 14, the first "still open" item above. Measured here (same hohlraum workload and
`bench3.jl`-style script, 4 threads, `tol` 1e-6, leaf 48; old = 1f20b2d). Two changes:

1. **All-zero blocks dropped.** A dense leaf pair or scattered block whose Float32
   entries are all exactly zero is not stored (a capsule's near field against itself
   across the sphere, rows cut at a terminator against a box's far side). Exact entry
   for entry: the product skipped `muladd(0, z, y)`. Only the split of blocks into
   product groups (by cost) moves, so threaded products sum in another order (round-off).
   At 17.6k: 6.6k of 25.7k dense leaves were zero, and 4.3M zero entries in scattered
   blocks.
2. **Facing-cut rows refined over the column tree** (`lowrank_index` blocks). The rows
   of a low-rank block that face only part of the other cluster (stored dense before)
   are re-tested against the column cluster's children, where most face wholly
   (indexed low-rank block, own ACA and sampled check) or not at all; dense only against
   a leaf or when ≤ 4 rows remain. These are new approximations at the same `tol`.

| N | 10,484 | 17,632 |
|---|---|---|
| stored, old → new | 99.4 → 81.8 MiB (−18 %) | 212.0 → 175.1 MiB (−17 %) |
| dense entries, old → new | 18.9M → 14.1M | 37.8M → 27.6M |
| indexed low-rank blocks | 52 | 215 |
| view factors (warm, 2 runs), old / new | 5.9–6.0 / 6.3–6.5 s | 11.5–11.6 / 11.7–11.8 s |
| F x, old → new | 6.6 → 5.6–5.9 ms | 16.6–17.1 → 13.5–13.7 ms |
| F X (16 columns), old → new | 93 → 80 ms | 256 → 212–217 ms |

New against old, balanced T-dependent solve at 2 ns (a₀ relative, P-ratios absolute):
a₀ −7.5e-12 / 1.4e-10, P2/P0 2.5e-11 / 3.8e-11, P4/P0 −1.5e-12 / 3.7e-11, H per element
max 1.4e-8 / 1.4e-8, moment response 1.4e-8 / 1.6e-8 (10.5k / 17.6k); against dense
Float32 at 10.5k unchanged (a₀ 1.4e-9, P2/P0 2.7e-9, H max 1.2e-7). The previous
session measured −15 to −27 % storage over 10.5k–28.6k. Assembly time is unchanged to
+8 %: the cut rows' extra ACA costs about what the dense entries did.

### Vectorized sector contour (H-matrix only)

A profile of the 10.5k assembly put ~45 % of its work in `_contour_vf` over polar-sector
loops (34 edges, one scalar `atan` each: wall and capsule points onto the end caps). The
H-matrix assembly now gives sector receivers their loop as x, y, z arrays
(`_SoAReceiver`) and takes the contour in one `@simd` loop with a branch-free
`atan(y ≥ 0, x)` (`_atan_nb`: fdlibm's polynomial after reduction to |u| ≤ tan(π/8);
≤ 2 ulp from `Base.atan` over 10⁶ random arguments, exact at 0, π/4, π/2, π), 2.4× faster
per contour (714 → 306 ns). Not bit for bit: the terms are summed in another order.
`view_factors` (dense) keeps `_Receiver` and is bit for bit (`==` on two hohlraums,
2.0k and 4.0k elements). Measured against 5173a3d (4 threads, tol 1e-6, leaf 48):

| N | 10,484 | 17,632 |
|---|---|---|
| H-matrix build, warm (3 runs), old → new | 6.07–6.20 → 4.38–4.45 s (−28 %) | 11.53–11.59 → 8.78–8.89 s (−24 %) |
| entries, sampled rows (400 / 300 rows, all columns) | 0.77M of 2.0M differ; max abs ΔF 3.9e-17, max rel 2.5e-12 | 0.85M of 2.6M; 4.2e-17, 3.6e-12 |
| blocks, ranks, storage | identical | identical |
| a₀ (rel), P2/P0, P4/P0 (abs) | 2.1e-11, 9.4e-12, 3.9e-13 | −5.6e-16, 3.1e-14, 9.8e-15 |
| H per element max (median); moment response | 7.4e-9 (2.0e-11); 1.2e-8 | 2.6e-12 (3.6e-15); 2.5e-12 |

(Balanced T-dependent solve at 2 ns, `bench3.jl`.) The largest relative entry
differences are small factors whose contour terms cancel (7e-7 with 1.1e-17 absolute).

**Contiguous quadrature layout: tried, dropped.** The per-element point rules (`qnear`,
`qmid`, `qfar`: three small arrays per element) were packed into three arrays with
offsets, read as views. Neither faster nor bit for bit: H-matrix build +2 % (4.47–4.55 s
at 10.5k, 9.03–9.05 s at 17.6k, on top of the contour change); dense Float32 11.79 against
11.82 s and 30.2 against 30.8 s (one run each). Dense F changed in 15 % of entries (max
rel 3.9e-13), from `_point_kernel` and `_self_kernel` with the weights as views:
StaticArrays' `dot` carries LLVM `reassoc contract` flags, so where FMAs are fused
depends on the surrounding code, and a change of data layout in a kernel is not bit for
bit in general.

## H-matrix near-field quadrature: where the time goes (2026-10-27) *derived*

ROADMAP 14, "the near-field quadrature itself". Same workload (`bench3.jl`, 4 threads
pinned to one CCX, `tol` 1e-6, leaf 48); old = d238f5d, timed interleaved in one sitting
(other jobs on the machine's other CCX, so absolute times run ~15 % above the previous
sections').

**Profile (10.5k, before).** Attribution by category of `_hassemble!` × class of the pair
kernel (tier × pair type; instrumented, ~1.7× the real time, so shares only):
near-field dense leaves 27 % of block assembly, ACA 37 %, rows cut by the facing test
(refined, split, sparse) 13 %, dense after a give-up 8 %, shadow-edge leaves 6 %. Of the
kernel time, a curved point (wall, capsule) onto an end-cap sector's 34-edge contour
("pc") is the largest class everywhere: far tier 6.3 s (3.7 s of it inside ACA), mid
2.1 s, near (adaptive rule) 1.3 s, against far point-to-point 2.3 s and the visibility
rays of occupied leaves 3.6 s. In the near-field leaves themselves: mid pc 1.6 s (0.5M
pairs, 3.2 µs), near pc 1.3 s (6.9k pairs, 194 µs each), mid point-to-point 0.8 s, far pc
0.5 s, rays 0.4 s, near point-to-point 0.3 s. Single pair costs (one thread, 20k random
pairs per class): far pc 1.24 µs, mid pc 2.84 µs, far sector–sector 2.37 µs, near pc
200 µs; of a near pc pair, ~half was `_distance_to` (32 point–triangle distances at
every cell the adaptive rule visits), a quarter the contours.

**Coarser rule: not viable.** The point rules' error is O((h/r)²): at the far tier's
12 element sizes a point-to-point or one-point-per-cell sector rule is off by ~1e-4 to
1e-3 relative, 100–1000× the block tolerance, and changes the discretization rather than
approximating the same matrix. An exact-enough quadrature of the same polygon (Gauss on a
sector to < 1e-9 needs ~16 points plus the 32 chord segments) costs as much as the
contour. Quadrature reuse between neighbouring pairs has nothing to share in this
workload: the dominant pairs are distinct points onto distinct receivers.

**What changed** (all on the receiver side; the cost is in evaluating the same kernel):

1. `_closer_than` on a sector receiver (`_Receiver` with a loop) takes the distance from
   the loop in its plane's coordinates (height inside the polygon, else the nearest of
   34 edges, one pass) and asks the 32 triangles only when the verdict is within 1e-9
   (relative) of the threshold; `_adapt_tri` and `_adapt_sector` ask it through
   `_within`. The refinement, and so every factor, is unchanged **bit for bit**: `==`
   on `view_factors` (dense) of `hohlraum(h = 0.6, capsule_h = 0.3)` and both bored
   cavities, and 0 wrong verdicts of 2.9M (points 1e-4 to 1 mm from sectors, in and off
   their plane, thresholds at D(1 ± 1e-6), D/2, 2D and D). Near pc pair 167 → 94 µs.
2. The vectorized contour counts the loop's vertices on or behind the point's tangent
   plane in the same `@simd` loop (an `Int` sum; a floating-point `min`/`ifelse`
   reduction did not vectorize and made the loop 2× slower) instead of two scalar
   `all` passes over the loop.
3. A loop the tangent plane cuts (1 % of points, 4 µs each through 32 clipped
   triangles) is clipped as one polygon (`_clipped_loop_vf`, streamed
   Sutherland–Hodgman; the coincident edges along the cutting line cancel), 2.5 µs.
   Against the triangles over 44.8k random cutting points: max abs 1.7e-14, within the
   test's 1e-12 rel + 1e-15 abs.
4. `_atan_nb` forms (t − 1)/(t + 1) of t = lo/hi as (lo − hi)/(lo + hi), one division
   instead of two; still ≤ 2 ulp from `Base.atan` (10⁶ arguments over 24 decades).

Pair costs after (same pairs): far pc 1.04 µs (−16 %), mid pc 2.43 µs (−14 %), far
sector–sector 1.83 µs (−23 %), near pc 114 µs (−43 %).

| | 10,484 | 17,632 |
|---|---|---|
| H-matrix build, warm, old → new (interleaved) | 4.30–4.43 → 3.89–4.06 s (−9 %) | 8.84–9.08 → 8.15–8.35 s (−8 %) |
| entries, 300 sampled rows (all columns) | 63 of 3.1M differ; max abs 8.9e-17 | 88 of 5.3M; max abs 1.9e-17 |
| blocks, ranks, storage | identical | identical |
| a₀ (rel), P2/P0, P4/P0 (abs) vs old | 0, −8.3e-17, 4.0e-17 | 0, 0, 5.2e-17 |
| H per element max (median); moment response vs old | 1.3e-15 (0); 2.5e-15 | 4.9e-15 (0); 3.9e-15 |

Against dense Float32 at 10.5k unchanged: a₀ −1.37e-9, P2/P0 2.67e-9, P4/P0 2.67e-9, H
max 1.19e-7 (median 4.7e-9), moment response 7.8e-7. (The largest relative entry
differences, up to 0.7, are factors of ~1e-17 whose contour terms cancel.) Dense
`view_factors` gains from item 1 alone, bit for bit: `bored_cavity(hw = 0.8, hc = 0.25)`
1.72–1.75 → 1.35–1.38 s, `bored_cavity()` 2.35–2.36 → 1.98–1.99 s (4 threads).

**What is left.** After this the H assembly is mostly the contour itself (≈ 9 ns per
edge: a square root and two divisions per edge, vectorized), point-to-point pairs whose
cost is the scattered loads of three small arrays per element (the contiguous layout
above did not help), and the rays of occupied leaves. None of the remaining near-field
ideas measured here is worth 10 %: a further step needs fewer kernel evaluations (rim
certificates for the rays; fewer ACA give-ups), not cheaper ones.

## Near flat pairs by the double contour (2026-10-27) *derived*

ROADMAP 6, "Adjacent flat CAD pairs". Scripts in the session scratchpad `adj/`
(`workloads.jl`: mock-ups only — a closed box of primitive quads, the flat and fitted
remeshed CAD hohlraum of the "Radiation in an imported CAD hohlraum" tutorial at
0.5/0.15 mm, the 2,000-strip sliver tube remeshed at 0.6; `profile.jl`, `compare.jl`,
`t2.jl`–`t5.jl`, `calib.jl`, `scope.jl`). Ryzen 9 4900HS, pinned to one CCX, other jobs on
the other.

**Share first (one thread, before; pair kernels classified by tier and by shared
vertices).** Pairs sharing an edge took 71 % of a box's view-factor time (4 × 4 quads
per face: 72 pairs at 3.5 ms each), 58 % of the flat CAD hohlraum's (1,959 pairs at
8 ms; all of them involve an element folded across the wall's strips) and 87 % of the
sliver tube's (1,482 pairs at 74 ms: 130-facet folded elements). The fitted CAD model
(walls on their cylinders) has almost none. Not a negligible share wherever flat
elements meet.

**Where the old rule's error was.** The cell rule refines toward the receiver to
maxdepth 8, but its one-point cells away from the shared edge set the error: one
element per plate, against Hottel's perpendicular rectangles, it is −5.6e-4 (1 × 1),
−4.4e-4 (1 × 0.6), −3.4e-4 (20:1 slivers), −5.0e-3 (shared edge 1/20 of the width);
maxdepth 13 changes none of these by more than 2e-5. On folded CAD elements the cell
rule emits from the element's mean plane: neighbors on the wall, folded by a few
degrees, exchanged 10–70 % off a facet-level converged quadrature (median 16 %, max 24 %
over 12 sampled edge-sharing pairs; 7 % / 10 % for vertex-sharing, 3 % / 53 % for
near pairs apart). Their factors are small (G ~ 4e-4 mm² against 0.25 mm² elements),
so rows moved by ≤ 2.2e-3.

**Adopted: the double contour integral** (Sparrow 1963; Mitalas & Stephenson 1966):
G = A_i F_ij = (1/2π) Σ_a Σ_b (â·b̂) ∫_a∫_b ln r ds dt over the directed outline edges,
exact for any two surfaces that see each other wholly from in front, folded or not
(`_contour_pair_vf`). Used for every near pair of flat elements (not curved, not
sectors) whose vertices all lie in front of the other's facet planes (to 1e-9 of the
size) and whose own facets lie in front of one another (no self-hiding fold); otherwise
the cell rule as before. Inner integral closed form; outer adaptive Gauss–Kronrod 7–15
with a round-off floor; edges meeting at a point have their s ln s subtracted
analytically, collinear edges are closed form (x² ln|x|/2 − 3x²/4), an edge is split
where the other ends on it, and edges far apart for their lengths take product
Gauss–Legendre (n = 6…2 for ρ = 2|Δmid|/(La + Lb) ≥ 16, 24, 40, 120, 2000: worst error
over 3,000 random segment pairs ≤ 1e-13 La Lb against 40-point product Gauss on 4 × 4
panels, the reference's own round-off). Coordinates are scaled to the pair's size; the
constant part of ln r cancels over closed outlines. The reverse direction is the same
integral, so reciprocity is exact without averaging.

Accuracy (relative to the independent reference):

| case | old (cell rule) | new |
|---|---|---|
| perpendicular rectangles, Hottel: 1 × 1 × 1, 1 × 1 × 0.6, 1 × 0.05 × 0.02, 0.05 × 1 × 1 | −5.6e-4 … −5.0e-3 | ≤ 2.2e-15 |
| unit squares at 45°, 135°, 170°, 1 × 0.5 × 2 at 60°: cell rule against the contour | +4.5e-4, −9.4e-4, −1.0e-3, −1.8e-3 (maxdepth 13 the same) | (the prism below checks non-right angles) |
| CAD flat: edge / vertex / apart near pairs vs facet-level quadrature (mindepth 6, maxdepth 12) | median 1.6e-1 / 6.8e-2 / 2.9e-2 | 1.5e-5 / 1.4e-5 / 7.4e-6 (the reference's error) |
| box: same classes | 5.5e-4 / 1.7e-3 / 1.5e-3 | 6e-6 / 1.1e-5 / 5.7e-6 |
| closed cube, triangular prism (one element per face): row sums | — | 1 to < 1e-12 |
| cube with bottom + front one folded element: its row sum | — | 1 − F(Hottel) to 1e-12 |
| closed box, 4 × 4 per face: closure error max / mean | 5.3e-3 / 2.9e-3 | 5.3e-3 / 2.8e-3 |
| perpendicular plates n × n, Hottel (H = 1; 0.6), n = 2 / 4 / 8 | +2.3e-4 / +2.4e-3 / +2.2e-3; −1.7e-4 / +1.2e-3 / +2.2e-3 | 2e-16 / +2.2e-3 / +2.1e-3; 4e-16 / +1.0e-3 / +2.1e-3 |

What is left in the plate sums is the mid tier (3–12 element sizes, one centroid point
onto the exact polygon, O((h/r)²)): taking the contour there too makes every case above
exact to 3e-15 and the closed box's rows one to 7.7e-15, at 3.4× the box's time (0.07 →
0.24 s at 8 per face), 2.3× the CAD hohlraum's (2.0 → 4.9 s) and 16× the sliver tube's
(132-edge outlines; 6.1 → 98 s). Not adopted (an open decision).

Timings, dense `view_factors`, 4 threads (old = the same build with the contour
disabled, interleaved; identical otherwise):

| model | N | old | new | entries changed | max ΔF / max Δ row sum |
|---|---|---|---|---|---|
| `hohlraum(h = 0.6, capsule_h = 0.3, sombrero = true)` | 1,888 | 0.45 s | 0.45 s | 0 (`==`) | — |
| `bored_cavity(hw = 0.8, hc = 0.25)` | 1,068 | 1.32–1.40 s | 1.32 s | 0 (`==`) | — |
| box, 8 × 8 quads per face | 832 | 0.24 s | 0.07 s | 7,392 | 1.4e-4 / 4.3e-4 |
| CAD hohlraum, flat 0.5/0.15 | 1,693 | 7.3–7.4 s | 2.1–2.3 s | 25,402 | 3.5e-3 / 2.2e-3 |
| CAD hohlraum, fitted 0.5/0.15 | 1,693 | 2.60 s | 2.59 s | 248 | 1.7e-3 / 1.6e-3 |
| sliver tube 0.6 (open) | 518 | 34.1 s | 7.2 s | 14,190 | 4.9e-4 / 6.4e-4 |

Per pair (one thread): an edge-sharing pair of quads 3.3–7.6 ms → 2–100 µs (slivers the
most), folded CAD wall pairs 8 ms → 0.4 ms, tube pairs 74 ms → 3.8 ms (132 × 132 outline
edges). Non-touching near pairs: the box's faster (4 × 4 per face: 0.013 → 0.006 s with them,
against edge-sharing pairs only), the CAD hohlraum's 25 % faster, the tube's ~15 % slower
(its total still 4.7× faster). Primitive hohlraums and bored cavities have
no flat non-curved elements and are unchanged bit for bit; the H-matrix assembly goes
through the same pair kernel. The flat CAD tutorial's radiation metrics move by little
(drive of the tutorial, balanced F): capsule flux +4.2e-5, a₂/a₀ 10.587 → 10.584 %,
a₄/a₀ −0.788 → −0.787 %, LEH and gap losses unchanged to 1e-5. The sliver tube's mean row
sum moves 0.72210 → 0.72197 (the open cylinder's exact 0.72301; the rest is the mid and
far tiers). Confidence: the method and the closed forms **H**; the ρ thresholds and the
1e-9 front tolerance engineering choices, **M**.

## H-matrix visibility certificates (2026-10-27) *derived*

ROADMAP 14, "rim certificates for non-convex occluders". Workloads: `bench3.jl`'s hohlraum
(10,484 and 17,632 elements) and a target where rims occlude, the bored cavity of
`test/partial_visibility.jl` (LEH bores, cap washers, a 20 µm waist gap with its washers;
latitude–longitude capsule) at `hw = 0.18, hc = 0.07` (11,360 elements). 4 threads pinned
to one CCX, `tol` 1e-6, leaf 48; old = abb3129 from a separate worktree, interleaved. Other
jobs shared the machine: process-to-process spread ~15 %, so ranges below are of 3 builds
per process, 2 processes each.

**Measured first.** Visibility rays timed by instrumenting `_vf_visible` (CPU-s inside
it over a warm build):

| | hohlraum 10.5k | bored cavity 11.4k |
|---|---|---|
| rays (old) | 2.9M, ~3.8 CPU-s (~23 % of assembly CPU) | 35.4M, 19–23 CPU-s (43 %) |
| of which occupied leaf pairs | ~0.4 % of assembly | ~32 % far leaves + 9 % near leaves |
| of which pairs only the facing slack admits, in empty shafts | ~23 % (wall–capsule 1.4M, wall–wall 1.0M, caps–capsule) | ~10 % |

So the rim question was a negative on the standard hohlraum (its occupied leaves are
~0.4 %) and a large positive where rims occlude; and on the standard hohlraum the rays
were nearly all of a different kind: pairs that face each other only within the slack of
the centroid facing test on curved elements, whose rays (lifted farther off the facet) are
cast even in an empty shaft.

**What changed** (`src/hmatrix.jl`; all exact, so `HMatrix` is **bit for bit**):

1. Certificates per component, built once (`_hrim_data`): a *planar* part (all facets
   within 1e-10 relative of one plane: washers, LEH caps, discs, foils) keeps its facets
   as 2-D edge lines binned on a polar grid; a *cylinder or cone* part keeps its axial and
   radial bounds about the surface axis; a convex part with balls keeps them. The shaft
   search of a leaf pair (`_shaft_occluder(…, full = true)`) now collects up to 8 such
   parts instead of giving up at the second one or at a non-convex one.
2. Each pair of a certified leaf is decided as `_vf_visible`'s ray (`any_hit` to
   1 − 1e-9, barycentric slack 1e-10) would: hidden if one certificate finds the segment
   through its facets, visible if all find it clear, else the ray. Planar: one side of
   the plane (margin 1e-9 relative plus twice the flatness), or the crossing outside every
   facet in its bin, or inside one with both ends off the plane — margins grown by
   1/sin of the crossing angle, and undecided when they exceed half a bin. Hull: the
   segment's part within the axial range (widened) inside the radial lower bound or
   outside the upper one (visible only). Ball: as before, with the per-pair condition that
   both ends lie outside the circumscribed ball (so a pair whose cluster box reaches into
   it no longer forces every pair to the ray), plus a separating plane: an end element
   lying on the convex part (its facets in one plane through its centroid) clears it when
   both ends are in front by 1e-8 × the target scale.
3. Pairs only the slack admits are judged on their own ray ends (`_hslack_ends`), when
   both ends stay inside their cluster boxes (lift ≤ the padding), every plane exclusion
   of the shaft held with the clusters' largest extra lift (`liftx`, recorded by the
   shaft search), and each convex part a cluster lies on (excluded from the search) is
   cleared by an end element's plane. This applies to empty shafts (the `_HOcc` of the
   shaft replaces `true` as the block's entry mode) and to certified masks; centroid
   visibility only.
4. Blocks above leaf level are unchanged (only the old single-ball masks), so the block
   structure, ranks and every stored entry stay as they were. (Superseded by the second
   pass below: certificates above leaf level, with a quadrature-tier check.)

| | hohlraum 10,484 | hohlraum 17,632 | bored cavity 11,360 |
|---|---|---|---|
| H-matrix build, warm, old | 4.46–4.69 s | 9.42–9.68 s | 9.79–11.40 s |
| new | 3.44–4.54 s (typ. 4.0) | 7.87–8.37 s | 7.74–7.92 s |
| change | ≈ −10 % | −16 % | −20 to −30 % |
| rays old → new | 2.9M → 0.98M | | 35.4M → 1.34M (0.9 CPU-s) |
| stored blocks, all entries | identical (`==` on every dense, low-rank, scattered block) | identical | identical |

The remaining hohlraum rays are capsule elements at their terminator facing slightly away
(c ≤ 0 within the slack), whose own plane does not separate; all visible. The wall-clock
gain is smaller than the ray CPU saved (bored: ~20 CPU-s of rays → ~10 CPU-s wall × 4):
the certificates cost ~140 ns per pair, and the shaft search no longer stops early on a
leaf.

**Exactness checks** (`test/hmatrix.jl`, "visibility certificates"): 57k random segments
across the parts of the bored cavity and a hohlraum, and segments crossing a washer's plane
within 1e-3 to 1e-12 (relative) of its inner and outer rims: every decided one (49k, 5,023
of them hidden) agrees with `any_hit` against that part's facets alone; masks of every
certified leaf pair equal the ray masks (6,193 leaf pairs over three targets, a
latitude–longitude capsule among them); the 833k pairs of empty-shaft leaves get
`_vf_pair`'s factor (`===`). At the benchmark sizes the certified masks equal ray masks on
all 9,889 certified leaf pairs of the bored cavity (12.4M pairs).

**Certificates above leaf level (adopted, second pass).** The certificates also decide
blocks above leaf level (`_shaft_occluder(…, full = true)` for every far block under
centroid visibility; a mask for blocks up to `mask_size`), so a block of uniform
visibility gets a cross approximation instead of being refined to leaves. A first try
also counted latitude–longitude capsule quads (two coplanar facets) as "on the convex
part" in the shaft search; that rule is **not** adopted (below). The first try lost row
accuracy against dense Float32 (2.4e-7 → 2.3e-6); the loss was diagnosed as follows.

- *One entry.* Capsule element 7631 is a polar-fan triangle of size 0.265 (3–4× its
  neighbours). Wall element 1380 lies at r = 3.1825 from it, just past the far boundary
  12 × 0.2649 = 3.1783. The dense matrix evaluates that pair with one point per triangle
  (far tier) and the other 42 pairs of the same column with the level-1 rule (mid tier).
  The rules differ by 0.53 % there: F = 4.2547e-4 (far) and 4.2319e-4 (mid). The
  converged value (level 4) is 4.2246e-4, so the dense entry is 0.71 % off and the mid
  value 0.17 %. In the low-rank block (43 wall × 107 capsule), the cross approximation's
  pivots and the check's 75 random samples all miss the isolated jump, and the block
  stores the mid-tier value. Its error against dense is 5.1e-5 of the block's norm; the
  other 4,728 low-rank blocks are all ≤ 7.6e-6. The row sum of 7631 is off by 2.28e-6.
  The dense row sum of that sliver is itself 1.00188 against 0.999998 converged, so the
  H-matrix's departure is ~1000× below the dense quadrature's own error, and toward the
  converged value.
- *Not the certificates.* The masks of every block that takes one are exact: 11,134
  blocks (8,763 leaves), 33.5M pairs (3.4M hidden), all equal to the ray masks. The 20.5M
  pairs of the shafts that the capsule rule empties get `_vf_pair`'s factor (`===`).
- *Not the capsule rule alone.* Above-leaf certificates alone reach the same block: its
  parent's shaft holds capsule facets, which the balls certify pair by pair. They give
  the same 2.28e-6. With both rules the blocks are identical to above-leaf alone, so the
  capsule rule adds nothing there. By itself it saves ~3 % (7.76 against 7.97 s), and it
  is dropped.
- *Fix.* The low-rank check (`_hcheck`) now adds the exact error of the entries where a
  row's or a column's quadrature tier switches (`_vf_tier`, `_htier_entries`). These are
  the pairs off their line's most common tier, taking the lines with the fewest such
  pairs first (an isolated jump is the one the pivots miss), at most m + n of them. There
  are none, and no pass over the pairs, when the centroid boxes are farther apart than
  `mid` × the largest element. The block above is rejected and refined. A block with no
  such entries is checked as before, bit for bit.

Bored cavity (11,360; 4 threads pinned, interleaved in one process, 2–3 builds each);
balanced solve as `bench3.jl` (NIF quads aimed at z = ±4 mm deposited by
`DepositionOperator`, `closed = true`, 2 ns), against the dense Float32 matrix:

| | committed | capsule rule | above-leaf | above-leaf + tier check (adopted) |
|---|---|---|---|---|
| build, warm | 7.0–8.2 s (typ. 7.9) | 7.8 s | 6.7 s | 5.7–6.6 s (typ. 6.4; ≈ −18 %) |
| stored | 117.8 MB | 115.7 MB | 104.9 MB | 104.9 MB (−11 %) |
| Frobenius error, 300 columns | 6.49e-8 | 7.07e-8 | 7.76e-8 | 7.75e-8 |
| row sums, max \|Δ\| | 2.42e-7 | 2.28e-6 | 2.28e-6 | 2.63e-7 |
| H per element, max (median) | 1.01e-7 (3.5e-9) | 5.5e-7 | 5.5e-7 (5.7e-9) | 1.10e-7 (5.7e-9) |
| capsule a0; P2/P0, P4/P0 | −1.7e-11; 1.1e-9, 6.6e-11 | | | −7.2e-10; 3.2e-9, −1.5e-9 |

Hohlraum: above-leaf alone leaves the blocks identical to the committed build (10,484).
The tier check rejects one more block at each size (10 → 11 at 10,484; 16 → 17 at
17,632), so the hohlraum is **not bit for bit**. At 10,484 the Frobenius error (2.679e-7)
and the row sums (3.214e-7) are unchanged to four digits, storage is 85.7 → 85.8 MB, and
the build stays at 4.0–4.2 s. At 17,632 storage stays 183.6 MB and the build 7.8–8.0 s.
The dense matrix is bit for bit (`_vf_tier` is `_vf_kernel`'s tier test, factored out).

Tests (`test/hmatrix.jl`): the mask-against-ray and empty-shaft checks of "visibility
certificates" now cover every block up to `mask_size`, not only leaf pairs. "Quadrature-
tier entries" uses a latitude–longitude hohlraum (3,254 far blocks): every entry returned
is off its row's or its column's majority tier; every pair alone off its line's majority
is returned when the cap is not reached; and a block with no entries has one tier.
Confidence: the diagnosis and the measurements **H**; the m + n cap and the
fewest-first order are engineering choices, **M**.

**`:sampled` visibility (third pass, adopted).** Same bored cavity with
`visibility = :sampled` (11,360 elements, 4 threads pinned, warm).

*Measured first* (old build 106 s with `_vf_visible` instrumented): 48.7M calls,
404 CPU-s, 93 % of the assembly CPU. By caller: 26.7M at occupied leaf pairs (227 CPU-s),
21.8M entries of near leaves and blocks without a certificate (176 CPU-s), 0.24M slack
pairs of empty shafts (1.3 CPU-s). 1.5M pairs disagreed with their probes and took all
256 sample rays. Under `:sampled` the old shaft search took no certificate but the
whole-ball one, so every occupied leaf took rays.

*What changed* (`src/hmatrix.jl`, `src/viewfactors.jl`):

1. Bundle certificates. All rays of a pair's `_visible_fraction` (the centroid ray, the
   probes, the 16 × 16 sample segments) run between points within `rvis[i]` of oi and
   `rvis[j]` of oj. Here `rvis` is the largest sample distance from the centroid. Each
   point of such a segment lies within ρ = max(rᵢ, rⱼ) of the centroid segment at the
   same parameter. So each certificate is applied to the centroid segment with margins
   grown by its tube:
   - ball: the ends outside the outer ball by rᵢ, rⱼ; clear past rout + ρ; hidden
     inside rin − ρ;
   - planar: the end heights shifted by rᵢ, rⱼ. The crossings lie within
     ρ (1 + |w|/|dh|) of the centroid crossing, since a point at height ≤ ρ on the
     axis lies that far along it from the crossing. The 1/sin margin uses
     (|w| + rᵢ + rⱼ)/(|dh| − rᵢ − rⱼ);
   - hull: axial range widened by ρ, rlo − ρ, rhi + ρ;
   - an end element's plane clearing a convex part: its own samples within `hdev` (≈ 0)
     of its plane, the other end's within rⱼ.

   Zero radii reproduce the old arithmetic bit for bit. A bundle is 0 if any part hides
   it whole, and 1 if all clear it whole. Otherwise `_visible_fraction` runs with each
   segment decided by the parts that did not clear the bundle (`_HSeen`), or by the ray.
2. Exclusions valid for the samples. In `:sampled` mode the shaft's plane tests take the
   box of every element's samples (and centroid) lifted by ε (`opts.qbox`) instead of the
   centroid box. A cluster's in-plane offset (`cplane`) also comes from those points.
   Per pair: the ends stay inside the boxes (rvis + lift ≤ the padding); a slack pair
   needs `occ.slack`; each excluded convex part must be cleared by an end element's plane
   for the bundle, else per segment, else the ray. Folded elements take `_vf_visible`.
3. Assembly: `:sampled` blocks take the full shaft search. Blocks up to `mask_size`
   with only certified parts get a fraction mask (`_HFrac`): wholly visible blocks get
   the cross approximation at any level, others are refined; near leaves are masked too.
4. Zero-weight sample pairs (facing away at an end) add nothing to either sum of the
   fraction, so their segments are tested only when no pair faces (the area-weighted
   fallback). The output is unchanged, and it removed 41M segments into the capsule.
5. One rounding for every caller. StaticArrays' `dot` is an `@simd` reduction, so it
   rounds as its caller is compiled: `ci` differed in the last bit between `_vf_visible`
   and the bundle code (a slack lift, so a fraction off by 6e-17). The ray ends
   (`_vf_ends`) and the sample weights (`_vis_weight`) are now `@noinline` with concrete
   arguments, one compiled instance for all.

| bored cavity 11,360, `:sampled` | old (7ec13bb) | new |
|---|---|---|
| H-matrix build, warm (interleaved) | 93.3, 93.6 s | 29.9, 29.7 s (3.1×) |
| stored | 121.0 MB | 108.3 MB |
| low-rank / masked blocks | 15,437 / 0 | 4,971 / 4,086 |
| `_vf_visible` calls | 48.7M | 0.61M |
| pairs decided by bundle (hidden / visible / per segment) | — | 1.41M / 23.5M / 11.3M |
| segments per segment (certified clear / hidden / ray) | — | 187M / 116M / 67M |
| Frobenius error, 300 dense columns | 5.98e-8 | 7.82e-8 |
| row sums, new − old | | ≤ 2.5e-7 |

The error is at the block tolerance (1e-6) as before; the comparable centroid change was
6.5e-8 → 7.8e-8. Most of the 67M rays left are segments through the 20 µm waist gap that
end on the wall. The wall's axisymmetric bounds cannot decide them, since the end lies
between the inscribed and circumscribed radii. A certificate for the inside of a convex
faceted wall (prism) would be needed: ~40M rays, roughly a quarter of the build.

*Exactness checks.* On the benchmark all 13,316 certified blocks (leaf pairs and blocks
up to `mask_size`) were checked. Their 33.2M facing pairs (2.92M hidden, 1.02M partial)
have `_hfrac` `===` `_vf_visible`'s fraction. The 7,644 empty-shaft blocks have
`_vf_pair_free` `===` `_vf_pair` per pair. No block of this target needs rays any more.
`test/hmatrix.jl` "sampled visibility bundles" (1,068 elements) checks:
- bundle verdicts against the rays of random segments inside the balls (some ends on
  their surfaces), against that part's facets alone;
- every certified block's fractions and every empty shaft's factors (`===`).

Centroid-visibility H-matrices are identical block by block (bored cavity 11,360,
hohlraum 10,484; 5.9 s and 3.5 s, unchanged). Dense centroid factors are bit for bit.
Dense `:sampled` factors are not: 3.8 % of entries move, by ≤ 4e-14 relative (item 5).
They are 5 % faster from item 4 (3,156 elements: 10.05 → 9.50 s). The dense path has no
shafts, so the bundle certificates do not apply to it. Confidence **H** for exactness and
timings; the gap-certificate estimate is **M**.

## Near sector pairs by the double contour with arcs (2026-10-27) *derived*

ROADMAP 14, "the bored cavity's gap pairs". Scripts in the session scratchpad `gap/`
(`share.jl`, `attr.jl`, `near.jl`, `ref.jl`, `cmp.jl`, `blocks.jl`, `nocap.jl`); 4 threads
pinned to one CCX, old = fb62e1b from an exported copy, run back to back (not
interleaved; other jobs on the other CCX).

**Share first (before).** The "gap pairs" are the near pairs with a gap component: the
slot's aperture (180 chord rectangles of the cylinder at R_out, 20 µm tall, facing in)
onto the annuli's sectors (50 µm wide, 9:1 at 0.5 mm), sharing their outer arc. Timed
one pair at a time on one thread:

| | `bored_cavity()` 2,244 | `(hw = 0.8, hc = 0.25)` 1,068 | `(0.18, 0.07)` 11,360 |
|---|---|---|---|
| aperture ↔ annulus, near | 1,520 pairs × 2.28 ms | 1,536 × 2.25 ms | 1,584 × 2.26 ms |
| annulus ↔ annulus | 200 × 0.48 ms | 120 × 1.79 ms | 520 × 0.52 ms |
| all gap near pairs | 3.58 CPU-s | 3.69 CPU-s | 3.86 CPU-s |
| build, 4 threads | dense 1.96 s, H 1.05 s | dense 1.33 s, H 0.78 s | H (centroid) 5.85 s |
| share (CPU-s / 4 over the build) | dense ~45 %, H ~85 % | dense ~69 % | H ~16 % |

The cost is fixed by the gap (its 180 aperture elements), so the share falls with N but
is still above 10 % at 11k. Other near pairs of the same kind: cap ↔ wall at the corner
(1,520 pairs, 184 µs), bore ↔ LEH disk (704, ~300 µs). All were the cell rule refined to
maxdepth 8 toward the shared arc onto a 34-edge sector receiver.

**Accuracy of the old rule.** Against a converged quadrature (2-D adaptive Gauss–Kronrod
over the emitter's true patch, (φ, z) on the cylinder or (r, φ) on the sector, of the
exact point-to-polygon factor onto the sector drawn with 1,024 chords per arc; tolerance
1e-8): aperture–annulus pairs −7.1e-4 (sharing the arc), −6.2e-4, −4.4e-2, +2.3e-3 at
0.04–0.9 mm; annulus pairs across the gap +8.4e-4 (opposed), −2.1e-2 (neighbors);
wall–cap −1.0e-3 … −4.9e-3; bore–LEH up to −1.9e-2. The 16-chord arcs of the receiver
leave a slit along the shared arc (the −0.38 % of "Deposition-operator rays").

**Adopted: the double contour with exact arcs** (`_arc_pair_vf`). The contour form of
"Near flat pairs by the double contour" depends only on the two outlines, not on the
surface spanning them, so it also holds for a cylinder patch seen from inside, as long
as each element lies wholly in front of the other (every point of a sector's partner on
its front side; every point of an inward cylinder patch's partner within its radius).
Outlines of sectors and chord rectangles are arcs about the axis and straight edges:

- arc × arc (coaxial): the double integral of ρ₁ρ₂ cos(φ − ψ) ln r depends on
  u = φ − ψ only, so it is one integral of it weighted by the overlap w(u) of the two
  ranges (linear between four breakpoints). ln r² = ln(c + 4ρ₁ρ₂ sin²(u/2)),
  c = Δρ² + Δz², is singular at u = 0 on coincident circles (the shared arc) and sharply
  peaked across the 20 µm gap; ½ ln(c + ρ₁ρ₂u²) times w is integrated in closed form on
  each piece and the smooth rest by adaptive Gauss–Kronrod (a first version without the
  subtraction took 500 µs per gap pair; an x⁴ grading toward the singular ends was
  slower still);
- arc × line: the line's integral in closed form (as for flat pairs), the arc's adaptive,
  split at the azimuths of the line's ends; axial lines are perpendicular to every arc
  and radial line, so their terms vanish;
- line × line: the flat pairs' `_lnr_edges`.

Used for near pairs of a sector and a sector or an inward chord rectangle of a cylinder,
coaxial (to 1e-12), each in front of the other to 1e-9 of their size; otherwise the cell
rule as before (the wall above the slot and the annulus below it, outward cylinders,
cones, spheres, flat facets). Reciprocity is exact (one value for both directions).

Accuracy. Summed over every pair of the bored cavity at 0.5 mm (thousands of pairs, all
tiers, all taking the new rule), against the coaxial-disk closed forms: aperture → annulus
7e-13 to 1.7e-12, annulus ↔ annulus −2.9e-13, wall → cap 1.7e-14, bore → LEH disk (center
fans) 5.8e-14, cap ↔ cap 1.8e-12, LEH disk → far cap −9.8e-13. Per pair against the
converged quadrature above: ≤ 1.6e-7 on every pair checked (16 pairs, 0.02–2.1 mm apart),
the reference's own error. **H** for the method and checks.

Cost per near pair, one thread: aperture–annulus 2.29 ms → 27 µs (85×), annuli 0.48 ms
→ 81 µs, cap–wall 184 → 18 µs, bore–LEH ~300 → 20 µs, LEH disk–wall 23 → 11 µs; all such
near pairs of `bored_cavity()` 4.6 → 0.19 CPU-s.

| 4 threads, warm | old | new | entries changed | max \|ΔF\| / max \|Δ row sum\| |
|---|---|---|---|---|
| `bored_cavity()` 2,244, dense | 1.96 s | 0.86 s (2.3×) | 14,000 | 5.5e-4 / 8.7e-4 |
| `bored_cavity(hw = 0.8, hc = 0.25)` 1,068, dense | 1.33 s | 0.27 s (5×) | 10,032 | 1.3e-3 / 2.7e-3 |
| same, H-matrix | 1.05 s / 0.78 s | 0.53 s / 0.18 s | | |
| bored 11,360, H-matrix, centroid | 5.85 s | 5.15 s (−12 %) | | |
| bored 11,360, H-matrix, `:sampled` | 33.0–34.1 s | 28.6 s (−14 %) | | |
| `hohlraum(h = 0.6, capsule_h = 0.3)` 1,120 | 0.43–0.45 s | 0.33–0.35 s | 6,400 | 7.1e-4 / 1.1e-3 |
| with `sombrero = true` 1,888 | 0.49 s | 0.38–0.39 s | 6,400 | 7.1e-4 / 1.1e-3 |
| `hohlraum()` 4,116 | 2.35–2.37 s | 2.18–2.21 s | 10,176 | 5.6e-4 / 7.6e-4 |

**Not bit for bit** wherever a sector meets a coaxial sector or inward cylinder near by:
every primitive hohlraum (cap–wall corner pairs; LEH disk–wall pairs), bored cavities,
coaxial disks. That is the point of the change. Whole-matrix component exchanges (all
tiers) against closed forms, old → new:

| | 2,244 | 1,068 |
|---|---|---|
| aperture → annulus | −7.9e-4 → −7.5e-6 | −3.2e-3 → −7.5e-6 |
| annulus ↔ annulus | +1.1e-4 → −1.4e-7 | −1.4e-4 → −5.6e-8 |
| wall → cap | +2.2e-4 → +4.6e-5 | +2.7e-4 → +6.9e-5 |
| bore → LEH disk | −2.7e-4 → −5.9e-5 | −8.3e-4 → −2.0e-5 |

(The rest is the mid and far tiers.) Row sums: the hohlraum without its capsule (closed
by its LEH apertures) closes better, max |1 − Σ| 8.7e-4 → 6.3e-4 and mean 3.4e-4 →
1.6e-4 at h = 0.6 (8.5e-4 → 5.8e-4, 2.5e-4 → 1.7e-4 at h = 0.3); its wall → cap exchange
+1.2e-4 → −1.2e-4 (+2.5e-4 → +1.3e-4 at h = 0.3). With the capsule, the changed rows of
`hohlraum(h = 0.6, capsule_h = 0.3)` move mean |1 − Σ| 4.0e-4 → 4.3e-4 and max 6.3e-4 →
1.7e-3 (cap rows +2.3e-4 → +4.5e-4: the old rule's error had cancelled part of the
other tiers' on those rows); `hohlraum()` 5.5e-4 → 4.0e-4 mean. Bored cavity gap rows
change by < 1e-3 (their 5 % closure error is the slot narrower than a sample cell,
"Partial visibility"). Tests: `test/runtests.jl` "view factors" (the four closed forms
summed over every pair to 1e-11, the out-of-front case, the near pair in
`view_factors`); `hmatrix.jl`, `partial_visibility.jl`, `corners.jl`, `meshops.jl`,
`cadfit.jl` pass unchanged.

Effect on the calibrated benchmark (coordinator, same day): `examples/n221204.jl 0.8 0.4
full`, main (d19d84f) against this change, every time of the CSV: P2/P0 within 1.4e-5,
P4/P0 within 4e-6, capsule Tr and Dante (both LEH areas) within 1.4e-3 eV, the P2
crossing unchanged to 1e-3 ns (6.241 ns), capsule energy unchanged (248.7 kJ, 14.0 %).
The calibrations fitted with the old rule (wall loss, CBET, the N221204 defaults) are not
moved by it; nothing is refitted.

## N221204 late-time drive (2026-10-24) *derived*

ROADMAP 16. Two results of the N221204 benchmark (`examples/n221204.jl`) disagreed about
the late drive. (a) Default energetics: capsule Tr peaked at 301.8 eV at 5.5 ns and was
flat to 8 ns; synthetic Dante 1 (open LEH area) peaked at 306 eV at 6.25 ns; measured
Dante 1 rises to 312 eV at 8.3 ns. (b) `energetics=new`: capsule Tr 334 eV at 8 ns and an
LEH share of 33 % against HYDRA's 25.2 %. Neither Dante nor Fig. 10 was used to fit
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
the published energetics it also left the wall too reflective (LEH share 33 %).
Disagreements (a) and (b) are the same error. The pulse was not the cause (below).

**Fix** (`drive_history(...; wall_history = true)`, the default): Lindl's wall energy is
carried to a varying T as E_w ∝ [∫ (T/100 eV)^(a/b) dt]^b (a = 3.3, b = 0.62; the form
of `HeatFront`'s drive integral, exact at constant T). The loss is then the constant-T
loss at each element's effective time t_eff = ∫ (T(t′)/T(t))^(a/b) dt′, integrated
trapezoidally from the element's own incident-flux temperature (constant T before the
first step). For T ∝ t^k it multiplies E_w by (1 + 5.32k)^-0.62 (*L*; the same form for
the heat front is within 2.2 % of Hammer & Rosen's exact solution at k = 0.18). It is
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
  flat into one that rises through the peak, 4 % (default) to 13 % (new energetics)
  from 5.5 to 8 ns, and moves ~10 % of the energy from the LEH and capsule into the wall.
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
split is wall 63.6 / LEH 24.8 / capsule 11.6 % of the coupled energy (HYDRA 58.6 with
gas / 25.2 / 16.2). With the default energetics (standard mesh): capsule Tr 255.2 / 278.4 / 279.9 / 283.1 / 263.7 eV at 5 / 6 / 7 / 8 / 8.5 ns (peak 283.1 at 8.0 ns), Dante 1 over the open area 255.4 / 281.1 / 288.8 / 292.3 / 271.9; split 67.0 / 17.5 / 15.6 %.

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
| F1: + backscatter, 11 % degradation | 246.3 / 234.4 / 247.1 | 264.5 / 253.9 / 272.6 | 273.5 / 261.4 / 282.3 | 254.8 / 243.1 / 262.6 | 17.3 / 15.2 | −0.192, 6.19, +0.147, 0.049 |
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

- **After the laser droop (from 8.05 ns)** the model falls with the power (−18 eV by
  8.5 ns, −49 by 9 ns) while Dante holds 311 eV to 8.6 ns and reads 301.5 at 9 ns. A
  quasi-static radiosity with a self-similar wall loss cannot give stored wall energy
  back; re-emission from the hot surface layer of a wall under a falling drive is the
  missing physics. Not addressed.
- **Early peak**: over the open area the model is 15–25 eV above Dante at 5.0–5.5 ns
  (initial area +11 / +2): the model's rise is ~0.15 ns ahead of the measured one.
- **Capsule share** 11.6 % against 16.2 % with α_c 0.22, which was calibrated on the
  old, too hot drive. Even α_c = 0 would give only ~15 %, so 16.2 % needs a larger
  effective radius than the ablation-front schedule, not a lower albedo. The
  recalibration is still open (ROADMAP 16) and is not done here.
- **Default energetics** (Au, η 0.75, the capsule at R₀) now under-drive the peak
  (capsule Tr 283 eV at 8 ns, Dante 270 / 292 eV). `energetics=new` is the configuration
  that matches Dante's level and shape through 8 ns and the LEH share. The example's
  default is unchanged (the 2026-10-24 benchmark); switching it is a decision for the
  next round.

## Wall re-emission and the capsule radius (2026-10-24) *derived*

ROADMAP 16, the two open items of "N221204 late-time drive": (1) after the laser droop
(8.05 ns) the model's drive fell with the power while Dante 1 holds 311 eV to 8.6 ns;
(2) the capsule share, 11.6 % against HYDRA's 16.2 %. Dante 1 and Fig. 10 are checks
here, never fit. All runs `examples/n221204.jl ... full energetics=new` (the full model:
calibrated CBET, wall, lips, bubble; DU wall, `nif_conversion`, backscatter, 11 %
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
it shares no code or assumption with the closure but the material power laws).

### 1. Wall re-emission (`drive_history(...; wall_release = true)`, default)

The self-similar loss with history (`wall_history`) never returns energy: when T falls,
t_eff grows and the loss only shrinks. The closure (*L*, in `_wall_split`):

- the history loss is the heating response, evaluated at each element's running maximum
  T_h = max T (bit for bit the old loss while T never falls below it);
- a drop is the deficit D = φ(T_h) − φ(T) in the wall-energy potential φ = (T/100 eV)^a
  (a = 4 + p = 3.3), and the wall answers it with the loss's own self-similar kernel,
  superposed (Duhamel, linear in φ): net flux in F = L(T_h, t_eff) −
  f_r b C ∫ (t − t′)^(b−1) dD(t′), C = L t_eff^(1−b)/(b φ(T_h)) (constant, c σ(100 eV)⁴/b,
  for gold's unclamped Lindl loss); a drop D at t_d gives back f_r C D (t − t_d)^b by t:
  fast at first (the hot surface layer), then as the heat wave's own t^b;
- D is piecewise linear between steps and the kernel integral exact for it. The step's
  own term is stiff (Δ^(b−1)); as an emission it made `solve`'s albedo iteration diverge
  (first try: albedos of 10⁴), so below T_h it goes into the albedo, 1 − f_r C Δ^(b−1)
  φ(T)/H ∝ T^-0.7 (iterates like Lindl's), and the T-independent rest into a source
  term (`solve` gained a per-element emission from the albedo state; `sol.source`
  includes it, so B = source + αH stays exact). The effective wall albedo exceeds 1
  where the wall gives back more than it takes.
- Bookkeeping: `drive_history` books `source` as ηP and `absorbed` net of the emission;
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
/ 30 % over 0.5 / 0.5 / 1.0 ns. Energy into the wall over 8–9.5 ns (MJ/cm²; 0.149 stored
at 8 ns):

| Drop | reference | no release | f_r = 1 | f_r = 0.77 | f_r the reference needs |
|---|---|---|---|---|---|
| 5 % | +0.0193 | +0.0234 | +0.0156 | +0.0191 | 0.759 |
| 10 % | +0.0086 | +0.0177 | +0.0023 | +0.0089 | 0.779 |
| 30 % | −0.0148 | +0.0094 | −0.0291 | −0.0153 | 0.762 |

The closure is linear in f_r for a prescribed T_D, so f_r is read off exactly; the same
0.76–0.78 serves all three drops, and with 0.77 the net flux follows the reference
point by point (10 % drop, F/σT_D⁴ at 8.1 / 8.5 / 9.0 ns: +0.011 / −0.010 / +0.015
against +0.013 / −0.009 / +0.014; 30 %: +0.003 / −0.035 / −0.095 against +0.006 /
−0.032 / −0.093). The linear response (f_r = 1) overshoots by ~1.6×: a cooled surface is
more opaque and throttles the release. Without the release the wall keeps absorbing
after the drop. Caveats: static wall (the ablating wall's b = 0.62 is not tested), one
material's exponents; the history loss itself is 8 % low against the reference early
in the peak (5 ns: 0.112 against 0.122 of σT_D⁴) and within 2 % from 6 ns.

**Checks.** Rising drives: bit for bit (`test/energetics.jl`, "wall re-emission"). Step
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

- **The published albedo raises the share to 12.3 %**, not 16.2 %: Young's ablator
  albedo is below 0.16 to 6 ns (more absorption early, −1.8 eV at 5 ns, which moves the
  early-peak excess by that much) and 0.3–0.4 late, where the small radius makes it
  matter little. With α_c = 0 the share is ~15 % (before). **16.2 % is not reachable
  with absorption at the ablation front** under this drive; HYDRA's "absorbed by the
  capsule" must include absorption outside the front (the blow-off, W-doped, which
  fills the low-fill hohlraum's centre) or a larger effective radius than the front.
  No published radius of that absorbing surface was found, so nothing is fitted.
- **The fixed mesh is an energetics-exact stand-in**: shrinking the mesh itself on the
  keyframes (`moving_hohlraum(...; capsule_shrink)`, new) gives the same E_cap (11.7 /
  12.3 %) and capsule Tr within 0.3 eV, and costs nothing (the keyframes already rebuild
  the target for the moving wall: 264 against 271 s per coarse run). The approximation
  `CapsuleAlbedo`'s radius makes is in the flux *pattern*: it evaluates the capsule
  flux on R₀ and re-emits the unabsorbed share from the R₀ sphere, as if light crossed
  the shell's surroundings undisturbed.
- **P2: the shrinking sphere sees less of the late pattern**: +0.073 → +0.047 at 8 ns,
  +0.047 → +0.026 at 8.5 ns, and the minimum −0.181 → −0.170 (Fig. 10 −0.170), the
  crossing 6.24 → 6.28 ns; rms 0.041 → 0.045 (the late values move away from Fig. 10's
  +0.095 / +0.109). Early (≤ 5 ns, R/R₀ ≥ 0.86) the two agree within 0.002. Fig. 10's
  flux radius is not stated (HYDRA's is the ablation front's own), so the shrinking-mesh
  P2 is the like-for-like comparison if Fig. 10 is the flux at the front.

### Standard mesh (0.28 / 0.16)

Capsule Tr / Dante 1 initial / open area (eV); "before" is the late-time-drive notes'
standard run (no release); coarse agrees within 0.6 eV and 0.008 in P2.

| t (ns) | 5.0 | 8.0 | 8.5 | 9.0 | LEH / capsule (%) | Fig. 10 min (ns), zero, 8 ns, rms |
|---|---|---|---|---|---|---|
| before | 260.9 / 248.7 / 262.2 | 310.5 / 295.0 / 318.6 | 292.5 / 277.2 / 299.5 | 261.6 / 246.7 / 266.4 | 24.8 / 11.6 | −0.189 (5.25), 6.23, +0.077, 0.041 |
| release (`energetics=new`, now) | 261.0 / 248.8 / 262.3 | 310.5 / 295.0 / 318.6 | 299.4 / 283.4 / 306.1 | 277.0 / 260.5 / 281.4 | 25.6 / 11.7 | −0.189 (5.25), 6.23, +0.077, 0.042 |
| + shrinking mesh, Young's α(t) | 259.3 / 247.2 / 260.7 | 310.9 / 295.4 / 319.0 | 300.0 / 283.9 / 306.7 | 277.5 / 261.1 / 282.1 | 25.2 / 12.3 | −0.176 (5.25), 6.24, +0.048, 0.047 |
| measured Dante 1 / HYDRA / Fig. 10 | 237.9 | 310.1 | 311.3 | 301.5 | 25.2 / 16.2 | −0.170 (5.6), 6.67, +0.095 |

Split with the release: wall 62.7 / LEH 25.6 / capsule 11.7 % of the coupled energy,
103 kJ re-emitted (gross).

### What remains

- The late drive: with the release the model falls 12 eV from 8 to 8.5 ns where Dante
  rises 1; the hohlraum plasma's stored energy and the LEH after the laser turns off are
  the candidates (neither modeled).
- The capsule share: 12.3 % with published radius and albedo against 16.2 %; the
  absorbing radius outside the ablation front is the open input.
- The early peak (+18 to +24 eV over the open area at 5–5.5 ns, +1 to +11 over the
  initial area) is not moved by either change (the published capsule albedo takes
  1.7 eV off).
- The release closure is calibrated on a static wall; an ablating-wall (hydro) reference
  would test b = 0.62.

## Computed LEH closure (2026-10-24) *calibrated*

ROADMAP 9 ("a computed lip closure") and the late-drive suspects of ROADMAP 16 ("Wall
re-emission and the capsule radius": the LEH closing faster after the laser turns off,
the energy stored in the hohlraum plasma). The prescribed lips (`GapClosure`, Young's
ramp shifted 0.8 ns) stay available and remain the example's default; the computed one
is `hohlraum(; leh_closure = :computed)`. Dante 1 and Fig. 10 are checks, never fit.
Confidence: **H** primary text or vector data, **L** our construction.

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| Effective LEH radius vs t, 5 basecamp shots (5.75 × 10.13 mm Au, 3.37 mm LEHs, r₀ 1.685 mm, ~400 TW from 4.2 to 6.05 ns) | 1.54–1.61 mm at 3.6 ns (0.92–0.95 r₀); 1.51–1.58 at 4.0–4.2; 1.40–1.50 through the peak (mean 1.44, 0.855 r₀; shot means ±0.035); after the laser turns off (6.05 ns) N171010 1.46 → 1.44 / 1.42 / 1.43 at 7.07 / 7.44 / 7.75 ns, N171212 1.40 at 6.6–6.9 | Chen et al., PoP 27, 022702 (2020) Fig. 27, vector paths (`chen2020_leh()`, `data/chen2020_fig27.csv`); analysis error 10–20 µm (Sec. II.2b) | H |
| LASNEX on the same shots (images analysed as the data) | closes at 0.07–0.11 mm/ns after laser-off; rms against the 31 points 0.045 (f = 0.15) to 0.056 mm (f = 0.03), after laser-off 0.09–0.12 mm, all too closed | same | H (rms: L) |
| Mechanism | a dense ridge forms on the rise where the outer beams cross the lip blow-off; "once the ridge is formed the effective LEH radius is fairly stable, closing only very slowly, until the laser turns off"; then (in LASNEX) "the ridge diffuses inward"; "heating of the plasma inside the LEH region by the lasers keeps the LEH open"; the edge is set by hydrodynamics, not opacity | Chen 2020 Sec. II.2a, Fig. 5b–6 | H |
| Closure at peak, other platforms | 1.83 → 1.57 mm (SXI end value, depth 0.26 mm; 3.64 mm LEH) | Young et al., PoP 27, 082702 (2020) Fig. 3a | H |
| | ~0.84 r₀ of 1.55 mm (depth 0.25), ≤ 4 % variation through peak; closes in the low-power part, reopens at the rise | Schneider et al., PoP 22, 122705 (2015) ("Closing entrance holes" above) | H |
| SXI | time-integrated (3–5 keV clear aperture, 870 eV channel); Dante's correction f (R_LEH/R_CA)² ≈ 1.15 (+10 eV at 300 eV) | Schneider et al., HTPD 2012 (OSTI 1043636); Schneider et al., RSI 81, 10E538 (2010) (OSTI 1009651) | H |
| Hohlraum gas | He 0.3 mg/cc in N221204; 4.5 % of the incident laser energy heats the gas | n221204_benchmark.md; PRE 109, 025204 App. 3 | H |
| Gas and window heat capacity | 25–30 % of the picket energy in low-fill hohlraums, a small fraction of the peak's | Callahan et al., PoP 27, 072704 (2020) Sec. II | H |
| Cooling phase | the 0-D model misses the Dante data "when the laser power drops and the hohlraum starts to cool ... the radiation lost out of the LEH becomes more important" | Callahan 2020 Sec. II | H |

OSTI searches ("laser entrance hole", "hole closure", "LEH closure", hohlraum filling)
found no other time-resolved LEH-radius data; Hall et al. 2017 (OSTI 1361586) and the
SXI papers are time-integrated. No published lip temperature or blow-off speed.

### Model (`LipClosure`, `lip_depth`; `drive_history` per lip) *L*

The lip's plasma edge sits a depth d inside r₀ and moves as

    dd/dt = (1 − h) v(T) + h (depth − d)/hold_time,    v(T) = speed (T/100 eV)^½,

T the radiation temperature on the lip (the sound speed of x-ray-ablated gold, c_s ∝
T^½ at fixed Z̄), h the smooth step of the laser intensity through the hole (beams whose
port is on that side, after `lpi`, over π r₀²) between `hold` = (6, 17) TW/mm². Below
`hold[1]` (the foot, and as the power falls) the blow-off expands freely; above it the
edge is held at `depth`, from either side (a lip that closed further in a long foot is
pushed back, Schneider's reopening). The lip is meshed in `nr_lip` = 6 equal-area rings
down to 0.6 r₀; each ring's closed fraction is the share of its area outside r₀ − d,
so the open area is exact (to 1e-10 against the chord polygon, `test/leh_closure.jl`).
`drive_history` integrates d per lip with the trapezoidal mean of v and h over each step
(predictor: last step's T; corrector: the step's own, re-solved until d moves < 1e-5 mm),
exactly in time for the hold's relaxation, and records `lip_radius`.

### Calibration (`examples`: none; scripts in the session scratchpad) *calibrated*

The lip temperature is the model's own: a coarse basecamp history (h 0.8 / 0.4; 10.13 ×
5.75 mm, 3.37 mm LEHs, 0.908 mm capsule, `WallLossAlbedo` wall, `nif_conversion()`, the
Izumi-trends pointing, each shot's Fig. 27 power split by N170827's cone fraction from
`young2020_pulses()`, prescribed lips): lip T 120–125 eV in the foot, 170–180 at 3–3.5 ns,
275–288 through the peak, 226 / 170 / 138 eV at 6.25 / 7.0 / 7.75 ns (capsule Tr peak
292–294 eV); the five shots' lip temperatures differ by < 3 eV.

Grid fit (speed 0–0.06 mm/ns, depth 0.18–0.32 mm, hold_time 0.05–0.8 ns) of the 31 GLEH
points at hold = (6, 17) TW/mm² (107–303 TW on two 3.37 mm LEHs, the rise):

| Parameter | Fit | Range | Note |
|---|---|---|---|
| speed (at 100 eV) | 0.022 mm/ns | 0.010–0.030 (rms within +10 %); 0.020–0.024 over the hold variants | ~0.55 c_s of gold at 100 eV (Z̄ ≈ 30, 39 km/s; *L*) |
| depth | 0.240 mm | 0.235–0.245 | Young 0.26 (3.64 mm LEH), Schneider 0.25 (3.1 mm) |
| hold_time | 0.40 ns | 0.2–0.6 (degenerate with `hold`) | the rise's 0.1 mm drop over ~0.5 ns |
| hold | (6, 17) TW/mm² | (3–10, 12–22) give rms 0.032–0.034 | free at 3.6 ns (110 TW), held at peak |

rms 0.033 mm over the 31 points (rise 0.034, peak 0.034, after laser-off 0.026), mean
−0.001; LASNEX 0.045–0.056 (after laser-off 0.09–0.12). The floor is the shot scatter:
at peak N171212 sits at 1.40 and N180117 at 1.49 mm with the same pulse and
temperatures. After the laser turns off the model's edge moves in at 0.02–0.03 mm/ns
(T falling from 280 to 140 eV), measured ≤ 0.02 (N171010: −0.03 mm in 1.7 ns), LASNEX
0.07–0.11: **the measured LEH does not close faster after laser-off**, and the model,
fitted to the foot and the peak, agrees with that within the analysis error. The
integrated `drive_history` on the same basecamp with `leh_closure = :computed`
reproduces the offline integration to 2e-5 mm and the data to rms 0.033 mm.

Depth or fraction: 0.24–0.26 mm on 3.1–3.64 mm LEHs is also 0.84–0.86 r₀; the three
platforms cannot tell an absolute depth from a fraction. The depth is the model's form
(the ridge stands where the beams meet the blow-off, which does not scale with the
hole). On the Zylstra 2022 pair (§5 of "Wall material, conversion and capsule
energetics") a common depth gives ×1.143 against ×1.122 for a common fraction, the
data ≥ ×1.19: the right direction, not enough.

The thermal-band radius is smaller than the GLEH (3–5 keV, > 5 keV) edge (longer mean
free path at higher photon energy; Chen Sec. II.2); calibrating on GLEH, like Young's
SXI, makes the hole somewhat too open for the x-ray loss (unquantified, *L*).

### N221204 (check, not fit)

`examples/n221204.jl 0.8 0.4 full energetics=new` (coarse; the full model with wall
re-emission), prescribed lips (Young's ramp shifted 0.8 ns, to 0.857 r₀) against
`lips=computed`. Edge radius of the lower LEH (mm; prescribed: the radius of equal open
area):

| t (ns) | 2.0 | 4.5 | 5.5 | 8.0 | 8.5 | 9.0 |
|---|---|---|---|---|---|---|
| prescribed | 1.513 | 1.415 | 1.375 | 1.328 | 1.328 | 1.328 |
| computed | 1.513 | 1.420 | 1.320 | 1.310 | 1.310 | 1.306 |

The computed lip follows the prescribed ramp through the foot (free expansion at
120–180 eV), is held 0.24 mm in from the rise on (prescribed: still closing to 8 ns), and
stays held through the droop: the intensity through each hole is 17.8 TW/mm² at 8.6 ns
and 8.8 at 9.2 ns against the hold's 6–17, so the edge is free again only after ~9 ns.

Capsule Tr / Dante 1 over the initial / open LEH area (eV):

| t (ns) | prescribed | computed | measured Dante 1 |
|---|---|---|---|
| 5.0 | 260.9 / 248.4 / 261.8 | 261.5 / 244.0 / 261.9 | 237.9 |
| 5.5 | 279.9 / 264.4 / 280.7 | 280.9 / 259.5 / 281.1 | 263.1 |
| 6.0 | 285.8 / 271.0 / 289.4 | 286.9 / 267.0 / 290.1 | 275.5 |
| 7.0 | 297.8 / 283.5 / 305.5 | 298.7 / 281.5 / 306.2 | 291.3 |
| 8.0 | 310.6 / 294.6 / 318.2 | 311.5 / 293.2 / 319.0 | 310.1 |
| 8.5 | 299.5 / 283.1 / 305.8 | 300.4 / 281.9 / 306.6 | 311.3 |
| 9.0 | 277.0 / 260.4 / 281.3 | 278.0 / 259.0 / 282.2 | 301.5 |
| LEH / capsule share (%) | 25.5 / 11.7 | 24.8 / 11.8 | HYDRA 25.2 / 16.2 |
| Fig. 10: min (ns), zero, 8 ns, rms | −0.181 (5.25), 6.24, +0.073, 0.041 | −0.176 (5.25), 6.23, +0.073, 0.041 | −0.170 (5.6), 6.67, +0.095 |

+1 eV through the peak (the hole is 0.02 mm smaller from 5.5 ns), the LEH share 0.7
points lower (HYDRA's 25.2 % lies between), Fig. 10 within 0.005. **The late gap is not
the LEH**: 8.0 → 8.5 ns the model still falls 11 eV (capsule), Dante rises 1. Bound: to
hold the drive the losses must fall by ~4 × 3.5 % ≈ 14 % (T^~3.5 at fixed source); with
the LEH 25 % of the loss that is −56 % of the LEH area, an edge 0.44 mm further in within
0.45 ns (~1 mm/ns): 25× the fitted free speed at 300 eV, 10× LASNEX's post-laser-off
closure, and against Chen's data, which show no faster closure after the laser turns
off. Even LASNEX's 0.1 mm/ns from 8.05 ns would add ~1.5 eV by 8.5 ns.

Standard mesh (0.28 / 0.16), computed lips: capsule Tr / Dante 1 initial / open area
261.6 / 244.4 / 262.2 at 5.0 ns, 311.4 / 293.4 / 319.2 at 8.0, 300.3 / 282.0 / 306.7 at
8.5, 278.0 / 258.9 / 282.2 at 9.0 (prescribed, "Wall re-emission and the capsule radius":
261.0 / 248.8 / 262.3, 310.5 / 295.0 / 318.6, 299.4 / 283.4 / 306.1, 277.0 / 260.5 /
281.4); LEH / capsule share 24.9 / 11.8 % (25.6 / 11.7); Fig. 10 −0.183 at 5.25 ns,
crossing 6.22, +0.077 at 8 ns, rms 0.042 (−0.189, 6.23, +0.077, 0.042). Coarse agrees
within 0.3 eV and 0.007 in P2.

### Plasma-stored energy (estimate; not implemented) *L*

The quasi-static radiosity converts each step's laser power at once; energy held in the
hohlraum plasma (He gas, gold blow-off) and given back after the droop has no store.

- **He gas** (the one published inventory). 0.3 mg/cc in 0.356 cm³ (6.4 × 11.24 mm less
  the capsule) is 0.107 mg, 1.6e19 ions; fully ionized at T_e 2.5 / 3 / 5 keV (the LEH
  T_e ≤ 5.5 keV of the two-zone CBET calibration; Chen's bubble interior 3–6 keV) it
  holds 4.5 n T_e = 29 / 35 / 58 kJ. PRE 109 App. 3 books 4.5 % of the incident energy
  (~92 kJ) into the gas over the pulse, consistent with a store of that size drained into
  the wall as it goes. He radiates nothing to speak of (bremsstrahlung ~0.5 TW at
  n_e 9e19 cm⁻³); it gives its energy to the gold by conduction, at most the
  flux-limited f n_e T_e v_e = 5e12 W/cm² (f 0.05, 3 keV) over ~2.75 cm², 14 TW, so
  ~35 kJ drains over ~2.5 ns: ≲ 12 TW of x-rays after conversion against the ~100 TW
  the x-ray source loses from 8.05 to 8.5 ns (laser 448 → ~320 TW × 0.89 × 0.87), worth
  ~1–2 eV at 8.5 ns. Not the gap.
- **A fast reservoir** (the laser-heated gold corona, which converts the laser at
  Te 3–5 keV) would act as a lag τ of the source with E = P τ. Sensitivity
  (`examples/n221204.jl ... lag=<ns>`, the delivered power replaced by its first-order
  lag; coarse, prescribed lips): capsule Tr / Dante 1 open area (eV), Fig. 10.

  | lag (ns) | stored at 395 TW coupled (kJ) | 4.5 ns | 5.0 ns | 8.0 ns | 8.5 ns | 9.0 ns | Fig. 10 min (5.25 ns), zero, 8 ns, rms |
  |---|---|---|---|---|---|---|---|
  | 0 | — | 216.5 / 222.1 | 260.9 / 261.8 | 310.6 / 318.2 | 299.5 / 305.8 | 277.0 / 281.3 | −0.181, 6.24, +0.073, 0.041 |
  | 0.08 | 32 | 209.9 / 216.0 | 254.9 / 256.6 | 310.4 / 318.1 | 302.7 / 309.4 | 280.8 / 285.4 | −0.169, 6.24, +0.073, 0.041 |
  | 0.15 | 59 | 204.8 / 211.1 | 249.5 / 251.9 | 310.2 / 317.9 | 304.8 / 311.7 | 284.5 / 289.4 | −0.156, 6.24, +0.074, 0.042 |
  | 0.25 | 99 | 198.9 / 205.5 | 242.3 / 245.4 | 309.7 / 317.4 | 306.6 / 313.7 | 289.3 / 294.6 | −0.138, 6.23, +0.075, 0.044 |
  | 0.15, computed lips | 59 | 204.5 / 210.9 | 249.9 / 251.8 | 311.0 / 318.6 | 305.7 / 312.5 | 285.5 / 290.3 | −0.152, 6.23, +0.074, 0.042 |
  | measured Dante 1 / Fig. 10 | | 203.4 | 237.9 | 310.1 | 311.3 | 301.5 | −0.170 (5.6), 6.67, +0.095 |

  A lag of ~0.15 ns (~60 kJ at 395 TW coupled) brings Dante 1 at 8.5 ns to the measured
  311 eV, lowers the early-peak excess at 5 ns from +24 to +14 eV (open area), puts the
  4.5 ns readings around the data (201.7 / 211.1 eV initial / open area against 203.4),
  and moves Fig. 10's minimum to −0.156 (Fig. 10 −0.170)
  without moving the crossing. One time constant explaining both the late rise and the
  early lag (the model's rise was ~0.15 ns ahead) points at a fast plasma store of
  ~60 kJ, the size of the gas inventory but needing a ~0.15 ns return, which only the
  radiating gold corona has. Its energy content and the conversion's response time are
  not published (PRE 109 gives integrals only; Callahan 2020's 0-D model notes the
  cooling phase as where it fails), so the lag stays a sensitivity in the example, not
  a model. The energy is conserved: the run's delivered energy drops by the 19 kJ the
  lag moves past the last step. 9.0 ns stays 7–12 eV low for lags of 0.15–0.25 ns (a
  slower store, the gas, would add there).

### What remains

- A measured HYBRID-E LEH history (3.1 mm, 440 TW, a longer peak) to test the hold
  depth's independence of the LEH size and power; the thermal-band radius against the
  GLEH edge.
- The hohlraum plasma's energy content and the conversion's response time (a HYDRA
  inventory, or Dante's delay behind a step in the laser power) to turn the lag into a
  model.
- Shot scatter at peak (±0.035 mm across the five shots) is not modeled: N180117's split
  quads sit 0.05 mm more open, N171212's 1.5× picket 0.05 mm more closed.

## Backlighter foil on the line of sight (2026-10-24) *derived*

ROADMAP Next 2, second half. `backlighter_foil(los; standoff, width, height, tilt,
tilt_toward, rotation)` and `line_of_sight(target, los; ...) -> LineOfSightReport`
(`src/backlighter.jl`, `test/backlighter.jl`), used in the radiography tutorial
(`docs/literate/dropped_bundle.jl`).
- **Defaults are engineering values.** The 2 mm square gold foil and the standoffs tried
  are generic: real backlighter dimensions and standoffs are shot-specific, and for THOR
  not public. The radiography windows (two open 3 mm disks in the wall on the line of
  sight) are likewise generic. **E**
- **Placement.** Centered `standoff` mm from TCC opposite the detector, front face (the
  laser-heated side) toward the detector, normal along the line of sight unless tilted;
  width along the line of sight's θ̂, height along φ̂.
- **Line-of-sight bundle.** Every ray from a foil sample point (mesh vertices plus
  `4^level` sub-triangle points) through a sample point of the object's projected disk
  (center and rings of 8k points from θ̂), continued without end: all foil-to-detector
  x-ray paths through the object, a conservative superset of one imager's. A ray counts
  as stopped if an opaque hardware triangle lies anywhere past the foil (apertures and
  open windows are transparent) or it would leave through the foil's back. `detector =
  L` (or `Inf`) takes the imaging rays to a point detector instead. `covered` is the
  fraction of disk points whose image ray meets the foil's front. Exact for obstruction
  edges that meet the bundle at sampled extremes (tests: a shield rim grazed by the
  steepest ray, x* = D(Z + w)/(a + w) − D, to 1e-7 relative; a foil tilted by β covers
  a disk iff (w/2) cos β ≥ a; a plate at a (1 − s/L) beside converging imaging rays).
- **Shields vs backlighter beams** (tutorial geometry: 10.1 mm hohlraum, `Sombrero()`
  defaults, untuned cone pointing; beams aimed at the foil center; worst case of
  `clearance` at the 1% contour, 50 µm pointing + 25 µm positioning, ±1 mm focus; no
  wings). *derived*
  - At 8 and 12 mm standoff every one of the six candidate bundles per line of sight
    (both quads within 70° of the foil normal) lands on a shield: worst-case clearance
    −0.99 to −1.12 mm, i.e. the whole envelope, the upper quad on the upper shield, the
    lower on the lower. The earlier selection without shields (12 mm: all six clear the
    hohlraum; 8 mm: only 46 for 90-78 and 11 for 90-348) is therefore not achievable.
    At 12 mm, Q46B / Q11B (lower quads) already pass the lower shield nominally
    (+0.05 / +0.08 mm) but not in the worst case; the upper quads are fully on the upper
    shield.
  - Smallest standoff with positive worst-case clearance from hohlraum, capsule and
    both shields (bisection to 0.05 mm): bundles 46 (90-78) and 11 (90-348) 13.24 mm;
    43/14 14.45; 12/22 15.18; 11 (90-78) and 23 (90-348) 15.39; 41/13 15.95; 16/21
    16.34. The closed-form standoff at which the quad axes cross the shield plane at the
    9.2 mm rim is 11.49–14.42 mm; the worst case lies 1.6–1.9 mm beyond it.
  - The bundle choice does not change: 46/11 remain the cheapest (50° cone, ε₁ 1.1%)
    and the first to clear, now at ≥ 13.2 mm instead of 8 mm.
  - The equatorial line of sight itself is never near the shields (> 3 mm clear); at
    14 mm, 6561 of 6561 rays pass with the 3 mm windows (all stopped by the wall without
    them). At 8 mm the conservative bundle reaches ~1.73 mm from the axis at the far
    wall and the 1.5 mm-radius detector window stops 356 rays.
- **Not modeled.** The shield's mounting sector, where the sheet runs past 9.2 mm (it
  can only push the standoffs out on that side), and its ±0.12 mm dimples (well below
  the clearance margins above). The radiography windows' radiation cost is not included
  in the tutorial's drive (see the diagnostic-windows tutorial). Backlighter beam
  pointing could also be offset from the foil center or the foil moved off the
  equator; neither is explored.

## Joint multi-shot CBET calibration (2026-10-24) *derived*

ROADMAP 10 (CBET bullets). `n221204_cbet()` was fitted to one HYDRA run of one shot; its
Te(t) reaches 7.4 keV and a Te ≤ 5.5 keV refit fitted Fig. 9b as well but moved the Δλ
null from 1.2 to 0.5 Å ("Hohlraum dynamics calibrated", above). Here the two-zone model is
fitted to three shots at three detunings at once, with the published constraints as
bounds and penalties. PRE 109 Fig. 10 is not used in any fit (blind check at the end).
Confidence: **H** primary text or figure, **L** our construction.

### Data extracted

All HYDRA post-processing, not measurements: no published per-shot *measured* transfer
exists for these shots (lit review §2). Figures are rasters (no vector paths, including
PRE 106 page 2, whose 534 drawings are glyphs and legend swatches). Axes were fitted to
the tick marks on all four frame sides with a shear term (the PRE 109 rasters are sheared
by −0.005 to −0.007 px/px: the right side of each panel sits 2–3 px higher); colours were
unmixed per pixel; centrelines are stroke centroids averaged over ±0.05 ns; samples were
overlaid on the figures and checked by eye.

| Data | Source | Stored | Reading error |
|---|---|---|---|
| N220919 (Δλ 2.5 Å) per-cone input/post power and cone fraction, 0.7–7.4 ns (the black curves end at 7.4 ns in every panel) | PRE 109, 025204 Fig. 9 (black) | `pre109_fig9_n220919.csv`, `n220919_fig9()` | 0.35 TW / 0.002 seen; 1.5 TW where hidden under N221204's red (flag `h`: inner inputs before 3.5 ns, much of the outer dashed foot); ±5 TW on the steep rise |
| N210808 (1.8 Å) incident/post-CBET per cone and cone fraction, 0–9 ns | Kritcher et al., PRE 106, 025201 (2022) Fig. 1(d, e) | `pre106_fig1_n210808.csv`, `n210808_fig1()` | 0.7 TW (native pixel), 0.002 |
| N210808 and N221204 cone fraction input/post, 23° cone power | PRE 109 Fig. 1(d, e) | `pre109_fig1_n210808.csv`, `pre109_fig1()` | 0.002–0.003, ~1 TW |
| Δλ, energy, peak power, hot-spot P2, yield, HYDRA clamps for N210808 / N220919 / N221204 / N230729 | K22, PRE 109, PoP 31 Table II, Pak PRE 109 025203 | `hybride_shots.csv`, `hybride_shots()` | as published |

Checks:
- *Method.* The N221204 red curves re-digitized the same way agree with
  `pre109_fig9_n221204.csv` to rms 0.004 (cf_in) and 0.006 (cf_post), 0.6–2.7 TW per cone
  (largest on steep edges). On plateaus the new values are 0.001–0.006 lower, growing with
  time, the pattern of the shear the 2026-10-23 digitization did not correct. The existing
  file is kept (the reference fit and its tests use it). PRE 109 Fig. 1e's N221204 curves
  (a different raster) agree with Fig. 9b to rms 0.002.
- *Two N210808 runs.* The incident curves of PRE 106 and PRE 109 agree (cone fraction rms
  0.004), the post-CBET curves do not: PRE 106 is +0.04 higher on average, +0.09 at
  3.4 ns (0.266 against 0.172) and +0.026 over 5.5–7.8 ns. They are different HYDRA runs
  (*L*: PRE 106 appeared before the clamp retune from N221106 shock data, PRE 109 App. 1). PRE 109's is
  the model of Fig. 9 and Fig. 10, so it is the joint fit's N210808 target; its incident
  per-cone powers come from PRE 106 Fig. 1d (Fig. 1e of PRE 109 gives only the cone
  fraction and the 23° cone).
- *N220919 against N221204* (Fig. 9). Foot cone fraction equal to 0.006 at 0.6–2.6 ns;
  0.02–0.036 lower at 3.0–4.0 ns (outer second step ~10 % higher); peak cone fraction
  equal (0.332), every cone ~2 % higher. No resolvable rise shift (< 0.05 ns).
- *Cross-shot transfer in HYDRA* (post-CBET cone fraction, input 0.33 for all):

  | t (ns) | N210808, 1.8 Å (PRE 109 / PRE 106) | N220919, 2.5 Å | N221204, 2.75 Å |
  |---|---|---|---|
  | 5.6 | 0.596 / 0.624 | 0.634 | 0.645–0.648 |
  | 7.0 | 0.465 / 0.508 | 0.564 | 0.571 |

  So HYDRA's peak dCF/dΔλ is 0.03–0.06 per Å between 2.5 and 2.75 Å (same pulse and
  inputs: 0.014 at 5.6 ns, 0.007 at 7 ns per 0.25 Å) and 0.05 (5.6 ns) to 0.11 (7 ns)
  per Å between 1.8 and 2.75 Å (different pulses). The lit-review target 0.15–0.18 per Å
  was derived from hot-spot P2 at 2.7 µm/%CF; HYDRA's own transfer is 2–3× flatter.
- *Not found.* No N230729 transfer curves (PoP 31 shows its pulse but no post-CBET
  powers): Δλ 2.9 Å and P2 +6 µm are in `hybride_shots()` as a trend check only. No
  published LEH Te or ne for any HYBRID-E shot (PRE 109, PoP 31, PRE 106 and HE21 give
  only the clamps and multipliers); the Te bound below is the high-fill measurement
  (B16: 4.2 ± 0.2 keV at peak) and the reduced-model value of Michel 2009 (5–6 keV).

### Joint fit (`calibrate_cbet(model, shots; penalty)`, `examples/hybride_cbet_calibration.jl`)

- **Setup.** One `TwoZoneCBET` for all shots: outside multiplier and Mach, inside
  multiplier, one Te(t) table (7 knots, 1–8 ns, both zones), CH at ne/nc 0.06, Ti 1.5 keV.
  Per shot: the HYDRA input powers of its own figure (`n221204_fig9_beams` on the mock
  HYBRID-E pointing; the inner total from the input cone fraction, split evenly 23°/30°),
  its Δλ, and its timing (`delay` = its `pulse_rise` minus N221204's: −0.03 ns for
  N220919, −0.06 for N210808; the model is evaluated at t − delay). Targets: post-CBET
  cone fraction every 0.2 ns over 0.8–8.4 ns where readable (107 rows: 39 / 31 / 37).
  Equal weights.
- **Constraints.** Te ≤ 5.5 keV (bound; B16's high-fill 4.2 ± 0.2 keV at peak, Michel
  2009's 5–6 keV). Δλ null: a squared-hinge penalty (weight 10) on N221204's beams at 5.6
  and 7.0 ns requiring the inner cones not to gain at 1.0 Å and to gain at 1.55 Å.
- **New model options** (item 4). `CBET(; saturation = :min)`: S17's hard limit
  δn = min(clamp, δn_lin), a factor min(1, clamp/δn), so HYDRA's clamp values carry over.
  `pulse_rise(pulse)`: the knee of the rise (the tangent at half maximum meets the power
  1 ns earlier; scale-free, so delivered and degraded inputs agree: 4.11 ns on Pak's
  delivered pulse, 3.98 on Fig. 9's inputs). `cbet_clamp_schedule(pulse; foot, peak)`:
  the two-knot clamp table switched there (the deck's 3.9–4.1 ns on Fig. 9's inputs).
  `cbet_null(model, beams, t; receivers)`: the Δλ of zero net transfer, by bisection.

| Fit | Te bound / null | clamps | loss | rms all (N221204 / N220919 / N210808) | Mach out; mult out / in | Te 1 → 8 ns (keV) | null at 5.6 ns | dCF/dΔλ at 5.6 ns, 1–1.8 / 2.5–2.75 Å |
|---|---|---|---|---|---|---|---|---|
| `n221204_cbet()` (one shot, reference) | — | deck, smooth | — | 0.044 (0.016 / 0.030 / 0.068) | 1.23; 0.30 / 1.23 | 2.8 → 7.4 | 1.18 Å | 0.20 / 0.19 |
| free | 20 keV / none | deck, smooth | 0.053 | 0.022 (0.021 / 0.017 / 0.027) | **0.0**; 0.07 / 5.6 (degenerate) | 6.1 → **20 (bound)** | **0 Å** | 0.12 / 0.10 |
| Te-bounded | 5.5 keV / none | deck, smooth | 0.066 (150 iterations, not converged) | 0.025 (0.023 / 0.017 / 0.032) | **0.0**; 0.09 / 0.37 | 1.8 → 5.5 (bound) | **0 Å** | 0.12 / 0.16 |
| constrained, clamps started at 0.002 / 0.0075 (not binding: zero gradient, they stay) | 5.5 keV / 1–1.55 Å | fitted, hard min | 0.136 | 0.036 (0.029 / 0.019 / 0.049) | 1.50; 0.058 / 0.60 | 2.0 → 5.5 | 0.98 Å | 0.20 / 0.18 |
| **constrained + clamps (`hybride_cbet()`)** | 5.5 keV / 1–1.55 Å | fitted, hard min | 0.120 | 0.033 (0.027 / 0.016 / 0.047) | 1.54; 0.142 / 0.729 | 2.1 → 5.5 | 0.97 ± 0.10 Å | 0.197 ± 0.016 / 0.173 ± 0.022 |
| HYDRA (data) | | | | | | | — | 0.05 (7 ns: 0.11) between 1.8 and 2.75 Å, other pulse / 0.06 (7 ns: 0.03) |

- **One parameter set does not fit all three shots under the constraints.** The
  constrained fit leaves a systematic cross-shot pattern at peak (5.4–7.4 ns mean
  residual): N221204 +0.033, N220919 −0.004, N210808 −0.065 (−0.10 at 5.8 ns). The model's
  Δλ response is linear with a slope 2–3× HYDRA's, so it cannot be right at both 1.8 and
  2.75 Å. Per-cone ratios (not fitted) stay as in the reference: 44° loses more than
  50° at peak (0.39 / 0.62 against HYDRA's 0.43 / 0.64), the transfer favours 30° (2.16
  against 1.88).
- **What the data alone want** (free): no flow-driven zone (outside Mach → 0, the
  multipliers then degenerate: σ_u ~10³–10⁵), the null at Δλ = 0, a flatter response
  (0.10–0.12 per Å at 5.6 ns) and a Te of 14–20 keV, i.e. the gain scale pushed down as far
  as the bound allows. That is HYDRA's cross-shot flatness bought with an unphysical
  plasma and the opposite of NC24's statement that HYBRID-E needs 1–1.55 Å to cancel
  flow-driven CBET. N210808 is still 0.05 low at 5.6 ns.
- **It is the null, not the temperature, that the data resist.** With Te ≤ 5.5 keV and
  no null constraint the loss rises only from 0.053 to 0.066 (rms 0.022 → 0.025), still
  with Mach 0 and the null at 0 Å; adding the null takes it to 0.120–0.136 (rms 0.033–0.036).
  So the 7.4 keV of `n221204_cbet()` is not needed by any of the data; the flow-driven
  outside zone that makes the null is what HYDRA's cross-shot data do not support.
- **The tension is the null against the slope.** With the null at ~1 Å and the
  cone fraction at 0.33 there, reaching N221204's 0.65 by 2.75 Å needs ~0.18 per Å on
  average, while HYDRA's two shots on the same pulse (2.5 and 2.75 Å) differ by 0.014 at
  5.6 ns and 0.007 at 7 ns. A response that steep near the null and that flat by 2.5 Å
  needs strong saturation in Δλ. The hard clamp does bind at peak in the fit (4.8×10⁻⁴,
  inside K18's 5×10⁻⁴–2×10⁻³; the linear δn/n of the strongest pairs is ~1×10⁻³ in the
  reference plasma), but a
  binding clamp scales the gain down without flattening it in Δλ: the saturated rate goes
  as Im Γ/|Γ|, which stays near one past the null (scan: clamp 2.5×10⁻⁴ flattens the
  slope only by also cutting the transfer at 2.75 Å to 0.47). Pump depletion of the 44°
  cone flattens it only near full transfer (CF > 0.8). The 2026-10-23 single-shot fit
  (start clamp 0.0075) never reached the binding region: a clamp above every pair's
  δn has zero gradient.
- **Possible reasons** (*L*, untested): HYDRA's N210808 and N220919 runs differ from
  N221204's in more than Δλ (multipliers "adjusted in a common way" but per-shot laser
  degradation, as-delivered powers, and for N210808 a different pulse);
  the cross-shot flatness may be a property of HYDRA's inline CBET in a filling hohlraum
  (the transfer region moves with Δλ), which a fixed two-plane model cannot have; or
  NC24's 1–1.55 Å is a hot-spot-symmetry statement (data-driven, with other pulse
  changes), not the null of the transfer itself.
- **Uncertainty (Laplace, Gauss–Newton at the constrained optimum; σ from the residuals,
  0.035, model error included; the Te knot at its bound and the penalty held fixed).**
  Outside multiplier ×/÷ 1.96, inside ±18 %, Mach ±0.38, Te knots ±6–14 %, peak clamp
  ±27 %, foot clamp unconstrained (σ_u 12: not binding, the foot is fixed by Te and the
  multipliers). Strong correlations: outside multiplier with Mach (−0.89), inside
  multiplier with the early Te knots (+0.75–0.83) and with the peak clamp (−0.83). Derived:
  null 0.97 ± 0.10 Å (the penalty holds it at the band's lower edge), dCF/dΔλ
  0.197 ± 0.016 (1–1.8 Å) and 0.173 ± 0.022 per Å (2.5–2.75 Å) at 5.6 ns. HYDRA's
  0.06 between 2.5 and 2.75 Å is > 5 σ away: the discrepancy is the model's structure,
  not its parameters.
  A profile in the null (a fit pinned at each null) is the next step; the two ends are in
  the table (null 0: loss 0.053, or 0.066 with the Te bound; null ~1 Å: loss 0.120), each
  fit 15–90 min on the shared 4 threads.
- **Not adopted as the default.** On its own targets the joint fit is better on N220919
  (rms 0.030 → 0.016) and N210808 (0.068 → 0.047) and over all rows (0.044 → 0.033), but
  worse on N221204 (0.016 → 0.027), the example's default shot. `n221204_cbet()` stays the
  example's default; `hybride_cbet()` is the `joint` transfer option.

### Blind check: PRE 109 Fig. 10 per shot

`examples/n221204.jl 0.8 0.4 full <transfer> times calibrated <shot>` (coarse mesh,
current defaults: calibrated dynamics, wall history and re-emission on, default
energetics), `calibrated` = `n221204_cbet()` (before), `joint` = `hybride_cbet()` (after;
N210808's proxy pulse rises 0.2 ns early, so its tables are shifted by −0.2 ns). The
proxies are as before: N220919 on N221204's delivered pulse, N210808 the DN24 pulse
shifted. Nothing was fitted to Fig. 10.

| Shot | Fig. 10: min (ns), zero, 8 ns | before: min, zero, 8 ns, rms | after: min, zero, 8 ns, rms | post-CBET CF at 5.5 / 7 ns, before → after |
|---|---|---|---|---|
| N221204 (2.75 Å) | −0.170 (5.60), 6.67, +0.095 | −0.225 (5.25), 6.26, +0.143, 0.046 | −0.209 (5.25), 6.31, +0.141, 0.043 | 0.687 / 0.587 → 0.663 / 0.586 |
| N220919 (2.5 Å) | −0.145 (5.60), 6.66, +0.094 | −0.195 (5.25), 6.22, +0.144, 0.045 | −0.185 (5.25), 6.24, +0.142, 0.044 | 0.634 / 0.546 → 0.619 / 0.551 |
| N210808 (1.8 Å) | −0.141 (5.50), 6.68, +0.084 | −0.126 (5.25), 6.08, +0.125, 0.071 | −0.137 (5.25), 6.13, +0.123, 0.071 | 0.478 / 0.428 → 0.491 / 0.440 |

- The joint calibration moves every shot's swing by ≤ 0.016 and the crossing by
  ≤ 0.05 ns, each toward Fig. 10 (minimum within 0.04 of Fig. 10 for all three now;
  N210808's −0.137 against −0.141). The cross-shot ordering of the minima is HYDRA's
  (N221204 deepest), and the spread shrinks from 0.099 to 0.072 (Fig. 10: 0.029).
- Fig. 10's remaining misfit (the crossing ~0.4–0.55 ns early, +0.03–0.05 at 8 ns, the
  N210808 foot +0.10–0.15 above at 4–4.5 ns) is not CBET's: it is the same with either
  calibration (bubble speed, proxy pulses; "Hohlraum dynamics calibrated", above).

### What remains

- ~~A per-shot profile of the null and an explicit Δλ-saturation mechanism (or a
  transfer region that moves with Δλ) to reconcile the null with HYDRA's flat cross-shot
  response; or evidence that NC24's 1–1.55 Å is not the transfer null.~~ Done: NC24's
  1–1.55 Å is the Δλ of round hot spots, not the null ("Joint multi-shot CBET
  calibration: the null", below).
- The 23°/30° split (the joint fit still favours 30°), an inner-quad SBS proxy
  (N220919's limit) and outer depletion reporting; N230729 (2.9 Å) has no published
  transfer curves, so its +6 µm overshoot stays a trend check.
- The N221204 Fig. 9 digitization could be replaced by the shear-corrected re-digitization
  (differences ≤ 0.006 in cone fraction); not done, to keep the reference fit reproducible.

## Joint multi-shot CBET calibration: the null (2026-10-25) *derived*

Follow-up to the previous section's open question: is NC24's 1–1.55 Å the Δλ null of
the transfer? Confidence: **H** primary text, **L** our construction.

### What NC24's 1–1.55 Å is

- **The Δλ of round hot spots, not the transfer null (H).** NC24 (Nat. Commun. 15, 2975,
  p. 3): "In experiments, 1 ≤ Δλ ≤ 1.55 Å was required to minimize P2 in HYBRID-E designs
  and 0.5 Å Δλ was required in I-Raum experiments mostly to negate the flow-induced
  CBET." P2 is that of the hot spot's 17 % self-emission contour. The same paper says
  HYBRID-E "principally uses wavelength separation … to compensate for the inner beam
  absorption in the bubble and the flow-induced CBET": at a round hot spot the inner
  cones must *gain*, so HYBRID-E's transfer null lies below 1 Å.
- **It is design-specific.** The round-hot-spot Δλ moved with the design: 1.55 Å
  (N210207/N210307, 3.64 mm LEH, 1050 µm), 1.8 (N210808, 3.1 mm LEH), 2.5–2.75 Å
  (N220919/N221204, longer pulse, thicker ablator); HE21: LEH 4.0 → 3.64 mm needed
  +0.5 Å. A transfer null would not track capsule and pulse changes.
- **The null itself (L for HYBRID-E).** The only published statement is I-Raum's:
  0.5 Å "mostly" negates the flow-induced CBET in a hohlraum built to need little
  bubble compensation, so the null is ≤ 0.5 Å there; flow-driven inner → outer transfer
  at Δλ = 0 exists (P20's HYDRA; NC24), so it is > 0. HYBRID-E has the same fill
  (0.3 mg/cc). Taken here as a band of 0.25–0.5 Å (the lower edge is our choice).
- `hybride_cbet()`'s constraint (null in 1–1.55 Å) therefore over-constrained the joint
  fit; the docstrings of `cbet_null` and `hybride_cbet`, and lit review §3, now say so.

### Is the cross-shot flatness physical?

- **The measured response is steep and unsaturated (H).** HE21: "no indication of
  transfer saturation", 40–50 µm P2 per Å. NC24: N201101 → N210307 (a repeat, 1.75 →
  1.55 Å) 15.7 ± 0.7 → 0.68 ± 2.8 µm, ≈ 75 µm/Å. PRE 109: N220919 → N221204, −12 → ~0 µm
  for +0.25 Å (≈ 48 µm/Å; the foot cone fraction also changed), Δλ set from "a
  data-driven model … together with an experimental sensitivity curve hotspot P2 vs Δλ"
  because "late-time laser beam propagation through the plasma filled hohlraum is
  difficult to model". Via 2.7 µm/%CF that is 0.15–0.18 per Å.
- **HYDRA's cross-shot transfer does not reproduce it.** Its post-CBET cone fractions of
  N220919 and N221204 differ by 0.014 at 5.6 ns (0.007 at 7 ns), ~4 µm of hot-spot P2 at
  2.7 µm/%CF against the ~12 µm measured. Read as a Δλ response, HYDRA's 0.05 per Å
  between 1.8 and 2.75 Å (5.6 ns) puts the null at 1.8 − (0.596 − 0.33)/0.05 ≈ −3.5 Å:
  no Δλ response with a null ≥ 0 passes through HYDRA's three shots. Their differences
  carry per-run settings (PRE 109: "the hohlraum model was adjusted … by modifying the
  CBET saturation" after N221106; K22: N210808's peak flux P2 × 0.7; the PRE 106 and
  PRE 109 N210808 runs differ by 0.04–0.09 in post-CBET cone fraction); *L* as to which
  run used which setting.
- So mechanisms that flatten the Δλ response (Δλ-dependent saturation, a Mach
  distribution, outer-cone depletion) would fit HYDRA's runs against the measured
  sensitivity, and none was implemented. A scan on `hybride_cbet()` at 5.6 ns
  (N221204's beams; dCF/dΔλ per Å) shows that none of them flattens the high end
  selectively anyway:

  | Variant | 1–1.8 Å | 2.5–2.75 Å |
  |---|---|---|
  | `hybride_cbet()` | 0.197 | 0.173 |
  | inside Mach −0.5 / +0.5 | 0.193 / 0.150 | 0.165 / 0.170 |
  | inside Mach uniform in 0–1 / 0–2 (5 sub-zones, multiplier / 5 each) | 0.151 / 0.067 | 0.171 / 0.138 |
  | Te × 0.7 / × 0.5 (ion-acoustic resonance nearer) | 0.264 / 0.345 | 0.209 / 0.183 |
  | peak clamp 2.5×10⁻⁴ / 1×10⁻⁴ (stronger saturation) | 0.103 / 0.041 | 0.099 / 0.039 |
  | inside Mach −1 (pump depletion: CF 1.00 at 2.5 Å) | 0.148 | 0.001 |

  A Mach distribution acts like its mean, saturation scales both ends together, a
  nearer resonance steepens the low end, and depletion flattens only where the outer
  cones are emptied (CF > 0.8, beyond every shot).

### Refit with the corrected constraint (`hybride_cbet_iraum()`)

`examples/hybride_cbet_calibration.jl iraum` (now the default variant): as `clamp`
(fitted hard-minimum clamps on the pulse-driven schedule, Te ≤ 5.5 keV, the same 107
rows and start) with the null penalty at 0.25–0.5 Å; `clampfree` has no null penalty.
2 threads each, run concurrently.

| Fit | null band | loss | rms all (N221204 / N220919 / N210808) | Mach out; mult out / in | Te (keV) at 1, 4, 5.5, 8 ns | clamps foot / peak | null at 5.6 / 7.0 ns | dCF/dΔλ at 5.6 ns, 1–1.8 / 2.5–2.75 Å |
|---|---|---|---|---|---|---|---|---|
| `hybride_cbet()` (`clamp`) | 1–1.55 Å | 0.120 | 0.033 (0.027 / 0.016 / 0.047) | 1.54; 0.142 / 0.729 | 2.1, 2.8, 3.4, 5.5 | 1.4e-3 / 4.8e-4 | 0.98 / 1.09 Å | 0.197 / 0.173 |
| **`hybride_cbet_iraum()`** (105 iterations, converged, 45 min) | 0.25–0.5 Å | **0.059** | **0.024 (0.024 / 0.016 / 0.028)** | 1.49; 0.072 / 1.60 | 3.2, 4.5, **2.5**, 5.5 | 9.1e-4 / 1.4e-4 | 0.24 / 0.29 Å | 0.134 / 0.137 |
| `clampfree` (150 iterations, not converged, 61 min) | none | 0.055 | 0.023 (0.023 / 0.016 / 0.026) | **0.0**; 0.011 / 1.39 | 3.1, 4.4, 2.4, 5.5 | 1.0e-3 / 1.5e-4 | 0 Å | 0.122 / 0.127 |
| `n221204_cbet()` (one shot) | — | — | 0.044 (0.016 / 0.030 / 0.068) | 1.23; 0.30 / 1.23 | 2.8, 3.7, 5.5, 7.4 | deck, smooth | 1.18 / 1.26 Å | 0.20 / 0.19 |
| HYDRA, cross-shot | | | | | | | — | 0.05 (1.8 → 2.75 Å, other pulse) / 0.06 |

- **The conflict was the constraint.** With the null in I-Raum's band the joint loss
  halves (0.120 → 0.059) and is within 0.004 of the unconstrained fit (0.055, still
  descending): the data are indifferent between a null of 0 and of 0.25 Å and resist
  1 Å. Every shot improves over `hybride_cbet()`, N210808 most (rms 0.047 → 0.028;
  peak-mean residual −0.065 → −0.030). The flow-driven outside zone survives (Mach 1.49;
  cone fraction 0.305 at Δλ = 0 against the input 0.33), unlike the unconstrained fit
  (Mach → 0).
- **Slope.** 0.13–0.14 per Å at 5.6 ns, linear from 1 to 2.75 Å (0.09 at 7 ns, 0.06 at
  8 ns): a little below the measured-P2 0.15–0.18 and 2× HYDRA's cross-shot 0.06.
  Model cone fractions at 5.6 ns on each shot's beams: 0.674 / 0.634 / 0.543 against
  HYDRA's 0.648 / 0.634 / 0.596, so N210808's peak stays 0.03–0.05 low.
- **Per-cone ratios** (not fitted): N221204 at 5.6 ns, 23/30/44/50, 1.92 / 2.09 /
  0.44 / 0.54 against HYDRA's 1.92 / 1.88 / 0.43 / 0.64 (the 23° ratio is now right;
  30° still favoured, 50° loses too much).
- **What is not physical about it.** Te dips to 2.5 keV at 5.5 ns between 4.5 keV at
  4 ns and 3.0 at 6.5 ns (B16: Te follows the drive and peaks late, 4.2 keV), and the
  peak clamp 1.4×10⁻⁴ is below HYDRA's deck values (5×10⁻⁴–7.5×10⁻³; our ion-wave
  coefficient is not HYDRA's, so the values compare only roughly). Both trade against
  the inside multiplier (correlations −0.92 to −0.97): with the hard clamp binding, the
  peak gain is set by clamp × multiplier and Te enters through Im Γ/|Γ|. The
  unconstrained fit lands in the same place, so the dip is the data's, not the null's.
- **Uncertainty** (Laplace, Gauss–Newton, σ from the residuals 0.025, penalty held
  fixed): outside multiplier ×/÷ 4.9, inside ±31 %, Mach ±0.83 (−0.84 with the outside
  multiplier), Te knots ±7–15 %, peak clamp ±29 %, foot clamp unconstrained. Derived:
  null 0.24 ± 0.20 Å (the data alone barely fix it; the penalty holds it at the band's
  lower edge), dCF/dΔλ 0.134 ± 0.011 (1–1.8 Å) and 0.137 ± 0.011 per Å (2.5–2.75 Å).

### Blind check: PRE 109 Fig. 10 per shot

As in the previous section (`examples/n221204.jl 0.8 0.4 full <transfer> times
calibrated <shot>`, coarse mesh, current defaults), `iraum` = `hybride_cbet_iraum()`.
Repeated on the current code, `calibrated` (N221204) and `joint` (N221204, N220919)
reproduce the previous section's rows to the last digit. Nothing was fitted to Fig. 10.

| Shot | Fig. 10: min (ns), zero, 8 ns | rms: `n221204_cbet()` / `hybride_cbet()` | `hybride_cbet_iraum()`: min, zero, 8 ns, rms | post-CBET CF at 5.5 / 7 ns: `hybride_cbet()` → iraum |
|---|---|---|---|---|
| N221204 (2.75 Å) | −0.170 (5.60), 6.67, +0.095 | 0.046 / 0.043 | −0.215 (5.50), 6.33, +0.139, **0.040** | 0.663 / 0.586 → 0.674 / 0.558 |
| N220919 (2.5 Å) | −0.145 (5.60), 6.66, +0.094 | 0.045 / 0.044 | −0.197 (5.50), 6.28, +0.140, **0.044** | 0.619 / 0.551 → 0.640 / 0.535 |
| N210808 (1.8 Å) | −0.141 (5.50), 6.68, +0.084 | 0.071 / 0.071 | −0.168 (5.25), 6.18, +0.120, **0.052** | 0.491 / 0.440 → 0.536 / 0.462 |

- rms against Fig. 10 improves on N221204 (0.043 → 0.040) and N210808 (0.071 → 0.052)
  and is unchanged on N220919; the upward crossing moves 0.02–0.05 ns later (toward
  Fig. 10) and the N221204/N220919 minima to 5.50 ns (Fig. 10: 5.60). The minima are
  deeper than `hybride_cbet()`'s by 0.006–0.031, now 0.03–0.05 below Fig. 10 (more
  transfer at 5.5 ns, less after 7 ns). Cross-shot spread of the minima 0.047
  (`hybride_cbet()` 0.072, Fig. 10 0.029), in HYDRA's order.
- The rest of Fig. 10's misfit (the early crossing, +0.03–0.05 at 8 ns) is the same with
  every calibration, as before.

### Adopted, and what remains

- `hybride_cbet_iraum()` (exported) and the `iraum` transfer option of
  `examples/n221204.jl`; `hybride_cbet()` and `n221204_cbet()` are unchanged. The
  calibration example's default variant is now `iraum` (loss 0.059 against `clamp`'s
  0.120 on the same targets). `examples/n221204.jl` keeps `n221204_cbet()` as its
  default: on N221204's own Fig. 9 it is still the best (rms 0.016 against 0.024).
- Left: the Te dip and the low peak clamp (a drive-following Te(t) prior, e.g. B16's
  shape, would test whether the data need them); the 23°/30° split; HYDRA's per-run
  settings cannot be separated from Δλ in the published curves.

## Backlighter conversion efficiency and the photometric chain (2026-10-24) *engineering*

ROADMAP Next 12. `kshell_line_energy`, `kshell_ce_data`, `kshell_ce_model`, `kshell_ce`,
`mass_attenuation`, `filter_transmission`, `backlighter_signal` (`src/backlighter_signal.jl`,
`test/backlighter_signal.jl`, `data/kshell_ce.csv`); used in the radiography tutorial
(`docs/literate/dropped_bundle.jl`, "The radiograph's signal") and `examples/backlighter.jl`.
All of it is public literature plus generic detector geometry; no facility-specific
backlighter design is used.

**Data** (`data/kshell_ce.csv`, 42 rows, 37 fitted). CE = He-α (w line plus satellites)
energy into 4π over the *incident* laser energy, 351 nm, planar foils; every source
converts J/sr at one view angle to 4π assuming isotropy.
- *Workman & Kyrala, LA-UR-01-4155 (2001), OSTI 783558* (OMEGA, 1 ns square, solid flat
  targets, no picket): Cr, Fe, Zn, Ge irradiance scans (Figs. 2–5) and Ni (Fig. 6), 24
  points **digitized** from the rendered figures (reading ~3% in CE). 550 eV band about
  He-α (may include part of Ly-α); crystal reflectivities estimated, so the absolute
  scale is uncertain; relative ±10%, reproducibility 20%; Ge's two lowest points ±50%.
  Peaks: Cr > 3% at ≤ 1e15 W/cm² (not resolved), Fe ~3% at 3.5e15, Zn ≥ 0.84% at
  ≤ 4.4e15, Ge 0.22% near 4e16 (unsmoothed tight-focus irradiances, uncertain). They
  report CE ∝ hν^-(2.5–4) across lines and I^-0.3 (Cr) to I^-0.5 (Zn) past the peak.
- *Barrios et al., HEDP 9, 626 (2013), OSTI 1237551*, Table 3 (exact): NIF, 1 ns picket
  3 ns before a 2–5 ns drive, 4–8 beams at 40°, CPP spots. Zn 0.62–1.08% (0.57–3.5e15),
  Ge 0.52–1.54% (0.57–2.96e15, peak 1.5e15), KBr 0.29% (half the emitters; not fitted),
  Zr 0.40, 0.32%.
- *Barrios et al., RSI 85 (2014), LLNL-JRNL-656005, OSTI 1162253*, text (exact): Ge
  one-sided 1.0 ± 0.3% at 2.9e15 (two-sided 1.3 ± 0.2, not fitted); Zr 0.18 ± 0.10 and
  0.15 ± 0.07% (1.94, 4.5e15; taken as a reanalysis of the 2013 Zr shots, whose rows
  are kept unfitted, *inferred*); Zr without picket 0.09 ± 0.05% (4.7e15, put in the
  single-pulse class); Nb 0.15, 0.18% (17, 23e15); Mo 0.16 ± 0.08% (9.5e15). Pre-pulsed
  foils follow hν^-(3–4) (their Fig. 3, with Workman's bounds).
- Context, not fitted: Kemp, APS-DPP 2025 slides (OSTI 3389300) compile foils (0.2–0.5%
  at 7–10 keV), pre-pulsed foils (1–1.5% at 9–10 keV, 0.3–0.4% at 16 keV), foil-lined
  cavities and under-dense targets (2× better). Pre-exploded Ti foils reached up to ~8%
  (Girard et al., PoP 12 (2005); seen only through abstracts and citing papers, not
  checked against the paper). Gas, aerogel, foam and cavity sources are out of scope.

**Model as first fitted** (2026-10-24; superseded by the weighted refit on the extended
compilation, "Weighted He-α CE compilation (2026-10-26)" below, which also replaces the
uncertainty and pulse-length bullets here): least squares on ln CE with equal weights, `ln CE = a_c + b L + c_c X + d X² + e X L`, L = ln(hν/10 keV),
X = ln(I/1e15 W/cm²), class c ∈ {foil (one pulse), prepulsed}. Coefficients: a = −5.694
(foil), −4.876 (prepulsed); b = −3.92; c = 0.244, 0.209; d = −0.098; e = 0.543. Residual
scatter **0.41 in ln CE (×1.51)** with 30 degrees of freedom (foil class rms 0.38,
prepulsed 0.35). The optimum irradiance per line, `exp(−(c_c + e L)/2d)` × 1e15 W/cm²:
foil Fe 1.2e15, Cu 2.1e15, Ge 3.8e15, Mo 1.8e16; prepulsed Zn 2.2e15, Ge 3.1e15, Zr 1.1e16.
Peak CE: foil Fe 1.6%, Cu 0.71%, Ge 0.36%; prepulsed Cu 1.6%, Ge 0.78%, Zr 0.20%,
Mo 0.15%. The peak falls as hν^-3.5 (foil Fe→Ge), inside the published hν^-(2.5–4).
- **Shortfalls, stated.** The peaks are flatter than the data: Workman's Fe peak (3%) and
  Barrios's Ge peak (1.5% at 1.5e15) sit ~2× above the fit, and the fitted prepulsed Ge
  optimum (3e15) is ~2× above Barrios's (1–1.5e15). The quadratic in ln I is the
  simplest form with an optimum; the scatter absorbs the rest. *confidence: medium
  inside the data (foil 5.7–16.3 keV, 0.35–76e15; prepulsed 9.0–18.0 keV,
  0.57–23e15), low beyond.*
- **Uncertainty.** `kshell_ce` returns σ = √(scatter² + vᵀ C v), C the coefficient
  covariance, v the basis at the point: 0.43–0.47 inside the data, 0.49 for Ti
  (single pulse), 0.59 for Ag prepulsed, 0.74 for Ag single pulse (all at 2e15), 1.0 for Ti
  prepulsed at its optimum (16% ×/÷ 2.7, i.e. unknown). `extrapolated` flags any point outside its
  class's line, irradiance or (if given) pulse range.
- **Pulse length.** Not in the model: the foil class is 1 ns (plus one 2–5 ns NIF Zr
  point), the prepulsed class 2–5 ns, so τ and the picket are confounded. `pulse` only
  sets the flag (foil 1–5 ns, prepulsed 2–5 ns). *E*
- **Line energies** (`kshell_line_energy`): Ti 4.750 … Ge 10.279, Kr 13.115 as commonly
  quoted; Br, Zr, Nb, Mo, Ag from `10.204 eV (Z − 0.42)² [1 + 0.0011 (Z − 22)]` fitted
  through Ti, Cu, Ge, Kr and Zr (≤ 0.5%). *derived*

**Attenuation** (`mass_attenuation`, `filter_transmission`): NIST X-Ray Mass Attenuation
Coefficients (Hubbell & Seltzer), μ/ρ with coherent scattering, 3–40 keV, for Be, Al, Ti,
Fe, Cu, Mo, Ta, Mylar and polyethylene; Kapton uses Mylar's μ/ρ at ρ = 1.42 (a few %).
Log-log interpolation between the tables' nodes (1–5 keV apart; a few % away from
edges), edges exact (Ti K 4.966, Fe K 7.112, Cu K 8.979, Mo K 20.00, Ta L 9.881/11.14/
11.68 keV). Beer–Lambert at normal incidence. A transmission estimate, not a channel
response.

**Chain** (`backlighter_signal(target, deposition, los; line, pulse, distance, ...)`):
1. Laser energy: each foil element's deposited power (Briareus deposits all incident
   power; CE is per incident energy, so backscatter is inside CE) × a flat-top `pulse`.
2. CE at the power-weighted mean irradiance (`irradiance = :element` evaluates each
   element at its own; σ stays the mean's), or a given `ce`.
3. Angular distribution: `:isotropic` by default (how the CE data were converted to 4π;
   half the energy toward the front), `:lambertian` (4 cos θ, front only), or `f(μ)`
   with ∫ f dΩ = 4π. Elements facing away or with hardware (any surface part but the
   foil and `object`) on the path contribute nothing.
4. Point projection: fluence Σ E_e f/(4π r_e²). Pinhole imager (diameter, distance p;
   detector at q = distance − p behind it): each image point gets R A_p/q² with
   R = E_e f/(4π A_e cos θ) the radiance of the foil element seen through the pinhole and
   a sample point of the object's projected disk; the mean and range over the disk.
5. Filters at the line energy; 6. photons = fluence/hν per mm² and per pixel; the CE's
   1σ gives `photons_lo`/`photons_hi`. The object does not attenuate (open field);
   continuum, Ly-α and detector response are not modeled.
- Tests (closed forms): point projection Σ E_e/(4π r²) by hand and E_x/(4π(L + D)²) to
  1e-5; Lambertian pinhole image E_x A_p/(π A q²) exactly at every image point; energy
  through the pinhole = fluence × image area (magnification q/(p + D)); hemispherical
  integrals of per_sr E_x/2 (isotropic) and E_x (Lambertian) to 1e-10; filters vs NIST
  nodes, edges, power-law interpolation, layer products; the fit's weighted χ²/dof = 1
  and weighted residual sum zero per class, scatter = τ ⊕ rms stated error, the
  held-out check of the next section, the optimum a maximum, the hν slope in 2.5–4.5.

**Tutorial numbers** (radiography tutorial; bundle 11's two quads, 16 TW for 2 ns with a
picket, aimed at the center of the 2 mm Ge foil at 14 mm; generic imager: 20 µm pinhole
at 200 mm, detector at 2 m, 50 µm pixels, 250 µm Be + 50 µm Kapton). *derived*; CE and
photon numbers from the saturating-peak form ("He-α CE: saturating peak (2026-10-26)",
tutorial rerun as a script), in brackets the weighted refit of the earlier form
(2026-10-26) and the 2026-10-24 fit.
- 13.6 of 16 TW land on the foil, mean irradiance 5.1e14 W/cm² (just below the NIF
  data; extrapolated), against Ge's fitted optimum 2.0e15 [2.5e15; 3.1e15]: the CPP spots
  at best focus fill most of the foil.
- CE 0.54% ×/÷ 1.47 [0.49% ×/÷ 1.47; 0.56% ×/÷ 1.62]: ~146 J of He-α into 4π [133; 154],
  ~11.6 J/sr toward the detector (= E/4π, asserted), filters 0.95, **~455 photons per
  pixel** averaged over the capsule's image (~310–670 at 1σ) [415, 280–610; 480,
  295–771]; across the image 2e3–5e5 photons/mm², the foil's radiance following the spots.
- At the same irradiance: Cu 1.33% → ~1330 photons/pixel [1320; 1340], Zn 0.99% → ~930
  [880; 940], Ge ~455, Zr 0.05% → ~26 [28; 43].
- `examples/backlighter.jl` (OMEGA Ti foil, 1 ns, point projection at 300 mm behind 25 µm
  Be), rerun: 8.0 kJ on the foil at 4.6e14 W/cm², below the fitted lines (Cr is the
  lowest), so extrapolated: CE 3.0% ×/÷ 1.52 [8.9% ×/÷ 1.53; 6.7% ×/÷ 1.78], ~240 J into
  4π, ~19 J/sr, ~2.6e11 photons/mm² [710 J, 56 J/sr, 7.8e11; 531 J, 42 J/sr, 5.8e11]. The
  held-out OMEGA Ti foil (~1.6% at 1.3e15) is 1.1σ below the fit there; the example's
  comment says to pass a measured `ce`.

**Not modeled / next.** Pulse-length and spot-size dependence (data confounded);
backside emission of thin foils; the line's opacity and limb brightening (the isotropic
assumption); continuum and higher lines in a filter's pass band; detector response
(image plate, CCD QE) and the object's transmission. The low-Z end: see the next
section.

## Weighted He-α CE compilation (2026-10-26) *engineering*

ROADMAP 12, "a weighted, larger low-Z compilation (Ti, V)". `data/kshell_ce.csv` grows
from 42 to 74 rows (68 fitted, from 37); `kshell_ce_model` becomes a weighted fit with a
`tau` field. The fit numbers follow from the CSV; `test/backlighter_signal.jl` asserts
the weighting and the C19 held-out check.

**Search** (2026-10-26; OSTI API, OSTI full texts, publisher abstracts, the local paper
folder). Ti/V/Sc/Ca/Cr/Mn He-α CE at 351 nm on planar foils is scarce in public, citable
form. Found and used:
- *Coppari et al., RSI 90, 125113 (2019), LLNL-JRNL-757490, OSTI 1632368* (C19): OMEGA and
  OMEGA EP, 12.5 µm Fe, Cu, Ge foils, 1 ns flat-top (Fe 0.6–2 ns), with or without a
  ~1 ns pre-pulse at ~2e14 W/cm² 3–4 ns ahead; He-α complex (continuum subtracted),
  calibrated flat-crystal spectrometer and image plate, J/sr × 4π (isotropic; Fe s77768:
  ~20 J of 2.2 kJ), over the **main-pulse** energy (the pre-pulse, ~20–25% of the total,
  not counted); irradiance on the elliptical footprint of the nominal spot. Stated
  errors: CE 20%, spot ±50 µm. 31 points **digitized** from Fig. 16 (600 dpi render;
  reading ±0.0003 in CE, ±0.02e15 in I); Fe irradiances and pulse lengths from Table I
  (printed with its facility-to-irradiance columns one row below the shot numbers;
  realigned by the paper's own checks, s72639 = 7.8 kJ at 2e15 and s77768 = 0.3e15 at
  0.6 ns); Table II's three values exact (Ge 0.006 at 1.3e15 with pre-pulse, Cu 0.01 at
  1.1e15, Fe 0.008 at 0.3e15). Peaks: Fe ~3% (0.8–2e15), Cu ~2% with pre-pulse, Ge ~1%.
  Omega shots irradiate both sides of the foil, EP shots one (Barrios 2014's Ge pair:
  two-sided ≈ 1.3× one-sided); kept in the fit, with a C19 offset tested below.
  *confidence: medium (digitized; two-sided; main-pulse normalization).*
- *Ti, OMEGA* (P16): 5–6 µm Ti foil, one side, 1 ns at 1.32e15 W/cm², no pre-pulse:
  1.8e-3 J/J/sr ±30% in 4–6 keV (DMX at 68.7° from the normal), as quoted in the CEA
  review, Phys. Plasmas 23, 040501 (2016), Sec. II (the original, Girard et al., PoP 12,
  092705 (2005), OSTI 20764547 abstract, gave 1.5% in 2π above 4 keV before the
  recalibration the review mentions). ×4π = 2.26% for the band; ×0.7 for the He-α share
  (the review's share for its Cu and Ge liner bands, ~70–75%, taken for Ti; *inferred*,
  ±0.2) = **1.58% ± 0.66** (30% ⊕ 29%). Held out (`fit = 0`) as the low-Z check.
  *confidence: low (one point; band share inferred).*
- Context, not in the CSV (conditions or convention not stated in what is public):
  Tanabe et al., APL 93, 051505 (2008) (Gekko XII, 351 nm): planar Ti foil 1.4 ± 0.9% in
  4–6 keV (TiO₂ nanofiber cotton 3.7 ± 0.5%); Hu et al., PoP 14, 033103 (2007)
  (Shenguang II): 5.8 µm Ti foil ~1.75% (half the 3.5% of 3.5 µm foils); Bullock, Landen
  & Bradley, UCRL-JC-139282 (OSTI 15006503) (OMEGA, Ti and Sc/Ti foils, 3 ns at 5e14):
  one- vs two-sided ratios only, no absolute CE; Do et al., RSI (2020), OSTI 1708986:
  continuum at 11–30 keV only; Babonneau et al., PoP 15, 092702 (2008): pre-exploded
  Ti 9%, Cu 1%, Ge 2.5–3% (2π, bands; pre-exploded, out of class). Phillion & Hailey,
  PRA 34, 4886 (1986) is 0.53 µm (another wavelength; not used). **No V, Sc or Ca He-α
  CE at 351 nm was found**; Mn/Co none either.

**Weights.** Inverse variance in ln CE with a random-effects term: σᵢ² = sᵢ² + τ², sᵢ =
u/ce (the paper's error bar, relative); W01's digitized points (no per-point bar) take
their stated relative 10% ⊕ 20% reproducibility = 0.22; C19 0.20; B13/B14 their bars
(0.20–0.58). τ (the scatter beyond the stated errors) is set so that χ²/dof = 1.
Predictive scatter = √(τ² + s̄²), s̄ = 0.273 the rms sᵢ; σ at a point adds vᵀCv,
C = (AᵀWA)⁻¹.

**Fit, before → after** (same form). *derived*
- 2026-10-24 (37 points, equal weights): a = −5.694/−4.876, b = −3.92, c = 0.244/0.209,
  d = −0.098, e = 0.543; scatter 0.41.
- Same 37 points, weighted: a = −5.736/−4.846, b = −4.16, c = 0.456/0.218, d = −0.145,
  e = 0.853; τ = 0.26.
- **New (68 points, weighted)**: a = −5.597 ± 0.10 (foil), −4.893 ± 0.08 (prepulsed);
  b = −4.27 ± 0.30; c = 0.571 ± 0.11, 0.315 ± 0.10; d = −0.188 ± 0.027; e = 0.961 ± 0.20;
  **τ = 0.243, scatter 0.366 (×1.44)** (was 0.41, ×1.51); coefficient σ's 1.5–3.5× smaller than the weighted 37-point fit's.
  Unweighted residual rms by source: W01 0.46 (mean −0.15), B13 0.44, B14 0.26, C19
  0.28; foil 0.41, prepulsed 0.30.
- Optimum / peak (old → new): foil Fe 1.1e15 / 1.6% → 1.6e15 / 2.2%; Cu 2.1e15 / 0.71%
  → 2.9e15 / 0.97%; Ge 3.8e15 / 0.36% → 4.9e15 / 0.53%; prepulsed Cu 1.8e15 / 1.57% →
  1.5e15 / 1.63%; Ge 3.1e15 / 0.78% → 2.5e15 / 0.78% (closer to Barrios's 1–1.5e15);
  Zr 1.1e16 / 0.20% → 8.1e15 / 0.21%. Peak slope Fe→Ge (foil) hν^-3.3, Cu→Ge
  (prepulsed) hν^-3.6.
- **CE at representative points, old → new** (%, ratio): foil Cr 1e15 3.09 → 4.14
  (×1.34); Fe 1e15 1.62 → 2.05 (×1.27); Fe 3e15 1.48 → 2.01 (×1.35); Cu 1e15 0.67 →
  0.78 (×1.17); Ge 3e15 0.36 → 0.51 (×1.42); prepulsed Cu 1e15 1.52 → 1.59 (×1.05); Zn
  2e15 1.22 → 1.25 (×1.02); Ge 5.1e14 0.56 → 0.49 (×0.86); Ge 1.5e15 0.74 → 0.74
  (×1.01); Zr 4.5e15 0.18 → 0.20 (×1.08); Mo 1e16 0.15 → 0.17 (×1.12); Zr 5.1e14 0.078
  → 0.050 (×0.65, extrapolated). σ 0.43–0.48 → 0.37–0.42 inside the data. Ti foil
  4.6e14 6.7 → 8.9 (×1.34, extrapolated; see below).

**Held-out checks** (independent data, not the points fitted).
- The form fitted to the 2026-10-24 rows alone (weighted, 37 points) predicts C19's 31
  points: mean residual (data − fit) +0.11, rms 0.35; prepulsed −0.04 / 0.26, foil
  +0.27 / 0.42; normalized by the prediction's σ ⊕ 20%, rms 0.71 (`test/
  backlighter_signal.jl` asserts < 1.2). The unweighted 2026-10-24 model does about the
  same (rms 0.44; prepulsed −0.01 / 0.25; Fe +0.48 / 0.63).
- Leave-one-source-out (new form, weighted): rms 0.43 (B13), 0.31 (B14), 0.40 (C19),
  0.57 (W01); normalized rms 0.5–1.15, so σ is about right for a new facility.
- Leave-one-element-out: Cr −0.42 mean (fit too high), Fe +0.41 (too low), Cu −0.22, Zn
  −0.18, Ge +0.15 (rms 0.27–0.57).
- **Ti (P16) fails**: the new fit gives 8.4% at 1.32e15 against 1.58% (He-α share) or
  2.26% (whole band): residual −1.67 in ln CE, **3.3σ**; the old fit gave 5.9% (−1.32
  unweighted fit; −1.45, 2.6σ, for the weighted fit to the old rows). Gekko (1.4%) and
  Shenguang II (~1.75%) Ti foils, in 4–6 keV bands, agree with the OMEGA point.

**The functional form is wrong at the low-Z end** (numbers above): the measured peaks
are flat from Cr to Fe (W01 Cr 3.0% at 0.93e15, unresolved; W01 and C19 Fe 3.0–3.25%)
and lower at Ti (≤ 2.3% for the whole band at 1.3e15), while `b L` with the `X L` drift
makes the peak rise steadily toward low hν (fitted Cr 4.2%, Ti 9.2%). The leave-one-
element-out signs (Cr over, Fe under) show the same tension inside the data. Tested and
rejected: a shared `L²` term (−0.91 ± 0.60 with Ti, Ti still −1.6 held out) and a
foil-only `L²` (−0.31 ± 0.84, Ti −1.59, 2.9σ). The form is kept; Ti stays outside the
fitted lines (flagged `extrapolated`), and the docstring and `examples/backlighter.jl`
say to pass a measured CE. Changing the form needs Ti/V data at stated conditions (an
open decision). *Superseded the same day by "He-α CE: saturating peak (2026-10-26)"
below: the convexity is in the form's parametrization, and a peak plateau below a fitted
hνₛ passes the Ti check without fitting it.*

**Pulse length and spot size** (the "if cheap" part). A `g ln τ` term on the foil class:
g = 0.63 ± 0.21 (0.58 ± 0.22 with a C19 offset of 0.10 ± 0.10; 0.60 ± 0.24 without the
0.6 ns point). It rests on C19's 2 ns Fe (2.8–3.2%) against W01's 1 ns Fe (2.4% at
0.9e15); within C19, the 1 ns Fe at 1.98e15 (3.25%) is above the 2 ns Fe at 1.75e15
(2.84–2.89%). So τ stays confounded with the source; not in the model. A C19 amplitude
offset alone: +0.19 ± 0.11 (two-sided irradiation or calibration; not significant).
Spot size: every source sets the irradiance with the spot (C19 and W01 change spot and
beam count together), so a spot term is degenerate with X here. *E*

## He-α CE: saturating peak (2026-10-26) *engineering*

The form of `kshell_ce_model` changes; the data (`data/kshell_ce.csv`, 68 fitted rows),
the weights and τ rule of the previous section are unchanged. Exploration script not in
the repository; the package fit reproduces its numbers.

**Diagnosis.** The previous form's peak over lines is `a_c + bL + (c_c + eL)²/(4|d|)`:
the `X L` drift makes it convex in `L`, so the peak rises faster than any power law
toward low hν (Cr 4.15%, Ti 9.2%). The shared `L²` term tried before only offsets part of
that curvature. Writing each line's parabola about its own optimum separates the peak
from the optimum's drift:
`ln CE = a_c + P(L) + c_c Y + d Y²`, `Y = X − sL` (optimum `I* ∝ hν^s`, peak `a_c + P(L)
− c_c²/4d`), linear in the rest for fixed `s`.

**Forms compared** (weighted, τ refitted per fit; nonlinear parameters refitted inside
every fold; held-out rms of data − prediction in `ln CE`, unweighted). LOSO = leave one
source out (B13, B14, C19, W01), LOEO = leave one element out (Cr … Mo); Ti = the held-
out P16 point, z with the prediction's σ ⊕ its stated 42%. *derived*

| form | params | −2 ln L | τ | C19 held out | LOSO pooled | LOEO pooled | Ti fit, z |
|---|---|---|---|---|---|---|---|
| previous (a4f1477) | 7 | −82.6 | 0.243 | 0.35 | 0.55 | 0.47 | 8.4%, −2.9 |
| about the optimum, `P = bL` | 7 | −87.4 | 0.226 | 0.34 | 0.44 | 0.44 | 6.1%, −2.4 |
| … `P = bL + gL²` | 8 | −88.4 | 0.226 | 0.37 | 0.60 | 0.45 | 7.1%, −2.6 |
| … `P = b max(L − Lₛ, 0)` (**adopted**) | 8 | −88.5 | 0.226 | 0.34 | 0.45 | 0.45 | 3.0%, −1.1 |
| … same, hinge fixed at Fe (6.7 keV) | 7 | −85.7 | 0.232 | 0.29 | 0.44 | 0.45 | 2.4%, −0.7 |
| … hinge smoothed (softplus, width 0.1 in L) | 8 | −87.3 | 0.230 | 0.33 | 0.45 | 0.45 | 4.9%, −2.0 |
| … hinge + `Y³` (asymmetric in ln I) | 9 | −101.3 | 0.185 | 0.45 | 2.26 | 0.47 | 3.1%, −1.3 |
| … hinge, curvature per class | 9 | −89.3 | 0.227 | 0.37 | 0.62 | 0.51 | 3.0%, −1.2 |
| previous form with `b·max(L − Lₛ, 0)` | 8 | −87.4 | 0.230 | 0.34 | 0.54 | 0.49 | 2.9%, −1.1 |

LOSO per source, previous → adopted: B13 0.44 → 0.43, B14 0.30 → 0.40, C19 0.35 → 0.34,
W01 0.80 → 0.58 (recomputed here; the previous section's LOSO line, 0.43/0.31/0.40/0.57,
was computed differently and its C19 value disagrees with its own C19 held-out 0.35).
LOEO, adopted (mean/rms): Cr see below, Fe +0.38/0.54, Ni +0.03, Cu −0.17/0.35, Zn
−0.29/0.33, Ge +0.10/0.54, Zr +0.12/0.23, Nb +0.44/0.47, Mo +0.30/0.30 (the previous
form: Cr −0.42, Fe +0.41, Nb −0.03, Mo −0.07). The cubic's in-sample gain is overfitting
(W01 left out: rms 3.8).

**What the comparison says.** (1) Most of the cross-validated gain is the
reparametrization at the same parameter count (LOSO 0.55 → 0.44). (2) The plateau adds
one parameter that the data inside Cr–Mo barely need: for `Lₛ` between Cr and Fe it is
exactly a free Cr peak level (Fe and above follow `hν^b`), and Cr sits 0.24 below the
power law, Δ(−2 ln L) = 1.1. Its profile: −2 ln L within 1 of the best for hνₛ 5.85–6.45
keV, within 3 up to 6.7; above 7 keV the fit degrades fast (hνₛ = 7.5: LOSO 0.87). Ti
across that range: 3.2–2.4% (z −1.3 to −0.7); with no plateau (hνₛ below Cr) 5.1–6.1%.
(3) The Ti check therefore selects the plateau, and only a sharp one: a softplus knee
of width 0.1 in `L` fits slightly worse and runs to the grid's edge (4.8 keV), back to
the power law. (4) With Cr left out, hνₛ is not identified (any value up to Fe fits the
rest equally); Cr's held-out residual then ranges from −0.21 (hinge below Cr) to +0.34
(at Fe; what the refined grid lands on). (5) The prepulsed high-Z end (B14 Zr/Nb/Mo)
loses: LOSO B14 0.30 → 0.40, Nb/Mo now underpredicted by 0.3–0.45, because the shared peak
slope `b` is now the Fe–Ge slope.

**Adopted fit** (`kshell_ce_model`, 68 points): `ln CE = a_c + b max(L − Lₛ, 0) + c_c Y
+ d Y²`, a = −3.955 ± 0.20 (foil), −3.175 ± 0.22 (prepulsed); b = −3.33 ± 0.25; c =
0.639 ± 0.10 (foil), 0.262 ± 0.12 (prepulsed); d = −0.203 ± 0.026; s = 2.48 ± 0.43;
Lₛ = −0.497 ± 0.065, **hνₛ = 6.09 keV** (5.7–6.5); τ = 0.23, scatter 0.35 (×1.43).
(s, Lₛ) by maximum likelihood on a refined grid, the other six linear; covariance
(JᵀWJ)⁻¹ over all eight, J the gradient (the hinge's one-sided), which σ uses.
Optimum: foil 4.8e15 (hν/10 keV)^2.48 W/cm² (Cr 1.2e15, Fe 1.8e15, Cu 3.1e15, Ge
5.2e15), prepulsed 1.9e15 × the same (Cu 1.2e15, Ge 2.0e15, Zr 6.4e15). Peak: foil 3.17%
for hν ≤ 6.09 keV, then ∝ hν^-3.33 (Fe 2.30%, Cu 1.09%, Ge 0.55%); prepulsed Cu 1.56%, Ge
0.79%, Zr 0.17%. The C19 held-out check (fit without C19, 37 points): mean +0.10, rms
0.34, normalized 0.69 (was +0.11, 0.35, 0.71).

**CE at representative points** (%, 2026-10-24 main fit → a4f1477 weighted refit →
adopted; ×/÷ the adopted 1σ factor). *derived*

| point | 2026-10-24 | a4f1477 | adopted |
|---|---|---|---|
| foil Ti 1e15 | 6.23 | 8.91 | 3.12 ×/÷ 1.49 |
| foil Ti 1.32e15 (held-out 1.58 ± 0.66) | 5.91 | 8.43 | 2.97 ×/÷ 1.48 |
| foil Cr 1e15 (W01 3.02 at 0.93e15) | 3.09 | 4.14 | 3.15 ×/÷ 1.50 |
| foil Fe 1e15 | 1.62 | 2.05 | 2.15 ×/÷ 1.44 |
| foil Fe 3e15 | 1.48 | 2.01 | 2.17 ×/÷ 1.44 |
| foil Cu 1e15 | 0.67 | 0.78 | 0.84 ×/÷ 1.44 |
| foil Ge 1e15 | 0.30 | 0.33 | 0.32 ×/÷ 1.45 |
| foil Ge 3e15 | 0.36 | 0.51 | 0.52 ×/÷ 1.44 |
| prepulsed Cu 1e15 | 1.52 | 1.59 | 1.55 ×/÷ 1.44 |
| prepulsed Zn 2e15 | 1.22 | 1.25 | 1.21 ×/÷ 1.44 |
| prepulsed Ge 5.1e14 | 0.56 | 0.49 | 0.54 ×/÷ 1.47 |
| prepulsed Ge 1.5e15 | 0.74 | 0.74 | 0.78 ×/÷ 1.43 |
| prepulsed Zr 4.5e15 | 0.18 | 0.20 | 0.17 ×/÷ 1.47 |
| prepulsed Mo 1e16 | 0.15 | 0.17 | 0.12 ×/÷ 1.48 |
| prepulsed Zr 5.1e14 (extrapolated) | 0.078 | 0.050 | 0.047 ×/÷ 1.66 |

**Cr–Fe flatness.** The fit does not make the Cr and Fe peaks equal (3.17% vs 2.30%).
The "flat at ~3%" reading pairs W01's Cr (3.02% at 0.93e15) with C19's two-sided Fe
(2.84–3.25%) and W01's single Fe point at 3.5e15 (3.0%); within W01 at ~0.9e15, Cr/Fe =
3.02/2.38 = 1.27, and the fit gives 3.13/2.14 = 1.46 there. Fe stays underpredicted
(in-sample mean +0.11, LOEO +0.38) in every form tried: the sources disagree by about the
two-sided factor (previous section). Forcing the plateau to start at Fe (7 parameters,
table) makes Cr and Fe equal at 2.42%, undershoots W01's Cr by 0.25 in `ln CE`, and costs
2.8 in −2 ln L.

**Mechanism** (not modeled; plausible, not tested here): the peak CE cannot grow without
bound toward low Z, because the absorbed fraction and the share of it radiated are
bounded, and at the lower temperatures that suffice for low-Z He-like ions more of the
energy goes into lower ionization stages and continuum. The Ti point (and the Gekko XII
and Shenguang II Ti foils, 1.4–1.8% in 4–6 keV) are below the plateau, which would mean
a turnover rather than a plateau; one point cannot separate them. *confidence: medium
inside the data (unchanged), low below Cr.*

**Tests** (`test/backlighter_signal.jl`): χ²/(68 − 8) = 1, weighted residual sums zero
per class, the C19 held-out check (normalized rms < 1.2), the held-out Ti within 1.5σ,
hνₛ between Cr and Fe, the Ti and Cr peaks equal (the plateau) and above Fe's, the
gradient behind σ against central differences, the optimum a maximum, the hν slope in
2.5–4.5.

## Dante channel responses and the unfold (2026-10-24) *derived*

ROADMAP 13, what remained of the synthetic Dante: the 18 Dante 1 channels built from
their hardware in place of Gaussians, the Dante unfold, and XRD lines of sight
(`src/xray_optics.jl`, `src/dante_response.jl`, `examples/dante_channels.jl`,
`test/dante_response.jl`). Read: Rubery et al., RSI (2023), OSTI 1975922 (Fig. 4
hardware table, Figs. 3 and 8, §II.C–D); Harris et al., "Spectral sensor error analysis
for measuring x-ray radiation drive using the DANTE diagnostic", RSI 92, 033505 (2021),
OSTI 1820019 (eq. 1, the standard unfold, §II); Barnak et al., "Soft x-ray spectrum
unfold of K-edge filtered x-ray diode arrays using cubic splines", RSI 91, 073102
(2020), OSTI 1642075; Kemp et al., "A compact filtered x-ray diode array spectrometer
for the NIF: SENTINEL", RSI (2020), OSTI 1736345; Rubery et al., "A 2–4 keV multilayer
mirrored channel for the NIF Dante system", RSI 93, 113502 (2022), OSTI 1959490 (not
needed beyond the 2023 review's Fig. 8). Confidence: **H** primary text/table, **H-img**
digitized raster figure, **L** our construction or assumption.

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| Dante 1 hardware, 18 channels | energies, filters (µm), 8.8 % W grids, mirrors (ch 2 2–4 keV multilayer, ch 3 3.5° SiO₂, ch 4 2.5° SiO₂), XRD cathodes (Al/Cr/Ni), solid angles (6.51e-7 sr; ch 2 1.90e-7, 3 5.74e-7, 4 5.56e-7, 16 and 18 1.77e-6), uncertainties 8–16 % | Rubery 2023 Fig. 4 (a table image, read at 300 dpi) | H |
| Channel 2 filter | "47 µm CH" in Fig. 4; 47 µm polyimide in §II.D and Fig. 8 — polyimide used | Rubery 2023 | H |
| Grids | 8.8 % open tungsten, behind the filters, to keep XRDs below the 500 V Child–Langmuir limit | Rubery 2023 §II.A | H |
| XRD-31 | reflection-mode planar vacuum diodes, Ni/Al/Cr/Au cathodes, 5 kV bias, 140 ps | Rubery 2023 §II.A | H |
| Cr XRD sensitivity | 0.07–1 keV, C K (277), O K (525), Cr L (573 eV) structure; axis labelled V/GW | Rubery 2023 Fig. 3b | H-img |
| Al XRD sensitivity | 1–4.5 keV, MV/GW, Al K edge | Rubery 2023 Fig. 8b | H-img |
| SiO₂ mirror "3.5°" | measured reflectivity 0–1 keV (0.84 → 0.06) | Rubery 2023 Fig. 3a | H-img |
| Ti 0.75 µm filter | measured "effective transmission" 0–1 keV | Rubery 2023 Fig. 3a | H-img |
| Cr/Sc multilayer, 1.9° | reflectivity 1.2–4.5 keV; system response (mirror × polyimide × Al XRD) in kV/GW | Rubery 2023 Fig. 8 | H-img |
| Channel 2 calibration | 0.97 TW/sr in 2–4 keV per volt at the XRD, 8.7 % | Rubery 2023 §II.D | H |
| Signal equation | V = Ω ∫ S(E) R(E) dE, S in GW/sr/eV, R in V/GW | Harris 2021 eq. 1 | H |
| Standard unfold | Planckian start from the channels < 1.5 keV (Rubery) / channels 1–9 (Harris), averaged per-channel brightness temperatures, plus an M-band Gaussian (1 GW/eV/sr, 2.5 keV, width 600 eV); each response replaced by a fitted Gaussian whose amplitude is adjusted channel by channel from the highest energy down, ~10 iterations to 5 % | Harris 2021 §II; Rubery 2023 §II.C | H |
| Sentinel | DIM-mounted, 16 MiniXRDs (6 Al, 5 Ni, 5 Cr) behind Dante filters, polar DIM 000-000 or equatorial 090-078 / 090-315, ~7 m | Kemp 2020 §II | H |
| Optical constants | Henke f₁, f₂, 10 eV–30 keV, 19 elements (`data/henke/`, CXRO, fetched 2026-10-24) | Henke, Gullikson & Davis, ADNDT 54, 181 (1993) + CXRO updates | H |

Digitization (`data/rubery2023_dante.csv`): the figures are raster images in the
accepted manuscript; the plot frames were located from the dark frame lines and the
axis ends taken at the frame (where the end tick labels sit), each curve masked by colour
and its median row read per pixel column. Overlaying the read points on the figures
agreed to about a pixel (~0.1 % of full scale); the narrow edge spikes (XANES) are
under-read by the median. Units: Fig. 3b's axis says V/GW; read as MV/GW like Fig. 8b's
Al XRD — with V/GW the quantum efficiency would be ~10⁻⁸, and with MV/GW channel 2
reproduces its published calibration (below).

### Model

- **Filters**: Henke photoabsorption (`henke_mass_attenuation`; `mass_attenuation(...;
  table = :henke)`, with `:auto` keeping NIST where it covered the backlighter filters
  before, so `backlighter_signal` is unchanged). Henke vs NIST at 3 keV: Be, Al, Ti, Fe,
  Cu, polyethylene within 3 % (NIST includes scattering). Compound densities CH 1.05,
  parylene (N) 1.11, Saran 1.70 are typical values (**L**).
- **Grids**: an energy-independent 0.088 (wires opaque) (**L**: the wires' few-µm W
  transmits nothing below ~10 keV).
- **SiO₂ mirrors**: Fresnel (s and p averaged) from Henke f₁, f₂ with a Névot–Croce
  roughness. At the nominal 3.5° the model's cutoff sits at ~0.5 keV, the measured one
  at ~0.9 keV (rms 0.24 in reflectivity); fitted, the measured curve is SiO₂ at 2.04°
  with 4.25 nm roughness (rms 0.037). The published angle is therefore not the grazing
  angle of a bulk SiO₂ surface; the 2.5° mirror takes the same 0.584 scale (1.46°)
  (**L**).
- **Multilayer (ch 2)**: Fig. 8a's reflectivity as read (held at its end values
  outside 1.2–4.5 keV).
- **XRDs**: Henke's total-yield form S = K (E/1 keV)^γ x/(1 + x), x = L ρ μ_s(E),
  with μ_s the mass attenuation of the bulk metal plus a surface layer (mass fractions
  w_C, w_O, from the C and O edges in Fig. 3b). Fitted in log to Fig. 3b (Cr: K 15.84
  MV/GW, L 33.6 nm, w_C 0.274, w_O 0.250, γ −0.244; rms 12 %) and to Fig. 8b (Al: K 35.12
  MV/GW, γ 0.102 on Cr's surface layer; rms 6.7 %; the K-edge region 1.54–1.75 keV
  excluded from the fit). Ni (no published curve) takes Cr's parameters on Ni (**L**);
  Cr above 1 keV and Al below 1 keV are extrapolations of the fitted form (**L**). As
  quantum efficiencies: Al 0.016 electrons/photon at 1 keV, Cr 0.03 at 70 eV.
- **Checks against independent published numbers**: channel 2 (multilayer × 47 µm
  polyimide × Al XRD, Ω 1.9e-7 sr) gives **0.968 TW/sr per volt** for power flat over
  2–4 keV against Rubery's **0.97**; its response is 0.942 × Fig. 8b's over 2–4 keV
  (rms 6 %). Not fitted. The Ti filter of Fig. 3a transmits 0.45 × Henke's 0.75 µm at
  0.45 keV and 0.83 × at 1 keV: measured filters absorb more than their nominal
  thickness at low energy (oxide, contamination, support), a caveat for every thin-filter
  channel here, all of which use nominal thicknesses.

Channel bands (half maximum, eV) against Rubery's nominal energies: 1 (1170) 1084–1305;
2 (3060) 2307–3625; 3 (250) 226–284; 4 (470) 444–513; 5 (710) 656–779; 6 (870) 867–934;
7/15 (940) 852–1024; 8 = 0.088 × 1; 9 (1430) 1314–1561; 10/14 (1710/1700) 1558–1840;
11 (3430) 1821–4467; 12 (2580) 2351–2825; 13 (3000) 2695–3357; 16 (14370) 8579–19247;
17 (690) 624–779; 18 (3660) 1771–4200. Each band's peak lies within 4–13 % above the
nominal energy (the absorption edge that closes it), except the broad 11, 16 and 18.
A Planckian at 300 eV through the LEH (π 1.55², 37°) gives 190 V on channel 1 and 14 V
on channel 3 (at the XRD, before attenuators), the scale of Harris's N190721 traces
(90–230 V on the low channels).

### Unfold (`dante_unfold`, `dante(view, ::SpectralSolution)`)

- `:traditional` — Harris's algorithm: Planckian start at the mean of the per-channel
  brightness temperatures of the channels below 1.5 keV, plus the M-band Gaussian; each
  response replaced by a unit Gaussian (centre between its half-maximum points, FWHM
  their distance: Harris fits one but does not say how, **L**); amplitudes adjusted
  channel by channel from the top. With Dante 1's duplicated hardware (1/8 differ only
  by the grid, 7/15 and 10/14 are identical) the Gaussians coincide and the sweep
  diverges; the same equations are then solved directly (least squares in the relative
  voltages, minimum-norm amplitudes).
- `:smooth` — positive: the start Planckian (no M-band) times a piecewise-linear factor
  on 48 log-spaced knots, nonnegative least squares on the relative voltages with a
  second-difference penalty (λ 1e-3). A Planckian is its own solution (factor 1), so
  it is exact; it is our construction (**L**), closest in spirit to Barnak 2020's
  spline unfold (which uses one knot per channel edge and no prior).
- Tr = (π I / (σ A cos θ))^¼ over the same area and angle as the synthetic Dante's
  `Tr`, I integrated over 50 eV – 20 keV.

**Closed-form check (Planckians)**: voltages computed by quadrature of the exact
Planckian through the responses, then unfolded. The exact answer is
T · F(T)^¼ with F the Planck fraction in 50 eV – 20 keV (the unfold's range):

| T (eV) | exact | traditional | smooth |
|---|---|---|---|
| 100 | 99.867 | 105.40 (misfit 95: fails) | 99.867 |
| 150 | 149.937 | 149.98 (misfit 1.3) | 149.937 |
| 200 | 199.964 | 199.84 (misfit 0.07) | 199.964 |
| 250 | 249.976 | 250.047 | 249.976 |
| 300 | 299.983 | 300.013 | 299.983 |
| 350 | 349.988 | 350.007 | 349.988 |

The smooth unfold and the traditional one without the M-band start agree with the
closed form to 1e-6. With Harris's absolute M-band start (1 GW/sr/eV, sized for a NIF
peak drive) the traditional unfold is within 0.07 eV from 250 eV up and fails below
~200 eV, where that start is larger than the whole signal of the M-band channels and the
Gaussians cannot cancel it within the channels' bands.

### N221204 (`examples/n221204.jl energetics=new unfold=on`)

Standard mesh (0.28 / 0.16), full model, published energetics; the spectrum Dante 1
sees from `spectral_solve` (42 groups to 20 keV), each element's source split as a
Planckian at its own brightness temperature ("Planckian"), and with the intensity-form
laser M-band (`MBandEmission()`, "+ M-band"). Within each group the voltages take the
spectrum as a Planckian at Tr_open (**L**). Tr over the initial LEH area (eV):

| t (ns) | spectrum | input Tr | traditional | smooth | color T, hardware / Gaussian channels | > 1.8 keV |
|---|---|---|---|---|---|---|
| 5.0 | Planckian | 248.8 | 249.7 (+0.90) | 248.8 (−0.02) | 272.2 / 271.9 | 0.098 |
| 5.0 | + M-band | 248.8 | 249.8 (+1.02) | 248.8 (−0.01) | 276.9 / 271.8 | 0.139 |
| 6.0 | Planckian | 271.4 | 272.3 (+0.90) | 271.3 (−0.02) | 301.1 / 300.5 | 0.144 |
| 6.0 | + M-band | 271.4 | 272.4 (+1.02) | 271.4 (+0.00) | 311.8 / 300.0 | 0.199 |
| 7.0 | Planckian | 283.8 | 285.0 (+1.13) | 283.8 (+0.02) | 328.3 / 326.8 | 0.192 |
| 7.0 | + M-band | 283.8 | 285.1 (+1.27) | 283.9 (+0.04) | 352.0 / 325.4 | 0.284 |
| 8.0 | Planckian | 295.0 | 296.2 (+1.18) | 295.0 (+0.04) | 346.6 / 345.2 | 0.223 |
| 8.0 | + M-band | 295.0 | 296.3 (+1.30) | 295.0 (+0.05) | 377.6 / 344.2 | 0.331 |
| 8.5 | Planckian | 283.4 | 284.2 (+0.82) | 283.3 (−0.02) | 313.2 / 312.1 | 0.165 |
| 8.5 | + M-band | 283.4 | 284.2 (+0.87) | 283.3 (−0.02) | 319.3 / 311.4 | 0.195 |

(coarse mesh 0.8 / 0.4: the same to 0.4 eV; every unfold matched all 18 voltages to
< 1e-7.) The gray history is unchanged by this work: Dante 1 295.0 eV at 8 ns over the
initial area, 318.6 over the open area, against 310.1 measured.

- **The unfold reads the input's Tr**: the traditional unfold 0.8–1.3 eV high (0.3–0.4 %
  in Tr, ~1.5 % in flux; its Gaussians spread the high channels' power into the
  thermal band), the smooth one within 0.05 eV. Both well inside Dante's quoted 1σ of
  1.5 % in Tr (Rubery 2023 §II.C). So comparing the synthetic Dante's `Tr` with the
  measured Dante Tr, as the N221204 tables do, was sound to ~1 eV; the brackets of
  initial vs open area (24 eV) and the late-time energetics remain the open terms.
- **Before/after for the channels**: with Gaussians there was no unfold, only the
  color temperature of the channels below 1.5 keV. That color temperature agrees
  between Gaussian and hardware channels for the Planckian spectra (≤ 1.5 eV), but
  with the laser M-band the hardware channels read 5–33 eV hotter: their high-energy
  leakage (nominal-energy bands with tails past the absorption edges) picks up the
  M-band, which Gaussians cannot. The color temperature is a start value, not a
  result, so Tr is unaffected; it is a reason not to quote a < 1.5 keV color
  temperature from Dante as the wall temperature.
- **Model finding — hot laser spots**: even without the laser M-band, 10–22 % of what
  Dante 1 sees lies above 1.8 keV and the color temperature is 23–52 eV above Tr: the
  per-element Planckian split puts the laser spots at their brightness temperature
  (B/σ)^¼, far above the wall's, and Dante 1 sees outer spots through the LEH. With the
  intensity-form M-band on top the share reaches 0.33 at 8 ns, against Dewald's
  measured 12–15 % (with the per-cone `nif_mband()` and α_M 0.5 the Dewald setup gave
  12 %). Not acted on here: a spot's thermal emission should be split at a plasma
  temperature model rather than its brightness temperature. (Resolved in "Laser-spot
  spectrum (2026-10-24)" below: spots now thermal at the local Tr plus the M-band; the
  0.33 was mostly the uncalibrated M-band pairing, 0.115–0.138 calibrated.)

### XRD lines of sight

`sentinel(:polar)` (Kemp 2020: polar DIM 000-000, θ 0°) views the hohlraum down its
axis through the upper LEH: on the test hohlraum (LEHs r 1.6 mm, capsule r 0.9 mm) it
sees the capsule and, around it, the lower LEH, so a blackbody cavity at 250 eV reads
the capsule's (0.9/1.6)² share (186 eV): the polar view measures the capsule, not the
wall. The equatorial DIMs (090-078, 090-315) see a hohlraum only through a window
(`through` required). Sentinel's own filter set is not tabulated (Dante 1's channels
stand in, **L**). The OMEGA Dante (port H16I, θ 142.62°, φ 126°; Rubery 2023, the Omega section) depends
on the hohlraum axis chosen per shot and is left to `DanteView`.

### What remains

- Ni XRD measured shape (Cr's parameters stand in); Al below 1 keV and Cr above 1 keV
  are extrapolations; measured thin-filter transmissions (Fig. 3a's Ti is 0.45–0.83 ×
  nominal) — real filters' oxide and contamination layers would move the low channels.
- The SiO₂ mirror angle scale (published 3.5° behaves as 2.04° bulk SiO₂) applied to
  the 2.5° mirror is an assumption.
- Within-group spectral shape: voltages assume a Planckian at Tr_open within each group;
  finer groups near the channel edges would remove it.
- Laser-spot spectra (above; done, "Laser-spot spectrum" below) and the M-band split
  per cone in time-dependent runs.
- Dante 2's hardware (11 mirrored channels; not tabulated in Rubery 2023) and the
  Omega Dante configuration (Rubery 2023 Fig. 13).

## Laser-spot spectrum (2026-10-24) *derived*

ROADMAP 13: the finding above ("hot laser spots"). `spectral_solve` split each
element's whole source as a Planckian at the brightness temperature (B/σ)^¼ of its total
radiosity. On a laser spot B = ηI + αH with ηI ≈ 1.1 × 10¹⁵ W/cm² at the outer cones'
1.5 × 10¹⁵, so a ~300 eV cavity put its spots at ~365 eV, whose Planck tail alone holds
25 % of the power above 1.8 keV (a 290 eV Planckian 12.4 %, 300 eV 14.0 %), and the
laser M-band was then carved out of that already hard spectrum. Read (cite only; in the
user's paper folder): Dewald et al., "X-ray conversion efficiency of high-Z hohlraum wall
materials for indirect drive ignition", PoP 15, 072706 (2008), LLNL-JRNL-402374 (OSTI
944359; abstract only, the full text is behind the LLNL library); Chen et al., "On the
deficiency of x-ray flux predictions of indirect-drive ICF experiments at the NIF",
PRE (2024), LLNL-JRNL-862277 (OSTI 2429360); Sacks et al., "A simple non-Planckian
radiation source for ICF and HED simulations", LA-UR-23-31376 rev. 1 (2024, OSTI
2377325); Dewald et al. PoP 27 (2020) (above). Confidence: **H** primary text, **M**
abstract or secondhand, **L** our construction.

### What the literature fixes

| Statement | Source | Conf. |
|---|---|---|
| Au spheres, direct drive, no cavity: x-ray conversion 95 % at 10¹⁴, 80 % at 10¹⁵ W/cm²; M-band (2–5 keV) ~1 % of the x-ray flux at 10¹⁴, 27 % at 10¹⁵ (LASNEX with XSN under-predicts it) | Dewald 2008 abstract | M |
| The > 1.8 keV emission in NIF hohlraums is "predominantly produced by hot, low-density NLTE gold bubble plasma" (the laser-heated corona); the M-band/total ratio is higher through the LEH, which sees the outer-beam bubble | Chen 2024 §III | H |
| A hohlraum drive spectrum = a Planckian at Tr plus a non-Planckian power law above 1.85 keV, normalized to the total; used in place of integrated hohlraum simulations and compared with them and with Dante-driven experiments | Sacks 2024 §2–3 | H |
| The recirculating M-band at the wall is the walls' Planckian thermal tail; the laser-produced part comes from the outer spots | Dewald 2020 §V.2 | H |

So the converted x-rays of a spot have two parts: the non-Planckian M-band, set by the
laser intensity (Dewald 2008's 27 % is the laser-produced hard fraction itself; an
extra Planck tail at an inflated temperature counts it twice), and a soft thermal part.
The soft part's temperature is not published for spots in a cavity (Ze et al., PoP 4,
1997, OSTI 467314, compare spot and disk x-ray yields, not spectra). The brightness
temperature is a bookkeeping number of the gray radiosity, not the temperature of an
emitting surface: the converted energy is radiated by a conversion layer bathed in, and
thermalized against, the cavity field. We take that layer's thermal emission at the
local radiation temperature (**L**: the Sacks form applied per element).

### Model (`spectral_solve(...; spot)`)

- `spot = :radiation` (new default): the converted x-rays `ηP/A` (`_xray_source`, the
  conversion model's own number) less the M-band at the element's radiation temperature
  `(H/σ)^¼`; the rest of the source (a re-emitting wall's stored energy) and, with
  `thermalize`, the thermal re-emission at the wall temperature `((B − ηP/A)/σ)^¼` (on a
  wall without laser light that is (B/σ)^¼, as before); the M-band flat over its band.
  Both temperatures follow the group totals when hard groups lose energy.
- `:wall`: the converted x-rays at the wall temperature too (a lower bracket).
- `:brightness`: the earlier rule, bit for bit (the Dewald scan and the N221204 rows
  below reproduce the earlier tables to the printed digit).
- Group sums are untouched: Planck fractions sum to one per element, so Σ_g S_g = S
  and every gray total of `test/spectral.jl` holds for all three rules (B, H to 1e-11,
  absorbed power to 1e-11). New tests: the source rows against the closed form of each
  rule (1e-12), spots softer than under the brightness rule, the rules identical on
  unlit elements, a stored-energy emission split at the wall temperature.
- `SpectralSolution` gains `spot_temperature`, `converted`, `spot`.

### Dewald's setup (`examples/mband_dewald2020.jl spot=…`, h 0.3/0.15, per-cone `nif_mband()`)

Cone fraction 1/3 (1/4 in brackets); Dante 1 Tr 291.9 (290.6) eV, Dewald 290.

| rule | α_M | F_M | capsule M-band P2/P0 | boundary P2/P0 | Dante-1 > 1.8 keV | Dante-1 Tc (eV) |
|---|---|---|---|---|---|---|
| brightness (before) | 0.5 | 2.15 (2.05) | +0.039 (+0.145) | −0.012 (+0.200) | 0.119 (0.125) | 294.8 (296.2) |
| brightness | gray | 7.63 | −0.036 | −0.140 | 0.300 | 336.7 |
| radiation | 0.5 | 1.82 (1.73) | +0.051 (+0.145) | +0.009 (+0.197) | 0.107 (0.112) | 278.8 (279.1) |
| **radiation** | **0.55** | **2.06 (1.95)** | **+0.041 (+0.127)** | −0.008 (+0.165) | **0.114 (0.120)** | 279.8 (280.2) |
| radiation | 0.6 | 2.33 (2.22) | +0.030 (+0.109) | −0.026 (+0.133) | 0.124 (0.130) | 281.0 (281.5) |
| radiation | gray | 6.86 | −0.034 | −0.134 | 0.273 | 309.5 |
| wall | 0.55 | 1.98 (1.88) | +0.049 (+0.134) | +0.008 (+0.178) | 0.111 (0.117) | 275.3 (275.7) |

- **The α_M calibration was compensating, slightly**: the brightness rule's spot tails
  counted as recirculating M-band, so F_M ≈ 2 took α_M 0.5; with thermal spots it takes
  **α_M ≈ 0.55** (thermalized). Recalibrated to F_M (Dewald's number, not to Dante):
  Dante 1 then sees 11.4–12.0 % above 1.8 keV against Dewald's measured 12 %
  (11.9–12.5 % before); the capsule M-band P2/P0 moves by +0.002 / −0.018.
- The colour temperature (channels < 1.5 keV) drops from 3–6 eV above Dante's Tr to
  ~11 eV below it: the walls emit at α^¼ Tr ≈ 0.95 Tr, the spots at Tr, and 11 % of
  the flux is M-band. Dewald does not publish a colour temperature, so this is no
  check either way.

### N221204 (`examples/n221204.jl energetics=new unfold=on`, standard mesh 0.28/0.16)

Dante 1 over the initial LEH area (Tr; open-area Tr_open in brackets), unchanged (gray):
5 ns 248.8 (262.3), 6 ns 271.4 (289.9), 7 ns 283.8 (305.9), 8 ns 295.0 (318.6), 8.5 ns
283.4 (306.1) eV. "> 1.8 keV" is the share of what Dante 1 sees; Tc the
hardware-channel colour temperature; "cones" the per-cone M-band (`nif_mband()`, each
element's conversion weighted by the cones' shares of its pre-CBET deposition; laser →
M-band 0.13–0.19 of the absorbed laser, the intensity form 0.04–0.13) with the
calibrated α_M, thermalized.

| t (ns) | brightness | thermal | brightness + `MBandEmission()` | thermal + `MBandEmission()` | brightness, cones, α_M 0.5 | **thermal, cones, α_M 0.55** |
|---|---|---|---|---|---|---|
| > 1.8 keV, 5.0 | 0.098 | 0.069 | 0.139 | 0.111 | 0.120 | **0.115** |
| 6.0 | 0.144 | 0.108 | 0.199 | 0.166 | 0.121 | **0.116** |
| 7.0 | 0.192 | 0.137 | 0.284 | 0.239 | 0.129 | **0.122** |
| 8.0 | 0.223 | 0.162 | 0.331 | 0.280 | 0.126 | **0.118** |
| 8.5 | 0.165 | 0.136 | 0.195 | 0.168 | 0.143 | **0.138** |
| Tc − Tr_open, 5.0 | +9.9 | −10.3 | +14.6 | −7.5 | +10.9 | **−8.3** |
| 6.0 | +11.2 | −9.9 | +21.9 | −1.1 | +5.6 | **−9.2** |
| 7.0 | +22.4 | −7.4 | +46.1 | +14.3 | +5.5 | **−10.2** |
| 8.0 | +28.0 | −5.7 | +59.0 | +23.9 | +0.8 | **−14.0** |
| 8.5 | +7.1 | −8.6 | +13.2 | −2.6 | +2.0 | **−7.8** |

(Tc − Tr over the initial area, the comparison of the finding: +23 to +52 eV before,
+3 to +18 thermal, without M-band.) Coarse mesh (0.8/0.4): the rule-only and calibrated
rows the same to 0.007 in share and 0.5 eV in Tc (thermal) / 2 eV (brightness); the
intensity-form rows 0.015–0.02 lower (its element-averaged intensity is mesh
dependent). Unfolds at 8 ns (input 295.0 eV): traditional 296.2 (brightness) → 295.5
(thermal), 295.4 (calibrated); smooth within 0.25 eV in every case; all 18 voltages
matched to < 5e-4.

- **The 0.33 at 8 ns was mostly the M-band pairing, not the spot split**: the example
  used the intensity-form M-band at the gray albedo in the M-band groups (α_M = α,
  F_M ≈ 7 on Dewald's setup), which no calibration supports. With the Dewald-calibrated
  pairing (per-cone M-band, α_M thermalized) both rules give 12–14 %: brightness
  0.120–0.143, thermal **0.115–0.138**, against Dewald's 12 % (Fig. 1b, at the 290 eV
  peak) to 15 % (Fig. 3). His shot is a different hohlraum, so agreement to ~2 points is
  all this supports; not fitted to it.
- **What the spot split changes**: without any M-band, Dante's > 1.8 keV share falls
  from 0.10–0.22 to 0.07–0.16 (the cavity's own Planck tail: 0.07 at 250 eV, 0.14 at
  300); the colour temperature goes from 7–28 eV above Tr_open to 6–10 eV below it; the
  traditional unfold's bias falls from +0.8–1.3 to +0.15–0.5 eV (its Gaussians have
  less hard power to spread). With the calibrated M-band the colour temperature sits
  8–14 eV below Tr_open (brightness rule: 1–11 above).
- Neither is a Dante check: published N221204 colour temperatures are not in hand, and
  the M-band share was compared, not fitted.

### What remains

- A measured spot spectrum (soft part) in a cavity: Dante or a spectrometer at a view
  dominated by spots, or a hohlraum simulation's spot emission, to replace the local-Tr
  assumption (**L**); Dewald 2008's full text (spectra vs intensity) was not reachable.
- The per-cone M-band in `drive_history` records: the example's pre-CBET cone weighting
  is a proxy; bubble-absorbed inner light radiating from the outer bubble (Dewald's
  clipping) is not split per cone.

## Zylstra 2022 LEH pair, time-dependent (2026-10-25) *derived*

ROADMAP 16: the LEH-size pair run with each shot's measured pulse, the remaining
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
cone fraction, with its rise and droop moved onto each shot's own by a piecewise-linear
time map through `pulse_rise` and the half-power fall), the pointing (the mock HYBRID-E
pointing; N191118's outer beams sat farther from the midplane), the capsule (the
N210808 design, 1.1285 mm, for both), the CBET (the joint fit `hybride_cbet`, its tables
shifted by the rise difference, −0.03 and +0.03 ns), and the capsule schedule (its end
moved to each shot's bang time). Backscatter is unpublished for the pair and is left
out; `energetics=new` otherwise (DU, `nif_conversion`, `CapsuleAlbedo` α_c 0.22, 11 %
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
- **The histories close most of the gap.** On the instrument's own reading (Dante) the
  model puts the 3.64 mm shot 8–10 eV above the 3.1 mm one at 5–6.5 ns, against 4–8 eV
  measured: 2–4 eV too much, ~1 % in T, ~3–4 % in power (*L*), against ~6 % before.
  The sign is the same as the static result's (the model's larger LEH loses slightly too
  little). The difference does not depend on the energetics (old and new agree within
  1 eV): it is geometry, not wall or conversion physics.
- **Dante's view widens the model's difference rather than closing it**: the capsule
  sees +5 eV, Dante +9. Through the larger hole Dante 1 sees deeper into the hohlraum,
  so the 3.64 mm shot reads hotter relative to its capsule. Comparing the model's capsule
  Tr with the measured Dante would hide the residual.
- **Computed lips widen it further** (+11–13 eV on the initial area): `LipClosure`'s
  absolute depth (0.24 mm) closes the 3.1 mm hole by a larger fraction, so Dante's
  initial-area reading darkens more on N210601. The prescribed relative closure (0.857
  r₀) does better here; Young's 0.26 mm on 3.64 mm and Schneider's 0.25 on 3.1 mm cannot
  separate the two forms (see "Computed LEH closure").
- **When the laser falls the difference flips sign too hard** (7.5 ns: −14 against
  −4.4): N191118's model Tr drops with the power while its Dante holds for ~0.3 ns more.
  A 0.15 ns lag halves that (−7.4). This is the late-time drive problem below.
- Absolute level, rise and peak: the measured Dante lies between the initial-area and
  open-area readings on both shots (N191118: +0.8 / +19 eV mean, initial / open,
  5–7 ns; N210601: −4.1 / +15 eV, 5–7.5 ns, coarse), as on N221204.

### 2. After the laser falls, and the drive at bang

On all three shots the model's drive falls far faster than Dante's after the power
drops. Model minus measured (initial / open area, `energetics=new`, standard mesh):
N191118 −41 / −21 eV at 8.0 ns, −73 / −58 at 8.5, −71 / −58 at 9.0; N210601 −27 / −5 at
8.0, −44 / −24 at 8.5, −68 / −51 at 9.0; N221204 −28 / −5 at 8.5, −41 / −20 at 9.0. At
bang time Dante reads 226.0 eV (N191118) and 253.8 eV (N210601), the model 154.7 / 167.1
and 185.5 / 200.4. The ordering Zylstra points to ("by bang time Trad is significantly
higher on N210601") is reproduced: model +31 eV (initial area), data +28. The level is
not.

Size of the gap (*L*): the LEH x-ray energy that Dante's Tr implies over the fall
(σT⁴ over both initial LEH areas, Lambertian) exceeds the model's by 125 kJ (N191118,
7.25–9.37 ns), 89 kJ (N210601, 7.75–9.48) and 71 kJ (N221204, 8.0–9.4). The LEHs carry
25–31 % of the model's losses, so holding the hohlraum at Dante's temperature would take
~250–450 kJ more x-ray source after the laser falls, 15–25 % of the laser energy. That is
several times the He gas's 30–60 kJ ("Computed LEH closure") and more than the wall
already gives back in the model (`wall_release`: 137–272 kJ). A 0.15 ns lag of the source
(~60 kJ stored) recovers only 10–15 eV of it (N191118 at 8.5 ns, coarse: −59 instead of −73). So the
tail is not a reservoir of the published size. What remains is either physics a gray,
quasi-static radiosity lacks over a 1–2 ns coast (a hot, optically thick LEH plasma or
gold blow-off in Dante's view; the wall's re-emission when the drive is fully cut, which
`release` = 0.77 was calibrated for only up to 30 % drops), or a Dante time base or
analysis that differs after laser-off. The pair, with N191118's laser off 1 ns before
bang, is the cleanest data for it. Open.

### 3. Bang time and the hot spot (ROADMAP 11 metrics)

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

- **Velocity**: zero to round-off. The mock pointing and the hohlraum are mirror
  symmetric and there are no windows, so the model has no mode 1. N221204's 36 km/s
  comes from sources that are not in it: per-quad power imbalance (DN24 Fig. 1 left,
  not attributable per quad), capsule thickness mode 1 (0.13 µm), the three equatorial
  diagnostic windows and the waist manufacturing gaps (H24), none of them published in a
  usable form. This check cannot pass with public inputs. It says the drive's own mode 1
  is not the cause only if those inputs are known.
- **Hot-spot P2**: with no shipped hot-spot kernel, the comparison is with HYDRA's flux
  history, whose design gave a near-round hot spot (0.3 µm measured). Over the whole
  window the published energetics match Fig. 10's mean to 0.002 (old: 0.013 more
  pole-hot). Over the peak window that the −9.2 µm per % transfer was derived on
  (Callahan 2018's line on another platform, good to a factor of ~2), `energetics=new`
  gives a −6.5 µm (oblate) shift and the old energetics +22 µm. Neither is zero within
  the transfer's uncertainty band; the published energetics are 3× closer.
- **Bang time as the end of the window**: extending the window from 8.9 ns to bang moves
  the mean flux P2 by 0.0014 (new): the coast's flux is low and its P2 small, so the
  window end hardly matters for P2. It matters for the drive. At bang the model's Tr is
  70 eV below Dante on the pair (§2), and on N221204 the capsule Tr is 234 eV at bang
  against ~276 eV that HYDRA gives N210808 there (Kr22 Table I; N221204 not tabulated).
  This is the same late-time deficit. Coast time: the measured pulse's half-power fall
  is at 8.75 ns, and 9.53 − 8.75 = 0.78 ns is AS24's coast time exactly (a consistency
  check of the pulse data, not of the model).

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
| HYDRA LEH / capsule share (25.2 / 16.2 %) | 25.9 / 11.7 | 18.2 / 16.3 |
| pair difference | same | same |

(The level rows use coarse runs; N221204's standard mesh is within 0.5 eV.)

The published energetics win on Dante's shape, the one check that does not depend on
how Dante's area is counted (they capture 64–78 % of N221204's rise against 39–51 %; the
pair alike). They also win on the level in Dante's own convention on all three shots, on
the late drive, on Fig. 10, on the hot-spot transfer and on HYDRA's LEH share. The old
energetics win on the level only if the lips are taken to glow like the wall (open
area), and on the capsule share. That share is HYDRA's, a model, and `energetics=new`
calibrated α_c on a drive that was then too hot, so it is a recalibration issue
(ROADMAP 16), not evidence for the old physics. **Decision**: `examples/n221204.jl` now
defaults to `energetics=new`; `energetics=old` restores the 2026-10-24 settings. Untagged
runs and CSVs before 2026-10-25 used the old energetics. No src default changed: the
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

## Deposition-operator rays and sliver view factors (2026-10-25) *derived*

ROADMAP 14. The test-suite profile found two costs: building a `DepositionOperator`
(~6 s per 48-quad and ~22–30 s per 192-beam operator on one thread, almost all in the
visibility rays of every triangle × cone direction × beam) and the view factors of
cavities with sliver elements (the bored cavity's 20 µm waist gap: 70–80 % of its
time). Scripts in the session scratchpad `depcost/` (`cmp.jl`, `harness.jl`, `vfcmp.jl`,
`closed.jl`, `conv2.jl`, `plates.jl`); timings on the Ryzen 9 4900HS with other jobs
sharing it (ratios are from interleaved runs).

### Deposition operator: packet visibility (bit for bit)

- **Screening on the profile's cut-off does not apply.** The operator's points and
  masks are reused as aims, spots and defocus move (by `reach` with full resolution,
  beyond it at `level`), so a point that receives an exact zero at the construction aim
  is still needed after re-aiming; dropping it changes `deposit` there. Only points whose
  mask is empty (no cone direction reaches them) are dispensable, and they are already
  dropped. For the record, the rays on triangles outside every beam's `near` region
  (zero for every aim within `reach`) are 2 % of a quad's (Q15T) and 22–26 % of a
  single NIF beam's: restricting the operator to `reach` would save at most a quarter on
  192 beams and lose the level-resolved deposition beyond it. Not adopted.
- **Where the rays go** (default `hohlraum()`, 6228 opaque triangles): ~0.5 M rays per
  quad, ~0.43 M per beam; 55–65 % blocked (wall points looking back at the far wall),
  the rest lit. Per triangle and direction the answer is uniform for 98 % of the rays
  (all points lit 3–7 %, all blocked 33–46 %, back sides ~50 %, mixed 1.5–2 %). A
  single `any_hit` costs ~300 ns (≈ 55 node boxes and 4.5 triangle tests; 1.7 of those
  on triangles the ray moves away from).
- **Adopted: one packet per triangle** (`_unblocked!`, all its quadrature points × the
  cone directions facing its lit side). Interior nodes are culled by a slab test that is
  conservative over the origins' box and the directions' inverse components, leaves by
  one such test per direction, and then each ray takes its own slab test and triangle
  tests. fp subtraction and products are monotone, so no ray that passes a box is
  culled, and the slab test is monotone in box containment, so a ray passes a leaf only
  if it passes all its ancestors: every ray tests exactly the triangles its own
  `any_hit` reaches until one blocks it. Each direction's last blocking leaf is tried
  first (a cache across a work unit's triangles; 38 % of blocked rays settle there).
  The operator is identical, bit for bit (points, weights, normals, masks, destinations
  and `deposit`), on the 48-quad and 192-beam NIF operators. The facing test d·n < 0
  is now summed in a fixed order (`_facing_mask`): StaticArrays' `dot` is an `@simd`
  reduction whose rounding depends on the calling code, which flipped a direction lying
  in a back side's plane (d·n = 0 or −3e-17, a zero-weight contribution) when the loop
  was restructured.
- Tried (same answers; replaying 37k packets of six quads, `any_hit` 0.95–1.0 s): the
  adopted packet 0.62–0.65 s; one packet per direction 0.68; direction masks carried down
  the tree 0.67–0.71; without the leaf cache 0.63 (the cache barely matters once the
  packet is shared); visiting the child farther along the rays first 0.63 (the current
  order), nearer first 0.92. The occluder cache with single rays alone gave 1.1×.
- Not exact, not tried: skipping the triangle tests a ray provably moves away from (the
  ray's own facet), safe by an error bound of nine orders but not by construction.

| NIF operator, `hohlraum()` | 1 thread | 4 threads |
|---|---|---|
| 48 quads, before → after | 6.05 → 4.26 s (1.42×) | 1.61 → 1.14 s (1.41×) |
| 192 beams, before → after | 22.2 → 13.3 s (1.66×) | 5.72 → 3.54 s (1.61×) |

### View factors: slivers cut across their length

`_adapt_tri` (and `_adapt_sector`) split cells at their midpoints, which keeps a cell's
shape: a sliver stays a sliver at every depth, and the test `distance < ratio × longest
edge` holds for every thin row of it near a receiver, so a 5:1 element took ~10 rows of
cells where a square one takes ~3, and an annular sector 15:1 refined radially 16× more
than needed. Now a triangle with h² > 4.5 × twice its area (a chord rectangle beyond
~4:1) is cut across its length instead: a needle parallel to its short edge (its tip and
the two halves of the base trapezoid), a flat one at its longest edge; a sector whose
arc is more than 4× its width (or the reverse) is halved one way only. Either halves the
longest edge, so `mindepth`/`maxdepth` keep their meaning. Primitive triangles are under
the bound (walls ≤ 2.5, geodesic capsules 1.5, corner bands 2.45; latitude–longitude
polar caps reach 5.2 but are never adaptive emitters), and their midpoint children are
similar to them: **`view_factors` is unchanged bit for bit** on `hohlraum()` (4116
elements), h 0.8/0.4, h 0.6/0.3 with `sombrero = true`, with `corner_radius = 0.5`, and
with `capsule_mesh = :latlong`. Reciprocity is exact as before (the pair is still
averaged over both directions); the bored cavities change by ≤ 3.5e-4 in F and ≤ 4.3e-4
in row sums.

Accuracy against closed forms (relative error of the summed exchange):

| case | before | after | converged quadrature |
|---|---|---|---|
| 1 × 0.05 plates 0.02 apart, one 20:1 element each (parallel opposed) | 1.0e-4 | 5.8e-4 | — |
| same, 5:1 elements | −8.6e-5 | 4.0e-4 | — |
| same, square elements (unchanged) | 8.8e-4 | 8.8e-4 | — |
| 1 × 0.05 plate to a 1 × 0.02 plate on its long edge, 20:1 | 8.8e-5 | −3.5e-4 | — |
| same, 5:1 | −2.6e-4 | −5.5e-4 | — |
| same, square elements (unchanged) | −1.6e-3 | −1.6e-3 | — |
| bored cavity (0.8/0.25): gap cylinder → gap annulus | −3.65e-3 | −3.22e-3 | −3.76e-3 |
| bored cavity: annulus ↔ annulus across the gap | −2.8e-4 | −1.4e-4 | −3.9e-5 |

Closed forms: directly opposed and perpendicular rectangles; coaxial parallel annuli
and a cylinder's inner surface to its end annulus from the disk–disk formula. The
"converged" column is the same kernel at maxdepth 12, ratio 8 (before and after agree
there): the −0.38 % of the gap exchange is the receiver's polygonized arcs (16 chords
per arc leave a ~0.1 µm slit along the shared corner), not the emitter quadrature.
Slivers now carry the error of square elements of their width or less (3–6e-4, where
they were over-resolved to ~1e-4); the gap–annulus exchange is 5e-4 from its converged
value (1e-4 before) and the annulus pair 1e-4 (2.4e-4 before).

| view factors | 1 thread | 4 threads |
|---|---|---|
| `bored_cavity(hw = 0.8, hc = 0.25)` (1068 el.), before → after | 12.2 → 7.0 s | 3.19 → 1.87 s |
| `bored_cavity()` (2244 el.) | 14.6 → 9.5 s | 3.68 → 2.43 s |
| `hohlraum(h = 0.6, capsule_h = 0.3, sombrero = true)` (1888 el., identical) | 2.07 → 1.98 s | 0.53 → 0.55 s |

The gap pairs are still ~60 % of the bored cavity's time: refinement toward the shared
long edge (maxdepth 8) onto sector receivers of ~34 edges each. A maxdepth relative to
the cell's width, or coarser arcs on receivers far from the point, would be the next
step; not tried. (Superseded 2026-10-27: those pairs take the double contour with exact
arcs, "Near sector pairs by the double contour with arcs".)

## Energetics defaults and the capsule share (2026-10-25) *calibrated*

ROADMAP 16. The user's decisions: the published ("new") energetics become the defaults
wherever they are the physics Briareus uses, the old paths go without deprecation, and
the constants calibrated on the old, too reflective wall (constant-T loss at the time
since the pulse start, before "N221204 late-time drive") are recalibrated. Dante 1,
PRE 109 Fig. 10 and HYDRA's LEH share are checks here, never fit. Confidence: **H**
primary text, **L** our construction.

### What became a default, and where

| Piece | Where now | Before |
|---|---|---|
| `nif_conversion()` (0.7 t^0.12 ≤ 0.95) | src: the default `conversion` of `gold()` and `uranium()` whose albedo follows the drive (any `AlbedoModel` but a constant) | 0.75 for every `gold()` / `uranium()` |
| static walls (`gold()`, `gold(albedo = 0.85)`) | unchanged, 0.75 (VisRad-style static solves; at `solve`'s default `time = 0` the rising conversion would read 0.53) | — |
| Au-lined DU wall | the material, `uranium()` (`du_wall_loss()`, now with the rising conversion); `examples/n221204.jl` builds HYBRID-E's wall with it | `gold(albedo = du_wall_loss())` behind `wall=du` |
| DU loss factor | `DU_LOSS_FACTOR` 0.785 → **0.836** (below) | calibrated at 13 ns, constant T |
| wall history, re-emission | already `drive_history`'s defaults; the example's `history=` / `release=` switches are gone. `wall_history = false` stays in src for comparing a march with standalone `solve`s (tests, the `diagnostic_windows` tutorial) | — |
| `CapsuleAlbedo` on HYDRA's ablation-front schedule, 11 % degradation, measured backscatter | `examples/n221204.jl` (shot-specific inputs); its `energetics=`, `wall=`, `conversion=`, `backscatter=` options and `capsule=ablator` are removed | `energetics=new` |
| `CapsuleAlbedo` default α_c 0.2 | unchanged (Swadling 2025's published f_cap, not calibrated on the old drive) | — |

*Superseded in part (2026-10-25, "Capsule absorption, conversion and the ablating-wall
release")*: `gold()`/`uranium()` now convert at a fixed 0.75 and the drive-history walls
are the presets `nif_gold()`/`nif_uranium()`; the example's α_c is Young's published
albedo on N221204's K24 schedule, not the fitted 0; `release` is 1.0.

### DU loss factor at the wall's effective time (*calibrated*)

The 6.5 % peak-power equivalence (Döppner 2015) was fitted on the high-foot hohlraum
with gold's loss at 13 ns, the time since the pulse start (α_Au 0.944): the
constant-T evaluation that "N221204 late-time drive" showed to leave the wall too
reflective. `drive_history` now uses the wall's effective age, t_eff = ∫(T/T_now)^5.32
dt. From N221204's measured Dante 1 (Pak 2024 Fig. 1; the foot stores almost nothing)
it is 0.76 / 1.04 / 1.64 / 2.04 ns at 5.5 / 6 / 7 / 8 ns (*L*). The high-foot pulse's
peak (~3–4 ns at 390 TW, ~1.7 MJ; the 90 eV foot negligible) gives t_eff ≈ 2–4 ns at
Dante's peak; 2.5 ns is taken (*L*). Same setup as before (h 0.4, 390 TW × 0.86,
conversion 0.9, two-band wall source, capsule albedo 0.3; session scratch `du_cal.jl`):

| t (ns) | f = 0.785: equivalence, ΔT at equal power | α_Au | f for 6.5 % | ΔT at that f |
|---|---|---|---|---|
| 13 (old) | 6.51 %, +4.9 eV | 0.944 | 0.785 | +4.9 |
| 2.0 | 9.16 %, +6.6 | 0.880 | 0.838 | +4.7 |
| 2.5 | 8.95 %, +6.5 | 0.891 | **0.836** | +4.8 |
| 3.0 | 8.74 %, +6.4 | 0.899 | 0.833 | +4.8 |
| 4.0 | 8.36 %, +6.2 | 0.910 | 0.827 | +4.8 |

- **f = 0.836, range 0.827–0.838** over t_eff 2–4 ns (the t_eff uncertainty dominates;
  the mesh moves the ratio < 0.3 %). Check, not fitted: Döppner's T-1 pair, +4.8 eV at
  equal power against the measured +5 ± 6.
- With a less reflective wall the wall loss is a larger share of the source, so a given
  loss reduction is worth more power and the fit needs less of it. The new factor
  agrees with the independent cocktail estimates, Rosen's "about 0.83" and Schein's
  0.825. The earlier reading ("the calibrated f below Rosen's cocktail value says the
  Lindl gold loss is low at late times") was an artefact of the constant-T evaluation.
- `test/energetics.jl` "wall material factor" checks the 6.5 % at 2.5 ns.

### Capsule albedo recalibrated on the current drive (*calibrated*, at its bound)

The target is the original fit's: E_cap = 16.2 % of the coupled energy (PRE 109 App. 3,
read as a fraction of the degraded laser energy), in the full model of
`examples/n221204.jl` (calibrated CBET, wall motion, lips, bubble; `uranium()` with
f = 0.836, wall history and re-emission, 11 % degradation, measured backscatter,
`CapsuleAlbedo` on HYDRA's ablation-front schedule). Coarse mesh 0.8 / 0.4:

| α_c | capsule share | LEH share | capsule Tr 8.0 ns | Dante 1 initial / open, 8.0 ns | Fig. 10 min, zero, 8 ns, rms |
|---|---|---|---|---|---|
| 0.22 (previous default; with f 0.785: 11.7 %) | 11.4 % | 24.8 % | 307.8 | 292.1 / 315.6 | −0.185, 6.24, +0.077, 0.042 |
| 0 | 13.9 % | 23.9 % | 305.2 | 289.9 / 313.1 | −0.183, 6.21, +0.085, 0.045 |

- The share is linear in α_c to the accuracy needed (−11.4 points per unit α_c), so
  16.2 % would need **α_c = −0.20 ± 0.03** (the range from the 9–13 % degradation and
  the mesh, *L*): unphysical. **The bounded fit is α_c = 0**, with 2.3 points of the
  share unexplained; the example's default is now 0. (*Superseded*: the 2.3 points came
  from the misread capsule anchor, "The P2 swing" side result; the example now uses
  Young's published albedo on the K24 schedule, 13.8 %.)
- So the missing capsule absorption is not an albedo. As found before ("Wall
  re-emission and the capsule radius"), HYDRA's share needs absorption outside the
  ablation front (the W-doped blow-off) or a larger effective radius than the front's.
  No published radius of that absorbing surface was found, so nothing else is fitted.
  α_c = 0 contradicts Young 2020's HYDRA ablator albedo (0.16 at 6 ns, 0.3–0.4 later;
  `capsule_albedo=young`, 12.3 % before this change): the calibrated value stands in
  for absorption the schedule leaves out; it is not a material property.
- The two recalibrations barely interact: at α_c 0.22 the share moved 11.7 → 11.4 %
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
- The LEH share falls to 23.9 % (HYDRA 25.2): the extra capsule sink takes from the
  LEH as well as from the wall, so an albedo alone cannot match both shares.
- Fig. 10: the late P2 rises by 0.012 at 8 ns (more capsule sink, less smoothing), the
  minimum is unchanged; rms 0.042 → 0.046.
- Not rerun: Zylstra's pair (its difference was independent of the energetics, within
  1 eV between old and new; its level should drop by about N221204's ~5 eV) and the
  implosion metrics.

### Other constants checked

- **Wall re-emission `release` 0.77**: calibrated against 1-D diffusion with its own
  constant-T scaling, independent of the drive. Unchanged. (Later 2026-10-25: 1.0, the ablating
  wall's value; "Capsule absorption, conversion and the ablating-wall release".)
- **`LipClosure`** (speed, depth, hold): fitted on a gold basecamp history with
  `nif_conversion()` (gold walls are unaffected by the DU factor); the speed scales as
  the lip temperature's square root, so the few eV these changes move a drive would move
  it ~1 %, inside its fitted range (0.010–0.030 mm/ns). Not refitted.
- **CBET (`n221204_cbet`, `hybride_cbet`), `BubbleAbsorption`, `WallMotion`**: fitted to
  cone fractions, HYDRA absorbed fractions and radii, not to the drive temperature.
  Unchanged.
- **M-band group albedo α_M** (0.55 with thermal spots, `examples/mband_dewald2020.jl`):
  calibrated on a static solve with gold's loss at 5.5 ns since the pulse start and
  η 0.75. With the wall at its effective age (0.8 ns, mid-plateau after a rise at
  ~4.4 ns, *L*) and `nif_conversion` at 5.5 ns (0.86), cone fraction 1/3, F_M at α_M
  0.55 / 0.6 is 1.90 / 2.16 (was 2.06 / 2.33): F_M ≈ 2 needs α_M ≈ 0.57. The move is
  below the scan's 0.05 step and F_M ≈ 2 is itself approximate, so α_M stays 0.55
  (negative result). The example now uses the consistent wall; its Dante 1 reads
  277 eV against Dewald's 290 eV peak (292 before, where the too reflective wall and the
  low conversion compensated), and the share above 1.8 keV at α_M 0.55 is 0.133
  (Dewald 0.12; 0.114 before).
- **`docs/literate/leh_power_balance.jl`**: its static wall albedo was gold's loss at
  6 ns since the start (0.925); it is now `du_wall_loss()` at the effective age from
  N221204's Dante (1.07 ns), 0.875. The predicted N210307 → N210808 power saving is
  8.3 % at that albedo (measured 9.4 %; 9.5 % at 0.90, 11 % at the old 0.925); the
  Monte Carlo agreement is unchanged (2×10⁻⁴ of the source).

### What remains

- The capsule share: 13.9 % against 16.2 % at the albedo's bound; the absorbing surface
  outside the ablation front is the open input (a published HYDRA absorption radius, or
  the blow-off's opacity). (Now 13.8 % with Young's albedo on the K24 schedule; the
  remaining 2.4–4.2 points are `CapsuleAlbedo`'s `outside`, unset: "Capsule
  absorption, conversion and the ablating-wall release" §2.)
- The late-time drive after the droop (unchanged by these defaults).

## Late-time drive: the wall, the plasma and Dante (2026-10-25) *derived*

ROADMAP 16, the late-time drive ("Zylstra 2022 LEH pair, time-dependent" §2): after the
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
| Corona energy | P_corona = (1 − η)(P_L − P_bs) "retained in the hot, under-dense coronal plasma"; 15–30 % of the energy on NIF (5–10 % Nova); η 85–90 % | Kline et al., PRL 106, 085003 (2011) | H |
| | "up to 40 %" stored in the corona at NIF scale (XSN); after N091204's laser peak the XSN, f = 0.05 model follows Dante's fall, DCA f = 0.15 falls ~20–25 % lower in flux by 1.5–2 ns | Town et al., PoP 18, 056302 (2011) (OSTI 1018446) Fig. 2 | H, H-img |
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
conserved to 1e-3; at constant 300 eV it stores E = 0.0039 T^3.48 t^0.67 MJ/cm², 9 %
below Lindl's 5.2e-3 T^3.3 t^0.62 at 1–2 ns and 2.4× the static wall); and
`drive_history`'s closure (Lindl's gold loss, `release` 0.77). Each column falls from its
own peak; Dante is shown for the shape only.

T (eV) after the laser falls, capsule at R₀ (static / ablating / closure / Dante):

| t (ns) | N191118 | N210601 | N221204 |
|---|---|---|---|
| 7.0 | 298 / 273 / 266 / 293 | 299 / 269 / 263 / 289 | 307 / 278 / 271 / 291 |
| 8.0 | 251 / 239 / 227 / 281 | 290 / 265 / 257 / 301 | 311 / 282 / 275 / 310 |
| 8.5 | 191 / 198 / 180 / 262 | 261 / 243 / 232 / 292 | 297 / 272 / 264 / 311 |
| 9.0 | 171 / 177 / 160 / 238 | 214 / 213 / 195 / 278 | 270 / 253 / 242 / 302 |
| 9.5 | 158 / 163 / 145 / 221 | 187 / 190 / 172 / 254 | 222 / 221 / 204 / — |

- **No gold wall holds the drive.** From 7.0 to 9.0 ns on N191118 the ablating wall falls
  35 %, the static 43 %, Dante 19 %; on N210601 from 7.5 to 9.5 ns 30 % against 15 %. With
  no capsule and the LEHs at 0.857 r₀ (the most favourable sinks) N191118's ablating wall
  still falls 30 % (298 → 210 eV) by 9.0 ns. The ablating wall stores 620–820 kJ at the
  fall but gives back only 30–45 kJ in the next ns: most of its energy is kinetic or deep.
  So the slow fall needs a source outside the wall, as HYDRA's internal Tr (−29 to −40 eV
  from maximum to bang, against −94 to −127 eV for Briareus's capsule Tr on the pair) says
  it has.
- **The closure is good to ~10 eV.** It sits 2–3 % below the ablating reference at the
  peak and 8–10 % below it 1–1.5 ns after the fall (N191118 160 against 177 eV at 9 ns).

**The release factor on an ablating wall** (`late_time_wall.jl release`): the closure
given the ablating wall's own constant-T scaling needs `release` 1.07 / 0.99 / 0.91 for
10 / 30 / 50 % drops (static wall 0.76–0.78): the cooled surface throttles the release
less when the heated layer is expanded. In the 0-D hohlraum itself (closure with the
ablating scaling against the ablating reference, T normalized at 7 ns), the fall to
8.5–9.5 ns is matched by 0.92–1.0 on N191118, 0.97–1.05 on N210601 and ~1.06 on N221204;
0.77 gives 6–15 % too low a T 1.5–2 ns after the fall (*L*). *Calibrated* candidate:
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

- Bit for bit until the power falls (rising drives are unchanged by construction).
- The Dante deficit (initial area) at bang falls from −75 to −43 eV (N191118) and −71 to
  −48 eV (N210601); at 8.5 ns from −76 to −56 and −48 to −43; N221204 at 9.0 ns from −45
  to −39. So ~30–45 % of the late gap at bang, 6–26 % at 8.5 ns. The rest is the
  plasma (§4).
- Fig. 10 (blind) unchanged to 8.5 ns (−0.183 at 5.25, 6.21 ns, +0.085, rms 0.045);
  mean flux P2 rise→bang −0.0081 against −0.0076 (−6.1 against −6.5 µm). Split (N221204,
  % of coupled): LEH 24.2 → 24.7, capsule 13.9 → 14.0; wall re-emitted 136 → 145 kJ
  (N191118: 270 → 295 kJ).

### 4. What remains, with its data

- **A plasma source after the fall.** The published evidence (Kline 2011's definition of
  1 − η as energy retained in the corona, 15–30 % of the energy on NIF; Town 2011 Fig. 2:
  the model with more coronal energy follows Dante's fall) points to the laser-heated
  gold corona/bubble, with the capsule's W-doped exhaust (~100 kJ thermal; it is also
  where HYDRA's extra 2 points of capsule absorption likely go) and the gas (30–60 kJ) as
  smaller reservoirs. Size of order the gap (1 − η = 10–25 % of 1.6–2 MJ); its return rate
  after laser-off is unpublished. Only the thermal part (~60 kJ, ~0.3 ns) is estimable,
  and that is the 0.15 ns lag already tested. **Not modeled**: no published number
  constrains the return without fitting Dante.
- **Data that would settle it**: HYDRA's energy inventory through the coast (radiation,
  wall internal and kinetic, corona/bubble, capsule exhaust, gas) at laser-off and at
  bang; a post-shot synthetic Dante 1 for N191118 or N210601; Dante spectra (color
  against brightness temperature, the >1.8 keV share) between laser-off and bang;
  Moody 2017's cooling constant (closed access); gated LEH/bubble images after
  laser-off.
- **Decision needed** (ROADMAP 16): whether `release` 1.0 (the ablating-wall
  calibration) becomes the default. Not changed here. *Decided: 1.0 is the default
  ("Capsule absorption, conversion and the ablating-wall release").*

## The P2 swing: capsule radius, backscatter and the collision (2026-10-25) *derived*

> **Superseded in part (2026-10-26).** K24 Fig. 1's radius is not the ice–gas surface:
> it leaves R(0) at ~3.6 ns, ~2 ns before that surface can move, so it follows the
> ablator–fuel interface or the fuel's centre of mass, and the lower bound on the front
> at 6 ns is 0.84 R₀, not 0.79 ("Capsule radius schedules: published histories and the
> API (2026-10-26)"). The runs and their numbers stand.

ROADMAP 9 ("what remains"). The N221204 flux P2/P0 swing turns positive ~0.45 ns before
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
| Which radius | the fuel's inner (ice–gas) surface: it holds until the first shock breaks out (keyhole N230110, the N221204 design: 3.86–3.91 ns, Ruof 2026 Table II), and its 0.906 at 6 ns matches the ice–gas interface drawn at 6 ns, 0.92 mm against 0.985 initially (63.0 px at 68.75 px/mm in K24 Fig. 1; PRE 109 Fig. 1a, the same picture, gave 0.92 mm before) | K24 Fig. 1 caption ("the position of the ice-gas interface (red)"); PRE 109, 025204 Fig. 1a caption (ablator light grey, DT ice blue) | H (caption), L (identification) |
| HYDRA's treatment of backscatter | not modeled explicitly, "but as part of an effective total laser power multiplier and δn/n" | PRE 109, 025204 §II | H |
| HYDRA's late-time symmetry | set with Callahan's data-based model and measured sensitivity curves; "artificial multipliers on the input laser power ... to match the in-flight symmetry"; late P2 "not currently predictive" (LEH closure, plasma near the bubble, CBET), HYBRID-E needing a 1.065 peak-CF multiplier and, with 4.0 mm LEHs, symmetrized drive from 7.3 ns | PRE 109 §II; Kritcher et al., PoP 28, 072706 (2021) §III | H |
| N221204 swing calibration | "a calibration experiment for N221204 had not yet been completed"; the swing level "is currently unknown" and was reduced by hand for N230729 | K24 §II | H |
| Inner-beam glint | < 8 TW (< 2 % of peak) escapes by specular reflection, at 0.3–0.6 mg/cc | Lemos et al., PoP 29, 092705 (2022), doi:10.1063/5.0099937 (abstract) | H |

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
platforms. The HYBRID-E bubble speed stays inferred from the picket energy (`energy`
1.1) and PRE 109 Fig. 1a's N210808 tip. Reading the outer picket from Pak's delivered
pulse (`n221204_pulse()` split by DN24's cone fraction) gives 43.9 kJ in the first ns
after the turn-on (the example's "~40" was by eye), i.e. `energy` ≈ 1.69 × 0.73 = 1.23
against the basecamp's ~26 kJ (by eye): a faster bubble, the wrong direction for Fig. 10.

### Effect of each mechanism (coarse; crossing, minimum)

| Run | upward zero (ns) | minimum (at ns) | at 8 ns | capsule share | rms vs Fig. 10 |
|---|---|---|---|---|---|
| baseline (`capsule=model`, Young-shaped schedule) | 6.21 | −0.183 (5.25) | +0.085 | 13.9 % | 0.045 |
| `capsule=mesh` (shrinking mesh, same schedule) | 6.22 (+0.01) | −0.171 (5.25) | +0.056 | 14.0 % | 0.047 |
| `capsule=mesh schedule=k24` | 6.18 (−0.03) | −0.176 (5.25) | +0.062 | 16.0 % | 0.050 |
| `schedule=k24` (model) | 6.17 (−0.04) | −0.179 (5.25) | +0.091 | 16.0 % | 0.049 |
| `backscatter=uniform` (then `hydra`) (48 kJ as a uniform multiplier) | 6.26 (+0.05) | −0.196 (5.25) | +0.084 | 14.0 % | 0.042 |
| `backscatter=off` | 6.27 (+0.06) | −0.199 (5.25) | +0.084 | 14.1 % | 0.042 |
| transfer `fig9b` (HYDRA's post-CBET cone fraction) | 6.32 (+0.11) | −0.203 (5.00) | +0.080 | — | 0.044 |
| `fig9b` + `backscatter=uniform` (then `hydra`) (uniform cut) | 6.37 (+0.16) | −0.203 (5.00) | +0.079 | — | 0.043 |
| `collision=6.39` (t_c + 0.3 ns) | 6.38 (+0.17) | −0.183 (5.25) | +0.085 | 13.9 % | 0.043 |
| `collision=5.79` (t_c − 0.3 ns) | 6.10 (−0.11) | −0.183 (5.25) | +0.085 | 13.9 % | 0.052 |
| `bubble_energy=0.8` (t_c 6.38) | 6.54 (+0.33) | −0.205 (5.25) | +0.077 | 14.1 % | 0.040 |
| `backscatter=uniform capsule=mesh schedule=k24` | 6.23 (+0.02) | −0.189 (5.25) | +0.062 | 16.2 % | 0.047 |
| PRE 109 Fig. 10 (HYDRA) | 6.67 | −0.170 (5.60) | +0.095 | 16.2 % (App. 3) | — |

1. **Capsule radius** (published: K24 Fig. 1): ≤ 0.04 ns either way, and the K24
   schedule moves the crossing *earlier* (a larger capsule sees more of the late
   pole-hot pattern). Minimum within 0.012. It does not close the gap. The old
   schedule's time map is confirmed: K24's minimum falls at bang time (9.54 ns), as
   the map assumed.
2. **Backscatter as HYDRA has it**: the measured inner-cone SBS (16.8 % of the inner
   power over 5.5–6.5 ns, the shape *L*) pulls the post-CBET cone fraction 0.06 below
   Fig. 9b at 5.5–6.5 ns (0.59 / 0.56 / 0.53 at 5.5 / 6.0 / 6.5 ns against 0.65 / 0.62 /
   0.59). Fig. 10's HYDRA books backscatter as a total-power multiplier, so
   `backscatter=uniform` (then `hydra`) is the like-for-like setting: +0.05 ns, minimum 0.013 deeper. The
   measured SBS stays the default (the physics); PRE 109 itself notes that "including
   the measured backscatter time history could impact implosion shape".
3. **CBET's time dependence**: forcing HYDRA's own post-transfer cone fraction (`fig9b`,
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
5. **Bubble speed**: `energy` 0.8 (PRE 109 Fig. 1a's N210808 tip alone) gives +0.33 ns,
   the largest single lever, as found before; Pak's pulse points the other way (1.23).
   No HYBRID-E datum separates them.
6. **Glint** (Lemos 2022, estimated, not run): < 8 TW out of the inner beams, ~1.8
   points of cone fraction at peak, −0.005 in P2 at the crossing (d(P2/P0)/dCF
   −0.29 %/%), i.e. ≤ 0.03 ns *earlier*. LEH lips (computed closure): ≤ 0.005 in P2
   (before). Gas fill: in the bubble speed and opacity scalings (0.3 mg/cc on both
   platforms).

So the published inputs move the crossing by −0.04 to +0.16 ns (HYDRA's CBET result and the
uniform backscatter cut together); ~0.3 ns remains, and only the two least constrained bubble inputs
(picket fluence, collision) can carry it. Fig. 10 is itself HYDRA with power and
cone-fraction multipliers set to in-flight symmetry data of other shots (PRE 109 §II),
on a design whose swing calibration had not been done (K24): the 0.4 ns is a difference
between two models whose late-time inner-beam propagation is calibrated on different
data, not a measured discrepancy. The minimum's timing (5.25 against 5.6 ns) moves with
none of these except the Fig. 9b transfer, which moves it the wrong way.

**Side result: the capsule share.** With K24's radius as the absorbing surface (the
upper bound), the capsule takes 16.0 % of the coupled energy at α_c = 0 (16.2 % with
the uniform backscatter cut): HYDRA's 16.2 % (App. 3). The 2.3-point "missing capsule
absorption" of "Energetics defaults and the capsule share" comes from the misread
anchor: the Young-shaped schedule shrinks the absorbing radius too early (0.81 against
≤ 0.906 R₀ at 6 ns, and from 2 ns instead of the shock breakout at ~3.6 ns). Capsule Tr
drops 2–3.5 eV through the peak (276.5 against 280.0 at 6 ns, 303.3 against 305.2 at
8 ns) and the LEH share falls 23.9 → 23.1 %. The front lies between the bounds, and with
the published HDC albedo (Young 2020, 0.16 at 6 ns) the k24 schedule would give ~14 %.
Not adopted: the schedule and α_c = 0 belong to item 16 (decision below).

### Figures of merit beside the swing (derived, not fits)

| Metric | N221204 | Meaning |
|---|---|---|
| Callahan 2018 x = √(E_pk,out/(A_out ρ)) τ/R · r_cap/R (A 53 mm², ρ 0.3 mg/cc, R 3.2, r_cap 1.135) | 1.61–2.11 (picket 43.9 kJ in the first ns or 65.4 kJ to the trough at 2.1 ns; τ 8.75 ns half-power fall or 9.40 ns end) | his fitted line (P2_hs − 2.7(CF − 33) = 101.9 − 85.0x µm) gives −35 to −78 µm at 33 % CF without CBET: +13 to +29 points of cone fraction needed from transfer for a round hot spot |
| Izumi 2021 FoM ∝ (E^0.2 t′/R)^2.5 (ρ/0.3)^0.13 (r_cap/R)^1.1, relative to the HDC basecamp (t′ 5 ns, R 2.875, r_cap 0.908) | 1.62 (E 1.1, t′ 6.33 ns: outer turn-on 0.4 to mid-peak 6.73 ns); 1.39–1.90 over E 0.8–1.5; 1.72 at Pak's E 1.23 | P2 = 90(3CF − FoM) (his Eq. 2, basecamp FoM 1): round at CF_eff ≈ FoM/3 = 0.54 (0.46–0.63) |
| post-CBET cone fraction, peak | 0.53–0.59 (Briareus, calibrated two-zone); 0.59–0.65 (Fig. 9b) | +20–26 (+26–32) points over 33 % |
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
- The N221204-design 2DConA in-flight shell P2 against radius (the calibration K24 says
  was pending; Ruof 2026's N240304 surrogate gives only 270–160 µm, at 2.9 Å):
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
- Whether `backscatter=uniform` (then `hydra`) (HYDRA's own bookkeeping) is the setting for the Fig. 10
  comparison (the measured SBS shape stays *L* either way).
- Whether to shrink the capsule mesh by default (`capsule=mesh`): crossing +0.01 ns,
  minimum 0.012 shallower, late P2 −0.03; the right choice depends on where Fig. 10's
  flux is evaluated, which is not stated.

## Capsule absorption, conversion and the ablating-wall release (2026-10-25) *derived*

> **Superseded in part (2026-10-26).** §2's bounds: the `young` schedule's 0.81 R₀ at
> 6 ns lies below the lower bound, now 0.84 R₀ (K24's radius is the ablator–fuel
> interface or fuel centre of mass, not the ice–gas surface); the schedules are
> `RadiusSchedule`s, bit for bit the same ("Capsule radius schedules: published
> histories and the API (2026-10-26)"). The numbers stand.

ROADMAP 16 (and 9): the user's four decisions of 2026-10-25 on "Energetics defaults and
the capsule share", "Late-time drive" and "The P2 swing", implemented together. Dante 1,
PRE 109 Fig. 10 and HYDRA's 16.2 % capsule share are checks here; nothing is fitted to
them. Runs: `examples/n221204.jl 0.8 0.4 full <transfer> x calibrated <shot> metrics=on`
(coarse mesh, full model) on the final defaults below. Confidence: **H** primary text,
**L** our construction.

### What changed

| Piece | Now | Before |
|---|---|---|
| capsule albedo α_c (example) | Young 2020's quasi-1D HYDRA ablator albedo on the Young time map (`capsule_albedo=young`, the default; a number still accepted) | 0, fitted to HYDRA's share and stopped at its bound |
| absorption outside the front | `CapsuleAlbedo(; outside = a_out)`: α_eff = 1 − (1 − α_c)s² − a_out(1 − s²), s = R/R₀; default 0, so the default is the front's alone (`capsule_outside=` in the example) | — |
| capsule radius schedule (example) | `schedule=k24`, N221204's own HYDRA history (Kritcher 2024 Fig. 1), the default; Zylstra's pair: K24's N210808 capsule, time × bang/9.18 (*L*) | Young's N170827 front through the misread 6 ns anchor (`schedule=young` kept it; removed 2026-10-26) |
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
SBS paper, PoP 31, 042707 (2024). Cited only: a per-cone, time-dependent backscatter
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
conservation verified to ε², ~8 %).

| case | t (ns) | E/A, 1-D over H&R | m_F, 1-D over H&R |
|---|---|---|---|
| T = 2 heV | 0.5 / 1 / 3 | 0.996 / 1.001 / 1.006 | 1.047 / 1.043 / 1.028 |
| T = 3 heV | 0.5 / 1 / 3 | 0.990 / 0.996 / 1.004 | 1.063 / 1.057 / 1.040 |
| T = 2.5 t^0.18 (to 3.05 heV) | 0.5 / 1 / 3 | 0.967 / 0.968 / 0.968 | 1.065 / 1.059 / 1.038 |

- **It agrees within the solution's own accuracy**: stored energy to 1–3 %, heated
  mass 3–7 % high, against ε² ≈ 8 %. Zoning: halving the zone masses moves E ≤ 0.3 % and
  m_F ≤ 2.5 % (`fine=on`).
- With the flux limiter the reference otherwise uses (`limiter=on`) E is 9–12 % lower at
  2–3 ns (0.71–0.77 at 0.25 ns): the limiter, absent from H&R, throttles the low-density
  blow-off. The fit E = 0.0039 T^3.48 t^0.67 MJ/cm² of "Late-time drive" refers to the
  incident T_D (Marshak boundary), not the surface T, hence its difference from H&R's.
- **The release** (`late_time_wall.jl release [limiter=off]`, the closure given the
  wall's own constant-T scaling, drops of 10 / 30 / 50 %): 1.07 / 0.99 / 0.91 flux
  limited (as before), 0.98 / 0.90 / 0.84 in pure diffusion; the 0-D hohlraums gave
  0.92–1.06 (flux limited). **Adopted: `release` = 1.0** (*calibrated*), range
  0.84–1.07 (drop size and limiter). The static wall's 0.76–0.78 (`wall_reemission.jl`)
  no longer applies: Lindl's loss, which `WallLossAlbedo` uses, is the ablative one.

### 2. The capsule share with Young's albedo on both radius bounds

N221204, full model, final defaults (release 1.0, `nif_uranium()`), coarse; Dante 1 over
the initial / open LEH area:

| schedule (front at 6 ns) | capsule share | LEH share | capsule Tr 6.0 / 8.0 / 8.5 | Dante 8.0 | Dante 8.5 | Dante 9.0 | Dante at bang 9.53 | Fig. 10 min, zero, 8 ns, rms: measured backscatter / uniform cut |
|---|---|---|---|---|---|---|---|---|
| `k24` (0.906 R₀, upper bound; **default**) | **13.8 %** (14.0 with the uniform cut) | 24.7 % | 278.8 / 306.3 / 298.4 | 290.9 / 314.2 | 282.2 / 304.8 | 264.6 / 285.8 | 233.6 / 252.4 | −0.181 (5.25), 6.19, +0.080, 0.047 / −0.195 (5.25), 6.24, +0.080, 0.044 |
| `young` (0.81 R₀, near the lower bound 0.79) | 12.0 % (12.1) | 25.4 % | 282.0 / 307.9 / 299.4 | 292.3 / 315.7 | 283.0 / 305.7 | 265.0 / 286.3 | 234.0 / 252.8 | −0.184 (5.25), 6.23, +0.076, 0.044 / −0.197 (5.25), 6.32, +0.075, 0.041 |
| before (α_c 0, `young`, release 0.77, `uranium()`) | 13.9 % | 23.9 % | — / 305.2 / — | 289.9 / 313.1 | — | — | — | −0.183 (5.25), 6.21, +0.085, 0.045 |
| HYDRA / measured / Fig. 10 | 16.2 % | 25.2 % | ~300 internal | 310.1 | 311.3 | 301.5 | — | −0.170 (5.6), 6.67, +0.095 |

- With the published ablator albedo, absorption at the front gives 12.0–13.8 % between
  the radius bounds, against HYDRA's 16.2 %: 2.4–4.2 points remain. That is the share
  `outside` would carry (a_out ≈ 0.1–0.2 by the share's near-linear response, *L*); it is
  not set, since no blow-off opacity or absorbing radius is published. Alternatively
  Young's N170827 albedo is too high for N221204's thicker W-doped ablator. Neither is
  fixed here.
- **Default schedule `k24`, on physical grounds**: it is N221204's own published HYDRA
  history with the right timing (the radius holds until the first shock breaks out,
  ~3.6–3.9 ns, and reaches its minimum at bang); `young` is another shot's front shape,
  anchored on a misidentified surface (the ice–gas interface) and shrinking from 2 ns.
  Its magnitude is the upper bound (a shell that keeps its relative thickness, whereas
  ablation and compression thin it), so its capsule share is an upper bound for
  absorption at the front. HYDRA's 16.2 % played no part in the choice (both bounds miss
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
  −54 / −37 eV. The release takes ~30–45 % of the late gap, as in "Late-time drive" §3;
  the capsule changes (Young's albedo, the K24-shaped schedules) add 1–2 eV. The rest is
  the plasma source after the fall, still not modeled.
- Capsule shares on the pair: 14.9 % (N191118), 14.7 % (N210601); wall re-emitted 290
  and 192 kJ (N221204 143 kJ).
- Tests: `test/energetics.jl` (presets, `outside`), `dynamics.jl`, `retune.jl`,
  `capsule_drive.jl`, `docrefs.jl` pass on the new defaults.

### Decisions needed

*Decided 2026-10-26 (user): `outside` stays 0 until a blow-off opacity or absorbing
radius is published (a candidate for the 13.8 vs 16.2 % capsule share); the pair's
K24-shaped schedule is accepted, labelled derived, with room in the API for other
schedules (research pending, ROADMAP 16).*

- Whether to set `outside` (a_out ≈ 0.1–0.2 would close HYDRA's share on the two
  bounds) without a published blow-off opacity or absorbing radius; it is 0.
- The pair's K24-shaped schedule (N210808 capsule, time scaled to bang) is *L*; the pair's
  own radius histories are unpublished.

## Compile latency (2026-10-25) *measured*

Where users and the suite pay compilation, and what removed it. Julia 1.13.2, Windows,
4 threads, on a machine shared with other runs (wall times ±10 %; compile times are
`@timed`'s `compile_time`). Base = 3950215, new = the same plus this change, each in its
own worktree and package image. Scripts in the session scratchpad (`compile/ttfr.jl`,
`cand.jl`, `bench_history.jl`, `bench_opt.jl`, SnoopCompile 3 in a throwaway
environment, not a dependency).

**Where it went.** `@snoop_inference` over the test files: `_drive_history` was the
largest single cost, 52 s of inference in energetics/retune/bubble/dynamics/capsule_drive
and 4 s more in uncertainty/parallel_gradient. Every call site's `beams_at` closure (and
each `lpi` type, `record` closure and `times` type) is a new type, so each recompiled the
whole march: ~2 s per call site, and no precompiled copy could serve it. Next:
`visrad_target(...; remesh)` (42 s inclusive in visrad.jl: `remesh_mesh` 10 s and
`fit_surfaces` 10 s of inference over ~3,600 concrete method instances; nothing
pathological, just never precompiled), `talis_check` 7.6 s, `write_visrad_kit` 7.4 s,
`dante` on a spectral solution 7.6 s, `spectral_solve` 3.7 s, and every new objective
closure re-inferring `optimize_pointing`'s body and Optim's L-BFGS-B (~1 s). Not
addressable in src: tests that differentiate through their own closures
(`ForwardDiff.gradient(fc, x)`) get a new dual tag each, so the deposition stack is
recompiled per test closure (uncertainty.jl: `_sum_points*` 15 s, `lpi_deposit` on
duals 13 s); `optimize_pointing`'s fixed `_PointingTag` already avoids this for users.

**Changes.** (1) `_drive_history` is `Base.@nospecializeinfer` with `@nospecialize` on
`model_at`, `operator_at`, `beams_at`, `times`, `lpi` and `record`: one inferred and
compiled instance for every closure. `@nospecialize` alone is not enough: the keyword
wrapper calls the body with concrete types, and inference specializes on them anyway.
The callables' values reach `solve`, `apply_lpi` and `lpi_deposit` by dynamic dispatch
once per step. (2) `optimize_pointing` likewise on `objective` and `callback`. (3) The
precompile workload adds a wall-loss drive history (with and without CBET), a spectral
solve and `dante`, `talis_check`, a CAD `cad_target` (cavity selection, remesh, cylinder
and sphere fits) with its view factors, and `write_visrad_kit`.

**Results (first call in a fresh process, wall s / compile s).**

| workflow | base | new |
|---|---|---|
| `drive_history`, first | 10.2 / 10.0 | 0.56 / 0.35 |
| `drive_history`, each new closure | 2.1 / 1.9 | 0.17 / 0.01 |
| `optimize_pointing` (3 it.), new objective | 1.65 / 1.07 | 0.45 / 0.10 |
| same, objective with keep-outs | 3.95 / 1.38 | 2.15 / 0.13 |
| `spectral_solve` | 4.6 / 4.6 | 0.04 / 0.04 |
| `dante(dante1(), sp)` | 8.9 / 8.4 | 0.21 / 0.03 |
| `cad_target(...; remesh)` (small can + capsule) | 34.2 / 34.3 | 0.06 / 0.02 |
| `talis_check` | 6.1 / 8.5 | 0.04 / 0.00 |
| `write_visrad_kit` | 13.5 / 11.4 | 1.7 / 0.03 |
| `DepositionOperator`, `view_factors`, `deposit` | ~0 compile | ~0 compile |

Test files visrad, talis, energetics, retune, dynamics, capsule_drive, spectral,
parallel_gradient (4 threads, run side by side): compile 197 → 113 s, wall 247 → 168 s.
The package image builds in 162 → 192 s on the loaded machine (95 → 143 s quiet):
the cost moves from every first call to every rebuild, which `BRIAREUS_NO_WORKLOAD=true`
still skips.

**Hot loops unchanged.** Drive history, 1,300 elements, 24 steps, wall loss and a
recorder, warmed: min of 5 runs, four alternating pairs, base 1.80/1.77/1.97/1.89 s, new
1.98/1.79/1.82/1.84 s (noise), allocations 196 → 182 MiB. Results bit for bit (every
field of the history, `isequal`). `optimize_pointing` (default NIF hohlraum, keep-outs,
8 iterations): 0.379 vs 0.389 s, `x`, value and evaluation count bitwise equal.

**Not done.** The GUI's 0.75 s first-import tick and 0.16 s first field-colouring frame
(ROADMAP 7) were not measured (they need the Qt window); the GUI's time panel and drives
call `drive_history`/`optimize_pointing` with their own closures and now reuse the
package's precompiled copies. `view_factors` on a CAD target still compiles ~2 s at
first call when its parts differ from the workload's (sampled visibility on other
element shapes). A robust objective (`uncertainty = ...`) is not in the workload.

### Update (2026-10-25): the workload cut to the common session

The user: the workload should cover what most sessions hit, and the VisRad kit is not
that. The common session is a GUI-less script: build a `hohlraum`, NIF quads, a
`Scenario` (deposition operator, view factors, radiation model), evaluate, moments,
`optimize_pointing` with keep-outs, solve, report, clearance, and a `drive_history`. All
nine tutorials and eleven examples build a `hohlraum`; four tutorials optimize pointing;
`drive_history` is in six examples and four tutorials. Dropped from the workload:
`write_visrad_kit` (GUI export only), `talis_check` (GUI only), `cad_target(...; remesh)`
and its view factors (one example, one tutorial), `spectral_solve`/`dante` (two/three
examples, the GUI's Dante panel) and the CBET drive history (three calibration
examples). The `@nospecializeinfer` changes stay: they are what make each new closure
cheap, workload or not. The pre-existing large-model, H-matrix and ray-trace items were
not revisited.

Package image (touch a src file, `julia --project=. -e 'using Briareus'`, same machine
and session, another agent's runs alongside): 138 s → 96 s (twice 96 s; a first,
colder build of the old workload took 338 s and is discarded). First call in a fresh
process, wall s / compile s (scratchpad `ttfr.jl`, small `hohlraum(nz = 10, nphi = 16)`
for all but the optimize row):

| workflow | old workload | new |
|---|---|---|
| import | 1.8 | 1.5 |
| kept: `Scenario` (default-resolution NIF hohlraum) | 4.58 / 0.02 | 4.57 / 0.02 |
| kept: `optimize_pointing` with keep-outs, 3 it. | 0.17 / 0.09 | 0.17 / 0.09 |
| kept: `clearance`, `solve`, `report` | 0.35 / 0 | 0.34 / 0 |
| kept: `drive_history` (wall loss) | 0.43 / 0.31 | 0.44 / 0.33 |
| dropped: `drive_history(...; lpi = CBET(...))` | 0.13 / 0.07 | 1.9 / 2.2 |
| dropped: `spectral_solve` | 0.05 / 0.04 | 3.2 / 3.2 |
| dropped: `dante(dante1(), sp)` | 0.17 / 0.01 | 8.1 / 7.9 |
| dropped: `talis_check` | 0.03 / 0 | 4.2 / 5.8 |
| dropped: `cad_target(...; remesh)` | 0.06 / 0.02 | 24.7 / 24.7 |
| dropped: `view_factors` on that CAD target | 0.51 / 1.1 | 0.56 / 1.3 |
| dropped: `write_visrad_kit` | 1.3 / 0.09 | 7.7 / 6.5 |

The kept paths are unchanged; each dropped path pays its compile once per session that
uses it (CAD import the most, 25 s), against 42 s off every rebuild for everyone.

## Inner-cone SBS time histories and the late-time refraction (2026-10-25) *derived*

ROADMAP 9 and 16. Literature: `reference/litreview/inner_cone_sbs.md` (Kemp et al., PoP
33, 082705 (2026), doi:10.1063/5.0320730, accepted manuscript OSTI 3441608; the AAC 2025
talk, OSTI 3388898; Higginson et al., PRL 134, 165101 (2025); Izumi et al., RSI 95,
103510 (2024); Kemp et al., PoP 31, 042707 (2024), outer cones). Fig. 10 and Dante are
checks; nothing is fitted to them. Confidence: **H** primary text, **H-img** digitized,
**L** ours.

### What is published

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| N210808 measured inner SBS, per quad (all light: in-lens × 3.4) | 23°: 2.8 kJ, peak 2.18 TW at 5.3 ns, gone by 6.95 ns; 30°: 0.9 kJ, shoulder 0.41 TW at 5.1, peak 0.90 TW at 5.75, gone by 6.7–7.0 ns; both start as the inner power reaches peak (half of its 9 TW per quad at 4.74 ns) | K26 Fig. 2 (`kemp2026_sbs()`, `data/kemp2026_fig2.csv`; integrates to 2.87 / 0.92 kJ) | H / H-img |
| Errors | single shot up to 100 % (lens fraction, calibration); 30–50 % between N210808 repeats; the 23° sharp rise may be a processing artifact; DrD timing ±0.2 ns | K26 §II–III | H |
| Inner SBS in total (N210808) | ~30 kJ (8 quads per cone; 23°:30° = 3.1:1) | arithmetic | L |
| Modeled (pF3D on LASNEX maps) | single-species: 30° 2.42 TW at 6.0 ns (2.7× the measured peak), R 0.08 / 0.21 / 0.003 at 5 / 6 / 7 ns; multispecies: 30° 1.20 TW at 5.5 (1.7 kJ, 1.9× measured), 23° 2.0 TW at 5.75 (2.8 kJ, as measured), R(30°) 0.06 / 0.09 / 10⁻⁶; the drop after peak is in both models | K26 Figs. 2, 4; K25s | H, H-img |
| Axis- vs equator-side beams | ~10:1 in peak in-lens SBS on HYBRID-E | K26 §III E | H |
| On-axis density spike (N210808 MSH) | at z ≈ ±4.5 mm between capsule and LEH: n_e/n_c 0.01 to 6.0 ns, 0.02 at 6.5, 0.08 at 6.75, 0.15 at 7.0, 0.14 at 7.75 ns; similar in SSH; spectra suggest it forms ~1 ns earlier in the experiment | K26 Fig. 6a inset, §III D | H-img, H |
| What the spike does to the inner beams | refracts them into the Au bubble and absorbs them on axis, "reducing the incident laser power in the region where SBS backscatter originates"; without it twice the 30° power reaches the ablator at 7 ns; the 30° (closer to the spike) loses its SBS 0.25 ns before the 23° | K26 §III C; K25s slide 11 | H |
| The spike measured | N230313-002 (DT-filled N210404 hohlraum, power cut after 6.15 ns): bubble–axis neutrons at z ≈ ±3.4 mm, burn peak 8.15 ns; DT compressed to 4.3–5.6 mg/cc on axis (n_e/n_c ≈ 0.12–0.15 at 3ω, L); MSH matches the yield, SSH 3.5× high | H25 | H |
| Outer 50° cone, N210808 (Q36B, FABS, one quad) | 170 J per quad (50 J into the lens; pF3D 40 J), one burst from ~6.8 ns to 0.20 TW at ~7.4 ns, gone by ~7.9 ns (pF3D peaks ~0.5 ns earlier); noise floor 0.02–0.05 TW; 44° not studied. HYBRID-E peak reflectivity 0.012–0.08 (BigFoot 0.04–0.35) | Kemp et al., PoP 31, 042707 (2024) Figs. 1, 3, 6 | H, H-img |
| N221204 per-cone, time-resolved backscatter | not published (only AS24's 42 ± 11 kJ SBS, "early in peak power and mostly on the inner cones", and 6 ± 1 kJ SRS) | — | — |

The refraction's published effect is less inner power reaching the ablator/waist region
(absorbed on axis or deposited in the bubble, nearer the LEH), not a displaced spot: the
transverse spectra broaden by a few hundredths of k₀ (a few degrees, a few tenths of a mm
over the remaining path, below a spot size; *L*). No energy, 23° fraction or symmetry
effect is given.

### Against the example's `measured` shapes

`measured` puts AS24's 42 kJ on a shape of our own: zero to 5.0 ns, flat 5.5–6.5, gone
by 7.5 ns, the same fraction on both inner cones (16.8 % of the inner power over the
plateau). The published N210808 history, on N221204's axis, is earlier, shorter and
uneven between the cones: it starts ~0.3 ns earlier (4.73 ns), peaks in the first
~0.7 ns of the peak (23° at 5.46 ns), is mostly over by 6.7 ns and gone by 7.1 (against
7.5), and the 23° cone carries 76 % of it. The two agree on the
energy scale (N210808 ~30 kJ, N221204 42 ± 11) and on "early in peak power".

### `backscatter=history` (new option)

`examples/n221204.jl backscatter=history`: each inner cone loses `kemp2026_sbs()`'s power
per quad over its share of the delivered inner power, plus the same SRS as `measured`.
Time map: N221204's inner power reaches half its peak 0.16 ns after N210808's 4.74 ns on
K26's axis, so the curves move +0.16 ns (*L*: the shut-off follows the hydrodynamics,
timed here by the rise). Energy: scaled ×1.39 to AS24's 42 kJ with N210808's 3.1:1 split
(`sbs_energy=<kJ>` changes it; *L* for N221204's split). Result: 23° 31.8 kJ (peak 35.8 %
of its delivered power), 30° 10.2 kJ (peak 13.9 %).

N221204, full model, final defaults, coarse (0.8 / 0.4); `measured` from "Capsule
absorption, conversion and the ablating-wall release" §2, reproduced to the digit by a
fresh `measured` run in this tree:

| backscatter | min (at) | upward zero | 8 ns | rms vs Fig. 10 | P2/P0 at 5.5 / 6.5 ns | capsule / LEH share | Dante 1 8.0 / 8.5 / 9.0 (initial area) |
|---|---|---|---|---|---|---|---|
| `measured` (default) | −0.181 (5.25) | 6.19 | +0.080 | 0.047 | — | 13.8 / 24.7 % | 290.9 / 282.2 / 264.6 |
| `history` | −0.152 (5.25) | 6.20 | +0.080 | 0.043 | −0.150 / +0.039 | 13.8 / 24.2 % | 290.9 / 282.2 / 264.6 |
| uniform 48 kJ cut (Fig. 10's like-for-like) | −0.195 (5.25) | 6.24 | +0.080 | 0.044 | −0.186 / +0.032 | 14.0 / 24.1 % | 289.7 / 281.1 / 263.8 |
| PRE 109 Fig. 10 (HYDRA) | −0.170 (5.60) | 6.67 | +0.095 | — | −0.169 / −0.045 | 16.2 / 25.2 % | measured 310.1 / 311.3 / 301.5 |

- **Swing timing: unchanged** (6.20 against 6.19 ns; Fig. 10 6.67). The published
  history removes the inner power earlier (mostly 4.9–6.3 ns), when the waist is most driven,
  so the minimum is 0.029 shallower (−0.152 against −0.181; Fig. 10 −0.170), but it is
  gone by ~7 ns either way, and the crossing at 6.2 ns is set by the bubble's absorption
  and the collision (ridge stop at 6.09 ns), not by the backscatter. The minimum stays at
  5.25 ns (Fig. 10 5.6).
- **Late-time drive: unchanged** (Dante 1 and capsule Tr identical to 0.1 eV from 7.5 ns
  on): in both shapes the SBS has ended before the late peak. The 6 kJ SRS is the same.
- rms against Fig. 10 0.043 (0.047 `measured`); the capsule share is the same (13.8 %).

### Refraction as a sensitivity (`refraction=<f>`, `refraction_shift=<ns>`)

Implemented as a time-dependent inner power loss: a fraction f of the 30° power at the
spike's peak, following the inset's density (normalized), the 23° cone 0.25 ns later, on
the same time map as `history`; `refraction_shift=-1` puts it ~1 ns earlier (K26's
hypothesis). K25s's "twice the power without the spike" gives f ≈ 0.5 for the 30° cone at
7 ns (on N221204's axis 7.16 ns). The light is booked as lost: an upper bound on its
effect on the wall drive, since part of it is absorbed in the bubble and converts there;
`BubbleAbsorption` already absorbs 0.57–0.77 of the inner light in the bubble over
6.5–7.0 ns on this shot, so the spike's "refraction into the bubble" partly double
counts. The spike's density is held at its last published value (0.14 at 7.75 ns) after that
(*L*), so the loss holds to the end of the pulse: a bound, not a model. Run on the `measured` default (coarse, full model):

| | min (at) | upward zero | 6.5 / 7.0 / 8.0 ns | rms vs Fig. 10 | Dante 1 8.0 / 8.5 / 9.0 | capsule Tr peak | lost (backscatter + degradation + refraction) |
|---|---|---|---|---|---|---|---|
| `measured` | −0.181 (5.25) | 6.19 | — / — / +0.080 | 0.047 | 290.9 / 282.2 / 264.6 | 306 (8.0) | 268 kJ |
| `refraction=0.5 refraction_shift=-1` (spike at 6.16 ns) | −0.180 (5.25) | **6.09** | +0.085 / +0.099 / +0.085 | 0.057 | 279.4 / 268.8 / 247.0 | 293.6 (8.0) | 458 kJ |
| the same, uniform cut | −0.194 (5.25) | 6.14 | +0.074 / +0.097 / +0.085 | 0.053 | 278.1 / 267.6 / 246.0 | 292.2 | 465 kJ |
| `refraction=0.5` (spike at 7.16 ns) | −0.181 (5.25) | 6.19 | +0.044 / +0.086 / +0.083 | 0.047 | 280.1 / 269.6 / 247.8 | 294.4 (8.0) | 407 kJ |
| the same, uniform cut | −0.195 (5.25) | 6.24 | +0.032 / +0.083 / +0.083 | 0.044 | 278.9 / 268.5 / 247.0 | 293.1 | 406 kJ |

- Booked as lost from ~6.2 ns to the end, 50 % of the 30° power (and of the 23°, 0.25 ns
  later) removes ~190 kJ: the crossing moves **0.10 ns earlier** (away from Fig. 10's
  6.67), P2 rises late (+0.085 to +0.099 at 6.5–7 ns against Fig. 10's −0.045 / +0.046),
  and Dante falls 11–18 eV through the late peak (the late deficit grows). So the
  refraction, if it removes inner light from the hohlraum, cannot be what delays HYDRA's
  swing; on its published numbers it acts the other way. If the refracted light is
  instead absorbed in the bubble and converted (the bubble model's ridge stop already
  does this from 6.09 ns), the drive loss shrinks and the symmetry effect stays.
- `refraction=0.5` at the hydro's own timing (spike at 7.16 ns; run 2026-10-26, 4.3 min):
  ~140 kJ booked as lost. The minimum and the crossing are unchanged (−0.181 at 5.25 ns,
  6.19 ns; uniform cut 6.24), P2 after it moves by ≤ +0.003 at 8 ns, and Dante falls
  11–17 eV over 8–9 ns (290.9 / 282.2 / 264.6 → 280.1 / 269.6 / 247.8): as expected, the
  late drive deficit grows and the swing timing does not move. Neither timing of the
  published refraction closes the 0.46 ns early crossing.

### Outer cones

Kemp 2024's outer-cone history is one 50° quad on N210808 (FABS; not cone-averaged, 44°
not studied, noise at a quarter of the peak): ~2.7 kJ over the sixteen 50° quads (*L*),
10–20 % of the inner SBS, late (7–8 ns). `backscatter=history` stays inner-only: the
outer history is not published cone-averaged, and ~2 % of the 50° power for ~0.7 ns moves
the post-CBET cone fraction by < 1 point (≲ 0.003 in P2/P0, after the crossing).

### Decisions (coordinator 2026-10-25; accepted by the user 2026-10-26)

- `measured` stays the default; `history` is an option; its 3.1:1 cone split is
  N210808's (Kemp 2026), unpublished for N221204.
- `refraction` stays a bound (booked as lost; overlaps `BubbleAbsorption`).

## Closure on bored, windowed and CAD cavities (2026-10-25) *derived*

Does `closed = true` (Sinkhorn balancing, `enforce_closure!`) apply beyond the flush-LEH
`hohlraum()` of "Balancing on the curved path against a continuum Monte Carlo"? ROADMAP
CAD decision 11. `examples/mc_cavities.jl` (4 threads; ~30 min, cached).

Geometry: the bored cavity of "Partial visibility" (R 2.875, 9.4 mm long, 0.1 mm LEH bores
of radius 1.55, a 20 µm waist gap, capsule 0.9; apertures on the LEHs and the gap exit);
the same cavity from CAD-like facets (each half a closed solid of 360 strips, test
`cad_half`) through `cad_target(...; inside, remesh, fit = true)` (`:sampled`); and
`hohlraum()` of the same size with an open 1.5 mm circular window at the waist and a
closed 1 × 2 mm window (α 0.5) at z = 1, each in an open 50 µm gap. Wall α 0.9, η 0.75,
capsule α 0.3. Every model's capsule is the same geodesic mesh (`icosphere_mesh`,
n = 15, 4,500 elements); source: the CAD tutorial's two Gaussian rings on the inner wall,
integrated over each element.

References: `montecarlo_incident(...; reemit = :hit)` (the continuum) on flat copies
faceted at 0.1 mm (capsule n = 29, 16,820 elements), watertight: bores, LEH disks and gap
aperture on the walls' polygon (`open_fraction` 0 and 3.4e-5). Bored cavity: 72.7M
capsule paths, a₀ 0.016511 ± 0.000001 per unit wall power, a₂/a₀ +5.334 ± 0.018 %,
a₄/a₀ −0.472 ± 0.021 %; the CAD cavity is the same geometry (its 360-gon lies within
1.1e-4 mm of the circle). Windows: 109M paths, a₀ 0.015232 ± 0.000001, a₂/a₀
+8.595 ± 0.015 %, a₄/a₀ −0.686 ± 0.019 %. The loss through the apertures (LEHs, gaps,
open window) per unit x-ray source, from the same Monte Carlo with the apertures as
receivers: 0.36880 ± 0.00005 and 0.40237 ± 0.00006.

**Is the enclosure closed?** In the continuum, yes: the apertures cover the openings
exactly, and the curved model integrates over the true surfaces. On facets, not as the
tests and examples build it: `open_fraction` (new; diffuse rays from the elements that
meet no facet, apertures counted as walls) is 9.2e-4 on the bored cavity (12–25% on
single bore elements: 24-sided LEH disks in 20-sided bores, bores of `nphi ÷ 2` sides on
caps of `nphi`, the 180-sided gap aperture on 40-sided gap faces), 7.6e-4 on the CAD
cavity (24-sided disks in 360-sided bores), 1–5e-4 on `hohlraum()` with or without
windows (LEH disks of `nphi ÷ 2` sides in caps of `nphi`), 2.7e-3 on `hohlraum(h = 1)`.
Openings without apertures: the CAD cavity without its LEH apertures 7e-2, the
40 µm waist gap of test/cadfit.jl 4.2e-3 (6.0e-4 with an aperture). Cracks and small
openings overlap, so the fraction locates leaks (`element`) but cannot decide closure.

Model − reference (a₀ and loss in %, moments in pp), raw → balanced:

| case | a₀ | a₂/a₀ | a₄/a₀ | loss |
|---|---|---|---|---|
| bored 0.5, centroid | +0.31 → −0.02 | −0.288 → +0.010 | +0.051 → −0.107 | +1.93 → +0.24 |
| bored 0.5, sampled | −0.01 → +0.06 | −0.153 → −0.168 | +0.050 → +0.023 | −0.29 → −0.01 |
| bored 0.25, centroid | −0.21 → −0.13 | +0.079 → +0.176 | +0.053 → −0.076 | +2.26 → +0.42 |
| bored 0.25, sampled | −0.05 → −0.06 | +0.028 → +0.022 | +0.050 → +0.028 | +0.60 → +0.16 |
| CAD 0.5, centroid | +0.34 → −0.10 | −0.262 → −0.012 | +0.022 → −0.092 | +2.82 → +0.33 |
| CAD 0.5, sampled | +0.16 → +0.03 | −0.254 → −0.176 | +0.048 → +0.039 | +0.73 → +0.09 |
| CAD 0.25, sampled | −0.06 → −0.06 | −0.014 → +0.020 | +0.045 → +0.032 | +0.76 → +0.19 |
| windows 0.45, centroid | +0.35 → +0.07 | −0.310 → −0.251 | +0.010 → +0.009 | +0.26 → +0.01 |
| windows 0.45, sampled | +0.17 → +0.07 | −0.237 → −0.249 | +0.006 → +0.011 | +0.12 → +0.00 |
| windows 0.28, centroid | +0.13 → −0.02 | −0.072 → −0.088 | +0.017 → +0.003 | +0.23 → +0.07 |
| windows 0.28, sampled | +0.06 → −0.01 | −0.052 → −0.089 | +0.016 → +0.004 | +0.16 → +0.06 |
| windows 0.17, centroid | −0.04 → −0.03 | −0.018 → −0.045 | −0.005 → −0.001 | +0.08 → +0.08 |
| windows 0.17, sampled | +0.04 → −0.03 | −0.045 → −0.042 | −0.004 → −0.002 | +0.16 → +0.08 |
| bored 0.5 flat, watertight, sampled | +0.71 → +0.23 | −0.071 → −0.015 | +0.043 → +0.040 | +0.53 → −0.08 |
| bored 0.5 flat, cracked, sampled | +0.71 → +0.49 | −0.072 → +0.128 | +0.043 → −0.051 | −0.63 → −0.61 |
| bored 0.5 flat, watertight, centroid | +1.04 → +0.18 | −0.207 → +0.158 | +0.044 → −0.060 | +2.51 → +0.10 |
| bored 0.5 flat, cracked, centroid | +1.04 → +0.41 | −0.207 → +0.329 | +0.044 → −0.220 | +1.59 → −0.37 |

(5,984 / 9,064 / 5,752 / 8,252 / 6,247 / 8,445 / 14,177 elements for the curved
models.) Largest |row sum − 1| before balancing: apertures 2–6% (centroid) and 1–3%
(sampled) on the rims, the 50 µm gap faces behind the 20 µm slot ~5% in every model
("Partial visibility"), 0.04–0.5% with windows. The reference's a₄/a₀ error (±0.02 pp)
is the size of most a₄ entries.

- **Balancing applies wherever apertures close the cavity, on the curved path.** It
  brings the aperture loss within 0.2% of the reference in every sampled case (raw
  −0.3 to +0.8%) and within 0.4% with centroid visibility (raw +1.9 to +2.8% on rims),
  and a₀ within 0.07% (sampled). With sampled visibility it moves a₂/a₀ by at most
  0.08 pp (CAD 0.5 mm, toward the reference); with centroid by 0.1–0.3 pp on rims,
  toward the reference at 0.5 mm and away at 0.25 mm (+0.08 → +0.18): sampled
  visibility is what gets the rims' a₂/a₀ right, and balancing then fixes the
  energetics. **H** for these geometries and albedos.
- **Windows and gaps in `hohlraum()`:** balanced loss within 0.08% and a₀ within 0.07%
  at every mesh; a₂/a₀ is set by the wall mesh (−0.25 pp at h 0.45, −0.04 at 0.17) and
  balancing changes it by ≤ 0.06 pp. `closed = true` stays the `hohlraum()` default with
  windows, gaps and lips (apertures) and sombreros (their rows are zero; test/sombrero.jl
  checks the hohlraum's solution unchanged to 1e-12).
- **Flat facets must be watertight.** On the cracked flat bored cavity balancing pushes
  a₂/a₀ up (+0.13 to +0.33 pp) and the loss low (−0.4 to −0.6%); the same mesh made
  watertight (matched polygons) balances to the loss within 0.1% and a₂/a₀ within
  0.02 pp (sampled). The upward bias of "Monte Carlo radiosity reference for partial
  visibility" is this (both its geometries are cracked). The curved path is immune: its
  apertures and rims are exact surfaces.
- **CAD:** the remeshed, fitted cavity behaves as the primitives (raw loss +0.7%,
  balanced +0.1–0.2%). `cad_target` cannot know whether the user's apertures close the
  cavity (a missing small opening and polygon cracks give the same `open_fraction`), so
  closure is declared: `cad_target(...; closed = true)` and `visrad_target(...; closed =
  true)` (new keyword, default `false`: existing targets are unchanged bit for bit).
  The GUI's CAD path still builds raw targets.
- Not covered: thin-film windows that transmit (a transmitting aperture is a partial
  loss, not a closure) and VisRad workspaces (no public data).

## Compact capsule rows at capsule_h = h/3 (2026-10-25) *derived*

ROADMAP CAD decision 12: on the compact capsule the LEH power-balance tutorial's raw
error grew from h 0.4 to 0.25 mm. The tutorial's N210808 hohlraum (6.4 × 11.24 mm,
3.1 mm LEHs, capsule 1.13 mm), capsule_h = h/3 unless stated, Float64 view factors,
4 threads (scripts in the session scratchpad). Reference for the capsule rows: the
enclosure is closed, so every capsule element's row sums to 1, and the capsule's exchange
with the wall, the caps and the LEHs is the solid-angle fraction from its center
(cos θ_c, cos θ_L − cos θ_c, 1 − cos θ_L), exact for a sphere.

| h / capsule_h | mesh | capsule rows − 1 (mean; largest) | capsule→wall / solid angle − 1 | wall rows − 1 (mean) |
|---|---|---|---|---|
| 0.40 / 0.133 | compact | −9.2e-5; 4.4e-4 | −1.0e-4 | −2.3e-4 |
| 0.40 / 0.133 | lat-long | +1.3e-5; 1.3e-3 | | −2.9e-4 |
| 0.25 / 0.083 | compact | −3.5e-4; 6.8e-4 | −4.1e-4 | −3.0e-4 |
| 0.25 / 0.083 | lat-long | −2.3e-4; 1.9e-3 | −2.3e-4 | −1.6e-4 |
| 0.40 / 0.083 | compact | −7.2e-5 | −8.1e-5 | −3.9e-4 |
| 0.25 / 0.133 | compact | −5.6e-4 | −6.5e-4 | −1.9e-4 |
| 0.25 / 0.083, `mid = 24` | compact | −6.2e-6; 8.7e-5 | −1.9e-6 | −3.2e-4 |
| 0.40 / 0.133, `mid = 24` | compact | +1.4e-5; 5.4e-5 | +1.8e-5 | −3.0e-4 |
| 0.25 / 0.083, `level = 2` | compact | −4.5e-4 | −5.2e-4 | −3.6e-4 |

The capsule's exchange with the caps and the LEHs is within 4e-5 of the solid angles
in every case.

- **The capsule rows' error is the wall mesh's, not the geodesic triangles'.** It is the
  capsule–wall exchange, and it grows when the *wall* is refined (0.4 → 0.25 mm at either
  capsule size), not the capsule (0.133 → 0.083 at h 0.25 shrinks it). Curved pairs
  farther than `mid` (12) sizes of the larger element use one point per triangle
  (`qfar`), whose O((h/r)²) bias does not cancel; at h 0.25 that is 3 mm, so most
  capsule–wall pairs (2–6.5 mm apart) are far, while at h 0.4 (4.8 mm) most were
  mid-range (four points per triangle). `mid = 24` removes it (−6e-6), `level` does not
  (near pairs only). Lat-long elements are two triangles, two far points each, hence
  their smaller error at the same size.
- **The source-weighted loss is the wall rows', mostly the capsule's shadow.** Removing
  the capsule changes the wall–cap–LEH exchange by its shadow, which in a closed
  enclosure equals the capsule's area exactly. The centroid rule's shadow is
  0.4–0.9% too large: +5.7e-3 (compact) and +3.7e-3 (lat-long) at h 0.25, +5.8e-3 and
  +7.2e-3 at 0.4, +8.7e-3 for both at capsule_h 0.067 (h 0.4); on flat facets
  +8.3e-3 / +1.0e-2, so it is the all-or-nothing sampling of the silhouette, not the
  surfaces. Sampled visibility overshoots the other way (−1.7e-3 / −2.3e-3 at 0.4).
  Compact and lat-long differ by its scatter (wall–wall exchange per unit area, compact −
  lat-long, −1.1e-4 at 0.25 with centroid rays, −1.2e-5 sampled).
- **Tutorial errors now** (fractions of the source; 2 × 10⁷ histories of the tutorial's
  exact-geometry walk, σ 7e-5 LEH and 6e-5 capsule): h 0.4 raw LEH −0.8e-4, capsule
  −1.7e-4 (lost 6.1e-4), balanced +0.9e-4 / −0.1e-4; h 0.25 raw −2.8e-4 / −1.8e-4
  (lost 1.1e-3), balanced −1.0e-4 / +1.3e-4 (1.4σ / 2.2σ), the same on lat-long
  (−1.1e-4 / +1.2e-4) and with `mid = 24` (−1.1e-4 / +1.3e-4); the tutorial's default
  Float32 storage at 10,928 elements agrees to 1e-6. Against the tutorial's own seed-2
  walk the balanced capsule share is +0.5e-4. The ROADMAP's −3.3e-4 / −5.0e-4 (lost
  1.4e-3) were of an earlier build.
- **Decision: documented, no default changed.** Balanced results at capsule_h = h/3 are
  within 2.2σ (1.3e-4 of the source) of the exact walk on either capsule mesh; the raw
  rows' 2–6e-4 come from the far-pair rule and the centroid shadow, neither specific to
  the compact mesh. `mid = 24` fixes the capsule rows for +30–60% view-factor time but
  not the loss (at h 0.4 it raises the raw loss, 6.1e-4 → 9.1e-4: the shadow's excess
  then stands alone), so it does not earn a default change. **H**

## Capsule radius schedules: published histories and the API (2026-10-26) *derived*

ROADMAP 16, the user's decision of 2026-10-26 on "Capsule absorption, conversion and the
ablating-wall release": the pair's K24-shaped schedule is accepted, labelled derived,
with room in the API for other schedules. What published radius histories exist for the
HYBRID-E shots in use (N221204, N210808, the Zylstra pair N191118 / N210601), which
surface each follows, and what the new `RadiusSchedule` types support. No default or
result changed. Confidence: **H** primary text or table, **H-img** digitized figure,
**M** known from citing papers only, **L** our construction.

### Published radius histories and anchors

| Source | Shot(s) | Which radius | Times | Uncertainty | Use | Conf. |
|---|---|---|---|---|---|---|
| Kritcher et al., PoP 31, 070502 (2024), doi:10.1063/5.0210904, Fig. 1 (`kritcher2024_capsule()`) | N221204 capsule, and the thinner N210808 capsule, on N221204's 2.05 MJ drive (HYDRA) | not stated; the ablator–fuel interface or the fuel's centre of mass (below), as R/R(0) | 0–9.5 ns (N210808 capsule to its minimum, 0.055 at 9.18 ns) | ±0.01 in R/R(0) (0.02 at 9.25 ns), ±0.04 ns | **schedule**: `kritcher2024_radius()`, the default `schedule=k24`; `:n210808` time-scaled for the pair | H-img (curve), L (surface) |
| Young et al., PoP 27, 082702 (2020) Fig. 3a (`young2020_radii()`) | N170827 (LLNL dynamic view-factor model) | the model's ablation front | 0–7.75 ns | ±0.03 mm | schedule on a piecewise time map (`schedule=young`, removed 2026-10-26); another shot's shape | H-img |
| Kritcher et al., PRE 106, 025201 (2022) Table I (HYDRA post-shot) | N210207 / N210307 / N210808 | limb-minimum radius at peak velocity (RPV) | point anchors | — | `TrajectoryRadius` anchors, checks: v_imp 393 / 387 / 391 km/s at 8.63 / 8.9 / 8.8 ns, RPV 195 / 211 / 185 µm, bang 9.05 / 9.35 / 9.26 ns (measured 9.09 / 9.31 / 9.26 ± 0.03), so peak velocity comes 0.42–0.46 ns before bang | H |
| Kritcher et al., PRE 109, 025204 (2024) Table I (as read in `litreview/n221204_benchmark.md`; OSTI 2474067 has no full text) | N210808 / N221204 | — | point anchors | — | anchors: v_imp 401 / 383 km/s (HYDRA, recalibrated on the surrogates; PRE 106 had 391 for N210808), ablator remaining 3.8 / 5.7 %, bang 9.28 / 9.53 ± 0.07 ns measured, 9.27 / 9.61 simulated; no time or radius of peak velocity | H |
| Ruof et al., PoP 33, 052706 (2026) (OSTI 3367476, LLNL-JRNL-2019015) Table II, keyholes | N221106 (N210808 design), N230110 (N221204 design), N231016 | first shock's breakout from the ablator into the liquid D2 (VISAR) | 3.6–6.1 ns | ±0.006–0.06 ns | the time the ablator–fuel interface starts to move: equator / pole 3.714 / 3.650 ns (N221106), 3.864 ± 0.012 / 3.906 ± 0.018 (N230110; HYDRA 3.899 / 3.886), 3.575 / 3.588 (N231016); shock 1 at 32.7–34.0 km/s, merged 1–2 at ~50 km/s from 4.85 ns (N230110) | H |
| Ruof et al. 2026 Table IV and Figs. 8–9, 2DConA | N220618, N221010 (N210808 design), N240304 (N221204 design, 85.9 µm, D2 gas, lower peak power) | limb minimum and outer slope of the backlit shell (inner and outer boundaries), Legendre P0 | from R < 300 µm to stagnation (the last ~0.6–0.8 ns) | figure; HYDRA trajectory offsets −50 / −35 / +25 µm | **check only**: shell velocity 392 ± 24 / 381 ± 15 / 355 ± 5 km/s (HYDRA 401 / 408 / 344), stagnation 9.170 ± 0.061 / 9.230 ± 0.029 / 9.876 ± 0.090 ns, mass remaining 8.8 ± 1.0 % (N240304; HYDRA 10.2). A surrogate without ice on its own pulse, and R/R₀ < 0.27 (R² < 0.07, where the capsule's absorption hardly matters) | H (the extracted table's rows realigned to their labels by us) |
| Olson et al., SAND2011-8633C (OSTI 1141159) | NIC CH ConA | centre of mass of the remaining ablator (ConA's radius) | ~19–22 ns, 900 → ~200 µm | — | context: baseline HYDRA over-predicted the ConA velocity by ~18 %, matched with 85–92 % peak-power multipliers; another platform | H (slides) |
| Hicks et al., PoP 19, 122702 (2012) (1-D streaked ConA) | NIC CH | remaining-mass centre and shell thickness | — | — | not read; another platform | M |
| Callahan et al., PoP 27, 072704 (2020) §II.B (LLNL-JRNL-809797) | scaling model | shell radius of a rocket model, m dv/dt = −4πr² P_abl, dm/dt = −4πr² ṁ, driven by their T_r(t) fit | — | — | an analytic r(t), but it needs an ablation rate and pressure fitted per design and the drive; not adopted (below) | H |
| Zylstra et al., PRE 106, 025202 (2022) | N191118, N210601 | none published | — | — | bang times 9.37 / 9.48 ns, which place the pair's scaled schedule; no v_imp or RPV | H |

**Which radius K24 Fig. 1 follows.** The curve holds 1.000 to 3.50 ns and reads 0.995 at
3.75 (the N210808 capsule 0.999 at 3.50): it starts to move at ~3.6 ns, 0.2–0.3 ns before
the N221204 design's keyhole breakout (3.86–3.91 ns), as the N210808 capsule (~3.4 ns)
does before its own (3.65–3.71 ns), a common offset between HYDRA on N221204's drive and
the keyholes' pulses plus reading. The ice–gas surface cannot move until the first shock
has crossed the 65 µm of ice: at the keyhole's speeds (32.7 km/s to the 1–2 merge at
4.85 ns, ~32 µm, then ~50 km/s) that is ~5.5 ns (*L*). So the curve follows the
ablator–fuel interface (1.050 mm initially) or the fuel's centre of mass, which start
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
in use, and only for N221204's drive (the N210808 capsule on it, not N210808's own
shot). Young 2020's front belongs to N170827. The 2DConA trajectories (Ruof 2026) cover
the last ~0.7 ns of surrogates, where R² < 0.07, and are figures. The rest are point
anchors (breakout, peak velocity and its time and radius, bang).

### The API (`src/target.jl`, `src/dynamics.jl`)

- `RadiusSchedule`: abstract; `radius_fraction(s, t)` gives R/R₀. `CapsuleAlbedo(;
  radius)` and `moving_hohlraum(; capsule_shrink)` take a schedule, a `(times, R/R₀)`
  table (now stored as a `TabulatedRadius`) or a function, as before.
- `TabulatedRadius(t, r; source)`: linear, held at its ends, times strictly increasing;
  `source` records provenance. `kritcher2024_radius(:n221204 | :n210808)` returns K24's
  two curves (the N210808 one closed at 0.055, 9.18 ns), labelled published.
- `ScaledRadius(base; time_scale, radius_scale, source)`: R = radius_scale ·
  base(t/time_scale); for a tabulated base the times are scaled once and interpolated as
  a table, which is what the example did, so the pair's schedule
  (`ScaledRadius(kritcher2024_radius(:n210808); time_scale = bang/9.18)`) is the same
  table bit for bit. Its `source` says derived.
- `TrajectoryRadius(; start, peak_time, peak_radius, velocity, bang, minimum, radius)`:
  a kinematic interpolant through published anchors, for a shot with anchors but no
  history: R(0) to `start`, a cubic Hermite segment to `peak_radius` at `peak_time` with
  speed `velocity` (at rest at `start`), a second to `minimum` at `bang` (at rest), held
  after. Both segments must be monotone (speed ≤ 3× the mean, Fritsch–Carlson), else an
  `ArgumentError`. *Derived*. Checked against K24 Fig. 1, which it is not built from:

  | anchors | max \|ΔR/R₀\| (3 ns to peak velocity) | max \|Δ(R/R₀)²\| | at 5 / 6 / 7 / 8 ns (K24) |
  |---|---|---|---|
  | N221204: start 3.864 (keyhole), bang 9.53, 383 km/s (PRE 109), peak velocity 0.44 ns before bang at 197 µm (the PRE 106 rows' means), R/1.050 mm | 0.019 (7.0 ns) | 0.029 | 0.972 / 0.892 / 0.748 / 0.528 (0.965 / 0.906 / 0.767 / 0.537) |
  | N210808 capsule: start 3.714, 391 km/s at 8.8 ns, 185 µm (its own drive, PRE 106) | 0.028 (6.75 ns) | 0.043 | 0.959 / 0.861 / 0.693 / 0.442 (0.959 / 0.881 / 0.720 / 0.452; on N221204's drive) |

  N221204's time and radius of peak velocity are borrowed from the N210808 family, so the
  first row is a check of the form, not of N221204's anchors. The trajectory runs inside
  K24 by up to 0.02 from 5.5 to 8 ns. Not used by the example; for the pair it would need
  v_imp and RPV, which are unpublished, so the pair keeps `ScaledRadius`.
- **Left out: a rocket model.** Callahan 2020's rocket model needs an ablation rate and
  pressure fitted per design and the drive (a coupled march with `drive_history`, since
  the capsule's absorption depends on R). CompactLES stage 0 on our full-model drive, at
  the effective rate fitted to the remaining mass (ṁ = 2.85e5 T³), stagnates at 10.2 ns
  against bang 9.53 ("Capsule drive hand-over and CompactLES stage 0"), and LLNL's own
  HYDRA times the surrogates with peak-power multipliers of 0.825–0.86, varying 5 % shot
  to shot (Ruof 2026; Olson 2011 needed 0.85–0.92 on NIC): a rocket trajectory from our drive would be a 0.5–0.7 ns
  timing error, not a schedule. A `RadiusSchedule` subtype can take one later, e.g.
  from CompactLES's trajectory as a `TabulatedRadius`.

### Bit for bit

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
- The pair's own radius histories, v_imp and RPV are unpublished.

## Hard-sphere direct-drive illumination (2026-10-27) *derived*

Polar direct drive, Stage 0 (ROADMAP item 15; `reference/litreview/polar_direct_drive.md`
§5), the part that does not depend on the pointing parameterization (open decision 5
there). Tests: `test/direct_drive.jl`. No pointing fields or types were added.

**What exists and is used.** A bare sphere is a one-component target,
`Target(Component(:capsule, icosphere_mesh(R; n), material; surface = SphereSurface(R)))`.
`DepositionOperator`/`deposit` on it give hard-sphere illumination: straight rays, the
incident `I cos ψ` absorbed where it first lands, the far side shadowed by the sphere
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
(bit for bit the same weights; n = 30 icosphere, lmax = 30: 5.3 s → 2.0 s).

**Checks** (R₀ = 0.5 mm; icospheres n = 20 and 30, edges 0.055 R and 0.037 R; σ_l
absolute in units of the mean):
- *Skupsky–Lee separability* (J. Appl. Phys. 54, 3662, 1983), OMEGA's 60 ports
  (`omega_ports`), SG4 spots (`omega_spot(:SG4)`), aimed at the centre: the reference is
  σ_l = |a_l/a₀| G_l/√(2l+1) with a_l from a 1-D Gauss quadrature of one beam and G_l from
  the ports. σ (l ≤ 30): reference 0.0821 %; Briareus 0.0792 % (n = 20), 0.0807 % (n = 30),
  Richardson-extrapolated 0.0819 %. max_l |Δσ_l| = 3.5e-5 and 1.6e-5 (ratio 2.2, the h²
  rate). Absorbed fraction equals the enclosed fraction of the spot inside the sphere's
  disk to 1e-4. (Exploration at R = 0.43 mm: 0.1874 % reference, 0.1869 % at n = 30.)
- *Facet quadrature at l ≤ 30*: the zonal spectrum of `ylm_weights` (level 2) and the
  Legendre spectrum of `legendre_weights` (level 3) agree to 3.8e-6.
- *Temporal et al. 2014* (arXiv:1406.7762), one ring at 49°/131°, m = 4: the direct
  quadrature reference reproduces the published values with incident I cos ψ fully
  absorbed, the shift normal to the beam axis toward the equator and the a axis
  azimuthal: σ = 11.93 % (circular, a = b = R₀) and 5.23 % (elliptical, a = 2b, b = R₀)
  centred (published 11.9, 5.2); minima 1.088 % at Δ/R₀ = 0.1374 and 1.341 % at 0.0743
  (published 1.1 % at 13.7 %, 1.3 % at 7.4 %). The conventions of the review's §4.2 are
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
  Briareus has no absorption model (η⊥ = 1), so it is checked per ring at the
  re-solved geometry against the same quadrature with η⊥ = 1: max_l≤30 |Δ(b_l/b₀)| =
  2.5e-3 / 5.1e-3 (n = 20), 1.1e-3 / 2.3e-3 (n = 30), 4.8e-5 / 7.4e-5 extrapolated.

**Accuracy of element metrics.** The error is the facet averaging of the element
representation, second order in the element size; the deposition quadrature does not
set it (operator `level` 1, 2, 3 at n = 20: OMEGA max |Δσ_l| 3.5e-5, 3.4e-5, 3.4e-5). Sharp spots (order 6) need finer meshes
than SG4. For σ (l ≤ 30) to ~2 % use n ≥ 30 (18 000 elements; the 961 × 18 000 weight
matrix is 140 MB), or extrapolate n = 20/30. Ba24 evaluates σ on HEALPix nside 64 (49 152
pixels, 0.9°); n = 30 elements are about 2.1°.

**Not done (Stage 0 remainder).** The absorption weight η_a(ψ) (MN17 Eq. 29) per
operator sample; the spherical pointing parameterization (open decision 5); the Ba24
OC σ₁ example; a Literate tutorial.

## GUI on CAD: import, parts and first frames, remeasured (2026-10-26) *measured*

Measured locally (8 threads pinned to one CCX; the NIF CAD workspace of the 2026-10-18
notes: 15 CAD objects and a sphere, the 133k-facet mount and two far parts drawn only,
12 parts in the calculations; numbers only, the workspace stays outside the
repository), coarse mesh (h = 0.6 mm).

**`remesh_mesh` on the CAD capsule (H).** Since the sliver merge of 2026-10-20
(`_merge_small!`/`_dissolve!`), remeshing the capsule part (38,200 facets, sphere and
fill tube) did not finish in 3.5 min at h = 0.6 or 0.28 mm (it took ~1 s before that
commit), so the GUI's import of the workspace never completed. The capsule's small
elements hold up to 4,000 facets; dissolving one tried every remaining facet against
every neighbor element after each move, each try a pass over the neighbor (≈2,400
facets): 1.9 M passes in 160 s, almost all rejected because the new facet itself lay
beyond `rmax` of the new centroid, and 267 of 298 dissolutions ended stuck and were
undone. Now a facet is retried only after a neighbor of it moved or an element it tried
grew, the element sums are extended by the new facet (the same left folds), and the
new facet's own flatness and reach are tested before the pass over the element (the
answers are those of the full test). Capsule: 4.1 s at 0.6 mm (663 elements), 8.4 s at
0.28 mm (978); before the merge commit ~1 s, so the merge pass still costs 3–7 s on it.
Bit for bit: the remeshes of the 13 other CAD parts and two synthetic meshes at 0.28
and 0.6 mm equal the old code's (`==` on vertices, triangles and patches);
test/meshops.jl and test/cadfit.jl pass. GUI import of the workspace (coarse, no window,
first call): 12 s.

**Large small elements are no longer dissolved (H; decision 16).** Most of what was
left was still `_dissolve!`: 10.8 of 13.9 s at 0.28 mm, 5.2 of 6.7 s at 0.6 (an
instrumented copy of `_merge_small!` logging each dissolution's triangle count n,
outcome and time). The capsule's are dominated by a few of n = 10³–4.6 × 10³ at ~n²
element tests each, mostly stuck (0.28 mm: 446 of 520 stuck, 233 of them over 100
triangles; 0.6 mm: 656 of 803). Elements of more than `_DISSOLVE_MAX` = 256 refined
triangles are now left whole. The threshold, from the largest dissolutions that succeed
elsewhere (every remesh with dissolutions logged, no skip): the test meshes of
test/meshops.jl and test/cadfit.jl 112 (two of 233 over 100), the remeshes of
test/hmatrix.jl and test/partial_visibility.jl 49, the cad_hohlraum tutorial (0.5/0.15
and 0.25/0.1 mm) 24, the 13 other CAD parts of the workspace 60, `sphere_mesh(0.9;
ntheta = 60, nphi = 120)` at 0.25–0.6 mm 66, and a 40k-facet sphere
(`sphere_mesh(1.1; ntheta = 100, nphi = 200)`) at 0.6 mm 237; the sliver tubes (test
meshops' 360 strips, the GUI's 2,000) and `cylinder_mesh(...; nphi = 400)` have none at
0.25–0.6 mm. A dissolution that ends stuck is undone, so the result changes only where
a dissolution over 256 would have succeeded: all of the above are `==` the old code's
(46 meshes compared directly, the test and tutorial meshes by that criterion), and the
test files pass. Capsule time against the threshold (one run each, s, 0.28 / 0.6 mm):
none 13.9 / 6.7, 512 2.28 / 1.21, 256 2.07 / 1.01, 128 2.06 / 0.97, 100 2.02 / 0.98; at
256 the dissolutions take 0.1 / 0.05 s and clustering is the rest. 100 would change two
test meshes (n 108, 112), 128 the 40k-facet sphere (237) and more of the capsule's
successful dissolutions (n 131–142 and 192–255); 512 costs 0.2 s more on the capsule.
256 keeps every public case at no measurable cost over 100. Timed
three times after a warm-up (8 threads, the machine shared): before 12.4–12.7 s at
0.28 mm and 4.1–4.2 s at 0.6 (the 8.4 s above was a less loaded run); after 2.0–3.3 s
and 1.0–1.5 s.

The capsule mesh changes where a dissolution over 256 succeeded in the old run (0.28 mm:
14, n 318–2,094; 0.6 mm: 8, n 260–1,767): 36 elements (0.17 % of the capsule's area) become
48 at 0.28 mm, 978 → 990 elements; 21 (0.18 %) become 28 at 0.6, 663 → 670; all in one
thin band of the sphere. Elements under 0.3 size²: 445 → 457 and 620 → 627 (at 0.6 mm
the angle, not the size, bounds the capsule's elements, so almost all are "small": the
merge cannot help there). Unchanged: max tilt 9.995° / 9.997°, max reach 0.9998 /
0.9996 size, every centroid on its element, max diam/√A 22.3 / 29.9 → 29.8, the area
to 3e-14, closed and conforming (no open or non-manifold edges). Output triangles
+7 % at 0.28 mm (59,958 → 64,280), +1.6 % at 0.6. The skip also applies to finer CAD
spheres: `sphere_mesh(1.1; ntheta = 140, nphi = 280)` at 0.6 mm goes 1.6 → 1.2 s and
238 → 256 elements, at 200 × 400 5.3 → 2.2 s and 244 → 274 (a facet-count threshold
skips more on finer faceting). Scripts: scratchpad `dsk/` (`instr.jl`,
`survey_cad.jl`, `survey_tests.jl`, `survey_more.jl`, `time_capsule.jl`,
`cap_stats.jl`, `mock.jl`).

**The window's GUI-thread stalls (H).** Scripted Qt window (scratchpad
`gui_timing.jl`: a coarse hohlraum's drive, coloring by a field, a selection, the
workspace import, the largest part left out and put back, the template and the import
again), the longest gap between ticks after each step, first and warm, two runs each;
"before" is the GUI as of 576587b with the remesh fix. Most of the first-call stalls
were not the GUI's own work but compilation on the GUI thread, of two kinds:
(1) `_run_job(job, spawn)` was specialized on each job's closure, and inferring its
inline branch inferred the whole job (the import's reading and meshing, a part's
remeshing, the drive's scenario) on the GUI thread before spawning it; it is now
`@nospecializeinfer` on `job`, so the job compiles on its worker. (2) A worker's result
is applied on the GUI thread by a closure compiled at its first call there; the worker
now compiles it first (`_deliver_compiled`, a `precompile` of the commit for the
result's type; a lookup once compiled). Compilation the GUI thread still needs waits
behind the workers' (the code generator's lock), so a few ms of it can stall for
hundreds while a worker compiles the CAD remesh.

| step (window) | before | after | warm |
|---|---|---|---|
| drive requested (call) | 0.39–0.40 s | 0.01 s | — |
| first coloring by a field | 0.36 s | 0.21–0.24 s | 0.03 s |
| first selection | 0.22–0.23 s | 0.20 s | 0.03 s |
| import requested (call) | 0.72–0.73 s | 0 | 0 |
| import's commit tick | 2.2–2.4 s (commit 0.95–1.0 s) | 1.5–1.6 s (commit 0.27–0.28 s) | 0.14–0.19 s |
| a part left out (call; tick) | 0.35–0.38 s; 0.36–0.38 s | 0.09 s; 0.17–0.21 s | 0.03–0.05 s; 0.07–0.16 s |
| a part put back (call; tick) | 0.86–0.87 s | 0; 0.13–0.17 s | 0; 0.12–0.15 s |

Without a window (scratchpad `gui_import_ticks.jl`, a session with a queue, first
call): import call 1.33 s → 8 ms, its commit 0.85 → 0.05 s; leave out 0.33 →
0.05–0.09 s, put back 1.16 s → 0. Also: `import_visrad_workspace!` is
`@nospecializeinfer` on `include` and `done` (each new `done` closure recompiled its
body, 0.85 s); the GUI workload leaves its disk part out and puts it back through a
window-style queue and starts an import of a missing file (the call's path), and
precompiles the import's `set_beams!` call (0.17 s of its commit). The field colorbar
warm-up at the window's first frame now sets the field's colormap, scale and label and
a second scale at the next frame (the colormap and gradient texture updates were
compiled at the first coloring). Left: the import's commit tick is Makie's first
drawing of the CAD scene (the ghost and highlight meshes, the beams' line segments,
camera fitting: compute-graph edges of 50–80 ms each) and the shell's first rows of a
VisRad part, ~1.3 s; the first field coloring and the first selection, ~0.2 s each, are
GL-side first drawings too. Drawing the highlight from a degenerate triangle at the
first frame, as the ghost is, did not remove the first selection's stall (0.57 s in
one run) and was dropped.

## Per-cone M-band in the drive history (2026-10-27) *derived*

ROADMAP 13. Until now the laser M-band entered only the `unfold=on` spectra of
`examples/n221204.jl`, with the per-cone conversions of `nif_mband()` (outer 0.20, inner
0.042; "Spectral radiosity, M-band and a synthetic Dante" §2) applied through a proxy:
each element's cone shares taken from the *pre-CBET* deposition of `beams_at(t)`, times
its post-CBET absorbed laser power.

**Change.** `mband_power(m, op, beams; lpi, time)` deposits each cone's *effective* beams
(after CBET, backscatter and the bubble) through `lpi_deposit`, so light the bubble's
ridge stop puts off the operator converts at its own cone's fraction. `drive_history(...;
mband, spectral)` does this at every step with that step's beams, splits the step's gray
solution by `spectral_solve` (default edges [0, 1.8 keV, ∞]) and records `mband_moments`
(Legendre moments of the observed flux above 1.8 keV) and `mband` (laser-produced M-band
power per part). Off by default; it is a diagnostic: the M-band is part of the converted
x-rays and, thermalized, the group sums are the gray solution's, so the gray history does
not depend on it.

**Checks.** (`test/spectral.jl`, "drive history M-band", coarse 10.1 × 5.75 mm hohlraum,
constant-albedo gold, an outer → inner transfer growing in time.) The gray records with
`mband` equal those without (`==`); the history's M-band moments equal a static solve at
each time with the M-band summed cone by cone by hand from the post-transfer beams to
1e-10, its M-band power to 1e-12; with unit fractions the per-cone sum is the deposition
(1e-12). Against the pre-change code (its `src/history.jl` and `src/spectral.jl`
re-included): a wall-loss history with a falling pulse, CBET, Dante record, `ylm` and
`keep_albedo` is `==` field by field with and without `mband`, and `mband_power` without
`lpi` is `==`. On the N221204 coarse model (`0.8 0.4 full calibrated table`, the comparison
times), the written P2/P0, P4/P0, capsule Tr and Dante columns are byte-identical with and
without `mband=on`.

**Effect** (N221204 full model, h 0.8 / 0.4, comparison times; α_M 0.55 thermalized,
spots at the local Tr, `nif_mband(materials = (:Au, :DU))`):

| t (ns) | capsule M-band P2/P0, per cone post-CBET | pre-CBET weighting | difference | M-band share of capsule flux | laser M-band / laser |
|---|---|---|---|---|---|
| 2.0 | +0.237 | +0.234 | +0.0026 | 0.123 | 0.168 |
| 4.0 | +0.331 | +0.330 | +0.0016 | 0.097 | 0.163 |
| 4.5 | +0.269 | +0.268 | +0.0011 | 0.102 | 0.155 |
| 5.5 | −0.014 | −0.011 | −0.0025 | 0.102 | 0.127 |
| 6.5 | +0.274 | +0.272 | +0.0020 | 0.105 | 0.163 |
| 7.0 | +0.342 | +0.341 | +0.0015 | 0.108 | 0.180 |
| 8.0 | +0.371 | +0.371 | −0.0000 | 0.103 | 0.194 |
| 8.5 | +0.209 | +0.210 | −0.0001 | 0.128 | 0.194 |

- **Negligible**: |Δ| ≤ 0.0026 in the M-band P2/P0, mean +0.0003 over 4.5–8.5 ns; the gray
  P2/P0 swing, its crossing, capsule Tr and Dante are unchanged by construction (0, not
  small). The proxy's premise holds: on this pointing nearly every element is lit by one
  cone, so its pre-CBET cone share is ~0 or 1 and the post-CBET laser power on it already
  carries the transfer; the bubble moves absorbed inner light to the outer footprint beams
  in both forms. Only the overlaps and the ridge stop differ.
- The M-band P2/P0 itself swings with the cone balance: +0.33 at 4 ns, −0.01 at 5.5 ns
  (the CBET-boosted inners at peak), +0.37 at 8 ns. Dewald 2020's Mo-sphere M-band P2/P0 on
  his shot, mostly > 0 and +0.15 at 6 ns (Fig. 10, **H-img**), is a check of sign and
  magnitude only (another shot and pointing; not fitted).
