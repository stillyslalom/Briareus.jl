# Research notes

These are literature findings gathered from public sources (October 2026) that set
Briareus models and defaults. Confidence levels:
- **H**: primary text read.
- **M**: abstract or secondary source only.
- **L**: inference.

Numbers marked *derived* are our own computations.

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
- JSON workspace files (VisRad 21.1.0).
- STEP CAD import.
- A "Laser-Target Capture Clearance Report" (capture fraction and maximum incidence
  angle) exportable as ASCII.
- Multi-column ASCII time tables.
- Per-element flux exports for rad-hydro input (MacFarlane 2003).

Not documented in the manual: the full JSON schema of `.vrw` workspaces. Sample files
are still needed for that.

## State of the art (pointing design and tools)

- **LLNL view-factor tools.**
  - A "dynamic view factor" with wall and spot motion taken from rad-hydro (Young, DPP
    2018–19). **M**
  - MacGowan's 3D view-factor metric (HEDP 40, 100944, 2021): 1% mode-1 flux gives
    about 100 km/s hot-spot velocity, and deliberate ±1.2% mode-1 drive was used to
    cancel capsule asymmetry. **H/M**
- **China.** IRAD3D adds time-dependent albedo and spot motion (arXiv:1606.00934,
  1703.01579). **M**
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

- **Workspaces:** `.vrw`, JSON since v21.0.0.
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

### Workspace `.vrw` (JSON) and `.bti` (not implemented)

- `.vrw` is JSON since v21.0.0 and "can be readily modified using scripting tools"
  (VR/files/workspace_files.html). **H**
- The only key the manual names is "Beam cone F number" (v21.1.0 notes). The schema is
  undocumented. **H** (that it is undocumented)
- `.bti` (File | Export Laser Beam Data → Beam-Target Intersection Points) gives each
  beam's intersection points in the component's and the chamber's frames, in cm and
  degrees (VR/files/laser_target_intpts_file.html). **H**
  - Format ID 2 adds the beam name and port angles.
  - An optional block adds a(z), b(z), n(z), fSQ, DLW, DLH, M.
  - The column layout is undocumented. Not implemented.

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

