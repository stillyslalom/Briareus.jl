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
  - Γ, ΔL(z), M(z) and z₀ are not given. a(z) and b(z) are negligible within 0.6 mm,
    but f_sq is not unless Γ is large: z/√(z² + ab) = 0.5 at 0.3 mm and 0.76 at 0.6 mm.
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
