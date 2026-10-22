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
  0.2 mm walls). **Visibility:** with bored LEHs and a
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
- **Tori (fillets), not done.** A torus fit (7 parameters) would turn the CAD fillets
  that stay unfitted or split into small cones/spheres into one region each; on the NIF
  halves they are ≤ 1.3 mm² apiece and their flat elements are within the faceting, so
  the gain is a smaller component count, which `minfraction` already gives. **M**
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
