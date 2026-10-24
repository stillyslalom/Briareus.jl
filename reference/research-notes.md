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
  constraint (ROADMAP 10) is what can fix the temperature scale.

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
- **f = 0.785, calibrated** to the 6.5 % equivalence: Döppner's high-foot hohlraum
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

- **Schedule (example, *L*)**: HYDRA's ablation front, 0.92 mm = 0.81 R₀ at 6 ns (PRE 109
  Fig. 1a), along the shape of Young 2020 Fig. 3a's N170827 front, its time axis mapped
  piecewise linearly (Young's 0.81 R₀ at 5.75 ns → 6.0 ns; its end, 0.23 R₀ at 7.75 ns →
  N221204's bang time 9.53 ns): R/R₀ = 0.94 at 4.2, 0.86 at 5.0, 0.81 at 6.0, 0.66 at 7.0, 0.51 at
  8.0, 0.31 at 9.0 ns.
- **Calibration**: α_c = 0.22 gives E_cap = 16.2 % of the coupled energy (285 kJ of
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
  independent agreement (Dante's M-band share was not fitted).
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
0.3 E (**L**, widths unpublished); `color_temperature` fits a Planckian to the channels
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
  block (G = AF symmetric) halves it again.

## N221204 late-time drive (2026-10-24) *derived*

ROADMAP 16. Two results of the N221204 benchmark (`examples/n221204.jl`) disagreed about
the late drive. (a) Default energetics: capsule Tr peaked at 301.8 eV at 5.5 ns and was
flat to 8 ns; synthetic Dante 1 (open LEH area) peaked at 306 eV at 6.25 ns; measured
Dante 1 rises to 312 eV at 8.3 ns. (b) `energetics=new`: capsule Tr 334 eV at 8 ns and an
LEH share of 33 % against HYDRA's 25.2 %. Neither Dante nor Fig. 10 was used to fit
anything below; both stay checks. Confidence: **H** primary text or figure vector data,
**L** our construction.

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
