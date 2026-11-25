# Research notes: Facility data, pointing conventions and VisRad

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

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

The manual does not document the `.vrw` layout. Two sample workspaces (2026-10-17)
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
- **How the 0.5° enters.** The document does not give a formula. Meaning (user, from facility
  practice, 2026-10-27): an insertion or assembly error of the target component. The
  beams must still clear with the component rotated up to 0.5° off nominal in the least
  favorable direction, whether the target is inserted rotated or the component is
  assembled off-square. It does not apply to the 70° incidence limit (user,
  2026-10-27). Briareus applies it to the clearances and, as a conservative extension
  of the rule text (which puts it on the clearance slides only), to the 0.5 mm capture
  edge margin. Briareus applies it as a lateral
  allowance r tan 0.5°, with r the distance from the rotation's pivot (TCC by default,
  or the component's mount); the pivot per component remains a modeling choice. **H**
  (meaning), **M** (pivot). VisRad's capture/clearance report has no angular term: a
  VisRad 21 probe (2026-10-27) clears at the fixed 0.5 mm at 20–120 mm from TCC. A
  VisRad user adds the 0.5° to "Min clearance" by hand, if at all. VisRad's incidence
  angle is the f-cone's steepest ray (axis + 1.38° at f/20.755), which `talis_check`
  follows. **H** (VisRad's behavior)

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
    2 mm / 5 mm. TaLIS states "< 2 mm / < 5 mm"; Briareus uses the strict form.
  - Cones and tori (*interpretation*, user's statement of the rule's purpose, 2026-11-25;
    confidence medium). The concern is specular 1ω reaching another beamline in the far
    field, so the limits apply to the unsigned radius, concave or convex. `talis_check`
    judges each sub-triangle in the footprint on its largest principal radii: it scatters
    in one dimension if the smaller is < 2 mm and in two if the larger is < 5 mm. A cone
    has the radius ρ/cos α at distance ρ from its axis, a torus the minor radius r and
    |R/cos ψ + r|. The document names only cylinders and spheres.
  - Free-form meshes (*engineering default*, 2026-11-25; confidence medium). The radii
    of a facet without a tagged surface come from the turn of its vertex normals along
    its edges (the per-face shape operator of Rusinkiewicz, 3DPVT 2004), where a vertex
    normal averages, by corner angle, the facets around the vertex within 40° of the
    facet. Measured on Briareus meshes: cylinders of 12 and 24 facets exact to rounding; a
    6 mm icosphere 5.77 to 6.29 mm; a 4 mm latitude–longitude sphere 3.15 to 4.57 mm, the
    extremes at its poles; a flat plate and the flange of a stalk-on-flange solid flat.
    A facet is taken as scattering only with a quarter added to its radii, a fit residual
    ≤ 0.1 and a normal spread ≤ 25°; the quarter covers the 20% seen on spheres. These
    three numbers are choices of this package and are not facility criteria.
  - Keep-out map (*engineering default*, 2026-11-25; confidence medium). `corridor_map`
    samples a sphere about a point, or a plane, with the distances and defaults of the
    rules above: the 1% contour with 0.5 mm for 3ω light, the 10% square with 0.2 mm for 1ω
    light, the same position errors and rotation tolerance. The volume of a line of sight
    is a choice of this package: a cylinder of the collimated field of the definition
    (60 mm for Dante 1), or for a pinhole imager the cone from the detector field at the
    object (24.6 mm for the SXI) to the pinhole, with no clearance added. The TaLIS
    document gives no keep-out volume for a line of sight. Checked against closed forms
    (a cylinder, a cone with its apex at the center, the 1% contour of a super-Gaussian
    spot with its f-cone, the near side of a 1ω square on the beamlet axis). Measured with
    192 NIF beams at their default aims, Dante 1 and both SXIs: nothing open at 20 mm from
    TCC, 33% of the solid angle at 50 mm, 60% at 150 mm.
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
  - Incidence does not take a rule tolerance (2026-10-27; it took +0.5° and the widened
    footprint before); the capture margin keeps the lateral allowance. Specular-reflection angles get −1° (a tilt by α turns a
    reflection by 2α).
- **1ω footprint model.** *Derived*, **M**.
  - Each beamlet's 10% 1ω region at its 3ω focal plane is a 22.7 mm square: the User
    Guide's 22 mm, grown 0.35 mm per side so that its near edge is the TaLIS 4.45 mm.
  - It is displaced along φ̂ to the side of the quad axis where the beamlet sits (User Guide
    Fig. 4-12 top view).
  - Away from that plane it is the cone converging to the 1ω focus 460 mm beyond the
    3ω focus, with its center on the chief ray from the lens center 7.7 m upstream. The
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
    `(:dimpled, r, half_side)`. The User Guide (§5.6.2) also needs dimple
    half-side/radius ≥ 0.25, so a dimpled surface without `half_side` is a `:review`
    (2026-10-28; before, it passed on the radius only).
- **Distances over whole triangles (2026-10-28).** An independent review found a 100 mm
  plate at z = 2 mm, crossed by a vertical beam but meshed as one quad, reported as
  passing 3ω clearance by 10.6 mm: the rules were evaluated only at triangle vertices
  and `4^level` sub-triangle points (substrate edges every 20 µm). Every distance in the
  rules is now minimized over each whole triangle (or edge) by a branch-and-bound on a
  Lipschitz bound (`envelope_distance` √(1 + cone²); 1ω footprint 1 + its taper; the
  tolerance allowance tan 0.5°). The reported value is a lower bound, within
  max(1e-4 mm, 1% of its distance from the rule threshold), so a status never
  flips on the bound. On the plate the check now fails at −0.719 mm, the analytic
  value. One analytic test was itself a sampling artifact: with the pivot on a plate
  edge point, r·tan 0.5° grows linearly along the edge and the distance quadratically,
  so the minimum is 14 µm off the vertex and 6e-5 mm lower. **H** (method; checked
  against closed forms and a 0.1 µm grid).
  - Cost: the 48-quad hohlraum test geometry takes 0.9 s (level 1) and 1.4 s (level 2)
    on 4 threads, against 0.2 s and 0.8 s with vertex sampling. Most of it is rims
    nearly concentric with a beam, where the distance is almost constant along an arc.
  - `clearance` minimizes its envelope distance the same way (the same routine). Its
    reachability rule is decided per `4^level` sub-triangle: one counts whole when any
    vertex or its centroid is reached, asked only for sub-triangles that could hold the
    minimum. Peak intensity stays a maximum over the reached vertices, plus the point of
    least clearance. On the 48-quad hohlraum (caps as hardware) the cost is unchanged
    (3.8 s quick, 11.6 s standard on 4 threads; the power ray trace dominates) and the
    worst cap clearance drops 7.5 µm, the minimum the vertices missed.

#### LEH edges and the capsule (2026-10-27)

The GUI's default NIF hohlraum (10.1 × 5.75 mm, 3.1 mm LEHs, 0.9 mm capsule, the
`NIF_CONE_POINTING` aims) failed the TaLIS check: 72 fails.
- 64 were `clearance_3w` on the end caps at the LEH rim. Worst: Q33T on `cap_top`,
  worst case 0.020 mm against the required 0.5 mm (margin −0.480 mm; nominal 0.205 mm).
  An independent pure-Julia check gives the same: the best-focus 1% ellipse from the HWHM
  table, grown by the f/7.9 cone, against the exact rim circle, nominal 0.207 mm, worst
  margin −0.477 mm. The difference is the 48-gon chords, which lie up to 3.3 µm inside
  the circle. **H** (geometry).
- 8 were `incidence_3w`: the 23.5° quads on the curved wall. Nominal 69.29° over the
  nominal 10% footprint; worst case 70.57° against 70°. That is 70.07° over the
  footprint grown by the errors and the lateral rule allowance, plus the 0.5° that
  Briareus adds to the incidence. With `tolerance = 0` these rows pass.
- A ray trace at nominal pointing deposits no power on the caps.

As a closed-form bound, the inner (23.5°) quad's best-focus 1% contour has a semi-axis of
1.145 mm along φ̂. With the axis through the LEH center, that leaves
1.55 − 1.145 = 0.405 mm to the lip along φ̂, before the f-cone, errors or tolerance.
No 3.1 mm-LEH NIF design meets "1% contour > 0.5 mm from the LEH lip".

The rule text does not separate the LEH of the target itself from hardware. It states
"clearance to target features", and every 3ω violation leads to the same consequence: the
configuration "must be assessed for blow-by by the backscatter working group". That is
an assessment, not a prohibition. Interpretation adopted (**assumption**, not stated by
TaLIS; confirm with the facility):
- `talis_check(...; leh)` names a target's entrance-hole edges (hohlraum end caps, LEH
  lips). Their 3ω clearance is checked, but a violation is reported as `:review` with an
  "LEH clearance" note, not as `:fail`. Hardware proper (shields or sombreros, stalks,
  diagnostics, closed windows) still fails.
- The GUI passes its checked hohlraum caps and lips as `leh` and leaves the capsule out
  of the default TaLIS hardware. Light on the capsule is drive physics (keep-outs such
  as `SphereKeepOut`, the deposition), and the rule text does not mention capsules.
  This is an engineering default; a user can mark the capsule as hardware again.

Incidence was left as it is. Slide 3 defines capture "by a target feature": the beam is
fully captured if it is "incident on the target substrate" with the 10% contour ≥ 0.5 mm
from the substrate edge, "and incident with an angle <70 degrees wrt the target normal".
It is not restricted to foils or capture substrates, so a hohlraum wall the beam lands on
falls under it. The rule-tolerance clause ("up to 0.5 degrees
metrology/alignment/positioning tolerance") appears only in the clearance definitions
(slides 4 and 5), not in the capture definition. Briareus added the 0.5° to incidence and
the lateral allowance to capture as a conservative extension until 2026-10-27, when the
user confirmed that the tolerance (a component rotated off nominal) does not enter the
70° incidence limit; the +0.5° was dropped, and the capture margin keeps the lateral
allowance.

All from the manual. **H**

- **Workspaces:** `.vrw`; the manual states JSON since v21.0.0, real 20.0.0/21.0.0 files
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
- **Capture/clearance report defaults:** equal to the TaLIS rules. The report applies the
  margins as fixed widths, with no 0.5° term, and measures incidence on the f-cone's
  steepest ray ("Capture/clearance report", 2026-10-27).

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
- The NIF 1ω prepulse level and shape (User Guide §4.4.3 gives only when it passes
  10⁸ W/cm²).
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
- Delimiter quoting, preamble lines, and whether the SST adds columns are not documented.
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
  - Spot size can be defined at 36.8% (1/e), 50% or 5% of peak; the conversion is
    r_f = a (ln 1/f)^{1/n}. **H**
  - NIF CPP radii have been HWHM since v10.0.0 (VR/intro/revision_versions_800to1020.html). **H**
  - Away from focus VisRad uses a(z) = a₀√((z/z_A)² + 1) and n(z) = n₀[1 + 0.2 z/z₀]. This
    differs from Briareus's f-cone convolution. **H**
- **Rotation.**
  - The rotation is that of the minor axis from vertical, positive clockwise; zero has
    the major axis "nearly horizontal", "coincident with the upper edge of the square
    beam as it enters the target chamber" (same page; revision notes). **H**
  - Which side "clockwise" is viewed from is not stated. VisRad 21's deposition settles
    it as the beam's-eye view (looking along propagation): ψ_major = 90° − ρ, Briareus's
    default `rotation_sense = :beam` (`:target` gives 90° + ρ). **H** (2026-10-27; see
    "VisRad probe kits (2026-10-27): results")
  - NIF spots (major axis along φ̂) have ρ = 0 either way.
- DPR ("On" splits an OMEGA beam into two cones 85 µm apart) is in the CSV. **H**
  Briareus ignores it with a warning.

### Phase-plate table

- Table at VR/power_sources/phase_plate_params.html. Radii are in µm: the "(mm)"
  headers are "µm" in the Symbol font. **H**
  - OMEGA, OMEGA EP: 1/e radii.
  - NIF CPPs: HWHM.
  - LMJ: radii at 3%.
  - GEKKO: definition not stated; 1/e assumed. **L**
- The OMEGA table has a "# radii at 5%" footnote that no row references.
- VisRad's NIF CPPs agree with Briareus's measured fits to ≤ 2% in 1/e radius:
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

- The manual states that `.vrw` is JSON since v21.0.0 and "can be readily modified using
  scripting tools" (VR/files/workspace_files.html). Real 20.0.0 and 21.0.0 workspaces
  are keyed text instead; `read_visrad_workspace` reads that format (next section). A
  JSON workspace (perhaps 21.1+) is rejected with an error. **H** for the text format
- The only key the manual names is "Beam cone F number" (v21.1.0 notes).
- `.bti` (File | Export Laser Beam Data → Beam-Target Intersection Points) gives each
  beam's intersection points in the component's and the chamber's frames, in cm and
  degrees (VR/files/laser_target_intpts_file.html). **H**
  - Format ID 2 adds the beam name and port angles.
  - An optional block adds a(z), b(z), n(z), fSQ, DLW, DLH, M.
  - The column layout is undocumented; Format ID 3 as VisRad 21 writes it is recorded in
    "VisRad's beam profile law" (2026-10-27). No reader.

### VisRad workspace files (2026-10-17)

The layout was read off two real workspaces: a NIF hohlraum with bundled STEP CAD
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
- **Generalized meshes (type 13, CAD).** `Gen. mesh object origin: x/y/z` (= −Position in
  every object), `Number of vertices`, `Number of surfaces`, `Node positions:`
  (`index x y z`, cm, centered on the object's bounding-box centroid, 0-based) and
  `Surface element connectivity:` (`index v1 v2 v3 -1`, all triangles). The manual states
  that a GM object's position defaults to its bounding-box centroid
  (VR/grid/import_export_cad_files.html). Placement x = R·node + Position is consistent with
  every object (all angles are zero in the file, so R is untested on meshes). **H** for zero
  angles
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
  named by the `!` header of the file, `Ref. Coord. Sys.:  Target Chamber`, `D toward Port
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
  `visrad_unconverted` (2026-11-25) gives the red (1ω) and green (2ω) cone powers of
  each beam: the color column as TW under model 1, times the power column under model 2
  (fraction of the total) and model 3 (fraction of the blue cone, the power column then
  being the blue power), 0 for a cone that is off. The multiplier of each color and the
  random multiplier are applied as to the blue cone. Models 1 and 3 and the multipliers
  follow VisRad's dialog text and are not checked against a computed workspace. **L**
  (model 2 with the cones off, measured: **H**)
- **Deposition check (foil workspace, 2026-10-17).** The `[Power Source Data]` section
  lists, per object, one row per element: index, direct laser power source (TW/cm²),
  uniform source, albedo, x-ray conversion efficiency, laser reflectivity, pre-reflection
  power, numbers of red and blue cones hitting, and the minimum overlap angles at 1ω and
  3ω. The rectangle's elements run with its local x fastest (the transpose of
  `rectangle_mesh`'s order). Values are element averages, not centroid values (centroid
  values differ from it by up to 13% next to the peak). Against Briareus, four beams with a 620/360 µm HWHM
  n = 4.5 spot aimed at the 4 mm foil, tilted ~40° (9 × 9 elements of 0.44 mm):
  - Totals: 3.6000 TW in both.
  - Elements above 20 TW/cm² (peak 371): within 1.4%, center within 0.1%.
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
    (fitted Γ ≈ 10, c ≈ 0.94, M ≈ 7.7) takes the minor-axis edge from −21% to −2% and
    keeps the core within 1.4%. It leaves the diagonals at +8% and +14%. Three
    parameters fitted to ~20 informative elements: suggestive, not identified. **L**
- **Ruled out:**
  - Centroid or few-point element quadrature (vertices, 2×2/3×3 Gauss, Simpson,
    p × p grids).
  - Averaging the profile at the ray hit points.
  - Fixed 600-cell beamlets.
  - An elliptical envelope cut by itself (1.5–2.0).
  - Global or per-beam spot rotations, including quad-parallel square beams (±1.6–1.9°).
  - Axis or exponent changes. These fit the wings only by moving the core by 3–6%.
- **Power on the tilted foil.** VisRad's foil total is 3.5999999 TW for 3.6 TW
  delivered, exact to the table's rounding (~1e-6). Normalizing per transverse plane
  as the manual describes, with the hyperbolic semi-axes and the envelope
  (`VisRadProfile`), gives 3.60041 TW on this ~40° foil. VisRad therefore also
  renormalizes each beam's deposited total, perhaps by its ray-grid capture fraction
  (1.00000 for all four beams here). **M**
- **What would settle it.** The VisRad runs of the probe kits ("VisRad probe kits
  (2026-10-27): results"). The intersection-point export
  (VR/files/laser_target_intpts_file.html) lists a(z), b(z), n(z), f_SQ(z), ΔL_W(z), ΔL_H(z)
  and M(z) at each hit, which gives the defocus law directly. Rerunning with the foil normal
  to the beams, where the distance from focus barely varies, separates the focal profile
  from that law.

#### VisRad JSON workspaces (2026-11-16)

VisRad 21.0.0 saves a workspace it has computed as JSON, whatever format it was opened
from. Read off computed workspaces of a NIF CAD hohlraum and of Briareus-built targets.
Two of these had been saved by VisRad in both formats from the same run. Implemented in
`src/visrad_json.jl`.

- **Layout.** CRLF, three-space indent, strict JSON. `"Workspace Format ID": 1001`, then one
  object per text section: `Header data`, `Target component data`, `Target comp list
  structure`, `Laser beam data`, `Diagnostic port data`, `Chamber component data`,
  `Dialog settings data`, `Snapshots data`, `Movie maker data`, `Memorized views data`,
  `Power source data`, `Radiosity data`, `Capture clearance data`, `Target pos viewer
  data`, `SRF beam group data`. Repeated items are numbered keys (`"Target component
  element[1]"`, …, 1-based). **H**
- **Objects.** `"Object name"`, then the text block's keys with the same spelling (Length
  k value, Position k, Tilt/Rotation/Spin angle, albedo, the Include flags, …).
  Exceptions:
  - `"Object color"` is `[r, g, b]`, and `"Gen. mesh object origin"` is `[x, y, z]`.
  - Time-dependent tables are `…-Num rows` / `…-Num cols`.
  - Meshes are `"Object node array-Num elements"` plus `"Object node array[k]"` (`[x, y,
    z]`, cm, object frame), and four connectivity columns `"Object connectivity
    array-Values Col c"`, the fourth all −1.
  - Numbers inside arrays are strings; scalars are JSON numbers. **H**
- **Folders.** `"Target comp list structure"` lists `"Folder name"` paths and
  `"Target comp name"` as `"Target Components/<folders>/<name>"`. **H**
- **Beams.** `"Laser beam element[k]"` with readable keys, one per text-row column. VisRad's
  spelling is kept (`"Beam cone F numuber"`).
  - Port angles and f-number are unrounded (52.382267335385706), where the text row has
    52.382. Colors have 6 digits.
  - The cone half-angles are not stored. The text row's cone_max is atan(1/2F) to its
    2 decimals, and cone_min is 0.
  - All three power tables are listed with their row counts (0 in every file seen). Row
    layout of a non-empty table: *assumption*, as the connectivity columns.
  - On the 192 beams of a NIF workspace every column agrees with the text row to its
    printed digits. **H**
- **Power sources.** `"Num pwr src up to date"` (0 or 1), then `"Sim time element[1]"`
  with one `"Power source element[k]"` per target component, in the components' order,
  with no names. Each has `"Num surface elements"` and, per text column, `"<label>-Num
  elems"` and `"<label>-Values"` (strings, `%e`). On the paired files the values match the
  text table to its 5 decimals: 5 × 10⁻⁷ relative at 1000 TW/cm². **H**
- **CAD windings and VisRad's deposition.** VisRad deposits laser light on a CAD mesh
  regardless of how its triangles are wound. On a NIF CAD workspace the hohlraum parts
  are wound inconsistently: `orient_mesh` turns 48–52% of the area of the wall halves,
  LEH washers and window patches. Every triangle of the inner wall sheet faces away from
  the cavity, yet VisRad puts 405 TW on those 5,500 triangles, matching Briareus with
  `two_sided = true` or `orient = true` to 0.2%. Read one-sided as stored, Briareus puts
  2.5 of 422.5 TW on front sides and the rest on back sides. `visrad_target` now warns
  when more than 10% of a mesh's area is wound against its surface (without `orient` or
  `two_sided`). **H**

#### VisRad deposition parity: hohlraum runs (2026-10-21)

The profile model below (growth γ, then the tables `_VR_MEASURED`) is superseded by
"VisRad's beam profile law" (2026-10-27); the runs and their data stand.

VisRad 21 ran laser deposition on three workspaces built from Briareus geometry
(below: `n210808`, `n210808_reference` and `repointed`, a re-pointing of the reference): an
N210808-like hohlraum (Au wall r = 3.2 mm, 11.24 mm long, open LEHs; 2388 wall
triangles of ~0.6 mm, one VisRad element each), capsule, two window patches, foam, a
4 × 4 mm foil tilted 20° (3200 elements), sombreros; 192 NIF beams with the built-in
CPP spots, f/20.755, port radius 769.5 cm. Each file's `[Power Source Data]` holds the
per-element direct laser source. The
fits below used the Briareus operator points with a free-form evaluator.

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
  VisRad/Briareus 0.93–0.97 at the spot center, 1.1–1.7 at q = 1–1.3 (q = elliptical
  radius in 1/e units), i.e. VisRad's spots are wider on the wall (3–8 mm past focus).
- **Free-form profile fit.** With intensity P/(a(z)b(z)) g(q) and g a hat-basis
  function per cone (17 nodes, q ≤ 2.4), least squares over wall and window elements
  reproduces VisRad to 1.09% rms: VisRad's profile is a function of q only, per cone,
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
  0.52, 0.49 on repointed. Times the envelope ratio ρ: 0.69, 0.67, 0.66, 0.62
  (n210808) and 0.68, 0.67, 0.66, 0.63 (repointed), one constant γ ≈ 0.67. The form
  γ(ln 2)^(1/n)/ρ spreads more (0.65–0.75); z_A ∝ a₀² fails (≥ 5.8%). One γ for all
  beams on n210808: 0.60 4.2%, 0.64 3.2%, 0.68 2.80%, 0.72 3.0%, 0.76 3.5%. Adopted
  γ = 0.67 (`VisRadProfile(; growth)`). Every beam has the default ρ = (ln 1000)^(1/n),
  so a dependence on ρ cannot be told from one on n. What a(FL) ≈ ρ R_lens/γ ≈ 1.5 ρ
  R_lens means in VisRad's terms (d(z), the lens ratio 0.4) is unknown.
  **M** (inferred from VisRad output; one constant, fitted on one workspace)
- **Validation** (`visrad_profile(ws)`: per-beam envelopes, γ = 0.67, through
  `DepositionOperator`; before → after):

  | case | component | rms diff | peak, Briareus vs VisRad (TW/cm²) | correlation |
  |---|---|---|---|---|
  | n210808 (fit) | wall | 10.21% → 2.81% | 858.8 → 774.4 vs 778.1 | 0.9969 → 0.9997 |
  | n210808 | windows | 2.46% → 0.94% | 315.2 → 308.6 vs 308.5 | 0.9974 → 0.9995 |
  | n210808_reference | wall | 10.21% → 2.81% | as n210808 | |
  | repointed | wall | 10.66% → 2.68% | 942.1 → 897.4 vs 915.9 | 0.9969 → 0.9997 |
  | repointed | windows | 3.80% → 1.41% | 386.5 → 384.6 vs 374.6 | 0.9899 → 0.9989 |
  | repointed | foil (±2.5 mm of focus) | 2.81% → 2.60% | 1039.1 → 1039.2 vs 1032.0 | 0.99991 → 0.99993 |
  | foil workspace, 4 beams | elements > 20 TW/cm² | max 1.43% → 1.31% | | |

  Capsule-drive departures of repointed from n210808_reference (Briareus
  radiosity on each deposition; drive ratio, then δ1..4 in %): VisRad 0.9261, 0.006,
  0.092, 0.018, 0.046; before 0.9261, 0.006, 0.054, 0.012, 0.024; after 0.9261, 0.007,
  0.090, 0.019, 0.037. n210808 δ2: VisRad 2.576, before 2.574, after 2.573. δ2 moves by
  ∓0.008 for γ ± 0.05. **M**
- **What remains.**
  - Wall rms 2.7–2.8%. Outer-dominated elements carry 82% of it (4.5% rms; the outer
    peak ring at |z| = 2.7 mm is 3–4% low), inner-dominated ones 1.2%. Per-cone γ
    would give 1.1–2.0% per cone. δ4 stays 0.009 pp below VisRad's.
  - Totals on tilted surfaces: per-plane normalization puts 0.15% too much on the 20°
    foil (VisRad gives 20.80008 TW for 8 × 2.6, exact) and 0.06% on the 40° foil of the
    foil workspace
    (0.01% with the manual z_A). VisRad renormalizes each beam's total; not modeled.
  - Window patches in repointed get 1.0% more than in VisRad (n210808: 0.02%).
  - Near-focus wings on the foil (10⁻³–10⁻² of peak): ±10–40%, as before.
  - Occultation by the ray grid: VisRad puts ≤ 7e-4 TW/cm² on the LEH washers'
    elements, Briareus nothing.
- **VisRad beam-profile probe (2026-10-21).** A VisRad 21
  workspace built for the purpose:
  - **Setup.** One beam per NIF quad (48, 1 TW each, the template's CPP rows), each
    focused 60 mm out along its port and deposited on a separate flat plate normal to the
    beam, at |z| = 0–12 mm (inner cones) or 0–15 mm (outer) from best focus, before and
    after it. The plates have 2 × 44² triangles, laser deposition only.
  - **Fit.** Each plate is fitted with a rotated, shifted super-Gaussian
    averaged over each element (16 sub-triangles), cut at ρ. The fit leaves 0.0–2.5% rms;
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
  - **Validation through `DepositionOperator`:**
    - **Probe plates:** element rms 0.1–1.6% (inner), 0.2–4.5% (outer, largest at
      15 mm), against 2.6–20% with γ/ρ.
    - **n210808 wall:** 0.9% rms (from 2.8%), peak 779.1 against 778.1 TW/cm²,
      correlation 1.0000.
    - **repointed wall:** 0.9% rms, peak 911.5 against 915.9.
    - **Totals:** wall 0.07% below VisRad (γ/ρ: 1e-4); windows 0.17% below; the 20°
      foil 0.37% above.
    - **Capsule modes:** absolute P2/P4 agree with VisRad's deposition to
      0.014/0.004 pp (repointed).
    - **Departures of repointed from the reference:** δ1..4 0.008, 0.104, 0.020,
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
  are plausible (for an LEH washer 6.2 mm across and 50 µm thick, 84 → 1.1 mm³), and the
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
  0.9–1.1 elements per size²; `rmax = size` (in place of 0.75 size) avoids splitting
  Lloyd cells into fragments (0.75 doubled the count). **H** (measured locally; the
  constants are engineering choices, **M**)
- **Cavity selection (`enclosure`, 2026-10-17).** On the same workspace (oriented,
  `coarsen = 0.3`: 155k elements, 212k facets), 64 diffuse rays per element and the
  gathered flux take 3.0–3.5 s on 8 threads (7.8 s on 3; the kept set and fluxes are
  identical). From the capsule center the kept set is only the capsule and the inside of
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
  weighted by area (the fill tube is ~7% of the area): center at TCC to 1e-4 mm,
  r = 1.1333 mm, the outer face of the shell (the median vertex radius within 1.5 mm
  is 1.1333). Cavity + 3,120-element capsule (h = 0.09): 9,726 elements, scenario
  built in 93 s on 8 threads, 3.1 GB live; an `optimize_pointing` evaluation ≈ 0.4 s.
  The x-ray view factors are refused above the element count whose three dense N×N
  Float64 matrices (view factors and two cached factorizations) fill half the memory
  (`radiation_size`: 26.8k elements on 32 GB, *derived*); the drive is then laser
  only.
- **The window's cavity and stalls (BriareusGUI, 2026-10-18, measured locally, H).** Cavity
  selection on the window's remeshed elements against the library's faceted two-pass
  (`visrad_target(...; remesh = 0.28, inside)`: select on the oriented CAD facets, remesh
  the kept facets, select again), same objects, tent left out: one selection 0.8 s (2.2 s
  with the pocket search) and 7,155 elements, 360 mm², against 5.3 s and 6,908 elements, 339
  mm². The window keeps each part remeshed whole (it draws it, deposits laser on it and
  checks clearance against it), so the faceted path would add the remeshing of the kept
  facets to every cavity build; the 6% extra area is whole elements straddling the cavity's
  edge (mostly a cap washer, 413 vs 295). For element-size edits, building the `Target` of
  the 12 parts and of the 200k-facet ghost took 1.2–1.5 s per edit (BVHs); the BVHs are now
  reused for identical components; a CAD remesh runs on workers with the part drawn from its
  facets meanwhile. The ~1 s stalls at the first drive result, the first coloring by a field
  and the first ghost were JIT compilation (Makie's heatmap colormapping, the colorbar's
  first drawing, the transparent mesh's first upload), not GC (≤ 0.26 s pauses at 3.3 GB
  live); the window now compiles them before it opens.
- **The window's target builds, capsule fit and drive end (BriareusGUI, 2026-10-18, measured
  locally in the scripted Qt window on a NIF CAD workspace, 8 threads, H).** `Target` of the
  12 calculation parts (275k facets, 20k elements) 0.39–0.40 s, of the 4 drawn-only parts
  (141k facets) 0.18–0.19 s; the components themselves 1–5 ms (cached parts). The CAD
  capsule's sphere fit 0.13 s (0.79 s at the first call, including compilation).
  `fit_surfaces` on that capsule (38k facets, sphere and fill tube) took 133 s
  single-threaded, so the window's whole-part curved fit first checks that every vertex lies
  on one algebraic sphere or cylinder (ms): none of the workspace's 14 CAD parts is one
  surface as a whole (a hohlraum half has a bore, walls and ends), so the fit is neutral
  there; a sliver tube (the synthetic workspace) gets its cylinder. The ~1 s GUI-thread gap
  at a drive's end was inside the window's refresh, not Qt: compilation of the drive's
  result path (selection, drive and tree refreshes: 1.0 s and 0.9 s with 1.3–1.5 s compiled,
  the evaluation compiling on a worker meanwhile) and a full reset of the model tree's rows
  (0.2–0.26 s, every delegate rebuilt) at every row-count change. A precompiled drive (the
  package's workload) and row inserts/removals where the tree changes bring the largest gap
  at the drive's end from 1.2 s to 0.31 s, of which 0.23 s is a GC pause.
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
    plus the facet deviation of the region of), 1.3 s: 74% of the area, the inner wall
    881 elements and 84.6 of its 86.7 mm² (elements straddling the wall's edges stay
    flat). **H** (measured locally)
- **Curved quadrature on remeshed elements (2026-10-18).** The cylinder/sphere quadrature
  (sub-triangle points projected radially, exact curved areas via `_curved_area`) is exact
  in area for any element whose facets are a radial graph over the surface, not only chord
  rectangles. Its cost was too high: per-facet points made a 12-facet element 6× (36×
  per pair) dearer than a primitive one. Elements of more than two facets now lump the
  facets' pieces into cells of size/(2^level √2), as many as a two-facet element has.
  Remeshed elements also need their centroid on their facets: on a convex surface (a CAD
  capsule) the area centroid of an element folded across facets lies inside the solid, and
  its visibility rays were blocked by the element itself (every capsule element saw
  nothing). **H** (synthetic, see below)
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
  separately (16 sub-triangles each, then adaptive). Now, (1) the receiver is the
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
  1.2 → 1.1 s (1.0 s for primitives at 0.5 mm, 2,244 elements); the synthetic
  thick-walled hohlraum of test/cadfit.jl at 0.5 mm (1,217 elements) flat 47.1 → 4.1 s,
  fitted 2.3 → 1.3 s; integrated results unchanged (capsule flux, a₂, a₄, LEH loss to
  ≤ 2e-6 relative or 1e-4 points; capsule rows identical to 4 digits). A flat CAD pair still costs ~5× a primitive pair: folded
  outlines have ~3× a quad's edges, and the cells' bookkeeping. Primitive targets
  (1–2 facets per element) are unchanged. **H** (synthetic)
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
  cone, are unchanged (seven primitive targets; the tutorial's half, washer,
  capsule and tent cup fitted raw and remeshed). **H** (synthetic)
- **Cones in the NIF CAD workspace (2026-10-18, local).** Raw facets, one thread; time
  with cones vs without: hohlraum halves 1.1/0.8 s vs 0.7/0.6 s. In each half, the inner
  end wall beside the LEH is a 84.47° cone facing the cavity (apex 5.789 mm from the
  waist plane, radius 2.153–2.363 mm, 3.04 mm², rms 0.25 µm; it meets the r = 2.15 mm
  LEH bore), and a 61.1° chamfer at the waist edge of the inner wall (r 3.21–3.25 mm,
  0.19 mm², axis within 1.7e-3 rad of z). The LEH washers have a 77.14° chamfer on the outer
  rim (r 3.065–3.082 mm, 0.35 mm², facing out). A waist sleeve has a 10.0° bore chamfer
  facing in (2.35 mm²) and two 45.0° chamfers facing out (1.72 mm² each), rms 2e-10 mm.
  Cones also absorb pieces of free-form faces (fillets) within tolerance: on the halves
  8–10 more cones of 0.001–0.12 mm², most with axes tilted 0.03–1.5 rad, rms 0.5–1 µm
  (previously unfitted or spheres). Two halves through `visrad_target(remesh = 0.3,
  fit = true)`: 116 components (102 cylinders, 12 spheres, 2 rest) before, 138 with 22
  cones; curved area 353.6 → 360.8 of 490 mm²; view factors no slower (10.8 s before,
  7.9 s with cones, 8 threads) and identical on the inner wall (self-view, wall-to-wall,
  6 digits). **H** (measured locally)
- **Near-flat and small fitted regions (2026-10-18, local).** `curved_components` now takes
  `maxradius` (relative to the part's bounding-box diagonal, 11 mm for a half) and
  `minfraction` (relative to the part's area, 245 mm²); both off by default. Two halves,
  remeshed at 0.3 mm, with cones and only `maxradius` set, maxradius 1, 2, 5 give 127, 128,
  131 components (the R = 23–670 mm spheres and cylinders are only ~10 of them; most of the
  138 are small holes and edge faces); minfraction 1e-3, 3e-3, 1e-2 → 25, 16, 10 components
  with 356.6, 353.2, 344.1 mm² curved (of 360.8); minfraction 1e-2 keeps the inner walls,
  outer sheets, steps and the two 84.5° cones. **H** (measured locally)
- **Tori (2026-10-24).** `fit_surfaces` fits tori (`:torus`, in the default kinds) and
  `TorusSurface` carries them into the radiation model. The fit takes candidate axes from the
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
  - neighboring faces: within the relative tolerance (rtol min(r, L), 5 µm for r = 10
    mm) a large torus (R 13.5, r 10 mm) took a NIF part's cone and 40% of a cylinder
    beside it at rms 16 nm, and a fillet torus 9 boundary strips of the cylinder it
    blends into;
  - band tessellations: strips with vertices only on their rims lie on any circle through
    three rims, and four rims placed symmetrically (two cone bands at either end of a
    cylinder) lie on one circle too: with tolerance `atol` the large torus above came back
    through such rims, at round-off, over the cone and the cylinder's strips. So tori are
    held to CAD round-off: tolerance `atol`, rms under `atol/16` (6 nm; the NIF CAD tori fit
    at 0.1–5 nm), at least four evenly spaced points along the tube circle (largest step
    under 3× the median past the empty arc), no plane, cylinder, cone or sphere fitting the
    region, and merges that keep that rms. A torus within 2 tol of a sphere (R), cylinder
    (L²/8R) or cone/plane (L²/8r) is left to it, and a torus growth that becomes one, or
    that a simpler surface fits within round-off, is given up at its checkpoint; the patches
    of a rejected torus region (half its radius) grow no torus again, and a torus starts
    from at most 64 patches (rings around CAD fan apexes hold hundreds). Without these the
    cylinder–sphere–torus fit of the NIF capsule took 103–123 s (torus growths retracing the
    sphere and the fill tube's fans); with them 21 s against 18 s. A boxed closure in the LM
    step (`ntuple` over a reassigned matrix) had made each torus fit 10× slower. Synthetic
    (`test/cadfit.jl`): tilted and moved quarter-to-150° sections of R = 2, r = 0.5 mm, fine
    (12.5k facets) and coarse (19° bands), through Float32 STL, facing out or in, raw and
    remeshed with `reference`: one torus, R, r, center to 1e-6, axis to 1e-10. Fine fillet
    (all kinds): 7.9k / 32.8k / 65.5k facets in 0.13 / 0.51 / 1.6 s, one region (without
    tori 4 s at 8k and ~80 s at 33k, quadratic; "Cost of fitting faces of none of the
    kinds", below). `TorusSurface` projection (rim-preserving since 2026-10-24): in the
    point's meridian half-plane, inside the tube circle the nearest point of the circle
    (radially from the spine point), on or outside it the circle's point at the same height
    on the point's side of the spine (radially from the axis, as on a cone; the pole when
    above or below the tube). The two agree on the circle. A rim (circle ψ = ψ0 about the
    axis) has its chords at its height, displaced toward the axis by their sag δ: outside
    the tube for ψ0 on the poles and the inner half (distance² r² + 2δr|cos ψ0| + δ²), so
    they map onto the rim; inside it on the outer half, where only the outer equator (ψ0 =
    0, the horizontal is the radius) maps exactly. Facet points outside the tube lie within
    the azimuthal sag of it horizontally (the meridional sag points inward), so the
    horizontal map moves them by at most that sag: second order like the nearest point (max
    distance from the nearest point on a quarter fillet 3.8e-2 / 9.6e-3 / 2.4e-3 / 6.0e-4 mm
    at 16 / 32 / 64 / 128 strips, ratio 4.00). Exact rims everywhere are not reachable by a
    projection of the points only: exactness on every outer-half rim forces the horizontal
    map throughout the inside of the outer half, whose tangential error ε tan ψ (depth ε) is
    first order next to the poles (on a quarter fillet's cap strip, ~0.2 r Δψ). A
    facet-aware map (barycentric in (φ, ψ)) would be exact everywhere but changes every
    `surface_point` caller. Patch areas r |∮ (Rψ + r sin ψ) dφ| by Gauss–Legendre on edge
    pieces, split by bisection where an edge crosses the tube circle (the integrand kinks
    there); a fine subdivision of the projected facet agrees within its discretization
    error. These cases are covered within round-off (|error| ≤ 7e-16): closed torus, ψ ∈ [0,
    π], [−π/2, π/2], [π/2, 3π/2], and the quarter fillet at 16 / 32 / 64 strips (before:
    +7.0e-3 / 1.8e-3 at 32 / 64, the cap-rim chords projecting past the rim); ψ ∈ [0, π/3]
    (rim at 60°) +4.6e-2 / 1.1e-2 / 2.8e-3 (second order). The rounded hohlraum (rc = 1 mm)
    with its corners as tori: fillet area exact (the frusta's chord polygon −15% / −3.8% /
    −1.0% with 1 / 2 / 4 bands per corner), raw row sums max |1 − ΣF| 5.6e-3 / 3.0e-3 with 1
    / 4 bands (frusta 7.1e-3 / 3.0e-3), on the fillets 6.0e-4 / 1.2e-3 (frusta 4.7e-3 /
    2.2e-3; with the earlier projection 4.5e-5 / 1.8e-3, the one-band value a cancellation
    of the cap rims' +1.8e-3 area against the quadrature). Exchange against the frusta (Σ
    A_i F_ij): capsule–fillets −2.9e-3 / +1.5e-3 / +7.1e-4, fillets–LEHs −9.6e-3 / −2.9e-3 /
    −1.5e-3 with 1 / 2 / 4 bands; `hohlraum` keeps its frusta (the default 4 bands at h =
    0.5; baselines). NIF CAD (14 parts other than the mount): the fits are unchanged
    (fitting uses a model separate from the projection); the 16 tori's 11.2k elements 177.64
    → 177.74 mm² (+6.0e-4): the full quarter fillets (r 0.9 mm) from +2.1e-3 of their band
    to −2e-7, inner-half fillets (r 1.5 mm) from −2.4e-3 to ≤ 1e-4, thin tori (r 25 µm,
    inner-half rims, azimuthal facets long against r) +8–13% (a full one from −11% of its
    band to −1.3e-5), the convex round between cones (rims −29°, +23°) unchanged to 1e-15,
    and a fillet (r 1.05 mm) whose region takes a 0.4° sliver of cap past the pole −0.30%:
    facets lying nearly in the pole's tangent plane are squeezed by the horizontal map (one
    1e-3 mm² facet −91%), though the band they tile is kept. Tori are not convex: their
    facets still occlude pairs inside them. NIF CAD workspace, its 14 parts other than the
    mount, raw facets, default kinds with and without `:torus` (and cylinder–cone–sphere,
    cylinder–sphere): 9 parts are unchanged; on the other 5 every region that changes is a
    fillet patchwork replaced whole, 29 small cones, cylinders and spheres (0.11–5.8 mm² per
    part) by 16 tori (R 1.6–3.4 mm, r 0.025–1.5 mm, 0.31–61 mm² per part, which also take
    the fillet facets that stayed unfitted), and on one part two spheres refitted (same
    patches, centers 1e-16 mm apart). The capsule (fill tube) is unchanged with every kind
    set, no tori: with the default kinds its tapered tube is two cones (0.147°, 0.110°) at
    round-off, with cylinders and spheres the 25 regions of "Cost of fitting faces of none
    of the kinds" (tolerance-stopped cylinder pieces). Fit time for the 14 parts 9.0 → 4.3 s
    (fillets no longer seed rejected growths); cylinder–cone–sphere 8.3 → 5.0 s;
    cylinder–sphere 21.6 → 25 s (capsule 18 → 21 s). The 133k-facet mount (default kinds,
    one thread): 1911 → 1948 regions, 1278 → 1100 s; 85 tori (347 mm², R 0.02–5.9 mm, r
    0.02–1 mm, rms ≤ 6 nm) take 54 of the 61 regions that change (39 cones, 19 cylinders, 3
    spheres; 7.4 mm²) whole, and every other region that changes or appears shares vertices
    with a torus: a 106 mm² cylinder (R 3.58 mm) that stayed unfitted before, and 0.3 mm² of
    small cylinders, cones and a sphere re-cut or newly fitted beside fillets. **H**
    (measured locally and synthetic)
- **Cost of `fit_surfaces` on finely faceted CAD (2026-10-18).** Not the sphere: a
  finely faceted sphere is one region grown from its first seed (lat-long spheres of
  10k and 40k facets, an icosphere of 20k: 0.01–0.07 s). The NIF CAD capsule part
  (38,200 facets, raw, one thread, local) has a fill tube tessellated in fans: vertex
  degrees up to 801, so the patch-neighbor lists (patches sharing a vertex) held 4.4M
  entries (116 per patch), and each seed's two-ring start walked 1,600 lists of 1,600:
  83% of the time. Neighbors are now visited through each vertex once per pass, with
  the same order of first occurrence, so results are unchanged (every region
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
  faceting and grows on. The growth steps themselves are
  unchanged, so a growth that is not given up ends exactly as before. Synthetic
  (`test/cadfit.jl` shapes): the 0.15° tapered tube (r = 50 µm, 2 mm) with cylinders and
  spheres, 2k / 8k / 18k facets 1.8 / 27 / 171 s before, after 1.5 s (8k) and
  3.3 / 6.5 / 30 s at 18k / 37k / 131k; patches tested 0.99M / 1.93M / 3.80M at 9k / 18k /
  37k facets (×1.96 per doubling; 41.8M at 9k without covers). A quarter-torus fillet
  (R = 2, r = 0.5 mm, all kinds) is unchanged, 4 s at 8k and ~80 s at 33k facets: its
  rejected growths stay under 512 patches (~200–500), so it is quadratic until they
  pass it (a torus kind would end that). NIF CAD workspace, its 14 parts other than the
  133k-facet mount, raw facets, kinds default / cylinder–cone–sphere / cylinder–sphere
  (42 fits, shared machine): 41 are unchanged (kind, surface, centroid, normal,
  patches, area, rms, max deviation, facing of every region); 8.3 → 8.6, 8.1 → 8.3,
  59 → 22 s (56 → 18 s for the capsule only); the 133k-facet mount with the default kinds
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
- **Cost of fitting large CAD parts (2026-10-25).** Profiled `fit_surfaces` (default kinds,
  one thread) on the 133k-facet NIF mount and on samples of it (the N patches nearest its
  centroid). Before the change, 40k facets took 26–30 s, 73% of it the merge of regions of
  one kind (every pair rescanned after each merge, each pair recomputing both regions'
  extents: quadratic in the regions times the merges). The merge now tries each pair of
  region versions once (a pair's outcome depends on its two regions only) and keeps the
  tolerances, and the final fit of a growth reuses the last fit when the members are
  unchanged. Results are unchanged (kind, surface, centroid, normal, patches, area, rms, max
  deviation, facing) on 90 fits: the 14 other NIF parts with four kind sets, the test shapes
  and a synthetic boss plate, and the mount samples; the full mount identical to the
  2026-10-24 run (1948 regions). Mount 1100 s → 29 s; samples 20k / 40k / 80k facets 7.4 /
  26 → 2.9 / 6.0 / 17 s (now linear). Of the remaining time 53% is the start fits on the
  seed rings (tori 29%, cones 13%), 37% refits while growing, 7% the merge. Comparisons must
  use one compile mode: the same code gives tori differing by ~1e-13 compiled in the package
  image and compiled in the session (likely FMA). Variable-radius fillet (spine ρ = 2 + r,
  tube r = 0.5(1 + 0.3 cos φ), all kinds): 16k / 26k / 37k / 50k facets 16 / 38 / 71 / 107
  s; cylinder and cone growths end tolerance-stopped (rms 0.41 of tol) at 200–370 patches,
  under the 512-patch covers, so the cost is quadratic in that window and drops to ~N^1.3
  once growths pass 512 (8,994 given up at 50k). No exact bound gives such a growth up early
  (the tolerance grows with the region's extent, so its rms in tolerance units can fall);
  smaller covers changed NIF parts before ("Cost of fitting faces of none of the kinds"), so
  none was adopted. **H** (measured locally and synthetic)
- **Fitting CAD by default (2026-10-24, local).** Should `fit = true` follow `remesh` in
  `cad_target`/`visrad_target`/the GUI (ROADMAP open question 1)? The NIF CAD workspace, its
  11 parts near the cavity (no mount, far parts or tent), `visrad_target(...; orient = true,
  remesh = 0.6 or 0.3 (capsule half that), inside = (0, 0, 0))`, 4 threads with other jobs
  running. As for the fit cost, `fit_surfaces` on the parts' raw facets 7 s in all (halves
  1.9–2.7 s; the 38k-facet capsule with its fill tube 0.29 s with the default kinds: the
  quadratic case of the previous item needs a face of none of the kinds), against ~220 s
  (0.6) and ~420 s (0.3) for the cavity selection and remeshing. Components: 10 parts → 158
  (0.6) / 159 (0.3) fitted, 68 with `minfraction = 1e-3` (40 cylinders, 9 cones, 9 spheres,
  9 rests, one part unsplit); the capsule's mean flux and moments with 1e-3 are within 3e-4
  and 0.04 points of fitting everything (`:sampled` runs, below), so 1e-3 became
  `curved_components`' default. Left in with `maxradius` off: cylinders of R 54–670 mm and
  spheres of R 0.94–6.4 mm besides the capsule's (all of which `SphereKeepOut` would take by
  default). View factors, default `:sampled` visibility: 47 → 49 s at 0.6, 177 → 161 s at
  0.3; with `:centroid` visibility at 0.6, 33 → 22 s. The 22× of the sliver tube (above) is
  a whole-part fit, which the GUI already makes (`_whole_surface`). Accuracy against
  `montecarlo_incident` on the flat facets (α = 0.7, a source on the inner walls, 4,000
  paths per receiver), `:centroid` visibility: outer capsule flux flat +0.03%, fitted +0.26%
  (MC `:element` 0.5096 ± 0.0003; the fitted capsule is the true sphere, and `:hit` gives
  +1.3% over `:element`); halves (400 elements) flat −0.81%, fitted −0.22%. **With the
  default `:sampled` visibility the outer capsule's row sums are 0.21 (0.0–0.99) against
  0.955 with `:centroid`, and its flux 7× too low, flat or fitted; the halves' rows 0.79
  against 0.91 (flux −15%).** `_visibility_samples` averages each grid cell's sub-triangle
  points, and on an element folded across the facets of a convex part (13.5 facets per
  capsule element) that mean lies behind the facets, where every ray is blocked (the failure
  `Target` avoids for the centroid by moving it onto the facets). Fix before relying on
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
  remaining blocked sample ray of the reduced capsule rows (561k of 561k) hit facets
  of the emitting element, rays toward the horizon of the facet under the sample that
  the element model (radiating from the mean plane) already counts. So an element of
  more than two facets no longer occludes the sample rays it emits (nor a receiving
  element those it receives). On primitives `:sampled` is unchanged (`==` on the
  bored cavity at 0.5 mm and `hohlraum(h = 0.8, capsule_h = 0.3)`). Raw factors
  (`closed = false`), a two-ring source on the driven walls, against `montecarlo_incident`
  (`:element`, 8,000 paths per receiver):

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
  the change in its capsule rows alone adds +0.04% at fixed radiosity, so the shift comes
  from the walls' rows. The driven inner wall (ρ 3.0–3.4, 3,000 paths, σ 0.07–0.08%) shows
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
  the `SphereKeepOut` default only the capsule part (9 spherical components),
  `show` "10 parts (68 components)". So fitting costs nothing measurable here and cuts
  the flux error 5–10×. **H** (measured locally)

### Result tables and fits

- **Surface Element Properties.**
  - Saved via File | Save Data for Object / All Objects / Selected Rows as
    "multi-column ascii". `*.sed` files from a run have the same contents
    (VR/results/display_surf_elem_table.html, VR/simulation_times/config_sim_params.html).
    **H**
  - Positions are in the coordinates of the component, (r, θ, φ) for spheres. **H**
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
  - The fit weights and the Condon–Shortley phase are not stated. `visrad_legendre_fit`
    and `visrad_ylm_fit` are unweighted and have no phase. **L**
- In VisRad, the grid size of Sphere, Cylinder, Cylindrical Hohlraum and similar
  objects is the outer wall; the inner wall lies inward by the Wall Thickness
  (VR/intro/revision_history.html, v21.0.0). **H** The benchmark kit's README states so.

### VisRad probe kits (2026-10-27): results

VisRad 21 ran the probe workspaces of 2026-10-23 and 2026-10-27, each built from
Briareus geometry with NIF beams at f/20.755, port radius 769.5 cm, and 1 TW 3ω per
beam. The workspaces, `.bti` exports, report export and readback scripts are kept
outside the repository with the other VisRad runs. Settings recorded by the user: AMR
segment length/d(z) 0.5 and distance/d(z) 1.5, beam cone display 24 × 1, and "write beam
spatial profile parameters" on.

#### VisRad's beam profile law

- **The `.bti` export (Format ID 3, profile data written).** Per beam with power:
  - index, name, port polar and azimuthal angles, intersection point in the object's
    frame, object type, and the point in chamber spherical coordinates (cm, degrees);
  - then a(z), b(z) (µm), n(z), ΔL_x(z), ΔL_y(z) (µm, always equal), a seventh column
    equal to a(z) in cm (to rounding; labeled m(z)), and f_SQ.
  **H** (four files, 299 rows)
- **Closed forms.** With z the hit's distance from best focus along the beam, focal 1/e
  semi-axes a₀ ≥ b₀, exponent n₀, f-number F and envelope ratio ρ, z_A = √2 a₀ F and
  z_B = √2 b₀ F, VisRad's exported parameters are:
  - a(z) = a₀ √(1 + (z/z_A)²), i.e. a(z)² = a₀² + z²/(2F²), and likewise b(z);
  - n(z) = n₀ (1 + 0.2 |z|/z_A), the manual's form with z₀ = z_A;
  - f_SQ = |z| / √(z² + z_A z_B), the manual's form with Γ = 1 and "a(z) b(z)" read as
    z_A z_B;
  - ΔL(z) = ρ a₀/√2 + 0.023995 |z| mm, a slope of 0.99605/(2F).

  Over 299 hits (four files: normal plates at z = −8…15 mm, tilted plates, envelope
  plates at ρ = 1.2 and the default, and a hohlraum's wall hits) the forms hold to
  5 × 10⁻⁵ (a, b), 1.4 × 10⁻⁴ (n), 5 × 10⁻⁴ absolute (f_SQ) and 1.3 × 10⁻³ (ΔL). The ΔL
  slope is constant to 0.04% across plates, cones and ρ. **H** (fit to the parameters
  VisRad exports). No dependence on the final-optic distance or on F other than through
  these forms was seen (every run had 769.5 cm and f/20.755). **L** for other F.
- **The profile.** The manual's blend, I = I₀(z) exp(−α) with
  α = f_SQ [(|x|/ΔL)^M + (|y|/ΔL)^M] + (1 − f_SQ) q^(n/2), x, y on the spot axes,
  evaluated with these parameters at the z of each point, leaves one unknown, M(z).
  - Fitted per plate on the 48 normal plates of the 2026-10-21 profile probe, M rises
    linearly in |z| with a slope that differs by cone. Every cone fits one form,
    M = M₀ (1 + |z|/√(z_A z_B)), with M₀ ≈ 9.9: the slope times √(z_A z_B) is 10.02,
    9.98, 9.90, 9.87 for the 50°, 44.5°, 30° and 23.5° plates. A free scale on ΔL comes
    back 1.001–1.004.
  - With M₀ = 10 fixed (a scan over 9–11 has its minimum there) and no per-plate
    parameter, the 48 plates are within 0.05–0.22% rms per element (elements above 1%
    of peak). The tabulated fit it replaces (`_VR_MEASURED`) left 0.5–4.5% on the same
    plates.
  - On the default-ρ plates an elliptical cut at q = ρ leaves 1.7–2.0% at
    15 mm, and a cut on α at ln 1000 (10⁻³ of the local peak) or none leaves 0.13–0.17%.
    On ρ = 1.2 plates (z = 0 and 10 mm) a cut at q = ρ leaves 4–6%; the fit improves up
    to a cut at ≈ 1.2ρ and is flat beyond, so ρ = 1.2 cuts nothing visible above
    ~0.5% of peak. Adopted: cut at α = ln 1000 for every ρ, so ρ enters ΔL only.
    **M** (the cut at ρ < default is constrained only to lie beyond ~1.45 in q).
  - I₀(z) normalizes each transverse plane (numerically). The plate totals are within
    10⁻⁴ of VisRad's on normal and tilted plates.
  **H** for M and the form on normal plates (fitted), **H** as checked below.
- **Checks on data not used in the fit** (absolute: Briareus's normalization, not
  rescaled; rms over elements above 1% of VisRad's peak):

  | Surfaces | Elements | Law | Tables it replaced |
  |---|---|---|---|
  | 32 plates at z = 4–5 mm, cos i = 1–0.6, tilted toward a or b | | 0.12–0.22% | 1.26–1.46% |
  | 23 normal plates, default ρ, z = 0–15 mm | | 0.10–0.27% | 0.24–4.6% |
  | normal plates at ρ = 1.2 (z = 0, 0, 10 mm) | | 0.06, 0.21, 0.14% | 5.1, 6.3, 10.7% |
  | N210808-like hohlraum wall (2388 el.) / windows | | 0.62% / 0.20% | 0.68% / 0.61% |
  | re-pointed hohlraum wall / foil (3200 el.) | | 0.69% / 0.55% | 0.77% / 0.54% |
  | N230729-like fine wall (38 272 el.) | | 0.45% | 0.96% |

  The hohlraum rows are Briareus's `DepositionOperator` at default resolution; the plate
  rows evaluate the law on 16 sub-triangles per element. **H**
- **Envelope ratio.** VisRad 21 replaced ratios of 2.0 and 3.0, written into the
  workspace, with the default (ln 1000)^(1/n₀) (1.53645 and 1.49576 in the saved file);
  1.2 was kept. `VisRadProfile` caps ρ at the default likewise. The lens radius ratio
  (0.2, 0.4, 0.8) does not change any profile parameter by more than 0.03%. **H**
- **Implemented** in `VisRadProfile` (2026-10-27), replacing the tables `_VR_MEASURED`,
  the hyperbolic semi-axes with growth γ = 0.67, and the `focal_length`, `growth` and
  `measured` options. The normalization ∫∫ exp(−α) is tabulated per spot in |z| (0.1 mm
  steps to 60 mm; Gauss–Legendre in scaled polar coordinates out to the cut-off along
  each direction, within 2 × 10⁻⁶ of a brute-force sum) and interpolated linearly. Its
  derivative with respect to the spot size enters gradients only through a(z) b(z).
- **A tapered cut for optimization (2026-11-18), `VisRadProfile(; taper)`.** The hard
  cut makes deposition discontinuous in the aim: an operator point crossing α = ln 1000
  adds or removes exactly 10⁻³ of the beam's local peak times its weight and |n·k|.
  On a public mock-up (the default `hohlraum()`, 192 `visrad_nif_spot` beams in the
  quad frame aimed at their entrance-hole centers, one 44.5° quad's aim moved along z in
  10⁻⁵ mm steps over 2 × 10⁻³ mm), 13 points crossed; each jump of the wall power
  matched the prediction from the crossing point to 10⁻⁴ (the ratio to W |n·k| I₀ was
  0.001000 every time), 5 × 10⁻⁶ TW or 6 × 10⁻⁷ of the quad's wall power, and 2.9 × 10⁻⁶
  of the wall's P2 moment. A central difference whose stencil straddles a jump is off
  from ForwardDiff by the jump over 2h: 0.36% of the P2 moment's derivative at h =
  10⁻⁴ mm, 36% at 10⁻⁶ mm. A constrained re-pointing run on a saved study case
  (scratch, not in the repository) stalled in its first inner run, where the augmented
  Lagrangian decreased along its descent direction only for steps ≤ 10⁻⁶ mm. **H**
  - The option multiplies exp(−α) by the quintic smoothstep of t = (ln 1000 − α)/taper,
    C² at both ends of the band, and renormalizes I₀(z) with the tapered integral (the
    normalization of the law is the integral inside the cut, so the taper's share is
    restored). The support, and with it the operator's footprint and the dual screen,
    is unchanged. Placing the band below the cut in α is where the law is least
    constrained: a cut at ln 1000 and no cut fit VisRad's plates equally well (above).
  - Its cost with taper = 1 (intensities between 10⁻³ and e × 10⁻³ of the local peak):
    the share of a plane's integral in the band is 4–12 × 10⁻⁵ (NIF plates, |z| ≤ 10
    mm), so the renormalized interior rises by that much; on the mock-up the absorbed
    total changes by +6 × 10⁻⁶ (redistributed from clipped light), wall elements by at
    most 6.9 × 10⁻⁵ of the wall peak (rms 2.5 × 10⁻⁵), bright elements (above 1% of
    the peak) by at most 0.24% of their own value, and the wall peak by 5 × 10⁻⁵.
    Against the VisRad deposition on the 2026-10-23 probe plates (16 points per
    element, rms over elements above 1% of VisRad's peak), the per-plate rms moves by
    −0.011…+0.010 percentage points: 0.111–0.181% → 0.112–0.189% on the 32 tilted
    plates, 0.041–0.460% → 0.030–0.458% on the 27 envelope plates; plate totals move
    by under 10⁻⁵; no element changes by more than 8.6 × 10⁻⁴ of its plate's peak. With
    taper = ln 2 the changes are about 60% of these. (The hard-cut statistics are recomputed
    here: the tilted plates agree with the row above, 0.11–0.18% against 0.12–0.22%,
    while some envelope plates, at 0.42–0.46%, exceed the normal-plate row's 0.27%;
    the difference was not traced.) **H**
  - What it buys, on the mock-up: no jumps (second differences of the wall power at
    10⁻⁵ mm steps 10⁻⁴ of the hard cut's jump and below), and ForwardDiff agrees with
    central differences as h²: 9.4 × 10⁻⁷, 1.5 × 10⁻⁸, 7.8 × 10⁻¹¹ (P2 moment) at h =
    10⁻³, 10⁻⁴, 10⁻⁵ mm across the former jump. The normalization table's linear
    interpolation in |z| (0.1 mm nodes) leaves kinks in the derivative, which did not
    show at these steps. On the stalled re-pointing run (scratch; 92 parameters, a wall
    ceiling and four mode gates, CBET on), started the same way with taper = 1: the
    gates' finite differences along a random direction agree with ForwardDiff to
    10⁻⁷ at h = 10⁻³ and 10⁻⁴ mm (2 × 10⁻³ to 8 × 10⁻³ with the hard cut); the first
    inner run ended after about 540 evaluations instead of stalling (1350 evaluations
    on its last 33 iterations with the hard cut), and in the same 7-minute limit the
    run took 8 outer iterations (1404 iterations, 1.24 evaluations each) against 1,
    bringing the gate ratios from 1.3–5.7 to 1.004–1.020. The case stays infeasible (its
    wall ceiling), and at the hard cut's stalled design the line search along the
    projected gradient still fails above 10⁻⁶ mm with the taper: what remains there is
    not smooth in the objective's keep-out penalty (its drive term agrees; the objective's
    ForwardDiff −1.229 against central differences −1.169, −1.193, −1.207 at h = 10⁻³,
    10⁻⁴, 10⁻⁵ mm), not in the deposition. **H**
  - The default stays VisRad's hard cut (the law reproduces VisRad); the taper is an
    engineering option for optimization. **M** for the choice of width (one e-fold,
    not tuned). Since 2026-11-19 `optimize_pointing` uses it on a hard-cut
    scenario by default and reports under the hard cut ("Defaults changed (2026-11-19)"
    §4); the default of the profile is unchanged.

#### Tilted plates and per-beam totals

- **No profile frozen at the axis hit.** On 32 tilted plates the near and far halves of
  each plate agree with Briareus within 1% in the same direction (0.986–1.001); VisRad
  evaluates the profile at the z of each element. The 44.5°/50° wall-shoulder deficit
  of the old tables (0.84–0.97 on n230729_fine) was the tables' a(z), not incidence:
  the shoulders are flat in cos i (0.995–0.997). **H**
- **Per-beam renormalization.** Every fully captured beam's total is its power to 10⁻⁵,
  at 0–78° incidence (steep foils: per-plane normalization gives +0.23, +0.55, +0.99,
  +2.2% at 40, 60, 70, 78°). Not implemented. **H**
- **Partly lost beams** (a foil alone in part of the spot, against the same foil with a
  plate catching the rest): the lost foil carries 0.9936, 1.0253 and 1.0029 of its twin
  for half-plates cut along the major and minor axes and a centered square. That is
  neither the beam's full power (2.49, 3.26, 1.73) nor no renormalization (1 within
  2 × 10⁻⁴). It fits a rescale of each beam to P × f_hit, with f_hit the share of some
  ray grid that hits an object. The grid is not identified: the 25 × 24 "cone power
  frac" grid of the capture report is ruled out (50 × 48 gives an identical
  deposition), and a search over polar and square grids found nothing beyond chance.
  **M**

#### Shadowing

Two 50° quads on foils at focus, one behind a washer (r 0.5–3 mm) and one behind a
sphere (r 0.4 mm, 0.3 mm off axis), both 5 mm upstream. VisRad conserves each quad's
4 TW over foil and occluder (3.99996). Neither shadowing along the beam axis (32% rms
per element) nor the f-cone (24–28%) reproduces VisRad, even rescaled to its totals.
A subagent's search found a model that does far better:
- each beam is rescaled after shadowing to P × f_hit (f_hit = 1 here: every ray hits
  something), which explains why the sphere foil's unshadowed side is 1.17 × the axis
  model's;
- the shadow ray from a point at spot coordinates (u, v) leans away from the axis:
  direction −k − tan θ_u (u/a) e₁ − tan θ_v (v/b) e₂, θ_u ≈ 2.0°, θ_v ≈ 1.75° (about
  1.4 × the cone half-angle per 1/e radius), as from a virtual source 19 mm (major axis)
  and 12 mm (minor) upstream of focus. Shadows and holes are magnified.

Under it the per-element rms is 5.6% (washer foil), 6.5% (washer), 11.6% (sphere foil,
mostly single-beam shadow edges) and 0.6% (sphere), against 20–33% for the axis model. Along
element cuts the sphere foil matches element by element. **M** (both occluders sit 5 mm
upstream, so an angle ∝ spot position and a fixed virtual point cannot be told apart;
occluders at ~2 and 10 mm would). The leaning shadows are not implemented: the shadow
direction depends on the point's spot coordinates, so the visibility mask would depend on
aim and spot, which a `DepositionOperator` built once per geometry does not allow.

#### Spot rotation

`spot_rotation`: ten beams (five 50°, five 23.5°) at ρ = 0, ±30, 60, 90° on plates at
best focus. The fitted major-axis angles match `rotation_sense = :beam` (ψ = 90° − ρ,
quad frame) within 0.02–0.04° on all ten, and `:target` is 60° off on the six
discriminating plates. The quad-frame reference holds to 0.04°. **H** The Laser Beam
Summary export was not saved, so `read_visrad_summary` on a real export is still
unchecked.

#### Capture/clearance report

`talis_report`: 18 beams; clearance ladders (strips at gaps g = 0.3125–1.7125 mm from
Briareus's 1% envelope, along the major or minor axis, at 20–120 mm from TCC, 10 mm
before and after focus, and on a 23.5° beam), capture plates with an edge m = 0.3–1.15
mm outside the 10% contour, and plates tilted 68.0–69.7° through focus. Stored setup:
0.1 / 0.5 mm / 70° (capture), 0.01 / 0.5 mm (3ω clearance).
- **Clearance.** Every ladder fails 0.3125–0.4625 and clears from 0.5125: a fixed
  0.5 mm, with no 0.5° term at any r (the r tan 0.5° reading predicts 0.17, 0.52,
  1.05 mm more at r = 20, 60, 120). The minor-axis ladder clears at the same rung, so
  the test volume is elliptical, and the rungs at ±10 mm from focus place VisRad's 1%
  cone within 0.05 mm of Briareus's envelope there. **H**
- **Capture.** Fractions 0.983 and 0.993 at m = 0.3 and 0.4 mm, 1 from 0.6 mm: the
  fixed 0.5 mm again. **H**
- **Incidence.** VisRad gives tilt + 1.384–1.390° on all five plates, toward either
  spot axis: the steepest ray of the f-cone (half-angle 1.380°, stored in each beam row),
  not the ±1.64° ray-grid model. On normal plates it gives 1.95° = √2 × 1.38°, the
  corner of a square cone of half-width 1.38°; a square cone also gives the tilted
  plates' extra 0.005–0.01°. **H** `talis_check` now adds atan(1/(2f#)) to the axis
  angle.

### To verify with real files at work

1. Export a Laser Beam Summary CSV from a NIF workspace with elliptical CPPs:
   - the exact headers after "Power", and whether exponent and spot definition appear
   - whether "TC Pt Pos" is the focus (as derived above)
   - cm units for the focus offset with mm preferences set
   (`spot_rotation.vrw` of the 2026-10-27 kit, step 4 of its README, is set up for it.)
2. Done (2026-10-27): rotation +30° turns the minor axis clockwise looking along the
   beam (`:beam`); see "VisRad probe kits (2026-10-27): results".
3. Export a NIF SST file. Check the header strings, any preamble or extra columns, and
   whether off beams are listed.
4. Save Data for Object on a sphere. Check the delimiter, the header spellings
   (especially the polar-angle column) and whether rows carry an index.
5. Obtain a sample CMT XML to implement the CMT reader.


## Example NIF cone pointing: unsourced (2026-10-27) *decision*

The cone aims used by the README quick start, `examples/nif_hohlraum.jl`, the tutorials
(Freeing a bundle for radiography and those built on it), several examples and tests,
and BriareusGUI's NIF presets (`NIF_CONE_POINTING`) are, in hohlraum (z, r) coordinates
with r toward the port azimuth of the quad (mm):

| cone | z | r |
|---|---|---|
| 23.5° | 3.09 | −1.09 |
| 30° | 3.63 | −0.60 |
| 44.5° | 5.20 | 0 |
| 50° | 5.01 | 0 |

They date from the initial commit (2026-10-15), labeled "published", with no citation.
A search of this file, `reference/litreview/`, the `data/` headers, the docstrings and
the history (`git log -S`) found no source: no paper, report or public file in the
repository gives them, and Build-A-Hohlraum (Swadling et al., PoP 32, 052707, 2025),
whose 5.75 mm diameter the example hohlraum shares, is recorded for its LEH closure and
calibration only, not for aims. They are therefore an **illustrative example pointing**,
not a facility or shot pointing, and are labeled so everywhere (**L** for any claim that
they resemble a fielded design). The form is the usual indirect-drive one: outer cones
aimed on the axis near the LEH planes, inner cones across the axis. Where a shot-like
pointing matters, `hybride_mock_pointing` (section "Mock HYBRID-E pointing
(2026-10-23)") is the documented stand-in. Earlier sections that said "published cone
aims" were relabeled "example" on this date.

Since 2026-11-02 BriareusGUI's NIF presets hold these aims as offsets from the default
10.1 mm hohlraum part's entrance holes (z − 5.05 mm in the hole's frame), so on a
halfraum or a hohlraum of another length they keep the same offsets from its holes
(`preset_pointing`; reference/assemblies-design.md, "Item 5 as built"). That carry-over
is an engineering default for illustrative designs (**L**), not a pointing rule.

**Quick start, re-posed** *derived*. The previous quick start optimized z, r and power
per cone (12 parameters) against P2, P4 and the mean flux only: a flat valley, stopped at
the 200-iteration limit (objective 1.9e-8, outer44 r at −0.98 against its −1 bound,
powers near their bounds), ~25 s on first run. The new one starts from the example aims
and moves each ring rigidly (`GroupPointing(quads; by = ring, fields = (:z, :power),
relative = true)`, inner = 23.5° + 30°, outer = 44.5° + 50°): one aim-height shift and
one power per ring (default hohlraum wall, 5.75 × 10.1 mm, 3.1 mm LEHs, 0.9 mm capsule,
`nif_quads(power = 8.0)`).
- Without keep-outs, the 4 parameters against P2, P4, P6 and the mean flux converge in
  31 iterations, none at a bound (inner −0.607 mm, outer +0.387 mm, 5.42 and 10.54 TW
  per quad), but the result clips: worst-case LEH clearance (`aperture_clearance`,
  default 1% envelope and errors) −0.265 mm on the 30° cone, −0.05 mm on the outer
  cones, and the 23.5° envelope 0.013 mm into the capsule keep-out. At the start the
  clearances are +0.066 mm (LEH, 23.5°) and +0.194 mm (capsule).
- Adopted: `ApertureKeepOut` + `SphereKeepOut` in `symmetry_objective`, P2 and P4 and the
  mean flux, power box ±60% (the inner power sits just below the default −50%). The
  LEH keep-out is active at the optimum (30° clearance +0.001 mm), so three targets and
  one active constraint fix the four parameters. 20 iterations, 57 evaluations, none at
  a bound; objective 0.02477 → 3.50e-8; P2/P0 −15.74% → −0.000%, P4/P0 +0.14% →
  +0.000%, P6/P0 +0.046% → +0.077%; inner −0.135 mm, outer +0.298 mm, 3.97 and
  11.21 TW per quad (inner-cone fraction 1/3 → 0.15), delivered power 384 → 422.1 TW at
  the same mean capsule flux (Tr 262.7 eV); capsule clearance +0.148 mm. Adding P6
  (weight 1) reaches the same point in 52 iterations and 221 evaluations, so it is left
  out. The ray trace of `examples/nif_hohlraum.jl` at the optimum gives P2/P0 −0.01%, P4/P0
  −0.00%, P6/P0 +0.08%.
- Timing, 16 threads (Ryzen, Julia 1.13.2): `using Briareus` 1.6 s; scenario build
  (deposition operator, view factors) 2.5 s; optimization 4.9 s on first call
  (compilation) and 1.0 s after.

**Starting P2 at equal power.** At 8 TW per quad (inner-cone fraction 1/3, near HYBRID-E's
input peak cone fraction: "N221204 benchmark (`n221204_benchmark.md`)" under "Literature
review against Kritcher et al. 2024 (2026-10-22)", **L**) the static capsule flux has
P2/P0 −15.7%. A static solve is not expected to be balanced at a shot's cone fraction:
in the N221204 benchmark the static P2/P0 at 5.5 ns is −0.26 even with the calibrated
CBET transfer (section "Two-zone CBET calibrated to PRE 109 Fig. 9 (2026-10-23)");
fielded balance relies on wavelength-detuned CBET and the hohlraum's dynamics (wall
motion, the gold bubble, LEH closure).

## VisRad benchmark beyond deposition: the saved runs (2026-11-18) *measured*

ROADMAP "Benchmark against VisRad". The saved VisRad runs (VisRad 21, outside the
repository, non-public; described generically here) were searched for VisRad's radiation
results, and Briareus's radiation side was measured on them with the VisRad deposition
as the source and with Briareus's.

- **What the saved runs hold.** Every computed workspace (hohlraum targets built from
  Briareus geometry as VisRad mesh objects, one VisRad element per triangle, at a ~0.6 mm
  and a 0.1 mm mesh; the probe kits; a TaLIS report kit; the primitive hohlraum + foil
  sample) holds VisRad's power sources at one time (t = 0, constant power) and nothing
  else of VisRad's solution: `[Radiosity Data]` and `[View Factor Data]` read "not up to
  date" or are empty in every file and in every version the study repository has
  committed. The JSON workspaces read for "VisRad JSON workspaces (2026-11-16)" are no
  longer on disk; the layout recorded from them (the test mock) has only the count
  `"Num radiosity up to date"` under `"Radiosity data"`. No radiosity, emission or radiation temperature, incident
  flux on the capsule, view factor, spectral or time-dependent result exists in them.
  The TaLIS kit also holds beam capture and clearance data, and the probe kits `.bti`
  intersection exports; neither is a radiation result. **H** (every file searched)
- **What the readers take.** `visrad_power_sources` reads every column of the power
  sources (direct laser source, uniform source, albedo, conversion, reflectivity,
  pre-reflection power, cone counts, overlap angles) in both formats. A populated
  radiosity section has not been seen in either format, so no reader was written for
  it; capture-clearance tables are not read (the writer resets them). **H**
- **The import's materials.** On all nine computed workspaces VisRad's per-element
  albedo and conversion columns equal the imported objects' values (to 10⁻⁵), the
  uniform source and laser reflectivity are zero, and the element counts match the
  meshes; on the primitive sample the only elements that differ are the hohlraum's
  entrance-hole ring at α = η = 0, which the import makes apertures. VisRad's
  source term is then ηI on every element, as in Briareus. **H**
- **Method** (scratch scripts kept with the runs). The target is the workspace's
  objects on VisRad's elements (`visrad_target(ws)`); laser deposition from VisRad's
  power sources (`visrad_deposition`), from Briareus with VisRad's profile law
  (`visrad_profile(ws)`) and with the default f-cone model, each on the whole target.
  The radiation solve takes the radiating parts inside the cavity (wall, caps, corner
  pieces, window patches, capsule): the parts outside it have α = η = 0, so B = 0 there,
  and a segment between two points inside the cylinder cannot cross them, so dropping
  them leaves the cavity's solution unchanged. One Briareus model per mesh (raw view
  factors, default `:sampled` visibility, constant albedo and conversion as imported)
  solves all three sources; differences are taken against the solve on VisRad's
  deposition. The regions of the gold are the elements whose laser source ηI exceeds B/5
  ("spots", 32–33% of the gold's area), the corners (within 0.6 mm of the wall–cap
  edge), the entrance-hole rim (within 0.6 mm of it) and the rest. Modes are the
  capsule's incident flux in real spherical harmonics by exact solid angles
  (`ylm_weights`), as Legendre ratios `aₗ/a₀` (m = 0) and per-l amplitudes
  `√((2l+1) Σₘ cₗₘ²)/c₀₀`, in percent.
- **Deposition propagated through the radiation, ~0.6 mm mesh** (four runs: one
  pointing with and without two window patches, two re-pointings; about 5.5k radiating
  elements). Briareus with VisRad's profile law against VisRad's deposition:
  - laser flux on the gold 0.86–0.94% rms per element (largest 4.0–4.8% of the mean);
    absorbed by the gold 0.07% below VisRad;
  - gold radiosity B 0.16–0.17% rms, largest 0.72–0.86% of the mean; window patches
    0.07–0.13%; 56–63% of Σ A ΔB² under the spots, 32–39% on the rest of the barrel,
    3–4% in the corners, 1% at the entrance-hole rim (each region 0.08–0.19% rms);
  - capsule incident flux −0.050 to −0.054% (Tr −0.04 eV); the power leaving through
    the openings −0.08%; mean gold B −0.07%;
  - capsule modes: ΔP2 −0.019 to −0.032 pp, ΔP4 −0.005 to −0.006 pp, ΔP1, ΔP3 ≤ 0.003 pp,
    the m ≠ 0 parts of l = 1–4 within 0.006 pp.

  Briareus's default f-cone model instead: laser flux 12.5–12.8% rms (largest
  66–75%), B 2.0–2.1% rms (largest 10.8–12.4%), 85–87% of it under the spots;
  capsule flux +0.075 to +0.078% (Tr +0.056 eV); openings −0.02 to −0.03%; ΔP2 −0.09 to
  −0.11 pp, ΔP4 +0.025 to +0.040 pp, the m ≠ 0 part of l = 2 up to 0.037 pp. The wall's
  high albedo smooths the deposition difference five- to sixfold into the
  radiosity, and again on its way to the capsule.
- **The sensitivities of the radiation solve itself** (same mesh, VisRad's deposition), the
  size of the known modeling differences a VisRad radiosity comparison would show:
  - visibility by one centroid ray (`:centroid`) instead of sampled partial visibility:
    capsule flux −0.12% (Tr −0.09 eV), openings +0.28%, gold B 0.14% rms, ΔP2 −0.003
    pp, ΔP4 +0.002 pp;
  - view factors by the centroid point kernel A cos θ cos θ′/(πR²) on every pair
    (`near = 0`, `exact_far = false`, `:centroid`): capsule flux +6.2% (Tr +4.5 eV),
    openings −16.7%, gold B 9.5% rms (20% in the corners, where neighbors share an
    edge), ΔP2 +1.5 pp, ΔP4 +0.07 pp. A code that integrates near pairs this crudely
    would show it at once in the drive level and the loss through the openings.
- **The same at the 0.1 mm mesh** (three runs on one geometry: a shot-like pointing and
  two re-pointings at a lower power; about 60k radiating elements; view factors as an
  `HMatrix` with `:centroid` visibility, since `:sampled` did not assemble within the
  run's five-minute budget on 16 threads; on the 0.6 mm mesh the deposition-propagated
  differences above are the same under either visibility rule to the digits quoted).
  Briareus with VisRad's profile law against VisRad's deposition:
  - laser flux on the gold 0.56–0.71% rms (largest 2.0–2.5%); absorbed by the gold
    0.08–0.09% below VisRad;
  - gold B 0.13–0.15% rms (largest 0.41–0.49%); 60–67% of Σ A ΔB² under the spots,
    12–22% in the rounded corners (38% of the gold's area here), 19–21% on the rest,
    2–3% at the rim;
  - capsule flux −0.065 to −0.077% (Tr −0.04 to −0.05 eV); openings −0.09 to −0.10%;
  - ΔP2 −0.009 to −0.013 pp, ΔP4 −0.009 to −0.010 pp, odd l and m ≠ 0 within 0.002 pp.

  The default f-cone model: laser flux 12.0–14.5% rms, B 2.0–2.4% rms (largest
  9.2–12.7%), 86–88% under the spots; capsule +0.046 to +0.055% (Tr +0.03 to +0.04 eV);
  openings −0.03 to −0.04%; ΔP2 −0.108 to −0.117 pp, ΔP4 +0.041 to +0.042 pp, the m ≠ 0
  part of l = 4 +0.02 pp.
- **Across the two meshes** (different cases, so a trend, not a convergence study):
  with VisRad's profile law the per-element differences shrink with the mesh (laser
  0.9 → 0.6–0.7% rms, B 0.16 → 0.13–0.15%, ΔP2 0.02–0.03 → 0.01 pp), while the drive
  level does not: the capsule flux and the loss through the openings follow the gold's
  total deposition (−0.07% → −0.08 to −0.09%), so the level difference is VisRad's
  per-beam total (ROADMAP item 6), not the spot shapes. The f-cone model's differences
  are mesh-independent (its spots are another shape).
- **Explained and open.** Explained: every radiation-side difference measured here is
  the deposition's, carried through one radiosity solve: the level by the beam totals,
  the shape by the spot profile (VisRad's law against the f-cone). Open, and not
  measurable from the saved runs: the VisRad radiosity, so whether its view factors,
  visibility, closure and albedo bookkeeping agree with Briareus's. With VisRad's
  deposition as the source, a VisRad radiosity solve that agrees with Briareus's
  should match to well under the f-cone differences above; a difference of the size of
  the `:centroid`/`:sampled` sensitivity (0.1% in drive, 0.3% in the openings) would
  point at partial visibility, one of the point-kernel size at the view-factor
  quadrature.
- **What would settle it** (no new kit needed: the saved hohlraum workspaces hold the
  geometry, beams and VisRad's up-to-date deposition). In VisRad, open a saved run per
  mesh (one at 0.6 mm, one at 0.1 mm), compute the radiosities, save, and export the
  Surface Element Properties (Flux In, Flux Out, Rad T, Emis T with the surface index)
  of the gold objects and the capsule, and the Power Summary. Back in Briareus:
  `compare_visrad(ws, sol; elements)` with `sol` solved on `visrad_target(ws)` from
  `visrad_deposition(ws, t)`. The saved file then also shows VisRad's populated
  `Radiosity data` layout, for a reader. `write_visrad_kit` cases remain only for what
  the Briareus-built meshes cannot show: VisRad's primitives with curved surfaces
  (flat facets renormalized to the nominal area), a cavity with a closed entrance hole
  (balanced view factors), and temperature-dependent albedo.
- **Code.** `compare_visrad(workspace, sol; elements, observe)` compares a radiation
  solution with a workspace's results element by element: laser deposition from the
  power sources, Flux In, Flux Out, Rad T and Emis T from the element exports, the
  incident and emitted powers and Tr per object, and the observed object's modes by
  exact solid angles on both sides. `visrad_element_values` maps one exported column
  onto the target's elements by the surface index (0- or 1-based, *assumption* until a
  real export is seen). Checked on the synthetic workspace with a CAD capsule whose
  exported field is f₀(1 + a₁P₁ + a₂P₂): P1 and P2 recovered to 0.003 and 0.02 pp, and a
  2% offset on another object recovered exactly. On a saved run with Briareus's
  deposition its laser rows reproduce the comparison above (gold 0.07% below VisRad).
  **H** (the comparison itself), **M** (the export layout: header names from the
  manual's screenshot, index base assumed)
