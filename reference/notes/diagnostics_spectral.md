# Research notes: Spectra, Dante, imagers and backlighters

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

## Spectral radiosity, M-band and a synthetic Dante (2026-10-24) *derived*

ROADMAP "Spectral radiosity" and the first absolute check of "Absolute energetics on N221204". Read (accepted
manuscripts): Dewald et al., "Direct measurement of
hohlraum M-band preheat asymmetry at the ICF capsule", PoP 27 (2020), LLNL-JRNL-813556;
Rubery et al., "Soft x-ray power diagnostics for fusion experiments at NIF, Omega and Z
facilities", RSI (2023), OSTI 1975922. Confidence: **H** primary text, **M** quoted
secondhand, **L** our construction, *calibrated* as marked.

### Inputs read

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| M-band definition | > 1.8 keV (Au); Mo fluorescence sensitive > 2.5 keV | Dewald 2020 §II | H |
| Laser → M-band, Au, normal incidence | ~1% at 10¹⁴, 20% at 10¹⁵ W/cm² | Dewald 2020 §I citing Dewald PoP 15, 072706 (2008) | M |
| Wall irradiance | inner 3 × 10¹⁴, outer 1.5 × 10¹⁵ W/cm²; inner turning point 2.5× lower density (M-band ∝ n_e² ∝ sin⁴θ) | Dewald 2020 §I | H |
| Dante M-band share | 12% of the flux at the 290 eV peak (Fig. 1b); 15% of the spectrum, 1/3 of it > 2.5 keV (Fig. 3) | Dewald 2020 §II–III | H |
| Setup | 10.1 × 5.75 mm Au, 3.37 mm LEHs, 360 TW / 1 MJ, inner/outer power 1/2, outer ring at z = ±3.3 mm (EHXI), θ_o ≈ 40°, R_o 4.4 mm, wall ~0.2 mm in, HDC front 0.815 mm at 4.35 ns | Dewald 2020 §II, V | H |
| Mo-sphere P2/P0 | swings, mostly > 0, +0.15 at 6 ns | Dewald 2020 Fig. 10 | H-img |
| Wall M-band P2/P0 | +0.53 (= 0.15/0.28); a pure ring at 40° gives 1.9, so F_M ≈ 2 (recirculating / laser-produced) with the LEH deficit (−25%) | Dewald 2020 §V.2 | H |
| Smoothing | P2: 0.4 at the HDC ablation front (CCR 5.4), 0.28 at the Mo sphere; inferred M-band P2/P0 at the inner HDC surface ≤ +0.28 | Dewald 2020 Figs. 6, 12–14 | H |
| Dante 1 channels | 18 channels, 0.25–3.66 keV and 14.4 keV (table of energies, filters, mirrors, XRD cathodes); the unfold fits a Planckian to channels < 1.5 keV then Gaussian responses | Rubery 2023 Fig. 4, §II.C | H |
| Dante lines of sight | Dante 1 θ 143°, φ 274°, 37° view of the LEH; Dante 2 θ 64°, φ 350°, shallow | Rubery 2023 §II | H |

### 1. Multigroup radiosity (`spectral_solve`)

View factors are spectrally independent, so each group g solves the gray system with
a separate source and albedo, `B_g = S_g + α_g H_g`. The source ηP/A minus the laser M-band
is split by Planck fractions `b_g(T)` at each element's emission brightness temperature
`(B/σ)^¼` (iterated when the totals change); the M-band is flat over its band. With
`thermalize = true` the rest of the gray albedo, `α − α_g`, is re-emitted with `b_g(T)`,
so the totals are the gray solution's whatever the group albedos (*energy as calibrated
gray*); without it, it is lost to the wall. Groups at the gray albedo reuse the model's
cached LU (or GMRES preconditioned with it; the GMRES became generic in the operator,
and the gray path is unchanged); hard groups are solved jointly by unpreconditioned GMRES.

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
  11.9–12.5% of its flux above 1.8 keV against Dewald's measured 12% at 290 eV: an
  independent agreement (Dante's M-band share was not fitted). (With the laser spots
  thermal at the local Tr, the 2026-10-24 default, F_M ≈ 2 takes α_M ≈ 0.55 and Dante
  sees 11.4–12.0%: "Laser-spot spectrum" below.)
- The M-band is pole-hot relative to the thermal drive by +0.15 to +0.21 in P2/P0 at the
  capsule (M-band minus gray), against Dewald's ≈ +0.21 at the ablation front (0.53 ×
  0.4, **L**). The absolute M-band P2 here (+0.04 / +0.15) sits lower because this mock
  pointing's thermal drive is waist-hot (−0.11 / −0.07) and the recirculating M-band
  (2/3 of it) carries that shape, where Dewald takes it as uniform. Likewise the
  boundary P2/P0 (−0.01 / +0.20) against his +0.53, which assumes a uniform background.
- Two effects are not modeled: the outer-spot motion and inner-beam clipping in the outer
  bubble that Dewald invokes for the late swing, and the angle dependence (sin⁴θ) of the
  conversion.

### 3. Synthetic Dante (`DanteView`, `dante1`, `dante`)

Parallel rays (400 × 400) across the LEH give each element's visible projected area,
partial views included; a ray through a closing lip is weighted by its open fraction.
`Tr = (π I / (σ A cos θ_v))^¼` over the full initial LEH area (Dante's convention),
`Tr_open` over the open area. Channels: Gaussians at Rubery's 18 Dante 1 energies, FWHM
0.3 E (**L**, widths unpublished; since replaced by the hardware responses, "Dante channel
responses and the unfold" below); `color_temperature` fits a Planckian to the channels
below 1.5 keV. Checks: a blackbody cavity gives T₀ to 5e-4 (after the caps' inscribed
24-gon opening, 1.1% below the equal-area hole on the coarse test mesh: on coarse
meshes the opening is smaller than π r²); lips at open fractions 0, ½, 1; a Planckian's
channel signals unfold to its T to 1e-6.

N221204 (`examples/n221204.jl` full model, standard mesh, default energetics):
| | peak | at 8.5 ns |
|---|---|---|
| capsule Tr (mean incident flux) | 301.8 eV at 5.50 ns | 270.3 |
| Dante 1, initial LEH area (Dante's convention) | 285.9 at 6.25 ns | 257.4 |
| Dante 1, open area (lips at their open fraction) | 306.2 at 6.25 ns | 278.1 |
| measured Dante 1 (Rubery 2024; Pak 2024 Fig. 1) | 313 ± 2 at ≈8.4 ns | ≈310 |

- Dante 1 sees only the wall (no capsule, no cap) at 37°. Over the open area it gives
  7 eV above the capsule at the Dante peak, the sign of PRE 109's "internal ~13 eV below
  Dante" (**L**: the simulated internal temperature is not exactly our capsule mean).
  Over the initial area it gives 13 eV *below* the capsule: the lips (closing to
  0.857 r₀) are dark to Dante here, while the real closing plasma radiates. Which of the
  two the measurement corresponds to depends on the LEH plasma's emission, which this
  model does not have; the two bracket the measured 313 only at the open-area end, and
  7 eV low there.
- The larger miss is the time history: the model's drive peaks at 5.5–6.25 ns and falls
  ~30 eV by 8.5 ns, where Dante's keeps rising to its 8.4 ns peak. The default
  energetics (Au wall loss, constant conversion 0.75, ablator at its initial radius, no
  backscatter) and the bubble's growing inner-beam absorption lose drive through the
  peak; this is the first absolute check of item 16. It points at the late-time
  energetics (wall loss vs time, conversion rising in time, the `energetics=new` options);
  the line of sight is not the cause. (Resolved in "N221204 late-time drive" below: the
  wall loss was evaluated at the time since the pulse start; these numbers are the old
  evaluation.)

## Backlighter foil on the line of sight (2026-10-24) *derived*

ROADMAP "Backlighter geometry", second half. `backlighter_foil(los; standoff, width, height, tilt,
tilt_toward, rotation)` and `line_of_sight(target, los; ...) -> LineOfSightReport`
(`src/backlighter.jl`, `test/backlighter.jl`), used in the radiography tutorial
(`docs/literate/dropped_bundle.jl`).
- **Defaults are engineering values.** The 2 mm square gold foil and the standoffs tried
  are generic: real backlighter dimensions and standoffs are shot-specific. The radiography
  windows (two open 3 mm disks in the wall on the line of sight) are likewise generic. **E**
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

## Backlighter conversion efficiency and the photometric chain (2026-10-24) *engineering*

ROADMAP "Backlighter physics". `kshell_line_energy`, `kshell_ce_data`, `kshell_ce_model`, `kshell_ce`,
`mass_attenuation`, `filter_transmission`, `backlighter_signal` (`src/backlighter_signal.jl`,
`test/backlighter_signal.jl`, `data/kshell_ce.csv`); used in the radiography tutorial
(`docs/literate/dropped_bundle.jl`, "The radiograph signal") and `examples/backlighter.jl`.
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
Mo 0.15%. The peak falls as hν^-3.5 (foil, from Fe to Ge), inside the published hν^-(2.5–4).
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
   element at its irradiance; σ stays the mean's), or a given `ce`.
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
  comment states that a measured `ce` should be passed.

**Not modeled / next.** Pulse-length and spot-size dependence (data confounded);
backside emission of thin foils; the line's opacity and limb brightening (the isotropic
assumption); continuum and higher lines in a filter's pass band; detector response
(image plate, CCD QE) and the object's transmission. For the low-Z end, see the next
section.

## Weighted He-α CE compilation (2026-10-26) *engineering*

ROADMAP "Backlighter physics", "a weighted, larger low-Z compilation (Ti, V)". `data/kshell_ce.csv` grows
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
  realigned by the checks in the paper, s72639 = 7.8 kJ at 2e15 and s77768 = 0.3e15 at
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
  Zr 1.1e16 / 0.20% → 8.1e15 / 0.21%. Peak slope from Fe to Ge (foil) hν^-3.3, from Cu to Ge
  (prepulsed) hν^-3.6.
- **CE at representative points, old → new** (%, ratio): foil Cr 1e15 3.09 → 4.14
  (×1.34); Fe 1e15 1.62 → 2.05 (×1.27); Fe 3e15 1.48 → 2.01 (×1.35); Cu 1e15 0.67 →
  0.78 (×1.17); Ge 3e15 0.36 → 0.51 (×1.42); prepulsed Cu 1e15 1.52 → 1.59 (×1.05); Zn
  2e15 1.22 → 1.25 (×1.02); Ge 5.1e14 0.56 → 0.49 (×0.86); Ge 1.5e15 0.74 → 0.74
  (×1.01); Zr 4.5e15 0.18 → 0.20 (×1.08); Mo 1e16 0.15 → 0.17 (×1.12); Zr 5.1e14 0.078
  → 0.050 (×0.65, extrapolated). σ 0.43–0.48 → 0.37–0.42 inside the data. Ti foil
  4.6e14 6.7 → 8.9 (×1.34, extrapolated; see below).

**Held-out checks** (independent data, not the points fitted).
- The form fitted only to the 2026-10-24 rows (weighted, 37 points) predicts C19's 31
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
offset by itself fits +0.19 ± 0.11 (two-sided irradiation or calibration; not significant).
On spot size, every source sets the irradiance with the spot (C19 and W01 change spot and
beam count together), so a spot term is degenerate with X here. *E*

## He-α CE: saturating peak (2026-10-26) *engineering*

The form of `kshell_ce_model` changes; the data (`data/kshell_ce.csv`, 68 fitted rows),
the weights and τ rule of the previous section are unchanged. Exploration script not in
the repository; the package fit reproduces its numbers.

**Diagnosis.** The previous form's peak over lines is `a_c + bL + (c_c + eL)²/(4|d|)`:
the `X L` drift makes it convex in `L`, so the peak rises faster than any power law
toward low hν (Cr 4.15%, Ti 9.2%). The shared `L²` term tried before only offsets part of
that curvature. Writing the parabola of each line about its optimum separates the peak
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
was computed differently and its C19 value disagrees with the C19 held-out value of that
line, 0.35). LOEO, adopted (mean/rms): Cr see below, Fe +0.38/0.54, Ni +0.03, Cu −0.17/0.35,
Zn −0.29/0.33, Ge +0.10/0.54, Zr +0.12/0.23, Nb +0.44/0.47, Mo +0.30/0.30 (the previous
form: Cr −0.42, Fe +0.41, Nb −0.03, Mo −0.07). The cubic's in-sample gain is overfitting
(rms 3.8 with W01 left out).

**What the comparison shows.** (1) Most of the cross-validated gain is the
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
(at Fe; the value the refined grid finds). (5) The prepulsed high-Z end (B14 Zr/Nb/Mo)
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
a turnover in place of a plateau; one point cannot separate them. *confidence: medium
inside the data (unchanged), low below Cr.*

**Tests** (`test/backlighter_signal.jl`): χ²/(68 − 8) = 1, weighted residual sums zero
per class, the C19 held-out check (normalized rms < 1.2), the held-out Ti within 1.5σ,
hνₛ between Cr and Fe, the Ti and Cr peaks equal (the plateau) and above Fe's, the
gradient behind σ against central differences, the optimum a maximum, the hν slope in
2.5–4.5.

## Dante channel responses and the unfold (2026-10-24) *derived*

> **Corrected in part (2026-10-29).** The grids are 25 µm tantalum (Kemp 2020, Fournier
> 2010, May 2015), not tungsten; as opaque wires the responses do not change. See
> "Dante-family hardware and geometry: what is public (2026-10-29)".

ROADMAP "Spectral radiosity", what remained of the synthetic Dante: the 18 Dante 1 channels built from
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
| Dante 1 hardware, 18 channels | energies, filters (µm), 8.8% W grids, mirrors (ch 2 2–4 keV multilayer, ch 3 3.5° SiO₂, ch 4 2.5° SiO₂), XRD cathodes (Al/Cr/Ni), solid angles (6.51e-7 sr; ch 2 1.90e-7, 3 5.74e-7, 4 5.56e-7, 16 and 18 1.77e-6), uncertainties 8–16% | Rubery 2023 Fig. 4 (a table image, read at 300 dpi) | H |
| Channel 2 filter | "47 µm CH" in Fig. 4; 47 µm polyimide in §II.D and Fig. 8 — polyimide used | Rubery 2023 | H |
| Grids | 8.8% open tungsten, behind the filters, to keep XRDs below the 500 V Child–Langmuir limit | Rubery 2023 §II.A | H |
| XRD-31 | reflection-mode planar vacuum diodes, Ni/Al/Cr/Au cathodes, 5 kV bias, 140 ps | Rubery 2023 §II.A | H |
| Cr XRD sensitivity | 0.07–1 keV, C K (277), O K (525), Cr L (573 eV) structure; axis labeled V/GW | Rubery 2023 Fig. 3b | H-img |
| Al XRD sensitivity | 1–4.5 keV, MV/GW, Al K edge | Rubery 2023 Fig. 8b | H-img |
| SiO₂ mirror "3.5°" | measured reflectivity 0–1 keV (0.84 → 0.06) | Rubery 2023 Fig. 3a | H-img |
| Ti 0.75 µm filter | measured "effective transmission" 0–1 keV | Rubery 2023 Fig. 3a | H-img |
| Cr/Sc multilayer, 1.9° | reflectivity 1.2–4.5 keV; system response (mirror × polyimide × Al XRD) in kV/GW | Rubery 2023 Fig. 8 | H-img |
| Channel 2 calibration | 0.97 TW/sr in 2–4 keV per volt at the XRD, 8.7% | Rubery 2023 §II.D | H |
| Signal equation | V = Ω ∫ S(E) R(E) dE, S in GW/sr/eV, R in V/GW | Harris 2021 eq. 1 | H |
| Standard unfold | Planckian start from the channels < 1.5 keV (Rubery) / channels 1–9 (Harris), averaged per-channel brightness temperatures, plus an M-band Gaussian (1 GW/eV/sr, 2.5 keV, width 600 eV); each response replaced by a fitted Gaussian whose amplitude is adjusted channel by channel from the highest energy down, ~10 iterations to 5% | Harris 2021 §II; Rubery 2023 §II.C | H |
| Sentinel | DIM-mounted, 16 MiniXRDs (6 Al, 5 Ni, 5 Cr) behind Dante filters, polar DIM 000-000 or equatorial 090-078 / 090-315, ~7 m | Kemp 2020 §II | H |
| Optical constants | Henke f₁, f₂, 10 eV–30 keV, 19 elements (`data/henke/`, CXRO, fetched 2026-10-24) | Henke, Gullikson & Davis, ADNDT 54, 181 (1993) + CXRO updates | H |

Digitization (`data/rubery2023_dante.csv`): the figures are raster images in the accepted
manuscript; the plot frames were located from the dark frame lines and the axis ends taken
at the frame (where the end tick labels sit), each curve masked by color and its median row
read per pixel column. Overlaying the read points on the figures agreed to about a pixel
(~0.1% of full scale); the narrow edge spikes (XANES) are under-read by the median. Units:
the axis of Fig. 3b is labeled V/GW; read as MV/GW like Fig. 8b's Al XRD — with V/GW the
quantum efficiency would be ~10⁻⁸, and with MV/GW channel 2 reproduces its published
calibration (below).

### Model

- **Filters**: Henke photoabsorption (`henke_mass_attenuation`; `mass_attenuation(...;
  table = :henke)`, with `:auto` keeping NIST where it covered the backlighter filters
  before, so `backlighter_signal` is unchanged). Henke vs NIST at 3 keV: Be, Al, Ti, Fe,
  Cu, polyethylene within 3% (NIST includes scattering). Compound densities CH 1.05,
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
  MV/GW, L 33.6 nm, w_C 0.274, w_O 0.250, γ −0.244; rms 12%) and to Fig. 8b (Al: K 35.12
  MV/GW, γ 0.102 on Cr's surface layer; rms 6.7%; the K-edge region 1.54–1.75 keV
  excluded from the fit). Ni (no published curve) takes Cr's parameters on Ni (**L**);
  Cr above 1 keV and Al below 1 keV are extrapolations of the fitted form (**L**). As
  quantum efficiencies: Al 0.016 electrons/photon at 1 keV, Cr 0.03 at 70 eV.
- **Checks against independent published numbers**: channel 2 (multilayer × 47 µm
  polyimide × Al XRD, Ω 1.9e-7 sr) gives **0.968 TW/sr per volt** for power flat over
  2–4 keV against Rubery's **0.97**; its response is 0.942 × Fig. 8b's over 2–4 keV
  (rms 6%). Not fitted. The Ti filter of Fig. 3a transmits 0.45 × Henke's 0.75 µm at
  0.45 keV and 0.83 × at 1 keV: measured filters absorb more than their nominal
  thickness at low energy (oxide, contamination, support), a caveat for every thin-filter
  channel here, all of which use nominal thicknesses.

Channel bands (half maximum, eV) against Rubery's nominal energies: 1 (1170) 1084–1305;
2 (3060) 2307–3625; 3 (250) 226–284; 4 (470) 444–513; 5 (710) 656–779; 6 (870) 867–934;
7/15 (940) 852–1024; 8 = 0.088 × 1; 9 (1430) 1314–1561; 10/14 (1710/1700) 1558–1840;
11 (3430) 1821–4467; 12 (2580) 2351–2825; 13 (3000) 2695–3357; 16 (14370) 8579–19247;
17 (690) 624–779; 18 (3660) 1771–4200. Each band's peak lies within 4–13% above the
nominal energy (the absorption edge that closes it), except the broad 11, 16 and 18.
A Planckian at 300 eV through the LEH (π 1.55², 37°) gives 190 V on channel 1 and 14 V
on channel 3 (at the XRD, before attenuators), the scale of Harris's N190721 traces
(90–230 V on the low channels).

### Interpolation of tabulated curves (2026-11-03) *derived*

`TabulatedCurve` (a measured filter transmission, mirror reflectivity, cathode
sensitivity or whole response read from CSV) interpolates linearly. The built-in
channels do not use it: their filters are Beer–Lambert from Henke's f₂, itself
interpolated log-log between the Henke energies, so the default Dante 1 responses and
every Dante Tr are independent of this choice. This was measured on the filters of the 18
Dante 1 channels, with each transmission tabulated on a log grid (50 eV – 30 keV), read
back by three rules, and the channel response (`R` × read/exact transmission) weighted
by Planckians at 100, 200 and 300 eV against the exact one on the 3001-point Dante
grid; the worst of the three temperatures:

| Table | linear | log-log | log(−ln T) vs log E |
|---|---|---|---|
| 167 per decade (CXRO's 500 points from 30 eV) | 0.1–4.0% (44% for 250 µm Al) | 0.1–8.8% | 0.05–8.0% |
| 20 per decade | 7–87% | 1.4–99% | 0.3–98% |

Neither rule wins: the error is set by the absorption edges falling between two points
(the edge that closes a channel's band), where linear reading is the less biased; the
power-law rules win on smooth steep curves (the thick Al filter of channel 16: 44% →
0.3% / 0.05%). Linear is kept for every kind (**L**, a decision); a measured table
should carry the energies of the edges, as CXRO's and Henke's do.

### Unfold (`dante_unfold`, `dante(view, ::SpectralSolution)`)

- `:traditional` — Harris's algorithm: Planckian start at the mean of the per-channel
  brightness temperatures of the channels below 1.5 keV, plus the M-band Gaussian; each
  response replaced by a unit Gaussian (center between its half-maximum points, FWHM
  their distance: Harris fits one but does not say how, **L**); amplitudes adjusted
  channel by channel from the top. With Dante 1's duplicated hardware (1/8 differ only
  by the grid, 7/15 and 10/14 are identical) the Gaussians coincide and the sweep
  diverges; the same equations are then solved directly (least squares in the relative
  voltages, minimum-norm amplitudes).
- `:smooth` — positive: the start Planckian (no M-band) times a piecewise-linear factor
  on 48 log-spaced knots, nonnegative least squares on the relative voltages with a
  second-difference penalty (λ 1e-3). A Planckian reproduces itself (factor 1), so
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
Planckian at the brightness temperature of the element ("Planckian"), and with the
intensity-form laser M-band (`MBandEmission()`, "+ M-band"). Within each group the voltages
take the spectrum as a Planckian at Tr_open (**L**). Tr over the initial LEH area (eV):

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

- **The unfold recovers the input's Tr**: the traditional unfold 0.8–1.3 eV high (0.3–0.4%
  in Tr, ~1.5% in flux; its Gaussians spread the high channels' power into the
  thermal band), the smooth one within 0.05 eV. Both well inside Dante's quoted 1σ of
  1.5% in Tr (Rubery 2023 §II.C). So comparing the synthetic Dante's `Tr` with the
  measured Dante Tr, as the N221204 tables do, was sound to ~1 eV; the brackets of
  initial vs open area (24 eV) and the late-time energetics remain the open terms.
- **Before/after for the channels**: with Gaussians there was no unfold, only the
  color temperature of the channels below 1.5 keV. That color temperature agrees
  between Gaussian and hardware channels for the Planckian spectra (≤ 1.5 eV), but
  with the laser M-band the hardware channels give 5–33 eV hotter values: their high-energy
  leakage (nominal-energy bands with tails past the absorption edges) picks up the
  M-band, which Gaussians cannot. The color temperature is a start value, not a
  result, so Tr is unaffected; it is a reason not to quote a < 1.5 keV color
  temperature from Dante as the wall temperature.
- **Model finding — hot laser spots**: even without the laser M-band, 10–22% of what Dante 1
  sees lies above 1.8 keV and the color temperature is 23–52 eV above Tr: the per-element
  Planckian split puts the laser spots at their brightness temperature (B/σ)^¼, far above
  the wall's, and Dante 1 sees outer spots through the LEH. With the intensity-form M-band
  on top the share reaches 0.33 at 8 ns, against Dewald's measured 12–15% (with the per-cone
  `nif_mband()` and α_M 0.5 the Dewald setup gave 12%). One conclusion was not acted on
  here: a spot's thermal emission should be split at a plasma temperature model in place of
  its brightness temperature. (Resolved in "Laser-spot spectrum (2026-10-24)" below: spots
  now thermal at the local Tr plus the M-band; the 0.33 was mostly the uncalibrated M-band
  pairing, 0.115–0.138 calibrated.)

### XRD lines of sight

`sentinel(:polar)` (Kemp 2020: polar DIM 000-000, θ 0°) views the hohlraum down its
axis through the upper LEH: on the test hohlraum (LEHs r 1.6 mm, capsule r 0.9 mm) it
sees the capsule and, around it, the lower LEH, so a blackbody cavity at 250 eV gives
the capsule's (0.9/1.6)² share (186 eV): the polar view measures the capsule, not the
wall. The equatorial DIMs (090-078, 090-315) see a hohlraum only through a window
(`through` required). The Sentinel filter set is not tabulated (Dante 1's channels
stand in, **L**). The OMEGA Dante (port H16I, θ 142.62°, φ 126°; Rubery 2023, the Omega section) depends
on the hohlraum axis chosen per shot and is left to `DanteView`.

### What remains

- Ni XRD measured shape (Cr's parameters stand in); Al below 1 keV and Cr above 1 keV
  are extrapolations; measured thin-filter transmissions (Fig. 3a's Ti is 0.45–0.83 ×
  nominal) — real filters' oxide and contamination layers would move the low channels.
- The SiO₂ mirror angle scale (published 3.5° behaves as 2.04° bulk SiO₂) applied to
  the 2.5° mirror is an assumption.
- Within each group, the voltages assume a Planckian spectral shape at Tr_open;
  finer groups near the channel edges would remove it.
- Laser-spot spectra (above; done, "Laser-spot spectrum" below) and the M-band split
  per cone in time-dependent runs.
- Dante 2's hardware (11 mirrored channels; not tabulated in Rubery 2023) and the
  Omega Dante configuration (Rubery 2023 Fig. 13).

## Laser-spot spectrum (2026-10-24) *derived*

ROADMAP "Spectral radiosity": the finding above ("hot laser spots"). `spectral_solve` split each
element's whole source as a Planckian at the brightness temperature (B/σ)^¼ of its total
radiosity. On a laser spot B = ηI + αH with ηI ≈ 1.1 × 10¹⁵ W/cm² at the outer cones'
1.5 × 10¹⁵, so a ~300 eV cavity put its spots at ~365 eV, whose Planck tail by itself contains
25% of the power above 1.8 keV (a 290 eV Planckian 12.4%, 300 eV 14.0%), and the
laser M-band was then carved out of that already hard spectrum. Read (cite
only): Dewald et al., "X-ray conversion efficiency of high-Z hohlraum wall
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
| Au spheres, direct drive, no cavity: x-ray conversion 95% at 10¹⁴, 80% at 10¹⁵ W/cm²; M-band (2–5 keV) ~1% of the x-ray flux at 10¹⁴, 27% at 10¹⁵ (LASNEX with XSN under-predicts it) | Dewald 2008 abstract | M |
| The > 1.8 keV emission in NIF hohlraums is "predominantly produced by hot, low-density NLTE gold bubble plasma" (the laser-heated corona); the M-band/total ratio is higher through the LEH, which sees the outer-beam bubble | Chen 2024 §III | H |
| A hohlraum drive spectrum = a Planckian at Tr plus a non-Planckian power law above 1.85 keV, normalized to the total; used in place of integrated hohlraum simulations and compared with them and with Dante-driven experiments | Sacks 2024 §2–3 | H |
| The recirculating M-band at the wall is the walls' Planckian thermal tail; the laser-produced part comes from the outer spots | Dewald 2020 §V.2 | H |

So the converted x-rays of a spot have two parts: the non-Planckian M-band, set by the
laser intensity (Dewald 2008's 27% is the laser-produced hard fraction itself; an
extra Planck tail at an inflated temperature counts it twice), and a soft thermal part.
The soft part's temperature is not published for spots in a cavity (Ze et al., PoP 4,
1997, OSTI 467314, compare spot and disk x-ray yields, not spectra). The brightness
temperature is a bookkeeping number of the gray radiosity, not the temperature of an
emitting surface: the converted energy is radiated by a conversion layer bathed in, and
thermalized against, the cavity field. We take that layer's thermal emission at the
local radiation temperature (**L**: the Sacks form applied per element).

### Model (`spectral_solve(...; spot)`)

- `spot = :radiation` (new default): the converted x-rays `ηP/A` (`_xray_source`, the
  number from the conversion model) less the M-band at the element's radiation temperature
  `(H/σ)^¼`; the rest of the source (a re-emitting wall's stored energy) and, with
  `thermalize`, the thermal re-emission at the wall temperature `((B − ηP/A)/σ)^¼` (on a
  wall without laser light that is (B/σ)^¼, as before); the M-band flat over its band.
  Both temperatures follow the group totals when hard groups lose energy.
- `:wall`: the converted x-rays at the wall temperature too (a lower bracket).
- `:brightness`: the earlier rule, unchanged (the Dewald scan and the N221204 rows
  below reproduce the earlier tables to the printed digit).
- Group sums are untouched: Planck fractions sum to one per element, so Σ_g S_g = S
  and every gray total of `test/spectral.jl` holds for all three rules (B, H to 1e-11,
  absorbed power to 1e-11). New tests check the source rows against the closed form of each
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
  Dante 1 then sees 11.4–12.0% above 1.8 keV against Dewald's measured 12%
  (11.9–12.5% before); the capsule M-band P2/P0 moves by +0.002 / −0.018.
- The color temperature (channels < 1.5 keV) drops from 3–6 eV above Dante's Tr to
  ~11 eV below it: the walls emit at α^¼ Tr ≈ 0.95 Tr, the spots at Tr, and 11% of
  the flux is M-band. Dewald does not publish a color temperature, so this is no
  check either way.

### N221204 (`examples/n221204.jl energetics=new unfold=on`, standard mesh 0.28/0.16)

Dante 1 over the initial LEH area (Tr; open-area Tr_open in brackets), unchanged (gray): 5
ns 243.4 (256.6), 6 ns 265.0 (283.1), 7 ns 277.9 (299.5), 8 ns 291.2 (314.5), 8.5 ns 282.4
(305.0) eV. "> 1.8 keV" is the share of what Dante 1 sees; Tc the hardware-channel color
temperature; "cones" the per-cone M-band (`nif_mband()`, each cone from its post-CBET,
post-bubble deposition, `mband_power`; the laser-to-M-band ratio is 0.125–0.194 of the
absorbed laser, the intensity form 0.040–0.130) with the calibrated α_M, thermalized.

| t (ns) | brightness | thermal | brightness + `MBandEmission()` | thermal + `MBandEmission()` | brightness, cones, α_M 0.5 | **thermal, cones, α_M 0.55** |
|---|---|---|---|---|---|---|
| > 1.8 keV, 5.0 | 0.091 | 0.060 | 0.132 | 0.103 | 0.123 | **0.116** |
| 6.0 | 0.134 | 0.097 | 0.190 | 0.156 | 0.124 | **0.118** |
| 7.0 | 0.183 | 0.127 | 0.276 | 0.229 | 0.134 | **0.125** |
| 8.0 | 0.218 | 0.155 | 0.326 | 0.274 | 0.130 | **0.122** |
| 8.5 | 0.160 | 0.133 | 0.183 | 0.158 | 0.164 | **0.158** |
| Tc − Tr_open, 5.0 | +10.3 | −11.7 | +14.4 | −9.7 | +11.9 | **−10.1** |
| 6.0 | +11.9 | −10.7 | +21.9 | −3.2 | +8.4 | **−8.9** |
| 7.0 | +23.6 | −7.5 | +46.6 | +12.7 | +9.5 | **−8.2** |
| 8.0 | +28.8 | −5.6 | +59.6 | +23.3 | +4.0 | **−12.0** |
| 8.5 | +5.0 | −9.1 | +9.6 | −4.6 | +5.5 | **−4.3** |

Rerun 2026-10-27 at 41e856f (`examples/n221204.jl 0.28 0.16 full calibrated x unfold=on`),
after the per-cone M-band switched to each cone's post-CBET deposition (decision 21) and
`MBandEmission()` gained the `:DU` liner. Before that fix the intensity-form columns were
equal to the columns without M-band on the DU wall (no element emitted); the table of
2026-10-24 predates the DU wall and the present energetics (its Tr was 248.8 eV at 5 ns).

(Tc − Tr over the initial area, the comparison of the finding: +23 to +52 eV before,
+3 to +18 thermal, without M-band.) Coarse mesh (0.8/0.4): the rule-only and calibrated
rows the same to 0.007 in share and 0.5 eV in Tc (thermal) / 2 eV (brightness); the
intensity-form rows 0.015–0.02 lower (its element-averaged intensity is mesh
dependent). Unfolds at 8 ns (input 291.2 eV): traditional 292.4 (brightness) → 291.7
(thermal), 291.7 (calibrated); smooth within 0.25 eV in every case; all 18 voltages
matched to < 5e-4.

- **The 0.33 at 8 ns was mostly the M-band pairing, not the spot split**: the example
  used the intensity-form M-band at the gray albedo in the M-band groups (α_M = α,
  F_M ≈ 7 on Dewald's setup), which no calibration supports. With the Dewald-calibrated
  pairing (per-cone M-band, α_M thermalized) both rules give 12–14%: brightness
  0.123–0.164, thermal **0.116–0.158** (the top of each range at 8.5 ns), against Dewald's 12% (Fig. 1b, at the 290 eV
  peak) to 15% (Fig. 3). His shot is a different hohlraum, so agreement to ~2 points is
  all this supports; not fitted to it.
- **What the spot split changes**: without any M-band, Dante's > 1.8 keV share falls
  from 0.09–0.22 to 0.06–0.16 (the Planck tail of the cavity: 0.07 at 250 eV, 0.14 at
  300); the color temperature goes from 5–29 eV above Tr_open to 6–12 eV below it; the
  traditional unfold's bias falls from +0.8–1.2 to +0.13–0.50 eV (its Gaussians have
  less hard power to spread). With the calibrated M-band the color temperature sits
  4–12 eV below Tr_open (4–12 above under the brightness rule).
- Neither is a Dante check: published N221204 color temperatures are not in hand, and
  the M-band share was compared, not fitted.

### What remains

- A measured spot spectrum (soft part) in a cavity: Dante or a spectrometer at a view
  dominated by spots, or a hohlraum simulation's spot emission, to replace the local-Tr
  assumption (**L**); Dewald 2008's full text (spectra vs intensity) was not reachable.
- Bubble-absorbed inner light radiating from the outer bubble (Dewald's clipping) is not
  split per cone. (The per-cone, post-CBET M-band is `drive_history`'s default since
  2026-10-26; the example's "cones" cases use it.)

## Per-cone M-band in the drive history (2026-10-27) *derived*

ROADMAP "Spectral radiosity". Until now the laser M-band entered only the `unfold=on` spectra of
`examples/n221204.jl`, with the per-cone conversions of `nif_mband()` (outer 0.20, inner
0.042; "Spectral radiosity, M-band and a synthetic Dante" §2) applied through a proxy:
each element's cone shares taken from the *pre-CBET* deposition of `beams_at(t)`, times
its post-CBET absorbed laser power.

**Change.** `mband_power(m, op, beams; lpi, time)` deposits each cone's *effective* beams
(after CBET, backscatter and the bubble) through `lpi_deposit`, so light the bubble's
ridge stop puts off the operator converts at the fraction of its cone. `drive_history(...;
mband, spectral)` does this at every step with that step's beams, splits the step's gray
solution by `spectral_solve` (default edges [0, 1.8 keV, ∞]) and records `mband_moments`
(Legendre moments of the observed flux above 1.8 keV) and `mband` (laser-produced M-band
power per part). (Off by default when added; on by default since, below.) It is a diagnostic: the M-band is part of the converted
x-rays and, thermalized, the group sums are the gray solution's, so the gray history does
not depend on it.

**Checks.** (`test/spectral.jl`, "drive history M-band", coarse 10.1 × 5.75 mm hohlraum,
constant-albedo gold, a transfer from outer to inner growing in time.) The gray records with
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

### On by default, and the unfold's cones (2026-10-26, decisions 21–22)

**Cost first** (warmed, minimum of ≥ 4 runs, 8 threads pinned to one CCX; N221204 = the
full model `0.8 0.4 full calibrated table`, 1088 elements, 8 steps; test-sized = the
`test/spectral.jl` hohlraum, 612 elements, 4 steps, CBET):

| `drive_history` | M-band off | on, as added | on, now |
|---|---|---|---|
| N221204 0.8/0.4, wall time | 18.3–18.5 s | +2.6 s (+12%) | ±0.2 s (below the ±0.5 s run-to-run spread) |
| N221204 0.8/0.4, allocated | 11 051 MiB | +62 MiB | +31 MiB (+0.3%) |
| test-sized, constant albedo | 0.334 s, 8.7 MiB | 0.724 s (+117%) | 0.354 s (+6%), 17.5 MiB |
| test-sized, wall loss | 0.357 s, 13.7 MiB | 0.837 s (+119%) | 0.378 s (+6%), 22.4 MiB |

As added, the per-cone M-band deposited every cone again with the other beams at zero
power: four full depositions per step (0.27 s of a 2.3 s N221204 step; on the small case
a deposition is most of a step). Two changes took it out:

- **The cones split in the deposition pass of the step.** `drive_history` deposits once
  with a per-group split (`_deposit(op, beams, split)`; the bubble's ridge stop through
  `_lpi_deposit`): each lit point also adds its value to its beam group's column. The
  deposition itself is exactly that of `deposit` (the same sums in the same order; the
  history's gray records `==` with and without the M-band), and each column is exactly
  the group deposited alone, `==` (zero-power beams add exact zeros; checked on the test
  hohlraum with CBET and on N221204 with the bubble's ridge stop at 8 ns, ridge fraction
  1). The split pass costs what a plain one does (0.061 / 0.074 s on N221204 0.8/0.4 /
  0.28/0.16). `mband_power(m, op, beams; lpi)` uses the same pass (from 0.27 to 0.06 s).
- **The two-group thermalized solve skips its gray-albedo group**: thermalized, the group
  sums are the gray radiosity, so the one group at the gray albedo is `B_gray − B_M`
  in place of another gray-type solve (a factorization or GMRES; 0.42 → 0.31 s per step on
  the standard mesh). The M-band group is untouched (N221204's M-band P2/P0, share and
  laser M-band columns are identical to the printed digits); the soft group matches the
  dense solve of the coupled system to 1e-9 (`test/spectral.jl`).

What remains is the M-band group's thermalized GMRES on the dense F: 0.015 s per step at
0.8/0.4 and 0.31 s at the standard mesh 0.28/0.16 (5528 elements, a gray solve there 0.39
s; a warmed step there takes ~10 s, so ~3%, an estimate from the parts: the whole history
is not yet timed on and off at that mesh).

**Default.** `drive_history`'s `mband` is `MBandSplit()`: `nif_mband()` (Dewald's
per-cone conversions) in the calibrated pairing, edges [0, 1.8 keV, ∞], `thermalize =
true`, α_M 0.55 (`mband_albedo(t, emission, edges)` on the emitting surfaces at and above
1.8 keV; the calibration above, "Laser-spot spectrum"). The split's keywords were first
`drive_history`'s `spectral` keywords, defaulted when unset; `MBandSplit` carries them
with the model since 2026-10-27, and `spectral` holds `spectral_solve`'s other keywords.
A plain call thus gives the per-cone, post-CBET M-band with the calibrated α_M; the gray
default before (α_M = α, no thermalization) gave F_M ≈ 7 on Dewald's setup, which no
calibration supports. `nif_mband()`'s materials are now `(:Au, :DU)`: the DU wall of a
NIF hohlraum is Au-lined and its liner is what the laser heats (the N221204 example
already passed `materials = (:Au, :DU)`). The M-band is on only where its inputs exist: a
target with a surface of the model's materials, lit (per-group model) by a beam of
`beams_at(times[1])` in its groups; otherwise `mband_moments` and `mband` are `nothing`
(e.g. an Omega target, whose beam groups `nif_mband()` does not convert; no calibrated
Omega M-band exists). `mband = nothing` turns it off. The example's `mband=on` option
(the per-cone M-band beside the pre-CBET proxy) is removed: the full model prints its
M-band P2/P0 at every run, and the comparison is the table above (commit d19d84f).

**The unfold's "cones" cases** (`examples/n221204.jl unfold=on`, decision 21) now use
`mband_power` with the step's effective beams and `full_lpi` (CBET + bubble) instead of
the pre-CBET cone weighting. Change on the coarse model (0.8/0.4, default times), the
calibrated rows (brightness α_M 0.5 / thermal α_M 0.55):

| t (ns) | Dante > 1.8 keV share | color temperature Tc | traditional unfold bias | laser → M-band |
|---|---|---|---|---|
| 5.0 | −0.0020 / −0.0021 | −0.26 / −0.08 eV | −0.008 / −0.009 eV | 0.1284 → 0.1255 |
| 6.0 | −0.0019 / −0.0020 | −0.34 / −0.22 eV | −0.006 / −0.006 eV | 0.1470 → 0.1442 |
| 7.0 | −0.0013 / −0.0014 | −0.29 / −0.22 eV | −0.004 / −0.003 eV | 0.1825 → 0.1804 |
| 8.0 | −0.0000 / −0.0000 | −0.01 / −0.00 eV | 0 | 0.1934 → 0.1935 |
| 8.5 | +0.0002 / +0.0002 | +0.04 / +0.04 eV | 0 | 0.1930 → 0.1937 |

So |Δ| ≤ 0.0021 in the share, ≤ 0.34 eV in Tc and ≤ 0.003 in the ratio of M-band to laser
power, largest at 5–7 ns, where CBET moves power between the cones; the gray Tr, Tr_open and
the CSV (P2, P4, Tr, Dante) are byte-identical. The standard-mesh table of "Laser-spot
spectrum" (N221204 section) was rerun after the switch (2026-10-27).

## Dante-family hardware and geometry: what is public (2026-10-29)

For the diagnostic definitions (`reference/diagnostics-design.md`): what a synthetic
Dante, Sentinel, Omega Dante or Z XRD can take from public sources, and what a user
supplies. Gathered by a literature pass and then checked by a second, independent pass
against the PDFs (tables read from page images); the field of view, standoffs, Chen's
flux pair, Kemp's solid angles and grids, the MLM channel's geometry and the BAH
occlusion were checked again against the text. Papers (OSTI IDs): Dewald et al., RSI 75,
3759 (2004), 15014120; Campbell et al., RSI 75 (2004), 15014117; May et al., RSI 81,
10E505 (2010), 1009652; Fournier et al. (2010), 1012430; Glenzer et al., PRL 104 (2010),
1119940; May et al. (2012), 1226204; May et al. (2015), 1234624; Marrs et al. (2015),
1366953; Kemp et al., RSI 91, 083507 (2020), 1736341 (Dante 1 configurations); Kemp et
al. (2020), 1736345 (Sentinel); Opachich et al. (2021), 1736342; Harris et al., RSI 92,
033505 (2021), 1820019; Rubery et al. (2022), 1959490 (the multilayer channel); Rubery
et al., RSI (2023), 1975922; Rubery et al., PRL 132, 065104 (2024); Chen et al. (2024),
2429360; arXiv 2501.10350 (the BAH platform); Schoenzeit, 3363346; Jones et al.,
SAND2014-4261J, 1145734 (Z); Fournier (2007), 920482 (Omega); NIF User Guide (2024-09).
Not obtained: Sorce et al., RSI 77, 10E518 (2006) (Omega Dante; no full text); Girard
et al., RSI 80 (2009), OSTI 22053575 (abstract only; first pass misattributed it to
Kline).

### Corrections to the earlier sections

- **Grids are tantalum.** 25 µm Ta, 8.8% open (Kemp 2020 Tables I–III; Fournier 2010;
  May 2015); only Rubery 2023 §II.A states tungsten ("Spectral radiosity, M-band and a
  synthetic Dante", "Dante channel responses and the unfold": "W grids"). The model
  treats a grid as opaque wires (transmission = open fraction), so nothing changes for
  the current channels; a grid on a channel above ~10 keV would need the Ta transmission.
- **Mirrored channels changed over time**, not a contradiction: Dante 1 had 4 in 2009,
  5 in 2012–15, 4 in 2017–18 and 3 from October 2021 (Rubery 2023); Dante 2 8 in 2015
  (User Guide, from a 2015 LLNL news item), 11 from 2021 (Rubery 2023). A channel set
  belongs to a date or shot.
- **Dante 1's port** is θ 143.25°, φ 274° (User Guide Table 5-1; 36.75° to the axis; the
  BAH paper uses 36.75°); Rubery, Kline, Marrs and Fournier give 37°, Chen 2024 37.4°.
  Since 2026-11-03 every Dante 1 definition names the port and takes the table's
  143.25° ("Named ports and per-channel fields (2026-11-03)" below, with the change
  measured). Dante 2: P64-350, θ 64°, φ 350° in the table too.
- **Tr's area is an analysis choice.** Rubery 2024 and the BAH paper (its eq. 12) use the
  initial LEH area; Glenzer 2010 used A_LEH(t) from the 3–5 keV SXI images (diameter to
  83%) and a 10–15 eV view-factor correction to the capsule; Fournier 2010 used 19.65 /
  19.89 mm² (Dante 1 / 2), the projected area of the cylinder wall and end cap as each
  Dante saw them.

### 1. Field of view, distances, apertures

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| NIF Dante standoff | 8.5 m average, collimators at 5.5 m (2004 design) | Dewald 2004 text, Fig. 2 | H |
| NIF Dante standoff (fielded) | 8.8 m | Fournier 2010 §II | H |
| Filter wheel; MLM filter | 8 m; 8.9 m | Rubery 2022 | H |
| Field of view | ±30 mm at TCC, set by collimators (2004 design) | Dewald 2004 text | H |
| Field of view (fielded) | "centimeter-wide"; it took in the target support lit by unconverted light, raising the foot Tr in the first two channels | Girard 2009 abstract | M |
| What a channel sees | the whole target, not the chamber | Marrs 2015 Fig. 1; Harris 2021 Fig. 1 | H |
| Channel solid angle | set by an aperture in front of each XRD | Rubery 2023 §II.A | H |
| Dante 1 Ω (2017–18) | 0.6514 µsr unmirrored; ch 18 1.773 µsr; mirrored 0.5558–0.5929 µsr | Kemp 2020 Tables I–III | H |
| Implied XRD apertures | 8.0 mm and 13.2 mm at 8.78 m give 0.6514 and 1.773 µsr | derived | L |
| MLM channel (ch 2) | 3.2 × 14 mm elliptical aperture at the mirror box, 2.6 × 8 mm at the XRD; FOV 12 mm diameter; Ω 1.91e-7 sr (GEANT4); acceptance ±0.05° | Rubery 2022 | H |
| Omega Dante port | H16I, θ 142.62°, φ 126° | Rubery 2023 §IV | H |
| Omega view angle | 37.37° to a P6–P7 or P5–P8 axis (37.3°, OSTI 920482) | May 2012 | H |
| Omega XRD apertures | ch 1 10.5 mm, ch 2–4 8.3 mm, ch 5–18 13.2 mm | Rubery 2023 Fig. 13 | H |
| Sentinel | ≈ 7 m (can be closer); entrance 12.9 mm² (1.45 × 8.89 mm), 18.1 mm² active, plates at 45° to the axis | Kemp 2020 Sentinel §II | H |
| Z LOS50 | 12° polar, 50° azimuth: XRD 23.20 m (11.6°), bolometer 18.84 m, PCD 18.36 m; sees the bottom 75% of a 10 mm array | Jones 2014 Tables 1–2 | H |

These are not public: the NIF collimators and field of view as fielded (only the 2004 design
and the MLM channel's 12 mm are), the large-aperture channels' geometry, Omega Dante's
standoff and solid angles, Dante 2's solid angles, current grid transmissions.

As a consequence, a NIF Dante sees the whole target, so the synthetic view is limited by
occlusion and the target's extent, not by the collimator, except the MLM channel (12 mm)
and large assemblies (a halfraum pair, a package off TCC).

### 2. Channel sets

Notation: cathode / filters (µm; "g" 8.8% Ta grid) / mirror / Ω (µsr; 0.6514 where
not given).

**NIF Dante 1, N170924-003** (Kemp 2020 Table I): 1 Al / Sc 0.8 + Ti 0.5 + g / SiO₂
3.5° / 0.5786; 2 Al / B 1.0 + CH 0.8 + g / a-C 5° / 0.5929; 3 Al / Lexan 4.0 + g / SiO₂
3.5° / 0.5786; 4 Ni / V 2.0 / SiO₂ 2.5° / 0.5603; 5 Cr / Co 1.4 + CH 0.2 + g; 6 Cr / Cu
1.1 + g; 7 Ni / Mg 5 + Zn 1.0 + CH 0.1 + g; 8 Ni / Mg 22 + g; 9 Ni / Al 20; 10 Al / Si 20;
11 Al / Ti 26; 12 Al / Saran 66; 13 Al / Ag 4.4; 14 Al / Lexan 2 + CaF₂ 21 + CH 3; 15 Ni /
Fe 30; 16 Ni / Ti 26; 17 Cr / Mn 26.4 + Ni 3.6; 18 Ni / Mn 26.4 + Ni 3.6 / 1.773.

**N180129-001** (Table II; mirrors and cathodes as Table I; Ω ch 1–4 0.5740, 0.5881,
0.5740, 0.5558): 1 Sc 0.8 + Ti 0.5; 2 B 1.5 + CH 1.2; 3 Lexan 4 + g; 4 Ni / V 1.75; 5 Co
1.4 + CH 0.2 + g; 6 Cu 2.0 + g; 7 Mg 5 + Zn 0.65 + CH 0.1 + g; 8 Mg 22 + g; 9 Al 20 + g;
10 Si 20; 11 Al / Fe 1 + Cr 0.65 + parylene 5 + g; 12 Saran 66; 13 Ag 4.4; 14 Al / Ti 13;
15 Ni / Zn 25; 16 Ni / Al 250; 17 Cr / Ni 10; 18 Ni / Fe 1 + Cr 0.65 + parylene 5 / 1.773.

**N180604-002** (Table III; Ω ch 1–4 0.5740, 0.5929, 0.5786, 0.5603): 1 Sc 0.8 + Ti 0.5
+ g; 2 B 1.0 + CH 0.8 + g; 3 Lexan 4 + g (grid printed twice); 4 V 2.0 + g; 5 Co 1.4 + CH
0.2 + g; 6 Cr / Al 250; 7 Ni / Mo 32; 8 Ni / Y 50 + Al 10; 9 Ni / Ge 42 (a "C/375" entry
as printed, unresolved); 10 Si 20; 11 Al / Mg 16.4 + Ni 3.6; 12 Saran 66; 13 Ag 4.4; 14
Al / Mg 22 + g; 15 Ni / Al 20; 16 Ni / Mg 5 + Zn 1.0 + CH 0.1 + g; 17 Cr / Cu 1.1 + g; 18
Ni / Al 250 / 1.773.

**From October 2021** (Rubery 2023 Fig. 4): `_DANTE1_HARDWARE`, checked channel by
channel against the figure (the one difference, ch 2's "47 µm CH" as polyimide, follows
Rubery 2022 and 2023 §II.D). Unresolved: Fig. 4 gives ch 16 and 18 both 1.77 µsr where
every earlier table has the large aperture on ch 18 only; Rubery's text calls ch 17 a
duplicate of ch 5, which the figure does not show.

**N091104 (2009)** (Fournier 2010 Table II; per-foil thickness × count, the count
multiplying as May 2015 Table IV shows; Ω 0.6490 on ch 5–17): 1 Cr / Al 0.75 × 2 / B₄C
7° / 0.598; 2 Al / (B 0.2 + Lexan 0.4) × 4 / a-C 5° / 0.588; 3 Al / Lexan 2 × 2 / SiO₂
3.5° / 0.574; 4 Ni / V 1 × 2 / SiO₂ 2.5° / 0.556; 5 Cr / Ni 12 × 2; 6 Cr / Cu 0.6 × 2 +
grid (T 0.091); 7 Ni / (Zn 0.65 + CH 0.1) × 2 + grid; 8 Ni / Mg 5.5 × 3 + grid; 9 Ni / Al
5 × 3 + grid; 10 Al / Si 11 × 2; 11 Al / (CaF₂ 7 + parylene 1) × 3; 12 Al / Saran 25 × 3;
13 Al / Ag 2.5 × 2; 14 Al / Ti 11 × 2; 15 Ni / Fe 15 × 2; 16 Ni / Zn 25; 17 Ni / Al 250;
18 Ni / Fe 1 + Cr 0.65 + parylene 5 + grid / 1.773.

**2012–14** (May 2015 Table IV, "typical" for N120627/28, N140220/21): five mirrors (ch
1 B₄C 7°, 2 a-C 5°, 3 SiO₂ 3.5°, 4 and 17 SiO₂ 2.5°); Ω ch 1 0.6026, 2 0.5929, 3 0.5786,
4 and 17 0.5603, 18 1.773, others 0.6514.

**Omega Dante, current** (Rubery 2023 Fig. 13; ch 1–2 offline; no solid angles):
cathode / filters / mirror / aperture mm, energy keV — 1 Cr / Al 1.5 / C 7° / 10.5, 0.06;
2 Al / B 0.4 + Lexan 0.8 / C 5° / 8.3, 0.17; 3 Al / Lexan 4 / Be 3.5° / 8.3, 0.27; 4 Al /
V 2 / C 2.5° / 8.3, 0.49; then 13.2 mm, no mirror: 5 Cr / Co 1.4 + Lexan 0.2, 0.79; 6 Cr /
Cu 1.2 + grid, 0.94; 7 Ni / Zn 1.3 + CH 0.2 + grid, 1.02; 8 Ni / Mg 11, 1.25; 9 Al / Al
10, 1.40; 10 Al / Si 20, 1.70; 11 Al / Fe 1 + Cr 0.65 + parylene 5, 3.50; 12 Al / Saran
75, 2.15; 13 Al / Ag 5, 2.92; 14 Ni / Ti 22, 4.00; 15 Ni / Mn 30, 5.30; 16 Ni / Fe 30,
6.15; 17 Ni / Zn 50, 8.40; 18 Ni / Al 250, 15.0. Channel counts disagree across sources
(13 in 2006, OSTI 920482; 15, May 2012; 18 on 15-position hardware, Rubery 2023).

**Omega, 2003–04** (Campbell 2004 Table 1; 12 channels, no thicknesses): bands, cathodes,
mirrors (C 2.5°, 5°, 7°; Be 3.5°) and filter elements only.

**Sentinel, N170924-003** (Kemp 2020 Sentinel Table II; filters chosen to emulate that
shot's Dante 2 without its mirror channels; "s" 35% sieve): 1 Ni / Mg 22; 2 Ni / Mn 13.2
+ Ni 1.8; 3 Cr / Al 5 + s; 4 Cr / B 5 + CH 1 + s; 5 Al / Lexan 2 + CaF₂ 14 + CH 18; 6 Cr
/ parylene 10 + s; 7 Al / Sc 1.6 + Ti 1 + s; 8 Al / Ti 11; 9 Ni / V 2; 10 Al / Fe 15; 11
Cr / Co 1.4 + CH 0.2; 12 Cr / Si 20; 13 Ni / Al 20; 14 Al / Saran 44; 15 Ni / Ti 26; 16 Al
/ Ag 2.2. Relative responsivities ±10–20%, absolute calibration crude.

**Dante 2:** no channel table is public (16 channels like Dante 1 in 2009, Fournier
2010; "like Dante 1 but missing several 3–5 keV channels" in 2012–14, May 2015;
normalized responses only in Opachich 2021 Fig. 5). Users supply it.

**Z** (Jones 2014): five XRD filters (µm): Kimfol 4; V 1; Zn 0.75 + parylene N 0.5; Be 8
+ parylene N 1; Be 10 + V 0.8.

### 3. Power to Tr; other surfaces; published synthetic views

| Item | Value | Source | Conf. |
|---|---|---|---|
| Relation | F = A cos θ σ Tr⁴; V_i = P Ω_i ∫ R_i S dE | May 2010 eq. 1; Harris 2021 eq. 1 | H |
| Flux band | to ≈ 13 keV (NIF); 0–2 or 0–3 keV (Omega) | Harris 2021; May 2010 | H |
| Area uncertainty | excluded from published error bars; emission outside the LEH is recorded | May 2010 §IV | H |
| Thin walls | A_eff(E) = A_LEH cos θ + T_wall(E) A_wall sin θ; halves an Fe K-shell yield (Ti ÷ 1.5) | May 2012 eq. 2 | H |
| Lambertian LEHs | Dante 1 / Dante 2 below 3 keV 1.8 = cos 37° / cos 64° (1.82) | May 2015 | H |
| Simulations through the instrument | LASNEX spectra through the channel responses and the same unfold as data | May 2015 | H |
| Synthetic Dante by view factors | VisRad with the Dante 2 line of sight occluded by an annular window-retaining ring: irradiance on a small element 6 m out, 27% blocked (≈ 20% by a simple estimate); I uncertainty 5% (Dante 1), 10% (Dante 2) | arXiv 2501.10350 §V.B | H |
| Side views | Dante 2 through a window onto the back wall, two LEH shields; in/out albedo correction; LASNEX synthetic Dante, method unpublished | Opachich 2021 | H |
| Obstructions | a gas-pipe stalk blocks Dante 2 | Schoenzeit | H |
| Orientation | through the LEH ≈ 15 TW/sr peak; through an open end ≈ 35 TW/sr (a larger, cooler area; higher M-band share through the LEH) | Chen 2024 pp. 3–4, Figs. 1–2 | H |
| Keeping clear of the views | a 10 µm polyethylene filter ≈ 5 cm from the target, outside both Dante views | Kemp 2020 Sentinel §III | H |

Checks for the field-of-view path (never fit targets): Chen 2024's LEH / open-end pair;
the BAH 27% occlusion; May 2015's cos 37° / cos 64° ratio.

### 4. Timing

| Item | Value | Source | Conf. |
|---|---|---|---|
| XRD-31 impulse response | ≈ 140 ps FWHM at 5 kV | Rubery 2023 Fig. 5 | H |
| NIF bandwidth | 2.5 GHz (54 m LMR-600); Wiener deconvolution; attenuation 0–62 dB + 19.1 dB fixed (Fig. 2 labels 12 dB) | Rubery 2023 §II.A | H |
| Digitizers | 10-bit Keysight, 10 GS/s: Dante 1 from summer 2016, Dante 2 from June 2017 (SCD5000 before) | LLNL news 2016-11-15; Opachich 2021 | M; H |
| Effective response; timing | ≈ 400 ps Gaussian; absolute 50–75 ps (±50 ps) | Rubery 2024; Rubery 2023 | H |
| Omega | 16 SCD5000, 5 GHz; 100–200 ps | Rubery 2023 §IV; May 2010 | H |
| Sentinel MiniXRD | ≈ 250 ps at 1.5 kV | Kemp 2020 Sentinel | H |

### 5. Uncertainty

| Item | Value | Source | Conf. |
|---|---|---|---|
| NIF at peak (1σ) | flux 4%, Tr 1.5% | Rubery 2023 §II.C | H |
| Dante 1 per channel | 8.4–16.4%; MLM 8.7% | Rubery 2023 Fig. 4, §II.D | H |
| Component calibration (2009) | 10–20% per channel | Fournier 2010 | H |
| Omega per channel | ≈ 20% at 200 eV to ≈ 5% at 2 keV; ΔF, ΔTr fits eqs. 3–6 | May 2010 | H |
| Omega 2004 | systematic 17.4% (60 eV) to 2.3%; Tr 2.7 / 1.6 / 1.2 / 0.7% at 50 / 100 / 200 / 300 eV | Campbell 2004 Table 2 | H |
| Opacity campaign | flux ±5%; Tr 2.1% at 100 eV, 1.5% at 200 eV | Opachich 2021 | H |
| Unfold with channels lost | 1–3% typical; 14 / 23 / 15% with ch 1–2, 6–9, 11 and 18 off (summer 2019) | Harris 2021 | H |
| XRD aging | response falls over years (C, O contamination) | Campbell 2004 Fig. 1 | H |

## X-ray imagers, radiography and streaked diagnostics: what is public (2026-10-29)

For the next diagnostic kinds (`reference/diagnostics-design.md`, phase 5). Gathered by a
literature pass and checked by a second, independent pass against the PDFs (tables from
page images); Gregor's Table II, Izumi 2018's spot-temperature statement and the
P18-123 port were checked again against the text. "D" design, "F" as fielded, "S" one
shot or campaign.

### Corrections to the first pass

- Izumi 2018 does **not** state that 3–5 keV hohlraum images are dominated by laser
  spots: its thin-wall camera images 7–12 keV, a photon energy 3–5 × Te (Te 2.2–2.7 keV
  simulated at 2–3 × 10¹⁴ W/cm²), so that signal rises faster than the fourth power of
  the laser intensity. Schneider 2010 §I and Chen 2020 say only that the SXI images show
  beam spots and expanding plasma. A thermal (Tr) model's adequacy for the 3–5 keV
  channel is open, to be measured, not assumed.
- The SXI soft channel is 870 eV in Schneider 2010/2012 and 900 eV in the User Guide.
- G-LEH-1's specifications are in Chen 2021 Table I, not Chen 2020; its filters are Ti
  5–25 µm (Chen 2018), not "Ti 10 or Al 30 µm".
- The 2015 K-mirror and SOP upgrade paper is Manuel et al. (SPIE 2015), OSTI 1240931.
- Omega ASBO is F/3.3, 15.3× to the interferometers, 25.5× to the slit.
- Lanier et al. is RSI 87 (2016).

### 1. Ports (NIF User Guide Tables 5-1, 6-1, App. E; H)

(θ, φ) in degrees: P18-123 (18, 123.75), the SXI-U port; P161-326 (161, 326.25), SXI-L
and G-LEH (Chen 2020 §II.1, Chen 2021 §II); DIMs P0-0 (0, 0), P90-315 (90, 315), P90-78
(90, 78.75); P90-123 (90, 123.75) and P90-348 (90, 348.75); P90-89, P90-100. VISAR and SOP
share the 90-315 DIM. SPIDER at (7, 90), through the upper LEH. That SXI-U's line of
sight is its port axis is inferred (M) from SXI-L's.

### 2. X-ray imagers

| Instrument | Values | Source | Conf. |
|---|---|---|---|
| SXI hard channel (F) | M = 2; detector 4.97 m (2010), 4.88 m (2012, upper); four 25 µm pinholes; 75 µm Be + 50 µm Ti, 3–5 keV; ≈ 60 µm | Schneider 2010 §II (OSTI 1009651); Schneider 2012 §II (1043636, a draft) | H |
| SXI detector | back-thinned Si CCD, 15 µm Si, 24 µm pixels, 2k × 2k, absolutely calibrated | Schneider 2010 §II | H |
| SXI soft channel | W/B₄C multilayer, 30 bilayers, 43.1 Å period (W 21.55 Å), 40 Å SiC cap, 10°, 0.7 mm Si substrate; 45 µm pinhole; Cu 2 µm + polyimide 0.5 µm (2010) / 1 µm (2012) | Schneider 2010, 2012 | H |
| G-LEH-1 (F) | Furi hCMOS 1024 × 448, 25 µm, 2 frames, ≥ 2 ns gate; 22 pinholes at 0.9–1.2×, up to 8 images by timing skew; 56 µm at 2×, 87 µm at 1×; Ti 5–25 µm | Chen 2021 Table I (1860726); Chen 2018 §III (1524285) | H |
| G-LEH-2 (F, from April 2020) | two Icarus sensors 1024 × 512, 25 µm, 4 frames each; gates 1, 2 or 38 ns, fastest "1 on / 2 off" (FWHM 1.4 / 1.0 / 1.4 / 1.0 ns); 0.7–10 keV; M 0.6–4×, at 1× pinholes 3.25 m and camera 6.5 m; 25 µm pinholes, ≈ 70 µm; Ti 5–40 µm; light-tight 0.5 µm polyimide + 1 µm Al + 0.5 µm polycarbonate; the halves' relative frame timing 2.86 / 2.60 / 6.84 / 6.60 / 10.90 / 10.65 / 14.98 / 14.74 ns | Chen 2021 Table I, text | H |
| G-LEH field | ≈ 80° in azimuth, a little over half the hohlraum length (N171227-001) | Chen 2020 §II.1 (1632379) | H |
| Gated response | S = G ∫ g(t − t₀(x)) I(t) dt per pixel, t₀ varying across the sensor | Chen 2018 §I | H |
| GXD (D) | Au photocathode, 30 mm, four 7.5 mm strips; gate 70–1500 ps; strip delays to 44 ns; ≈ 50 µm; P43; 4k × 4k CCD | Oertel, UCRL-TR-224211 (896014) | H |
| Camera standard | one, two or four strips, up to 50 ns apart; Kodak 16801E 4096², 9 µm | Kimbrough (1124976) | H |
| Implosion GXD | 10 µm pinholes in 70 µm Ta, 200 µm pitch, M = 15, array 80 mm from TCC; yield limits CCD ≈ 10¹³, film ≈ 10¹⁵ | Bell (1124981) | H |
| Detector MTF | ≈ 0 near 15 mm⁻¹ (detector plane) | Trosseille 2018 (1479070) | H |
| Izumi's three cameras (S: N171010-001 … N181010-001; 5.75 × 10.13 mm, LEH 3.37 mm; ≈ 100 ps gates, 12 frames) | axis: 2×, 100 µm × 24, 380 mm to detector, Fe 12 µm + Kapton 650 µm (< 7.1 and ≈ 14 keV); φ 315°: 0.75×, 80 µm × 16, 250 mm, Kapton 1.15 mm (7–12 keV); φ 78.75°: 12×, 10 µm × 186, 1200 mm, Au 2.5 µm + HDC 162 µm + polycarbonate 1.5 mm + Kapton 1.025 mm (≈ 15 keV) | Izumi, RSI 89, 10K111 (2018) Table I | H |
| Omega XRFC (S: CEA 2013, 2015) | 4-strip MCP, Au, 4096²; 4 images per strip 50 ps apart, strips 200 ps apart | Poujade, arXiv 1705.10147 §IV | H |
| Omega SLOS-TRXI | port H4F (45.2°, 234°); 10 µm pinholes at 166.5 mm; cathode 2238 mm, M 12.4; 4–9 keV; Icarus, 2 ns; ≈ 40 ps | Theobald (1477184) | H |

### 3. Radiography

| Platform | Values | Source | Conf. |
|---|---|---|---|
| NIF in-flight radiography (S) | Ge foil 12 mm from the capsule, equatorial; 8 beams, 52 kJ; Ge He-α 10.25 keV; ≈ 1 mm source, ≈ 2 ns; hohlraum holes 0.9 × 0.83 mm, HDC windows; GXD 8×, two strips, ≈ 90 ps, 20–25 µm pinholes | Rygg 2014 (1184746) | H |
| NIF area foils | 3 × 3 or 1 × 1 mm, 10–15 µm; 4–8 beams, 16–60 kJ at 40° to the normal; 2–5 ns after a 1 ns picket 3 ns earlier | Barrios (1162253) | H |
| ARC | 4 beamlets of Q35T (8 planned); 1 kJ at 38 ps … 0.25 kJ at ≈ 1 ps; B353 delayable 3.6 ns after B354; pointing ±50, ±50, +10/−45 mm | User Guide §4.10 | H |
| ARC dual backlighter (no shot numbers) | two Au wires (25 or 10 µm) tipped 3–6° to the line of sight, two 900 J / 30 ps beamlets each; wires 7 and 8 mm from the target, ≈ 5° apart; MS image plate ≈ 600 mm, M ≈ 75 (87 and 76 recomputed); Cu 500 + Al 500 µm; frames at bang time −50, −10, +160 ps, ±40 ps | Tommasini, PRL 125, 155003 (2020) (1718992) | H |
| AXIS | two gated thick-MCP imagers at DIMs 90-78 and 90-124; 40–200 keV; 8° pore bias | User Guide Table 6-1; Hall 2016 (1258550) | H |
| ARC Dy flag | 300 × 300 × 5 µm edge-on, ≈ 2.5 kJ / 10 ps, ≈ 100 keV effective; object 41 mm, plates 500 mm beyond, M 13.2; 10 Fuji SR | Hill 2021 (1811762) | H |
| Omega Compton | Au wire 10 µm on a 300 µm square; plate ≈ 400 mm; ≈ 10 µm, ≈ 10 ps | Tommasini 2011 (1016984) | H |
| Omega streaked radiography (S) | V 5 µm foil (5.2 keV), 8 beams ≈ 320 J; slots 1000 × 300 / 1000 × 100 µm; slit M 20; Au cathode, 4 ns, ≈ 40 ps; P6–P7 hohlraum, H7–H14 axis | Hicks 2010 (1113950) | H |
| Image plates (NIF) | Fuji MS/SR/TR; FLA-7000 at 25/50/100 µm; S1272/S4000/S10000 = 500/588/678 V; fade correction | Holder (1479074) | H |
| Image-plate response | PSL = (1/W) ∫ e^(−z/L) (dE/dz) dz, SR W 2.1–2.3 MeV, L 40–50 µm; mPSL/photon at S1000/S4000/S10000 — SR 1.80/2.00/2.10 (5.9 keV), 2.03/2.33/2.44 (22.4), 0.72/1.14/1.29 (59.5); MS 4.24/4.50/4.72, 7.91/8.42/8.82, 4.73/5.61/6.07; S labels are scanner-specific; MS fades 58% in 105 min; ±8.7% | Rosenberg, RSI 90, 013506 (2019) (1497010) | H |

### 4. Streaked optical and x-ray

| Instrument | Values | Source | Conf. |
|---|---|---|---|
| VISAR fringe constant | τ = 2h(n − 1/n)/c; VPF₀ = λ/[2τ(1 + δ)], δ 0.0318 (fused silica, 532 nm); window and in-dielectric corrections; 2 mm fused silica 24.9 km/s, 30 mm 1.66 (recomputed) | Celliers 2004 eqs. 1–5 (926409) | H |
| NIF VISAR (D) | 90-315 DIM; probe 659.5 nm; f/3 fused-silica triplet ≥ 500 mm from TCC; objects 1–5 mm; two interferometers; Dove prism, later K-mirror; 19 mm slit | Malone 2004 (15014809, 15014425); Manuel 2015 (1240931) | H |
| NIF VISAR fields | 1 and 2 mm; 660 nm with a VPF correction | Celliers & Millot, RSI 94, 011101 (2023) (2005084) | H |
| NIF VISAR performance | < 0.35 km/s; < 50 ps; keyhole field ≈ 350 µm | Celliers 2011 (1035299); Robey 2013 (1113420) | H |
| NIF SOP | 540–645 nm, f/5, 2 mm (2005, D); 560–640 nm, ≈ 12 µm, gated cathode (2015) | Malone 2005 (15016768); Manuel 2015 | H |
| NIF optical streak camera | P510 tube, P43; windows 2–90 ns, four per camera; 4096², 9 µm, binned 3 × 3; one calibration: 3/6/16/40 ns windows, 19.1/28.5/54.7/159.3 ps | Datte 2013 slides (1092512) | H |
| Omega ASBO | 532 nm; F/3.3, 146 mm; 15.3× and 25.5×; ≈ 760 µm field; ≈ 4 µm; etalons 2–30 mm (25–1.7 km/s per fringe) | Celliers 2004 | H |
| Omega SOP (from 2011) | OG590 to S20 (≈ 590–850 nm); P510; windows 6/17/46/96 ns; slit ≤ 80 µm apparent (≈ 20 / 50 ps at 6 / 17 ns); M ≈ 21; 800 µm slit recommended; ≥ 0.5 eV | Gregor, RSI 87, 114903 (2016) (1341873, accepted manuscript) | H |
| Omega SOP calibration | T = T₀ / ln(1 + A/I), A = A₀ X(WS₀)/η; ND 0 → T₀ 1.909 eV, A₀ 481,000 ADU/ns … ND 2.5 → 1.789, 3016 (Table II, 15 rows); for shots ≥ 74528; A₀ × 0.434 / 0.831 / 0.886 for shots 68276–70551 / 70552–71572 / 71573–74527; ±3.5% (A₀), ±0.2% (T₀) | same, eqs. 7–8, Tables I–II | H |
| NIF DISC | 263 nm fiducials, 45 ps, ≈ 300 ps apart; 500 Å Al fiducial cathode | Opachich 2012 (1044234) | H |
| X-ray streak | windows 1–100 ns; a hardened camera records 680 ps at ≈ 30 ps, CsI, cathode 1.2 m from target | Kimbrough; MacPhee (1891708) | H |
| VISAR on hohlraum witnesses | Al mirrors tamped with Be on the axis (−9° to +13°) and side-on (+37° to +47°) | Lanier, RSI 87 (2016) (1394969) | H |

### Not public, or not found

NIF VISAR sweep windows and etalon sets; NIF SOP's absolute calibration; Omega TIM port
assignments (PI portal); ARC radiography shot numbers; the AXIS 128× / 300 µm figure
(news item only); G-LEH-1's original filter set; fused silica's δ at 659.5 nm; the
Meadowcroft (2008) and Maddox (2011) image-plate papers (quoted by Rosenberg only).

## Dante field of view and reference area (2026-10-29) *derived*

`DanteView(...; field_of_view, center, area)` (`src/dante.jl`; design in
`reference/diagnostics-design.md`). A field of view casts the parallel rays over a disk of
the given diameter normal to the line of sight about `center`, cut to the target's extent;
the view through a hole (no field of view) is the square about the hole, as before. Every
surface facing the instrument counts, occluded by any part; the outsides of one-sided
surfaces are dark. The reference area for Tr is the aperture's initial area (default), a
number, or a part's visible projected area (`area = :part`, Fournier 2010's convention). No
value here is fitted.

Checks (`test/spectral.jl`, "Dante field of view"):
- The through-hole path is unchanged: `dante1()`, `dante2()`, an area override and
  `sentinel()` on the coarse hohlraum with random radiosity give `isequal` intensity, Tr,
  area, cos, Tr_open and visible areas before and after the change.
- The blackbody cavity of "synthetic Dante" seen with a 60 mm field (800 pixels): the
  intensity of the view through the hole to 2e-3 and the same visible parts (the
  cavity is seen only through its hole).
- A flat Lambertian 32-gon facing the detector: I = (B/π) A to 2e-3 (1000 pixels); a
  field narrower than the disk sees π(D/2)² to 2e-3; a dark sheet over half the disk
  halves I to 1e-12 (the pixel grid and the polygon share the symmetry); Tr over the
  emitter's visible projected area is the emitter temperature to 1e-12, shielded or not.

## Several cavities: wall temperatures and the whole-target Dante (2026-11-02) *derived*

`wall_temperature(sol, names...)` and `drive_history(; cavities)` (`cavity_Tr`): each
cavity's Tr from the area-weighted mean incident flux on its wall components
(assemblies design, decision 8). No value is fitted. The built-in `nif_dante1` (and so
`dante1()`) takes the field of view of the 2004 design, ±30 mm at TCC (Dewald 2004;
section 1 of "Dante-family hardware and geometry" above, **H** for the design; the
fielded collimators are not public, and Girard 2009's "centimeter-wide" field is the one
fielded hint), with Dante 1's channels and Tr over the lower entrance hole. (C2 first
added it as `nif_dante1_field` beside a through-hole `nif_dante1`; since 2026-11-02 the
one definition is the field view, what the instrument sees, and so are the dated
configurations, `nif_dante1_2021` and older, which share the 2004 design; Dante 2 and
Omega's Dante still look through the hole, their fields being unpublished.) Cost and difference on the
default GUI session (4796 elements, N221204's shape, 49 steps, 8 threads): the history
takes 8.35–8.52 s with the field view and 8.35–8.37 s through the hole (two runs each;
the footprint is cast once per target); Dante 1's Tr differs by 7.3e-5 at its peak
(289.563 against 289.542 eV) and at most, the intensity by 2.9e-4 (400 pixels across
the target's extent against 400 across the hole). On the coarse hohlraum with a closing
lip (`test/spectral.jl`), the field's 400 pixels put Tr 4e-4 off the converged view.

Checks (`test/cavities.jl`; GUI `test_halfraum.jl`, "per-cavity wall temperature"):
- Each cavity's Tr equals Σ H A / Σ A over its wall elements (selected by component name
  in the test) to 1e-12, per step, for two hohlraums side by side driven unequally and
  for the shock-tube starter with half its lower beams off; one cavity holding every wall
  gives the history's `wall_Tr` exactly, and a history without `cavities` is unchanged.
- An 8 mm field about one hohlraum receives exactly zero from a blackbody in the other;
  it matches the view through the hole of that hohlraum to 1.2e-4 (400 pixels), a 60 mm
  field the sum of both holes' views to 1.4e-4 (1200 pixels); on the shock tube the
  whole-target view matches the view through the lower hole to 3.5e-4 (the package is black
  and the outsides of one-sided surfaces dark, so only that hole shows radiance).

## Named ports and per-channel fields (2026-11-03) *derived*

**Named ports.** `NIF_DIAGNOSTIC_PORTS` holds the NIF User Guide's (2024-09, a public
LLNL document) Table 5-1: 110 diagnostic ports, name, θ, φ and diameter, transcribed
from page images at 110 and 200 dpi (the text layer of the table's second page is
scrambled), **H**. As printed, P116-250 sits at θ 116.8°, φ 250.5° and P121-83 at θ
120.5°; every other name truncates its angles. A definition's `[line_of_sight]` may give
`port = "P143-274"` (with `facility = "NIF"`) in place of `theta_deg` and `phi_deg`,
resolved by `port_angles(:nif, port)`. Omega's diagnostic ports (H16I and the rest) are
not tabulated in the repository, so `omega_dante` keeps Rubery's angles.

Adopted: every Dante 1 definition (2009 to 2021) names P143-274, θ 143.25° (from the
papers' 143°); Dante 2 P64-350, SXI-L and G-LEH-2 P161-326, SXI-U P18-123 carry the
table's angles, which they already had. Change on the default GUI session (5228
elements, 49 steps of the default pulse shape, 8 threads, one history recording both
lines of sight): Dante 1's intensity rises by 2.00e-3 and its Tr falls by 3.18e-4 at
every step (289.564 → 289.472 eV at the peak, 7.8 ns). The reference area's projection
grows by cos 36.75° / cos 37° = 1.0033, which alone would lower Tr by 8.2e-4; the
intensity, the hole seen a little more face-on, grows by less than that projection, and
Tr ∝ (I / A cos θ)^¼ falls by the difference. Checks (`test/diagnostics.jl`,
"named ports"): 110 unique names, each within 1° of its angles; the built-ins' angles;
the errors for a port with angles, a port without a facility, an Omega port and an
unknown port; angles in an extending file replace an inherited port.

**A field of view per channel.** A channel's `field_of_view_mm` (`DanteResponse(…;
field_of_view)`) replaces the view's for that channel's signal: a disk of that diameter
about the view's center, or, in a view through a hole, the hole's square cut to that
disk. One footprint is cast per distinct field and target; channels on the view's field
take the view's footprint and signals unchanged. The intensity, Tr and the in-group
spectral shape (a Planckian at the view's `Tr_open`) stay the view's. Applied: Dante 1's
channel 2, the 2–4 keV multilayer, 12 mm (Rubery et al. 2022, OSTI 1959490, **H**; the
one public per-channel field) in `nif_dante1_2021`, so in `nif_dante1` and `dante1()`.
`dante2()` and `sentinel()` take Dante 1's channels without it (their fields are not
public). Nothing changes on a hohlraum of N221204's size (10.1 × 5.75 mm, LEH 3.2 mm; a
blackbody cavity and a driven radiosity solution, 42 groups): its extent across Dante
1's line of sight is under 6 mm, so the 12 mm and 60 mm fields cast the same rays
(`DanteView`'s rays cover the field cut to the target's extent) and channel 2 gives the
same voltage exactly. A channel's field changes a signal only on a target wider than
the field (a halfraum pair, a package off TCC). Also, `sentinel(:equatorial_78)` moves
to the table's φ 78.75° (P90-78) from 78°.

**`dante2()` from its definition.** `dante2()` builds `nif_dante2` (port P64-350, the
upper hole, 400 pixels) with Dante 1's channels without their per-channel fields: an
assumption in code (**L**), since Dante 2's channel set, solid angles and fields are not
public (section 2 above); the definition has no channels and stays unusable from the
search path until a file that extends it supplies them. Checks (`test/diagnostics.jl`,
`test/dante_response.jl`): its geometry equals the bare `DanteView(64, 350; through =
:leh_top)` it replaced, and on the blackbody hohlraum its intensity, Tr, channel
voltages and unfolded Tr equal that view's exactly (with Dante 1's channels, no fields);
against the bare view before the per-channel field, nothing changes.

Checks (`test/cavities.jl`, two hohlraums 12 mm apart as blackbodies, 400 pixels):
a channel given the 60 mm field of the view gives exactly what it gives without one
(every output, the unfold included); channel 1 with 8 mm about one cavity, in the 60 mm
view, receives exactly 0 from the other while the other channels receive it, and matches the
8 mm view's channel 1 to 1.5e-6 (the same rays; the in-group shape is each view's
Tr_open); the other channels equal those of the 60 mm view exactly; through one hole, a
60 mm channel field changes nothing and an 8 mm field about TCC (the hole 6 mm off it)
receives 0.

## Public diagnostic definitions and the time-dependent Tr area (2026-10-29) *derived*

Built-in definitions in `data/diagnostics` (docs page "Diagnostic definitions"), each
transcribed from the source table as read from the page image (not from text
extraction): `nif_dante1_2009` (Fournier 2010 Table II), `nif_dante1_2012` (May 2015
Table IV), `nif_dante1_N170924`, `_N180129`, `_N180604` (Kemp 2020 Tables I–III),
`nif_dante1_2021` (Rubery 2023 Fig. 4, phase 1), `nif_dante2` (line of sight only),
`nif_sentinel_N170924` (Kemp 2020 Sentinel Table II), `omega_dante` (Rubery 2023 Fig. 13).

Transcription decisions (each stated in its file):
- **Stacked foils.** Fournier's "Number" and May's "Qty" multiply the listed thickness:
  May's Mn 13.2 + Ni 1.8 × 2, Ti 13 × 2, Fe 15 × 2 and Cu 0.55 × 2 are Kemp's single
  Mn 26.4 + Ni 3.6, Ti 26, Fe 30 and Cu 1.1. Each foil is listed once per count.
- **Kemp N180604 ch 3** "Lexan/4.0 + Ta grid + Ta grid": two grids as printed (0.088²).
  **Ch 9** "Ge/42 C/375": read as Ge 42 µm + C 375 µm (C as graphite, 2.2 g/cm³,
  assumed); its modeled half-maximum midpoint is 1.00 × Kemp's peak energy.
- **Grids** "Yes" without a transmission (May 2015): 0.088 (Kemp 2020), assumed for
  2012–14; Fournier's 0.091 kept for 2009.
- **Mirrors.** B₄C 7°, amorphous C 5° (NIF 2009–18), C 7°, 5°, 2.5° and Be 3.5° (Omega) as
  `GrazingMirror`: Fresnel from Henke's constants at bulk densities (B₄C 2.52, C 2.2, Be
  1.848 g/cm³), no roughness, **L** (uncalibrated; the SiO₂ model's calibration to Rubery
  Fig. 3a does not carry over). The coatings' densities and roughness are not public.
- **Sentinel.** Equatorial DIM not published (φ left missing); MiniXRD responses not
  public (the NIF cross-calibration needed the Al MiniXRDs 3–4 × as responsive as an
  XRD-31, Kemp 2020 Sentinel §III B), so no cathode is given; solid angle 12.9 mm² /
  (7 m)² = 2.633e-7 sr (derived, **L**); "35% transmission sieve" as grid 0.35; no
  nominal energies published (each channel takes its response's peak).
- **Omega.** Standoff not public, so no solid angles; XRD cathodes as NIF's XRD-31 model
  (assumed); the hohlraum axis on Omega is P6–P7 or P5–P8, not the chamber z axis.

Materials added: Henke f₁, f₂ for B, F, Ca, Sc, Mn, Ge, Y, Mo (CXRO, fetched 2026-10-29);
densities B 2.34, Sc 2.985, Mn 7.3, Ge 5.323, Y 4.472, C 2.2 (graphite) g/cm³ (handbook
values, **M**); CaF₂ 3.18, B₄C 2.52 g/cm³.

Check (`test/diagnostics.jl`, "public definitions"): each channel's modeled half-maximum
midpoint against the energy its source tabulates (Kemp: peak sensitivity from the
calibrated responses; Fournier, May: band centers; Rubery: nominal energies), 108
channels in six NIF Dante 1 configurations: 97 within 5%, all within 20%; the outliers
are broad multi-edge channels (Fe + Cr + parylene near 3.5 keV, 0.82–0.92; Lexan + CaF₂ +
CH, 1.08), where a half-maximum midpoint is ill-defined. Kemp's tabulated widths are 20%
of the energy on most broad channels (a fitting default), not a response measure, and are
not compared. Nothing is fitted.

Time-dependent Tr area (`DanteView(...; area)`): `area = :open` is the hole's open area
at the solution (aperture plus each closing-lip element times its open fraction, the
entrance-hole closure of `GapClosure` / `LipClosure`), so `Tr` equals `Tr_open`; `area =
TabulatedCurve(t, A)` is a measured history at the solution's time (Glenzer 2010 took
A_LEH(t) from SXI). Checks (`test/spectral.jl`, "Dante field of view"): the open area
equals aperture + ½ lip at open fraction ½ to 1e-12; a history's Tr scales as A^(-¼) to
1e-12.

## The x-ray pinhole imager model and its public definitions (2026-10-29) *derived*

`XrayImager` (`src/imager.jl`, kind `xray_imager`; design in
`reference/diagnostics-design.md`, phase 5). No value here is fitted.

Model (each an approximation, stated in the docstring):
- **Response** R(E): the filters (Henke), an optional mirror, and the detector's
  absorbing layer, 1 − e^(−μρt), with no dead layer.
- **Spectrum** per element: the Planckian of its brightness temperature, (B/σ)^¼, for a
  gray solution (the reflected part of a radiosity really carries the incident
  spectrum), or the wall-temperature Planckian within each group of a spectral one.
  ⟨R⟩(T) is tabulated on 600 temperatures, 1 eV–5 keV, and interpolated in log ⟨R⟩ and
  log T: linear in ⟨R⟩ it was 0.4% off at 200 eV for the SXI hard channel, where
  ⟨R⟩ ∝ e^(−E/T) with E/T ≈ 22.
- **Point-spread function**: the penumbra disk, d(1 + 1/M), convolved with one Gaussian
  of the variances of the Fraunhofer FWHM, 1.03 λ p / d at the definition's energy, and
  of the detector pixel's box, (pixel/M)²/12. Projection is orthographic; the SXI
  pinhole is 1.66 m away, against a hohlraum of ≈ 10 mm.

Definitions (`data/diagnostics`):
- **`nif_sxi_lower`** (Schneider 2010 §II, read from the text): 2× magnification, the
  detector 4.97 m from the target, so the pinhole at 1656.7 mm (derived). One channel of
  the four 25 µm pinholes: Be 75 + Ti 50 µm. CCD of 15 µm Si with 24 µm pixels, 2k × 2k.
  Port P161-326. Diffraction energy 4 keV, the middle of the paper's 3–5 keV (derived).
- **`nif_sxi_upper`**: P18-123, the detector 4.88 m away (Schneider 2012, a draft with
  open "(CHECK)" marks: M); the filters are taken to be the lower one's (assumed, M).
- **`nif_gleh2`** (Chen 2021 Table I; P161-326, the 3.25 m pinhole distance at 1× and
  the filters checked against the PDF): what the paper fixes. The Ti thickness
  (5–40 µm per campaign), the Icarus sensor's absorbing layer, the frame times and the
  gate are left to a file that extends it.

Consistency, not a check of the model: Schneider 2010 calls diffraction and geometry
"about equal" for the 25 µm pinhole at 3–5 keV, with ≈ 60 µm resolution. The model's
penumbra is 37.5 µm and its diffraction FWHM 21.2 µm at 4 keV; their variances are
equal to 5%, and their sum is 59 µm. The quoted resolution's definition is not stated.

Checks (`test/imager.jl`):
- The detected radiance of a 200 eV Lambertian disk against a direct trapezoid of
  R × Planck on a separate 10⁵-point grid, to 2e-3.
- The blur conserves an image's total to 1e-12, and the rendered area matches the
  32-gon's to 1e-2.
- The point-spread function's variance per axis equals D²/16 + σ² to 1e-2 (1 µm
  sampling).
- Integrated and Gaussian-gated records of a signal linear in time match their closed
  forms to 1e-12 and 1e-6 (a square gate to 1e-3, the sampling of its edges).

## Backlit radiography: the model and its public definitions (2026-11-01) *derived*

`Radiograph` and `radiograph` (`src/radiography.jl`, kind `radiography`; design in
`reference/diagnostics-design.md`, phase 5). The source is the chain of "Backlighter
conversion efficiency and the photometric chain (2026-10-24)" (the laser energy on each
foil element, the CE fit, the angular distribution, hardware occlusion), and the
point-spread function and gates are those of "The x-ray pinhole imager model and its
public definitions (2026-10-29)"; both are shared code and unchanged. No value
here is fitted.

Model (each an approximation, stated in the docstrings):
- **Line only.** The He-α line at `kshell_line_energy` (Ge 10.279 keV; Rygg 2014 quotes
  10.25 keV), no continuum or satellites; the filters by `filter_transmission` at the
  line (NIST μ/ρ with coherent scattering where tabulated, as `backlighter_signal`), the
  detector's absorbing layer by Henke photoabsorption, 1 − e^(−μρt).
- **Object.** Each object part is a closed surface filled with one material (μ/ρ from
  `mass_attenuation`, `:auto`); the transmission along a ray is exp(−Σ μρL), L the ray's
  length inside each part from the signed sum of its crossings (a crossing on a shared
  edge or vertex counted once). Overlapping parts add, so a hollow shell is its outer
  surface filled and its inner surface at the negative density. The object is static
  over the emission.
- **Pinhole.** Each object-plane point sees the foil through the pinhole's center and
  that point; the foil's radiance there times the path's transmission, A_p/q² at the
  detector; then the imager's point-spread function (penumbra d(1 + 1/M), Gaussian of
  the diffraction at the line and of the detector pixel). The blur is applied to the
  product, so the foil's out-of-focus structure is blurred as the object's is (an
  approximation).
- **Point projection.** Each foil element is a point source at its centroid; each
  detector point sums them along their own paths through the object, so the source's
  extent sets the blur, w L/(s + L) on the object plane. The image point of a detector
  point is the crossing of its line to the foil's centroid; M = (s + L)/s.
- **Time.** The He-α emission is uniform over the laser's flat-top pulse (an assumption;
  Rygg 2014's ≈ 2 ns emission against a 2–5 ns drive is not resolved). A gated frame
  records the emission against its gate (Gaussian or square FWHM, peak 1); the weights
  come from the imager's gate integral over 400 intervals of the pulse.

Built-in definition (`data/diagnostics/nif_inflight_radiography.toml`):
- **`nif_inflight_radiography`** (Rygg et al. 2014, LLNL-JRNL-646555, OSTI 1184746, text
  read from the PDF, H): Ge foil 12 mm from the capsule in the equatorial plane (the
  target's geometry, not the definition's); eight beams, 52 kJ, CE ≈ 1%, emission
  ≈ 1 mm across and ≈ 2 ns; pinhole array imaging at 8× (20–25 µm pinholes) onto a gated
  MCP detector (Oertel 2006, Rygg's Ref. 18) of two strips with ≈ 90 ps integration
  (taken as the gate's FWHM, M); viewing holes 0.9 × 0.83 mm with HDC windows. The
  picket class for the CE (`prepulse = true`) follows Rygg's Ref. 17 (Barrios 2013: a
  1 ns picket 3 ns before a 2–5 ns drive), M.
- **Not public, so absent** (`check_diagnostics` reports them): the equatorial port
  (φ); the pinholes' diameter within 20–25 µm and their distance from the capsule; the
  filters; the MCP's detection efficiency; the backlighter pulse and its timing; the
  frame times; the object-plane field. The test completes it with a labeled mock-up.
- **Not built.** ARC radiography (Tommasini 2020; Hill 2021) uses 40–200 keV
  bremsstrahlung and Compton sources: the attenuation tables stop at 40 keV (NIST) and
  30 keV (Henke), and the source is not a He-α line, so neither the source nor the
  object's attenuation is modeled. Image-plate detectors (Rosenberg 2019's PSL response),
  streaked radiography (Hicks 2010, the streaked mode) and opacity (point-projection
  spectroscopy) radiographs are not built either.

Checks (`test/radiography.jl`, independent references):
- A 0.1 mm Al slab (a closed box) seen through a pinhole: the radiograph's transmission
  at every pixel against exp(−(μ/ρ)ρ t / cos α), α the angle of that pixel's ray, to
  1e-7 (the blur's mixing of neighboring rays); nested boxes, the outer filled and the
  inner at −ρ, against the 0.1 mm between them, to 1e-7; a box not listed as object
  obstructs every ray.
- The open field of a uniform Lambertian foil through the pinhole against the chain by
  hand, E_x/(πA) · A_p/q² · T_filters · η_Si / hν, to 1e-12.
- Point projection: every pixel against a hand sum over the foil elements, E_e/(4π r²)
  times exp(−μρL), L from slab clipping of an axis-aligned box (closed form), to 1e-10;
  M = (s + L)/s; the center's open field against a point source at s + L to 1e-4 (the
  foil's extent); partial transmission only within w L/(s + L) of the edge.
- A pinhole radiograph of an opaque edge: partial transmission only within D/2 = d(1 +
  1/M)/2 of the edge (plus half a pixel) and reaching within 1.5 pixels of it; the
  row's blocked length conserved (1 mm) to 1e-6.
- The open field at TCC (no filters or detector) equals `backlighter_signal`'s photons
  per mm² to 1e-12, with its CE, σ and extrapolation flag (a consistency check of the
  shared source).
- Gates: a Gaussian gate inside the emission records w √(π/4 ln 2)/τ of it and one at
  the pulse's start half that, to 1e-6; a square gate w/τ to 1e-3.
