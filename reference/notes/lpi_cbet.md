# Research notes: Laser-plasma interaction and CBET

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

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
- **Briareus convention.** `Beam.dlambda` is at the (3ω) wavelength of the beam.
- **N170706 example (1ω).** 23.5° = 10530.5 Å, 30° = 10529.8 Å, 44.5°/50° = 10524.3 Å.

Published Δλ values and sensitivities:

| Design | Δλ | Sensitivity |
|---|---|---|
| NIC 2009 | 1.7 Å (3ω) for round | Inner +60%, outer −30% |
| HiFoot | 6.2 Å (1ω) | ~9 µm P2/Å |
| Low-fill CH | — | 1 Å ≈ +24 µm P2 ≈ cone fraction 33→41.5% |
| HYBRID-E (0.3 mg/cc) | 1–2 Å (1ω) | 40–50 µm P2/Å |
| N210808 / N220919 / N221204 | 1.8 / 2.5 / 2.75 Å (1ω) | Peak cone fraction 33% |

A rule of thumb is about 2.7 µm of hot-spot P2 per 1% of cone fraction (Pickworth 2020).
**H**

**Split Δλ.** One Be shot ran the 30° cone 0.7 Å shorter than 23.5° to correct
azimuthal asymmetry (arXiv:1701.09179). **M** The m=4 "cloverleaf" is attributed to
the four-fold geometry of the inner cone (Michel DPP 2024). **M**

## Backscatter

In practice, designers subtract the measured, time-dependent FABS/NBI backscatter from
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
  right of the 3ω spot along ±φ, the two columns of a quad to opposite sides. Sources: NIF
  User Guide Fig. 4-12 (top/side views) and VisRad's "NIF Unconverted Light in 3ω
  Focus Plane" sample. **H/M**
- **Rechecked against the User Guide (revision 2024-09, 2026-10-28).**
  - §4.5 and §5.6.1: 1ω ≈ 22 mm square at ≥ 4.8 mm, 2ω ≈ 12 mm square at ≥ 2.9 mm,
    an envelope of about 32 mm from TCC. Most of the unconverted light is 1ω, and
    shaped pulses leave more of it. No fraction is given.
  - Fig. 4-12: 1ω focus 460 mm beyond the 3ω focus, 2ω focus between the two. The two
    columns of a quad put their footprints on opposite sides. The figure's 4.9 mm edge and
    31 mm outer dimension are schematic (its squares are not to scale), so the text's
    4.8 mm and 32 mm are used. **H** for the text, **M** for reading the figure.
- **Deposition model (`unconverted_deposition`, 2026-10-28).** *Derived*, **M**.
  - Each beamlet's color is the uniform square beam of the TaLIS footprint geometry
    ("1ω footprint model" under TaLIS), at the User Guide's nominal size (22 mm with a
    4.8 mm edge; 12 mm with a 2.9 mm edge), converging to its focus. The 2ω focus is
    the derived +270 mm.
  - Flat top: the light at the target is the near-field aperture's image. The phase
    plate's angular spread blurs the edge by ~0.1 mm, which is not modeled.
  - Without shadowing, the power per flat triangle is exact: its projection from the
    focus onto the 3ω focal plane, clipped to the square. Shadowing is resolved by
    rays toward the lens, refined where a shadow edge crosses a piece.
  - Checks:
    - A plate in the focal plane: the totals are exact, and every element inside the
      square gets P·A/22² (or /12²).
    - A 40°-tilted plate 10 mm upstream, shadowed by a 3 mm sphere: against a 10⁶-ray
      Monte Carlo of the same beam, every lit element is within 5σ (4·10⁶ rays: max
      z = 3.5 over 447 elements).
    - A quad splits evenly between its two sides.
  - Cost: 48 quads, both colors, on the 4116-element NIF hohlraum: 0.06 s on 4
    threads. No light reaches it except 1.8 mW of 2ω on the wall's outer side.
  - **Power fraction:** an input with no default (`ratio`, per color, relative to the
    beam's 3ω power). No measured unconverted fraction vs. pulse shape was found (see
    "Still not found"). VisRad workspaces seen so far carry 5% / 5% of the total.

### Unconverted light over time (2026-11-01)

`unconverted_history`, `irradiance_history`, `threshold_crossing`, `specular_check`,
`preplasma_check` (`src/unconverted.jl`).
- **Prepulse and plasma onset.** NIF User Guide (rev. 2024-09) §4.4.3, reread from the
  text: the primary prepulse is 1ω light, kept off the 3ω focus by the wedged lens (its
  footprint 4.8 mm clear). "As a general guideline", intensities above about 10⁸ W/cm²
  begin to form weakly ionized plasma on a metal surface (the guide cites Spaeth et
  al., "Description of the NIF laser", Fusion Sci. Technol. 69, 25 (2016)). The NIF
  prepulse exceeds 10⁸ W/cm² of 1ω about 3–5 ns before the nominal start of the pulse.
  **H** for what the guide states. The value is a guideline, not a criterion. The guide
  does not give the prepulse level or shape, or where the 10⁸ W/cm² is reached (the 1ω
  footprint is the interpretation used here).
- **Derived scale.** The 22 mm 1ω footprint is 4.84 cm² in the 3ω focal plane, so it is
  at 10⁸ W/cm² when a beamlet carries 0.48 GW (4.8e-4 TW) of 1ω. *Derived*, **M** (flat
  top at the nominal size).
- **Thresholds in the checks.**
  - `preplasma_check`: 10⁸ W/cm², the guideline above, for metal foils. **M**.
  - `specular_check`: the same 10⁸ W/cm², taken as the end of cold-mirror (specular)
    reflection. *Engineering default*, **L**: no facility criterion for specularity was
    found, and a surface's reflection can change below plasma onset (melting and
    roughening, a fluence effect not modeled) or persist above it. Metals only;
    silicon and plastics have no value here.
  - For comparison, `clearance`'s 3ω tiers (10⁹ ablation and debris, 10¹² W/cm² plasma
    formation) are engineering defaults for 3ω beam clipping, a different question.
- **Fraction in time.** Still an input (`ratio`, a number or a function of time per
  color, plus an optional `prepulse` power): no measured unconverted fraction vs. pulse
  shape was found ("Still not found").
- **Cost.** The geometry and pointing are fixed and the light is linear in each beam's
  power, so a history is one deposition at unit power per (pulse, prepulse) schedule,
  times scalar powers in time; the schedules partition the beams, so the geometric work
  is one `unconverted_deposition` in all. Measured on the 4116-element NIF hohlraum with
  48 quads, 241 times, 4 threads: one deposition 0.08 s; the history with one pulse 0.09
  s, with 48 per-quad pulses 0.26 s (per-call overhead and fewer jobs per thread call);
  241 separate `unconverted_deposition` calls 17.5 s. The 3ω `irradiance_history` costs
  one `deposit` per distinct pulse. A pointing or target that changes in time needs one
  deposition per configuration (not covered).
- **Checks** (`test/unconverted_history.jl`): a plate in the 3ω focal plane under one
  beamlet, against the closed form 10¹⁴ P(t)/footprint² W/cm² inside each square for a
  shaped pulse and a time-varying 1ω fraction (1e-12); a constant fraction against
  `unconverted_deposition` of `with_pulses` at every time, for a quad (1e-12); per-beam
  pulses against the sum of single-beam histories; a linear ramp's crossing of 3e9
  W/cm² against τ₀ + Δ I* 22²/(10¹⁴ r P₀) (1e-12); a prepulse ramp's 10⁸ W/cm² onset
  against its hand value; the 3ω history against `deposit` at every time (1e-12).

### Reflected 1ω power of a TaLIS reflector (2026-11-26)

`Omega1Reflection.reflected_fraction` and `omega1_reflected_power` (`src/talis.jl`).
*Derived*, **M**, for the geometry; *engineering default*, **L**, for the reflectivity.
- **Geometry.** Each beamlet is the uniform square of the deposition model above (22 mm
  at 4.8 mm in the 3ω focal plane, converging to the 1ω focus 460 mm beyond). A
  sub-triangle of the TaLIS samples counts whole toward a beamlet when its centroid lies
  in that square and its specular reflection of that beamlet comes within the required
  20° (less twice the rule tolerance) of an open port. Its share of the beamlet power is
  its area projected along the beamlet over the area of the square at its distance along
  the beam. The fraction of the beam is the sum over beamlets divided by their number.
  Shadowing is left out, so the fraction is an upper bound for a reflector behind other
  hardware.
- **Reflectivity.** `reflectivity = 1.0` by default, an upper bound. Cold metals reflect
  1ω almost entirely (Au ≈0.98, Al ≈0.94 at normal incidence; "Unconverted light:
  absorption and conversion vs. wavelength" below), and a plasma, formed above about
  10⁸ W/cm² ("Unconverted light over time"), absorbs most of it. No facility number was
  found for the reflected power that an opposing beamline tolerates, so the estimate is
  for ranking reflectors, not a pass or fail.
- **1ω power.** An input: the `omega1` power of `visrad_unconverted` for a workspace, or
  a fraction of the 3ω power (no measured fraction against pulse shape was found).
- **What is left out.** The part of the reflected light that reaches the optics of the
  port depends on the divergence of the reflected cone and the aperture of the port at
  the 7.7 m lens distance, which the rule's 20° replaces.
- **Checks** (`test/talis.jl`): a 4 mm plate normal to one beamlet in the focal plane
  gives 16/22² of its power (1e-12); 100 mm upstream it gives 16/(22·(1 + 100/460))²
  (1e-12); a 30 × 26 mm plate over the whole square, with mesh lines on the square's
  edges, gives 1 (1e-12); the plate tilted 45° gives 0.

### Unconverted light: absorption and conversion vs. wavelength (2026-10-28)

*Engineering default*, **L**. `nif_unconverted` scales the 1ω and 2ω light relative to
3ω. The `absorption` factor scales the absorbed share of the incident light (3ω
deposition takes all of it as absorbed). The `conversion` factor scales the x-ray
conversion of the absorbed energy (it multiplies the material's 3ω conversion in a
drive).

**Measured (gold disks, Argus/Novette, ~0.5–1 ns):**

| λ (µm) | I (W/cm²) | Absorption | CE of absorbed | Source |
|---|---|---|---|---|
| 1.06 | 2.7e13 / 3e14 / 1.6e15 | 0.82 / 0.52 / 0.38 | – | Ze et al., UCRL-86414 (1981; OSTI 5947742), Fig. 4 **H** |
| 0.53 | 3e13 / 3e14 / 1e15 | ≈0.93 / 0.80 / 0.72 | – | same |
| 0.35 | 3e13 / 1.5e14 / 3e14 | ≈0.93 / 0.95 / 0.91 | – | same |
| 1.06 | 3.5e13 / 3e14 / 1e15 | – | 0.50 / 0.41 / 0.35 | Mead et al., LA-UR-86-300 (1985; OSTI 5918786), Fig. 16 **H** |
| 0.53 Argus (small spots) | 3.5e13 / 3e14 | – | 0.41 / 0.49 | same |
| 0.53 Novette (large spots) | 4e13 / 2.3e14 | – | 0.66 / 0.59 | same |
| 0.35 | 2e13 / 1e14 / 4e14 | – | 0.57 / 0.59 / 0.51 | same |

- Mead Fig. 21: on an Iλ² axis every wavelength plateaus at 0.4–0.6 below Iλ² ≈ 10¹⁴.
  The LASNEX fit CE ∝ I^-0.15 λ^-0.3 holds above that. **M**
- The Argus low-intensity conversion deficit disappears with larger spots (Mead, pp.
  9–10). NIF footprints are 22 mm, the large-spot limit. **M**
- At 2.7e13 W/cm² and 1ω, absorption falls with Z: Au 0.82, Ti ≈0.76, Be 0.68 (Ze,
  Fig. 7). **H**
- Hogan, UCRL-94321 (1986), Fig. 6, at 3e14: absorption 0.35 / 0.66 / 0.91 and CE
  0.32 / 0.51 / 0.65 for 1ω / 2ω / 3ω. These are high-intensity ratios, not used.

**The rule:** take the ratio to 3ω at the lowest intensity measured and hold it below
that.
- 1ω: absorption 0.82/0.93 = 0.88, conversion 0.50/0.57 = 0.88 (same laboratory).
- 2ω: absorption 1.0. Conversion is the geometric mean of 0.72 (Argus/Argus) and 1.16
  (Novette/Argus), 0.9.
- Ranges: 1ω absorption 0.6–1.0 and conversion 0.65–1.0; 2ω absorption 0.8–1.0 and
  conversion 0.72–1.16. In total the 1ω x-ray source is 0.77 of 3ω's, the 2ω one 0.9.

**Why low confidence:**
- Unconverted light reaches hardware at ~1e10–1e12 W/cm², 1.5–3.5 decades below every
  measurement found.
- The ratios rise toward 1 as intensity falls, and inverse-bremsstrahlung absorption is
  nearly complete in that regime (a tens-of-eV gold corona). That favors values nearer 1.
- Against that, before a plasma forms cold surfaces reflect or transmit 1ω: Au ≈0.98,
  Al ≈0.94 and Cu ≈0.97 reflectance at normal incidence, and Si is nearly transparent
  (Johnson & Christy / Palik optical constants; approximate, **M**). The NIF User Guide
  (§4.4) puts plasma formation on metals above ~1e8 W/cm², with the 1ω prepulse past
  that 3–5 ns before the pulse. The time before ignition at 1e10–1e11 was not found.
- The absolute conversion is lower than the 3ω value Briareus uses for gold (0.75,
  calibrated in hohlraums). Open 3ω disks give 0.55–0.6 of the absorbed energy at
  2e13–1e14. At tens of eV the emission is also softer than 0.1 keV, below the
  measured bands (0.1–1.5 keV).

**In `solve` (2026-11-01):** *engineering default*, **L**. With `unconverted`, each
color's x-ray source is the material's 3ω conversion at the absorbed flux of that color,
times the color's `conversion` ratio. The data above were taken one wavelength at a
time, so they say nothing about overlapping colors. Mead's Iλ² plateau would suggest
evaluating a 1ω or 2ω conversion at (λ/0.351 µm)² times its flux. Briareus does not do
that: below 10¹⁴ W/cm² the 3ω table (Dewald 2008) is held at its low-intensity end
either way, and unconverted light lies below that.

**In `drive_history` and the flux response (2026-11-02):** *engineering default*, **L**
(the rule above, applied at each step). Each step's colors are deposited by
`unconverted_deposition` from the beams' 3ω powers before the laser–plasma model
(unconverted light is taken not to take part in 3ω CBET or backscatter), at that step's
fractions. The M-band of 1ω and 2ω light is not modeled: the M-band conversions are 3ω
data (Dewald 2020), and no 1ω or 2ω value was looked for.
`flux_response`'s per-color response fixes the conversion of each color at its power (a
secant per color), exact there and for conversions independent of intensity. On
Dewald's table the secant's source exceeds the per-color one by 0.065 δ/η for a relative
flux rise δ (0.15/ln 10 = 0.065, the table's slope), and the two agree where both
fluxes lie below 10¹⁴ W/cm², as unconverted light does.
Checks (`test/unconverted_source.jl`): a two-plate history solved by hand over two steps
with a time-varying 1ω fraction (1e-12); NIF hohlraum histories against `solve` with
`unconverted_deposition` at every step, for a constant fraction and for a falling 1ω
fraction on per-cone pulses (1e-12); zero fractions reproduce the 3ω history exactly;
the secant's error on the two-plate system against the hand source change (1e-9).
Cost and effect, measured on the 4116-element NIF hohlraum (`nif_gold`, 48 quads at
8 TW on two per-cone pulses, 16 steps to 4 ns, 1ω fraction 0.08 falling to 0.02, 2ω
0.01, 4 threads): 16.4 s without and 16.6 s with the unconverted light (one
`unconverted_deposition` 0.68 s). Of the 23.2 kJ of 1ω and 2ω delivered, 0.03 kJ is
absorbed: the squares, 2.9 mm or more off the 3ω axis, miss the 1.55 mm entrance holes
and almost all of the case (the absorbed share includes the case's outside). The capsule's Tr and P2/P0 agree to the printed 1e-3 eV and
1e-4.

**Prepulse and retune (2026-11-02):** no new physics; the rule above with a power not
proportional to the 3ω pulse. `drive_history`'s `prepulse` adds, per color, a power per
beam (a function of time, one for all beams or keyed by group, quad or beam) to the
fraction times the 3ω power, as `unconverted_history`'s does; the prepulse level is an
input (the NIF User Guide gives none; "Unconverted light over time"). Beams share a
unit-power deposition when their 3ω powers are equal at every step and their prepulse
schedules are the same object. `retune_drive` with `unconverted` adds each group's
per-color response at its 1 TW-per-beam light (times the step's fraction) to the 3ω
response, and the prepulse's flux as a fixed offset of the least squares (it is not
retuned). As for the 3ω response, the conversions must not depend on intensity, so the
per-color secant does not enter. Checks: `test/unconverted_source.jl`, a 1ω prepulse
checked by itself on the two-plate system by hand at a step before and one during the 3ω pulse
(1e-12), and per-cone prepulses on per-cone pulses against `solve` with each color's
`unconverted_deposition` at r_c P₃ω + prepulse_c (1e-12; mock squares near the axis,
since NIF's miss the hohlraum's inside); a zero prepulse reproduces the history without one
exactly. `test/retune.jl` (constant albedo and conversion, mock squares that raise
the drive 5%): the reference retuned with its powers pinned gives its drive and
shape (9e-16, 2e-16); the free retune returns those powers to 1e-9 (the bounded solver's
tolerance; 1e-8 for the 3ω-only reference); a dropped bundle's retuned powers, run
through `drive_history`, give the retune's predicted drive and shape (1e-15, 2e-16);
zero fractions reproduce the retune without unconverted light exactly.

**The prepulse in the drive history (2026-11-02):** a drive history steps only at its
`times`, so a prepulse before the 3ω pulse reached the GUI's irradiance checks only. The
GUI's history now also steps through the prepulse window every 0.5 ns (ends included),
one step 0.5 ns after it, and every time between the first and last firing time,
zero-power ones included (`_history_times`): the wall cools through a gap and re-emits
what it stored. 0.5 ns is *engineering*: on the coarse hohlraum (`nz = 8`, Lindl wall,
48 quads, a flat 5 ns 1ω prepulse of 0.1 TW per beam on mock squares near the axis, main
pulse from 0 in 0.2 ns steps), the Tr at the window's end and over the first 1 ns of the
main pulse moves by ≤ 0.0053 eV (2.4e-4) from 0.125 ns steps (1 ns steps: 0.022 eV,
1e-3). Three changes
in the wall's memory, each exact where it did not apply before:
- the first step's wall-energy integral is `max(t, 0) · r(T)` (constant T since t = 0);
  a first step at negative time, a prepulse's, starts with no history (it was `t r(T)`,
  negative);
- at a step where nothing fires, a heated element without flux starts `solve`'s albedo
  iteration at half its running maximum (it took `reference_temperature`, 250 eV, which
  put the dark wall on the heating branch at 250 eV and gave negative fluxes);
- a cooled wall heated again (the step term of the release, `Fc < 0`, the surface layer
  reabsorbed) could take more than its incident flux, a negative radiosity (the capsule's Tr
  went to 0 at the first main-pulse step after a 1.25 ns gap); the radiosity is kept ≥ 0
  (albedo 0, no emission, the re-emitted flux booked as `L(T_h) − H`). Checks
  (`test/unconverted_source.jl`, "the prepulse in the drive history"): a prepulse-only
  history at -3, -2, -1 ns (constant albedo 0.8, an intensity-dependent conversion, 1ω and
  2ω prepulses on mock squares) against `solve` with each color's `unconverted_deposition`
  at that step's prepulse power: Tr equal to the last bit, H to 3e-16. On the Lindl wall, a
  1ω prepulse of 0.1 TW per beam from -3 to 0 ns (prepulse Tr 17.0 → 20.3 eV) raises the
  main pulse's capsule Tr by 9.5 eV at 0.25 ns (8.9 → 18.3), falling every step to 0.96 eV
  at 3 ns (28.9 eV); the same steps with a zero prepulse reproduce the history without one
  exactly. With the prepulse stopping at -1 ns, the zero-power steps at -0.5 and 0 ns give
  18.9 and 17.6 eV (only the re-emission of the wall, from 19.1), and no radiosity is
  negative after the gap. GUI (`test_time_unconverted.jl`): without a prepulse the history's
  times are the pulses' firing times, unchanged; with one, the window's steps are added and
  match the core's history at the same times to 1e-12.

**Not found:**
- Garban-Labaune et al., PRL 48, 1018 (1982): CH absorption at 1.06/0.53/0.26 µm down
  to 1e10 W/cm²; only the abstract was available.
- Eidmann & Schwanda, Laser Part. Beams 9, 551 (1991): Au/Cu/Al at 0.53 µm,
  1e12–1e14; only the abstract was available.
- Celliers & Eidmann, PRA 41, 3270 (1990).
- Any absorption or conversion of metals or Si at ≤ 1e12 W/cm² at 1ω or 2ω.
- A NIF measurement of unconverted-light coupling to hardware.

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
    peak. Total power is conserved post/input to < 1%.
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
    moves power from inner to outer even at Δλ = 0.
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
    1.92 and 1.88. The model favors 30°, as I23 saw in x-ray spots, but HYDRA splits
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
  - The response is linear within 10% over 0.5–2.75 Å at 5.6 ns.
  - At Δλ = 0 the net transfer is from inner to outer (cone fraction 0.336 → 0.090 at peak).
    The null is at ≈1.2 Å, inside NC24's 1–1.55 Å for HYBRID-E.
  - Caveat: HYDRA's peak post-CBET cone fraction differs by only ≈0.07 between N210808
    (1.8 Å) and N221204 (2.75 Å): 0.55–0.6 against 0.65 in Fig. 1e. The calibrated model
    is therefore ~2.5× more sensitive to Δλ than HYDRA's cross-shot difference
    (different pulses, L). Beyond the clamp, the model has no saturation of transfer
    with Δλ.
- **Benchmark rerun.** `examples/n221204.jl` with transfer `calibrated`, now the
  example's default. Standard mesh (h 0.28 / capsule_h 0.16), 20 min on 8 threads.
  - The calibrated model acts on the example's delivered powers (DN24 Fig. 1, 442 TW),
    ~19% above Fig. 9's.
  - Resulting post-transfer cone fraction: 0.21 (2 ns), 0.25 (4), 0.29 (4.5), 0.56 (5),
    0.69 (5.5), 0.66 (6), 0.63 (6.5), 0.59 (7), 0.52 (8).
  - At 4–4.5 ns this is below Fig. 9b's 0.35–0.40, because DN24's input cone fraction
    there (0.13–0.14) is below Fig. 9b's (0.17–0.18). At 5.5 ns it is 0.04 above,
    from the higher intensity.
  - The example also changes the bubble's `capsule_scale`, now 1.135/0.91 = 1.25
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
       less transfer there than Fig. 9b shows.
  2. `capsule_scale` 1.25 moves the ridge stop 0.5 ns earlier (g = 0.83 already at
     5.5 ns).
     - The 5.5 ns minimum rises from −0.138 (scale 1) to −0.086, against Fig. 10's
       −0.168. The values at 6.5–8.5 ns are unchanged, since the ridge was already fully
       on there.
     - So the swing turns up ≈0.5 ns too early and ends 0.06–0.08 high.
     - The open items are the same as before: the bubble's inner absorption (48% at
       6.5 ns and 85% at 8 ns, above Ralph's 51–62% at peak) and the wall motion
       (constant 0.10 mm/ns from the turn-on).
     - The collision time is new. Izumi's 0.4 ns per 20% of capsule radius is
       extrapolated here to 25% and to a 3.2 mm hohlraum. Whether r_cap or r_cap/R_hohl
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

## Joint multi-shot CBET calibration (2026-10-24) *derived*

ROADMAP "Shot-data calibration" (CBET bullets). `n221204_cbet()` was fitted to one HYDRA run of one shot; its
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
by −0.005 to −0.007 px/px: the right side of each panel sits 2–3 px higher); colors were
unmixed per pixel; centerlines are stroke centroids averaged over ±0.05 ns; samples were
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
  0.02–0.036 lower at 3.0–4.0 ns (outer second step ~10% higher); peak cone fraction
  equal (0.332), every cone ~2% higher. No resolvable rise shift (< 0.05 ns).
- *Cross-shot transfer in HYDRA* (post-CBET cone fraction, input 0.33 for all):

  | t (ns) | N210808, 1.8 Å (PRE 109 / PRE 106) | N220919, 2.5 Å | N221204, 2.75 Å |
  |---|---|---|---|
  | 5.6 | 0.596 / 0.624 | 0.634 | 0.645–0.648 |
  | 7.0 | 0.465 / 0.508 | 0.564 | 0.571 |

  So HYDRA's peak dCF/dΔλ is 0.03–0.06 per Å between 2.5 and 2.75 Å (same pulse and
  inputs: 0.014 at 5.6 ns, 0.007 at 7 ns per 0.25 Å) and 0.05 (5.6 ns) to 0.11 (7 ns)
  per Å between 1.8 and 2.75 Å (different pulses). The lit-review target 0.15–0.18 per Å
  was derived from hot-spot P2 at 2.7 µm/%CF; the HYDRA transfer is 2–3× flatter.
- *Not found.* No N230729 transfer curves (PoP 31 shows its pulse but no post-CBET
  powers): Δλ 2.9 Å and P2 +6 µm are in `hybride_shots()` as a trend check only. No
  published LEH Te or ne for any HYBRID-E shot (PRE 109, PoP 31, PRE 106 and HE21 give
  only the clamps and multipliers); the Te bound below is the high-fill measurement
  (B16: 4.2 ± 0.2 keV at peak) and the reduced-model value of Michel 2009 (5–6 keV).

### Joint fit (`calibrate_cbet(model, shots; penalty)`, `examples/hybride_cbet_calibration.jl`)

- **Setup.** One `TwoZoneCBET` for all shots: outside multiplier and Mach, inside
  multiplier, one Te(t) table (7 knots, 1–8 ns, both zones), CH at ne/nc 0.06, Ti 1.5 keV.
  Per shot: the HYDRA input powers of the figure for that shot (`n221204_fig9_beams` on the
  mock HYBRID-E pointing; the inner total from the input cone fraction, split evenly
  23°/30°), its Δλ, and its timing (`delay` = its `pulse_rise` minus N221204's: −0.03 ns for
  N220919, −0.06 for N210808; the model is evaluated at t − delay). Targets: post-CBET
  cone fraction every 0.2 ns over 0.8–8.4 ns where readable (107 rows, 39 / 31 / 37).
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
  50° at peak (0.39 / 0.62 against HYDRA's 0.43 / 0.64), the transfer favors 30° (2.16
  against 1.88).
- **What the data by themselves favor** (free): no flow-driven zone (outside Mach → 0, the
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
  Outside multiplier ×/÷ 1.96, inside ±18%, Mach ±0.38, Te knots ±6–14%, peak clamp
  ±27%, foot clamp unconstrained (σ_u 12: not binding, the foot is fixed by Te and the
  multipliers). Strong correlations: outside multiplier with Mach (−0.89), inside
  multiplier with the early Te knots (+0.75–0.83) and with the peak clamp (−0.83). Derived:
  null 0.97 ± 0.10 Å (the penalty holds it at the band's lower edge), dCF/dΔλ
  0.197 ± 0.016 (1–1.8 Å) and 0.173 ± 0.022 per Å (2.5–2.75 Å) at 5.6 ns. HYDRA's
  0.06 between 2.5 and 2.75 Å is > 5 σ away: the discrepancy is the model's structure,
  not its parameters.
  A profile in the null (a fit pinned at each null) is the next step; the two ends are in
  the table (null 0: loss 0.053, or 0.066 with the Te bound; null ~1 Å: loss 0.120), each
  fit 15–90 min on the shared 4 threads.
- **Not adopted as the default.** On the shots it was fitted to, the joint fit is better on
  N220919 (rms 0.030 → 0.016) and N210808 (0.068 → 0.047) and over all rows (0.044 → 0.033),
  but worse on N221204 (0.016 → 0.027), the example's default shot. `n221204_cbet()` stays
  the example's default; `hybride_cbet()` is the `joint` transfer option.

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
- The 23°/30° split (the joint fit still favors 30°), an inner-quad SBS proxy
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
  CBET." P2 is that of the hot spot's 17% self-emission contour. The same paper states
  HYBRID-E "principally uses wavelength separation … to compensate for the inner beam
  absorption in the bubble and the flow-induced CBET": at a round hot spot the inner
  cones must *gain*, so HYBRID-E's transfer null lies below 1 Å.
- **It is design-specific.** The round-hot-spot Δλ moved with the design: 1.55 Å
  (N210207/N210307, 3.64 mm LEH, 1050 µm), 1.8 (N210808, 3.1 mm LEH), 2.5–2.75 Å
  (N220919/N221204, longer pulse, thicker ablator); HE21: LEH 4.0 → 3.64 mm needed
  +0.5 Å. A transfer null would not track capsule and pulse changes.
- **The null itself (L for HYBRID-E).** The only published statement is I-Raum's:
  0.5 Å "mostly" negates the flow-induced CBET in a hohlraum built to need little
  bubble compensation, so the null is ≤ 0.5 Å there; flow-driven transfer from inner to outer
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
  30° still favored, 50° loses too much).
- **What is not physical about it.** Te dips to 2.5 keV at 5.5 ns between 4.5 keV at
  4 ns and 3.0 at 6.5 ns (B16: Te follows the drive and peaks late, 4.2 keV), and the
  peak clamp 1.4×10⁻⁴ is below HYDRA's deck values (5×10⁻⁴–7.5×10⁻³; our ion-wave
  coefficient is not HYDRA's, so the values compare only roughly). Both trade against
  the inside multiplier (correlations −0.92 to −0.97): with the hard clamp binding, the
  peak gain is set by clamp × multiplier and Te enters through Im Γ/|Γ|. The
  unconstrained fit ends in the same place, so the dip is the data's, not the null's.
- **Uncertainty** (Laplace, Gauss–Newton, σ from the residuals 0.025, penalty held
  fixed): outside multiplier ×/÷ 4.9, inside ±31%, Mach ±0.83 (−0.84 with the outside
  multiplier), Te knots ±7–15%, peak clamp ±29%, foot clamp unconstrained. Derived:
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
  default: on the N221204 Fig. 9 it is still the best (rms 0.016 against 0.024).
- Left: the Te dip and the low peak clamp (a drive-following Te(t) prior, e.g. B16's
  shape, would test whether the data need them); the 23°/30° split; HYDRA's per-run
  settings cannot be separated from Δλ in the published curves.

## Inner-cone SBS time histories and the late-time refraction (2026-10-25) *derived*

ROADMAP "The P2 swing" and "Absolute energetics on N221204". Literature: `reference/litreview/inner_cone_sbs.md` (Kemp et al., PoP
33, 082705 (2026), doi:10.1063/5.0320730, accepted manuscript OSTI 3441608; the AAC 2025
talk, OSTI 3388898; Higginson et al., PRL 134, 165101 (2025); Izumi et al., RSI 95,
103510 (2024); Kemp et al., PoP 31, 042707 (2024), outer cones). Fig. 10 and Dante are
checks; nothing is fitted to them. Confidence: **H** primary text, **H-img** digitized,
**L** ours.

### What is published

| Quantity | Value | Source | Conf. |
|---|---|---|---|
| N210808 measured inner SBS, per quad (all light: in-lens × 3.4) | 23°: 2.8 kJ, peak 2.18 TW at 5.3 ns, gone by 6.95 ns; 30°: 0.9 kJ, shoulder 0.41 TW at 5.1, peak 0.90 TW at 5.75, gone by 6.7–7.0 ns; both start as the inner power reaches peak (half of its 9 TW per quad at 4.74 ns) | K26 Fig. 2 (`kemp2026_sbs()`, `data/kemp2026_fig2.csv`; integrates to 2.87 / 0.92 kJ) | H / H-img |
| Errors | single shot up to 100% (lens fraction, calibration); 30–50% between N210808 repeats; the 23° sharp rise may be a processing artifact; DrD timing ±0.2 ns | K26 §II–III | H |
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
by 7.5 ns, the same fraction on both inner cones (16.8% of the inner power over the
plateau). The published N210808 history, on N221204's axis, is earlier, shorter and
uneven between the cones: it starts ~0.3 ns earlier (4.73 ns), peaks in the first
~0.7 ns of the peak (23° at 5.46 ns), is mostly over by 6.7 ns and gone by 7.1 (against
7.5), and the 23° cone carries 76% of it. The two agree on the
energy scale (N210808 ~30 kJ, N221204 42 ± 11) and on "early in peak power".

### `backscatter=history` (new option)

`examples/n221204.jl backscatter=history`: each inner cone loses `kemp2026_sbs()`'s power
per quad over its share of the delivered inner power, plus the same SRS as `measured`.
Time map: N221204's inner power reaches half its peak 0.16 ns after N210808's 4.74 ns on
K26's axis, so the curves move +0.16 ns (*L*: the shut-off follows the hydrodynamics,
timed here by the rise). Energy: scaled ×1.39 to AS24's 42 kJ with N210808's 3.1:1 split
(`sbs_energy=<kJ>` changes it; *L* for N221204's split). Result: 23° 31.8 kJ (peak 35.8%
of its delivered power), 30° 10.2 kJ (peak 13.9%).

N221204, full model, final defaults, coarse (0.8 / 0.4); `measured` from "Capsule
absorption, conversion and the ablating-wall release" §2, reproduced to the digit by a
fresh `measured` run in this tree:

| backscatter | min (at) | upward zero | 8 ns | rms vs Fig. 10 | P2/P0 at 5.5 / 6.5 ns | capsule / LEH share | Dante 1 8.0 / 8.5 / 9.0 (initial area) |
|---|---|---|---|---|---|---|---|
| `measured` (default) | −0.181 (5.25) | 6.19 | +0.080 | 0.047 | — | 13.8 / 24.7% | 290.9 / 282.2 / 264.6 |
| `history` | −0.152 (5.25) | 6.20 | +0.080 | 0.043 | −0.150 / +0.039 | 13.8 / 24.2% | 290.9 / 282.2 / 264.6 |
| uniform 48 kJ cut (the setting matching Fig. 10) | −0.195 (5.25) | 6.24 | +0.080 | 0.044 | −0.186 / +0.032 | 14.0 / 24.1% | 289.7 / 281.1 / 263.8 |
| PRE 109 Fig. 10 (HYDRA) | −0.170 (5.60) | 6.67 | +0.095 | — | −0.169 / −0.045 | 16.2 / 25.2% | measured 310.1 / 311.3 / 301.5 |

- **Swing timing: unchanged** (6.20 against 6.19 ns; Fig. 10 6.67). The published
  history removes the inner power earlier (mostly 4.9–6.3 ns), when the waist is most driven,
  so the minimum is 0.029 shallower (−0.152 against −0.181; Fig. 10 −0.170), but it is
  gone by ~7 ns either way, and the crossing at 6.2 ns is set by the bubble's absorption
  and the collision (ridge stop at 6.09 ns), not by the backscatter. The minimum stays at
  5.25 ns (Fig. 10 5.6).
- **Late-time drive: unchanged** (Dante 1 and capsule Tr identical to 0.1 eV from 7.5 ns
  on): in both shapes the SBS has ended before the late peak. The 6 kJ SRS is the same.
- rms against Fig. 10 0.043 (0.047 `measured`); the capsule share is the same (13.8%).

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

- Booked as lost from ~6.2 ns to the end, 50% of the 30° power (and of the 23°, 0.25 ns
  later) removes ~190 kJ: the crossing moves **0.10 ns earlier** (away from Fig. 10's
  6.67), P2 rises late (+0.085 to +0.099 at 6.5–7 ns against Fig. 10's −0.045 / +0.046),
  and Dante falls 11–18 eV through the late peak (the late deficit grows). So the
  refraction, if it removes inner light from the hohlraum, cannot be what delays HYDRA's
  swing; on its published numbers it acts the other way. If the refracted light is
  instead absorbed in the bubble and converted (the bubble model's ridge stop already
  does this from 6.09 ns), the drive loss shrinks and the symmetry effect stays.
- `refraction=0.5` at the hydrodynamic timing (spike at 7.16 ns; run 2026-10-26, 4.3 min):
  ~140 kJ booked as lost. The minimum and the crossing are unchanged (−0.181 at 5.25 ns,
  6.19 ns; uniform cut 6.24), P2 after it moves by ≤ +0.003 at 8 ns, and Dante falls
  11–17 eV over 8–9 ns (290.9 / 282.2 / 264.6 → 280.1 / 269.6 / 247.8): as expected, the
  late drive deficit grows and the swing timing does not move. Neither timing of the
  published refraction closes the 0.46 ns early crossing.

### Outer cones

Kemp 2024's outer-cone history is one 50° quad on N210808 (FABS; not cone-averaged, 44°
not studied, noise at a quarter of the peak): ~2.7 kJ over the sixteen 50° quads (*L*),
10–20% of the inner SBS, late (7–8 ns). `backscatter=history` stays inner-only: the
outer history is not published cone-averaged, and ~2% of the 50° power for ~0.7 ns moves
the post-CBET cone fraction by < 1 point (≲ 0.003 in P2/P0, after the crossing).

### Decisions (coordinator 2026-10-25; accepted by the user 2026-10-26)

- `measured` stays the default; `history` is an option; its 3.1:1 cone split is
  N210808's (Kemp 2026), unpublished for N221204.
- `refraction` stays a bound (booked as lost; overlaps `BubbleAbsorption`).

## CBET without new data: the Δλ slope and the Te history (2026-11-18) *derived*

ROADMAP item 11, the part that needs no N230729 transfer curves or per-shot plasma
states: why the calibrated models' cone-fraction response to Δλ is 2–3× steeper than
HYDRA's cross-shot response (§1–2), and whether a drive-following Te(t) removes the
oddities of the joint fit's Te history (§3). Only published numbers already in these
notes and `litreview/cbet_calibration.md`, `litreview/n221204_benchmark.md` are used
(one web search found nothing newer). No default or calibrated constant changed; every
run is a scratch script on the joint-fit setup of `examples/hybride_cbet_calibration.jl`
(the three shots' HYDRA input powers on the mock pointing, each shot's Δλ and rise
delay; N221204's beams where a single shot is named). Confidence: **H** primary text or
figure, **L** ours.

### 1. The sequence: what is published, and what the models give

| Shot | Δλ (1ω) | What else changed (**H**) | HYDRA post-CBET CF at 5.6 / 7.0 ns, input 0.33 (**H-img**) | Hot-spot P2 (µm, **H**) |
|---|---|---|---|---|
| N210808 | 1.8 Å | reference pulse; 3.1 mm LEH; 1.92 MJ (1.89 corrected), 441 TW | 0.596 / 0.465 (PRE 109 Fig. 1e; the PRE 106 run 0.624 / 0.508) | −2.3 ± 0.5 |
| N220919 | 2.5 Å | peak +~300 ps, longer trough, thicker ablator; Δλ capped by SBS risk; 2.05 MJ, 440 TW | 0.634 / 0.564 | −10.8 ± 1.6 |
| N221204 | 2.75 Å | N220919's pulse with re-tuned powers (foot CF raised, 0.09 → 0.13 at 3.4 ns) | 0.648 / 0.571 | 0.3 ± 0.1 |
| N230729 | 2.9 Å | early pulse +~30 ps, ablator +6.9 µm vs N210808 | not published | ≈ +6 |

Ruof 2026 (2DConA, N210808 design, Δλ 1.8 Å; **H**): rise cone fraction 0.21 → 0.24
raised the in-flight shell P2 from 9.1 to 27.5 µm; 0.21 → 0.16 left it at 10.9 µm (HYDRA
had predicted −10 µm).

HYDRA's slope (5.6 ns): 0.055 per Å over 1.8–2.75 Å (other pulse), 0.056 over 2.5–2.75 Å
(same pulse); at 7.0 ns 0.112 and 0.030. As a shape: (CF(1.8) − CF_in) / (CF(2.75) − CF_in)
= 0.83 at 5.6 ns, 0.54 at 7.0 ns (a response proportional to Δλ gives 0.65).

The three calibrated models on the beams and delay of each shot (N230729: N221204's
beams at 2.9 Å, tables 0.03 ns later), post-CBET CF at 5.6 ns (7.0 ns), and the hot-spot
P2 steps they imply at Pickworth's 2.7 µm per % of cone fraction (**L**: a peak-CF rule
of thumb applied to the instantaneous CF):

| Model | N210808 | N220919 | N221204 | N230729 (proxy) | ΔP2 N220919 → N221204 | ΔP2 N221204 → N230729 | ΔP2 N210808 → N221204 |
|---|---|---|---|---|---|---|---|
| `n221204_cbet()` | 0.453 (0.409) | 0.591 (0.514) | 0.641 (0.561) | 0.671 (0.585) | +13.6 µm | +8.1 µm | +50.7 µm |
| `hybride_cbet()` | 0.496 (0.441) | 0.618 (0.547) | 0.666 (0.594) | 0.692 (0.618) | +13.0 | +7.1 | +45.8 |
| `hybride_cbet_iraum()` | 0.543 (0.469) | 0.634 (0.533) | 0.674 (0.568) | 0.697 (0.586) | +10.7 | +6.1 | +35.5 |
| HYDRA | 0.596 (0.465) | 0.634 (0.564) | 0.648 (0.571) | — | +3.8 | — | +14 |
| measured hot-spot P2 step | | | | | **+11.1** | **+5.7** (confounded) | +2.6 (other pulse) |

- **On the same pulse the models' slope is the measured one.** N220919 → N221204 is
  the only step with the pulse held (powers re-tuned): the models give +10.7 to +13.6 µm
  against the measured +11.1 µm, HYDRA's runs +3.8 µm. N221204 → N230729 (foot +30 ps,
  thicker ablator, so only a trend) is +6.1 to +8.1 against +5.7 µm. The N210808 step
  crosses the pulse change, which by itself moved the hot spot oblate (N220919 at 2.5 Å
  is −10.8 µm against N210808's −2.3 µm at 1.8 Å), so it does not isolate Δλ.
- So "2–3× steeper than HYDRA" is a statement about HYDRA's runs, and where measured P2
  can arbitrate it sides with the models. HYDRA's same-pulse step is a third of the
  measured one ("Joint multi-shot CBET calibration: the null" found the same from the
  P2 sensitivity curves).

### 2. Which model term carries the slope

N221204's beams, post-CBET CF against the inner cones' Δλ (1ω), dCF/dΔλ per Å over
1.8–2.75 Å and 2.5–2.75 Å. Each variant's inside multiplier is rescaled (bisection) so
that its CF at 2.75 Å equals that of the model; the variants differ only in shape.

| Variant (5.6 ns) | `n221204_cbet()` | `hybride_cbet()` | `hybride_cbet_iraum()` |
|---|---|---|---|
| model | 0.189 / 0.186 | 0.176 / 0.173 | 0.135 / 0.137 |
| outside zone alone | 0.031 / 0.030 | 0.028 / 0.026 | 0.006 / 0.005 |
| inside zone alone | 0.114 / 0.120 | 0.118 / 0.122 | 0.125 / 0.127 |
| clamps removed (both zones) | 0.190 / 0.186 (clamp not binding) | 0.288 / 0.274 | 0.242 / 0.248 |
| overlap-limited depletion, f = 0.8 | 0.168 / 0.156 | 0.158 / 0.150 | 0.128 / 0.126 |
| overlap-limited depletion, f = 0.6 | not reachable | 0.125 / 0.100 | 0.110 / 0.095 |
| Te × 0.7 / × 1.3 (both zones) | 0.255 / 0.162 | 0.196 / 0.167 | 0.143 / 0.129 |
| HYDRA (cross-shot) | 0.055 / 0.056 | | |

At 7.0 ns (1.8–2.75 Å): models 0.143 / 0.143 / 0.089; with f = 0.6 0.112 / 0.113 /
0.083; HYDRA 0.112 (2.5–2.75 Å: 0.030).

- **The linear gain's Δλ dependence carries the slope, through zero.** In the inside
  zone (Mach 0) the pair gains are proportional to Δλ within 15%: g(1.8)/g(2.75) =
  0.56–0.66 and g(2.5)/g(2.75) = 0.86–0.91 for 23°/30° × 44°/50° pairs (θ 21–33°),
  against 0.655 and 0.909 for strict proportionality. The ion-acoustic resonance lies
  at ≳ 5 Å (1ω) for these crossing angles and Te ≥ 2.5 keV (estimate from k c_s, **L**),
  so the response is sub-resonant and linear (M13: once heated, the gain is linear in
  Δλ). The whole inside transfer (CF after the inside zone minus after the outside one)
  is 0.28–0.36 of its 2.75 Å value at 1 Å, 0.61–0.66 at 1.8 Å, 0.90–0.91 at 2.5 Å:
  proportional to Δλ. The slope is therefore the transfer that HYDRA's 0.65 at 2.75 Å
  demands, spread over 0–2.75 Å: (0.674 − 0.305) / 2.75 = 0.134 per Å for `iraum`.
- **The clamp flattens, already.** Where it binds (the two joint fits) the hard clamp
  makes the response *less* steep: without it the same CF at 2.75 Å needs a 1.8–2.75 Å
  slope of 0.24–0.29. The binding clamp turns the coupling of a pair into
  K·clamp/D·√(P_p P_q), which keeps the Δλ dependence of Im Γ/|Γ| (linear) while
  removing the superlinear part of the exchange.
- **Pump depletion is present; what is absent is local depletion, and it does not
  flatten enough.** The exchange dP/dτ = P ⊙ (K P) conserves power and drains the donor
  cones, so beam-averaged depletion is in every model. A model in which only a fraction
  f of each beam's power crosses the others is, in closed form, P_post = (1 − f) P_in +
  f · P_model (the interacting part evolves as the full model; energy is conserved
  between the cones exactly). With f = 0.6 the inside gain must rise ×1.7–1.8 to hold
  CF(2.75 Å), and the slope falls by 19–29% at 5.6 ns (`iraum` 0.135 → 0.110) and
  7–21% at 7 ns; on `n221204_cbet()` f ≤ 0.7 cannot reach the target at all (the
  exchange diverges first). HYDRA's 0.055 is not reached by any f.
- **The plasma state moves the slope by ≤ 0.02 at fixed transfer** on the joint fits
  (Te ×0.7 / ×1.3, gain rescaled: `iraum` 0.143 / 0.129, `hybride_cbet()` 0.196 / 0.167); a
  colder plasma brings the resonance nearer and steepens it (`n221204_cbet()` 0.255).
- **The outside zone adds 0.005–0.03 per Å** (its flow-driven transfer from inner to outer
  shrinks as the inner cones redden).

What HYDRA's cross-shot numbers would require instead: the inside-zone gain that puts
each shot on HYDRA's curve, at the beams and delay of that shot (inside multiplier factor):

| Model | 5.6 ns: N221204 / N220919 / N210808 | 7.0 ns |
|---|---|---|
| `n221204_cbet()` | 1.02 / 1.12 / 1.57 | 1.03 / 1.18 / 1.30 |
| `hybride_cbet()` | 0.95 / 1.05 / 1.38 | 0.93 / 1.06 / 1.11 |
| `hybride_cbet_iraum()` | 0.92 / 1.00 / 1.24 | 1.01 / 1.14 / 0.97 |

- A plasma state that differs between the shots would have to give N210808 35–54% more
  gain than N221204 at 5.6 ns but −4 to +26% at 7 ns, and N220919 13–15% more than
  N221204 at 7 ns on the same pulse. Nothing published differs between those shots at
  peak that would do this (same LEH, peak power, fill; the N210808 pulse is shorter but
  equal at 5.6 ns); the published per-run differences are HYDRA settings (PRE 109: the
  saturation retuned after N221106; PRE 106: N210808's peak flux P2 × 0.7; the two
  N210808 runs differ by 0.03–0.09 in CF).

**Diagnosis.** No term of the model carries an *excess*: the slope is that of a
sub-resonant linear gain through a null near 0 Å, with the clamp already flattening it.
To bend a response that passes the I-Raum null into HYDRA's shape (0.83 of the 2.75 Å
transfer already at 1.8 Å) a saturation in Δλ is needed that neither the clamp, nor
local depletion, nor the plasma state supplies; the null itself would have to sit near
−3.5 Å ("the null" section). The one measured same-pulse step agrees with the models
(+10.7 to +13.6 µm against +11.1 µm), not with HYDRA (+3.8 µm). HYDRA's flat cross-shot
response is best read as per-run settings, not as a Δλ response to reproduce (**L**).
No depletion keyword was added: the measurement shows it cannot close the gap, and the
measured P2 does not ask for it.

**Ruof 2026's rise response.** The models amplify a change of the rise input CF
(4.0–4.2 ns, N221204's beams, total fixed) almost symmetrically: post-CBET CF 0.32–0.34
/ 0.40–0.42 / 0.43–0.46 for inputs 0.16 / 0.21 / 0.24, dCF_post/dCF_in 1.25–1.55 down
and 1.18–1.41 up. A cut to 0.16 is amplified at least as much as a raise to 0.24, so the
reduced CBET cannot produce Ruof's asymmetry (+18.4 µm shell P2 up, +1.8 µm down, where
HYDRA predicted −10 µm); its cause lies outside the transfer at the LEH (**L**).

### 3. A drive-following Te(t) prior (a measurement; not adopted)

`hybride_cbet_iraum()`'s Te(t) dips to 2.5 keV at 5.5 ns between 4.5 keV at 4 ns and 3.0
at 6.5 ns, its peak clamp (1.4×10⁻⁴) lies below HYDRA's deck values, and Te is at the
5.5 keV bound at 8 ns. B16 (**H**, high fill): LEH Te "follows the drive" and peaks late.

- **Prior.** Te at the 4, 4.75, 5.5, 6.5 and 8 ns knots tied to N221204's measured
  Dante 1 radiation temperature (Pak et al., PRE 109, 025203 Fig. 1, `n221204_pulse()`:
  178.8, 219.9, 263.1, 282.8, 310.0 eV): Te(t) = Te(8 ns) · (Tr(t)/Tr(8 ns))^β, with
  Te(8 ns) ≤ 5.5 keV and β ∈ [0, 4] fitted; the 1 and 3 ns knots stay free (Dante 1
  starts at 3.3 ns and B16 does not discuss the foot). The form is ours (**L**); B16 gives
  only the shape. 9 parameters instead of 12.
- **Data and everything else as `iraum`.** The same 107 rows of HYDRA post-CBET cone
  fraction (N221204 PRE 109 Fig. 9b, N220919 Fig. 9b black, N210808 PRE 109 Fig. 1e with
  PRE 106 Fig. 1d's incident powers), the same null penalty (0.25–0.5 Å at 5.6 and 7 ns),
  bounds, hard-minimum clamps on the 3.88–4.08 ns schedule, start at `iraum`'s values
  (Te(8 ns) 5.4 keV, β 1). PRE 109 Fig. 10 is not used. The same loss gives 0.05944 for
  `hybride_cbet_iraum()` (the notes' 0.059). Bounded L-BFGS-B with forward-mode gradients
  on 8 threads, two legs of 196 + 237 s (77 iterations, not converged; the second leg
  moved the loss by 0.0004).

| | `hybride_cbet_iraum()` | Dante prior |
|---|---|---|
| loss (107 rows + penalty) | 0.059 | 0.122 |
| CF rms all (N221204 / N220919 / N210808) | 0.024 (0.024 / 0.016 / 0.028) | 0.034 (0.029 / 0.026 / 0.043) |
| peak-mean residual 5.4–7.4 ns, N210808 | −0.030 | −0.049 (worst −0.100 at 5.8 ns) |
| Te (keV) at 1, 3, 4, 4.75, 5.5, 6.5, 8 ns | 3.2, 3.9, 4.5, 3.8, 2.5, 3.0, 5.5 | 2.1, 2.5, 3.2, 3.9, 4.7, 5.0, 5.5 (β 0.97) |
| clamp foot / peak | 9.1e-4 / 1.4e-4 | 4.3e-3 (not binding) / 6.5e-4 |
| outside Mach; multipliers out / in | 1.49; 0.072 / 1.60 | 0.46; 0.091 / 0.59 (outside pair drifting along its valley: 0.86; 0.048 after the first leg) |
| null at 5.6 / 7.0 ns | 0.24 / 0.29 Å | 0.23 / 0.24 Å |
| dCF/dΔλ at 5.6 ns, 1–1.8 / 2.5–2.75 Å | 0.134 / 0.137 | 0.117 / 0.116 |
| CF at 5.6 ns, each shot's beams (HYDRA 0.648 / 0.634 / 0.596) | 0.674 / 0.634 / 0.543 | 0.624 / 0.592 / 0.504 |
| N221204 5.6 ns post/input 23/30/44/50 (HYDRA 1.92/1.88/0.43/0.64) | 1.92 / 2.09 / 0.44 / 0.54 | 1.75 / 1.96 / 0.52 / 0.61 |

Blind check, PRE 109 Fig. 10 (not fitted): `examples/n221204.jl 0.8 0.4 full iraum x
<shot>` on a scratch copy that substitutes the prior model (same `leh_z` and delays;
coarse mesh, current defaults, measured backscatter; the uniform 48 kJ cut in brackets).
Both columns run today on the same code; the `iraum` rms differs from the 2026-10-25
section's (0.040 / 0.044 / 0.052) through the defaults changed since.

| Shot | Fig. 10: min (ns), zero | `iraum`: min (ns), zero, 8 ns, rms | Dante prior: min (ns), zero, 8 ns, rms | post-CBET CF at 5.5 / 7 ns, `iraum` → prior |
|---|---|---|---|---|
| N221204 | −0.170 (5.60), 6.67 | −0.182 (5.50), 6.27, +0.077, 0.038 [0.037] | −0.197 (5.00), 6.22, +0.076, 0.046 [0.043] | 0.631 / 0.535 → 0.573 / 0.555 |
| N220919 | −0.145 (5.60), 6.66 | −0.164 (5.50), 6.23, +0.077, 0.037 [0.037] | −0.175 (5.00), 6.19, +0.076, 0.043 [0.040] | 0.595 / 0.512 → 0.542 / 0.529 |
| N210808 | −0.141 (5.50), 6.68 | −0.147 (5.25), 6.14, +0.063, 0.056 [0.054] | −0.147 (5.00), 6.09, +0.063, 0.068 [0.065] | 0.491 / 0.439 → 0.456 / 0.453 |

- **The oddities go, at a price.** Tied to the drive, Te rises monotonically
  (3.2 → 5.5 keV over 4–8 ns, Te ∝ Tr^0.97) and the peak clamp moves up to 6.5×10⁻⁴,
  inside the 5×10⁻⁴–7.5×10⁻³ of HYDRA decks (K18, PRE 109; Ruof 2026's fits 0.005–0.007;
  the ion-wave coefficients differ, so the comparison is rough). So the low clamp was
  the partner of the dip, as the correlations of −0.92 to −0.97 suggested. Te at 8 ns
  stays at the 5.5 keV bound: the late transfer still needs a weaker coupling (hotter
  electrons) than the bound allows.
- The fit loss doubles (0.059 → 0.122), back to `hybride_cbet()`'s 0.120. Every shot
  fits worse, N210808 most (rms 0.028 → 0.043; its peak CF 0.504 against HYDRA's 0.596).
  Under this model the HYDRA cone-fraction histories require the dip; which feature of
  the three histories it fits was not isolated.
- The Δλ slope falls 13–15% (0.117 per Å), still twice HYDRA's cross-shot 0.055.
- **Blind check: worse on every shot.** rms against Fig. 10 +0.006 to +0.012 (N221204
  0.038 → 0.046), the minimum 0.25–0.5 ns earlier (5.00 ns against Fig. 10's 5.5–5.6),
  the upward crossing 0.04–0.05 ns earlier; 8 ns unchanged. The prior transfers less at
  5.5 ns (CF 0.573 on the delivered powers, against 0.631) and the swing's minimum moves
  into the rise.
- **Not adopted.** The calibrated constants are unchanged. The prior removes what looks
  unphysical in the Te history but pays for it in the fit to the data it was fitted to and
  in the blind check, so the dip is not an artifact the data are indifferent to. Since §2
  interprets HYDRA's cross-shot differences as per-run settings, the dip is best interpreted
  as the model absorbing those settings into the one place that varies in time (**L**).

Not run (lock time): the laser-power form of the prior (Te ∝ P_laser^α, flat through the
peak), a Laplace covariance of the prior fit, a third optimizer leg, the standard mesh.

### What remains

- HYDRA's flat cross-shot Δλ response is not a model target: on the one same-pulse
  step the measured hot-spot P2 agrees with the calibrated models' slope, and no model
  term (clamp, local depletion, plasma state, the outside zone) bends a response
  through the I-Raum null into HYDRA's shape. A calibration against measured P2 steps
  (N220919 → N221204, +11.1 µm; the N230729 trend) instead of HYDRA's per-shot cone
  fractions would need a hydro transfer from flux P2 to hot-spot P2.
- N230729 transfer curves and per-shot plasma states (unchanged); Ruof's asymmetric
  rise-CF response lies outside the reduced CBET.
