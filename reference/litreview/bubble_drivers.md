# What drives the gold bubble: irradiance, fluence and pulse epoch (2026-11-18)

ROADMAP items 10 (the P2 swing and the gold bubble) and 36 (constrained re-pointing). A
re-pointing study holds the peak laser irradiance over the wall's elements at or below a
reference pointing's, on the ground that high irradiance drives the growth of the gold
bubble that absorbs the inner beams. This review asks which irradiance measure and which
pulse epoch the published bubble data tie to its growth, and with what scaling. Applied
in the research notes, "The swing and the bubble without new data (2026-11-18)".

Confidence: **H** primary text read; **M** abstract or search-index text only; **L** our
inference. Cite-only papers: numbers with citations, no text or figures reproduced.

## 1. Sources

| Key | Reference | Access | Used for |
|---|---|---|---|
| Ca18 | D. A. Callahan et al., PoP 25, 056305 (2018), doi:10.1063/1.5020057 | full text: LLNL-CONF-742209, OSTI 1788342 (read 2026-10-22, `hohlraum_dynamics.md` §1) | the 39-shot bubble metric |
| Ra18 | J. E. Ralph et al., "The influence of hohlraum dynamics on implosion symmetry in indirect drive inertial confinement fusion experiments", PoP 25, 082701 (2018), doi:10.1063/1.5023008 | full text: LLNL-JRNL-738440, OSTI 1469454 | the speed model and its irradiance exponent |
| Iz21 | N. Izumi et al., PoP 28, 022706 (2021), doi:10.1063/5.0030826 | full text: OSTI 1773252 (read 2026-10-22) | measured tip trajectories, the E^0.2 scaling |
| Iz23 | N. Izumi et al., "Control of low-mode drive asymmetry in an efficient long pulse low-density filled hohlraum", PoP 30 (2023), doi:10.1063/5.0136145 | full text: LLNL-JRNL-842205, OSTI 1971298 | stall and relaunch; ablation scalings; the 15× contrast criterion |
| Ra21 | J. E. Ralph et al., "The effects of multispecies hohlraum walls on stimulated Brillouin scattering, hohlraum dynamics, and beam propagation", PoP 28, 072704 (2021), doi:10.1063/5.0044404 | full text: LLNL-JRNL-832507, OSTI 1860800 | Au and Ta₂O₅ bubble speeds; outer-cone intensities |
| Sc15 | M. B. Schneider et al., "The evolution of the gold bubble in NIF ignition gas-filled hohlraums", APS DPP 2015, post-deadline abstract UO7.5 | search-index text of the abstract only | 50° against 44.5° bubbles |

## 2. Which irradiance, which epoch

- **The picket's fluence, not the peak's intensity (Ra18, Iz21, Ca18; H).** Every
  published scaling of the bubble's speed or size uses the outer cones' *picket*: Ra18's
  picket energy over its FWHM net of outer backscatter, divided by the outer spot area of
  the phase plates (a fluence, kJ/cm²: the outer irradiance integrated over the
  picket); Iz21's picket fluence relative to the HDC basecamp's; Ca18's outer picket
  energy over the one-side outer spot area (53 mm²). The tip then moves at a constant
  speed through the rise and the peak, extrapolating back to the outer turn-on (Iz21:
  seven HDC SymCap shots, intercept 0.5 ± 0.1 ns for a 0.3 ns turn-on; Ra18: N160504,
  182 ± 11 µm/ns from 2.47 ± 0.24 ns). No published trajectory shows the peak's higher
  irradiance accelerating the tip in a pulse that rises gently.
- **The exponent is weak.** Speed ∝ (fluence)^0.2 ± 0.03 (Iz21, fitted to the seven
  shots; a 1.5× picket +160 µm, a 27 % lower fluence from split quads −120 µm at
  x̄ ≈ 1.8 mm); ∝ I^0.36 (Ra18, from an isothermal-expansion model with the laser and
  x-ray terms taken equal at 0.6 mg/cc; the data fit equally over I^1/3 to I^2/3, with I
  measured as the picket fluence); size ∝ (E/A)^0.5 in Ca18's 39-shot metric (R² 84 %;
  exponents 1/3 and 2/3 give 76 % and 79 %). Iz23 adds that the wall's mass ablation rate
  and exhaust speed go as the laser intensity to the 1/3 and the x-ray drive as
  T_r^1.8 and T_r^0.7 (about I_x^0.5 and I_x^0.2), and that 10 % less laser-ablated mass
  would extend adequate inner-beam propagation by ~0.4 ns.
- **The fluence on the wall, not only per beam (Iz21, H; the reading, L).** Splitting the outer quads
  lowers the wall fluence by 27 % at fixed beam energy and slows the tip as E^0.2 predicts
  (−120 µm): the overlap of spots on the wall, which a re-pointing changes, enters as
  fluence does (**L** for the re-pointing case: no published re-pointing study).
- **The peak drives the bubble when the picket does not (Iz23; H).** In a low-trough
  4-shock CH design (1 TW trough, 280 TW peak) the picket-launched bubble stalls in the
  trough and is relaunched at ~6 ns, at most 300 µm/ns, slowing again near the axis on
  the gas's back-pressure; extrapolated ballistically it would appear launched ~7 ns
  late. Iz23's criterion: the trajectory is ballistic (set by the picket) for a
  foot-to-peak intensity contrast below ~15×.
- **Cone and material (Ra21, Sc15).** Ra21 (H): on the 5.75 mm Au basecamp hohlraum
  (0.3 mg/cc, N190325 and repeats) the Au bubble moves in at 364 km/s averaged from t = 0,
  a Ta₂O₅-lined wall's 399 km/s (~10 % faster); the 50° beams start at more than twice the
  inner beams' intensity. Sc15 (M): the bubble under the 50° quads grows faster than
  under the interleaved 44.5° quads, varying with laser intensity and CBET (no numbers in
  the indexed text).

## 3. What it means for a re-pointing constraint (L)

- A ceiling on the *peak* irradiance over wall elements during the *peak* is not the
  measure the bubble data scale with. The published measure is the outer cones' picket
  fluence on the wall (energy per area during the picket, or per phase-plate spot area),
  with exponent 0.2–0.5 on the speed or size; the epoch is the picket for pulses below
  ~15× foot-to-peak contrast, and the main drive as well for low-trough pulses.
- At fixed beam powers a re-pointing changes the wall irradiance only through incidence
  angle and spot overlap; at exponent 0.2 a 20 % rise in local picket fluence is a 3.7 %
  faster tip (~0.07 mm at 1.9 mm), which on N221204 moves the swing's crossing by about
  −0.15 ns if the collision time follows E as the model has it (E^-0.16; scaled from
  the research notes' E = 1.23 run, −0.09 ns for +12 %).
- What Briareus already models: `BubbleAbsorption` takes the picket fluence ratio as a
  scalar input (`picket_ratio`, speed ∝ E^0.2, Iz21) and computes the clipping from the
  actual outer pointing (the bubble plane where the 50° axes reach the wall) and the
  actual inner beams. Re-pointing therefore already moves the bubble's effect through
  geometry when the model is in the laser–plasma chain; it does not see the wall
  irradiance.

## 4. Open

- No published bubble speed against outer-spot irradiance at fixed picket energy (e.g. a
  phase-plate or spot-size change), and no HYBRID-E bubble trajectory (searches
  2026-11-18, as 2026-10-25).
- Whether the peak's irradiance matters in a gently rising pulse beyond the I^1/3 of
  laser ablation is not separated from the picket in any published data.
