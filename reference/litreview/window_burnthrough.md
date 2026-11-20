# Backed diagnostic windows: burnthrough, emergence and symmetry (2026-11-18)

ROADMAP item 28 ("Diagnostic windows"). `HeatFront` puts a 2.5 µm gold window through
about 2.7 ns before the leading edge of the window-emergent flux that Lester et al.
measured on an ignition-scale shot, and their eq. 3 estimates a drive P2 that had not been
compared with Briareus in a common convention. This review collects what the public
record states about the window, its backing, low-Z burnthrough and the eq. 3 estimate. The
calculations built on it are in the research notes, "Backed windows: what delays the
emergent flux (2026-11-18)", and in `examples/tamped_foil.jl`.

Confidence: **H** primary text read; **H-img** read from a rendered figure (by eye);
**M** abstract or secondary; **L** inference. Papers are cite-only; downloaded texts and
rendered pages stayed in the session scratchpad.

## 1. Sources read

| Key | Reference | Access | Used for |
|---|---|---|---|
| Le26 | R. S. Lester et al., "Demonstration of ignition-driven radiation transport through a THOR hohlraum", Phys. Rev. E 114, L023201 (2026), doi:10.1103/k8dz-pjdx (CC BY 4.0) | full text and figures | window build, eqs. 1–3, Dante 2 extraction, leading edge, core shape |
| Ha26 | B. M. Haines et al., "Relating metrology to performance in recent high-yield implosions on the National Ignition Facility", Phys. Plasmas 33, 112712 (2026), doi:10.1063/5.0349151 (CC BY) | Sec. V C, VII and Fig. 9 (text extracted by the coordinator) | the window shot's gold thickness and gap; applied flux P2 per shot from the neutron images |
| Ol07 | R. E. Olson et al., "X-Ray Ablation Rates of ICF Capsule Materials", APS-DPP 2007, SAND2007-7414C (OSTI 1147322); same data at IFSA07, SAND2007-5659C (OSTI 1147367) | slides | burnthrough ablation rates of Be, Cu-doped Be, CH(Ge), polyimide, HDC vs Tr |
| Ol11 | R. E. Olson, G. A. Rochau, O. L. Landen, R. J. Leeper, "X-ray ablation rates in inertial confinement fusion capsule materials", Phys. Plasmas 18, 032706 (2011), doi:10.1063/1.3566009 | **abstract only** (OSTI 21537306) | rate ranges 160–260 eV, ~T³ |
| HR03 | J. H. Hammer, M. D. Rosen, Phys. Plasmas 10, 1829 (2003) | read before (research notes, "Gold wall loss and thin-foil burnthrough (2026-10-16)") | gold power laws, the ablative front, its density dependence |
| Li04 | J. D. Lindl et al., Phys. Plasmas 11, 339 (2004), Fig. 4-6 | read before (same section) | Nova gold-foil burnthrough |
| FTI | UW Fusion Technology Institute report on BUCKY simulations of Nova gold-foil burnthrough (fti.neep.wisc.edu, fdm1017) | full text | the Nova foils: 1–3 µm gold in the wall of a 1.6 mm hohlraum, rear emission seen by a streaked x-ray imager |
| Pa24, Zy22 | Pak et al., PRE 109, 025203 (2024) Fig. 1; Zylstra et al., PRE 106, 025202 (2022) Fig. 2 | already in `n221204_pulse()`, `zylstra2022_fig2()` | a public ignition-scale drive |
| Mc25 | M. J. MacDonald et al., "Experimental assessment of planar and rippled hohlraum window designs ...", APS DPP 2025 TP13.52 | title and teaser only | an Omega window-burnthrough campaign exists (Dante); no numbers public |

Not found in open form: Olson 2011's full text (its CH points beyond the 2007 slides);
Lester et al.'s DPP 2024 (TO04.10) and DPP 2025 (JI02.6) abstracts, whose archive pages
give only a teaser; MacGowan et al., HEDP 40, 100944 (2021) and Milovich et al., PPCF 63,
025012 (2021), the sources of eq. 3's smoothing factor S₂ and recirculation factor F;
whether Nova's burnthrough foils had a backing.

## 2. The window and the measurement (Le26) **H**

- Two windows at (90°, 150°) and (90°, 330°) in chamber coordinates, opposed to cancel
  mode 1; 2 mm in diameter; 2.5–3.6 µm of pure gold on the inner face, flush with the
  wall and curved with it. The Letter's estimates use 2.5 µm; Haines et al., Phys.
  Plasmas 33, 112712 (2026), Sec. VII, describe the same shot's windows as "≈ 3 µm"
  thick and the gap around them as at most 400 µm as built, expected to close to at most
  200 µm, the shot performing as the nominal pre-shot simulation (read as the gap having
  closed). **H**
- The gold is "tamped" by 35 µm of EPON resin of 1.185 g/cm³ (4.15 mg/cm², about the
  areal mass of the 2.5 µm gold, 4.83 mg/cm²). Its stated purpose is mechanical: to stop
  the thin membrane from moving or detaching, as earlier windows did. The modeling they
  quote found the total transmission "minimally affected" by this thickness.
- Eq. 1–2: the gray optical depth τ = L/λ_Au(T_h) (τ ∼ 5 at 2.5 µm, gold assumed near
  1 g/cm³, T_h 300–350 eV) and the Eddington–Marshak transmission 4/(3τ + 2). These are
  Briareus's `eddington_transmission` and `optical_depth` defaults.
- Dante 2, near (64°, 350°), looks at the (90°, 330°) window with the upper LEH, arm and
  backlighter shielded. The window signal is the shot's Dante 2 minus a no-window shot's
  scaled to match in the foot (their Fig. 4): it is near zero to about 7 ns and rises
  steeply after. "The leading edge at 7.5 ns marks the initial x-ray emergence"; the
  walk-up shot's excess appears "only after the predicted burnthrough time (8 ns)". The
  extracted window flux is ~1.3 TW/sr near 8.5–9 ns, then the reheating burst at 10 ns
  (H-img).
- The leading edge is therefore the emergence through gold and backing together, as
  seen from outside, at the detection level of a background-subtracted difference.
- Eq. 3: ΔP₂/P₀ ≈ (2.5 n A S₂ / 4πr²)(Δα/α)(F/(F + 1)), Δα/α the window's fractional
  albedo change, F a "recirculation factor", S₂ a geometric smoothing factor, n windows
  of area A; for two 2.5 µm windows "a drive ΔP₂ of approximately −0.9%" and ΔM₂/M₀
  below 1.5%. S₂, F and r are not given.
- Core shape: reference shot P2 9% prolate, the window shot −1%; neutron M2/M0 12%
  along the window axis against 0% (their Table/Fig. 5–6). They attribute part of the M2
  to gaps around the windows.

## 3. Low-Z burnthrough rates (Ol07, Ol11)

- Method **H**: halfraum on Omega, two sample thicknesses side by side on the wall,
  burnthrough timed with a streaked x-ray imager at 0.1 of the peak signal ("by
  convention"), Dante for Tr; the ablation rate is the areal-mass difference over the
  burnthrough-time difference. Samples of ~3.5–4.5 mg/cm², the areal mass of the
  windows' backings.
- Rates **H-img** (read by eye, mg/cm²/ns at Dante Tr): polyimide (C₂₂H₁₀N₂O₅) 2.4 at
  170 eV, 4.3 at 188, 4.6 at 193; CH + 1% Ge 3.5–3.8 at 156–158 eV; pure and Cu-doped Be
  3.1–4.5 at 161–170, 5.5 at 188, 5.8–7.3 at 202–215, 7–9 at 225–230, 8.8–11.4 at
  239–243; HDC 2.3 at 175, 3.7–3.8 at 201–207, 4.8 at 220, 7–8.8 at 237–254.
- Fits ṁ = A (Tr/200 eV)³ (**L**, ours): polyimide A = 4.8 (3.9–5.4 over its three
  points), Be 5.9 (4.9–7.7 over sixteen), HDC 3.4–4.7, CH(Ge) 7.5 (one temperature).
  Ol11's abstract (**M**): Be 3–12, Ge-doped CH 2–6, HDC 2–9 mg/cm²/ns over 160–260 eV,
  "an approximate T³ dependence", Be faster than HDC and CH below 230 eV.
- No epoxy was measured. Polyimide is the closest CHNO polymer; the scatter and the
  CH(Ge) point span ×0.8–1.5 of its fit.

## 4. The gold front and density (HR03, Li04, FTI)

- In mass coordinates the diffusion flux is F = −(c/3κ) ∂(aT⁴)/∂m, so density enters the
  front's mass only through κ(ρ) and e(ρ). With HR03's fits (κ ∝ ρ^0.2, e ∝ ρ^−0.14) a
  diffusion front at fixed density heats m ∝ (κe)^(−1/2) ∝ ρ^(−0.03): holding the gold
  at 1 g/cm³ instead of letting it expand to 0.1 g/cm³ changes the heated mass by 7%
  (*L*, derived). The ablative solution has no initial-density dependence at all (its
  m₀ is a mass per area).
- The Nova foils (FTI **H**, Li04 **H-img**): gold foils of 1–3 µm in the wall of a
  1.6 mm hohlraum, rear emission timed with a streaked imager. No backing is mentioned in
  the sources read (**L**: bare, or on a backing thin enough to be ignored).

## 5. Eq. 3 and the static view-factor smoothing (**L**, derived)

- A small window of solid angle Ω (from the capsule), emitting a fraction f w less than
  the mean wall radiance, at polar angle θ_w changes the axial Legendre moment by
  a_l/a₀ = −(2l + 1) P_l(cos θ_w) (Ω/4π) f w S_l. For l = 2 at θ_w = 90°, P₂(0) = −1/2:
  a₂/a₀ = +2.5 (nΩ/4π) S₂ f w. This is eq. 3 with Ω = A/r², Δα/α = −f and
  F/(F + 1) = w, the wall's re-emitted share of the radiance. Eq. 3 as printed carries
  the sign of Δα/α; the deficit on the waist is a positive P₂ of the drive about the
  hohlraum axis (more flux at the poles), the sign that pushes a prolate core rounder,
  as Le26 found.
- S_l for a point capsule in a sphere is the clamped-cosine ratio 1, 2/3, 1/4, 0, −1/24
  (l = 0–4); for a capsule of radius x R it follows from the kernel
  K(μ) = (μ − x)(1 − xμ)/(1 + x² − 2xμ)², μ > x.
- About the hohlraum axis two waist windows are not axisymmetric: at 90° |Y₂₀|² is a
  quarter of Σ_m |Y₂m|², so the m = ±2 part carries three quarters of the l = 2 power,
  and about the window axis the same pattern is an axial P₂ of −5 (nΩ/4π) S₂ f w, twice
  eq. 3's amplitude and of the opposite sign (less drive along the window axis, the
  elongation Le26's neutron image shows along it).
