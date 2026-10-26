# Polar direct drive: literature, benchmarks and a staged plan (2026-10-26)

ROADMAP item 15: "Polar direct drive mode, benchmarked against published GA optimizers."
This review asks what a PDD mode needs that Briareus lacks, which published optimizers
and objectives exist, which published numbers a Briareus PDD mode could be checked
against (and whether their inputs are complete), and in what order to build it. Research
and design only; no code was changed.

Confidence: **H** primary text read; **M** abstract, slides or secondary citation only;
**L** our inference. Cite-only papers: numbers with citations, no text or figures
reproduced.

## 1. Sources

| Key | Reference | Access | Used for |
|---|---|---|---|
| Ba24 | D. E. M. Barlow, A. Colaïtis, D. Viala, M. J. Rosenberg, I. Igumenshchev, V. Goncharov, L. Ceurvorst, P. B. Radha, W. Theobald, R. S. Craxton, M. J. V. Streeter, R. H. H. Scott, K. Glize, T. Chapman, J. Mathiaud, "Optimization methodology of polar direct-drive illumination for the National Ignition Facility", PRL 133, 175101 (2024), doi:10.1103/PhysRevLett.133.175101 | full text: accepted manuscript OSTI 2475934; preprint arXiv:2311.18342 (v2) | the only published GA PDD optimizer; Table I (optimized pointing), Table II (σ₁, σ₂, P̄abl, f), appendices (defocus model, pressure proxy) |
| Pol | D. Barlow, "Polaris" (github.com/DuncanBarlow/Polaris, MIT licence; successor of the "PDDOptimisation" repository cited in Ba24) | code read (README, beam tables, deck generator, Ifriit base input) | how Ba24's σ₁ is evaluated; NIF/LMJ/OMEGA beam tables; a 1-D-averaged Ifriit plasma profile (`plasma_profiles/ifriit_1davg_input.nc`) |
| Ye21 | C. B. Yeamans, G. E. Kemp, Z. B. Walters, H. D. Whitley, P. W. McKenty, E. M. Garcia, Y. Yang, R. S. Craxton, B. E. Blue, "High yield polar direct drive fusion neutron sources at the National Ignition Facility", Nucl. Fusion 61, 046031 (2021) | full text: LLNL-JRNL-814327, OSTI 1772693 | complete 24-parameter pointing tables (Table 1); SAGE optimization; Eq. 3 azimuthal model; V_rms 2.22 % |
| Co08 | A. M. Cok, R. S. Craxton, P. W. McKenty, "Polar-drive designs for optimizing neutron yields on the National Ignition Facility", PoP 15, 082705 (2008), doi:10.1063/1.2975213 | full text (author copy, pas.rochester.edu) | Table I pointing/defocus/split quads; CPP 50 % ellipses; geometric-optics defocus model; Table III energy fits |
| Ho15 | M. Hohenberger et al., "Polar-direct-drive experiments on the National Ignition Facility", PoP 22, 056308 (2015) | full text: OSTI 1182632 | the standard NIF PDD pointing (Table I); absorption with/without CBET; measured P₂–P₆ |
| CM18 | T. J. B. Collins, J. A. Marozas, "Mitigation of cross-beam energy transfer in ignition-scale polar-direct-drive target designs for the National Ignition Facility", PoP 25 (2018), doi:10.1063/1.5039513 | full text: OSTI 1467080 | ignition-scale pointing (Table II), parametric spots with spot-masking apodization, detuning (Table I) |
| MN17 | M. Murakami, D. Nishi, "Optimization of laser illumination configuration for directly driven inertial confinement fusion", Matter Radiat. Extremes 2, 55 (2017), doi:10.1016/j.mre.2016.12.002 (open access; arXiv:2512.04722) | full text | analytic PDD model (Eqs. 28–31, 37) and a fully specified 2-cone-per-hemisphere solution |
| Te14 | M. Temporal, B. Canaud, W. J. Garbett, R. Ramis, "Comparison between illumination model and hydrodynamic simulation for a direct drive laser irradiated target", Laser Part. Beams (published version not checked), arXiv:1406.7762 (2014) | full text (arXiv) | single-ring PDD σ vs shift for circular and elliptical super-Gaussians |
| TM11 | Y. Tsuda, M. Murakami, "Optimization of illumination configuration for polar direct drive", PLASMA2011 poster 23P167-B | full text (2 pp.) | potential-gradient PDD optimizer on ring geometry (qualitative) |
| SL83 | S. Skupsky, K. Lee, J. Appl. Phys. 54, 3662 (1983) | cited in MN17, Te14, TM11 (**M**) | separability σ_l = |a_l|·G_l for beams aimed at the centre |
| Sk04 | S. Skupsky et al., "Polar direct drive on the National Ignition Facility", PoP 11, 2763 (2004), doi:10.1063/1.1689665 | abstract (OSTI 823148) | the concept: repointing, equatorial power, elliptical spots |
| Cr05 | R. S. Craxton et al., "Polar direct drive: proof-of-principle experiments on OMEGA and prospects for ignition on the NIF", PoP 12, 056304 (2005), doi:10.1063/1.1876252 | abstract (OSTI 20736596) | OMEGA 40-beam PDD emulation; Saturn ring |
| Ma06 | J. A. Marozas et al., "Polar-direct-drive simulations and experiments", PoP 13, 056311 (2006) | abstract (OSTI 882305) | spot shapes, pulse shapes, Saturn designs |
| Co12 | T. J. B. Collins et al., "A polar-drive-ignition design for the NIF", PoP 19, 056308 (2012) | abstract (OSTI 22072412) | gain 32 in 2-D |
| Ma18 | J. A. Marozas et al., PRL 120, 085001 (2018) and the companion PoP 25 (2018), doi:10.1063/1.5022181 (OSTI 1423123, 1444067) | cited (**M**) | wavelength-detuning CBET mitigation on NIF PDD |

Not found: any other peer-reviewed GA PDD optimizer (searched NIF, LMJ, OMEGA, SG-III);
any published gradient or adjoint PDD optimizer (the research notes' "Optimization"
bullet already records this for illumination in general). Polaris's README lists Bayesian
and gradient-descent drivers besides the GA, unpublished (**M**). The LLE high-school
reports (Craxton's programme) contain many PDD pointing studies; none was found that
uses a GA (**M**).

## 2. What PDD needs that Briareus lacks

Briareus today deposits straight-line, phase-plate-shaped, f-cone-blurred beams on
surface patches (`DepositionOperator`, `deposit`), measures moments on a capsule
(`legendre_moments`, `ylm_moments`, exact facet solid angles), and optimizes per-group
aims (`GroupPointing`: `:z, :r, :phi, :defocus, :power, :spot, :dlambda`) with ForwardDiff.
For a bare sphere this already gives the *hard-sphere* (vacuum, no plasma) illumination
that Ba24 calls σ₁ (**L**). What is missing, in order of need:

1. **Target-surface pointing.** PDD papers give aims as points on the initial sphere:
   (θ, φ) on the surface (Ho15 Table I; Ye21 Table 1, with Δφ split left/right within a
   quad; Ba24 Table I as polar angle plus an azimuthal offset from the port), as a shift
   perpendicular to the beam (Co08 Fig. 1b), or as a far-field translation (CM18
   Table II). Briareus's `aim` is a translation of the axis through a point, and
   `GroupPointing`'s aim fields are cylindrical hohlraum coordinates. Needed: a mapping
   surface (θ_s, φ_s; R) → `aim` (the axis passes through that surface point, first
   intersection along k) and a parameterization in those coordinates, including per-beam
   split-quad offsets (Co08: split-quad x/y shifts; Ye21: four aim points per quad) (**L**).
2. **Capsule-frame metrics at high l.** Ba24 uses σ over l ≤ 30 (HEALPix nside 64 in
   Ifriit; Pol); Ye21 splits σ into the m = 0 part and the m ≠ 0 part. Briareus's
   objective cache stops at l = 6 (`_CACHED_YLM_LMAX`), and the capsule meshes are sized
   for hohlraum low modes. l ≈ 30 needs ~60+ elements around a great circle (**L**) and a
   check of the exact-solid-angle weights' cost at that l.
3. **Absorption versus incidence angle.** Without plasma, every model in §3 uses either
   the incident projected flux I·cos ψ (Te14, Ba24 σ₁; **L** for Ba24) or a closed-form
   absorption fraction: MN17 Eq. 29, η_a(ψ) = 1 − (1 − η)^{cos³ψ} with
   η = 1 − exp(−8ν_ei L / 3c) (exponential corona of scale length L). The operator knows
   each sample's direction and normal, so a per-sample weight η_a(ψ) is cheap and
   differentiable (**L**).
4. **Refraction and inverse-bremsstrahlung absorption in the corona.** All the
   hydro-based designs (Co08, Ye21, Ho15, CM18, Ba24's σ₂) deposit along refracted rays:
   oblique equatorial rays turn at lower density (n_e = n_c cos²ψ for a planar
   gradient), deposit less and spread their energy over a wider angle (Co08 §III
   explains SiO₂ vs CH this way; Ye21 reports 94.5 % absorption, all inverse
   bremsstrahlung). Briareus has no corona. In a spherically symmetric profile each ray
   stays in the plane through the centre, so its deposition depends only on its impact
   parameter b and rotates with its plane: a 1-D family of kernels K(b; r, Δα) can be
   tabulated once per profile and summed over the existing beam samples, keeping the
   operator structure and differentiability (**L**). Ba24's approximation 2 (an
   angularly averaged plasma serves for 3-D beam evaluation) is the same assumption (**H**).
5. **A drive proxy.** Ba24 Eq. 2 turns absorbed intensity into ablation pressure,
   P_abl = 24.7 Mbar·(F/10¹⁴ W cm⁻²)^{2/3}, with F a radial sum of r²·n_r^{2/3}·I_r
   normalized at the 351-nm critical radius and density (the n^{2/3} weight is empirical,
   fitted to their 3-D runs; the 2/3 power follows Manheimer et al. 1982) (**H**). Ye21
   Eq. 3 scales shell displacement with deposited energy per steradian to the 0.72 power
   (from implosion simulations) (**H**). Both are engineering proxies, not facility
   criteria.
6. **Time dependence.** PDD drive uniformity swings as the corona grows and CBET
   increases (Ba24 discussion; Ho15: P₂ prolate without CBET, oblate with it, in 2-D
   DRACO) (**H**). Needed: corona snapshots (Ba24 uses one, at 4 ns, plus the no-plasma
   state), per-cone pulse shapes (Ho15 Fig. 4; `with_pulses` exists) and a multi-snapshot
   objective. Briareus's `Keyframes` (dynamics.jl) is the hohlraum analogue.
7. **Phase plates for PDD.** Indirect-drive CPPs used defocused (Co08, Ye21, Ho15, Ba24;
   Briareus has the LMF-measured NIF CPP fits, `nif_spot`); the planned PDD plates with a
   "secondary ellipse" (Ho15 §IV) and skewed ellipses for the equatorial beams (Ho15 §I);
   circular super-Gaussians enveloped by a higher-order masking super-Gaussian
   (spot-masking apodization: CM18 uses orders 3 polar / 2.2 others, 5 % contour at the
   initial radius R₀, mask order 10 with 5 % contour at 1.15 R₀) (**H**).
   `SuperGaussianSpot` has an elliptical core plus an additive halo; a multiplicative mask
   and a two-ellipse sum are new (**L**).
8. **Defocus model.** Co08 computes out-of-focus spots by geometric optics from the
   best-focus profile (rays from a lens-plane grid with angles drawn from the best-focus
   distribution), checked against an OMEGA SG4 beam at 0, 5 and 10 mm (Fig. 4) (**H**).
   This is Briareus's f-cone convolution in substance (**L**). Two gaps: NIF's near field
   is square, so large defocus tends to a square spot (Co08 Fig. 3c), while
   `cone_directions` uses a disk quadrature; and Ba24 instead scales the super-Gaussian
   1/e radii linearly, w(d) = w₀(1 + κd), κ = {0.30, 0.40, 0.80, 0.95}/50 mm⁻¹ for the
   23.5°/30°/44.5°/50° cones, fitted to pF3D, with Ifriit treating each beam as a plane
   wave (**H**).
9. **Direct-drive CBET.** The dominant loss and symmetry driver in NIF PDD (Ho15: 89 %
   absorption without CBET vs 67 % with it in a DRACO example, Fig. 5; CM18: losses
   concentrate over the equator where repointed beams from both hemispheres overlap,
   mitigated by hemispheric or banded detuning, Table I) (**H**). `cbet.jl` is a
   two-zone LEH model calibrated to N221204 and does not transfer (**L**).

## 3. Published optimizers and their objectives

| Work | Facility / target | Free parameters | Objective | Evaluator | Search |
|---|---|---|---|---|---|
| Ba24 | NIF, solid CD/CH sphere R = 1.1 mm (N190204-003) | 16: per cone pair (θ_s, φ_s) surface aim in a disk of radius R about the port's nearest surface point (area-preserving square-to-disk map), defocus 0–10 mm, power 2–4 TW/quad | maximize f = 10·exp(−√(σ₁²/9 + σ₂²/18))·(P̄abl/50 Mbar)², σ₁ = vacuum surface-intensity rms, σ₂ = ablation-pressure rms in a 4 ns angle-averaged plasma with CBET, both in % of mean over l ≤ 30 (form reconstructed from the garbled extraction and checked: it reproduces Table II's f = 0.18 and 0.35 to two digits) | Ifriit inverse ray trace on a fixed 1-D plasma, ~1000× cheaper than rad-hydro per evaluation | PyGAD GA (185 evaluations) then coordinate descent (365); < 2 M CPU h; > 10,000 evaluations needed to converge σ₁ alone **H** |
| Ye21 | NIF, CH shells 3 and 4 mm diameter, 1.5 mg/cc DT | 24: per sub-ring (8 per hemisphere: upper/lower beams of each quad ring) θ, Δφ (left/right split), defocus | minimize velocity nonuniformity V: m = 0 from SAGE 2-D centre-of-mass radius, m ≠ 0 from the 3-D time-integrated deposited energy via Eq. 3 (exponent 0.72), combined in quadrature | SAGE 2-D hydro + 3-D ray trace | manual expert iteration; integer degrees only (1° ≈ 35 µm at 2 mm) **H** |
| Co08 | NIF, 1.1 mm SiO₂/CH thin shell, 350 kJ–1.5 MJ | 16: per ring pointing shift, defocus, split-quad x and y | minimize centre-of-mass radius rms (and V_rms) at 2.8 ns | SAGE 2-D | manual; quadratic fits in energy (Table III) **H** |
| Ma06, Co12, CM18 | NIF ignition designs | pointing to 3 rings/hemisphere, spot shapes, per-cone pulse shapes, detuning | shell uniformity, gain | DRACO 2-D + 3-D ray trace (Mazinisin), CBET (Adaawam) | manual; Ma06 automated pulse-shape tuning (**M**) |
| MN17 | idealized N_C cones (continuous rings) | per cone: axis polar angle j, axis–polar-axis intersection z (stand-off), energy E | zero the Legendre amplitudes B₂…B₁₀ of the absorbed intensity | analytic (Eqs. 28–31: parallel rays, SG profile, η_a(ψ)) | simultaneous nonlinear equations **H** |
| Te14 | one ring at 49°/131° (LMJ second cone), R₀ = 0.5 mm | shift Δ toward the equator; SG widths a, b and order m | σ_2D of the illumination (Eq. 1) | Skupsky–Lee illumination model, compared with 2-D hydro shock-front rms | parameter scans **H** |
| TM11, MN17 §5 | ring systems / symmetric DD | ring angles, axis offsets | rms absorption nonuniformity | separable model (SL83) | potential gradient / Coulomb self-organization **H** |

All published PDD optimizers are derivative-free or manual. Only Ba24 uses a GA, and
Ba24 states that the search procedure was not its focus (**H**).

## 4. Benchmark candidates, ranked by reproducibility from published inputs

**1. MN17 §4, analytic PDD solution (fully specified; H).** Two cones per hemisphere,
each an axisymmetric continuous ring of parallel-ray beams whose axes make angle j with
the polar axis and cross it at z (units of R₀); SG profile exp[−(x/aR₀)^b]; absorption
η_a(ψ). Published solution: j₁ = 46°, j₂ = 53°, z₁ = −0.65, z₂ = 0.80, E₁ = 1,
E₂ = 0.40 for a = 0.6, b = 6, η = 0.9, giving B₂ = … = B₁₀ = 0. The 4-cone centred
solution (θ₁ = 70°, θ₂ = 30°, E₂ = 0.53, independent of a, b, η) is also given; checking
it with the separable formula G_l = ΣE_i P_l(cos θ_i) leaves |G₂| ≈ 0.007 and
|G₄| ≈ 0.009 (**L**), so the published values are rounded and a test must re-solve
them, not assert B_l = 0 at the printed digits. Reference: direct numerical integration
of MN17 Eqs. 28–31 (independent of Briareus's operator). Briareus emulation: N beams per
ring (N → many), f-number → ∞, defocus 0, η_a weight. *Gap*: none beyond §2 items 1–3.

**2. Te14 single-ring σ_2D (nearly complete; H).** Ring at 49°/131°, R₀ = 0.5 mm, SG
I₀·exp{−[(x/a)² + (y/b)²]^{m/2}} (same form as `SuperGaussianSpot`), m = 4, y in the
meridional plane. Centred beams: circular a = b = R₀ → σ_2D = 11.9 %; elliptical
a = 2b, b = R₀ → 5.2 %. Optimum shift: Δ/R₀ = 13.7 % → 1.1 % (circular), 7.4 % → 1.3 %
(elliptical). Hydro absorption 85 % vs 61 % (circular vs elliptical, 2-D) is model
output, not a benchmark. *Ambiguities*: whether the illumination is I·cos ψ with full
absorption (likely, **L**), whether Δ is perpendicular to the beam axis, and the ring is
continuous (Eq. 1 is axisymmetric). Four numbers, two conventions to try; a mismatch
under both would be informative.

**3. Ba24 σ₁, NIF 192 beams, vacuum (most relevant; partial; H).** OC parameters
(Table I, per cone 1–4: polar angle 6.36/14.95/52.48/85.62°, azimuthal offset
52.02/160.75/−12.14/6.90°, defocus 8.22/9.07/4.67/6.83 mm, power 66.26/66.26/95.05/82.88 %
of 4 TW/quad), all four beams of a quad overlapping at the aim on R = 1.1 mm; defocus
spot growth κ as in §2.8; l ≤ 30. Published: σ₁ = 6.03 % (OC), 8.56 % (N190204-003);
σ₂ = 12.56 % / 14.83 %; P̄abl = 56.0 / 63.5 Mbar; f = 0.35 / 0.18; absorption
77 % / 82 % and other 3-D ASTER results in Table II. *Missing or ambiguous*: the
best-focus CPP shapes are Ifriit internal models ("NIF-inner-23" … "NIF-outer-50" in
Pol's deck generator), not tabulated (Briareus's LMF-measured fits are the public
substitute); whether "polar angle" in Table I is measured from the pole or about the
port (the caption and the method text differ); whether σ₁ is cos-projected incident
intensity (Pol's `DIAGNOSE_INPUT_BEAMS` at 1100 µm suggests incident intensity, **L**);
the power rescaling after optimization. N190204-003's own pointing is only drawn
(Fig. 1a; quad splitting, stalk and backlighter adjustments), so its σ₁ is not
reproducible. σ₂ needs the plasma and CBET (Pol ships a 1-D-averaged profile under MIT;
which shot and time it is was not verified). Expected use: a documented consistency
check (σ₁ within ~1 % absolute, **L**), and an apples-to-apples optimizer comparison:
Ba24's f with σ₂ dropped (or replaced at Stage 1) evaluated by Briareus at OC and at
Briareus's own optimum, with evaluation counts.

**4. Ye21, NIF 3-mm and 4-mm designs (inputs complete, metric not; H).** Table 1 gives
all 24 parameters for both targets (e.g. 4 mm "Orange": ring-1 A/B θ = 13°/24°,
ring-2 38°/46°, ring-3 53°/71°, ring-4 79°/88°; Δφ 22.5°/11.25°/15°; defocus 2.5, 2.0,
2.0, 2.0 cm), the 750 kJ pulse (linear rise to 357 TW at 0.7 ns, fall 2.3–2.6 ns), 25.4 µm
CH. Published: V = 2.22 % (V_{m=0} = 1.23 % from SAGE, V_{m≠0} = 1.85 % from the 3-D
deposition map via Eq. 3), R_CM rms 11.4 µm after 923 µm travel at 3.1 ns, 94.5 %
absorption. Best-focus spots are Munro's model (private communication), defocused by
Co08's method. Only the m ≠ 0 part is in Briareus's reach, and only once deposition
includes refraction (time-integrated to 3.1 ns). Good as a realistic pointing fixture.

**5. Co08 350 kJ design (inputs nearly complete, metric hydro-only; H).** Table I:
rings 1–4 shifts 50/200/250/594 µm (SiO₂; CH: 100/300/350/694), defocus 2.6/2.6/2.6/1.3 cm,
inner/inner/outer/outer CPPs, split quads ring 3 (200, 400) µm and ring 4 (75, 300) µm;
CPP 50 % contours 739 × 636 µm (inner) and 593 × 343 µm (outer) with an unspecified
radial fall-off; R = 1.1 mm. Published: R_CM rms 6.5 µm and V_rms 1.7 % at 2.8 ns
(SAGE), ±50 µm ring-4 pointing tolerance (Fig. 8). No static illumination number.

**6. Ho15 standard NIF PDD pointing (incomplete; H).** Table I: cone 23.5° unrepointed;
30° → 35°; 44.5° split to 46° and 69° (one beam per focus); 50° → 83° (two per focus);
foci in 22.5° azimuthal steps; R = 1.1 mm. Defocus and spots (indirect-drive CPPs) are
not tabulated; pulse shapes only in Fig. 4. Measured self-emission P₂/P₄/P₆ and DRACO
comparisons are hydro-level. Use: a named fixture and qualitative mode signs.

**7. CM18 ignition-scale designs (spots fully parametric, metric hydro-only; H).**
R₀ = 1482 µm, 1.8 MJ; repointed to rings at 23°, 49° and 84.4° (northern hemisphere;
far-field translations −12.9, 482.5, 950.6/837.3 µm; the table extraction is partly
garbled, **M** on the cone-to-ring assignment), spots as in §2.7, equatorial pulses 50 %
higher. Published metrics are DRACO absorption and implosion speed with and without
detuning. Use: a PDD-plate fixture.

Also useful as *validation* (not benchmarks): the SL83 separability, exact for beams
aimed at the centre — σ_l = |a_l|·G_l with a_l the single-beam Legendre coefficients and
G_l the geometric factor of the port set — on OMEGA's 60 beams (`omega_ports`) with SG4
spots; and converged `raytrace_deposition` (nspot ≥ 121) on the bare sphere.

## 5. Staged implementation plan

**Stage 0 — static hard-sphere illumination (no plasma).**
- Bare-sphere target (icosphere at R) and a surface-aim helper; spherical pointing fields
  (θ_s, Δφ_s per group, optional per-beam split offsets) feeding the existing `aim`.
- Optional absorption weight η_a(ψ) (MN17 Eq. 29), per operator sample.
- Capsule-frame metrics: σ_rms over l ≤ l_max (30), the σ_l spectrum, the m = 0 / m ≠ 0
  split, as `Evaluation`-level functions on the exact-solid-angle Y_lm weights.
- Tests: SL83 separability on centred OMEGA beams (closed form); MN17 re-solved B_l
  against direct integration of Eqs. 28–31; Te14's four σ_2D values; converged ray trace.
- Example/tutorial: Ba24 OC σ₁ with Briareus's LMF CPPs (and with Ba24's linearly scaled
  spots emulated as zero-defocus spots), a 16-parameter σ₁ optimization from the Ho15
  pointing, and the comparison of optimum, evaluations and wall time with Ba24's.

**Stage 1 — static corona.** Spherically symmetric n_e(r), T_e(r), Z̄(r) (analytic
exponential/linear, or tabulated from a 1-D hydro snapshot); refraction and inverse
bremsstrahlung per impact parameter, tabulated as kernels and summed over the beam
samples; absorbed power per steradian and radius, total absorption, and the Ba24 P_abl
and Ye21 0.72-power proxies. Tests: closed-form absorption for linear (cos⁵) and
exponential (cos³) profiles, turning radius n_e = n_c cos²ψ, energy conservation, a
direct 3-D ray integration in the same profile. Benchmarks: MN17 (thin-corona limit),
Ye21's m ≠ 0 part once a comparable profile is available.

**Stage 2 — time dependence.** A sequence of corona snapshots (from the user's 1-D runs,
e.g. CompactLES, or tabulated), per-cone pulse shapes and time-varying power balance,
snapshot-weighted objectives (Ba24 suggests σ < 2 % at all times needs power balance
varied between snapshots), and an optional in-memory 1-D coupling.

**Stage 3 — direct-drive LPI and hardware.** Ray-based CBET with outgoing rays and
wavelength detuning (Ma18, CM18), Saturn rings (Cr05), PDD plates (secondary ellipse,
SMA), square near field, robustness to ±50 µm pointing and power imbalance through
`mode_uncertainty`/`sample_modes`. Ba24's σ₂ becomes reachable only here.

## 6. Open decisions

1. **What Stage 0 is validated against.**
   - (a) The analytic cases (MN17, Te14, SL83) as tests. + Exact, cheap, independent;
     inputs complete. − Idealized (continuous rings, parallel rays); no facility spots.
   - (b) Ba24 σ₁ as the headline check. + The facility case the roadmap names; NIF ports
     and public CPPs exist. − Spot models differ (Ifriit internal vs LMF fits); Table I
     convention ambiguous; agreement can only be approximate.
   - *Recommendation*: (a) as tests, (b) as a documented example comparison with the
     differences measured, never as a fit target.
2. **Defocus model.**
   - (a) Keep Briareus's f-cone geometric optics (Co08's model in substance). + Physical,
     already validated, smooth in defocus. − Disagrees with Ba24's empirical scaling;
     ignores NIF's square near field.
   - (b) Add Ba24's linear radius scaling. + Matches Ba24's numbers more closely.
     − Empirical pF3D fit for one set of plates; a second defocus convention.
   - *Recommendation*: (a); emulate (b) in the Ba24 comparison with explicitly scaled
     zero-defocus spots (no new API). Add a square-aperture option only if Ye21 or Co08
     comparisons show it matters.
3. **Corona model for Stage 1.**
   - (a) Analytic profiles only. + Closed forms, trivial to test. − Not a real implosion.
   - (b) Tabulated 1-D spherical profiles with tabulated ray kernels. + Covers real
     snapshots and the analytic ones as special cases; keeps the operator structure and
     ForwardDiff. − Kernel tables add memory and an interpolation error to control.
   - (c) 3-D ray tracing through 3-D hydro. − Duplicates Ifriit/Mazinisin; breaks the
     precomputed-operator design.
   - *Recommendation*: (b), with (a) as its tests.
4. **Drive metric and objective.**
   - (a) Absorbed (or incident) intensity rms only. + Model-free, matches σ₁. − Ignores
     that oblique light drives less (the reason for PDD equatorial power).
   - (b) Add the P_abl (Ba24, n^{2/3}-weighted, 2/3 power) and Ye21 (0.72 power) proxies.
     + Published, cheap. − Empirical exponents fitted to other codes' runs.
   - *Recommendation*: absorbed power per steradian as the primary output; (b) as
     optional objective transforms with the exponent a parameter and the provenance in
     the docstring.
5. **Parameterization.**
   - (a) New spherical fields in `GroupPointing` (e.g. `:theta_s`, `:dphi_s`, split
     offsets). + One parameterization; existing bounds, inference guards and optimizer.
     − More cases in an already broad type.
   - (b) A separate `SpherePointing`. + Clean separation. − Duplicates mirroring,
     grouping, bounds and the inference tests.
   - *Recommendation*: (a).
6. **Optimizer comparison with "GA optimizers".**
   - (a) Gradient + multistart (existing ForwardDiff/Optim path). + Briareus's
     advantage; reports evaluations to converge. − Not the same algorithm as Ba24.
   - (b) Add a GA dependency. − New dependency for a comparison only.
   - (c) Also run a derivative-free baseline already in Optim (e.g. particle swarm or
     Nelder–Mead) on the same objective. + Algorithm comparison without a new
     dependency.
   - *Recommendation*: (a) + (c); state that Ba24 is the only published GA PDD optimizer
     and that its full objective needs Stage 3 (CBET).
7. **First facility.**
   - (a) NIF 192 beams. + Every quantitative benchmark (Ba24, Ye21, Co08, Ho15, CM18) is
     NIF; ports and CPPs are in the package. − Larger problem.
   - (b) OMEGA 40-beam PDD emulation (`omega_beams` cones 1–3). + Small; Cr05 heritage.
     − No numerical benchmark found (pointing shifts not in the abstracts).
   - *Recommendation*: (a).
8. **Direct-drive CBET.**
   - (a) Defer to Stage 3. + Keeps Stages 0–2 clean and validated. − PDD designs without
     CBET are optimistic (Ho15: 89 % → 67 % absorption).
   - (b) A crude equatorial loss factor sooner. − An uncalibrated number with no data
     to fit it to.
   - *Recommendation*: (a); flag CBET-free results as such in docs and reports.
