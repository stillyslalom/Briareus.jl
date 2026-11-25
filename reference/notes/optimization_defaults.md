# Research notes: Optimizer internals, defaults and cost

This file is one part of the research notes. The index, the confidence levels and the
table of current defaults are in [research-notes.md](../research-notes.md).

## Mode gates in physical units (2026-11-19) *measured*

Roadmap item 36 (the gates sub-bullet). Mode gates set from the pointing-noise floor
(`√(2l+1) × rms_mode` of `mode_uncertainty` at 50 µm rms per quad) state what pointing error
by itself produces: a lower bound for a gate. Laser delivery and the hohlraum plasma move
the shape from shot to shot by more (§ Drive asymmetry inferred from post-shot capsule
simulations: Haines 2026: no yield effect for a fixed P2/P0 of −0.01 to +0.015 on one
design, nominal repeats spanning 0.046, the modelers' repeat uncertainty 0.005), and the
user's judgment (2026-11-18) is that a study's tolerances can be relaxed within that
uncertainty without a loss of performance. This section measures what a relaxed gate buys.

**Unit.** A gate bounds `mode_departure`, `δₗ = √((2l+1) Σₘ (c_lm/c_00 − r_lm/r_00)²)`
against the 48-quad shape; for an axisymmetric difference it is `|Δ(aₗ/a₀)|` in the
Legendre expansion of the capsule flux, so a gate of 0.5% at l = 2 is ±0.005 in P2/P0.
"0.5%" below means 0.005 in every degree l = 1–4.

**Case.** The dropped-bundle tutorial's (bundle 11 freed; the 46 other quads re-aimed in
`:phi` and `:z` within ±0.3 mm; its mesh and coarse deposition operator; scratch
scripts reusing its setup verbatim), with the tutorial's limits as constraints: the
wall-irradiance ceiling at the 48-quad peak (`irradiance_ceiling`, exact) and the
end-cap keep-out (`clearance_floor` of the `HardwareKeepOut` at 10 µm, one row per
quad). Gates at k × the floor (floor 0.346, 0.395, 0.313, 0.095% for l = 1–4) and at
0.25, 0.5 and 1% in every degree. Per gate set: (a) `shift_objective` (the least
re-pointing), (b) `drive_objective` (the largest drive), one call each from the
unshifted start, `iterations = 1000`, `evaluations = 3000`, a time limit of 60–110 s
and `from` to continue a run that reached it (three runs, marked ²). Every run
converged (`:interior`, or `:box` where marked ᵇ) with every constraint within
`feasibility`. The start (the 46 quads at the 48-quad aims) keeps 0.96247 of the
48-quad drive, departs 1.091, 1.075, 0.769, 0.121% in l = 1–4, has its peak at the
ceiling and clears the caps by 66 µm. 8 threads, 40 ms per evaluation (another run may
have shared the machine).

| gates | start δ/gate, l = 1–4 | (a) rms / max shift (mm) | (a) drive | (a) evals | (b) drive | (b) δ/gate, l = 1–4 | (b) rms / max shift (mm) | (b) evals | (b) active besides gates |
|---|---|---|---|---|---|---|---|---|---|
| 0.8 × floor | 3.95 3.40 3.08 1.59 | 0.108 / 0.201 | 0.9586 | 1,924 | 0.9614 | 0.98 1.00 0.95 0.85 | 0.181 / 0.319 | 694 | ceiling |
| 1 × floor | 3.16 2.72 2.46 1.27 | 0.093 / 0.178 | 0.9591 | 712 | 0.9624ᵇ² | 1.00 1.00 1.00 0.88 | 0.204 / 0.331 | 1,558 | ceiling, 1 keep-out |
| 2 × floor | 1.58 1.36 1.23 0.64 | 0.036 / 0.060 | 0.9612 | 642 | 0.9660ᵇ | 0.82 1.00 0.89 1.00 | 0.169 / 0.308 | 1,495 | 6 keep-out |
| 3 × floor | 1.05 0.91 0.82 0.42 | 0.0024 / 0.0036 | 0.9625 | 91 | 0.9694 | 0.67 1.00 0.64 1.00 | 0.175 / 0.319 | 1,102 | 4 keep-out |
| 5 × floor | 0.63 0.54 0.49 0.25 | 0 / 0 | 0.9625 | 7 | 0.9749² | 0.77 1.00 0.22 0.83 | 0.200 / 0.337 | 1,645 | ceiling, 10 keep-out |
| 0.25% | 4.36 4.30 3.08 0.48 | 0.104 / 0.172 | 0.9575 | 705 | 0.9632ᵇ | 1.00 1.00 1.00 0.89 | 0.217 / 0.326 | 1,721 | 6 keep-out |
| 0.5% | 2.18 2.15 1.54 0.24 | 0.063 / 0.105 | 0.9594 | 636 | 0.9674ᵇ² | 1.00 1.00 1.00 0.97 | 0.217 / 0.334 | 1,496 | 9 keep-out |
| 1% | 1.09 1.08 0.77 0.12 | 0.0058 / 0.010 | 0.9622 | 287 | 0.9709 | 0.84 1.00 0.71 0.59 | 0.218 / 0.345 | 882 | 9 keep-out |

Drives are relative to the 48 quads; δ/gate is the departure over the gate (1.00:
active); "evals" are totals over the segments of a continued run. The least re-pointing
ends with every violated degree at its gate (and the ceiling active within
`feasibility` in most rows); at the floor it also holds l = 4 there.

**What it shows.** The least re-pointing falls fast with the gate: 93 µm rms at the
floor, 63 µm at 0.5%, 36 µm at twice the floor, and a few micrometers once the gates
reach the departures of the start itself (1%, three times the floor), none at five times; the
gates a start already meets need no re-pointing, and at 1% only l = 1 and 2 are
violated, by under 10%. That design does not seek drive and loses up to 0.5% of the
start's. The largest drive keeps rising across the whole range, from 0.9614 at 0.8 ×
the floor to 0.9749 at five times, with the l = 2 gate active in every run (and l = 4
at two and three times the floor): the drive is bought with P2, about 0.8–1% of the
48-quad drive per percent of l = 2 departure on this untuned pointing. No return stops
within the range measured: besides l = 2 the drive runs press against the cap keep-out
(4–10 quads at 10 µm from twice the floor and 0.25% on) and the ±0.3 mm box, and move
the quads 0.17–0.22 mm rms whatever the gate. A tolerance of 0.5% per
degree, the modelers' P2 repeat uncertainty, buys 0.5% of drive over the floor gates
and needs a third less re-pointing to meet. Gate sets that are not nested (0.25%
against 0.8 × the floor, whose l = 4 gate is a third as wide) are not ordered in drive.
Confidence: one public case and one start per run; the drive runs are local optima of a
landscape where the active gates leave a flat family (differences of a few 1e-4 between
neighboring designs are not significant).

**The reported multipliers of the drive runs.** Every largest-drive run ends with all
gate multipliers zero (`PointingResult.constraints`, `state.λ`), the gates active at
normalized violations of −1e-5 to −6e-4 and the penalty parameter at 1e3–1e4; the
least-shift runs end with positive gate multipliers (2e-4 to 6e-3) at μ of 7–1e3. With
μ that large, the last update `λ ← max(0, λ + μg)` at a design a few 1e-4 inside a gate
sets its multiplier to zero, and a drive objective's last inner run, flat along the
limits, converges on the relative value change. The run of the tutorial within limits
(ceiling and floor gates, no keep-out) shows it: it reports multipliers 0, 0, 0, 1.1e-3,
0 (ceiling, l = 1–4) at μ = 204; at its design the drive objective's gradient over the
90 free parameters has norm 1.75e-2, the reported multipliers leave a stationarity
residual `|∇f + Jᵀλ|` of 1.80e-2, and least-squares multipliers from the same
ForwardDiff Jacobian (ceiling 5.7e-3, l = 2 2.4e-3, l = 4 1.2e-3, the others under
1e-4) leave 8.8e-4. The design is near a KKT point and the run's drive is the
tutorial's, but the multiplier column of a drive run does not show which limit costs
drive, and `bound_sensitivity` of such a run is close to `bound_gradient`. Not changed
here (a change to the multiplier estimate the augmented Lagrangian reports).

## Multipliers from the converged active set (2026-11-19) *measured*

Roadmap item 36 (a drive run under gates reports zero gate multipliers). The augmented
Lagrangian of `optimize_pointing` updates its multipliers after each inner run by
`λ ← max(0, λ + μ g)`. That update is the multiplier of the inner problem when the
inner run ends at a stationary point of the augmented Lagrangian; an inner run that ends
elsewhere leaves an error of μ times its error in each `g`. The inner runs end on a
projected gradient below the tolerance or on a relative value change below
`_AL_F_RELTOL` = 1e-10 per iteration (§ "Constrained re-pointing": without it a drive
objective's inner run creeps for hundreds of iterations), and the outer loop counts
either as converged. At the penalties a drive run reaches (μ of 1e2 to 1e6, grown
because the infeasibility–complementarity measure stalls with the inner runs), a limit
the design ends a few 1e-5 inside gets the multiplier max(0, λ − μ · 1e-5) = 0.

**Closed form** (scratch scripts; `test/constraints.jl`, "multipliers from the
converged active set", runs the first two rows below): 48 quad powers, objective
`-c·P + ε/2 |P − 8|²` with `cᵢ = (1 + 0.3 sin i)/384` (the inner cones' raised by 50%
where stated) and ε = 2e-4, under ceilings on the total power
(1.01 × 384 TW) and on the inner cones' (130 TW) and an inactive floor on the 50° cone's;
the Karush–Kuhn–Tucker point and multipliers in closed form. Inner runs from the third
outer iteration on ended on the value change after 1–3 L-BFGS-B iterations at
projected gradients of 3e-5 to 6e-5, against tolerances of 1e-5 to 1e-8.

| case | design error | f − f* | iteration's λ (total, inner) | reported λ | KKT λ |
|---|---|---|---|---|---|
| inner cones' `c` +50%, from the start | 7.0e-5 | 1.0e-5 | 1.00893, 0.17363 | 1.00892, 0.17364 | 1.00890, 0.17364 |
| same, restarted at μ = 10⁶ | 6.8e-5 | 3.1e-7 | 1.22289, **0** | 1.00890, 0.17364 | 1.00890, 0.17364 |
| plain `c`, inner cap 129.2 TW | 1.3e-3 | 1.6e-6 | 1.01011, **0** | 1.00696, 0.00389 | 1.00696, 0.00388 |
| plain `c`, inner cap 129.1 TW | 1.8e-3 | −2.8e-6 | 1.01661, **0** | 1.00665, 0.00419 | 1.00672, 0.00412 |
| plain `c`, inner cap 130 TW | 7.9e-3 | 3.4e-5 | 1.01352, 0 | 1.01085, 0 | 1.00890, 0.00195 |

With the inner cap at 129.1, 129.2 and 129.5 TW the run ends 2e-5 to 2e-4 inside the
inner limit, and the iteration's multiplier for it is zero where the closed form's is
3–4e-3 (at 129.3 TW it ends 4e-6 outside, with a multiplier 8% high).
At 130 TW the run ends 0.66% inside that limit (outside the `feasibility` band) with
the objective 3.4e-5 above the optimum: the design itself, not only the multiplier,
is short of the optimum, which no multiplier estimate at that design repairs.

**Fix: the report** (`_active_set_multipliers`). At the end of a constrained run the
reported multipliers (`PointingResult.constraints`) are the nonnegative least-squares
solution of `∂f/∂xᵢ + Σⱼ λⱼ ∂gⱼ/∂xᵢ = 0` over the parameters inside the box, on the
components within `feasibility` of their (working) limits or beyond them, zero for the
others (Lawson–Hanson; a weight 1e-6 × the largest column norm on |λ| picks the
least-norm solution where several fit equally: the twelve symmetric beams of the
keep-out test, whose plain NNLS solution put 7.49 on one beam, each get a positive
multiplier, as the iteration's twelve of 0.624).
One Jacobian of the objective and the constraints, a gradient's cost, after the run.
`bound_sensitivity` is `∂L/∂xᵢ` with these multipliers, and the new field
`PointingResult.stationarity` is the largest `|∂L/∂xᵢ|` over the parameters inside the
box. The multipliers of the iteration stay in `state.λ` for `from`. Chosen over a change
to the iteration because the iteration is ALGENCAN's and its update is exact for inner
runs that end stationary; what departs is the inner value tolerance, a default (below).
Designs, values, evaluation counts and constraint values are unchanged: all 27 runs
of `test/constraints.jl` (its version before this change, every call logged) and the
tutorial's run within limits (`x` compared with `==`).

**The tutorial's run within limits** (ceiling and pointing-floor gates; 1,079
evaluations, drive 0.962575 of the 48 quads, both before and after). Normalized
violations at the end: ceiling (stand-in at its working limit) −1.3e-4, l = 1–4
−4.3e-4, −3.2e-5, +5.3e-6, −5.4e-3; μ = 204.

| multipliers | ceiling | l = 1 | l = 2 | l = 3 | l = 4 | residual ‖·‖₂ (free) | ‖·‖∞ |
|---|---|---|---|---|---|---|---|
| none (the drive's gradient) | | | | | | 1.75e-2 | 3.05e-3 |
| iteration's (reported before) | 0 | 0 | 0 | 1.08e-3 | 0 | 1.80e-2 | 3.44e-3 |
| reported now (band-active NNLS) | 3.13e-2 | 8.3e-5 | 2.68e-3 | 1.89e-4 | 0 | 5.65e-3 | 1.21e-3 |
| least squares, all five, signs free | 5.7e-3 | 1.6e-5 | 2.38e-3 | 6.5e-5 | 1.22e-3 | 8.8e-4 | 2.1e-4 |

The reported multipliers name the limits that cost drive (the ceiling and l = 2) where
the iteration's named none, but they leave a third of the drive's gradient: the design
is short of stationarity, and the five-constraint fit does better only by giving l = 4,
0.54% inside its gate, a multiplier. The bound sensitivities move from 2.0e-6 and
−3.6e-5 (Q12T.phi at its upper and Q14B.phi at its lower bound; both of the wrong sign)
to −1.5e-3 and +1.4e-3. Restarts through `from` (scratch script): in the same box the
run leaves both bounds in 46 evaluations for a gain of 1.6e-7; with either bound
widened by 0.01 or 0.02 mm it gains under 6e-9 (the same design); with a bound
narrowed, the drive falls by amounts that are not monotonic in the step (Q12T.phi:
1.4e-4 at 0.01 mm, 1.1e-5 at 0.02 mm), the differences between neighboring designs of
the flat family. The bounds do not hold these parameters, and no first-order
derivative at this design predicts that: the overstatement of a wider box's gain
recorded on a re-pointing study (§ "Constrained re-pointing", "On the study") has the
same cause, a design that ends short of stationarity, and the corrected multipliers do
not remove it. `stationarity` now shows it on every result.

**User decision: the inner value tolerance.** On the closed form (plain `c`, inner cap
130 TW, default options otherwise):

| `_AL_F_RELTOL` | evaluations | outer | design error | f − f* | inner limit (normalized) | `stationarity` | iteration's inner-cap λ |
|---|---|---|---|---|---|---|---|
| 1e-10 (default) | 219 | 7 | 7.9e-3 | 3.4e-5 | −6.6e-3 | 3.5e-5 | 0 |
| 1e-12 | 204 | 7 | 2.3e-4 | 5.5e-7 | −1.6e-5 | 5.7e-7 | 0 |
| 1e-14 | 206 | 7 | 3.5e-6 | −2.5e-8 | −4.4e-9 | 8.5e-9 | 1.95391e-3 |
| 0 (off) | 206 | 7 | 3.5e-6 | −2.5e-8 | −4.4e-9 | 8.5e-9 | 1.95391e-3 |

(KKT: 1.95409e-3.) The negative f − f* is the run ending 2.5e-8 outside the total-power
cap, within `feasibility`.

On the tutorial's run within limits (8 threads, budget 3,000 evaluations; another run
shared the machine, so evaluation counts are given in place of times):

| `_AL_F_RELTOL` | evaluations | end | drive / 48 quads | `stationarity` | reported λ (ceiling, l = 1–4) | iteration's λ |
|---|---|---|---|---|---|---|
| 1e-10 (default) | 1,079 | converged (box) | 0.962575 | 1.2e-3 | 3.1e-2, 8e-5, 2.68e-3, 1.9e-4, 0 | 0, 0, 0, 1.08e-3, 0 |
| 1e-12 | 2,308 | converged (box) | 0.962609 | 1.9e-5 | 2.58e-3, 1e-5, 2.38e-3, 4.0e-5, 1.27e-3 | 2.63e-3, 0, 2.27e-3, 8.2e-5, 1.28e-3 |
| 1e-14 | 3,000 | evaluation budget | 0.962611 | – | 2.40e-3, 1e-5, 2.38e-3, 3.9e-5, 1.27e-3 | 2.53e-3, 0, 2.40e-3, 3.9e-5, 1.27e-3 |

At 1e-12 every gate ends active (the l = 4 gate, 0.54% inside at the default, at its
limit), the exact peak at 1.0010 of the ceiling (within `feasibility`), the fit leaves
0.3% of the drive's gradient (2-norm over the free parameters) and the multipliers of
the iteration come within 5% of the fit's on the three limits that bind (ceiling, l = 2,
l = 4), at 2.1 times the evaluations for 3.5e-5 more drive (inside the 1e-4
differences of the flat family). Five parameters end at a bound; for four of them
restarts with the bound narrowed by 0.01 and 0.02 mm cost 1.4e-4 to 6.4e-4 of the
objective per mm, not monotonic in the step, against reported sensitivities of 1.1e-4
to 1.9e-4 (the fifth's restarts reached the budget); with a bound widened the restart
(38–60 evaluations, on the same value tolerance) leaves the parameter at the old bound
and gains under 6e-7 per mm. A
first-order sensitivity on this flat family is an order-of-magnitude estimate of the
cost of narrowing, and a restart is the measure of what widening gains.

Recommendation: keep 1e-10. The designs agree within the flat family's spread at half
the evaluations, and the reported multipliers with `stationarity` now state how far a
result is from a Karush–Kuhn–Tucker point. A run whose multipliers or bound
sensitivities matter can be continued with `from` at a tighter inner tolerance if the
tolerance becomes a keyword (not one now; the user's decision).

Confidence: the closed forms are exact; the tutorial's case is one start and one mesh.

## The multi-phase study's hand-written pieces as API (2026-11-19) *derived*

Roadmap item 34 (the bullet on the quad stage rerun on `PhasedScenario`). A multi-phase
re-pointing study (six phases, 184 beams in 46 quads, 276 parameters) still wrote eight
pieces by hand on the phased scenario; each is now public, and the study's start
rewritten on them (scratch runs, one evaluation, no optimization) reproduces the
hand-built quantities.

**Choices.**

- *Drive band*: `drive_ceiling(reference; fraction)`, the upper side of `drive_floor`
  with the same reference, scale and laser–plasma-model and time check. Chosen over an
  upper limit inside `drive_floor` because a `Constraint` has one side: a band is two
  constraints, and two builders compose with `in_phase` and the reports without a
  vector-of-constraints return type.
- *Gates*: `mode_gate(...; fraction)` scales a `mode_uncertainty` floor; with `floor`,
  the gate per degree is `max(gates[l], fraction × floor[l])` (the floor as a lower bound
  of physical gates, § Mode gates in physical units). `fraction` without a floor is an
  error, so that `fraction` has one meaning; it is never a scale of physical gates.
- *Shifted target*: `mode_gate(...; shift = l => Δ)` (or `(l, m) => Δ`) measures the
  departure from the reference's coefficients with `c_lm → c_lm + Δ c_00/√(2l+1)`, the
  reference's model and time still checked; `shift_modes` gives the coefficients for an
  objective. For m = 0, Δ is the change of aₗ/a₀ in the Legendre expansion: with
  orthonormal harmonics aₗ/a₀ = √(2l+1) c_l0/c_00, checked to 1e-10 against the Legendre
  moments integrated on the same sub-triangles (the cached `flux_moments` and `flux_ylm`
  integrate at different subdivision levels and differ by 0.2% in P2/P0 on a capsule
  of 0.9 mm facets; the relation itself is exact). `drive_objective` has no reference
  shape, so the shift has no place there.
- *Keep-outs over phases*: `clearance_floor` evaluated on a `PerPhase` returns per beam
  the least clearance over the phases. The phases share the aims; a keep-out at a
  `peak_fraction` gives each phase the same clearance (the value is exactly that of one
  phase), one at an absolute `intensity` widens with the power, and the brightest
  phase then binds.
- *Power scales*: `power_scales(ps, ev | beams | x)` reads each (key, phase) scale as a
  beam's power over its phase's power (the first beam of the key with positive phase
  power; 1 for a key dark in its phase, as `parameters` reads it), differentiable;
  `phase_powers` and `power_scale_table` for reports. A scaled phase on the powers of
  the pointing has nothing to divide by: an error, the scales then come from `x`.
- *History phases*: `remake(scenario; history, time)` and `Phase(...; history, time)`
  are exactly `history_state` followed by `remake(...; time, albedo)`, for a phase
  and its reference alike.
- *Summary*: `pointing_summary(...; gates, floor)`; a `mode_gate` passed as the gate
  reports the value of that gate (the departure from a shifted target where it has one), and
  the run's constraints passed as `gates` select each phase's `in_phase` gate.
- *Ceiling on another target*: `picket_ceiling(...; reference_names)` accepts a
  reference on another target for the three integral measures. `:average` is absorbed
  power over a fixed area, independent of the mesh up to the deposition operator's
  sampling; `:weighted` (∫I² dA/∫I dA) converges as O(h²), about 1% low at h = 0.2 mm
  (§ The published picket fluence), a bias that cancels between targets meshed at the
  same h; `:corridor` is an area mean over a band. `:peak` stays refused: it follows one
  element (8% between 0.2 and 0.45 mm elements, § Constrained re-pointing) and its
  smooth stand-in is scaled by the reference mesh's peak-to-norm ratio.
  `irradiance_ceiling` holds only the peak and keeps its refusal. With `lpi = true`
  the reference's model and time are checked on another target (the same target's path
  is unchanged: it has never checked the model; extending the check there would make
  existing calls throw and is left to the user).

**The study's start on the API** (scratch, coarse deposition operator, 8 threads).
Variant with windows from a drive history: the phase and reference scenarios from
`remake(...; history, time)` have the hand construction's albedos, times and moment
responses exactly; the floors from `mode_uncertainty` equal the study's cache; the
objective with `power_scales` and the floor vectors is 378.899827494, the study's
construction's value exactly; the 212 constraint rows (184 of them the keep-out's,
now without a phase in their names) agree with the hand-built constraints' within
1.2e-16 relative in value and limit, with equal statuses (`drive_ceiling` reads the
drive itself where the study's cap read its ratio: the ratio of value to limit is
compared), and match the study's logged start to every printed digit; the picket
ceilings with `reference_names` on the design's windowed target reproduce the study's
`limit = picket_irradiance(reference, ...)` route. Variant without windows (one target
per phase, no history): objective 330.041515363 in both constructions, the same 212
rows within 1.2e-16 and equal statuses; the phased summary with the run's constraints
as `gates` and the `mode_uncertainty` results as `floor` reproduces the study's
logged departures, multiples of the floor and of the gate in every phase.

Confidence: the API's arithmetic is exact or checked against closed forms in
`test/phased_study.jl`; the study comparison is one start per variant.

## Inner-run options, a wall-loss scenario warning and laser delivery figures (2026-11-19) *measured*

Seven changes the user approved on 2026-11-19 from the recommendations of the sections
above, each opt-in, a warning, a check on misuse or a display: no default and no
existing result changes.

### 1. The inner value tolerance as a keyword

`optimize_pointing(...; inner_reltol = 1e-10)` exposes the relative value tolerance
that ends a constrained run's inner runs (formerly the internal `_AL_F_RELTOL`; phase 1
uses it too; 0 turns it off). Default runs are unchanged: the closed forms of
`test/constraints.jl`, "multipliers from the converged active set", give the same
design with the keyword given as 1e-10 (tested). The measurement behind it is
§ "Multipliers from the converged active set": on the dropped-bundle tutorial's run
within limits 1e-12 reaches a stationarity of 1.9e-5 (1.2e-3 at the default) for 2.1
times the evaluations, the drive the same within the flat family's 1e-4.

A converged result continues tighter through `from` (scratch script, `-t 1`; the same
closed form, 48 quad powers under a total and an inner-cone cap; design error relative
to the KKT point):

| case | run | evaluations | design error | `stationarity` | reported inner-cap λ (KKT) |
|---|---|---|---|---|---|
| inner cones' `c` +50% | default | 162 | 7.0e-5 | 2.1e-7 | 0.173637 (0.173643) |
| | `from` it, default | 6 | 6.8e-5 | 2.1e-7 | 0.173643 |
| | `from` it, `inner_reltol = 0` | 19 | 3.5e-8 | 5.4e-12 | 0.173643 |
| plain `c`, cap 130 TW | default | 219 | 7.9e-3 | 3.5e-5 | 0 (1.95409e-3) |
| | fresh, 1e-12 | 204 | 2.3e-4 | 5.7e-7 | 1.95903e-3 |
| | `from` it, default | 121 | 3.1e-3 | 7.6e-6 | 1.95390e-3 |
| | `from` it, 1e-12 | 121 | 3.1e-3 | 7.6e-6 | 1.95390e-3 |
| | `from` it, 1e-14 or 0 | 189 | 6.2e-7 | 1.6e-11 | 1.95388e-3 |

A continuation runs one more inner run at the final tolerance from the result's
multipliers and penalty, so even the default tolerance moves a result that ended short
of stationarity (the plain case's error falls from 7.9e-3 to 3.1e-3); with the value
tolerance off it reaches the KKT point in 189 evaluations, fewer than a fresh run takes.
`test/constraints.jl` checks the continuation with the tolerance off against the closed
form (design within 1e-5, multipliers within 1e-3), which any converged run of that
strictly convex problem meets. Confidence: exact closed forms; one start.

### 2. The inner budget: an evaluation cap per inner run

`iterations` caps an inner run's L-BFGS-B iterations, not its evaluations: a line
search can take dozens, and an inner run of the multi-phase study spent hundreds of
evaluations on 1e-4 of the objective while limits stayed violated (§ "Re-pointing
tools: the inner budget, the smallest shift, a pulse-shape budget"). The opt-in
keyword `inner_evaluations` (by default `nothing`, so no cap and unchanged results) ends an
inner run at that many evaluations, at its best design, as the iteration cap ends one:
the outer loop updates the multipliers and the penalty, and a capped run never counts
as converged. Phase 1 is not capped (it has no multipliers to update); the
relaxations' constrained runs are.

**Coarse dropped-bundle case** (`test/constraints.jl`, "dropped bundle: gates under a
ceiling in one call": 46 quads, 92 parameters, the exact wall-irradiance ceiling and
four pointing-floor gates, largest drive, `iterations = 1000`, `evaluations = 3000`;
8 threads, another run may have shared the machine):

| cap | evaluations | outer | end | drive / 48 quads | `stationarity` | peak / limit | evaluations per inner run |
|---|---|---|---|---|---|---|---|
| none (default) | 1,044 | 7 | converged, interior | 0.962592 | 3.3e-5 | 0.99735 | 132, 54, **795**, 52, 4, 3, 4 |
| 100 | 766 | 10 | converged, interior | 0.961987 | 8.2e-4 | 0.99922 | 100, 82, 100, 100, 43, 65, 72, 51, 100, 53 |
| 30 | 895 | 30 | `outer` reached | 0.960894 | 1.5e-3 | 1.00141 | 30 (27–30) each |

The default run's third inner run took 795 evaluations; a cap of 100 converged with
27% fewer evaluations at 6e-4 less drive (more than the 1e-4 spread of the flat family
measured in § "Multipliers from the converged active set", so the cap moved the design,
not only the path), and a cap of 30 ran out of outer iterations with the ceiling 0.14%
over. A small cap therefore needs a larger `outer`; the cap's effect is case-dependent,
and the default stays none.

**Closed form** (the multipliers testset's 48 quad powers; `-t 1`; `outer = 100`):

| case | cap | evaluations | outer | end | design error | longest inner run |
|---|---|---|---|---|---|---|
| inner cones' `c` +50% | none | 162 | 7 | converged | 7.0e-5 | 43 |
| | 50 | 162 | 7 | converged (unchanged; no run reaches it) | 7.0e-5 | 43 |
| | 20 | 97 | 7 | converged | 7.0e-5 | 20 |
| | 10 | 75 | 8 | converged | 4.8e-4 | 10 |
| | 5 | 500 | 100 | `outer` reached | 0.34 | 5 |
| | 3 | 300 | 100 | `outer` reached | 0.37 | 3 |
| plain `c`, cap 130 TW | none | 219 | 7 | converged | 7.9e-3 | 43 |
| | 20 | 202 | 12 | converged | 7.8e-3 | 20 |
| | 10 | 164 | 17 | converged | 5.9e-3 | 10 |
| | 5 | 500 | 100 | `outer` reached | 0.27 | 5 |

A cap of 3 or 5 evaluations never lets an inner run reach the final tolerance (the first
L-BFGS-B iterations take that many), so the run cannot converge whatever `outer`; from
10 up the capped runs converge, here with fewer evaluations than the default.

`test/constraints.jl` checks on the closed form that each inner run stays within the
cap (the callback's evaluations per outer iteration), that the capped run converges to
the KKT point, and that a cap above every inner run of the default run leaves it
unchanged. Not run: the multi-phase study itself under a cap (its runs take hours; scratch).

### 3. A warning for a scenario on a history-dependent wall without albedos

A `Scenario` with a response (observed elements and a model) built without `albedo` on
a target whose materials include a `WallLossAlbedo`, or a `BurnthroughAlbedo` whose
`before` is one (the elements whose loss `drive_history` follows through their own
temperature history), takes the loss at the model's reference temperature (250 eV by
default) at `time`, as a wall held at it since the pulse start. On the public mock-up
of § "Drive histories: sliver steps and the albedo of a phase" that overstates the
capsule drive against the history by up to 1.64 times mid-pulse (the second step) and
understates P2/P0's change under a pointing perturbation by 11–37%. The construction
now warns once (one `@warn` per scenario built, naming the component, the reference
temperature and `time`, and pointing to `history_state` and `remake(scenario; history,
time)`). `albedo = :reference` marks the reference temperature as intended: the same
scenario (moment responses identical, tested), no warning; `remake` takes it too.
Chosen over a refusal (the user's option of 2026-11-18): the reference temperature is
a legitimate fixed-time model (a peak scenario, a test), and a warning keeps existing
scripts running. No warning is given for a scenario without a response (deposition only), one
with explicit albedos, a `remake` (the scenario it copies warned at construction), and
`GapClosure`/`LipClosure` closed albedos (`drive_history` does not follow those
through a wall history either).

To find the reach of the warning, every file in `src/`, `test/`, `docs/literate/`,
`examples/` and `BriareusGUI/` was searched that both builds a `Scenario` (or a `Phase`
from a target) and uses a wall-loss albedo, a `BurnthroughAlbedo` or a preset built on
one (`nif_gold`, `uranium`, `nif_uranium`, `du_wall_loss`). Three tests use the reference
temperature by design and now pass `albedo = :reference`: `test/history_state.jl` (a
scenario before `remake` with a state), `test/large_models.jl` (a scenario on the large
path) and `test/phased_study.jl` (the shared scenario whose phases take the history's
albedos). The scenarios of the diagnostic-windows tutorial already pass the albedos of a
solution. The GUI builds its materials from constant albedos, and the precompile
workload does not build a scenario on its wall-loss target. A full run of the
dropped-bundle tutorial as a script gives no warning.

### 4. The laser–plasma model check of `picket_ceiling` on the reference target

`picket_ceiling(...; lpi = true)` deposits the evaluation's effective beams, so its
limit depends on the laser–plasma model and the time at which the reference was
evaluated. Until now the model was checked only on an evaluation of another target
(§ "The multi-phase study's hand-written pieces as API"); on the reference target
a design evaluated with another model was measured against the reference's limit
silently. The check now runs on either target, as for `irradiance_ceiling` and
`mode_gate`: another model, none against one, or a model at another time throws an
`ArgumentError`; a model built alike is the same model. `lpi = false` (the default)
deposits the setup beams and checks nothing. Every existing call in the tests, the
tutorials and the examples evaluates the ceiling with the model of the reference and
passes; a call that mixed models now throws (the user's decision of 2026-11-19).
`test/picket_fluence.jl` covers each case on the reference's target.

### 5. A note in the result display for an exact ceiling's multiplier

An exact ceiling (`irradiance_ceiling`, `picket_ceiling` with `exact = true`) reports
the exact peak in its row while the multiplier fit uses its smooth stand-in at the
working limit (§ "Multipliers from the converged active set"): the row can read
`:satisfied` with a positive multiplier (1.2% inside the limit under one CI rounding,
the case of `test/constraints.jl`'s dropped-bundle testset). The `PointingResult`
display now adds one line naming such constraints and saying that the multiplier is
the stand-in's at its working limit; `ConstraintRow`'s fields and the result are
unchanged. Display only.

### 6. A tolerance for the number of pulse shapes

`pulse_shape_groups(...; shapes, tolerance)` returns the smallest number of shapes,
from one per class up to `shapes`, whose best split between the classes fits the
per-key multipliers with an rms in log multiplier at or below `tolerance` (`shapes`
when none does); `counts` lists the best rms at each number tried. Without `tolerance`
the grouping is the fixed count as before (unchanged, with the same groupings per class and
the same splits). Each class's best grouping per count is computed once for the
largest count already, so the search adds only the splits of the smaller counts.
Checked in `test/pulse_shapes.jl` against exhaustive enumeration of every assignment
and split at each count.

**Scale, not a default** (public facility figures; the model's parameters are not
public): the NIF User Guide, section 4.8, "Laser Performance Operations Model"
(https://nifuserguide.llnl.gov/home/4-laser-system/48-laser-performance-operations-model,
read 2026-11-19), states that the laser performance model predicts total energies
"within 2–3 percent" and provides "beam-to-beam energy balance better than 5 percent
for full system shots" (the page does not state whether the balance is rms or extreme).
A relative error ε is ε in the natural log to first order, so tolerances of 0.02–0.05
span those figures: shapes whose multipliers differ by less than the laser delivers
distinctly need not be separate pulses. The tolerance is the study's choice.
Confidence: the figures are a public statement of the facility itself, **M**; how they map
onto a per-shape rms is an engineering reading.

### 7. Laser delivery in the noise floor

`mode_gate`'s floor from `mode_uncertainty` counts only pointing errors by default (50
µm rms per pointing unit, NIF's specification); `mode_uncertainty` also takes `power`,
an independent relative rms per pointing unit. The public balance figure of § 6
("better than 5 percent" beam to beam) puts a number on what the floor leaves out, the
variation in laser delivery that the user's review of 2026-11-18 named as moving the
shape more than pointing does. On the dropped-bundle tutorial's case (48 quads at
8 TW, its coarse deposition operator; the floor per degree √(2l+1) × `rms_mode`, as
`mode_gate` takes it; a scratch run of the coordinator's script, rerun here to the
digits shown; 8 threads):

| errors per quad | floor l = 1 | l = 2 | l = 3 | l = 4 (%) |
|---|---|---|---|---|
| pointing 50 µm (the tutorial's floor) | 0.346 | 0.395 | 0.313 | 0.095 |
| power 2% alone | 0.227 | 0.179 | 0.107 | 0.024 |
| power 3% alone | 0.340 | 0.268 | 0.161 | 0.036 |
| power 5% alone | 0.567 | 0.447 | 0.268 | 0.060 |
| pointing 50 µm + power 2% | 0.413 | 0.434 | 0.330 | 0.098 |
| pointing 50 µm + power 3% | 0.485 | 0.477 | 0.351 | 0.102 |
| pointing 50 µm + power 5% | 0.664 | 0.596 | 0.412 | 0.112 |

The pointing unit here is a quad: 5% rms per beam over four independent beams is 2.5%
per quad (the sum's relative rms falls as 1/√4), between the 2% and 3% rows; 5% per
quad is the figure read as the balance within one quad. The power term adds in quadrature,
mostly at l = 1 and 2: at 3% per quad the floor rises by 40% at l = 1 and 21% at
l = 2, at 5% by 92% and 51%; l = 4 moves by at most 18%. No default changes; the
floor stays pointing-only, and the `mode_gate` docstring, where it presents the floor
as a lower bound, points here for what the delivery figure adds (`power` in
`mode_uncertainty`). Confidence: the propagation is linear and exact for the model's
sensitivities (`mode_uncertainty`); the balance figure's statistic is not stated.

## Defaults changed (2026-11-19) *measured*

The user's decisions of 2026-11-19 on defaults listed as open in roadmap items 17, 34, 35
and 36 ("Agreed with all recs"), one change each, in the order made. Each subsection
gives the before and after numbers on the cases that move and states which do not.

### 1. `Conversion()`'s default efficiency 0.9

`Conversion(; efficiency)` defaulted to 0.75 after `gold()` and `uranium()` had moved to
0.9 ("Material defaults and their sources", "Decisions of 2026-11-18"). It now defaults
to 0.9, `gold()`'s conversion and the input of Young 2020's static VisRad model of a NIF
hohlraum (**H**); 0.75 stays available as `Conversion(efficiency = 0.75)` for 1 ns Nova
and Omega hohlraums. `nif_conversion()` sets a separate efficiency (0.7 t^0.12) and an
intensity table replaces the constant, so neither moves.

No call site in `src/`, `test/`, `docs/literate/`, `examples/` or `BriareusGUI/` built a
`Conversion` on the default efficiency: every one passes `efficiency` or an intensity
table, so no existing result moves and no call site is pinned. Check by the closed
form (*derived*): the default-size hohlraum (`hohlraum(length = 10.1, radius = 2.875,
leh_radius = 1.55, capsule_radius = 0.9)`, h 0.8 / capsule_h 0.3) with a wall of albedo
0.8 and 384 TW absorbed uniformly over the wall, the wall's conversion a constant
`Conversion` at the old and the new default:

| quantity | `Conversion(efficiency = 0.75)` (old default) | `Conversion()` | check |
|---|---|---|---|
| capsule mean flux (TW/mm²) | 4.581977 | 5.498372 | ratio 1.200000000000000 = 0.9/0.75 |
| capsule Tr (eV) | 258.3648 | 270.4136 | (0.9/0.75)^¼ × 258.3648, difference 0 eV in Float64 |

`test/energetics.jl` checks that the default equals `gold()`'s conversion and that a wall
on `Conversion()` solves exactly as one on the number 0.9.

### 2. `subset` keeps a `GroupPointing`'s groups

`subset(scenario, keep)` without a parameterization rebuilt a `GroupPointing` scenario's
parameterization on the kept beams grouped by `:group` (the cone), whatever the
grouping of the pointing; the composed parameterizations (`HeldParameters`,
`PowerScales`, a `PhasedScenario`'s pointing) already kept theirs ("Re-pointing studies
on the public API"). Now every default subset keeps the groups: a group that keeps a
beam keeps its key, so the subset's parameter names are those of the full pointing, with
its fields, mirroring, `relative`, `unit` and sphere (one function, `_subset_param`,
where there were two). A pointing grouped by `:group` gets the same groups in the same
order as before, exactly.

What moves (every case found): a pointing grouped otherwise and subset without a
parameterization. In `test/pointing_summary.jl` the per-quad pointing (`by = :quad`, `:z`
and `:phi`, unmirrored, relative) of 44 kept quads runs one iteration on 88 parameters
instead of 8 (one pair per quad, not per cone); the test compares two summaries of the same
result and passes. In `test/basics.jl` (`by = :quad`, `:power`) and
`BriareusGUI/test/test_halfraum.jl` the subset is evaluated on given beams, so its
parameterization does not enter. The tutorials pass a parameterization to every `subset`
they optimize on, and `examples/n221204.jl` does not build a subset: neither can move.
`test/parameterizations.jl` checks the default against a `GroupPointing` built by hand on
the kept beams by quad (names and beams equal).

### 3. A smooth minimum over the rim in `ApertureKeepOut` (*engineering*)

The aperture keep-out's clearance was the exact minimum over the rim (the samples and
the four sample intervals beside the closest ones, refined by golden section), with
genuine kinks where the worst rim point jumps between separate local minima of the
distance ("Optimization defaults decided", "The aperture keep-out's slope along the
aim"). It is now a smooth minimum, `ApertureKeepOut(...; softness = 1e-4)` (mm):

- every sample interval whose nearer end lies within W = max(3 µm, 30 softness) of the
  closest sample is refined by golden section; its value is the least of the golden
  result and its two ends (on an interval where the distance is not unimodal, golden
  section can end at the higher end and switch ends as the beam moves: a jump, seen
  at 4 nm before the ends were included);
- the clearance is −λ log Σₖ σₖ exp(−gₖ/λ) (λ = `softness`), the interval's weight σₖ 1
  inside W and falling C¹ (smoothstep) to 0 at 2W, so intervals enter and leave the sum
  continuously;
- it never exceeds the exact minimum (the sum holds the interval of the exact minimum at
  weight 1), lies below it by at most λ log(number of intervals counted), and is the
  exact minimum for `softness = 0`, which reproduces the former clearance exactly
  on the tutorial's 48 quads and the 480 random ones below.

W: the distance dips between two samples by at most 0.51 µm below the nearer one on
the cases below (96 rim samples), so 3 µm keeps every interval that can hold the
minimum; in the taper band the weights are below e⁻²⁵.

**Superseded in part** ("The aperture keep-out's exact interval minima (2026-11-20)"):
each interval's minimum is now found exactly from the distance's slopes at the samples,
with no window, and the weights σₖ are anchored on the least interval minimum (1 within
10 λ, 0 at 20 λ). The smooth minimum over intervals, the softness and every number below
are unchanged: the clearances move by at most 1.2e-7 µm.

Measured on the dropped-bundle tutorial's hohlraum and keep-out (1% contour, minimum
10 µm, as "Optimization defaults decided"): its 48 quads as pointed there, the same
quads moved to within ±20 µm of the minimum along random directions, and 480 quads with
random aim shifts up to ±0.3 mm; exact minus smooth clearance, in µm (*measured*):

| softness (mm) | 48 quads | 48 near the minimum | 480 random |
|---|---|---|---|
| 3e-5 | 0.003–0.041 | 0.000–0.040 | 0.000–0.027 |
| **1e-4** (default) | 0.059–0.171 | 0.003–0.176 | 0.001–0.131 |
| 3e-4 | 0.29–0.62 | 0.11–0.63 | 0.08–0.52 |
| 1e-3 | 1.37–2.57 | 0.67–2.46 | 0.78–2.23 |

The departure is largest on broad minima, where many intervals lie within λ of the
minimum (up to about 8 at the default: the rim nearly concentric with the envelope over
half a millimeter). The 144 scans of "Optimization defaults decided" (48 quads near the
minimum, 3 random directions each, ±20 µm at 1e-5 mm steps, then ±1e-6 mm at 1e-9 mm
steps around each scan's largest slope change):

| | largest slope change per 1e-5 mm step (median) | per 1e-9 mm step at it |
|---|---|---|
| exact (`softness = 0`) | 0.21 (4.5e-6) | 0.22: the kink |
| smooth, 1e-4 | 0.0043 (7.0e-6) | 0.0018 |

The kink of the exact clearance is rounded off over about a micrometer of aim; what
remains, a slope change of 0.002 on the same scan, is an interval that holds two local
minima of the distance (its value, the least of the golden result and its ends, is a
minimum of two branches, weighted by its share of the sum). Elsewhere every scan is
smooth to the rounding of the steps. The forward-mode derivative against central
differences at h = 1e-6 mm along 240 random aim directions (5 per quad): largest
difference 3.6e-8, median 3.4e-9. At the former kink the smooth slope is −0.202
(ForwardDiff) against −0.204 (central differences at 1e-4 to 1e-6 mm), the residual
kink; on the other two scans examined there they agree to 1e-6.

The tutorial's comparison of the keep-out with `clearance` ("a few micrometers") is
untouched by 0.17 µm. Cost: about 2.4 times the exact clearance's (28 against 12 µs per
beam, one thread, while other runs shared the machine): more intervals are refined.
The other keep-outs (`SphereKeepOut`, `HardwareKeepOut`, the TaLIS rules) are
unchanged, and so is `keepout_clearance`'s meaning for them. Tests: `test/keepout.jl`
checks the geometry on `softness = 0` (as before, to 1e-9 mm and 1e-12 relative) and the
default against it: within the bound on three holes and three beams, and across the
switch between two edges of a rectangular hole (slopes 0 and −1 on either side, between
them at the switch, ForwardDiff against central differences to 1e-6).

### 4. The tapered cut by default inside `optimize_pointing`

On a scenario whose deposition uses `VisRadProfile` with its hard cut (`taper = 0`,
the default of the profile, unchanged), `optimize_pointing` now runs on the same profile
tapered over one e-fold below the cut (`taper = 1`; "VisRad's beam profile law", "A
tapered cut for optimization") and reports the design under the cut of the scenario.
`smooth_cut = false` optimizes the profile of the scenario as before.

- *Reuse.* The taper keeps the profile's support, so the operator's footprint, its
  points, normals, visibility masks, work plan and envelope ratios are shared by both
  (`_operator_beam` sizes the footprint from the envelope, not the taper; the screens
  `_vr_beyond` test the support): the tapered operator is the scenario's with only its
  profile replaced, and the tapered scenario shares the view factors, the responses
  (`moment_response`, `ylm_response`, the lazy full response) and the albedos. Nothing
  is rebuilt; the per-spot normalization tables of the tapered profile are built once and
  cached, as for any taper. A `PhasedScenario` is tapered phase by phase.
- *Report.* `value`, `initial_value`, `constraints` (with the tapered run's
  multipliers) and `outcome` are the scenario's under its hard cut at the run's design
  and start; `PointingResult.tapered_value` holds the tapered run's objective (`NaN` for
  a run on the profile of the scenario); everything else (parameters, beams, counts,
  `state`, stationarity, bound sensitivities, phase 1, relaxations) is the tapered run's.
- *Exact ceilings.* Decided: report, do not force. The tapered run ends with an exact
  ceiling within `feasibility` (1e-3 of the limit) of its limit; the taper moves wall
  elements by at most 7e-5 of the wall's peak, so under the hard cut the same design is
  within 1.07e-3 of the limit at worst. Rather than continuing the run under the hard cut
  (where the jumps the taper removes return), the rows are computed under the hard cut,
  and a row over its limit by more than `feasibility` shows `:violated` with the outcome
  `:infeasible`: the report states so plainly. On the case below the hard cut's peak came
  out 6e-5 of the limit *below* the tapered one.

Measured (*measured*; 8 threads, another run may have shared the machine): the
dropped-bundle case of `test/constraints.jl` (h 0.45, its coarse operator, bundle 11
dropped, 44 quads in `:phi` and `:z`, the exact wall ceiling at the 48 quads' peak and
the four mode gates at the pointing-noise floor, drive objective, 1,000 iterations,
3,000 evaluations) with the deposition on `VisRadProfile()`:

| run | evaluations (outer) | outcome | drive objective, hard cut (tapered) | peak / limit, hard cut | gates / limit, l = 1…4 |
|---|---|---|---|---|---|
| `smooth_cut = true` (default) | 1,010 (14) | `:interior` | −0.961760 (−0.961772) | 1.00091 (tapered 1.00097) | 0.9994, 1.0005, 0.9998, 0.9996 |
| `smooth_cut = false` | 1,177 (8) | `:box` | −0.961765 | 0.99155 | 0.9642, 0.9991, 0.9940, 0.9970 |

Both end feasible within `feasibility`, at the same drive within 5e-6 (a flat family of
designs); the tapered run took 14% fewer evaluations and ends with every limit active
and no parameter at a bound. The tapered and hard-cut objectives at the tapered run's
design differ by 1.2e-5.

As for existing results, no test, tutorial, example or GUI default optimizes a scenario on
`VisRadProfile` with its hard cut (the tutorials and the GUI use the f-cone model by
default; `examples/n221204.jl` optimizes nothing), so none moves; every other scenario takes
the code path it took, with unchanged results (`tapered_value` `NaN`). The GUI's optimizer
moves when its VisRad beam model is chosen. `test/constraints.jl` checks the default run
against the same run on a scenario built with `taper = 1` (identical in design, beams,
evaluations and multipliers), the value and rows against the hard-cut evaluation, and a
`PhasedScenario`.

### 5. Phase 1 on by default

`optimize_pointing(...; phase1 = true)` is the default (`relaxation` stays off): a
constrained run whose start violates a constraint first minimizes only the
infeasibility and starts its main run from the first design within the limits, or ends with a
verdict ("Constrained re-pointing", "A verdict on conflicting constraints"). A start
that meets every constraint, a restart with `from` and an unconstrained run are
unchanged; `phase1 = false` is the former path.

Measured on every constrained run of `test/constraints.jl` and `test/phased.jl`
(*measured*; each run logged at the step-4 tree, where phase 1 was off, and with the
default on; 8 threads). Runs whose start meets the constraints, runs that set `phase1`
or `relaxation` themselves, and runs stopped at their first evaluation are
identical (22 of 39 runs). The others:

| case | evaluations, off → on | objective, off → on | outcome, off → on |
|---|---|---|---|
| closed form, cap and floor (feasibility 1e-8) | 58 → 56 | 3.335025864e-4 → 3.335025436e-4 | interior |
| the same, restarted at an outer iteration (`from`) | 47 → 36 | as above | interior |
| the same, two outer iterations | 11 → 20 | 3.697e-4 → 4.016e-4 | iterations |
| the same, a 20-evaluation budget, then restarted | 20 → 20, 24 → 37 | 3.3360e-4 → 4.0162e-4, then 3.335025689e-4 → 3.335025436e-4 | evaluations, interior |
| the same, a 10-evaluation budget | 10 → 10 | 7.18e-4 → 1.81e-2 | evaluations |
| a power above the box's cap share (verdict) | 23 → 6 | 1.0007e-2 (same) | infeasible_box |
| bound sensitivities, three runs | 69 → 101, 93 → 113, 80 → 93 | equal within 5e-7 relative | box |
| dropped bundle: ceiling and four gates (drive) | 1,044 → 1,037 | −0.962591 → −0.962284 | interior |
| keep-out, light and power on held power scales (infeasible) | 204 → 1 | −1.04574 → −1.00000 (the start) | infeasible_box → infeasible |
| hard-cut VisRad case (§4), default and taper 1 | 88 → 106 (each) | 8.527433 → 8.527574 | interior |
| the same, `smooth_cut = false` (300-evaluation budget) | 300 → 125 | 8.527446 → 8.527786 | evaluations → interior |
| phased: three phases, gates and a ceiling | 474 → 470 | −1.0000071 → −0.9997388 | interior |
| the same, 10 evaluations | 10 → 10 | −1.0000073 → −0.9997394 | evaluations |

- On the closed forms, the KKT point is reached either way (the designs agree within the
  runs' 1e-8 feasibility); phase 1 costs −23% to +54% of the evaluations. Runs cut short
  by an iteration or evaluation budget end elsewhere, as any change of path does.
- Drive objectives on a flat family of designs end at another member: the dropped-bundle
  case at 3.2e-4 less drive in 7 fewer evaluations, the phased case at 2.7e-4 less, both
  within every limit.
- In the verdicts, the conflicting cap is found in 6 evaluations instead of 23. The keep-out
  case of held power scales starts with a clearance violation (5.5e-3 of its scale) that no
  parameter of that parameterization moves (the aims are held): phase 1 stops at the start,
  a local minimum of the infeasibility, with the verdict `:infeasible` after one evaluation,
  where the augmented Lagrangian spent 204 to end at `:infeasible_box` with the same
  violation.

Tests were re-derived: `test/constraints.jl`'s augmented-Lagrangian-only comparison passes
`phase1 = false` (it is the reference for phase 1's evaluation count), and its checks of
`res.phase1 === nothing` now read phase 1 from the result; every other assertion holds
on both paths. The tutorials' constrained runs are in "6." below.

### 6. The tutorials on the new defaults; the capped polish stays

The user decided that the dropped-bundle tutorial's per-beam polish stays capped at 120
evaluations (a converged one costs minutes of the docs build); its text states what the
cap leaves and that a design to hand off runs to convergence. No code change.

Both tutorials run as scripts at the step-4 tree (phase 1 off) and at the final tree (8
threads; another run may have shared the machine), every assertion holding. The
dropped-bundle tutorial's constrained runs (drive as a fraction of the 48 quads'):

| run | step-4 tree (phase 1 off) | final (phase 1 on) |
|---|---|---|
| within limits | 1,079 evaluations, drive 0.9626, gates 1.000 1.000 1.000 0.995, peak 1.0000 of the ceiling, `:box`; cap review 2 beams failing, worst −60 µm | 978 (phase 1: 15), drive 0.9620, gates 1.000 0.998 0.997 0.954, peak 0.9986, `:interior`; no beam failing, worst +2 µm |
| stage 1 (quads, four limits) | 1,526, drive 0.9574, gates 0.976 1.000 0.998 1.000 | 1,591 (12), drive 0.9550, gates 0.973 0.999 1.000 1.000; cones 1.000 / 0.981 / 1.000 / 0.994 |
| stage 2 (per-beam polish, 120 evaluations) | drive 0.9593; l = 2 0.1% over its gate, peak 0.44% over; review worst 10 µm | drive 0.9567; l = 2 0.19% over, peak 0.28% over, keep-out 1.3 µm short on one beam; review worst 9 µm |
| least re-pointing at 0.5% per degree | 636, 63 µm rms | 508 (12), 63 µm rms, 0.105 mm at most, drive 0.9594 of the 48 quads' |

Phase 1 moves the first run to another member of the flat family at 0.06% less drive,
one that happens to clear the end caps by 2 µm (the text no longer states that the run reaches
into them; it states that the margin is unheld and gives the phase-1-off result). With step
3's smooth keep-out the cap review of the penalty re-pointing gives 8 µm where the text had
9 µm. The diagnostic-windows tutorial has no constrained run; its pair re-pointing
under the aperture keep-out's penalty is 0.105–0.142 mm rms, largest 0.26 mm, where
the text had 0.11–0.15 mm and 0.29 mm, with l = 1–3 still restored to 0.02%; the other
numbers its text quotes are unchanged.

## Evaluation cost, measured (2026-11-19) *measured*

Whether the cost of one objective-plus-gradient evaluation inside `optimize_pointing`
rose across versions, on a public workload, and what made it rise. Every timing is the
best of several runs on one machine (8 fast cores; no other job; the thread count
given), with `BRIAREUS_NO_WORKLOAD=true`, from `git archive` exports of each version
with the same package environment.

**The workload.** The dropped-bundle tutorial's re-pointing: 46 quads (bundle 11
freed), φ and z per quad (92 parameters), `asymmetry_objective` with the tutorial's
goals and an `ApertureKeepOut` penalty (1% contour, minimum 10 µm, weight 1e3).
Case A: the tutorial's default mesh and its coarse f-cone operator (resolution 0.6,
reach 0.5, level 0). Case B: a study's size, the hohlraum at h = 0.2 and VisRad's
profile (resolution 0.6, reach 1.0, level 0). Steps: the deposition operator, the
radiation model, one `evaluate`, the objective's value, its gradient by ForwardDiff
(chunk 12, the chunks serial), the wall time per evaluation of a six-iteration
`optimize_pointing` run (its chunks threaded), and the keep-out's clearance and its
gradient. Versions: the development tree of 2026-10-22 (0faede7), v0.1.0, v0.2.0,
v0.2.1, faa37b8 (2026-11-18), 69761ff (0.3.0) and 0.3.0 with the fixes below.

The same physics in every version: the operator's points agree from v0.1.0 on (328,791
in A, 644,585 in B; 652,649 in B on 2026-10-22, before the VisRad profile law); the
hohlraum's default mesh gained elements at v0.2.0 (4,116 to 4,308 in A, 8,628 to 9,396
in B: a model change, not bisected), and the objective's value agrees to 3e-5 across
versions (VisRad's law changes it by 1.5e-3 in B).

Case B, per evaluation (ms):

| step, threads | 10-22 | 0.1.0 | 0.2.0 | 0.2.1 | faa37b8 | 0.3.0 | fixed |
|---|---|---|---|---|---|---|---|
| `evaluate`, 1 | 4.0 | 22.3 | 22.3 | 22.2 | 22.1 | 22.3 | 3.2 |
| gradient (serial chunks), 1 | 53 | 223 | 229 | 230 | 235 | 264 | 54 |
| optimizer, per evaluation, 1 | 54 | 234 | 244 | 240 | 237 | 273 | 56 |
| optimizer, per evaluation, 8 | 17.7 | 31.7, 36.0 | 34.5, 42.3 | 34.8, 39.8 | 53.5, 39.9 | 50.6, 47.6 | 9.0, 9.2 |
| keep-out gradient (92 parameters), 1 | 2.8 | 2.9 | 11.8 | 11.6 | 11.9 | 30.3 | 17.4 |
| gradient's allocations (MB) | 17.3 | 17.2 | 18.4 | 18.4 | 18.4 | 45.9 | 19.3 |

At 8 threads a short run's time per evaluation varies by about 20% between runs: the
row gives a single run and the best of three where both were taken. Since fix 5 below
(2026-11-20) the keep-out's gradient in case B is 1.8 ms on one thread (v0.1.0: 3.4,
"fixed" above: 18.2, rerun beside it) and its allocations 1.3 MB.

The "fixed" column is fixes 1–9 below, measured together on 2026-11-20 (fix 7 does
not enter case B), except the keep-out's row, fixes 1–4 (fix 5's number is in the
sentence above). The first four rows with only fixes 1–4: 12.3, 162, 168, and 24.6
and 25.6 ms at 8 threads; with fixes 1–4 and 6: 5.9, 112, 117, and 18.7 and 18.5 ms;
with fixes 1–4, 6 and 8 (before fix 5): 3.1, 69, 71, and 14.0 and 12.0 ms.

Case A, per evaluation (ms): `evaluate` 9.0–9.4 on one thread in every version; the
serial-chunk gradient 152, 149, 157, 158, 162, 181 and 165 (fixed); the optimizer per
evaluation on one thread 163, 163, 171, 172, 173, 194 and 178 (fixes 1–9: 151 and
160, at 8 threads 22–23 ms; fix 11, measured beside fixes 1–10 in another session: the
gradient 164 to 141 ms, the optimizer 187 to 156 ms, at 8 threads 21.6 to 19.8 ms); the keep-out's gradient
2.7, 2.7, 10.0, 9.7, 9.7, 28.2 and 14.1 (1.5 since fix 5). Building the operator takes 2.1–2.4 s and the
radiation model 5.7 s (10-22), 7.2 s (0.1.0) and 7.7–8.1 s since (one thread; build
time, not per evaluation; the step at 0.1.0 not bisected: the near-field view-factor
quadrature changed in those days).

**What moved, by commit** (bisected on exports):

- b5173bc (2026-10-27, `VisRadProfile` takes the VisRad profile law in place of the
  tabulated fit): a VisRad-profile deposition costs 5.6 times as much (case B's 46 quads, one
  thread: 3.8 to 22.5 ms on Float64, 5.4 to 27.0 ms on Dual{12}). By design, the law
  takes three non-integer powers per point, and every point of the operator was
  evaluated on Float64 and screened by the full law on duals.
- ed46747 (2026-11-03, the rim of the hole itself, refined): the aperture keep-out 3.6 times
  (0.1 to 0.6 ms for 46 quads; gradient 2.7 to 10 ms). By design, the exact rim and a
  golden-section search instead of 96 samples.
- 29bb958 (2026-11-19, a smooth minimum over the rim): the keep-out a further 2.4–2.9
  times (gradient 10 to 28 ms, allocations 0.4 to 28 MB): every interval within 3 µm
  of the closest sample is refined, not the four beside the closest. `softness = 0`
  costs the same (the window, not the smooth minimum, sets the cost).
- e8e3887 (2026-11-18, bound sensitivities): one gradient at the end of an
  unconstrained run with a parameter at a bound; a constrained run takes its gradients
  from the reported multipliers' Jacobian (one gradient's cost per run). Per run, not
  per evaluation (in the six-iteration runs above, one evaluation in nine).

**The fixes** (2026-11-19, each measured on one thread unless stated):

1. `clearance_floor` over a `PhasedScenario`'s phases computes the keep-out once for
   phases whose beams are alike (same aims and spots; for an `ApertureKeepOut` the
   same envelope radius), not once per phase; results are unchanged.
2. `ApertureKeepOut` on a dual beam runs its samples and golden sections on the
   beam's values and evaluates in duals only the points that enter the result:
   66.5 to 38.8 µs per beam (Float64: 23.4 to 21.2), identical in values and partials
   on 5,760 beam cases (quads and beams in two hohlraums, aims moved up to 0.4 mm,
   some defocused, softness 1e-4 and 0).
3. `deposit` skips beams at exactly zero power (value and partials; a picket measure's other
   cones), with results unchanged (48 quads with the inner cones at zero, f-cone profile:
   9.6 to 4.4 ms on Float64, 19.1 to 8.5 ms on Dual{12}).
4. VisRad's deposition skips the points that a bound without powers puts beyond the
   cut-off (`_vr_far`: α ≥ min(sq, gauss), sq ≥ (max(|u|, |v|)/ΔL)^10 and gauss ≥ q^(n₀/2)
   for q ≥ 1): on 12 cases (48 quads and 192 beams at h = 0.2, hard cut and taper,
   moved and defocused aims) 714 to 370 ms on Float64 and 863 to 495 ms on Dual{12},
   with results unchanged (elements and per-beam totals, values and partials, against 0.3.0 on
   those cases). No skipped point has a nonzero intensity (every point of the 48
   quads' operator checked, and a grid across and along four spots in
   `test/basics.jl`). The screen first sat behind a `@noinline` barrier, which kept the
   depositions unchanged.
5. (2026-11-20) `ApertureKeepOut` finds each sample interval's minimum exactly from the
   distance's slopes at the samples, without a window of refined intervals, and takes a
   dual beam's partials by the chain rule from the clearance's sensitivities to the focus
   and the envelope's semi-axes ("The aperture keep-out's exact interval minima
   (2026-11-20)"). Case B, 46 quads: value 1.49 to 0.20 ms (v0.1.0: 0.11), gradient
   with serial chunks 18.2 to 1.8 ms (v0.1.0: 3.4), the optimizer's threaded gradient at
   8 threads 2.96 to 0.50 ms (v0.1.0: 0.47). Results move by at most 1.2e-7 µm.
6. The far screen inlined (2026-11-19 evening, after the user's ruling that a
   rounding-level change is judged against physics tolerances): the compiler then
   inlines the law's `exp` and `log` in the loop differently, about half the cost
   again. On the 12 cases: 370 to 173 ms on Float64 and 493 to 292 ms on Dual{12} (one
   thread, the two measured in one session). Result change, Float64 depositions only:
   10,905 element entries differ, by at most 1.24e-13 of the value of that element (3.7e-16
   of the largest element), per-beam totals by at most 3.3e-16; the Dual{12}
   depositions, values and partials, are unchanged. Against the tolerances of
   delivered power (2–5%) and the drive's modes (tenths of a percent), 1e-11 of the
   smallest. Case B's objective prints the same 15 digits and its six-iteration run
   takes the same 8 evaluations to the same value; the dropped-bundle tutorial, run as
   a script, prints the same output.
7. The picket ceilings take the deposition of the evaluation (2026-11-19 evening). A
   ceiling at the powers of the beams (no `pulses` or `power`: a phase's evaluation at
   that phase's powers) deposited the outer cones again on every evaluation. Inside a
   constrained `optimize_pointing` run the scenario now keeps, in the deposition
   pass itself, the element powers of each group set such a ceiling reads (`_GroupSplit`, a
   `Scenario` field set per phase from the constraints), and the ceiling reads them;
   outside a run, with other powers, or when the laser–plasma model deposits light
   itself (or, with `lpi = false`, changes the beams) it deposits as before. The kept
   sums are those of the deposition by the ceiling itself, in the same order
   (`test/picket_fluence.jl`: every measure, overlapping group sets, with and without
   backscatter, to 1e-12). On the multi-phase study (scratch runs: six phases, 552
   parameters, the three picket ceilings on `:average`, plain and tapered cut) the
   Dual{12} depositions and every term's violations, values and partials, are
   unchanged; on Float64 the loop with the kept sums compiles the law differently and the
   depositions move by at most 5.9e-16 of the value of an element, the terms' violations
   by 6.8e-14 of their scale. A constrained run of 20 evaluations takes 53.4 s without
   and 40.5 s with (8 threads), the same 20 evaluations, and ends 6.6e-8 (mm, or scale
   units) away, its objective 4.7e-7 relative (the round-off moves the line searches).
   The optimizer's augmented-Lagrangian gradient falls from 2.39 to 1.87 s at 8 threads
   (1.57 s for the objective only): the constraints from 0.52 to 0.19 of the objective.
8. Cheaper screens and coordinates in VisRad's deposition (2026-11-20). (a) `_vr_far`
   adds to its test of the smaller term a bound on the weighted sum itself, α ≥ f r^10 +
   (1 − f) q^k (q^(1/4))^j (r ≥ 1, q ≥ 1, n₀/2 ≥ k + j/4), with multiplications and two
   square roots: on case B's operator the old bound skipped 84.5% of the points and left
   60,600 beyond the cut beside 39,182 lit ones; the new one skips 92.4% of those 60,600
   (99.2% of all points beyond the cut). (b) A run of 8 consecutive points (in the
   operator's order, by destination) is skipped when its bounding ball lies beyond the
   cut by the same bound over the ball (78.5% of case B's points, one test per run).
   (c) A dual beam no longer screens the remaining points on values first: 90% of them
   are lit, and the law gives an exact zero, value and partials, beyond its cut. (d) A
   point's coordinates in the beam frame are taken once, on values, and a dual beam's
   carry the partials every point shares (those of the origin's coordinates: the frame's
   axes carry none). Screens (a)–(c) skip only exact zeros. Results: on the 12 cases of
   fix 4, Float64 depositions move by at most 1.6e-13 of the value of an element
   (3.8e-16 of the largest), Dual{12} values and partials by at most 1.5e-15 of the
   largest element's, per-beam totals by 2.2e-16 and their partials by 1.2e-15 of the
   beam's power; the f-cone model is unchanged. Case B's objective prints the same 15
   digits and its six-iteration run takes the same 8 evaluations to the same value; the
   dropped-bundle tutorial prints the same output and `examples/n221204.jl` at 0.8/0.4
   writes a byte-identical CSV (both on the f-cone model). Timings (one thread unless
   stated; another run shared the machine's second lock slot): the 12 cases'
   depositions 171 to 85 ms on Float64 and 297 to 155 ms on Dual{12} (in one session);
   case B's `evaluate` 5.9 to 3.1 ms, its serial-chunk gradient 112 to 69 ms, the
   optimizer per evaluation 117 to 71 ms and at 8 threads 18.5–18.7 to 12.0–14.0 ms; the
   multi-phase study's gradients at 8 threads 1.57 to 0.88 s (objective) and 1.87 to
   1.32 s (the optimizer's augmented Lagrangian), on one thread 6.45 and 7.44 s (10.8
   and 15.7 s with fixes 1–4).

   The same coordinates on the f-cone model, measured here and held back while the
   dropped-bundle tutorial quoted the digits of its long runs, landed as fix 11.
   Building VisRad's frames once per beam instead of per work unit showed no difference
   beyond the machine's noise (about 10%) and allocated 8% more per gradient: not
   landed.
9. An evaluation computes its moment and Y_lm products once (2026-11-20): `flux_moments` (on
   the cached degrees, axis and center), `mean_flux` through it, and `flux_ylm` read the
   scenario's response product with the deposition from the evaluation after the first read,
   instead of each forming it again (a dual matrix product, about 0.4 ms per read per dual
   pass on the multi-phase study, where a phase's gate, drive floor, drive ceiling and the
   objective read the same products: 16 products a pass before, at most two a phase after).
   Results are unchanged (the same product, read as a copy). On the study the optimizer's
   augmented-Lagrangian gradient allocates 2.14 GB instead of 2.34 GB; its time at 8 threads
   moved within the machine's noise that session (1.17 to 1.69 s before, 1.51 to 1.62 s
   after, alternating, another job on the second slot).
10. Phases on one mesh deposit in one pass (2026-11-20). A `PhasedScenario` keeps, per
   phase, the first phase whose deposition operator has the same points, weights,
   normals, visibility, destinations, ports and profile (the same operator, or one built
   alike on a target that differs in its materials; compared once at construction), and
   `evaluate` deposits together the phases of one such class whose effective beams
   differ only in power (an exact test on aims, focus, spots and ports, values and
   partials) and whose laser–plasma models deposit no light themselves: each point is
   evaluated once, for the beam at unit power, and added to every phase times that
   phase's beam power, into its elements, back sides and kept group sets (fix 7). A
   beam's deposition is linear in its power (the profile's shape, cut-off and screens do
   not read it); a model that changed a spot by power would fail the exact test and
   deposit alone, as phases on other points do. Results change at rounding: on the
   multi-phase study (scratch runs: six phases on two meshes, three each; Float64 and two
   Dual{12} chunks, plain and tapered with the kept group sets) the depositions move by
   at most 2.2e-13 of the value of an element (6.4e-16 of the largest), per-beam totals by
   8.5e-16, back sides by 3.1e-16, the kept group sets by 4.6e-16, every term's
   violations by at most 3.2e-13 of its scale and the objective by 4.5e-16. A constrained
   run of 20 evaluations takes the same 20 evaluations and ends 5.9e-8 mm away (1.6e-8
   mm rms in the aims, 5.9e-8 in the scales), its objective 4.2e-7 relative. Timings,
   the former path being the same build with the grouping switched off, alternated in one
   session (another job held the second lock slot): at 8 threads a dual pass's
   evaluation 37–40 to 15–17 ms, the objective's gradient 1.21–1.26 to 0.46–0.49 s and the
   optimizer's augmented-Lagrangian gradient 1.58–1.59 to 0.74 s (allocations 2.14 to
   1.87 GB), the 20-evaluation run 27.0 to 13.3 s; on one thread the dual pass's
   evaluation 138 to 52 ms (Float64: 81 to 27 ms), the objective's gradient 7.11 to
   2.65 s and the augmented Lagrangian's 7.85 to 3.29 s.
11. The f-cone model's point coordinates taken once (2026-11-20), as fix 8 (d) does for
   VisRad's: a point's (s, u₀, v₀) in the beam frame are computed on values, are used by the
   screens, and a dual beam's are built from them with the partials of the origin's
   coordinates, which every point shares (`_cone_coords`). Results change at rounding:
   on the f-cone model's six cases of fix 4's set (48 quads and 192 beams at h = 0.2,
   moved and defocused aims) Float64 depositions move by at most 4.8e-14 of the value
   of an element (3.4e-16 of the largest), Dual{12} values and partials by 1.3e-15 of the
   largest element's, per-beam totals by 2.2e-16 and their partials by 1.6e-15 of the
   beam's power. Case A's objective prints the same 15 digits and its six-iteration run
   takes the same 8 evaluations to the same value (12 digits). Timings (one thread unless
   stated; another job held the second lock slot): the six cases' Dual{12} depositions
   779 to 721 ms (Float64 386 and 384 ms), case A's serial-chunk gradient 164 to 141 ms,
   the optimizer per evaluation 187 to 156 ms and at 8 threads 21.6 to 19.8 ms. The
   dropped-bundle tutorial's long constrained runs end elsewhere in their families, as
   any rounding moves them, within the ranges its text gives ("The dropped-bundle
   tutorial's long runs under rounding (2026-11-20)", where both trees' runs are among
   the 13 starts); every assertion holds. The beam-errors tutorial's tuning (200
   iterations) ends nearer its optimum (P₂/P₀ +0.000% against −0.023%, the l = 1
   budget 0.236% against 0.233%): its text now states 0.23–0.24%. The
   diagnostic-windows tutorial's pair re-pointing ends elsewhere on its flat family, 134
   and 98 µm rms in z and along the wall, the largest move 0.267 mm (before 141, 93 and
   0.247 mm), l = 1–3 still restored to 0.02%: its text now states 0.25–0.27 mm; its
   power trim converges in 106 iterations instead of 111, its printed cone factors moving
   by at most 0.01 (none quoted in the text), and the pulse table's entries by at most
   0.01 (within the text's ranges). The backlighter-sweep tutorial prints the same output.

Together, on case B the optimizer's evaluation costs 0.21 of 0.3.0's on one thread (56
against 273 ms) and 0.18–0.19 at 8 threads (9.0–9.2 against 47.6–50.6 ms), 0.24–0.29 of
v0.1.0's, and about the 2026-10-22 tree's on one thread (1.04; 0.5 at 8 threads), all nine
fixes measured together on 2026-11-20 (case A: 160 ms on one thread, 22–23 ms at 8). Fixes
1–4 by themselves (0.52–0.62 of 0.3.0's) left every result identical to that of 0.3.0; fixes
6–8, 10 and 11 change depositions at the rounding level only, and fix 9 nothing. The
keep-out's gradient, 6 times v0.1.0's after fixes 1–4 (the refined rim and the window), is
0.53 of it after fix 5. On a multi-phase study's evaluation (scratch runs: six phases, 552
parameters, 17 constraints, 8 threads) the augmented Lagrangian's gradient fell from 7.9 to
3.6 s (0.45) with fixes 1–4, and the constraints from 0.98 to 0.40 of the objective's
gradient: before the fixes the keep-out over the six phases and the three picket ceilings'
depositions each cost about a quarter of the gradient. Fix 6 then takes the objective's
gradient from 2.83 to 1.62 s and the augmented Lagrangian's from 3.80 to 2.44 s (the two
measured in one session; another run may have shared the machine, whose second lock slot is
open to other jobs), fix 7 takes the augmented Lagrangian's to 1.87 s, and fix 8 takes the
two to 0.88 and 1.32 s; with all nine, on one thread, 7.1 and 7.7 s. The six phases'
depositions were then more than half of a dual pass (56% before fix 9), on two meshes (three
phases each, the same operator points), the law evaluated six times where two would do: fix
10 evaluates it twice, and takes the two gradients to 0.37 and 0.42 of fix 9's on one thread
(2.65 and 3.29 s) and to 0.38 and 0.47 at 8 threads (0.46–0.49 and 0.74 s; fix 10 does not
enter cases A and B, one scenario each). The gradient's chunks run threaded on the phased
path (0.3.0: 6.3 times faster at 8 threads than on one, 7.2 at 16, fast and slow cores).

**The end state, on a quiet machine** (2026-11-20, on the tree at the end of that day's
work, whose commits after fix 11 do not enter an evaluation; the timings of fixes 6–11 above
were taken while another job held the machine's second slot). The multi-phase study's
gradients (552 parameters, the same scratch harness, best of several): on one thread the
objective's 30.3 s (0.3.0) to 2.29 s and the augmented Lagrangian's 49.7 to 4.98 s; at 8
threads 0.58 s and 0.78 s, the augmented Lagrangian's against 7.9 s on 0.3.0. The harness
builds the augmented Lagrangian outside `optimize_pointing`, so its three picket ceilings
still deposit (15 ms each per dual pass on one thread, the whole of the constraints' cost
there beside the keep-out's 4.5 ms): inside a run fix 7 removes them, and the run's gradient
is nearer the objective's. One dual pass of `evaluate` is 47 ms on one thread and 11 ms at
8. A study script's penalty re-pointing with VisRad's profile (scratch runs, 30 iterations
at 8 threads, v0.1.0 beside the end state in one session): 0.123 to 0.036 s per evaluation
at 92 parameters and 0.478 to 0.145 s at 368 (0.148 and 0.711 for 0.3.0, measured the day
before).

**By design, not changed.** Phase 1 (on by default in constrained runs since 2026-11-19)
adds separate evaluations before the main run, not cost per evaluation; each lowering of
an exact ceiling's working limit continues the run; the reported multipliers cost one
Jacobian at the end of a constrained run (as much as one gradient, whatever the number
of constraints); the optimizer's tapered VisRad cut (`smooth_cut`) costs as much as the
hard cut per deposition (within 15%). A picket ceiling outside `optimize_pointing`, or
one given the picket's powers by `pulses` or `power`, deposits the outer cones once per
evaluation (fix 7 supplies the powers of the beams inside a run): a hand-written sum of
`ev.deposition.per_beam` over the outer beams costs nothing but measures another
quantity (it counts the light on the capsule and the entrance holes' outside).

## The aperture keep-out's exact interval minima (2026-11-20) *measured*

Why `ApertureKeepOut`'s gradient cost 6 times v0.1.0's, and the design that replaces the
window of refined intervals ("Evaluation cost, measured (2026-11-19)", fix 5).

**Before.** Every sample interval whose nearer end lay within 3 µm of the closest sample
was refined by a 40-step golden section ("Defaults changed (2026-11-19)" §3): on NIF
quads through a hohlraum's holes 10 such intervals per beam (median; up to 36 on a broad
minimum), about 1,500 distance evaluations per beam, repeated in every chunk of a
gradient, and the terms of the smooth minimum evaluated in the gradient's duals.

**Now** (the smooth minimum over sample intervals, the softness and the rim unchanged):

- *Each interval's minimum, exactly.* The rim points and their derivatives along the rim
  are stored at the samples, and the search computes the distance and its slope there
  (on Float64, a dual beam's values). An interval whose slope rises through zero holds
  a minimum, found by regula falsi with the Illinois step to 1e-12 of the interval (6
  steps on average); otherwise its minimum is its lower end. An interval is first split
  where the f-cone term |s| has its kink (the rim crossing the plane of best focus:
  closed form on an edge and an arc; the kink is concave, never a minimum); where the
  rim crosses an axis of an elongated envelope, once the direction in envelope units
  turns over the interval through a quarter of 1/√(e² − 1) (e the ratio of the
  semi-axes: the angle over which the Sampson distance's direction factor peaks); where the
  chord through the ends passes nearest the beam axis, if it reaches inside the
  envelope (there the Sampson distance (q − 1)·h falls to −h within the axis's
  distance); and where the cubic Hermite interpolant of the ends' distances and slopes
  turns twice (a minimum and a maximum between two samples). Each split point is a
  candidate.
- *Weights on the least interval minimum.* σₖ is 1 within 10 λ of the least interval
  minimum and falls (C¹) to 0 at 20 λ (a term at 10 λ weighs e⁻¹⁰ = 4.5e-5 of the least;
  before: within 3 µm of the closest sample, 0 at 6 µm). Only the terms within 20 λ enter
  in the beam's number type. *Engineering* constants.
- *Partials by the chain rule.* On a dual beam whose values are Float64 the terms are
  evaluated in 5-partial duals of the keep-out, seeded on the focus (three) and the
  envelope's semi-axes, the only quantities of the beam the distance reads, and the result's
  partials are those sensitivities times the beam frame's partials: the cost no longer grows
  with the gradient's chunk. Nested duals evaluate the terms in the number type of the beam.
  The clearance of a dual beam has exactly the value of the clearance of the Float64 beam.
- The envelope radius (a root search) is computed once per distinct spot and power in
  a call.

**Validation** (*measured*, softness 0, against an independent reference: the distance
sampled at 401 points per sample interval along the analytic rim pieces, the least
refined by golden section between its neighbors). 2,820 beams (NIF quads and beams through
two hohlraums' holes, aims moved up to ±0.4 mm, defocused up to ±2 mm, elliptical spots
turned at random, and 900 oblique beams through the three holes of `test/keepout.jl`):
within 9e-16 mm. 18,000 random oblique, defocused beams with spots of semi-axes 0.05 to
0.3 mm (aspect up to 6) through those three holes, 5,058 with positive clearance and
12,942 clipping: within 8e-15 mm; the former golden section missed the minimum on 37
clipping beams by up to 0.29 mm (the dip near the axis inside the envelope). An
intermediate form without the axis splits and the Hermite test missed it on 3 of the
positive-clearance beams by up to 0.2 µm; `test/keepout.jl` keeps those beams and two
whose axes pass 10 µm inside and 2 µm outside a rim, against 10⁵ samples per piece.

**Change of results** (the 5,760 beam cases of "Evaluation cost, measured (2026-11-19)"
fix 2: quads and beams in two hohlraums, aims moved up to 0.4 mm, some defocused):

| | softness 1e-4 (default) | softness 0 |
|---|---|---|
| clearance, new − former (µm) | −2e-13 to +1.2e-7 (median 1e-8) | −4e-13 to +8e-13 |
| partials, largest relative change | 5.6e-7 (median 1.1e-8) | 1.1e-8 (median), except at exact ties |

The partials move by the former golden section's position error (about 1e-9 mm along
the rim). At softness 0 the partials differ on 31 of the 2,880 beams, all at the unmoved
aims, whose symmetric pointing ties two rim points: the exact minimum's derivative takes
the other one-sided slope there.
Differentiability, on the 144 scans of §3 (48 quads near the minimum, ±20 µm of aim):
largest slope change per 1e-5 mm step 4.3e-3 (as before), per 1e-9 mm step at it 1.76e-3
(1.77e-3: the interval that holds two local minima); forward mode against central
differences at h = 1e-6 mm along 240 directions: largest difference 4.9e-9 (3.6e-8
before), median 1.1e-10 (3.4e-9).

**Cost** (*measured*, best of 50, one machine, the dropped-bundle workload of "Evaluation
cost": 46 quads, 92 parameters; v0.1.0, the former tree 3a89a32 and this one):

| case B (h = 0.2) | v0.1.0 | former | now |
|---|---|---|---|
| value, 1 thread (ms) | 0.114 | 1.49 | 0.197 |
| gradient, serial chunks of 12, 1 thread (ms) | 3.41 | 18.2 | 1.81 |
| the optimizer's threaded gradient, 8 threads (ms) | 0.469 | 2.96 | 0.498 |
| allocations of a gradient (MB) | 0.41 | 13.1 | 1.26 |

Case A (the tutorial's mesh): value 0.114, 1.26 and 0.164 ms; gradient 3.40, 15.5 and
1.51 ms; threaded at 8 threads 0.38, 1.92 and 0.50 ms. Float64 per beam on the 5,760 cases
18.9 to 4.2 µs, Dual{12} 30.7 to 4.5 µs. The value costs 1.7 times v0.1.0's (96 samples
and no refinement then; 231 samples and exact minima now); what remains of a dual pass
is that same search on the values, repeated in every chunk of a gradient (about 4 µs a
beam a chunk): a cache across the chunks would cut the gradient to about a fifth, at
the price of mutable state in the keep-out, not built.

**Constrained runs** (*measured*, 8 threads; another run may have shared the machine).
A study-sized re-pointing in scratch (a dropped bundle, the kept quads shifted in z and φ, the study's
hohlraum at h = 0.2 with VisRad's profile, the mode-gate and peak-irradiance hinges, the
entrance-hole keep-out as a `clearance_floor` at the least clearance of the dropped
design, active on about 50 beams; 3,000-evaluation budget, which both runs reach):
the designs differ by 0.35 µm rms in aim (1.0 µm at most), the objective by 3.4e-5
relative, the least clearance by 3e-4 µm; 291 s before, 264 s now.
The dropped-bundle tutorial's re-pointing under the keep-out's penalty (`second_try`, 20
iterations): 87 evaluations either way, aims within 0.002 µm (0.0009 µm rms), the
objective within 1.7e-7 relative, the least clearance within 2e-5 µm.
The whole dropped-bundle tutorial prints the former tree's output to every digit, and
the beam-errors tutorial's printed numbers agree with its text. The diagnostic-windows
tutorial's pair re-pointing under the keep-out's penalty (an unconstrained run of up to
200 iterations on a flat family of designs) ends elsewhere in that family, at 0.093–0.141 mm rms and 0.247 mm at most where the former
tree ends at 0.105–0.142 mm and 0.257 mm, with l = 1–3 within 0.02% either way; its text
now states 0.09–0.14 mm and 0.25 mm. The repository's constrained runs with an aperture
keep-out floor (`test/phased.jl`) hold it far below the designs' clearance; their tests
pass, plain and with bounds checks.

## The dropped-bundle tutorial's long runs under rounding (2026-11-20) *measured*

Why the documentation build failed on CI's Linux runner at the dropped-bundle tutorial's
per-beam polish (`@assert` on every limit's normalized violation at most 0.01) while the
same commit printed the expected output on a Windows machine, and what the tutorial's
long constrained runs determine.

**Method** (scratch runs, 8 threads; another job held the machine's second lock slot).
The tutorial's runs repeated from 13 starts: the start of the tutorial on two trees (the tree
before the f-cone model's shared coordinates and the tree with them, whose depositions
differ at 1e-15), and starts moved by 10⁻¹² or 10⁻⁹ mm in a random direction (seeds 1–5,
on both trees). On one tree and one start the runs are identical at 4 and 8 threads
(every printed digit and evaluation count): the threading is deterministic,
and a CPU that rounds differently (CI's) moves a run as a start moved by 10⁻¹² mm does.

**What is determined** (ranges over the 13 starts):

| run | determined | not determined |
|---|---|---|
| re-pointing within limits (drive, gates, ceiling; ~1,000 evaluations) | outcome converged (`:interior` or `:box`), no limit violated; drive 0.9616–0.9625 of the 48-quad drive; peak 0.992–1.0006 of the ceiling | which gates are active (departure/gate l = 1 0.78–1.000, l = 2 0.994–1.000, l = 3 0.91–1.000, l = 4 0.66–1.000); evaluations 843–1,500; the caps' worst-case clearance, −38 to +4 µm, 0–2 beams failing |
| the same limits as penalty terms | drive 0.959208, departures/gate 0.002 0.120 0.005 0.066 at the printed digits, peak 1.0089–1.0098 | evaluations 1,072–1,142; caps (not printed) |
| stage 1 (quads, four limits) | outcome `:box`, no limit violated, no beam failing the caps' check | drive 0.950–0.957; cones 0.990–1.000 (23.5°), 0.967–0.985 (30°), 1.000 (44.5°), 0.990–1.000 (50°); gates active 1–3; keep-out least 13–66 µm; evaluations 656–1,653 |
| stage 2 (per-beam polish, 120 evaluations) | outcome `:evaluations`; gates within 0.2% and the keep-out within 1.3 µm of their limits (largest normalized violation of the gates, keep-out and light ceiling 0.0019); the ceiling's 64-norm 0.944–0.994 of its limit; drive 0.15–0.36% above the first stage's aims on beams; no beam failing the caps' check, worst 8–49 µm | the wall's peak, 0.966–1.020 of the ceiling (over it by 1.2–2.0% on 3 of 13 starts); drive 0.952–0.959 |
| least re-pointing at 0.5% per degree | outcome `:interior`; 62.6–62.8 µm rms, 0.1046–0.1047 mm at most, drive 0.95936; departures at the printed digits | evaluations 508–768 |

**The failure.** The budget-stopped polish holds the irradiance ceiling's smooth
stand-in, the area-weighted 64-norm, not the peak over the elements: the exact ceiling
compares the peak with the limit and lowers the stand-in's working limit only once a run
has converged on the stand-in (`irradiance_ceiling(...; exact = true)`). Stopped after
120 evaluations, the polish's peak ends anywhere in a 5% band about the ceiling, and on
3 of the 13 starts above the 1% the assert allowed. The other limits the polish works
on are held well inside 0.01.

**The tutorial now** asserts for the polish what every start meets with margin: the
outcome, the gates', keep-out's and light ceiling's normalized violations at most 0.01
(largest 0.0019), the stand-in at most 1% over its limit (largest 0.994 of it), and the
clearance check passing; the peak over the elements is printed, not asserted, and the
text gives its range. The within-limits run's and the first stage's text quotes ranges
for what the family does not determine (drive, active gates, the caps' clearance)
instead of the digits of the run. The penalty run and the least re-pointing are
determined, and their text keeps its numbers. Making the polish's peak determined would
take a converged polish (861 evaluations on an earlier first stage, minutes of the docs
build), left as the user decided ("Defaults changed (2026-11-19)" §6). *Superseded in part* ("Least violation, capped inner runs and phase 1's budget", § 8):
a budget-stopped run now returns its least-violating iterate, and the polish's peak
ends within the ceiling on every start measured.

## A ladder study's hand-written pieces as API (2026-11-20) *derived*

Roadmap item 34. The newest script of a multi-phase re-pointing study (a design ladder:
six phases, 184 beams in 46 quads, one aim shift per quad and one power scale per cone
and early phase, 100 parameters, 212 constraint rows, 184 of them the entrance-hole
keep-out's) still wrote pieces by hand that the earlier pass (§ "The multi-phase
study's hand-written pieces as API") had not covered. Those below are public now. A
scratch copy of the script with those lines rewritten on the API, run beside a copy of
the original on the same tree (8 threads, the start only, with the start's objective,
every constraint row and the gates' limits printed to 17 digits), prints the same log
as the original, identical in every character: the start report, the objective, the 212
rows' values, limits and statuses and every phase's gate limits. Each subsection names
the lines it replaced and what of the log checks them.

### 1. The pointing-noise floor and per-degree gate multiples

The script wrote `[sqrt(2l + 1) * r for (l, r) in enumerate(u.rms_mode)]` from a
`mode_uncertainty` result where it needed the floor (the gates, the objective's
normalization, the report), the expression of the internal `_gate_values`, and its
gates as a multiple of the floor per degree (`GATES .* floor`, e.g. 1, 2, 2, 2).
`mode_uncertainty` now returns `floor` beside `rms_mode`, `√(2l+1) × rms_mode[l]`, the
rms of `mode_departure` under the errors (the long Monte Carlo check of
`test/uncertainty.jl` compares the sampled departures' rms with it, within 12% at 400
samples), and `mode_gate(...; fraction)` takes one multiple per degree as well as one
number, for a floor given as the gates or as `floor`. Chosen over a function
(`mode_floor(u)`): the floor is a quantity of the result, read where the result is, and
a field does not add a name. `shift` already takes several modes, and `floor` a vector.
Results are unchanged: the floor is the same expression as before, and a vector of multiples
the same product as the script's. In the script the floor of each reference and its
gates (`mode_gate(reference, GATES .* floor)`) became `u.floor` and
`mode_gate(reference, u; fraction = GATES)`: the 16 gate rows and the 24 gate limits
agree to 17 digits.

### 2. The beam-weighted regularization of the power scales

The script's objective held the early phases' power scales near a uniform start, `Σ_j n_j
(s_j/s0_j − 1)² / Σ_j n_j`, each scale weighted by the number of beams it governs, which it
counted by reading the internal field `ps.param.by`. `power_scale_objective(ps; reference,
weight)` is that term: the scales of `power_scales` (so generic in the number type and
differentiable), `reference` one scale, a vector in the order of the keys or a function
`(key, phase) -> s0`, and per scale the summed `weight` of its beams (default 1 each: the
mean over beams and scaled phases). It was chosen over a public accessor of the counts: it
is placed beside `shift_objective` (the mean squared aim shift, also a mean over beams),
composes as a term of an objective function of the `PerPhase` evaluation as every other term
does, and keeps the grouping internal; the script uses the counts for nothing else. A plain
`Scenario` with `PowerScales` has no phase powers to read its scales against from an
evaluation (its scales are parameters), so the term is for phased scenarios only.
`test/phased_study.jl` checks it against the hand count (16 inner and 32 outer quads per
scale of a cone), a per-beam weight and the closed-form gradient `2 W_j (s_j/s0_j − 1)/s0_j
/ Σ W`. In the script the objective's three lines of counts and reference scales became
`power_scale_objective(ps; reference = (key, phase) -> S0[phase])`: the start's objective
agrees to 17 digits.

### 3. A phased scenario at another deposition resolution

The script's check evaluated a finished design on a finer deposition operator by
building every phase's target, reference and scenario again. `remake(::Scenario;
deposition_kw)` already rebuilt one scenario's operator; `remake(::PhasedScenario;
deposition_kw)` now rebuilds every phase's for the pointing's beams and keeps each
phase's target, radiation model, response, laser–plasma model and time. Phases that
shared one operator share the rebuilt one, and a phase on a target of an earlier
phase's mesh (the same triangles, sides, components, opacity and elements; materials
may differ) takes that phase's rebuilt operator on the target of that phase (an operator
reads the geometry and the ports, not the materials), so one operator is built per mesh
instead of per phase and the phases of one mesh keep depositing in one pass (the
`_mesh` classes of "Evaluation cost, measured (2026-11-19)", fix 10). On the mock-up of
`test/phased_study.jl` (three phases on targets that differ in albedo) the remade
scenario's depositions and modes equal those of the phased scenario built afresh at the
new resolution to 1e-12, its three phases in one class and one pass. The `deposition_kw`
of each phase is replaced too: one set of options for every phase. In the script the
check of a finished design on the finer operator (the study's default resolution, for
the designed windows) built the six references, the six phases and the ceiling's
scenario again; rewritten, it remakes them (`remake(ps; deposition_kw)` for the phases,
`remake(scenario; deposition_kw)` for the references and the ceiling), with two
operators built for the six phases instead of six. Run on a finished design of the
ladder (100 parameters, quads moved by tenths of a millimeter), the rewritten check prints the same
log as the original, identical in every character, the drive ratios and departures of every
phase printed to 17 digits included; the leg took 3.6 instead of 5.5 minutes, setup
included (8 threads; another run may have shared the machine).

### 4. The departure per degree split by order

The script reported each degree's departure split into its `m = 0` and `m ≠ 0` parts
(`√((2l+1) Σ_{m ∈ M} (c_lm/c_00 − r_lm/r_00)²)`), from `flux_ylm`'s exact facet integrals.
`mode_departure(...; m = :zero | :nonzero)` gives the parts, the split and the names of
`mode_rms`'s `m` (`:all`, the default, is exactly the departure as before), and
`PointingSummary` records them as `departure_zero` and `departure_nonzero` whenever it has a
`drive_reference`. The zonal part depends on the coefficients' polar axis (the hohlraum axis
for `flux_ylm`). `test/pointing_summary.jl` checks the zonal part against
`√(2l+1) |Δ(c_l0/c_00)|` by hand, the parts' squares against the whole, and two single-mode
changes against their closed forms. In the script the report's hand-written sum over orders
became `mode_departure(e, r; m = :zero)` and `m = :nonzero`: its four split columns agree at
their printed digits in every phase, at the start and at the finished design of § 3's check
(where they are not zero).

### 5. Quads moved past given shifts

The script counted the quads whose largest aim shift from the as-shot beams exceeded
0.1, 0.2 and 0.3 mm. `pointing_summary(...; moved = (0.1, 0.2, 0.3), moved_by = :quad)`
records the counts (`moved`, of `ngroups` groups) and shows them in the one-line
display; `moved_by` is a beam field (`:group` for cones, `:name` for beams) or a
function of a beam. `test/pointing_summary.jl` checks the counts by hand on rigid
shifts per quad and per cone. In the script the report's six lines of shifts, per-quad
maxima and clearance became one `pointing_summary(ps[:picket], beams[:picket];
reference, moved, keepout)`: the aims line agrees at its printed digits, at the start
(no quad moved) and at the finished design of § 3's check (most quads moved past
0.1 mm, some past 0.3 mm).

**Not built.** The script's uniform start scales (each early phase's scale set so that
the start's drive is a fraction of the reference's,
`fraction × mean_flux(reference) / mean_flux(start)`) are exact only where the drive is
linear in the beams' powers; with CBET or another power-dependent laser–plasma model a
single ratio is a first step of an iteration. Left as a few lines of a study (no design
obvious enough for a public name).

## Least violation, capped inner runs and phase 1's budget (2026-11-20) *measured*

A variant of a multi-phase re-pointing study (scratch runs) whose limits admit no design
exposed three defects of `optimize_pointing(...; constraints)`. Its start met the
entrance-hole keep-out (a `clearance_floor`) and the peak-irradiance ceiling and lay far
over the picket's mode gates; the run, phase 1 off and every inner run stopped by the
caller's `iterations`, ended "infeasible at the largest penalty" and returned a design
inside the keep-out on many beams and over the ceiling: worse than its start on the two
limits the start met.

### 1. The iterate of least violation

The returned design was the last outer iterate (the best augmented-Lagrangian point of
the last inner run, at the largest penalty); nothing else was kept. Now every iterate
(the start, phase 1's design, the design after each inner run, a stopped inner run's
best design) is recorded with its largest normalized violation over the reported rows,
`max(0, maxⱼ (±(value − limit)/scale))`, the measure of `ConstraintRow.violation` and of
`feasibility` (an exact ceiling's exact value against its limit). `least_violation` is
the least of them, violations within `feasibility` counting as equal and the lower
objective then preferred; `last_feasible` is the latest within `feasibility`. A run that
ends above its limits (at the penalty's ceiling, `outer`, the evaluation budget, the
time limit or the callback) returns `least_violation`; `state.x` keeps the last iterate,
whose multipliers `state` holds, and a restart with `from` continues there. Chosen over
ψ = ½ Σ max(0, gᵢ)² (phase 1's measure) because the largest violation is what the rows,
the outcome and `feasibility` judge; a run that met the limits at some iterate returns a
design within them. Cost: one constraint evaluation per outer iteration (no deposition).

Checks (`test/constraints.jl`): the conflicting closed form without phase 1 (a 120 TW
cap against a drive floor 5% above the cap's largest drive) returns an iterate whose
largest violation is at most the last iterate's (0.0252364 either way; outer iteration
10 of 11, the objective 1.3e-8 relative from the last's); a budget stop of the quadratic
program started on its cap (8 evaluations) returns a design within the cap where the
last iterate exceeds it. Every other constrained run of `test/constraints.jl`,
`test/phased.jl`, `test/picket_fluence.jl` and `test/pointing_summary.jl` (54 runs
logged at c40cb98 and after) is unchanged in design and evaluation count: they end
within their limits, or at a verdict whose design is the least violating.

### 2. Inner runs ended by their cap

ALGENCAN raises the penalty tenfold when V = maxᵢ |max(gᵢ, −λᵢ/μ)| has not halved
(Andreani, Birgin, Martínez and Schuverdt, SIAM J. Optim. 18, 1286, 2008; Birgin and
Martínez 2014, Algorithm 4.1), a test of the outer method that presumes the subproblem
solved to its tolerance (step 1 there; the ω_k of Conn, Gould and Toint, SIAM J. Numer.
Anal. 28, 545, 1991). An inner run stopped by `iterations` or `inner_evaluations` has
not solved it, and on the study every such run raised the penalty, by 10⁸ in nine
capped runs, then the verdict. Neither source states a rule for a capped subproblem
(the theory assumes step 1 succeeds); four were measured on the closed forms of
`test/constraints.jl` (8 threads; distance relative to the closed-form KKT point):

- *hold*: after a capped run, the multipliers updated (λ ← max(0, λ + μ g)), the
  penalty, the inner tolerance and the test's reference V held;
- *hold μ*: the same, the tolerance tightened and the reference updated;
- *separate*: the penalty test as before, but the infeasible verdict only after a
  completed inner run (the penalty held at its ceiling meanwhile);
- *hold above the limits* (chosen): *hold* after a capped run that ends above the
  limits, the usual test after one that ends within them.

| case | before (c40cb98) | hold | hold μ | separate | hold above the limits |
|---|---|---|---|---|---|
| KKT problem, `inner_evaluations = 20` | `:interior`, 97, 7.0e-5 | `outer` (100) reached, 1,958, 5.8e-5 | `outer` reached, 1,960 | as before | as before |
| the same, 10 | `:interior`, 75, 4.8e-4 | `:interior`, 273, 2.6e-6 | 147, 6.5e-5 | as before | as before |
| the same, `iterations = 3` | `:interior`, 254, 0.24 | 524, 1.6e-2 | 270, 3.2e-2 | as before | 187, 0.23 |
| quadratic program, `iterations = 2`, phase 1 off | `:infeasible` (false verdict), 213 | `:interior`, 194, 4.1e-6 | 187, 4.1e-6 | 1,121, 2.4e-4 | `:interior`, 170, 4.0e-5 |
| the same, `inner_evaluations = 5` | `:infeasible` (false verdict), 105 | `:interior`, 283, 4.0e-6 | 251, 3.9e-6 | `outer` reached, 500 | `outer` reached, 500, 1.8e-4 |
| conflict, `iterations = 3`, phase 1 off | `:infeasible_box`, 54 | 268 | 403 | 54 | 268 |
| conflict, `inner_evaluations = 10` | `:infeasible_box`, 58 | 94 | 91 | 58 | 94 |

The false verdicts are the study's defect on a feasible problem: capped runs that end
above the limits raise the penalty to its ceiling. Holding the penalty there removes
them; holding it also within the limits leaves the dual ascent at the first penalty,
too slow for the run to converge where every inner run is capped (the KKT problem at
20 evaluations), while the usual test there converges as before: within the limits V is
the complementarity, which the penalty's growth sharpens, and no verdict can follow from
a design within the limits. *Separate* keeps the verdict away but not the growth (10⁹ on
the quadratic program, its design 60 times farther from the KKT point). Chosen: hold
above the limits; the cost is a slower verdict on a true conflict (268 against 54
evaluations at `iterations = 3`), which phase 1 (on by default) gives in 20. With the
caps above every inner run (the default `iterations` and no `inner_evaluations` on every
other test) nothing changes. Confidence: closed forms, one start each; the rule is an
engineering reading of the method's premise, not a published rule.

### 3. Phase 1's budget

Phase 1's single L-BFGS-B run took the inner runs' cap, `iterations`. On the study
variant with `iterations = 100` it ended on that cap, neither converged nor stalled: no
verdict, the display read "no design within the limits", and the main run started from
the unchanged start (the code adopted phase 1's design only when it met the limits; a
capped phase 1 was not foreseen). With a cap ten times larger the same phase 1 reached
its verdict. `iterations` caps the augmented
Lagrangian's subproblems, so that the multipliers change often; phase 1 has none, and
stops only by its convergence, its stall rule (1% over 100 evaluations) and the
run's `evaluations`, time limit and callback. Phase 1 is no longer capped by
`iterations`; it now ends at a verdict, within the limits, or stopped by the budget, the
time limit or the callback, which ends the whole run at its design of least violation
(§ 1), so the main run never starts from a phase 1 short of both. `phase1.reason` names
the end (`:feasible`, `:converged`, `:stalled`, or the stop), and the display shows when
a stop came before a verdict.

Closed form (the 120 TW cap against the drive floor, `feasibility = 1e-8`; scratch
script, one thread): with `iterations = 3`, c40cb98's phase 1 stopped after 7
evaluations without a verdict and the augmented Lagrangian reached `:infeasible_box`
after 12 outer iterations (61 evaluations); now phase 1 gives the verdict in 20, the
evaluations and the design of the default `iterations` exactly (violations 0.0237667 and
0.0252364 of their scales, as the main run's). Runs whose phase 1 ended within its
iteration cap are unchanged (every constrained run of the four test files above).
`test/constraints.jl` checks the verdict at `iterations = 3` and the stop at a budget of
5 evaluations, with its line in the display.

### 4. The laser's limit between the phase samples

`bounds(phased.param; power_max)` caps each power scale so that its beams stay within
the laser's limit at the phase's time; the pulse a design is delivered as,
`scaled_pulse` of the beam's base pulse with the scales as multipliers at the scaled
phases' times, interpolates the multipliers between those times. Where the base pulse
rises between two samples the product exceeds the limit with every sample inside it: on
the mock-up of `test/pulse_shapes.jl` (a base pulse at its 2.0 TW limit at 5.2 ns,
phases at 1, 2, 3 and 5 ns, the multiplier ramped to 1 by 5.5 ns) multipliers of 1.05,
1.05, 1.05, 1.03 give 2.036 TW at 5.2 ns, and the corners of the box itself give 2.13 TW
(between the first two phases, where the base pulse is linear, the product of two rising and
falling linear pieces, 2.08 TW). On the multi-phase study a scale of the last scaled
phase carried a cone's pulse past its laser limit between that phase and the peak
(scratch runs).

The layers considered were a clip in `scaled_pulse` or the table writer (changes the pulse's
energy, so the delivered drive is no longer the optimized one: rejected); a check that
reports (kept: `scaled_pulse(...; power_max)` and `shape_pulses(...; power_max)` throw
an `ArgumentError` naming the time, the power and the excess); and the design's limits.

*Superseded the same day* (§ 4a): the first version bounded each scale by
`mᵢ ≤ L / max(base on [tᵢ₋₁, tᵢ₊₁])` (the pulse's start before the first phase, `peak`
or its end after the last), a box that is sufficient but excludes designs whose written
pulse is within the limit.

### 4a. The exact condition as a constraint

The written pulse (`scaled_pulse`) is tabulated at the base pulse's points and the phase
times and is linear between them, so its maximum lies at a point. At a point τ between
the phases at tᵢ and tᵢ₊₁ the multiplier is `(1 − w) mᵢ + w mᵢ₊₁`, `w = (τ − tᵢ)/(tᵢ₊₁ −
tᵢ)`, so the condition `base(τ) ((1 − w) mᵢ + w mᵢ₊₁) ≤ L` is linear in two neighboring
scales; before the first phase it bounds only m₁, and after the last only mₙ (on the
`peak` ramp `(1 − w) mₙ + w ≤ L/base(τ)`; the pulse from `peak` on has no scale).
Between two phases the line `(1 − w) mᵢ + w mᵢ₊₁` must lie below every point
`(w, L/base(τ))`, which holds if and only if it lies below the vertices of their lower
convex hull: those points carry the whole condition.

The box excluded what it should not. On a multi-phase re-pointing study (scratch runs) a
cone's base pulse rises steeply between the last two scaled phases and peaks at the
laser's limit inside the `peak` ramp. The box capped the last scale at 1, which is exact
(the ramp then crosses the base pulse's maximum), and also the scale before it at the
limit over the base pulse at the end of its interval, although by then the delivered
multiplier is interpolated down to the last scale's: every run under the box ended
infeasible, a drive floor short, while a box that applied the interval rule to the last
scales only converged with every limit met and the written pulses within the limit. No
box is exact here: the condition couples two scales, and a box excludes the corner of
the feasible polygon that a design trading one scale against its neighbor needs.

Options weighed:

- *Only a tighter box*: the exact caps on one scale (the sample times, the pulse before
  the first phase, the ramp or the hold after the last) are necessary, not sufficient,
  for a pulse that rises between two phases away from both; making them sufficient puts
  back the exclusion. Kept as the box part.
- *One component per key and interval*, its value the largest written power over the
  interval: exact, but a maximum over points, with a kink wherever the point of the
  maximum changes, which the augmented Lagrangian's L-BFGS-B handles poorly.
- *A smooth stand-in* (a p-norm over the interval, the exact value reported, as
  `irradiance_ceiling`): smooth, but the working limit has to be lowered after
  convergence, and the stand-in is not needed when the exact condition is linear.
- *One component per hull vertex* (chosen): each linear in two scales, so smooth and
  exact, and few: on the mock-ups of `test/pulse_shapes.jl` 5 and 6 components for 9
  table points, and 16 for a 141-point random pulse (keys sharing a pulse share the
  points).

`pulse_ceiling(scenario; pulses, times, power_max, peak)` builds it: a `Constraint`
named `:pulse`, a component per hull vertex and per key, the value the written power at
the point (TW) against `power_max`, labeled by key and time; it takes the scales from the
evaluation through `power_scales` (generic in the number type; a case in
`test/inference.jl`). `bounds(...; pulses, times, peak)` now holds the conditions on a single
scale, and on each interval the cap that follows with the neighboring scale at its
least value in the box (both exact, since they exclude no design within the limit). Without
`pulses` the box is unchanged. The augmented Lagrangian holds the coupled components to
within the run's `feasibility` (1e-3 of the limit by default); a ceiling at
`(1 − feasibility) L` leaves the written pulse within `L`.

The table is the deliverable: the product of the two linear pieces (base and multiplier)
between table points is quadratic and can exceed the table, which `scaled_pulse` does
not write (on the rise mock-up below, by at most 0.0004 TW at 2 TW over 300 random
designs; at the corners of the former box on the mock-up above, 0.08 TW).

Checks (`test/pulse_shapes.jl`, against the written pulse sampled every 1e-4 ns, every
table point on that grid): on 300 random designs for each of three mock-ups (a steep
rise between the last two phases to a peak at the limit inside the ramp, a hump between
two phases away from both, a 141-point random pulse), a design meets the constraint if
and only if its sampled pulse stays within the limit (239, 240 and 99 within), and the
box contains every design within the limit; the analog of the study's case (the
earlier scales at 1.06, the last at 0.99, where the former box capped the third at
1.031) is admitted with its pulse within the limit, and a design over the limit at the
hump, inside the box, is reported at that point. `test/phased.jl` runs the largest
scales under the constraint and the box (a linear program, closed form `m_picket = 7/3`,
`m_step = 1.5`, above the former box's 2.0 for the picket).

On the study (scratch runs, 8 threads; values not recorded here): the per-beam rung
rerun from its converged design with `pulse_ceiling` and the new box in place of the
hand-made box converged inside the box with every limit met, in 3% fewer evaluations
than the hand-made run, to the same scales at three decimals and the same written
pulses; the constraint (48 components for two cone keys on measured pulses) ended
satisfied, its closest component 0.6% below the limit (the last scale set below its cap
by its drive floor and the objective, not by the laser's limit). The former box had capped the scale
before the rise below what the design needs, and every run under it ended infeasible.

### 5. On the study (scratch runs, 8 threads; values not recorded here)

- The infeasible variant as first run (phase 1 off, `iterations = 100`, an inner
  evaluation cap): every inner run but one ended on its cap above the limits and held
  the penalty, so the run went to its evaluation budget instead of the verdict at the
  penalty's ceiling, and returned an iterate from early in the run (outer iteration 3 of
  28) whose ceiling excess and largest violation are below the former return's, the
  gates' about as far over. It still lies inside the keep-out on many beams: by the
  normalized measure the start, far over its gates, is the more violating design (§ 6).
- The same with phase 1 (the default) and `iterations = 100`: the verdict in a few
  hundred evaluations, as the former run reached only with a cap ten times larger.
- The converging variant continued from its converged design: the same design and
  evaluation count as on the base tree, exactly (11 outer iterations; the 9 capped
  inner runs ended within the limits, where the usual penalty test applies), within
  0.01 mm rms of the design it started from.

### 6. The clearance's unit

`clearance_floor`'s `scale` (1 mm) makes 0.1 mm of clearance weigh as a gate exceeded
by 10% of itself, each beam a component. Measured with the scale at 0.05 mm (the
pointing tolerance), everything else as in §§ 1–3:

- the infeasible variant: phase 1's verdict trades the other way, the keep-out
  shortfall smaller and the gates and the ceiling several times further over (the
  largest normalized violation, in gate units, ten times the 1 mm run's); the phase-1-off
  run likewise. Neither weighting is the better design: which shortfall to accept is
  the study's.
- the dropped-bundle tutorial (a script run): its first stage converges in 767
  evaluations instead of 1,590 but at 0.9499 of the 48-quad drive instead of 0.9574 (the
  keep-out held 50 µm clear of the caps instead of 15 µm, against a 10 µm minimum), and
  the budget-stopped polish then ends at its start's drive (0.9489); the runs without a
  keep-out do not move.
- the constrained runs of `test/constraints.jl`, `test/phased.jl`,
  `test/picket_fluence.jl` and `test/pointing_summary.jl`: two move, the one-parameter
  keep-out run (the same design and value in 99 evaluations instead of 75) and the
  held-scales verdict (the same design, its violation 0.111 in the new unit instead of
  0.0055); every other run is unchanged (`test/phased.jl`'s floors lie far below the
  designs' clearance), and two assertions on the scale itself fail, as expected.

Decision: the default stays 1 mm (nothing in the notes argues for it, and nothing
measured argues for 0.05 mm: on the conflicting case the weighting moves the violation
between the limits, and on the tutorial it costs 0.8% of drive for clearance the
keep-out does not require). `clearance_floor`'s docstring now states the weighing.

### 7. An exact ceiling in a run stopped by its budget

The question from the dropped-bundle tutorial's budget-stopped polish ("The
dropped-bundle tutorial's long runs under rounding"): its wall peak ended over the
ceiling on some starts while the run held the ceiling's smooth stand-in. The reading
holds in part. The reported rows carry the exact peak: `_constraint_rows` reports
`_exact_value` (the peak over the elements) for an exact ceiling, and `in_phase`
forwards it, so `PointingResult.constraints` gives the peak and its violation as they
stand. The stand-in's working limit is lowered only after a run has converged on the
stand-in (`irradiance_ceiling(...; exact = true)`), so a run stopped earlier never
corrects for the stand-in's offset from the peak. Since § 1, a run that ends above its
limits returns its iterate of least violation measured on those exact rows.

Lowering the working limit at every outer iteration from the exact value at the iterate was
not built. The lowering only ratchets (a peak below the limit leaves the working limit, §
"Constrained re-pointing: an irradiance ceiling and mode gates", exact peak, to avoid
cycling where the peak moves between elements), and the stand-in-to-peak ratio of an early,
unconverged iterate is not the converged design's, so early rounds would hold the converged
design below the limit by the difference, a drive cost that later rounds cannot return.

### 8. The dropped-bundle tutorial's polish under § 1

The per-beam polish (120 evaluations, always stopped by the budget) from six starts:
the start of the tutorial and five of the perturbed starts of "The dropped-bundle tutorial's
long runs under rounding" (10⁻⁹ mm, seeds 1 and 2; 10⁻¹² mm, seeds 3 to 5), on the tree
with the f-cone model's shared coordinates (8 threads):

| start | returned | drive gained over the first stage on beams | departure / gate, largest | peak / ceiling | caps, closest |
|---|---|---|---|---|---|
| the tutorial's | last iterate (within the limits) | +0.17% | 0.998 | 0.966 | 34 µm |
| 10⁻⁹, seed 1 | an earlier one within the limits | −0.002% | 0.954 | 0.931 | 42 µm |
| 10⁻⁹, seed 2 | an earlier one | 0.000% | 1.000 | 0.959 | 49 µm |
| 10⁻¹², seed 3 | the start (within the limits) | 0 | 0.948 | 0.927 | 49 µm |
| 10⁻¹², seed 4 | last iterate | +0.32% | 1.000 | 0.995 | 17 µm |
| 10⁻¹², seed 5 | an earlier one | −0.002% | 0.992 | 0.942 | 49 µm |

Before § 1 the same starts ended 0.966–1.020 of the ceiling with gates up to 0.2%
over; now every returned design meets every limit, the peak at 0.93–0.995, and where
the last iterate exceeded a gate the run gives up the polish's drive gain for a design
within the limits. The tutorial's text now states so, and the tutorial asserts the peak
again, at most 1% over the ceiling (its largest over the six starts 0.995). The
start of the tutorial prints the former output to every digit (its last iterate met the
limits).
