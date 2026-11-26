# Scripting reference for language models

This page is a compact map of the package for a language model that writes Briareus
scripts, and for anyone who reviews such a script. Every name below is exported and
links to its docstring. The same text without the links is the file `llms.txt` at the
root of the repository.

## Units and conventions

- Lengths are in mm, times in ns, powers in TW, fluxes and intensities in TW/mm²
  (multiply by 1e14 for W/cm²), temperatures in eV. Angles are in degrees at the API.
- The chamber frame has z up. A port is a unit vector from the chamber center, and a
  beam travels along `k = -port` ([`direction`](@ref)).
- A [`Beam`](@ref) has fixed `port`, and `aim` translates its axis. Best focus lies at
  `aim + defocus * k`. Beams are immutable: `remake(beam; aim = (x, y, z))` returns a
  changed copy.
- A [`Target`](@ref) is a list of [`Component`](@ref)s, each a mesh with a
  [`Material`](@ref). Elements are radiating patch sides. A one-sided surface radiates
  from the front, and its back absorbs laser light as a loss.
- `Component.part` names the part of the user, which a CAD fit may split into pieces.
  Functions that take names take part names.
- A [`Scenario`](@ref) holds everything precomputed for one target and one set of ports.
  It stays valid while the ports are unchanged, for every aim, power and spot, so a
  sweep builds it once and calls [`evaluate`](@ref) many times.

## Workflows

A capsule in a hohlraum, with the drive, the modes and a pointing optimization:

```julia
using Briareus
target = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
beams = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.1)))
         for b in nif_quads(power = 8.0)]
scenario = Scenario(target, beams; observe = :capsule, viewfactor_kw = (; closed = true))
ev = evaluate(scenario, beams)
a = flux_moments(ev; lmax = 4)
Tr = radiation_temperature(mean_flux(ev))        # eV
p2, p4 = a[3] / a[1], a[5] / a[1]
v = hotspot_velocity(ev)                         # km/s, chamber x, y, z

pointing = GroupPointing(beams; by = :group, fields = (:z,))
tuned = Scenario(target, pointing; observe = :capsule, viewfactor_kw = (; closed = true))
result = optimize_pointing(tuned, symmetry_objective())
result.beams                                      # the optimized beams
```

The names used are [`hohlraum`](@ref), [`nif_quads`](@ref), [`flux_moments`](@ref),
[`mean_flux`](@ref), [`radiation_temperature`](@ref), [`hotspot_velocity`](@ref),
[`GroupPointing`](@ref), [`optimize_pointing`](@ref) and [`symmetry_objective`](@ref).
A shaped pulse goes through [`PhasedScenario`](@ref) and [`drive_history`](@ref).

The NIF target-configuration rules and the 1ω review, on the 192 beams:

```julia
nb = nif_beams(beams)                             # the four beams of each quad
report = talis_check(target, nb; capture = (:wall,))
bad = violations(report)                          # failing and review rows, worst first
row = first(x for x in report.rows if x.rule === :clearance_1w)
approach = omega1_closest_approach(target, nb, report, row)
approach.margin, approach.x_hw, approach.x_beam   # mm, the hardware and beam points
review = omega1_review(target, nb; capture = (:wall,), threshold = 1.0)
write_omega1_review("review", review)             # a self-contained HTML packet
```

The names used are [`nif_beams`](@ref), [`talis_check`](@ref), [`violations`](@ref),
[`omega1_closest_approach`](@ref), [`omega1_review`](@ref) and
[`write_omega1_review`](@ref). [`omega1_disposition`](@ref) gives the verdict per mesh
triangle, [`corridor_map`](@ref) the directions that beams and lines of sight leave
open, and [`clearance`](@ref) the worst-case clearance check with engineering defaults.

Which beams a shot can lose, for any target:

```julia
choice = roe_default_criteria(scenario, beams)    # a starting set, engineering defaults
roe = roe_scan(scenario, beams, choice.criteria)
write_roe_packet("roe", roe, beams; choice.notes)

custom = ROECriterion(:wall_power, ev -> sum(ev.deposition.element); green = 0.03, amber = 0.06,
                      measure = :loss, label = "fraction")
roe_scan(scenario, beams, [roe_drive(), roe_mode(2), custom])
```

The names used are [`roe_default_criteria`](@ref), [`roe_scan`](@ref),
[`write_roe_packet`](@ref), [`ROECriterion`](@ref), [`roe_drive`](@ref) and
[`roe_mode`](@ref). The other built-in criteria are [`roe_hotspot_velocity`](@ref),
[`roe_yield`](@ref), [`roe_face_flux`](@ref), [`roe_face_nonuniformity`](@ref),
[`roe_temperature`](@ref), [`roe_flux_balance`](@ref), [`roe_foil_irradiance`](@ref),
[`roe_foil_nonuniformity`](@ref), [`roe_foil_centroid`](@ref),
[`roe_backlighter_signal`](@ref), [`roe_illumination_rms`](@ref) and
[`roe_corona_rms`](@ref).

A backlighter foil with no hohlraum, and direct drive on a sphere:

```julia
ft = Target(foil(:backlighter, (0, 0, 0), (0, 0, 1); width = 2.0, n = 16))
fb = [remake(b; aim = (0.0, 0.0, 0.0)) for b in omega_beams(; cones = (:cone1,))]
fs = Scenario(ft, fb)                             # laser only: no observe, no radiation
absorbed = sum(evaluate(fs, fb).deposition.element)          # TW
roe_scan(fs, fb, [roe_foil_irradiance(:backlighter), roe_foil_centroid(:backlighter)])

ct = Target(Component(:capsule, icosphere_mesh(0.5; n = 4), Material(:CH); surface = SphereSurface(0.5)))
ob = omega_beams(; cones = :all)
dep = evaluate(Scenario(ct, ob), ob).deposition.element
W = ylm_weights(ct, 1:nelements(ct); lmax = 8)
sigma = mode_rms(W * (dep ./ ct.elem_area), 8)    # RMS nonuniformity of the illumination
```

The names used are [`foil`](@ref), [`omega_beams`](@ref), [`icosphere_mesh`](@ref),
[`SphereSurface`](@ref), [`ylm_weights`](@ref) and [`mode_rms`](@ref).
[`backlighter_signal`](@ref) gives the photons of a line on a detector, and
[`CoronaOperator`](@ref) the absorption in a direct-drive corona.

A VisRad workspace, read and written back:

```julia
ws = read_visrad_workspace("shot.vrw")
target = visrad_target(ws)
beams = ws.beams
moved = [remake(b; aim = b.aim + [0.0, 0.0, 0.1]) for b in beams]
write_visrad_workspace("shot_moved.vrw", "shot.vrw", moved)
diff = design_diff(target, beams, target, moved; talis = (; capture = (:wall,)))
write_design_diff("diff", diff)
```

The names used are [`read_visrad_workspace`](@ref), [`visrad_target`](@ref),
[`write_visrad_workspace`](@ref), [`design_diff`](@ref) and
[`write_design_diff`](@ref). [`visrad_unconverted`](@ref) gives the 1ω and 2ω powers of
the workspace, and [`compare_visrad`](@ref) and [`visrad_residuals`](@ref) compare a
solution with the results that VisRad saved.

The desktop application saves a design as a session file, and File ▸ Export Julia script
writes it as a script that defines `target`, `beams`, `talis`, `scenario` and
`evaluation`. `briareus run script.jl` runs a script with the package loaded, and
`briareus design.briar --tab talis` opens a session written by a script on a results
tab.

## Traps

- `Beam.dlambda` is in Å at the wavelength of the beam, 3ω on NIF. A NIF detuning quoted
  at 1ω is divided by 3 before it is passed.
- A [`Scenario`](@ref) at a time partway through a pulse needs the state of the history:
  with `st = history_state(target, h, t)`, the keywords `st.time` and `st.albedo` go to
  `Scenario` ([`history_state`](@ref), on a [`drive_history`](@ref) run with
  `keep_albedo = true`). The default albedos evaluate a [`WallLossAlbedo`](@ref) as if
  the wall had been at one temperature since `t = 0`, and the capsule drive then differs
  from the history by tens of percent.
- [`VisRadProfile`](@ref) cuts the beam off sharply by default, as VisRad does. An
  optimization passes `VisRadProfile(; taper = 1)`, since a gradient misses the step of
  the sharp cut.
- A [`Scenario`](@ref) is rebuilt when a port, the target or the observed part changes.
  Aims, powers, defocus and spots change through [`evaluate`](@ref) alone.
- [`talis_check`](@ref) counts every surface as hardware unless `hardware` names the
  parts, so a capsule or an end cap in the beam path fails the 3ω clearance.
- The rules-of-engagement thresholds, the clearance thresholds and the line-of-sight
  volumes of [`corridor_map`](@ref) are engineering defaults. A shot sheet needs the
  criteria of the experiment.
- Checks are validated against independent references. A script that checks a result
  against the function that produced it checks nothing.
- Objectives passed to [`optimize_pointing`](@ref) run on several threads and on dual
  numbers: they write no shared state and return the number type of their input.
