# Shot review

Every check is judged on the worst case over pointing, target positioning, hardware
tolerance and focus errors, never on the nominal beams. [`clearance`](@ref) gives the
worst-case beam-to-hardware clearance, power and intensity, with intensity-contour
envelopes, clearance-only spot wings ([`nif_wings`](@ref)) and pass, review and fail
tiers; its thresholds and wing models are engineering defaults, which the
review requirements and measured far fields of a facility should replace.
[`talis_check`](@ref) applies NIF's published target-configuration rules (3ω capture and
incidence, 3ω and 1ω clearance, unconverted-light management, the TAS stay-in zone), and
[`unconverted_light_check`](@ref) covers the 1ω and 2ω footprints and the 32 mm zone, and
[`unconverted_deposition`](@ref) puts their power on the target (a uniform square beam
per beamlet and color, with shadowing; the unconverted fraction is an input).
[`unconverted_history`](@ref) follows that light through the pulses, with a fraction and
a prepulse that may vary in time, and [`irradiance_history`](@ref) the 3ω light;
[`specular_check`](@ref) and [`preplasma_check`](@ref) compare them with the NIF User
Guide's 10⁸ W/cm² plasma-onset guideline for metals.
[`ApertureKeepOut`](@ref) puts the worst-case clearance of each beam from its entrance hole
into an optimization objective, [`SphereKeepOut`](@ref) its clearance from the capsule
and [`HardwareKeepOut`](@ref) its clearance from named hardware; each clearance is also a
separate quantity ([`keepout_clearance`](@ref)) for a hand-written objective or
constraint, as are the fractions of the laser power on named parts
([`clipped_fraction`](@ref) beside [`power_fraction`](@ref)). The TaLIS rules that carry
a margin are such quantities too, each checked against `talis_check`: the 3ω clearance
is a `HardwareKeepOut` with the rule tolerance, [`CaptureKeepIn`](@ref) gives the 3ω
capture margin ([`capture_margin`](@ref)), [`IncidenceRule`](@ref) the 3ω incidence
([`incidence_angle`](@ref)) and [`Omega1KeepOut`](@ref) the 1ω clearance
([`omega1_clearance`](@ref)). [`omega1_closest_approach`](@ref) gives the geometry behind
a 1ω clearance: the hardware point, the closest point of the footprint square, the
beamlet, focal plane and side of the square that hold it, and the terms between the
nominal distance and the margin ([`Omega1Approach`](@ref)). The square itself comes
from [`omega1_contour`](@ref) for the rules, the keep-out and the application alike.
[`omega1_reflection`](@ref) gives the specular reflection behind a 1ω management row: the
reflected ray, the open port it threatens and the projected area that reflects toward it
([`Omega1Reflection`](@ref)), and [`omega1_reflected_power`](@ref) turns the share of the
1ω light on that area into a reflected power. [`omega1_disposition`](@ref) gives the same rules as a field
over the hardware, one verdict per mesh triangle over all beams: clear, tight, curved
enough to scatter, clear by its reflection, to review, or failing
([`Omega1Disposition`](@ref)). [`omega1_slice`](@ref) gives the plane across the beamlet
through the hardware point of an approach as 2-D geometry, and
[`write_omega1_slice_svg`](@ref) draws it as a vector file without a graphics context.
[`corridor_map`](@ref) turns the same keep-outs around: it samples a sphere about a point
or a plane and gives the margin of each sample to the 3ω envelopes, the 1ω footprints and
the lines of sight of the diagnostics ([`SightLine`](@ref), [`sight_line`](@ref)), so the
open region shows where new hardware may be placed. [`corridor_margin`](@ref) gives the
margins of one point, and [`write_corridor_map_svg`](@ref) draws the map.
[`backlighter_foil`](@ref) and [`line_of_sight`](@ref)
place a foil on a diagnostic line of sight and check it against shields and beams, and
[`backlighter_signal`](@ref) follows the deposition through a He-α conversion efficiency
fitted to published foil data ([`kshell_ce`](@ref)) to the photons at the detector.

```@autodocs
Modules = [Briareus]
Pages = ["src/clearance.jl", "src/backlighter.jl", "src/backlighter_signal.jl", "src/keepout.jl", "src/unconverted.jl", "src/talis.jl", "src/omega1_slice.jl", "src/corridor.jl", "src/design_rules.jl"]
Private = false
```

## Beam rules of engagement

[`roe_scan`](@ref) plans which beams and quads a shot can lose before the shot day. An
outage is a set of beams at zero power with every other beam as planned, and each outage
is scored against a list of criteria ([`ROECriterion`](@ref)), each a function of an
evaluation with a green and an amber threshold. The scan covers every beam, every quad,
pairs of them and pre-planned companion drops, and sorts each unit into a tier: green
units may drop without a request, amber ones need approval, and red ones hold the shot.
The physics lies in the criteria, so the same scan applies to a capsule in a hohlraum
(`roe_drive`, `roe_mode`, `roe_hotspot_velocity`, `roe_yield`), a package face or the
wall of one cavity (`roe_face_flux`, `roe_face_nonuniformity`, `roe_temperature`,
`roe_flux_balance`), a backlighter foil (`roe_foil_irradiance`, `roe_foil_nonuniformity`,
`roe_foil_centroid`, `roe_backlighter_signal`) and direct drive (`roe_illumination_rms`,
`roe_corona_rms`), and any other function of an evaluation enters as a criterion of the
same kind. Several cases in one scan cover the epochs of a pulse and several observed
faces. The default thresholds of the built-in criteria are engineering defaults.
[`roe_table`](@ref), [`write_roe_csv`](@ref) and [`roe_markdown`](@ref) export the sheet
by role and tier, with the scores behind each tier and the pairs that are worse together.
[`port_map_svg`](@ref) draws the ports of the beams of any facility on two polar disks,
one per hemisphere, colored by the tiers of a report or by any value per beam or quad.
[`roe_default_criteria`](@ref) selects a starting set of criteria from what a design
contains.

## Review packets

A packet is a folder that a reviewer opens without the package. [`omega1_review`](@ref)
runs the TaLIS check and collects the 1ω disposition and the closest approach of every
failing or tight pair, and [`write_omega1_review`](@ref) writes it as one self-contained
HTML page with the slice drawings inline, the CSV of the check, one SVG per slice and
optionally JSON. [`write_roe_packet`](@ref) writes the rules-of-engagement sheet in the
same forms, with the chamber port map. The pages refer to no other file and to nothing
on the network, and the same input gives the same text.

[`design_diff`](@ref) compares two designs: the parts that moved, the beams that were
repointed, the TaLIS rows that newly fail or were fixed, and the physics of every role
that either design has. [`write_design_diff`](@ref) writes the page of the result.

```@autodocs
Modules = [Briareus]
Pages = ["src/roe.jl", "src/roe_criteria.jl", "src/roe_export.jl", "src/port_map.jl", "src/review_packets.jl",
         "src/design_diff.jl"]
Private = false
```
