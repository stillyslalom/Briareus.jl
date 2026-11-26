# Outages and review packets

This page covers the beam rules of engagement and the self-contained pages that carry a
review to readers without the package. The clearance and TaLIS checks that the 1ω packet runs are on the
[Shot review](review.md) page.

## Beam rules of engagement

[`roe_scan`](@ref) plans which beams and quads a shot can lose before the shot day. An
outage is a set of beams at zero power with every other beam as planned, and each outage
is scored against a list of criteria ([`ROECriterion`](@ref)), each a function of an
evaluation with a green and an amber threshold. The scan covers every beam, every quad,
pairs of them, pre-planned companion drops and, when asked, triples of units that drop
safely in pairs, and sorts each unit into a tier: green
units may drop without a request, amber ones need approval, and red ones hold the shot.
The physics lies in the criteria, so the same scan applies to a capsule in a hohlraum
(`roe_drive`, `roe_mode`, `roe_hotspot_velocity`, `roe_yield`), a package face or the
wall of one cavity (`roe_face_flux`, `roe_face_nonuniformity`, `roe_temperature`,
`roe_flux_balance`), a backlighter foil (`roe_foil_irradiance`, `roe_foil_nonuniformity`,
`roe_foil_centroid`, `roe_backlighter_signal`) and direct drive (`roe_illumination_rms`,
`roe_corona_rms`, `roe_ablation_pressure`), and any other function of an evaluation enters as a criterion of the
same kind. Several cases in one scan cover the epochs of a pulse and several observed
faces. The default thresholds of the built-in criteria are engineering defaults.
A beam that stands for several facility beams, such as a NIF quad held as one beam, is
scanned as a quad and as one of its beams (`stands_for` of [`roe_scan`](@ref)).
[`roe_outage`](@ref) scores any other set of units, such as three beams or a quad with a
beam of another quad, against the criteria and baseline of a report.
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
