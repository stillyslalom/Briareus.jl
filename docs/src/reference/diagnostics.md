# Diagnostic definitions

A synthetic diagnostic is described by a definition file: one TOML file per instrument,
holding its line of sight, its hardware and its recording, each value with its source.
The public definitions ship with the package (`data/diagnostics`); a user group keeps
its own, including details that are not public, in directories on a search path, and
modifies a public definition by extending it rather than copying it.
[`diagnostic`](@ref) builds the instrument of a definition (for the soft x-ray power
kind, a [`DanteView`](@ref) for [`dante`](@ref)); [`check_diagnostics`](@ref) reports
every problem in the files on the search path, by file and key.

## The search path

[`diagnostic_paths`](@ref) lists the directories read, in order of precedence: the
built-in definitions, the directories saved in the user configuration file
([`diagnostic_config_path`](@ref), key `diagnostics.paths`), those of the environment
variable `BRIAREUS_DIAGNOSTICS` (for batch jobs; `;`-separated on Windows,
`:`-separated elsewhere) and those added for the current session. Directories are read
recursively. A definition later on the path shadows an earlier one with the same `id`,
with a warning naming both files. A directory is added once and kept for later
sessions with

```julia
add_diagnostic_path("//groupserver/briareus/diagnostics")
```

and [`remove_diagnostic_path`](@ref) removes it. The built-in constructors
([`dante1`](@ref), [`dante1_responses`](@ref)) read the built-in definitions only, so
user files do not change them; `diagnostic("nif_dante1")` takes the search path.

## A definition file

Every file has a header: `schema` (the format version, [`DIAGNOSTIC_SCHEMA`](@ref)),
`id`, and either `kind` or `extends`; optionally `name`, `facility`, `public`
(`false` for non-public detail), `valid_from` and `valid_to` (dates), `shots`, and the
provenance keys `source`, `confidence` (H, M or L, as in the research notes) and
`note`, which every table may also carry. Keys are checked: an unknown key or a value
of the wrong type is an error naming the file and the key. Units are in the key names
(`theta_deg`, `energy_eV`, `thickness_um`, `solid_angle_sr`).

The kind `xray_power` (Dante, Sentinel) takes:

| Table | Keys |
|---|---|
| `[line_of_sight]` (required) | `theta_deg`, `phi_deg`: the instrument's direction from TCC |
| `[view]` | `field_of_view_mm` (the collimated field, a disk of this diameter normal to the line of sight; without it the view looks through `through` alone), `center_part` or `center_mm` (the field's center, a part's centroid or a point; default TCC), `through` (the aperture looked through, and the reference for Tr; default `leh_bottom` without a field of view, none with one), `pixels` (rays per side; default 400), `area_mm2` (the area Tr is quoted over; default the aperture's initial area) or `area_part` (Tr over that part's visible projected area) |
| `[recording]` | `mode = "history"` |
| `[[channel]]` (one per channel) | `number`, `energy_eV` (nominal band energy) and `solid_angle_sr` (required); `filters` (an array of `{ material, thickness_um }`, optionally `density_g_cm3`, or `{ transmission_csv }`); `grid` (open fraction; default 1); `mirror` (`{ kind = "SiO2", angle_deg }`, `{ kind = "mlm", angle_deg = 1.9 }` or `{ reflectivity_csv }`); `cathode` (`"Al"`, `"Cr"` or `"Ni"`, [`xrd_sensitivity`](@ref)) or `sensitivity_csv`; or a measured `response_csv` in place of all of these |
| top level | `remove_channels` (numbers dropped from the definition extended) |

A field of view takes in every surface of the target facing the instrument, occluded by
any part: a NIF Dante's is ±30 mm at TCC (the 2004 design), wider than any hohlraum, so
windows, packages and shields enter the signal as they would on the shot. Within an
extended definition, giving one of `area_mm2` / `area_part` or `center_part` /
`center_mm` drops the other inherited from the base.

The response of a channel is the product of its parts ([`DanteResponse`](@ref)).
Measured curves are CSV files with a header line and two columns, photon energy (eV)
and value (a fraction, or V/GW for a sensitivity or a response), at a path relative to
the definition file ([`TabulatedCurve`](@ref)); outside a table a factor holds its end
values and a whole response is zero, and [`check_diagnostics`](@ref) warns when a table
does not cover 50 eV – 20 keV.

## Extending a definition

A file with `extends = "<id>"` takes that definition and replaces what it gives:
tables key by key, channels by `number` (a channel's listed keys replace the base
channel's; a new number adds a channel), and `remove_channels` drops channels. A
campaign that replaced Dante 1's channel 6 filter, removed channel 16 and added a
calibrated channel is, in full:

```toml
schema = 1
id = "group_dante1_campaign"
extends = "nif_dante1"
name = "Dante 1, campaign configuration"
public = false
remove_channels = [16]

[[channel]]
number = 6
filters = [{ material = "Cu", thickness_um = 1.1 }]
source = "campaign shot sheet"

[[channel]]
number = 19
energy_eV = 1000
solid_angle_sr = 6.5e-7
response_csv = "calibrations/ch19.csv"
```

`DiagnosticDefinition.chain` lists the files a definition was resolved from, and
`public` is `false` when any of them sets `public = false`. Values a definition lacks
to be built (Dante 2's channels are not public, for example) are reported as
`:missing` by [`check_diagnostics`](@ref), with the key to supply.

## New kinds

[`register_diagnostic_kind`](@ref) adds a kind: a schema of
[`DiagnosticField`](@ref)s and a builder. The format already reserves the recording
modes of the kinds planned next (time-integrated, gated and streaked imagers,
radiography, streaked optical diagnostics); see `reference/diagnostics-design.md` in
the repository.

```@autodocs
Modules = [Briareus]
Pages = ["src/diagnostics.jl"]
Private = false
```
