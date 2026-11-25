"""
    Briareus.GUIAPI

The internals of Briareus that the desktop application BriareusGUI calls, under names
without the underscore. The application takes everything else from the exported API.
This list is the whole boundary between the two packages: a name is added here when the
application needs another internal, and a test of the application checks that its source
has no other `Briareus._name` reference. The names are for that package. They are not
exported and can change with any release.

- Frames and placement: `step_references`, `step_axis`, `is_relative_step`,
  `axis_rotation`.
- Targets, beams and pulses: `profile_radius`, `beam_key`, `beam_pulse`, `is_nif`,
  `nif_beam_table`, `omega1_beamlets`, `omega1_square`.
- Drive and unconverted light: `has_radiation`, `observed_flux`, `UnconvertedLight`,
  `unconverted_at!`, `part_series`.
- Diagnostics and files: `BUILTIN_DIAGNOSTICS`, `read_definitions`,
  `checked_definitions`, `definitions`, `resolve_definition`, `saved_paths`,
  `imager_view`, `map_csv!`.
"""
module GUIAPI

using ..Briareus: _step_refs as step_references, _step_axis as step_axis, _relative as is_relative_step,
                  _axis_rotation as axis_rotation,
                  _profile_radius as profile_radius, _beam_key as beam_key, _beam_pulse as beam_pulse,
                  _is_nif as is_nif, _nif_beam_table as nif_beam_table, _omega1_beamlets as omega1_beamlets,
                  _omega1_square as omega1_square,
                  _radiative as has_radiation, _observed_flux as observed_flux,
                  _UnconvertedLight as UnconvertedLight, _unconverted_at! as unconverted_at!,
                  _part_series as part_series,
                  _BUILTIN_DIAGNOSTICS as BUILTIN_DIAGNOSTICS, _read_all as read_definitions,
                  _checked as checked_definitions, _definitions as definitions, _resolve as resolve_definition,
                  _saved_paths as saved_paths, _imager_view as imager_view, _map_csv! as map_csv!

end # module GUIAPI
