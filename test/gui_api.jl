# The boundary with the desktop application (`Briareus.GUIAPI`): every listed name is
# bound, to the internal that its docstring line names.

using Briareus, Test

@testset "GUIAPI" begin
    api = Briareus.GUIAPI
    listed = (:step_references, :step_axis, :is_relative_step, :axis_rotation, :profile_radius, :beam_key,
              :beam_pulse, :is_nif, :nif_beam_table, :omega1_beamlets, :omega1_square, :has_radiation,
              :observed_flux, :UnconvertedLight, :unconverted_at!, :part_series, :BUILTIN_DIAGNOSTICS,
              :read_definitions, :checked_definitions, :definitions, :resolve_definition, :saved_paths,
              :imager_pixels, :map_csv!)
    @test all(n -> isdefined(api, n), listed)
    @test all(n -> occursin(string(n), string(@doc(Briareus.GUIAPI))), listed)
    @test api.step_references === Briareus._step_refs && api.has_radiation === Briareus._radiative
    @test isdir(api.BUILTIN_DIAGNOSTICS)
    @test api.is_nif(first(nif_beams())) && !api.is_nif(first(omega_beams()))
end
