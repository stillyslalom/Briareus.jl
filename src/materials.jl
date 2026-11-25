# This file defines `Material` and the material presets.

"""
    Material(name; albedo = 0.0, conversion = 0.0)

Surface properties for the radiation model. `albedo` (a number or `AlbedoModel`) is the
x-ray re-emission fraction; `conversion` (a number or a [`Conversion`](@ref)) is the
fraction of absorbed laser power re-emitted as x-rays.
"""
struct Material
    name::Symbol
    albedo::AlbedoModel
    conversion::Union{Float64,Conversion}
    Material(name, albedo, conversion) =
        new(name, albedo, conversion isa Conversion ? conversion : Float64(conversion))
end

Material(name::Symbol; albedo = 0.0, conversion = 0.0) =
    Material(name, albedo isa AlbedoModel ? albedo : ConstantAlbedo(albedo), conversion)

albedo(m::Material, t = 0.0) = albedo(m.albedo, t)

"""
    gold(; albedo = 0.8, conversion = 0.9)

Gold hohlraum wall for static solves. Its laser-to-x-ray `conversion` is a constant 0.9
whatever the `albedo`. Both defaults are the inputs of a published static VisRad model of
a NIF hohlraum (Young et al., PoP 27, 082702 (2020)). The conversion lies within the
0.85–0.95 published for NIF hohlraums near peak power, above the 0.68–0.85 of 1 ns Nova
and Omega hohlraums, which `gold(conversion = 0.75)` suits. The research notes section
"Material defaults and their sources (2026-11-18)" gives the sources. For a drive history
on NIF, use [`nif_gold`](@ref), whose albedo and conversion follow the drive. For a
directly irradiated foil, use [`foil_material`](@ref).
"""
gold(; albedo = 0.8, conversion = 0.9) = Material(:Au; albedo, conversion)
"""
    uranium(; albedo = du_wall_loss(), conversion = 0.9)

Au-lined depleted-uranium hohlraum wall: [`du_wall_loss`](@ref)`()` and the constant
conversion 0.9 of [`gold`](@ref) (the conversion happens in the gold liner). For a
drive history on NIF use [`nif_uranium`](@ref).
"""
uranium(; albedo = du_wall_loss(), conversion = 0.9) = Material(:DU; albedo, conversion)
"""
    nif_gold(; albedo = WallLossAlbedo(), conversion = nif_conversion())

The gold wall of a NIF drive history: Lindl's Marshak wall loss, evaluated by
[`drive_history`](@ref) at the temperature history of each element, and the rising
[`nif_conversion`](@ref)`()`, `0.7 t^0.12 ≤ 0.95` with `t` the time since the pulse
start (at `solve`'s default `time = 0` it is 0.53, so a static solve needs
[`gold`](@ref)).
"""
nif_gold(; albedo = WallLossAlbedo(), conversion = nif_conversion()) = Material(:Au; albedo, conversion)
"""
    nif_uranium(; albedo = du_wall_loss(), conversion = nif_conversion())

The Au-lined depleted-uranium wall of a NIF drive history (HYBRID-E's):
[`du_wall_loss`](@ref)`()` and the rising [`nif_conversion`](@ref)`()` of its gold liner.
"""
nif_uranium(; albedo = du_wall_loss(), conversion = nif_conversion()) = Material(:DU; albedo, conversion)
"Capsule ablator (defaults: albedo 0.3, conversion 0.3)."
ablator(; albedo = 0.3, conversion = 0.3) = Material(:ablator; albedo, conversion)
"Perfect absorber that emits nothing."
const BLACK = Material(:black; albedo = 0.0, conversion = 0.0)
