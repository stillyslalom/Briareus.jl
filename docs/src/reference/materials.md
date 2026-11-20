# Material defaults

A `Material` contains the two surface properties the radiation model uses: the x-ray albedo
(re-emitted over incident flux, a number or an albedo model) and the laser-to-x-ray
conversion of the absorbed laser power (a number or a [`Conversion`](@ref)). The
constructors [`gold`](@ref), [`uranium`](@ref), [`nif_gold`](@ref),
[`nif_uranium`](@ref), [`ablator`](@ref) and [`BLACK`](@ref) build the package
defaults.

[`material_library`](@ref) lists named material defaults whose every quantity has a
published source or a derivation recorded in the research notes ("Material defaults and
their sources (2026-11-18)"). Each entry, a [`MaterialDefault`](@ref), contains the
`Material` it builds and, per quantity, a [`SourcedValue`](@ref): the value, how it was
obtained (one of [`MATERIAL_BASES`](@ref)), its source and its uncertainty. Displaying an
entry prints that record:

```julia
d = material_default(:gold)
display(d)                       # albedo and conversion, each with basis, source, uncertainty
t = hohlraum(; wall = d.material)
```

The entries named after a constructor build exactly the material that constructor builds,
identical in every field, so a model built from the library gives the same results. The foil
entries are the materials of [`foil_material`](@ref), the default of [`foil`](@ref) and
[`backlighter_foil`](@ref): gold from published disk and sphere data, the mid-Z
backlighter elements with the basis `:suggested` (a suggested range that no measurement
backs). Package defaults with no source found, such as the capsule conversion 0.3 in
[`ablator`](@ref), are not in the library.

```@autodocs
Modules = [Briareus]
Pages = ["src/material_library.jl"]
Private = false
```
