# Input and output

Beam tables in a canonical CSV ([`write_beam_table`](@ref), [`read_beam_table`](@ref)),
with column and unit mapping for foreign tables; STL ([`read_stl`](@ref),
[`read_stl_solids`](@ref), [`write_stl`](@ref)); and [`write_results`](@ref), which
exports per-element laser, radiosity and flux, per-component and per-beam power, and the
capsule's Legendre and ``Y_{lm}`` modes as CSV for postprocessing. VisRad files are on
[VisRad interchange](@ref).

```@autodocs
Modules = [Briareus]
Pages = ["src/io.jl"]
Private = false
```
