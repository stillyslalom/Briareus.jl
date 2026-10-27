# Mode analysis

Legendre moments ([`legendre_moments`](@ref)) and real spherical harmonics
([`ylm_moments`](@ref), [`flux_ylm`](@ref)) are integrated over each facet's exact solid
angle; [`mode_amplitude`](@ref) gives the fractional amplitude of each mode (the mode-1
drive asymmetry, for example). [`energy_balance`](@ref) accounts for the source, wall and
capsule absorption, LEH loss and leakage, and [`radiation_temperature`](@ref) converts
flux to temperature.

```@autodocs
Modules = [Briareus]
Pages = ["src/analysis.jl"]
Private = false
```
