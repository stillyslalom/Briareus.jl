# Shot review

Every check is judged on the worst case over pointing, target positioning, hardware
tolerance and focus errors, never on the nominal beams. [`clearance`](@ref) gives the
worst-case beam-to-hardware clearance, power and intensity, with intensity-contour
envelopes, clearance-only spot wings ([`nif_wings`](@ref)) and pass, review and fail
tiers; its thresholds and wing models are engineering defaults, not facility criteria,
and a facility's own review requirements and measured far fields should replace them.
[`talis_check`](@ref) applies NIF's published target-configuration rules (3ω capture and
incidence, 3ω and 1ω clearance, unconverted-light management, the TAS stay-in zone), and
[`unconverted_light_check`](@ref) covers the 1ω and 2ω footprints and the 32 mm zone.
[`ApertureKeepOut`](@ref) puts each beam's worst-case clearance from its entrance hole
into an optimization objective. [`backlighter_foil`](@ref) and [`line_of_sight`](@ref)
place a foil on a diagnostic line of sight and check it against shields and beams, and
[`backlighter_signal`](@ref) follows the deposition through a He-α conversion efficiency
fitted to published foil data ([`kshell_ce`](@ref)) to the photons at the detector.

```@autodocs
Modules = [Briareus]
Pages = ["src/clearance.jl", "src/backlighter.jl", "src/backlighter_signal.jl", "src/keepout.jl", "src/unconverted.jl", "src/talis.jl"]
Private = false
```
