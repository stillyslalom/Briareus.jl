# Beams and facilities

Beams enter from facility ports: [`nif_quads`](@ref) and [`nif_beams`](@ref) (the LMF
beam table, with measured phase-plate fits), NIF bundles ([`nif_bundle`](@ref)), OMEGA,
or any port table ([`beams_from_table`](@ref)). Each beam has an elliptical
super-Gaussian spot, an f-number cone, pointing and defocus.

A beam's direction is fixed by its port (``\mathbf{k} = -\hat{\mathbf{p}}``), and
pointing translates its axis through `aim`. Best focus is at `aim + defocus·k`; positive
defocus moves the focus along the beam, as in VisRad and the NIF CMT. The spot axes are
given in the port's (``\hat\theta``, ``\hat\phi``) frame. [`VisRadProfile`](@ref)
reproduces VisRad's beam away from focus, by the law VisRad's own exported parameters
follow (research notes, "VisRad's beam profile law").

```@autodocs
Modules = [Briareus]
Pages = ["src/beams.jl", "src/facilities.jl"]
Private = false
```
