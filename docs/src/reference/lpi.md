# Laser–plasma models

Corrections applied between the beam setup and the deposition, differentiable in power
and Δλ: per-cone, time-dependent [`Backscatter`](@ref); an empirical
[`CalibratedTransfer`](@ref) (transfer as a function of Δλ); a reduced kinetic
[`CBET`](@ref) model at the LEHs (multi-species Maxwellian ion-acoustic response, flow
Doppler shift, polarization smoothing and a smooth δn/n saturation clamp), whose
multiplier, flow and clamp [`calibrate_cbet`](@ref) fits to transfer data;
[`TwoZoneCBET`](@ref), which separates the flowing plasma outside the LEH from the
stagnant one inside, with calibrations to published HYDRA transfer for HYBRID-E shots;
and inner-beam absorption in the gold bubble ([`BubbleAbsorption`](@ref)), whose picket
fluence ratio [`picket_ratio`](@ref) can take from a pointing's own outer picket fluence
on the wall ([`picket_fluence`](@ref)). [`PicketLinkedBubble`](@ref) takes it from the
beams inside every evaluation, so that a re-pointed design grows its own bubble, and
[`bubble_prediction`](@ref) reports what a re-pointing does to the picket, the bubble
and the P2 swing; [`picket_irradiance`](@ref) and [`picket_ceiling`](@ref) give the
irradiance during the picket and hold it in a constrained optimization.
[`lpi_deposit`](@ref) applies them together. The CBET model assumes uniform LEH plasma
and full pairwise overlap, and its unsaturated coupling overpredicts transfer: its
magnitudes hold only after calibration against shot data or reference transfer points.

```@autodocs
Modules = [Briareus]
Pages = ["src/lpi.jl", "src/bubble.jl", "src/picket_fluence.jl", "src/cbet.jl"]
Private = false
```
