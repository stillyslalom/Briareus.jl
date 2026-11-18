# Capture and clearance rules as quantities of an evaluation, for objectives and
# constraints written by hand: each a pure function of the evaluation (or of its beams),
# generic in the number type and thread-safe, with a natural limit (a fraction at most a
# ceiling or at least a floor, a clearance at least a minimum).

"""
    direct_fraction(ev) -> Real
    direct_fraction(ev, names) -> Real

Fraction of the delivered laser power absorbed by the scenario's observed elements (the
capsule whose drive is measured), or with `names` (a collection of part or component
names) by those parts instead: the fraction whose square the `direct` penalty of
[`symmetry_objective`](@ref), [`asymmetry_objective`](@ref) and
[`face_flux_objective`](@ref) weighs (`direct_names` there is `names` here). A name that
is not in the target is an `ArgumentError`. Relative to the effective (post-LPI) power,
as every deposition fraction; see [`clipped_fraction`](@ref) and
[`missed_fraction`](@ref) for the other two penalized fractions.
"""
direct_fraction(ev::Evaluation, names = nothing) = _direct_fraction(ev, names)
