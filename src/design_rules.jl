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

"""
    clipped_fraction(ev, names...) -> Real

Fraction of the delivered laser power striking the outside (back side) of the named
parts or components (a part stands for all its pieces, see [`part_pieces`](@ref)); with
no names, of every component. Beside [`power_fraction`](@ref)`(ev, names...)`, the power
absorbed on their elements, the sum

    power_fraction(ev, names...) + clipped_fraction(ev, names...)

is all the laser light the parts receive, from either side: the quantity of a ceiling
on hardware meant to stay dark (a closed window, a shield, a stalk: at most the value
of an accepted design) or of a floor on capture by a part meant to be hit. Both are
totals over the beams, in the evaluation's number type (differentiable), and relative
to the effective (post-LPI) power, as [`component_power`](@ref) reports the same powers
in TW (`name` and `Symbol(name, "_outside")`). An open window is an aperture: light
passes through it, so it receives none.
"""
clipped_fraction(ev::Evaluation, names::Symbol...) =
    sum(ev.deposition.backside[ci] for n in names for ci in part_pieces(ev.scenario.target, n);
        init = zero(ev.deposition.delivered)) / ev.deposition.delivered
