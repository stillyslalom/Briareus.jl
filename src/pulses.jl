# This file defines the pulse shape of a beam group and the beams of a pulse at a time.

"""
    PulseShape(times, powers)

Piecewise-linear power history (ns, TW per beam), zero outside the table. Callable:
`pulse(t)`. Two pulse shapes are equal (`==`, `isequal`, `hash`) when their times and
powers are.
"""
struct PulseShape
    t::Vector{Float64}
    p::Vector{Float64}
    function PulseShape(t, p)
        length(t) == length(p) || throw(ArgumentError("times and powers differ in length"))
        issorted(t) || throw(ArgumentError("times must be sorted"))
        new(collect(Float64, t), collect(Float64, p))
    end
end

Base.:(==)(a::PulseShape, b::PulseShape) = a.t == b.t && a.p == b.p
Base.isequal(a::PulseShape, b::PulseShape) = isequal(a.t, b.t) && isequal(a.p, b.p)
Base.hash(s::PulseShape, h::UInt) = hash(s.p, hash(s.t, hash(:PulseShape, h)))

(s::PulseShape)(t) = (t < s.t[1] || t > s.t[end]) ? 0.0 : _interp(s.t, s.p, t)

# The trapezoid rule for ∫ y dt through the samples (t, y), of which there are at least two.
_trapz(t, y) = sum((t[k+1] - t[k]) * (y[k] + y[k+1]) / 2 for k in 1:length(t)-1)

"Laser energy (kJ = TW·ns) of one beam following the pulse."
pulse_energy(s::PulseShape) = length(s.t) < 2 ? 0.0 : _trapz(s.t, s.p)

"""
    with_pulses(beams, pulses, t; by = :group) -> Vector{Beam}

Beams with power set to `pulses[key](t)`. `key` is the `group` of each beam (pulses
given as a `Dict` or NamedTuple keyed by group), its `quad` as a `Symbol` (`by = :quad`),
its `name` as a `Symbol` (`by = :beam`) or `by(beam)` for a function. `pulses` may also
be a single `PulseShape` for all beams. Every pulse reader and writer keys by `Symbol`:
[`read_pulse_table`](@ref) by column header, [`cone_fraction_pulses`](@ref) by group and
[`read_visrad_powers`](@ref) by beam name (`by = :beam`).
"""
function with_pulses(beams, pulses, t; by = :group)
    key = _beam_key(by)
    return [remake(b; power = _beam_pulse(pulses, key, b)(t)) for b in beams]
end

_beam_key(by) = by === :group ? (b -> b.group) : by === :quad ? (b -> Symbol(b.quad)) :
                by === :beam ? (b -> Symbol(b.name)) : by

# The pulse of beam `b`: one `PulseShape` for all, or `pulses[key(b)]`.
_beam_pulse(p::PulseShape, key, b) = p
function _beam_pulse(p, key, b)
    k = key(b)
    haskey(p, k) || throw(KeyError(k))
    return p[k]
end
