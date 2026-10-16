# Time-dependent drive: pulse shapes and quasi-static radiosity histories.

"""
    PulseShape(times, powers)

Piecewise-linear power history (ns, TW per beam), zero outside the table. Callable:
`pulse(t)`.
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

(s::PulseShape)(t) = (t < s.t[1] || t > s.t[end]) ? 0.0 : _interp(s.t, s.p, t)

"Laser energy (kJ = TW·ns) of one beam following the pulse."
pulse_energy(s::PulseShape) =
    sum((s.t[k+1] - s.t[k]) * (s.p[k] + s.p[k+1]) / 2 for k in 1:length(s.t)-1; init = 0.0)

"""
    with_pulses(beams, pulses, t) -> Vector{Beam}

Beams with power set to `pulses[key](t)`, where `key` is each beam's `group` (pulses
given as a `Dict` or NamedTuple keyed by group), or a single `PulseShape` for all.
"""
with_pulses(beams, pulse::PulseShape, t) = [remake(b; power = pulse(t)) for b in beams]
with_pulses(beams, pulses, t) = [remake(b; power = pulses[b.group](t)) for b in beams]

"""
    drive_history(model, op, beams_at, times; observe = :capsule, lmax = 6,
                  wall = (:wall,)) -> NamedTuple

Quasi-static drive history: at each time, deposit `beams_at(t)` with `op`, solve the
radiosity at that time (time- and temperature-dependent albedos included), and record

- `delivered` — laser power (TW)
- `flux`, `Tr` — mean incident flux on `observe` (TW/mm²) and its radiation
  temperature (eV)
- `moments` — `(lmax+1) × nt` Legendre moments of the observed flux (`aₗ/a₀` rows
  are `moments[l+1, :] ./ moments[1, :]`)
- `wall_Tr`, `wall_albedo` — area-weighted mean incident-flux temperature and albedo
  of the `wall` components.
"""
function drive_history(m::RadiationModel, op::DepositionOperator, beams_at, times;
                       observe = :capsule, lmax::Integer = 6, wall = (:wall,))
    t = m.target
    obs = elements(t, (observe isa Symbol ? (observe,) : observe)...)
    we = elements(t, wall...)
    Aw = t.elem_area[we]
    nt = length(times)
    delivered, flux, Tr = zeros(nt), zeros(nt), zeros(nt)
    wall_Tr, wall_albedo = zeros(nt), zeros(nt)
    moments = zeros(lmax + 1, nt)
    for (k, τ) in enumerate(times)
        beams = beams_at(τ)
        dep = deposit(op, beams)
        sol = solve(m, dep; time = τ)
        a = legendre_moments(sol.H, t, obs; lmax)
        delivered[k] = dep.delivered
        moments[:, k] = a
        flux[k] = a[1]
        Tr[k] = radiation_temperature(a[1])
        Hw = sum(sol.H[we] .* Aw) / sum(Aw)
        wall_Tr[k] = radiation_temperature(Hw)
        wall_albedo[k] = sum(sol.albedo[we] .* Aw) / sum(Aw)
    end
    return (; time = collect(Float64, times), delivered, flux, Tr, moments, wall_Tr, wall_albedo)
end
