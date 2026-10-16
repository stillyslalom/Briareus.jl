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
                  wall = (:wall,), lpi = nothing, holes = nothing, ylm = false) -> NamedTuple

Quasi-static drive history: at each time, deposit `beams_at(t)` with `op` (after the
laser–plasma model `lpi`, evaluated at that time), solve the
radiosity at that time (time- and temperature-dependent albedos included), and record

- `delivered` — laser power (TW)
- `flux`, `Tr` — mean incident flux on `observe` (TW/mm²) and its radiation
  temperature (eV)
- `moments` — `(lmax+1) × nt` Legendre moments of the observed flux (`aₗ/a₀` rows
  are `moments[l+1, :] ./ moments[1, :]`)
- `ylm` — with `ylm = true`, the `(lmax+1)² × nt` real spherical-harmonic
  coefficients of the observed flux (`ylm_weights`, in `mode_index` order; columns
  work with `mode_amplitude`), otherwise `nothing`
- `wall_Tr`, `wall_albedo` — area-weighted mean incident-flux temperature and albedo
  of the `wall` components
- `escaped` — power (TW) through each aperture, and `transmitted` through each
  thin-foil (`BurnthroughAlbedo`) component, as in `energy_balance`
- `window_loss` — power lost through diagnostic windows: `transmitted` plus `escaped`
  through every aperture except the entrance `holes` (default: the apertures not on a
  `CylinderSurface`, i.e. a hohlraum's LEHs)
- `front` — per thin-foil component, its area-mean heat-front areal mass as a fraction
  of the foil's (`≥ 1` once burned through), from its own flux history even when the
  model has an explicit `burnthrough`
- `burnthrough` — per thin-foil component, each element's burnthrough time (ns; `Inf`
  if it has not burned through).

Each thin-foil element without an explicit `burnthrough` time tracks its own heat
front: the drive integral of `heat_front` over that element's incident-flux temperature
history (trapezoidal between `times`, from `times[1]`), so its burnthrough is
consistent with the flux the windowed hohlraum puts on it. A step is solved with the
front extended from the previous step at constant temperature, and solved again
until the step's own temperatures no longer move a burnthrough that affects its
albedo. Resolve the
burnthrough with steps short against the transition `width`.
"""
function drive_history(m::RadiationModel, op::DepositionOperator, beams_at, times;
                       observe = :capsule, lmax::Integer = 6, wall = (:wall,),
                       lpi = nothing, holes = nothing, ylm::Bool = false)
    t = m.target
    obs = elements(t, (observe isa Symbol ? (observe,) : observe)...)
    we = elements(t, wall...)
    Aw = t.elem_area[we]
    nt = length(times)
    W = ylm ? ylm_weights(t, obs; lmax) : nothing
    clm = ylm ? zeros((lmax + 1)^2, nt) : nothing
    N = nelements(t)
    delivered, flux, Tr = zeros(nt), zeros(nt), zeros(nt)
    wall_Tr, wall_albedo = zeros(nt), zeros(nt)
    moments = zeros(lmax + 1, nt)
    # thin-foil components and the elements whose burnthrough is tracked
    foils = [ci for (ci, c) in enumerate(t.components)
             if c.kind === :surface && c.material.albedo isa BurnthroughAlbedo]
    foil_elems = [e for ci in foils for e in t.comp_elems[ci]]
    model(e) = element_material(t, e).albedo
    tracked = [e for e in foil_elems if model(e).burnthrough === nothing]
    need = zeros(N)
    for e in foil_elems
        need[e] = _front_integral(model(e).front, areal_mass(model(e)))
    end
    tb = fill(NaN, N)              # burnthrough times (NaN: not burned)
    for e in foil_elems
        model(e).burnthrough === nothing || (tb[e] = model(e).burnthrough)
    end
    I = zeros(N)                   # drive integral of each thin-foil element
    rate = zeros(N)                # its integrand at the previous time
    holes = holes === nothing ?
            [c.name for c in t.components if c.kind === :aperture && !(c.surface isa CylinderSurface)] :
            collect(Symbol, holes)
    apertures = [c.name for c in t.components if c.kind === :aperture]
    escaped = Dict(n => zeros(nt) for n in apertures)
    transmitted = Dict(t.components[ci].name => zeros(nt) for ci in foils)
    front = Dict(t.components[ci].name => zeros(nt) for ci in foils)
    window_loss = zeros(nt)
    for (k, τ) in enumerate(times)
        beams = apply_lpi(lpi, beams_at(τ), τ)
        dep = deposit(op, beams)
        Δ = k == 1 ? 0.0 : τ - times[k-1]
        # predictor: each front advanced through the step at its last temperature
        pred = copy(tb)
        for e in tracked
            isnan(tb[e]) && rate[e] > 0 && I[e] + Δ * rate[e] >= need[e] &&
                (pred[e] = times[k-1] + (need[e] - I[e]) / rate[e])
        end
        sol = solve(m, dep; time = τ, burnthrough = pred)
        # corrector: the trapezoidal integral with this step's temperatures
        function advance(sol)
            Ik, tk, rk = copy(I), copy(tb), copy(rate)
            for e in foil_elems
                rk[e] = _front_rate(model(e).front, radiation_temperature(ForwardDiff.value(sol.H[e])))
                k == 1 && continue
                dI = Δ * (rate[e] + rk[e]) / 2
                if isnan(tk[e]) && Ik[e] + dI >= need[e]   # (explicit times are never NaN)
                    tk[e] = times[k-1] + _linear_crossing(rate[e], rk[e], Δ, need[e] - Ik[e])
                end
                Ik[e] += dI
            end
            return Ik, tk, rk
        end
        Ik, tk, rk = advance(sol)
        # re-solve while the step's own temperatures move a burnthrough that its albedo
        # sees (a fixed point; the transition is smooth, so it converges in a few passes)
        used = pred
        for _ in 1:20
            moved(e) = !isequal(tk[e], used[e]) && min(_nan_inf(tk[e]), _nan_inf(used[e])) <= τ &&
                       !(abs(tk[e] - used[e]) <= 1e-12 * max(abs(τ), 1.0))
            any(moved, tracked) || break
            used = tk
            sol = solve(m, dep; time = τ, burnthrough = used)
            Ik, tk, rk = advance(sol)
        end
        I, tb, rate = Ik, tk, rk
        a = legendre_moments(sol.H, t, obs; lmax)
        delivered[k] = dep.delivered
        moments[:, k] = a
        ylm && (clm[:, k] = W * sol.H[obs])
        flux[k] = a[1]
        Tr[k] = radiation_temperature(a[1])
        Hw = sum(sol.H[we] .* Aw) / sum(Aw)
        wall_Tr[k] = radiation_temperature(Hw)
        wall_albedo[k] = sum(sol.albedo[we] .* Aw) / sum(Aw)
        bal = energy_balance(sol)
        for n in apertures
            escaped[n][k] = bal.escaped[n]
            n in holes || (window_loss[k] += bal.escaped[n])
        end
        for ci in foils
            c = t.components[ci]
            transmitted[c.name][k] = bal.transmitted[c.name]
            window_loss[k] += bal.transmitted[c.name]
            es = t.comp_elems[ci]
            front[c.name][k] = sum(_front_mass(model(e).front, I[e]) / areal_mass(model(e)) * t.elem_area[e]
                                   for e in es) / sum(t.elem_area[es])
        end
    end
    burnthrough = Dict(t.components[ci].name => [_nan_inf(tb[e]) for e in t.comp_elems[ci]] for ci in foils)
    return (; time = collect(Float64, times), delivered, flux, Tr, moments, ylm = clm, wall_Tr,
            wall_albedo, escaped, transmitted, window_loss, front, burnthrough)
end

_nan_inf(x) = isnan(x) ? Inf : x
