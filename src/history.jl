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
    with_pulses(beams, pulses, t; by = :group) -> Vector{Beam}

Beams with power set to `pulses[key](t)`, where `key` is each beam's `group` (pulses
given as a `Dict` or NamedTuple keyed by group), its `quad` (`by = :quad`), its `name`
(`by = :beam`) or `by(beam)` for a function, or a single `PulseShape` for all.
"""
with_pulses(beams, pulse::PulseShape, t; by = :group) = [remake(b; power = pulse(t)) for b in beams]
function with_pulses(beams, pulses, t; by = :group)
    key = _beam_key(by)
    return [remake(b; power = pulses[key(b)](t)) for b in beams]
end

_beam_key(by) = by === :group ? (b -> b.group) : by === :quad ? (b -> Symbol(b.quad)) :
                by === :beam ? (b -> Symbol(b.name)) : by

"""
    drive_history(model, op, beams_at, times; observe = :capsule, lmax = 6,
                  wall = (:wall,), lpi = nothing, holes = nothing, ylm = false,
                  transition = nothing, flux_change = nothing, min_step = 1e-3,
                  tol = 1e-4, keep_albedo = false) -> NamedTuple

Quasi-static drive history: at each time, deposit `beams_at(t)` with `op` (after the
laser–plasma model `lpi`, evaluated at that time), solve the
radiosity at that time (time- and temperature-dependent albedos included), and record

- `time` — the times solved (ns): `times`, plus any steps inserted by `transition` or
  `flux_change`
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
  thin-foil (`BurnthroughAlbedo`) and closing-gap (`GapClosure`) component, as in
  `energy_balance` (by default the fraction `s T` of the flux on each burned-through
  element, and `1 - s` on a gap a fraction `s` closed)
- `window_loss` — power lost through diagnostic windows: `transmitted` plus `escaped`
  through every aperture except the entrance `holes` (default: the apertures not on a
  `CylinderSurface`, i.e. a hohlraum's LEHs)
- `front` — per thin-foil component, its area-mean heat-front areal mass as a fraction
  of the foil's (`≥ 1` once burned through), from its own flux history even when the
  model has an explicit `burnthrough`
- `burnthrough` — per thin-foil component, each element's burnthrough time (ns; `Inf`
  if it has not burned through)
- `albedo` — with `keep_albedo = true`, each step's converged albedo per element (the
  state a flux response linearizes about, e.g. for `retune_drive`), otherwise `nothing`
- `setup` — the power of `beams_at(t)` (TW, before `lpi`); `laser` — absorbed laser
  power per component (TW, back side included); `source` — x-ray source power per
  component (Σ ηP); `absorbed` — x-ray power absorbed per surface component (and per
  closing lip), as in `energy_balance`; `lost` — x-ray power lost to open space.
  `energy_budget` integrates them.

Each step's temperature-dependent albedos are iterated to `tol` (see `solve`), starting
from the previous step's temperatures.

Each thin-foil element without an explicit `burnthrough` time tracks its own heat
front: the drive integral of `heat_front` over that element's incident-flux temperature
history (trapezoidal between the steps, from `times[1]`), so its burnthrough is
consistent with the flux the windowed hohlraum puts on it. A step is solved with the
front extended from the previous step at constant temperature, and solved again
until the step's own temperatures no longer move a burnthrough that affects its
albedo.

A burnthrough transition takes the foil's `width` (0.4 ns by default), and a step
much longer than that jumps across it. Two optional criteria insert steps between
`times`:

- `transition` — no thin-foil element's transition (the smooth step `s` of
  `BurnthroughAlbedo`, 0 intact to 1 burned through), and no closing gap's (`GapClosure`,
  0 open to 1 closed), advances by more than this in one step.
  Steps are planned from each element's known or predicted burnthrough time, and
  halved if the step's own temperatures move a burnthrough so that the advance
  exceeds 1.5 `transition`.
- `flux_change` — the observed flux relative to its mean, `H/⟨H⟩` on each element of
  `observe`, changes by no more than this in one step; a step that exceeds it is
  halved and solved again.

Neither shortens a step below `min_step` (ns). Without them (the default) the steps
are exactly `times`.
"""
drive_history(m::RadiationModel, op::DepositionOperator, beams_at, times; kwargs...) =
    _drive_history(m.target, τ -> m, τ -> op, beams_at, times; kwargs...)

# The history on the models `model_at(τ)` and deposition operators `operator_at(τ)`, whose
# elements are those of `t` (a fixed target, or the keyframes of a moving one).
function _drive_history(t::Target, model_at, operator_at, beams_at, times;
                        observe = :capsule, lmax::Integer = 6, wall = (:wall,),
                        lpi = nothing, holes = nothing, ylm::Bool = false,
                        transition = nothing, flux_change = nothing, min_step::Real = 1e-3,
                        tol = 1e-4, keep_albedo::Bool = false)
    transition === nothing || transition > 0 || throw(ArgumentError("transition must be positive"))
    flux_change === nothing || flux_change > 0 || throw(ArgumentError("flux_change must be positive"))
    min_step > 0 || throw(ArgumentError("min_step must be positive"))
    issorted(times) || throw(ArgumentError("times must be sorted"))
    obs = elements(t, (observe isa Symbol ? (observe,) : observe)...)
    we = elements(t, wall...)
    weights = _TargetCache(tg -> ylm_weights(tg, obs; lmax))
    N = nelements(t)
    # thin-foil components and the elements whose burnthrough is tracked
    foils = [ci for (ci, c) in enumerate(t.components)
             if c.kind === :surface && c.material.albedo isa BurnthroughAlbedo]
    foil_elems = [e for ci in foils for e in t.comp_elems[ci]]
    # components that transmit part of their flux out: thin foils and closing gaps
    transmitters = [ci for (ci, c) in enumerate(t.components)
                    if c.kind === :surface && _transmits(c.material.albedo)]
    model(e) = element_material(t, e).albedo
    tracked = [e for e in foil_elems if model(e).burnthrough === nothing]
    # transitions a step can resolve (closing entrance-hole lips included)
    smooth = [e for (ci, c) in enumerate(t.components) if _transmits(c.material.albedo)
              for e in t.comp_elems[ci] if model(e).width > 0]
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
    # records, one entry per step taken
    time = Float64[]
    delivered, flux, Tr = Float64[], Float64[], Float64[]
    wall_Tr, wall_albedo = Float64[], Float64[]
    moments = Vector{Float64}[]
    clm = Vector{Float64}[]
    escaped = Dict(n => Float64[] for n in apertures)
    transmitted = Dict(t.components[ci].name => Float64[] for ci in transmitters)
    front = Dict(t.components[ci].name => Float64[] for ci in foils)
    window_loss = Float64[]
    albedo_steps = Vector{Float64}[]
    names = [c.name for c in t.components]
    setup, lost = Float64[], Float64[]
    laser = Dict(n => Float64[] for n in names)
    source = Dict(n => Float64[] for n in names)
    absorbed = Dict(c.name => Float64[] for c in t.components
                    if c.kind !== :aperture || _transmits(c.material.albedo))

    # The step from τ0 (`nothing` at the first time) to τ: its solution and the
    # thin-foil state at τ. The state at τ0 is left untouched, so a step can be retried.
    function step(τ0, τ, dep, guess)
        Δ = τ0 === nothing ? 0.0 : τ - τ0
        # predictor: each front advanced through the step at its last temperature
        pred = copy(tb)
        for e in tracked
            isnan(tb[e]) && rate[e] > 0 && I[e] + Δ * rate[e] >= need[e] &&
                (pred[e] = τ0 + (need[e] - I[e]) / rate[e])
        end
        mτ = model_at(τ)
        sol = solve(mτ, dep; time = τ, burnthrough = pred, guess, tol)
        # corrector: the trapezoidal integral with this step's temperatures
        function advance(sol)
            Ik, tk, rk = copy(I), copy(tb), copy(rate)
            for e in foil_elems
                rk[e] = _front_rate(model(e).front, radiation_temperature(ForwardDiff.value(sol.H[e])))
                τ0 === nothing && continue
                dI = Δ * (rate[e] + rk[e]) / 2
                if isnan(tk[e]) && Ik[e] + dI >= need[e]   # (explicit times are never NaN)
                    tk[e] = τ0 + _linear_crossing(rate[e], rk[e], Δ, need[e] - Ik[e])
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
            sol = solve(mτ, dep; time = τ, burnthrough = used, guess = sol, tol)
            Ik, tk, rk = advance(sol)
        end
        return sol, Ik, tk, rk
    end
    # the largest advance of a smooth transition from τ0 to τ, given burnthrough times
    advance_of(τ0, τ, tbs) =
        maximum((_transition_step(model(e), τ, _nan_inf(tbs[e])) -
                 _transition_step(model(e), τ0, _nan_inf(tbs[e])) for e in smooth); init = 0.0)
    shape(sol, a0) = ForwardDiff.value.(sol.H[obs]) ./ ForwardDiff.value(a0)

    τ0 = nothing                   # time of the last step taken
    last = nothing                 # its solution, which starts the next step's iteration
    shape0 = nothing               # its observed flux shape
    for τgoal in times
        while true
            τ = τgoal
            if τ0 !== nothing && transition !== nothing && !isempty(smooth)
                # plan the step from the burnthrough times known or predicted at τ0
                tc = copy(tb)
                for e in tracked
                    isnan(tb[e]) && rate[e] > 0 && (tc[e] = τ0 + (need[e] - I[e]) / rate[e])
                end
                if advance_of(τ0, τ, tc) > transition
                    lo, hi = τ0 + min_step, τ       # the longest step within the tolerance
                    if advance_of(τ0, lo, tc) <= transition
                        while hi - lo > 1e-3 * min_step
                            mid = (lo + hi) / 2
                            advance_of(τ0, mid, tc) > transition ? (hi = mid) : (lo = mid)
                        end
                    end
                    τ = min(lo, τgoal)
                end
            end
            local sol, Ik, tk, rk, a, dep, P0
            while true
                b0 = beams_at(τ)
                P0 = sum(b -> ForwardDiff.value(b.power), b0; init = 0.0)
                dep = lpi_deposit(lpi, operator_at(τ), apply_lpi(lpi, b0, τ), τ)
                sol, Ik, tk, rk = step(τ0, τ, dep, last)
                a = legendre_moments(sol.H, sol.model.target, obs; lmax)
                τ0 === nothing && break
                h = τ - τ0
                h / 2 >= min_step || break
                too_far = (transition !== nothing && advance_of(τ0, τ, tk) > 1.5 * transition) ||
                          (flux_change !== nothing &&
                           maximum(abs.(shape(sol, a[1]) .- shape0)) > flux_change)
                too_far || break
                τ = τ0 + h / 2
            end
            I, tb, rate = Ik, tk, rk
            τ0, last = τ, sol
            shape0 = shape(sol, a[1])
            push!(time, τ)
            push!(delivered, dep.delivered)
            push!(moments, a)
            ylm && push!(clm, weights(sol.model.target) * sol.H[obs])
            Aw = sol.model.target.elem_area[we]
            push!(flux, a[1])
            push!(Tr, radiation_temperature(a[1]))
            Hw = sum(sol.H[we] .* Aw) / sum(Aw)
            push!(wall_Tr, radiation_temperature(Hw))
            push!(wall_albedo, sum(sol.albedo[we] .* Aw) / sum(Aw))
            keep_albedo && push!(albedo_steps, copy(sol.albedo))
            bal = energy_balance(sol)
            loss = 0.0
            for n in apertures
                push!(escaped[n], bal.escaped[n])
                n in holes || (loss += bal.escaped[n])
            end
            for ci in transmitters
                c = t.components[ci]
                push!(transmitted[c.name], bal.transmitted[c.name])
                loss += bal.transmitted[c.name]
            end
            for ci in foils
                c = t.components[ci]
                es = t.comp_elems[ci]
                push!(front[c.name], sum(_front_mass(model(e).front, I[e]) / areal_mass(model(e)) * t.elem_area[e]
                                         for e in es) / sum(t.elem_area[es]))
            end
            push!(window_loss, loss)
            push!(setup, P0)
            push!(lost, bal.lost)
            tg = sol.model.target
            v = ForwardDiff.value
            for (ci, c) in enumerate(tg.components)
                es = tg.comp_elems[ci]
                push!(laser[c.name], sum(v(sol.laser[e]) for e in es; init = 0.0) + v(dep.backside[ci]))
                push!(source[c.name], sum(v(sol.source[e]) * tg.elem_area[e] for e in es; init = 0.0))
                haskey(absorbed, c.name) && push!(absorbed[c.name], bal.absorbed[c.name])
            end
            τ == τgoal && break
        end
    end
    burnthrough = Dict(t.components[ci].name => [_nan_inf(tb[e]) for e in t.comp_elems[ci]] for ci in foils)
    return (; time, delivered, flux, Tr, moments = _columns(moments, lmax + 1),
            ylm = ylm ? _columns(clm, (lmax + 1)^2) : nothing, wall_Tr,
            wall_albedo, escaped, transmitted, window_loss, front, burnthrough,
            albedo = keep_albedo ? albedo_steps : nothing, setup, laser, source, absorbed, lost)
end

"""
    energy_budget(history; from = nothing, to = nothing) -> NamedTuple

Energies (kJ = TW·ns) of a `drive_history`, the trapezoidal integrals of its power
records over the solved times (or over `[from, to]`, with the powers held at zero
outside the solved times only through the trapezoids they bound):

- `setup` — the laser power `beams_at(t)` delivers (before the laser–plasma model);
  `coupled` — after it (`delivered`); `backscattered` — their difference (what
  `Backscatter` and other lossy models remove; CBET conserves power);
- `laser` — absorbed laser energy per component (front elements plus back side),
  `unabsorbed` — `coupled − Σ laser`, the light that leaves through the holes;
- `source` — x-ray energy converted per component (Σ η P); `absorbed` — x-ray energy
  absorbed per surface component (`(1 − α) H A`; a `CapsuleAlbedo` capsule's at its
  effective radius: E_cap); `escaped` per aperture; `transmitted` per thin-foil or
  closing-gap component; `lost` to open space;
- `heating` — per surface component, what stays in it: unconverted laser plus absorbed
  x-rays, `laser − source + absorbed`.

`setup ≈ backscattered + unabsorbed + Σ heating + Σ escaped + Σ transmitted + lost`, to
the closure of the view factors.
"""
function energy_budget(h; from = nothing, to = nothing)
    t = h.time
    keep = [(from === nothing || τ >= from) && (to === nothing || τ <= to) for τ in t]
    tt = t[keep]
    ∫(p) = (q = p[keep]; sum((tt[k+1] - tt[k]) * (q[k] + q[k+1]) / 2 for k in 1:length(tt)-1; init = 0.0))
    d(x) = Dict(n => ∫(v) for (n, v) in x)
    setup, coupled = ∫(h.setup), ∫(h.delivered)
    laser, source, absorbed = d(h.laser), d(h.source), d(h.absorbed)
    heating = Dict(n => laser[n] - source[n] + absorbed[n] for n in keys(absorbed))
    return (; setup, coupled, backscattered = setup - coupled, laser,
            unabsorbed = coupled - sum(values(laser); init = 0.0), source, absorbed,
            escaped = d(h.escaped), transmitted = d(h.transmitted), lost = ∫(h.lost), heating)
end

_columns(v, n) = isempty(v) ? zeros(n, 0) : reduce(hcat, v)
_nan_inf(x) = isnan(x) ? Inf : x

# A value computed from a target and kept while the target stays the same object.
mutable struct _TargetCache{F}
    f::F
    target::Union{Nothing,Target}
    value::Any
end
_TargetCache(f) = _TargetCache(f, nothing, nothing)
function (c::_TargetCache)(t::Target)
    c.target === t || (c.value = c.f(t); c.target = t)
    return c.value
end
