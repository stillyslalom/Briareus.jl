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
                  tol = 1e-4, keep_albedo = false, record = nothing,
                  wall_history = true, wall_release = true, mband = MBandSplit(),
                  spectral = nothing) -> NamedTuple

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
  state a flux response linearizes about, e.g. for `retune_drive`; where a wall re-emits
  stored energy, the effective albedo that carries its emission, so that `B = ηP/A + αH`),
  otherwise `nothing`
- `setup` — the power of `beams_at(t)` (TW, before `lpi`); `laser` — absorbed laser
  power per component (TW, back side included); `source` — x-ray source power per
  component (Σ ηP); `absorbed` — x-ray power absorbed per surface component (and per
  closing lip), as in `energy_balance` (net of `released`, the power wall components
  re-emit from stored energy; see below); `lost` — x-ray power lost to open space.
  `energy_budget` integrates them. Like `energy_balance`, `escaped`, `transmitted`,
  `laser`, `source` and `absorbed` are keyed by part, the pieces of a part split by
  `curved_components` summed into it; `front` and `burnthrough` stay per component.
- `lip_radius` — per computed entrance-hole lip (`LipClosure`, `hohlraum(; leh_closure =
  :computed)`), the radius of its plasma edge at each step (mm)
- `record` — with `record = f`, the vector of `f(sol)` for each step's
  `RadiationSolution` (e.g. `sol -> dante(dante1(), sol)` for a synthetic Dante, or a
  `spectral_solve`), otherwise `nothing`.
- `mband_moments`, `mband` — the `(lmax+1) × nt` Legendre moments of the observed flux
  at photon energies above the M-band's lower edge (`band_flux(s, mband.emission.band[1])`
  of the step's spectral solution `s`: laser-produced and recirculating M-band and the
  thermal tail) and the laser-produced M-band power per component (TW, keyed by part like
  `laser`); `nothing` with `mband = nothing`, or when the target has no surface of the
  model's materials or (per-group model) no beam of `beams_at(times[1])` is in its
  groups.

**M-band.** `mband` is an [`MBandSplit`](@ref): the M-band emission model and the
spectral split it is solved with. By default (`MBandSplit()`: `nif_mband()`, Dewald
2020's per-cone conversions on gold and Au-lined uranium, in the calibrated two-group
split) each step's laser-produced M-band comes from that step's own effective beams, as
`mband_power(mband.emission, op, beams; lpi, time)` gives it with the beams after `lpi`
(CBET, backscatter, the bubble): each beam group's (cone's) share of every element,
split off in the step's own deposition pass, so the M-band follows the cone balance in
time (an intensity-form `MBandSplit(MBandEmission())` converts the step's deposition
itself). The step's gray solution is then split by `spectral_solve` with the split's
`edges`, `thermalize` and group albedos; `spectral` holds `spectral_solve`'s other
keywords (`spot`, `tol`, ...). A diagnostic: the gray history (every other record) is
the same bit for bit with or without it. It costs ~6 % of a small history (600
elements) and less than the run-to-run spread on N221204 (research notes, "Per-cone
M-band in the drive history (2026-10-27)"); `mband = nothing` turns it off.

Each step's temperature-dependent albedos are iterated to `tol` (see `solve`), starting
from the previous step's temperatures.

**Wall-loss history.** With `wall_history = true` (the default) each element whose
albedo is a [`WallLossAlbedo`](@ref) loses energy as its own temperature history
requires, not as if it had been at today's temperature since `t = 0`. Lindl's
constant-T wall energy `E_w ∝ T^a t^b` (a = 4 + p = 3.3, b = 1 + q = 0.62 for the
defaults) is carried to a varying `T(t)` as `E_w ∝ [∫ (T/100 eV)^(a/b) dt]^b`, the form
of [`HeatFront`](@ref)'s drive integral (exact at constant T); its rate is the
constant-T loss at the effective time

    t_eff = ∫₀ᵗ (T(t′)/T(t))^(a/b) dt′,

integrated trapezoidally between the steps from each element's incident-flux
temperature (constant T before `times[1]`). A wall heated in a low foot holds little
energy at the peak temperature, so its loss early in the peak is several times the
constant-T value at `t` (`t_eff` ≈ the time since the rise), and the drive rises
through a flat peak. `wall_history = false` evaluates the albedo at `t` (the behavior
before 2026-10-24). Research notes, "N221204 late-time drive (2026-10-24)".

**Wall re-emission.** The self-similar loss never gives energy back: when `T` falls,
`t_eff` grows and the loss only shrinks. With `wall_release = true` (the default, with
`wall_history`) a falling drive makes the wall re-emit stored energy. The loss above is
taken as the heating response to the element's highest temperature so far,
`T_h = max T`, and a drop below it, the deficit `D = φ(T_h) − φ(T)` in the wall-energy
potential `φ = (T/100 eV)^a`, is answered by the same self-similar kernel, superposed
(Duhamel's principle for the wall's diffusion, linear in `φ`, *engineering*):

    F = L(T_h, t_eff) − f_r b C ∫ (t − t′)^(b−1) dD(t′),    C = L t_eff^(1−b) / (b φ(T_h)),

the net flux into the wall (`C`, the loss's own coefficient, is constant for gold's
Lindl loss; `f_r` is `WallLossAlbedo`'s `release`, 1.0, calibrated against 1-D
radiation hydrodynamics of an ablating gold wall, `examples/late_time_wall.jl`). A constant-T wall loses `C b φ t^(b−1)`; a drop `D` at `t_d` gives back
`f_r C D (t − t_d)^b` by `t`, quickly at first (the hot surface layer) and then as the
heat wave's own `t^b`. The effective albedo `1 − F/H` exceeds 1 where the wall gives
back more than it takes; `solve` carries the release as an emission added to the
element's source. While no element's temperature falls below its running maximum the
history is the heating loss's bit for bit. `released` records the re-emitted power
per component (`energy_budget` integrates it); `absorbed` is net of it. Research notes,
"Wall re-emission and the capsule radius (2026-10-24)" (the closure) and "Capsule
absorption, conversion and the ablating-wall release (2026-10-25)" (`release` 1.0).

**Computed lip closure.** Each entrance-hole lip with a [`LipClosure`](@ref) albedo
carries its closure depth `d` (open at `times[1]`), advanced over each step from the
trapezoidal mean of its free speed `v(T)` (`T` the lip's area-mean incident-flux
temperature, at the step's start and, corrected, at its end) and of the laser's hold
(the power of the beams in `beams_at`, after `lpi`, whose ports lie on the hole's side,
over the hole's initial area), exactly for the hold's relaxation. The step is solved
with each lip ring closed outside the edge `r₀ − d`, and again until the corrected depth
moves by less than 1e-5 mm.

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
# Compiled and inferred once for every closure, laser–plasma model and recorder: each
# call site's `beams_at` closure used to recompile the whole march (~2 s). Their values
# reach the solves by dynamic dispatch, at no measurable run-time cost, and the march
# is in the package's precompile workload (research notes, "Compile latency").
Base.@nospecializeinfer function _drive_history(
                        t::Target, @nospecialize(model_at), @nospecialize(operator_at),
                        @nospecialize(beams_at), @nospecialize(times);
                        observe = :capsule, lmax::Integer = 6, wall = (:wall,),
                        @nospecialize(lpi = nothing), holes = nothing, ylm::Bool = false,
                        transition = nothing, flux_change = nothing, min_step::Real = 1e-3,
                        tol = 1e-4, keep_albedo::Bool = false, @nospecialize(record = nothing),
                        wall_history::Bool = true, wall_release::Bool = true,
                        mband = MBandSplit(), @nospecialize(spectral = nothing))
    # the M-band model of each step (`nothing`: off) and its `spectral_solve` keywords
    emission = nothing
    if mband === nothing
        spectral === nothing ||
            throw(ArgumentError("`spectral` sets the M-band's spectral_solve; with mband = nothing it has no effect"))
    else
        mband isa MBandSplit ||
            throw(ArgumentError("mband must be an MBandSplit (e.g. MBandSplit(nif_mband())) or nothing"))
        spectral = spectral === nothing ? (;) : NamedTuple(spectral)
        split_keys = filter(in((:edges, :thermalize, :albedo, :mband, :mband_power)), collect(keys(spectral)))
        isempty(split_keys) || throw(ArgumentError("`spectral` takes spectral_solve's other keywords; " *
                                                   join(split_keys, ", ") * " are the MBandSplit's (`mband`)"))
        albedo = mband.albedo isa Float64 ?
                 mband_albedo(t, mband.emission, mband.edges; albedo = mband.albedo) : mband.albedo
        spectral = (; spectral..., edges = mband.edges, thermalize = mband.thermalize, albedo)
        # on only where there is M-band to make: a target with an emitting surface, lit by
        # a beam group the model converts (any beam, for the intensity form)
        em = mband.emission
        emits = any(e -> _emits_mband(em, t, e), 1:nelements(t))
        lit = em.groups === nothing || isempty(times) ||
              any(b -> haskey(em.groups, b.group), beams_at(first(times)))
        emits && lit && (emission = em)
    end
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
    smooth = [e for (ci, c) in enumerate(t.components)
              if _transmits(c.material.albedo) && !(c.material.albedo isa LipClosure)
              for e in t.comp_elems[ci] if model(e).width > 0]
    # entrance-hole lips whose closure is computed (`LipClosure`): their geometry, and
    # the closure depth, lip temperature and hole intensity at the last step taken
    lips = _lip_geometry(t)
    ldepth, lT, lI = zeros(length(lips)), zeros(length(lips)), zeros(length(lips))
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
    # wall-loss elements whose loss follows their own temperature history
    lossy = wall_history ?
            [e for e in 1:N if !is_aperture(t, e) && model(e) isa WallLossAlbedo] : Int[]
    Iw, ratew = zeros(N), zeros(N) # their wall-energy integral and integrand
    # re-emission (`wall_release`): each element's highest temperature, the deficit
    # φ(T_max) − φ(T) at every step taken, and the times of those steps
    Tmax, Dw = zeros(N), zeros(N)
    deficits, dtimes = Vector{Float64}[], Float64[]
    cooled = falses(N)             # elements whose deficit has ever been nonzero
    holes = holes === nothing ?
            [c.name for c in t.components if c.kind === :aperture && !(c.surface isa CylinderSurface)] :
            collect(Symbol, holes)
    apertures = [c.name for c in t.components if c.kind === :aperture]
    partof = Dict(c.name => c.part for c in t.components)
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
    recorded = Any[]
    names = [c.name for c in t.components]
    setup, lost = Float64[], Float64[]
    laser = Dict(n => Float64[] for n in names)
    source = Dict(n => Float64[] for n in names)
    released = Dict(n => Float64[] for n in names)
    absorbed = Dict(c.name => Float64[] for c in t.components
                    if c.kind !== :aperture || _transmits(c.material.albedo))
    lip_radius = Dict(L.name => Float64[] for L in lips)
    mband_moments = Vector{Float64}[]
    mband_laser = Dict(n => Float64[] for n in names)

    # The step from τ0 (`nothing` at the first time) to τ: its solution and the
    # thin-foil state at τ. The state at τ0 is left untouched, so a step can be retried.
    # `Il`: the laser intensity through each computed lip's hole at τ.
    function step(τ0, τ, dep, guess, Il)
        Δ = τ0 === nothing ? 0.0 : τ - τ0
        # predictor: each front advanced through the step at its last temperature
        pred = copy(tb)
        for e in tracked
            isnan(tb[e]) && rate[e] > 0 && I[e] + Δ * rate[e] >= need[e] &&
                (pred[e] = τ0 + (need[e] - I[e]) / rate[e])
        end
        mτ = model_at(τ)
        # the per-element state the albedos see: burnthrough times, and the wall-loss
        # history of the `lossy` elements (`_WallHistory`)
        hist = Dict(e => (τ0 === nothing ? _WallHistory(0.0, τ) :
                          !wall_release ? _WallHistory(Iw[e] + Δ * ratew[e] / 2, Δ / 2) :
                          _WallHistory(Iw[e] + Δ * ratew[e] / 2, Δ / 2, Tmax[e],
                                       cooled[e] ? _deficit_sum(model(e), deficits, dtimes, e, τ) : 0.0,
                                       Dw[e], Δ)) for e in lossy)
        # the lips' closure depths at τ, from the lip temperatures `T1` at τ
        depths(T1) = [τ0 === nothing ? 0.0 :
                      clamp(_lip_advance(L.model, ldepth[k], Δ,
                                         (_lip_speed(L.model, lT[k]) + _lip_speed(L.model, T1[k])) / 2,
                                         (_lip_hold(L.model, lI[k]) + _lip_hold(L.model, Il[k])) / 2),
                            0.0, L.r0 - L.rin) for (k, L) in enumerate(lips)]
        ds = depths(lT)            # predictor: the last step's lip temperatures
        function state(bt)
            isempty(lossy) && isempty(lips) && return bt
            s = Vector{Any}(bt)
            for e in lossy
                s[e] = hist[e]
            end
            for (k, L) in enumerate(lips), (j, e) in enumerate(L.elems)
                s[e] = _ring_closed(L.ra[j], L.rb[j], L.r0 - ds[k])
            end
            return s
        end
        sol = solve(mτ, dep; time = τ, burnthrough = state(pred), guess, tol)
        # corrector: the lips' depths with this step's own lip temperatures
        for _ in 1:5
            (isempty(lips) || τ0 === nothing) && break
            dn = depths(_lip_temperatures(lips, sol))
            maximum(abs.(dn .- ds)) <= 1e-5 && break
            ds = dn
            sol = solve(mτ, dep; time = τ, burnthrough = state(pred), guess = sol, tol)
        end
        # corrector: the trapezoidal integral with this step's temperatures
        function advance(sol)
            # `local`: the step's own `Ik, tk, rk` below would otherwise be captured
            # (boxed, and possibly undefined here; JET's report_call)
            local Ik, tk, rk = copy(I), copy(tb), copy(rate)
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
            sol = solve(mτ, dep; time = τ, burnthrough = state(used), guess = sol, tol)
            Ik, tk, rk = advance(sol)
        end
        # the wall-energy integrals at τ, with this step's temperatures (the heating
        # response at the highest temperature so far when the release is on)
        Iwk, rwk = copy(Iw), copy(ratew)
        for e in lossy
            T = radiation_temperature(ForwardDiff.value(sol.H[e]))
            wall_release && (T = max(T, Tmax[e]))
            rwk[e] = _loss_rate(model(e), T)
            Iwk[e] = τ0 === nothing ? τ * rwk[e] : hist[e].A + Δ * rwk[e] / 2
        end
        return sol, Ik, tk, rk, Iwk, rwk, hist, ds
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
            local sol, Ik, tk, rk, Iwk, rwk, hist, ds, Il, a, dep, P0, bl, msplit
            while true
                b0 = beams_at(τ)
                P0 = sum(b -> ForwardDiff.value(b.power), b0; init = 0.0)
                bl = apply_lpi(lpi, b0, τ)
                # per-group M-band: the cones split in the same deposition pass
                msplit = emission === nothing || emission.groups === nothing ? nothing :
                         _mband_split(emission, t, bl)
                dep = _lpi_deposit(lpi, operator_at(τ), bl, τ, msplit)
                Il = _lip_intensities(lips, bl)
                sol, Ik, tk, rk, Iwk, rwk, hist, ds = step(τ0, τ, dep, last, Il)
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
            Iw, ratew = Iwk, rwk
            if !isempty(lips)
                ldepth .= ds
                lT .= _lip_temperatures(lips, sol)
                lI .= Il
                for (k, L) in enumerate(lips)
                    push!(lip_radius[L.name], L.r0 - ds[k])
                end
            end
            # the release: power re-emitted from stored energy per component, then each
            # element's running maximum and deficit at τ
            rel = zeros(N)
            if wall_release
                for e in lossy
                    T = radiation_temperature(ForwardDiff.value(sol.H[e]))
                    rel[e] = _wall_split(model(e), hist[e], T)[3] * sol.model.target.elem_area[e]
                    Th = max(T, Tmax[e])
                    Dw[e] = _wall_potential(model(e), Th) - _wall_potential(model(e), T)
                    Tmax[e] = Th
                    Dw[e] != 0 && (cooled[e] = true)
                end
                push!(deficits, copy(Dw))
                push!(dtimes, τ)
            end
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
            # the emission a re-emitting wall adds to the solution's source (TW/mm²); the
            # wall albedo is the effective one, reflected over incident flux
            xe = wall_release && !isempty(lossy) ?
                 ForwardDiff.value.(sol.source) .- _xray_source(sol.model.target, ForwardDiff.value.(sol.laser), τ) :
                 nothing
            αw = xe === nothing || all(iszero, xe) ? sol.albedo[we] :
                 sol.albedo[we] .+ xe[we] ./ ForwardDiff.value.(sol.H[we])
            push!(wall_albedo, sum(αw .* Aw) / sum(Aw))
            keep_albedo && push!(albedo_steps, xe === nothing || all(iszero, xe) ? copy(sol.albedo) :
                                   sol.albedo .+ xe ./ max.(ForwardDiff.value.(sol.H), eps()))
            record === nothing || push!(recorded, record(sol))
            if emission !== nothing
                # the step's own effective beams, each cone's share of the deposition
                # (its beams alone, = `mband_power(emission, op, bl; lpi, time)`); the
                # intensity form converts the step's deposition itself
                tm = sol.model.target
                PM = msplit === nothing ? mband_power(emission, tm, dep.element) :
                     _mband_from_split(emission, tm, msplit)
                sm = spectral_solve(sol; spectral..., mband = emission, mband_power = PM)
                push!(mband_moments, legendre_moments(band_flux(sm, emission.band[1]), tm, obs; lmax))
                for (ci, c) in enumerate(tm.components)
                    push!(mband_laser[c.name], sum(sm.mband[e] * tm.elem_area[e] for e in tm.comp_elems[ci]; init = 0.0))
                end
            end
            bal =energy_balance(sol; by = :component)
            loss = 0.0
            for n in apertures
                push!(escaped[n], bal.escaped[n])
                n in holes || partof[n] in holes || (loss += bal.escaped[n])
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
            # the wall's re-emission rides on the solution's source: book the x-ray
            # source as converted laser (ηP) and the absorbed power net of the emission
            emitted = xe === nothing ? nothing : xe .* tg.elem_area
            for (ci, c) in enumerate(tg.components)
                es = tg.comp_elems[ci]
                Ec = emitted === nothing ? 0.0 : sum(emitted[e] for e in es; init = 0.0)
                push!(laser[c.name], sum(v(sol.laser[e]) for e in es; init = 0.0) + v(dep.backside[ci]))
                push!(source[c.name], sum(v(sol.source[e]) * tg.elem_area[e] for e in es; init = 0.0) - Ec)
                push!(released[c.name], sum(rel[e] for e in es; init = 0.0))
                haskey(absorbed, c.name) && push!(absorbed[c.name], bal.absorbed[c.name] - Ec)
            end
            τ == τgoal && break
        end
    end
    burnthrough = Dict(t.components[ci].name => [_nan_inf(tb[e]) for e in t.comp_elems[ci]] for ci in foils)
    if any(c -> c.part !== c.name, t.components)       # pieces of split parts summed into them
        escaped, transmitted, laser, source, absorbed, released =
            _by_part(t, escaped), _by_part(t, transmitted), _by_part(t, laser), _by_part(t, source),
            _by_part(t, absorbed), _by_part(t, released)
        emission === nothing || (mband_laser = _by_part(t, mband_laser))
    end
    return (; time, delivered, flux, Tr, moments = _columns(moments, lmax + 1),
            ylm = ylm ? _columns(clm, (lmax + 1)^2) : nothing, wall_Tr,
            wall_albedo, escaped, transmitted, window_loss, front, burnthrough,
            albedo = keep_albedo ? albedo_steps : nothing, setup, laser, source, absorbed, lost,
            released, lip_radius, record = record === nothing ? nothing : [r for r in recorded],
            mband_moments = emission === nothing ? nothing : _columns(mband_moments, lmax + 1),
            mband = emission === nothing ? nothing : mband_laser)
end

# The computed entrance-hole lips of a target (`LipClosure` apertures): per lip its name,
# model, hole center and normal, initial (outer) radius `r0`, inner radius `rin` (the
# floor of the closure), elements and each element's radial span [ra, rb] from the
# center (its vertices' extreme radii in the hole's plane).
function _lip_geometry(t::Target)
    lips = []
    for (ci, c) in enumerate(t.components)
        c.kind === :aperture && c.material.albedo isa LipClosure || continue
        s = c.surface
        s isa DiskSurface || throw(ArgumentError("a LipClosure lip needs a DiskSurface (component $(c.name))"))
        es = t.comp_elems[ci]
        radius(v) = (d = v - s.center; norm(d - dot(d, s.normal) * s.normal))
        span(e) = extrema(radius(v) for k in t.elem_tris[e]
                          for v in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k]))
        sp = span.(es)
        ra, rb = first.(sp), last.(sp)
        push!(lips, (; name = c.name, model = c.material.albedo, center = s.center, normal = s.normal,
                     r0 = maximum(rb), rin = minimum(ra), elems = es, ra, rb,
                     area = t.elem_area[es]))
    end
    return lips
end

# Each lip's incident-flux radiation temperature (area mean of its elements' H).
_lip_temperatures(lips, sol) =
    [radiation_temperature(sum(ForwardDiff.value(sol.H[e]) * A for (e, A) in zip(L.elems, L.area)) / sum(L.area))
     for L in lips]

# The laser intensity through each lip's hole (TW/mm² over its initial area): the power
# of the beams entering through it (ports on the hole's side of the chamber center).
_lip_intensities(lips, beams) =
    [sum((ForwardDiff.value(b.power) for b in beams if dot(b.port, L.center) > 0); init = 0.0) / (π * L.r0^2)
     for L in lips]

# The wall-release Duhamel sum of element e at time τ over the steps taken: the deficit
# D (piecewise linear between the steps) convolved with the kernel's integral,
# Σₖ (Dₖ − Dₖ₋₁)/(tₖ − tₖ₋₁) [(τ − tₖ₋₁)^b − (τ − tₖ)^b].
function _deficit_sum(m::WallLossAlbedo, D, ts, e, τ)
    b = 1 + m.time_exponent
    G = 0.0
    for k in 2:length(ts)
        dD = D[k][e] - D[k-1][e]
        dD == 0 && continue
        G += dD / (ts[k] - ts[k-1]) * ((τ - ts[k-1])^b - (τ - ts[k])^b)
    end
    return G
end

# Per-step series keyed by component summed into the components' parts.
function _by_part(t::Target, d::AbstractDict)
    out = Dict{Symbol,Vector{Float64}}()
    for c in t.components
        haskey(d, c.name) || continue
        v = d[c.name]
        out[c.part] = haskey(out, c.part) ? out[c.part] .+ v : copy(v)
    end
    return out
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
  effective radius: ``E_{cap}``); `escaped` per aperture; `transmitted` per thin-foil or
  closing-gap component; `lost` to open space;
- `heating` — per surface component, what stays in it: unconverted laser plus absorbed
  x-rays, `laser − source + absorbed`;
- `released` — per component, the stored wall energy re-emitted under a falling drive
  (`drive_history`'s `wall_release`). `absorbed` is net of it: a wall component took in
  `absorbed + released` and holds `absorbed`.

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
            escaped = d(h.escaped), transmitted = d(h.transmitted), lost = ∫(h.lost), heating,
            released = haskey(h, :released) ? d(h.released) : Dict{Symbol,Float64}())
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
