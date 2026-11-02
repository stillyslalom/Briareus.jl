# Time in the window: the pulses the beams fire (a pulse table, or a pulse shape scaled
# to the table's powers), the drive history they give on the drive's scenario (capsule
# Tᵣ and its Legendre modes, or a flat observed face's Tᵣ and nonuniformity, the
# synthetic diagnostics ticked in the list of diagnostic
# definitions and the energy split, `drive_history` and `energy_budget`), the latest
# retune's power histories, and the time cursor that colors the 3D view by a step of the
# history. Each step's element fields are kept with the history, so moving the cursor
# recomputes nothing.

# ---------------------------------------------------------------- pulses

"""
    load_pulses!(sess::DesignSession, path) -> sess

Fire the beams with the pulses of a pulse table (`read_pulse_table`: power per beam in
TW against `time_ns`, one column per cone group, quad or beam of the table; disabled
beams stay off). The time panel plots them, and [`start_history!`](@ref) computes the
drive they give.
"""
function load_pulses!(sess::DesignSession, path::AbstractString)
    table = read_pulse_table(path)
    by = _pulse_key(sess.beams.beams[], table)
    nodes = sort(unique(reduce(vcat, [p.t for p in values(table)])))
    _set_pulses!(sess, (; kind = :table, table, by, shape = nothing, nodes, times = nodes,
                        label = basename(path)))
    return sess
end

"""
    use_pulse_shape!(sess::DesignSession; step = 0.2) -> sess

Fire the beams with N221204's delivered pulse shape (Pak et al., PRE 109, 025203 (2024)
Fig. 1, `n221204_pulse`) scaled so that each beam reaches its table power at the peak: a
NIF pulse for a first look when no pulse table is at hand (the published curve is the
total of the 192 beams; every beam follows its shape here). The history is solved every
`step` ns.
"""
function use_pulse_shape!(sess::DesignSession; step::Real = 0.2)
    step > 0 || throw(ArgumentError("the step must be positive"))
    p = n221204_pulse()
    shape = PulseShape(p.t, p.power ./ maximum(p.power))
    times = collect(range(first(p.t), last(p.t); step))
    _set_pulses!(sess, (; kind = :shape, table = nothing, by = :beam, shape, nodes = collect(p.t), times,
                        label = "N221204's pulse shape × the table's powers"))
    return sess
end

function _set_pulses!(sess::DesignSession, src)
    tm = sess.time
    tm.pulses[] = src
    tm.history[] === nothing && _assign!(tm.chart, :pulses)
    tm.status[] = "pulses: " * src.label
    return sess
end

# The table's beams at time τ of the pulses `src` (disabled beams off).
function _pulse_beams(src, beams, enabled, τ)
    bs = src.kind === :table ? with_pulses(beams, src.table, τ; by = src.by) :
         [remake(b; power = b.power * src.shape(τ)) for b in beams]
    return [on ? b : remake(b; power = 0.0) for (b, on) in zip(bs, enabled)]
end

# ---------------------------------------------------------------- the history job

struct _HistoryCanceled <: Exception end

# What a history is computed for: the scenario, the beams (pointing and table power),
# which are on, the pulses, the walls and the unconverted light (fractions and prepulse).
_history_key(sess::DesignSession) = (sess.drive.scenario[], sess.beams.beams[], sess.beams.enabled[],
                                     sess.time.pulses[], wall_names(sess.target), _history_unconverted(sess),
                                     _history_prepulse(sess))

# The unconverted light of a drive history of the table's beams `idx`, as
# `drive_history`'s and `retune_drive`'s `unconverted`: per color, the drive's fraction
# of the 3ω power (`unconverted_ratios`), `(; omega1, omega2)`, one number when it is
# the same on every beam of `idx`, else a vector with one per beam of the table (zero on
# the beams whose light does not take part); `nothing` when no unconverted light of the
# beams `idx` takes part (the history is then the 3ω one, bit for bit).
function _history_unconverted(sess::DesignSession, idx = findall(sess.beams.enabled[]))
    u = unconverted_ratios(sess)
    u === nothing && return nothing
    f = map((u.omega1, u.omega2)) do o
        v = unique(o[idx])
        length(v) > 1 ? copy(o) : isempty(v) ? 0.0 : v[1]
    end
    all(x -> x isa Real && iszero(x), f) && return nothing
    return (; omega1 = f[1], omega2 = f[2])
end

"""
    set_prepulse!(tm::TimeState, omega1, omega2; start = current, stop = current) -> tm

The unconverted light's prepulse: a 1ω and a 2ω power per beam (TW, ≥ 0; 0 turns that
color's prepulse off), flat from `start` to `stop` (ns) and off outside, on the enabled
beams whose light of that color takes part ([`set_participation!`](@ref)). It is light
not proportional to the 3ω pulse (`drive_history`'s, `retune_drive`'s and
`unconverted_history`'s `prepulse`), so it adds to the drive's 1ω and 2ω ratios
([`set_unconverted_ratio!`](@ref)): a beam's power in a color is `ratio · P₃ω(t) +
prepulse(t)`. The drive histories, the retune and the irradiance checks
([`start_irradiance!`](@ref)) take it; the drive at the table's powers, which has no
time, does not. A drive history solves at the pulses' times where they fire, so a
prepulse before them reaches the irradiance histories only. There is no default level:
the NIF User Guide (§4.4) puts the 1ω prepulse above 10⁸ W/cm² 3–5 ns before the nominal
start of the pulse and gives no power (as a *derived* scale, 4.8e-4 TW of 1ω on a
beamlet puts its 22 mm footprint at 10⁸ W/cm²; research notes, "Unconverted light over
time"). A change makes the drive history and the checks out of date.
"""
function set_prepulse!(tm::TimeState, omega1::Real, omega2::Real; start::Real = tm.prepulse[][3],
                       stop::Real = tm.prepulse[][4])
    all(isfinite, (omega1, omega2, start, stop)) && omega1 >= 0 && omega2 >= 0 ||
        throw(ArgumentError("the prepulse powers must be finite and ≥ 0"))
    start < stop || throw(ArgumentError("the prepulse must start before it stops"))
    _assign!(tm.prepulse, (Float64(omega1), Float64(omega2), Float64(start), Float64(stop)))
    return tm
end

# The prepulse of the table's beams `idx` (`set_prepulse!`): per color, each table
# beam's power (TW; zero where the beam is not in `idx` or its light of that color does
# not take part) and the window (ns), `(; omega1, omega2, window)`; `nothing` when no
# beam of `idx` carries one.
function _history_prepulse(sess::DesignSession, idx = findall(sess.beams.enabled[]))
    p1, p2, a, b = sess.time.prepulse[]
    u = sess.beams.unconverted[]
    on = falses(length(u))
    on[collect(Int, idx)] .= true
    v1 = [on[i] && u[i][1] ? p1 : 0.0 for i in eachindex(u)]
    v2 = [on[i] && u[i][2] ? p2 : 0.0 for i in eachindex(u)]
    any(>(0), v1) || any(>(0), v2) || return nothing
    return (; omega1 = v1, omega2 = v2, window = (a, b))
end

# `drive_history`'s `prepulse` for the prepulse `pp` (`_history_prepulse`) of `beams`
# (in the order of `pp`'s vectors), keyed by `by`: per color with power, each key's flat
# pulse, one `PulseShape` per power, so that the beams on equal powers share a
# deposition. The beams of one key must carry the same power.
function _prepulse_schedules(pp, beams, by)
    key = Briareus._beam_key(by)
    a, b = pp.window
    out = Dict{Symbol,Any}()
    for c in (:omega1, :omega2)
        v = pp[c]
        any(>(0), v) || continue
        shapes = Dict{Float64,PulseShape}()
        d = Dict{Any,PulseShape}()
        for (i, bm) in enumerate(beams)
            k = key(bm)
            s = get!(() -> PulseShape([a, b], [v[i], v[i]]), shapes, v[i])
            haskey(d, k) && d[k] !== s &&
                throw(ArgumentError("the beams of $k carry different prepulses (the light of some of them " *
                                    "only takes part): retune by beam"))
            d[k] = s
        end
        out[c] = d
    end
    return out
end

# `drive_history`'s keywords for the unconverted light: the fractions `u`
# (`_history_unconverted`, zero for both colors when only a prepulse is set) and the
# prepulse `pp` (`_history_prepulse`) of `beams`, keyed by `by`. None without either,
# so that the history is the 3ω one bit for bit.
_light_kw(::Nothing, ::Nothing, beams; by = :beam) = (;)
function _light_kw(u, pp, beams; by = :beam)
    kw = (; unconverted = something(u, (omega1 = 0.0, omega2 = 0.0)), colors = nif_unconverted())
    pp === nothing && return kw
    return (; kw..., prepulse = _prepulse_schedules(pp, beams, by), by)
end

"""
    history_stale(sess::DesignSession) -> Bool

Whether the drive history was computed for another scenario, other beams or on/off
states, other pulses, other walls ([`wall_names`](@ref)) or other unconverted light
(the drive's ratios, [`set_unconverted_ratio!`](@ref), the prepulse,
[`set_prepulse!`](@ref), and the beams' 1ω and 2ω participation) (`false` without a
history).
"""
function history_stale(sess::DesignSession)
    h = sess.time.history[]
    h === nothing && return false
    k = _history_key(sess)
    return drive_stale(sess) || h.key[1] !== k[1] || !isequal(h.key[2:end], k[2:end])
end

"""
    start_history!(sess::DesignSession) -> sess

Compute the drive history of the table's beams on the pulses ([`load_pulses!`](@ref),
[`use_pulse_shape!`](@ref)) with the drive's deposition operator and radiation model
(`drive_history`, quasi-static: one radiosity solve per time of the pulses where they
fire), recording each step's element fields for the time cursor, each ticked diagnostic
([`set_diagnostic!`](@ref); `dante` for the soft x-ray power kind), and the energy split
(`energy_budget`). The unconverted 1ω and 2ω light takes part as in the drive
([`set_unconverted_ratio!`](@ref), `drive_history`'s `unconverted`), with the prepulse
([`set_prepulse!`](@ref), `drive_history`'s `prepulse`), on the beams whose 1ω or 2ω
light takes part ([`participates`](@ref)).
It runs on a worker when `spawn[]` is set, with progress per step and
[`cancel_history!`](@ref). A drive that is missing or stale is computed first, as a
coloring by a drive field computes it.
"""
function start_history!(sess::DesignSession)
    tm, dr = sess.time, sess.drive
    tm.running[] && return sess
    src = tm.pulses[]
    src === nothing && (tm.status[] = "load a pulse table first (or use N221204's pulse shape)"; return sess)
    if !_drive_current(sess)
        tm.queued[] = true
        tm.status[] = "computing the drive first…"
        dr.building[] || build_scenario!(sess)
        return sess
    end
    tm.queued[] = false
    s = dr.scenario[]
    s.model === nothing && (tm.status[] = "no drive history: the drive is laser only" * _radiation_note(s.target);
                            return sess)
    isempty(s.observed) &&
        (tm.status[] = "the drive history follows the observed surface: give a sphere (a capsule) or a flat " *
                       "face the role Observed"; return sess)
    beams, enabled = sess.beams.beams[], sess.beams.enabled[]
    times = filter(τ -> sum(b -> b.power, _pulse_beams(src, beams, enabled, τ)) > 0, src.times)
    isempty(times) && (tm.status[] = "the pulses fire no enabled beam"; return sess)
    light = _light_kw(_history_unconverted(sess), _history_prepulse(sess), beams)
    key = _history_key(sess)
    g = (tm.generation[] += 1)
    tm.cancel[] = false
    n = length(times)
    tm.running[] = true
    tm.progress[] = (0, n)
    tm.started[] = time()
    tm.status[] = "computing the drive history: $n steps…"
    diags = _recorded_definitions(sess)
    deliver = sess.deliver[]
    progress(k) = deliver(() -> g == tm.generation[] && _assign!(tm.progress, (k, n)))
    cancel = tm.cancel
    runner = _runner(sess.spawn, sess.deliver)
    wall = Tuple(wall_names(sess.target))
    cavities = cavity_walls(sess.target)
    runner(() -> _history_job(s, beams, enabled, src, times; diags, progress, cancel, key, wall, light,
                              cavities),
           function (out)
        g == tm.generation[] || return
        tm.running[] = false
        if out.err isa _HistoryCanceled
            tm.status[] = "drive history canceled"
        elseif out.err !== nothing
            tm.status[] = "drive history failed: " * _errmsg(out.err)
        else
            h = out.value
            tm.history[] = h
            _assign!(tm.chart, :drive)
            k = argmax(h.Tr)
            tm.status[] = @sprintf("%d steps in %s · peak Tr %.1f eV at %.2f ns", length(h.time),
                                   _seconds(h.seconds), h.Tr[k], h.time[k])
        end
    end)
    return sess
end

"""
    cancel_history!(sess::DesignSession) -> sess

Stop the drive history in progress at its next step (nothing is kept).
"""
cancel_history!(sess::DesignSession) = (sess.time.running[] && (sess.time.cancel[] = true); sess)

# The drive history of `beams` on the pulses `src` at `times` on scenario `s` (worker).
# `wall`: the wall-role holders (`wall_names`); those in the scenario's target give the
# history's wall temperature, else the observed surface does. `light`: the unconverted
# light's keywords of `drive_history` (`_light_kw`; empty: none). `cavities`: the hohlraum
# parts' walls (`cavity_walls`), each cavity's wall temperature per step.
function _history_job(s::Scenario, beams, enabled, src, times; diags = NamedTuple[],
                        progress = k -> nothing, cancel = Threads.Atomic{Bool}(false), key = nothing,
                        wall = (), light = (;), cavities = ())
    t0 = time()
    t = s.target
    names = component_names(t)
    observe = _observed_parts(s)
    walls = Tuple(n for n in wall if n in names)
    wall = isempty(walls) ? observe : walls
    views = _diagnostic_views(diags, names)
    A = t.elem_area
    aperture = [t.components[t.elem_comp[e]].kind === :aperture for e in 1:nelements(t)]
    steps = Ref(0)
    flat = _observed_flat(s)
    record = function (sol)
        cancel[] && throw(_HistoryCanceled())
        laser = Float32.(to_Wcm2.(sol.laser ./ A))
        laser[aperture] .= NaN32
        diag = [_record_diagnostic!(v, sol) for v in views]
        progress(steps[] += 1)
        # (per component: the history sums a part's pieces, a hohlraum's holes into its walls)
        escaped = energy_balance(sol; by = :component).escaped
        rec = (; laser, radiosity = Float32.(to_Wcm2.(sol.B)), incident = Float32.(to_Wcm2.(sol.H)), diag,
               escaped)
        # a flat observed face: its statistics over the face, as the drive's (`_face_stats`)
        return flat ? merge(rec, (; face = _face_stats(s, sol))) : rec
    end
    cavs = _drive_cavities(s, cavities)
    h = drive_history(s.model, s.deposition, τ -> _pulse_beams(src, beams, enabled, τ), times;
                      observe, wall, lmax = 4, record, cavities = cavs, light...)
    fields = Dict{Symbol,Matrix{Float32}}(f => reduce(hcat, [r[f] for r in h.record]) for f in DRIVE_FIELDS)
    ranges = Dict{Symbol,Tuple{Float64,Float64}}()
    for (f, M) in fields
        fin = filter(isfinite, M)
        lo, hi = isempty(fin) ? (0.0, 1.0) : Float64.(extrema(fin))
        f === :laser && (lo = 0.0)
        ranges[f] = (lo, hi > lo ? hi : lo + 1)
    end
    recorded = NamedTuple[_recorded(v, j, h, t) for (j, v) in enumerate(views)]
    if flat
        # the face's Tr from its area-mean incident flux, as the drive's `mean_flux` (the
        # core history's `Tr` is from the Legendre moment P0, meaningless on a plane)
        rs = [r.face for r in h.record]
        face = (; mean = [f.mean for f in rs], peak = [f.peak for f in rs], low = [f.low for f in rs],
                rms = [f.rms for f in rs])
        Tr = [radiation_temperature(m / to_Wcm2(1.0)) for m in face.mean]
        legendre = Matrix{Float64}(undef, 0, length(h.time))
    else
        face, Tr = nothing, h.Tr
        legendre = h.moments[2:5, :] ./ h.moments[1:1, :]
    end
    walls = NamedTuple[(; name = first(c), Tr = h.cavity_Tr[first(c)]) for c in cavs]
    return DriveHistory(h.time, Tr, legendre, face, recorded, _energy_series(h, t, observe),
                        fields, ranges, h, key, time() - t0, walls)
end

# Where the laser energy has gone by each time of a history (kJ, `energy_budget` up to
# that time): the laser delivered, the heat left in the walls (every surface but the
# capsule: unconverted laser and absorbed x-rays), in the capsule (the observed parts
# `observed`), the x-rays out of the entrance holes, and out of windows (other apertures
# and thin foils), and the unconverted 1ω and 2ω light delivered (part of the laser's).
function _energy_series(h, t::Target, observed)
    holes = Set(c.name for c in t.components if c.kind === :aperture && !(c.surface isa CylinderSurface))
    capsule = Set(c.part for c in t.components if c.name in observed || c.part in observed)
    eb = [energy_budget(h; to = τ) for τ in h.time]
    sumof(d, keep) = sum((v for (k, v) in d if keep(k)); init = 0.0)
    # the escape through each aperture up to each time, from the steps' records (the
    # budget's are by part), integrated as `energy_budget` integrates
    tt = h.time
    function upto(p)
        out, acc = zeros(length(p)), 0.0
        for k in 2:length(p)
            acc += (tt[k] - tt[k-1]) * (p[k-1] + p[k]) / 2
            out[k] = acc
        end
        return out
    end
    esc = isempty(h.record) ? Dict{Symbol,Vector{Float64}}() :
          Dict(n => upto([r.escaped[n] for r in h.record]) for n in keys(first(h.record).escaped))
    out(keep) = [sum((v[k] for (n, v) in esc if keep(n)); init = 0.0) for k in eachindex(tt)]
    series = NamedTuple[
        (; name = :laser, values = [e.setup for e in eb]),
        (; name = :wall, values = [sumof(e.heating, k -> !(k in capsule)) for e in eb]),
        (; name = :capsule, values = [sumof(e.heating, in(capsule)) for e in eb]),
        (; name = :holes, values = out(in(holes))),
        (; name = :windows, values = out(n -> !(n in holes)) .+ [sumof(e.transmitted, _ -> true) for e in eb])]
    # the unconverted light delivered, per color (part of the laser's)
    uc = h.unconverted
    uc === nothing || append!(series, [(; name = c, values = upto(uc.delivered[c])) for c in (:omega1, :omega2)])
    # (a hohlraum without windows, or a target without an observed surface, has no such series)
    return filter(s -> s.name === :laser || any(>(1e-9), s.values), series)
end

const ENERGY_LABELS = Dict(:laser => "Laser delivered", :wall => "Walls (heat)", :capsule => "Observed surface (heat)",
                           :holes => "Out of the entrance holes", :windows => "Out of windows",
                           :omega1 => "1ω delivered", :omega2 => "2ω delivered")

# ---------------------------------------------------------------- irradiance checks

"""
    set_irradiance_threshold!(tm::TimeState, threshold) -> tm

The irradiance (W/cm², > 0) the checks of [`start_irradiance!`](@ref) compare with
([`DEFAULT_THRESHOLD`](@ref) in a new session). The checks follow from the histories
kept: nothing is computed again.
"""
function set_irradiance_threshold!(tm::TimeState, x::Real)
    isfinite(x) && x > 0 || throw(ArgumentError("the threshold must be positive"))
    _assign!(tm.threshold, Float64(x))
    r = tm.irradiance[]
    r === nothing || r.threshold == x || (tm.irradiance[] = _with_checks(r, Float64(x)))
    return tm
end

# What the irradiance checks are computed for: the scenario, the beams, which are on,
# the pulses, the unconverted light (fractions and prepulse) and the active foil.
_irradiance_key(sess::DesignSession) = (sess.drive.scenario[], sess.beams.beams[], sess.beams.enabled[],
                                        sess.time.pulses[], _history_unconverted(sess), _history_prepulse(sess),
                                        active_foil(sess.target))

"""
    irradiance_stale(sess::DesignSession) -> Bool

Whether the irradiance checks were computed for another drive, other beams or on/off
states, other pulses, other unconverted light (ratios, prepulse, participation) or
another active foil (`false` without them).
"""
function irradiance_stale(sess::DesignSession)
    r = sess.time.irradiance[]
    r === nothing && return false
    k = _irradiance_key(sess)
    return drive_stale(sess) || r.key[1] !== k[1] || !isequal(r.key[2:end], k[2:end])
end

"""
    start_irradiance!(sess::DesignSession) -> sess

Compute the irradiance histories of the unconverted 1ω and 2ω light on the drive's
target through the pulses (`unconverted_history`: the drive's ratios on the beams whose
light takes part, [`set_unconverted_ratio!`](@ref) and [`participates`](@ref), with the
prepulse, [`set_prepulse!`](@ref); one `unconverted_deposition` per set of beams on one
schedule) and check them against the threshold ([`set_irradiance_threshold!`](@ref)):
`specular_check` per surface part, the first crossing on any part
(`threshold_crossing`), and, on a target with a foil, `preplasma_check` on the active
foil against the arrival of the 3ω light there (the first crossing of its
`irradiance_history`, on the drive's deposition operator). The histories are sampled at
the pulses' times and knots, at the prepulse's edges, and at the instants just outside
each pulse's and the prepulse's ends (where the power jumps from zero, so that a
crossing there is not drawn out over the step before it); they are linear in between,
as the pulses are. The result goes to `sess.time.irradiance` (see
[`irradiance_checks`](@ref)), from a worker when `spawn[]` is set; a drive that is
missing or stale is computed first.
"""
function start_irradiance!(sess::DesignSession)
    tm, dr = sess.time, sess.drive
    tm.checking[] && return sess
    src = tm.pulses[]
    src === nothing && (tm.status[] = "load a pulse table first (or use N221204's pulse shape)"; return sess)
    ratio, pp = _history_unconverted(sess), _history_prepulse(sess)
    ratio === nothing && pp === nothing &&
        (tm.status[] = "no unconverted light takes part: set a 1ω or 2ω ratio or prepulse (Drive properties) " *
                       "and turn on the beams' 1ω or 2ω dots"; return sess)
    if !_drive_current(sess)
        tm.check_queued[] = true
        tm.status[] = "computing the drive first…"
        dr.building[] || build_scenario!(sess)
        return sess
    end
    tm.check_queued[] = false
    s = dr.scenario[]
    beams, enabled = sess.beams.beams[], sess.beams.enabled[]
    foil = active_foil(sess.target)
    key = _irradiance_key(sess)
    threshold = tm.threshold[]
    g = (tm.check_generation[] += 1)
    tm.checking[] = true
    tm.status[] = "computing the irradiance histories…"
    t0 = time()
    runner = _runner(sess.spawn, sess.deliver)
    runner(() -> _irradiance_job(s, beams, enabled, src, ratio, pp, foil; key, threshold), function (out)
        g == tm.check_generation[] || return
        tm.checking[] = false
        if out.err !== nothing
            tm.status[] = "irradiance checks failed: " * _errmsg(out.err)
            return
        end
        r = out.value
        tm.threshold[] == r.threshold || (r = _with_checks(r, tm.threshold[]))
        tm.irradiance[] = r
        _assign!(tm.chart, :irradiance)
        tm.status[] = @sprintf("irradiance histories at %d times in %s", length(r.history.time),
                               _seconds(time() - t0))
    end)
    return sess
end

# Each table beam's 3ω power history (TW per beam), keyed by its name (`by = :beam`):
# the pulse table's function of its key (shared by the beams of that key), or the shape
# scaled to the beam's table power (one per power), and one zero pulse for the disabled
# beams; beams on the same function share a deposition.
function _beam_pulses(src, beams, enabled)
    allunique(b.name for b in beams) || throw(ArgumentError("the irradiance checks need unique beam names"))
    off = PulseShape([0.0], [0.0])
    scaled = Dict{Float64,PulseShape}()
    key = src.kind === :table ? Briareus._beam_key(src.by) : nothing
    out = Dict{Symbol,Any}()
    for (b, on) in zip(beams, enabled)
        out[Symbol(b.name)] = !on ? off :
                              src.kind === :table ? Briareus._beam_pulse(src.table, key, b) :
                              get!(() -> PulseShape(src.shape.t, b.power .* src.shape.p), scaled, b.power)
    end
    return out
end

# The times of the irradiance histories: the pulses' times and knots, and the ends of
# each pulse (zero outside its knots) and of the prepulse with the instants just outside.
function _irradiance_times(src, pp)
    ts = Float64[src.times; src.nodes]
    shapes = src.kind === :table ? collect(values(src.table)) : [src.shape]
    for p in shapes
        isempty(p.t) || append!(ts, (prevfloat(p.t[1]), nextfloat(p.t[end])))
    end
    if pp !== nothing
        a, b = pp.window
        append!(ts, (prevfloat(a), a, b, nextfloat(b)))
    end
    return sort!(unique!(ts))
end

# The irradiance histories (worker): the unconverted light's on the scenario's target,
# the 3ω light's when `foil` is one of its parts, and per surface part the peak of the
# unconverted irradiance at each time (all colors summed, W/cm²), with the checks at
# `threshold` (`_with_checks`).
function _irradiance_job(s::Scenario, beams, enabled, src, ratio, pp, foil; key = nothing,
                         threshold = DEFAULT_THRESHOLD)
    t = s.target
    times = _irradiance_times(src, pp)
    pulses = _beam_pulses(src, beams, enabled)
    h = unconverted_history(t, beams, pulses, times; ratio = something(ratio, (omega1 = 0.0, omega2 = 0.0)),
                            prepulse = pp === nothing ? (;) : _prepulse_schedules(pp, beams, :beam), by = :beam)
    peak(hh, p) = vec(maximum(Briareus._part_series(hh, (p,), nothing)[1]; dims = 1, init = 0.0))
    parts = unique!([c.part for c in t.components if c.kind === :surface])
    peaks = [(; part = p, values = peak(h, p)) for p in parts]
    f = foil !== nothing && any(c -> c.part === foil && c.kind === :surface, t.components) ? foil : nothing
    arrival = f === nothing ? nothing : irradiance_history(s.deposition, beams, pulses, times; by = :beam)
    r = (; history = h, arrival, foil = f, peaks, arrival_peak = arrival === nothing ? nothing : peak(arrival, f),
         key, threshold = NaN, specular = NamedTuple[], crossing = (; time = Inf, part = :none, element = 0),
         preplasma = nothing)
    return _with_checks(r, threshold)
end

# The irradiance histories `r` with their checks at `threshold`.
function _with_checks(r, threshold::Float64)
    h = r.history
    specular = specular_check(h; threshold)
    crossing = threshold_crossing(h, threshold)
    preplasma = r.arrival === nothing ? nothing : preplasma_check(h, r.arrival; foil = r.foil, threshold)
    return merge(r, (; threshold, specular, crossing, preplasma))
end

"""
    irradiance_checks(sess::DesignSession) -> Vector{NamedTuple}

The irradiance checks of the latest [`start_irradiance!`](@ref), one row `(; check,
pass, text, tip)` each: `:specular` (every surface part below the threshold through the
pulses, `specular_check`; `text` the first crossing, `threshold_crossing`, when one is
not), and on a target with a foil `:preplasma` (the unconverted light below the
threshold on the active foil until the 3ω light arrives, `preplasma_check`). `tip`
lists the parts' peaks. Empty without checks.
"""
function irradiance_checks(sess::DesignSession)
    r = sess.time.irradiance[]
    r === nothing && return NamedTuple[]
    thr = _irr_text(r.threshold)
    rows = NamedTuple[]
    c = r.crossing
    ok = isinf(c.time)
    hot = sort(r.specular; by = x -> -x.peak)
    tip = "Peak unconverted irradiance per part:\n" *
          join((@sprintf("%s: %s at %.2f ns%s", x.part, _irr_text(x.peak), x.peak_time,
                         x.specular ? "" : @sprintf(", above %s from %.2f ns", thr, x.crossing))
                for x in hot if x.peak > 0), "\n")
    text = ok ? "Specular: every part below " * thr :
           @sprintf("Specular: %s reaches %s at %.2f ns", c.part, thr, c.time) *
           (count(x -> !x.specular, r.specular) > 1 ? @sprintf(" (%d parts)", count(x -> !x.specular, r.specular)) : "")
    push!(rows, (; check = :specular, pass = ok, text, tip))
    p = r.preplasma
    if p !== nothing
        arrive = isinf(p.arrival) ? "the 3ω light stays below it" : @sprintf("3ω arrives at %.2f ns", p.arrival)
        text = if p.preplasma
            isinf(p.arrival) ? @sprintf("Pre-plasma on %s: from %.2f ns (3ω below %s)", r.foil, p.onset, thr) :
            @sprintf("Pre-plasma on %s: from %.2f ns, %.2f ns before the 3ω light", r.foil, p.onset, p.lead)
        else
            isinf(p.onset) ? "Pre-plasma on $(r.foil): none (below $thr)" :
            @sprintf("Pre-plasma on %s: none before the 3ω light (%s from %.2f ns)", r.foil, thr, p.onset)
        end
        push!(rows, (; check = :preplasma, pass = !p.preplasma,
                     text, tip = "The unconverted light on $(r.foil) against the arrival of the 3ω light: " *
                                 "the first time each reaches $thr there; " * arrive))
    end
    return rows
end

# An irradiance (W/cm²) as text, three significant digits.
_irr_text(x::Real) = @sprintf("%.3g W/cm²", x)

# ---------------------------------------------------------------- the time cursor

"""
    set_time_cursor!(sess::DesignSession, t::Real) -> sess

Show the drive history at time `t` (ns) in the 3D view: the history step nearest `t`
([`time_step`](@ref)) colors the target by the chosen drive field (by the incident x-ray
flux when the target was colored by component), on one color scale for the whole
history. `NaN` returns to the drive at the table's powers. The fields were recorded
with the history: nothing is computed.
"""
function set_time_cursor!(sess::DesignSession, t::Real)
    tm = sess.time
    _assign!(tm.cursor, Float64(t))
    isnan(t) || tm.history[] === nothing || sess.view.color[] !== :component || set_color!(sess, :incident)
    return sess
end

"""
    set_time_step!(sess::DesignSession, k::Integer) -> sess

Put the time cursor on step `k` of the drive history.
"""
function set_time_step!(sess::DesignSession, k::Integer)
    h = sess.time.history[]
    h === nothing && return sess
    return set_time_cursor!(sess, h.time[clamp(k, 1, length(h.time))])
end

"""
    time_step(tm::TimeState) -> Union{Nothing,Int}

The history step the cursor shows: the solved time nearest it (the first or last step
outside the history); `nothing` without a history or a cursor.
"""
function time_step(tm::TimeState)
    h, c = tm.history[], tm.cursor[]
    (h === nothing || isnan(c)) && return nothing
    return argmin(abs.(h.time .- c))
end

"""
    cursor_values(sess::DesignSession, field::Symbol) -> Union{Nothing,Vector{Float64}}

Per-element values (W/cm², on the drive's elements) of `field` at the cursor's history
step, or `nothing` when the cursor is off or the history is out of date.
"""
function cursor_values(sess::DesignSession, field::Symbol)
    k = time_step(sess.time)
    (k === nothing || history_stale(sess)) && return nothing
    return Float64.(view(sess.time.history[].fields[field], :, k))
end

"""
    cursor_range(sess::DesignSession, field::Symbol) -> Union{Nothing,Tuple}

The color scale of `field` over the whole history (W/cm²) while the cursor shows it.
"""
cursor_range(sess::DesignSession, field::Symbol) =
    cursor_values(sess, field) === nothing ? nothing : sess.time.history[].range[field]

"""
    time_label(sess::DesignSession) -> String

What the 3D view shows of the history: `"t = 6.20 ns"` at the cursor's step, `""` when
it shows the drive at the table's powers.
"""
function time_label(sess::DesignSession)
    k = time_step(sess.time)
    (k === nothing || history_stale(sess)) && return ""
    return @sprintf("t = %.2f ns", sess.time.history[].time[k])
end

# ---------------------------------------------------------------- charts

const _INK = "#333333"
const _GRAY = "#8C8C8C"
# the irradiance chart: the parts drawn, at most, and the floor of its log axis (a
# fraction of the threshold; display choices)
const _IRRADIANCE_PARTS = 12
const _IRRADIANCE_FLOOR = 1e-4

_series(label, x, y, color; style = :solid) =
    (; label = String(label), x = collect(Float64, x), y = collect(Float64, y), color, style)

"""
    time_chart(sess::DesignSession, chart = sess.time.chart[]) -> NamedTuple

What the time panel plots: `(; chart, title, xlabel, left, right, legend, cursor,
empty)`. `left` and `right` are the y axes, `(; label, series)` (`right` is `nothing`
for one axis), each series `(; label, x, y, color, style)` with `style` `:solid`,
`:dash` or `:dot`; `legend` the key, `(; label, color, style)` rows; `cursor` the time
of the cursor's history step (`NaN` without); `empty` a message when there is nothing to
plot (then no series); `scale` (`:log`), on a chart whose left axis is logarithmic only.

- `:pulses`: the laser power of each beam group (the sum over its enabled beams, TW;
  groups colored as in the model tree) and the total.
- `:drive`: the capsule's Tr and each recorded diagnostic's (eV, dashed), and a₂/a₀ and
  a₄/a₀ (%); for a flat observed face, its Tr (from its area-mean flux) and, on the
  right, its RMS nonuniformity and its peak and lowest flux relative to its mean (%,
  dashed), as [`DriveHistory`](@ref)'s `face`.
- `:diagnostics`: each recorded diagnostic's Tr (eV) and the intensity toward it
  (TW/sr, dashed).
- `:energy`: the laser energy delivered by each time and where it went (kJ), and the
  unconverted 1ω and 2ω light's part of the delivered energy (dotted) when it takes part.
- `:irradiance`: from the latest [`start_irradiance!`](@ref), the peak irradiance of the
  unconverted light (1ω and 2ω summed, W/cm²) on each surface part it reaches (the
  12 highest), the 3ω light's on the active foil (dotted) and the threshold (dashed), on
  a log axis (`scale = :log`; values below 10⁻⁴ of the threshold drawn at that level).
- `:retune`: each group's power on the reference pulses (all beams, dashed) and on the
  retuned ones (enabled beams), and the capsule's mode 1 ε₁ before and after (%).
"""
function time_chart(sess::DesignSession, chart::Symbol = sess.time.chart[])
    chart in TIME_CHARTS || throw(ArgumentError("chart must be one of $TIME_CHARTS"))
    tm = sess.time
    h = tm.history[]
    k = time_step(tm)
    cursor = k === nothing ? NaN : h.time[k]
    stale = h !== nothing && history_stale(sess) ? " (out of date: run it again)" : ""
    none(title, msg) = (; chart, title, xlabel = "t (ns)", left = (; label = "", series = NamedTuple[]),
                        right = nothing, legend = NamedTuple[], cursor = NaN, empty = msg)
    if chart === :pulses
        src = tm.pulses[]
        src === nothing && return none("Pulses", "Load a pulse table (power per beam against time), or use " *
                                                 "N221204's pulse shape")
        series = _group_series(sess.beams, src.nodes, τ -> _pulse_beams(src, sess.beams.beams[],
                                                                       sess.beams.enabled[], τ))
        return (; chart, title = "Laser pulses: " * src.label, xlabel = "t (ns)",
                left = (; label = "Laser power (TW)", series), right = nothing,
                legend = [_key(s) for s in series], cursor, empty = "")
    elseif chart === :drive
        h === nothing && return none("Drive history", tm.pulses[] === nothing ?
                                                      "Load pulses, then Run the drive history" :
                                                      "Run the drive history on these pulses")
        f = h.face
        left = [_series(f === nothing ? "Capsule Tr" : "Face Tr", h.time, h.Tr, _INK)]
        # with several cavities, each one's wall temperature (colored as its wall); one
        # cavity's wall is in the drive's summary, its observed surface the headline
        if length(h.cavities) > 1
            cc = component_colors(sess.drive.scenario[] === nothing ? sess.target.target[] :
                                  sess.drive.scenario[].target)
            walls = Dict(cavity_walls(sess.target))
            for (k, c) in enumerate(h.cavities)
                w = get(walls, c.name, Symbol[])
                color = isempty(w) ? palette_color(k) : get(cc, first(w), palette_color(k))
                push!(left, _series(String(c.name) * " wall Tr", h.time, c.Tr, color; style = :dot))
            end
        end
        for (j, d) in enumerate(h.diagnostics)
            any(isfinite, d.Tr) && push!(left, _series(d.name, h.time, d.Tr, _diagnostic_color(j); style = :dash))
        end
        if f === nothing
            right = [_series("a₂/a₀", h.time, 100 .* h.legendre[2, :], "#0072B2"),
                     _series("a₄/a₀", h.time, 100 .* h.legendre[4, :], "#009E73")]
            title, rlabel = "Capsule drive history", "aₗ/a₀ (%)"
        else
            # a flat observed face: its RMS nonuniformity, and its peak and lowest flux
            # relative to its mean
            right = [_series("RMS", h.time, 100 .* f.rms, "#0072B2"),
                     _series("peak", h.time, 100 .* (f.peak ./ f.mean .- 1), "#D55E00"; style = :dash),
                     _series("lowest", h.time, 100 .* (f.low ./ f.mean .- 1), "#009E73"; style = :dash)]
            title, rlabel = "Observed face drive history", "Face nonuniformity (%)"
        end
        return (; chart, title = title * stale, xlabel = "t (ns)",
                left = (; label = "Tr (eV)", series = left), right = (; label = rlabel, series = right),
                legend = [_key(s) for s in vcat(left, right)], cursor, empty = "")
    elseif chart === :diagnostics
        h === nothing && return none("Diagnostics", "Tick diagnostics below, then Run the drive history")
        isempty(h.diagnostics) && return none("Diagnostics", "No diagnostic was recorded: tick one, then Run the drive history")
        left = [_series(d.name * " Tr", h.time, d.Tr, _diagnostic_color(j))
                for (j, d) in enumerate(h.diagnostics) if any(isfinite, d.Tr)]
        right = [_series(d.name * (d.kind === :xray_imager ? " detected" : " intensity"), h.time, d.intensity,
                         _diagnostic_color(j); style = :dash)
                 for (j, d) in enumerate(h.diagnostics) if any(isfinite, d.intensity)]
        notes = [d.name * ": " * d.note for d in h.diagnostics if !isempty(d.note)]
        isempty(left) && isempty(right) &&
            return none("Diagnostics", "No recorded diagnostic sees the target: " * join(notes, "; "))
        return (; chart, title = "Synthetic diagnostics" * stale, xlabel = "t (ns)",
                left = (; label = "Tr (eV)", series = left),
                right = isempty(right) ? nothing : (; label = "Intensity (TW/sr)", series = right),
                legend = [_key(s) for s in vcat(left, right)], cursor, empty = "")
    elseif chart === :energy
        h === nothing && return none("Energy", "Run the drive history to see where the energy goes")
        cc = component_colors(sess.drive.scenario[] === nothing ? sess.target.target[] :
                              sess.drive.scenario[].target)
        # the colors of the first wall and of the observed surface
        ws = wall_names(sess.target)
        w, o = isempty(ws) ? nothing : ws[1], observed_name(sess.target)
        color = Dict(:laser => _INK, :wall => get(cc, w, palette_color(1)),
                     :capsule => get(cc, o, palette_color(4)), :holes => _GRAY, :windows => "#D55E00",
                     :omega1 => "#CC79A7", :omega2 => "#E69F00")
        series = [_series(ENERGY_LABELS[s.name], h.time, s.values, color[s.name];
                          style = s.name === :laser ? :dash : s.name in (:omega1, :omega2) ? :dot : :solid)
                  for s in h.energy]
        return (; chart, title = "Energy by each time" * stale, xlabel = "t (ns)",
                left = (; label = "Energy (kJ)", series), right = nothing,
                legend = [_key(s) for s in series], cursor, empty = "")
    elseif chart === :irradiance
        r = tm.irradiance[]
        r === nothing && return none("Irradiance", tm.pulses[] === nothing ? "Load pulses, then Run checks" :
                                                   "Run checks on these pulses")
        lit = sort(filter(p -> maximum(p.values; init = 0.0) > 0, r.peaks); by = p -> -maximum(p.values))
        isempty(lit) && r.arrival_peak === nothing &&
            return none("Irradiance", "The unconverted light reaches no surface of the drive's target")
        ts = r.history.time
        lo = r.threshold * _IRRADIANCE_FLOOR
        cc = component_colors(r.history.target)
        pc = Dict{Symbol,String}()
        for c in r.history.target.components
            c.kind === :surface && haskey(cc, c.name) && get!(pc, c.part, cc[c.name])
        end
        series = [_series(String(p.part), ts, max.(p.values, lo), get(pc, p.part, palette_color(k)))
                  for (k, p) in enumerate(first(lit, _IRRADIANCE_PARTS))]
        r.arrival_peak === nothing ||
            push!(series, _series("3ω on " * String(r.foil), ts, max.(r.arrival_peak, lo), _INK; style = :dot))
        push!(series, _series("threshold " * _irr_text(r.threshold), [ts[1], ts[end]], [r.threshold, r.threshold],
                              _GRAY; style = :dash))
        st = irradiance_stale(sess) ? " (out of date: run it again)" : ""
        return (; chart, title = "Peak irradiance of the unconverted light" * st, xlabel = "t (ns)",
                left = (; label = "Irradiance (W/cm²)", series), right = nothing,
                legend = [_key(s) for s in series], cursor, empty = "", scale = :log)
    else
        r = tm.retune[]
        r === nothing && return none("Retune", "Retune to a reference drive (Compute ▸ Retune to a reference drive) to see its pulses")
        ref = _group_series(r.table, r.nodes, τ -> with_pulses(r.beams, r.reference_pulses, τ; by = r.reference_by);
                            total = false, style = :dash)
        new = _group_series(r.table, r.nodes,
                            τ -> [on ? remake(b; power = r.pulses[Briareus._beam_key(r.by)(b)](τ)) : remake(b; power = 0.0)
                                  for (b, on) in zip(r.beams, r.enabled)]; total = false)
        ε1(hh) = [100 * mode_amplitude(hh.ylm[:, j], 1) for j in eachindex(hh.time)]
        right = [_series("ε₁ before", r.before.time, ε1(r.before), _GRAY; style = :dot),
                 _series("ε₁ retuned", r.after.time, ε1(r.after), _INK)]
        legend = vcat([_key(s) for s in new], [(; label = "reference (dashed)", color = _GRAY, style = :dash)],
                      [_key(s) for s in right])
        return (; chart, title = "Retune to " * r.label, xlabel = "t (ns)",
                left = (; label = "Laser power (TW)", series = vcat(ref, new)),
                right = (; label = "ε₁ (%)", series = right), legend, cursor, empty = "")
    end
end

_key(s) = (; s.label, s.color, s.style)

# Power per beam group (the sum over its beams, TW) at the times `x` of the beams
# `beams_at(τ)` (the table's order), groups colored as in the model tree; with the total
# when there is more than one group. `table` is a `BeamTable`, or a NamedTuple of a
# table's `cluster` and `colors` at the time (the retune's).
function _group_series(table, x, beams_at; total::Bool = true, style::Symbol = :solid)
    cluster, colors = table isa BeamTable ? (table.cluster[], cluster_colors(table)) : (table.cluster, table.colors)
    groups = unique(cluster)
    P = [b.power for τ in x for b in beams_at(τ)]
    nb = length(cluster)
    M = reshape(P, nb, length(x))
    series = NamedTuple[]
    for g in groups
        rows = findall(==(g), cluster)
        y = vec(sum(M[rows, :]; dims = 1))
        any(>(0), y) || continue
        push!(series, _series(String(g), x, y, colors[g]; style))
    end
    total && length(series) > 1 && push!(series, _series("total", x, vec(sum(M; dims = 1)), _INK; style = :dash))
    return series
end

# ---------------------------------------------------------------- diagnostics

# Colors of the recorded diagnostics on the charts (Okabe–Ito, after those of Tr and
# the modes).
const _DIAGNOSTIC_COLORS = ("#D55E00", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
_diagnostic_color(j) = _DIAGNOSTIC_COLORS[mod1(j, length(_DIAGNOSTIC_COLORS))]

_builtin_file(f) = startswith(normpath(f), normpath(Briareus._BUILTIN_DIAGNOSTICS))

"""
    refresh_diagnostics!(sess::DesignSession) -> Vector{NamedTuple}

Read the diagnostic definitions on the search path again (`diagnostic_paths`: the
built-ins, the saved folders, `BRIAREUS_DIAGNOSTICS`) with their problems
(`check_diagnostics`, in `sess.time.issues`), and list them with those saved with the
session in `sess.time.catalog`, one row `(; id, name, kind, origin, file, chain,
public, usable, problem, changed)` per id: `origin` is `:builtin`, `:user` (a folder on
the search path) or `:session` (saved with the session and not on the search path);
`usable` whether it builds, else `problem` says why; `changed` the keys in which a
definition saved with the session differs from the search path's of that id
(`diagnostic_diff`; the session's is used). The built-ins first, then by id.
"""
function refresh_diagnostics!(sess::DesignSession)
    tm = sess.time
    issues = check_diagnostics()
    raws = Briareus._read_all(diagnostic_paths(), DiagnosticIssue[])
    found = Briareus._checked(raws, DiagnosticIssue[])
    problem(d) = join([i.message * (isempty(i.key) ? "" : " [" * i.key * "]") for i in issues
                       if i.file == d.file && i.level !== :warning], "; ")
    # listed too: those that lack values (not public) but have no error
    for id in setdiff(keys(raws), keys(found))
        any(i -> i.file == raws[id].file && i.level === :error, issues) && continue
        d = Briareus._resolve(id, raws, DiagnosticIssue[])
        d === nothing || (found[id] = d)
    end
    rows = NamedTuple[]
    for id in union(keys(found), keys(tm.embedded[]))
        e = get(tm.embedded[], id, nothing)
        p = get(found, id, nothing)
        d = e === nothing ? p : e.def
        chain = e === nothing ? d.chain : e.chain
        msg = e === nothing ? problem(d) : ""
        push!(rows, (; id, d.name, d.kind, origin = p === nothing ? :session : _builtin_file(p.file) ? :builtin : :user,
                     file = e === nothing ? d.file : e.file, chain, d.public, usable = isempty(msg), problem = msg,
                     changed = e === nothing || p === nothing ? String[] : diagnostic_diff(e.def, p)))
    end
    sort!(rows; by = r -> (r.origin !== :builtin, r.id))
    tm.issues[] = issues
    tm.catalog[] = rows
    return rows
end

"""
    diagnostic_catalog(sess::DesignSession) -> Vector{NamedTuple}

The rows of [`refresh_diagnostics!`](@ref), read from the search path the first time.
"""
diagnostic_catalog(sess::DesignSession) =
    isempty(sess.time.catalog[]) ? refresh_diagnostics!(sess) : sess.time.catalog[]

# The definition a session uses for `id`: the one saved with it, else the search path's
# (`nothing` without either).
function _definition(sess::DesignSession, id::AbstractString)
    e = get(sess.time.embedded[], id, nothing)
    e === nothing || return e.def
    r = findfirst(r -> r.id == id, diagnostic_catalog(sess))
    r === nothing && return nothing
    row = diagnostic_catalog(sess)[r]
    return row.origin === :session ? nothing : Briareus._definitions(diagnostic_paths(), DiagnosticIssue[])[id]
end

"""
    set_diagnostic!(sess::DesignSession, id, on::Bool = true) -> sess

Record the diagnostic of definition `id` in the next drive history (`on`), or not.
"""
function set_diagnostic!(sess::DesignSession, id::AbstractString, on::Bool = true)
    tm = sess.time
    ids = filter(!=(String(id)), tm.diagnostics[])
    on && push!(ids, String(id))
    ids == tm.diagnostics[] || (tm.diagnostics[] = ids)
    return sess
end

"""
    use_path_diagnostic!(sess::DesignSession, id) -> sess

Use the search path's definition of `id` instead of the one saved with the session
(when the two differ, see [`refresh_diagnostics!`](@ref)).
"""
function use_path_diagnostic!(sess::DesignSession, id::AbstractString)
    e = copy(sess.time.embedded[])
    pop!(e, String(id), nothing) === nothing && return sess
    sess.time.embedded[] = e
    refresh_diagnostics!(sess)
    return sess
end

# The ticked definitions, in the order ticked: `(; id, name, def, note)`, `def` nothing
# (and `note` why) for one the session cannot build.
function _recorded_definitions(sess::DesignSession)
    rows = Dict(r.id => r for r in diagnostic_catalog(sess))
    out = NamedTuple[]
    for id in sess.time.diagnostics[]
        r = get(rows, id, nothing)
        if r === nothing
            push!(out, (; id, name = id, def = nothing, note = "not on the diagnostic search path"))
        elseif !r.usable
            push!(out, (; id, r.name, def = nothing, note = r.problem))
        else
            push!(out, (; id, r.name, def = _definition(sess, id), note = ""))
        end
    end
    return out
end

# The instrument of each recorded definition on a target with components `names`
# (worker): `(; id, name, view, note)`, `view` nothing (and `note` why) when it records
# nothing in a history: a kind recorded neither as a power history nor as images, a
# view through an aperture the target lacks, or a field centered on a part it lacks.
function _diagnostic_views(diags, names)
    out = NamedTuple[]
    for d in diags
        view, note = nothing, d.note
        if d.def !== nothing
            if !(d.def.kind in (:xray_power, :xray_imager))
                note = "a $(d.def.kind) diagnostic is not recorded in a drive history"
            else
                v = _view_on(diagnostic(d.def), names)
                if v isa DanteView && !(v.through === nothing || v.through in names)
                    note = "the target has no $(v.through) to look through"
                elseif v.center isa Symbol && !(v.center in names)
                    note = "the target has no $(v.center) to center the field on"
                else
                    view = v
                end
            end
        end
        push!(out, (; d.id, d.name, view, note = Ref(note)))
    end
    return out
end

# One step of a recorded diagnostic: `(; Tr, intensity, L)`, the brightness temperature
# (eV) and intensity (TW/sr) of a power diagnostic, or for an imager the detected
# radiance per element (`detected_radiance`, Float32) and NaN; NaN, with the reason kept,
# where the view gives none.
function _record_diagnostic!(v, sol)
    v.view === nothing && return (; Tr = NaN, intensity = NaN, L = nothing)
    v.view isa XrayImager &&
        return (; Tr = NaN, intensity = NaN, L = Float32.(detected_radiance(v.view, sol)))
    try
        r = dante(v.view, sol)
        return (; Tr = Float64(r.Tr), intensity = Float64(r.intensity), L = nothing)
    catch err
        err isa ArgumentError || rethrow()
        isempty(v.note[]) && (v.note[] = sprint(showerror, err))
        return (; Tr = NaN, intensity = NaN, L = nothing)
    end
end

# What history `h` keeps of recorded diagnostic `j` (worker): `(; id, name, kind, Tr,
# intensity, note, view, L, frames)`. For an imager, `L` the detected radiance per
# element and step, `intensity` the detected intensity toward the pinhole (TW/sr: the
# image's total over the object plane) and `frames` its recording
# (`imager_record`: one integrated image, or a frame per gate).
function _recorded(v, j, h, t::Target)
    rows = [r.diag[j] for r in h.record]
    if v.view isa XrayImager
        L = reduce(hcat, [r.L for r in rows])
        P = Briareus._imager_pixels(v.view, t)
        a = (v.view.camera.width / v.view.camera.nx)^2
        intensity = [a * sum((p > 0 ? Float64(L[p, k]) : 0.0 for p in P); init = 0.0) for k in axes(L, 2)]
        frames = imager_record(v.view, t, h.time, Float64.(L))
        return (; v.id, v.name, kind = :xray_imager, Tr = fill(NaN, length(rows)), intensity,
                note = v.note[], v.view, L, frames)
    end
    return (; v.id, v.name, kind = :xray_power, Tr = [r.Tr for r in rows], intensity = [r.intensity for r in rows],
            note = v.note[], v.view, L = nothing, frames = nothing)
end

"""
    diagnostic_images(sess::DesignSession) -> Vector{NamedTuple}

The images the latest drive history holds of its recorded imagers, one row `(; id,
label, frame)` per image: `frame` 0 for the image at the time cursor, else the
recording's frame (the integrated image, or each gate's). Empty without a history or
an imager.
"""
function diagnostic_images(sess::DesignSession)
    h = sess.time.history[]
    h === nothing && return NamedTuple[]
    out = NamedTuple[]
    for d in h.diagnostics
        d.frames === nothing && continue
        push!(out, (; d.id, label = d.name * " at the cursor", frame = 0))
        if isempty(d.frames.times)
            push!(out, (; d.id, label = d.name * ", time-integrated", frame = 1))
        else
            for (k, τ) in enumerate(d.frames.times)
                push!(out, (; d.id, label = @sprintf("%s, frame %d (%.2f ns)", d.name, k, τ), frame = k))
            end
        end
    end
    return out
end

"""
    show_diagnostic_image!(sess::DesignSession, id, frame::Integer = 0) -> sess

Show an image of recorded imager `id` of the latest drive history in the image view
(as Compute ▸ X-ray image does, [`set_view!`](@ref) `:image`): `frame` 0 the image at
the time cursor's step (detected radiance, W/cm²/sr; the last step without a cursor),
else frame `frame` of its recording (detected fluence per steradian, J/cm²/sr: the
integral over the shot or against a gate). The view opens on the region with signal
(`window`, or the whole field when that covers most of it); Ctrl + Shift + click on the
image shows the whole field. The image stays until another is asked
for: a drive edit does not replace it.
"""
function show_diagnostic_image!(sess::DesignSession, id::AbstractString, frame::Integer = 0)
    h = sess.time.history[]
    h === nothing && return sess
    j = findfirst(d -> d.id == id && d.frames !== nothing, h.diagnostics)
    j === nothing && return sess
    d = h.diagnostics[j]
    t = sess.drive.scenario[].target
    if frame == 0
        k = something(time_step(sess.time), length(h.time))
        img = to_Wcm2.(imager_image(d.view, t, Float64.(view(d.L, :, k))))
        title, unit = @sprintf("%s at %.2f ns", d.name, h.time[k]), "W/cm²/sr (detected)"
    else
        f = d.frames
        1 <= frame <= length(f.images) || throw(ArgumentError("no frame $frame"))
        img = 1e5 .* f.images[frame]                          # kJ/mm²/sr to J/cm²/sr
        title = isempty(f.times) ? d.name * ", time-integrated" :
                @sprintf("%s, frame %d (%.2f ns, gate %.2f ns)", d.name, frame, f.times[frame], d.view.gate)
        unit = "J/cm²/sr (detected)"
    end
    cam = d.view.camera
    w = cam.width / 2
    x = range(-w, w; length = cam.nx)
    values = permutedims(img[end:-1:1, :])
    sess.drive.image_request[] = nothing                     # not the drive's camera image
    sess.drive.image[] = (; x, y = x, values, title, unit, window = _signal_window(x, x, values))
    set_view!(sess.drive, :image)
    return sess
end

# The region of an image with signal, for the image view's first limits: `(x0, x1, y0,
# y1)` (mm) bounding the pixels above `frac` of the peak, widened by `margin` of its
# extent on each side (at least 4 pixels) and clamped to the field; `nothing` when that
# region covers most of the field or the image is empty. `values[i, j]` is at `(x[i],
# y[j])`. The threshold and margin are display defaults, not instrument criteria.
function _signal_window(x, y, values; frac = 1e-3, margin = 0.15)
    fin = filter(isfinite, values)
    peak = isempty(fin) ? 0.0 : maximum(fin)
    peak > 0 || return nothing
    hot = findall(v -> isfinite(v) && v > frac * peak, values)
    is, js = extrema(i[1] for i in hot), extrema(i[2] for i in hot)
    dx, dy = step(x), step(y)
    function span(r, lo, hi, d)
        pad = max(margin * (r[hi] - r[lo]), 4 * abs(d))
        return max(first(r) - abs(d) / 2, r[lo] - pad), min(last(r) + abs(d) / 2, r[hi] + pad)
    end
    x0, x1 = span(x, is..., dx)
    y0, y1 = span(y, js..., dy)
    (x1 - x0) * (y1 - y0) > 0.5 * (last(x) - first(x)) * (last(y) - first(y)) && return nothing
    return (x0, x1, y0, y1)
end

# The group of a diagnostic kind in the time panel's definition picker.
const _DIAGNOSTIC_GROUPS = (:xray_power => "X-ray power", :xray_imager => "X-ray imagers",
                            :radiography => "Radiography")
_diagnostic_group(kind::Symbol) =
    something(get(Dict(_DIAGNOSTIC_GROUPS), kind, nothing), uppercasefirst(replace(String(kind), '_' => ' ')))

# The catalog's rows grouped by kind for the picker (the groups of `_DIAGNOSTIC_GROUPS`
# in that order, then the other kinds by name; the catalog's order within a group):
# `(; row, group, head)`, `head` for the first row of its group.
function _diagnostic_groups(rows)
    rank(k) = something(findfirst(p -> p.first === k, _DIAGNOSTIC_GROUPS), length(_DIAGNOSTIC_GROUPS) + 1)
    order = sort(collect(eachindex(rows)); by = i -> (rank(rows[i].kind), _diagnostic_group(rows[i].kind)))
    out = NamedTuple[]
    for i in order
        g = _diagnostic_group(rows[i].kind)
        push!(out, (; row = rows[i], group = g, head = isempty(out) || out[end].group != g))
    end
    return out
end

"""
    diagnostic_folders() -> Vector{NamedTuple}

The folders on the diagnostic search path besides the built-ins, `(; path, saved)`:
`saved` for those kept in the configuration file (Diagnostics ▸ Folders), not for those
of `BRIAREUS_DIAGNOSTICS` or added for this Julia session only.
"""
function diagnostic_folders()
    saved = Set(normpath.(Briareus._saved_paths()))
    return [(; path, saved = normpath(path) in saved) for path in diagnostic_paths(; builtin = false)]
end

"""
    add_diagnostic_folder!(sess::DesignSession, dir) -> Bool

Add the folder `dir` to the diagnostic search path and keep it for later sessions
(`add_diagnostic_path`), then read the definitions again. Returns whether it was added
(the status says why not).
"""
function add_diagnostic_folder!(sess::DesignSession, dir::AbstractString)
    isdir(dir) || (sess.time.status[] = "no such folder: $dir"; return false)
    add_diagnostic_path(dir; save = true)
    rows = refresh_diagnostics!(sess)
    n = count(r -> r.origin === :user && startswith(normpath(r.file), normpath(abspath(dir))), rows)
    sess.time.status[] = "diagnostic folder added: $n definition$(n == 1 ? "" : "s") in $(basename(normpath(dir)))"
    return true
end

"""
    remove_diagnostic_folder!(sess::DesignSession, dir) -> sess

Take the folder `dir` off the saved diagnostic search path and read the definitions
again.
"""
function remove_diagnostic_folder!(sess::DesignSession, dir::AbstractString)
    remove_diagnostic_path(dir; save = true)
    refresh_diagnostics!(sess)
    return sess
end
