# This file covers the retune of the pulses of a session to a reference drive.

# ---------------------------------------------------------------- retune

# The key a pulse table's columns are: the beams' groups, quads or names.
function _pulse_key(beams, pulses)
    for by in (:group, :quad, :beam)
        key = Briareus.GUIAPI.beam_key(by)
        all(b -> haskey(pulses, key(b)), beams) && return by
    end
    throw(ArgumentError("the pulse table's columns are not the beams' groups, quads or names"))
end

"""
    retune_pulses(s::Scenario, beams, enabled, pulses; by = :quad, cap = 1.3, stay = 0.0,
                  target = :reference, level = nothing, wall = (), unconverted = nothing,
                  prepulse = nothing, foil = nothing, spot_radius = 0.5, on = :auto) -> NamedTuple

Retune the enabled beams' powers to the drive of a reference shot in which every beam of
`beams` (the scenario's, in its order) fired `pulses` (power per beam against time, keyed
by group, quad or beam name, e.g. `read_pulse_table`): at each step of a drive history
of the reference ([`start_history!`](@ref)'s steps: the pulses' times, the prepulse
window every 0.5 ns and one step after it, and the zero-power steps between, from the
first time anything fires to the last), where the reference's 3ω light fires (elsewhere
the retuned powers are zero), the drive is held (`Briareus.retune_drive`) by one power
per `by` key (`:quad`, `:group` or `:beam`) within `cap` times its reference power,
preferring the reference powers with weight `stay`. `target = :reference` holds the
reference's drive shape, `:uniform` a uniform drive (no modes; a flat flux), the least
nonuniformity the powers allow; `level` sets the drive's level at every step (TW/mm²:
the face's or the foil spot's mean, or the capsule's mean incident flux, whose `c_00` is
`√(4π)` times it), and `nothing` holds the reference's level at each step. Where the
power limits cannot reach the level the retune reports the level they allow
(`retune.drive`). Disabled beams are the dropped ones. A cone's or a quad's column of
`pulses` is the power per `unit` ([`PULSE_UNITS`](@ref): per quad by default, shared
among the quad's beams), and so are the new pulses' columns (the result's `unit`); a key
whose quads are held some as one beam and some as four is refused (its beams would not
share one power). What is held (`kind`):

- `:capsule`, an observed sphere: its modes through l = 4 (`Briareus.CapsuleModes`);
- `:face`, an observed flat face ([`observed_flat`](@ref)): its incident x-ray flux on
  each element, so its RMS nonuniformity, peak-to-mean and mean (`Briareus.FaceFlux`);
- `:foil`, the backlighter foil `foil` (a part or component of the scenario's target,
  e.g. [`active_foil`](@ref)): its absorbed laser irradiance (3ω, and the unconverted
  light below) within `spot_radius` (mm) of the foil's area centroid on its front, as the
  drive's spot statistics (`Briareus.FoilIrradiance`); laser only, no radiation model
  needed.

`on = :auto` holds the observed surface when the scenario has a radiation model and an
observed surface, else the foil; `:observed` or `:foil` chooses. Returns the new `pulses`
(by `by`, on the table's times and the steps, zero outside the steps and for the keys of
the dropped beams, so that every beam of `beams` has one), the `kind`, the `retune`, and the `reference`,
the enabled beams on the reference pulses (`before`) and on the new ones (`after`): for
a capsule, drive histories, whose wall temperature is that of the names `wall` in the
scenario's target (a session passes its [`wall_names`](@ref)), else the capsule's; for
a face or a foil, `(; time, mean, peak, low, rms)` per time (W/cm², the RMS relative), as
the drive's face and spot statistics, the face's from drive histories; and the `target`
and `level` asked for. With
`unconverted` (per color, the unconverted power as a fraction of each beam's 3ω power,
`(; omega1, omega2)`: a number for every beam, or a vector with one per beam of
`beams`), every history and the retune take the unconverted 1ω and 2ω light
(`drive_history`'s and `retune_drive`'s `unconverted`, NIF's colors; the enabled beams'
fractions in the retune); for a foil, each color's absorbed power adds to the spot's
irradiance in the reference, the retune and the statistics, as `drive_history` deposits
it. `prepulse` adds light not proportional to the 3ω pulse, `(; omega1, omega2,
window)`: per color a power per beam of `beams` (TW, a vector), flat over `window =
(start, stop)` (ns); the histories and the retune take it as `drive_history`'s and
`retune_drive`'s `prepulse`, keyed by `by` (the beams of one key must carry the same
prepulse), and the retune does not retune it. Without either they are identical to the
3ω ones.
"""
function retune_pulses(s::Scenario, beams::AbstractVector{<:Beam}, enabled::AbstractVector{Bool},
                       pulses::AbstractDict; by::Symbol = :quad, cap::Real = 1.3, stay::Real = 0.0,
                       unit::Symbol = :quad,
                       target::Symbol = :reference, level::Union{Nothing,Real} = nothing,
                       wall = (), unconverted::Union{Nothing,NamedTuple} = nothing,
                       prepulse::Union{Nothing,NamedTuple} = nothing, foil::Union{Nothing,Symbol} = nothing,
                       spot_radius::Real = 0.5, on::Symbol = :auto)
    kind = _retune_kind(s, foil, on)
    kind === nothing && throw(ArgumentError(_RETUNE_NONE))
    by in (:quad, :group, :beam) || throw(ArgumentError("retune by :quad, :group or :beam, not :$by"))
    target in (:reference, :uniform) || throw(ArgumentError("the retune's target is :reference or :uniform, not :$target"))
    level === nothing || (isfinite(level) && level > 0) || throw(ArgumentError("the level must be positive"))
    cap > 0 || throw(ArgumentError("the power limit must be positive"))
    any(enabled) || throw(ArgumentError("every beam is disabled"))
    unit in PULSE_UNITS || throw(ArgumentError("a pulse table's unit is one of $PULSE_UNITS"))
    ref_by = _pulse_key(beams, pulses)
    refkey, key = Briareus.GUIAPI.beam_key(ref_by), Briareus.GUIAPI.beam_key(by)
    # each beam's share of its column, in the reference table and in the retuned one
    # (`_pulse_shares`); a retuned key gives one power to its beams, so they must share
    # alike: a cone whose quads are held some as one beam and some as four cannot
    rsh, osh = _pulse_shares(beams, ref_by, unit), _pulse_shares(beams, by, unit)
    keyshare = Dict{Any,Float64}()
    for (b, sh) in zip(beams, osh)
        get!(keyshare, key(b), sh) == sh ||
            throw(ArgumentError("$(key(b)) holds quads as one beam and as four: split or merge them all, " *
                                "or retune with one power per beam"))
    end
    nodes = sort(unique(reduce(vcat, [p.t for p in values(pulses)])))
    # the steps of a drive history of the reference (`_history_times`): the table's times,
    # the prepulse window and one step after it, gaps where nothing fires; a step without
    # 3ω light in the reference has none in the retune (its limits are zero)
    times = _history_times((; kind = :table, table = pulses, by = ref_by, times = nodes, unit), prepulse, beams,
                           fill(true, length(beams)))
    any(τ -> sum(b -> pulses[refkey(b)](τ), beams) > 0, times) ||
        throw(ArgumentError("the reference pulses have no power"))
    op = s.deposition
    keep = findall(enabled)
    kept = beams[keep]
    opk = op[keep]
    # a key's reference power: its beams' mean (the table may be finer than the retune)
    members = Dict(k => [(refkey(beams[i]), rsh[i]) for i in keep if key(beams[i]) == k] for k in unique(key.(kept)))
    nominal(k, τ) = sum(sh * pulses[r](τ) for (r, sh) in members[k]) / length(members[k])
    limits(k, τ) = (0.0, cap * nominal(k, τ))
    # the level in the quantity's units: a capsule's c_00 = √(4π) × its mean flux
    lv = level === nothing ? nothing : kind === :capsule ? sqrt(4π) * Float64(level) : Float64(level)
    lkw = (; by, limits, stay, nominal = stay > 0 ? nominal : nothing, target, level = lv)
    refbeams(τ) = _table_beams(beams, pulses, τ; by = ref_by, shares = rsh)
    oldbeams(τ) = _table_beams(kept, pulses, τ; by = ref_by, shares = rsh[keep])
    # The unconverted light of every beam (for the reference) and of the kept beams, as
    # keywords of `drive_history`; the light of the kept beams also goes to the retune,
    # with a separate `by`.
    uc = _light_kw(unconverted, prepulse, beams; by)
    uck = _light_kw(_kept_unconverted(unconverted, keep), _kept_unconverted(prepulse, keep), kept; by)
    ucr = Base.structdiff(uck, NamedTuple{(:by,)})
    if kind === :foil
        t = s.target
        es = _spot_elements(t, foil, spot_radius)
        es === nothing && throw(ArgumentError("the foil $foil has no front elements within $spot_radius mm of its center"))
        ref = (; time = times, deposition = _absorbed_light(op, refbeams, times, uc, by))
        r = retune_drive(FoilIrradiance(foil; radius = spot_radius), opk, kept, ref; lkw..., ucr...)
        newbeams(τ) = with_pulses(kept, r.pulses, τ; by)
        series(o, at, kw) = _spot_series(t, es, times, [d.element for d in _absorbed_light(o, at, times, kw, by)])
        reference, before, after = series(op, refbeams, uc), series(opk, oldbeams, uck), series(opk, newbeams, uck)
    else
        observe = _observed_parts(s)
        # the histories' wall: the walls in the scenario's target (`wall_names`), else the
        # observed surface
        names = component_names(s.target)
        walls = Tuple(n for n in wall if n in names)
        wall = isempty(walls) ? observe : walls
        face = kind === :face
        # a face's statistics per step, as the drive's (`_face_stats`)
        hkw = face ? (; record = sol -> _face_stats(s, sol)) : (; ylm = true)
        hist(o, at, kw) = drive_history(s.model, o, at, times; observe, wall, lmax = 4, hkw..., kw...)
        ref = hist(op, refbeams, (; keep_albedo = true, observed_flux = face, uc...))
        q = face ? FaceFlux(s.model, observe) : CapsuleModes(s.model; observe)
        r = retune_drive(q, opk, kept, ref; lkw..., ucr...)
        newbeams2(τ) = with_pulses(kept, r.pulses, τ; by)
        reference, before, after = ref, hist(opk, oldbeams, uck), hist(opk, newbeams2, uck)
        if face
            reference, before, after = _face_series(reference), _face_series(before), _face_series(after)
        end
    end
    # the new pulses on the table's times and the steps, off outside the steps, and off
    # for the keys of the dropped beams (so that the table keys every beam); in the
    # table's unit: a key's power per beam over its beams' share of the column
    at(g, τ) = (k = findfirst(==(τ), times); k === nothing ? 0.0 : r.power[g, k])
    ts = sort(unique([nodes; times]))
    out = Dict{Symbol,PulseShape}(Symbol(k) => PulseShape(ts, [at(g, τ) / keyshare[k] for τ in ts])
                                  for (g, k) in enumerate(r.keys))
    for k in unique(key.(beams))
        haskey(out, k) || (out[k] = PulseShape(ts, zeros(length(ts))))
    end
    return (; pulses = out, unit = by === :beam ? :beam : unit, kind, target, level, retune = r, reference, before,
            after)
end

# What a retune of scenario `s` holds (`retune_pulses`'s `kind`; `on`: `:auto`,
# `:observed` or `:foil`), or `nothing` when it has neither an observed surface under a
# radiation model nor the foil `foil`.
function _retune_kind(s::Scenario, foil, on::Symbol)
    on in (:auto, :observed, :foil) || throw(ArgumentError("retune on :auto, :observed or :foil, not :$on"))
    observed = s.model !== nothing && !isempty(s.observed) ? (_observed_flat(s) ? :face : :capsule) : nothing
    hasfoil = _has_part(s.target, foil) ? :foil : nothing
    on === :observed && return observed
    on === :foil && return hasfoil
    return something(observed, hasfoil, Some(nothing))
end

# What a retune of the drive's scenario can hold: `(; on, kind, label)` per choice (`on`,
# `retune_pulses`'s: `:observed` or `:foil`), the one `on = :auto` takes first; none
# without a scenario.
function _retune_quantities(sess::DesignSession)
    s = sess.drive.scenario[]
    s === nothing && return NamedTuple[]
    foil = active_foil(sess.target)
    out = NamedTuple[]
    obs = _retune_kind(s, foil, :observed)
    name = String(something(observed_name(sess.target), :observed))
    obs === nothing ||
        push!(out, (; on = :observed, kind = obs, label = obs === :capsule ? "$name's modes" : "$name's flux"))
    _retune_kind(s, foil, :foil) === nothing || push!(out, (; on = :foil, kind = :foil, label = "$foil's irradiance"))
    return out
end

# A retune's `level` from the dialog's text in the units of the quantity `kind`: a
# capsule's Tᵣ (eV) as its mean flux σTᵣ⁴, a face's or a foil spot's mean (TW/mm²); an
# empty text is the reference's level (`nothing`).
function _retune_level(kind::Symbol, text::AbstractString)
    isempty(strip(text)) && return nothing
    x = _parse_number(text, kind === :capsule ? "the level's Tᵣ" : "the level"; positive = true)
    return kind === :capsule ? Briareus.SIGMA_SB * x^4 : x
end

# The status that explains why a drive cannot be retuned when there is nothing to hold.
const _RETUNE_NONE = "the retune holds an observed surface's x-ray drive or a backlighter foil's laser " *
                     "irradiance: observe a part (Properties ▸ Role), or make one a backlighter foil"

# A face's statistics per step of a drive history recorded with `_face_stats`.
_finite_max(x) = maximum(filter(isfinite, x); init = 0.0)

_face_series(h) = (; time = h.time, mean = [f.mean for f in h.record], peak = [f.peak for f in h.record],
                   low = [f.low for f in h.record], rms = [f.rms for f in h.record])

# A foil spot's statistics (W/cm²; the RMS relative) at each of `times` from each time's
# absorbed laser power per element `P`, over the spot elements `es`, as `_spot_stats`.
function _spot_series(t::Target, es, times, Ps)
    A = t.elem_area[es]
    I = [P[es] ./ A for P in Ps]
    mean = [sum(i .* A) / sum(A) for i in I]
    return (; time = collect(times), mean = to_Wcm2.(mean), peak = to_Wcm2.(maximum.(I)),
            low = to_Wcm2.(minimum.(I)), rms = [m > 0 ? rms_nonuniformity(i, A) : NaN for (i, m) in zip(I, mean)])
end

# The absorbed laser power on the target of `op` from the beams `at(τ)` (3ω setup powers)
# at each of `times`: the 3ω light's `Deposition`, and with the unconverted light `kw`
# (`_light_kw`'s; empty: none, the 3ω deposition itself) each color's absorbed power
# added per element, `(; element)`, as `drive_history` deposits it (one
# `_UnconvertedLight` per call: its unit-power patterns serve every time).
function _absorbed_light(op, at, times, kw, by)
    isempty(kw) && return [deposit(op, at(τ)) for τ in times]
    ul = Briareus.GUIAPI.UnconvertedLight(kw.unconverted, get(kw, :prepulse, nothing), get(kw, :by, by), kw.colors, (;))
    return map(times) do τ
        bs = at(τ)
        E = copy(deposit(op, bs).element)
        u = Briareus.GUIAPI.unconverted_at!(ul, op.target, bs, τ)
        for c in ul.colors
            E .+= u[c.name].element
        end
        (; element = E)
    end
end

# A retune's face or foil statistics `f` (`retune_pulses`'s `before`, `after`): the mean's
# relative departure at each step from the reference's mean `ref.mean`, or with a
# `level` (TW/mm²) from it where the retune solved (`solved`; NaN elsewhere).
_mean_departure(f, level, solved, ref) =
    level === nothing ? f.mean ./ ref.mean .- 1 : [s ? m / to_Wcm2(level) - 1 : NaN for (m, s) in zip(f.mean, solved)]

# The retune's level as text: the reference's, or the set one.
_level_text(v) = v.level === nothing ? "the reference's level" :
                 v.kind === :capsule ? @sprintf("Tᵣ %.1f eV", radiation_temperature(v.level)) :
                 @sprintf("%.3g TW/mm²", v.level)

# The unconverted light `u` (`retune_pulses`'s `unconverted`) of the beams `keep`: each
# color's fraction per beam, if it has one, at those beams.
_kept_unconverted(::Nothing, keep) = nothing
_kept_unconverted(u::NamedTuple, keep) = map(f -> f isa AbstractVector ? f[keep] : f, u)

"""
    retune_to_reference!(sess::DesignSession, table, out; by = :quad, cap = 1.3, stay = 0.0,
                         unit = :quad, target = :reference, level = nothing, on = :auto, load = true)

Retune the enabled beams to the reference shot whose pulses are the table `table`
(`read_pulse_table`: a `time_ns` column and one per group, quad or beam; a cone's or a
quad's column the power per `unit`, [`PULSE_UNITS`](@ref), as the table written), as
[`retune_pulses`](@ref) with the drive's scenario (which must be current), its spot
radius and the active foil ([`active_foil`](@ref)), to the reference's drive shape
(`target = :reference`) or a uniform drive (`:uniform`), at the reference's level or at
`level` (TW/mm², `retune_pulses`'s), and write the new pulses to `out`
(`write_pulse_table`), on a worker when `spawn[]` is set. With `load` the beams then fire
the new pulses: they become the Time panel's (as [`load_pulses!`](@ref) of `out`, before
the CSV's rounding), an edit that undo takes back. The disabled beams are the ones
the shot loses. The retune holds the observed surface's drive (a capsule's modes, a flat
face's flux) when the drive has a radiation model, else the active foil's irradiance;
`on = :observed` or `:foil` chooses (`retune_pulses`'s `on`) when the design has both.
The unconverted 1ω and 2ω light takes part as in the drive
([`set_unconverted_ratio!`](@ref); `retune_pulses`'s `unconverted`), with the prepulse
([`set_prepulse!`](@ref); `retune_pulses`'s `prepulse`), in the reference shot and in
the retune, on the beams whose 1ω or 2ω light takes part ([`participates`](@ref)), for a
foil as each color's absorbed irradiance on its spot. Progress and the outcome go to
`drive.retune_running`/`retune_status`: for a capsule the worst-time mode 1 before and
after, for a face or a foil the worst RMS nonuniformity and the worst departure of the
mean from the reference's (from `level` when it is set; for a capsule the Tᵣ's likewise).
A drive with neither is refused with a status message.
"""
function retune_to_reference!(sess::DesignSession, table::AbstractString, out::AbstractString;
                              by::Symbol = :quad, cap::Real = 1.3, stay::Real = 0.0, unit::Symbol = :quad,
                              target::Symbol = :reference, level::Union{Nothing,Real} = nothing,
                              on::Symbol = :auto, load::Bool = true)
    dr = sess.drive
    dr.retune_running[] && return sess
    if !_drive_current(sess)
        dr.retune_status[] = "compute the drive first: the retune reuses its operator and radiation model"
        return sess
    end
    s = dr.scenario[]
    foil = active_foil(sess.target)
    kind = _retune_kind(s, foil, on)
    kind === nothing && (dr.retune_status[] = on === :auto ? _RETUNE_NONE :
                                              "the design has no $(on === :foil ? "foil" : "observed surface") to hold"; return sess)
    # the unconverted light of every beam (the reference fires them all)
    unconverted = _history_unconverted(sess, eachindex(sess.beams.beams[]))
    prepulse = _history_prepulse(sess, eachindex(sess.beams.beams[]))
    pulses = read_pulse_table(table)
    beams = sess.beams.beams[]
    enabled = copy(sess.beams.enabled[])
    out = String(out)
    spot_radius = dr.spot_radius[]
    dr.retune_running[] = true
    dr.retune_status[] = "retuning to $(basename(table))…"
    runner = _runner(sess.spawn, sess.deliver)
    wall = Tuple(wall_names(sess.target))
    runner(() -> retune_pulses(s, beams, enabled, pulses; by, cap, stay, unit, target, level, wall, unconverted,
                               prepulse, foil, spot_radius, on), function (res)
        dr.retune_running[] = false
        if res.err !== nothing
            dr.retune_status[] = "retune failed: " * _errmsg(res.err)
            return
        end
        v = res.value
        write_pulse_table(out, v.pulses)
        # the time panel plots the retune
        bt = sess.beams
        nodes = sort(unique(reduce(vcat, [p.t for p in values(pulses)])))
        sess.time.retune[] = (; label = basename(table), v.kind, v.target, v.level,
                              solved = isfinite.(v.retune.residual), by, unit, pulses = v.pulses, reference_pulses = pulses,
                              reference_by = _pulse_key(beams, pulses), nodes, beams, enabled,
                              table = (; cluster = copy(bt.cluster[]), colors = cluster_colors(bt)),
                              before = v.before, after = v.after, reference = v.reference)
        # the beams fire the new pulses (the table's beams unchanged while the worker ran)
        loaded = load && length(bt.beams[]) == length(beams) &&
                 all(b -> haskey(v.pulses, Briareus.GUIAPI.beam_key(by)(b)), bt.beams[])
        loaded && _set_pulses!(sess, _table_pulses(v.pulses, by, basename(out); unit))
        set_chart!(sess.time, :retune)
        n = length(v.retune.keys)
        aim = v.target === :uniform ? "uniform, " : ""
        solved = isfinite.(v.retune.residual)
        if v.kind === :capsule
            ε1(h) = 100 * maximum(mode_amplitude(h.ylm[:, k], 1) for k in eachindex(h.time))
            # the Tᵣ against the reference's, or against the level's where the retune solves
            ΔT(h) = v.level === nothing ? maximum(abs.(h.Tr .- v.reference.Tr)) :
                    _finite_max(abs.(h.Tr[solved] .- radiation_temperature(v.level)))
            dr.retune_status[] = @sprintf("%d %s powers (%s%s): worst ε₁ %.2f%% → %.2f%%, |ΔTᵣ| %.1f → %.1f eV; wrote %s",
                                          n, by, aim, _level_text(v), ε1(v.before), ε1(v.after), ΔT(v.before),
                                          ΔT(v.after), basename(out))
        else
            # (over the steps where the reference lights the face or the spot)
            rms(f) = 100 * _finite_max(f.rms)
            Δ(f) = 100 * _finite_max(abs.(_mean_departure(f, v.level, solved, v.reference)))
            dr.retune_status[] = @sprintf("%d %s powers (%s%s): %s worst RMS %.2f%% → %.2f%% (reference %.2f%%), |Δmean| %.2f%% → %.2f%%; wrote %s",
                                          n, by, aim, _level_text(v),
                                          v.kind === :face ? "face" : "foil spot", rms(v.before), rms(v.after),
                                          rms(v.reference), Δ(v.before), Δ(v.after), basename(out))
        end
        sess.status[] = "wrote the retuned pulses to " * basename(out) *
                        (loaded ? "; the beams fire them (Undo takes them back)" : "")
    end)
    return sess
end
