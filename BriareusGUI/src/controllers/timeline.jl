# Time in the window: the pulses the beams fire (a pulse table, or a pulse shape scaled
# to the table's powers), the drive history they give on the drive's scenario (capsule
# Tᵣ, its Legendre modes, the synthetic diagnostics ticked in the list of diagnostic
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
# which are on, and the pulses.
_history_key(sess::DesignSession) = (sess.drive.scenario[], sess.beams.beams[], sess.beams.enabled[],
                                     sess.time.pulses[])

"""
    history_stale(sess::DesignSession) -> Bool

Whether the drive history was computed for another scenario, other beams or on/off
states, or other pulses (`false` without a history).
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
(`energy_budget`). It runs on a worker when `spawn[]` is set, with progress per step and
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
        (tm.status[] = "the drive history follows the capsule: a hohlraum with a capsule, or a sphere part " *
                       "with the role Capsule"; return sess)
    beams, enabled = sess.beams.beams[], sess.beams.enabled[]
    times = filter(τ -> sum(b -> b.power, _pulse_beams(src, beams, enabled, τ)) > 0, src.times)
    isempty(times) && (tm.status[] = "the pulses fire no enabled beam"; return sess)
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
    runner(() -> _history_job(s, beams, enabled, src, times; diags, progress, cancel, key), function (out)
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
function _history_job(s::Scenario, beams, enabled, src, times; diags = NamedTuple[],
                        progress = k -> nothing, cancel = Threads.Atomic{Bool}(false), key = nothing)
    t0 = time()
    t = s.target
    names = component_names(t)
    observe = Tuple(unique([t.components[t.elem_comp[e]].part for e in s.observed]))
    wall = :wall in names ? (:wall,) : observe
    views = _diagnostic_views(diags, names)
    A = t.elem_area
    aperture = [t.components[t.elem_comp[e]].kind === :aperture for e in 1:nelements(t)]
    steps = Ref(0)
    record = function (sol)
        cancel[] && throw(_HistoryCanceled())
        laser = Float32.(to_Wcm2.(sol.laser ./ A))
        laser[aperture] .= NaN32
        diag = [_record_diagnostic!(v, sol) for v in views]
        progress(steps[] += 1)
        return (; laser, radiosity = Float32.(to_Wcm2.(sol.B)), incident = Float32.(to_Wcm2.(sol.H)), diag)
    end
    h = drive_history(s.model, s.deposition, τ -> _pulse_beams(src, beams, enabled, τ), times;
                      observe, wall, lmax = 4, record)
    fields = Dict{Symbol,Matrix{Float32}}(f => reduce(hcat, [r[f] for r in h.record]) for f in DRIVE_FIELDS)
    ranges = Dict{Symbol,Tuple{Float64,Float64}}()
    for (f, M) in fields
        fin = filter(isfinite, M)
        lo, hi = isempty(fin) ? (0.0, 1.0) : Float64.(extrema(fin))
        f === :laser && (lo = 0.0)
        ranges[f] = (lo, hi > lo ? hi : lo + 1)
    end
    legendre = h.moments[2:5, :] ./ h.moments[1:1, :]
    recorded = NamedTuple[(; v.id, v.name, Tr = [r.diag[j][1] for r in h.record],
                           intensity = [r.diag[j][2] for r in h.record], note = v.note[])
                          for (j, v) in enumerate(views)]
    return DriveHistory(h.time, h.Tr, legendre, recorded, _energy_series(h, t),
                        fields, ranges, h, key, time() - t0)
end

# Where the laser energy has gone by each time of a history (kJ, `energy_budget` up to
# that time): the laser delivered, the heat left in the walls (every surface but the
# capsule: unconverted laser and absorbed x-rays), in the capsule, the x-rays out of
# the entrance holes, and out of windows (other apertures and thin foils).
function _energy_series(h, t::Target)
    holes = Set(c.part for c in t.components if c.kind === :aperture && !(c.surface isa CylinderSurface))
    capsule = Set(c.part for c in t.components if c.name === :capsule || c.part === :capsule)
    eb = [energy_budget(h; to = τ) for τ in h.time]
    sumof(d, keep) = sum((v for (k, v) in d if keep(k)); init = 0.0)
    series = NamedTuple[
        (; name = :laser, values = [e.setup for e in eb]),
        (; name = :wall, values = [sumof(e.heating, k -> !(k in capsule)) for e in eb]),
        (; name = :capsule, values = [sumof(e.heating, in(capsule)) for e in eb]),
        (; name = :holes, values = [sumof(e.escaped, in(holes)) for e in eb]),
        (; name = :windows, values = [sumof(e.escaped, k -> !(k in holes)) + sumof(e.transmitted, _ -> true)
                                      for e in eb])]
    # (a hohlraum without windows, or a target without a capsule, has no such series)
    return filter(s -> s.name === :laser || any(>(1e-9), s.values), series)
end

const ENERGY_LABELS = Dict(:laser => "Laser delivered", :wall => "Walls (heat)", :capsule => "Capsule (heat)",
                           :holes => "Out of the entrance holes", :windows => "Out of windows")

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

_series(label, x, y, color; style = :solid) =
    (; label = String(label), x = collect(Float64, x), y = collect(Float64, y), color, style)

"""
    time_chart(sess::DesignSession, chart = sess.time.chart[]) -> NamedTuple

What the time panel plots: `(; chart, title, xlabel, left, right, legend, cursor,
empty)`. `left` and `right` are the y axes, `(; label, series)` (`right` is `nothing`
for one axis), each series `(; label, x, y, color, style)` with `style` `:solid`,
`:dash` or `:dot`; `legend` the key, `(; label, color, style)` rows; `cursor` the time
of the cursor's history step (`NaN` without); `empty` a message when there is nothing to
plot (then no series).

- `:pulses`: the laser power of each beam group (the sum over its enabled beams, TW;
  groups colored as in the model tree) and the total.
- `:drive`: the capsule's Tr and each recorded diagnostic's (eV, dashed), and a₂/a₀ and
  a₄/a₀ (%).
- `:diagnostics`: each recorded diagnostic's Tr (eV) and the intensity toward it
  (TW/sr, dashed).
- `:energy`: the laser energy delivered by each time and where it went (kJ).
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
        left = [_series("Capsule Tr", h.time, h.Tr, _INK)]
        for (j, d) in enumerate(h.diagnostics)
            any(isfinite, d.Tr) && push!(left, _series(d.name, h.time, d.Tr, _diagnostic_color(j); style = :dash))
        end
        right = [_series("a₂/a₀", h.time, 100 .* h.legendre[2, :], "#0072B2"),
                 _series("a₄/a₀", h.time, 100 .* h.legendre[4, :], "#009E73")]
        return (; chart, title = "Capsule drive history" * stale, xlabel = "t (ns)",
                left = (; label = "Tr (eV)", series = left), right = (; label = "aₗ/a₀ (%)", series = right),
                legend = [_key(s) for s in vcat(left, right)], cursor, empty = "")
    elseif chart === :diagnostics
        h === nothing && return none("Diagnostics", "Tick diagnostics below, then Run the drive history")
        isempty(h.diagnostics) && return none("Diagnostics", "No diagnostic was recorded: tick one, then Run the drive history")
        left = [_series(d.name * " Tr", h.time, d.Tr, _diagnostic_color(j))
                for (j, d) in enumerate(h.diagnostics) if any(isfinite, d.Tr)]
        right = [_series(d.name * " intensity", h.time, d.intensity, _diagnostic_color(j); style = :dash)
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
        color = Dict(:laser => _INK, :wall => get(cc, :wall, palette_color(1)),
                     :capsule => get(cc, :capsule, palette_color(4)), :holes => _GRAY, :windows => "#D55E00")
        series = [_series(ENERGY_LABELS[s.name], h.time, s.values, color[s.name];
                          style = s.name === :laser ? :dash : :solid) for s in h.energy]
        return (; chart, title = "Energy by each time" * stale, xlabel = "t (ns)",
                left = (; label = "Energy (kJ)", series), right = nothing,
                legend = [_key(s) for s in series], cursor, empty = "")
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
# nothing in a history: a kind other than the soft x-ray power one, or a view through an
# aperture the target lacks.
function _diagnostic_views(diags, names)
    out = NamedTuple[]
    for d in diags
        view, note = nothing, d.note
        if d.def !== nothing
            if d.def.kind !== :xray_power
                note = "a $(d.def.kind) diagnostic is not recorded in a drive history"
            else
                v = diagnostic(d.def)
                v.through === nothing || v.through in names ? (view = v) :
                    (note = "the target has no $(v.through) to look through")
            end
        end
        push!(out, (; d.id, d.name, view, note = Ref(note)))
    end
    return out
end

# One step of a recorded diagnostic: (Tr in eV, intensity in TW/sr); NaN, with the
# reason kept, where the view gives none.
function _record_diagnostic!(v, sol)
    v.view === nothing && return (NaN, NaN)
    try
        r = dante(v.view, sol)
        return (Float64(r.Tr), Float64(r.intensity))
    catch err
        err isa ArgumentError || rethrow()
        isempty(v.note[]) && (v.note[] = sprint(showerror, err))
        return (NaN, NaN)
    end
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
