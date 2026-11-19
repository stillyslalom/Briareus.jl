# Several phases of one pulse under one pointing: a scenario per phase (its own target,
# radiation model, deposition settings, observation and laser–plasma model) and beam
# powers per phase, the aims of one pointing parameterization shared by every phase, and
# a power scale per (group, phase) beside them; one parameter vector and one evaluation
# per phase for `optimize_pointing`.

"""
    PerPhase{V}
    PerPhase(:picket => a, :peak => b, ...)

One value per phase of a [`PhasedScenario`](@ref), in the phases' order: the beams
(`apply` of a [`PhasedPointing`](@ref), `PointingResult.beams`), the evaluations
([`evaluate`](@ref)) or the summaries ([`pointing_summary`](@ref)) of a phased design.
Indexed by the phase's name (`v[:picket]`) or position (`v[1]`); `keys` gives the names,
`values` the values, `pairs` the `name => value` pairs, and iteration, `length`,
`haskey` and `map` (a function applied to each value, the names kept) work as for an
ordered dictionary whose iteration yields the values.
"""
struct PerPhase{V}
    names::Vector{Symbol}
    values::Vector{V}
    function PerPhase(names::Vector{Symbol}, values::Vector{V}) where {V}
        length(names) == length(values) || throw(DimensionMismatch("$(length(names)) names for $(length(values)) values"))
        allunique(names) || throw(ArgumentError("phase names must be unique, got $(names)"))
        return new{V}(names, values)
    end
end

PerPhase(ps::Pair{Symbol}...) = PerPhase(Symbol[first(p) for p in ps], [last(p) for p in ps])
PerPhase(itr) = PerPhase(collect(itr)...)

Base.keys(v::PerPhase) = v.names
Base.values(v::PerPhase) = v.values
Base.pairs(v::PerPhase) = (n => x for (n, x) in zip(v.names, v.values))
Base.length(v::PerPhase) = length(v.values)
Base.iterate(v::PerPhase, i...) = iterate(v.values, i...)
Base.eltype(::Type{PerPhase{V}}) where {V} = V
Base.haskey(v::PerPhase, k::Symbol) = k in v.names
Base.getindex(v::PerPhase, i::Integer) = v.values[i]
function Base.getindex(v::PerPhase, k::Symbol)
    i = findfirst(==(k), v.names)
    i === nothing && throw(KeyError(k))
    return v.values[i]
end
Base.map(f, v::PerPhase) = PerPhase(v.names, map(f, v.values))
Base.:(==)(a::PerPhase, b::PerPhase) = a.names == b.names && a.values == b.values

# A value in a `PerPhase`'s display: a summary of arrays and evaluations, else compact.
_show_phase_value(io::IO, x) = x isa Union{AbstractArray,Evaluation} ? print(io, summary(x)) :
                               show(IOContext(io, :compact => true, :limit => true), x)

function Base.show(io::IO, v::PerPhase)
    print(io, "PerPhase(")
    for (i, (n, x)) in enumerate(pairs(v))
        i > 1 && print(io, ", ")
        print(io, ":", n, " => ")
        _show_phase_value(io, x)
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", v::PerPhase)
    print(io, "PerPhase with ", length(v), length(v) == 1 ? " phase:" : " phases:")
    w = maximum(n -> length(string(n)), v.names; init = 0)
    for (n, x) in pairs(v)
        print(io, "\n  ", rpad(string(n), w), "  ")
        _show_phase_value(io, x)
    end
end

"""
    Phase(name, target; power = nothing, pulses = nothing, time = nothing, pulse_by = :group,
          scaled = true, history = nothing, kwargs...)
    Phase(name, scenario::Scenario; power = nothing, pulses = nothing, time = scenario.time,
          pulse_by = :group, scaled = true, history = nothing)

One phase of a [`PhasedScenario`](@ref): its name, what its beams deposit in and radiate
from, and the beams' powers in it. With a `target`, the phased scenario builds the
phase's [`Scenario`](@ref) for the shared pointing, with the phased scenario's
`Scenario` keywords and this phase's `kwargs` over them (`observe`, `deposition_kw`,
`viewfactor_kw`, `lpi`, `albedo`, `model`, `radiation`, `time`); a phase with its own
wall albedos is a target with its own materials, which can share the view factors of
another phase's model through `viewfactor_kw = (; F = model.F)`. A ready `scenario` is
used with its operators, radiation model and options, its parameterization replaced by
the shared pointing (its deposition operator must serve the pointing's beams, as for
`remake(scenario; parameterization)`). One target, view-factor matrix and deposition
operator can serve every phase, windows intact or burned through, as `remake`s of a
shared scenario with each phase's albedos and time. `history`, a [`drive_history`](@ref)
of the phase's target run with `keep_albedo = true`, gives them at the phase's `time`:
`Phase(:rise, shared; history = h, time = τ, pulses)` is `Phase(:rise, remake(shared;
history = h, time = τ); pulses)`, the albedos of [`history_state`](@ref)`(target, h, τ)`
(with a `target`, they are the built scenario's `albedo`). A reference scenario at the
same time is `remake(reference_scenario; history = reference_history, time = τ)`, so
that the reference checks of [`mode_gate`](@ref) and [`drive_floor`](@ref) find the
phase's time.

The phase's beam powers (TW per beam) replace the pointing's:

- `power`: a number for every beam, a vector (one per beam of the pointing, in its
  order) or a function of a beam; `nothing` keeps the powers the pointing gives (its
  template's, or a `:power` field's).
- `pulses` at `time` (ns): each beam's power is its pulse's at that time, as
  [`with_pulses`](@ref)`(beams, pulses, time; by = pulse_by)` sets it.

`time` is also the built scenario's time (the materials' albedos and the laser–plasma
model at that time) unless `kwargs` give another; a ready scenario keeps its own.
`scaled = false` fixes the phase's powers (a peak at the laser's limit): the phase has
no power scales among the parameters.
"""
struct Phase
    name::Symbol
    source::Union{Target,Scenario}
    power::Any
    pulses::Any
    time::Union{Nothing,Float64}
    pulse_by::Any
    scaled::Bool
    kwargs::NamedTuple
end

function Phase(name, source::Union{Target,Scenario}; power = nothing, pulses = nothing,
               time::Union{Nothing,Real} = nothing, pulse_by = :group, scaled::Bool = true,
               history = nothing, kwargs...)
    power === nothing || pulses === nothing ||
        throw(ArgumentError("phase $name: give `power` or `pulses`, not both"))
    source isa Scenario && !isempty(kwargs) &&
        throw(ArgumentError("phase $name: a ready scenario takes no Scenario keywords ($(join(keys(kwargs), ", "))); " *
                            "remake the scenario instead"))
    if history !== nothing
        time === nothing && throw(ArgumentError("phase $name: a `history` needs the phase's `time`"))
        haskey(kwargs, :albedo) && throw(ArgumentError("phase $name: give the albedos by `history` or by `albedo`, not both"))
        if source isa Scenario
            source = remake(source; history, time)
        else
            st = history_state(source, history, time)
            kwargs = pairs(merge(values(kwargs), (; st.albedo)))
        end
    end
    tm = time === nothing ? (source isa Scenario ? source.time : nothing) : Float64(time)
    pulses === nothing || tm !== nothing ||
        throw(ArgumentError("phase $name: pulses need the phase's `time`"))
    return Phase(Symbol(name), source, power, pulses, tm, pulse_by, scaled, values(kwargs))
end

Base.show(io::IO, ph::Phase) =
    print(io, "Phase(:", ph.name, ph.time === nothing ? "" : string(", ", ph.time, " ns"),
          ph.scaled ? "" : ", fixed powers", ")")

# The phase's powers per beam of `beams` (the pointing's at its initial parameters), or
# `nothing` for the pointing's own.
function _phase_power(ph::Phase, beams::AbstractVector{<:Beam})
    if ph.pulses !== nothing
        return Float64[ForwardDiff.value(b.power) for b in with_pulses(beams, ph.pulses, ph.time; by = ph.pulse_by)]
    end
    P = ph.power
    P === nothing && return nothing
    P isa Real && return fill(Float64(P), length(beams))
    if P isa AbstractVector
        length(P) == length(beams) ||
            throw(DimensionMismatch("phase $(ph.name): $(length(P)) powers for $(length(beams)) beams"))
        return collect(Float64, P)
    end
    return Float64[P(b) for b in beams]
end

"""
    PhasedPointing(pointing, phases, power; scaled = phases, by = :group)

The parameterization of a [`PhasedScenario`](@ref) (its `param`): the parameters of
`pointing` (any parameterization with `apply` and `initial_parameters`: a
[`GroupPointing`](@ref) by cone, quad or beam, `relative` or not, a
[`FramePointing`](@ref), a [`HeldParameters`](@ref) of either), which set the aims of
every phase, followed by a power scale per (group, phase) for the phases in `scaled`.
`phases` are the phases' names in order, and `power[k]` the powers of phase `k` (a
vector per beam of the pointing, or `nothing` for the pointing's own). Groups are keyed
by `by` on the pointing's beams as for [`PowerScales`](@ref) (`:group`, `:quad`,
`:beam` or a function of a beam), so a function that gives several quads one key ties
their scales (one pulse shape for them, a multiplier per phase). The scales' keys
`keys` are `(key, phase)` tuples, phase by phase in order of first appearance, and the
scales are named `"<key>.power_scale.<phase>"` after the pointing's parameter names.

`apply(p, x)` gives a [`PerPhase`](@ref) of beam vectors: the pointing's beams at
`x[1:n]`, in each phase with the phase's power (times its group's scale in a scaled
phase). Generic in the number type.
"""
struct PhasedPointing{P,K}
    pointing::P
    n::Int                                          # parameters of `pointing`
    phases::Vector{Symbol}
    power::Vector{Union{Nothing,Vector{Float64}}}   # per phase: per beam, or the pointing's own
    keys::Vector{Tuple{K,Symbol}}                   # the scales' (group key, phase)
    group::Vector{Vector{Int}}                      # per phase: per beam, its scale (empty: fixed)
    by::Any
end

function PhasedPointing(pointing, phases::AbstractVector{Symbol}, power::AbstractVector;
                        scaled = phases, by = :group)
    length(power) == length(phases) || throw(DimensionMismatch("$(length(power)) powers for $(length(phases)) phases"))
    allunique(phases) || throw(ArgumentError("phase names must be unique, got $(phases)"))
    for k in scaled
        k in phases || throw(ArgumentError("no phase $(repr(k)) to scale; the phases are $(join(repr.(phases), ", "))"))
    end
    x0 = initial_parameters(pointing)
    beams = apply(pointing, x0)
    pw = Union{Nothing,Vector{Float64}}[P === nothing ? nothing : collect(Float64, P) for P in power]
    for (k, P) in zip(phases, pw)
        P === nothing || length(P) == length(beams) ||
            throw(DimensionMismatch("phase $k: $(length(P)) powers for $(length(beams)) beams"))
    end
    by === :quad && any(b -> isempty(b.quad), beams) &&
        throw(ArgumentError("by = :quad needs quad names on every beam (OMEGA beams have none)"))
    key = by === :group ? (b -> b.group) : by === :quad ? (b -> b.quad) : by === :beam ? (b -> b.name) : by
    ks = map(key, beams)
    K = eltype(ks)
    sk = Tuple{K,Symbol}[]
    index = Dict{Tuple{K,Symbol},Int}()
    group = Vector{Int}[]
    for ph in phases
        if ph in scaled
            g = Vector{Int}(undef, length(beams))
            for i in eachindex(beams)
                kk = (ks[i], ph)
                g[i] = get!(index, kk) do
                    push!(sk, kk)
                    length(sk)
                end
            end
            push!(group, g)
        else
            push!(group, Int[])
        end
    end
    return PhasedPointing{typeof(pointing),K}(pointing, length(x0), collect(Symbol, phases), pw, sk, group, by)
end

"Number of parameters of a `PhasedPointing`: the pointing's plus one per power scale."
nparams(p::PhasedPointing) = p.n + length(p.keys)

"""
Names of a `PhasedPointing`'s parameters: the pointing's (`"x[i]"` for one without
names), then `"<key>.power_scale.<phase>"` per scale.
"""
parameter_names(p::PhasedPointing) =
    vcat(_inner_names(p.pointing, p.n), [string(k, ".power_scale.", ph) for (k, ph) in p.keys])

"""
    initial_parameters(p::PhasedPointing) -> Vector{Float64}

The pointing's initial parameters, then a scale of 1 per (group, phase).
"""
initial_parameters(p::PhasedPointing) = vcat(initial_parameters(p.pointing), ones(length(p.keys)))

_nparams_shown(p::PhasedPointing) = nparams(p)

"""
    apply(p::PhasedPointing, x) -> PerPhase

The beams of every phase at `x`: the pointing's beams at `x[1:n]` (one set of aims), in
each phase with its powers, each times its group's scale in `x[n+1:end]` in a scaled
phase (element type follows `x`).
"""
function apply(p::PhasedPointing, x::AbstractVector{T}) where {T}
    length(x) == nparams(p) || throw(DimensionMismatch("expected $(nparams(p)) parameters"))
    beams = apply(p.pointing, x[1:p.n])
    return PerPhase(p.phases, [_phase_beams(beams, p.power[k], p.group[k], x, p.n) for k in eachindex(p.phases)])
end

# The pointing's beams with phase powers `P` (`nothing`: their own), each times its
# scale `x[n + group[i]]` (no `group`: a fixed phase). A function barrier: the beams'
# type is known here.
function _phase_beams(beams::Vector{Beam{T}}, P::Union{Nothing,Vector{Float64}}, group::Vector{Int},
                      x::AbstractVector, n::Int) where {T}
    if isempty(group)
        P === nothing && return beams
        return Beam{T}[remake(b; power = T(P[i])) for (i, b) in enumerate(beams)]
    end
    P === nothing && return Beam{T}[remake(b; power = b.power * x[n+group[i]]) for (i, b) in enumerate(beams)]
    return Beam{T}[remake(b; power = P[i] * x[n+group[i]]) for (i, b) in enumerate(beams)]
end

"""
    phase_beams(p, beams; scale = 1) -> PerPhase

The beams `beams` (the pointing's beams, in its order: their aims, spots and focus) in
every phase of `p` (a [`PhasedScenario`](@ref) or its [`PhasedPointing`](@ref)) at the
phase's powers times `scale` in the scaled phases (a number, or a function `(beam,
phase) -> scale`): a reference pointing evaluated over the same phases,
`evaluate(scenario, phase_beams(scenario, reference_beams))`. A phase whose powers are
the pointing's own keeps the beams' powers (times the scale).
"""
function phase_beams(p::PhasedPointing, beams::AbstractVector{<:Beam}; scale = 1)
    tmpl = apply(p.pointing, initial_parameters(p.pointing))
    length(beams) == length(tmpl) && all(i -> beams[i].name == tmpl[i].name, eachindex(beams)) ||
        throw(ArgumentError("phase_beams: the beams are not the pointing's (its $(length(tmpl)) beams, by name, in order)"))
    sc(b, ph) = scale isa Real ? scale : scale(b, ph)
    vals = map(eachindex(p.phases)) do k
        ph, P, s = p.phases[k], p.power[k], !isempty(p.group[k])
        [remake(b; power = (P === nothing ? b.power : P[i]) * (s ? sc(b, ph) : 1)) for (i, b) in enumerate(beams)]
    end
    return PerPhase(p.phases, vals)
end

"""
    power_scales(p, ev::PerPhase) -> Vector
    power_scales(p, beams::PerPhase) -> Vector
    power_scales(p, x::AbstractVector) -> Vector

The power scales of a phased design, one per `(key, phase)` of `p.param.keys` in that
order (`p` a [`PhasedScenario`](@ref) or its [`PhasedPointing`](@ref)), the parameters
`x[n+1:end]` that `apply` multiplies the phases' powers by. From a [`PerPhase`](@ref) of
evaluations (their setup beams, `ev[phase].beams`) or of beam vectors (in the pointing's
order, as `apply` and `evaluate` give them), each scale is a beam's power over its
phase's power ([`phase_powers`](@ref)), for the first beam of the key with positive
phase power; a key whose beams are all at zero power in the phase has no scale to read
and gives 1, as [`parameters`](@ref) does. A scaled phase whose powers are the
pointing's own (no `power` or `pulses` in its [`Phase`](@ref)) has no phase power to
divide by: an `ArgumentError` (read the scales from `x`).

Generic in the number type and thread-safe: inside an objective of the `PerPhase`
evaluation the scales carry the optimizer's duals, e.g. a regularization of the scales
toward `s0` (a vector in the order of the keys),
`ev -> J(ev) + w * sum(abs2, power_scales(ps, ev) ./ s0 .- 1) / length(s0)`.
[`power_scale_table`](@ref) gives them per key and phase for a report.
"""
function power_scales(p::PhasedPointing, x::AbstractVector{<:Real})
    length(x) == nparams(p) || throw(DimensionMismatch("expected $(nparams(p)) parameters"))
    return x[p.n+1:end]
end
power_scales(p::PhasedPointing, ev::PerPhase{<:Evaluation}) =
    _power_scales(p, ev.names, map(e -> e.beams, ev.values))
power_scales(p::PhasedPointing, beams::PerPhase{<:AbstractVector{<:Beam}}) =
    _power_scales(p, beams.names, beams.values)

function _power_scales(p::PhasedPointing, names::Vector{Symbol}, bs::AbstractVector{<:AbstractVector{<:Beam}})
    names == p.phases ||
        throw(ArgumentError("the phases ($(join(names, ", "))) are not the scenario's ($(join(p.phases, ", ")))"))
    return _read_scales(p, bs)
end

# A function barrier: the beams' number type is known here.
function _read_scales(p::PhasedPointing, bs::AbstractVector{<:AbstractVector{Beam{T}}}) where {T}
    s = ones(T, length(p.keys))
    done = falses(length(p.keys))
    for k in eachindex(p.phases)
        g = p.group[k]
        isempty(g) && continue
        P = p.power[k]
        P === nothing &&
            throw(ArgumentError("phase $(repr(p.phases[k])) is scaled on the pointing's own powers: " *
                                "its scales are read from the parameters, `power_scales(scenario, x)`"))
        b = bs[k]
        length(b) == length(g) ||
            throw(DimensionMismatch("phase $(repr(p.phases[k])): $(length(b)) beams for the pointing's $(length(g))"))
        for i in eachindex(g)
            j = g[i]
            (done[j] || !(P[i] > 0)) && continue
            s[j] = b[i].power / P[i]
            done[j] = true
        end
    end
    return s
end

"""
    phase_powers(p, phase) -> Union{Vector{Float64},Nothing}

The powers (TW per beam of the pointing, in its order) of the phase named `phase` of `p`
(a [`PhasedScenario`](@ref) or its [`PhasedPointing`](@ref)) before its power scales:
the [`Phase`](@ref)'s `power`, or its `pulses` at its time; `nothing` for a phase whose
powers are the pointing's own. A copy.
"""
function phase_powers(p::PhasedPointing, phase::Symbol)
    k = findfirst(==(phase), p.phases)
    k === nothing && throw(KeyError(phase))
    P = p.power[k]
    return P === nothing ? nothing : copy(P)
end

"""
    PowerScaleTable{K}

The power scales of a phased design per key and phase, from
[`power_scale_table`](@ref): a row per key (`keys`, the scales' group keys in order of
first appearance) and a column per scaled phase (`phases`), with `scales[i, j]` the
scale of `keys[i]` in `phases[j]` (`NaN` where that key has no scale in the phase). The
text display prints the table and each phase's least, mean and largest scale, for a
study's log; the fields serve other reports.
"""
struct PowerScaleTable{K}
    keys::Vector{K}
    phases::Vector{Symbol}
    scales::Matrix{Float64}
end

"""
    power_scale_table(p, result::PointingResult) -> PowerScaleTable
    power_scale_table(p, x::AbstractVector) / power_scale_table(p, beams_or_ev::PerPhase)

The power scales of a phased design (`p` a [`PhasedScenario`](@ref) or its
[`PhasedPointing`](@ref); the design a result, a parameter vector, or the evaluations
or beams that [`power_scales`](@ref) reads) as a [`PowerScaleTable`](@ref), in Float64
(the values of duals).
"""
power_scale_table(p::PhasedPointing, r::PointingResult) = power_scale_table(p, r.x)
function power_scale_table(p::PhasedPointing{<:Any,K}, v) where {K}
    s = Float64.(ForwardDiff.value.(power_scales(p, v)))
    keys = unique(first.(p.keys))
    phases = unique(last.(p.keys))
    S = fill(NaN, length(keys), length(phases))
    for (g, (k, ph)) in enumerate(p.keys)
        S[findfirst(isequal(k), keys), findfirst(==(ph), phases)] = s[g]
    end
    return PowerScaleTable{K}(keys, phases, S)
end

Base.show(io::IO, t::PowerScaleTable) =
    print(io, "PowerScaleTable(", length(t.keys), " keys × ", length(t.phases), " phases (", join(t.phases, ", "), "))")

function Base.show(io::IO, ::MIME"text/plain", t::PowerScaleTable)
    print(io, "Power scales of ", length(t.keys), " keys in ", length(t.phases), length(t.phases) == 1 ? " phase:" : " phases:")
    w = max(maximum(k -> length(string(k)), t.keys; init = 0), 4)
    cw = [max(length(string(ph)), 7) for ph in t.phases]
    print(io, "\n  ", rpad("", w))
    for (j, ph) in enumerate(t.phases)
        print(io, "  ", lpad(string(ph), cw[j]))
    end
    row(label, vals) = (print(io, "\n  ", rpad(label, w));
                        foreach(j -> print(io, "  ", lpad(isnan(vals[j]) ? "–" : @sprintf("%.4f", vals[j]), cw[j])), eachindex(vals)))
    for (i, k) in enumerate(t.keys)
        row(string(k), t.scales[i, :])
    end
    col(f) = [(c = filter(!isnan, t.scales[:, j]); isempty(c) ? NaN : f(c)) for j in eachindex(t.phases)]
    row("min", col(minimum))
    row("mean", col(c -> sum(c) / length(c)))
    row("max", col(maximum))
end

"""
    parameters(p::PhasedPointing, beams::PerPhase) -> Vector{Float64}
    parameters(p::PhasedPointing, values::Pair...; group = (field = value, ...), ...)
        -> Vector{Float64}

The parameter vector of `p` that reproduces `beams` (a [`PerPhase`](@ref) of beam
vectors, such as a phased result's `beams` in another parameterization: staging from
quads to beams, or to tied scales): the pointing's parameters from the first phase's
beams (`parameters(pointing, beams)`, matched by name; the aims are the same in every
phase), and each scale the mean over its group's beams of the beam's power over the
phase's power (beams at zero phase power left out; a scale with none is 1). Every
scaled phase must be among `beams`.

With named values, as for [`PowerScales`](@ref): the pointing's fields pass through,
and the field `power_scale` of a key sets its scales, a number for every scaled phase
or a named tuple by phase, `"Q31T" => (z = 0.1, power_scale = (picket = 1.1, rise = 0.95))`.
"""
function parameters(p::PhasedPointing, beams::PerPhase)
    k1 = findfirst(k -> haskey(beams, k), p.phases)
    k1 === nothing && throw(ArgumentError("the beams have none of the phases $(join(repr.(p.phases), ", "))"))
    xp = parameters(p.pointing, beams[p.phases[k1]])
    tmpl = apply(p.pointing, xp)
    acc, cnt = zeros(length(p.keys)), zeros(Int, length(p.keys))
    for (k, ph) in enumerate(p.phases)
        isempty(p.group[k]) && continue
        haskey(beams, ph) || throw(ArgumentError("no beams for the scaled phase $(repr(ph))"))
        byname = Dict(b.name => b for b in beams[ph])
        for (i, t) in enumerate(tmpl)
            haskey(byname, t.name) || throw(ArgumentError("no beam named $(t.name) in phase $(repr(ph))"))
            P0 = p.power[k] === nothing ? Float64(ForwardDiff.value(t.power)) : p.power[k][i]
            P0 > 0 || continue
            g = p.group[k][i]
            acc[g] += Float64(ForwardDiff.value(byname[t.name].power)) / P0
            cnt[g] += 1
        end
    end
    return vcat(xp, [c > 0 ? a / c : 1.0 for (a, c) in zip(acc, cnt)])
end

parameters(p::PhasedPointing, values::Pair...; kwargs...) = _with_named(p, (values..., pairs(kwargs)...))

function _named_values(p::PhasedPointing, kv)
    out = Pair{Int,Any}[]
    rest = Pair[]
    for (key, fv) in kv
        others = Pair{Symbol,Any}[]
        for (f, v) in pairs(fv)
            if f === :power_scale
                hit = false
                for (g, (k, ph)) in enumerate(p.keys)
                    isequal(k, key) || continue
                    if v isa Real
                        push!(out, p.n + g => v); hit = true
                    elseif haskey(v, ph)
                        push!(out, p.n + g => v[ph]); hit = true
                    end
                end
                if !(v isa Real)
                    for ph in keys(v)
                        any(kp -> isequal(kp[1], key) && kp[2] === ph, p.keys) ||
                            throw(ArgumentError("no power scale for $(repr(key)) in phase $(repr(ph))"))
                    end
                end
                hit || throw(ArgumentError("no power scale for $(repr(key)); the scales' keys are " *
                                           join(unique(repr.(first.(p.keys))), ", ")))
            else
                push!(others, f => v)
            end
        end
        isempty(others) || push!(rest, key => (; others...))
    end
    return append!(out, _named_values(p.pointing, rest))
end

"""
    bounds(p::PhasedPointing, x0 = initial_parameters(p); scale = 0.5, scale_min = 0.0,
           scale_max = Inf, power_max = Inf, kwargs...) -> (lower, upper)

The pointing's box (its `bounds` about `x0[1:n]` with `kwargs`) for the first
parameters, and each power scale within `x0 (1 ± scale)` (`scale = Inf`: no relative
limit), clipped to `[scale_min, scale_max]` and to the laser's limit: `power_max` (TW
per beam, a number or a function of a beam) caps every beam in every scaled phase, so a
scale is at most the least `power_max / P` over its group's beams in its phase (`P` the
phase's power; the pointing's own at `x0` in a phase without powers). A fixed phase with
a beam above `power_max` is an `ArgumentError`, and so is an empty box (`scale_min`
above a scale's cap). The defaults are engineering defaults, not facility limits.
"""
function bounds(p::PhasedPointing, x0 = initial_parameters(p); scale = 0.5, scale_min = 0.0,
                scale_max = Inf, power_max = Inf, kwargs...)
    length(x0) == nparams(p) || throw(DimensionMismatch("expected $(nparams(p)) parameters"))
    xp = collect(Float64, x0[1:p.n])
    lo, hi = bounds(p.pointing, xp; kwargs...)
    tmpl = apply(p.pointing, xp)
    lim(b) = power_max isa Real ? Float64(power_max) : Float64(power_max(b))
    cap = fill(Inf, length(p.keys))
    for (k, ph) in enumerate(p.phases)
        P = p.power[k]
        for (i, b) in enumerate(tmpl)
            Pi = P === nothing ? Float64(ForwardDiff.value(b.power)) : P[i]
            L = lim(b)
            if isempty(p.group[k])
                P === nothing || Pi <= L * (1 + 1e-12) ||
                    throw(ArgumentError("phase $(repr(ph)) has fixed powers and beam $(b.name) is at $Pi TW, above power_max = $L"))
            elseif Pi > 0
                g = p.group[k][i]
                cap[g] = min(cap[g], L / Pi)
            end
        end
    end
    s = collect(Float64, x0[p.n+1:end])
    slo = isinf(scale) ? fill(Float64(scale_min), length(s)) : max.(s .* (1 - scale), scale_min)
    shi = isinf(scale) ? min.(scale_max, cap) : min.(s .* (1 + scale), scale_max, cap)
    for g in eachindex(s)
        slo[g] <= shi[g] ||
            throw(ArgumentError("empty power-scale bounds for $(repr(p.keys[g])): $(slo[g]) to $(shi[g])" *
                                (cap[g] < slo[g] ? " (the laser's limit caps the scale at $(cap[g]))" : "")))
    end
    return vcat(lo, slo), vcat(hi, shi)
end

function Base.show(io::IO, p::PhasedPointing)
    print(io, "PhasedPointing(", length(p.phases), " phases (", join(p.phases, ", "), "), ",
          p.n, " pointing parameters, ", length(p.keys), " power scales)")
end

"""
    PhasedScenario(pointing, phases; by = :group, kwargs...)
    PhasedScenario(pointing, phase::Phase...; by = :group, kwargs...)

Several phases of one pulse (a picket, a foot, the peak before and after windows burn
through) under one pointing: a [`Scenario`](@ref) per [`Phase`](@ref), each with its
own target (wall albedos, windows), radiation model, deposition options, observed
elements and laser–plasma model, and its own beam powers, with the ports and spots of
`pointing`'s beams shared and its parameters setting the aims of every phase. The
parameterization (`scenario.param`, a [`PhasedPointing`](@ref)) adds a power scale
per (group, phase) for the phases not fixed (`scaled = false`), groups keyed by `by`
(`:group`, `:quad`, `:beam` or a function of a beam; one key for several quads ties
their scales). `kwargs` are `Scenario` keywords shared by the phases built from a
target (`observe`, `deposition_kw`, `viewfactor_kw`, …; a phase's own override them).

[`evaluate`](@ref)`(scenario, x)` gives a [`PerPhase`](@ref) of [`Evaluation`](@ref)s,
each computed once per parameter vector: [`optimize_pointing`](@ref) passes it to the
objective and every constraint. A builder of an evaluation's objective or constraint
([`mode_gate`](@ref), [`drive_floor`](@ref), [`irradiance_ceiling`](@ref),
[`picket_ceiling`](@ref), [`light_ceiling`](@ref), [`drive_ceiling`](@ref),
[`drive_objective`](@ref), a user [`Constraint`](@ref)) applies to one phase through
[`in_phase`](@ref); a function of the whole `PerPhase` combines phases
(`ev -> w1 * f1(ev[:intact]) + w2 * f2(ev[:burned])`), and a keep-out's
[`clearance_floor`](@ref), which reads the beams alone, takes every phase at once.
[`power_scales`](@ref) reads the power scales inside such a function. The result's
`beams` are a `PerPhase` of beam vectors, and [`pointing_summary`](@ref) reports every
phase.

`scenario[:picket]` is a phase's `Scenario` and `keys(scenario)` the phases' names.
[`remake`](@ref)`(scenario; pointing, by, scaled)` replaces the pointing (staging per
quad, then per beam), the scales' grouping (tied scales) or the scaled phases, keeping
every phase's operators; [`subset`](@ref)`(scenario, keep)` keeps some beams in every
phase. The box of a run: [`bounds`](@ref)`(scenario.param, x0; power_max)` caps every
beam in every phase at the laser's limit. A scenario of one phase whose powers are the
pointing's own and not scaled (`Phase(name, scenario; scaled = false)`) evaluates and
optimizes exactly as that `Scenario` does.

```julia
pointing = GroupPointing(beams; by = :beam, fields = (:z, :phi), relative = true)
ps = PhasedScenario(pointing,
                    Phase(:picket, target_picket; pulses, time = 1.0),
                    Phase(:peak, target_peak; power = 2.0, scaled = false);
                    by = :quad, observe = :capsule)
gates = [in_phase(:picket, mode_gate(ref[:picket], g1)), in_phase(:peak, mode_gate(ref[:peak], g2))]
lower, upper = bounds(ps.param; power_max = 2.0)
r = optimize_pointing(ps, in_phase(:peak, drive_objective(; drive)); constraints = gates, lower, upper)
r.beams[:picket]
```
"""
struct PhasedScenario{P,S<:Scenario,K} <: _AbstractScenario
    phases::Vector{Phase}
    scenarios::Vector{S}
    param::PhasedPointing{P,K}
end

PhasedScenario(pointing, phases::Phase...; kwargs...) = PhasedScenario(pointing, collect(Phase, phases); kwargs...)

function PhasedScenario(pointing, phases::AbstractVector{Phase}; by = :group, kwargs...)
    isempty(phases) && throw(ArgumentError("give at least one phase"))
    names = [ph.name for ph in phases]
    allunique(names) || throw(ArgumentError("phase names must be unique, got $(names)"))
    shared = values(kwargs)
    scen = [_phase_scenario(ph, pointing, shared) for ph in phases]
    return _phased(pointing, collect(Phase, phases), scen, by)
end

# The phased scenario of the built phase scenarios (one concrete `Scenario` type: each is
# on the pointing) and their phases' powers.
function _phased(pointing, phases::Vector{Phase}, scen::AbstractVector, by)
    beams = apply(pointing, initial_parameters(pointing))
    pw = [_phase_power(ph, beams) for ph in phases]
    pp = PhasedPointing(pointing, [ph.name for ph in phases], pw;
                        scaled = [ph.name for ph in phases if ph.scaled], by)
    S = typeof(first(scen))
    all(s -> s isa S, scen) || throw(ArgumentError("the phases' scenarios differ in type"))
    return PhasedScenario(phases, collect(S, scen), pp)
end

# A phase's scenario on the shared pointing: built from its target with the shared and
# its own keywords, or its ready scenario with the pointing as parameterization.
function _phase_scenario(ph::Phase, pointing, shared::NamedTuple)
    ph.source isa Scenario && return remake(ph.source; parameterization = pointing)
    kw = merge(shared, ph.kwargs)
    ph.time === nothing || haskey(kw, :time) || (kw = merge(kw, (; time = ph.time)))
    return Scenario(ph.source, pointing; kw...)
end

Base.keys(ps::PhasedScenario) = ps.param.phases
function Base.getindex(ps::PhasedScenario, k::Symbol)
    i = findfirst(==(k), ps.param.phases)
    i === nothing && throw(KeyError(k))
    return ps.scenarios[i]
end

_nparams_shown(ps::PhasedScenario) = nparams(ps.param)

function Base.show(io::IO, ps::PhasedScenario)
    nb = length(first(ps.scenarios).deposition.ports)
    print(io, "PhasedScenario(", length(ps.phases), " phases (", join(ps.param.phases, ", "), "), ", nb,
          " beams, ", nparams(ps.param), " parameters: ", ps.param.n, " pointing, ", length(ps.param.keys),
          " power scales)")
end

"""
    evaluate(scenario::PhasedScenario, x) / evaluate(scenario::PhasedScenario, beams::PerPhase)
        -> PerPhase

One [`Evaluation`](@ref) per phase: each phase's beams (`apply(scenario.param, x)`, or
`beams`, a `PerPhase` of beam vectors with the scenario's phases in order) in the
phase's scenario. Generic in the number type; the evaluations share nothing mutable, so
concurrent gradient chunks each evaluate their own.
"""
evaluate(ps::PhasedScenario, x::AbstractVector{<:Real}) = evaluate(ps, apply(ps.param, x))

function evaluate(ps::PhasedScenario, beams::PerPhase)
    beams.names == ps.param.phases ||
        throw(ArgumentError("the beams' phases ($(join(beams.names, ", "))) are not the scenario's ($(join(ps.param.phases, ", ")))"))
    return PerPhase(ps.param.phases, _phase_evaluations(ps.scenarios, beams.values))
end

_phase_evaluations(scen::Vector{<:Scenario}, beams::AbstractVector) =
    [evaluate(scen[k], beams[k]) for k in eachindex(scen, beams)]

phase_beams(ps::PhasedScenario, beams::AbstractVector{<:Beam}; kwargs...) = phase_beams(ps.param, beams; kwargs...)
power_scales(ps::PhasedScenario, v) = power_scales(ps.param, v)
phase_powers(ps::PhasedScenario, phase::Symbol) = phase_powers(ps.param, phase)
power_scale_table(ps::PhasedScenario, v) = power_scale_table(ps.param, v)

"""
    remake(scenario::PhasedScenario; pointing, by, scaled) -> PhasedScenario

Copy of a phased scenario with another `pointing` (a parameterization of the same beams,
whose ports the phases' deposition operators serve: per quad, then per beam), another
grouping `by` of the power scales (a function giving several quads one key ties their
scales) or other `scaled` phases (their names); every phase's operators, radiation
model and response are kept, and the phases' powers are set again on the new
pointing's beams.
"""
function remake(ps::PhasedScenario; pointing = _UNCHANGED, by = _UNCHANGED, scaled = _UNCHANGED)
    p = pointing === _UNCHANGED ? ps.param.pointing : pointing
    b = by === _UNCHANGED ? ps.param.by : by
    phases = ps.phases
    if scaled !== _UNCHANGED
        for k in scaled
            k in ps.param.phases || throw(ArgumentError("no phase $(repr(k))"))
        end
        phases = [Phase(ph.name, ph.source, ph.power, ph.pulses, ph.time, ph.pulse_by, ph.name in scaled, ph.kwargs)
                  for ph in phases]
    end
    scen = pointing === _UNCHANGED ? ps.scenarios : [remake(s; parameterization = p) for s in ps.scenarios]
    return _phased(p, phases, scen, b)
end

"""
    subset(scenario::PhasedScenario, keep; pointing = the pointing's subset that keeps its groups)

The phased scenario of the beams `keep` (indices into the pointing's beams) in every
phase, each phase's scenario restricted as by [`subset`](@ref) (its view factors and
response shared), the phases' powers restricted to the kept beams and the scales'
grouping kept. The default pointing keeps the groups (and so the parameter names) of
the scenario's.
"""
function subset(ps::PhasedScenario, keep::AbstractVector{<:Integer}; pointing = nothing)
    p = pointing === nothing ? _subset_param(ps.param.pointing, keep) : pointing
    beams = apply(ps.param.pointing, initial_parameters(ps.param.pointing))
    phases = map(ps.phases) do ph
        P = _phase_power(ph, beams)
        Phase(ph.name, ph.source, P === nothing ? nothing : P[keep], nothing, ph.time, ph.pulse_by, ph.scaled, ph.kwargs)
    end
    scen = [subset(s, keep; parameterization = p) for s in ps.scenarios]
    return _phased(p, phases, scen, ps.param.by)
end

"""
    in_phase(phase, constraint::Constraint) -> Constraint
    in_phase(phase, f) -> function

A constraint or an objective of one [`Evaluation`](@ref) applied to the phase `phase`
(its name) of a [`PhasedScenario`](@ref)'s evaluations: the constraint's function and
labels read `ev[phase]` of the [`PerPhase`](@ref) evaluation, and its name becomes
`"<phase>.<name>"` (e.g. `picket.mode_gate[l = 2]` in the reports); an exact ceiling
([`irradiance_ceiling`](@ref), [`picket_ceiling`](@ref)) keeps its exact value. Any
other function `f` of an evaluation (an objective) becomes `ev -> f(ev[phase])`.

```julia
gate = in_phase(:picket, mode_gate(reference[:picket], gates))
objective = in_phase(:peak, drive_objective(; drive))
```
"""
in_phase(phase::Symbol, c::Constraint) =
    Constraint(Symbol(phase, ".", c.name), _InPhase(phase, c.f), c.limit, c.scale, c.upper,
               c.labels isa Vector{String} ? c.labels : _InPhase(phase, c.labels))
in_phase(phase::Symbol, f) = _InPhase(phase, f)

# A function of one phase's evaluation, applied to a `PerPhase` of them.
struct _InPhase{F} <: Function
    phase::Symbol
    f::F
end
(p::_InPhase)(ev::PerPhase) = p.f(ev[p.phase])
_exact_value(p::_InPhase, ev::PerPhase) = _exact_value(p.f, ev[p.phase])

# A keep-out's clearance over every phase (`clearance_floor` on a `PerPhase` of
# evaluations): per beam the least over the phases. The phases share the aims, so a
# keep-out that does not read the beams' power gives every phase the first's values.
function (c::_ClearanceValue)(ev::PerPhase)
    v = c(ev.values[1])
    for k in 2:length(ev.values)
        v = min.(v, c(ev.values[k]))
    end
    return v
end

_beam_names(ev::PerPhase) = _beam_names(ev.values[1])

check_constraints(ev::PerPhase, c::Constraint; kwargs...) = check_constraints(ev, (c,); kwargs...)
check_constraints(ev::PerPhase, cs; feasibility::Real = 1e-3) =
    _constraint_rows(ev, Tuple(cs), nothing, Float64(feasibility))

"""
    pointing_summary(scenario::PhasedScenario, beams::PerPhase; reference = the start's beams,
                     keepout = nothing, inner = (:inner23, :inner30), drive_reference = nothing,
                     lmax = 4) -> PerPhase
    pointing_summary(scenario::PhasedScenario, result::PointingResult; kwargs...)

A [`PointingSummary`](@ref) per phase: each phase's beams in its scenario against the
`reference` beams (a `PerPhase` of beam vectors, by default the scenario's start
`apply(scenario.param, initial_parameters(scenario.param))`, or one beam vector for
every phase), with the drive against `drive_reference[phase]` where given (a
`PerPhase`, a named tuple or a dictionary by phase of evaluations or coefficient
vectors; a phase without one has no drive ratio). The display is one line per phase.

`gates` and `floor` show each phase's departure against its gates and its pointing
floor, as for one [`Scenario`](@ref): per phase (a `PerPhase`, a named tuple or a
dictionary by phase; a phase without an entry shows none), or one value for every
phase (a vector, a [`mode_uncertainty`](@ref) result, a [`mode_gate`](@ref)). `gates`
may also be the run's constraints (a vector or tuple, e.g. the `constraints` passed to
[`optimize_pointing`](@ref)): each phase then shows its `in_phase(phase, mode_gate(...))`,
from that gate's own reference shape.

```julia
pointing_summary(ps, r; drive_reference = refs, gates = constraints, floor = floors)
```
"""
function pointing_summary(ps::PhasedScenario, beams::PerPhase;
                          reference = apply(ps.param, initial_parameters(ps.param)),
                          keepout = nothing, inner = (:inner23, :inner30),
                          drive_reference = nothing, lmax::Integer = 4, gates = nothing,
                          floor = nothing)
    beams.names == ps.param.phases ||
        throw(ArgumentError("the beams' phases are not the scenario's ($(join(ps.param.phases, ", ")))"))
    out = map(eachindex(ps.scenarios)) do k
        ph = ps.param.phases[k]
        ref = reference isa PerPhase ? reference[ph] : reference
        dr = drive_reference === nothing ? nothing : haskey(drive_reference, ph) ? drive_reference[ph] : nothing
        pointing_summary(ps.scenarios[k], beams[ph]; reference = ref, keepout, inner, drive_reference = dr, lmax,
                         gates = _phase_gate(gates, ph), floor = _phase_entry(floor, ph))
    end
    return PerPhase(ps.param.phases, out)
end

# One phase's entry of a summary's per-phase option: a value for every phase as it is,
# else the phase's entry (`nothing` without one).
_phase_entry(::Nothing, ph) = nothing
_phase_entry(x::Union{AbstractVector{<:Real},Constraint}, ph) = x
_phase_entry(x, ph) = x isa NamedTuple && haskey(x, :rms_mode) ? x : haskey(x, ph) ? x[ph] : nothing

# The gates of a phase: from a run's constraints, the phase's `in_phase` mode gate.
_phase_gate(x, ph) = _phase_entry(x, ph)
function _phase_gate(cs::Union{Tuple{Vararg{Constraint}},AbstractVector{<:Constraint}}, ph)
    i = findfirst(c -> c.f isa _InPhase && c.f.phase === ph && c.f.f isa _GateValue, cs)
    return i === nothing ? nothing : cs[i]
end

_gate_value(f::_InPhase, ev::Evaluation) = _gate_value(f.f, ev)

pointing_summary(ps::PhasedScenario, r::PointingResult; kwargs...) = pointing_summary(ps, r.beams; kwargs...)
