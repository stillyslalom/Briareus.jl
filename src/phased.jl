# This file models several phases of one pulse under one pointing. Each phase has a
# scenario, with its target, radiation model, deposition settings, observation and
# laser–plasma model, and its beam powers. One pointing parameterization sets the aims of
# every phase, and a power scale per (group, phase) is added to its parameters. The
# result is one parameter vector, and one evaluation per phase, for `optimize_pointing`.

"""
    PerPhase{V}
    PerPhase(:picket => a, :peak => b, ...)

One value per phase of a [`PhasedScenario`](@ref), in the order of the phases. The
values are the beams (`apply` of a [`PhasedPointing`](@ref), `PointingResult.beams`), the
evaluations ([`evaluate`](@ref)) or the summaries ([`pointing_summary`](@ref)) of a
phased design. A `PerPhase` is indexed by phase name (`v[:picket]`) or position (`v[1]`).
`keys` gives the names, `values` the values and `pairs` the `name => value` pairs.
Iteration, `length`, `haskey` and `map` (which applies a function to each value and keeps
the names) work as for an ordered dictionary whose iteration yields the values.
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

# Shows a value in the display of a `PerPhase`: arrays and evaluations as a summary, and
# other values in compact form.
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

One phase of a [`PhasedScenario`](@ref). A phase consists of its name, the target that
its beams deposit in and radiate from, and the beam powers in the phase. With a `target`,
the phased scenario builds the [`Scenario`](@ref) of the phase for the shared pointing.
The build uses the `Scenario` keywords of the phased scenario, overridden by the
`kwargs` of this phase (`observe`, `deposition_kw`, `viewfactor_kw`, `lpi`, `albedo`,
`model`, `radiation`, `time`). A phase with separate wall albedos is a target with
separate materials. It can share the view factors of the model of another phase through
`viewfactor_kw = (; F = model.F)`. A ready `scenario` is used with its operators,
radiation model and options, and its parameterization is replaced by the shared
pointing. Its deposition operator must then serve the beams of the pointing, as for
`remake(scenario; parameterization)`.

One target, view-factor matrix and deposition operator can be used for every phase,
with windows intact or burned through, as `remake`s of a shared scenario with the
albedos and time of each phase. `history`, a [`drive_history`](@ref) of the target of
the phase run with `keep_albedo = true`, gives the albedos at the `time` of the phase.
`Phase(:rise, shared; history = h, time = τ, pulses)` is `Phase(:rise, remake(shared;
history = h, time = τ); pulses)`, with the albedos of
[`history_state`](@ref)`(target, h, τ)`. With a `target`, they are the `albedo` of the
built scenario. A reference scenario at the same time is `remake(reference_scenario;
history = reference_history, time = τ)`, so that the reference checks of
[`mode_gate`](@ref) and [`drive_floor`](@ref) use the time of the phase.

The beam powers of the phase (TW per beam) replace those of the pointing:

- `power`: a number for every beam, a vector (one per beam of the pointing, in its
  order) or a function of a beam. With `nothing`, the powers given by the pointing are
  kept, which are those of its template or of a `:power` field.
- `pulses` at `time` (ns): the power of each beam is the power of its pulse at that
  time, as set by [`with_pulses`](@ref)`(beams, pulses, time; by = pulse_by)`.

`time` is also the time of the built scenario, which sets the albedos of the materials
and the laser–plasma model, unless `kwargs` give another. The time of a ready scenario
is kept. `scaled = false` fixes the powers of the phase, for example for a peak at the
laser limit. The phase then has no power scales among the parameters.
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

# Returns the powers of the phase per beam of `beams`, the beams of the pointing at its
# initial parameters. Returns `nothing` for the powers of the pointing.
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

The parameterization of a [`PhasedScenario`](@ref), its `param`. The parameters of
`pointing` set the aims of every phase. `pointing` is any parameterization with `apply`
and `initial_parameters`, such as a [`GroupPointing`](@ref) by cone, quad or beam
(`relative` or not), a [`FramePointing`](@ref) or a [`HeldParameters`](@ref) of either.
A power scale per (group, phase) follows for the phases in `scaled`. `phases` are the
names of the phases in order, and `power[k]` the powers of phase `k`, a vector per beam
of the pointing or `nothing` for the powers of the pointing. Groups are keyed by `by` on
the beams of the pointing as for [`PowerScales`](@ref) (`:group`, `:quad`, `:beam` or a
function of a beam). A function that gives several quads one key ties their scales, with
one pulse shape for them and a multiplier per phase. The keys of the scales, `keys`, are
`(key, phase)` tuples, phase by phase in order of first appearance. The scales are named
`"<key>.power_scale.<phase>"` and follow the parameter names of the pointing.

`apply(p, x)` gives a [`PerPhase`](@ref) of beam vectors. They are the beams of the
pointing at `x[1:n]` with the power of each phase, multiplied by the scale of the group
in a scaled phase. Generic in the number type.
"""
struct PhasedPointing{P,K}
    pointing::P
    n::Int                                          # parameters of `pointing`
    phases::Vector{Symbol}
    power::Vector{Union{Nothing,Vector{Float64}}}   # per-beam powers of each phase, or `nothing`
    keys::Vector{Tuple{K,Symbol}}                   # (group key, phase) of each scale
    group::Vector{Vector{Int}}                      # scale of each beam per phase (empty if fixed)
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
Names of the parameters of a `PhasedPointing`. The names of the pointing come first
(`"x[i]"` for a pointing without names), then `"<key>.power_scale.<phase>"` per scale.
"""
parameter_names(p::PhasedPointing) =
    vcat(_inner_names(p.pointing, p.n), [string(k, ".power_scale.", ph) for (k, ph) in p.keys])

"""
    initial_parameters(p::PhasedPointing) -> Vector{Float64}

The initial parameters of the pointing, then a scale of 1 per (group, phase).
"""
initial_parameters(p::PhasedPointing) = vcat(initial_parameters(p.pointing), ones(length(p.keys)))

_nparams_shown(p::PhasedPointing) = nparams(p)

"""
    apply(p::PhasedPointing, x) -> PerPhase

The beams of every phase at `x`. They are the beams of the pointing at `x[1:n]`, one set
of aims, with the powers of each phase. In a scaled phase each power is multiplied by
the scale of its group in `x[n+1:end]`. The element type follows `x`.
"""
function apply(p::PhasedPointing, x::AbstractVector{T}) where {T}
    length(x) == nparams(p) || throw(DimensionMismatch("expected $(nparams(p)) parameters"))
    beams = apply(p.pointing, x[1:p.n])
    return PerPhase(p.phases, [_phase_beams(beams, p.power[k], p.group[k], x, p.n) for k in eachindex(p.phases)])
end

# Returns the beams of the pointing with phase powers `P`, or with their present powers
# when `P` is `nothing`. Each power is multiplied by its scale `x[n + group[i]]`; an empty
# `group` marks a fixed phase. This is a function barrier, inside which the beam type is
# known.
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

The beams `beams` in every phase of `p`, a [`PhasedScenario`](@ref) or its
[`PhasedPointing`](@ref). `beams` are beams of the pointing in its order, which set the
aims, spots and focus. Their powers are those of each phase, multiplied by `scale` in
the scaled phases; `scale` is a number or a function `(beam, phase) -> scale`. A
reference pointing is evaluated over the same phases with `evaluate(scenario,
phase_beams(scenario, reference_beams))`. In a phase whose powers are those of the
pointing, the beam powers are kept and multiplied by the scale.
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
order, where `p` is a [`PhasedScenario`](@ref) or its [`PhasedPointing`](@ref). They are
the parameters `x[n+1:end]` by which `apply` multiplies the powers of the phases. From a
[`PerPhase`](@ref) of evaluations (their setup beams, `ev[phase].beams`) or of beam
vectors (in the order of the pointing, as `apply` and `evaluate` give them), each scale
is the power of a beam over its phase power ([`phase_powers`](@ref)). The beam used is
the first beam of the key with positive phase power. A key whose beams are all at zero
power in the phase has no scale to recover and gives 1, as in [`parameters`](@ref). A
scaled phase whose powers are the pointing powers (no `power` or `pulses` in its
[`Phase`](@ref)) has no phase power to divide by, and an `ArgumentError` is thrown; the
scales of such a phase come from `x`.

Generic in the number type and thread-safe. Inside an objective of the `PerPhase`
evaluation the scales are duals of the optimizer, for example in a regularization of the
scales toward `s0` (a vector in the order of the keys),
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

# A function barrier, inside which the number type of the beams is known.
function _read_scales(p::PhasedPointing, bs::AbstractVector{<:AbstractVector{Beam{T}}}) where {T}
    s = ones(T, length(p.keys))
    done = falses(length(p.keys))
    for k in eachindex(p.phases)
        g = p.group[k]
        isempty(g) && continue
        P = p.power[k]
        P === nothing &&
            throw(ArgumentError("phase $(repr(p.phases[k])) is scaled on the powers of the pointing: " *
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

The powers of the phase named `phase` of `p` before its power scales, in TW per beam of
the pointing and in its order. `p` is a [`PhasedScenario`](@ref) or its
[`PhasedPointing`](@ref). The powers are the `power` of the [`Phase`](@ref), or its
`pulses` at its time, and `nothing` for a phase whose powers are those of the pointing.
The result is a copy.
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
[`power_scale_table`](@ref). The table has a row per key (`keys`, the group keys of the
scales in order of first appearance) and a column per scaled phase (`phases`).
`scales[i, j]` is the scale of `keys[i]` in `phases[j]`, and `NaN` where that key has no
scale in the phase. The text display prints the table and the least, mean and largest
scale of each phase for the log of a study; other reports use the fields.
"""
struct PowerScaleTable{K}
    keys::Vector{K}
    phases::Vector{Symbol}
    scales::Matrix{Float64}
end

"""
    power_scale_table(p, result::PointingResult) -> PowerScaleTable
    power_scale_table(p, x::AbstractVector) / power_scale_table(p, beams_or_ev::PerPhase)

The power scales of a phased design as a [`PowerScaleTable`](@ref) in Float64, with the
values of any duals. `p` is a [`PhasedScenario`](@ref) or its [`PhasedPointing`](@ref).
The design is given as a result, a parameter vector, or the evaluations or beams that
[`power_scales`](@ref) accepts.
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

The parameter vector of `p` that reproduces `beams`, a [`PerPhase`](@ref) of beam
vectors. An example is the `beams` of a phased result in another parameterization, when
staging from quads to beams or to tied scales. The parameters of the pointing come from
the beams of the first phase (`parameters(pointing, beams)`, matched by name), since the
aims are the same in every phase. Each scale is the mean, over the beams of its group,
of the beam power over the phase power. Beams at zero phase power are left out, and a
scale with no other beams is 1. Every scaled phase must be among `beams`.

With named values, as for [`PowerScales`](@ref), the fields of the pointing are passed
to the pointing, and the field `power_scale` of a key sets its scales. Its value is a
number for every scaled phase or a named tuple by phase,
`"Q31T" => (z = 0.1, power_scale = (picket = 1.1, rise = 0.95))`.
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
           scale_max = Inf, power_max = Inf, pulses = nothing, times = nothing,
           peak = nothing, kwargs...) -> (lower, upper)

The box of the pointing (its `bounds` about `x0[1:n]` with `kwargs`) for the first
parameters, and a range for each power scale. Each scale lies within `x0 (1 ± scale)`,
where `scale = Inf` removes the relative limit, and is clipped to `[scale_min, scale_max]`
and to the laser limit. `power_max` (TW per beam, a number or a function of a beam) caps
every beam in every scaled phase, so a scale is at most the least `power_max / P` over
the beams of its group in its phase. `P` is the phase power, or the pointing power at
`x0` in a phase without powers. A fixed phase with a beam above `power_max` is an
`ArgumentError`, and so is an empty box (`scale_min` above the cap of a scale). The
defaults are engineering defaults; no facility limit sets them.

**Between the phases.** The scale bounds limit the beams only at the times of the
phases. A design is delivered as the pulse [`scaled_pulse`](@ref) (or
[`shape_pulses`](@ref)) of the base pulse of each beam, with the scales of its key as
multipliers at the times of the scaled phases. The multipliers are interpolated between
those times, so a base pulse that rises between two phase samples can exceed
`power_max` where neither sample does. Given `pulses`, `times` and `peak`, the box also
includes every condition of [`pulse_ceiling`](@ref) that bounds only one scale. `pulses`
is the base pulse of each beam, a [`PulseShape`](@ref) for every beam or a function of a
beam. `times` is the time of each phase in ns by name, as a named tuple or dictionary,
or a vector in the order of the phases. `peak` is the value given to `scaled_pulse`.
The conditions on one scale are:

- the pulse before the first scaled phase, with its multiplier held at the first scale;
- the pulse after the last scaled phase, on the `peak` ramp to 1, or held at the last
  scale without `peak`;
- on each interval between two scaled phases, the cap that follows with the neighboring
  scale at its least value in the box.

No design whose pulses stay within `power_max` is excluded by these caps. The conditions
that couple two neighboring scales, for the multiplier between their times, are linear
in both, and no box can express them; [`pulse_ceiling`](@ref) imposes them as a
constraint. The power of each scaled phase must be the base pulse at the phase time.
Otherwise an `ArgumentError` is thrown, since the multipliers would scale another pulse.
Without `pulses` the caps are those at the phase times, and
`scaled_pulse(...; power_max)` checks a pulse table.
"""
function bounds(p::PhasedPointing, x0 = initial_parameters(p); scale = 0.5, scale_min = 0.0,
                scale_max = Inf, power_max = Inf, pulses = nothing, times = nothing, peak = nothing,
                kwargs...)
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
    rows = pulses === nothing ? nothing : _pulse_rows(p, tmpl, pulses, times, peak, lim)
    if rows !== nothing
        for r in rows                      # the conditions on one scale
            r.j == 0 && (cap[r.i] = min(cap[r.i], (r.limit - r.c) / r.a))
        end
    end
    s = collect(Float64, x0[p.n+1:end])
    slo = isinf(scale) ? fill(Float64(scale_min), length(s)) : max.(s .* (1 - scale), scale_min)
    if rows !== nothing
        for r in rows                      # conditions on two scales, the other at its least value
            r.j == 0 && continue
            cap[r.i] = min(cap[r.i], (r.limit - r.b * slo[r.j]) / r.a)
            cap[r.j] = min(cap[r.j], (r.limit - r.a * slo[r.i]) / r.b)
        end
    end
    shi = isinf(scale) ? min.(scale_max, cap) : min.(s .* (1 + scale), scale_max, cap)
    for g in eachindex(s)
        slo[g] <= shi[g] ||
            throw(ArgumentError("empty power-scale bounds for $(repr(p.keys[g])): $(slo[g]) to $(shi[g])" *
                                (cap[g] < slo[g] ? " (the laser's limit caps the scale at $(cap[g]))" : "")))
    end
    return vcat(lo, slo), vcat(hi, shi)
end

"""
    pulse_ceiling(p; pulses, times, power_max, peak = nothing, name = :pulse) -> Constraint

A ceiling on the pulses of a phased design between the phase samples as well as at them. The
[`Constraint`](@ref) requires every point of the [`scaled_pulse`](@ref) of each beam to stay
at or below `power_max`. That pulse is the base pulse of the beam, with the power scales of
its key as multipliers at the times of the scaled phases, ramped to 1 at `peak`. `p` is a
[`PhasedScenario`](@ref) or its [`PhasedPointing`](@ref). `pulses` is the base pulse of each
beam, a [`PulseShape`](@ref) for every beam or a function of a beam. `times` is the time of
each phase in ns by name, as a named tuple or dictionary, or a vector in the order of the
phases. `power_max` is the limit in TW per beam (a number, or a function of a beam), and
`peak` the time given to `scaled_pulse`. The power of each scaled phase must be the base
pulse at the phase time, or an `ArgumentError` is thrown. Every scaled phase needs powers, a
`power` or `pulses` in the [`Phase`](@ref), since the constraint takes the scales from the
beams of the evaluation ([`power_scales`](@ref)). The constraint is generic in the number
type.

The written pulse is linear between its points, which are the points of the base pulse
and the phase times, so its largest power is at a point. At a point τ between the scaled
phases at `tᵢ` and `tᵢ₊₁`, the multiplier is `(1 − w) sᵢ + w sᵢ₊₁` with
`w = (τ − tᵢ)/(tᵢ₊₁ − tᵢ)`. The condition `P(τ) ((1 − w) sᵢ + w sᵢ₊₁) ≤ power_max` is
therefore linear in the two scales. Before the first phase and after the last, the
condition bounds one scale; on the `peak` ramp the multiplier runs to 1. Of the points
between two phases, only those on the lower convex hull of `(w, power_max/P(τ))` are
kept, since their conditions imply those at the other points. The constraint is
therefore exact, with a few components per key and interval: a design meets it if and
only if its written pulses stay within `power_max`. The pulse from `peak` on is the base
pulse, which lies outside the constraint; `scaled_pulse(...; power_max)` checks it.

The value of a component is the written power at its point (TW), and its limit is
`power_max`. The label gives the key and the time (`pulse[outer, 4.35 ns]`), and the
default `scale` makes the violation relative to the limit.
[`bounds`](@ref)`(p.param; power_max, pulses, times, peak)` imposes exactly, as a box,
the components that bound only one scale. With both, the augmented Lagrangian includes
only the coupled components, each to within the `feasibility` of the run (1e-3 of the
limit by default), so a written pulse can end that far over the limit. With a limit
lowered by the `feasibility`, `power_max = (1 - 1e-3) * L`, the pulse stays within `L`.
"""
function pulse_ceiling(p::PhasedPointing; pulses, times, power_max, peak = nothing, name = :pulse)
    for k in eachindex(p.phases)
        isempty(p.group[k]) || p.power[k] !== nothing ||
            throw(ArgumentError("pulse_ceiling: the scaled phase $(repr(p.phases[k])) has the pointing's powers; " *
                                "give it a `power` or `pulses` so that its scales can be read"))
    end
    tmpl = apply(p.pointing, initial_parameters(p.pointing))
    lim(b) = power_max isa Real ? Float64(power_max) : Float64(power_max(b))
    rows = _pulse_rows(p, tmpl, pulses, times, peak, lim)
    isempty(rows) && throw(ArgumentError("pulse_ceiling: no point of the pulses depends on a power scale"))
    # `power_scales` returns 1 for a scale whose key has no power in its phase, whatever
    # the parameter
    readable = falses(length(p.keys))
    for k in eachindex(p.phases), (i, g) in enumerate(p.group[k])
        p.power[k][i] > 0 && (readable[g] = true)
    end
    for r in rows
        readable[r.i] && (r.j == 0 || readable[r.j]) ||
            throw(ArgumentError("pulse_ceiling: $(repr(p.keys[readable[r.i] ? r.j : r.i])) has no power in its phase, " *
                                "so its scale cannot be read from the evaluation, but its pulse depends on it at $(r.label)"))
    end
    value = _PulseCeilingValue(p, [r.i for r in rows], [r.j for r in rows], [r.a for r in rows],
                               [r.b for r in rows], [r.c for r in rows])
    return Constraint(name, value; upper = [r.limit for r in rows], labels = [r.label for r in rows])
end

# The written power at the point of each component, `a sᵢ + b sⱼ + c` (`j = 0` for one
# scale), from the scales of the evaluation or of a `PerPhase` of beams. Generic in the
# number type.
struct _PulseCeilingValue{P}
    param::P
    i::Vector{Int}
    j::Vector{Int}
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
end

(v::_PulseCeilingValue)(ev) = _pulse_powers(v, power_scales(v.param, ev))

function _pulse_powers(v::_PulseCeilingValue, s::AbstractVector{T}) where {T}
    out = Vector{T}(undef, length(v.i))
    for r in eachindex(out)
        x = v.a[r] * s[v.i[r]] + v.c[r]
        v.j[r] == 0 || (x += v.b[r] * s[v.j[r]])
        out[r] = x
    end
    return out
end

# One condition on a written pulse at a point: `a sᵢ + b sⱼ + c ≤ limit` (TW), with sᵢ and
# sⱼ the scales on either side of the point (`j = 0` for one scale), and `c` the share of
# the `peak` ramp.
struct _PulseRow
    i::Int
    j::Int
    a::Float64
    b::Float64
    c::Float64
    limit::Float64
    label::String
    time::Float64
end

# Returns the conditions that keep the written pulse of every beam at or below its limit.
# The written pulse is the `scaled_pulse` of the base pulse of the beam, with the scales
# of its key at the times of the scaled phases and `peak`. A check first confirms that the
# power of each scaled phase is the base pulse at its time. The conditions are, per key,
# the least bound on each single scale, and between two phases the points on the lower
# convex hull of (w, limit/P), which imply the others. `tmpl` contains the beams of the
# pointing and, for a phase without powers, their powers.
function _pulse_rows(p::PhasedPointing, tmpl::AbstractVector{<:Beam}, pulses, times, peak, lim)
    times === nothing && throw(ArgumentError("`pulses` need the phases' `times`"))
    sc = [k for k in eachindex(p.phases) if !isempty(p.group[k])]
    isempty(sc) && return _PulseRow[]
    ts = Float64[times isa AbstractVector ? times[k] : times[p.phases[k]] for k in sc]
    n = length(ts)
    issorted(ts; lt = <=) || throw(ArgumentError("the scaled phases' times must increase in their order"))
    tp = peak === nothing ? Inf : Float64(peak)
    tp > ts[end] || throw(ArgumentError("peak ($peak ns) must follow the last scaled phase ($(ts[end]) ns)"))
    knots = isfinite(tp) ? [ts; tp] : ts
    single = Dict{Int,Tuple{Float64,_PulseRow}}()           # scale => (its bound, the row)
    inner = Dict{Tuple{Int,Int},Vector{Tuple{Float64,Float64,_PulseRow}}}()  # (key, interval) => (w, limit/P, row)
    function one!(g, a, c, L, label, τ)
        bound = (L - c) / a
        (haskey(single, g) && single[g][1] <= bound) || (single[g] = (bound, _PulseRow(g, 0, a, 0.0, c, L, label, τ)))
    end
    for (i, b) in enumerate(tmpl)
        q = pulses isa PulseShape ? pulses : pulses(b)
        g = [p.group[k][i] for k in sc]                     # scale of the beam in each scaled phase
        for (j, k) in enumerate(sc)
            P = p.power[k] === nothing ? Float64(ForwardDiff.value(b.power)) : p.power[k][i]
            isapprox(q(ts[j]), P; rtol = 1e-9) ||
                throw(ArgumentError("phase $(repr(p.phases[k])): beam $(b.name)'s power, $P TW, is not its base " *
                                    "pulse's at $(ts[j]) ns ($(q(ts[j])) TW)"))
        end
        L = lim(b)
        key = p.keys[g[1]][1]
        lo, hi = q.t[1], q.t[end]
        for τ in sort!(unique!([q.t; filter(τ -> lo < τ < hi, knots)]))
            P = q(τ)
            (P > 0 && τ < tp) || continue
            label = string(key, ", ", round(τ; sigdigits = 6), " ns")
            τ <= ts[1] && one!(g[1], P, 0.0, L, label, τ)
            if τ >= ts[n]
                w = isfinite(tp) ? (τ - ts[n]) / (tp - ts[n]) : 0.0
                one!(g[n], P * (1 - w), P * w, L, label, τ)
            end
            for s in 1:n-1
                ts[s] <= τ <= ts[s+1] || continue
                w = (τ - ts[s]) / (ts[s+1] - ts[s])
                push!(get!(() -> Tuple{Float64,Float64,_PulseRow}[], inner, (g[1], s)),
                      (w, L / P, _PulseRow(g[s], g[s+1], P * (1 - w), P * w, 0.0, L, label, τ)))
            end
        end
    end
    rows = _PulseRow[]
    for ((_, s), pts) in inner
        for (w, _, r) in _lower_hull(pts)
            if w == 0
                one!(r.i, r.a, 0.0, r.limit, r.label, r.time)
            elseif w == 1
                one!(r.j, r.b, 0.0, r.limit, r.label, r.time)
            else
                push!(rows, r)
            end
        end
    end
    append!(rows, last.(values(single)))
    order = Dict{Any,Int}()                                 # the place of each key, in order of first appearance
    for (g, kk) in enumerate(p.keys)
        get!(order, kk[1], g)
    end
    return sort!(rows; by = r -> (order[p.keys[r.i][1]], r.time, r.j))
end

# Returns the lower convex hull of points (w, c, row). Of points tied in w, the one with
# the least c is kept. The conditions `(1 − w) sᵢ + w sⱼ ≤ c` at the hull vertices imply
# those at every other point, because an affine function of w that lies below the
# vertices lies below the hull.
function _lower_hull(pts::Vector{T}) where {T<:Tuple}
    sorted = sort(pts; by = x -> (x[1], x[2]))
    h = T[]
    for x in sorted
        isempty(h) || h[end][1] != x[1] || continue          # the least c at each w
        while length(h) >= 2 && _cross(h[end-1], h[end], x) <= 0
            pop!(h)
        end
        push!(h, x)
    end
    return h
end

_cross(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])

function Base.show(io::IO, p::PhasedPointing)
    print(io, "PhasedPointing(", length(p.phases), " phases (", join(p.phases, ", "), "), ",
          p.n, " pointing parameters, ", length(p.keys), " power scales)")
end

"""
    PhasedScenario(pointing, phases; by = :group, kwargs...)
    PhasedScenario(pointing, phase::Phase...; by = :group, kwargs...)

Several phases of one pulse under one pointing, such as a picket, a foot, and the peak
before and after the windows burn through. Each [`Phase`](@ref) has a
[`Scenario`](@ref) with its target (wall albedos, windows), radiation model, deposition
options, observed elements and laser–plasma model, and its beam powers. The ports and
spots of the beams of `pointing` are shared, and its parameters set the aims of every
phase. The parameterization (`scenario.param`, a [`PhasedPointing`](@ref)) adds a power
scale per (group, phase) for the phases not fixed with `scaled = false`. Groups are
keyed by `by` (`:group`, `:quad`, `:beam` or a function of a beam), and one key for
several quads ties their scales. `kwargs` are `Scenario` keywords shared by the phases
built from a target (`observe`, `deposition_kw`, `viewfactor_kw`, …), and the keywords
of a phase override them.

[`evaluate`](@ref)`(scenario, x)` gives a [`PerPhase`](@ref) of [`Evaluation`](@ref)s,
computed once per parameter vector, and [`optimize_pointing`](@ref) passes it to the
objective and every constraint. A builder of an objective or constraint of one
evaluation is applied to one phase through [`in_phase`](@ref). Such builders are
[`mode_gate`](@ref), [`drive_floor`](@ref), [`irradiance_ceiling`](@ref),
[`picket_ceiling`](@ref), [`light_ceiling`](@ref), [`drive_ceiling`](@ref),
[`drive_objective`](@ref) and user [`Constraint`](@ref)s. A function of the whole
`PerPhase` combines phases (`ev -> w1 * f1(ev[:intact]) + w2 * f2(ev[:burned])`). The
[`clearance_floor`](@ref) of a keep-out depends only on the beams and takes every phase
at once. [`power_scales`](@ref) extracts the power scales inside such a function. The
`beams` of the result are a `PerPhase` of beam vectors, and [`pointing_summary`](@ref)
reports every phase.

`scenario[:picket]` is the `Scenario` of a phase, and `keys(scenario)` gives the names of
the phases. [`remake`](@ref)`(scenario; pointing, by, scaled)` replaces the pointing
(staging per quad, then per beam), the grouping of the scales (tied scales) or the
scaled phases; the operators of every phase are kept. [`subset`](@ref)`(scenario, keep)`
keeps some beams in every phase. For the box of a run,
[`bounds`](@ref)`(scenario.param, x0; power_max)` caps every beam in every phase at the
laser limit. A scenario of one phase whose powers are the pointing powers and not scaled
(`Phase(name, scenario; scaled = false)`) evaluates and optimizes exactly as that
`Scenario` does.

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
    _mesh::Vector{Int}   # per phase, the first phase with an operator on the same points (`_operator_classes`)
end

PhasedScenario(phases::Vector{Phase}, scen::Vector{S}, param::PhasedPointing{P,K}) where {P,S<:Scenario,K} =
    PhasedScenario{P,S,K}(phases, scen, param, _operator_classes(scen))

Base.propertynames(::PhasedScenario) = (:phases, :scenarios, :param)

# Returns, for each phase, the first phase whose deposition operator has the same points,
# weights, normals, visibility, destinations, cone directions, ports and profile. That is
# the same operator, or one built alike on a target that differs only in its materials.
# Phases of one class whose beams differ only in power are deposited in one pass
# (`_deposit_groups`).
function _operator_classes(scen::AbstractVector{<:Scenario})
    cls = collect(eachindex(scen))
    for k in eachindex(scen)
        for j in 1:k-1
            cls[j] == j && _same_points(scen[j].deposition, scen[k].deposition) && (cls[k] = j; break)
        end
    end
    return cls
end

function _same_points(a::DepositionOperator, b::DepositionOperator)
    a === b && return true
    same(x, y) = x === y || x == y
    nelements(a.target) == nelements(b.target) && length(a.target.components) == length(b.target.components) ||
        return false
    (a.profile === nothing) == (b.profile === nothing) || return false
    a.profile === nothing || (a.profile.taper == b.profile.taper && a.envelope == b.envelope) || return false
    return a.ports == b.ports && same(a.dirs, b.dirs) && same(a.pts, b.pts) && same(a.w, b.w) &&
           same(a.normal, b.normal) && same(a.mask, b.mask) && same(a.dest, b.dest)
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

# Returns the phased scenario of the built phase scenarios and the powers of their
# phases. The phase scenarios share one concrete `Scenario` type, since each is on the
# pointing.
function _phased(pointing, phases::Vector{Phase}, scen::AbstractVector, by)
    beams = apply(pointing, initial_parameters(pointing))
    pw = [_phase_power(ph, beams) for ph in phases]
    pp = PhasedPointing(pointing, [ph.name for ph in phases], pw;
                        scaled = [ph.name for ph in phases if ph.scaled], by)
    S = typeof(first(scen))
    all(s -> s isa S, scen) || throw(ArgumentError("the phases' scenarios differ in type"))
    return PhasedScenario(phases, collect(S, scen), pp)
end

# Returns the scenario of a phase on the shared pointing. It is built from the target of
# the phase with the shared and the per-phase keywords, or it is the ready scenario of the
# phase with the pointing as parameterization.
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

One [`Evaluation`](@ref) per phase, of the beams of each phase in the scenario of that
phase. The beams are `apply(scenario.param, x)`, or `beams`, a `PerPhase` of beam vectors
with the phases of the scenario in order. Generic in the number type. No mutable state
is shared between the evaluations, so concurrent gradient chunks can evaluate them
independently.

Some phases are deposited in one pass: those whose deposition operators have the same
points (one operator, or operators built alike on targets that differ only in their
materials) and whose effective beams differ only in power. The beam profile is then
evaluated once per point and scaled by the beam powers of each phase, and the deposition
of each phase equals a separate `deposit` of that phase within rounding. A phase is
deposited separately when its effective beams differ from those of the other phases in
aim, spot or focus, when the points of its operator differ, or when its laser–plasma
model deposits light itself.
"""
evaluate(ps::PhasedScenario, x::AbstractVector{<:Real}) = evaluate(ps, apply(ps.param, x))

function evaluate(ps::PhasedScenario, beams::PerPhase)
    beams.names == ps.param.phases ||
        throw(ArgumentError("the beams' phases ($(join(beams.names, ", "))) are not the scenario's ($(join(ps.param.phases, ", ")))"))
    return PerPhase(ps.param.phases, _phase_evaluations(ps, beams.values))
end

# Returns the evaluations of the phases. Phases are deposited in one pass
# (`_deposit_sets`) when they share operator points (`ps._mesh`), their effective beams
# differ only in power, and their laser–plasma models do not deposit light themselves.
# The profile is then evaluated once for them, and the deposition of each phase equals a
# separate `deposit` of that phase within rounding. Other phases are deposited
# separately, including phases whose aims, spots or focus differ because a model changes
# them with power.
function _phase_evaluations(ps::PhasedScenario, beams::AbstractVector)
    scen = ps.scenarios
    lpi = [_resolve_lpi(scen[k].lpi, beams[k]) for k in eachindex(scen, beams)]
    eff = [apply_lpi(lpi[k], beams[k], scen[k].time) for k in eachindex(scen, beams)]
    evs = Vector{Any}(undef, length(scen))
    for g in _deposit_groups(ps, eff)
        if length(g) == 1
            k = only(g)
            evs[k] = _evaluate(scen[k], beams[k], eff[k], lpi[k])
        else
            _shared_evaluations!(evs, scen, beams, eff, g)
        end
    end
    return [ev for ev in evs]
end

# Returns the groups of phases that are deposited together. A phase is added to the
# group of the first earlier phase on the same operator points whose effective beams
# differ from those of the phase only in power (an exact test, `_same_but_power`),
# provided that neither laser–plasma model deposits light.
function _deposit_groups(ps::PhasedScenario, eff::AbstractVector)
    n = length(eff)
    lead = collect(1:n)
    for k in 2:n
        _plain_deposit(ps.scenarios[k].lpi) || continue
        for j in 1:k-1
            lead[j] == j && ps._mesh[j] == ps._mesh[k] && _plain_deposit(ps.scenarios[j].lpi) &&
                _same_shape(eff[j], eff[k]) || continue
            lead[k] = j
            break
        end
    end
    return [findall(==(j), lead) for j in 1:n if lead[j] == j]
end

_same_shape(a::Vector{Beam{T}}, b::Vector{Beam{T}}) where {T} =
    length(a) == length(b) && all(i -> _same_but_power(a[i], b[i]), eachindex(a, b))
_same_shape(a, b) = false

# Returns the evaluations of the phases `g`, one group of `_deposit_groups`, from one pass
# over the operator of the first phase. The group sets of each phase (the split of its
# `Scenario`) are kept. This is a function barrier, inside which the number type of the
# beams is known.
function _shared_evaluations!(evs, scen, beams, eff, g)
    T = typeof(first(eff[g[1]]).power)
    sets = Vector{Beam{T}}[eff[k] for k in g]
    splits = [_GroupSplit(getfield(scen[k], :_split), eff[k], nelements(scen[k].target), T) for k in g]
    deps = _deposit_sets(scen[g[1]].deposition, sets, splits)
    for (i, k) in enumerate(g)
        sp = isempty(splits[i].sets) ? nothing : splits[i]
        evs[k] = _evaluation(scen[k], beams[k], eff[k], deps[i], sp)
    end
    return evs
end

phase_beams(ps::PhasedScenario, beams::AbstractVector{<:Beam}; kwargs...) = phase_beams(ps.param, beams; kwargs...)
power_scales(ps::PhasedScenario, v) = power_scales(ps.param, v)
phase_powers(ps::PhasedScenario, phase::Symbol) = phase_powers(ps.param, phase)
power_scale_table(ps::PhasedScenario, v) = power_scale_table(ps.param, v)
pulse_ceiling(ps::PhasedScenario; kwargs...) = pulse_ceiling(ps.param; kwargs...)

"""
    power_scale_objective(scenario::PhasedScenario; reference = 1, weight = b -> 1.0)

An objective term for [`optimize_pointing`](@ref) on a phased scenario. The term is the
beam-weighted mean squared relative change of the power scales from `reference`,

    Σ_j W_j (s_j / s0_j - 1)² / Σ_j W_j,

over the scales `s_j` of [`power_scales`](@ref), one per `(key, phase)` of
`scenario.param.keys`. `W_j` is the summed `weight` of the beams that the scale applies
to, by default their number. With unit weights the term is the mean, over every beam and
scaled phase, of the squared relative change of the beam scale. A scale shared by a cone
therefore counts as many times as the cone has beams, and a scale per quad as many times
as the quad has beams. The term regularizes a re-pointing whose powers should stay near
a reference drive. It is added to an objective, `ev -> J(ev) + w * reg(ev)` with
`reg = power_scale_objective(ps; reference = s0)`.

`reference` is one scale for every key, a vector in the order of the keys (for example
`power_scales(ps, x0)` of a start), or a function `(key, phase) -> scale`. Its values
must be positive and finite. `weight` is a function of a beam of the pointing (of its
group, quad or name), evaluated once at construction on the beams of the pointing. The
weights must be finite and non-negative with a positive sum. The scales are taken from
the evaluations as in [`power_scales`](@ref), so every scaled phase needs explicit powers
(a `power` or `pulses` in its [`Phase`](@ref)). Generic in the number type and
thread-safe, since it does not modify its tables.
"""
function power_scale_objective(ps::PhasedScenario; reference = 1, weight = b -> 1.0)
    p = ps.param
    nk = length(p.keys)
    nk > 0 || throw(ArgumentError("the phased scenario has no power scales (no scaled phase)"))
    s0 = reference isa Real ? fill(Float64(reference), nk) :
         reference isa AbstractVector ? collect(Float64, reference) :
         Float64[reference(k, ph) for (k, ph) in p.keys]
    length(s0) == nk || throw(DimensionMismatch("$(length(s0)) reference scales for $nk keys"))
    all(s -> isfinite(s) && s > 0, s0) || throw(ArgumentError("the reference scales must be positive and finite"))
    beams = apply(p.pointing, initial_parameters(p.pointing))
    wb = Float64[weight(b) for b in beams]
    all(w -> isfinite(w) && w >= 0, wb) || throw(ArgumentError("the beams' weights must be finite and non-negative"))
    W = zeros(nk)
    for g in p.group, i in eachindex(g)
        W[g[i]] += wb[i]
    end
    sum(W) > 0 || throw(ArgumentError("the beams' weights sum to zero"))
    return _ScaleObjective(p, s0, W, sum(W))
end

struct _ScaleObjective{P<:PhasedPointing}
    param::P
    s0::Vector{Float64}
    w::Vector{Float64}   # summed beam weight of each scale
    W::Float64
end

(o::_ScaleObjective)(ev::PerPhase) = sum(o.w .* abs2.(power_scales(o.param, ev) ./ o.s0 .- 1)) / o.W

"""
    remake(scenario::PhasedScenario; pointing, by, scaled, deposition_kw) -> PhasedScenario

Copy of a phased scenario with another `pointing`, another grouping `by` of the power
scales, or other `scaled` phases (by name). The new `pointing` is a parameterization of
the same beams, whose ports the deposition operators of the phases serve, as when
staging per quad, then per beam. A function that gives several quads one key ties their
scales. The operators, radiation model and response of every phase are kept, and the
powers of the phases are set again on the beams of the new pointing.

With `deposition_kw`, the deposition operator of every phase is rebuilt with these
options for the beams of the pointing, as `remake(scenario; deposition_kw)` does for one
[`Scenario`](@ref). Options not given take their defaults, including those in the
`deposition_kw` of a phase. The target, radiation model, response, laser–plasma model
and time of each phase are kept. An example is a finished design checked on a finer
operator, `evaluate(remake(ps; deposition_kw = (; resolution = 0.3)), r.x)`. With
`pointing` as well, the operators are built for the beams of the new pointing, so its
ports need not be the old ones. Phases that shared one operator share the rebuilt one.
A phase on a target with the same mesh as an earlier phase (a target that differs only
in its materials) is given the rebuilt operator of that phase on its target, so that
phases on one mesh are still deposited in one pass.
"""
function remake(ps::PhasedScenario; pointing = _UNCHANGED, by = _UNCHANGED, scaled = _UNCHANGED,
                deposition_kw = _UNCHANGED)
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
    scen = deposition_kw !== _UNCHANGED ? _rebuilt_operators(ps.scenarios, p, deposition_kw) :
           pointing === _UNCHANGED ? ps.scenarios : [remake(s; parameterization = p) for s in ps.scenarios]
    return _phased(p, phases, scen, b)
end

# Returns the scenarios of the phases on the pointing `p`, with deposition operators
# built with `kw` for its beams and everything else unchanged. One operator is built per
# old operator. A phase on a target with the mesh of an earlier phase (`_same_mesh`) is
# given the operator of that phase on its target, since an operator depends only on the
# geometry and the ports.
function _rebuilt_operators(scen::AbstractVector{<:Scenario}, p, kw)
    beams = _scenario_beams(p)
    built = Tuple{DepositionOperator,DepositionOperator}[]   # (old, rebuilt)
    return map(scen) do s
        i = findfirst(o -> first(o) === s.deposition, built)
        op = if i !== nothing
            last(built[i])
        else
            j = findfirst(o -> _same_mesh(last(o).target, s.target), built)
            new = j === nothing ? DepositionOperator(s.target, beams; kw...) : _on_target(last(built[j]), s.target)
            push!(built, (s.deposition, new))
            new
        end
        Scenario(s.target, p, op, s.model, s.observed, s.moment_weights, s.ylm_weights, s.moment_response,
                 s.ylm_response, s.lpi, s.time, getfield(s, :_response), s.albedo)
    end
end

# Whether two targets have the same mesh, that is, the same triangles, the same elements
# and components of the triangle sides, the same opacity and the same elements. Materials
# may differ. These are all the parts of a target that the construction of a deposition
# operator uses.
_same_mesh(a::Target, b::Target) =
    a === b || nelements(a) == nelements(b) && length(a.components) == length(b.components) &&
    a.tri_v0 == b.tri_v0 && a.tri_e1 == b.tri_e1 && a.tri_e2 == b.tri_e2 && a.tri_normal == b.tri_normal &&
    a.tri_area == b.tri_area && a.tri_comp == b.tri_comp && a.tri_front == b.tri_front &&
    a.tri_back == b.tri_back && a.tri_opaque == b.tri_opaque && a.elem_comp == b.elem_comp &&
    a.elem_side == b.elem_side && a.elem_tris == b.elem_tris

# The operator `op` on the target `t` of the same mesh.
_on_target(op::DepositionOperator, t::Target) =
    op.target === t ? op :
    DepositionOperator(t, op.ports, op.dirs, op.pts, op.w, op.normal, op.mask, op.dest, op.plan, op.profile,
                       op.envelope)

# Returns the phased scenario that the optimizer uses on the VisRad profile with its hard
# cut, with the scenario of each phase tapered (`_tapered(::Scenario)`). Returns `nothing`
# when no phase has a tapered form.
function _tapered(ps::PhasedScenario)
    ts = map(_tapered, ps.scenarios)
    all(isnothing, ts) && return nothing
    return PhasedScenario(ps.phases, [something(t, sc) for (t, sc) in zip(ts, ps.scenarios)], ps.param)
end

"""
    subset(scenario::PhasedScenario, keep; pointing = the pointing subset with the same groups)

The phased scenario of the beams `keep` (indices into the beams of the pointing) in
every phase. The scenario of each phase is restricted as by [`subset`](@ref), with its
view factors and response shared. The powers of the phases are restricted to the kept
beams, and the grouping of the scales is kept. The default pointing has the same groups,
and so the same parameter names, as the pointing of the scenario.
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

A constraint or an objective of one [`Evaluation`](@ref), applied to the phase `phase`
(by name) of the evaluations of a [`PhasedScenario`](@ref). The function and labels of
the constraint use `ev[phase]` of the [`PerPhase`](@ref) evaluation, and its name
becomes `"<phase>.<name>"` (e.g. `picket.mode_gate[l = 2]` in the reports). The exact
value of an exact ceiling ([`irradiance_ceiling`](@ref), [`picket_ceiling`](@ref)) is
kept. Any other function `f` of an evaluation (an objective) becomes
`ev -> f(ev[phase])`.

```julia
gate = in_phase(:picket, mode_gate(reference[:picket], gates))
objective = in_phase(:peak, drive_objective(; drive))
```
"""
in_phase(phase::Symbol, c::Constraint) =
    Constraint(Symbol(phase, ".", c.name), _InPhase(phase, c.f), c.limit, c.scale, c.upper,
               c.labels isa Vector{String} ? c.labels : _InPhase(phase, c.labels))
in_phase(phase::Symbol, f) = _InPhase(phase, f)

# A function of the evaluation of one phase, applied to a `PerPhase` of evaluations.
struct _InPhase{F} <: Function
    phase::Symbol
    f::F
end
(p::_InPhase)(ev::PerPhase) = p.f(ev[p.phase])
_exact_value(p::_InPhase, ev::PerPhase) = _exact_value(p.f, ev[p.phase])
_split_sets(p::_InPhase) = _SplitRequest([something(k, p.phase) => g for (k, g) in _split_sets(p.f)])

# `_with_split` per phase. The group sets that the constraints of each phase use are kept
# in the scenario of that phase.
function _with_split(ps::PhasedScenario, req)
    scen = [_with_split(sc, _SplitRequest([nothing => g for (k, g) in req if k === name]))
            for (sc, name) in zip(ps.scenarios, ps.param.phases)]
    all(i -> scen[i] === ps.scenarios[i], eachindex(scen)) && return ps
    return PhasedScenario(ps.phases, scen, ps.param)
end

# Returns the clearance of a keep-out over every phase (`clearance_floor` on a `PerPhase`
# of evaluations), which is the least over the phases for each beam. The phases share the
# aims, so a keep-out that does not depend on the beam power gives every phase the values
# of the first phase. A phase whose beams give the same clearances as those of the first
# phase (`_same_clearance`, an exact test) is skipped, since its clearances are exactly
# equal to those of the first phase.
function (c::_ClearanceValue)(ev::PerPhase)
    b1 = ev.values[1].beams
    v = c(ev.values[1])
    for k in 2:length(ev.values)
        bk = ev.values[k].beams
        length(bk) == length(b1) && all(i -> _same_clearance(c.keepout, b1[i], bk[i]), eachindex(b1)) &&
            continue
        v = min.(v, c(ev.values[k]))
    end
    return v
end

# Whether the keep-out `k` gives beams `a` and `b` exactly the same clearance. By default
# the test holds when the beams are identical in values and partials. An
# `ApertureKeepOut` uses the beam power only through the envelope radius, computed from
# the values, so for it the test holds when the beams differ at most in power and have
# equal envelope radii.
_same_clearance(k, a::Beam, b::Beam) = _same_but_power(a, b) && _identical(a.power, b.power)
_same_clearance(k::ApertureKeepOut, a::Beam, b::Beam) =
    _same_but_power(a, b) &&
    (ForwardDiff.value(a.power) == ForwardDiff.value(b.power) ||
     envelope_radius(a.spot, a.power; intensity = k.intensity, peak_fraction = k.peak_fraction) ==
     envelope_radius(b.spot, b.power; intensity = k.intensity, peak_fraction = k.peak_fraction))

_same_but_power(a::Beam, b::Beam) =
    a.name == b.name && a.port == b.port && _identical(a.aim, b.aim) && _identical(a.defocus, b.defocus) &&
    _identical(a.spot, b.spot) && a.fnumber == b.fnumber && a.group === b.group && a.quad == b.quad &&
    a.wavelength == b.wavelength && _identical(a.dlambda, b.dlambda)

# Whether two numbers of one type have equal values and partials. `==` on duals compares
# only values, and `isequal` distinguishes -0.0 from 0.0.
_identical(x::D, y::D) where {D<:ForwardDiff.Dual} =
    _identical(ForwardDiff.value(x), ForwardDiff.value(y)) && _identical(ForwardDiff.partials(x).values, ForwardDiff.partials(y).values)
_identical(x::Real, y::Real) = !(x isa ForwardDiff.Dual) && !(y isa ForwardDiff.Dual) && isequal(x, y)
_identical(x::Tuple, y::Tuple) = length(x) == length(y) && all(i -> _identical(x[i], y[i]), eachindex(x, y))
_identical(x::SVector{3}, y::SVector{3}) = _identical(x[1], y[1]) && _identical(x[2], y[2]) && _identical(x[3], y[3])
_identical(x::SuperGaussianSpot, y::SuperGaussianSpot) =
    _identical(x.a, y.a) && _identical(x.b, y.b) && x.order == y.order && x.angle == y.angle && x.halo == y.halo &&
    x.halo_scale == y.halo_scale && x.halo_order == y.halo_order

_beam_names(ev::PerPhase) = _beam_names(ev.values[1])

check_constraints(ev::PerPhase, c::Constraint; kwargs...) = check_constraints(ev, (c,); kwargs...)
check_constraints(ev::PerPhase, cs; feasibility::Real = 1e-3) =
    _constraint_rows(ev, Tuple(cs), nothing, Float64(feasibility))

"""
    pointing_summary(scenario::PhasedScenario, beams::PerPhase; reference = the beams of the start,
                     keepout = nothing, inner = (:inner23, :inner30), drive_reference = nothing,
                     lmax = 4) -> PerPhase
    pointing_summary(scenario::PhasedScenario, result::PointingResult; kwargs...)

A [`PointingSummary`](@ref) per phase, of the beams of each phase in its scenario against
the `reference` beams. `reference` is a `PerPhase` of beam vectors, by default the start
of the scenario `apply(scenario.param, initial_parameters(scenario.param))`, or one beam
vector for every phase. The drive is compared with `drive_reference[phase]` where given,
as a `PerPhase`, a named tuple or a dictionary by phase of evaluations or coefficient
vectors; a phase without one has no drive ratio. The display is one line per phase.

`gates` and `floor` show the departure of each phase against its gates and its pointing
floor, as for one [`Scenario`](@ref). They are given per phase (a `PerPhase`, a named
tuple or a dictionary by phase), or as one value for every phase (a vector, a
[`mode_uncertainty`](@ref) result, a [`mode_gate`](@ref)). Nothing is shown for a phase
without an entry. `gates` may also be the constraints of the run, a vector or tuple such
as the `constraints` passed to [`optimize_pointing`](@ref). The summary of each phase
then shows its `in_phase(phase, mode_gate(...))`, from the reference shape of that gate.
`moved` and `moved_by` count the groups moved past each shift threshold, as for one
scenario. The phases share the aims, so the counts are the same in every phase.

```julia
pointing_summary(ps, r; drive_reference = refs, gates = constraints, floor = floors)
```
"""
function pointing_summary(ps::PhasedScenario, beams::PerPhase;
                          reference = apply(ps.param, initial_parameters(ps.param)),
                          keepout = nothing, inner = (:inner23, :inner30),
                          drive_reference = nothing, lmax::Integer = 4, gates = nothing,
                          floor = nothing, moved = (), moved_by = :quad)
    beams.names == ps.param.phases ||
        throw(ArgumentError("the beams' phases are not the scenario's ($(join(ps.param.phases, ", ")))"))
    out = map(eachindex(ps.scenarios)) do k
        ph = ps.param.phases[k]
        ref = reference isa PerPhase ? reference[ph] : reference
        dr = drive_reference === nothing ? nothing : haskey(drive_reference, ph) ? drive_reference[ph] : nothing
        pointing_summary(ps.scenarios[k], beams[ph]; reference = ref, keepout, inner, drive_reference = dr, lmax,
                         gates = _phase_gate(gates, ph), floor = _phase_entry(floor, ph), moved, moved_by)
    end
    return PerPhase(ps.param.phases, out)
end

# Returns the entry of one phase in a per-phase option of a summary. A value for every
# phase is returned as it is; otherwise the entry of the phase is returned, or `nothing`
# without one.
_phase_entry(::Nothing, ph) = nothing
_phase_entry(x::Union{AbstractVector{<:Real},Constraint}, ph) = x
_phase_entry(x, ph) = x isa NamedTuple && haskey(x, :rms_mode) ? x : haskey(x, ph) ? x[ph] : nothing

# Returns the gates of a phase. From the constraints of a run, the gates are the
# `in_phase` mode gate of the phase.
_phase_gate(x, ph) = _phase_entry(x, ph)
function _phase_gate(cs::Union{Tuple{Vararg{Constraint}},AbstractVector{<:Constraint}}, ph)
    i = findfirst(c -> c.f isa _InPhase && c.f.phase === ph && c.f.f isa _GateValue, cs)
    return i === nothing ? nothing : cs[i]
end

_gate_value(f::_InPhase, ev::Evaluation) = _gate_value(f.f, ev)

pointing_summary(ps::PhasedScenario, r::PointingResult; kwargs...) = pointing_summary(ps, r.beams; kwargs...)
