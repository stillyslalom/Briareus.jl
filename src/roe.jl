# Beam rules of engagement: which beams and quads a shot can lose. An outage is a set of
# beams at zero power. Every outage is scored against a list of criteria, each a function
# of an `Evaluation` with a green and an amber threshold, so the engine holds no model of
# the target: a hohlraum, a halfraum, a foil and a bare sphere differ only in the criteria.

const _ROE_MEASURES = (:value, :change, :relative, :loss, :rise)
const _ROE_TIERS = (:green, :amber, :red)
_roe_rank(tier::Symbol) = tier === :green ? 1 : tier === :amber ? 2 : 3

"""
    ROECriterion(name, f; green, amber, measure = :change, case = :default, label = "")

One criterion of a beam rules-of-engagement scan ([`roe_scan`](@ref)). `f` maps an
[`Evaluation`](@ref) of the shot with some beams at zero power to a number or a vector,
such as `ev -> flux_moments(ev)[3] / flux_moments(ev)[1]`. The score of an outage is
taken from `v = f(outage)` and `v₀ = f(nominal)` by `measure`:

- `:value` is `v`, or its norm for a vector;
- `:change` is `|v − v₀|`, the norm of the difference for a vector;
- `:relative` is `|v − v₀| / |v₀|`;
- `:loss` is `(v₀ − v) / |v₀|`, the fractional loss of a scalar;
- `:rise` is `v − v₀` of a scalar.

An outage is green on the criterion when the score is at most `green`, amber when it is
at most `amber`, and red above that or when the score is not a number. `green ≤ amber`
is required. `case` names the case of the scan that the criterion is evaluated on, by
default the first. `label` is the unit of the score in the exported tables.

`f` runs on several threads at once and must not write shared state. Any function of an
evaluation serves, so quantities outside the built-in constructors (`roe_drive`,
`roe_mode`, `roe_hotspot_velocity`, `roe_face_flux`, `roe_foil_irradiance`, …) enter the
same scan, tiers and tables.
"""
struct ROECriterion
    name::Symbol
    f::Any
    green::Float64
    amber::Float64
    measure::Symbol
    case::Symbol
    label::String
end

function ROECriterion(name, @nospecialize(f); green::Real, amber::Real, measure::Symbol = :change,
                      case::Symbol = :default, label::AbstractString = "")
    measure in _ROE_MEASURES ||
        throw(ArgumentError("measure must be one of $(join(_ROE_MEASURES, ", ")), got :$measure"))
    green <= amber || throw(ArgumentError("the green threshold $green exceeds the amber threshold $amber"))
    return ROECriterion(Symbol(name), f, Float64(green), Float64(amber), measure, case, String(label))
end

Base.show(io::IO, c::ROECriterion) =
    print(io, "ROECriterion(:", c.name, ", ", c.measure, ", green ≤ ", c.green, ", amber ≤ ", c.amber,
          isempty(c.label) ? "" : " " * c.label, c.case === :default ? "" : ", case :$(c.case)", ")")

_roe_size(v::Real) = abs(float(v))
_roe_size(v) = float(norm(v))

function _roe_score(c::ROECriterion, @nospecialize(v), @nospecialize(v0))
    m = c.measure
    m === :value && return Float64(v isa Real ? float(v) : norm(v))
    m === :change && return Float64(_roe_size(v - v0))
    m === :relative && return Float64(_roe_size(v - v0) / _roe_size(v0))
    (v isa Real && v0 isa Real) ||
        throw(ArgumentError("criterion :$(c.name): measure :$m needs a scalar value"))
    m === :loss && return Float64((v0 - v) / abs(v0))
    return Float64(v - v0)
end

_roe_tier(c::ROECriterion, score::Float64) =
    score <= c.green ? :green : score <= c.amber ? :amber : :red

"""
    ROEUnit

A droppable unit of a rules-of-engagement scan: `name`, `kind` (`:beam` or `:quad`),
`beams` (indices into the beam list of the scan), `role` and `share`, the fraction of the
power of those beams that the unit drops. The share is 1 except for one facility beam of
a beam that stands for several (`stands_for` of [`roe_scan`](@ref)), such as one beam of
a NIF quad held as one beam, which drops a quarter of it.
"""
struct ROEUnit
    name::String
    kind::Symbol
    beams::Vector{Int}
    role::Symbol
    share::Float64
end

ROEUnit(name, kind, beams, role) = ROEUnit(name, kind, beams, role, 1.0)

"""
    ROEOutage

The result of one outage of a rules-of-engagement scan, with fields:

- `units`, the indices of the dropped units in `report.units`;
- `beams`, the indices of the beams that lose power, and `shares`, the fraction of the
  power that each of them loses (1 for a beam at zero power);
- `scores` and `tiers`, one per criterion of the scan;
- `tier`, the worst of `tiers`;
- `binding`, the index of the criterion that sets `tier`, the one farthest past its
  thresholds.
"""
struct ROEOutage
    units::Vector{Int}
    beams::Vector{Int}
    shares::Vector{Float64}
    scores::Vector{Float64}
    tiers::Vector{Symbol}
    tier::Symbol
    binding::Int
end

"""
    ROEBudget

The multiple-outage rule of one beam role, from the pair and triple scans: `role`,
`kind`, the number `green` of green units, the number `pairs` of scanned pairs of green
units, the number `safe` of those pairs that stay green, the number `triples` of scanned
triples of green units and the number `safe_triples` of those that stay green, `limit`
(how many green units may drop together without a request: 0 to 3), `exceptions` and
`triple_exceptions` (the pairs and triples of green units that leave green, by name) and
`text`, the rule as a sentence. Triples are scanned only with `triples = true` in
[`roe_scan`](@ref), and four or more simultaneous outages are not scanned, so `limit` is
at most 3.
"""
struct ROEBudget
    role::Symbol
    kind::Symbol
    green::Int
    pairs::Int
    safe::Int
    triples::Int
    safe_triples::Int
    limit::Int
    exceptions::Vector{Tuple{String,String}}
    triple_exceptions::Vector{NTuple{3,String}}
    text::String
end

"""
    ROEReport

The result of [`roe_scan`](@ref), with fields:

- `criteria`, the criteria of the scan, and `baseline`, their values on the full shot;
- `cases`, the case names, and `linear`, whether each case was scanned through the
  linear deposition tables;
- `beams`, the beam names;
- `units`, the droppable beams and quads ([`ROEUnit`](@ref));
- `singles`, the outage of each unit alone ([`ROEOutage`](@ref)), in the order of `units`;
- `pairs`, the scanned outages of two units, and `triples`, those of three units (empty
  unless the scan was asked for them);
- `companions`, per unit the outage of the unit together with its companion, or `nothing`;
- `tiers`, the tier of each unit: the tier of its single outage, raised from green to
  amber when a pair with another green unit is red;
- `caveats`, per unit the indices into `pairs` of the pairs that are worse than both of
  their units alone;
- `budgets`, the multiple-outage rule of each role ([`ROEBudget`](@ref)).

[`roe_table`](@ref), [`write_roe_csv`](@ref) and [`roe_markdown`](@ref) export it.
"""
struct ROEReport
    criteria::Vector{ROECriterion}
    cases::Vector{Symbol}
    linear::Vector{Bool}
    beams::Vector{String}
    units::Vector{ROEUnit}
    baseline::Vector{Any}
    singles::Vector{ROEOutage}
    pairs::Vector{ROEOutage}
    triples::Vector{ROEOutage}
    companions::Vector{Union{Nothing,ROEOutage}}
    tiers::Vector{Symbol}
    caveats::Vector{Vector{Int}}
    budgets::Vector{ROEBudget}
end

function Base.show(io::IO, r::ROEReport)
    n(k, t) = count(i -> r.units[i].kind === k && r.tiers[i] === t, eachindex(r.units))
    print(io, "ROEReport(", length(r.criteria), " criteria")
    for k in unique(u.kind for u in r.units)
        print(io, "; ", count(u -> u.kind === k, r.units), " ", k, "s: ", n(k, :red), " red, ",
              n(k, :amber), " amber, ", n(k, :green), " green")
    end
    print(io, "; ", length(r.pairs), " pairs")
    isempty(r.triples) || print(io, ", ", length(r.triples), " triples")
    print(io, ")")
end

# One case of a scan: a scenario with the beams of its epoch, the nominal evaluation and
# the deposition of each beam by itself. The tables are the columns of the linear path.
struct _ROECase
    name::Symbol
    scenario::Scenario
    beams::Vector{Beam{Float64}}
    effective::Vector{Beam{Float64}}
    nominal::Evaluation
    element::Matrix{Float64}    # elements × beams: absorbed power of each beam alone
    backside::Matrix{Float64}   # components × beams
    linear::Bool
end

# Builds a case. The tables hold the deposition of each effective beam with every other
# beam at zero power, which `deposit` skips, so the columns cost one deposition together.
function _roe_case(name::Symbol, s::Scenario, beams, method::Symbol)
    bs = collect(Beam{Float64}, beams)
    lpi = _resolve_lpi(s.lpi, bs)
    eff = collect(Beam{Float64}, apply_lpi(lpi, bs, s.time))
    nominal = evaluate(s, bs)
    ne, nc, nb = nelements(s.target), length(s.target.components), length(bs)
    E, L = zeros(ne, nb), zeros(nc, nb)
    off = [remake(b; power = 0.0) for b in eff]
    for b in 1:nb
        eff[b].power > 0 || continue
        only_b = copy(off)
        only_b[b] = eff[b]
        d = lpi_deposit(lpi, s.deposition, only_b, s.time)
        E[:, b] = d.element
        L[:, b] = d.backside
    end
    linear = method === :linear || (method === :auto && s.lpi === nothing)
    return _ROECase(name, s, bs, eff, nominal, E, L, linear)
end

# The evaluation of a case with the beams `drop` less the fractions `share` of their
# power (1: at zero power). On the linear path the deposition is the nominal one less
# those fractions of the columns of the dropped beams, floored at zero where the
# difference of rounded sums would leave a negative power. Otherwise the shot is
# evaluated again, laser–plasma model included.
_roe_left(p, share) = share == 1 ? zero(p) : p * (1 - share)
function _roe_evaluation(c::_ROECase, drop::Vector{Int}, share::Vector{Float64})
    beams = copy(c.beams)
    for (b, f) in zip(drop, share)
        beams[b] = remake(beams[b]; power = _roe_left(beams[b].power, f))
    end
    c.linear || return evaluate(c.scenario, beams)
    eff = copy(c.effective)
    d0 = c.nominal.deposition
    E, L, per_beam = copy(d0.element), copy(d0.backside), copy(d0.per_beam)
    for (b, f) in zip(drop, share)
        eff[b] = remake(eff[b]; power = _roe_left(eff[b].power, f))
        if f == 1
            E .-= @view c.element[:, b]
            L .-= @view c.backside[:, b]
        else
            E .-= f .* @view c.element[:, b]
            L .-= f .* @view c.backside[:, b]
        end
        per_beam[b] = _roe_left(per_beam[b], f)
    end
    E .= max.(E, 0.0)
    L .= max.(L, 0.0)
    dep = Deposition(E, L, per_beam, sum(b.power for b in eff))
    return _evaluation(c.scenario, beams, eff, dep, nothing)
end

function _roe_outage(cases::Vector{_ROECase}, criteria::Vector{ROECriterion}, case_of::Vector{Int},
                     base::Vector{Any}, units::Vector{Int}, beams::Vector{Int}, share::Vector{Float64})
    evs = Vector{Union{Nothing,Evaluation}}(nothing, length(cases))
    scores = Vector{Float64}(undef, length(criteria))
    for (k, c) in enumerate(criteria)
        j = case_of[k]
        evs[j] === nothing && (evs[j] = _roe_evaluation(cases[j], beams, share))
        scores[k] = _roe_score(c, c.f(evs[j]::Evaluation), base[k])
    end
    return _roe_result(criteria, units, beams, share, scores)
end

function _roe_result(criteria::Vector{ROECriterion}, units, beams, share, scores::Vector{Float64})
    tiers = Symbol[_roe_tier(c, s) for (c, s) in zip(criteria, scores)]
    binding, worst = 0, (0, -Inf)
    for (k, c) in enumerate(criteria)
        scale = isfinite(c.amber) && c.amber > 0 ? c.amber : 1.0
        key = (_roe_rank(tiers[k]), isnan(scores[k]) ? Inf : scores[k] / scale)
        key > worst && ((binding, worst) = (k, key))
    end
    tier = isempty(tiers) ? :green : _ROE_TIERS[maximum(_roe_rank, tiers)]
    return ROEOutage(collect(Int, units), collect(Int, beams), collect(Float64, share), scores, tiers, tier, binding)
end

# This function is the barrier for the threaded loop over outages.
function _roe_outages(cases::Vector{_ROECase}, criteria::Vector{ROECriterion}, case_of::Vector{Int},
                      base::Vector{Any}, list::Vector{Tuple{Vector{Int},Vector{Int},Vector{Float64}}})
    out = Vector{ROEOutage}(undef, length(list))
    Threads.@threads for i in eachindex(list)
        out[i] = _roe_outage(cases, criteria, case_of, base, list[i]...)
    end
    return out
end

# The role of each beam. `:part` is the part that absorbs most of the beam in the first
# case in which the beam has power, and `:none` for a beam that lands on nothing.
function _roe_roles(@nospecialize(roles), cases::Vector{_ROECase}, names::Vector{String})
    nb = length(names)
    auto() = Symbol[_roe_part_role(cases, b) for b in 1:nb]
    roles === :part && return auto()
    roles === :group && return Symbol[cases[1].beams[b].group for b in 1:nb]
    if roles isa AbstractVector{Symbol}
        length(roles) == nb || throw(DimensionMismatch("roles has $(length(roles)) entries for $nb beams"))
        return collect(Symbol, roles)
    end
    if roles isa AbstractDict
        a = auto()
        for (k, v) in roles
            i = findfirst(==(String(k)), names)
            i === nothing && throw(ArgumentError("roles names the unknown beam $k"))
            a[i] = Symbol(v)
        end
        return a
    end
    roles isa Symbol && throw(ArgumentError("roles must be :part, :group, a vector, a dictionary or a function"))
    return Symbol[Symbol(roles(b)) for b in cases[1].beams]
end

function _roe_part_role(cases::Vector{_ROECase}, b::Int)
    for c in cases
        c.effective[b].power > 0 || continue
        t = c.scenario.target
        P = Dict{Symbol,Float64}()
        for e in axes(c.element, 1)
            p = c.element[e, b]
            p > 0 || continue
            part = t.components[t.elem_comp[e]].part
            P[part] = get(P, part, 0.0) + p
        end
        isempty(P) && return :none
        best = first(parts(t))
        for p in parts(t)   # in target order, so ties resolve the same way every run
            get(P, p, 0.0) > get(P, best, 0.0) && (best = p)
        end
        return best
    end
    return :none
end

function _roe_mode(v::Vector{Symbol})
    best, n = v[1], 0
    for x in v
        k = count(==(x), v)
        k > n && ((best, n) = (x, k))
    end
    return best
end

# The beams that the units `idx` drop together and the share of each, by beam index.
function _roe_drop(us::Vector{ROEUnit}, idx)
    bs, fs = Int[], Float64[]
    for i in idx
        append!(bs, us[i].beams)
        append!(fs, fill(us[i].share, length(us[i].beams)))
    end
    p = sortperm(bs)
    return bs[p], fs[p]
end

# The number of facility beams that each beam stands for (`stands_for` of `roe_scan`).
function _roe_stands_for(@nospecialize(stands_for), beams::Vector{Beam{Float64}})
    stands_for === nothing && return ones(Int, length(beams))
    n = if stands_for isa AbstractDict
        names = Set(b.name for b in beams)
        for k in keys(stands_for)
            String(k) in names || throw(ArgumentError("stands_for names the unknown beam $k"))
        end
        d = Dict(String(k) => v for (k, v) in stands_for)
        Int[get(d, b.name, 1) for b in beams]
    else
        Int[stands_for(b) for b in beams]
    end
    all(>=(1), n) || throw(ArgumentError("a beam stands for at least one facility beam"))
    return n
end

# The companion of each unit as a unit index, or 0.
function _roe_companions(@nospecialize(companions), units::Vector{ROEUnit}, ports::Vector{Vec3}, axis::Vec3)
    comp = zeros(Int, length(units))
    companions === nothing && return comp
    mean_port(u) = normalize(sum(ports[b] for b in u.beams))
    index(name) = (i = findfirst(u -> u.name == String(name), units);
                   i === nothing ? throw(ArgumentError("companions names the unknown unit $name")) : i)
    if companions === :opposite || companions === :mirror
        P = [mean_port(u) for u in units]
        for (i, u) in enumerate(units)
            want = companions === :opposite ? -P[i] : P[i] - 2 * dot(P[i], axis) * axis
            best, d = 0, -Inf
            for (j, v) in enumerate(units)
                (j == i || v.kind !== u.kind || v.role !== u.role) && continue
                dot(want, P[j]) > d && ((best, d) = (j, dot(want, P[j])))
            end
            comp[i] = best
        end
    elseif companions isa Symbol
        throw(ArgumentError("companions must be :opposite, :mirror, pairs of unit names or a function"))
    elseif companions isa Union{AbstractDict,AbstractVector,Tuple}
        for (a, b) in companions
            comp[index(a)] = index(b)
        end
    else
        for (i, u) in enumerate(units)
            c = companions(u.name)
            c === nothing || (comp[i] = index(c))
        end
    end
    any(i -> comp[i] == i, eachindex(comp)) && throw(ArgumentError("a unit cannot be its companion"))
    return comp
end

_roe_cases(c::NamedTuple) = [k => v for (k, v) in pairs(c)]
_roe_cases(c::AbstractDict) = sort!([Symbol(k) => v for (k, v) in c]; by = first)
_roe_cases(c) = [Symbol(first(p)) => last(p) for p in c]

"""
    roe_scan(scenario, beams, criteria; kwargs...) -> ROEReport
    roe_scan(cases, criteria; units = :auto, pairs = :auto, pair_scope = :all,
             extra_pairs = (), triples = false, roles = :part, companions = nothing,
             axis = (0, 0, 1), stands_for = nothing, method = :auto)

Beam rules of engagement: scores the loss of every beam, every quad and pairs of them
against `criteria` (a vector of [`ROECriterion`](@ref)) and sorts each beam and quad into
a tier. Green units may drop without a request, amber ones need the approval of the
responsible individual for the shot, and red ones hold the shot. The scan takes criteria
and assumes nothing about the target, so one call covers a capsule in a hohlraum, the
face of a package, a backlighter foil and every other quantity that some function of an
[`Evaluation`](@ref) returns.

`cases` are pairs `name => (scenario, beams)`, as a vector, a named tuple or a
dictionary. Each case is one view of the shot: an epoch of the pulse with the beam
powers of that epoch, or a scenario that observes another face. The beams of every case
are the same beams in the same order, with powers that may differ, and a beam with no
power in any case is left out. The first form is the single case `:default`.

An outage is a set of beams at zero power, with every other beam as planned (no
repointing and no power change). The scan covers:

- each beam with power, and each quad (`beam.quad`) of more than one beam, as selected
  by `units`: `:auto` for both, `:beams` or `:quads`. A beam that stands for several
  facility beams, by `stands_for` (a dictionary from beam name to count, or a function
  of the beam), is a quad unit named by its quad (by the beam when `beam.quad` is
  empty), and its beam unit, named such as `"Q31B (1 of 4)"`, drops one of those beams
  as that fraction of its power: a NIF quad held as one beam (`nif_quads`) stands for 4.
  The beams of a quad leave from ports within a few degrees of the quad port, so this
  beam unit approximates the loss of each of them by a share of the quad;
- pairs of units, by `pairs`: `:quads`, `:beams`, `:both` or `:none`, and `:auto` for the
  quad pairs when the beams form quads and the beam pairs otherwise. `pair_scope` keeps
  `:all` pairs, those `:within` one role or those `:across` roles. `extra_pairs` adds
  pairs of unit names;
- with `triples = true`, the triples of green units of one role and kind in a role whose
  green units may all drop in pairs (every pair of them scanned and green), so that the
  budget can allow a third drop. A role of `g` green units adds `g(g - 1)(g - 2)/6`
  outages, about 20,000 for 50 units;
- with `companions`, each unit together with its companion, a pre-planned second drop
  meant to restore the balance. `:opposite` pairs a unit with the unit of the same kind
  and role whose mean port direction is nearest the antipode, `:mirror` with the one nearest the
  reflection through the plane normal to `axis`, and pairs of names or a function of the
  unit name give the companions explicitly.

`roles` groups the beams in the tables and the budgets. `:part` gives each beam the part
that absorbs most of its power, so the drive beams of a hohlraum, the beams of a second
halfraum and the beams of a backlighter foil separate without input. `:group` uses
`beam.group`. A vector of symbols, a dictionary from beam name to role (other beams keep
their part) or a function of the beam assigns roles explicitly. A quad has the role of
most of its beams.

Deposition without a laser–plasma model is linear in the beam powers. With `method =
:auto`, a case whose scenario has no such model is scanned through tables of the
deposition of each beam alone: the deposition of an outage is the nominal deposition
less the columns of the dropped beams, and the tables cost one deposition together. The
criteria are evaluated in full on that deposition, so nonlinear quantities (a ratio of
moments, a temperature, a conversion efficiency that depends on irradiance) need no
linearization. A case with a laser–plasma model is evaluated again for every outage
(`method = :full`), since a model such as `CBET` couples the surviving beams to the
dropped ones. `method = :linear` forces the tables with the effective beams of the full
shot, which is exact for models that act on each beam separately and freezes the
coupling otherwise; `method = :full` on every case is the independent check of a linear
scan.

The tier of a unit is the tier of its single outage. A green unit whose pair with
another green unit is red is raised to amber, because a drop without a request must stay
safe after another such drop. The other pairs that are worse than both of their units
are recorded as caveats, and the budget of each role ([`ROEBudget`](@ref)) gives the
number of green units that may drop together: two when the pairs stay green, and three
when the triples are scanned and stay green.

The outages run on all threads, so the criteria must be safe to run concurrently.
"""
roe_scan(cases, criteria::AbstractVector{ROECriterion}; kwargs...) = _roe_scan(cases, criteria; kwargs...)

# The scan compiles once: the criteria are untyped fields, and the role and companion
# functions are not specialized on.
Base.@nospecializeinfer function _roe_scan(@nospecialize(cases), criteria::AbstractVector{ROECriterion};
                                           units::Symbol = :auto, pairs::Symbol = :auto,
                                           pair_scope::Symbol = :all, @nospecialize(extra_pairs = ()),
                                           triples::Bool = false,
                                           @nospecialize(roles = :part), @nospecialize(companions = nothing),
                                           axis = (0.0, 0.0, 1.0), @nospecialize(stands_for = nothing),
                                           method::Symbol = :auto)
    units in (:auto, :beams, :quads) || throw(ArgumentError("units must be :auto, :beams or :quads"))
    pairs in (:auto, :quads, :beams, :both, :none) ||
        throw(ArgumentError("pairs must be :auto, :quads, :beams, :both or :none"))
    pair_scope in (:all, :within, :across) ||
        throw(ArgumentError("pair_scope must be :all, :within or :across"))
    method in (:auto, :linear, :full) || throw(ArgumentError("method must be :auto, :linear or :full"))
    isempty(criteria) && throw(ArgumentError("a scan needs at least one criterion"))
    allunique(c.name for c in criteria) || throw(ArgumentError("criterion names must be unique"))
    crit = collect(ROECriterion, criteria)
    cs = _ROECase[]
    for (name, sb) in _roe_cases(cases)
        (sb isa Tuple && length(sb) == 2 && sb[1] isa Scenario) ||
            throw(ArgumentError("case :$name must be a (scenario, beams) tuple"))
        push!(cs, _roe_case(name, sb[1], sb[2], method))
    end
    isempty(cs) && throw(ArgumentError("a scan needs at least one case"))
    allunique(c.name for c in cs) || throw(ArgumentError("case names must be unique"))
    names = String[b.name for b in cs[1].beams]
    allunique(names) || throw(ArgumentError("beam names must be unique"))
    for c in cs
        String[b.name for b in c.beams] == names ||
            throw(ArgumentError("case :$(c.name) has other beams than case :$(cs[1].name)"))
    end
    case_of = Int[]
    for c in crit
        j = c.case === :default ? 1 : findfirst(x -> x.name === c.case, cs)
        j === nothing && throw(ArgumentError("criterion :$(c.name) names the unknown case :$(c.case)"))
        push!(case_of, j)
    end
    base = Any[c.f(cs[j].nominal) for (c, j) in zip(crit, case_of)]

    # The units: beams with power in some case, and quads of more than one such beam.
    # A beam that stands for several facility beams is a quad, and its beam unit drops
    # one of them.
    nb = length(names)
    active = [b for b in 1:nb if any(c -> c.beams[b].power > 0, cs)]
    isempty(active) && throw(ArgumentError("no beam has power"))
    role = _roe_roles(roles, cs, names)
    held = _roe_stands_for(stands_for, cs[1].beams)
    us = ROEUnit[]
    if units !== :quads
        for b in active
            n = held[b]
            push!(us, n == 1 ? ROEUnit(names[b], :beam, [b], role[b]) :
                      ROEUnit(string(names[b], " (1 of ", n, ")"), :beam, [b], role[b], 1 / n))
        end
    end
    quads = String[]
    for b in active
        q = cs[1].beams[b].quad
        held[b] == 1 && !isempty(q) && !(q in quads) && push!(quads, q)
    end
    nquad = 0
    if units !== :beams
        for b in active
            held[b] > 1 || continue
            q = cs[1].beams[b].quad
            push!(us, ROEUnit(isempty(q) ? names[b] : String(q), :quad, [b], role[b]))
            nquad += 1
        end
        for q in quads
            bs = [b for b in active if held[b] == 1 && cs[1].beams[b].quad == q]
            length(bs) > 1 || continue
            push!(us, ROEUnit(q, :quad, bs, _roe_mode(role[bs])))
            nquad += 1
        end
    end
    isempty(us) && throw(ArgumentError("units = :$units leaves nothing to drop"))
    allunique(u.name for u in us) || throw(ArgumentError("a quad and a beam share a name"))

    kinds = pairs === :none ? Symbol[] : pairs === :both ? [:beam, :quad] :
            pairs === :auto ? [nquad > 0 ? :quad : :beam] : [pairs === :quads ? :quad : :beam]
    list = Tuple{Vector{Int},Vector{Int},Vector{Float64}}[([i], _roe_drop(us, (i,))...) for i in eachindex(us)]
    npair = 0
    seen = Set{Tuple{Int,Int}}()
    function add_pair(i, j)
        (i, j) = minmax(i, j)
        (i == j || (i, j) in seen || !isempty(intersect(us[i].beams, us[j].beams))) && return
        push!(seen, (i, j))
        push!(list, ([i, j], _roe_drop(us, (i, j))...))
        npair += 1
    end
    for i in eachindex(us), j in i+1:length(us)
        (us[i].kind === us[j].kind && us[i].kind in kinds) || continue
        same = us[i].role === us[j].role
        (pair_scope === :all || (pair_scope === :within) == same) && add_pair(i, j)
    end
    uindex(n) = (i = findfirst(u -> u.name == String(n), us);
                 i === nothing ? throw(ArgumentError("extra_pairs names the unknown unit $n")) : i)
    for (a, b) in extra_pairs
        add_pair(uindex(a), uindex(b))
    end
    comp = _roe_companions(companions, us, [b.port for b in cs[1].beams], normalize(Vec3(axis...)))
    cslot = zeros(Int, length(us))
    for (i, j) in enumerate(comp)
        j == 0 && continue
        isempty(intersect(us[i].beams, us[j].beams)) || continue
        push!(list, ([i, j], _roe_drop(us, (i, j))...))
        cslot[i] = length(list)
    end

    out = _roe_outages(cs, crit, case_of, base, list)
    nu = length(us)
    singles, pairouts = out[1:nu], out[nu+1:nu+npair]
    companion = Union{Nothing,ROEOutage}[k == 0 ? nothing : out[k] for k in cslot]
    tiers, caveats = _roe_tiers(singles, pairouts)
    budgets = _roe_budgets(us, tiers, singles, pairouts)
    triouts = ROEOutage[]
    if triples
        tlist = [(g, _roe_drop(us, g)...) for g in _roe_triple_sets(us, tiers, budgets)]
        triouts = _roe_outages(cs, crit, case_of, base, tlist)
        budgets = _roe_triple_budgets(budgets, us, tiers, triouts)
    end
    return ROEReport(crit, [c.name for c in cs], [c.linear for c in cs], names, us, base, singles,
                     pairouts, triouts, companion, tiers, caveats, budgets)
end

roe_scan(s::Scenario, beams::AbstractVector{<:Beam}, criteria::AbstractVector{ROECriterion}; kwargs...) =
    _roe_scan(Pair{Symbol,Any}[:default => (s, beams)], criteria; kwargs...)

"""
    roe_outage(report, scenario, beams, drop) -> ROEOutage
    roe_outage(report, cases, drop) -> ROEOutage

A hypothetical outage of any number of units, scored against the criteria and the
baseline of `report` (from [`roe_scan`](@ref) on the same `scenario` and `beams`, or the
same `cases`). `drop` lists the names of units of the report (beams and quads) and of
beams of the report that are not units, such as beams with no power. The beams of the
dropped units lose the shares of their power that the units drop (all of it, except for
the beam unit of a beam that stands for several), and the shares of units on one beam
add up to at most all of it. Every case of a criterion is evaluated
again in full, laser–plasma model included, so the result is that of a scan with
`method = :full`. The scan covers single units, pairs and, when asked, triples; this
function covers any other set of units.

The returned outage lists the dropped units and beams by their indices in the report.
"""
function roe_outage(r::ROEReport, s::Scenario, beams::AbstractVector{<:Beam}, drop)
    return roe_outage(r, Pair{Symbol,Any}[:default => (s, beams)], drop)
end

function roe_outage(r::ROEReport, @nospecialize(cases), @nospecialize(drop))
    cs = _roe_cases(cases)
    Symbol[first(p) for p in cs] == r.cases ||
        throw(ArgumentError("the cases must be those of the report, " * join((":" * String(c) for c in r.cases), ", ")))
    units = Int[]
    share = Dict{Int,Float64}()
    for name in unique(String(n) for n in drop)
        i = findfirst(u -> u.name == name, r.units)
        if i === nothing
            b = findfirst(==(name), r.beams)
            b === nothing && throw(ArgumentError("$name is neither a unit nor a beam of the report"))
            share[b] = 1.0
        else
            push!(units, i)
            for b in r.units[i].beams
                share[b] = min(1.0, get(share, b, 0.0) + r.units[i].share)
            end
        end
    end
    isempty(share) && throw(ArgumentError("an outage drops at least one beam"))
    sort!(unique!(units))
    bs = sort!(collect(keys(share)))
    fs = Float64[share[b] for b in bs]
    evs = Vector{Union{Nothing,Evaluation}}(nothing, length(cs))
    function ev(j)
        evs[j] === nothing || return evs[j]::Evaluation
        sb = last(cs[j])
        (sb isa Tuple && length(sb) == 2 && sb[1] isa Scenario) ||
            throw(ArgumentError("case :$(first(cs[j])) must be a (scenario, beams) tuple"))
        b2 = collect(Beam{Float64}, sb[2])
        String[b.name for b in b2] == r.beams ||
            throw(ArgumentError("case :$(first(cs[j])) has other beams than the report"))
        for (b, f) in zip(bs, fs)
            b2[b] = remake(b2[b]; power = _roe_left(b2[b].power, f))
        end
        return evs[j] = evaluate(sb[1], b2)
    end
    scores = Vector{Float64}(undef, length(r.criteria))
    for (k, c) in enumerate(r.criteria)
        j = c.case === :default ? 1 : findfirst(==(c.case), r.cases)
        scores[k] = _roe_score(c, c.f(ev(j)), r.baseline[k])
    end
    return _roe_result(r.criteria, units, bs, fs, scores)
end

# The final tier of each unit and its pair caveats.
function _roe_tiers(singles::Vector{ROEOutage}, pairouts::Vector{ROEOutage})
    tiers = Symbol[o.tier for o in singles]
    caveats = [Int[] for _ in singles]
    for (k, p) in enumerate(pairouts)
        i, j = p.units
        ri, rj, rp = _roe_rank(singles[i].tier), _roe_rank(singles[j].tier), _roe_rank(p.tier)
        rp > max(ri, rj) || continue
        push!(caveats[i], k)
        push!(caveats[j], k)
        if ri == 1 && rj == 1 && rp == 3
            tiers[i] = tiers[j] = :amber
        end
    end
    return tiers, caveats
end

# The budget of each role and kind, from the pairs of units that are green alone.
function _roe_budgets(us::Vector{ROEUnit}, tiers::Vector{Symbol}, singles::Vector{ROEOutage},
                      pairouts::Vector{ROEOutage})
    out = ROEBudget[]
    for kind in unique(u.kind for u in us), role in unique(u.role for u in us if u.kind === kind)
        member(i) = us[i].kind === kind && us[i].role === role
        green = [i for i in eachindex(us) if member(i) && tiers[i] === :green]
        scanned = [p for p in pairouts if all(member, p.units) && all(i -> singles[i].tier === :green, p.units)]
        ok = [p for p in scanned if all(i -> tiers[i] === :green, p.units)]
        safe = count(p -> p.tier === :green, ok)
        exceptions = [(us[p.units[1]].name, us[p.units[2]].name) for p in ok if p.tier !== :green]
        word = string(kind, "s")
        npossible = length(green) * (length(green) - 1) ÷ 2
        limit, text = if isempty(green)
            0, "no $kind drops without a request"
        elseif length(green) == 1
            1, "the one green $kind may drop without a request"
        elseif length(ok) < npossible
            1, "1 green $kind at a time ($(length(ok)) of $npossible pairs scanned)"
        elseif safe == 0
            1, "1 green $kind at a time (every pair of green $word leaves green)"
        elseif safe == length(ok)
            2, "up to 2 green $word together (all $safe pairs stay green; 3 or more not scanned)"
        else
            2, "up to 2 green $word together, except $(length(ok) - safe) of $(length(ok)) pairs (3 or more not scanned)"
        end
        push!(out, ROEBudget(role, kind, length(green), length(ok), safe, 0, 0, limit, exceptions,
                             NTuple{3,String}[], text))
    end
    return out
end

# The triples of green units of each role and kind whose budget allows any two of them
# together, with every pair scanned.
function _roe_triple_sets(us::Vector{ROEUnit}, tiers::Vector{Symbol}, budgets::Vector{ROEBudget})
    out = Vector{Int}[]
    for b in budgets
        (b.limit == 2 && b.safe == b.pairs == b.green * (b.green - 1) ÷ 2 && b.green >= 3) || continue
        g = [i for i in eachindex(us) if us[i].kind === b.kind && us[i].role === b.role && tiers[i] === :green]
        for x in 1:length(g), y in x+1:length(g), z in y+1:length(g)
            push!(out, [g[x], g[y], g[z]])
        end
    end
    return out
end

# The budgets with the triples: three green units together when the triples stay green.
function _roe_triple_budgets(budgets::Vector{ROEBudget}, us::Vector{ROEUnit}, tiers::Vector{Symbol},
                             triouts::Vector{ROEOutage})
    return map(budgets) do b
        ts = [o for o in triouts if all(i -> us[i].kind === b.kind && us[i].role === b.role, o.units)]
        isempty(ts) && return b
        n, safe = length(ts), count(o -> o.tier === :green, ts)
        exc = [Tuple(us[i].name for i in o.units)::NTuple{3,String} for o in ts if o.tier !== :green]
        word = string(b.kind, "s")
        limit, text = if safe == n
            3, "up to 3 green $word together (all $n triples stay green; 4 or more not scanned)"
        elseif safe == 0
            2, "up to 2 green $word together (all $(b.safe) pairs stay green; each of the $n triples leaves green)"
        else
            3, "up to 3 green $word together, except $(n - safe) of $n triples (4 or more not scanned)"
        end
        ROEBudget(b.role, b.kind, b.green, b.pairs, b.safe, n, safe, limit, b.exceptions, exc, text)
    end
end
