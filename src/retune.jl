# Retuning a drive: per-group power histories that hold a drive quantity through time
# (a capsule's spherical-harmonic shape and level, a face's x-ray flux, a foil's laser
# irradiance), to a reference shot's or uniform at a prescribed level, under power
# limits, e.g. after dropping a bundle. Plus the pulse tables they come from.

"""
    cone_fraction_pulses(beams, times, total, fraction; inner = (:inner23, :inner30))
        -> Dict{Symbol,PulseShape}

Per-group pulses (power per beam, TW) from a total laser power history `total` (TW, at
`times` ns) and a cone fraction history `fraction` (inner-cone share of the total), the
way NIF pulses are specified: the inner groups' beams share `fraction × total` equally,
the other groups' beams the rest. Keyed by each beam's `group`, for `with_pulses`.
"""
function cone_fraction_pulses(beams::AbstractVector{<:Beam}, times, total, fraction;
                              inner = (:inner23, :inner30))
    length(times) == length(total) == length(fraction) ||
        throw(ArgumentError("times, total and fraction differ in length"))
    all(f -> 0 <= f <= 1, fraction) || throw(ArgumentError("the cone fraction must lie in [0, 1]"))
    nin = count(b -> b.group in inner, beams)
    nout = length(beams) - nin
    nin > 0 && nout > 0 || throw(ArgumentError("the beams need both inner ($inner) and outer groups"))
    out = Dict{Symbol,PulseShape}()
    for g in unique(b.group for b in beams)
        out[g] = g in inner ? PulseShape(times, total .* fraction ./ nin) :
                 PulseShape(times, total .* (1 .- fraction) ./ nout)
    end
    return out
end

"""
    write_pulse_table(path, pulses) -> path

Pulses (power per beam, TW, against time, ns) per group, quad or beam as CSV: a
`time_ns` column and one column per key, on the union of the pulses' times. Read back
with [`read_pulse_table`](@ref).
"""
function write_pulse_table(path::AbstractString, pulses::AbstractDict)
    ks = sort(collect(keys(pulses)); by = string)
    ts = sort(unique(reduce(vcat, [pulses[k].t for k in ks])))
    open(path, "w") do io
        println(io, join(["time_ns"; string.(ks)], ","))
        for τ in ts
            println(io, join(_csvfield.([τ; [pulses[k](τ) for k in ks]]), ","))
        end
    end
    return path
end

"""
    read_pulse_table(path) -> Dict{Symbol,PulseShape}

Pulses per group, quad or beam from a CSV table in the layout of
[`write_pulse_table`](@ref): a `time_ns` column and one column per key (power per beam,
TW). Read a reference shot's measured or requested powers this way.
"""
function read_pulse_table(path::AbstractString)
    header, rows = read_csv(path)
    ti = findfirst(==("time_ns"), strip.(header))
    ti === nothing && throw(ArgumentError("$(basename(path)) has no time_ns column"))
    t = [parse(Float64, r[ti]) for r in rows]
    return Dict{Symbol,PulseShape}(Symbol(strip(h)) => PulseShape(t, [parse(Float64, r[j]) for r in rows])
                                   for (j, h) in enumerate(header) if j != ti)
end

"""
    RetuneQuantity

What [`retune_drive`](@ref) holds through time: a vector `y` of the drive, linear in the
beam groups' powers at fixed albedos, with a level `L = μ'y` and a shape `y/L`. The kinds
are [`CapsuleModes`](@ref) (a capsule's x-ray flux, by its real spherical-harmonic
coefficients), [`FaceFlux`](@ref) (the x-ray flux on each element of a face, e.g. a
package's face in a halfraum) and [`FoilIrradiance`](@ref) (the absorbed laser irradiance
on a foil's spot, without a radiation model).
"""
abstract type RetuneQuantity end

"""
    CapsuleModes(model; observe = :capsule, modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0),
                 components = ())

A capsule's drive for [`retune_drive`](@ref): the real spherical-harmonic coefficients
`c_lm` of the incident x-ray flux on `observe` in the radiation model `model`
(`ylm_weights` about the observed elements' center, `mode_index` order), with level
`c_00`. The shape terms weight each degree `l => w` of `modes` over its `2l + 1`
coefficients and each `(l, m) => w` of `components` alone (e.g. the `l = 1, m = ±1`
drift of a dropped bundle). The reference is a [`drive_history`](@ref) run with
`ylm = true` (through the highest degree weighted) and `keep_albedo = true`.
"""
struct CapsuleModes{R<:RadiationModel} <: RetuneQuantity
    model::R
    observe::Any
    modes::Vector{Pair{Int,Float64}}
    components::Vector{Pair{Tuple{Int,Int},Float64}}
end
CapsuleModes(model::RadiationModel; observe = :capsule, modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0),
             components = ()) =
    CapsuleModes(model, observe, Pair{Int,Float64}[Int(l) => Float64(w) for (l, w) in modes],
                 Pair{Tuple{Int,Int},Float64}[(Int(lm[1]), Int(lm[2])) => Float64(w) for (lm, w) in components])

"""
    FaceFlux(model, face; uniformity = 1.0)

A face's x-ray drive for [`retune_drive`](@ref): the incident flux `H_e` (TW/mm²) on
each element of `face` (a part or component name, or a tuple of them) in the radiation
model `model`, with level the area mean `H̄ = Σ A_e H_e / Σ A_e`. The shape terms are
`uniformity × Σ_e (A_e/ΣA) (H_e/H̄ - ŝ_e)²` (times `(H̄/L_t)²`, see `retune_drive`): for the
uniform target (`ŝ = 1`) the square of the face's RMS nonuniformity,
`rms_nonuniformity(H, A)`. Meant for a flat face (a package's face in a halfraum, an
opacity sample), whose flux has no modes to restore. The reference is a
[`drive_history`](@ref) run with `observe = face`, `observed_flux = true` and
`keep_albedo = true`.
"""
struct FaceFlux{R<:RadiationModel} <: RetuneQuantity
    model::R
    face::Any
    uniformity::Float64
end
FaceFlux(model::RadiationModel, face; uniformity::Real = 1.0) = FaceFlux(model, face, Float64(uniformity))

"""
    FoilIrradiance(foil; center = nothing, radius = Inf, side = 1, uniformity = 1.0)

A foil's laser drive for [`retune_drive`](@ref) (a backlighter): the absorbed laser
irradiance `I_e = P_e/A_e` (TW/mm²) of each element of the part or component `foil` on
side `side` (1 front, -1 back) whose centroid lies within `radius` (mm) of `center`
(default: the area centroid of the foil's elements on that side), with level the mean
`Σ P_e / Σ A_e` and shape terms as [`FaceFlux`](@ref)'s. These are the uniformity and
intensity terms of [`irradiance_objective`](@ref) (its spill has no counterpart here:
the pointing is fixed). No radiation model takes part: the deposition is linear in the
powers. The reference is `(; time, deposition)`, the reference shot's absorbed laser
`Deposition` at each time (e.g. `[deposit(op, with_pulses(beams, pulses, τ)) for τ in
times]`, every color it fires included); with `target = :uniform` and a `level`,
`(; time)` alone.
"""
struct FoilIrradiance <: RetuneQuantity
    foil::Symbol
    center::Union{Nothing,Vec3}
    radius::Float64
    side::Int
    uniformity::Float64
end
FoilIrradiance(foil::Symbol; center = nothing, radius::Real = Inf, side::Integer = 1, uniformity::Real = 1.0) =
    FoilIrradiance(foil, center === nothing ? nothing : Vec3(center...), Float64(radius), Int(side),
                   Float64(uniformity))

# `y = X[es, :] ./ A`: a foil spot's irradiance from per-element absorbed powers, the
# response of `FoilIrradiance` (to the 3ω light and to each unconverted color alike).
struct _SpotRows
    es::Vector{Int}
    A::Vector{Float64}
end
Base.:*(s::_SpotRows, X::AbstractMatrix) = X[s.es, :] ./ s.A
Base.:*(s::_SpotRows, x::AbstractVector) = x[s.es] ./ s.A

# Per quantity: the target, the length `ny` of `y`, the level functional `μ`, the shape
# rows (index, weight), the uniform shape, the reference's `y` (`ny × nt`, or nothing
# when the target does not need it), the area weights of the RMS (nothing: none), and
# `response(τ, k, colors)`, a Dict of the responses (`:omega3` and each unconverted
# color's) from per-element absorbed power to `y` at reference step `k`.
function _retune_setup(q::CapsuleModes, ref, need_ref::Bool)
    m = q.model
    t = m.target
    obs = elements(t, (q.observe isa Symbol ? (q.observe,) : q.observe)...)
    lmax = maximum(vcat([l for (l, _) in q.modes], [lm[1] for (lm, _) in q.components], 0))
    nm = (lmax + 1)^2
    _need_albedo(ref)
    yref = nothing
    if need_ref
        get(ref, :ylm, nothing) === nothing && throw(ArgumentError("the reference history needs ylm = true"))
        size(ref.ylm, 1) >= nm || throw(ArgumentError("the reference's Y_lm reach l = $(isqrt(size(ref.ylm, 1)) - 1)"))
        yref = ref.ylm[1:nm, :]
    end
    Y = ylm_weights(t, obs; lmax)
    rows = Tuple{Int,Float64}[]
    for (l, w) in q.modes, mm in -l:l
        push!(rows, (mode_index(l, mm), w))
    end
    for ((l, mm), w) in q.components
        push!(rows, (mode_index(l, mm), w))
    end
    e1 = zeros(nm)
    e1[1] = 1.0
    return (; target = t, ny = nm, μ = e1, rows, uniform = e1, yref, area = nothing,
            response = _xray_response(m, obs, ref, Y))
end

function _retune_setup(q::FaceFlux, ref, need_ref::Bool)
    m = q.model
    t = m.target
    obs = elements(t, (q.face isa Symbol ? (q.face,) : q.face)...)
    isempty(obs) && throw(ArgumentError("the face $(repr(q.face)) has no elements"))
    q.uniformity >= 0 || throw(ArgumentError("uniformity must not be negative"))
    _need_albedo(ref)
    yref = nothing
    if need_ref
        get(ref, :observed_flux, nothing) === nothing &&
            throw(ArgumentError("the reference history needs observed_flux = true (observing the face)"))
        size(ref.observed_flux, 1) == length(obs) ||
            throw(ArgumentError("the reference history observed $(size(ref.observed_flux, 1)) elements, " *
                                "the face has $(length(obs))"))
        yref = ref.observed_flux
    end
    A = t.elem_area[obs]
    μ = A ./ sum(A)
    return (; target = t, ny = length(obs), μ, rows = [(i, q.uniformity * μ[i]) for i in eachindex(obs)],
            uniform = ones(length(obs)), yref, area = A, response = _xray_response(m, obs, ref, nothing))
end

function _retune_setup(q::FoilIrradiance, ref, need_ref::Bool, t::Target)
    q.uniformity >= 0 || throw(ArgumentError("uniformity must not be negative"))
    q.radius > 0 || throw(ArgumentError("the spot radius must be positive"))
    on = [e for e in elements(t, q.foil) if t.elem_side[e] == q.side]
    isempty(on) && throw(ArgumentError("the foil $(q.foil) has no elements on side $(q.side)"))
    Aon = t.elem_area[on]
    c0 = q.center === nothing ? sum(Aon .* t.elem_centroid[on]) / sum(Aon) : q.center
    es = [e for e in on if norm(t.elem_centroid[e] - c0) <= q.radius]
    isempty(es) && throw(ArgumentError("no elements of $(q.foil) within $(q.radius) mm of the spot center"))
    A = t.elem_area[es]
    S = _SpotRows(es, A)
    yref = nothing
    if need_ref
        get(ref, :deposition, nothing) === nothing &&
            throw(ArgumentError("the foil's reference needs `deposition`, the reference shot's at each time"))
        length(ref.deposition) == length(ref.time) ||
            throw(ArgumentError("the reference has $(length(ref.deposition)) depositions for $(length(ref.time)) times"))
        yref = reduce(hcat, [S * d.element for d in ref.deposition])
    end
    μ = A ./ sum(A)
    response(τ, k, colors) = Dict{Symbol,Any}(:omega3 => S, (c.name => S for c in colors)...)
    return (; target = t, ny = length(es), μ, rows = [(i, q.uniformity * μ[i]) for i in eachindex(es)],
            uniform = ones(length(es)), yref, area = A, response)
end

_need_albedo(ref) = get(ref, :albedo, nothing) === nothing &&
    throw(ArgumentError("the reference history needs keep_albedo = true"))

# The x-ray responses of `obs` (weighted by `W`, or per element) at reference step `k`,
# about its converged albedos; with unconverted colors, one per color as well, each
# color's conversion independent of intensity (as the 3ω response requires), so the
# secant's powers do not enter.
function _xray_response(m, obs, ref, W)
    ne = nelements(m.target)
    kw(k) = W === nothing ? (; albedo = ref.albedo[k]) : (; albedo = ref.albedo[k], weights = W)
    return function (τ, k, colors)
        isempty(colors) && return Dict{Symbol,Any}(:omega3 => flux_response(m, obs; time = τ, kw(k)...))
        R = flux_response(m, obs; time = τ, kw(k)...,
                          unconverted = _with_colors(Dict(c.name => zeros(ne) for c in colors), colors))
        return Dict{Symbol,Any}(R)
    end
end

"""
    retune_drive(quantity, op, beams, reference; target = :reference, level = nothing,
                 by = :group, limits = nothing, smooth = 0.0, nominal = nothing, stay = 0.0,
                 unconverted = nothing, prepulse = nothing, colors = nif_unconverted(),
                 unconverted_kw = (;)) -> NamedTuple

Power histories for `beams` (grouped `by` `:group`, `:quad`, `:beam` or a function
`beam -> key`; one power per beam of a group at each time) that hold a drive `quantity`
through the times of `reference`, e.g. after dropping a bundle: a capsule's modes
([`CapsuleModes`](@ref)), a face's x-ray flux ([`FaceFlux`](@ref)) or a foil's laser
irradiance ([`FoilIrradiance`](@ref)), whose docstrings list what `reference` must hold.
`op` is the deposition operator of `beams` (e.g. the reference's restricted to the beams
that remain, `op[keep]`), whose pointing holds through the pulse.

**Target.** The quantity's vector `y` (a capsule's `c_lm`, a face's `H_e`, a foil's `I_e`)
has a level `L = μ'y` (`c_00`; the area mean) and a shape `y/L`. At each time the
target is a shape `ŝ` and a level `L_t`:

- `target = :reference` (the default): the reference's shape, `ŝ = y_ref/L_ref`, all of
  it: a capsule's modes, a face's or a foil's flux distribution (so its RMS
  nonuniformity and peak-to-mean ratio as well);
- `target = :uniform`: a uniform drive (`c_lm = 0` for `l ≥ 1`; a flat flux), i.e. the
  least nonuniformity that the powers allow;

and `L_t` is `level` (TW/mm²; a number or a function of time in ns) when given, else the
reference's `L_ref`. `target = :uniform, level = I₀` thus gives a foil's most uniform
irradiance at the mean `I₀`.

**Solve.** At each reference time `y` is linear in the groups' powers `P`: `y = M P + y₀`,
for the x-ray quantities with the albedos fixed at the reference's converged ones
(`flux_response(...; albedo)`), `y₀` the prepulse's part (below). The powers minimize

    Σ_i w_i ((y_i - ŝ_i L)/L_t)² + smooth Σ_g ((P_g - P_g,prev)/P_scale)²
        + stay Σ_g ((P_g - P_g,nominal)/P_scale)²

at the level `L = L_t`, with the quantity's weights `w_i` (its modes and components for
a capsule; `uniformity × A_e/ΣA` for a face or a foil, whose first sum is then the
square of its RMS departure from `ŝ`), within the bounds `limits` (TW per beam):
`nothing` (from 0, unbounded above), a `Dict` key => `(lo, hi)` (missing keys:
`(0, Inf)`), or a function `(key, t) -> (lo, hi)`, e.g. a quad's power capability, or
`(0, 0)` for a dropped quad. `smooth` ties each step's powers to the previous step's,
and `stay` to `nominal` powers (a `Dict` key => `PulseShape`, or a function
`(key, t) -> P`, e.g. the reference shot's): among the many power sets that restore the
drive it selects the one closest to them, instead of one that runs many quads to their
limits. This is a box-bounded linear least-squares problem per time, the level held
by a penalty row whose target is shifted until the solution's level is `L_t` (the
secant method on that level, piecewise linear in the shift: one step while the active
bounds hold). A penalty alone would pull the level down by about the square of the
nonuniformity (`L/L_t ≈ 1/(1 + RMS²)`), a percent at an RMS of 10 %. Where the bounds
cannot reach `L_t` the level is the nearest they allow, and `drive` reports it.

**Unconverted light.** A reference run with unconverted (1ω, 2ω) light takes that
history's `unconverted`, `prepulse`, `colors` and `unconverted_kw` here (see
[`drive_history`](@ref)), so that the retune solves with the same per-color sources:
`y = (R₃ E + Σ_c R_c E_c(t)) P + Σ_c R_c Q_c p_c(t)`, with `R₃` and `R_c` the 3ω and
per-color responses (`flux_response(...; unconverted)` for the x-ray quantities; for a
foil, each color's absorbed power adds to the irradiance), `E` each group's 3ω absorbed
power at 1 TW per beam, `E_c(t)` its color-`c` absorbed power at 1 TW per beam times each
beam's fraction `r_c(t)` (the sum, over its beams on each fraction schedule, of
[`unconverted_deposition`](@ref) of those beams at that fraction), and `p_c(t)` the
color's prepulse, fixed powers deposited by `Q_c`. The fractions and prepulses are keyed
by `by`, as the retuned groups, and a color's fraction may also be a vector with one
entry per beam of `beams`: the beams of a retuned group may carry different fractions
(e.g. 1ω on some beams of a quad only). The prepulse is not retuned: it shifts the
target of the least squares. As for the 3ω response, the conversions must not depend on
intensity. A reference history with unconverted light and no `unconverted` here is an
error.

Returns `(; time, keys, power, pulses, shape, target_shape, level, drive, rms, peak,
residual)`: the per-beam powers (`length(keys) × nt`), the same as `PulseShape`s for
`with_pulses` (`by` the same key), the predicted shape `y/L` and the target `ŝ`
(`ny × nt`), the predicted level `L` (TW/mm²) and `L/L_t`, for a face or a foil the
predicted RMS nonuniformity and peak-to-mean ratio `max(y)/L` (`nothing` for a capsule),
and the weighted residual norm per time. For the x-ray quantities the prediction is
linear about the reference's albedos: confirm with a `drive_history` of the pulses,
where the albedos follow the new drive.
"""
Base.@nospecializeinfer function retune_drive(@nospecialize(q::RetuneQuantity), op::DepositionOperator,
                      beams::AbstractVector{<:Beam}, @nospecialize(ref); target::Symbol = :reference,
                      @nospecialize(level = nothing), @nospecialize(by = :group),
                      @nospecialize(limits = nothing), smooth = 0.0,
                      @nospecialize(nominal = nothing), stay = 0.0, @nospecialize(unconverted = nothing),
                      @nospecialize(prepulse = nothing), colors = nif_unconverted(),
                      @nospecialize(unconverted_kw = (;)))
    target in (:reference, :uniform) || throw(ArgumentError("target must be :reference or :uniform, not :$target"))
    stay > 0 && nominal === nothing && throw(ArgumentError("stay needs nominal powers"))
    level isa Real && !(level > 0) && throw(ArgumentError("the level must be positive"))
    nom(k, τ) = nominal isa AbstractDict ? nominal[k](τ) : nominal(k, τ)
    # the unconverted light, as the reference history's (`nothing`: none)
    ul = unconverted === nothing ? nothing : _UnconvertedLight(unconverted, prepulse, by, colors, unconverted_kw)
    if unconverted === nothing
        (prepulse === nothing && isempty(unconverted_kw)) ||
            throw(ArgumentError("prepulse and unconverted_kw set the unconverted light; they need `unconverted`"))
        get(ref, :unconverted, nothing) === nothing ||
            throw(ArgumentError("the reference history has unconverted light: pass its `unconverted` " *
                                "(and `prepulse`, `colors`, `unconverted_kw`)"))
    end
    # the reference's own drive: for its shape, or for its level when none is prescribed
    need_ref = target === :reference || level === nothing
    st = q isa FoilIrradiance ? _retune_setup(q, ref, need_ref, op.target) : _retune_setup(q, ref, need_ref)
    t = st.target
    ny, μ = st.ny, st.μ
    key = _beam_key(by)
    ks = unique(key.(beams))
    ng = length(ks)
    # each group's deposition at 1 TW per beam: the pointing holds through the pulse
    E = reduce(hcat, [deposit(op, [remake(b; power = key(b) == k ? 1.0 : 0.0) for b in beams]).element
                      for k in ks])
    # with unconverted light, each color's absorbed power per group at 1 TW per beam
    # (`Eu`), and per set of beams on one prepulse schedule at 1 TW per beam (`Q`, with
    # the schedules `qs`); the beams' setup powers carry it, as in `drive_history`
    ne = nelements(t)
    ucolors = ul === nothing ? UnconvertedColor[] : ul.colors
    unit_light(bs) = unconverted_deposition(t, [remake(b; power = 1.0) for b in bs];
                                            ratio = Dict(c.name => 1.0 for c in ul.colors),
                                            colors = ul.colors, ul.kw...)
    Eu = Dict{Symbol,Matrix{Float64}}()
    Q = Dict{Symbol,Matrix{Float64}}()
    qs = Any[]
    # the groups' beams on each fraction schedule (`subs`: group index, fraction per
    # color, beam indices; one per group when the fractions hold for every beam), with
    # each one's absorbed power per color in the columns of `Eu`
    subs = Tuple{Int,Any,Vector{Int}}[]
    if ul !== nothing
        fr = _fractions_of(ul.ratio, Symbol[c.name for c in ul.colors], key, beams)
        for (g, k) in enumerate(ks)
            idx = [i for (i, b) in enumerate(beams) if key(b) == k]
            sg, gg = _schedule_groups(idx, i -> fr[i])
            append!(subs, [(g, f, idx[p]) for (f, p) in zip(sg, gg)])
        end
        du = [unit_light(beams[idx]) for (_, _, idx) in subs]
        for c in ul.colors
            Eu[c.name] = reduce(hcat, [d[c.name].element for d in du]; init = zeros(ne, 0))
        end
        if prepulse !== nothing
            qs, qgroups = _schedule_groups(beams, b -> _prepulse_of(ul, b))
            dq = [unit_light(beams[g]) for g in qgroups]
            for c in ul.colors
                Q[c.name] = reduce(hcat, [d[c.name].element for d in dq]; init = zeros(ne, 0))
            end
        end
    end
    lim(k, τ) = limits === nothing ? (0.0, Inf) :
                limits isa AbstractDict ? get(limits, k, (0.0, Inf)) : limits(k, τ)
    nt = length(ref.time)
    power = zeros(ng, nt)
    shape = zeros(ny, nt)
    target_shape = zeros(ny, nt)
    lev = zeros(nt)
    drive = zeros(nt)
    residual = zeros(nt)
    rms = st.area === nothing ? nothing : zeros(nt)
    peak = st.area === nothing ? nothing : zeros(nt)
    nr = length(st.rows) + 1 + (stay > 0 ? ng : 0) + (smooth > 0 ? ng : 0)
    prev = nothing
    for (k, τ) in enumerate(ref.time)
        # y at powers P: y = M P (+ y0, the prepulse's)
        R = st.response(τ, k, ucolors)
        M = R[:omega3] * E                                       # ny × ng
        y0 = nothing
        for (j, c) in enumerate(ucolors)
            X = R[c.name] * Eu[c.name]
            for (si, (g, f, _)) in enumerate(subs)
                rc = _fraction_at(f[j], τ)
                rc >= 0 || throw(ArgumentError("the unconverted fractions must not be negative (at $τ ns)"))
                @views M[:, g] .+= rc .* X[:, si]
            end
        end
        if !isempty(qs)
            y0 = zeros(ny)
            for (j, c) in enumerate(ucolors)
                w = Float64[qq[j] === nothing ? 0.0 : Float64(qq[j](τ)) for qq in qs]
                all(>=(0), w) || throw(ArgumentError("the prepulse power must not be negative (at $τ ns)"))
                y0 .+= R[c.name] * (Q[c.name] * w)
            end
        end
        Lt = level === nothing ? dot(μ, view(st.yref, :, k)) : level isa Real ? Float64(level) : Float64(level(τ))
        Lt > 0 || throw(ArgumentError("the target level must be positive (at $τ ns)"))
        s = target === :reference ? st.yref[:, k] ./ dot(μ, view(st.yref, :, k)) : st.uniform
        target_shape[:, k] = s
        LM = vec(μ' * M)                                         # the level per group
        L0 = y0 === nothing ? 0.0 : dot(μ, y0)
        A = zeros(nr, ng)
        b = zeros(nr)
        for (r, (i, w)) in enumerate(st.rows)
            sw = sqrt(w)
            @views A[r, :] .= sw .* (M[i, :] .- s[i] .* LM) ./ Lt
            y0 === nothing || (b[r] = -sw * (y0[i] - s[i] * L0) / Lt)
        end
        r = rl = length(st.rows) + 1                             # the level row
        A[r, :] .= LM ./ Lt
        b[r] = (Lt - L0) / Lt
        bounds = [lim(kk, τ) for kk in ks]
        lo, hi = first.(bounds), last.(bounds)
        all(lo .<= hi) || throw(ArgumentError("empty power limits at $τ ns"))
        if stay > 0
            Pn = [nom(kk, τ) for kk in ks]
            scale = max(sum(Pn) / ng, eps())
            for g in 1:ng
                A[r+g, g] = sqrt(stay) / scale
                b[r+g] = sqrt(stay) * Pn[g] / scale
            end
            r += ng
        end
        if smooth > 0
            if prev === nothing
                A = A[1:r, :]                                    # no previous step to tie to
                b = b[1:r]
            else
                scale = max(maximum(abs, prev), eps())
                for g in 1:ng
                    A[r+g, g] = sqrt(smooth) / scale
                    b[r+g] = sqrt(smooth) * prev[g] / scale
                end
            end
        end
        x0 = prev === nothing ? clamp.(fill(Lt / max(sum(LM), eps()), ng), lo, hi) : clamp.(prev, lo, hi)
        P = _level_lsq(A, b, rl, Lt, L0, LM, lo, hi, x0)
        power[:, k] = P
        y = y0 === nothing ? M * P : M * P .+ y0
        L = dot(μ, y)
        shape[:, k] = y ./ L
        lev[k] = L
        drive[k] = L / Lt
        if st.area !== nothing
            rms[k] = rms_nonuniformity(y, st.area)
            peak[k] = maximum(y) / L
        end
        residual[k] = norm(A * P - b)
        prev = P
    end
    pulses = Dict{Any,PulseShape}(k => PulseShape(ref.time, power[g, :]) for (g, k) in enumerate(ks))
    return (; time = collect(ref.time), keys = ks, power, pulses, shape, target_shape, level = lev, drive,
            rms, peak, residual)
end

# argmin |A P - b|² for lo ≤ P ≤ hi at the level `LM'P + L0 = Lt` where the bounds allow:
# row `r` of `A, b` is the level's penalty, `(LM'P + L0 - T)/Lt`, whose target `T` the
# secant method shifts until the solution's level is `Lt` (the level is piecewise linear
# and nondecreasing in `T`, so a step is exact while the active bounds hold). Stops where
# the bounds hold the level. `b[r]` is left at the unshifted target, for the residual.
function _level_lsq(A, b, r, Lt, L0, LM, lo, hi, x0)
    T, Tp, Lp = Lt, NaN, NaN
    P = x0
    for _ in 1:30
        b[r] = (T - L0) / Lt
        P = _bounded_lsq(A, b, lo, hi, P)
        L = dot(LM, P) + L0
        abs(L - Lt) <= 1e-12 * Lt && break
        if isnan(Tp)
            Tn = T + (Lt - L)                       # the first shift: the deficit
        else
            # the level moves by at most the shift; a flatter secant is the bounds' (a
            # slope under 1e-6: they hold the level), kept from extrapolating far past it
            s = (L - Lp) / (T - Tp)
            s < 1e-6 && break
            Tn = T + (Lt - L) / clamp(s, 0.05, 1.0)
        end
        Tp, Lp, T = T, L, Tn
    end
    b[r] = (Lt - L0) / Lt
    return P
end

# argmin |A x - b|² for lo ≤ x ≤ hi (small, dense): L-BFGS-B on the quadratic.
function _bounded_lsq(A, b, lo, hi, x0)
    AtA, Atb = A' * A, A' * b
    f(x) = dot(x, AtA * x) - 2dot(Atb, x)
    g!(G, x) = (G .= 2 .* (AtA * x .- Atb))
    any(lo .== hi) && all(lo .== hi) && return copy(lo)
    res = Optim.optimize(f, g!, lo, hi, clamp.(x0, lo, hi), Optim.LBFGSB(),
                         Optim.Options(; g_abstol = 1e-14, iterations = 2000))
    return _polish_lsq(A, b, lo, hi, Optim.minimizer(res))
end

# Active-set Gauss–Newton steps from a bounded least-squares point `x`: on the free
# variables (off their bounds, or at one with the gradient pointing inward), the
# least-norm correction of the residual, cut where it meets the box (the variable that
# meets it is held there at the next step), taken while it lowers the residual.
# L-BFGS-B stops at a projected gradient near 1e-6 of the gradient's scale on the
# retune's rank-deficient problems; the steps take it to rounding.
function _polish_lsq(A, b, lo, hi, x)
    x = clamp.(x, lo, hi)
    for _ in 1:(2length(x) + 5)
        r = b .- A * x
        G = -(A' * r)
        free = [lo[i] < x[i] < hi[i] || (x[i] <= lo[i] && G[i] < 0) || (x[i] >= hi[i] && G[i] > 0)
                for i in eachindex(x)]
        any(free) || break
        F = findall(free)
        δ, α = Float64[], 0.0
        while !isempty(F)
            δ = pinv(A[:, F]) * r
            all(isfinite, δ) || return x
            # variables at a bound that the step would leave the box from are held
            stuck = [(δ[j] > 0 && x[i] >= hi[i]) || (δ[j] < 0 && x[i] <= lo[i]) for (j, i) in enumerate(F)]
            any(stuck) || break
            F = F[.!stuck]
        end
        isempty(F) && break
        α = 1.0
        for (j, i) in enumerate(F)
            δ[j] > 0 && hi[i] < Inf && (α = min(α, (hi[i] - x[i]) / δ[j]))
            δ[j] < 0 && (α = min(α, (lo[i] - x[i]) / δ[j]))
        end
        y = copy(x)
        y[F] .= clamp.(x[F] .+ α .* δ, lo[F], hi[F])
        norm(b .- A * y) < norm(r) || break
        x = y
    end
    return x
end
