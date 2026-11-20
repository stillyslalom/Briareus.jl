# This file retunes a drive. A retune computes per-group power histories that hold a
# drive quantity through time under power limits, for example after a bundle is dropped.
# The quantity is the spherical-harmonic shape and level of the drive on a capsule, the
# x-ray flux on a face, or the laser irradiance on a foil. It is held at the value of a
# reference shot, or uniform at a prescribed level. The file also reads and writes the
# pulse tables of the histories.

"""
    cone_fraction_pulses(beams, times, total, fraction; inner = (:inner23, :inner30))
        -> Dict{Symbol,PulseShape}

Per-group pulses (power per beam, TW) from a total laser power history `total` (TW, at
`times` ns) and a cone fraction history `fraction` (the inner-cone share of the total),
in the way NIF pulses are specified. The beams of the inner groups share
`fraction × total` equally, and the beams of the other groups share the rest. The pulses
are keyed by the `group` of each beam, for `with_pulses`.
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

Write pulses (power per beam, TW, against time, ns) per group, quad or beam as CSV. The
table has a `time_ns` column and one column per key, on the union of the times of the
pulses. [`read_pulse_table`](@ref) reads the table back.
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
[`write_pulse_table`](@ref), with a `time_ns` column and one column per key (power per
beam, TW). The measured or requested powers of a reference shot are read this way.
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

The quantity that [`retune_drive`](@ref) holds through time. It is a vector `y` of the
drive, linear in the powers of the beam groups at fixed albedos, with a level `L = μ'y`
and a shape `y/L`. The kinds are [`CapsuleModes`](@ref) (the x-ray flux on a capsule, by
its real spherical-harmonic coefficients), [`FaceFlux`](@ref) (the x-ray flux on each
element of a face, such as the face of a package in a halfraum) and
[`FoilIrradiance`](@ref) (the absorbed laser irradiance on the spot of a foil, without a
radiation model).
"""
abstract type RetuneQuantity end

"""
    CapsuleModes(model; observe = :capsule, modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0),
                 components = ())

The drive of a capsule for [`retune_drive`](@ref). The quantity is the real
spherical-harmonic coefficients `c_lm` of the incident x-ray flux on `observe` in the
radiation model `model` (`ylm_weights` about the center of the observed elements, in
`mode_index` order), with level `c_00`. The shape terms weight each degree `l => w` of
`modes` over its `2l + 1` coefficients, and each `(l, m) => w` of `components`
separately, for example the `l = 1, m = ±1` drift of a dropped bundle. The reference is
a [`drive_history`](@ref) run with `ylm = true` (through the highest degree weighted) and
`keep_albedo = true`.
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

The x-ray drive of a face for [`retune_drive`](@ref). The quantity is the incident flux
`H_e` (TW/mm²) on each element of `face` (a part or component name, or a tuple of them)
in the radiation model `model`, and the level is the area mean
`H̄ = Σ A_e H_e / Σ A_e`. The shape terms are `uniformity × Σ_e (A_e/ΣA) (H_e/H̄ - ŝ_e)²`,
times `(H̄/L_t)²` (see `retune_drive`). For the uniform target (`ŝ = 1`) this is the
square of the RMS nonuniformity of the face, `rms_nonuniformity(H, A)`. The quantity is
meant for a flat face, such as the face of a package in a halfraum or an opacity sample,
whose flux has no modes to restore. The reference is a [`drive_history`](@ref) run with
`observe = face`, `observed_flux = true` and `keep_albedo = true`.
"""
struct FaceFlux{R<:RadiationModel} <: RetuneQuantity
    model::R
    face::Any
    uniformity::Float64
end
FaceFlux(model::RadiationModel, face; uniformity::Real = 1.0) = FaceFlux(model, face, Float64(uniformity))

"""
    FoilIrradiance(foil; center = nothing, radius = Inf, side = 1, uniformity = 1.0)

The laser drive of a foil, such as a backlighter, for [`retune_drive`](@ref). The
quantity is the absorbed laser irradiance `I_e = P_e/A_e` (TW/mm²) of each element of
the part or component `foil` on side `side` (1 front, -1 back) whose centroid lies
within `radius` (mm) of `center`. By default `center` is the area centroid of the foil
elements on that side. The level is the mean `Σ P_e / Σ A_e`, and the shape terms are
those of [`FaceFlux`](@ref). These are the uniformity and intensity terms of
[`irradiance_objective`](@ref). Its spill term has no counterpart here, since the
pointing is fixed. No radiation model is used, since the deposition is linear in the
powers. The reference is `(; time, deposition)`, the absorbed laser `Deposition` of the
reference shot at each time, with every color that it fires (e.g.
`[deposit(op, with_pulses(beams, pulses, τ)) for τ in times]`). With
`target = :uniform` and a `level`, the reference is only `(; time)`.
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

# `y = X[es, :] ./ A` is the irradiance of a foil spot from per-element absorbed powers.
# It is the response of `FoilIrradiance`, to the 3ω light and to each unconverted color
# alike.
struct _SpotRows
    es::Vector{Int}
    A::Vector{Float64}
end
Base.:*(s::_SpotRows, X::AbstractMatrix) = X[s.es, :] ./ s.A
Base.:*(s::_SpotRows, x::AbstractVector) = x[s.es] ./ s.A

# Returns, per quantity, the target, the length `ny` of `y`, the level functional `μ`, the
# shape rows (index, weight) and the uniform shape. It also returns the `y` of the
# reference (`ny × nt`, or nothing when the target does not need it), the area weights of
# the RMS (nothing when there are none), and `response(τ, k, colors)`. That function
# returns a Dict of the responses from per-element absorbed power to `y` at reference
# step `k`, under `:omega3` and the name of each unconverted color.
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

# Returns the x-ray responses of `obs` (weighted by `W`, or per element) at reference step
# `k`, about its converged albedos. With unconverted colors there is also one response per
# color. The conversion of each color is independent of intensity, as the 3ω response
# requires, so the powers of the secant do not enter.
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

Power histories for `beams` that hold a drive `quantity` through the times of
`reference`, for example after a bundle is dropped. Beams are grouped `by` `:group`,
`:quad`, `:beam` or a function `beam -> key`, with one power per beam of a group at each
time. The quantity is the modes of a capsule ([`CapsuleModes`](@ref)), the x-ray flux on
a face ([`FaceFlux`](@ref)) or the laser irradiance on a foil
([`FoilIrradiance`](@ref)); their docstrings list what `reference` must contain. `op`
is the deposition operator of `beams`, for example that of the reference restricted to
the remaining beams, `op[keep]`. Its pointing is fixed through the pulse.

**Target.** The vector `y` of the quantity (the `c_lm` of a capsule, the `H_e` of a
face, the `I_e` of a foil) has a level `L = μ'y` (`c_00`, or the area mean) and a shape
`y/L`. At each time the target is a shape `ŝ` and a level `L_t`. The shape is set by
`target`:

- `target = :reference` (the default): the full shape of the reference,
  `ŝ = y_ref/L_ref`. This is the modes of a capsule, or the flux distribution of a face
  or a foil, and so also its RMS nonuniformity and peak-to-mean ratio;
- `target = :uniform`: a uniform drive (`c_lm = 0` for `l ≥ 1`, or a flat flux), which
  is the least nonuniformity that the powers allow.

`L_t` is `level` (TW/mm², a number or a function of time in ns) when given, and
otherwise the level `L_ref` of the reference. `target = :uniform, level = I₀` thus gives
the most uniform irradiance of a foil at the mean `I₀`.

**Solve.** At each reference time `y` is linear in the powers `P` of the groups,
`y = M P + y₀`, where `y₀` is the part of the prepulse (below). For the x-ray quantities
the albedos are fixed at the converged albedos of the reference
(`flux_response(...; albedo)`). The powers minimize

    Σ_i w_i ((y_i - ŝ_i L)/L_t)² + smooth Σ_g ((P_g - P_g,prev)/P_scale)²
        + stay Σ_g ((P_g - P_g,nominal)/P_scale)²

at the level `L = L_t`, with the weights `w_i` of the quantity. For a capsule they are
its modes and components. For a face or a foil they are `uniformity × A_e/ΣA`, and the
first sum is then the square of the RMS departure from `ŝ`. The powers lie within the
bounds `limits` (TW per beam), which are `nothing` (from 0, unbounded above), a `Dict`
key => `(lo, hi)` (with `(0, Inf)` for missing keys), or a function
`(key, t) -> (lo, hi)`. Examples are the power capability of a quad, or `(0, 0)` for a
dropped quad. `smooth` ties the powers of each step to those of the previous step, and
`stay` ties them to `nominal` powers, such as those of the reference shot. `nominal` is
a `Dict` key => `PulseShape` or a function `(key, t) -> P`. Among the many power sets
that restore the drive, `stay` selects the one closest to the nominal powers, instead of
one that runs many quads to their limits. The problem is a box-bounded linear
least-squares problem per time. The level is held by a penalty row whose target is
shifted until the level of the solution is `L_t`, by the secant method on that level.
The level is piecewise linear in the shift, so one step suffices while the active bounds
are unchanged. A penalty without the shift would pull the level down by about the square
of the nonuniformity (`L/L_t ≈ 1/(1 + RMS²)`), a percent at an RMS of 10%. Where the
bounds cannot reach `L_t`, the level is the nearest they allow, and `drive` reports it.

**Unconverted light.** For a reference run with unconverted (1ω, 2ω) light, the
`unconverted`, `prepulse`, `colors` and `unconverted_kw` of that history are passed here
(see [`drive_history`](@ref)), so that the retune uses the same per-color sources:

    y = (R₃ E + Σ_c R_c E_c(t)) P + Σ_c R_c Q_c p_c(t).

`R₃` and `R_c` are the 3ω and per-color responses, from
`flux_response(...; unconverted)` for the x-ray quantities; for a foil, the absorbed
power of each color adds to the irradiance. `E` is the 3ω absorbed power of each group
at 1 TW per beam. `E_c(t)` is its absorbed power in color `c` at 1 TW per beam, times the
fraction `r_c(t)` of each beam. It is the sum, over the beams of the group on each
fraction schedule, of [`unconverted_deposition`](@ref) of those beams at that fraction.
`p_c(t)` is the prepulse of the color, fixed powers deposited by `Q_c`. The fractions and
prepulses are keyed by `by`, as the retuned groups are. The fraction of a color may also
be a vector with one entry per beam of `beams`, so the beams of a retuned group may have
different fractions (e.g. 1ω only on some beams of a quad). The prepulse is not
retuned, and it shifts the target of the least squares. As for the 3ω response, the
conversions must not depend on intensity. A reference history with unconverted light
and no `unconverted` here is an error.

Returns `(; time, keys, power, pulses, shape, target_shape, level, drive, rms, peak,
residual)`. These are the per-beam powers (`length(keys) × nt`), the same powers as
`PulseShape`s for `with_pulses` (with `by` the same key), the predicted shape `y/L` and
the target `ŝ` (`ny × nt`), and the predicted level `L` (TW/mm²) and `L/L_t`. For a face
or a foil they include the predicted RMS nonuniformity and peak-to-mean ratio
`max(y)/L` (`nothing` for a capsule), and finally the weighted residual norm per time. A
time whose limits fix every power (`lo == hi`) is not solved. An example is `(0, 0)`
where the reference has no 3ω light, as with only a prepulse or in a gap between
pulses. The powers at such a time are the limits, the target level may be zero, and the
residual is `NaN`. For the x-ray quantities the prediction is linear about the albedos
of the reference, and a `drive_history` of the pulses, where the albedos follow the new
drive, confirms it.
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
    # the unconverted light, as in the reference history (`nothing` without it)
    ul = unconverted === nothing ? nothing : _UnconvertedLight(unconverted, prepulse, by, colors, unconverted_kw)
    if unconverted === nothing
        (prepulse === nothing && isempty(unconverted_kw)) ||
            throw(ArgumentError("prepulse and unconverted_kw set the unconverted light; they need `unconverted`"))
        get(ref, :unconverted, nothing) === nothing ||
            throw(ArgumentError("the reference history has unconverted light: pass its `unconverted` " *
                                "(and `prepulse`, `colors`, `unconverted_kw`)"))
    end
    # the drive of the reference, used for its shape, or for its level when none is prescribed
    need_ref = target === :reference || level === nothing
    st = q isa FoilIrradiance ? _retune_setup(q, ref, need_ref, op.target) : _retune_setup(q, ref, need_ref)
    t = st.target
    ny, μ = st.ny, st.μ
    key = _beam_key(by)
    ks = unique(key.(beams))
    ng = length(ks)
    # deposition of each group at 1 TW per beam (the pointing is fixed through the pulse)
    E = reduce(hcat, [deposit(op, [remake(b; power = key(b) == k ? 1.0 : 0.0) for b in beams]).element
                      for k in ks])
    # With unconverted light, the absorbed power of each color per group at 1 TW per beam
    # (`Eu`), and per set of beams on one prepulse schedule at 1 TW per beam (`Q`, with
    # the schedules `qs`). The setup powers of the beams carry the light, as in
    # `drive_history`.
    ne = nelements(t)
    ucolors = ul === nothing ? UnconvertedColor[] : ul.colors
    unit_light(bs) = unconverted_deposition(t, [remake(b; power = 1.0) for b in bs];
                                            ratio = Dict(c.name => 1.0 for c in ul.colors),
                                            colors = ul.colors, ul.kw...)
    Eu = Dict{Symbol,Matrix{Float64}}()
    Q = Dict{Symbol,Matrix{Float64}}()
    qs = Any[]
    # The beams of the groups on each fraction schedule are in `subs`, as (group index,
    # fraction per color, beam indices), with one entry per group when the fractions
    # apply to every beam. The absorbed power per color of each entry is in the columns
    # of `Eu`.
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
        # y at powers P is y = M P, plus y0 from the prepulse
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
        bounds = [lim(kk, τ) for kk in ks]
        lo, hi = first.(bounds), last.(bounds)
        all(lo .<= hi) || throw(ArgumentError("empty power limits at $τ ns"))
        Lt = level === nothing ? dot(μ, view(st.yref, :, k)) : level isa Real ? Float64(level) : Float64(level(τ))
        if lo == hi
            # The limits fix every power, for example at zero where the reference fires
            # no light. There is nothing to solve, and the target may be empty.
            y = y0 === nothing ? M * lo : M * lo .+ y0
            power[:, k] = lo
            target_shape[:, k] .= target === :reference ? st.yref[:, k] ./ dot(μ, view(st.yref, :, k)) : st.uniform
            lev[k] = L = dot(μ, y)
            shape[:, k] = y ./ L
            drive[k] = L / Lt
            if st.area !== nothing
                rms[k] = L > 0 ? rms_nonuniformity(y, st.area) : NaN
                peak[k] = maximum(y) / L
            end
            residual[k] = NaN
            continue
        end
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
                A = A[1:r, :]                                    # no smoothing rows without a previous step
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

# Returns argmin |A P - b|² for lo ≤ P ≤ hi at the level `LM'P + L0 = Lt` where the bounds
# allow. Row `r` of `A, b` is the penalty of the level, `(LM'P + L0 - T)/Lt`. The secant
# method shifts its target `T` until the level of the solution is `Lt`. The level is
# piecewise linear and nondecreasing in `T`, so a step is exact while the active bounds
# are unchanged. The iteration stops where the bounds fix the level. `b[r]` is left at the
# unshifted target, for the residual.
function _level_lsq(A, b, r, Lt, L0, LM, lo, hi, x0)
    T, Tp, Lp = Lt, NaN, NaN
    P = x0
    for _ in 1:30
        b[r] = (T - L0) / Lt
        P = _bounded_lsq(A, b, lo, hi, P)
        L = dot(LM, P) + L0
        abs(L - Lt) <= 1e-12 * Lt && break
        if isnan(Tp)
            Tn = T + (Lt - L)                       # the first shift is the deficit
        else
            # The level moves by at most the shift. A flatter secant comes from the bounds,
            # and a slope under 1e-6 means that they fix the level. Such a secant is kept
            # from extrapolating far past the shift.
            s = (L - Lp) / (T - Tp)
            s < 1e-6 && break
            Tn = T + (Lt - L) / clamp(s, 0.05, 1.0)
        end
        Tp, Lp, T = T, L, Tn
    end
    b[r] = (Lt - L0) / Lt
    return P
end

# Returns argmin |A x - b|² for lo ≤ x ≤ hi (small, dense) by L-BFGS-B on the quadratic.
function _bounded_lsq(A, b, lo, hi, x0)
    AtA, Atb = A' * A, A' * b
    f(x) = dot(x, AtA * x) - 2dot(Atb, x)
    g!(G, x) = (G .= 2 .* (AtA * x .- Atb))
    any(lo .== hi) && all(lo .== hi) && return copy(lo)
    res = Optim.optimize(f, g!, lo, hi, clamp.(x0, lo, hi), Optim.LBFGSB(),
                         Optim.Options(; g_abstol = 1e-14, iterations = 2000))
    return _polish_lsq(A, b, lo, hi, Optim.minimizer(res))
end

# Takes active-set Gauss–Newton steps from a bounded least-squares point `x`. Each step
# is the least-norm correction of the residual on the free variables, which are those off
# their bounds or at a bound with the gradient pointing inward. The step is cut where it
# meets the box, and the variable that meets the box is held there at the next step.
# Steps are taken while they lower the residual. On the rank-deficient problems of the
# retune, L-BFGS-B stops at a projected gradient near 1e-6 of the gradient scale, and the
# steps reduce it to rounding.
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
            # variables at a bound from which the step would leave the box are held
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
