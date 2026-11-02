# Retuning a drive to a reference: per-group power histories that reproduce a reference
# shot's capsule drive (its spherical-harmonic shape and level) through time, under
# power limits, e.g. after dropping a bundle. Plus the pulse tables they come from.

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
    retune_drive(model, op, beams, reference; by = :group, limits = nothing,
                 modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0), components = (),
                 drive_weight = 1.0, smooth = 0.0, nominal = nothing, stay = 0.0,
                 observe = :capsule, unconverted = nothing, prepulse = nothing,
                 colors = nif_unconverted(), unconverted_kw = (;)) -> NamedTuple

Power histories for `beams` (grouped `by` `:group`, `:quad`, `:beam` or a function
`beam -> key`; one power per beam of a group at each time) that reproduce a reference
drive: `reference` is a `drive_history` of the reference shot run with `ylm = true` and
`keep_albedo = true`, whose capsule-flux coefficients `c_lm(t)` and albedos are the
target. `op` is the deposition operator of `beams` (e.g. the reference's restricted to
the beams that remain, `op[keep]`), whose pointing holds through the pulse.

At each reference time the capsule flux is linear in the groups' powers, with the
albedos fixed at the reference's converged ones (`flux_response(...; albedo)`). The
powers minimize

    Σ_(l => w) w Σ_m ((c_lm - ŝ_lm c_00)/r_00)² + Σ_((l,m) => w) w ((c_lm - ŝ_lm c_00)/r_00)²
        + drive_weight ((c_00 - r_00)/r_00)² + smooth Σ_g ((P_g - P_g,prev)/P_scale)²
        + stay Σ_g ((P_g - P_g,nominal)/P_scale)²,

where `ŝ = r/r_00` is the reference shape, within the bounds `limits` (TW per beam):
`nothing` (from 0, unbounded above), a `Dict` key => `(lo, hi)` (missing keys: `(0, Inf)`),
or a function `(key, t) -> (lo, hi)`, e.g. a quad's power capability, or `(0, 0)` for a
dropped quad. `smooth` ties each step's powers to the previous step's, and `stay` to
`nominal` powers (a `Dict` key => `PulseShape`, or a function `(key, t) -> P`, e.g. the
reference shot's): among the many power sets that restore the drive it selects the one
closest to them, instead of one that runs many quads to their limits. This is a
box-bounded linear least-squares problem per time; the drive term keeps the level, so
`drive_weight` must be positive.

**Unconverted light.** A reference run with unconverted (1ω, 2ω) light takes that
history's `unconverted`, `prepulse`, `colors` and `unconverted_kw` here (see
[`drive_history`](@ref)), so that the retune solves with the same per-color sources:
the flux is then `c = (R₃ E + Σ_c r_c(t) R_c E_c) P + Σ_c R_c Q_c p_c(t)`, with `R₃` and
`R_c` the 3ω and per-color responses (`flux_response(...; unconverted)`), `E` and `E_c`
each group's 3ω and color-`c` absorbed power at 1 TW per beam
([`unconverted_deposition`](@ref) of the group's beams), `r_c(t)` the color's fraction
and `p_c(t)` its prepulse, fixed powers deposited by `Q_c` (keyed by `by`, as the
retuned groups). The prepulse is not retuned: it shifts the target of the least squares.
As for the 3ω response, the conversions must not depend on intensity. A reference with
unconverted light and no `unconverted` here is an error.

Returns `(; time, keys, power, pulses, shape, reference_shape, drive, residual)`: the
per-beam powers (`length(keys) × nt`), the same as `PulseShape`s for `with_pulses`
(`by` the same key), the predicted and reference `c_lm/c_00` (`(lmax+1)² × nt`), the
predicted `c_00/r_00` and the weighted residual norm per time. The prediction is
linear about the reference's albedos: confirm with a `drive_history` of the pulses,
where the albedos follow the new drive.
"""
Base.@nospecializeinfer function retune_drive(m::RadiationModel, op::DepositionOperator,
                      beams::AbstractVector{<:Beam}, ref; @nospecialize(by = :group),
                      @nospecialize(limits = nothing),
                      modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0), components = (),
                      drive_weight = 1.0, smooth = 0.0, @nospecialize(nominal = nothing), stay = 0.0,
                      observe = :capsule, @nospecialize(unconverted = nothing),
                      @nospecialize(prepulse = nothing), colors = nif_unconverted(),
                      @nospecialize(unconverted_kw = (;)))
    drive_weight > 0 || throw(ArgumentError("drive_weight must be positive (it keeps the drive level)"))
    stay > 0 && nominal === nothing && throw(ArgumentError("stay needs nominal powers"))
    nom(k, τ) = nominal isa AbstractDict ? nominal[k](τ) : nominal(k, τ)
    ref.ylm === nothing && throw(ArgumentError("the reference history needs ylm = true"))
    get(ref, :albedo, nothing) === nothing &&
        throw(ArgumentError("the reference history needs keep_albedo = true"))
    # the unconverted light, as the reference history's (`nothing`: none)
    ul = unconverted === nothing ? nothing : _UnconvertedLight(unconverted, prepulse, by, colors, unconverted_kw)
    if unconverted === nothing
        (prepulse === nothing && isempty(unconverted_kw)) ||
            throw(ArgumentError("prepulse and unconverted_kw set the unconverted light; they need `unconverted`"))
        get(ref, :unconverted, nothing) === nothing ||
            throw(ArgumentError("the reference history has unconverted light: pass its `unconverted` " *
                                "(and `prepulse`, `colors`, `unconverted_kw`)"))
    end
    t = m.target
    obs = elements(t, (observe isa Symbol ? (observe,) : observe)...)
    lmax = maximum(vcat([l for (l, _) in modes], [lm[1] for (lm, _) in components], 0))
    nm = (lmax + 1)^2
    size(ref.ylm, 1) >= nm || throw(ArgumentError("the reference's Y_lm reach l = $(isqrt(size(ref.ylm, 1)) - 1)"))
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
    unit_light(bs) = unconverted_deposition(t, [remake(b; power = 1.0) for b in bs];
                                            ratio = Dict(c.name => 1.0 for c in ul.colors),
                                            colors = ul.colors, ul.kw...)
    Eu = Dict{Symbol,Matrix{Float64}}()
    Q = Dict{Symbol,Matrix{Float64}}()
    qs = Any[]
    if ul !== nothing
        du = [unit_light([b for b in beams if key(b) == k]) for k in ks]
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
    Y = ylm_weights(t, obs; lmax)
    # weighted rows: (index, weight) of the shape terms
    rows = Tuple{Int,Float64}[]
    for (l, w) in modes, mm in -l:l
        push!(rows, (mode_index(l, mm), w))
    end
    for ((l, mm), w) in components
        push!(rows, (mode_index(l, mm), w))
    end
    lim(k, τ) = limits === nothing ? (0.0, Inf) :
                limits isa AbstractDict ? get(limits, k, (0.0, Inf)) : limits(k, τ)
    nt = length(ref.time)
    power = zeros(ng, nt)
    shape = zeros(nm, nt)
    drive = zeros(nt)
    residual = zeros(nt)
    rshape = ref.ylm[1:nm, :] ./ ref.ylm[1:1, :]
    prev = nothing
    for (k, τ) in enumerate(ref.time)
        # the flux coefficients at powers P: c = M P (+ c0, the prepulse's)
        c0 = nothing
        if ul === nothing
            M = flux_response(m, obs; time = τ, albedo = ref.albedo[k], weights = Y) * E   # nm × ng
        else
            # each color's response; its conversion independent of intensity (as the 3ω
            # response requires), so the secant's powers do not enter
            R = flux_response(m, obs; time = τ, albedo = ref.albedo[k], weights = Y,
                              unconverted = Dict(c.name => zeros(ne) for c in ul.colors), colors = ul.colors)
            M = R[:omega3] * E
            for c in ul.colors
                rc = _fraction_at(ul.ratio[c.name], τ)
                rc >= 0 || throw(ArgumentError("the unconverted fractions must not be negative (at $τ ns)"))
                M .+= rc .* (R[c.name] * Eu[c.name])
            end
            if !isempty(qs)
                c0 = zeros(nm)
                for (j, c) in enumerate(ul.colors)
                    w = Float64[q[j] === nothing ? 0.0 : Float64(q[j](τ)) for q in qs]
                    all(>=(0), w) || throw(ArgumentError("the prepulse power must not be negative (at $τ ns)"))
                    c0 .+= R[c.name] * (Q[c.name] * w)
                end
            end
        end
        r00 = ref.ylm[1, k]
        s = rshape[:, k]
        A = zeros(0, ng)
        b = Float64[]
        for (i, w) in rows
            A = [A; sqrt(w) .* (M[i, :] .- s[i] .* M[1, :])' ./ r00]
            push!(b, c0 === nothing ? 0.0 : -sqrt(w) * (c0[i] - s[i] * c0[1]) / r00)
        end
        A = [A; sqrt(drive_weight) .* M[1, :]' ./ r00]
        push!(b, c0 === nothing ? sqrt(drive_weight) : sqrt(drive_weight) * (r00 - c0[1]) / r00)
        bounds = [lim(kk, τ) for kk in ks]
        lo, hi = first.(bounds), last.(bounds)
        all(lo .<= hi) || throw(ArgumentError("empty power limits at $τ ns"))
        if stay > 0
            Pn = [nom(kk, τ) for kk in ks]
            scale = max(sum(Pn) / ng, eps())
            A = [A; sqrt(stay) .* Matrix{Float64}(I, ng, ng) ./ scale]
            append!(b, sqrt(stay) .* Pn ./ scale)
        end
        if smooth > 0 && prev !== nothing
            scale = max(maximum(abs, prev), eps())
            A = [A; sqrt(smooth) .* Matrix{Float64}(I, ng, ng) ./ scale]
            append!(b, sqrt(smooth) .* prev ./ scale)
        end
        x0 = prev === nothing ? clamp.(fill(r00 / max(sum(M[1, :]), eps()), ng), lo, hi) : clamp.(prev, lo, hi)
        P = _bounded_lsq(A, b, lo, hi, x0)
        power[:, k] = P
        c = c0 === nothing ? M * P : M * P .+ c0
        shape[:, k] = c ./ c[1]
        drive[k] = c[1] / r00
        residual[k] = norm(A * P - b)
        prev = P
    end
    pulses = Dict{Any,PulseShape}(k => PulseShape(ref.time, power[g, :]) for (g, k) in enumerate(ks))
    return (; time = collect(ref.time), keys = ks, power, pulses, shape, reference_shape = rshape,
            drive, residual)
end

# argmin |A x - b|² for lo ≤ x ≤ hi (small, dense): L-BFGS-B on the quadratic.
function _bounded_lsq(A, b, lo, hi, x0)
    AtA, Atb = A' * A, A' * b
    f(x) = dot(x, AtA * x) - 2dot(Atb, x)
    g!(G, x) = (G .= 2 .* (AtA * x .- Atb))
    any(lo .== hi) && all(lo .== hi) && return copy(lo)
    res = Optim.optimize(f, g!, lo, hi, clamp.(x0, lo, hi), Optim.LBFGSB(),
                         Optim.Options(; g_abstol = 1e-14, iterations = 2000))
    return Optim.minimizer(res)
end
