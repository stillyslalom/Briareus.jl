# The Drive step: the Scenario (deposition operator, view factors and capsule flux
# response; built once per target and beam ports, on a worker) and its evaluation for
# the current pointing (deposition, radiosity solve, capsule modes, power accounting).
# Pointing, power and on/off edits only re-evaluate; the scenario is rebuilt when the
# target or the beam ports change.

"""
Element fields the 3D view colors the target by (all W/cm²): `:laser` (absorbed
laser intensity), `:radiosity` (x-ray radiosity B), `:incident` (incident x-ray flux H).
"""
const DRIVE_FIELDS = (:laser, :radiosity, :incident)

const DRIVE_FIELD_LABELS = Dict(:laser => "Absorbed laser intensity (W/cm²)",
                                :radiosity => "X-ray radiosity (W/cm²)",
                                :incident => "Incident x-ray flux (W/cm²)")

"Degrees of the mode table: ε_l and P_l/P0 for l = 1…4."
const MODE_DEGREES = 1:4

"""
    DriveResult

One evaluation of the drive: `evaluation` (core `Evaluation` of the display
beams: every beam of the table, disabled ones at zero power), `solution` (the
radiosity `RadiationSolution`, or `nothing` without a radiation model), `modes` (rows
`(; l, epsilon, legendre)` of the capsule flux, empty without a capsule), `mean_flux`
(TW/mm², `NaN` without a capsule), `power` (rows `(; name, power, fraction)`), `spot`
(foil irradiance statistics, or `nothing`), `map` (see [`flux_map`](@ref)), and
`seconds`.
"""
struct DriveResult
    evaluation::Evaluation
    solution::Any
    modes::Vector{NamedTuple}
    mean_flux::Float64
    power::Vector{NamedTuple}
    spot::Any
    map::Any
    seconds::Float64
end

"""
    DriveState()

State of the drive analysis: `scenario` (the core `Scenario`, `nothing` before the first
build), `building`/`started`/`build_seconds`, `result` (the latest
[`DriveResult`](@ref)), `evaluating`, `field` (one of [`DRIVE_FIELDS`](@ref)), `view`
(`:scene`, the 3D target, or `:map`, the capsule θ–φ map or foil face), `spot_radius`
(mm; the foil spot over which irradiance statistics and the foil objective are
taken), `status`, and the VisRad benchmark kit's `kit_running`/`kit_status` (see
[`write_visrad_kit!`](@ref)).
"""
struct DriveState
    scenario::Observable{Union{Nothing,Scenario}}
    building::Observable{Bool}
    started::Observable{Float64}
    build_seconds::Observable{Float64}
    result::Observable{Union{Nothing,DriveResult}}
    evaluating::Observable{Bool}
    field::Observable{Symbol}
    view::Observable{Symbol}
    spot_radius::Observable{Float64}
    status::Observable{String}
    kit_running::Observable{Bool}
    kit_status::Observable{String}
    key::Base.RefValue{Any}                  # (target, ports) of the scenario
    lookup::Base.RefValue{Any}               # flux-map raster lookup of the scenario
    generation::Base.RefValue{Int}           # build requests
    eval_generation::Base.RefValue{Int}      # evaluation requests
end

DriveState() = DriveState(Observable{Union{Nothing,Scenario}}(nothing), Observable(false),
                          Observable(0.0), Observable(0.0),
                          Observable{Union{Nothing,DriveResult}}(nothing), Observable(false),
                          Observable(:laser), Observable(:scene), Observable(0.5),
                          Observable(""), Observable(false), Observable(""),
                          Ref{Any}(nothing), Ref{Any}(nothing), Ref(0), Ref(0))

"""
    set_field!(dr::DriveState, field::Symbol)

Color the target by one of [`DRIVE_FIELDS`](@ref).
"""
function set_field!(dr::DriveState, field::Symbol)
    field in DRIVE_FIELDS || throw(ArgumentError("field must be one of $DRIVE_FIELDS"))
    _assign!(dr.field, field)
    return dr
end

"""
    set_view!(dr::DriveState, view::Symbol)

Show the 3D target (`:scene`) or the flux map (`:map`).
"""
function set_view!(dr::DriveState, view::Symbol)
    view in (:scene, :map) || throw(ArgumentError("view must be :scene or :map"))
    _assign!(dr.view, view)
    return dr
end

# Build a scenario for a target and the table's beams (all of them: disabled beams are
# evaluated at zero power, so on/off needs no rebuild). Hohlraums observe the capsule;
# `closed` (a hohlraum and nothing outside it) balances the view factors of the
# enclosure. Every target gets a radiation model for the B/H fields.
function _build_scenario(t::Target, beams::Vector{Beam{Float64}}; closed::Bool = false)
    names = component_names(t)
    observe = :capsule in names ? :capsule : nothing
    s = Scenario(t, beams; observe, radiation = true, viewfactor_kw = (; closed))
    return s, flux_lookup(s)
end

# ---------------------------------------------------------------- flux map

"""
    flux_lookup(s::Scenario; ntheta = 90, nphi = 180, npix = 100)

Raster lookup for [`flux_map`](@ref): for a scenario observing a capsule, the observed
element seen in each direction of a θ–φ grid about the capsule center; for a foil
target, the front element nearest each pixel of the foil face (y–z plane). `nothing`
otherwise.
"""
function flux_lookup(s::Scenario; ntheta::Integer = 90, nphi::Integer = 180, npix::Integer = 100)
    t = s.target
    if !isempty(s.observed)
        obs = s.observed
        A = t.elem_area[obs]
        c0 = sum(A .* t.elem_centroid[obs]) / sum(A)
        dirs = [normalize(t.elem_centroid[e] - c0) for e in obs]
        θs = range(180 / 2ntheta, 180 - 180 / 2ntheta; length = ntheta)
        φs = range(360 / 2nphi, 360 - 360 / 2nphi; length = nphi)
        idx = Matrix{Int}(undef, nphi, ntheta)
        for (j, θ) in enumerate(θs), (i, φ) in enumerate(φs)
            d = port_vector(θ, φ)
            idx[i, j] = argmax(k -> dot(dirs[k], d), eachindex(dirs))
        end
        return (; kind = :capsule, x = collect(φs), y = collect(θs), index = idx, elems = obs,
                xlabel = "φ (°)", ylabel = "θ (°)",
                title = "Capsule incident x-ray flux (W/cm²)")
    end
    FOIL_NAME in component_names(t) || return nothing
    es = [e for e in elements(t, FOIL_NAME) if t.elem_side[e] == 1]
    us = [t.elem_centroid[e][2] for e in es]
    vs = [t.elem_centroid[e][3] for e in es]
    ylo, yhi = extrema(us)
    zlo, zhi = extrema(vs)
    pad = (yhi - ylo) / (2 * sqrt(length(es)))        # half an element
    ys = range(ylo - pad, yhi + pad; length = npix)
    zs = range(zlo - pad, zhi + pad; length = npix)
    idx = Matrix{Int}(undef, npix, npix)
    for (j, z) in enumerate(zs), (i, y) in enumerate(ys)
        idx[i, j] = argmin(k -> (us[k] - y)^2 + (vs[k] - z)^2, eachindex(es))
    end
    return (; kind = :foil, x = collect(ys), y = collect(zs), index = idx, elems = es,
            xlabel = "y (mm)", ylabel = "z (mm)", title = "Foil laser irradiance (W/cm²)")
end

"""
    flux_map(dr::DriveState) -> Union{Nothing,NamedTuple}

The flux map of the latest evaluation: `(; kind, x, y, values, xlabel, ylabel,
title)` with `values[i, j]` at `(x[i], y[j])` in W/cm² — capsule incident x-ray flux
over polar angle θ and azimuth φ (degrees) for a hohlraum, absorbed laser irradiance
over the foil face for a foil. `nothing` before an evaluation or for targets with
neither.
"""
flux_map(dr::DriveState) = dr.result[] === nothing ? nothing : dr.result[].map

function _flux_map(lk, ev::Evaluation)
    lk === nothing && return nothing
    t = ev.scenario.target
    vals = lk.kind === :capsule ? to_Wcm2.(ev.flux) :
           to_Wcm2.(ev.deposition.element[lk.elems] ./ t.elem_area[lk.elems])
    return (; lk.kind, lk.x, lk.y, values = vals[lk.index], lk.xlabel, lk.ylabel, lk.title)
end

# ---------------------------------------------------------------- evaluation

# Evaluate beams on a scenario: deposition, radiosity, modes, power accounting, map.
function _evaluate_drive(s::Scenario, beams, lookup, spot_radius::Real)
    t0 = time()
    ev = evaluate(s, beams)
    sol = s.model === nothing ? nothing : solve(s.model, ev.deposition)
    modes = NamedTuple[]
    mean = NaN
    if ev.flux !== nothing
        c = flux_ylm(ev; lmax = maximum(MODE_DEGREES))
        a = flux_moments(ev; lmax = maximum(MODE_DEGREES))
        mean = a[1]
        modes = [(; l, epsilon = mode_amplitude(c, l), legendre = a[l+1] / a[1]) for l in MODE_DEGREES]
    end
    return DriveResult(ev, sol, modes, mean, power_rows(ev), _spot_stats(ev, spot_radius),
                       _flux_map(lookup, ev), time() - t0)
end

"""
    power_rows(ev::Evaluation) -> Vector{NamedTuple}

Laser power accounting rows `(; name, power, fraction)` (TW, fraction of delivered):
each component's absorbed power, back-side (outside) losses, and missed power.
"""
function power_rows(ev::Evaluation)
    t = ev.scenario.target
    d = ev.deposition
    P = d.delivered
    cp = component_power(d, t)
    rows = NamedTuple[]
    for c in t.components
        c.kind === :aperture && continue
        v = get(cp, c.name, 0.0)
        push!(rows, (; name = String(c.name), power = v, fraction = P > 0 ? v / P : 0.0))
    end
    for c in t.components
        k = Symbol(c.name, :_outside)
        haskey(cp, k) && push!(rows, (; name = "$(c.name) (outside)", power = cp[k],
                                      fraction = P > 0 ? cp[k] / P : 0.0))
    end
    m = get(cp, :missed, 0.0)
    abs(m) > 1e-6 * P &&
        push!(rows, (; name = m > 0 ? "missed" : "imbalance", power = m, fraction = P > 0 ? m / P : 0.0))
    return rows
end

# Foil irradiance within `radius` of the foil center: mean, peak (W/cm²) and RMS
# nonuniformity; `nothing` without a foil.
function _spot_stats(ev::Evaluation, radius::Real)
    t = ev.scenario.target
    FOIL_NAME in component_names(t) || return nothing
    front = [e for e in elements(t, FOIL_NAME) if t.elem_side[e] == 1]
    A = t.elem_area[front]
    c0 = sum(A .* t.elem_centroid[front]) / sum(A)
    es = [e for e in front if norm(t.elem_centroid[e] - c0) <= radius]
    isempty(es) && return nothing
    I = ev.deposition.element[es] ./ t.elem_area[es]
    Ae = t.elem_area[es]
    mean = sum(I .* Ae) / sum(Ae)
    return (; radius = Float64(radius), mean = to_Wcm2(mean), peak = to_Wcm2(maximum(I)),
            rms = mean > 0 ? rms_nonuniformity(I, Ae) : NaN,
            fraction = sum(ev.deposition.element[es]) / ev.deposition.delivered)
end

"""
    element_values(res::DriveResult, field::Symbol) -> Vector{Float64}

Per-element values (W/cm²) of one of [`DRIVE_FIELDS`](@ref); `NaN` where the field
is undefined (apertures for the laser field; every element without a radiation
model).
"""
function element_values(res::DriveResult, field::Symbol)
    ev = res.evaluation
    t = ev.scenario.target
    if field === :laser
        v = to_Wcm2.(laser_flux(ev.deposition, t))
        for e in 1:nelements(t)
            t.components[t.elem_comp[e]].kind === :aperture && (v[e] = NaN)
        end
        return v
    end
    field in DRIVE_FIELDS || throw(ArgumentError("field must be one of $DRIVE_FIELDS"))
    res.solution === nothing && return fill(NaN, nelements(t))
    return to_Wcm2.(field === :radiosity ? res.solution.B : res.solution.H)
end

"""
    drive_summary(res::DriveResult) -> String

The drive in a line: capsule mean flux, radiation temperature, and ε₂/ε₄; or the foil
spot irradiance and nonuniformity.
"""
function drive_summary(res::DriveResult)
    if !isempty(res.modes)
        m = Dict(r.l => r.epsilon for r in res.modes)
        return @sprintf("capsule %s W/cm² (Tr %.0f eV) · ε₂ %.2f %% · ε₄ %.2f %%",
                        _wcm2(res.mean_flux), radiation_temperature(res.mean_flux),
                        100 * get(m, 2, NaN), 100 * get(m, 4, NaN))
    elseif res.spot !== nothing
        s = res.spot
        return @sprintf("spot %.3g W/cm² · rms %.1f %% · %.0f %% of the power in r ≤ %.2g mm",
                        s.mean, 100 * s.rms, 100 * s.fraction, s.radius)
    end
    d = res.evaluation.deposition
    return @sprintf("%.4g TW delivered", d.delivered)
end
