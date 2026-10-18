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
(`:scene`, the 3D target, `:map`, the capsule θ–φ map or foil face, or `:image`, an
x-ray image), `spot_radius` (mm; the foil spot over which irradiance statistics and
the foil objective are taken), `status`, the VisRad benchmark kit's
`kit_running`/`kit_status` (see [`write_visrad_kit!`](@ref)), the retune's
`retune_running`/`retune_status` (see [`retune_to_reference!`](@ref)), the x-ray image's
`image_request` (its camera, see [`xray_image!`](@ref)) and `image` (the latest
rendering, `(; x, y, values, title)` in mm and W/cm²/sr, or `nothing`), `side`
(show the map or the image beside the 3D target instead of in its place), and the
cavity (see [`set_cavity!`](@ref)): `cavity` (compute the drive on the elements that
bound the free space around `cavity_point` only, `Briareus.enclosure`: an imported CAD
assembly's hohlraum cavity without its outer faces and mounting hardware), `films`
(components the cavity leaves out: thin films a CAD file models as solids, such as a
capsule tent or an entrance-hole window, which would seal the point into a pocket) and
`enclosure` (the latest cavity selection, with `cavity_index`, the target element of
each drive element).
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
    retune_running::Observable{Bool}
    retune_status::Observable{String}
    image_request::Observable{Any}           # the x-ray image's camera, or nothing
    image::Observable{Any}                   # its latest rendering, or nothing
    side::Observable{Bool}                   # map or image beside the 3D view
    cavity::Observable{Bool}                 # drive on the cavity around cavity_point
    cavity_point::Observable{Vec3}
    films::Observable{Set{Symbol}}           # left out of the cavity (tents, windows)
    enclosure::Base.RefValue{Any}            # the latest cavity selection, or nothing
    cavity_index::Base.RefValue{Any}         # target element of each drive element
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
                          Observable(false), Observable(""), Observable{Any}(nothing), Observable{Any}(nothing), Observable(false),
                          Observable(false), Observable(zero(Vec3)), Observable(Set{Symbol}()),
                          Ref{Any}(nothing), Ref{Any}(nothing),
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
    set_side_view!(dr::DriveState, on::Bool)

Show the flux map (or the x-ray image) beside the 3D target (`true`), or in its place
when chosen (`false`).
"""
set_side_view!(dr::DriveState, on::Bool) = (_assign!(dr.side, on); dr)

"""
    set_view!(dr::DriveState, view::Symbol)

Show the 3D target (`:scene`), the flux map (`:map`) or the x-ray image (`:image`, see
[`xray_image!`](@ref)).
"""
function set_view!(dr::DriveState, view::Symbol)
    view in (:scene, :map, :image) || throw(ArgumentError("view must be :scene, :map or :image"))
    _assign!(dr.view, view)
    return dr
end

# The x-ray image of the latest drive for the requested camera (`nothing` without one).
function _render_image(dr::DriveState)
    req, res = dr.image_request[], dr.result[]
    (req === nothing || res === nothing || res.solution === nothing) && return nothing
    sol = res.solution
    t = sol.model.target
    cam = OrthoCamera(req.direction; up = req.up, center = req.center, width = req.width,
                      nx = req.nx)
    img = render(t, radiance(sol), cam)
    h = cam.width / 2
    x = range(-h, h; length = cam.nx)
    y = range(-h, h; length = cam.ny)
    # rows from the top: as a matrix over (x, y) with y increasing
    values = permutedims(to_Wcm2.(img)[end:-1:1, :])
    return (; x, y, values, title = "X-ray radiance seen along " *
                                    join((@sprintf("%.2f", c) for c in cam.dir), ", "))
end

"""
    set_cavity!(dr::DriveState, on::Bool; point = dr.cavity_point[], films = dr.films[]) -> dr

Compute the drive on the cavity around `point` (mm) only (`on`), or on the whole target.
The cavity is what `Briareus.enclosure` keeps: the elements bounding the free space the
point is in, without the outer faces, end faces and mounting hardware of an imported CAD
assembly (a NIF CAD workspace's 10⁴–10⁵ elements become a few thousand), so its x-ray
view factors are within reach. `films` (component names) are left out first: thin
films a CAD file models as closed solids, such as a capsule tent or an entrance-hole
window, seal the point into a pocket (they transmit x-rays at their thickness). The
drive's laser deposition then sees the cavity only (light that would clip outer faces
passes); the clearance and TaLIS checks keep using the whole target.
"""
function set_cavity!(dr::DriveState, on::Bool; point = dr.cavity_point[], films = dr.films[])
    _assign!(dr.cavity_point, Vec3(point...))
    _assign!(dr.films, Set{Symbol}(films))
    _assign!(dr.cavity, on)
    return dr
end

# The cavity of a target around `point`, `films` left out: the enclosure, and the
# target element of each of its elements.
function _cavity(t::Target, point, films)
    comps = [c for c in t.components if !(c.name in films)]
    isempty(comps) && throw(ArgumentError("every component is a film left out of the cavity"))
    enc = enclosure(comps, point)
    src = enc.source
    pos = Dict{Symbol,Vector{Int}}(c.name => collect(elements(t, c.name)) for c in comps)
    rank = zeros(Int, nelements(src))
    for c in comps
        for (k, e) in enumerate(elements(src, c.name))
            rank[e] = k
        end
    end
    index = [pos[src.components[src.elem_comp[e]].name][rank[e]] for e in enc.index]
    return enc, index
end

"""
    cavity_rows(dr::DriveState) -> Vector{NamedTuple}

What the latest cavity selection kept of each component: rows `(; name, kept, elements,
area_kept, area)` (element counts and mm²), largest kept area first; empty before one.
"""
function cavity_rows(dr::DriveState)
    enc = dr.enclosure[]
    enc === nothing && return NamedTuple[]
    src = enc.source
    rows = map(src.components) do c
        es = elements(src, c.name)
        k = [e for e in es if enc.kept[e]]
        (; name = c.name, kept = length(k), elements = length(es),
         area_kept = sum(src.elem_area[e] for e in k; init = 0.0),
         area = sum(src.elem_area[e] for e in es; init = 0.0))
    end
    return sort(rows; by = r -> -r.area_kept)
end

"""
    cavity_summary(dr::DriveState) -> String

The latest cavity in a few lines: elements and area kept, the open area, and per
component what was kept; with a warning when the point sits in a sealed pocket.
"""
function cavity_summary(dr::DriveState)
    enc = dr.enclosure[]
    enc === nothing && return ""
    rows = cavity_rows(dr)
    n = sum(r.kept for r in rows; init = 0)
    A = sum(r.area_kept for r in rows; init = 0.0)
    head = @sprintf("Cavity around (%.3g, %.3g, %.3g) mm: %d of %d elements, %.4g mm², open area %.3g mm²",
                    enc.point..., n, nelements(enc.source), A, enc.open_area)
    sealed = cavity_sealed(dr) ?
             "\nSealed: nothing escapes, and few components bound it. A thin film modeled as a solid " *
             "(a capsule tent, an entrance-hole window) closes the point in: mark it as a thin film " *
             "(Properties)." : ""
    lines = [@sprintf("%s: %d of %d elements, %.3g of %.3g mm²", r.name, r.kept, r.elements,
                      r.area_kept, r.area) for r in rows if r.kept > 0]
    dropped = [String(r.name) for r in rows if r.kept == 0]
    isempty(dropped) || push!(lines, "dropped whole: " * join(dropped, ", "))
    return head * sealed * "\n" * join(lines, "\n")
end

"Whether the latest cavity is a sealed pocket (no open area: a film closes the point in)."
cavity_sealed(dr::DriveState) = dr.enclosure[] !== nothing && _sealed(dr.enclosure[])

# A pocket: nothing escapes from the point, and at most two components bound it (a
# capsule in its tent, say).
_sealed(enc) = enc.escape < 1e-3 && length(enc.target.components) <= 2

# The cavity around a point, and when films seal it into a pocket, the component of the
# pocket whose leaving out opens the largest cavity (the film: a tent or a window, tried
# among the pocket's three largest components). Returns (enclosure, index, film or nothing).
function _open_cavity(t::Target, point, films)
    enc, index = _cavity(t, point, films)
    _sealed(enc) || return enc, index, nothing
    best = nothing
    for c in sort(enc.target.components; by = c -> -sum(enc.target.elem_area[e] for e in elements(enc.target, c.name)))[1:min(3, end)]
        e2, i2 = _cavity(t, point, union(films, (c.name,)))
        _sealed(e2) && continue
        A = sum(e2.target.elem_area)
        (best === nothing || A > best[4]) && (best = (e2, i2, c.name, A))
    end
    best === nothing && return enc, index, nothing
    return best[1], best[2], best[3]
end

# Build a scenario for a target and the table's beams (all of them: disabled beams are
# evaluated at zero power, so on/off needs no rebuild). Hohlraums observe the capsule;
# `closed` (a hohlraum and nothing outside it) balances the view factors of the
# enclosure. Every target gets a radiation model for the B/H fields, unless it is too
# large for its view factors (`radiation_size`: laser only). With `cavity = (point,
# films)` the scenario is built on the cavity around the point (see `set_cavity!`); a
# pocket sealed by a film is opened by leaving that film out too (`_open_cavity`), and
# a point in the open (most of its rays escape: a backlighter foil and its shields, say)
# keeps the whole target. Returns `(; s, lookup, enc, index, film, escape)`: the
# scenario, its flux-map lookup, the cavity's enclosure and element index, the film
# found (`nothing`s without a cavity), and the point's escaping fraction.
function _build_scenario(t::Target, beams::Vector{Beam{Float64}}; closed::Bool = false,
                         radiation::Union{Bool,Nothing} = nothing, cavity = nothing)
    enc = index = film = nothing
    escape = NaN
    if cavity !== nothing
        # the parts are meshed to their element sizes already, so the selection runs on
        # elements (`visrad_target(...; inside)` selects on CAD facets before remeshing)
        enc, index, film = _open_cavity(t, cavity...)
        escape = enc.escape
        if escape > CAVITY_MAX_ESCAPE
            enc = index = film = nothing
        else
            t = enc.target
            closed = false
        end
    end
    radiation === nothing && (radiation = radiation_size(t).feasible)
    names = component_names(t)
    observe = radiation && :capsule in names ? :capsule : nothing
    s = Scenario(t, beams; observe, radiation, viewfactor_kw = (; closed))
    return (; s, lookup = flux_lookup(s), enc, index, film, escape)
end

"""
The largest fraction of the rays from the cavity point that may escape for the drive to
use the cavity: above it the point is in the open (a backlighter foil's side of the
chamber) and the drive keeps the whole target.
"""
const CAVITY_MAX_ESCAPE = 0.5

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

# The capsule flux comes from the drive's radiosity solution (`ev.flux` would solve for
# the full flux response, which only repeated evaluations repay).
function _flux_map(lk, ev::Evaluation, sol)
    lk === nothing && return nothing
    t = ev.scenario.target
    vals = lk.kind === :capsule ? to_Wcm2.(sol.H[ev.scenario.observed]) :
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
    if _has_flux(ev)
        c = flux_ylm(ev; lmax = maximum(MODE_DEGREES))
        a = flux_moments(ev; lmax = maximum(MODE_DEGREES))
        mean = a[1]
        modes = [(; l, epsilon = mode_amplitude(c, l), legendre = a[l+1] / a[1]) for l in MODE_DEGREES]
    end
    return DriveResult(ev, sol, modes, mean, power_rows(ev), _spot_stats(ev, spot_radius),
                       _flux_map(lookup, ev, sol), time() - t0)
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
    return @sprintf("%.4g TW delivered", d.delivered) *
           (res.solution === nothing ? " (laser only: no x-ray drive)" : "")
end
