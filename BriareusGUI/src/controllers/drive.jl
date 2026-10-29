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

"""
Beam models of the drive's laser deposition: `:fcone` (Briareus's: the focal spot
convolved over the f-cone) and `:visrad` (VisRad's, `VisRadProfile`, so that a design
exported to VisRad deposits there as here).
"""
const BEAM_MODELS = (:fcone, :visrad)

const BEAM_MODEL_LABELS = Dict(:fcone => "Briareus (f-cone)", :visrad => "VisRad profile")

"Degrees of the mode table: ε_l and P_l/P0 for l = 1…4."
const MODE_DEGREES = 1:4

"""
    DriveResult

One evaluation of the drive: `evaluation` (core `Evaluation` of the display
beams: every beam of the table, disabled ones at zero power), `solution` (the
radiosity `RadiationSolution`, or `nothing` without a radiation model), `modes` (rows
`(; l, epsilon, legendre)` of the capsule flux, empty without a capsule), `mean_flux`
(TW/mm², `NaN` without a capsule), `power` (rows `(; name, power, fraction)`), `spot`
(foil irradiance statistics, or `nothing`), `map` (see [`flux_map`](@ref)),
`seconds`, and `unconverted` (`nothing`, or the unconverted light's share when it takes
part: `(; omega1, omega2)`, each `(; absorbed, delivered, conversion)` (TW, and its x-ray
conversion relative to 3ω), and `skipped`, the beams without a footprint; see
[`set_unconverted_ratio!`](@ref)). With unconverted light, `evaluation`'s deposition is
the absorbed power of the 3ω light plus the 1ω and 2ω `unconverted_deposition`, and the
x-ray source (`solution`, the capsule flux) counts each unconverted color's power at the
material's conversion times the color's `conversion` (`UnconvertedColor`).
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
    unconverted::Any
end

"""
    DriveState()

State of the drive analysis: `scenario` (the core `Scenario`, `nothing` before the first
build), `building`/`started`/`build_seconds`, `result` (the latest
[`DriveResult`](@ref)), `evaluating`, `field` (one of [`DRIVE_FIELDS`](@ref)), `view`
(`:scene`, the 3D target, `:map`, the capsule θ–φ map or foil face, or `:image`, an
x-ray image), `spot_radius` (mm; the foil spot over which irradiance statistics and
the foil objective are taken), `status`, the retune's
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
each drive element), and the beam model (see [`set_beam_model!`](@ref)): `beam_model`
(one of [`BEAM_MODELS`](@ref)) and `profile` (the `VisRadProfile` of `:visrad`), and
`unconverted_ratio` (the 1ω and 2ω power as fractions of each beam's 3ω power; 0: that
color off; see [`set_unconverted_ratio!`](@ref)).
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
    retune_running::Observable{Bool}
    retune_status::Observable{String}
    image_request::Observable{Any}           # the x-ray image's camera, or nothing
    image::Observable{Any}                   # its latest rendering, or nothing
    side::Observable{Bool}                   # map or image beside the 3D view
    cavity::Observable{Bool}                 # drive on the cavity around cavity_point
    cavity_point::Observable{Vec3}
    films::Observable{Set{Symbol}}           # left out of the cavity (tents, windows)
    beam_model::Observable{Symbol}           # :fcone or :visrad
    profile::Observable{VisRadProfile}       # the deposition profile of :visrad
    enclosure::Base.RefValue{Any}            # the latest cavity selection, or nothing
    cavity_index::Base.RefValue{Any}         # target element of each drive element
    key::Base.RefValue{Any}                  # (target, ports) of the scenario
    lookup::Base.RefValue{Any}               # flux-map raster lookup of the scenario
    generation::Base.RefValue{Int}           # build requests
    eval_generation::Base.RefValue{Int}      # evaluation requests
    unconverted_ratio::Observable{NTuple{2,Float64}}   # 1ω, 2ω power / 3ω power
end

DriveState() = DriveState(Observable{Union{Nothing,Scenario}}(nothing), Observable(false),
                          Observable(0.0), Observable(0.0),
                          Observable{Union{Nothing,DriveResult}}(nothing), Observable(false),
                          Observable(:laser), Observable(:scene), Observable(0.5),
                          Observable(""),
                          Observable(false), Observable(""), Observable{Any}(nothing), Observable{Any}(nothing), Observable(false),
                          Observable(false), Observable(zero(Vec3)), Observable(Set{Symbol}()),
                          Observable(:fcone), Observable(VisRadProfile()),
                          Ref{Any}(nothing), Ref{Any}(nothing),
                          Ref{Any}(nothing), Ref{Any}(nothing), Ref(0), Ref(0), Observable((0.0, 0.0)))

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
    set_beam_model!(dr::DriveState, model::Symbol; profile = dr.profile[]) -> dr

Deposit the laser with one of [`BEAM_MODELS`](@ref): `:visrad` with `profile` (a
`VisRadProfile`; a VisRad workspace import sets the workspace's, `visrad_profile`),
`:fcone` with Briareus's model. `profile` is kept for a later switch back to
`:visrad`. A change makes the scenario stale.
"""
function set_beam_model!(dr::DriveState, model::Symbol; profile::VisRadProfile = dr.profile[])
    model in BEAM_MODELS || throw(ArgumentError("beam model must be one of $BEAM_MODELS"))
    _profile_state(profile) == _profile_state(dr.profile[]) || (dr.profile[] = profile)
    _assign!(dr.beam_model, model)
    return dr
end

"""
    set_unconverted_ratio!(dr::DriveState, omega1, omega2) -> dr

The unconverted 1ω and 2ω power of every beam as a fraction of its 3ω power (≥ 0; 0
turns the color off). The beams whose 1ω or 2ω light takes part
([`set_participation!`](@ref)) then add `unconverted_deposition`'s absorbed power to the
drive. There is no default: the NIF User Guide gives no fraction (most of the light is
1ω, more for shaped pulses), so the ratio is an input, as in VisRad.
"""
function set_unconverted_ratio!(dr::DriveState, omega1::Real, omega2::Real)
    (isfinite(omega1) && isfinite(omega2) && omega1 >= 0 && omega2 >= 0) ||
        throw(ArgumentError("the unconverted power ratios must be finite and ≥ 0"))
    _assign!(dr.unconverted_ratio, (Float64(omega1), Float64(omega2)))
    return dr
end

# The deposition profile of the drive's scenario (`nothing`: the f-cone model).
_deposition_profile(dr::DriveState) = dr.beam_model[] === :visrad ? dr.profile[] : nothing

# A `VisRadProfile` as plain data (session files, and the scenario's key: profiles
# compare by value). TOML has no `nothing`: an absent envelope is VisRad's default
# ratio.
function _profile_state(p::VisRadProfile)
    d = Dict{String,Any}()
    p.envelope === nothing ||
        (d["envelope"] = p.envelope isa AbstractDict ? Dict{String,Any}(p.envelope) : p.envelope)
    return d
end
_profile_state(::Nothing) = nothing

function _state_profile(d::AbstractDict)
    e = get(d, "envelope", nothing)
    envelope = e isa AbstractDict ? Dict{String,Float64}(String(k) => Float64(v) for (k, v) in e) : e
    return VisRadProfile(; envelope)
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
    return (; x, y, values, title = "X-ray radiance viewed from " * _view_angles(-cam.dir))
end

# A direction as polar and azimuth angles (degrees), as ports are given; along ±z the
# azimuth is left out.
function _view_angles(d)
    θ = acosd(clamp(d[3] / norm(d), -1.0, 1.0))
    s = @sprintf("θ %.0f°", θ)
    hypot(d[1], d[2]) < 1e-3 * norm(d) && return s
    return s * @sprintf(", φ %.0f°", mod(atand(d[2], d[1]), 360))
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
passes); the clearance and TaLIS checks keep using the whole target. A CAD cavity is
open wherever its file is (entrance holes), so its view factors are left raw unless
aperture parts close its openings: the drive then balances them (`closed = true`,
`enforce_closure!`), when at most [`CAVITY_MAX_LEAK`](@ref) of the cavity's rays leak,
and its status says which.
"""
function set_cavity!(dr::DriveState, on::Bool; point = dr.cavity_point[], films = dr.films[])
    _assign!(dr.cavity_point, Vec3(point...))
    _assign!(dr.films, Set{Symbol}(films))
    _assign!(dr.cavity, on)
    return dr
end

# The cavity of a target around `point`, `films` (parts) left out: the enclosure, and the
# target element of each of its elements.
function _cavity(t::Target, point, films)
    comps = [c for c in t.components if !(c.part in films || c.name in films)]
    isempty(comps) && throw(ArgumentError("every component is a film left out of the cavity"))
    enc = enclosure(comps, point)
    src = enc.source
    own(tg, c) = tg.comp_elems[component_index(tg, c.name)]     # a component's own elements
    pos = Dict{Symbol,Vector{Int}}(c.name => own(t, c) for c in comps)
    rank = zeros(Int, nelements(src))
    for c in comps
        for (k, e) in enumerate(own(src, c))
            rank[e] = k
        end
    end
    index = [pos[src.components[src.elem_comp[e]].name][rank[e]] for e in enc.index]
    return enc, index
end

"""
Whether the drive fits the curved faces of CAD parts (cylinders, cones, spheres:
`curved_components` on each part's file facets) in its own target, so the x-ray view
factors integrate on the true surfaces. The window keeps showing and keying the parts;
the drive's elements map back to theirs. Set it to `false` for flat CAD elements.
"""
const FIT_CAD = Ref(true)

# The drive's target of a design's target `t`: CAD components without a surface split by
# `curved_components` on their file facets `refs` (component name => facets, chamber
# frame), and the element of `t` each of its elements is (`nothing` when nothing split:
# then `t` itself). The pieces keep their part's facets, so the elements match by side
# and centroid.
const _FIT_MEMO = Dict{Tuple{UInt,UInt},Pair{Component,Vector{Component}}}()

# The file facets of a design's CAD parts (STL, VisRad CAD) in the target `t` (one flat
# component named after the part), by name, in the chamber frame: what `_fitted_target`
# fits them on.
function _fit_references(t::Target, kind::Symbol, params, mesh, parts)
    refs = Dict{Symbol,Mesh}()
    FIT_CAD[] || return refs
    flat = Set(c.name for c in t.components if c.kind === :surface && c.surface === nothing)
    frames = isempty(parts) ? Dict{Symbol,Frame}() : part_frames(kind, params, parts)
    for p in parts
        p.kind in (:visrad, :mesh) && p.name in flat || continue
        try
            p.kind === :visrad && _visrad_object(p.file, p.object).kind !== :mesh && continue
            cs = part_components(p, mesh.h; frame = frames[p.name], raw = true)
            length(cs) == 1 && only(cs).name === p.name && (refs[p.name] = only(cs).mesh)
        catch err
            err isa ArgumentError || rethrow()
        end
    end
    return refs
end

function _fitted_target(t::Target, refs::AbstractDict)
    isempty(refs) && return t, nothing
    comps = Component[]
    split = false
    for c in t.components
        ref = get(refs, c.name, nothing)
        if ref === nothing || c.surface !== nothing || c.kind !== :surface
            push!(comps, c)
            continue
        end
        key = (objectid(c), objectid(ref))
        lock(() -> length(_FIT_MEMO) > 256 && empty!(_FIT_MEMO), _STL_LOCK)
        e = _cached(() -> c => curved_components(c; reference = ref), _FIT_MEMO, key)
        pieces = first(e) === c ? last(e) : curved_components(c; reference = ref)
        split |= !(length(pieces) == 1 && only(pieces) === c)
        append!(comps, pieces)
    end
    split || return t, nothing
    tf = Target(comps; closed = t.closed, visibility = t.visibility)
    at = Dict{Tuple{Symbol,Int8,Vec3},Int}()
    for e in 1:nelements(t)
        at[(t.components[t.elem_comp[e]].part, t.elem_side[e], t.elem_centroid[e])] = e
    end
    fmap = [get(at, (tf.components[tf.elem_comp[e]].part, tf.elem_side[e], tf.elem_centroid[e]), 0)
            for e in 1:nelements(tf)]
    all(>(0), fmap) && allunique(fmap) || return t, nothing      # (not expected: keep the flat target)
    return tf, fmap
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
    rows = map(Briareus.parts(src)) do p                  # a part split by the drive's fit as one
        es = elements(src, p)
        k = [e for e in es if enc.kept[e]]
        (; name = p, kept = length(k), elements = length(es),
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

# A pocket: nothing escapes from the point, and at most two parts bound it (a capsule in
# its tent, say).
_sealed(enc) = enc.escape < 1e-3 && length(Briareus.parts(enc.target)) <= 2

# The cavity around a point, and when films seal it into a pocket, the part of the
# pocket whose leaving out opens the largest cavity (the film: a tent or a window, tried
# among the pocket's three largest parts). Returns (enclosure, index, film or nothing).
function _open_cavity(t::Target, point, films)
    enc, index = _cavity(t, point, films)
    _sealed(enc) || return enc, index, nothing
    best = nothing
    tg = enc.target
    for p in sort(Briareus.parts(tg); by = p -> -sum(tg.elem_area[e] for e in elements(tg, p)))[1:min(3, end)]
        e2, i2 = _cavity(t, point, union(films, (p,)))
        _sealed(e2) && continue
        A = sum(e2.target.elem_area)
        (best === nothing || A > best[4]) && (best = (e2, i2, p, A))
    end
    best === nothing && return enc, index, nothing
    return best[1], best[2], best[3]
end

# Whether `t` has a part or component named `name` (`false` for `nothing`).
_has_part(t::Target, name) =
    name !== nothing && (name in component_names(t) || any(c -> c.part === name, t.components))

# The foil a flux-map lookup was made for (`nothing`: none, or a capsule map).
_lookup_foil(lk) = lk !== nothing && lk.kind === :foil ? lk.foil : nothing

# Build a scenario for a target and the table's beams (all of them: disabled beams are
# evaluated at zero power, so on/off needs no rebuild). Hohlraums observe the capsule;
# `closed` (a hohlraum and nothing outside it) balances the view factors of the
# enclosure. Every target gets a radiation model for the B/H fields, unless it is too
# large for its view factors (`radiation_size`: laser only). With `cavity = (point,
# films)` the scenario is built on the cavity around the point (see `set_cavity!`); a
# pocket sealed by a film is opened by leaving that film out too (`_open_cavity`), and
# a point in the open (most of its rays escape: a backlighter foil and its shields, say)
# keeps the whole target. The cavity is balanced when apertures close it
# (`_cavity_closure`), whatever `closed` says. `profile` is the deposition's beam model
# (`nothing`: the f-cone; a `VisRadProfile`: VisRad's, see `set_beam_model!`).
# `references` (CAD part name => file facets, `_fit_references`) splits those parts into
# their fitted curved pieces in the drive's target (`_fitted_target`). `foil`: the
# backlighter foil of the flux map (`active_foil`). Returns `(; s,
# lookup, enc, index, film, escape, closed, leak)`: the scenario, its flux-map lookup,
# the cavity's enclosure and the target element of each drive element (`nothing` when
# the drive's elements are the target's), the film found (`nothing` without a cavity),
# the point's escaping fraction, whether the view factors were balanced, and the
# cavity's leaking fraction (`open_fraction`; `NaN` without a cavity) and whether it
# has apertures, as `(; fraction, apertures)`.
function _build_scenario(t::Target, beams::Vector{Beam{Float64}}; closed::Bool = false,
                         radiation::Union{Bool,Nothing} = nothing, cavity = nothing,
                         profile::Union{Nothing,VisRadProfile} = nothing,
                         references::AbstractDict = Dict{Symbol,Mesh}(),
                         foil::Union{Nothing,Symbol} = nothing)
    enc = index = film = nothing
    escape = NaN
    leak = (; fraction = NaN, apertures = false)
    # (the fit serves the x-ray view factors only)
    t, fmap = radiation === false || !radiation_size(t).feasible && cavity === nothing ?
              (t, nothing) : _fitted_target(t, references)
    if cavity !== nothing
        # the parts are meshed to their element sizes already, so the selection runs on
        # elements (`visrad_target(...; inside)` selects on CAD facets before remeshing)
        enc, index, film = _open_cavity(t, cavity...)
        escape = enc.escape
        if escape > CAVITY_MAX_ESCAPE
            enc = index = film = nothing
        else
            t = enc.target
            closed, leak = _cavity_closure(t)
        end
    end
    fmap === nothing || (index = index === nothing ? fmap : fmap[index])
    radiation === nothing && (radiation = radiation_size(t).feasible)
    names = component_names(t)
    observe = radiation && :capsule in names ? :capsule : nothing
    s = Scenario(t, beams; observe, radiation, viewfactor_kw = (; closed), deposition_kw = (; profile))
    return (; s, lookup = flux_lookup(s; foil), enc, index, film, escape, closed, leak)
end

"""
The largest fraction of the rays from the cavity point that may escape for the drive to
use the cavity: above it the point is in the open (a backlighter foil's side of the
chamber) and the drive keeps the whole target.
"""
const CAVITY_MAX_ESCAPE = 0.5

"""
The largest fraction of a cavity's diffuse rays that may leak (`Briareus.open_fraction`
of the cavity selected, apertures stopping the rays) for the drive to balance the view
factors of a cavity that has apertures (`closed = true`, `enforce_closure!`). Adding
apertures is the user's declaration that they close the cavity's openings; this guard
keeps a cavity with an opening left without one, or with an aperture facing out of the
cavity (dropped by the selection), unbalanced. An engineering default, not a facility
criterion, set from the fractions measured in the research notes ("Closure on bored,
windowed and CAD cavities"): cracks where meshes meet 2–7 × 10⁻⁴ (2.7 × 10⁻³ on a coarse
hohlraum) and a 40 µm waist gap 4 × 10⁻³ lie below it; a CAD cavity's two open LEHs,
7 × 10⁻², lie well above it.
"""
const CAVITY_MAX_LEAK = 1e-2

# Whether the drive balances the view factors of the cavity `t`: when it has apertures
# (the user's declaration that they close it) and they do (`CAVITY_MAX_LEAK`). Returns
# `(closed, (; fraction, apertures))`.
function _cavity_closure(t::Target)
    apertures = any(c -> c.kind === :aperture, t.components)
    fraction = open_fraction(t).fraction
    return apertures && fraction <= CAVITY_MAX_LEAK, (; fraction, apertures)
end

# The drive status's note on the cavity's closure: balanced, or how much leaks and why.
function _closure_note(leak, closed::Bool)
    isnan(leak.fraction) && return ""
    pct = @sprintf("%.2g %%", 100leak.fraction)
    closed && return " · closed by its apertures: view factors balanced ($pct of the rays leak)"
    leak.apertures ||
        return " · open: $pct of the cavity's rays leak: add apertures to its LEHs to balance its view factors"
    return " · open: $pct of the cavity's rays leak past its apertures (above " *
           @sprintf("%.2g %%", 100CAVITY_MAX_LEAK) * "): an opening has no aperture, or one faces out of " *
           "the cavity (Properties ▸ Inward); view factors not balanced"
end

# ---------------------------------------------------------------- flux map

"""
    flux_lookup(s::Scenario; foil = nothing, ntheta = 90, nphi = 180, npix = 100)

Raster lookup for [`flux_map`](@ref): for a scenario observing a capsule, the observed
element seen in each direction of a θ–φ grid about the capsule center; otherwise, with
a backlighter foil `foil` (a part or component name of the scenario's target, see
[`active_foil`](@ref)), the front element nearest each pixel of the foil face (y–z
plane; the lookup records the foil as `foil`). `nothing` otherwise.
"""
function flux_lookup(s::Scenario; foil::Union{Nothing,Symbol} = nothing, ntheta::Integer = 90,
                     nphi::Integer = 180, npix::Integer = 100)
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
    _has_part(t, foil) || return nothing
    es = [e for e in elements(t, foil) if t.elem_side[e] == 1]
    isempty(es) && return nothing
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
    return (; kind = :foil, foil, x = collect(ys), y = collect(zs), index = idx, elems = es,
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
# `ratio`: `nothing`, or the unconverted power fractions per beam, `(; omega1, omega2)`
# (vectors), deposited with `colors` (`UnconvertedColor`s: footprints, and absorption and
# x-ray conversion relative to 3ω).
#
# With unconverted light, the evaluation's deposition is the absorbed power of every
# color (what the laser coloring, power accounting and back-side losses show:
# `unconverted_deposition` has applied each color's `absorption`). The x-ray source is
# kept per color: the 3ω power at the material's conversion, each unconverted color's at
# the material's conversion times the color's `conversion`. The radiosity solve and the
# capsule flux take their sum as an equivalent 3ω power, P₃ + Σ c P_c (exact for a
# conversion independent of intensity; with an intensity-dependent conversion, it is
# evaluated at that equivalent power's intensity).
function _evaluate_drive(s::Scenario, beams, lookup, spot_radius::Real; ratio = nothing,
                         colors = nif_unconverted(), foil = _lookup_foil(lookup))
    t0 = time()
    ev = evaluate(s, beams)
    uc = nothing
    source = ev.deposition.element                   # the x-ray source's equivalent 3ω power
    if ratio !== nothing
        u = unconverted_deposition(s.target, beams; ratio, colors)
        d = ev.deposition
        tot = Deposition(d.element .+ u.omega1.element .+ u.omega2.element,
                         d.backside .+ u.omega1.backside .+ u.omega2.backside,
                         d.per_beam .+ u.omega1.per_beam .+ u.omega2.per_beam,
                         d.delivered + u.omega1.delivered + u.omega2.delivered)
        conv = Dict(c.name => c.conversion for c in colors)
        source = d.element .+ conv[:omega1] .* u.omega1.element .+ conv[:omega2] .* u.omega2.element
        flux = Briareus._radiative(s) ? Briareus._observed_flux(s, source) : nothing
        ev = Evaluation(s, ev.beams, ev.effective, tot, Ref{Union{Nothing,Vector{Float64}}}(flux))
        share(x, h) = (; absorbed = sum(x.per_beam), delivered = x.delivered, conversion = conv[h])
        uc = (; omega1 = share(u.omega1, :omega1), omega2 = share(u.omega2, :omega2), u.skipped)
    end
    sol = s.model === nothing ? nothing : solve(s.model, source)
    modes = NamedTuple[]
    mean = NaN
    if _has_flux(ev)
        c = flux_ylm(ev; lmax = maximum(MODE_DEGREES))
        a = flux_moments(ev; lmax = maximum(MODE_DEGREES))
        mean = a[1]
        modes = [(; l, epsilon = mode_amplitude(c, l), legendre = a[l+1] / a[1]) for l in MODE_DEGREES]
    end
    return DriveResult(ev, sol, modes, mean, power_rows(ev), _spot_stats(ev, spot_radius, foil),
                       _flux_map(lookup, ev, sol), time() - t0, uc)
end

"""
    unconverted_summary(res::DriveResult) -> String

The unconverted light's share of the drive in a line (absorbed and delivered 1ω and 2ω
power), `""` when none takes part.
"""
function unconverted_summary(res::DriveResult)
    u = res.unconverted
    u === nothing && return ""
    parts = String[]
    for (name, x) in (("1ω", u.omega1), ("2ω", u.omega2))
        x.delivered > 0 && push!(parts, @sprintf("%s %.4g TW absorbed of %.4g TW", name, x.absorbed, x.delivered))
    end
    isempty(u.skipped) || push!(parts, "$(length(u.skipped)) beams without a footprint left out")
    return _nobreak("unconverted light: " * join(parts, " · "))
end

"""
    power_rows(ev::Evaluation; minor = 1e-3) -> Vector{NamedTuple}

Laser power accounting rows `(; name, power, fraction)` (TW, fraction of delivered), as
`Briareus.report` prints them: each part's absorbed power and back-side (outside) loss,
in the target's order and named as in the model tree, with the parts below `minor` × the
delivered power summed on one row; then the residual (delivered less absorbed and
clipped), when not below that fraction, as `missed` (positive: light that strikes nothing) or `quadrature excess`
(negative: the deposition quadrature counts more power than the beams deliver).
"""
function power_rows(ev::Evaluation; minor::Real = 1e-3)
    t = ev.scenario.target
    d = ev.deposition
    P = d.delivered
    frac(v) = P > 0 ? v / P : 0.0
    cp = component_power(d, t)                  # by part: the drive's fitted pieces summed
    ps = [p for p in Briareus.parts(t) if t.components[first(part_pieces(t, p))].kind !== :aperture]
    parts = Tuple{String,Float64}[(component_label(p), get(cp, p, 0.0)) for p in ps]
    for p in ps
        k = Symbol(p, :_outside)
        haskey(cp, k) && push!(parts, (component_label(p) * " (outside)", cp[k]))
    end
    rows = NamedTuple[]
    small = 0
    rest = 0.0
    for (name, v) in parts
        if abs(v) < minor * P
            small += 1
            rest += v
        else
            push!(rows, (; name, power = v, fraction = frac(v)))
        end
    end
    small > 0 && push!(rows, (; name = string(small, small == 1 ? " part" : " parts", " < ",
                                              _pct_text(minor)), power = rest, fraction = frac(rest)))
    m = get(cp, :missed, 0.0)
    abs(m) < minor * P ||
        push!(rows, (; name = m > 0 ? "missed" : "quadrature excess", power = m, fraction = frac(m)))
    return rows
end

# A fraction as a percentage without trailing zeros: 1e-3 → "0.1 %".
_pct_text(x::Real) = _nobreak(string(round(100x; sigdigits = 3), " %"))

# Text with a no-break space before each unit and after each mode symbol (ε₂), so that a
# wrapped line keeps a number with its unit and its symbol.
_nobreak(s::AbstractString) = replace(s, r"(?<=ε[₀-₉]) | (?=%|W/cm²|eV|mm|TW)" => string(Char(0xa0)))

# Irradiance of foil `foil` (a part or component name) within `radius` of its center:
# mean, peak (W/cm²) and RMS nonuniformity; `nothing` without that foil.
function _spot_stats(ev::Evaluation, radius::Real, foil::Union{Nothing,Symbol})
    t = ev.scenario.target
    _has_part(t, foil) || return nothing
    front = [e for e in elements(t, foil) if t.elem_side[e] == 1]
    isempty(front) && return nothing
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
    drive_summary(res::DriveResult; short = false) -> String

The drive in a line: capsule mean flux, radiation temperature, and ε₂/ε₄; or the foil
spot irradiance and nonuniformity. `short`: the model tree's detail, the capsule's
radiation temperature and ε₂, or the foil spot's nonuniformity.
"""
function drive_summary(res::DriveResult; short::Bool = false)
    if short && (!isempty(res.modes) || res.spot !== nothing)
        isempty(res.modes) && return _nobreak(@sprintf("spot rms %.1f %%", 100 * res.spot.rms))
        e2 = only(r.epsilon for r in res.modes if r.l == 2)
        return _nobreak(@sprintf("Tr %.0f eV · ε₂ %.2f %%", radiation_temperature(res.mean_flux), 100 * e2))
    end
    if !isempty(res.modes)
        m = Dict(r.l => r.epsilon for r in res.modes)
        return _nobreak(@sprintf("capsule %s W/cm² (Tr %.0f eV) · ε₂ %.2f %% · ε₄ %.2f %%",
                                 _wcm2(res.mean_flux), radiation_temperature(res.mean_flux),
                                 100 * get(m, 2, NaN), 100 * get(m, 4, NaN)))
    elseif res.spot !== nothing
        s = res.spot
        return _nobreak(@sprintf("spot %.3g W/cm² · rms %.1f %% · %.0f %% of the power in r ≤ %.2g mm",
                                 s.mean, 100 * s.rms, 100 * s.fraction, s.radius))
    end
    d = res.evaluation.deposition
    return _nobreak(@sprintf("%.4g TW delivered", d.delivered)) *
           (res.solution === nothing ? " (laser only: no x-ray drive)" : "")
end

# ---------------------------------------------------------------- time

# The time panel's state (its logic is in timeline.jl, after the session).

"""
Charts of the time panel: `:pulses` (laser power per beam group), `:drive` (capsule
Tᵣ, Dante 1 and a_l/a_0), `:energy` (where the laser energy has gone by each time) and
`:retune` (the latest retune's reference and retuned powers and its mode 1).
"""
const TIME_CHARTS = (:pulses, :drive, :energy, :retune)

const TIME_CHART_LABELS = Dict(:pulses => "Pulses", :drive => "Drive history", :energy => "Energy",
                               :retune => "Retune")

"""
    DriveHistory

A drive history of the session ([`start_history!`](@ref)): `time` (ns), `Tr` (capsule
radiation temperature, eV), `legendre` (rows l = 1…4 of a_l/a_0 per step), `dante`
(Dante 1 brightness temperature, eV, `NaN` without one), `energy` (series `(; name,
values)` in kJ, cumulative), `fields` (per [`DRIVE_FIELDS`](@ref) an elements × steps
matrix in W/cm², on the drive's elements), `range` (each field's finite extremes over
the history), `history` (the core `drive_history`), `key` (what it was computed for)
and `seconds`.
"""
struct DriveHistory
    time::Vector{Float64}
    Tr::Vector{Float64}
    legendre::Matrix{Float64}
    dante::Vector{Float64}
    energy::Vector{NamedTuple}
    fields::Dict{Symbol,Matrix{Float32}}
    range::Dict{Symbol,Tuple{Float64,Float64}}
    history::Any
    key::Any
    seconds::Float64
end

"""
    TimeState()

State of the time panel: `pulses` (the pulses the beams fire, see
[`load_pulses!`](@ref) and [`use_pulse_shape!`](@ref); `nothing` before one), `history`
(the latest [`DriveHistory`](@ref)), `retune` (the latest retune's histories, see
[`retune_to_reference!`](@ref)), `chart` (one of [`TIME_CHARTS`](@ref)), `cursor` (the
time shown in the 3D view, ns; `NaN`: the drive at the table's powers), `dante` (record
a synthetic Dante 1), and the history job's `running`, `progress` (steps done, steps),
`started` and `status`.
"""
struct TimeState
    pulses::Observable{Any}
    history::Observable{Any}
    retune::Observable{Any}
    chart::Observable{Symbol}
    cursor::Observable{Float64}
    dante::Observable{Bool}
    running::Observable{Bool}
    progress::Observable{Tuple{Int,Int}}
    started::Observable{Float64}
    status::Observable{String}
    generation::Base.RefValue{Int}
    cancel::Threads.Atomic{Bool}
    queued::Base.RefValue{Bool}               # run once the drive is computed
end

TimeState() = TimeState(Observable{Any}(nothing), Observable{Any}(nothing), Observable{Any}(nothing),
                        Observable(:pulses), Observable(NaN), Observable(true), Observable(false),
                        Observable((0, 0)), Observable(0.0), Observable(""), Ref(0),
                        Threads.Atomic{Bool}(false), Ref(false))

"""
    set_chart!(tm::TimeState, chart::Symbol) -> tm

Show one of [`TIME_CHARTS`](@ref).
"""
function set_chart!(tm::TimeState, chart::Symbol)
    chart in TIME_CHARTS || throw(ArgumentError("chart must be one of $TIME_CHARTS"))
    _assign!(tm.chart, chart)
    return tm
end
