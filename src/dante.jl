# This file computes a synthetic Dante signal: the x-ray power per steradian, its
# spectrum and the brightness temperature along the line of sight of a broadband
# spectrometer through an entrance hole.

"""
    AbstractDanteChannel

The spectral response of a Dante channel: a [`DanteChannel`](@ref) (Gaussian approximation)
or a [`DanteResponse`](@ref) (built from filters, mirror and XRD).
"""
abstract type AbstractDanteChannel end

"""
    DanteChannel(energy, fwhm)

The spectral response of one Dante channel, approximated by a unit-peak Gaussian centered on
`energy` (eV) with full width `fwhm` (eV): the fallback of `dante1_channels(model =
:gaussian)`.
"""
struct DanteChannel <: AbstractDanteChannel
    energy::Float64
    fwhm::Float64
end

# These are the NIF Dante 1 channel energies (eV), channels 1–18, as configured for ICF
# experiments (Rubery et al., "Soft x-ray power diagnostics for fusion experiments at NIF,
# Omega and Z facilities", RSI 2023, OSTI 1975922, Fig. 4).
const _DANTE1_ENERGIES = [1170.0, 3060.0, 250.0, 470.0, 710.0, 870.0, 940.0, 1170.0, 1430.0,
                          1710.0, 3430.0, 2580.0, 3000.0, 1700.0, 940.0, 14370.0, 690.0, 3660.0]

"""
    dante1_channels(; model = :hardware, resolution = 0.3) -> Vector

The 18 NIF Dante 1 channels (Rubery et al., RSI 2023, Fig. 4). `model = :hardware`
(default): [`dante1_responses`](@ref), each response built from the channel filters,
grid, mirror and XRD cathode. `model = :gaussian`: the fallback approximation, a
[`DanteChannel`](@ref) Gaussian at each published band energy (0.25–3.66 keV and one at
14.4 keV) of FWHM `resolution × energy`; it has no K- or L-edge structure, mirror
cutoffs or high-energy tails, and `resolution` is an engineering choice.
"""
function dante1_channels(; model::Symbol = :hardware, resolution::Real = 0.3)
    model === :hardware && return dante1_responses()
    model === :gaussian || throw(ArgumentError("model must be :hardware or :gaussian"))
    return [DanteChannel(E, resolution * E) for E in _DANTE1_ENERGIES]
end

"""
    DanteView(theta, phi; through = :leh_bottom, pixels = 400, channels = dante1_channels(),
              area = nothing, field_of_view = nothing, center = (0, 0, 0))

A broadband spectrometer at port (θ, φ) (degrees). The visible surface is found by
casting `pixels × pixels` parallel rays toward the target, since the detector is distant.
Every element then counts with its visible projected area, and partial views and
occlusion by any part are included. A ray through a transmitting aperture, such as a
closing lip, is weighted by its open fraction (`sol.transmission`). The rays cover one of
two regions:

- `field_of_view = nothing` (default) gives the view through one entrance hole. The rays
  cover a square just enclosing the aperture component `through` and its closing lip
  `<through>_lip`, if the target has one.
- `field_of_view = D` (mm) gives the collimated view. The rays cover a disk of diameter
  `D` normal to the line of sight about `center`, cut to the target extent, with
  `pixels` rays across it. `center` is a point (mm) or a part name, which stands for the
  area-weighted centroid of the part. The field of view of a NIF Dante is ±30 mm at TCC
  in the 2004 design (research notes, "Dante-family hardware and geometry"). Every
  surface of a target that faces the detector therefore lies inside it, such as walls
  seen through windows, packages and shields.

A channel with a separate `field_of_view` ([`DanteResponse`](@ref)) records the signal over
that field in place of the view field. That field is a disk of its diameter about `center`;
in a view through a hole, it is the square about the hole cut to that disk. The signal of
the channel comes from the footprint of that field, which is cast once per distinct field
and target. Channels that share the field of the view therefore cost nothing more and record
exactly what they record without a separate field. The intensity, `Tr` and the visible areas
belong to the view.

The brightness temperature is

    T = (π I / (σ A cos θ_v))^¼,

with `I` the power per steradian toward the detector and `A cos θ_v` the reference area
seen from the detector. By default `A` is the full initial area of `through`, aperture
plus lip, as in the Dante convention of the pre-shot LEH size. `θ_v` is the angle between
the line of sight and the normal of `through`. The keyword `area` sets `A` otherwise:

- `area = A` (mm²) replaces `A`, with the aperture `θ_v`, or `cos θ_v = 1` without
  `through`.
- `area = :part` takes the visible projected area of the part as `A cos θ_v`. This is the
  convention of Fournier et al. 2010, who used the projected wall and cap.
- `area = :open` takes the open area of the hole at the solution time. That area is the
  aperture plus the area of each closing-lip element times its open fraction. It comes
  from the entrance-hole closure of the model, so `Tr` then equals `Tr_open`.
- `area = TabulatedCurve(t, A)` takes a measured area history, with time in ns, area in
  mm² and the end values held outside. Glenzer et al. 2010 used the hole area from SXI
  images.

With a field of view, `through = nothing` and no `area` give `Tr = NaN`, though the power
is still computed. The footprint is computed once per target, or once per keyframe of a
moving target. `dante1()` and `dante2()` give the two NIF Dantes.
"""
struct DanteView
    theta::Float64
    phi::Float64
    direction::Vec3
    through::Union{Nothing,Symbol}
    pixels::Int
    channels::Vector{<:AbstractDanteChannel}
    area::Union{Nothing,Float64,Symbol,TabulatedCurve}
    field_of_view::Union{Nothing,Float64}
    center::Union{Vec3,Symbol}
    cache::WeakKeyDict{Target,Any}
    lock::ReentrantLock
end

function DanteView(theta::Real, phi::Real; through::Union{Nothing,Symbol} = :leh_bottom,
                   pixels::Integer = 400, channels = dante1_channels(), area = nothing,
                   field_of_view = nothing, center = zero(Vec3))
    pixels >= 1 || throw(ArgumentError("pixels must be positive"))
    field_of_view === nothing && through === nothing &&
        throw(ArgumentError("a view without a field of view looks through an aperture (`through`)"))
    field_of_view === nothing || field_of_view > 0 || throw(ArgumentError("field_of_view must be positive"))
    a = area === nothing || area isa Union{Symbol,TabulatedCurve} ? area : Float64(area)
    a === :open && through === nothing && throw(ArgumentError("area = :open is the open area of `through`"))
    c = center isa Symbol ? center : Vec3(center...)
    return DanteView(Float64(theta), Float64(phi), port_vector(theta, phi), through, Int(pixels),
                     collect(channels), a, field_of_view === nothing ? nothing : Float64(field_of_view),
                     c, WeakKeyDict{Target,Any}(), ReentrantLock())
end

"""
    sentinel(dim = :polar; through, kwargs...)

Returns the NIF Sentinel filtered-XRD array (Kemp et al., RSI 2020, LLNL-JRNL-813698), a
DIM-mounted broadband spectrometer. It can be fielded on the polar DIM or on either
equatorial DIM (NIF User Guide 2024-09, Table 5-1):

- `:polar` is port 000-000 at θ 0°, on the hohlraum axis, with the view straight down
  through the upper entrance hole.
- `:equatorial_78` is P90-78 at φ 78.75°, and `:equatorial_315` is P90-315. From either
  port a hohlraum is seen only through a diagnostic window, so `through` must name it.

The Sentinel filter set of 16 Al, Ni and Cr MiniXRDs behind Dante filters is not
tabulated in Briareus. The channels default to those of Dante 1 (`dante1_channels()`,
*assumed*) without their per-channel fields of view, since the Sentinel fields are not
public. The channels of the shot configuration should replace them.
"""
function sentinel(dim::Symbol = :polar; through = nothing, channels = _with_field.(dante1_channels(), nothing),
                  kwargs...)
    dim === :polar && return DanteView(0, 0; through = something(through, :leh_top), channels, kwargs...)
    φ = dim === :equatorial_78 ? 78.75 : dim === :equatorial_315 ? 315 :
        throw(ArgumentError("dim must be :polar, :equatorial_78 or :equatorial_315"))
    through === nothing && throw(ArgumentError("an equatorial view needs `through` (a window)"))
    return DanteView(90, φ; through, channels, kwargs...)
end

function Base.show(io::IO, v::DanteView)
    print(io, "DanteView(θ = ", v.theta, "°, φ = ", v.phi, "°",
          v.field_of_view === nothing ? "" : ", field of view $(v.field_of_view) mm",
          v.through === nothing ? "" : " through :$(v.through)", ")")
end

# Returns the elements, area and centroid of the reference aperture `through` (with its
# lip) and the cosine of the viewing angle. It returns `nothing` when the view has no
# aperture (a field-of-view view only) or the target lacks it.
function _dante_hole(t::Target, v::DanteView)
    v.through === nothing && return nothing
    if !any(c -> c.name === v.through, t.components)
        v.field_of_view === nothing && throw(ArgumentError("the target has no component :$(v.through)"))
        return nothing
    end
    names = [v.through]
    lip = Symbol(v.through, :_lip)
    any(c -> c.name === lip, t.components) && push!(names, lip)
    open = elements(t, v.through)
    hole = elements(t, names...)
    c = sum(t.elem_centroid[e] * t.elem_area[e] for e in open) / sum(t.elem_area[open])
    return (; open, hole, area = sum(t.elem_area[hole]), centroid = c,
            cosv = abs(dot(v.direction, t.elem_normal[first(hole)])))
end

# Casts the rays of a view over the field `fov` (mm), which is the view field or the
# field of a channel. The rays cover a disk about the view center. When the view has no
# field of view, they cover the square about the hole instead, cut to the disk of `fov`
# when a channel gives one. The result holds:
# - `elems`, `lips` and `area`, the visible projected area (mm², normal to the line of
#   sight) per element and the lip element that its ray crossed, 0 for none;
# - `hole_area`, `cosv` and `open_area`, the area of the reference hole, the cosine of
#   the viewing angle and the open area of the hole;
# - `lip_elems` and `lip_area`, the lip elements and their areas;
# - `ref_area` and `ref_cos`, the reference area and its cosine (`DanteView`).
function _dante_footprint(t::Target, v::DanteView, fov = v.field_of_view)
    h = _dante_hole(t, v)
    d = v.direction
    cf = fov === nothing ? zero(Vec3) : v.center isa Symbol ? part_centroid(t, v.center) : v.center
    # In a view through a hole, only the hole rays within that disk count for a channel
    # field.
    Rc2 = v.field_of_view === nothing && fov !== nothing ? (fov / 2)^2 : Inf
    if v.field_of_view === nothing
        c = h.centroid
        pts = [p for e in h.hole for k in t.elem_tris[e]
               for p in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])]
        half = 1.05 * maximum(norm(p - c) for p in pts)
        R2 = Inf
    else
        c = cf
        # Measure the target extent across the line of sight about c.
        ext = 0.0
        for k in eachindex(t.tri_v0), p in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])
            q = p - c
            ext = max(ext, norm(q - dot(q, d) * d))
        end
        R = fov / 2
        half = min(R, 1.0001 * ext)
        R2 = R^2
    end
    r = normalize(cross(d, abs(d[3]) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0)))
    u = cross(r, d)
    far = 4 * (t.scale + norm(c)) + 1
    np = v.pixels
    lips = _lip_bvh(t)      # of transmitting apertures, which are closing lips
    acc = _footprint_rays(t, lips, c, cf, r, u, d, Float64(half), Float64(R2), Float64(Rc2), far, np)
    ks = collect(keys(acc))
    elems, area = first.(ks), [acc[k] for k in ks]
    if h === nothing
        hole_area = cosv = open_area = NaN
        lip_elems = Int[]
    else
        hole_area, cosv, open_area = h.area, h.cosv, sum(t.elem_area[h.open])
        lip_elems = setdiff(h.hole, h.open)
    end
    # Take the reference area A cos θ_v seen from the detector (`DanteView`). `:open` and
    # an area history depend on the solution, so `_reference` handles them.
    A, cref = if v.area isa Symbol && v.area !== :open
        cis = Set(part_pieces(t, v.area))
        sum((a for (e, a) in zip(elems, area) if t.elem_comp[e] in cis); init = 0.0), 1.0
    elseif v.area isa Float64
        v.area, h === nothing ? 1.0 : cosv
    else
        hole_area, cosv
    end
    return (; elems, lips = last.(ks), area, hole_area, cosv, open_area, lip_elems,
            lip_area = t.elem_area[lip_elems], ref_area = A, ref_cos = cref)
end

# Casts the np² rays of a footprint toward the detector and returns the area per pair of
# first element and first lip element (0 for none). The loop is behind this barrier since
# the center and the extents of the caller depend on the kind of view.
function _footprint_rays(t::Target, lips::BVH, c::Vec3, cf::Vec3, r::Vec3, u::Vec3, d::Vec3,
                         half::Float64, R2::Float64, Rc2::Float64, far::Float64, np::Int)
    px = 2half / np
    ray = -d
    acc = Dict{Tuple{Int,Int},Float64}()
    for i in 1:np, j in 1:np
        x, y = -half + (i - 0.5) * px, -half + (j - 0.5) * px
        x^2 + y^2 > R2 && continue
        if Rc2 < Inf
            q = c + x * r + y * u - cf
            sum(abs2, q - dot(q, d) * d) > Rc2 && continue
        end
        o = c + x * r + y * u + far * d
        th, tri = closest_hit(t.bvh, o, ray)
        tri == 0 && continue
        e = Int(dot(ray, t.tri_normal[tri]) < 0 ? t.tri_front[tri] : t.tri_back[tri])
        e > 0 || continue
        l = _first_lip(t, lips, o, ray, th)
        acc[(e, l)] = get(acc, (e, l), 0.0) + px^2
    end
    return acc
end

# Returns the footprint over the field `fov`, the view field by default. Each footprint
# is cast once per target and field and then cached.
_footprint(v::DanteView, t::Target, fov = v.field_of_view) =
    lock(v.lock) do
        byfield = get!(() -> Dict{Union{Nothing,Float64},Any}(), v.cache, t)
        get!(() -> _dante_footprint(t, v, fov), byfield, fov)
    end

# Returns the field that a channel records, which is the channel field if it has one and
# the view field otherwise.
_channel_field(v::DanteView, c) =
    c isa DanteResponse && c.field_of_view !== nothing ? c.field_of_view : v.field_of_view

# Returns the power per steradian toward the detector from the radiosity `B` (TW/mm²).
_los_intensity(fp, B, τ) =
    sum(fp.area[k] * ForwardDiff.value(B[fp.elems[k]]) / π * (fp.lips[k] == 0 ? 1.0 : τ[fp.lips[k]])
        for k in eachindex(fp.elems); init = 0.0)

"""
    visible_area(view::DanteView, target) -> Dict{Symbol,Float64}

Projected area (mm², normal to the line of sight) of each part that the view sees,
through its hole or within its field of view (closing lips counted fully open); the
pieces of a part split for computation are summed into it.
"""
function visible_area(v::DanteView, t::Target)
    fp = _footprint(v, t)
    out = Dict{Symbol,Float64}()
    for (e, a) in zip(fp.elems, fp.area)
        n = t.components[t.elem_comp[e]].part
        out[n] = get(out, n, 0.0) + a
    end
    return out
end

"""
    dante(view::DanteView, sol::RadiationSolution) -> NamedTuple
    dante(view::DanteView, s::SpectralSolution; unfold = :traditional, kwargs...) -> NamedTuple

Computes the synthetic Dante signal. The result holds:

- `intensity` (TW/sr toward the detector);
- `Tr` (eV), the brightness temperature over the reference area (see `DanteView`);
- `area` (mm²) and `cos` of the reference area, so that `area × cos` is its projection
  toward the detector;
- `Tr_open`, the brightness temperature over the open area of the hole, which is the
  aperture plus the area of each closing-lip element times its open fraction. It is the
  radiation temperature behind the hole that the view samples. It equals `Tr` without
  lips and is `NaN` without `through`.

Without a reference area `Tr` is `NaN`, and the in-group spectrum and the unfold take the
visible projected area instead, which gives a brightness temperature over what is seen.
For a `SpectralSolution` the result also holds:

- `spectrum` (TW/sr per group) and `edges`;
- `channels`, the signal of each channel over its field of view where it has one, with
  nothing from the open-ended top group;
- `color_temperature` (eV), the Planckian whose channel signals best match, in log, those
  of the channels below 1.5 keV, as the Dante unfold starts;
- `mband_fraction`, the fraction of the intensity above 1.8 keV, `NaN` without that edge.

With [`DanteResponse`](@ref) channels, the default of `dante1_channels()`, `channels` are
the voltages `Ω ∫ (dI/dE) R dE` (V). The spectrum within each group is taken as a
Planckian at `Tr_open`, the temperature behind the hole. The voltages are then unfolded
as data are, by [`dante_unfold`](@ref) with `method = unfold` and `kwargs` over the area
and angle of the view. The result then also holds `unfold`, the result of the unfold,
and `Tr_unfold`, the unfolded radiation temperature to compare with `Tr`.
`unfold = nothing` skips the unfold. With Gaussian [`DanteChannel`](@ref)s, `channels`
is ∫ (dI/dE) R(E) dE (TW/sr) with the spectrum flat within each group, and there is no
unfold.
"""
function dante(v::DanteView, sol::RadiationSolution)
    fp = _footprint(v, sol.model.target)
    I = _los_intensity(fp, sol.B, sol.transmission)
    A, c = _reference(v, fp, sol.transmission, sol.time)
    return (; intensity = I, Tr = radiation_temperature(π * I / (A * c)), area = A,
            cos = c, Tr_open = _open_temperature(fp, I, sol.transmission))
end

# Returns the open area of the hole, the aperture plus each closing-lip element times its
# open fraction.
_open_area(fp, τ) = fp.open_area + sum(fp.lip_area[k] * τ[e] for (k, e) in enumerate(fp.lip_elems); init = 0.0)

_open_temperature(fp, I, τ) = radiation_temperature(π * I / (_open_area(fp, τ) * fp.cosv))

# Returns the reference area A and cos θ_v for a solution with lip transmissions τ at
# `time` (ns).
function _reference(v::DanteView, fp, τ, time)
    v.area === :open && return (_open_area(fp, τ), fp.cosv)
    v.area isa TabulatedCurve && return (v.area(time), isnan(fp.cosv) ? 1.0 : fp.cosv)
    return (fp.ref_area, fp.ref_cos)
end

function dante(v::DanteView, s::SpectralSolution; unfold = :traditional, kwargs...)
    sol = s.gray
    fp = _footprint(v, sol.model.target)
    G = ngroups(s)
    spectrum = [_los_intensity(fp, view(s.B, :, g), sol.transmission) for g in 1:G]
    I = sum(spectrum)
    A, c = _reference(v, fp, sol.transmission, sol.time)
    Tr_open = _open_temperature(fp, I, sol.transmission)
    # Without a reference area, the in-group spectrum and the unfold use the visible
    # projected area, which gives a brightness temperature over what is seen.
    seen = sum(fp.area; init = 0.0)
    Tshape = isfinite(Tr_open) ? Tr_open : radiation_temperature(π * I / seen)
    Au, cu = isfinite(A * c) ? (A, c) : (seen, 1.0)
    hw = eltype(v.channels) <: DanteResponse
    W = hw ? _channel_matrix(v.channels, s.edges; T = Tshape) : _channel_matrix(v.channels, s.edges)
    signals = hw ? W * (1e3 .* spectrum) : W * spectrum
    # For channels with a separate field, compute the footprint and spectrum of each
    # distinct field.
    fields = [_channel_field(v, c) for c in v.channels]
    for f in unique(fields)
        isequal(f, v.field_of_view) && continue
        fpf = _footprint(v, sol.model.target, f)
        sf = [_los_intensity(fpf, view(s.B, :, g), sol.transmission) for g in 1:G]
        sig = hw ? W * (1e3 .* sf) : W * sf
        for c in eachindex(signals)
            isequal(fields[c], f) && (signals[c] = sig[c])
        end
    end
    k = findfirst(==(1800.0), s.edges)
    mband_fraction = k === nothing ? NaN : sum(spectrum[k:end]) / I
    base = (; intensity = I, Tr = radiation_temperature(π * I / (A * c)), area = A,
            cos = c, Tr_open, spectrum, edges = s.edges, channels = signals,
            color_temperature = _color_temperature(v.channels, s.edges, W, signals),
            mband_fraction)
    (hw && unfold !== nothing) || return base
    u = dante_unfold(v.channels, signals; area = Au, cos = cu, method = unfold, kwargs...)
    return (; base..., unfold = u, Tr_unfold = u.Tr)
end

# Fills W[c, g] with the signal of channel c per unit power in group g, with the power
# spread flat over the group.
function _channel_matrix(ch::Vector{DanteChannel}, E::Vector{Float64}; T = nothing)
    G = length(E) - 1
    W = zeros(length(ch), G)
    for (c, h) in enumerate(ch), g in 1:G
        isfinite(E[g+1]) || continue
        σ = h.fwhm / (2 * sqrt(2 * log(2)))
        a, b = (E[g] - h.energy) / (sqrt(2) * σ), (E[g+1] - h.energy) / (sqrt(2) * σ)
        W[c, g] = σ * sqrt(π / 2) * (erf(b) - erf(a)) / (E[g+1] - E[g])
    end
    return W
end

function _color_temperature(ch, E, W, signals; below = 1500.0)
    use = [c for c in eachindex(ch) if ch[c].energy < below && signals[c] > 0]
    length(use) >= 2 || return NaN
    ls = log.(signals[use])
    function misfit(lT)
        T = exp(lT)
        model = W[use, :] * planck_fractions(T, E)
        any(<=(0), model) && return Inf
        r = ls .- log.(model)
        return sum(abs2, r .- sum(r) / length(r))
    end
    grid = range(log(10.0), log(3000.0), 241)
    k = argmin(misfit.(grid))
    lo, hi = grid[max(k - 1, 1)], grid[min(k + 1, end)]
    φ = (sqrt(5) - 1) / 2
    for _ in 1:60                                   # golden section
        a, b = hi - φ * (hi - lo), lo + φ * (hi - lo)
        misfit(a) < misfit(b) ? (hi = b) : (lo = a)
    end
    return exp((lo + hi) / 2)
end
