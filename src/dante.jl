# Synthetic Dante: x-ray power per steradian, its spectrum and the brightness
# temperature along a broadband spectrometer's line of sight through an entrance hole.

"""
    AbstractDanteChannel

A Dante channel's spectral response: a [`DanteChannel`](@ref) (Gaussian approximation)
or a [`DanteResponse`](@ref) (built from filters, mirror and XRD).
"""
abstract type AbstractDanteChannel end

"""
    DanteChannel(energy, fwhm)

One Dante channel's spectral response, approximated by a unit-peak Gaussian centred on
`energy` (eV) with full width `fwhm` (eV): the fallback of `dante1_channels(model =
:gaussian)`.
"""
struct DanteChannel <: AbstractDanteChannel
    energy::Float64
    fwhm::Float64
end

# NIF Dante 1 channel energies (eV), channels 1–18, as configured for ICF experiments
# (Rubery et al., "Soft x-ray power diagnostics for fusion experiments at NIF, Omega and
# Z facilities", RSI 2023, OSTI 1975922, Fig. 4).
const _DANTE1_ENERGIES = [1170.0, 3060.0, 250.0, 470.0, 710.0, 870.0, 940.0, 1170.0, 1430.0,
                          1710.0, 3430.0, 2580.0, 3000.0, 1700.0, 940.0, 14370.0, 690.0, 3660.0]

"""
    dante1_channels(; model = :hardware, resolution = 0.3) -> Vector

The 18 NIF Dante 1 channels (Rubery et al., RSI 2023, Fig. 4). `model = :hardware`
(default): [`dante1_responses`](@ref), each response built from the channel's filters,
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
              area = nothing)

A broadband spectrometer at port (θ, φ) (degrees) viewing a target through the aperture
component `through` (with its closing lip `<through>_lip`, if the target has one). The
visible surface is found by casting `pixels × pixels` parallel rays (the detector is
distant) across the hole, so every element counts with its visible projected area,
partial views past the rim included; a ray through a closing lip is weighted by the
lip's open fraction (`sol.transmission`). The brightness temperature is

    T = (π I / (σ A cos θ_v))^¼,

with `I` the power per steradian toward the detector, `θ_v` the angle between the line
of sight and the hole's normal and `A` the hole's full initial area (aperture plus lip;
Dante's convention of the pre-shot LEH size) unless `area` (mm²) is given. The footprint is
computed once per target (each keyframe of a moving target). `dante1()` and `dante2()`
give NIF's two.
"""
struct DanteView
    theta::Float64
    phi::Float64
    direction::Vec3
    through::Symbol
    pixels::Int
    channels::Vector{<:AbstractDanteChannel}
    area::Union{Nothing,Float64}
    cache::IdDict{Target,Any}
    lock::ReentrantLock
end

function DanteView(theta::Real, phi::Real; through::Symbol = :leh_bottom, pixels::Integer = 400,
                   channels = dante1_channels(), area = nothing)
    pixels >= 1 || throw(ArgumentError("pixels must be positive"))
    return DanteView(Float64(theta), Float64(phi), port_vector(theta, phi), through, Int(pixels),
                     collect(channels), area === nothing ? nothing : Float64(area),
                     IdDict{Target,Any}(), ReentrantLock())
end

"""
    dante1(; kwargs...)

NIF Dante 1: θ 143°, φ 274°, viewing the lower entrance hole at 37° to the hohlraum axis
(Rubery et al., RSI 2023; Rubery et al., PRL 132, 065104, 2024).
"""
dante1(; kwargs...) = DanteView(143, 274; through = :leh_bottom, kwargs...)

"""
    dante2(; kwargs...)

NIF Dante 2: θ 64°, φ 350°, a shallow view of the upper entrance hole (Rubery et al.,
RSI 2023).
"""
dante2(; kwargs...) = DanteView(64, 350; through = :leh_top, kwargs...)

"""
    sentinel(dim = :polar; through, kwargs...)

The NIF Sentinel filtered-XRD array (Kemp et al., RSI 2020, LLNL-JRNL-813698), a
DIM-mounted broadband spectrometer that can be fielded on the polar DIM (`:polar`,
000-000: θ 0°, straight down the hohlraum axis through the upper entrance hole) or on
either equatorial DIM (`:equatorial_78`, 090-078; `:equatorial_315`, 090-315), where a
hohlraum is seen only through a diagnostic window, so `through` must name it. Its own
filter set (16 Al, Ni and Cr MiniXRDs behind Dante filters) is not tabulated here: the
channels default to Dante 1's (`dante1_channels()`, *assumed*), to be replaced by the
shot's configuration.
"""
function sentinel(dim::Symbol = :polar; through = nothing, kwargs...)
    dim === :polar && return DanteView(0, 0; through = something(through, :leh_top), kwargs...)
    φ = dim === :equatorial_78 ? 78 : dim === :equatorial_315 ? 315 :
        throw(ArgumentError("dim must be :polar, :equatorial_78 or :equatorial_315"))
    through === nothing && throw(ArgumentError("an equatorial view needs `through` (a window)"))
    return DanteView(90, φ; through, kwargs...)
end

function Base.show(io::IO, v::DanteView)
    print(io, "DanteView(θ = ", v.theta, "°, φ = ", v.phi, "° through :", v.through, ")")
end

# Visible projected area (mm², normal to the line of sight) per (element, lip element
# the ray crossed or 0), the hole's area and the cosine of the viewing angle.
function _dante_footprint(t::Target, v::DanteView)
    names = [v.through]
    lip = Symbol(v.through, :_lip)
    any(c -> c.name === lip, t.components) && push!(names, lip)
    open = elements(t, v.through)
    hole = elements(t, names...)
    hole_area = sum(t.elem_area[hole])
    c = sum(t.elem_centroid[e] * t.elem_area[e] for e in open) / sum(t.elem_area[open])
    cosv = abs(dot(v.direction, t.elem_normal[first(hole)]))
    pts = [p for e in hole for k in t.elem_tris[e]
           for p in (t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k])]
    half = 1.05 * maximum(norm(p - c) for p in pts)
    d = v.direction
    r = normalize(cross(d, abs(d[3]) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0)))
    u = cross(r, d)
    far = 4 * (t.scale + norm(c)) + 1
    np = v.pixels
    px = 2half / np
    ray = -d
    # triangles of transmitting apertures (closing lips)
    lips = [k for (ci, comp) in enumerate(t.components)
            if comp.kind === :aperture && _transmits(comp.material.albedo)
            for e in t.comp_elems[ci] for k in t.elem_tris[e]]
    acc = Dict{Tuple{Int,Int},Float64}()
    for i in 1:np, j in 1:np
        o = c + (-half + (i - 0.5) * px) * r + (-half + (j - 0.5) * px) * u + far * d
        th, tri = closest_hit(t.bvh, o, ray)
        tri == 0 && continue
        e = Int(dot(ray, t.tri_normal[tri]) < 0 ? t.tri_front[tri] : t.tri_back[tri])
        e > 0 || continue
        l, tl = 0, th
        for k in lips
            s = _intersect(o, ray, t.tri_v0[k], t.tri_e1[k], t.tri_e2[k], tl)
            s < tl && ((tl, l) = (s, Int(max(t.tri_front[k], t.tri_back[k]))))
        end
        acc[(e, l)] = get(acc, (e, l), 0.0) + px^2
    end
    ks = collect(keys(acc))
    lip_elems = setdiff(hole, open)
    return (; elems = first.(ks), lips = last.(ks), area = [acc[k] for k in ks], hole_area, cosv,
            open_area = sum(t.elem_area[open]), lip_elems, lip_area = t.elem_area[lip_elems])
end

_footprint(v::DanteView, t::Target) = lock(() -> get!(() -> _dante_footprint(t, v), v.cache, t), v.lock)

# Power per steradian toward the detector from radiosity `B` (TW/mm²).
_los_intensity(fp, B, τ) =
    sum(fp.area[k] * ForwardDiff.value(B[fp.elems[k]]) / π * (fp.lips[k] == 0 ? 1.0 : τ[fp.lips[k]])
        for k in eachindex(fp.elems); init = 0.0)

"""
    visible_area(view::DanteView, target) -> Dict{Symbol,Float64}

Projected area (mm², normal to the line of sight) of each part that the view sees
through its hole (closing lips counted fully open); the pieces of a part split for
computation are summed into it.
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

Synthetic Dante signal: `intensity` (TW/sr toward the detector), `Tr` (eV, the
brightness temperature over the hole; see `DanteView`), `area` (mm²) and `cos` of the
viewing angle; `Tr_open`, the brightness temperature over the hole's open area (the
aperture plus each closing-lip element's area times its open fraction), the radiation
temperature behind the hole that the view samples (equal to `Tr` without lips). For a
`SpectralSolution` also `spectrum` (TW/sr per group), `edges`, `channels` (each
channel's signal; the open-ended top group contributes nothing), `color_temperature`
(eV, the Planckian whose channel signals best match, in log, those of the channels
below 1.5 keV, as Dante's unfold starts) and `mband_fraction` (of the intensity above
1.8 keV, `NaN` without that edge).

With [`DanteResponse`](@ref) channels (the default, `dante1_channels()`) `channels` are
the voltages `Ω ∫ (dI/dE) R dE` (V), the spectrum taken within each group as a
Planckian at `Tr_open` (the temperature behind the hole), and the voltages are unfolded
as data are ([`dante_unfold`](@ref) with `method = unfold` and `kwargs`, over the
view's area and angle): `unfold` (its result) and `Tr_unfold`, the unfolded radiation
temperature, to compare with `Tr`. `unfold = nothing` skips it. With Gaussian
[`DanteChannel`](@ref)s `channels` is ∫ (dI/dE) R(E) dE (TW/sr) with the spectrum flat
within each group, and there is no unfold.
"""
function dante(v::DanteView, sol::RadiationSolution)
    fp = _footprint(v, sol.model.target)
    I = _los_intensity(fp, sol.B, sol.transmission)
    A = something(v.area, fp.hole_area)
    return (; intensity = I, Tr = radiation_temperature(π * I / (A * fp.cosv)), area = A,
            cos = fp.cosv, Tr_open = _open_temperature(fp, I, sol.transmission))
end

function _open_temperature(fp, I, τ)
    A = fp.open_area + sum(fp.lip_area[k] * τ[e] for (k, e) in enumerate(fp.lip_elems); init = 0.0)
    return radiation_temperature(π * I / (A * fp.cosv))
end

function dante(v::DanteView, s::SpectralSolution; unfold = :traditional, kwargs...)
    sol = s.gray
    fp = _footprint(v, sol.model.target)
    G = ngroups(s)
    spectrum = [_los_intensity(fp, view(s.B, :, g), sol.transmission) for g in 1:G]
    I = sum(spectrum)
    A = something(v.area, fp.hole_area)
    Tr_open = _open_temperature(fp, I, sol.transmission)
    hw = eltype(v.channels) <: DanteResponse
    W = hw ? _channel_matrix(v.channels, s.edges; T = Tr_open) : _channel_matrix(v.channels, s.edges)
    signals = hw ? W * (1e3 .* spectrum) : W * spectrum
    k = findfirst(==(1800.0), s.edges)
    mband_fraction = k === nothing ? NaN : sum(spectrum[k:end]) / I
    base = (; intensity = I, Tr = radiation_temperature(π * I / (A * fp.cosv)), area = A,
            cos = fp.cosv, Tr_open, spectrum, edges = s.edges, channels = signals,
            color_temperature = _color_temperature(v.channels, s.edges, W, signals),
            mband_fraction)
    (hw && unfold !== nothing) || return base
    u = dante_unfold(v.channels, signals; area = A, cos = fp.cosv, method = unfold, kwargs...)
    return (; base..., unfold = u, Tr_unfold = u.Tr)
end

# W[c, g]: channel c's signal per unit power in group g, spread flat over the group.
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
