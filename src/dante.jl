# Synthetic Dante: x-ray power per steradian, its spectrum and the brightness
# temperature along a broadband spectrometer's line of sight through an entrance hole.

"""
    DanteChannel(energy, fwhm)

One Dante channel's spectral response, approximated by a unit-peak Gaussian centred on
`energy` (eV) with full width `fwhm` (eV).
"""
struct DanteChannel
    energy::Float64
    fwhm::Float64
end

# NIF Dante 1 channel energies (eV), channels 1–18, as configured for ICF experiments
# (Rubery et al., "Soft x-ray power diagnostics for fusion experiments at NIF, Omega and
# Z facilities", RSI 2023, OSTI 1975922, Fig. 4).
const _DANTE1_ENERGIES = [1170.0, 3060.0, 250.0, 470.0, 710.0, 870.0, 940.0, 1170.0, 1430.0,
                          1710.0, 3430.0, 2580.0, 3000.0, 1700.0, 940.0, 14370.0, 690.0, 3660.0]

"""
    dante1_channels(; resolution = 0.3) -> Vector{DanteChannel}

The 18 NIF Dante 1 channels at their published band energies (Rubery et al., RSI 2023,
Fig. 4: 0.25–3.66 keV and one at 14.4 keV), each a Gaussian of FWHM `resolution × energy`.
An **approximation**: the real responses combine filter transmission (with K- and
L-edges), grazing-incidence mirror cutoffs and x-ray diode sensitivity, and have
high-energy tails; Dante's own unfold approximates each response by a Gaussian, whose
widths are not published there, so `resolution` is an engineering choice.
"""
dante1_channels(; resolution::Real = 0.3) =
    [DanteChannel(E, resolution * E) for E in _DANTE1_ENERGIES]

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
    channels::Vector{DanteChannel}
    area::Union{Nothing,Float64}
    cache::IdDict{Target,Any}
    lock::ReentrantLock
end

function DanteView(theta::Real, phi::Real; through::Symbol = :leh_bottom, pixels::Integer = 400,
                   channels = dante1_channels(), area = nothing)
    pixels >= 1 || throw(ArgumentError("pixels must be positive"))
    return DanteView(Float64(theta), Float64(phi), port_vector(theta, phi), through, Int(pixels),
                     collect(DanteChannel, channels), area === nothing ? nothing : Float64(area),
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

function Base.show(io::IO, v::DanteView)
    print(io, "DanteView(θ = ", v.theta, "°, φ = ", v.phi, "° through :", v.through, ")")
end

# Visible projected area (mm², normal to the line of sight) per (element, lip element
# the ray crossed or 0), the hole's area and the cosine of the viewing angle.
function _dante_footprint(t::Target, v::DanteView)
    names = [v.through]
    lip = Symbol(v.through, :_lip)
    any(c -> c.name === lip, t.components) && push!(names, lip)
    open = t.comp_elems[component_index(t, v.through)]
    hole = reduce(vcat, [t.comp_elems[component_index(t, nm)] for nm in names])
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
    return (; elems = first.(ks), lips = last.(ks), area = [acc[k] for k in ks], hole_area, cosv)
end

_footprint(v::DanteView, t::Target) = lock(() -> get!(() -> _dante_footprint(t, v), v.cache, t), v.lock)

# Power per steradian toward the detector from radiosity `B` (TW/mm²).
_los_intensity(fp, B, τ) =
    sum(fp.area[k] * ForwardDiff.value(B[fp.elems[k]]) / π * (fp.lips[k] == 0 ? 1.0 : τ[fp.lips[k]])
        for k in eachindex(fp.elems); init = 0.0)

"""
    visible_area(view::DanteView, target) -> Dict{Symbol,Float64}

Projected area (mm², normal to the line of sight) of each component that the view sees
through its hole (closing lips counted fully open).
"""
function visible_area(v::DanteView, t::Target)
    fp = _footprint(v, t)
    out = Dict{Symbol,Float64}()
    for (e, a) in zip(fp.elems, fp.area)
        n = t.components[t.elem_comp[e]].name
        out[n] = get(out, n, 0.0) + a
    end
    return out
end

"""
    dante(view::DanteView, sol::RadiationSolution) -> NamedTuple
    dante(view::DanteView, s::SpectralSolution) -> NamedTuple

Synthetic Dante signal: `intensity` (TW/sr toward the detector), `Tr` (eV, the
brightness temperature over the hole; see `DanteView`), `area` (mm²) and `cos` of the
viewing angle. For a `SpectralSolution` also `spectrum` (TW/sr per group), `edges`,
`channels` (each channel's signal, ∫ (dI/dE) R(E) dE in TW/sr, the spectrum taken flat
within each group; the open-ended top group contributes nothing), `color_temperature`
(eV, the Planckian whose channel signals best match, in log, those of the channels
below 1.5 keV, as Dante's unfold starts) and `mband_fraction` (of the intensity above
1.8 keV, `NaN` without that edge).
"""
function dante(v::DanteView, sol::RadiationSolution)
    fp = _footprint(v, sol.model.target)
    I = _los_intensity(fp, sol.B, sol.transmission)
    A = something(v.area, fp.hole_area)
    return (; intensity = I, Tr = radiation_temperature(π * I / (A * fp.cosv)), area = A,
            cos = fp.cosv)
end

function dante(v::DanteView, s::SpectralSolution)
    sol = s.gray
    fp = _footprint(v, sol.model.target)
    G = ngroups(s)
    spectrum = [_los_intensity(fp, view(s.B, :, g), sol.transmission) for g in 1:G]
    I = sum(spectrum)
    A = something(v.area, fp.hole_area)
    W = _channel_matrix(v.channels, s.edges)
    signals = W * spectrum
    k = findfirst(==(1800.0), s.edges)
    mband_fraction = k === nothing ? NaN : sum(spectrum[k:end]) / I
    return (; intensity = I, Tr = radiation_temperature(π * I / (A * fp.cosv)), area = A,
            cos = fp.cosv, spectrum, edges = s.edges, channels = signals,
            color_temperature = _color_temperature(v.channels, s.edges, W, signals),
            mband_fraction)
end

# W[c, g]: channel c's signal per unit power in group g, spread flat over the group.
function _channel_matrix(ch::Vector{DanteChannel}, E::Vector{Float64})
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
