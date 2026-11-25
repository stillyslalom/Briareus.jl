# This file builds Dante channel responses from their parts, which are filters, grids,
# grazing-incidence mirrors and x-ray diode cathodes. It also holds the Dante unfold,
# which turns channel voltages back into a spectrum, a flux and a radiation temperature.
# See the research notes, "Dante channel responses and the unfold (2026-10-24)".

# The energy grid of the responses and of the unfold spans 50 eV – 30 keV in 0.21% steps.
const _DANTE_GRID = exp.(range(log(50.0), log(29_990.0), 3001))

# The digitized curves of Rubery et al., RSI 2023 (OSTI 1975922), Figs. 3 and 8, are read
# from `data/rubery2023_dante.csv` as curve name => (energy eV, value).
const _RUBERY_DANTE = Ref{Union{Nothing,Dict{String,Tuple{Vector{Float64},Vector{Float64}}}}}(nothing)
const _RUBERY_LOCK = ReentrantLock()

function _rubery_dante()
    lock(_RUBERY_LOCK) do
        _RUBERY_DANTE[] === nothing || return _RUBERY_DANTE[]
        h, rows = read_csv(joinpath(@__DIR__, "..", "data", "rubery2023_dante.csv"))
        d = Dict{String,Tuple{Vector{Float64},Vector{Float64}}}()
        for r in rows
            E, v = get!(() -> (Float64[], Float64[]), d, r[1])
            push!(E, parse(Float64, r[2])); push!(v, parse(Float64, r[3]))
        end
        _RUBERY_DANTE[] = d
    end
end

# ---- x-ray diode cathodes ---------------------------------------------------------------
#
# The sensitivity in V/GW into 50 Ω is S(E) = K (E/1 keV)^γ x/(1 + x), with
# x = L ρ μ_s(E). This is Henke's total-electron-yield form, with a surface layer and a
# power-law correction E^γ (Henke, Knauer & Premaratne, J. Appl. Phys. 52, 1509, 1981).
# The yield is proportional to the energy absorbed within an escape depth L, so S ∝ μ
# while μL ≪ 1. The surface layer has mass fractions w_C of carbon and w_O of oxygen over
# the bulk metal; μ_s is the μ/ρ of the mixture, and ρ is the density of the metal.
#
# The parameters are fitted in log (`examples/dante_channels.jl`) to the measured curves
# digitized from Rubery 2023:
# - For Cr, L, w_C, w_O, γ and K are fitted to Fig. 3b over 0.075–1 keV. The rms is 12%
#   against the median-read curve, and the fit does not reproduce the edge spikes.
# - For Al, K and γ are fitted to Fig. 8b over 1–4.5 keV, across the Al K edge, with the
#   surface layer of Cr. The rms is 7%.
# - Ni has no published curve and takes the Cr parameters on the Ni bulk (*assumed*).
# Above 1 keV for Cr and below 1 keV for Al the model is an extrapolation, since the
# "representative sensitivity shapes" of Rubery are above 2.2 keV.
const _XRD_SURFACE = (L = 33.64, wC = 0.2742, wO = 0.2504)     # nm, mass fractions
const _XRD_FIT = Dict(:Cr => (K = 15.84e6, γ = -0.2442), :Al => (K = 35.12e6, γ = 0.1018),
                      :Ni => (K = 15.84e6, γ = -0.2442))
const _XRD_BULK_DENSITY = Dict(:Al => 2.699, :Cr => 7.19, :Ni => 8.908, :Au => 19.32)

function _xrd_model(bulk::Symbol, E::Real, K, γ; L = _XRD_SURFACE.L, wC = _XRD_SURFACE.wC,
                    wO = _XRD_SURFACE.wO)
    μ = (1 - wC - wO) * henke_mass_attenuation(bulk, E) + wC * henke_mass_attenuation(:C, E) +
        wO * henke_mass_attenuation(:O, E)
    x = 1e-7 * L * _XRD_BULK_DENSITY[bulk] * μ
    return K * (E / 1000)^γ * x / (1 + x)
end

"""
    xrd_sensitivity(cathode, E) -> V/GW

Returns the sensitivity of a NIF Dante XRD-31 x-ray diode with an `:Al`, `:Cr` or `:Ni`
cathode at photon energy `E` (eV). The sensitivity is in volts into 50 Ω per GW of x-rays
on the cathode, where 1 MV/GW = 2 × 10⁻⁵ A/W. The model is a Henke total-yield model
fitted to the curves digitized from Rubery et al., RSI 2023. Cr is fitted to Fig. 3b over
0.075–1 keV and Al to Fig. 8b over 1–4.5 keV. Ni takes the Cr parameters as an
*assumption*, since no curve is published for it. See `data/rubery2023_dante.csv` and the
research notes, "Dante channel responses and the unfold".
"""
function xrd_sensitivity(cathode::Symbol, E::Real)
    haskey(_XRD_FIT, cathode) || throw(ArgumentError("no XRD model for :$cathode (have :Al, :Cr, :Ni)"))
    p = _XRD_FIT[cathode]
    return _xrd_model(cathode, E, p.K, p.γ)
end

# ---- mirrors ------------------------------------------------------------------------------
#
# SiO2 grazing-incidence mirrors take the Fresnel reflectivity (Henke f1, f2) at the
# nominal angle × 0.584 with 4.25 nm rms roughness. These values are fitted to the measured
# "3.5° SiO2" curve of Rubery 2023 Fig. 3a, with an rms of 0.037 in reflectivity over
# 0.06–1 keV. At the nominal 3.5° the Fresnel cutoff lies at ~0.5 keV, where the measured
# one is at ~0.9 keV. Either the published angle is not the grazing angle of a bulk SiO2
# surface, or the surface is not bulk SiO2. The same scale is applied to the 2.5° mirror
# (*assumed*).
const _SIO2_ANGLE_SCALE = 2.044 / 3.5
const _SIO2_ROUGHNESS = 4.25

# The Cr/Sc multilayer of channel 2 (1.9°) is digitized from Rubery 2023 Fig. 8a over
# 1.2–4.5 keV and held at its end values outside. Below 1.2 keV the 47 µm polyimide
# transmits < 1e-5.
function _mlm_reflectivity(E::Real)
    Et, R = _rubery_dante()["mlm_reflectivity"]
    E <= Et[1] && return R[1] / 100
    E >= Et[end] && return R[end] / 100
    i = searchsortedlast(Et, E)
    return (R[i] + (E - Et[i]) / (Et[i+1] - Et[i]) * (R[i+1] - R[i])) / 100
end

"""
    GrazingMirror(material, angle; roughness = 0, density = nothing)

A grazing-incidence mirror of `material` at `angle` (degrees from the surface). The material
is an element or compound of [`henke_mass_attenuation`](@ref), such as `:C`, `:Be` and
`:B4C`. Its reflectivity is [`mirror_reflectivity`](@ref) with `roughness` (nm rms) and
`density` (g/cm³), which defaults to that of the material. The reflectivity is Fresnel from
Henke's constants and is *uncalibrated*. The SiO₂ mirrors of NIF Dante take `(:SiO2, angle)`
instead, the model calibrated to Rubery et al. 2023 Fig. 3a ([`DanteResponse`](@ref)).
"""
struct GrazingMirror
    material::Symbol
    angle::Float64
    roughness::Float64
    density::Union{Nothing,Float64}
end
GrazingMirror(material::Symbol, angle::Real; roughness::Real = 0.0, density = nothing) =
    GrazingMirror(material, Float64(angle), Float64(roughness), density === nothing ? nothing : Float64(density))

Base.show(io::IO, m::GrazingMirror) = print(io, m.material, " ", m.angle, "°")

_mirror_reflectivity(m::TabulatedCurve, E) = m(E)
_mirror_reflectivity(m::GrazingMirror, E) =
    mirror_reflectivity(m.material, m.angle, E; roughness = m.roughness, density = m.density)
function _mirror_reflectivity(m, E)
    m === nothing && return 1.0
    kind, angle = m
    kind === :mlm && return _mlm_reflectivity(E)
    return mirror_reflectivity(kind, angle * _SIO2_ANGLE_SCALE, E; roughness = _SIO2_ROUGHNESS)
end

# ---- channels -----------------------------------------------------------------------------

"""
    DanteResponse(number, energy; filters = [], grid = 1, mirror = nothing, cathode,
                  solid_angle, tables = [], field_of_view = nothing)
    DanteResponse(number, energy; solid_angle, response::TabulatedCurve, field_of_view = nothing)

One Dante channel built from its hardware. `filters` are `material => thickness` pairs
(thickness in mm) or `material => (thickness, density)`, with transmission from Henke's
tables ([`filter_transmission`](@ref) with `table = :henke`); `tables` are measured
filter transmissions ([`TabulatedCurve`](@ref)s); `grid` is the open fraction of a 25 µm
tantalum mesh (opaque wires). The `mirror` is `(:SiO2, angle°)` for grazing incidence,
calibrated to Rubery 2023 Fig. 3a; `(:mlm, 1.9)` for the Cr/Sc multilayer of channel 2;
a [`GrazingMirror`](@ref) of another material; or a measured reflectivity as a
`TabulatedCurve`. The `cathode` is an XRD, `:Al`, `:Cr` or `:Ni`
([`xrd_sensitivity`](@ref)), or a measured sensitivity as a `TabulatedCurve` in V/GW. The
second form takes a measured whole-channel `response` (V/GW). `energy` (eV) is the nominal
band energy and `solid_angle` (sr) that of the channel. Fields `E` (eV) and `R` (V/GW)
tabulate the response on 50 eV – 30 keV; a spectrum `dI/dE` (GW/sr/eV) gives the voltage
`V = Ω ∫ (dI/dE) R dE` (Harris et al., RSI 92, 033505, 2021, eq. 1).

`field_of_view` (mm) is the collimated field of the channel, the diameter of the disk
seen at the view center in place of the field of the view ([`DanteView`](@ref)). With
`nothing`, the default, the channel takes the field of the view.
"""
struct DanteResponse <: AbstractDanteChannel
    number::Int
    energy::Float64
    filters::Vector{Pair{Symbol,Union{Float64,Tuple{Float64,Float64}}}}
    tables::Vector{TabulatedCurve}
    grid::Float64
    mirror::Union{Nothing,Tuple{Symbol,Float64},GrazingMirror,TabulatedCurve}
    cathode::Union{Symbol,TabulatedCurve}
    response::Union{Nothing,TabulatedCurve}
    solid_angle::Float64
    E::Vector{Float64}
    R::Vector{Float64}
    field_of_view::Union{Nothing,Float64}
end

_filter_spec(v::Real) = Float64(v)
_filter_spec(v::Tuple) = (Float64(v[1]), Float64(v[2]))

function DanteResponse(number::Integer, energy::Real; filters = (), grid::Real = 1.0, mirror = nothing,
                       cathode = nothing, solid_angle::Real, tables = (), response = nothing,
                       field_of_view = nothing)
    E = _DANTE_GRID
    fov = field_of_view === nothing ? nothing : Float64(field_of_view)
    fov === nothing || fov > 0 || throw(ArgumentError("field_of_view must be positive"))
    if response !== nothing
        response isa TabulatedCurve || throw(ArgumentError("response must be a TabulatedCurve"))
        (isempty(filters) && isempty(tables) && grid == 1 && mirror === nothing && cathode === nothing) ||
            throw(ArgumentError("a measured response replaces filters, grid, mirror and cathode"))
        return DanteResponse(Int(number), Float64(energy), [], TabulatedCurve[], 1.0, nothing, response,
                             response, Float64(solid_angle), E, [response(e) for e in E], fov)
    end
    cathode === nothing && throw(ArgumentError("a cathode (an XRD symbol or a TabulatedCurve) is required"))
    f = Pair{Symbol,Union{Float64,Tuple{Float64,Float64}}}[k => _filter_spec(v) for (k, v) in filters]
    m = mirror === nothing || mirror isa Union{TabulatedCurve,GrazingMirror} ? mirror :
        (Symbol(mirror[1]), Float64(mirror[2]))
    m isa Tuple && !(m[1] in (:SiO2, :mlm)) &&
        throw(ArgumentError("a (kind, angle) mirror is :SiO2 or :mlm; use GrazingMirror(:$(m[1]), $(m[2]))"))
    tb = collect(TabulatedCurve, tables)
    R = map(E) do e
        t = filter_transmission(f, e / 1000; table = :henke)
        for c in tb
            t *= c(e)
        end
        grid * t * _mirror_reflectivity(m, e) * _sensitivity(cathode, e)
    end
    return DanteResponse(Int(number), Float64(energy), f, tb, Float64(grid), m, cathode, nothing,
                         Float64(solid_angle), E, R, fov)
end

_sensitivity(c::Symbol, E) = xrd_sensitivity(c, E)
_sensitivity(c::TabulatedCurve, E) = c(E)

# Returns a channel with field of view `f` (mm), or with the field of the view for
# `nothing`. Other channel types have no field of view.
_with_field(c::DanteResponse, f) =
    DanteResponse(c.number, c.energy, c.filters, c.tables, c.grid, c.mirror, c.cathode, c.response,
                  c.solid_angle, c.E, c.R, f === nothing ? nothing : Float64(f))
_with_field(c, f) = c

_show_filter(k, v::Float64) = "$(k) $(round(1000v; digits = 3)) µm"
_show_filter(k, v::Tuple) = "$(k) $(round(1000v[1]; digits = 3)) µm ($(v[2]) g/cm³)"
_show_part(x::TabulatedCurve) = isempty(x.file) ? "measured" : basename(x.file)

function Base.show(io::IO, c::DanteResponse)
    print(io, "DanteResponse(", c.number, ", ", c.energy, " eV, ")
    fov = c.field_of_view === nothing ? "" : ", field of view $(c.field_of_view) mm"
    c.response === nothing || return print(io, "measured response ", _show_part(c.response), fov, ")")
    parts = [[_show_filter(k, v) for (k, v) in c.filters]; ["filter " * _show_part(t) for t in c.tables]]
    print(io, join(parts, " + "), c.grid < 1 ? " + $(round(100c.grid; digits = 1))% grid" : "",
          c.mirror === nothing ? "" : c.mirror isa TabulatedCurve ? ", mirror " * _show_part(c.mirror) :
          c.mirror isa GrazingMirror ? ", " * sprint(show, c.mirror) :
          ", $(c.mirror[1]) $(c.mirror[2])°",
          ", ", c.cathode isa Symbol ? "$(c.cathode) XRD" : "detector " * _show_part(c.cathode), fov, ")")
end

# Returns the response at E by linear interpolation on the grid, 0 outside it.
function _response(c::DanteResponse, E::Real)
    (E < c.E[1] || E > c.E[end]) && return 0.0
    i = min(searchsortedlast(c.E, E), length(c.E) - 1)
    return c.R[i] + (E - c.E[i]) / (c.E[i+1] - c.E[i]) * (c.R[i+1] - c.R[i])
end

"""
    dante1_responses() -> Vector{DanteResponse}

Returns the 18 NIF Dante 1 channels as configured for ICF experiments from October 2021
(Rubery et al., RSI 2023, Fig. 4). The configuration covers the filters, the 8.8% grids,
two SiO₂ mirrors and the 2–4 keV multilayer, the Al, Cr and Ni XRDs, and the solid
angles. Each channel is built as a [`DanteResponse`](@ref) from the built-in definition
`nif_dante1_2021` ([`diagnostic`](@ref)). The transcription and its assumptions are
recorded in `data/diagnostics/nif_dante1_2021.toml`. The channels are built once.
"""
dante1_responses() = _channels(_builtin_definition("nif_dante1_2021").data["channel"])

# Fills W[c, g] with the voltage of channel c per GW/sr in group g, including the solid
# angle. The power is spread over the group as a Planckian at `T` (eV), or flat for
# `nothing`.
function _channel_matrix(ch::Vector{<:DanteResponse}, E::Vector{Float64}; T = nothing)
    G = length(E) - 1
    W = zeros(length(ch), G)
    for g in 1:G
        isfinite(E[g+1]) || continue
        a, b = max(E[g], 1.0), E[g+1]
        n = max(16, ceil(Int, 1000 * log(b / a)))
        x = exp.(range(log(a), log(b), n + 1))
        wt = T === nothing ? ones(n + 1) : _planck_shape.(x, T)
        dx = diff(x)
        den = sum(dx .* (wt[1:end-1] .+ wt[2:end])) / 2
        den > 0 || continue
        for (c, h) in enumerate(ch)
            r = [_response(h, e) for e in x] .* wt
            W[c, g] = h.solid_angle * sum(dx .* (r[1:end-1] .+ r[2:end])) / 2 / den
        end
    end
    return W
end

# Returns the Planck spectrum per eV, normalized to 1, as x³/(eˣ − 1) · 15/π⁴ / T.
_planck_shape(E, T) = (x = E / T; x > 700 ? 0.0 : _PLANCK_NORM * x^3 / expm1(x) / T)

# ---- unfold --------------------------------------------------------------------------------

"""
    dante_unfold(channels, voltages; area, cos = 1, method = :traditional,
                 below = 1500, mband = (1.0, 2500.0, 600.0), λ = 1e-3, knots = 48,
                 tol = 1e-4, maxiter = 50) -> NamedTuple

Unfolds Dante channel `voltages` (V) into a spectrum, as the data are processed. The
channels are of type [`DanteResponse`](@ref), and nonpositive or non-finite entries are
dropped. `area` (mm²) and `cos`, the cosine of the angle between the line of sight and
the emitter normal, turn the power per steradian into a radiation temperature,
`Tr = (π I / (σ A cos))^¼`.

Both methods start as the Dante analysis does (Rubery et al., RSI 2023 §II.C; Harris et al.,
RSI 92, 033505, 2021 §II). Each channel below `below` eV in nominal energy is matched by a
separate Planckian, which gives its `channel_temperatures` entry, and the start temperature
is their mean.

- `method = :traditional` is Harris's algorithm. The start spectrum is that Planckian plus a
  Gaussian M-band `mband = (amplitude GW/sr/eV, center eV, FWHM eV)`. Harris uses amplitude
  1 at 2.5 keV and a "width" of 600 eV, which `dante_unfold` takes as the FWHM. The response
  of each channel is replaced by a unit-peak Gaussian, centered between its half-maximum
  points with their distance as the FWHM, and the Gaussian amplitudes are adjusted until the
  computed voltage of every channel equals the measured one. In Harris's sweeps the
  amplitudes are set from the highest-energy channel down, each to match its channel, for up
  to `maxiter` sweeps and to `tol`. Duplicate channels (e.g. Dante 1's 1/8, 10/14 and 7/15)
  have identical Gaussians, and the sweep then diverges; where the sweeps do not converge,
  the same equations are solved directly by least squares in the relative voltages, with
  minimum-norm amplitudes. The `solver` field of the result is the method that ran. The
  spectrum may go negative locally, as in the original.
- `method = :smooth` gives a positive, smooth unfold. The spectrum is the start
  Planckian, without an M-band, times a piecewise-linear factor on `knots` log-spaced
  knots over 50 eV – 20 keV. The factor is fitted by nonnegative least squares to the
  relative voltages with a second-difference penalty `λ`, so a Planckian is recovered
  exactly. The cubic-spline unfold of Barnak et al. (RSI 91, 073102, 2020) is similar
  but does not use a prior.

Returns a named tuple with fields:

- `energy` (eV grid), `spectrum` (GW/sr/eV) and `intensity` (TW/sr, 50 eV – 20 keV);
- `Tr` (eV) and `start_temperature`;
- `channel_temperatures`, per channel, `NaN` above `below`;
- `fit` (V), the voltages of the unfolded spectrum;
- `misfit`, the largest |fit/V − 1| over the channels used;
- `solver`, which is `:sweep`, `:direct` or `:nnls`;
- `iterations`, the number of sweeps.
"""
function dante_unfold(ch::AbstractVector{<:DanteResponse}, V::AbstractVector{<:Real}; area::Real,
                      cos::Real = 1.0, method::Symbol = :traditional, below::Real = 1500.0,
                      mband = (1.0, 2500.0, 600.0), λ::Real = 1e-3, knots::Integer = 48,
                      tol::Real = 1e-4, maxiter::Integer = 50)
    length(ch) == length(V) || throw(ArgumentError("one voltage per channel"))
    method in (:traditional, :smooth) || throw(ArgumentError("method must be :traditional or :smooth"))
    E = [e for e in _DANTE_GRID if e <= 20_000.0]
    nE = length(E)
    w = zeros(nE)                                  # trapezoid weights
    for k in 1:nE-1
        h = (E[k+1] - E[k]) / 2
        w[k] += h; w[k+1] += h
    end
    use = [c for c in eachindex(ch) if isfinite(V[c]) && V[c] > 0]
    length(use) >= 2 || throw(ArgumentError("need at least two positive voltages"))
    Rw = [ch[c].solid_angle * _response(ch[c], E[k]) * w[k] for c in eachindex(ch), k in 1:nE]
    scale = 1e3 * area * cos / π * SIGMA_SB          # GW/sr per (eV⁴ · normalized shape)
    planck(T) = scale * T^4 .* _planck_shape.(E, T)
    volts(T, c) = dot(view(Rw, c, :), planck(T))
    function brightness(c)
        lo, hi = log(5.0), log(5000.0)
        for _ in 1:100
            m = (lo + hi) / 2
            volts(exp(m), c) < V[c] ? (lo = m) : (hi = m)
        end
        return exp((lo + hi) / 2)
    end
    Tc = fill(NaN, length(ch))
    low = [c for c in use if ch[c].energy < below]
    isempty(low) && throw(ArgumentError("no channel below $below eV"))
    for c in low
        Tc[c] = brightness(c)
    end
    T0 = sum(Tc[low]) / length(low)
    iterations = 0
    if method === :traditional
        amp, Em, fw = mband
        S0 = planck(T0) .+ amp .* exp.(-4log(2) .* ((E .- Em) ./ fw) .^ 2)
        Gs = zeros(nE, length(ch))
        for c in use
            r = [_response(ch[c], e) for e in E]
            k = argmax(r)
            i = k; while i > 1 && r[i] > r[k] / 2; i -= 1; end
            j = k; while j < nE && r[j] > r[k] / 2; j += 1; end
            μ, f = (E[i] + E[j]) / 2, max(E[j] - E[i], 2 * (E[k] - E[max(k - 1, 1)]))
            Gs[:, c] .= exp.(-4log(2) .* ((E .- μ) ./ f) .^ 2)
        end
        D = Rw * Gs                                 # D[c, j]: V of channel c per unit amplitude j
        V0 = Rw * S0
        a = zeros(length(ch))
        centre(c) = sum(Gs[:, c] .* E) / sum(Gs[:, c])
        order = sort(use; by = centre, rev = true)
        miss(a) = maximum(abs((V0[c] + dot(view(D, c, :), a)) / V[c] - 1) for c in use)
        solver = :sweep
        for it in 1:maxiter
            iterations = it
            for c in order
                a[c] += (V[c] - (V0[c] + dot(view(D, c, :), a))) / D[c, c]
            end
            m = miss(a)
            (m < tol || !isfinite(m)) && break
        end
        if !(miss(a) < tol)
            solver = :direct
            a = zeros(length(ch))
            a[use] = pinv(D[use, use] ./ V[use]; rtol = 1e-12) * ((V[use] .- V0[use]) ./ V[use])
        end
        S = S0 .+ Gs * a
    else
        S0 = planck(T0)
        kn = exp.(range(log(50.0), log(20_000.0), knots))
        Φ = zeros(nE, knots)                         # hat functions in log E
        for k in 1:nE
            t = log(E[k])
            j = clamp(searchsortedlast(kn, E[k]), 1, knots - 1)
            s = (t - log(kn[j])) / (log(kn[j+1]) - log(kn[j]))
            Φ[k, j] += 1 - s; Φ[k, j+1] += s
        end
        P = Φ .* S0                                  # basis: the start spectrum × hats
        A = (Rw * P)[use, :] ./ V[use]
        L = zeros(knots - 2, knots)
        for j in 1:knots-2
            L[j, j] = 1; L[j, j+1] = -2; L[j, j+2] = 1
        end
        y = _nnls(vcat(A, sqrt(λ) .* L), vcat(ones(length(use)), zeros(knots - 2)))
        S = P * y
        solver = :nnls
        iterations = 0
    end
    fit = Rw * S
    I = dot(w, S) / 1e3
    return (; energy = E, spectrum = S, intensity = I,
            Tr = radiation_temperature(π * I / (area * cos)), start_temperature = T0,
            channel_temperatures = Tc, fit,
            misfit = maximum(abs(fit[c] / V[c] - 1) for c in use), solver, iterations)
end
