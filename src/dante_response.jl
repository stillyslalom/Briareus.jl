# Dante channel responses built from their parts (filters, grids, grazing-incidence
# mirrors, x-ray diode cathodes) and the Dante unfold: channel voltages back to a
# spectrum, a flux and a radiation temperature. Research notes, "Dante channel responses
# and the unfold (2026-10-24)".

# Energy grid of the responses and of the unfold (eV): 50 eV – 30 keV, 0.21 % steps.
const _DANTE_GRID = exp.(range(log(50.0), log(29_990.0), 3001))

# Digitized curves of Rubery et al., RSI 2023 (OSTI 1975922), Figs. 3 and 8
# (`data/rubery2023_dante.csv`): curve name => (energy eV, value).
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
# Sensitivity (V/GW into 50 Ω) S(E) = K (E/1 keV)^γ x/(1 + x), x = L ρ μ_s(E): Henke's
# total-electron-yield form (yield ∝ the energy absorbed within an escape depth L, so
# S ∝ μ while μL ≪ 1; Henke, Knauer & Premaratne, J. Appl. Phys. 52, 1509, 1981) with a
# surface layer of mass fractions w_C carbon and w_O oxygen over the bulk metal (μ_s the
# mixture's μ/ρ, ρ the metal's density) and a power-law correction E^γ. Fitted in log
# (`examples/dante_channels.jl`) to the measured curves digitized from Rubery 2023: L,
# w_C, w_O, γ and K for Cr to Fig. 3b (0.075–1 keV, rms 12 % against the median-read
# curve; it does not reproduce the edge spikes); K and γ for Al to Fig. 8b (1–4.5 keV,
# Al K edge, rms 7 %) with Cr's surface layer. Ni (no published curve) takes Cr's
# parameters on the Ni bulk (*assumed*), and above 1 keV for Cr and below 1 keV for Al
# the model is an extrapolation, as Rubery's own "representative sensitivity shapes"
# are above 2.2 keV.
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

Sensitivity of a NIF Dante XRD-31 x-ray diode with an `:Al`, `:Cr` or `:Ni` cathode at
photon energy `E` (eV): volts into 50 Ω per GW of x-rays on the cathode (1 MV/GW =
2 × 10⁻⁵ A/W). A Henke total-yield model fitted to the curves digitized from Rubery et
al., RSI 2023 (Cr: Fig. 3b, 0.075–1 keV; Al: Fig. 8b, 1–4.5 keV); Ni takes Cr's
parameters (*assumed*: no published curve). See `data/rubery2023_dante.csv` and the
research notes, "Dante channel responses and the unfold".
"""
function xrd_sensitivity(cathode::Symbol, E::Real)
    haskey(_XRD_FIT, cathode) || throw(ArgumentError("no XRD model for :$cathode (have :Al, :Cr, :Ni)"))
    p = _XRD_FIT[cathode]
    return _xrd_model(cathode, E, p.K, p.γ)
end

# ---- mirrors ------------------------------------------------------------------------------
#
# SiO2 grazing-incidence mirrors: Fresnel reflectivity (Henke f1, f2) at the nominal angle
# × 0.584 with 4.25 nm rms roughness, fitted to Rubery 2023 Fig. 3a's measured "3.5° SiO2"
# curve (rms 0.037 in reflectivity over 0.06–1 keV). At the nominal 3.5° the Fresnel
# cutoff sits at ~0.5 keV where the measured one is at ~0.9 keV, so the published angle
# is not the grazing angle of a bulk SiO2 surface (or the surface is not bulk SiO2);
# the same scale is applied to the 2.5° mirror (*assumed*).
const _SIO2_ANGLE_SCALE = 2.044 / 3.5
const _SIO2_ROUGHNESS = 4.25

# the Cr/Sc multilayer of channel 2 (1.9°): Rubery 2023 Fig. 8a digitized (1.2–4.5 keV),
# held at its end values outside (below 1.2 keV the 47 µm polyimide passes < 1e-5)
function _mlm_reflectivity(E::Real)
    Et, R = _rubery_dante()["mlm_reflectivity"]
    E <= Et[1] && return R[1] / 100
    E >= Et[end] && return R[end] / 100
    i = searchsortedlast(Et, E)
    return (R[i] + (E - Et[i]) / (Et[i+1] - Et[i]) * (R[i+1] - R[i])) / 100
end

function _mirror_reflectivity(m, E)
    m === nothing && return 1.0
    kind, angle = m
    kind === :mlm && return _mlm_reflectivity(E)
    return mirror_reflectivity(kind, angle * _SIO2_ANGLE_SCALE, E; roughness = _SIO2_ROUGHNESS)
end

# ---- channels -----------------------------------------------------------------------------

"""
    DanteResponse(number, energy; filters, grid = 1, mirror = nothing, cathode,
                  solid_angle)

One Dante channel built from its hardware: `filters` (`material => thickness` pairs,
thickness in mm, transmission from Henke's tables, [`filter_transmission`](@ref) with
`table = :henke`), a 25 µm tantalum mesh of open fraction `grid` (opaque wires), a `mirror`
(`(:SiO2, angle°)` grazing-incidence, or `(:mlm, 1.9)` for channel 2's Cr/Sc multilayer)
and an XRD `cathode` ([`xrd_sensitivity`](@ref)); `energy` (eV) is the nominal band
energy and `solid_angle` (sr) the channel's. Fields `E` (eV) and `R` (V/GW) tabulate the
response on 50 eV – 30 keV; a spectrum `dI/dE` (GW/sr/eV) gives the voltage
`V = Ω ∫ (dI/dE) R dE` (Harris et al., RSI 92, 033505, 2021, eq. 1).
"""
struct DanteResponse <: AbstractDanteChannel
    number::Int
    energy::Float64
    filters::Vector{Pair{Symbol,Float64}}
    grid::Float64
    mirror::Union{Nothing,Tuple{Symbol,Float64}}
    cathode::Symbol
    solid_angle::Float64
    E::Vector{Float64}
    R::Vector{Float64}
end

function DanteResponse(number::Integer, energy::Real; filters, grid::Real = 1.0, mirror = nothing,
                       cathode::Symbol, solid_angle::Real)
    f = Pair{Symbol,Float64}[k => Float64(v) for (k, v) in filters]
    m = mirror === nothing ? nothing : (Symbol(mirror[1]), Float64(mirror[2]))
    E = _DANTE_GRID
    R = [grid * filter_transmission(f, e / 1000; table = :henke) * _mirror_reflectivity(m, e) *
         xrd_sensitivity(cathode, e) for e in E]
    return DanteResponse(Int(number), Float64(energy), f, Float64(grid), m, cathode,
                         Float64(solid_angle), E, R)
end

function Base.show(io::IO, c::DanteResponse)
    print(io, "DanteResponse(", c.number, ", ", c.energy, " eV, ",
          join(("$(k) $(round(1000v; digits = 3)) µm" for (k, v) in c.filters), " + "),
          c.grid < 1 ? " + $(round(100c.grid; digits = 1)) % grid" : "",
          c.mirror === nothing ? "" : ", $(c.mirror[1]) $(c.mirror[2])°", ", ", c.cathode, " XRD)")
end

# response at E (linear interpolation on the grid, 0 outside it)
function _response(c::DanteResponse, E::Real)
    (E < c.E[1] || E > c.E[end]) && return 0.0
    i = min(searchsortedlast(c.E, E), length(c.E) - 1)
    return c.R[i] + (E - c.E[i]) / (c.E[i+1] - c.E[i]) * (c.R[i+1] - c.R[i])
end

# NIF Dante 1 hardware for ICF experiments, Rubery et al., RSI 2023, Fig. 4 (thicknesses
# µm; "grid" the 8.8 % open mesh, §II.A: 25 µm Ta per Kemp 2020, "tungsten" in
# Rubery). Channel 2's "47 µm CH" is the polyimide of the 2–4 keV multilayer channel
# (§II.D, Fig. 8); "CH" elsewhere is taken as polystyrene and "parylene" as parylene N
# (*assumed*).
const _DANTE1_HARDWARE = [
    (1, 1170, [:Mg => 22], false, nothing, :Al, 6.51e-7),
    (2, 3060, [:polyimide => 47], false, (:mlm, 1.9), :Al, 1.90e-7),
    (3, 250, [:lexan => 4], true, (:SiO2, 3.5), :Al, 5.74e-7),
    (4, 470, [:V => 1.75], false, (:SiO2, 2.5), :Al, 5.56e-7),
    (5, 710, [:Co => 1.4, :CH => 0.2], true, nothing, :Al, 6.51e-7),
    (6, 870, [:Cu => 2], true, nothing, :Cr, 6.51e-7),
    (7, 940, [:Mg => 5, :Zn => 0.65, :CH => 0.1], true, nothing, :Ni, 6.51e-7),
    (8, 1170, [:Mg => 22], true, nothing, :Al, 6.51e-7),
    (9, 1430, [:Al => 20], true, nothing, :Ni, 6.51e-7),
    (10, 1710, [:Si => 20], false, nothing, :Al, 6.51e-7),
    (11, 3430, [:Fe => 1, :Cr => 0.65, :parylene => 5], true, nothing, :Al, 6.51e-7),
    (12, 2580, [:saran => 66], false, nothing, :Al, 6.51e-7),
    (13, 3000, [:Ag => 4.4], false, nothing, :Al, 6.51e-7),
    (14, 1700, [:Si => 20], false, nothing, :Al, 6.51e-7),
    (15, 940, [:Mg => 5, :Zn => 0.65, :CH => 0.1], true, nothing, :Ni, 6.51e-7),
    (16, 14370, [:Al => 250], false, nothing, :Ni, 1.77e-6),
    (17, 690, [:Co => 1], true, nothing, :Cr, 6.51e-7),
    (18, 3660, [:Fe => 1, :Cr => 0.65, :parylene => 5], false, nothing, :Ni, 1.77e-6)]

const _DANTE1_RESPONSES = Ref{Union{Nothing,Vector{DanteResponse}}}(nothing)
const _DANTE1_RESPONSES_LOCK = ReentrantLock()

"""
    dante1_responses() -> Vector{DanteResponse}

The 18 NIF Dante 1 channels as configured for ICF experiments (Rubery et al., RSI 2023,
Fig. 4: filters, 8.8 % grids, two SiO₂ mirrors and the 2–4 keV multilayer, Al/Cr/Ni
XRDs, solid angles), each built as a [`DanteResponse`](@ref). Computed once.
"""
function dante1_responses()
    lock(_DANTE1_RESPONSES_LOCK) do
        _DANTE1_RESPONSES[] === nothing || return _DANTE1_RESPONSES[]
        _DANTE1_RESPONSES[] = [DanteResponse(n, E; filters = [k => 1e-3t for (k, t) in f],
                                             grid = g ? 0.088 : 1.0, mirror = m, cathode = x,
                                             solid_angle = Ω)
                               for (n, E, f, g, m, x, Ω) in _DANTE1_HARDWARE]
    end
end

# W[c, g]: channel c's voltage per GW/sr in group g, spread over the group as a
# Planckian at `T` (eV; flat for `nothing`), solid angle included
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

# Planck spectrum per eV, normalized to 1 (x³/(eˣ − 1) · 15/π⁴ / T)
_planck_shape(E, T) = (x = E / T; x > 700 ? 0.0 : _PLANCK_NORM * x^3 / expm1(x) / T)

# ---- unfold --------------------------------------------------------------------------------

"""
    dante_unfold(channels, voltages; area, cos = 1, method = :traditional,
                 below = 1500, mband = (1.0, 2500.0, 600.0), λ = 1e-3, knots = 48,
                 tol = 1e-4, maxiter = 50) -> NamedTuple

Unfold Dante channel `voltages` (V; channels of type [`DanteResponse`](@ref); nonpositive
or non-finite entries are dropped) into a spectrum, as the data are processed. `area`
(mm²) and `cos` (of the angle between the line of sight and the emitter's normal) turn
the power per steradian into a radiation temperature, `Tr = (π I / (σ A cos))^¼`.

Both methods start as Dante's analysis does (Rubery et al., RSI 2023 §II.C; Harris et
al., RSI 92, 033505, 2021 §II): each channel below `below` eV (nominal energy) is
matched alone by a Planckian (its `channel_temperatures`), and the start temperature is
their mean.

- `method = :traditional` (Harris's algorithm): the start spectrum is that Planckian
  plus a Gaussian M-band `mband = (amplitude GW/sr/eV, center eV, FWHM eV)` (Harris:
  amplitude 1 at 2.5 keV, "width" 600 eV, taken here as the FWHM). Each channel's
  response is replaced by a unit-peak Gaussian (here centered between its half-maximum
  points, FWHM their distance), and the Gaussians' amplitudes are adjusted until every
  channel's computed voltage equals the measured one: by Harris's sweeps (from the
  highest-energy channel down, each amplitude set to match its channel; up to `maxiter`
  sweeps, to `tol`), and where the sweeps do not converge (duplicate channels, e.g.
  Dante 1's 1/8, 10/14 and 7/15, give identical Gaussians, and the sweep then diverges)
  by solving the same equations directly (least squares in the relative voltages,
  minimum-norm amplitudes; `solver` selects which). The spectrum may go negative locally,
  as in the original.
- `method = :smooth`: a positive, smooth unfold: the start Planckian (no M-band) times a
  piecewise-linear factor on `knots` log-spaced knots (50 eV – 20 keV), fitted by
  nonnegative least squares to the relative voltages with a second-difference penalty
  `λ` on the factor, so a Planckian is recovered exactly (cf. the cubic-spline unfold of
  Barnak et al., RSI 91, 073102, 2020, which uses no prior).

Returns `energy` (eV grid), `spectrum` (GW/sr/eV), `intensity` (TW/sr, 50 eV – 20 keV),
`Tr` (eV), `start_temperature`, `channel_temperatures` (per channel, `NaN` above
`below`), `fit` (V, the voltages of the unfolded spectrum), `misfit` (largest
|fit/V − 1| over the channels used), `solver` (`:sweep`, `:direct` or `:nnls`) and
`iterations` (sweeps).
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

# Lawson–Hanson nonnegative least squares: min ‖Ax − b‖, x ≥ 0
function _nnls(A::AbstractMatrix, b::AbstractVector; tol = 1e-12, maxiter = 10 * size(A, 2))
    n = size(A, 2)
    x = zeros(n)
    P = falses(n)
    r = b - A * x
    w = A' * r
    for _ in 1:maxiter
        cand = [j for j in 1:n if !P[j] && w[j] > tol * norm(w, Inf)]
        isempty(cand) && break
        P[cand[argmax(w[cand])]] = true
        for _ in 1:10n
            idx = findall(P)
            z = zeros(n)
            z[idx] = A[:, idx] \ b
            if all(z[idx] .> 0)
                x = z
                break
            end
            neg = [i for i in idx if z[i] <= 0]
            α = minimum(x[i] / (x[i] - z[i]) for i in neg)
            x .+= α .* (z .- x)
            for i in idx
                x[i] <= 1e-15 && (P[i] = false; x[i] = 0.0)
            end
        end
        w = A' * (b - A * x)
    end
    return x
end
