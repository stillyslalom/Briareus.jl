# Soft x-ray optical constants: Henke atomic scattering factors (f1, f2), photoabsorption
# of elements and compounds from 10 eV to 30 keV, and grazing-incidence mirror
# reflectivity. Used by the Dante channel responses (`dante.jl`) and, below 3 keV or for
# materials outside NIST's table, by `mass_attenuation`.
#
# Data: `data/henke/<element>.nff`, verbatim from the CXRO database
# (https://henke.lbl.gov/optical_constants/sf/, fetched 2026-10-24): B. L. Henke,
# E. M. Gullikson and J. C. Davis, "X-ray interactions: photoabsorption, scattering,
# transmission, and reflection at E = 50–30000 eV, Z = 1–92", At. Data Nucl. Data Tables
# 54, 181 (1993), with CXRO's later updates (Mg, Si, Be; Cr from Delmotte et al., J. Appl.
# Phys. 124, 035107, 2018). Columns E (eV), f1, f2; f1 = -9999 where it is not tabulated.

const _HENKE_ELEMENTS = (:H, :Be, :C, :N, :O, :Mg, :Al, :Si, :Cl, :Ti, :V, :Cr, :Fe, :Co,
                         :Ni, :Cu, :Zn, :Ag, :Au)
# standard atomic weights (g/mol)
const _ATOMIC_WEIGHT = Dict(:H => 1.008, :Be => 9.0122, :C => 12.011, :N => 14.007, :O => 15.999,
                            :Mg => 24.305, :Al => 26.982, :Si => 28.085, :Cl => 35.45, :Ti => 47.867,
                            :V => 50.942, :Cr => 51.996, :Fe => 55.845, :Co => 58.933, :Ni => 58.693,
                            :Cu => 63.546, :Zn => 65.38, :Ag => 107.87, :Au => 196.97)
# compounds: atoms per formula unit and the density (g/cm³) used for filters; the
# densities marked *assumed* are typical values, not the Dante parts' own.
const _COMPOUNDS = Dict{Symbol,Tuple{Vector{Pair{Symbol,Int}},Float64}}(
    :CH => ([:C => 1, :H => 1], 1.05),                          # polystyrene, assumed
    :polystyrene => ([:C => 1, :H => 1], 1.05),
    :lexan => ([:C => 16, :H => 14, :O => 3], 1.20),            # polycarbonate
    :parylene => ([:C => 8, :H => 8], 1.11),                    # parylene N, assumed
    :parylene_c => ([:C => 8, :H => 7, :Cl => 1], 1.289),
    :saran => ([:C => 2, :H => 2, :Cl => 2], 1.70),             # PVDC, assumed (1.63–1.78)
    :polyimide => ([:C => 22, :H => 10, :N => 2, :O => 5], 1.42),
    :kapton => ([:C => 22, :H => 10, :N => 2, :O => 5], 1.42),
    :mylar => ([:C => 10, :H => 8, :O => 4], 1.38),
    :polyethylene => ([:C => 2, :H => 4], 0.93),
    :SiO2 => ([:Si => 1, :O => 2], 2.2))                        # fused silica

const _HENKE_CACHE = Dict{Symbol,NTuple{3,Vector{Float64}}}()
const _HENKE_LOCK = ReentrantLock()

# (E, f1, f2) of an element, read once
function _henke(el::Symbol)
    el in _HENKE_ELEMENTS || throw(ArgumentError("no Henke table for :$el"))
    lock(_HENKE_LOCK) do
        get!(_HENKE_CACHE, el) do
            E, f1, f2 = Float64[], Float64[], Float64[]
            for line in eachline(joinpath(@__DIR__, "..", "data", "henke", lowercase(String(el)) * ".nff"))
                w = split(line)
                length(w) == 3 || continue
                v = tryparse.(Float64, w)
                any(isnothing, v) && continue
                push!(E, v[1]); push!(f1, v[2]); push!(f2, v[3])
            end
            (E, f1, f2)
        end
    end
end

# log-log (f2) and log-linear (f1) interpolation at E (eV, 10–30000)
function _henke_f(el::Symbol, E::Real)
    Et, f1, f2 = _henke(el)
    Et[1] <= E <= Et[end] || throw(ArgumentError("the Henke tables cover 10 eV – 30 keV"))
    i = clamp(searchsortedlast(Et, E), 1, length(Et) - 1)
    s = log(E / Et[i]) / log(Et[i+1] / Et[i])
    g2 = exp(log(f2[i]) + s * (log(f2[i+1]) - log(f2[i])))
    g1 = f1[i] == -9999 || f1[i+1] == -9999 ? NaN : f1[i] + s * (f1[i+1] - f1[i])
    return g1, g2
end

_composition(mat::Symbol) = haskey(_COMPOUNDS, mat) ? _COMPOUNDS[mat][1] :
                            mat in _HENKE_ELEMENTS ? [mat => 1] :
                            throw(ArgumentError("no composition for :$mat"))

const _R_E = 2.8179403262e-13        # classical electron radius (cm)
const _HC = 1.239841984e-4           # hc (eV cm)
const _N_A = 6.02214076e23

"""
    henke_mass_attenuation(material, E) -> cm²/g

Photoabsorption mass attenuation coefficient of an element (`:Al`, `:Cr`, …) or compound
(`:CH`, `:lexan`, `:parylene`, `:saran`, `:polyimide`, `:mylar`, `:SiO2`, …) at photon
energy `E` (**eV**, 10–30 000) from Henke's f₂: μ/ρ = 2 ``r_e`` λ ``N_A`` Σ nᵢ f₂ᵢ / Σ nᵢ Aᵢ
(CXRO tables, `data/henke/`; scattering is not included, < 1 % of μ below ~10 keV for
the filter materials used here). Absorption edges are tabulated as point pairs 0.2 eV
apart, so log–log interpolation keeps them sharp.
"""
function henke_mass_attenuation(mat::Symbol, E::Real)
    comp = _composition(mat)
    λ = _HC / E
    s = sum(n * _henke_f(el, E)[2] for (el, n) in comp)
    M = sum(n * _ATOMIC_WEIGHT[el] for (el, n) in comp)
    return 2 * _R_E * λ * _N_A * s / M
end

"""
    mirror_reflectivity(material, angle, E; roughness = 0, density = nothing) -> Float64

Specular reflectivity of a thick, flat `material` mirror (an element or compound, e.g.
`:SiO2`) at grazing `angle` (degrees from the surface) and photon energy `E` (eV), for
unpolarized light: the mean of the s and p Fresnel reflectivities with the refractive
index n = 1 − δ − iβ from Henke's f₁, f₂ (`data/henke/`), and a Névot–Croce factor
exp(−2 ``k_{z0} k_{z1}`` σ²) for an rms `roughness` σ (nm). `density` (g/cm³) defaults
to the material's (`SiO2` 2.2).
"""
function mirror_reflectivity(mat::Symbol, angle::Real, E::Real; roughness::Real = 0.0,
                             density = nothing)
    comp = _composition(mat)
    ρ = density === nothing ? (haskey(_COMPOUNDS, mat) ? _COMPOUNDS[mat][2] : _FILTER_DENSITY[mat]) :
        Float64(density)
    λ = _HC / E                                           # cm
    M = sum(n * _ATOMIC_WEIGHT[el] for (el, n) in comp)
    N = ρ * _N_A / M                                      # formula units per cm³
    f = sum(n * complex(_henke_f(el, E)...) for (el, n) in comp)
    isnan(real(f)) && throw(ArgumentError("f1 is not tabulated at $E eV"))
    χ = _R_E * λ^2 * N * f / π                            # 2(δ + iβ)
    ε = 1 - conj(χ)                                       # n² ≈ 1 − 2δ + 2iβ (e^{−iωt})
    θ = deg2rad(angle)
    kz0 = sin(θ)
    kz1 = sqrt(ε - cos(θ)^2)                              # principal branch: Im ≥ 0, decaying
    k = 2π / λ
    nc = exp(-2 * k^2 * kz0 * kz1 * (1e-7 * roughness)^2)
    rs = (kz0 - kz1) / (kz0 + kz1) * nc
    rp = (ε * kz0 - kz1) / (ε * kz0 + kz1) * nc
    return (abs2(rs) + abs2(rp)) / 2
end

# ---- tabulated curves ------------------------------------------------------------------

"""
    TabulatedCurve(E, values; extrapolate = :hold, file = "")
    TabulatedCurve(file; extrapolate = :hold)

A measured curve against photon energy `E` (eV, increasing), interpolated linearly in
`E`: a filter's transmission, a mirror's reflectivity, a detector's sensitivity (V/GW)
or a whole channel's response (V/GW). Outside the table it holds its end values
(`extrapolate = :hold`, for factors) or is zero (`:zero`, for a whole response). The
file form reads a CSV with a header line and two columns, energy (eV) and value.
"""
struct TabulatedCurve
    E::Vector{Float64}
    v::Vector{Float64}
    extrapolate::Symbol
    file::String
    function TabulatedCurve(E::AbstractVector{<:Real}, v::AbstractVector{<:Real};
                            extrapolate::Symbol = :hold, file::AbstractString = "")
        length(E) == length(v) || throw(ArgumentError("one value per energy"))
        length(E) >= 2 || throw(ArgumentError("a curve needs at least two points"))
        all(>(0), diff(E)) || throw(ArgumentError("energies must increase"))
        extrapolate in (:hold, :zero) || throw(ArgumentError("extrapolate must be :hold or :zero"))
        all(isfinite, v) || throw(ArgumentError("values must be finite"))
        return new(Float64.(E), Float64.(v), extrapolate, String(file))
    end
end

function TabulatedCurve(file::AbstractString; extrapolate::Symbol = :hold)
    _, rows = read_csv(file)
    all(r -> length(r) >= 2, rows) || throw(ArgumentError("$file: two columns (energy eV, value) expected"))
    E = [parse(Float64, r[1]) for r in rows]
    v = [parse(Float64, r[2]) for r in rows]
    return TabulatedCurve(E, v; extrapolate, file)
end

function (c::TabulatedCurve)(E::Real)
    E < c.E[1] && return c.extrapolate === :zero ? 0.0 : c.v[1]
    E > c.E[end] && return c.extrapolate === :zero ? 0.0 : c.v[end]
    i = min(searchsortedlast(c.E, E), length(c.E) - 1)
    return c.v[i] + (E - c.E[i]) / (c.E[i+1] - c.E[i]) * (c.v[i+1] - c.v[i])
end

Base.show(io::IO, c::TabulatedCurve) =
    print(io, "TabulatedCurve(", length(c.E), " points, ", c.E[1], "–", c.E[end], " eV",
          isempty(c.file) ? "" : ", " * basename(c.file), ")")
