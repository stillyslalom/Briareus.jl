# Backlighter physics: He-α conversion efficiency CE(hν, I) with its uncertainty, filter
# attenuation, and the photometric chain from the laser power on a foil to photons at a
# detector. Research notes, "Backlighter conversion efficiency and the photometric chain
# (2026-10-24)".

const _HEALPHA_KEV = Dict(:Ti => 4.750, :V => 5.181, :Cr => 5.682, :Mn => 6.181, :Fe => 6.700,
                          :Co => 7.242, :Ni => 7.806, :Cu => 8.392, :Zn => 9.000, :Ga => 9.629,
                          :Ge => 10.279, :Br => 12.35, :Kr => 13.115, :Zr => 16.30, :Nb => 17.16,
                          :Mo => 18.03, :Ag => 22.75)

"""
    kshell_line_energy(element) -> keV

Photon energy of the He-α resonance line (He-like 1s²–1s2p ¹P₁, the "w" line) of a
backlighter element given as a `Symbol` (`:Ti`, `:V`, `:Cr`, `:Mn`, `:Fe`, `:Co`, `:Ni`,
`:Cu`, `:Zn`, `:Ga`, `:Ge`, `:Br`, `:Kr`, `:Zr`, `:Nb`, `:Mo`, `:Ag`). Ti–Ge and Kr are
the commonly quoted line energies (to ~0.01 keV); Br, Zr, Nb, Mo and Ag come from a
screened-hydrogenic fit anchored on them,
`hν ≈ 10.204 eV · (Z − 0.42)² · [1 + 0.0011 (Z − 22)]`, good to ~0.5% (Zr: 16.30 keV,
Barrios et al. 2013 quote 16.3). *derived* for Z ≥ 35.
"""
function kshell_line_energy(el::Symbol)
    haskey(_HEALPHA_KEV, el) || throw(ArgumentError("no He-α line energy for :$el; pass the energy in keV"))
    return _HEALPHA_KEV[el]
end
_line_keV(line::Symbol) = kshell_line_energy(line)
_line_keV(line::Real) = (line > 0 || throw(ArgumentError("the line energy must be positive")); Float64(line))

"""
    KShellCEModel

The fitted He-α conversion-efficiency model returned by [`kshell_ce_model`](@ref):
`coef` (the eight parameters `(a_foil, a_prepulsed, b, c_foil, c_prepulsed, d, s, Lₛ)` of
the form described there; the saturation energy is `hνₛ = 10 keV · exp(Lₛ)`), `cov` (their
covariance), `tau` (the scatter of `ln CE` beyond the stated measurement errors),
`scatter` (the predictive scatter of a measured `ln CE`, `√(tau² + s̄²)` with `s̄` the rms
stated error of the fitted points), `n` (data points fitted), and the data's ranges per
class: `hnu` (keV), `irradiance` (W/cm²) and `pulse` (ns), each a `Dict` from `:foil`
and `:prepulsed` to `(min, max)`.
"""
struct KShellCEModel
    coef::SVector{8,Float64}
    cov::SMatrix{8,8,Float64,64}
    tau::Float64
    scatter::Float64
    n::Int
    hnu::Dict{Symbol,Tuple{Float64,Float64}}
    irradiance::Dict{Symbol,Tuple{Float64,Float64}}
    pulse::Dict{Symbol,Tuple{Float64,Float64}}
end

function Base.show(io::IO, ::MIME"text/plain", m::KShellCEModel)
    @printf(io, "KShellCEModel: %d points, scatter %.2f in ln CE (×%.2f; τ = %.2f); peak ∝ hν^%.2f above %.2f keV, flat below; optimum ∝ hν^%.2f",
            m.n, m.scatter, exp(m.scatter), m.tau, m.coef[3], 10exp(m.coef[8]), m.coef[7])
    for c in (:foil, :prepulsed)
        @printf(io, "\n  %-9s hν %.1f–%.1f keV, I %.2g–%.2g W/cm², τ %g–%g ns", c, m.hnu[c]...,
                m.irradiance[c]..., m.pulse[c]...)
    end
end

# ln CE = a_c + b max(L − Lₛ, 0) + c_c Y + d Y², Y = X − s L, with L = ln(hν/10 keV),
# X = ln(I/1e15 W/cm²); for fixed (s, Lₛ) it is linear in the first six parameters on this
# basis
_ce_basis(L, X, prepulsed::Bool, s, Ls) = (Y = X - s * L;
    SVector{6,Float64}(!prepulsed, prepulsed, max(L - Ls, 0.0), Y * !prepulsed, Y * prepulsed, Y^2))
_ce_lncE(p, L, X, prepulsed::Bool) = dot(_ce_basis(L, X, prepulsed, p[7], p[8]), p[SVector{6}(1:6)])
# gradient of ln CE in all eight parameters (the hinge's derivative taken one-sided)
function _ce_grad(p, L, X, prepulsed::Bool)
    v = _ce_basis(L, X, prepulsed, p[7], p[8])
    Y = X - p[7] * L
    return SVector{8,Float64}(v..., -L * ((prepulsed ? p[5] : p[4]) + 2p[6] * Y), L > p[8] ? -p[3] : 0.0)
end

"""
    kshell_ce_data() -> NamedTuple

The He-α conversion-efficiency compilation behind [`kshell_ce_model`](@ref)
(`data/kshell_ce.csv`, sources and provenance in its header): columns `element`, `Z`,
`hnu` (keV), `irradiance` (10¹⁵ W/cm²), `pulse` (ns), `class` (`"foil"` or
`"prepulsed"`), `ce` and `u` (%; `u` is NaN where the paper gives none), `fit` (1 if
fitted), `facility`, `source` and `note`. Workman & Kyrala (2001; OMEGA, 1 ns,
digitized from the figures), Barrios et al. (2013, 2014; NIF, a picket ahead of a
2–5 ns drive, from the tables and text), Coppari et al. (2019; OMEGA and OMEGA EP Fe,
Cu and Ge foils, 0.6–2 ns with or without a pre-pulse, digitized from their Fig. 16),
and one OMEGA Ti foil (Girard et al. 2005, as recalibrated in the 2016 review; held
out, `fit = 0`, as a check below the fitted lines).
"""
function kshell_ce_data()
    h, rows = read_csv(joinpath(@__DIR__, "..", "data", "kshell_ce.csv"))
    num = ("Z", "hnu", "irradiance", "pulse", "ce", "u", "fit")
    col(i) = h[i] in num ? [isempty(r[i]) ? NaN : parse(Float64, r[i]) for r in rows] :
             String[r[i] for r in rows]
    return NamedTuple{Tuple(Symbol.(h))}(Tuple(col(i) for i in eachindex(h)))
end

const _KSHELL_MODEL = Ref{Union{Nothing,KShellCEModel}}(nothing)
const _KSHELL_LOCK = ReentrantLock()

"""
    kshell_ce_model() -> KShellCEModel

Weighted least-squares fit of `ln CE` to the rows of [`kshell_ce_data`](@ref) with
`fit = 1` (68 points). CE is the He-α line with its satellites, into 4π, over the
incident laser energy, for planar foils driven at 351 nm:

```math
\\ln \\mathrm{CE} = a_c + b\\,\\max(L - L_s, 0) + c_c\\,Y + d\\,Y^2,\\qquad Y = X - s\\,L,
```

with ``L = \\ln(h\\nu/10\\ \\mathrm{keV})`` and ``X = \\ln(I/10^{15}\\ \\mathrm{W/cm^2})``.
Each line's CE is a parabola in `ln I` of shared curvature `d`, peaking at the optimum
irradiance ``I^* \\propto (h\\nu)^s`` (`Y* = −c_c/2d`). The peak CE falls as
``(h\\nu)^b`` above the saturation energy ``h\\nu_s = 10\\ \\mathrm{keV}\\, e^{L_s}`` and is
flat below it. The class `c`, `:foil` (one main pulse) or `:prepulsed` (a pre-pulse or
picket ahead of the main drive), carries its own amplitude `a_c` and `c_c` (so its own
peak and optimum); `b`, `d`, `s` and `hνₛ` are shared. Fitted: `b = −3.33 ± 0.25`
(Workman & Kyrala report `hν^-(2.5–4)`, Barrios et al. `hν^-(3–4)`), `s = 2.5 ± 0.4`,
`hνₛ = 6.1 ± 0.4 keV`, between Cr and Fe; foil optimum 1.2e15 W/cm² for Cr, 1.8e15 for
Fe, 5.2e15 for Ge, peak 3.2% at and below hνₛ; prepulsed optimum 2.0e15 for Ge.

Why this form (research notes, "He-α CE: saturating peak (2026-10-26)"): the earlier
form, `a_c + bL + c_c X + dX² + eXL`, has a peak `a_c + bL + (c_c + eL)²/(4|d|)` that is
convex in `L`, so it rose steeply toward low `hν` (Cr 4.1%, Ti 9.2%) against the Cr data
(3.0% at 0.93e15) and the held-out Ti foil. Writing the parabola about its optimum
(`Y`) removes that convexity at the same number of parameters and lowers the
cross-validated error (leave-one-source-out rms 0.55 → 0.44 in `ln CE`); the
saturation, one more parameter, is what the Cr points require (they sit 0.24 below the
power law through Fe–Mo) and what predicts the held-out Ti foil within 1.1σ. Inside the
data the saturation is weakly determined (it moves the likelihood by 1.1 in `−2 ln L`),
and with Cr left out `hνₛ` is not identified. Plausible reasons for a plateau below Fe
(not modeled): absorption and the fraction of it radiated cannot grow without bound, and
at the lower temperatures that suffice for low-Z lines the energy goes increasingly into
lower ionization stages and the continuum.

**Weights.** Each point carries `σᵢ² = sᵢ² + τ²`: `sᵢ = u/ce`, the paper's error bar as a
relative error in `ln CE`; where a row gives none (Workman & Kyrala's digitized points),
`sᵢ = √(0.1² + 0.2²)`, their stated relative error and shot-to-shot reproducibility.
`τ`, the scatter the stated errors do not cover (calibrations, spot definitions, bands
and the model's form), is set so that the weighted χ² per degree of freedom is 1
(random effects): 0.23. For given `(s, Lₛ)` the other six parameters are linear;
`(s, Lₛ)` maximize the likelihood on a refined grid, and the covariance of all eight is
`(JᵀWJ)⁻¹`, `J` the gradient of `ln CE` at the data. The model's `scatter`, the spread of
a measured CE about the fit, is `√(τ² + s̄²)` with `s̄` the rms `sᵢ`: 0.35 in `ln CE` (a
factor of 1.43); outside the data the parameters' covariance adds to it (see
[`kshell_ce`](@ref)).

**Pulse length and spot size** do not enter. The foil data are 1 ns (OMEGA) except
Coppari et al.'s Fe at 0.6 and 2 ns; a `ln τ` term on the foil class fits
`CE ∝ τ^(0.6 ± 0.2)`, but only through Coppari's 2 ns Fe against Workman's 1 ns Fe (within
Coppari's own shots the 1 ns Fe point lies above the 2 ns ones), so the pulse length
stays confounded with the source and, in the prepulsed class, with the picket. Every
source varies the irradiance through the spot size, so the spot size is confounded with
the irradiance.

**Below Cr** the plateau is an extrapolation. The one held-out Ti foil (OMEGA, 1 ns,
1.3e15 W/cm²; `fit = 0`) measured 1.6% (the He-α share of a 4–6 keV band, *inferred*;
2.3% for the whole band) where the fit gives 3.0% ×/÷ 1.48: −0.63 in `ln CE`, 1.1σ with
the point's own 42%. Ti foils on Gekko XII and Shenguang II gave 1.4–1.8% in the band, so
the data may turn over below Cr rather than level off; the fit does not separate the
two. Ti is
flagged `extrapolated`; pass a measured `ce` for Ti backlighters when one exists. The
fitted Fe peak (2.3%) sits below Coppari et al.'s two-sided Fe (2.8–3.25%) and
Workman's 3.0% at 3.5e15, but agrees with Workman's Fe near 1e15 (2.4%): the sources
disagree by about the two-sided factor.

*Engineering model*: confidence medium inside the data, low beyond it. The compilation
mixes facilities, spectrometer calibrations, spot definitions, one- and two-sided
irradiation and bands (Workman's 550 eV band may include part of Ly-α). The fit runs
once per session, from the CSV.
"""
function kshell_ce_model()
    m = _KSHELL_MODEL[]
    m === nothing || return m
    lock(_KSHELL_LOCK) do
        _KSHELL_MODEL[] === nothing && (_KSHELL_MODEL[] = _fit_kshell_ce(kshell_ce_data()))
    end
    return _KSHELL_MODEL[]
end

# relative error in ln CE taken for a row without its own error bar: Workman & Kyrala's
# stated ±10% relative error and 20% reproducibility in quadrature
const _KSHELL_UNSTATED = sqrt(0.1^2 + 0.2^2)

# inverse-variance weights 1/(sᵢ² + τ²), τ set so that χ² = dof (zero when the stated
# errors already cover the scatter); returns the coefficients, χ², τ and the weights
function _kshell_wls(A, y, s, dof)
    function solve(t)
        w = 1 ./ (s .^ 2 .+ t^2)
        β = (A' * (w .* A)) \ (A' * (w .* y))
        return β, sum(w .* abs2.(y - A * β)), w
    end
    lo, hi = 0.0, 4.0
    if solve(0.0)[2] <= dof
        hi = 0.0
    else
        for _ in 1:60
            mid = (lo + hi) / 2
            solve(mid)[2] > dof ? (lo = mid) : (hi = mid)
        end
    end
    β, χ², w = solve(hi)
    return β, χ², hi, w
end

# For each (s, Lₛ) the other six parameters are a weighted linear least-squares fit with
# its own τ; (s, Lₛ) minimize the Gaussian −2 ln L, Σ ln(sᵢ² + τ²) + χ², on a grid refined
# four times (s in 0–5, hνₛ between the data's lowest and highest line). The covariance is
# (JᵀWJ)⁻¹, J the gradient of ln CE in all eight parameters at the data.
function _fit_kshell_ce(d)
    k = findall(==(1.0), d.fit)
    pre = [d.class[i] == "prepulsed" for i in k]
    L = [log(d.hnu[i] / 10) for i in k]
    X = [log(d.irradiance[i]) for i in k]
    y = [log(d.ce[i] / 100) for i in k]
    s = [isnan(d.u[i]) ? _KSHELL_UNSTATED : d.u[i] / d.ce[i] for i in k]
    dof = length(k) - 8
    function profile(sl, Ls)
        A = zeros(length(k), 6)
        for j in eachindex(k)
            A[j, :] = _ce_basis(L[j], X[j], pre[j], sl, Ls)
        end
        all(iszero, view(A, :, 3)) && return Inf, nothing, 0.0, nothing   # hinge above every line
        β, χ², τ, w = _kshell_wls(A, y, s, dof)
        return sum(log, 1 ./ w) + χ², β, τ, w
    end
    Lmin, Lmax = extrema(L)
    best = (Inf, 0.0, 0.0)
    for sl in range(0, 5; length = 21), Ls in range(Lmin, Lmax; length = 41)
        f = profile(sl, Ls)[1]
        f < best[1] && (best = (f, sl, Ls))
    end
    hs, hL = 0.25, (Lmax - Lmin) / 40
    for _ in 1:4
        hs, hL = hs / 4, hL / 4
        c = best
        for sl in c[2] .+ hs .* (-4:4), Ls in c[3] .+ hL .* (-4:4)
            (0 <= sl <= 5 && Lmin <= Ls <= Lmax) || continue
            f = profile(sl, Ls)[1]
            f < best[1] && (best = (f, sl, Ls))
        end
    end
    _, sl, Ls = best
    _, β, τ, w = profile(sl, Ls)
    p = SVector{8,Float64}(β..., sl, Ls)
    J = reduce(vcat, [_ce_grad(p, L[j], X[j], pre[j])' for j in eachindex(k)])
    C = inv(Symmetric(J' * (w .* J)))
    scatter = sqrt(τ^2 + sum(abs2, s) / length(s))
    rng(f, c) = extrema(f[i] for (i, q) in zip(k, pre) if q == (c === :prepulsed))
    cls = (:foil, :prepulsed)
    # pulse lengths: the data's shortest up to 5 ns (the NIF rows' "2" stands for 2–5 ns)
    return KShellCEModel(p, SMatrix{8,8}(Matrix(C)), τ, scatter, length(k),
                         Dict(c => rng(d.hnu, c) for c in cls),
                         Dict(c => 1e15 .* rng(d.irradiance, c) for c in cls),
                         Dict(c => (rng(d.pulse, c)[1], 5.0) for c in cls))
end

"""
    kshell_ce(line; irradiance, prepulse = false, pulse = nothing,
              model = kshell_ce_model()) -> NamedTuple

He-α conversion efficiency (a fraction, into 4π, over the incident laser energy) of a
planar foil whose `line` is an element `Symbol` (see [`kshell_line_energy`](@ref)) or a
line energy in keV, driven at `irradiance` (W/cm², the spot-average incident value),
from the fit [`kshell_ce_model`](@ref).

Returns `(ce, sigma, lo, hi, optimum, peak, extrapolated)`:
- `ce`: the median CE.
- `sigma`: the 1σ uncertainty of `ln CE`, the model's `scatter` and the fit's
  coefficient uncertainty at this point in quadrature. CE is lognormal here, uncertain by a factor
  `e^σ` either way.
- `lo = ce·e^-σ` and `hi = ce·e^σ`.
- `optimum` and `peak`: the irradiance at which this line's CE peaks, and that CE.
- `extrapolated`: whether the point lies outside its class's data in line energy or
  irradiance, or in pulse length when `pulse` (ns) is given. `pulse` changes nothing
  else.
"""
function kshell_ce(line; irradiance::Real, prepulse::Bool = false, pulse = nothing,
                   model::KShellCEModel = kshell_ce_model())
    irradiance > 0 || throw(ArgumentError("the irradiance must be positive"))
    hν = _line_keV(line)
    L, X = log(hν / 10), log(irradiance / 1e15)
    p = model.coef
    g = _ce_grad(p, L, X, prepulse)
    σ = sqrt(model.scatter^2 + dot(g, model.cov * g))
    Xs = p[7] * L - (prepulse ? p[5] : p[4]) / (2p[6])
    cls = prepulse ? :prepulsed : :foil
    inside(x, r) = r[1] * (1 - 1e-9) <= x <= r[2] * (1 + 1e-9)
    ext = !(inside(hν, model.hnu[cls]) && inside(irradiance, model.irradiance[cls]) &&
            (pulse === nothing || inside(pulse, model.pulse[cls])))
    ce = exp(_ce_lncE(p, L, X, prepulse))
    return (ce = ce, sigma = σ, lo = ce * exp(-σ), hi = ce * exp(σ), optimum = 1e15 * exp(Xs),
            peak = exp(_ce_lncE(p, L, Xs, prepulse)), extrapolated = ext)
end

# Mass attenuation coefficients μ/ρ (cm²/g, with coherent scattering) from NIST's X-Ray
# Mass Attenuation Coefficients (Hubbell & Seltzer, Tables 3 and 4), 3–40 keV, each
# absorption edge listed twice (below, above). Interpolated log-log between the
# tabulated energies.
const _MU_RHO = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}(
    :Be => ([3, 4, 5, 6, 8, 10, 15, 20, 30, 40],
            [21.27, 8.685, 4.369, 2.527, 1.124, 0.6466, 0.307, 0.2251, 0.1792, 0.164]),
    :Al => ([3, 4, 5, 6, 8, 10, 15, 20, 30, 40],
            [788, 360.5, 193.4, 115.3, 50.33, 26.23, 7.955, 3.441, 1.128, 0.5685]),
    :Ti => ([3, 4, 4.966, 4.966, 5, 6, 8, 10, 15, 20, 30, 40],
            [332.3, 151.7, 83.8, 687.8, 683.8, 432.3, 202.3, 110.7, 35.87, 15.85, 4.972, 2.214]),
    :Fe => ([3, 4, 5, 6, 7.112, 7.112, 8, 10, 15, 20, 30, 40],
            [557.6, 256.7, 139.8, 84.84, 53.19, 407.6, 305.6, 170.6, 57.08, 25.68, 8.176, 3.629]),
    :Cu => ([3, 4, 5, 6, 8, 8.979, 8.979, 10, 15, 20, 30, 40],
            [748.8, 347.3, 189.9, 115.6, 52.55, 38.29, 278.4, 215.9, 74.05, 33.79, 10.92, 4.862]),
    :Mo => ([3, 4, 5, 6, 8, 10, 15, 20, 20, 30, 40],
            [2011, 970.3, 545, 337.3, 156.5, 85.76, 28.54, 13.08, 80.55, 28.1, 12.94]),
    :Ta => ([3, 4, 5, 6, 8, 9.881, 9.881, 10, 11.14, 11.14, 11.41, 11.68, 11.68, 15, 20, 30, 40],
            [1838, 922.2, 532.8, 338.2, 163.9, 95.99, 244.3, 237.9, 179.1, 244.9, 239.3, 218.2,
             251.8, 134, 63.34, 21.87, 10.25]),
    :mylar => ([3, 4, 5, 6, 8, 10, 15, 20, 30, 40],
               [128.8, 54.66, 27.92, 16.08, 6.75, 3.481, 1.132, 0.5798, 0.3009, 0.2304]),
    :polyethylene => ([3, 4, 5, 6, 8, 10, 15, 20, 30, 40],
                      [77.43, 32.42, 16.43, 9.432, 3.975, 2.088, 0.7452, 0.4315, 0.2706, 0.2275]))
# densities (g/cm³); under the NIST table Kapton (polyimide) takes Mylar's μ/ρ (close in
# composition). The elements added with the Henke tables (`xray_optics.jl`) carry
# handbook densities; compounds carry their own (`_COMPOUNDS`).
const _FILTER_DENSITY = Dict(:Be => 1.848, :Al => 2.699, :Ti => 4.54, :Fe => 7.874, :Cu => 8.96,
                             :Mo => 10.22, :Ta => 16.65, :mylar => 1.38, :kapton => 1.42,
                             :polyethylene => 0.93, :Mg => 1.738, :Si => 2.329, :V => 6.11,
                             :Cr => 7.19, :Co => 8.90, :Ni => 8.908, :Zn => 7.134, :Ag => 10.49,
                             :Au => 19.32)

"""
    mass_attenuation(material, hν; table = :auto) -> cm²/g

Mass attenuation coefficient μ/ρ of a filter `material` at photon energy `hν` (keV).

- `table = :nist`: NIST's X-Ray Mass Attenuation Coefficients (Hubbell & Seltzer;
  coherent scattering included), 3–40 keV, for `:Be`, `:Al`, `:Ti`, `:Fe`, `:Cu`, `:Mo`,
  `:Ta`, `:mylar`, `:kapton` (Mylar's μ/ρ, a few-percent approximation) and
  `:polyethylene`, interpolated log-log with the absorption edges exact (an energy
  exactly on an edge takes the value above it); between the tables' 1–5 keV spaced
  energies good to a few percent away from edges.
- `table = :henke`: photoabsorption from Henke's f₂ ([`henke_mass_attenuation`](@ref),
  0.01–30 keV), for the elements H, Be, C, N, O, Mg, Al, Si, Cl, Ti, V, Cr, Fe, Co, Ni,
  Cu, Zn, Ag, Au and the compounds `:CH`, `:lexan`, `:parylene`, `:parylene_c`,
  `:saran`, `:polyimide` (= `:kapton`), `:mylar`, `:polyethylene`, `:SiO2`.
- `table = :auto` (default): NIST where it covers the material and energy, Henke
  elsewhere (so the backlighter filters keep NIST's values).
"""
function mass_attenuation(mat::Symbol, hν::Real; table::Symbol = :auto)
    table in (:auto, :nist, :henke) || throw(ArgumentError("table must be :auto, :nist or :henke"))
    key = mat === :kapton ? :mylar : mat
    if table === :henke || (table === :auto && !(haskey(_MU_RHO, key) && 3 <= hν <= 40))
        return henke_mass_attenuation(mat, 1000 * hν)
    end
    haskey(_MU_RHO, key) || throw(ArgumentError("no attenuation table for :$mat"))
    3 <= hν <= 40 || throw(ArgumentError("the attenuation tables cover 3–40 keV"))
    E, μ = _MU_RHO[key]
    i = clamp(searchsortedlast(E, hν), 1, length(E) - 1)
    E[i] == E[i+1] && (i += 1)                     # exactly on an edge: above it
    s = log(hν / E[i]) / log(E[i+1] / E[i])
    return exp(log(μ[i]) + s * (log(μ[i+1]) - log(μ[i])))
end

"""
    filter_transmission(filters, hν; table = :auto) -> Real

Transmission `exp(−Σ (μ/ρ) ρ t)` at normal incidence through `filters`, a collection of
`material => thickness` pairs (thickness in mm), at photon energy `hν` (keV). The
materials and μ/ρ are those of [`mass_attenuation`](@ref) (`table` passed on); the
densities are Be 1.848, Al 2.699, Ti 4.54, Fe 7.874, Cu 8.96, Mo 10.22, Ta 16.65, Mylar
1.38, Kapton 1.42, polyethylene 0.93, Mg 1.738, Si 2.329, V 6.11, Cr 7.19, Co 8.90, Ni
8.908, Zn 7.134, Ag 10.49 and Au 19.32 g/cm³; the other compounds of
[`henke_mass_attenuation`](@ref) carry CH 1.05, Lexan 1.20, parylene (N) 1.11, Saran
1.70 and SiO₂ 2.2 g/cm³ (CH, parylene and Saran densities assumed). Pass
`material => (thickness, density)` to override a density.
"""
function filter_transmission(filters, hν::Real; table::Symbol = :auto)
    τ = 0.0
    for (mat, spec) in filters
        t, ρ = spec isa Tuple ? (Float64(spec[1]), Float64(spec[2])) :
               (Float64(spec), get(_FILTER_DENSITY, mat,
                                   haskey(_COMPOUNDS, mat) ? _COMPOUNDS[mat][2] : NaN))
        isnan(ρ) && throw(ArgumentError("no density for :$mat; pass $mat => (thickness, density)"))
        t >= 0 || throw(ArgumentError("filter thicknesses must be non-negative"))
        τ += mass_attenuation(mat, hν; table) * ρ * 0.1t       # mm → cm
    end
    return exp(-τ)
end

# angular distribution f(μ), μ = cos θ from the foil normal, normalized so ∫ f dΩ = 4π
# over the sphere: isotropic f = 1 (half of it toward the front), Lambertian front face
# f = 4μ (all of it toward the front)
_angular(a::Symbol, μ) = a === :isotropic ? one(μ) :
                         a === :lambertian ? 4 * max(μ, zero(μ)) :
                         throw(ArgumentError("angular must be :isotropic, :lambertian or a function"))
_angular(f, μ) = f(μ)

"""
    BacklighterSignal

Result of [`backlighter_signal`](@ref).

Source fields:
- `line` (keV) and `laser_energy` (J incident on all the foil's elements).
- `irradiance` (W/cm², the power-weighted mean over the foil).
- `ce`, `sigma` (1σ of `ln CE`) and `extrapolated`.
- `xray_energy` (J of He-α into 4π).
- `facing` (fraction of the x-ray energy from foil elements that face the detector).
- `per_sr` (J/sr toward a distant detector along the line of sight, from unobstructed
  elements only).

Detector fields:
- `mode` (`:projection` or `:pinhole`).
- `solid_angle` (sr, seen from the foil center: of 1 mm² of detector, or of the
  pinhole).
- `transmission` (through the filters).
- `fluence` (J/mm² of He-α at the detector, after the filters; for a pinhole, the mean
  over the object's image).
- `photons` (per mm² of detector), with `photons_lo` and `photons_hi` (the CE's 1σ
  range).
- `photons_range` (the minimum and maximum per mm² over the object's image, pinhole
  mode only).
- `pixel` (mm, or `nothing`) and `photons_pixel`.

The object does not attenuate here: these are open-field (bright) levels, to be
multiplied by the object's transmission.
"""
struct BacklighterSignal
    line::Float64
    laser_energy::Float64
    irradiance::Float64
    ce::Float64
    sigma::Float64
    extrapolated::Bool
    xray_energy::Float64
    facing::Float64
    per_sr::Float64
    mode::Symbol
    solid_angle::Float64
    transmission::Float64
    fluence::Float64
    photons::Float64
    photons_lo::Float64
    photons_hi::Float64
    photons_range::Tuple{Float64,Float64}
    pixel::Union{Nothing,Float64}
    photons_pixel::Float64
end

function Base.show(io::IO, ::MIME"text/plain", s::BacklighterSignal)
    @printf(io, "BacklighterSignal (%s): %.4g keV line, %.3g J of laser at %.2g W/cm²\n", s.mode,
            s.line, s.laser_energy, s.irradiance)
    @printf(io, "  CE %.3g%% (×%.2f at 1σ)%s: %.3g J into 4π, %.3g J/sr toward the detector\n",
            100s.ce, exp(s.sigma), s.extrapolated ? ", extrapolated" : "", s.xray_energy, s.per_sr)
    @printf(io, "  filter transmission %.3g: %.3g photons/mm² (1σ %.3g–%.3g)", s.transmission,
            s.photons, s.photons_lo, s.photons_hi)
    s.pixel === nothing || @printf(io, ", %.3g per %.3g mm pixel", s.photons_pixel, s.pixel)
end

const _KEV_J = 1.602176634e-16

"""
    backlighter_signal(target, deposition, los; line, pulse, distance, foil = :backlighter,
                       prepulse = false, pinhole = nothing, pinhole_distance = nothing,
                       pixel = nothing, filters = (), angular = :isotropic,
                       ce = nothing, irradiance = :mean, object = :capsule,
                       center = (0, 0, 0), radius = nothing, nrings = 4,
                       model = kshell_ce_model()) -> BacklighterSignal

The photometric chain from the laser power on a backlighter foil to He-α photons per
unit area at a detector `distance` mm from `center` on the line of sight `los` ((θ, φ)
in degrees, or a direction vector, from `center` toward the detector).

1. **Laser energy.** The power that each foil element receives in `deposition` (TW; all of
   the incident power, as Briareus deposits it), over a flat-top pulse of `pulse` ns.
2. **Conversion.** The CE comes from [`kshell_ce`](@ref) for `line` (an element
   `Symbol` or keV). With `irradiance = :mean`, it is evaluated at the foil's
   power-weighted mean irradiance. With `:element`, each element gets the CE at its own
   irradiance, but σ is still the mean's. `ce = value` overrides the model (σ = 0).
   The He-α energy into 4π is `E = CE · E_laser`.
3. **Angular distribution** (an assumption):
   - `:isotropic` (the default): the CE data were converted to 4π this way, each from
     one view angle.
   - `:lambertian`: an optically thick emitting face, `4 cos θ`, toward the front only.
   - A function `f(cos θ)` normalized to `∫ f dΩ = 4π`.

   Each element radiates `E_e f(cos θ)/4π` per steradian. An element contributes
   nothing if its front faces away from the detector, or if its path to the detector
   crosses hardware (every surface part except the foil and `object`).
4. **Geometry.**
   - With `pinhole = nothing` (point projection), the fluence at the detector point is
     `Σ_e E_e f(cos θ_e) / (4π r_e²)`, where `r_e` is the element–detector distance.
   - With a `pinhole` (diameter, mm) at `pinhole_distance` p mm from `center`, and the
     detector at `distance` (q = distance − p behind the pinhole), this is an imager of
     the object. Each image point receives `R · A_pinhole / q²`, where
     `R = E_e f / (4π A_e cos θ)` is the radiance of the foil element seen through the
     pinhole and the object's disk. The disk has radius `radius`, by default the
     `object`'s projected radius, and is sampled as in [`line_of_sight`](@ref).
     `fluence` is the mean over the disk's image points.
5. **Filters.** The [`filter_transmission`](@ref) of `filters` (`material => mm`
   pairs) at the line energy, at normal incidence.
6. **Photons.** `fluence / hν` per mm², and per `pixel` (mm, in the detector plane) if
   given.

The object does not attenuate here (the open-field level). Only the He-α line is
counted: no continuum, no Ly-α, no detector response. The CE's lognormal 1σ
propagates to `photons_lo` and `photons_hi`; filter and geometry uncertainties are not
included.
"""
function backlighter_signal(t::Target, dep::Deposition, los; line, pulse::Real, distance::Real,
                            foil::Symbol = :backlighter, prepulse::Bool = false,
                            pinhole = nothing, pinhole_distance = nothing, pixel = nothing,
                            filters = (), angular = :isotropic, ce = nothing,
                            irradiance::Symbol = :mean, object::Union{Nothing,Symbol} = :capsule,
                            center = (0.0, 0.0, 0.0), radius = nothing, nrings::Integer = 4,
                            model::KShellCEModel = kshell_ce_model())
    pulse > 0 || throw(ArgumentError("the pulse length must be positive"))
    distance > 0 || throw(ArgumentError("the detector distance must be positive"))
    irradiance in (:mean, :element) || throw(ArgumentError("irradiance must be :mean or :element"))
    d = los isa Tuple{Real,Real} ? port_vector(los...) : normalize(Vec3(los...))
    c = Vec3(center...)
    hν = _line_keV(line)
    fe = elements(t, foil)
    isempty(fe) && throw(ArgumentError("the foil :$foil has no elements"))
    P = [Float64(ForwardDiff.value(dep.element[e])) for e in fe]        # TW
    Ptot = sum(P)
    Ptot > 0 || throw(ArgumentError("no laser power on the foil :$foil"))
    Ie = [1e14 * P[i] / t.elem_area[e] for (i, e) in enumerate(fe)]    # W/cm²
    Ibar = sum(P .* Ie) / Ptot
    if ce === nothing
        est = kshell_ce(hν; irradiance = Ibar, prepulse, pulse, model)
        cee = irradiance === :mean ? fill(est.ce, length(fe)) :
              [p > 0 ? kshell_ce(hν; irradiance = I, prepulse, model).ce : 0.0 for (p, I) in zip(P, Ie)]
        σ, ext = est.sigma, est.extrapolated
    else
        0 <= ce <= 1 || throw(ArgumentError("ce must be a fraction in [0, 1]"))
        cee, σ, ext = fill(Float64(ce), length(fe)), 0.0, false
    end
    EL = 1e3 * pulse * Ptot                                   # J (TW·ns = kJ)
    Ee = 1e3 .* pulse .* P .* cee                             # J into 4π
    Ex = sum(Ee)
    # hardware: surface components other than the foil and the object
    fi = Set(part_pieces(t, foil))
    oi = object === nothing ? Set{Int}() : Set(part_pieces(t, object))
    hw = Set(ci for (ci, cp) in enumerate(t.components)
             if cp.kind === :surface && !(ci in fi) && !(ci in oi))
    skip(id) = !(t.tri_comp[id] in hw)
    ε = _ray_eps(t)
    facing = sum((Ee[i] for (i, e) in enumerate(fe) if dot(t.elem_normal[e], d) > 0); init = 0.0) /
             max(Ex, floatmin())
    per_sr = 0.0
    for (i, e) in enumerate(fe)
        μ = dot(t.elem_normal[e], d)
        (μ > 0 && Ee[i] > 0) || continue
        any_hit(t.bvh, t.elem_centroid[e] + ε * t.elem_normal[e], d, Inf, skip) && continue
        per_sr += Ee[i] * _angular(angular, μ) / 4π
    end
    T = filter_transmission(filters, hν)
    fc = sum(t.elem_area[e] * t.elem_centroid[e] for e in fe) / sum(t.elem_area[e] for e in fe)
    if pinhole === nothing
        Pd = c + Float64(distance) * d
        fl = 0.0
        for (i, e) in enumerate(fe)
            Ee[i] > 0 || continue
            r = Pd - t.elem_centroid[e]
            μ = dot(t.elem_normal[e], r) / norm(r)
            μ > 0 || continue
            any_hit(t.bvh, t.elem_centroid[e] + ε * t.elem_normal[e], r, 1.0, skip) && continue
            fl += Ee[i] * _angular(angular, μ) / (4π * dot(r, r))
        end
        mode, Ω, rng = :projection, 1 / sum(abs2, Pd - fc), (NaN, NaN)
    else
        pinhole_distance === nothing && throw(ArgumentError("give pinhole_distance with pinhole"))
        p = Float64(pinhole_distance)
        0 < p < distance || throw(ArgumentError("need 0 < pinhole_distance < distance"))
        q = distance - p
        Ap = π * (Float64(pinhole) / 2)^2
        Pp = c + p * d
        if radius === nothing
            object === nothing && throw(ArgumentError("give `radius` when there is no `object`"))
            a = maximum(x -> norm((x - c) - dot(x - c, d) * d),
                        (x for ci in part_pieces(t, object) for x in t.components[ci].mesh.vertices))
        else
            a = Float64(radius)
        end
        ftris = [k for k in eachindex(t.tri_comp) if t.tri_comp[k] in fi]
        idx = Dict(e => i for (i, e) in enumerate(fe))
        vals = Float64[]
        for qd in _disk_samples(c, d, a, nrings)
            u = qd - Pp                                       # from the pinhole through the disk
            s, k = _first_hit(t, ftris, Pp, u)
            R = 0.0
            if k > 0
                e = dot(u, t.tri_normal[k]) < 0 ? t.tri_front[k] : t.tri_back[k]
                x = Pp + s * u
                if e > 0 && haskey(idx, e) &&
                   !any_hit(t.bvh, x + ε * t.elem_normal[e], Pp - x, 1.0, skip)
                    μ = -dot(normalize(u), t.elem_normal[e])
                    R = Ee[idx[e]] * _angular(angular, μ) / (4π * t.elem_area[e] * μ)   # J/sr/mm²
                end
            end
            push!(vals, R * Ap / q^2)
        end
        fl = sum(vals) / length(vals)
        mode, Ω = :pinhole, Ap / sum(abs2, Pp - fc)
        rng = (T * minimum(vals) / (hν * _KEV_J), T * maximum(vals) / (hν * _KEV_J))
    end
    fl *= T
    nph = fl / (hν * _KEV_J)
    pix = pixel === nothing ? nothing : Float64(pixel)
    return BacklighterSignal(hν, EL, Ibar, Ex / EL, σ, ext, Ex, facing, per_sr, mode, Ω, T, fl,
                             nph, nph * exp(-σ), nph * exp(σ), rng, pix,
                             pix === nothing ? NaN : nph * pix^2)
end
