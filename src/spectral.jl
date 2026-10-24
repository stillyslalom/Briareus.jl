# Multigroup (spectral) radiosity: the gray solution split into photon-energy groups.
#
# View factors do not depend on photon energy, so each group g obeys the gray system
# with its own source and albedo,
#
#     B_g = S_g + α_g H_g,      H_g = F B_g,
#
# where S_g splits the gray source ηP/A by a Planckian at each element's emission
# temperature, less a laser-produced non-Planckian M-band part (`MBandEmission`) placed
# in its own band. With `thermalize = true` an element instead re-emits everything it
# re-emits thermally: B_g = S_g + b_g(T) Σ_g′ α_g′ H_g′, with b_g(T) the Planck fraction
# of group g. With equal group albedos either form sums to the gray system exactly
# (Σ_g S_g = S, Σ_g b_g = 1), so the group sum of every flux is the gray one.
#
# A diagnostic, not on the optimization path: everything here is Float64 (the values of
# a dual-number gray solution).

const _PLANCK_NORM = 15 / π^4
# Bernoulli-series coefficients of ∫₀ˣ t³/(eᵗ - 1) dt = Σ cₙ x^(n+3): B_n / ((n + 3) n!)
const _PLANCK_SERIES = let B = Dict(0 => 1 // 1, 1 => -1 // 2, 2 => 1 // 6, 4 => -1 // 30,
                                     6 => 1 // 42, 8 => -1 // 30, 10 => 5 // 66,
                                     12 => -691 // 2730, 14 => 7 // 6, 16 => -3617 // 510)
    [(n, Float64(B[n] / ((n + 3) * factorial(big(n))))) for n in sort(collect(keys(B)))]
end

# Fraction of σT⁴ above x = hν/T: the Planck tail, 1 at x = 0 and 0 at x = ∞, computed
# directly (not as 1 - cdf) so that the hard groups keep their relative accuracy.
function _planck_tail(x::Float64)
    x <= 0 && return 1.0
    isinf(x) && return 0.0
    if x < 0.5                       # 1 - (15/π⁴) ∫₀ˣ, the series converging as (x/2π)ⁿ
        s = 0.0
        for (n, c) in _PLANCK_SERIES
            s += c * x^(n + 3)
        end
        return 1 - _PLANCK_NORM * s
    end
    tail = 0.0                       # Σₖ e^(-kx) (x³/k + 3x²/k² + 6x/k³ + 6/k⁴)
    for k in 1:200
        term = exp(-k * x) * (x^3 / k + 3x^2 / k^2 + 6x / k^3 + 6 / k^4)
        tail += term
        term < 1e-17 * tail && break
    end
    return _PLANCK_NORM * tail
end

"""
    planck_fraction(T, E1, E2) -> Float64

Fraction of a Planckian's power (σT⁴) emitted between photon energies `E1` and `E2`
(eV; `E2 = Inf` for the whole tail) at temperature `T` (eV). `planck_fraction(T, 0, Inf)`
is 1 exactly.
"""
function planck_fraction(T::Real, E1::Real, E2::Real)
    E1 <= E2 || throw(ArgumentError("need E1 ≤ E2"))
    T > 0 || return E1 == 0 ? 1.0 : 0.0     # a cold element: everything in the softest group
    return _planck_tail(Float64(E1 / T)) - _planck_tail(Float64(E2 / T))
end

"""
    planck_fractions(T, edges) -> Vector

Planck fractions of the groups `[edges[g], edges[g+1]]` at temperature `T` (eV). With
`edges[1] = 0` and `edges[end] = Inf` they sum to one (to rounding).
"""
planck_fractions(T::Real, edges::AbstractVector) =
    [planck_fraction(T, edges[g], edges[g+1]) for g in 1:length(edges)-1]

"""
    spectral_edges(; lo = 50, hi = 20_000, n = 40, include = (1800,)) -> Vector

Group edges (eV) for `spectral_solve`: `0`, `n` log-spaced edges from `lo` to `hi`, the
energies in `include` (default the M-band threshold, 1.8 keV) and `Inf`. The default
spans the NIF Dante channels (50 eV – 20 keV; Rubery et al., RSI 2023).
"""
function spectral_edges(; lo::Real = 50.0, hi::Real = 20_000.0, n::Integer = 40,
                        include = (1800.0,))
    0 < lo < hi || throw(ArgumentError("need 0 < lo < hi"))
    e = exp.(range(log(lo), log(hi), n))
    return sort(unique(vcat(0.0, e, collect(Float64, include), Inf)))
end

"""
    MBandEmission(; intensity = [1e14, 1e15], fraction = [0.01, 0.20], groups = nothing,
                  band = (1800, 2850), materials = (:Au,))

Non-Planckian gold M-band emission from laser-heated spots, an **engineering model**:
elements of a material in `materials` emit part of their absorbed laser power `P` as
M-band x-rays, flat in photon energy over `band` (eV). The fraction is

- by default `fraction(I)`, with `I` the element's absorbed laser flux `P/A` (W/cm²),
  interpolated in log–log between the `(intensity, fraction)` points and held constant
  beyond them; or
- with `groups` (a `Dict` from beam group to fraction, e.g. `nif_mband()`), fixed per
  beam group: an element emits `Σ_g fraction_g P_g`, `P_g` the power it absorbs from
  group g's beams (`mband_power` deposits the groups separately).

The M-band is part of the converted x-rays, not added to them: the element's Planckian
source is `ηP/A` less the M-band (which is capped at `ηP/A`).

Defaults (provenance in the research notes, "Spectral radiosity, M-band and a synthetic
Dante"): laser → M-band (> 1.8 keV) conversion of gold at normal incidence ~1 % at 10¹⁴
and ~20 % at 10¹⁵ W/cm² (Dewald et al., PoP 15, 072706, 2008, as quoted in Dewald et
al., PoP 27, 2020, **M**); the band's upper edge 2.85 keV puts 1/3 of a flat band above
2.5 keV, Dewald 2020's share of the hohlraum M-band above the Mo edge (**L**). The
element-averaged intensity of the default form depends on the mesh across a spot (the
fraction goes as I^1.3); the per-group form does not.
"""
struct MBandEmission
    intensity::Vector{Float64}
    fraction::Vector{Float64}
    groups::Union{Nothing,Dict{Symbol,Float64}}
    band::NTuple{2,Float64}
    materials::Vector{Symbol}
    function MBandEmission(intensity, fraction, groups, band, materials)
        length(intensity) == length(fraction) >= 1 ||
            throw(ArgumentError("intensity and fraction differ in length"))
        issorted(intensity) && all(>(0), intensity) ||
            throw(ArgumentError("intensities must be positive and sorted"))
        all(f -> 0 <= f <= 1, fraction) || throw(ArgumentError("fractions must lie in [0, 1]"))
        groups === nothing || all(f -> 0 <= f <= 1, values(groups)) ||
            throw(ArgumentError("group fractions must lie in [0, 1]"))
        0 <= band[1] < band[2] < Inf || throw(ArgumentError("need 0 ≤ band[1] < band[2] < ∞"))
        new(collect(Float64, intensity), collect(Float64, fraction),
            groups === nothing ? nothing : Dict{Symbol,Float64}(groups),
            (Float64(band[1]), Float64(band[2])), collect(Symbol, materials))
    end
end
MBandEmission(; intensity = [1e14, 1e15], fraction = [0.01, 0.20], groups = nothing,
              band = (1800.0, 2850.0), materials = (:Au,)) =
    MBandEmission(intensity, fraction, groups, band, materials)

"""
    nif_mband(; outer = 1.5e15, inner = 3e14, kwargs...) -> MBandEmission

Per-cone NIF M-band emission: each cone's laser → M-band conversion is the default
`fraction(I)` at its mean wall irradiance, outer 1.5 × 10¹⁵ and inner 3 × 10¹⁴ W/cm²
(Dewald et al., PoP 27, 2020, **M**): 0.20 (held beyond 10¹⁵) and 0.042. Dewald also
notes the inner beams' 2.5× lower turning-point density, which lowers their M-band
further (∝ n_e², ~5×); that factor is not applied (**L**: the inner M-band is an upper
bound). `kwargs` pass to `MBandEmission`.
"""
function nif_mband(; outer::Real = 1.5e15, inner::Real = 3e14, kwargs...)
    m = MBandEmission(; kwargs...)
    fo, fi = mband_conversion(m, outer), mband_conversion(m, inner)
    return MBandEmission(; kwargs...,
                         groups = Dict(:outer44 => fo, :outer50 => fo, :inner23 => fi, :inner30 => fi))
end

"""
    mband_conversion(m::MBandEmission, I) -> Float64

Laser → M-band conversion at absorbed laser flux `I` (W/cm²).
"""
function mband_conversion(m::MBandEmission, I::Real)
    I > 0 || return 0.0
    x, y = m.intensity, m.fraction
    I <= x[1] && return y[1]
    I >= x[end] && return y[end]
    k = searchsortedlast(x, I)
    w = log(I / x[k]) / log(x[k+1] / x[k])
    (y[k] > 0 && y[k+1] > 0) || return y[k] + w * (y[k+1] - y[k])
    return exp(log(y[k]) + w * (log(y[k+1]) - log(y[k])))
end

"""
    mband_power(m::MBandEmission, target, laser) -> Vector
    mband_power(m::MBandEmission, op::DepositionOperator, beams) -> Vector

Laser-produced M-band power per element (TW): from the absorbed laser power per element
`laser` (TW; the default intensity form), or by depositing `beams` with `op` (either
form; the per-group form deposits each group's beams separately). Elements of other
materials and apertures emit none. `spectral_solve` takes the result as `mband_power`.
"""
function mband_power(m::MBandEmission, t::Target, P::AbstractVector)
    m.groups === nothing ||
        throw(ArgumentError("per-group M-band: pass the deposition operator and the beams"))
    out = zeros(nelements(t))
    for e in 1:nelements(t)
        _emits_mband(m, t, e) || continue
        p = ForwardDiff.value(P[e])
        p > 0 && (out[e] = mband_conversion(m, to_Wcm2(p / t.elem_area[e])) * p)
    end
    return out
end

function mband_power(m::MBandEmission, op::DepositionOperator, beams::AbstractVector{<:Beam})
    t = op.target
    m.groups === nothing && return mband_power(m, t, deposit(op, beams).element)
    out = zeros(nelements(t))
    for (g, f) in m.groups
        any(b -> b.group === g, beams) || continue
        P = deposit(op, [b.group === g ? b : remake(b; power = 0.0) for b in beams]).element
        for e in 1:nelements(t)
            _emits_mband(m, t, e) && (out[e] += f * ForwardDiff.value(P[e]))
        end
    end
    return out
end

_emits_mband(m::MBandEmission, t::Target, e) =
    !is_aperture(t, e) && element_material(t, e).name in m.materials

# share of a flat band [a, b] in each group
_band_shares(edges, (a, b)) =
    [max(0.0, min(b, edges[g+1]) - max(a, edges[g])) / (b - a) for g in 1:length(edges)-1]

"""
    SpectralSolution

Result of `spectral_solve`. Fields: `gray` (the `RadiationSolution` it splits),
`edges` (group edges, eV); per element × group (TW/mm²): `source`, `B` (radiosity), `H`
(incident flux), `albedo`; per element: `temperature` (eV, the wall temperature: the
Planckian of the thermal re-emission and of any stored-energy emission; with `spot =
:brightness` the one temperature of the whole split), `spot_temperature` (eV, the
Planckian of the converted laser x-rays less the M-band), `converted` (the converted
laser source ηP/A, TW/mm²), `mband` (the laser-produced M-band source, TW/mm²);
`mband_shares` (its split over the groups); `thermalize`; `spot` (the spot-spectrum
rule, see `spectral_solve`).
"""
struct SpectralSolution
    gray::RadiationSolution
    edges::Vector{Float64}
    source::Matrix{Float64}
    B::Matrix{Float64}
    H::Matrix{Float64}
    albedo::Matrix{Float64}
    temperature::Vector{Float64}
    spot_temperature::Vector{Float64}
    converted::Vector{Float64}
    mband::Vector{Float64}
    mband_shares::Vector{Float64}
    thermalize::Bool
    spot::Symbol
end
# a single-temperature split (the `:brightness` rule): every source at `temperature`
SpectralSolution(gray::RadiationSolution, edges, source, B, H, albedo, temperature::AbstractVector,
                 mband::AbstractVector, mband_shares, thermalize::Bool) =
    SpectralSolution(gray, edges, source, B, H, albedo, temperature, copy(temperature),
                     zeros(length(temperature)), mband, mband_shares, thermalize, :brightness)

"Number of photon-energy groups of a `SpectralSolution`."
ngroups(s::SpectralSolution) = length(s.edges) - 1

function Base.show(io::IO, s::SpectralSolution)
    print(io, "SpectralSolution(", size(s.B, 1), " elements, ", ngroups(s), " groups",
          s.thermalize ? ", thermalized" : "", ", spots at ", s.spot, ")")
end

"""
    spectral_solve(sol::RadiationSolution; edges = [0, Inf], albedo = nothing,
                   mband = nothing, mband_power = nothing, thermalize = false,
                   spot = :radiation, temperature = nothing, spot_temperature = nothing,
                   tol = 1e-6, maxiter = 20) -> SpectralSolution

Multigroup radiosity about a gray solution `sol` (same model, laser deposition and
time). The groups are `[edges[g], edges[g+1]]` (eV; `edges[1] = 0`, `edges[end] = Inf`,
e.g. `spectral_edges()`).

- Sources: an element's gray source is its converted laser x-rays `ηP/A` plus, on a
  re-emitting wall of a `drive_history`, the emission of its stored energy. The
  laser-produced M-band (`mband`, an `MBandEmission`; its power per element
  `mband_power`, TW, defaults to `mband_power(mband, target, sol.laser)`, and must be
  given for a per-group model) is part of the converted x-rays and goes into its band.
  The rest is split by Planck fractions, by `spot`:
  - `:radiation` (default): the converted x-rays less the M-band at the element's
    radiation temperature `(H/σ)^¼` (its incident flux; the coronal emission of a laser
    spot is thermal at about the cavity's temperature, its hard part non-Planckian M-band:
    research notes, "Laser-spot spectrum"); the stored-energy emission (and, with
    `thermalize`, the thermal re-emission) at the wall temperature `((B − ηP/A)/σ)^¼`,
    the brightness of the element's radiosity without its converted laser part;
  - `:wall`: the converted x-rays also at the wall temperature;
  - `:brightness`: everything at the brightness temperature `(B/σ)^¼` of the total
    radiosity (the rule before 2026-10-24's spot spectrum; a laser spot then emits as a
    Planckian well above the wall's temperature, with a hard tail of its own).
  The temperatures follow the group totals when those differ from the gray ones
  (iterated to `tol`). A given `temperature` vector (eV) fixes the wall (with
  `:brightness` the only) temperature, `spot_temperature` the converted x-rays' one.
  On elements without laser light the three rules agree.
- Group albedos `α_g`: the fraction of an element's incident flux in group g that it
  re-emits in group g. The gray converged albedo `sol.albedo` in every group, except
  for the components in `albedo`, a `Dict` (or pairs) from component name to a vector
  of group albedos (`NaN` keeps the gray albedo in that group).
- `thermalize`: what happens to the rest of the gray albedo, `α − α_g`. Lost to the
  wall by default (`B_g = S_g + α_g H_g`, each group a gray-type solve); with
  `thermalize = true` re-emitted thermally, spread over the groups by the element's
  Planck fractions `b_g(T)` at its wall temperature:
  `B_g = S_g + α_g H_g + b_g(T) Σ_g′ (α − α_g′) H_g′` (group
  albedos above `α` are lowered to it). Thermalized, the group sums are the gray
  solution's whatever the group albedos; with every group albedo 0, all re-emission is
  thermal.

With every group albedo equal to the gray one (either way), the group sums of `source`,
`B` and `H` equal the gray solution's (absorbed and emitted power per element too),
whatever the split; with one group and no M-band the solution is the gray one exactly.
Groups at the gray albedo reuse the model's cached factorizations: solved directly when
that albedo is factored, by GMRES preconditioned with a cached factorization otherwise
(fresh factorizations are not cached, so the gray `solve` path is unaffected). Groups
with their own albedos are solved by GMRES (jointly when thermalized), which converges
fast for the low albedos of hard groups.
"""
function spectral_solve(sol::RadiationSolution; edges = [0.0, Inf], albedo = nothing,
                        mband::Union{Nothing,MBandEmission} = nothing, mband_power = nothing,
                        thermalize::Bool = false, spot::Symbol = :radiation, temperature = nothing,
                        spot_temperature = nothing, tol = 1e-6, maxiter::Integer = 20)
    E = collect(Float64, edges)
    length(E) >= 2 && E[1] == 0 && E[end] == Inf && issorted(E) && allunique(E) ||
        throw(ArgumentError("edges must increase from 0 to Inf"))
    spot in (:radiation, :wall, :brightness) ||
        throw(ArgumentError("spot must be :radiation, :wall or :brightness, got :$spot"))
    m = sol.model
    t = m.target
    n = nelements(t)
    G = length(E) - 1
    v = ForwardDiff.value
    S = Float64[v(x) for x in sol.source]
    Bgray = Float64[v(x) for x in sol.B]
    Hgray = Float64[v(x) for x in sol.H]
    α = sol.albedo
    A = _group_albedos(t, α, albedo, G)
    thermalize && (A .= min.(A, α))
    # laser-produced M-band
    SM = zeros(n)
    shares = zeros(G)
    if mband !== nothing
        shares = _band_shares(E, mband.band)
        PM = mband_power === nothing ? Briareus.mband_power(mband, t, sol.laser) : mband_power
        length(PM) == n || throw(ArgumentError("mband_power must have one entry per element"))
        SM .= min.(PM ./ t.elem_area, S)
    end
    # converted laser x-rays (the rest of the source is a re-emitting wall's stored energy)
    SC = spot === :brightness ? copy(S) :
         clamp.(Float64[v(x) for x in _xray_source(t, v.(sol.laser), sol.time)], SM, S)
    # wall temperature (the only one with :brightness) and the converted x-rays'
    wallT(Btot) = spot === :brightness ? radiation_temperature.(Btot) :
                  radiation_temperature.(max.(Btot .- SC, 0.0))
    spotT(Btot, Htot) = spot === :radiation ? radiation_temperature.(Htot) : wallT(Btot)
    T = temperature === nothing ? wallT(Bgray) : collect(Float64, temperature)
    length(T) == n || throw(ArgumentError("temperature must have one entry per element"))
    Tc = spot_temperature !== nothing ? collect(Float64, spot_temperature) :
         spot === :brightness ? T : spot === :wall ? T : spotT(Bgray, Hgray)
    length(Tc) == n || throw(ArgumentError("spot_temperature must have one entry per element"))
    soft = [g for g in 1:G if view(A, :, g) == α]     # groups at the gray albedo
    hard = setdiff(1:G, soft)
    gray = G == 1 && isempty(hard) && all(iszero, SM)
    local B, Sg
    for it in 1:maxiter
        b = _planck_matrix(T, E)
        bc = Tc == T ? b : _planck_matrix(Tc, E)
        Sg = (SC .- SM) .* bc .+ (S .- SC) .* b .+ SM .* shares'
        if gray
            Sg = reshape(copy(S), n, 1)         # the gray system itself
            B = reshape(copy(Bgray), n, 1)
        elseif thermalize && all(iszero, A)
            B = Sg .+ b .* (α .* Hgray)         # all re-emission thermal: R = αH exactly
        else
            B = similar(Sg)
            R = zeros(n)
            if !isempty(hard)
                X = thermalize ? _thermalized_groups(m.F, A[:, hard], b[:, hard], Sg[:, hard], α) :
                    _group_radiosity(m, A[:, hard], Sg[:, hard], α)
                B[:, hard] .= X
                thermalize && (R = vec(sum((α .- A[:, hard]) .* (m.F * X); dims = 2)))
            end
            isempty(soft) ||
                (B[:, soft] .= _group_radiosity(m, A[:, soft], Sg[:, soft] .+ b[:, soft] .* R, α))
        end
        # thermalized or at the gray albedo, the totals (and so T) are the gray ones
        fixedT = temperature !== nothing
        fixedTc = spot_temperature !== nothing || (spot !== :radiation && fixedT)
        (fixedT && fixedTc || isempty(hard) || thermalize) && break
        Btot = vec(sum(B; dims = 2))
        Tn = fixedT ? T : wallT(Btot)
        Tcn = fixedTc ? Tc : spot === :radiation ? radiation_temperature.(m.F * Btot) : Tn
        rel(x, y) = maximum(abs.(x .- y) ./ max.(y, 1.0))
        max(rel(Tn, T), rel(Tcn, Tc)) < tol && break   # keep the T B was solved with
        it == maxiter && @warn "spectral temperature iteration did not converge"
        T, Tc = Tn, Tcn
    end
    H = gray ? reshape(copy(Hgray), n, 1) : m.F * B
    return SpectralSolution(sol, E, Sg, B, H, A, T, Tc, SC, SM, shares, thermalize, spot)
end

# Planck fractions per element (rows) and group (columns).
function _planck_matrix(T::AbstractVector, E::AbstractVector)
    G = length(E) - 1
    b = Matrix{Float64}(undef, length(T), G)
    for e in eachindex(T)
        if G == 1
            b[e, 1] = 1.0
            continue
        end
        # tails at the edges: 1 at E = 0, 0 at Inf (a cold element: all in the first group)
        x = T[e] > 0 ? [_planck_tail(E[g] / T[e]) for g in 1:G+1] : [i == 1 ? 1.0 : 0.0 for i in 1:G+1]
        for g in 1:G
            b[e, g] = x[g] - x[g+1]
        end
    end
    return b
end

# Group albedos (elements × groups): the gray albedo, overridden per component.
function _group_albedos(t::Target, α::Vector{Float64}, spec, G::Integer)
    A = repeat(α, 1, G)
    spec === nothing && return A
    for (name, a) in spec
        length(a) == G || throw(ArgumentError("$name: $(length(a)) group albedos for $G groups"))
        for e in elements(t, name), g in 1:G
            isnan(a[g]) || (A[e, g] = a[g])
        end
    end
    return A
end

# A cached factorization of the model to solve or precondition with: the one for `α`
# when cached (`exact = true`), else the iterative preconditioner or the most recent one.
function _cached_lu(m::RadiationModel, α::Vector{Float64})
    c = m.cache
    lock(c.lock) do
        i = findfirst(==(α), c.keys)
        i !== nothing && return (c.facts[i], true)
        c.precond[] !== nothing && return (c.precond[], false)
        isempty(c.facts) || return (c.facts[1], false)
        return (nothing, false)
    end
end

# In-group radiosity, B_g = (I - diag(α_g) F)⁻¹ S_g, one solve per distinct albedo vector.
function _group_radiosity(m::RadiationModel, A::Matrix{Float64}, S::Matrix{Float64},
                          αgray::Vector{Float64})
    n, G = size(S)
    B = similar(S)
    done = falses(G)
    for g in 1:G
        done[g] && continue
        αg = A[:, g]
        cols = [k for k in g:G if !done[k] && view(A, :, k) == αg]
        done[cols] .= true
        if m.iterative
            B[:, cols] .= _iterative_radiosity(m, αg, S[:, cols], nothing)
            continue
        end
        f, exact = _cached_lu(m, αg)
        if f !== nothing && exact
            B[:, cols] .= f \ S[:, cols]
            continue
        end
        local_lu = nothing
        f === nothing && (f, _ = _cached_lu(m, αgray))
        for k in cols
            x = f === nothing ? nothing :
                _gmres_radiosity(m.F, αg, f, S[:, k], f \ S[:, k])
            if x === nothing
                local_lu === nothing && (local_lu = lu!(_radiosity_matrix(m.F, αg)))
                x = local_lu \ S[:, k]
            end
            B[:, k] .= x
        end
    end
    return B
end

# Groups h with their own albedos, thermalized: X_h - a_h∘(F X_h) - b_h∘R = S_h with
# R = Σ_h (α - a_h)∘(F X_h), solved jointly by (unpreconditioned) GMRES; the iteration
# matrix has spectral radius about max(a_h) + Σ_h b_h (α - a_h), small for hard groups.
function _thermalized_groups(F::AbstractMatrix{Float64}, Ah::Matrix{Float64}, bh::Matrix{Float64},
                             Sh::Matrix{Float64}, α::Vector{Float64})
    n, Gh = size(Sh)
    loss = α .- Ah
    function op!(y, x)
        X = reshape(x, n, Gh)
        FX = F * X
        R = vec(sum(loss .* FX; dims = 2))
        y .= vec(X .- Ah .* FX .- bh .* R)
        return y
    end
    rhs = vec(Sh)
    x = _gmres(op!, nothing, rhs, copy(rhs); maxiter = 100)
    x === nothing && error("thermalized group solve did not converge (group albedos close to the gray albedo: use thermalize = false)")
    return reshape(x, n, Gh)
end

"""
    band_flux(s::SpectralSolution, lo, hi = Inf; field = :H) -> Vector

Per-element sum of the groups between photon energies `lo` and `hi` (eV, both group
edges) of `field` (`:H` incident flux, `:B` radiosity, `:source`).
"""
function band_flux(s::SpectralSolution, lo::Real, hi::Real = Inf; field::Symbol = :H)
    g = _band_groups(s.edges, lo, hi)
    X = getfield(s, field)
    return vec(sum(view(X, :, g); dims = 2))
end

function _band_groups(E, lo, hi)
    i, j = findfirst(==(lo), E), findfirst(==(hi), E)
    (i === nothing || j === nothing || j <= i) &&
        throw(ArgumentError("band [$lo, $hi] eV must start and end on group edges"))
    return i:j-1
end

"""
    mband_recirculation(s::SpectralSolution, elems; lo = 1800) -> Float64

Dewald's F_M on the elements `elems`: the ratio of recirculating to laser-produced
M-band flux (photon energies ≥ `lo`, eV), where the laser-produced part is the first
flight of the laser-heated spots' non-Planckian M-band (`F · S_M`) and the recirculating
part everything else incident in the band (re-emission and the walls' thermal tail),
both area-weighted over `elems`. Dewald et al., PoP 27 (2020), infer F_M ≈ 2 on the
wall through the peak.
"""
function mband_recirculation(s::SpectralSolution, elems::AbstractVector{<:Integer}; lo::Real = 1800.0)
    g = _band_groups(s.edges, lo, Inf)
    t = s.gray.model.target
    direct = s.gray.model.F * (s.mband .* sum(s.mband_shares[g]))
    HM = band_flux(s, lo)
    w = t.elem_area[elems]
    d = sum(direct[elems] .* w)
    d > 0 || return Inf
    return (sum(HM[elems] .* w) - d) / d
end
