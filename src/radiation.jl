# Gray diffuse radiosity: laser-heated surfaces re-emit x-rays, walls re-emit by albedo.
#
# For element i with albedo αᵢ and laser→x-ray conversion ηᵢ, the radiosity (emitted
# flux) Bᵢ and incident flux Hᵢ satisfy
#
#     Bᵢ = ηᵢ Pᵢ / Aᵢ + αᵢ Hᵢ,      Hᵢ = Σⱼ Fᵢⱼ Bⱼ,
#
# where Pᵢ is absorbed laser power. Apertures have α = η = 0, so H on an aperture is
# the flux escaping through it.

"Stefan–Boltzmann constant in TW mm⁻² eV⁻⁴ (σT⁴ with T in eV gives TW/mm²)."
const SIGMA_SB = 1.0283e-9

"Convert a flux in TW/mm² to W/cm²."
to_Wcm2(x) = 1e14 * x

"""
    radiation_temperature(H)

Radiation temperature (eV) of incident flux `H` (TW/mm²): `(H/σ)^(1/4)`.
"""
radiation_temperature(H) = (max(H, zero(H)) / SIGMA_SB)^(1 / 4)

# LU factorizations of I - diag(α) F, most recently used first, keyed on α.
struct _FactorCache
    lock::ReentrantLock
    keys::Vector{Vector{Float64}}
    facts::Vector{LU{Float64,Matrix{Float64},Vector{Int}}}
end
_FactorCache() = _FactorCache(ReentrantLock(), Vector{Float64}[],
                              LU{Float64,Matrix{Float64},Vector{Int}}[])

"Number of radiosity factorizations a `RadiationModel` keeps (each is n² Float64)."
const FACTOR_CACHE_SIZE = 2

"""
    RadiationModel(target; F = view_factors(target), closed = false, kwargs...)

Target plus view factors: everything needed to solve for x-ray fluxes given laser
deposition. `kwargs` pass through to `view_factors`. Set `closed = true` when the
elements (apertures included) form a closed enclosure — e.g. a hohlraum with its LEHs
— to remove discretization closure error with `enforce_closure!`; radiation that would
otherwise leak through the mesh is amplified by multiple bounces at high albedo.

The model caches LU factorizations of the radiosity matrix `I - diag(α) F` keyed on
the albedo vector (the two most recently used), so `solve` and `flux_response` with
unchanged albedos skip the factorization. The cache is thread-safe. Do not mutate `F`
after the first solve (or call `Briareus.empty_cache!(model)` if you do).
"""
struct RadiationModel
    target::Target
    F::Matrix{Float64}
    cache::_FactorCache
    RadiationModel(t::Target, F::AbstractMatrix) = new(t, F, _FactorCache())
end

function RadiationModel(t::Target; F = nothing, closed::Bool = false, kwargs...)
    F = F === nothing ? view_factors(t; kwargs...) : copy(F)
    closed && enforce_closure!(F, t.elem_area)
    return RadiationModel(t, F)
end

"""
    enforce_closure!(F, areas; tol = 1e-12, maxiter = 1000) -> F

Rescale view factors so every row sums to one while keeping reciprocity
`Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ`: the exchange matrix `Gᵢⱼ = Aᵢ Fᵢⱼ` becomes `sᵢ Gᵢⱼ sⱼ` with `s`
from symmetric (Sinkhorn) balancing. Only valid for closed enclosures. Rows that see
nothing are left untouched.
"""
function enforce_closure!(F::AbstractMatrix, A::AbstractVector; tol = 1e-12,
                          maxiter::Integer = 1000)
    G = A .* F
    s = ones(length(A))
    active = vec(sum(G; dims = 2)) .> 0
    for _ in 1:maxiter
        r = G * s
        snew = [active[i] ? sqrt(s[i] * A[i] / r[i]) : 1.0 for i in eachindex(s)]
        δ = maximum(abs.(snew .- s))
        s = snew
        δ < tol && break
    end
    F .= (s .* G .* s') ./ A
    return F
end

function Base.show(io::IO, m::RadiationModel)
    print(io, "RadiationModel(", nelements(m.target), " elements)")
end

# LU factorization of I - diag(α) F, from the model's cache when α was seen recently.
# The factorization runs under the lock, so concurrent solves with the same new albedo
# wait for one factorization instead of each repeating it; the returned LU is only read.
function _radiosity_lu(m::RadiationModel, α::Vector{Float64})
    c = m.cache
    lock(c.lock) do
        i = findfirst(==(α), c.keys)
        if i !== nothing
            if i != 1
                pushfirst!(c.keys, popat!(c.keys, i))
                pushfirst!(c.facts, popat!(c.facts, i))
            end
            return c.facts[1]
        end
        f = lu!(I - α .* m.F)
        pushfirst!(c.keys, copy(α))
        pushfirst!(c.facts, f)
        if length(c.keys) > FACTOR_CACHE_SIZE
            pop!(c.keys)
            pop!(c.facts)
        end
        return f
    end
end

"""
    empty_cache!(model::RadiationModel) -> model

Drop the cached radiosity factorizations (frees memory; required after mutating
`model.F`).
"""
function empty_cache!(m::RadiationModel)
    lock(m.cache.lock) do
        empty!(m.cache.keys)
        empty!(m.cache.facts)
    end
    return m
end

"""
    albedos(target, time = 0.0, H = nothing) -> Vector

Albedo per element at `time`; temperature-dependent models use each element's
radiation temperature from incident flux `H` (or their reference temperature when `H`
is `nothing`).
"""
function albedos(t::Target, time = 0.0, H = nothing)
    return [is_aperture(t, e) ? 0.0 :
            albedo(element_material(t, e).albedo, time,
                   H === nothing ? nothing : radiation_temperature(ForwardDiff.value(H[e])))
            for e in 1:nelements(t)]
end
conversions(t::Target) =
    [is_aperture(t, e) ? 0.0 : element_material(t, e).conversion for e in 1:nelements(t)]

_temperature_dependent(t::Target) =
    any(c -> c.kind === :surface && c.material.albedo isa WallLossAlbedo, t.components)

"""
    RadiationSolution

Fields: `laser` (absorbed laser power per element, TW), `source` (x-ray source flux
ηP/A), `B` (radiosity), `H` (incident flux), all per element in TW/mm²; `albedo` per
element; `time`.
"""
struct RadiationSolution{T}
    model::RadiationModel
    laser::Vector{T}
    source::Vector{T}
    B::Vector{T}
    H::Vector{T}
    albedo::Vector{Float64}
    time::Float64
end

"""
    solve(model, dep::Deposition; time = 0.0, albedo = nothing, tol = 1e-4, maxiter = 50)
    solve(model, laser_power::AbstractVector; ...)

Solve the radiosity system for the given absorbed laser power per element. Albedos come
from the materials at `time` unless an explicit per-element `albedo` vector is passed.
Temperature-dependent albedos (`WallLossAlbedo`) are iterated to self-consistency with
each element's incident-flux temperature (to `tol` in albedo).

The factorization of the system for the initial albedo (the materials' at `time`, or the
explicit `albedo`) is cached on the model, so repeated solves with unchanged albedos
cost only a pair of triangular solves and a matrix–vector product. Iterates of a
temperature-dependent albedo differ from call to call and are factored afresh.
"""
solve(m::RadiationModel, d::Deposition; kwargs...) = solve(m, d.element; kwargs...)

function solve(m::RadiationModel, P::AbstractVector; time = 0.0, albedo = nothing,
               tol = 1e-4, maxiter::Integer = 50)
    t = m.target
    S = conversions(t) .* P ./ t.elem_area
    α = albedo === nothing ? albedos(t, time) : collect(Float64, albedo)
    iterate = albedo === nothing && _temperature_dependent(t)
    B = H = nothing
    for it in 1:(iterate ? maxiter : 1)
        B = (it == 1 ? _radiosity_lu(m, α) : lu!(I - α .* m.F)) \ S
        H = m.F * B
        iterate || break
        αn = albedos(t, time, H)
        δ = maximum(abs.(αn .- α))
        α = αn
        δ < tol && break
        it == maxiter && @warn "albedo iteration did not converge" δ
    end
    return RadiationSolution(m, collect(P), S, B, H, α, Float64(time))
end

"""
    flux_response(model, observed; time = 0.0, albedo = nothing) -> Matrix

Linear operator `R` with `H[observed] = R * P`: incident x-ray flux on the observed
elements (TW/mm²) per unit absorbed laser power on every element (TW). It is fixed
for a given geometry and albedo, so optimizers that only change the laser source
apply it as a matrix product instead of re-solving. Pass the converged `sol.albedo` of
a reference solution to linearize about a temperature-dependent albedo state. Shares
the model's factorization cache with `solve`.
"""
function flux_response(m::RadiationModel, observed::AbstractVector{<:Integer}; time = 0.0,
                       albedo = nothing)
    t = m.target
    α = albedo === nothing ? albedos(t, time) : collect(Float64, albedo)
    η = conversions(t)
    # R = F[obs,:] M⁻¹ diag(η ./ A)  ⇒  Rᵀ = diag(η ./ A) M⁻ᵀ F[obs,:]ᵀ
    X = _radiosity_lu(m, α)' \ Matrix(m.F[observed, :]')
    return Matrix((X .* (η ./ t.elem_area))')
end

"""
    energy_balance(sol) -> NamedTuple

Global x-ray power accounting (TW): `source` (Σ ηPᵢ), `absorbed` by each component
(`(1-α)HA`), `escaped` through each aperture, and `lost` to open space (what leaves
the geometry without reaching any element). `source ≈ Σabsorbed + Σescaped + lost`.
"""
function energy_balance(s::RadiationSolution)
    m = s.model
    t = m.target
    α = s.albedo
    absorbed = Dict{Symbol,Float64}()
    escaped = Dict{Symbol,Float64}()
    for (ci, c) in enumerate(t.components)
        es = t.comp_elems[ci]
        if c.kind === :aperture
            escaped[c.name] = sum(ForwardDiff.value(s.H[e]) * t.elem_area[e] for e in es; init = 0.0)
        else
            absorbed[c.name] = sum((1 - α[e]) * ForwardDiff.value(s.H[e]) * t.elem_area[e]
                                   for e in es; init = 0.0)
        end
    end
    lost = sum(ForwardDiff.value(s.B[e]) * t.elem_area[e] * c for (e, c) in
               enumerate(closure_error(m.F)); init = 0.0)
    source = sum(ForwardDiff.value(s.source[e]) * t.elem_area[e] for e in 1:nelements(t))
    return (; source, absorbed, escaped, lost)
end
