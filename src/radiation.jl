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

const _LUF = LU{Float64,Matrix{Float64},Vector{Int}}

# LU factorizations of I - diag(α) F, most recently used first, keyed on α; and the
# preconditioner of `solve`'s iterative path: the factorization of a nearby system (an
# earlier albedo iterate, or another model with the same elements, e.g. the previous
# time of a moving target).
struct _FactorCache
    lock::ReentrantLock
    keys::Vector{Vector{Float64}}
    facts::Vector{_LUF}
    precond::Base.RefValue{Union{Nothing,_LUF}}
end
_FactorCache() = _FactorCache(ReentrantLock(), Vector{Float64}[], _LUF[],
                              Ref{Union{Nothing,_LUF}}(nothing))

"Number of radiosity factorizations a `RadiationModel` keeps (each is n² Float64)."
const FACTOR_CACHE_SIZE = 2

"""
    RadiationModel(target; F = view_factors(target), closed, kwargs...)

Target plus view factors: everything needed to solve for x-ray fluxes given laser
deposition. `kwargs` pass through to `view_factors`. `closed = true`, for elements
(apertures included) that form a closed enclosure, removes discretization closure
error with `enforce_closure!`: view-factor rows sum to one only to ~1e-3, and multiple
bounces at high albedo amplify that leak or excess (at α = 0.925, a 5e-4 shift in the
wall's mean row sum moves the capsule drive by 0.3% and its P2 by ~0.05%). It defaults
to the target's own declaration (`hohlraum` makes it) when the view factors are
computed here, and to `false` for a given `F`, which is used as is; pass
`closed = false` for the raw view factors of a closed target.

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

function RadiationModel(t::Target; F = nothing, closed::Union{Bool,Nothing} = nothing, kwargs...)
    closed = something(closed, F === nothing && t.closed)
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
    # with G = diag(A) F, (G s)ᵢ = Aᵢ (F s)ᵢ and sᵢ Gᵢⱼ sⱼ / Aᵢ = sᵢ Fᵢⱼ sⱼ: no copy of F
    s = ones(length(A))
    active = vec(sum(F; dims = 2)) .> 0
    r = similar(s)
    for _ in 1:maxiter
        r .= A .* mul!(r, F, s)
        snew = [active[i] ? sqrt(s[i] * A[i] / r[i]) : 1.0 for i in eachindex(s)]
        δ = maximum(abs.(snew .- s))
        s = snew
        δ < tol && break
    end
    F .= s .* F .* s'
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
        f = lu!(_radiosity_matrix(m.F, α))
        pushfirst!(c.keys, copy(α))
        pushfirst!(c.facts, f)
        if length(c.keys) > FACTOR_CACHE_SIZE
            pop!(c.keys)
            pop!(c.facts)
        end
        return f
    end
end

# I - diag(α) F, in one allocation.
function _radiosity_matrix(F::Matrix{Float64}, α::Vector{Float64})
    M = F .* (.-α)
    for i in axes(M, 1)
        M[i, i] += 1
    end
    return M
end

# Iterative radiosity solves. An albedo iterate differs from the last by a small δα, and
# a history step from the previous one by a small change of α (and, on a moving target,
# of F): GMRES on (I - diag(α)F) B = S, right-preconditioned with the LU of such a nearby
# system, reaches a relative residual of 1e-12 in a few iterations, each a product with
# F and a pair of triangular solves (3n² flops against the LU's 2n³/3). It falls back to
# a fresh factorization, which becomes the preconditioner, when it does not converge in
# `_GMRES_MAXITER` iterations.

"Use the preconditioned iterative path for albedo iterates (`false`: factor each afresh)."
const _ITERATIVE_SOLVE = Ref(true)
const _GMRES_RTOL = 1e-12
const _GMRES_MAXITER = 30
# counters for benchmarks: fresh factorizations, iterative solves, GMRES iterations
const _SOLVE_STATS = (lu = Threads.Atomic{Int}(0), gmres = Threads.Atomic{Int}(0),
                      iterations = Threads.Atomic{Int}(0))
_reset_solve_stats!() = foreach(a -> a[] = 0, _SOLVE_STATS)

# y = (I - diag(α) F) x
function _radiosity_mul!(y, F, α, x)
    mul!(y, F, x)
    @. y = x - α * y
    return y
end

# Solve (I - diag(α)F) x = b by right-preconditioned GMRES with the LU `P` of a nearby
# system, starting from x0. Returns `nothing` when the true residual does not reach
# rtol‖b‖ within maxiter iterations (restarting from the last iterate if roundoff in
# the recurrence stops short of it).
function _gmres_radiosity(F::Matrix{Float64}, α::Vector{Float64}, P::_LUF,
                          b::Vector{Float64}, x0::Vector{Float64};
                          rtol = _GMRES_RTOL, maxiter::Integer = _GMRES_MAXITER)
    n = length(b)
    bn = norm(b)
    bn == 0 && return zeros(n)
    x = copy(x0)
    r = similar(b)
    z = similar(b)
    w = similar(b)
    V = Matrix{Float64}(undef, n, maxiter + 1)
    R = zeros(maxiter + 1, maxiter)
    cs, sn = zeros(maxiter), zeros(maxiter)
    g = zeros(maxiter + 1)
    used = 0
    while true
        r .= b .- _radiosity_mul!(r, F, α, x)
        β = norm(r)
        β <= rtol * bn && return x
        used >= maxiter && return nothing
        m = maxiter - used
        fill!(g, 0.0)
        g[1] = β
        V[:, 1] .= r ./ β
        k = 0
        for j in 1:m
            ldiv!(P, copyto!(z, view(V, :, j)))
            _radiosity_mul!(w, F, α, z)
            # classical Gram–Schmidt, twice (BLAS-2 and as stable as modified)
            Vj = view(V, :, 1:j)
            h = Vj' * w
            mul!(w, Vj, h, -1.0, 1.0)
            h2 = Vj' * w
            mul!(w, Vj, h2, -1.0, 1.0)
            h .+= h2
            hn = norm(w)
            for i in 1:j-1                       # apply the earlier rotations
                hi = cs[i] * h[i] + sn[i] * h[i+1]
                h[i+1] = -sn[i] * h[i] + cs[i] * h[i+1]
                h[i] = hi
            end
            ρ = hypot(h[j], hn)
            cs[j], sn[j] = h[j] / ρ, hn / ρ
            h[j] = ρ
            g[j+1] = -sn[j] * g[j]
            g[j] = cs[j] * g[j]
            R[1:j, j] .= view(h, 1:j)
            k = j
            hn > 0 && (V[:, j+1] .= w ./ hn)
            (abs(g[j+1]) <= 0.5 * rtol * bn || hn == 0) && break
        end
        used += k
        Threads.atomic_add!(_SOLVE_STATS.iterations, k)
        y = UpperTriangular(view(R, 1:k, 1:k)) \ g[1:k]
        ldiv!(P, mul!(z, view(V, :, 1:k), y))
        x .+= z
    end
end

# Solution of (I - diag(α)F) B = S for an albedo iterate: GMRES from `x0` with the
# model's preconditioner when it has one, else (or when GMRES fails) a fresh LU, which
# becomes the preconditioner. Dual-number sources are solved column by column (value and
# partials) with the same Float64 operator.
function _iterate_solve(m::RadiationModel, α::Vector{Float64}, S::AbstractVector, x0)
    if !_ITERATIVE_SOLVE[]
        Threads.atomic_add!(_SOLVE_STATS.lu, 1)
        return lu!(_radiosity_matrix(m.F, α)) \ S
    end
    cols, rebuild = _real_columns(S)
    if cols === nothing                                  # other number types: direct
        Threads.atomic_add!(_SOLVE_STATS.lu, 1)
        return lu!(_radiosity_matrix(m.F, α)) \ S
    end
    x0c = x0 === nothing || length(x0) != length(S) ? nothing : first(_real_columns(x0))
    P = lock(() -> m.cache.precond[], m.cache.lock)
    X = similar(cols)
    if P !== nothing
        ok = true
        for c in axes(cols, 2)
            Threads.atomic_add!(_SOLVE_STATS.gmres, 1)
            start = x0c !== nothing && c <= size(x0c, 2) ? x0c[:, c] : P \ cols[:, c]
            x = _gmres_radiosity(m.F, α, P, cols[:, c], start)
            x === nothing && (ok = false; break)
            X[:, c] .= x
        end
        ok && return rebuild(X)
    end
    Threads.atomic_add!(_SOLVE_STATS.lu, 1)
    f = lu!(_radiosity_matrix(m.F, α))
    lock(() -> m.cache.precond[] = f, m.cache.lock)
    return rebuild(ldiv!(f, cols))
end

# A real vector as a one-column matrix, or a vector of first-order Float64 duals as the
# columns (value, partials...), with the function that maps such columns back.
_real_columns(S::AbstractVector{Float64}) = (reshape(collect(S), :, 1), X -> vec(X))
function _real_columns(S::AbstractVector{ForwardDiff.Dual{T,Float64,N}}) where {T,N}
    cols = Matrix{Float64}(undef, length(S), N + 1)
    for (i, s) in enumerate(S)
        cols[i, 1] = ForwardDiff.value(s)
        p = ForwardDiff.partials(s)
        for k in 1:N
            cols[i, k+1] = p[k]
        end
    end
    rebuild(X) = [ForwardDiff.Dual{T}(X[i, 1], ntuple(k -> X[i, k+1], Val(N))...)
                  for i in axes(X, 1)]
    return cols, rebuild
end
_real_columns(S::AbstractVector) = (nothing, nothing)

"""
    empty_cache!(model::RadiationModel) -> model

Drop the cached radiosity factorizations (frees memory; required after mutating
`model.F`).
"""
function empty_cache!(m::RadiationModel)
    lock(m.cache.lock) do
        empty!(m.cache.keys)
        empty!(m.cache.facts)
        m.cache.precond[] = nothing
    end
    return m
end

"""
    albedos(target, time = 0.0, H = nothing; burnthrough = nothing) -> Vector

Albedo per element at `time`; temperature-dependent models use each element's
radiation temperature from incident flux `H` (or their reference temperature when `H`
is `nothing`). `burnthrough` (ns, per element; `NaN` where not given) sets the
burnthrough time of `BurnthroughAlbedo` elements.
"""
albedos(t::Target, time = 0.0, H = nothing; burnthrough = nothing) =
    first(_albedos(t, time, H, burnthrough))

# Albedo per element and, alongside, the albedo it would have intact (a thin-foil
# element's `before` albedo at the same temperature, a closing gap's closed albedo, the
# albedo itself elsewhere) and the fraction of its incident flux it transmits out of
# the enclosure (0 but on a burned-through thin foil or an open gap).
function _albedos(t::Target, time, H, burnthrough)
    n = nelements(t)
    α, intact, transmission = zeros(n), zeros(n), zeros(n)
    for e in 1:n
        is_aperture(t, e) && !_transmits(element_material(t, e).albedo) && continue
        a = element_material(t, e).albedo
        T = H === nothing ? nothing : radiation_temperature(ForwardDiff.value(H[e]))
        α[e], intact[e], transmission[e] =
            _albedo_split(a, time, T, burnthrough === nothing ? NaN : burnthrough[e])
    end
    return α, intact, transmission
end
"""
    conversions(target, time = 0.0, laser = nothing) -> Vector

Laser → x-ray conversion efficiency per element (0 on apertures) at `time` (ns), for
the absorbed laser power `laser` (TW per element; needed by intensity-dependent
[`Conversion`](@ref) models, which see the absorbed flux `laser/area`). Generic in the
number type of `laser`.
"""
function conversions(t::Target, time = 0.0, P = nothing)
    _constant_conversion(t) &&
        return Float64[is_aperture(t, e) ? 0.0 : element_material(t, e).conversion::Float64
                       for e in 1:nelements(t)]
    T = P === nothing || !_intensity_dependent(t) ? Float64 : promote_type(Float64, eltype(P))
    η = zeros(T, nelements(t))
    for e in 1:nelements(t)
        is_aperture(t, e) && continue
        m = element_material(t, e).conversion
        if _intensity_dependent(m)
            P === nothing && throw(ArgumentError("the conversion of $(element_material(t, e).name) depends on the laser intensity: pass the absorbed laser power"))
            η[e] = conversion_efficiency(m, time, to_Wcm2(P[e] / t.elem_area[e]))
        else
            η[e] = conversion_efficiency(m, time)
        end
    end
    return η
end
_constant_conversion(t::Target) = all(c -> c.material.conversion isa Real, t.components)
_intensity_dependent(t::Target) = any(c -> _intensity_dependent(c.material.conversion), t.components)

# The x-ray source flux ηP/A per element at `time`: with constant conversions exactly
# the product of the materials' numbers.
_xray_source(t::Target, P, time) = _constant_conversion(t) ? conversions(t) .* P ./ t.elem_area :
                                   conversions(t, time, P) .* P ./ t.elem_area

_temperature_dependent(t::Target) =
    any(c -> _temperature_dependent(c.material.albedo), t.components)
_temperature_dependent(::AlbedoModel) = false
_temperature_dependent(::WallLossAlbedo) = true
_temperature_dependent(::TemperatureAlbedo) = true
_temperature_dependent(m::CapsuleAlbedo) = _temperature_dependent(m.albedo)
_temperature_dependent(m::BurnthroughAlbedo) =
    _temperature_dependent(m.before) || m.transmission isa RosselandTransmission ||
    (m.after !== nothing && _temperature_dependent(m.after))
_temperature_dependent(m::GapClosure) = _temperature_dependent(m.after)

"""
    RadiationSolution

Fields: `laser` (absorbed laser power per element, TW), `source` (x-ray source flux
ηP/A), `B` (radiosity), `H` (incident flux), all per element in TW/mm²; `albedo` per
element (the albedos the radiosity was solved with, so `B = source + albedo .* H`
exactly); `intact`, the albedo each element would have intact (a thin-foil element's
`before` albedo, the albedo elsewhere); `transmission`, the fraction of its incident
flux each element transmits out of the enclosure (nonzero only on a burned-through
`BurnthroughAlbedo` foil; see `energy_balance`); `time`.
"""
struct RadiationSolution{T}
    model::RadiationModel
    laser::Vector{T}
    source::Vector{T}
    B::Vector{T}
    H::Vector{T}
    albedo::Vector{Float64}
    intact::Vector{Float64}
    transmission::Vector{Float64}
    time::Float64
end

"""
    solve(model, dep::Deposition; time = 0.0, albedo = nothing, burnthrough = nothing,
          tol = 1e-4, maxiter = 50, guess = nothing)
    solve(model, laser_power::AbstractVector; ...)

Solve the radiosity system for the given absorbed laser power per element. Albedos come
from the materials at `time` unless an explicit per-element `albedo` vector is passed.
Temperature-dependent albedos (`WallLossAlbedo`) are iterated to self-consistency with
each element's incident-flux temperature (to `tol` in albedo). `burnthrough` (ns, per
element, `NaN` where not given) sets the burnthrough times of `BurnthroughAlbedo`
elements; without it a thin-foil element uses its model's explicit `burnthrough`, or
is intact (`drive_history` tracks it from the flux history instead).

The factorization of the system for the initial albedo (the materials' at `time`, or the
explicit `albedo`) is cached on the model, so repeated solves with unchanged albedos
cost only a pair of triangular solves and a matrix–vector product. Iterates of a
temperature-dependent albedo differ from call to call; they are solved by GMRES
preconditioned with the factorization of a nearby system (the first iterate's, or the
`guess` model's), to a relative residual of 1e-12, and factored afresh only when that
does not converge in a few tens of iterations (the new factorization then preconditions
later iterates). On the fine 12.5k-element hohlraum a 50-step `drive_history` took
32 s per step with a fresh factorization per iterate (127 in all) and takes 3.6 s with
none (127 GMRES solves, 819 iterations); 3.1 → 1.3 s per step on the standard mesh.
Results agree with the factored path to ~1e-12 eV in Tr and ~1e-11 percentage points
in P2/P0.

`guess`, a previous solution of the same model or one with the same elements (e.g. the
last step of a history, or a keyframe of a moving target),
starts the albedo iteration from its incident-flux temperatures instead of the
materials' reference temperatures, and its radiosity starts the first iterate's GMRES
(its model's preconditioner serves when this model has none): fewer iterations when
it is close. The result then
differs from a cold start by less than the iteration's tolerance. The first albedo is
not cached in that case.
"""
solve(m::RadiationModel, d::Deposition; kwargs...) = solve(m, d.element; kwargs...)

function solve(m::RadiationModel, P::AbstractVector; time = 0.0, albedo = nothing,
               burnthrough = nothing, tol = 1e-4, maxiter::Integer = 50, guess = nothing)
    t = m.target
    S = _xray_source(t, P, time)
    iterate = albedo === nothing && _temperature_dependent(t)
    warm = iterate && guess !== nothing
    warm && nelements(guess.model.target) != nelements(t) &&
        throw(ArgumentError("guess is a solution of a model with other elements"))
    α, intact, tr = albedo === nothing ? _albedos(t, time, warm ? guess.H : nothing, burnthrough) :
                    (a = collect(Float64, albedo); (a, a, zero(a)))
    B = H = nothing
    if warm && nelements(guess.model.target) == nelements(t)
        # a nearby system's factorization (another model's: the previous time of a
        # moving target) preconditions this one's iterates until it factors its own
        c = m.cache
        lock(c.lock) do
            c.precond[] === nothing && (c.precond[] = guess.model.cache.precond[])
        end
    end
    for it in 1:(iterate ? maxiter : 1)
        if it == 1 && !warm
            f = _radiosity_lu(m, α)
            iterate && lock(() -> m.cache.precond[] = f, m.cache.lock)
            B = f \ S
        else
            B = _iterate_solve(m, α, S, it == 1 ? guess.B : B)
        end
        H = m.F * B
        iterate || break
        αn, intactn, trn = _albedos(t, time, H, burnthrough)
        δ = maximum(abs.(αn .- α))
        δ < tol && break          # keep the albedo B was solved with, so B = S + αH exactly
        α, intact, tr = αn, intactn, trn
        it == maxiter && @warn "albedo iteration did not converge" δ
    end
    return RadiationSolution(m, collect(P), S, B, H, α, intact, tr, Float64(time))
end

"""
    flux_response(model, observed; time = 0.0, albedo = nothing, weights = nothing) -> Matrix

Linear operator `R` with `H[observed] = R * P`: incident x-ray flux on the observed
elements (TW/mm²) per unit absorbed laser power on every element (TW). It is fixed
for a given geometry and albedo, so optimizers that only change the laser source
apply it as a matrix product instead of re-solving. Pass the converged `sol.albedo` of
a reference solution to linearize about a temperature-dependent albedo state. Shares
the model's factorization cache with `solve`.

With `weights` (a matrix with one column per observed element, e.g. `legendre_weights`),
returns `weights * R` directly, solving for its rows only: far cheaper than `R` when
only a few moments of the observed flux are wanted.

A time-dependent [`Conversion`](@ref) is evaluated at `time`. An intensity-dependent one
makes the source nonlinear in `P`: pass the absorbed laser power `laser` (TW per
element) to linearize about, i.e. to fix each element's conversion at its value for
that deposition (the secant, exact at `laser` itself).
"""
function flux_response(m::RadiationModel, observed::AbstractVector{<:Integer}; time = 0.0,
                       albedo = nothing, weights = nothing, laser = nothing)
    t = m.target
    α = albedo === nothing ? albedos(t, time) : collect(Float64, albedo)
    η = _constant_conversion(t) ? conversions(t) :
        conversions(t, time, laser === nothing ? nothing : ForwardDiff.value.(laser))
    # R = F[obs,:] M⁻¹ diag(η ./ A)  ⇒  Rᵀ = diag(η ./ A) M⁻ᵀ F[obs,:]ᵀ
    rhs = if weights === nothing
        Matrix(transpose(view(m.F, observed, :)))
    else
        # (weights F[obs,:])ᵀ as one BLAS product over all elements
        Wn = zeros(size(weights, 1), nelements(t))
        Wn[:, observed] .= weights
        Matrix(transpose(Wn * m.F))
    end
    X = ldiv!(_radiosity_lu(m, α)', rhs)
    return Matrix(transpose(X .* (η ./ t.elem_area)))
end

"""
    energy_balance(sol) -> NamedTuple

Global x-ray power accounting (TW): `source` (Σ ηPᵢ), `absorbed` by each component
(`(1-α)HA`), `escaped` through each aperture, `transmitted` through each burned-through
thin-foil component (`BurnthroughAlbedo`) and each open or closing gap (`GapClosure`,
`(1-s)HA` while a fraction `s` is closed), and `lost` to open space (what leaves the
geometry without reaching any element). A thin-foil component's loss `(1-α)HA` splits
into what it transmits, `τHA` with `τ` its `sol.transmission`, and what it absorbs,
`(1-α-τ)HA`: with the default `transmission` `T`, `τ = sT` and `α = (1-sT)α_before`
(`s` the transition's progress), so it absorbs `(1-sT)(1-α_before)HA`; with an
explicit `after`, `τ = α_before-α` and it absorbs `(1-α_before)HA`, as if intact.
A capsule with a `CapsuleAlbedo` absorbs `(1-α_c)(R/R₀)² HA`, its effective sphere's
share (E_cap's power; `energy_budget` integrates it over a `drive_history`).
`source ≈ Σabsorbed + Σescaped + Σtransmitted + lost`.
"""
function energy_balance(s::RadiationSolution)
    m = s.model
    t = m.target
    α, τ = s.albedo, s.transmission
    absorbed = Dict{Symbol,Float64}()
    escaped = Dict{Symbol,Float64}()
    transmitted = Dict{Symbol,Float64}()
    for (ci, c) in enumerate(t.components)
        es = t.comp_elems[ci]
        P(e) = ForwardDiff.value(s.H[e]) * t.elem_area[e]
        if c.kind === :aperture && _transmits(c.material.albedo)     # a closing lip
            escaped[c.name] = sum(τ[e] * P(e) for e in es; init = 0.0)
            absorbed[c.name] = sum((1 - α[e] - τ[e]) * P(e) for e in es; init = 0.0)
        elseif c.kind === :aperture
            escaped[c.name] = sum(P, es; init = 0.0)
        elseif _transmits(c.material.albedo)
            absorbed[c.name] = sum((1 - α[e] - τ[e]) * P(e) for e in es; init = 0.0)
            transmitted[c.name] = sum(τ[e] * P(e) for e in es; init = 0.0)
        else
            absorbed[c.name] = sum((1 - α[e]) * P(e) for e in es; init = 0.0)
        end
    end
    lost = sum(ForwardDiff.value(s.B[e]) * t.elem_area[e] * c for (e, c) in
               enumerate(closure_error(m.F)); init = 0.0)
    source = sum(ForwardDiff.value(s.source[e]) * t.elem_area[e] for e in 1:nelements(t))
    return (; source, absorbed, escaped, transmitted, lost)
end
