# This file solves gray diffuse radiosity. Laser-heated surfaces re-emit x-rays, and walls
# re-emit by albedo.
#
# For element i with albedo αᵢ and laser-to-x-ray conversion ηᵢ, the radiosity (emitted
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

# LU factorizations of I - diag(α) F, most recently used first, keyed on α, and the
# preconditioner of the iterative path of `solve`. The preconditioner is the
# factorization of a nearby system, from an earlier albedo iterate or from another model
# with the same elements, such as the previous time of a moving target.
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

Target plus view factors, which is everything needed to solve for x-ray fluxes given
laser deposition. `kwargs` pass through to `view_factors`. `closed = true`, for elements
(apertures included) that form a closed enclosure, removes discretization closure
error with `enforce_closure!`. View-factor rows sum to one only to about 1e-3, and
multiple bounces at high albedo amplify that leak or excess. At α = 0.925, a 5e-4 shift
in the mean row sum of the wall moves the capsule drive by 0.3% and its P2 by about
0.05%. `closed` defaults to the declaration of the target (which `hohlraum` makes) when
the view factors are computed here, and to `false` for a given `F`, which is used as
is. `closed = false` gives the raw view factors of a closed target.

For a target of several cavities, such as two halfraums across a package or two
hohlraums side by side, `closed` is a collection of cavities. Each cavity is the
elements of one closed enclosure, given as a vector of element indices or of component
or part names. A part name selects every piece of the part, entrance holes included,
which is what the cavity of a hohlraum contains. The rows of each cavity are balanced
to sum to one over its elements (`enforce_closure!(F, A; cavities)`). Elements in no
cavity, such as the outer face of a package or free parts in the open, keep their raw
rows apart from the reciprocal scaling of their entries with cavity elements. One
cavity that contains every element is exactly `closed = true`. The caller declares
which elements form a cavity, and `open_fraction(target; cavity)` checks the
declaration by rays.

The model caches LU factorizations of the radiosity matrix `I - diag(α) F` keyed on
the albedo vector (the two most recently used), so `solve` and `flux_response` with
unchanged albedos skip the factorization. The cache is thread-safe. Do not mutate `F`
after the first solve (or call `Briareus.empty_cache!(model)` after mutating it).

Large models, with more than `Briareus.LARGE_MODEL[]` = 10,000 elements by default,
store the view factors they compute in a [`Float32Matrix`](@ref). Above
`Briareus.HMATRIX_MODEL[]` = 15,000 they use a hierarchical [`HMatrix`](@ref)
(`storage = :auto`; `Float64`, `Float32` or `:hmatrix` select a storage). Large models
are solved by GMRES without factorizations (`solver = :auto`; `:direct` or `:iterative`
select a solver). The storage of a given `F` is kept unless `storage` names another, and an
`HMatrix` is solved only iteratively. Each GMRES iteration is one product with `F`
accumulated in Float64, so iterative solves reach the same relative residual (1e-12) as
the preconditioned ones. Float32 storage changes the capsule moments by about 1e-8
relative (research notes, "Large radiation models"). `RadiationModel(target, F;
solver = :auto)` wraps a given matrix as is.
"""
struct RadiationModel{M<:AbstractMatrix{Float64}}
    target::Target
    F::M
    cache::_FactorCache
    iterative::Bool
    function RadiationModel(t::Target, F::AbstractMatrix; solver::Symbol = :auto)
        solver in (:auto, :direct, :iterative) ||
            throw(ArgumentError("solver must be :auto, :direct or :iterative, got :$solver"))
        F = F isa Union{Float32Matrix,_OperatorMatrix} ? F : convert(Matrix{Float64}, F)
        F isa _OperatorMatrix && solver === :direct &&
            throw(ArgumentError("hierarchical view factors (HMatrix) are solved iteratively only"))
        iterative = solver === :iterative ||
                    (solver === :auto && (F isa Union{Float32Matrix,_OperatorMatrix} || size(F, 1) > LARGE_MODEL[]))
        return new{typeof(F)}(t, F, _FactorCache(), iterative)
    end
end

"""
    Briareus.LARGE_MODEL[]

Element count above which `RadiationModel` stores the view factors it computes in
Float32 (`storage = :auto`) and solves iteratively (`solver = :auto`), 10,000 by
default. A dense Float64 model stores the view factors and up to two LU
factorizations, 24N² bytes (6 GB at 16k elements, 38 GB at 40k). The large path needs
4N² bytes and a Krylov basis (research notes, "Large radiation models").
"""
const LARGE_MODEL = Ref(10_000)

"""
    Briareus.HMATRIX_MODEL[]

Element count above which `RadiationModel` computes hierarchical view factors
([`HMatrix`](@ref), `storage = :auto`), 15,000 by default. At 17.6k elements (the fine
hohlraum preset), hierarchical view factors are assembled in half the time of the
Float32 matrix and use a fifth of its memory. They change the capsule moments by about
1e-8 (research notes, "Hierarchical view factors").
"""
const HMATRIX_MODEL = Ref(15_000)

function RadiationModel(t::Target; F = nothing, closed = nothing,
                        storage = :auto, solver::Symbol = :auto, kwargs...)
    storage in (:auto, Float64, Float32, :hmatrix) ||
        throw(ArgumentError("storage must be :auto, Float64, Float32 or :hmatrix, got $storage"))
    closed = something(closed, F === nothing && t.closed)
    cavities = _cavities(t, closed)
    closed = cavities !== nothing
    if F === nothing
        N = nelements(t)
        st = storage === :auto ? (N > HMATRIX_MODEL[] ? :hmatrix : N > LARGE_MODEL[] ? Float32 : Float64) : storage
        F = view_factors(t; storage = st, kwargs...)
    elseif storage === :hmatrix && !(F isa HMatrix)
        throw(ArgumentError("storage = :hmatrix computes the view factors: pass no F"))
    elseif storage === Float32 && !(F isa Float32Matrix)
        F = Float32Matrix(F)
    elseif storage === Float64 && F isa Union{Float32Matrix,HMatrix}
        F = Matrix{Float64}(F)
    elseif closed
        F = copy(F)                 # balanced below; the matrix of the caller stays unchanged
    end
    closed && enforce_closure!(F, t.elem_area; cavities)
    return RadiationModel(t, F; solver)
end

# Returns the cavities that `closed` declares on `t` (`RadiationModel`). The result is
# `nothing` for none (`false`), `[1:N]` for the whole target (`true`), and otherwise one
# element vector per cavity (from indices, or component or part names). The vectors are
# checked to be nonempty, in range and disjoint.
function _cavities(t::Target, closed)
    closed === false && return nothing
    N = nelements(t)
    closed === true && return [collect(1:N)]
    closed isa Union{Symbol,AbstractString,Integer} &&
        throw(ArgumentError("closed is true, false or a collection of cavities (element indices or names each)"))
    sets = Vector{Int}[]
    for c in closed
        es = c isa AbstractVector{<:Integer} ? collect(Int, c) :
             c isa Symbol ? elements(t, c) : elements(t, (Symbol(n) for n in c)...)
        isempty(es) && throw(ArgumentError("a cavity has no elements"))
        all(e -> 1 <= e <= N, es) || throw(ArgumentError("a cavity's element lies outside 1:$N"))
        push!(sets, es)
    end
    seen = falses(N)
    for es in sets, e in es
        seen[e] && throw(ArgumentError("element $e lies in two cavities (or twice in one)"))
        seen[e] = true
    end
    return sets
end

"""
    enforce_closure!(F, areas; cavities = nothing, tol = 1e-12, maxiter = 1000) -> F

Rescale view factors so every row sums to one while keeping reciprocity
`Aᵢ Fᵢⱼ = Aⱼ Fⱼᵢ`. The exchange matrix `Gᵢⱼ = Aᵢ Fᵢⱼ` becomes `sᵢ Gᵢⱼ sⱼ`, with `s`
from symmetric (Sinkhorn) balancing. Only valid for closed enclosures. Rows with zero
view factors are left untouched.

With `cavities` (disjoint vectors of element indices), each cavity is balanced
separately. The rows of its elements sum to one over its elements, from the Sinkhorn
balancing of its block of `G`, and elements in no cavity keep `sᵢ = 1`. The one scaling
`s` applies to the whole matrix, so reciprocity holds everywhere. An entry between a
cavity and another element changes only by the factors of the cavity; the other element
may belong to another cavity, or lie outside every cavity and be seen through an
aperture or a crack. A single cavity of every element gives exactly the whole-matrix
balancing.
"""
function enforce_closure!(F::AbstractMatrix, A::AbstractVector; cavities = nothing, tol = 1e-12,
                          maxiter::Integer = 1000)
    cavities === nothing || (length(cavities) == 1 && length(only(cavities)) == length(A)) ||
        return _balance!(F, _cavity_scaling(F, A, cavities, tol, maxiter))
    # with G = diag(A) F, (G s)ᵢ = Aᵢ (F s)ᵢ and sᵢ Gᵢⱼ sⱼ / Aᵢ = sᵢ Fᵢⱼ sⱼ, so F is not copied
    s = ones(length(A))
    active = _rowsums(F) .> 0
    r = similar(s)
    for _ in 1:maxiter
        r .= A .* mul!(r, F, s)
        snew = [active[i] ? sqrt(s[i] * A[i] / r[i]) : 1.0 for i in eachindex(s)]
        δ = maximum(abs.(snew .- s))
        s = snew
        δ < tol && break
    end
    return _balance!(F, s)
end

# Returns the Sinkhorn scaling of the block of G = diag(A) F for each cavity
# (`enforce_closure!`, with disjoint element vectors). (G s)ᵢ over the cavity of i comes
# from one product with s masked to that cavity.
function _cavity_scaling(F, A, cavities, tol, maxiter)
    N = length(A)
    s, u, q, r = ones(N), zeros(N), zeros(N), zeros(N)
    active = falses(N)
    for c in cavities
        fill!(u, 0.0)
        u[c] .= 1.0
        mul!(q, F, u)
        active[c] .= q[c] .> 0
    end
    for _ in 1:maxiter
        for c in cavities
            fill!(u, 0.0)
            u[c] .= s[c]
            mul!(q, F, u)
            r[c] .= A[c] .* q[c]
        end
        δ = 0.0
        for c in cavities, i in c
            active[i] || continue
            snew = sqrt(s[i] * A[i] / r[i])
            δ = max(δ, abs(snew - s[i]))
            s[i] = snew
        end
        δ < tol && break
    end
    return s
end

"""
    open_fraction(target; cavity = nothing, rays = 64, seed = 1) -> (; fraction, element)

Whether the facets of the target close an enclosure, apertures counted. The result
gives the area-weighted `fraction` of diffuse rays leaving the elements that meet no
facet at all, and that fraction per `element` (zero on apertures). `rays`
cosine-distributed rays leave from uniform points on each radiating side except those of
the apertures. Apertures stop the rays here, as walls would. A target is closed, and
`closed = true` balancing (`enforce_closure!`) applies, when every opening has an
aperture. The fraction is then zero up to cracks where meshes meet without sharing
vertices, as for a disk aperture of 24 sides in a bore of 20 or 360. Such cracks give
2–7 × 10⁻⁴ on the cavities of "Closure on bored, windowed and CAD cavities" in the
research notes. The holes of `hohlraum` meet its caps at shared vertices, and the
hohlraum does not leak; the fraction was 2.7 × 10⁻³ at h = 1 when the holes had half the
count of the caps. An opening left without an aperture shows as the share of the view
that it takes, which may be no larger: 4 × 10⁻³ for a 40 µm waist gap, and 7 × 10⁻² for
the two 3.1 mm LEHs of a CAD cavity. The fraction therefore indicates where to look
(`element`), not whether to balance, which follows from the apertures given. Curved
surfaces are ignored: cracks between the facets of a curved part and an aperture do not
leak in the view factors, which integrate over the true surfaces.

`cavity` (element indices, or component or part names) checks one cavity of a target
of several, among the cavities that `RadiationModel(...; closed)` balances. The rays
leave only from the elements of the cavity. A ray leaks when it misses every facet, or when
it first meets an element side outside the cavity or a one-sided surface from behind.
An element that the declaration of the cavity misses, or one that belongs to another
cavity, therefore shows as a leak where it is seen. `element` is zero outside the
cavity.

Each element draws from a separate stream (from `seed` and the element), so the result
does not depend on the thread count.
"""
function open_fraction(t::Target; cavity = nothing, rays::Integer = 64, seed::Integer = 1)
    rays > 0 || throw(ArgumentError("rays must be positive"))
    N = nelements(t)
    nt = ntriangles(t)
    bvh = _sah_bvh(BVH(t.tri_v0, t.tri_e1, t.tri_e2, Int32.(1:nt)))   # apertures included
    cum = [cumsum(t.tri_area[t.elem_tris[e]]) for e in 1:N]
    within = cavity === nothing ? nothing : (b = falses(N); b[only(_cavities(t, [cavity]))] .= true; b)
    open = [t.components[t.elem_comp[e]].kind !== :aperture && (within === nothing || within[e]) for e in 1:N]
    miss = zeros(N)
    Threads.@threads for e in 1:N
        open[e] && (miss[e] = _open_rays(t, bvh, cum, e, Int(rays), Random.Xoshiro(1_000_003 * seed + e), within))
    end
    A = t.elem_area
    w = sum(A[open])
    return (; fraction = w > 0 ? sum(miss[e] * A[e] for e in 1:N if open[e]) / w : 0.0, element = miss)
end

# Returns the fraction of `rays` diffuse rays from element e that meet no facet or, with
# `within` (a mask of the elements of a cavity), first meet a side outside the cavity.
function _open_rays(t, bvh, cum, e, rays, rng, within)
    eps = _ray_eps(t)
    n = 0
    for _ in 1:rays
        x, nrm = _mc_point(t, cum, e, rng)
        d = _cosine_direction(nrm, rng)
        _, k = closest_hit(bvh, x + eps * nrm, d)
        if k == 0
            n += 1
        elseif within !== nothing
            h = dot(d, t.tri_normal[k]) < 0 ? t.tri_front[k] : t.tri_back[k]
            (h == 0 || !within[h]) && (n += 1)
        end
    end
    return n / rays
end

function Base.show(io::IO, m::RadiationModel)
    print(io, "RadiationModel(", nelements(m.target), " elements",
          m.F isa Float32Matrix ? ", Float32" : m.F isa _OperatorMatrix ? ", hierarchical" : "", m.iterative ? ", iterative" : "", ")")
end

# Returns the LU factorization of I - diag(α) F, from the cache of the model when α was
# used recently. The factorization runs under the lock, so concurrent solves with the
# same new albedo wait for one factorization and do not repeat it. The returned LU is
# only read.
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
# of F). GMRES on (I - diag(α)F) B = S, right-preconditioned with the LU of such a nearby
# system, reaches a relative residual of 1e-12 in a few iterations. Each iteration is a
# product with F and a pair of triangular solves, 3n² flops against 2n³/3 for the LU.
# When GMRES does not converge in `_GMRES_MAXITER` iterations, a fresh factorization is
# used and becomes the preconditioner.

# Test seams, kept in the solve path on purpose. `_ITERATIVE_SOLVE[] = false` factors
# each albedo iterate afresh. That is the reference against which the preconditioned
# path is checked (test/solve_cache.jl), and no public option reproduces it, since
# `solver = :direct` still reuses the cached LU as preconditioner. `_SOLVE_STATS` counts
# fresh factorizations, iterative solves and GMRES iterations, which the tests use to
# check that caches are reused. The cost is one relaxed atomic add per solve, against at
# least O(n²) work for the solve.
"Use the preconditioned iterative path for albedo iterates (`false`: factor each afresh)."
const _ITERATIVE_SOLVE = Ref(true)
const _GMRES_RTOL = 1e-12
const _GMRES_MAXITER = 30
# counters of fresh factorizations, iterative solves and GMRES iterations
const _SOLVE_STATS = (lu = Threads.Atomic{Int}(0), gmres = Threads.Atomic{Int}(0),
                      iterations = Threads.Atomic{Int}(0))
_reset_solve_stats!() = foreach(a -> a[] = 0, _SOLVE_STATS)

# y = (I - diag(α) F) x
function _radiosity_mul!(y, F, α, x)
    mul!(y, F, x)
    @. y = x - α * y
    return y
end

# Solves (I - diag(α)F) x = b by right-preconditioned GMRES with the LU `P` of a nearby
# system, starting from x0. Returns `nothing` when the true residual does not reach
# rtol‖b‖ within maxiter iterations. The iteration restarts from the last iterate if
# roundoff in the recurrence stops it short of that residual.
_gmres_radiosity(F::Matrix{Float64}, α::Vector{Float64}, P::_LUF, b::Vector{Float64},
                 x0::Vector{Float64}; kwargs...) =
    _gmres((y, x) -> _radiosity_mul!(y, F, α, x), P, b, x0; kwargs...)

# Right-preconditioned restarted GMRES for A x = b, with `op!(y, x)` computing y = A x.
# The preconditioner is the LU `P` of a nearby matrix for the radiosity solves above, or
# none (`P = nothing`) for the thermalized group solves of `spectral_solve`.
_precondition!(P::_LUF, z) = ldiv!(P, z)
_precondition!(::Nothing, z) = z
function _gmres(op!, P::Union{_LUF,Nothing}, b::Vector{Float64}, x0::Vector{Float64};
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
        r .= b .- op!(r, x)
        β = norm(r)
        β <= rtol * bn && return x
        used >= maxiter && return nothing
        m = maxiter - used
        fill!(g, 0.0)
        g[1] = β
        V[:, 1] .= r ./ β
        k = 0
        for j in 1:m
            _precondition!(P, copyto!(z, view(V, :, j)))
            op!(w, z)
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
        _precondition!(P, mul!(z, view(V, :, 1:k), y))
        x .+= z
    end
end

# Returns the solution of (I - diag(α)F) B = S for an albedo iterate. GMRES starts from
# `x0` with the preconditioner of the model when it has one. Otherwise, or when GMRES
# fails, a fresh LU is used and becomes the preconditioner. Dual-number sources are
# solved column by column (value and partials) with the same Float64 operator.
function _iterate_solve(m::RadiationModel, α::Vector{Float64}, S::AbstractVector, x0)
    m.iterative && return _iterative_radiosity(m, α, S, x0)
    if !_ITERATIVE_SOLVE[]
        Threads.atomic_add!(_SOLVE_STATS.lu, 1)
        return lu!(_radiosity_matrix(m.F, α)) \ S
    end
    cols, rebuild = _real_columns(S)
    if cols === nothing                                  # other number types are solved directly
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

# Converts a real vector to a one-column matrix, or a vector of first-order Float64
# duals to the columns (value, partials...), and returns the function that maps such
# columns back.
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

# Unpreconditioned solves of a large model (`m.iterative`), which use only products with
# F and no factorization. The operator I - diag(α)F has its eigenvalues in a disk of
# radius max α about 1, since F is similar to a symmetric matrix of spectral radius ≤ 1.
# The eigenvalues of F cluster at 0, so restarted GMRES converges in a few tens of
# iterations. The columns of a block (several right-hand sides, or the value and partials
# of a dual) are iterated in lockstep, each with a separate Krylov basis, and share one
# product with F per iteration.

const _BLOCK_RESTART = 40
const _BLOCK_MAXITER = 400
# the Krylov bases of a block are kept under this many bytes (columns are batched)
const _BLOCK_BYTES = 1 << 28

# Restarted GMRES for the columns of A X = B, `op!(Y, X)` computing Y = A X for an
# n × k block. Returns `nothing` if the true residual of some column does not reach
# rtol‖bᵢ‖ within `maxiter` iterations.
function _gmres_block(op!, B::Matrix{Float64}, X0::Matrix{Float64}; rtol = _GMRES_RTOL,
                      restart::Integer = _BLOCK_RESTART, maxiter::Integer = _BLOCK_MAXITER)
    n, k = size(B)
    bn = [norm(view(B, :, c)) for c in 1:k]
    X = copy(X0)
    R = similar(B)
    V = Array{Float64,3}(undef, n, restart + 1, k)
    Hs = zeros(restart + 1, restart, k)
    cs, sn = zeros(restart, k), zeros(restart, k)
    g = zeros(restart + 1, k)
    Z, W = similar(B), similar(B)
    steps = zeros(Int, k)
    used = 0
    while true
        op!(R, X)
        R .= B .- R
        active = [c for c in 1:k if norm(view(R, :, c)) > rtol * bn[c]]
        isempty(active) && return X
        used >= maxiter && return nothing
        fill!(g, 0.0)
        fill!(steps, 0)
        for c in active
            β = norm(view(R, :, c))
            g[1, c] = β
            V[:, 1, c] .= view(R, :, c) ./ β
        end
        live = copy(active)
        for j in 1:min(restart, maxiter - used)
            nl = length(live)
            Zl, Wl = view(Z, :, 1:nl), view(W, :, 1:nl)
            for (q, c) in enumerate(live)
                Zl[:, q] .= view(V, :, j, c)
            end
            op!(Wl, Zl)
            Threads.atomic_add!(_SOLVE_STATS.iterations, 1)
            done = falses(nl)
            for (q, c) in enumerate(live)
                w = view(Wl, :, q)
                Vj = view(V, :, 1:j, c)
                h = Vj' * w                      # classical Gram–Schmidt, twice
                mul!(w, Vj, h, -1.0, 1.0)
                h2 = Vj' * w
                mul!(w, Vj, h2, -1.0, 1.0)
                h .+= h2
                hn = norm(w)
                for i in 1:j-1
                    hi = cs[i, c] * h[i] + sn[i, c] * h[i+1]
                    h[i+1] = -sn[i, c] * h[i] + cs[i, c] * h[i+1]
                    h[i] = hi
                end
                ρ = hypot(h[j], hn)
                cs[j, c], sn[j, c] = h[j] / ρ, hn / ρ
                h[j] = ρ
                g[j+1, c] = -sn[j, c] * g[j, c]
                g[j, c] = cs[j, c] * g[j, c]
                Hs[1:j, j, c] .= view(h, 1:j)
                steps[c] = j
                hn > 0 && (V[:, j+1, c] .= w ./ hn)
                done[q] = abs(g[j+1, c]) <= 0.5 * rtol * bn[c] || hn == 0
            end
            live = live[.!done]
            isempty(live) && break
        end
        used += maximum(steps)
        U = view(Z, :, 1:length(active))
        for (q, c) in enumerate(active)
            kc = steps[c]
            y = UpperTriangular(view(Hs, 1:kc, 1:kc, c)) \ g[1:kc, c]
            mul!(view(U, :, q), view(V, :, 1:kc, c), y)
        end
        for (q, c) in enumerate(active)
            X[:, c] .+= view(U, :, q)
        end
    end
end

# Solves op X = B one batch of columns at a time, with the Krylov bases under
# `_BLOCK_BYTES`.
function _solve_columns(op!, B::Matrix{Float64}, X0::Matrix{Float64}, what)
    n, k = size(B)
    kb = max(1, _BLOCK_BYTES ÷ (8 * n * (_BLOCK_RESTART + 1)))
    X = similar(B)
    for c0 in 1:kb:k
        cols = c0:min(c0 + kb - 1, k)
        Threads.atomic_add!(_SOLVE_STATS.gmres, length(cols))
        x = _gmres_block(op!, B[:, cols], X0[:, cols])
        x === nothing && error("iterative radiosity solve ($what) did not converge in $_BLOCK_MAXITER iterations")
        X[:, cols] .= x
    end
    return X
end

# (I - diag(α)F) X and its transpose Xᵀ ↦ X - Fᵀ(α ∘ X), on blocks.
function _radiosity_block!(Y, F, α, X)
    mul!(Y, F, X)
    Y .= X .- α .* Y
    return Y
end
function _radiosity_tblock!(Y, F, α, X)
    mul!(Y, transpose(F), α .* X)
    Y .= X .- Y
    return Y
end

# B = (I - diag(α)F)⁻¹ S on a large model, from `x0` (an earlier radiosity) or from S;
# Float64 or first-order dual sources (value and partials as columns).
function _iterative_radiosity(m::RadiationModel, α::Vector{Float64}, S::AbstractVector, x0)
    cols, rebuild = _real_columns(S)
    cols === nothing &&
        throw(ArgumentError("iterative radiosity solves take Float64 or first-order dual sources, got $(eltype(S))"))
    X0 = copy(cols)
    if x0 !== nothing && length(x0) == length(S)
        c0 = first(_real_columns(x0))
        if c0 !== nothing
            nc = min(size(c0, 2), size(X0, 2))
            X0[:, 1:nc] .= view(c0, :, 1:nc)
        end
    end
    F = m.F
    return rebuild(_solve_columns((Y, X) -> _radiosity_block!(Y, F, α, X), cols, X0, "radiosity"))
end

# ... and for the columns of a matrix of sources
function _iterative_radiosity(m::RadiationModel, α::Vector{Float64}, S::Matrix{Float64}, ::Nothing)
    F = m.F
    return _solve_columns((Y, X) -> _radiosity_block!(Y, F, α, X), S, copy(S), "radiosity")
end

# X = (I - diag(α)F)⁻ᵀ R on a large model.
function _iterative_transposed(m::RadiationModel, α::Vector{Float64}, R::Matrix{Float64})
    F = m.F
    return _solve_columns((Y, X) -> _radiosity_tblock!(Y, F, α, X), R, copy(R), "flux response")
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
        m.cache.precond[] = nothing
    end
    return m
end

"""
    albedos(target, time = 0.0, H = nothing; burnthrough = nothing) -> Vector

Albedo per element at `time`. Temperature-dependent models use the radiation
temperature of each element from incident flux `H`, or their reference temperature when
`H` is `nothing`. `burnthrough` (ns, per element; `NaN` where not given) sets the
burnthrough time of `BurnthroughAlbedo` elements.
"""
albedos(t::Target, time = 0.0, H = nothing; burnthrough = nothing) =
    first(_albedos(t, time, H, burnthrough))

# Returns the albedo per element and, alongside, the albedo that the element would have
# intact. The intact albedo is the `before` albedo at the same temperature for a
# thin-foil element, the closed albedo for a closing gap, and the albedo itself
# elsewhere. The third output is the fraction of the incident flux that the element
# transmits out of the enclosure, which is zero except on a burned-through thin foil or
# an open gap.
#
# The fourth output is the emission per element (TW/mm²) that the state of a model adds
# to the x-ray source. Such an element is a wall-loss element of a `drive_history`
# re-emitting stored energy (`_WallHistory`), or a thin foil whose `before` follows such
# a history (`_FoilHistory`). The emission is zero elsewhere, and the output is `nothing`
# when no element has an emission.
function _albedos(t::Target, time, H, burnthrough)
    n = nelements(t)
    α, intact, transmission = zeros(n), zeros(n), zeros(n)
    emission = nothing
    for e in 1:n
        is_aperture(t, e) && !_transmits(element_material(t, e).albedo) && continue
        a = element_material(t, e).albedo
        T = H === nothing ? nothing : radiation_temperature(ForwardDiff.value(H[e]))
        state = _element_value(burnthrough, e)
        wall = _element_wall(burnthrough, e)
        if a isa WallLossAlbedo && wall !== nothing
            α[e], E, _ = _wall_split(a, wall, T)
            intact[e] = α[e]
            if E != 0
                emission === nothing && (emission = zeros(n))
                emission[e] = E
            end
        elseif a isa BurnthroughAlbedo && wall !== nothing
            # a thin foil whose `before` follows the wall history of the element
            α[e], intact[e], transmission[e], E, _ = _foil_history_split(a, time, T, _FoilHistory(state, wall))
            if E != 0
                emission === nothing && (emission = zeros(n))
                emission[e] = E
            end
        else
            α[e], intact[e], transmission[e] = _albedo_split(a, time, T, state)
        end
    end
    return α, intact, transmission, emission
end
"""
    conversions(target, time = 0.0, laser = nothing) -> Vector

Laser-to-x-ray conversion efficiency per element (0 on apertures) at `time` (ns), for
the absorbed laser power `laser` (TW per element). Intensity-dependent
[`Conversion`](@ref) models need `laser`, and they use the absorbed flux `laser/area`.
Generic in the number type of `laser`.
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

# Returns the x-ray source flux ηP/A per element at `time`. With constant conversions it
# is exactly the product of the numbers of the materials.
_xray_source(t::Target, P, time) = _constant_conversion(t) ? conversions(t) .* P ./ t.elem_area :
                                   conversions(t, time, P) .* P ./ t.elem_area

_temperature_dependent(t::Target) =
    any(c -> _temperature_dependent(c.material.albedo)::Bool, t.components)
_temperature_dependent(::AlbedoModel) = false
_temperature_dependent(::WallLossAlbedo) = true
_temperature_dependent(::TemperatureAlbedo) = true
_temperature_dependent(m::CapsuleAlbedo) = _temperature_dependent(m.albedo)
_temperature_dependent(m::BurnthroughAlbedo) =
    _temperature_dependent(m.before) || m.transmission isa Union{RosselandTransmission,DiffusionTransmission} ||
    (m.after !== nothing && _temperature_dependent(m.after))
_temperature_dependent(m::GapClosure) = _temperature_dependent(m.after)
_temperature_dependent(m::LipClosure) = _temperature_dependent(m.after)

"""
    RadiationSolution

Fields:

- `laser`: absorbed laser power per element (TW), of every color with the `unconverted`
  of `solve`;
- `source`: x-ray source flux ηP/A, plus the emission of wall elements re-emitting
  stored energy in a `drive_history` step;
- `B` (radiosity) and `H` (incident flux); `source`, `B` and `H` are per element in
  TW/mm²;
- `albedo`: per element, the albedos with which the radiosity was solved, so that
  `B = source + albedo .* H` exactly;
- `intact`: the albedo that each element would have intact, which is the `before` albedo
  of a thin-foil element and the albedo elsewhere;
- `transmission`: the fraction of the incident flux that each element transmits out of
  the enclosure, nonzero only on a burned-through `BurnthroughAlbedo` foil (see
  `energy_balance`);
- `time`;
- `converted`: the part of `source` converted from laser light (TW/mm²), from every
  color and without the emission of a re-emitting wall;
- `laser3`: the 3ω part of `laser` (TW per element; `laser` itself without
  `unconverted`), from which [`spectral_solve`](@ref) forms its default M-band.

Without `converted`, the constructor takes it as the 3ω conversion of `laser`, `ηP/A`,
and without `laser3` it takes `laser` as all 3ω. `iterations`, `residual` and
`converged` report the albedo fixed point of a temperature-dependent albedo (`solve`).
They are the number of radiosity solves, the largest albedo change after the last solve
(an emission change counted as the albedo change it equals), and whether that change
fell below `tol`. An unconverged solution is the radiosity of the last iterate, with the
albedos with which it was solved; it is not self-consistent. A solution without a fixed
point (explicit or constant albedos, or one built directly) has `iterations = 1`,
`residual = 0` and `converged = true`.
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
    converted::Vector{T}
    laser3::Vector{T}
    iterations::Int
    residual::Float64
    converged::Bool
end
RadiationSolution(m, laser, source, B, H, albedo, intact, transmission, time) =
    RadiationSolution(m, laser, source, B, H, albedo, intact, transmission, time,
                      convert(Vector{eltype(source)}, _xray_source(m.target, laser, time)))
RadiationSolution(m, laser, source, B, H, albedo, intact, transmission, time, converted) =
    RadiationSolution(m, laser, source, B, H, albedo, intact, transmission, time, converted, laser)
RadiationSolution(m, laser, source, B, H, albedo, intact, transmission, time, converted, laser3) =
    RadiationSolution(m, laser, source, B, H, albedo, intact, transmission, time, converted,
                      laser3, 1, 0.0, true)

"""
    solve(model, dep::Deposition; time = 0.0, albedo = nothing, burnthrough = nothing,
          tol = 1e-4, maxiter = 50, guess = nothing, unconverted = nothing)
    solve(model, laser_power::AbstractVector; ...)

Solve the radiosity system for the given absorbed laser power per element. Albedos come
from the materials at `time` unless an explicit per-element `albedo` vector is passed.
Temperature-dependent albedos (`WallLossAlbedo`) are iterated to self-consistency with
the incident-flux temperature of each element (to `tol` in albedo), in at most `maxiter`
radiosity solves. The first ten iterates use the new albedos in full, and so does every
later iterate until the iteration oscillates. An albedo that falls steeply enough with
`T` overshoots, so that the albedo change reverses direction (a negative dot product
with the previous change) and shrinks by less than half. Each such iterate halves the
step from then on, `α ← α + ω (α_new − α)` with `ω` down to 1/64. Results that converge
in ten iterates are therefore those of the plain iteration. An albedo that rises
steeply with `T` runs away monotonically toward its `range` instead, and relaxation does
not help there. The `iterations`, `residual` and `converged` of the result report the
fixed point ([`RadiationSolution`](@ref)), and an unconverged solve also warns.
`burnthrough` (ns, per element, `NaN` where not given) sets the burnthrough times of
`BurnthroughAlbedo` elements. Without it, a thin-foil element uses the explicit
`burnthrough` of its model, or is intact; `drive_history` tracks the burnthrough from
the flux history instead.

A single solve has no history. It iterates a `WallLossAlbedo` to the present
temperature of each element at `time`, as if the wall had been at that temperature
since `t = 0`, and with no re-emission of stored energy. Partway through a shaped pulse,
that wall differs from the wall found by a [`drive_history`](@ref) of the pulse. On a
NIF-like mock-up the wall albedo is 0.74 against 0.67 at the picket, 0.81 against 0.86
in the trough and 0.94 against 0.87 at the peak, and the capsule drive is 15–44% off. A
fixed-time scenario of such a phase takes its albedos from the history
([`history_state`](@ref); research notes, "Drive histories: sliver steps and the albedo
of a phase (2026-11-18)").

`unconverted` adds the x-ray source of unconverted (1ω, 2ω) light absorbed by the
target. It is the result of [`unconverted_deposition`](@ref), which contains the
[`UnconvertedColor`](@ref)s with which it was deposited (its `colors`). It may also be
any collection (a `NamedTuple`, or a `Dict` keyed by `Symbol`) with an entry `colors`,
the `UnconvertedColor`s, and, under the name of each color, a `Deposition` or a
per-element absorbed power (TW). The conversions used are those of these colors, so a
deposition and the source made from it share one set of colors. The source flux of
element i is then

    Sᵢ = [ηᵢ(t, Pᵢ/Aᵢ) Pᵢ + Σ_c r_c ηᵢ(t, P_c,ᵢ/Aᵢ) P_c,ᵢ] / Aᵢ,

with ηᵢ the [`conversion_efficiency`](@ref) of its material, `P` the 3ω power, `P_c` the
absorbed power of the color and `r_c` the `conversion` of the color relative to 3ω. The
conversion of each color is evaluated at the absorbed flux of that color, an
*engineering default*. The relative conversions come from gold disks irradiated at one
wavelength at a time, and no measurement of the conversion of overlapping colors was
found (research notes, "Unconverted light: absorption and conversion vs. wavelength").
With a conversion independent of intensity, this is the 3ω solve of the equivalent
power `P + Σ_c r_c P_c`, within rounding. `sol.laser` is then the absorbed power of all
colors, `P + Σ_c P_c`, and `sol.converted` the source above, which
[`spectral_solve`](@ref) takes as the converted share. Without `unconverted` the source
is the unchanged 3ω source. Generic in the number type of the powers of every color.

The factorization of the system for the initial albedo (that of the materials at
`time`, or the explicit `albedo`) is cached on the model. Repeated solves with unchanged
albedos then cost only a pair of triangular solves and a matrix–vector product. Iterates
of a temperature-dependent albedo differ from call to call. They are solved by GMRES
preconditioned with the factorization of a nearby system (that of the first iterate, or
of the `guess` model), to a relative residual of 1e-12. They are factored afresh only
when GMRES does not converge in a few tens of iterations, and the new factorization then
preconditions later iterates. On the fine 12.5k-element hohlraum, a 50-step
`drive_history` ran about nine times faster per step with no fresh factorization (127
GMRES solves, 819 iterations) than with a fresh factorization per iterate (127 in all).
On the standard mesh it ran about 2.4 times faster. Results agree with the factored path
to about 1e-12 eV in Tr and 1e-11 percentage points in P2/P0.

`guess` is a previous solution of the same model or of one with the same elements, such
as the last step of a history or a keyframe of a moving target. The albedo iteration
then starts from its incident-flux temperatures instead of the reference temperatures of
the materials, and its radiosity starts the GMRES of the first iterate. The
preconditioner of its model is used when this model has none. Fewer iterations are
needed when the guess is close. The result then differs from a cold start by less than
the tolerance of the iteration. The first albedo is not cached in that case.

A large model (`model.iterative`, see [`RadiationModel`](@ref)) is never factored.
Every iterate is solved by unpreconditioned GMRES to the same 1e-12 relative residual,
from the previous iterate or from `guess`. Dual-number sources are solved as one block,
with the value and partials sharing each product with `F`.
"""
solve(m::RadiationModel, d::Deposition; kwargs...) = solve(m, d.element; kwargs...)

function solve(m::RadiationModel, P::AbstractVector; time = 0.0, albedo = nothing,
               burnthrough::Union{Nothing,AbstractVector{<:Real},_StepState} = nothing,
               tol = 1e-4, maxiter::Integer = 50, guess = nothing, unconverted = nothing)
    t = m.target
    S = _xray_source(t, P, time)
    laser = collect(P)
    laser3 = laser                          # the 3ω part (all of it without `unconverted`)
    if unconverted !== nothing
        Sc, Pc = _unconverted_source(t, unconverted, time)
        S = S .+ Sc
        laser = laser .+ Pc
    end
    iterate = albedo === nothing && _temperature_dependent(t)
    warm = iterate && guess !== nothing
    warm && nelements(guess.model.target) != nelements(t) &&
        throw(ArgumentError("guess is a solution of a model with other elements"))
    α, intact, tr, E = albedo === nothing ? _albedos(t, time, warm ? guess.H : nothing, burnthrough) :
                       (a = collect(Float64, albedo); (a, a, zero(a), nothing))
    S0 = S
    E === nothing || (S = S0 .+ E)
    B = H = nothing
    δ, δprev, ω, its, Δprev = 0.0, Inf, 1.0, 0, α
    if warm && nelements(guess.model.target) == nelements(t)
        # The factorization of a nearby system preconditions the iterates of this system
        # until it is factored. The nearby system may belong to another model, such as the
        # previous time of a moving target.
        c = m.cache
        lock(c.lock) do
            c.precond[] === nothing && (c.precond[] = guess.model.cache.precond[])
        end
    end
    for it in 1:(iterate ? maxiter : 1)
        its = it
        if it == 1 && !warm && m.iterative
            B = _iterative_radiosity(m, α, S, nothing)
        elseif it == 1 && !warm
            f = _radiosity_lu(m, α)
            iterate && lock(() -> m.cache.precond[] = f, m.cache.lock)
            B = f \ S
        else
            B = _iterate_solve(m, α, S, it == 1 ? guess.B : B)
        end
        H = m.F * B
        iterate || break
        αn, intactn, trn, En = _albedos(t, time, H, burnthrough)
        δ = maximum(abs.(αn .- α))
        # a change of emission counts as the albedo change it equals (ΔE/H)
        if E !== nothing || En !== nothing
            ev(x, e) = x === nothing ? 0.0 : x[e]
            δ = max(δ, maximum(abs(ev(En, e) - ev(E, e)) / max(ForwardDiff.value(H[e]), eps())
                               for e in eachindex(H)))
        end
        # keep the albedo B was solved with, so B = S + αH exactly (also when unconverged)
        (δ < tol || it == maxiter) && break
        # In an overshooting iteration the change reverses and shrinks by less than half,
        # and the iteration is under-relaxed from here on.
        Δ = αn .- α
        it > _UNDAMPED && δ > δprev / 2 && dot(Δ, Δprev) < 0 && (ω = max(ω / 2, 1 / 64))
        δprev, Δprev = δ, Δ
        if ω == 1
            α, intact, tr, E = αn, intactn, trn, En
        else
            α, intact, tr = _relax(α, αn, ω), _relax(intact, intactn, ω), _relax(tr, trn, ω)
            E = _relax(E, En, ω)
        end
        S = E === nothing ? S0 : S0 .+ E
    end
    converged = !iterate || δ < tol
    converged || @warn "albedo iteration did not converge" iterations = its δ tol
    return RadiationSolution(m, laser, S, B, H, α, intact, tr, Float64(time), S0,
                             laser3 === laser ? laser : convert(Vector{eltype(laser)}, laser3),
                             its, iterate ? Float64(δ) : 0.0, converged)
end

_relax(x, xn, ω) = x .+ ω .* (xn .- x)
_relax(::Nothing, ::Nothing, ω) = nothing                 # for emissions, `nothing` means none
_relax(x::AbstractVector, ::Nothing, ω) = (1 - ω) .* x
_relax(::Nothing, xn::AbstractVector, ω) = ω .* xn
# The number of albedo iterates used in full before `solve` may relax an oscillating
# iteration. The clamp of a re-emitting wall switches between iterates and converges
# undamped in a few iterates; relaxation slowed one such step from 6 to 17 solves.
const _UNDAMPED = 10

# Returns the x-ray source flux Σ_c r_c η(t, P_c/A) P_c/A and absorbed power Σ_c P_c of
# the unconverted colors per element, with each color converted at its absorbed flux.
_color_power(d::Deposition) = d.element
_color_power(p::AbstractVector{<:Number}) = p
_color_power(x) =      # (an entry of a NamedTuple read by a runtime name may be `colors`)
    throw(ArgumentError("an unconverted color's power is a Deposition or a per-element vector, not $(typeof(x))"))
# Returns the colors of the `unconverted` argument of `solve`, its entry `colors`, which
# are the `UnconvertedColor`s with which its light was deposited. A check first confirms
# that it contains a power per element for each color and nothing else but `skipped`.
function _unconverted_colors(t::Target, u)
    haskey(u, :colors) ||
        throw(ArgumentError("unconverted needs the colors its light was deposited with: the result of " *
                            "`unconverted_deposition`, or a collection with an entry `colors` " *
                            "(`UnconvertedColor`s) beside each color's power"))
    colors = u[:colors]
    for c in colors
        c isa UnconvertedColor || throw(ArgumentError("unconverted's `colors` must be `UnconvertedColor`s"))
        haskey(u, c.name) || throw(ArgumentError("unconverted has no power for $(c.name)"))
        length(_color_power(u[c.name])) == nelements(t) ||
            throw(DimensionMismatch("unconverted power of $(c.name): one entry per element"))
    end
    for k in keys(u)
        k in (:colors, :skipped) || any(c -> c.name === k, colors) ||
            throw(ArgumentError("unconverted has an entry $(repr(k)), which is not one of its colors"))
    end
    return colors
end
function _unconverted_source(t::Target, u, time)
    n = nelements(t)
    colors = _unconverted_colors(t, u)
    T = mapreduce(c -> eltype(_color_power(u[c.name])), promote_type, colors; init = Float64)
    S = zeros(T, n)
    Pc = zeros(T, n)
    for c in colors
        p = _color_power(u[c.name])
        S .+= c.conversion .* _xray_source(t, p, time)
        Pc .+= p
    end
    return S, Pc
end

"""
    flux_response(model, observed; time = 0.0, albedo = nothing, weights = nothing,
                  laser = nothing, unconverted = nothing) -> Matrix

Linear operator `R` with `H[observed] = R * P`, the incident x-ray flux on the observed
elements (TW/mm²) per unit absorbed laser power on every element (TW). It is fixed for
a given geometry and albedo, so optimizers that only change the laser source apply it
as a matrix product instead of solving again. The converged `sol.albedo` of a reference
solution, passed as `albedo`, linearizes about a temperature-dependent albedo state. The
factorization cache of the model is shared with `solve`.

With `weights` (a matrix with one column per observed element, e.g. `legendre_weights`),
`weights * R` is returned directly, and only its rows are solved for. This is far
cheaper than `R` when only a few moments of the observed flux are wanted. On a large
model (`model.iterative`) the rows are solved by GMRES on the transposed system, in
blocks that share each product with `F`, and the full response costs a solve per
observed element.

A time-dependent [`Conversion`](@ref) is evaluated at `time`. An intensity-dependent
conversion makes the source nonlinear in `P`. The absorbed laser power `laser` (TW per
element) then sets the state about which the response is linearized: the conversion of
each element is fixed at its value for that deposition (the secant, exact at `laser`
itself).

With `unconverted` (as in [`solve`](@ref), with its `colors` and, per color, a
`Deposition` or a per-element absorbed power), the result is a `Dict` of responses keyed
by color. `:omega3` is the 3ω response above. For each unconverted color `c`, the
response to the absorbed power of that color is the 3ω response with the conversion
`η(t, P/A)` of each element replaced by `r_c η(t, P_c/A)`, where `r_c` is the
`conversion` of the color. Then `H[observed] ≈ R[:omega3] * P + Σ_c R[c] * P_c`. The
conversion of each color is fixed at its value for the power of that color in
`unconverted`, a secant per color. The sum then equals the per-color result of `solve`
within rounding at those powers, and at any powers for conversions independent of
intensity. With an intensity-dependent conversion, other powers `P_c′` leave the error
of the source change `r_c [η(P_c/A) − η(P_c′/A)] P_c′/A`, carried through the radiosity
solve. On the gold table of Dewald (`dη/dlog₁₀I = −0.15` between 10¹⁴ and
10¹⁵ W/cm²), the secant source of the color exceeds the per-color source by `0.065 δ/η`
of it for a small relative rise `δ` of the flux of the color. The two are equal where
both fluxes lie beyond the same end of the table, which is held at its end values, as
for unconverted light below 10¹⁴ W/cm².
"""
function flux_response(m::RadiationModel, observed::AbstractVector{<:Integer}; time = 0.0,
                       albedo = nothing, weights = nothing, laser = nothing,
                       unconverted = nothing)
    t = m.target
    α = albedo === nothing ? albedos(t, time) : collect(Float64, albedo)
    η = _constant_conversion(t) ? conversions(t) :
        conversions(t, time, laser === nothing ? nothing : ForwardDiff.value.(laser))
    colors = unconverted === nothing ? nothing : _unconverted_colors(t, unconverted)
    # R = F[obs,:] M⁻¹ diag(η ./ A)  ⇒  Rᵀ = diag(η ./ A) M⁻ᵀ F[obs,:]ᵀ
    rhs = if weights === nothing
        _rows_transposed(m.F, observed)
    else
        # (weights F[obs,:])ᵀ as one BLAS product over all elements
        Wn = zeros(size(weights, 1), nelements(t))
        Wn[:, observed] .= weights
        Matrix(transpose(Wn * m.F))
    end
    X = m.iterative ? _iterative_transposed(m, α, rhs) : ldiv!(_radiosity_lu(m, α)', rhs)
    unconverted === nothing && return Matrix(transpose(X .* (η ./ t.elem_area)))
    R = Dict{Symbol,Matrix{Float64}}(:omega3 => Matrix(transpose(X .* (η ./ t.elem_area))))
    for c in colors
        ηc = _constant_conversion(t) ? η :
             conversions(t, time, Float64[ForwardDiff.value(x) for x in _color_power(unconverted[c.name])])
        R[c.name] = Matrix(transpose(X .* (c.conversion .* ηc ./ t.elem_area)))
    end
    return R
end

# F[observed, :]ᵀ (an HMatrix by products with unit vectors)
_rows_transposed(F::AbstractMatrix, observed) = Matrix(transpose(view(F, observed, :)))
function _rows_transposed(F::_OperatorMatrix, observed)
    E = zeros(size(F, 1), length(observed))
    for (k, i) in enumerate(observed)
        E[i, k] = 1.0
    end
    return transpose(F) * E
end

"""
    energy_balance(sol; by = :part) -> NamedTuple

Global x-ray power accounting (TW). The terms are `source` (Σ ηPᵢ), `absorbed` by each
component (`(1-α)HA`), `escaped` through each aperture, `transmitted` through each
burned-through thin-foil component (`BurnthroughAlbedo`) and each open or closing gap
(`GapClosure`, `(1-s)HA` while a fraction `s` is closed), and `lost` to open space,
which is the power that leaves the geometry without reaching any element. The loss
`(1-α)HA` of a thin-foil component splits into the transmitted part `τHA`, with `τ` its
`sol.transmission`, and the absorbed part `(1-α-τ)HA`. With the default `transmission`
`T`, `τ = sT` and `α = (1-sT)α_before`, with `s` the progress of the transition, so the
foil absorbs `(1-sT)(1-α_before)HA`. With an explicit `after`, `τ = α_before-α`, and the
foil absorbs `(1-α_before)HA`, as if intact. A capsule with a `CapsuleAlbedo` absorbs
`(1-α_c)(R/R₀)² HA`, the share of its effective sphere, plus its `outside` share
`a_out (1-(R/R₀)²) HA`. This is the power of ``E_{cap}``, which `energy_budget`
integrates over a `drive_history`. `source ≈ Σabsorbed + Σescaped + Σtransmitted + lost`.

The dictionaries are keyed by part: the pieces of a part split for computation
(`curved_components`) are summed into it. `by = :component` keys them by component.
"""
function energy_balance(s::RadiationSolution; by::Symbol = :part)
    by in (:part, :component) || throw(ArgumentError("by must be :part or :component"))
    m = s.model
    t = m.target
    α, τ = s.albedo, s.transmission
    absorbed = Dict{Symbol,Float64}()
    escaped = Dict{Symbol,Float64}()
    transmitted = Dict{Symbol,Float64}()
    add!(d, k, v) = (d[k] = get(d, k, 0.0) + v)
    for (ci, c) in enumerate(t.components)
        es = t.comp_elems[ci]
        n = by === :part ? c.part : c.name
        P(e) = ForwardDiff.value(s.H[e]) * t.elem_area[e]
        if c.kind === :aperture && _transmits(c.material.albedo)     # a closing lip
            add!(escaped, n, sum(τ[e] * P(e) for e in es; init = 0.0))
            add!(absorbed, n, sum((1 - α[e] - τ[e]) * P(e) for e in es; init = 0.0))
        elseif c.kind === :aperture
            add!(escaped, n, sum(P, es; init = 0.0))
        elseif _transmits(c.material.albedo)
            add!(absorbed, n, sum((1 - α[e] - τ[e]) * P(e) for e in es; init = 0.0))
            add!(transmitted, n, sum(τ[e] * P(e) for e in es; init = 0.0))
        else
            add!(absorbed, n, sum((1 - α[e]) * P(e) for e in es; init = 0.0))
        end
    end
    lost = sum(ForwardDiff.value(s.B[e]) * t.elem_area[e] * c for (e, c) in
               enumerate(closure_error(m.F)); init = 0.0)
    source = sum(ForwardDiff.value(s.source[e]) * t.elem_area[e] for e in 1:nelements(t))
    return (; source, absorbed, escaped, transmitted, lost)
end

"""
    wall_temperature(sol::RadiationSolution, names...) -> Float64

The radiation temperature (eV) of the area-weighted mean incident flux `ΣHᵢAᵢ / ΣAᵢ` on
the elements of the components or parts `names`. When `names` are the wall and caps of a
cavity, it is the wall temperature of the cavity. With several cavities in one target
(two halfraums across a package), the wall pieces of each cavity give a separate
temperature. A part name selects every component of the part, entrance holes included,
so the wall of a hohlraum is named piece by piece. The `cavity_Tr` of `drive_history`
records it per step.
"""
wall_temperature(sol::RadiationSolution, names::Symbol...) =
    _wall_temperature(sol, elements(sol.model.target, names...))

function _wall_temperature(sol::RadiationSolution, es::AbstractVector{<:Integer})
    isempty(es) && throw(ArgumentError("no wall elements"))
    A = sol.model.target.elem_area[es]
    return radiation_temperature(sum(ForwardDiff.value.(sol.H[es]) .* A) / sum(A))
end

# Returns the wall elements of each cavity of `cavities` as `name => elements`. The
# cavities are pairs `name => wall names`, with a name or a collection of names, and the
# cavity names must be unique.
function _cavity_elements(t::Target, cavities)
    out = Pair{Symbol,Vector{Int}}[]
    for c in cavities
        c isa Pair || throw(ArgumentError("a cavity is a pair name => wall names"))
        ns = last(c) isa Symbol ? (last(c),) : Tuple(last(c))
        push!(out, Symbol(first(c)) => elements(t, ns...))
    end
    allunique(first.(out)) || throw(ArgumentError("cavity names must be unique"))
    return out
end
