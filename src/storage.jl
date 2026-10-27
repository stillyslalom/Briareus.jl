# Single-precision storage of large dense matrices (view factors), read as Float64.

"""
    Float32Matrix(A::AbstractMatrix)
    Float32Matrix(undef, m, n)

A dense matrix stored in Float32 that reads as Float64: indexing returns the stored
value widened, and products with Float64 vectors and matrices (`*`, `mul!`, also
transposed) accumulate in Float64. A product is therefore exactly the Float64 product
with the matrix `Float64.(A.data)`; the only error against a Float64 original is the
rounding of each entry (relative ≤ 2⁻²⁴ ≈ 6e-8). It halves the memory of a view-factor
matrix: `view_factors(target; storage = Float32)`, and `RadiationModel` by default for
more than `Briareus.LARGE_MODEL[]` elements (research notes, "Large radiation models").

Matrix–vector products are threaded over row (or, transposed, column) blocks and
memory-bound; matrix–matrix products widen column panels of the matrix to Float64 and
hand them to BLAS.
"""
struct Float32Matrix <: AbstractMatrix{Float64}
    data::Matrix{Float32}
end
Float32Matrix(A::AbstractMatrix) = Float32Matrix(Matrix{Float32}(A))
Float32Matrix(::UndefInitializer, m::Integer, n::Integer) =
    Float32Matrix(Matrix{Float32}(undef, m, n))
Float32Matrix(A::Float32Matrix) = A

Base.size(A::Float32Matrix) = size(A.data)
Base.IndexStyle(::Type{Float32Matrix}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(A::Float32Matrix, i::Int) = Float64(A.data[i])
Base.@propagate_inbounds Base.setindex!(A::Float32Matrix, v, i::Int) = (A.data[i] = v; A)
Base.similar(A::Float32Matrix) = Float32Matrix(similar(A.data))
Base.copy(A::Float32Matrix) = Float32Matrix(copy(A.data))
Base.zero(A::Float32Matrix) = Float32Matrix(zero(A.data))
Base.:(==)(A::Float32Matrix, B::Float32Matrix) = A.data == B.data

# Row sums (Float64 accumulation).
_rowsums(F::AbstractMatrix) = vec(sum(F; dims = 2))
_rowsums(F::Float32Matrix) = F * ones(size(F, 2))

# Blocks of `1:n` for `k` workers.
_blocks(n, k) = (k = clamp(k, 1, max(n, 1)); [((b - 1) * n ÷ k + 1):(b * n ÷ k) for b in 1:k])

# y = D x (β = 0), Float64 accumulation; threaded over row blocks.
function _f32_mul!(y::AbstractVector{Float64}, D::Matrix{Float32}, x::AbstractVector{Float64})
    n = size(D, 1)
    Threads.@threads for r in _blocks(n, n < 2048 ? 1 : Threads.nthreads())
        _f32_rows!(y, D, x, r)
    end
    return y
end
function _f32_rows!(y, D, x, r)
    @inbounds for i in r
        y[i] = 0.0
    end
    @inbounds for j in axes(D, 2)
        xj = x[j]
        xj == 0 && continue
        @simd for i in r
            y[i] = muladd(Float64(D[i, j]), xj, y[i])
        end
    end
    return nothing
end

# y = Dᵀ x (β = 0); threaded over column blocks.
function _f32_tmul!(y::AbstractVector{Float64}, D::Matrix{Float32}, x::AbstractVector{Float64})
    m = size(D, 2)
    Threads.@threads for r in _blocks(m, size(D, 1) < 2048 ? 1 : Threads.nthreads())
        _f32_cols!(y, D, x, r)
    end
    return y
end
function _f32_cols!(y, D, x, r)
    @inbounds for j in r
        s = 0.0
        @simd for i in axes(D, 1)
            s = muladd(Float64(D[i, j]), x[i], s)
        end
        y[j] = s
    end
    return nothing
end

# Fused block products for a few columns (k ≤ `_FUSED_MAX`): one pass over D per chunk
# of `_CHUNK` columns of X, four columns of D at a time against a row (or, transposed,
# column) block of the output that stays in cache. Memory-bound like the vector kernels
# up to k ≈ 8 (3–7× faster than the panels at k = 2–16); past ~50 columns the BLAS
# panels win.
const _CHUNK = 8
const _FUSED_MAX = 48
const _RB = 512           # rows per block of D X
const _IB = 1024          # rows of D per pass of Dᵀ X

# Y = D X (β = 0)
function _f32_bmul!(Y::AbstractMatrix{Float64}, D::Matrix{Float32}, X::AbstractMatrix{Float64})
    n, k = size(D, 1), size(X, 2)
    rbs = [r0:min(r0 + _RB - 1, n) for r0 in 1:_RB:n]
    for c0 in 1:_CHUNK:k
        cc = c0:min(c0 + _CHUNK - 1, k)
        Threads.@threads for r in rbs
            _f32_brows!(Y, D, X, r, cc)
        end
    end
    return Y
end
function _f32_brows!(Y, D, X, r, cc)
    m = size(D, 2)
    @inbounds for c in cc, i in r
        Y[i, c] = 0.0
    end
    j = 1
    @inbounds while j + 3 <= m
        for c in cc
            x1, x2, x3, x4 = X[j, c], X[j+1, c], X[j+2, c], X[j+3, c]
            @simd for i in r
                y = muladd(Float64(D[i, j]), x1, Y[i, c])
                y = muladd(Float64(D[i, j+1]), x2, y)
                y = muladd(Float64(D[i, j+2]), x3, y)
                Y[i, c] = muladd(Float64(D[i, j+3]), x4, y)
            end
        end
        j += 4
    end
    @inbounds while j <= m
        for c in cc
            xj = X[j, c]
            @simd for i in r
                Y[i, c] = muladd(Float64(D[i, j]), xj, Y[i, c])
            end
        end
        j += 1
    end
    return nothing
end

# Y = Dᵀ X (β = 0)
function _f32_btmul!(Y::AbstractMatrix{Float64}, D::Matrix{Float32}, X::AbstractMatrix{Float64})
    m, k = size(D, 2), size(X, 2)
    jbs = _blocks(m, 4 * Threads.nthreads())
    for c0 in 1:_CHUNK:k
        cc = c0:min(c0 + _CHUNK - 1, k)
        Threads.@threads for jr in jbs
            _f32_bcols!(Y, D, X, jr, cc)
        end
    end
    return Y
end
function _f32_bcols!(Y, D, X, jr, cc)
    n = size(D, 1)
    @inbounds for c in cc, j in jr
        Y[j, c] = 0.0
    end
    for i0 in 1:_IB:n
        ib = i0:min(i0 + _IB - 1, n)
        j = first(jr)
        @inbounds while j + 3 <= last(jr)
            for c in cc
                s1 = s2 = s3 = s4 = 0.0
                @simd for i in ib
                    x = X[i, c]
                    s1 = muladd(Float64(D[i, j]), x, s1)
                    s2 = muladd(Float64(D[i, j+1]), x, s2)
                    s3 = muladd(Float64(D[i, j+2]), x, s3)
                    s4 = muladd(Float64(D[i, j+3]), x, s4)
                end
                Y[j, c] += s1
                Y[j+1, c] += s2
                Y[j+2, c] += s3
                Y[j+3, c] += s4
            end
            j += 4
        end
        @inbounds while j <= last(jr)
            for c in cc
                s = 0.0
                @simd for i in ib
                    s = muladd(Float64(D[i, j]), X[i, c], s)
                end
                Y[j, c] += s
            end
            j += 1
        end
    end
    return nothing
end

# Columns `cols` of D widened into the first columns of `P`.
function _widen!(P::Matrix{Float64}, D::Matrix{Float32}, cols::UnitRange{Int})
    Threads.@threads for c in eachindex(cols)
        j = cols[c]
        @inbounds @simd for i in axes(D, 1)
            P[i, c] = D[i, j]
        end
    end
    return P
end

# Panel width: about 64 MB of Float64 per panel.
_panel(n, m) = clamp(cld(1 << 23, max(n, 1)), 16, max(m, 1))

# Y = op(D) X (β = 0) through Float64 panels of D and BLAS: `tr = false` D X,
# `tr = true` Dᵀ X, `left = true` X D (X on the left).
function _f32_panel_mul!(Y::AbstractMatrix{Float64}, D::Matrix{Float32}, X::AbstractMatrix{Float64};
                         tr::Bool = false, left::Bool = false)
    n, m = size(D)
    b = _panel(n, m)
    P = Matrix{Float64}(undef, n, min(b, m))
    !tr && !left && fill!(Y, 0.0)
    for j0 in 1:b:m
        cols = j0:min(j0 + b - 1, m)
        _widen!(P, D, cols)
        Pc = view(P, :, 1:length(cols))
        if left
            mul!(view(Y, :, cols), X, Pc)                    # (X D)[:, cols] = X D[:, cols]
        elseif tr
            mul!(view(Y, cols, :), transpose(Pc), X)         # (Dᵀ X)[cols, :]
        else
            mul!(Y, Pc, view(X, cols, :), 1.0, 1.0)          # D X = Σ D[:, cols] X[cols, :]
        end
    end
    return Y
end

const _F32T = Union{Transpose{Float64,Float32Matrix},Adjoint{Float64,Float32Matrix}}

# y = α op(A) x + β y in terms of the β = 0 kernels
function _axpby!(y, z, α, β)
    if iszero(β)
        isone(α) ? copyto!(y, z) : (y .= α .* z)
    else
        y .= α .* z .+ β .* y
    end
    return y
end

function LinearAlgebra.mul!(y::StridedVector{Float64}, A::Float32Matrix, x::StridedVector{Float64},
                            α::Number, β::Number)
    size(A, 2) == length(x) && size(A, 1) == length(y) || throw(DimensionMismatch("matrix–vector sizes"))
    z = iszero(β) && isone(α) ? y : similar(y)
    _f32_mul!(z, A.data, Base.mightalias(z, x) ? copy(x) : x)
    z === y || _axpby!(y, z, α, β)
    return y
end
function LinearAlgebra.mul!(y::StridedVector{Float64}, At::_F32T, x::StridedVector{Float64},
                            α::Number, β::Number)
    A = parent(At)
    size(A, 1) == length(x) && size(A, 2) == length(y) || throw(DimensionMismatch("matrix–vector sizes"))
    z = iszero(β) && isone(α) ? y : similar(y)
    _f32_tmul!(z, A.data, Base.mightalias(z, x) ? copy(x) : x)
    z === y || _axpby!(y, z, α, β)
    return y
end
function LinearAlgebra.mul!(Y::StridedMatrix{Float64}, A::Float32Matrix, X::StridedMatrix{Float64},
                            α::Number, β::Number)
    size(A, 2) == size(X, 1) && size(Y) == (size(A, 1), size(X, 2)) ||
        throw(DimensionMismatch("matrix sizes"))
    Z = iszero(β) && isone(α) && !Base.mightalias(Y, X) ? Y : similar(Y)
    if size(X, 2) == 1
        _f32_mul!(view(Z, :, 1), A.data, view(X, :, 1))
    elseif size(X, 2) <= _FUSED_MAX
        _f32_bmul!(Z, A.data, X)
    else
        _f32_panel_mul!(Z, A.data, X)
    end
    Z === Y || _axpby!(Y, Z, α, β)
    return Y
end
function LinearAlgebra.mul!(Y::StridedMatrix{Float64}, At::_F32T, X::StridedMatrix{Float64},
                            α::Number, β::Number)
    A = parent(At)
    size(A, 1) == size(X, 1) && size(Y) == (size(A, 2), size(X, 2)) ||
        throw(DimensionMismatch("matrix sizes"))
    Z = iszero(β) && isone(α) && !Base.mightalias(Y, X) ? Y : similar(Y)
    if size(X, 2) == 1
        _f32_tmul!(view(Z, :, 1), A.data, view(X, :, 1))
    elseif size(X, 2) <= _FUSED_MAX
        _f32_btmul!(Z, A.data, X)
    else
        _f32_panel_mul!(Z, A.data, X; tr = true)
    end
    Z === Y || _axpby!(Y, Z, α, β)
    return Y
end
function LinearAlgebra.mul!(Y::StridedMatrix{Float64}, X::StridedMatrix{Float64}, A::Float32Matrix,
                            α::Number, β::Number)
    size(X, 2) == size(A, 1) && size(Y) == (size(X, 1), size(A, 2)) ||
        throw(DimensionMismatch("matrix sizes"))
    Z = iszero(β) && isone(α) && !Base.mightalias(Y, X) ? Y : similar(Y)
    _f32_panel_mul!(Z, A.data, X; left = true)
    Z === Y || _axpby!(Y, Z, α, β)
    return Y
end

# Products with strided vectors of first-order Float64 duals (the radiosities of a dual
# solve): the values and partials as the columns of one matrix product (generic dual
# products would index element by element). Strided rather than abstract vectors, so
# FillArrays zero vectors keep their own methods without an ambiguity.
for T in (:Float32Matrix, :_F32T)
    @eval function Base.:*(A::$T, x::StridedVector{ForwardDiff.Dual{Tag,Float64,N}}) where {Tag,N}
        X = Matrix{Float64}(undef, length(x), N + 1)
        for (i, d) in enumerate(x)
            X[i, 1] = ForwardDiff.value(d)
            p = ForwardDiff.partials(d)
            for k in 1:N
                X[i, k+1] = p[k]
            end
        end
        Y = A * X
        return [ForwardDiff.Dual{Tag}(Y[i, 1], ntuple(k -> Y[i, k+1], Val(N))...) for i in axes(Y, 1)]
    end
end

# F .= s .* F .* s' (symmetric balancing), threaded over columns.
function _balance!(F::Float32Matrix, s::AbstractVector{Float64})
    D = F.data
    Threads.@threads for j in axes(D, 2)
        sj = s[j]
        @inbounds @simd for i in axes(D, 1)
            D[i, j] = s[i] * D[i, j] * sj
        end
    end
    return F
end
_balance!(F::AbstractMatrix, s) = (F .= s .* F .* s'; F)
