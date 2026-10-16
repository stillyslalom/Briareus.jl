# Mode decompositions of fields on (approximately) spherical element sets.

"""
    solid_angle(target, e, center) -> Float64

Exact solid angle subtended by element `e` at `center` (Van Oosterom–Strackee formula
summed over its triangles). Mode decompositions use it as the quadrature weight, so the
weights of a closed surface around `center` sum to exactly 4π.
"""
function solid_angle(t::Target, e::Integer, center = zero(Vec3))
    Ω = 0.0
    for k in t.elem_tris[e]
        a = t.tri_v0[k] - center
        b = a + t.tri_e1[k]
        c = a + t.tri_e2[k]
        la, lb, lc = norm(a), norm(b), norm(c)
        num = abs(dot(a, cross(b, c)))
        den = la * lb * lc + dot(a, b) * lc + dot(a, c) * lb + dot(b, c) * la
        Ω += 2 * atan(num, den)
    end
    return Ω
end

"""
    legendre(l, μ)

Legendre polynomial `Pₗ(μ)`.
"""
function legendre(l::Integer, μ)
    l == 0 && return one(μ)
    p0, p1 = one(μ), μ
    for n in 1:l-1
        p0, p1 = p1, ((2n + 1) * μ * p1 - n * p0) / (n + 1)
    end
    return p1
end

# Sub-triangles of element e at the given refinement, as (solid angle, unit direction)
# about `c`; used to integrate angular functions over each facet's projection.
function _facet_samples(t::Target, e::Integer, c::Vec3, level::Integer)
    out = Tuple{Float64,Vec3}[]
    for k in t.elem_tris[e]
        a = t.tri_v0[k] - c
        _subtri_solid!(out, a, a + t.tri_e1[k], a + t.tri_e2[k], level)
    end
    return out
end

function _subtri_solid!(out, a, b, c, level)
    if level == 0
        la, lb, lc = norm(a), norm(b), norm(c)
        num = abs(dot(a, cross(b, c)))
        den = la * lb * lc + dot(a, b) * lc + dot(a, c) * lb + dot(b, c) * la
        push!(out, (2 * atan(num, den), normalize(a / la + b / lb + c / lc)))
        return
    end
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    _subtri_solid!(out, a, ab, ca, level - 1)
    _subtri_solid!(out, ab, b, bc, level - 1)
    _subtri_solid!(out, ca, bc, c, level - 1)
    _subtri_solid!(out, ab, bc, ca, level - 1)
end

"""
    legendre_weights(target, elems; lmax = 8, axis = ẑ, center = 0, level = 3) -> Matrix

Matrix `W` with `legendre_moments(values, …) == W * values[elems]`:
`W[l+1, k] = (2l+1)/Ω ∫ Pₗ(cosθ) dΩ` over element `elems[k]`, integrated over `4^level`
sub-triangles per facet with exact solid angles, and `Ω` the total solid angle of
`elems` (4π for a closed surface around `center`). Treating the field as constant on
each element makes a uniform field have no `l > 0` moments, however the facets tile
the sphere.
"""
function legendre_weights(t::Target, elems::AbstractVector{<:Integer}; lmax::Integer = 8,
                          axis = Vec3(0, 0, 1), center = zero(Vec3), level::Integer = 3)
    ax = normalize(Vec3(axis))
    c = Vec3(center)
    W = zeros(lmax + 1, length(elems))
    Ω = 0.0
    for (k, e) in enumerate(elems)
        for (dΩ, u) in _facet_samples(t, e, c, level)
            μ = dot(u, ax)
            Ω += dΩ
            p0, p1 = 1.0, μ
            W[1, k] += dΩ
            lmax >= 1 && (W[2, k] += μ * dΩ)
            for l in 1:lmax-1
                p0, p1 = p1, ((2l + 1) * μ * p1 - l * p0) / (l + 1)
                W[l+2, k] += p1 * dΩ
            end
        end
    end
    for l in 0:lmax
        W[l+1, :] .*= (2l + 1) / Ω
    end
    return W
end

"""
    legendre_moments(values, target, elems; lmax = 8, axis = ẑ, center = 0, level = 3) -> Vector

Legendre coefficients `aₗ` (index `l+1`) of the azimuthally averaged field,
`aₗ = (2l+1)/(4π) ∫ f Pₗ(cosθ) dΩ`, where `θ` is measured from `axis` about `center`
and `values[e]` is indexed by element (taken constant over each element). `aₗ/a₀` is
the usual "Pₗ/P₀" asymmetry. `legendre_weights` precomputes the linear operator.
"""
legendre_moments(values::AbstractVector, t::Target, elems::AbstractVector{<:Integer}; kwargs...) =
    legendre_weights(t, elems; kwargs...) * values[elems]

"""
    ylm_moments(values, target, elems; lmax = 8, axis = ẑ, center = 0, level = 2) -> Matrix

Coefficients of real orthonormal spherical harmonics, `c[l+1, l+1+m]` for
`-l ≤ m ≤ l`, with the polar axis along `axis` (field constant on each element,
harmonics integrated over each facet). Use `mode_power` for rotation-invariant per-`l`
amplitudes.
"""
function ylm_moments(values::AbstractVector, t::Target, elems::AbstractVector{<:Integer};
                     lmax::Integer = 8, axis = Vec3(0, 0, 1), center = zero(Vec3),
                     level::Integer = 2)
    ax = normalize(Vec3(axis))
    ex = normalize(abs(ax[1]) < 0.9 ? cross(Vec3(1, 0, 0), ax) : cross(Vec3(0, 1, 0), ax))
    ey = cross(ax, ex)
    c = Vec3(center)
    T = eltype(values)
    C = zeros(T, lmax + 1, 2lmax + 1)
    Ω = 0.0
    Y = zeros(lmax + 1, 2lmax + 1)
    Ye = zeros(lmax + 1, 2lmax + 1)
    for e in elems
        fill!(Ye, 0)
        for (dΩ, u) in _facet_samples(t, e, c, level)
            Ω += dΩ
            _real_ylm!(Y, lmax, dot(u, ax), atan(dot(u, ey), dot(u, ex)))
            Ye .+= Y .* dΩ
        end
        C .+= values[e] .* Ye
    end
    return C .* (4π / Ω)       # normalize the covered solid angle to 4π
end

# Real orthonormal spherical harmonics via normalized associated Legendre recursion.
function _real_ylm!(Y, lmax, μ, φ)
    s = sqrt(max(0.0, 1 - μ^2))
    # P̄ₗᵐ fully normalized such that Yₗᵐ = P̄ₗᵐ(μ) * (cos mφ or sin mφ) [* √2 for m≠0]
    P = zeros(lmax + 1, lmax + 1)
    P[1, 1] = sqrt(1 / (4π))
    for m in 1:lmax
        P[m+1, m+1] = -sqrt((2m + 1) / (2m)) * s * P[m, m]
    end
    for m in 0:lmax-1
        P[m+2, m+1] = sqrt(2m + 3) * μ * P[m+1, m+1]
    end
    for m in 0:lmax, l in m+2:lmax
        a = sqrt((4l^2 - 1) / (l^2 - m^2))
        b = sqrt(((l - 1)^2 - m^2) / (4(l - 1)^2 - 1))
        P[l+1, m+1] = a * (μ * P[l, m+1] - b * P[l-1, m+1])
    end
    for l in 0:lmax
        Y[l+1, lmax+1] = P[l+1, 1]
        for m in 1:l
            Y[l+1, lmax+1+m] = sqrt(2) * P[l+1, m+1] * cos(m * φ)
            Y[l+1, lmax+1-m] = sqrt(2) * P[l+1, m+1] * sin(m * φ)
        end
    end
    return Y
end

"""
    mode_power(C) -> Vector

Per-`l` amplitude `√(Σₘ cₗₘ²) / c₀₀` from `ylm_moments`, a rotation-invariant
nonuniformity measure (index `l+1`). Equals `|aₗ/a₀|/√(2l+1)` for an axisymmetric field.
"""
function mode_power(C::AbstractMatrix)
    lmax = size(C, 1) - 1
    return [sqrt(sum(abs2, C[l+1, lmax+1-l:lmax+1+l])) / abs(C[1, lmax+1]) for l in 0:lmax]
end

"""
    rms_nonuniformity(values, weights) -> Real

Weighted RMS deviation over weighted mean.
"""
function rms_nonuniformity(v::AbstractVector, w::AbstractVector)
    W = sum(w)
    μ = sum(v .* w) / W
    return sqrt(sum(w .* (v .- μ) .^ 2) / W) / μ
end
