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
    mode_center(target, elems) -> Vec3

Center about which the mode analysis of `elems` (`legendre_weights`, `ylm_weights` and
the scenario caches built from them) measures directions by default: the center of
their `SphereSurface` when every element lies on sphere components with one center (a
capsule, also when it is offset from chamber center), otherwise the origin. Analyzing
an offset capsule about the origin mixes its P2 into a spurious mode 1 at first order
in the offset (0.74 % of mode 1 per 100 µm on a 0.83 mm capsule with P2/P0 ≈ −0.1).
"""
function mode_center(t::Target, elems::AbstractVector{<:Integer})
    c = nothing
    for e in elems
        s = t.components[t.elem_comp[e]].surface
        s isa SphereSurface || return zero(Vec3)
        c === nothing ? (c = s.center) : (s.center == c || return zero(Vec3))
    end
    return c === nothing ? zero(Vec3) : c
end

_mode_center(t, elems, center) = center === nothing ? mode_center(t, elems) : Vec3(center...)

"""
    legendre_weights(target, elems; lmax = 8, axis = ẑ, center = mode_center(target, elems), level = 3) -> Matrix

Matrix `W` with `legendre_moments(values, …) == W * values[elems]`:
`W[l+1, k] = (2l+1)/Ω ∫ Pₗ(cosθ) dΩ` over element `elems[k]`, integrated over `4^level`
sub-triangles per facet with exact solid angles, and `Ω` the total solid angle of
`elems` (4π for a closed surface around `center`). Treating the field as constant on
each element makes a uniform field have no `l > 0` moments, however the facets tile
the sphere.
"""
function legendre_weights(t::Target, elems::AbstractVector{<:Integer}; lmax::Integer = 8,
                          axis = Vec3(0, 0, 1), center = nothing, level::Integer = 3)
    ax = normalize(Vec3(axis))
    c = _mode_center(t, elems, center)
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
    legendre_moments(values, target, elems; lmax = 8, axis = ẑ, center = mode_center(target, elems), level = 3) -> Vector

Legendre coefficients `aₗ` (index `l+1`) of the azimuthally averaged field,
`aₗ = (2l+1)/(4π) ∫ f Pₗ(cosθ) dΩ`, where `θ` is measured from `axis` about `center`
and `values[e]` is indexed by element (taken constant over each element). `aₗ/a₀` is
the usual "Pₗ/P₀" asymmetry. `legendre_weights` precomputes the linear operator.
"""
legendre_moments(values::AbstractVector, t::Target, elems::AbstractVector{<:Integer}; kwargs...) =
    legendre_weights(t, elems; kwargs...) * values[elems]

"""
    mode_index(l, m) -> Int

Position of the real spherical harmonic `(l, m)` in the flat ordering used by
`ylm_weights`: `l² + l + m + 1`, i.e. `(0,0), (1,-1), (1,0), (1,1), (2,-2), …`.
"""
mode_index(l::Integer, m::Integer) = l^2 + l + m + 1

"""
    ylm_frame(axis) -> (ẑ, x̂, ŷ)

Right-handed frame for spherical harmonics with the polar axis along `axis`. For the
default `axis = ẑ` this is the chamber frame, so azimuthal mode phases are relative to
`+x` (NIF: the direction of Switchyard 1).
"""
function ylm_frame(axis)
    ax = normalize(Vec3(axis))
    if abs(ax[3]) > 1 - 1e-12
        ex = Vec3(1, 0, 0)
        ey = cross(ax, ex)
        return ax, ex, ey
    end
    ex = normalize(abs(ax[1]) < 0.9 ? cross(Vec3(1, 0, 0), ax) : cross(Vec3(0, 1, 0), ax))
    return ax, ex, cross(ax, ex)
end

"""
    ylm_weights(target, elems; lmax = 8, axis = ẑ, center = mode_center(target, elems), level = 2) -> Matrix

Matrix `W` with `W * values[elems]` the coefficients of real orthonormal spherical
harmonics in the flat `mode_index` ordering: `Y_l0 = P̄_l0`, `Y_lm ∝ cos mφ` for `m > 0`
and `∝ sin |m|φ` for `m < 0`, without the Condon–Shortley phase, so `Y_11 ∝ +x` and
`Y_1,-1 ∝ +y`. Each element's field is taken as constant over the
element and the harmonics are integrated over its facets with exact sub-triangle solid
angles (`4^level` per triangle); the covered solid angle is normalized to 4π.
"""
function ylm_weights(t::Target, elems::AbstractVector{<:Integer}; lmax::Integer = 8,
                     axis = Vec3(0, 0, 1), center = nothing, level::Integer = 2)
    ax, ex, ey = ylm_frame(axis)
    c = _mode_center(t, elems, center)
    W = zeros((lmax + 1)^2, length(elems))
    Y = zeros(lmax + 1, 2lmax + 1)
    P = zeros(lmax + 1, lmax + 1)
    Ω = 0.0
    for (k, e) in enumerate(elems)
        for (dΩ, u) in _facet_samples(t, e, c, level)
            Ω += dΩ
            _real_ylm!(Y, lmax, dot(u, ax), atan(dot(u, ey), dot(u, ex)), P)
            for l in 0:lmax, m in -l:l
                W[mode_index(l, m), k] += Y[l+1, lmax+1+m] * dΩ
            end
        end
    end
    return W .* (4π / Ω)
end

"""
    ylm_matrix(c, lmax) -> Matrix

Flat coefficients (`mode_index` order) as the `(lmax+1) × (2lmax+1)` matrix
`C[l+1, lmax+1+m]` used by `ylm_moments` and `mode_power` (unused entries zero).
"""
function ylm_matrix(c::AbstractVector, lmax::Integer)
    C = zeros(eltype(c), lmax + 1, 2lmax + 1)
    for l in 0:lmax, m in -l:l
        C[l+1, lmax+1+m] = c[mode_index(l, m)]
    end
    return C
end

"""
    ylm_moments(values, target, elems; lmax = 8, axis = ẑ, center = mode_center(target, elems), level = 2) -> Matrix

Coefficients of real orthonormal spherical harmonics, `c[l+1, l+1+m]` for
`-l ≤ m ≤ l` (`m > 0`: cos mφ, `m < 0`: sin |m|φ), with the polar axis along `axis`
(field constant on each element, harmonics integrated over each facet; see
`ylm_weights`). Use `mode_power` for rotation-invariant per-`l` amplitudes.
"""
ylm_moments(values::AbstractVector, t::Target, elems::AbstractVector{<:Integer};
            lmax::Integer = 8, kwargs...) =
    ylm_matrix(ylm_weights(t, elems; lmax, kwargs...) * values[elems], lmax)

"""
    relative_mode_power(c, l) -> Real

`Σₘ (c_lm / c_00)²` from flat coefficients `c` (`mode_index` order): the squared
rotation-invariant amplitude of degree `l` relative to the mean. Smooth (no square
root), so suited to objectives; `mode_power` is its square root. For an axisymmetric
field it equals `(aₗ/a₀)² / (2l+1)`.
"""
relative_mode_power(c::AbstractVector, l::Integer) =
    sum(abs2, (c[mode_index(l, m)] for m in -l:l)) / c[1]^2

"""
    mode_amplitude(c, l) -> Real

Rotation-invariant amplitude of degree `l` relative to the mean, in the normalization
of Legendre moments: `√((2l+1) Σₘ (c_lm/c_00)²)`. It equals `|aₗ/a₀|` for a field
axisymmetric about any axis, so for `l = 1` it is the fractional amplitude `ε` of a
dipole drive `1 + ε cos γ`, the quantity hot-spot velocity sensitivities refer to.
"""
mode_amplitude(c::AbstractVector, l::Integer) = sqrt((2l + 1) * relative_mode_power(c, l))

# Real orthonormal spherical harmonics via normalized associated Legendre recursion; `P`
# is an (lmax+1)² work array (every entry read is written first).
function _real_ylm!(Y, lmax, μ, φ, P = zeros(lmax + 1, lmax + 1))
    s = sqrt(max(0.0, 1 - μ^2))
    # P̄ₗᵐ fully normalized such that Yₗᵐ = P̄ₗᵐ(μ) * (cos mφ or sin mφ) [* √2 for m≠0]
    P[1, 1] = sqrt(1 / (4π))
    for m in 1:lmax
        P[m+1, m+1] = sqrt((2m + 1) / (2m)) * s * P[m, m]     # no Condon–Shortley phase
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
    end
    for m in 1:lmax
        cm, sm = cos(m * φ), sin(m * φ)
        for l in m:lmax
            Y[l+1, lmax+1+m] = sqrt(2) * P[l+1, m+1] * cm
            Y[l+1, lmax+1-m] = sqrt(2) * P[l+1, m+1] * sm
        end
    end
    return Y
end

"""
    mode_power(C; m = :all) -> Vector

Per-`l` amplitude `σₗ = √(Σₘ cₗₘ²) / c₀₀` from `ylm_moments` (index `l+1`; `ylm_matrix`
turns flat coefficients into `C`), a rotation-invariant nonuniformity spectrum: the rms
of the field over the sphere relative to its mean is `√(Σ_{l≥1} σₗ²)` (`mode_rms`).
Equals `|aₗ/a₀|/√(2l+1)` for an axisymmetric field. `m = :zero` keeps only the zonal
harmonic `m = 0` of each degree (`|c_l0|/c₀₀`, the spectrum of the azimuthally averaged
field, `|aₗ/a₀|/√(2l+1)` from its `legendre_moments`) and `m = :nonzero` the rest; their
squares add to the `:all` spectrum's. The zonal part depends on the polar `axis` the
coefficients were taken about; `:all` does not.
"""
function mode_power(C::AbstractMatrix; m::Symbol = :all)
    lmax = size(C, 1) - 1
    m === :all &&
        return [sqrt(sum(abs2, C[l+1, lmax+1-l:lmax+1+l])) / abs(C[1, lmax+1]) for l in 0:lmax]
    m === :zero && return [abs(C[l+1, lmax+1]) / abs(C[1, lmax+1]) for l in 0:lmax]
    m === :nonzero ||
        throw(ArgumentError("m must be :all, :zero or :nonzero, got :$m"))
    return [sqrt(sum(abs2, C[l+1, lmax+1-l:lmax]; init = zero(eltype(C))) +
                 sum(abs2, C[l+1, lmax+2:lmax+1+l]; init = zero(eltype(C)))) / abs(C[1, lmax+1])
            for l in 0:lmax]
end

"""
    mode_rms(c, lmax; m = :all) -> Real

RMS nonuniformity over degrees `1 ≤ l ≤ lmax`, `σ = √(Σₗ Σₘ c_lm²) / |c₀₀|`, from flat
real-harmonic coefficients `c` (`mode_index` order: `ylm_weights(…) * values`, or
`flux_ylm`): the rms deviation of the field over the sphere relative to its mean, with
the degrees above `lmax` filtered out (`σ² = Σₗ σₗ²` with `σₗ` from `mode_power`). With
`m = :zero` only the zonal harmonics `m = 0` are summed (the rms of the azimuthally
averaged field about the coefficients' polar axis), with `m = :nonzero` only the rest;
`mode_rms(c, l; m = :zero)^2 + mode_rms(c, l; m = :nonzero)^2 == mode_rms(c, l)^2` up to
rounding. The split is the m = 0 / m ≠ 0 decomposition of polar-direct-drive designs
(Yeamans et al., Nucl. Fusion 61, 046031, 2021); the filtered `σ` over `l ≤ 30` is the
illumination nonuniformity of Barlow et al. (PRL 133, 175101, 2024).

Generic in the number type (it sits on the optimization path). It has a square root,
so its derivative is singular for an exactly uniform field; objectives can sum
`relative_mode_power` over `l` instead, which is `σ²` without the root.

For hard-sphere laser illumination (no plasma: every ray straight, the incident
`I cos ψ` absorbed where it lands) the field is the absorbed intensity of a `deposit` on
a bare sphere, `t = Target(Component(:capsule, icosphere_mesh(R; n), material; surface =
SphereSurface(R)))`, with `dep.element ./ t.elem_area` projected by
`W = ylm_weights(t, 1:nelements(t); lmax = 30)`, built once (`DepositionOperator(t,
beams; ncone = 1)` for collimated beams). The element values average the illumination
over each facet, which lowers `σₗ` at second order in the element size: on OMEGA's 60
beams with SG4 spots on R = 0.5 mm, `σ` (l ≤ 30) is 3.5 % low at `n = 20` (edges
0.055 R) and 1.7 % low at `n = 30`; a Richardson extrapolation of the two is 0.3 % low
(`reference/research-notes.md`, "Hard-sphere direct-drive illumination").
"""
function mode_rms(c::AbstractVector, lmax::Integer; m::Symbol = :all)
    m in (:all, :zero, :nonzero) ||
        throw(ArgumentError("m must be :all, :zero or :nonzero, got :$m"))
    length(c) >= (lmax + 1)^2 ||
        throw(ArgumentError("coefficients reach l = $(isqrt(length(c)) - 1), not lmax = $lmax"))
    s = zero(c[1])^2
    for l in 1:lmax, k in -l:l
        (m === :all || (m === :zero) == (k == 0)) && (s += c[mode_index(l, k)]^2)
    end
    return sqrt(s) / abs(c[1])
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
