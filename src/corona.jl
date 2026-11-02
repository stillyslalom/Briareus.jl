# Laser absorption in a static, spherically symmetric corona (polar direct drive, Stage 1;
# reference/litreview/polar_direct_drive.md, decision 3 (b); research notes "Corona
# absorption operator"). Rays are refracted and absorbed by inverse bremsstrahlung; the
# ray paths depend only on the profile and the impact parameter, so they are integrated
# once per profile, by radial quadrature, as harmonic kernels. A pointing change then
# only reweights a fixed lattice of rays per beam and f-cone direction.

# Kruer, The Physics of Laser Plasma Interactions, ch. 5: the electron–ion collision
# frequency ν_ei = 4√(2π) Z e⁴ n_e lnΛ / (3 m_e^{1/2} T_e^{3/2}) (Gaussian units), written
# ν_ei = _NU_EI Z n_e lnΛ / T_e^{3/2} with n_e in cm⁻³ and T_e in eV (s⁻¹; 2.91e-6, the
# NRL Plasma Formulary's electron collision rate), and the critical density
# n_c = π m_e c² / (e² λ²) (cm⁻³). CODATA 2018 constants, cgs.
const _E_ESU = 4.80320471e-10                 # elementary charge (statC)
const _ME_G = 9.1093837015e-28                # electron mass (g)
const _EV_ERG = 1.602176634e-12               # 1 eV (erg)
const _C_CM = 2.99792458e10                   # speed of light (cm/s)
const _NU_EI = 4sqrt(2π) / 3 * _E_ESU^4 / (sqrt(_ME_G) * _EV_ERG^1.5)
_critical_density(λ_µm) = π * _ME_G * _C_CM^2 / (_E_ESU^2 * (1e-4λ_µm)^2)

"""
    CoronaProfile(r, ne_over_nc, Te; Z, lnΛ, wavelength = 0.351)
    CoronaProfile(; critical_radius, scale_length, absorption,
                  outer_radius = critical_radius + 14scale_length)

A static, spherically symmetric laser corona for [`CoronaOperator`](@ref), given by the
electron density as a fraction ``n_e/n_c`` of the light's critical density (not in cm⁻³):
that ratio and the inverse-bremsstrahlung absorption rate as functions of the radius from
the target center.

The first form tabulates a snapshot, such as a 1-D hydrodynamics output: radii `r` (mm,
increasing, not necessarily uniform), `ne_over_nc` (``n_e/n_c``, positive), electron
temperature `Te` (eV), and the mean ionization `Z` and Coulomb logarithm `lnΛ` (each a
number or a vector over `r`; no defaults). `wavelength` is the laser's vacuum wavelength
(µm), which sets ``n_c``. The absorption rate per unit length of a ray, divided by its
refractive index, is `(ν_ei/c)(n_e/n_c)` with ``ν_{ei}`` the electron–ion collision
frequency of Kruer (*The Physics of Laser Plasma Interactions*, ch. 5),
`2.91e-6 Z n_e lnΛ / T_e^{3/2}` s⁻¹ (``n_e`` in cm⁻³, ``T_e`` in eV). Both `ln(n_e/n_c)` and the logarithm of
that rate are interpolated by cubic Hermite polynomials (three-point slopes), so an
exponential profile is represented exactly; the profile is zero beyond `r[end]`.
`ne_over_nc[1] ≥ 1` is required: every ray must turn before the innermost radius (the
light never reaches a hard surface, whose absorption the corona model does not include).

The second form is the analytic isothermal corona `n_e/n_c = exp(-(r - r_c)/L)` with
`critical_radius` ``r_c`` and `scale_length` L (mm), truncated at `outer_radius`, whose
absorption rate is set by `absorption`, the fraction η⊥ that a planar profile of the
same scale length absorbs at normal incidence: η⊥ = 1 - exp(-8``ν_c`` L/3c) (Kruer; the
η⊥ of `deposit(op, beams; absorption)`, Murakami & Nishi 2017, Eq. 29). It is also the exact
absorption of a radial ray when the truncation is far out (e⁻¹⁴ of ``n_c`` by default).

Not calibrated: a profile is the user's input (a snapshot or an engineering estimate),
and `Z`, `lnΛ` and `absorption` have no defaults.
"""
struct CoronaProfile
    r::Vector{Float64}
    lne::Vector{Float64}       # ln(n_e/n_c)
    dlne::Vector{Float64}      # its slope at the nodes (mm⁻¹)
    lrate::Vector{Float64}     # ln((ν_ei/c)(n_e/n_c)), the rate in mm⁻¹
    dlrate::Vector{Float64}
end

function CoronaProfile(r::AbstractVector{<:Real}, ne_over_nc::AbstractVector{<:Real},
                       Te::AbstractVector{<:Real}; Z, lnΛ, wavelength::Real = 0.351)
    n = length(r)
    n >= 2 && length(ne_over_nc) == n && length(Te) == n ||
        throw(DimensionMismatch("r, ne_over_nc and Te need the same length (at least 2)"))
    for (name, v) in (("Z", Z), ("lnΛ", lnΛ))
        v isa Real || length(v) == n ||
            throw(DimensionMismatch("$name must be a number or a vector over r"))
    end
    all(>(0), ne_over_nc) && all(>(0), Te) && all(>(0), Z) && all(>(0), lnΛ) ||
        throw(ArgumentError("ne_over_nc, Te, Z and lnΛ must be positive"))
    wavelength > 0 || throw(ArgumentError("wavelength must be positive"))
    nc = _critical_density(wavelength)
    rate = @. _NU_EI * Z * lnΛ * nc * ne_over_nc^2 / Te^1.5 / (10 * _C_CM)     # mm⁻¹
    return _corona_profile(collect(Float64, r), log.(Float64.(ne_over_nc)), log.(Float64.(rate)))
end

function CoronaProfile(; critical_radius::Real, scale_length::Real, absorption::Real,
                       outer_radius::Real = critical_radius + 14scale_length)
    rc, L = Float64(critical_radius), Float64(scale_length)
    rc > L > 0 || throw(ArgumentError("need critical_radius > scale_length > 0"))
    0 < absorption < 1 || throw(ArgumentError("absorption must be in (0, 1), got $absorption"))
    outer_radius > rc || throw(ArgumentError("outer_radius must exceed critical_radius"))
    νc = -3log1p(-absorption) / (8L)                     # ν_c/c (mm⁻¹)
    r = [rc - L, Float64(outer_radius)]
    lne = @. -(r - rc) / L                               # linear: the table is exact
    return _corona_profile(r, lne, log(νc) .+ 2 .* lne)
end

function _corona_profile(r, lne, lrate)
    all(>(0), diff(r)) || throw(ArgumentError("radii must increase"))
    r[1] > 0 || throw(ArgumentError("radii must be positive"))
    lne[1] >= 0 ||
        throw(ArgumentError("the profile must reach the critical density at its first radius " *
                            "(ne_over_nc[1] = $(exp(lne[1])) < 1): light would reach a hard surface"))
    return CoronaProfile(r, lne, _hermite_slopes(r, lne), lrate, _hermite_slopes(r, lrate))
end

function Base.show(io::IO, p::CoronaProfile)
    @printf(io, "CoronaProfile(%d radii, %.4g–%.4g mm, n_e/n_c %.3g–%.3g)", length(p.r),
            p.r[1], p.r[end], exp(p.lne[1]), exp(p.lne[end]))
end

# Three-point slopes on a non-uniform grid (exact for quadratics), one-sided at the ends.
function _hermite_slopes(r, y)
    n = length(r)
    d = similar(y)
    d[1] = (y[2] - y[1]) / (r[2] - r[1])
    d[n] = (y[n] - y[n-1]) / (r[n] - r[n-1])
    for k in 2:n-1
        h0, h1 = r[k] - r[k-1], r[k+1] - r[k]
        d[k] = (h0^2 * (y[k+1] - y[k]) + h1^2 * (y[k] - y[k-1])) / (h0 * h1 * (h0 + h1))
    end
    return d
end

@inline function _hermite(r, y, dy, x)
    i = clamp(searchsortedlast(r, x), 1, length(r) - 1)
    h = r[i+1] - r[i]
    t = (x - r[i]) / h
    t2 = t * t
    t3 = t2 * t
    return (2t3 - 3t2 + 1) * y[i] + (t3 - 2t2 + t) * h * dy[i] + (3t2 - 2t3) * y[i+1] +
           (t3 - t2) * h * dy[i+1]
end

# n_e/n_c and the absorption rate (ν_ei/c)(n_e/n_c) (mm⁻¹) at radius x; zero outside.
_corona_ne(p::CoronaProfile, x) = x > p.r[end] ? 0.0 : exp(_hermite(p.r, p.lne, p.dlne, x))
_corona_rate(p::CoronaProfile, x) = x > p.r[end] ? 0.0 : exp(_hermite(p.r, p.lrate, p.dlrate, x))

# ---------------------------------------------------------------------------------------
# Ray kernels: one ray of impact parameter b by radial quadrature (no ray integration).
# Geometric optics in a cold plasma, refractive index n² = 1 - n_e/n_c: the ray stays in
# its plane, turns where r n(r) = b, sweeps the angle ∫ b dr / (r √(r²n² - b²)) and
# accumulates the optical depth ∫ κ r dr / √(r²n² - b²) per pass (κ the rate above).
# ---------------------------------------------------------------------------------------

const _GL4 = ((-0.8611363115940526, 0.3478548451374538), (-0.3399810435848563, 0.6521451548625461),
              (0.3399810435848563, 0.6521451548625461), (0.8611363115940526, 0.3478548451374538))

# The outermost turning radius in [rmin, R): r² (1 - n_e/n_c) = b² (n_e(rmin) ≥ n_c).
function _turning_radius(ne::F, R, rmin, b) where {F}
    h(r) = r^2 * (1 - ne(r)) - b^2
    hi = R * (1 - 1e-12)
    lo = hi
    step = (R - rmin) / 200
    while h(lo) > 0
        hi = lo
        lo = max(lo - step, rmin)
        lo == rmin && break
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        h(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    return hi
end

# Deposits of a ray with impact parameter b entering the sphere of radius R: angles α
# (from the incoming side, -d, toward the impact vector) and absorbed fractions f of 2N
# panels (N in and N out, on r = r_t + (R - r_t) t² with uniform t, which removes the
# 1/√ singularity at the turning point; 4-point Gauss each), each panel's deposit at its
# mid angle; the escaping fraction (Σ f + escaped = 1); and the exit angle. With `edges`
# (increasing radii, covering [rmin, R]), every edge inside (r_t, R) is also a panel
# boundary, so each panel lies in one shell, and the fifth result is each deposit's shell
# (`edges[s] ≤ r < edges[s+1]`; zeros without edges, whose panels are as above).
function _ray_deposits(ne::F, rate::G, R, rmin, b, N, edges = nothing) where {F,G}
    rt = _turning_radius(ne, R, rmin, b)
    Δ = R - rt
    tb = [k / N for k in 0:N]
    if edges !== nothing
        append!(tb, (sqrt((e - rt) / Δ) for e in edges if rt < e < R))
        sort!(unique!(tb))
    end
    N = length(tb) - 1
    A = zeros(N)
    D = zeros(N)
    for k in 1:N
        t0, t1 = tb[k], tb[k+1]
        a = 0.0
        τ = 0.0
        for (ξ, w) in _GL4
            t = (t0 + t1) / 2 + (t1 - t0) / 2 * ξ
            r = rt + Δ * t^2
            dr = Δ * t * (t1 - t0) * w
            q = r^2 * (1 - ne(r)) - b^2
            q <= 0 && continue
            s = sqrt(q)
            a += b * dr / (r * s)
            τ += rate(r) * r * dr / s
        end
        A[k], D[k] = a, τ
    end
    α = zeros(2N)
    f = zeros(2N)
    shell = zeros(Int, 2N)
    θ = asin(min(b / R, 1.0))
    τ = 0.0
    for (i, k) in enumerate(Iterators.flatten((N:-1:1, 1:N)))
        α[i] = θ + A[k] / 2
        f[i] = exp(-τ) * -expm1(-D[k])
        θ += A[k]
        τ += D[k]
        if edges !== nothing
            rmid = rt + Δ * ((tb[k] + tb[k+1]) / 2)^2
            shell[i] = clamp(searchsortedlast(edges, rmid), 1, length(edges) - 1)
        end
    end
    return α, f, exp(-τ), θ, shell
end

_ray_deposits(p::CoronaProfile, b, N, edges = nothing) =
    _ray_deposits(x -> _corona_ne(p, x), x -> _corona_rate(p, x), p.r[end], p.r[1], b, N, edges)

# The critical radius of a profile: the outermost radius where n_e/n_c = 1 (ln = 0).
function _critical_radius(p::CoronaProfile)
    i = something(findlast(>=(0), p.lne), 1)    # p.lne[1] ≥ 0 by construction
    i == length(p.r) && return p.r[end]
    lo, hi = p.r[i], p.r[i+1]
    for _ in 1:200
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        _hermite(p.r, p.lne, p.dlne, mid) >= 0 ? (lo = mid) : (hi = mid)
    end
    return lo
end

# Position of (l, m ≥ 0) in the zonal-kernel ordering (m-major, l = m…L).
@inline _zonal_index(l, m, L) = m * (L + 1) - (m * (m - 1)) ÷ 2 + (l - m) + 1

# Harmonic kernel of a ray in its own frame (pole -d, azimuth 0 along the impact vector):
# κ[(l, m)] = Σ f Ȳ_lm(α), so that a ray at impact azimuth ψ contributes κ cos mψ to the
# frame coefficient (l, m) and κ sin mψ to (l, -m). A deposit past the antipode (α > π)
# lies at polar angle 2π - α and azimuth π.
function _ray_kernel!(κ, Y, P, L, α, f)
    for (a, w) in zip(α, f)
        w > 0 || continue
        a = mod(a, 2π)
        _real_ylm!(Y, L, cos(a), 0.0, P)
        flip = a > π
        for m in 0:L
            s = flip && isodd(m) ? -w : w
            for l in m:L
                κ[_zonal_index(l, m, L)] += s * Y[l+1, L+1+m]
            end
        end
    end
    return κ
end

# Gauss–Legendre nodes and weights on [-1, 1] (Newton on the three-term recurrence).
function _gauss_legendre(n)
    x, w = zeros(n), zeros(n)
    for i in 1:n
        z = cos(π * (i - 0.25) / (n + 0.5))
        dp = 0.0
        for _ in 1:100
            p0, p1 = 1.0, z
            for k in 2:n
                p0, p1 = p1, ((2k - 1) * z * p1 - (k - 1) * p0) / k
            end
            dp = n * (z * p1 - p0) / (z^2 - 1)
            δ = p1 / dp
            z -= δ
            abs(δ) < 1e-15 && break
        end
        x[i], w[i] = z, 2 / ((1 - z^2) * dp^2)
    end
    return x, w
end

# Blocks D_l of the rotation Q acting on real-harmonic coefficients: a field g in frame
# coordinates becomes f(x) = g(Qᵀx) in the lab, with lab coefficients D_l c_l (m = -l…l).
# D_ab = ∫ Y_a(Qx) Y_b(x) dΩ by a product quadrature exact to degree 2L.
function _rotation_blocks(Q, L)
    xg, wg = _gauss_legendre(L + 2)
    nφ = 2L + 2
    np = length(xg) * nφ
    Yf = [zeros(np, 2l + 1) for l in 0:L]
    Yl = [zeros(np, 2l + 1) for l in 0:L]
    wts = zeros(np)
    Y = zeros(L + 1, 2L + 1)
    P = zeros(L + 1, L + 1)
    p = 0
    for (μ, w) in zip(xg, wg), j in 1:nφ
        p += 1
        φ = 2π * (j - 1) / nφ
        wts[p] = w * 2π / nφ
        s = sqrt(1 - μ^2)
        _real_ylm!(Y, L, μ, φ, P)
        for l in 0:L
            Yf[l+1][p, :] .= @view Y[l+1, L+1-l:L+1+l]
        end
        x = Q * Vec3(s * cos(φ), s * sin(φ), μ)
        _real_ylm!(Y, L, x[3], atan(x[2], x[1]), P)
        for l in 0:L
            Yl[l+1][p, :] .= @view Y[l+1, L+1-l:L+1+l]
        end
    end
    return [Yl[l+1]' * (wts .* Yf[l+1]) for l in 0:L]
end

"""
    CoronaOperator(profile, beams; center = (0, 0, 0), lmax = 30, ncone = 13,
                   spacing = 0.1 * smallest spot semi-axis, panels = 200, shells = nothing)

Precomputed laser absorption in the static, spherically symmetric corona `profile`
([`CoronaProfile`](@ref)) about `center`, for polar and symmetric direct drive on a bare
sphere: `deposit(op, beams)` gives the absorbed power per steradian on the sphere about
`center` as real spherical harmonics up to degree `lmax`, the primary output of
direct-drive designs (`reference/litreview/polar_direct_drive.md`, decision 4).

Each beam is `ncone` collimated sub-beams along its f-cone directions
(`cone_directions`), each carrying the best-focus profile, as in
[`DepositionOperator`](@ref). Their rays are refracted by the corona (cold-plasma
geometric optics, ``n^2 = 1 - n_e/n_c``) and absorbed by inverse bremsstrahlung along their
paths (Kruer, *The Physics of Laser Plasma Interactions*, ch. 5); light that is not
absorbed leaves the corona. Absorbed power is assigned to the direction from `center` of
the point where it is absorbed (projected radially). In a spherical corona every ray
stays in its plane and is fixed by its impact parameter, so the paths are integrated
once, by radial quadrature (`panels` per pass), into harmonic kernels at the impact
parameters of a square lattice of rays (`spacing`, mm) across each sub-beam, out to the
profile's outer radius. These depend on the profile, the lattice and the ports only, so
the operator serves every aim, defocus, spot and power: the beams' ports and f-numbers
must not change after construction.

`shells`, increasing radii (mm) from at most the profile's first radius to at least its
outer radius, resolves the absorption in radius as well: `deposit` then also gives the
absorbed power per steradian in each shell `shells[s] ≤ r < shells[s+1]` and the
absorbed power per shell (a radial deposition profile), which sum to the totals. Each
shell edge is then also a panel boundary of every ray, so every panel lies in one shell
and the partition is exact; the totals then come from these refined panels (on an OMEGA
case the absorbed power equals the operator's without shells to round-off and the
harmonics by 1e-7 to 1e-6 of the mean, the panels' angular placement). A shell's
harmonics cost as much as the totals' in each evaluation, so the evaluation time grows
in proportion to the number of shells.
[`ablation_pressure_proxy`](@ref) sums the shells with a density weight.

Scope and limits:
- **Spherical corona only.** A non-spherical corona is outside the model: a P₂ of 4 % in
  the critical radius moves σ₂ by about 3 % of its value and the absorption by 1 point
  (measured by 3-D ray tracing; research notes, "Corona model prototype").
- **Bare sphere only.** Nothing but the corona interacts with the light: no other target
  parts, stalks or shields, which would need their own shadowing; the profile must reach
  ``n_c`` (`CoronaProfile`), so no light reaches the surface inside.
- **No cross-beam energy transfer** (decision 8 of the review; deferred): results are
  CBET-free, and optimistic where CBET matters (the PDD equator).

`deposit(op, beams; absorption = η⊥)` on a hard sphere is this model's limit for a
thin corona; oblique and grazing rays depart from it at finite scale lengths (research
notes, "Corona model prototype").
"""
struct CoronaOperator
    profile::CoronaProfile
    center::Vec3
    lmax::Int
    ports::Vector{Vec3}
    dirs::Vector{Vector{Vec3}}                 # f-cone directions per beam
    axes::Vector{Vector{NTuple{2,Vec3}}}       # per sub-beam: lattice axes e₁′, e₂′ ⟂ d
    euler::Vector{Vector{NTuple{2,Float64}}}   # per sub-beam: (α, β), frame = R_z(α) R_y(β)
    spacing::Float64
    n::Int                                     # lattice indices -n:n on each axis
    class::Matrix{Int32}                       # impact-parameter class of each lattice ray (0: misses)
    b::Vector{Float64}                         # impact parameter of each class (mm)
    kernel::Matrix{Float64}                    # zonal kernels (l, m ≥ 0) × class
    absorbed::Vector{Float64}                  # absorbed fraction of each class
    J::Vector{Matrix{Float64}}                 # rotation blocks of R_x(-90°)
    shells::Vector{Float64}                    # shell edges (mm); empty: no shells
    shell_kernel::Array{Float64,3}             # zonal kernels × class × shell
    radial::Matrix{Float64}                    # absorbed fraction of each class per shell
end

function CoronaOperator(profile::CoronaProfile, beams::AbstractVector{<:Beam};
                        center = (0.0, 0.0, 0.0), lmax::Integer = 30, ncone::Integer = 13,
                        spacing = nothing, panels::Integer = 200, shells = nothing)
    isempty(beams) && throw(ArgumentError("no beams"))
    lmax >= 0 && 1 <= ncone <= 16 && panels >= 1 ||
        throw(ArgumentError("need lmax ≥ 0, 1 ≤ ncone ≤ 16 and panels ≥ 1"))
    c = Vec3(center...)
    R = profile.r[end]
    edges = shells === nothing ? nothing : collect(Float64, shells)
    edges === nothing || (length(edges) >= 2 && all(>(0), diff(edges)) &&
                          edges[1] <= profile.r[1] && edges[end] >= R) ||
        throw(ArgumentError("shells must increase from at most the profile's first radius " *
                            "$(profile.r[1]) to at least its outer radius $R"))
    h = spacing === nothing ?
        0.1 * minimum(b -> ForwardDiff.value(min(b.spot.a, b.spot.b)), beams) : Float64(spacing)
    0 < h < R || throw(ArgumentError("spacing must be in (0, outer radius)"))
    n = floor(Int, R / h)
    # lattice rays that enter the corona (b < R), classed by i² + j²
    ks = sort!(unique!([i^2 + j^2 for i in -n:n for j in -n:n if h^2 * (i^2 + j^2) < R^2]))
    rank = Dict(k => Int32(q) for (q, k) in enumerate(ks))
    class = [h^2 * (i^2 + j^2) < R^2 ? rank[i^2+j^2] : Int32(0) for i in -n:n, j in -n:n]
    bs = h .* sqrt.(ks)
    L = Int(lmax)
    nz = (L + 1) * (L + 2) ÷ 2
    K = zeros(nz, length(bs))
    absorbed = zeros(length(bs))
    ns = edges === nothing ? 0 : length(edges) - 1
    SK = zeros(nz, length(bs), ns)
    radial = zeros(ns, length(bs))
    _tforeach(length(bs)) do q
        α, f, esc, _, sh = _ray_deposits(profile, bs[q], panels, edges)
        Y, P = zeros(L + 1, 2L + 1), zeros(L + 1, L + 1)
        _ray_kernel!(view(K, :, q), Y, P, L, α, f)
        absorbed[q] = 1 - esc
        for s in 1:ns
            sel = sh .== s
            any(sel) || continue
            _ray_kernel!(view(SK, :, q, s), Y, P, L, α[sel], f[sel])
            radial[s, q] = sum(f[sel])
        end
    end
    dirs = [cone_directions(b, ncone)[1] for b in beams]
    axes = [NTuple{2,Vec3}[] for _ in beams]
    euler = [NTuple{2,Float64}[] for _ in beams]
    for (bi, ds) in enumerate(dirs), d in ds
        z = -d
        α, β = atan(z[2], z[1]), acos(clamp(z[3], -1.0, 1.0))
        e1 = Vec3(cos(α) * cos(β), sin(α) * cos(β), -sin(β))   # R_z(α) R_y(β) x̂
        e2 = Vec3(-sin(α), cos(α), 0.0)                         # R_z(α) R_y(β) ŷ
        push!(axes[bi], (e1, e2))
        push!(euler[bi], (α, β))
    end
    J = _rotation_blocks(SMatrix{3,3,Float64}(1, 0, 0, 0, 0, -1, 0, 1, 0), L)   # R_x(-90°): ẑ ↦ ŷ
    return CoronaOperator(profile, c, L, [b.port for b in beams], dirs, axes, euler, h, n,
                          class, bs, K, absorbed, J, edges === nothing ? Float64[] : edges,
                          SK, radial)
end

function Base.show(io::IO, op::CoronaOperator)
    nr = count(>(0), op.class)
    @printf(io, "CoronaOperator(%d beams × %d directions, lattice %.4g mm: %d rays and %d impact parameters per direction, l ≤ %d)",
            length(op.ports), length(op.dirs[1]), op.spacing, nr, length(op.b), op.lmax)
end

"""
    deposit(op::CoronaOperator, beams) -> (; ylm, per_beam, delivered, shell_ylm, radial)

Laser absorption in the corona of `op` for `beams` (same ports and order as at
construction): `ylm`, the flat real-harmonic coefficients (`mode_index` order, l ≤
`op.lmax`, polar axis ẑ, about the operator's `center`) of the absorbed power per
steradian `F(Ω) = Σ ylm[mode_index(l, m)] Y_lm(Ω)` (TW/sr); `per_beam`, each beam's
absorbed power (TW); and `delivered`, the beams' total power. The absorbed fraction is
`sum(per_beam) / delivered` (equal to `√(4π) ylm[1] / delivered`), the nonuniformity
`mode_rms(ylm, lmax)` and its spectrum `mode_power(ylm_matrix(ylm, lmax))`. No
cross-beam energy transfer: results are CBET-free.

For an operator with `shells`, `shell_ylm[:, s]` holds the same coefficients for the
power absorbed in shell `s` (TW/sr; a `(lmax + 1)² × nshells` matrix, the absorbed power
per steradian and radius) and `radial[s]` the power absorbed in it (TW, the radial
deposition profile); `sum(eachcol(shell_ylm)) ≈ ylm` and `sum(radial) ≈ sum(per_beam)`
to round-off, and `√(4π) shell_ylm[1, s] ≈ radial[s]`. Without shells both are empty.

Generic in the beams' number type (aims, powers, defocus and spot axes may be
ForwardDiff duals); threads share the beams and the sums are taken in a fixed order, so
the result does not depend on the number of threads.
"""
function deposit(op::CoronaOperator, beams::AbstractVector{<:Beam})
    length(beams) == length(op.ports) || throw(DimensionMismatch("beam count changed"))
    for (b, p) in zip(beams, op.ports)
        b.port == p || throw(ArgumentError("beam $(b.name) port differs from the operator's"))
    end
    T = mapreduce(b -> typeof(b.power), promote_type, beams)
    nb = length(beams)
    frames = [BeamFrame(b) for b in beams]
    screens = [_screen_frame(f, b) for (f, b) in zip(frames, beams)]
    slopes = [_cone_slopes(op.dirs[bi], frames[bi]) for bi in 1:nb]
    nlm = (op.lmax + 1)^2
    ns = size(op.radial, 1)
    C = zeros(T, nlm, nb)
    A = zeros(T, nb)
    Cs = zeros(T, nlm, ns, nb)
    As = zeros(T, ns, nb)
    _corona_beams!(C, A, Cs, As, op, frames, screens, slopes)
    ylm = zeros(T, nlm)
    for bi in 1:nb, k in 1:nlm
        ylm[k] += C[k, bi]
    end
    shell_ylm = zeros(T, nlm, ns)
    radial = zeros(T, ns)
    for bi in 1:nb, s in 1:ns
        radial[s] += As[s, bi]
        for k in 1:nlm
            shell_ylm[k, s] += Cs[k, s, bi]
        end
    end
    return (; ylm, per_beam = A, delivered = convert(T, sum(b.power for b in beams)),
            shell_ylm, radial)
end

# Function barrier for the threads: each takes beams in turn and writes their columns.
function _corona_beams!(C::Matrix{T}, A, Cs, As, op, frames, screens, slopes) where {T}
    nb = length(frames)
    nt = min(Threads.nthreads(), nb)
    next = Threads.Atomic{Int}(1)
    if nt <= 1
        _corona_worker!(C, A, Cs, As, op, frames, screens, slopes, next)
    else
        @sync for _ in 1:nt
            Threads.@spawn _corona_worker!(C, A, Cs, As, op, frames, screens, slopes, next)
        end
    end
    return nothing
end

function _corona_worker!(C::Matrix{T}, A, Cs, As, op, frames, screens, slopes, next) where {T}
    L = op.lmax
    nlm = (L + 1)^2
    ws = (; Sc = zeros(T, L + 1, length(op.b)), Ss = zeros(T, L + 1, length(op.b)),
          cf = zeros(T, nlm), cfs = zeros(T, nlm, size(Cs, 2)), t1 = zeros(T, 2L + 1), t2 = zeros(T, 2L + 1),
          cosα = zeros(L + 1), sinα = zeros(L + 1), cosβ = zeros(L + 1), sinβ = zeros(L + 1))
    while (bi = Threads.atomic_add!(next, 1)) <= length(frames)
        A[bi] = _corona_beam!(view(C, :, bi), view(Cs, :, :, bi), view(As, :, bi), ws, op, bi,
                              frames[bi], screens[bi], slopes[bi])
    end
    return nothing
end

# Lattice index range on one axis that can reach the spot: a ray through lattice point
# (i, j) crosses the focal plane at (u, v), with h i = c₁ + u (e₁·e₁′) + v (e₂·e₁′), and
# the profile is zero beyond (u/a)² + (v/b)² = qcut, so |h i - c₁| ≤ √qcut ‖(a e₁·e₁′,
# b e₂·e₁′)‖ (Cauchy–Schwarz) for every live ray; one more index on each side for
# rounding. An exact screen: the rays outside contribute exact zeros.
function _lattice_range(f::BeamFrame{Float64}, c1, ep, h, n)
    p1, p2 = dot(f.e1, ep), dot(f.e2, ep)
    reach(s) = sqrt(s.qcut * (s.a2 * p1^2 + s.b2 * p2^2))
    ρ = f.has_halo ? max(reach(f.core), reach(f.halo)) : reach(f.core)
    return max(-n, floor(Int, (c1 - ρ) / h) - 1), min(n, ceil(Int, (c1 + ρ) / h) + 1)
end

# One beam: for each f-cone direction, the lattice rays' weights (best-focus profile at
# their focal-plane crossing × spacing² / (k·d) / ncone, the power of the collimated
# sub-beam per lattice cell normal to d) binned with cos mψ, sin mψ by impact-parameter
# class, contracted with the kernels into frame coefficients and rotated to the lab
# (D = Z(α) J Z(β) Jᵀ). Adds into `c`, and per shell into the columns of `cs` and the
# entries of `as`; returns the beam's absorbed power.
function _corona_beam!(c::AbstractVector{T}, cs, as, ws, op::CoronaOperator, bi, f::BeamFrame,
                       fv::BeamFrame{Float64}, (κ, α, β)) where {T}
    (; Sc, Ss, cf, cfs, t1, t2, cosα, sinα, cosβ, sinβ) = ws
    ns = size(cs, 2)
    L, n, h = op.lmax, op.n, op.spacing
    ds = op.dirs[bi]
    nd = length(ds)
    cls, K, ctr = op.class, op.kernel, op.center
    absorbed = zero(T)
    for j in eachindex(ds)
        e1p, e2p = op.axes[bi][j]
        ilo, ihi = _lattice_range(fv, dot(fv.focus - ctr, e1p), e1p, h, n)
        jlo, jhi = _lattice_range(fv, dot(fv.focus - ctr, e2p), e2p, h, n)
        (ilo > ihi || jlo > jhi) && continue
        scale = h^2 / (κ[j] * nd)
        lo, hi = typemax(Int32), Int32(0)
        @inbounds for jj in jlo:jhi, ii in ilo:ihi
            q = cls[ii+n+1, jj+n+1]
            q == 0 && continue
            x = ctr + (h * ii) * e1p + (h * jj) * e2p
            if T !== Float64           # screen on values: rays off the spot add exact zeros
                d = x - fv.focus
                s = dot(d, fv.k)
                u, v = dot(d, fv.e1) - s * α[j], dot(d, fv.e2) - s * β[j]
                _sg_inside(fv.core, u, v) || (fv.has_halo && _sg_inside(fv.halo, u, v)) || continue
            end
            d = x - f.focus
            s = dot(d, f.k)
            u, v = dot(d, f.e1) - s * α[j], dot(d, f.e2) - s * β[j]
            I = _sg_focal(f.core, u, v)
            f.has_halo && (I += _sg_focal(f.halo, u, v))
            iszero(I) && continue
            w = I * scale
            ρ = hypot(ii, jj)
            c1, s1 = ρ > 0 ? (ii / ρ, jj / ρ) : (1.0, 0.0)
            cm, sm = 1.0, 0.0
            for m in 0:L
                Sc[m+1, q] += w * cm
                Ss[m+1, q] += w * sm
                cm, sm = cm * c1 - sm * s1, sm * c1 + cm * s1
            end
            lo, hi = min(lo, q), max(hi, q)
        end
        lo > hi && continue
        # per shell (before the totals, which clear the moments)
        fill!(cfs, zero(T))
        @inbounds for s in 1:ns, q in lo:hi
            as[s] += Sc[1, q] * op.radial[s, q]
            for m in 0:L
                a, b = Sc[m+1, q], Ss[m+1, q]
                z = _zonal_index(m, m, L) - 1
                for l in m:L
                    k = op.shell_kernel[z+l-m+1, q, s]
                    cfs[mode_index(l, m), s] += k * a
                    m > 0 && (cfs[mode_index(l, -m), s] += k * b)
                end
            end
        end
        # frame coefficients: Σ over classes of kernel × azimuthal moments
        fill!(cf, zero(T))
        @inbounds for q in lo:hi
            absorbed += Sc[1, q] * op.absorbed[q]
            for m in 0:L
                a, b = Sc[m+1, q], Ss[m+1, q]
                z = _zonal_index(m, m, L) - 1
                for l in m:L
                    k = K[z+l-m+1, q]
                    cf[mode_index(l, m)] += k * a
                    m > 0 && (cf[mode_index(l, -m)] += k * b)
                end
                Sc[m+1, q] = zero(T)
                Ss[m+1, q] = zero(T)
            end
        end
        # to the lab frame
        αe, βe = op.euler[bi][j]
        for m in 0:L
            sinα[m+1], cosα[m+1] = sincos(m * αe)
            sinβ[m+1], cosβ[m+1] = sincos(m * βe)
        end
        _rotate_add!(c, cf, op.J, cosα, sinα, cosβ, sinβ, L, t1, t2)
        for s in 1:ns
            _rotate_add!(view(cs, :, s), view(cfs, :, s), op.J, cosα, sinα, cosβ, sinβ, L, t1, t2)
        end
    end
    return absorbed
end

# c += Z(α) J Z(β) Jᵀ cf, block by block (Z: rotation about ẑ, which mixes the
# cos mφ and sin mφ coefficients of each m).
function _rotate_add!(c::AbstractVector{T}, cf, J, cosα, sinα, cosβ, sinβ, L, t1, t2) where {T}
    @inbounds for l in 0:L
        o = l^2
        Jl = J[l+1]
        nl = 2l + 1
        for a in 1:nl                                  # t1 = Jᵀ cf
            acc = zero(T)
            for b in 1:nl
                acc += Jl[b, a] * cf[o+b]
            end
            t1[a] = acc
        end
        _zrotate!(t1, l, cosβ, sinβ)
        for a in 1:nl                                  # t2 = J t1
            t2[a] = zero(T)
        end
        for b in 1:nl
            tb = t1[b]
            for a in 1:nl
                t2[a] += Jl[a, b] * tb
            end
        end
        _zrotate!(t2, l, cosα, sinα)
        for a in 1:nl
            c[o+a] += t2[a]
        end
    end
    return c
end

@inline function _zrotate!(v, l, cs, sn)
    @inbounds for m in 1:l
        a, b = v[l+1+m], v[l+1-m]
        v[l+1+m] = cs[m+1] * a - sn[m+1] * b
        v[l+1-m] = sn[m+1] * a + cs[m+1] * b
    end
    return v
end

# ---------------------------------------------------------------------------------------
# Drive proxies (review decision 4): published, empirical transforms of the absorbed power
# per steradian into a drive. Each is nonlinear in the field, so it is applied on a
# product quadrature and projected back to harmonics.
# ---------------------------------------------------------------------------------------

# Applies `g(values, μ)` ring by ring to the band-limited field with coefficients `c`
# (l ≤ L) and returns the coefficients (l ≤ L) of the result. Gauss–Legendre in μ with
# 2L + 2 nodes and 4L + 4 uniform azimuths: the projection is exact for results up to
# degree 3L + 3 (higher degrees alias), and a ring's azimuthal mean is exact.
function _sphere_transform(g::G, c::AbstractVector, L::Integer) where {G}
    T = eltype(c)
    nθ, nφ = 2L + 2, 4L + 4
    μs, wμ = _gauss_legendre(nθ)
    cs = [cos(2π * m * (j - 1) / nφ) for m in 0:L, j in 1:nφ]
    sn = [sin(2π * m * (j - 1) / nφ) for m in 0:L, j in 1:nφ]
    Y, P = zeros(L + 1, 2L + 1), zeros(L + 1, L + 1)
    A, B, v = zeros(T, L + 1), zeros(T, L + 1), zeros(T, nφ)
    out = zeros(T, (L + 1)^2)
    for i in 1:nθ
        _real_ylm!(Y, L, μs[i], 0.0, P)        # Y[l+1, L+1+m]: the μ part of Y_l,±m
        fill!(A, zero(T))
        fill!(B, zero(T))
        for m in 0:L, l in m:L
            q = Y[l+1, L+1+m]
            A[m+1] += c[mode_index(l, m)] * q
            m > 0 && (B[m+1] += c[mode_index(l, -m)] * q)
        end
        for j in 1:nφ
            s = A[1]
            for m in 1:L
                s += A[m+1] * cs[m+1, j] + B[m+1] * sn[m+1, j]
            end
            v[j] = s
        end
        gv = g(v, μs[i])
        for m in 0:L
            a, b = zero(T), zero(T)
            for j in 1:nφ
                a += gv[j] * cs[m+1, j]
                b += gv[j] * sn[m+1, j]
            end
            a *= 2π / nφ * wμ[i]
            b *= 2π / nφ * wμ[i]
            for l in m:L
                q = Y[l+1, L+1+m]
                out[mode_index(l, m)] += a * q
                m > 0 && (out[mode_index(l, -m)] += b * q)
            end
        end
    end
    return out
end

"""
    ablation_pressure_proxy(op::CoronaOperator, dep; lmax = op.lmax, exponent = 2/3,
                            density_exponent = 2/3) -> Vector

The ablation-pressure proxy of Barlow et al. (PRL 133, 175101, 2024, Eqs. 2a–2b of the
published Letter's appendix; Eqs. 1a–1b of arXiv:2311.18342v2) from a corona deposition
`dep = deposit(op, beams)` of an operator with `shells`: the flat real-harmonic
coefficients (`mode_index` order, l ≤ `lmax`) of the pressure map (Mbar)

``P_{abl}(Ω) = 24.7\\,\\mathrm{Mbar}\\, (F(Ω) / 10^{14}\\,\\mathrm{W\\,cm^{-2}})^{2/3}``,
``F(Ω) = r_c^{-2} Σ_s (n_s/n_c)^{2/3}\\, dP_s/dΩ``,

with ``dP_s/dΩ`` the absorbed power per steradian in shell `s` (`dep.shell_ylm`), ``n_s``
the profile's electron density at the middle of the part of the shell inside the
profile, and ``r_c`` the profile's critical radius (outermost ``n_e = n_c``). Ba24 write
the sum over their hydrodynamic cells, ``Σ_r r^2 n_r^{2/3} I_r / (R_{351}^2 n_{351}^{2/3})``
with ``r^2 I_r`` the cell's absorbed power per steradian, so with a snapshot's cell edges
as `shells` this is their sum; finer shells converge to the radial integral, with an
error of second order in the shell width (on an OMEGA case, shells of an eighth of the
density scale length L out to 4L past the critical radius put σ₂ within 0.1 % of its
converged value, a quarter of L within 1 %; research notes). ``n_c`` is the critical
density of the profile's wavelength, Ba24's ``n_{351}`` at 351 nm. `exponent` replaces the outer 2/3
and `density_exponent` the density weight's 2/3.

The angular mean, Ba24's ``\\bar P_{abl}``, is `c[1] / √(4π)` and their ``σ_2`` is
`mode_rms(c, 30)`. The map is formed on a Gauss–Legendre product grid and projected back
to degree `lmax`, so it is ``F^{2/3}`` filtered to l ≤ `lmax` (as Ba24's σ₂ is), and
``F`` must be positive everywhere (an `ArgumentError` otherwise).

Provenance: an engineering proxy, not a facility criterion. The 2/3 power follows
Manheimer et al. (Phys. Fluids 25, 1644, 1982; Ba24's Ref. 33); the density weight is
empirical, chosen by Ba24 to match pressures and asymmetries of their 3-D
radiation-hydrodynamics runs; the prefactor 24.7 Mbar is Ba24's. CBET-free, as the
deposition. Generic in the number type (for ForwardDiff through `deposit`).
"""
function ablation_pressure_proxy(op::CoronaOperator, dep; lmax::Integer = op.lmax,
                                 exponent::Real = 2 / 3, density_exponent::Real = 2 / 3)
    isempty(op.shells) &&
        throw(ArgumentError("the operator has no shells: build it with CoronaOperator(...; shells)"))
    0 <= lmax <= op.lmax || throw(ArgumentError("need 0 ≤ lmax ≤ op.lmax = $(op.lmax)"))
    S = dep.shell_ylm
    ns = length(op.shells) - 1
    size(S, 2) == ns || throw(DimensionMismatch("dep has $(size(S, 2)) shells, op has $ns"))
    p = op.profile
    rc = _critical_radius(p)
    nlm = (lmax + 1)^2
    F = zeros(eltype(S), nlm)
    for s in 1:ns
        a, b = max(op.shells[s], p.r[1]), min(op.shells[s+1], p.r[end])
        a < b || continue                    # outside the profile: nothing absorbed
        w = _corona_ne(p, (a + b) / 2)^density_exponent / rc^2
        for k in 1:nlm
            F[k] += w * S[k, s]
        end
    end
    return _sphere_transform(F, lmax) do v, _
        minimum(v) > 0 ||
            throw(ArgumentError("the weighted absorbed intensity F is not positive everywhere"))
        24.7 .* v .^ exponent
    end
end

"""
    shell_displacement_proxy(ylm, lmax; exponent = 0.72, mean_displacement = nothing) -> Vector

The azimuthal shell-displacement model of Yeamans et al. (Nucl. Fusion 61, 046031, 2021,
Eq. 3) from the absorbed power per steradian `ylm` (flat real-harmonic coefficients,
`mode_index` order, polar axis ẑ; e.g. `deposit(op::CoronaOperator, beams).ylm`): the
coefficients (l ≤ `lmax`) of

``d(θ, φ) = \\bar d(θ)\\, (E(θ, φ) / \\bar E(θ))^{0.72}``,

with ``E`` the deposited energy per steradian (in a static deposition proportional to the
absorbed power per steradian), ``\\bar E(θ)`` its mean over φ at each polar angle θ, and
``\\bar d(θ)`` the mean distance the shell has moved at θ, `mean_displacement(θ)` (θ in
degrees) or 1 when `nothing`. Ye21 take ``\\bar d(θ)`` from a 2-D (axisymmetric)
hydrodynamics run; with ``\\bar d = 1`` the m ≠ 0 part, `mode_rms(d, lmax; m = :nonzero)`,
is the azimuthal nonuniformity of the shell velocity (their ``ΔV_{m≠0}``, which adds in
quadrature to the 2-D run's ``ΔV_{m=0}``). `exponent` replaces 0.72.

The printed Eq. 3 sets the factor 0.72 beside ``E/\\bar E`` without marking it as an
exponent; the text defines it as the power with which the shell velocity scales with
laser energy, which this function follows. The linear reading
``\\bar d\\,(1 + 0.72 (E/\\bar E - 1))`` agrees with it to first order in ``E/\\bar E - 1``,
and Ye21's statement that their velocity modes are 0.72 times the energy's holds
exactly for that reading and to first order for this one (research notes, "Corona
deposition per radius and drive proxies"). Formed on a Gauss–Legendre product grid and projected back to degree
`lmax`; ``E`` must be positive everywhere (an `ArgumentError` otherwise).

Provenance: an engineering proxy, not a facility criterion: the 0.72 power is from
implosion simulations (Ye21's Ref. 47, an LLE report), close to an experimental 0.71
that Ye21 cite.
Generic in the number type of `ylm`.
"""
function shell_displacement_proxy(ylm::AbstractVector, lmax::Integer; exponent::Real = 0.72,
                                  mean_displacement = nothing)
    lmax >= 0 && length(ylm) >= (lmax + 1)^2 ||
        throw(ArgumentError("coefficients reach l = $(isqrt(length(ylm)) - 1), not lmax = $lmax"))
    c = ylm[1:(lmax+1)^2]
    return _sphere_transform(c, lmax) do v, μ
        minimum(v) > 0 ||
            throw(ArgumentError("the absorbed power per steradian is not positive everywhere"))
        Ē = sum(v) / length(v)
        d̄ = mean_displacement === nothing ? 1.0 : mean_displacement(acosd(μ))
        d̄ .* (v ./ Ē) .^ exponent
    end
end
