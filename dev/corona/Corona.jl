# Corona prototype for polar direct drive, Stage 1 (decision 3 of
# reference/litreview/polar_direct_drive.md): two ways to deposit laser light in a static
# corona, to be compared on evidence. Research code, not package API: nothing here is
# exported or documented in the manual.
#
#   (b) `KernelTable` + `CoronaOperator`: a spherically symmetric profile n_e(r), with the
#       ray paths and inverse-bremsstrahlung absorption tabulated against impact parameter
#       by 1-D radial quadrature (no ray integration), then mapped through the samples of
#       a Briareus `DepositionOperator` built on the corona's outer sphere. Per pointing
#       change only the samples' beam weights are recomputed (ForwardDiff-generic).
#   (c) `trace_ray` / `trace_beams`: RK4 integration of the ray equations in 3-D through
#       any n_e(x) (here spherical, planar, or P2-perturbed), rays launched on a lattice
#       in each beam's transverse plane, retraced on every evaluation.
#
# Physics (both): geometric optics in a cold unmagnetized plasma, refractive index
# n² = 1 - n_e/n_c. With the ray parameter τ (ds = n dτ) the ray equations are
#     dx/dτ = v,   dv/dτ = -∇(n_e/n_c)/2,   |v|² = n²,
# and inverse-bremsstrahlung absorption (Kruer, The Physics of Laser Plasma Interactions,
# ch. 5: κ = ν_ei (n_e/n_c) / (c n)) is d ln P/dτ = -g(x), g = (ν_ei/c)(n_e/n_c), with
# ν_ei the local electron–ion collision frequency (∝ n_e Z lnΛ / T_e^{3/2}). For an
# isothermal profile g = (ν_c/c)(n_e/n_c)², ν_c the frequency at critical density.
# Lengths are mm, so ν/c is in mm⁻¹. Absorbed energy is projected radially onto the unit
# sphere (absorbed power per steradian, decision 4's primary output).

module Corona

using Briareus
using Briareus: Vec3, BeamFrame, _sg_focal, _cone_slopes, mode_index, _real_ylm!, cone_directions
using LinearAlgebra
using StaticArrays
using ForwardDiff

const V3 = SVector{3,Float64}

# ---------------------------------------------------------------------------------------
# Profiles. Radial ones give ne(p, r), dne(p, r), g(p, r) (n_e/n_c, its r-derivative and
# the absorption rate per unit τ, mm⁻¹), zero beyond p.Rout; 3-D ones give ne3, grad3, g3.
# ---------------------------------------------------------------------------------------

abstract type RadialProfile end

"Isothermal exponential corona n_e/n_c = exp(-(r - rc)/L), truncated at Rout."
struct ExpProfile <: RadialProfile
    rc::Float64
    L::Float64
    nuc::Float64     # ν_c / c (mm⁻¹)
    Rout::Float64
end
ne(p::ExpProfile, r) = r > p.Rout ? 0.0 : exp(-(r - p.rc) / p.L)
dne(p::ExpProfile, r) = r > p.Rout ? 0.0 : -exp(-(r - p.rc) / p.L) / p.L
g(p::ExpProfile, r) = p.nuc * ne(p, r)^2

"n_e/n_c = a/r² (a in mm²), isothermal: closed-form ray paths and absorption (validation)."
struct InvSqProfile <: RadialProfile
    a::Float64
    nuc::Float64
    Rout::Float64
end
ne(p::InvSqProfile, r) = r > p.Rout ? 0.0 : p.a / r^2
dne(p::InvSqProfile, r) = r > p.Rout ? 0.0 : -2p.a / r^3
g(p::InvSqProfile, r) = p.nuc * ne(p, r)^2

"""
Tabulated 1-D profile: ln(n_e/n_c) and ln(g) on a uniform radius grid, cubic Hermite
with centred-difference slopes (C¹, so ∇n_e is continuous for the ray tracer), zero
beyond the last node. This is the form a 1-D hydro snapshot would take.
"""
struct TableProfile <: RadialProfile
    r0::Float64
    dr::Float64
    y::Vector{Float64}     # ln ne
    z::Vector{Float64}     # ln g
    Rout::Float64
end
function TableProfile(f_ne, f_g, r0, Rout, n)
    r = range(r0, Rout; length = n)
    return TableProfile(r0, step(r), log.(f_ne.(r)), log.(f_g.(r)), Rout)
end
@inline function _hermite(y, r0, dr, r)
    n = length(y)
    s = (r - r0) / dr
    i = clamp(floor(Int, s), 0, n - 2)
    t = s - i
    i1 = i + 1
    m(k) = k == 1 ? (y[2] - y[1]) : k == n ? (y[n] - y[n-1]) : (y[k+1] - y[k-1]) / 2
    y0, y1, m0, m1 = y[i1], y[i1+1], m(i1), m(i1 + 1)
    h00, h10, h01, h11 = 2t^3 - 3t^2 + 1, t^3 - 2t^2 + t, -2t^3 + 3t^2, t^3 - t^2
    v = h00 * y0 + h10 * m0 + h01 * y1 + h11 * m1
    dv = ((6t^2 - 6t) * y0 + (3t^2 - 4t + 1) * m0 + (-6t^2 + 6t) * y1 + (3t^2 - 2t) * m1) / dr
    return v, dv
end
ne(p::TableProfile, r) = r > p.Rout ? 0.0 : exp(_hermite(p.y, p.r0, p.dr, r)[1])
function dne(p::TableProfile, r)
    r > p.Rout && return 0.0
    v, dv = _hermite(p.y, p.r0, p.dr, r)
    return exp(v) * dv
end
g(p::TableProfile, r) = r > p.Rout ? 0.0 : exp(_hermite(p.z, p.r0, p.dr, r)[1])

"General radial profile from functions (n_e/n_c and g of r); dn_e/dr by ForwardDiff."
struct FuncProfile{F,G} <: RadialProfile
    fne::F
    fg::G
    Rout::Float64
end
ne(p::FuncProfile, r) = r > p.Rout ? 0.0 : p.fne(r)
dne(p::FuncProfile, r) = r > p.Rout ? 0.0 : ForwardDiff.derivative(p.fne, r)
g(p::FuncProfile, r) = r > p.Rout ? 0.0 : p.fg(r)

# 3-D views of a radial profile
ne3(p::RadialProfile, x) = ne(p, norm(x))
grad3(p::RadialProfile, x) = (r = norm(x); dne(p, r) / r * x)
g3(p::RadialProfile, x) = g(p, norm(x))
outer(p::RadialProfile) = p.Rout

"""
Exponential corona whose critical radius carries a P2 perturbation about z:
n_e/n_c = exp(-(r - rc (1 + ε P₂(cos θ)))/L), isothermal. Not representable by (b).
"""
struct P2ExpProfile
    rc::Float64
    L::Float64
    nuc::Float64
    Rout::Float64
    ε::Float64
end
function ne3(p::P2ExpProfile, x)
    r = norm(x)
    r > p.Rout && return 0.0
    μ = x[3] / r
    return exp(-(r - p.rc * (1 + p.ε * (1.5μ^2 - 0.5))) / p.L)
end
function grad3(p::P2ExpProfile, x)
    r = norm(x)
    r > p.Rout && return zero(V3)
    μ = x[3] / r
    n = exp(-(r - p.rc * (1 + p.ε * (1.5μ^2 - 0.5))) / p.L)
    rhat = x / r
    gμ = (V3(0, 0, 1) - μ * rhat) / r
    return -n / p.L * (rhat - p.rc * p.ε * 3μ * gμ)
end
g3(p::P2ExpProfile, x) = p.nuc * ne3(p, x)^2
outer(p::P2ExpProfile) = p.Rout

"Planar profiles along x (validation): linear ramp n_e/n_c = x/L (x > 0) and exp(x/L)."
struct LinearSlab
    L::Float64
    nuc::Float64
end
ne3(p::LinearSlab, x) = x[1] > 0 ? x[1] / p.L : 0.0
grad3(p::LinearSlab, x) = x[1] > 0 ? V3(1 / p.L, 0, 0) : zero(V3)
g3(p::LinearSlab, x) = p.nuc * ne3(p, x)^2

struct ExpSlab
    L::Float64
    nuc::Float64
end
ne3(p::ExpSlab, x) = exp(x[1] / p.L)
grad3(p::ExpSlab, x) = V3(exp(x[1] / p.L) / p.L, 0, 0)
g3(p::ExpSlab, x) = p.nuc * ne3(p, x)^2

# ---------------------------------------------------------------------------------------
# (c) 3-D ray tracer: RK4 in τ on (x, v, ℓ = ∫ g dτ).
# ---------------------------------------------------------------------------------------

"""
    trace_ray(p, x, v, h, done; maxsteps) -> (x, v, ℓ, steps)

Integrate one ray from `x` with `v` (|v| = local n) by RK4 steps `h` in τ until
`done(x, v)`; calls `deposit(xmid, frac)` for the power fraction absorbed in each step
(fraction of the power at entry). Returns the final state; `exp(-ℓ)` is the fraction
that escapes.
"""
function trace_ray(deposit::F, p, x::V3, v::V3, h, done::D; maxsteps = 200_000) where {F,D}
    ℓ = 0.0
    P = 1.0
    n = 0
    while n < maxsteps
        n += 1
        k1x, k1v, k1l = v, -0.5 * grad3(p, x), g3(p, x)
        x2 = x + 0.5h * k1x
        k2x, k2v, k2l = v + 0.5h * k1v, -0.5 * grad3(p, x2), g3(p, x2)
        x3 = x + 0.5h * k2x
        k3x, k3v, k3l = v + 0.5h * k2v, -0.5 * grad3(p, x3), g3(p, x3)
        x4 = x + h * k3x
        k4x, k4v, k4l = v + h * k3v, -0.5 * grad3(p, x4), g3(p, x4)
        xn = x + h / 6 * (k1x + 2k2x + 2k3x + k4x)
        vn = v + h / 6 * (k1v + 2k2v + 2k3v + k4v)
        ℓn = ℓ + h / 6 * (k1l + 2k2l + 2k3l + k4l)
        Pn = exp(-ℓn)
        P - Pn > 0 && deposit(0.5 * (x + xn), P - Pn)
        x, v, ℓ, P = xn, vn, ℓn, Pn
        done(x, v) && break
        if P < 1e-10           # the rest is absorbed here (≤ 1e-10 of the ray)
            deposit(x, P)
            ℓ = Inf
            break
        end
    end
    return x, v, ℓ, n
end

# Entry point of the straight line x0 + t d on the sphere |x| = R (nothing if it misses).
function sphere_entry(x0::V3, d::V3, R)
    bq = dot(x0, d)
    c = dot(x0, x0) - R^2
    disc = bq^2 - c
    disc <= 0 && return nothing
    return x0 + (-bq - sqrt(disc)) * d
end

# Exit from a sphere-bounded corona: outside Rout and moving outward.
sphere_done(R) = (x, v) -> dot(x, x) > R^2 && dot(x, v) > 0

# ---------------------------------------------------------------------------------------
# Angular bins and spherical-harmonic spectrum (shared by (b) and (c)).
# ---------------------------------------------------------------------------------------

"Equal-angle (θ, φ) bins on the unit sphere about z, and Ylm tables at bin centres."
struct AngGrid
    nθ::Int
    nφ::Int
    lmax::Int
    Y0::Array{Float64,3}     # [l+1, lmax+1+m, iθ] real Ylm at φ = 0 (√2 included for m ≠ 0)
    cosmφ::Matrix{Float64}   # [m+1, iφ]
    sinmφ::Matrix{Float64}
end
function AngGrid(nθ, lmax)
    nφ = 2nθ
    Y = zeros(lmax + 1, 2lmax + 1)
    Y0 = zeros(lmax + 1, 2lmax + 1, nθ)
    for i in 1:nθ
        θ = (i - 0.5) * π / nθ
        _real_ylm!(Y, lmax, cos(θ), 0.0)
        Y0[:, :, i] .= Y
    end
    φ = [(j - 0.5) * 2π / nφ for j in 1:nφ]
    return AngGrid(nθ, nφ, lmax, Y0, [cos(m * f) for m in 0:lmax, f in φ],
                   [sin(m * f) for m in 0:lmax, f in φ])
end
nbins(G::AngGrid) = G.nθ * G.nφ
# Cloud-in-cell deposit of `P` at direction x onto the four nearest bin centres
# (bilinear in θ, φ): unlike nearest-bin assignment it has no first-order bias in the
# harmonic projection (nearest-bin was 4e-4 of the mean on a coarse (c) run).
@inline function deposit_cic!(E, G::AngGrid, x, P)
    r = norm(x)
    θ = acos(clamp(x[3] / r, -1.0, 1.0))
    φ = mod(atan(x[2], x[1]), 2π)
    sθ = θ / π * G.nθ - 0.5
    sφ = φ / 2π * G.nφ - 0.5
    i0 = floor(Int, sθ); tθ = sθ - i0
    j0 = floor(Int, sφ); tφ = sφ - j0
    i1 = i0 + 1
    i0 < 0 && (i0 = 0); i1 > G.nθ - 1 && (i1 = G.nθ - 1)
    j1 = mod(j0 + 1, G.nφ); j0 = mod(j0, G.nφ)
    @inbounds begin
        E[i0*G.nφ+j0+1] += P * (1 - tθ) * (1 - tφ)
        E[i0*G.nφ+j1+1] += P * (1 - tθ) * tφ
        E[i1*G.nφ+j0+1] += P * tθ * (1 - tφ)
        E[i1*G.nφ+j1+1] += P * tθ * tφ
    end
    return E
end

@inline function binindex(G::AngGrid, x)
    r = norm(x)
    θ = acos(clamp(x[3] / r, -1.0, 1.0))
    φ = mod(atan(x[2], x[1]), 2π)
    i = clamp(floor(Int, θ / π * G.nθ), 0, G.nθ - 1)
    j = clamp(floor(Int, φ / 2π * G.nφ), 0, G.nφ - 1)
    return i * G.nφ + j + 1
end

"Flat real-harmonic coefficients (Briareus `mode_index` order) of binned power `E`."
function ylm_coeffs(E::AbstractVector{T}, G::AngGrid) where {T}
    L = G.lmax
    c = zeros(T, (L + 1)^2)
    Fc = zeros(T, L + 1)
    Fs = zeros(T, L + 1)
    for i in 1:G.nθ
        fill!(Fc, zero(T)); fill!(Fs, zero(T))
        off = (i - 1) * G.nφ
        for j in 1:G.nφ
            e = E[off+j]
            iszero(e) && continue
            for m in 0:L
                Fc[m+1] += e * G.cosmφ[m+1, j]
                Fs[m+1] += e * G.sinmφ[m+1, j]
            end
        end
        for l in 0:L
            c[mode_index(l, 0)] += G.Y0[l+1, L+1, i] * Fc[1]
            for m in 1:l
                c[mode_index(l, m)] += G.Y0[l+1, L+1+m, i] * Fc[m+1]
                c[mode_index(l, -m)] += G.Y0[l+1, L+1+m, i] * Fs[m+1]
            end
        end
    end
    return c
end

"σ_l (l = 0…lmax, index l+1) from flat coefficients."
sigma_l(c, L) = [sqrt(sum(abs2, c[mode_index(l, m)] for m in -l:l)) / abs(c[1]) for l in 0:L]
sigma_rms(c, L) = sqrt(sum(abs2, c[mode_index(l, m)] for l in 1:L for m in -l:l)) / abs(c[1])

# Direct (unbinned) accumulation: c_lm += P Y_lm(x̂). Exact for point deposits.
function accumulate_ylm!(c, Y, Pw, L, x, P)
    r = norm(x)
    _real_ylm!(Y, L, x[3] / r, atan(x[2], x[1]), Pw)
    for l in 0:L, m in -l:l
        c[mode_index(l, m)] += P * Y[l+1, L+1+m]
    end
    return c
end

# ---------------------------------------------------------------------------------------
# (b) Impact-parameter kernels by radial quadrature.
# ---------------------------------------------------------------------------------------

const _GL4 = ((-0.8611363115940526, 0.3478548451374538), (-0.3399810435848563, 0.6521451548625461),
              (0.3399810435848563, 0.6521451548625461), (0.8611363115940526, 0.3478548451374538))

# Turning radius: the largest r < Rout with r² n²(r) = b².
function turning_radius(p::RadialProfile, b; scan = 200)
    R = p.Rout
    h(r) = r^2 * (1 - ne(p, r)) - b^2
    hi = R * (1 - 1e-12)
    lo = hi
    step = R / scan
    while h(lo) > 0
        hi = lo
        lo -= step
        lo <= 0 && (lo = 0.0; break)
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        h(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    return hi
end

"""
    ray_kernel(p, b; N) -> (; rt, αin, A, T)

Planar path of a ray with impact parameter `b` in a radial profile: turning radius,
entry angle on the Rout sphere (from -k toward the impact vector), and per panel
(r = rt + (Rout - rt) t², N uniform panels in t ∈ [0, 1], 4-point Gauss each) the swept
angle A and the absorption depth T. Integrable 1/√ singularity at rt removed by the t²
substitution.
"""
function ray_kernel(p::RadialProfile, b; N = 400)
    R = p.Rout
    rt = turning_radius(p, b)
    Δ = R - rt
    A = zeros(N)
    T = zeros(N)
    for k in 1:N
        t0, t1 = (k - 1) / N, k / N
        a = 0.0; τ = 0.0
        for (ξ, w) in _GL4
            t = (t0 + t1) / 2 + (t1 - t0) / 2 * ξ
            r = rt + Δ * t^2
            dr = 2Δ * t * (t1 - t0) / 2 * w
            q = r^2 * (1 - ne(p, r)) - b^2
            q <= 0 && continue
            s = sqrt(q)
            a += b * dr / (r * s)
            τ += g(p, r) * r * dr / s
        end
        A[k] = a; T[k] = τ
    end
    return (; rt, αin = asin(min(b / R, 1.0)), A, T)
end

"""
Deposit nodes of a kernel: 2N (α, fraction) pairs, inward panels N…1 then outward 1…N,
each panel's absorbed fraction at its mid angle; `escaped` = exp(-total depth).
"""
function kernel_nodes(k)
    N = length(k.A)
    α = zeros(2N); f = zeros(2N)
    a = k.αin; τ = 0.0
    for (i, kk) in enumerate(N:-1:1)
        α[i] = a + k.A[kk] / 2
        f[i] = exp(-τ) * (-expm1(-k.T[kk]))
        a += k.A[kk]; τ += k.T[kk]
    end
    for kk in 1:N
        i = N + kk
        α[i] = a + k.A[kk] / 2
        f[i] = exp(-τ) * (-expm1(-k.T[kk]))
        a += k.A[kk]; τ += k.T[kk]
    end
    return α, f, exp(-τ), a     # a: total swept angle (exit angle from -k)
end

"Kernel table: deposit nodes on a uniform b grid in [0, Rout)."
struct KernelTable
    b::Vector{Float64}
    α::Matrix{Float64}    # [node, ib]
    f::Matrix{Float64}
    escaped::Vector{Float64}
    Rout::Float64
end
function KernelTable(p::RadialProfile; Nb = 1024, N = 200)
    R = p.Rout
    bs = collect(range(0, R * (1 - 1e-9); length = Nb))
    α = zeros(2N, Nb); f = zeros(2N, Nb); esc = zeros(Nb)
    Threads.@threads for i in 1:Nb
        a, ff, e, _ = kernel_nodes(ray_kernel(p, bs[i]; N))
        α[:, i] .= a; f[:, i] .= ff; esc[i] = e
    end
    return KernelTable(bs, α, f, esc, R)
end
# Linear interpolation weights in b.
@inline function _bslot(K::KernelTable, b)
    s = b / (K.b[2] - K.b[1])
    i = clamp(floor(Int, s), 0, length(K.b) - 2)
    return i + 1, s - i
end
"Absorbed fraction for impact parameter b (interpolated)."
function absorbed(K::KernelTable, b)
    b >= K.Rout && return 0.0
    i, t = _bslot(K, b)
    return 1 - ((1 - t) * K.escaped[i] + t * K.escaped[i+1])
end

# ---------------------------------------------------------------------------------------
# (b) operator, harmonic form (the one used): every ray of a sub-beam (one beam, one
# f-cone direction d) deposits along the meridian of its impact azimuth ψ in the frame
# z' = -d, at polar angles α (the kernel nodes). Its real-harmonic coefficients in that
# frame are κ_lm(b) cos mψ and κ_lm(b) sin mψ (m ≥ 0), κ tabulated against b once per
# profile. A sub-beam's coefficients are Σ_b κ_lm(b) S_m(b), with S_m(b) the beam
# weights' azimuthal moments binned (linearly) in b, and a fixed rotation per sub-beam
# takes them to the lab frame. Exact harmonics (no angular bins); memory per sample is a
# few numbers. Generic in the number type of the beams.
# ---------------------------------------------------------------------------------------

"κ[l+1, m+1, ib]: zonal-frame harmonic coefficients of the kernel nodes at each b node."
function kappa_table(K::KernelTable, lmax)
    Nb = length(K.b)
    κ = zeros(lmax + 1, lmax + 1, Nb)
    Threads.@threads for i in 1:Nb
        Y = zeros(lmax + 1, 2lmax + 1); P = zeros(lmax + 1, lmax + 1)
        for k in axes(K.α, 1)
            f = K.f[k, i]
            f > 0 || continue
            a = mod(K.α[k, i], 2π)
            θ, φ = a <= π ? (a, 0.0) : (2π - a, π)
            _real_ylm!(Y, lmax, cos(θ), φ, P)
            for l in 0:lmax, m in 0:l
                κ[l+1, m+1, i] += f * Y[l+1, lmax+1+m]
            end
        end
    end
    return κ
end
# Single-node models (hard sphere): node (α(b), f(b)) on the same b grid.
function kappa_table(bs::AbstractVector, model, lmax)
    κ = zeros(lmax + 1, lmax + 1, length(bs))
    Y = zeros(lmax + 1, 2lmax + 1); P = zeros(lmax + 1, lmax + 1)
    for (i, b) in enumerate(bs)
        α, f = model(b)
        for (a, w) in zip(α, f)
            _real_ylm!(Y, lmax, cos(a), 0.0, P)
            for l in 0:lmax, m in 0:l
                κ[l+1, m+1, i] += w * Y[l+1, lmax+1+m]
            end
        end
    end
    return κ
end

# Gauss–Legendre nodes on [-1, 1]
function gauss(n)
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

"""
Rotation of real-harmonic coefficients from the frame with axes Q[:,1..3] to the lab:
blocks D[l+1] with c_lab[l] = D[l+1] * c_frame[l] (m = -l…l order), by a product
quadrature exact to degree 2lmax.
"""
function rotation_blocks(Q::AbstractMatrix, lmax)
    xg, wg = gauss(lmax + 2)
    nφ = 2lmax + 2
    np = length(xg) * nφ
    Yf = [zeros(np, 2l + 1) for l in 0:lmax]   # Y_lm'(x') (frame)
    Yl = [zeros(np, 2l + 1) for l in 0:lmax]   # Y_lm(Q x') (lab)
    wts = zeros(np)
    Y = zeros(lmax + 1, 2lmax + 1); P = zeros(lmax + 1, lmax + 1)
    p = 0
    for (μ, w) in zip(xg, wg), j in 1:nφ
        p += 1
        φ = 2π * (j - 1) / nφ
        wts[p] = w * 2π / nφ
        s = sqrt(1 - μ^2)
        xp = V3(s * cos(φ), s * sin(φ), μ)
        _real_ylm!(Y, lmax, μ, φ, P)
        for l in 0:lmax
            Yf[l+1][p, :] .= @view Y[l+1, lmax+1-l:lmax+1+l]
        end
        x = Q * xp
        _real_ylm!(Y, lmax, x[3], atan(x[2], x[1]), P)
        for l in 0:lmax
            Yl[l+1][p, :] .= @view Y[l+1, lmax+1-l:lmax+1+l]
        end
    end
    return [Yl[l+1]' * (wts .* Yf[l+1]) for l in 0:lmax]
end

struct HarmonicOperator
    op::DepositionOperator
    lmax::Int
    db::Float64                       # b-grid spacing (grid from 0)
    Nb::Int
    κ::Array{Float64,3}
    # per sub-beam s = (bi, j): pairs as (point q, b slot ib, weight t, cos ψ, sin ψ)
    sub::Vector{Tuple{Int,Int}}
    q::Vector{Vector{Int32}}
    ib::Vector{Vector{Int32}}
    t::Vector{Vector{Float64}}
    cψ::Vector{Vector{Float64}}
    sψ::Vector{Vector{Float64}}
    D::Vector{Vector{Matrix{Float64}}}
end

function HarmonicOperator(op::DepositionOperator, bs::AbstractVector, κ::Array{Float64,3}, lmax)
    db = bs[2] - bs[1]
    Nb = length(bs)
    sub = [(bi, j) for bi in eachindex(op.ports) for j in eachindex(op.dirs[bi])]
    ns = length(sub)
    Qv = Vector{Vector{Int32}}(undef, ns); IB = Vector{Vector{Int32}}(undef, ns)
    TT = Vector{Vector{Float64}}(undef, ns); CS = Vector{Vector{Float64}}(undef, ns)
    SN = Vector{Vector{Float64}}(undef, ns); D = Vector{Vector{Matrix{Float64}}}(undef, ns)
    Threads.@threads for s in 1:ns
        bi, j = sub[s]
        d = V3(op.dirs[bi][j])
        z = -d
        x1 = normalize(abs(z[3]) < 0.9 ? cross(V3(0, 0, 1), z) : cross(V3(1, 0, 0), z))
        y1 = cross(z, x1)
        Q = hcat(x1, y1, z)
        D[s] = rotation_blocks(Q, lmax)
        qs = Int32[]; ibs = Int32[]; ts = Float64[]; cs = Float64[]; sn = Float64[]
        for q in eachindex(op.pts[bi])
            op.dest[bi][q] > 0 || continue
            (op.mask[bi][q] >> (j - 1)) & 0x1 == 0x1 || continue
            x = V3(op.pts[bi][q])
            bv = x - dot(x, d) * d
            b = norm(bv)
            sl = b / db
            i = min(floor(Int, sl), Nb - 2)
            ψ = atan(dot(bv, y1), dot(bv, x1))
            push!(qs, q); push!(ibs, i + 1); push!(ts, sl - i); push!(cs, cos(ψ)); push!(sn, sin(ψ))
        end
        Qv[s], IB[s], TT[s], CS[s], SN[s] = qs, ibs, ts, cs, sn
    end
    return HarmonicOperator(op, lmax, db, Nb, κ, sub, Qv, IB, TT, CS, SN, D)
end

npairs(H::HarmonicOperator) = sum(length, H.q)
bytes(H::HarmonicOperator) = sizeof(H.κ) + sum(v -> sum(sizeof, v), (H.q, H.ib, H.t, H.cψ, H.sψ)) +
                             sum(s -> sum(sizeof, s), H.D)

"Flat lab-frame coefficients (mode_index order) of absorbed power, and incident power."
function apply(H::HarmonicOperator, beams)
    op, L = H.op, H.lmax
    f1 = BeamFrame(beams[1])
    T = promote_type(eltype(f1.focus), typeof(f1.core.peak))
    ns = length(H.sub)
    nt = Threads.nthreads()
    cbuf = [zeros(T, (L + 1)^2) for _ in 1:nt]
    ibuf = [zero(T) for _ in 1:nt]
    frames = [BeamFrame(b) for b in beams]
    slopes = [_cone_slopes(op.dirs[bi], frames[bi]) for bi in eachindex(beams)]
    next = Threads.Atomic{Int}(1)
    Threads.@threads for tid in 1:nt
        Sc = zeros(T, H.Nb, L + 1); Ss = zeros(T, H.Nb, L + 1)
        cm = zeros(T, L + 1); sm = zeros(T, L + 1)
        cf = zeros(T, (L + 1)^2)
        cc = cbuf[tid]
        while true
            s = Threads.atomic_add!(next, 1)
            s > ns && break
            bi, j = H.sub[s]
            f = frames[bi]
            κc, α, β = slopes[bi]
            dj = op.dirs[bi][j]
            P, W, Nr = op.pts[bi], op.w[bi], op.normal[bi]
            fill!(Sc, zero(T)); fill!(Ss, zero(T))
            lo, hi = H.Nb, 1
            inc = zero(T)
            for k in eachindex(H.q[s])
                q = H.q[s][k]
                d = P[q] - f.focus
                sd = dot(d, f.k)
                u = dot(d, f.e1) - sd * α[j]
                v = dot(d, f.e2) - sd * β[j]
                I = _sg_focal(f.core, u, v)
                f.has_halo && (I += _sg_focal(f.halo, u, v))
                iszero(I) && continue
                w = W[q] * I * abs(dot(Nr[q], dj)) / κc[j] / length(op.dirs[bi])
                inc += w
                ib, t, c1, s1 = H.ib[s][k], H.t[s][k], H.cψ[s][k], H.sψ[s][k]
                lo = min(lo, ib); hi = max(hi, ib + 1)
                # cos mψ, sin mψ by recurrence
                cprev, sprev, ccur, scur = 1.0, 0.0, 1.0, 0.0
                w0, w1 = (1 - t) * w, t * w
                for m in 0:L
                    if m == 1
                        ccur, scur = c1, s1
                    elseif m > 1
                        ccur, scur = ccur * c1 - scur * s1, scur * c1 + ccur * s1
                    end
                    Sc[ib, m+1] += w0 * ccur; Sc[ib+1, m+1] += w1 * ccur
                    Ss[ib, m+1] += w0 * scur; Ss[ib+1, m+1] += w1 * scur
                end
            end
            ibuf[tid] += inc
            lo > hi && continue
            # frame coefficients
            fill!(cf, zero(T))
            for i in lo:hi, m in 0:L
                a, bsn = Sc[i, m+1], Ss[i, m+1]
                (iszero(a) && iszero(bsn)) && continue
                for l in m:L
                    κ = H.κ[l+1, m+1, i]
                    cf[mode_index(l, m)] += κ * a
                    m > 0 && (cf[mode_index(l, -m)] += κ * bsn)
                end
            end
            # rotate to the lab frame
            for l in 0:L
                Dl = H.D[s][l+1]
                r = mode_index(l, -l):mode_index(l, l)
                for a in 1:2l+1
                    acc = zero(T)
                    for b in 1:2l+1
                        acc += Dl[a, b] * cf[r[b]]
                    end
                    cc[r[a]] += acc
                end
            end
        end
    end
    return sum(cbuf), sum(ibuf)
end

# ---------------------------------------------------------------------------------------
# (b) operator, binned form (kept for the record: mapping each sample's kernel onto
# 0.5° bins costs ~130 nonzeros per sample, 6.5 GB on the OMEGA case).
# ---------------------------------------------------------------------------------------

struct CoronaOperator
    op::DepositionOperator
    ncone::Int
    pair_q::Vector{Vector{Int32}}       # per beam: point index of each pair
    pair_j::Vector{Vector{UInt8}}       # per beam: cone direction of each pair
    ptr::Vector{Vector{Int}}            # per beam: CSC-like column pointers (length npairs+1)
    bin::Vector{Vector{Int32}}
    val::Vector{Vector{Float64}}
end

# Bins of one ray (entry x on the Rout sphere, direction d) through the kernel table;
# `nodes(b) -> (α, f)` lets the same mapper serve the hard-sphere η_a model.
function ray_bins!(acc::Dict{Int32,Float64}, G::AngGrid, x::V3, d::V3, α, f)
    bv = x - dot(x, d) * d
    b = norm(bv)
    eb = b > 0 ? bv / b : (abs(d[3]) < 0.9 ? normalize(cross(d, V3(0, 0, 1))) : normalize(cross(d, V3(1, 0, 0))))
    for (a, w) in zip(α, f)
        w > 0 || continue
        dir = -cos(a) * d + sin(a) * eb
        k = Int32(binindex(G, dir))
        acc[k] = get(acc, k, 0.0) + w
    end
    return b
end

function kernel_at(K::KernelTable, b)
    b >= K.Rout && return Float64[], Float64[]
    i, t = _bslot(K, b)
    return (1 - t) .* view(K.α, :, i) .+ t .* view(K.α, :, i + 1),
           (1 - t) .* view(K.f, :, i) .+ t .* view(K.f, :, i + 1)
end

"""
    CoronaOperator(op, ncone, G, model)

Map every (point, live cone direction) sample of `op` (built on the Rout sphere) to
angular bins. `model(b) -> (α, f)` gives the deposit nodes for impact parameter b:
`b -> kernel_at(K, b)` for the corona, or a hard-sphere model.
"""
function CoronaOperator(op::DepositionOperator, ncone::Integer, G::AngGrid, model)
    nb = length(op.ports)
    PQ = Vector{Vector{Int32}}(undef, nb); PJ = Vector{Vector{UInt8}}(undef, nb)
    PTR = Vector{Vector{Int}}(undef, nb); BIN = Vector{Vector{Int32}}(undef, nb)
    VAL = Vector{Vector{Float64}}(undef, nb)
    Threads.@threads for bi in 1:nb
        pq = Int32[]; pj = UInt8[]; ptr = [1]; bin = Int32[]; val = Float64[]
        acc = Dict{Int32,Float64}()
        ds = op.dirs[bi]
        for q in eachindex(op.pts[bi])
            op.dest[bi][q] > 0 || continue
            m = op.mask[bi][q]
            x = V3(op.pts[bi][q])
            for j in eachindex(ds)
                (m >> (j - 1)) & 0x1 == 0x1 || continue
                d = V3(ds[j])
                bv = x - dot(x, d) * d
                α, f = model(norm(bv))
                empty!(acc)
                ray_bins!(acc, G, x, d, α, f)
                push!(pq, q); push!(pj, j)
                for (k, w) in acc
                    push!(bin, k); push!(val, w)
                end
                push!(ptr, length(bin) + 1)
            end
        end
        PQ[bi], PJ[bi], PTR[bi], BIN[bi], VAL[bi] = pq, pj, ptr, bin, val
    end
    return CoronaOperator(op, ncone, PQ, PJ, PTR, BIN, VAL)
end

npairs(C::CoronaOperator) = sum(length, C.pair_q)
nnz(C::CoronaOperator) = sum(length, C.bin)
bytes(C::CoronaOperator) = sum(sizeof, C.bin) + sum(sizeof, C.val) + sum(sizeof, C.ptr) +
                           sum(sizeof, C.pair_q) + sum(sizeof, C.pair_j)

"Incident power of each pair of beam bi (generic in the beams' number type)."
function pair_weights(C::CoronaOperator, bi, b::Beam)
    op = C.op
    f = BeamFrame(b)
    κ, α, β = _cone_slopes(op.dirs[bi], f)
    P, W, Nr = op.pts[bi], op.w[bi], op.normal[bi]
    T = promote_type(eltype(f.focus), typeof(f.core.peak))
    w = zeros(T, length(C.pair_q[bi]))
    for (i, (q, j)) in enumerate(zip(C.pair_q[bi], C.pair_j[bi]))
        d = P[q] - f.focus
        s = dot(d, f.k)
        u = dot(d, f.e1) - s * α[j]
        v = dot(d, f.e2) - s * β[j]
        I = _sg_focal(f.core, u, v)
        f.has_halo && (I += _sg_focal(f.halo, u, v))
        w[i] = W[q] * I * abs(dot(Nr[q], op.dirs[bi][j])) / κ[j] / length(op.dirs[bi])
    end
    return w
end

"Absorbed power per bin (and incident power) for `beams`; ForwardDiff-generic."
function apply(C::CoronaOperator, G::AngGrid, beams)
    T = (f1 = BeamFrame(beams[1]); promote_type(eltype(f1.focus), typeof(f1.core.peak)))
    nb = length(beams)
    nt = Threads.nthreads()
    bufs = [zeros(T, nbins(G)) for _ in 1:nt]
    inc = zeros(T, nb)
    chunks = collect(Iterators.partition(1:nb, cld(nb, nt)))
    Threads.@threads for c in eachindex(chunks)
        E = bufs[c]
        for bi in chunks[c]
            w = pair_weights(C, bi, beams[bi])
            inc[bi] = sum(w)
            ptr, bin, val = C.ptr[bi], C.bin[bi], C.val[bi]
            for i in eachindex(w)
                wi = w[i]
                iszero(wi) && continue
                for k in ptr[i]:ptr[i+1]-1
                    E[bin[k]] += wi * val[k]
                end
            end
        end
    end
    E = bufs[1]
    for c in 2:length(chunks)
        E .+= bufs[c]
    end
    return E, sum(inc)
end

# ---------------------------------------------------------------------------------------
# (c) on beams: a lattice of parallel rays per cone direction in the plane normal to the
# beam through the origin (fixed by the port; aim-independent positions, so the weights
# alone carry the beams' parameters and ForwardDiff goes through them).
# ---------------------------------------------------------------------------------------

"""
    trace_beams(p, beams, G; ncone, hlat, hτ) -> (E, incident, nrays)

Absorbed power per bin. Rays: lattice spacing `hlat` (mm), the f-cone sampled by
Briareus's `cone_directions(beam, ncone)`; RK4 step `hτ`. Each ray's weight is the
best-focus profile at its focal-plane crossing times hlat²/ncone. `E` has the beams'
number type.
"""
function trace_beams(p, beams, G::AngGrid; ncone = 7, hlat = 0.02, hτ = 0.004,
                     direct_ylm = false)
    R = outer(p)
    done = sphere_done(R)
    T = (f1 = BeamFrame(beams[1]); promote_type(eltype(f1.focus), typeof(f1.core.peak)))
    # work list: (beam, direction, lattice point) on values
    work = Tuple{Int,Int,V3,V3}[]
    for (bi, b) in enumerate(beams)
        f = BeamFrame(b)
        ds, _ = cone_directions(b, ncone)
        c0 = ForwardDiff.value.(f.focus)
        reach = sqrt(f.core.qcut) * sqrt(max(ForwardDiff.value(f.core.a2), ForwardDiff.value(f.core.b2))) +
                R * tan(atan(1 / (2b.fnumber))) + norm(c0)
        n = ceil(Int, reach / hlat)
        for (j, d) in enumerate(ds), iu in -n:n, iv in -n:n
            x0 = V3(iu * hlat * f.e1 + iv * hlat * f.e2)
            sphere_entry(x0, V3(d), R) === nothing && continue
            push!(work, (bi, j, x0, V3(d)))
        end
    end
    nt = Threads.nthreads()
    bufs = [zeros(T, nbins(G)) for _ in 1:nt]
    cl = [zeros(T, (G.lmax + 1)^2) for _ in 1:nt]
    incs = [zero(T) for _ in 1:nt]
    live = zeros(Int, nt)
    frames = [BeamFrame(b) for b in beams]
    slopes = [_cone_slopes(cone_directions(b, ncone)[1], f) for (b, f) in zip(beams, frames)]
    chunks = collect(Iterators.partition(eachindex(work), cld(length(work), 8nt)))
    lk = ReentrantLock()
    next = Threads.Atomic{Int}(1)
    Threads.@threads for t in 1:nt
        E = bufs[t]; cc = cl[t]
        Y = zeros(G.lmax + 1, 2G.lmax + 1); Pw = zeros(G.lmax + 1, G.lmax + 1)
        while true
            c = Threads.atomic_add!(next, 1)
            c > length(chunks) && break
            for w in chunks[c]
                bi, j, x0, d = work[w]
                f = frames[bi]
                κ, α, β = slopes[bi]
                dd = x0 - f.focus
                s = dot(dd, f.k)
                u = dot(dd, f.e1) - s * α[j]
                v = dot(dd, f.e2) - s * β[j]
                I = _sg_focal(f.core, u, v)
                f.has_halo && (I += _sg_focal(f.halo, u, v))
                iszero(I) && continue
                wt = I * hlat^2 / ncone
                incs[t] += wt
                live[t] += 1
                xe = sphere_entry(x0, d, R)
                trace_ray(p, xe, d, hτ, done) do xm, fr
                    deposit_cic!(E, G, xm, wt * fr)
                    direct_ylm && accumulate_ylm!(cc, Y, Pw, G.lmax, xm, ForwardDiff.value(wt) * fr)
                end
            end
        end
    end
    E = bufs[1]
    for t in 2:nt
        E .+= bufs[t]
    end
    return E, sum(incs), sum(live), sum(cl)
end

end # module
