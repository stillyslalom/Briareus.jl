# 3-D ray tracing through a static corona: the prototype's model (c), kept as an offline
# study tool for coronas the package's `CoronaOperator` (model (b), spherical only)
# cannot represent, such as a P2-perturbed critical surface or a non-spherical hydro
# snapshot (decision 3 of reference/litreview/polar_direct_drive.md; research notes,
# "Corona model prototype" and "Corona absorption operator"). Research code, not package
# API. The prototype's (b) moved into src/corona.jl; a spherical-exponential copy of this
# tracer guards it in test/corona_raytrace.jl. p2.jl runs the P2 study.
#
#   `trace_ray` / `trace_beams`: RK4 integration of the ray equations in 3-D through
#   any n_e(x) (spherical, tabulated, planar, or P2-perturbed), rays launched on a
#   lattice in each beam's transverse plane, retraced on every evaluation.
#
# Physics: geometric optics in a cold unmagnetized plasma, refractive index
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
