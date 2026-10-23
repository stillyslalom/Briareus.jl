# Implosion metrics: the drive's mode 1 as a hot-spot velocity, the piston model's
# yield factors, P2/P4 swing metrics over a drive history, and the mode uncertainty
# carried to a velocity covariance, a velocity inverse and a velocity penalty.
# Sources and confidence: reference/research-notes.md, "Implosion metrics".

"""
    HOTSPOT_VELOCITY_PER_PERCENT = 100.0

Hot-spot velocity (km/s, burn-off) per percent of capsule flux mode 1 (`a₁/a₀`, the
norm of the `l = 1` vector, `mode_amplitude(c, 1)`), applied over the rise and peak.
An engineering default, not a facility criterion: Young et al., Phys. Plasmas 27,
082702 (2020) (N170601: 0.44 % → 37.6 km/s in 3D HYDRA, 44.7 km/s from view factors
and their Green's function), MacGowan et al., HEDP 40, 100944 (2021) (~100 km/s per
%), and Kritcher et al. 2024 (N230729: ~0.7 % → 70–75 km/s burn-off). It is platform
specific and linear only below about 0.75 % (Young; Spears et al. 2014 find 120 km/s at
2 %, saturated).
"""
const HOTSPOT_VELOCITY_PER_PERCENT = 100.0

"""
    mode1_vector(c) -> SVector{3}

The flux's mode-1 vector from real spherical-harmonic coefficients `c` (`mode_index`
order, e.g. `flux_ylm`): `√3 (c₁₁, c₁,₋₁, c₁₀)/c₀₀` along (x, y, z), pointing toward
the high-flux side. Its norm is `mode_amplitude(c, 1)`, the fractional amplitude `ε` of
a dipole drive `1 + ε cos γ` (not a percentage). Generic in the number type.
"""
function mode1_vector(c::AbstractVector)
    length(c) >= 4 || throw(ArgumentError("mode 1 needs the coefficients through l = 1"))
    s = sqrt(3) / c[1]
    return SVector(s * c[mode_index(1, 1)], s * c[mode_index(1, -1)], s * c[mode_index(1, 0)])
end

"""
    hotspot_velocity(c; per_percent = HOTSPOT_VELOCITY_PER_PERCENT, amplification = 1) -> SVector{3}
    hotspot_velocity(ev::Evaluation; ...)
    hotspot_velocity(history; kernel = nothing, window = nothing, per_percent, amplification)

Hot-spot velocity (km/s, chamber x, y, z) that the capsule flux's mode 1 drives: the
hot spot moves away from the high-flux side (Casey et al., PRL 126, 025002, 2021), at
`per_percent` km/s per percent of `mode1_vector` (an engineering default, see
`HOTSPOT_VELOCITY_PER_PERCENT`). `amplification` multiplies the burn-off velocity, e.g.
by `burn_amplification` to compare with a burn-on nToF measurement. Generic in the number type (usable in objectives).

From a `drive_history` result (run with `ylm = true`; without it only the z component,
from the Legendre moments), the velocity is the convolution of Young et al. (Phys.
Plasmas 27, 082702, 2020, Eq. 9), `v = −∫ G₁(t) · 100 p₁(t) dt`, with `G₁` in km/s per
% per ns: `kernel` (a `ResponseKernel` or any function of time). The default is
`flat_kernel(window...; total = per_percent)`, a flat weight over `window` (default:
the whole history) that gives `per_percent` for a constant 1 %; Young's published
kernels are negative in the foot and dominated by the rise and peak, so pass `window`
as the rise and peak, or a kernel of your own. The mode 1 is taken as linear between
the history's steps and integrated exactly against a tabulated kernel (two-point Gauss
on each interval between their breakpoints; a function kernel gets 16 pieces per
interval); the kernel's weight outside the history's span is dropped.
"""
function hotspot_velocity(c::AbstractVector; per_percent::Real = HOTSPOT_VELOCITY_PER_PERCENT,
                          amplification::Real = 1.0)
    return -(100 * per_percent * amplification) * mode1_vector(c)
end

hotspot_velocity(ev::Evaluation; kwargs...) = hotspot_velocity(flux_ylm(ev; lmax = 1); kwargs...)

function hotspot_velocity(h::NamedTuple; kernel = nothing, window = nothing,
                          per_percent::Real = HOTSPOT_VELOCITY_PER_PERCENT, amplification::Real = 1.0)
    t = h.time
    w = window === nothing ? (t[1], t[end]) : window
    K = kernel === nothing ? flat_kernel(w...; total = per_percent) : kernel
    p = _mode1_history(h)
    v = -100 * amplification .* _kernel_integral(t, p, K, (t[1], t[end]))
    return SVector(v[1], v[2], v[3])
end

# 3 × nt mode-1 vectors of a drive history (z only from Legendre moments)
function _mode1_history(h)
    if h.ylm !== nothing
        size(h.ylm, 1) >= 4 || throw(ArgumentError("history's ylm does not reach l = 1"))
        return reduce(hcat, [mode1_vector(view(h.ylm, :, k)) for k in axes(h.ylm, 2)])
    end
    size(h.moments, 1) >= 2 || throw(ArgumentError("history's moments do not reach l = 1"))
    z = h.moments[2, :] ./ h.moments[1, :]
    return vcat(zeros(2, length(z)), transpose(z))
end

"""
    ResponseKernel(times, values)

A response (Green's) function of time, piecewise linear through the table and zero
outside it; callable, `k(t)`. A time may repeat once to make a step. Units are the
observable per percent of drive mode per ns (km/s per % per ns for the hot-spot
velocity, Young et al. 2020; µm per % per ns for a hot-spot P2 kernel, Masse et al.
2019). Briareus ships no published kernel: those are platform specific and printed only
as figures.
"""
struct ResponseKernel
    t::Vector{Float64}
    g::Vector{Float64}
    function ResponseKernel(t, g)
        length(t) == length(g) >= 2 || throw(ArgumentError("a kernel table needs at least two times and as many values"))
        issorted(t) || throw(ArgumentError("kernel times must be sorted"))
        new(collect(Float64, t), collect(Float64, g))
    end
end

function (k::ResponseKernel)(t)
    (t < k.t[1] || t > k.t[end]) && return 0.0
    i = clamp(searchsortedlast(k.t, t), 1, length(k.t) - 1)
    a, b = k.t[i], k.t[i+1]
    b == a && return k.g[i+1]
    w = (t - a) / (b - a)
    return (1 - w) * k.g[i] + w * k.g[i+1]
end

"""
    flat_kernel(t0, t1; total = HOTSPOT_VELOCITY_PER_PERCENT) -> ResponseKernel

A constant kernel on `[t0, t1]` whose integral is `total`: a constant 1 % mode over the
window gives `total`. The placeholder for Young's mode-1 Green's function (*derived*:
it keeps the published static sensitivity and puts all weight on the window, which
should be the rise and peak).
"""
function flat_kernel(t0::Real, t1::Real; total::Real = HOTSPOT_VELOCITY_PER_PERCENT)
    t1 > t0 || throw(ArgumentError("the kernel window must have t1 > t0"))
    g = total / (t1 - t0)
    return ResponseKernel([t0, t1], [g, g])
end

_breakpoints(::Any) = Float64[]
_breakpoints(k::ResponseKernel) = k.t

# ∫ K(t) y(t) dt over `span`, y piecewise linear through (t, columns of Y) and K a kernel:
# two-point Gauss on each interval between the merged breakpoints (exact for a
# tabulated kernel; a function kernel gets 16 pieces per interval). Y is a vector or a matrix (one integral per row).
function _kernel_integral(t, Y, K, span)
    a, b = span
    pts = sort!(unique!(vcat(a, b, filter(x -> a < x < b, vcat(t, _breakpoints(K))))))
    y(τ) = Y isa AbstractVector ? _interp(t, Y, τ) : [_interp(t, view(Y, i, :), τ) for i in axes(Y, 1)]
    out = zero(y(a))
    g = 1 / (2sqrt(3))
    n = K isa ResponseKernel ? 1 : 16        # a smooth function: 16 pieces per interval
    for k in 1:length(pts)-1, j in 1:n
        lo, hi = pts[k] + (j - 1) * (pts[k+1] - pts[k]) / n, pts[k] + j * (pts[k+1] - pts[k]) / n
        hi > lo || continue
        m, h = (lo + hi) / 2, hi - lo
        for τ in (m - g * h, m + g * h)
            out = out .+ (h / 2 * K(τ)) .* y(τ)
        end
    end
    return out
end

"""
    yield_factor(f; a = 3.3) -> Real

No-α yield relative to 1D for a mode-1 (piston) asymmetry `f = v_HS/v_imp`:
`(1 − f²)ᵃ` (Hurricane et al., Phys. Plasmas 27, 062704, 2020, Eq. 20), with `a` the
temperature exponent of ⟨σv⟩ (3.3 for 4–7 keV; `reactivity_exponent`). Generic.
"""
yield_factor(f; a::Real = 3.3) = max(1 - f^2, zero(f))^a

"""
    alpha_factor(f, chi; a = 3.3, q = a - 1) -> Real

The α-heating part of the yield degradation, `Y_amp,3D/Y_amp,1D = exp(−q χ^1.2 f²)`,
with `χ` the 1D generalized Lawson parameter (`Y_amp ≈ exp(χ^1.2)`). `q = a − 1` is
Hurricane et al. 2020 (Eq. 21 and the tolerance Eq. 26, from which their Table I comes);
Casey et al., PRE 108, L053203 (2023), Eq. 5, has `q = 1.2(a − 1)`. Generic.
"""
alpha_factor(f, chi; a::Real = 3.3, q::Real = a - 1) = exp(-q * chi^1.2 * f^2)

"""
    rke_fraction(f; c = 0.66) -> Real

Residual kinetic energy over the shell's implosion kinetic energy, `RKE/KE = c f²`
(Hurricane et al. 2020, Eq. 23: `f = v_HS/v_imp = √(RKE/(c KE))`; c ≈ 0.66 in NIF-like
simulations, 1 for incompressible pistons).
"""
rke_fraction(f; c::Real = 0.66) = c * f^2

"""
    reactivity_exponent(T) -> Float64

The exponent `a` of `⟨σv⟩ ∝ Tᵃ` at ion temperature `T` (keV): `6 − 1.5 ln T` (Casey et
al., PRE 108, L053203, 2023). 3.3 at 6 keV, Hurricane's 4–7 keV value.
"""
reactivity_exponent(T::Real) = 6 - 1.5 * log(T)

"""
    burn_amplification(Yamp) -> Float64

Factor by which α heating raises the measured (burn-on) hot-spot velocity over the
burn-off one, as a function of the yield amplification `Y_amp`: `0.324 ln(Y_amp +
2.239) + 0.653` (Nora et al., Phys. Plasmas 30, 092705, 2023, Eq. 4, fitted for
`Y_amp < 40`; refused above). It is about 1.8 at `Y_amp = 30`, against burn-on/burn-off
ratios of 1.27–1.39 in Kritcher et al. 2024 (Table II): treat it as an upper estimate.
"""
function burn_amplification(Yamp::Real)
    1 <= Yamp <= 40 || throw(DomainError(Yamp, "Nora et al.'s fit covers 1 ≤ Y_amp ≤ 40"))
    return 0.324 * log(Yamp + 2.239) + 0.653
end

"""
    implosion_metrics(v; vimp = 380, a = 3.3, c = 0.66, chi = nothing, q = a - 1, ke = nothing)
    implosion_metrics(ev::Evaluation; per_percent, amplification, kwargs...)
    implosion_metrics(history; kernel, window, per_percent, amplification, kwargs...)

The piston model's metrics for a hot-spot velocity `v` (km/s, a 3-vector; or computed
with `hotspot_velocity` from an evaluation or a drive history). Returns `(; velocity,
speed, f, rke_fraction, rke, yield_factor, alpha_factor, total)`:
- `f = |v|/vimp`, the shell asymmetry `(ρR_max − ρR_min)/(ρR_max + ρR_min)` (Hurricane
  et al. 2020, Eq. 23); `vimp` defaults to Hurricane's Table I value, 380 km/s (an
  engineering default: use the design's implosion velocity);
- `rke_fraction = c f²` (`rke_fraction`), and `rke` in the units of `ke` when given;
- `yield_factor`, the no-α yield over 1D, `(1 − f²)ᵃ`;
- `alpha_factor`, with `chi` (the 1D generalized Lawson parameter) given, otherwise 1;
- `total = yield_factor × alpha_factor`.
Generic in the number type.
"""
function implosion_metrics(v::AbstractVector; vimp::Real = 380.0, a::Real = 3.3, c::Real = 0.66,
                           chi = nothing, q::Real = a - 1, ke = nothing)
    length(v) == 3 || throw(ArgumentError("a hot-spot velocity is a 3-vector (km/s)"))
    vimp > 0 || throw(ArgumentError("vimp must be positive"))
    velocity = SVector(v[1], v[2], v[3])
    speed = norm(velocity)
    f = speed / vimp
    Y = yield_factor(f; a)
    A = chi === nothing ? one(Y) : alpha_factor(f, chi; a, q)
    r = rke_fraction(f; c)
    return (; velocity, speed, f, rke_fraction = r, rke = ke === nothing ? nothing : r * ke,
            yield_factor = Y, alpha_factor = A, total = Y * A)
end

function implosion_metrics(ev::Evaluation; per_percent::Real = HOTSPOT_VELOCITY_PER_PERCENT,
                           amplification::Real = 1.0, kwargs...)
    return implosion_metrics(hotspot_velocity(ev; per_percent, amplification); kwargs...)
end

function implosion_metrics(h::NamedTuple; kernel = nothing, window = nothing,
                           per_percent::Real = HOTSPOT_VELOCITY_PER_PERCENT,
                           amplification::Real = 1.0, kwargs...)
    return implosion_metrics(hotspot_velocity(h; kernel, window, per_percent, amplification); kwargs...)
end

"""
    mode1_tolerance(loss; vimp = 380, a = 3.3, chi = 0.0, q = a - 1) -> NamedTuple

The largest mode-1 asymmetry that keeps the yield loss `1 − (1 − f²)ᵃ exp(−q χ^1.2 f²)`
(Hurricane et al. 2020, Eq. 26) at `loss`: `(; f, velocity)`, the shell asymmetry and
the hot-spot speed (km/s) at `vimp`. The inverse of `implosion_metrics`'s `total`.
"""
function mode1_tolerance(loss::Real; vimp::Real = 380.0, a::Real = 3.3, chi::Real = 0.0,
                         q::Real = a - 1)
    0 < loss < 1 || throw(ArgumentError("loss is a fraction in (0, 1)"))
    kept(f) = yield_factor(f; a) * alpha_factor(f, chi; a, q)
    lo, hi = 0.0, 1.0
    for _ in 1:60
        mid = (lo + hi) / 2
        kept(mid) > 1 - loss ? (lo = mid) : (hi = mid)
    end
    f = (lo + hi) / 2
    return (; f, velocity = f * vimp)
end

"""
    swing_metrics(history; l = 2, window = nothing, picket = nothing, limit = 0.075,
                  kernels = (;)) -> NamedTuple

Time metrics of the capsule flux's `aₗ/a₀` over a `drive_history` result (taken as
linear between its steps). Returns `(; mean, max, min, swing, deviation, picket_mean,
picket_max, picket_ok, response)`:
- over `window` (default: the whole history): the time average `mean`, the extremes and
  `swing = max − min`, and `deviation`, the time average of `|aₗ/a₀ − mean|`;
- over `picket` (default: the first 2 ns of the history), the time average and the
  largest `|aₗ/a₀|`; `picket_ok` compares the average with `limit`, Milovich et al.'s
  (Phys. Plasmas 23, 032701, 2016) 7.5 % early-time P2/P0 for under 10 % margin loss
  (whether their limit is on the average or the peak is not stated: *derived*);
- `response`: for each `name => kernel` in `kernels` (`ResponseKernel`s or functions,
  per % per ns), `∫ K(t) · 100 aₗ/a₀ dt` in the kernel's units, e.g. a hot-spot or ρR
  P2 Green's function (Masse et al., Phys. Plasmas 26, 2019). Without kernels it is an
  empty NamedTuple: the flux history alone does not give a hot-spot shape. Masse's
  hot-spot kernel changes sign before about 4.5 ns while the ρR one does not, so a
  time-integrated flux P2 says neither, and Briareus reports no µm by default.
All fractions (not percent).
"""
function swing_metrics(h::NamedTuple; l::Integer = 2, window = nothing, picket = nothing,
                       limit::Real = 0.075, kernels = (;))
    t = h.time
    size(h.moments, 1) > l || throw(ArgumentError("the history's moments do not reach l = $l"))
    y = h.moments[l+1, :] ./ h.moments[1, :]
    w = window === nothing ? (t[1], t[end]) : window
    p = picket === nothing ? (t[1], min(t[1] + 2.0, t[end])) : picket
    tw, yw = _clip_series(t, y, w)
    T = tw[end] - tw[1]
    mean = _trapz(tw, yw) / T
    deviation = sum(_abs_integral(tw[k+1] - tw[k], yw[k] - mean, yw[k+1] - mean)
                    for k in 1:length(tw)-1) / T
    tp, yp = _clip_series(t, y, p)
    picket_mean = _trapz(tp, yp) / (tp[end] - tp[1])
    response = map(K -> 100 * _kernel_integral(t, y, K, (t[1], t[end])), kernels)
    return (; mean, max = maximum(yw), min = minimum(yw), swing = maximum(yw) - minimum(yw),
            deviation, picket_mean, picket_max = maximum(abs, yp),
            picket_ok = abs(picket_mean) <= limit, response)
end

# the piecewise-linear series restricted to [a, b], with interpolated end values
function _clip_series(t, y, (a, b))
    t[1] <= a < b <= t[end] || throw(ArgumentError("window ($a, $b) must lie within the history ($(t[1]), $(t[end]))"))
    ts = vcat(a, [τ for τ in t if a < τ < b], b)
    return ts, [_interp(t, y, τ) for τ in ts]
end

_trapz(t, y) = sum((t[k+1] - t[k]) * (y[k] + y[k+1]) / 2 for k in 1:length(t)-1)

# ∫|y| over an interval of length h on which y goes linearly from u to v
function _abs_integral(h, u, v)
    (u >= 0) == (v >= 0) && return h * abs(u + v) / 2
    return h * (u^2 + v^2) / (2 * (abs(u) + abs(v)))
end

"""
    velocity_uncertainty(scenario, beams; per_percent, amplification = 1, vimp = 380,
                         a = 3.3, chi = nothing, q = a - 1, kwargs...) -> NamedTuple
    velocity_uncertainty(u; ...)

The beam errors of `mode_uncertainty` (its keywords: `pointing`, `power`,
`group_pointing`, `group_power`, `groups`, `positioning`), or its result `u`, carried to
the hot-spot velocity at `per_percent` km/s per % of mode 1 (`hotspot_velocity`).
Returns `(; velocity, covariance, sigma, rms, rms_by, speed_rms, yield_factor,
expected_yield, total, expected_total)`:
- `velocity`, the nominal (km/s), `covariance` (3 × 3, km²/s²) and `sigma` its per-axis
  rms; `rms = √tr Σ`, split by source in `rms_by` (`beam`, `group`, `positioning`; they
  add in quadrature) and `speed_rms = √E|v|²`;
- `yield_factor`, the no-α factor at the nominal velocity, and `expected_yield`, its
  expectation over the errors, from the Gaussian closed form of `E[exp(−a f²)]`
  (`(1 − f²)ᵃ ≈ exp(−a f²)`, within 1 % at f = 0.25); `total` and `expected_total`
  include the α factor (`chi` given) the same way, `exp(−(a + q χ^1.2) f²)`.
Linear in the errors, like `mode_uncertainty`.
"""
function velocity_uncertainty(s::Scenario, beams::AbstractVector{<:Beam};
                              per_percent::Real = HOTSPOT_VELOCITY_PER_PERCENT,
                              amplification::Real = 1.0, vimp::Real = 380.0, a::Real = 3.3,
                              chi = nothing, q::Real = a - 1, kwargs...)
    u = mode_uncertainty(s, beams; lmax = 1, kwargs...)
    return velocity_uncertainty(u; per_percent, amplification, vimp, a, chi, q)
end

function velocity_uncertainty(u::NamedTuple; per_percent::Real = HOTSPOT_VELOCITY_PER_PERCENT,
                              amplification::Real = 1.0, vimp::Real = 380.0, a::Real = 3.3,
                              chi = nothing, q::Real = a - 1)
    M = _velocity_map(per_percent * amplification)
    velocity = SVector{3}(M * u.c[2:4])
    C = M * u.covariance[2:4, 2:4] * transpose(M)
    C = (C + transpose(C)) / 2
    K = 100 * per_percent * amplification * sqrt(3)
    rms_by = map(r -> K * r[1], u.rms_by)
    f = norm(velocity) / vimp
    β = chi === nothing ? Float64(a) : a + q * chi^1.2
    E(b) = _expected_exp(b / vimp^2, velocity, C)
    Y = yield_factor(f; a)
    return (; velocity, covariance = C, sigma = sqrt.(max.(diag(C), 0.0)),
            rms = sqrt(max(tr(C), 0.0)), rms_by, speed_rms = sqrt(sum(abs2, velocity) + max(tr(C), 0.0)),
            yield_factor = Y, expected_yield = E(a),
            total = chi === nothing ? Y : Y * alpha_factor(f, chi; a, q), expected_total = E(β))
end

# v = M (c₁,₋₁, c₁₀, c₁₁)/c₀₀ for a sensitivity k (km/s per %)
function _velocity_map(k)
    K = 100 * k * sqrt(3)
    return [0.0 0.0 -K; -K 0.0 0.0; 0.0 -K 0.0]
end

# E[exp(−b |v|²)] for v ~ N(μ, Σ): det(I + 2bΣ)^(−1/2) exp(−b μᵀ (I + 2bΣ)⁻¹ μ)
function _expected_exp(b, μ, Σ)
    A = I + 2b * Σ
    return exp(-b * dot(μ, A \ Vector(μ))) / sqrt(det(A))
end

"""
    ntof_covariance(v; speed = 14, theta = 13, phi = 29) -> Matrix

Covariance (km²/s²) of an nToF hot-spot velocity measurement `v` (km/s): rms `speed` in
magnitude and `theta`, `phi` degrees in direction (NIF: ±14 km/s, 13°, 29°; Young et al.
2020), along the local radial, polar and azimuthal directions of `v`.
"""
function ntof_covariance(v; speed::Real = 14.0, theta::Real = 13.0, phi::Real = 29.0)
    u = SVector{3,Float64}(v...)
    r = norm(u)
    r > 0 || return Matrix(Diagonal(fill(Float64(speed)^2, 3)))
    er = u / r
    θ = acos(clamp(er[3], -1.0, 1.0))
    φ = atan(er[2], er[1])
    eθ = SVector(cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ))
    eφ = SVector(-sin(φ), cos(φ), 0.0)
    σ = (speed, r * deg2rad(theta), r * sin(θ) * deg2rad(phi))
    return sum(σ[i]^2 * e * transpose(e) for (i, e) in enumerate((er, eθ, eφ)))
end

"""
    infer_errors(scenario, beams, v::AbstractVector{<:Real}; sigma = ntof_covariance(v),
                 per_percent, amplification = 1, thickness = 0.002,
                 thickness_per_percent = 110, kwargs...) -> NamedTuple

The most probable beam errors behind a measured hot-spot velocity `v` (km/s, chamber
x, y, z): `infer_errors` on the mode-1 components, with the measurement mapped to them
through `hotspot_velocity` (`per_percent`, `amplification`; divide out a burn
amplification here) and its covariance `sigma` (3 × 3, km²/s², or one rms number; default
`ntof_covariance(v)`). The priors are `infer_errors`'s keywords.

A capsule-thickness mode 1 competes with the drive (Casey et al., PRL 126, 025002,
2021): the hot spot moves toward the thicker side, at `thickness_per_percent` km/s per %
of the shell thickness (110 simulated, 140 ± 30 measured), so the opposite sign of a
flux mode 1 on that side. It enters as a nuisance error, a vector `t` (fractions of the
thickness, toward the thicker side) with rms `thickness` per axis (default 0.2 %, about
0.15 µm metrology on a 78 µm shell: *derived*; 0 leaves it out), and is inferred with
the beam errors.

Returns `infer_errors`'s fields, with `measured`, `predicted` and `residual` still in
the mode-1 components, plus `thickness` (the inferred `t`) and `velocity = (; measured,
nominal, predicted, residual, chi2, chi2_prior)` in km/s: `predicted` is the nominal plus
the inferred errors' and thickness's velocity; `chi2` is MacGowan's χ²_ν of the nominal
(the forward prediction) against the measurement, `(v − v₀)ᵀ Σₘ⁻¹ (v − v₀)/3` (Casey et
al. 2023), and `chi2_prior` the same with the errors' prior spread added to `Σₘ`
(about 1 when the measurement is what the priors expect).
"""
function infer_errors(s::Scenario, beams::AbstractVector{<:Beam}, v::AbstractVector{<:Real};
                      sigma = ntof_covariance(v), per_percent::Real = HOTSPOT_VELOCITY_PER_PERCENT,
                      amplification::Real = 1.0, thickness::Real = 0.002,
                      thickness_per_percent::Real = 110.0, kwargs...)
    length(v) == 3 || throw(ArgumentError("a hot-spot velocity is a 3-vector (km/s)"))
    thickness >= 0 || throw(ArgumentError("thickness must be a non-negative rms"))
    Σv = sigma isa Real ? Matrix(Diagonal(fill(Float64(sigma)^2, 3))) : Matrix{Float64}(sigma)
    size(Σv) == (3, 3) || throw(ArgumentError("sigma is one rms or a 3 × 3 covariance"))
    M = _velocity_map(per_percent * amplification)
    Mi = inv(M)
    eff = apply_lpi(s.lpi, beams, s.time)
    c = _ylm_product(s.ylm_response, deposit(s.deposition, eff).element, 4)
    c = c ./ c[1]
    v0 = M * c[2:4]
    Δc = Mi * (Vector{Float64}(v) - v0)
    Σc = Mi * Σv * transpose(Mi)
    # thickness: v = 100 k_t A t (toward the thicker side), as mode-1 columns
    nuisance = thickness > 0 ?
        (; columns = Mi * (100 * thickness_per_percent * amplification * Matrix(1.0I, 3, 3)),
           variances = fill(Float64(thickness)^2, 3), names = ["thickness.x", "thickness.y", "thickness.z"],
           label = "capsule thickness") : nothing
    comps = [(1, -1), (1, 0), (1, 1)]
    r = _infer_errors(s, beams, comps, Δc, (Σc + transpose(Σc)) / 2; nuisance, kwargs...)
    t = thickness > 0 ? SVector{3}(r.x[end-2:end]) : zero(Vec3)
    predicted = v0 + M * r.predicted
    # the prior predictive spread of the velocity: the errors' (and thickness's) plus Σₘ
    Sv = M * (r.J * (r.prior_sigma .^ 2 .* transpose(r.J))) * transpose(M) + Σv
    d = Vector{Float64}(v) - v0
    velocity = (; measured = SVector{3,Float64}(v...), nominal = SVector{3}(v0),
                predicted = SVector{3}(predicted), residual = SVector{3}(Vector{Float64}(v) - predicted),
                chi2 = dot(d, Symmetric(Σv) \ d) / 3, chi2_prior = dot(d, Symmetric(Sv) \ d) / 3)
    return (; r..., thickness = t, velocity)
end

"""
    VelocityPenalty(; scale = 10, weight = 1, per_percent, amplification = 1,
                    uncertainty = nothing)

An objective term on the hot-spot velocity of an evaluation's capsule flux
(`hotspot_velocity`): `weight × |v|²/scale²` (`scale` in km/s), the mode 1 weighted by
its velocity sensitivity rather than equally with the other modes. Callable on an
`Evaluation`; add it to an objective, e.g.
`obj = let J = asymmetry_objective(; modes = (2 => 1.0, 4 => 1.0)), P = VelocityPenalty(); ev -> J(ev) + P(ev) end`.
With `uncertainty` (`mode_uncertainty`'s error keywords, as for the objectives) it adds
the errors' expected `|δv|²`. Generic in the number type (ForwardDiff).
"""
struct VelocityPenalty{E}
    scale::Float64
    weight::Float64
    per_percent::Float64
    amplification::Float64
    errors::E
end

function VelocityPenalty(; scale::Real = 10.0, weight::Real = 1.0,
                         per_percent::Real = HOTSPOT_VELOCITY_PER_PERCENT,
                         amplification::Real = 1.0, uncertainty = nothing)
    scale > 0 || throw(ArgumentError("scale must be positive"))
    return VelocityPenalty(Float64(scale), Float64(weight), Float64(per_percent),
                           Float64(amplification), _uncertainty_spec(uncertainty))
end

function (p::VelocityPenalty)(ev::Evaluation)
    v = hotspot_velocity(ev; per_percent = p.per_percent, amplification = p.amplification)
    J = sum(abs2, v)
    if p.errors !== nothing
        var = _error_variance(ev.scenario, ev.effective, ev.scenario.ylm_response, 4, p.errors)
        J += (100 * p.per_percent * p.amplification)^2 * 3 * (var[2] + var[3] + var[4])
    end
    return p.weight * J / p.scale^2
end
