# Reduced cross-beam energy transfer (CBET) at the laser entrance holes.
#
# The coupling is the kinetic ray-pair gain of Nguyen et al., arXiv:2304.13268, Eqs. 1–5,
# which is identical to Michel et al. 2009's γ = (Δk²/8k₀)Γ. In practical units it is
#
#     g_bc = −C λ_µm sin²(θ/2) Im Γ(Ω − k·u, k) · P_pol,   C = 4.59×10⁻¹⁸ µm⁻¹ (W/cm²)⁻¹,
#
# the spatial gain rate of beam b per unit intensity of beam c, with
# Γ = χe(1+χi)/(1+χe+χi), Maxwellian χ_s = (kλ_Ds)⁻²[1 + ξ_s Z(ξ_s)],
# ξ_s = (Ω − k·u)/(√2 k v_Ts), Ω = ω_b − ω_c, k = k_b − k_c. Power flows from the beam
# of higher frequency in the flow frame to the lower one.

const _C_LIGHT = 2.99792458e8      # m/s
const _QE = 1.602176634e-19        # C
const _ME = 9.1093837015e-31       # kg
const _MP = 1.67262192369e-27      # kg
const _EPS0 = 8.8541878128e-12     # F/m

"""
    plasma_dispersion(ξ) -> Complex

Fried–Conte plasma dispersion function for real argument, `Z(ξ) = −2 D(ξ) + i√π exp(−ξ²)`,
where D is Dawson's integral. The function is differentiable.
"""
plasma_dispersion(ξ::Real) = complex(-2 * dawson(ξ), sqrt(π) * exp(-ξ^2))

"""
    cbet_coupling(Ω, k, kdotu; ne, Te, Ti, ions) -> Γ

Ion-acoustic response `Γ = χe(1+χi)/(1+χe+χi)` for beat frequency `Ω` (rad/s), beat
wavenumber `k` (1/m), Doppler shift `k·u` (rad/s), electron density `ne` (m⁻³),
temperatures `Te`, `Ti` (eV), and ion species `ions` — `(fraction, Z, A)` tuples with
number fractions summing to one. `χi` sums over species, each with its Landau damping;
the damping matters for CH, where the protons dominate.
"""
function cbet_coupling(Ω, k, kdotu; ne, Te, Ti, ions)
    Zbar = sum(f * Z for (f, Z, _) in ions)
    ni = ne / Zbar
    ω = Ω - kdotu
    λDe2 = _EPS0 * Te * _QE / (ne * _QE^2)
    ξe = ω / (sqrt(2) * k * sqrt(Te * _QE / _ME))
    χe = (1 + ξe * plasma_dispersion(ξe)) / (k^2 * λDe2)
    χi = zero(χe)
    for (f, Z, A) in ions
        λD2 = _EPS0 * Ti * _QE / (f * ni * Z^2 * _QE^2)
        ξ = ω / (sqrt(2) * k * sqrt(Ti * _QE / (A * _MP)))
        χi += (1 + ξ * plasma_dispersion(ξ)) / (k^2 * λD2)
    end
    return χe * (1 + χi) / (1 + χe + χi)
end

"Fully ionized CH (1:1): `[(0.5, 6, 12.011), (0.5, 1, 1.008)]`."
const CH_IONS = [(0.5, 6.0, 12.011), (0.5, 1.0, 1.008)]

"""
    CBET(; leh_z, ne_nc = 0.06, Te = 5000, Ti = 1500, ions = CH_IONS, mach = 1.0,
         flow = nothing, length = 1.5, polarization = :smoothed, multiplier = 1.0,
         nsteps = 20)

Reduced steady-state CBET among the beams that share each laser entrance hole, applied
as a laser–plasma model (`apply_lpi`, `Scenario(...; lpi)`).

- `leh_z` — z (mm) of the LEH planes; beams whose ports have `z > 0` cross at `+leh_z`,
  the rest at `-leh_z`.
- Plasma in the crossing region: `ne_nc` (electron density over the critical density at
  the nominal wavelength), `Te`, `Ti` (eV), and ion species `ions` (`(fraction, Z, A)`
  number fractions; default fully ionized CH). Defaults are the NIF LEH conditions of
  Michel et al. 2009; `Ti` is not given in that paper and is an assumption.
- `flow` — outflow speed along the hohlraum axis (km/s, away from the hohlraum), or `mach`
  times the ion sound speed `√((Z̄ Te + 3 Ti)/(Ā m_p))` when `flow` is nothing. The value
  is signed: a negative `mach` or `flow` points into the hohlraum. With outflow the inner
  (steeper) cones get the larger blue shift, so flow alone moves power from the inner to
  the outer cones; red-detuned inner cones (`dlambda > 0`) take power back.
  [`TwoZoneCBET`](@ref) puts a supersonic zone outside the LEH and a slow one inside.
- `length` — maximum interaction length (mm). A pair crossing at angle θ interacts over
  `min(length, w/sin θ)`, with `w` the mean beam width at the LEH plane.
- `polarization` — `:smoothed` for NIF polarization smoothing with uncorrelated
  speckle, which multiplies the gain by (1+cos²θ)/4, or `:parallel` for co-polarized
  beams, which leave the gain unchanged.
- `multiplier` — calibration factor on all gains. Fit it to shot data: published NIF
  designs tune CBET with saturation clamps and power multipliers.
- `clamp` — ion-wave saturation level δn/n (`nothing` = linear theory, which grossly
  overpredicts transfer at NIF intensities). The coupling of each pair is multiplied by
  `1/√(1 + (δn/clamp)²)`, a smooth version of capping the linear δn/n at `clamp`.
  Published NIF calibrations range from 3×10⁻⁴ (2011) to 0.002 (foot) / 0.0075 (peak)
  for N221204 (Kritcher et al., PRE 109, 025204, 2024); the default is 0.002.
  [`cbet_clamp_schedule`](@ref) builds the foot/peak schedule from a `PulseShape`.
- `saturation` — how `clamp` acts: `:smooth` (default, the factor above) or `:min`,
  the hard limit `δn = min(clamp, δn_linear)` of Strozzi et al., PRL 118, 025002 (2017),
  i.e. a factor `min(1, clamp/δn)`, as in HYDRA's inline model, so its clamp values carry
  over one to one. At δn = clamp the smooth factor is √2 lower.

Any plasma parameter and `multiplier` may be a number, a function of time (ns), or a
`(times, values)` table. The exchange of each pair uses the geometric-mean beam area at
the LEH, so total power is conserved. The coupled exchange is integrated with `nsteps` RK4
steps. Beam detuning comes from `Beam.dlambda` (Å at the beam wavelength — NIF
quotes Δλ at 1ω, so divide those values by 3 for 3ω beams).

The pair coefficients and the exchanges at the two LEHs are computed on threads; the
results are identical for every thread count.
"""
Base.@kwdef struct CBET <: LaserPlasmaModel
    leh_z::Any
    ne_nc::Any = 0.06
    Te::Any = 5000.0
    Ti::Any = 1500.0
    ions::Vector{NTuple{3,Float64}} = CH_IONS
    mach::Any = 1.0
    flow::Any = nothing
    length::Float64 = 1.5
    polarization::Symbol = :smoothed
    multiplier::Any = 1.0
    clamp::Any = 0.002
    saturation::Symbol = :smooth
    nsteps::Int = 20
end

"Plasma parameters of `m` at time `t` as a NamedTuple (SI densities, eV, m/s flow)."
cbet_plasma(m::CBET, λ0_µm, t) = _cbet_plasma(_cbet_state(m, t), λ0_µm)

# Returns the state of the model at time `t`: the schedules of the plasma (`ne_nc`, `Te`,
# `Ti`, the flow) and the `multiplier` evaluated once, with the sound speed, the flow
# speed `u` (m/s) and the fields that the pair coefficients use. The schedules of a `CBET`
# are untyped (numbers, functions, tables), so this evaluation dispatches at run time.
# The result is a NamedTuple of concrete types, and the pair loop, which takes it through
# a function barrier, is inferred for it.
function _cbet_state(m::CBET, t)
    ne_nc, Te, Ti = _tvalue(m.ne_nc, t), _tvalue(m.Te, t), _tvalue(m.Ti, t)
    Z = sum(f * Z for (f, Z, _) in m.ions)
    A = sum(f * A for (f, _, A) in m.ions)
    cs = sqrt((Z * Te + 3Ti) * _QE / (A * _MP))
    u = m.flow === nothing ? _tvalue(m.mach, t) * cs : 1e3 * _tvalue(m.flow, t)
    return (; ne_nc, Te, Ti, ions = m.ions, u, cs, multiplier = _tvalue(m.multiplier, t),
            polarization = m.polarization)
end

function _cbet_plasma(st::NamedTuple, λ0_µm)
    ω0 = 2π * _C_LIGHT / (λ0_µm * 1e-6)
    nc = _EPS0 * _ME * ω0^2 / _QE^2
    return (; ne = st.ne_nc * nc, nc, st.Te, st.Ti, st.ions, st.u, st.cs)
end

"""
    cbet_gain(m, b, c, t) -> g

Linear (unsaturated) spatial gain rate (1/mm per TW/mm² of beam `c`) of beam `b` due to
beam `c` crossing at the LEH, at time `t`, including the polarization factor and
`multiplier`. The gain is antisymmetric: `cbet_gain(m, c, b, t) == -cbet_gain(m, b, c, t)`.
"""
cbet_gain(m::CBET, b::Beam, c::Beam, t = 0.0) = _cbet_pair(m, b, c, t)[1]

# Linear gain g and ion-wave amplitude coefficient D for a pair, with δn/n = D √(I_b I_c)
# for intensities in TW/mm². The amplitude is δn/n = (k²c²/2ω_pe²) |Γ| |a_b||a_c| with
# |a|² = 7.31e-19 I λ² for linear polarization, I in W/cm² and λ in µm. Polarization
# smoothing reduces the rms beat amplitude by √pol.
_cbet_pair(m::CBET, b::Beam, c::Beam, t) = _cbet_pair(_cbet_state(m, t), b, c)

function _cbet_pair(st::NamedTuple, b::Beam, c::Beam)
    pl = _cbet_plasma(st, b.wavelength)
    db, dc = direction(b), direction(c)
    cosθ = clamp(dot(db, dc), -1.0, 1.0)
    sin2h = (1 - cosθ) / 2                       # sin²(θ/2)
    z = zero(b.power * c.power * b.dlambda * c.dlambda)
    sin2h < 1e-12 && return (z, z)
    λb, λc = wavelength(b), wavelength(c)         # µm
    ωb, ωc = 2π * _C_LIGHT / (λb * 1e-6), 2π * _C_LIGHT / (λc * 1e-6)
    kb = ωb / _C_LIGHT * sqrt(1 - pl.ne / (_EPS0 * _ME * ωb^2 / _QE^2))
    kc = ωc / _C_LIGHT * sqrt(1 - pl.ne / (_EPS0 * _ME * ωc^2 / _QE^2))
    kvec = kb * db - kc * dc
    kmag = norm(kvec)
    # The outflow is along the hohlraum axis, away from the hohlraum at this LEH.
    uvec = Vec3(0, 0, b.port[3] > 0 ? 1.0 : -1.0) * pl.u
    Γ = cbet_coupling(ωb - ωc, kmag, dot(kvec, uvec); pl.ne, pl.Te, pl.Ti, pl.ions)
    pol = st.polarization === :smoothed ? (1 + cosθ^2) / 4 :
          st.polarization === :parallel ? 1.0 :
          throw(ArgumentError("polarization must be :smoothed or :parallel"))
    # 4.59e-18 µm⁻¹ per W/cm²  →  ×1e3 (per mm) ×1e14 (per TW/mm²)
    g = -0.459 * b.wavelength * sin2h * imag(Γ) * pol * st.multiplier
    ωpe2 = pl.ne * _QE^2 / (_EPS0 * _ME)
    D = kmag^2 * _C_LIGHT^2 / (2ωpe2) * abs(Γ) * 7.31e-19 * 1e14 * b.wavelength^2 * sqrt(pol)
    return (g, D)
end

# Width (mm) and effective area (mm², power / peak intensity) of a beam where its axis
# crosses the plane z = zL, using the moment-matched defocused core profile.
function _beam_at_plane(b::Beam, zL)
    k = direction(b)
    τ = (zL - b.aim[3]) / k[3]
    s = τ - b.defocus
    sp = b.spot
    c = (1 / (2b.fnumber))^2 / (4 * _sg_moment(sp.order))
    a = sqrt(sp.a^2 + c * s^2)
    bb = sqrt(sp.b^2 + c * s^2)
    return 2 * sqrt(a * bb), _sg_norm(sp.order) * a * bb
end

apply_lpi(m::CBET, beams::AbstractVector{<:Beam}, time = 0.0) =
    _apply_cbet(m, beams, time, Threads.nthreads() > 1)

# `threaded = false` is the serial reference. Both paths give identical results
# for any thread count: each pair is written to separate slots of K and D, the two LEHs
# share no data, and no sum is taken in a thread-dependent order.
function _apply_cbet(m::CBET, beams::AbstractVector{<:Beam}, time, threaded::Bool)
    T = promote_type(mapreduce(b -> typeof(b.power), promote_type, beams),
                     mapreduce(b -> typeof(b.dlambda), promote_type, beams),
                     mapreduce(b -> eltype(b.aim), promote_type, beams))
    # Model parameters may be dual numbers in a calibration.
    for f in (:ne_nc, :Te, :Ti, :mach, :flow, :multiplier, :clamp)
        v = getfield(m, f)
        v === nothing || (T = promote_type(T, typeof(_tvalue(v, time))))
    end
    # every pair coefficient is built from T-typed inputs and Float64 constants
    Tk = promote_type(T, Float64)
    P = T[b.power for b in beams]
    zL = _tvalue(m.leh_z, time)
    clampv = m.clamp === nothing ? nothing : _tvalue(m.clamp, time)
    m.saturation in (:smooth, :min) ||
        throw(ArgumentError("saturation must be :smooth or :min (got $(m.saturation))"))
    sat = clampv === nothing ? nothing : (clampv, m.saturation === :min)
    idxs = Vector{Int}[]
    for top in (true, false)
        idx = [i for (i, b) in enumerate(beams) if (b.port[3] > 0) == top]
        length(idx) < 2 || push!(idxs, idx)
    end
    isempty(idxs) && return [remake(b; power = P[i]) for (i, b) in enumerate(beams)]
    geos = [[_beam_at_plane(beams[i], beams[i].port[3] > 0 ? zL : -zL) for i in idx]
            for idx in idxs]
    # The exchange coefficients K[p,q] = g_pq ℓ_pq / Ā_pq (antisymmetric) and the ion-wave
    # amplitude coefficients D[p,q] (symmetric) form one pair of matrices per LEH.
    Ks = [zeros(Tk, length(idx), length(idx)) for idx in idxs]
    Ds = [zeros(Tk, length(idx), length(idx)) for idx in idxs]
    tasks = [(l, p, q) for l in eachindex(idxs) for p in eachindex(idxs[l])
             for q in p+1:length(idxs[l])]
    _cbet_couplings!(Ks, Ds, _cbet_state(m, time), m.length, beams, idxs, geos, tasks, threaded)
    Pls = [P[idx] for idx in idxs]
    Aeffs = [[g[2] for g in geo] for geo in geos]
    _cbet_exchanges!(Pls, Ks, Ds, Aeffs, sat, m.nsteps, threaded)
    for (idx, Pl) in zip(idxs, Pls)
        P[idx] = Pl
    end
    return [remake(b; power = P[i]) for (i, b) in enumerate(beams)]
end

# This function is a barrier: the coefficient eltype is a type parameter of this method,
# so the threaded closure does not box it, and the state `st` of the model at this time
# (`_cbet_state`) is concretely typed in it.
function _cbet_couplings!(Ks::Vector{Matrix{Tk}}, Ds::Vector{Matrix{Tk}}, st::NamedTuple, len::Float64, beams,
                          idxs, geos, tasks, threaded::Bool) where {Tk}
    function pair!(j)
        l, p, q = tasks[j]
        idx, geo, K, D = idxs[l], geos[l], Ks[l], Ds[l]
        b, c = beams[idx[p]], beams[idx[q]]
        cosθ = clamp(dot(direction(b), direction(c)), -1.0, 1.0)
        sinθ = sqrt(max(1 - cosθ^2, 1e-12))
        ℓ = min(len, (geo[p][1] + geo[q][1]) / 2 / sinθ)
        Ā = sqrt(geo[p][2] * geo[q][2])
        g, d = _cbet_pair(st, b, c)
        kpq = g * ℓ / Ā
        K[p, q] = kpq
        K[q, p] = -kpq
        D[p, q] = D[q, p] = d
        return nothing
    end
    if threaded && length(tasks) > 1
        Threads.@threads for j in eachindex(tasks)
            pair!(j)
        end
    else
        foreach(pair!, eachindex(tasks))
    end
    return nothing
end

function _cbet_exchanges!(Pls, Ks, Ds, Aeffs, clampv, nsteps, threaded::Bool)
    if threaded && length(Pls) > 1
        Threads.@threads for l in eachindex(Pls)
            _cbet_exchange!(Pls[l], Ks[l], Ds[l], Aeffs[l], clampv, nsteps)
        end
    else
        for l in eachindex(Pls)
            _cbet_exchange!(Pls[l], Ks[l], Ds[l], Aeffs[l], clampv, nsteps)
        end
    end
    return nothing
end

# RK4 integration of dP_p/dτ = P_p Σ_q K_pq S_pq P_q over τ ∈ [0, 1] in place. The
# saturation factor is S_pq = 1/√(1 + (δn_pq/clamp)²) when smooth and min(1, clamp/δn_pq)
# when hard. It is symmetric in p and q, so power is conserved. `clampv` is `nothing` or
# (clamp, hard). Work buffers are allocated once per call.
function _cbet_exchange!(x::AbstractVector, K, D, Aeff, clampv, nsteps)
    k1, k2, k3, k4, y = similar(x), similar(x), similar(x), similar(x), similar(x)
    I = similar(x, promote_type(eltype(x), eltype(Aeff)))
    h = 1 / nsteps
    for _ in 1:nsteps
        _cbet_rate!(k1, x, K, D, Aeff, clampv, I)
        y .= x .+ (h / 2) .* k1
        _cbet_rate!(k2, y, K, D, Aeff, clampv, I)
        y .= x .+ (h / 2) .* k2
        _cbet_rate!(k3, y, K, D, Aeff, clampv, I)
        y .= x .+ h .* k3
        _cbet_rate!(k4, y, K, D, Aeff, clampv, I)
        x .= x .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end
    return x
end

function _cbet_rate!(out, x, K, D, Aeff, ::Nothing, I)
    mul!(out, K, x)
    out .= x .* out
    return out
end

function _cbet_rate!(out, x, K, D, Aeff, sat::Tuple, I)
    clampv, hard = sat
    n = length(x)
    I .= x ./ Aeff
    for p in 1:n
        acc = zero(eltype(x))
        for q in 1:n
            q == p && continue
            δn = D[p, q] * sqrt(max(I[p] * I[q], zero(I[p])))
            s = hard ? (δn <= clampv ? one(δn) : clampv / δn) : 1 / sqrt(1 + (δn / clampv)^2)
            acc += K[p, q] * x[q] * s
        end
        out[p] = x[p] * acc
    end
    return out
end

"""
    ion_wave_amplitude(m::CBET, b, c, t = 0.0) -> δn/n

Linear ion-acoustic density perturbation driven by beams `b` and `c` at their LEH
plane intensities, for comparison with `m.clamp`.
"""
function ion_wave_amplitude(m::CBET, b::Beam, c::Beam, t = 0.0)
    zL = _tvalue(m.leh_z, t) * (b.port[3] > 0 ? 1 : -1)
    _, d = _cbet_pair(m, b, c, t)
    Ib = b.power / _beam_at_plane(b, zL)[2]
    Ic = c.power / _beam_at_plane(c, zL)[2]
    return d * sqrt(Ib * Ic)
end

"""
    TwoZoneCBET(outside::CBET, inside::CBET)
    TwoZoneCBET(; leh_z, offset = 0.5, outside = (;), inside = (;))

CBET in two crossing zones per entrance hole, applied in the order in which the beams meet
them. The `outside` zone lies outside the LEH, in the supersonic outflow, where the
Doppler shift of the flow moves power from the inner cones to the outer ones even without
detuning. The `inside` zone lies inside the LEH, in slow flow, where the transfer is set
by the wavelength detuning; red-detuned inner cones take power from the outer cones. Each
zone is a full [`CBET`](@ref) with a separate crossing plane, flow (signed `mach` or
`flow`), plasma state (`ne_nc`, `Te`, `Ti`, `ions`, as numbers, functions of time or
`(times, values)` tables), `clamp`, `multiplier` and `length`, and stays generic in the
number type.

The keyword form places the zones `offset` (mm) outside and inside the LEH planes
`±leh_z` and passes `outside`/`inside` (NamedTuples) to `CBET`; their defaults are
Mach 2 outflow outside (the near-sonic LEH flow turns supersonic outside, Michel et al.
2013) and Mach 0 inside. Both are engineering defaults; neither is measured. In a low-fill
hohlraum the net transfer is the balance of the two. In Pickworth et al. 2020's HYDRA
runs, power is transferred from the inner to the outer cones outside the LEH and from the
outer to the inner cones inside. For this reason Δλ = 0 does not give zero transfer in
such a hohlraum; I-Raum needed ~0.5 Å to cancel flow-driven CBET (NC24).
[`n221204_cbet`](@ref Briareus.Benchmarks.n221204_cbet) is a version calibrated to the N221204 HYDRA transfer.
"""
struct TwoZoneCBET <: LaserPlasmaModel
    outside::CBET
    inside::CBET
end

_zshift(z::Real, d) = z + d
_zshift(z, d) = t -> _tvalue(z, t) + d

function TwoZoneCBET(; leh_z, offset::Real = 0.5, outside = (;), inside = (;))
    return TwoZoneCBET(CBET(; leh_z = _zshift(leh_z, offset), mach = 2.0, outside...),
                       CBET(; leh_z = _zshift(leh_z, -offset), mach = 0.0, inside...))
end

apply_lpi(m::TwoZoneCBET, beams::AbstractVector{<:Beam}, time = 0.0) =
    apply_lpi(m.inside, apply_lpi(m.outside, beams, time), time)

"""
    pulse_rise(pulse; lookback = 1.0, delta = 0.15) -> t

Start of the rise to peak power (ns) of `pulse`, a [`PulseShape`](@ref) or
`(times, powers)`. The start is the knee where the tangent at the half-maximum crossing
t½ (slope over t½ ± `delta`) meets the power `lookback` ns before t½, the step the rise
leaves. It is
invariant to a scale factor on the power, so delivered and HYDRA-degraded inputs give
nearly the same time: 4.11 ns on N221204's delivered pulse ([`n221204_pulse`](@ref Briareus.Benchmarks.n221204_pulse)),
4.0 on PRE 109 Fig. 9's degraded input powers, where the HYDRA deck switches its CBET
clamp ("until the beginning of the rise to peak power", 3.9–4.1 ns). This is an
engineering definition for pulses with one main rise.
"""
pulse_rise(p::PulseShape; kw...) = pulse_rise((p.t, p.p); kw...)
function pulse_rise(p::Tuple{<:AbstractVector,<:AbstractVector}; lookback = 1.0, delta = 0.15)
    t, P = p
    half = maximum(P) / 2
    k = findfirst(>=(half), P)
    (k === nothing || k == 1) && throw(ArgumentError("the pulse has no rise to its maximum"))
    th = t[k-1] + (half - P[k-1]) / (P[k] - P[k-1]) * (t[k] - t[k-1])
    slope = (_interp(t, P, th + delta) - _interp(t, P, th - delta)) / 2delta
    return th - (half - _interp(t, P, th - lookback)) / slope
end

"""
    cbet_clamp_schedule(pulse; foot = 0.002, peak = 0.0075, width = 0.2) -> (times, values)

CBET saturation clamp δn/n stepping from `foot` to `peak` over `width` (ns) centered on
the start of the rise of `pulse` to peak power ([`pulse_rise`](@ref)), as a
two-knot `(times, values)` table for `CBET`'s `clamp` (constant beyond the knots, so
fitting `:clamp` fits the foot and peak values). The defaults are PRE 109's HYDRA setting
for N221204 (0.002 until the rise, 0.0075 after; PoP 31, 070502 quotes 0.008), which
on PRE 109 Fig. 9's input powers puts the switch at 3.9–4.1 ns as in the deck (4.0–4.2 ns
on the delivered pulse). The schedule is a HYDRA-equivalent setting; it is not derived
from physics. With `saturation = :min` the values mean what they mean in HYDRA.
"""
function cbet_clamp_schedule(pulse; foot = 0.002, peak = 0.0075, width = 0.2)
    tr = pulse_rise(pulse)
    return ([tr - width / 2, tr + width / 2], [foot, peak])
end

"""
    cbet_null(model, beams, t; receivers, lo = 0.0, hi = 1.5, tol = 1e-4) -> dlambda

The detuning (Å at the beam wavelength; × 3 for NIF's 1ω convention) of the `receivers`
groups, all other beams at zero, for which `model` leaves the total power of the receivers
unchanged at time `t`. This Δλ cancels flow-driven transfer ([`TwoZoneCBET`](@ref)).
NC24's I-Raum needed ~0.5 Å at 1ω "mostly to negate" it. Their 1–1.55 Å for HYBRID-E is
the Δλ of round hot spots, which also compensates bubble absorption, so it lies above the
null. The null is found by bisection on [`receiver_ratio`](@ref) over [`lo`, `hi`] and is
`NaN` when the ratio does not cross 1 on that interval.
"""
function cbet_null(m::LaserPlasmaModel, beams::AbstractVector{<:Beam}, t; receivers,
                   lo = 0.0, hi = 1.5, tol = 1e-4)
    f(x) = receiver_ratio(m, beams, x; receivers, time = t) - 1
    a, b = float(lo), float(hi)
    fa, fb = f(a), f(b)
    sign(fa) == sign(fb) && return NaN
    while b - a > tol
        c = (a + b) / 2
        fc = f(c)
        sign(fc) == sign(fa) ? ((a, fa) = (c, fc)) : ((b, fb) = (c, fc))
    end
    return a - fa * (b - a) / (fb - fa)
end

"""
    receiver_ratio(model, beams, dlambda; receivers, time = 0.0)

Ratio of the `receivers` groups' total power after/before `model` when those groups are
detuned by `dlambda` (Å, beam wavelength) and all other beams by zero.
"""
function receiver_ratio(m::LaserPlasmaModel, beams::AbstractVector{<:Beam}, dλ; receivers,
                        time = 0.0)
    rs = Set(Symbol[receivers...])
    bs = [remake(b; dlambda = (b.group in rs ? one(dλ) : zero(dλ)) * dλ) for b in beams]
    eff = apply_lpi(m, bs, time)
    return sum(b.power for b in eff if b.group in rs) / sum(b.power for b in bs if b.group in rs)
end

"""
    calibrate_cbet(model, beams, data; receivers, fit = (:multiplier,), time = 0.0,
                   bounds = Dict(:multiplier => (1e-3, 10.0), :mach => (0.0, 3.0),
                                 :clamp => (1e-4, 0.05))) -> (model, result)

Fit CBET calibration parameters to measured or reference transfer data. `data` is a list
of `(dlambda, ratio)` pairs: with the `receivers` groups (e.g. the inner cones) detuned by
`dlambda` (Å, beam wavelength — 1ω-quoted NIF values ÷ 3), their total power should be
multiplied by `ratio` (e.g. 1.6 for "inner +60%"). The constant-valued fields in `fit`
(any of `:multiplier`, `:mach`, `:clamp`) are fitted by bounded L-BFGS-B on squared
log-ratio residuals with automatic-differentiation gradients. Returns the calibrated model
and the Optim result.
"""
function calibrate_cbet(m::Union{CBET,TwoZoneCBET}, beams::AbstractVector{<:Beam}, data;
                        receivers, fit = (:multiplier,), time = 0.0,
                        bounds = _CBET_BOUNDS, iterations = 200)
    loss(mm) = sum((log(receiver_ratio(mm, beams, dλ; receivers, time)) - log(r))^2
                   for (dλ, r) in data)
    return _calibrate(loss, m, fit, bounds, iterations)
end

"""
    calibrate_cbet(model, beams_at, data; receivers, fit = (:multiplier,), bounds,
                   iterations = 200) -> (model, result)

Multi-epoch calibration, which fits `model` (a [`CBET`](@ref) or [`TwoZoneCBET`](@ref))
to a time series of post-transfer power fractions. `beams_at(t)` gives the setup beams at
time `t` (ns); `data` lists `(t, fraction)` or `(t, fraction, weight)` rows, `fraction`
being the `receivers` groups' share of the total power after `model` (e.g. the inner
cone fraction of PRE 109 Fig. 9b). The fit minimizes the weighted sum of squared
fraction residuals.

`fit` names the parameters: fields of a `CBET` (`:multiplier`, `:mach`, `:flow`,
`:clamp`, `:Te`, `:Ti`, `:ne_nc`), or for a `TwoZoneCBET` `(zone, field)` pairs with
zone `:outside`, `:inside` or `:both` (one value shared by the two zones; a bare
field name means `:both`). A field holding a number is one parameter; a field holding
a `(times, values)` table fits every value at its fixed time, so plasma-state
histories are calibration inputs. `bounds` maps fields to `(lo, hi)`; parameters
with `lo > 0` are fitted in log space. The optimizer is bounded L-BFGS-B with
forward-mode gradients.
"""
function calibrate_cbet(m::Union{CBET,TwoZoneCBET}, beams_at::Function, data; receivers,
                        fit = (:multiplier,), bounds = _CBET_BOUNDS, iterations = 200)
    rs = Set(Symbol[receivers...])
    rows = [(Float64(r[1]), Float64(r[2]), length(r) > 2 ? Float64(r[3]) : 1.0) for r in data]
    beams = [beams_at(t) for (t, _, _) in rows]
    function loss(mm)
        acc = 0.0
        for ((t, f, w), bs) in zip(rows, beams)
            eff = apply_lpi(mm, bs, t)
            acc += w * (sum(b.power for b in eff if b.group in rs) / sum(b.power for b in eff) - f)^2
        end
        return acc
    end
    return _calibrate(loss, m, fit, bounds, iterations)
end

"""
    calibrate_cbet(model, shots::AbstractVector; receivers, fit = (:multiplier,), bounds,
                   penalty = nothing, iterations = 200) -> (model, result)

Joint multi-shot calibration, in which one `model` (shared plasma physics) is fitted to
the post-transfer power fractions of several shots at once. Each entry of `shots` is a
NamedTuple `(; beams_at, data, delay = 0.0)`: `beams_at(t)` gives the setup beams of
that shot (its powers and detuning) and `data` its `(t, fraction[, weight])` rows, as
in the multi-epoch method. `delay` (ns) shifts the shot against the time tables of the
model: the model is evaluated at `t - delay`, so a shot whose rise to peak comes 0.2 ns
earlier than that of the reference pulse has `delay = -0.2` and is evaluated with the
plasma history (and clamp schedule) of the reference 0.2 ns early. `penalty(model)`
(optional) is added to the loss, e.g. a constraint on the Δλ null; it must be generic in
the number type. `fit` and `bounds` are as in the multi-epoch method.
"""
function calibrate_cbet(m::Union{CBET,TwoZoneCBET}, shots::AbstractVector; receivers,
                        fit = (:multiplier,), bounds = _CBET_BOUNDS, penalty = nothing,
                        iterations = 200)
    loss = _joint_loss(shots, receivers, penalty)
    return _calibrate(loss, m, fit, bounds, iterations)
end

# Loss of the joint calibration, the weighted squared fraction residuals over all shots
# plus the penalty. Beams are built once.
function _joint_loss(shots, receivers, penalty)
    rs = Set(Symbol[receivers...])
    sets = map(shots) do s
        rows = [(Float64(r[1]), Float64(r[2]), length(r) > 2 ? Float64(r[3]) : 1.0) for r in s.data]
        (rows, [s.beams_at(t) for (t, _, _) in rows], Float64(get(s, :delay, 0.0)))
    end
    function loss(mm)
        acc = penalty === nothing ? 0.0 : penalty(mm)
        for (rows, beams, delay) in sets, ((t, f, w), bs) in zip(rows, beams)
            eff = apply_lpi(mm, bs, t - delay)
            acc += w * (sum(b.power for b in eff if b.group in rs) / sum(b.power for b in eff) - f)^2
        end
        return acc
    end
    return loss
end

const _CBET_BOUNDS = Dict(:multiplier => (1e-3, 10.0), :mach => (0.0, 3.0),
                          :flow => (-1000.0, 1000.0), :clamp => (1e-4, 0.05),
                          :Te => (300.0, 2e4), :Ti => (100.0, 1e4), :ne_nc => (1e-3, 0.25))
const _CBET_FITTABLE = (:multiplier, :mach, :flow, :clamp, :Te, :Ti, :ne_nc)

# normalized (zone, field) specs of a fit
_fit_specs(::CBET, fit) = [(:self, f isa Symbol ? f : Symbol(last(f))) for f in fit]
_fit_specs(::TwoZoneCBET, fit) =
    [f isa Symbol ? (:both, f) : (Symbol(first(f)), Symbol(last(f))) for f in fit]
_zone(m::CBET, ::Symbol) = m
function _zone(m::TwoZoneCBET, z::Symbol)
    z === :inside && return m.inside
    z in (:outside, :both) && return m.outside
    throw(ArgumentError("zone must be :outside, :inside or :both (got $z)"))
end

# Parameter values of a field, which are a number or the values of a (times, values) table.
function _field_values(v, f)
    v isa Real && return [Float64(v)]
    v isa Tuple{<:AbstractVector,<:AbstractVector} && return Float64.(collect(v[2]))
    throw(ArgumentError("$f must be a number or a (times, values) table to be fitted (got $(typeof(v)))"))
end

function _calibrate(loss, m, fit, bounds, iterations)
    u0, ulo, uhi, with = _calibration_space(m, fit, bounds)
    f(u) = loss(with(u))
    g!(G, u) = ForwardDiff.gradient!(G, f, u)
    res = Optim.optimize(f, g!, ulo, uhi, u0, Optim.LBFGSB(), Optim.Options(; iterations))
    return with(Optim.minimizer(res)), res
end

# Parameter space of the fit, with the start, the bounds (log space where lo > 0) and the
# map from u to the model. It is also used for Laplace (curvature) uncertainties of a fit.
function _calibration_space(m, fit, bounds)
    specs = _fit_specs(m, fit)
    for (_, f) in specs
        f in _CBET_FITTABLE || throw(ArgumentError("cannot fit $f"))
    end
    x0, lo, hi, sizes = Float64[], Float64[], Float64[], Int[]
    for (z, f) in specs
        v = _field_values(getfield(_zone(m, z), f), f)
        append!(x0, v)
        append!(lo, fill(bounds[f][1], length(v)))
        append!(hi, fill(bounds[f][2], length(v)))
        push!(sizes, length(v))
    end
    logs = lo .> 0
    to_u(x) = [l ? log(xi) : xi for (xi, l) in zip(x, logs)]
    from_u(u) = [l ? exp(ui) : ui for (ui, l) in zip(u, logs)]
    with(u) = _cbet_set(m, specs, sizes, from_u(u))
    return to_u(clamp.(x0, lo, hi)), to_u(lo), to_u(hi), with
end

# the model with the fitted fields set from x (numbers, or the values of tables)
function _cbet_set(m, specs, sizes, x)
    vals = Dict(z => Dict{Symbol,Any}() for z in (:self, :outside, :inside))
    k = 0
    for ((z, f), n) in zip(specs, sizes)
        old = getfield(_zone(m, z), f)
        v = old isa Real ? x[k+1] : (old[1], x[k+1:k+n])
        k += n
        for zz in (z === :both ? (:outside, :inside) : (z,))
            vals[zz][f] = v
        end
    end
    m isa CBET && return _cbet_with(m, vals[:self])
    return TwoZoneCBET(_cbet_with(m.outside, vals[:outside]), _cbet_with(m.inside, vals[:inside]))
end

function _cbet_with(m::CBET, vals::AbstractDict)
    kw = Dict{Symbol,Any}(f => getfield(m, f) for f in fieldnames(CBET))
    merge!(kw, vals)
    return CBET(; kw...)
end
