# Reduced cross-beam energy transfer (CBET) at the laser entrance holes.
#
# Coupling: kinetic ray-pair gain (Nguyen et al., arXiv:2304.13268, Eqs. 1–5; identical
# to Michel et al. 2009's γ = (Δk²/8k₀)Γ), in practical units
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

Fried–Conte plasma dispersion function for real argument,
`Z(ξ) = −2 D(ξ) + i√π exp(−ξ²)` (D: Dawson's integral). Differentiable.
"""
plasma_dispersion(ξ::Real) = complex(-2 * dawson(ξ), sqrt(π) * exp(-ξ^2))

"""
    cbet_coupling(Ω, k, kdotu; ne, Te, Ti, ions) -> Γ

Ion-acoustic response `Γ = χe(1+χi)/(1+χe+χi)` for beat frequency `Ω` (rad/s), beat
wavenumber `k` (1/m), Doppler shift `k·u` (rad/s), electron density `ne` (m⁻³),
temperatures `Te`, `Ti` (eV), and ion species `ions` — `(fraction, Z, A)` tuples with
number fractions summing to one; `χi` sums over species (each with its own Landau
damping, which matters for CH: the protons dominate).
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
  Michel et al. 2009; `Ti` is not given there and is an assumption.
- `flow` — outflow speed along the hohlraum axis (km/s, away from the hohlraum), or
  `mach` times the ion sound speed `√((Z̄ Te + 3 Ti)/(Ā m_p))` when `flow` is nothing.
- `length` — maximum interaction length (mm). A pair crossing at angle θ interacts over
  `min(length, w/sin θ)`, with `w` the mean beam width at the LEH plane.
- `polarization` — `:smoothed` (NIF polarization smoothing, uncorrelated speckle:
  gain × (1+cos²θ)/4) or `:parallel` (co-polarized, gain × 1).
- `multiplier` — calibration factor on all gains. Fit it to shot data: published NIF
  designs tune CBET with saturation clamps and power multipliers.
- `clamp` — ion-wave saturation level δn/n (`nothing` = linear theory, which grossly
  overpredicts transfer at NIF intensities). Each pair's coupling is multiplied by
  `1/√(1 + (δn/clamp)²)`, a smooth version of capping the linear δn/n at `clamp`.
  Published NIF calibrations range from 3×10⁻⁴ (2011) to 0.002 (foot) / 0.0075 (peak)
  for N221204 (Kritcher et al., PRE 109, 025204, 2024); default 0.002.

Any plasma parameter and `multiplier` may be a number, a function of time (ns), or a
`(times, values)` table. Each pair's exchange uses the geometric-mean beam area at the
LEH, so total power is conserved. The coupled exchange is integrated with `nsteps` RK4
steps. Beam detuning comes from `Beam.dlambda` (Å at the beam's own wavelength — NIF
quotes Δλ at 1ω, so divide those values by 3 for 3ω beams).

The pair coefficients, and the two LEHs' exchanges, run on threads; results are
bitwise independent of the thread count.
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
    nsteps::Int = 20
end

"Plasma parameters of `m` at time `t` as a NamedTuple (SI densities, eV, m/s flow)."
function cbet_plasma(m::CBET, λ0_µm, t)
    ne_nc, Te, Ti = _tvalue(m.ne_nc, t), _tvalue(m.Te, t), _tvalue(m.Ti, t)
    Z = sum(f * Z for (f, Z, _) in m.ions)
    A = sum(f * A for (f, _, A) in m.ions)
    ω0 = 2π * _C_LIGHT / (λ0_µm * 1e-6)
    nc = _EPS0 * _ME * ω0^2 / _QE^2
    cs = sqrt((Z * Te + 3Ti) * _QE / (A * _MP))
    u = m.flow === nothing ? _tvalue(m.mach, t) * cs : 1e3 * _tvalue(m.flow, t)
    return (; ne = ne_nc * nc, nc, Te, Ti, ions = m.ions, u, cs)
end

"""
    cbet_gain(m, b, c, t) -> g

Linear (unsaturated) spatial gain rate (1/mm per TW/mm² of beam `c`) of beam `b` due to
beam `c` crossing at the LEH, at time `t`, including the polarization factor and
`multiplier`. Antisymmetric: `cbet_gain(m, c, b, t) == -cbet_gain(m, b, c, t)`.
"""
cbet_gain(m::CBET, b::Beam, c::Beam, t = 0.0) = _cbet_pair(m, b, c, t)[1]

# (linear gain g, ion-wave amplitude coefficient D) for a pair: δn/n = D √(I_b I_c),
# intensities in TW/mm². δn/n = (k²c²/2ω_pe²) |Γ| |a_b||a_c| with |a|² = 7.31e-19 I λ²
# (I in W/cm², λ in µm; linear polarization), and polarization smoothing reducing the
# rms beat amplitude by √pol.
function _cbet_pair(m::CBET, b::Beam, c::Beam, t)
    pl = cbet_plasma(m, b.wavelength, t)
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
    # outflow along the hohlraum axis, away from the hohlraum at this LEH
    uvec = Vec3(0, 0, b.port[3] > 0 ? 1.0 : -1.0) * pl.u
    Γ = cbet_coupling(ωb - ωc, kmag, dot(kvec, uvec); pl.ne, pl.Te, pl.Ti, pl.ions)
    pol = m.polarization === :smoothed ? (1 + cosθ^2) / 4 :
          m.polarization === :parallel ? 1.0 :
          throw(ArgumentError("polarization must be :smoothed or :parallel"))
    # 4.59e-18 µm⁻¹ per W/cm²  →  ×1e3 (per mm) ×1e14 (per TW/mm²)
    g = -0.459 * b.wavelength * sin2h * imag(Γ) * pol * _tvalue(m.multiplier, t)
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

# `threaded = false` is the serial reference. Both paths give bitwise identical results
# for any thread count: every pair writes its own slots of K and D, the two LEHs never
# share data, and no sum runs in a thread-dependent order.
function _apply_cbet(m::CBET, beams::AbstractVector{<:Beam}, time, threaded::Bool)
    T = promote_type(mapreduce(b -> typeof(b.power), promote_type, beams),
                     mapreduce(b -> typeof(b.dlambda), promote_type, beams),
                     mapreduce(b -> eltype(b.aim), promote_type, beams))
    # model parameters may carry dual numbers (calibration)
    for f in (:ne_nc, :Te, :Ti, :mach, :flow, :multiplier, :clamp)
        v = getfield(m, f)
        v === nothing || (T = promote_type(T, typeof(_tvalue(v, time))))
    end
    # every pair coefficient is built from T-typed inputs and Float64 constants
    Tk = promote_type(T, Float64)
    P = T[b.power for b in beams]
    zL = _tvalue(m.leh_z, time)
    clampv = m.clamp === nothing ? nothing : _tvalue(m.clamp, time)
    idxs = Vector{Int}[]
    for top in (true, false)
        idx = [i for (i, b) in enumerate(beams) if (b.port[3] > 0) == top]
        length(idx) < 2 || push!(idxs, idx)
    end
    isempty(idxs) && return [remake(b; power = P[i]) for (i, b) in enumerate(beams)]
    geos = [[_beam_at_plane(beams[i], beams[i].port[3] > 0 ? zL : -zL) for i in idx]
            for idx in idxs]
    # exchange coefficients K[p,q] = g_pq ℓ_pq / Ā_pq (antisymmetric) and ion-wave
    # amplitude coefficients D[p,q] (symmetric), one pair of matrices per LEH
    Ks = [zeros(Tk, length(idx), length(idx)) for idx in idxs]
    Ds = [zeros(Tk, length(idx), length(idx)) for idx in idxs]
    tasks = [(l, p, q) for l in eachindex(idxs) for p in eachindex(idxs[l])
             for q in p+1:length(idxs[l])]
    _cbet_couplings!(Ks, Ds, m, beams, idxs, geos, tasks, time, threaded)
    Pls = [P[idx] for idx in idxs]
    Aeffs = [[g[2] for g in geo] for geo in geos]
    _cbet_exchanges!(Pls, Ks, Ds, Aeffs, clampv, m.nsteps, threaded)
    for (idx, Pl) in zip(idxs, Pls)
        P[idx] = Pl
    end
    return [remake(b; power = P[i]) for (i, b) in enumerate(beams)]
end

# Function barrier: the coefficient eltype is a type parameter here, so the threaded
# closure does not box it.
function _cbet_couplings!(Ks::Vector{Matrix{Tk}}, Ds::Vector{Matrix{Tk}}, m::CBET, beams,
                          idxs, geos, tasks, time, threaded::Bool) where {Tk}
    function pair!(j)
        l, p, q = tasks[j]
        idx, geo, K, D = idxs[l], geos[l], Ks[l], Ds[l]
        b, c = beams[idx[p]], beams[idx[q]]
        cosθ = clamp(dot(direction(b), direction(c)), -1.0, 1.0)
        sinθ = sqrt(max(1 - cosθ^2, 1e-12))
        ℓ = min(m.length, (geo[p][1] + geo[q][1]) / 2 / sinθ)
        Ā = sqrt(geo[p][2] * geo[q][2])
        g, d = _cbet_pair(m, b, c, time)
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

# RK4 integration of dP_p/dτ = P_p Σ_q K_pq S_pq P_q over τ ∈ [0, 1] in place, with the
# smooth saturation factor S_pq = 1/√(1 + (δn_pq/clamp)²) symmetric in p, q so power
# is conserved. Work buffers are allocated once per call.
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

function _cbet_rate!(out, x, K, D, Aeff, clampv, I)
    n = length(x)
    I .= x ./ Aeff
    for p in 1:n
        acc = zero(eltype(x))
        for q in 1:n
            q == p && continue
            δn = D[p, q] * sqrt(max(I[p] * I[q], zero(I[p])))
            acc += K[p, q] * x[q] / sqrt(1 + (δn / clampv)^2)
        end
        out[p] = x[p] * acc
    end
    return out
end

"""
    ion_wave_amplitude(m::CBET, b, c, t = 0.0) -> δn/n

Linear ion-acoustic density perturbation driven by beams `b` and `c` at their LEH
plane intensities — compare with `m.clamp`.
"""
function ion_wave_amplitude(m::CBET, b::Beam, c::Beam, t = 0.0)
    zL = _tvalue(m.leh_z, t) * (b.port[3] > 0 ? 1 : -1)
    _, d = _cbet_pair(m, b, c, t)
    Ib = b.power / _beam_at_plane(b, zL)[2]
    Ic = c.power / _beam_at_plane(c, zL)[2]
    return d * sqrt(Ib * Ic)
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
of `(dlambda, ratio)` pairs: with the `receivers` groups (e.g. the inner cones) detuned
by `dlambda` (Å, beam wavelength — 1ω-quoted NIF values ÷ 3), their total power should be
multiplied by `ratio` (e.g. 1.6 for "inner +60%"). Fits the constant-valued fields in
`fit` (any of `:multiplier`, `:mach`, `:clamp`) by bounded L-BFGS-B on squared log-ratio
residuals with automatic-differentiation gradients. Returns the calibrated model and the
Optim result.
"""
function calibrate_cbet(m::CBET, beams::AbstractVector{<:Beam}, data; receivers,
                        fit = (:multiplier,), time = 0.0,
                        bounds = Dict(:multiplier => (1e-3, 10.0), :mach => (0.0, 3.0),
                                      :clamp => (1e-4, 0.05)))
    fields = Symbol[fit...]
    for f in fields
        f in (:multiplier, :mach, :clamp) || throw(ArgumentError("cannot fit $f"))
        v = getfield(m, f)
        v isa Real || throw(ArgumentError("$f must be a number to be fitted (got $(typeof(v)))"))
    end
    with(x) = _cbet_with(m, Dict(f => x[i] for (i, f) in enumerate(fields)))
    loss(x) = sum((log(receiver_ratio(with(x), beams, dλ; receivers, time)) - log(r))^2
                  for (dλ, r) in data)
    x0 = Float64[getfield(m, f) for f in fields]
    lo = Float64[bounds[f][1] for f in fields]
    hi = Float64[bounds[f][2] for f in fields]
    g!(G, x) = ForwardDiff.gradient!(G, loss, x)
    res = Optim.optimize(loss, g!, lo, hi, clamp.(x0, lo, hi), Optim.LBFGSB(),
                         Optim.Options(; iterations = 200))
    return with(Optim.minimizer(res)), res
end

function _cbet_with(m::CBET, vals::Dict)
    return CBET(; leh_z = m.leh_z, ne_nc = m.ne_nc, Te = m.Te, Ti = m.Ti, ions = m.ions,
                mach = get(vals, :mach, m.mach), flow = m.flow, length = m.length,
                polarization = m.polarization, multiplier = get(vals, :multiplier, m.multiplier),
                clamp = get(vals, :clamp, m.clamp), nsteps = m.nsteps)
end
