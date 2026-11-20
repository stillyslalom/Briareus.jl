# Why a backed diagnostic window lets flux out later than the bare heat front predicts:
# 1-D planar Lagrangian radiation hydrodynamics of a gold foil, bare or backed, and the
# burnthrough of a low-Z backing by the flux the gold transmits. Research notes,
# "Backed windows: what delays the emergent flux (2026-11-18)";
# reference/litreview/window_burnthrough.md.
#
#   julia --project=. examples/tamped_foil.jl nova            # ~2 min
#   julia --project=. examples/tamped_foil.jl table [scale]   # ~3 min; scale of the drive (1.0)
#
# The gold follows the power laws of Hammer & Rosen, Phys. Plasmas 10, 1829 (2003), with e =
# 3.4 T^1.6 ρ^-0.14 MJ/g, κ_R = 7200 T^-1.5 ρ^0.2 cm²/g and T in heV, and P = 0.25 e/V. The
# model uses flux-limited diffusion and a Marshak boundary at the incident temperature T_D
# on the front, as in `examples/late_time_wall.jl`, which validates this wall against the
# self-similar solution of Hammer & Rosen. In this example the wall is a foil of finite
# areal mass. Its rear radiates into a cold, black absorber, which is vacuum or a backing
# before the backing burns through. The rear is either free (bare) or pushes on a backing of
# given areal mass. The backing is an ideal gas (γ = 5/3) that does not absorb radiation, so
# it acts on the gold only through its hydrodynamic hold.
#
# The flux F that leaves the rear of the gold burns through the backing at the measured
# x-ray ablation rate of the backing material, ṁ = A (F/σ(200 eV)⁴)^(3/4). This is the
# approximate T³ scaling of the rate with the drive temperature, taken as the temperature
# whose blackbody flux is F. A, in mg/cm²/ns at 200 eV, is fitted in this example to the
# burnthrough ablation rates of Olson et al. (SAND2007-7414C, OSTI 1147322; Phys. Plasmas
# 18, 032706 (2011)), read from their figure. It is 4.8 for polyimide (C₂₂H₁₀N₂O₅, the
# measured material closest to an epoxy) and 5.9 for beryllium, and ×0.7–1.3 brackets the
# scatter and the material difference. The flux emerges when the areal mass of the backing
# has burned through.
#
# `nova` runs bare and backed foils of 0.77–3 µm on the Nova drive of 260 eV (t/1 ns)^0.15
# against the measured half-signal times (Lindl et al., Phys. Plasmas 11, 339 (2004), Fig.
# 4-6) and the arrival of `HeatFront`. `table` runs 1.5–3 µm with 35 µm of epoxy (1.185
# g/cm³; Lester et al., Phys. Rev. E 114, L023201 (2026)) or 50 µm of Be on a public
# ignition-scale drive. The drive is the Dante 1 temperature of N221204 from 3.3 ns (Pak et
# al., Phys. Rev. E 109, 025203 (2024) Fig. 1). Before that it is the temperature of N210601
# (Zylstra et al., Phys. Rev. E 106, 025202 (2022) Fig. 2) scaled by 0.958 to join it, after
# a linear rise over 0.2–0.8 ns to the first value of the foot. `scale` multiplies the
# drive, since the incident temperature at the window is a few percent below that of Dante.
# The columns are:
# - the arrival of `HeatFront`;
# - two calibrations of it to the 7.5 ns leading edge that Lester et al. measured through
#   2.5 µm of gold and 35 µm of epoxy, a heat-front coefficient and a constant offset,
#   each fitted on this drive;
# - the transmission onset (1%) of the 1-D gold;
# - the emergence with the arrival of `HeatFront` and the transition and transmission of
#   `BurnthroughAlbedo` (`RosselandTransmission`) driving the backing;
# - the 1-D emergence.

using Briareus, Printf, LinearAlgebra

const σ1 = SIGMA_SB * 1e8 * 1e2 * 1e-3       # σ (1 heV)^4: TW/mm² → MJ/(ns cm²)
const r_eos = 0.25
const c2 = 1e-5                               # (cm/ns)² per MJ/g
eT(T, V) = 3.4 * T^1.6 * V^0.14
Te(e, V) = (max(e, 1e-12) / (3.4 * V^0.14))^(1 / 1.6)
κR(T, V) = 7200 * T^-1.5 * V^-0.2
cvT(T, V) = 3.4 * 1.6 * T^0.6 * V^0.14

# zones by areal mass from the front: `ng` gold zones graded from `dm0` (capped), then `nt`
# equal backing zones; nodes 1 (front) to n+1 (rear), both free surfaces
mutable struct Foil
    dm::Vector{Float64}; x::Vector{Float64}; u::Vector{Float64}
    V::Vector{Float64}; T::Vector{Float64}; e::Vector{Float64}
    ng::Int; γ::Float64
end
function Foil(; mass, ng = 160, dm0 = 2e-6, backing = 0.0, ρb = 1.185, nt = 60, γ = 5 / 3)
    grade(q) = [min(dm0 * q^(i - 1), 4mass / ng) for i in 1:ng]
    lo, hi = 1.0, 1.2
    for _ in 1:200                            # the ratio that fills the gold's areal mass
        q = (lo + hi) / 2
        sum(grade(q)) > mass ? (hi = q) : (lo = q)
    end
    dmg = grade((lo + hi) / 2)
    dmg .*= mass / sum(dmg)
    dmb = backing > 0 ? fill(backing / nt, nt) : Float64[]
    dm = [dmg; dmb]
    V = [fill(1 / 19.3, ng); fill(1 / ρb, length(dmb))]
    T = fill(0.02, length(dm))
    e = [eT.(T[1:ng], V[1:ng]); zeros(length(dmb))]
    return Foil(dm, [0.0; cumsum(dm .* V)], zeros(length(dm) + 1), V, T, e, ng, γ)
end
pressure(f, i) = i <= f.ng ? r_eos * f.e[i] / f.V[i] : (f.γ - 1) * max(f.e[i], 0.0) / f.V[i]
sound2(f, i) = i <= f.ng ? (1 + r_eos) * r_eos * f.e[i] * c2 : f.γ * (f.γ - 1) * max(f.e[i], 0.0) * c2
hydro_dt(f) = minimum(0.2 * (f.x[i+1] - f.x[i]) / (sqrt(sound2(f, i)) + abs(f.u[i+1] - f.u[i]) + 1e-30)
                      for i in eachindex(f.dm))

# explicit Lagrangian step (quadratic viscosity), P = 0 beyond both surfaces
function hydro!(f, dt)
    n = length(f.dm)
    Pt = [pressure(f, i) + (f.u[i+1] < f.u[i] ? 2.0 * (f.u[i+1] - f.u[i])^2 / f.V[i] / c2 : 0.0) for i in 1:n]
    for j in 1:n+1
        mj = j == 1 ? f.dm[1] / 2 : j == n + 1 ? f.dm[n] / 2 : (f.dm[j-1] + f.dm[j]) / 2
        f.u[j] -= c2 * ((j == n + 1 ? 0.0 : Pt[j]) - (j == 1 ? 0.0 : Pt[j-1])) / mj * dt
    end
    f.x .+= f.u .* dt
    for i in 1:n
        Vn = (f.x[i+1] - f.x[i]) / f.dm[i]
        f.e[i] -= Pt[i] * (Vn - f.V[i])
        f.V[i] = Vn
        i <= f.ng && (f.T[i] = Te(f.e[i], Vn))
    end
end

# backward-Euler diffusion in y = T⁴ over the gold (Picard; sum flux limiter): a Marshak
# boundary at T_D in front, a vacuum boundary behind; returns the flux out of the rear
function rad!(f, TD, dt; tol = 1e-8)
    n = f.ng
    T, V, dm = view(f.T, 1:n), view(f.V, 1:n), view(f.dm, 1:n)
    eold = f.e[1:n]
    Fr = 0.0
    for _ in 1:60
        y = T .^ 4
        A = Tridiagonal(zeros(n - 1), zeros(n), zeros(n - 1))
        rhs = zeros(n)
        for i in 1:n
            dedy = cvT(T[i], V[i]) / (4T[i]^3)
            A.d[i] += dm[i] * dedy / dt
            rhs[i] += dm[i] * (dedy * y[i] - (eT(T[i], V[i]) - eold[i])) / dt
        end
        g0 = 2 / (1 + 3κR(T[1], V[1]) * dm[1] / 4)
        A.d[1] += σ1 * g0
        rhs[1] += σ1 * g0 * TD^4
        gn = 2 / (1 + 3κR(T[n], V[n]) * dm[n] / 4)
        A.d[n] += σ1 * gn
        for i in 1:n-1
            τ = (κR(T[i], V[i]) * dm[i] + κR(T[i+1], V[i+1]) * dm[i+1]) / 2
            c = (4σ1 / 3) / (τ + (4 / 3) * abs(y[i+1] - y[i]) / ((y[i] + y[i+1]) / 2))
            A.d[i] += c; A.du[i] -= c; A.d[i+1] += c; A.dl[i] -= c
        end
        Tn = max.(A \ rhs, 1e-12) .^ 0.25
        δ = maximum(abs.(Tn .- T) ./ Tn)
        T .= Tn
        Fr = σ1 * gn * T[n]^4
        δ < tol && break
    end
    for i in 1:n
        f.e[i] = eT(T[i], V[i])
    end
    return Fr
end

# a foil on the drive TD(t) (heV) to `tend`: times, rear flux, transmitted fraction, and the
# density of the gold's rear quarter (by mass), every `out` ns
function foil_run(TD, tend; out = 0.01, kw...)
    f = Foil(; kw...)
    t, nxt = 0.0, out
    ts, Fr, tr, ρ = Float64[], Float64[], Float64[], Float64[]
    while t < tend - 1e-12
        dt = min(2e-3, hydro_dt(f), tend - t)
        hydro!(f, dt)
        F = rad!(f, max(TD(t + dt), 0.02), dt)
        t += dt
        if t >= nxt - 1e-9
            cm = cumsum(f.dm[1:f.ng])
            ks = findall(>=(0.75cm[end]), cm)
            Fi = σ1 * TD(t)^4
            push!(ts, t); push!(Fr, F); push!(tr, Fi > σ1 * 0.5^4 ? F / Fi : 0.0)
            push!(ρ, sum(f.dm[ks]) / sum(f.dm[ks] .* f.V[ks]))
            nxt += out
        end
    end
    return (; t = ts, Fr, tr, ρ)
end
crossing(t, y, level) = (k = findfirst(>=(level), y); k === nothing ? Inf :
                         k == 1 ? t[1] : t[k-1] + (level - y[k-1]) / (y[k] - y[k-1]) * (t[k] - t[k-1]))

# the time a backing of areal mass `mb` (mg/cm²) burns through at ṁ = A (F/σ(200 eV)⁴)^(3/4),
# F(t) the flux on it sampled at `ts` (MJ/(ns cm²), trapezoidal in ṁ)
function backing_burnthrough(ts, F, mb, A)
    rate(x) = A * (max(x, 0.0) / (σ1 * 2.0^4))^0.75
    I = 0.0
    for k in 2:length(ts)
        dI = (ts[k] - ts[k-1]) * (rate(F[k]) + rate(F[k-1])) / 2
        I + dI >= mb && return ts[k-1] + (mb - I) / dI * (ts[k] - ts[k-1])
        I += dI
    end
    return Inf
end

# measured x-ray ablation rate coefficients (mg/cm²/ns at 200 eV) and the backings
const RATE = (epoxy = 4.8, Be = 5.9)
const BACKINGS = (("epoxy 35 µm", 1.185, 35.0, RATE.epoxy), ("Be 50 µm", 1.85, 50.0, RATE.Be))

interp(ts, ys, x) = x <= ts[1] ? ys[1] : x >= ts[end] ? ys[end] :
    (k = searchsortedlast(ts, x); ys[k] + (ys[k+1] - ys[k]) * (x - ts[k]) / (ts[k+1] - ts[k]))

function public_drive(scale)
    p, z = n221204_pulse(), zylstra2022_fig2()
    k = 0.958                                 # N210601's foot to N221204's 164.4 eV at 3.3 ns
    foot(τ) = k * z.dante_n210601[argmin(abs.(z.t .- τ))]
    t = collect(0.0:0.05:10.0)
    T = map(t) do τ
        τ >= 3.3 ? p.dante1[argmin(abs.(p.t .- τ))] : τ >= 1.65 ? foot(τ) :
        τ < 0.2 ? 0.0 : τ < 0.8 ? foot(1.65) * (τ - 0.2) / 0.6 : foot(1.65)
    end
    return t, scale .* T
end

mode = get(ARGS, 1, "table")
if mode == "nova"
    nova = [(0.77, 0.33), (1.07, 0.565), (1.55, 0.78), (1.93, 0.86), (2.05, 1.03), (2.27, 1.14), (3.0, 1.52)]
    TD(t) = 2.6 * max(t, 1e-3)^0.15
    tn = range(0, 2.5; length = 2501)
    println("Nova, 260 eV (t/1 ns)^0.15: measured half signal, HeatFront arrival, 1-D transmission")
    println("onset (1%), bare and backed by 35 µm of epoxy, and the rear gold's density then (g/cm³)")
    println("| Au (µm) | measured | HeatFront | 1-D bare | ρ | 1-D backed | ρ |")
    for (x, measured) in nova
        hr = burnthrough_time(BurnthroughAlbedo(thickness = x), tn, 260 .* tn .^ 0.15)
        b = foil_run(TD, 2.0; mass = 19.3e-4x)
        e = foil_run(TD, 2.0; mass = 19.3e-4x, backing = 1.185 * 35e-4)
        on(r) = (t1 = crossing(r.t, r.tr, 0.01); (t1, r.ρ[something(findfirst(>=(t1), r.t), length(r.t))]))
        @printf("| %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |\n", x, measured, hr, on(b)..., on(e)...)
    end
else
    scale = parse(Float64, get(ARGS, 2, "1.0"))
    td, Td = public_drive(scale)
    TDh(t) = interp(td, Td, t) / 100
    t = collect(0.0:0.002:9.0)
    T = [interp(td, Td, x) for x in t]
    hf = HeatFront()
    # Calibration 1 is the coefficient that puts 2.5 µm through at 7.5 ns on this drive.
    k75 = findlast(<=(7.5 + 1e-9), t)
    rate(x) = (max(x, 0.0) / 100)^(hf.temperature_exponent / hf.time_exponent)
    J = sum((t[k] - t[k-1]) * (rate(T[k]) + rate(T[k-1])) / 2 for k in 2:k75)
    kc = areal_mass(BurnthroughAlbedo(thickness = 2.5)) / (hf.coefficient * J^hf.time_exponent)
    # Calibration 2 is the offset.
    off = 7.5 - burnthrough_time(BurnthroughAlbedo(thickness = 2.5), t, T)
    @printf("drive ×%.2f: %s eV at 2.5, 4.5, 5.5, 6.5, 7.5, 8.5 ns; calibrated coefficient ×%.3f, offset %+.2f ns\n",
            scale, join((@sprintf("%.0f", interp(td, Td, x)) for x in (2.5, 4.5, 5.5, 6.5, 7.5, 8.5)), ", "), kc, off)
    println("| Au (µm) | backing | HeatFront | ×coefficient | +offset | 1-D onset | front + backing | 1-D emergence (A ×1.3–×0.7) |")
    for x in (1.5, 2.0, 2.5, 3.0)
        m = BurnthroughAlbedo(thickness = x, transmission = RosselandTransmission())
        tb = burnthrough_time(m, t, T)
        tk = burnthrough_time(BurnthroughAlbedo(thickness = x, front = HeatFront(coefficient = kc * hf.coefficient)), t, T)
        # the flux through the gold as `BurnthroughAlbedo` has it: s f(T) σT⁴ after the front
        s(τ) = (u = clamp((τ - tb) / m.width, 0, 1); u^2 * (3 - 2u))
        Fm = [s(t[k]) * foil_transmission(m, T[k]) * σ1 * (T[k] / 100)^4 for k in eachindex(t)]
        for (name, ρb, d, A) in BACKINGS
            mb = ρb * d * 0.1                 # mg/cm²
            r = foil_run(TDh, 9.0; mass = 19.3e-4x, backing = mb * 1e-3, ρb)
            @printf("| %.1f | %s | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f (%.2f–%.2f) |\n", x, name, tb, tk, tb + off,
                    crossing(r.t, r.tr, 0.01), backing_burnthrough(t, Fm, mb, A),
                    backing_burnthrough(r.t, r.Fr, mb, A), backing_burnthrough(r.t, r.Fr, mb, 1.3A),
                    backing_burnthrough(r.t, r.Fr, mb, 0.7A))
        end
    end
end
