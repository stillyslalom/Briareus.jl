# Can the gold wall hold the drive after the laser falls? An independent reference for
# `drive_history`'s wall closure (`wall_history`, `wall_release`): 1-D planar Lagrangian
# radiation hydrodynamics of an ablating gold wall with the power laws of Hammer & Rosen,
# PoP 10, 1829 (2003) (e = 3.4 T^1.6 ρ^-0.14 MJ/g, κ_R = 7200 T^-1.5 ρ^0.2 cm²/g, T in heV)
# and P = r e/V, r = 0.25 (Rosen, "Fundamentals of ICF hohlraums", 2005, OSTI 885380,
# Lecture 1 eqs. 1.26, 1.30), flux-limited radiation diffusion and a Marshak boundary.
# `hydro = false` freezes the wall (the static reference of `examples/wall_reemission.jl`).
# Research notes, "Late-time drive: the wall, the plasma and Dante (2026-10-25)";
# reference/litreview/late_time_drive.md.
#
#   julia --project=. examples/late_time_wall.jl release            # ~3 min
#   julia --project=. examples/late_time_wall.jl hohlraum [shot]    # ~15 min per shot
#
# `release`: the release factor (`WallLossAlbedo`'s `release`) the ablating wall needs
# under prescribed drive drops, the closure given the ablating wall's own constant-T
# scaling (as `wall_reemission.jl` does for the static wall).
# `hohlraum`: a 0-D hohlraum on the shot's measured pulse (`n191118`, `n210601`: Zylstra
# 2022 Fig. 2; `n221204`: Pak 2024 Fig. 1), ηP = A_w F_wall(T) + A_open σT⁴ solved for T
# every step, with the static wall, the ablating wall and `drive_history`'s closure
# (Lindl's gold loss, the default release). The hohlraum is HYBRID-E's 6.4 × 11.24 mm,
# η = 0.75, the capsule a 0.7-absorbing sphere of its initial radius (`capsule=0` drops
# it and closes the LEHs to 0.857 r₀, the most favourable case for the wall). Dante is
# printed beside it for the shape of the fall only: nothing here is fitted to it.

using Briareus, Printf, LinearAlgebra

const σ1 = SIGMA_SB * 1e8 * 1e2 * 1e-3       # σ (1 heV)^4: TW/mm² → MJ/(ns cm²)
const r_eos = 0.25
const c2 = 1e-5                               # (cm/ns)² per MJ/g
eT(T, V) = 3.4 * T^1.6 * V^0.14               # MJ/g (ρ^-0.14 = V^0.14)
Te(e, V) = (max(e, 1e-12) / (3.4 * V^0.14))^(1 / 1.6)
κR(T, V) = 7200 * T^-1.5 * V^-0.2
cvT(T, V) = 3.4 * 1.6 * T^0.6 * V^0.14

# zones by areal mass from the surface (zone 1), nodes 1 (surface, vacuum) to n+1 (rigid)
mutable struct Wall
    dm::Vector{Float64}
    x::Vector{Float64}
    u::Vector{Float64}
    V::Vector{Float64}
    T::Vector{Float64}
    e::Vector{Float64}
    hydro::Bool
end
function Wall(; n = 220, dm0 = 2e-6, q = 1.04, ρ0 = 19.3, T0 = 0.02, hydro = true)
    dm = [dm0 * q^(i - 1) for i in 1:n]
    V = fill(1 / ρ0, n)
    T = fill(T0, n)
    return Wall(dm, [0.0; cumsum(dm .* V)], zeros(n + 1), V, T, eT.(T, V), hydro)
end
Base.copy(w::Wall) = Wall(copy(w.dm), copy(w.x), copy(w.u), copy(w.V), copy(w.T), copy(w.e), w.hydro)
# internal + kinetic energy above the cold wall (MJ/cm²)
energy(w) = sum(w.dm .* (w.e .- eT(0.02, 1 / 19.3))) +
            sum(((w.u[i] + w.u[i+1]) / 2)^2 / (2c2) * w.dm[i] for i in eachindex(w.dm))

function hydro_dt(w)
    w.hydro || return Inf
    return minimum(0.25 * (w.x[i+1] - w.x[i]) /
                   (sqrt((1 + r_eos) * r_eos * w.e[i] * c2) + abs(w.u[i+1] - w.u[i]) + 1e-30)
                   for i in eachindex(w.dm))
end

# explicit Lagrangian step (quadratic viscosity), P = 0 outside the surface
function hydro!(w, dt)
    w.hydro || return
    n = length(w.dm)
    Pt = [r_eos * w.e[i] / w.V[i] +
          (w.u[i+1] < w.u[i] ? 2.0 * (w.u[i+1] - w.u[i])^2 / w.V[i] / c2 : 0.0) for i in 1:n]
    for j in 1:n
        mj = j == 1 ? w.dm[1] / 2 : (w.dm[j-1] + w.dm[j]) / 2
        w.u[j] -= c2 * (Pt[j] - (j == 1 ? 0.0 : Pt[j-1])) / mj * dt
    end
    w.x .+= w.u .* dt
    for i in 1:n
        Vn = (w.x[i+1] - w.x[i]) / w.dm[i]
        w.e[i] -= Pt[i] * (Vn - w.V[i])
        w.V[i] = Vn
        w.T[i] = Te(w.e[i], Vn)
    end
end

# backward-Euler diffusion in y = T⁴ (Picard; sum flux limiter), Marshak boundary at the
# incident T_D; returns the net flux into the wall (MJ/(ns cm²))
function rad!(w, TD, dt; tol = 1e-9)
    n = length(w.dm)
    T, V, dm = w.T, w.V, w.dm
    eold = copy(w.e)
    F0 = 0.0
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
        for i in 1:n-1
            τ = (κR(T[i], V[i]) * dm[i] + κR(T[i+1], V[i+1]) * dm[i+1]) / 2
            c = (4σ1 / 3) / (τ + (4 / 3) * abs(y[i+1] - y[i]) / ((y[i] + y[i+1]) / 2))
            A.d[i] += c; A.du[i] -= c; A.d[i+1] += c; A.dl[i] -= c
        end
        Tn = max.(A \ rhs, 1e-12) .^ 0.25
        δ = maximum(abs.(Tn .- T) ./ Tn)
        T .= Tn
        F0 = σ1 * g0 * (TD^4 - T[1]^4)
        δ < tol && break
    end
    w.e .= eT.(T, V)
    return F0
end

# monotone scalar root (Illinois)
function root(g, a, b; tol = 1e-7)
    fa, fb = g(a), g(b)
    while fa * fb > 0
        b *= 1.5; fb = g(b)
    end
    side = 0
    for _ in 1:100
        c = (a * fb - b * fa) / (fb - fa)
        fc = g(c)
        abs(b - a) < tol * c && return c
        if fc * fb > 0
            b, fb = c, fc; side == -1 && (fa /= 2); side = -1
        else
            a, fa = c, fc; side == 1 && (fb /= 2); side = 1
        end
    end
    return (a + b) / 2
end

# prescribed incident T_D(t) (heV): times, net flux, stored energy every `out` ns
function driven(TD, tend; hydro = true, dtmax = 2e-3, out = 0.05)
    w = Wall(; hydro)
    t, nxt = 0.0, out
    ts, Fs, Es = Float64[], Float64[], Float64[]
    while t < tend - 1e-12
        dt = min(dtmax, hydro_dt(w), tend - t)
        hydro!(w, dt)
        F = rad!(w, TD(t + dt), dt)
        t += dt
        t >= nxt - 1e-9 && (push!(ts, t); push!(Fs, F); push!(Es, energy(w)); nxt += out)
    end
    return ts, Fs, Es
end

# `drive_history`'s per-element wall bookkeeping (as in `wall_reemission.jl`): the net
# flux into the wall (MJ/(ns cm²)) at T (eV), and the state update after the step
mutable struct Closure
    m::WallLossAlbedo
    Iw::Float64; rate::Float64; Tmax::Float64; Dw::Float64
    D::Vector{Vector{Float64}}; ts::Vector{Float64}; τ0::Float64
end
Closure(m) = Closure(m, 0.0, 0.0, 0.0, 0.0, Vector{Float64}[], Float64[], NaN)
function history(c::Closure, τ)
    isnan(c.τ0) && return Briareus._WallHistory(0.0, τ)
    Δ = τ - c.τ0
    return Briareus._WallHistory(c.Iw + Δ * c.rate / 2, Δ / 2, c.Tmax,
                                 Briareus._deficit_sum(c.m, c.D, c.ts, 1, τ), c.Dw, Δ)
end
flux(c::Closure, h, T) = (s = Briareus._wall_split(c.m, h, T); ((1 - s[1]) * SIGMA_SB * T^4 - s[2]) * 1e2 * 1e-3)
function advance!(c::Closure, h, τ, T)
    Th = max(T, c.Tmax)
    c.rate = Briareus._loss_rate(c.m, Th)
    c.Iw = isnan(c.τ0) ? τ * c.rate : h.A + (τ - c.τ0) * c.rate / 2
    c.Dw = Briareus._wall_potential(c.m, Th) - Briareus._wall_potential(c.m, T)
    c.Tmax = Th
    push!(c.D, [c.Dw]); push!(c.ts, τ); c.τ0 = τ
end

∫(t, f) = sum((t[k+1] - t[k]) * (f[k] + f[k+1]) / 2 for k in 1:length(t)-1)

mode = get(ARGS, 1, "release")
if mode == "release"
    # the ablating wall's constant-T scaling E = K T^a t^b (1–2 ns)
    ts0, _, E3 = driven(t -> 3.0, 2.0)
    _, _, E2 = driven(t -> 2.0, 2.0)
    k1 = findfirst(t -> t >= 1 - 1e-9, ts0)
    b = log(E3[end] / E3[k1]) / log(2)
    a = log(E3[end] / E2[end]) / log(1.5)
    K = E3[end] / (3.0^a * 2.0^b)
    @printf("ablating wall at constant T: E = %.4f T^%.3f t^%.3f MJ/cm² (Lindl: 5.2e-3 T^3.3 t^0.62; HR k=0: 5.8e-3 T^3.34 t^0.59)\n", K, a, b)
    wall(f) = WallLossAlbedo(; coefficient = b * K / σ1, temperature_exponent = a - 4,
                             time_exponent = b - 1, range = (-Inf, Inf), t0 = 1e-6, release = f)
    grid = collect(0.05:0.05:9.5)
    k8 = findfirst(≈(8.0), grid)
    closure_flux(f, TD) = (c = Closure(wall(f));
                           [(h = history(c, τ); F = flux(c, h, 100TD(τ)); advance!(c, h, τ, 100TD(τ)); F) for τ in grid])
    println("100 eV foot to 4 ns, 300 eV by 5 ns, flat to 8 ns, then down: energy into the wall 8–9.5 ns (MJ/cm²)")
    for (drop, fall) in ((0.10, 0.5), (0.30, 1.0), (0.50, 0.5))
        TD(t) = t < 4 ? 1.0 : t < 5 ? 1.0 + 2.0 * (t - 4) : t < 8 ? 3.0 :
                t < 8 + fall ? 3.0 * (1 - drop * (t - 8) / fall) : 3.0 * (1 - drop)
        ts, Fr, Es = driven(TD, 9.5)
        i8 = findfirst(t -> t >= 8 - 1e-9, ts)
        ref = ∫(ts[i8:end], Fr[i8:end])
        E0, E1, Ed = (∫(grid[k8:end], closure_flux(f, TD)[k8:end]) for f in (0.0, 1.0, WallLossAlbedo().release))
        @printf("  %2.0f %% over %.1f ns: stored at 8 ns %.4f; reference %+.4f; closure, release 0 %+.4f, %.2f %+.4f; release the reference needs %.3f\n",
                100drop, fall, Es[i8], ref, E0, WallLossAlbedo().release, Ed, (E0 - ref) / (E0 - E1))
    end
else
    shot = get(ARGS, 2, "n191118")
    capsule = get(Dict(String.(split(x, '=', limit = 2)) for x in ARGS if occursin('=', x)), "capsule", "1")
    Lh, R = 11.24, 3.2
    r = shot == "n191118" ? 1.82 : 1.55
    rc, rfac = capsule == "0" ? (0.0, 0.857) : (shot == "n221204" ? 1.135 : 1.1285, 1.0)
    Aw = (2π * R * Lh + 2π * (R^2 - r^2)) * 1e-2                  # cm²
    Ao = (2π * (rfac * r)^2 + 0.7 * 4π * rc^2) * 1e-2
    d = shot == "n221204" ? n221204_pulse() : zylstra2022_fig2()
    tp = d.t
    pp = shot == "n221204" ? d.power : getfield(d, Symbol("power_", shot))
    dd = shot == "n221204" ? d.dante1 : getfield(d, Symbol("dante_", shot))
    P(t) = t <= tp[1] ? pp[1] : t >= tp[end] ? 0.0 : (k = searchsortedlast(tp, t);
           pp[k] + (pp[k+1] - pp[k]) * (t - tp[k]) / (tp[k+1] - tp[k]))
    S(t) = 0.75 * P(t) * 1e-3                                        # MJ/ns
    tend = 10.0
    function reference(hydro)
        w = Wall(; hydro)
        t, nxt, TD = 0.0, 0.05, 0.3
        ts, Ts, Es = Float64[], Float64[], Float64[]
        while t < tend - 1e-12
            dt = min(2e-3, hydro_dt(w), tend - t)
            hydro!(w, dt)
            s = S(t + dt)
            TD = root(x -> Aw * rad!(copy(w), x, dt; tol = 1e-8) + Ao * σ1 * x^4 - s,
                      max(TD * 0.7, 0.03), TD * 1.2 + 0.1)
            rad!(w, TD, dt)
            t += dt
            t >= nxt - 1e-9 && (push!(ts, t); push!(Ts, 100TD); push!(Es, energy(w)); nxt += 0.05)
        end
        return ts, Ts, Es
    end
    function closure(m; dt = 0.05)
        c = Closure(m)
        ts, Ts, T = Float64[], Float64[], 30.0
        for τ in dt:dt:tend
            h = history(c, τ)
            T = root(x -> Aw * flux(c, h, x) + Ao * σ1 * (x / 100)^4 - S(τ), 1.0, max(1.3T, 50.0))
            advance!(c, h, τ, T)
            push!(ts, τ); push!(Ts, T)
        end
        return ts, Ts
    end
    tsd, Tsd, _ = reference(false)
    ta, Ta, Ea = reference(true)
    tc, Tc = closure(WallLossAlbedo())
    at(ts, Ts, t) = Ts[argmin(abs.(ts .- t))]
    @printf("%s: wall %.3f cm², open %.3f cm² (capsule %s), η 0.75\n", shot, Aw, Ao, capsule == "0" ? "none, LEH 0.857 r₀" : "at R₀")
    println("  t (ns)  P (TW)   T_D (eV): static / ablating / closure    Dante   ablating wall energy (kJ)")
    for t in 5.0:0.25:10.0
        @printf("  %5.2f  %6.1f   %6.1f / %6.1f / %6.1f   %6.1f   %5.0f\n", t, P(t), at(tsd, Tsd, t),
                at(ta, Ta, t), at(tc, Tc, t), dd[argmin(abs.(tp .- t))], 1e3 * Aw * at(ta, Ea, t))
    end
end
