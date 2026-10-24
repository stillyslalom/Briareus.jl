# The wall re-emission closure of `drive_history` (`wall_release`) against an independent
# reference: 1-D nonlinear radiation diffusion in a static gold wall with the power laws
# of Hammer & Rosen, PoP 10, 1829 (2003) (e = 3.4 T^1.6 ρ^-0.14 MJ/g, κ_R = 7200 T^-1.5
# ρ^0.2 cm²/g, T in heV) and a Marshak boundary driven by the incident flux σT_D(t)^4.
# The closure is given the same T_D(t) and the medium's own constant-T scaling (a static
# medium: E ∝ T^((4+α+β)/2) t^0.5, a = 3.55 and b = 0.5; coefficient from the
# reference's constant-T energy). Prints the release factor each drive needs (the
# closure is linear in it for a prescribed T_D) and the check at `WallLossAlbedo`'s
# default. Research notes, "Wall re-emission and the capsule radius (2026-10-24)".
#
#   julia --project=. examples/wall_reemission.jl [ρ g/cc (19.3)]      # ~10 min

using Briareus, Printf, LinearAlgebra

const σ1 = SIGMA_SB * 1e8 * 1e2 * 1e-3      # σ (1 heV)^4: TW/mm² → TW/cm² → MJ/(ns cm²)
ρ = parse(Float64, get(ARGS, 1, "19.3"))
cv(T) = 3.4 * 1.6 * T^0.6 * ρ^-0.14          # ∂e/∂T, MJ/(g heV)
ee(T) = 3.4 * T^1.6 * ρ^-0.14
κ(T) = 7200 * T^-1.5 * ρ^0.2
const dm = [1e-8 * 1.04^(i - 1) for i in 1:300] # areal-mass zones from the surface (0.032 g/cm²)
const Tcold = 0.02

# One backward-Euler step: Picard in y = T⁴ (κ and ∂e/∂y from the iterate, energy
# conserved); returns the net flux into the wall through the Marshak boundary.
function step!(T, Told, TD, dt)
    n = length(T)
    F0 = 0.0
    for _ in 1:40
        y = T .^ 4
        A = Tridiagonal(zeros(n - 1), zeros(n), zeros(n - 1))
        rhs = zeros(n)
        for i in 1:n
            dedy = cv(T[i]) / (4T[i]^3)
            A.d[i] += dm[i] * dedy / dt
            rhs[i] += dm[i] * (dedy * y[i] - (ee(T[i]) - ee(Told[i]))) / dt
        end
        g0 = 2 / (1 + 3κ(T[1]) * dm[1] / 4)       # Marshak: F = 2(σT_D⁴ − σT_B⁴)
        A.d[1] += σ1 * g0
        rhs[1] += σ1 * g0 * TD^4
        for i in 1:n-1
            km = (κ(T[i]) * dm[i] + κ(T[i+1]) * dm[i+1]) / (dm[i] + dm[i+1])
            c = (4σ1 / 3) / (km * (dm[i] + dm[i+1]) / 2)
            A.d[i] += c; A.du[i] -= c; A.d[i+1] += c; A.dl[i] -= c
        end
        Tn = max.(A \ rhs, 1e-12) .^ 0.25
        δ = maximum(abs.(Tn .- T) ./ Tn)
        T .= Tn
        F0 = σ1 * g0 * (TD^4 - T[1]^4)
        δ < 1e-9 && break
    end
    return F0
end

# The reference's net flux (MJ/(ns cm²)) and stored energy every `out` ns.
function reference(TD, tend; dt = 1e-3, out = 0.05)
    T = fill(Tcold, length(dm))
    ts, Fs, Es = Float64[], Float64[], Float64[]
    every = round(Int, out / dt)
    for k in 1:round(Int, tend / dt)
        F = step!(T, copy(T), TD(k * dt), dt)
        if k % every == 0
            push!(ts, k * dt); push!(Fs, F); push!(Es, sum(dm .* (ee.(T) .- ee(Tcold))))
        end
    end
    return ts, Fs, Es
end

# the medium's constant-T scaling
_, _, E2 = reference(t -> 2.0, 1.0; out = 1.0)
_, _, E3 = reference(t -> 3.0, 1.0; out = 1.0)
a, b = log(E3[1] / E2[1]) / log(1.5), 0.5
K = E3[1] / 3.0^a                                # E = K (T/1 heV)^a t^b
@printf("constant T: a = %.3f (similarity 3.55), E(300 eV, 1 ns) = %.4f MJ/cm²\n", a, E3[1])
wall(release) = WallLossAlbedo(; coefficient = b * K / σ1, temperature_exponent = a - 4,
                               time_exponent = b - 1, range = (-Inf, Inf), t0 = 1e-6, release)

# The closure's net flux along a prescribed T_D history: `drive_history`'s per-element
# bookkeeping (the wall-energy integral at the running maximum, the deficit history).
function closure(m, TD, times; release = true)
    Iw = rate = Tmax = Dw = 0.0
    D, ts, F = Vector{Float64}[], Float64[], Float64[]
    τ0 = nothing
    for τ in times
        T = 100 * TD(τ)
        Δ = τ0 === nothing ? 0.0 : τ - τ0
        h = τ0 === nothing ? Briareus._WallHistory(0.0, τ) :
            !release ? Briareus._WallHistory(Iw + Δ * rate / 2, Δ / 2) :
            Briareus._WallHistory(Iw + Δ * rate / 2, Δ / 2, Tmax, Briareus._deficit_sum(m, D, ts, 1, τ), Dw, Δ)
        α, E, _ = Briareus._wall_split(m, h, T)
        push!(F, ((1 - α) * SIGMA_SB * T^4 - E) * 1e2 * 1e-3)
        Th = release ? max(T, Tmax) : T
        rate = Briareus._loss_rate(m, Th)
        Iw = τ0 === nothing ? τ * rate : h.A + Δ * rate / 2
        if release
            Dw = Briareus._wall_potential(m, Th) - Briareus._wall_potential(m, T)
            Tmax = Th
            push!(D, [Dw]); push!(ts, τ)
        end
        τ0 = τ
    end
    return F
end

∫(t, f) = sum((t[k+1] - t[k]) * (f[k] + f[k+1]) / 2 for k in 1:length(t)-1)
default = WallLossAlbedo().release
grid = collect(0.1:0.1:9.5)
k8 = findfirst(≈(8.0), grid)
# a NIF-like drive: 100 eV foot to 4 ns, 300 eV by 5 ns, flat to 8 ns, then down by
# `drop` over `fall` ns and held
for (drop, fall) in ((0.05, 0.5), (0.10, 0.5), (0.30, 1.0))
    TD(t) = t < 4 ? 1.0 : t < 5 ? 1.0 + 2.0 * (t - 4) : t < 8 ? 3.0 :
            t < 8 + fall ? 3.0 * (1 - drop * (t - 8) / fall) : 3.0 * (1 - drop)
    ts, Fr, _ = reference(TD, 9.5)
    i8 = findfirst(≈(8.0), ts)
    ref = ∫(ts[i8:end], Fr[i8:end])
    F0, F1 = closure(wall(0.0), TD, grid), closure(wall(1.0), TD, grid)
    E0, E1 = ∫(grid[k8:end], F0[k8:end]), ∫(grid[k8:end], F1[k8:end])
    Fd, Fn = closure(wall(default), TD, grid), closure(wall(default), TD, grid; release = false)
    @printf("\nT_D 300 → %.0f eV over %.1f ns from 8 ns: energy into the wall 8–9.5 ns (MJ/cm²)\n",
            300(1 - drop), fall)
    @printf("  reference %+.4f; closure: no release %+.4f, release %.2f %+.4f; the release the reference needs %.3f\n",
            ref, ∫(grid[k8:end], Fn[k8:end]), default, ∫(grid[k8:end], Fd[k8:end]), (E0 - ref) / (E0 - E1))
    println("  t (ns)  F/σT_D⁴: reference / release / no release")
    for tq in (5.0, 7.0, 8.1, 8.25, 8.5, 9.0, 9.5)
        i, j = findfirst(≈(tq), ts), findfirst(≈(tq), grid)
        (i === nothing || j === nothing) && continue
        s = σ1 * TD(tq)^4
        @printf("  %5.2f   %+.3f / %+.3f / %+.3f\n", tq, Fr[i] / s, Fd[j] / s, Fn[j] / s)
    end
end
