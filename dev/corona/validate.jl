# Validation of the corona prototype (dev/corona/Corona.jl) against closed forms.
#   julia --project=. -t 8 dev/corona/validate.jl
# References (independent of both methods):
#   V1 planar linear ramp n_e/n_c = x/L: parabolic path, turning depth x_t = L cos²θ,
#      re-exit at z = 2L sin 2θ, absorption 1 - exp(-(32/15)(ν_c L/c) cos⁵θ) (Kruer ch. 5).
#   V2 planar exponential n_e/n_c = e^{x/L}: turning where n_e/n_c = cos²θ, absorption
#      1 - exp(-(8/3)(ν_c L/c) cos³θ) (Kruer; MN17 Eq. 29, the Stage-0 η_a weight).
#   V3 spherical n_e/n_c = a/r²: r_t = √(a + b²); swept angle inside radius R,
#      2(b/r_t) arccos(r_t/R); optical depth 2(ν_c/c) a² [√(R²-r_t²)/(2r_t²R²) +
#      arcsec(R/r_t)/(2r_t³)] (elementary integrals; derived here).
#   V4 spherical exponential, L/r_c → 0: tends to V2 with sin ψ = b/r_c; at finite L/r_c
#      this measures how far the Stage-0 η_a weight is from the spherical corona.
#   V5 kernel-table interpolation in b and the tabulated-profile error.
include(joinpath(@__DIR__, "Corona.jl"))
using .Corona, LinearAlgebra, Printf, StaticArrays
const C = Corona
const V3 = SVector{3,Float64}

println("== V1 planar linear ramp (c tracer) ==")
let L = 0.1, nuc = 5.0
    @printf("%6s %10s %10s %10s %10s %12s %12s %9s\n", "θ°", "x_t", "ref", "z_exit", "ref", "absorbed", "ref", "rel err")
    worst = 0.0
    for θd in (0.0, 20.0, 40.0, 60.0, 75.0), h in (0.004, 0.001)
        θ = deg2rad(θd)
        xmax = Ref(0.0)
        x, v, ℓ, _ = C.trace_ray(C.LinearSlab(L, nuc), V3(1e-15, 0, 0), V3(cos(θ), sin(θ), 0), h,
                                 (x, v) -> x[1] < 0 && v[1] < 0) do xm, f
            xmax[] = max(xmax[], xm[1])
        end
        zex = x[3] == 0 ? x[2] - x[1] * v[2] / v[1] : NaN
        ref = 1 - exp(-32 / 15 * nuc * L * cos(θ)^5)
        A = 1 - exp(-ℓ)
        h == 0.001 && (worst = max(worst, abs(A - ref) / ref))
        @printf("%6.1f %10.6f %10.6f %10.6f %10.6f %12.9f %12.9f %9.1e  (hτ=%g)\n", θd, xmax[],
                L * cos(θ)^2, zex, 2L * sin(2θ), A, ref, abs(A - ref) / ref, h)
    end
    @printf("V1 worst relative absorption error at hτ = 0.001: %.2e\n", worst)
end

println("\n== V2 planar exponential (c tracer) ==")
let L = 0.05, nuc = 7.5   # (8/3) ν L/c = 1.0
    x0 = -25L
    worst = 0.0
    for θd in (0.0, 30.0, 50.0, 70.0, 80.0), h in (0.002, 0.0005)
        θ = deg2rad(θd)
        n0 = sqrt(1 - exp(x0 / L))
        xmax = Ref(-Inf)
        x, v, ℓ, _ = C.trace_ray(C.ExpSlab(L, nuc), V3(x0, 0, 0), n0 * V3(cos(θ), sin(θ), 0), h,
                                 (x, v) -> x[1] < x0 && v[1] < 0) do xm, f
            xmax[] = max(xmax[], xm[1])
        end
        ref = 1 - exp(-8 / 3 * nuc * L * cos(θ)^3)
        A = 1 - exp(-ℓ)
        h == 0.0005 && (worst = max(worst, abs(A - ref) / ref))
        @printf("θ=%5.1f° hτ=%-7g x_t=%10.6f (ref %10.6f)  A=%.9f (ref %.9f) rel %.1e\n", θd, h,
                xmax[], L * log(cos(θ)^2), A, ref, abs(A - ref) / ref)
    end
    @printf("V2 worst relative absorption error at hτ = 0.0005: %.2e\n", worst)
end

println("\n== V3 spherical a/r² profile: (b) radial quadrature and (c) tracer vs closed form ==")
let a = 0.04, nuc = 2.0, R = 2.0
    p = C.InvSqProfile(a, nuc, R)
    wb = (0.0, 0.0, 0.0); wc = (0.0, 0.0, 0.0)
    for b in (0.0, 0.1, 0.3, 0.6, 1.0, 1.5)
        rt = sqrt(a + b^2)
        Θ = 2 * (b / rt) * acos(rt / R)
        τ = 2nuc * a^2 * (sqrt(R^2 - rt^2) / (2rt^2 * R^2) + acos(rt / R) / (2rt^3))
        k = C.ray_kernel(p, b; N = 400)
        Θb, τb = 2sum(k.A), 2sum(k.T)
        # (c): start on r = R inside the medium, Bouguer-consistent direction
        nR = sqrt(1 - a / R^2)
        sφ = b / (R * nR)
        v0 = nR * V3(-sqrt(1 - sφ^2), sφ, 0)
        x, v, ℓ, _ = C.trace_ray((xm, f) -> nothing, p, V3(R, 0, 0), v0, 0.0005,
                                 (x, v) -> dot(x, x) > R^2 && dot(x, v) > 0)
        δτ = (norm(x) - R) / dot(x / norm(x), v)        # overshoot past r = R
        Θc = mod(atan(x[2], x[1]), 2π) - b / R^2 * δτ
        b == 0 && (Θc = 0.0)
        ℓc = ℓ - C.g(p, R) * δτ
        wb = max.(wb, (abs(k.rt - rt), abs(Θb - Θ), abs(τb - τ) / τ))
        wc = max.(wc, (0.0, abs(Θc - Θ), abs(ℓc - τ) / τ))
        @printf("b=%.2f r_t=%.12f (b: %.12f)  Θ=%.9f b:%.9f c:%.9f  τ=%.9f b:%.9f c:%.9f\n",
                b, rt, k.rt, Θ, Θb, Θc, τ, τb, ℓc)
    end
    @printf("V3 (b) max |Δr_t| %.1e mm, |ΔΘ| %.1e rad, rel Δτ %.1e\n", wb...)
    @printf("V3 (c) max |ΔΘ| %.1e rad, rel Δτ %.1e (hτ = 5e-4)\n", wc[2], wc[3])
end

println("\n== V4 spherical exponential → planar cos³ (Stage-0 η_a) as L/r_c → 0 ==")
let rc = 0.4
    D = 2.4      # (8/3) ν_c L / c, η⊥ = 1 - e^{-2.4} = 0.909
    bs = (0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95)
    @printf("%8s", "L/r_c"); foreach(s -> @printf("  %8s", "b/rc=$(s)"), bs); println()
    @printf("%8s", "η_a(ψ)"); foreach(s -> @printf("  %8.5f", 1 - exp(-D * (1 - s^2)^1.5)), bs); println()
    for Lr in (0.0025, 0.01, 0.04, 0.1, 0.2)
        L = Lr * rc
        p = C.ExpProfile(rc, L, D / (8 / 3 * L), rc + 25L)
        @printf("%8.4f", Lr)
        for s in bs
            k = C.ray_kernel(p, s * rc; N = 800)
            @printf("  %8.5f", 1 - exp(-sum(2k.T)))
        end
        println()
    end
    # the (b) kernel and the (c) tracer on the same rays (L/r_c = 0.2): kernel-level agreement
    L = 0.08; p = C.ExpProfile(rc, L, D / (8 / 3 * L), rc + 14L)
    R = p.Rout
    worstA, worstΘ = 0.0, 0.0
    for s in (0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 2.0)
        b = s * rc
        k = C.ray_kernel(p, b; N = 400)
        _, _, esc, αout = C.kernel_nodes(k)
        xe = V3(-sqrt(R^2 - b^2), b, 0)
        x, v, ℓ, _ = C.trace_ray((xm, f) -> nothing, p, xe, V3(1, 0, 0), 0.001, C.sphere_done(R))
        # exit angle of position from -k (= -x̂), back-corrected for the overshoot
        δτ = (norm(x) - R) / dot(x / norm(x), v)
        αc = mod(atan(x[2], -x[1]), 2π) - b / R^2 * δτ
        Ac = isinf(ℓ) ? 1.0 : 1 - exp(-ℓ)
        worstA = max(worstA, abs(Ac - (1 - esc))); worstΘ = max(worstΘ, abs(αc - αout))
        @printf("L/rc=0.2 b/rc=%.1f: absorbed b %.7f c %.7f; exit angle b %.6f c %.6f rad\n",
                s, 1 - esc, Ac, αout, αc)
    end
    @printf("V4 (b) vs (c) per ray: max |ΔA| %.1e, max |Δα_exit| %.1e rad\n", worstA, worstΘ)
end

println("\n== V5 kernel table and tabulated profile ==")
let rc = 0.4, L = 0.08
    p = C.ExpProfile(rc, L, 2.4 / (8 / 3 * L), rc + 14L)
    for Nb in (256, 1024, 4096)
        K = C.KernelTable(p; Nb, N = 200)
        bs = [(i + 0.5) * (K.b[2] - K.b[1]) for i in 0:Nb-2]
        err = maximum(abs(C.absorbed(K, b) - (1 - C.kernel_nodes(C.ray_kernel(p, b; N = 200))[3])) for b in bs[1:max(1, Nb ÷ 256):end])
        @printf("Nb=%5d: max interpolation error of the absorbed fraction (mid-cells) %.2e\n", Nb, err)
    end
    for N in (50, 100, 200, 400)
        e = maximum(abs(sum(C.ray_kernel(p, b; N).T) - sum(C.ray_kernel(p, b; N = 3200).T)) / sum(C.ray_kernel(p, b; N = 3200).T) for b in (0.0, 0.2, 0.4, 0.6))
        @printf("N=%4d panels: max relative depth error %.2e\n", N, e)
    end
    # an exponential is exact in the table (ln n_e linear); use a non-trivial profile:
    # scale length growing outward and T_e varying (a stand-in for a 1-D hydro snapshot)
    Lr(r) = 0.06 * (r / rc)^2
    fne(r) = exp(-(r - rc) / Lr(r))
    Te(r) = 1 + 0.3 * tanh((r - rc) / 0.2)
    fg(r) = 15.0 * fne(r)^2 / Te(r)^1.5
    pf = C.FuncProfile(fne, fg, 1.6)
    for dr in (0.02, 0.01, 0.005, 0.0025)
        n = round(Int, (pf.Rout - 0.2) / dr) + 1
        tp = C.TableProfile(r -> C.ne(pf, r), r -> C.g(pf, r), 0.2, pf.Rout, n)
        rs = range(0.3, pf.Rout - 0.01; length = 2001)
        e1 = maximum(abs(C.ne(tp, r) / C.ne(pf, r) - 1) for r in rs)
        e2 = maximum(abs((1 - C.kernel_nodes(C.ray_kernel(tp, b))[3]) - (1 - C.kernel_nodes(C.ray_kernel(pf, b))[3])) for b in 0:0.05:0.9)
        @printf("table Δr=%5.3f mm: max rel n_e error %.1e, max absorbed-fraction error %.1e\n", dr, e1, e2)
    end
end
