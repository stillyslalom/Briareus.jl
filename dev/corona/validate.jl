# Validation of the 3-D ray tracer (dev/corona/Corona.jl) against closed forms in planar
# profiles (the spherical checks, V3–V5 of the research notes, now guard the package's
# operator in test/corona.jl).
#   julia --project=. -t 8 dev/corona/validate.jl
#   V1 planar linear ramp n_e/n_c = x/L: parabolic path, turning depth x_t = L cos²θ,
#      re-exit at z = 2L sin 2θ, absorption 1 - exp(-(32/15)(ν_c L/c) cos⁵θ) (Kruer ch. 5).
#   V2 planar exponential n_e/n_c = e^{x/L}: turning where n_e/n_c = cos²θ, absorption
#      1 - exp(-(8/3)(ν_c L/c) cos³θ) (Kruer; MN17 Eq. 29, the Stage-0 η_a weight).
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

