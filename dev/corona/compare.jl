# (b) tabulated impact-parameter kernels vs (c) 3-D ray tracing on an OMEGA 40-beam PDD
# case (public inputs; reference/litreview/polar_direct_drive.md §4.1):
#   rings 21.4°/42.0°/58.9° of `omega_beams` (5, 5, 10 beams per hemisphere), shifted
#   Δr = 90/150/150 µm normal to the beam toward the PDD equator (Mn14 Table 1, Ma16's
#   "spherical-shell" pointing), SG4 spots exp[-(r/r₀)^n], r₀ = 308 µm, n = 3.66 (Ma16),
#   f/6.7, equal powers.
# Corona (an engineering snapshot, not from a simulation): isothermal exponential,
#   n_e/n_c = exp(-(r - r_c)/L), r_c = 0.40 mm, L = 0.08 mm, ν_c/c = 11.25 mm⁻¹ (planar
#   normal-incidence absorption 1 - e^{-2.4} = 0.909; ≈ T_e 3.6 keV, Z 3.5, lnΛ 8 at
#   351 nm), truncated at r_c + 14 L.
#
#   julia --project=. -t 8 dev/corona/compare.jl <stage> [args]
# stages: check (rotation and harmonic-vs-binned checks), b (build, evaluate, gradients
#         of (b) and the hard-sphere models), c <hlat> <hτ> (one (c) run),
#         cpart <hlat> <hτ> <k> <K> then cmerge <hlat> <hτ> <K> (a (c) run split over
#         processes, beams k:K:40), cgrad <hlat> <hτ> (ForwardDiff through (c)),
#         p2 <ε> <hlat> <hτ> ((c) on the P2-perturbed corona), ylm <hlat> (grid bias of (c)).
# The numbers are in reference/research-notes.md, "Corona model prototype".
# Results are appended to ENV["CORONA_OUT"] (a Serialization file).
include(joinpath(@__DIR__, "Corona.jl"))
using .Corona, Briareus, LinearAlgebra, Printf, StaticArrays, ForwardDiff, Serialization, Random
const C = Corona

const OUT = get(ENV, "CORONA_OUT", joinpath(tempdir(), "corona_results.jls"))
const LMAX = 30
const rc, Lsc = 0.40, 0.08
const nuc = 2.4 / (8 / 3 * Lsc)
const ETA_PERP = 1 - exp(-2.4)
const prof = C.ExpProfile(rc, Lsc, nuc, rc + 14Lsc)
const Rout = prof.Rout
const NCONE = 7
const SPOT = SuperGaussianSpot(0.308; order = 3.66)
const RING = Dict(:cone1 => 1, :cone2 => 2, :cone3 => 3)
const BASE = omega_beams(; cones = (:cone1, :cone2, :cone3), spot = SPOT, power = 1.0, fnumber = 6.7)

# toward the PDD equator, normal to the beam, azimuth unchanged
shiftdir(b) = (eθ = Briareus.port_basis(b.port)[1]; sign(b.port[3]) * eθ)
pdd_beams(Δr) = [remake(b; aim = Δr[RING[b.group]] * shiftdir(b)) for b in BASE]
const DR0 = [0.090, 0.150, 0.150]

function record(key, val)
    d = isfile(OUT) ? deserialize(OUT) : Dict{String,Any}()
    d[key] = val
    serialize(OUT, d)
end

summ(c) = (; σl = C.sigma_l(c, LMAX), σ = C.sigma_rms(c, LMAX),
           a20 = c[Briareus.mode_index(2, 0)] / c[1], a40 = c[Briareus.mode_index(4, 0)] / c[1])

function show_spec(name, s, absorbed, incident)
    @printf("%-30s absorbed %.5f  σ(l≤30) %.4f %%  c20/c00 %+.5f  c40/c00 %+.5f\n", name,
            absorbed / incident, 100s.σ, s.a20, s.a40)
    print("   σ_l (%), l=1..12: ")
    foreach(l -> @printf("%.4f ", 100s.σl[l+1]), 1:12); println()
end

opsphere(n) = Target(Component(:corona, icosphere_mesh(Rout; n), Material(:CH); surface = SphereSurface(Rout)))
hs_model(b) = b < rc ? ([asin(b / rc)], [1.0]) : (Float64[], Float64[])
eta_model(b) = b < rc ? ([asin(b / rc)], [1 - (1 - ETA_PERP)^((1 - (b / rc)^2)^1.5)]) : (Float64[], Float64[])

function main(stage)

if stage == "check"
    # 1. rotation blocks: Y(Q x') = D Y(x') at random points
    Random.seed!(1)
    Q = Matrix(qr(randn(3, 3)).Q); det(Q) < 0 && (Q[:, 1] .*= -1)
    D = C.rotation_blocks(Q, LMAX)
    Y = zeros(LMAX + 1, 2LMAX + 1); Y2 = similar(Y)
    err = 0.0
    for _ in 1:50
        x = normalize(randn(3))
        Briareus._real_ylm!(Y, LMAX, x[3], atan(x[2], x[1]))
        xl = Q * x
        Briareus._real_ylm!(Y2, LMAX, xl[3], atan(xl[2], xl[1]))
        for l in 0:LMAX
            # field = δ at x' in the frame has frame coefficients Y_lm(x'); in the lab Y_lm(Qx')
            err = max(err, maximum(abs, D[l+1] * Y[l+1, LMAX+1-l:LMAX+1+l] - Y2[l+1, LMAX+1-l:LMAX+1+l]))
        end
    end
    @printf("rotation blocks l ≤ 30: max |D Y(x') - Y(Qx')| = %.1e\n", err)
    # 2. harmonic operator vs per-ray exact harmonics on the hard sphere (no binning in
    #    either: every pair's single deposit point summed with its exact Y_lm)
    beams = pdd_beams(DR0)
    op = DepositionOperator(opsphere(12), beams; ncone = NCONE)
    K = C.KernelTable(prof; Nb = 2048, N = 100)
    H = C.HarmonicOperator(op, K.b, C.kappa_table(K.b, hs_model, LMAX), LMAX)
    ch, inc = C.apply(H, beams)
    CO = C.CoronaOperator(op, NCONE, C.AngGrid(16, 2), hs_model)    # just for the pairs
    cref = zeros((LMAX + 1)^2); Pw = zeros(LMAX + 1, LMAX + 1)
    for bi in eachindex(beams)
        w = C.pair_weights(CO, bi, beams[bi])
        for (k, (q, j)) in enumerate(zip(CO.pair_q[bi], CO.pair_j[bi]))
            x = C.V3(op.pts[bi][q]); d = C.V3(op.dirs[bi][j])
            bv = x - dot(x, d) * d; b = norm(bv)
            b < rc || continue
            hit = -sqrt(rc^2 - b^2) * d + bv
            C.accumulate_ylm!(cref, Y, Pw, LMAX, hit, w[k])
        end
    end
    # the harmonic operator interpolates the node angle α(b) linearly in b (Nb = 2048)
    @printf("harmonic operator vs exact per-ray harmonics (hard sphere): max |Δc_lm|/c00 = %.1e; σ %.6f vs %.6f %%\n",
            maximum(abs.(ch .- cref)) / cref[1], 100C.sigma_rms(ch, LMAX), 100C.sigma_rms(cref, LMAX))
elseif stage == "b"
    beams = pdd_beams(DR0)
    n = parse(Int, get(ARGS, 2, "20"))
    Nb = parse(Int, get(ARGS, 3, "2048"))
    DepositionOperator(opsphere(4), beams[1:1]; ncone = NCONE)       # compile
    t_op = @elapsed op = DepositionOperator(opsphere(n), beams; ncone = NCONE)
    hs = deposit(op, beams)
    @printf("operator on the Rout = %.2f mm sphere (icosphere n = %d): %d points, %.2f s, %.1f MB; intercepted %.6f of %.1f\n",
            Rout, n, sum(length, op.w), t_op, Base.summarysize(op) / 2^20, sum(hs.per_beam), hs.delivered)
    C.kappa_table(C.KernelTable(prof; Nb = 16, N = 100), LMAX)
    t_K = @elapsed K = C.KernelTable(prof; Nb, N = 100)
    t_κ = @elapsed κb = C.kappa_table(K, LMAX)
    @printf("kernel table Nb = %d (2N = 200 nodes) %.2f s; κ table %.2f s, %.1f MB\n", Nb, t_K, t_κ, sizeof(κb) / 2^20)
    res = Dict{String,Any}("t_op" => t_op, "t_K" => t_K + t_κ, "n" => n, "Nb" => Nb,
                           "opbytes" => Base.summarysize(op), "points" => sum(length, op.w))
    for (name, κ) in (("hard sphere r_c, all absorbed", C.kappa_table(K.b, hs_model, LMAX)),
                      ("hard sphere r_c, eta_a (MN17)", C.kappa_table(K.b, eta_model, LMAX)),
                      ("(b) kernel table", κb))
        t_H = @elapsed H = C.HarmonicOperator(op, K.b, κ, LMAX)
        c, inc = C.apply(H, beams)
        t_ev = minimum(@elapsed(C.apply(H, beams)) for _ in 1:3)
        s = summ(c)
        absd = c[1] * sqrt(4π)
        show_spec(name, s, absd, inc)
        @printf("   map %.2f s; %d pairs, %.1f MB (+ operator); evaluate %.3f s\n", t_H, C.npairs(H), C.bytes(H) / 2^20, t_ev)
        r = Dict{String,Any}("σl" => s.σl, "σ" => s.σ, "a20" => s.a20, "a40" => s.a40, "absorbed" => absd / inc,
                             "t_H" => t_H, "t_ev" => t_ev, "bytes" => C.bytes(H), "pairs" => C.npairs(H))
        if name == "(b) kernel table"
            fb(Δr) = C.sigma_rms(C.apply(H, pdd_beams(Δr))[1], LMAX)
            ForwardDiff.gradient(fb, DR0)
            t_g3 = minimum(@elapsed(ForwardDiff.gradient(fb, DR0)) for _ in 1:2)
            gr = ForwardDiff.gradient(fb, DR0)
            h = 1e-5
            fd = [(fb(DR0 .+ h .* (1:3 .== i)) - fb(DR0 .- h .* (1:3 .== i))) / 2h for i in 1:3]
            @printf("   dσ/dΔr (3 rings): AD %s  FD %s  max rel %.1e; gradient %.2f s\n", string(round.(gr; sigdigits = 7)),
                    string(round.(fd; sigdigits = 7)), maximum(abs.(gr .- fd) ./ abs.(fd)), t_g3)
            U = [Briareus.port_basis(b.port) for b in beams]
            beams80(x) = [remake(b; aim = b.aim + x[2i-1] * U[i][1] + x[2i] * U[i][2]) for (i, b) in enumerate(beams)]
            f80(x) = C.sigma_rms(C.apply(H, beams80(x))[1], LMAX)
            ForwardDiff.gradient(f80, zeros(80))
            t_g80 = @elapsed ForwardDiff.gradient(f80, zeros(80))
            @printf("   gradient in 80 aim coordinates: %.2f s (%.1f evaluations)\n", t_g80, t_g80 / t_ev)
            r["gr"] = gr; r["fd"] = fd; r["t_g3"] = t_g3; r["t_g80"] = t_g80
        end
        res[name] = r
    end
    record("b_n$(n)_Nb$(Nb)", res)
elseif stage == "c"
    hlat = parse(Float64, get(ARGS, 2, "0.02"))
    hτ = parse(Float64, get(ARGS, 3, "0.002"))
    nθ = parse(Int, get(ARGS, 4, "720"))
    G = C.AngGrid(nθ, LMAX)
    beams = pdd_beams(DR0)
    C.trace_beams(prof, beams[1:1], G; ncone = 1, hlat = 0.2, hτ)   # compile
    t = @elapsed (E, inc, nr, _) = C.trace_beams(prof, beams, G; ncone = NCONE, hlat, hτ)
    s = summ(C.ylm_coeffs(E, G))
    show_spec("(c) hlat=$hlat hτ=$hτ", s, sum(E), inc)
    @printf("   %d live rays, %.1f s (%d threads), incident %.6f of 40\n", nr, t, Threads.nthreads(), inc)
    record("c_$(hlat)_$(hτ)_$(nθ)", Dict("σl" => s.σl, "σ" => s.σ, "a20" => s.a20, "a40" => s.a40,
                                    "absorbed" => sum(E) / inc, "t" => t, "nr" => nr, "inc" => inc))
elseif stage == "cpart"
    # a subset of the beams (k:K:40), absorbed bins saved for "cmerge" (E is linear in beams)
    hlat = parse(Float64, get(ARGS, 2, "0.01"))
    hτ = parse(Float64, get(ARGS, 3, "0.004"))
    k = parse(Int, get(ARGS, 4, "1")); K = parse(Int, get(ARGS, 5, "2"))
    G = C.AngGrid(720, LMAX)
    beams = pdd_beams(DR0)[k:K:end]
    C.trace_beams(prof, beams[1:1], G; ncone = 1, hlat = 0.2, hτ)
    t = @elapsed (E, inc, nr, _) = C.trace_beams(prof, beams, G; ncone = NCONE, hlat, hτ)
    @printf("part %d/%d: %d rays, %.1f s\n", k, K, nr, t)
    serialize(OUT * ".part_$(hlat)_$(k)_$(K)", (; E, inc, nr, t))
elseif stage == "cmerge"
    hlat = parse(Float64, get(ARGS, 2, "0.01"))
    hτ = parse(Float64, get(ARGS, 3, "0.004"))
    K = parse(Int, get(ARGS, 4, "2"))
    parts = [deserialize(OUT * ".part_$(hlat)_$(k)_$(K)") for k in 1:K]
    G = C.AngGrid(720, LMAX)
    E = sum(p.E for p in parts); inc = sum(p.inc for p in parts)
    nr = sum(p.nr for p in parts); t = sum(p.t for p in parts)
    s = summ(C.ylm_coeffs(E, G))
    show_spec("(c) hlat=$hlat hτ=$hτ", s, sum(E), inc)
    @printf("   %d live rays, %.1f s in %d parts (8 threads), incident %.6f of 40\n", nr, t, K, inc)
    record("c_$(hlat)_$(hτ)_720", Dict("σl" => s.σl, "σ" => s.σ, "a20" => s.a20, "a40" => s.a40,
                                      "absorbed" => sum(E) / inc, "t" => t, "nr" => nr, "inc" => inc))
elseif stage == "cgrad"
    hlat = parse(Float64, get(ARGS, 2, "0.04"))
    hτ = parse(Float64, get(ARGS, 3, "0.004"))
    G = C.AngGrid(720, LMAX)
    fc(Δr) = C.sigma_rms(C.ylm_coeffs(C.trace_beams(prof, pdd_beams(Δr), G; ncone = NCONE, hlat, hτ)[1], G), LMAX)
    t0 = @elapsed fc(DR0)
    t1 = @elapsed fc(DR0)
    tg = @elapsed gr = ForwardDiff.gradient(fc, DR0)
    @printf("(c) hlat=%g hτ=%g: value %.2f s (first %.2f s), gradient in 3 ring shifts %.2f s; dσ/dΔr = %s\n",
            hlat, hτ, t1, t0, tg, string(round.(gr; sigdigits = 6)))
    record("cgrad_$(hlat)_$(hτ)", Dict("t1" => t1, "tg" => tg, "gr" => gr))
elseif stage == "p2"
    ε = parse(Float64, get(ARGS, 2, "0.04"))
    hlat = parse(Float64, get(ARGS, 3, "0.02"))
    hτ = parse(Float64, get(ARGS, 4, "0.002"))
    G = C.AngGrid(720, LMAX)
    beams = pdd_beams(DR0)
    pp = C.P2ExpProfile(rc, Lsc, nuc, Rout, ε)
    C.trace_beams(pp, beams[1:1], G; ncone = 1, hlat = 0.2, hτ)
    t = @elapsed (E, inc, nr, _) = C.trace_beams(pp, beams, G; ncone = NCONE, hlat, hτ)
    s = summ(C.ylm_coeffs(E, G))
    show_spec("(c) P2 ε=$ε hlat=$hlat", s, sum(E), inc)
    @printf("   %d rays, %.1f s\n", nr, t)
    record("p2_$(ε)_$(hlat)_$(hτ)", Dict("σl" => s.σl, "σ" => s.σ, "a20" => s.a20, "a40" => s.a40,
                                        "absorbed" => sum(E) / inc, "t" => t))
elseif stage == "ylm"
    hlat = parse(Float64, get(ARGS, 2, "0.08"))
    beams = pdd_beams(DR0)
    for nθ in (360, 720)
        G = C.AngGrid(nθ, LMAX)
        E, inc, nr, cdirect = C.trace_beams(prof, beams, G; ncone = NCONE, hlat, hτ = 0.004, direct_ylm = true)
        s = summ(C.ylm_coeffs(E, G))
        σd = C.sigma_l(cdirect, LMAX)
        @printf("grid %d×%d vs exact per-deposit Ylm (%d rays): max_l |Δσ_l| = %.2e; σ %.5f vs %.5f %%; rel at l = 10, 20, 30: %.1e %.1e %.1e\n",
                nθ, 2nθ, nr, maximum(abs.(s.σl[2:end] .- σd[2:end])), 100s.σ, 100C.sigma_rms(cdirect, LMAX),
                (abs(s.σl[l+1] - σd[l+1]) / σd[l+1] for l in (10, 20, 30))...)
    end
end
end

main(get(ARGS, 1, "b"))
