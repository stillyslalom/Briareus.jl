# Retuning a drive to a reference: pulses from a cone fraction, pulse tables, and the
# per-time bounded least-squares retune against a full drive history.

@testset "retune to a reference drive" begin
    t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                 wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8),
                 nz = 16, nphi = 24, capsule_ntheta = 8, capsule_nphi = 16)
    rz = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20), :inner30 => (0.60, 3.63),
              :inner23 => (1.09, 3.09))
    quads = map(nif_quads()) do b
        r, z = rz[b.group]
        φ = port_angles(b)[2]
        remake(b; aim = (-r * cosd(φ), -r * sind(φ), sign(b.port[3]) * z))
    end
    # pulses from the total power and the cone fraction
    tt = [0.0, 1.0, 3.0, 4.0]
    total = [0.0, 100.0, 300.0, 0.0]
    cf = [0.3, 0.3, 0.25, 0.25]
    pul = cone_fraction_pulses(quads, tt, total, cf)
    inner = count(b -> b.group in (:inner23, :inner30), quads)
    @test pul[:inner23](3.0) ≈ 0.25 * 300 / inner                  # at the given times
    @test sum(pul[b.group](2.0) for b in quads) ≈ 200
    @test_throws ArgumentError cone_fraction_pulses(quads, tt, total, [0.3, 0.3, 1.2, 0.2])
    # pulse tables round trip
    f = joinpath(mktempdir(), "pulses.csv")
    write_pulse_table(f, pul)
    back = read_pulse_table(f)
    @test Set(keys(back)) == Set(keys(pul))
    @test all(back[k](τ) ≈ pul[k](τ) for k in keys(pul) for τ in (0.5, 2.0, 3.5))
    # keyed pulses
    qp = Dict(Symbol(b.quad) => PulseShape([0.0, 1.0], [1.0, 2.0]) for b in quads)
    @test with_pulses(quads, qp, 0.5; by = :quad)[1].power ≈ 1.5

    # the reference shot, then a dropped bundle, then the retune
    m = RadiationModel(t; closed = true)
    op = DepositionOperator(t, quads)
    times = 0.5:0.5:3.5
    ref = drive_history(m, op, τ -> with_pulses(quads, pul, τ), times; lmax = 4, ylm = true,
                        keep_albedo = true)
    @test length(ref.albedo) == length(ref.time) && length(ref.albedo[1]) == nelements(t)
    keep = findall(b -> nif_bundle(b) != "31", quads)
    kept = quads[keep]
    opk = op[keep]
    amp1(h) = [mode_amplitude(h.ylm[:, k], 1) for k in eachindex(h.time)]
    dropped = drive_history(m, opk, τ -> with_pulses(kept, pul, τ), ref.time; lmax = 4, ylm = true)
    cap(k, τ) = (0.0, 1.5 * pul[first(b.group for b in kept if Symbol(b.quad) == k)](τ) + 1e-9)
    rt = retune_drive(m, opk, kept, ref; by = :quad, limits = cap)
    @test size(rt.power) == (length(kept), length(ref.time))
    @test all(rt.power[g, k] <= cap(rt.keys[g], rt.time[k])[2] + 1e-9 for g in eachindex(rt.keys)
              for k in eachindex(rt.time))
    @test all(>=(-1e-12), rt.power)
    @test maximum(abs.(rt.drive .- 1)) < 1e-3
    tuned = drive_history(m, opk, τ -> with_pulses(kept, rt.pulses, τ; by = :quad), ref.time;
                          lmax = 4, ylm = true)
    # the per-time linear model is the history's own: same albedos to the solve's tolerance
    @test maximum(abs.(rt.shape .- tuned.ylm[1:25, :] ./ tuned.ylm[1:1, :])) < 1e-4
    # and the dropped bundle's mode 1 is gone, the drive restored
    @test maximum(amp1(dropped)) > 0.01
    @test maximum(amp1(tuned)) < maximum(amp1(dropped)) / 10
    @test maximum(abs.(tuned.Tr .- ref.Tr)) < 0.5
    # staying near nominal powers: a smaller departure, a larger residual
    nominal(k, τ) = pul[first(b.group for b in kept if Symbol(b.quad) == k)](τ)
    near = retune_drive(m, opk, kept, ref; by = :quad, limits = cap, nominal, stay = 1e-2)
    dep(r) = maximum(abs(r.power[g, k] / nominal(r.keys[g], r.time[k]) - 1)
                     for g in eachindex(r.keys) for k in eachindex(r.time))
    @test dep(near) < dep(rt)
    @test sum(near.residual) > sum(rt.residual)
    # fixed powers (a quad held at zero) are honored
    q0 = rt.keys[1]
    fixed = retune_drive(m, opk, kept, ref; by = :quad,
                         limits = (k, τ) -> k == q0 ? (0.0, 0.0) : cap(k, τ))
    @test all(iszero, fixed.power[1, :])
    @test_throws ArgumentError retune_drive(m, opk, kept, dropped)            # no albedos
    @test_throws ArgumentError retune_drive(m, opk, kept, ref; drive_weight = 0)
    @test_throws ArgumentError retune_drive(m, opk, kept, ref; stay = 1.0)
end
