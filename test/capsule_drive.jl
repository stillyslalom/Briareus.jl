# CapsuleDrive: the drive handed to a capsule-only hydrodynamics run (ROADMAP 17).

@testset "capsule drive" begin
    t = hohlraum(nz = 10, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12,
                 wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
    m = RadiationModel(t)
    q = nif_quads()
    pulse = PulseShape([0.0, 1.0, 5.0, 6.0], [0.0, 2.0, 8.0, 0.0])
    times = [0.5, 2.0, 4.0, 5.5]
    # every beam crosses the axis at z = ±4 and lands toward the waist (equator-hot), or
    # at z = ±1 and lands on the far cap and entrance-hole end (pole-hot)
    run(z) = (qz = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * z)) for b in q];
              op = DepositionOperator(t, qz);
              (drive_history(m, op, τ -> with_pulses(qz, pulse, τ), times; lmax = 8), op, qz))
    h, op, qz = run(4.0)
    d = capsule_drive(h, m; provenance = (; pulse = "test ramp", hohlraum = "coarse NIF"))
    # shapes: a uniform grid of ≤ 50 ps over the solved times, l = 0…8
    @test d.time[1] == 0.5 && d.time[end] == 5.5
    @test maximum(diff(d.time)) <= 0.05 + 1e-12 && length(d.time) == 101
    @test size(d.modes) == (9, 101) && all(==(1), d.modes[1, :])
    @test length(d.Tr) == length(d.flux) == length(d.absorbed) == 101
    @test d.radius == 0.9 && d.area ≈ sum(t.elem_area[elements(t, :capsule)])
    @test d.provenance[:pulse] == "test ramp" && d.provenance[:code] == "Briareus"
    @test d.provenance[:observe] == [:capsule] && haskey(d.provenance, :git)
    # round trip: on the history's own times the drive is the history
    dh = capsule_drive(h, m; times = h.time)
    @test dh.flux == h.flux && dh.Tr ≈ h.Tr
    @test dh.modes ≈ h.moments ./ h.moments[1:1, :]
    @test dh.absorbed == h.absorbed[:capsule]
    @test d(2.0).flux ≈ h.flux[2] && d(2.0).modes ≈ h.moments[:, 2] ./ h.moments[1, 2]
    @test d(3.0).flux ≈ (h.flux[2] + h.flux[3]) / 2           # linear between steps
    # resampling the resampled drive at the history times recovers it (piecewise linear)
    @test capsule_drive(h, m; times = [2.0, 4.0]).flux ≈ h.flux[2:3]
    @test capsule_drive(h, m; dt = 1.0).time == collect(0.5:1.0:5.5)
    @test size(capsule_drive(h, m; lmax = 4).modes, 1) == 5
    @test_throws ArgumentError capsule_drive(h, m; lmax = 10)
    @test_throws ArgumentError capsule_drive(h, m; times = [0.0, 1.0])
    @test_throws ArgumentError capsule_drive(h, m; observe = :wall)   # not one sphere
    # a₀ is the capsule flux energy_balance books: (1 − α) a₀ A = absorbed (constant albedo)
    α = albedos(t)[elements(t, :capsule)[1]]
    @test all(k -> isapprox(dh.absorbed[k], (1 - α) * dh.flux[k] * d.area; rtol = 2e-3), 1:4)
    sol = solve(m, deposit(op, with_pulses(qz, pulse, 4.0)); time = 4.0)
    @test energy_balance(sol).absorbed[:capsule] ≈ (1 - α) * d(4.0).flux * d.area rtol = 2e-3
    # sign: P2 > 0 is pole-hot. A field that is larger toward the poles…
    es = elements(t, :capsule)
    H = zeros(nelements(t))
    for e in es
        c = t.elem_centroid[e]
        H[e] = 1 + 0.5 * (c[3] / 0.9)^2
    end
    a = legendre_moments(H, t, es; lmax = 8)
    hs = (; time = [0.0, 1.0], moments = hcat(a, a), absorbed = Dict{Symbol,Vector{Float64}}())
    ds = capsule_drive(hs, t)
    @test ds.modes[3, 1] > 0.1 && all(isnan, ds.absorbed)
    # …and the beams: landing on the caps is pole-hot, toward the waist equator-hot
    @test all(d.modes[3, :] .< 0)
    hp, _, _ = run(1.0)
    @test all(capsule_drive(hp, m).modes[3, :] .> 0)
    @test occursin("CapsuleDrive(101 times", sprint(show, d))
end
