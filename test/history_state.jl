# `history_state`: a drive history's albedos at a time of the pulse, for a linear
# scenario at that time (a phase), and the state of its thin-foil windows.

using Briareus, Test

@testset "history_state" begin
    kw = (; nz = 12, nphi = 16, capsule_ntheta = 8, capsule_nphi = 12)
    wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8)
    foil = Material(:Au; albedo = BurnthroughAlbedo(; thickness = 1.0), conversion = 0.8)
    t = hohlraum(; kw..., wall, windows = [WallWindow(phi = 0, diameter = 2.0, material = foil)])
    m = RadiationModel(t; closed = true)
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * 4.0)) for b in nif_quads()]
    op = DepositionOperator(t, q)
    # a picket, a trough and a peak
    pulse = PulseShape([0.0, 0.2, 0.6, 0.8, 1.6, 2.0, 3.0], [0.0, 3.0, 3.0, 0.5, 0.5, 8.0, 8.0])
    beams_at(τ) = with_pulses(q, pulse, τ)
    times = collect(0.25:0.25:3.0)
    h = drive_history(m, op, beams_at, times; keep_albedo = true, lmax = 4)
    @test h.transmission isa Vector{Vector{Float64}} && length(h.transmission) == length(h.time)
    @test drive_history(m, op, beams_at, times[1:2]).transmission === nothing
    obs, we, ew = elements(t, :capsule), elements(t, :wall), elements(t, :window_1)
    tb = h.burnthrough[:window_1]
    @test all(isfinite, tb) && times[1] < minimum(tb) && maximum(tb) + 0.4 < times[end]

    # at a step's time: the step's own albedos, and a scenario with them gives the
    # history's capsule flux moments at the history's beams (a linear solve against the
    # history's iterated one)
    for τ in (0.5, 1.25, 2.5)                       # picket, trough, peak
        k = findfirst(==(τ), h.time)
        st = history_state(t, h, τ)
        @test st.time == τ && st.albedo == h.albedo[k] && st.transmission == h.transmission[k]
        s = Scenario(t, q; observe = :capsule, model = m, st.time, st.albedo)
        @test flux_moments(evaluate(s, beams_at(τ)); lmax = 4) ≈ h.moments[:, k] rtol = 1e-9
    end
    # between steps: linear in time
    st = history_state(t, h, 1.375)
    k = findfirst(==(1.25), h.time)
    @test st.albedo ≈ (h.albedo[k] .+ h.albedo[k+1]) ./ 2 rtol = 1e-14
    # `mean`: the wall's area mean, the history's `wall_albedo`; the rest as it was
    sm = history_state(t, h, 2.5; mean = :wall)
    k = findfirst(==(2.5), h.time)
    @test all(≈(h.wall_albedo[k]; rtol = 1e-12), sm.albedo[we])
    rest = setdiff(1:nelements(t), we)
    @test sm.albedo[rest] == h.albedo[k][rest]
    # `remake` takes the state as a new scenario does
    s0 = Scenario(t, q; observe = :capsule, model = m)
    s1 = remake(s0; sm.time, sm.albedo)
    @test s1.albedo == sm.albedo && s1.time == 2.5
    @test s1.moment_response == Scenario(t, q; observe = :capsule, model = m, sm.time, sm.albedo).moment_response

    # the window: intact before its burnthrough, burned through after its transition,
    # where it transmits the foil's fraction 𝒯 of the flux on it
    early, late = history_state(t, h, times[1]), history_state(t, h, times[end])
    @test early.windows[:window_1] == (; transition = 0.0, transmission = 0.0, burned = false)
    @test late.windows[:window_1].transition == 1 && late.windows[:window_1].burned
    @test late.windows[:window_1].transmission ≈ foil_transmission(foil.albedo) rtol = 1e-12
    @test all(≈(foil_transmission(foil.albedo); rtol = 1e-12), late.transmission[ew])
    mid = history_state(t, h, minimum(tb) + 0.2)
    @test 0 < mid.windows[:window_1].transition < 1 && !mid.windows[:window_1].burned

    @test_throws ArgumentError history_state(t, drive_history(m, op, beams_at, times[1:2]), 0.25)
    @test_throws ArgumentError history_state(t, h, 0.1)
    @test_throws ArgumentError history_state(t, h, 3.5)

    # a window whose transmission follows its temperature profile (`DiffusionTransmission`):
    # once burned through, each element transmits the fraction for the flux on it at that step
    fd = Material(:Au; albedo = BurnthroughAlbedo(; thickness = 1.0, transmission = DiffusionTransmission()),
                  conversion = 0.8)
    td = hohlraum(; kw..., wall, windows = [WallWindow(phi = 0, diameter = 2.0, material = fd)])
    hd = drive_history(RadiationModel(td; F = m.F), op, beams_at, times; keep_albedo = true, lmax = 4,
                       record = sol -> radiation_temperature.(sol.H[ew]))
    ld = history_state(td, hd, times[end])
    @test ld.windows[:window_1].burned
    @test ld.transmission[ew] ≈ foil_transmission.(Ref(fd.albedo), hd.record[end]) rtol = 1e-3
    # below the gray estimate at the same temperature (`RosselandTransmission`)
    gray = BurnthroughAlbedo(; thickness = 1.0, transmission = RosselandTransmission())
    @test all(ld.transmission[ew] .< foil_transmission.(Ref(gray), hd.record[end]))
end
