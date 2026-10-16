# The Target and Beams steps' controllers: no scenario, no GL.

@testset "TargetSetup (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    t0 = ts.target[]
    @test component_names(t0) == [:wall, :cap_top, :cap_bottom, :leh_top, :leh_bottom, :capsule]
    @test surface_names(t0) == [:wall, :cap_top, :cap_bottom, :capsule]
    @test occursin("hohlraum", target_summary(ts)) && occursin("$(nelements(t0)) elements", target_summary(ts))

    # an accepted edit rebuilds; a rejected one keeps the target and reports why
    @test edit_target!(ts, :length, "9.4")
    @test ts.target[] !== t0 && ts.params[][:length] == 9.4 && ts.error[] == ""
    t1 = ts.target[]
    @test !edit_target!(ts, :leh_radius, "3.5")               # wider than the hohlraum
    @test ts.target[] === t1 && ts.params[][:leh_radius] == 1.55 && occursin("leh_radius", ts.error[])
    @test !edit_target!(ts, :wall_albedo, "1.2")
    @test occursin("between 0 and 1", ts.error[])
    @test !edit_target!(ts, :radius, "abc")
    @test occursin("must be a number", ts.error[])
    @test edit_target!(ts, :radius, "2.875") && ts.error[] == ""      # unchanged value
    @test_throws ArgumentError edit_target!(ts, :height, "1")

    # no capsule
    @test edit_target!(ts, :capsule_radius, "0")
    @test !(:capsule in component_names(ts.target[]))
    edit_target!(ts, :capsule_radius, "0.9")

    # resolution and kind
    n = nelements(ts.target[])
    set_resolution!(ts, :standard)
    @test nelements(ts.target[]) > 4n
    set_resolution!(ts, :coarse)
    @test_throws ArgumentError set_resolution!(ts, :ultra)
    set_target_kind!(ts, :foil)
    @test component_names(ts.target[]) == [FOIL_NAME, :main]
    @test foil_center(ts) == Briareus.Vec3(6.0, 0, 0)
    @test !edit_target!(ts, :distance, "1.2")                 # inside the sphere
    @test edit_target!(ts, :main_radius, "0")
    @test component_names(ts.target[]) == [FOIL_NAME]
    @test_throws ArgumentError set_target_kind!(ts, :sphere)
    set_target_kind!(ts, :hohlraum)
    @test :capsule in component_names(ts.target[])            # the hohlraum edits were kept
    @test ts.params[][:length] == 9.4
end

@testset "BeamTable (no GL)" begin
    bt = BeamTable(preset_beams(:nif_quads); preset = :nif_quads)
    @test length(bt.beams[]) == 48 && all(bt.enabled[])
    rows = beam_rows(bt)
    @test [r.key for r in rows] == ["inner23", "inner30", "outer44", "outer50"]
    @test all(r -> r.count == 8 || r.count == 16, rows) && sum(r.count for r in rows) == 48
    r23 = rows[1]
    @test r23.z ≈ 3.09 && r23.r ≈ -1.09 && r23.power ≈ 8.0 && r23.enabled === true
    # published pointing, mirrored for the lower hemisphere
    @test all(b -> abs(b.aim[3]) ≈ 3.09, bt.beams[][r23.members])
    @test any(b -> b.aim[3] < 0, bt.beams[][r23.members])
    @test occursin("48 of 48 beams on", beams_summary(bt)) && occursin("384 TW", beams_summary(bt))

    # a group edit sets one field on every member and nothing else
    before = copy(bt.beams[])
    @test edit_beams!(bt, 1, :z, "3.3")
    after = bt.beams[]
    m = r23.members
    @test all(b -> abs(b.aim[3]) ≈ 3.3, after[m])
    @test [b.power for b in after] == [b.power for b in before]
    others = setdiff(1:48, m)
    @test after[others] == before[others]
    @test beam_rows(bt)[1].r ≈ -1.09 atol = 1e-12                   # r kept
    @test !edit_beams!(bt, 1, :power, "-2")
    @test occursin("must not be negative", bt.error[]) && bt.beams[] === after
    @test !edit_beams!(bt, 1, :z, "far")
    @test edit_beams!(bt, 1, :z, "3.3") && bt.error[] == ""
    @test_throws ArgumentError edit_beams!(bt, 1, :spot, "1")

    # per quad and per beam
    set_grouping!(bt, :quad)
    @test length(beam_rows(bt)) == 48
    q = beam_rows(bt)[5]
    @test edit_beams!(bt, 5, :defocus, "1.5")
    @test bt.beams[][only(q.members)].defocus == 1.5
    set_grouping!(bt, :beam)
    @test beam_rows(bt)[5].z ≈ bt.beams[][5].aim[3]                  # not mirrored per beam
    @test_throws ArgumentError set_grouping!(bt, :bundle)
    set_grouping!(bt, :group)

    # on/off: disabled beams keep their place and deliver no power
    set_enabled!(bt, 2, false)
    @test beam_rows(bt)[2].enabled === false
    @test length(active_beams(bt)) == 48 - beam_rows(bt)[2].count
    @test length(display_beams(bt)) == 48
    @test all(b -> b.power == 0, display_beams(bt)[beam_rows(bt)[2].members])
    m2 = beam_rows(bt)[2].members
    set_grouping!(bt, :quad)                                          # one quad back on
    set_enabled!(bt, findfirst(r -> r.members == [m2[1]], beam_rows(bt)), true)
    set_grouping!(bt, :group)
    @test beam_rows(bt)[2].enabled === missing
    set_enabled!(bt, 2, false)

    # CSV round trip
    path = joinpath(mktempdir(), "beams.csv")
    export_beams(bt, path)
    @test bt.path[] == path
    bt2 = BeamTable()
    @test beams_summary(bt2) == "no beams"
    import_beams!(bt2, path)
    @test bt2.preset[] === :custom && length(bt2.beams[]) == 48
    @test all(((a, b),) -> a.port ≈ b.port && a.aim ≈ b.aim && a.power ≈ b.power &&
                           a.defocus ≈ b.defocus && a.name == b.name,
              zip(bt.beams[], bt2.beams[]))
    export_beams(bt, path; enabled_only = true)
    @test length(read_beam_table(path)) == length(active_beams(bt))

    # presets
    for p in BEAM_PRESETS
        @test !isempty(preset_beams(p))
    end
    @test length(preset_beams(:nif_beams)) == 192
    om = preset_beams(:omega_hohlraum; params = TargetSetup(resolution = :coarse).params[])
    @test length(om) == 40 && all(b -> abs(b.aim[3]) ≈ 10.1 / 2, om)
    foilp = TargetSetup(kind = :foil, resolution = :coarse).params[]
    bl = preset_beams(:omega_backlighter; target_kind = :foil, params = foilp)
    @test all(b -> b.aim ≈ foil_center(foilp) && b.port[1] < -cosd(70), bl)
    @test_throws ArgumentError preset_beams(:lmj)
    set_grouping!(bt, :quad)
    load_preset!(bt, :omega_hohlraum)
    @test bt.preset[] === :omega_hohlraum && length(bt.enabled[]) == 40 && all(bt.enabled[])
    @test bt.by[] === :group                         # OMEGA beams have no quad names
    @test_throws ArgumentError set_grouping!(bt, :quad)
    @test_throws ArgumentError Controllers._check_grouping(bt.beams[], :quad)   # the optimizer's check
end
