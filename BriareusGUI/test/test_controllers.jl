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

@testset "diagnostic windows (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    t0 = ts.target[]
    n0 = nelements(t0)
    @test add_window!(ts)                                   # closed, 1 × 1 mm at φ = 90°
    w1 = only(ts.windows[])
    @test w1.phi == 90 && w1.z == 0 && !w1.open && w1.albedo == ts.params[][:wall_albedo]
    @test :window_1 in surface_names(ts.target[]) && window_name(1) === :window_1
    @test occursin("with 1 window", target_summary(ts))
    @test add_window!(ts) && ts.windows[][2].phi == 270     # the next free default place
    t = ts.target[]

    # the geometry's rejections surface as messages and keep the target
    @test !edit_window!(ts, 2, :phi, "90.2")
    @test occursin("overlap", ts.error[]) && ts.target[] === t && ts.windows[][2].phi == 270
    @test !edit_window!(ts, 1, :z, "5.0")                   # past the end of the 10.1 mm wall
    @test occursin("past the end", ts.error[]) && ts.windows[][1].z == 0
    @test startswith(ts.error[], "window 1")                  # which window, for the user
    @test !edit_window!(ts, 1, :width, "0") && occursin("positive", ts.error[])
    @test !edit_window!(ts, 1, :albedo, "1.5") && occursin("between 0 and 1", ts.error[])
    @test !edit_window!(ts, 1, :height, "tall") && occursin("must be a number", ts.error[])
    @test !add_window!(ts; phi = 92.0) && occursin("overlap", ts.error[])
    @test length(ts.windows[]) == 2 && ts.target[] === t
    @test_throws ArgumentError edit_window!(ts, 1, :depth, "1")
    @test_throws ArgumentError edit_window!(ts, 3, :z, "1")
    @test_throws ArgumentError remove_window!(ts, 3)

    # accepted edits rebuild; a closed window re-emits with its own material
    @test edit_window!(ts, 1, :z, "1.5") && ts.error[] == "" && ts.target[] !== t
    @test edit_window!(ts, 1, :albedo, "0.3") && edit_window!(ts, 1, :conversion, "0.1")
    tw = ts.target[]
    c1 = tw.components[component_index(tw, :window_1)]
    @test c1.kind === :surface && albedo(c1.material) ≈ 0.3 && c1.material.conversion ≈ 0.1
    zs = [tw.elem_centroid[e][3] for e in elements(tw, :window_1)]
    @test all(z -> 1.0 < z < 2.0, zs)                       # the mesh conforms to the window
    # an open window is a hole (an aperture)
    @test set_window_open!(ts, 2, true)
    @test ts.target[].components[component_index(ts.target[], :window_2)].kind === :aperture
    @test !(:window_2 in surface_names(ts.target[]))
    # a dimension edit that leaves a window off the wall is rejected
    @test !edit_target!(ts, :length, "3.0") && occursin("past the end", ts.error[])
    @test ts.params[][:length] == 10.1
    # resolution and kind keep the windows; a foil has none
    set_resolution!(ts, :standard)
    @test :window_1 in component_names(ts.target[]) && nelements(ts.target[]) > 4n0
    set_resolution!(ts, :coarse)
    set_target_kind!(ts, :foil)
    @test !any(n -> startswith(String(n), "window"), component_names(ts.target[]))
    @test !occursin("window", target_summary(ts))
    @test !add_window!(ts) && occursin("hohlraum", ts.error[])
    set_target_kind!(ts, :hohlraum)
    @test component_names(ts.target[])[end-1:end] == [:window_1, :window_2]
    # removing renumbers: the open window becomes window_1
    remove_window!(ts, 1)
    @test length(ts.windows[]) == 1 && ts.windows[][1].open
    @test ts.target[].components[component_index(ts.target[], :window_1)].kind === :aperture
    remove_window!(ts, 1)
    @test isempty(ts.windows[]) && nelements(ts.target[]) == n0
    # windows given up front
    ts2 = TargetSetup(resolution = :coarse, windows = [WindowSpec(phi = 0.0, height = 2.0, open = true)])
    @test :window_1 in component_names(ts2.target[])
end

@testset "VisRad interchange through the beam table (no GL)" begin
    bt = BeamTable(preset_beams(:nif_quads); preset = :nif_quads)
    edit_beams!(bt, 1, :z, "3.3")
    edit_beams!(bt, 2, :defocus, "0.5")
    set_enabled!(bt, 4, false)                                 # the 50° cone off
    ref, en = copy(bt.beams[]), copy(bt.enabled[])
    dir = mktempdir()

    # SST pointing: one row per NIF beam of the enabled quads
    sst = export_visrad_sst(bt, joinpath(dir, "sst.csv"))
    @test bt.path[] == sst
    @test countlines(sst) == 1 + 4 * count(en)
    fresh() = BeamTable(preset_beams(:nif_quads); preset = :nif_quads)
    b2 = fresh()
    import_visrad_sst!(b2, sst)
    @test all(((a, b),) -> a.aim ≈ b.aim && a.defocus ≈ b.defocus,
              zip(bt.beams[][en], b2.beams[][en]))
    @test beam_rows(b2)[1].z ≈ 3.3 && beam_rows(b2)[2].defocus ≈ 0.5
    # the quads absent from the file come back disabled, at their power
    @test b2.enabled[] == en && all(b -> b.power == 8.0, b2.beams[])
    @test b2.preset[] === :nif_quads && b2.path[] == sst
    @test occursin("disabled", b2.warnings[][1])
    # the reader's warnings are collected: a beam the table does not have
    open(io -> println(io, "B999,1.0,2.0,3.0,0.0"), sst, "a")
    b2b = fresh()
    import_visrad_sst!(b2b, sst)
    @test any(w -> occursin("ignored", w) && occursin("B999", w), b2b.warnings[])

    # Laser Beam Summary: pointing, power and spot round trip; off beams disabled
    summ = export_visrad_summary(bt, joinpath(dir, "summary.csv"))
    b3 = fresh()
    edit_beams!(b3, 1, :power, "5")                       # the file restores it
    import_visrad_summary!(b3, summ)
    @test b3.enabled[] == en
    @test all(((a, b),) -> a.aim ≈ b.aim && a.defocus ≈ b.defocus && a.power ≈ b.power &&
                           a.spot.a ≈ b.spot.a && a.spot.b ≈ b.spot.b,
              zip(ref, b3.beams[]))
    # a file that is not VisRad's leaves the table unchanged
    bad = joinpath(dir, "bad.csv")
    write(bad, "Beam,X\nB115,1\n")
    b4 = fresh()
    before = b4.beams[]
    @test_throws ArgumentError import_visrad_sst!(b4, bad)
    @test b4.beams[] === before
    @test_throws ArgumentError import_visrad_sst!(BeamTable(), sst)       # nothing to apply to
    # loading a preset clears the warnings
    load_preset!(b2, :nif_quads)
    @test isempty(b2.warnings[])
end
