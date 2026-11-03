# Diagnostic definitions in the window: the list of definitions found (built-in, a
# folder's, the session's), several recorded per drive history, the folders preference,
# and the definitions saved with a session.

# A definition file in `dir`.
function _defn(dir, id, body)
    write(joinpath(dir, id * ".toml"), "schema = 1\nid = \"$id\"\n" * body)
    return joinpath(dir, id * ".toml")
end

@testset "diagnostic definitions: list, histories, folders, sessions (no GL)" begin
    root = mktempdir()
    withenv("BRIAREUS_CONFIG" => joinpath(root, "config.toml"), "BRIAREUS_DIAGNOSTICS" => nothing) do
        sess = DesignSession(resolution = :coarse)
        tm = sess.time
        @test tm.diagnostics[] == DEFAULT_DIAGNOSTICS
        rows = diagnostic_catalog(sess)
        byid = Dict(r.id => r for r in rows)
        @test byid["nif_dante1"].origin === :builtin && byid["nif_dante1"].usable && isempty(byid["nif_dante1"].changed)
        @test !byid["nif_dante2"].usable && occursin("channel", byid["nif_dante2"].problem)
        @test all(r -> r.origin === :builtin, rows) && issorted([r.id for r in rows])
        @test isempty(diagnostic_folders())
        # the picker's groups: every row once, by kind in the fixed order, the catalog's
        # order within a group, a heading on each group's first row
        gs = Controllers._diagnostic_groups(rows)
        @test sort([g.row.id for g in gs]) == sort([r.id for r in rows])
        @test unique([g.group for g in gs]) == ["X-ray power", "X-ray imagers", "Radiography"]
        @test [g.group for g in gs if g.head] == unique([g.group for g in gs])
        @test all(g -> Controllers._diagnostic_group(g.row.kind) == g.group, gs)
        @test [g.row.id for g in gs if g.row.kind === :xray_power] == [r.id for r in rows if r.kind === :xray_power]
        @test Controllers._diagnostic_group(:neutron_yield) == "Neutron yield"
        # the image view's first limits: the pixels above 1e-3 of the peak plus 15 % (at
        # least 4 pixels) each side, none for an image filling the field or an empty one
        xs = range(-10.0, 10.0; length = 201)
        img = zeros(201, 201)
        img[101:121, 91:101] .= 1.0                               # x 0..2, y -1..0
        img[150, 150] = 1e-4                                      # below the threshold
        @test all(Controllers._signal_window(xs, xs, img) .≈ (-0.4, 2.4, -1.4, 0.4))
        @test Controllers._signal_window(xs, xs, ones(201, 201)) === nothing
        @test Controllers._signal_window(xs, xs, zeros(201, 201)) === nothing

        # a group's folder: Dante 1 with a 12 mm field of view (Tr over the lower hole), a
        # measured channel from a CSV, and a view through a hole the target lacks
        dir = mkpath(joinpath(root, "group"))
        write(joinpath(dir, "ch19.csv"), "energy_eV,response_V_per_GW\n50,0.0\n1000,0.2\n20000,0.0\n")
        _defn(dir, "group_dante1_fov", """
            extends = "nif_dante1"
            name = "Dante 1, 12 mm field"
            public = false
            [view]
            field_of_view_mm = 12
            through = "leh_bottom"
            [[channel]]
            number = 19
            energy_eV = 1000
            solid_angle_sr = 6.5e-7
            response_csv = "ch19.csv"
            """)
        _defn(dir, "group_blind", "extends = \"nif_dante1\"\n[view]\nthrough = \"no_such_hole\"\n")
        @test add_diagnostic_folder!(sess, dir)
        @test occursin("2 definitions", tm.status[])
        @test diagnostic_folders() == [(; path = normpath(dir), saved = true)]
        r = Dict(r.id => r for r in tm.catalog[])["group_dante1_fov"]
        @test r.origin === :user && r.usable && !r.public && length(r.chain) == 3
        @test !add_diagnostic_folder!(sess, joinpath(root, "nope"))

        # several diagnostics in one history, in the order ticked
        set_diagnostic!(sess, "group_dante1_fov")
        set_diagnostic!(sess, "group_blind")
        set_diagnostic!(sess, "nif_dante1", false)
        set_diagnostic!(sess, "nif_dante1")
        set_diagnostic!(sess, "nif_sxi_lower")                     # an imager, time-integrated
        @test tm.diagnostics[] == ["group_dante1_fov", "group_blind", "nif_dante1", "nif_sxi_lower"]
        use_pulse_shape!(sess; step = 2.0)
        start_history!(sess)
        h = tm.history[]
        @test [d.id for d in h.diagnostics] == tm.diagnostics[]
        fov, blind, d1, sxi = h.diagnostics
        @test all(isnan, blind.Tr) && blind.note == "the target has no no_such_hole to look through"
        @test isempty(fov.note) && isempty(d1.note) && all(isfinite, d1.Tr)
        # independent of the definition path: a 12 mm field about TCC and Dante 1's 60 mm
        # field over the whole target see the same emission (the walls radiate inward only)
        k = argmax(d1.intensity)
        @test fov.intensity[k] ≈ d1.intensity[k] rtol = 2e-2
        @test fov.Tr[k] ≈ d1.Tr[k] rtol = 1e-2
        ch = time_chart(sess, :diagnostics)
        @test isempty(ch.empty) && length(ch.left.series) == 2 && length(ch.right.series) == 3
        # the imager: its integrated image is the time integral of the detected intensity
        # (the blur keeps the total of a target inside the field)
        @test sxi.kind === :xray_imager && isempty(sxi.note) && all(isnan, sxi.Tr) && maximum(sxi.intensity) > 0
        a = (sxi.view.camera.width / sxi.view.camera.nx)^2
        @test sum(sxi.frames.images[1]) * a ≈
              sum((h.time[k+1] - h.time[k]) * (sxi.intensity[k] + sxi.intensity[k+1]) / 2
                  for k in 1:length(h.time)-1) rtol = 1e-6
        @test [r.frame for r in diagnostic_images(sess)] == [0, 1]
        show_diagnostic_image!(sess, "nif_sxi_lower", 1)
        im1 = sess.drive.image[]
        @test sess.drive.view[] === :image && im1.unit == "J/cm²/sr (detected)" && size(im1.values) == (1024, 1024)
        @test sum(im1.values) ≈ 1e5 * sum(sxi.frames.images[1]) rtol = 1e-12
        # it opens on its region with signal, inside the field
        w = im1.window
        @test w !== nothing && first(im1.x) - step(im1.x) <= w[1] < w[2] <= last(im1.x) + step(im1.x) &&
              first(im1.y) - step(im1.y) <= w[3] < w[4] <= last(im1.y) + step(im1.y)
        set_time_step!(sess, argmax(sxi.intensity))
        show_diagnostic_image!(sess, "nif_sxi_lower")
        @test sess.drive.image[].unit == "W/cm²/sr (detected)"
        @test sum(sess.drive.image[].values) * a ≈ to_Wcm2(maximum(sxi.intensity)) rtol = 1e-5
        @test [s.label for s in time_chart(sess, :drive).left.series] ==
              ["Capsule Tr", "Dante 1, 12 mm field", "NIF Dante 1"]

        # saved with the session: the definitions and the CSV they read
        file = joinpath(root, "s.briar")
        save_session(file, sess)
        @test occursin("from outside the built-ins: group_blind, group_dante1_fov", sess.status[]) &&
              occursin("not public: group_dante1_fov", sess.status[])
        names = zip_names(ZipReader(read(file)))
        h19 = bytes2hex(open(sha256, joinpath(dir, "ch19.csv")))
        @test "diagnostics/$h19.csv" in names
        state = TOML.parse(zip_readentry(ZipReader(read(file)), "session.toml", String))
        @test state["time"]["diagnostics"] == tm.diagnostics[]
        @test Set(keys(state["diagnostic_definitions"])) == Set(tm.diagnostics[])
        @test state["diagnostic_definitions"]["group_dante1_fov"]["definition"]["channel"][end]["response_csv"] ==
              "diagnostics/$h19.csv"

        # a plain-text session references the CSV by path
        plain = joinpath(root, "s.briar.toml")
        save_session(plain, sess)
        st = TOML.parsefile(plain)
        @test st["diagnostic_definitions"]["group_dante1_fov"]["definition"]["channel"][end]["response_csv"] ==
              normpath(joinpath(dir, "ch19.csv"))
        # the folder's definition changes: the session's is used, and the change reported
        orig = Controllers._definition(sess, "group_dante1_fov")
        f = joinpath(dir, "group_dante1_fov.toml")
        write(f, replace(read(f, String), "field_of_view_mm = 12" => "field_of_view_mm = 40"))
        s2 = DesignSession(resolution = :coarse)
        load_session!(s2, file)
        @test s2.time.diagnostics[] == tm.diagnostics[]
        @test occursin("differ from the search path's (the saved ones are used): group_dante1_fov", s2.status[])
        r2 = Dict(r.id => r for r in s2.time.catalog[])
        @test r2["group_dante1_fov"].changed == ["view.field_of_view_mm"] && isempty(r2["nif_dante1"].changed)
        @test isempty(diagnostic_diff(Controllers._definition(s2, "group_dante1_fov"), orig))
        use_path_diagnostic!(s2, "group_dante1_fov")
        @test Controllers._definition(s2, "group_dante1_fov").data["view"]["field_of_view_mm"] == 40
        @test isempty(Dict(r.id => r for r in s2.time.catalog[])["group_dante1_fov"].changed)

        # without the folder: the session's definitions alone, still recorded
        remove_diagnostic_folder!(sess, dir)
        @test isempty(diagnostic_folders())
        s3 = DesignSession(resolution = :coarse)
        load_session!(s3, file)
        r3 = Dict(r.id => r for r in s3.time.catalog[])
        @test r3["group_dante1_fov"].origin === :session && r3["group_dante1_fov"].usable
        @test [d.id for d in Controllers._recorded_definitions(s3) if d.def !== nothing] == tm.diagnostics[]
        @test isempty(diagnostic_diff(Controllers._definition(s3, "group_dante1_fov"), orig))

        # sessions before version 6 record Dante 1
        old = session_state(sess)
        delete!(old, "time")
        restore!(s3, old)
        @test s3.time.diagnostics[] == DEFAULT_DIAGNOSTICS
    end
end
