# Quads as the unit of pointing and power, facilities, instruments and their
# configurations, published pulse shapes, the drive history's export, the TaLIS notes
# and the window's rows for them (no GL).

@testset "quads split and merged together; power per quad (no GL)" begin
    s = DesignSession(resolution = :coarse)
    bt = s.beams
    r0 = beam_rows(bt)
    @test all(r -> r.power_unit === :quad, r0) && all(r -> r.power == 8, r0)
    @test quads_to_split(bt) == collect(1:48) && isempty(quads_to_merge(bt))
    P0 = sum(b.power for b in bt.beams[])
    aims0 = Dict(b.quad => b.aim for b in bt.beams[])
    # every quad into its four beams: the rows' powers stay the quads'
    @test split_quads!(bt) == 48 && length(bt.beams[]) == 192
    @test sum(b.power for b in bt.beams[]) ≈ P0 && all(b -> b.power == 2, bt.beams[])
    @test all(b -> b.aim == aims0[b.quad], bt.beams[])
    r1 = beam_rows(bt)
    @test [r.key for r in r1] == [r.key for r in r0] && [r.power for r in r1] == [r.power for r in r0]
    @test all(r -> r.power_unit === :quad, r1) && all(r -> r.power_unit === :quad, beam_rows(bt; by = :quad))
    @test all(r -> r.power_unit === :beam && r.power == 2, beam_rows(bt; by = :beam))
    @test split_quads!(bt) == 0 && occursin("no quad", bt.error[])
    # an edit of a cone's power sets each quad's, shared by its beams
    @test edit_beams!(bt, 1, :power, "6") && beam_rows(bt)[1].power == 6
    @test all(i -> bt.beams[][i].power == 1.5, beam_rows(bt)[1].members)
    # a quad's node: per quad; one of its beams: per beam; half a quad: per beam
    q = bt.beams[][beam_rows(bt)[1].members[1]].quad                   # a quad of that cone
    qi = findall(b -> b.quad == q, bt.beams[])
    @test power_unit(bt, qi) === :quad && beam_set_values(bt, qi).power == 6
    @test power_unit(bt, qi[1:1]) === :beam && beam_set_values(bt, qi[1:1]).power == 1.5
    @test power_unit(bt, qi[1:2]) === :beam
    @test edit_beam_set!(bt, qi, :power, "8") && all(i -> bt.beams[][i].power == 2, qi)
    @test edit_beam_set!(bt, qi[1:1], :power, "3") && bt.beams[][qi[1]].power == 3 && beam_set_values(bt, qi).power == 9
    # a click on a beam of a split quad selects the quad while its node is closed
    g = String(bt.cluster[][qi[1]])
    @test picked_beam_id(s, qi[1]) == "quad:$q@$g"
    toggle_expanded!(s.view, "quad:$q@$g")
    @test picked_beam_id(s, qi[1]) == "beam:" * bt.beams[][qi[1]].name
    # a group's quads merged, then the rest: the first group's are whole again
    gi = node_beams(bt, "group:" * g)
    @test length(quads_to_merge(bt, gi)) == 8 && merge_quads!(bt, gi) == 8 && length(bt.beams[]) == 192 - 24
    @test quads_to_split(bt, node_beams(bt, "group:" * g)) == node_beams(bt, "group:" * g)
    @test picked_beam_id(s, 1) == "beam:" * bt.beams[][1].name      # a quad held as one beam
    @test merge_quads!(bt) == 40 && length(bt.beams[]) == 48 && merge_quads!(bt) == 0
    @test sum(b.power for b in bt.beams[]) ≈ P0 - 16 + 3             # the cone's quads at 6, one of them then at 9
    # single quads still split and merge, with their errors
    @test split_quad!(bt, q) && count(b -> b.quad == q, bt.beams[]) == 4
    @test !split_quad!(bt, q) && occursin("already split", bt.error[])
    @test merge_quad!(bt, q) && !merge_quad!(bt, q) && occursin("not split", bt.error[])
    # beams without quad names: per beam
    o = DesignSession(starter = :foil, resolution = :coarse)
    @test all(r -> r.power_unit === :beam, beam_rows(o.beams)) && isempty(quads_to_split(o.beams))
end

@testset "a quad's beams share one aim point (no GL)" begin
    V = Briareus.Vec3
    s = DesignSession(resolution = :coarse)
    bt = s.beams
    spread(q) = (a = [b.aim for b in bt.beams[] if b.quad == q]; maximum(Briareus.norm(x - y) for x in a, y in a))
    quads = unique(b.quad for b in bt.beams[])
    r0 = beam_rows(bt)
    q0 = beam_rows(bt; by = :quad)
    aim0 = Dict(b.quad => b.aim for b in bt.beams[])
    @test split_quads!(bt) == 48 && beams_label(bt) == "NIF quads (48) as 192 beams"
    # the split leaves every row as it was: aim z, r, the tangential offset, defocus, power
    same(a, b) = all(isapprox(getfield(a, f), getfield(b, f); atol = 5e-4) for f in (:z, :r, :phi_offset, :defocus)) &&
                 a.power == b.power
    @test all(same(a, b) for (a, b) in zip(r0, beam_rows(bt))) &&
          all(same(a, b) for (a, b) in zip(q0, beam_rows(bt; by = :quad)))
    @test all(q -> spread(q) == 0, quads)
    # a cone's aim r and tangential offset move each quad's aim point, its beams together
    @test edit_beams!(bt, 1, :r, "-1.4") && edit_beams!(bt, 1, :phi, "0.3")
    cone = [bt.beams[][i].quad for i in beam_rows(bt)[1].members]
    @test all(q -> spread(q) < 1e-12, quads) && beam_rows(bt)[1].r ≈ -1.4 && beam_rows(bt)[1].phi_offset ≈ 0.3
    @test all(q -> Briareus.norm(first(b.aim for b in bt.beams[] if b.quad == q) - aim0[q]) > 0.3, unique(cone))
    # as the same edits move the quads held as one beam (their ports agree to 0.004°)
    u = DesignSession(resolution = :coarse)
    @test edit_beams!(u.beams, 1, :r, "-1.4") && edit_beams!(u.beams, 1, :phi, "0.3")
    ua = Dict(b.quad => b.aim for b in u.beams.beams[])
    @test maximum(Briareus.norm(b.aim - ua[b.quad]) for b in bt.beams[]) < 1e-3
    # a quad's node, and one of its beams read in the quad's azimuth
    q = cone[1]
    qi = findall(b -> b.quad == q, bt.beams[])
    @test edit_beam_set!(bt, qi, :r, "-0.8") && spread(q) < 1e-12 && beam_set_values(bt, qi).r ≈ -0.8
    @test beam_set_values(bt, qi[1:1]).r ≈ -0.8 && beam_set_values(bt, qi[1:1]).phi ≈ 0.3
    # pointing one beam apart is still possible, from its own row
    @test edit_beam_set!(bt, qi[1:1], :phi, "0.5") && isapprox(spread(q), 0.2; atol = 1e-9)
    # the 192-beam preset aims each quad's beams where the 48 quads aim
    load_preset!(s, :nif_beams)
    @test all(q -> spread(q) < 1e-12, quads) && beams_label(bt) == "NIF beams (192)"
    @test maximum(Briareus.norm(b.aim - aim0[b.quad]) for b in bt.beams[]) < 1e-3
    # the optimizer's parameterization per cone or per quad keeps them together
    p = Controllers._parameterization(bt.beams[], nothing, b -> b.group, :group, (:z, :r, :phi))
    @test p.unit === :quad
    opt = Briareus.apply(p, Briareus.initial_parameters(p) .+ 0.2)
    @test all(q -> (a = [b.aim for b in opt if b.quad == q]; maximum(Briareus.norm(x - y) for x in a, y in a) < 1e-12), quads)
    # the Beams label reads what the table holds after a merge of some quads, or of all
    @test merge_quads!(bt, findall(b -> b.quad == quads[1], bt.beams[])) == 1 &&
          beams_label(bt) == "NIF quads (48) as 189 beams"
    @test merge_quads!(bt) == 47 && beams_label(bt) == "NIF quads (48)"
end

@testset "pulse tables per quad (no GL)" begin
    s = DesignSession(resolution = :coarse)
    bt = s.beams
    times = [0.0, 1.0, 2.0, 3.0]
    pul = cone_fraction_pulses(bt.beams[], times, [0.0, 100.0, 300.0, 0.0], [0.3, 0.3, 0.25, 0.25])
    table = write_pulse_table(joinpath(mktempdir(), "pulses.csv"), pul)
    load_pulses!(s, table)
    src = s.time.pulses[]
    @test src.by === :group && src.unit === :quad && PULSE_UNITS == (:quad, :beam)
    total(τ) = sum(b.power for b in Controllers._pulse_beams(s.time.pulses[], bt.beams[], bt.enabled[], τ))
    P0 = [total(τ) for τ in times]
    b0 = Controllers._pulse_beams(src, bt.beams[], bt.enabled[], 2.0)
    @test P0 ≈ [0.0, 100.0, 300.0, 0.0] && b0 == with_pulses(bt.beams[], pul, 2.0)     # quads as one beam: as read
    # split: each beam fires a quarter of its quad's power, the totals as before
    split_quads!(bt)
    b1 = Controllers._pulse_beams(src, bt.beams[], bt.enabled[], 2.0)
    @test [total(τ) for τ in times] ≈ P0 && length(b1) == 192
    quadpower(bs) = (d = Dict{String,Float64}(); foreach(b -> d[b.quad] = get(d, b.quad, 0.0) + b.power, bs); d)
    @test all(isapprox(quadpower(b1)[q], P; rtol = 1e-12) for (q, P) in quadpower(b0))
    # the irradiance checks' per-beam pulses likewise
    bp = Controllers._beam_pulses(src, bt.beams[], bt.enabled[])
    @test sum(p -> p(2.0), values(bp)) ≈ 300.0
    # read per beam instead: every beam fires its column's power (four times the total)
    set_pulse_unit!(s, :beam)
    @test s.time.pulses[].unit === :beam && [total(τ) for τ in times] ≈ 4 .* P0
    @test Controllers._pulse_beams(s.time.pulses[], bt.beams[], bt.enabled[], 2.0) == with_pulses(bt.beams[], pul, 2.0)
    @test_throws ArgumentError set_pulse_unit!(s, :cone)
    # in a session; one saved before the unit read every beam at its column's power
    st = Controllers._pulse_state(s.time.pulses[])
    @test st["unit"] == "beam" && Controllers._state_pulses(st).unit === :beam
    set_pulse_unit!(s, :quad)
    st = Controllers._pulse_state(s.time.pulses[])
    @test st["unit"] == "quad" && Controllers._state_pulses(st).unit === :quad
    @test Controllers._state_pulses(delete!(copy(st), "unit")).unit === :beam
    # a table per beam name is per beam, whatever is asked; a shape has no unit to set
    perbeam = Dict(Symbol(b.name) => PulseShape(times, [0.0, 1.0, 2.0, 0.0]) for b in bt.beams[])
    @test Controllers._table_pulses(perbeam, :beam, "per beam").unit === :beam
    use_pulse_shape!(s; step = 1.0)
    @test set_pulse_unit!(s, :beam) === s && s.time.pulses[].kind === :shape
    # a retune writes its pulses in the table's unit: per quad, four times the per-beam
    # powers of the same retune on the table read per beam (a quarter of it per beam)
    s2 = DesignSession(resolution = :coarse)
    split_quads!(s2.beams)
    build_scenario!(s2)
    sc, beams = s2.drive.scenario[], s2.beams.beams[]
    enabled = fill(true, length(beams))
    drop = findall(b -> b.quad == beams[1].quad, beams)
    enabled[drop] .= false
    quarter = Dict(k => PulseShape(p.t, p.p ./ 4) for (k, p) in pul)
    rq = retune_pulses(sc, beams, enabled, pul; by = :quad, cap = 1.5)
    rb = retune_pulses(sc, beams, enabled, quarter; by = :quad, cap = 1.5, unit = :beam)
    @test rq.unit === :quad && rb.unit === :beam && keys(rq.pulses) == keys(rb.pulses)
    @test all(isapprox(rq.pulses[k].p, 4 .* rb.pulses[k].p; rtol = 1e-9, atol = 1e-9) for k in keys(rq.pulses))
    @test all(==(0), rq.pulses[Symbol(beams[1].quad)].p) && rq.retune.power ≈ rb.retune.power
    # per cone on a table whose cones hold quads both ways: refused before any solve
    merge_quads!(s2.beams, drop)
    mixed = s2.beams.beams[]
    @test_throws ArgumentError retune_pulses(sc, mixed, fill(true, length(mixed)), pul; by = :group)
    # and so is the optimizer's power per cone (one power per beam of the cone); per quad
    # the check passes on to the next (the drive is out of date after the merge)
    st = s2.optimize
    st.fields[] = [:z, :power]
    st.by[] = :group
    start_optimize!(s2; spawn = false)
    @test occursin("split or merge them all", st.status[]) && occursin(String(mixed[drop[1]].group), st.status[])
    st.by[] = :quad
    start_optimize!(s2; spawn = false)
    @test occursin("compute the drive first", st.status[]) && !st.running[]
end

@testset "facilities and the starters' beams (no GL)" begin
    s = DesignSession(resolution = :coarse)
    bt = s.beams
    @test beams_label(bt) == "NIF quads (48)" && FACILITY_LABELS[:custom] == "Other"
    # a facility with presets loads its first; one without takes the preset's beams out
    set_facility!(bt, :omega; target = s.target)
    @test bt.preset[] === :omega_hohlraum && bt.facility[] === :omega && length(bt.beams[]) == 40
    set_facility!(bt, :omega_ep; target = s.target)
    @test isempty(bt.beams[]) && bt.facility[] === :omega_ep && bt.preset[] === :custom && beams_label(bt) == "OMEGA EP"
    @test any(r -> r.id == "beams" && r.label == "Beams: OMEGA EP", model_tree(s))
    build_scenario!(s)
    @test occursin("no beams", s.drive.status[])
    # imported beams stay when the facility changes, and say whose they are
    set_beams!(bt, nif_quads()[1:4]; preset = :custom)
    @test beams_label(bt) == "OMEGA EP (imported)"
    set_facility!(bt, :custom)
    @test length(bt.beams[]) == 4 && beams_label(bt) == "imported table"
    set_facility!(bt, :nif; target = s.target)
    @test bt.preset[] === :nif_quads && length(bt.beams[]) == 48
    # a new target takes its facility's beams; imported beams and an empty target keep theirs
    @test new_target!(s, :foil) && bt.preset[] === :omega_backlighter && bt.facility[] === :omega
    @test new_target!(s, :frustraum) && bt.preset[] === :nif_quads && bt.facility[] === :nif
    load_preset!(s, :nif_beams)
    @test new_target!(s, :shock_tube) && bt.preset[] === :nif_beams          # the same facility's preset stays
    load_preset!(s, :omega_hohlraum)
    @test new_target!(s, :hohlraum) && bt.preset[] === :nif_quads            # a NIF mock-up gets NIF beams
    @test new_target!(s, :none) && isempty(s.target.parts[]) && bt.preset[] === :nif_quads
    set_beams!(bt, nif_quads()[1:4]; preset = :custom)
    @test new_target!(s, :foil) && bt.preset[] === :custom && length(bt.beams[]) == 4
    @test Set(keys(STARTER_PRESETS)) == Set(STARTERS)
end

@testset "instruments and their configurations (no GL)" begin
    s = DesignSession(resolution = :coarse)
    tm = s.time
    rows = Dict(r.id => r for r in diagnostic_catalog(s))
    @test rows["nif_dante1"].instrument == "NIF Dante 1" && rows["nif_dante1"].config == "default"
    @test rows["nif_dante1_2012"].instrument == "NIF Dante 1" && rows["nif_dante1_2012"].config == "2012–14"
    @test rows["nif_sxi_upper"].instrument == "NIF SXI upper" && rows["nif_sxi_upper"].config == "hard channel"
    @test rows["omega_dante"].instrument == "Omega Dante" && rows["omega_dante"].facility == "Omega"
    @test Controllers._config_label("Dante (2012–14)", "Dante") == "2012–14" &&
          Controllers._config_label("Dante, soft", "Dante") == "soft" && Controllers._config_label("Other", "Dante") == "Other"
    # one row per instrument of the session's facility, by kind; the recorded configuration
    ins = diagnostic_instruments(s)
    names = [r.instrument for r in ins]
    @test allunique(names) && "NIF Dante 1" in names && !("Omega Dante" in names)
    d1 = ins[findfirst(==("NIF Dante 1"), names)]
    @test length(d1.configs) == 7 && d1.configs[1].id == "nif_dante1" && d1.recorded == ["nif_dante1"] && d1.here
    @test [r.group for r in ins if r.head] == unique(r.group for r in ins)
    # the window records an instrument once and swaps its configuration in place
    set_instrument!(s, "NIF SXI lower")
    set_configuration!(s, "nif_dante1", "nif_dante1_2012")
    @test tm.diagnostics[] == ["nif_dante1_2012", "nif_sxi_lower"]
    set_instrument!(s, "NIF Dante 1")                                  # recorded already: as it is
    @test tm.diagnostics[] == ["nif_dante1_2012", "nif_sxi_lower"]
    # an instrument that does not build is not recorded, and the status says why
    set_instrument!(s, "NIF Dante 2")
    @test tm.diagnostics[] == ["nif_dante1_2012", "nif_sxi_lower"] && occursin("NIF Dante 2 cannot be recorded", s.status[])
    @test Controllers._definition(s, "nif_dante2") === nothing
    # a script may record several definitions of one instrument side by side
    set_diagnostic!(s, "nif_dante1")
    @test diagnostic_instruments(s)[1].recorded == ["nif_dante1_2012", "nif_dante1"]
    set_configuration!(s, "nif_dante1", "nif_dante1_2012")              # recorded already: dropped
    @test tm.diagnostics[] == ["nif_dante1_2012", "nif_sxi_lower"]
    set_instrument!(s, "NIF Dante 1", false)
    @test tm.diagnostics[] == ["nif_sxi_lower"]
    # another facility's list; an instrument recorded for the other facility stays listed
    set_facility!(s.beams, :omega; target = s.target)
    ins = diagnostic_instruments(s)
    @test "Omega Dante" in [r.instrument for r in ins] && !("NIF Dante 1" in [r.instrument for r in ins])
    @test only(r for r in ins if r.instrument == "NIF SXI lower").here == false
    @test diagnostic_here(s, "") && diagnostic_here(s, "OMEGA") && !diagnostic_here(s, "NIF")
    set_facility!(s.beams, :custom)
    @test all(r -> r.here, diagnostic_instruments(s)) && "Omega Dante" in [r.instrument for r in diagnostic_instruments(s)]
end

@testset "published pulse shapes (no GL)" begin
    s = DesignSession(resolution = :coarse)
    @test PULSE_SHAPES == (:n221204, :n191118, :n210601) && Set(keys(PULSE_SHAPE_LABELS)) == Set(PULSE_SHAPES)
    for key in PULSE_SHAPES
        use_pulse_shape!(s, key; step = 1.0)
        src = s.time.pulses[]
        t, P = Controllers._published_pulse(key)
        @test src.key === key && maximum(src.shape.(t)) ≈ 1 && src.nodes == t
        @test occursin(uppercase(String(key)), src.label)
        # in a session: the shape by name
        st = Controllers._pulse_state(src)
        @test st["shape"] == String(key) && Controllers._state_pulses(st).key === key
    end
    @test Controllers._state_pulses(Dict("kind" => "shape", "step" => 0.5)).key === :n221204   # sessions before the name
    @test_throws ArgumentError use_pulse_shape!(s, :n000000)
    use_pulse_shape!(s; step = 1.0)
    @test s.time.pulses[].key === :n221204
end

@testset "the drive history in the results export (no GL)" begin
    s = DesignSession(resolution = :coarse)
    set_instrument!(s, "NIF SXI lower")
    use_pulse_shape!(s; step = 3.0)
    build_scenario!(s)
    dir = mktempdir()
    # without a history: the drive's files only
    files = export_results(s, dir)
    @test !any(f -> startswith(basename(f), "history"), files) && !occursin("history", s.status[])
    @test_throws ArgumentError export_history(s, dir)
    start_history!(s)
    h = s.time.history[]
    files = basename.(export_results(s, dir))
    @test occursin("with the drive history", s.status[])
    @test issubset(["history.csv", "history_laser_Wcm2.csv", "history_radiosity_Wcm2.csv", "history_incident_Wcm2.csv",
                    "image_nif_sxi_lower_integrated.csv", "session.briar.toml", "elements.csv"], files)
    lines = readlines(joinpath(dir, "history.csv"))
    head = split(lines[1], ',')
    @test head[1:2] == ["time_ns", "Tr_eV"] && "a2_over_a0" in head && "energy_laser_kJ" in head &&
          "nif_dante1_Tr_eV" in head && "nif_sxi_lower_TW_per_sr" in head
    @test length(lines) == length(h.time) + 1
    vals = [parse.(Float64, split(l, ',')) for l in lines[2:end]]
    col(name) = [v[findfirst(==(name), head)] for v in vals]
    @test col("time_ns") ≈ h.time && isapprox(col("Tr_eV"), h.Tr; rtol = 1e-7) &&
          isapprox(col("a2_over_a0"), h.legendre[2, :]; rtol = 1e-6, atol = 1e-12)
    # the fields: a row per element of the drive's target, a column per step
    t = s.drive.scenario[].target
    fl = readlines(joinpath(dir, "history_incident_Wcm2.csv"))
    @test length(fl) == nelements(t) + 1 && length(split(fl[1], ',')) == 2 + length(h.time)
    row = split(fl[2], ',')
    @test row[1] == "1" && Symbol(row[2]) === t.components[t.elem_comp[1]].name &&
          isapprox(parse.(Float64, row[3:end]), Float64.(h.fields[:incident][1, :]); rtol = 1e-5)
    # the imager's image: a comment line, the x coordinates, a row per y
    im = readlines(joinpath(dir, "image_nif_sxi_lower_integrated.csv"))
    d = only(x for x in h.diagnostics if x.id == "nif_sxi_lower")
    n = d.view.camera.nx
    @test startswith(im[1], "# NIF SXI lower") && length(im) == n + 2 && length(split(im[2], ',')) == n + 1
    @test maximum(maximum(parse.(Float64, split(l, ',')[2:end])) for l in im[3:end]) ≈ 1e5 * maximum(d.frames.images[1]) rtol = 1e-5
    # an image as the Diagnostics tab draws it, and in the image view (then exported too)
    di = diagnostic_image(s, "nif_sxi_lower", 1)
    @test di !== nothing && size(di.values) == (n, n) && occursin("time-integrated", di.title)
    @test diagnostic_image(s, "nif_dante1", 0) === nothing && diagnostic_image(s, "nope", 0) === nothing
    show_diagnostic_image!(s, "nif_sxi_lower", 1)
    @test s.drive.view[] === :image && s.drive.image[].title == di.title
    @test "image.csv" in basename.(export_results(s, dir))
    # a history out of date is not exported, and the status says so
    edit_beam_set!(s.beams, [1], :power, "7")
    files = basename.(export_results(s, mktempdir()))
    @test !("history.csv" in files) && occursin("out of date", s.status[])
end

@testset "TaLIS notes: short lines (no GL)" begin
    tl = TaLISState()
    notes = talis_notes(tl)
    @test notes[1].kind === :tolerance && !occursin(TALIS_TOLERANCE_NOTE, notes[1].text)
    scope = [n for n in notes if n.kind === :scope]
    @test length(scope) == 1 && scope[1].text == TALIS_SCOPE_NOTE && !occursin("blow-by", TALIS_SCOPE_NOTE)
    @test all(n -> length(n.text) < 130, notes)
    @test !any(n -> occursin("blow-by", n), Briareus.TALIS_OUT_OF_SCOPE)
end

@testset "the window's rows for these (no window)" begin
    s = DesignSession(resolution = :coarse)
    sh = BriareusGUI.BriareusShell(s)
    B = BriareusGUI
    B._SHELL[] = sh
    try
        refresh() = B._refresh!(sh)
        # the row menu's numbers: shown, included, quads to split and merge, piece roles
        @test B.br_node_menu("group:inner23") == [1, 1, 8, 0, 0, 0, 0, 0, 0]
        m = B.br_node_menu("comp:hohlraum_cap_top")
        @test m[1:4] == [1, 1, 0, 0] && m[5] == 1 && m[6] == findfirst(==(:wall), PIECE_ROLES) - 1 && m[7:9] == [1, 1, 1]
        # split from the Beams row: the selected quad beams become their quads' nodes
        q = s.beams.beams[][1].quad
        g = String(s.beams.cluster[][1])
        @test s.beams.beams[][1].name == q                             # a quad held as one beam is named as its quad
        B.br_select("beam:" * q)
        B.br_split_quads("beams")
        refresh()
        @test length(s.beams.beams[]) == 192 && s.view.selection[] == "quad:$q@$g"
        @test sh.shown["beamPowerUnit"] == "quad" && sh.shown["bnPowerUnit"] == "quad" && sh.shown["bn_power"] == "8"
        @test sh.shown["bnSplit"] == 0 && sh.shown["bnMerge"] == 1 && B.br_node_menu("beams")[3:4] == [0, 48]
        B.br_merge_quads("")                                           # the selection's quad
        refresh()
        @test length(s.beams.beams[]) == 192 - 3 && s.view.selection[] == "beam:" * q && sh.shown["bnSplit"] == 1
        B.br_merge_quads("beams")
        @test length(s.beams.beams[]) == 48
        # several nodes removed at once; the status says what stays where it is
        B.br_add_part("disk")
        B.br_select("part:hohlraum")
        B.br_select_add("part:disk_1")
        refresh()
        @test sh.shown["selectionCount"] == 2
        B.br_remove_selected()
        refresh()
        @test [p.name for p in s.target.parts[]] == [:capsule] && occursin("capsule keeps its place", s.status[])
        B.br_remove_selected()                                         # the target is selected: nothing to remove
        @test occursin("nothing to remove", s.status[])
        # an empty target, a bare hohlraum part inserted, a starter with its facility's beams
        B.br_new_target("empty")
        refresh()
        @test isempty(s.target.parts[]) && sh.shown["partCount"] == 0 && !sh.shown["hasHohlraum"]
        B.br_add_part("hohlraum")
        refresh()
        @test [p.name for p in s.target.parts[]] == [:hohlraum_1] && sh.shown["hasHohlraum"] &&
              s.view.selection[] == "part:hohlraum_1"
        B.br_new_target("foil")
        refresh()
        @test sh.shown["facility"] == "omega" && s.view.selection[] == "target" && occursin("OMEGA backlighter", s.status[])
        # a facility without beams, and back
        B.br_set_facility("omega_ep")
        refresh()
        @test sh.shown["beamCount"] == 0 && sh.shown["presetLabel"] == "OMEGA EP" && occursin("import a beam table", s.status[])
        B.br_new_target("hohlraum")
        B.br_set_facility("nif")
        refresh()
        @test sh.shown["beamCount"] == 48 && sh.shown["facility"] == "nif"
        # the Diagnostics tab's rows: an instrument once, its configurations, the categories
        rows = sh.rows[:diagModel]
        d = only(r for r in rows if r.name == "NIF Dante 1")
        @test d.ticked && d.first && d.config == "default" && d.kind == "xray_power" && !d.foreign &&
              split(d.configIds, "|")[1] == "nif_dante1" && length(split(d.configLabels, "|")) == 7
        @test sh.shown["diagCounts"] == [1, 0, 0] && sh.shown["diagKinds"] == ["xray_power", "xray_imager", "radiography"]
        B.br_diag_config("nif_dante1", "nif_dante1_2021")
        B.br_diag_record("NIF SXI upper", true)
        B.br_diag_category("xray_imager")
        refresh()
        d = only(r for r in sh.rows[:diagModel] if r.name == "NIF Dante 1")
        u = only(r for r in sh.rows[:diagModel] if r.name == "NIF SXI upper")
        @test d.id == "nif_dante1_2021" && d.config == "from October 2021" && u.ticked && u.config == "hard channel"
        @test sh.shown["diagCounts"] == [1, 1, 0] && sh.shown["diagCategory"] == "xray_imager"
        B.br_diag_record("NIF Dante 2", true)                          # does not build: not recorded
        refresh()
        @test !only(r for r in sh.rows[:diagModel] if r.name == "NIF Dante 2").ticked
        B.br_time_diagnostic("nif_dante1", true)                       # a second configuration, as a script records it
        refresh()
        ds = [r for r in sh.rows[:diagModel] if r.name == "NIF Dante 1"]
        @test length(ds) == 2 && all(r -> r.ticked, ds) && [r.first for r in ds] == [true, false]
        B.br_time_diagnostic("nif_dante1", false)
        # pulses from the menu, the history, what the 3D view shows, the tab's image
        B.br_time_chart("diagnostics")                                 # the chart of the Diagnostics tab:
        refresh()
        @test s.time.chart[] === :drive                                # the Time tab falls back to the drive's
        B.br_time_shape("n191118")
        refresh()
        @test s.time.pulses[].key === :n191118 && occursin("N191118", sh.shown["pulsesLabel"]) &&
              sh.shown["pulseShapeKeys"] == String.(collect(PULSE_SHAPES))
        use_pulse_shape!(s; step = 3.0)
        build_scenario!(s)
        start_history!(s)
        refresh()
        @test sh.shown["timeMode"] == "beams" && sh.shown["diagEmpty"] == "" && sh.shown["diagImageIndex"] == 0
        @test length(sh.rows[:diagImageModel]) == 2 && !isempty(sh.rows[:diagLegendModel])
        @test sh.diag_image.ax.title[] != "No image yet"
        B.br_time_mode("history")
        refresh()
        n = length(s.time.history[].time)
        @test sh.shown["timeMode"] == "history" && sh.shown["timeStep"] == n
        B.br_time_step(2)
        refresh()
        B.br_time_mode("beams")
        refresh()
        @test sh.shown["timeMode"] == "beams" && time_label(s) == ""
        B.br_time_mode("history")                                      # back at the step shown last
        refresh()
        @test sh.shown["timeStep"] == 2
        B.br_diag_pick_image("nif_sxi_upper", 1)
        refresh()
        @test sh.shown["diagImageIndex"] == 1 && occursin("time-integrated", sh.diag_image.ax.title[])
        # the tick's bits: the Diagnostics tab's canvases
        sh.diag.dirty[] = true
        @test B.br_tick() & 16 == 16
    finally
        B._SHELL[] = nothing
    end
end
