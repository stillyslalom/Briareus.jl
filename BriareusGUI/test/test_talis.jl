# This file tests the TaLIS checks through a DesignSession, first on a small foil with one
# beam and a plate upstream (the analytic case of the core tests), then on the coarse
# hohlraum with the NIF quads and a diagnostic window.

using LinearAlgebra: dot, cross

@testset "TaLIS checks: foil and one beam (no GL)" begin
    spot = SuperGaussianSpot(0.3; order = 4.0)
    b = Beam("b"; theta = 0, phi = 0, power = 1.0, spot, fnumber = 8.0)
    x0 = 1.6                                         # plate edge, 2 mm upstream of the foil
    plate = Component(:plate, translate(rectangle_mesh(2.0, 2.0; nu = 4, nv = 4), (x0 + 1, 0, 2.0)),
                      gold(); two_sided = true)
    t = Target(foil(:f, (0, 0, 0), (0, 0, 1); width = 4.0, n = 8), plate)
    tl, rv = TaLISState(), ReviewState()
    default_selection!(tl, t)
    @test tl.capture[] == [:f] && tl.hardware[] == [:plate]      # the first component captures
    @test_throws ArgumentError set_talis_role!(tl, :plate, :target)

    # before a check, the notes already carry the tolerance caveat
    n0 = talis_notes(tl)
    @test n0[1].kind === :tolerance && occursin("lateral allowance", n0[1].text) && occursin("TCC", n0[1].text)
    @test talis_brief(tl) == "not run" && isempty(talis_rows(tl))

    kw = talis_kwargs(tl, rv)
    @test kw.capture == (:f,) && kw.hardware == (:f, :plate) && kw.omega1_components == kw.hardware
    @test kw.tolerance == 0.5 && kw.pivot === :tcc && kw.pointing == 0.05 && kw.level == 1
    start_talis!(tl, rv, t, [b]; spawn = false)
    r = tl.report[]
    @test r isa TaLISReport && !tl.running[] && startswith(tl.status[], "done")
    ref = talis_check(t, [b]; kw...)
    key(rep) = [(x.beam, x.rule, x.component, x.status, x.margin) for x in rep.rows]
    @test key(r) == key(ref)
    c3 = only(x for x in r.rows if x.rule === :clearance_3w && x.component === :plate)
    @test c3.status === :pass && c3.unit === :mm
    rr = talis_rule_rows(tl)
    @test [x.rule for x in rr] == [:capture_3w, :incidence_3w, :clearance_3w]
    @test all(x -> x.pass + x.review + x.fail == count(y -> y.rule === x.rule, r.rows), rr)
    @test occursin("contour", rr[1].requirement) && rr[2].label == "3ω incidence"
    rows = talis_rows(tl)
    rank(s) = s === :fail ? 0 : s === :review ? 1 : 2
    @test length(rows) == length(r.rows) && issorted([(rank(x.status), x.margin) for x in rows])
    notes = talis_notes(tl)
    @test notes[1].kind === :tolerance && occursin("lateral allowance", notes[1].text)
    @test any(n -> n.kind === :skipped && occursin("1ω", n.text), notes)    # not a NIF beam
    @test [n.text for n in notes if n.kind === :scope] == [TALIS_SCOPE_NOTE]      # one short line
    @test notes[end].kind === :source && occursin(TALIS_SOURCE, notes[end].text)

    # The tolerance interpretation reaches talis_check. With no tolerance, r·tan(0.5°) is
    # removed at the edge point of the plate, at distance hypot(x0, 2) from TCC.
    @test edit_talis!(tl, :tolerance, "0") && tl.tolerance[] == 0
    start_talis!(tl, rv, t, [b]; spawn = false)
    c0 = only(x for x in tl.report[].rows if x.rule === :clearance_3w && x.component === :plate)
    # Each worst value is a lower bound within max(1e-4, 0.01 |worst|), the defaults.
    @test c0.worst - c3.worst ≈ tand(0.5) * hypot(x0, 2.0) atol = 0.01 * (c0.worst + c3.worst) + 2e-4
    @test edit_talis!(tl, :tolerance, "0.5")
    # A pivot at the edge point gives zero allowance at that point.
    @test edit_talis!(tl, :pivot_point, "1.6, 0, 2")
    set_talis_pivot!(tl, :point)
    @test talis_kwargs(tl, rv).pivot == (1.6, 0.0, 2.0)
    start_talis!(tl, rv, t, [b]; spawn = false)
    cP = only(x for x in tl.report[].rows if x.rule === :clearance_3w && x.component === :plate)
    @test cP.worst ≈ c0.worst atol = 0.01 * c0.worst + 1e-4
    @test occursin("(1.6, 0, 2)", talis_notes(tl)[1].text)
    set_talis_pivot!(tl, :aim)
    @test talis_kwargs(tl, rv).pivot === :aim
    @test_throws ArgumentError set_talis_pivot!(tl, :lens)
    @test !edit_talis!(tl, :pivot_point, "1, 2") && occursin("three coordinates", tl.error[])
    @test !edit_talis!(tl, :tolerance, "-1") && !edit_talis!(tl, :focus, "1")
    # the Review step's error stack and sampling are used
    edit_review!(rv, :pointing, "0.1")
    rv.quality[] = :standard
    @test talis_kwargs(tl, rv).pointing == 0.1 && talis_kwargs(tl, rv).level == 2

    start_talis!(tl, rv, t, [b]; spawn = false)                    # on the current settings
    path = export_talis(tl, joinpath(mktempdir(), "talis.csv"))
    lines = readlines(path)
    @test count(!startswith("#"), lines) == 1 + length(tl.report[].rows) && tl.status[] == "wrote talis.csv"
    @test "# setting pivot = aim" in lines && any(startswith("# out of scope: "), lines)
    @test_throws ArgumentError export_talis(TaLISState(), path)
    # stale on the target, the beams, or any setting the check used
    @test !talis_stale(tl, rv, t, [b]) && talis_stale(tl, rv, t, [remake(b; aim = (0.1, 0.0, 0.0))])
    changes = (() -> edit_talis!(tl, :tolerance, "0.3"), () -> set_talis_pivot!(tl, :tcc),
               () -> (tl.omega1[] = false), () -> (rv.quality[] = :quick),
               () -> edit_review!(rv, :focus, "2"), () -> set_talis_role!(tl, :plate, :none))
    stale = map(changes) do change!
        start_talis!(tl, rv, t, [b]; spawn = false)
        fresh = !talis_stale(tl, rv, t, [b])
        change!()
        fresh && talis_stale(tl, rv, t, [b])
    end
    @test all(stale)

    # diagnostic hardware: 3ω-checked as hardware, given to the 1ω diagnostic rule
    # instead of the 0.2 mm 1ω clearance
    set_talis_role!(tl, :plate, :diagnostic)
    @test talis_role(tl, :plate) === :diagnostic && talis_role(tl, :f) === :capture
    kd = talis_kwargs(tl, rv)
    @test kd.diagnostics == (:plate,) && :plate in kd.hardware && !(:plate in kd.omega1_components)
    @test_throws ArgumentError set_talis_role!(tl, :plate, :diagnostics)

    # no capture component, no beams
    tl2 = TaLISState()
    default_selection!(tl2, t)
    set_talis_role!(tl2, :f, :none)
    start_talis!(tl2, rv, t, [b]; spawn = false)
    @test occursin("capture", tl2.status[]) && tl2.report[] === nothing
    start_talis!(tl2, rv, t, Beam{Float64}[]; spawn = false)
    @test tl2.status[] == "no beams are enabled"

    # On a worker the result arrives through the queue, and a cancel drops it.
    q = Channel{Any}(Inf)
    tl3 = TaLISState()
    default_selection!(tl3, t)
    start_talis!(tl3, rv, t, [b]; deliver = f -> put!(q, f), spawn = true)
    @test tl3.running[]
    cancel_talis!(tl3)
    @test !tl3.running[] && tl3.status[] == "canceled"
    t0 = time()
    while !isready(q) && time() - t0 < 120
        sleep(0.02)
    end
    take!(q)()
    @test tl3.report[] === nothing && tl3.status[] == "canceled"
    start_talis!(tl3, rv, t, [b]; deliver = f -> put!(q, f), spawn = true)
    drain_until(() -> !tl3.running[], q)
    @test tl3.report[] isa TaLISReport
end

@testset "1ω disposition of the TaLIS check (no GL)" begin
    # one NIF beamlet, a plate that faces it in its 1ω footprint, a small plate just
    # outside the near edge of the footprint and one on the dark side
    bl = first(nif_beams())
    k = direction(bl)
    pq = Briareus.nif_quads()[findfirst(q -> q.quad == bl.quad, Briareus.nif_quads())].port
    e1, e2 = Briareus.port_basis(bl.port)
    ô = dot(bl.port - pq, e2) >= 0 ? e2 : -e2
    Rm = hcat(ô, cross(k, ô), k)
    plate(name, w, c) = Component(name, transform(rectangle_mesh(w, w; nu = 2, nv = 2), Rm, c), gold(); two_sided = true)
    t = Target(plate(:mirror, 4.0, 8.0 * ô), plate(:near, 0.1, 3.9 * ô), plate(:far, 0.1, -3.9 * ô))
    tl, rv = TaLISState(), ReviewState()
    @test disposition_of(tl, t) === nothing && disposition_status(tl, t) === nothing
    @test disposition_info(tl, t, [bl], 1) == ""
    _assign = Controllers._assign!
    _assign(tl.capture, [:mirror]); _assign(tl.hardware, [:near, :far])
    start_talis!(tl, rv, t, [bl]; spawn = false)
    d = disposition_of(tl, t)
    kw = talis_kwargs(tl, rv)
    @test d isa Omega1Disposition && d.status == omega1_disposition(t, [bl]; disposition_kwargs(t, kw)...).status
    @test disposition_kwargs(t, kw).level == kw.level && disposition_kwargs(t, kw).margin_range == DISPOSITION_RANGE
    of(name) = [j for j in eachindex(t.tri_area) if t.components[t.tri_comp[j]].name === name]
    st = disposition_status(tl, t)
    @test Set(st[of(:mirror)]) ⊆ Set([:fail, :specular_clear]) && Set(st[of(:near)]) == Set([:tight])
    @test Set(st[of(:far)]) == Set([:clear])
    # the tight threshold and the curved-safe filter change the drawing, not the check
    set_disposition_filter!(tl; tight = 0.05)
    @test Set(disposition_status(tl, t)[of(:near)]) == Set([:clear]) && tl.disposition[] === d
    set_disposition_filter!(tl; tight = DISPOSITION_RANGE)
    @test Set(disposition_status(tl, t)[of(:near)]) == Set([:tight])
    @test_throws ArgumentError set_disposition_filter!(tl; tight = 6.0)
    set_disposition_filter!(tl; tight = 1.0, hide_curved = true)
    @test tl.hide_curved[] && !(:curved_safe in disposition_status(tl, t))
    # what a facet reports
    info = disposition_info(tl, t, [bl], first(of(:near)))
    @test startswith(info, "tight; 1ω margin +") && occursin(bl.name, info) && occursin("beamlet 1", info)
    im = disposition_info(tl, t, [bl], first(of(:mirror)))
    @test occursin("radii of curvature flat, flat", im) && occursin("from open port", im)
    @test occursin("above 5 mm", disposition_info(tl, t, [bl], first(of(:far))))
    # a row under inspection: the approach of a clearance row with what the view draws of
    # it, the reflection and cone of a management row decided on reflection
    rows = talis_rows(tl)
    @test tl.focus[] === nothing && talis_focus_text(tl) == ""
    focus_talis_row!(tl, findfirst(x -> x.rule === :clearance_1w && x.component === :near, rows))
    f = tl.focus[]
    @test f.row.component === :near && f.approach.component === :near && f.approach.x_hw == f.row.point
    @test length(f.surface.outline) == 4 && f.surface.caliper == (f.approach.x_hw, f.approach.x_beam)
    @test f.reflection === nothing && f.cone === nothing
    @test length(split(talis_focus_text(tl), '\n')) == 8 && startswith(talis_focus_text(tl), "nominal")
    # the slice of the row, as the tab draws it: the numbers of the core slice in order
    sl = tl.focus[].slice
    @test sl.distance == tl.focus[].approach.distance && sl.component === :near
    js = slice_drawing(sl)
    nums(key) = parse.(Float64, split(match(Regex("\"$key\":\\[([^\\]]*)\\]"), js).captures[1], ','; keepempty = false))
    @test startswith(js, "{\"window\":[") && endswith(js, "}") && !occursin("NaN", js) && !occursin(" ", js)
    @test maximum(abs.(nums("window") .- collect(sl.window))) < 1e-4
    @test maximum(abs.(nums("nominal") .- [sl.nominal.center..., sl.nominal.half])) < 1e-4
    @test maximum(abs.(nums("worst") .- [sl.worst.center..., sl.worst.half])) < 1e-4
    @test length(nums("hardware")) == 4 * length(sl.hardware) && length(nums("caliper")) == 4
    @test occursin("1ω margin", slice_text(sl)) && startswith(slice_text(sl), tl.focus[].approach.beam)
    im = findfirst(x -> x.rule === :management_1w && x.component === :mirror, rows)
    focus_talis_row!(tl, im)
    @test tl.focus[].index == im && tl.focus[].reflection.component === :mirror
    @test tl.focus[].cone.half_angle == 20 + 2 * 0.5 && occursin("reflection", talis_focus_text(tl)) &&
          occursin("ringed in the view", talis_focus_text(tl)) && length(tl.focus[].cone.ring) == 32
    focus_talis_row!(tl, findfirst(x -> x.rule === :capture_3w, rows))
    @test tl.focus[].approach === nothing && tl.focus[].surface === nothing && talis_focus_text(tl) == ""
    @test tl.focus[].slice === nothing
    focus_talis_row!(tl, length(rows) + 1)
    @test tl.focus[] === nothing
    # the pinch mode: every beam within the threshold of a part. The whole quad of the
    # beamlet lights both sides of the 3ω axis, so the quad and the beamlet both come
    # within 1 mm of the near plate, and only the beamlet is near in a check of it alone.
    @test isempty(pinch_approaches(TaLISState(), :near))
    pa = pinch_approaches(tl, :near)
    @test [a.beam for a in pa] == [bl.name] && abs(pa[1].margin - f.approach.margin) < 1e-4
    @test isempty(pinch_approaches(tl, :near; threshold = 0.05)) && isempty(pinch_approaches(tl, :far))
    quad = only(q for q in nif_quads() if q.quad == bl.quad)
    tl2 = TaLISState()
    _assign(tl2.capture, [:mirror]); _assign(tl2.hardware, [:near, :far])
    start_talis!(tl2, rv, t, [bl, quad]; spawn = false)
    p2 = pinch_approaches(tl2, :near)
    @test Set(a.beam for a in p2) == Set([bl.name, quad.name]) && issorted([a.margin for a in p2])
    r2 = talis_rows(tl2)
    i2 = findfirst(x -> x.rule === :clearance_1w && x.component === :near && x.beam == bl.name, r2)
    focus_talis_row!(tl2, i2)
    @test isempty(tl2.focus[].pinch)
    set_talis_pinch!(tl2, true)
    @test tl2.focus[].index == i2 && [p.approach.beam for p in tl2.focus[].pinch] == [quad.name]
    @test length(tl2.focus[].pinch[1].surface.outline) >= 4 && occursin(quad.name, talis_focus_text(tl2))
    set_talis_pinch!(tl2, false)
    @test isempty(tl2.focus[].pinch) && !tl2.pinch[]
    focus_talis_row!(tl, im)
    # another target has no disposition, and the 1ω rules off leave none
    @test disposition_of(tl, Target(plate(:mirror, 4.0, 8.0 * ô))) === nothing
    _assign(tl.omega1, false)
    start_talis!(tl, rv, t, [bl]; spawn = false)
    @test tl.disposition[] === nothing && tl.report[] isa TaLISReport
    @test tl.focus[] === nothing                                      # a new check clears the focus
    @test Set(DISPOSITION_MODES) ⊆ Set(COLOR_MODES) && all(m -> haskey(COLOR_LABELS, m), COLOR_MODES)
end

@testset "TaLIS checks: hohlraum with a window (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    @test add_window!(sess.target)                            # closed, at φ = 90°
    tl = sess.talis
    @test tl.capture[] == [:hohlraum_wall] && :hohlraum_window_1 in tl.hardware[] && :hohlraum_cap_top in tl.hardware[]
    @test !(:capsule in tl.hardware[]) && !(:capsule in tl.capture[])   # drive physics; the capsule is not hardware
    @test talis_kwargs(tl, sess.review).leh == (:hohlraum_cap_top, :hohlraum_cap_bottom)
    @test :hohlraum_window_1 in sess.review.hardware[]                 # clearance checks it too
    @test task_status(sess, :talis)[1] === :todo
    start_talis!(sess; spawn = false)
    r = tl.report[]
    # The pointing of the template grazes the LEH lips, which is an LEH-clearance review and
    # does not count as a failure. Apart from the closed window, which is hardware in the 1%
    # contours of beams that land beside it, the only failures are the incidence of quads on
    # the wall. Their axes meet the wall below 70°, and the steepest rays of their f-cones
    # (axis + 3.6° at f/7.9) above it.
    lip = [x for x in r.rows if x.rule === :clearance_3w && x.component in (:hohlraum_cap_top, :hohlraum_cap_bottom)]
    @test any(x -> x.status === :review, lip) && all(x -> x.status in (:pass, :review), lip)
    @test all(x -> occursin("LEH clearance", x.note), (x for x in lip if x.status === :review))
    bad = [x for x in r.rows if x.status === :fail && x.component !== :hohlraum_window_1]
    @test !isempty(bad) && all(x -> x.rule === :incidence_3w && x.component === :hohlraum_wall &&
                                   x.nominal - atand(1 / (2 * 7.9)) < 70, bad)   # quads are f/7.9
    @test any(x -> x.status === :fail && x.rule === :clearance_3w, (x for x in r.rows if x.component === :hohlraum_window_1))
    @test !any(x -> x.component === :capsule && x.rule !== :tas_stay_in, r.rows)   # stay-in covers every part
    rules = Set(x.rule for x in r.rows)
    @test :tas_stay_in in rules && :clearance_1w in rules && :clearance_3w in rules
    @test any(x -> x.component === :hohlraum_window_1, r.rows)
    @test length(unique(x.beam for x in r.rows if x.rule !== :tas_stay_in)) == 48
    @test r.settings.capture == [:hohlraum_wall]
    # Every quad is deposited on the wall, so capture rows have no free substrate edge to
    # measure.
    @test all(x -> x.component === :hohlraum_wall, (x for x in r.rows if x.rule === :capture_3w))
    @test occursin("48 beams", talis_brief(tl))
    st, txt = task_status(sess, :talis)
    @test st === :attention && occursin(r"\d+ fail, \d+ review", txt)  # LEH reviews, incidence fails
    # The roles follow the target. A removed window leaves the selection, and a new one is
    # checked.
    set_talis_role!(tl, :capsule, :hardware)
    edit_part!(sess.target, 1, :length, "10.0")
    @test :capsule in tl.hardware[] && !(:capsule in tl.capture[])        # kept as chosen
    @test !(:capsule in talis_kwargs(tl, sess.review).leh)
    @test talis_stale(tl, sess.review, sess.target.target[], active_beams(sess.beams))
    @test task_status(sess, :talis)[1] === :attention
    remove_window!(sess.target, 1)
    @test !(:hohlraum_window_1 in tl.hardware[]) && !(:hohlraum_window_1 in sess.review.hardware[])
    add_window!(sess.target)
    @test :hohlraum_window_1 in tl.hardware[] && :hohlraum_window_1 in sess.review.hardware[]
    # switching to a foil chooses afresh
    @test new_target!(sess, :foil)
    @test tl.capture[] == [:foil] && tl.hardware[] == [:main]
    set_review_view!(sess.review, :talis)
    @test sess.review.view[] === :talis
    @test_throws ArgumentError set_review_view!(sess.review, :tally)
end

@testset "keep-out map (no GL)" begin
    spot = SuperGaussianSpot(0.3; order = 4.0)
    b = Beam("b"; theta = 0, phi = 0, power = 1.0, spot, fnumber = 8.0)
    cs, rv, tl = CorridorState(), ReviewState(), TaLISState()
    @test cs.radius[] == 50.0 && !cs.shown[] && cs.map[] === nothing
    kw = corridor_kwargs(rv, tl)
    @test kw == (; tolerance = 0.5, pivot = :tcc, pointing = 0.05, target_offset = 0.025, hardware_tolerance = 0.0,
                 focus = 1.0, omega1 = true)
    @test !corridor_stale(cs, rv, tl, [b], SightLine[])
    sl = [SightLine("cam", 90, 0; field_of_view = 10.0)]
    update_corridor!(cs, rv, tl, [b], sl)
    m = cs.map[]
    ref = corridor_map([b]; sight_lines = sl, radius = 50.0, kw...)
    @test m.margin == ref.margin && m.sources == [(:omega3, "b"), (:omega1, "b"), (:sight, "cam")]
    @test startswith(cs.status[], "sphere of 50 mm about TCC: ") && occursin("1 beams, 1 lines of sight", cs.status[])
    @test occursin("widest margin", cs.status[]) && !occursin(" %", cs.status[])
    @test !corridor_stale(cs, rv, tl, [b], sl) && corridor_stale(cs, rv, tl, Beam[], sl) &&
          corridor_stale(cs, rv, tl, [b], SightLine[])
    tl.tolerance[] = 1.0
    @test corridor_stale(cs, rv, tl, [b], sl)
    tl.tolerance[] = 0.5
    # the surface: the corners of the cells, two triangles per cell but one at each pole
    g = corridor_surface(m)
    nθ, nφ = size(m.margin)
    @test length(g.vertices) == (nθ + 1) * nφ && length(g.faces) == 2 * nθ * nφ - 2 * nφ == length(g.fills)
    @test all(abs(sqrt(sum(abs2, v)) - 50.0) < 1e-9 for v in g.vertices)
    @test all(1 <= i <= length(g.vertices) for f in g.faces for i in f) && issubset(g.fills, 1:5)
    fl = corridor_fills(m)
    @test count(==(3), g.fills) == 2 * count(==(3), fl[2:end-1, :]) + count(==(3), fl[[1, end], :])
    @test Set(g.fills) == Set(fl) && 1 in g.fills && 3 in g.fills && 5 in g.fills
    # a plane: one square cell about each sample, two triangles each, in the plane
    mp = corridor_map([b]; plane = ((0, 0, 3), (0, 0, 1)), extent = 5.0, n = 5)
    gp = corridor_surface(mp)
    @test length(gp.vertices) == 36 && length(gp.faces) == 50 == length(gp.fills)
    @test all(v -> abs(v[3] - 3) < 1e-12 && maximum(abs, v[1:2]) <= 6.25 + 1e-12, gp.vertices)
    @test gp.fills[1:2:end] == vec(permutedims(corridor_fills(mp)))
    @test corridor_extent(mp) ≈ 3 + sqrt(2) * 5 && corridor_extent(m) == 50
    # the drawing of the Clearance tab: the classes as digits, the first index fastest
    dm = corridor_drawing(m)
    @test startswith(dm, "{\"surface\":\"sphere\",\"nu\":$nθ,\"nv\":$nφ,") &&
          occursin("\"fills\":\"" * join(vec(fl)) * "\"", dm)
    @test occursin("\"surface\":\"plane\",\"nu\":5,\"nv\":5,\"extent\":5.0000,", corridor_drawing(mp))
    clear_corridor!(cs)
    @test cs.map[] === nothing && cs.status[] == "" && !corridor_stale(cs, rv, tl, [b], sl)

    # through a session: the enabled beams and the ticked diagnostics
    sess = DesignSession(resolution = :coarse)
    sess.time.diagnostics[] = ["nif_dante1"]
    lines = corridor_sight_lines(sess)
    @test length(lines) == 1 && lines[1].radius == 30.0 && lines[1].axis ≈ sight_line("nif_dante1").axis
    show_corridor!(sess)
    cm = sess.corridor
    @test cm.shown[] && cm.map[] isa CorridorMap && cm.map[].radius == 50.0
    @test count(s -> s[1] === :sight, cm.map[].sources) == 1 &&
          count(s -> s[1] === :omega3, cm.map[].sources) == length(active_beams(sess.beams))
    @test any(<=(0), cm.map[].omega1) && any(<=(0), cm.map[].sight) && any(>(0), cm.map[].margin)
    set_corridor_radius!(sess, 120.0)
    @test cm.map[].radius == 120.0 && !corridor_stale(cm, sess.review, sess.talis, active_beams(sess.beams), lines)
    # a ticked diagnostic that gives no line of sight is named in the status line
    @test isempty(corridor_omitted(sess)) && !occursin("left out", cm.status[])
    sess.time.diagnostics[] = ["nif_dante1", "no_such_camera"]
    @test length(corridor_sight_lines(sess)) == 1
    @test corridor_omitted(sess) == ["no_such_camera (not on the diagnostic search path)"]
    show_corridor!(sess)
    @test endswith(cm.status[], "; left out, no line of sight: no_such_camera (not on the diagnostic search path)")
    @test count(s -> s[1] === :sight, cm.map[].sources) == 1
    sess.time.diagnostics[] = String[]
    @test corridor_stale(cm, sess.review, sess.talis, active_beams(sess.beams), corridor_sight_lines(sess))
    @test_throws ArgumentError set_corridor_radius!(sess, 0.0)
    show_corridor!(sess, false)
    @test !cm.shown[] && cm.map[] === nothing
    set_corridor_radius!(sess, 80.0)
    @test cm.radius[] == 80.0 && cm.map[] === nothing
    @test_throws ArgumentError export_corridor_svg(sess, tempname())
    # a plane through the mounting zone, set from text; a shown map follows the settings
    @test edit_corridor!(sess, :point, "0, 0, 8") && edit_corridor!(sess, :normal, "0 0 2") &&
          edit_corridor!(sess, :extent, "6")
    @test cm.point[] == (0.0, 0.0, 8.0) && cm.normal[] == (0.0, 0.0, 2.0) && cm.extent[] == 6.0 && cm.error[] == ""
    @test !edit_corridor!(sess, :normal, "0, 0, 0") && occursin("must not be zero", cm.error[])
    @test !edit_corridor!(sess, :extent, "-1") && !edit_corridor!(sess, :point, "1, 2") && cm.extent[] == 6.0
    set_corridor_surface!(sess, :plane)
    show_corridor!(sess)
    ref = corridor_map(active_beams(sess.beams); sight_lines = corridor_sight_lines(sess), plane = ((0, 0, 8), (0, 0, 2)),
                       extent = 6.0, corridor_kwargs(sess.review, sess.talis)...)
    @test cm.map[].surface === :plane && cm.map[].margin == ref.margin
    @test startswith(cm.status[], "plane through (0, 0, 8) mm, normal (0, 0, 1), ±6 mm: ")
    @test edit_corridor!(sess, :extent, "4") && cm.map[].u[end] == 4.0
    @test !corridor_stale(cm, sess.review, sess.talis, active_beams(sess.beams), corridor_sight_lines(sess))
    path = export_corridor_svg(sess, joinpath(mktempdir(), "map.svg"))
    @test read(path, String) == corridor_map_svg(cm.map[]; tight = sess.talis.tight[])
    @test_throws ArgumentError set_corridor_surface!(sess, :cylinder)
    set_corridor_surface!(sess, :sphere)
    @test cm.map[].surface === :sphere && cm.map[].radius == 80.0
    show_corridor!(sess, false)
end

@testset "the view from a diagnostic (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    ds = diagnostic_sights(sess)
    ids = [d.id for d in ds]
    @test "nif_dante1" in ids && "nif_sxi_upper" in ids && !("omega_dante" in ids) && issorted([d.name for d in ds])
    v = look_from_diagnostic(sess, "nif_sxi_upper")
    @test v.line == sight_line("nif_sxi_upper")
    @test v.intrusions == sight_intrusions(sess.target.target[], v.line; tolerance = 0.025)
    # the hohlraum and capsule hold the center of the view and are observed; nothing else
    # enters it
    @test any(r -> r.part === :hohlraum && r.intrudes, v.intrusions) && Set(v.observed) == Set([:hohlraum, :capsule])
    @test isempty(v.hardware) && startswith(v.text, "looking from NIF SXI upper") &&
          occursin("observing ", v.text) && occursin("no other part within 0.025 mm", v.text)
    # a shield across the view, 8 mm toward the instrument, enters it by 1 mm and more
    a = v.line.axis
    c = 8.0 * a
    add_part!(sess.target, :disk; name = :shield, two_sided = true, x = c[1], y = c[2], z = c[3],
              theta = acosd(a[3]), phi = atand(a[2], a[1]), r1 = 1.0)
    w = look_from_diagnostic(sess, "nif_sxi_upper")
    @test [r.part for r in w.hardware] == [:shield] && w.hardware[1].margin < -1.0 &&
          occursin("shield (", w.text) && occursin(" enters the view within 0.025 mm", w.text)
    @test_throws ArgumentError look_from_diagnostic(sess, "nif_inflight_radiography")
end

@testset "beam rules of engagement (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    rs = sess.roe
    @test rs.report[] === nothing && isempty(roe_tiers(rs)) && roe_beam_text(rs, "Q11T") == "" && !roe_stale(sess)
    start_roe!(sess)
    r = rs.report[]
    @test r !== nothing && sess.view.beam_color[] === :roe && !roe_stale(sess)
    # the scan is the core scan of the scenario of the session with the default criteria,
    # each quad of the preset standing for its four beams: a unit for the quad and one
    # for a quarter of it
    sc = sess.drive.scenario[]
    beams = display_beams(sess.beams)
    held = roe_stands_for(sess.beams)
    @test length(held) == 48 && all(==(4), values(held))
    direct = roe_scan(sc, beams, roe_default_criteria(sc, beams).criteria; stands_for = held)
    @test r.tiers == direct.tiers && [c.name for c in r.criteria] == [:drive, :P2, :P4, :hotspot_velocity]
    @test count(u -> u.kind === :quad, r.units) == 48 && count(u -> u.kind === :beam && u.share == 0.25, r.units) == 48
    @test r.units[1].name == r.beams[1] * " (1 of 4)" && r.units[49].name == r.beams[1]
    # one beam of a quad is not critical where the loss of the quad is
    @test all(==(:red), r.tiers[49:96]) && !any(==(:red), r.tiers[1:48])
    @test occursin("48 single beams (one of four in a quad): 0 critical", rs.status[]) &&
          occursin("48 quads: 48 critical", rs.status[])
    tiers = roe_tiers(rs)
    @test Set(keys(tiers)) == Set(r.beams) && all(i -> tiers[r.beams[i]] === r.tiers[i], 1:48)
    i = argmax([o.scores[o.binding] for o in r.singles])
    u, o = r.units[i], r.singles[i]
    text = roe_beam_text(rs, u.name)
    @test occursin(string(r.tiers[i]), text) && occursin(string(r.criteria[o.binding].name), text) && occursin("amber ≤", text)
    @test roe_beam_text(rs, "nonesuch") == ""
    @test occursin("units against drive", rs.status[]) && occursin("engineering defaults", rs.status[])
    # The port map of the Outages tab, by quad, in the colors of the tiers of the quads
    # and of their single beams.
    ms = roe_marks(rs)
    @test [m.name for m in ms] == r.beams && all(m -> m.tier === :red && m.beam === tiers[m.name] && !m.picked, ms)
    pd = roe_port_drawing(rs)
    @test startswith(pd, "[{\"name\":\"" * ms[1].name * "\",") && count("\"tier\":", pd) == length(ms) &&
          count("\"beam\":", pd) == length(ms)
    # A hypothetical outage picked on the map: a pair of quads agrees with the linear pair
    # of the scan, and a third quad joins it.
    a, b, c = r.units[49].name, r.units[50].name, r.units[51].name
    @test roe_outage_text(rs) == "" && isempty(roe_outage_rows(rs))
    pick_outage!(sess, a)
    pick_outage!(sess, b)
    pair = only(p for p in r.pairs if Set(p.units) == Set((49, 50)))
    @test rs.picked[] == [a, b] && isapprox(rs.outage[].scores, pair.scores; rtol = 1e-9)
    @test count(m -> m.picked, roe_marks(rs)) == 2
    pick_outage!(sess, c)
    o3 = rs.outage[]
    @test o3.scores == roe_outage(r, sc, beams, [a, b, c]).scores && o3.beams == [1, 2, 3]
    rows = roe_outage_rows(rs)
    @test [x.criterion for x in rows] == ["drive", "P2", "P4", "hotspot_velocity"] && count(x -> x.binding, rows) == 1
    @test startswith(roe_outage_text(rs), "outage of $a, $b, $c (12 beams): $(o3.tier), bound by ")
    pick_outage!(sess, a)
    @test rs.picked[] == [b, c] && rs.outage[].beams == [2, 3]
    # one beam of a quad, as a quarter of its power
    clear_outage!(sess)
    pick_outage!(sess, r.units[1].name)
    @test rs.outage[].shares == [0.25] && startswith(roe_outage_text(rs), "outage of $(r.units[1].name) (1 beam): ")
    @test [m.name for m in roe_marks(rs) if m.picked] == [r.beams[1]]
    @test_throws ArgumentError pick_outage!(sess, "nonesuch")
    clear_outage!(sess)
    @test isempty(rs.picked[]) && rs.outage[] === nothing
    # a changed beam marks the scan out of date
    bs = sess.beams.beams[]
    set_beams!(sess.beams, [k == 1 ? remake(b; power = 0.5 * b.power) : b for (k, b) in enumerate(bs)])
    @test roe_stale(sess)
end
