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
    # another target has no disposition, and the 1ω rules off leave none
    @test disposition_of(tl, Target(plate(:mirror, 4.0, 8.0 * ô))) === nothing
    _assign(tl.omega1, false)
    start_talis!(tl, rv, t, [bl]; spawn = false)
    @test tl.disposition[] === nothing && tl.report[] isa TaLISReport
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
