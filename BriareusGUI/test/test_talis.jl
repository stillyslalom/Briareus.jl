# The TaLIS checks: a small foil with one beam and a plate upstream
# (the core's analytic case), then the coarse hohlraum with the NIF quads and a
# diagnostic window, through a DesignSession.

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
    @test n0[1].kind === :assumed && occursin("ASSUMED", n0[1].text) && occursin("TCC", n0[1].text)
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
    @test notes[1].kind === :assumed && occursin("ASSUMED", notes[1].text)
    @test any(n -> n.kind === :skipped && occursin("1ω", n.text), notes)    # not a NIF beam
    @test count(n -> n.kind === :scope, notes) == length(TALIS_OUT_OF_SCOPE)
    @test notes[end].kind === :source && occursin(TALIS_SOURCE, notes[end].text)

    # the tolerance interpretation reaches talis_check: no tolerance removes r·tan(0.5°)
    # at the plate's edge point (distance from TCC hypot(x0, 2))
    @test edit_talis!(tl, :tolerance, "0") && tl.tolerance[] == 0
    start_talis!(tl, rv, t, [b]; spawn = false)
    c0 = only(x for x in tl.report[].rows if x.rule === :clearance_3w && x.component === :plate)
    @test c0.worst - c3.worst ≈ tand(0.5) * hypot(x0, 2.0) atol = 1e-9
    @test edit_talis!(tl, :tolerance, "0.5")
    # a pivot at the edge point gives no allowance there
    @test edit_talis!(tl, :pivot_point, "1.6, 0, 2")
    set_talis_pivot!(tl, :point)
    @test talis_kwargs(tl, rv).pivot == (1.6, 0.0, 2.0)
    start_talis!(tl, rv, t, [b]; spawn = false)
    cP = only(x for x in tl.report[].rows if x.rule === :clearance_3w && x.component === :plate)
    @test cP.worst ≈ c0.worst atol = 1e-9
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

    path = export_talis(tl, joinpath(mktempdir(), "talis.csv"))
    @test countlines(path) == 1 + length(tl.report[].rows) && tl.status[] == "wrote talis.csv"
    @test_throws ArgumentError export_talis(TaLISState(), path)
    @test !talis_stale(tl, t, [b]) && talis_stale(tl, t, [remake(b; aim = (0.1, 0.0, 0.0))])

    # no capture component, no beams
    tl2 = TaLISState()
    default_selection!(tl2, t)
    set_talis_role!(tl2, :f, :none)
    start_talis!(tl2, rv, t, [b]; spawn = false)
    @test occursin("capture", tl2.status[]) && tl2.report[] === nothing
    start_talis!(tl2, rv, t, Beam{Float64}[]; spawn = false)
    @test tl2.status[] == "no beams are enabled"

    # on a worker: the result arrives through the queue; a cancel drops it
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

@testset "TaLIS checks: hohlraum with a window (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    @test add_window!(sess.target)                            # closed, at φ = 90°
    tl = sess.talis
    @test tl.capture[] == [:wall] && :window_1 in tl.hardware[] && :cap_top in tl.hardware[]
    @test !(:capsule in tl.hardware[]) && !(:capsule in tl.capture[])   # drive physics, not hardware
    @test talis_kwargs(tl, sess.review).leh == (:cap_top, :cap_bottom)
    @test :window_1 in sess.review.hardware[]                 # clearance checks it too
    @test task_status(sess, :talis)[1] === :todo
    start_talis!(sess; spawn = false)
    r = tl.report[]
    # the template's pointing grazes the LEH lips (an LEH-clearance review, not a
    # failure); apart from the closed window (hardware, in the 1% contours of beams
    # landing beside it), the only failures are the 23.5° quads' incidence on the wall,
    # below 70° at nominal and over 70° only with the 0.5° rule tolerance
    lip = [x for x in r.rows if x.rule === :clearance_3w && x.component in (:cap_top, :cap_bottom)]
    @test any(x -> x.status === :review, lip) && all(x -> x.status in (:pass, :review), lip)
    @test all(x -> occursin("LEH clearance", x.note), (x for x in lip if x.status === :review))
    bad = [x for x in r.rows if x.status === :fail && x.component !== :window_1]
    @test !isempty(bad) && all(x -> x.rule === :incidence_3w && x.component === :wall && x.nominal < 70, bad)
    @test any(x -> x.status === :fail && x.rule === :clearance_3w, (x for x in r.rows if x.component === :window_1))
    @test !any(x -> x.component === :capsule && x.rule !== :tas_stay_in, r.rows)   # stay-in: every part
    rules = Set(x.rule for x in r.rows)
    @test :tas_stay_in in rules && :clearance_1w in rules && :clearance_3w in rules
    @test any(x -> x.component === :window_1, r.rows)
    @test length(unique(x.beam for x in r.rows if x.rule !== :tas_stay_in)) == 48
    @test r.settings.capture == [:wall]
    # every quad lands on the wall: capture rows have no free substrate edge to measure
    @test all(x -> x.component === :wall, (x for x in r.rows if x.rule === :capture_3w))
    @test occursin("48 beams", talis_brief(tl))
    st, txt = task_status(sess, :talis)
    @test st === :attention && occursin(r"\d+ fail, \d+ review", txt)  # LEH reviews, incidence fails
    # roles follow the target: a removed window leaves the selection, a new one is checked
    set_talis_role!(tl, :capsule, :hardware)
    edit_target!(sess.target, :length, "10.0")
    @test :capsule in tl.hardware[] && !(:capsule in tl.capture[])        # kept as chosen
    @test !(:capsule in talis_kwargs(tl, sess.review).leh)
    @test talis_stale(tl, sess.target.target[], active_beams(sess.beams))
    @test task_status(sess, :talis)[1] === :attention
    remove_window!(sess.target, 1)
    @test !(:window_1 in tl.hardware[]) && !(:window_1 in sess.review.hardware[])
    add_window!(sess.target)
    @test :window_1 in tl.hardware[] && :window_1 in sess.review.hardware[]
    # switching to a foil chooses afresh
    set_target_kind!(sess.target, :foil)
    @test tl.capture[] == [FOIL_NAME] && tl.hardware[] == [:main]
    set_review_view!(sess.review, :talis)
    @test sess.review.view[] === :talis
    @test_throws ArgumentError set_review_view!(sess.review, :tally)
end
