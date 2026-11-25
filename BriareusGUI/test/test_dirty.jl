# The dirty domains of the window (qt/dirty.jl): the masks, the domains each observable of
# a session marks, and a refresh of the marked domains against a full refresh. A shell is
# built without a window, as in the precompile workload.

@testset "dirty domains: masks (no window)" begin
    B = BriareusGUI
    @test B.DIRTY_DOMAINS == (:tree, :selection, :beams, :drive, :time, :optimize, :review, :talis, :sweep, :view)
    @test B._dirty_mask() == B._DIRTY_EDIT != B._DIRTY_NONE && B._dirty_domains(B._dirty_mask()) == ()
    @test B._dirty_domains(B._DIRTY_ALL) == B.DIRTY_DOMAINS && B._dirty_mask(B.DIRTY_DOMAINS...) == B._DIRTY_ALL
    @test B._dirty_domains(B._dirty_mask(:view, :time)) == (:time, :view)
    @test all(d -> B._is_dirty(B._dirty_mask(d), d) && B._dirty_domains(B._dirty_mask(d)) == (d,), B.DIRTY_DOMAINS)
    @test_throws ArgumentError B._dirty_mask(:nowhere)
    r = Ref(B._DIRTY_NONE)
    B._mark!(r, B._dirty_mask(:beams))
    B._mark!(r, B._dirty_mask(:tree))
    @test B._dirty_domains(r[]) == (:tree, :beams)
    @test B._take_dirty!(r) == B._dirty_mask(:tree, :beams) && r[] == B._DIRTY_NONE
    B._mark!(r)
    @test r[] == B._DIRTY_ALL
end

@testset "dirty domains: a refresh of the marked domains matches a full refresh (no window)" begin
    B = BriareusGUI
    s = DesignSession(resolution = :coarse)
    build_scenario!(s)
    use_pulse_shape!(s; step = 3.0)
    start_history!(s)
    sh = B.BriareusShell(s)
    B._SHELL[] = sh
    try
        # every observable the window shows marks its domains, once
        src = B._dirty_sources(s)
        @test allunique(objectid.(first.(src))) && all(m -> m & B._DIRTY_EDIT != 0, last.(src))
        mask_of(obs) = only(m for (o, m) in src if o === obs)
        @test mask_of(s.target.parts) == B._DIRTY_ALL && mask_of(s.view.selection) == B._DIRTY_ALL
        @test B._dirty_domains(mask_of(s.time.cursor)) == (:time,)
        @test B._dirty_domains(mask_of(s.view.cutaway)) == (:view,)
        @test B._dirty_domains(mask_of(s.review.progress)) == (:tree, :review, :view)

        select!(s.view, "comp:hohlraum_wall")
        s.view.expanded[] = Set(r.id for r in model_tree(s))
        B._refresh!(sh)
        B._take_dirty!(sh.dirty)
        # what QML shows: the properties (figures aside) and the rows of every model
        fields(r) = Any[getfield(r, f) for f in fieldnames(typeof(r))]
        snapshot() = (Dict{String,Any}(k => v for (k, v) in sh.shown if !(v isa Figure)),
                      Dict{Symbol,Any}(k => fields.(rows) for (k, rows) in sh.rows))
        # After `edit`, the domains marked, and whether refreshing them shows what a full
        # refresh shows.
        function narrow(edit)
            edit()
            mask = B._take_dirty!(sh.dirty)
            B._refresh!(sh, mask)
            a = snapshot()
            B._refresh!(sh)
            return B._dirty_domains(mask), mask != B._DIRTY_NONE && isequal(a, snapshot())
        end
        same(edit) = last(narrow(edit))
        # an edit of the design marks every domain
        @test narrow(() -> B.br_edit_mesh("h", "0.55")) == (B.DIRTY_DOMAINS, true)
        build_scenario!(s)
        start_history!(s)
        B._refresh!(sh)
        B._take_dirty!(sh.dirty)
        # The first step of the cursor colors the view by a field of the history, which
        # sets the field of the drive too (every domain); the next ones only move the cursor.
        @test same(() -> B.br_time_step(2)) && time_label(s) != "" && s.view.color[] !== :component
        @test narrow(() -> B.br_time_step(3)) == ((:time,), true)
        @test narrow(() -> set_time_cursor!(s, s.time.history[].time[1])) == ((:time,), true)
        @test narrow(() -> B.br_time_chart(s.time.chart[] === :energy ? "pulses" : "energy")) == ((:time,), true)
        d, ok = narrow(() -> B.br_time_now())
        @test ok && !(:tree in d) && !(:selection in d)
        # the viewer's options
        @test narrow(() -> B.br_set_cutaway(!s.view.cutaway[])) == ((:view,), true)
        @test narrow(() -> B.br_cut_width(s.view.cut_width[] - 25.0)) == ((:view,), true)
        @test narrow(() -> B.br_cut_azimuth(s.view.cut_azimuth[] + 45.0)) == ((:view,), true)
        @test narrow(() -> B.br_cut_axis("x")) == ((:view,), true)
        @test narrow(() -> B.br_cut_kind("plane")) == ((:view,), true)
        @test narrow(() -> B.br_cut_offset(0.3)) == ((:view,), true)
        @test narrow(() -> B.br_show_off(!s.view.show_off[])) == ((:view,), true)
        @test all(st -> narrow(() -> B.br_beam_style(String(st))) == ((:view,), true), filter(!=(s.view.beams[]), BEAM_STYLES))
        @test narrow(() -> (s.view.beam_length[] = 7.5)) == ((:view,), true)
        @test narrow(() -> (s.view.beam_color[] = BEAM_COLORINGS[end])) == ((:view,), true)
        @test narrow(() -> (s.view.camera[] = (; eye = (9.0, 2.0, 4.0), lookat = (0.0, 0.0, 0.0), up = (0.0, 0.0, 1.0),
                                              ortho = true, turntable = false))) == ((:view,), true)
        # the coloring by component; a coloring by a field sets the field of the drive too
        @test narrow(() -> B.br_set_color("component")) == ((:drive, :time, :view), true)
        @test all(m -> same(() -> B.br_set_color(String(m))), COLOR_MODES)
        # the beam table's sort and the tree's expanded nodes
        @test narrow(() -> B.br_sort_beams(String(BEAM_SORT_COLUMNS[end]))) == ((:beams,), true)
        @test narrow(() -> B.br_toggle_expand("target")) == ((:tree,), true)
        # a job's progress: its panel, the tree and the running jobs
        @test narrow(() -> (s.review.progress[] = (3, 7))) == ((:tree, :review, :view), true)
        @test narrow(() -> (s.sweep.progress[] = (1, 5))) == ((:tree, :sweep, :view), true)
        @test narrow(() -> (s.time.progress[] = (2, 9))) == ((:tree, :time, :view), true)
        @test narrow(() -> (s.optimize.evaluations[] = 12)) == ((:tree, :optimize, :view), true)
        # with nothing marked, a tick refreshes nothing; a failing callback marks every domain
        @test sh.dirty[] == B._DIRTY_NONE
        B.br_time_chart("no such chart")
        @test sh.dirty[] == B._DIRTY_ALL && startswith(s.status[], "error")
    finally
        B._SHELL[] = nothing
    end
end
