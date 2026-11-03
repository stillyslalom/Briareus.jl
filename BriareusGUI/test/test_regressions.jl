# Fixes from the code survey (ROADMAP item 3): each analysis starts its own job, a renamed
# part keeps what was keyed by its name, undo keeps a current optimization result, the
# default cutaway follows the hohlraum parts, and one rule decides a flat observed face
# (no GL).

isdefined(@__MODULE__, :synthetic_workspace) || include(joinpath(pkgdir(Briareus), "test", "vrw_fixture.jl"))

@testset "Run starts each analysis's own job (no GL)" begin
    # with no beams every start refuses at once, on its own panel's status
    sess = DesignSession(resolution = :coarse, beams = Briareus.Beam[])
    status = Dict(:drive => sess.drive.status, :optimize => sess.optimize.status,
                  :clearance => sess.review.status, :sweep => sess.sweep.status,
                  :talis => sess.talis.status)
    @test Set(keys(status)) == Set(TASKS)
    changed = Dict{Symbol,Vector{Symbol}}()
    for task in TASKS
        foreach(o -> o[] = "", values(status))
        run_task!(sess, task; spawn = false)
        changed[task] = [k for (k, o) in status if !isempty(o[])]
    end
    @test all(task -> changed[task] == [task], TASKS)
    @test sess.view.beam_color[] === :talis          # the last check run colors the beams
    @test_throws ArgumentError run_task!(sess, :results)
end

@testset "a renamed part keeps what is keyed by its name (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    ts, v, dr, rv, tl, sw = sess.target, sess.view, sess.drive, sess.review, sess.talis, sess.sweep
    leh = Symbol("hohlraum.leh_top")
    @test set_excluded!(ts, :hohlraum_cap_top, true)
    set_display!(v, :hohlraum_wall, :wireframe)
    set_cavity!(dr, true; films = [:hohlraum_leh_bottom])
    @test set_drive_cavity!(sess, :hohlraum)
    set_role!(rv, :hohlraum_wall, :hardware)                  # not the default roles
    set_talis_role!(tl, :hohlraum_cap_bottom, :diagnostic)
    set_sweep!(sw; subject = :hohlraum, frame = leh, los_object = :hohlraum_wall)
    set_triad!(v, leh, true)
    select!(v, "comp:hohlraum_wall")
    record!(sess)
    @test edit_part!(ts, part_index(ts, :hohlraum), :name, "hr")
    # the cap stays out of the calculations, under its new name
    @test :hr_cap_top in ts.excluded[] && !(:hohlraum_cap_top in ts.excluded[])
    @test !(:hr_cap_top in component_names(ts.target[])) && :hr_cap_top in all_component_names(ts)
    @test v.display[] == Dict(:hr_wall => :wireframe)
    @test dr.films[] == Set([:hr_leh_bottom]) && dr.wall_cavity[] === :hr
    @test :hr_wall in rv.hardware[] && !(:hr_wall in rv.targets[]) && !any(n -> startswith(String(n), "hohlraum"), rv.targets[])
    @test tl.diagnostics[] == [:hr_cap_bottom] && :hr_wall in tl.capture[]
    @test sw.subject[] === :hr && sw.frame[] === Symbol("hr.leh_top") && sw.los_object[] === :hr_wall
    @test v.triads[] == Set([Symbol("hr.leh_top")]) && v.selection[] == "comp:hr_wall"
    @test v.cut_components[] == Set([:hr_wall, :hr_cap_top, :hr_cap_bottom])
    # the rename undone: the view keyed to the old names comes back with them (the view is
    # in the edit history), and follows the redo
    @test record!(sess) && undo!(sess)
    @test v.display[] == Dict(:hohlraum_wall => :wireframe) && v.triads[] == Set([leh])
    @test v.cut_components[] == Set([:hohlraum_wall, :hohlraum_cap_top, :hohlraum_cap_bottom])
    @test redo!(sess) && v.display[] == Dict(:hr_wall => :wireframe) && v.triads[] == Set([Symbol("hr.leh_top")])

    # a VisRad part's pieces carry its name as their part's: they follow it too
    dir = mktempdir()
    src = joinpath(dir, "work.vrw")
    write(src, synthetic_workspace())
    s = DesignSession(resolution = :coarse)
    import_visrad_workspace!(s, src)
    ts = s.target
    set_display!(s.view, :Cyl_Hohlraum_0_wall, :transparent)
    n0 = length(node_components(s, "part:Cyl_Hohlraum_0"))
    @test n0 > 1
    @test edit_part!(ts, part_index(ts, :Cyl_Hohlraum_0), :name, "hohl")
    cs = node_components(s, "part:hohl")
    @test length(cs) == n0 && :hohl_wall in cs && all(n -> startswith(String(n), "hohl_"), cs)
    @test s.view.display[] == Dict(:hohl_wall => :transparent)
    @test edit_part!(ts, part_index(ts, :hohl), :name, "cyl")
    @test length(node_components(s, "part:cyl")) == n0 && s.view.display[] == Dict(:cyl_wall => :transparent)
end

@testset "undo keeps a current optimization result (no GL)" begin
    sess = DesignSession(resolution = :coarse)
    st = sess.optimize
    build_scenario!(sess)
    st.iterations[] = 1
    run_task!(sess, :optimize; spawn = false)
    r = st.result[]
    @test r isa OptimizeResult && !optimize_stale(sess)
    record!(sess)
    # an edit the optimization does not depend on, undone: the result stays current
    @test edit_review!(sess.review, :peak_fraction, "0.001") && record!(sess)
    @test undo!(sess) && st.result[] === r && !optimize_stale(sess)
    @test task_status(sess, :optimize)[1] === :attention          # not applied
    # applied, then the apply undone and redone
    accept_optimum!(sess)
    @test record!(sess) && r.applied && task_status(sess, :optimize)[1] === :ok
    @test undo!(sess) && st.result[] === r && !r.applied && !optimize_stale(sess)
    @test redo!(sess) && r.applied && !optimize_stale(sess)
    # a design change: the result is out of date, kept, and current again on the way back
    @test edit_part!(sess.target, 1, :length, "10.2") && record!(sess)
    @test st.result[] === r && optimize_stale(sess) && task_status(sess, :optimize)[1] === :attention
    @test undo!(sess) && st.result[] === r && !optimize_stale(sess)
    @test redo!(sess) && st.result[] === r && optimize_stale(sess)
    # a session loaded from a file clears it
    path = joinpath(mktempdir(), "s.briar.toml")
    save_session(path, sess)
    load_session!(sess, path)
    @test st.result[] === nothing
end

@testset "the default cutaway cuts every hohlraum part's shell (no GL)" begin
    shells(ts, parts...) = Set(c.name for c in ts.target[].components if c.part in parts && c.kind === :surface)
    s = DesignSession(starter = :frustraum, resolution = :coarse)
    @test s.view.cut_components[] == default_cut(s.target) == shells(s.target, :frustraum)
    @test length(s.view.cut_components[]) > 3                     # sections and caps
    s = DesignSession(starter = :shock_tube, resolution = :coarse)
    @test s.view.cut_components[] == shells(s.target, :halfraum_top, :halfraum_bottom)
    @test length(s.view.cut_components[]) == 4                    # walls and outer caps
    s = DesignSession(resolution = :coarse)
    ts = s.target
    @test insert_starter!(ts, :hohlraum)
    @test s.view.cut_components[] == shells(ts, :hohlraum, :hohlraum_2) && :hohlraum_2_wall in s.view.cut_components[]
    @test edit_part!(ts, part_index(ts, :hohlraum), :name, "shell")
    @test s.view.cut_components[] == shells(ts, :shell, :hohlraum_2)
    # a cut chosen by hand stays as chosen
    set_cut!(s.view; components = [:capsule])
    @test insert_starter!(ts, :frustraum) && s.view.cut_components[] == Set([:capsule])
end

@testset "one rule for a flat observed face (no GL)" begin
    # the design's rule (on the built target) and the drive's (on its scenario) agree
    function agree(sess)
        build_scenario!(sess)
        s = sess.drive.scenario[]
        return !isempty(s.observed) && Controllers._observed_flat(s) == observed_flat(sess.target)
    end
    sess = DesignSession(resolution = :coarse)
    @test !observed_flat(sess.target) && agree(sess)                      # the capsule
    @test set_piece_role!(sess.target, :hohlraum_cap_top, :observed) && observed_flat(sess.target) && agree(sess)
    sess = DesignSession(starter = :opacity, resolution = :coarse)
    ts = sess.target
    i = part_index(ts, :sample)
    @test observed_flat(ts) && agree(sess)
    # a two-sided face faces both ways (no flat face): an observed face stays one-sided
    @test !set_part_flag!(ts, i, :two_sided, true) && occursin("one-sided", ts.error[])
    @test set_part_role!(ts, i, :none) && set_part_flag!(ts, i, :two_sided, true)
    @test !set_part_role!(ts, i, :observed) && occursin("one-sided", ts.error[])
    @test !Controllers._flat_face(ts.target[], elements(ts.target[], :sample))
    # a VisRad rectangle
    dir = mktempdir()
    src = joinpath(dir, "work.vrw")
    write(src, synthetic_workspace())
    s = DesignSession(resolution = :coarse)
    import_visrad_workspace!(s, src)
    @test set_part_role!(s.target, part_index(s.target, :Rectangle_0), :observed) && observed_flat(s.target)
end
