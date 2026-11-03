# Hohlraum profiles (assemblies phase D): a hohlraum part's frustraum and section-list
# profiles against core `hohlraum(; profile)`, its wall pieces and their roles, the
# frustraum starter and its drive, the section table's commands, window frames on a
# section, and sessions with and without profiles.

@testset "hohlraum profiles: frustraum and sections against core (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    h = ts.mesh[].h
    i = part_index(ts, :hohlraum)
    wall = gold(albedo = 0.85, conversion = 0.8)
    core(prof) = hohlraum(; profile = prof, leh_radius = 1.55, capsule_radius = 0.0, wall, h,
                          leh_open_radius = 0.84 * 1.55)
    # a frustraum: the end radius filled in on the switch (60 % of the waist radius, at
    # least 1.2 times the LEH radius)
    @test edit_part!(ts, i, :profile, "1")
    p = ts.parts[][i]
    @test p.values[:end_radius] == 1.86                    # 1.2 times the LEH radius here
    @test hohlraum_pieces(p) == [:wall_1, :wall_2, :cap_top, :cap_bottom, :leh_top, :leh_bottom]
    ref = core(frustraum_profile(; length = 10.1, waist_radius = 2.875, end_radius = p.values[:end_radius]))
    cs = part_components(p, h)
    @test [c.name for c in cs] == [Symbol(:hohlraum_, c.name) for c in ref.components]
    @test all(_same_component(a, b) for (a, b) in zip(cs, ref.components))
    @test all(c -> c.part === :hohlraum, cs)
    @test wall_names(ts) == [:hohlraum_wall_1, :hohlraum_wall_2, :hohlraum_cap_top, :hohlraum_cap_bottom]
    @test closed_design(ts)
    @test edit_part!(ts, i, :bottom_end_radius, "2.0")
    t = ts.target[]
    @test sum(t.elem_area[elements(t, :hohlraum_cap_bottom, :hohlraum_leh_bottom)]) ≈ 4π rtol = 1e-12
    @test !edit_part!(ts, i, :end_radius, "1.5") && occursin("narrower than the end", ts.error[])
    @test !edit_part!(ts, i, :profile, "3")

    # sections: from the frustraum's breaks; a cylinder starts as four equal sections
    @test edit_part!(ts, i, :profile, "2")
    p = ts.parts[][i]
    @test p.sections == [(-5.05, 2.0), (0.0, 2.875), (5.05, 1.86)]
    @test (hohlraum_profile(p).z, hohlraum_profile(p).r) == ([-5.05, 0.0, 5.05], [2.0, 2.875, 1.86])
    ts2 = TargetSetup(resolution = :coarse)
    @test edit_part!(ts2, 1, :profile, "2")
    q = ts2.parts[][1]
    @test q.sections == [(z, 2.875) for z in (-5.05, -2.525, 0.0, 2.525, 5.05)]
    @test hohlraum_pieces(q)[1:4] == [:wall_1, :wall_2, :wall_3, :wall_4]
    t = ts2.target[]
    @test sum(t.elem_area[elements(t, (Symbol(:hohlraum_wall_, k) for k in 1:4)...)]) ≈ 2π * 2.875 * 10.1 rtol = 1e-12
    # the table: a radius, an inner height, a break added (the shape kept) and removed
    @test edit_hohlraum_break!(ts2, 1, 3, :r, "3.1")
    @test edit_hohlraum_break!(ts2, 1, 2, :z, "-3.0")
    @test ts2.parts[][1].sections[2:3] == [(-3.0, 2.875), (0.0, 3.1)]
    @test !edit_hohlraum_break!(ts2, 1, 1, :z, "-4.0") && occursin("length", ts2.error[])
    @test !edit_hohlraum_break!(ts2, 1, 2, :z, "0.5") && occursin("go up", ts2.error[])
    @test !edit_hohlraum_break!(ts2, 1, 2, :r, "-1")
    @test add_hohlraum_break!(ts2, 1)
    q = ts2.parts[][1]
    @test length(q.sections) == 6 && hohlraum_pieces(q)[5] === :wall_5
    @test remove_hohlraum_break!(ts2, 1, 2)
    @test length(ts2.parts[][1].sections) == 5
    @test !remove_hohlraum_break!(ts2, 1, 1)
    # the length moves the end breaks
    @test edit_part!(ts2, 1, :length, "11.0")
    @test hohlraum_profile(ts2.parts[][1]).z[[1, end]] == [-5.5, 5.5]
    # a window on a cylindrical section, its frame on that section's wall
    @test edit_part!(ts2, 1, :profile, "0") && edit_part!(ts2, 1, :profile, "2")   # (the breaks kept)
    @test length(ts2.parts[][1].sections) == 5
    ts3 = TargetSetup(resolution = :coarse)
    @test edit_part!(ts3, 1, :profile, "2") && edit_hohlraum_break!(ts3, 1, 3, :r, "3.0")
    @test !add_window!(ts3)                                       # z = 0 is a break between cones
    @test add_window!(ts3; phi = 90.0, z = -4.0)
    F = part_frames(ts3)[Symbol("hohlraum.window_1")]
    @test F.origin ≈ Briareus.Vec3(0, 2.875, -4.0) atol = 1e-12
end

@testset "frustraum starter, its drive and sessions (no GL)" begin
    ps = starter_parts(:frustraum)
    @test [q.name for q in ps] == [:frustraum, :capsule] && ps[2].observed && ps[2].frame === :frustraum
    sess = DesignSession(starter = :frustraum, resolution = :coarse)
    ts = sess.target
    @test observed_name(ts) === :capsule && closed_design(ts)
    t = ts.target[]
    zs = extrema(v[3] for c in t.components if c.name === :frustraum_wall_2 for v in c.mesh.vertices)
    @test zs == (0.0, 6.12)
    # the presets aim in its entrance holes
    @test Set(sess.beams.aim_frame[]) == Set([Symbol("frustraum.leh_top"), Symbol("frustraum.leh_bottom")])
    build_scenario!(sess)
    res = sess.drive.result[]
    @test !isempty(res.modes) && res.mean_flux > 0
    eb = energy_balance(res.solution)
    @test sum(values(eb.absorbed)) + sum(values(eb.escaped)) + eb.lost ≈ eb.source rtol = 1e-10
    # the keep-out sees the two entrance holes
    @test ApertureKeepOut(sess.drive.scenario[].target).names == [:frustraum_leh_top, :frustraum_leh_bottom]

    # sessions: a frustraum and a section list round trip; a session saved before
    # profiles (no profile values, no sections) loads as before
    s = DesignSession(resolution = :coarse)
    @test edit_part!(s.target, 1, :profile, "2") && edit_hohlraum_break!(s.target, 1, 3, :r, "3.0")
    path = joinpath(mktempdir(), "profile.briar.toml")
    save_session(path, s)
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, path)
    @test session_state(s2) == session_state(s)
    @test s2.target.parts[][1].sections == s.target.parts[][1].sections
    @test nelements(s2.target.target[]) == nelements(s.target.target[])
    s0 = DesignSession(resolution = :coarse)
    state = session_state(s0)
    for p in state["target"]["parts"]
        delete!(p, "sections")
        foreach(k -> delete!(p["values"], k), ("profile", "end_radius", "bottom_end_radius"))
    end
    s3 = DesignSession(resolution = :coarse)
    restore!(s3, state)
    @test session_state(s3) == session_state(s0)
    @test all(_same_component(a, b) for (a, b) in zip(s3.target.target[].components, s0.target.target[].components))
end

@testset "hohlraum breaks: a new break splits the longest section, the wall kept (no GL)" begin
    ts = TargetSetup(resolution = :coarse)
    i = part_index(ts, :hohlraum)
    @test edit_part!(ts, i, :profile, "1") && edit_part!(ts, i, :bottom_end_radius, "2.0")
    @test edit_part!(ts, i, :profile, "2")
    @test ts.parts[][i].sections == [(-5.05, 2.0), (0.0, 2.875), (5.05, 1.86)]
    wall(ts) = (t = ts.target[]; sum(t.elem_area[elements(t, (c.name for c in t.components
                                                              if startswith(String(c.name), "hohlraum_wall"))...)]))
    A0 = wall(ts)
    # the frustum sections' lateral areas, π (r₀ + r₁) × slant
    @test A0 ≈ π * (2.0 + 2.875) * hypot(5.05, 0.875) + π * (2.875 + 1.86) * hypot(5.05, 1.015) rtol = 1e-12
    # the upper section is the longer (slant 5.151 against 5.125): its midpoint is the new break
    @test add_hohlraum_break!(ts, i)
    p = ts.parts[][i]
    @test p.sections == [(-5.05, 2.0), (0.0, 2.875), (2.525, (2.875 + 1.86) / 2), (5.05, 1.86)]
    @test hohlraum_pieces(p)[1:3] == [:wall_1, :wall_2, :wall_3]
    @test wall(ts) ≈ A0 rtol = 1e-12
    @test hohlraum_profile(p).r[3] ≈ (2.875 + 1.86) / 2
    # again: now the lower section is the longest
    @test add_hohlraum_break!(ts, i)
    @test ts.parts[][i].sections[2] == (-2.525, (2.0 + 2.875) / 2) && isapprox(wall(ts), A0; rtol = 1e-12)
    # a part that is not a hohlraum has no sections
    @test !add_hohlraum_break!(ts, part_index(ts, :capsule)) && !isempty(ts.error[])
end
