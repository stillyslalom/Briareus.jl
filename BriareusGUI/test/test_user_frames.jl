# User frames: the Frames folder (add, rename, remove), parts placed and beams aimed in
# a user frame, port frames, the coordinate transform against hand computations, the
# drawn axes, and sessions (round trip, version 8).

using LinearAlgebra: I, norm

const U3 = Briareus.Vec3

@testset "user frames: add, place, aim, rename, remove (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts, bt = s.target, s.beams
    n0, sum0 = nelements(ts.target[]), target_summary(ts)
    # a frame at (1, 0, 0) with z along chamber +x (θ = 90°): x is chamber -z, y chamber +y
    @test add_frame!(ts; x = 1.0, theta = 90.0)
    @test [p.name for p in user_frames(ts)] == [:frame_1] && nelements(ts.target[]) == n0
    @test part_frames(ts)[:frame_1] ≈ Frame([0.0 0 1; 0 1 0; -1 0 0], U3(1, 0, 0))
    @test target_summary(ts) == sum0                               # not a part of the target
    # the Frames folder lists it once expanded
    rows = model_tree(s)
    @test only(r for r in rows if r.id == "frames").detail == "1"
    toggle_expanded!(s.view, "frames")
    @test only(r for r in model_tree(s) if r.id == "frame:frame_1").detail == "in TCC"
    # a part placed in it: 0.5 mm along its z is chamber (1.5, 0, 0)
    @test add_part!(ts, :disk; name = :shield, frame = :frame_1, z = 0.5, r1 = 0.3)
    @test part_frames(ts)[:shield].origin ≈ U3(1.5, 0, 0) atol = 1e-15
    # beams aimed in it, at local (0, 0.2, 0): chamber (1, 0.2, 0)
    idx = findall(==(:outer50), bt.cluster[])[1:2]
    set_aim_frame!(s, idx, :frame_1)
    edit_beam_set!(bt, idx, :aim_x, "0")
    edit_beam_set!(bt, idx, :aim_y, "0.2")
    edit_beam_set!(bt, idx, :aim_z, "0")
    @test all(i -> bt.beams[][i].aim ≈ U3(1, 0.2, 0), idx)
    # moved: the part and the beams follow exactly
    i = part_index(ts, :frame_1)
    @test edit_part!(ts, i, :x, "1.25")
    @test all(i -> isapprox(bt.beams[][i].aim, U3(1.25, 0.2, 0); atol = 1e-15), idx)
    @test part_frames(ts)[:shield].origin ≈ U3(1.75, 0, 0) atol = 1e-15
    # renamed: the part placed in it, the beams aimed in it and its drawn axes follow
    set_triad!(s.view, :frame_1, true)
    @test edit_part!(ts, i, :name, "mount")
    @test ts.parts[][part_index(ts, :shield)].frame === :mount && all(==(:mount), bt.aim_frame[][idx])
    @test s.view.triads[] == Set([:mount]) && drawn_triads(s) == [:mount]
    @test !edit_part!(ts, i, :name, "shield") && occursin("defined twice", ts.error[])
    # a frame has no surface fields and no role
    @test !edit_part!(ts, i, :albedo, "0.5") && occursin("no field", ts.error[])
    @test !set_part_role!(ts, i, :foil)
    # removal is refused while a part is placed in it, then allowed; its beams keep their aims
    @test frame_users(ts, :mount) == [:shield]
    @test !remove_part!(ts, i) && occursin("in use", ts.error[]) && occursin("shield", ts.error[])
    @test remove_part!(ts, part_index(ts, :shield))
    a = [bt.beams[][i].aim for i in idx]
    @test remove_part!(ts, part_index(ts, :mount))
    @test isempty(user_frames(ts)) && all(==(:tcc), bt.aim_frame[][idx]) && [bt.beams[][i].aim for i in idx] == a
    @test drawn_triads(s) == Symbol[]
    # a frame in another frame (a hohlraum piece's), and the selection's axes drawn
    @test add_frame!(ts; name = :lip, frame = Symbol("hohlraum.leh_top"), z = 0.4)
    L = ts.parts[][part_index(ts, :hohlraum)].values[:length]
    @test part_frames(ts)[:lip].origin ≈ U3(0, 0, L / 2 + 0.4) atol = 1e-14
    select!(s.view, "frame:lip")
    @test drawn_triads(s) == [:lip]
end

@testset "port frames (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts, bt = s.target, s.beams
    ports = port_choices(bt)
    @test length(ports) == 48
    q = only(p for p in ports if p.name == "Q31B")
    @test (q.theta, q.phi) == (150.0, 236.25)                     # the NIF port table's port 67
    @test add_port_frame!(ts, q.name, q.theta, q.phi)
    p = ts.parts[][end]
    @test p.name === :port_Q31B && is_port_frame(p) && p.object == "Q31B"
    @test part_frames(ts)[:port_Q31B] == port_frame(:nif, "Q31B")
    i = part_index(ts, :port_Q31B)
    # its origin, axis, frame and steps are its port's
    @test !edit_part!(ts, i, :theta, "10") && occursin("port frame", ts.error[])
    @test !edit_part!(ts, i, :x, "1") && occursin("port frame", ts.error[])
    @test !set_part_frame!(ts, i, :hohlraum) && occursin("port frame", ts.error[])
    @test !add_part_step!(ts, i, :offset) && occursin("port frame", ts.error[])
    @test !nudge_part!(ts, i, 1.0) && occursin("port frame", ts.error[])
    # the roll turns x about the port axis
    @test edit_part!(ts, i, :rotation, "30")
    @test part_frames(ts)[:port_Q31B] == port_frame(:nif, "Q31B"; roll = 30.0)
    # a frame 5 mm out along the port, and the port's beam aimed on its axis there
    @test add_frame!(ts; name = :out5, frame = :port_Q31B, z = 5.0)
    @test part_frames(ts)[:out5].origin ≈ 5 * port_vector(150.0, 236.25) atol = 1e-14
    b = findfirst(x -> x.quad == "Q31B", bt.beams[])
    set_aim_frame!(s, [b], :port_Q31B)
    edit_beam_set!(bt, [b], :aim_x, "0")
    edit_beam_set!(bt, [b], :aim_y, "0")
    edit_beam_set!(bt, [b], :aim_z, "2")
    @test bt.beams[][b].aim ≈ 2 * port_vector(150.0, 236.25) atol = 1e-14
    @test occursin("port Q31B · roll 30°", only(r for r in (toggle_expanded!(s.view, "frames"); model_tree(s))
                                               if r.id == "frame:port_Q31B").detail)
    # beams without quads (OMEGA): each beam's port, by the beam's name, at its angles
    o = DesignSession(starter = :foil, resolution = :coarse)
    bo = o.beams.beams[]
    @test all(c -> (c.theta, c.phi) == port_angles(bo[findfirst(x -> x.name == c.name, bo)]), port_choices(o.beams))
end

@testset "coordinate transform (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts = s.target
    # the frame at (1, 0, 0) with z = +x, x = -z, y = +y
    @test add_frame!(ts; name = :fx, x = 1.0, theta = 90.0)
    # by hand: (2, 3, 4) - (1, 0, 0) = (1, 3, 4) → (x·d, y·d, z·d) = (-4, 3, 1)
    @test all(isapprox.(transform_coordinates(ts, :tcc, :fx, (2, 3, 4)), (-4.0, 3.0, 1.0); atol = 1e-15))
    # spherical: r = √26, θ = acos(1/√26), φ = 180° - atan(3/4)
    @test all(isapprox.(transform_coordinates(ts, :tcc, :fx, (2, 3, 4); output = :spherical),
                        (√26, acosd(1 / √26), 180 - atand(3 / 4)); atol = 1e-12))
    # cylindrical in, Cartesian out: r = 2, φ = 90°, z = 1 in fx is (0, 2, 1) local, (2, 2, 0) chamber
    @test all(isapprox.(transform_coordinates(ts, :fx, :tcc, (2, 90, 1); input = :cylindrical), (2.0, 2.0, 0.0);
                        atol = 1e-15))
    # a direction: fx's z is chamber +x, whatever the origins
    @test all(isapprox.(transform_coordinates(ts, :fx, :tcc, (0, 0, 1); direction = true), (1.0, 0.0, 0.0);
                        atol = 1e-15))
    @test_throws ArgumentError transform_coordinates(ts, :tcc, :nowhere, (0, 0, 0))
    # a frame from a point: in fx, at the point, with TCC's axes
    @test add_frame_from_transform!(ts, :tcc, :fx, (2, 3, 4); name = :at)
    p = ts.parts[][part_index(ts, :at)]
    @test p.frame === :fx && part_frames(ts)[:at] ≈ Frame(Matrix(1.0I, 3, 3), U3(2, 3, 4))
    # from a direction: in fx, at its origin, z along chamber +y
    @test add_frame_from_transform!(ts, :tcc, :fx, (0, 1, 0); direction = true, name = :up)
    F = part_frames(ts)[:up]
    @test F.origin ≈ U3(1, 0, 0) && F.R[:, 3] ≈ U3(0, 1, 0)
    @test !add_frame_from_transform!(ts, :tcc, :fx, (0, 0, 0); direction = true) && occursin("nonzero", ts.error[])
end

@testset "user frames in sessions (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts = s.target
    @test add_frame!(ts; name = :mount, x = 2.0, phi = 30.0, theta = 45.0, rotation = 10.0)
    @test add_port_frame!(ts, "Q31B", port_angles(:nif, "Q31B")...; roll = 15.0)
    @test add_part!(ts, :disk; name = :shield, frame = :mount, z = 0.5, r1 = 0.3)
    set_aim_frame!(s, [1, 2], :port_Q31B)
    set_triad!(s.view, :mount, true)
    f = joinpath(mktempdir(), "frames.briar.toml")
    save_session(f, s)
    saved = TOML.parsefile(f)
    @test saved["version"] == 9 && saved["view"]["triads"] == ["mount"]
    @test count(p -> p["kind"] == "frame", saved["target"]["parts"]) == 2
    s2 = DesignSession(resolution = :coarse)
    load_session!(s2, f)
    @test session_state(s2) == session_state(s)
    @test part_frames(s2.target) == part_frames(ts) && s2.beams.aim_frame[][1:2] == [:port_Q31B, :port_Q31B]
    @test is_port_frame(s2.target.parts[][part_index(s2.target, :port_Q31B)])
    # a version 8 session (no user frames, no drawn axes) loads as it is
    v = DesignSession(resolution = :coarse)
    st = session_state(v)
    old = deepcopy(st)
    old["version"] = 8
    delete!(old["view"], "triads")
    v8 = DesignSession(resolution = :coarse)
    @test add_frame!(v8.target)                                    # replaced by the loaded design
    set_triad!(v8.view, :frame_1, true)
    restore!(v8, old)
    @test session_state(v8) == st && isempty(user_frames(v8.target)) && isempty(v8.view.triads[])
end
