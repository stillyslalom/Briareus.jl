# The placement triad: a drag commits the Move-along row's edit (`nudge_part!`) or its
# rotational counterpart (`turn_part!`), checked against that edit and against a
# `Rotate` placement step; the 3D view's handles driven by synthetic mouse events
# through the camera's gesture handlers (offscreen GL).

using LinearAlgebra: norm

# Parts equal field by field (a `PartSpec` holds a Dict: `==` would compare identities).
_same_parts(a, b) = length(a) == length(b) &&
    all(isequal(getfield(p, f), getfield(q, f)) for (p, q) in zip(a, b) for f in fieldnames(PartSpec))

# A design with a tilted disk at x = 3 mm, a plate placed in its frame and four beams
# aimed in the disk's frame.
function _triad_design()
    s = DesignSession(starter = :none, resolution = :coarse)
    ts = s.target
    @assert add_part!(ts, :disk; name = :shield, x = 3.0, theta = 90, phi = 20, r1 = 0.5)
    @assert add_part!(ts, :rectangle; name = :plate, frame = :shield, z = 0.4, theta = 30, r1 = 0.3, r2 = 0.3)
    set_aim_frame!(s, 1:4, :shield; keep = false)
    edit_beam_set!(s.beams, 1:4, :aim_y, "0.2")
    return s
end

@testset "placement triad: edits (no GL)" begin
    a, b = _triad_design(), _triad_design()
    i = part_index(a.target, :shield)
    @test a.target.parts[][i].frame === :tcc
    # an arrow drag of three 0.1 mm steps is the Move-along row's 0.3 mm move, exactly
    Δ = triad_snap(0.287, 0.1)
    @test Δ === 0.3 && triad_snap(0.04, 0.1) === 0.0 && triad_snap(-0.26, 0.05) === -0.25
    @test commit_triad!(a.target, i, (:move, 2), Δ)
    @test nudge_part!(b.target, i, 0.3; frame = :tcc, axis = :y)
    @test _same_parts(a.target.parts[], b.target.parts[])
    # the beams aimed in the moved frame follow, as the Move-along row's move takes them
    @test a.beams.beams[] == b.beams.beams[] && a.beams.aim_frame[] == b.beams.aim_frame[]
    aims0 = [x.aim for x in _triad_design().beams.beams[][1:4]]
    @test all(isapprox(x.aim, y + Briareus.Vec3(0, 0.3, 0); atol = 1e-12) for (x, y) in zip(a.beams.beams[][1:4], aims0))
    # the preview places the parts and the beams where the edit puts them (to round-off)
    c = _triad_design()
    parts = triad_parts(c.target.parts[], i, (:move, 2), 0.3)
    @test _same_parts(parts, b.target.parts[])
    pv = triad_preview(c, parts)
    @test sort(unique(x.name for x in pv.components)) == [:plate, :shield]
    @test first.(pv.aims) == 1:4
    @test all(isapprox(aim, b.beams.beams[][k].aim; atol = 1e-12) for (k, aim) in pv.aims)
    byname(t, n) = t.components[findfirst(x -> x.name === n, t.components)]
    moved = Dict(x.name => x for x in pv.components)
    @test all(isapprox(u, v; atol = 1e-12) for n in (:plate, :shield)
              for (u, v) in zip(moved[n].mesh.vertices, byname(b.target.target[], n).mesh.vertices))
    # a zero drag makes no edit
    @test commit_triad!(c.target, i, (:turn, 1), 0.0) && _same_parts(c.target.parts[], _triad_design().target.parts[])

    # a ring drag of θ about axis k equals a `Rotate` step about the same chamber axis
    θ, k = 30.0, 1
    d, r = _triad_design(), _triad_design()
    fr = part_frames(d.target)
    e = fr[:tcc].R[:, k]                                       # the "Placed in" frame's x
    @test commit_triad!(d.target, i, (:turn, k), θ)
    ps = copy(r.target.parts[])
    ps[i] = Controllers._with(ps[i]; steps = Briareus.PlacementStep[Briareus.Rotate(fr[:shield].R' * e, θ)])
    @test Controllers._try_build!(r.target; parts = ps)
    Fd, Fr = part_frames(d.target), part_frames(r.target)
    @test all(n -> isapprox(Fd[n].R, Fr[n].R; atol = 1e-12) && isapprox(Fd[n].origin, Fr[n].origin; atol = 1e-12),
              [:shield, :plate])
    @test Fd[:shield].origin == fr[:shield].origin                       # it turns about its origin
    @test d.target.parts[][i].values[:x] == 3.0                          # its center as typed
    @test all(isapprox(u.aim, v.aim; atol = 1e-12) for (u, v) in zip(d.beams.beams[], r.beams.beams[]))
    # the triad's commit is turn_part! itself
    t2 = _triad_design()
    @test turn_part!(t2.target, i, θ; frame = :tcc, axis = :x) && _same_parts(t2.target.parts[], d.target.parts[])
    # a part placed in a tilted frame turns about that frame's axes, through its own origin
    j = part_index(d.target, :plate)
    F0 = part_frames(d.target)
    @test commit_triad!(d.target, j, (:turn, 3), -45.0)
    F1 = part_frames(d.target)
    Q = Briareus._axis_rotation(F0[:shield].R[:, 3], deg2rad(-45))
    @test isapprox(F1[:plate].R, Q * F0[:plate].R; atol = 1e-12) && isapprox(F1[:plate].origin, F0[:plate].origin; atol = 1e-12)
    @test triad_frame(d.target, j).R == F1[:shield].R && triad_frame(d.target, j).origin == F1[:plate].origin

    # one drag, one undo: undo returns the design exactly
    u = _triad_design()
    record!(u)
    parts0, beams0 = u.target.parts[], u.beams.beams[]
    n0 = length(u.history.undo)
    @test commit_triad!(u.target, i, (:turn, 2), 15.0)
    @test record!(u) && length(u.history.undo) == n0 + 1
    @test !record!(u)
    @test undo!(u)
    @test _same_parts(u.target.parts[], parts0) && u.beams.beams[] == beams0

    # a square's drag: one nudge along the sum of its two components, as the two
    # Move-along moves make it (to round-off), its preview the edit's parts, one undo
    q, w = _triad_design(), _triad_design()
    record!(q)
    n0 = length(q.history.undo)
    pq = triad_parts(q.target.parts[], i, (:plane, 3), (0.3, -0.2))
    @test commit_triad!(q.target, i, (:plane, 3), (0.3, -0.2)) && _same_parts(q.target.parts[], pq)
    @test nudge_part!(w.target, i, 0.3; frame = :tcc, axis = :x) && nudge_part!(w.target, i, -0.2; frame = :tcc, axis = :y)
    Fq, Fw = part_frames(q.target), part_frames(w.target)
    @test all(n -> isapprox(Fq[n].origin, Fw[n].origin; atol = 1e-12) && Fq[n].R == Fw[n].R, [:shield, :plate])
    @test all(isapprox(u.aim, v.aim; atol = 1e-12) for (u, v) in zip(q.beams.beams[], w.beams.beams[]))
    @test record!(q) && length(q.history.undo) == n0 + 1 && undo!(q) &&
          _same_parts(q.target.parts[], _triad_design().target.parts[])
    @test commit_triad!(q.target, i, (:plane, 1), (0.0, 0.0)) && _same_parts(q.target.parts[], _triad_design().target.parts[])

    # refused: a missing frame; a port frame (its axis is its port's)
    @test !turn_part!(u.target, i, 1.0; frame = :nowhere) && occursin("no frame", u.target.error[])
    @test add_port_frame!(u.target, "Q31T", 23.5, 78.75)
    pf = findfirst(is_port_frame, u.target.parts[])
    @test !turn_part!(u.target, pf, 10.0) && occursin("port", u.target.error[])
end

@testset "placement triad: the 3D view's handles (offscreen)" begin
    sess = DesignSession(resolution = :coarse)
    sc = B.scene_canvas(sess)
    c = sc.camera
    nplots = length(sc.lscene.scene.plots)
    V(x) = B.Makie.Vec3d(x...)
    # hidden with nothing selected and with a beam selected; on a piece, its part's
    select!(sess.view, "")
    @test sc.triad_at[] === nothing && !sc.place_arrows.visible[] && !sc.place_rings.visible[]
    select!(sess.view, beam_id(sess.beams.beams[][1]))
    @test sc.triad_at[] === nothing && !sc.place_arrows.visible[]
    select!(sess.view, "comp:hohlraum_wall")
    @test sc.triad_at[] !== nothing && sc.triad_at[][2] === :hohlraum && sc.place_arrows.visible[]
    select!(sess.view, "part:hohlraum")
    i, _, F = sc.triad_at[]
    @test F == triad_frame(sess.target, i)
    B.orient_camera!(c, :front)                    # looking along +y: x to the right, z up
    L = B._triad_length(c, F.origin)
    O, ex, ez = V(F.origin), V(F.R[:, 1]), V(F.R[:, 3])
    record!(sess)
    parts0, beams0, n0 = sess.target.parts[], sess.beams.beams[], length(sess.history.undo)
    nfollow = count(f -> startswith(String(f), "hohlraum"), sess.beams.aim_frame[])
    @test nfollow > 0                              # the preset aims in the hohlraum's frames

    # a click on a handle (no drag) is a click: no edit
    p_z = B.project_point(c, O + 0.85L * ez)
    @test B.triad_handle_at(c, F, L, p_z) == (:move, 3)
    c.last_click = -Inf
    B.camera_press!(c, p_z, :left)
    @test c.mode === :handle_pending && sc.placing[] !== nothing
    B.camera_release!(c, p_z)
    @test sc.placing[] === nothing && sess.target.parts[] === parts0
    select!(sess.view, "part:hohlraum")

    # the z arrow dragged 0.31 mm: the preview, snapped to three 0.1 mm steps, edits nothing
    B.camera_press!(c, p_z, :left)
    B.camera_move!(c, B.project_point(c, O + (0.85L + 0.12) * ez))
    @test c.mode === :handle && sc.placing[].amount == 0.1
    p_end = B.project_point(c, O + (0.85L + 0.31) * ez)
    B.camera_move!(c, p_end)
    g = sc.placing[]
    @test g.amount === 0.3 && g.handle == (:move, 3)
    @test sc.place_ghost.visible[] && sc.place_beams.visible[] && sc.place_label.visible[]
    @test sc.place_label.text[] == ["+0.3 mm"]
    @test length(sc.place_beams[1][]) == 2nfollow
    @test sess.target.parts[] === parts0 && !record!(sess)      # a preview is not an edit
    @test !isempty(colorbuffer(sc.fig; px_per_unit = 1))
    # released: one edit, the Move-along row's 0.3 mm along z of the frame it is placed in
    B.camera_release!(c, p_end)
    @test sc.placing[] === nothing && !sc.place_ghost.visible[]
    ref = DesignSession(resolution = :coarse)
    @test nudge_part!(ref.target, i, 0.3; frame = :tcc, axis = :z)
    @test _same_parts(sess.target.parts[], ref.target.parts[])
    @test sess.beams.beams[] == ref.beams.beams[]               # the followed beams, exactly
    @test sess.beams.beams[] != beams0
    @test record!(sess) && length(sess.history.undo) == n0 + 1 && !record!(sess)
    @test sc.triad_at[][3].origin ≈ F.origin + 0.3 * F.R[:, 3]  # the triad moved with it
    @test undo!(sess)
    @test _same_parts(sess.target.parts[], parts0) && sess.beams.beams[] == beams0

    # the y ring (facing the view) dragged 31° clockwise on screen: +30° about y
    select!(sess.view, "part:hohlraum")
    F = sc.triad_at[][3]
    ρ = B._TRIAD_RING * L
    at(t) = B.project_point(c, O + ρ * (cos(t) * ez + sin(t) * ex))
    @test B.triad_handle_at(c, F, L, at(π / 4)) == (:turn, 2)
    B.camera_press!(c, at(π / 4), :left)
    B.camera_move!(c, at(π / 4 + deg2rad(20)))
    @test sc.placing[].amount == 15.0
    B.camera_move!(c, at(π / 4 + deg2rad(31)))
    @test sc.placing[].amount == 30.0 && sc.place_label.text[] == ["+30°"]
    B.camera_release!(c, at(π / 4 + deg2rad(31)))
    rot = DesignSession(resolution = :coarse)
    @test turn_part!(rot.target, i, 30.0; frame = :tcc, axis = :y)
    @test _same_parts(sess.target.parts[], rot.target.parts[]) && sess.beams.beams[] == rot.beams.beams[]
    # and the Rotate step's frames and beams, to round-off
    stp = DesignSession(resolution = :coarse)
    ps = copy(stp.target.parts[])
    ps[i] = Controllers._with(ps[i]; steps = Briareus.PlacementStep[Briareus.Rotate(F.R' * [0, 1, 0], 30.0)])
    @test Controllers._try_build!(stp.target; parts = ps)
    Fa, Fb = part_frames(sess.target), part_frames(stp.target)
    @test all(n -> isapprox(Fa[n].R, Fb[n].R; atol = 1e-12) && isapprox(Fa[n].origin, Fb[n].origin; atol = 1e-12),
              collect(keys(Fb)))
    @test all(isapprox(u.aim, v.aim; atol = 1e-12) for (u, v) in zip(sess.beams.beams[], stp.beams.beams[]))
    @test undo!(sess) && _same_parts(sess.target.parts[], parts0)

    # Ctrl at the press: a tenth of the step; focus lost mid-drag: nothing changes
    select!(sess.view, "part:hohlraum")
    B.camera_press!(c, p_z, :left; ctrl = true)
    B.camera_move!(c, B.project_point(c, O + (0.85L + 0.123) * ez))
    @test sc.placing[].step ≈ 0.01 && sc.placing[].amount === 0.12
    B.release_input!(sc)
    @test sc.placing[] === nothing && c.mode === :none && _same_parts(sess.target.parts[], parts0)
    # the arrows snap to the Move-along row's step
    sc.triad_step[] = 0.25
    B.camera_press!(c, p_z, :left)
    B.camera_move!(c, B.project_point(c, O + (0.85L + 0.4) * ez))
    @test sc.placing[].amount === 0.5
    B.release_input!(sc)
    sc.triad_step[] = 0.1

    # hover: the arrow under the pointer drawn brighter and thicker, plain again off it
    arrow(k) = (6k-5):6k                           # its points in `place_arrows` (handle order)
    square(k) = 6 * B._TRIAD_RING_SEGMENTS .+ ((8k-7):8k)     # in `place_rings`, after the rings
    col(p, r) = p.color[][r]
    wid(p, r) = p.linewidth[][r]
    plain = B.TRIAD_COLORS[3]
    @test c.mode === :none && sc.triad_hover[] === nothing
    B.camera_move!(c, p_z)
    @test sc.triad_hover[] == (:move, 3)
    @test all(==(B._triad_hover_color(plain)), col(sc.place_arrows, arrow(3))) &&
          all(==(B._TRIAD_ARROW_WIDTH + B._TRIAD_HOT_WIDTH), wid(sc.place_arrows, arrow(3)))
    @test all(==(B.TRIAD_COLORS[1]), col(sc.place_arrows, arrow(1))) &&
          all(==(B._TRIAD_ARROW_WIDTH), wid(sc.place_arrows, arrow(1)))
    B.camera_move!(c, (5.0, 5.0))                  # off every handle
    @test sc.triad_hover[] === nothing && all(==(plain), col(sc.place_arrows, arrow(3))) &&
          all(==(B._TRIAD_ARROW_WIDTH), wid(sc.place_arrows, arrow(3)))
    # a drag keeps the dragged arrow highlighted, the pointer off its line included
    B.camera_move!(c, p_z)
    B.camera_press!(c, p_z, :left)
    B.camera_move!(c, B.project_point(c, O + (0.85L + 0.2) * ez + 0.3L * ex))
    @test c.mode === :handle && sc.triad_hover[] == (:move, 3)
    @test all(==(B._TRIAD_ACTIVE), col(sc.place_arrows, arrow(3))) &&
          all(==(B._TRIAD_ARROW_WIDTH + B._TRIAD_HOT_WIDTH), wid(sc.place_arrows, arrow(3)))
    B.release_input!(sc)
    @test _same_parts(sess.target.parts[], parts0)

    # the square of the zx plane (normal y, facing the view): hovered, then dragged 0.31 mm
    # along z and 0.22 mm along x, snapped to (0.3, 0.2): one move, one undo
    F = sc.triad_at[][3]
    mid = O + sum(B._TRIAD_PLANE) / 2 * L * (ez + ex)
    p_sq = B.project_point(c, mid)
    @test B.triad_handle_at(c, F, L, p_sq) == (:plane, 2)
    B.camera_move!(c, p_sq)
    @test sc.triad_hover[] == (:plane, 2) &&
          all(==(B._TRIAD_RING_WIDTH + B._TRIAD_HOT_WIDTH), wid(sc.place_rings, square(2))) &&
          all(==(B._TRIAD_RING_WIDTH), wid(sc.place_rings, square(3)))
    record!(sess)
    n1 = length(sess.history.undo)
    B.camera_press!(c, p_sq, :left)
    p_end = B.project_point(c, mid + 0.31 * ez + 0.22 * ex)
    B.camera_move!(c, p_end)
    @test sc.placing[].amount == (0.3, 0.2) && sc.place_label.text[] == ["+0.3, +0.2 mm"]
    @test all(==(B._TRIAD_ACTIVE), col(sc.place_rings, square(2)))
    B.camera_release!(c, p_end)
    pl = DesignSession(resolution = :coarse)
    @test nudge_part!(pl.target, i, 0.3; frame = :tcc, axis = :z) && nudge_part!(pl.target, i, 0.2; frame = :tcc, axis = :x)
    Fa, Fb = part_frames(sess.target), part_frames(pl.target)
    @test all(n -> isapprox(Fa[n].origin, Fb[n].origin; atol = 1e-12) && Fa[n].R == Fb[n].R, collect(keys(Fb)))
    @test all(isapprox(u.aim, v.aim; atol = 1e-12) for (u, v) in zip(sess.beams.beams[], pl.beams.beams[]))
    @test record!(sess) && length(sess.history.undo) == n1 + 1 && !record!(sess)
    @test undo!(sess) && _same_parts(sess.target.parts[], parts0)
    @test length(sc.lscene.scene.plots) == nplots
end
