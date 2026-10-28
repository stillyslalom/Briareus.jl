# Frames in the GUI: part frames through `resolve_frames`, beams held in frames, the
# optimizer's frame fields, the placement sweep, and sessions before version 3.

using LinearAlgebra: norm, I
import TOML

const V3 = Briareus.Vec3

# The pre-frames GUI's part rotation, written out: Rz(φ) Ry(θ) Rz(ψ) through radians.
function _old_rotation(θ, φ, ψ)
    θ, φ, ψ = deg2rad.((θ, φ, ψ))
    Rz(a) = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]
    Ry(a) = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)]
    return Rz(φ) * Ry(θ) * Rz(ψ)
end

# Whether every beam keeps `aim ≈ to_chamber(frame, aim_local)`.
_aims_held(bt; atol = 1e-12) =
    all(i -> isapprox(bt.beams[][i].aim, to_chamber(bt.frames[][bt.aim_frame[][i]], bt.aim_local[][i]); atol),
        eachindex(bt.beams[]))

@testset "part frames and the foil frame (no GL)" begin
    # the foil template built from its frame, bit for bit with the core's `foil`
    s = DesignSession(kind = :foil, resolution = :coarse)
    ts = s.target
    p = ts.params[]
    n = max(4, ceil(Int, p[:width] / ts.mesh[].capsule_h))
    ref = foil(FOIL_NAME, foil_center(p), (-1, 0, 0); width = p[:width], n,
               material = Material(:Ti; conversion = p[:conversion]), two_sided = false)
    c = ts.target[].components[findfirst(c -> c.name === FOIL_NAME, ts.target[].components)]
    @test c.mesh.vertices == ref.mesh.vertices && c.mesh.triangles == ref.mesh.triangles
    F = backlighter_frame(p)
    @test F.origin == V3(p[:distance], 0, 0) && F.R[:, 3] == V3(-1, 0, 0)    # front: toward TCC
    @test keys(template_frames(:foil, p)) == Set([:tcc, :backlighter])
    @test keys(template_frames(:hohlraum, p)) == Set([:tcc, :leh_top, :leh_bottom])

    # `euler_rotation` against the old radians formula: the last bits only
    angles = [(θ, φ, ψ) for θ in (0.0, 12.5, 90.0, 133.0, 180.0) for φ in (0.0, 47.0, 270.0)
              for ψ in (0.0, 33.0, 211.0)]
    δ = maximum(maximum(abs, euler_rotation(a...) - _old_rotation(a...)) for a in angles)
    @test δ <= 1e-15

    # the exact Euler angles of a rotation, also with the axis along ±z
    Rs = [euler_rotation(a...) for a in angles]
    @test all(R -> maximum(abs, euler_rotation(Controllers._euler_angles(R)...) - R) < 1e-14, Rs)

    # a grandchild follows its grandparent exactly; vertices against the frames written out
    h = DesignSession(resolution = :coarse)
    th = h.target
    @test add_part!(th, :disk; name = :shield, frame = :leh_top, x = 0.1, z = 0.5, theta = 12,
                    phi = 77, rotation = 31, r1 = 2.0)
    @test add_part!(th, :rectangle; name = :plate, frame = :shield, y = 0.4, z = 0.3, theta = 47.5,
                    phi = 333, rotation = 211, r1 = 0.6, r2 = 0.9)
    L = th.params[][:length]
    pl = th.parts[][part_index(th, :plate)]
    local_mesh = part_component(PartSpec(:rectangle; name = :plate, values = Dict(:r1 => 0.6, :r2 => 0.9)),
                                th.mesh[].h).mesh.vertices
    R1, t1 = _old_rotation(12, 77, 31), [0.1, 0.0, 0.5 + L / 2]
    R2 = R1 * _old_rotation(47.5, 333, 211)
    t2 = R1 * [0.0, 0.4, 0.3] + t1
    got = th.target[].components[findfirst(c -> c.name === :plate, th.target[].components)].mesh.vertices
    @test maximum(norm(g - (R2 * v + t2)) for (g, v) in zip(got, local_mesh)) < 1e-12
    @test part_frames(th)[:plate] ≈ compose(part_frames(th)[:shield], euler_frame((0, 0.4, 0.3), 47.5, 333, 211))
    # a frame that is not there, and a loop through the frames, are refused
    @test !add_part!(th, :sphere; name = :lost, frame = :nowhere) && occursin("no frame named nowhere", th.error[])
    @test !set_part_frame!(th, part_index(th, :shield), :plate) && occursin("loop", th.error[])
    # changing frames keeps a tilted part where it is (the exact angles)
    v0 = copy(got)
    @test set_part_frame!(th, part_index(th, :plate), :leh_bottom)
    got = th.target[].components[findfirst(c -> c.name === :plate, th.target[].components)].mesh.vertices
    @test maximum(norm(a - b) for (a, b) in zip(got, v0)) < 1e-12
end

@testset "beams follow a nudged or turned part (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts, bt = s.target, s.beams
    @test add_part!(ts, :disk; name = :shield, x = 3.0, r1 = 0.5)
    @test add_part!(ts, :sphere; name = :bead, frame = :shield, z = 1.0, r1 = 0.1)
    idx = findall(==(:outer50), bt.cluster[])[1:3]
    set_aim_frame!(s, idx, :shield; keep = false)                 # at the disk's center
    @test all(i -> bt.beams[][i].aim ≈ V3(3, 0, 0), idx) && all(==(:shield), bt.aim_frame[][idx])
    edit_beam_set!(bt, idx, :aim_x, "0.2")
    edit_beam_set!(bt, idx, :aim_z, "0.1")
    aim(i) = bt.beams[][i].aim
    @test all(i -> aim(i) ≈ V3(3.2, 0, 0.1), idx)
    # turned: θ = 90° tilts the disk's z to chamber +x and its x to chamber -z
    @test edit_part!(ts, part_index(ts, :shield), :theta, "90")
    @test all(i -> isapprox(aim(i), V3(3.1, 0, -0.2); atol = 1e-12), idx)
    bead = part_centroid(ts.target[], :bead)
    @test isapprox(bead, V3(4.0, 0, 0); atol = 1e-12)             # the child turned with it
    # moved along its own axis (now chamber +x), then along the chamber's y: one edit each
    record!(s)
    others = [aim(i) for i in eachindex(bt.beams[]) if !(i in idx)]
    @test nudge_part!(ts, part_index(ts, :shield), 0.5; frame = :shield, axis = :z)
    @test all(i -> isapprox(aim(i), V3(3.6, 0, -0.2); atol = 1e-12), idx)
    @test isapprox(part_centroid(ts.target[], :bead), V3(4.5, 0, 0); atol = 1e-12)
    @test [aim(i) for i in eachindex(bt.beams[]) if !(i in idx)] == others     # the rest stay
    record!(s)                                                    # (a window records each edit)
    @test nudge_part!(ts, part_index(ts, :shield), -0.25; frame = :tcc, axis = :y)
    @test all(i -> isapprox(aim(i), V3(3.6, -0.25, -0.2); atol = 1e-12), idx)
    @test !nudge_part!(ts, part_index(ts, :shield), 1.0; frame = :nowhere) && occursin("no frame", ts.error[])
    record!(s)
    @test undo!(s)
    @test ts.parts[][part_index(ts, :shield)].values[:y] == 0 &&
          ts.parts[][part_index(ts, :shield)].values[:x] == 3.5
    @test all(i -> isapprox(aim(i), V3(3.6, 0, -0.2); atol = 1e-12), idx)
    # a part removed: its beams keep their aims, held in the chamber
    a = [aim(i) for i in idx]
    @test remove_part!(ts, part_index(ts, :bead)) && remove_part!(ts, part_index(ts, :shield))
    @test all(==(:tcc), bt.aim_frame[][idx]) && [aim(i) for i in idx] == a
    @test _aims_held(bt)
end

@testset "every aim writer keeps the local aims (no GL)" begin
    s = DesignSession(resolution = :coarse)
    ts, bt = s.target, s.beams
    @test add_part!(ts, :disk; name = :lid, frame = :leh_top, z = 0.5, theta = 20, phi = 30, r1 = 2.0)
    held = Ref(true)
    calls = Ref(0)
    on(_ -> (calls[] += 1; held[] &= _aims_held(bt)), bt.beams)
    record!(s)
    set_aim_frame!(s, 1:8, :lid)
    set_aim_frame!(s, 9:12, :leh_bottom; keep = false)
    edit_beams!(bt, 1, :z, "3.3")
    edit_beam_set!(bt, 1:8, :aim_x, "0.4")
    edit_beam_set!(bt, [1, 9], :aim_y, "0.1")                       # two frames: chamber y
    edit_beam_set!(bt, 1:4, :power, "6.0")
    nudge_part!(ts, part_index(ts, :lid), 0.3; frame = :tcc, axis = :x)
    edit_part!(ts, part_index(ts, :lid), :rotation, "40")
    edit_target!(ts, :length, "10.4")                               # the entrance holes move
    q = bt.beams[][1].quad
    @test split_quad!(bt, q) && merge_quad!(bt, q)
    f = joinpath(mktempdir(), "beams.csv")
    export_beams(bt, f)
    record!(s)
    @test undo!(s) && redo!(s)
    sst = joinpath(mktempdir(), "pointing.csv")
    export_visrad_sst(bt, sst)
    import_visrad_sst!(bt, sst)
    import_beams!(bt, f)
    load_preset!(s, :nif_quads)
    @test held[] && calls[] > 10
    # the optimizer's result is written through the same path
    build_scenario!(s)
    st = s.optimize
    st.fields[] = [:r]
    st.iterations[] = 1
    st.max_evaluations[] = 3
    set_aim_frame!(s, 1:8, :lid)
    start_optimize!(s; spawn = false)
    accept_optimum!(s)
    @test held[] && st.result[].applied

    # and no writer bypasses it: in the controllers, only `_set_aims!` assigns the table's beams
    writes = r"(?:\bbt|\btable|\bbeams|sess\.beams|s\.beams)\.beams(?:\[\]|\.val)\s*=[^=]|_assign!\(\s*(?:bt|[\w.]*\.beams)\.beams\b|notify\(\s*(?:bt|[\w.]*\.beams)\.beams\)"
    src = joinpath(pkgdir(BriareusGUI), "src")
    bad = String[]
    for (root, _, files) in walkdir(src), file in files
        endswith(file, ".jl") || continue
        fn = ""
        for (k, line) in enumerate(eachline(joinpath(root, file)))
            m = match(r"^(?:Base\.@nospecializeinfer )?function ([\w!.]+)", line)
            m === nothing || (fn = m[1])
            occursin(writes, line) && fn != "_set_aims!" && push!(bad, "$file:$k")
        end
    end
    @test isempty(bad)
end

@testset "frame fields of the optimizer (no GL)" begin
    s = DesignSession(kind = :foil, resolution = :coarse)
    bt, st = s.beams, s.optimize
    @test all(==(:backlighter), bt.aim_frame[]) && all(iszero, bt.aim_local[])  # at the foil's center
    @test st.fields[] == [:z, :r]
    set_free!(st, :u, true)
    set_free!(st, :v, true)
    @test st.fields[] == [:u, :v]                                   # the chamber fields fixed
    set_free!(st, :z, true)
    @test st.fields[] == [:z]
    set_free!(st, :z, false)
    set_free!(st, :u, true)
    set_free!(st, :v, true)
    @test edit_optimize!(st, :u, "0.3") && st.deltas[][:u] == 0.3
    set_objective!(st, :irradiance)
    st.by[] = :group
    st.iterations[] = 4
    build_scenario!(s)
    start_optimize!(s; spawn = false)
    r = st.result[]
    @test r !== nothing && r.value <= r.initial_value
    @test all(n -> endswith(n, ".u") || endswith(n, ".v"), r.names) && all(abs.(r.x) .<= 0.3 + 1e-12)
    accept_optimum!(s)
    # each group's common point moved across the foil's face (w = 0: chamber x = distance)
    d = s.target.params[][:distance]
    groups = unique(bt.cluster[])
    loc = bt.aim_local[]
    @test all(g -> (m = findall(==(g), bt.cluster[]); all(i -> isapprox(loc[i], loc[m[1]]; atol = 1e-12), m)), groups)
    @test all(l -> abs(l[3]) < 1e-12, loc) && all(b -> isapprox(b.aim[1], d; atol = 1e-12), bt.beams[])
    @test any(l -> hypot(l[1], l[2]) > 1e-3, loc)
    @test _aims_held(bt)
end

@testset "placement sweep (no GL)" begin
    s = DesignSession(kind = :foil, resolution = :coarse)
    ts, bt, sw = s.target, s.beams, s.sweep
    @test sweep_subject(s) === :backlighter && sweep_frame(s) === :backlighter
    @test sweep_los_object(s) === :main
    θ, φ = sweep_los(s)                                              # foil (+x) → TCC → detector
    @test isapprox(θ, 90; atol = 1e-9) && isapprox(φ, 180; atol = 1e-9)
    set_sweep!(sw; clearance = true, los = true, capture = true)
    @test edit_sweep!(sw, :start, "-1") && edit_sweep!(sw, :stop, "1") && edit_sweep!(sw, :step, "1")
    @test sweep_values(sw) == [-1.0, 0.0, 1.0]
    @test !edit_sweep!(sw, :step, "0") && occursin("positive", sw.error[])
    aims0 = [b.aim for b in bt.beams[]]
    start_sweep!(s; spawn = false)
    r = sw.result[]
    @test r !== nothing && length(r.rows) == 3 && !r.canceled
    @test task_status(s, :sweep)[1] in (:ok, :attention) && !sweep_stale(s)
    # the beams follow the foil: +1 along its front normal is 1 mm toward TCC (-x)
    @test [b.aim for b in r.rows[3].beams] == [a + V3(-1, 0, 0) for a in aims0]
    @test all(w -> w.los_clear === true, r.rows)                    # nothing between foil and detector
    # the captured fraction at zero shift, from the deposition directly
    t = ts.target[]
    op = DepositionOperator(t, active_beams(bt); resolution = SWEEP_QUALITY[:quick].resolution)
    M = beam_component_power(op, active_beams(bt))
    fc = [sum(M[i, part_pieces(t, FOIL_NAME)]) / b.power for (i, b) in enumerate(active_beams(bt))]
    @test r.rows[2].min_captured ≈ minimum(fc) rtol = 1e-12
    ch = sweep_chart(s)
    @test ch.chart === :sweep && length(ch.left.series) == 1 && ch.left.series[1].x == [-1.0, 0.0, 1.0]
    select_sweep_step!(sw, 3)
    @test sweep_chart(s).cursor == 1.0 && length(sweep_rows(sw)) == 3
    # the step applied as one edit, and undone as one
    record!(s)
    before = (copy(ts.params[]), bt.beams[])
    @test apply_sweep_step!(s, 3)
    @test ts.params[][:distance] == before[1][:distance] - 1
    @test all(k -> isapprox(bt.beams[][k].aim, aims0[k] + V3(-1, 0, 0); atol = 1e-12), eachindex(aims0))
    @test sweep_stale(s) && !apply_sweep_step!(s, 1)                # the design moved since
    @test undo!(s)
    @test ts.params[] == before[1] && bt.beams[] == before[2]
    # a part as the subject: it and what is placed in its frame move
    @test add_part!(ts, :disk; name = :shield, frame = :backlighter, z = -0.5, r1 = 0.4)
    set_sweep!(sw; subject = :shield, frame = :tcc, axis = :y, los = false, capture = false)
    start_sweep!(s; spawn = false)
    @test length(sw.result[].rows) == 3 && apply_sweep_step!(s, 1)
    @test isapprox(part_centroid(ts.target[], :shield), V3(6.5, -1, 0); atol = 1e-12)
    # in a window: on a worker, canceled at once
    queue = Channel{Any}(Inf)
    s.deliver[] = f -> put!(queue, f)
    set_sweep!(sw; clearance = false, capture = true)
    edit_sweep!(sw, :stop, "3")
    start_sweep!(s; spawn = true)
    cancel_sweep!(s)
    drain_until(() -> !sw.running[], queue)
    @test startswith(sw.status[], "canceled")
    s.deliver[] = f -> f()
    # settings the design cannot serve start nothing
    set_sweep!(sw; clearance = false, capture = false, los = false, optimize = false)
    start_sweep!(s; spawn = false)
    @test occursin("at least one check", sw.status[])
end

@testset "sessions: version 2 converted, version 3 round trip (no GL)" begin
    file = joinpath(@__DIR__, "data", "session_v2.briar.toml")
    saved = TOML.parsefile(file)
    @test saved["version"] == 2
    s = DesignSession(resolution = :coarse)
    load_session!(s, file)
    ts, bt = s.target, s.beams
    # the saved chamber aims are kept; each beam that followed a component aims in its frame
    @test [b.aim for b in bt.beams[]] == [V3(r["aim"]...) for r in saved["beams"]["beam"]]
    @test bt.aim_frame[] == [:backlighter, :backlighter, :substrate, :tcc, :tcc]
    @test _aims_held(bt)
    # the old foil frame re-expressed in the backlighter frame: the vertices where they were,
    # against the old frames written out (axes [0 0 1; 0 1 0; -1 0 0] at the foil's center)
    @test [p.frame for p in ts.parts[]] == [:backlighter, :substrate, :tcc, :backlighter]
    d = ts.params[][:distance]
    Rfoil, tfoil = [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0], [d, 0.0, 0.0]
    old = Dict(p["name"] => p for p in saved["target"]["parts"])
    vals(n) = old[n]["values"]
    function old_frame(n)
        p = old[n]
        R0, t0 = p["anchor"] == "foil" ? (Rfoil, tfoil) : p["anchor"] == "tcc" ? (Matrix(1.0I, 3, 3), zeros(3)) :
                 old_frame(p["anchor"])
        v = vals(n)
        return R0 * _old_rotation(v["theta"], v["phi"], v["rotation"]), R0 * [v["x"], v["y"], v["z"]] + t0
    end
    worst = 0.0
    for n in keys(old)
        p = ts.parts[][part_index(ts, Symbol(n))]
        sizes = Dict(k => p.values[k] for (k, _, _, _) in PART_SIZE_FIELDS[p.kind])
        loc = part_component(PartSpec(p.kind; name = p.name, values = sizes), ts.mesh[].h).mesh.vertices
        got = ts.target[].components[findfirst(c -> c.name === p.name, ts.target[].components)].mesh.vertices
        R, t = old_frame(n)
        worst = max(worst, maximum(norm(g - (R * v + t)) for (g, v) in zip(got, loc)))
    end
    @test worst < 1e-12
    # saved again: version 3, read back the same
    f = joinpath(mktempdir(), "v3.briar.toml")
    set_sweep!(s.sweep; subject = :substrate, axis = :x, los = true)
    edit_sweep!(s.sweep, :step, "0.25")
    save_session(f, s)
    @test TOML.parsefile(f)["version"] == 3
    s3 = DesignSession(resolution = :coarse)
    load_session!(s3, f)
    @test session_state(s3) == session_state(s)
    @test s3.beams.aim_local[] == bt.aim_local[] && s3.beams.beams[] == bt.beams[]
    @test s3.target.target[].tri_v0 == ts.target[].tri_v0
    @test s3.sweep.subject[] === :substrate && s3.sweep.step[] == 0.25 && s3.sweep.los[]
    @test isnan(s3.sweep.los_theta[])
end
