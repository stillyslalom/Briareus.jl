# The design saved in data/session_v13.briar.toml, the session file of the current schema
# that test_workspace.jl loads. It has a hohlraum with a window and closing entrance holes,
# a shield at the top entrance hole in an assembly, a user frame, four quads with two aimed
# in the frame of the bottom entrance hole, and settings changed from their defaults in
# every section. No diagnostic is recorded, since a recorded definition is saved with the
# path of its file. After a change of the schema (and of `SESSION_VERSION`), write the file
# again with `save_session(path, _fixture_design())` and rename it for the new version.
function _fixture_design()
    s = DesignSession(resolution = :coarse)
    ts, bt = s.target, s.beams
    keep = 1:4
    set_beams!(bt, bt.beams[][keep]; preset = :custom, cluster = bt.cluster[][keep])
    @assert add_window!(ts; phi = 60.0, z = 0.5)
    @assert set_leh_closure!(ts, true)
    @assert add_frame!(ts)
    @assert add_part!(ts, :disk; name = :shield, frame = Symbol("hohlraum.leh_top"), z = 0.8, r1 = 2.0, two_sided = true)
    @assert add_part_step!(ts, part_index(ts, :shield), :offset)
    @assert add_assembly!(ts; name = :mount, members = [:shield])
    @assert set_piece_role!(ts, :hohlraum_cap_bottom, :none)
    set_excluded!(ts, :shield, true)
    set_aim_frame!(s, 3:4, Symbol("hohlraum.leh_bottom"))
    edit_beams!(bt, 1, :z, "3.2")
    set_participation!(bt, [2], 1, true)
    set_set_enabled!(bt, [4], false)
    set_cavity!(s.drive, true; films = [:hohlraum_window_1])
    set_beam_model!(s.drive, :visrad; profile = VisRadProfile(; envelope = Dict("Q15T" => 1.2)))
    set_unconverted_ratio!(s.drive, 0.03, 0.01)
    s.drive.spot_radius[] = 0.6
    s.optimize.weights[] = Dict(2 => 1.0, 4 => 0.5)
    s.optimize.iterations[] = 37
    edit_review!(s.review, :peak_fraction, "0.2")
    edit_talis!(s.talis, :tolerance, "0.3")
    set_sweep!(s.sweep; subject = :shield, axis = :x, los = true)
    edit_sweep!(s.sweep, :step, "0.25")
    use_pulse_shape!(s; step = 0.5)
    set_prepulse!(s.time, 0.02, 0.0; start = -4.0, stop = 0.5)
    set_irradiance_threshold!(s.time, 3e8)
    v = s.view
    v.beams[] = :volumes
    set_display!(v, :hohlraum_wall, :wireframe)
    set_cut!(v; kind = :plane, axis = :x, range = (-1.0, Inf), components = [:hohlraum_wall])
    set_cut_plane!(v; normal = (0.0, 1.0, 0.0), offset = 0.4)
    set_triad!(v, :frame_1, true)
    v.hidden_beams[] = Set([bt.beams[][2].name])
    v.camera[] = (; eye = (12.5, -3.25, 6.0), lookat = (0.0, 0.5, 1.0), up = (0.0, 0.0, 1.0), ortho = true,
                  turntable = false)
    set_chart!(s.time, :energy)
    foreach(id -> set_diagnostic!(s, id, false), copy(s.time.diagnostics[]))
    return s
end

# Whether two session states agree: numbers within `tol` (relative or absolute), the rest equal.
_same_state(a::AbstractDict, b::AbstractDict; tol = 1e-12) =
    keys(a) == keys(b) && all(k -> _same_state(a[k], b[k]; tol), keys(a))
_same_state(a::AbstractVector, b::AbstractVector; tol = 1e-12) =
    length(a) == length(b) && all(_same_state(x, y; tol) for (x, y) in zip(a, b))
_same_state(a::Bool, b::Bool; tol = 1e-12) = a == b
_same_state(a::Real, b::Real; tol = 1e-12) = a == b || isapprox(a, b; rtol = tol, atol = tol)
_same_state(a, b; tol = 1e-12) = a == b
