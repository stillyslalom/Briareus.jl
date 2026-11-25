# Saving and restoring a design. The session is kept as plain data (`session_state`) and
# written to and read from a session file, either a zip of the TOML and the CAD files it
# reads or only the TOML, with the diagnostic definitions the design records. The undo/redo
# history is built on the same snapshots. Results, such as the drive, an optimization or a
# clearance report, are not saved; they are computed again from the restored design.

const SESSION_FORMAT = "BriareusGUI session"
const SESSION_VERSION = 14              # the one schema `restore!` accepts (`session_state` writes it)

_strs(xs) = [String(x) for x in xs]
_syms(xs) = Symbol[Symbol(x) for x in xs]
_fdict(d) = Dict{String,Any}(String(k) => Float64(v) for (k, v) in d)

function _beam_state(b::Beam, enabled::Bool, cluster::Symbol, frame::Symbol, loc, uc)
    sp = b.spot
    return Dict{String,Any}("name" => b.name, "port" => collect(b.port), "aim" => collect(b.aim),
                            "power" => b.power, "defocus" => b.defocus, "fnumber" => b.fnumber,
                            "group" => String(b.group), "quad" => b.quad,
                            "wavelength" => b.wavelength, "dlambda" => b.dlambda,
                            "spot" => Dict{String,Any}("a" => sp.a, "b" => sp.b, "order" => sp.order,
                                                       "angle" => sp.angle, "halo" => sp.halo,
                                                       "halo_scale" => sp.halo_scale,
                                                       "halo_order" => sp.halo_order),
                            "enabled" => enabled, "cluster" => String(cluster),
                            "aim_frame" => String(frame), "aim_local" => collect(loc),
                            "omega1" => uc[1], "omega2" => uc[2])
end

function _state_beam(d)
    s = d["spot"]
    spot = SuperGaussianSpot(Float64(s["a"]), Float64(s["b"]); order = s["order"], angle = s["angle"],
                             halo = s["halo"], halo_scale = s["halo_scale"], halo_order = s["halo_order"])
    return Beam{Float64}(d["name"], Vec3(d["port"]...), Vec3(d["aim"]...), d["power"],
                         d["defocus"], spot, d["fnumber"], Symbol(d["group"]), d["quad"],
                         d["wavelength"], d["dlambda"])
end

"""
    session_state(sess::DesignSession) -> Dict{String,Any}

The design as plain data, ready for TOML: the target (element sizes, parts with their
windows, what is left out, the active foil), the beams (each beam's port, aim, power, focus, spot, cone and
quad, group, whether it is on, and whether its unconverted 1ω and 2ω light takes part),
the frame each beam aims in and its aim there, the settings of the drive (beam model
included), the optimization, the clearance check, the TaLIS checks and the placement
sweep, the diagnostics a drive history records (their ids; a session file carries their
definitions too, see [`save_session`](@ref)), the pulses the beams fire (a pulse table's
data itself, its keys' times and powers, so the session opens without the CSV file, or
the pulse shape and its step; see [`load_pulses!`](@ref)), and the 3D view's coloring, beam
drawing (which beams' 1ω, 2ω and 3ω light is drawn), cutaway, component display modes and
frame axes. It also holds the window as it was left (`window`): the 3D view's camera
pose and projection (`ViewState`'s `camera`), whether the drive shows the 3D view, the
map or the x-ray image (and the image's camera), and the Time tab's chart, so a session
reopens as it was; the edit history leaves these out (they are not undo steps).
Results are not part of it; they are computed again, and neither are the selection,
the model tree's expanded nodes or the beam table's sort. Restore it with
[`restore!`](@ref); the edit history ([`record!`](@ref)) holds the same snapshots.
"""
function session_state(sess::DesignSession)
    ts, bt, dr, op, rv, tl = sess.target, sess.beams, sess.drive, sess.optimize, sess.review, sess.talis
    target = Dict{String,Any}(
        "mesh" => Dict{String,Any}("h" => ts.mesh[].h, "capsule_h" => ts.mesh[].capsule_h),
        "excluded" => sort(_strs(ts.excluded[])),
        "active_foil" => String(something(active_foil(ts), Symbol(""))),
        "parts" => _part_state.(ts.parts[]))
    beams = Dict{String,Any}("preset" => String(bt.preset[]), "by" => String(bt.by[]),
                             "facility" => String(bt.facility[]),
                             "beam" => [_beam_state(b, e, c, r, p, u) for (b, e, c, r, p, u) in
                                        zip(bt.beams[], bt.enabled[], bt.cluster[], bt.aim_frame[],
                                            bt.aim_local[], bt.unconverted[])])
    optimize = Dict{String,Any}(
        "by" => String(op.by[]), "fields" => _strs(op.fields[]), "deltas" => _fdict(op.deltas[]),
        "objective" => String(op.objective[]),
        "weights" => Dict{String,Any}(string(l) => w for (l, w) in op.weights[]),
        "hold_drive" => op.hold_drive[], "keepout" => op.keepout[], "keepout_min" => op.keepout_min[],
        "capsule_keepout" => op.capsule_keepout[], "capsule_min" => op.capsule_min[],
        "spill" => op.spill[], "uniformity" => op.uniformity[], "iterations" => op.iterations[],
        "max_evaluations" => op.max_evaluations[])
    review = Dict{String,Any}(
        "targets" => _strs(rv.targets[]), "hardware" => _strs(rv.hardware[]),
        "peak_fraction" => rv.peak_fraction[], "wings" => rv.wings[], "errors" => _fdict(rv.errors[]),
        "quality" => String(rv.quality[]))
    talis = Dict{String,Any}(
        "capture" => _strs(tl.capture[]), "hardware" => _strs(tl.hardware[]),
        "diagnostics" => _strs(tl.diagnostics[]),
        "tolerance" => tl.tolerance[], "pivot" => String(tl.pivot[]),
        "pivot_point" => collect(tl.pivot_point[]), "omega1" => tl.omega1[])
    state = Dict{String,Any}("format" => SESSION_FORMAT, "version" => SESSION_VERSION,
                             "visrad_source" => sess.visrad_source[],
                             "target" => target, "beams" => beams,
                             "drive" => Dict{String,Any}("spot_radius" => dr.spot_radius[],
                                                         "cavity" => dr.cavity[],
                                                         "cavity_point" => collect(dr.cavity_point[]),
                                                         "films" => sort(_strs(dr.films[])),
                                                         "beam_model" => String(dr.beam_model[]),
                                                         "profile" => _profile_state(dr.profile[]),
                                                         "unconverted_ratio" => collect(dr.unconverted_ratio[])),
                             "optimize" => optimize, "review" => review, "talis" => talis,
                             "sweep" => _sweep_state(sess.sweep),
                             "time" => Dict{String,Any}("diagnostics" => copy(sess.time.diagnostics[]),
                                                        "prepulse" => collect(sess.time.prepulse[]),
                                                        "threshold" => sess.time.threshold[]))
    pulses = _pulse_state(sess.time.pulses[])
    pulses === nothing || (state["time"]["pulses"] = pulses)
    # the parts the diagnostics' views name in this session (after renames), when any
    vp = sess.time.view_parts[]
    isempty(vp) || (state["time"]["view_parts"] = Dict{String,Any}(id => Dict{String,Any}(d) for (id, d) in vp))
    state["window"] = _window_state(sess)
    v = sess.view
    state["view"] = Dict{String,Any}("color" => String(v.color[]), "beams" => String(v.beams[]),
                                     "beam_color" => String(v.beam_color[]),
                                     "beam_length" => v.beam_length[], "cutaway" => v.cutaway[],
                                     "cut_azimuth" => v.cut_azimuth[], "cut_width" => v.cut_width[],
                                     "cut_components" => sort(_strs(v.cut_components[])),
                                     "cut_kind" => String(v.cut_kind[]),
                                     "cut_axis" => String(v.cut_axis[]),
                                     "cut_range" => collect(v.cut_range[]),
                                     "cut_point" => collect(v.cut_point[]),
                                     "cut_normal" => collect(v.cut_normal[]),
                                     "cut_shift" => v.cut_shift[],
                                     "cut_line" => collect(v.cut_line[]),
                                     "cut_rotation" => v.cut_rotation[],
                                     "display" => Dict{String,Any}(String(k) => String(m)
                                                                   for (k, m) in v.display[]),
                                     "hidden_beams" => sort(collect(v.hidden_beams[])),
                                     "shown_1w" => sort(collect(v.shown_1w[])),
                                     "shown_2w" => sort(collect(v.shown_2w[])),
                                     "side" => dr.side[], "show_off" => v.show_off[],
                                     "triads" => sort(_strs(v.triads[])))
    return state
end

"""
    restore!(sess::DesignSession, state) -> sess

Set the session to a [`session_state`](@ref), view included. The target is rebuilt once; analyses whose inputs changed become stale
as after any edit, and an optimization result stays, current when the restored design
is the one it ran on ([`optimize_stale`](@ref)). Throws `ArgumentError` for a state of another format, or one the
geometry rejects (the session is then left as it was before the target, if the
target fails, or partly restored otherwise).

The state must have the current schema (`version` equal to `SESSION_VERSION`): a state of
another version is refused with an `ArgumentError` that names the version found and the
version supported, and nothing is converted. A state with no `window` section (an undo
step) leaves the window as it is.
"""
restore!(sess::DesignSession, state::AbstractDict) = _restore!(sess, state, false)

# `keep_target`: the state's target section is the design's (an undo step's): the target
# is kept unless a build, remeshing, capsule fit or rejected edit is in flight.
function _restore!(sess::DesignSession, state::AbstractDict, keep_target::Bool)
    get(state, "format", "") == SESSION_FORMAT ||
        throw(ArgumentError("not a $SESSION_FORMAT (format \"$(get(state, "format", ""))\")"))
    version = get(state, "version", nothing)
    version == SESSION_VERSION ||
        throw(ArgumentError("session version $(something(version, "none")) is not supported: " *
                            "this program loads version $SESSION_VERSION only"))
    ts, bt = sess.target, sess.beams
    t = state["target"]
    mesh = MeshSpec(t["mesh"]["h"], t["mesh"]["capsule_h"])
    parts = PartSpec[_state_part(p) for p in t["parts"]]
    # (an undo step whose target is the design's, built and settled, keeps it as it is)
    if !(keep_target && !ts.building[] && isempty(ts.pending[]) && isempty(ts.fitting[]) && isempty(ts.error[]))
        _cancel_fits!(ts)                 # a capsule fit in flight was for another design
        _try_build!(ts; mesh, parts, excluded = Set(_syms(t["excluded"]))) ||
            throw(ArgumentError("the saved target was rejected: $(ts.error[])"))
    end
    a = t["active_foil"]
    _assign!(ts.active_foil, isempty(a) ? nothing : Symbol(a))
    _assign!(sess.visrad_source, state["visrad_source"])
    b = state["beams"]
    rows = b["beam"]
    # Aims are held in their frames. A beam whose frame the target lacks aims in `:tcc`,
    # and the saved chamber aims are then kept and re-expressed in their frames.
    frames = part_frames(ts)
    function frame_of(r)
        f = Symbol(r["aim_frame"])
        return f === :none || !haskey(frames, f) ? :tcc : f
    end
    af = Symbol[frame_of(r) for r in rows]
    local_saved = all(r -> haskey(frames, Symbol(r["aim_frame"])), rows)
    set_beams!(bt, [_state_beam(r) for r in rows]; preset = Symbol(b["preset"]),
               enabled = Bool[r["enabled"] for r in rows], cluster = Symbol[Symbol(r["cluster"]) for r in rows],
               aim_frame = af,
               aim_local = local_saved ? [Vec3(Float64.(r["aim_local"])...) for r in rows] : nothing,
               unconverted = [(Bool(r["omega1"]), Bool(r["omega2"])) for r in rows])
    _assign!(bt.by, Symbol(b["by"]))
    _assign!(bt.facility, Symbol(b["facility"]))
    d = state["drive"]
    _assign!(sess.drive.spot_radius, Float64(d["spot_radius"]))
    set_cavity!(sess.drive, d["cavity"]; point = Float64.(d["cavity_point"]), films = _syms(d["films"]))
    set_beam_model!(sess.drive, Symbol(d["beam_model"]); profile = _state_profile(d["profile"]))
    set_unconverted_ratio!(sess.drive, Float64.(d["unconverted_ratio"])...)
    # An optimization result stays, since it is current while the design is the one it
    # ran on.
    o, op = state["optimize"], sess.optimize
    _assign!(op.by, Symbol(o["by"]))
    _assign!(op.fields, _syms(o["fields"]))
    _assign!(op.deltas, Dict{Symbol,Float64}(Symbol(k) => v for (k, v) in o["deltas"]))
    _assign!(op.objective, Symbol(o["objective"]))
    _assign!(op.weights, Dict{Int,Float64}(parse(Int, k) => v for (k, v) in o["weights"]))
    for k in (:hold_drive, :keepout, :keepout_min, :capsule_keepout, :capsule_min, :spill,
              :uniformity, :iterations, :max_evaluations)
        obs = getfield(op, k)
        _assign!(obs, convert(typeof(obs[]), o[String(k)]))
    end
    # the checks' component roles once the saved design's target is in place (a target
    # built on a worker would otherwise give its new components their default roles)
    r, c = state["review"], state["talis"]
    roles = function ()
        _assign!(sess.review.targets, _syms(r["targets"]))
        _assign!(sess.review.hardware, _syms(r["hardware"]))
        _assign!(sess.talis.capture, _syms(c["capture"]))
        _assign!(sess.talis.hardware, _syms(c["hardware"]))
        _assign!(sess.talis.diagnostics, _syms(c["diagnostics"]))
    end
    roles()
    _after_build!(roles, ts)
    rv, tl = sess.review, sess.talis
    _assign!(rv.peak_fraction, Float64(r["peak_fraction"]))
    _assign!(rv.wings, r["wings"])
    _assign!(rv.errors, Dict{Symbol,Float64}(Symbol(k) => v for (k, v) in r["errors"]))
    _assign!(rv.quality, Symbol(r["quality"]))
    _assign!(tl.tolerance, Float64(c["tolerance"]))
    _assign!(tl.pivot, Symbol(c["pivot"]))
    _assign!(tl.pivot_point, Tuple(Float64.(c["pivot_point"])))
    _assign!(tl.omega1, c["omega1"])
    _restore_sweep!(sess.sweep, state["sweep"])
    tms = state["time"]
    ids = String[String(x) for x in tms["diagnostics"]]
    ids == sess.time.diagnostics[] || (sess.time.diagnostics[] = ids)
    # (the parts the diagnostics' views name are saved when there are any)
    vp = Dict{String,Dict{String,String}}(String(id) => Dict{String,String}(String(k) => String(v) for (k, v) in d)
                                          for (id, d) in get(tms, "view_parts", Dict()))
    _assign!(sess.time.view_parts, vp)
    p1, p2, t0, t1 = Float64.(tms["prepulse"])
    set_prepulse!(sess.time, p1, p2; start = t0, stop = t1)
    set_irradiance_threshold!(sess.time, Float64(tms["threshold"]))
    # The pulses are saved when the session has any. Pulses equal to the loaded ones stay
    # as they are.
    p = get(tms, "pulses", nothing)
    if p != _pulse_state(sess.time.pulses[])
        p === nothing ? (sess.time.pulses[] = nothing) : _set_pulses!(sess, _state_pulses(p))
    end
    _restore_view!(sess, state["view"])
    haskey(state, "window") && _restore_window!(sess, state["window"])
    _restored_optimum!(sess)
    return sess
end

# The view section of a session state (`session_state`): the coloring, the beam drawing,
# the cut, the display modes, the frame axes.
function _restore_view!(sess::DesignSession, v::AbstractDict)
    vs = sess.view
    _assign!(vs.color, Symbol(v["color"]))
    _assign!(vs.beams, Symbol(v["beams"]))
    _assign!(vs.beam_color, Symbol(v["beam_color"]))
    _assign!(vs.beam_length, Float64(v["beam_length"]))
    _assign!(vs.cutaway, v["cutaway"])
    set_cut!(vs; azimuth = v["cut_azimuth"], width = v["cut_width"], components = _syms(v["cut_components"]),
             kind = Symbol(v["cut_kind"]), axis = Symbol(v["cut_axis"]), range = Tuple(Float64.(v["cut_range"])))
    # The plane is set without making the cut a plane cut on the way, since a sector cut
    # that briefly became a plane cut redrew the 3D view twice on every undo.
    _set_plane!(vs, Vec3(Float64.(v["cut_normal"])...), Float64.(v["cut_point"]), Float64(v["cut_shift"]))
    _assign!(vs.cut_line, Tuple(Float64.(v["cut_line"])))
    _assign!(vs.cut_rotation, Float64(v["cut_rotation"]))
    _assign!(vs.cut_kind, Symbol(v["cut_kind"]))
    disp = Dict{Symbol,Symbol}(Symbol(k) => Symbol(m) for (k, m) in v["display"])
    filter!(p -> last(p) in DISPLAY_MODES && last(p) !== :opaque, disp)
    _assign!(vs.display, disp)
    _assign!(vs.hidden_beams, Set{String}(v["hidden_beams"]))
    _assign!(vs.shown_1w, Set{String}(v["shown_1w"]))
    _assign!(vs.shown_2w, Set{String}(v["shown_2w"]))
    _assign!(sess.drive.side, Bool(v["side"]))
    _assign!(vs.show_off, Bool(v["show_off"]))
    _assign!(vs.triads, Set{Symbol}(_syms(v["triads"])))
    return sess
end

# The window as it was left (`session_state`'s `window`): the camera pose (none for the
# default view), the view of the drive (`:scene`, `:map`, `:image`) with the camera of
# the x-ray image, and the Time chart.
function _window_state(sess::DesignSession)
    w = Dict{String,Any}("drive_view" => String(sess.drive.view[]), "time_chart" => String(sess.time.chart[]))
    c = sess.view.camera[]
    c === nothing ||
        (w["camera"] = Dict{String,Any}("eye" => collect(Float64, c.eye), "lookat" => collect(Float64, c.lookat),
                                        "up" => collect(Float64, c.up), "ortho" => c.ortho,
                                        "turntable" => c.turntable))
    r = sess.drive.image_request[]
    r === nothing ||
        (w["image"] = Dict{String,Any}("direction" => collect(Float64, r.direction), "up" => collect(Float64, r.up),
                                       "center" => collect(Float64, r.center), "width" => r.width, "nx" => r.nx))
    return w
end

# Set the window from `_window_state`'s data; what is missing takes its default (the
# isometric camera, the 3D view, the Pulses chart). An image view without its camera (a
# diagnostic's image, which is not saved) opens on the 3D view; the image itself is
# rendered when the drive is computed.
function _restore_window!(sess::DesignSession, w::AbstractDict)
    t3(x) = Tuple(Float64.(x))
    c = get(w, "camera", nothing)
    pose = c === nothing ? nothing :
           (; eye = t3(c["eye"]), lookat = t3(c["lookat"]), up = t3(c["up"]), ortho = Bool(c["ortho"]),
            turntable = Bool(c["turntable"]))
    isequal(pose, sess.view.camera[]) || (sess.view.camera[] = pose)
    dr = sess.drive
    im = get(w, "image", nothing)
    req = im === nothing ? nothing :
          (; direction = Vec3(t3(im["direction"])...), up = Vec3(t3(im["up"])...),
           center = Vec3(t3(im["center"])...), width = Float64(im["width"]), nx = Int(im["nx"]))
    isequal(req, dr.image_request[]) || (dr.image_request[] = req; dr.image[] = nothing)
    view = Symbol(get(w, "drive_view", "scene"))
    (view in (:scene, :map, :image) && !(view === :image && req === nothing)) || (view = :scene)
    set_view!(dr, view)
    chart = Symbol(get(w, "time_chart", "pulses"))
    set_chart!(sess.time, chart in TIME_CHARTS ? chart : :pulses)
    return sess
end

# A part as saved: what every node of the model tree has and, for a part with geometry,
# the table `geometry` (`PartGeometry`). A frame or an assembly has no such table.
function _part_state(p::PartSpec)
    d = Dict{String,Any}("kind" => String(p.kind), "name" => String(p.name), "values" => _fdict(p.values),
                         "frame" => String(p.frame), "coords" => String(p.coords), "object" => p.object,
                         "steps" => _step_state.(p.steps), "assembly" => String(p.assembly))
    g = p.geometry
    g === nothing ||
        (d["geometry"] = Dict{String,Any}(
            "two_sided" => g.two_sided, "inward" => g.inward, "aperture" => g.aperture, "file" => g.file,
            "foil" => g.foil, "observed" => g.observed, "wall" => g.wall,
            "windows" => _window_state.(g.windows), "sections" => [[z, r] for (z, r) in g.sections],
            "piece_roles" => Dict{String,Any}(String(k) => String(r) for (k, r) in g.piece_roles)))
    return d
end

# A saved part.
function _state_part(p)
    kind = Symbol(p["kind"])
    bare = _no_geometry(kind)
    bare == !haskey(p, "geometry") ||
        throw(ArgumentError("part $(p["name"]): a $kind " * (bare ? "has no geometry" : "needs its geometry")))
    g = get(p, "geometry", nothing)
    geometry = g === nothing ? (;) :
               (; two_sided = g["two_sided"], inward = g["inward"], aperture = g["aperture"], file = g["file"],
                foil = g["foil"], observed = g["observed"], wall = g["wall"],
                windows = WindowSpec[_state_window(w) for w in g["windows"]],
                sections = NTuple{2,Float64}[(s[1], s[2]) for s in g["sections"]],
                piece_roles = Dict{Symbol,Symbol}(Symbol(k) => Symbol(r) for (k, r) in g["piece_roles"]))
    return PartSpec(kind; name = Symbol(p["name"]),
                    values = Dict{Symbol,Float64}(Symbol(k) => v for (k, v) in p["values"]),
                    frame = Symbol(p["frame"]), coords = Symbol(p["coords"]), object = p["object"],
                    steps = _state_step.(p["steps"]), assembly = Symbol(p["assembly"]), geometry...)
end

_window_state(w::WindowSpec) =
    Dict{String,Any}("phi" => w.phi, "z" => w.z, "width" => w.width, "height" => w.height,
                     "open" => w.open, "albedo" => w.albedo, "conversion" => w.conversion, "round" => w.round)

_state_window(w) = WindowSpec(; phi = w["phi"], z = w["z"], width = w["width"], height = w["height"],
                              open = w["open"], albedo = w["albedo"], conversion = w["conversion"],
                              round = w["round"])

# ---------------------------------------------------------------- CAD files in sessions

_sha256(file) = bytes2hex(open(sha256, file))

# The folder saved CAD copies are extracted to: one subfolder per content hash, so a
# copy's path (and the remesh caches keyed on it) is the same in every session that has it.
_cad_cache(hash) = joinpath(tempdir(), "BriareusGUI-cad", hash)

"""
    cad_source(sess::DesignSession, file) -> String

The file a CAD part's `file` was saved from: for a copy extracted from a `.briar`
session, the path it had when the session was saved; otherwise `file` itself.
"""
cad_source(sess::DesignSession, file::AbstractString) = get(sess.cad_sources, file, String(file))

# The CAD files the design reads, from its STL and VisRad parts and the VisRad export
# template.
function _cad_files(sess::DesignSession)
    files = String[p.file for p in sess.target.parts[] if p.kind in (:mesh, :visrad) && !isempty(p.file)]
    isempty(sess.visrad_source[]) || push!(files, sess.visrad_source[])
    return unique!(files)
end

# The state as written to a file: each CAD file by the path it was saved from and its
# content hash, and the definition of each recorded diagnostic (`_diagnostic_state`).
# Returns the state, the hash => file of each CAD file found and of each CSV file of a
# definition (stored in an archive when `embed`).
function _portable_state(sess::DesignSession; embed::Bool = true)
    state = session_state(sess)
    hashes = Dict{String,String}()
    for f in _cad_files(sess)
        isfile(f) && (hashes[f] = _sha256(f))
    end
    blobs = Dict{String,String}(h => f for (f, h) in hashes)
    for p in state["target"]["parts"]
        g = get(p, "geometry", nothing)
        (g === nothing || isempty(g["file"])) && continue
        f = g["file"]
        haskey(hashes, f) && (g["sha256"] = hashes[f])
        g["file"] = cad_source(sess, f)
    end
    v = sess.visrad_source[]
    if !isempty(v)
        haskey(hashes, v) && (state["visrad_sha256"] = hashes[v])
        state["visrad_source"] = cad_source(sess, v)
    end
    csvs = Dict{String,String}()
    defs = _diagnostic_state(sess, csvs; embed)
    isempty(defs) || (state["diagnostic_definitions"] = defs)
    return state, blobs, csvs
end

# ---------------------------------------------------------------- diagnostic definitions in sessions

# The definition of each recorded diagnostic the session can build, by id: `file` (the
# file that defines it), `chain` (the files it was resolved from), `origin` (builtin,
# user, session) and `definition` (`diagnostic_table`). Each CSV file is named
# `diagnostics/<sha256>.csv` and listed in `csvs` (hash => file) when `embed`, else
# referenced by its path.
function _diagnostic_state(sess::DesignSession, csvs; embed::Bool)
    rows = Dict(r.id => r for r in diagnostic_catalog(sess))
    out = Dict{String,Any}()
    for id in sess.time.diagnostics[]
        r = get(rows, id, nothing)
        (r === nothing || !r.usable) && continue
        def = Controllers._definition(sess, id)
        def === nothing && continue
        csv = function (f)
            (embed && isfile(f)) || return f
            h = _sha256(f)
            csvs[h] = f
            return _csv_entry(h)
        end
        out[id] = Dict{String,Any}("file" => r.file, "chain" => copy(r.chain), "origin" => String(r.origin),
                                   "definition" => diagnostic_table(def; csv))
    end
    return out
end

_csv_entry(hash) = "diagnostics/" * hash * ".csv"
_csv_cache() = joinpath(tempdir(), "BriareusGUI-diagnostics")

# The definitions saved in a session file, by id: `(; def, file, chain)`; an embedded
# CSV is extracted to a cache folder. A definition that no longer builds is left out
# and named in `bad`.
function _read_diagnostic_state(state, zr, bad::Vector{String})
    out = Dict{String,Any}()
    for (id, e) in get(state, "diagnostic_definitions", Dict())
        table = Briareus.GUIAPI.map_csv!(deepcopy(e["definition"]), function (v)
            (zr === nothing || !startswith(v, "diagnostics/")) && return v
            zip_findlast_entry(zr, v) === nothing && return v
            f = joinpath(_csv_cache(), basename(v))
            if !isfile(f)
                mkpath(dirname(f))
                write(f * ".part", zip_readentry(zr, v))
                mv(f * ".part", f; force = true)
            end
            return f
        end)
        def = try
            diagnostic_definition(table; file = e["file"])
        catch err
            err isa ArgumentError || rethrow()
            push!(bad, id)
            continue
        end
        out[id] = (; def, file = String(e["file"]), chain = String.(e["chain"]))
    end
    return out
end

# The zip entry of a CAD file, named by its hash and the extension of its source.
_cad_entry(hash, source) = "cad/" * hash * lowercase(splitext(source)[2])

function _write_toml(io::IO, state)
    println(io, "# ", SESSION_FORMAT, ": a Briareus design (target, beams, analysis settings).")
    println(io, "# Lengths in mm, power in TW, angles in degrees; results are not saved.")
    TOML.print(io, state; sorted = true)
end

"""
    save_session(path, sess::DesignSession; embed = !endswith(path, ".toml")) -> path

Write the design ([`session_state`](@ref), view included) to a session file.

- `embed = true` (a `.briar` file): a zip archive of the design as TOML
  (`session.toml`) and a copy of every CAD file it reads (`cad/<sha256>.<ext>`: STL
  parts, VisRad workspaces), so the file opens on another computer. A file used by
  several parts is stored once. The definition of each diagnostic the design records
  is saved in the TOML (`diagnostic_table`) with its CSV files
  (`diagnostics/<sha256>.csv`), so the session records the same instruments
  anywhere; the status names those from outside the built-ins and those not public,
  which the file then carries.
- `embed = false` (a `.briar.toml` file): only the TOML, plain text that can be read
  and edited by hand. CAD parts are referenced by path, so the file travels without
  its geometry (for CAD that may not leave the computer it is on); diagnostic
  definitions are saved in it, their CSV files referenced by path.

Either way each CAD part records the path it was read from and its content's SHA-256,
which [`load_session!`](@ref) checks.
"""
function save_session(path::AbstractString, sess::DesignSession; embed::Bool = !endswith(path, ".toml"))
    state, blobs, csvs = _portable_state(sess; embed)
    if embed
        tmp = path * ".part"
        ZipWriter(tmp) do w
            zip_newfile(w, "session.toml"; compress = true)
            _write_toml(w, state)
            zip_commitfile(w)
            for (h, f) in sort!(collect(blobs))
                zip_newfile(w, _cad_entry(h, cad_source(sess, f)); compress = true)
                write(w, read(f))
                zip_commitfile(w)
            end
            for (h, f) in sort!(collect(csvs))
                zip_newfile(w, _csv_entry(h); compress = true)
                write(w, read(f))
                zip_commitfile(w)
            end
        end
        mv(tmp, path; force = true)            # an interrupted save leaves the old file
        mb = sum(filesize, values(blobs); init = 0) / 2^20
        sess.status[] = "Saved the session to $(basename(path))" *
                        (isempty(blobs) ? "" : @sprintf(" with %d CAD file%s (%.1f MB)", length(blobs),
                                                        length(blobs) == 1 ? "" : "s", mb)) *
                        _diagnostics_note(state)
    else
        open(io -> _write_toml(io, state), path, "w")
        sess.status[] = "Saved the session to $(basename(path))" *
                        (isempty(blobs) ? "" : " (CAD files by path, not embedded)") *
                        _diagnostics_note(state)
    end
    return path
end

# What a saved session carries of diagnostic definitions beyond the built-ins.
function _diagnostics_note(state)
    defs = get(state, "diagnostic_definitions", Dict())
    other = sort([id for (id, e) in defs if e["origin"] != "builtin"])
    private = sort([id for (id, e) in defs if !get(e["definition"], "public", true)])
    return (isempty(other) ? "" : "; with diagnostic definitions from outside the built-ins: " * join(other, ", ")) *
           (isempty(private) ? "" : "; not public: " * join(private, ", "))
end

# Point a saved CAD reference at a file the session reads: the embedded copy when the
# archive has one (extracted to the cache), else the path itself. Returns that file and
# whether its source on disk differs from what was saved.
function _resolve_cad!(sess::DesignSession, zr, source::AbstractString, hash)
    hash === nothing && return String(source), false
    changed = isfile(source) && _sha256(source) != hash
    entry = zr === nothing ? nothing : _cad_entry(hash, source)
    if entry === nothing || zip_findlast_entry(zr, entry) === nothing
        return String(source), changed      # by reference: the design reads the source
    end
    dir = _cad_cache(hash)
    file = joinpath(dir, basename(source))
    if !(isfile(file) && _sha256(file) == hash)
        mkpath(dir)
        tmp = file * ".part"
        write(tmp, zip_readentry(zr, entry))
        mv(tmp, file; force = true)
    end
    sess.cad_sources[file] = source
    return file, changed
end

"""
    load_session!(sess::DesignSession, path) -> sess

Read a session file written by [`save_session`](@ref) (a `.briar` archive or a
`.briar.toml` file) into the session ([`restore!`](@ref)). The
edit history starts over at the loaded design.

A CAD file embedded in the archive is the one the design reads (extracted to a cache
folder; [`cad_source`](@ref) gives the path it was saved from). When that source still
exists but its content changed since the save, the copy is kept and the source is listed
in `sess.cad_changed`; [`reload_cad_sources!`](@ref) switches to the sources and
[`relink_cad!`](@ref) to another file. A session saved without its CAD reads each
file from its path; the status names those whose content changed since the save.

The diagnostic definitions saved with the session are the ones its drive histories
use; the status names those that differ from the search path's definition of the same
id (the keys are in `diagnostic_catalog`), and [`use_path_diagnostic!`](@ref) switches
to the search path's.
"""
function load_session!(sess::DesignSession, path::AbstractString)
    bytes = read(path)
    zr = length(bytes) >= 4 && bytes[1:4] == b"PK\x03\x04" ? ZipReader(bytes) : nothing
    state = TOML.parse(zr === nothing ? String(bytes) : zip_readentry(zr, "session.toml", String))
    changed = String[]
    if haskey(state, "target") && haskey(state["target"], "parts")
        for p in state["target"]["parts"]
            g = get(p, "geometry", nothing)
            (g === nothing || isempty(g["file"])) && continue
            f, c = _resolve_cad!(sess, zr, g["file"], get(g, "sha256", nothing))
            g["file"] = f
            c && push!(changed, f)
        end
    end
    if !isempty(get(state, "visrad_source", ""))
        f, c = _resolve_cad!(sess, zr, state["visrad_source"], get(state, "visrad_sha256", nothing))
        state["visrad_source"] = f
        c && push!(changed, f)
    end
    bad = String[]
    embedded = _read_diagnostic_state(state, zr, bad)
    clear_optimize!(sess.optimize)                        # (a result of the design before)
    restore!(sess, state)
    sess.time.embedded[] = embedded
    rows = refresh_diagnostics!(sess)
    differ = [r.id for r in rows if !isempty(r.changed)]
    reset_history!(sess)
    unique!(changed)
    copies = filter(f -> haskey(sess.cad_sources, f), changed)
    ondisk = setdiff(changed, copies)
    listed(fs) = join(basename.(cad_source.(Ref(sess), fs)), ", ")
    sess.cad_changed[] = copies
    sess.status[] = "Loaded the session from $(basename(path))" *
                    (isempty(copies) ? "" : "; CAD changed on disk since the save (the design reads " *
                                             "the saved copies): " * listed(copies)) *
                    (isempty(ondisk) ? "" : "; CAD read from disk, changed since the save: " * listed(ondisk)) *
                    (isempty(differ) ? "" : "; diagnostics saved with the session differ from the search path's " *
                                            "(the saved ones are used): " * join(differ, ", ")) *
                    (isempty(bad) ? "" : "; saved diagnostic definitions that no longer build: " * join(bad, ", "))
    return sess
end

"""
    relink_cad!(sess::DesignSession, file, newfile) -> Bool

Read every CAD part that reads `file` (and the VisRad export template, if it is `file`)
from `newfile` instead, and rebuild the target. Returns whether the target accepted it
(an error is in `sess.target.error` otherwise, and the design is unchanged).
"""
function relink_cad!(sess::DesignSession, file::AbstractString, newfile::AbstractString)
    ts = sess.target
    newfile = String(newfile)
    isfile(newfile) || (ts.error[] = "file not found: $newfile"; return false)
    ps = [p.file == file ? _with(p; file = newfile) : p for p in ts.parts[]]
    ps == ts.parts[] || _try_build!(ts; parts = ps) || return false
    sess.visrad_source[] == file && (sess.visrad_source[] = newfile)
    sess.cad_changed[] = filter(!=(file), sess.cad_changed[])
    return true
end

"""
    reload_cad_sources!(sess::DesignSession, files = sess.cad_changed[]) -> Bool

Read each of `files` (saved CAD copies) from the path it was saved from
([`cad_source`](@ref)), for a source that changed since the save. Returns whether every
one was accepted ([`relink_cad!`](@ref)).
"""
function reload_cad_sources!(sess::DesignSession, files = copy(sess.cad_changed[]))
    ok = all([relink_cad!(sess, f, cad_source(sess, f)) for f in files])
    ok && !isempty(files) && (sess.status[] = "Read the CAD from its source: " *
                                              join(unique(basename.(cad_source.(Ref(sess), files))), ", "))
    return ok
end

# ---------------------------------------------------------------- undo and redo

"""
    record!(sess::DesignSession) -> Bool

Note the design's current state in the edit history: if it differs from the last one
noted, the last one becomes an undo step (and the redo steps are dropped). Returns
whether a step was added. A window calls this after every event, so any edit (a field,
a preset, an applied optimization, loaded pulses) can be undone.

The history holds what a session saves ([`session_state`](@ref)), view included: a
change of display mode, cutaway, beam drawing, hidden beams, drawn colors or frame axes
is a step of the one history, so undo takes back the last thing done, whichever it was,
and view state keyed to part names follows the undo of a rename. The camera, the
selection and the model tree's expanded nodes are not saved and not steps, and the
window as it was left (the camera, the drive's view, the Time chart) is saved but not a
step: an undo leaves it as it is. An undo
that changes only the view leaves the target as it is (no rebuild). Consecutive changes
of the same numeric view settings and nothing else (a cutaway slider dragged) merge
into one step, so the step returns to where the drag started. The diagnostics a drive
history records (ticked in the list) are saved but not steps: an undo keeps them as
they are.
"""
function record!(sess::DesignSession)
    h = sess.history
    h.current !== nothing && h.changes[] == h.seen && return false      # nothing changed since
    state = _snapshot!(sess)
    if h.current === nothing
        h.current = state
        return false
    end
    isequal(state, h.current) && return false
    if !isempty(h.undo) && isempty(h.redo)
        now, last = _numeric_view_change(h.current, state), _numeric_view_change(h.undo[end], h.current)
        if now !== nothing && now == last
            h.current = state                 # the same drag, so the step keeps its start
            return false
        end
    end
    push!(h.undo, h.current)
    length(h.undo) > h.limit && popfirst!(h.undo)
    h.current = state
    empty!(h.redo)
    return true
end

"""
    undo!(sess::DesignSession) -> Bool

Return the design to its state before the last edit (`false` if there is none).
"""
function undo!(sess::DesignSession)
    record!(sess)
    h = sess.history
    isempty(h.undo) && return false
    push!(h.redo, h.current)
    _step_to!(sess, pop!(h.undo))
    sess.status[] = "Undid an edit"
    return true
end

"""
    redo!(sess::DesignSession) -> Bool

Repeat the last undone edit (`false` if there is none).
"""
function redo!(sess::DesignSession)
    record!(sess)
    h = sess.history
    isempty(h.redo) && return false
    push!(h.undo, h.current)
    _step_to!(sess, pop!(h.redo))
    sess.status[] = "Redid an edit"
    return true
end

# A snapshot of the edit history: the session state without the diagnostics recorded
# and the window (camera, views and charts), which are not undo steps.
function _history_state(sess::DesignSession)
    state = session_state(sess)
    delete!(state["time"], "diagnostics")
    delete!(state, "window")
    return state
end

# A snapshot of the design as it is, which `current` matches once it is noted: the
# count of changes is read first, so a change made while the snapshot is taken counts
# as one after it.
function _snapshot!(sess::DesignSession)
    n = sess.history.changes[]
    state = _history_state(sess)
    sess.history.seen = n
    return state
end

# The view keys whose values changed from state `a` to `b` when only numeric view
# settings changed (a slider's), else `nothing`.
function _numeric_view_change(a, b)
    all(k -> k == "view" || isequal(get(a, k, nothing), get(b, k, nothing)), union(keys(a), keys(b))) ||
        return nothing
    va, vb = a["view"], b["view"]
    ks = Set(k for k in union(keys(va), keys(vb)) if !isequal(get(va, k, nothing), get(vb, k, nothing)))
    (isempty(ks) || !all(k -> get(va, k, nothing) isa Real && get(vb, k, nothing) isa Real &&
                              !(get(va, k, nothing) isa Bool), ks)) && return nothing
    return ks
end

# A step of the history: only the view when nothing else differs from the current state
# (the design and its target stay as they are), else the whole state, the target kept
# when its section is the current one (a beam or setting edit undone). `current` is the
# design as it is: undo and redo record it first.
function _step_to!(sess::DesignSession, state)
    cur = sess.history.current
    if cur !== nothing && all(k -> k == "view" || isequal(get(state, k, nothing), get(cur, k, nothing)),
                              union(keys(state), keys(cur)))
        _restore_view!(sess, state["view"])
    else
        # The ticked diagnostics stay as they are, since they are not steps.
        tm = merge(get(state, "time", Dict{String,Any}()),
                   Dict{String,Any}("diagnostics" => copy(sess.time.diagnostics[])))
        _restore!(sess, merge(state, Dict{String,Any}("time" => tm)),
                  cur !== nothing && isequal(get(state, "target", nothing), get(cur, "target", nothing)))
    end
    sess.history.current = _snapshot!(sess)
    return sess
end

"""
    can_undo(sess::DesignSession) -> Bool
    can_redo(sess::DesignSession) -> Bool

Whether there is an edit to undo (an undo step, or an edit not recorded yet) / to redo.
Between edits the answer comes from the history's count of changes
([`EditHistory`](@ref)), without a snapshot of the design.
"""
function can_undo(sess::DesignSession)
    h = sess.history
    isempty(h.undo) || return true
    n = h.changes[]
    (h.current === nothing || n == h.seen) && return false
    isequal(_history_state(sess), h.current) || return true
    h.seen = n                                          # (the design is `current` again)
    return false
end
can_redo(sess::DesignSession) = !isempty(sess.history.redo)

"""
    reset_history!(sess::DesignSession)

Forget the undo and redo steps; the current design becomes the starting point.
"""
function reset_history!(sess::DesignSession)
    h = sess.history
    empty!(h.undo)
    empty!(h.redo)
    h.current = _snapshot!(sess)
    return sess
end
