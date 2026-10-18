# Saving and restoring a design: the session as plain data (`session_state`), written to
# and read from a TOML session file, and the undo/redo history built on the same
# snapshots. Results (the drive, an optimization, a clearance report) are not saved:
# they are computed again from the restored design.

const SESSION_FORMAT = "BriareusGUI session"
const SESSION_VERSION = 1

_strs(xs) = [String(x) for x in xs]
_syms(xs) = Symbol[Symbol(x) for x in xs]
_fdict(d) = Dict{String,Any}(String(k) => Float64(v) for (k, v) in d)

function _beam_state(b::Beam, enabled::Bool, cluster::Symbol, ref::Symbol, refpt)
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
                            "aim_ref" => String(ref), "ref_point" => collect(refpt))
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
    session_state(sess::DesignSession; view = true) -> Dict{String,Any}

The design as plain data, ready for TOML: the target (template, dimensions, element
sizes, windows, parts), the beams (each beam's port, aim, power, focus, spot, cone and
quad, group, and whether it is on), the settings of the drive, the optimization, the
clearance check and the TaLIS checks, and (with `view`) the 3D view's coloring, beam
drawing, cutaway and hidden components. Results are not part of it; they are
computed again. Restore it with [`restore!`](@ref).
"""
function session_state(sess::DesignSession; view::Bool = true)
    ts, bt, dr, op, rv, tl = sess.target, sess.beams, sess.drive, sess.optimize, sess.review, sess.talis
    target = Dict{String,Any}(
        "kind" => String(ts.kind[]), "params" => _fdict(ts.params[]),
        "mesh" => Dict{String,Any}("h" => ts.mesh[].h, "capsule_h" => ts.mesh[].capsule_h),
        "windows" => [Dict{String,Any}("phi" => w.phi, "z" => w.z, "width" => w.width,
                                       "height" => w.height, "open" => w.open, "albedo" => w.albedo,
                                       "conversion" => w.conversion, "round" => w.round)
                      for w in ts.windows[]],
        "excluded" => sort(_strs(ts.excluded[])),
        "parts" => [Dict{String,Any}("kind" => String(p.kind), "name" => String(p.name),
                                     "values" => _fdict(p.values), "two_sided" => p.two_sided,
                                     "inward" => p.inward, "aperture" => p.aperture, "file" => p.file,
                                     "anchor" => String(p.anchor), "coords" => String(p.coords),
                                     "object" => p.object)
                    for p in ts.parts[]])
    beams = Dict{String,Any}("preset" => String(bt.preset[]), "by" => String(bt.by[]),
                             "facility" => String(bt.facility[]),
                             "beam" => [_beam_state(b, e, c, r, p) for (b, e, c, r, p) in
                                        zip(bt.beams[], bt.enabled[], bt.cluster[], bt.aim_ref[],
                                            bt.ref_point[])])
    optimize = Dict{String,Any}(
        "by" => String(op.by[]), "fields" => _strs(op.fields[]), "deltas" => _fdict(op.deltas[]),
        "objective" => String(op.objective[]),
        "weights" => Dict{String,Any}(string(l) => w for (l, w) in op.weights[]),
        "hold_drive" => op.hold_drive[], "keepout" => op.keepout[], "keepout_min" => op.keepout_min[],
        "capsule_keepout" => op.capsule_keepout[], "capsule_min" => op.capsule_min[],
        "spill" => op.spill[], "iterations" => op.iterations[],
        "max_evaluations" => op.max_evaluations[])
    review = Dict{String,Any}(
        "targets" => _strs(rv.targets[]), "hardware" => _strs(rv.hardware[]),
        "peak_fraction" => rv.peak_fraction[], "wings" => rv.wings[], "errors" => _fdict(rv.errors[]),
        "quality" => String(rv.quality[]))
    talis = Dict{String,Any}(
        "capture" => _strs(tl.capture[]), "hardware" => _strs(tl.hardware[]),
        "tolerance" => tl.tolerance[], "pivot" => String(tl.pivot[]),
        "pivot_point" => collect(tl.pivot_point[]), "omega1" => tl.omega1[])
    state = Dict{String,Any}("format" => SESSION_FORMAT, "version" => SESSION_VERSION,
                             "visrad_source" => sess.visrad_source[],
                             "target" => target, "beams" => beams,
                             "drive" => Dict{String,Any}("spot_radius" => dr.spot_radius[]),
                             "optimize" => optimize, "review" => review, "talis" => talis)
    if view
        v = sess.view
        state["view"] = Dict{String,Any}("color" => String(v.color[]), "beams" => String(v.beams[]),
                                         "beam_color" => String(v.beam_color[]),
                                         "beam_length" => v.beam_length[], "cutaway" => v.cutaway[],
                                         "cut_azimuth" => v.cut_azimuth[], "cut_width" => v.cut_width[],
                                         "cut_components" => sort(_strs(v.cut_components[])),
                                         "hidden" => sort(_strs(v.hidden[])),
                                         "side" => dr.side[])
    end
    return state
end

"""
    restore!(sess::DesignSession, state) -> sess

Set the session to a [`session_state`](@ref) (a view section, when present, restores
the view too). The target is rebuilt once; analyses whose inputs changed become stale
as after any edit. Throws `ArgumentError` for a state of another format, or one the
geometry rejects (the session is then left as it was before the target, if the
target fails, or partly restored otherwise).
"""
function restore!(sess::DesignSession, state::AbstractDict)
    get(state, "format", "") == SESSION_FORMAT ||
        throw(ArgumentError("not a $SESSION_FORMAT (format \"$(get(state, "format", ""))\")"))
    get(state, "version", 0) <= SESSION_VERSION ||
        throw(ArgumentError("session version $(state["version"]) is newer than this program's $SESSION_VERSION"))
    ts, bt = sess.target, sess.beams
    t = state["target"]
    params = copy(ts.params[])
    for (k, v) in t["params"]
        params[Symbol(k)] = v
    end
    mesh = MeshSpec(t["mesh"]["h"], t["mesh"]["capsule_h"])
    windows = WindowSpec[WindowSpec(; phi = w["phi"], z = w["z"], width = w["width"], height = w["height"],
                                    open = w["open"], albedo = w["albedo"], conversion = w["conversion"],
                                    round = w["round"]) for w in t["windows"]]
    parts = PartSpec[PartSpec(Symbol(p["kind"]); name = Symbol(p["name"]),
                              values = Dict{Symbol,Float64}(Symbol(k) => v for (k, v) in p["values"]),
                              two_sided = p["two_sided"], inward = p["inward"], aperture = p["aperture"],
                              file = get(p, "file", ""), anchor = Symbol(get(p, "anchor", "tcc")),
                              coords = Symbol(get(p, "coords", "cartesian")),
                              object = get(p, "object", ""))
                     for p in t["parts"]]
    _try_build!(ts; kind = Symbol(t["kind"]), params, mesh, windows, parts,
                excluded = Set(_syms(get(t, "excluded", String[])))) ||
        throw(ArgumentError("the saved target was rejected: $(ts.error[])"))
    sess.visrad_source[] = get(state, "visrad_source", "")
    b = state["beams"]
    rows = b["beam"]
    set_beams!(bt, [_state_beam(r) for r in rows]; preset = Symbol(b["preset"]),
               enabled = Bool[r["enabled"] for r in rows], cluster = Symbol[Symbol(r["cluster"]) for r in rows],
               aim_ref = Symbol[Symbol(get(r, "aim_ref", "none")) for r in rows],
               ref_point = [Vec3(get(r, "ref_point", [0.0, 0.0, 0.0])...) for r in rows])
    _assign!(bt.by, Symbol(b["by"]))
    haskey(b, "facility") && _assign!(bt.facility, Symbol(b["facility"]))
    haskey(state, "drive") && _assign!(sess.drive.spot_radius, Float64(state["drive"]["spot_radius"]))
    if haskey(state, "optimize")
        o, op = state["optimize"], sess.optimize
        _assign!(op.by, Symbol(o["by"]))
        _assign!(op.fields, _syms(o["fields"]))
        _assign!(op.deltas, Dict{Symbol,Float64}(Symbol(k) => v for (k, v) in o["deltas"]))
        _assign!(op.objective, Symbol(o["objective"]))
        _assign!(op.weights, Dict{Int,Float64}(parse(Int, k) => v for (k, v) in o["weights"]))
        for k in (:hold_drive, :keepout, :keepout_min, :capsule_keepout, :capsule_min, :spill,
                  :iterations, :max_evaluations)
            obs = getfield(op, k)
            _assign!(obs, convert(typeof(obs[]), o[String(k)]))
        end
    end
    if haskey(state, "review")
        r, rv = state["review"], sess.review
        _assign!(rv.targets, _syms(r["targets"]))
        _assign!(rv.hardware, _syms(r["hardware"]))
        _assign!(rv.peak_fraction, Float64(r["peak_fraction"]))
        _assign!(rv.wings, r["wings"])
        _assign!(rv.errors, Dict{Symbol,Float64}(Symbol(k) => v for (k, v) in r["errors"]))
        _assign!(rv.quality, Symbol(r["quality"]))
    end
    if haskey(state, "talis")
        c, tl = state["talis"], sess.talis
        _assign!(tl.capture, _syms(c["capture"]))
        _assign!(tl.hardware, _syms(c["hardware"]))
        _assign!(tl.tolerance, Float64(c["tolerance"]))
        _assign!(tl.pivot, Symbol(c["pivot"]))
        _assign!(tl.pivot_point, Tuple(Float64.(c["pivot_point"])))
        _assign!(tl.omega1, c["omega1"])
    end
    if haskey(state, "view")
        v, vs = state["view"], sess.view
        _assign!(vs.color, Symbol(v["color"]))
        _assign!(vs.beams, Symbol(v["beams"]))
        _assign!(vs.beam_color, Symbol(v["beam_color"]))
        _assign!(vs.beam_length, Float64(v["beam_length"]))
        _assign!(vs.cutaway, v["cutaway"])
        haskey(v, "cut_width") && set_cut!(vs; azimuth = v["cut_azimuth"], width = v["cut_width"],
                                           components = _syms(v["cut_components"]))
        _assign!(vs.hidden, Set(_syms(v["hidden"])))
        haskey(v, "side") && _assign!(sess.drive.side, Bool(v["side"]))
    end
    return sess
end

"""
    save_session(path, sess::DesignSession) -> path

Write the design ([`session_state`](@ref), view included) to a TOML file (by
convention `.briareus.toml`). The file is plain text: beams, windows and parts can be
read and edited by hand.
"""
function save_session(path::AbstractString, sess::DesignSession)
    state = session_state(sess)
    open(path, "w") do io
        println(io, "# ", SESSION_FORMAT, ": a Briareus design (target, beams, analysis settings).")
        println(io, "# Lengths in mm, power in TW, angles in degrees; results are not saved.")
        TOML.print(io, state; sorted = true)
    end
    sess.status[] = "Saved the session to $(basename(path))"
    return path
end

"""
    load_session!(sess::DesignSession, path) -> sess

Read a session file written by [`save_session`](@ref) into the session
([`restore!`](@ref)). The edit history starts over at the loaded design.
"""
function load_session!(sess::DesignSession, path::AbstractString)
    state = TOML.parsefile(path)
    restore!(sess, state)
    reset_history!(sess)
    sess.status[] = "Loaded the session from $(basename(path))"
    return sess
end

# ---------------------------------------------------------------- undo and redo

"""
    record!(sess::DesignSession) -> Bool

Note the design's current state in the edit history: if it differs from the last one
noted, the last one becomes an undo step (and the redo steps are dropped). Returns
whether a step was added. A window calls this after every event, so any edit (a field,
a preset, an applied optimization) can be undone.
"""
function record!(sess::DesignSession)
    h = sess.history
    state = session_state(sess; view = false)
    if h.current === nothing
        h.current = state
        return false
    end
    isequal(state, h.current) && return false
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

function _step_to!(sess::DesignSession, state)
    restore!(sess, state)
    sess.history.current = session_state(sess; view = false)
    return sess
end

"Whether there is an edit to undo / to redo."
can_undo(sess::DesignSession) = !isempty(sess.history.undo) ||
                                (sess.history.current !== nothing &&
                                 !isequal(session_state(sess; view = false), sess.history.current))
can_redo(sess::DesignSession) = !isempty(sess.history.redo)

"""
    reset_history!(sess::DesignSession)

Forget the undo and redo steps; the current design becomes the starting point.
"""
function reset_history!(sess::DesignSession)
    h = sess.history
    empty!(h.undo)
    empty!(h.redo)
    h.current = session_state(sess; view = false)
    return sess
end
