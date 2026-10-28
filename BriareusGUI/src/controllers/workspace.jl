# Saving and restoring a design: the session as plain data (`session_state`), written to
# and read from a session file (a zip of the TOML and the CAD files it reads, or the TOML
# alone), and the undo/redo history built on the same snapshots. Results (the drive, an
# optimization, a clearance report) are not saved: they are computed again from the
# restored design.

const SESSION_FORMAT = "BriareusGUI session"
const SESSION_VERSION = 2               # 2: `.briar` archives with the CAD files

_strs(xs) = [String(x) for x in xs]
_syms(xs) = Symbol[Symbol(x) for x in xs]
_fdict(d) = Dict{String,Any}(String(k) => Float64(v) for (k, v) in d)

function _beam_state(b::Beam, enabled::Bool, cluster::Symbol, ref::Symbol, refpt, uc)
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
                            "aim_ref" => String(ref), "ref_point" => collect(refpt),
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
    session_state(sess::DesignSession; view = true) -> Dict{String,Any}

The design as plain data, ready for TOML: the target (template, dimensions, element
sizes, windows, parts), the beams (each beam's port, aim, power, focus, spot, cone and
quad, group, whether it is on, and whether its unconverted 1ω and 2ω light takes part),
the settings of the drive (beam model included), the optimization, the clearance check
and the TaLIS checks, and (with `view`) the 3D view's coloring, beam drawing (which
beams' 1ω, 2ω and 3ω light is drawn), cutaway and component display modes. Results are
not part of it; they are
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
                             "beam" => [_beam_state(b, e, c, r, p, u) for (b, e, c, r, p, u) in
                                        zip(bt.beams[], bt.enabled[], bt.cluster[], bt.aim_ref[],
                                            bt.ref_point[], bt.unconverted[])])
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
                             "optimize" => optimize, "review" => review, "talis" => talis)
    if view
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
                                         "side" => dr.side[], "show_off" => v.show_off[])
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
    _cancel_fits!(ts)                     # a capsule fit in flight was for another design
    _try_build!(ts; kind = Symbol(t["kind"]), params, mesh, windows, parts,
                excluded = Set(_syms(get(t, "excluded", String[])))) ||
        throw(ArgumentError("the saved target was rejected: $(ts.error[])"))
    sess.visrad_source[] = get(state, "visrad_source", "")
    b = state["beams"]
    rows = b["beam"]
    set_beams!(bt, [_state_beam(r) for r in rows]; preset = Symbol(b["preset"]),
               enabled = Bool[r["enabled"] for r in rows], cluster = Symbol[Symbol(r["cluster"]) for r in rows],
               aim_ref = Symbol[Symbol(get(r, "aim_ref", "none")) for r in rows],
               ref_point = [Vec3(get(r, "ref_point", [0.0, 0.0, 0.0])...) for r in rows],
               unconverted = [(Bool(get(r, "omega1", false)), Bool(get(r, "omega2", false))) for r in rows])
    _assign!(bt.by, Symbol(b["by"]))
    haskey(b, "facility") && _assign!(bt.facility, Symbol(b["facility"]))
    if haskey(state, "drive")
        d = state["drive"]
        _assign!(sess.drive.spot_radius, Float64(d["spot_radius"]))
        set_cavity!(sess.drive, get(d, "cavity", false); point = Float64.(get(d, "cavity_point", [0.0, 0.0, 0.0])),
                    films = _syms(get(d, "films", String[])))
        # sessions before the beam-model choice deposited with the f-cone
        set_beam_model!(sess.drive, Symbol(get(d, "beam_model", "fcone"));
                        profile = haskey(d, "profile") ? _state_profile(d["profile"]) : VisRadProfile())
        # sessions before the unconverted-light deposition: off
        set_unconverted_ratio!(sess.drive, Float64.(get(d, "unconverted_ratio", [0.0, 0.0]))...)
    end
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
    # the checks' component roles once the saved design's target is in place (a target
    # built on a worker would otherwise give its new components their default roles)
    roles = function ()
        if haskey(state, "review")
            r = state["review"]
            _assign!(sess.review.targets, _syms(r["targets"]))
            _assign!(sess.review.hardware, _syms(r["hardware"]))
        end
        if haskey(state, "talis")
            c = state["talis"]
            _assign!(sess.talis.capture, _syms(c["capture"]))
            _assign!(sess.talis.hardware, _syms(c["hardware"]))
            _assign!(sess.talis.diagnostics, _syms(get(c, "diagnostics", String[])))
        end
    end
    roles()
    _after_build!(roles, ts)
    if haskey(state, "review")
        r, rv = state["review"], sess.review
        _assign!(rv.peak_fraction, Float64(r["peak_fraction"]))
        _assign!(rv.wings, r["wings"])
        _assign!(rv.errors, Dict{Symbol,Float64}(Symbol(k) => v for (k, v) in r["errors"]))
        _assign!(rv.quality, Symbol(r["quality"]))
    end
    if haskey(state, "talis")
        c, tl = state["talis"], sess.talis
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
        # sessions before the plane cut and the axis and range of the sector: a sector
        # about z, unbounded
        haskey(v, "cut_width") && set_cut!(vs; azimuth = v["cut_azimuth"], width = v["cut_width"],
                                           components = _syms(v["cut_components"]),
                                           kind = Symbol(get(v, "cut_kind", "sector")),
                                           axis = Symbol(get(v, "cut_axis", "z")),
                                           range = Tuple(Float64.(get(v, "cut_range", [-Inf, Inf]))))
        # sessions before the plane's offset and line of sight kept the shifted point and
        # no line
        if haskey(v, "cut_normal")
            set_cut_plane!(vs; normal = Float64.(v["cut_normal"]), point = Float64.(v["cut_point"]),
                           offset = Float64(get(v, "cut_shift", 0.0)))
            _assign!(vs.cut_line, Tuple(Float64.(get(v, "cut_line", [0.0, 0.0, 0.0]))))
            _assign!(vs.cut_rotation, Float64(get(v, "cut_rotation", 0.0)))
            _assign!(vs.cut_kind, Symbol(get(v, "cut_kind", "sector")))
        end
        # sessions before display modes listed hidden components: invisible
        disp = Dict{Symbol,Symbol}(Symbol(k) => Symbol(m) for (k, m) in get(v, "display", Dict()))
        for n in get(v, "hidden", String[])
            disp[Symbol(n)] = :invisible
        end
        filter!(p -> last(p) in DISPLAY_MODES && last(p) !== :opaque, disp)
        _assign!(vs.display, disp)
        _assign!(vs.hidden_beams, Set{String}(get(v, "hidden_beams", String[])))
        _assign!(vs.shown_1w, Set{String}(get(v, "shown_1w", String[])))
        _assign!(vs.shown_2w, Set{String}(get(v, "shown_2w", String[])))
        haskey(v, "side") && _assign!(sess.drive.side, Bool(v["side"]))
        haskey(v, "show_off") && _assign!(vs.show_off, Bool(v["show_off"]))
    end
    return sess
end

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

# The CAD files the design reads: its STL and VisRad parts and the VisRad export template.
function _cad_files(sess::DesignSession)
    files = String[p.file for p in sess.target.parts[] if p.kind in (:mesh, :visrad) && !isempty(p.file)]
    isempty(sess.visrad_source[]) || push!(files, sess.visrad_source[])
    return unique!(files)
end

# The state as written to a file: each CAD file by the path it was saved from and its
# content hash. Returns the state and the hash => file of each file found.
function _portable_state(sess::DesignSession)
    state = session_state(sess)
    hashes = Dict{String,String}()
    for f in _cad_files(sess)
        isfile(f) && (hashes[f] = _sha256(f))
    end
    blobs = Dict{String,String}(h => f for (f, h) in hashes)
    for p in state["target"]["parts"]
        f = p["file"]
        isempty(f) && continue
        haskey(hashes, f) && (p["sha256"] = hashes[f])
        p["file"] = cad_source(sess, f)
    end
    v = sess.visrad_source[]
    if !isempty(v)
        haskey(hashes, v) && (state["visrad_sha256"] = hashes[v])
        state["visrad_source"] = cad_source(sess, v)
    end
    return state, blobs
end

# The zip entry of a CAD file: its hash and its source's extension.
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
  several parts is stored once.
- `embed = false` (a `.briar.toml` file): the TOML alone, plain text that can be read
  and edited by hand. CAD parts are referenced by path, so the file travels without
  its geometry (for CAD that may not leave the computer it is on).

Either way each CAD part records the path it was read from and its content's SHA-256,
which [`load_session!`](@ref) checks.
"""
function save_session(path::AbstractString, sess::DesignSession; embed::Bool = !endswith(path, ".toml"))
    state, blobs = _portable_state(sess)
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
        end
        mv(tmp, path; force = true)            # an interrupted save leaves the old file
        mb = sum(filesize, values(blobs); init = 0) / 2^20
        sess.status[] = "Saved the session to $(basename(path))" *
                        (isempty(blobs) ? "" : @sprintf(" with %d CAD file%s (%.1f MB)", length(blobs),
                                                        length(blobs) == 1 ? "" : "s", mb))
    else
        open(io -> _write_toml(io, state), path, "w")
        sess.status[] = "Saved the session to $(basename(path))" *
                        (isempty(blobs) ? "" : " (CAD files by path, not embedded)")
    end
    return path
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
`.briar.toml` file; older sessions too) into the session ([`restore!`](@ref)). The
edit history starts over at the loaded design.

A CAD file embedded in the archive is the one the design reads (extracted to a cache
folder; [`cad_source`](@ref) gives the path it was saved from). When that source still
exists but its content changed since the save, the copy is kept and the source is listed
in `sess.cad_changed`; [`reload_cad_sources!`](@ref) switches to the sources and
[`relink_cad!`](@ref) to another file. A session saved without its CAD reads each
file from its path; the status names those whose content changed since the save.
"""
function load_session!(sess::DesignSession, path::AbstractString)
    bytes = read(path)
    zr = length(bytes) >= 4 && bytes[1:4] == b"PK\x03\x04" ? ZipReader(bytes) : nothing
    state = TOML.parse(zr === nothing ? String(bytes) : zip_readentry(zr, "session.toml", String))
    changed = String[]
    if haskey(state, "target") && haskey(state["target"], "parts")
        for p in state["target"]["parts"]
            isempty(get(p, "file", "")) && continue
            f, c = _resolve_cad!(sess, zr, p["file"], get(p, "sha256", nothing))
            p["file"] = f
            c && push!(changed, f)
        end
    end
    if !isempty(get(state, "visrad_source", ""))
        f, c = _resolve_cad!(sess, zr, state["visrad_source"], get(state, "visrad_sha256", nothing))
        state["visrad_source"] = f
        c && push!(changed, f)
    end
    restore!(sess, state)
    reset_history!(sess)
    unique!(changed)
    copies = filter(f -> haskey(sess.cad_sources, f), changed)
    ondisk = setdiff(changed, copies)
    listed(fs) = join(basename.(cad_source.(Ref(sess), fs)), ", ")
    sess.cad_changed[] = copies
    sess.status[] = "Loaded the session from $(basename(path))" *
                    (isempty(copies) ? "" : "; CAD changed on disk since the save (the design reads " *
                                             "the saved copies): " * listed(copies)) *
                    (isempty(ondisk) ? "" : "; CAD read from disk, changed since the save: " * listed(ondisk))
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
