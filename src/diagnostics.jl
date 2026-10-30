# Diagnostic definitions: one TOML file per diagnostic, read from a search path (the
# built-ins in data/diagnostics, then the directories of a user configuration file and of
# BRIAREUS_DIAGNOSTICS), checked against the schema of its kind, resolved through
# `extends` and built into an instrument. Design: reference/diagnostics-design.md; public
# values: research notes, "Dante-family hardware and geometry: what is public".

"""
    DIAGNOSTIC_SCHEMA

The version of the diagnostic definition format that this release reads (the `schema`
key of every definition file).
"""
const DIAGNOSTIC_SCHEMA = 1

const _BUILTIN_DIAGNOSTICS = joinpath(@__DIR__, "..", "data", "diagnostics")

# ---- schema ------------------------------------------------------------------------------

"""
    DiagnosticField(type; required = false, doc = "")

One key of a diagnostic kind's schema ([`register_diagnostic_kind`](@ref)): `type` is
`Real`, `Integer`, `String`, `Bool`, `:date`, `Vector{String}`, `Vector{Real}`, a
schema table (a `Dict{String,DiagnosticField}`) or a one-element vector holding a
schema table (an array of tables). Keys ending in `_csv` name a CSV file, a path
relative to the definition file.
"""
struct DiagnosticField
    type::Any
    required::Bool
    doc::String
end
DiagnosticField(type; required::Bool = false, doc::AbstractString = "") =
    DiagnosticField(type, required, String(doc))

# keys every table may carry: provenance
const _PROVENANCE_KEYS = Dict("source" => DiagnosticField(String; doc = "citation"),
                              "confidence" => DiagnosticField(String; doc = "H, M or L"),
                              "note" => DiagnosticField(String))

# top-level keys of every definition; those below the line are not inherited through
# `extends` (they belong to the file that states them)
const _HEADER_KEYS = Dict(
    "schema" => DiagnosticField(Integer; required = true, doc = "format version"),
    "id" => DiagnosticField(String; required = true, doc = "unique name, e.g. nif_dante1"),
    "kind" => DiagnosticField(String; doc = "a registered kind; inherited through extends"),
    "name" => DiagnosticField(String; doc = "display name"),
    "facility" => DiagnosticField(String),
    "extends" => DiagnosticField(String; doc = "the id this definition modifies"),
    "public" => DiagnosticField(Bool; doc = "false for non-public detail"),
    "valid_from" => DiagnosticField(:date),
    "valid_to" => DiagnosticField(:date),
    "shots" => DiagnosticField(Vector{String}))
const _OWN_KEYS = ("schema", "id", "name", "extends", "public", "valid_from", "valid_to", "shots",
                   "source", "confidence", "note")

struct _Kind
    schema::Dict{String,DiagnosticField}
    build::Any
end

const _DIAGNOSTIC_KINDS = Dict{Symbol,_Kind}()

"""
    register_diagnostic_kind(kind::Symbol, schema::Dict{String,DiagnosticField}, build)

Make `kind` a diagnostic kind: definition files with `kind = "<kind>"` are checked
against `schema` (the keys besides the header ones: `schema`, `id`, `kind`, `name`,
`facility`, `extends`, `public`, `valid_from`, `valid_to`, `shots`, and `source`,
`confidence`, `note` in every table) and built by `build(def::DiagnosticDefinition;
kwargs...)`. A kind registered again replaces the earlier one.
"""
function register_diagnostic_kind(kind::Symbol, schema::Dict{String,DiagnosticField}, build)
    _DIAGNOSTIC_KINDS[kind] = _Kind(schema, build)
    return kind
end

"""
    DiagnosticIssue

A problem found in a definition file ([`check_diagnostics`](@ref)): `file`, `key` (the
path to the value, e.g. `channel[6].filters[1].thickness_um`), `message`, and `level`:
`:error` (the file cannot be used), `:missing` (a value the definition lacks to be
built, e.g. one that is not public) or `:warning`.
"""
struct DiagnosticIssue
    file::String
    key::String
    message::String
    level::Symbol
end

Base.show(io::IO, i::DiagnosticIssue) =
    print(io, i.level, ": ", i.file, isempty(i.key) ? "" : " [" * i.key * "]", ": ", i.message)

_typename(t) = t === :date ? "a date" : t isa Dict ? "a table" :
               t isa Vector ? "an array of tables" : t === Real ? "a number" :
               t === Integer ? "an integer" : t === String ? "a string" : t === Bool ? "true or false" :
               t === Vector{String} ? "an array of strings" : t === Vector{Real} ? "an array of numbers" :
               string(t)

function _check_value(v, t)
    t === :date && return v isa Dates.Date
    t === Real && return v isa Real && !(v isa Bool)
    t === Integer && return v isa Integer && !(v isa Bool)
    t === String && return v isa AbstractString
    t === Bool && return v isa Bool
    t === Vector{String} && return v isa AbstractVector && all(x -> x isa AbstractString, v)
    t === Vector{Real} && return v isa AbstractVector && all(x -> x isa Real && !(x isa Bool), v)
    t isa Dict && return v isa AbstractDict
    t isa Vector && return v isa AbstractVector && all(x -> x isa AbstractDict, v)
    return false
end

# Check table `d` against `spec`; push issues. `required` is checked only when asked
# (a file that extends another need not repeat what it inherits).
function _check_table!(issues, file, path, d::AbstractDict, spec::Dict{String,DiagnosticField};
                       required::Bool = true)
    for (k, v) in d
        f = get(spec, k, get(_PROVENANCE_KEYS, k, nothing))
        key = isempty(path) ? k : path * "." * k
        if f === nothing
            push!(issues, DiagnosticIssue(file, key, "unknown key (allowed: " *
                                          join(sort(collect(keys(spec))), ", ") * ")", :error))
        elseif !_check_value(v, f.type)
            push!(issues, DiagnosticIssue(file, key, "must be " * _typename(f.type), :error))
        elseif f.type isa Dict
            _check_table!(issues, file, key, v, f.type; required)
        elseif f.type isa Vector
            for (i, x) in enumerate(v)
                tag = haskey(x, "number") ? "[$(x["number"])]" : "[$i]"
                _check_table!(issues, file, key * tag, x, only(f.type); required)
            end
        end
    end
    if required
        for (k, f) in spec
            f.required && !haskey(d, k) &&
                push!(issues, DiagnosticIssue(file, isempty(path) ? k : path * "." * k, "required", :missing))
        end
    end
    return issues
end

# ---- definitions ------------------------------------------------------------------------

"""
    DiagnosticDefinition

A diagnostic definition read from a file and resolved through `extends`: `id`, `name`,
`kind`, `file` (the file that defines it), `chain` (the files it was resolved from,
the base first), `public` (`false` when any file in the chain sets `public = false`) and `data` (the
resolved table, CSV paths made absolute). [`diagnostic`](@ref) builds the instrument.
"""
struct DiagnosticDefinition
    id::String
    name::String
    kind::Symbol
    file::String
    chain::Vector{String}
    public::Bool
    data::Dict{String,Any}
end

Base.show(io::IO, d::DiagnosticDefinition) =
    print(io, "DiagnosticDefinition(\"", d.id, "\", ", d.kind, ", ", basename(d.file),
          length(d.chain) > 1 ? " extending " * join(basename.(d.chain[1:end-1]), " ← ") : "", ")")

# a file as read: its table with CSV paths made absolute
struct _Raw
    file::String
    data::Dict{String,Any}
end

function _absolute_csv!(d::AbstractDict, dir)
    for (k, v) in d
        if endswith(k, "_csv") && v isa AbstractString
            d[k] = normpath(isabspath(v) ? v : joinpath(dir, v))
        elseif v isa AbstractDict
            _absolute_csv!(v, dir)
        elseif v isa AbstractVector
            foreach(x -> x isa AbstractDict && _absolute_csv!(x, dir), v)
        end
    end
    return d
end

function _read_raw(file::AbstractString, issues)
    d = try
        TOML.parsefile(file)
    catch e
        push!(issues, DiagnosticIssue(file, "", "not valid TOML: " * sprint(showerror, e), :error))
        return nothing
    end
    for (k, f) in _HEADER_KEYS
        haskey(d, k) && !_check_value(d[k], f.type) &&
            (push!(issues, DiagnosticIssue(file, k, "must be " * _typename(f.type), :error)); return nothing)
        f.required && !haskey(d, k) &&
            (push!(issues, DiagnosticIssue(file, k, "required", :error)); return nothing)
    end
    d["schema"] == DIAGNOSTIC_SCHEMA ||
        (push!(issues, DiagnosticIssue(file, "schema", "this release reads schema $DIAGNOSTIC_SCHEMA", :error));
         return nothing)
    haskey(d, "kind") || haskey(d, "extends") ||
        (push!(issues, DiagnosticIssue(file, "kind", "required (or extends)", :error)); return nothing)
    return _Raw(String(file), _absolute_csv!(d, dirname(abspath(file))))
end

# merge child into parent: tables key by key, channels by number, other values replaced
function _merge(parent::AbstractDict, child::AbstractDict)
    out = deepcopy(Dict{String,Any}(parent))
    for k in _OWN_KEYS
        delete!(out, k)
    end
    for (k, v) in child
        k in ("extends", "remove_channels") && continue
        if k == "channel" && haskey(out, "channel")
            # a channel of the same number takes the child's keys (each value replaced
            # whole: filters, mirror); a new number is added
            byn = Dict(c["number"] => c for c in out["channel"] if haskey(c, "number"))
            for c in v
                n = get(c, "number", 0)
                byn[n] = haskey(byn, n) ? merge(byn[n], deepcopy(Dict{String,Any}(c))) : deepcopy(c)
            end
            out[k] = Any[byn[n] for n in sort(collect(keys(byn)))]
        elseif v isa AbstractDict && get(out, k, nothing) isa AbstractDict
            out[k] = _merge_table(out[k], v)
        else
            out[k] = deepcopy(v)
        end
    end
    if haskey(child, "remove_channels") && haskey(out, "channel")
        drop = Set(child["remove_channels"])
        out["channel"] = Any[c for c in out["channel"] if !(get(c, "number", 0) in drop)]
    end
    return out
end

function _merge_table(a::AbstractDict, b::AbstractDict)
    out = deepcopy(Dict{String,Any}(a))
    for (k, v) in b
        out[k] = v isa AbstractDict && get(out, k, nothing) isa AbstractDict ? _merge_table(out[k], v) : deepcopy(v)
    end
    return out
end

function _resolve(id, raws::Dict{String,_Raw}, issues, seen = String[])
    r = raws[id]
    if id in seen
        push!(issues, DiagnosticIssue(r.file, "extends", "a cycle: " * join([seen; id], " → "), :error))
        return nothing
    end
    d = r.data
    if haskey(d, "extends")
        p = d["extends"]
        if !haskey(raws, p)
            push!(issues, DiagnosticIssue(r.file, "extends", "no definition \"$p\" on the search path", :error))
            return nothing
        end
        base = _resolve(p, raws, issues, [seen; id])
        base === nothing && return nothing
        data = _merge(base.data, d)
        chain = [base.chain; r.file]
        public = base.public && get(d, "public", true)
        name = get(d, "name", base.name)
    else
        data = deepcopy(d)
        delete!(data, "remove_channels")
        chain = [r.file]
        public = get(d, "public", true)
        name = get(d, "name", d["id"])
    end
    for k in ("id", "schema")
        data[k] = d[k]
    end
    return DiagnosticDefinition(d["id"], name, Symbol(data["kind"]), r.file, chain, public, data)
end

# ---- search path and configuration -----------------------------------------------------

"""
    diagnostic_config_path() -> String

The user configuration file that stores the diagnostic search path (key
`diagnostics.paths`): `BRIAREUS_CONFIG` if set, else `%APPDATA%\\Briareus\\config.toml` on
Windows, `~/Library/Application Support/Briareus/config.toml` on macOS and
`\$XDG_CONFIG_HOME/briareus/config.toml` (default `~/.config`) elsewhere.
"""
function diagnostic_config_path()
    haskey(ENV, "BRIAREUS_CONFIG") && return ENV["BRIAREUS_CONFIG"]
    Sys.iswindows() && return joinpath(get(ENV, "APPDATA", homedir()), "Briareus", "config.toml")
    Sys.isapple() && return joinpath(homedir(), "Library", "Application Support", "Briareus", "config.toml")
    return joinpath(get(ENV, "XDG_CONFIG_HOME", joinpath(homedir(), ".config")), "briareus", "config.toml")
end

_read_config(file = diagnostic_config_path()) = isfile(file) ? TOML.parsefile(file) : Dict{String,Any}()

function _saved_paths()
    c = _read_config()
    p = get(get(c, "diagnostics", Dict{String,Any}()), "paths", String[])
    return String[x for x in p if x isa AbstractString]
end

const _SESSION_PATHS = String[]
const _PATHS_LOCK = ReentrantLock()

"""
    diagnostic_paths(; builtin = true) -> Vector{String}

The diagnostic search path, in order of precedence (a later definition of an `id`
shadows an earlier one): the built-in definitions (`data/diagnostics`, unless `builtin =
false`), the directories saved in the configuration file
([`diagnostic_config_path`](@ref)), those of the environment variable
`BRIAREUS_DIAGNOSTICS` (separated by `;` on Windows, `:` elsewhere), and those added for
this Julia session ([`add_diagnostic_path`](@ref) with `save = false`).
"""
function diagnostic_paths(; builtin::Bool = true)
    env = filter(!isempty, split(get(ENV, "BRIAREUS_DIAGNOSTICS", ""), Sys.iswindows() ? ';' : ':'))
    p = String[]
    builtin && push!(p, normpath(_BUILTIN_DIAGNOSTICS))
    append!(p, _saved_paths())
    append!(p, String.(env))
    lock(() -> append!(p, _SESSION_PATHS), _PATHS_LOCK)
    return unique(p)
end

function _write_config(c)
    file = diagnostic_config_path()
    mkpath(dirname(file))
    open(io -> TOML.print(io, c; sorted = true), file, "w")
    return file
end

"""
    add_diagnostic_path(dir; save = true) -> Vector{String}

Add the directory `dir` to the diagnostic search path, after those already on it:
`save = true` stores it in the configuration file ([`diagnostic_config_path`](@ref)),
so later sessions find it; `save = false` adds it for this Julia session only. Returns
the search path ([`diagnostic_paths`](@ref)).
"""
function add_diagnostic_path(dir::AbstractString; save::Bool = true)
    isdir(dir) || throw(ArgumentError("no such directory: $dir"))
    d = normpath(abspath(dir))
    if save
        c = _read_config()
        t = get!(c, "diagnostics", Dict{String,Any}())
        p = String[x for x in get(t, "paths", String[])]
        d in p || push!(p, d)
        t["paths"] = p
        _write_config(c)
    else
        lock(() -> (d in _SESSION_PATHS || push!(_SESSION_PATHS, d)), _PATHS_LOCK)
    end
    return diagnostic_paths()
end

"""
    remove_diagnostic_path(dir; save = true) -> Vector{String}

Remove `dir` from the saved search path (`save = true`) or from this session's
([`add_diagnostic_path`](@ref)). Returns the search path.
"""
function remove_diagnostic_path(dir::AbstractString; save::Bool = true)
    d = normpath(abspath(dir))
    if save
        c = _read_config()
        t = get(c, "diagnostics", nothing)
        if t !== nothing && haskey(t, "paths")
            t["paths"] = String[x for x in t["paths"] if normpath(x) != d]
            _write_config(c)
        end
    else
        lock(() -> filter!(!=(d), _SESSION_PATHS), _PATHS_LOCK)
    end
    return diagnostic_paths()
end

# the definition files of a path entry: a file, or the .toml files under a directory
function _definition_files(p::AbstractString)
    isfile(p) && return [p]
    isdir(p) || return String[]
    out = String[]
    for (root, dirs, files) in walkdir(p)
        sort!(dirs)
        append!(out, joinpath(root, f) for f in sort(files) if endswith(lowercase(f), ".toml"))
    end
    return out
end

# read every file on `paths`; later ids shadow earlier ones (a warning names both),
# two files in one entry with one id are an error
function _read_all(paths, issues)
    raws = Dict{String,_Raw}()
    for p in paths
        isfile(p) || isdir(p) ||
            (push!(issues, DiagnosticIssue(p, "", "no such file or directory on the search path", :warning));
             continue)
        here = Dict{String,String}()
        for f in _definition_files(p)
            r = _read_raw(f, issues)
            r === nothing && continue
            id = r.data["id"]
            if haskey(here, id)
                push!(issues, DiagnosticIssue(f, "id", "\"$id\" is also defined in $(here[id])", :error))
                continue
            end
            here[id] = f
            haskey(raws, id) &&
                push!(issues, DiagnosticIssue(f, "id", "\"$id\" shadows $(raws[id].file)", :warning))
            raws[id] = r
        end
    end
    return raws
end

"""
    diagnostics(; paths = diagnostic_paths()) -> Dict{String,DiagnosticDefinition}

Every diagnostic definition on the search path (directories, read recursively, or
single files), resolved through `extends`, by `id`. A later definition of an `id`
shadows an earlier one, with a warning naming both files; a file with errors is left
out with a warning ([`check_diagnostics`](@ref) reports them all).
"""
function diagnostics(; paths = diagnostic_paths())
    issues = DiagnosticIssue[]
    out = _definitions(paths, issues)
    for i in issues
        i.level === :missing && continue
        @warn sprint(show, i)
    end
    return out
end

function _definitions(paths, issues)
    raws = _read_all(paths, issues)
    out = Dict{String,DiagnosticDefinition}()
    for id in sort(collect(keys(raws)))
        n = length(issues)
        d = _resolve(id, raws, issues)
        d === nothing && continue
        if !haskey(_DIAGNOSTIC_KINDS, d.kind)
            push!(issues, DiagnosticIssue(d.file, "kind", "unknown kind \"$(d.kind)\" (registered: " *
                                          join(sort(string.(keys(_DIAGNOSTIC_KINDS))), ", ") * ")", :error))
            continue
        end
        spec = _DIAGNOSTIC_KINDS[d.kind].schema
        body = Dict(k => v for (k, v) in d.data if !haskey(_HEADER_KEYS, k) && !haskey(_PROVENANCE_KEYS, k) &&
                    k != "remove_channels")
        # the file's own keys first (so an error names it), then the resolved whole
        own = Dict(k => v for (k, v) in raws[id].data if !haskey(_HEADER_KEYS, k) &&
                   !haskey(_PROVENANCE_KEYS, k))
        _check_table!(issues, d.file, "", own, spec; required = false)
        any(i -> i.level === :error, @view issues[n+1:end]) && continue
        _check_table!(issues, d.file, "", body, spec; required = true)
        any(i -> i.level in (:error, :missing), @view issues[n+1:end]) && continue
        out[id] = d
    end
    return out
end

"""
    diagnostic(id; paths = diagnostic_paths(), kwargs...)
    diagnostic(def::DiagnosticDefinition; kwargs...)

Build the instrument of definition `id` found on `paths` ([`diagnostics`](@ref)), or of
`def`: for `kind = "xray_power"` a [`DanteView`](@ref). `kwargs` go to the kind's
builder (for `xray_power`, those of `DanteView`, e.g. `pixels`), replacing the
definition's values. A definition that lacks a value the instrument needs (e.g. a
channel set that is not public) throws, naming the key to supply.
"""
function diagnostic(id::AbstractString; paths = diagnostic_paths(), kwargs...)
    issues = DiagnosticIssue[]
    defs = _definitions(paths, issues)
    haskey(defs, id) && return diagnostic(defs[id]; kwargs...)
    mine = [i for i in issues if i.level !== :warning]
    throw(ArgumentError("no usable diagnostic \"$id\" on the search path" *
                        (isempty(mine) ? "" : ":\n" * join(sprint.(show, mine), "\n"))))
end

diagnostic(def::DiagnosticDefinition; kwargs...) = _DIAGNOSTIC_KINDS[def.kind].build(def; kwargs...)

"""
    check_diagnostics(; paths = diagnostic_paths()) -> Vector{DiagnosticIssue}

Read every definition on `paths`, resolve `extends`, check each against its kind's
schema and build it, and return the problems found ([`DiagnosticIssue`](@ref)): errors
(unknown keys, wrong types, unreadable CSV files, unknown materials), values a
definition lacks to be built (`:missing`) and warnings (shadowed ids, CSV tables that do
not cover a channel's energy range). An empty result means every definition builds.
"""
function check_diagnostics(; paths = diagnostic_paths())
    issues = DiagnosticIssue[]
    defs = _definitions(paths, issues)
    for id in sort(collect(keys(defs)))
        d = defs[id]
        try
            diagnostic(d)
        catch e
            e isa _MissingValue ? push!(issues, DiagnosticIssue(d.file, e.key, e.msg, :missing)) :
                push!(issues, DiagnosticIssue(d.file, e isa _KeyedError ? e.key : "",
                                              e isa _KeyedError ? e.msg : sprint(showerror, e), :error))
        end
        _coverage!(issues, d)
    end
    return issues
end

# errors a builder raises with the key that caused them
struct _MissingValue <: Exception
    key::String
    msg::String
end
struct _KeyedError <: Exception
    key::String
    msg::String
end
Base.showerror(io::IO, e::_MissingValue) = print(io, "[", e.key, "] ", e.msg)
Base.showerror(io::IO, e::_KeyedError) = print(io, "[", e.key, "] ", e.msg)

# CSV tables that do not cover 50 eV – 20 keV (the unfold's range)
function _coverage!(issues, d::DiagnosticDefinition)
    for c in get(d.data, "channel", Any[])
        for (k, f) in _csv_keys(c)
            isfile(f) || continue
            t = try TabulatedCurve(f) catch; continue end
            (t.E[1] > 50 || t.E[end] < 20_000) &&
                push!(issues, DiagnosticIssue(d.file, "channel[$(get(c, "number", "?"))].$k",
                                              "covers $(t.E[1])–$(t.E[end]) eV, not 50 eV – 20 keV; " *
                                              (k == "response_csv" ? "zero outside" : "end values held outside"),
                                              :warning))
        end
    end
    return issues
end

function _csv_keys(d::AbstractDict, path = "")
    out = Pair{String,String}[]
    for (k, v) in d
        key = isempty(path) ? k : path * "." * k
        if endswith(k, "_csv") && v isa AbstractString
            push!(out, key => v)
        elseif v isa AbstractDict
            append!(out, _csv_keys(v, key))
        elseif v isa AbstractVector
            for (i, x) in enumerate(v)
                x isa AbstractDict && append!(out, _csv_keys(x, key * "[$i]"))
            end
        end
    end
    return out
end

# ---- kind "xray_power": broadband x-ray power diagnostics (Dante, Sentinel) ------------

const _FILTER_SPEC = Dict(
    "material" => DiagnosticField(String; doc = "a filter_transmission material, e.g. Mg, polyimide"),
    "thickness_um" => DiagnosticField(Real; doc = "µm"),
    "density_g_cm3" => DiagnosticField(Real; doc = "g/cm³; default the material's"),
    "transmission_csv" => DiagnosticField(String; doc = "measured transmission (energy eV, fraction)"))

const _MIRROR_SPEC = Dict(
    "kind" => DiagnosticField(String; doc = "SiO2 (grazing incidence) or mlm (Rubery 2023's Cr/Sc multilayer)"),
    "angle_deg" => DiagnosticField(Real; doc = "grazing angle, degrees"),
    "reflectivity_csv" => DiagnosticField(String; doc = "measured reflectivity (energy eV, fraction)"))

const _CHANNEL_SPEC = Dict(
    "number" => DiagnosticField(Integer; required = true),
    "energy_eV" => DiagnosticField(Real; required = true, doc = "nominal band energy"),
    "solid_angle_sr" => DiagnosticField(Real; required = true),
    "filters" => DiagnosticField([_FILTER_SPEC]),
    "grid" => DiagnosticField(Real; doc = "open fraction of the grid or sieve; default 1"),
    "mirror" => DiagnosticField(_MIRROR_SPEC),
    "cathode" => DiagnosticField(String; doc = "XRD cathode: Al, Cr or Ni"),
    "sensitivity_csv" => DiagnosticField(String; doc = "measured detector sensitivity (energy eV, V/GW)"),
    "response_csv" => DiagnosticField(String; doc = "measured whole-channel response (energy eV, V/GW)"))

const _XRAY_POWER_SPEC = Dict(
    "line_of_sight" => DiagnosticField(Dict(
        "theta_deg" => DiagnosticField(Real; required = true),
        "phi_deg" => DiagnosticField(Real; required = true)); required = true),
    "view" => DiagnosticField(Dict(
        "through" => DiagnosticField(String; doc = "the aperture component the view looks through"),
        "pixels" => DiagnosticField(Integer; doc = "rays across the aperture, per side"),
        "area_mm2" => DiagnosticField(Real; doc = "the area Tr is quoted over; default the aperture's"))),
    "recording" => DiagnosticField(Dict(
        "mode" => DiagnosticField(String; doc = "history (a voltage history per channel)"))),
    "channel" => DiagnosticField([_CHANNEL_SPEC]; required = true),
    "remove_channels" => DiagnosticField(Vector{Real}; doc = "channel numbers to drop from the base"))

function _curve(path, key; extrapolate = :hold)
    try
        return TabulatedCurve(path; extrapolate)
    catch e
        throw(_KeyedError(key, "cannot read $path: " * sprint(showerror, e)))
    end
end

function _channel(c::AbstractDict)
    n = c["number"]
    key(k) = "channel[$n].$k"
    Ω = Float64(c["solid_angle_sr"])
    E0 = c["energy_eV"]
    if haskey(c, "response_csv")
        for k in ("filters", "grid", "mirror", "cathode", "sensitivity_csv")
            haskey(c, k) && throw(_KeyedError(key(k), "a channel with response_csv takes no $k"))
        end
        return DanteResponse(n, E0; solid_angle = Ω,
                             response = _curve(c["response_csv"], key("response_csv"); extrapolate = :zero))
    end
    filters = Pair{Symbol,Any}[]
    tables = TabulatedCurve[]
    for (i, f) in enumerate(get(c, "filters", Any[]))
        if haskey(f, "transmission_csv")
            (haskey(f, "material") || haskey(f, "thickness_um")) &&
                throw(_KeyedError(key("filters[$i]"), "give material and thickness_um, or transmission_csv"))
            push!(tables, _curve(f["transmission_csv"], key("filters[$i].transmission_csv")))
        else
            (haskey(f, "material") && haskey(f, "thickness_um")) ||
                throw(_KeyedError(key("filters[$i]"), "give material and thickness_um, or transmission_csv"))
            m = Symbol(f["material"])
            t = 1e-3 * f["thickness_um"]
            push!(filters, haskey(f, "density_g_cm3") ? m => (t, Float64(f["density_g_cm3"])) : m => t)
            try
                filter_transmission([filters[end]], 1.0; table = :henke)
            catch e
                throw(_KeyedError(key("filters[$i].material"), sprint(showerror, e)))
            end
        end
    end
    mirror = nothing
    if haskey(c, "mirror")
        m = c["mirror"]
        if haskey(m, "reflectivity_csv")
            (haskey(m, "kind") || haskey(m, "angle_deg")) &&
                throw(_KeyedError(key("mirror"), "give kind and angle_deg, or reflectivity_csv"))
            mirror = _curve(m["reflectivity_csv"], key("mirror.reflectivity_csv"))
        else
            (haskey(m, "kind") && haskey(m, "angle_deg")) ||
                throw(_KeyedError(key("mirror"), "give kind and angle_deg, or reflectivity_csv"))
            k = Symbol(m["kind"])
            k in (:SiO2, :mlm) || throw(_KeyedError(key("mirror.kind"), "SiO2 or mlm (or reflectivity_csv)"))
            k === :mlm && m["angle_deg"] != 1.9 &&
                throw(_KeyedError(key("mirror.angle_deg"),
                                  "the mlm reflectivity is Rubery 2023 Fig. 8a at 1.9°; give reflectivity_csv for another"))
            mirror = (k, Float64(m["angle_deg"]))
        end
    end
    if haskey(c, "sensitivity_csv")
        haskey(c, "cathode") && throw(_KeyedError(key("cathode"), "give cathode or sensitivity_csv, not both"))
        cathode = _curve(c["sensitivity_csv"], key("sensitivity_csv"))
    else
        haskey(c, "cathode") || throw(_MissingValue(key("cathode"), "a cathode or sensitivity_csv is required"))
        cathode = Symbol(c["cathode"])
        haskey(_XRD_FIT, cathode) ||
            throw(_KeyedError(key("cathode"), "no XRD model for $cathode (Al, Cr, Ni); give sensitivity_csv"))
    end
    g = Float64(get(c, "grid", 1.0))
    0 < g <= 1 || throw(_KeyedError(key("grid"), "an open fraction in (0, 1]"))
    return DanteResponse(n, E0; filters, grid = g, mirror, cathode, solid_angle = Ω, tables)
end

# channels built once per content (the responses take ~0.1 s per channel set)
const _CHANNEL_CACHE = Dict{UInt,Vector{DanteResponse}}()
const _CHANNEL_CACHE_LOCK = ReentrantLock()

function _channels(cs::AbstractVector)
    h = hash(cs)
    hit = lock(() -> get(_CHANNEL_CACHE, h, nothing), _CHANNEL_CACHE_LOCK)
    hit === nothing || return hit
    numbers = [c["number"] for c in cs]
    allunique(numbers) || throw(_KeyedError("channel", "channel numbers must be unique"))
    ch = [_channel(c) for c in sort(cs; by = c -> c["number"])]
    lock(() -> (_CHANNEL_CACHE[h] = ch), _CHANNEL_CACHE_LOCK)
    return ch
end

function _build_xray_power(def::DiagnosticDefinition; kwargs...)
    d = def.data
    los = d["line_of_sight"]
    v = get(d, "view", Dict{String,Any}())
    mode = get(get(d, "recording", Dict{String,Any}()), "mode", "history")
    mode == "history" || throw(_KeyedError("recording.mode", "xray_power records a history"))
    cs = get(d, "channel", Any[])
    isempty(cs) && throw(_MissingValue("channel", "no channels: supply [[channel]] tables (extend \"$(def.id)\")"))
    kw = (; through = Symbol(get(v, "through", "leh_bottom")), pixels = Int(get(v, "pixels", 400)),
          channels = _channels(cs), area = haskey(v, "area_mm2") ? Float64(v["area_mm2"]) : nothing)
    return DanteView(los["theta_deg"], los["phi_deg"]; merge(kw, values(kwargs))...)
end

register_diagnostic_kind(:xray_power, _XRAY_POWER_SPEC, _build_xray_power)

# the built-in definitions alone (the public ones; user files cannot change these)
_builtin_diagnostic(id; kwargs...) = diagnostic(_builtin_definition(id); kwargs...)
_builtin_definition(id) = _definitions([normpath(_BUILTIN_DIAGNOSTICS)], DiagnosticIssue[])[id]
