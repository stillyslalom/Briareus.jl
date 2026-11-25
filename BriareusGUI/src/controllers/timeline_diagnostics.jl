# This file covers the diagnostics a drive history records, with their catalog, definitions
# and folders, the recorded signals and images, and the history export.

# ---------------------------------------------------------------- diagnostics

# Colors of the recorded diagnostics on the charts (Okabe–Ito, after those of Tr and
# the modes).
const _DIAGNOSTIC_COLORS = ("#D55E00", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
_diagnostic_color(j) = _DIAGNOSTIC_COLORS[mod1(j, length(_DIAGNOSTIC_COLORS))]

_builtin_file(f) = startswith(normpath(f), normpath(Briareus.GUIAPI.BUILTIN_DIAGNOSTICS))

"""
    refresh_diagnostics!(sess::DesignSession) -> Vector{NamedTuple}

Read the diagnostic definitions on the search path again (`diagnostic_paths`: the
built-ins, the saved folders, `BRIAREUS_DIAGNOSTICS`) with their problems
(`check_diagnostics`, in `sess.time.issues`), and list them with those saved with the
session in `sess.time.catalog`, one row `(; id, name, kind, origin, file, chain,
public, usable, problem, changed, facility, instrument, config)` per id: `origin` is
`:builtin`, `:user` (a folder on the search path) or `:session` (saved with the session
and not on the search path); `usable` whether it builds, else `problem` says why;
`changed` the keys in which a definition saved with the session differs from the search
path's of that id (`diagnostic_diff`; the session's is used); `facility` the definition's
(`""` without one); `instrument` the instrument it is a configuration of
(`diagnostic_instrument`) and `config` what tells it from the instrument's other
configurations (its name without the instrument's: "2012–14" for "NIF Dante 1
(2012–14)", "default" for the definition named as the instrument). The built-ins
first, then by id.
"""
function refresh_diagnostics!(sess::DesignSession)
    tm = sess.time
    issues = check_diagnostics()
    raws = Briareus.GUIAPI.read_definitions(diagnostic_paths(), DiagnosticIssue[])
    found = Briareus.GUIAPI.checked_definitions(raws, DiagnosticIssue[])
    problem(d) = join([i.message * (isempty(i.key) ? "" : " [" * i.key * "]") for i in issues
                       if i.file == d.file && i.level !== :warning], "; ")
    # Definitions that lack values (not public) but have no error are listed too.
    for id in setdiff(keys(raws), keys(found))
        any(i -> i.file == raws[id].file && i.level === :error, issues) && continue
        d = Briareus.GUIAPI.resolve_definition(id, raws, DiagnosticIssue[])
        d === nothing || (found[id] = d)
    end
    rows = NamedTuple[]
    for id in union(keys(found), keys(tm.embedded[]))
        e = get(tm.embedded[], id, nothing)
        p = get(found, id, nothing)
        d = e === nothing ? p : e.def
        chain = e === nothing ? d.chain : e.chain
        msg = e === nothing ? problem(d) : ""
        inst = diagnostic_instrument(d)
        push!(rows, (; id, d.name, d.kind, origin = p === nothing ? :session : _builtin_file(p.file) ? :builtin : :user,
                     file = e === nothing ? d.file : e.file, chain, d.public, usable = isempty(msg), problem = msg,
                     changed = e === nothing || p === nothing ? String[] : diagnostic_diff(e.def, p),
                     facility = String(get(d.data, "facility", "")), instrument = inst,
                     config = _config_label(d.name, inst)))
    end
    sort!(rows; by = r -> (r.origin !== :builtin, r.id))
    tm.issues[] = issues
    tm.catalog[] = rows
    return rows
end

# What tells a definition named `name` from the other configurations of `instrument`: the
# name without the instrument's and the brackets around the rest; "default" for the
# definition named as the instrument.
function _config_label(name::AbstractString, instrument::AbstractString)
    name == instrument && return "default"
    s = startswith(name, instrument) ? name[nextind(name, ncodeunits(instrument)):end] : name
    s = strip(s, [' ', ',', ':', '·', '-', '–'])
    startswith(s, '(') && endswith(s, ')') && (s = strip(s[nextind(s, 1):prevind(s, end)]))
    return isempty(s) ? "default" : String(s)
end

# A facility name as written in a definition or as the session's symbol, for comparison.
_facility_key(f) = replace(lowercase(String(f)), r"[^a-z0-9]" => "")

"""
    diagnostic_here(sess::DesignSession, facility) -> Bool

Whether a diagnostic of `facility` (a definition's `facility` key) belongs to the
session's: the same facility, a definition without one, or a session whose facility is
Other (`:custom`), which lists every definition.
"""
function diagnostic_here(sess::DesignSession, facility::AbstractString)
    f = sess.beams.facility[]
    return f === :custom || isempty(facility) || _facility_key(facility) == _facility_key(f)
end

"""
    diagnostic_instruments(sess::DesignSession) -> Vector{NamedTuple}

The diagnostics list of the window: one row per instrument of the session's facility
([`diagnostic_here`](@ref)), and per recorded instrument of another, grouped by kind as
[`refresh_diagnostics!`](@ref)'s catalog. Each row is `(; instrument, kind, group, head,
facility, here, usable, recorded, configs)`: `recorded` the ids of the configurations
the next drive history records, in the order recorded (empty: not recorded; one, as the
window records them, unless a script or an earlier session recorded several), `configs`
the instrument's definitions as catalog rows (the one named as the instrument first,
then by id), `usable` whether one of them builds, `head` whether the row is the first of
its `group`.
"""
function diagnostic_instruments(sess::DesignSession)
    rows = diagnostic_catalog(sess)
    ticked = sess.time.diagnostics[]
    by = Dict{Tuple{Symbol,String},Vector{Any}}()
    order = Tuple{Symbol,String}[]
    for g in _diagnostic_groups(rows)
        r = g.row
        k = (r.kind, r.instrument)
        haskey(by, k) || (by[k] = Any[]; push!(order, k))
        push!(by[k], r)
    end
    out = NamedTuple[]
    for k in order
        configs = sort(by[k]; by = r -> (r.config != "default", r.id))
        recorded = String[id for id in ticked if any(r -> r.id == id, configs)]
        here = any(r -> diagnostic_here(sess, r.facility), configs)
        (here || !isempty(recorded)) || continue
        group = _diagnostic_group(k[1])
        push!(out, (; instrument = k[2], kind = k[1], group, head = isempty(out) || out[end].group != group,
                    facility = configs[1].facility, here, usable = any(r -> r.usable, configs), recorded, configs))
    end
    return out
end

"""
    set_instrument!(sess::DesignSession, instrument, on::Bool = true) -> sess

Record `instrument` (a row of [`diagnostic_instruments`](@ref)) in the next drive
history, in its first configuration that builds (the default one first), or stop
recording it. [`set_diagnostic!`](@ref) chooses another configuration. An instrument
none of whose configurations builds is not recorded; the time status says why.
"""
function set_instrument!(sess::DesignSession, instrument::AbstractString, on::Bool = true)
    rows = [r for r in diagnostic_catalog(sess) if r.instrument == instrument]
    isempty(rows) && return sess
    configs = sort(rows; by = r -> (r.config != "default", r.id))
    if on
        any(r -> r.id in sess.time.diagnostics[], configs) && return sess
        j = findfirst(r -> r.usable, configs)
        j === nothing && (sess.status[] = "$instrument cannot be recorded: " * configs[1].problem; return sess)
        return set_diagnostic!(sess, configs[j].id, true)
    end
    foreach(r -> set_diagnostic!(sess, r.id, false), configs)
    return sess
end

"""
    set_configuration!(sess::DesignSession, id, new) -> sess

Record definition `new`, another configuration of the same instrument chosen from the
menu of a chip in the window, in place of the recorded definition `id`, where `id` was in
the list. When `new` is recorded already, `id` is dropped.
"""
function set_configuration!(sess::DesignSession, id::AbstractString, new::AbstractString)
    id, new = String(id), String(new)
    ids = copy(sess.time.diagnostics[])
    k = findfirst(==(id), ids)
    (k === nothing || id == new) && return sess
    new in ids ? deleteat!(ids, k) : (ids[k] = new)
    sess.time.diagnostics[] = ids
    w = _view_part_warning(sess)
    isempty(w) || (sess.status[] = w)
    return sess
end

"""
    diagnostic_catalog(sess::DesignSession) -> Vector{NamedTuple}

The rows of [`refresh_diagnostics!`](@ref), read from the search path the first time.
"""
diagnostic_catalog(sess::DesignSession) =
    isempty(sess.time.catalog[]) ? refresh_diagnostics!(sess) : sess.time.catalog[]

# The definition a session uses for `id`: the one saved with it, else the search path's
# (`nothing` without either).
function _definition(sess::DesignSession, id::AbstractString)
    e = get(sess.time.embedded[], id, nothing)
    e === nothing || return e.def
    r = findfirst(r -> r.id == id, diagnostic_catalog(sess))
    r === nothing && return nothing
    row = diagnostic_catalog(sess)[r]
    # (`nothing` too for one that does not build: it lacks values, e.g. ones that are not public)
    return row.origin === :session ? nothing : get(Briareus.GUIAPI.definitions(diagnostic_paths(), DiagnosticIssue[]), id, nothing)
end

"""
    set_diagnostic!(sess::DesignSession, id, on::Bool = true) -> sess

Record the diagnostic of definition `id` in the next drive history (`on`), or not.
Several definitions of one instrument can be recorded side by side (two fields of view
of one Dante); the window records an instrument once and swaps its configuration
([`set_instrument!`](@ref), [`set_configuration!`](@ref)). A ticked definition whose
view names a part the design lacks (`center_part`, `area_part`; e.g. a part renamed
before it was ticked, its definition being the user's file) is warned about in the
status bar, as it is when a drive history starts.
"""
function set_diagnostic!(sess::DesignSession, id::AbstractString, on::Bool = true)
    tm = sess.time
    ids = filter(!=(String(id)), tm.diagnostics[])
    on && push!(ids, String(id))
    ids == tm.diagnostics[] || (tm.diagnostics[] = ids)
    w = on ? _view_part_warning(sess) : ""
    isempty(w) || (sess.status[] = w)
    return sess
end

"""
    use_path_diagnostic!(sess::DesignSession, id) -> sess

Use the search path's definition of `id` instead of the one saved with the session
(when the two differ, see [`refresh_diagnostics!`](@ref)).
"""
function use_path_diagnostic!(sess::DesignSession, id::AbstractString)
    e = copy(sess.time.embedded[])
    pop!(e, String(id), nothing) === nothing && return sess
    sess.time.embedded[] = e
    refresh_diagnostics!(sess)
    return sess
end

# Definition `def` with the session's parts `vp` (`time.view_parts[id]`; `nothing`: none)
# in its view in place of the file's: a copy, the definition itself unchanged.
_with_view_parts(def, ::Nothing) = def
_with_view_parts(::Nothing, vp::AbstractDict) = nothing
function _with_view_parts(def, vp::AbstractDict)
    data = copy(def.data)
    data["view"] = merge(Dict{String,Any}(get(def.data, "view", Dict{String,Any}())), vp)
    return Briareus.DiagnosticDefinition(def.id, def.name, def.kind, def.file, def.chain, def.public, data)
end

# The ticked diagnostics whose view names a part the design lacks: `(id, key, name)` per
# such name (the session name, `time.view_parts`, where it has one).
function _missing_view_parts(sess::DesignSession)
    names = Symbol[component_names(sess.target.target[])...; [p.name for p in sess.target.parts[]]...]
    out = Tuple{String,String,String}[]
    for id in sess.time.diagnostics[]
        def = _definition(sess, id)
        def === nothing && continue
        for (k, name) in _view_parts(def, get(sess.time.view_parts[], id, nothing))
            _piece_named(names, Symbol(name)) in names || push!(out, (id, k, name))
        end
    end
    return out
end

# A status warning for the ticked diagnostics whose view names a part the design lacks,
# or an empty string when there are none.
function _view_part_warning(sess::DesignSession)
    miss = _missing_view_parts(sess)
    isempty(miss) && return ""
    return "the design has no part " * join(["$name ($id's $k)" for (id, k, name) in miss], ", ") *
           ": the diagnostic records nothing (rename the part back, or edit its definition)"
end

# The ticked definitions, in the order ticked: `(; id, name, def, note)`, `def` nothing
# (and `note` why) for one the session cannot build.
function _recorded_definitions(sess::DesignSession)
    rows = Dict(r.id => r for r in diagnostic_catalog(sess))
    out = NamedTuple[]
    for id in sess.time.diagnostics[]
        r = get(rows, id, nothing)
        if r === nothing
            push!(out, (; id, name = id, def = nothing, note = "not on the diagnostic search path"))
        elseif !r.usable
            push!(out, (; id, r.name, def = nothing, note = r.problem))
        else
            def = _with_view_parts(_definition(sess, id), get(sess.time.view_parts[], id, nothing))
            push!(out, (; id, r.name, def, note = ""))
        end
    end
    return out
end

# The instrument of each recorded definition on a target with components `names`
# (worker): `(; id, name, view, note)`, `view` nothing (and `note` why) when it records
# nothing in a history: a kind recorded neither as a power history nor as images, a
# view through an aperture the target lacks, or a field centered on a part it lacks.
function _diagnostic_views(diags, names)
    out = NamedTuple[]
    for d in diags
        view, note = nothing, d.note
        if d.def !== nothing
            if !(d.def.kind in (:xray_power, :xray_imager))
                note = "a $(d.def.kind) diagnostic is not recorded in a drive history"
            else
                v = _view_on(diagnostic(d.def), names)
                if v isa DanteView && !(v.through === nothing || v.through in names)
                    note = "the target has no $(v.through) to look through"
                elseif v.center isa Symbol && !(v.center in names)
                    note = "the target has no $(v.center) to center the field on"
                else
                    view = v
                end
            end
        end
        push!(out, (; d.id, d.name, view, note = Ref(note)))
    end
    return out
end

# One step of a recorded diagnostic: `(; Tr, intensity, L)`, the brightness temperature
# (eV) and intensity (TW/sr) of a power diagnostic, or for an imager the detected
# radiance per element (`detected_radiance`, Float32) and NaN; NaN, with the reason kept,
# where the view gives none.
function _record_diagnostic!(v, sol)
    v.view === nothing && return (; Tr = NaN, intensity = NaN, L = nothing)
    v.view isa XrayImager &&
        return (; Tr = NaN, intensity = NaN, L = Float32.(detected_radiance(v.view, sol)))
    try
        r = dante(v.view, sol)
        return (; Tr = Float64(r.Tr), intensity = Float64(r.intensity), L = nothing)
    catch err
        err isa ArgumentError || rethrow()
        isempty(v.note[]) && (v.note[] = sprint(showerror, err))
        return (; Tr = NaN, intensity = NaN, L = nothing)
    end
end

# Whether target `t` has a closing entrance-hole lip: an aperture whose albedo is a
# `GapClosure` or a `LipClosure`, the apertures whose open fraction the core images
# weight a line of sight by (`render`, `imager_image`, `imager_record`).
_closing_lips(t::Target) =
    any(c -> c.kind === :aperture && c.material.albedo isa Union{GapClosure,LipClosure}, t.components)

# What history `h` keeps of recorded diagnostic `j` (worker): `(; id, name, kind, Tr,
# intensity, note, view, L, transmission, frames)`. For an imager, `L` the detected
# radiance per element and step, `transmission` the open fraction per element and step
# where the target has closing lips (`nothing` without them), `intensity` the detected
# intensity toward the pinhole (TW/sr: the total of the image over the object plane, the
# lines of sight through a lip weighted by its open fraction) and `frames` its recording (`imager_record`: one integrated image, or a
# frame per gate, the lines of sight through a lip weighted by its open fraction).
function _recorded(v, j, h, t::Target)
    rows = [r.diag[j] for r in h.record]
    if v.view isa XrayImager
        L = reduce(hcat, [r.L for r in rows])
        transmission = first(h.record).open === nothing ? nothing : reduce(hcat, [r.open for r in h.record])
        P, Q = Briareus.GUIAPI.imager_view(v.view, t)      # the element and the lip element of each pixel
        a = (v.view.camera.width / v.view.camera.nx)^2
        seen(i, k) = transmission === nothing || Q[i] == 0 ? 1.0 : Float64(transmission[Q[i], k])
        intensity = [a * sum((P[i] > 0 ? Float64(L[P[i], k]) * seen(i, k) : 0.0 for i in eachindex(P)); init = 0.0)
                     for k in axes(L, 2)]
        frames = imager_record(v.view, t, h.time, Float64.(L);
                               transmission = transmission === nothing ? nothing : Float64.(transmission))
        return (; v.id, v.name, kind = :xray_imager, Tr = fill(NaN, length(rows)), intensity,
                note = v.note[], v.view, L, transmission, frames)
    end
    return (; v.id, v.name, kind = :xray_power, Tr = [r.Tr for r in rows], intensity = [r.intensity for r in rows],
            note = v.note[], v.view, L = nothing, transmission = nothing, frames = nothing)
end

"""
    diagnostic_images(sess::DesignSession) -> Vector{NamedTuple}

The images the latest drive history holds of its recorded imagers, one row `(; id,
label, frame)` per image: `frame` 0 for the image at the time cursor, else the
recording's frame (the integrated image, or each gate's). Empty without a history or
an imager.
"""
function diagnostic_images(sess::DesignSession)
    h = sess.time.history[]
    h === nothing && return NamedTuple[]
    out = NamedTuple[]
    for d in h.diagnostics
        d.frames === nothing && continue
        push!(out, (; d.id, label = d.name * " at the cursor", frame = 0))
        if isempty(d.frames.times)
            push!(out, (; d.id, label = d.name * ", time-integrated", frame = 1))
        else
            for (k, τ) in enumerate(d.frames.times)
                push!(out, (; d.id, label = @sprintf("%s, frame %d (%.2f ns)", d.name, k, τ), frame = k))
            end
        end
    end
    return out
end

"""
    show_diagnostic_image!(sess::DesignSession, id, frame::Integer = 0) -> sess

Show an image of recorded imager `id` of the latest drive history in the image view
(as Compute ▸ X-ray image does, [`set_view!`](@ref) `:image`): `frame` 0 the image at
the time cursor's step (detected radiance, W/cm²/sr; the last step without a cursor),
else frame `frame` of its recording (detected fluence per steradian, J/cm²/sr: the
integral over the shot or against a gate). Where the target has closing entrance-hole
lips, a line of sight through a lip is weighted by the open fraction of the lip at each
step (`Briareus.imager_image`, `Briareus.imager_record`). The view opens on the region with signal
(`window`, or the whole field when that covers most of it); Ctrl + Shift + click on the
image shows the whole field. The image stays until another is asked
for: a drive edit does not replace it.
"""
function show_diagnostic_image!(sess::DesignSession, id::AbstractString, frame::Integer = 0)
    im = diagnostic_image(sess, id, frame)
    im === nothing && return sess
    sess.drive.image_request[] = nothing                     # not the drive's camera image
    sess.drive.image[] = im
    set_view!(sess.drive, :image)
    return sess
end

"""
    diagnostic_image(sess::DesignSession, id, frame::Integer = 0) -> Union{Nothing,NamedTuple}

An image of recorded imager `id` of the latest drive history, as
[`show_diagnostic_image!`](@ref) shows it: `(; x, y, values, title, unit, window)` with
`values[i, j]` at `(x[i], y[j])` (mm on the object plane) and `window` the region with
signal. `nothing` without a history or such an imager.
"""
function diagnostic_image(sess::DesignSession, id::AbstractString, frame::Integer = 0)
    h = sess.time.history[]
    h === nothing && return nothing
    j = findfirst(d -> d.id == id && d.frames !== nothing, h.diagnostics)
    j === nothing && return nothing
    d = h.diagnostics[j]
    t = sess.drive.scenario[].target
    if frame == 0
        k = something(time_step(sess.time), length(h.time))
        transmission = d.transmission === nothing ? nothing : Float64.(view(d.transmission, :, k))
        img = to_Wcm2.(imager_image(d.view, t, Float64.(view(d.L, :, k)); transmission))
        title, unit = @sprintf("%s at %.2f ns", d.name, h.time[k]), "W/cm²/sr (detected)"
    else
        f = d.frames
        1 <= frame <= length(f.images) || throw(ArgumentError("no frame $frame"))
        img = 1e5 .* f.images[frame]                          # kJ/mm²/sr to J/cm²/sr
        title = isempty(f.times) ? d.name * ", time-integrated" :
                @sprintf("%s, frame %d (%.2f ns, gate %.2f ns)", d.name, frame, f.times[frame], d.view.gate)
        unit = "J/cm²/sr (detected)"
    end
    cam = d.view.camera
    w = cam.width / 2
    x = range(-w, w; length = cam.nx)
    values = permutedims(img[end:-1:1, :])
    return (; x, y = x, values, title, unit, window = _signal_window(x, x, values))
end

# An image as CSV: a comment line with `title`, then a header of the x coordinates (mm)
# and one row per y coordinate (mm) with the values; `values[i, j]` is at `(x[i], y[j])`.
function _write_image(path::AbstractString, x, y, values, title::AbstractString)
    open(path, "w") do io
        println(io, "# ", title, "; rows: y_mm, columns: x_mm")
        println(io, "y_mm\\x_mm,", join((@sprintf("%.6g", a) for a in x), ","))
        for j in eachindex(y)
            println(io, @sprintf("%.6g", y[j]), ",", join((@sprintf("%.6g", values[i, j]) for i in eachindex(x)), ","))
        end
    end
    return path
end

"""
    export_history(sess::DesignSession, dir) -> Vector{String}

Write the latest drive history to the folder `dir` as CSV files and return their paths:

- `history.csv`: against `time_ns`, the observed surface's `Tr_eV` and the capsule's
  `a1_over_a0` … `a4_over_a0` (or a flat face's mean, peak and lowest flux in W/cm² and
  its RMS nonuniformity), each cavity's wall temperature (`<cavity>_wall_Tr_eV`), each
  energy series (`energy_<name>_kJ`, cumulative) and each recorded diagnostic's
  brightness temperature and intensity toward it (`<id>_Tr_eV`, `<id>_TW_per_sr`; an
  imager's detected intensity);
- `history_laser_Wcm2.csv`, `history_radiosity_Wcm2.csv`, `history_incident_Wcm2.csv`:
  per element of the drive's target (`element`, `component`) one column per step
  (`t_<time>ns`), the fields the time cursor shows;
- `image_<id>_<frame>.csv`: each recorded imager's recording (`integrated`, or
  `frame<k>` per gate; detected fluence per steradian, J/cm²/sr), as [`export_results`](@ref)
  writes the image view's image.

Throws without a history; [`export_results`](@ref) calls it for a current one.
"""
function export_history(sess::DesignSession, dir::AbstractString)
    h = sess.time.history[]
    h === nothing && throw(ArgumentError("run the drive history first"))
    mkpath(dir)
    files = String[]
    cols = Pair{String,Vector{Float64}}["time_ns" => h.time, "Tr_eV" => h.Tr]
    for l in axes(h.legendre, 1)
        push!(cols, "a$(l)_over_a0" => h.legendre[l, :])
    end
    f = h.face
    f === nothing || append!(cols, ["face_mean_Wcm2" => f.mean, "face_peak_Wcm2" => f.peak,
                                    "face_lowest_Wcm2" => f.low, "face_rms" => f.rms])
    for c in h.cavities
        push!(cols, String(c.name) * "_wall_Tr_eV" => c.Tr)
    end
    for s in h.energy
        push!(cols, "energy_" * String(s.name) * "_kJ" => s.values)
    end
    for d in h.diagnostics
        any(isfinite, d.Tr) && push!(cols, d.id * "_Tr_eV" => d.Tr)
        any(isfinite, d.intensity) && push!(cols, d.id * "_TW_per_sr" => d.intensity)
    end
    path = joinpath(dir, "history.csv")
    open(path, "w") do io
        println(io, join(first.(cols), ","))
        for k in eachindex(h.time)
            println(io, join((@sprintf("%.8g", last(c)[k]) for c in cols), ","))
        end
    end
    push!(files, path)
    # the fields per element and step, on the drive's elements
    s = sess.drive.scenario[]
    t = s === nothing ? nothing : s.target
    for field in DRIVE_FIELDS
        M = get(h.fields, field, nothing)
        (M === nothing || t === nothing || size(M, 1) != nelements(t)) && continue
        path = joinpath(dir, "history_$(field)_Wcm2.csv")
        open(path, "w") do io
            println(io, "element,component,", join((@sprintf("t_%.4gns", τ) for τ in h.time), ","))
            for e in axes(M, 1)
                println(io, e, ",", t.components[t.elem_comp[e]].name, ",",
                        join((@sprintf("%.6g", M[e, k]) for k in axes(M, 2)), ","))
            end
        end
        push!(files, path)
    end
    # the imagers' recordings
    for d in h.diagnostics
        d.frames === nothing && continue
        cam = d.view.camera
        x = range(-cam.width / 2, cam.width / 2; length = cam.nx)
        for (k, img) in enumerate(d.frames.images)
            name = isempty(d.frames.times) ? "integrated" : "frame$k"
            title = isempty(d.frames.times) ? d.name * ", time-integrated" :
                    @sprintf("%s, frame %d (%.2f ns)", d.name, k, d.frames.times[k])
            values = permutedims(1e5 .* img[end:-1:1, :])               # kJ/mm²/sr to J/cm²/sr
            push!(files, _write_image(joinpath(dir, "image_$(d.id)_$name.csv"), x, x, values,
                                      title * ", J/cm²/sr (detected)"))
        end
    end
    return files
end

# The region of an image with signal, for the image view's first limits: `(x0, x1, y0,
# y1)` (mm) bounding the pixels above `frac` of the peak, widened by `margin` of its
# extent on each side (at least 4 pixels) and clamped to the field; `nothing` when that
# region covers most of the field or the image is empty. `values[i, j]` is at `(x[i],
# y[j])`. The threshold and margin are display defaults; no instrument criterion sets them.
function _signal_window(x, y, values; frac = 1e-3, margin = 0.15)
    fin = filter(isfinite, values)
    peak = isempty(fin) ? 0.0 : maximum(fin)
    peak > 0 || return nothing
    hot = findall(v -> isfinite(v) && v > frac * peak, values)
    is, js = extrema(i[1] for i in hot), extrema(i[2] for i in hot)
    dx, dy = step(x), step(y)
    function span(r, lo, hi, d)
        pad = max(margin * (r[hi] - r[lo]), 4 * abs(d))
        return max(first(r) - abs(d) / 2, r[lo] - pad), min(last(r) + abs(d) / 2, r[hi] + pad)
    end
    x0, x1 = span(x, is..., dx)
    y0, y1 = span(y, js..., dy)
    (x1 - x0) * (y1 - y0) > 0.5 * (last(x) - first(x)) * (last(y) - first(y)) && return nothing
    return (x0, x1, y0, y1)
end

# The group of a diagnostic kind in the time panel's definition picker.
const _DIAGNOSTIC_GROUPS = (:xray_power => "X-ray power", :xray_imager => "X-ray imagers",
                            :radiography => "Radiography")
_diagnostic_group(kind::Symbol) =
    something(get(Dict(_DIAGNOSTIC_GROUPS), kind, nothing), uppercasefirst(replace(String(kind), '_' => ' ')))

# The catalog's rows grouped by kind for the picker (the groups of `_DIAGNOSTIC_GROUPS`
# in that order, then the other kinds by name; the catalog's order within a group):
# `(; row, group, head)`, `head` for the first row of its group.
function _diagnostic_groups(rows)
    rank(k) = something(findfirst(p -> p.first === k, _DIAGNOSTIC_GROUPS), length(_DIAGNOSTIC_GROUPS) + 1)
    order = sort(collect(eachindex(rows)); by = i -> (rank(rows[i].kind), _diagnostic_group(rows[i].kind)))
    out = NamedTuple[]
    for i in order
        g = _diagnostic_group(rows[i].kind)
        push!(out, (; row = rows[i], group = g, head = isempty(out) || out[end].group != g))
    end
    return out
end

"""
    diagnostic_folders() -> Vector{NamedTuple}

The folders on the diagnostic search path besides the built-ins, `(; path, saved)`:
`saved` is true for those kept in the configuration file (Diagnostics ▸ Folders) and
false for those of `BRIAREUS_DIAGNOSTICS` or added for this Julia session only.
"""
function diagnostic_folders()
    saved = Set(normpath.(Briareus.GUIAPI.saved_paths()))
    return [(; path, saved = normpath(path) in saved) for path in diagnostic_paths(; builtin = false)]
end

"""
    add_diagnostic_folder!(sess::DesignSession, dir) -> Bool

Add the folder `dir` to the diagnostic search path and keep it for later sessions
(`add_diagnostic_path`), then read the definitions again. Returns whether it was added
(the status says why not).
"""
function add_diagnostic_folder!(sess::DesignSession, dir::AbstractString)
    isdir(dir) || (sess.time.status[] = "no such folder: $dir"; return false)
    add_diagnostic_path(dir; save = true)
    rows = refresh_diagnostics!(sess)
    n = count(r -> r.origin === :user && startswith(normpath(r.file), normpath(abspath(dir))), rows)
    sess.time.status[] = "diagnostic folder added: $n definition$(n == 1 ? "" : "s") in $(basename(normpath(dir)))"
    return true
end

"""
    remove_diagnostic_folder!(sess::DesignSession, dir) -> sess

Take the folder `dir` off the saved diagnostic search path and read the definitions
again.
"""
function remove_diagnostic_folder!(sess::DesignSession, dir::AbstractString)
    remove_diagnostic_path(dir; save = true)
    refresh_diagnostics!(sess)
    return sess
end
