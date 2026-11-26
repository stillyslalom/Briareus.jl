# The headless commands of `briareus`: a design file in, a report or a packet out, with
# no window and no graphics context. Each command loads the design into a session as the
# window does and calls one core function of Briareus on it, so a script gets the same
# numbers from `Briareus.talis_check`, `Briareus.omega1_review`, `Briareus.roe_scan` and
# `Briareus.design_diff`. The optimization runs as the Optimize page of the window runs
# it, with the settings saved in the design.

"The results tabs that `briareus design.briar --tab <name>` opens."
const TABS = ("beams", "drive", "time", "diagnostics", "optimize", "clearance", "sweep", "outages", "talis")

const USAGE = """
usage:
  briareus [design.briar] [--tab <name>]      open the window (tabs: $(join(TABS, ", ")))
  briareus run script.jl [args...]            run a Julia script with Briareus and the controllers loaded
  briareus talis <design> --out report.csv [--json]
                                              the TaLIS checks as CSV (and JSON); exit code 2 when a check fails
  briareus review-1w <design> --out dir/ [--threshold 1.0] [--json] [--images]
                                              the 1ω review packet: HTML, CSV and one SVG slice per tight pair;
                                              --images adds a rendered Daylight view of each pair
  briareus roe <design> --out dir/ [--json]   the beam rules of engagement: HTML, CSV, Markdown and port map
  briareus diff <reference> <design> --out dir/ [--json]
                                              what changed between two designs: parts, beams, TaLIS margins, physics
  briareus optimize <design> --out dir/ [--json]
                                              the optimization saved in the design: the optimized session, its
                                              pointing as a VisRad SST file and, for a workspace, a copy of it
  briareus --shortcut                         add Briareus to the desktop and the applications menu
<design> is a saved session (.briar, .briar.toml) or a VisRad workspace (.vrw).
"""

"""
    parse_command(args) -> (; command, files, out, threshold, json, images, tab, rest)

The arguments of the `briareus` command. `command` is `"run"`, `"talis"`, `"review-1w"`,
`"roe"`, `"diff"` or `"optimize"`, or `""` for the window. `files` are the arguments that are not options,
`out`, `threshold` and `tab` the values of `--out`, `--threshold` and `--tab` (`nothing`
when absent), `json` and `images` whether `--json` and `--images` were given, and `rest`
every argument after the
script of `run`. An option without its value or a `--threshold` that is not a number
throws an `ArgumentError`.
"""
function parse_command(args::AbstractVector{<:AbstractString})
    command = !isempty(args) && args[1] in ("run", "talis", "review-1w", "roe", "diff", "optimize") ? String(args[1]) : ""
    as = command == "" ? collect(String, args) : collect(String, args[2:end])
    if command == "run"
        return (; command, files = first(as, 1), out = nothing, threshold = nothing, json = false, images = false,
                tab = nothing, rest = as[2:end])
    end
    files, out, threshold, tab, json, images = String[], nothing, nothing, nothing, false, false
    i = 1
    while i <= length(as)
        a = as[i]
        if a in ("--out", "--threshold", "--tab")
            i < length(as) || throw(ArgumentError("$a needs a value"))
            v = as[i += 1]
            a == "--out" ? (out = v) : a == "--tab" ? (tab = v) :
            (threshold = something(tryparse(Float64, v), Some(nothing));
             threshold === nothing && throw(ArgumentError("--threshold needs a number of mm, got \"$v\"")))
        elseif a == "--json"
            json = true
        elseif a == "--images"
            images = true
        elseif !startswith(a, "--")
            push!(files, a)
        end
        i += 1
    end
    return (; command, files, out, threshold, json, images, tab, rest = String[])
end

"""
    load_design(path) -> DesignSession

The session of a design file, without a window: a saved session (`load_session!`) or,
for a `.vrw` file, a VisRad workspace (`import_visrad_workspace!`). Jobs of the session
run in place.
"""
function load_design(path::AbstractString)
    isfile(path) || throw(ArgumentError("no such file: $path"))
    sess = DesignSession()
    endswith(lowercase(path), ".vrw") ? import_visrad_workspace!(sess, path; spawn = false) : load_session!(sess, path)
    await_target!(sess.target)
    return sess
end

# The target, enabled beams and TaLIS keywords of a session, or a message when the checks
# do not apply to it.
function _talis_inputs(sess::DesignSession)
    :talis in session_tasks(sess) || return "the TaLIS rules apply to NIF beams, and this design is not on NIF"
    beams = active_beams(sess.beams)
    isempty(beams) && return "no beams are enabled"
    isempty(sess.talis.capture[]) && return "no capture component: choose the components the beams land on and save the design"
    return (; target = sess.target.target[], beams, kw = talis_kwargs(sess.talis, sess.review))
end

"""
    talis_command(sess::DesignSession, out; json = false) -> (; report, paths)

The `briareus talis` command on a session: `Briareus.talis_check` with the settings of
the session, written to `out` as CSV (`Briareus.write_talis_csv`), and with `json` also
as JSON beside it. A session that the rules do not apply to throws an `ArgumentError`.
"""
function talis_command(sess::DesignSession, out::AbstractString; json::Bool = false)
    inp = _talis_inputs(sess)
    inp isa String && throw(ArgumentError(inp))
    report = talis_check(inp.target, inp.beams; inp.kw...)
    mkpath(dirname(abspath(out)))
    paths = String[write_talis_csv(out, report)]
    json && push!(paths, write_talis_json(first(splitext(out)) * ".json", report))
    return (; report, paths)
end

"""
    review_images(sess::DesignSession, review; max_images = 24) -> Dict{Int,String}

A rendered Daylight view of each approach of a 1ω review (`Briareus.omega1_review` of the
session), as PNG data URIs keyed by the index in `review.items`, for the `images` of
`Briareus.write_omega1_review`. Each picture is the 3D view of the session with the row
under inspection ([`focus_talis_row!`](@ref)) and the camera of
[`daylight_view`](@ref), drawn offscreen. At most `max_images` approaches are drawn, the
tightest first.

The pictures need a graphics context of GLMakie, which cannot share a process with the
window: this function serves the headless command and scripts, and it throws an
`ArgumentError` while the window is open. The window grabs the same views from its 3D
view instead (File ▸ Export 1ω review report).
"""
function review_images(sess::DesignSession, review::Omega1Review; max_images::Integer = 24)
    _SHELL[] === nothing ||
        throw(ArgumentError("the Daylight pictures are drawn offscreen, which the open window does not allow"))
    images = Dict{Int,String}()
    isempty(review.items) && return images
    tl = sess.talis
    start_talis!(sess; spawn = false)
    rows = talis_rows(tl)
    sc = scene_canvas(sess)
    mktempdir() do dir
        for (i, it) in enumerate(review.items)
            length(images) >= max_images && break
            k = findfirst(x -> x.beam == it.row.beam && x.rule === it.row.rule && x.component === it.row.component, rows)
            k === nothing && continue
            focus_talis_row!(tl, k)
            talis_view!(sess, :daylight)
            png = joinpath(dir, "view.png")
            GLMakie.save(png, sc.fig; px_per_unit = 1)
            images[i] = "data:image/png;base64," * base64encode(read(png))
        end
    end
    clear_talis_focus!(tl)
    return images
end

"""
    review_1w_command(sess::DesignSession, out; threshold = 1.0, json = false, images = false,
                      title = "1ω review") -> (; review, paths)

The `briareus review-1w` command on a session: `Briareus.omega1_review` with the settings
of the session, written into the folder `out` (`Briareus.write_omega1_review`). With
`images` the page also shows a rendered Daylight view of each approach
([`review_images`](@ref)).
"""
function review_1w_command(sess::DesignSession, out::AbstractString; threshold::Real = 1.0, json::Bool = false,
                           images::Bool = false, title::AbstractString = "1ω review")
    inp = _talis_inputs(sess)
    inp isa String && throw(ArgumentError(inp))
    review = omega1_review(inp.target, inp.beams; threshold, inp.kw...)
    pics = images ? review_images(sess, review) : Dict{Int,String}()
    return (; review, paths = write_omega1_review(out, review; json, title, images = pics))
end

"""
    roe_command(sess::DesignSession, out; json = false, title = "Beam rules of engagement")
        -> (; report, criteria, notes, paths)

The `briareus roe` command on a session: the drive scenario of the session
([`build_scenario!`](@ref)), the criteria that `Briareus.roe_default_criteria` selects
from it and the foils of the design, `Briareus.roe_scan` over the beams with power (a NIF
quad held as one beam stands for its four beams, [`roe_stands_for`](@ref)), and the
packet written into the folder `out` (`Briareus.write_roe_packet`). `notes` gives
the criteria used and why. A design whose scenario cannot be built, or whose beams land
on nothing, throws an `ArgumentError`.
"""
function roe_command(sess::DesignSession, out::AbstractString; json::Bool = false,
                     title::AbstractString = "Beam rules of engagement")
    sess.spawn[] = false
    build_scenario!(sess)
    s = sess.drive.scenario[]
    s === nothing && throw(ArgumentError("the drive could not be computed: " * sess.drive.status[]))
    beams = display_beams(sess.beams)
    foils = [f for f in foil_names(sess.target) if f in parts(s.target)]
    choice = roe_default_criteria(s, beams; foils)
    isempty(choice.criteria) && throw(ArgumentError(join(choice.notes, " ")))
    report = roe_scan(s, beams, choice.criteria; stands_for = roe_stands_for(sess.beams))
    paths = write_roe_packet(out, report, beams; json, title, choice.notes)
    return (; report, choice.criteria, choice.notes, paths)
end

"""
    diff_command(a::DesignSession, b::DesignSession, out; names = ("reference", "current"), json = false,
                 title = "Design diff") -> (; diff, paths)

The `briareus diff` command on two sessions: `Briareus.design_diff` of the targets and
enabled beams, written into the folder `out` (`Briareus.write_design_diff`). The TaLIS
rows are compared when the rules apply to both designs, each checked with the settings
saved in it, and the physics when the drive of both can be computed. A comparison that
is left out is named in `notes` of the result.
"""
function diff_command(a::DesignSession, b::DesignSession, out::AbstractString; names = ("reference", "current"),
                      json::Bool = false, title::AbstractString = "Design diff")
    notes = String[]
    ia, ib = _talis_inputs(a), _talis_inputs(b)
    reports = nothing
    if ia isa String || ib isa String
        push!(notes, "TaLIS rows not compared: " * (ia isa String ? ia : ib))
    else
        reports = (talis_check(ia.target, ia.beams; ia.kw...), talis_check(ib.target, ib.beams; ib.kw...))
    end
    evs = map((a, b)) do s
        s.spawn[] = false
        build_scenario!(s)
        sc = s.drive.scenario[]
        sc === nothing ? nothing : evaluate(sc, display_beams(s.beams))
    end
    evaluations = any(isnothing, evs) ? nothing : evs
    evaluations === nothing && push!(notes, "physics not compared: the drive of one design could not be computed")
    foils = [f for f in foil_names(a.target) if f in foil_names(b.target)]
    diff = design_diff(a.target.target[], active_beams(a.beams), b.target.target[], active_beams(b.beams);
                       names, reports, evaluations, foils)
    return (; diff, notes, paths = write_design_diff(out, diff; json, title))
end

"""
    optimize_command(sess::DesignSession, out; json = false) -> (; result, paths)

The `briareus optimize` command on a session: the drive of the design
([`build_scenario!`](@ref)), the optimization of the enabled beams with the settings that
the design saves (the objective, the free fields and their bounds, the weights and the
keep-outs of the Optimize page, [`start_optimize!`](@ref)), and the optimum applied to
the beams ([`accept_optimum!`](@ref)). Written into the folder `out`: the optimized
session (`optimized.briar`), its pointing as a VisRad SST file (`pointing.csv`), a copy
of the VisRad workspace with the optimized beams (`optimized.vrw`) when the design was
imported from one, and with `json` a summary (`optimize.json`) of the objective before
and after and the beams with their aims and powers. A design whose drive cannot be
computed, or whose optimization cannot run or fails, throws an `ArgumentError` with the
reason.
"""
function optimize_command(sess::DesignSession, out::AbstractString; json::Bool = false)
    sess.spawn[] = false
    build_scenario!(sess)
    sess.drive.scenario[] === nothing && throw(ArgumentError("the drive could not be computed: " * sess.drive.status[]))
    st = sess.optimize
    st.result[] = nothing
    start_optimize!(sess; spawn = false)
    r = st.result[]
    r === nothing && throw(ArgumentError("the optimization did not run: " * (isempty(st.error[]) ? st.status[] : st.error[])))
    accept_optimum!(sess)
    mkpath(out)
    paths = String[save_session(joinpath(out, "optimized.briar"), sess), export_visrad_sst(sess.beams, joinpath(out, "pointing.csv"))]
    isempty(sess.visrad_source[]) || push!(paths, export_visrad_workspace(sess, joinpath(out, "optimized.vrw")))
    if json
        num(x) = isfinite(x) ? repr(Float64(x)) : "null"
        str(x) = "\"" * escape_string(string(x)) * "\""
        beam(b) = "    {\"name\": " * str(b.name) * ", \"aim_mm\": [" * join((num(v) for v in b.aim), ", ") *
                  "], \"power_TW\": " * num(b.power) * "}"
        text = "{\n  \"objective\": " * str(st.objective[]) * ",\n  \"initial_value\": " * num(r.initial_value) *
               ",\n  \"value\": " * num(r.value) * ",\n  \"evaluations\": " * string(r.evaluations) *
               ",\n  \"iterations\": " * string(r.iterations) * ",\n  \"reason\": " * str(r.reason) *
               ",\n  \"completed\": " * string(r.completed) * ",\n  \"status\": " * str(st.status[]) *
               ",\n  \"beams\": [\n" * join((beam(b) for b in r.beams), ",\n") * "\n  ]\n}\n"
        p = joinpath(out, "optimize.json")
        write(p, text)
        push!(paths, p)
    end
    return (; result = r, paths)
end

# Run a script as `briareus run` does: in a fresh module with Briareus and the
# controllers in scope, and its arguments in `ARGS`.
function run_script(path::AbstractString, args::AbstractVector{<:AbstractString} = String[])
    isfile(path) || throw(ArgumentError("no such file: $path"))
    m = Module(:BriareusScript)
    Core.eval(m, :(const Briareus = $Briareus; const BriareusGUI = $(@__MODULE__);
                   using .Briareus; using .BriareusGUI.Controllers;
                   include(p::AbstractString) = Base.include($m, p)))
    old = copy(Base.ARGS)
    empty!(Base.ARGS)
    append!(Base.ARGS, args)
    try
        Base.include(m, abspath(path))
    finally
        empty!(Base.ARGS)
        append!(Base.ARGS, old)
    end
    return m
end

# One headless command, with its messages on `io`. Returns the exit code.
function _run_command(c, io::IO = stdout)
    if c.command == "run"
        isempty(c.files) && (println(stderr, "briareus run: give a script"); return 1)
        run_script(c.files[1], c.rest)
        return 0
    end
    if c.command == "diff"
        (length(c.files) == 2 && c.out !== nothing) ||
            (println(stderr, "briareus diff: give two design files and --out\n", USAGE); return 1)
        names = (basename(c.files[1]), basename(c.files[2]))
        r = diff_command(load_design(c.files[1]), load_design(c.files[2]), c.out; names, c.json,
                         title = "Design diff: " * names[1] * " against " * names[2])
        foreach(x -> println(io, x), r.notes)
        println(io, r.diff)
        foreach(p -> println(io, "wrote ", p), r.paths)
        return 0
    end
    (length(c.files) == 1 && c.out !== nothing) ||
        (println(stderr, "briareus ", c.command, ": give one design file and --out\n", USAGE); return 1)
    sess = load_design(c.files[1])
    name = basename(c.files[1])
    if c.command == "talis"
        r = talis_command(sess, c.out; c.json)
        n(st) = count(x -> x.status === st, r.report.rows)
        println(io, "TaLIS: ", n(:fail), " fail, ", n(:review), " review, ", n(:pass), " pass")
        foreach(p -> println(io, "wrote ", p), r.paths)
        return n(:fail) > 0 ? 2 : 0
    elseif c.command == "optimize"
        r = optimize_command(sess, c.out; c.json)
        println(io, sess.optimize.status[])
        foreach(p -> println(io, "wrote ", p), r.paths)
    elseif c.command == "review-1w"
        r = review_1w_command(sess, c.out; threshold = something(c.threshold, 1.0), c.json, c.images, title = "1ω review: " * name)
        println(io, r.review)
        foreach(p -> println(io, "wrote ", p), first(r.paths, 2))
        length(r.paths) > 2 && println(io, "wrote ", length(r.paths) - 2, " more files in ", c.out)
    else
        r = roe_command(sess, c.out; c.json, title = "Beam rules of engagement: " * name)
        foreach(x -> println(io, x), r.notes)
        println(io, r.report)
        foreach(p -> println(io, "wrote ", p), r.paths)
    end
    return 0
end

"""
    compare_design(sess::DesignSession, reference::DesignSession) -> (; diff, ghost, arrows, text)

The design of `sess` against a `reference` design (a session loaded from a saved file or
a VisRad workspace) for the 3D view to overlay: `diff`, the core `DesignDiff` of the
parts, beams and, when the TaLIS rules apply to both designs, the TaLIS rows
(`Briareus.design_diff`, the reference first); `ghost`, the components of the reference
whose part moved, changed or is gone in `sess`, drawn where the reference had them;
`arrows`, one pair (reference aim, current aim) per beam whose aim moved; and `text`, a
line about the differences.
"""
function compare_design(sess::DesignSession, reference::DesignSession)
    ta, tb = reference.target.target[], sess.target.target[]
    ba, bb = active_beams(reference.beams), active_beams(sess.beams)
    ia, ib = _talis_inputs(reference), _talis_inputs(sess)
    reports = ia isa String || ib isa String ? nothing :
              (talis_check(ia.target, ia.beams; ia.kw...), talis_check(ib.target, ib.beams; ib.kw...))
    d = design_diff(ta, ba, tb, bb; names = ("reference", "current"), reports)
    differ = Set(p.part for p in d.parts if p.status in (:moved, :changed, :removed))
    ghost = Component[c for c in ta.components if c.part in differ && c.kind !== :aperture]
    aims = Dict(b.name => Briareus.Vec3(b.aim) for b in bb)
    arrows = Tuple{Briareus.Vec3,Briareus.Vec3}[]
    for x in d.beams
        (x.status === :changed && x.distance > 0 && haskey(aims, x.beam)) || continue
        b = aims[x.beam]
        push!(arrows, (b - x.aim_shift, b))
    end
    n(f) = count(f, d.parts)
    io = IOBuffer()
    @printf(io, "against the reference: %d parts moved, %d changed, %d removed, %d added; %d beams changed",
            n(p -> p.status === :moved), n(p -> p.status === :changed), n(p -> p.status === :removed),
            n(p -> p.status === :added), count(x -> x.status !== :same, d.beams))
    isempty(d.talis) || @printf(io, "; TaLIS: %d newly failing, %d fixed, %d worse, %d improved",
                                count(x -> x.change === :newly_failing, d.talis), count(x -> x.change === :fixed, d.talis),
                                count(x -> x.change === :worse, d.talis), count(x -> x.change === :improved, d.talis))
    return (; diff = d, ghost, arrows, text = String(take!(io)))
end
