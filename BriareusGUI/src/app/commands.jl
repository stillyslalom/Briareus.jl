# The headless commands of `briareus`: a design file in, a report or a packet out, with
# no window and no graphics context. Each command loads the design into a session as the
# window does and calls one core function of Briareus on it, so a script gets the same
# numbers from `Briareus.talis_check`, `Briareus.omega1_review` and `Briareus.roe_scan`.

"The results tabs that `briareus design.briar --tab <name>` opens."
const TABS = ("beams", "drive", "time", "diagnostics", "optimize", "clearance", "sweep", "talis")

const USAGE = """
usage:
  briareus [design.briar] [--tab <name>]      open the window (tabs: $(join(TABS, ", ")))
  briareus run script.jl [args...]            run a Julia script with Briareus and the controllers loaded
  briareus talis <design> --out report.csv [--json]
                                              the TaLIS checks as CSV (and JSON); exit code 2 when a check fails
  briareus review-1w <design> --out dir/ [--threshold 1.0] [--json]
                                              the 1ω review packet: HTML, CSV and one SVG slice per tight pair
  briareus roe <design> --out dir/ [--json]   the beam rules of engagement: HTML, CSV, Markdown and port map
  briareus --shortcut                         add Briareus to the desktop and the applications menu
<design> is a saved session (.briar, .briar.toml) or a VisRad workspace (.vrw).
"""

"""
    parse_command(args) -> (; command, files, out, threshold, json, tab, rest)

The arguments of the `briareus` command. `command` is `"run"`, `"talis"`, `"review-1w"`
or `"roe"`, or `""` for the window. `files` are the arguments that are not options,
`out`, `threshold` and `tab` the values of `--out`, `--threshold` and `--tab` (`nothing`
when absent), `json` whether `--json` was given, and `rest` every argument after the
script of `run`. An option without its value or a `--threshold` that is not a number
throws an `ArgumentError`.
"""
function parse_command(args::AbstractVector{<:AbstractString})
    command = !isempty(args) && args[1] in ("run", "talis", "review-1w", "roe") ? String(args[1]) : ""
    as = command == "" ? collect(String, args) : collect(String, args[2:end])
    if command == "run"
        return (; command, files = first(as, 1), out = nothing, threshold = nothing, json = false, tab = nothing,
                rest = as[2:end])
    end
    files, out, threshold, tab, json = String[], nothing, nothing, nothing, false
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
        elseif !startswith(a, "--")
            push!(files, a)
        end
        i += 1
    end
    return (; command, files, out, threshold, json, tab, rest = String[])
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
    review_1w_command(sess::DesignSession, out; threshold = 1.0, json = false, title = "1ω review")
        -> (; review, paths)

The `briareus review-1w` command on a session: `Briareus.omega1_review` with the settings
of the session, written into the folder `out` (`Briareus.write_omega1_review`).
"""
function review_1w_command(sess::DesignSession, out::AbstractString; threshold::Real = 1.0, json::Bool = false,
                           title::AbstractString = "1ω review")
    inp = _talis_inputs(sess)
    inp isa String && throw(ArgumentError(inp))
    review = omega1_review(inp.target, inp.beams; threshold, inp.kw...)
    return (; review, paths = write_omega1_review(out, review; json, title))
end

"""
    roe_command(sess::DesignSession, out; json = false, title = "Beam rules of engagement")
        -> (; report, criteria, notes, paths)

The `briareus roe` command on a session: the drive scenario of the session
([`build_scenario!`](@ref)), the criteria that `Briareus.roe_default_criteria` selects
from it and the foils of the design, `Briareus.roe_scan` over the beams with power, and
the packet written into the folder `out` (`Briareus.write_roe_packet`). `notes` gives
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
    report = roe_scan(s, beams, choice.criteria)
    paths = write_roe_packet(out, report, beams; json, title, choice.notes)
    return (; report, choice.criteria, choice.notes, paths)
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
    elseif c.command == "review-1w"
        r = review_1w_command(sess, c.out; threshold = something(c.threshold, 1.0), c.json, title = "1ω review: " * name)
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
