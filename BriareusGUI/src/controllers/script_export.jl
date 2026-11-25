# This file writes a session as a Julia script: the state of the design as a literal, and
# the calls that turn it into a target, beams, the TaLIS check and the drive.

# A value of a session state as Julia source that evaluates to the same value.
function _julia_literal(io::IO, x, indent::Int = 0)
    pad = "    "^(indent + 1)
    if x isa AbstractDict
        isempty(x) && return print(io, "Dict{String,Any}()")
        println(io, "Dict{String,Any}(")
        for k in sort!(collect(keys(x)); by = string)
            print(io, pad, repr(String(k)), " => ")
            _julia_literal(io, x[k], indent + 1)
            println(io, ",")
        end
        print(io, "    "^indent, ")")
    elseif x isa AbstractVector
        if any(v -> v isa AbstractDict || v isa AbstractVector, x)
            println(io, "Any[")
            for v in x
                print(io, pad)
                _julia_literal(io, v, indent + 1)
                println(io, ",")
            end
            print(io, "    "^indent, "]")
        else
            T = isempty(x) ? "Any" : all(v -> v isa AbstractString, x) ? "String" :
                all(v -> v isa Bool, x) ? "Bool" : all(v -> v isa AbstractFloat, x) ? "Float64" :
                all(v -> v isa Integer, x) ? "Int" : "Any"
            print(io, T, "[")
            for (i, v) in enumerate(x)
                i > 1 && print(io, ", ")
                _julia_literal(io, v, indent)
            end
            print(io, "]")
        end
    elseif x isa AbstractString
        print(io, repr(String(x)))
    elseif x isa Symbol
        print(io, repr(String(x)))
    elseif x isa Bool || x isa Integer || x isa AbstractFloat || x === nothing
        print(io, repr(x))
    else
        throw(ArgumentError("a session state holds a $(typeof(x)), which has no literal form"))
    end
    return io
end

"""
    session_script(sess::DesignSession; title = "") -> String

The design of the session as a Julia script, for sweeps and batch work outside the
window. The script holds the [`session_state`](@ref) as a literal, without the window
layout, and rebuilds it with [`restore!`](@ref). It then defines `target` (a
`Briareus.Target`), `beams` (the enabled beams, `Briareus.Beam`s), `talis` (the report
of `Briareus.talis_check` with the settings of the design, or `nothing` where the rules
do not apply), `scenario` and `evaluation` (the drive of the design, or `nothing` when
it cannot be computed). Every value after `target` and `beams` comes from a core
function that takes them, so a changed beam list or target goes through the same calls.

The script runs with `briareus run script.jl`, or in Julia with Briareus and BriareusGUI
in the environment. CAD files are named by path and are not embedded: a design with CAD
parts needs them where the session read them.
"""
function session_script(sess::DesignSession; title::AbstractString = "")
    state = session_state(sess)
    delete!(state, "window")
    io = IOBuffer()
    println(io, "# Briareus design script", isempty(title) ? "" : ": " * title)
    println(io, "#")
    println(io, "# Run with `briareus run <this file>`, or in Julia with Briareus and BriareusGUI in the")
    println(io, "# environment. The design is the `state` below, the same data as a saved session.")
    println(io, "# `target` and `beams` are core Briareus values, and every result after them comes")
    println(io, "# from a core function, so a sweep changes `beams` (for example with `remake(b; aim = …)`)")
    println(io, "# and calls the same functions again.")
    println(io)
    println(io, "using Briareus, BriareusGUI.Controllers")
    println(io)
    print(io, "state = ")
    _julia_literal(io, state)
    println(io)
    println(io)
    println(io, "session = DesignSession()")
    println(io, "session.spawn[] = false                     # jobs run in place")
    println(io, "restore!(session, state)")
    println(io, "await_target!(session.target)")
    println(io)
    println(io, "target = session.target.target[]            # Briareus.Target")
    println(io, "beams = active_beams(session.beams)         # the enabled beams, Vector{Briareus.Beam}")
    println(io)
    println(io, "# The TaLIS check with the settings of the design, where the rules apply.")
    println(io, "talis = :talis in session_tasks(session) && !isempty(session.talis.capture[]) && !isempty(beams) ?")
    println(io, "    talis_check(target, beams; talis_kwargs(session.talis, session.review)...) : nothing")
    println(io)
    println(io, "# The drive: the scenario of the design and its evaluation with the beams that have power.")
    println(io, "build_scenario!(session)")
    println(io, "scenario = session.drive.scenario[]")
    println(io, "evaluation = scenario === nothing ? nothing : evaluate(scenario, display_beams(session.beams))")
    return String(take!(io))
end

"""
    write_session_script(path, sess::DesignSession) -> path

Writes [`session_script`](@ref) to `path` (File ▸ Export Julia script).
"""
function write_session_script(path::AbstractString, sess::DesignSession)
    write(path, session_script(sess; title = basename(path)))
    sess.status[] = "wrote the design as a Julia script to $(basename(path))"
    return path
end
