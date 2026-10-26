# Small helpers shared by the controllers: error messages, text-field parsing,
# number display, and the background-job runner. (Long jobs cancel through the core's
# own callbacks: `optimize_pointing` and `clearance` stop when theirs returns `true`.)

# First line of an exception's message, for status lines.
_errmsg(err) = first(split(sprint(showerror, err), '\n'))

# Run a background computation, capturing its outcome as `(; value, err)`.
_try_job(job) = try
    (; value = job(), err = nothing)
catch err
    (; value = nothing, err)
end

# Run `job` on a worker thread (`spawn`) or inline. Not specialized on `job`: inferring
# the inline branch for a job's closure inferred the whole job on the calling (GUI)
# thread before it was spawned (0.6–0.9 s at a first CAD import or part put back).
Base.@nospecializeinfer _run_job(@nospecialize(job), spawn::Bool) =
    spawn ? errormonitor(Threads.@spawn job()) : job()

# Hand a worker's outcome `out` to the GUI thread (`deliver`) as `apply(out)`, compiled
# here first: the first call of a commit otherwise compiles on the GUI thread (most of a
# second for a workspace import's). Once compiled, `precompile` only looks it up.
function _deliver_compiled(deliver, apply, out)
    precompile(apply, (typeof(out),))
    deliver(() -> apply(out))
end

# Where background computations run: inline, or (with `spawn[]`) on a worker whose
# result is applied through `deliver[]` on the GUI thread.
function _runner(spawn::Base.RefValue{Bool}, deliver::Base.RefValue{Any})
    return function (job, apply)
        d = deliver[]
        _run_job(() -> (out = _try_job(job); d(() -> apply(out))), spawn[])
        return nothing
    end
end

"""
    display_number(x) -> String

A number as text-field text: as printed when that is short, otherwise rounded to 6
significant digits (the stored value stays exact until the field is edited).
"""
function display_number(x::Real)
    x isa Integer && return string(x)
    isinteger(x) && abs(x) < 1e15 && return string(Int(x))
    s = string(x)
    return length(s) <= 8 ? s : @sprintf("%.6g", x)
end

# Parse a number text field entry (finite; positive / nonnegative when asked).
function _parse_number(str::AbstractString, what::AbstractString; positive::Bool = false,
                       nonnegative::Bool = false)
    v = tryparse(Float64, strip(str))
    (v === nothing || !isfinite(v)) && throw(ArgumentError("$what must be a number, got \"$str\""))
    positive && v <= 0 && throw(ArgumentError("$what must be positive, got $str"))
    nonnegative && v < 0 && throw(ArgumentError("$what must not be negative, got $str"))
    return v
end
_parse_number(x::Real, what::AbstractString; kwargs...) = _parse_number(string(x), what; kwargs...)

# Flux in TW/mm² as W/cm² text.
_wcm2(x::Real) = @sprintf("%.3g", to_Wcm2(x))

# Seconds as elapsed-time text.
_seconds(t::Real) = t < 90 ? @sprintf("%.0f s", t) : @sprintf("%.1f min", t / 60)

# Write an observable only when its value changes (Observables notify on every write).
_assign!(obs, v) = (isequal(obs[], v) || (obs[] = v); obs)
