# Small helpers shared by the controllers: error messages, text-field parsing,
# number display, the background-job runner, and the cancellation signal.

"""
    JobCancelled()

Thrown inside a background job (from an objective or between beam chunks) to stop it
at the next check. What the job finished before stays available (the best design an
optimization found so far; the clearance rows already computed).
"""
struct JobCancelled <: Exception end

# First line of an exception's message, for status lines.
_errmsg(err) = first(split(sprint(showerror, err), '\n'))

# Run a background computation, capturing its outcome as `(; value, err)`.
_try_job(job) = try
    (; value = job(), err = nothing)
catch err
    (; value = nothing, err)
end

# Run `job` on a worker thread (`spawn`) or inline.
_run_job(job, spawn::Bool) = spawn ? errormonitor(Threads.@spawn job()) : job()

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
