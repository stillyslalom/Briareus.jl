# Placement sweeps: step a design through a sequence of placements (a foil along its line
# of sight, a part along a frame axis) and check each step for beam clearance, the line
# of sight and beam capture, optionally re-optimizing the pointing at each step.

"""
    SweepRow

One step of a [`placement_sweep`](@ref):

- `value` — the swept quantity (the convenience method: the shift along the axis, mm);
- `beams` — the step's beams (after re-optimization, when on);
- `min_clearance`, `worst_beam`, `worst_component`, `fails`, `reviews` — the least
  worst-case clearance of the [`clearance`](@ref) check (mm), its beam and hardware part,
  and the numbers of `:fail` and `:review` pairs (`NaN`, `""`, `:none`, 0 when not
  checked);
- `los_clear`, `blocked_fraction`, `covered` — the [`line_of_sight`](@ref) check: clear,
  the fraction of its rays obstructed, the fraction of the object covered (`missing`,
  `NaN` when not checked);
- `min_captured` — the least fraction of a beam's power absorbed by the capture parts
  (`NaN` when not checked);
- `objective_before`, `objective_after`, `x` — the objective at the step's starting
  parameters and after re-optimization, and the optimized parameters (`NaN`, empty when
  not optimized);
- `clearance`, `los` — the full reports (`nothing` when not checked).
"""
struct SweepRow
    value::Float64
    beams::Vector{Beam{Float64}}
    min_clearance::Float64
    worst_beam::String
    worst_component::Symbol
    fails::Int
    reviews::Int
    los_clear::Union{Missing,Bool}
    blocked_fraction::Float64
    covered::Float64
    min_captured::Float64
    objective_before::Float64
    objective_after::Float64
    x::Vector{Float64}
    clearance::Union{Nothing,ClearanceReport}
    los::Union{Nothing,LineOfSightReport}
end

"""
    SweepResult

Result of [`placement_sweep`](@ref): `values` (the requested steps), `rows` (one
[`SweepRow`](@ref) per completed step, in order) and `canceled` (a callback stopped the
sweep; `rows` then hold the steps completed before it). [`write_sweep_csv`](@ref) writes
the rows.
"""
struct SweepResult
    values::Vector{Float64}
    rows::Vector{SweepRow}
    canceled::Bool
end

"""
    placement_sweep(design, values; clearance_kw = nothing, los = nothing, los_kw = (;),
                    capture = nothing, capture_kw = (;), optimize = nothing,
                    callback = nothing) -> SweepResult
    placement_sweep(target, beams, parts, values; frame = CHAMBER, axis = :z,
                    follow = false, kwargs...) -> SweepResult

Steps a design through `values` and checks each step. `design(value)` returns a named
tuple with the step's `target` and `beams`, and with `pointing`, the parameterization
whose initial parameters give those beams, when the step is re-optimized.

Each check runs when its keyword is given:
- `clearance_kw`, the keywords of [`clearance`](@ref) (`targets` among them): the least
  worst-case clearance, its beam and part, and the numbers of failing and review pairs;
- `los`, a line of sight (θ, φ in degrees or a direction) for [`line_of_sight`](@ref)
  with keywords `los_kw`: whether it is clear, the fraction of its rays blocked and the
  object's coverage;
- `capture`, a part name or a tuple of them (a backlighter foil): the least fraction of
  a beam's power those parts absorb, from [`beam_component_power`](@ref) on a
  [`DepositionOperator`](@ref) built for the step with keywords `capture_kw`.

`optimize`, a named tuple, re-optimizes the pointing at each step before the checks:
`objective` (required), `scenario_kw` (keywords of [`Scenario`](@ref)), `bounds_kw`
(keywords of `bounds`, a box about the step's own initial parameters), `warm` (default
`true`: start from the previous step's optimum when it has the same number of
parameters) and any keywords of [`optimize_pointing`](@ref). With a
[`FramePointing`](@ref) in a frame that moves with the part, the parameters keep their
meaning from step to step, which is what makes the warm start useful.

`callback((; done, total, row))` is called once with `done = 0` and `row = nothing`
before any work, then after each step with that step's row; returning `true` cancels,
and the result holds the completed rows with `canceled = true` (a `true` after the last
step is ignored).

**Moving parts along a frame axis.** The second form moves the named `parts` (a
`Symbol` or a collection) of `target` by `value · a` at each step, `a` the unit `axis`
of `frame` in the chamber (`:x`, `:y`, `:z`, or a local vector). Beams selected by
`follow` (`true` for all, `false` for none, beam indices, or a predicate `beam -> Bool`)
move with them: their aims shift by the same vector. `beams` may be a
[`FramePointing`](@ref) instead: its template beams are swept, followed beams' frames
move with them, and it is the `pointing` that `optimize` re-optimizes. `value` is then the
shift (mm) and every other keyword passes to the first form.

Each step costs a target rebuild, the clearance check and a deposition operator for the
capture check, and the optimization when on: give these checks their quick settings
(`clearance_kw = (; nangles = 8, nspot = 31, ...)`, a coarse `capture_kw`) for long sweeps.
"""
function placement_sweep(design, values::AbstractVector{<:Real}; clearance_kw = nothing,
                         los = nothing, los_kw = (;), capture = nothing, capture_kw = (;),
                         optimize = nothing, callback = nothing)
    opts = (; clearance_kw, los, los_kw, capture, capture_kw, optimize)
    return _sweep(design, collect(Float64, values), opts, callback)
end

function placement_sweep(t::Target, beams, partnames, values::AbstractVector{<:Real};
                         frame::Frame = CHAMBER, axis = :z, follow = false, kwargs...)
    pnames = partnames isa Symbol ? Symbol[partnames] : collect(Symbol, partnames)
    for p in pnames
        part_pieces(t, p)                 # a KeyError for an unknown part
    end
    a = axis === :x ? Vec3(1, 0, 0) : axis === :y ? Vec3(0, 1, 0) :
        axis === :z ? Vec3(0, 0, 1) : normalize(Vec3(axis...))
    dir = Vec3(to_chamber(frame, a; direction = true))
    pointing = beams isa FramePointing ? beams : nothing
    bs = pointing === nothing ? collect(beams) : pointing.template
    moved = _followed(bs, follow)
    design = _MovedParts(t, Set(pnames), bs, moved, pointing, dir)
    return placement_sweep(design, values; kwargs...)
end

_followed(bs, f::Bool) = fill(f, length(bs))
function _followed(bs, f::AbstractVector{<:Integer})
    m = fill(false, length(bs))
    m[f] .= true
    return m
end
_followed(bs, f) = Bool[f(b) for b in bs]

# The design of the convenience method: parts moved by value · dir, followed beams with them.
struct _MovedParts{P}
    target::Target
    parts::Set{Symbol}
    beams::Vector
    moved::Vector{Bool}
    pointing::P
    dir::Vec3
end

function (d::_MovedParts)(value)
    Δ = value * d.dir
    t = d.target
    comps = Component[c.part in d.parts ? translate(c, Δ) : c for c in t.components]
    target = Target(comps; closed = t.closed, visibility = t.visibility)
    beams = Beam{Float64}[]
    for (b, m) in zip(d.beams, d.moved)
        bv = _value_beam(b)
        push!(beams, m ? remake(bv; aim = bv.aim + Δ) : bv)
    end
    pointing = d.pointing === nothing ? nothing : _translated(d.pointing, Δ, d.moved)
    return (; target, beams, pointing)
end

# The sweep, compiled once for every design, objective and callback.
Base.@nospecializeinfer function _sweep(@nospecialize(design), values::Vector{Float64},
                                        @nospecialize(opts::NamedTuple), @nospecialize(callback))
    rows = SweepRow[]
    nv = length(values)
    finish(canceled) = SweepResult(values, rows, canceled)
    callback === nothing || callback((; done = 0, total = nv, row = nothing)) === true &&
        nv > 0 && return finish(true)
    opt = opts.optimize
    if opt !== nothing
        haskey(opt, :objective) || throw(ArgumentError("optimize needs an `objective`"))
    end
    xprev = nothing
    for (k, v) in enumerate(values)
        d = design(v)
        t::Target = d.target
        beams = Beam{Float64}[_value_beam(b) for b in d.beams]
        x, J0, J1 = Float64[], NaN, NaN
        if opt !== nothing
            p = get(d, :pointing, nothing)
            p === nothing && throw(ArgumentError("re-optimizing needs the design's `pointing`"))
            s = Scenario(t, p; get(opt, :scenario_kw, (;))...)
            xi = initial_parameters(p)
            lo, hi = bounds(p, xi; get(opt, :bounds_kw, (;))...)
            warm = get(opt, :warm, true)
            x0 = warm && xprev !== nothing && length(xprev) == length(xi) ? xprev : xi
            rest = Base.structdiff(opt, NamedTuple{(:objective, :scenario_kw, :bounds_kw, :warm)})
            res = optimize_pointing(s, opt.objective; x0, lower = lo, upper = hi, rest...)
            beams = Beam{Float64}[_value_beam(b) for b in res.beams]
            x, J0, J1 = res.x, res.initial_value, res.value
            xprev = x
        end
        row = _sweep_row(t, v, beams, x, J0, J1, opts)
        push!(rows, row)
        callback === nothing ||
            callback((; done = k, total = nv, row)) === true && k < nv && return finish(true)
    end
    return finish(false)
end

function _sweep_row(t::Target, v, beams, x, J0, J1, opts)
    cr = nothing
    cmin, wb, wc, nf, nr = NaN, "", :none, 0, 0
    if opts.clearance_kw !== nothing
        cr = clearance(t, beams; opts.clearance_kw...)
        if !isempty(cr.rows)
            w = argmin(r -> r.clearance, cr.rows)
            cmin, wb, wc = w.clearance, w.beam, w.component
        else
            cmin = Inf
        end
        nf = count(r -> r.status === :fail, cr.rows)
        nr = count(r -> r.status === :review, cr.rows)
    end
    lr = nothing
    lc, bf, cov = missing, NaN, NaN
    if opts.los !== nothing
        lr = line_of_sight(t, opts.los; opts.los_kw...)
        lc, bf, cov = lr.clear, lr.blocked / max(lr.rays, 1), lr.covered
    end
    cap = NaN
    if opts.capture !== nothing && !isempty(beams)
        names = opts.capture isa Symbol ? (opts.capture,) : Tuple(opts.capture)
        cols = unique!(Int[ci for n in names for ci in part_pieces(t, n)])
        op = DepositionOperator(t, beams; opts.capture_kw...)
        M = beam_component_power(op, beams)
        fr = [sum(M[i, cols]) / beams[i].power for i in eachindex(beams) if beams[i].power > 0]
        cap = isempty(fr) ? NaN : minimum(fr)
    end
    return SweepRow(v, beams, cmin, wb, wc, nf, nr, lc, bf, cov, cap, J0, J1, x, cr, lr)
end

function Base.show(io::IO, r::SweepResult)
    print(io, "SweepResult(", length(r.rows), " of ", length(r.values), " steps",
          r.canceled ? ", canceled" : "", ")")
end

function Base.show(io::IO, ::MIME"text/plain", r::SweepResult)
    @printf(io, "SweepResult: %d of %d steps%s\n", length(r.rows), length(r.values),
            r.canceled ? " (canceled: partial)" : "")
    print(io, "     value  clearance mm  worst beam / part        fail review  LOS      blocked  covered  captured  objective")
    f(x, fmt) = isnan(x) ? lpad("–", length(fmt(0.0))) : fmt(x)
    for w in r.rows
        who = w.worst_beam == "" ? "–" : string(w.worst_beam, " / ", w.worst_component)
        los = w.los_clear === missing ? "–" : w.los_clear ? "clear" : "blocked"
        J = isnan(w.objective_after) ? "–" :
            @sprintf("%.4g → %.4g", w.objective_before, w.objective_after)
        @printf(io, "\n  %+8.4g  %12s  %-22s %5d %6d  %-7s  %7s  %7s  %8s  %s", w.value,
                f(w.min_clearance, x -> @sprintf("%+12.4f", x)), who, w.fails, w.reviews, los,
                f(w.blocked_fraction, x -> @sprintf("%7.4f", x)),
                f(w.covered, x -> @sprintf("%7.4f", x)),
                f(w.min_captured, x -> @sprintf("%8.5f", x)), J)
    end
end

"""
    write_sweep_csv(path, result) -> path

Write the rows of a [`SweepResult`](@ref) as CSV, one line per step: the value, the least
clearance (mm), its beam and part, the failing and review pairs, the line of sight
(clear, blocked fraction, covered fraction), the least captured fraction and the
objective before and after re-optimization. A check that did not run leaves its cells
empty.
"""
function write_sweep_csv(path::AbstractString, r::SweepResult)
    g(x) = isnan(x) ? "" : @sprintf("%.8g", x)
    open(path, "w") do io
        println(io, "value,min_clearance_mm,worst_beam,worst_component,fails,reviews,los_clear,",
                "blocked_fraction,covered,min_captured,objective_before,objective_after")
        for w in r.rows
            checked = w.clearance !== nothing
            println(io, join((g(w.value), g(w.min_clearance), w.worst_beam,
                              checked ? string(w.worst_component) : "",
                              checked ? string(w.fails) : "", checked ? string(w.reviews) : "",
                              w.los_clear === missing ? "" : string(w.los_clear),
                              g(w.blocked_fraction), g(w.covered), g(w.min_captured),
                              g(w.objective_before), g(w.objective_after)), ","))
        end
    end
    return path
end
