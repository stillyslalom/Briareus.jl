# Beam clearance: geometric beam envelopes, hardware clearance distances, and
# per-beam/per-component intercepted power.

"""
    enclosed_radius(order, fraction)

Radius (in units of the 1/e semi-axes) of the super-Gaussian contour enclosing
`fraction` of the power: the fraction inside `ρ` is `P(2/n, ρⁿ)` (regularized lower
incomplete gamma function).
"""
enclosed_radius(n::Real, f::Real) = gamma_inc_inv(2 / n, f, 1 - f)^(1 / n)

"""
    BeamEnvelope(beam; enclosed = 0.99, margin = 0.0)

Geometric envelope of a beam: at axial distance `s` from best focus, the ellipse that
encloses `enclosed` of the focal-spot power, grown by the f-cone radius `|s|/(2F)` and
a uniform `margin` (mm; e.g. three times the rms pointing error). Clearance is the
signed distance to this envelope; negative means the hardware is inside it.
"""
struct BeamEnvelope
    focus::Vec3
    aim::Vec3
    k::Vec3
    e1::Vec3
    e2::Vec3
    a::Float64
    b::Float64
    cone::Float64          # tan(half-angle)
    margin::Float64
end

function BeamEnvelope(b::Beam; enclosed = 0.99, margin = 0.0)
    ρ = enclosed_radius(b.spot.order, enclosed)
    e1, e2 = spot_axes(b)
    k = direction(b)
    aim = Vec3(ForwardDiff.value.(b.aim))
    focus = aim + ForwardDiff.value(b.defocus) * k
    return BeamEnvelope(focus, aim, k, e1, e2, ρ * ForwardDiff.value(b.spot.a),
                        ρ * ForwardDiff.value(b.spot.b), 1 / (2b.fnumber), Float64(margin))
end

"""
    envelope_distance(env, x) -> (distance, s)

Signed distance (mm) from point `x` to the beam envelope in the plane through `x`
normal to the beam, and the axial position `s` of that plane relative to best focus.
"""
function envelope_distance(env::BeamEnvelope, x::Vec3)
    d = x - env.focus
    s = dot(d, env.k)
    u, v = dot(d, env.e1), dot(d, env.e2)
    return _ellipse_distance(u, v, env.a, env.b) - env.cone * abs(s) - env.margin, s
end

# Signed distance from (u, v) to the ellipse with semi-axes (a, b): negative inside.
# Eberly's robust point–ellipse distance (bisection on the Lagrange multiplier), valid
# inside and outside.
function _ellipse_distance(u, v, a, b)
    y0, y1 = abs(u), abs(v)
    e0, e1 = a, b
    if e0 < e1
        e0, e1, y0, y1 = e1, e0, y1, y0
    end
    inside = (y0 / e0)^2 + (y1 / e1)^2 < 1
    if y1 > 0
        if y0 > 0
            z0, z1 = y0 / e0, y1 / e1
            g = z0^2 + z1^2 - 1
            g == 0 && return 0.0
            r0 = (e0 / e1)^2
            n0 = r0 * z0
            s0, s1 = z1 - 1, g < 0 ? 0.0 : hypot(n0, z1) - 1
            s = s0
            for _ in 1:200
                s = (s0 + s1) / 2
                (s == s0 || s == s1) && break
                gs = (n0 / (s + r0))^2 + (z1 / (s + 1))^2 - 1
                gs > 0 ? (s0 = s) : gs < 0 ? (s1 = s) : break
            end
            x0, x1 = r0 * y0 / (s + r0), y1 / (s + 1)
            d = hypot(x0 - y0, x1 - y1)
        else
            d = abs(y1 - e1)
        end
    else
        numer, denom = e0 * y0, e0^2 - e1^2
        if numer < denom
            ξ = numer / denom
            d = hypot(e0 * ξ - y0, e1 * sqrt(1 - ξ^2))
        else
            d = abs(y0 - e0)
        end
    end
    return inside ? -d : d
end

"""
    ClearanceRow

One beam/component pair: minimum envelope `clearance` (mm, negative = overlap), the
hardware `point` where it occurs and its axial position `s` from best focus, and the
fraction of the beam's power the component `intercepts` (from the deposition
operator, both sides).
"""
struct ClearanceRow
    beam::String
    component::Symbol
    clearance::Float64
    point::Vec3
    s::Float64
    intercepted::Float64
end

"""
    ClearanceReport

Result of `clearance`: `rows` (one per beam and hardware component) plus the settings
used. `violations(report, required)` lists rows closer than `required` mm.
"""
struct ClearanceReport
    rows::Vector{ClearanceRow}
    enclosed::Float64
    margin::Float64
    targets::Vector{Symbol}
end

"""
    clearance(target, beams; targets, hardware = all others, enclosed = 0.99,
              margin = 0.0, level = 2, op = nothing) -> ClearanceReport

Beam-to-hardware clearance check. Each beam's `BeamEnvelope` is compared with every
`hardware` component, sampled at triangle vertices and `4^level` sub-triangle points.
Only hardware upstream of the beam's landing is considered: the landing is where the
beam axis first meets one of the `targets` components (the surfaces the beam is meant
to hit, e.g. `(:wall,)` in a hohlraum or the backlighter foil), or the aim plane if the
axis misses them.

The intercepted fraction per component includes the far wings of the spot that a
contour-based envelope leaves out. It comes from the reference ray trace (exact power
accounting at any hardware mesh resolution), or from `op` — a `DepositionOperator` for
the same beams — when given (faster; accurate when hardware meshes resolve the spot).
"""
function clearance(t::Target, beams::AbstractVector{<:Beam}; targets,
                   hardware = nothing, enclosed = 0.99, margin = 0.0, level::Integer = 2,
                   op::Union{Nothing,DepositionOperator} = nothing)
    tgt = Symbol[targets...]
    hw = hardware === nothing ?
         [c.name for c in t.components if c.kind === :surface && !(c.name in tgt)] :
         Symbol[hardware...]
    tgt_idx = Set(component_index(t, n) for n in tgt)
    hw_idx = [component_index(t, n) for n in hw]
    # sample points per hardware component
    samples = Dict{Int,Vector{Vec3}}()
    for ci in hw_idx
        pts = Vec3[]
        for k in eachindex(t.tri_area)
            t.tri_comp[k] == ci || continue
            v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
            _subtri_vertices!(pts, v0, v0 + e1, v0 + e2, level)
        end
        samples[ci] = unique(pts)
    end
    frac = op === nothing ? raytrace_component_power(t, beams; nspot = 61) :
           beam_component_power(op, beams)
    rows = Vector{ClearanceRow}(undef, length(beams) * length(hw_idx))
    Threads.@threads for bi in eachindex(beams)
        b = beams[bi]
        env = BeamEnvelope(b; enclosed, margin)
        s_stop = _landing(t, env, tgt_idx)
        P = ForwardDiff.value(b.power)
        for (j, ci) in enumerate(hw_idx)
            best, bp, bs = Inf, zero(Vec3), NaN
            for x in samples[ci]
                d, s = envelope_distance(env, x)
                s <= s_stop || continue
                d < best && ((best, bp, bs) = (d, x, s))
            end
            rows[(bi-1)*length(hw_idx)+j] = ClearanceRow(b.name, t.components[ci].name, best, bp, bs,
                                                         P > 0 ? frac[bi, ci] / P : 0.0)
        end
    end
    return ClearanceReport(rows, Float64(enclosed), Float64(margin), tgt)
end

# Axial position (relative to best focus) where the beam axis first meets a target
# component; the aim plane if it never does.
function _landing(t::Target, env::BeamEnvelope, tgt_idx)
    far = 4 * (t.scale + norm(env.aim)) + 1
    o = env.aim - far * env.k
    traveled = 0.0
    while true
        τ, tri = closest_hit(t.bvh, o, env.k)
        tri == 0 && break
        x = o + τ * env.k
        t.tri_comp[tri] in tgt_idx && return dot(x - env.focus, env.k)
        o = x + _ray_eps(t) * env.k
        traveled += τ
        traveled > 2far && break
    end
    return dot(env.aim - env.focus, env.k)
end

function _subtri_vertices!(pts, a, b, c, level)
    if level == 0
        push!(pts, a, b, c)
        return
    end
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    _subtri_vertices!(pts, a, ab, ca, level - 1)
    _subtri_vertices!(pts, ab, b, bc, level - 1)
    _subtri_vertices!(pts, ca, bc, c, level - 1)
    _subtri_vertices!(pts, ab, bc, ca, level - 1)
end

"""
    beam_component_power(op, beams) -> Matrix

Absorbed power (TW) per beam (rows) and component (columns), both sides of each
component, from the deposition operator.
"""
function beam_component_power(op::DepositionOperator, beams::AbstractVector{<:Beam})
    t = op.target
    nc = length(t.components)
    M = zeros(length(beams), nc)
    Threads.@threads for bi in eachindex(beams)
        E = zeros(nelements(t))
        L = zeros(nc)
        _deposit_beam!(E, L, op, bi, BeamFrame(remake(beams[bi]; aim = ForwardDiff.value.(beams[bi].aim),
                                                       power = ForwardDiff.value(beams[bi].power),
                                                       defocus = ForwardDiff.value(beams[bi].defocus),
                                                       spot = _value_spot(beams[bi].spot))))
        for e in eachindex(E)
            M[bi, t.elem_comp[e]] += E[e]
        end
        M[bi, :] .+= L
    end
    return M
end

_value_spot(s::SuperGaussianSpot) =
    SuperGaussianSpot(ForwardDiff.value(s.a), ForwardDiff.value(s.b); order = s.order, angle = s.angle)

"""
    violations(report, required = 0.0) -> Vector{ClearanceRow}

Rows whose clearance is below `required` mm, closest first.
"""
violations(r::ClearanceReport, required = 0.0) =
    sort(filter(x -> x.clearance < required, r.rows); by = x -> x.clearance)

"""
    min_clearance(report) -> Dict{Symbol,ClearanceRow}

Closest beam for each hardware component.
"""
function min_clearance(r::ClearanceReport)
    out = Dict{Symbol,ClearanceRow}()
    for row in r.rows
        (!haskey(out, row.component) || row.clearance < out[row.component].clearance) &&
            (out[row.component] = row)
    end
    return out
end

function Base.show(io::IO, ::MIME"text/plain", r::ClearanceReport)
    nb = length(unique(x.beam for x in r.rows))
    @printf(io, "ClearanceReport: %d beams, envelope = %.3g%% enclosed + f-cone + %.3g mm margin\n",
            nb, 100r.enclosed, r.margin)
    println(io, "  targets: ", join(r.targets, ", "))
    println(io, "  component        closest beam     clearance (mm)   at s (mm)   max intercepted")
    mc = min_clearance(r)
    for (c, row) in sort(collect(mc); by = x -> x[2].clearance)
        maxi = maximum(x.intercepted for x in r.rows if x.component === c)
        @printf(io, "  %-16s %-16s %+12.3f   %+10.2f   %10.3g%%\n", c, row.beam, row.clearance,
                row.s, 100maxi)
    end
    nv = length(violations(r))
    nv > 0 && println(io, "  ", nv, " beam/component pairs overlap (clearance < 0)")
end

"""
    write_clearance_csv(path, report)

Write every row of a clearance report as CSV.
"""
function write_clearance_csv(path::AbstractString, r::ClearanceReport)
    open(path, "w") do io
        println(io, "beam,component,clearance_mm,x_mm,y_mm,z_mm,s_mm,intercepted_fraction")
        for x in r.rows
            @printf(io, "%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", x.beam, x.component,
                    x.clearance, x.point..., x.s, x.intercepted)
        end
    end
    return path
end
