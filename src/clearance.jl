# Beam clearance: beam envelopes at an intensity contour, worst-case hardware clearance
# over pointing and focus tolerances, and worst-case power and intensity on hardware.

"""
    enclosed_radius(order, fraction)

Radius (in units of the 1/e semi-axes) of the super-Gaussian contour enclosing
`fraction` of the power: the fraction inside `ρ` is `P(2/n, ρⁿ)` (regularized lower
incomplete gamma function). For spots with a halo use `envelope_radius`.
"""
enclosed_radius(n::Real, f::Real) = gamma_inc_inv(2 / n, f, 1 - f)^(1 / n)

"""
    envelope_radius(spot, power; intensity = nothing, peak_fraction = nothing,
                    enclosed = nothing) -> ρ

Best-focus contour radius in units of the core 1/e semi-axes (halo included), set by
exactly one of: an absolute `intensity` (W/cm²), a fraction of the peak focal intensity
`peak_fraction`, or the `enclosed` power fraction. Returns 0 if the whole spot is below
an intensity threshold.
"""
function envelope_radius(sp::SuperGaussianSpot, P; intensity = nothing, peak_fraction = nothing,
                         enclosed = nothing)
    count(!isnothing, (intensity, peak_fraction, enclosed)) == 1 ||
        throw(ArgumentError("give exactly one of intensity, peak_fraction, enclosed"))
    sp = _value_spot(sp)
    P = ForwardDiff.value(P)
    if enclosed !== nothing
        0 < enclosed < 1 || throw(ArgumentError("enclosed must be in (0, 1)"))
        g = ρ -> enclosed_fraction(sp, ρ) - enclosed            # increasing
    else
        thr = intensity !== nothing ? intensity / 1e14 : peak_fraction * focal_profile(sp, P, 0.0)
        focal_profile(sp, P, 0.0) <= thr && return 0.0
        g = ρ -> thr - focal_profile(sp, P, ρ)                  # increasing
    end
    lo, hi = 0.0, 1.0
    while g(hi) < 0
        hi *= 2
        hi > 1e6 && throw(ArgumentError("envelope contour not found"))
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        g(mid) < 0 ? (lo = mid) : (hi = mid)
        hi - lo < 1e-12 * hi && break
    end
    return hi
end

"""
    BeamEnvelope(beam; intensity, peak_fraction, enclosed, margin = 0.0)

Geometric beam envelope: at axial distance `s` from best focus, the best-focus contour
ellipse (see `envelope_radius`) grown by the f-cone radius `|s|/(2F)` and a uniform
`margin` (mm). In geometric optics with a phase plate (every lens point fills the
whole spot), the intensity at a point is the mean of the focal intensity over a disk of
radius `|s|/(2F)`, so outside the envelope it is below the contour level at every `s`.
Clearance is the signed distance to the envelope; negative means inside.
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

function BeamEnvelope(b::Beam; intensity = nothing, peak_fraction = nothing,
                      enclosed = nothing, margin = 0.0)
    if all(isnothing, (intensity, peak_fraction, enclosed))
        peak_fraction = 1e-4
    end
    ρ = envelope_radius(b.spot, b.power; intensity, peak_fraction, enclosed)
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
# inside and outside. A degenerate ellipse (a = b = 0) is the point at the origin.
function _ellipse_distance(u, v, a, b)
    (a == 0 || b == 0) && return hypot(u, v)
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

One beam/hardware-component pair. All "worst" values are over the pointing and focus
tolerance:

- `clearance` — worst-case envelope clearance (mm; negative = overlap)
- `nominal` — clearance at nominal pointing and focus
- `point`, `s` — hardware point of minimum clearance and its axial position from best
  focus (mm)
- `power`, `intercepted` — worst-case power (TW) landing on the component and its
  fraction of the beam power
- `peak_intensity` — worst-case beam intensity (W/cm², normal to the beam) at the
  component's upstream surface points, ignoring shadowing (conservative)
"""
struct ClearanceRow
    beam::String
    component::Symbol
    clearance::Float64
    nominal::Float64
    point::Vec3
    s::Float64
    power::Float64
    intercepted::Float64
    peak_intensity::Float64
end

"""
    ClearanceReport

Result of `clearance`: one `ClearanceRow` per beam and hardware component plus the
envelope and tolerance settings. See `violations`.
"""
struct ClearanceReport
    rows::Vector{ClearanceRow}
    contour::String
    pointing::Float64
    focus::Float64
    margin::Float64
    targets::Vector{Symbol}
end

"""
    clearance(target, beams; targets, hardware = all other surfaces,
              intensity = nothing, peak_fraction = 1e-4, enclosed = nothing,
              pointing = 0.05, focus = 0.0, margin = 0.0, nangles = 8,
              level = 2, op = nothing, nspot = 61) -> ClearanceReport

Beam-to-hardware clearance check.

Envelope: the best-focus contour at an absolute `intensity` (W/cm²), a `peak_fraction`
of the peak, or an `enclosed` power fraction (give one), grown by the f-cone and
`margin` (see `BeamEnvelope`).

Tolerances: `pointing` is the radius (mm) of the transverse aim error disk, `focus` the
± focal-position error (mm) along the beam. Worst-case clearance is exact for these —
the nominal clearance minus `pointing` minus `focus/(2F)`. (For an independent ±δ box
in two transverse axes, pass `pointing = √2·δ`.) Worst-case power and peak intensity
are maxima over aims displaced to `nangles` points on the tolerance circle (plus
nominal) at focus offsets `-focus, 0, +focus`.

Only hardware upstream of the beam's landing is checked: where the beam axis first
meets one of the `targets` components (the surfaces the beam is meant to hit, e.g.
`(:wall,)`), or the aim plane if it misses them. Hardware is sampled at triangle
vertices and `4^level` sub-triangle points.

Power on hardware comes from the reference ray trace (`nspot` grid per spot component,
exact power accounting at any hardware mesh resolution, halo included), or from `op` —
a `DepositionOperator` for these beam ports — when given (faster; accurate when the
hardware mesh resolves the spot).
"""
function clearance(t::Target, beams::AbstractVector{<:Beam}; targets,
                   hardware = nothing, intensity = nothing, peak_fraction = nothing,
                   enclosed = nothing, pointing = 0.05, focus = 0.0, margin = 0.0,
                   nangles::Integer = 8, level::Integer = 2,
                   op::Union{Nothing,DepositionOperator} = nothing, nspot::Integer = 61)
    all(isnothing, (intensity, peak_fraction, enclosed)) && (peak_fraction = 1e-4)
    contour = intensity !== nothing ? @sprintf("%.3g W/cm²", intensity) :
              peak_fraction !== nothing ? @sprintf("%.3g of peak", peak_fraction) :
              @sprintf("%.6g enclosed", enclosed)
    tgt = Symbol[targets...]
    hw = hardware === nothing ?
         [c.name for c in t.components if c.kind === :surface && !(c.name in tgt)] :
         Symbol[hardware...]
    tgt_idx = Set(component_index(t, n) for n in tgt)
    hw_idx = [component_index(t, n) for n in hw]
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

    # beam variants over the tolerance region
    nb = length(beams)
    variants = Beam{Float64}[]
    owner = Int[]
    for (bi, b) in enumerate(beams)
        bv = remake(b; aim = ForwardDiff.value.(b.aim), power = ForwardDiff.value(b.power),
                    defocus = ForwardDiff.value(b.defocus), spot = _value_spot(b.spot),
                    dlambda = ForwardDiff.value(b.dlambda))
        e1, e2 = spot_axes(bv)
        shifts = pointing > 0 ?
                 vcat([zero(Vec3)], [pointing * (cos(2π * j / nangles) * e1 + sin(2π * j / nangles) * e2)
                                     for j in 0:nangles-1]) : [zero(Vec3)]
        for δ in shifts, df in (focus > 0 ? (0.0, -focus, focus) : (0.0,))
            push!(variants, remake(bv; aim = Vec3(bv.aim) + δ, defocus = bv.defocus + df))
            push!(owner, bi)
        end
    end
    Pc = op === nothing ? raytrace_component_power(t, variants; nspot) :
         beam_component_power(op, variants)

    rows = Vector{ClearanceRow}(undef, nb * length(hw_idx))
    Threads.@threads for bi in 1:nb
        vids = findall(==(bi), owner)
        b = variants[vids[1]]                       # nominal variant comes first
        env = BeamEnvelope(b; intensity, peak_fraction, enclosed, margin)
        s_stop = _landing(t, env, tgt_idx)
        frames = [BeamFrame(variants[v]) for v in vids]
        P = b.power
        for (j, ci) in enumerate(hw_idx)
            best, bp, bs = Inf, zero(Vec3), NaN
            Imax = 0.0
            for x in samples[ci]
                d, s = envelope_distance(env, x)
                s <= s_stop || continue
                d < best && ((best, bp, bs) = (d, x, s))
                for f in frames
                    Imax = max(Imax, intensity_at(f, x))
                end
            end
            worst = best - pointing - focus * env.cone
            pw = maximum(Pc[v, ci] for v in vids)
            rows[(bi-1)*length(hw_idx)+j] = ClearanceRow(b.name, t.components[ci].name, worst, best,
                                                         bp, bs, pw, P > 0 ? pw / P : 0.0,
                                                         to_Wcm2(Imax))
        end
    end
    return ClearanceReport(rows, contour, Float64(pointing), Float64(focus), Float64(margin), tgt)
end

intensity_at(f::BeamFrame, x) = intensity(f, x)

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
component, from the deposition operator. `beams` may differ from the operator's in
pointing, focus, spot and power but must share its ports, in the same order — or be
any list whose ports all appear among the operator's (each is matched by port).
"""
function beam_component_power(op::DepositionOperator, beams::AbstractVector{<:Beam})
    t = op.target
    nc = length(t.components)
    M = zeros(length(beams), nc)
    Threads.@threads for i in eachindex(beams)
        b = beams[i]
        bi = findfirst(==(b.port), op.ports)
        bi === nothing && throw(ArgumentError("beam $(b.name) port is not in the operator"))
        E = zeros(nelements(t))
        L = zeros(nc)
        bv = remake(b; aim = ForwardDiff.value.(b.aim), power = ForwardDiff.value(b.power),
                    defocus = ForwardDiff.value(b.defocus), spot = _value_spot(b.spot),
                    dlambda = ForwardDiff.value(b.dlambda))
        _deposit_beam!(E, L, op, bi, BeamFrame(bv))
        for e in eachindex(E)
            M[i, t.elem_comp[e]] += E[e]
        end
        M[i, :] .+= L
    end
    return M
end

_value_spot(s::SuperGaussianSpot) =
    with_axes(s, ForwardDiff.value(s.a), ForwardDiff.value(s.b))

"""
    violations(report; clearance = 0.0, power = Inf, intensity = Inf) -> Vector{ClearanceRow}

Rows failing any criterion: worst-case clearance below `clearance` mm, worst-case power
on the component above `power` TW, or worst-case intensity above `intensity` W/cm².
Sorted by clearance.
"""
violations(r::ClearanceReport; clearance = 0.0, power = Inf, intensity = Inf) =
    sort(filter(x -> x.clearance < clearance || x.power > power || x.peak_intensity > intensity,
                r.rows); by = x -> x.clearance)

"""
    min_clearance(report) -> Dict{Symbol,ClearanceRow}

Closest beam (worst-case clearance) for each hardware component.
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
    @printf(io, "ClearanceReport: %d beams; envelope at %s + f-cone + %.3g mm margin\n",
            nb, r.contour, r.margin)
    @printf(io, "  tolerances: pointing %.3g mm (transverse radius), focus ±%.3g mm; targets: %s\n",
            r.pointing, r.focus, join(r.targets, ", "))
    println(io, "  component        closest beam   worst (nominal) clearance mm   max power on it   max intensity")
    mc = min_clearance(r)
    for (c, row) in sort(collect(mc); by = x -> x[2].clearance)
        rs = [x for x in r.rows if x.component === c]
        pmax = maximum(x.power for x in rs)
        imax = maximum(x.peak_intensity for x in rs)
        @printf(io, "  %-16s %-14s %+9.3f (%+9.3f)          %9.3g GW      %9.3g W/cm²\n", c,
                row.beam, row.clearance, row.nominal, 1e3pmax, imax)
    end
    nv = count(x -> x.clearance < 0, r.rows)
    nv > 0 && println(io, "  ", nv, " beam/component pairs overlap the envelope in the worst case")
end

"""
    write_clearance_csv(path, report)

Write every row of a clearance report as CSV.
"""
function write_clearance_csv(path::AbstractString, r::ClearanceReport)
    open(path, "w") do io
        println(io, "beam,component,worst_clearance_mm,nominal_clearance_mm,x_mm,y_mm,z_mm,s_mm,",
                "worst_power_TW,worst_intercepted_fraction,worst_peak_intensity_Wcm2")
        for x in r.rows
            println(io, join((x.beam, x.component,
                              (@sprintf("%.6g", v) for v in (x.clearance, x.nominal, x.point..., x.s,
                                                              x.power, x.intercepted, x.peak_intensity))...), ","))
        end
    end
    return path
end
