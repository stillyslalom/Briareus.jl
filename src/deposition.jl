# Laser deposition: a differentiable quadrature operator and an energy-conserving ray trace.

"""
    DepositionOperator(target, beams; level = 1, ncone = 7)

Precomputed laser-deposition quadrature. For each beam, every opaque triangle side that
faces the beam is sampled at `4^level` points; each point stores its area times the
incidence cosine and a bitmask of which of the `ncone` f-cone directions reach it from
the lens without hitting another surface. These depend only on the target and the beam
*directions*, so the operator stays valid as pointing, defocus, spot size, and power
change — deposition is then a smooth function of those parameters (see `deposit`).

Fully visible points use the moment-matched defocused profile (see `intensity`). At
partially shadowed points each visible cone direction is traced back to the best-focus
plane and contributes the focal profile there (widened by its own share of the
aperture), so shadows cast by objects near focus — a capsule at TCC, an LEH edge — are
correlated with the part of the spot they block.

The beams' ports must not change after construction.
"""
struct DepositionOperator
    target::Target
    ports::Vector{Vec3}
    dirs::Vector{Vector{Vec3}}     # cone directions per beam
    pts::Vector{Vector{Vec3}}
    w::Vector{Vector{Float64}}     # area × |cos| (mm²)
    mask::Vector{Vector{UInt16}}   # visible cone directions (bit j-1 set: direction j)
    dest::Vector{Vector{Int32}}    # > 0: element; < 0: -(component) for an unmodeled back side
end

function DepositionOperator(t::Target, beams::AbstractVector{<:Beam}; level::Integer = 1,
                            ncone::Integer = 7)
    1 <= ncone <= 16 || throw(ArgumentError("ncone must be in 1:16"))
    nb = length(beams)
    dirs = Vector{Vector{Vec3}}(undef, nb)
    pts = Vector{Vector{Vec3}}(undef, nb)
    w = Vector{Vector{Float64}}(undef, nb)
    mask = Vector{Vector{UInt16}}(undef, nb)
    dest = Vector{Vector{Int32}}(undef, nb)
    # sub-triangle quadrature in barycentric form, shared by all triangles
    bary, bw = _bary_quadrature(level)
    ε = _ray_eps(t)
    Threads.@threads for bi in 1:nb
        b = beams[bi]
        k = direction(b)
        ds, _ = cone_directions(b, ncone)
        P, W, M, D = Vec3[], Float64[], UInt16[], Int32[]
        for tri in eachindex(t.tri_area)
            t.tri_opaque[tri] || continue
            n = t.tri_normal[tri]
            c = dot(k, n)
            c == 0 && continue
            side_n = c < 0 ? n : -n            # normal of the lit side
            e = c < 0 ? t.tri_front[tri] : t.tri_back[tri]
            target_id = e > 0 ? e : -t.tri_comp[tri]
            A = t.tri_area[tri]
            v0, e1, e2 = t.tri_v0[tri], t.tri_e1[tri], t.tri_e2[tri]
            for (λ, wq) in zip(bary, bw)
                x = v0 + λ[1] * e1 + λ[2] * e2
                o = x + ε * side_n
                m = UInt16(0)
                for (j, d) in enumerate(ds)
                    dot(d, side_n) < 0 || continue
                    any_hit(t.bvh, o, -d) || (m |= UInt16(1) << (j - 1))
                end
                m == 0 && continue
                push!(P, x); push!(W, A * wq * abs(c)); push!(M, m); push!(D, target_id)
            end
        end
        dirs[bi], pts[bi], w[bi], mask[bi], dest[bi] = ds, P, W, M, D
    end
    return DepositionOperator(t, [b.port for b in beams], dirs, pts, w, mask, dest)
end

function _bary_quadrature(level)
    pts = Vec3[]
    ws = Float64[]
    _subtri!(pts, ws, Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), level)
    return [(p[1], p[2]) for p in pts], ws ./ sum(ws)
end

_fullmask(ds) = UInt16((1 << length(ds)) - 1)

function Base.show(io::IO, op::DepositionOperator)
    np = sum(length, op.w; init = 0)
    npart = sum(count(!=(_fullmask(ds)), M) for (M, ds) in zip(op.mask, op.dirs); init = 0)
    print(io, "DepositionOperator(", length(op.ports), " beams, ", np,
          " lit quadrature points, ", npart, " partially shadowed)")
end

"""
    Deposition

Absorbed laser power by element (`element`, TW), power striking the unmodeled back side
of each component (`backside`, TW — e.g. the outside of a hohlraum, i.e. LEH clipping),
and per-beam totals (`per_beam`, TW, power that landed on any opaque surface).
"""
struct Deposition{T}
    element::Vector{T}
    backside::Vector{T}
    per_beam::Vector{T}
    delivered::T
end

"""
    deposit(op, beams) -> Deposition

Evaluate laser deposition for `beams` (same ports and order as at construction). The
result is differentiable in every beam's `aim`, `power`, `defocus` and spot axes.
Power that strikes no opaque surface (e.g. passes through a hohlraum and out the
other LEH) is `delivered - sum(per_beam)`.
"""
function deposit(op::DepositionOperator, beams::AbstractVector{<:Beam})
    length(beams) == length(op.ports) || throw(DimensionMismatch("beam count changed"))
    for (b, p) in zip(beams, op.ports)
        b.port == p || throw(ArgumentError("beam $(b.name) port differs from the operator's"))
    end
    T = mapreduce(b -> typeof(b.power), promote_type, beams)
    t = op.target
    ne, nc = nelements(t), length(t.components)
    nb = length(beams)
    nchunk = min(nb, Threads.nthreads())
    chunks = [round(Int, (i - 1) * nb / nchunk) + 1:round(Int, i * nb / nchunk) for i in 1:nchunk]
    E = [zeros(T, ne) for _ in 1:nchunk]
    L = [zeros(T, nc) for _ in 1:nchunk]
    per_beam = zeros(T, nb)
    Threads.@threads for ic in 1:nchunk
        for bi in chunks[ic]
            per_beam[bi] = _deposit_beam!(E[ic], L[ic], op, bi, BeamFrame(beams[bi]))
        end
    end
    return Deposition(sum(E), sum(L), per_beam, sum(b.power for b in beams))
end

function _deposit_beam!(Ec::Vector{T}, Lc::Vector{T}, op::DepositionOperator, bi,
                        f::BeamFrame) where {T}
    ds = op.dirs[bi]
    full = _fullmask(ds)
    P, W, M, D = op.pts[bi], op.w[bi], op.mask[bi], op.dest[bi]
    acc = zero(T)
    @inbounds for q in eachindex(W)
        I = M[q] == full ? intensity(f, P[q]) : _partial_intensity(f, ds, M[q], P[q])
        val = W[q] * I
        d = D[q]
        d > 0 ? (Ec[d] += val) : (Lc[-d] += val)
        acc += val
    end
    return acc
end

# Sum over visible cone directions of the focal profile at each ray's best-focus
# crossing, each widened by its 1/N share of the aperture's second moment.
@inline function _partial_intensity(f::BeamFrame, ds, m::UInt16, x)
    N = length(ds)
    d = x - f.focus
    s = dot(d, f.k)
    c = f.defocus_c / N
    acc = zero(f.a2 * s)
    for j in 1:N
        (m >> (j - 1)) & 0x1 == 0x1 || continue
        dj = ds[j]
        sj = s / dot(dj, f.k)
        xf = d - sj * dj                      # crossing of the best-focus plane
        u = dot(xf, f.e1)
        v = dot(xf, f.e2)
        a2 = f.a2 + c * sj^2
        b2 = f.b2 + c * sj^2
        acc += exp(-(u^2 / a2 + v^2 / b2)^f.halfn) / sqrt(a2 * b2)
    end
    return f.power_norm * acc / N
end

"""
    laser_flux(dep, target) -> Vector

Absorbed laser intensity per element (TW/mm²; multiply by 1e14 for W/cm²).
"""
laser_flux(d::Deposition, t::Target) = d.element ./ t.elem_area

"""
    component_power(dep, target) -> Dict{Symbol,Float64}

Absorbed laser power (TW) per component, plus back-side (outside) losses keyed as
`Symbol(name, "_outside")` and `:missed` for power that strikes nothing.
"""
function component_power(d::Deposition, t::Target)
    out = Dict{Symbol,Float64}()
    for (ci, c) in enumerate(t.components)
        c.kind === :aperture && continue
        out[c.name] = sum(ForwardDiff.value(d.element[e]) for e in t.comp_elems[ci]; init = 0.0)
        bs = ForwardDiff.value(d.backside[ci])
        bs > 0 && (out[Symbol(c.name, :_outside)] = bs)
    end
    out[:missed] = ForwardDiff.value(d.delivered - sum(d.per_beam))
    return out
end

# ---------------------------------------------------------------------------------------
# Ray trace (reference)
# ---------------------------------------------------------------------------------------

"""
    raytrace_deposition(target, beams; nspot = 41, ncone = 7, extent = 3.5) -> Deposition

Reference deposition by geometric ray tracing. Each beam is sampled with an
`nspot × nspot` grid of best-focus points spanning `±extent` 1/e semi-axes, each
launched along `ncone` f-cone directions and weighted by the focal intensity. Rays
travel straight (no refraction) to the first opaque surface. Power is conserved
exactly; defocus is exact within the sampling. Not differentiable.
"""
function raytrace_deposition(t::Target, beams::AbstractVector{<:Beam}; nspot::Integer = 41,
                             ncone::Integer = 7, extent = 3.5)
    ne, nc, nb = nelements(t), length(t.components), length(beams)
    E = [zeros(ne) for _ in 1:nb]
    L = [zeros(nc) for _ in 1:nb]
    per_beam = zeros(nb)
    Threads.@threads for bi in 1:nb
        per_beam[bi] = _raytrace_beam!(E[bi], L[bi], t, beams[bi], nspot, ncone, extent)
    end
    return Deposition(sum(E), sum(L), per_beam, sum(ForwardDiff.value(b.power) for b in beams))
end

"""
    raytrace_component_power(target, beams; nspot = 41, ncone = 7, extent = 3.5) -> Matrix

Power (TW) each beam (rows) deposits on each component (columns, both sides), by the
reference ray trace — exact power accounting independent of mesh resolution.
"""
function raytrace_component_power(t::Target, beams::AbstractVector{<:Beam}; nspot::Integer = 41,
                                  ncone::Integer = 7, extent = 3.5)
    nc = length(t.components)
    M = zeros(length(beams), nc)
    Threads.@threads for bi in eachindex(beams)
        E = zeros(nelements(t))
        L = zeros(nc)
        _raytrace_beam!(E, L, t, beams[bi], nspot, ncone, extent)
        for e in eachindex(E)
            M[bi, t.elem_comp[e]] += E[e]
        end
        M[bi, :] .+= L
    end
    return M
end

function _raytrace_beam!(E, L, t::Target, b::Beam, nspot, ncone, extent)
    us = range(-extent, extent; length = nspot)
    val = ForwardDiff.value
    a, bb, n = val(b.spot.a), val(b.spot.b), b.spot.order
    e1, e2 = spot_axes(b)
    k = direction(b)
    focus = Vec3(val.(b.aim)) + val(b.defocus) * k
    dirs, dw = cone_directions(b, ncone)
    wsum = sum(exp(-(u^2 + w^2)^(n / 2)) for u in us, w in us)
    far = 4 * (t.scale + norm(focus)) + 1
    P = val(b.power)
    total = 0.0
    for u in us, w in us
        wp = exp(-(u^2 + w^2)^(n / 2)) / wsum
        wp < 1e-12 && continue
        pf = focus + (u * a) * e1 + (w * bb) * e2
        for (d, wd) in zip(dirs, dw)
            pw = P * wp * wd
            _, tri = closest_hit(t.bvh, pf - far * d, d)
            tri == 0 && continue
            total += pw
            e = dot(d, t.tri_normal[tri]) < 0 ? t.tri_front[tri] : t.tri_back[tri]
            e > 0 ? (E[e] += pw) : (L[t.tri_comp[tri]] += pw)
        end
    end
    return total
end
