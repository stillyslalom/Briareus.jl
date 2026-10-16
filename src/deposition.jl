# Laser deposition: a differentiable quadrature operator and an energy-conserving ray trace.

"""
    DepositionOperator(target, beams; level = 1, ncone = 13, resolution = 0.3,
                       reach = 2.0, maxlevel = 5, shadow = 0)

Precomputed laser-deposition quadrature. For each beam, every opaque triangle side that
faces the beam is sampled at `4^level` points, or more where the beam can land: a
triangle within `reach` (mm) of the beam's footprint (2.5 spot semi-axes grown by the
f-cone along the beam axis through the aim) is divided until its sub-triangles are no
longer than `resolution` times the spot's smaller semi-axis, up to `4^maxlevel` points;
with `shadow > 0`, such a triangle crossed by a shadow edge (its points see different
cone directions) is divided `shadow` times more. `reach` is the aim motion the operator
stays accurate for; beyond it, deposition is resolved only at `level`. The defaults
conserve each beam's power to 1e-4 when nothing clips and to 2e-3 when an entrance-hole
lip clips several percent of it, on fine and coarse hohlraum meshes alike.

Each point stores its area weight, the normal of the lit side and a bitmask of which of
the `ncone` f-cone directions reach it from the lens without hitting another surface.
These depend only on the target and the beam *directions*, so the operator stays valid
as pointing, defocus, spot size, and power change — deposition is then a smooth
function of those parameters (see `deposit`).

The beam is the superposition of `ncone` collimated sub-beams, one per f-cone
direction `dⱼ` (`cone_directions`), each carrying the best-focus profile. A point
receives from each visible direction the focal profile where the ray along `dⱼ`
through it crosses the best-focus plane, times that ray's incidence factor
`|n·dⱼ|/(k·dⱼ)`. This is the geometric-optics convolution of the focal spot with the
f-cone, sampled at the cone directions as the reference ray trace samples it, so the
power on any surface that intercepts the whole beam is the beam power, at any
distance from focus and any incidence; and shadows cast by objects near focus — a
capsule at TCC, an LEH edge — are correlated with the part of the spot they block.

The beams' ports must not change after construction.
"""
struct DepositionOperator
    target::Target
    ports::Vector{Vec3}
    dirs::Vector{Vector{Vec3}}     # cone directions per beam
    pts::Vector{Vector{Vec3}}
    w::Vector{Vector{Float64}}     # area weight (mm²)
    normal::Vector{Vector{Vec3}}   # normal of the lit side
    mask::Vector{Vector{UInt16}}   # visible cone directions (bit j-1 set: direction j)
    dest::Vector{Vector{Int32}}    # > 0: element; < 0: -(component) for an unmodeled back side
end

function DepositionOperator(t::Target, beams::AbstractVector{<:Beam}; level::Integer = 1,
                            ncone::Integer = 13, resolution = 0.3, reach = 2.0,
                            maxlevel::Integer = 5, shadow::Integer = 0)
    1 <= ncone <= 16 || throw(ArgumentError("ncone must be in 1:16"))
    maxlevel = max(maxlevel, level)
    nb = length(beams)
    dirs = Vector{Vector{Vec3}}(undef, nb)
    pts = Vector{Vector{Vec3}}(undef, nb)
    w = Vector{Vector{Float64}}(undef, nb)
    normal = Vector{Vector{Vec3}}(undef, nb)
    mask = Vector{Vector{UInt16}}(undef, nb)
    dest = Vector{Vector{Int32}}(undef, nb)
    # sub-triangle quadratures in barycentric form, shared by all triangles
    rules = [_bary_quadrature(l) for l in 0:maxlevel]
    ε = _ray_eps(t)
    Threads.@threads for bi in 1:nb
        b = beams[bi]
        k = direction(b)
        ds, _ = cone_directions(b, ncone)
        sp = b.spot
        amin = ForwardDiff.value(min(sp.a, sp.b))
        extent = 2.5 * ForwardDiff.value(max(sp.a, sp.b)) * (sp.halo > 0 ? max(1, sp.halo_scale) : 1)
        aim = Vec3(ForwardDiff.value.(b.aim))
        focus = aim + ForwardDiff.value(b.defocus) * k
        cone = 1 / (2b.fnumber)
        P, W, N, M, D = Vec3[], Float64[], Vec3[], UInt16[], Int32[]
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
            lv = level
            L = max(norm(e1), norm(e2), norm(e2 - e1))
            cen = v0 + (e1 + e2) / 3
            rad = max(norm(v0 - cen), norm(v0 + e1 - cen), norm(v0 + e2 - cen))
            d = cen - focus
            s = dot(d, k)
            near = norm(d - s * k) - rad <= extent + cone * abs(s) + reach
            if near && L > resolution * amin * 2^level
                lv = clamp(ceil(Int, log2(L / (resolution * amin))), level, maxlevel)
            end
            xs, ms = _visible_directions(t, rules[lv+1][1], v0, e1, e2, side_n, ds, ε)
            # a shadow edge crosses the triangle: resolve it more finely
            if near && shadow > 0 && lv < maxlevel && any(!=(ms[1]), ms)
                lv = min(lv + shadow, maxlevel)
                xs, ms = _visible_directions(t, rules[lv+1][1], v0, e1, e2, side_n, ds, ε)
            end
            for (x, m, wq) in zip(xs, ms, rules[lv+1][2])
                m == 0 && continue
                push!(P, x); push!(W, A * wq); push!(N, side_n); push!(M, m); push!(D, target_id)
            end
        end
        dirs[bi], pts[bi], w[bi], normal[bi], mask[bi], dest[bi] = ds, P, W, N, M, D
    end
    return DepositionOperator(t, [b.port for b in beams], dirs, pts, w, normal, mask, dest)
end

# Quadrature points of a triangle and, for each, the bitmask of cone directions that
# reach it unobstructed.
function _visible_directions(t::Target, bary, v0, e1, e2, side_n, ds, ε)
    xs = Vec3[]
    ms = UInt16[]
    for λ in bary
        x = v0 + λ[1] * e1 + λ[2] * e2
        o = x + ε * side_n
        m = UInt16(0)
        for (j, d) in enumerate(ds)
            dot(d, side_n) < 0 || continue
            any_hit(t.bvh, o, -d) || (m |= UInt16(1) << (j - 1))
        end
        push!(xs, x)
        push!(ms, m)
    end
    return xs, ms
end

function _bary_quadrature(level)
    pts = Vec3[]
    ws = Float64[]
    _subtri!(pts, ws, Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), level)
    return [(p[1], p[2]) for p in pts], ws ./ sum(ws)
end

_fullmask(ds) = UInt16((1 << length(ds)) - 1)

"Operator restricted to the beams `idx` (same order as `idx`)."
Base.getindex(op::DepositionOperator, idx::AbstractVector{<:Integer}) =
    DepositionOperator(op.target, op.ports[idx], op.dirs[idx], op.pts[idx], op.w[idx],
                       op.normal[idx], op.mask[idx], op.dest[idx])

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
    # transverse slopes of the cone directions in the beam frame (fixed by the port)
    κ = [dot(d, f.k) for d in ds]
    α = [dot(d, f.e1) / dot(d, f.k) for d in ds]
    β = [dot(d, f.e2) / dot(d, f.k) for d in ds]
    P, W, Nr, M, D = op.pts[bi], op.w[bi], op.normal[bi], op.mask[bi], op.dest[bi]
    acc = zero(T)
    @inbounds for q in eachindex(W)
        val = W[q] * _cone_intensity(f, ds, α, β, κ, M[q], P[q], Nr[q])
        d = D[q]
        d > 0 ? (Ec[d] += val) : (Lc[-d] += val)
        acc += val
    end
    return acc
end

# Absorbed power per unit area at `x` (lit-side normal `n`): the mean over the visible
# cone directions dⱼ of the best-focus profile where the ray along dⱼ through `x`
# crosses the focal plane, times the ray's incidence factor |n·dⱼ|/(k·dⱼ). Each
# direction's sub-beam is a bundle of parallel rays, so the visibility mask that
# selects it is exactly consistent with where its light goes, and power is conserved
# at shadow edges. (Giving each sub-beam a share of the f-cone spread would smooth the
# footprint but break that consistency: 0.6% of a beam clipped by an LEH edge is lost.)
@inline function _cone_intensity(f::BeamFrame, ds, α, β, κ, m::UInt16, x, n)
    d = x - f.focus
    s = dot(d, f.k)
    u0 = dot(d, f.e1)
    v0 = dot(d, f.e2)
    acc = zero(f.core.power_norm * s)
    for j in eachindex(ds)
        (m >> (j - 1)) & 0x1 == 0x1 || continue
        u = u0 - s * α[j]                     # crossing of the best-focus plane
        v = v0 - s * β[j]
        I = _sg_focal(f.core, u, v)
        f.has_halo && (I += _sg_focal(f.halo, u, v))
        acc += I * (abs(dot(n, ds[j])) / κ[j])
    end
    return acc / length(ds)
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
    raytrace_deposition(target, beams; nspot = 41, ncone = 13, extent = nothing) -> Deposition

Reference deposition by geometric ray tracing. Each spot component (core, and halo if
present) is sampled with an `nspot × nspot` grid of best-focus points spanning
`±extent` 1/e semi-axes (default: out to 1e-12 of the component's peak), each
launched along `ncone` f-cone directions and weighted by the focal intensity. Rays
travel straight (no refraction) to the first opaque surface. Power is conserved
exactly; defocus is exact within the sampling. Not differentiable.
"""
function raytrace_deposition(t::Target, beams::AbstractVector{<:Beam}; nspot::Integer = 41,
                             ncone::Integer = 13, extent = nothing)
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
    raytrace_component_power(target, beams; nspot = 41, ncone = 13, extent = nothing) -> Matrix

Power (TW) each beam (rows) deposits on each component (columns, both sides), by the
reference ray trace — exact power accounting independent of mesh resolution.
"""
function raytrace_component_power(t::Target, beams::AbstractVector{<:Beam}; nspot::Integer = 41,
                                  ncone::Integer = 13, extent = nothing)
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

# Each spot component (core, halo) is sampled on its own nspot × nspot grid spanning
# `extent` 1/e radii (default: out to where the component falls to 1e-12 of its peak).
function _raytrace_beam!(E, L, t::Target, b::Beam, nspot, ncone, extent)
    val = ForwardDiff.value
    sp = b.spot
    a, bb = val(sp.a), val(sp.b)
    e1, e2 = spot_axes(b)
    k = direction(b)
    focus = Vec3(val.(b.aim)) + val(b.defocus) * k
    dirs, dw = cone_directions(b, ncone)
    far = 4 * (t.scale + norm(focus)) + 1
    P = val(b.power)
    comps = [(1 - sp.halo, 1.0, sp.order)]
    sp.halo > 0 && push!(comps, (sp.halo, sp.halo_scale, sp.halo_order))
    total = 0.0
    for (frac, sc, n) in comps
        ext = extent === nothing ? 1.02 * log(1e12)^(1 / n) : extent
        us = range(-ext, ext; length = nspot)
        wsum = sum(exp(-(u^2 + w^2)^(n / 2)) for u in us, w in us)
        for u in us, w in us
            wp = frac * exp(-(u^2 + w^2)^(n / 2)) / wsum
            wp < 1e-15 && continue
            pf = focus + (sc * u * a) * e1 + (sc * w * bb) * e2
            for (d, wd) in zip(dirs, dw)
                pw = P * wp * wd
                _, tri = closest_hit(t.bvh, pf - far * d, d)
                tri == 0 && continue
                total += pw
                e = dot(d, t.tri_normal[tri]) < 0 ? t.tri_front[tri] : t.tri_back[tri]
                e > 0 ? (E[e] += pw) : (L[t.tri_comp[tri]] += pw)
            end
        end
    end
    return total
end
