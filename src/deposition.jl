# Laser deposition: a differentiable quadrature operator and an energy-conserving ray trace.

# Run `f(i)` for `i in 1:n` on all threads, handing out indices one at a time in order
# (dynamic load balancing; with heterogeneous cores a static split waits on the slowest).
# Each `f(i)` must write only outputs owned by `i` for results to be deterministic.
function _tforeach(f, n::Integer)
    nt = min(Threads.nthreads(), n)
    if nt <= 1
        for i in 1:n
            f(i)
        end
        return nothing
    end
    next = Threads.Atomic{Int}(1)
    @sync for _ in 1:nt
        Threads.@spawn while (i = Threads.atomic_add!(next, 1)) <= n
            f(i)
        end
    end
    return nothing
end

# Each beam's points are sorted by destination (elements ascending, then back sides by
# component). `deposit` splits the work into units that own their outputs: a block of
# consecutive elements takes every beam's points in it (and is the only writer of those
# elements), and each beam's back-side points are cut into pieces with their own partial
# sums. Units are handed to threads dynamically, largest first; every output is summed
# in a fixed order (an element over beams, then points, as a serial loop would), so the
# result is bitwise independent of thread count and scheduling.
struct _DepositPlan
    boff::Matrix{Int}                     # (nblock+1) × nb: first point of each element block; last row: first back-side point
    back::Vector{Tuple{Int,UnitRange{Int}}}  # back-side pieces (beam, points)
    order::Vector{Int}                    # units by decreasing size; ≤ nblock: element block, else back-side piece
end

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
    plan::_DepositPlan             # work units of the threaded `deposit`
end

const _DEPOSIT_UNIT = 8192   # points per work unit

function _deposit_plan(ne::Integer, dest::Vector{Vector{Int32}}; unit::Integer = _DEPOSIT_UNIT)
    nb = length(dest)
    cnt = zeros(Int, ne)
    for D in dest, d in D
        d > 0 && (cnt[d] += 1)
    end
    bounds = [1]
    acc = 0
    for e in 1:ne
        acc += cnt[e]
        if acc >= unit && e < ne
            push!(bounds, e + 1)
            acc = 0
        end
    end
    push!(bounds, ne + 1)
    nblock = length(bounds) - 1
    boff = Matrix{Int}(undef, nblock + 1, nb)
    back = Tuple{Int,UnitRange{Int}}[]
    for (bi, D) in enumerate(dest)
        q = 1
        for j in 1:nblock+1
            while q <= length(D) && 0 < D[q] < bounds[j]
                q += 1
            end
            boff[j, bi] = q
        end
        issorted(@view D[1:q-1]) && all(<(0), @view D[q:end]) ||
            throw(ArgumentError("operator points are not sorted by destination"))
        nq = length(D) - q + 1
        npiece = cld(nq, unit)
        for p in 1:npiece
            push!(back, (bi, q + (p - 1) * nq ÷ npiece:q + p * nq ÷ npiece - 1))
        end
    end
    cost = vcat([sum(boff[j+1, bi] - boff[j, bi] for bi in 1:nb; init = 0) for j in 1:nblock],
                [length(r) for (_, r) in back])
    return _DepositPlan(boff, back, sortperm(cost; rev = true))
end

DepositionOperator(t::Target, ports, dirs, pts, w, normal, mask, dest) =
    DepositionOperator(t, ports, dirs, pts, w, normal, mask, dest, _deposit_plan(nelements(t), dest))

function DepositionOperator(t::Target, beams::AbstractVector{<:Beam}; level::Integer = 1,
                            ncone::Integer = 13, resolution = 0.3, reach = 2.0,
                            maxlevel::Integer = 5, shadow::Integer = 0)
    1 <= ncone <= 16 || throw(ArgumentError("ncone must be in 1:16"))
    maxlv = max(maxlevel, level)
    nb = length(beams)
    # sub-triangle quadratures in barycentric form, shared by all triangles
    rules = [_bary_quadrature(l) for l in 0:maxlv]
    ε = _ray_eps(t)
    geo = [_operator_beam(b, ncone) for b in beams]
    # Work units are (beam, triangle range): beams differ in cost and the triangles near
    # a beam's footprint are refined, so units are small and handed out dynamically.
    # Each unit's points are stored apart and joined in order: the result does not
    # depend on scheduling.
    tris = findall(t.tri_opaque)
    nchunk = clamp(cld(length(tris), 256), 1, 16)
    ranges = [tris[(i-1)*length(tris)÷nchunk+1:i*length(tris)÷nchunk] for i in 1:nchunk]
    parts = Matrix{_OperatorPoints}(undef, nchunk, nb)
    _tforeach(nchunk * nb) do u
        ic, bi = mod1(u, nchunk), cld(u, nchunk)
        parts[ic, bi] = _operator_points(t, geo[bi], ranges[ic], rules, ε, level, maxlv,
                                         resolution, reach, shadow)
    end
    ne = nelements(t)
    key(d) = d > 0 ? Int(d) : ne - Int(d)
    pts = Vector{Vector{Vec3}}(undef, nb)
    w = Vector{Vector{Float64}}(undef, nb)
    normal = Vector{Vector{Vec3}}(undef, nb)
    mask = Vector{Vector{UInt16}}(undef, nb)
    dest = Vector{Vector{Int32}}(undef, nb)
    _tforeach(nb) do bi
        col = @view parts[:, bi]
        D = reduce(vcat, [p.D for p in col])
        perm = sortperm(key.(D))                      # stable: keeps point order per element
        pts[bi] = reduce(vcat, [p.P for p in col])[perm]
        w[bi] = reduce(vcat, [p.W for p in col])[perm]
        normal[bi] = reduce(vcat, [p.N for p in col])[perm]
        mask[bi] = reduce(vcat, [p.M for p in col])[perm]
        dest[bi] = D[perm]
    end
    return DepositionOperator(t, [b.port for b in beams], [g.ds for g in geo], pts, w, normal,
                              mask, dest)
end

const _OperatorPoints = @NamedTuple{P::Vector{Vec3}, W::Vector{Float64}, N::Vector{Vec3},
                                    M::Vector{UInt16}, D::Vector{Int32}}

# Beam quantities the operator's sampling depends on (Float64 values).
function _operator_beam(b::Beam, ncone)
    k = direction(b)
    ds, _ = cone_directions(b, ncone)
    sp = b.spot
    amin = ForwardDiff.value(min(sp.a, sp.b))
    extent = 2.5 * ForwardDiff.value(max(sp.a, sp.b)) * (sp.halo > 0 ? max(1, sp.halo_scale) : 1)
    aim = Vec3(ForwardDiff.value.(b.aim))
    focus = aim + ForwardDiff.value(b.defocus) * k
    return (; k, ds, amin, extent, focus, cone = 1 / (2b.fnumber))
end

# Lit quadrature points of one beam on the triangles `tris`.
function _operator_points(t::Target, g, tris, rules, ε, level, maxlevel, resolution, reach, shadow)
    (; k, ds, amin, extent, focus, cone) = g
    P, W, N, M, D = Vec3[], Float64[], Vec3[], UInt16[], Int32[]
    for tri in tris
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
    return (; P, W, N, M, D)
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
Threads share the work dynamically, and every sum is taken in a fixed order: the result
is the same, bit for bit, on any number of threads.
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
    plan = op.plan
    nblock, nback = size(plan.boff, 1) - 1, length(plan.back)
    frames = [BeamFrame(b) for b in beams]
    slopes = [_cone_slopes(op.dirs[bi], frames[bi]) for bi in 1:nb]
    E = zeros(T, ne)
    PB = zeros(T, nb, nblock)     # per-beam power by element block
    LB = zeros(T, nc, nback)      # back-side power by piece
    SB = zeros(T, nback)          # piece totals
    _deposit_units!(E, PB, LB, SB, op, frames, slopes)
    # partial sums reduced in a fixed order
    per_beam = [sum(@view PB[bi, :]; init = zero(T)) for bi in 1:nb]
    L = zeros(T, nc)
    for (u, (bi, _)) in enumerate(plan.back)
        per_beam[bi] += SB[u]
        L .+= @view LB[:, u]
    end
    return Deposition(E, L, per_beam, sum(b.power for b in beams))
end

# Function barrier for the threaded loop over `op.plan`'s work units.
function _deposit_units!(E, PB, LB, SB, op::DepositionOperator, frames, slopes)
    plan = op.plan
    nblock = size(plan.boff, 1) - 1
    _tforeach(length(plan.order)) do k
        u = plan.order[k]
        if u <= nblock
            # an element block: every point lands on an element (E is never indexed < 0)
            for bi in axes(plan.boff, 2)
                r = plan.boff[u, bi]:plan.boff[u+1, bi]-1
                isempty(r) || (PB[bi, u] = _deposit_points!(E, E, op, bi, frames[bi], slopes[bi], r))
            end
        else
            v = u - nblock
            bi, r = plan.back[v]
            SB[v] = _deposit_points!(E, view(LB, :, v), op, bi, frames[bi], slopes[bi], r)
        end
    end
end

# Transverse slopes of the cone directions in the beam frame (fixed by the port).
function _cone_slopes(ds, f::BeamFrame)
    κ = [dot(d, f.k) for d in ds]
    α = [dot(d, f.e1) / dot(d, f.k) for d in ds]
    β = [dot(d, f.e2) / dot(d, f.k) for d in ds]
    return (κ, α, β)
end

# Deposit beam `bi`'s points `r` into elements `Ec` and back sides `Lc`; returns their sum.
function _deposit_points!(Ec::AbstractVector{T}, Lc::AbstractVector{T}, op::DepositionOperator,
                          bi, f::BeamFrame, (κ, α, β), r) where {T}
    ds = op.dirs[bi]
    P, W, Nr, M, D = op.pts[bi], op.w[bi], op.normal[bi], op.mask[bi], op.dest[bi]
    acc = zero(T)
    @inbounds for q in r
        val = W[q] * _cone_intensity(f, ds, α, β, κ, M[q], P[q], Nr[q])
        d = D[q]
        d > 0 ? (Ec[d] += val) : (Lc[-d] += val)
        acc += val
    end
    return acc
end

_deposit_beam!(Ec, Lc, op::DepositionOperator, bi, f::BeamFrame) =
    _deposit_points!(Ec, Lc, op, bi, f, _cone_slopes(op.dirs[bi], f), eachindex(op.w[bi]))

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
    _tforeach(nb) do bi
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
    _tforeach(length(beams)) do bi
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
