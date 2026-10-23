# Moving targets: geometry that changes through the pulse, solved on keyframes.

"""
    WallMotion(; speed = 0.10, start = 0.0)

Inward motion of a hohlraum's radiating wall: the x-ray emission surface sits at radius
`R₀ - speed * max(0, t - start)` (mm, ns), where `R₀` is the initial wall radius and
`start` the outer cones' turn-on.

The default speed is an engineering default for low gas fill (0.3 mg/cc), *derived* from
Chen et al. 2020 (the emitting wall 500 µm in by 5.5 ns at 0.3 mg/cc, 400 µm at 0.45;
about 0.07 mm/ns) and the HYBRID-E wall of PRE 109, 025204 Fig. 1a (0.35–0.4 mm in by
6 ns at the waist); use ~0.05 mm/ns at 0.6 mg/cc. See research notes, "Hohlraum
dynamics".
"""
struct WallMotion
    speed::Float64
    start::Float64
end
WallMotion(; speed::Real = 0.10, start::Real = 0.0) = WallMotion(speed, start)

"Wall radius (mm) at time `t` (ns) of a wall that starts at radius `R0`."
wall_radius(w::WallMotion, R0, t) = R0 - w.speed * max(0.0, t - w.start)

"""
    moving_hohlraum(motion::WallMotion; radius = 2.875, kwargs...) -> (τ -> Target)

The `hohlraum(; radius, kwargs...)` at time `τ` with its wall moved in by `motion`: the
cylinder has radius `wall_radius(motion, radius, τ)`, the caps keep their planes and
the entrance holes their radius (the caps' annuli shrink with the wall), rounded corners
follow the wall, and the capsule stays. The mesh counts are those of the initial radius
(from `h` when given), so every time has the same elements, as [`Keyframes`](@ref)
needs.
"""
function moving_hohlraum(motion::WallMotion; radius::Real = 2.875, length::Real = 9.4,
                         leh_radius::Real = 1.55, h = nothing, corner_radius::Real = 0.0,
                         nz::Integer = 36, nphi::Integer = 48, nr_cap::Integer = 6,
                         nr_leh::Integer = 4, nr_corner = nothing, kwargs...)
    grid = _hohlraum_grid(; length, radius, leh_radius, h, rc = Float64(corner_radius),
                          nz, nphi, nr_cap, nr_leh, nr_corner)
    return τ -> hohlraum(; radius = wall_radius(motion, radius, τ), length, leh_radius,
                         corner_radius, grid..., kwargs...)
end

"""
    Keyframes(build, times, beams; closed = nothing, keep = 2, deposition_kw = (;),
              viewfactor_kw = (;))

A target that changes through the pulse. `build(τ)` makes its `Target` at time `τ`
(ns), with the same elements at every time (e.g. [`moving_hohlraum`](@ref)); the view
factors and the deposition operator of `beams` are computed at the keyframe `times`.
Between two keyframes, [`model_at`](@ref) interpolates the exchange matrix `Aᵢ Fᵢⱼ`
linearly in time and divides by the areas of `build(τ)`, which keeps reciprocity, and
rebalances it when `closed` (default: the target's declaration, as in
`RadiationModel`). Deposition between keyframes uses an operator built at that time
([`deposition_at`](@ref)). Before the first keyframe and after the last, the nearest one
holds. The error against rebuilding at every time falls as the square of the spacing: on
a NIF hohlraum whose wall moves at 0.1 mm/ns, keyframes 2, 1 and 0.5 ns apart keep P2/P0
within 0.035, 0.012 and 0.009 percentage points (research notes, "Moving wall").

Keyframes are built when first needed and only the `keep` most recently used models are
held (a fine hohlraum's view factors take over 1 GB each); deposition operators are
kept. A history walks forward in time, so each keyframe is built once.

`drive_history(kf, beams_at, times; ...)` runs a drive history on the moving target.
"""
struct Keyframes{B}
    build::B
    time::Vector{Float64}
    beams::Vector{Beam}
    closed::Bool
    keep::Int
    deposition_kw::Any
    viewfactor_kw::Any
    reference::Target                              # build(time[1]): the elements
    models::Dict{Int,RadiationModel}
    used::Vector{Int}                              # keyframes, most recently used first
    ops::Dict{Int,DepositionOperator}
    between::Base.RefValue{Any}                    # (τ, model) last interpolated
    between_op::Base.RefValue{Any}                 # (τ, operator) last built between keyframes
    lock::ReentrantLock
end

function Keyframes(build, times, beams::AbstractVector{<:Beam}; closed = nothing,
                   keep::Integer = 2, deposition_kw = (;), viewfactor_kw = (;))
    ts = collect(Float64, times)
    isempty(ts) && throw(ArgumentError("need at least one keyframe"))
    issorted(ts) && allunique(ts) || throw(ArgumentError("keyframe times must increase"))
    keep >= 2 || throw(ArgumentError("keep at least two keyframes"))
    ref = build(ts[1])
    return Keyframes(build, ts, collect(Beam, beams), something(closed, ref.closed), Int(keep),
                     deposition_kw, viewfactor_kw, ref, Dict{Int,RadiationModel}(), Int[],
                     Dict{Int,DepositionOperator}(), Ref{Any}(nothing), Ref{Any}(nothing),
                     ReentrantLock())
end

function Base.show(io::IO, kf::Keyframes)
    print(io, "Keyframes(", length(kf.time), " times from ", kf.time[1], " to ", kf.time[end],
          " ns, ", nelements(kf.reference), " elements)")
end

# A target with the reference's elements: the same components and element ranges.
function _check_elements(kf::Keyframes, t::Target, τ)
    r = kf.reference
    nelements(t) == nelements(r) && t.comp_elems == r.comp_elems &&
        [c.name for c in t.components] == [c.name for c in r.components] ||
        throw(ArgumentError("the target at t = $τ ns has other elements than at $(kf.time[1]) ns"))
    return t
end

# The model and deposition operator of keyframe k, built on first use.
function _keyframe(kf::Keyframes, k::Int)
    lock(kf.lock) do
        m = get(kf.models, k, nothing)
        if m === nothing
            τ = kf.time[k]
            t = k == 1 ? kf.reference : _check_elements(kf, kf.build(τ), τ)
            m = RadiationModel(t; closed = kf.closed, kf.viewfactor_kw...)
            kf.models[k] = m
            haskey(kf.ops, k) || (kf.ops[k] = DepositionOperator(t, kf.beams; kf.deposition_kw...))
        end
        filter!(!=(k), kf.used)
        pushfirst!(kf.used, k)
        while length(kf.used) > kf.keep
            delete!(kf.models, pop!(kf.used))
        end
        return m, kf.ops[k]
    end
end

# Keyframes bracketing τ and the weight of the later one.
function _bracket(ts::Vector{Float64}, τ)
    τ <= ts[1] && return 1, 1, 0.0
    τ >= ts[end] && return length(ts), length(ts), 0.0
    k = searchsortedlast(ts, τ)
    w = (τ - ts[k]) / (ts[k+1] - ts[k])
    return k, k + 1, w
end

"""
    model_at(kf::Keyframes, τ) -> RadiationModel

The radiation model of the moving target at time `τ`: a keyframe's, or between two the
target `build(τ)` with interpolated exchange factors (see [`Keyframes`](@ref)). The last
interpolated model is kept, so repeated calls at one time build it once.
"""
function model_at(kf::Keyframes, τ::Real)
    k0, k1, w = _bracket(kf.time, τ)
    w == 0 && return first(_keyframe(kf, k0))
    b = kf.between[]
    b !== nothing && b[1] == τ && return b[2]
    m0, _ = _keyframe(kf, k0)
    m1, _ = _keyframe(kf, k1)
    t = _check_elements(kf, kf.build(τ), τ)
    A0, A1, A = m0.target.elem_area, m1.target.elem_area, t.elem_area
    F = similar(m0.F)
    Threads.@threads for j in axes(F, 2)
        @inbounds for i in axes(F, 1)
            F[i, j] = ((1 - w) * A0[i] * m0.F[i, j] + w * A1[i] * m1.F[i, j]) / A[i]
        end
    end
    kf.closed && enforce_closure!(F, A)
    m = RadiationModel(t, F)
    kf.between[] = (τ, m)
    return m
end

"""
    deposition_at(kf::Keyframes, τ) -> DepositionOperator

The deposition operator of the keyframes' beams on the moving target at time `τ`: a
keyframe's, or between two one built on `build(τ)` (the last one is kept). Unlike the
view factors it is not interpolated: an operator costs a tenth of the view factors, and
interpolating depositions blurs footprints that slide with the wall (on a NIF hohlraum
at 0.1 mm/ns with six keyframes, P2/P0 off by 0.24 points instead of 0.035).
"""
function deposition_at(kf::Keyframes, τ::Real)
    k0, _, w = _bracket(kf.time, τ)
    w == 0 && return last(_keyframe(kf, k0))
    lock(kf.lock) do
        b = kf.between_op[]
        b !== nothing && b[1] == τ && return b[2]
        op = DepositionOperator(model_at(kf, τ).target, kf.beams; kf.deposition_kw...)
        kf.between_op[] = (τ, op)
        return op
    end
end

"""
    deposit(kf::Keyframes, beams, τ) -> Deposition

Laser deposition of `beams` (the keyframes' ports, in order) on the moving target at time
`τ`, with [`deposition_at`](@ref).
"""
deposit(kf::Keyframes, beams, τ::Real) = deposit(deposition_at(kf, τ), beams)

"""
    drive_history(kf::Keyframes, beams_at, times; kwargs...) -> NamedTuple

[`drive_history`](@ref) on a moving target: each time is solved on `model_at(kf, τ)`
with the deposition `deposit(kf, beams_at(τ), τ)`. Keywords and results are those of the
fixed-target history; areas, moments and spherical harmonics are taken on each time's
target.
"""
drive_history(kf::Keyframes, beams_at, times; kwargs...) =
    _drive_history(kf.reference, τ -> model_at(kf, τ), (τ, beams) -> deposit(kf, beams, τ),
                   beams_at, times; kwargs...)
