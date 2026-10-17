# Materials, target components, and the assembled `Target` (elements + ray-casting BVH).

"""
    AlbedoModel

Supertype of x-ray albedo models: `ConstantAlbedo`, `TabulatedAlbedo`,
`WallLossAlbedo`, `BurnthroughAlbedo`. Evaluate with `albedo`.
"""
abstract type AlbedoModel end

"""
    ConstantAlbedo(α)

Time-independent x-ray albedo (re-emitted / incident flux).
"""
struct ConstantAlbedo <: AlbedoModel
    value::Float64
end

"""
    TabulatedAlbedo(times, values)

Albedo interpolated linearly in time (ns) and held constant outside the table.
"""
struct TabulatedAlbedo <: AlbedoModel
    t::Vector{Float64}
    value::Vector{Float64}
    function TabulatedAlbedo(t, v)
        length(t) == length(v) || throw(ArgumentError("times and values differ in length"))
        issorted(t) || throw(ArgumentError("times must be sorted"))
        new(collect(Float64, t), collect(Float64, v))
    end
end

"""
    WallLossAlbedo(; coefficient = 0.32, temperature_exponent = -0.7, time_exponent = -0.38,
                   reference_temperature = 250.0, t0 = 0.1, range = (0.0, 0.95))

Marshak-wave wall-loss albedo, `1 - α = c (T/100 eV)^p (t/1 ns)^q`, with defaults for
gold from Lindl's scaling (T in heV, t in ns). `t` is clamped below at `t0`, and α to
`range`. The radiation temperature is each element's incident-flux temperature, so
`solve` iterates it to self-consistency; `reference_temperature` (eV) is used where
no flux is known.
"""
Base.@kwdef struct WallLossAlbedo <: AlbedoModel
    coefficient::Float64 = 0.32
    temperature_exponent::Float64 = -0.7
    time_exponent::Float64 = -0.38
    reference_temperature::Float64 = 250.0
    t0::Float64 = 0.1
    range::Tuple{Float64,Float64} = (0.0, 0.95)
end

"""
    albedo(model, t, T = nothing)
    albedo(material, t = 0.0)

X-ray albedo at time `t` (ns); temperature-dependent models use the radiation
temperature `T` (eV) when given.
"""
albedo(m::ConstantAlbedo, t, T = nothing) = m.value
albedo(m::TabulatedAlbedo, t, T = nothing) = _interp(m.t, m.value, t)
function albedo(m::WallLossAlbedo, t, T = nothing)
    Te = (T === nothing || T <= 0) ? m.reference_temperature : T
    loss = m.coefficient * (Te / 100)^m.temperature_exponent * max(t, m.t0)^m.time_exponent
    return clamp(1 - loss, m.range...)
end

"""
    HeatFront(; coefficient = 9.51e-4, temperature_exponent = 1.914, time_exponent = 0.5156)

Areal mass (g/cm²) behind the ablative (subsonic) Marshak heat front that radiation of
temperature `T` drives into a gold wall: at constant `T`,

    m_F = c (T/100 eV)^a (t/1 ns)^b,

with the defaults from the self-similar solution of Hammer & Rosen, Phys. Plasmas 10,
1829 (2003), eq. 112 for a constant drive (q₀ = 0), built on their power-law fits to
gold's opacity and equation of state. For a drive `T(t)` the front is integrated as

    m_F(t) = c [∫₀ᵗ (T(t′)/100 eV)^(a/b) dt′]^b,

which is exact at constant `T`, depends only on the drive history, and for a ramp
`T ∝ tᵏ` gives `m_F = c (1 + k a/b)^(-b) T(t)^a t^b`: within 2% of Hammer & Rosen's
exact coefficient for `k = 0.182` and within 1% of Shussman & Heizler's (Phys. Plasmas
22, 082109 (2015)) for `k = 0.123` (see `heat_front` and the diagnostic-windows
tutorial). Lindl et al. (Phys. Plasmas 11, 339 (2004), eq. 4-24) give
`c = 1.20e-3, a = 1.86, b = 0.54`, the scaling the `WallLossAlbedo` defaults come
from; it heats about 20% more mass at 300 eV, 1 ns. `T` is the incident-flux
radiation temperature (Lindl's convention; Hammer & Rosen drive with the surface
temperature, a few percent lower for gold late in a pulse).
"""
Base.@kwdef struct HeatFront
    coefficient::Float64 = 9.51e-4
    temperature_exponent::Float64 = 1.914
    time_exponent::Float64 = 0.5156
end

# The integrand (T/100 eV)^(a/b) of the front's drive integral, and the areal mass for
# a value of the integral.
_front_rate(f::HeatFront, T) = (max(T, zero(T)) / 100)^(f.temperature_exponent / f.time_exponent)
_front_mass(f::HeatFront, I) = f.coefficient * max(I, zero(I))^f.time_exponent
_front_integral(f::HeatFront, m) = (m / f.coefficient)^(1 / f.time_exponent)

"""
    heat_front(front::HeatFront, times, Tr) -> Vector

Areal mass (g/cm²) behind the heat front at each of `times` (ns) for the radiation
temperature history `Tr` (eV) sampled at those times, piecewise linear between them
(trapezoidal drive integral), from zero before `times[1]`.
"""
function heat_front(f::HeatFront, times, Tr)
    length(times) == length(Tr) || throw(ArgumentError("times and temperatures differ in length"))
    I = 0.0
    m = zeros(length(times))
    for k in 2:length(times)
        I += (times[k] - times[k-1]) * (_front_rate(f, Tr[k-1]) + _front_rate(f, Tr[k])) / 2
        m[k] = _front_mass(f, I)
    end
    return m
end

"""
    eddington_transmission(τ)

Fraction of the incident flux that a gray, planar foil of optical depth `τ` transmits
in the Eddington approximation with a Marshak boundary condition,

    F_trans / F_inc = (4/3) / (τ + 2/3) = 4 / (3τ + 2),

eq. 2 of R. S. Lester et al., Phys. Rev. E 114, L023201 (2026), doi:10.1103/k8dz-pjdx.
It is a diffusion result for optically thick foils: it exceeds 1 below `τ = 2/3`,
where it is capped at 1. `BurnthroughAlbedo` takes its default transmission from it
(see `optical_depth`). Generic in the number type.
"""
function eddington_transmission(τ)
    τ >= 0 || throw(DomainError(τ, "the optical depth must not be negative"))
    T = 4 / (3τ + 2)
    return min(one(T), T)
end

# The point the default foil optical depth is scaled from: τ ∼ 5 for a gold foil
# 2.5 µm thick (Lester et al. 2026, after eq. 1), i.e. 19.3 g/cm³ × 2.5 µm of gold.
const _FOIL_REFERENCE_DEPTH = 5.0
const _FOIL_REFERENCE_AREAL_MASS = 19.3 * 2.5e-4      # g/cm²

_foil_optical_depth(areal_mass) = _FOIL_REFERENCE_DEPTH * areal_mass / _FOIL_REFERENCE_AREAL_MASS

"""
    BurnthroughAlbedo(; thickness, density = 19.3, before = WallLossAlbedo(),
                      transmission = eddington_transmission(τ), after = nothing,
                      width = 0.4, burnthrough = nothing, front = HeatFront())

Albedo of a thin wall layer (a thin-foil diagnostic window: `thickness` µm of a
material of `density` g/cm³, gold by default) that re-emits like the wall until the
Marshak heat front (`front`, a `HeatFront`) has heated its whole areal mass
`ρ × thickness`, and then lets part of the flux that reaches it through.

Once burned through, a fraction `T = transmission` of the incident flux passes through
the foil and leaves the hohlraum. The rest meets the foil as before and is re-emitted
with the albedo of `before`. Across the transition

    α(t) = α_before (1 - s T),    transmitted fraction s T,

with `s = s((t - t_b)/width)` the smooth step `3x² - 2x³` on [0, 1] and `t_b` the
burnthrough time, so the burned-through albedo is `(1 - T) α_before`. Of the flux
`H A` on the foil, `s T H A` is transmitted and `(1 - s T)(1 - α_before) H A` absorbed
(`energy_balance`, `drive_history`). `before` is an albedo model or a number, by
default the wall's Marshak loss (`WallLossAlbedo`).

The default `T` is the gray-foil estimate of R. S. Lester et al., Phys. Rev. E 114,
L023201 (2026), doi:10.1103/k8dz-pjdx: `T = 4/(3τ + 2)` (their eq. 2,
`eddington_transmission`). The optical depth is their eq. 1, `τ = L/λ_Au(T_h)` at a
fixed Rosseland mean free path λ_Au, scaled linearly in areal mass (`optical_depth`)
from their point `τ ∼ 5` at `L = 2.5 µm`:

    τ = 5 × ρ thickness / (19.3 g/cm³ × 2.5 µm),

so `T ≈ 0.33`, `0.235` and `0.17` for 1.7, 2.5 and 3.6 µm of gold. Limitations:
- λ_Au is held at their evaluation, for a hohlraum radiation temperature `T_h` of
  typically 300–350 eV and an assumed gold density of about 1 g/cm³ (the expanded
  foil). Its dependence on temperature and density is not modeled, and `T` is
  constant once burned through.
- Pairing that density with the 2.5 µm as-built thickness is theirs. Taken literally,
  eq. 1 at 1 g/cm³ describes a layer of 0.25 mg/cm², not the foil's 4.8. The scaling
  with thickness (areal mass) is ours, from eq. 1.
- The estimate is gray, planar and quasi-static, and for gold: pass `transmission`
  for another material. Below `τ = 2/3` (0.33 µm of gold) it is capped at 1.

`transmission = 1` makes a burned-through foil an open hole whose whole incident flux
is transmitted.

An explicit `after` (an albedo model or number; exclusive with `transmission`) sets
the burned-through albedo directly instead,

    α(t) = α_before + s (α_after - α_before),

and counts as transmitted only the loss beyond the intact foil's, `(α_before - α) H A`;
the intact foil's own `(1 - α_before) H A` stays absorbed. `after = 0` is the open-hole
bound on the loss with this split (the default before `transmission` existed).

The burnthrough time `t_b` is, in order:
- the per-element time passed to `solve`/`albedos` (`burnthrough = ...`), which
  `drive_history` integrates from each element's own incident-flux temperature
  history (`heat_front`);
- the explicit `burnthrough` (ns), when given;
- otherwise the foil is intact (`t_b = ∞`): a single solve has no history to estimate
  it from. `burnthrough_time` gives the estimate for a temperature history.

Before `t_b` the foil loses exactly what a thick wall loses and transmits nothing: in
the diffusion picture the cold material ahead of a nonlinear heat front does not
change the flux into it, so no early loss is modeled. Measured burnthrough of 0.8–3 µm
gold foils on a 260 eV hohlraum wall (Nova; Lindl et al., Phys. Plasmas 11, 339
(2004), Fig. 4-6) reaches half signal 0.1–0.25 ns after the computed front arrival,
0.19 ns on average, which sets the default `width` (half transmission at
`t_b + width/2`; see the diagnostic-windows tutorial).
"""
struct BurnthroughAlbedo <: AlbedoModel
    thickness::Float64
    density::Float64
    before::AlbedoModel
    after::Union{Nothing,AlbedoModel}       # an explicit burned-through albedo, or
    transmission::Union{Nothing,Float64}    # the transmitted fraction (exactly one is set)
    width::Float64
    burnthrough::Union{Nothing,Float64}
    front::HeatFront
    function BurnthroughAlbedo(thickness, density, before, after, transmission, width,
                               burnthrough, front)
        thickness > 0 && density > 0 || throw(ArgumentError("thickness and density must be positive"))
        width >= 0 || throw(ArgumentError("the transition width must not be negative"))
        after === nothing || transmission === nothing ||
            throw(ArgumentError("give `after` or `transmission`, not both"))
        if after === nothing
            transmission === nothing &&
                (transmission = eddington_transmission(_foil_optical_depth(density * thickness * 1e-4)))
            0 <= transmission <= 1 || throw(ArgumentError("transmission must lie in [0, 1]"))
        end
        new(thickness, density, _albedo_model(before),
            after === nothing ? nothing : _albedo_model(after),
            transmission === nothing ? nothing : Float64(transmission), width,
            burnthrough === nothing ? nothing : Float64(burnthrough), front)
    end
end
_albedo_model(x) = x isa AlbedoModel ? x : ConstantAlbedo(x)

BurnthroughAlbedo(; thickness, density = 19.3, before = WallLossAlbedo(), after = nothing,
                  transmission = nothing, width = 0.4, burnthrough = nothing,
                  front = HeatFront()) =
    BurnthroughAlbedo(thickness, density, before, after, transmission, width, burnthrough, front)

"Areal mass (g/cm²) of a burnthrough layer."
areal_mass(m::BurnthroughAlbedo) = m.density * m.thickness * 1e-4

"""
    optical_depth(model::BurnthroughAlbedo)

Gray optical depth of a burned-through gold foil, scaled linearly in areal mass from
`τ ∼ 5` at 2.5 µm (Lester et al., Phys. Rev. E 114, L023201 (2026), eq. 1 at a fixed
Rosseland mean free path): `5 × areal_mass(model) / (19.3 g/cm³ × 2.5 µm)`. The
default `transmission` is `eddington_transmission(optical_depth(model))`. An explicit
`transmission` or `after` does not change this estimate.
"""
optical_depth(m::BurnthroughAlbedo) = _foil_optical_depth(areal_mass(m))

# Transmitting fraction of the transition at time t for burnthrough time tb.
function _burnthrough_step(m::BurnthroughAlbedo, t, tb)
    isfinite(tb) || return 0.0
    m.width == 0 && return t >= tb ? 1.0 : 0.0
    x = clamp((t - tb) / m.width, 0.0, 1.0)
    return x^2 * (3 - 2x)
end

# The albedo at time t, the intact (`before`) albedo, and the fraction of the incident
# flux transmitted out of the enclosure.
function _burnthrough_split(m::BurnthroughAlbedo, t, T, burnthrough)
    tb = burnthrough !== nothing && !isnan(burnthrough) ? burnthrough :
         m.burnthrough !== nothing ? m.burnthrough : Inf
    αb = albedo(m.before, t, T)
    s = _burnthrough_step(m, t, tb)
    s == 0 && return αb, αb, zero(αb)
    if m.after === nothing
        τ = s * m.transmission
        return αb * (1 - τ), αb, τ + zero(αb)
    end
    α = αb + s * (albedo(m.after, t, T) - αb)
    return α, αb, αb - α
end

albedo(m::BurnthroughAlbedo, t, T = nothing; burnthrough = nothing) =
    first(_burnthrough_split(m, t, T, burnthrough))

"""
    burnthrough_time(model::BurnthroughAlbedo, times, Tr) -> Float64

Time (ns) at which the heat front driven by the radiation temperature history `Tr`
(eV, at `times`, piecewise linear; `heat_front`) reaches the back of the layer, or `Inf`
if it does not within `times`. This is the estimate alone; an explicit `burnthrough`
of the model is not consulted.
"""
function burnthrough_time(m::BurnthroughAlbedo, times, Tr)
    length(times) == length(Tr) || throw(ArgumentError("times and temperatures differ in length"))
    f = m.front
    need = _front_integral(f, areal_mass(m))
    I = 0.0
    for k in 2:length(times)
        r0, r1 = _front_rate(f, Tr[k-1]), _front_rate(f, Tr[k])
        Δ = times[k] - times[k-1]
        dI = Δ * (r0 + r1) / 2
        if I + dI >= need
            # the integrand is linear across the step: solve the quadratic for the crossing
            return times[k-1] + _linear_crossing(r0, r1, Δ, need - I)
        end
        I += dI
    end
    return Inf
end

# x in [0, Δ] with ∫₀ˣ (r0 + (r1 - r0) s/Δ) ds = J (0 ≤ J ≤ Δ(r0 + r1)/2).
function _linear_crossing(r0, r1, Δ, J)
    a = (r1 - r0) / 2Δ
    abs(a) * Δ^2 <= 1e-12 * max(J, eps()) && return r0 > 0 ? J / r0 : Δ
    # a x² + r0 x - J = 0, root in [0, Δ], written to avoid cancellation
    return 2J / (r0 + sqrt(max(r0^2 + 4a * J, 0.0)))
end

function _interp(ts, vs, t)
    t <= ts[1] && return vs[1]
    t >= ts[end] && return vs[end]
    k = searchsortedlast(ts, t)
    w = (t - ts[k]) / (ts[k+1] - ts[k])
    return (1 - w) * vs[k] + w * vs[k+1]
end

"""
    Material(name; albedo = 0.0, conversion = 0.0)

Surface properties for the radiation model. `albedo` (a number or `AlbedoModel`) is the
x-ray re-emission fraction; `conversion` is the fraction of absorbed laser power
re-emitted as x-rays.
"""
struct Material
    name::Symbol
    albedo::AlbedoModel
    conversion::Float64
end

Material(name::Symbol; albedo = 0.0, conversion = 0.0) =
    Material(name, albedo isa AlbedoModel ? albedo : ConstantAlbedo(albedo), conversion)

albedo(m::Material, t = 0.0) = albedo(m.albedo, t)

"Gold hohlraum wall (defaults: albedo 0.8, laser→x-ray conversion 0.75)."
gold(; albedo = 0.8, conversion = 0.75) = Material(:Au; albedo, conversion)
"Capsule ablator (defaults: albedo 0.3, conversion 0.3)."
ablator(; albedo = 0.3, conversion = 0.3) = Material(:ablator; albedo, conversion)
"Perfect absorber that emits nothing."
const BLACK = Material(:black; albedo = 0.0, conversion = 0.0)

"""
    Component(name, mesh, material; kind = :surface, two_sided = false, surface = nothing)

A named piece of a target. `kind = :surface` components are opaque: they block laser
rays and lines of sight, absorb laser light, and radiate from their front side (and
back side when `two_sided`). `kind = :aperture` components (e.g. laser entrance holes)
are transparent diagnostic surfaces that record the x-ray flux crossing them and
re-emit nothing. `surface` is the exact `CurvedSurface` the mesh approximates, if any;
the radiation model then uses the curved geometry (see `view_factors`).
"""
struct Component
    name::Symbol
    mesh::Mesh
    material::Material
    kind::Symbol
    two_sided::Bool
    surface::Union{Nothing,CurvedSurface}
    function Component(name, mesh, material; kind::Symbol = :surface, two_sided::Bool = false,
                       surface::Union{Nothing,CurvedSurface} = nothing)
        kind in (:surface, :aperture) || throw(ArgumentError("kind must be :surface or :aperture"))
        kind === :aperture && two_sided && throw(ArgumentError("apertures are one-sided"))
        new(Symbol(name), mesh, material, kind, two_sided, surface)
    end
end

"""
    aperture(name, mesh)

Transparent diagnostic surface (e.g. a laser entrance hole): rays pass through it and
its incident x-ray flux is the flux escaping through it.
"""
aperture(name, mesh; surface = nothing) = Component(name, mesh, BLACK; kind = :aperture, surface)

"Component with its mesh rigidly moved (`x ↦ R*x + t`)."
transform(c::Component, R::AbstractMatrix, t = zero(Vec3)) =
    Component(c.name, transform(c.mesh, R, Vec3(t...)), c.material; kind = c.kind,
              two_sided = c.two_sided, surface = _transform_surface(c.surface, R, Vec3(t...)))
translate(c::Component, t) = transform(c, SMatrix{3,3,Float64}(I), t)

"""
    Target(components)

Assembled target geometry. Each patch side that radiates becomes an *element* (the
unknowns of the radiation model): one element per patch for one-sided surfaces and
apertures, two for two-sided surfaces. The back side of a one-sided surface absorbs
laser light as a loss but is not modeled radiatively.
"""
struct Target
    components::Vector{Component}
    # triangles
    tri_v0::Vector{Vec3}
    tri_e1::Vector{Vec3}
    tri_e2::Vector{Vec3}
    tri_normal::Vector{Vec3}
    tri_area::Vector{Float64}
    tri_comp::Vector{Int32}
    tri_front::Vector{Int32}       # element on the front side (0 = none)
    tri_back::Vector{Int32}        # element on the back side (0 = none)
    tri_opaque::BitVector
    # elements
    elem_centroid::Vector{Vec3}
    elem_normal::Vector{Vec3}      # unit normal of the radiating side
    elem_area::Vector{Float64}
    elem_comp::Vector{Int32}
    elem_side::Vector{Int8}        # +1 front, -1 back
    elem_tris::Vector{UnitRange{Int}}
    elem_curved::Vector{Int8}      # radiating side relative to the curved surface's outward normal (0: flat)
    comp_elems::Vector{Vector{Int}}
    bvh::BVH
    scale::Float64                 # bounding-sphere radius about the origin (mm)
end

function Target(components::AbstractVector{Component})
    names = [c.name for c in components]
    allunique(names) || throw(ArgumentError("component names must be unique: $names"))
    tv0, te1, te2, tn, ta = Vec3[], Vec3[], Vec3[], Vec3[], Float64[]
    tcomp, tfront, tback = Int32[], Int32[], Int32[]
    topaque = BitVector()
    ec, en, ea, ecomp, eside, ecurv = Vec3[], Vec3[], Float64[], Int32[], Int8[], Int8[]
    etris = UnitRange{Int}[]
    comp_elems = Vector{Int}[]
    for (ci, c) in enumerate(components)
        m = c.mesh
        t_offset = length(tv0)
        for (a, b, d) in m.triangles
            v0, v1, v2 = m.vertices[a], m.vertices[b], m.vertices[d]
            cr = cross(v1 - v0, v2 - v0)
            push!(tv0, v0); push!(te1, v1 - v0); push!(te2, v2 - v0)
            ar = norm(cr) / 2
            push!(ta, ar); push!(tn, ar > 0 ? cr / (2ar) : zero(Vec3))
            push!(tcomp, ci); push!(tfront, 0); push!(tback, 0)
            push!(topaque, c.kind === :surface)
        end
        mine = Int[]
        for side in (c.two_sided ? (1, -1) : (1,))
            for r in m.patches
                gr = r .+ t_offset
                A = sum(ta[gr])
                A > 0 || continue
                cen = sum(ta[k] * (tv0[k] + (te1[k] + te2[k]) / 3) for k in gr) / A
                nrm = normalize(sum(ta[k] * tn[k] for k in gr))
                curv = Int8(0)
                if c.surface !== nothing
                    # radiation uses the curved patch: its area, and normals on the side
                    # the flat element radiates from
                    A = c.surface isa DiskSurface ?
                        _sector_area(_sector(c.surface, _patch_vertices(tv0, te1, te2, gr))) :
                        sum(_curved_area(c.surface, tv0[k], tv0[k] + te1[k], tv0[k] + te2[k]) for k in gr)
                    curv = Int8(sign(dot(side * nrm, last(surface_point(c.surface, cen)))))
                end
                push!(ec, cen); push!(en, side * nrm); push!(ea, A); push!(ecurv, curv)
                push!(ecomp, ci); push!(eside, side); push!(etris, gr)
                e = length(ec)
                push!(mine, e)
                for k in gr
                    side == 1 ? (tfront[k] = e) : (tback[k] = e)
                end
            end
        end
        push!(comp_elems, mine)
    end
    opaque_idx = findall(topaque)
    bvh = _sah_bvh(BVH(tv0[opaque_idx], te1[opaque_idx], te2[opaque_idx], Int32.(opaque_idx)))
    scale = isempty(tv0) ? 1.0 :
            maximum(max(norm(v), norm(v + a), norm(v + b)) for (v, a, b) in zip(tv0, te1, te2))
    return Target(collect(components), tv0, te1, te2, tn, ta, tcomp, tfront, tback, topaque,
                  ec, en, ea, ecomp, eside, etris, ecurv, comp_elems, bvh, scale)
end

Target(cs::Component...) = Target(collect(cs))

"Number of radiating elements (patch sides) of a target."
nelements(t::Target) = length(t.elem_area)
ntriangles(t::Target) = length(t.tri_area)

"""
    component_index(target, name) -> Int
"""
function component_index(t::Target, name::Symbol)
    i = findfirst(c -> c.name === name, t.components)
    i === nothing && throw(KeyError(name))
    return i
end

"""
    elements(target, names...) -> Vector{Int}

Element indices of the named components (both sides of two-sided components).
"""
elements(t::Target, names::Symbol...) =
    reduce(vcat, (t.comp_elems[component_index(t, n)] for n in names); init = Int[])

"""
    element_material(target, e)
"""
element_material(t::Target, e::Integer) = t.components[t.elem_comp[e]].material
is_aperture(t::Target, e::Integer) = t.components[t.elem_comp[e]].kind === :aperture

function Base.show(io::IO, ::MIME"text/plain", t::Target)
    println(io, "Target with ", length(t.components), " components, ", nelements(t),
            " elements, ", ntriangles(t), " triangles")
    for (i, c) in enumerate(t.components)
        A = sum(t.elem_area[e] for e in t.comp_elems[i] if t.elem_side[e] == 1; init = 0.0)
        @printf(io, "  %-14s %-8s %-9s %6d elements  area %9.4f mm²%s\n", c.name, c.kind,
                c.material.name, length(t.comp_elems[i]), A, c.two_sided ? "  (two-sided)" : "")
    end
end
Base.show(io::IO, t::Target) = print(io, "Target(", length(t.components), " components, ",
                                     nelements(t), " elements)")

# Small offset used to start rays just off a surface (mm).
_ray_eps(t::Target) = 1e-7 * max(t.scale, 1.0)

"""
    visible(target, a, b) -> Bool

Whether the open segment between points `a` and `b` is free of opaque surfaces.
"""
visible(t::Target, a::Vec3, b::Vec3) = !any_hit(t.bvh, a, b - a, 1 - 1e-9)

"""
    element_quadrature(target, e, level) -> (points, weights, normals)

Quadrature over element `e` from `4^level` sub-triangles per triangle: sub-triangle
centroids weighted by area, with the element normal. On a component with a
`CurvedSurface`, the points are projected onto the surface and carry its normals
(on the element's radiating side) and curved-patch areas. Weights sum to the element
area.
"""
function element_quadrature(t::Target, e::Integer, level::Integer)
    pts = Vec3[]
    ws = Float64[]
    for k in t.elem_tris[e]
        _subtri!(pts, ws, t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k], t.tri_v0[k] + t.tri_e2[k], level)
    end
    s = t.components[t.elem_comp[e]].surface
    if s === nothing || t.elem_curved[e] == 0
        return pts, ws, fill(t.elem_normal[e], length(pts))
    end
    if s isa DiskSurface
        r0, r1, a, b = element_sector(t, e)
        nr, nφ = 2^level, 2^(level + 1)
        rs = sqrt.(range(r0^2, r1^2; length = nr + 1))        # equal-area rings
        φs = range(a, b; length = nφ + 1)
        cpts, cws = Vec3[], Float64[]
        for i in 1:nr, k in 1:nφ
            push!(cpts, _sector_point(s, sqrt((rs[i]^2 + rs[i+1]^2) / 2), (φs[k] + φs[k+1]) / 2))
            push!(cws, (rs[i+1]^2 - rs[i]^2) / 2 * (φs[k+1] - φs[k]))
        end
        return cpts, cws, fill(t.elem_curved[e] * s.normal, length(cpts))
    end
    cpts, cws, cns = Vec3[], Float64[], Vec3[]
    for k in t.elem_tris[e]
        _curved_subtri!(cpts, cws, cns, s, t.elem_curved[e], t.tri_v0[k], t.tri_v0[k] + t.tri_e1[k],
                        t.tri_v0[k] + t.tri_e2[k], level)
    end
    return cpts, cws, cns
end

_patch_vertices(tv0, te1, te2, gr) =
    unique(vcat([[tv0[k], tv0[k] + te1[k], tv0[k] + te2[k]] for k in gr]...))

"""
    element_sector(target, e) -> (r0, r1, φa, φb)

Polar sector of an element of a component with a `DiskSurface`.
"""
element_sector(t::Target, e::Integer) =
    _sector(t.components[t.elem_comp[e]].surface,
            _patch_vertices(t.tri_v0, t.tri_e1, t.tri_e2, t.elem_tris[e]))

function _curved_subtri!(pts, ws, ns, s, sgn, a, b, c, level)
    if level == 0
        p, n = surface_point(s, (a + b + c) / 3)
        push!(pts, p); push!(ns, sgn * n); push!(ws, _curved_area(s, a, b, c))
        return
    end
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    _curved_subtri!(pts, ws, ns, s, sgn, a, ab, ca, level - 1)
    _curved_subtri!(pts, ws, ns, s, sgn, ab, b, bc, level - 1)
    _curved_subtri!(pts, ws, ns, s, sgn, ca, bc, c, level - 1)
    _curved_subtri!(pts, ws, ns, s, sgn, ab, bc, ca, level - 1)
end

function _subtri!(pts, ws, a, b, c, level)
    if level == 0
        push!(pts, (a + b + c) / 3)
        push!(ws, norm(cross(b - a, c - a)) / 2)
        return
    end
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    _subtri!(pts, ws, a, ab, ca, level - 1)
    _subtri!(pts, ws, ab, b, bc, level - 1)
    _subtri!(pts, ws, ca, bc, c, level - 1)
    _subtri!(pts, ws, ab, bc, ca, level - 1)
end
