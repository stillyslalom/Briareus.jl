# This file defines the pointing parameterizations: `AbstractParameterization` and
# `GroupPointing`, which maps parameters to beams per cone, quad or beam.

const _CACHED_LMAX = 12

const _FIELDS = (:z, :r, :phi, :split, :defocus, :power, :spot, :dlambda,
                 :theta_s, :dphi_s, :split_theta, :split_phi)
const _CYL_FIELDS = (:z, :r, :phi, :split)
const _SPH_FIELDS = (:theta_s, :dphi_s, :split_theta, :split_phi)

"""
    AbstractParameterization

Supertype of the maps from a parameter vector to beams that a [`Scenario`](@ref)
optimizes: [`GroupPointing`](@ref), [`FramePointing`](@ref), [`HeldParameters`](@ref),
[`PowerScales`](@ref) and [`PhasedPointing`](@ref). A subtype has the methods
`apply(p, x)`, which returns the beams of the parameters `x` and is generic in their
number type, and `initial_parameters(p)`. `parameter_names(p)` labels the parameters in
results and selects them by name in [`HeldParameters`](@ref), and `nparams(p)` defaults
to the length of the initial parameters. An object of another type with the first two
methods is also accepted by `Scenario`, and its parameters are then labeled `x[i]`.
"""
abstract type AbstractParameterization end

nparams(p::AbstractParameterization) = length(initial_parameters(p))

# Returns the names of the `n` parameters of `p`: its `parameter_names`, or `x[i]` for an
# object that is not an `AbstractParameterization`.
_parameter_names(p::AbstractParameterization, n) = collect(String, parameter_names(p))
_parameter_names(p, n) = ["x[$i]" for i in 1:n]

"""
    GroupPointing(beams; by = :group, fields = (:z, :r, :defocus), mirror = true,
                  relative = false, unit = :quad, radius = nothing, center = (0, 0, 0))

Maps a parameter vector to beams. Beams are grouped by `by`: `:group` (e.g. a cone),
`:quad`, `:beam`, or a function `beam -> key`. Each group gets one value per field:

- `:z` — axial aim coordinate (mm). With `mirror`, beams from the lower hemisphere
  (port `z < 0`) aim at `-z`, so one value describes a top/bottom-symmetric cone.
- `:r` — radial aim offset (mm) toward the quad port azimuth. With `unit = :beam`
  (below), the offset is toward the beam port azimuth.
- `:phi` — tangential aim offset (mm), counterclockwise about `z`.
- `:split` — tangential split of the beams of a quad about its aim point (mm), along the
  quad φ̂. Each beam aim moves by `±split`, with the sign of the beam port azimuth
  measured from the mean port of the quad; beams counterclockwise of it take `+split`.
  A positive split therefore opens the aims of a quad in the sense of its ports, and a
  negative one crosses them. The `s` of [`hybride_mock_pointing`](@ref Briareus.Benchmarks.hybride_mock_pointing) is
  `split = -s`. The field sets the tangential offset of each beam from the mean
  template aim of its quad, or shifts that offset with `relative`. A quad held as one
  beam is not split. The field needs `unit = :quad` and `quad` names.
- `:defocus` — focal-plane offset along the beam (mm, see `Beam`).
- `:power` — power per beam (TW).
- `:spot` — multiplier on the best-focus spot semi-axes.
- `:dlambda` — wavelength detuning (Å), used by laser–plasma models.

Fields not listed keep the template value of each beam.

**Pointing on a sphere** (polar direct drive). With `radius = R` (mm), the aim can
instead be a point on the sphere of radius `R` about `center`. Direct-drive pointings
are published on this target surface (Barlow et al., PRL 133, 175101, 2024, Table I;
Yeamans et al., Nucl. Fusion 61, 046031, 2021, Table 1; Hohenberger et al., PoP 22,
056308, 2015, Table I). The spherical fields are:

- `:theta_s` — polar angle of the aim point (degrees from `+z`). With `mirror`, beams
  from the lower hemisphere aim at `180° - theta_s`, so one value describes a
  top/bottom-symmetric cone pair, and a larger value moves both toward the equator.
- `:dphi_s` — azimuth of the aim point relative to the beam port azimuth, in degrees
  counterclockwise about `z` and the same for both hemispheres.
- `:split_theta`, `:split_phi` — split of the beams of a quad (degrees). Each beam aim
  moves by `±split_theta` in polar angle and `±split_phi` in azimuth, with the sign of
  the beam port offset from the mean port of the quad. Beams at a larger polar angle
  than the quad port take `+split_theta`, with the angle counted from the pole of their
  hemisphere when mirrored. Beams counterclockwise of the quad port take `+split_phi`.
  Positive values therefore open the aims of a quad in the sense of its ports. This is
  the split-quad pointing of Cok et al. (PoP 15, 082705, 2008) and Yeamans et al. 2021
  (left/right `Δφ`) as angles on the sphere. A quad given as one beam, or a beam level
  with the mean port of its quad, is not split. The fields need `quad` names; other
  per-beam aims come from `by`, e.g. one group per sub-ring.

The aim is the surface point itself, `aim = center + R * port_vector(θ, φ)`, which the
axis meets first when it faces the port (`Beam`). Best focus is then at
`aim + defocus * k`, so defocus is measured from the surface aim point, as VisRad
measures the focus offset from its reference point (research notes, "Pointing and focus
conventions"). The angles of each template beam are those of the first point where its
axis meets the sphere; an `ArgumentError` is thrown if the axis misses. Listed spherical
values are shared by the members of a group. For an angle not set by a listed field, the
value of each beam is kept. The aim moves to the surface along the template axis, so the
focus moves with it unless the template aim was already on the sphere. The spherical fields
and the hohlraum fields `:z`, `:r`, `:phi` exclude each other.

By default (`unit = :quad`), the unit of pointing for the hohlraum fields is a quad.
`:r`, `:phi` and `:split` are then measured in the azimuth of the quad port, which is
the mean port of the beams that share its `quad` name, instead of the port of each beam.
The beams of a quad that share their values then share one aim point, whether the table
contains the quad as one beam or as four. A NIF pointing given per cone or per quad thus
stays one aim per quad while the deposition is computed per beam. The beams of a quad
are pointed apart about its aim point as follows:

```julia
p = GroupPointing(nif_beams(); by = :quad, fields = (:z, :r, :phi, :split))
beams = apply(p, parameters(p, "Q31T" => (z = 5.8, split = 0.5)))
```

`by = :beam` gives each beam a separate aim, with its offsets still in the quad azimuth.
With `relative = true` (below), the differences between the beams of a per-beam table
are kept under one shift per quad. An as-shot SST export,
`read_visrad_sst(path; template = nif_beams(quads))`, is such a table. `unit = :beam`
measures `:r` and `:phi` of each beam in the beam port azimuth instead. The four beams
of a quad with the same `r` then aim `r sin Δφ` apart, where Δφ is the azimuth of a beam
from that of its quad, up to 5.7° on the NIF 23.5° cone. A quad held as one beam, and a
beam without a quad name, give identical results either way.

`relative = true` makes the aim fields (`:z`, `:r`, `:phi`, `:split`, or the spherical
ones) shifts added to the template aim of each beam, instead of values shared by the
beams of a group. A group then moves rigidly, and the differences between its
members are kept. This re-points an as-shot table, whose beams within a quad aim at
different points (the per-beam offsets of a facility), by one shift per quad. The initial
parameters are then zero; the other fields stay absolute.
"""
struct GroupPointing{K} <: AbstractParameterization
    template::Vector{Beam{Float64}}
    keys::Vector{K}
    members::Vector{Vector{Int}}
    fields::Vector{Symbol}
    mirror::Bool
    relative::Bool
    radius::Float64                      # sphere of the spherical fields (NaN: none)
    center::Vec3
    sph::Vector{NTuple{5,Float64}}       # per beam: θ' (mirrored frame), Δφ, sθ, sφ, φ_port
    sph0::Vector{NTuple{4,Float64}}      # per group: mean θ', mean Δφ, θ' split, Δφ split
    unit::Symbol                         # :beam or :quad, the port azimuth of :r and :phi
    radial::Vector{NTuple{2,Vec3}}       # per beam: those radial and tangential directions
    split::Vector{NTuple{2,Float64}}     # per beam with `:split`: its sign, its template offset
end

function GroupPointing(beams::AbstractVector{<:Beam}; by = :group,
                       fields = (:z, :r, :defocus), mirror::Bool = true, relative::Bool = false,
                       unit::Symbol = :quad,
                       radius::Union{Nothing,Real} = nothing, center = (0.0, 0.0, 0.0))
    by === :quad && any(b -> isempty(b.quad), beams) &&
        throw(ArgumentError("by = :quad needs quad names on every beam (OMEGA beams have none)"))
    key = by === :group ? (b -> b.group) : by === :quad ? (b -> b.quad) :
          by === :beam ? (b -> b.name) : by
    return _group_pointing(beams, map(key, beams); fields, mirror, relative, unit, radius, center)
end

# Builds a `GroupPointing` of `beams` with group key `ks[i]` for beam i. The keywords are
# those of the constructor, except `by`.
function _group_pointing(beams::AbstractVector{<:Beam}, ks::AbstractVector; fields, mirror::Bool,
                         relative::Bool, unit::Symbol, radius::Union{Nothing,Real}, center)
    for f in fields
        f in _FIELDS || throw(ArgumentError("unknown field $f; choose from $_FIELDS"))
    end
    unit in (:beam, :quad) || throw(ArgumentError("unit is :beam or :quad, got :$unit"))
    if :split in fields
        unit === :quad || throw(ArgumentError("the :split field is in the quad's frame: it needs unit = :quad"))
        any(b -> isempty(b.quad), beams) && throw(ArgumentError("the :split field needs quad names on every beam"))
    end
    spherical = any(in(_SPH_FIELDS), fields)
    spherical && any(in(_CYL_FIELDS), fields) &&
        throw(ArgumentError("the spherical fields $_SPH_FIELDS exclude $_CYL_FIELDS"))
    spherical && radius === nothing &&
        throw(ArgumentError("the spherical fields $_SPH_FIELDS need the sphere's `radius`"))
    radius === nothing || radius > 0 || throw(ArgumentError("radius must be positive"))
    any(in((:split_theta, :split_phi)), fields) && any(b -> isempty(b.quad), beams) &&
        throw(ArgumentError("split fields need quad names on every beam"))
    uk = unique(ks)
    members = [findall(==(k), ks) for k in uk]
    tmpl = [remake(b; aim = Float64.(ForwardDiff.value.(b.aim)),
                   power = Float64(ForwardDiff.value(b.power)),
                   defocus = Float64(ForwardDiff.value(b.defocus)),
                   spot = with_axes(b.spot, Float64(ForwardDiff.value(b.spot.a)),
                                    Float64(ForwardDiff.value(b.spot.b))),
                   dlambda = Float64(ForwardDiff.value(b.dlambda)))
            for b in beams]
    R = radius === nothing ? NaN : Float64(radius)
    c = Vec3(center...)
    sph = spherical ? _sphere_angles(tmpl, R, c, mirror) : NTuple{5,Float64}[]
    sph0 = spherical ? [_sphere_group0(sph, mem) for mem in members] : NTuple{4,Float64}[]
    radial = unit === :quad ? _quad_radial(tmpl) : map(_radial, tmpl)
    split = :split in fields ? _quad_split(tmpl, radial) : NTuple{2,Float64}[]
    return GroupPointing(tmpl, uk, members, collect(Symbol, fields), mirror, relative, R, c,
                         sph, sph0, unit, radial, split)
end

# Returns, per template beam, the side of its quad that it splits to and its tangential
# offset. The side is the sign of the beam port azimuth measured from the mean port of
# the quad. It is 0 for a beam level with that port and for a quad held as one beam. The
# offset is that of the beam aim from the mean aim of the quad beams, along the quad φ̂
# (the second vector of `radial`).
function _quad_split(tmpl, radial)
    qsum, asum, n = Dict{String,Vec3}(), Dict{String,Vec3}(), Dict{String,Int}()
    for b in tmpl
        qsum[b.quad] = get(qsum, b.quad, zero(Vec3)) + b.port
        asum[b.quad] = get(asum, b.quad, zero(Vec3)) + b.aim
        n[b.quad] = get(n, b.quad, 0) + 1
    end
    return map(eachindex(tmpl)) do i
        b = tmpl[i]
        q = qsum[b.quad]
        δ = mod(atan(b.port[2], b.port[1]) - atan(q[2], q[1]) + π, 2π) - π
        side = abs(δ) < 1e-9 ? 0.0 : sign(δ)
        (side, dot(b.aim - asum[b.quad] / n[b.quad], radial[i][2]))
    end
end

# Returns the tangential offset of each entry of `aims` from the mean of its quad, along
# the quad φ̂. `aims` contains one aim per template beam, or the shifts of the template
# aims.
function _split_offsets(p::GroupPointing, aims::AbstractVector{Vec3})
    asum, n = Dict{String,Vec3}(), Dict{String,Int}()
    for (b, a) in zip(p.template, aims)
        asum[b.quad] = get(asum, b.quad, zero(Vec3)) + a
        n[b.quad] = get(n, b.quad, 0) + 1
    end
    return [dot(a - asum[b.quad] / n[b.quad], p.radial[i][2]) for (i, (b, a)) in enumerate(zip(p.template, aims))]
end

# Returns the split value of a group from the tangential offsets `d` of its members from
# the mean aims of their quads. The value is Σ side × d over the members with a side,
# divided by their number. It is 0 when no member has a side.
function _group_split(p::GroupPointing, mem, d)
    k = count(bi -> p.split[bi][1] != 0, mem)
    k == 0 && return 0.0
    return sum(p.split[bi][1] * d[bi] for bi in mem if p.split[bi][1] != 0) / k
end

# Returns, per template beam, the radial and tangential directions of its quad port. The
# quad port is the mean port of the beams that share the quad name. For a quad held as
# one beam, it is exactly the port of that beam. For a beam without a quad name, the
# directions of its port are used.
function _quad_radial(tmpl)
    qsum = Dict{String,Vec3}()
    for b in tmpl
        isempty(b.quad) || (qsum[b.quad] = get(qsum, b.quad, zero(Vec3)) + b.port)
    end
    return [isempty(b.quad) ? _radial(b) : _radial(qsum[b.quad]) for b in tmpl]
end

# Returns five values per template beam. The first two locate the first point where the
# beam axis meets the sphere: its polar angle θ' in degrees, counted from the pole of the
# beam hemisphere when mirrored, and its azimuth offset Δφ from the beam port. The others
# are the split signs (sθ, sφ) of the beam relative to the mean port of its quad, and the
# port azimuth.
function _sphere_angles(tmpl, R, c, mirror)
    qsum = Dict{String,Vec3}()
    for b in tmpl
        isempty(b.quad) || (qsum[b.quad] = get(qsum, b.quad, zero(Vec3)) + b.port)
    end
    return map(tmpl) do b
        zs = (mirror && b.port[3] < 0) ? -1.0 : 1.0
        k = -b.port
        d = b.aim - c
        dk = dot(d, k)
        disc = dk^2 - (dot(d, d) - R^2)
        disc >= 0 || throw(ArgumentError("beam $(b.name)'s axis misses the sphere of radius $R"))
        s = d + (-dk - sqrt(disc)) * k                  # first hit, relative to the center
        θ = acosd(clamp(zs * s[3] / R, -1, 1))
        φp = atand(b.port[2], b.port[1])
        Δφ = mod(atand(s[2], s[1]) - φp + 180, 360) - 180
        sθ = sφ = 0.0
        if !isempty(b.quad)
            q = qsum[b.quad]
            δθ = acosd(clamp(zs * b.port[3], -1, 1)) - acosd(clamp(zs * q[3] / norm(q), -1, 1))
            δφ = mod(φp - atand(q[2], q[1]) + 180, 360) - 180
            sθ = abs(δθ) < 1e-9 ? 0.0 : sign(δθ)
            sφ = abs(δφ) < 1e-9 ? 0.0 : sign(δφ)
        end
        (θ, Δφ, sθ, sφ, φp)
    end
end

# Returns the initial spherical values of a group: the mean θ' and Δφ, and the mean
# split Σ s (x - x̄) / #(s ≠ 0), which is zero without split members.
function _sphere_group0(sph, mem)
    θm = sum(sph[i][1] for i in mem) / length(mem)
    φm = sum(sph[i][2] for i in mem) / length(mem)
    split(j, xm) = (n = count(i -> sph[i][j+2] != 0, mem);
                    n == 0 ? 0.0 : sum(sph[i][j+2] * (sph[i][j] - xm) for i in mem) / n)
    return (θm, φm, split(1, θm), split(2, φm))
end

"Number of parameters of a `GroupPointing` (groups × fields)."
nparams(p::GroupPointing) = length(p.keys) * length(p.fields)

"Names `\"<group>.<field>\"` for each parameter."
parameter_names(p::GroupPointing) = [string(k, ".", f) for k in p.keys for f in p.fields]

_pidx(p::GroupPointing, g, fi) = (g - 1) * length(p.fields) + fi

_zsign(p::GroupPointing, b::Beam) = (p.mirror && b.port[3] < 0) ? -1.0 : 1.0
_radial(b::Beam) = _radial(b.port)
_radial(port::Vec3) = (φ = atan(port[2], port[1]); (Vec3(cos(φ), sin(φ), 0), Vec3(-sin(φ), cos(φ), 0)))
_spherical(p::GroupPointing) = !isempty(p.sph)

"""
    initial_parameters(p::GroupPointing) -> Vector{Float64}

Parameters reproducing the template beams (group means where members differ; zero
shifts for the aim fields of a `relative` parameterization). Spherical fields take the
angles where each template axis first meets the sphere (see `GroupPointing`).
"""
function initial_parameters(p::GroupPointing)
    x = zeros(nparams(p))
    for (g, mem) in enumerate(p.members), (fi, f) in enumerate(p.fields)
        p.relative && f in _CYL_FIELDS && continue
        if f in _SPH_FIELDS
            p.relative || (x[_pidx(p, g, fi)] = p.sph0[g][findfirst(==(f), _SPH_FIELDS)])
            continue
        end
        # The group split averages the offset of each beam from its quad mean, signed by
        # its side.
        sg = :split in p.fields ? _group_split(p, mem, last.(p.split)) : 0.0
        if f === :split
            x[_pidx(p, g, fi)] = sg
            continue
        end
        vals = map(mem) do bi
            b = p.template[bi]
            er, eφ = p.radial[bi]
            f === :z ? _zsign(p, b) * b.aim[3] :
            f === :r ? dot(b.aim, er) :
            f === :phi ? (isempty(p.split) ? dot(b.aim, eφ) : dot(b.aim, eφ) - p.split[bi][1] * sg) :
            f === :defocus ? b.defocus :
            f === :power ? b.power :
            f === :dlambda ? b.dlambda : 1.0
        end
        x[_pidx(p, g, fi)] = sum(vals) / length(vals)
    end
    return x
end

"""
    parameters(p::GroupPointing, values::Pair...; group = (field = value, ...), ...)
        -> Vector{Float64}

The parameter vector of `p` with named values: each keyword (for `Symbol` group keys,
such as NIF's cones `:inner23` … `:outer50`) or `key => (; field = value, ...)` pair (for
any key, such as a quad name) sets those fields of that group. Fields and groups not
named keep [`initial_parameters`](@ref), so only the values that differ from the
template beams need to be given. An unknown group or a field that is not a parameter of
`p` is an `ArgumentError`.

```julia
pointing = GroupPointing(nif_quads(); fields = (:z, :r, :power))
x0 = parameters(pointing; inner23 = (z = 3.09, r = -1.09), outer50 = (z = 5.01,))
```
"""
parameters(p::GroupPointing, values::Pair...; kwargs...) = _with_named(p, (values..., pairs(kwargs)...))

# Returns the initial parameters of `p` with the named values `kv` set, given as
# `key => (; field = value)` pairs. `_named_values` (parameterizations.jl) resolves the
# names to indices.
function _with_named(p, kv)
    x = initial_parameters(p)
    for (i, v) in _named_values(p, kv)
        x[i] = v
    end
    return x
end

"""
    parameters(p::GroupPointing, beams::AbstractVector{<:Beam}) -> Vector{Float64}

The parameter vector of `p` that reproduces `beams`, such as the result of an
optimization in another parameterization of the same beams. Each template beam of `p`
is matched to the beam of its name in `beams`. An `ArgumentError` is thrown for a
missing name, and extra beams are ignored. The value of a field for each group is the mean
over its members of what `apply` would need for that beam:

- For the aim fields, it is the axial coordinate of the aim (`:z`, mirrored as in
  `apply`) and the components of the aim along the radial and tangential directions of
  the beam (`:r`, `:phi`, in the frame of `unit`). For a `relative` parameterization, it
  is the same components of the shift from the template aim.
- For `:split`, it is the group mean of the signed tangential offset of each beam from
  the mean of its quad, as in `initial_parameters`.
- For `:defocus`, `:power` and `:dlambda`, it is the value of the beam. For `:spot`, it
  is the ratio of the spot semi-axis `a` of the beam to that of the template.

Spherical fields are not supported and throw an `ArgumentError`.

`apply(p, parameters(p, beams))` reproduces the beams exactly when `p` can express
them, as when the groups of `p` split those of the parameterization that produced
`beams` and share its aim frame: a re-pointing by quad handed to a per-beam polish on
the same template,

```julia
quad = GroupPointing(beams; by = :quad, fields = (:z, :phi), relative = true)
s = remake(scenario; parameterization = quad)
r = optimize_pointing(s, objective; ...)
per_beam = GroupPointing(beams; by = :beam, fields = (:z, :phi), relative = true)
x0 = parameters(per_beam, r.beams)       # the box `bounds(per_beam)` stays about `beams`
polish = optimize_pointing(remake(s; parameterization = per_beam), objective; x0, ...)
```

Otherwise the result is the group means of those values: per-beam aims read by quad
give each quad its mean shift.
"""
function parameters(p::GroupPointing, beams::AbstractVector{<:Beam})
    _spherical(p) && throw(ArgumentError("parameters from beams: spherical fields are not supported"))
    byname = Dict(b.name => b for b in beams)
    for t in p.template
        haskey(byname, t.name) || throw(ArgumentError("no beam named $(t.name) among the beams"))
    end
    # These are the tangential offsets of the beam aims, or of their shifts, from the mean
    # aims of their quads.
    d = :split in p.fields ?
        _split_offsets(p, [Vec3(ForwardDiff.value.(byname[t.name].aim)) - (p.relative ? t.aim : zero(Vec3))
                           for t in p.template]) : Float64[]
    x = zeros(nparams(p))
    for (g, mem) in enumerate(p.members), (fi, f) in enumerate(p.fields)
        sg = isempty(d) ? 0.0 : _group_split(p, mem, d)
        if f === :split
            x[_pidx(p, g, fi)] = sg
            continue
        end
        s = 0.0
        for bi in mem
            t = p.template[bi]
            b = _value_beam(byname[t.name])
            er, eφ = p.radial[bi]
            rel = p.relative && f in (:z, :r, :phi)
            a = rel ? b.aim - t.aim : b.aim
            s += f === :z ? _zsign(p, t) * a[3] :
                 f === :r ? dot(a, er) :
                 f === :phi ? (isempty(d) ? dot(a, eφ) : dot(a, eφ) - p.split[bi][1] * sg) :
                 f === :defocus ? b.defocus :
                 f === :power ? b.power :
                 f === :dlambda ? b.dlambda : b.spot.a / t.spot.a
        end
        x[_pidx(p, g, fi)] = s / length(mem)
    end
    return x
end

"""
    apply(p::GroupPointing, x) -> Vector{Beam}

Beams for parameter vector `x` (element type follows `x`, so this is differentiable).
"""
function apply(p::GroupPointing, x::AbstractVector{T}) where {T}
    length(x) == nparams(p) || throw(DimensionMismatch("expected $(nparams(p)) parameters"))
    out = Vector{Beam{T}}(undef, length(p.template))
    for (g, mem) in enumerate(p.members)
        val(f) = (fi = findfirst(==(f), p.fields); fi === nothing ? nothing : x[_pidx(p, g, fi)])
        z, r, ϕ, dfc, pw, sp = val(:z), val(:r), val(:phi), val(:defocus), val(:power), val(:spot)
        dl, spl = val(:dlambda), val(:split)
        # The spherical values are held as concrete types: a flag for whether each is a
        # parameter, and its value.
        sv = (val(:theta_s), val(:dphi_s), val(:split_theta), val(:split_phi))
        has = map(!isnothing, sv)
        vals = map(v -> v === nothing ? zero(T) : v, sv)
        for bi in mem
            b = p.template[bi]
            if _spherical(p)
                aim = _sphere_aim(p, g, bi, has, vals)
            else
                er, eφ = p.radial[bi]
                aim0 = b.aim
                if p.relative
                    rr = r === nothing ? dot(aim0, er) : dot(aim0, er) + r
                    tt = ϕ === nothing ? dot(aim0, eφ) : dot(aim0, eφ) + ϕ
                    zz = z === nothing ? aim0[3] : aim0[3] + _zsign(p, b) * z
                    spl === nothing || (tt += p.split[bi][1] * spl)
                else
                    rr = r === nothing ? dot(aim0, er) : r
                    tt = ϕ === nothing ? dot(aim0, eφ) : ϕ
                    zz = z === nothing ? aim0[3] : _zsign(p, b) * z
                    if spl !== nothing
                        # The split sets the offset of the beam from the mean of its
                        # quad. Where the aim is the template aim, the template offset
                        # is taken out first.
                        side, d0 = p.split[bi]
                        tt = (ϕ === nothing ? tt - d0 : tt) + side * spl
                    end
                end
                aim = rr * er + tt * eφ + SVector{3,T}(0, 0, zz)
            end
            out[bi] = _pointed(b, aim, T, sp, pw, dfc, dl)
        end
    end
    return out
end

# Returns template beam `b` in element type `T` with aim `aim` and the non-aim fields of a
# parameterization: spot scale `sp`, power `pw`, defocus `dfc` and detuning `dl`. For a
# field given as `nothing`, the template value is used. `GroupPointing` and
# `FramePointing` both call this function.
@inline function _pointed(b::Beam, aim, ::Type{T}, sp, pw, dfc, dl) where {T}
    spot = sp === nothing ? _retype(b.spot, T) :
           _retype(with_axes(b.spot, sp * b.spot.a, sp * b.spot.b), T)
    return remake(b; aim, spot,
                  power = pw === nothing ? T(b.power) : pw,
                  defocus = dfc === nothing ? T(b.defocus) : dfc,
                  dlambda = dl === nothing ? T(b.dlambda) : dl)
end

# Returns the surface aim of template beam `bi` in group `g`. `has` flags which of
# `:theta_s`, `:dphi_s`, `:split_theta` and `:split_phi` are parameters, and `v` contains
# their values. The polar angle in the mirrored frame of the beam is a group value plus the
# beam split, and so is the azimuth offset from the beam port. When neither of the two
# fields for an angle is a parameter, the template value of the angle is used. With
# `relative`, the values are shifts of the template value.
function _sphere_aim(p::GroupPointing, g, bi, has::NTuple{4,Bool}, v::NTuple{4,T}) where {T}
    θb, Δφb, sθ, sφ, φp = p.sph[bi]
    Θ0, Φ0, Sθ0, Sφ0 = p.sph0[g]
    if p.relative
        θ = θb + v[1] + sθ * v[3]
        Δφ = Δφb + v[2] + sφ * v[4]
    else
        θ = !(has[1] || has[3]) ? T(θb) :
            (has[1] ? v[1] : T(Θ0)) + sθ * (has[3] ? v[3] : T(Sθ0))
        Δφ = !(has[2] || has[4]) ? T(Δφb) :
             (has[2] ? v[2] : T(Φ0)) + sφ * (has[4] ? v[4] : T(Sφ0))
    end
    θ = _zsign(p, p.template[bi]) < 0 ? 180 - θ : θ
    φ = φp + Δφ
    R, c = p.radius, p.center
    return SVector{3,T}(c[1] + R * sind(θ) * cosd(φ), c[2] + R * sind(θ) * sind(φ),
                        c[3] + R * cosd(θ))
end

"""
    bounds(p::GroupPointing, x0; z = 1.0, r = 1.0, phi = 1.0, split = 1.0, defocus = 2.0,
           power = 0.5, spot = 0.5, dlambda = 2.0, theta_s = 5.0, dphi_s = 10.0,
           split_theta = 5.0, split_phi = 5.0, power_min = 0.0, power_max = Inf)
           -> (lower, upper)

Box bounds around `x0`: `±z`, `±r`, `±phi`, `±split`, `±defocus` mm, relative `±power`, `±spot`
fractions, `±dlambda` Å, and `±theta_s`, `±dphi_s`, `±split_theta`, `±split_phi`
degrees (engineering defaults that are not facility limits). Power bounds are further
clipped to the absolute limits `power_min ≤ P ≤ power_max` (TW per beam, e.g. the per-quad
capability of a facility). A group at zero power in `x0` (a dropped quad) gets `±power` of
the mean of the groups' nonzero powers, `[0, power · P̄]` before the clipping, so that it
can come back; with every group at zero the upper bound is `power_max`, which must then be
finite.
"""
bounds(p::GroupPointing, x0 = initial_parameters(p); z = 1.0, r = 1.0, phi = 1.0,
       split = 1.0, defocus = 2.0, power = 0.5, spot = 0.5, dlambda = 2.0, theta_s = 5.0,
       dphi_s = 10.0, split_theta = 5.0, split_phi = 5.0, power_min = 0.0,
       power_max = Inf) =
    _box_bounds(p, x0, (; z, r, phi, split, defocus, dlambda, theta_s, dphi_s, split_theta, split_phi),
                power, spot, power_min, power_max)

# Returns the box of `bounds` for a parameterization with groups `p.members`, fields
# `p.fields` and `_pidx`. The half-widths `widths` are additive and given by field name;
# `power` and `spot` are relative.
function _box_bounds(p, x0::AbstractVector, widths::NamedTuple, power, spot, power_min, power_max)
    lo, hi = copy(x0), copy(x0)
    fp = findfirst(==(:power), p.fields)
    # A group at zero power is scaled by the mean of the nonzero powers.
    P = fp === nothing ? Float64[] : [x0[_pidx(p, g, fp)] for g in eachindex(p.members)]
    Pbar = any(>(0), P) ? sum(filter(>(0), P)) / count(>(0), P) : NaN
    for g in eachindex(p.members), (fi, f) in enumerate(p.fields)
        i = _pidx(p, g, fi)
        if f === :power
            if x0[i] > 0
                lo[i] = max(x0[i] * (1 - power), power_min)
                hi[i] = min(x0[i] * (1 + power), power_max)
            else
                isnan(Pbar) && !isfinite(power_max) &&
                    throw(ArgumentError("every group is at zero power: give a finite power_max"))
                lo[i] = max(x0[i], power_min)
                hi[i] = isnan(Pbar) ? power_max : min(x0[i] + power * Pbar, power_max)
            end
            lo[i] <= hi[i] || throw(ArgumentError("empty power bounds for group $(p.keys[g])"))
        elseif f === :spot
            lo[i] = x0[i] * (1 - spot); hi[i] = x0[i] * (1 + spot)
        else
            haskey(widths, f) || throw(ArgumentError("no bound width for the field $f"))
            δ = widths[f]
            lo[i] = x0[i] - δ; hi[i] = x0[i] + δ
        end
    end
    return lo, hi
end
