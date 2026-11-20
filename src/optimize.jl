# Pointing parameterizations, objectives, and gradient-based optimization.

const _CACHED_LMAX = 12

const _FIELDS = (:z, :r, :phi, :split, :defocus, :power, :spot, :dlambda,
                 :theta_s, :dphi_s, :split_theta, :split_phi)
const _CYL_FIELDS = (:z, :r, :phi, :split)
const _SPH_FIELDS = (:theta_s, :dphi_s, :split_theta, :split_phi)

"""
    GroupPointing(beams; by = :group, fields = (:z, :r, :defocus), mirror = true,
                  relative = false, unit = :quad, radius = nothing, center = (0, 0, 0))

Maps a parameter vector to beams. Beams are grouped by `by` (`:group` — e.g. cone —,
`:quad`, `:beam`, or a function `beam -> key`); each group gets one value per field:

- `:z` — axial aim coordinate (mm). With `mirror`, beams from the lower hemisphere
  (port `z < 0`) aim at `-z`, so one value describes a top/bottom-symmetric cone.
- `:r` — radial aim offset (mm) toward the port azimuth of the beam's quad (its own with
  `unit = :beam`, below).
- `:phi` — tangential aim offset (mm), counterclockwise about `z`.
- `:split` — tangential split of a quad's beams about its aim point (mm), along the
  quad's φ̂: each beam's aim moves by `±split`, with the sign of its own port's azimuth
  from its quad's mean port (beams counterclockwise of it `+split`), so a positive split
  opens a quad's aims in the sense of its ports and a negative one crosses them (the `s`
  of [`hybride_mock_pointing`](@ref) is `split = -s`). It sets each beam's tangential
  offset from its quad's mean template aim (a shift of it with `relative`); a quad held
  as one beam gets none. Needs `unit = :quad` and `quad` names.
- `:defocus` — focal-plane offset along the beam (mm, see `Beam`).
- `:power` — power per beam (TW).
- `:spot` — multiplier on the best-focus spot semi-axes.
- `:dlambda` — wavelength detuning (Å), used by laser–plasma models.

Fields not listed keep each beam's template value.

**Pointing on a sphere** (polar direct drive). With `radius = R` (mm) the aim can instead
be a point on the sphere of radius `R` about `center`, the target surface on which
direct-drive pointings are published (Barlow et al., PRL 133, 175101, 2024, Table I;
Yeamans et al., Nucl. Fusion 61, 046031, 2021, Table 1; Hohenberger et al., PoP 22,
056308, 2015, Table I):

- `:theta_s` — polar angle of the aim point (degrees from `+z`). With `mirror`, beams
  from the lower hemisphere aim at `180° - theta_s`, so one value describes a
  top/bottom-symmetric cone pair and a larger value moves both toward the equator.
- `:dphi_s` — azimuth of the aim point relative to the beam's own port azimuth
  (degrees, counterclockwise about `z`; the same for both hemispheres).
- `:split_theta`, `:split_phi` — split of a quad's beams (degrees): each beam's aim
  moves by `±split_theta` in polar angle and `±split_phi` in azimuth, with the sign of
  its own port's offset from its quad's mean port (beams at a larger polar angle than
  their quad's, counted from the beam's own pole when mirrored, take `+split_theta`;
  beams counterclockwise of it `+split_phi`), so positive values open a quad's aims in
  the sense of its ports. This is the split-quad pointing of Cok et al. (PoP 15,
  082705, 2008) and Yeamans et al. 2021 (left/right `Δφ`) as angles on the sphere.
  A quad given as one beam, or a beam level with its quad's mean, gets no split. Needs
  `quad` names; other per-beam aims come from `by` (e.g. one group per sub-ring).

The aim is the surface point itself, `aim = center + R * port_vector(θ, φ)`, which the
axis meets first when it faces the port (`Beam`). Best focus is then at
`aim + defocus * k`: defocus is measured from the surface aim point, as VisRad measures
the focus offset from its reference point (research notes, "Pointing and focus
conventions"). Each template beam's angles are those of the first point where its axis
meets the sphere (an `ArgumentError` if it misses); listed spherical values are shared
by a group's members, and an angle that no listed field governs keeps each beam's own
value. The aim moves to the surface along the template's axis, so the focus moves with
it unless the template aim was already on the sphere. The spherical fields and the
hohlraum fields `:z`, `:r`, `:phi` exclude each other.

A quad is the unit of the hohlraum fields' pointing by default (`unit = :quad`): `:r`,
`:phi` and `:split` are measured in the azimuth of the quad's port, the mean port of the
beams that share its `quad` name, instead of each beam's own. The beams of a quad that
share their values then share one aim point, whether the table holds the quad as one
beam or as its four, so a NIF pointing given per cone or per quad stays one aim per quad
while the deposition is computed per beam. A quad's beams pointed apart about its aim
point:

```julia
p = GroupPointing(nif_beams(); by = :quad, fields = (:z, :r, :phi, :split))
beams = apply(p, parameters(p, "Q31T" => (z = 5.8, split = 0.5)))
```

`by = :beam` gives each beam its own aim (its offsets still in its quad's azimuth), and
a per-beam table (an as-shot SST export, `read_visrad_sst(path; template =
nif_beams(quads))`) keeps its beams' differences under one shift per quad with
`relative = true` (below). `unit = :beam` measures each beam's `:r` and `:phi` in its own
port azimuth instead, so a quad's four beams with the same `r` aim `r sin Δφ` apart (Δφ:
a beam's azimuth from its quad's, up to 5.7° on NIF's 23.5° cone). A quad held as one
beam, and a beam without a quad name, are the same either way, bit for bit.

`relative = true` makes the aim fields (`:z`, `:r`, `:phi`, `:split`, or the spherical ones) shifts
added to each beam's own template aim instead of values the group's beams share, so a
group moves rigidly and keeps its members' differences: re-pointing an as-shot table,
whose beams within a quad aim at different points (a facility's per-beam offsets), by
one shift per quad. Its initial parameters are then zero; the other fields stay
absolute.
"""
struct GroupPointing{K}
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
    unit::Symbol                         # :beam or :quad: whose port azimuth :r and :phi refer to
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

# A `GroupPointing` of `beams` with group key `ks[i]` for beam i (the constructor's
# keywords but `by`).
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

# Per template beam, the side of its quad it splits to: the sign of its port's azimuth
# from the quad's mean port (0 for a beam level with it, and for a quad held as one
# beam), and its aim's tangential offset (along the quad's φ̂, `radial`'s second vector)
# from the mean aim of its quad's beams.
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

# The tangential offsets of `aims` (one per template beam, or their shifts from the
# template's) from the mean of their quad's, along each quad's φ̂.
function _split_offsets(p::GroupPointing, aims::AbstractVector{Vec3})
    asum, n = Dict{String,Vec3}(), Dict{String,Int}()
    for (b, a) in zip(p.template, aims)
        asum[b.quad] = get(asum, b.quad, zero(Vec3)) + a
        n[b.quad] = get(n, b.quad, 0) + 1
    end
    return [dot(a - asum[b.quad] / n[b.quad], p.radial[i][2]) for (i, (b, a)) in enumerate(zip(p.template, aims))]
end

# A group's split value from its members' tangential offsets `d` from their quads' mean
# aims: Σ side × d over the members with a side, over their number (0 without one).
function _group_split(p::GroupPointing, mem, d)
    k = count(bi -> p.split[bi][1] != 0, mem)
    k == 0 && return 0.0
    return sum(p.split[bi][1] * d[bi] for bi in mem if p.split[bi][1] != 0) / k
end

# Per template beam, the radial and tangential directions of its quad's port: of the mean
# port of the beams that share its quad name (a quad held as one beam: its own port,
# exactly); a beam without a quad name keeps its own.
function _quad_radial(tmpl)
    qsum = Dict{String,Vec3}()
    for b in tmpl
        isempty(b.quad) || (qsum[b.quad] = get(qsum, b.quad, zero(Vec3)) + b.port)
    end
    return [isempty(b.quad) ? _radial(b) : _radial(qsum[b.quad]) for b in tmpl]
end

# Per template beam: the polar angle θ' (degrees; from the beam's own pole when mirrored)
# and the azimuth offset Δφ from its port of the first point where its axis meets the
# sphere, its split signs (sθ, sφ) relative to its quad's mean port, and its port azimuth.
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

# A group's initial spherical values: mean θ' and Δφ, and the mean split
# Σ s (x - x̄) / #(s ≠ 0) (zero without split members).
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
        # the group's split: each beam's offset from its quad's mean, by its side
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

# The initial parameters of `p` with the named values `kv` (`key => (; field = value)`
# pairs) set; `_named_values` (parameterizations.jl) resolves them to indices.
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
optimization in another parameterization of the same beams: each template beam of `p`
takes the beam of its name in `beams` (an `ArgumentError` when there is none; extra
beams are ignored), and each group's value of a field is the mean over its members of
what `apply` would need for that beam. For the aim fields that is the aim's axial
coordinate (`:z`, mirrored as in `apply`) and its components along the beam's radial
and tangential directions (`:r`, `:phi`, in the frame of `unit`), or, for a `relative`
parameterization, the same components of the shift from the template beam's aim
(`:split`: the group's mean of each beam's signed tangential offset from its quad's
mean, as in `initial_parameters`);
`:defocus`, `:power` and `:dlambda` are the beams' values and `:spot` the ratio of
the beam's spot semi-axis `a` to the template's. Spherical fields are not supported (an
`ArgumentError`).

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
    # the beams' tangential offsets from their quads' mean aims (or of their shifts)
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
        # spherical values as concrete types (whether each is a parameter, and its value)
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
                        # the split sets the beam's offset from its quad's mean: the
                        # template's own offset out where the aim is the template's
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

# Template beam `b` with aim `aim` and the parameterizations' non-aim fields: spot scale
# `sp`, power `pw`, defocus `dfc` and detuning `dl` (`nothing`: the template's value),
# in element type `T`. Shared by `GroupPointing` and `FramePointing`.
@inline function _pointed(b::Beam, aim, ::Type{T}, sp, pw, dfc, dl) where {T}
    spot = sp === nothing ? _retype(b.spot, T) :
           _retype(with_axes(b.spot, sp * b.spot.a, sp * b.spot.b), T)
    return remake(b; aim, spot,
                  power = pw === nothing ? T(b.power) : pw,
                  defocus = dfc === nothing ? T(b.defocus) : dfc,
                  dlambda = dl === nothing ? T(b.dlambda) : dl)
end

# The surface aim of template beam `bi` in group `g`: `has` tells which of
# (:theta_s, :dphi_s, :split_theta, :split_phi) are parameters and `v` holds their values.
# The polar angle in the beam's mirrored frame and the azimuth offset from its port are
# each a group value plus the beam's split, or the template beam's own when neither of
# the two is a parameter (`relative`: shifts of the template beam's own).
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
degrees (engineering defaults, not facility limits). Power bounds are further clipped to
the absolute limits `power_min ≤ P ≤ power_max` (TW per beam, e.g. a facility's per-quad
capability). A group at zero power in `x0` (a dropped quad) gets `±power` of the mean
of the groups' nonzero powers, `[0, power · P̄]` before the clipping, so that it can
come back; with every group at zero the upper bound is `power_max`, which must then be
finite.
"""
bounds(p::GroupPointing, x0 = initial_parameters(p); z = 1.0, r = 1.0, phi = 1.0,
       split = 1.0, defocus = 2.0, power = 0.5, spot = 0.5, dlambda = 2.0, theta_s = 5.0,
       dphi_s = 10.0, split_theta = 5.0, split_phi = 5.0, power_min = 0.0,
       power_max = Inf) =
    _box_bounds(p, x0, (; z, r, phi, split, defocus, dlambda, theta_s, dphi_s, split_theta, split_phi),
                power, spot, power_min, power_max)

# The box of `bounds` for a parameterization with groups `p.members`, fields `p.fields`
# and `_pidx`: additive half-widths `widths` by field name; `power` and `spot` relative.
function _box_bounds(p, x0::AbstractVector, widths::NamedTuple, power, spot, power_min, power_max)
    lo, hi = copy(x0), copy(x0)
    fp = findfirst(==(:power), p.fields)
    # the scale of a group at zero power: the mean of the nonzero ones
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

# ---------------------------------------------------------------------------------------
# Scenario: precomputed operators + evaluation
# ---------------------------------------------------------------------------------------

# The full flux response of a scenario, computed when first read.
mutable struct _LazyResponse
    const lock::ReentrantLock
    R::Union{Nothing,Matrix{Float64}}
end
_LazyResponse() = _LazyResponse(ReentrantLock(), nothing)

# What `optimize_pointing` runs on: a `Scenario`, or a `PhasedScenario` (phased.jl).
abstract type _AbstractScenario end

"""
    Scenario(target, pointing::GroupPointing; observe = nothing,
             radiation = observe !== nothing, time = 0.0, albedo = nothing, lpi = nothing,
             deposition_kw = (;), viewfactor_kw = (;), model = nothing)
    Scenario(target, beams; kwargs...)
    Scenario(target, parameterization; kwargs...)

The parameterization (`scenario.param`) maps the parameter vectors of `evaluate` and
`optimize_pointing` to beams. The first form builds the deposition operator for the
template beams of `pointing` (a [`FramePointing`](@ref) too); the second is `Scenario(target, GroupPointing(beams))`
(fields `:z`, `:r`, `:defocus` per group). The third takes any other parameterization
with methods of `apply(parameterization, x)` and `initial_parameters(parameterization)`,
and builds the operator for the beams at its initial parameters. The operator resolves
the beams' footprints where those beams aim (`reach`, `DepositionOperator`).

Everything precomputed for fast repeated evaluation of a pointing design:
the `DepositionOperator`, and — when `radiation` and `observe` (component names whose
incident x-ray flux is the figure of merit, e.g. `:capsule`) are given — the view
factors and the linear flux response of the observed elements. `lpi` is a laser–plasma
model (`Backscatter`, `CalibratedTransfer`, `CBET`, or a tuple of them) applied at `time`
between the setup beams and deposition.

The flux response is linear in the laser source at fixed albedos: the materials' at
`time` (temperature-dependent ones at their reference temperatures, thin foils
intact), or an explicit per-element `albedo`, such as the converged `sol.albedo` of a
solution at the design point (e.g. a drive history's peak, windows burned through).
For a time partway through a shaped pulse (a phase), the materials' albedos at `time`
are not the wall's: a `WallLossAlbedo`'s follows each element's temperature history,
which neither the reference temperature nor a standalone `solve` at `time` knows. Pass
the history's own, `st = history_state(target, history, τ)` as `st.time, st.albedo`
([`history_state`](@ref)). A scenario with a response (`observe` and a model) built
without `albedo` on a target whose materials include a [`WallLossAlbedo`](@ref) (or a
[`BurnthroughAlbedo`](@ref) whose `before` is one) warns at construction: its loss is
taken at the reference temperature (250 eV by default), which mid-pulse overstated the
drive against a drive history by up to 1.64 times on a public mock-up (research notes,
"Drive histories: sliver steps and the albedo of a phase"). `albedo = :reference`
marks the reference temperature as intended: the same scenario, without the warning.

Objectives need only the moments of the observed flux, so the scenario keeps the
response contracted with the Legendre (`lmax = 12`) and real-harmonic (`lmax = 6`)
weights: per evaluation a 13 × n and a 49 × n product, and a 62-column solve at setup.
The full observed × n response `scenario.response` (and an evaluation's `flux`) is
computed when first read; on a large (iterative) model an evaluation's `flux` is one
solve instead, until the response has been read.

The operators' own options pass through two named tuples:
- `deposition_kw` to [`DepositionOperator`](@ref): `profile` (a
  [`VisRadProfile`](@ref): VisRad's beam model instead of the f-cone superposition),
  `ncone` (f-cone directions, 13), `reach` (the aim motion, in mm, that the operator
  stays accurate for, 2.0), `resolution`, `level`, `maxlevel` and `shadow` (the
  quadrature's refinement);
- `viewfactor_kw` to [`RadiationModel`](@ref) and through it to
  [`view_factors`](@ref): `closed` (the closure correction; by default the target's own
  declaration, which `hohlraum` makes), `F` (given view factors, e.g. those of a target
  that differs only in its materials), `storage` (`Float64`, `Float32` or `:hmatrix`),
  `solver`, `visibility` (`:centroid` or `:sampled`), `visibility_grid`, `level` and
  `near`.

`model`, a `RadiationModel` of the target, is used instead of building one; it excludes
`viewfactor_kw` (an `ArgumentError`).

[`remake`](@ref)`(scenario; ...)` replaces the laser–plasma model, the parameterization,
the deposition options or the response's observed elements, time and albedos, and
rebuilds only what depends on them.
"""
struct Scenario{P} <: _AbstractScenario
    target::Target
    param::P
    deposition::DepositionOperator
    model::Union{Nothing,RadiationModel}
    observed::Vector{Int}
    moment_weights::Matrix{Float64}   # legendre_weights of the observed elements (lmax = 12)
    ylm_weights::Matrix{Float64}      # ylm_weights of the observed elements (lmax = 6)
    moment_response::Matrix{Float64}  # moment_weights * response: moments per element power
    ylm_response::Matrix{Float64}     # ylm_weights * response
    lpi::Any                          # laser–plasma model (see apply_lpi)
    time::Float64
    _response::_LazyResponse          # the full response, once read (shared by subsets)
    albedo::Union{Nothing,Vector{Float64}}   # explicit albedos of the response
    _split::Vector{Vector{Symbol}}    # beam-group sets whose element powers `evaluate` keeps
end

# A scenario whose evaluations keep no group sets' powers (`_with_split` sets them).
Scenario(t::Target, param, op::DepositionOperator, model, obs, W, Y, WR, YR, lpi, time, response::_LazyResponse,
         albedo) = Scenario(t, param, op, model, obs, W, Y, WR, YR, lpi, time, response, albedo, Vector{Symbol}[])

const _CACHED_YLM_LMAX = 6

Base.getproperty(s::Scenario, name::Symbol) =
    name === :response ? _response(s) : getfield(s, name)
Base.propertynames(::Scenario) =
    (:target, :param, :deposition, :model, :observed, :response, :moment_weights,
     :ylm_weights, :moment_response, :ylm_response, :lpi, :time, :albedo)

# Whether the scenario has a radiation response (observed elements and a model).
_radiative(s::Scenario) = s.model !== nothing && !isempty(s.observed)

function _response(s::Scenario)
    _radiative(s) || return nothing
    c, m = getfield(s, :_response), s.model::RadiationModel
    lock(c.lock) do
        c.R === nothing && (c.R = flux_response(m, s.observed; time = s.time, s.albedo))
        c.R
    end
end

Scenario(t::Target, p::GroupPointing; kwargs...) = _scenario(t, p.template, p; kwargs...)
Scenario(t::Target, beams::AbstractVector{<:Beam}; kwargs...) =
    _scenario(t, beams, GroupPointing(beams); kwargs...)
Scenario(t::Target, p; kwargs...) = _scenario(t, apply(p, initial_parameters(p)), p; kwargs...)

function _scenario(t::Target, beams::AbstractVector{<:Beam}, parameterization;
                   observe = nothing, radiation::Bool = observe !== nothing, time = 0.0,
                   albedo = nothing, deposition_kw = (;), viewfactor_kw = (;), lpi = nothing,
                   model = nothing)
    model === nothing || model.target === t ||
        throw(ArgumentError("model is a radiation model of another target"))
    model === nothing || isempty(viewfactor_kw) ||
        throw(ArgumentError("give a model or viewfactor_kw, not both (the model's view factors are fixed)"))
    # `:reference`: the materials' own albedos at `time`, stated as intended
    meant = albedo === :reference
    albedo isa Symbol && !meant &&
        throw(ArgumentError("albedo takes albedos per element, `nothing` or `:reference`, got :$albedo"))
    meant && (albedo = nothing)
    op = DepositionOperator(t, beams; deposition_kw...)
    obs = _observed(t, observe)
    model = !radiation ? nothing : model === nothing ? RadiationModel(t; viewfactor_kw...) : model
    meant || albedo !== nothing || model === nothing || isempty(obs) || _warn_reference_albedo(t, time)
    W, Y, WR, YR = _moment_responses(t, model, obs, time, albedo)
    return Scenario(t, parameterization, op, model, obs, W, Y, WR, YR, lpi, Float64(time),
                    _LazyResponse(), albedo === nothing ? nothing : collect(Float64, albedo))
end

# A scenario's response without albedos on a target whose wall loss follows a
# temperature history in a drive history (a `WallLossAlbedo`, or a thin foil's
# `BurnthroughAlbedo` over one): the loss is taken at the model's reference temperature,
# which mid-pulse overstates the drive against the history (research notes, "Drive
# histories: sliver steps and the albedo of a phase"). One warning per construction.
_history_albedo(a) = a isa WallLossAlbedo || a isa BurnthroughAlbedo && a.before isa WallLossAlbedo

function _warn_reference_albedo(t::Target, time)
    i = findfirst(c -> c.kind === :surface && _history_albedo(c.material.albedo), t.components)
    i === nothing && return nothing
    c = t.components[i]
    w = c.material.albedo isa BurnthroughAlbedo ? c.material.albedo.before : c.material.albedo
    @warn "Scenario: $(c.part)'s WallLossAlbedo is evaluated at its reference temperature " *
          "($(w.reference_temperature) eV) at time = $(Float64(time)) ns, as a wall held at it since the " *
          "pulse start; mid-pulse that can overstate the drive against a drive_history (by up to 1.64 times " *
          "on a public mock-up). For a phase of a pulse pass the history's albedos, " *
          "`st = history_state(target, history, time)` with `time = st.time, albedo = st.albedo` " *
          "(or `remake(scenario; history, time)`); `albedo = :reference` marks the reference " *
          "temperature is intended."
    return nothing
end

# The observed elements of `observe` (component or part names; `nothing`: none).
_observed(t::Target, observe) =
    observe === nothing ? Int[] : elements(t, (observe isa Symbol ? (observe,) : observe)...)

# The mode weights of the observed elements and, with a model, their responses: one
# flux-response solve for all weights.
function _moment_responses(t::Target, model, obs, time, albedo)
    W = isempty(obs) ? zeros(0, 0) : legendre_weights(t, obs; lmax = _CACHED_LMAX)
    Y = isempty(obs) ? zeros(0, 0) : ylm_weights(t, obs; lmax = _CACHED_YLM_LMAX)
    WR, YR = zeros(0, 0), zeros(0, 0)
    if model !== nothing && !isempty(obs)
        WY = flux_response(model, obs; time, albedo, weights = vcat(W, Y))
        WR, YR = WY[1:size(W, 1), :], WY[size(W, 1)+1:end, :]
    end
    return W, Y, WR, YR
end

# The beams a scenario's deposition operator is built for, from its parameterization.
_scenario_beams(p::GroupPointing) = p.template
_scenario_beams(p) = apply(p, initial_parameters(p))

# Marks a keyword of `remake(::Scenario)` that was not given.
struct _Unchanged end
const _UNCHANGED = _Unchanged()

"""
    remake(scenario; parameterization, lpi, deposition_kw, observe, time, albedo, model,
           history) -> Scenario

Copy of `scenario` with the given parts replaced. What depends on a replaced part is
rebuilt and everything else is shared:

- `lpi` (a laser–plasma model or `nothing`): nothing is rebuilt.
- `parameterization`: the deposition operator is kept when the new parameterization's
  beams (a `GroupPointing`'s or `FramePointing`'s template, another's beams at its
  initial parameters) enter through the operator's ports in the same order with the same
  f-cone directions (and, under a [`VisRadProfile`](@ref), the same envelope ratios);
  otherwise the change is refused with an `ArgumentError`, unless `deposition_kw` is
  given too. A kept operator resolves the footprints where the original beams aimed
  (its `reach`), as for [`subset`](@ref).
- `deposition_kw`: the operator is rebuilt for the parameterization's beams with these
  options. Options not given take their defaults: a scenario does not record the
  options its operator was built with, so pass all that differ from the defaults (a
  `profile`, `reach`, `resolution`).
- `observe`, `time`, `albedo` and `model` (a `RadiationModel` of the same target, or
  `nothing` for none): the mode weights and moment responses are recomputed (one
  flux-response solve), and the full `response` is computed again when read. `albedo =
  nothing` or `:reference` takes the materials' own albedos at `time` (no warning:
  the warning, if any, came at the scenario's construction).

The target, `radiation` and `viewfactor_kw` do not change through `remake` (an
`ArgumentError`): they need a new `Scenario`, or a `model` built for the purpose.

`history` (a [`drive_history`](@ref) of the scenario's target run with `keep_albedo =
true`) with `time` sets the albedos and time of the history's state at `time`:
`remake(scenario; history = h, time = τ)` is `st = `[`history_state`](@ref)`(scenario.target,
h, τ)` and `remake(scenario; st.time, st.albedo)`, bit for bit. Remade so, a design's
phase and its reference at the same time carry the same `time`, which the reference
checks of [`mode_gate`](@ref) and [`drive_floor`](@ref) compare.
"""
function remake(s::Scenario; parameterization = _UNCHANGED, lpi = _UNCHANGED,
                deposition_kw = _UNCHANGED, observe = _UNCHANGED, time = _UNCHANGED,
                albedo = _UNCHANGED, model = _UNCHANGED, history = nothing, kwargs...)
    isempty(kwargs) ||
        throw(ArgumentError("remake(::Scenario) cannot change $(join(keys(kwargs), ", ")): build a new Scenario " *
                            "(or pass a `model`); it changes parameterization, lpi, deposition_kw, observe, " *
                            "time, albedo and model"))
    if history !== nothing
        time === _UNCHANGED && throw(ArgumentError("a `history` needs the `time` (ns) of its state"))
        albedo === _UNCHANGED ||
            throw(ArgumentError("give the albedos by `history` or by `albedo`, not both"))
        st = history_state(s.target, history, time)
        time, albedo = st.time, st.albedo
    end
    t = s.target
    param = parameterization === _UNCHANGED ? s.param : parameterization
    op = s.deposition
    if deposition_kw !== _UNCHANGED
        op = DepositionOperator(t, _scenario_beams(param); deposition_kw...)
    elseif parameterization !== _UNCHANGED
        _check_operator(op, _scenario_beams(param))
    end
    lp = lpi === _UNCHANGED ? s.lpi : lpi
    if all(v -> v === _UNCHANGED, (observe, time, albedo, model))
        return Scenario(t, param, op, s.model, s.observed, s.moment_weights, s.ylm_weights,
                        s.moment_response, s.ylm_response, lp, s.time, getfield(s, :_response),
                        s.albedo)
    end
    m = model === _UNCHANGED ? s.model : model
    m === nothing || m.target === t ||
        throw(ArgumentError("model is a radiation model of another target"))
    obs = observe === _UNCHANGED ? s.observed : _observed(t, observe)
    tm = time === _UNCHANGED ? s.time : Float64(time)
    alb = albedo === _UNCHANGED ? s.albedo : albedo === nothing || albedo === :reference ? nothing :
          collect(Float64, albedo)
    W, Y, WR, YR = _moment_responses(t, m, obs, tm, alb)
    return Scenario(t, param, op, m, obs, W, Y, WR, YR, lp, tm, _LazyResponse(), alb)
end

# Whether the operator `op` serves `beams`: their ports and f-cone directions in order,
# and their envelope ratios under VisRad's profile (an `ArgumentError` otherwise).
function _check_operator(op::DepositionOperator, beams)
    fail(why) = throw(ArgumentError("the parameterization's beams $why the deposition operator's: " *
                                    "give deposition_kw to rebuild it"))
    length(beams) == length(op.ports) || fail("are $(length(beams)), not the $(length(op.ports)) of")
    all(i -> beams[i].port == op.ports[i], eachindex(beams)) || fail("enter through other ports than")
    all(i -> cone_directions(beams[i], length(op.dirs[i]))[1] == op.dirs[i], eachindex(beams)) ||
        fail("have other f-cone directions (f-numbers) than")
    op.profile === nothing || [_vr_envelope(op.profile, b) for b in beams] == op.envelope ||
        fail("have other VisRad envelope ratios than")
    return nothing
end

"""
    subset(scenario, keep; parameterization = the scenario's, restricted to the kept beams)

Scenario for the beams `keep` (indices into the scenario's beams) — e.g. after
dropping a bundle — reusing the view factors, flux response and mode weights; only the
deposition operator is restricted. The kept beams' template values (pointing, power)
come from the scenario's parameterization. The default new parameterization keeps the
groups of the scenario's: a [`GroupPointing`](@ref)'s groups that keep a beam keep their
keys, so the parameter names of the subset are those of the full pointing, with its
fields, mirroring, `relative`, `unit` and sphere (with `unit = :quad`, the default, a
quad's azimuth and split are those of its kept beams); pass a parameterization to group
the kept beams otherwise. For a [`FramePointing`](@ref) scenario, the default keeps the
kept beams' frames, groups, fields and `relative`.

For a [`HeldParameters`](@ref) or [`PowerScales`](@ref) scenario the default wraps the
inner parameterization's subset, so the kept parameters keep their names, and is what
the wrappers built by hand on the kept beams give: a parameter whose group keeps some of
its beams keeps its value and stays held or free (the value applies to the kept beams as
it applied to all, an absolute aim or a shift); one whose group keeps none is dropped,
held or free, and so is the power scale of a group with no kept beam. Any other
parameterization needs one.
"""
function subset(s::Scenario, keep::AbstractVector{<:Integer}; parameterization = nothing)
    param = parameterization === nothing ? _subset_param(s.param, keep) : parameterization
    return Scenario(s.target, param, s.deposition[keep], s.model, s.observed,
                    s.moment_weights, s.ylm_weights, s.moment_response, s.ylm_response,
                    s.lpi, s.time, getfield(s, :_response), s.albedo)
end

# The beams `keep` of `p` in their own groups (the groups that keep a beam keep their
# keys, so their parameter names), with its fields, mirroring, `relative`, `unit` and
# sphere.
function _subset_param(p::GroupPointing, keep)
    key = similar(p.keys, length(p.template))
    for (g, mem) in enumerate(p.members), bi in mem
        key[bi] = p.keys[g]
    end
    return _group_pointing(p.template[keep], key[keep]; fields = Tuple(p.fields), mirror = p.mirror,
                           relative = p.relative, unit = p.unit, radius = _spherical(p) ? p.radius : nothing,
                           center = p.center)
end
_subset_param(p, keep) =
    throw(ArgumentError("give a parameterization for a scenario that is not on a GroupPointing or FramePointing"))

_nparams_shown(p::GroupPointing) = nparams(p)
_nparams_shown(p) = "?"

function Base.show(io::IO, s::Scenario)
    print(io, "Scenario(", nelements(s.target), " elements, ", length(s.deposition.ports),
          " beams, ", _nparams_shown(s.param), " parameters, ",
          length(s.observed), " observed elements)")
end

"""
    Evaluation

Result of `evaluate`: the setup `beams`, the `effective` beams after the scenario's
laser–plasma model, their `deposition`, and `flux` — incident x-ray flux on the
scenario's observed elements (TW/mm²; `nothing` without radiation). Deposition
fractions are relative to effective (post-LPI) power. The flux is computed when first
read; `flux_moments` and `flux_ylm` do not need it.
"""
struct Evaluation{T,B,E}
    scenario::Scenario
    beams::Vector{B}
    effective::Vector{E}
    deposition::Deposition{T}
    _flux::Base.RefValue{Union{Nothing,Vector{T}}}   # observed flux, once computed
    _split::Union{Nothing,_GroupSplit{T}}            # the scenario's group sets' powers
    _moments::Base.RefValue{Union{Nothing,Vector{T}}}  # cached Legendre moments, once computed
    _ylm::Base.RefValue{Union{Nothing,Vector{T}}}      # cached Y_lm coefficients, once computed
end

Evaluation(s::Scenario, beams::Vector, effective::Vector, d::Deposition{T},
           flux::Base.RefValue{Union{Nothing,Vector{T}}}, split = nothing) where {T} =
    Evaluation(s, beams, effective, d, flux, split, Ref{Union{Nothing,Vector{T}}}(nothing),
               Ref{Union{Nothing,Vector{T}}}(nothing))

# The scenario's cached moment (`moment_response`) or Y_lm (`ylm_response`) products with
# the deposition, each computed once per evaluation: a gate, a drive floor and a ceiling
# and the objective on one phase read them alike.
function _cached_product(c::Base.RefValue{Union{Nothing,Vector{T}}}, R::Matrix{Float64},
                         P::Vector{T}) where {T}
    v = c[]
    v === nothing || return v
    w = _mul(R, P)
    c[] = w
    return w
end

Base.getproperty(ev::Evaluation, name::Symbol) =
    name === :flux ? _flux(ev) : getfield(ev, name)
Base.propertynames(::Evaluation) = (:scenario, :beams, :effective, :deposition, :flux)

function _flux(ev::Evaluation)
    s = getfield(ev, :scenario)
    _radiative(s) || return nothing
    c = getfield(ev, :_flux)
    c[] === nothing && (c[] = _observed_flux(s, getfield(ev, :deposition).element))
    return c[]
end

# The observed flux for absorbed laser power P: the full response times P, or on a large
# model whose response has not been read, one iterative solve (the response would take
# a solve per observed element).
function _observed_flux(s::Scenario, P)
    m = s.model
    (m.iterative && getfield(s, :_response).R === nothing) || return _mul(s.response, P)
    t = m.target
    α = s.albedo === nothing ? albedos(t, s.time) : s.albedo
    η = _constant_conversion(t) ? conversions(t) : conversions(t, s.time)
    B = _iterative_radiosity(m, α, η .* P ./ t.elem_area, nothing)
    return _mul(m.F, B)[s.observed]
end

"""
    evaluate(scenario, x) / evaluate(scenario, beams) -> Evaluation
"""
evaluate(s::Scenario, x::AbstractVector{<:Real}) = evaluate(s, apply(s.param, x))

evaluate(s::Scenario, beams::AbstractVector{<:Beam}) = _evaluate(s, beams, apply_lpi(s.lpi, beams, s.time))

# `evaluate` of the setup `beams` and their effective beams `eff`.
function _evaluate(s::Scenario, beams::AbstractVector{<:Beam}, eff)
    sets = getfield(s, :_split)
    (isempty(sets) || !_plain_deposit(s.lpi)) && return _evaluation(s, beams, eff, lpi_deposit(s.lpi, s.deposition, eff, s.time), nothing)
    # the element powers of the group sets the run's constraints read (`_with_split`), in
    # the same pass
    T = mapreduce(b -> typeof(b.power), promote_type, eff)
    split = _GroupSplit(sets, eff, nelements(s.target), T)
    return _evaluation(s, beams, eff, _deposit(s.deposition, eff, split), split)
end

_evaluation(s, beams, eff, dep::Deposition{T}, split) where {T} =
    Evaluation(s, beams, eff, dep, Ref{Union{Nothing,Vector{T}}}(nothing), split)

# Whether the laser–plasma model deposits the effective beams as `deposit` does (no model
# that deposits light off the operator, `_deposits`).
_plain_deposit(m) = !_deposits(m)
_plain_deposit(ms::Union{Tuple,AbstractVector}) = !any(_deposits, ms)

# `A * x` for a Float64 matrix and a vector of duals through BLAS: the values and each
# partial are the columns of one matrix product (a generic dual mat-vec runs ~10×
# slower). Other element types use the generic product.
_mul(A::AbstractMatrix, x::AbstractVector) = A * x
function _mul(A::Matrix{Float64}, x::Vector{ForwardDiff.Dual{Tag,Float64,N}}) where {Tag,N}
    X = Matrix{Float64}(undef, length(x), N + 1)
    for (i, d) in enumerate(x)
        X[i, 1] = ForwardDiff.value(d)
        p = ForwardDiff.partials(d)
        for j in 1:N
            X[i, j+1] = p[j]
        end
    end
    Y = A * X
    return [ForwardDiff.Dual{Tag}(Y[i, 1], ntuple(j -> Y[i, j+1], Val(N))) for i in axes(Y, 1)]
end

"Fraction of delivered laser power striking the outside (back side) of any component."
clipped_fraction(ev::Evaluation) = sum(ev.deposition.backside) / ev.deposition.delivered
"Fraction of delivered power that strikes no surface at all."
missed_fraction(ev::Evaluation) = 1 - sum(ev.deposition.per_beam) / ev.deposition.delivered
"Fraction of delivered power absorbed by the named parts or components (see `elements`)."
power_fraction(ev::Evaluation, names::Symbol...) =
    sum(ev.deposition.element[e] for e in elements(ev.scenario.target, names...); init = zero(ev.deposition.delivered)) /
    ev.deposition.delivered

"""
    flux_moments(ev; lmax = 8, axis = ẑ, center = mode_center(target, observed)) -> Vector

Legendre moments `aₗ` of the observed incident flux (index `l+1`), about the observed
capsule's center by default (`mode_center`).
"""
function flux_moments(ev::Evaluation; lmax::Integer = 8, axis = Vec3(0, 0, 1),
                      center = nothing)
    s = ev.scenario
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    if lmax <= _CACHED_LMAX && Vec3(axis...) == Vec3(0, 0, 1) &&
       (center === nothing || Vec3(center...) == mode_center(s.target, s.observed))
        return _cached_product(getfield(ev, :_moments), s.moment_response, ev.deposition.element)[1:lmax+1]
    end
    return legendre_moments(_full(ev), s.target, s.observed; lmax, axis, center)
end

"""
    mean_flux(ev) -> Real

Mean incident x-ray flux on the observed elements (TW/mm²): the solid-angle mean seen
from the mode center, `flux_moments(ev)[1]`. It is the drive level that the `drive`
keyword of [`symmetry_objective`](@ref) and [`asymmetry_objective`](@ref) holds, so
`drive = mean_flux(evaluate(scenario, x0))` keeps the drive of the starting design.
"""
mean_flux(ev::Evaluation) = flux_moments(ev; lmax = 0)[1]

function _full(ev::Evaluation)
    ev.flux === nothing && throw(ArgumentError("scenario has no radiation response"))
    v = zeros(eltype(ev.flux), nelements(ev.scenario.target))
    v[ev.scenario.observed] .= ev.flux
    return v
end

"""
    symmetry_objective(; modes = (2 => 1.0, 4 => 1.0), odd = 1.0, lmax = 6,
                       clip = 10.0, missed = 10.0, direct = 10.0, direct_names = nothing,
                       drive = nothing, drive_weight = 1.0, center = nothing, keepout = nothing,
                       uncertainty = nothing)

Objective for indirect-drive symmetry: `Σ wₗ (aₗ/a₀)²` over `modes`, plus `odd ×` the
sum of squared odd moments (top/bottom imbalance), plus `clip`, `missed` and `direct`
times the squared fractions of laser power clipped on the outside of the target,
missing it, or hitting the observed elements directly (the scenario's `observe`, the
capsule whose drive is measured), or with `direct_names` the named parts or components
instead (a name that is not in the target is an error). When `drive`
(TW/mm²) is given, adds `drive_weight × (a₀/drive - 1)²`, with `a₀` the
[`mean_flux`](@ref): `drive = mean_flux(evaluate(scenario, x0))` holds the drive of the
design `x0`. A `keepout` (an
`ApertureKeepOut` for the entrance-hole rims, a `SphereKeepOut` for the capsule, or a
tuple of them) adds its penalty on the beams' worst-case clearance. The `direct`
penalty only sees power already on the capsule; keep a margin with a `SphereKeepOut`.

With `uncertainty = (; pointing, power)` (rms beam-to-beam pointing in mm and relative
power, as in `mode_uncertainty`) the mode terms are their expected values over those
errors, linearized: each `(aₗ/a₀)²` gains the variance the errors give it, so the
optimum also makes the drive shape insensitive to them (robust design). This costs
about three gradients per evaluation. The shared errors of `mode_uncertainty`
(`group_pointing`, `group_power`, `groups`, `positioning`) are accepted too; they add
a fourth per-beam sensitivity (about four gradients per evaluation).
"""
function symmetry_objective(; modes = (2 => 1.0, 4 => 1.0), odd = 1.0, lmax::Integer = 6,
                            clip = 10.0, missed = 10.0, direct = 10.0,
                            direct_names = nothing, drive = nothing, drive_weight = 1.0,
                            center = nothing, keepout = nothing, uncertainty = nothing)
    errs = _uncertainty_spec(uncertainty)
    errs === nothing || (center === nothing && lmax <= _CACHED_LMAX) ||
        throw(ArgumentError("uncertainty needs the default center and lmax ≤ $_CACHED_LMAX"))
    return function (ev::Evaluation)
        a = flux_moments(ev; lmax, center)
        J = sum(w * (a[l+1] / a[1])^2 for (l, w) in modes)
        odd > 0 && (J += odd * sum((a[l+1] / a[1])^2 for l in 1:2:lmax))
        if errs !== nothing
            v = _error_variance(ev.scenario, ev.effective, ev.scenario.moment_response, lmax + 1, errs)
            J += sum(w * v[l+1] for (l, w) in modes)
            odd > 0 && (J += odd * sum(v[l+1] for l in 1:2:lmax))
        end
        J += _penalties(ev; clip, missed, direct, direct_names)
        J += keepout_penalty(keepout, ev.beams)
        drive === nothing || (J += drive_weight * (a[1] / drive - 1)^2)
        return J
    end
end

# The errors of an objective's `uncertainty` keyword (`mode_uncertainty`'s keywords), or
# nothing.
_uncertainty_spec(::Nothing) = nothing
function _uncertainty_spec(u)
    e = (; pointing = Float64(get(u, :pointing, 0.0)), power = Float64(get(u, :power, 0.0)),
         group_pointing = Float64(get(u, :group_pointing, 0.0)),
         group_power = Float64(get(u, :group_power, 0.0)), groups = get(u, :groups, :quad),
         positioning = _positioning_rms(get(u, :positioning, 0.0)))
    min(e.pointing, e.power, e.group_pointing, e.group_power) >= 0 ||
        throw(ArgumentError("uncertainty must be nonnegative rms errors"))
    return e
end

function _penalties(ev::Evaluation; clip, missed, direct, direct_names)
    J = zero(ev.deposition.delivered)
    clip > 0 && (J += clip * clipped_fraction(ev)^2)
    missed > 0 && (J += missed * missed_fraction(ev)^2)
    direct > 0 && (J += direct * _direct_fraction(ev, direct_names)^2)
    return J
end

# The fraction of the delivered laser power absorbed by the observed elements, or by the
# named parts or components (a fitted CAD capsule is the part of its pieces).
_direct_fraction(ev::Evaluation, ::Nothing) =
    sum(ev.deposition.element[e] for e in ev.scenario.observed; init = zero(ev.deposition.delivered)) /
    ev.deposition.delivered
function _direct_fraction(ev::Evaluation, names)
    cs = ev.scenario.target.components
    for n in names
        any(c -> c.name === n || c.part === n, cs) ||
            throw(ArgumentError("direct_names: the target has no part or component $(repr(n)) " *
                                "(parts: $(join(repr.(parts(ev.scenario.target)), ", ")))"))
    end
    return power_fraction(ev, names...)
end

"""
    flux_ylm(ev; lmax = 6) -> Vector

Real spherical-harmonic coefficients of the observed incident flux in the flat
`mode_index` order (polar axis ẑ, azimuth from +x, about the observed capsule's center:
`mode_center`, TCC unless the capsule is offset), from the scenario's cached
`ylm_weights`. Use `ylm_matrix(c, lmax)` for the `C[l+1, lmax+1+m]` layout.
"""
function flux_ylm(ev::Evaluation; lmax::Integer = _CACHED_YLM_LMAX)
    _radiative(ev.scenario) || throw(ArgumentError("scenario has no radiation response"))
    lmax <= _CACHED_YLM_LMAX || throw(ArgumentError("cached Y_lm weights reach lmax = $_CACHED_YLM_LMAX"))
    return _cached_product(getfield(ev, :_ylm), ev.scenario.ylm_response, ev.deposition.element)[1:(lmax+1)^2]
end

"""
    asymmetry_objective(; modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0), components = (),
                        clip = 10.0, missed = 10.0, direct = 10.0, direct_names = nothing,
                        drive = nothing, drive_weight = 1.0, keepout = nothing,
                        reference = nothing, uncertainty = nothing)

Objective for three-dimensional drive symmetry. For each `l => w` in `modes` it adds
`w × relative_mode_power(c, l)` — the rotation-invariant power of all `2l+1` real
harmonics of degree `l` over the mean (`(aₗ/a₀)²/(2l+1)` for an axisymmetric field) — and
for each `(l, m) => w` in `components` adds `w × (c_lm/c_00)²`, so a particular
azimuthal mode (e.g. the `l = 1, m = ±1` drift of a dropped bundle) can be weighted
separately. Penalties for clipping, missing and direct capsule illumination and the
`drive` target and `keepout` are as in `symmetry_objective`, with `drive` compared to
the mean flux `c_00/√(4π)`.

With `reference` (a coefficient vector from `flux_ylm`, e.g. of a tuned design before a
change to the target or the beams), the modes and components measure the departure of
`c_lm/c_00` from the reference's `c_lm/c_00` instead of from zero, so the objective
restores a known-good drive shape.

With `uncertainty = (; pointing, power)` the mode and component terms are their
expected values over independent beam errors of those rms sizes (linearized; see
`mode_uncertainty` and `symmetry_objective`): a robust design. Shared errors
(`group_pointing`, `group_power`, `groups`, `positioning`) are accepted as in
`mode_uncertainty`.
"""
function asymmetry_objective(; modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0),
                             components = (), clip = 10.0, missed = 10.0, direct = 10.0,
                             direct_names = nothing, drive = nothing, drive_weight = 1.0,
                             keepout = nothing, reference = nothing, uncertainty = nothing)
    errs = _uncertainty_spec(uncertainty)
    lmax = maximum(vcat([l for (l, _) in modes], [lm[1] for (lm, _) in components], 0))
    nref = (max(lmax, 1) + 1)^2
    reference === nothing || length(reference) >= nref ||
        throw(ArgumentError("reference needs coefficients through l = $(max(lmax, 1))"))
    shape = reference === nothing ? zeros(nref) : reference[1:nref] ./ reference[1]
    return function (ev::Evaluation)
        c = flux_ylm(ev; lmax = max(lmax, 1))
        d = c .- c[1] .* shape                 # departure from the reference shape
        J = zero(c[1])
        for (l, w) in modes
            J += w * sum(abs2, view(d, l^2+1:(l+1)^2)) / c[1]^2   # relative_mode_power
        end
        for ((l, m), w) in components
            J += w * (d[mode_index(l, m)] / c[1])^2
        end
        if errs !== nothing
            v = _error_variance(ev.scenario, ev.effective, ev.scenario.ylm_response, nref, errs)
            for (l, w) in modes
                J += w * sum(@view v[l^2+1:(l+1)^2])
            end
            for ((l, m), w) in components
                J += w * v[mode_index(l, m)]
            end
        end
        J += _penalties(ev; clip, missed, direct, direct_names)
        J += keepout_penalty(keepout, ev.beams)
        drive === nothing || (J += drive_weight * (c[1] / sqrt(4π) / drive - 1)^2)
        return J
    end
end

"""
    influence_matrix(scenario, beams; lmax = 6) -> Matrix

Mode response of the observed flux to each beam: column `b` is the real-harmonic
coefficient vector (`mode_index` order) produced by beam `b` alone at 1 TW, with the
beam's pointing, spot and focus. Without a laser–plasma model the flux is linear in
beam powers, so the coefficients of any power assignment `P` are `M * P` — the basis of
power-rebalancing solves and of tables of what losing each beam costs. `beams` are the
effective beams (after `apply_lpi`); a laser–plasma model that deposits light itself
(the ridge stop of `BubbleAbsorption`, linear in each beam's power) is included through
`lpi_deposit` at the scenario's time.
"""
function influence_matrix(s::Scenario, beams::AbstractVector{<:Beam};
                          lmax::Integer = _CACHED_YLM_LMAX)
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    nb = length(beams)
    M = zeros((lmax + 1)^2, nb)
    lmax <= _CACHED_YLM_LMAX || throw(ArgumentError("cached Y_lm weights reach lmax = $_CACHED_YLM_LMAX"))
    Wy = view(s.ylm_response, 1:(lmax+1)^2, :)
    for i in 1:nb
        unit = [remake(b; power = j == i ? 1.0 : 0.0) for (j, b) in enumerate(beams)]
        dep = lpi_deposit(s.lpi, s.deposition, unit, s.time)
        M[:, i] = Wy * dep.element
    end
    return M
end

"""
    irradiance_objective(name; center, radius, intensity = nothing, uniformity = 1.0,
                         spill = 1.0, intensity_weight = 1.0, side = 1)

Objective for a laser-irradiated foil (e.g. a backlighter): RMS nonuniformity of the
absorbed laser intensity over the elements of component `name` within `radius` of
`center` (on side `side`: 1 front, -1 back), plus `spill ×` the squared fraction of
delivered power landing outside that region, plus — if `intensity` (TW/mm²) is given —
`intensity_weight × (mean/intensity - 1)²`.
"""
function irradiance_objective(name::Symbol; center, radius, intensity = nothing,
                              uniformity = 1.0, spill = 1.0, intensity_weight = 1.0, side = 1)
    c0 = Vec3(center...)
    return function (ev::Evaluation)
        t = ev.scenario.target
        es = [e for e in elements(t, name) if t.elem_side[e] == side &&
              norm(t.elem_centroid[e] - c0) <= radius]
        isempty(es) && throw(ArgumentError("no elements of $name within radius"))
        P = ev.deposition.element[es]
        A = t.elem_area[es]
        J = uniformity * rms_nonuniformity(P ./ A, A)^2
        spill > 0 && (J += spill * (1 - sum(P) / ev.deposition.delivered)^2)
        if intensity !== nothing
            J += intensity_weight * (sum(P) / sum(A) / intensity - 1)^2
        end
        return J
    end
end

"""
    face_flux(ev) -> Vector

Incident x-ray flux `H_e` (TW/mm²) on each of the scenario's observed elements: the
full linear flux response (`scenario.response`, computed when first read) times the
absorbed laser power, generic in the number type (the optimizer's duals). On a dense
model, the same numbers as the evaluation's `flux`; on a large (iterative) model it
builds the response (a solve per observed element) instead of one iterative solve, so
meant for an observed face of moderate size.
"""
function face_flux(ev::Evaluation)
    s = ev.scenario
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    return _mul(_response(s)::Matrix{Float64}, ev.deposition.element)
end

"""
    face_mean_flux(ev) -> Real

Area mean `H̄ = Σ A_e H_e / Σ A_e` (TW/mm²) of the incident x-ray flux on the scenario's
observed elements ([`face_flux`](@ref)): the level of [`FaceFlux`](@ref) and the level
that the `flux` keyword of [`face_flux_objective`](@ref) holds, so
`flux = face_mean_flux(evaluate(scenario, x0))` keeps the drive of the design `x0`.
"""
face_mean_flux(ev::Evaluation) = _face_terms(ev)[1]

"""
    face_nonuniformity(ev) -> Real

Area-weighted RMS nonuniformity `σ = rms_nonuniformity(H, A)` (relative to the mean) of
the incident x-ray flux on the scenario's observed elements ([`face_flux`](@ref)): the
shape term of [`face_flux_objective`](@ref), and the scale of its `rms_scale` keyword,
so `rms_scale = face_nonuniformity(evaluate(scenario, x0))` weighs the start's
nonuniformity as 1.
"""
face_nonuniformity(ev::Evaluation) = _face_terms(ev)[2]

# The area mean and the RMS nonuniformity of the observed face's flux.
function _face_terms(ev::Evaluation)
    H = face_flux(ev)
    A = ev.scenario.target.elem_area[ev.scenario.observed]
    return sum(A .* H) / sum(A), rms_nonuniformity(H, A)
end

"""
    face_flux_objective(; uniformity = 1.0, rms_scale = nothing, flux = nothing,
                        flux_weight = 1.0, clip = 10.0, missed = 10.0, direct = 10.0,
                        direct_names = nothing, keepout = nothing)

Objective for the x-ray drive of a flat observed face (a package's face in a halfraum,
an opacity sample), which has no capsule modes: `uniformity × (σ/σ₀)²` with `σ =
rms_nonuniformity(H, A)` the area-weighted RMS nonuniformity of the incident flux on
the scenario's observed elements ([`face_nonuniformity`](@ref), through the scenario's
linear flux response, so differentiable), the shape term of [`FaceFlux`](@ref) for a
uniform target, and `σ₀ = rms_scale` (`nothing`: 1, the absolute `σ²`). When `flux`
(TW/mm²) is given, adds `flux_weight × (H̄/flux - 1)²` with `H̄` the
[`face_mean_flux`](@ref). The `clip`, `missed` and `direct` penalties (`direct`: laser
light on the face itself, or with `direct_names` on the named parts or components) and
the `keepout` penalty are as in [`symmetry_objective`](@ref).

`rms_scale = face_nonuniformity(evaluate(scenario, x0))` makes the term 1 at the design
`x0` (times `uniformity`), so `uniformity = 1` balances the shape against the penalties
on any face. Without it the term is the absolute `σ²`, which is as small as the face is
uniform: on a uniform 1 mm sample (σ ~ 1e-4, σ² ~ 1e-8) the penalties (10 × a squared
light fraction, ~1e-2) set the optimum, while a capsule's modes (`(Pₗ/P₀)²` ~ 1e-2
before tuning) are of their size. The scale is a number fixed at construction, so the
objective stays a function of the beams alone (thread-safe, differentiable).
"""
function face_flux_objective(; uniformity = 1.0, rms_scale = nothing, flux = nothing,
                             flux_weight = 1.0, clip = 10.0, missed = 10.0, direct = 10.0,
                             direct_names = nothing, keepout = nothing)
    uniformity >= 0 || throw(ArgumentError("uniformity must not be negative"))
    rms_scale === nothing || (isfinite(rms_scale) && rms_scale > 0) ||
        throw(ArgumentError("rms_scale must be a positive RMS nonuniformity"))
    flux === nothing || flux > 0 || throw(ArgumentError("the flux level must be positive"))
    w = rms_scale === nothing ? Float64(uniformity) : uniformity / Float64(rms_scale)^2
    return function (ev::Evaluation)
        H̄, σ = _face_terms(ev)
        J = w * σ^2
        J += _penalties(ev; clip, missed, direct, direct_names)
        J += keepout_penalty(keepout, ev.beams)
        flux === nothing || (J += flux_weight * (H̄ / flux - 1)^2)
        return J
    end
end

"""
    optimize_pointing(scenario, objective; from = nothing,
                      x0 = initial_parameters(scenario.param),
                      lower = bounds(scenario.param, x0)[1], upper = bounds(scenario.param, x0)[2],
                      iterations = 200, evaluations = 20 * iterations, show_trace = false,
                      g_tol = 1e-8, callback = nothing, constraints = (),
                      feasibility = 1e-3, outer = 30, time_limit = Inf, phase1 = true,
                      relaxation = false, inner_reltol = 1e-10,
                      inner_evaluations = nothing, smooth_cut = true) -> PointingResult

Minimize `objective(evaluate(scenario, x))` within the box `lower ≤ x ≤ upper` by
L-BFGS-B with forward-mode automatic-differentiation gradients, from `x0` (clamped to
the box). Returns a [`PointingResult`](@ref), whose display lists the parameters that
end at a bound.

The default box is `bounds(scenario.param, x0)` with its default widths: aims within
±1 mm of `x0` (`z`, `r`, `phi`), defocus ±2 mm, power and spot scale ±50 %, Δλ ±2 Å and
the spherical angles ±5–10° (engineering defaults, not facility limits). A parameter at
its bound usually means the box, not the objective, set it: widen it with keywords of
[`bounds`](@ref), e.g. `lower, upper = bounds(pointing, x0; z = 1.5, power = 0.4,
power_max = 12.0)`, or pass any `lower` and `upper` vectors. A parameterization other
than `GroupPointing` and [`FramePointing`](@ref) has no `bounds` method and needs
`lower` and `upper`.

`scenario` may be a [`PhasedScenario`](@ref): one pointing over several phases of a
pulse, with power scales per (group, phase). The objective and every constraint then
receive the [`PerPhase`](@ref) of the phases' evaluations at each parameter vector,
evaluated once, and the result's `beams` are a `PerPhase` of beam vectors; everything
else (constraints, `from`, `phase1`, `relaxation`, the time limit) is the same.

The run ends after `iterations` L-BFGS-B iterations or `evaluations` objective
evaluations (each gradient counts one), whichever comes first: a line search may take
many evaluations in one iteration, on a kinked objective (keep-out hinges, a coarse
mesh) dozens, so the iteration limit alone does not bound the run time. `time_limit`
(seconds of wall time, default none) ends the run at the first evaluation after it has
passed. `reason` is `:converged` (projected gradient below `g_tol`, or no further
progress), `:iterations`, `:evaluations`, `:time_limit` or `:callback`. On
`:evaluations`, `:time_limit` and `:callback` the result is the best design evaluated
so far, with `stopped = true` and `result = nothing`; `iterations` is the number of
iterations completed in every case.

`callback`, if given, is called after every (non-derivative) objective evaluation with
`(; evaluations, iteration, value, best_value, best_x)`; returning `true` stops the
optimization. Use it for progress reports and cancellation.

**Constraints.** `constraints`, a tuple or vector of [`Constraint`](@ref)s (ceilings
`f(ev) ≤ limit` and floors `f(ev) ≥ limit`, such as [`irradiance_ceiling`](@ref) and
[`mode_gate`](@ref)), are handled as constraints rather than as penalty terms of the
objective: by an augmented Lagrangian, whose inner problems are box-constrained
L-BFGS-B runs on `objective + Σ μ/2 [max(0, gᵢ + λᵢ/μ)² - (λᵢ/μ)²]` (`gᵢ` the
normalized violations, `±(value - limit)/scale`), each started from the last, with
the multipliers `λ` and the penalty `μ` updated between them (ALGENCAN's rules,
Andreani et al., SIAM J. Optim. 18, 1286, 2008). The run converges when every
violation is at most `feasibility` (in the constraints' scaled units: 1e-3 is 0.1 % of
a limit by default), the multipliers are complementary to the same tolerance, and the
last inner run has converged at the final tolerance: the inner runs' tolerance on the
projected gradient tightens tenfold per outer iteration from 1e-3 to `g_tol`, and the
run ends only once an inner run at `g_tol` has converged (its projected gradient below
it, or its value changing by less than `inner_reltol` relative per iteration: a drive
objective is flat along the limits), so a start that already meets the limits, as on a
restart or a hand-off to a finer parameterization, is still optimized. A constraint that cannot
be met ends the run at phase 1's verdict (below), or with `phase1 = false` once the
penalty has grown by 10⁸ without progress. `outer` caps the outer
iterations, `iterations` caps each inner run and `evaluations` is the budget of the
whole run (the result counts both over all inner runs). The multipliers change only
between inner runs, so the cap sets how long the run works with each set of them. An
inner run ended by its cap has not solved its subproblem, which the penalty's test
(has the infeasibility–complementarity measure halved?) presumes: after one that ends
above the limits the multipliers are updated, while the penalty, the inner tolerance
and the test's reference are held, so that a cap never grows the penalty toward the
infeasible verdict; one that ends within the limits takes the usual test (the result's
`capped` counts both; research notes, "Least violation, capped inner runs and phase
1's budget"). A run that ends above its limits (infeasible, or stopped
by a budget, the time limit or the callback) returns its iterate of least violation,
not its last ([`PointingResult`](@ref): `least_violation`, `last_feasible`). The
public cases of up to about a hundred parameters converge with a large cap (the
dropped-bundle tutorial passes 1,000) in about a thousand evaluations; small early
caps, measured as an alternative, saved evaluations on some of them and cost drive on
others. On a problem with hundreds of parameters, whose evaluations take seconds and
whose start is far from the limits, a cap of about 100 updates the multipliers sooner
(a 276-parameter study with `iterations = 1000` completed one outer iteration in an
hour, with 100 six in 40 minutes; research notes, "The inner budget"). The cap on
iterations does not bound an inner run's evaluations (line searches); an inner run can
spend hundreds of them on 1e-4 of the objective while limits stay violated.
`inner_evaluations` (default none) caps the evaluations of each inner run: one that
reaches it ends as one ended by `iterations` does, at its best design, and the outer
loop updates the multipliers and the penalty from there (phase 1 is not capped; the
relaxations' runs are). A capped inner run never counts as converged, so the run
still ends only after an inner run converges within the cap at the final tolerance,
and a small cap needs a larger `outer`: on a coarse version of the dropped-bundle
tutorial's run within limits (92 parameters), whose default run spent 795 of its 1,044
evaluations in one inner run, a cap of 100 converged in 10 outer iterations and 766
evaluations at 0.06 % less drive, and a cap of 30
ended at `outer = 30` unconverged (research notes, "Inner-run options, a wall-loss
scenario warning and laser delivery figures").
`inner_reltol` (default 1e-10; 0 turns it off) is the inner runs' relative value
tolerance above, also phase 1's. The default ends a drive objective's inner runs where
the objective is flat along the limits, short of stationarity (the result's
`stationarity` shows how far); a tighter one costs evaluations: on the dropped-bundle
tutorial's run within limits, 1e-12 reaches a stationarity of 1.9e-5 against 1.2e-3 for
2.1 times the evaluations, the drive the same within the flat family's 1e-4 (research
notes, "Multipliers from the converged active set"). A result whose multipliers or
bound sensitivities matter continues tighter through `from` with a smaller
`inner_reltol`. The result's
`constraints` and `outcome` report each constraint against its limit, with the
multipliers that best satisfy stationarity at the final design, and whether the box,
the constraints or the budget ended the run ([`PointingResult`](@ref)). The
callback's named tuple then also carries `outer`, and its `value` is the inner
problem's (the objective with the constraint terms). Without constraints the run is
the plain L-BFGS-B above.

**Infeasible constraints.** `phase1 = true` (the default) first minimizes, from a start that violates
a constraint, the infeasibility `½ Σ max(0, gᵢ)²` over every component's normalized
violation (the augmented Lagrangian's limit as its penalty grows), ignoring the
objective, by the same L-BFGS-B. It stops at the first design where every violation is
within `feasibility`, and the main run starts there. If it converges, or stalls (its
least infeasibility falls by less than 1 % over 100 evaluations, an engineering
default), with a violation above `feasibility` instead, the run ends there with
`outcome` `:infeasible` (or
`:infeasible_box`, parameters at a bound) and no outer iteration: that design is the
least infeasible one the search reached from this start (a local verdict), and its
`constraints` rows show each constraint's violation there. `relaxation = true` (which
implies `phase1`) adds, after such a verdict, one constrained run per constraint
violated or active there: the smallest largest violation of its components with every
other constraint met, which is how far that constraint alone must be relaxed
([`Relaxation`](@ref), in `PointingResult.relaxations`); these runs share the budget
and the time limit. `phase1` is on by default and `relaxation` off: a run whose start
violates a constraint starts its main run from phase 1's design and can end at another
design of a flat family than one started at `x0`; `phase1 = false` starts the main run
at `x0`. A start that already meets every constraint, and a restart with `from`, skip
phase 1. The
callback's `outer` is 0 during phase 1. On the dropped-bundle test case made
infeasible, phase 1 reached the verdict in 577 evaluations where the augmented
Lagrangian had none after 2,195, and on feasible starts it cost 8 % fewer to 29 % more
evaluations (research notes, "A verdict on conflicting constraints"); on the test
suite's constrained cases that start outside their limits, −23 % to +54 %, with drive
designs 3e-4 apart on a flat family (research notes, "Defaults changed (2026-11-19)").

**VisRad's hard cut.** On a scenario whose deposition uses [`VisRadProfile`](@ref)
with its hard cut (`taper = 0`, its default), every deposition-derived quantity jumps
where an operator point crosses the cut, and line searches park on the jumps. With
`smooth_cut = true` (the default) the run optimizes the same scenario with the profile
tapered over one e-fold below the cut (`taper = 1`: the operator's points and screens
and the scenario's responses are shared, nothing is rebuilt) and reports the design
under the scenario's own cut: the result's `value`, `initial_value`, `constraints` and
`outcome` are the hard cut's at the run's design and start, and `tapered_value` holds
the objective the tapered run reached. The taper moves wall elements by at most 7e-5 of
the wall's peak on NIF hohlraums, well inside `feasibility`'s default 1e-3, so a
design within its limits under the taper is normally within them under the cut; where
it is not, its row reads `:violated` and `outcome` `:infeasible`. The multipliers
(fitted under the taper), stationarity, bound sensitivities, phase 1, relaxations and
`state` are the tapered run's; a restart with `from` continues it. `smooth_cut = false`
optimizes the scenario's own profile; other scenarios (the f-cone model, a profile
already tapered) are optimized as they are. A [`PhasedScenario`](@ref) is tapered
phase by phase.

**Restarts.** `from`, a previous run's [`PointingResult`](@ref), continues that run:
`x0`, `lower` and `upper` default to its last iterate (`state.x` of a constrained run,
its `x` otherwise and whenever the run returned its last iterate) and box, and a constrained run starts from
its multipliers, penalty parameter, inner tolerance and working limits (`result.state`)
instead of learning them again from zero. A run cut short by `evaluations` (or by the
callback) and continued with `from` resumes the interrupted inner problem; the counts of
the new result (`iterations`, `evaluations`, `outer`) are its own, and its
`initial_value` is the objective at its start. The constraints must have the same
components as the previous run's (a `DimensionMismatch` otherwise). An unconstrained
`from` result only sets the start and the box.

Compared with penalty terms at fixed weights, the multipliers carry the constraints'
weights from one inner run to the next: the limits are met without stiff weights
(which stall the line search), without stages of loosened limits set by hand, and the
design ends at its limits rather than where the penalty terms balance the objective
(on the dropped-bundle tutorial, 0.35 % more drive in a comparable number of
evaluations, 997 against the penalty form's 1,072: research notes, "Constrained
re-pointing").
"""
Base.@nospecializeinfer function optimize_pointing(
                           s::_AbstractScenario, @nospecialize(objective);
                           @nospecialize(from = nothing),
                           x0 = from === nothing ? initial_parameters(s.param) : _restart_x(from),
                           lower = from === nothing ? bounds(s.param, x0)[1] : from.lower,
                           upper = from === nothing ? bounds(s.param, x0)[2] : from.upper,
                           iterations::Integer = 200, evaluations::Integer = 20 * iterations,
                           show_trace::Bool = false, g_tol = 1e-8, @nospecialize(callback = nothing),
                           @nospecialize(constraints = ()), feasibility::Real = 1e-3,
                           outer::Integer = 30, time_limit::Real = Inf, phase1::Bool = true,
                           relaxation::Bool = false, inner_reltol::Real = 1e-10,
                           @nospecialize(inner_evaluations = nothing), smooth_cut::Bool = true)
    # Inferred and compiled once for every objective and callback (each objective
    # closure type used to recompile this body and Optim's L-BFGS-B, ~1 s); the objective
    # is reached by dynamic dispatch once per evaluation.
    deadline = time() + Float64(time_limit)
    time_limit > 0 || throw(ArgumentError("time_limit must be positive (seconds)"))
    # VisRad's hard cut: the run on the profile tapered below the cut, reported under the
    # scenario's own cut
    st = smooth_cut ? _tapered(s) : nothing
    if st !== nothing
        r = optimize_pointing(st, objective; from, x0, lower, upper, iterations, evaluations,
                              show_trace, g_tol, callback, constraints, feasibility, outer,
                              time_limit, phase1, relaxation, inner_reltol, inner_evaluations,
                              smooth_cut = false)
        return _hard_cut_result(s, r, objective, constraints, clamp.(x0, lower, upper),
                                Float64(feasibility))
    end
    from === nothing || from isa PointingResult ||
        throw(ArgumentError("`from` takes a PointingResult, got $(typeof(from))"))
    xs = clamp.(x0, lower, upper)
    inner_evaluations === nothing || inner_evaluations isa Integer ||
        throw(ArgumentError("inner_evaluations takes an integer or `nothing`"))
    inner_evaluations === nothing || !isempty(constraints) ||
        throw(ArgumentError("inner_evaluations caps the inner runs of a constrained run; without " *
                            "constraints `evaluations` is the run's budget"))
    isempty(constraints) ||
        return _constrained_pointing(s, objective, Tuple(constraints), x0, xs, lower, upper;
                                     iterations, evaluations, show_trace, g_tol, callback,
                                     feasibility = Float64(feasibility), outer = Int(outer),
                                     state = from === nothing ? nothing : from.state, deadline,
                                     phase1, relaxation, inner_reltol = Float64(inner_reltol),
                                     inner_evaluations = inner_evaluations === nothing ? typemax(Int) :
                                                         Int(inner_evaluations))
    r = _lbfgsb(s, objective, x0, xs, lower, upper; iterations, evaluations, show_trace, g_tol,
                callback, deadline)
    return _pointing_result(s, r.x, r.value, r.initial_value, r.iterations, r.evaluations,
                            r.stopped, r.reason, r.result, lower, upper; objective)
end

# A restart's start: the iterate whose multipliers the result's `state` holds (`state.x`;
# the result's `x` unless the run ended above its limits and returned its iterate of
# least violation), else the result's `x`.
# Declared `Vector{Float64}`: `from` is not specialized on, and an inferred `Any` here would
# compile `optimize_pointing`'s body a second time.
function _restart_x(r)::Vector{Float64}
    st = r.state
    return st isa NamedTuple && haskey(st, :x) ? st.x : r.x
end

# The gradient of `F(evaluate(s, x))` at `x`, through the optimizer's dual tag and chunks
# (bitwise ForwardDiff's chunk-mode gradient; see `_pointing_gradient!`).
Base.@nospecializeinfer function _scenario_gradient(s::_AbstractScenario, @nospecialize(F), x::Vector{Float64})
    J(x) = F(evaluate(s, x))
    cfg = ForwardDiff.GradientConfig(J, x, ForwardDiff.Chunk(min(length(x), 12)),
                                     ForwardDiff.Tag(_PointingTag(), Float64))
    out = ForwardDiff.DiffResults.GradientResult(similar(x))
    return copy(ForwardDiff.DiffResults.gradient(_pointing_gradient!(out, J, x, cfg)))
end

# The Jacobian of the vector `F(evaluate(s, x))` at `x` (one row per entry), through the
# optimizer's dual tag and chunks: each chunk's dual pass carries every entry, so the
# Jacobian costs one gradient whatever its number of rows. The passes run concurrently
# as in `_pointing_gradient!`; each row is bitwise ForwardDiff's chunk-mode gradient of
# that entry.
Base.@nospecializeinfer function _scenario_jacobian(s::_AbstractScenario, @nospecialize(F), x::Vector{Float64})
    J(x) = F(evaluate(s, x))
    cfg = ForwardDiff.GradientConfig(J, x, ForwardDiff.Chunk(min(length(x), 12)),
                                     ForwardDiff.Tag(_PointingTag(), Float64))
    return _pointing_jacobian(J, x, cfg)
end

function _pointing_jacobian(f, x::Vector{Float64}, cfg::ForwardDiff.GradientConfig{T,Float64,N}) where {T,N}
    nc = cld(length(x), N)
    blocks = if nc <= 1 || Threads.nthreads() == 1 || !_PARALLEL_CHUNKS[]
        [_jacobian_chunk(f, x, cfg.seeds, T, c) for c in 1:nc]
    else
        tasks = [Threads.@spawn _jacobian_chunk(f, x, cfg.seeds, T, c) for c in 1:nc]
        map(tasks) do t
            try
                fetch(t)
            catch err
                err isa TaskFailedException ? throw(err.task.exception) : rethrow()
            end
        end
    end
    return reduce(hcat, blocks)
end

# Dual pass `c` of `_pointing_jacobian`: the columns of its own parameters.
function _jacobian_chunk(f, x, seeds::NTuple{N,ForwardDiff.Partials{N,Float64}}, ::Type{T}, c) where {N,T}
    lo = (c - 1) * N
    zp = zero(ForwardDiff.Partials{N,Float64})
    xd = [ForwardDiff.Dual{T}(x[i], lo < i <= lo + N ? seeds[i-lo] : zp) for i in eachindex(x)]
    y = f(xd)
    k = min(N, length(x) - lo)
    return Float64[y[r] isa ForwardDiff.Dual ? ForwardDiff.partials(y[r], j) : 0.0
                   for r in eachindex(y), j in 1:k]
end

# One L-BFGS-B run of `objective(evaluate(s, x))` in the box from `xs` (`x0` sizes the
# gradient's configuration): the parameters, objective value, iterations, evaluations
# and stop of `optimize_pointing`, and Optim's result (`nothing` when stopped).
# `f_reltol` (Optim's; 0, its default, off) ends the run on a small relative change of
# the value, for the inner runs of the constrained path.
Base.@nospecializeinfer function _lbfgsb(s::_AbstractScenario, @nospecialize(objective), x0, xs, lower,
                                         upper; iterations::Integer, evaluations::Integer,
                                         show_trace::Bool, g_tol, @nospecialize(callback),
                                         f_reltol::Float64 = 0.0, deadline::Float64 = Inf,
                                         cap::Int = typemax(Int))
    best_x, best_value, count, iter = copy(xs), Inf, Ref(0), Ref(0)
    J(x) = objective(evaluate(s, x))
    # every objective value the optimizer sees (bookkeeping, callback, budget)
    function seen(x, v)
        count[] += 1
        v < best_value && ((best_value, best_x) = (v, copy(x)))
        callback === nothing ||
            callback((; evaluations = count[], iteration = iter[], value = v, best_value,
                      best_x)) === true && throw(_OptimizationStop(:callback))
        count[] >= evaluations && throw(_OptimizationStop(:evaluations))
        time() > deadline && throw(_OptimizationStop(:time_limit))
        count[] >= cap && throw(_OptimizationStop(:cap))
        return v
    end
    f(x) = seen(x, J(x))
    # One dual tag for every objective, so the compiled deposition and objective code
    # (and the package's precompiled copy) serve any objective, not only the first.
    # ForwardDiff's tag check (`Val(false)` in `_pointing_gradient!`) compares the tag's
    # function with `J` and would refuse this shared tag; it guards against a config
    # reused for another function, not against nesting. An objective that differentiates
    # internally makes its own tag, which `_PointingTag`'s ordering (`≺`) nests inside
    # this one whatever the tags' counts (test/nested_gradient.jl).
    cfg = ForwardDiff.GradientConfig(J, x0, ForwardDiff.Chunk(min(length(x0), 12)),
                                     ForwardDiff.Tag(_PointingTag(), Float64))
    # The chunks' dual passes run concurrently (`_pointing_gradient!`); the result is
    # bitwise ForwardDiff's serial chunk-mode gradient.
    out = ForwardDiff.DiffResults.GradientResult(similar(xs))
    g!(G, x) = copyto!(G, ForwardDiff.DiffResults.gradient(_pointing_gradient!(out, J, x, cfg)))
    # L-BFGS-B asks for value and gradient together: the dual pass carries the value
    function fg!(G, x)
        _pointing_gradient!(out, J, x, cfg)
        copyto!(G, ForwardDiff.DiffResults.gradient(out))
        return seen(x, ForwardDiff.DiffResults.value(out))
    end
    initial_value = J(xs)
    track(st) = (iter[] = st.iteration; false)
    res = try
        Optim.optimize(Optim.OnceDifferentiable(f, g!, fg!, xs, 0.0), lower, upper, xs,
                       Optim.LBFGSB(),
                       Optim.Options(; iterations, show_trace, g_abstol = g_tol, f_reltol,
                                     callback = track))
    catch err
        err isa _OptimizationStop || rethrow()
        err
    end
    if res isa _OptimizationStop
        # the cap on one run's evaluations ends it as the iteration limit does: not a stop
        # of the whole optimization
        capped = res.reason === :cap
        return (; x = best_x, value = best_value, initial_value, iterations = iter[],
                evaluations = count[], stopped = !capped, reason = capped ? :iterations : res.reason,
                result = nothing)
    end
    n = Optim.iterations(res)
    return (; x = Optim.minimizer(res), value = Optim.minimum(res), initial_value,
            iterations = n, evaluations = count[], stopped = false,
            reason = Optim.converged(res) || n < iterations ? :converged : :iterations,
            result = res)
end

"""
    ConstraintRow

One component of a [`Constraint`](@ref) at a design, from [`check_constraints`](@ref) or
in `PointingResult.constraints`. Fields: `name` (the constraint's name, with the
component's label for a vector constraint, e.g. `"mode_gate[l = 2]"`), `value` and
`limit` (in the constraint's units), `upper` (`true` for a ceiling `value ≤ limit`,
`false` for a floor), `violation` (`±(value - limit)/scale`, the normalized excess:
positive when violated, negative inside), `multiplier` (the component's Lagrange
multiplier at the design, in the normalized units: in a constrained run's result, the
nonnegative least-squares solution of the stationarity condition on the parameters
inside the box, over the components within `feasibility` of their limits or beyond
them, zero for the other components ([`PointingResult`](@ref)); 0 when not computed,
as in [`check_constraints`](@ref)) and `status`:
`:violated` (`violation > feasibility`), `:active` (within `feasibility` of the limit)
or `:satisfied`.
"""
struct ConstraintRow
    name::String
    value::Float64
    limit::Float64
    upper::Bool
    violation::Float64
    multiplier::Float64
    status::Symbol
end

"""
    Relaxation

How far one constraint must be relaxed, alone, for the others to hold: an entry of
`PointingResult.relaxations` after a run that [`optimize_pointing`](@ref) found
infeasible with `relaxation = true`. Fields: `name` (the constraint's), `violation`
(the largest normalized violation of its components at a design that meets every
other constraint: the relaxation in units of its `scale`, so a multiple `1 +
violation` of the limit where the scale is the limit, the default), `limit` and
`relaxed_limit` (one per component, in the constraint's units: a ceiling raised and a
floor lowered by `violation × scale`), `upper`, `status`, `x` (that design),
`constraints` (every constraint's [`ConstraintRow`](@ref) at `x`) and `evaluations`.
`status` is `:relaxed` (the others met at `x`), `:others_conflict` (the other
constraints cannot all be met either, so relaxing this one alone does not suffice;
`violation` is `NaN`) or `:stopped` (the budget, the time limit or the callback ended
the search first; `violation` is the best design's, the others possibly unmet).

The relaxation is achievable (the design `x` meets the others with this one at
`relaxed_limit`) and the least one only locally: the run that finds it minimizes a
smooth stand-in for the largest violation (a 64-norm of the positive violations) from
the infeasible design, a local search.
"""
struct Relaxation
    name::String
    violation::Float64
    limit::Vector{Float64}
    relaxed_limit::Vector{Float64}
    upper::Bool
    status::Symbol
    x::Vector{Float64}
    constraints::Vector{ConstraintRow}
    evaluations::Int
end

"""
    PointingResult

Result of [`optimize_pointing`](@ref). Fields: `x` (the best parameters), `beams`
(`apply(scenario.param, x)`; a [`PerPhase`](@ref) of beam vectors for a phased
scenario), `value` and `initial_value` (the objective at `x` and at
the clamped `x0`), `iterations` (L-BFGS-B iterations completed), `evaluations`
(objective evaluations, each gradient counting one), `reason` (`:converged`,
`:iterations`, `:evaluations`, `:time_limit` or `:callback`), `stopped` (`true` for the
last three, which
return the best design evaluated so far), `result` (Optim's result, `nothing` when
`stopped`), `names` (`parameter_names` of the parameterization, `x[i]` for one without
a method), the box `lower`, `upper`, and `at_bounds`: the names of the parameters at a
bound (within 1e-8 of the box width; parameters fixed by `lower == upper` excluded). Its
display lists those parameters with their values and bounds.

Two more fields record what ended the run. `constraints` holds one [`ConstraintRow`](@ref)
per constraint component at `x` (empty without `constraints`; the display shows a
constraint of more than eight components, such as a keep-out per beam, as its worst row
and a count by status), and `outcome` is one of
- `:interior`: converged with every constraint met and no parameter at a bound;
- `:box`: converged with every constraint met and parameters at a bound (the box is
  active: a wider box may lower the objective);
- `:infeasible_box`: a constraint is still violated and parameters are at a bound (the
  box keeps the run from meeting it);
- `:infeasible`: a constraint is still violated and no parameter is at a bound (a local
  point where the constraints cannot be met; another start may help);
- `:iterations`, `:evaluations`, `:time_limit` or `:callback`: the budget, the time
  limit or the callback ended the run first (its constraints, if any, are in
  `constraints`).
`reason` reads `:converged` for the first four. `outer` is the number of outer
(multiplier) iterations of a constrained run, 0 without constraints; `iterations` and
`evaluations` are then totals over them, and `value` and `initial_value` are the
objective's own values, without the constraint terms.

`state` is what `optimize_pointing(...; from = result)` continues from: for a
constrained run, a named tuple of the multipliers `λ` (normalized units, one per
constraint component in the order of `constraints`), the penalty parameter `μ` and its
ceiling `μmax`, the inner runs' tolerance `tol`, the last infeasibility–complementarity
measure `Vold`, the working limits `limits` (an exact ceiling's lowered stand-in
limit; see [`irradiance_ceiling`](@ref)) and the iterate `x` they belong to (the
result's `x` unless the run returned its iterate of least violation, below); `nothing`
without constraints. A run ended by
its budget, time limit or callback keeps the state of the interrupted inner problem,
so that the restart resumes it.

**Multipliers.** The multipliers in `constraints` are not the iteration's `state.λ`
but those that best satisfy stationarity at `x`: the nonnegative least-squares
solution of `∂f/∂xᵢ + Σⱼ λⱼ ∂gⱼ/∂xᵢ = 0` over the parameters inside the box (`gⱼ` the
normalized violations, with an exact ceiling's stand-in at its working limit) and the
components within `feasibility` of their limits or beyond them, zero for the other
components (complementary). The iteration's own estimate, `λ ← max(0, λ + μ g)`
after each inner run, is exact only when the inner run ended at a stationary point; a
drive objective's inner runs end on their value tolerance instead, and at a large
penalty `μ` the update then zeroes the multipliers of active limits at which the design
ends a few 1e-5 inside (research notes, "Multipliers from the converged active set").
The fit costs one Jacobian of the objective and the constraints after the run (a
gradient's cost, not counted in `evaluations`). Where several sets of multipliers fit
equally (components with the same gradient, such as the clearances of beams related
by the ports' symmetry) the fit returns the one of least norm. `stationarity` is what
the fit leaves: the largest `|∂L/∂xᵢ|` over the parameters inside the box with
`L = f + Σⱼ λⱼ gⱼ` (objective units per parameter unit, the measure the inner runs
compare with `g_tol`; `NaN` without constraints or when phase 1 ended the run). Near
zero, the design is a Karush–Kuhn–Tucker point and the multipliers are its own; a
drive run whose inner runs ended on their value tolerance can end with it a sizable
fraction of the objective's gradient, and the multipliers and bound sensitivities are
then estimates at a design short of stationarity (on the dropped-bundle tutorial's
case, 1.2e-3 against the drive objective's largest derivative of 3.0e-3 per mm).
An exact ceiling ([`irradiance_ceiling`](@ref), [`picket_ceiling`](@ref) with `exact =
true`) reports the exact value in its row but takes part in the fit through its smooth
stand-in at the working limit: its row can read `:satisfied` (the exact value inside
the limit by more than `feasibility`) with a positive multiplier (the stand-in active).
The display notes such rows; the multiplier is the stand-in's.

**Bound sensitivities.** For each parameter in `at_bounds` (same order),
`bound_gradient` holds the objective's derivative `∂f/∂xᵢ` at `x` and
`bound_sensitivity` the Lagrangian's, `∂L/∂xᵢ` with `L = f + Σⱼ λⱼ gⱼ` (the
multipliers of `constraints`, above; without constraints the two are equal). By the
envelope theorem `bound_sensitivity` is the derivative of the optimal objective with
respect to the bound's value, the active set held: positive at a lower bound and
negative at an upper bound when the run has converged, so widening the box by δ in
that parameter (lowering a lower bound, raising an upper one) lowers the optimum by
about `|bound_sensitivity| δ` (to first order, in the objective's units per parameter
unit). `bound_gradient` is what the objective alone would gain, before the constraints
take their share. Without constraints they cost one gradient after the run; with
constraints they come from the multipliers' Jacobian (neither counted in
`evaluations`); they are meaningful at convergence. Where the active constraints and
bounds leave a flat family of designs (a drive objective under several active gates),
and where the run ended short of stationarity (its inner runs on their value
tolerance), the derivative can differ from what a wider box gains by its whole size: a
sign that would move the parameter inward shows a bound that does not hold it
(research notes, "Multipliers from the converged active set"). A restart with `from`
in the wider box measures it.

**VisRad's hard cut.** `tapered_value` is the objective that a run on VisRad's
profile tapered below its cut reached, when `optimize_pointing` optimized a scenario
with the hard cut that way (`smooth_cut`, its default); `value`, `initial_value`,
`constraints` and `outcome` are then the scenario's own, under its hard cut. `NaN`
otherwise.

**Phase 1 and relaxations** (`phase1 = true`, `relaxation = true`). `phase1` is
`nothing` when phase 1 did not run, else a named tuple of its `evaluations` and
`iterations` (included in the totals), `feasible` (it reached a design within the
limits, where the main run started), `verdict` (it converged above `feasibility`: the
run ended there, `outer` 0) and its design `x`. `relaxations` holds one
[`Relaxation`](@ref) per constraint violated or active at a verdict's design, with
`relaxation = true` (empty otherwise); their evaluations count in `evaluations`. The
display lists them under the constraints.

**Least violation and the last design within the limits** (constrained runs; `nothing`
without constraints). The run records each iterate's largest normalized violation over
the reported rows, `max(0, maxⱼ violation)` of [`ConstraintRow`](@ref) (an exact
ceiling's exact value against its limit; zero when every component holds), at the
start, at phase 1's design and after every inner run. `least_violation` is the iterate
of least violation, violations within `feasibility` counting as equal and the lower
objective then preferred, and `last_feasible` the latest iterate with every violation
within `feasibility` (`nothing` when none had), each a named tuple `(; x, violation,
value, outer)`: its parameters, its largest violation, its objective and its outer
iteration (0 for the start or phase 1's design). A run that ends with a violation above
`feasibility` (`outcome` `:infeasible` or `:infeasible_box`, or a budget, time limit or
callback stop above the limits) returns `least_violation` as `x`, its beams, value and
rows, rather than its last iterate, which can be worse than the start on the limits it
met; `state.x` is then the last iterate, the one the multipliers of `state` belong to,
where a restart with `from` continues. When an earlier iterate met the limits, the
returned design meets them and a run that ended at the penalty's ceiling reads `reason`
and `outcome` `:iterations`. Under VisRad's hard cut (`tapered_value` finite) the
records are the tapered run's.

`capped` counts the outer iterations whose inner run ended on its cap (`iterations`, or
`inner_evaluations`) instead of converging: after each that ended above the limits the
multipliers were updated and the penalty, the inner tolerance and the penalty test's
reference held ([`optimize_pointing`](@ref)); the display gives the count.
"""
struct PointingResult{B}
    x::Vector{Float64}
    beams::B
    value::Float64
    initial_value::Float64
    iterations::Int
    evaluations::Int
    stopped::Bool
    reason::Symbol
    result::Any
    names::Vector{String}
    lower::Vector{Float64}
    upper::Vector{Float64}
    at_bounds::Vector{String}
    constraints::Vector{ConstraintRow}
    outcome::Symbol
    outer::Int
    state::Union{Nothing,NamedTuple}
    bound_gradient::Vector{Float64}
    bound_sensitivity::Vector{Float64}
    phase1::Union{Nothing,NamedTuple}
    relaxations::Vector{Relaxation}
    stationarity::Float64
    tapered_value::Float64
    least_violation::Union{Nothing,NamedTuple}
    last_feasible::Union{Nothing,NamedTuple}
    capped::Int
end

_names(p::GroupPointing, n) = parameter_names(p)
_names(p, n) = hasmethod(parameter_names, Tuple{typeof(p)}) ?
               collect(String, parameter_names(p)) : ["x[$i]" for i in 1:n]

# 0: inside the box, -1 at the lower bound, 1 at the upper (within 1e-8 of the width;
# a parameter fixed by lower == upper is inside)
function _bound_side(x, lo, hi)
    lo < hi || return 0
    tol = isfinite(hi - lo) ? 1e-8 * (hi - lo) : 0.0
    return x - lo <= tol ? -1 : hi - x <= tol ? 1 : 0
end

Base.@nospecializeinfer function _pointing_result(
        s::_AbstractScenario, x, value, initial_value, iterations, evaluations, stopped, reason, result,
        lower, upper; constraints::Vector{ConstraintRow} = ConstraintRow[], outer::Int = 0,
        @nospecialize(state = nothing), @nospecialize(objective = nothing),
        gradients::Union{Nothing,NTuple{2,Vector{Float64}}} = nothing, stationarity::Float64 = NaN,
        @nospecialize(phase1 = nothing), relaxations::Vector{Relaxation} = Relaxation[],
        @nospecialize(least_violation = nothing), @nospecialize(last_feasible = nothing), capped::Int = 0)
    x = collect(Float64, x)
    names = _names(s.param, length(x))
    side = map(_bound_side, x, lower, upper)
    atb = names[side .!= 0]
    outcome = _outcome(reason, constraints, atb)
    # the derivatives at the bound parameters: from the constrained run's gradients of the
    # objective and the Lagrangian at `x`, else one gradient of the objective, only when
    # a parameter is at a bound
    gf, gL = Float64[], Float64[]
    at = findall(!=(0), side)
    if gradients !== nothing
        gf, gL = gradients[1][at], gradients[2][at]
    elseif !isempty(atb) && objective !== nothing
        gf = _scenario_gradient(s, objective, x)[at]
        gL = copy(gf)
    end
    return PointingResult(x, apply(s.param, x), Float64(value), Float64(initial_value),
                          Int(iterations), Int(evaluations), stopped, reason, result, names,
                          collect(Float64, lower), collect(Float64, upper), atb,
                          constraints, outcome, outer, state, gf, gL, phase1, relaxations, stationarity,
                          NaN, least_violation, last_feasible, capped)
end

# What ended a run: the stop's reason, or at convergence whether the constraints are met
# and parameters at a bound (`atb`, their names).
function _outcome(reason::Symbol, rows::Vector{ConstraintRow}, atb)
    reason !== :converged && return reason
    met = all(r -> r.status !== :violated, rows)
    return met ? (isempty(atb) ? :interior : :box) : (isempty(atb) ? :infeasible : :infeasible_box)
end

# The optimizer's view of a scenario on VisRad's profile with its hard cut: the same
# operator (points, normals, visibility, work plan, envelope ratios) with the cut
# tapered over `_OPT_TAPER` e-folds, the responses and caches shared; `nothing` for any
# other scenario. The taper keeps the profile's support, so the operator's points and
# the screens that skip points outside it serve both.
const _OPT_TAPER = 1.0

function _tapered(op::DepositionOperator)
    prof = op.profile
    (prof === nothing || prof.taper > 0) && return nothing
    return DepositionOperator(op.target, op.ports, op.dirs, op.pts, op.w, op.normal, op.mask, op.dest,
                              op.plan, VisRadProfile(prof.envelope, _OPT_TAPER), op.envelope)
end

_tapered(::_AbstractScenario) = nothing

function _tapered(s::Scenario)
    op = _tapered(s.deposition)
    op === nothing && return nothing
    return Scenario(s.target, s.param, op, s.model, s.observed, s.moment_weights, s.ylm_weights,
                    s.moment_response, s.ylm_response, s.lpi, s.time, getfield(s, :_response), s.albedo,
                    getfield(s, :_split))
end

# The optimizer's view of a scenario whose constraints read the element powers of sets of
# beam groups (a picket ceiling's outer cones, `_split_sets`): the same scenario, its
# evaluations keeping those powers from their own deposition (`evaluate`), which the
# constraints then read instead of depositing the groups again. `req` pairs a phase
# (`nothing`: the scenario itself) with a group set; `s` itself when none applies.
const _SplitRequest = Vector{Pair{Union{Nothing,Symbol},Vector{Symbol}}}

_split_sets(f) = _SplitRequest()

function _split_request(cs::Tuple)
    req = _SplitRequest()
    for c in cs
        append!(req, _split_sets(c.f))
    end
    return req
end

_with_split(s::_AbstractScenario, req) = s

function _with_split(s::Scenario, req)
    old = getfield(s, :_split)
    sets = copy(old)
    for (k, g) in req
        k === nothing && !(g in sets) && push!(sets, g)
    end
    length(sets) == length(old) && return s
    return Scenario(s.target, s.param, s.deposition, s.model, s.observed, s.moment_weights, s.ylm_weights,
                    s.moment_response, s.ylm_response, s.lpi, s.time, getfield(s, :_response), s.albedo, sets)
end

# A run on the tapered profile (`_tapered`) reported under the scenario `s`'s hard cut:
# the objective at the run's design and at its start `xs`, the constraint rows there
# (with the run's multipliers) and the outcome, the tapered run's value kept as
# `tapered_value`; everything else is the run's.
Base.@nospecializeinfer function _hard_cut_result(s::_AbstractScenario, r::PointingResult,
                                                  @nospecialize(objective), @nospecialize(constraints),
                                                  xs::AbstractVector, feasibility::Float64)
    ev = evaluate(s, r.x)
    value = Float64(ForwardDiff.value(objective(ev)))
    initial = Float64(ForwardDiff.value(objective(evaluate(s, collect(Float64, xs)))))
    rows = isempty(r.constraints) ? r.constraints :
           _constraint_rows(ev, Tuple(constraints), [c.multiplier for c in r.constraints], feasibility)
    return PointingResult(r.x, r.beams, value, initial, r.iterations, r.evaluations, r.stopped,
                          r.reason, r.result, r.names, r.lower, r.upper, r.at_bounds, rows,
                          _outcome(r.reason, rows, r.at_bounds), r.outer, r.state, r.bound_gradient,
                          r.bound_sensitivity, r.phase1, r.relaxations, r.stationarity, r.value,
                          r.least_violation, r.last_feasible, r.capped)
end

function Base.show(io::IO, r::PointingResult)
    print(io, "PointingResult(", r.reason, ", objective ", round(r.initial_value; sigdigits = 4),
          " → ", round(r.value; sigdigits = 4), ", ", length(r.at_bounds), " of ",
          length(r.x), " parameters at bounds")
    isempty(r.constraints) ||
        print(io, ", ", count(c -> c.status === :violated, r.constraints), " of ",
              length(r.constraints), " constraints violated")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", r::PointingResult)
    how = r.outcome in (:infeasible, :infeasible_box) ? "ended with constraints violated" :
          r.reason === :converged ? "converged" :
          r.reason === :iterations ? "reached the iteration limit" :
          r.reason === :evaluations ? "reached the evaluation limit (best design kept)" :
          r.reason === :time_limit ? "reached the time limit (best design kept)" :
          "stopped by the callback (best design kept)"
    @printf(io, "Pointing optimization %s: objective %.4g → %.4g\n", how, r.initial_value,
            r.value)
    isnan(r.tapered_value) ||
        @printf(io, "  optimized on VisRad's profile tapered below its cut (objective %.4g there); reported under the scenario's hard cut\n",
                r.tapered_value)
    @printf(io, "  %d iterations, %d evaluations, %d parameters", r.iterations,
            r.evaluations, length(r.x))
    if isempty(r.at_bounds)
        print(io, ", none at a bound")
    else
        print(io, "; at a bound (widen the box if the objective should set them):")
        w = maximum(length, r.at_bounds)
        for (j, n) in enumerate(r.at_bounds)
            i = findfirst(==(n), r.names)
            lo = _bound_side(r.x[i], r.lower[i], r.upper[i]) < 0
            @printf(io, "\n    %s = %.6g  (%s bound; box %.6g to %.6g)", rpad(n, w), r.x[i],
                    lo ? "lower" : "upper", r.lower[i], r.upper[i])
            j <= length(r.bound_sensitivity) &&
                @printf(io, "  d objective/d bound %.3g", r.bound_sensitivity[j])
        end
    end
    isempty(r.constraints) && return nothing
    what = r.outcome === :interior ? "constraints met inside the box" :
           r.outcome === :box ? "constraints met, parameters at a bound" :
           r.outcome === :infeasible_box ? "constraints violated with parameters at a bound (the box blocks them)" :
           r.outcome === :infeasible ? "constraints violated inside the box" :
           "ended by the budget, the time limit or the callback"
    p1 = r.phase1
    if p1 !== nothing && p1.verdict
        @printf(io, "\n  phase 1 (%d evaluations): no design meets every constraint from this start; the least infeasibility found:",
                p1.evaluations)
    else
        p1 === nothing ||
            @printf(io, "\n  phase 1: %s in %d evaluations", p1.feasible ? "a design within the limits" :
                    "no design within the limits", p1.evaluations)
        @printf(io, "\n  %d outer iterations", r.outer)
        r.capped > 0 && @printf(io, " (%d ended on the inner cap)", r.capped)
        print(io, ": ", what)
        isnan(r.stationarity) ||
            @printf(io, "\n  stationarity: largest |∂L/∂x| inside the box %.2e (objective units per parameter unit)",
                    r.stationarity)
    end
    lv, lf = r.least_violation, r.last_feasible
    if lv !== nothing && r.state !== nothing && haskey(r.state, :x) && r.x == lv.x && r.state.x != r.x
        @printf(io, "\n  returned the iterate of least violation (%s, largest violation %.3g), not the last (`state.x`)",
                lv.outer == 0 ? "the start or phase 1's design" : "outer iteration $(lv.outer)", lv.violation)
    end
    if lf === nothing
        r.outcome in (:interior, :box) || print(io, "\n  no iterate met every limit")
    elseif r.x != lf.x
        @printf(io, "\n  the last iterate within the limits: %s, objective %.4g (`last_feasible`)",
                lf.outer == 0 ? "the start or phase 1's design" : "outer iteration $(lf.outer)", lf.value)
    end
    _show_constraint_rows(io, r.constraints)
    # an exact ceiling's row reads the exact value (inside its limit) while its multiplier
    # is the stand-in's, active at its working limit
    stand = unique(first(split(x.name, '[')) for x in r.constraints if x.status === :satisfied && x.multiplier > 0)
    isempty(stand) ||
        print(io, "\n  ", join(stand, ", "), ": :satisfied with a positive multiplier, the smooth stand-in's at its ",
              "working limit, not the reported value's (an exact ceiling; PointingResult, \"Multipliers\")")
    isempty(r.relaxations) && return nothing
    print(io, "\n  each violated constraint relaxed alone, the others met:")
    for x in r.relaxations
        print(io, "\n    ", rpad(x.name, 12), " ", _relaxation_text(x))
    end
    return nothing
end

# A relaxation in words: the relaxed limits in the constraint's units, and as a multiple
# of the limit where the scale is the limit's magnitude.
function _relaxation_text(x::Relaxation)
    x.status === :others_conflict && return "not enough alone: the others conflict too"
    v = x.violation
    s = x.status === :stopped ? " (stopped first: the others possibly unmet)" : ""
    # relative: the scale is the limit's magnitude, so each limit moves by the same multiple
    m = x.relaxed_limit ./ x.limit
    rel = all(l -> l != 0, x.limit) &&
          all(i -> isapprox(abs(x.relaxed_limit[i] - x.limit[i]), v * abs(x.limit[i]); rtol = 1e-9),
              eachindex(x.limit))
    lims = length(x.limit) == 1 ? @sprintf("limit %.5g → %.5g", x.limit[1], x.relaxed_limit[1]) :
           "limits " * join((@sprintf("%.4g", l) for l in x.relaxed_limit), ", ")
    return @sprintf("by %.4g × its scale%s: %s%s (%d evaluations)", v,
                    rel ? @sprintf(" (the limit × %.4g)", m[1]) : "", lims, s, x.evaluations)
end

# A constraint's rows (its components) in a display: all of them up to this many, beyond
# it the worst (largest violation) and a count by status (a keep-out per beam has
# hundreds; the rows themselves stay in the vector).
const _ROWS_SHOWN = 8

function _show_constraint_rows(io::IO, rows::AbstractVector{ConstraintRow}; indent = "\n    ")
    # consecutive rows of one constraint: the name before the component's label
    stem(r) = first(split(r.name, '['))
    groups = Vector{ConstraintRow}[]
    for r in rows
        isempty(groups) || stem(groups[end][1]) != stem(r) ? push!(groups, [r]) : push!(groups[end], r)
    end
    shown = [length(g) <= _ROWS_SHOWN ? g : [g[argmax([r.violation for r in g])]] for g in groups]
    w = maximum(r -> length(r.name), Iterators.flatten(shown); init = 10)
    @printf(io, "%s%s %12s %12s %10s %10s  %s", indent, rpad("constraint", w), "value", "limit",
            "violation", "multiplier", "status")
    for (g, rs) in zip(groups, shown), r in rs
        @printf(io, "\n    %s %12.5g %s%11.5g %10.2e %10.3g  %s", rpad(r.name, w), r.value,
                r.upper ? "≤" : "≥", r.limit, r.violation, r.multiplier, r.status)
        length(g) > _ROWS_SHOWN &&
            print(io, " (the worst of ", length(g), ": ", count(x -> x.status === :violated, g), " violated, ",
                  count(x -> x.status === :active, g), " active)")
    end
end

Base.show(io::IO, ::MIME"text/plain", rows::AbstractVector{ConstraintRow}) =
    _show_constraint_rows(io, rows; indent = "    ")

struct _OptimizationStop <: Exception
    reason::Symbol
end

struct _PointingTag end

# The optimizer's dual tag precedes every other tag. ForwardDiff orders the tags of
# nested derivatives by `tagcount`, a counter read when a tag type is first compiled,
# so that a derivative taken inside the objective (its tag `Tag{typeof(g), <:Dual}`)
# nests inside the optimizer's. The precompile workload compiles this tag, so its count
# is the one from the precompile session, while the counter restarts on loading: a tag
# made at run time can then get an equal or smaller count, and its duals would not nest
# in the optimizer's ("Cannot determine ordering of Dual tags"). The gradient of
# `optimize_pointing` is always the outermost derivative (its parameters and scenario
# are Float64, and the tag is not public), so stating the order is exact.
const _PTAG = ForwardDiff.Tag{_PointingTag,Float64}
ForwardDiff.:≺(::Type{_PTAG}, ::Type{ForwardDiff.Tag{F,V}}) where {F,V} = true
ForwardDiff.:≺(::Type{ForwardDiff.Tag{F,V}}, ::Type{_PTAG}) where {F,V} = false
ForwardDiff.:≺(::Type{_PTAG}, ::Type{_PTAG}) = false

# ForwardDiff's chunk-mode `gradient!(out, f, x, cfg)` with the chunks' dual passes run
# concurrently, one task each (their inner threaded loops nest on the same pool). A
# pass is seeded exactly as ForwardDiff seeds it (the last chunk carries N partials, its
# trailing seeds zero) and every pass is deterministic, so the gradient and value are
# bitwise those of the serial `ForwardDiff.gradient!`; the value is taken from the last
# chunk, as ForwardDiff takes it. Serial on one thread, for a single chunk, or with
# `_PARALLEL_CHUNKS[] = false` (the reference path).
const _PARALLEL_CHUNKS = Ref(true)

function _pointing_gradient!(out, f, x::AbstractVector{Float64},
                             cfg::ForwardDiff.GradientConfig{T,Float64,N}) where {T,N}
    nc = cld(length(x), N)
    (nc <= 1 || Threads.nthreads() == 1 || !_PARALLEL_CHUNKS[]) &&
        return ForwardDiff.gradient!(out, f, x, cfg, Val(false))
    G = ForwardDiff.DiffResults.gradient(out)
    vals = Vector{Float64}(undef, nc)
    tasks = [Threads.@spawn _gradient_chunk!(G, vals, f, x, cfg.seeds, T, c) for c in 1:nc]
    for t in tasks
        try
            wait(t)
        catch err
            # rethrow the chunk's own exception (e.g. a stop request) unwrapped
            err isa TaskFailedException ? throw(err.task.exception) : rethrow()
        end
    end
    return ForwardDiff.DiffResults.value!(out, vals[nc])
end

# Dual pass `c` of `_pointing_gradient!`: writes only its own gradient slots and value.
function _gradient_chunk!(G, vals, f, x, seeds::NTuple{N,ForwardDiff.Partials{N,Float64}},
                          ::Type{T}, c) where {N,T}
    lo = (c - 1) * N
    zp = zero(ForwardDiff.Partials{N,Float64})
    xd = [ForwardDiff.Dual{T}(x[i], lo < i <= lo + N ? seeds[i-lo] : zp) for i in eachindex(x)]
    y = f(xd)
    for k in 1:min(N, length(x) - lo)
        G[lo+k] = y isa ForwardDiff.Dual ? ForwardDiff.partials(y, k) : 0.0
    end
    vals[c] = ForwardDiff.value(y)
    return nothing
end

"""
    report([io = stdout,] ev; lmax = 6, minor = 1e-3) -> nothing

Print a power accounting and flux-moment summary of an evaluation: the setup and
effective power per group when the scenario has a laser–plasma model; the delivered
laser power and its split over the parts (`component_power`), with the parts below
`minor` × the delivered power summed on one line; the residual, if not below that
fraction; and the observed flux's mean ([`mean_flux`](@ref)), radiation temperature and
Legendre moments `Pₗ/P₀` up to `lmax`.

The residual is the delivered power minus the power absorbed and clipped: light that
strikes no surface plus the deposition quadrature's error (`DepositionOperator`
conserves each beam's power to 1e-4 when nothing clips). It prints as `missed` when
positive and as `quadrature excess` when negative (the quadrature counts more power
than the beams deliver).
"""
function report(io::IO, ev::Evaluation; lmax::Integer = 6, minor::Real = 1e-3)
    d = ev.deposition
    P = ForwardDiff.value(d.delivered)
    if ev.scenario.lpi !== nothing
        Ps = sum(ForwardDiff.value(b.power) for b in ev.beams)
        @printf(io, "Setup laser power: %.4g TW; after the laser–plasma model:\n", Ps)
        for (g, v) in sort(collect(lpi_summary(ev.beams, ev.effective)); by = first)
            @printf(io, "  %-20s %9.4g → %9.4g TW  (×%.3f)\n", g, v.setup, v.effective, v.ratio)
        end
    end
    @printf(io, "Delivered laser power: %.4g TW\n", P)
    cp = component_power(d, ev.scenario.target)
    missed = pop!(cp, :missed)
    rows = sort(collect(cp); by = first)
    small = filter(kv -> abs(kv[2]) < minor * P, rows)
    for (k, v) in rows
        abs(v) < minor * P || @printf(io, "  %-20s %9.4g TW  (%6.2f%%)\n", k, v, 100v / P)
    end
    if !isempty(small)
        v = sum(last, small)
        @printf(io, "  %-20s %9.4g TW  (%6.2f%%)\n",
                string(length(small), length(small) == 1 ? " part" : " parts", " < ",
                       100minor, "%"), v, 100v / P)
    end
    abs(missed) < minor * P ||
        @printf(io, "  %-20s %9.4g TW  (%6.2f%%)\n", missed > 0 ? "missed" : "quadrature excess",
                missed, 100missed / P)
    if _radiative(ev.scenario)
        a = ForwardDiff.value.(flux_moments(ev; lmax))
        @printf(io, "Observed flux: mean %.4g TW/mm² (%.4g W/cm², Tr = %.1f eV)\n", a[1],
                to_Wcm2(a[1]), radiation_temperature(a[1]))
        for l in 1:lmax
            @printf(io, "  P%d/P0 = %+7.3f%%\n", l, 100 * a[l+1] / a[1])
        end
    end
    return nothing
end
report(ev::Evaluation; kwargs...) = report(stdout, ev; kwargs...)
