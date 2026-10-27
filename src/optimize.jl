# Pointing parameterizations, objectives, and gradient-based optimization.

const _CACHED_LMAX = 12

const _FIELDS = (:z, :r, :phi, :defocus, :power, :spot, :dlambda,
                 :theta_s, :dphi_s, :split_theta, :split_phi)
const _CYL_FIELDS = (:z, :r, :phi)
const _SPH_FIELDS = (:theta_s, :dphi_s, :split_theta, :split_phi)

"""
    GroupPointing(beams; by = :group, fields = (:z, :r, :defocus), mirror = true,
                  relative = false, radius = nothing, center = (0, 0, 0))

Maps a parameter vector to beams. Beams are grouped by `by` (`:group` — e.g. cone —,
`:quad`, `:beam`, or a function `beam -> key`); each group gets one value per field:

- `:z` — axial aim coordinate (mm). With `mirror`, beams from the lower hemisphere
  (port `z < 0`) aim at `-z`, so one value describes a top/bottom-symmetric cone.
- `:r` — radial aim offset (mm) toward the beam's own port azimuth.
- `:phi` — tangential aim offset (mm), counterclockwise about `z`.
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

`relative = true` makes the aim fields (`:z`, `:r`, `:phi`, or the spherical ones) shifts
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
end

function GroupPointing(beams::AbstractVector{<:Beam}; by = :group,
                       fields = (:z, :r, :defocus), mirror::Bool = true, relative::Bool = false,
                       radius::Union{Nothing,Real} = nothing, center = (0.0, 0.0, 0.0))
    for f in fields
        f in _FIELDS || throw(ArgumentError("unknown field $f; choose from $_FIELDS"))
    end
    by === :quad && any(b -> isempty(b.quad), beams) &&
        throw(ArgumentError("by = :quad needs quad names on every beam (OMEGA beams have none)"))
    spherical = any(in(_SPH_FIELDS), fields)
    spherical && any(in(_CYL_FIELDS), fields) &&
        throw(ArgumentError("the spherical fields $_SPH_FIELDS exclude $_CYL_FIELDS"))
    spherical && radius === nothing &&
        throw(ArgumentError("the spherical fields $_SPH_FIELDS need the sphere's `radius`"))
    radius === nothing || radius > 0 || throw(ArgumentError("radius must be positive"))
    any(in((:split_theta, :split_phi)), fields) && any(b -> isempty(b.quad), beams) &&
        throw(ArgumentError("split fields need quad names on every beam"))
    key = by === :group ? (b -> b.group) : by === :quad ? (b -> b.quad) :
          by === :beam ? (b -> b.name) : by
    ks = map(key, beams)
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
    return GroupPointing(tmpl, uk, members, collect(Symbol, fields), mirror, relative, R, c,
                         sph, sph0)
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
_radial(b::Beam) = (φ = atan(b.port[2], b.port[1]); (Vec3(cos(φ), sin(φ), 0), Vec3(-sin(φ), cos(φ), 0)))
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
        p.relative && f in (:z, :r, :phi) && continue
        if f in _SPH_FIELDS
            p.relative || (x[_pidx(p, g, fi)] = p.sph0[g][findfirst(==(f), _SPH_FIELDS)])
            continue
        end
        vals = map(mem) do bi
            b = p.template[bi]
            er, eφ = _radial(b)
            f === :z ? _zsign(p, b) * b.aim[3] :
            f === :r ? dot(b.aim, er) :
            f === :phi ? dot(b.aim, eφ) :
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
function parameters(p::GroupPointing, values::Pair...; kwargs...)
    x = initial_parameters(p)
    for (key, fv) in (values..., pairs(kwargs)...)
        g = findfirst(==(key), p.keys)
        g === nothing && throw(ArgumentError("no group $(repr(key)); the groups are " *
                                             join(repr.(p.keys), ", ")))
        for (f, v) in pairs(fv)
            fi = findfirst(==(f), p.fields)
            fi === nothing && throw(ArgumentError("$f is not a parameter; the fields are " *
                                                  join(p.fields, ", ")))
            x[_pidx(p, g, fi)] = v
        end
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
        dl = val(:dlambda)
        # spherical values as concrete types (whether each is a parameter, and its value)
        sv = (val(:theta_s), val(:dphi_s), val(:split_theta), val(:split_phi))
        has = map(!isnothing, sv)
        vals = map(v -> v === nothing ? zero(T) : v, sv)
        for bi in mem
            b = p.template[bi]
            if _spherical(p)
                aim = _sphere_aim(p, g, bi, has, vals)
            else
                er, eφ = _radial(b)
                aim0 = b.aim
                if p.relative
                    rr = r === nothing ? dot(aim0, er) : dot(aim0, er) + r
                    tt = ϕ === nothing ? dot(aim0, eφ) : dot(aim0, eφ) + ϕ
                    zz = z === nothing ? aim0[3] : aim0[3] + _zsign(p, b) * z
                else
                    rr = r === nothing ? dot(aim0, er) : r
                    tt = ϕ === nothing ? dot(aim0, eφ) : ϕ
                    zz = z === nothing ? aim0[3] : _zsign(p, b) * z
                end
                aim = rr * er + tt * eφ + SVector{3,T}(0, 0, zz)
            end
            spot = sp === nothing ? _retype(b.spot, T) :
                   _retype(with_axes(b.spot, sp * b.spot.a, sp * b.spot.b), T)
            out[bi] = remake(b; aim, spot,
                             power = pw === nothing ? T(b.power) : pw,
                             defocus = dfc === nothing ? T(b.defocus) : dfc,
                             dlambda = dl === nothing ? T(b.dlambda) : dl)
        end
    end
    return out
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
    bounds(p::GroupPointing, x0; z = 1.0, r = 1.0, phi = 1.0, defocus = 2.0,
           power = 0.5, spot = 0.5, dlambda = 2.0, theta_s = 5.0, dphi_s = 10.0,
           split_theta = 5.0, split_phi = 5.0, power_min = 0.0, power_max = Inf)
           -> (lower, upper)

Box bounds around `x0`: `±z`, `±r`, `±phi`, `±defocus` mm, relative `±power`, `±spot`
fractions, `±dlambda` Å, and `±theta_s`, `±dphi_s`, `±split_theta`, `±split_phi`
degrees (engineering defaults, not facility limits). Power bounds are further clipped to
the absolute limits `power_min ≤ P ≤ power_max` (TW per beam, e.g. a facility's per-quad
capability).
"""
function bounds(p::GroupPointing, x0 = initial_parameters(p); z = 1.0, r = 1.0, phi = 1.0,
                defocus = 2.0, power = 0.5, spot = 0.5, dlambda = 2.0, theta_s = 5.0,
                dphi_s = 10.0, split_theta = 5.0, split_phi = 5.0, power_min = 0.0,
                power_max = Inf)
    lo, hi = copy(x0), copy(x0)
    for g in eachindex(p.members), (fi, f) in enumerate(p.fields)
        i = _pidx(p, g, fi)
        if f === :power
            lo[i] = max(x0[i] * (1 - power), power_min)
            hi[i] = min(x0[i] * (1 + power), power_max)
            lo[i] <= hi[i] || throw(ArgumentError("empty power bounds for group $(p.keys[g])"))
        elseif f === :spot
            lo[i] = x0[i] * (1 - spot); hi[i] = x0[i] * (1 + spot)
        else
            δ = f === :z ? z : f === :r ? r : f === :phi ? phi : f === :dlambda ? dlambda :
                f === :theta_s ? theta_s : f === :dphi_s ? dphi_s :
                f === :split_theta ? split_theta : f === :split_phi ? split_phi : defocus
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

"""
    Scenario(target, pointing::GroupPointing; observe = nothing, radiation = true,
             time = 0.0, albedo = nothing, lpi = nothing, deposition_kw = (;),
             viewfactor_kw = (;), model = nothing)
    Scenario(target, beams; kwargs...)
    Scenario(target, parameterization; kwargs...)

The parameterization (`scenario.param`) maps the parameter vectors of `evaluate` and
`optimize_pointing` to beams. The first form builds the deposition operator for the
template beams of `pointing`; the second is `Scenario(target, GroupPointing(beams))`
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

Objectives need only the moments of the observed flux, so the scenario keeps the
response contracted with the Legendre (`lmax = 12`) and real-harmonic (`lmax = 6`)
weights: per evaluation a 13 × n and a 49 × n product, and a 62-column solve at setup.
The full observed × n response `scenario.response` (and an evaluation's `flux`) is
computed when first read; on a large (iterative) model an evaluation's `flux` is one
solve instead, until the response has been read. `model`, a `RadiationModel` of the
target, is used instead of building one from `viewfactor_kw`.
"""
struct Scenario{P}
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
end

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
    op = DepositionOperator(t, beams; deposition_kw...)
    obs = observe === nothing ? Int[] : elements(t, (observe isa Symbol ? (observe,) : observe)...)
    model = !radiation ? nothing : model === nothing ? RadiationModel(t; viewfactor_kw...) : model
    W = isempty(obs) ? zeros(0, 0) : legendre_weights(t, obs; lmax = _CACHED_LMAX)
    Y = isempty(obs) ? zeros(0, 0) : ylm_weights(t, obs; lmax = _CACHED_YLM_LMAX)
    WR, YR = zeros(0, 0), zeros(0, 0)
    if model !== nothing && !isempty(obs)
        WY = flux_response(model, obs; time, albedo, weights = vcat(W, Y))
        WR, YR = WY[1:size(W, 1), :], WY[size(W, 1)+1:end, :]
    end
    return Scenario(t, parameterization, op, model, obs, W, Y, WR, YR, lpi, Float64(time),
                    _LazyResponse(), albedo === nothing ? nothing : collect(Float64, albedo))
end

"""
    subset(scenario, keep; parameterization = GroupPointing(kept beams; same fields))

Scenario for the beams `keep` (indices into the scenario's beams) — e.g. after
dropping a bundle — reusing the view factors, flux response and mode weights; only the
deposition operator is restricted. The kept beams' template values (pointing, power)
come from the scenario's parameterization; the default new parameterization keeps its
fields, mirroring, `relative` and sphere and groups by `group` (pass one to group by quad
or beam).
"""
function subset(s::Scenario, keep::AbstractVector{<:Integer}; parameterization = nothing)
    s.param isa GroupPointing || parameterization !== nothing ||
        throw(ArgumentError("give a parameterization for a non-GroupPointing scenario"))
    beams = s.param.template[keep]
    param = parameterization === nothing ?
            GroupPointing(beams; fields = Tuple(s.param.fields), mirror = s.param.mirror,
                          relative = s.param.relative,
                          radius = _spherical(s.param) ? s.param.radius : nothing,
                          center = s.param.center) :
            parameterization
    return Scenario(s.target, param, s.deposition[keep], s.model, s.observed,
                    s.moment_weights, s.ylm_weights, s.moment_response, s.ylm_response,
                    s.lpi, s.time, getfield(s, :_response), s.albedo)
end

function Base.show(io::IO, s::Scenario)
    print(io, "Scenario(", nelements(s.target), " elements, ", length(s.deposition.ports),
          " beams, ", s.param isa GroupPointing ? nparams(s.param) : "?", " parameters, ",
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

function evaluate(s::Scenario, beams::AbstractVector{<:Beam})
    eff = apply_lpi(s.lpi, beams, s.time)
    dep = lpi_deposit(s.lpi, s.deposition, eff, s.time)
    T = eltype(dep.element)
    return Evaluation(s, beams, eff, dep, Ref{Union{Nothing,Vector{T}}}(nothing))
end

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
        return _mul(s.moment_response, ev.deposition.element)[1:lmax+1]
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
                       clip = 10.0, missed = 10.0, direct = 10.0, direct_names = (:capsule,),
                       drive = nothing, drive_weight = 1.0, center = nothing, keepout = nothing,
                       uncertainty = nothing)

Objective for indirect-drive symmetry: `Σ wₗ (aₗ/a₀)²` over `modes`, plus `odd ×` the
sum of squared odd moments (top/bottom imbalance), plus `clip`, `missed` and `direct`
times the squared fractions of laser power clipped on the outside of the target,
missing it, or hitting `direct_names` (e.g. the capsule) directly. When `drive`
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
                            direct_names = (:capsule,), drive = nothing, drive_weight = 1.0,
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
        J += _keepout_penalty(keepout, ev.beams)
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
    # parts or components (a fitted CAD capsule is the part of its pieces)
    if direct > 0 && all(n -> any(c -> c.name === n || c.part === n, ev.scenario.target.components), direct_names)
        J += direct * power_fraction(ev, direct_names...)^2
    end
    return J
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
    return _mul(ev.scenario.ylm_response, ev.deposition.element)[1:(lmax+1)^2]
end

"""
    asymmetry_objective(; modes = (1 => 1.0, 2 => 1.0, 3 => 1.0, 4 => 1.0), components = (),
                        clip = 10.0, missed = 10.0, direct = 10.0, direct_names = (:capsule,),
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
                             direct_names = (:capsule,), drive = nothing, drive_weight = 1.0,
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
        J += _keepout_penalty(keepout, ev.beams)
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
    optimize_pointing(scenario, objective; x0 = initial_parameters(scenario.param),
                      lower = bounds(scenario.param, x0)[1], upper = bounds(scenario.param, x0)[2],
                      iterations = 200, evaluations = 20 * iterations, show_trace = false,
                      g_tol = 1e-8, callback = nothing) -> PointingResult

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
than `GroupPointing` has no `bounds` method and needs `lower` and `upper`.

The run ends after `iterations` L-BFGS-B iterations or `evaluations` objective
evaluations (each gradient counts one), whichever comes first: a line search may take
many evaluations in one iteration, on a kinked objective (keep-out hinges, a coarse
mesh) dozens, so the iteration limit alone does not bound the run time. `reason` is
`:converged` (projected gradient below `g_tol`, or no further progress), `:iterations`,
`:evaluations` or `:callback`. On `:evaluations` and `:callback` the result is the best
design evaluated so far, with `stopped = true` and `result = nothing`; `iterations` is
the number of iterations completed in every case.

`callback`, if given, is called after every (non-derivative) objective evaluation with
`(; evaluations, iteration, value, best_value, best_x)`; returning `true` stops the
optimization. Use it for progress reports and cancellation.
"""
Base.@nospecializeinfer function optimize_pointing(
                           s::Scenario, @nospecialize(objective); x0 = initial_parameters(s.param),
                           lower = bounds(s.param, x0)[1], upper = bounds(s.param, x0)[2],
                           iterations::Integer = 200, evaluations::Integer = 20 * iterations,
                           show_trace::Bool = false, g_tol = 1e-8, @nospecialize(callback = nothing))
    # Inferred and compiled once for every objective and callback (each objective
    # closure type used to recompile this body and Optim's L-BFGS-B, ~1 s); the objective
    # is reached by dynamic dispatch once per evaluation.
    xs = clamp.(x0, lower, upper)
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
        return v
    end
    f(x) = seen(x, J(x))
    # One dual tag for every objective, so the compiled deposition and objective code
    # (and the package's precompiled copy) serve any objective, not only the first.
    # The tag check is off: objectives do not differentiate internally.
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
                       Optim.Options(; iterations, show_trace, g_abstol = g_tol,
                                     callback = track))
    catch err
        err isa _OptimizationStop || rethrow()
        err
    end
    if res isa _OptimizationStop
        return _pointing_result(s, best_x, best_value, initial_value, iter[], count[], true,
                                res.reason, nothing, lower, upper)
    end
    x = Optim.minimizer(res)
    n = Optim.iterations(res)
    return _pointing_result(s, x, Optim.minimum(res), initial_value, n, count[], false,
                            Optim.converged(res) || n < iterations ? :converged : :iterations,
                            res, lower, upper)
end

"""
    PointingResult

Result of [`optimize_pointing`](@ref). Fields: `x` (the best parameters), `beams`
(`apply(scenario.param, x)`), `value` and `initial_value` (the objective at `x` and at
the clamped `x0`), `iterations` (L-BFGS-B iterations completed), `evaluations`
(objective evaluations, each gradient counting one), `reason` (`:converged`,
`:iterations`, `:evaluations` or `:callback`), `stopped` (`true` for the last two, which
return the best design evaluated so far), `result` (Optim's result, `nothing` when
`stopped`), `names` (`parameter_names` of the parameterization, `x[i]` for one without
a method), the box `lower`, `upper`, and `at_bounds`: the names of the parameters at a
bound (within 1e-8 of the box width; parameters fixed by `lower == upper` excluded). Its
display lists those parameters with their values and bounds.
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

function _pointing_result(s::Scenario, x, value, initial_value, iterations, evaluations,
                          stopped, reason, result, lower, upper)
    x = collect(Float64, x)
    names = _names(s.param, length(x))
    side = map(_bound_side, x, lower, upper)
    return PointingResult(x, apply(s.param, x), Float64(value), Float64(initial_value),
                          Int(iterations), Int(evaluations), stopped, reason, result, names,
                          collect(Float64, lower), collect(Float64, upper),
                          names[side .!= 0])
end

function Base.show(io::IO, r::PointingResult)
    print(io, "PointingResult(", r.reason, ", objective ", round(r.initial_value; sigdigits = 4),
          " → ", round(r.value; sigdigits = 4), ", ", length(r.at_bounds), " of ",
          length(r.x), " parameters at bounds)")
end

function Base.show(io::IO, ::MIME"text/plain", r::PointingResult)
    how = r.reason === :converged ? "converged" :
          r.reason === :iterations ? "reached the iteration limit" :
          r.reason === :evaluations ? "reached the evaluation limit (best design kept)" :
          "stopped by the callback (best design kept)"
    @printf(io, "Pointing optimization %s: objective %.4g → %.4g\n", how, r.initial_value,
            r.value)
    @printf(io, "  %d iterations, %d evaluations, %d parameters", r.iterations,
            r.evaluations, length(r.x))
    if isempty(r.at_bounds)
        print(io, ", none at a bound")
        return nothing
    end
    print(io, "; at a bound (widen the box if the objective should set them):")
    w = maximum(length, r.at_bounds)
    for n in r.at_bounds
        i = findfirst(==(n), r.names)
        lo = _bound_side(r.x[i], r.lower[i], r.upper[i]) < 0
        @printf(io, "\n    %s = %.6g  (%s bound; box %.6g to %.6g)", rpad(n, w), r.x[i],
                lo ? "lower" : "upper", r.lower[i], r.upper[i])
    end
    return nothing
end

struct _OptimizationStop <: Exception
    reason::Symbol
end

struct _PointingTag end

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
