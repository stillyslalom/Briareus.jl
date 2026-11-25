# This file propagates beam-to-beam errors in pointing and power to the drive shape. The
# propagation is linear, through exact per-beam sensitivities, and a Monte Carlo
# evaluation checks it. Robust objectives add the expected mode power that the errors
# cause.

struct _ErrorTag end

"""
    mode_sensitivity(scenario, beams; lmax = 4) -> NamedTuple

Computes the sensitivity of the real spherical-harmonic shape `c_lm/c_00` of the observed
flux to the errors of each beam. The coefficients are in `mode_index` order, through
`lmax ≤ 6`. Each beam has two kinds of error. One is a pointing displacement of its aim
transverse to the beam, in mm along the spot `a` and `b` axes. The other is a relative
power error. Returns `(; J, c, names)`:

- `J` is `(lmax+1)² × 3nb`, and column `3(b-1) + k` is the derivative with respect to
  error `k` of beam `b`. Errors `k` = 1 and 2 are the two transverse directions, and
  `k` = 3 is the relative power.
- `c` is the nominal `c_lm/c_00`.
- `names` contains the column labels.

The deposition of each beam depends on its errors only, so the derivatives come from
one forward-mode pass per beam with three dual partials: exact (within round-off), and
about the cost of one gradient of an objective. The beams are the setup beams of the
scenario. Errors move the aim, and the beam direction, set by the port, is fixed.

The errors are those of the effective beams. A laser–plasma model of the scenario is
applied to the nominal setup beams first, and the errors then perturb the effective
beams after it, with the model held at its nominal state. The response of the model to
the errors is excluded: the transfer of a CBET model and the clipping of a bubble do not
change with the perturbed aims and powers. `sample_modes(...; through_lpi = true)`
measures that response.
"""
function mode_sensitivity(s::Scenario, beams::AbstractVector{<:Beam}; lmax::Integer = 4)
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    lmax <= _CACHED_YLM_LMAX || throw(ArgumentError("cached Y_lm weights reach lmax = $_CACHED_YLM_LMAX"))
    eff = apply_lpi(s.lpi, beams, s.time)
    dc, c = _beam_mode_derivatives(s, eff, s.ylm_response, (lmax + 1)^2)
    names = [string(b.name, ".", k) for b in eff for k in ("a", "b", "power")]
    return (; J = dc, c, names)
end

# Returns the derivatives of the shape c_i/c_1 with respect to the three errors of each
# beam, and the shape. The shape holds the first `nm` coefficients of the contracted
# response `R`, which is the `ylm_response` or `moment_response` of the scenario. The
# function is generic in the number type of the beams, so robust objectives can pass
# nested duals. With `axial`, each beam gets a fourth column for the aim moved along the
# beam, with the focus moving too. Shared displacements of a quad or of the target have
# this component.
function _beam_mode_derivatives(s::Scenario, beams, R::Matrix{Float64}, nm; axial::Bool = false)
    T = mapreduce(b -> typeof(b.power), promote_type, beams)
    nb = length(beams)
    N = axial ? 4 : 3
    # Compute the nominal coefficients from all beams together.
    dep = lpi_deposit(s.lpi, s.deposition, beams, s.time)
    craw = _ylm_product(R, dep.element, nm)
    dcraw = Array{T}(undef, nm, N, nb)
    _tforeach(nb) do bi
        dcraw[:, :, bi] = _beam_coefficient_partials(s, beams[bi], bi, R, nm, Val(N))
    end
    c00 = craw[1]
    c = craw ./ c00
    # d(c/c00) = (dc - (c/c00) dc00) / c00
    J = Matrix{T}(undef, nm, N * nb)
    for bi in 1:nb, k in 1:N
        col = N * (bi - 1) + k
        d00 = dcraw[1, k, bi]
        for i in 1:nm
            J[i, col] = (dcraw[i, k, bi] - c[i] * d00) / c00
        end
    end
    return J, c
end

# Returns the first `nm` coefficients of an element deposition through a contracted
# response, Σₑ R[:, e] Pₑ over the elements that the deposition reaches.
function _ylm_product(YR, P::AbstractVector{T}, nm) where {T}
    out = zeros(T, nm)
    @inbounds for e in eachindex(P)
        p = P[e]
        iszero(p) && continue
        for i in 1:nm
            out[i] += YR[i, e] * p
        end
    end
    return out
end

# Returns ∂c/∂(errors of one beam), nm × N, from a pass over only that beam with N dual
# partials. The partials move the aim along the spot axes and scale the power. With
# N = 4, the fourth partial moves the aim along the beam.
function _beam_coefficient_partials(s::Scenario, b::Beam{V}, bi, R, nm, ::Val{N} = Val(3)) where {V,N}
    op = s.deposition
    D = ForwardDiff.Dual{ForwardDiff.Tag{_ErrorTag,V},V,N}
    one3(k) = D(zero(V), ForwardDiff.Partials(ntuple(j -> V(j == k), N)))
    e1, e2 = spot_axes(b)
    aim = SVector{3,D}(b.aim) + one3(1) * e1 + one3(2) * e2
    N == 4 && (aim += one3(4) * _beam_k(b))
    bd = remake(b; aim, power = b.power * (one(V) + one3(3)), defocus = D(b.defocus),
                spot = _retype(b.spot, D), dlambda = D(b.dlambda))
    # With a ridge stop in the laser–plasma model (`lpi_deposit`), a fraction g of the beam
    # is deposited directly on the ridge elements, and the rest goes through the operator.
    rs = _ridge_split(s.lpi, op, bd, s.time)
    Pr = rs === nothing ? zero(bd.power) : bd.power * rs.g
    rs === nothing || (bd = remake(bd; power = bd.power * (1 - rs.g)))
    f = BeamFrame(bd)
    # Accumulate the coefficients point by point. An element vector per beam, of nested
    # duals in a robust objective, would be most of the cost.
    craw = zeros(D, nm)
    _sum_points(op, bi, f, _screen_frame(f, bd), _cone_slopes(op.dirs[bi], f),
                eachindex(op.w[bi]), D) do d, val
        if d > 0
            @inbounds for i in 1:nm
                craw[i] += R[i, d] * val
            end
        end
    end
    if rs !== nothing
        for (j, e) in enumerate(rs.es), i in 1:nm
            craw[i] += R[i, e] * Pr * rs.w[j]
        end
    end
    return [ForwardDiff.partials(craw[i], k) for i in 1:nm, k in 1:N]
end

_beam_k(b::Beam) = -normalize(Vec3(b.port...))

"""
    mode_uncertainty(scenario, beams; pointing = 0.05, power = 0.0, lmax = 4,
                     group_pointing = 0.0, group_power = 0.0, groups = :quad,
                     positioning = 0.0) -> NamedTuple

Computes the expected spread of the observed drive shape from beam errors. The errors are
propagated linearly through exact per-beam sensitivities (`mode_sensitivity`).

Independent beam-to-beam errors displace the aim of each beam at random transverse to the
beam, with rms `pointing` mm. That rms is of the displacement length, so each axis gets
`pointing/√2`; NIF specifies 50 µm. The power of each beam is scaled by an rms fraction
`power`.

Correlated errors are zero by default:

- `group_pointing` displaces every beam of a group by one shared vector transverse to the
  mean beam direction of the group, with rms length `group_pointing` mm. `groups` is a
  `Beam` field such as `:quad` or a function of the beam, which can select bundles or the
  cone in one hemisphere, `b -> (b.group, b.port[3] > 0)`. Beams at an angle to that mean
  also move along their axes, with the focus moving with the aim. A group whose beams have
  no common direction, such as `:group` for a cone over both hemispheres, is refused.
- `group_power` gives one shared relative power error per group, with rms `group_power`.
- `positioning` displaces the target from chamber center, which moves every aim by the
  same vector relative to the target. It is an isotropic rms length in mm or a 3-tuple of
  rms values along x, y and z.

Each source adds to the covariance independently. Returns
`(; c, sigma, covariance, rms_mode, floor, rms_by, J)`:

- `c` is the nominal `c_lm/c_00`, `sigma` the rms of each coefficient and `covariance`
  their covariance.
- `rms_mode[l]`, for each degree `l = 1:lmax`, is the rms of the mode amplitude over the
  errors, `√E[Σ_m (δc_lm/c_00)²]`. It is the expected mode-`l` budget used by the errors.
  `mode_amplitude` gives `√(2l+1)` times the nominal amplitude for comparison.
- `floor[l]` is `√(2l+1) × rms_mode[l]`, the same spread in the units of
  [`mode_departure`](@ref), which is `|Δaₗ/a₀|` for an axisymmetric change. It is the
  noise floor of a drive shape under these errors. [`mode_gate`](@ref) takes its gate
  from this result, with `fraction` a multiple of the floor, one number or one per degree.
- `rms_by` splits `rms_mode` by source into `beam`, `group` and `positioning`, which add in
  quadrature.
- `J` is the matrix of `mode_sensitivity`.

The propagation is linear in the errors. Large errors can be checked with `sample_modes`,
which takes the same keywords.
"""
function mode_uncertainty(s::Scenario, beams::AbstractVector{<:Beam}; pointing::Real = 0.05,
                          power::Real = 0.0, lmax::Integer = 4, group_pointing::Real = 0.0,
                          group_power::Real = 0.0, groups = :quad, positioning = 0.0)
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    lmax <= _CACHED_YLM_LMAX || throw(ArgumentError("cached Y_lm weights reach lmax = $_CACHED_YLM_LMAX"))
    min(pointing, power, group_pointing, group_power) >= 0 ||
        throw(ArgumentError("error rms values must be non-negative"))
    σpos = _positioning_rms(positioning)
    nm = (lmax + 1)^2
    eff = apply_lpi(s.lpi, beams, s.time)
    shared = group_pointing > 0 || group_power > 0 || any(>(0), σpos)
    N = shared ? 4 : 3
    J4, c = _beam_mode_derivatives(s, eff, s.ylm_response, nm; axial = shared)
    J = shared ? J4[:, [N * (bi - 1) + k for bi in eachindex(eff) for k in 1:3]] : J4
    cov(V, σ2) = isempty(σ2) ? zeros(nm, nm) : V * (σ2 .* transpose(V))
    Cb = cov(J, _error_variances(length(eff), pointing, power))
    Cg = cov(_group_sources(J4, N, eff, groups, group_pointing, group_power)...)
    Cp = cov(_positioning_sources(J4, N, eff, σpos)...)
    C = Cb + Cg + Cp
    sigma = sqrt.(max.(diag(C), 0.0))
    budget(C) = [sqrt(max(sum(@view(diag(C)[l^2+1:(l+1)^2])), 0.0)) for l in 1:lmax]
    rms_by = (; beam = budget(Cb), group = budget(Cg), positioning = budget(Cp))
    rms = budget(C)
    return (; c, sigma, covariance = C, rms_mode = rms, floor = _mode_floor(rms), rms_by, J)
end

# Converts the rms of each mode-l budget to the noise floor per degree in the units of
# `mode_departure`, √(2l+1) times that rms.
_mode_floor(rms_mode::AbstractVector{<:Real}) = [sqrt(2l + 1) * r for (l, r) in enumerate(rms_mode)]

_error_variances(nb, pointing, power) =
    repeat([pointing^2 / 2, pointing^2 / 2, power^2], nb)

function _positioning_rms(p)
    σ = p isa Real ? ntuple(_ -> p / sqrt(3), 3) : Tuple(Float64.(p))
    length(σ) == 3 && all(>=(0), σ) ||
        throw(ArgumentError("positioning is an rms length or three non-negative rms values (x, y, z)"))
    return σ
end

_group_key(b::Beam, groups::Symbol) = getfield(b, groups)
_group_key(b::Beam, groups) = groups(b)

# Returns a unit transverse basis for the mean beam direction k of a group.
function _transverse_basis(k)
    u1 = normalize(cross(abs(k[3]) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0), k))
    return u1, cross(k, u1)
end

# Builds the transverse basis of the shared displacement of a group. The mean direction
# of its beams must be well defined; the beams of a cone from both hemispheres nearly
# cancel and have none.
function _group_basis(beams, idx, g)
    k = sum(_beam_k(beams[i]) for i in idx)
    norm(k) >= length(idx) / 2 ||
        throw(ArgumentError("group $(repr(g)): its beams have no common direction (|mean k| = " *
                            "$(round(norm(k) / length(idx); digits = 2))), so a shared transverse " *
                            "displacement is undefined; group by hemisphere, e.g. " *
                            "groups = b -> (b.group, b.port[3] > 0)"))
    return _transverse_basis(normalize(k))
end

# Returns δc when the aim of beam bi moves by the vector d, from the a, b and axial
# columns of J4.
function _displaced(J4, N, bi, b, d)
    e1, e2 = spot_axes(b)
    o = N * (bi - 1)
    return J4[:, o+1] .* dot(d, e1) .+ J4[:, o+2] .* dot(d, e2) .+ J4[:, o+4] .* dot(d, _beam_k(b))
end

# Returns the shared group errors as δc columns and their variances. Each group has two
# transverse displacements and one power scale.
function _group_sources(J4, N, beams, groups, pointing, power)
    cols = Vector{eltype(J4)}[]
    σ2 = Float64[]
    keys = [_group_key(b, groups) for b in beams]
    for g in unique(keys)
        idx = findall(==(g), keys)
        if pointing > 0
            for u in _group_basis(beams, idx, g)
                push!(cols, sum(_displaced(J4, N, i, beams[i], u) for i in idx))
                push!(σ2, pointing^2 / 2)
            end
        end
        if power > 0
            push!(cols, sum(J4[:, N*(i-1)+3] for i in idx))
            push!(σ2, power^2)
        end
    end
    return isempty(cols) ? zeros(eltype(J4), size(J4, 1), 0) : reduce(hcat, cols), σ2
end

# Returns the target positioning errors, which move every aim by one vector, as columns
# for its x, y and z components.
function _positioning_sources(J4, N, beams, σ)
    any(>(0), σ) || return zeros(eltype(J4), size(J4, 1), 0), Float64[]
    cols = [sum(_displaced(J4, N, i, beams[i], u) for i in eachindex(beams))
            for u in (Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1))]
    return reduce(hcat, cols), [σ[i]^2 for i in 1:3]
end

"""
    sample_modes(scenario, beams; pointing = 0.05, power = 0.0, n = 100,
                 rng = Random.default_rng(), lmax = 4, group_pointing = 0.0,
                 group_power = 0.0, groups = :quad, positioning = 0.0,
                 through_lpi = false) -> Matrix

Runs a Monte Carlo of the drive shape under the errors of `mode_uncertainty`, with the
same keywords. Returns `n` samples of `c_lm/c_00` as columns. Each sample comes from a
full evaluation with the shared group and target-positioning errors and with the
independent errors of every beam. The aim of each beam is displaced along its spot axes
by normal errors of rms `pointing/√2` each, and its power is scaled by a normal fraction
of rms `power`. The samples check the linear propagation independently. They also give
the answer when the errors are large enough for the drive to respond nonlinearly.

`through_lpi` selects the beams that carry the errors when the scenario has a
laser–plasma model. By default the errors perturb the effective beams after the model,
as in `mode_sensitivity`, and the samples check the linear propagation. With
`through_lpi = true` the errors perturb the setup beams, which are the beams that the
laser delivers, and the model is applied to every sample. The spread then includes the
response of the model to the errors, such as the transfer of a CBET model and the
clipping of a bubble at the perturbed aims and powers. The two settings draw the same
random errors, so the difference of their samples is that response. Without a model
they are identical.
"""
function sample_modes(s::Scenario, beams::AbstractVector{<:Beam}; pointing::Real = 0.05,
                      power::Real = 0.0, n::Integer = 100, rng = Random.default_rng(),
                      lmax::Integer = 4, group_pointing::Real = 0.0, group_power::Real = 0.0,
                      groups = :quad, positioning = 0.0, through_lpi::Bool = false)
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    σpos = _positioning_rms(positioning)
    nm = (lmax + 1)^2
    # the beams that carry the errors: the setup beams, or the effective beams
    eff = through_lpi ? beams : apply_lpi(s.lpi, beams, s.time)
    keys = [_group_key(b, groups) for b in eff]
    gs = unique(keys)
    basis = Dict(g => group_pointing > 0 ? _group_basis(eff, findall(==(g), keys), g) : (zero(Vec3), zero(Vec3))
                 for g in gs)
    out = zeros(nm, n)
    for k in 1:n
        # Draw the shared errors only when present, so that runs with only independent
        # errors keep the same random sequence.
        D = any(>(0), σpos) ? Vec3(σpos[1] * randn(rng), σpos[2] * randn(rng), σpos[3] * randn(rng)) : zero(Vec3)
        gd = Dict(g => group_pointing > 0 ?
                       group_pointing / sqrt(2) * (randn(rng) * basis[g][1] + randn(rng) * basis[g][2]) : zero(Vec3)
                  for g in gs)
        gp = Dict(g => group_power > 0 ? group_power * randn(rng) : 0.0 for g in gs)
        perturbed = map(eachindex(eff)) do i
            b = eff[i]
            e1, e2 = spot_axes(b)
            remake(b; aim = Vec3(b.aim) + pointing / sqrt(2) * (randn(rng) * e1 + randn(rng) * e2) +
                            gd[keys[i]] + D,
                   power = b.power * (1 + power * randn(rng)) * (1 + gp[keys[i]]))
        end
        through_lpi && (perturbed = apply_lpi(s.lpi, perturbed, s.time))
        c = _ylm_product(s.ylm_response, lpi_deposit(s.lpi, s.deposition, perturbed, s.time).element, nm)
        out[:, k] = c ./ c[1]
    end
    return out
end

"""
    infer_errors(scenario, beams, measured; sigma, pointing = 0.05, power = 0.0,
                 group_pointing = 0.0, group_power = 0.0, groups = :quad,
                 positioning = 0.0) -> NamedTuple

Infers the most probable beam errors behind a measured drive-shape departure. This
function is the inverse of `mode_uncertainty` and takes the same error priors. A measured
low mode can trend from shot to shot; MacGowan et al. (HEDP 2021) trend mode 1 of NIF
implosions against beam and target data. Such a trend raises the question of which errors
could have caused it, and `infer_errors` gives a linear answer for one shot.

`measured` contains the departures `Δc = (c_lm/c_00)_measured − (c_lm/c_00)_nominal` as
pairs. A pair `(l, m) => Δc` gives one component, and `l => [Δc_l,-l, …, Δc_l,l]` gives a
whole degree in `mode_index` order. `sigma` is the rms error of the measurement, one
number or one per component in the order that `measured` expands to. The nominal shape
is that of the scenario for `beams`. Components are allowed through `l = 6`, the limit of
the cached ``Y_{lm}`` weights.

The model is `Δc = J x + ε`, with `x` the errors and `J` their exact sensitivities
(`mode_sensitivity`, with the shared columns). The errors are zero-mean normal with the
variances of `mode_uncertainty`: per beam, its aim along the spot axes `a`, `b` and its
relative power; per group, a shared displacement and power; and the target positioning.
The posterior mean, which is also the most probable `x`, and its covariance come from the
data-space (Woodbury) form, a solve the size of the measurement:

    x̂ = Σₚ Jᵀ (J Σₚ Jᵀ + Σₘ)⁻¹ Δc,   Σ = Σₚ − Σₚ Jᵀ (J Σₚ Jᵀ + Σₘ)⁻¹ J Σₚ.

Returns `(; beam, group, positioning, x, sigma, prior_sigma, covariance, names,
components, measured, predicted, residual, sources, resolved, beams, J, c)`:
- `beam`: `nb × 3`, the inferred `a`, `b` displacement (mm) and relative power of each
  beam;
- `group`: a `Dict` from group key to `(; displacement, power)`, with the displacement a
  `Vec3` in mm;
- `positioning`: the inferred common displacement of every aim relative to the target.
  The target therefore lies at minus this vector from where the beams are aimed.
- `x`, `sigma`, `prior_sigma`, `covariance`, `names`: all errors as one vector, with
  their posterior and prior rms, the posterior covariance and labels. The vector holds
  the beam `a`, `b` and `power`, then the group errors, then x, y and z.
- `predicted = J x̂` and `residual = Δc − J x̂`, over the `(l, m)` list `components`;
- `sources`: the `fraction` of `Δc` reproduced by the share of `J x̂` from each error
  source, `dot(J_s x̂_s, Δc)/|Δc|²`, sorted descending. A source is the pointing or power
  of a beam or of a group, or the positioning. The fractions sum to
  `dot(predicted, Δc)/|Δc|²`, and the first entries are the sources that best explain
  the measurement.
- `resolved`: Σⱼ (1 − Σⱼⱼ/σⱼ²), the number of error combinations that the data
  determine, at most the number of components;
- `beams`: `beams` with `x̂` applied, for a full `evaluate` as the nonlinear check. The
  errors are applied before a laser–plasma model, as in the sensitivities.

The name describes the result: the errors that are most probable under a prior, which are
not a unique explanation. Many errors map to a few modes, so a handful of combinations are
determined and every other direction stays at its prior. `sigma ./ prior_sigma` therefore
stays near 1 for most errors, and `resolved` gives the number of determined combinations.
The result is only as good as the prior. An error left out of it, such as
`group_pointing = 0`, cannot be blamed, and its effect is attributed to the sources that are
allowed. The inference is linear in the errors. [`infer_velocity_errors`](@ref) is the same
inverse for a measured hot-spot velocity.
"""
function infer_errors(s::Scenario, beams::AbstractVector{<:Beam}, measured; sigma, kwargs...)
    comps, Δc = _measured_components(measured)
    σm = sigma isa Real ? fill(Float64(sigma), length(Δc)) : Float64.(collect(sigma))
    length(σm) == length(Δc) && all(>=(0), σm) ||
        throw(ArgumentError("sigma is one non-negative number or one per measured component"))
    return _infer_errors(s, beams, comps, Δc, Matrix(Diagonal(σm .^ 2)); kwargs...)
end

# Computes the inverse for measured departures Δc of the components `comps`, with
# measurement covariance Σm. A `nuisance` of the form (; columns, variances, names, label)
# adds error sources after the positioning. These sources act on the measured components
# directly, and `columns` holds their δc per unit error.
function _infer_errors(s::Scenario, beams::AbstractVector{<:Beam}, comps, Δc, Σm;
                       pointing::Real = 0.05, power::Real = 0.0, group_pointing::Real = 0.0,
                       group_power::Real = 0.0, groups = :quad, positioning = 0.0,
                       nuisance = nothing)
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    min(pointing, power, group_pointing, group_power) >= 0 ||
        throw(ArgumentError("error rms values must be non-negative"))
    lmax = maximum(first, comps)
    lmax <= _CACHED_YLM_LMAX || throw(ArgumentError("cached Y_lm weights reach lmax = $_CACHED_YLM_LMAX"))
    σpos = _positioning_rms(positioning)
    nm = (lmax + 1)^2
    eff = apply_lpi(s.lpi, beams, s.time)
    nb = length(eff)
    shared = group_pointing > 0 || group_power > 0 || any(>(0), σpos)
    N = shared ? 4 : 3
    J4, c = _beam_mode_derivatives(s, eff, s.ylm_response, nm; axial = shared)
    Jb = shared ? J4[:, [N * (bi - 1) + k for bi in 1:nb for k in 1:3]] : J4
    Vg, σg = _group_sources(J4, N, eff, groups, group_pointing, group_power)
    Vp, σp = _positioning_sources(J4, N, eff, σpos)
    layout = _group_layout(eff, groups, group_pointing, group_power)
    rows = [mode_index(l, m) for (l, m) in comps]
    A = hcat(Jb, Vg, Vp)[rows, :]
    v = vcat(_error_variances(nb, pointing, power), σg, σp)
    if nuisance !== nothing
        A = hcat(A, nuisance.columns)
        v = vcat(v, nuisance.variances)
    end
    # Form the data-space matrix S = A Σₚ Aᵀ + Σₘ, of size k × k.
    F = cholesky(Symmetric(A * (v .* transpose(A)) + Σm); check = false)
    issuccess(F) || throw(ArgumentError("the measured components are not constrained: " *
                                        "give a measurement sigma or a prior that moves them"))
    K = transpose(F \ (A .* transpose(v)))                   # Σₚ Aᵀ S⁻¹ (n × k)
    x = K * Δc
    post = Diagonal(v) - K * (A .* transpose(v))
    post = (post + transpose(post)) / 2
    sd = sqrt.(max.(diag(post), 0.0))
    resolved = sum((1 - post[j, j] / v[j] for j in eachindex(v) if v[j] > 0); init = 0.0)
    predicted = A * x
    # Collect the column ranges and labels of the error sources.
    srcs = Tuple{String,UnitRange{Int}}[]
    for (bi, b) in enumerate(eff)
        push!(srcs, (string(b.name, " pointing"), 3bi-2:3bi-1), (string(b.name, " power"), 3bi:3bi))
    end
    o = 3nb
    names = [string(b.name, ".", k) for b in eff for k in ("a", "b", "power")]
    gd = Dict(g => (; displacement = zero(Vec3), power = 0.0) for g in unique(_group_key(b, groups) for b in eff))
    for (g, kind, u) in layout
        o += 1
        if kind === :power
            push!(srcs, (string(g, " power"), o:o))
            push!(names, string(g, ".power"))
            gd[g] = (; gd[g]..., power = x[o])
        else
            kind === :u1 && push!(srcs, (string(g, " pointing"), o:o+1))
            push!(names, string(g, ".", kind))
            gd[g] = (; gd[g]..., displacement = gd[g].displacement + x[o] * u)
        end
    end
    D = zero(Vec3)
    if !isempty(σp)
        push!(srcs, ("positioning", o+1:o+3))
        append!(names, ("positioning.x", "positioning.y", "positioning.z"))
        D = Vec3(x[o+1], x[o+2], x[o+3])
        o += 3
    end
    if nuisance !== nothing
        push!(srcs, (nuisance.label, o+1:o+length(nuisance.variances)))
        append!(names, nuisance.names)
    end
    n2 = dot(Δc, Δc)
    sources = [(; source, fraction = n2 > 0 ? dot(A[:, r] * x[r], Δc) / n2 : 0.0, delta = A[:, r] * x[r])
               for (source, r) in srcs if any(>(0), @view v[r])]
    sort!(sources; by = t -> -t.fraction)
    beam = permutedims(reshape(x[1:3nb], 3, nb))
    keys = [_group_key(b, groups) for b in eff]
    applied = map(eachindex(beams)) do i
        b = beams[i]
        e1, e2 = spot_axes(eff[i])
        g = gd[keys[i]]
        remake(b; aim = Vec3(b.aim) + beam[i, 1] * e1 + beam[i, 2] * e2 + g.displacement + D,
               power = b.power * (1 + beam[i, 3]) * (1 + g.power))
    end
    return (; beam, group = gd, positioning = D, x, sigma = sd, prior_sigma = sqrt.(v),
            covariance = post, names, components = comps, measured = Δc, predicted,
            residual = Δc - predicted, sources, resolved, beams = applied, J = A, c = c[rows])
end

# Expands the measured departures, pairs (l, m) => Δc or l => vector, into an (l, m)
# list and its values.
function _measured_components(measured)
    comps = Tuple{Int,Int}[]
    vals = Float64[]
    for kv in measured
        key, val = first(kv), last(kv)
        if key isa Tuple
            l, m = key
            1 <= l && abs(m) <= l || throw(ArgumentError("measured component $key is not (l ≥ 1, |m| ≤ l)"))
            push!(comps, (l, m)); push!(vals, val)
        else
            key >= 1 && length(val) == 2key + 1 ||
                throw(ArgumentError("measured degree $key takes 2l + 1 = $(2key + 1) values"))
            append!(comps, [(key, m) for m in -key:key]); append!(vals, val)
        end
    end
    isempty(comps) && throw(ArgumentError("no measured components"))
    allunique(comps) || throw(ArgumentError("a measured component is given twice"))
    return comps, vals
end

# Lists the columns of `_group_sources` in order, each as (group, :u1 | :u2 | :power,
# displacement unit).
function _group_layout(beams, groups, pointing, power)
    out = Tuple{Any,Symbol,Vec3}[]
    keys = [_group_key(b, groups) for b in beams]
    for g in unique(keys)
        idx = findall(==(g), keys)
        if pointing > 0
            u1, u2 = _group_basis(beams, idx, g)
            push!(out, (g, :u1, u1), (g, :u2, u2))
        end
        power > 0 && push!(out, (g, :power, zero(Vec3)))
    end
    return out
end

# Returns the expected squared departure that the errors add to each of the first `nm`
# shape coefficients of the response `R`, linearized as Σ_j σ_j² (∂(c_i/c_1)/∂e_j)².
function _error_variance(s::Scenario, beams, R, nm, e)
    shared = e.group_pointing > 0 || e.group_power > 0 || any(>(0), e.positioning)
    N = shared ? 4 : 3
    J, _ = _beam_mode_derivatives(s, beams, R, nm; axial = shared)
    σ2 = _error_variances(length(beams), e.pointing, e.power)
    v = [sum(σ2[3(bi-1)+k] * J[i, N*(bi-1)+k]^2 for bi in eachindex(beams) for k in 1:3) for i in 1:nm]
    shared || return v
    for (V, σ) in (_group_sources(J, N, beams, e.groups, e.group_pointing, e.group_power),
                   _positioning_sources(J, N, beams, e.positioning))
        for j in eachindex(σ), i in 1:nm
            v[i] += σ[j] * V[i, j]^2
        end
    end
    return v
end
