# Drive uncertainty from beam-to-beam errors (pointing and power), by linear propagation
# through exact per-beam sensitivities, with a Monte Carlo check, and robust objectives
# that add the expected mode power the errors cause.

struct _ErrorTag end

"""
    mode_sensitivity(scenario, beams; lmax = 4) -> NamedTuple

Sensitivity of the observed flux's real spherical-harmonic shape `c_lm/c_00`
(`mode_index` order, through `lmax ≤ 6`) to each beam's own errors: a pointing
displacement of its aim transverse to the beam (mm, along the spot's `a` and `b` axes)
and a relative power error. Returns `(; J, c, names)`: `J` is `(lmax+1)² × 3nb`,
column `3(b-1) + k` the derivative with respect to error `k` of beam `b` (`k` = 1, 2
the two transverse directions, 3 the relative power), `c` the nominal `c_lm/c_00`, and
`names` the column labels.

Each beam's deposition depends on its own errors only, so the derivatives come from
one forward-mode pass per beam with three dual partials: exact (to round-off), and
about the cost of one gradient of an objective. The beams are the scenario's setup
beams; a laser–plasma model is applied before the errors (the errors are in the
delivered beams). Errors move the aim; the beam direction (the port) is fixed.
"""
function mode_sensitivity(s::Scenario, beams::AbstractVector{<:Beam}; lmax::Integer = 4)
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    lmax <= _CACHED_YLM_LMAX || throw(ArgumentError("cached Y_lm weights reach lmax = $_CACHED_YLM_LMAX"))
    eff = apply_lpi(s.lpi, beams, s.time)
    dc, c = _beam_mode_derivatives(s, eff, s.ylm_response, (lmax + 1)^2)
    names = [string(b.name, ".", k) for b in eff for k in ("a", "b", "power")]
    return (; J = dc, c, names)
end

# Derivatives of the shape c_i/c_1 (first `nm` coefficients of the contracted response
# `R`: the scenario's `ylm_response` or `moment_response`) with respect to each beam's
# three errors, and the shape, generic in the beams' number type (nested duals for
# robust objectives).
# With `axial`, a fourth column per beam: the aim moved along the beam (its focus moves
# with it), which shared displacements (a quad's, the target's) have.
function _beam_mode_derivatives(s::Scenario, beams, R::Matrix{Float64}, nm; axial::Bool = false)
    T = mapreduce(b -> typeof(b.power), promote_type, beams)
    nb = length(beams)
    N = axial ? 4 : 3
    # the nominal coefficients, from all beams together
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

# Coefficients (first `nm`) of an element deposition through a contracted response:
# Σₑ R[:, e] Pₑ, over the elements it reaches.
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

# ∂c/∂(errors of one beam), nm × N, from a pass of that beam alone with N dual
# partials: the aim moved along the spot axes, the power scaled, and (N = 4) the aim
# moved along the beam.
function _beam_coefficient_partials(s::Scenario, b::Beam{V}, bi, R, nm, ::Val{N} = Val(3)) where {V,N}
    op = s.deposition
    D = ForwardDiff.Dual{ForwardDiff.Tag{_ErrorTag,V},V,N}
    one3(k) = D(zero(V), ForwardDiff.Partials(ntuple(j -> V(j == k), N)))
    e1, e2 = spot_axes(b)
    aim = SVector{3,D}(b.aim) + one3(1) * e1 + one3(2) * e2
    N == 4 && (aim += one3(4) * _beam_k(b))
    bd = remake(b; aim, power = b.power * (one(V) + one3(3)), defocus = D(b.defocus),
                spot = _retype(b.spot, D), dlambda = D(b.dlambda))
    # a laser–plasma model's ridge stop (`lpi_deposit`) takes g of the beam off the operator
    rs = _ridge_split(s.lpi, op, bd, s.time)
    Pr = rs === nothing ? zero(bd.power) : bd.power * rs.g
    rs === nothing || (bd = remake(bd; power = bd.power * (1 - rs.g)))
    f = BeamFrame(bd)
    # the coefficients accumulated point by point (an element vector per beam, of nested
    # duals in a robust objective, is most of the cost)
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

Expected spread of the observed drive shape from beam errors, by linear propagation
through exact per-beam sensitivities (`mode_sensitivity`).

Independent beam-to-beam errors: each beam's aim is displaced at random transverse to
the beam with rms `pointing` mm (the rms of the displacement's length, `pointing/√2`
along each axis; NIF specifies 50 µm), and its power by an rms fraction `power`.

Correlated errors (zero by default):
- `group_pointing`: every beam of a group (`groups`: a `Beam` field such as `:quad`, or
  a function of the beam, e.g. its bundle, or its cone in one hemisphere,
  `b -> (b.group, b.port[3] > 0)`) is displaced by one shared vector transverse to the
  group's mean beam direction, rms length `group_pointing` mm. Beams at an angle to that
  mean also move along themselves (their focus moves with the aim). A group whose beams
  have no common direction (`:group`, a cone over both hemispheres) is refused.
- `group_power`: one shared relative power error per group, rms `group_power`.
- `positioning`: the target is displaced from chamber center, which moves every aim by
  the same vector relative to the target: rms length `positioning` mm, isotropic, or a
  3-tuple of rms values along x, y, z.

Each source adds to the covariance independently. Returns `(; c, sigma, covariance,
rms_mode, rms_by, J)`: the nominal `c_lm/c_00`, the rms of each coefficient, their
covariance, and for each degree `l = 1:lmax` the rms of the mode's amplitude over the
errors, `√E[Σ_m (δc_lm/c_00)²]` (the expected mode-`l` budget the errors use: compare
`mode_amplitude`, which is `√(2l+1)` times the nominal one). `rms_by` splits `rms_mode`
by source (`beam`, `group`, `positioning`; they add in quadrature), and `J` is
`mode_sensitivity`'s. Linear in the errors: check large errors with `sample_modes`,
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
    return (; c, sigma, covariance = C, rms_mode = budget(C), rms_by, J)
end

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

# the unit transverse basis of a group's mean beam direction k
function _transverse_basis(k)
    u1 = normalize(cross(abs(k[3]) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0), k))
    return u1, cross(k, u1)
end

# The transverse basis of a group's shared displacement: its beams' mean direction must be
# well defined (a cone's beams from both hemispheres nearly cancel and have none).
function _group_basis(beams, idx, g)
    k = sum(_beam_k(beams[i]) for i in idx)
    norm(k) >= length(idx) / 2 ||
        throw(ArgumentError("group $(repr(g)): its beams have no common direction (|mean k| = " *
                            "$(round(norm(k) / length(idx); digits = 2))), so a shared transverse " *
                            "displacement is undefined; group by hemisphere, e.g. " *
                            "groups = b -> (b.group, b.port[3] > 0)"))
    return _transverse_basis(normalize(k))
end

# δc when beam bi's aim moves by the vector d (columns a, b and axial of J4)
function _displaced(J4, N, bi, b, d)
    e1, e2 = spot_axes(b)
    o = N * (bi - 1)
    return J4[:, o+1] .* dot(d, e1) .+ J4[:, o+2] .* dot(d, e2) .+ J4[:, o+4] .* dot(d, _beam_k(b))
end

# Shared group errors as (δc columns, variances): per group two transverse
# displacements and one power scale
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

# Target positioning: every aim moves by one vector (its x, y, z components)
function _positioning_sources(J4, N, beams, σ)
    any(>(0), σ) || return zeros(eltype(J4), size(J4, 1), 0), Float64[]
    cols = [sum(_displaced(J4, N, i, beams[i], u) for i in eachindex(beams))
            for u in (Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1))]
    return reduce(hcat, cols), [σ[i]^2 for i in 1:3]
end

"""
    sample_modes(scenario, beams; pointing = 0.05, power = 0.0, n = 100,
                 rng = Random.default_rng(), lmax = 4, group_pointing = 0.0,
                 group_power = 0.0, groups = :quad, positioning = 0.0) -> Matrix

Monte Carlo of the drive shape under the errors of `mode_uncertainty` (same keywords):
`n` samples of `c_lm/c_00` (columns), each from a full evaluation with every beam's aim
displaced along its spot axes by normal errors of rms `pointing/√2` each and its power
by a normal fraction of rms `power`, plus the shared group and target-positioning
errors. The independent check of the linear propagation, and the answer when errors
are large enough for the drive to respond nonlinearly.
"""
function sample_modes(s::Scenario, beams::AbstractVector{<:Beam}; pointing::Real = 0.05,
                      power::Real = 0.0, n::Integer = 100, rng = Random.default_rng(),
                      lmax::Integer = 4, group_pointing::Real = 0.0, group_power::Real = 0.0,
                      groups = :quad, positioning = 0.0)
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    σpos = _positioning_rms(positioning)
    nm = (lmax + 1)^2
    eff = apply_lpi(s.lpi, beams, s.time)
    keys = [_group_key(b, groups) for b in eff]
    gs = unique(keys)
    basis = Dict(g => group_pointing > 0 ? _group_basis(eff, findall(==(g), keys), g) : (zero(Vec3), zero(Vec3))
                 for g in gs)
    out = zeros(nm, n)
    for k in 1:n
        # shared errors (drawn only when present, so independent-only runs keep their sequence)
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
        c = _ylm_product(s.ylm_response, lpi_deposit(s.lpi, s.deposition, perturbed, s.time).element, nm)
        out[:, k] = c ./ c[1]
    end
    return out
end

"""
    infer_errors(scenario, beams, measured; sigma, pointing = 0.05, power = 0.0,
                 group_pointing = 0.0, group_power = 0.0, groups = :quad,
                 positioning = 0.0) -> NamedTuple

The most probable beam errors behind a measured drive-shape departure: the inverse of
`mode_uncertainty`, which takes the same error priors. A measured low mode that trends
from shot to shot (MacGowan et al., HEDP 2021, trend mode 1 of NIF implosions against
beam and target data) raises the question of which errors could have caused it; this
answers it for one
shot, linearly.

`measured` holds the departures `Δc = (c_lm/c_00)_measured − (c_lm/c_00)_nominal` as
pairs: `(l, m) => Δc` for a component, `l => [Δc_l,-l, …, Δc_l,l]` for a whole degree
(`mode_index` order). `sigma` is the measurement's rms error, one number or one per
component (in the order `measured` expands to). The nominal shape is the scenario's
for `beams`; components through `l = 6` (the cached ``Y_{lm}`` weights).

Model: `Δc = J x + ε`, with `x` the errors (zero-mean normal, the variances of
`mode_uncertainty`: per beam its aim along the spot axes `a`, `b` and its relative
power; per group a shared displacement and power; the target positioning) and `J`
their exact sensitivities (`mode_sensitivity`, with the shared columns). The posterior
mean, which is also the most probable `x`, and its covariance come from the data-space
(Woodbury) form, a solve the size of the measurement:

    x̂ = Σₚ Jᵀ (J Σₚ Jᵀ + Σₘ)⁻¹ Δc,   Σ = Σₚ − Σₚ Jᵀ (J Σₚ Jᵀ + Σₘ)⁻¹ J Σₚ.

Returns `(; beam, group, positioning, x, sigma, prior_sigma, covariance, names,
components, measured, predicted, residual, sources, resolved, beams, J, c)`:
- `beam`: `nb × 3`, each beam's inferred `a`, `b` displacement (mm) and relative power;
- `group`: a `Dict` from group key to `(; displacement, power)` (a `Vec3` in mm);
- `positioning`: the inferred common displacement of every aim relative to the target
  (the target sits at minus this from where the beams are aimed);
- `x`, `sigma`, `prior_sigma`, `covariance`, `names`: all errors as one vector (beams'
  `a`, `b`, `power`, then the groups', then x, y, z), their posterior and prior rms,
  the posterior covariance, and labels;
- `predicted = J x̂` and `residual = Δc − J x̂`, over `components` (the `(l, m)` list);
- `sources`: per error source (a beam's pointing, a beam's power, a group's pointing or
  power, the positioning), the `fraction` of `Δc` reproduced by its share of `J x̂`,
  `dot(J_s x̂_s, Δc)/|Δc|²`, sorted descending (the fractions sum to
  `dot(predicted, Δc)/|Δc|²`): the few sources that best explain the measurement;
- `resolved`: Σⱼ (1 − Σⱼⱼ/σⱼ²), the number of error combinations the data determine
  (at most the number of components);
- `beams`: `beams` with `x̂` applied, for a full `evaluate` (the nonlinear check; the
  errors are applied before a laser–plasma model, which the sensitivities follow).

The name reflects what comes back: errors, under a prior, not a unique explanation. Many
errors map to a few modes, so a handful of combinations are determined and every other
direction keeps its prior: `sigma ./ prior_sigma` stays near 1 for most errors, and
`resolved` counts what was learned. What the inverse reports is only as good as the
prior: an error left out of it (`group_pointing = 0`, say) cannot be blamed, and its
effect is attributed to the sources that are allowed. Linear in the errors.
[`infer_velocity_errors`](@ref) is the same inverse for a measured hot-spot velocity.
"""
function infer_errors(s::Scenario, beams::AbstractVector{<:Beam}, measured; sigma, kwargs...)
    comps, Δc = _measured_components(measured)
    σm = sigma isa Real ? fill(Float64(sigma), length(Δc)) : Float64.(collect(sigma))
    length(σm) == length(Δc) && all(>=(0), σm) ||
        throw(ArgumentError("sigma is one non-negative number or one per measured component"))
    return _infer_errors(s, beams, comps, Δc, Matrix(Diagonal(σm .^ 2)); kwargs...)
end

# The inverse for measured departures Δc of the components `comps`, with measurement
# covariance Σm. A `nuisance` (; columns, variances, names, label) adds error sources
# that act on the measured components directly (columns: their δc per unit error),
# after the positioning.
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
    # data space: S = A Σₚ Aᵀ + Σₘ (k × k)
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
    # the error sources: column ranges and labels
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

# The measured departures as ((l, m) list, values): pairs (l, m) => Δc or l => vector
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

# The order of `_group_sources`' columns: (group, :u1 | :u2 | :power, displacement unit)
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

# The expected squared departure the errors add to each shape coefficient (first `nm`
# of the response `R`): Σ_j σ_j² (∂(c_i/c_1)/∂e_j)², linearized.
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
