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
function _beam_mode_derivatives(s::Scenario, beams, R::Matrix{Float64}, nm)
    T = mapreduce(b -> typeof(b.power), promote_type, beams)
    nb = length(beams)
    # the nominal coefficients, from all beams together
    dep = deposit(s.deposition, beams)
    craw = _ylm_product(R, dep.element, nm)
    dcraw = Array{T}(undef, nm, 3, nb)
    _tforeach(nb) do bi
        dcraw[:, :, bi] = _beam_coefficient_partials(s, beams[bi], bi, R, nm)
    end
    c00 = craw[1]
    c = craw ./ c00
    # d(c/c00) = (dc - (c/c00) dc00) / c00
    J = Matrix{T}(undef, nm, 3nb)
    for bi in 1:nb, k in 1:3
        col = 3(bi - 1) + k
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

# ∂c/∂(errors of one beam), nm × 3, from a pass of that beam alone with three dual
# partials: the aim moved along the spot axes, the power scaled.
function _beam_coefficient_partials(s::Scenario, b::Beam{V}, bi, R, nm) where {V}
    op = s.deposition
    D = ForwardDiff.Dual{ForwardDiff.Tag{_ErrorTag,V},V,3}
    one3(k) = D(zero(V), ForwardDiff.Partials(ntuple(j -> V(j == k), 3)))
    e1, e2 = spot_axes(b)
    aim = SVector{3,D}(b.aim) + one3(1) * e1 + one3(2) * e2
    bd = remake(b; aim, power = b.power * (one(V) + one3(3)), defocus = D(b.defocus),
                spot = _retype(b.spot, D), dlambda = D(b.dlambda))
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
    return [ForwardDiff.partials(craw[i], k) for i in 1:nm, k in 1:3]
end

"""
    mode_uncertainty(scenario, beams; pointing = 0.05, power = 0.0, lmax = 4) -> NamedTuple

Expected spread of the observed drive shape from independent beam-to-beam errors, by
linear propagation through `mode_sensitivity`: each beam's aim is displaced at random
transverse to the beam with rms `pointing` mm (the rms of the displacement's length,
`pointing/√2` along each axis; NIF specifies 50 µm), and its power by an rms fraction
`power`. Returns `(; c, sigma, covariance, rms_mode, J)`: the nominal `c_lm/c_00`,
the rms of each coefficient, their covariance, and for each degree `l = 1:lmax` the
rms of the mode's amplitude over the errors, `√E[Σ_m (δc_lm/c_00)²]` (the expected
mode-`l` budget the errors use: compare `mode_amplitude`, which is `√(2l+1)` times
the nominal one). Linear in the errors: check large errors with `sample_modes`.
"""
function mode_uncertainty(s::Scenario, beams::AbstractVector{<:Beam}; pointing::Real = 0.05,
                          power::Real = 0.0, lmax::Integer = 4)
    (; J, c) = mode_sensitivity(s, beams; lmax)
    σ2 = _error_variances(size(J, 2) ÷ 3, pointing, power)
    C = J * (σ2 .* transpose(J))
    sigma = sqrt.(max.(diag(C), 0.0))
    rms_mode = [sqrt(sum(@view(diag(C)[l^2+1:(l+1)^2]))) for l in 1:lmax]
    return (; c, sigma, covariance = C, rms_mode, J)
end

_error_variances(nb, pointing, power) =
    repeat([pointing^2 / 2, pointing^2 / 2, power^2], nb)

"""
    sample_modes(scenario, beams; pointing = 0.05, power = 0.0, n = 100,
                 rng = Random.default_rng(), lmax = 4) -> Matrix

Monte Carlo of the drive shape under the errors of `mode_uncertainty`: `n` samples of
`c_lm/c_00` (columns), each from a full evaluation with every beam's aim displaced
along its spot axes by normal errors of rms `pointing/√2` each and its power by a
normal fraction of rms `power`. The independent check of the linear propagation, and
the answer when errors are large enough for the drive to respond nonlinearly.
"""
function sample_modes(s::Scenario, beams::AbstractVector{<:Beam}; pointing::Real = 0.05,
                      power::Real = 0.0, n::Integer = 100, rng = Random.default_rng(),
                      lmax::Integer = 4)
    _radiative(s) || throw(ArgumentError("scenario has no radiation response"))
    nm = (lmax + 1)^2
    eff = apply_lpi(s.lpi, beams, s.time)
    out = zeros(nm, n)
    for k in 1:n
        perturbed = map(eff) do b
            e1, e2 = spot_axes(b)
            remake(b; aim = Vec3(b.aim) + pointing / sqrt(2) * (randn(rng) * e1 + randn(rng) * e2),
                   power = b.power * (1 + power * randn(rng)))
        end
        c = _ylm_product(s.ylm_response, deposit(s.deposition, perturbed).element, nm)
        out[:, k] = c ./ c[1]
    end
    return out
end

# The expected squared departure the errors add to each shape coefficient (first `nm`
# of the response `R`): Σ_j σ_j² (∂(c_i/c_1)/∂e_j)², linearized.
function _error_variance(s::Scenario, beams, R, nm, (pointing, power))
    J, _ = _beam_mode_derivatives(s, beams, R, nm)
    σ2 = _error_variances(size(J, 2) ÷ 3, pointing, power)
    return [sum(σ2[j] * J[i, j]^2 for j in axes(J, 2)) for i in 1:nm]
end
