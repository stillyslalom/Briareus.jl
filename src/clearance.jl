# This file checks beam clearance. It builds beam envelopes at an intensity contour and
# computes the worst-case hardware clearance over pointing and focus tolerances and the
# worst-case power and intensity on hardware.

"""
    enclosed_radius(order, fraction)

Returns the radius, in units of the 1/e semi-axes, of the super-Gaussian contour that
encloses `fraction` of the power. The fraction inside `ρ` is `P(2/n, ρⁿ)`, the
regularized lower incomplete gamma function. For spots with a halo use `envelope_radius`.
"""
enclosed_radius(n::Real, f::Real) = gamma_inc_inv(2 / n, f, 1 - f)^(1 / n)

"""
    envelope_radius(spot, power; intensity = nothing, peak_fraction = nothing,
                    enclosed = nothing) -> ρ

Returns the best-focus contour radius in units of the core 1/e semi-axes, with the halo
included. Exactly one of three keywords sets the contour: an absolute `intensity`
(W/cm²), a fraction `peak_fraction` of the peak focal intensity, or the `enclosed` power
fraction. Returns 0 if the whole spot is below
an intensity threshold.
"""
function envelope_radius(sp::SuperGaussianSpot, P; intensity = nothing, peak_fraction = nothing,
                         enclosed = nothing)
    count(!isnothing, (intensity, peak_fraction, enclosed)) == 1 ||
        throw(ArgumentError("give exactly one of intensity, peak_fraction, enclosed"))
    sp = _value_spot(sp)
    P = ForwardDiff.value(P)
    if enclosed !== nothing
        0 < enclosed < 1 || throw(ArgumentError("enclosed must be in (0, 1)"))
        g = ρ -> enclosed_fraction(sp, ρ) - enclosed            # increasing
    else
        thr = intensity !== nothing ? intensity / 1e14 : peak_fraction * focal_profile(sp, P, 0.0)
        focal_profile(sp, P, 0.0) <= thr && return 0.0
        g = ρ -> thr - focal_profile(sp, P, ρ)                  # increasing
    end
    lo, hi = 0.0, 1.0
    while g(hi) < 0
        hi *= 2
        hi > 1e6 && throw(ArgumentError("envelope contour not found"))
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        g(mid) < 0 ? (lo = mid) : (hi = mid)
        hi - lo < 1e-12 * hi && break
    end
    return hi
end

"""
    BeamEnvelope(beam; intensity, peak_fraction, enclosed, margin = 0.0)

A geometric beam envelope. At axial distance `s` from best focus, the envelope is the
best-focus contour ellipse (see `envelope_radius`) grown by the f-cone radius `|s|/(2F)`
and a uniform `margin` (mm). In geometric optics with a phase plate, every lens point
fills the whole spot. The intensity at a point is then the mean of the focal intensity
over a disk of radius `|s|/(2F)`, so outside the envelope it is below the contour level
at every `s`.
Clearance is the signed distance to the envelope; negative means inside.
"""
struct BeamEnvelope
    focus::Vec3
    aim::Vec3
    k::Vec3
    e1::Vec3
    e2::Vec3
    a::Float64
    b::Float64
    cone::Float64          # tan(half-angle)
    margin::Float64
end

function BeamEnvelope(b::Beam; intensity = nothing, peak_fraction = nothing,
                      enclosed = nothing, margin = 0.0)
    if all(isnothing, (intensity, peak_fraction, enclosed))
        peak_fraction = 1e-4
    end
    ρ = envelope_radius(b.spot, b.power; intensity, peak_fraction, enclosed)
    e1, e2 = spot_axes(b)
    k = direction(b)
    aim = Vec3(ForwardDiff.value.(b.aim))
    focus = aim + ForwardDiff.value(b.defocus) * k
    return BeamEnvelope(focus, aim, k, e1, e2, ρ * ForwardDiff.value(b.spot.a),
                        ρ * ForwardDiff.value(b.spot.b), 1 / (2b.fnumber), Float64(margin))
end

"""
    envelope_distance(env, x) -> (distance, s)

Returns the signed distance (mm) from point `x` to the beam envelope in the plane through
`x` normal to the beam, and the axial position `s` of that plane relative to best focus.
"""
function envelope_distance(env::BeamEnvelope, x::Vec3)
    d = x - env.focus
    s = dot(d, env.k)
    u, v = dot(d, env.e1), dot(d, env.e2)
    return _ellipse_distance(u, v, env.a, env.b) - env.cone * abs(s) - env.margin, s
end

# Returns the signed distance from (u, v) to the ellipse with semi-axes (a, b), negative
# inside. The distance is the projection of (u, v) minus the closest point onto the
# outward normal at that point (`_ellipse_closest`). It is valid inside and outside, to a
# few units in the last place also next to the major axis, as checked against a BigFloat
# minimum (research notes, "Optimization defaults decided"). A degenerate ellipse, with
# a = 0 or b = 0, is the point at the origin.
function _ellipse_distance(u::Float64, v::Float64, a::Float64, b::Float64)
    (a == 0 || b == 0) && return hypot(u, v)
    X0, X1, _ = _ellipse_closest(u, v, a, b)
    n1, n2 = _ellipse_unit_normal(u, v, a, b, X0, X1)
    return n1 * (u - X0) + n2 * (v - X1)
end

# Returns the point (X0, X1) of the ellipse (u/a)² + (v/b)² = 1 closest to (u, v), and
# whether (u, v) lies on the medial axis. The medial axis is the part of the major axis,
# inside the ellipse, between its centers of curvature. A point on it has two closest
# points, each the mirror image of the other across the major axis. The function uses
# Eberly's robust point–ellipse method, with the Lagrange multiplier found by Newton's
# method, and returns the center for a degenerate ellipse. Only a point on the major axis
# takes the on-axis branch. A point next to the axis keeps its closest point, which moves
# fast along the ellipse next to a center of curvature. `HardwareKeepOut` marks the minima
# next to the medial axis by their mirror point (`_medial_tie`).
function _ellipse_closest(u::Float64, v::Float64, a::Float64, b::Float64)
    (a == 0 || b == 0) && return 0.0, 0.0, false
    swap = a < b
    e0, e1, y0, y1 = swap ? (b, a, abs(v), abs(u)) : (a, b, abs(u), abs(v))
    medial = false
    if y1 > 0 && y1 / e1 > 0                                  # (z1 = 0 only by underflow)
        if y0 > 0
            z0, z1 = y0 / e0, y1 / e1
            if z0^2 + z1^2 == 1
                x0, x1 = y0, y1
            else
                r0 = (e0 / e1)^2
                q = (e0 - e1) * (e0 + e1) / e1^2               # r0 − 1
                n0 = r0 * z0
                # Find the root σ = s + 1 of F = (n0/(σ + r0 − 1))² + (z1/σ)² − 1, which is
                # decreasing and convex for σ > 0, so Newton's method from F(z1) ≥ 0 climbs
                # to it. The iteration is in σ. Inside, next to the major axis, σ ~ y1, and
                # iterating in the multiplier s = σ − 1 would give σ only to an absolute
                # 1e-16.
                σ = z1
                for _ in 1:200
                    α, β = n0 / (σ + q), z1 / σ
                    F = α^2 + β^2 - 1
                    F <= 0 && break
                    σn = σ + F / (2 * (α^2 / (σ + q) + β^2 / σ))
                    σn <= σ && break
                    σ = σn
                end
                x0, x1 = r0 * y0 / (σ + q), y1 / σ
            end
        else
            x0, x1 = 0.0, e1
        end
    else
        numer, denom = e0 * y0, e0^2 - e1^2
        if numer < denom
            ξ = numer / denom
            x0, x1 = e0 * ξ, e1 * sqrt(1 - ξ^2)
            medial = true
        else
            x0, x1 = e0, 0.0
        end
    end
    X0, X1 = swap ? (x1, x0) : (x0, x1)
    return copysign(X0, u), copysign(X1, v), medial
end

# Returns the outward unit normal of the ellipse at its point (X0, X1) closest to (u, v).
# For a degenerate ellipse, whose closest point is its center, it returns the direction
# of (u, v).
function _ellipse_unit_normal(u::Float64, v::Float64, a::Float64, b::Float64, X0, X1)
    if a == 0 || b == 0
        r = hypot(u, v)
        return r > 0 ? (u / r, v / r) : (1.0, 0.0)
    end
    n1, n2 = X0 / a^2, X1 / b^2
    m = hypot(n1, n2)
    return n1 / m, n2 / m
end

"""
    SpotWing

Supertype of low-level circular far-field components outside the phase-plate core of a
beam. Only clearance checks use them, since they carry ≲1% of the energy and do not
matter for drive symmetry. Each wing is centered on the beam axis at best focus, and
radii are in mm.
"""
abstract type SpotWing end

"""
    PowerLawWing(fraction; rh = :r95, rmax = 1.5)

Halo with `I ∝ [1 + (r/rh)²]^(-3/2)` for `r ≤ rmax`, carrying `fraction` of the beam power.
`rh = :r95` uses the radius enclosing 95% of the core, the geometric mean of the axes.
`rmax` is the spatial-filter cutoff, ≈ ±200 µrad × 7.7 m ≈ 1.5 mm on NIF.
"""
struct PowerLawWing <: SpotWing
    fraction::Float64
    rh::Union{Float64,Symbol}
    rmax::Float64
end
PowerLawWing(fraction; rh = :r95, rmax = 1.5) = PowerLawWing(fraction, rh isa Symbol ? rh : Float64(rh), rmax)

"""
    UniformScatter(fraction; radius = 23.0)

A far-angle scatter floor, with `fraction` of the beam power spread uniformly over a disk
of `radius` mm at the target plane. The default covers NIF final-optics scatter up to
≈3 mrad × 7.7 m.
"""
struct UniformScatter <: SpotWing
    fraction::Float64
    radius::Float64
end
UniformScatter(fraction; radius = 23.0) = UniformScatter(fraction, radius)

"""
    nif_wings(; halo = 0.01, rmax = 1.5, scatter = 1e-3, scatter_radius = 23.0)

Returns the default NIF clearance wing model. Its values are engineering defaults (see
reference/research-notes.md), and no facility criterion sets them; override them with
measured far-field data. The model has two wings:

- a 1% power-law halo, truncated at the spatial-filter angle. OMEGA measured < 1%
  outside the super-Gaussian, and the NIF/LMJ design specs allowed 2% outside the
  nominal spot.
- a 10⁻³ uniform scatter floor out to 3 mrad.
"""
nif_wings(; halo = 0.01, rmax = 1.5, scatter = 1e-3, scatter_radius = 23.0) =
    SpotWing[PowerLawWing(halo; rmax), UniformScatter(scatter; radius = scatter_radius)]

_wing_rh(w::PowerLawWing, b::Beam) =
    w.rh === :r95 ? envelope_radius(b.spot, 1.0; enclosed = 0.95) * sqrt(ForwardDiff.value(b.spot.a) * ForwardDiff.value(b.spot.b)) :
    w.rh isa Symbol ? throw(ArgumentError("unknown rh $(w.rh)")) : w.rh

# Returns the intensity (TW/mm²) of a wing at transverse radius r for beam power P.
function _wing_intensity(w::PowerLawWing, rh, P, r)
    r > w.rmax && return 0.0
    norm = 2π * rh^2 * (1 - 1 / sqrt(1 + (w.rmax / rh)^2))
    return w.fraction * P / norm * (1 + (r / rh)^2)^(-1.5)
end
_wing_intensity(w::UniformScatter, _, P, r) = r > w.radius ? 0.0 : w.fraction * P / (π * w.radius^2)

# Returns the radius (mm) of the contour of a wing at intensity thr (TW/mm²), or 0 if the
# wing is never that bright.
function _wing_contour(w::PowerLawWing, rh, P, thr)
    I0 = _wing_intensity(w, rh, P, 0.0)
    I0 <= thr && return 0.0
    return min(w.rmax, rh * sqrt((I0 / thr)^(2 / 3) - 1))
end
_wing_contour(w::UniformScatter, _, P, thr) = _wing_intensity(w, nothing, P, 0.0) > thr ? w.radius : 0.0

# Returns the cumulative power fraction of a wing inside radius r.
_wing_cdf(w::PowerLawWing, rh, r) = (1 - 1 / sqrt(1 + (min(r, w.rmax) / rh)^2)) / (1 - 1 / sqrt(1 + (w.rmax / rh)^2))
_wing_cdf(w::UniformScatter, _, r) = min(r, w.radius)^2 / w.radius^2
_wing_rmax(w::PowerLawWing) = w.rmax
_wing_rmax(w::UniformScatter) = w.radius

# Ray traces a circular wing on a polar grid of equal-power annuli into per-component
# power.
function _raytrace_wing!(L, t::Target, b::Beam, w::SpotWing, rh; nr = 48, nθ = 64, ncone = 7)
    val = ForwardDiff.value
    e1, e2 = spot_axes(b)
    k = direction(b)
    focus = Vec3(val.(b.aim)) + val(b.defocus) * k
    dirs, dw = cone_directions(b, ncone)
    far = 4 * (t.scale + norm(focus)) + 1
    P = val(b.power) * w.fraction
    rmax = _wing_rmax(w)
    # Place the annulus edges at equal power and sample each annulus at its median radius.
    cdf_inv(u) = begin
        lo, hi = 0.0, rmax
        for _ in 1:60
            mid = (lo + hi) / 2
            _wing_cdf(w, rh, mid) < u ? (lo = mid) : (hi = mid)
        end
        (lo + hi) / 2
    end
    for i in 1:nr
        r = cdf_inv((i - 0.5) / nr)
        for j in 1:nθ
            ϕ = 2π * (j - 0.5) / nθ
            pf = focus + r * (cos(ϕ) * e1 + sin(ϕ) * e2)
            for (d, wd) in zip(dirs, dw)
                _, tri = closest_hit(t.bvh, pf - far * d, d)
                tri == 0 && continue
                L[t.tri_comp[tri]] += P * wd / (nr * nθ)
            end
        end
    end
    return L
end

"""
    ClearanceRow

One pair of a beam and a hardware component. "Worst" values are taken over the error set
of pointing, target offset, hardware tolerance and focus:

- `clearance` — worst-case signed distance to the beam envelope (mm), negative for
  overlap
- `nominal` — the same at nominal pointing and focus
- `point`, `s` — hardware point of minimum clearance and its axial position from best
  focus (mm)
- `power`, `intercepted` — worst-case power (TW) on the component (core + wings) and
  its fraction of the beam power
- `peak_intensity` — worst-case beam intensity (W/cm², normal to the beam, core +
  wings) at the upstream surface points of the component, ignoring shadowing
- `status` — `:fail`, `:review` or `:pass` against the thresholds of the report
"""
struct ClearanceRow
    beam::String
    component::Symbol
    clearance::Float64
    nominal::Float64
    point::Vec3
    s::Float64
    power::Float64
    intercepted::Float64
    peak_intensity::Float64
    status::Symbol
end

"""
    ClearanceReport

The result of `clearance`. It holds one `ClearanceRow` per beam and hardware component
and the envelope, error-set and threshold settings; see `violations` and
`min_clearance`. `canceled` is true when a `clearance` callback stopped the check early.
`rows` then hold only the beams completed before the cancel, which are the leading beams
in input order.
"""
struct ClearanceReport
    rows::Vector{ClearanceRow}
    contour::String
    transverse::Float64      # total transverse error radius (mm)
    errors::NamedTuple
    margin::Float64
    targets::Vector{Symbol}
    thresholds::NamedTuple
    canceled::Bool
end
ClearanceReport(rows, contour, transverse, errors, margin, targets, thresholds) =
    ClearanceReport(rows, contour, transverse, errors, margin, targets, thresholds, false)

"""
    DEFAULT_CLEARANCE_THRESHOLDS

The status tiers of `clearance`. They are engineering defaults, and no facility
criterion sets them:

- `:fail` above 1e12 W/cm², for plasma formation, or above 1e-3 of the beam power on a
  component;
- `:review` above 1e9 W/cm², for the onset of ablation and debris, or above 1e-5.
"""
const DEFAULT_CLEARANCE_THRESHOLDS = (fail_intensity = 1e12, fail_fraction = 1e-3,
                                      review_intensity = 1e9, review_fraction = 1e-5)

function _status(I, frac, th)
    (I > th.fail_intensity || frac > th.fail_fraction) && return :fail
    (I > th.review_intensity || frac > th.review_fraction) && return :review
    return :pass
end

"""
    clearance(target, beams; targets, hardware = all other surfaces,
              intensity = nothing, peak_fraction = 1e-4, enclosed = nothing,
              wings = nif_wings(), pointing = 0.05, target_offset = 0.025,
              hardware_tolerance = 0.0, focus = 1.0, margin = 0.0, ssd = 0.0,
              thresholds = DEFAULT_CLEARANCE_THRESHOLDS, nangles = 16, level = 2,
              resolution = 1e-4, rtol = 0.01,
              op = nothing, nspot = 61, callback = nothing,
              batch = Threads.nthreads()) -> ClearanceReport

Checks beam-to-hardware clearance, judged on the worst case over an error set. `targets`
and `hardware` name parts or components (see `elements`). A part split for computation
stands for all its pieces and gets one row per beam, the worst case over its pieces.

**Spot.** The spot is the `SuperGaussianSpot` core of each beam, with its optional
super-Gaussian halo, plus clearance-only circular `wings`. The default wings,
`nif_wings()`, are a 1% power-law halo and a 10⁻³ scatter floor. They are engineering
defaults, to be overridden with measured far fields. `ssd` (mm) is the half-width of the
SSD smear along the dispersion axis, added isotropically to the envelope as a
conservative choice.

**Envelope.** The envelope is the best-focus contour grown by the f-cone, `ssd` and
`margin`. The contour is set by one of an absolute `intensity` (W/cm²), a `peak_fraction` of
the core peak, or an `enclosed` core fraction. With wings, the contour of each spot
component is taken at the threshold divided by the number of components, and the envelope is
their union, so the total intensity outside it is below the threshold.

**Error set.** Transverse errors add linearly, not in quadrature, as
`R = pointing + target_offset + hardware_tolerance` (mm). The three terms are:

- the per-beam pointing, with a NIF spec of 50 µm rms; pass ≈ 0.12 for a 3σ statistical
  radius;
- the rigid target offset, with a NIF mean beams-to-target offset of ≈ 25 µm;
- the as-built hardware tolerance, from metrology, or 0.1–0.25 mm without it.

`focus` is the ± focal-position error in mm, ≤ 1 mm in the 1994 NIF spec. The worst-case
clearance is exact, nominal − `R` − `focus/(2F)`. Worst-case power and intensity are
maxima over aims displaced to `nangles` points on the radius-`R` circle, plus the
nominal aim, at focus offsets `0, ±focus`.

**Status.** `thresholds` are in W/cm² and in fractions of the beam power. A pair is
`:fail` above `fail_intensity` (1e12, plasma formation) or `fail_fraction` (1e-3). It is
`:review` above `review_intensity` (1e9, onset of ablation and debris) or
`review_fraction` (1e-5). These are engineering defaults, and no facility review
criterion sets them; `talis_check` applies the NIF rules.

**Reach.** Only hardware upstream of the landing point of each beam is checked. The
landing point is where the beam axis first meets one of the `targets` components, the
surfaces that the beam should hit, or the aim plane if the axis misses them. Hardware is
split into `4^level` sub-triangles per triangle. A point counts only if beam light can
reach it past the opaque surfaces of the target. Apertures are transparent, and the
triangles of the component itself never block it. Beam light reaches a point in one of
three cases:

- the straight segment to the point from the beam axis at the same axial position is
  unobstructed;
- the point is outside the envelope, and the segment to it from the nearest envelope
  point along the same transverse ray is unobstructed, provided that envelope point is
  lit. An envelope point is lit when some f-cone direction reaches it from upstream, so
  envelope behind a surface that the beam edge already struck is dark.
- the point is inside the envelope and is itself lit.

Hardware outside a hohlraum wall is therefore not judged against a beam inside it, while
hardware beside the beam path before an entrance hole still is. A sub-triangle counts
whole when any of its vertices or its centroid is reached. This choice is conservative,
since a sub-triangle is dropped only when every path to each of those points is blocked.

**Minimum.** The clearance is minimized over the whole of each counted sub-triangle,
and not only at its vertices. The envelope distance is Lipschitz in the point, with
constant √(1 + cone²), so a branch-and-bound subdivision bounds the minimum from below.
The bound is within `resolution` (mm) or `rtol` times its distance from zero worst-case
clearance, whichever is larger. A plate crossed by a beam between its vertices is
therefore found. The reported clearance is that lower bound. The peak intensity is the
largest over the reached vertices and the point of least clearance. The power on
hardware comes from the reference ray trace, with an `nspot` grid per core component and
equal-power polar grids for the wings. When `op`, a `DepositionOperator` for these ports,
is given, the core power comes from it.

**Progress and cancel.** With a `callback`, beams are checked in batches of `batch` in
input order. `callback((; done, total))` is called with the number of beams completed,
once with `done = 0` before any work and then after each batch. Returning `true`
cancels the check. `clearance` then returns at once with the rows of the completed beams
and `canceled = true` in the report. A `true` returned after the last batch is
ignored. The
callback always runs on the calling task, never on a worker thread, so it may touch
task-local state; it must not block for long, since no work runs while it does. Rows
do not depend on `callback` or `batch`.
"""
function clearance(t::Target, beams::AbstractVector{<:Beam}; targets,
                   hardware = nothing, intensity = nothing, peak_fraction = nothing,
                   enclosed = nothing, wings = nif_wings(), pointing = 0.05,
                   target_offset = 0.025, hardware_tolerance = 0.0, focus = 1.0,
                   margin = 0.0, ssd = 0.0, thresholds = DEFAULT_CLEARANCE_THRESHOLDS,
                   nangles::Integer = 16, level::Integer = 2, resolution = 1e-4, rtol = 0.01,
                   op::Union{Nothing,DepositionOperator} = nothing, nspot::Integer = 61,
                   callback = nothing, batch::Integer = Threads.nthreads())
    all(isnothing, (intensity, peak_fraction, enclosed)) && (peak_fraction = 1e-4)
    contour = intensity !== nothing ? @sprintf("%.3g W/cm²", intensity) :
              peak_fraction !== nothing ? @sprintf("%.3g of core peak", peak_fraction) :
              @sprintf("%.6g of core enclosed", enclosed)
    R = pointing + target_offset + hardware_tolerance
    tgt = Symbol[targets...]
    tgt_idx = Set(ci for n in tgt for ci in part_pieces(t, n))
    # List the hardware as parts, with all their pieces, or components, one row each.
    hw = hardware === nothing ?
         unique!([c.part for (ci, c) in enumerate(t.components) if c.kind === :surface && !(ci in tgt_idx)]) :
         Symbol[hardware...]
    hw_idx = [hardware === nothing ? setdiff(part_pieces(t, n), tgt_idx) : part_pieces(t, n)
              for n in hw]
    resolution > 0 || throw(ArgumentError("resolution must be > 0"))
    0 <= rtol < 1 || throw(ArgumentError("rtol must be in [0, 1)"))
    # Collect the sub-triangles of each hardware item, which are the roots of the distance
    # search, and their vertices.
    samples = Dict{Int,Tuple{Vector{NTuple{3,Vec3}},Vector{Vec3}}}()
    for (j, cis) in enumerate(hw_idx)
        pts = Vec3[]
        for k in eachindex(t.tri_area)
            t.tri_comp[k] in cis || continue
            v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
            _subtri_vertices!(pts, v0, v0 + e1, v0 + e2, level)
        end
        samples[j] = ([(pts[i], pts[i+1], pts[i+2]) for i in 1:3:length(pts)], unique(pts))
    end
    errs = (; pointing = Float64(pointing), target_offset = Float64(target_offset),
            hardware_tolerance = Float64(hardware_tolerance), focus = Float64(focus),
            ssd = Float64(ssd))
    finish(rows, canceled) = ClearanceReport(rows, contour, R, errs, Float64(margin), tgt,
                                             thresholds, canceled)
    rowsof(bs) = _clearance_rows(t, bs, samples, hw, hw_idx, tgt_idx, wings, op, nspot, R,
                                 focus, nangles, thresholds, (Float64(resolution), Float64(rtol));
                                 intensity, peak_fraction, enclosed, margin = margin + ssd)
    callback === nothing && return finish(rowsof(beams), false)

    # Check the beams in batches, reporting progress and polling for a cancel between them.
    nb = length(beams)
    rows = ClearanceRow[]
    callback((; done = 0, total = nb)) === true && nb > 0 && return finish(rows, true)
    step = max(1, batch)
    for lo in 1:step:nb
        hi = min(lo + step - 1, nb)
        append!(rows, rowsof(view(beams, lo:hi)))
        callback((; done = hi, total = nb)) === true && hi < nb && return finish(rows, true)
    end
    return finish(rows, false)
end

# Returns the rows of the clearance check for `beams`, which are all the beams of a call
# without batches.
function _clearance_rows(t::Target, beams, samples, hw, hw_idx, tgt_idx, wings, op, nspot, R,
                         focus, nangles, thresholds, acc; intensity, peak_fraction, enclosed,
                         margin)
    nc = length(t.components)
    # Build the beam variants over the error set, with the nominal variant first for each
    # beam.
    nb = length(beams)
    variants = Beam{Float64}[]
    owner = Int[]
    for (bi, b) in enumerate(beams)
        bv = _value_beam(b)
        e1, e2 = spot_axes(bv)
        shifts = R > 0 ?
                 vcat([zero(Vec3)], [R * (cos(2π * j / nangles) * e1 + sin(2π * j / nangles) * e2)
                                     for j in 0:nangles-1]) : [zero(Vec3)]
        for δ in shifts, df in (focus > 0 ? (0.0, -focus, focus) : (0.0,))
            push!(variants, remake(bv; aim = Vec3(bv.aim) + δ, defocus = bv.defocus + df))
            push!(owner, bi)
        end
    end
    Pc = op === nothing ? raytrace_component_power(t, variants; nspot) :
         beam_component_power(op, variants)
    if !isempty(wings)
        Threads.@threads for v in eachindex(variants)
            b = variants[v]
            L = zeros(nc)
            for w in wings
                _raytrace_wing!(L, t, b, w, w isa PowerLawWing ? _wing_rh(w, b) : nothing)
            end
            Pc[v, :] .+= L
        end
    end

    rows = Vector{ClearanceRow}(undef, nb * length(hw_idx))
    wt = Tuple(wings)        # concrete wing types in the per-sample loop
    Threads.@threads for bi in 1:nb
        vids = findall(==(bi), owner)
        b = variants[vids[1]]
        envs, wingc = _envelopes(b, wings; intensity, peak_fraction, enclosed, margin)
        s_stop = _landing(t, envs[1], tgt_idx)
        frames = [BeamFrame(variants[v]) for v in vids]
        cone = cone_directions(b, _CLEARANCE_NCONE)[1]
        rhs = map(w -> w isa PowerLawWing ? _wing_rh(w, b) : 0.0, wt)
        P = b.power
        shift = R + focus * envs[1].cone
        for (j, cis) in enumerate(hw_idx)
            # The worst clearance, best − shift, is judged against 0.
            best, bp, bs, Imax = _scan_hardware(t, cis, samples[j]..., envs[1], wingc, s_stop,
                                                frames, cone, wt, rhs, P, (acc..., shift))
            worst = best - shift
            pw = maximum(sum(ci -> Pc[v, ci], cis; init = 0.0) for v in vids)
            frac = P > 0 ? pw / P : 0.0
            Iw = to_Wcm2(Imax)
            rows[(bi-1)*length(hw_idx)+j] = ClearanceRow(b.name, hw[j], worst, best,
                                                         bp, bs, pw, frac, Iw,
                                                         _status(Iw, frac, thresholds))
        end
    end
    return rows
end

const _CLEARANCE_NCONE = 37

# Checks the hardware components `ci` against the envelope upstream of `s_stop`. Returns
# the smallest distance over the reachable sub-triangles `leaves`, with its point and
# axial position. The distance is a lower bound with accuracy `acc` (see
# `_lipschitz_minimum`). Also returns the largest intensity over the error-set variants
# `frames` at the reachable vertices `verts` and at that point. A sub-triangle is
# reachable if one of its vertices or its centroid is reachable (`_reachable`).
function _scan_hardware(t::Target, ci, leaves, verts, env::BeamEnvelope, wingc, s_stop,
                        frames, cone, wings::Tuple, rhs, P, acc)
    Imax = 0.0
    # Compute the transverse slopes of the cone directions, which all variants share
    # through the port.
    f1 = frames[1]
    αmax = maximum(dj -> abs(dot(dj, f1.e1) / dot(dj, f1.k)), cone)
    βmax = maximum(dj -> abs(dot(dj, f1.e2) / dot(dj, f1.k)), cone)
    reach = Dict{Vec3,Bool}()
    for x in verts
        d, s = _envelope_union_distance((env,), wingc, x)
        s <= s_stop || continue
        reach[x] = _reachable(t, ci, env, cone, x, d, s)
        reach[x] && (Imax = max(Imax, _variant_intensity(frames, αmax, βmax, cone, wings, rhs, P, x)))
    end
    dist = x -> _envelope_union_distance((env,), wingc, x)[1]
    axial = x -> dot(x - env.focus, env.k)
    admit = function (leaf)
        any(x -> get(reach, x, false), leaf) && return true
        m = sum(leaf) / 3
        d, s = _envelope_union_distance((env,), wingc, m)
        return s <= s_stop && _reachable(t, ci, env, cone, m, d, s)
    end
    best, bp = _lipschitz_minimum(dist, leaves, sqrt(1 + env.cone^2), acc, zero(Vec3), axial, s_stop;
                                  admit)
    isfinite(best) || return best, bp, NaN, Imax
    Imax = max(Imax, _variant_intensity(frames, αmax, βmax, cone, wings, rhs, P, bp))
    return best, bp, axial(bp), Imax
end

# Returns the largest intensity at `x` over the error-set variants `frames`. The core and
# halo enter through the f-cone, and the wings are added.
function _variant_intensity(frames, αmax, βmax, cone, wings::Tuple, rhs, P, x)
    Imax = 0.0
    for f in frames
        I = _beyond(f, αmax, βmax, x) ? 0.0 : _cone_mean_intensity(f, cone, x)
        if !isempty(wings)
            dv = x - f.focus
            r = hypot(dot(dv, f.e1), dot(dv, f.e2))
            I += _wing_sum(wings, rhs, P, r)
        end
        Imax = max(Imax, I)
    end
    return Imax
end

# Checks whether beam light can reach the hardware sample `x` of components `ci` at axial
# position `s` and signed envelope distance `d`. The triangles of that component are
# ignored, and apertures are not in the BVH. The sample is reachable in these cases:
# - The segment to `x` from the beam axis at the same `s` is unobstructed.
# - Outside the envelope (d > 0), the envelope point `p` nearest `x` along the transverse
#   ray from the axis is lit, and the segment from `p` to `x` is unobstructed. A point is
#   lit when it sees upstream along some f-cone direction, so envelope beyond a surface
#   that the beam edge already struck is dark.
# - Inside the envelope (d ≤ 0), `x` itself is lit.
# A point is skipped only when every one of these paths is blocked.
function _reachable(t::Target, ci, env::BeamEnvelope, cone, x, d, s)
    skip = let comp = t.tri_comp, ci = ci
        id -> comp[id] in ci
    end
    ε = _ray_eps(t)
    a = env.focus + s * env.k
    _unblocked(t, a, x, skip, ε) && return true
    d > 0 || return _lit(t, x, cone, skip, ε)
    ρ = norm(x - a)
    p = a + (clamp(ρ - d, 0.0, ρ) / ρ) * (x - a)
    return _unblocked(t, p, x, skip, ε) && _lit(t, p, cone, skip, ε)
end

# Checks whether some f-cone direction reaches `x` from upstream without crossing an
# opaque triangle.
function _lit(t::Target, x, cone, skip, ε)
    for dj in cone
        any_hit(t.bvh, x - ε * dj, -dj, Inf, skip) || return true
    end
    return false
end

# Checks that no opaque triangle, apart from the `skip`ped ones, lies on the open segment
# from `a` to `x`.
function _unblocked(t::Target, a, x, skip, ε)
    L = norm(x - a)
    L <= 2ε && return true
    return !any_hit(t.bvh, a, (x - a) / L, L - ε, skip)
end

_wing_sum(wings::Tuple, rhs::Tuple, P, r) =
    sum(map((w, rh) -> _wing_intensity(w, rh, P, r), wings, rhs); init = 0.0)

# Returns the beam intensity at `x`, the mean over cone directions of the best-focus
# profile where each ray through `x` crosses the focal plane. This is the f-cone
# convolution that deposition uses. The moment-matched `intensity` is not used, because
# it sharpens the edge of a defocused flat-topped spot, where clearance is judged.
function _cone_mean_intensity(f::BeamFrame, ds, x)
    d = x - f.focus
    s = dot(d, f.k)
    acc = 0.0
    for dj in ds
        sj = s / dot(dj, f.k)
        xf = d - sj * dj
        u, v = dot(xf, f.e1), dot(xf, f.e2)
        acc += _sg_focal(f.core, u, v)
        f.has_halo && (acc += _sg_focal(f.halo, u, v))
    end
    return acc / length(ds)
end

# Returns the core envelope, the contour at thr/N for N spot components, and the wing
# contour radii.
function _envelopes(b::Beam, wings; intensity, peak_fraction, enclosed, margin)
    N = 1 + length(wings)
    P = b.power
    if enclosed !== nothing
        core = BeamEnvelope(b; enclosed, margin)
        thr = focal_profile(b.spot, P, enclosed_radius(b.spot.order, enclosed))
    else
        thr = intensity !== nothing ? intensity / 1e14 : peak_fraction * focal_profile(b.spot, P, 0.0)
        core = BeamEnvelope(b; intensity = 1e14 * thr / N, margin)
    end
    wingc = [_wing_contour(w, w isa PowerLawWing ? _wing_rh(w, b) : nothing, P, thr / N) for w in wings]
    return [core], wingc
end

# Returns the signed distance to the union of the core envelope and the wing disks. Each
# wing disk is grown by the f-cone and by the margin of the core envelope.
function _envelope_union_distance(envs, wingc, x)
    env = envs[1]
    d, s = envelope_distance(env, x)
    if !isempty(wingc)
        dv = x - env.focus
        r = hypot(dot(dv, env.e1), dot(dv, env.e2))
        for rw in wingc
            rw > 0 || continue
            d = min(d, r - rw - env.cone * abs(s) - env.margin)
        end
    end
    return d, s
end

# Returns the axial position, relative to best focus, where the beam axis first meets a
# target component. If the axis never meets one, it returns the aim plane, or `miss`
# when given.
_landing(t::Target, env::BeamEnvelope, tgt_idx; miss = nothing) = first(_landing_hit(t, env, tgt_idx; miss))

# Returns `_landing` and the component met, or 0 if none.
function _landing_hit(t::Target, env::BeamEnvelope, tgt_idx; miss = nothing)
    far = 4 * (t.scale + norm(env.aim)) + 1
    o = env.aim - far * env.k
    traveled = 0.0
    while true
        τ, tri = closest_hit(t.bvh, o, env.k)
        tri == 0 && break
        x = o + τ * env.k
        t.tri_comp[tri] in tgt_idx && return Float64(dot(x - env.focus, env.k)), Int(t.tri_comp[tri])
        o = x + _ray_eps(t) * env.k
        traveled += τ
        traveled > 2far && break
    end
    return Float64(miss === nothing ? dot(env.aim - env.focus, env.k) : miss), 0
end

function _subtri_vertices!(pts, a, b, c, level)
    if level == 0
        push!(pts, a, b, c)
        return
    end
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    _subtri_vertices!(pts, a, ab, ca, level - 1)
    _subtri_vertices!(pts, ab, b, bc, level - 1)
    _subtri_vertices!(pts, ca, bc, c, level - 1)
    _subtri_vertices!(pts, ab, bc, ca, level - 1)
end

# Returns a lower bound of the minimum of `f` over the simplices `roots` and the point of
# the smallest value found, or `x0` if none is found. The simplices are segments or
# triangles given as tuples of points. With `acc = (ε, rtol, thr)`, the bound is within
# max(ε, rtol·|value − thr|). The bound never crosses the threshold `thr` that the value
# is judged against, so precision is spent only near it.
#
# `f` is Lipschitz with constant `lip`, so on a simplex with centroid m and vertices
# within ρ of it, f ≥ f(m) − lip·ρ. The search is best first. The simplex of lowest bound
# is bisected until that bound is within the accuracy of the smallest value found, or
# until lip·ρ is. Only points with `axial(x) ≤ s_stop` count. Because `axial` is linear,
# a simplex whose vertices all lie beyond `s_stop` lies beyond it entirely. With `admit`,
# a root counts only if `admit(root)` holds. The test is made when the root first comes
# up, so only roots that could hold the minimum are tested. Values found on roots not yet
# admitted are not used. The result is `(Inf, x0)` when no point counts.
function _lipschitz_minimum(f, roots::AbstractVector{S}, lip, acc, x0, axial = nothing,
                            s_stop = Inf; admit = nothing) where {S}
    ε0, rtol, thr = acc
    inside = x -> axial === nothing || axial(x) <= s_stop
    ub, xb = Inf, x0
    if admit === nothing
        for r in roots, x in r
            inside(x) || continue
            v = f(x)
            v < ub && ((ub, xb) = (v, x))
        end
    end
    heap = Tuple{Float64,Float64,S,Bool}[]      # (bound, lip·ρ, simplex, admitted), a min-heap
    for r in roots
        ub, xb = _simplex_bound!(heap, f, inside, lip, r, ub, xb, admit === nothing)
    end
    lb = Inf
    while !isempty(heap)
        low, lρ, sx, ok = _heap_pop!(heap)
        ε = isfinite(ub) ? max(ε0, rtol * abs(ub - thr)) : ε0
        if low >= ub - ε                        # every simplex left is bounded by `low`
            lb = min(lb, low)
            break
        elseif !ok                              # a root not yet admitted
            (admit === nothing || admit(sx)) || continue
            ub, xb = _simplex_bound!(heap, f, inside, lip, sx, ub, xb, true)
            for x in sx
                inside(x) || continue
                v = f(x)
                v < ub && ((ub, xb) = (v, x))
            end
        elseif lρ <= max(ε0, rtol * abs(low + lρ - thr))    # relative to f(m) = low + lρ
            lb = min(lb, low)
            for x in sx                         # a settled simplex across s_stop
                inside(x) || continue
                v = f(x)
                v < ub && ((ub, xb) = (v, x))
            end
        else
            for c in _simplex_split(sx)
                ub, xb = _simplex_bound!(heap, f, inside, lip, c, ub, xb, true)
            end
        end
    end
    return isfinite(ub) ? (lb, xb) : (Inf, x0)
end

# Pushes the simplex `sx` with its bound f(m) − lip·ρ, unless it lies beyond `s_stop`.
# Returns the smallest value found and its point, updated with f(m) on an admitted
# simplex.
function _simplex_bound!(heap, f, inside, lip, sx, ub, xb, ok::Bool)
    any(inside, sx) || return ub, xb
    m = sum(sx) / length(sx)
    lρ = lip * maximum(x -> norm(x - m), sx)
    fm = f(m)
    ok && inside(m) && fm < ub && ((ub, xb) = (fm, m))
    _heap_push!(heap, (fm - lρ, lρ, sx, ok))
    return ub, xb
end

function _heap_push!(h, x)
    push!(h, x)
    i = length(h)
    while i > 1 && h[i>>1][1] > h[i][1]
        h[i>>1], h[i] = h[i], h[i>>1]
        i >>= 1
    end
    return h
end

function _heap_pop!(h)
    top, last = h[1], pop!(h)
    isempty(h) && return top
    h[1] = last
    i, n = 1, length(h)
    while 2i <= n
        c = 2i < n && h[2i+1][1] < h[2i][1] ? 2i + 1 : 2i
        h[i][1] <= h[c][1] && break
        h[i], h[c] = h[c], h[i]
        i = c
    end
    return top
end

_simplex_split(s::NTuple{2,Vec3}) = (m = (s[1] + s[2]) / 2; ((s[1], m), (m, s[2])))
function _simplex_split(s::NTuple{3,Vec3})
    a, b, c = s
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    return ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))
end

"""
    beam_component_power(op, beams) -> Matrix

Returns the absorbed power (TW) per beam (rows) and component (columns), on both sides of
each component, from the deposition operator. Beams may differ from the operator's in
pointing, focus, spot and power; each is matched to the operator by its port.
"""
function beam_component_power(op::DepositionOperator, beams::AbstractVector{<:Beam})
    t = op.target
    nc = length(t.components)
    M = zeros(length(beams), nc)
    # Match the ports before the threaded loop, so that a beam without a port throws the
    # ArgumentError itself and not the CompositeException of a task.
    idx = map(beams) do b
        bi = findfirst(==(b.port), op.ports)
        bi === nothing && throw(ArgumentError("beam $(b.name) port is not in the operator"))
        bi
    end
    Threads.@threads for i in eachindex(beams)
        b, bi = beams[i], idx[i]
        E = zeros(nelements(t))
        L = zeros(nc)
        _deposit_beam!(E, L, op, bi, BeamFrame(_value_beam(b)))
        for e in eachindex(E)
            M[i, t.elem_comp[e]] += E[e]
        end
        M[i, :] .+= L
    end
    return M
end

# Returns the Float64 value of a number, through nested duals.
_fvalue(x::ForwardDiff.Dual) = _fvalue(ForwardDiff.value(x))
_fvalue(x::Real) = x

_value_spot(s::SuperGaussianSpot) = with_axes(s, _fvalue(s.a), _fvalue(s.b))

# Returns a Float64 copy of a beam, with ForwardDiff duals stripped, nested ones included.
_value_beam(b::Beam) =
    remake(b; aim = _fvalue.(b.aim), power = _fvalue(b.power), defocus = _fvalue(b.defocus),
           spot = _value_spot(b.spot), dlambda = _fvalue(b.dlambda))

"""
    violations(report; status = (:fail, :review), clearance = -Inf) -> Vector{ClearanceRow}

Returns the rows whose status is in `status`, or whose worst-case clearance is below
`clearance` mm, with failures first and then sorted by clearance.
"""
function violations(r::ClearanceReport; status = (:fail, :review), clearance = -Inf)
    rows = filter(x -> x.status in status || x.clearance < clearance, r.rows)
    rank(x) = x.status === :fail ? 0 : x.status === :review ? 1 : 2
    return sort(rows; by = x -> (rank(x), x.clearance))
end

"""
    min_clearance(report) -> Dict{Symbol,ClearanceRow}

Returns the closest beam, by worst-case clearance, for each hardware component.
"""
function min_clearance(r::ClearanceReport)
    out = Dict{Symbol,ClearanceRow}()
    for row in r.rows
        (!haskey(out, row.component) || row.clearance < out[row.component].clearance) &&
            (out[row.component] = row)
    end
    return out
end

function Base.show(io::IO, ::MIME"text/plain", r::ClearanceReport)
    nb = length(unique(x.beam for x in r.rows))
    e = r.errors
    @printf(io, "ClearanceReport: %d beams%s; envelope at %s + f-cone + %.3g mm margin + %.3g mm SSD\n",
            nb, r.canceled ? " (canceled: partial)" : "", r.contour, r.margin, e.ssd)
    @printf(io, "  error set: transverse %.3g mm (pointing %.3g + target %.3g + hardware %.3g), focus ±%.3g mm\n",
            r.transverse, e.pointing, e.target_offset, e.hardware_tolerance, e.focus)
    th = r.thresholds
    @printf(io, "  status: fail > %.0e W/cm² or %.0e of beam; review > %.0e W/cm² or %.0e\n",
            th.fail_intensity, th.fail_fraction, th.review_intensity, th.review_fraction)
    println(io, "  component        status  closest beam  worst (nominal) clearance mm  max power   max intensity")
    mc = min_clearance(r)
    for (c, row) in sort(collect(mc); by = x -> x[2].clearance)
        rs = [x for x in r.rows if x.component === c]
        pmax = maximum(x.power for x in rs)
        imax = maximum(x.peak_intensity for x in rs)
        st = any(x -> x.status === :fail, rs) ? :fail : any(x -> x.status === :review, rs) ? :review : :pass
        @printf(io, "  %-16s %-7s %-13s %+9.3f (%+9.3f)       %9.3g GW  %9.3g W/cm²\n", c, st,
                row.beam, row.clearance, row.nominal, 1e3pmax, imax)
    end
    nf = count(x -> x.status === :fail, r.rows)
    nr = count(x -> x.status === :review, r.rows)
    println(io, "  ", nf, " failing and ", nr, " review beam/component pairs")
end

"""
    write_clearance_csv(path, report)

Writes every row of a clearance report as CSV.
"""
function write_clearance_csv(path::AbstractString, r::ClearanceReport)
    open(path, "w") do io
        println(io, "beam,component,status,worst_clearance_mm,nominal_clearance_mm,x_mm,y_mm,z_mm,s_mm,",
                "worst_power_TW,worst_intercepted_fraction,worst_peak_intensity_Wcm2")
        for x in r.rows
            println(io, join((x.beam, x.component, x.status,
                              (@sprintf("%.6g", v) for v in (x.clearance, x.nominal, x.point..., x.s,
                                                              x.power, x.intercepted, x.peak_intensity))...), ","))
        end
    end
    return path
end
