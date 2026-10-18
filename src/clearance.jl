# Beam clearance: beam envelopes at an intensity contour, worst-case hardware clearance
# over pointing and focus tolerances, and worst-case power and intensity on hardware.

"""
    enclosed_radius(order, fraction)

Radius (in units of the 1/e semi-axes) of the super-Gaussian contour enclosing
`fraction` of the power: the fraction inside `ρ` is `P(2/n, ρⁿ)` (regularized lower
incomplete gamma function). For spots with a halo use `envelope_radius`.
"""
enclosed_radius(n::Real, f::Real) = gamma_inc_inv(2 / n, f, 1 - f)^(1 / n)

"""
    envelope_radius(spot, power; intensity = nothing, peak_fraction = nothing,
                    enclosed = nothing) -> ρ

Best-focus contour radius in units of the core 1/e semi-axes (halo included), set by
exactly one of: an absolute `intensity` (W/cm²), a fraction of the peak focal intensity
`peak_fraction`, or the `enclosed` power fraction. Returns 0 if the whole spot is below
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

Geometric beam envelope: at axial distance `s` from best focus, the best-focus contour
ellipse (see `envelope_radius`) grown by the f-cone radius `|s|/(2F)` and a uniform
`margin` (mm). In geometric optics with a phase plate (every lens point fills the
whole spot), the intensity at a point is the mean of the focal intensity over a disk of
radius `|s|/(2F)`, so outside the envelope it is below the contour level at every `s`.
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

Signed distance (mm) from point `x` to the beam envelope in the plane through `x`
normal to the beam, and the axial position `s` of that plane relative to best focus.
"""
function envelope_distance(env::BeamEnvelope, x::Vec3)
    d = x - env.focus
    s = dot(d, env.k)
    u, v = dot(d, env.e1), dot(d, env.e2)
    return _ellipse_distance(u, v, env.a, env.b) - env.cone * abs(s) - env.margin, s
end

# Signed distance from (u, v) to the ellipse with semi-axes (a, b): negative inside.
# Eberly's robust point–ellipse distance (bisection on the Lagrange multiplier), valid
# inside and outside. A degenerate ellipse (a = b = 0) is the point at the origin.
function _ellipse_distance(u, v, a, b)
    (a == 0 || b == 0) && return hypot(u, v)
    y0, y1 = abs(u), abs(v)
    e0, e1 = a, b
    if e0 < e1
        e0, e1, y0, y1 = e1, e0, y1, y0
    end
    inside = (y0 / e0)^2 + (y1 / e1)^2 < 1
    if y1 > 0
        if y0 > 0
            z0, z1 = y0 / e0, y1 / e1
            g = z0^2 + z1^2 - 1
            g == 0 && return 0.0
            r0 = (e0 / e1)^2
            n0 = r0 * z0
            s0, s1 = z1 - 1, g < 0 ? 0.0 : hypot(n0, z1) - 1
            s = s0
            for _ in 1:200
                s = (s0 + s1) / 2
                (s == s0 || s == s1) && break
                gs = (n0 / (s + r0))^2 + (z1 / (s + 1))^2 - 1
                gs > 0 ? (s0 = s) : gs < 0 ? (s1 = s) : break
            end
            x0, x1 = r0 * y0 / (s + r0), y1 / (s + 1)
            d = hypot(x0 - y0, x1 - y1)
        else
            d = abs(y1 - e1)
        end
    else
        numer, denom = e0 * y0, e0^2 - e1^2
        if numer < denom
            ξ = numer / denom
            d = hypot(e0 * ξ - y0, e1 * sqrt(1 - ξ^2))
        else
            d = abs(y0 - e0)
        end
    end
    return inside ? -d : d
end

"""
    SpotWing

Supertype of low-level circular far-field components outside a beam's phase-plate
core, used only by clearance checks (they carry ≲1% of the energy and do not matter for
drive symmetry). Each wing is centered on the beam axis at best focus; radii in mm.
"""
abstract type SpotWing end

"""
    PowerLawWing(fraction; rh = :r95, rmax = 1.5)

Halo with `I ∝ [1 + (r/rh)²]^(-3/2)` for `r ≤ rmax`, carrying `fraction` of the beam
power. `rh = :r95` uses the radius enclosing 95% of the core (geometric mean of the
axes). `rmax` is the spatial-filter cutoff (NIF ≈ ±200 µrad × 7.7 m ≈ 1.5 mm).
"""
struct PowerLawWing <: SpotWing
    fraction::Float64
    rh::Union{Float64,Symbol}
    rmax::Float64
end
PowerLawWing(fraction; rh = :r95, rmax = 1.5) = PowerLawWing(fraction, rh isa Symbol ? rh : Float64(rh), rmax)

"""
    UniformScatter(fraction; radius = 23.0)

Far-angle scatter floor: `fraction` of the beam power spread uniformly over a disk of
`radius` mm at the target plane (NIF final-optics scatter up to ≈3 mrad × 7.7 m).
"""
struct UniformScatter <: SpotWing
    fraction::Float64
    radius::Float64
end
UniformScatter(fraction; radius = 23.0) = UniformScatter(fraction, radius)

"""
    nif_wings(; halo = 0.01, rmax = 1.5, scatter = 1e-3, scatter_radius = 23.0)

Default NIF clearance wing model. These are engineering defaults (see
reference/research-notes.md), not facility criteria; override them with measured
far-field data:
- a 1% power-law halo (OMEGA measured < 1% outside the super-Gaussian; the NIF/LMJ
  design specs allowed 2% outside the nominal spot), truncated at the spatial-filter
  angle;
- a 10⁻³ uniform scatter floor out to 3 mrad.
"""
nif_wings(; halo = 0.01, rmax = 1.5, scatter = 1e-3, scatter_radius = 23.0) =
    SpotWing[PowerLawWing(halo; rmax), UniformScatter(scatter; radius = scatter_radius)]

_wing_rh(w::PowerLawWing, b::Beam) =
    w.rh === :r95 ? envelope_radius(b.spot, 1.0; enclosed = 0.95) * sqrt(ForwardDiff.value(b.spot.a) * ForwardDiff.value(b.spot.b)) :
    w.rh isa Symbol ? throw(ArgumentError("unknown rh $(w.rh)")) : w.rh

# intensity (TW/mm²) of a wing at transverse radius r for beam power P
function _wing_intensity(w::PowerLawWing, rh, P, r)
    r > w.rmax && return 0.0
    norm = 2π * rh^2 * (1 - 1 / sqrt(1 + (w.rmax / rh)^2))
    return w.fraction * P / norm * (1 + (r / rh)^2)^(-1.5)
end
_wing_intensity(w::UniformScatter, _, P, r) = r > w.radius ? 0.0 : w.fraction * P / (π * w.radius^2)

# radius (mm) of a wing's contour at intensity thr (TW/mm²); 0 if never that bright
function _wing_contour(w::PowerLawWing, rh, P, thr)
    I0 = _wing_intensity(w, rh, P, 0.0)
    I0 <= thr && return 0.0
    return min(w.rmax, rh * sqrt((I0 / thr)^(2 / 3) - 1))
end
_wing_contour(w::UniformScatter, _, P, thr) = _wing_intensity(w, nothing, P, 0.0) > thr ? w.radius : 0.0

# cumulative power fraction of a wing inside radius r
_wing_cdf(w::PowerLawWing, rh, r) = (1 - 1 / sqrt(1 + (min(r, w.rmax) / rh)^2)) / (1 - 1 / sqrt(1 + (w.rmax / rh)^2))
_wing_cdf(w::UniformScatter, _, r) = min(r, w.radius)^2 / w.radius^2
_wing_rmax(w::PowerLawWing) = w.rmax
_wing_rmax(w::UniformScatter) = w.radius

# Ray trace a circular wing (polar grid, equal-power annuli) into per-component power.
function _raytrace_wing!(L, t::Target, b::Beam, w::SpotWing, rh; nr = 48, nθ = 64, ncone = 7)
    val = ForwardDiff.value
    e1, e2 = spot_axes(b)
    k = direction(b)
    focus = Vec3(val.(b.aim)) + val(b.defocus) * k
    dirs, dw = cone_directions(b, ncone)
    far = 4 * (t.scale + norm(focus)) + 1
    P = val(b.power) * w.fraction
    rmax = _wing_rmax(w)
    # annulus edges at equal power, sample at each annulus' median radius
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

One beam/hardware-component pair. "Worst" values are over the error set (pointing,
target offset, hardware tolerance, focus):

- `clearance` — worst-case signed distance to the beam envelope (mm; negative = overlap)
- `nominal` — the same at nominal pointing and focus
- `point`, `s` — hardware point of minimum clearance and its axial position from best
  focus (mm)
- `power`, `intercepted` — worst-case power (TW) on the component (core + wings) and
  its fraction of the beam power
- `peak_intensity` — worst-case beam intensity (W/cm², normal to the beam, core +
  wings) at the component's upstream surface points, ignoring shadowing
- `status` — `:fail`, `:review` or `:pass` against the report's thresholds
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

Result of `clearance`: one `ClearanceRow` per beam and hardware component plus the
envelope, error-set and threshold settings. See `violations`, `min_clearance`.
`canceled` is true when a `clearance` callback stopped the check early; `rows` then
hold only the beams completed before the cancel (the leading beams, in input order).
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

Status tiers of `clearance` (engineering defaults, not facility criteria): `:fail`
above 1e12 W/cm² (plasma formation) or 1e-3 of the beam power on a component,
`:review` above 1e9 W/cm² (ablation and debris onset) or 1e-5.
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
              op = nothing, nspot = 61, callback = nothing,
              batch = Threads.nthreads()) -> ClearanceReport

Beam-to-hardware clearance check, judged on the worst case over an error set.

**Spot.** Each beam's `SuperGaussianSpot` core (with its optional super-Gaussian halo)
plus clearance-only circular `wings` (default `nif_wings()`: a 1% power-law halo and a
10⁻³ scatter floor — engineering defaults, override with measured far fields). `ssd`
(mm) is the half-width of the SSD smear along the dispersion axis, added isotropically
to the envelope (conservative).

**Envelope.** The best-focus contour at an absolute `intensity` (W/cm²), a
`peak_fraction` of the core peak, or an `enclosed` core fraction (give one), grown by
the f-cone, `ssd` and `margin`. With wings, each spot component's contour is taken at
the threshold divided by the number of components, and the envelope is their union,
so the total intensity outside it is below the threshold.

**Error set.** Transverse errors add linearly (not in quadrature):
`R = pointing + target_offset + hardware_tolerance` (mm) — per-beam pointing (NIF
50 µm rms spec; pass ≈ 0.12 for a 3σ statistical radius), rigid target offset (NIF
mean beams-to-target ≈ 25 µm), and as-built hardware tolerance (metrology, or
0.1–0.25 mm without it). `focus` is the ± focal-position error (mm; 1994 NIF spec ≤ 1 mm).
Worst-case clearance is exact: nominal − `R` − `focus/(2F)`. Worst-case power and
intensity are maxima over aims displaced to `nangles` points on the radius-`R` circle
(plus nominal) at focus offsets `0, ±focus`.

**Status.** `thresholds` (W/cm² and fraction of beam power): `:fail` above
`fail_intensity` (1e12, plasma formation) or `fail_fraction` (1e-3), `:review` above
`review_intensity` (1e9, ablation/debris onset) or `review_fraction` (1e-5). These are
engineering defaults; NIF (TaLIS/BLIP) and LLE review criteria are not public.

Only hardware upstream of each beam's landing is checked: where the beam axis first
meets one of the `targets` components (the surfaces the beam should hit), or the aim
plane if it misses them. Hardware is sampled at triangle vertices and `4^level`
sub-triangle points. Power on hardware comes from the reference ray trace (`nspot`
grid per core component; wings on equal-power polar grids), or for the core from `op`,
a `DepositionOperator` for these ports, when given.

**Progress and cancel.** With a `callback`, beams are checked in batches of `batch`
(input order) and `callback((; done, total))` is called with the number of beams
completed — once with `done = 0` before any work, then after each batch. Returning
`true` cancels: `clearance` returns at once with the rows of the completed beams and
`canceled = true` in the report (a `true` after the last batch is ignored). The
callback always runs on the calling task, never on a worker thread, so it may touch
task-local state; it must not block for long, since no work runs while it does. Rows
do not depend on `callback` or `batch`.
"""
function clearance(t::Target, beams::AbstractVector{<:Beam}; targets,
                   hardware = nothing, intensity = nothing, peak_fraction = nothing,
                   enclosed = nothing, wings = nif_wings(), pointing = 0.05,
                   target_offset = 0.025, hardware_tolerance = 0.0, focus = 1.0,
                   margin = 0.0, ssd = 0.0, thresholds = DEFAULT_CLEARANCE_THRESHOLDS,
                   nangles::Integer = 16, level::Integer = 2,
                   op::Union{Nothing,DepositionOperator} = nothing, nspot::Integer = 61,
                   callback = nothing, batch::Integer = Threads.nthreads())
    all(isnothing, (intensity, peak_fraction, enclosed)) && (peak_fraction = 1e-4)
    contour = intensity !== nothing ? @sprintf("%.3g W/cm²", intensity) :
              peak_fraction !== nothing ? @sprintf("%.3g of core peak", peak_fraction) :
              @sprintf("%.6g of core enclosed", enclosed)
    R = pointing + target_offset + hardware_tolerance
    tgt = Symbol[targets...]
    hw = hardware === nothing ?
         [c.name for c in t.components if c.kind === :surface && !(c.name in tgt)] :
         Symbol[hardware...]
    tgt_idx = Set(component_index(t, n) for n in tgt)
    hw_idx = [component_index(t, n) for n in hw]
    samples = Dict{Int,Vector{Vec3}}()
    for ci in hw_idx
        pts = Vec3[]
        for k in eachindex(t.tri_area)
            t.tri_comp[k] == ci || continue
            v0, e1, e2 = t.tri_v0[k], t.tri_e1[k], t.tri_e2[k]
            _subtri_vertices!(pts, v0, v0 + e1, v0 + e2, level)
        end
        samples[ci] = unique(pts)
    end
    errs = (; pointing = Float64(pointing), target_offset = Float64(target_offset),
            hardware_tolerance = Float64(hardware_tolerance), focus = Float64(focus),
            ssd = Float64(ssd))
    finish(rows, canceled) = ClearanceReport(rows, contour, R, errs, Float64(margin), tgt,
                                             thresholds, canceled)
    rowsof(bs) = _clearance_rows(t, bs, samples, hw_idx, tgt_idx, wings, op, nspot, R,
                                 focus, nangles, thresholds;
                                 intensity, peak_fraction, enclosed, margin = margin + ssd)
    callback === nothing && return finish(rowsof(beams), false)

    # In batches of beams, reporting progress (and polling for a cancel) between them.
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

# Rows of the clearance check for `beams` (all beams of a non-batched call).
function _clearance_rows(t::Target, beams, samples, hw_idx, tgt_idx, wings, op, nspot, R,
                         focus, nangles, thresholds; intensity, peak_fraction, enclosed,
                         margin)
    nc = length(t.components)
    # beam variants over the error set (nominal first for each beam)
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
        for (j, ci) in enumerate(hw_idx)
            best, bp, bs, Imax = _scan_hardware(samples[ci], envs[1], wingc, s_stop, frames,
                                                cone, wt, rhs, P)
            worst = best - R - focus * envs[1].cone
            pw = maximum(Pc[v, ci] for v in vids)
            frac = P > 0 ? pw / P : 0.0
            Iw = to_Wcm2(Imax)
            rows[(bi-1)*length(hw_idx)+j] = ClearanceRow(b.name, t.components[ci].name, worst, best,
                                                         bp, bs, pw, frac, Iw,
                                                         _status(Iw, frac, thresholds))
        end
    end
    return rows
end

const _CLEARANCE_NCONE = 37

# Nearest hardware sample to the envelope (distance, point, axial position) upstream of
# `s_stop`, and the largest intensity over the error-set variants `frames` there.
function _scan_hardware(samples, env::BeamEnvelope, wingc, s_stop, frames, cone, wings::Tuple,
                        rhs, P)
    best, bp, bs, Imax = Inf, zero(Vec3), NaN, 0.0
    # transverse slopes of the cone directions (all variants share the port)
    f1 = frames[1]
    αmax = maximum(dj -> abs(dot(dj, f1.e1) / dot(dj, f1.k)), cone)
    βmax = maximum(dj -> abs(dot(dj, f1.e2) / dot(dj, f1.k)), cone)
    for x in samples
        d, s = _envelope_union_distance((env,), wingc, x)
        s <= s_stop || continue
        d < best && ((best, bp, bs) = (d, x, s))
        for f in frames
            I = _beyond(f, αmax, βmax, x) ? 0.0 : _cone_mean_intensity(f, cone, x)
            if !isempty(wings)
                dv = x - f.focus
                r = hypot(dot(dv, f.e1), dot(dv, f.e2))
                I += _wing_sum(wings, rhs, P, r)
            end
            Imax = max(Imax, I)
        end
    end
    return best, bp, bs, Imax
end

_wing_sum(wings::Tuple, rhs::Tuple, P, r) =
    sum(map((w, rh) -> _wing_intensity(w, rh, P, r), wings, rhs); init = 0.0)

# Beam intensity at `x`: the mean over cone directions of the best-focus profile where
# each ray through `x` crosses the focal plane, i.e. the f-cone convolution that
# deposition uses. (The moment-matched `intensity` sharpens the edge of a defocused
# flat-topped spot, where clearance is judged.)
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

# Core envelope (contour at thr/N for N spot components) and wing contour radii.
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

# Signed distance to the union of the core envelope and the wing disks (each grown by
# the f-cone and the core envelope's margin).
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

# Axial position (relative to best focus) where the beam axis first meets a target
# component; the aim plane if it never does (or `miss`, when given).
function _landing(t::Target, env::BeamEnvelope, tgt_idx; miss = nothing)
    far = 4 * (t.scale + norm(env.aim)) + 1
    o = env.aim - far * env.k
    traveled = 0.0
    while true
        τ, tri = closest_hit(t.bvh, o, env.k)
        tri == 0 && break
        x = o + τ * env.k
        t.tri_comp[tri] in tgt_idx && return dot(x - env.focus, env.k)
        o = x + _ray_eps(t) * env.k
        traveled += τ
        traveled > 2far && break
    end
    return miss === nothing ? dot(env.aim - env.focus, env.k) : miss
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

"""
    beam_component_power(op, beams) -> Matrix

Absorbed power (TW) per beam (rows) and component (columns), both sides of each
component, from the deposition operator. Beams may differ from the operator's in
pointing, focus, spot and power; each is matched to the operator by its port.
"""
function beam_component_power(op::DepositionOperator, beams::AbstractVector{<:Beam})
    t = op.target
    nc = length(t.components)
    M = zeros(length(beams), nc)
    Threads.@threads for i in eachindex(beams)
        b = beams[i]
        bi = findfirst(==(b.port), op.ports)
        bi === nothing && throw(ArgumentError("beam $(b.name) port is not in the operator"))
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

# The Float64 value of a number, through nested duals.
_fvalue(x::ForwardDiff.Dual) = _fvalue(ForwardDiff.value(x))
_fvalue(x::Real) = x

_value_spot(s::SuperGaussianSpot) = with_axes(s, _fvalue(s.a), _fvalue(s.b))

# Float64 copy of a beam (strips ForwardDiff duals, nested ones included).
_value_beam(b::Beam) =
    remake(b; aim = _fvalue.(b.aim), power = _fvalue(b.power), defocus = _fvalue(b.defocus),
           spot = _value_spot(b.spot), dlambda = _fvalue(b.dlambda))

"""
    violations(report; status = (:fail, :review), clearance = -Inf) -> Vector{ClearanceRow}

Rows whose status is in `status`, or whose worst-case clearance is below `clearance`
mm; failures first, then by clearance.
"""
function violations(r::ClearanceReport; status = (:fail, :review), clearance = -Inf)
    rows = filter(x -> x.status in status || x.clearance < clearance, r.rows)
    rank(x) = x.status === :fail ? 0 : x.status === :review ? 1 : 2
    return sort(rows; by = x -> (rank(x), x.clearance))
end

"""
    min_clearance(report) -> Dict{Symbol,ClearanceRow}

Closest beam (worst-case clearance) for each hardware component.
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

Write every row of a clearance report as CSV.
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
