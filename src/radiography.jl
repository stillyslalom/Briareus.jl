# Backlit x-ray radiography: a backlighter foil's He-α line through an object, a pinhole
# or point projection, filters and a detector, recorded time-integrated or in gated
# frames. The source is `backlighter_signal`'s chain (deposition, CE, angular
# distribution); the point-spread function and the gates are the x-ray imager's. Design:
# reference/diagnostics-design.md (phase 5); public values: research notes, "Backlit
# radiography: the model and its public definitions".

"""
    Radiograph(theta, phi; line, distance, pinhole = nothing, pinhole_distance = nothing,
               width, pixels = 128, center = 0, up = ẑ, filters = (), tables = (),
               detector = nothing, detector_pixel = 0.0, foil = :backlighter,
               angular = :isotropic, prepulse = false, ce = nothing, irradiance = :mean,
               pulse, start = 0.0, mode = :integrated, times = Float64[], gate = 0.0,
               gate_shape = :gaussian)

A backlit radiograph on the line of sight (θ, φ) (degrees, from the object toward the
detector): the He-α `line` (an element `Symbol` of [`kshell_line_energy`](@ref), or keV)
of the backlighter `foil` (a part of the target, behind the object), recorded on a
detector `distance` mm from `center` (a point, or a part's name). The object plane, a
`width` × `width` mm field about `center` normal to the line of sight sampled on
`pixels` × `pixels` points, is the image's coordinate system; `up` (a chamber
direction) is the image's up.

The optic:
- **Pinhole imaging** (`pinhole`, diameter mm, at `pinhole_distance` p mm from
  `center`): each image point sees the foil through the pinhole and the object point,
  at magnification `M = (distance − p)/p`. The point-spread function is the x-ray
  imager's ([`XrayImager`](@ref)): the penumbra `d(1 + 1/M)` and a Gaussian of the
  diffraction at the line energy and the detector pixel, on the object plane.
- **Point projection** (`pinhole = nothing`): each detector point receives every foil
  element as a point source at its centroid, along its own path through the object, so
  the source's extent blurs the image (by the source width × `L/(s + L)` on the object
  plane, `s` the foil's distance behind the object and `L = distance`). The image point
  of a detector point is the object-plane crossing of its line to the foil's centroid,
  and `M = (s + L)/s`. Only the detector pixel blurs further.

The spectral response at the line is the product of the filters' transmission
(`filter_transmission` materials with thicknesses in mm, NIST's tables where they cover
the material and energy, as [`backlighter_signal`](@ref)), measured `tables`
([`TabulatedCurve`](@ref)s, eV) and the detector's efficiency: `detector = (:Si, 0.015)`
absorbs `1 − e^{−μρt}` (Henke photoabsorption, no dead layer), a `TabulatedCurve` is a
measured efficiency, `nothing` absorbs every photon (`rg.filter_transmission` and
`rg.efficiency` hold the two at the line). Images are of detected photons.

The source is that of [`backlighter_signal`](@ref) (`prepulse`, `ce`, `irradiance`,
`angular` as there): the He-α energy of each foil element, emitted uniformly over the
flat-top laser pulse of `pulse` ns from `start` ns (an assumption: the emission follows
the laser). `mode` is `:integrated` (the whole emission) or `:gated` (one frame per
time of `times`, ns, each recording the emission against a gate of FWHM `gate` ns,
`:gaussian` or `:square`, of peak 1, as [`imager_record`](@ref)). The object is static
over the emission: for a moving object, one radiograph per frame of its geometry.
"""
struct Radiograph
    theta::Float64
    phi::Float64
    camera::OrthoCamera
    center::Union{Vec3,Symbol}
    line::Float64
    distance::Float64
    pinhole::Union{Nothing,Float64}
    pinhole_distance::Float64
    detector_pixel::Float64
    filter_transmission::Float64
    efficiency::Float64
    foil::Symbol
    angular::Any
    prepulse::Bool
    ce::Union{Nothing,Float64}
    irradiance::Symbol
    pulse::Float64
    start::Float64
    mode::Symbol
    times::Vector{Float64}
    gate::Float64
    gate_shape::Symbol
end

function Radiograph(theta::Real, phi::Real; line, distance::Real, pinhole = nothing, pinhole_distance = nothing,
                    width::Real, pixels::Integer = 128, center = zero(Vec3), up = Vec3(0, 0, 1),
                    filters = (), tables = (), detector = nothing, detector_pixel::Real = 0.0,
                    foil::Symbol = :backlighter, angular = :isotropic, prepulse::Bool = false, ce = nothing,
                    irradiance::Symbol = :mean, pulse::Real, start::Real = 0.0, mode::Symbol = :integrated,
                    times = Float64[], gate::Real = 0.0, gate_shape::Symbol = :gaussian)
    distance > 0 && width > 0 && pixels > 0 && pulse > 0 ||
        throw(ArgumentError("distance, width, pixels and pulse must be positive"))
    detector_pixel >= 0 || throw(ArgumentError("detector_pixel must be nonnegative"))
    if pinhole === nothing
        pinhole_distance === nothing || throw(ArgumentError("pinhole_distance needs a pinhole"))
        p = NaN
    else
        pinhole > 0 || throw(ArgumentError("the pinhole diameter must be positive"))
        pinhole_distance === nothing && throw(ArgumentError("give pinhole_distance with pinhole"))
        p = Float64(pinhole_distance)
        0 < p < distance || throw(ArgumentError("need 0 < pinhole_distance < distance"))
    end
    ce === nothing || 0 <= ce <= 1 || throw(ArgumentError("ce must be a fraction in [0, 1]"))
    irradiance in (:mean, :element) || throw(ArgumentError("irradiance must be :mean or :element"))
    mode in (:integrated, :gated) || throw(ArgumentError("mode is :integrated or :gated"))
    gate_shape in (:gaussian, :square) || throw(ArgumentError("gate_shape is :gaussian or :square"))
    mode === :gated && (isempty(times) || gate <= 0) &&
        throw(ArgumentError("a gated radiograph needs times and a positive gate"))
    hν = _line_keV(line)
    f = Pair{Symbol,Union{Float64,Tuple{Float64,Float64}}}[k => _filter_spec(v) for (k, v) in filters]
    T = filter_transmission(f, hν)
    for c in tables
        T *= c(1000hν)
    end
    cam = OrthoCamera(-port_vector(theta, phi); up, width, nx = pixels)
    return Radiograph(Float64(theta), Float64(phi), cam, center isa Symbol ? center : Vec3(center...), hν,
                      Float64(distance), pinhole === nothing ? nothing : Float64(pinhole), p,
                      Float64(detector_pixel), T, _detection(detector, 1000hν), foil, angular, prepulse,
                      ce === nothing ? nothing : Float64(ce), irradiance, Float64(pulse), Float64(start), mode,
                      collect(Float64, times), Float64(gate), gate_shape)
end

# ---- the object's attenuation -------------------------------------------------------

# Every crossing (s, ±1) of the segment o + s r, 0 < s < 1, with the triangles of `b`,
# the sign that of r · (e1 × e2), into `buf`
function _crossings!(buf::Vector{Tuple{Float64,Int}}, b::BVH, o::Vec3, r::Vec3)
    empty!(buf)
    isempty(b.nodes) && return buf
    inv = _safe_inv(r)
    stack = MVector{64,Int32}(undef)
    sp = 1
    stack[1] = 1
    @inbounds while sp > 0
        ni = stack[sp]; sp -= 1
        node = b.nodes[ni]
        _slab(node, o, inv, 1.0) || continue
        if node.b < 0
            for k in node.a:(node.a - node.b - 1)
                it = b.items[k]
                s = _intersect(o, r, b.v0[it], b.e1[it], b.e2[it], 1.0)
                s < Inf && push!(buf, (s, dot(r, cross(b.e1[it], b.e2[it])) > 0 ? 1 : -1))
            end
        else
            stack[sp+1] = node.a
            stack[sp+2] = node.b
            sp += 2
        end
    end
    return buf
end

# The length (mm) of the segment from `o` to `o + r` inside the closed surface `b`: the
# signed sum of its crossings, Σ ±s |r|, whose magnitude is the inside length for a
# consistently oriented surface (either way) and ends outside it. A crossing on an edge
# or vertex shared by several triangles counts once (crossings of one sign within 1e-9
# of each other in s).
function _inside_length(b::BVH, o::Vec3, r::Vec3, buf)
    _crossings!(buf, b, o, r)
    sort!(buf)
    S, last, sgn = 0.0, -Inf, 0
    for (s, σ) in buf
        (s - last <= 1e-9 && σ == sgn) && continue
        S += σ * s
        last, sgn = s, σ
    end
    return abs(S) * norm(r)
end

# The object `parts => material` as closed-surface BVHs and attenuation coefficients
# (1/mm) at `hν` keV: `part => :Al` (the material's density) or `part => (:Al, ρ)`
function _object(t::Target, object, hν)
    bvhs, κ, parts = BVH[], Float64[], Symbol[]
    for (part, spec) in object
        mat, ρ = spec isa Tuple ? (Symbol(spec[1]), Float64(spec[2])) : (Symbol(spec), _density(Symbol(spec)))
        isnan(ρ) && throw(ArgumentError("no density for :$mat; pass $part => (:$mat, density)"))
        cis = Set(part_pieces(t, Symbol(part)))
        ks = [k for k in eachindex(t.tri_comp) if t.tri_comp[k] in cis]
        isempty(ks) && throw(ArgumentError("the object part :$part has no triangles"))
        push!(bvhs, BVH(t.tri_v0[ks], t.tri_e1[ks], t.tri_e2[ks], Int32.(ks)))
        push!(κ, mass_attenuation(mat, hν) * ρ * 0.1)
        push!(parts, Symbol(part))
    end
    return bvhs, κ, parts
end

# transmission of the segment o → o + r through the object
function _object_transmission(bvhs, κ, o::Vec3, r::Vec3, buf)
    τ = 0.0
    for (b, k) in zip(bvhs, κ)
        τ += k * _inside_length(b, o, r, buf)
    end
    return exp(-τ)
end

# ---- the radiograph -----------------------------------------------------------------

"""
    radiograph(rg::Radiograph, target, deposition; object = (),
               model = kshell_ce_model()) -> NamedTuple

The radiograph `rg` records of `target`, its backlighter foil driven by `deposition`
(the laser power on each element, TW, as [`backlighter_signal`](@ref) takes it).

`object` lists the parts that attenuate, `part => material` (the material's density) or
`part => (material, density)` (g/cm³), materials as [`mass_attenuation`](@ref)'s: each
part is a closed surface filled with its material, and along each ray the transmission
is `exp(−Σ μ/ρ · ρ · L)`, `L` the ray's length inside each part. Where parts overlap their
attenuations add, so a hollow shell is its outer surface with the shell's material and
its inner surface with the negative density. The object's parts do not obstruct;
every other surface part except the foil does (as hardware in
[`backlighter_signal`](@ref)).

Returns `(; images, open, transmission, weights, times, magnification, line, ce,
sigma, extrapolated)`:
- `images`: one matrix per frame (one for `:integrated`), detected photons per mm² of
  detector at each object-plane point (row 1 at the top, `up`); per detector pixel,
  times its area.
- `open`: the same over the whole emission without the object (the open field).
- `transmission`: the whole emission's image over `open` (NaN where `open` is 0): the
  object's transmission as the radiograph measures it, blurred.
- `weights`: the fraction of the emission each frame records (1 for `:integrated`), so
  `images[k] = weights[k] × transmission × open`; `times`: the frames' times.
- `magnification`, `line` (keV), and from the source `ce` (the He-α energy over the
  laser energy), `sigma` (1σ of `ln CE`, by which every image is uncertain) and
  `extrapolated`.
"""
function radiograph(rg::Radiograph, t::Target, dep::Deposition; object = (),
                    model::KShellCEModel = kshell_ce_model())
    src = _backlighter_source(t, dep; line = rg.line, pulse = rg.pulse, foil = rg.foil, prepulse = rg.prepulse,
                              ce = rg.ce, irradiance = rg.irradiance, model)
    bvhs, κ, parts = _object(t, object, rg.line)
    _, skip = _hardware(t, rg.foil, parts)
    k = rg.camera
    c = rg.center isa Symbol ? part_centroid(t, rg.center) : rg.center
    d = -k.dir                                                  # toward the detector
    fe = src.fe
    fc = sum(t.elem_area[e] * t.elem_centroid[e] for e in fe) / sum(t.elem_area[e] for e in fe)
    gain = rg.filter_transmission * rg.efficiency / (rg.line * _KEV_J)     # photons per J at the detector
    h = k.width / k.nx
    if rg.pinhole === nothing
        s = dot(c - fc, d)
        s > 0 || throw(ArgumentError("the foil :$(rg.foil) must lie behind the object, away from the detector"))
        M = (s + rg.distance) / s
        psf = _psf(h, 0.0, 0.0, rg.detector_pixel / M)
    else
        p = rg.pinhole_distance
        M = (rg.distance - p) / p
        λ = 1.239841984e-6 / rg.line                            # mm
        psf = _psf(h, rg.pinhole * (1 + 1 / M), 1.03λ * p / rg.pinhole, rg.detector_pixel / M)
    end
    raw, opn = _radiograph_pixels(rg, t, src, bvhs, κ, skip, c, d, fc, rg.angular)
    img = _blur(psf, gain .* raw)
    op = _blur(psf, gain .* opn)
    w = rg.mode === :integrated ? [1.0] : _emission_weights(rg)
    return (; images = [wk .* img for wk in w], open = op,
            transmission = [o > 0 ? i / o : NaN for (i, o) in zip(img, op)], weights = w,
            times = rg.mode === :integrated ? Float64[] : copy(rg.times), magnification = M, line = rg.line,
            ce = src.Ex / src.EL, sigma = src.σ, extrapolated = src.ext)
end

# the fluence (J/mm² of He-α at the detector, before the filters) with and without the
# object at each object-plane point (`angular` passed through a function barrier)
function _radiograph_pixels(rg::Radiograph, t::Target, src, bvhs, κ, skip, c::Vec3, d::Vec3, fc::Vec3,
                            angular::A) where {A}
    k = rg.camera
    raw = zeros(k.ny, k.nx)
    opn = zeros(k.ny, k.nx)
    ε = _ray_eps(t)
    fe, Ee = src.fe, src.Ee
    fi = Set(part_pieces(t, rg.foil))
    ftris = [q for q in eachindex(t.tri_comp) if t.tri_comp[q] in fi]
    idx = Dict(e => i for (i, e) in enumerate(fe))
    Threads.@threads for j in 1:k.nx
        buf = Tuple{Float64,Int}[]
        x = (j - 0.5) / k.nx * k.width - k.width / 2
        for i in 1:k.ny
            y = k.height / 2 - (i - 0.5) / k.ny * k.height
            q = c + x * k.right + y * k.up
            if rg.pinhole === nothing
                # the detector point on the line from the foil's centroid through q
                Pd = fc + (q - fc) * (dot(c + rg.distance * d - fc, d) / dot(q - fc, d))
                tr = (o, r) -> _object_transmission(bvhs, κ, o, r, buf)
                opn[i, j] = _point_fluence(t, fe, Ee, angular, Pd, skip, ε)
                raw[i, j] = _point_fluence(t, fe, Ee, angular, Pd, skip, ε, tr)
            else
                Pp = c + rg.pinhole_distance * d
                R, X = _foil_radiance(t, ftris, idx, Ee, angular, Pp, q - Pp, skip, ε)
                R > 0 || continue
                f = R * π * (rg.pinhole / 2)^2 / (rg.distance - rg.pinhole_distance)^2
                opn[i, j] = f
                raw[i, j] = f * _object_transmission(bvhs, κ, X, Pp - X, buf)
            end
        end
    end
    return raw, opn
end

# the fraction of the emission (uniform over the pulse) each gate records: the gates of
# `imager_record` integrated over the pulse on 400 intervals of 64 points
function _emission_weights(rg::Radiograph)
    tg = collect(range(rg.start, rg.start + rg.pulse; length = 401))
    return [sum(_gate_weights(tg, _gate(rg.gate_shape, t0, rg.gate))) / rg.pulse for t0 in rg.times]
end
