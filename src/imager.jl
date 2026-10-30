# Synthetic x-ray pinhole imagers: the target's x-ray emission seen through a pinhole,
# filters and a detector, recorded time-integrated or in gated frames. Design:
# reference/diagnostics-design.md (phase 5); public values: research notes, "X-ray
# imagers, radiography and streaked diagnostics: what is public".

"""
    XrayImager(theta, phi; pinhole, distance, magnification, width, pixels = 256,
               center = 0, up = ẑ, filters = (), tables = (), mirror = nothing,
               detector = nothing, detector_pixel = 0.0, energy = nothing,
               mode = :integrated, times = Float64[], gate = 0.0, gate_shape = :gaussian)

A pinhole x-ray imager on the line of sight (θ, φ) (degrees, from target chamber
center toward the instrument): a pinhole of diameter `pinhole` (mm) at `distance` (mm)
from the object, imaging at `magnification` onto a detector of pixel `detector_pixel`
(mm). The object plane is sampled over a `width` × `width` mm field about `center` (a
point, or a part's name) on `pixels` × `pixels` rays; `up` (a chamber direction) is the
image's up.

The spectral response `R(E)` is the product of the filters (`filter_transmission`
materials with thicknesses in mm, as [`DanteResponse`](@ref)'s, and measured `tables`),
an optional `mirror` (as [`DanteResponse`](@ref)'s) and the detector's absorption:
`detector = (:Si, 0.015)` absorbs `1 − e^{−μρt}` of the photons in 15 µm of Si (no
dead layer), a [`TabulatedCurve`](@ref) is a measured efficiency, and `nothing`
absorbs every photon. `R` is a fraction, so images are of detected x-ray energy.

The point-spread function on the object plane is the pinhole's geometric penumbra, a
uniform disk of diameter `penumbra = d(1 + 1/M)`, convolved with a Gaussian of the
summed variances of its Fraunhofer diffraction (FWHM `diffraction = 1.03 λ
distance / d` at photon `energy`, eV) and of the detector pixel (a box of
`detector_pixel / M`). The Gaussian is an engineering approximation: the Airy pattern
and the box are not Gaussian, and a pinhole whose Fresnel number `d² / λz` is near 1 is
not in the far field. `energy` defaults to the mean photon energy of what the imager
detects of a 300 eV Planckian (an engineering default); a definition states it.

`mode` is `:integrated` (one image of the time integral) or `:gated` (one frame per
time of `times`, ns, each the integral of the image against a gate of FWHM `gate` ns,
`:gaussian` or `:square`, of peak 1).
"""
struct XrayImager
    theta::Float64
    phi::Float64
    camera::OrthoCamera
    center::Union{Vec3,Symbol}
    pinhole::Float64
    distance::Float64
    magnification::Float64
    detector_pixel::Float64
    E::Vector{Float64}
    R::Vector{Float64}
    energy::Float64
    penumbra::Float64
    diffraction::Float64
    mode::Symbol
    times::Vector{Float64}
    gate::Float64
    gate_shape::Symbol
    cache::IdDict{Target,Any}
    lock::ReentrantLock
end

function XrayImager(theta::Real, phi::Real; pinhole::Real, distance::Real, magnification::Real,
                    width::Real, pixels::Integer = 256, center = zero(Vec3), up = Vec3(0, 0, 1),
                    filters = (), tables = (), mirror = nothing, detector = nothing,
                    detector_pixel::Real = 0.0, energy = nothing, mode::Symbol = :integrated, times = Float64[],
                    gate::Real = 0.0, gate_shape::Symbol = :gaussian)
    pinhole > 0 && distance > 0 && magnification > 0 && width > 0 && pixels > 0 ||
        throw(ArgumentError("pinhole, distance, magnification, width and pixels must be positive"))
    detector_pixel >= 0 || throw(ArgumentError("detector_pixel must be nonnegative"))
    mode in (:integrated, :gated) || throw(ArgumentError("mode is :integrated or :gated"))
    gate_shape in (:gaussian, :square) || throw(ArgumentError("gate_shape is :gaussian or :square"))
    mode === :gated && (isempty(times) || gate <= 0) &&
        throw(ArgumentError("a gated imager needs times and a positive gate"))
    E = _DANTE_GRID
    f = Pair{Symbol,Union{Float64,Tuple{Float64,Float64}}}[k => _filter_spec(v) for (k, v) in filters]
    tb = collect(TabulatedCurve, tables)
    m = mirror === nothing || mirror isa Union{TabulatedCurve,GrazingMirror} ? mirror :
        (Symbol(mirror[1]), Float64(mirror[2]))
    R = map(E) do e
        t = filter_transmission(f, e / 1000; table = :henke)
        for c in tb
            t *= c(e)
        end
        t * _mirror_reflectivity(m, e) * _detection(detector, e)
    end
    if energy === nothing
        P = [e^3 / expm1(e / 300) for e in E]
        energy = _trapz(E, E .* R .* P) / _trapz(E, R .* P)
    end
    energy > 0 || throw(ArgumentError("energy must be positive"))
    d, M = Float64(pinhole), Float64(magnification)
    λ = 1.239841984e-6 / (energy / 1000)      # mm (hc = 1.2398 keV nm)
    cam = OrthoCamera(-port_vector(theta, phi); up, width, nx = pixels)
    return XrayImager(Float64(theta), Float64(phi), cam, center isa Symbol ? center : Vec3(center...), d,
                      Float64(distance), M, Float64(detector_pixel), E, R, Float64(energy), d * (1 + 1 / M),
                      1.03λ * distance / d, mode,
                      collect(Float64, times), Float64(gate), gate_shape, IdDict{Target,Any}(), ReentrantLock())
end

_detection(::Nothing, e) = 1.0
_detection(c::TabulatedCurve, e) = c(e)
_detection(d::Tuple{Symbol,<:Real}, e) = 1 - filter_transmission([d[1] => Float64(d[2])], e / 1000; table = :henke)


"""
    pinhole_solid_angle(im::XrayImager) -> Float64

The solid angle (sr) of the pinhole seen from the object, `π d² / 4 distance²`: a
detector pixel of side `p` records `L × (p/M)² × Ω` of an image of detected radiance `L`.
"""
pinhole_solid_angle(im::XrayImager) = π * im.pinhole^2 / 4 / im.distance^2

# The element each pixel sees first (0: none), on target `t` (cached per target).
function _imager_pixels(im::XrayImager, t::Target)
    lock(im.lock) do
        get!(im.cache, t) do
            c = im.center isa Symbol ? part_centroid(t, im.center) : im.center
            k = im.camera
            cam = OrthoCamera(k.dir; up = k.up, center = c, width = k.width, nx = k.nx)
            img = render(t, Float64.(1:nelements(t)), cam; background = 0.0)
            Int32.(img)
        end
    end
end

# Planck-weighted means of R over [lo, hi) eV at temperatures `T` (eV)
function _planck_means(im::XrayImager, lo, hi, T)
    E, R = im.E, im.R
    out = similar(T, Float64)
    for (j, θ) in enumerate(T)
        num = den = 0.0
        for i in 1:length(E)-1
            e = (E[i] + E[i+1]) / 2
            (lo <= e < hi) || continue
            w = (E[i+1] - E[i]) * e^3 / expm1(e / θ)
            num += w * (R[i] + R[i+1]) / 2
            den += w
        end
        out[j] = den > 0 ? num / den : 0.0
    end
    return out
end

const _IMAGER_T = exp.(range(log(1.0), log(5000.0), 600))

# ⟨R⟩(T) over [lo, hi) as a function, by a table linear in log ⟨R⟩ and log T (⟨R⟩ of a
# band far above T falls as e^{−E/T}: linear in ⟨R⟩ it is 1 % off between nodes)
function _mean_response(im::XrayImager, lo = 0.0, hi = Inf)
    tab = _planck_means(im, lo, hi, _IMAGER_T)
    lT = log.(_IMAGER_T)
    return function (T)
        T > 0 || return 0.0
        x = clamp(log(T), lT[1], lT[end])
        i = clamp(searchsortedlast(lT, x), 1, length(lT) - 1)
        f = (x - lT[i]) / (lT[i+1] - lT[i])
        a, b = tab[i], tab[i+1]
        return a > 0 && b > 0 ? exp((1 - f) * log(a) + f * log(b)) : (1 - f) * a + f * b
    end
end

"""
    detected_radiance(im::XrayImager, B::AbstractVector) -> Vector{Float64}
    detected_radiance(im::XrayImager, sol::RadiationSolution) -> Vector{Float64}
    detected_radiance(im::XrayImager, s::SpectralSolution) -> Vector{Float64}

The radiance of each element that the imager detects (TW mm⁻² sr⁻¹ of detected x-ray
energy): `B ⟨R⟩ / π`, the element's Lambertian radiance times the imager's response
averaged over its spectrum. For a gray radiosity `B` (TW/mm²) the spectrum is the
Planckian at the element's brightness temperature `(B/σ)^{1/4}` (an assumption: the
reflected part of the radiosity carries the incident spectrum); for a
`SpectralSolution`, each group's radiosity with the Planckian of the element's wall
temperature within the group.
"""
function detected_radiance(im::XrayImager, B::AbstractVector)
    mean = _mean_response(im)
    return [b > 0 ? b / π * mean(radiation_temperature(b)) : 0.0 for b in B]
end
detected_radiance(im::XrayImager, sol::RadiationSolution) = detected_radiance(im, sol.B)

function detected_radiance(im::XrayImager, s::SpectralSolution)
    L = zeros(size(s.B, 1))
    for g in 1:ngroups(s)
        mean = _mean_response(im, s.edges[g], s.edges[g+1])
        for e in eachindex(L)
            s.B[e, g] > 0 && (L[e] += s.B[e, g] / π * mean(s.temperature[e]))
        end
    end
    return L
end

# The point-spread function on the object-plane pixels: the penumbra's disk (each
# pixel's covered fraction, 8 × 8 subsamples) convolved with the Gaussian of diffraction
# and the detector pixel; normalized to sum 1
function _psf(im::XrayImager)
    h = im.camera.width / im.camera.nx
    σ = sqrt((im.diffraction / 2.3548)^2 + (im.detector_pixel / im.magnification)^2 / 12) / h
    ρ = im.penumbra / 2 / h
    r = ceil(Int, ρ + 4σ + 1)
    disk = zeros(2r + 1, 2r + 1)
    for i in -r:r, j in -r:r
        disk[i+r+1, j+r+1] = count(((i + (a - 4.5) / 8)^2 + (j + (b - 4.5) / 8)^2 <= ρ^2 for a in 1:8, b in 1:8)) / 64
    end
    sum(disk) == 0 && (disk[r+1, r+1] = 1.0)              # a penumbra below the pixel
    g = σ < 1e-3 ? [x == 0 ? 1.0 : 0.0 for x in -r:r] : [exp(-x^2 / (2σ^2)) for x in -r:r]
    k = _convolve(_convolve(disk, g; dims = 1), g; dims = 2)
    return k ./ sum(k)
end

# `img` convolved with the 1D kernel `w` (odd length, centered) along `dims`, the same size
function _convolve(img::AbstractMatrix, w::AbstractVector; dims)
    r = length(w) ÷ 2
    out = zeros(size(img))
    n1, n2 = size(img)
    for j in 1:n2, i in 1:n1, q in -r:r
        a, b = dims == 1 ? (i + q, j) : (i, j + q)
        (1 <= a <= n1 && 1 <= b <= n2) && (out[i, j] += w[q+r+1] * img[a, b])
    end
    return out
end

# Blur an object-plane image by the point-spread function (zero outside the field).
function _blur(im::XrayImager, img::AbstractMatrix)
    k = _psf(im)
    r = size(k, 1) ÷ 2
    n1, n2 = size(img)
    out = zeros(n1, n2)
    Threads.@threads for j in 1:n2
        for i in 1:n1
            s = 0.0
            for q in -r:r, p in -r:r
                a, b = i + p, j + q
                (1 <= a <= n1 && 1 <= b <= n2) || continue
                s += k[p+r+1, q+r+1] * img[a, b]
            end
            out[i, j] = s
        end
    end
    return out
end

"""
    imager_image(im::XrayImager, t::Target, L::AbstractVector) -> Matrix{Float64}
    imager_image(im::XrayImager, sol) -> Matrix{Float64}

The image on the object plane (row 1 at the top, `up`) of the per-element detected
radiance `L` ([`detected_radiance`](@ref)), blurred by the imager's resolution:
TW mm⁻² sr⁻¹ of detected x-ray energy. For a solution, its target and
`detected_radiance(im, sol)`.
"""
function imager_image(im::XrayImager, t::Target, L::AbstractVector)
    P = _imager_pixels(im, t)
    return _blur(im, [p > 0 ? Float64(L[p]) : 0.0 for p in P])
end
imager_image(im::XrayImager, sol::RadiationSolution) = imager_image(im, sol.model.target, detected_radiance(im, sol))
imager_image(im::XrayImager, s::SpectralSolution) = imager_image(im, s.gray.model.target, detected_radiance(im, s))

# The weight of each history step in the integral of the piecewise-linear signal
# against gate `g` (the trapezoid rule's weights for g ≡ 1)
function _gate_weights(time::AbstractVector, g)
    n = length(time)
    w = zeros(n)
    for i in 1:n-1
        a, b = time[i], time[i+1]
        m = 64
        for q in 1:m
            τ = a + (q - 0.5) / m * (b - a)
            f = (τ - a) / (b - a)
            gv = g(τ) * (b - a) / m
            w[i] += (1 - f) * gv
            w[i+1] += f * gv
        end
    end
    return w
end

_gate(shape, t0, fwhm) = shape === :square ? (τ -> abs(τ - t0) <= fwhm / 2 ? 1.0 : 0.0) :
                         (τ -> exp(-4log(2) * (τ - t0)^2 / fwhm^2))

"""
    imager_record(im::XrayImager, t::Target, time, L::AbstractMatrix) -> NamedTuple

What the imager records of a history: `L` the detected radiance per element (rows) at
each time of `time` (ns, columns; [`detected_radiance`](@ref) of each step), taken as
linear between the times. `(; images, times)`: for `mode = :integrated` one image, the
time integral (kJ mm⁻² sr⁻¹, as TW·ns), `times` empty; for `:gated` one frame per gate
time, the integral against the gate (peak 1; kJ mm⁻² sr⁻¹), `times` the gate times.
"""
function imager_record(im::XrayImager, t::Target, time::AbstractVector, L::AbstractMatrix)
    size(L, 2) == length(time) || throw(ArgumentError("one column of L per time"))
    issorted(time) || throw(ArgumentError("times must increase"))
    gates = im.mode === :integrated ? [τ -> 1.0] : [_gate(im.gate_shape, t0, im.gate) for t0 in im.times]
    P = _imager_pixels(im, t)
    images = map(gates) do g
        w = _gate_weights(time, g)
        Lw = L * w
        _blur(im, [p > 0 ? Lw[p] : 0.0 for p in P])
    end
    return (; images, times = im.mode === :integrated ? Float64[] : copy(im.times))
end
