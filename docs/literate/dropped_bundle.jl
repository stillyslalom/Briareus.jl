# # Freeing a bundle for radiography
#
# A radiograph of an imploding capsule needs a backlighter, and on NIF the
# backlighter lasers come out of the drive. Freeing one bundle (two quads of four
# beams, one quad in each hemisphere, sharing an amplifier chain) leaves 46 quads to
# drive the hohlraum, and the missing ports leave an asymmetry in the capsule drive. On a
# layered DT ignition shot there is little margin for it. A mode-1 asymmetry of the drive,
# a flux difference between opposite sides of the capsule, pushes the hot spot sideways at
# about 100 km/s for each percent of mode-1 flux (MacGowan et al., High Energy Density
# Physics 40, 100944, 2021), and the drive itself is the scarcest resource of the shot.
#
# Which bundle can be freed is set first by the radiograph: its quads must reach a
# backlighter foil placed across the target from the detector, perpendicular to the line
# of sight, without their light touching the hohlraum or its entrance-hole shields on the
# way, and the x rays from the foil must reach the detector unobstructed. This tutorial
# finds the bundles that can serve each of two equatorial lines of sight, ranks them by
# the symmetry they cost, and compares three ways to win that symmetry back: freeing the
# symmetric partner of the bundle as well, rebalancing the powers of the remaining quads,
# and re-pointing them within the clearance the entrance holes allow. It introduces
# three-dimensional mode analysis with real spherical harmonics, the symmetry of the NIF
# ports, scenarios on a subset of the beams, the influence matrix, optimization with
# automatic-differentiation gradients and under constraints, beam-clearance checks, and a
# re-pointing study in two stages, from quads to their beams.

using Briareus
using CairoMakie
using LinearAlgebra
using Printf
CairoMakie.activate!(type = "png", px_per_unit = 2)
include(joinpath(pkgdir(Briareus), "docs", "literate", "_plotting.jl"))   # figure helpers

# ## Modes on the capsule
#
# Expand the incident flux on the capsule in real spherical harmonics,
#
# ```math
# H(\theta, \phi) = \sum_{l, m} c_{lm}\, Y_{lm}(\theta, \phi),
# ```
#
# where ``Y_{lm}`` with ``m > 0`` varies as ``\cos m\phi`` and with ``m < 0`` as
# ``\sin |m|\phi``; ``Y_{11}`` points along ``+x`` and ``Y_{1,-1}`` along ``+y``.
# [`flux_ylm`](@ref) returns the coefficients, in the order of [`mode_index`](@ref).
# The size of mode ``l`` relative to the mean is
#
# ```math
# \epsilon_l = \frac{1}{c_{00}}\sqrt{(2l + 1)\sum_m c_{lm}^2},
# ```
#
# computed by [`mode_amplitude`](@ref). For a field with azimuthal symmetry about any
# axis this is the familiar Legendre ratio ``|a_l/a_0|``; for ``l = 1`` it is the
# fractional amplitude ``\epsilon`` of a flux ``1 + \epsilon\cos\gamma``, where
# ``\gamma`` is the angle from the direction of the asymmetry.

# ## The hohlraum and its pointing
#
# The hohlraum is 10.1 mm long with a radius of 2.875 mm and 3.1 mm laser entrance
# holes, around a capsule of radius 0.9 mm. Its wall, end caps, entrance holes and
# capsule carry their exact cylindrical, planar and spherical surfaces (see the
# [Integrating sphere](@ref) tutorial). [`nif_quads`](@ref) gives the 48 quads, at 8 TW
# each, with their measured phase-plate spots. The outer quads are aimed on the
# hohlraum axis near the entrance-hole planes, at ``z = \pm 5.01`` mm (50° cone) and
# ``\pm 5.20`` mm (44.5°); the inner quads at ``z = \pm 3.63`` mm (30°) and
# ``\pm 3.09`` mm (23.5°), displaced 0.60 and 1.09 mm from the axis away from their own
# port, across the hohlraum. This pointing is an illustrative example and was not
# fielded on any shot.

target = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  wall = gold(albedo = 0.85, conversion = 0.8))
pointing = Dict(:outer50 => (0.0, 5.01), :outer44 => (0.0, 5.20),
                :inner30 => (0.60, 3.63), :inner23 => (1.09, 3.09))
function aimed(b)
    r, z = pointing[b.group]
    φ = port_angles(b)[2]
    return remake(b; aim = (-r * cosd(φ), -r * sind(φ), sign(b.port[3]) * z))
end
quads = aimed.(nif_quads(power = 8.0))

# A [`Scenario`](@ref) holds what does not change while beams are re-aimed and their
# powers adjusted: the target, the deposition operator for these ports, the view factors,
# the linear response of the incident flux on the capsule to the power each element
# absorbs, and the mode weights of the capsule. Its parameterization here gives each quad
# a separate power.
#
# The deposition operator is sampled more coarsely than its defaults, to keep the
# optimizations below quick. Its points lie at most 0.6 of the smaller spot semi-axis
# apart (the default is 0.3) and are refined only within 0.5 mm of the footprint of each
# beam (the default is 2 mm; no aim moves by more than 0.3 mm here); elsewhere each
# triangle gets one point. That cuts the operator to a fifth of its points and changes the
# modes below in the fourth digit.

scenario = Scenario(target, GroupPointing(quads; by = :quad, fields = (:power,));
                    observe = :capsule, deposition_kw = (; resolution = 0.6, reach = 0.5, level = 0))
full = flux_ylm(evaluate(scenario, quads); lmax = 4)
show_modes(label, c) = @printf("%-24s ε₁ %6.3f%%  ε₂ %6.3f%%  ε₃ %6.3f%%  ε₄ %6.3f%%  drive %.4f\n",
                               label, (100 * mode_amplitude(c, l) for l in 1:4)..., c[1] / full[1])
show_modes("48 quads", full)

# With all 48 quads the drive has the large ``l = 2`` asymmetry of this untuned
# pointing (the cone powers and pointing would tune it away) and no mode 1.

# ## The symmetry of the ports
#
# The NIF ports are not symmetric under every rotation that leaves the hohlraum
# unchanged. Rotating every port by 180° about the ``x``, ``y`` or ``z`` axis maps the
# set of ports onto itself, but rotating by 90° about ``z``, inverting through the
# center, or reflecting through the equatorial plane does not:

ports = [b.port for b in quads]
maps_onto_itself(f) = all(any(norm(f(p) - q) < 1e-6 for q in ports) for p in ports)
operations = ["180° about x" => p -> [p[1], -p[2], -p[3]],
              "180° about y" => p -> [-p[1], p[2], -p[3]],
              "180° about z" => p -> [-p[1], -p[2], p[3]],
              "90° about z" => p -> [-p[2], p[1], p[3]],
              "inversion" => p -> -p,
              "reflection in z = 0" => p -> [p[1], p[2], -p[3]]]
for (name, f) in operations
    @printf("%-20s %s\n", name, maps_onto_itself(f) ? "symmetry" : "not a symmetry")
end

# The three rotations form the group ``D_2``. A harmonic survives the average over
# ``D_2`` only if each rotation leaves it unchanged. A rotation by 180° about ``z``
# multiplies ``Y_{lm}`` by ``(-1)^m``, and one about ``x`` multiplies it by ``(-1)^{l+m}``
# for ``m \ge 0`` and by ``(-1)^{l+m+1}`` for ``m < 0``. A drive with the symmetry of the
# ports can therefore contain only even ``m``, with ``l`` even for the cosine harmonics
# and ``l`` odd for the sine harmonics. Every coefficient the rule forbids is zero within
# round-off:

allowed(l, m) = iseven(m) && (m >= 0 ? iseven(l) : isodd(l))
forbidden = [abs(full[mode_index(l, m)]) for l in 0:4 for m in -l:l if !allowed(l, m)]
@printf("largest forbidden |c_lm|/c₀₀ %.1e\n", maximum(forbidden) / full[1])
for l in 1:4, m in -l:l
    allowed(l, m) && @printf("c(%d,%2d)/c₀₀ = %+.5f\n", l, m, full[mode_index(l, m)] / full[1])
end

# The cosine ``l = 2, m = 2`` and sine ``l = 3, m = -2`` harmonics are allowed by the
# symmetry of the ports, and the ports excite them weakly.

# ## Which bundles can drive the backlighter
#
# Two lines of sight at the equator (``\theta = 90°``) suit radiography of the
# capsule, at ``\phi = 78°`` and ``\phi = 348°``. For each, [`backlighter_foil`](@ref)
# places a 2 mm square foil on the far side of the target, a distance ``d`` (the
# standoff) from its center, facing the detector across the capsule; the radiography
# beams strike its face toward the target, so the x rays they generate cross the
# capsule to the detector. The size and standoff of the foil here are generic engineering
# values: real backlighter dimensions are shot-specific. The radiograph looks through a
# pair of 3 mm windows in the hohlraum wall, one on each side of the capsule along the
# line of sight.
#
# The scene around the foil is the hohlraum with its entrance-hole shields
# (`sombrero = true`), the dimpled annular sheets around each end that keep
# unconverted 1ω light off the hohlraum and its thermo-mechanical package. Each
# reaches 9.2 mm from the axis, 0.62 mm inside the entrance-hole plane (see
# [`Sombrero`](@ref)), so it stands across much of the space between the ports and a
# foil near the waist. The shields lie outside the cavity and do not change its
# radiation, and the cost of the windows to the drive is the subject of the
# [Diagnostic windows in the hohlraum wall](@ref) tutorial, so the scenario above
# uses the closed hohlraum without them; the scene is for the geometry only.

shielded = Dict(φ => hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                              wall = gold(albedo = 0.85, conversion = 0.8), sombrero = true,
                              windows = (WallWindow(; phi = φ, diameter = 3.0, open = true, name = :detector_window),
                                         WallWindow(; phi = mod(φ + 180, 360), diameter = 3.0, open = true,
                                                    name = :foil_window)))
                for φ in (78.0, 348.0))
scene(φ_los, d; hohl = shielded[φ_los]) =
    Target([hohl.components..., backlighter_foil((90.0, φ_los); standoff = d)])

# A quad can drive the foil if it strikes the foil face at less than 70° from its
# normal, the LMF limit for laser capture on a target surface, and if its beam, aimed
# at the foil center, clears the hardware on the way: the hohlraum and capsule, and
# the shields. The clearance is judged as a shot review judges it, on the worst case:
# [`clearance`](@ref) grows the envelope of each beam (its contour at 1% of the peak,
# spread by the f-cone) by 50 µm of pointing error and 25 µm of target positioning,
# adds ±1 mm of focus error, and reports the smallest distance to each piece of
# hardware. That distance is exact; only the power on the hardware, not used here,
# depends on how finely the error set and the spot are sampled, so both are sampled
# coarsely. The standoff is a design choice: for each bundle, the table gives the
# worst-case clearance from the hohlraum and from the shields at 8 and 12 mm, and the
# smallest standoff at which both are positive, found by bisection to 0.05 mm.

hohlraum_parts, shields = (:wall, :cap_top, :cap_bottom, :capsule), (:sombrero_top, :sombrero_bottom)
function worst_clearance(φ_los, k, d)
    n = port_vector(90.0, φ_los)                           # toward the detector
    report = clearance(scene(φ_los, d), [remake(b; aim = Tuple(-d .* n)) for b in quads[bundles[k]]];
                       targets = (:backlighter,), hardware = (hohlraum_parts..., shields...),
                       peak_fraction = 0.01, wings = [], nangles = 4, nspot = 15)
    worst(parts) = minimum(r -> r.clearance, filter(r -> r.component in parts, report.rows))
    return (hohlraum = worst(hohlraum_parts), shield = worst(shields))
end
clears(φ_los, k, d) = (c = worst_clearance(φ_los, k, d); min(c.hohlraum, c.shield) > 0)
function min_standoff(φ_los, k; lo = 8.0, hi = 30.0, tol = 0.05)
    @assert !clears(φ_los, k, lo) && clears(φ_los, k, hi)
    while hi - lo > tol
        mid = (lo + hi) / 2
        clears(φ_los, k, mid) ? (hi = mid) : (lo = mid)
    end
    return hi
end

# A closed form bounds the bisection from below. The axis of a quad runs from the foil
# center ``-d\,\hat n`` back toward its port ``\hat p``, crossing the shield plane of
# the same hemisphere, ``|z| = z_s`` (4.427 mm), a distance ``s = z_s/|p_z|`` along
# ``\hat p``. It passes outside the rim, ``\rho = 9.2`` mm, once
# ``|-d\,\hat n + s\,\hat p_\perp| > \rho``, that is beyond the larger root
#
# ```math
# d_\mathrm{axis} = s\,\hat n\cdot\hat p_\perp
#   + \sqrt{\rho^2 - s^2\left(|\hat p_\perp|^2 - (\hat n\cdot\hat p_\perp)^2\right)},
# ```
#
# taken here as the larger over the two quads of the bundle. The beam envelope and the
# error set add to it.

z_shield = 10.1 / 2 - Sombrero().offset
function axis_standoff(φ_los, k)
    n = port_vector(90.0, φ_los)
    return maximum(quads[bundles[k]]) do b
        s = z_shield / abs(b.port[3])
        p⊥ = [b.port[1], b.port[2], 0.0]
        s * dot(n, p⊥) + sqrt(Sombrero().outer_radius^2 - s^2 * (dot(p⊥, p⊥) - dot(n, p⊥)^2))
    end
end

bundles = nif_bundles(quads)
standoff = Dict{Tuple{Float64,String},Float64}()
for φ_los in (78.0, 348.0)
    n = port_vector(90.0, φ_los)
    candidates = [k for k in sort(collect(keys(bundles)))
                  if all(b -> acosd(dot(b.port, n)) < 70, quads[bundles[k]])]
    @printf("line of sight 90-%g: %d bundles within 70°\n", φ_los, length(candidates))
    @printf("  %-6s %-16s %-13s %-34s %s\n", "bundle", "cones", "incidence",
            "worst case, hohlraum / shields (mm)", "standoff: axis, worst case (mm)")
    for k in candidates
        cones = join(sort(unique(string(quads[i].group) for i in bundles[k])), "+")
        incidence = [acosd(dot(b.port, n)) for b in quads[bundles[k]]]
        c8, c12 = worst_clearance(φ_los, k, 8.0), worst_clearance(φ_los, k, 12.0)
        standoff[(φ_los, k)] = min_standoff(φ_los, k)
        @printf("  %-6s %-16s %4.1f°, %4.1f°  d = 8: %5.2f / %5.2f, 12: %5.2f / %5.2f   %5.2f, %5.2f\n",
                k, cones, incidence..., c8.hohlraum, c8.shield, c12.hohlraum, c12.shield,
                axis_standoff(φ_los, k), standoff[(φ_los, k)])
        @assert standoff[(φ_los, k)] > axis_standoff(φ_los, k)
    end
end

# Each line of sight can be served by six bundles. Without the shields, a foil 12 mm from
# the center would let all six clear the hohlraum, and at 8 mm the hohlraum would leave
# only bundle 46 for 90-78 and bundle 11 for 90-348, whose quads pass above and below it
# at steep incidence. The shields rule out both designs: at 8 and 12 mm the beam of every
# candidate quad reaches a shield (a worst-case clearance near −1 mm is the whole
# envelope on it), the beam of the upper quad on the upper shield and that of the lower
# quad on the lower. Bundles 46 and 11 are still the first to clear, at a standoff of
# 13.2 mm; the other candidates need 14.5–16.3 mm. Each worst-case standoff lies
# 1.6–1.9 mm beyond the closed-form one at which the axes of the bundle graze the rim. The
# difference is the half-width of the envelope plus the errors, stretched by the shallow
# angle at which the beams cross the shield plane (40° for the 50° cone).
#
# The beams are half of it; the x rays must also reach the detector, and the same hardware
# may stand between the foil and the detector. [`line_of_sight`](@ref) traces a bundle of
# rays from every sample point of the emitting face of the foil through every sample point
# of the projected disk of the capsule, on to the detector, and reports any that hardware
# stops. The traced paths include every x-ray path from the foil through the capsule
# toward a detector in the cone of the bundle, more paths than one imager accepts. Here
# the foil is placed at the standoff that bundle 11 needs, rounded up to a millimeter:

d_los = ceil(standoff[(348.0, "11")])
los = line_of_sight(scene(348.0, d_los), (90.0, 348.0))

# Every one of the 81 × 81 rays passes: the shields stand more than 3 mm above and
# below the equatorial line of sight, the windows open the wall, and the foil covers
# the capsule disk. Without the windows the hohlraum wall stops every ray:

closed_wall = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                       wall = gold(albedo = 0.85, conversion = 0.8), sombrero = true)
line_of_sight(scene(348.0, d_los; hohl = closed_wall), (90.0, 348.0))

# and a foil closer in spreads the bundle more beyond the capsule. At 8 mm the rays from
# the corners of the foil through the rim of the capsule reach about 1.7 mm from the line
# of sight at the far wall, past the 1.5 mm radius of the detector-side window, which
# stops a twentieth of them:

s_far = (8.0 + 2.875) / 8.0           # distance from the foil to the far wall, in standoffs
@printf("bundle radius at the far wall %.2f mm
", s_far * 0.9 + (s_far - 1) * sqrt(2))
line_of_sight(scene(348.0, 8.0), (90.0, 348.0))

# ## The radiograph signal
#
# A clear line of sight still needs enough photons. The freed quads heat the foil. A
# fraction of their energy, the conversion efficiency (CE), comes out in the He-α line
# of the foil material, and the imager records what passes its pinhole and filters.
# [`kshell_ce`](@ref) estimates the CE from published measurements on planar foils
# ([`kshell_ce_data`](@ref)): OMEGA foils driven by a single 1 ns pulse (Workman and
# Kyrala, 2001), NIF foils with a picket 3 ns ahead of a 2–5 ns drive (Barrios et al.,
# 2013 and 2014), and OMEGA foils with and without a pre-pulse (Coppari et al., 2019).
# [`kshell_ce_model`](@ref) fits `ln CE`, weighting each point by its stated error, as a
# quadratic in `ln I` for each line, with an optimum irradiance that rises with photon
# energy and a peak that falls as ``(h\nu)^{-3.3}`` above about 6 keV and levels off
# below. The scatter of the data about the fit, a factor of 1.43, is the CE's
# uncertainty. The data do not separate the pulse
# length from the picket or the source.
#
# [`backlighter_signal`](@ref) computes the chain in six steps:
# 1. the laser energy on each foil element, from a deposition calculation;
# 2. the CE at the mean irradiance of the foil;
# 3. isotropic emission (the assumption the CE measurements were converted with);
# 4. a pinhole imager on the line of sight;
# 5. the attenuation of the filters ([`filter_transmission`](@ref), from NIST mass
#    attenuation coefficients);
# 6. the photons per detector pixel, with the CE's uncertainty.
#
# The capsule itself does not attenuate in this calculation: the result is the open-field
# level that its transmission multiplies. Here bundle 11's two quads, 16 TW for 2 ns
# with a picket, are aimed at the center of a germanium foil (He-α at 10.3 keV) at the
# 14 mm standoff. The imager is generic: a 20 µm pinhole 200 mm from the capsule, the
# detector 2 m away (magnification 9), behind 250 µm of beryllium and 50 µm of Kapton.

foil_center = Tuple(-d_los .* port_vector(90.0, 348.0))
radiography = [remake(b; aim = foil_center) for b in quads[bundles["11"]]]
bl_scene = scene(348.0, d_los)
bl_deposition = raytrace_deposition(bl_scene, radiography; nspot = 61)
imager = (; distance = 2000.0, pinhole = 0.02, pinhole_distance = 200.0, pixel = 0.05,
          filters = (:Be => 0.25, :kapton => 0.05))
signal = backlighter_signal(bl_scene, bl_deposition, (90.0, 348.0); line = :Ge, pulse = 2.0,
                            prepulse = true, imager...)

# The phase-plate spots at best focus spread the light of the quads over most of the foil,
# so the mean irradiance, about ``5\times10^{14}`` W/cm², lies just below the lowest
# measured NIF point, which flags the result as extrapolated. It is also well below the
# fitted optimum for germanium, about ``2\times10^{15}`` W/cm². About 450 photons per
# 50 µm pixel (6 µm at the capsule), averaged over the capsule image, is uncertain by a
# factor of 1.5 either way. Across the image the level varies a hundredfold
# (`signal.photons_range`), because the radiance of the foil follows the spots of the
# quads: the rim of the capsule is backlit by the weakly irradiated edges of the foil. All
# the emitting area of the foil faces the detector, and through the open windows every
# element sees it, so the power per steradian toward the detector is the isotropic
# ``E_x/4\pi``. Lower-Z lines convert better at this irradiance; the transmission of the
# capsule at each energy determines the choice between them:

@assert signal.facing ≈ 1 && signal.per_sr ≈ signal.xray_energy / 4π
println("line   hν (keV)  CE (%, 1σ range)      filters  photons per pixel (1σ range)")
for line in (:Cu, :Zn, :Ge, :Zr)
    s = backlighter_signal(bl_scene, bl_deposition, (90.0, 348.0); line, pulse = 2.0,
                           prepulse = true, imager...)
    @printf("%-4s   %6.2f   %5.2f (%4.2f–%4.2f)%s   %5.3f    %6.0f (%.0f–%.0f)\n", line, s.line,
            100s.ce, 100s.ce * exp(-s.sigma), 100s.ce * exp(s.sigma), s.extrapolated ? "*" : " ",
            s.transmission, s.photons_pixel, s.photons_lo * 0.05^2, s.photons_hi * 0.05^2)
end

# (A `*` marks a value outside the measured range of its class, here in irradiance.)

# ## What each bundle costs
#
# Without a laser–plasma model the capsule flux is linear in the beam powers.
# [`influence_matrix`](@ref) returns the harmonic coefficients produced by each quad
# alone at 1 TW, so the coefficients for any assignment of powers ``P`` are ``M P``.
# [`nif_bundles`](@ref) groups the quads into bundles.

M = influence_matrix(scenario, quads; lmax = 4)
P_full = [b.power for b in quads]
without(k) = setdiff(eachindex(quads), bundles[k])

# The rotations of ``D_2`` carry each bundle onto another, so the 24 bundles fall into
# sets of equivalent bundles, and freeing any member of a set costs the same. Each
# quad maps under a rotation onto a quad of the image bundle:

image(i, f) = findfirst(q -> norm(f(quads[i].port) - q.port) < 1e-6, quads)
rotations = last.(operations[1:3])
image_bundle(k, f) = only(unique(nif_bundle(quads[image(i, f)]) for i in bundles[k]))
orbits = unique(sort(unique([k; [image_bundle(k, f) for f in rotations]])) for k in keys(bundles))

cost(k) = M[:, without(k)] * P_full[without(k)]
sort!(orbits; by = o -> mode_amplitude(cost(o[1]), 1))
@printf("%-13s %-15s %7s %7s %7s %7s %7s\n", "bundles", "cones", "ε₁ (%)", "ε₂ (%)",
        "ε₃ (%)", "ε₄ (%)", "drive")
for o in orbits
    cones = join(sort(unique(string(quads[i].group) for i in bundles[o[1]])), "+")
    ε = [[mode_amplitude(cost(k), l) for l in 1:4] for k in o]
    @assert all(e -> isapprox(e, ε[1]; atol = 1e-9), ε)
    @printf("%-13s %-15s %7.3f %7.3f %7.3f %7.3f %7.4f\n", join(o, " "), cones,
            (100 .* ε[1])..., cost(o[1])[1] / full[1])
end

# The two quads of a bundle need not share a cone. The bundles whose quads both sit in the
# 50° cone cost least, in mode 1 and in drive, and those of the 44.5° cone next. A
# bundle with an inner quad costs three times as much mode 1: the inner quads deposit
# their power around the waist of the hohlraum, closest to the capsule. Among the
# bundles that can drive a backlighter, the 50° bundles are therefore the ones to free:
# 11 or 23 for the line of sight at 90-348, 46 or 11 for 90-78. Of these, 11 (for
# 90-348) and 46 (for 90-78) are also the ones the shields let through closest in, at
# 13.2 mm; the others need 15.4 mm. This tutorial frees bundle 11, the pair Q11T and Q11B,
# for the line of sight at 90-348; bundle 46 is its image under a rotation of ``D_2``
# and costs the same.

freed = "11"
keep = without(freed)
remaining = quads[keep]

# ## Freeing bundle 11
#
# A [`subset`](@ref) of the scenario contains the 46 remaining quads and reuses everything
# except the deposition operator, which it restricts to their ports. The same coefficients
# follow from the full scenario with the quads of the bundle at zero power, and from the
# influence matrix:

c_free = flux_ylm(evaluate(subset(scenario, keep), remaining); lmax = 4)
zeroed = [i in bundles[freed] ? remake(b; power = 0.0) : b for (i, b) in enumerate(quads)]
c_zeroed = flux_ylm(evaluate(scenario, zeroed); lmax = 4)
c_matrix = M[:, keep] * P_full[keep]
@printf("subset vs zero power %.1e   subset vs influence matrix %.1e\n",
        norm(c_free - c_zeroed) / norm(c_free), norm(c_free - c_matrix) / norm(c_free))
show_modes("bundle 11 freed", c_free)

# The drive falls by 4% and a mode-1 asymmetry of 1.1% appears, about 110 km/s of
# hot-spot velocity. Its direction lies between the azimuths of the freed quads'
# ports:

dipole = normalize([c_free[mode_index(1, 1)], c_free[mode_index(1, -1)], c_free[mode_index(1, 0)]])
@printf("mode-1 direction θ = %.1f°, φ = %.1f°\n", acosd(dipole[3]),
        mod(atand(dipole[2], dipole[1]), 360))
for i in bundles[freed]
    @printf("%s port θ = %.2f°, φ = %.2f°\n", quads[i].name, port_angles(quads[i])...)
end

# An outer quad crosses the hohlraum axis near its entrance hole and heats the wall on
# the far side, so without it the side of the capsule facing away from its port is
# colder.
#
# The maps below show the flux on the capsule as seen from outside, in an equal-area
# (Hammer) projection, so that equal areas of the map are equal solid angles of the
# capsule. Polar angle θ is marked at the left edge and azimuth φ along the equator,
# which runs across the middle. The flux is computed per element, on the 1,620
# triangles of the geodesic capsule mesh; the maps interpolate it linearly from the
# element values averaged onto their corners. The left map shows the flux relative to the
# mean flux, ``H/\langle H\rangle - 1``, where ``\langle H\rangle = c_{00}/\sqrt{4\pi}``.
# The right map subtracts the same quantity for all 48 quads, which isolates what
# freeing the bundle changed. The figure helpers `sphere_axis`, `sphere_map` and
# `mark_directions!` come from `docs/literate/_plotting.jl`.

capsule = elements(target, :capsule)
deviation(H, c) = 100 .* (H ./ (c[1] / sqrt(4π)) .- 1)       # percent of ⟨H⟩
H_full = evaluate(scenario, quads).flux
H_free = evaluate(subset(scenario, keep), remaining).flux
port_marker = (marker = :xcross, markersize = 12, color = :black)
los_marker = (marker = :circle, markersize = 11, color = :transparent, strokecolor = :black,
              strokewidth = 1.5)
function capsule_panel(pos, field, range; title = "")
    ax = sphere_axis(pos; title)
    hm = sphere_heatmap!(ax, sphere_map(target, field, :capsule); colormap = :balance,
                         colorrange = (-range, range))
    mark_directions!(ax, [port_angles(quads[i]) for i in bundles[freed]]; port_marker...)
    mark_directions!(ax, [(90, 348)]; los_marker...)
    return hm
end
capsule_legend(pos) = Legend(pos, [MarkerElement(; port_marker...), MarkerElement(; los_marker...)],
                             ["ports of the freed quads", "line of sight 90-348"];
                             orientation = :horizontal, framevisible = false)

fig = Figure(size = (800, 330))
hm = capsule_panel(fig[1, 1], deviation(H_free, c_free), 15; title = "bundle 11 freed")
Colorbar(fig[2, 1], hm; vertical = false, flipaxis = false, width = Relative(0.7),
         label = "H/⟨H⟩ − 1 (%)")
hm = capsule_panel(fig[1, 2], deviation(H_free, c_free) - deviation(H_full, full), 2;
                   title = "change from 48 quads")
Colorbar(fig[2, 2], hm; vertical = false, flipaxis = false, width = Relative(0.7),
         label = "difference in H/⟨H⟩ − 1 (%)")
capsule_legend(fig[3, 1:2])
rowgap!(fig.layout, 6)
fig

# The untuned ``l = 2`` asymmetry, a hot waist and cold poles, dominates the flux and
# barely changes. What freeing the bundle changed is mostly mode 1: the side facing the
# freed ports, around ``\phi = 45°``, is warmer relative to the mean by about 1%, and
# the far side, where the freed quads deposited their power, colder by up to 2%.

# ## Remedy 1: free the partner as well
#
# The rotation by 180° about ``z`` carries bundle 11 onto bundle 34. Freeing both
# leaves ports with that symmetry, which removes every harmonic of odd ``m``, among
# them the two transverse components of mode 1. The axial component ``c_{10}`` vanishes
# for another reason: the hohlraum is axisymmetric, so the contribution of a quad to
# ``c_{10}`` depends only on its cone and hemisphere, and each freed bundle takes one
# quad from each hemisphere of the same cone.

partner = image_bundle(freed, rotations[3])
keep_pair = setdiff(keep, bundles[partner])
c_pair = M[:, keep_pair] * P_full[keep_pair]
show_modes("bundles 11 and $partner freed", c_pair)

# Mode 1 is gone, but so is another 4% of the drive, which an ignition shot cannot
# spare.

# ## Remedy 2: rebalance the powers
#
# Keeping the other 46 quads at their pointing, find the smallest change ``\delta P``
# of their powers that restores the full-set coefficients of every mode up to degree
# ``L``. With ``M_L`` the rows of ``M`` through degree ``L``, this is the
# minimum-norm solution
#
# ```math
# \delta P = M_L^{+}\,\left(c^{\mathrm{full}}_L - M_L P\right),
# ```
#
# where ``M_L^{+}`` is the pseudoinverse:

Mk = M[:, keep]
P_keep = P_full[keep]
rebalanced = Dict{Int,Vector{Float64}}()
for L in 1:3
    rows = 1:(L+1)^2
    rebalanced[L] = P_keep + pinv(Mk[rows, :]) * (full[rows] - Mk[rows, :] * P_keep)
    show_modes("rebalanced through l = $L", Mk * rebalanced[L])
    @printf("%24s quad powers %.2f – %.2f TW\n", "", extrema(rebalanced[L])...)
end

# Restoring the drive and mode 1 needs some quads 11% above the nominal 8 TW, and
# restoring ``l = 3`` as well needs 24%. An ignition pulse runs the laser close to its
# energy and peak-power limits, so this headroom is rarely available. The quads near
# the freed bundle carry the extra power. Each hemisphere is shown separately, with
# the nominal 8 TW as a gray line and the azimuth of its freed quad dashed; the three
# solutions are offset by ±1.5° in azimuth so that they do not hide one another:

φk = [port_angles(b)[2] for b in remaining]
upper = [b.port[3] > 0 for b in remaining]
function rebalance_panel(pos, hemisphere, sel; bottom)
    ax = Axis(pos; ylabel = "power (TW)", title = "quads in the $hemisphere hemisphere",
              xticks = 0:45:360, limits = ((0, 360), (6.8, 10.3)),
              xlabel = "port azimuth φ (degrees)", xlabelvisible = bottom,
              xticklabelsvisible = bottom)
    hlines!(ax, [8.0]; color = :gray, linewidth = 1)
    for i in bundles[freed]
        b = quads[i]
        (b.port[3] > 0) == (hemisphere == "upper") || continue
        vlines!(ax, [port_angles(b)[2]]; color = :gray, linestyle = :dash)
        text!(ax, port_angles(b)[2] + 2, 10.2; text = "$(b.name) freed", align = (:left, :top),
              fontsize = 12, color = :gray30)
    end
    for (L, marker, dx) in ((1, :circle, -1.5), (2, :utriangle, 0.0), (3, :rect, 1.5))
        scatter!(ax, φk[sel] .+ dx, rebalanced[L][sel]; marker, markersize = 9,
                 label = "restored through l = $L")
    end
    return ax
end
fig = Figure(size = (800, 470))
rebalance_panel(fig[1, 1], "upper", upper; bottom = false)
ax = rebalance_panel(fig[2, 1], "lower", .!upper; bottom = true)
Legend(fig[3, 1], ax; orientation = :horizontal, framevisible = false)
rowgap!(fig.layout, 8)
fig

# The same answer comes from the optimizer. [`optimize_pointing`](@ref) minimizes an
# objective of an [`Evaluation`](@ref) over the parameters of the scenario, here the 46
# quad powers, with gradients from forward-mode automatic differentiation. With the
# misfit of the ``l \le 2`` coefficients as the objective and the bounds out of reach,
# L-BFGS-B starts from ``P`` and moves only within the span of the gradients, which is
# the row space of ``M_2``; it therefore converges to the minimum-norm solution:

power_only = subset(scenario, keep; parameterization = GroupPointing(remaining; by = :quad,
                                                                     fields = (:power,)))
misfit(ev) = sum(abs2, flux_ylm(ev; lmax = 2) .- full[1:9]) / full[1]^2
P0 = initial_parameters(power_only.param)
fit = optimize_pointing(power_only, misfit; x0 = P0, lower = fill(0.0, length(P0)),
                        upper = fill(20.0, length(P0)), g_tol = 1e-14)
@printf("misfit %.1e → %.1e,  max |P_optimizer − P_pseudoinverse| %.1e TW\n",
        fit.initial_value, fit.value, maximum(abs.(fit.x - rebalanced[2])))

# ## Remedy 3: re-point
#
# Re-aiming quads moves drive around the hohlraum at no cost in power. Let each of the
# remaining quads move its aim by up to 0.3 mm along the hohlraum axis and sideways
# around it, and minimize mode 1 strongly and modes 2–4 weakly while holding the drive
# at its value after the bundle is freed. [`asymmetry_objective`](@ref) builds that
# objective and adds penalties for light clipped by the entrance holes and light
# striking the capsule. Each optimization here stops after a fixed number of L-BFGS-B
# iterations, 40 and then 20, short of full convergence but well past the point where
# mode 1 is gone.

repoint = subset(scenario, keep; parameterization = GroupPointing(remaining; by = :quad,
                                                                  fields = (:phi, :z)))
x0 = initial_parameters(repoint.param)
lower, upper = bounds(repoint.param, x0; phi = 0.3, z = 0.3)
goals = (modes = (1 => 200.0, 2 => 1.0, 3 => 1.0, 4 => 1.0), drive = c_free[1] / sqrt(4π),
         drive_weight = 10.0)
first_try = optimize_pointing(repoint, asymmetry_objective(; goals...); x0, lower, upper,
                              iterations = 40)
show_modes("bundle 11 freed", c_free)
show_modes("re-pointed", flux_ylm(evaluate(repoint, first_try.x); lmax = 4))
@printf("clipped %.1e of the delivered power\n", clipped_fraction(evaluate(repoint, first_try.x)))

# Mode 1 is gone and nothing is clipped. But the clipping penalty sees only the
# nominal beams, and a shot review judges a pointing on its worst case.
# [`clearance`](@ref) grows the envelope of each beam (here its contour at 1% of the peak,
# the level of LMF's target-configuration rules, spread by the f-cone) by 50 µm of
# pointing error and 25 µm of target positioning, adds ±1 mm of focus error, and reports
# the smallest distance to the end caps, together with the worst-case power and
# intensity the caps receive. The scattered-light wings, which reach the
# entrance-hole lip in any design, are left out:

function cap_review(beams)
    report = clearance(target, beams; targets = (:wall, :capsule),
                       hardware = (:cap_top, :cap_bottom), peak_fraction = 0.01, wings = [],
                       nangles = 8)
    return (nominal = minimum(r -> r.nominal, report.rows),
            worst = minimum(r -> r.clearance, report.rows),
            power = maximum(r -> r.intercepted, report.rows),
            fail = count(r -> r.status === :fail, report.rows))
end
review_line(label, r) = @printf("%-22s nominal %6.3f mm   worst case %6.3f mm   cap power ≤ %.1e   %d fail\n",
                                label, r.nominal, r.worst, r.power, r.fail)
review_48 = cap_review(quads)
review_line("48 quads", review_48)
review_line("re-pointed", cap_review(first_try.beams))

# The original pointing clears the caps by 66 µm in the worst case. The re-pointed
# quads have used up that margin: in the worst case their 1% contours overlap the
# entrance-hole lips, several of them fail, and up to several times 10⁻⁴ of the
# power of a beam is deposited on a cap.
#
# The optimizer has to see the same worst case. An [`ApertureKeepOut`](@ref) computes
# the worst-case clearance of each beam from the rim of its entrance hole with the same
# envelope and error set, smoothly enough to differentiate. The objective adds a penalty
# when the clearance falls below a minimum, here 10 µm, with the weight given beside the
# keep-out, here 10³. On the original pointing the keep-out agrees with
# [`clearance`](@ref):

keepout = ApertureKeepOut(target; peak_fraction = 0.01, minimum = 0.01)
from_keepout = aperture_clearance(keepout, quads)
@printf("worst-case clearance from the keep-out %.3f mm, from clearance %.3f mm\n",
        minimum(from_keepout), review_48.worst)

second_try = optimize_pointing(repoint, asymmetry_objective(; goals..., keepout = keepout => 1e3); x0,
                               lower, upper, iterations = 20)
repointed = evaluate(repoint, second_try.x)
c_repointed = flux_ylm(repointed; lmax = 4)
show_modes("re-pointed, keep-out", c_repointed)
review_line("re-pointed, keep-out", cap_review(second_try.beams))

# Re-pointing within the clearance still removes nearly all of mode 1, to a few kilometers
# per second of hot-spot velocity, and the drive stays within 0.3% of its value after the
# bundle is freed. The worst-case clearance comes out at 8 µm, short of the
# 10 µm minimum: the penalty is soft, and the distance from the keep-out differs from that
# of the clearance check by a few micrometers. The maps below show the flux on the capsule
# before and after re-pointing, on the same color scales as above:

fig = Figure(size = (800, 470))
for (j, (title, H, c)) in enumerate((("bundle 11 freed", H_free, c_free),
                                     ("re-pointed within clearance", repointed.flux, c_repointed)))
    capsule_panel(fig[1, j], deviation(H, c), 15; title)
    capsule_panel(fig[2, j], deviation(H, c) - deviation(H_full, full), 2)
end
Colorbar(fig[1, 3]; colormap = :balance, colorrange = (-15, 15), height = Relative(0.8),
         label = "H/⟨H⟩ − 1 (%)")
Colorbar(fig[2, 3]; colormap = :balance, colorrange = (-2, 2), height = Relative(0.8),
         label = "change from 48 quads (%)")
capsule_legend(fig[3, 1:2])
rowgap!(fig.layout, 6)
fig

# Re-pointing removes the warm side around the freed ports. What remains of the change
# from the 48-quad drive is in the higher modes, which the objective weighted two
# hundred times less than mode 1: the re-pointed drive trades the dipole for a
# stronger ``l = 3``, while ``l = 2`` is as it was with all 48 quads.

# ### Independent check
#
# The optimizer worked through the deposition operator and the linear flux response.
# [`raytrace_deposition`](@ref) recomputes the deposition of the re-pointed quads by
# tracing rays from each spot through its f-cone, and [`solve`](@ref) recomputes the
# radiosity directly:

traced = raytrace_deposition(target, second_try.beams; nspot = 61)
H = solve(scenario.model, traced).H[capsule]
show_modes("re-pointed, ray-traced", ylm_weights(target, capsule; lmax = 4) * H)

# ## Re-pointing within limits
#
# The objective above weighs mode 1 against the other modes and the drive by hand-set
# weights. A design review states limits instead: every mode within a tolerance of the
# 48-quad shape, the laser irradiance on the wall no higher than under the 48-quad
# pointing, and within those limits as much drive as possible.
#
# The first limit is a gate per degree on the departure from the 48-quad shape,
# ``\delta_l = \sqrt{(2l+1)\sum_m (c_{lm}/c_{00} - c'_{lm}/c'_{00})^2}``
# ([`mode_departure`](@ref)), for an axisymmetric difference the change in the Legendre
# ratio ``a_l/a_0`` (``P_2/P_0`` at ``l = 2``). [`mode_uncertainty`](@ref) gives the rms
# amplitude that random pointing errors of 50 µm per quad leave in each degree, and
# [`mode_gate`](@ref) turns it into the gate ``\sqrt{2l+1}`` times that rms. That gate
# is a lower bound: it is what pointing error alone produces, while the delivery of the
# laser and the hohlraum plasma move the shape from shot to shot by more, in ways a
# view-factor model does not quantify. The usual gate is a tolerance per degree in
# physical units, a choice of the design; the floor is used here as the strictest case,
# and the research notes ("Mode gates in physical units") measure what relaxed gates buy
# on this case in aim shift and drive. The second limit is a ceiling on the laser
# irradiance over the wall elements ([`irradiance_ceiling`](@ref)): a hot spot on the gold
# drives the wall plasma inward. The peak under the 48-quad pointing serves as the limit,
# an engineering choice; no published criterion sets it. The peak over the elements
# ([`peak_irradiance`](@ref)) has a kink wherever the hottest element changes, so the
# optimizer works on a smooth stand-in, the area-weighted 64-norm of the irradiance
# ([`irradiance_norm`](@ref)), scaled to equal the peak at the reference. Once a run has
# converged on the stand-in, the ceiling compares the peak itself with the limit and
# lowers the working limit of the stand-in while the peak exceeds it (`exact = true`, its
# default).

ref48 = evaluate(scenario, quads)
gates = mode_gate(full, mode_uncertainty(scenario, quads; pointing = 0.05, lmax = 4))
ceiling = irradiance_ceiling(ref48, :wall)
check_constraints(evaluate(repoint, x0), (ceiling, gates))

# Freeing the bundle puts every degree outside its gate, l = 1 by a factor of three,
# while the hottest wall element, lit by quads that stay, remains at the 48-quad peak:
# the ceiling is active from the start.
#
# [`optimize_pointing`](@ref) takes the limits as `constraints`. It handles them by an
# augmented Lagrangian: a sequence of L-BFGS-B runs on the objective plus
# multiplier-weighted terms, each run started where the last ended, with the multipliers
# updated in between until every limit is met. The objective is the drive
# ([`drive_objective`](@ref)), with the usual penalties on light clipped by the
# entrance holes.

most_drive = drive_objective(; drive = mean_flux(ref48))
limited = optimize_pointing(repoint, most_drive; x0, lower, upper, constraints = (ceiling, gates),
                            iterations = 1000, evaluations = 3000)

# One call from the start meets every gate and the ceiling. The start violates the gates,
# so the run first searches for a design within the limits by minimizing only the
# violations (phase 1, `phase1` on by default) and starts its main run there. The table
# lists the value of each limit, its normalized excess and its multiplier; the outcome
# records that the run converged with every constraint met inside the box. The multipliers
# are those that best satisfy stationarity at the design, and the stationarity line gives
# the residual they leave. Here the residual is not small against the drive gradient: the
# drive is flat along the limits and the run ends on its value tolerance, so the
# multipliers rank the cost of the limits only roughly.
#
# The same limits as penalty terms are what a hand-written objective would add: 100
# times the squared excess of each degree over its gate and of the irradiance over the
# ceiling, with a small reward, 0.01 times the squared ratio of each degree to its gate, for
# symmetry inside the gates.

function penalized(ev)
    J = most_drive(ev)
    ratio = mode_departure(ev, full) ./ gates.limit
    J += sum(100 * max(0, r^2 - 1)^2 + 0.01 * r^2 for r in ratio)
    excess = ceiling.f(ev) / only(ceiling.limit) - 1
    return J + 100 * max(0, excess)^2
end
by_penalty = optimize_pointing(repoint, penalized; x0, lower, upper, iterations = 1000, evaluations = 3000)
for (label, r) in (("within limits", limited), ("penalty terms", by_penalty))
    ev = evaluate(repoint, r.x)
    @printf("%-14s drive %.4f  δ/gate %s  peak/limit %.4f  %4d evaluations\n", label,
            flux_ylm(ev; lmax = 0)[1] / full[1],
            join((@sprintf("%.3f", d) for d in mode_departure(ev, full) ./ gates.limit), " "),
            peak_irradiance(ev, :wall) / peak_irradiance(ref48, :wall), r.evaluations)
end
@assert limited.outcome in (:interior, :box) && all(r -> r.status !== :violated, limited.constraints)
@assert peak_irradiance(evaluate(repoint, limited.x), :wall) <= only(ceiling.limit) * (1 + 1e-3)

# Both designs meet the gates. The penalty terms keep rewarding symmetry inside the
# gates, beyond what the limits require, and pay for it with drive: modes 1 and 3 end
# far inside their gates and three tenths of a percent of the drive is gone. With the
# limits as constraints the gates end at or near their limits, with more drive left
# over. The drive is nearly flat along the limits, so where the run ends in that family
# of designs depends on its path, down to rounding: runs started 10⁻¹² to 10⁻⁹ mm from
# this one end with other gates active, after 840 to 1,500 evaluations, at 0.9616 to
# 0.9625 of the 48-quad drive (research notes, "The dropped-bundle tutorial's long runs
# under rounding (2026-11-20)"). The printed departures of this run are one member of
# the family; its drive is determined to about 0.1%. The penalty design is determined:
# the same perturbed starts reproduce its drive and departures to the printed digits and
# its peak to 0.1%, and only its evaluation count varies. Staging the penalties by hand
# (gates loosened four-, two- and 1.4-fold, then the final ones) reaches the penalty
# design again, in more evaluations (research notes, "Constrained re-pointing").
#
# The constrained run ends with the peak over the elements at the ceiling or under it
# (from 0.8% under to 0.06% over, within the run tolerance of 0.1%). The penalty form
# bounds the 64-norm, a mean over the hot spots, and leaves the hottest element a percent
# over the ceiling. Without the ceiling the same gates are met with the wall peak 27–38%
# above the peak under the 48-quad pointing (research notes), because the re-aimed quads
# are stacked on the same spots.

# ## A re-pointing study, end to end
#
# The design within limits meets the shape and wall-irradiance limits, but no limit in its
# run covered the entrance holes. The clearance check of the end caps shows what that run
# left unconstrained:

review_line("within limits", cap_review(limited.beams))

# In the worst case the 1% contours of the quads end within a few tens of micrometers of
# the caps, on either side: no limit in the run held them, and over the perturbed starts
# above, the worst case ranges from 38 µm into a cap to 4 µm clear of it, with up to two
# beams failing the check. A study repeats the run with every limit a design review
# checks, on the parameters the design may change, in two stages: the quads first, then
# each of their beams.
#
# ### The parameters
#
# Here the inner quads keep their pointing (their spots set the drive at the waist and
# pass close to the capsule), the outer quads move by up to 0.3 mm as before, and the
# power of each cone may fall by up to 20% but not rise: an ignition pulse runs the laser
# at its limit. [`PowerScales`](@ref) adds a power scale per cone after the pointing
# parameters, and [`HeldParameters`](@ref) fixes at their values the parameters that a
# predicate on their names selects, here those of the inner quads. With `relative = true`
# the parameters of each quad are shifts from its present aim. [`remake`](@ref) puts the
# composed parameterization in the scenario and reuses its deposition operator, since the
# beams enter through the same ports.

inner_quads = Set(b.quad for b in remaining if b.group in (:inner23, :inner30))
by_quad = GroupPointing(remaining; by = :quad, fields = (:phi, :z), relative = true)
study = remake(repoint; parameterization = HeldParameters(PowerScales(by_quad),
                                                          n -> first(split(n, '.')) in inner_quads))
lower_s, upper_s = bounds(study.param; phi = 0.3, z = 0.3, scale = 0.2, scale_max = 1.0)
println(nparams(study.param), " free parameters, ending ", join(parameter_names(study.param)[end-4:end], ", "))

# ### The limits
#
# The study has four limits: the gates and the wall-irradiance ceiling above, a keep-out
# from the end caps, and a ceiling on the light the caps receive.
# [`HardwareKeepOut`](@ref) gives the worst-case clearance of each beam from named parts,
# here the caps, exactly and with the envelope and error set of the clearance check, and
# [`clearance_floor`](@ref) makes the minimum of the keep-out, 10 µm, a lower bound on the
# clearance of each beam, as a constraint with one component and one multiplier per beam.
# The keep-out checks the 1% contour; [`light_ceiling`](@ref) bounds all the nominal light
# on the caps, the wings of the spots included ([`light_fraction`](@ref), absorbed plus
# striking their outside), here at 10⁻⁵ of the delivered power. On the 48-quad pointing
# the keep-out and the clearance check agree:

cap_keepout = HardwareKeepOut(target, (:cap_top, :cap_bottom); targets = (:wall, :capsule),
                              peak_fraction = 0.01, minimum = 0.01)
@printf("worst-case clearance from the caps: keep-out %.4f mm, clearance check %.4f mm\n",
        minimum(hardware_clearance(cap_keepout, quads)), review_48.worst)
limits = (ceiling, gates, clearance_floor(cap_keepout), light_ceiling(1e-5, :cap_top, :cap_bottom))
stage1 = optimize_pointing(study, most_drive; lower = lower_s, upper = upper_s, constraints = limits,
                           iterations = 1000, evaluations = 3000)

# The run ends with every limit met and the box active: the 44.5° cone at full power,
# the 30° cone 1.5–3.3% below it and the other two within 1%, and one to three of
# the four gates active. Every quad stays above the 10 µm minimum of the keep-out (the
# display shows the worst of the 46 rows, one per quad, each with a separate multiplier,
# and counts the active ones), and the light ceiling is met with room to spare: the gates
# and the box are the binding limits in this run, and the cap limits are not. This run too
# ends in a flat family: from the perturbed starts above it ends at 0.950 to 0.957 of the
# 48-quad drive, with the least keep-out clearance at 13 to 66 µm, after 650 to 1,650
# evaluations.

ev1 = evaluate(study, stage1.x)
@printf("%s after %d evaluations; cone powers %s of nominal\n", stage1.outcome, stage1.evaluations,
        join((@sprintf("%s %.3f", n[1:7], x) for (n, x) in zip(stage1.names[end-3:end], stage1.x[end-3:end])), ", "))
gate_line(label, ev) = @printf("%-22s drive %.4f, departure / gate for l = 1–4: %s\n", label,
                               mean_flux(ev) / mean_flux(ref48),
                               join((@sprintf("%.3f", d) for d in mode_departure(ev, full) ./ gates.limit), " "))
gate_line("stage 1 (quads)", ev1)
review_line("stage 1 (quads)", cap_review(stage1.beams))

# ### From quads to beams
#
# So far each quad has been one beam: one aim, one power and the deposition of its four
# beams as a single f/7.9 cone. [`nif_beams`](@ref) gives the four beams of each quad,
# each through its port at f/20.7, at the aim of the quad and with a quarter of its power.
# A scenario on them needs a deposition operator for their ports: `remake` builds one,
# given the options of the first. The polish lets each outer beam move independently by up
# to 0.1 mm, along the hohlraum axis and along the tangent of its quad (a quad is the unit
# of [`GroupPointing`](@ref)'s frame by default); the beams of the inner quads are held as
# before.

beams = nif_beams(stage1.beams)
inner_beams = Set(b.name for b in beams if b.quad in inner_quads)
by_beam = GroupPointing(beams; by = :beam, fields = (:phi, :z), relative = true)
polish = remake(study; parameterization = HeldParameters(by_beam, n -> first(split(n, '.')) in inner_beams),
                deposition_kw = (; resolution = 0.6, reach = 0.5, level = 0))
polish_start = evaluate(polish, initial_parameters(polish.param))
gate_line("stage 1 on beams", polish_start)
check_constraints(polish_start, limits)

# The same aims, deposited per beam, give a slightly different drive: the four beams of a
# quad enter at angles up to a few degrees apart, each with a separate, narrower f-cone,
# and a departure that was at its gate can now exceed it by up to a couple of percent
# (l = 1 or l = 3, depending on where the first stage ended). The limits are the same
# objects as in the first stage. The components of the keep-out follow the beams of the
# evaluation, 184 of them now, so the constraints carry over unchanged. The polish starts
# from the aims of the first stage (`relative = true`, all shifts zero) with a capped
# budget:

lower_p, upper_p = bounds(polish.param; phi = 0.1, z = 0.1)
stage2 = optimize_pointing(polish, most_drive; lower = lower_p, upper = upper_p, constraints = limits,
                           iterations = 1000, evaluations = 120)

# ### The result
#
# The report of the run above gives the value of each constrained quantity against its
# limit, with its multiplier and status. [`pointing_summary`](@ref) puts the polished
# beams beside the beams of the quads before re-pointing (`nif_beams(remaining)`, matched
# by name): the aim shifts, the power and cone fraction, the least keep-out clearance
# against that of the reference, and the light clipped and missed. Together they record
# the design:

pointing_summary(polish, stage2; reference = nif_beams(remaining), keepout = cap_keepout)

#-

ev2 = evaluate(polish, stage2.x)
@printf("%s after %d evaluations\n", stage2.outcome, stage2.evaluations)
gate_line("polished, per beam", ev2)
@printf("peak wall irradiance %.4f of the ceiling\n", peak_irradiance(ev2, :wall) / only(ceiling.limit))
polished_review = cap_review(stage2.beams)
review_line("polished", polished_review)
held = filter(r -> r.name != "irradiance", stage2.constraints)   # gates, keep-out, light
@assert stage2.outcome in (:interior, :box, :evaluations) && all(r -> r.violation <= 0.01, held)
@assert ceiling.f(ev2) <= only(ceiling.limit) * 1.01   # the smooth stand-in of the ceiling
@assert peak_irradiance(ev2, :wall) <= only(ceiling.limit) * 1.01   # the peak itself
@assert polished_review.fail == 0 && polished_review.worst > 0

# The polish brings the departures back to their gates. It stops at its budget, and a run
# stopped by its budget with a limit beyond its tolerance returns its design of least
# violation, here its best design within every limit (`least_violation` in the result).
# Six starts were compared: the start of this tutorial and five moved by 10⁻¹² or 10⁻⁹ mm,
# which stand for another CPU's rounding. Two runs end within the limits and keep the
# drive they gained, 0.17 and 0.32% over the aims of the first stage on beams. On the
# other four the last design exceeds a gate by more than the tolerance, and the run
# returns an earlier one within the limits, at the drive it started from. Every limit then
# holds: the gates, the keep-out, and the wall peak at 0.93 to 0.995 of the ceiling. The
# returned design is chosen on the reported rows, which contain the peak itself, while the
# run constrains the smooth stand-in of the ceiling. The ceiling constrains the peak
# itself only once a run has converged; before this choice the peak ended 0.966 to 1.020
# of the ceiling over thirteen such starts. The budget is kept small here to keep the
# tutorial short. A design to hand off runs to convergence, which on an earlier first
# stage took 861 evaluations (five and a half minutes on 8 threads) and met every limit,
# with the peak at 0.997 of the ceiling and the same drive (research notes, "Optimization
# defaults decided"). The clearance check confirms the keep-out on the worst case: no beam
# fails, and the closest stays 17 to 49 µm from a cap.
#
# Each stage was one call of [`optimize_pointing`](@ref) on exported names: a
# parameterization, its box and the limits, with the scenario of the second stage built by
# `remake` on the beams of the first stage. Both results stay at hand for a report or a
# further stage.

# ## A tolerance in physical units
#
# The gates above sit at the pointing-noise floor, the strictest case. A study states a
# tolerance per degree instead, a choice made for the design: 0.5% in ``P_2/P_0``
# is the modelers' repeat uncertainty for one ignition design (research notes, "Drive
# asymmetry inferred from post-shot capsule simulations"), and degrees 1, 3 and 4 have no
# published number. Here, as an example, every degree is held within 0.5% of the
# 48-quad shape. [`shift_objective`](@ref), the mean squared aim shift from the start,
# with the ceiling and the keep-out from the end caps as further constraints, gives the
# least re-pointing that meets the tolerance, in one call:

tolerance = mode_gate(full, fill(0.005, 4))
smallest = optimize_pointing(repoint, shift_objective(repoint); x0, lower, upper,
                             constraints = (ceiling, tolerance, clearance_floor(cap_keepout)),
                             iterations = 1000, evaluations = 3000)
smallest_summary = pointing_summary(repoint, smallest; keepout = cap_keepout, drive_reference = full)

#-

@assert smallest.outcome in (:interior, :box) && all(r -> r.status !== :violated, smallest.constraints)
@assert sqrt(smallest.value) ≈ smallest_summary.shift_rms

# Degrees 1–3 end at the tolerance and l = 4 well inside it; the quads move by 63 µm rms
# and 0.105 mm at most, and the drive ends 0.3% below that of the start (the objective
# does not reward drive). The same call with the gates at the noise floor moves the quads
# by 93 µm rms and up to 0.18 mm; with 1% per degree, which the departures of the start in
# l = 1 and 2 barely exceed, it moves them by 6 µm. The largest drive within each
# tolerance, from
# [`drive_objective`](@ref) under the same constraints, rises from 0.962 of the 48-quad
# drive at the floor to 0.967 at 0.5% and 0.971 at 1% (research notes, "Mode gates in
# physical units", measured before phase 1 was the default).

# ## What this checks
#
# - The NIF ports have the symmetry ``D_2`` and no more, and every harmonic that
#   ``D_2`` forbids vanishes from the computed drive within round-off.
# - Bundles in the same ``D_2`` orbit cost the same in every mode, as the symmetry
#   requires.
# - Three routes to the drive without a bundle agree within round-off: a scenario on the
#   remaining beams, the full scenario with the bundle at zero power, and the
#   influence matrix.
# - Freeing a bundle together with its image under the 180° rotation about the
#   hohlraum axis removes mode 1.
# - The power rebalance from the pseudoinverse and from the automatic-differentiation
#   optimizer agree, as they must for a linear least-squares problem solved from the
#   same start.
# - The worst-case standoff of every candidate bundle from the shields lies beyond the
#   closed-form standoff at which the axes of its quads graze the shield rim (asserted).
# - The reach of the line-of-sight bundle at the far wall with the foil at 8 mm, from the
#   closed form for rays from the corners of the foil through the rim of the capsule,
#   exceeds the radius of the detector-side window, and [`line_of_sight`](@ref) finds rays
#   stopped there and none at 14 mm. Its closed-form tests (`test/backlighter.jl`) cover a
#   bundle exactly grazing the rim of a shield, a foil rotated off the line of sight, and
#   imaging rays converging on a point detector.
# - The power per steradian of the backlighter toward the detector equals the isotropic
#   ``E_x/4\pi``: every emitting element faces the detector and sees it through the
#   windows (asserted). The closed forms of the chain are tested in
#   `test/backlighter_signal.jl`: point-projection fluence ``E_x/(4\pi r^2)``, the uniform
#   image ``E_x A_p/(\pi A q^2)`` of a Lambertian foil, energy through the pinhole equal
#   to the fluence times the image area, the hemispherical integrals of the angular
#   distributions (``E_x/2`` isotropic, ``E_x`` Lambertian), Beer–Lambert filters against
#   the NIST nodes and edges, and the weighted residuals, optimum and ``h\nu`` scaling of
#   the CE fit against the published ``(h\nu)^{-(2.5\text{–}4)}``. The form fitted without
#   Coppari et al.'s 31 points predicts them within their stated 20% and the scatter of
#   the fit, and the fit predicts an OMEGA Ti foil below the fitted lines within 1.5σ
#   (both held out).
# - The worst-case clearance from the keep-out agrees with that from the clearance check,
#   and re-pointing against it passes that check.
# - The re-pointed drive, recomputed by ray tracing and a direct radiosity solve,
#   reproduces the modes from the optimizer.
# - Re-pointing within limits meets every gate and the ceiling in one call (asserted, with
#   no constraint violated, the run converged and the peak over the wall elements within
#   0.1% of the ceiling or below it). In its closed-form test (`test/constraints.jl`),
#   with the coefficients of the drive linear in the quad powers, a fit under a ceiling on
#   the total power matches the solution of the Karush–Kuhn–Tucker equations to 1e-4 and
#   its multiplier to 1e-3. On a drive-like objective linear in the powers under two
#   active ceilings, the reported multipliers match those of the closed form to 1e-3, also
#   where the estimate of the iteration has zeroed one, and a bound sensitivity matches
#   the closed form and restarts with the bound moved. The same file tests the power-mean
#   ordering of the irradiance norm, its area-mean limit at p = 1, and the gradients of
#   the constraint values against central differences.
# - The keep-out of the study from the end caps agrees with the worst case of the
#   clearance check on the 48-quad pointing to within a micrometer, and the polished beams
#   pass that check with no beam failing (asserted). The budget-stopped polish meets the
#   limits it works on: the normalized violation of every gate, the keep-out and the light
#   ceiling is at most 0.01 (1% of a gate, 10 µm of the keep-out), and the smooth stand-in
#   of the ceiling is at most 1% over its limit (asserted; over 13 starts, perturbed by
#   10⁻¹² to 10⁻⁹ mm or by rounding, the largest of those violations was 0.002 and the
#   stand-in at most 0.994 of its limit). Only a converged run constrains the wall peak
#   itself, so the peak is not asserted here. In the closed-form tests of the keep-out and
#   light constraints (`test/constraints.jl`), one free parameter, the aim of the 50°
#   cone, is pulled by the objective into the caps, and a constrained run ends with the
#   constraint active at the aim where its value reaches the limit, found by bisection, to
#   1e-5 mm. The tests also check the values of the constraints against the quantities
#   that define them (the light on the caps against `component_power` to 1e-12) and their
#   gradients against central differences to 1e-6.
# - The least re-pointing within a tolerance of 0.5% per degree meets the tolerance,
#   the ceiling and the keep-out, and converges (asserted, with no constraint violated), and
#   the square root of its objective equals the rms aim shift of the summary, computed
#   from the aims directly (asserted). The closed-form test of the objective
#   (`test/pointing_summary.jl`) checks the least weighted re-pointing for one linear
#   constraint against its solution in closed form, to 1e-5 mm.
