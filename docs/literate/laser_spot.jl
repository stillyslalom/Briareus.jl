# # Laser spot on a foil
#
# A phase plate placed before the focusing lens of a high-power laser breaks the beam
# into many small, mutually phased elements. Each element illuminates the whole focal
# region, and their sum, averaged over the speckle, is a smooth spot whose shape is set
# by the plate rather than by the aberrations of the beam. The spots of the continuous
# phase plates on NIF and the distributed phase plates on OMEGA are flat-topped and
# elliptical or round, and an elliptical super-Gaussian describes them well. A foil
# placed in such a beam, as in an area backlighter, absorbs the spot's power
# distribution as it is projected onto the foil's surface.
#
# Three properties of that distribution have exact forms: the total power, the growth
# of the spot as the foil moves away from best focus, and the stretching of the
# footprint at oblique incidence. This tutorial computes the deposition on a foil with
# both of Briareus's methods and compares each property with its exact value. It
# introduces beams, spot profiles, the deposition operator, the reference ray trace and
# rendered images.

using Briareus
using CairoMakie
using LinearAlgebra
using Printf
CairoMakie.activate!(type = "png")

# ## The spot and its defocus
#
# A spot with semi-axes ``a`` and ``b`` and super-Gaussian order ``n`` has the
# best-focus intensity
#
# ```math
# I(u, v) = \frac{P}{N_n\,a b}\,\exp\!\left[-\left(\frac{u^2}{a^2} + \frac{v^2}{b^2}\right)^{n/2}\right],
# \qquad N_n = \frac{2\pi\,\Gamma(2/n)}{n},
# ```
#
# where ``u`` and ``v`` are transverse coordinates and ``P`` is the beam power. Its
# second moment along ``u`` is ``\langle u^2\rangle = a^2 m_n`` with
# ``m_n = \Gamma(4/n)/[2\Gamma(2/n)]``.
#
# Because every point of the lens illuminates the whole spot, the light that reaches a
# point a distance ``s`` from best focus arrives from a disk of focal points of radius
# ``|s|/2F``, where ``F`` is the f-number. In geometric optics the defocused profile is
# therefore the focal profile convolved with a uniform disk. Second moments of a
# convolution add, and a uniform disk of radius ``\rho`` has ``\langle u^2 \rangle =
# \rho^2/4``, so for a round spot the mean squared radius grows as
#
# ```math
# \langle r^2 \rangle(s) = 2 a^2 m_n + \frac{s^2}{8F^2}.
# ```
#
# Briareus evaluates deposition in two ways. A [`DepositionOperator`](@ref) samples
# every lit surface with a quadrature that is fixed by the target and the beam
# directions. It represents the defocused spot by a super-Gaussian of the same order
# whose axes reproduce these second moments, which keeps deposition differentiable in
# pointing, focus and power. [`raytrace_deposition`](@ref) launches rays from a grid
# of focal points along a quadrature of the f-cone and follows them to the first
# surface, which carries out the convolution itself.
#
# The beam here comes from directly above, through OMEGA's SG4 phase plate, which gives
# a round spot with a 1/e radius of 0.36 mm and order 4.2 at f/6.7.

spot = omega_spot(:SG4)
F = 6.7
P = 0.5                                  # TW
m_n = Briareus._sg_moment(spot.order)
moment(s) = 2 * spot.a^2 * m_n + s^2 / (8F^2)
nothing #hide

# The foil is a 6 mm square in the plane ``z = 0``, facing up, divided into
# ``120 \times 120`` patches. [`foil`](@ref) builds the [`Component`](@ref), and
# [`Target`](@ref) assembles components into the geometry the calculations use.

W, N = 6.0, 120
target = Target(foil(:foil, (0, 0, 0), (0, 0, 1); width = W, n = N, two_sided = false))

# ## Power and centroid
#
# A [`Beam`](@ref) is given by its port angles, power, spot, f-number and aim point.
# Here it is aimed off center, at ``(0.3, -0.2, 0)`` mm.

beam = Beam("B"; theta = 0, phi = 0, power = P, spot, fnumber = F, aim = (0.3, -0.2, 0.0))
op = DepositionOperator(target, [beam])
fast = deposit(op, [beam])
ray = raytrace_deposition(target, [beam]; nspot = 61, ncone = 61)
nothing #hide

# The operator's deposition is a vector of absorbed power per element, and its sum is
# the power that reached the foil. The ray trace conserves power exactly; the operator
# integrates the profile by quadrature.

centroid(d) = sum(d.element .* target.elem_centroid) / sum(d.element)
for (label, d) in (("operator", fast), ("ray trace", ray))
    c = centroid(d)
    @printf("%-9s  power %.6f TW   centroid (%.4f, %.4f) mm\n", label, sum(d.element), c[1], c[2])
end

# ## Spot size through focus
#
# Moving the focus a distance ``s`` along the beam leaves the foil at distance ``s``
# from best focus. The second moment of the deposited power is computed from the
# patch centroids. Assigning each patch's power to its centroid adds the variance of a
# uniform distribution over the patch, ``\Delta^2/12`` per axis for a patch of side
# ``\Delta``, which is removed (Sheppard's correction).

Δ = W / N
function spread(d)
    c = centroid(d)
    return sum(d.element .* [sum(abs2, x - c) for x in target.elem_centroid]) /
           sum(d.element) - 2Δ^2 / 12
end

shifts = range(-4.0, 4.0; length = 9)
fast_spread, ray_spread = Float64[], Float64[]
for s in shifts
    b = Beam("B"; theta = 0, phi = 0, power = P, spot, fnumber = F, defocus = -s)
    o = DepositionOperator(target, [b])
    push!(fast_spread, spread(deposit(o, [b])))
    push!(ray_spread, spread(raytrace_deposition(target, [b]; nspot = 61, ncone = 61)))
end

fig = Figure(size = (720, 420))
ax = Axis(fig[1, 1], xlabel = "distance from best focus s (mm)",
          ylabel = "⟨r²⟩ (mm²)", title = "Spot growth through focus")
ss = range(-4, 4; length = 200)
lines!(ax, ss, moment.(ss), color = :black, label = "exact: 2a²mₙ + s²/8F²")
scatter!(ax, shifts, fast_spread, marker = :rect, markersize = 11, label = "deposition operator")
scatter!(ax, shifts, ray_spread, marker = :circle, markersize = 7, label = "ray trace")
axislegend(ax, position = :ct)
fig

# The largest relative errors over the scan:

@printf("operator  %.1e\nray trace %.1e\n",
        maximum(abs.(fast_spread ./ moment.(shifts) .- 1)),
        maximum(abs.(ray_spread ./ moment.(shifts) .- 1)))

# The operator reproduces the exact growth to 0.1% and the ray trace to 0.3%. The ray
# trace's residual comes from sampling the focal plane on a grid; its quadrature of
# the f-cone is constructed to carry the disk's second moment exactly. The operator's
# residual comes from integrating a flat-topped profile over patches and from the
# approximation behind Sheppard's correction.

# ## Oblique incidence
#
# A beam arriving at polar angle ``\theta`` onto the same horizontal foil projects its
# spot onto the surface. Along the plane of incidence the footprint is stretched by
# ``1/\cos\theta``, across it the footprint is unchanged, and the absorbed power per
# unit area is the beam intensity times ``\cos\theta``. The second moments of the
# footprint at best focus are therefore ``\langle x^2 \rangle = a^2 m_n/\cos^2\theta``
# and ``\langle y^2 \rangle = a^2 m_n``, apart from one geometric correction: a tilted
# foil lies above the focal plane on one side of the aim point and below it on the
# other, where the spot is wider. That correction grows as ``1/F^2``, so it appears at
# f/6.7 and vanishes for a nearly collimated f/50 beam. The ray trace here samples the
# spot on a 241-point grid, a sixth of a patch apart: the rays of a nearly collimated
# beam land on a regular lattice, and a lattice comparable to the patch spacing beats
# against it and biases the second moment by about 2%.

θ = 50.0
function footprint(d)
    c = centroid(d)
    w = d.element ./ sum(d.element)
    x2 = sum(w .* [(x[1] - c[1])^2 for x in target.elem_centroid]) - Δ^2 / 12
    y2 = sum(w .* [(x[2] - c[2])^2 for x in target.elem_centroid]) - Δ^2 / 12
    return x2 * cosd(θ)^2 / (spot.a^2 * m_n), y2 / (spot.a^2 * m_n)
end
println("                     ⟨x²⟩cos²θ/a²mₙ   ⟨y²⟩/a²mₙ")
for fnum in (F, 50.0)
    b = Beam("B"; theta = θ, phi = 0, power = P, spot, fnumber = fnum)
    for (label, d) in (("operator", deposit(DepositionOperator(target, [b]), [b])),
                       ("ray trace", raytrace_deposition(target, [b]; nspot = 241, ncone = 61)))
        x, y = footprint(d)
        @printf("f/%-4.1f %-10s      %.4f          %.4f\n", fnum, label, x, y)
    end
end

# The transverse moment matches the projection. The in-plane moment exceeds it by
# 0.7% at f/6.7, identically in the two methods, and agrees with it at f/50.

tilted = Beam("B"; theta = θ, phi = 0, power = P, spot, fnumber = F)
tilted_dep = deposit(DepositionOperator(target, [tilted]), [tilted])
nothing #hide

#
# The peak absorbed flux is the peak beam intensity times ``\cos\theta``:

flux_peak(d) = maximum(laser_flux(d, target))
I0 = focal_profile(spot, P, 0.0)
@printf("peak absorbed flux: normal %.4f, oblique/cos θ %.4f, exact %.4f TW/mm²\n",
        flux_peak(fast), flux_peak(tilted_dep) / cosd(θ), I0)

# ## Images
#
# [`render`](@ref) draws any per-element field as seen by an [`OrthoCamera`](@ref).
# Here the camera looks down onto the foil, and the field is the absorbed laser flux
# from [`laser_flux`](@ref), converted to W/cm² by [`to_Wcm2`](@ref).

camera = OrthoCamera((0, 0, -1); up = (0, 1, 0), width = W, nx = 240)
defocused = Beam("B"; theta = 0, phi = 0, power = P, spot, fnumber = F, defocus = -4.0)
cases = ("best focus" => (beam, fast),
         "4 mm from focus" => (defocused, deposit(DepositionOperator(target, [defocused]), [defocused])),
         "θ = 50°" => (tilted, tilted_dep))
fig = Figure(size = (900, 330))
extent = (-W / 2, W / 2)
for (i, (label, (b, d))) in enumerate(cases)
    img = render(target, to_Wcm2.(laser_flux(d, target)), camera)
    panel = Axis(fig[1, i], title = label, aspect = DataAspect(), xlabel = "x (mm)",
                 ylabel = i == 1 ? "y (mm)" : "")
    heatmap!(panel, range(extent..., length = 240), range(extent..., length = 240),
             reverse(img; dims = 1)'; colormap = :inferno, colorrange = (0, 1.3e15))
end
Colorbar(fig[1, 4], colormap = :inferno, colorrange = (0, 1.3e15), label = "absorbed flux (W/cm²)")
fig

# ## What this checks
#
# - The deposition operator and the ray trace both deliver the beam's power to the
#   foil and center it at the aim point.
# - Through focus, the spot's second moment follows the exact geometric-optics growth
#   ``2a^2m_n + s^2/8F^2``, to 0.1% with the operator and 0.3% with the ray trace. This
#   checks the operator's moment-matched defocus and the ray trace's sampling of the
#   f-cone.
# - At oblique incidence the footprint is stretched by ``1/\cos\theta`` along the
#   plane of incidence, the transverse width is unchanged, and the absorbed flux is
#   reduced by ``\cos\theta``. The defocus across the tilted footprint adds 0.7% to the
#   in-plane moment at f/6.7, the same in both methods, and nothing at f/50.
