# # Integrating sphere
#
# An integrating sphere is a hollow sphere with a diffusely reflecting interior and a
# small port. Light admitted through the port and scattered from any point of the wall
# illuminates every other point of the wall equally: for two points on a sphere of
# radius ``R`` a distance ``r`` apart, ``\cos\theta_i = \cos\theta_j = r/2R``, so the
# view-factor kernel ``\cos\theta_i\cos\theta_j/\pi r^2`` is the constant
# ``1/4\pi R^2``. Ulbricht built photometers on this property in 1900, and it makes the
# sphere an exact test of a radiosity calculation.
#
# A hohlraum is a cavity of the same kind: its gold wall re-emits most of the x rays it
# absorbs, and its laser entrance holes let radiation escape. This tutorial heats a spot
# on the wall of a spherical cavity with a laser admitted through a port, solves for the
# x-ray flux, and compares it with the integrating-sphere solution. The comparison
# shows why Briareus integrates radiation over the true curved surface rather than over
# the flat facets of the mesh. The tutorial introduces materials, apertures, curved
# surfaces, view factors, the radiosity solve and the energy balance.

using Briareus
using CairoMakie
using LinearAlgebra
using Printf
CairoMakie.activate!(type = "png", px_per_unit = 2)
include(joinpath(pkgdir(Briareus), "docs", "literate", "_plotting.jl"))   # figure helpers

# ## The sphere multiplier
#
# Let the wall absorb laser power ``P`` and convert a fraction ``\eta`` of it into x
# rays, and let it re-emit a fraction ``\alpha`` (the albedo) of the x rays that reach
# it. Every point of the wall sees the whole interior with the same kernel, so the
# incident flux ``H`` is the same everywhere: the total power radiated by the wall,
# ``\Phi``, divided by the area ``A`` of the whole sphere, port included. The wall
# radiates its source and re-emits what it absorbs, ``\Phi = \eta P + \alpha H A_w``,
# where ``A_w = (1 - f)A`` is the wall area and ``f`` the port's fraction. Hence
#
# ```math
# H = \frac{\eta P}{A}\,\frac{1}{1 - \alpha (1 - f)} ,
# ```
#
# the direct flux ``\eta P/A`` times the sphere multiplier ``1/[1 - \alpha(1 - f)]``,
# independent of where the laser deposits its power; ``H f A`` escapes through the port.

# ## Geometry
#
# The sphere has a radius of 2 mm. [`cubed_sphere_mesh`](@ref) divides it into
# triangles projected from the faces of a cube, with inward normals because the interior
# radiates. [`split_mesh`](@ref) cuts out the patches within 30° of ``+z``, and
# [`aperture`](@ref) makes them a transparent port whose incident flux is the flux that
# escapes. The wall is a [`Material`](@ref) with albedo 0.5 and conversion efficiency
# 0.7.
#
# A component can carry the exact surface its mesh approximates. With a
# [`SphereSurface`](@ref), the radiation model integrates over the curved patches (true
# points, normals and areas); without one, it integrates over the flat triangles. Ray
# casting uses the triangles in both cases. Both versions are built here.

R = 2.0
α, η = 0.5, 0.7
function cavity(n; curved = true)
    sphere = cubed_sphere_mesh(R; n, outward = false)
    port, wall = split_mesh(sphere, c -> c[3] > R * cosd(30))
    surface = curved ? SphereSurface(R) : nothing
    return Target(Component(:wall, wall, Material(:wall; albedo = α, conversion = η); surface),
                  aperture(:port, port; surface))
end
curved_target = cavity(16)
flat_target = cavity(16; curved = false)

# ## Heating a spot
#
# A beam enters through the port from directly above and heats a spot on the lower
# hemisphere. The reference ray trace delivers it to the wall, apart from a few rays of
# the far wings that graze the rim of the port and strike the outside of the sphere.

beam = Beam("B"; theta = 0, phi = 0, power = 1.0, spot = omega_spot(:SG4), fnumber = 6.7,
            aim = (0.3, 0.0, -1.0))

function sphere_solution(target)
    dep = raytrace_deposition(target, [beam]; nspot = 61)
    wall = elements(target, :wall)
    A = sum(target.elem_area)
    f = sum(target.elem_area[elements(target, :port)]) / A
    H_exact = η * sum(dep.element[wall]) / A / (1 - α * (1 - f))
    model = RadiationModel(target)
    return (; target, dep, wall, A, f, H_exact, model, solution = solve(model, dep))
end
curved = sphere_solution(curved_target)
flat = sphere_solution(flat_target)
@printf("laser power on the wall %.6f TW\n", sum(curved.dep.element[curved.wall]))

# ## Uniformity
#
# [`RadiationModel`](@ref) computes the view factors between all elements, and
# [`solve`](@ref) returns the radiosity ``B`` and the incident flux ``H`` of every
# element. The relative error of the incident flux on each element of the wall, against
# its angular distance from the center of the heated spot, on a logarithmic scale that
# holds both calculations. Elements where the curved calculation is exact to the last
# bit are drawn at ``10^{-17}``:

function angle_from_spot(s)
    center = normalize(sum(s.dep.element .* s.target.elem_centroid))
    return [acosd(clamp(dot(normalize(s.target.elem_centroid[e]), center), -1, 1)) for e in s.wall]
end
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], xlabel = "angle from the heated spot (degrees)",
          ylabel = "|H/H_exact − 1|", title = "Incident flux on the wall", yscale = log10,
          yticks = LogTicks(-17:2:-1), xticks = 0:30:180, limits = ((0, 180), (4e-18, 0.4)))
for (s, label, marker, markersize) in ((flat, "flat facets", :circle, 4),
                                       (curved, "curved surface", :diamond, 5))
    rel_error = abs.(s.solution.H[s.wall] ./ s.H_exact .- 1)
    scatter!(ax, angle_from_spot(s), max.(rel_error, 1e-17); marker, markersize, label)
end
hlines!(ax, [eps()], color = :gray40, linestyle = :dash, label = "machine ε")
Legend(fig[1, 2], ax; framevisible = false)
fig

# On the flat facets the flux is uniform far from the spot but scattered by several
# percent near it. Next to a flat facet, a neighboring facet is tilted by about half an
# element's angular size, while on the sphere the cosines of nearby points are only of
# order ``r/2R``; the facets therefore misstate the exchange between nearby elements by
# an amount of order one, and the spot, which holds most of the source in a few
# elements, makes that error visible in the flux of its neighbors. Refining the facets
# reduces the scatter only in proportion to the element size. On the curved surface
# the flux is uniform to round-off, within 10⁻¹⁴:

for (label, s) in (("flat facets", flat), ("curved surface", curved))
    h = s.solution.H[s.wall] ./ s.H_exact
    @printf("%-15s H/H_exact %.6f – %.6f\n", label, extrema(h)...)
end

# ## Global quantities
#
# The mean flux, the power escaping through the port, and the energy balance from
# [`energy_balance`](@ref):

for (label, s) in (("flat facets", flat), ("curved surface", curved))
    Aw = s.target.elem_area[s.wall]
    b = energy_balance(s.solution)
    @printf("%-15s mean H/H_exact %.6f   escape/(H f A) %.6f   source − sinks %.1e TW\n", label,
            sum(s.solution.H[s.wall] .* Aw) / sum(Aw) / s.H_exact,
            b.escaped[:port] / (s.H_exact * s.f * s.A),
            b.source - sum(values(b.absorbed)) - sum(values(b.escaped)) - b.lost)
end

# The faceted cavity has the right mean to 0.1% because the misstated local exchanges
# nearly cancel in sums. The sphere multiplier grows rapidly as the albedo approaches
# one, and so does the amplification of any radiation the view factors lose. Over a
# range of albedos:

function mean_ratio(s, a)
    sol = solve(s.model, s.dep; albedo = [e in s.wall ? a : 0.0 for e in 1:nelements(s.target)])
    Aw = s.target.elem_area[s.wall]
    exact = η * sum(s.dep.element[s.wall]) / s.A / (1 - a * (1 - s.f))
    return sum(sol.H[s.wall] .* Aw) / sum(Aw) / exact
end
albedos = 0.0:0.05:0.95
ratios = [[mean_ratio(s, a) for a in albedos] for s in (flat, curved)]
multiplier(a) = 1 / (1 - a * (1 - curved.f))
nothing #hide

# The upper panel shows the computed multiplier, the mean wall flux over the direct
# flux ``\eta P/A``, against the exact one; the lower panel shows the relative error of
# the mean flux, with exact results drawn at ``10^{-17}``:

fig = Figure(size = (800, 520))
ax = Axis(fig[1, 1], ylabel = "⟨H⟩ / (ηP/A)", title = "Sphere multiplier",
          xticklabelsvisible = false)
αs = range(0, 0.96; length = 200)
lines!(ax, αs, multiplier.(αs), color = :black, label = "exact 1/[1 − α(1 − f)]")
scatter!(ax, albedos, ratios[1] .* multiplier.(albedos), marker = :rect, markersize = 11,
         label = "flat facets")
scatter!(ax, albedos, ratios[2] .* multiplier.(albedos), marker = :circle, markersize = 6,
         label = "curved surface")
axislegend(ax, position = :lt)
err = Axis(fig[2, 1], xlabel = "albedo α", ylabel = "|⟨H⟩/⟨H⟩_exact − 1|", yscale = log10,
           yticks = LogTicks(-16:4:-2), limits = (nothing, (4e-18, 1e-2)))
scatterlines!(err, albedos, max.(abs.(ratios[1] .- 1), 1e-17), marker = :rect, markersize = 9)
scatterlines!(err, albedos, max.(abs.(ratios[2] .- 1), 1e-17), marker = :circle, markersize = 7)
hlines!(err, [eps()], color = :gray40, linestyle = :dash)
text!(err, 0.5, eps(); text = "machine ε", align = (:center, :top), fontsize = 12,
      color = :gray30)
linkxaxes!(ax, err)
rowsize!(fig.layout, 2, Relative(0.4))
rowgap!(fig.layout, 8)
fig

# The multiplier reaches ``1/[1 - 0.95(1-f)] \approx 9`` at the right of the plot, with
# the port taking ``f = 6.3\%`` of the sphere. The faceted cavity's error grows with it,
# tenfold over these albedos. The curved calculation follows the multiplier to
# round-off; the rows of its view-factor matrix sum to one to round-off, because the
# sphere's kernel is integrated exactly:

for (label, s) in (("flat facets", flat), ("curved surface", curved))
    @printf("%-15s largest |closure error| %.1e\n", label, maximum(abs.(closure_error(s.model.F))))
end

# ## Images
#
# Looking into the sphere from above through the port, an imager sees the lower
# hemisphere of the wall. [`radiance`](@ref) converts radiosity to the Lambertian
# radiance it records, and [`render`](@ref) draws any per-element field as an
# [`OrthoCamera`](@ref) sees it. The radiance shows the heated spot, about four
# elements across. The radiance is computed per element, and the left image smooths
# the element values with a Gaussian 90 µm wide, about half an element, so that the
# facets of this coarse mesh do not hide the spot's shape; the smoothing lowers the
# peak by 7%. The other two images show the error of the incident flux of each element,
# unsmoothed, on one color scale. With flat facets it alternates in sign from triangle
# to triangle and is largest beside the heated spot; on the curved surface it is too
# small to see.

camera = OrthoCamera((0, 0, -1); up = (0, 1, 0), width = 2.1, nx = 420)
B = render_smoothed(curved.target, to_Wcm2.(radiance(curved.solution)) ./ 1e12, :wall, camera;
                    σ = 0.09)
flux_error(s) = render(s.target, 100 .* (s.solution.H ./ s.H_exact .- 1), camera)
fig = Figure(size = (800, 350))
ax = Axis(fig[1, 1]; title = "radiance B/π", aspect = DataAspect())
hidedecorations!(ax)
hm = camera_image!(ax, B, camera; colormap = :inferno, colorrange = (0, maximum(filter(isfinite, B))))
Colorbar(fig[2, 1], hm; vertical = false, flipaxis = false, width = Relative(0.85),
         label = "B/π (10¹² W/cm²/sr)")
for (j, (title, s)) in enumerate((("H, flat facets", flat), ("H, curved surface", curved)))
    local ax = Axis(fig[1, j + 1]; title, aspect = DataAspect())
    hidedecorations!(ax)
    camera_image!(ax, flux_error(s), camera; colormap = :balance, colorrange = (-12.5, 12.5))
end
Colorbar(fig[2, 2:3]; colormap = :balance, colorrange = (-12.5, 12.5), vertical = false,
         flipaxis = false, width = Relative(0.6), label = "H/H_exact − 1 (%)")
rowgap!(fig.layout, 6)
fig

# ## What this checks
#
# - On the curved surface, the radiosity solve reproduces the integrating-sphere flux
#   on every element to round-off, wherever the laser deposits its power, and the
#   sphere multiplier over albedos from 0 to 0.95.
# - Energy is conserved to round-off: the x-ray source equals the power absorbed by
#   the wall, escaping through the port, and lost.
# - Flat facets standing in for the sphere get the mean right to 0.1% but misstate
#   the local exchange near a concentrated source by several percent. This is why
#   hohlraum walls, end caps, entrance holes and capsules carry their exact surfaces.
