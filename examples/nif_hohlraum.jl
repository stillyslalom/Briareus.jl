# Capsule drive symmetry in a NIF hohlraum: evaluate an illustrative cone pointing, then
# move the inner and outer rings along the axis and re-balance their power to null P2
# and P4 at the same mean flux, clear of the entrance holes and the capsule (the
# README's first calculation, with a ray-trace check, an image and the energy balance).
#
#   julia --project=. -t auto examples/nif_hohlraum.jl

using Briareus

# 5.75 mm × 10.1 mm Au hohlraum (the Build-A-Hohlraum diameter), 3.1 mm LEHs, capsule
# at TCC.
target = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9)
display(target)

# 48 quads at 8 TW; cone pointing in hohlraum (r, z) coordinates (mm). Positive r points
# toward the quad's own port azimuth, so the inner cones' aims across the axis are
# negative. These aims are an illustrative example, not a facility or shot pointing
# (research notes, "Example NIF cone pointing: unsourced (2026-10-27)").
cones = GroupPointing(nif_quads(power = 8.0); fields = (:z, :r))
quads = apply(cones, parameters(cones; inner23 = (z = 3.09, r = -1.09),
                                inner30 = (z = 3.63, r = -0.60),
                                outer44 = (z = 5.20,), outer50 = (z = 5.01,)))

# One axial aim shift and one power per ring (inner: 23.5° and 30°; outer: 44.5° and
# 50°), mirrored top to bottom: four parameters, fixed at the optimum by P2, P4, the mean
# flux and the entrance-hole keep-out.
ring(b) = b.group in (:inner23, :inner30) ? :inner : :outer
pointing = GroupPointing(quads; by = ring, fields = (:z, :power), relative = true)
scenario = Scenario(target, pointing; observe = :capsule)   # closed: the hohlraum default
x0 = initial_parameters(pointing)
println("\nExample pointing, 8 TW per quad:")
report(evaluate(scenario, x0))

# Hold the drive at its initial level, without clipping the LEH or hitting the capsule.
a0 = mean_flux(evaluate(scenario, x0))
keepout = (ApertureKeepOut(target), SphereKeepOut(target))
objective = symmetry_objective(modes = (2 => 1.0, 4 => 1.0); drive = a0, keepout)
lower, upper = bounds(pointing, x0; power = 0.6)     # ±1 mm of aim, ±60 % of power
result = optimize_pointing(scenario, objective; x0, lower, upper)
display(result)   # lists any parameter left at a bound

println("\nOptimized ($(round(result.initial_value; sigdigits = 3)) → ",
        "$(round(result.value; sigdigits = 3))):")
for (name, a, b) in zip(parameter_names(pointing), x0, result.x)
    println(rpad(name, 16), round(a; digits = 3), "  →  ", round(b; digits = 3))
end
report(evaluate(scenario, result.x))

# Independent check: brute-force ray trace + direct radiosity solve at the optimum.
rt = raytrace_deposition(target, result.beams)
sol = solve(scenario.model, rt)
a = legendre_moments(sol.H, target, elements(target, :capsule); lmax = 6)
println("\nRay-trace check: ", join(("P$l/P0 = $(round(100a[l+1]/a[1]; digits = 2))%" for l in 2:2:6), ", "))

# X-ray image down the hohlraum axis through the upper LEH, and the energy balance.
img = render(target, radiance(sol), OrthoCamera((0, 0, -1); up = (0, 1, 0), width = 6.0, nx = 200))
println("LEH image: ", count(isfinite, img), " lit pixels, peak radiance ",
        round(maximum(filter(isfinite, img)); sigdigits = 3), " TW/mm²/sr")
eb = energy_balance(sol)
println("X-ray power (TW): source ", round(eb.source; digits = 1),
        ", wall ", round(eb.absorbed[:wall]; digits = 1),
        ", capsule ", round(eb.absorbed[:capsule]; digits = 1),
        ", LEHs ", round(eb.escaped[:leh_top] + eb.escaped[:leh_bottom]; digits = 1))
