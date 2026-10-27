# Capsule drive symmetry in a NIF hohlraum: evaluate a published pointing, then tune
# cone pointing and power balance to null P2 and P4.
#
#   julia --project=. -t auto examples/nif_hohlraum.jl

using Briareus

# Build-A-Hohlraum geometry (5.75 mm × 10.1 mm Au, 3.1 mm LEHs), capsule at TCC.
target = hohlraum(length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                  wall = gold(albedo = 0.85, conversion = 0.8))
display(target)

# 48 quads; cone pointing in hohlraum (r, z) coordinates (mm). Positive r points toward
# the quad's own port azimuth, so the published "back-azimuth" offsets are negative.
quads = nif_quads(power = 8.0)
pointing = GroupPointing(quads; fields = (:z, :r, :power))
x0 = parameters(pointing; inner23 = (z = 3.09, r = -1.09), inner30 = (z = 3.63, r = -0.60),
                outer44 = (z = 5.20, r = 0.0), outer50 = (z = 5.01, r = 0.0))

scenario = Scenario(target, pointing; observe = :capsule)   # closed: the hohlraum default
println("\nPublished pointing:")
report(evaluate(scenario, x0))

# Null P2, P4 (and keep P6 small) while holding the drive at its initial level, without
# clipping the LEH or hitting the capsule.
a0 = flux_moments(evaluate(scenario, x0))[1]
objective = symmetry_objective(modes = (2 => 1.0, 4 => 1.0, 6 => 0.25), drive = a0)
lower, upper = bounds(pointing, x0; z = 1.5, r = 1.0, power = 0.4)
result = optimize_pointing(scenario, objective; x0, lower, upper)

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
