# Area backlighter: point OMEGA beams onto a foil for a uniform irradiance spot without
# clipping the main target, then report the x-ray power toward the imager.
#
#   julia --project=. -t auto examples/backlighter.jl

using Briareus
using LinearAlgebra

# Foil 6 mm from TCC along +x, driven on the face toward TCC; a 1 mm-radius sphere at
# TCC stands in for the main target. The imager looks along +x through TCC at the foil.
foil_center = Vec3(6.0, 0.0, 0.0)
backlighter = foil(:backlighter, foil_center, (-1, 0, 0); width = 2.0, n = 40,
                   material = Material(:Ti; conversion = 0.02), two_sided = false)
main = Component(:main, sphere_mesh(1.0; ntheta = 16, nphi = 24), BLACK)
target = Target(backlighter, main)

# OMEGA beams that can reach the driven face (they propagate toward +x) from ports at
# least 40° away from the imager line of sight, each initially aimed at the foil center.
beams = [remake(b; aim = foil_center, power = 0.4)
         for b in omega_beams(cones = :all, spot = omega_spot(:SG4))
         if b.port[1] < -cosd(70)]
println(length(beams), " beams drive the backlighter")

# Each beam's aim point is free in 3D within ±0.6 mm.
pointing = GroupPointing(beams; by = :beam, fields = (:r, :phi, :z), mirror = false)
scenario = Scenario(target, beams; parameterization = pointing)
x0 = initial_parameters(pointing)
objective = irradiance_objective(:backlighter; center = foil_center, radius = 0.5, spill = 0.2)
lower, upper = bounds(pointing, x0; r = 0.6, phi = 0.6, z = 0.6)
result = optimize_pointing(scenario, objective; x0, lower, upper)

spot = [e for e in elements(target, :backlighter)
        if norm(target.elem_centroid[e] - foil_center) < 0.5]
A = target.elem_area[spot]
for (label, x) in (("initial", x0), ("optimized", result.x))
    ev = evaluate(scenario, x)
    I = laser_flux(ev.deposition, target)[spot]
    println(rpad(label, 10), ": spot irradiance ",
            round(to_Wcm2(sum(I .* A) / sum(A)); sigdigits = 3), " W/cm², rms nonuniformity ",
            round(100 * rms_nonuniformity(I, A); digits = 1), "%, on main target ",
            round(100 * power_fraction(ev, :main); digits = 2), "%")
end

# X-ray emission (conversion × absorbed laser) toward the imager, which sits at -x and
# sees the driven face through TCC — the main target blocks part of it.
ev = evaluate(scenario, result.x)
model = RadiationModel(target)
sol = solve(model, ev.deposition)
println("X-ray power toward the imager: ",
        round(line_of_sight_power(target, sol.B, (-1, 0, 0)); sigdigits = 3), " TW/sr")
