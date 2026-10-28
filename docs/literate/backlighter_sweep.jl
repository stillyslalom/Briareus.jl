# # Backlighter placement sweep
#
# An x-ray radiograph needs a backlighter: a foil heated by a few laser beams, placed
# behind the object on the diagnostic's line of sight, whose x rays cross the object to
# the detector. Where the foil goes is a trade. Close to the object, the foil subtends a
# wide cone, and the x-ray bundle from it can clip hardware near the line of sight;
# the beams that heat it also pass closer to that hardware. Far from it, the foil's
# emission spreads over a larger solid angle before it reaches the object.
#
# This tutorial designs the pointing of an OMEGA-style backlighter in the foil's own
# coordinates, then moves the foil along the line of sight with its beams following
# and checks every step for beam clearance, an unobstructed line of sight and the power
# the foil captures. It introduces frames ([`Frame`](@ref), [`foil_frame`](@ref)),
# pointing parameterized in a frame ([`FramePointing`](@ref)) and placement sweeps
# ([`placement_sweep`](@ref)).
#
# The target is a mock-up with generic dimensions, not a fielded design: a 2 mm square
# gold foil, a 0.43 mm radius plastic sphere at target-chamber center (TCC) as the
# object, and between them an aperture plate, a gold washer whose bore limits what the
# detector sees of the foil's surroundings.

using Briareus
using CairoMakie
using LinearAlgebra
using Printf
CairoMakie.activate!(type = "png", px_per_unit = 2)

# ## The foil and its frame
#
# The line of sight runs from TCC toward the detector at polar angle 90° and azimuth 0°,
# in the frame of [`omega_ports`](@ref). [`backlighter_foil`](@ref) places the foil on
# the far side, a standoff ``D`` from TCC, facing the detector. The same placement as a
# frame, [`foil_frame`](@ref), has its origin at the foil's center, its ``x`` and ``y``
# axes along the foil's edges and its ``z`` axis along the front normal, toward the
# detector. Local coordinates ``(u, v, 0)`` are points on the foil's face.

los = (90.0, 0.0)
D0 = 8.0                       # standoff of the design (mm)
foil_edge = 2.0
F = foil_frame(los; standoff = D0)
foil_part = backlighter_foil(los; standoff = D0, width = foil_edge, n = 8)
F

# The aperture plate sits 2.5 mm from TCC toward the foil, normal to the line of sight,
# with a 1.0 mm bore radius and a 1.5 mm outer radius. Its frame is the line of sight's
# [`port_frame`](@ref) moved to that point, so the plate's mesh, built about the
# ``z`` axis, lands across the line of sight.

d = port_vector(los...)
X0, bore, outer = 2.5, 1.0, 1.2
P = port_frame(los...; origin = -X0 .* d)
washer = Component(:washer, transform(annulus_mesh(bore, outer, 0.0; nr = 2, nphi = 64), P.R, P.origin),
                   gold(); two_sided = true)
capsule = Component(:capsule, sphere_mesh(0.43; ntheta = 12, nphi = 24), ablator())
target = Target([foil_part, washer, capsule])

# ## Beams pointed in the foil's frame
#
# The beams come from the ports on the detector's side, 20° to 40° from the line of
# sight: three symmetric pairs, at 25.0°, 33.2° and 38.5°. Each delivers 0.4 TW in an
# SG4 phase-plate spot ([`omega_spot`](@ref)); the powers are illustrative. They start
# aimed at a ring of radius 0.3 mm on the foil, written in the foil's frame and mapped
# to the chamber by [`to_chamber`](@ref).

beams = [b for b in omega_beams(cones = :all, spot = omega_spot(:SG4), power = 0.4)
         if 20 < acosd(dot(b.port, d)) < 40]
ring(i) = (0.3 * cospi(2i / length(beams)), 0.3 * sinpi(2i / length(beams)), 0.0)
beams = [remake(b; aim = to_chamber(F, ring(i))) for (i, b) in enumerate(beams)]
[(b.name, round(acosd(dot(b.port, d)); digits = 1)) for b in beams]

# A [`FramePointing`](@ref) makes each beam's aim a pair of parameters ``(u, v)`` in
# the foil's frame. Every beam is its own group (`by = :beam`), so there are twelve
# parameters, and each maps to the chamber as ``R\,(u, v, 0) + o``: the parameters keep
# their meaning wherever the foil goes.

pointing = FramePointing(beams, F; by = :beam, fields = (:u, :v))
scenario = Scenario(target, pointing)

# The objective is the uniformity of the absorbed flux over the foil's front, the
# squared rms nonuniformity, plus the squared fraction of the beams' power that misses
# it. It refers to the foil's elements only, not to chamber positions, so it serves at
# any placement.

function foil_objective(ev)
    t = ev.scenario.target
    front = [e for e in elements(t, :backlighter) if t.elem_side[e] == 1]
    absorbed = ev.deposition.element[front]
    area = t.elem_area[front]
    return rms_nonuniformity(absorbed ./ area, area)^2 +
           (1 - sum(absorbed) / ev.deposition.delivered)^2
end

result = optimize_pointing(scenario, foil_objective; iterations = 100)

# The optimizer spreads the spots across the face, each within the default ±1 mm box
# about its start:

for (name, u, v) in zip(beams, result.x[1:2:end], result.x[2:2:end])
    @printf("%s  u = %+.3f mm  v = %+.3f mm\n", name.name, u, v)
end
designed = FramePointing(result.beams, F; by = :beam, fields = (:u, :v))

# ## Sweeping the standoff
#
# [`placement_sweep`](@ref) moves the foil along its frame's ``z`` axis, toward TCC for
# a positive shift, with every beam following (`follow = true`): each step translates
# the foil and the beams' aims by the same vector, rebuilds the target and runs the
# checks. Here they are the beam clearance of [`clearance`](@ref) to every part but
# the foil (with its quick settings and without the NIF clearance wings, which are
# engineering defaults for NIF), the [`line_of_sight`](@ref) check through the
# sphere's disk, and the fraction of each beam's power the foil absorbs. Standoffs run
# from 4 to 9 mm.

standoffs = 4.0:0.5:9.0
checks = (; clearance_kw = (; targets = (:backlighter,), wings = (), nangles = 8, nspot = 31),
          los, los_kw = (; object = :capsule, radius = 0.43),
          capture = :backlighter, capture_kw = (; resolution = 0.6))
sweep = placement_sweep(target, designed, :backlighter, D0 .- standoffs; frame = F, axis = :z,
                        follow = true, checks...)

# The rows list, for each shift ``D_0 - D``, the least worst-case clearance between a
# beam's envelope and the hardware (negative: the envelope overlaps it under the error
# set), the beam and part it belongs to, the failing and review pairs, the line of
# sight, and the least captured fraction:

sweep

# Two limits appear. At 4 mm the line of sight is blocked: the x-ray bundle from the
# foil's corners widens past the aperture plate's bore. The beams set the tighter
# limit. A beam's axis crosses the plate's plane about ``(D - X_0) \tan\alpha`` from
# the line of sight, ``\alpha`` the angle of its port from it, shifted by its aim on the
# foil; a 25° beam aimed 0.68 mm off the foil's center (`B31`) passes closest, and its worst-case
# envelope clears the plate only beyond about 8.8 mm. At the design standoff of 8 mm it
# still overlaps under the error set (a review pair: some light of its displaced
# envelope reaches the plate). The captured fraction stays the same from 7 to 9 mm,
# where the nominal beams miss the plate, since the foil and its beams move as one;
# closer in, the plate takes part of the light.

fig = Figure(size = (800, 360))
ax1 = Axis(fig[1, 1], xlabel = "standoff D (mm)", ylabel = "least clearance (mm)",
           title = "Beams and the aperture plate")
clear = [r.min_clearance for r in sweep.rows]
hlines!(ax1, [0.0], color = :gray60, linestyle = :dash)
scatterlines!(ax1, collect(standoffs), clear, color = :black)
ax2 = Axis(fig[1, 2], xlabel = "standoff D (mm)", ylabel = "fraction",
           title = "Line of sight and capture")
scatterlines!(ax2, collect(standoffs), [r.blocked_fraction for r in sweep.rows],
              label = "rays blocked")
scatterlines!(ax2, collect(standoffs), [r.min_captured for r in sweep.rows],
              label = "least captured")
axislegend(ax2, position = :rc)
fig

# The rows can be written for a spreadsheet with [`write_sweep_csv`](@ref).
# `placement_sweep` also re-optimizes the pointing at each step when given an
# `optimize` named tuple, starting from the previous step's optimum; with the beams
# following the foil rigidly, as here, the design does not change, so that is left off.

# ## The line-of-sight threshold
#
# The sweep's steps bracket the standoff where the line of sight opens; its exact
# value follows from similar triangles. The widest ray of the bundle runs from a corner
# of the foil, at ``w\sqrt{2}`` from the line of sight (``w`` the half-width), through
# the opposite rim of the object's disk of radius ``a``. At the aperture plate,
# a fraction ``s = X_0/D`` of the way from the disk to the foil, it lies
# ``(1 - s)\,a + s\,w\sqrt{2}`` from the axis, and it reaches the bore's rim ``\rho``
# at
#
# ```math
# D^* = X_0\,\frac{w\sqrt{2} - a}{\rho - a}.
# ```

w, a = foil_edge / 2, 0.43
Dstar = X0 * (w * sqrt(2) - a) / (bore - a)
edge = placement_sweep(target, designed, :backlighter, D0 .- [Dstar + 1e-6, Dstar - 1e-6];
                       frame = F, axis = :z, follow = true, los, los_kw = checks.los_kw)
@printf("D* = %.6f mm: %s just outside, %s just inside\n", Dstar,
        edge.rows[1].los_clear ? "clear" : "blocked", edge.rows[2].los_clear ? "clear" : "blocked")

# The followed aims move by exactly the shift along the foil's normal:

shifts = [maximum(norm(b.aim - (c.aim + r.value * F.R[:, 3])) for (b, c) in zip(r.beams, result.beams))
          for r in sweep.rows]
@printf("largest aim error over the sweep: %.1e mm\n", maximum(shifts))

# ## What this checks
#
# - The line of sight opens at the similar-triangles standoff ``D^*``: the check
#   through the sphere's disk is clear a nanometer outside it and blocked a nanometer
#   inside, because the bundle's extreme ray (a foil corner through the disk's rim, both
#   sampled) meets a vertex of the bore's polygon.
# - The followed beams' aims equal the optimized aims plus the shift along the foil's
#   normal at every step, exactly.
# - The captured fraction is the same, to the five digits shown, at every step where
#   the nominal beams miss the aperture plate, since the foil and its beams translate
#   together.
