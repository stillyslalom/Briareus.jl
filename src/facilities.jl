# Facility beam configurations.

"""
    beams_from_table(rows; spot, fnumber, power) -> Vector{Beam}

Build beams from rows (NamedTuples or Dicts) with keys `name`, `theta`, `phi`
(degrees) and optionally `group`, `quad`, `power`, `aim`, `defocus`, `fnumber`, `spot`.
The keyword arguments supply defaults for rows that lack them; each may be a value or a
function `row -> value`.
"""
function beams_from_table(rows; spot = SuperGaussianSpot(0.5), fnumber = 20.0, power = 1.0)
    _get(r, k, default) = haskey(r, k) ? r[k] : (default isa Function ? default(r) : default)
    return [Beam(string(r[:name]); theta = r[:theta], phi = r[:phi],
                 group = _get(r, :group, :none), quad = string(_get(r, :quad, "")),
                 power = _get(r, :power, power), aim = _get(r, :aim, (0.0, 0.0, 0.0)),
                 defocus = _get(r, :defocus, 0.0), fnumber = _get(r, :fnumber, fnumber),
                 spot = _get(r, :spot, spot))
            for r in rows]
end

# ---------------------------------------------------------------------------------------
# NIF
# ---------------------------------------------------------------------------------------

# Indirect-drive quad ports (port number, θ, φ in degrees), NIF User Guide §4.1,
# Table 4-1. Port 63 is printed as "62" there. Lower ports are (180-θ, 360-φ) of upper.
const NIF_PORTS = (
    (1, 23.5, 78.75), (2, 23.5, 168.75), (3, 23.5, 258.75), (4, 23.5, 348.75),
    (5, 30.58, 34.33), (6, 30.0, 123.75), (7, 30.58, 214.33), (8, 30.0, 303.75),
    (9, 44.5, 16.29), (10, 44.5, 62.46), (11, 44.5, 106.29), (12, 44.5, 152.46),
    (13, 44.5, 196.29), (14, 44.5, 242.46), (15, 44.5, 286.29), (16, 44.5, 332.46),
    (17, 50.0, 39.67), (18, 50.0, 84.38), (19, 50.0, 129.38), (20, 50.0, 174.38),
    (21, 50.0, 219.67), (22, 50.0, 264.38), (23, 50.0, 309.38), (24, 50.0, 354.38),
    (49, 130.0, 5.62), (50, 130.0, 50.62), (51, 130.0, 95.62), (52, 130.0, 140.33),
    (53, 130.0, 185.62), (54, 130.0, 230.62), (55, 130.0, 275.62), (56, 130.0, 320.33),
    (57, 135.5, 27.54), (58, 135.5, 73.71), (59, 135.5, 117.54), (60, 135.5, 163.71),
    (61, 135.5, 207.54), (62, 135.5, 253.71), (63, 135.5, 297.54), (64, 135.5, 343.71),
    (65, 150.0, 56.25), (66, 149.42, 145.67), (67, 150.0, 236.25), (68, 149.42, 325.67),
    (69, 156.5, 11.25), (70, 156.5, 101.25), (71, 156.5, 191.25), (72, 156.5, 281.25),
)

function _nif_cone(θ)
    c = min(θ, 180 - θ)
    c < 27 ? :inner23 : c < 37 ? :inner30 : c < 47 ? :outer44 : :outer50
end

# Best-focus CPP spots (1/e semi-axes in mm). 30°, 44.5°, 50°: super-Gaussian fits of
# Kemp et al., Phys. Plasmas 23, 101204 (2016), Table III. 23.5°: NIF User Guide
# Table 4-3 FWHM (1764 × 1262 µm) converted at the 30° fit's order (n = 8.7). Minor axis
# along the port's θ̂ (in the plane of the beam and the chamber axis), per Kemp.
const NIF_SPOTS = Dict(
    :inner23 => (major = 0.920, minor = 0.658, order = 8.7),
    :inner30 => (major = 0.897, minor = 0.6206, order = 8.7),
    :outer44 => (major = 0.6898, minor = 0.3936, order = 4.3),
    :outer50 => (major = 0.6744, minor = 0.3613, order = 3.9),
)

"""
    nif_spot(cone; scale = 1.0)

Best-focus CPP spot for a NIF cone (`:inner23`, `:inner30`, `:outer44`, `:outer50`),
major axis along the port's φ̂ and minor along θ̂.
"""
function nif_spot(cone::Symbol; scale = 1.0)
    s = NIF_SPOTS[cone]
    return SuperGaussianSpot(scale * s.major, scale * s.minor; order = s.order, angle = π / 2)
end

"""
    nif_quads(; power = 8.0, fnumber = 7.9, spot = nif_spot)

The 48 indirect-drive NIF quads as single beams (`power` in TW per quad, i.e. four
beamlines). Port angles are the NIF User Guide table; `group` is the cone (`:inner23`,
`:inner30`, `:outer44`, `:outer50`, shared by top and bottom) and `quad` is
`"P<port>"`. `spot` maps a cone to its `SuperGaussianSpot`.
"""
function nif_quads(; power = 8.0, fnumber = 7.9, spot = nif_spot)
    rows = [(name = "P$p", theta = θ, phi = φ, group = _nif_cone(θ), quad = "P$p",
             spot = spot(_nif_cone(θ))) for (p, θ, φ) in NIF_PORTS]
    return beams_from_table(rows; fnumber, power)
end

"""
    nif_beams(; power = 2.0, fnumber = 20.7, spot = nif_spot, dtheta = 2.35, darc = 2.07)

All 192 indirect-drive NIF beams, four per quad port. Within a quad the beams sit at
`±dtheta` degrees in polar angle and `±darc` degrees of arc in azimuth about the port
axis — derived from the published 63.2 cm (polar) × 55.7 cm (azimuthal) beam spacing and
7.7 m focal length; per-beam tables are not public. Names are `"P<port><a|b|c|d>"`.
"""
function nif_beams(; power = 2.0, fnumber = 20.7, spot = nif_spot, dtheta = 2.35, darc = 2.07)
    rows = NamedTuple[]
    for (p, θ, φ) in NIF_PORTS
        cone = _nif_cone(θ)
        for (tag, sθ, sφ) in (("a", -1, -1), ("b", -1, 1), ("c", 1, -1), ("d", 1, 1))
            θb = θ + sθ * dtheta
            φb = φ + sφ * darc / sind(θ)
            push!(rows, (name = "P$p$tag", theta = θb, phi = φb, group = cone, quad = "P$p",
                         spot = spot(cone)))
        end
    end
    return beams_from_table(rows; fnumber, power)
end

# ---------------------------------------------------------------------------------------
# OMEGA
# ---------------------------------------------------------------------------------------

"""
    omega_ports() -> Vector{NamedTuple}

The 60 OMEGA beam directions in a hohlraum frame whose z axis is a pentagon (P6–P7)
axis. They are the vertices of a truncated icosahedron with hexagon:pentagon edge
ratio 1.2 (cut at 1.2/3.4 of each icosahedron edge), which gives cones at 21.415°,
42.020°, 58.852° and 81.249° and their supplements (Sci. Rep. 13, 9883, 2023).
The azimuthal origin puts cone 1 at φ = 0°, 72°, …; it is not tied to chamber
coordinates. Groups are `:cone1`–`:cone4`; lower-hemisphere beams use the same groups.
"""
function omega_ports()
    ϕ = (1 + sqrt(5)) / 2
    V = [Vec3(0, 1, ϕ), Vec3(0, -1, ϕ), Vec3(0, 1, -ϕ), Vec3(0, -1, -ϕ), Vec3(1, ϕ, 0),
         Vec3(-1, ϕ, 0), Vec3(1, -ϕ, 0), Vec3(-1, -ϕ, 0), Vec3(ϕ, 0, 1), Vec3(-ϕ, 0, 1),
         Vec3(ϕ, 0, -1), Vec3(-ϕ, 0, -1)]
    V = normalize.(V)
    R = rotation_to(V[1])'                           # takes V[1] to +z
    V = [R * v for v in V]
    L = minimum(norm(V[i] - V[j]) for i in 1:12 for j in 1:12 if i != j)
    t = 1.2 / 3.4
    pts = [normalize(V[i] + t * (V[j] - V[i])) for i in 1:12 for j in 1:12
           if i != j && norm(V[i] - V[j]) < 1.01L]
    θs = [acosd(clamp(p[3], -1, 1)) for p in pts]
    φs = [mod(atand(p[2], p[1]), 360) for p in pts]
    # rotate azimuths so that cone 1 sits at φ = 0
    φ1 = minimum(φs[i] for i in eachindex(pts) if θs[i] < 30)
    φs = mod.(φs .- φ1, 360)
    cones = (21.415, 42.020, 58.852, 81.249)
    cone(θ) = Symbol(:cone, argmin(abs.(min(θ, 180 - θ) .- cones)))
    order = sortperm(collect(zip(round.(θs; digits = 3), round.(φs; digits = 3))))
    return [(name = @sprintf("B%02d", n), theta = θs[i], phi = φs[i], group = cone(θs[i]))
            for (n, i) in enumerate(order)]
end

"""
    omega_spot(name = :EIDI300)

OMEGA phase-plate spots (1/e semi-axes, mm): `:EIDI300` hohlraum elliptical plate
(n = 4.1, 1/e diameters 290 × 212 µm), `:SG4` (r₀ = 0.36 mm, n = 4.2 with SSD+PS),
`:SG8` (r₀ = 0.412 mm, n = 4.7). OMEGA Users' Guide (2014).
"""
function omega_spot(name::Symbol = :EIDI300)
    name === :EIDI300 && return SuperGaussianSpot(0.145, 0.106; order = 4.1, angle = π / 2)
    name === :SG4 && return SuperGaussianSpot(0.36; order = 4.2)
    name === :SG8 && return SuperGaussianSpot(0.412; order = 4.7)
    throw(ArgumentError("unknown OMEGA spot $name"))
end

"""
    omega_beams(; cones = (:cone1, :cone2, :cone3), power = 0.5, fnumber = 6.7,
                spot = omega_spot(:EIDI300))

OMEGA beams in the hohlraum (P6–P7) frame (see `omega_ports`), restricted to `cones`
(the 40 indirect-drive beams by default; pass `cones = :all` for all 60).
"""
function omega_beams(; cones = (:cone1, :cone2, :cone3), power = 0.5, fnumber = 6.7,
                     spot = omega_spot(:EIDI300))
    rows = [r for r in omega_ports() if cones === :all || r.group in cones]
    return beams_from_table(rows; spot, fnumber, power)
end
