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

# Best-focus CPP spots (1/e semi-axes in mm) from the measured profiles of LMF's
# phase-plate-table-updated-10-15-2020.xlsx (sheet "Measured profiles"): super-Gaussian
# fits `exp(-½ρⁿ)` to beam-normal TAS images of rod shot N151022-002-006, averaged per
# cone, whose tabulated radii are half widths at half maximum. Converted here to the 1/e
# semi-axes of `exp(-ρⁿ)`: a = r_HWHM / (ln 2)^(1/n). CPPs: Rev3 Inner (23.5°),
# 1.0 Rev1a Inner (30°), Rev3 Outer (44.5°), Rev1 Outer (50°). Major axis horizontal
# (along the port's φ̂), minor axis in the polar plane.
const NIF_SPOTS = Dict(
    :inner23 => (major = 0.88585 / log(2)^(1 / 7.37167), minor = 0.65579 / log(2)^(1 / 7.37167),
                 order = 7.37167),
    :inner30 => (major = 0.85838 / log(2)^(1 / 8.046), minor = 0.60591 / log(2)^(1 / 8.046),
                 order = 8.046),
    :outer44 => (major = 0.63727 / log(2)^(1 / 4.7925), minor = 0.37234 / log(2)^(1 / 4.7925),
                 order = 4.7925),
    :outer50 => (major = 0.61910 / log(2)^(1 / 4.51833), minor = 0.36187 / log(2)^(1 / 4.51833),
                 order = 4.51833),
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

# Quad names of the NIF ports, read from the LMF chamber-geometry figures 3 (top quads)
# and 4 (bottom quads), both top views with φ counterclockwise from Switchyard 1
# (https://lmf.llnl.gov/for-users/experimental-capabilities/chamber-geometry, LLNL
# 2012-040636 / 2012-040642). Each label was located by its projected radius (∝ sin θ,
# separating the cones) and azimuth and assigned to the nearest port of that cone; all
# 48 assignments are distinct and within 9° of the tabulated azimuth. The figure's
# "FABS (130-186)_36B" label fixes the bottom view's orientation (Q36B = port 53).
"""
    NIF_QUAD_NAMES

NIF quad name of each indirect-drive chamber port (`Dict(port => "Q11T", …)`), read from
the labeled LMF chamber-geometry figures; see `nif_port`, `nif_bundle`.
"""
const NIF_QUAD_NAMES = Dict(
    1 => "Q15T", 2 => "Q42T", 3 => "Q33T", 4 => "Q24T",
    5 => "Q13T", 6 => "Q41T", 7 => "Q36T", 8 => "Q26T",
    9 => "Q14T", 10 => "Q12T", 11 => "Q43T", 12 => "Q45T",
    13 => "Q35T", 14 => "Q32T", 15 => "Q25T", 16 => "Q22T",
    17 => "Q11T", 18 => "Q16T", 19 => "Q46T", 20 => "Q44T",
    21 => "Q34T", 22 => "Q31T", 23 => "Q23T", 24 => "Q21T",
    49 => "Q13B", 50 => "Q11B", 51 => "Q41B", 52 => "Q46B",
    53 => "Q36B", 54 => "Q34B", 55 => "Q26B", 56 => "Q23B",
    57 => "Q14B", 58 => "Q12B", 59 => "Q43B", 60 => "Q45B",
    61 => "Q35B", 62 => "Q32B", 63 => "Q25B", 64 => "Q22B",
    65 => "Q16B", 66 => "Q44B", 67 => "Q31B", 68 => "Q21B",
    69 => "Q15B", 70 => "Q42B", 71 => "Q33B", 72 => "Q24B",
)

"""
    nif_port(quad) -> Int

Chamber port number of a NIF quad name such as `"Q31B"` (see `NIF_QUAD_NAMES`).
"""
function nif_port(quad::AbstractString)
    for (p, q) in NIF_QUAD_NAMES
        q == quad && return p
    end
    throw(KeyError(quad))
end

"""
    nif_bundle(beam_or_quad) -> String

Bundle of a NIF quad or beam: the cluster and bundle digits of its quad name, e.g.
`"31"` for Q31T and Q31B. A bundle's two quads sit in opposite hemispheres and need
not share a cone.
"""
nif_bundle(quad::AbstractString) = quad[2:3]
nif_bundle(b::Beam) = nif_bundle(b.quad)

"""
    nif_bundles(beams) -> Dict{String,Vector{Int}}

Indices of `beams` (NIF quads or beams) belonging to each bundle.
"""
function nif_bundles(beams::AbstractVector{<:Beam})
    out = Dict{String,Vector{Int}}()
    for (i, b) in enumerate(beams)
        push!(get!(out, nif_bundle(b), Int[]), i)
    end
    return out
end

"""
    nif_quads(; power = 8.0, fnumber = 7.9, spot = nif_spot)

The 48 indirect-drive NIF quads as single beams (`power` in TW per quad, i.e. four
beamlines). Port angles are the NIF User Guide table; `group` is the cone (`:inner23`,
`:inner30`, `:outer44`, `:outer50`, shared by top and bottom), and `name` and `quad` are
the quad name (e.g. `"Q31B"`, see `NIF_QUAD_NAMES`; `nif_port`, `nif_bundle`). `spot`
maps a cone to its `SuperGaussianSpot`.
"""
function nif_quads(; power = 8.0, fnumber = 7.9, spot = nif_spot)
    rows = [(name = NIF_QUAD_NAMES[p], theta = θ, phi = φ, group = _nif_cone(θ),
             quad = NIF_QUAD_NAMES[p],
             spot = spot(_nif_cone(θ))) for (p, θ, φ) in NIF_PORTS]
    return beams_from_table(rows; fnumber, power)
end

"""
    nif_beams(; power = 2.0, fnumber = 20.7, spot = nif_spot)

All 192 indirect-drive NIF beams with their individual directions from the LMF beam
table (`data/nif_beams.csv`, transcribed from beam-port-locations-updated-10-15-2020.xlsx
of the LMF Shot RI Resources). Names are the facility beam names (`"B111"`–`"B468"`:
bundle `xy`, beam `z`; beams 1–4 form quad QxyT, 5–8 QxyB); `quad` is the quad name and
`group` the cone. The mean direction of each quad's four beams is its port direction to
within 0.004°.
"""
function nif_beams(; power = 2.0, fnumber = 20.7, spot = nif_spot)
    rows = NamedTuple[]
    for r in _nif_beam_table()
        cone = _nif_cone(r.theta)
        push!(rows, (name = r.beam, theta = r.theta, phi = r.phi, group = cone,
                     quad = "Q" * r.bundle * r.hemisphere, spot = spot(cone)))
    end
    return beams_from_table(rows; fnumber, power)
end

function _nif_beam_table()
    _, rows = read_csv(joinpath(@__DIR__, "..", "data", "nif_beams.csv"))
    return [(beam = r[1], bundle = r[2], hemisphere = r[3], theta = parse(Float64, r[5]),
             phi = parse(Float64, r[6])) for r in rows]
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
