# This file handles interchange with VisRad (Prism Computational Sciences, VISRAD 21.1.1).
# It reads and writes pointing, beam-summary and power files, reads element tables, and
# writes a benchmark kit for rebuilding a Briareus case in VisRad.
#
# Every format detail cites a page of the online manual, https://prism-cs.com/Manuals/VisRad/
# (abbreviated `VR/` below). Details that the manual does not fix are labeled
# *assumption*. The sources and confidence of each detail are collected in
# `reference/research-notes.md`, section "VisRad interchange formats".
#
# Some conventions are shared with VisRad (VR/power_sources/set_laser_beam_pointing.html).
# A beam is pointed at a Reference Point, and its focus lies a signed "focus offset" from
# it along the propagation direction, with a positive offset moving the focus along the
# propagation. These are exactly the `aim` and `defocus` of Briareus. VisRad draws the
# beam axis from the port position on the chamber wall through the Reference Point, and
# Briareus keeps the axis parallel to the port direction. The two differ by
# |aim⊥|/R_chamber (≲ 1e-3 rad for mm-scale aims), that is, by ≲ 1 µm at the focus for
# defocus ≲ 1 mm.

const _VR_COMMENT = "VisRad manual (VISRAD 21.1.1): https://prism-cs.com/Manuals/VisRad/"

# ---------------------------------------------------------------------------------------
# Low-level text helpers
# ---------------------------------------------------------------------------------------

_vr_v(x) = Float64(ForwardDiff.value(x))
_vr_fmt(x::Real) = @sprintf("%.10g", _vr_v(x))
_vr_field(x::Real) = _vr_fmt(x)
_vr_field(x) = _csvfield(string(x))

function _vr_lines(path::AbstractString)
    s = read(path, String)
    startswith(s, "\ufeff") && (s = s[nextind(s, 1):end])
    return [rstrip(l, '\r') for l in split(s, '\n')]
end

_vr_delim(line) = occursin('\t', line) ? '\t' :
                  count(==(';'), line) > count(==(','), line) ? ';' : ','

# Splits one delimited line, honoring double quotes and "" escapes (as `read_csv` does).
function _vr_split(line::AbstractString, delim::Char)
    out = String[]
    buf = IOBuffer()
    q = false
    cs = collect(line)
    i = 1
    while i <= length(cs)
        c = cs[i]
        if q
            if c == '"'
                (i < length(cs) && cs[i+1] == '"') ? (write(buf, '"'); i += 1) : (q = false)
            else
                write(buf, c)
            end
        elseif c == '"'
            q = true
        elseif c == delim
            push!(out, strip(String(take!(buf))))
        else
            write(buf, c)
        end
        i += 1
    end
    push!(out, strip(String(take!(buf))))
    return out
end

function _vr_num(s::AbstractString)
    t = strip(s)
    (isempty(t) || t == "-") && return nothing
    return tryparse(Float64, t)
end

# Returns the length scale to mm from the unit in a header, e.g. "X (cm)", or `default`
# when the header gives no unit.
function _vr_length_scale(header::AbstractString, default::Real)
    m = match(r"\(\s*(mm|cm|um|µm|μm|microns?|m)\s*\)"i, header)
    m === nothing && return Float64(default)
    u = lowercase(m[1])
    return u == "mm" ? 1.0 : u == "cm" ? 10.0 : u == "m" ? 1000.0 : 1e-3
end

# Returns the canonical beam key. NIF beams are "B115", and VisRad also accepts "BL115"
# (VR/files/laser_beam_data_files.html). OMEGA beams are "Beam_10" (same page), also
# written "BL-10" or as a bare beam number "10" in older summary tables with a "Beam #"
# column (VR/intro/images/laser_summary.jpg). Other names are compared
# case-insensitively.
function _vr_key(name::AbstractString)
    u = uppercase(strip(name))
    m = match(r"^BL-?(\d{3})$", u)
    m !== nothing && return "B" * m[1]
    m = match(r"^(?:BEAM|BL)?[ _-]*(\d{1,2})$", u)
    m !== nothing && return "Beam_" * m[1]
    m = match(r"^BEAM[ _-]*(\d{1,2})(_.*)?$", u)
    m !== nothing && return "Beam_" * m[1] * (m[2] === nothing ? "" : m[2])
    return u
end

_is_nif_quad(name::AbstractString) = occursin(r"^Q[1-4][1-8][TB]$", name)

# Beam names of each NIF quad, from the LMF beam table behind `nif_beams`.
function _nif_quad_members()
    d = Dict{String,Vector{String}}()
    for r in _nif_beam_table()
        push!(get!(d, "Q" * r.bundle * r.hemisphere, String[]), r.beam)
    end
    return d
end

# Returns the VisRad beams of a beam set. NIF quads (names "Q11T"…) become their four
# facility beams with the aim, defocus and spot of the quad and a quarter of its power,
# and other beams pass through unchanged. Returns `(beam, source_index, power_fraction)`
# per VisRad beam.
function _vr_expand(beams::AbstractVector{<:Beam})
    any(b -> _is_nif_quad(b.name), beams) || return [(b, i, 1.0) for (i, b) in enumerate(beams)]
    nb = Dict(b.name => b for b in nif_beams())
    members = _nif_quad_members()
    out = Tuple{Beam,Int,Float64}[]
    for (i, b) in enumerate(beams)
        if _is_nif_quad(b.name)
            for n in members[b.name]
                # the spot of the quad, oriented as on the quad, since the VisRad rotation is
                # quad-relative
                spot = _vr_turn(b.spot, _vr_reference_angle(nb[n]))
                push!(out, (remake(nb[n]; aim = b.aim, defocus = b.defocus, power = b.power / 4,
                                   spot, group = b.group), i, 0.25))
            end
        else
            push!(out, (b, i, 1.0))
        end
    end
    return out
end

# ---------------------------------------------------------------------------------------
# Spots and phase plates
# ---------------------------------------------------------------------------------------

"""
    VISRAD_PHASE_PLATES

The built-in phase-plate profiles of VisRad (`VR/power_sources/phase_plate_params.html`),
keyed by the names in the manual, as `(facility, major, minor, order, level)`. The radii
are in µm; the "(mm)" header in the manual is a µ rendered in the Symbol font. They are
measured where the intensity is `level` × peak: `exp(-1)` (1/e) for the OMEGA, OMEGA EP
and GEKKO tables, `0.5` (HWHM) for NIF CPPs and `0.03` for LMJ. Circular plates have
`minor == major`. `order` is the super-Gaussian exponent `n` of
`exp(-[(x/a)²+(y/b)²]^(n/2))` (`VR/power_sources/set_laser_profile_props.html`). The NIF
entries with n = 2.0 are placeholders in the manual where no data were available. The
GEKKO table does not state its radius definition, and 1/e is an *assumption*. See
`visrad_spot`.
"""
const VISRAD_PHASE_PLATES = Dict{String,NamedTuple{(:facility, :major, :minor, :order, :level),
                                                    Tuple{Symbol,Float64,Float64,Float64,Float64}}}(
    # OMEGA DPPs (1/e radii)
    "SG5-650" => (facility = :OMEGA, major = 287.0, minor = 287.0, order = 5.11, level = exp(-1)),
    "SG5" => (facility = :OMEGA, major = 357.0, minor = 357.0, order = 4.77, level = exp(-1)),
    "SG8" => (facility = :OMEGA, major = 438.0, minor = 438.0, order = 4.5, level = exp(-1)),
    "SG8-Flat" => (facility = :OMEGA, major = 438.0, minor = 438.0, order = 5.8, level = exp(-1)),
    "PD-SG2-600" => (facility = :OMEGA, major = 300.0, minor = 300.0, order = 2.2, level = exp(-1)),
    "E-IDI-300" => (facility = :OMEGA, major = 144.0, minor = 106.0, order = 4.3, level = exp(-1)),
    "E-SG4-865" => (facility = :OMEGA, major = 430.0, minor = 396.0, order = 4.7, level = exp(-1)),
    "E-SG10-300" => (facility = :OMEGA, major = 171.0, minor = 157.0, order = 2.5, level = exp(-1)),
    "100 um (s/n 083)" => (facility = :OMEGA, major = 100.0, minor = 100.0, order = 2.0, level = exp(-1)),
    "100 um (s/n 086)" => (facility = :OMEGA, major = 99.0, minor = 99.0, order = 2.2, level = exp(-1)),
    "100 um (SSD off, DPR out)" => (facility = :OMEGA, major = 50.0, minor = 50.0, order = 4.0, level = exp(-1)),
    "200 um (s/n 120)" => (facility = :OMEGA, major = 107.0, minor = 107.0, order = 2.0, level = exp(-1)),
    "200 um (s/n 018)" => (facility = :OMEGA, major = 105.0, minor = 105.0, order = 2.1, level = exp(-1)),
    "500 um (s/n 006)" => (facility = :OMEGA, major = 150.0, minor = 150.0, order = 1.9, level = exp(-1)),
    "500 um (s/n 026)" => (facility = :OMEGA, major = 149.0, minor = 149.0, order = 1.9, level = exp(-1)),
    "700 um (s/n 056)" => (facility = :OMEGA, major = 195.0, minor = 195.0, order = 1.9, level = exp(-1)),
    "700 um (s/n 058)" => (facility = :OMEGA, major = 202.0, minor = 202.0, order = 1.8, level = exp(-1)),
    "LLNL-3w-150-328" => (facility = :OMEGA, major = 75.0, minor = 75.0, order = 2.0, level = exp(-1)),
    "LLNL-3w-250-024" => (facility = :OMEGA, major = 125.0, minor = 125.0, order = 2.0, level = exp(-1)),
    "LLNL-3w-250-027" => (facility = :OMEGA, major = 125.0, minor = 125.0, order = 2.0, level = exp(-1)),
    # OMEGA EP (1/e radii)
    "SG8-0750" => (facility = :OMEGA_EP, major = 336.0, minor = 336.0, order = 8.92, level = exp(-1)),
    "SG8-0800" => (facility = :OMEGA_EP, major = 397.3, minor = 397.3, order = 11.86, level = exp(-1)),
    "SG8-1100" => (facility = :OMEGA_EP, major = 483.0, minor = 483.0, order = 8.57, level = exp(-1)),
    "SG8-1800" => (facility = :OMEGA_EP, major = 956.0, minor = 956.0, order = 24.94, level = exp(-1)),
    "SG8-2000" => (facility = :OMEGA_EP, major = 909.6, minor = 909.6, order = 15.68, level = exp(-1)),
    "SG10-0400" => (facility = :OMEGA_EP, major = 233.0, minor = 233.0, order = 16.66, level = exp(-1)),
    # NIF CPPs (HWHM radii)
    "0.7-Inner-2w" => (facility = :NIF, major = 577.0, minor = 413.0, order = 2.0, level = 0.5),
    "0.7-Outer-2w" => (facility = :NIF, major = 415.0, minor = 240.0, order = 2.0, level = 0.5),
    "1.0-Inner-1w" => (facility = :NIF, major = 824.0, minor = 590.0, order = 2.0, level = 0.5),
    "1.0-Inner-2w" => (facility = :NIF, major = 824.0, minor = 590.0, order = 2.0, level = 0.5),
    "1.07-Inner-1w" => (facility = :NIF, major = 882.0, minor = 631.0, order = 2.0, level = 0.5),
    "1.07-Inner-2w" => (facility = :NIF, major = 882.0, minor = 631.0, order = 2.0, level = 0.5),
    "1.0-Outer-2w" => (facility = :NIF, major = 593.0, minor = 343.0, order = 2.0, level = 0.5),
    "1.07-Outer-2w" => (facility = :NIF, major = 635.0, minor = 367.0, order = 2.0, level = 0.5),
    "1.0-Outer-1w" => (facility = :NIF, major = 593.0, minor = 343.0, order = 2.0, level = 0.5),
    "1.07-Outer-1w" => (facility = :NIF, major = 635.0, minor = 367.0, order = 2.0, level = 0.5),
    "1.175-Inner-2w" => (facility = :NIF, major = 968.0, minor = 693.0, order = 2.0, level = 0.5),
    "1.175-Outer-2w" => (facility = :NIF, major = 697.0, minor = 403.0, order = 2.0, level = 0.5),
    "Rev3-Inner-2w" => (facility = :NIF, major = 890.0, minor = 660.0, order = 7.6, level = 0.5),
    "1.0 Rev1a-Inner-2w" => (facility = :NIF, major = 860.0, minor = 610.0, order = 8.0, level = 0.5),
    "Rev3-Outer-2w" => (facility = :NIF, major = 640.0, minor = 370.0, order = 4.8, level = 0.5),
    "Rev1-Outer-2w" => (facility = :NIF, major = 620.0, minor = 360.0, order = 4.5, level = 0.5),
    "No CPP" => (facility = :NIF, major = 75.0, minor = 75.0, order = 2.0, level = 0.5),
    "400 um CPP" => (facility = :NIF, major = 189.0, minor = 189.0, order = 7.145, level = 0.5),
    "3 mm 2W CPP" => (facility = :NIF, major = 1575.0, minor = 1190.0, order = 13.1, level = 0.5),
    # LMJ (radii at 3% of peak)
    "Type A" => (facility = :LMJ, major = 757.0, minor = 328.0, order = 2.87, level = 0.03),
    "Type B" => (facility = :LMJ, major = 904.0, minor = 584.0, order = 3.83, level = 0.03),
    "Type D" => (facility = :LMJ, major = 440.0, minor = 419.0, order = 3.01, level = 0.03),
    "Type E" => (facility = :LMJ, major = 796.0, minor = 769.0, order = 4.28, level = 0.03),
    "Type F" => (facility = :LMJ, major = 315.0, minor = 315.0, order = 2.0, level = 0.03),
    # GEKKO (radius definition not stated; 1/e assumed)
    "DPP-1" => (facility = :GEKKO, major = 353.6, minor = 353.6, order = 6.0, level = exp(-1)),
)

# Phase-plate names as they appear in other places, in CMT CPP names
# (VR/power_sources/import_nif_cmt_data_file.html) and in summary-table model strings.
const _VR_PLATE_ALIASES = Dict("1.0 Rev1 Outer" => "Rev1-Outer-2w")

function _vr_plate_norm(s::AbstractString)
    t = lowercase(s)
    for w in ("omega", "nif", "lmj", "gekko", "dpp", "cpp", "cone", "phase plate", "profile")
        t = replace(t, w => "")
    end
    return replace(t, r"[^a-z0-9]" => "")
end

function _vr_plate(model::AbstractString)
    n = _vr_plate_norm(model)
    isempty(n) && return nothing
    for (k, v) in VISRAD_PHASE_PLATES
        _vr_plate_norm(k) == n && return v
    end
    for (k, v) in _VR_PLATE_ALIASES
        _vr_plate_norm(k) * "2w" == n && return VISRAD_PHASE_PLATES[v]
        _vr_plate_norm(k) == n && return VISRAD_PHASE_PLATES[v]
    end
    return nothing
end

# The rotation of the minor axis of a spot in the VisRad convention is in degrees from
# vertical, positive clockwise, and zero with the major axis "nearly horizontal"
# (VR/power_sources/set_laser_profile_props.html). The manual does not state from which
# side "clockwise" is seen. `sense = :beam` (default, *assumption*) views the spot along
# the propagation direction, for the beam "as it enters the target chamber". `:target`
# views it from TCC back toward the port. In the beam view, screen-up is -θ̂ and
# screen-right φ̂, so the minor axis is cos ρ θ̂ - sin ρ φ̂. From TCC, screen-right is -φ̂
# and the minor axis is cos ρ θ̂ + sin ρ φ̂.
function _vr_sense(sense::Symbol)
    sense === :beam && return -1.0
    sense === :target && return 1.0
    throw(ArgumentError("rotation_sense must be :beam or :target"))
end

# VisRad measures the rotation of a NIF beam from the horizontal of its quad. At rotation
# 0 the major axis lies along φ̂ of the quad port, which is turned about the beam axis by
# 1.6–5.3° from φ̂ of the beam. This was measured from VisRad 21 deposition on plates
# normal to 48 beams, where the fitted spot rotation equals this angle to 0.02° on all of
# them (reference/research-notes.md, "VisRad deposition parity"). Returns the angle of φ̂
# of the quad from φ̂ of the beam, in radians in the Briareus sense (from θ̂ toward φ̂).
# Returns 0 for beams of no NIF quad, such as a quad itself or OMEGA beams.
function _vr_reference_angle(b::Beam)
    q = get(_NIF_QUAD_PORTS, b.quad, nothing)
    q === nothing && return 0.0
    eθ, eφ = port_basis(b.port)
    _, qφ = port_basis(q)
    return atan(dot(qφ, eφ), dot(qφ, eθ)) - π / 2
end

const _NIF_QUAD_PORTS = Dict(NIF_QUAD_NAMES[p] => port_vector(θ, φ) for (p, θ, φ) in NIF_PORTS)

# Returns the NIF quad of a facility beam name, B<cluster><bundle><beam>, where beams 1–4
# form quad QxyT and beams 5–8 quad QxyB. Returns "" for other names.
function _nif_beam_quad(name::AbstractString)
    m = match(r"^B([1-4][1-6])([1-8])$", name)
    m === nothing && return ""
    return "Q" * m[1] * (parse(Int, m[2]) <= 4 ? "T" : "B")
end

# Completes a beam built from a file without a template. Its NIF quad comes from its
# name, and its spot, read relative to φ̂ of the beam, is turned to the horizontal of the
# quad.
function _vr_quad_framed(b::Beam)
    q = _nif_beam_quad(b.name)
    isempty(q) && return b
    b = Beam(b.name, b.port, b.aim, b.power, b.defocus, b.spot, b.fnumber, b.group, q, b.wavelength, b.dlambda)
    return remake(b; spot = _vr_turn(b.spot, _vr_reference_angle(b)))
end

_vr_turn(sp::SuperGaussianSpot, δ) = δ == 0 ? sp :
    SuperGaussianSpot(sp.a, sp.b, sp.order, sp.angle + δ, sp.halo, sp.halo_scale, sp.halo_order)

# Returns the VisRad rotation (degrees) of a Briareus spot, with `ref` the
# `_vr_reference_angle` of its beam.
function _vr_rotation(sp::SuperGaussianSpot, sense::Symbol, ref = 0.0)
    _vr_v(sp.a) == _vr_v(sp.b) && return 0.0
    ψ = rad2deg(sp.angle - ref)
    ψminor = _vr_v(sp.a) > _vr_v(sp.b) ? ψ + 90 : ψ
    ρ = _vr_sense(sense) * ψminor
    return mod(ρ + 90, 180) - 90
end

# Returns the Briareus spot from VisRad radii (µm) at intensity `level` × peak and
# rotation ρ (deg), with `ref` the `_vr_reference_angle` of the beam.
function _vr_make_spot(major, minor, ρ, order, level, sense::Symbol;
                       halo = 0.0, halo_scale = 3.0, halo_order = 2.0, ref = 0.0)
    s = log(1 / level)^(1 / order)
    a = major / 1000 / s
    b = minor / 1000 / s
    angle = a == b ? 0.0 : deg2rad(mod(90 + _vr_sense(sense) * ρ, 180)) + ref
    return SuperGaussianSpot(a, b; order, angle, halo, halo_scale, halo_order)
end

"""
    visrad_spot(plate; rotation = 0.0, rotation_sense = :beam) -> SuperGaussianSpot

Briareus spot equivalent to a VisRad built-in phase plate, given by its
`VISRAD_PHASE_PLATES` name or a CMT or summary-table spelling such as
`"Rev3 Inner Cone 2W CPP"` or `"OMEGA E-IDI-300 DPP"`. The tabulated radii, defined at
the intensity level of the table, are converted to 1/e semi-axes (mm) of the same
super-Gaussian, `r_1/e = r_f / (ln(1/f))^(1/n)`. `rotation` is the minor-axis rotation
of VisRad (degrees, clockwise from vertical, 0 with the major axis horizontal); see
`read_visrad_summary` for `rotation_sense`.
"""
function visrad_spot(plate::AbstractString; rotation = 0.0, rotation_sense::Symbol = :beam)
    p = _vr_plate(plate)
    p === nothing && throw(ArgumentError("unknown VisRad phase plate \"$plate\""))
    return _vr_make_spot(p.major, p.minor, rotation, p.order, p.level, rotation_sense)
end

"""
    VISRAD_NIF_PLATES

The VisRad phase plate (a `VISRAD_PHASE_PLATES` name) for each NIF cone: Rev3 Inner on
the 23.5° cone, 1.0 Rev1a Inner on the 30° cone, Rev3 Outer on the 44.5° cone and Rev1
Outer on the 50° cone. The assignment is that of the measured CPP profiles behind
[`nif_spot`](@ref) (LMF phase-plate table, 2020-10-15). VisRad lists the four plates
among its NIF CPPs (`VR/power_sources/import_nif_cmt_data_file.html`) and tabulates
their radii and exponents (`VR/power_sources/phase_plate_params.html`). Other campaigns
field other plates. This is a default for VisRad parity and is not a facility record.
"""
const VISRAD_NIF_PLATES = Dict(:inner23 => "Rev3-Inner-2w", :inner30 => "1.0 Rev1a-Inner-2w",
                               :outer44 => "Rev3-Outer-2w", :outer50 => "Rev1-Outer-2w")

"""
    visrad_spot(cone::Symbol) -> SuperGaussianSpot

The best-focus spot of the VisRad phase plate for a NIF cone (`VISRAD_NIF_PLATES`, the
plate at rotation 0), with the major axis along φ̂ of the port as for
[`nif_spot`](@ref). VisRad measures the rotation of the spot of a NIF beam from the
horizontal of its quad, so the spots of the VisRad beams are
`nif_beams(; spot = visrad_spot, spot_frame = :quad)`. With
`deposition_kw = (; profile = VisRadProfile())`, the deposition follows the VisRad beam
model (research notes, "VisRad deposition parity").
"""
visrad_spot(cone::Symbol) = visrad_spot(VISRAD_NIF_PLATES[cone])

# ---------------------------------------------------------------------------------------
# Applying per-beam records to a template beam set
# ---------------------------------------------------------------------------------------

# A record is (name, on, theta, phi, aim, defocus, power, spot), with `nothing` for a
# value not given. `spot` is raw VisRad data (model, major, minor, rot, order, level) in
# µm and degrees.
const _VRRecord = NamedTuple{(:name, :on, :theta, :phi, :aim, :defocus, :power, :spot)}

function _vr_merge(recs::Vector)
    aims = [r.aim for r in recs if r.aim !== nothing]
    aim = isempty(aims) ? nothing : sum(aims) / length(aims)
    if length(aims) > 1
        spread = maximum(norm(a - aim) for a in aims)
        spread > 1e-3 && @warn "beams of a quad are pointed differently; using their mean aim" quad_beams = [r.name for r in recs] spread_mm = spread
    end
    dfs = [r.defocus for r in recs if r.defocus !== nothing]
    pws = [r.power for r in recs if r.power !== nothing]
    sps = [r.spot for r in recs if r.spot !== nothing]
    return (name = join((r.name for r in recs), "+"), on = any(r -> r.on, recs),
            theta = nothing, phi = nothing, aim = aim,
            defocus = isempty(dfs) ? nothing : sum(dfs) / length(dfs),
            power = length(pws) == length(recs) ? sum(pws) : nothing,
            spot = isempty(sps) ? nothing : first(sps))
end

# The "Uniform" profile of VisRad is flat across R0 at focus (n(z) = 0,
# VR/power_sources/set_laser_profile_props.html). Briareus approximates it by a
# super-Gaussian of this order with 1/e radius R0 (*approximation*).
const _VR_UNIFORM_ORDER = 50.0

function _vr_spot(raw, tb::Union{Nothing,Beam}, spot_definition, sense)
    raw === nothing && return tb === nothing ? nothing : tb.spot
    raw.major === nothing && return tb === nothing ? nothing : tb.spot
    plate = _vr_plate(raw.model)
    order = raw.order !== nothing ? raw.order :
            occursin(r"uniform"i, raw.model) ? _VR_UNIFORM_ORDER :
            plate !== nothing ? plate.order :
            tb !== nothing ? tb.spot.order :
            (@warn "no super-Gaussian exponent for beam spot; using n = 2" model = raw.model; 2.0)
    level = raw.level !== nothing ? raw.level :
            spot_definition === :auto ? (plate !== nothing ? plate.level : exp(-1)) :
            _vr_level(spot_definition)
    minor = raw.minor === nothing ? raw.major : raw.minor
    h = tb === nothing ? (halo = 0.0, halo_scale = 3.0, halo_order = 2.0) :
        (halo = tb.spot.halo, halo_scale = tb.spot.halo_scale, halo_order = tb.spot.halo_order)
    return _vr_make_spot(raw.major, minor, something(raw.rot, 0.0), order, level, sense; h...,
                         ref = tb === nothing ? 0.0 : _vr_reference_angle(tb))
end

_vr_level(x::Real) = Float64(x)
function _vr_level(s::Symbol)
    s === :e && return exp(-1)
    s === :hwhm && return 0.5
    s === :pct5 && return 0.05
    throw(ArgumentError("spot_definition must be :auto, :e, :hwhm, :pct5 or a fraction of peak"))
end

function _vr_update(b::Beam, r, angle_tol, spot_definition, sense)
    if r.theta !== nothing && r.phi !== nothing && angle_tol < Inf
        δ = acosd(clamp(dot(port_vector(r.theta, r.phi), b.port), -1, 1))
        δ > angle_tol && @warn "VisRad port angles differ from the template beam's" beam = b.name visrad = (r.theta, r.phi) briareus = port_angles(b) deviation_deg = δ
    end
    aim = r.aim === nothing ? b.aim : r.aim
    defocus = r.defocus === nothing ? b.defocus : r.defocus
    power = !r.on ? zero(b.power) : r.power === nothing ? b.power : r.power
    spot = r.spot === nothing ? b.spot : _vr_spot(r.spot, b, spot_definition, sense)
    return remake(b; aim, defocus, power, spot)
end

function _vr_apply(template::AbstractVector{<:Beam}, recs::Dict, absent::Symbol;
                   angle_tol = 0.05, spot_definition = :auto, sense = :beam)
    absent in (:off, :keep, :error) || throw(ArgumentError("absent must be :off, :keep or :error"))
    used = Set{String}()
    members = any(b -> _is_nif_quad(b.name), template) ? _nif_quad_members() : nothing
    out = Beam[]
    for b in template
        k = _vr_key(b.name)
        r = get(recs, k, nothing)
        if r !== nothing
            push!(used, k)
        elseif _is_nif_quad(b.name)
            ks = [_vr_key(n) for n in members[b.name] if haskey(recs, _vr_key(n))]
            if !isempty(ks)
                length(ks) < 4 && @warn "only $(length(ks)) of 4 beams of quad $(b.name) in the file" beams = ks
                union!(used, ks)
                r = _vr_merge([recs[x] for x in ks])
            end
        end
        if r === nothing
            absent === :error && throw(ArgumentError("beam $(b.name) is not in the VisRad file"))
            push!(out, absent === :keep ? b : remake(b; power = zero(b.power)))
        else
            push!(out, _vr_update(b, r, _is_nif_quad(b.name) ? Inf : angle_tol, spot_definition, sense))
        end
    end
    extra = sort!(collect(setdiff(keys(recs), used)))
    isempty(extra) || @warn "VisRad beams without a template beam were ignored" beams = extra
    return identity.(out)
end

function _vr_from_records(recs::Dict, order::Vector{String}, spot_definition, sense)
    out = Beam[]
    for k in order
        r = recs[k]
        (r.theta === nothing || r.phi === nothing) &&
            throw(ArgumentError("beam $(r.name) has no port angles; pass a template beam set"))
        sp = _vr_spot(r.spot, nothing, spot_definition, sense)
        name = occursin(r"^\s*\d+\s*$", r.name) ? k : r.name     # bare OMEGA beam numbers
        push!(out, _vr_quad_framed(Beam(name; theta = r.theta, phi = r.phi,
                                        power = r.on ? something(r.power, 1.0) : 0.0,
                                        aim = something(r.aim, Vec3(0, 0, 0)),
                                        defocus = something(r.defocus, 0.0),
                                        spot = something(sp, SuperGaussianSpot(0.5)))))
    end
    return identity.(out)
end

# ---------------------------------------------------------------------------------------
# NIF SST pointing CSV
# ---------------------------------------------------------------------------------------

"""
    VISRAD_SST_COLUMNS

Column headers of the NIF SST beam-pointing CSV, as VisRad maps them when it imports an
SST file (screenshot in `VR/intro/revision_history.html`, v21.0.0, and
`VR/power_sources/import_laser_beam_data_csv.html`). The columns are the beam name, the
aim point x, y, z (mm, target-chamber Cartesian about TCC) and the defocus (mm).
"""
const VISRAD_SST_COLUMNS = ("Beam", "Pointing: Beam Aim Point X (mm)",
                            "Pointing: Beam Aim Point Y (mm)", "Pointing: Beam Aim Point Z (mm)",
                            "Pointing: Beam Aim Point Beam Defocus (mm)")

"""
    write_visrad_sst(path, beams; include_off = false) -> path

Write beam pointing as a NIF SST CSV (`VISRAD_SST_COLUMNS`). The file has one row per NIF
beam (`"B111"`–`"B468"`), with the aim point in target-chamber Cartesian coordinates
about TCC (mm) *without* the defocus applied, and the defocus (mm; a positive defocus
moves the focus along the propagation direction). This is the export format of VisRad
(`VR/files/laser_beam_data_files.html`, "NIF SST Data Files"). It imports into VisRad
through File | Import Laser Beam Data | Pointing Data, "Custom" CSV, which recognizes
these headers (`VR/power_sources/import_laser_beam_data_csv.html`).

NIF quads (names `"Q11T"`…, e.g. `nif_quads()`) are written as their four beams with
the pointing of the quad. VisRad turns on exactly the beams that have pointing data, so
beams with zero power are left out unless `include_off`.

The file does not contain powers or spots. It is paired with `write_visrad_powers`, and the
spot model is set in VisRad (see `write_visrad_kit`).
"""
function write_visrad_sst(path::AbstractString, beams::AbstractVector{<:Beam};
                          include_off::Bool = false)
    open(path, "w") do io
        println(io, join(VISRAD_SST_COLUMNS, ","))
        for (b, _, _) in _vr_expand(beams)
            (include_off || _vr_v(b.power) > 0) || continue
            println(io, join((_vr_field(b.name), _vr_fmt.(b.aim)..., _vr_fmt(b.defocus)), ","))
        end
    end
    return path
end

"""
    read_visrad_sst(path; template = nif_beams(), absent = :off, units = nothing)
        -> Vector{Beam}

Read a NIF SST pointing CSV and apply it to `template`. The file may come from File |
Export Laser Beam Data | NIF SST Beam Data in VisRad, from the SST, or from
`write_visrad_sst`. Each beam found in the file gets its aim (mm, chamber Cartesian about
TCC, defocus not applied) and defocus (mm), and power and spot are kept. As in the
VisRad import (`VR/power_sources/import_laser_beam_data_csv.html`), template beams
absent from the file are turned off (`absent = :off`, power 0). `:keep` leaves them
unchanged, and `:error` throws.

Beam names match case-insensitively, with `BL115` read as `B115`. With a quad template
(`nif_quads()`), each quad takes the mean aim and defocus of its four beams, with a
warning if they differ by more than 1 µm. The header row is found by its
"Beam Aim Point" columns. Lines before it are skipped, and extra columns are ignored.
Units come from the headers ("(mm)", "(cm)", "(um)"), and `units` (mm per file unit)
overrides them.
"""
function read_visrad_sst(path::AbstractString; template = nif_beams(), absent::Symbol = :off,
                         units = nothing)
    lines = _vr_lines(path)
    hi = findfirst(l -> occursin(r"beam\s+aim\s+point"i, l), lines)
    hi === nothing && throw(ArgumentError("$path has no \"Pointing: Beam Aim Point\" header (NIF SST format)"))
    delim = _vr_delim(lines[hi])
    header = _vr_split(lines[hi], delim)
    h = lowercase.(header)
    find(re) = findfirst(x -> occursin(re, x), h)
    ib = findfirst(x -> occursin(r"^beam\b", x) && !occursin("aim point", x), h)
    ix, iy, iz = find(r"aim\s+point\s+x\b"), find(r"aim\s+point\s+y\b"), find(r"aim\s+point\s+z\b")
    idf = find(r"defocus")
    any(isnothing, (ib, ix, iy, iz)) && throw(ArgumentError("$path: missing beam or aim-point columns"))
    sx = units === nothing ? _vr_length_scale(header[ix], 1.0) : Float64(units)
    sd = units === nothing ? (idf === nothing ? 1.0 : _vr_length_scale(header[idf], 1.0)) : Float64(units)
    recs = Dict{String,Any}()
    for l in lines[hi+1:end]
        isempty(strip(l)) && continue
        f = _vr_split(l, delim)
        length(f) >= max(ib, ix, iy, iz) || continue
        name = f[ib]
        isempty(name) && continue
        xyz = _vr_num.((f[ix], f[iy], f[iz]))
        any(isnothing, xyz) && continue                 # a beam without pointing stays off
        d = idf === nothing || idf > length(f) ? 0.0 : something(_vr_num(f[idf]), 0.0)
        recs[_vr_key(name)] = (name = name, on = true, theta = nothing, phi = nothing,
                               aim = sx * Vec3(xyz...), defocus = sd * d, power = nothing,
                               spot = nothing)
    end
    template === nothing && throw(ArgumentError("SST files carry no port angles: pass a template beam set"))
    return _vr_apply(template, recs, absent)
end

# ---------------------------------------------------------------------------------------
# Laser Beam Summary CSV
# ---------------------------------------------------------------------------------------

# Columns written by `write_visrad_summary`. Through "TC Pt Pos-z" they are the headers
# of the v21 summary dialog (VR/power_sources/images/laser_summary.jpg). The remaining
# headers are hidden in that screenshot and are *assumptions*. The manual lists the
# quantities as power, DPR, profile model, R0 major/minor and minor-axis rotation
# (VR/power_sources/laser_summary_table.html). The last two columns are Briareus
# extensions.
const _VR_SUMMARY_COLUMNS = ("Cluster", "Beam ID", "On/Off", "Polar Ang", "Azm Ang",
                             "Pt Ref-CS", "Pt Geom", "Pt Param-1", "Pt Param-2", "Pt Param-3",
                             "Focus offset", "TC Pt Pos-x", "TC Pt Pos-y", "TC Pt Pos-z",
                             "Power", "DPR", "Spatial Profile", "R0 Major", "R0 Minor",
                             "Rot Angle", "SG Exponent", "R0 Def")

# Header recognition uses the first matching rule not yet used, in this order. The rules
# cover the v21 headers and the older headers of VR/intro/images/laser_summary.jpg
# ("Beam #", "Pt CS", "Pt X-1", "Sp Prof", "Ro").
const _VR_SUMMARY_RULES = (
    :cluster => r"^cluster",
    :beam => r"^beam\b",
    :onoff => r"on\s*/\s*off",
    :polar => r"polar",
    :azim => r"\baz(m|im)",
    :refcs => r"\bref",
    :geom => r"^pt\s*(geom|cs)\b",
    :p1 => r"^pt\s*(param|x)\s*-?\s*1",
    :p2 => r"^pt\s*(param|x)\s*-?\s*2",
    :p3 => r"^pt\s*(param|x)\s*-?\s*3",
    :tcx => r"^tc\s*pt\s*pos\s*-?\s*x",
    :tcy => r"^tc\s*pt\s*pos\s*-?\s*y",
    :tcz => r"^tc\s*pt\s*pos\s*-?\s*z",
    :focus => r"focus",
    :power => r"^(power|pwr)\b",
    :dpr => r"dpr",
    :level => r"\bdef",
    :major => r"maj",
    :minor => r"min",
    :r0 => r"^r[o0]\b",
    :rot => r"rot",
    :order => r"expon|^sg\s*n\b|super",
    :profile => r"prof|model",
)

function _vr_summary_columns(header)
    col = Dict{Symbol,Int}()
    for (j, x) in enumerate(lowercase.(strip.(header)))
        for (k, re) in _VR_SUMMARY_RULES
            haskey(col, k) && continue
            if occursin(re, x)
                col[k] = j
                break
            end
        end
    end
    return col
end

# Returns the pointing reference point from summary-table parameters, when it is given in
# the target-chamber frame (cm, with angles in degrees since v9.0.0,
# VR/intro/revision_versions_800to1020.html).
function _vr_ref_point(cs::AbstractString, geom::AbstractString, p, scale)
    c = lowercase(strip(cs))
    c == "tcc" && return Vec3(0, 0, 0)
    occursin("target chamber", c) || return nothing
    any(isnothing, p) && return nothing
    g = lowercase(strip(geom))
    v = startswith(g, "cyl") ? Vec3(p[1] * cosd(p[3]), p[1] * sind(p[3]), p[2]) :
        startswith(g, "sph") ? p[1] * port_vector(p[2], p[3]) : Vec3(p...)
    return scale * v                       # lengths are scaled, and angles were used above
end

"""
    read_visrad_summary(path; template = nothing, absent = :off, spot_definition = :auto,
                        rotation_sense = :beam, angle_tol = 0.05) -> Vector{Beam}

Read a VisRad Laser Beam Summary CSV (Export Data in the Laser Beams Summary dialog,
`VR/power_sources/laser_summary_table.html`; also written by `write_visrad_summary`).
Columns are recognized by header, so older layouts and extra columns are accepted.

Conversions:
- **Pointing.** "TC Pt Pos-x/y/z" (cm) is the *focus* of the beam in the target-chamber
  frame, which is the Reference Point plus the focus offset along the beam. This follows
  from the example table in the manual (`VR/power_sources/images/laser_summary.jpg`).
  Beams pointed at the same foil-frame Reference Point with focus offset -0.215 cm have
  TC positions that differ by exactly 0.215 cm × the difference of their port vectors.
  So `aim = 10·TC + defocus·port` (mm) with `defocus = 10·focus offset` (cm → mm). When
  the TC columns are absent, the Reference Point is taken from "Pt Param-1/2/3" (cm, as
  Cartesian, cylindrical (r, z, φ) or spherical (r, θ, φ) coordinates, with angles in
  degrees). This works only for references in the target-chamber frame ("Target
  Chamber" or "TCC"). When both are present they are cross-checked.
- **Power.** "Power" is in TW, and for NIF it is the sum of the 3ω, 2ω and 1ω cones
  (`VR/power_sources/nif_beam_params.html`). Time-dependent beams show their peak as
  "T-<peak>", which is read as the peak with a warning. Off beams get power 0.
- **Spot.** "R0" major/minor are radii (µm) at the spot-size definition level, "36.8%
  (1/e)", "50%" or "5%" of peak (`VR/power_sources/set_laser_profile_props.html`). They
  become the 1/e semi-axes `a` (major) and `b` (minor) of the same super-Gaussian
  `exp(-[(x/a)²+(y/b)²]^(n/2))`, which is the focal profile of VisRad. The definition
  comes from an "R0 Def" column if present, and otherwise from `spot_definition`. With
  `:auto`, the level of a recognized built-in phase plate is used (`VISRAD_PHASE_PLATES`;
  NIF CPP radii are HWHM), and 1/e otherwise. 1/e is the VisRad default for Custom
  Supergaussian and an *assumption* for other models. `:e`, `:hwhm`, `:pct5` or a number
  set the level explicitly. The exponent `n` comes from an "SG Exponent" column, else
  from the phase-plate table, else from the template beam. A "Uniform" profile (flat
  across R0 at focus) is approximated by n = 50 with 1/e radius R0.
- **Rotation.** "Rot Angle" is the minor-axis rotation ρ (degrees, clockwise from
  vertical; ρ = 0 puts the major axis "nearly horizontal"). The Briareus `angle` of the
  major axis from θ̂ toward φ̂ is `90° − ρ` for `rotation_sense = :beam` (the spot seen
  along the beam) and `90° + ρ` for `:target` (seen from TCC). The manual does not state
  which. The deposition of VisRad 21 settles it as `:beam`: at ρ = ±30°, 60° and 90° on
  ten NIF beams, its spot axes match `:beam` within 0.04° and differ from `:target` by
  60°. Horizontal is φ̂ of the NIF quad port of the beam (`quad` of the template beam),
  which is turned 1.6–5.3° from φ̂ of the beam, as measured in VisRad 21 to 0.02° on 48
  beams. Beams of no NIF quad use their φ̂.
- DPR, halo, f-number and wavelength are not in the file. DPR "On" is ignored with a
  warning; VisRad splits such OMEGA beams into two cones 85 µm apart. The other values
  come from the template.

With a `template`, beams are matched by name, case-insensitively, with "BL115" read as
"B115" and "10" and "BL-10" as `"Beam_10"`. Port angles are checked against the "Polar
Ang"/"Azm Ang" of the file, with a warning beyond `angle_tol` degrees, which indicates a
frame or naming mismatch. Template beams absent from the file follow `absent` (`:off`,
`:keep`, `:error`), since VisRad omits hidden unused beams. NIF quads aggregate their
four beams, with the mean aim and defocus and the summed power. Without a template,
beams are built from the port angles of the file.
"""
function read_visrad_summary(path::AbstractString; template = nothing, absent::Symbol = :off,
                             spot_definition = :auto, rotation_sense::Symbol = :beam,
                             angle_tol = 0.05)
    _vr_sense(rotation_sense)
    lines = _vr_lines(path)
    hi = findfirst(l -> occursin(r"(^|[,\t;])\s*\"?beam(\s*(id|#))?\s*\"?\s*([,\t;]|$)"i, l), lines)
    hi === nothing && throw(ArgumentError("$path has no Beam ID column header"))
    delim = _vr_delim(lines[hi])
    header = _vr_split(lines[hi], delim)
    col = _vr_summary_columns(header)
    haskey(col, :beam) || throw(ArgumentError("$path: no beam column"))
    if !haskey(col, :major) && haskey(col, :r0)
        col[:major] = col[:r0]
    end
    unit(k, default) = haskey(col, k) ? _vr_length_scale(header[col[k]], default) : default
    s_tc, s_p, s_f = unit(:tcx, 10.0), unit(:p1, 10.0), unit(:focus, 10.0)
    s_r = unit(:major, 1e-3) * 1000                 # radii to µm
    recs = Dict{String,Any}()
    order = String[]
    timedep = String[]
    for l in lines[hi+1:end]
        isempty(strip(l)) && continue
        f = _vr_split(l, delim)
        fld(k) = haskey(col, k) && col[k] <= length(f) ? f[col[k]] : ""
        num(k) = _vr_num(fld(k))
        name = fld(:beam)
        isempty(name) && continue
        on = haskey(col, :onoff) ? lowercase(fld(:onoff)) in ("on", "yes", "true", "1") : true
        θ, φ = num(:polar), num(:azim)
        pstr = fld(:power)
        if occursin(r"^t-"i, pstr)
            push!(timedep, name)
            pstr = pstr[3:end]
        end
        P = _vr_num(pstr)
        d = num(:focus)
        defocus = d === nothing ? nothing : s_f * d
        tc = (num(:tcx), num(:tcy), num(:tcz))
        ref = _vr_ref_point(fld(:refcs), fld(:geom), (num(:p1), num(:p2), num(:p3)), s_p)
        port = θ === nothing || φ === nothing ? nothing : port_vector(θ, φ)
        aim = nothing
        if !any(isnothing, tc) && port !== nothing
            aim = s_tc * Vec3(tc...) + something(defocus, 0.0) * port
            if ref !== nothing && norm(aim - ref) > 2e-3 + 1e-3 * norm(ref)
                @warn "summary row: TC Pt Pos disagrees with Pt Param + focus offset (check units)" beam = name tc_aim = aim param_aim = ref
            end
        elseif ref !== nothing
            aim = ref
        else
            throw(ArgumentError("$path: beam $name is pointed relative to \"$(fld(:refcs))\"; " *
                                "its chamber position needs the TC Pt Pos columns"))
        end
        lowercase(fld(:dpr)) == "on" && @warn "DPR is on for beam $name; Briareus does not model DPR"
        rot = num(:rot)
        lv = fld(:level)
        pm = match(r"([\d.]+)\s*%", lv)
        level = occursin(r"1/e", lv) ? exp(-1) : occursin(r"hwhm|half"i, lv) ? 0.5 :
                pm !== nothing ? parse(Float64, pm[1]) / 100 : nothing
        major = num(:major)
        spot = major === nothing ? nothing :
               (model = fld(:profile), major = s_r * major,
                minor = (mi = num(:minor)) === nothing ? nothing : s_r * mi,
                rot = rot, order = num(:order), level = level)
        k = _vr_key(name)
        haskey(recs, k) || push!(order, k)
        recs[k] = (name = name, on = on, theta = θ, phi = φ, aim = aim, defocus = defocus,
                   power = P === nothing ? nothing : (on ? P : 0.0), spot = spot)
    end
    isempty(timedep) || @warn "time-dependent beams: \"Power\" is the peak of their table (import the power file for the history)" beams = timedep
    template === nothing && return _vr_from_records(recs, order, spot_definition, rotation_sense)
    return _vr_apply(template, recs, absent; angle_tol, spot_definition, sense = rotation_sense)
end

"""
    write_visrad_summary(path, beams; rotation_sense = :beam, cluster = b -> string(b.group),
                         profile = "Custom Supergaussian", extras = true) -> path

Write beams in the layout of the VisRad Laser Beam Summary CSV, one row per VisRad beam.
NIF quads expand to their four beams with a quarter of the power each. The columns are
the cluster, beam ID, on/off (power > 0), port polar angle and azimuth (degrees), the
pointing reference "Target Chamber"/"Cartesian" with the aim as Pt Param-1/2/3 (cm), the
focus offset (cm), the focus position "TC Pt Pos" (cm, aim + defocus·k), the power (TW),
DPR "Off", the spatial profile model, R0 major/minor as 1/e radii (µm), and the
minor-axis rotation (degrees; see `read_visrad_summary` for `rotation_sense`). `extras`
appends "SG Exponent" and "R0 Def" ("36.8% (1/e)") so that the file is self-contained.

Columns through "TC Pt Pos-z" reproduce the headers of the v21 dialog. The names of the
later columns are not shown in the manual (*assumption*), and halo parameters are not
written, since VisRad has none. The file serves for records and for
`read_visrad_summary`. Pointing is set in VisRad by importing the output of
`write_visrad_sst`, whose format is documented.
"""
function write_visrad_summary(path::AbstractString, beams::AbstractVector{<:Beam};
                              rotation_sense::Symbol = :beam, cluster = b -> string(b.group),
                              profile = "Custom Supergaussian", extras::Bool = true)
    _vr_sense(rotation_sense)
    ncol = extras ? length(_VR_SUMMARY_COLUMNS) : length(_VR_SUMMARY_COLUMNS) - 2
    open(path, "w") do io
        println(io, join(_VR_SUMMARY_COLUMNS[1:ncol], ","))
        for (b, _, _) in _vr_expand(beams)
            θ, φ = port_angles(b)
            aim = _vr_v.(b.aim)
            d = _vr_v(b.defocus)
            focus = aim + d * direction(b)
            sp = b.spot
            a, bb = _vr_v(sp.a), _vr_v(sp.b)
            sp.halo > 0 && @warn "spot halo of beam $(b.name) is not representable in VisRad; dropped" maxlog = 1
            row = Any[cluster(b), b.name, _vr_v(b.power) > 0 ? "On" : "Off", θ, φ,
                      "Target Chamber", "Cartesian", (aim ./ 10)..., d / 10, (focus ./ 10)...,
                      _vr_v(b.power), "Off", profile, 1000 * max(a, bb), 1000 * min(a, bb),
                      _vr_rotation(sp, rotation_sense, _vr_reference_angle(b))]
            extras && append!(row, (sp.order, "36.8% (1/e)"))
            println(io, join(_vr_field.(row), ","))
        end
    end
    return path
end

# ---------------------------------------------------------------------------------------
# Power files
# ---------------------------------------------------------------------------------------

_vr_scale(p::PulseShape, f) = f == 1 ? p : PulseShape(p.t, f .* p.p)

"""
    write_visrad_powers(path, beams, pulses; by = :group, times = nothing) -> path
    write_visrad_powers(path, pulses::AbstractDict{Symbol,PulseShape}; times = nothing) -> path

Write a VisRad time-dependent beam-power CSV (`VR/files/laser_beam_data_files.html`).
The header is `time (ns)`, followed by one column per beam named as VisRad names beams
(`B115` for NIF, `Beam_10` for OMEGA), and the rows hold time (ns) and power (TW per
beam). The file is imported with File | Import Laser Beam Data | Time-Dependent Powers
(`VR/power_sources/import_tdep_beam_powers.html`). For NIF, VisRad assigns the powers to
the 3ω (blue) cones.

`pulses` takes the forms of [`with_pulses`](@ref): a `PulseShape` for every beam, or a
`Dict`/NamedTuple keyed `by` group (default), quad, beam name (`by = :beam`, as
[`read_visrad_powers`](@ref) returns them) or a function of the beam. The second form
writes a `Dict` keyed by beam name as it is. NIF quads expand to their four beams at a
quarter of the quad power each. Times default to the union of all pulse knots, so
piecewise-linear pulses are reproduced exactly when each starts and ends at zero power.
A `PulseShape` is zero outside its table, so a pulse that starts or ends at a finite
power gets a ramp to the knots of its neighbors instead. Briareus beam names that are
not facility names, such as the hohlraum-frame `B01`… of `omega_beams()`, are written
as they are and should be renamed first.
"""
function write_visrad_powers(path::AbstractString, beams::AbstractVector{<:Beam}, pulses;
                             by = :group, times = nothing)
    key = _beam_key(by)
    cols = Pair{String,PulseShape}[]
    for (b, i, f) in _vr_expand(beams)
        push!(cols, b.name => _vr_scale(_beam_pulse(pulses, key, beams[i]), f))
    end
    return _vr_write_powers(path, cols, times)
end

write_visrad_powers(path::AbstractString, pulses::AbstractDict{Symbol,PulseShape};
                    times = nothing) =
    _vr_write_powers(path, [String(k) => pulses[k] for k in sort!(collect(keys(pulses)))], times)

function _vr_write_powers(path, cols, times)
    ts = times === nothing ? sort!(unique!(reduce(vcat, (p.t for (_, p) in cols); init = Float64[]))) :
         collect(Float64, times)
    open(path, "w") do io
        println(io, join(("time (ns)", (_vr_field(n) for (n, _) in cols)...), ","))
        for t in ts
            println(io, join((_vr_fmt(t), (_vr_fmt(p(t)) for (_, p) in cols)...), ","))
        end
    end
    return path
end

"""
    read_visrad_powers(path; template = nothing) -> Dict{Symbol,PulseShape}

Read a VisRad time-dependent beam-power CSV. The first column is time in ns, and a
"(s)"/"(sec)" or "(ps)" unit in its header is converted. One power column (TW) per beam
follows, headed by the beam name (`VR/files/laser_beam_data_files.html`). Negative
powers are set to zero, as VisRad does on import
(`VR/power_sources/import_tdep_beam_powers.html`).

Without a template the keys are the normalized beam names as `Symbol`s, so that `BL115`
becomes `:B115`, and OMEGA beams become keys such as `:Beam_10`. With a `template` beam
set the keys are the names of the template. Beams match by name, and a NIF quad gets the
sum of the pulses of its four beams, on the union of their knots. Template beams absent
from the file have no key. The beams at time `t` are
`with_pulses(beams, P, t; by = :beam)` ([`with_pulses`](@ref)).
"""
function read_visrad_powers(path::AbstractString; template = nothing)
    lines = filter(l -> !isempty(strip(l)), _vr_lines(path))
    hi = findfirst(l -> occursin(r"^\s*\"?time"i, l), lines)
    hi = hi === nothing ? 1 : hi
    delim = _vr_delim(lines[hi])
    header = _vr_split(lines[hi], delim)
    th = lowercase(header[1])
    ts = occursin(r"\(\s*(s|sec)\s*\)", th) ? 1e9 : occursin(r"\(\s*ps\s*\)", th) ? 1e-3 : 1.0
    rows = [_vr_split(l, delim) for l in lines[hi+1:end]]
    rows = [r for r in rows if _vr_num(r[1]) !== nothing]
    t = [ts * _vr_num(r[1]) for r in rows]
    p = sortperm(t)
    out = Dict{String,PulseShape}()   # by normalized file name
    for j in 2:length(header)
        isempty(header[j]) && continue
        P = [j <= length(r) ? max(something(_vr_num(r[j]), 0.0), 0.0) : 0.0 for r in rows]
        out[_vr_key(header[j])] = PulseShape(t[p], P[p])
    end
    template === nothing && return Dict{Symbol,PulseShape}(Symbol(k) => v for (k, v) in out)
    res = Dict{Symbol,PulseShape}()
    members = any(b -> _is_nif_quad(b.name), template) ? _nif_quad_members() : nothing
    for b in template
        k = _vr_key(b.name)
        if haskey(out, k)
            res[Symbol(b.name)] = out[k]
        elseif _is_nif_quad(b.name)
            ps = [out[_vr_key(n)] for n in members[b.name] if haskey(out, _vr_key(n))]
            isempty(ps) && continue
            length(ps) < 4 && @warn "only $(length(ps)) of 4 beams of quad $(b.name) in the power file"
            tt = sort!(unique!(reduce(vcat, (q.t for q in ps))))
            res[Symbol(b.name)] = PulseShape(tt, [sum(q(τ) for q in ps) for τ in tt])
        end
    end
    return res
end

"""
    read_visrad_pulse_shape(path) -> PulseShape

Read one of the standard laser pulse-shape files of VisRad (`data/laser_pulse_shapes` of
the installation), a two-column ASCII table of power (W) versus time (s)
(`VR/files/laser_pulse_shape_files.html`), and return it in ns and TW. The column order
is not stated. Time first is an *assumption*, and the columns are swapped when the first
column is not sorted while the second is. Non-numeric lines are skipped.
"""
function read_visrad_pulse_shape(path::AbstractString)
    xs, ys = Float64[], Float64[]
    for l in _vr_lines(path)
        f = split(strip(l), r"[\s,;]+")
        length(f) >= 2 || continue
        a, b = tryparse(Float64, f[1]), tryparse(Float64, f[2])
        (a === nothing || b === nothing) && continue
        push!(xs, a); push!(ys, b)
    end
    isempty(xs) && throw(ArgumentError("$path has no numeric rows"))
    t, P = (!issorted(xs) && issorted(ys)) ? (ys, xs) : (xs, ys)
    return PulseShape(1e9 .* t, 1e-12 .* P)
end
