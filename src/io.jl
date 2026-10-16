# Beam-table and result I/O (plain CSV, no dependencies).
#
# The canonical beam table has one row per beam with the columns in `BEAM_COLUMNS`
# (lengths mm, power TW, angles degrees, wavelength µm, detuning Å). Foreign tables
# (e.g. exported from VisRad or a shot setup sheet) are read by mapping their column
# names and units onto these fields.

"Column names of the canonical beam table (see `write_beam_table`)."
const BEAM_COLUMNS = (:name, :quad, :group, :theta, :phi, :power, :aim_x, :aim_y, :aim_z,
                      :defocus, :spot_a, :spot_b, :spot_order, :spot_angle, :spot_halo,
                      :spot_halo_scale, :spot_halo_order, :fnumber,
                      :wavelength, :dlambda)

"""
    write_beam_table(path, beams) -> path

Write beams as CSV with the canonical columns `$(join(BEAM_COLUMNS, ", "))`
(mm, TW, degrees, µm, Å). Hohlraum-cylindrical aim coordinates `aim_r` (toward the
port azimuth), `aim_phi_offset` (tangential) are appended for reference.
"""
function write_beam_table(path::AbstractString, beams::AbstractVector{<:Beam})
    v = ForwardDiff.value
    open(path, "w") do io
        println(io, join((BEAM_COLUMNS..., :aim_r, :aim_tangential), ","))
        for b in beams
            θ, φ = port_angles(b)
            er, eφ = _radial(b)
            aim = Vec3(v.(b.aim))
            row = (b.name, b.quad, b.group, θ, φ, v(b.power), aim..., v(b.defocus),
                   v(b.spot.a), v(b.spot.b), b.spot.order, rad2deg(b.spot.angle), b.spot.halo,
                   b.spot.halo_scale, b.spot.halo_order, b.fnumber,
                   b.wavelength, v(b.dlambda), dot(aim, er), dot(aim, eφ))
            println(io, join((_csvfield(x) for x in row), ","))
        end
    end
    return path
end

_csvfield(x::AbstractFloat) = @sprintf("%.10g", x)
_csvfield(x::Symbol) = _csvfield(string(x))
_csvfield(x::AbstractString) = occursin(r"[,\"]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : x
_csvfield(x) = string(x)

"""
    read_csv(path; comment = '#') -> (header::Vector{String}, rows::Vector{Vector{String}})

Minimal CSV reader (quoted fields, `""` escapes, comment and blank lines skipped).
Tab- or semicolon-separated files are detected from the header line.
"""
function read_csv(path::AbstractString; comment = '#')
    lines = [l for l in eachline(path) if !isempty(strip(l)) && !startswith(lstrip(l), comment)]
    isempty(lines) && throw(ArgumentError("$path has no rows"))
    delim = occursin('\t', lines[1]) ? '\t' : (count(==(';'), lines[1]) > count(==(','), lines[1]) ? ';' : ',')
    split_line(l) = begin
        out = String[]; buf = IOBuffer(); q = false; i = 1
        cs = collect(l)
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
        out
    end
    return split_line(lines[1]), [split_line(l) for l in lines[2:end]]
end

"""
    read_beam_table(path; columns = Dict(), scale = Dict(), defaults = (;), template = nothing)
        -> Vector{Beam}

Read beams from a CSV. Canonical column names (see `write_beam_table`) are recognized
case-insensitively; `columns` maps a canonical field to the file's column name for
foreign tables, e.g. `Dict(:aim_x => "X (cm)", :theta => "Theta")`, and `scale`
multiplies a field after reading, e.g. `Dict(:aim_x => 10, :aim_y => 10, :aim_z => 10)`
for cm → mm. Fields missing from the file come from `defaults` (a NamedTuple of
canonical fields), or — when `template` is a beam vector — from the template beam with
the same `name` (so a pointing table can update a facility beam set), else from the
`Beam` defaults.
"""
function read_beam_table(path::AbstractString; columns = Dict{Symbol,String}(),
                         scale = Dict{Symbol,Float64}(), defaults = (;), template = nothing)
    header, rows = read_csv(path)
    lower = lowercase.(header)
    colidx = Dict{Symbol,Int}()
    for f in BEAM_COLUMNS
        name = lowercase(get(columns, f, string(f)))
        i = findfirst(==(name), lower)
        i === nothing || (colidx[f] = i)
    end
    haskey(colidx, :name) || throw(ArgumentError("no beam name column in $path (map it with columns)"))
    tmpl = template === nothing ? Dict{String,Beam}() : Dict(b.name => b for b in template)
    beams = Beam{Float64}[]
    for r in rows
        name = r[colidx[:name]]
        tb = get(tmpl, name, nothing)
        get_num(f, fallback) = begin
            x = haskey(colidx, f) && !isempty(r[colidx[f]]) ? parse(Float64, r[colidx[f]]) :
                haskey(defaults, f) ? Float64(defaults[f]) : fallback
            x === nothing ? nothing : x * get(scale, f, 1.0)
        end
        get_str(f, fallback) = haskey(colidx, f) && !isempty(r[colidx[f]]) ? r[colidx[f]] :
                               haskey(defaults, f) ? string(defaults[f]) : fallback
        θ0, φ0 = tb === nothing ? (nothing, nothing) : port_angles(tb)
        θ = get_num(:theta, θ0)
        φ = get_num(:phi, φ0)
        (θ === nothing || φ === nothing) && throw(ArgumentError("beam $name has no port angles"))
        tv(f, g) = tb === nothing ? g : f(tb)
        a = get_num(:spot_a, tv(b -> b.spot.a, 0.5))
        spot = SuperGaussianSpot(a, get_num(:spot_b, tv(b -> b.spot.b, a));
                                 order = get_num(:spot_order, tv(b -> b.spot.order, 2.0)),
                                 angle = deg2rad(get_num(:spot_angle, tv(b -> rad2deg(b.spot.angle), 0.0))),
                                 halo = get_num(:spot_halo, tv(b -> b.spot.halo, 0.0)),
                                 halo_scale = get_num(:spot_halo_scale, tv(b -> b.spot.halo_scale, 3.0)),
                                 halo_order = get_num(:spot_halo_order, tv(b -> b.spot.halo_order, 2.0)))
        aim0 = tv(b -> b.aim, Vec3(0, 0, 0))
        push!(beams, Beam(name; theta = θ, phi = φ,
                          power = get_num(:power, tv(b -> b.power, 1.0)),
                          aim = (get_num(:aim_x, aim0[1]), get_num(:aim_y, aim0[2]), get_num(:aim_z, aim0[3])),
                          defocus = get_num(:defocus, tv(b -> b.defocus, 0.0)), spot,
                          fnumber = get_num(:fnumber, tv(b -> b.fnumber, 20.0)),
                          group = Symbol(get_str(:group, string(tv(b -> b.group, :none)))),
                          quad = get_str(:quad, tv(b -> b.quad, "")),
                          wavelength = get_num(:wavelength, tv(b -> b.wavelength, 0.351)),
                          dlambda = get_num(:dlambda, tv(b -> b.dlambda, 0.0))))
    end
    return beams
end

"""
    write_element_table(path, target, fields...; elems = all) -> path

Per-element CSV: index, component, side, centroid (mm), normal, area (mm²), polar and
azimuthal angle of the centroid (degrees, about TCC), then one column per
`name => values` pair in `fields` (e.g. `:H => sol.H`, `:laser => laser_flux(dep, t)`).
"""
function write_element_table(path::AbstractString, t::Target, fields::Pair...;
                             elems = 1:nelements(t))
    open(path, "w") do io
        println(io, join(("element", "component", "side", "x", "y", "z", "nx", "ny", "nz", "area",
                          "theta_deg", "phi_deg", (string(first(f)) for f in fields)...), ","))
        for e in elems
            c, n = t.elem_centroid[e], t.elem_normal[e]
            θ = acosd(clamp(c[3] / max(norm(c), eps()), -1, 1))
            φ = mod(atand(c[2], c[1]), 360)
            vals = (ForwardDiff.value(last(f)[e]) for f in fields)
            println(io, join((e, t.components[t.elem_comp[e]].name, t.elem_side[e],
                              _csvfield.((c..., n..., t.elem_area[e], θ, φ))..., _csvfield.(vals)...), ","))
        end
    end
    return path
end

"""
    write_moments(path, moments; label = "a") -> path

Write Legendre moments `aₗ` and `aₗ/a₀` (vector) or a time history (matrix
`(lmax+1) × nt` with `times`) as CSV.
"""
function write_moments(path::AbstractString, a::AbstractVector)
    open(path, "w") do io
        println(io, "l,a_l,a_l_over_a0")
        for (i, x) in enumerate(a)
            println(io, i - 1, ",", _csvfield(ForwardDiff.value(x)), ",",
                    _csvfield(ForwardDiff.value(x / a[1])))
        end
    end
    return path
end

function write_moments(path::AbstractString, times::AbstractVector, M::AbstractMatrix)
    open(path, "w") do io
        println(io, join(("time_ns", "a0", ("P$(l)_over_P0" for l in 1:size(M, 1)-1)...), ","))
        for (k, τ) in enumerate(times)
            println(io, join(_csvfield.((τ, M[1, k], (M[l+1, k] / M[1, k] for l in 1:size(M, 1)-1)...)), ","))
        end
    end
    return path
end
