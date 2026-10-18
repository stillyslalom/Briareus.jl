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

"""
    read_stl(path; scale = 1.0) -> Mesh

Triangle mesh from an STL file (binary or ASCII), one patch per facet, with vertex
coordinates multiplied by `scale` (STL carries no units: 1 for millimetres, 10 for
centimetres, 25.4 for inches). Facet normals in the file are ignored; the winding
decides the front side, as STL specifies (counterclockwise seen from outside), but
exported CAD often mixes windings: pass the mesh through `orient_mesh`, and through
`coarsen_mesh` before view factors.

All facets of the file go into one mesh, also when an ASCII file holds several
`solid … endsolid` blocks (one per CAD body, as many exporters write them);
[`read_stl_solids`](@ref) keeps them apart. A binary file is recognized by its size
(84 + 50n bytes for the n facets its header announces), so binary files whose header
begins with "solid" (as SolidWorks writes them) read as binary. ASCII keywords are read
in any case, with any whitespace and line endings.
"""
function read_stl(path::AbstractString; scale::Real = 1.0)
    solids = read_stl_solids(path; scale)
    return length(solids) == 1 ? solids[1][2] : merge_meshes(last.(solids)...)
end

"""
    read_stl_solids(path; scale = 1.0) -> Vector{Pair{String,Mesh}}

The solids of an STL file as `name => mesh` pairs in file order, each mesh as
[`read_stl`](@ref) reads it (one patch per facet, coordinates multiplied by `scale`).
An ASCII file gives one pair per `solid <name> … endsolid` block, named by the rest
of its `solid` line (`"solid<k>"` for the k-th solid when that is blank; repeated
names are kept as they are), and an empty solid gives an empty mesh. A binary file has
no solids: it gives one pair named after the file. To build a target from the parts,
e.g. `Target([Component(n, orient_mesh(m), gold()) for (n, m) in read_stl_solids(path)])`.
"""
function read_stl_solids(path::AbstractString; scale::Real = 1.0)
    data = read(path)
    n = length(data) >= 84 ? Int(_stl_u32(data, 81)) : -1
    if n >= 0 && length(data) == 84 + 50n
        return [splitext(basename(path))[1] => _stl_binary(data, n, Float64(scale))]
    end
    i = findfirst(b -> !_stl_space(b), data)
    ascii = i !== nothing && i + 4 <= length(data) &&
            all(k -> lowercase(Char(data[i+k])) == "solid"[k+1], 0:4)
    if !ascii || findfirst(iszero, data) !== nothing     # ASCII STL never holds a NUL byte
        throw(ArgumentError(length(data) >= 84 ?
            "$path is not an STL file: as binary STL its $(length(data)) bytes do not match " *
            "the 84 + 50×$n its header announces (truncated?), and it is not ASCII STL" :
            "$path is not an STL file (neither binary nor ASCII)"))
    end
    return _stl_ascii(String(data), Float64(scale))
end

_stl_u32(d::Vector{UInt8}, i::Int) =
    UInt32(d[i]) | UInt32(d[i+1]) << 8 | UInt32(d[i+2]) << 16 | UInt32(d[i+3]) << 24
_stl_f32(d::Vector{UInt8}, i::Int) = Float64(reinterpret(Float32, _stl_u32(d, i)))
_stl_space(b::UInt8) = b == UInt8(' ') || UInt8('\t') <= b <= UInt8('\r')

# One patch per facet: vertices 3k-2:3k, triangle k, patch k:k.
function _stl_mesh(vertices::Vector{Vec3})
    n = length(vertices) ÷ 3
    return Mesh(vertices, [(3k - 2, 3k - 1, 3k) for k in 1:n], [k:k for k in 1:n])
end

function _stl_binary(data::Vector{UInt8}, n::Integer, scale::Float64)
    vs = Vector{Vec3}(undef, 3n)
    @inbounds for k in 1:n, j in 1:3
        o = 84 + 50(k - 1) + 12j + 1               # skip the facet normal
        vs[3(k-1)+j] = scale * Vec3(_stl_f32(data, o), _stl_f32(data, o + 4), _stl_f32(data, o + 8))
    end
    return _stl_mesh(vs)
end

# ASCII STL by a byte-level tokenizer over the file's text: no per-line or per-token
# allocation (numbers parse from SubStrings).
function _stl_token(d, i::Int, line::Int)           # next token d[i0:i1] from i; i0 > end at EOF
    len = length(d)
    while i <= len && _stl_space(d[i])
        d[i] == UInt8('\n') && (line += 1)
        i += 1
    end
    j = i
    while j <= len && !_stl_space(d[j])
        j += 1
    end
    return i, j - 1, line
end

_stl_iskey(d, i0::Int, i1::Int, key::String) =     # case-insensitive keyword match
    i1 - i0 + 1 == ncodeunits(key) && all(k -> d[i0+k] | 0x20 == codeunit(key, k + 1), 0:ncodeunits(key)-1)

function _stl_eol(d, i::Int)                        # index of the end of i's line
    while i <= length(d) && d[i] != UInt8('\n')
        i += 1
    end
    return i
end

const _STL_POW10 = ntuple(k -> Float64(big(10)^(k - 1)), 23)   # 1e0 … 1e22, exact

# Clinger's fast path: a decimal with an integer mantissa m ≤ 2⁵³ and a power of ten
# |e| ≤ 22 is m·10^e or m/10^-e of two exact doubles, so one IEEE operation rounds it
# correctly. The usual `%e` output of CAD exporters takes it; anything else (more
# digits, inf, nan, malformed) goes to `tryparse` (strtod: ~0.4 µs per number on
# Windows, ten times the fast path).
function _stl_float(text::String, i0::Int, i1::Int)
    d = codeunits(text)
    i = i0
    neg = d[i] == UInt8('-')
    (neg || d[i] == UInt8('+')) && (i += 1)
    m = UInt64(0)
    nd = 0                                          # mantissa digits (leading zeros skipped)
    e = 0
    dot = false
    any = false
    while i <= i1
        c = d[i]
        if UInt8('0') <= c <= UInt8('9')
            any = true
            (m > 0 || c != UInt8('0')) && (nd += 1)
            nd > 18 && return tryparse(Float64, SubString(text, i0, i1))
            m = 10m + (c - UInt8('0'))
            dot && (e -= 1)
        elseif c == UInt8('.') && !dot
            dot = true
        else
            break
        end
        i += 1
    end
    any || return tryparse(Float64, SubString(text, i0, i1))
    if i <= i1                                      # exponent
        d[i] | 0x20 == UInt8('e') && i < i1 || return tryparse(Float64, SubString(text, i0, i1))
        i += 1
        eneg = d[i] == UInt8('-')
        (eneg || d[i] == UInt8('+')) && (i += 1)
        i <= i1 || return tryparse(Float64, SubString(text, i0, i1))
        x = 0
        while i <= i1
            c = d[i]
            UInt8('0') <= c <= UInt8('9') && x < 10000 || return tryparse(Float64, SubString(text, i0, i1))
            x = 10x + (c - UInt8('0'))
            i += 1
        end
        e += eneg ? -x : x
    end
    m <= UInt64(2)^53 && -22 <= e <= 22 || return tryparse(Float64, SubString(text, i0, i1))
    v = e >= 0 ? Float64(m) * _STL_POW10[e+1] : Float64(m) / _STL_POW10[-e+1]
    return neg ? -v : v
end

function _stl_number(text::String, i::Int, line::Int)
    i0, i1, line = _stl_token(codeunits(text), i, line)
    v = i0 > i1 ? nothing : _stl_float(text, i0, i1)
    v === nothing && throw(ArgumentError("STL line $line: malformed vertex"))
    return v, i1 + 1, line
end

function _stl_close!(solids, name, vs, nv::Int, line::Int)
    nv == 0 || throw(ArgumentError("STL line $line: a facet must have 3 vertices, not $nv"))
    name === nothing || push!(solids, name => _stl_mesh(vs))
    return nothing
end

function _stl_ascii(text::String, scale::Float64)
    d = codeunits(text)
    len = length(d)
    solids = Pair{String,Mesh}[]
    vs = Vec3[]                                     # vertices of the current solid
    nv = 0                                          # vertices of the current facet
    name::Union{Nothing,String} = nothing           # current solid (nothing: between solids)
    line = 1
    i = 1
    while true
        i0, i1, line = _stl_token(d, i, line)
        i0 > len && break
        i = i1 + 1
        if _stl_iskey(d, i0, i1, "vertex")
            x, i, line = _stl_number(text, i, line)
            y, i, line = _stl_number(text, i, line)
            z, i, line = _stl_number(text, i, line)
            name === nothing && (name = "solid$(length(solids) + 1)")   # facets outside a solid
            push!(vs, scale * Vec3(x, y, z))
            nv += 1
        elseif _stl_iskey(d, i0, i1, "endfacet")
            nv == 3 || throw(ArgumentError("STL line $line: a facet must have 3 vertices, not $nv"))
            nv = 0
        elseif _stl_iskey(d, i0, i1, "solid")
            _stl_close!(solids, name, vs, nv, line)
            vs = Vec3[]
            j = _stl_eol(d, i)                      # the name: the rest of the line
            nm = strip(SubString(text, i, prevind(text, j)))
            name = isempty(nm) ? "solid$(length(solids) + 1)" : String(nm)
            i = j
        elseif _stl_iskey(d, i0, i1, "endsolid")
            _stl_close!(solids, name, vs, nv, line)
            vs = Vec3[]
            name = nothing
            i = _stl_eol(d, i)                      # skip the repeated name
        end                                         # facet normal, outer loop, endloop: skipped
    end
    (name !== nothing || !isempty(vs)) && _stl_close!(solids, name, vs, nv, line)  # no endsolid
    return solids
end

"""
    write_stl(path, x; binary = true, scale = 1.0, apertures = true, name = "briareus") -> path

Write `x` as an STL file for CAD or another code: a `Mesh`, a `Component`, a
`Target`, or a vector of components or of `name => mesh` pairs. Every triangle
becomes a facet whose normal is computed from its winding (the front side); vertex
coordinates are divided by `scale` (STL carries no units: 1 writes millimetres, 10
centimetres, 25.4 inches), so `read_stl(path; scale)` reads the same geometry back.

`binary = true` writes binary STL (Float32 coordinates, so a coordinate of size L
moves by up to ~6e-8 L; no names: the header lists them). `binary = false` writes ASCII
STL with one `solid <name> … endsolid <name>` block per component (or pair), in order,
at full Float64 precision; [`read_stl_solids`](@ref) reads them back as named parts
(whitespace in a name is written as `_`). A mesh alone is named `name`.

A two-sided component is written once, wound toward its front side. Apertures (laser
entrance holes and other open windows, transparent to light) are written as their own
solids so a CAD user sees them; pass `apertures = false` to leave them out, e.g. for a
code that would take every facet for a wall. A component's exact `CurvedSurface` is
not written, only its facets, and the patch structure is lost (`read_stl` returns one
patch per triangle).
"""
function write_stl(path::AbstractString, x; binary::Bool = true, scale::Real = 1.0,
                   apertures::Bool = true, name::AbstractString = "briareus")
    solids = _stl_solids(x, apertures, String(name))
    s = Float64(scale)
    s > 0 || throw(ArgumentError("scale must be positive, not $scale"))
    io = IOBuffer()
    binary ? _write_stl_binary(io, solids, s) : _write_stl_ascii(io, solids, s)
    write(path, take!(io))
    return path
end

_stl_solids(m::Mesh, _, name) = [name => m]
_stl_solids(c::Component, apertures, name) = _stl_solids([c], apertures, name)
_stl_solids(t::Target, apertures, name) = _stl_solids(t.components, apertures, name)
_stl_solids(cs::AbstractVector{Component}, apertures, _) =
    [String(c.name) => c.mesh for c in cs if apertures || c.kind !== :aperture]
_stl_solids(ps::AbstractVector{<:Pair}, _, _) = [string(n) => m::Mesh for (n, m) in ps]

function _stl_facet(m::Mesh, t::Integer, s::Float64)
    a, b, c = triangle_vertices(m, t)
    n = cross(b - a, c - a)
    l = norm(n)
    return (l > 0 ? n / l : zero(Vec3)), a / s, b / s, c / s
end

function _write_stl_binary(io::IO, solids, s::Float64)
    # 80 header bytes, not starting with "solid": the units and the solids' names
    header = Vector{UInt8}("Briareus.jl binary STL (units: mm/$s): " * join(first.(solids), ", "))
    nh = length(header)
    resize!(header, 80)
    header[min(nh, 80)+1:80] .= UInt8(' ')
    write(io, header)
    n = sum(ntriangles(m) for (_, m) in solids; init = 0)
    write(io, htol(UInt32(n)))
    for (_, m) in solids, t in eachindex(m.triangles)
        for v in _stl_facet(m, t, s), x in v
            write(io, htol(Float32(x)))
        end
        write(io, UInt16(0))
    end
end

function _write_stl_ascii(io::IO, solids, s::Float64)
    for (name, m) in solids
        nm = replace(name, r"\s+" => "_")         # a name is one token on the endsolid line
        println(io, "solid ", nm)
        for t in eachindex(m.triangles)
            n, a, b, c = _stl_facet(m, t, s)
            println(io, "  facet normal ", n[1], ' ', n[2], ' ', n[3])
            println(io, "    outer loop")
            for v in (a, b, c)
                println(io, "      vertex ", v[1], ' ', v[2], ' ', v[3])
            end
            println(io, "    endloop")
            println(io, "  endfacet")
        end
        println(io, "endsolid ", nm)
    end
end

"""
    write_results(dir, target, dep; sol = nothing, beams = nothing, observe = nothing,
                  lmax = 8) -> Vector{String}

Write a calculation's results to the folder `dir` as CSV files, for postprocessing or
comparison with another code, and return their paths:

- `elements.csv` — per element (`write_element_table`): absorbed laser power (TW) and
  intensity (W/cm²), and with a radiosity solution `sol` the x-ray source, radiosity
  and incident flux (W/cm²), albedo and transmitted fraction;
- `components.csv` — per component: absorbed laser power (TW, front and the unmodeled
  back side), and with `sol` the x-ray power absorbed, escaped (apertures) and
  transmitted (thin foils, gaps) from `energy_balance`;
- `beams.csv` — with `beams` (in the deposition's order): each beam's name, port angles,
  delivered and landed power (TW);
- `legendre.csv` and `ylm.csv` — with `observe` (component names) and `sol`: the
  Legendre moments `aₗ` (`write_moments`) and the real spherical-harmonic coefficients
  `c_lm` of the observed incident flux, with their ratios to the mean term;
- `summary.txt` — totals (delivered, landed, missed; x-ray source, lost).
"""
function write_results(dir::AbstractString, t::Target, dep::Deposition; sol = nothing,
                       beams = nothing, observe = nothing, lmax::Integer = 8)
    mkpath(dir)
    files = String[]
    v(x) = ForwardDiff.value.(x)
    P = v(dep.element)
    fields = Pair[:laser_TW => P, :laser_Wcm2 => to_Wcm2.(P ./ t.elem_area)]
    if sol !== nothing
        append!(fields, [:source_Wcm2 => to_Wcm2.(v(sol.source)), :radiosity_Wcm2 => to_Wcm2.(v(sol.B)),
                         :incident_Wcm2 => to_Wcm2.(v(sol.H)), :albedo => sol.albedo,
                         :transmission => sol.transmission])
    end
    push!(files, write_element_table(joinpath(dir, "elements.csv"), t, fields...))
    cp = component_power(dep, t)
    bal = sol === nothing ? nothing : energy_balance(sol)
    path = joinpath(dir, "components.csv")
    open(path, "w") do io
        println(io, "component,kind,laser_TW,laser_outside_TW",
                bal === nothing ? "" : ",xray_absorbed_TW,xray_escaped_TW,xray_transmitted_TW")
        for c in t.components
            row = Any[c.name, c.kind, get(cp, c.name, 0.0), get(cp, Symbol(c.name, :_outside), 0.0)]
            if bal !== nothing
                push!(row, get(bal.absorbed, c.name, 0.0), get(bal.escaped, c.name, 0.0),
                      get(bal.transmitted, c.name, 0.0))
            end
            println(io, join((_csvfield(x) for x in row), ","))
        end
    end
    push!(files, path)
    if beams !== nothing
        length(beams) == length(dep.per_beam) ||
            throw(DimensionMismatch("one beam per deposited beam expected"))
        path = joinpath(dir, "beams.csv")
        open(path, "w") do io
            println(io, "name,quad,group,theta_deg,phi_deg,delivered_TW,landed_TW")
            for (b, landed) in zip(beams, dep.per_beam)
                θ, φ = port_angles(b)
                println(io, join(_csvfield.((b.name, b.quad, b.group, θ, φ, ForwardDiff.value(b.power),
                                             ForwardDiff.value(landed))), ","))
            end
        end
        push!(files, path)
    end
    if observe !== nothing && sol !== nothing
        obs = elements(t, (observe isa Symbol ? (observe,) : observe)...)
        H = v(sol.H)
        push!(files, write_moments(joinpath(dir, "legendre.csv"), legendre_moments(H, t, obs; lmax)))
        c = ylm_weights(t, obs; lmax) * H[obs]
        path = joinpath(dir, "ylm.csv")
        open(path, "w") do io
            println(io, "l,m,c_lm,c_lm_over_c00")
            for l in 0:lmax, m in -l:l
                x = c[mode_index(l, m)]
                println(io, l, ",", m, ",", _csvfield(x), ",", _csvfield(x / c[1]))
            end
        end
        push!(files, path)
    end
    path = joinpath(dir, "summary.txt")
    open(path, "w") do io
        delivered = ForwardDiff.value(dep.delivered)
        landed = sum(v(dep.per_beam))
        println(io, "elements = ", nelements(t))
        println(io, "laser_delivered_TW = ", delivered)
        println(io, "laser_landed_TW = ", landed)
        println(io, "laser_missed_TW = ", delivered - landed)
        println(io, "laser_outside_TW = ", sum(v(dep.backside)))
        if bal !== nothing
            println(io, "xray_source_TW = ", bal.source)
            println(io, "xray_lost_TW = ", bal.lost)
            println(io, "time_ns = ", sol.time)
        end
    end
    push!(files, path)
    return files
end
