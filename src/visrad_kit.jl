# This file reads VisRad result tables, fits their Legendre and Y_lm modes, writes the
# benchmark kit for rebuilding a Briareus case in VisRad and compares the two results.

# ---------------------------------------------------------------------------------------
# Result tables
# ---------------------------------------------------------------------------------------

# Unit-less column names of the VisRad element table. VR/results/images/surf_elem_table.jpg
# shows "Surf Index", and the others are *assumed* spellings of the listed quantities.
const _VR_UNITLESS = ("Surf Index", "Albedo", "Laser Reflectivity", "X-ray Conv Eff", "Object")

# Splits a whitespace-separated header into `n` names, on runs of 2+ spaces, or else by
# attaching each "(unit)" token to the words before it. Known unit-less names stand
# alone.
function _vr_header_names(line::AbstractString, n::Integer)
    s = strip(line)
    parts = split(s, r"\s{2,}")
    length(parts) == n && return String.(strip.(parts))
    toks = split(s)
    groups = String[]
    cur = String[]
    i = 1
    while i <= length(toks)
        known = findfirst(k -> (m = length(split(k)); i + m - 1 <= length(toks) &&
                                lowercase(join(toks[i:i+m-1], " ")) == lowercase(k)), _VR_UNITLESS)
        if isempty(cur) && known !== nothing
            m = length(split(_VR_UNITLESS[known]))
            push!(groups, join(toks[i:i+m-1], " "))
            i += m
            continue
        end
        t = toks[i]
        push!(cur, t)
        if endswith(t, ")")
            push!(groups, join(cur, " ")); empty!(cur)
        end
        i += 1
    end
    isempty(cur) || push!(groups, join(cur, " "))
    length(groups) == n && return groups
    length(toks) == n && return String.(toks)
    return String[]
end

"""
    read_visrad_table(path; header = nothing) -> (; header, labels, data)

Read a multi-column ASCII table exported by VisRad. The table is the Surface Element
Properties table (File | Save Data for Object / All Objects / Selected Rows, and the
`*.sed` files of a run; `VR/results/display_surf_elem_table.html`,
`VR/simulation_times/config_sim_params.html`) or a Power Summary (`*.psd`). The manual
shows the columns, for example "Z (cm)", "Azim. ang (deg)", "Radius (cm)", "Surf Index",
"Emis T (eV)", "Rad T (eV)", "Flux In (TW/cm2)" and "Flux Out (TW/cm2)"
(`VR/results/images/surf_elem_table.jpg`), but not the file layout. The reader is
therefore tolerant (*assumption* throughout):

- the table may be comma-, tab- or whitespace-separated;
- data rows are the most common shape of lines ending in numbers, and any leading
  non-numeric text is kept as the label of the row (an object name);
- the header is the last non-data line before the first data row. It is split on the
  delimiter, on runs of 2+ spaces, or by attaching "(unit)" tokens to the words before
  them.

When that fails, `header` (a vector of names for the numeric columns) sets the names.

Returns the numeric column names, row labels (empty strings when none) and an
`nrows × ncols` matrix.
"""
function read_visrad_table(path::AbstractString; header = nothing)
    lines = [l for l in _vr_lines(path) if !isempty(strip(l))]
    delim = any(l -> occursin('\t', l), lines) ? '\t' :
            count(l -> occursin(',', l), lines) > length(lines) ÷ 2 ? ',' : ' '
    function toks(l)
        delim == ' ' ? String.(split(strip(l))) : _vr_split(l, delim)
    end
    function shape(l)
        f = toks(l)
        k = 0
        for x in Iterators.reverse(f)
            _vr_num(x) === nothing && break
            k += 1
        end
        return k, f
    end
    shapes = [shape(l) for l in lines]
    ks = [k for (k, _) in shapes if k > 0]
    isempty(ks) && throw(ArgumentError("$path has no numeric rows"))
    kmode = argmax(k -> count(==(k), ks), unique(ks))
    first_row = findfirst(s -> s[1] == kmode, shapes)
    labels = String[]
    data = Vector{Float64}[]
    for (k, f) in shapes[first_row:end]
        k == kmode || continue
        push!(labels, join(f[1:end-k], " "))
        push!(data, [_vr_num(x) for x in f[end-k+1:end]])
    end
    names = if header !== nothing
        String.(collect(header))
    else
        hl = findlast(i -> shapes[i][1] != kmode, 1:first_row-1)
        hs = hl === nothing ? String[] :
             delim == ' ' ? _vr_header_names(lines[hl], kmode + any(!isempty, labels)) :
             _vr_split(lines[hl], delim)
        if length(hs) == kmode + 1 && any(!isempty, labels)
            hs = hs[2:end]
        elseif length(hs) == kmode - 1
            hs = ["index"; hs]
        end
        length(hs) == kmode ? hs : ["col$j" for j in 1:kmode]
    end
    length(names) == kmode || throw(ArgumentError("header has $(length(names)) names for $kmode columns"))
    M = isempty(data) ? zeros(0, kmode) : permutedims(reduce(hcat, data))
    return (; header = names, labels, data = M)
end

# Flux unit scale to TW/cm² from a header.
_vr_flux_scale(h) = occursin(r"mm\s*\^?2|mm²"i, h) ? 100.0 : occursin(r"w/cm"i, h) && !occursin(r"tw"i, h) ? 1e-12 : 1.0

"""
    visrad_legendre_fit(theta, values; terms = 0:4, weights = nothing) -> Vector

Least-squares fit `f(θ) = Σₗ Aₗ Pₗ(cos θ)` over the given `terms`, as in the Surface
Element Statistics of VisRad, with coefficients up to A4 and each term selected by a
check box (`VR/results/display_surf_elem_table.html`). `theta` is in degrees. The result
`A` is indexed `l+1`, with zeros for unused terms, so `A[3]/A[1]` is P2/P0. The fit is
unweighted by default, since the VisRad weighting is not documented (*assumption*);
element areas passed as `weights` weight it. The same fit applied to Briareus and VisRad
element values gives a like-for-like comparison.
"""
function visrad_legendre_fit(θ::AbstractVector, values::AbstractVector; terms = 0:4,
                             weights = nothing)
    ls = collect(terms)
    X = [legendre(l, cosd(t)) for t in θ, l in ls]
    w = weights === nothing ? ones(length(values)) : sqrt.(collect(Float64, weights))
    c = (w .* X) \ (w .* collect(Float64, values))
    A = zeros(maximum(ls) + 1)
    for (j, l) in enumerate(ls)
        A[l+1] = c[j]
    end
    return A
end

# Associated Legendre P_l^m(x) without the Condon–Shortley phase, P[l+1, m+1].
function _vr_plm(lmax, x)
    s = sqrt(max(0.0, 1 - x^2))
    P = zeros(lmax + 1, lmax + 1)
    P[1, 1] = 1.0
    for m in 1:lmax
        P[m+1, m+1] = (2m - 1) * s * P[m, m]
    end
    for m in 0:lmax-1
        P[m+2, m+1] = x * (2m + 1) * P[m+1, m+1]
    end
    for m in 0:lmax, l in m+2:lmax
        P[l+1, m+1] = ((2l - 1) * x * P[l, m+1] - (l + m - 1) * P[l-1, m+1]) / (l - m)
    end
    return P
end

"""
    visrad_ylm_fit(theta, phi, values; lmax = 4, weights = nothing) -> Vector

Least-squares fit to the real spherical harmonics of VisRad
(`VR/results/display_surf_elem_table.html`). These are `Y_LM = N_LM P_L^M(cos θ)` for
M = 0 and `√2 N_LM P_L^|M|(cos θ)` × `cos Mφ` (M > 0) or `sin |M|φ` (M < 0), with
`N_LM = √((2L+1)(L−M)!/(L+M)!)`. Then `Y_00 = 1`, and `A_00` is the mean; the example in
the manual gives a mean of 9.736 and `A_00 = 9.739`. `A_L0 = A_L/√(2L+1)` for the
Legendre `A_L`, and the same example gives 1.143·√5 = 2.556 against `A_2 = 2.563`.
Angles are in degrees, and the coefficients are in the flat `mode_index` order of
Briareus. They are the orthonormal `ylm_moments` of Briareus times `√(4π)`. The manual
does not state whether its `P_L^M` includes the Condon–Shortley phase. No phase is
assumed, so the signs of odd-M coefficients may flip; `√(Σ_M A_LM²)/A_00` is
unaffected.
"""
function visrad_ylm_fit(θ::AbstractVector, φ::AbstractVector, values::AbstractVector;
                        lmax::Integer = 4, weights = nothing)
    0 <= lmax <= 10 || throw(ArgumentError("lmax must be in 0:10"))
    n = length(values)
    X = zeros(n, (lmax + 1)^2)
    for i in 1:n
        P = _vr_plm(lmax, cosd(θ[i]))
        for l in 0:lmax, m in -l:l
            am = abs(m)
            N = sqrt((2l + 1) * factorial(l - am) / factorial(l + am))
            X[i, mode_index(l, m)] = m == 0 ? N * P[l+1, 1] :
                                     m > 0 ? sqrt(2) * N * P[l+1, am+1] * cosd(m * φ[i]) :
                                     sqrt(2) * N * P[l+1, am+1] * sind(am * φ[i])
        end
    end
    w = weights === nothing ? ones(n) : sqrt.(collect(Float64, weights))
    return (w .* X) \ (w .* collect(Float64, values))
end

_vr_mode_amp(A, l) = sqrt(sum(abs2, (A[mode_index(l, m)] for m in -l:l))) / abs(A[1])

# Fit-based quantities shared by the kit and the comparison.
function _vr_fit_quantities(θ, φ, flux, w; terms, lmax)
    q = Pair{String,Tuple{Float64,String}}[]
    A = visrad_legendre_fit(θ, flux; terms, weights = w)
    push!(q, "Flux In mean" => (sum(flux) / length(flux), "TW/cm2"))
    push!(q, "Legendre A0" => (A[1], "TW/cm2"))
    for l in terms
        l == 0 && continue
        push!(q, "Legendre A$l/A0" => (A[l+1] / A[1], ""))
    end
    Y = visrad_ylm_fit(θ, φ, flux; lmax, weights = w)
    for l in 1:lmax
        push!(q, "Ylm amplitude L=$l" => (_vr_mode_amp(Y, l), ""))
    end
    return q
end

# ---------------------------------------------------------------------------------------
# Benchmark kit
# ---------------------------------------------------------------------------------------

const _VR_POWER_COLUMNS = ("Surface Area", "Incident Rad Pwr", "Emitted Rad Pwr", "Pwr Src Term",
                           "Laser Pwr Dep")

# Returns the conversion of a component as one number. That is the value itself, or for a
# `Conversion` model the laser-weighted mean of the solved source. Without laser, it is
# the model at the time of the solution.
_mean_conversion(η::Real, sol, t, es) = η
function _mean_conversion(m::Conversion, sol, t, es)
    v = ForwardDiff.value
    P = sum(v(sol.laser[e]) for e in es; init = 0.0)
    P > 0 && return sum(v(sol.source[e]) * t.elem_area[e] for e in es) / P
    return _intensity_dependent(m) ? NaN : conversion_efficiency(m, sol.time)
end

function _vr_describe(io, t::Target, sol)
    for (ci, c) in enumerate(t.components)
        vs = c.mesh.vertices
        es = t.comp_elems[ci]
        lo = reduce((a, b) -> min.(a, b), vs)
        hi = reduce((a, b) -> max.(a, b), vs)
        A = sum(t.elem_area[es]; init = 0.0) / (c.two_sided ? 2 : 1)
        α = isempty(es) ? 0.0 : sum(sol.albedo[es] .* t.elem_area[es]) / sum(t.elem_area[es])
        lip = c.kind === :aperture && _transmits(c.material.albedo)
        kind = lip ? "closing aperture lip (transparent to laser light; re-emits the x-ray albedo below, " *
                     "the rest of its incident flux escapes; VisRad: an aperture while open, else a " *
                     "laser-transparent surface of that albedo)" :
               c.kind === :aperture ? "aperture (transparent; its incident flux is the escaping flux)" :
               c.two_sided ? "opaque, two-sided" : "opaque, radiates from its front side"
        @printf(io, "- %s: %s, %d elements, area %.6g mm² (%.6g cm²)\n", c.name, kind,
                length(es), A, A / 100)
        (c.kind === :aperture && !lip) ||
            @printf(io, "    material %s: x-ray albedo %.4g (area mean as solved), laser→x-ray conversion %.4g\n",
                    c.material.name, α, _mean_conversion(c.material.conversion, sol, t, es))
        s = c.surface
        p3(x) = @sprintf("(%.6g, %.6g, %.6g)", x[1], x[2], x[3])
        if s isa SphereSurface
            @printf(io, "    sphere: radius %.6g mm (%.6g µm), center %s mm\n",
                    s.radius, 1000 * s.radius, p3(s.center))
        elseif s isa CylinderSurface
            zs = [dot(v - s.origin, s.axis) for v in vs]
            @printf(io, "    cylinder: radius %.6g mm (%.6g µm), axis %s through %s mm, from %.6g to %.6g mm along it\n",
                    s.radius, 1000 * s.radius, p3(s.axis), p3(s.origin), minimum(zs), maximum(zs))
        elseif s isa ConeSurface
            zs = [dot(v - s.apex, s.axis) for v in vs]
            @printf(io, "    cone: half-angle %.6g°, apex %s mm, axis %s, from %.6g to %.6g mm along it
",
                    s.half_angle, p3(s.apex), p3(s.axis), minimum(zs), maximum(zs))
        elseif s isa TorusSurface
            @printf(io, "    torus: major radius %.6g mm, minor radius %.6g mm, axis %s through %s mm\n",
                    s.major_radius, s.minor_radius, p3(s.axis), p3(s.center))
        elseif s isa DiskSurface
            rs = [norm((v - s.center) - dot(v - s.center, s.normal) * s.normal) for v in vs]
            @printf(io, "    planar disk/annulus: center %s mm, normal %s, radii %.6g to %.6g mm\n",
                    p3(s.center), p3(s.normal), minimum(rs), maximum(rs))
        end
        @printf(io, "    bounding box x [%.6g, %.6g], y [%.6g, %.6g], z [%.6g, %.6g] mm\n",
                lo[1], hi[1], lo[2], hi[2], lo[3], hi[3])
    end
end

"""
    write_visrad_kit(dir, target, beams; model = nothing, op = nothing, pulses = nothing, by = :group,
                     time = 0.0, observe = :capsule, terms = (0, 2, 4), lmax = 4,
                     rotation_sense = :beam, notes = "") -> NamedTuple

Write everything needed to rebuild a Briareus case in VisRad, together with the Briareus
results for comparison (`compare_visrad`), into directory `dir`:

- `pointing_sst.csv` — per-beam pointing in the NIF SST format (`write_visrad_sst`),
  imported with File | Import Laser Beam Data | Pointing Data, Custom CSV, units mm.
- `beam_summary.csv` — the full beam set in Laser Beam Summary layout
  (`write_visrad_summary`): port angles, pointing, focus, power, 1/e spot radii,
  exponents and rotations.
- `powers.csv` — per-beam time-dependent powers (`write_visrad_powers`), when `pulses` is
  given (keyed `by`, as in `write_visrad_powers`); beam powers are then evaluated at `time`.
- `briareus_beams.csv` — the canonical Briareus beam table (`write_beam_table`).
- `README.txt` — target geometry (component shapes and dimensions in mm, cm and µm,
  materials, albedos and conversion efficiencies), beam settings, and the VisRad steps.
- `briareus_results.csv` — `quantity,value,unit` rows: per-component "Laser Pwr Dep",
  "Pwr Src Term" (conversion × laser), "Incident Rad Pwr" (Σ H·A), "Emitted Rad Pwr"
  (Σ B·A), all TW, and "Surface Area" (cm²), named as the Power Summary columns of
  VisRad; the laser power that misses or strikes outer surfaces; and for `observe` the mean
  incident flux (TW/cm²), the VisRad-style Legendre fit (`visrad_legendre_fit` over
  `terms`) and per-L amplitudes of the VisRad spherical-harmonic fit (`visrad_ylm_fit`).
- `briareus_<observe>_elements.csv` — per-element polar and azimuthal angle about the
  center of the component (degrees), area (cm²), incident flux and absorbed laser flux
  (TW/cm²).

The deposition uses `op` (default `DepositionOperator(target, beams)`) and the radiation
`model` (default `RadiationModel(target)`). Returns `(; dir, files, results)` with
`results` a `Dict` of the quantities in `briareus_results.csv`.
"""
function write_visrad_kit(dir::AbstractString, target::Target, beams::AbstractVector{<:Beam};
                          model = nothing, op = nothing, pulses = nothing, by = :group, time = 0.0,
                          observe = :capsule, terms = (0, 2, 4), lmax::Integer = 4,
                          rotation_sense::Symbol = :beam, notes::AbstractString = "")
    mkpath(dir)
    files = Dict{Symbol,String}()
    bt = pulses === nothing ? beams : with_pulses(beams, pulses, time; by)
    files[:sst] = write_visrad_sst(joinpath(dir, "pointing_sst.csv"), bt)
    files[:summary] = write_visrad_summary(joinpath(dir, "beam_summary.csv"), bt; rotation_sense)
    pulses === nothing || (files[:powers] = write_visrad_powers(joinpath(dir, "powers.csv"), beams, pulses; by))
    files[:beams] = write_beam_table(joinpath(dir, "briareus_beams.csv"), bt)

    op = op === nothing ? DepositionOperator(target, bt) : op
    model = model === nothing ? RadiationModel(target) : model
    dep = deposit(op, bt)
    sol = solve(model, dep; time)
    t = target
    v = ForwardDiff.value
    results = Pair{String,Tuple{Float64,String}}[]
    for p in parts(t)                    # a part split for computation is treated as one
        cis = part_pieces(t, p)
        c = t.components[first(cis)]
        es = elements(t, p)
        n = string(p)
        A = sum(t.elem_area[es]; init = 0.0) / (c.two_sided ? 2 : 1)
        laser = sum(v(dep.element[e]) for e in es; init = 0.0)
        push!(results, "Surface Area [$n]" => (A / 100, "cm2"))
        push!(results, "Incident Rad Pwr [$n]" => (sum(v(sol.H[e]) * t.elem_area[e] for e in es; init = 0.0), "TW"))
        push!(results, "Emitted Rad Pwr [$n]" => (sum(v(sol.B[e]) * t.elem_area[e] for e in es; init = 0.0), "TW"))
        c.kind === :aperture && continue
        η = c.material.conversion
        src = η isa Real ? η * laser : sum(v(sol.source[e]) * t.elem_area[e] for e in es; init = 0.0)
        push!(results, "Pwr Src Term [$n]" => (src, "TW"))
        push!(results, "Laser Pwr Dep [$n]" => (laser, "TW"))
        bs = sum(ci -> v(dep.backside[ci]), cis)
        bs > 0 && push!(results, "Laser Pwr Dep [$(n)_outside]" => (bs, "TW"))
    end
    push!(results, "Laser delivered" => (v(dep.delivered), "TW"))
    push!(results, "Laser missed" => (v(dep.delivered - sum(dep.per_beam)), "TW"))

    obs_ok = any(c -> c.name === observe || c.part === observe, t.components)
    if obs_ok
        ci = first(part_pieces(t, observe))     # the largest piece of a split part comes first
        es = elements(t, observe)
        s = t.components[ci].surface
        ctr = s isa SphereSurface ? s.center :
              sum(t.elem_centroid[es] .* t.elem_area[es]) / sum(t.elem_area[es])
        θ = [acosd(clamp(normalize(t.elem_centroid[e] - ctr)[3], -1, 1)) for e in es]
        φ = [mod(atand((t.elem_centroid[e]-ctr)[2], (t.elem_centroid[e]-ctr)[1]), 360) for e in es]
        H = [100 * v(sol.H[e]) for e in es]
        L = [100 * v(dep.element[e]) / t.elem_area[e] for e in es]
        fpath = joinpath(dir, "briareus_$(observe)_elements.csv")
        open(fpath, "w") do io
            println(io, "element,theta_deg,phi_deg,area_cm2,flux_in_TWcm2,laser_TWcm2")
            for (k, e) in enumerate(es)
                println(io, join((e, _vr_fmt(θ[k]), _vr_fmt(φ[k]), _vr_fmt(t.elem_area[e] / 100),
                                  _vr_fmt(H[k]), _vr_fmt(L[k])), ","))
            end
        end
        files[:elements] = fpath
        for (k, val) in _vr_fit_quantities(θ, φ, H, nothing; terms, lmax)
            push!(results, "$k [$observe]" => val)
        end
        push!(results, "Tr mean [$observe]" =>
              (radiation_temperature(sum(v(sol.H[e]) * t.elem_area[e] for e in es) / sum(t.elem_area[es])), "eV"))
    end
    files[:results] = joinpath(dir, "briareus_results.csv")
    open(files[:results], "w") do io
        println(io, "quantity,value,unit")
        for (k, (x, u)) in results
            println(io, _csvfield(k), ",", _vr_fmt(x), ",", u)
        end
    end

    files[:readme] = joinpath(dir, "README.txt")
    open(files[:readme], "w") do io
        _vr_readme(io, t, bt, sol, dep, pulses, time, observe, obs_ok, terms, lmax,
                   rotation_sense, notes, files)
    end
    return (; dir, files, results = Dict(k => first(x) for (k, x) in results))
end

function _vr_readme(io, t, beams, sol, dep, pulses, time, observe, obs_ok, terms, lmax,
                    sense, notes, files)
    v = ForwardDiff.value
    println(io, "Briareus -> VisRad benchmark kit")
    println(io, "================================")
    println(io)
    println(io, "Units here: mm, ns, TW, degrees unless stated. VisRad works in cm, s, TW and")
    println(io, "TW/cm2 (VR/intro/conventions.html); spot radii in VisRad dialogs are in µm.")
    println(io, "Frame: target-chamber Cartesian about TCC, z up; for NIF +x toward Switchyard 1")
    println(io, "(φ = 0), +y at φ = 90°. VisRad's NIF frame agrees (Q15T at θ = 23.5°, φ = 78.75°")
    println(io, "in VR/intro/nif/sample_nif_3wrect_lasers.html).")
    isempty(notes) || (println(io); println(io, notes))
    println(io)
    println(io, "Target components")
    println(io, "-----------------")
    _vr_describe(io, t, sol)
    println(io)
    println(io, "VisRad note: for Cylindrical Hohlraum, Sphere, Cylinder and similar objects the")
    println(io, "grid size defines the OUTER wall and the inner wall lies inward by the Wall")
    println(io, "Thickness (VR/intro/revision_history.html, v21.0.0). Briareus dimensions above are")
    println(io, "the radiating surfaces: use zero wall thickness or add the thickness to the size.")
    println(io, "Name VisRad objects after the Briareus components so compare_visrad can match them")
    println(io, "(or pass `names`).")
    println(io)
    println(io, "Beams")
    println(io, "-----")
    on = count(b -> v(b.power) > 0, beams)
    @printf(io, "%d beams (%d on), total power %.6g TW at t = %.6g ns.\n", length(beams), on,
            sum(b -> v(b.power), beams), time)
    groups = unique(b.group for b in beams)
    for g in groups
        bs = [b for b in beams if b.group === g]
        b = first(bs)
        sp = b.spot
        @printf(io, "- group %s: %d beams, %.6g TW each (first), f/%.4g, defocus %.6g mm\n",
                g, length(bs), v(b.power), b.fnumber, v(b.defocus))
        _is_nif_quad(b.name) &&
            @printf(io, "    (NIF quads: each of the 4 VisRad beams of a quad carries %.6g TW)\n", v(b.power) / 4)
        rmaj, rmin = 1000 * max(v(sp.a), v(sp.b)), 1000 * min(v(sp.a), v(sp.b))
        hw = log(2)^(1 / sp.order)
        @printf(io, "    spot: super-Gaussian n = %.6g, 1/e radii %.6g × %.6g µm (major × minor), rotation %.6g° (VisRad minor-axis convention, sense :%s)\n",
                sp.order, rmaj, rmin, _vr_rotation(sp, sense, _vr_reference_angle(b)), sense)
        @printf(io, "    = radii %.6g × %.6g µm at half maximum (VisRad NIF CPPs are HWHM)\n",
                rmaj * hw, rmin * hw)
        sp.halo > 0 && @printf(io, "    halo: %.3g of the power at %.3g× the core radii (not representable in VisRad)\n",
                               sp.halo, sp.halo_scale)
    end
    println(io)
    println(io, "Steps in VisRad")
    println(io, "---------------")
    println(io, "1. New workspace for the facility; build the target components listed above.")
    println(io, "2. File | Import Laser Beam Data | Pointing Data -> Custom *.csv file ->")
    println(io, "   pointing_sst.csv. Columns map to Beam / Pointing: Beam Aim Point X, Y, Z (mm) /")
    println(io, "   Beam Defocus (mm); Units: mm. VisRad sets the reference frame to Target Chamber")
    println(io, "   and turns on exactly the beams in the file")
    println(io, "   (VR/power_sources/import_laser_beam_data_csv.html).")
    println(io, "3. Spatial profile of every beam: Custom Supergaussian, elliptical, spot size")
    println(io, "   defined at 36.8% (1/e), exponent, major/minor radii and rotation per group as")
    println(io, "   listed above (beam_summary.csv has them per beam). Check the rotation sense")
    println(io, "   (see research notes) on one elliptical beam with a spot overlay.")
    if pulses === nothing
        println(io, "4. Powers: constant, per beam as in beam_summary.csv (\"Power\", TW).")
    else
        println(io, "4. Powers: File | Import Laser Beam Data | Time-Dependent Powers -> powers.csv")
        @printf(io, "   (time in ns, TW per beam). Compare at t = %.6g ns.\n", time)
    end
    println(io, "   NIF: VisRad splits a beam into 3ω/2ω/1ω cones; turn the 2ω and 1ω cones off")
    println(io, "   (or set their power fraction to zero) so all power is 3ω like Briareus")
    println(io, "   (VR/power_sources/nif_beam_params.html).")
    println(io, "5. Materials: albedo and x-ray conversion efficiency per component as above.")
    println(io, "   Use many beamlets (e.g. 1000; Lasers | Beam Profile) to limit sampling noise.")
    println(io, "6. Run; then export View | Surface Elements Data -> File | Save Data for Object")
    obs_ok && println(io, "   for `$observe` (incident flux \"Flux In\" with polar/azimuthal angles), and")
    println(io, "   the Power Summary table values. Note the Legendre fit (terms $(collect(terms)))")
    println(io, "   and spherical-harmonic table from Analysis | Statistics for comparison.")
    println(io, "7. At home: compare_visrad(<this directory>; elements = <element file>,")
    println(io, "   power_summary = <file or Dict>) in Briareus.")
    println(io)
    println(io, "Expected differences")
    println(io, "--------------------")
    println(io, "- Away from focus VisRad widens the spot as a(z) = a0 [ (z/zA)^2 + 1 ]^(1/2) with a")
    println(io, "  z-dependent exponent (VR/power_sources/set_laser_profile_props.html); Briareus")
    println(io, "  convolves the focal spot with the f-cone. Defocused beams differ in footprint shape.")
    println(io, "- VisRad's beam axis runs from the port on the chamber wall through the aim point;")
    println(io, "  Briareus keeps it parallel to the port direction (≲ 1 µm at the focus).")
    println(io, "- Briareus integrates radiosity on the exact curved surfaces; VisRad uses flat facets")
    println(io, "  renormalized to the nominal area (VR/results/display_surf_elem_table.html).")
    println(io)
    println(io, "Briareus results (briareus_results.csv)")
    println(io, "---------------------------------------")
    @printf(io, "Laser delivered %.6g TW; absorbed per component:\n", v(dep.delivered))
    for (k, x) in component_power(dep, t)
        @printf(io, "  %-20s %.6g TW\n", k, x)
    end
    println(io)
    println(io, _VR_COMMENT)
end

"""
    compare_visrad(kit_dir; elements = nothing, power_summary = nothing, observe = :capsule,
                   names = Dict(), terms = (0, 2, 4), lmax = 4, io = nothing)
        -> Vector{NamedTuple}

Compare a kit written by `write_visrad_kit` with the VisRad results for the same case.

- `elements`: the VisRad Surface Element Properties export for the observed component
  (`read_visrad_table`). The polar and azimuthal columns are those whose headers contain
  "polar"/"theta" and "azim"/"phi", in degrees and in the component frame, which for an
  unrotated sphere is the chamber frame. They and "Flux In" (TW/cm²; a TW/mm² header is
  converted) are fitted exactly as the Briareus elements were (`visrad_legendre_fit`
  over `terms`, `visrad_ylm_fit` to `lmax`).
- `power_summary`: the VisRad Power Summary, as a file (`read_visrad_table`, with object
  names as row labels) or a `Dict(object => Dict(column => value))` typed from the
  dialog. The columns are as in VisRad ("Laser Pwr Dep", "Incident Rad Pwr", "Emitted
  Rad Pwr", "Pwr Src Term", "Surface Area"; TW and cm²). `names` maps VisRad object
  names to Briareus component names, e.g. `Dict("hohlraum_barrel" => "wall")`.

Returns rows `(; quantity, unit, briareus, visrad, difference, relative)` and prints
them to `io` when given.
"""
function compare_visrad(kit_dir::AbstractString; elements = nothing, power_summary = nothing,
                        observe = :capsule, names = Dict{String,String}(), terms = (0, 2, 4),
                        lmax::Integer = 4, io = nothing)
    hdr, rows = read_csv(joinpath(kit_dir, "briareus_results.csv"))
    bri = Dict(r[1] => (parse(Float64, r[2]), length(r) >= 3 ? r[3] : "") for r in rows)
    out = NamedTuple{(:quantity, :unit, :briareus, :visrad, :difference, :relative),
                     Tuple{String,String,Float64,Float64,Float64,Float64}}[]
    add(q, u, x, y) = push!(out, (quantity = q, unit = u, briareus = x, visrad = y,
                                  difference = y - x, relative = x == 0 ? NaN : (y - x) / abs(x)))
    if elements !== nothing
        tbl = elements isa AbstractString ? read_visrad_table(elements) : elements
        jθ = findfirst(h -> occursin(r"polar|theta|θ"i, h), tbl.header)
        jφ = findfirst(h -> occursin(r"azim|phi|φ"i, h), tbl.header)
        jH = findfirst(h -> occursin(r"flux\s*in"i, h), tbl.header)
        any(isnothing, (jθ, jφ, jH)) &&
            throw(ArgumentError("element table needs polar, azimuthal and Flux In columns; has $(tbl.header)"))
        H = tbl.data[:, jH] .* _vr_flux_scale(tbl.header[jH])
        for (k, (y, u)) in _vr_fit_quantities(tbl.data[:, jθ], tbl.data[:, jφ], H, nothing; terms, lmax)
            key = "$k [$observe]"
            haskey(bri, key) && add(key, u, first(bri[key]), y)
        end
    end
    if power_summary !== nothing
        ps = power_summary isa AbstractString ? _vr_power_summary(power_summary) : power_summary
        for (obj, cols) in ps
            comp = get(names, String(obj), String(obj))
            for (cname, y) in cols
                c = findfirst(p -> occursin(lowercase(p), lowercase(String(cname))), collect(_VR_POWER_COLUMNS))
                c === nothing && continue
                key = "$(_VR_POWER_COLUMNS[c]) [$comp]"
                haskey(bri, key) && add(key, last(bri[key]), first(bri[key]), Float64(y))
            end
        end
    end
    if io !== nothing
        @printf(io, "%-36s %12s %12s %12s %9s\n", "quantity", "Briareus", "VisRad", "V - B", "rel")
        for r in out
            @printf(io, "%-36s %12.5g %12.5g %12.4g %8.2f%% %s\n", r.quantity, r.briareus, r.visrad,
                    r.difference, 100r.relative, r.unit)
        end
    end
    return out
end

function _vr_power_summary(path)
    tbl = read_visrad_table(path)
    return Dict(tbl.labels[i] => Dict(tbl.header[j] => tbl.data[i, j] for j in eachindex(tbl.header))
                for i in eachindex(tbl.labels) if !isempty(tbl.labels[i]))
end
