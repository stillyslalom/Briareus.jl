# This file maps VisRad's per-element results onto a target built from its workspace and
# compares them per element with a Briareus radiation solution. The comparison covers:
# - laser deposition, from the power sources of the workspace;
# - incident and emitted flux and their temperatures, from VisRad's Surface Element
#   Properties export;
# - the modes of the observed component.

# Returns the column of a VisRad element table that matches `re`, or `nothing`.
_vr_find_col(tbl, re::Regex) = findfirst(h -> occursin(re, h), tbl.header)

# Checks whether a header names a flux in TW/cm², TW/mm² or W/cm².
_vr_is_flux(h) = occursin(r"w\s*/\s*(cm|mm)"i, h)

"""
    visrad_element_values(workspace, target, tables, column; index_base = nothing) -> Vector{Float64}

Maps one column of VisRad's per-element results onto the elements of `target`. The target
must be built on the VisRad elements, from `visrad_target(workspace)` or a selection of
its objects without `remesh`, `coarsen`, `fit` or `inside`. [`visrad_deposition`](@ref)
describes the element orders that result. `tables` maps object names to VisRad's Surface
Element Properties export of that object, which is saved with File | Save Data for
Object or as the `*.sed` file of a run. Each entry is a path, or a table as
[`read_visrad_table`](@ref) returns it. `column` is a
`Regex` matched against the table header (e.g. `r"flux\\s*out"i`, `r"rad\\s*t"i`).

Rows are matched to VisRad's surface elements by the "Surf Index" column of the table,
which may list a subset in any order. Without that column the rows must be the elements
of the object in order. VisRad does not document whether the index counts from 0, as in
the power-source tables of the workspace, or from 1 (*assumption*). With
`index_base = nothing`, the base is 0 when an index is 0 and 1 when the largest index
equals the element count; otherwise an `ArgumentError` is thrown. Pass `index_base = 0`
or `1` to fix it. Fluxes, with headers in TW/cm², TW/mm² or W/cm², are converted to the
Briareus unit TW/mm², the unit of `H` and `B` in a `RadiationSolution`. Other columns,
such as temperatures in eV, are returned as they are. A two-sided component takes the
values on its front elements. Elements without a value are `NaN`. An object named in
`tables` that is not in the workspace throws an error, and one without a component in
`target` is skipped with a warning.
"""
function visrad_element_values(ws::VisRadWorkspace, t::Target, tables::AbstractDict, column::Regex;
                               index_base::Union{Nothing,Integer} = nothing)
    out = fill(NaN, nelements(t))
    objs = Dict(o.name => o for o in ws.objects)
    for (name, tb) in tables
        name = String(name)
        o = get(objs, name, nothing)
        o === nothing && throw(ArgumentError("element table for \"$name\": the workspace has no such object"))
        tbl = tb isa AbstractString ? read_visrad_table(tb) : tb
        j = _vr_find_col(tbl, column)
        j === nothing && throw(ArgumentError("element table for $name has no column matching $column (columns: $(tbl.header))"))
        ci = findfirst(c -> String(c.name) == name, t.components)
        if ci === nothing
            @warn "VisRad object $name has an element table but no component in the target; not mapped"
            continue
        end
        n = o.kind === :mesh ? npatches(o.mesh) : npatches(t.components[ci].mesh)
        es = _vrw_row_elements(o, t, ci, n)
        idx = _vr_element_rows(tbl, n, index_base, name)
        h = tbl.header[j]
        f = _vr_is_flux(h) ? _vr_flux_scale(h) / 100 : 1.0             # to TW/mm²
        for (r, k) in enumerate(idx)
            e = es[k]
            e == 0 || (out[e] = tbl.data[r, j] * f)
        end
    end
    return out
end

# Returns the VisRad element (1-based) of each row of an element table of an object with
# `n` elements.
function _vr_element_rows(tbl, n, base, name)
    ji = _vr_find_col(tbl, r"surf\w*\s*ind|^index$"i)
    nr = size(tbl.data, 1)
    if ji === nothing
        nr == n || throw(ArgumentError("element table for $name has $nr rows and no \"Surf Index\" column; $name has $n elements"))
        return collect(1:n)
    end
    raw = tbl.data[:, ji]
    all(isinteger, raw) || throw(ArgumentError("element table for $name: the surface indices are not integers"))
    ix = Int.(raw)
    if base === nothing
        base = any(==(0), ix) ? 0 : maximum(ix) == n ? 1 :
               throw(ArgumentError("element table for $name: cannot tell whether its surface indices count from 0 or 1; pass index_base"))
    end
    k = ix .- base .+ 1
    all(i -> 1 <= i <= n, k) || throw(ArgumentError("element table for $name: surface indices outside 0:$(n - 1) (base $base)"))
    allunique(k) || throw(ArgumentError("element table for $name lists a surface element twice"))
    return k
end

# Returns the area-weighted agreement of the VisRad values `v` with the Briareus values
# `b` on elements with area `a`.
function _vr_element_stats(object, quantity, unit, b, v, a)
    A = sum(a)
    mb, mv = sum(a .* b) / A, sum(a .* v) / A
    d = v .- b
    rms = sqrt(sum(a .* d .^ 2) / A)
    cov = sum(a .* (b .- mb) .* (v .- mv))
    sb, sv = sqrt(sum(a .* (b .- mb) .^ 2)), sqrt(sum(a .* (v .- mv) .^ 2))
    return (; object, quantity, unit, n = length(b), briareus = mb, visrad = mv, rms = rms / abs(mv),
            max = maximum(abs, d) / abs(mv), correlation = sb > 0 && sv > 0 ? cov / (sb * sv) : NaN)
end

const _VR_ELEMENT_QUANTITIES = (
    (r"flux\s*in"i, "Flux In", "TW/cm2", :H), (r"flux\s*out"i, "Flux Out", "TW/cm2", :B),
    (r"rad\.?\s*t"i, "Rad T", "eV", :TH), (r"emis\.?\s*t"i, "Emis T", "eV", :TB))

"""
    compare_visrad(workspace, sol; elements = Dict(), sources = visrad_power_sources(workspace; stale = :nothing),
                   observe = :capsule, lmax = 4, center = nothing, index_base = nothing, io = nothing)
        -> (; rows, elements)

Compares a Briareus radiation solution `sol` with VisRad's results for the same
workspace, per element. `sol` is solved on a target with the VisRad elements, from
`visrad_target(workspace)` or a selection of its objects (see
[`visrad_deposition`](@ref)). The solve can start from the VisRad deposition, as in
`solve(RadiationModel(t), visrad_deposition(workspace, t))`, to isolate the radiation
solve, or from the Briareus deposition. Three comparisons are made:

- **Laser deposition** comes from the power sources of the workspace, `sources`, and
  `nothing` skips it. Per object, the rows give the absorbed power ("Laser Pwr Dep", TW)
  as Σ source × area against `sol.laser`. They also give the per-element agreement of the
  absorbed flux (TW/cm²) on the elements that either side lights.
- **`elements`** holds VisRad's Surface Element Properties export per object, as a
  `Dict` from object names to paths or [`read_visrad_table`](@ref) tables, mapped as in
  [`visrad_element_values`](@ref). Each column of the table among "Flux In" (`sol.H`),
  "Flux Out" (`sol.B`), "Rad T" and "Emis T" is compared per element on the elements
  that the table lists; "Rad T" and "Emis T" are the radiation temperatures of `H` and
  `B`. Per object, rows are added for the incident and emitted powers ("Incident Rad
  Pwr", "Emitted Rad Pwr"), Σ flux × area in TW. A row is also added for "Rad T", the
  radiation temperature (eV) of the area-mean incident flux. Each of these rows uses the
  elements that both sides have.
- **`observe`** names the component whose incident flux is analyzed in spherical
  harmonics, such as the capsule. Both sides use the same exact solid-angle weights
  ([`ylm_weights`](@ref) to `lmax`, about `center`). `center` defaults to the object
  position, the origin of the VisRad angles for that object; for a CAD mesh it is the
  bounding-box center. The rows are "P2/P0" and so on, the m = 0 coefficient as a
  Legendre ratio `aₗ/a₀`. Rows "mode amplitude l=2" and so on give
  `√((2l+1) Σₘ cₗₘ²) / c₀₀` for the orthonormal coefficients `cₗₘ`, which is `|aₗ/a₀|`
  for an axisymmetric field. All modes for l = 1, …, `lmax` are in percent. The analysis
  needs the "Flux In" of that object over all its elements. A table that VisRad fits
  itself (Analysis | Statistics) fits element centroids by least squares
  (`visrad_ylm_fit`). The weights in `compare_visrad` integrate the solid angle of each
  element exactly, so both sides are analyzed alike.

Returns `rows` and `elements`. `rows` holds
`(; quantity, unit, briareus, visrad, difference, relative)` entries, as the kit form of
`compare_visrad` gives them, with modes in percent. `elements` holds one
`(; object, quantity, unit, n, briareus, visrad, rms, max, correlation)` per object and
quantity. Its fields are the area-weighted means of both sides, the area-weighted rms
and the largest absolute difference as fractions of VisRad's mean, and the area-weighted
correlation.
Both are printed to `io` when given.
"""
function compare_visrad(ws::VisRadWorkspace, sol::RadiationSolution; elements::AbstractDict = Dict{String,Any}(),
                        sources = visrad_power_sources(ws; stale = :nothing), observe = :capsule,
                        lmax::Integer = 4, center = nothing, index_base::Union{Nothing,Integer} = nothing, io = nothing)
    t = sol.model.target
    v = ForwardDiff.value
    A = t.elem_area
    rows = NamedTuple{(:quantity, :unit, :briareus, :visrad, :difference, :relative),
                      Tuple{String,String,Float64,Float64,Float64,Float64}}[]
    add(q, u, x, y) = push!(rows, (quantity = q, unit = u, briareus = x, visrad = y,
                                   difference = y - x, relative = x == 0 ? NaN : (y - x) / abs(x)))
    stats = NamedTuple[]
    comp(name) = findfirst(c -> String(c.name) == name, t.components)
    if sources !== nothing
        PV = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do   # objects outside `t` are skipped
            visrad_deposition(ws, t; sources)
        end
        for name in sort!(collect(keys(sources.tables)))
            ci = comp(name)
            ci === nothing && continue
            es = t.comp_elems[ci]
            add("Laser Pwr Dep [$name]", "TW", sum(e -> v(sol.laser[e]), es; init = 0.0), sum(PV[es]; init = 0.0))
            lit = [e for e in es if PV[e] != 0 || v(sol.laser[e]) != 0]
            isempty(lit) || push!(stats, _vr_element_stats(name, "Laser", "TW/cm2", [100v(sol.laser[e]) / A[e] for e in lit],
                                                           [100PV[e] / A[e] for e in lit], A[lit]))
        end
    end
    H = v.(sol.H)
    B = v.(sol.B)
    field = Dict(:H => H, :B => B, :TH => radiation_temperature.(H), :TB => radiation_temperature.(B))
    for (name, tb) in sort!([String(k) => x for (k, x) in elements]; by = first)
        tbl = tb isa AbstractString ? read_visrad_table(tb) : tb
        tabs = Dict(name => tbl)
        for (re, label, unit, k) in _VR_ELEMENT_QUANTITIES
            _vr_find_col(tbl, re) === nothing && continue
            V = visrad_element_values(ws, t, tabs, re; index_base)
            es = findall(!isnan, V)
            isempty(es) && continue
            sc = k in (:H, :B) ? 100.0 : 1.0                                  # TW/mm² → TW/cm²
            push!(stats, _vr_element_stats(name, label, unit, sc .* field[k][es], sc .* V[es], A[es]))
            if k === :H
                add("Incident Rad Pwr [$name]", "TW", sum(H[es] .* A[es]), sum(V[es] .* A[es]))
                add("Rad T [$name]", "eV", radiation_temperature(sum(H[es] .* A[es]) / sum(A[es])),
                    radiation_temperature(sum(V[es] .* A[es]) / sum(A[es])))
                if String(observe) == name
                    ci = comp(name)
                    oes = t.comp_elems[ci]
                    length(es) == length(oes) ||
                        throw(ArgumentError("observed $name: the element table lists $(length(es)) of its $(length(oes)) elements"))
                    ctr = center !== nothing ? center : only(o.position for o in ws.objects if o.name == name)
                    Wy = ylm_weights(t, oes; lmax, center = ctr)
                    cb, cv = Wy * H[oes], Wy * V[oes]
                    for l in 1:lmax
                        add("P$l/P0 [$name]", "%", 100sqrt(2l + 1) * cb[mode_index(l, 0)] / cb[1],
                            100sqrt(2l + 1) * cv[mode_index(l, 0)] / cv[1])
                    end
                    amp(c, l) = 100sqrt((2l + 1) * sum(abs2, (c[mode_index(l, m)] for m in -l:l))) / c[1]
                    for l in 1:lmax
                        add("mode amplitude l=$l [$name]", "%", amp(cb, l), amp(cv, l))
                    end
                end
            elseif k === :B
                add("Emitted Rad Pwr [$name]", "TW", sum(B[es] .* A[es]), sum(V[es] .* A[es]))
            end
        end
    end
    stats = [s for s in stats]
    if io !== nothing
        @printf(io, "%-36s %12s %12s %12s %9s\n", "quantity", "Briareus", "VisRad", "V - B", "rel")
        for r in rows
            @printf(io, "%-36s %12.5g %12.5g %12.4g %8.2f%% %s\n", r.quantity, r.briareus, r.visrad,
                    r.difference, 100r.relative, r.unit)
        end
        println(io)
        @printf(io, "%-20s %-9s %6s %12s %12s %8s %8s %8s\n", "object", "quantity", "n", "Briareus", "VisRad", "rms", "max", "corr")
        for s in stats
            @printf(io, "%-20s %-9s %6d %12.5g %12.5g %7.2f%% %7.2f%% %8.5f %s\n", s.object, s.quantity, s.n,
                    s.briareus, s.visrad, 100s.rms, 100s.max, s.correlation, s.unit)
        end
    end
    return (; rows, elements = stats)
end

const _VR_RESIDUAL_COLUMNS = Dict(:flux_in => (r"flux\s*in"i, :H), :flux_out => (r"flux\s*out"i, :B),
                                  :rad_t => (r"rad\.?\s*t"i, :TH), :emis_t => (r"emis\.?\s*t"i, :TB))

"""
    visrad_residuals(workspace, sol; quantity = :laser, elements = Dict(),
                     sources = visrad_power_sources(workspace; stale = :nothing), relative = false,
                     index_base = nothing) -> (; residual, briareus, visrad, unit)
    visrad_residuals(workspace, target, laser; sources = …, relative = false, resample = false)

Element residuals `Briareus − VisRad` of one quantity, on the elements of the target of
`sol`, for coloring the mesh by the disagreement. [`compare_visrad`](@ref) gives the
summary rows of the same comparison, and the target must be built on the VisRad elements
as described there.

- `quantity = :laser` is the absorbed laser flux (TW/mm²), `sol.laser` over the element
  area against the power sources that the workspace saves (`sources`).
- `:flux_in`, `:flux_out`, `:rad_t` and `:emis_t` are the incident and emitted x-ray
  flux (TW/mm²) and their radiation temperatures (eV), against VisRad's Surface Element
  Properties export per object in `elements`, mapped as in
  [`visrad_element_values`](@ref).

The second form takes the absorbed laser power per element (TW), such as the `element`
field of a deposition, where no radiation solution exists. With `resample = true` the
target may be meshed differently from VisRad, as after `remesh`: the VisRad value of a
front element is that of the VisRad element nearest to its centroid on the object named
by the part of the element, by the distance to the triangles of the VisRad elements. The
comparison is then exact where the centroid lies on a VisRad element, and it assumes that
the object stands where the workspace places it. An object whose VisRad element order is
not established, such as a hohlraum, is left out and named in the field `skipped` of the
result. The field `objects` gives the summary of each compared object, as the laser rows
of [`compare_visrad`](@ref) give it: `object`, the absorbed power (TW) of both sides,
`power_briareus` on the front elements of the part and `power_visrad` on the elements of
the object in the workspace (Σ source × area, independent of the mesh), and the
agreement of the absorbed flux on the front elements (`n`, the area-weighted means
`briareus` and `visrad` in TW/cm², and `rms`, `max` and `correlation` as in
`compare_visrad`, with VisRad resampled).

`residual`, `briareus` and `visrad` have one entry per element. An element without a
VisRad value is `NaN` in `visrad` and `residual`. With `relative = true` the residual is
divided by the largest absolute VisRad value of the quantity, so the entries are
fractions of the peak.
"""
function visrad_residuals(ws::VisRadWorkspace, sol::RadiationSolution; quantity::Symbol = :laser,
                          elements::AbstractDict = Dict{String,Any}(),
                          sources = visrad_power_sources(ws; stale = :nothing), relative::Bool = false,
                          index_base::Union{Nothing,Integer} = nothing)
    t = sol.model.target
    v = ForwardDiff.value
    quantity === :laser && return visrad_residuals(ws, t, sol.laser; sources, relative)
    begin
        haskey(_VR_RESIDUAL_COLUMNS, quantity) ||
            throw(ArgumentError("quantity must be :laser, :flux_in, :flux_out, :rad_t or :emis_t, got :$quantity"))
        re, k = _VR_RESIDUAL_COLUMNS[quantity]
        V = fill(NaN, nelements(t))
        for (name, tb) in elements
            tbl = tb isa AbstractString ? read_visrad_table(tb) : tb
            _vr_find_col(tbl, re) === nothing && continue
            Vn = visrad_element_values(ws, t, Dict(String(name) => tbl), re; index_base)
            for e in eachindex(V)
                isnan(Vn[e]) || (V[e] = Vn[e])
            end
        end
        H = k in (:H, :TH) ? v.(sol.H) : v.(sol.B)
        B = k in (:TH, :TB) ? radiation_temperature.(H) : H
        unit = k in (:TH, :TB) ? "eV" : "TW/mm2"
    end
    return _vr_residual(B, V, unit, relative)
end

function visrad_residuals(ws::VisRadWorkspace, t::Target, laser::AbstractVector;
                          sources = visrad_power_sources(ws; stale = :nothing), relative::Bool = false,
                          resample::Bool = false)
    sources === nothing && throw(ArgumentError("the laser residual needs the power sources of the workspace"))
    length(laser) == nelements(t) || throw(ArgumentError("the laser power has $(length(laser)) entries for $(nelements(t)) elements"))
    v = ForwardDiff.value
    if resample
        B = [v(laser[e]) / t.elem_area[e] for e in 1:nelements(t)]
        V = fill(NaN, nelements(t))
        skipped = String[]
        objects = NamedTuple[]
        for name in sort!(collect(String, keys(sources.tables)))
            es = [e for e in 1:nelements(t) if t.components[t.elem_comp[e]].part === Symbol(name) && t.elem_side[e] == 1]
            isempty(es) && continue
            tv, PV = try
                Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                    tv = visrad_target(ws; select = o -> o.name == name)
                    tv, visrad_deposition(ws, tv; sources)
                end
            catch err
                err isa ArgumentError || rethrow()
                push!(skipped, name)
                continue
            end
            tris = [k for k in eachindex(tv.tri_front) if tv.tri_front[k] > 0]
            for e in es
                c = t.elem_centroid[e]
                k = argmin(k -> _point_triangle_distance(c, tv.tri_v0[k], tv.tri_v0[k] + tv.tri_e1[k], tv.tri_v0[k] + tv.tri_e2[k]), tris)
                ve = tv.tri_front[k]
                V[e] = PV[ve] / tv.elem_area[ve]
            end
            st = _vr_element_stats(name, "Laser", "TW/cm2", 100 .* B[es], 100 .* V[es], t.elem_area[es])
            push!(objects, (; object = name, power_briareus = sum(e -> v(laser[e]), es), power_visrad = sum(PV),
                            st.n, st.briareus, st.visrad, st.rms, st.max, st.correlation))
        end
        return (; _vr_residual(B, V, "TW/mm2", relative)..., skipped, objects = [o for o in objects])
    end
    names = Set(keys(sources.tables))
    PV = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        visrad_deposition(ws, t; sources)
    end
    B = [v(laser[e]) / t.elem_area[e] for e in eachindex(PV)]
    V = [String(t.components[t.elem_comp[e]].name) in names ? PV[e] / t.elem_area[e] : NaN for e in eachindex(PV)]
    return _vr_residual(B, V, "TW/mm2", relative)
end

function _vr_residual(B, V, unit, relative::Bool)
    residual = B .- V
    if relative
        peak = maximum(abs, (x for x in V if !isnan(x)); init = 0.0)
        peak > 0 && (residual ./= peak)
        unit = "of peak"
    end
    return (; residual, briareus = B, visrad = V, unit)
end
