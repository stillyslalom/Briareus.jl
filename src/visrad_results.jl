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
