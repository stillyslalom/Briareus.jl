# The difference of two designs: which parts moved, which beams were repointed, how the
# TaLIS margins changed and how the physics of each role changed. The page is written by
# html_report.jl, as the review packets are.

const _DiffPart = @NamedTuple{part::Symbol, status::Symbol, shift::Vec3, distance::Float64,
                              displacement::Float64, material::Bool}
const _DiffBeam = @NamedTuple{beam::String, status::Symbol, aim_shift::Vec3, distance::Float64,
                              power::NTuple{2,Float64}, defocus::NTuple{2,Float64}, changes::Vector{Symbol}}
const _DiffTaLIS = @NamedTuple{rule::Symbol, beam::String, component::Symbol, unit::Symbol,
                               status::NTuple{2,Symbol}, margin::NTuple{2,Float64}, delta::Float64,
                               change::Symbol}
const _DiffPhysics = @NamedTuple{role::String, quantity::String, unit::String, a::Float64, b::Float64,
                                 delta::Float64}

"""
    DesignDiff

The difference of two designs, from [`design_diff`](@ref). `names` labels the reference
and the changed design. Every other field is a vector of rows:

- `parts`: `(; part, status, shift, distance, displacement, material)`. `status` is
  `:same`, `:moved` for a rigid translation, `:changed` for any other change of the
  mesh, `:added` or `:removed`. `shift` (mm) is the translation of a moved part and the
  shift of the area centroid otherwise, `distance` its length, `displacement` the
  largest vertex displacement (`NaN` when the vertex counts differ), and `material`
  whether the material changed.
- `beams`: `(; beam, status, aim_shift, distance, power, defocus, changes)`. `status` is
  `:same`, `:changed`, `:added` or `:removed`, `power` and `defocus` hold the two values,
  and `changes` lists what differs among `:aim`, `:power`, `:defocus`, `:spot` and
  `:port`.
- `talis`: `(; rule, beam, component, unit, status, margin, delta, change)` per row of
  either report, with the two statuses and margins and `delta = margin[2] − margin[1]`.
  `change` is `:newly_failing`, `:fixed`, `:worse`, `:improved`, `:same`, `:added` or
  `:removed`.
- `physics`: `(; role, quantity, unit, a, b, delta)` with `delta = b − a`.
"""
struct DesignDiff
    names::NTuple{2,String}
    parts::Vector{_DiffPart}
    beams::Vector{_DiffBeam}
    talis::Vector{_DiffTaLIS}
    physics::Vector{_DiffPhysics}
end

function Base.show(io::IO, d::DesignDiff)
    n(v, f) = count(f, v)
    print(io, "DesignDiff(", d.names[1], " → ", d.names[2], "; parts: ", n(d.parts, x -> x.status !== :same),
          " of ", length(d.parts), " differ; beams: ", n(d.beams, x -> x.status !== :same), " of ", length(d.beams),
          " differ; TaLIS: ", n(d.talis, x -> x.change === :newly_failing), " newly failing, ",
          n(d.talis, x -> x.change === :fixed), " fixed; ", length(d.physics), " physics rows)")
end

# Vertices and the area centroid of every piece of a part, in the order of the target.
function _diff_part_geometry(t::Target, part::Symbol)
    V = Vec3[]
    c, A = zero(Vec3), 0.0
    for comp in t.components
        comp.part === part || continue
        m = comp.mesh
        append!(V, m.vertices)
        for (i, j, k) in m.triangles
            a = norm(cross(m.vertices[j] - m.vertices[i], m.vertices[k] - m.vertices[i])) / 2
            c += a * (m.vertices[i] + m.vertices[j] + m.vertices[k]) / 3
            A += a
        end
    end
    return V, (A > 0 ? c / A : zero(Vec3))
end

_diff_materials(t::Target, part::Symbol) = [(c.material.name, c.kind, c.two_sided) for c in t.components if c.part === part]

function _diff_parts(ta::Target, tb::Target, tol::Real)
    rows = _DiffPart[]
    pa, pb = parts(ta), parts(tb)
    nanv = Vec3(NaN, NaN, NaN)
    for p in pa
        if !(p in pb)
            push!(rows, (; part = p, status = :removed, shift = nanv, distance = NaN, displacement = NaN, material = false))
            continue
        end
        (Va, ca), (Vb, cb) = _diff_part_geometry(ta, p), _diff_part_geometry(tb, p)
        material = _diff_materials(ta, p) != _diff_materials(tb, p)
        if length(Va) == length(Vb)
            D = Vb .- Va
            disp = maximum(norm, D; init = 0.0)
            rigid = all(d -> norm(d - D[1]) <= tol, D)
            status = disp <= tol ? :same : rigid ? :moved : :changed
            shift = status === :moved ? D[1] : status === :same ? zero(Vec3) : cb - ca
            push!(rows, (; part = p, status, shift, distance = norm(shift), displacement = disp, material))
        else
            push!(rows, (; part = p, status = :changed, shift = cb - ca, distance = norm(cb - ca), displacement = NaN, material))
        end
    end
    for p in pb
        p in pa || push!(rows, (; part = p, status = :added, shift = nanv, distance = NaN, displacement = NaN, material = false))
    end
    return rows
end

function _diff_beams(ba, bb, tol::Real)
    rows = _DiffBeam[]
    ib = Dict(b.name => b for b in bb)
    na = Set(b.name for b in ba)
    nanv = Vec3(NaN, NaN, NaN)
    for a in ba
        b = get(ib, a.name, nothing)
        if b === nothing
            push!(rows, (; beam = a.name, status = :removed, aim_shift = nanv, distance = NaN,
                         power = (Float64(a.power), NaN), defocus = (Float64(a.defocus), NaN), changes = Symbol[]))
            continue
        end
        shift = Vec3(b.aim - a.aim)
        changes = Symbol[]
        norm(shift) > tol && push!(changes, :aim)
        abs(b.power - a.power) > tol * max(1.0, abs(a.power)) && push!(changes, :power)
        abs(b.defocus - a.defocus) > tol && push!(changes, :defocus)
        a.spot == b.spot || push!(changes, :spot)
        norm(b.port - a.port) > tol && push!(changes, :port)
        push!(rows, (; beam = a.name, status = isempty(changes) ? :same : :changed, aim_shift = shift,
                     distance = norm(shift), power = (Float64(a.power), Float64(b.power)),
                     defocus = (Float64(a.defocus), Float64(b.defocus)), changes))
    end
    for b in bb
        b.name in na || push!(rows, (; beam = b.name, status = :added, aim_shift = nanv, distance = NaN,
                                     power = (NaN, Float64(b.power)), defocus = (NaN, Float64(b.defocus)), changes = Symbol[]))
    end
    return rows
end

# The rows of a report by (rule, beam, component). The tightest row is kept when a key repeats.
function _diff_talis_index(r::TaLISReport)
    d = Dict{Tuple{Symbol,String,Symbol},TaLISRow}()
    order = Tuple{Symbol,String,Symbol}[]
    for x in r.rows
        k = (x.rule, x.beam, x.component)
        if !haskey(d, k)
            d[k] = x
            push!(order, k)
        elseif x.margin < d[k].margin
            d[k] = x
        end
    end
    return d, order
end

function _diff_talis(ra::TaLISReport, rb::TaLISReport, tol::Real)
    (da, oa), (db, ob) = _diff_talis_index(ra), _diff_talis_index(rb)
    rows = _DiffTaLIS[]
    for k in oa
        x = da[k]
        y = get(db, k, nothing)
        if y === nothing
            push!(rows, (; rule = k[1], beam = k[2], component = k[3], unit = x.unit, status = (x.status, :absent),
                         margin = (x.margin, NaN), delta = NaN, change = :removed))
            continue
        end
        delta = x.margin == y.margin ? 0.0 : y.margin - x.margin
        fa, fb = x.status === :fail, y.status === :fail
        change = fb && !fa ? :newly_failing : fa && !fb ? :fixed : delta < -tol ? :worse : delta > tol ? :improved : :same
        push!(rows, (; rule = k[1], beam = k[2], component = k[3], unit = x.unit, status = (x.status, y.status),
                     margin = (x.margin, y.margin), delta, change))
    end
    for k in ob
        haskey(da, k) && continue
        y = db[k]
        push!(rows, (; rule = k[1], beam = k[2], component = k[3], unit = y.unit, status = (:absent, y.status),
                     margin = (NaN, y.margin), delta = NaN, change = y.status === :fail ? :newly_failing : :added))
    end
    return rows
end

# Role-aware quantities of one evaluation, as (role, quantity, unit) => value in order.
function _diff_quantities(ev::Evaluation, foils)
    s = ev.scenario
    t = s.target
    q = Pair{NTuple{3,String},Float64}[]
    v = ForwardDiff.value
    dep = ev.deposition.element
    push!(q, ("laser", "absorbed power", "TW") => v(sum(dep)))
    for p in parts(t)
        es = elements(t, p)
        P = v(sum(dep[e] for e in es; init = zero(eltype(dep))))
        P > 0 && push!(q, ("part $p", "absorbed power", "TW") => P)
    end
    for f in foils
        es = [e for e in elements(t, Symbol(f)) if t.elem_side[e] == 1]
        isempty(es) && continue
        A = t.elem_area[es]
        I = v.(dep[es]) ./ A
        m = sum(I .* A) / sum(A)
        push!(q, ("foil $f", "absorbed intensity", "TW/mm2") => m)
        push!(q, ("foil $f", "intensity RMS over mean", "%") => (m > 0 ? 100 * sqrt(sum(A .* (I .- m) .^ 2) / sum(A)) / m : NaN))
    end
    obs = s.observed
    if !isempty(obs) && s.model !== nothing
        A = sum(t.elem_area[e] for e in obs)
        names = join(unique(t.components[t.elem_comp[e]].part for e in obs), ", ")
        if norm(sum(t.elem_area[e] * t.elem_normal[e] for e in obs)) < 0.05 * A
            a = v.(flux_moments(ev; lmax = 4))
            role = "capsule ($names)"
            push!(q, (role, "mean flux", "TW/mm2") => a[1])
            push!(q, (role, "Tr", "eV") => radiation_temperature(a[1]))
            push!(q, (role, "P2/P0", "%") => 100a[3] / a[1])
            push!(q, (role, "P4/P0", "%") => 100a[5] / a[1])
            vel = v.(hotspot_velocity(ev))
            push!(q, (role, "mode 1", "%") => 100 * norm(v.(mode1_vector(flux_ylm(ev; lmax = 1)))))
            for (i, ax) in enumerate(("x", "y", "z"))
                push!(q, (role, "hot-spot velocity $ax", "km/s") => vel[i])
            end
        else
            role = "face ($names)"
            H = v(face_mean_flux(ev))
            push!(q, (role, "mean flux", "TW/mm2") => H)
            push!(q, (role, "Tr", "eV") => radiation_temperature(H))
            push!(q, (role, "flux RMS over mean", "%") => 100v(face_nonuniformity(ev)))
        end
    end
    return q
end

function _diff_physics(eva::Evaluation, evb::Evaluation, foils)
    qa, qb = _diff_quantities(eva, foils), _diff_quantities(evb, foils)
    db = Dict(qb)
    da = Dict(qa)
    rows = _DiffPhysics[]
    for (k, x) in qa
        y = get(db, k, NaN)
        push!(rows, (; role = k[1], quantity = k[2], unit = k[3], a = x, b = y, delta = y - x))
    end
    for (k, y) in qb
        haskey(da, k) || push!(rows, (; role = k[1], quantity = k[2], unit = k[3], a = NaN, b = y, delta = NaN))
    end
    return rows
end

"""
    design_diff(target_a, beams_a, target_b, beams_b; names = ("reference", "current"),
                reports = nothing, talis = nothing, evaluations = nothing, foils = (),
                tolerance = 1e-9, margin_tolerance = 1e-6) -> DesignDiff

Compares design `b` with the reference design `a`. Parts are matched by part name and
beams by beam name, so a renamed part is one part removed and one added.

- Parts are compared at every vertex. A part whose vertices all move by the same vector
  (within `tolerance`, mm) is `:moved` and any other change is `:changed`, with the
  shift of the area centroid.
- Beams are compared in aim, power, defocus, spot and port.
- TaLIS rows are compared when `reports = (report_a, report_b)` is given, or when
  `talis` gives the keywords of [`talis_check`](@ref), which is then run on both
  designs. Rows are matched by rule, beam and component. A row that fails only in `b` is
  `:newly_failing`, one that fails only in `a` is `:fixed`, and any other change of the
  margin beyond `margin_tolerance` is `:worse` or `:improved`.
- The physics is compared when `evaluations = (ev_a, ev_b)` is given, over the roles
  that each design has: the absorbed laser power in total and per part, the absorbed
  intensity and its nonuniformity on each part of `foils`, and on the observed elements
  of each scenario the mean flux, radiation temperature, `P2/P0`, `P4/P0`, mode 1 and
  hot-spot velocity when the elements close on themselves, as a capsule does, or the
  mean flux, temperature and nonuniformity of an open face. A quantity that only one
  design has is `NaN` on the other side.

[`design_diff_html`](@ref) and [`write_design_diff`](@ref) give the page of the result.
"""
function design_diff(ta::Target, ba::AbstractVector{<:Beam}, tb::Target, bb::AbstractVector{<:Beam};
                     names = ("reference", "current"), reports = nothing, talis = nothing,
                     evaluations = nothing, foils = (), tolerance::Real = 1e-9, margin_tolerance::Real = 1e-6)
    if reports === nothing && talis !== nothing
        reports = (talis_check(ta, ba; talis...), talis_check(tb, bb; talis...))
    end
    fs = foils isa Symbol || foils isa AbstractString ? (Symbol(foils),) : Tuple(Symbol(f) for f in foils)
    return DesignDiff((String(names[1]), String(names[2])), _diff_parts(ta, tb, tolerance), _diff_beams(ba, bb, tolerance),
                      reports === nothing ? _DiffTaLIS[] : _diff_talis(reports[1], reports[2], margin_tolerance),
                      evaluations === nothing ? _DiffPhysics[] : _diff_physics(evaluations[1], evaluations[2], fs))
end

const _DIFF_CHANGE_RANK = Dict(:newly_failing => 0, :worse => 1, :removed => 2, :added => 3, :improved => 4, :fixed => 5, :same => 6)
const _DIFF_CHANGE_CLASS = Dict(:newly_failing => "fail", :worse => "review", :fixed => "pass", :improved => "pass")

"""
    design_diff_html(diff::DesignDiff; title = "Design diff", generated = nothing, unchanged = false) -> String

The difference of two designs as one self-contained HTML page: the parts and beams that
differ, the TaLIS rows that changed with the newly failing ones first, and the physics
of each role. With `unchanged = true` the rows without a change are listed too. The page
refers to no other file and the same input gives the same text; `generated` is a line
for the foot of the page, such as a date.
"""
function design_diff_html(d::DesignDiff; title::AbstractString = "Design diff", generated = nothing,
                          unchanged::Bool = false)
    io = IOBuffer()
    keep(v, f) = unchanged ? v : filter(f, v)
    n(v, f) = count(f, v)
    println(io, "<h2>Summary</h2>")
    _html_table(io, ("", "differ", "of"),
                [("parts", n(d.parts, x -> x.status !== :same), length(d.parts)),
                 ("beams", n(d.beams, x -> x.status !== :same), length(d.beams)),
                 ("TaLIS rows", n(d.talis, x -> x.change !== :same), length(d.talis)),
                 ("TaLIS rows newly failing", n(d.talis, x -> x.change === :newly_failing), length(d.talis)),
                 ("TaLIS rows fixed", n(d.talis, x -> x.change === :fixed), length(d.talis))]; sortable = false)
    ps = keep(d.parts, x -> x.status !== :same || x.material)
    println(io, "<h2>Parts</h2>")
    if isempty(ps)
        println(io, "<p>No part differs.</p>")
    else
        _html_table(io, ("part", "change", "shift (mm)", "distance (mm)", "largest vertex displacement (mm)", "material"),
                    [(string(x.part), string(x.status), _packet_vec(x.shift), _html_num(x.distance, _packet_mm(x.distance)),
                      _html_num(x.displacement, _packet_mm(x.displacement)), x.material ? "changed" : "") for x in ps]; id = "parts")
    end
    bs = keep(d.beams, x -> x.status !== :same)
    println(io, "<h2>Beams</h2>")
    if isempty(bs)
        println(io, "<p>No beam differs.</p>")
    else
        _html_table(io, ("beam", "change", "aim shift (mm)", "distance (mm)", "power (TW)", "defocus (mm)"),
                    [(x.beam, x.status === :changed ? join(x.changes, ", ") : string(x.status), _packet_vec(x.aim_shift),
                      _html_num(x.distance, _packet_mm(x.distance)),
                      @sprintf("%.4g to %.4g", x.power...), @sprintf("%.3f to %.3f", x.defocus...)) for x in bs]; id = "beams")
    end
    if !isempty(d.talis)
        ts = sort(keep(d.talis, x -> x.change !== :same); by = x -> (_DIFF_CHANGE_RANK[x.change], x.delta))
        println(io, "<h2>TaLIS margins</h2>")
        if isempty(ts)
            println(io, "<p>No row changed.</p>")
        else
            _html_table(io, ("change", "rule", "beam", "component", d.names[1], d.names[2], "margin " * d.names[1],
                             "margin " * d.names[2], "delta", "unit"),
                        [(_html_chip(get(_DIFF_CHANGE_CLASS, x.change, "none"), replace(string(x.change), "_" => " ");
                                     rank = _DIFF_CHANGE_RANK[x.change]),
                          string(x.rule), x.beam, string(x.component), string(x.status[1]), string(x.status[2]),
                          _html_num(x.margin[1], _packet_mm(x.margin[1])), _html_num(x.margin[2], _packet_mm(x.margin[2])),
                          _html_num(x.delta, _packet_mm(x.delta)), string(x.unit)) for x in ts]; id = "talis")
        end
    end
    if !isempty(d.physics)
        println(io, "<h2>Physics</h2>")
        g(x) = isfinite(x) ? @sprintf("%.5g", x) : string(x)
        _html_table(io, ("role", "quantity", d.names[1], d.names[2], "delta", "unit"),
                    [(x.role, x.quantity, _html_num(x.a, g(x.a)), _html_num(x.b, g(x.b)), _html_num(x.delta, g(x.delta)), x.unit)
                     for x in d.physics]; id = "physics")
    end
    return _html_document(title, String(take!(io)); lead = string(d.names[1], " against ", d.names[2]), generated)
end

function _diff_data(d::DesignDiff)
    return Dict("names" => collect(d.names),
                "parts" => [Dict("part" => x.part, "status" => x.status, "shift_mm" => _json_vec(x.shift),
                                 "distance_mm" => _json_number(x.distance), "displacement_mm" => _json_number(x.displacement),
                                 "material_changed" => x.material) for x in d.parts],
                "beams" => [Dict("beam" => x.beam, "status" => x.status, "aim_shift_mm" => _json_vec(x.aim_shift),
                                 "distance_mm" => _json_number(x.distance), "power_TW" => _json_vec(x.power),
                                 "defocus_mm" => _json_vec(x.defocus), "changes" => x.changes) for x in d.beams],
                "talis" => [Dict("rule" => x.rule, "beam" => x.beam, "component" => x.component, "unit" => x.unit,
                                 "status" => collect(x.status), "margin" => _json_vec(x.margin),
                                 "delta" => _json_number(x.delta), "change" => x.change) for x in d.talis],
                "physics" => [Dict("role" => x.role, "quantity" => x.quantity, "unit" => x.unit, "a" => _json_number(x.a),
                                   "b" => _json_number(x.b), "delta" => _json_number(x.delta)) for x in d.physics])
end

"""
    write_design_diff(dir, diff::DesignDiff; title = "Design diff", generated = nothing, json = false,
                      unchanged = false) -> Vector{String}

Writes the difference of two designs into the folder `dir` and returns the paths:
`diff.html` ([`design_diff_html`](@ref)) and, with `json = true`, `diff.json` with every
row of the four tables, the unchanged ones included.
"""
function write_design_diff(dir::AbstractString, d::DesignDiff; title::AbstractString = "Design diff", generated = nothing,
                           json::Bool = false, unchanged::Bool = false)
    mkpath(dir)
    paths = String[joinpath(dir, "diff.html")]
    write(paths[1], design_diff_html(d; title, generated, unchanged))
    json && push!(paths, _write_json(joinpath(dir, "diff.json"), _diff_data(d)))
    return paths
end
