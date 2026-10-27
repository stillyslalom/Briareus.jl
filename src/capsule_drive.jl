# The drive handed to a capsule-only hydrodynamics run (ROADMAP "Briareus → CompactLES", LLNL's two-step
# method; reference/litreview/hydro_coupling.md §(c)).

"""
    CapsuleDrive

The radiation drive on the capsule over time, for a capsule-only hydrodynamics run
(LLNL's two-step method: an integrated hohlraum step, then capsule-only runs driven by
its capsule flux; Clark et al., PoP 26, 050601 (2019); Masse et al., PoP 26 (2019)).
Build one with [`capsule_drive`](@ref). Fields:

- `time` — uniform time grid (ns)
- `Tr` — radiation temperature of the capsule-averaged incident flux (eV),
  `radiation_temperature(flux)`
- `flux` — capsule-averaged incident flux `a₀` (TW/mm²; ×1e14 → W/cm²)
- `modes` — `(lmax+1) × nt` Legendre ratios `aₗ/a₀` (row `l+1`, so row 1 is 1), of
  the incident flux `F(θ,t) = a₀(t) Σₗ (aₗ/a₀)(t) Pₗ(cos θ)`, `θ` from the hohlraum axis
  (+z) about the capsule center. **`P2 > 0` is pole-hot** (Masse 2019's convention).
- `absorbed` — x-ray power absorbed by the capsule (TW), as in `energy_balance`
  (`NaN` where the history does not record it)
- `radius` — the radius (mm) of the sphere the flux was evaluated on (the capsule mesh)
- `area` — the observed elements' area (mm²)
- `provenance` — `Dict{Symbol,Any}`: code, version, git commit, the history's source
  times, observed components, capsule material, mesh size, and whatever the caller adds
- `ylm` — `nothing` (axisymmetric), or the drive with every `m`: the `(lmax+1)² × nt`
  ratios `bₗₘ` (rows in [`mode_index`](@ref) order, so row 1 is 1) of
  `F(θ,φ,t) = a₀(t) Σₗₘ bₗₘ(t) Ŷₗₘ(θ,φ)`, where `Ŷₗₘ` are the real Schmidt
  semi-normalized harmonics, `Ŷₗ₀ = Pₗ(cos θ)` and
  `Ŷₗ,±ₘ = √(2(l−m)!/(l+m)!) Pₗᵐ(cos θ) (cos mφ, sin mφ)` for `m > 0`, without the
  Condon–Shortley phase, `φ` from +x toward +y (NIF: Switchyard 1). That is the
  normalization of `modes`: the `m = 0` rows *are* `modes` (`bₗ₀ = aₗ/a₀`), and
  `√(Σₘ bₗₘ²)` is `mode_amplitude`, the rotation-invariant `|aₗ/a₀|` of degree `l`.
  From the orthonormal coefficients `cₗₘ` of `ylm_moments` (or
  `drive_history(...; ylm = true)`), `bₗₘ = √((2l+1)/4π) cₗₘ / a₀`.

The angular structure is that of the flux on the initial capsule radius; applying it
at a moving ablation front far inside that radius changes the geometric smoothing
(`hydro_coupling.md` §2, pitfall 2). Gray: no spectral (M-band) split yet (ROADMAP "Spectral radiosity").
"""
struct CapsuleDrive
    time::Vector{Float64}
    Tr::Vector{Float64}
    flux::Vector{Float64}
    modes::Matrix{Float64}
    absorbed::Vector{Float64}
    radius::Float64
    area::Float64
    provenance::Dict{Symbol,Any}
    ylm::Union{Nothing,Matrix{Float64}}
    function CapsuleDrive(time, Tr, flux, modes, absorbed, radius, area, provenance,
                          ylm = nothing)
        nt = length(time)
        (length(Tr) == nt && length(flux) == nt && length(absorbed) == nt && size(modes, 2) == nt) ||
            throw(ArgumentError("time, Tr, flux, absorbed and the columns of modes differ in length"))
        size(modes, 1) >= 1 || throw(ArgumentError("modes needs at least the l = 0 row"))
        issorted(time) || throw(ArgumentError("times must be sorted"))
        radius > 0 || throw(ArgumentError("radius must be positive"))
        if ylm !== nothing
            L = size(modes, 1) - 1
            size(ylm) == ((L + 1)^2, nt) ||
                throw(ArgumentError("ylm must be (lmax+1)² × nt = $((L + 1)^2) × $nt, as modes"))
            all(l -> view(ylm, mode_index(l, 0), :) == view(modes, l + 1, :), 0:L) ||
                throw(ArgumentError("the m = 0 rows of ylm must equal modes (bₗ₀ = aₗ/a₀)"))
        end
        new(collect(Float64, time), collect(Float64, Tr), collect(Float64, flux),
            Matrix{Float64}(modes), collect(Float64, absorbed), Float64(radius),
            Float64(area), Dict{Symbol,Any}(provenance),
            ylm === nothing ? nothing : Matrix{Float64}(ylm))
    end
end

"Highest Legendre order a `CapsuleDrive` carries."
_lmax(d::CapsuleDrive) = size(d.modes, 1) - 1

function Base.show(io::IO, d::CapsuleDrive)
    print(io, "CapsuleDrive(", length(d.time), " times on [", first(d.time), ", ", last(d.time),
          "] ns, lmax = ", _lmax(d), d.ylm === nothing ? "" : " (all m)", ", R = ", d.radius, " mm, peak Tr = ",
          round(maximum(d.Tr); digits = 1), " eV)")
end

"""
    (d::CapsuleDrive)(t) -> (; Tr, flux, modes, ylm)

The drive at time `t` (ns), linear in the moments `aₗ` (and `a₀ bₗₘ`) between grid
times and held at the end values outside them; `ylm` is `nothing` for an axisymmetric
drive.
"""
function (d::CapsuleDrive)(t::Real)
    a = [_interp(d.time, view(d.modes, l + 1, :) .* d.flux, t) for l in 0:_lmax(d)]
    flux = a[1]
    ratio(v) = flux > 0 ? v ./ flux : [1; zeros(length(v) - 1)]
    ylm = d.ylm === nothing ? nothing :
          ratio([_interp(d.time, view(d.ylm, q, :) .* d.flux, t) for q in axes(d.ylm, 1)])
    return (; Tr = radiation_temperature(flux), flux, modes = ratio(a), ylm)
end

"""
    capsule_drive(history, target; observe = :capsule, dt = 0.05, times = nothing,
                  lmax = size(history.moments, 1) - 1, radius = nothing,
                  ylm = history.ylm !== nothing, provenance = (;)) -> CapsuleDrive

The capsule drive of a [`drive_history`](@ref) result for a capsule-only hydrodynamics
run, on a uniform grid of spacing `≤ dt` (ns; 50 ps by default, `hydro_coupling.md`
§(c)) spanning the history's solved times, or on the given `times` (within them). The
Legendre moments `aₗ` are interpolated linearly between the history's steps and
divided by the interpolated `a₀`; `flux` is `a₀` and `Tr` its radiation temperature,
exactly the history's `flux` and `Tr` at its own times. Run `drive_history(...; lmax = 8)`
for the moments through `l = 8` the hand-over calls for; `lmax` here may only truncate.

With `ylm` (the default when the history recorded them,
`drive_history(...; ylm = true)`), the drive carries every `m` as well: the history's
orthonormal coefficients `cₗₘ`, interpolated linearly like `aₗ`, become the Schmidt ratios
`bₗₘ = √((2l+1)/4π) cₗₘ / a₀` of the `ylm` field, whose `m = 0` rows are then set to
`modes`, so the two agree exactly. (The history integrates its Legendre moments one
sub-triangle level finer than its `Y_lm`; the replaced rows differ by less than 1e-3
of `a₀` on the coarse test hohlraum.)

`target` (a `Target` or `RadiationModel`; the keyframes' reference target for a moving
hohlraum) supplies the radius of the observed sphere (its `SphereSurface`; pass
`radius` for any other surface), the observed area, and the provenance of the capsule
material and mesh. `observe` must name the components the history observed.
`provenance` (NamedTuple or Dict: hohlraum, pulse, albedo model, …) is merged into the
recorded provenance.
"""
function capsule_drive(h::NamedTuple, target; observe = :capsule, dt::Real = 0.05,
                       times = nothing, lmax::Integer = size(h.moments, 1) - 1,
                       radius = nothing, ylm::Bool = get(h, :ylm, nothing) !== nothing,
                       provenance = (;))
    t = target isa RadiationModel ? target.target : target isa Keyframes ? target.reference : target
    t isa Target || throw(ArgumentError("target must be a Target, RadiationModel or Keyframes"))
    names = observe isa Symbol ? (observe,) : Tuple(observe)
    elems = elements(t, names...)
    isempty(elems) && throw(ArgumentError("no elements in $(names)"))
    0 <= lmax <= size(h.moments, 1) - 1 ||
        throw(ArgumentError("the history carries moments through l = $(size(h.moments, 1) - 1); " *
                            "run drive_history(...; lmax = $lmax)"))
    th = h.time
    length(th) >= 2 || throw(ArgumentError("the history needs at least two times"))
    if times === nothing
        dt > 0 || throw(ArgumentError("dt must be positive"))
        n = max(ceil(Int, (th[end] - th[1]) / dt - 1e-9), 1)
        grid = collect(range(th[1], th[end]; length = n + 1))
    else
        grid = collect(Float64, times)
        issorted(grid) || throw(ArgumentError("times must be sorted"))
        (grid[1] >= th[1] - 1e-12 && grid[end] <= th[end] + 1e-12) ||
            throw(ArgumentError("times must lie within the history's [$(th[1]), $(th[end])] ns"))
    end
    nt = length(grid)
    A = zeros(lmax + 1, nt)
    for l in 0:lmax, k in 1:nt
        A[l+1, k] = _interp(th, view(h.moments, l + 1, :), grid[k])
    end
    flux = A[1, :]
    modes = A ./ reshape(map(x -> x > 0 ? x : Inf, flux), 1, :)
    modes[1, :] .= 1
    B = nothing
    if ylm
        C = get(h, :ylm, nothing)
        C === nothing && throw(ArgumentError("the history has no Y_lm; run drive_history(...; ylm = true)"))
        size(C, 1) >= (lmax + 1)^2 || throw(ArgumentError("the history's Y_lm stop below l = $lmax"))
        B = zeros((lmax + 1)^2, nt)
        for l in 0:lmax, m in -l:l
            q, s = mode_index(l, m), sqrt((2l + 1) / 4π)
            for k in 1:nt
                B[q, k] = m == 0 ? modes[l+1, k] :
                          flux[k] > 0 ? s * _interp(th, view(C, q, :), grid[k]) / flux[k] : 0.0
            end
        end
    end
    absorbed = length(names) == 1 && haskey(h.absorbed, names[1]) ?
               [_interp(th, h.absorbed[names[1]], τ) for τ in grid] : fill(NaN, nt)
    R = radius === nothing ? _observed_radius(t, elems) : Float64(radius)
    area = sum(t.elem_area[elems])
    prov = Dict{Symbol,Any}(
        :code => "Briareus", :version => string(pkgversion(@__MODULE__)), :git => _git_commit(),
        :quantity => "incident x-ray flux on the capsule (gray), Legendre about +z, P2 > 0 pole-hot" *
                     (B === nothing ? "" : "; ylm: Schmidt semi-normalized real Y_lm ratios, phi from +x"),
        :units => "time ns, Tr eV, flux TW/mm^2, absorbed TW, radius mm, area mm^2",
        :observe => collect(Symbol, names), :history_times => (th[1], th[end], length(th)),
        :interpolation => "linear in a_l between history steps",
        :capsule_material => string(element_material(t, elems[1]).albedo),
        :elements => length(elems), :target_elements => nelements(t),
        :components => [c.name for c in t.components])
    for (k, v) in pairs(provenance)
        prov[Symbol(k)] = v
    end
    return CapsuleDrive(grid, radiation_temperature.(flux), flux, modes, absorbed, R, area, prov, B)
end

# The radius of the sphere the observed elements lie on.
function _observed_radius(t::Target, elems)
    s = t.components[t.elem_comp[elems[1]]].surface
    for e in elems
        se = t.components[t.elem_comp[e]].surface
        (se isa SphereSurface && se.radius == s.radius && se.center == s.center) ||
            throw(ArgumentError("the observed elements are not on one sphere; pass `radius`"))
    end
    return s.radius
end

# The package checkout's commit (with "-dirty" for uncommitted changes), or "unknown".
function _git_commit()
    dir = pkgdir(@__MODULE__)
    try
        sha = readchomp(pipeline(`git -C $dir rev-parse --short HEAD`; stderr = devnull))
        dirty = !isempty(readchomp(pipeline(`git -C $dir status --porcelain --untracked-files=no`;
                                            stderr = devnull)))
        return dirty ? sha * "-dirty" : sha
    catch
        return "unknown"
    end
end
