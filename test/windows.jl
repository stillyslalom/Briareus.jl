# Diagnostic windows in the hohlraum wall: conforming mesh, symmetry, energy balance,
# an integrating-cavity estimate of the drive loss, and laser transmission. First, the
# exact curved area of cylinder patches that window (and wall) areas rest on.

@testset "cylinder patch areas" begin
    # Reference: the radially projected area as an integral over the flat triangle,
    # ∬ (R/ρ)|ρ̂·n| dA (ρ the distance from the axis, n the triangle normal), Duffy-mapped
    # to the unit square and integrated by 40×40-point Gauss–Legendre in 256-bit
    # arithmetic from the exact Float64 vertices. It shares nothing with the closed form
    # (Green's theorem along the projected edges) but the geometry.
    setprecision(BigFloat, 256) do
        n = 40
        xs, ws = BigFloat[], BigFloat[]
        for i in 1:n                                  # Legendre roots by Newton
            x = BigFloat(cos(π * (i - 0.25) / (n + 0.5)))
            local dp
            for _ in 1:60
                p0, p1 = one(x), x
                for k in 2:n
                    p0, p1 = p1, ((2k - 1) * x * p1 - (k - 1) * p0) / k
                end
                dp = n * (x * p1 - p0) / (x^2 - 1)
                δ = p1 / dp
                x -= δ
                abs(δ) < big(1e-70) && break
            end
            push!(xs, (x + 1) / 2)
            push!(ws, 1 / ((1 - x^2) * dp^2))
        end
        function reference(s, a, b, c)
            a, b, c = big.(a), big.(b), big.(c)
            ax, o = big.(s.axis), big.(s.origin)
            e1, e2 = b - a, c - b
            N = cross(e1, e2)
            J = norm(N)
            acc = zero(BigFloat)
            for i in 1:n, j in 1:n
                d = a + xs[i] * e1 + xs[i] * xs[j] * e2 - o
                u = d - dot(d, ax) * ax
                ρ = norm(u)
                acc += ws[i] * ws[j] * xs[i] * (s.radius / ρ) * abs(dot(u, N) / (ρ * J))
            end
            return acc * J
        end
        s = CylinderSurface(2.875)
        tilt = CylinderSurface(1.7; origin = (0.3, -0.2, 0.5), axis = (0.2, 0.1, 1.0))
        e1 = normalize(cross(tilt.axis, Vec3(1, 0, 0)))
        e2 = cross(tilt.axis, e1)
        on(φ, z; r = 2.875) = Vec3(r * cos(φ), r * sin(φ), z)
        ontilt(φ, z; r = 1.7) = tilt.origin + r * (cos(φ) * e1 + sin(φ) * e2) + z * tilt.axis
        cases = [
            (s, on(0, 0), on(2π / 48, 0), on(2π / 48, 9.4 / 36)),               # default wall half-cell
            (s, on(0, -1), on(2π / 3, -1), on(2π / 3, 2)),                       # 120° patch
            (s, on(0, -2), on(π / 2, -2), on(π / 2, 1)),                         # 90°
            (s, on(0.1, -0.7), on(1.3, 0.4), on(0.5, 1.9)),                      # general, on the surface
            (s, on(0.3, 0), on(0.3 + 1e-6, 0), on(0.3 + 1e-6, 1)),               # azimuthal sliver
            (s, on(1.1, 0), on(1.1 + 1e-9, 0), on(1.1 + 1e-9, 2)),               # thinner, longer
            (s, on(0.3, 0), on(1.0, 0), on(1.0, 1e-6)),                          # axial sliver
            (s, on(0.7, 4), on(0.7 + 1e-7, 4), on(0.7, 4 + 1e-7)),               # tiny, far from z = 0
            (s, on(0.1, 0; r = 2.86), on(0.15, 0.02; r = 2.87), on(0.12, 0.05; r = 2.8749)),  # sub-triangle, off the surface
            (s, on(0, 0; r = 2.0), on(0.01, 0; r = 2.5), on(0, 0.01; r = 2.2)),  # short edges seen obliquely
            (tilt, ontilt(0.2, 0.1), ontilt(0.4, 0.1), ontilt(0.4, 0.6)),        # tilted, offset cylinder
            (tilt, ontilt(0.2, 0.1; r = 1.6), ontilt(0.25, 0.12; r = 1.69), ontilt(0.22, 0.2; r = 1.65)),
        ]
        for (sf, a, b, c) in cases
            @test Briareus._curved_area(sf, a, b, c) ≈ reference(sf, a, b, c) rtol = 1e-14
        end
        # a sliver on a tilted cylinder is limited by its inputs: its 1.5 mm axial edge,
        # rounded at ε·1.5 in global coordinates, spans 1.7e-6 mm tangentially, so the
        # area itself is defined only to ~ε·1.5/1.7e-6 ≈ 2e-10
        a, b, c = ontilt(1.0, 0), ontilt(1.0 + 1e-6, 0), ontilt(1.0, 1.5)
        @test Briareus._curved_area(tilt, a, b, c) ≈ reference(tilt, a, b, c) rtol = 1e-10
    end
    # areas are additive under subdivision (vertices off the surface), so element
    # quadrature weights sum to the element area; a whole cylinder is 2πRL
    t = hohlraum()
    @test all(e -> sum(element_quadrature(t, e, 2)[2]) ≈ t.elem_area[e], elements(t, :wall)[1:97:end])
    @test sum(t.elem_area[elements(t, :wall)]) ≈ 2π * 2.875 * 9.4 rtol = 1e-14
end

@testset "diagnostic windows" begin
    R, L, rl, rc = 2.875, 9.4, 1.55, 0.9
    h = L / 2

    # conforming grids contain every break exactly and are uniform between breaks
    g = conforming_grid(0.0, 1.0, 10, [0.3])
    @test g ≈ 0:0.1:1 atol = 1e-15
    @test 0.3 in g
    g = conforming_grid(-1.0, 1.0, 8, [0.2, 0.2 + 1e-12, 5.0])     # merged, out of range dropped
    @test g[1] == -1 && g[end] == 1 && count(x -> abs(x - 0.2) < 1e-9, g) == 1
    gp = conforming_grid(0.0, 2π, 24, [-0.3, 0.3]; periodic = true)   # wraps, starts at a break
    @test gp[1] == 0.3 && gp[end] ≈ 0.3 + 2π && any(x -> x ≈ 2π - 0.3, gp)
    @test all(>(0), diff(gp)) && length(gp) - 1 == 24
    @test_throws ArgumentError WallWindow(phi = 0, width = 1, height = 1, open = true, material = gold())
    @test_throws ArgumentError hohlraum(windows = [WallWindow(phi = 0, width = 2, height = 1),
                                                   WallWindow(phi = 20, width = 2, height = 1)])
    @test_throws ArgumentError hohlraum(windows = [WallWindow(phi = 0, z = 4.5, width = 1, height = 1)])

    # --- geometry: the window is tiled exactly by whole patches on the curved surface
    w, ht = 2.0, 1.5
    tw = hohlraum(windows = [WallWindow(phi = 30, width = w, height = ht),
                             WallWindow(phi = 150, z = -1.0, width = 1.0, height = 0.8, open = true, name = :hole)])
    @test [c.name for c in tw.components][end-1:end] == [:window_1, :hole]
    @test tw.components[end].kind === :aperture && tw.components[end-1].kind === :surface
    @test all(c -> c.surface isa CylinderSurface, tw.components[[1, 7, 8]])
    # (φ, z) extent of each patch from its vertices; the curved patch is the cylindrical
    # rectangle R Δφ × Δz exactly (radial projection preserves z and maps chords to arcs)
    function patch_rects(m)
        map(m.patches) do r
            vs = unique(vcat([[m.vertices[i] for i in m.triangles[k]] for k in r]...))
            φ = [atan(v[2], v[1]) for v in vs]
            φ = [φ[1] + mod(x - φ[1] + π, 2π) - π for x in φ]
            (extrema(φ), extrema(v[3] for v in vs))
        end
    end
    rect = patch_rects(tw.components[end-1].mesh)
    φa, φb = deg2rad(30) - w / 2R, deg2rad(30) + w / 2R
    @test minimum(first(first(r)) for r in rect) ≈ φa atol = 1e-14
    @test maximum(last(first(r)) for r in rect) ≈ φb atol = 1e-14
    @test minimum(first(last(r)) for r in rect) ≈ -ht / 2 atol = 1e-14
    @test maximum(last(last(r)) for r in rect) ≈ ht / 2 atol = 1e-14
    exact_area(rects) = sum(R * (φ[2] - φ[1]) * (z[2] - z[1]) for (φ, z) in rects)
    @test exact_area(rect) ≈ w * ht rtol = 1e-13
    @test exact_area(patch_rects(tw.components[end].mesh)) ≈ 1.0 * 0.8 rtol = 1e-13
    @test exact_area(rect) + exact_area(patch_rects(tw.components[1].mesh)) +
          exact_area(patch_rects(tw.components[end].mesh)) ≈ 2π * R * L rtol = 1e-13
    # and the radiation model's element areas are those exact curved areas
    e1 = elements(tw, :window_1)
    @test sum(tw.elem_area[e1]) ≈ w * ht rtol = 1e-14
    @test sum(tw.elem_area[elements(tw, :hole)]) ≈ 1.0 * 0.8 rtol = 1e-14
    @test sum(tw.elem_area[elements(tw, :wall, :window_1, :hole)]) ≈ 2π * R * L rtol = 1e-14
    # caps share the wall's azimuthal grid: the rims meet vertex for vertex
    rim(m, z) = sort([atan(v[2], v[1]) for v in unique(m.vertices) if abs(v[3] - z) < 1e-12 && abs(hypot(v[1], v[2]) - R) < 1e-12])
    @test rim(tw.components[1].mesh, h) ≈ rim(tw.components[2].mesh, h) atol = 1e-12
    # watertight: rays from inside toward the wall rims and the window's edge lines all
    # hit an opaque patch when the windows are closed; through an open window they escape
    tr = hohlraum(capsule_radius = 0, windows = [WallWindow(phi = 30, width = w, height = ht)])
    o = Vec3(0.3, -0.2, 0.1)
    rng = MersenneTwister(7)
    targets_ = [Vec3(R * cos(φ), R * sin(φ), z) for φ in 2π .* rand(rng, 400) for z in (h, -h, ht / 2, -ht / 2)]
    @test all(p -> Briareus.closest_hit(tr.bvh, o, p - o)[2] != 0, targets_)
    to = hohlraum(capsule_radius = 0, windows = [WallWindow(phi = 30, width = w, height = ht, open = true)])
    toward = Vec3(R * cosd(30), R * sind(30), 0.2)
    @test Briareus.closest_hit(to.bvh, o, toward - o)[2] == 0
    hit = Briareus.closest_hit(tr.bvh, o, toward - o)[2]
    @test tr.components[tr.tri_comp[hit]].name === :window_1

    # --- symmetry: a pair of windows related by the 180° rotation about z adds no odd-m
    # harmonics (NIF quads share that symmetry); a single window drives mode 1 toward
    # itself as a flux deficit
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.0)))
         for b in nif_quads()]
    kw = (; nz = 16, nphi = 24, capsule_ntheta = 12, capsule_nphi = 24)
    capsule_modes(t) = flux_ylm(evaluate(Scenario(t, q; observe = :capsule,
                                                  viewfactor_kw = (; closed = true)), q); lmax = 4)
    oddm(c) = maximum(abs(c[mode_index(l, m)]) for l in 1:4 for m in -l:l if isodd(m)) / c[1]
    c0 = capsule_modes(hohlraum(; kw...))
    global windowless_modes = c0                  # (the same beams and mesh as below)
    pair = [WallWindow(phi = 30, width = 2.0, height = 2.0, open = true),
            WallWindow(phi = 210, width = 2.0, height = 2.0, open = true)]
    cp2 = capsule_modes(hohlraum(; kw..., windows = pair))
    @test oddm(cp2) < 1e-10
    @test cp2[1] < 0.9 * c0[1]                                   # two 4 mm² holes cost > 10%
    @test abs(cp2[mode_index(2, -2)] / cp2[1]) > 0.01            # m = 2 at the waist (sin 2φ: φ = 30°)
    cpair_closed = capsule_modes(hohlraum(; kw..., windows = [WallWindow(phi = 30, width = 2.0, height = 2.0,
                                                                         material = ablator()),
                                                              WallWindow(phi = 210, width = 2.0, height = 2.0,
                                                                         material = ablator())]))
    @test oddm(cpair_closed) < 1e-10
    # a closed window of the wall's own material is just wall (up to the remeshing)
    cw = capsule_modes(hohlraum(; kw..., windows = [WallWindow(phi = 0, width = 2.0, height = 2.0)]))
    @test mode_amplitude(cw, 1) < 1e-3 && abs(cw[1] / c0[1] - 1) < 1e-3
    for (φw, wopen) in ((0.0, true), (0.0, false), (135.0, true))
        mat = wopen ? nothing : Material(:lowZ; albedo = 0.3, conversion = 0.3)
        c1 = capsule_modes(hohlraum(; kw..., windows = [WallWindow(phi = φw, width = 2.0, height = 2.0,
                                                                   open = wopen, material = mat)]))
        d = Vec3(c1[mode_index(1, 1)], c1[mode_index(1, -1)], c1[mode_index(1, 0)]) / c1[1]
        u = Vec3(cosd(φw), sind(φw), 0)
        @test dot(d, u) < -0.005                                 # deficit toward the window
        @test norm(d - dot(d, u) * u) < 0.02 * norm(d)           # and along it (within ~1°)
        φw == 0 && @test abs(d[2]) < 1e-10                       # exact for a mirror-symmetric setup
    end

    # --- energy balance with open windows, and an integrating-cavity estimate of the drive
    # loss. A uniform source of x rays, ηs per unit area of wall and caps (area A_s), and
    # a window of area A_w at the waist with wall albedo α when closed. If every surface
    # saw the same incident flux H (the integrating-sphere idealization), global balance
    # would give S = η s A_s = H D, with
    #   D  = (1-α)(A_s + A_w) + A_LEH + (1-α_c) A_c          (window closed),
    #   D' = D + α A_w                                        (open: it re-emits nothing),
    # the sphere multiplier 1/D. The capsule at the center sees the source walls (solid
    # angle Ω_s = 4π - Ω_LEH - Ω_w) with radiosity ηs + αH, and the window with αH
    # when closed and nothing when open, so
    #   H_c ∝ Ω_s (ηs + αH) + Ω_w αH (closed),   Ω_s (ηs + αH') (open).
    # Seen from the axis, Ω_LEH = 4π(1 - h/√(h² + r²)) for both holes and the window
    # (angular width w/R, height t) has Ω_w = (w/R) t/√(R² + t²/4). The open/closed ratio
    # is then a closed form independent of the view factors. Its error is that of the
    # uniform-flux idealization: the window's fractional effect 1 - ratio is misstated by
    # at most the largest relative deviation of the solved wall flux from H.
    α, αc, η = 0.8, 0.3, 0.75
    function cavity(wopen; ww, wt)
        t = hohlraum(; kw..., windows = [WallWindow(phi = 0, width = ww, height = wt, open = wopen)])
        src = elements(t, :wall, :cap_top, :cap_bottom)
        laser = zeros(nelements(t))
        laser[src] .= t.elem_area[src]                           # s = 1 TW/mm² absorbed
        sol = solve(RadiationModel(t; closed = true), laser)
        c = elements(t, :capsule)
        return (; t, sol, Hc = sum(sol.H[c] .* t.elem_area[c]) / sum(t.elem_area[c]))
    end
    for (ww, wt) in ((2.0, 2.0), (1.0, 3.0))
        op, cl = cavity(true; ww, wt), cavity(false; ww, wt)
        # energy balance: source = absorbed + escaped + lost, and the window's escape
        # H_w A_w is the power all emitters send it (Σⱼ Aⱼ Fⱼ,w Bⱼ, by reciprocity)
        bal = energy_balance(op.sol)
        @test bal.source ≈ sum(values(bal.absorbed)) + sum(values(bal.escaped)) + bal.lost rtol = 1e-12
        ew = elements(op.t, :window_1)
        F = op.sol.model.F
        @test bal.escaped[:window_1] ≈ sum(op.t.elem_area[j] * F[j, i] * op.sol.B[j]
                                           for i in ew, j in 1:nelements(op.t)) rtol = 1e-10
        @test bal.escaped[:window_1] ≈ sum(op.sol.H[ew] .* op.t.elem_area[ew]) rtol = 1e-12
        # integrating-cavity estimate
        Aw = ww * wt
        As = 2π * R * L - Aw + 2π * (R^2 - rl^2)
        D = (1 - α) * (As + Aw) + 2π * rl^2 + (1 - αc) * 4π * rc^2
        H, Hp = η * As / D, η * As / (D + α * Aw)
        ΩL = 4π * (1 - h / sqrt(h^2 + rl^2))
        Ωw = (ww / R) * wt / sqrt(R^2 + wt^2 / 4)
        Ωs = 4π - ΩL - Ωw
        ratio = Ωs * (η + α * Hp) / (Ωs * (η + α * H) + Ωw * α * H)
        dev = maximum(abs.(cl.sol.H[elements(cl.t, :wall)] ./ H .- 1))
        @test abs(op.Hc / cl.Hc - ratio) < (1 - ratio) * dev
        @test 0.03 < 1 - op.Hc / cl.Hc < 0.1
    end

    # --- laser: an open window transmits, a closed one absorbs, and the operator conserves
    # power. One outer quad aimed through the LEH lands near z = 2.07 mm at φ = -163.7°;
    # a 1.5 mm window there passes most of it.
    b = q[findfirst(x -> x.name == "Q14T", q)]
    φw, zw = -163.71, 2.0744
    pw = Dict{Bool,Any}()
    for wopen in (true, false)
        t = hohlraum(; kw..., windows = [WallWindow(phi = φw, z = zw, width = 1.5, height = 1.5, open = wopen)])
        pw[wopen] = (op = component_power(deposit(DepositionOperator(t, [b]), [b]), t),
                     rt = component_power(raytrace_deposition(t, [b]; nspot = 241), t))
    end
    P = b.power
    # closed: every bit of the beam lands (operator to its documented 1e-4)
    @test abs(pw[false].op[:missed]) < 2e-4 * P
    @test abs(pw[false].rt[:missed]) < 1e-12 * P
    @test pw[false].op[:window_1] > 0.9 * P
    # open: what the closed window absorbed now leaves through it, exactly (same
    # quadrature points, same visibility elsewhere)
    for m in (:op, :rt)
        @test pw[true][m][:missed] - pw[false][m][:missed] ≈ pw[false][m][:window_1] rtol = 1e-10
        @test pw[true][m][:wall] ≈ pw[false][m][:wall] rtol = 1e-12
    end
    # and the operator's transmitted power matches the converged ray trace
    @test abs(pw[true].op[:missed] - pw[true].rt[:missed]) < 1e-3 * P
    # keep-out ignores wall windows (they are not entrance holes)
    @test Set(ApertureKeepOut(tw).names) == Set([:leh_top, :leh_bottom])
end

@testset "circular windows and gaps" begin
    R, L = 2.875, 9.4
    h = L / 2
    @test_throws ArgumentError WallWindow(phi = 0, diameter = 1, width = 1, height = 1)
    @test_throws ArgumentError WallWindow(phi = 0, width = 1)
    @test_throws ArgumentError WallWindow(phi = 0, diameter = 1, gap = -0.1)
    @test_throws ArgumentError WallWindow(phi = 0, diameter = 1, segments = 6)
    @test_throws ArgumentError hohlraum(windows = [WallWindow(phi = 0, diameter = 2),
                                                   WallWindow(phi = 40, diameter = 2)])
    @test_throws ArgumentError hohlraum(windows = [WallWindow(phi = 0, z = 3.9, diameter = 1.5)])
    # the equal-area polygon has the disk's area
    @test let n = 32, r = 1.3, ρ = equal_area_radius(r, n)
        n / 2 * ρ^2 * sin(2π / n) ≈ π * r^2
    end

    # --- geometry: exact areas, every vertex on the cylinder, triangles only
    d, g = 2.0, 0.05
    pair = [WallWindow(phi = 78, diameter = d, gap = g),
            WallWindow(phi = 258, diameter = d, gap = g, name = :south)]
    t = hohlraum(windows = pair)
    names = [c.name for c in t.components]
    @test names[end-3:end] == [:window_1, :window_1_gap, :south, :south_gap]
    @test t.components[end-3].kind === :surface && t.components[end-2].kind === :aperture
    for n in (:window_1, :south)
        @test sum(t.elem_area[elements(t, n)]) ≈ π * d^2 / 4 rtol = 1e-12
    end
    for n in (:window_1_gap, :south_gap)
        @test sum(t.elem_area[elements(t, n)]) ≈ π * ((d / 2 + g)^2 - (d / 2)^2) rtol = 1e-12
    end
    @test sum(t.elem_area[elements(t, :wall, :window_1, :window_1_gap, :south, :south_gap)]) ≈
          2π * R * L rtol = 1e-13
    for c in t.components[[1; length(t.components)-3:length(t.components)]]
        @test all(v -> abs(hypot(v[1], v[2]) - R) < 1e-12, c.mesh.vertices)
    end
    @test all(r -> length(r) == 1, t.components[end].mesh.patches)
    # the window's edge: n sides on the equal-area circle, in the unrolled coordinates
    uv(v, φc) = (R * (mod(atan(v[2], v[1]) - deg2rad(φc) + π, 2π) - π), v[3])
    ρs = [hypot(uv(v, 78)...) for v in t.components[end-3].mesh.vertices]
    n = Briareus._window_segments(pair[1], min(2π * R / 48, L / 36))
    @test maximum(ρs) ≈ equal_area_radius(d / 2, n) rtol = 1e-12
    @test count(ρ -> abs(ρ - maximum(ρs)) < 1e-12, unique(ρs)) <= n
    # explicit segment count; a gap-free window is just the disk
    t1 = hohlraum(windows = [WallWindow(phi = 10, z = -1.0, diameter = 1.2, segments = 64, open = true)])
    @test sum(t1.elem_area[elements(t1, :window_1)]) ≈ π * 0.6^2 rtol = 1e-12
    @test maximum(hypot(uv(v, 10)[1], v[3] + 1) for v in t1.components[end].mesh.vertices) ≈
          equal_area_radius(0.6, 64) rtol = 1e-12

    # --- watertight: with the window and its gap closed every ray from inside toward the
    # window's region hits an opaque patch; with them open, rays toward the window and
    # gap centerlines escape
    closed = hohlraum(capsule_radius = 0, windows = [WallWindow(phi = 78, diameter = d, gap = g,
                                                                gap_material = gold())])
    opened = hohlraum(capsule_radius = 0, windows = [WallWindow(phi = 78, diameter = d, gap = g,
                                                                open = true)])
    o = Vec3(0.3, -0.2, 0.1)
    rng = MersenneTwister(11)
    onwall(u, v) = (φ = deg2rad(78) + u / R; Vec3(R * cos(φ), R * sin(φ), v))
    pts = [(ρ * cos(θ), ρ * sin(θ)) for θ in 2π .* rand(rng, 500)
           for ρ in (d / 2, d / 2 + g, d / 2 + g / 2, 1.4 * rand(rng), d / 2 + 1e-4 * randn(rng))]
    @test all(p -> Briareus.closest_hit(closed.bvh, o, onwall(p...) - o)[2] != 0, pts)
    # well inside the polygons (their edge is within ~3 µm of the circles)
    inner = [(ρ * cos(θ), ρ * sin(θ)) for θ in 2π .* rand(rng, 200) for ρ in (0.9 * rand(rng) * d / 2, d / 2 + g / 2)]
    @test all(p -> Briareus.closest_hit(opened.bvh, o, onwall(p...) - o)[2] == 0, inner)

    # --- view factors: in the capsule-free hohlraum (a convex enclosure, so the single
    # visibility ray is exact) rows close as well as without windows, reciprocity is
    # exact by construction, and closure enforcement keeps it
    kw = (; nz = 16, nphi = 24)
    tw = hohlraum(; kw..., capsule_radius = 0,
                  windows = [WallWindow(phi = 78, diameter = d, gap = g, gap_material = gold()),
                             WallWindow(phi = 258, diameter = d, gap = g)])
    t0 = hohlraum(; kw..., capsule_radius = 0)
    F = view_factors(tw)
    err0 = maximum(abs, closure_error(view_factors(t0)))
    @test err0 < 2e-3
    @test maximum(abs, closure_error(F)) < 1.2 * err0
    for n in (:window_1, :window_1_gap, :window_2, :window_2_gap)
        @test maximum(abs, closure_error(F)[elements(tw, n)]) < 1e-3
    end
    G = tw.elem_area .* F
    @test maximum(abs, G - G') <= 1e-14 * maximum(G)
    m = RadiationModel(tw; closed = true)
    @test maximum(abs, closure_error(m.F)) < 1e-10
    Gc = tw.elem_area .* m.F
    @test maximum(abs, Gc - Gc') <= 1e-12 * maximum(Gc)

    # --- a pair related by the 180° rotation about z adds no odd-m harmonics; the open
    # gap's escape is its share of the open area
    q = [remake(b; aim = (0.0, 0.0, sign(b.port[3]) * (b.group in (:inner23, :inner30) ? 3.3 : 5.0)))
         for b in nif_quads()]
    kwc = (; nz = 16, nphi = 24, capsule_ntheta = 12, capsule_nphi = 24)
    capsule_modes(t) = flux_ylm(evaluate(Scenario(t, q; observe = :capsule,
                                                  viewfactor_kw = (; closed = true)), q); lmax = 4)
    oddm(c) = maximum(abs(c[mode_index(l, m)]) for l in 1:4 for m in -l:l if isodd(m)) / c[1]
    c0 = windowless_modes                         # diagnostic windows' (same q, mesh)
    cp = capsule_modes(hohlraum(; kwc..., windows = [WallWindow(phi = 78, diameter = d, open = true),
                                                     WallWindow(phi = 258, diameter = d, open = true)]))
    @test oddm(cp) < 1e-10
    @test 0.9 * c0[1] < cp[1] < 0.97 * c0[1]
    # a closed window of wall material is just wall, up to the remeshing
    cw = capsule_modes(hohlraum(; kwc..., windows = [WallWindow(phi = 0, diameter = d)]))
    @test mode_amplitude(cw, 1) < 1e-3 && abs(cw[1] / c0[1] - 1) < 1e-3
    tg = hohlraum(; kwc..., windows = [WallWindow(phi = 0, diameter = d, gap = g)])
    src = elements(tg, :wall, :cap_top, :cap_bottom)
    laser = zeros(nelements(tg))
    laser[src] .= tg.elem_area[src]
    sol = solve(RadiationModel(tg; closed = true), laser)
    bal = energy_balance(sol)
    @test bal.source ≈ sum(values(bal.absorbed)) + sum(values(bal.escaped)) + bal.lost rtol = 1e-12
    eg = elements(tg, :window_1_gap)
    @test bal.escaped[:window_1_gap] ≈ sum(sol.H[eg] .* tg.elem_area[eg]) rtol = 1e-12
    # the gap sees about the wall's mean flux: its escape is its area times that flux
    ew = elements(tg, :wall)
    Hw = sum(sol.H[ew] .* tg.elem_area[ew]) / sum(tg.elem_area[ew])
    @test bal.escaped[:window_1_gap] / (Hw * π * ((d / 2 + g)^2 - (d / 2)^2)) ≈ 1 atol = 0.05
end
