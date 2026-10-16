# Diagnostic windows in the hohlraum wall: conforming mesh, symmetry, energy balance,
# an integrating-cavity estimate of the drive loss, and laser transmission.

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
    # the radiation areas carry the curved-area quadrature of every cylinder patch
    # (16 flat sub-triangles each: (Δφ/4)²/24 ≈ 5e-5 low on these patches)
    e1 = elements(tw, :window_1)
    @test sum(tw.elem_area[e1]) ≈ w * ht rtol = 1e-4
    @test sum(tw.elem_area[elements(tw, :wall, :window_1, :hole)]) ≈
          sum(hohlraum().elem_area[elements(hohlraum(), :wall)]) rtol = 1e-5
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
