# Diagnostic definitions: the built-in Dante 1 against the hardware table it replaced
# (bit for bit), extends and shadowing on a search path, measured CSV curves against the
# models they tabulate, the schema checks, the configuration file and a registered kind.
# Design: reference/diagnostics-design.md.

using Briareus, Test
using Briareus: _DANTE_GRID, _XRD_FIT

# NIF Dante 1 as `src/dante_response.jl` held it before the definitions (Rubery 2023
# Fig. 4): number, energy, filters (µm), grid, mirror, cathode, solid angle
const OLD_DANTE1 = [
    (1, 1170, [:Mg => 22], false, nothing, :Al, 6.51e-7),
    (2, 3060, [:polyimide => 47], false, (:mlm, 1.9), :Al, 1.90e-7),
    (3, 250, [:lexan => 4], true, (:SiO2, 3.5), :Al, 5.74e-7),
    (4, 470, [:V => 1.75], false, (:SiO2, 2.5), :Al, 5.56e-7),
    (5, 710, [:Co => 1.4, :CH => 0.2], true, nothing, :Al, 6.51e-7),
    (6, 870, [:Cu => 2], true, nothing, :Cr, 6.51e-7),
    (7, 940, [:Mg => 5, :Zn => 0.65, :CH => 0.1], true, nothing, :Ni, 6.51e-7),
    (8, 1170, [:Mg => 22], true, nothing, :Al, 6.51e-7),
    (9, 1430, [:Al => 20], true, nothing, :Ni, 6.51e-7),
    (10, 1710, [:Si => 20], false, nothing, :Al, 6.51e-7),
    (11, 3430, [:Fe => 1, :Cr => 0.65, :parylene => 5], true, nothing, :Al, 6.51e-7),
    (12, 2580, [:saran => 66], false, nothing, :Al, 6.51e-7),
    (13, 3000, [:Ag => 4.4], false, nothing, :Al, 6.51e-7),
    (14, 1700, [:Si => 20], false, nothing, :Al, 6.51e-7),
    (15, 940, [:Mg => 5, :Zn => 0.65, :CH => 0.1], true, nothing, :Ni, 6.51e-7),
    (16, 14370, [:Al => 250], false, nothing, :Ni, 1.77e-6),
    (17, 690, [:Co => 1], true, nothing, :Cr, 6.51e-7),
    (18, 3660, [:Fe => 1, :Cr => 0.65, :parylene => 5], false, nothing, :Ni, 1.77e-6)]

old_dante1() = [DanteResponse(n, E; filters = [k => 1e-3t for (k, t) in f], grid = g ? 0.088 : 1.0,
                              mirror = m, cathode = x, solid_angle = Ω)
                for (n, E, f, g, m, x, Ω) in OLD_DANTE1]

# a definition file in `dir`
put(dir, name, text) = (write(joinpath(dir, name), text); joinpath(dir, name))

# a whole xray_power definition: header, then `body` (channels), then a line of sight
defn(id, body) = "schema = 1\nid = \"$id\"\nkind = \"xray_power\"\n" * body *
                 "\n[line_of_sight]\ntheta_deg = 90.0\nphi_deg = 0.0\n"
# one channel with a Cr XRD, plus `extra` lines
chan(extra = "") = "[[channel]]\nnumber = 1\nenergy_eV = 500\nsolid_angle_sr = 1e-6\n" * extra
csv(dir, name, E, v) = (open(joinpath(dir, name), "w") do io
                            println(io, "energy_eV,value")
                            foreach((e, x) -> println(io, e, ",", x), E, v)
                        end; name)

@testset "diagnostic definitions" begin
    @testset "built-in Dante 1, bit for bit" begin
        new, old = dante1_responses(), old_dante1()
        @test [c.number for c in new] == 1:18
        @test all(new[i].R == old[i].R for i in 1:18)
        @test all(new[i].energy === old[i].energy && new[i].solid_angle === old[i].solid_angle &&
                  new[i].filters == old[i].filters && new[i].grid === old[i].grid &&
                  new[i].mirror == old[i].mirror && new[i].cathode === old[i].cathode for i in 1:18)
        @test dante1_responses() === dante1_responses()            # built once
        v = dante1()
        @test (v.theta, v.phi, v.through, v.pixels, v.area) == (143.0, 274.0, :leh_bottom, 400, nothing)
        @test all(v.channels .=== dante1_responses())
        @test dante1(pixels = 50).pixels == 50
        @test dante1(channels = dante1_channels(model = :gaussian)).channels[1] isa DanteChannel
        # every built-in builds without issues
        @test isempty(check_diagnostics(; paths = diagnostic_paths(builtin = true)[1:1]))
        defs = diagnostics(; paths = diagnostic_paths()[1:1])
        @test defs["nif_dante1"].chain == [defs["nif_dante1_2021"].file, defs["nif_dante1"].file]
        @test defs["nif_dante1"].name == "NIF Dante 1" && defs["nif_dante1"].public
    end

    @testset "extends, shadowing, CSV curves" begin
        mktempdir() do dir
            builtin = diagnostic_paths()[1]
            a, b = mkpath(joinpath(dir, "a")), mkpath(joinpath(dir, "b"))
            # a campaign's channel 6 (Cu 1.1 µm + grid, as N170924-003), channel 16 removed,
            # a measured channel 19 added
            E = [100.0, 1000.0, 10_000.0, 20_000.0]
            csv(a, "ch19.csv", E, [0.0, 2.0e6, 1.0e6, 0.0])
            put(a, "campaign.toml", """
                schema = 1
                id = "group_dante1"
                extends = "nif_dante1"
                public = false
                remove_channels = [16]
                [[channel]]
                number = 6
                filters = [{ material = "Cu", thickness_um = 1.1 }]
                [[channel]]
                number = 19
                energy_eV = 1000
                solid_angle_sr = 6.5e-7
                response_csv = "ch19.csv"
                """)
            d = diagnostic("group_dante1"; paths = [builtin, a])
            base = dante1_responses()
            @test [c.number for c in d.channels] == [1:15; 17:19]
            ref6 = DanteResponse(6, 870; filters = [:Cu => 1e-3 * 1.1], grid = 0.088, cathode = :Cr,
                                 solid_angle = 6.51e-7)
            @test d.channels[6].R == ref6.R                        # filters replaced, the rest kept
            @test all(d.channels[i].R == base[i].R for i in [1:5; 7:15])
            r19 = d.channels[end]
            @test r19.response isa TabulatedCurve
            @test r19.R[1] == 0.0 && r19.R[end] == 0.0           # zero outside the table
            k = findfirst(>=(1000.0), _DANTE_GRID)
            e = _DANTE_GRID[k]
            @test r19.R[k] ≈ 2.0e6 + (e - 1000) / 9000 * (1.0e6 - 2.0e6) rtol = 1e-12
            defs = diagnostics(; paths = [builtin, a])
            @test !defs["group_dante1"].public
            @test length(defs["group_dante1"].chain) == 3
            # a later definition of an id shadows an earlier one, with a warning naming both
            put(b, "dante1.toml", defn("nif_dante1", chan("cathode = \"Cr\"\n")))
            issues = check_diagnostics(; paths = [builtin, a, b])
            @test any(i -> i.level === :warning && occursin("shadows", i.message), issues)
            @test length(diagnostic("nif_dante1"; paths = [builtin, b]).channels) == 1
            # group_dante1 now extends the one-channel nif_dante1: its partial channel 6
            # lacks the values the built-in gave it, reported and not built
            @test any(i -> basename(i.file) == "campaign.toml" && i.key == "channel[6].energy_eV" &&
                      i.level === :missing, issues)
            @test_throws ArgumentError diagnostic("group_dante1"; paths = [builtin, a, b])
            @test length(dante1().channels) == 18                  # the built-in one is unchanged
            # measured factors that tabulate the models on the grid reproduce them
            g = _DANTE_GRID
            c1 = dante1_responses()[1]
            T = [filter_transmission([:Mg => 0.022], e / 1000; table = :henke) for e in g]
            X = [xrd_sensitivity(:Al, e) for e in g]
            put(a, "tab.toml", defn("tabulated", """
                [[channel]]
                number = 1
                energy_eV = 1170
                solid_angle_sr = 6.51e-7
                filters = [{ transmission_csv = "$(csv(a, "mg.csv", g, T))" }]
                sensitivity_csv = "$(csv(a, "al.csv", g, X))"
                """))
            t1 = diagnostic("tabulated"; paths = [a]).channels[1]
            @test t1.R ≈ c1.R rtol = 1e-12
            @test t1.cathode isa TabulatedCurve && length(t1.tables) == 1
            @test isempty(filter(i -> i.level !== :warning, check_diagnostics(; paths = [builtin, a])))
        end
    end

    @testset "schema checks" begin
        mktempdir() do dir
            put(dir, "typo.toml", defn("typo", chan("cathode = \"Cr\"\nthicknes_um = 3\n")))
            put(dir, "types.toml", defn("types", replace(chan("cathode = \"Cr\"\n"), "500" => "\"500\"")))
            put(dir, "material.toml", defn("material",
                chan("cathode = \"Cr\"\nfilters = [{ material = \"Unobtainium\", thickness_um = 1 }]\n")))
            put(dir, "nocathode.toml", defn("nocathode", chan()))
            put(dir, "ok.toml", defn("ok", chan("cathode = \"Cr\"\n")))
            put(dir, "nochannels.toml", "schema = 1\nid = \"nochannels\"\nextends = \"ok\"\nremove_channels = [1]\n")
            put(dir, "schema2.toml", "schema = 2\nid = \"future\"\nkind = \"xray_power\"\n")
            put(dir, "cycle.toml", "schema = 1\nid = \"x\"\nextends = \"y\"\n")
            put(dir, "cycle2.toml", "schema = 1\nid = \"y\"\nextends = \"x\"\n")
            put(dir, "kind.toml", "schema = 1\nid = \"k\"\nkind = \"spectrometer\"\n")
            put(dir, "csv.toml", defn("csv", chan("sensitivity_csv = \"missing.csv\"\n")))
            issues = check_diagnostics(; paths = [dir])
            has(file, key, level) = any(i -> basename(i.file) == file && occursin(key, i.key) && i.level === level,
                                        issues)
            @test has("typo.toml", "channel[1].thicknes_um", :error)
            @test has("types.toml", "channel[1].energy_eV", :error)
            @test has("material.toml", "channel[1].filters[1].material", :error)
            @test has("nocathode.toml", "channel[1].cathode", :missing)
            @test has("nochannels.toml", "channel", :missing)
            @test !any(i -> basename(i.file) == "ok.toml", issues)
            @test has("schema2.toml", "schema", :error)
            @test has("cycle.toml", "extends", :error) || has("cycle2.toml", "extends", :error)
            @test has("kind.toml", "kind", :error)
            @test has("csv.toml", "channel[1].sensitivity_csv", :error)
            @test_throws ArgumentError diagnostic("typo"; paths = [dir])
            @test_throws ArgumentError diagnostic("absent"; paths = [dir])
        end
    end

    @testset "configuration file and search path" begin
        mktempdir() do dir
            cfg = joinpath(dir, "cfg", "config.toml")
            g1, g2, g3 = (mkpath(joinpath(dir, s)) for s in ("g1", "g2", "g3"))
            withenv("BRIAREUS_CONFIG" => cfg, "BRIAREUS_DIAGNOSTICS" => g3) do
                @test diagnostic_config_path() == cfg
                p = add_diagnostic_path(g1)
                @test isfile(cfg)
                @test p == [diagnostic_paths()[1], normpath(g1), normpath(g3)]
                add_diagnostic_path(g1)                                 # once only
                @test count(==(normpath(g1)), diagnostic_paths()) == 1
                add_diagnostic_path(g2; save = false)
                @test diagnostic_paths()[end] == normpath(g2)
                @test !occursin("g2", read(cfg, String))
                remove_diagnostic_path(g2; save = false)
                remove_diagnostic_path(g1)
                @test diagnostic_paths() == [diagnostic_paths()[1], normpath(g3)]
                @test diagnostic_paths(builtin = false) == [normpath(g3)]
                @test_throws ArgumentError add_diagnostic_path(joinpath(dir, "none"))
            end
        end
    end

    @testset "a registered kind" begin
        spec = Dict("gain" => DiagnosticField(Real; required = true))
        register_diagnostic_kind(:test_counter, spec, (d; scale = 1) -> scale * d.data["gain"])
        mktempdir() do dir
            put(dir, "c.toml", "schema = 1\nid = \"counter\"\nkind = \"test_counter\"\ngain = 3.0\n")
            @test diagnostic("counter"; paths = [dir]) == 3.0
            @test diagnostic("counter"; paths = [dir], scale = 2) == 6.0
            put(dir, "d.toml", "schema = 1\nid = \"counter2\"\nextends = \"counter\"\ngain = 5.0\n")
            @test diagnostic("counter2"; paths = [dir]) == 5.0
        end
        delete!(Briareus._DIAGNOSTIC_KINDS, :test_counter)
    end
end
