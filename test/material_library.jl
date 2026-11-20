# The tests of the material library cover entries against the package defaults, sources
# and printing.

# Field-by-field equality, recursing into structs (`Material`'s albedo models and
# conversions hold vectors, which `===` compares by identity).
function _same_fields(a, b)
    typeof(a) === typeof(b) || return false
    a isa Union{Number,Symbol,AbstractString,Nothing,AbstractArray,Tuple} && return a == b
    return all(_same_fields(getfield(a, f), getfield(b, f)) for f in fieldnames(typeof(a)))
end

@testset "material library" begin
    lib = material_library()
    names = [d.name for d in lib]
    @test allunique(names)

    # the entries named after a constructor build that constructor's material exactly
    @test _same_fields(material_default(:gold).material, gold())
    @test _same_fields(material_default(:uranium).material, uranium())
    @test _same_fields(material_default(:nif_gold).material, nif_gold())
    @test _same_fields(material_default(:nif_uranium).material, nif_uranium())
    @test _same_fields(material_default(:black).material, BLACK)
    # The foil entries are the default of foil and backlighter_foil and one entry per foil
    # element.
    @test _same_fields(material_default(:gold_foil).material, foil_material())
    @test _same_fields(foil(:f, (0, 0, 0), (0, 0, 1); width = 1.0, n = 1).material, foil_material())
    @test _same_fields(backlighter_foil((90.0, 0.0); standoff = 5.0, n = 1).material, foil_material())
    @test all(_same_fields(material_default(Symbol(n, :_foil)).material, foil_material(e)) for (e, n, _) in FOIL_ELEMENTS)
    @test_throws ArgumentError foil_material(:Kr)
    # the helper sees a difference in a nested field
    @test !_same_fields(gold(), gold(conversion = 0.05))
    @test !_same_fields(nif_gold(), Material(:Au; albedo = WallLossAlbedo(), conversion = nif_conversion(efficiency = 0.71)))

    # every entry's material is built from its sourced values, and every value has a source
    @test all(_same_fields(d.material, Material(d.material.name; albedo = d.albedo.value,
                                                 conversion = d.conversion.value)) for d in lib)
    values = [v for d in lib for v in (d.albedo, d.conversion, last.(d.other)...)]
    @test all(!isempty(strip(v.source)) && !isempty(strip(v.uncertainty)) for v in values)
    @test all(v.basis in MATERIAL_BASES for v in values)

    # The published values are checked: gold() is the static model of Young 2020, and
    # the alternatives follow.
    y, s, disk = material_default(:gold), material_default(:gold_swadling2025), material_default(:gold_disk)
    @test (albedo(y.material), conversion_efficiency(y.material.conversion)) == (0.8, 0.9)
    @test y.albedo.basis === y.conversion.basis === :published
    @test material_default(:uranium).conversion.basis === :published
    @test (albedo(s.material), conversion_efficiency(s.material.conversion)) == (0.77, 0.95)
    @test conversion_efficiency(disk.material.conversion, 0.0, 1e14) ≈ 0.95 rtol = 1e-12
    @test conversion_efficiency(disk.material.conversion, 0.0, 1e15) ≈ 0.80 rtol = 1e-12
    @test first.(disk.other) == [:absorption, :cold_reflectance_1w]
    # the gold foil is the disk table's value at its 1e15 W/cm² end, held beyond it
    gf = material_default(:gold_foil)
    @test conversion_efficiency(gf.material.conversion) == conversion_efficiency(disk.material.conversion, 0.0, 1e15)
    @test conversion_efficiency(disk.material.conversion, 0.0, 3e15) ≈ gf.conversion.value rtol = 1e-12
    # Mid-Z foils have suggested values, each the midpoint of the range its uncertainty
    # states.
    midz = [material_default(Symbol(n, :_foil)) for (_, n, _) in FOIL_ELEMENTS[2:end]]
    @test all(d.albedo.basis === d.conversion.basis === :suggested for d in midz)
    @test all(d.conversion.value ≈ sum(Briareus._MIDZ_CONVERSION_RANGE) / 2 &&
              d.albedo.value ≈ sum(Briareus._MIDZ_ALBEDO_RANGE) / 2 for d in midz)
    @test all(occursin("suggested range 0.1-0.6", d.conversion.uncertainty) &&
              occursin("not a measurement", d.conversion.uncertainty) for d in midz)
    @test gf.conversion.basis === :published && !any(d -> d.conversion.basis === :suggested,
                                                       filter(d -> !endswith(String(d.name), "_foil") || d.name === :gold_foil, lib))

    # defaults without a source are left out
    @test !(:ablator in names)

    # Printing gives one source line and one uncertainty line per quantity.
    for d in (material_default(:gold), disk)
        txt = sprint(show, MIME"text/plain"(), d)
        n = 2 + length(d.other)
        @test count("source: ", txt) == n && count("uncertainty: ", txt) == n
        @test occursin(String(d.name), txt)
    end
    @test sprint(show, lib) == "MaterialDefault[" * join(("MaterialDefault(:$n)" for n in names), ", ") * "]"
    @test occursin("WallLossAlbedo(1 - α = 0.32", sprint(show, MIME"text/plain"(), material_default(:nif_gold)))
    @test occursin("x 0.836 (T/300 eV)^-0.2", sprint(show, MIME"text/plain"(), material_default(:uranium)))

    @test_throws ArgumentError material_default(:unobtainium)
    @test_throws ArgumentError SourcedValue(0.5, :published, " ", "none")
    @test_throws ArgumentError SourcedValue(0.5, :guess, "a source", "none")
end
