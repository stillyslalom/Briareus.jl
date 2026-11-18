# A registry of named material defaults, each quantity with its data source and
# uncertainty (research notes, "Material defaults and their sources (2026-11-18)").

"""
    SourcedValue(value, basis, source, uncertainty)

One quantity of a [`MaterialDefault`](@ref) with its provenance: `value` (a number, an
`AlbedoModel` or a [`Conversion`](@ref)), `basis` (how the value was obtained, one of
[`MATERIAL_BASES`](@ref)), `source` (the citation, and the research-notes section that
records it) and `uncertainty` (the stated or estimated spread, or the published values
around it). `source` must not be empty: a value without a source does not enter the
library.
"""
struct SourcedValue
    value::Any
    basis::Symbol
    source::String
    uncertainty::String
    function SourcedValue(value, basis::Symbol, source::AbstractString, uncertainty::AbstractString)
        basis in MATERIAL_BASES || throw(ArgumentError("basis must be one of $MATERIAL_BASES, got :$basis"))
        isempty(strip(source)) && throw(ArgumentError("a sourced value needs a source"))
        return new(value, basis, String(source), String(uncertainty))
    end
end

"""
    MATERIAL_BASES

How a [`SourcedValue`](@ref) was obtained:
- `:published`: a measured value, a published scaling, or the input of a published model,
  as printed;
- `:range`: a round value inside a published range, not itself a published number;
- `:calibrated`: fitted in Briareus to the published data the source names;
- `:derived`: Briareus's own arithmetic on published inputs;
- `:definition`: true by construction (a perfect absorber).
"""
const MATERIAL_BASES = (:published, :range, :calibrated, :derived, :definition)

"""
    MaterialDefault

A named entry of [`material_library`](@ref): `name`, a one-line `description`, the
`material` it builds (a `Material`), the [`SourcedValue`](@ref)s `albedo` and
`conversion` the material is built from, and `other`, recorded quantities the radiation
model does not use (laser absorption, cold reflectance; `quantity => SourcedValue`
pairs). Printed with `show` (`display(entry)`), it lists every quantity with its basis,
source and uncertainty.
"""
struct MaterialDefault
    name::Symbol
    description::String
    material::Material
    albedo::SourcedValue
    conversion::SourcedValue
    other::Vector{Pair{Symbol,SourcedValue}}
end

function MaterialDefault(name::Symbol, description, material_name::Symbol, albedo::SourcedValue,
                         conversion::SourcedValue, other = Pair{Symbol,SourcedValue}[])
    m = Material(material_name; albedo = albedo.value, conversion = conversion.value)
    return MaterialDefault(name, String(description), m, albedo, conversion, collect(Pair{Symbol,SourcedValue}, other))
end

# Shared provenance strings. Each names the citation and the research-notes section that
# records the number.
const _SRC_YOUNG_STATIC = "Young et al., Phys. Plasmas 27, 082702 (2020): the static VisRad model of " *
    "N170827 (albedo 0.8 Au, 0.3 ablator, conversion 0.9, radii at 85 %); research notes, " *
    "\"Young 2020 validation (2026-10-24)\" and \"Literature review against Kritcher et al. 2024 (2026-10-22)\""
const _UNC_GOLD_ALBEDO = "published static and 0-D wall albedos 0.77 (Swadling 2025, VisRad calibrated to " *
    "Dante at 400 TW) to 0.84 (Swadling 2025 0-D, 288 eV); Lindl's loss gives 0.88-0.94 at a NIF peak; " *
    "a single static value stands for an albedo that rises in time"
const _SRC_CONV_075 = "Lindl et al., Phys. Plasmas 11, 339 (2004) Sec. IV (about 0.7 in 1 ns Nova " *
    "hohlraums, up to 0.85 for longer shaped pulses); Rosen, OSTI 885380 (2005) (0.68 at 1 ns, Omega); " *
    "research notes, \"Energetics (energetics.md)\" under \"Literature review against Kritcher et al. 2024\" " *
    "and \"Capsule absorption, conversion and the ablating-wall release (2026-10-25)\""
const _UNC_CONV_075 = "published hohlraum values 0.68-0.85 at 1 ns (Lindl, Rosen, Schein 0.85 t^0.2); " *
    "0.85-0.95 at a NIF peak (Olson 2012, Swadling 2025): 0.75 is low for NIF"
const _SRC_LINDL_LOSS = "Lindl et al., Phys. Plasmas 11, 339 (2004), eq. 4-31: 1 - α = 0.32 T^-0.7 t^-0.38 " *
    "(T in heV, t in ns); re-emission `release` 1.0 calibrated against 1-D ablating-wall " *
    "radiation hydrodynamics; research notes, \"Gold wall loss and thin-foil burnthrough (2026-10-16)\" " *
    "and \"Capsule absorption, conversion and the ablating-wall release (2026-10-25)\" sec. 1"
const _UNC_LINDL_LOSS = "Hammer & Rosen 2003 and Shussman & Heizler 2015 losses within about 10 % of " *
    "Lindl's over 1-3 heV and 1-10 ns (derived); release 1.0, range 0.84-1.07"
const _SRC_DU = "du_wall_loss(): Rosen, OSTI 885380 (2005) sec. VI and Schein et al., PRL 98, 175003 " *
    "(2007) for the T^-0.2 scaling; factor 0.836 at 300 eV calibrated to the 6.5 % peak-power " *
    "equivalence of Döppner et al., PRL 115, 055001 (2015); 0.7 µm gold liner from the same paper; " *
    "research notes, \"Energetics defaults and the capsule share (2026-10-25)\" and \"Wall material, " *
    "conversion and capsule energetics (2026-10-24)\" sec. 1"
const _UNC_DU = "factor 0.827-0.838 over effective wall times 2-4 ns; independent cocktail ratios " *
    "0.83 (Rosen) and 0.825 (Schein); check: +4.8 eV at equal power against the measured +5 ± 6 eV"
const _SRC_NIF_CONV = "nif_conversion(): Lindl 2004 Fig. 4-4 caption (0.7 t^0.12), capped at 0.95 " *
    "(Dewald et al., PoP 15, 072706 (2008) plateau at 1e14 W/cm²; Swadling et al., PoP 32, 052707 " *
    "(2025)); research notes, \"Wall material, conversion and capsule energetics (2026-10-24)\" sec. 2"
const _UNC_NIF_CONV = "0.80/0.85/0.90 at 3/5/8 ns against Olson et al., PoP 19, 053301 (2012): " *
    "0.85-0.90 in NIF vacuum hohlraums at 235-345 eV; published time exponents 0.12-0.2"

"""
    material_library() -> Vector{MaterialDefault}

The named material defaults whose every quantity has a published source or a
derivation recorded in the research notes ("Material defaults and their sources
(2026-11-18)"), freshly built on each call. Each entry is a [`MaterialDefault`](@ref);
`material_default(name)` returns one, and its `material` is the `Material` to use. The
entries named after a constructor build exactly the material that constructor builds
today:

- `:gold`: [`gold`](@ref)`()`, albedo 0.8 and conversion 0.75 (static solves);
- `:uranium`: [`uranium`](@ref)`()`, the Au-lined DU wall, conversion 0.75;
- `:nif_gold`, `:nif_uranium`: [`nif_gold`](@ref)`()`, [`nif_uranium`](@ref)`()`, the
  drive-history walls (Lindl's wall loss, the rising `nif_conversion()`);
- `:black`: [`BLACK`](@ref), the perfect absorber.

The others are published alternatives:
- `:gold_young2020`: the gold wall of a published static VisRad model of a NIF hohlraum
  (albedo 0.8, conversion 0.9; Young et al., PoP 27, 082702 (2020));
- `:gold_swadling2025`: the only published VisRad energetics calibration (albedo 0.77,
  conversion 0.95 on a 5.75 mm gold hohlraum at 400 TW; Swadling et al., PoP 32, 052707
  (2025));
- `:gold_disk`: directly irradiated gold in open geometry (a foil or disk): Lindl's
  wall loss and Dewald et al.'s (PoP 15, 072706 (2008)) conversion plateau, 0.95 at
  10¹⁴ and 0.80 at 10¹⁵ W/cm² (an intensity table), with the 3ω laser absorption and
  the 1ω cold reflectance recorded.

Package defaults without a source are not in the library: the capsule's laser → x-ray
conversion 0.3 ([`ablator`](@ref)), the backlighter foil's conversion 0.05
([`foil`](@ref), [`backlighter_foil`](@ref)) and `Material`'s own 0/0. `ablator`'s
albedo 0.3 has a source (Young 2020's static model) but its conversion has none, so the
material is left out whole.
"""
function material_library()
    young_albedo = SourcedValue(0.8, :published, _SRC_YOUNG_STATIC, _UNC_GOLD_ALBEDO)
    conv075 = SourcedValue(0.75, :range, _SRC_CONV_075, _UNC_CONV_075)
    lindl = SourcedValue(WallLossAlbedo(), :published, _SRC_LINDL_LOSS, _UNC_LINDL_LOSS)
    du = SourcedValue(du_wall_loss(), :calibrated, _SRC_DU, _UNC_DU)
    nifconv = SourcedValue(nif_conversion(), :derived, _SRC_NIF_CONV, _UNC_NIF_CONV)
    swadling = "Swadling et al., Phys. Plasmas 32, 052707 (2025), Fig. 7 caption: VisRad of a " *
               "5.75 mm Au hohlraum at a constant 400 TW, conversion set to 0.95 and the wall albedo " *
               "adjusted to 0.77 to match the Dante 1 and 2 peak intensity; research notes, \"Energetics " *
               "(energetics.md)\" under \"Literature review against Kritcher et al. 2024 (2026-10-22)\""
    black = SourcedValue(0.0, :definition, "a perfect absorber that emits nothing (by construction)", "none")
    return [
        MaterialDefault(:gold, "gold hohlraum wall for static solves (gold())", :Au, young_albedo, conv075),
        MaterialDefault(:uranium, "Au-lined depleted-uranium hohlraum wall for static solves (uranium())", :DU,
                        du, conv075),
        MaterialDefault(:nif_gold, "gold wall of a NIF drive history (nif_gold())", :Au, lindl, nifconv),
        MaterialDefault(:nif_uranium, "Au-lined DU wall of a NIF drive history (nif_uranium())", :DU, du, nifconv),
        MaterialDefault(:gold_young2020, "gold wall of a published static VisRad model of a NIF hohlraum", :Au,
                        young_albedo,
                        SourcedValue(0.9, :published, _SRC_YOUNG_STATIC,
                                     "agrees with 0.85-0.90 in NIF vacuum hohlraums (Olson 2012)")),
        MaterialDefault(:gold_swadling2025, "gold wall of the published VisRad energetics calibration", :Au,
                        SourcedValue(0.77, :published, swadling,
                                     "a fitted value: an albedo error of 0.02 is a 10 % capsule-drive error " *
                                     "in the same paper's 0-D model, which gives 0.84 at 288 eV"),
                        SourcedValue(0.95, :published, swadling,
                                     "set, not fitted; at the 1e14 W/cm² plateau of Dewald 2008")),
        MaterialDefault(:gold_disk, "directly irradiated gold in open geometry (foil, disk)", :Au, lindl,
                        SourcedValue(Conversion(intensities = (1e14, 1e15), efficiencies = (0.95, 0.80)),
                                     :published,
                                     "Dewald et al., Phys. Plasmas 15, 072706 (2008), abstract: the " *
                                     "backscatter-corrected plateau after ~0.5 ns on Au, U and cocktail " *
                                     "spheres; research notes, \"Wall material, conversion and capsule " *
                                     "energetics (2026-10-24)\" (Inputs read)",
                                     "Mead et al., LA-UR-86-300 (1985) Fig. 16: 0.51-0.59 of the absorbed " *
                                     "energy on open 3ω disks at 2e13-4e14 W/cm² (0.5-1 ns, older " *
                                     "small spots); Dewald's values read from the abstract only (confidence M)"),
                        [:absorption => SourcedValue((intensities = (3e13, 1.5e14, 3e14), values = (0.93, 0.95, 0.91)),
                                                     :published,
                                                     "Ze et al., UCRL-86414 (1981), Fig. 4, 0.35 µm on gold disks; " *
                                                     "research notes, \"Unconverted light: absorption and " *
                                                     "conversion vs. wavelength (2026-10-28)\"",
                                                     "read from a figure, two of three marked approximate; Briareus deposits " *
                                                     "all incident light (absorption 1) and books backscatter " *
                                                     "separately"),
                         :cold_reflectance_1w => SourcedValue(0.98, :published,
                                                              "Johnson & Christy / Palik optical constants at 1.06 µm, " *
                                                              "normal incidence; research notes, \"Unconverted light: " *
                                                              "absorption and conversion vs. wavelength (2026-10-28)\"",
                                                              "approximate, confidence M; holds below plasma onset " *
                                                              "(about 1e8 W/cm², NIF User Guide sec. 4.4)")]),
        MaterialDefault(:black, "perfect absorber (BLACK; apertures)", :black, black, black),
    ]
end

"""
    material_default(name::Symbol) -> MaterialDefault

The entry `name` of [`material_library`](@ref); its `material` field is the `Material`.
Throws an `ArgumentError` listing the names for an unknown one.
"""
function material_default(name::Symbol)
    lib = material_library()
    i = findfirst(d -> d.name === name, lib)
    i === nothing && throw(ArgumentError("no material default :$name; the library has " *
                                         join((":" * String(d.name) for d in lib), ", ")))
    return lib[i]
end

_value_summary(x::Real) = string(x)
_value_summary(x::ConstantAlbedo) = string(x.value)
_value_summary(x::WallLossAlbedo) =
    @sprintf("WallLossAlbedo(1 - α = %g T^%g t^%g%s)", x.coefficient, x.temperature_exponent, x.time_exponent,
             x.factor == 1 && x.factor_exponent == 0 ? "" :
             @sprintf(" x %g (T/%g eV)^%g, liner %.3g g/cm²", x.factor, x.factor_temperature, x.factor_exponent, x.liner))
function _value_summary(c::Conversion)
    η0 = isempty(c.table) ? string(c.efficiency) :
         "table " * join((@sprintf("%.3g at %.0e W/cm²", v, exp10(l)) for (l, v) in zip(c.log_intensity, c.table)), ", ")
    s = c.time_exponent == 0 ? η0 : "$η0 x t^$(c.time_exponent)"
    return c.range == (0.0, 1.0) ? "Conversion($s)" : "Conversion($s, range $(c.range))"
end
_value_summary(x::NamedTuple) = join(("$(v) at $(@sprintf("%.1e", i)) W/cm²" for (i, v) in zip(x.intensities, x.values)), ", ")
_value_summary(x) = repr(x)

Base.show(io::IO, d::MaterialDefault) = print(io, "MaterialDefault(:", d.name, ")")

function Base.show(io::IO, ::MIME"text/plain", d::MaterialDefault)
    println(io, "MaterialDefault :", d.name, " — ", d.description, " (Material :", d.material.name, ")")
    rows = [:albedo => d.albedo, :conversion => d.conversion, d.other...]
    for (k, (q, v)) in enumerate(rows)
        println(io, "  ", q, ": ", _value_summary(v.value), " [", v.basis, "]")
        println(io, "    source: ", v.source)
        print(io, "    uncertainty: ", v.uncertainty)
        k < length(rows) && println(io)
    end
end
