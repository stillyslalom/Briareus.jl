# This file defines the starters, which insert assemblies of parts.

# ---------------------------------------------------------------- starters

"""
Starter assemblies, inserted with [`insert_starter!`](@ref): `:hohlraum` (*Hohlraum with
capsule*: a hohlraum part and a capsule sphere with the observed role), `:foil`
(*Backlighter foil*: a foil rectangle facing TCC and an absorbing sphere standing in for
the main target), `:shock_tube` (*Halfraum pair with shock tube*: two halfraums facing
each other across a cylinder package whose end face has the observed role) and
`:opacity` (*Opacity halfraum*: a halfraum with a sample foil, observed, at its open
end) and `:frustraum` (*Frustraum with capsule*: a hohlraum part with a frustraum
profile and a capsule). Their dimensions are illustrative (the NIF-scale defaults of the
former templates; mock-ups for the halfraums; the frustraum at a published platform's
stated size); none is a facility design.
"""
const STARTERS = (:hohlraum, :foil, :shock_tube, :opacity, :frustraum)

const STARTER_LABELS = Dict(:hohlraum => "Hohlraum with capsule", :foil => "Backlighter foil",
                            :shock_tube => "Halfraum pair with shock tube", :opacity => "Opacity halfraum",
                            :frustraum => "Frustraum with capsule")

"""
    starter_parts(kind; frame = :tcc, taken = (), capsule_h = 0.3) -> Vector{PartSpec}

The parts of starter `kind` (one of [`STARTERS`](@ref)), placed in the frame named
`frame`, their names made unique against `taken` (`hohlraum`, then `hohlraum_2`, …).

- `:hohlraum`: a `:hohlraum` part `hohlraum` (length 10.1 mm, radius 2.875 mm, entrance
  holes of 1.55 mm radius, wall albedo 0.8 and conversion 0.9, core `gold()`) and,
  placed in its frame at its center, the sphere `capsule` (radius 0.9 mm, albedo and
  conversion 0.3: the ablator's, meshed as a capsule, `capsule_mesh`) with the observed
  role. At TCC it meshes as the former hohlraum template did, at every element.
- `:foil`: the rectangle `foil` (2 mm square, albedo and conversion 0.8, core
  `foil_material()`, a gold foil; its element size the target's capsule size
  `capsule_h` at the insert), 6 mm from the frame's origin
  along its `x` and facing it (its edges along the frame's `z` and `y`), with the foil
  role, and the absorbing sphere `main` (radius 1 mm, albedo and conversion 0) at the
  origin. At TCC both mesh as the former foil template's did.
- `:shock_tube`: two halfraums (hohlraum parts `halfraum_top` and `halfraum_bottom`,
  each 2.5 mm long, radius 1.5 mm, of the hohlraum part's wall), each with an entrance
  hole of 1 mm radius at its outer end and its inner end open, facing each other across
  the package `package` (a cylinder of radius 1.5 mm and length 2 mm at the frame's
  origin, along its `z`, absorbing: albedo and conversion 0) whose end faces, the disks
  `package_top` (with the observed role) and `package_bottom` (radius 1.5 mm, albedo and
  conversion 0.3, placed in the package's frame), close each halfraum's open end.
- `:opacity`: a halfraum (the hohlraum part `halfraum`, 3 mm long, radius 1.5 mm, an
  entrance hole of 1 mm radius at its top, its bottom open) and the sample `sample` (a
  1 mm square rectangle, albedo 0.3, conversion 0, with the observed role) in its open
  end's plane, facing into it; the radiation that misses the sample leaves through the
  open end.
- `:frustraum`: a `:hohlraum` part `frustraum` with the frustraum profile (length
  12.24 mm, waist radius 4.88 mm, end radius 2.282 mm, entrance holes of 1.55 mm radius)
  and, at its center, the capsule sphere `capsule` (radius 1.281 mm, as the hohlraum
  starter's capsule otherwise) with the observed role. The length, waist and holes are
  the stated dimensions of the NIF frustraum of Baker et al., Phys. Plasmas 30, 092708
  (2023); the end radius is derived from its 23° wall angle taken from the axis, and the
  capsule is its 1.2 mm inner radius plus the ~81 µm shell. The starter is an
  illustrative mock-up, and its wall area differs from the one the paper quotes. The NIF
  presets aim in its holes with a pointing derived for it ([`preset_pointing`](@ref));
  the pointing of the platform is not used.

The two halfraum starters are mock-ups at an illustrative scale; reference/research-notes.md,
"Halfraum ends and starters" and "Hohlraum profiles", has the sources.
"""
function starter_parts(kind::Symbol; frame::Symbol = :tcc, taken = (), capsule_h::Real = 0.3)
    kind in STARTERS || throw(ArgumentError("a starter is one of $STARTERS"))
    used = Set{Symbol}(taken)
    function unique_name(base::Symbol)
        name, k = base, 1
        while name in used || any(n -> startswith(String(n), String(name) * "."), used)
            k += 1
            name = Symbol(base, :_, k)
        end
        push!(used, name)
        return name
    end
    if kind === :hohlraum
        h = unique_name(:hohlraum)
        c = unique_name(:capsule)
        return [PartSpec(:hohlraum; name = h, frame),
                PartSpec(:sphere; name = c, frame = h, observed = true,
                         values = Dict(:r1 => 0.9, :albedo => 0.3, :conversion => 0.3, :capsule_mesh => 1.0))]
    end
    if kind === :frustraum
        h = unique_name(:frustraum)
        c = unique_name(:capsule)
        return [PartSpec(:hohlraum; name = h, frame,
                         values = Dict(:profile => 1.0, :length => 12.24, :radius => 4.88, :end_radius => 2.282,
                                       :leh_radius => 1.55)),
                PartSpec(:sphere; name = c, frame = h, observed = true,
                         values = Dict(:r1 => 1.281, :albedo => 0.3, :conversion => 0.3, :capsule_mesh => 1.0))]
    end
    if kind === :shock_tube
        a, b = unique_name(:halfraum_top), unique_name(:halfraum_bottom)
        t, ft, fb = unique_name(:package), unique_name(:package_top), unique_name(:package_bottom)
        R, Lh, Lp = 1.5, 2.5, 2.0
        half(name, z, top, bottom) =
            PartSpec(:hohlraum; name, frame, values = Dict(:z => z, :length => Lh, :radius => R, :leh_radius => 1.0,
                                                           :top_end => top, :bottom_end => bottom))
        face(name, z, theta, observed) =
            PartSpec(:disk; name, frame = t, observed,
                     values = Dict(:z => z, :theta => theta, :r1 => R, :albedo => 0.3, :conversion => 0.3))
        return [half(a, (Lp + Lh) / 2, 1.0, 2.0), half(b, -(Lp + Lh) / 2, 2.0, 1.0),
                PartSpec(:cylinder; name = t, frame, values = Dict(:r1 => R, :length => Lp)),
                face(ft, Lp / 2, 0.0, true), face(fb, -Lp / 2, 180.0, false)]
    elseif kind === :opacity
        h, s = unique_name(:halfraum), unique_name(:sample)
        return [PartSpec(:hohlraum; name = h, frame,
                         values = Dict(:length => 3.0, :radius => 1.5, :leh_radius => 1.0, :bottom_end => 2.0)),
                PartSpec(:rectangle; name = s, frame = h, observed = true,
                         values = Dict(:z => -1.5, :r1 => 1.0, :r2 => 1.0, :albedo => 0.3))]
    end
    f = unique_name(:foil)
    m = unique_name(:main)
    return [PartSpec(:rectangle; name = f, frame, foil = true,
                     values = Dict(:x => 6.0, :theta => 90.0, :phi => 180.0, :rotation => 180.0,
                                   :r1 => 2.0, :r2 => 2.0,
                                   :albedo => albedo(foil_material()), :conversion => foil_material().conversion,
                                   :h => Float64(capsule_h))),
            PartSpec(:sphere; name = m, frame, values = Dict(:r1 => 1.0))]
end

"""
    insert_starter!(ts::TargetSetup, kind; frame = :tcc) -> Bool

Add the parts of starter `kind` ([`starter_parts`](@ref)) to the design, placed in the
frame named `frame` (a part's, a piece's or TCC), their names made unique against the
design's parts and components (`hohlraum_2`, …). The first foil a starter brings
becomes the active foil when the design had none; a starter's observed surface takes
the role from the part that held it. A design the geometry rejects sets `error` and
returns `false`.
"""
function insert_starter!(ts::TargetSetup, kind::Symbol; frame::Symbol = :tcc)
    try
        haskey(part_frames(ts), frame) || throw(ArgumentError("no frame named $frame"))
        taken = union(Set(all_component_names(ts)), (p.name for p in ts.parts[]), keys(part_frames(ts)))
        new = starter_parts(kind; frame, taken, capsule_h = ts.mesh[].capsule_h)
        ps = copy(ts.parts[])
        any(p -> p.observed, new) && _drop_observed!(ps)
        first_foil = active_foil(ts) === nothing ? findfirst(p -> p.foil, new) : nothing
        _try_build!(ts; parts = [ps; new]) || return false
        first_foil === nothing || (ts.active_foil[] = new[first_foil].name)
        return true
    catch err
        err isa ArgumentError || rethrow()
        ts.error[] = _errmsg(err)
        return false
    end
end
