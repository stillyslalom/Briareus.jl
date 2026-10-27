# Local only: regenerate the Briareus icon files from one drawing.
#
#     julia docs/make_icons.jl
#
# Writes the SVG (docs sidebar logo) and a multi-size ICO (docs favicon, Windows title
# bar and taskbar icon of the GUI window) to docs/src/assets/ and
# BriareusGUI/src/qml/icons/, and to the latter a 256 px PNG (the launcher splash and the
# Linux menu entry) and an ICNS (the macOS app bundle). The drawing: a backlighter shot, three focusing beam
# cones in Julia's green, red and purple converging on one spot of a blue foil. From
# 64 px up, x-rays fan out of the back of the foil; smaller sizes leave them out.
using Pkg
Pkg.activate(; temp = true)
Pkg.add("Librsvg_jll"; io = devnull)
using Librsvg_jll

const BLUE, RED, GREEN, PURPLE = "#4063D8", "#CB3C33", "#389826", "#9558B2"

fmt(v) = string(round(v; digits = 2))

# A cone 39 units long arriving at `deg` from the -x axis, base half-width 7, its apex
# 1 unit short of the focal spot on the foil face.
function cone(deg, color; focus = (44.0, 32.0), len = 39.0, half = 7.0, stop = 1.0)
    u = (cosd(deg), sind(deg))
    p = focus .- len .* u
    n = (-u[2], u[1])
    pts = (p .+ half .* n, focus .- stop .* u, p .- half .* n)
    return """<polygon points="$(join((fmt(x) * "," * fmt(y) for (x, y) in pts), " "))" fill="$color"/>"""
end

function drawing(; detail)
    xrays = detail ? join("""<line x1="56" y1="$y1" x2="$x2" y2="$y2" stroke="$BLUE" stroke-width="2.5" stroke-linecap="round"/>"""
                          for (y1, x2, y2) in ((32, 62, 32), (27, 61, 20), (37, 61, 44))) : ""
    return """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">\
           $xrays$(cone(38, GREEN))$(cone(-38, PURPLE))$(cone(0, RED))\
           <rect x="44" y="6" width="8" height="52" rx="2.5" fill="$BLUE"/></svg>
           """
end

# An ICO holding one PNG per size (Windows Vista and later, and browsers).
function write_ico(path, pngs::Vector{Pair{Int,Vector{UInt8}}})
    open(path, "w") do io
        write(io, UInt16(0), UInt16(1), UInt16(length(pngs)))
        offset = 6 + 16 * length(pngs)
        for (px, data) in pngs
            b = UInt8(px >= 256 ? 0 : px)
            write(io, b, b, UInt8(0), UInt8(0), UInt16(1), UInt16(32),
                  UInt32(length(data)), UInt32(offset))
            offset += length(data)
        end
        foreach(p -> write(io, p.second), pngs)
    end
end

# An ICNS holding one PNG per size (macOS 10.7 and later): big-endian type and length
# (8-byte header included) per entry, after the file header.
const ICNS_TYPES = Dict(16 => "icp4", 32 => "icp5", 64 => "icp6", 128 => "ic07", 256 => "ic08",
                        512 => "ic09", 1024 => "ic10")
function write_icns(path, pngs::Vector{Pair{Int,Vector{UInt8}}})
    entries = [(ICNS_TYPES[px], data) for (px, data) in pngs if haskey(ICNS_TYPES, px)]
    open(path, "w") do io
        write(io, "icns", hton(UInt32(8 + sum(8 + length(d) for (_, d) in entries))))
        for (t, d) in entries
            write(io, t, hton(UInt32(8 + length(d))), d)
        end
    end
end

root = dirname(@__DIR__)
dirs = (joinpath(root, "docs", "src", "assets"), joinpath(root, "BriareusGUI", "src", "qml", "icons"))
mktempdir() do tmp
    svgs = Dict(d => joinpath(tmp, "$d.svg") for d in (false, true))
    for (d, path) in svgs
        write(path, drawing(; detail = d))
    end
    pngs = map([16, 24, 32, 48, 64, 128, 256, 512, 1024]) do px
        out = joinpath(tmp, "$px.png")
        run(`$(rsvg_convert()) -w $px -h $px -o $out $(svgs[px >= 64])`)
        px => read(out)
    end
    for d in dirs
        mkpath(d)
        write(joinpath(d, d == dirs[1] ? "logo.svg" : "briareus.svg"), drawing(; detail = true))
        write_ico(joinpath(d, d == dirs[1] ? "favicon.ico" : "briareus.ico"), filter(p -> p.first <= 256, pngs))
    end
    write(joinpath(dirs[2], "briareus.png"), Dict(pngs)[256])
    write_icns(joinpath(dirs[2], "briareus.icns"), pngs)
end
println("wrote the icons to ", join(dirs, " and "))
