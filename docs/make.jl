const T_MAKE_START = time()

using Briareus
using Documenter
using Literate
using Printf

# Tutorials are Literate scripts in `literate/`, converted to markdown in
# `src/tutorials/` (generated, gitignored) and executed by Documenter through their
# `@example` blocks, so the documentation build runs every tutorial end to end. They
# are registered here rather than discovered, so that adding a script never adds an
# expensive page to every build. Each page reports its own execution time: the
# conversion wraps it in a pair of hidden `@setup` blocks that run in the page's module
# and record the elapsed time in `Main.PAGE_TIMES`.
const LITERATE_DIR = joinpath(@__DIR__, "literate")
const TUTORIAL_DIR = joinpath(@__DIR__, "src", "tutorials")
const TUTORIALS = [
    "laser_spot.jl",
    "integrating_sphere.jl",
    "dropped_bundle.jl",
    "diagnostic_windows.jl",
]

const PAGE_TIMES = Dict{String,Float64}()
const PHASES = Tuple{String,Float64}[]

macro phase(name, ex)
    quote
        local t0 = time()
        local val = $(esc(ex))
        push!(PHASES, ($(esc(name)), time() - t0))
        val
    end
end

"Wrap a Literate-generated page so it reports its own execution time."
function time_tutorial(name)
    start = "\n```@setup $name\nMain.PAGE_TIMES[\"$name\"] = -time()\n```\n"
    stop = "\n```@setup $name\nMain.PAGE_TIMES[\"$name\"] += time()\n```\n"
    return function (markdown)
        cut = findfirst("```\n", markdown)
        at = startswith(markdown, "```@meta") && cut !== nothing ? last(cut) : 0
        return markdown[1:at] * start * markdown[at+1:end] * stop
    end
end

push!(PHASES, ("load", time() - T_MAKE_START))

@phase "literate" begin
    mkpath(TUTORIAL_DIR)
    for file in TUTORIALS
        name = splitext(file)[1]
        Literate.markdown(joinpath(LITERATE_DIR, file), TUTORIAL_DIR;
                          documenter = true, credit = true,
                          postprocess = time_tutorial(name))
    end
end

tutorial_page(file) = joinpath("tutorials", splitext(file)[1] * ".md")

@phase "makedocs" makedocs(;
    sitename = "Briareus.jl",
    modules = [Briareus],
    authors = "Alex Ames",
    repo = Remotes.GitHub("stillyslalom", "Briareus.jl"),
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true", edit_link = "main",
                             size_threshold = 400 * 1024, example_size_threshold = 0),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Laser spot on a foil" => tutorial_page("laser_spot.jl"),
            "Integrating sphere" => tutorial_page("integrating_sphere.jl"),
            "Freeing a bundle for radiography" => tutorial_page("dropped_bundle.jl"),
            "Diagnostic windows in the hohlraum wall" => tutorial_page("diagnostic_windows.jl"),
        ],
        "Explanation" => [
            "Verification and validation" => "explanation/verification-validation.md",
        ],
        "Reference" => [
            "Geometry and targets" => "reference/geometry.md",
            "Beams and facilities" => "reference/beams.md",
            "Laser deposition" => "reference/deposition.md",
            "Radiation" => "reference/radiation.md",
            "Mode analysis" => "reference/analysis.md",
            "Optimization" => "reference/optimization.md",
            "Laser–plasma models" => "reference/lpi.md",
            "Shot review" => "reference/review.md",
            "Input and output" => "reference/io.md",
            "VisRad interchange" => "reference/visrad.md",
        ],
    ],
    checkdocs = :exports,
    warnonly = false,
)

if get(ENV, "CI", "false") == "true"
    @phase "deploy" deploydocs(; repo = "github.com/stillyslalom/Briareus.jl.git",
                               devbranch = "main", push_preview = false)
end

println("\nDocumentation build times")
for (name, seconds) in PHASES
    @printf("  %-10s %7.1f s\n", name, seconds)
end
for name in sort(collect(keys(PAGE_TIMES)); by = n -> -PAGE_TIMES[n])
    @printf("  tutorial  %-22s %7.1f s\n", name, PAGE_TIMES[name])
end
