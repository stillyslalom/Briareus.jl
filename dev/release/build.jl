# Builds the relocatable `briareus` application with PackageCompiler (USERBASE_ROADMAP.md,
# E.1): one system image holding Briareus, the controllers, the window and the headless
# commands, and a launcher `bin/briareus` (`bin/briareus.exe` on Windows) that needs no
# Julia installation on the machine that runs it.
#
#     julia --project=dev/release -e 'using Pkg; Pkg.instantiate()'
#     julia --project=dev/release dev/release/build.jl [out] [--no-window] [--cpu-target=<targets>] [--incremental]
#
# `out` defaults to build/Briareus. `--cpu-target` replaces PackageCompiler's default
# targets (a generic build with code for several x86 generations, which runs on any
# machine of the architecture). The default needs about 14 GB of memory to compile the
# system image on Linux x64; `--cpu-target=native` needs less and runs only on machines
# like the one that built it. `--incremental` builds on the system image of the Julia
# that runs the build, which compiles less at once. The precompile run (precompile.jl) runs every headless
# command on a starter design, and opens and closes the window unless --no-window is given
# or no display is available (it needs an OpenGL context). Signing and installers are
# separate steps (README.md in this folder).

using PackageCompiler

const ROOT = dirname(dirname(@__DIR__))
const GUI = joinpath(ROOT, "BriareusGUI")
const POSITIONAL = filter(a -> !startswith(a, "--"), ARGS)
const OUT = abspath(isempty(POSITIONAL) ? joinpath(ROOT, "build", "Briareus") : POSITIONAL[1])

# The GUI environment finds the core package in this checkout, as in CI. Developing it
# writes a [sources] entry into BriareusGUI/Project.toml, which must not be committed
# (CI checks that it has none), so the file is put back after the build.
project = joinpath(GUI, "Project.toml")
saved = read(project, String)
try
    run(`$(Base.julia_cmd()) --project=$GUI -e "using Pkg; Pkg.develop(PackageSpec(path = raw\"$ROOT\")); Pkg.instantiate()"`)
    ENV["BRIAREUS_RELEASE_WINDOW"] = "--no-window" in ARGS ? "false" : get(ENV, "BRIAREUS_RELEASE_WINDOW", "true")
    cpu = [a[length("--cpu-target=")+1:end] for a in ARGS if startswith(a, "--cpu-target=")]
    kw = isempty(cpu) ? (;) : (; cpu_target = only(cpu))
    "--incremental" in ARGS && (kw = (; kw..., incremental = true))
    t = @elapsed create_app(GUI, OUT; kw...,
                            executables = ["briareus" => "julia_main"],
                            precompile_execution_file = joinpath(@__DIR__, "precompile.jl"),
                            include_lazy_artifacts = true,
                            force = true)
    println("built ", OUT, " in ", round(t / 60; digits = 1), " min")
finally
    write(project, saved)
end
