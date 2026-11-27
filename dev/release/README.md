# Release builds

The files in this folder build `briareus` as a relocatable application with
PackageCompiler (item E.1 of `USERBASE_ROADMAP.md`). The application is a folder holding
one system image with Briareus, the controllers, the window and the headless commands,
the Julia runtime and every artifact, and a launcher `bin/briareus` (`bin/briareus.exe`
on Windows) that needs no Julia installation.

```bash
julia --project=dev/release -e 'using Pkg; Pkg.instantiate()'
julia --project=dev/release dev/release/build.jl build/Briareus      # add --no-window without a display
dev/release/smoke.sh build/Briareus                                   # the launcher, with no Julia on PATH
```

`build.jl` develops the core package of this checkout into the GUI environment and puts
`BriareusGUI/Project.toml` back afterward, since that step writes a `[sources]` entry
that CI refuses. `precompile.jl` runs every headless command on a starter design and
opens and closes the window, so a first call compiles nothing. `smoke.sh` writes a
design with `briareus run fixture.jl`, runs `talis`, `review-1w` and `roe` on it, checks
their files and prints the time of each command.

`.github/workflows/Release.yml` runs these on Windows x64, macOS arm64 and Linux x64 for
a tag `BriareusGUI-v*` or by hand, and keeps each archive as a workflow artifact for 30
days. It publishes nothing.

## Steps left to the maintainers

- Signing secrets. Windows Authenticode signing runs when the repository secrets
  `WINDOWS_CERTIFICATE` (a .pfx, base64) and `WINDOWS_CERTIFICATE_PASSWORD` are set.
  macOS signing and notarization (`sign_macos.sh`) run when `MACOS_CERTIFICATE` (a .p12,
  base64), `MACOS_CERTIFICATE_PASSWORD`, `MACOS_IDENTITY`, `APPLE_ID`, `APPLE_TEAM_ID`
  and `APPLE_APP_PASSWORD` are set. Without them the archives are unsigned, and
  SmartScreen and Gatekeeper warn on first launch; an IT allowlist of the folder is the
  alternative.
- Installers. No InnoSetup or MSI script exists for Windows, and the macOS disk image
  holds the folder, not a `Briareus.app` bundle with an `Info.plist` and icon.
- The GitHub release itself is drafted by hand from the workflow artifacts.

## What was measured

All runs are on Linux x64 (Julia 1.13.4, PackageCompiler 2.4.5) in a container whose
processes share a memory limit of about 13.6 GB, with `JULIA_IMAGE_THREADS=1`,
`--cpu-target=native` and `--no-window` unless the table says otherwise. The memory was
sampled every 2 s and the compiling process was stopped by a native stack sample with
gdb, so each phase below is named from the stack, not inferred.

| Flags of the compiling process | Peak resident memory | Result |
|---|---|---|
| default (`-g1`) | over 13.2 GB | killed while writing the data of the image |
| `-g0 --strip-ir --strip-metadata` | 12.1 GB | built; the application stops at startup with `Core.MissingCodeError` |
| `-g0` (the default of build.jl) | 12.6 GB | built in 29 min; `smoke.sh` passes |

The compiling process goes through four phases. Loading the packages takes it to about
2 GB. Selecting and inferring the methods to compile then runs for six to seven minutes
at that level: for a nonincremental system image Julia 1.13 compiles every method of
every loaded module whose signature is concrete, and every specialization already
inferred in the precompile caches of the packages, Makie included. The precompile run of
precompile.jl adds 1,456 statements, about 500 of them in Briareus, a small part of the
whole. LLVM code for all these methods is then generated and merged into one module, from
2 GB to 10.5 GB in about two minutes; debug information is part of it. Last, the data of
the image is written while that module is still held, and the module is optimized and
emitted. With debug information the data phase reached the limit; without it the
optimization phase holds about 11.7 GB for about 15 minutes, with a short peak of 12.6 GB
as the object file is written. `--strip-ir` removes the inferred code that the compiler
of the running application still needs, so it cannot be used.

Earlier runs with debug information (default targets, `--incremental`, the window in the
precompile run, and `BRIAREUS_NO_WORKLOAD=true`) were all stopped at the same limit, so
they show only that each needs more than about 13.5 GB.

The application built with `-g0` is 2.2 GB, of which the system image is 1.2 GB.
`briareus --help` takes about 1.0 s; the commands of `smoke.sh` take 2 to 6 s each,
startup included. The window, from a build whose precompile run left it out, draws its
first frame about 15 s after launch under Xvfb with software OpenGL. Headless commands on
a Linux machine without a display stop when QML starts (the module initializer of QML
creates a Qt application as the system image loads); `QT_QPA_PLATFORM=offscreen` avoids
it, and `smoke.sh` needs it on a runner without a display.

The 12.6 GB peak fits the 16 GB standard Linux and Windows runners with little headroom,
for native code without the window; the release builds compile for the default CPU
targets and trace the window, which adds code, and their peak has not been measured.
The 7 GB macOS arm64 runner needs a larger runner. The other route is to ship Julia with
a depot of package images (no custom system image), which builds in far less memory but
starts more slowly.
