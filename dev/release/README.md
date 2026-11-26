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

No build has finished yet. The precompile run (every headless command, with and without
the window) completes, but compiling the system image is killed for memory on a Linux
x64 machine whose processes share a limit of about 14.5 GB. The table gives the runs of
2026-11-26 (Julia 1.13.4, PackageCompiler 2.4.5), each killed by the kernel when the
compiling process reached the limit, so the peaks are lower bounds:

| Settings | Resident memory when killed |
|---|---|
| `JULIA_IMAGE_THREADS=1`, `--cpu-target=native` | 13.5 GB |
| `JULIA_IMAGE_THREADS=1`, `--cpu-target=native --no-window` | 13.5 GB |
| `JULIA_IMAGE_THREADS=1`, `--cpu-target=native --no-window --incremental` | 13.5 GB |
| `JULIA_IMAGE_THREADS=1`, `--cpu-target=native --no-window`, `BRIAREUS_NO_WORKLOAD=true` | 13.6 GB |

Earlier runs with the default targets, with `--cpu-target=native` and with
`--cpu-target=native --incremental`, all with the default number of image threads, were
killed at about 14 GB on the same machine. None of these settings brings the compile
under the limit: one code-generation thread, leaving the window out of the precompile
run, building on the system image of Julia, and skipping the precompile workload of
Briareus. The memory is therefore taken by compiling the package stack itself
(GLMakie, QML and their dependencies with Briareus), not by the choice of settings, and
a build needs a machine with well over 14 GB, such as 32 GB. The standard GitHub runners
have 16 GB (Linux and Windows) and 7 GB (macOS arm64), so the workflow needs larger
runners. The other route is to ship Julia with a depot of package images (no custom
system image), which builds in far less memory but starts more slowly. Neither the
archives nor `smoke.sh` on a built application have been run.
