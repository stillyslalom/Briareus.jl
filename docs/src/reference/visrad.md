# VisRad interchange

Readers and writers for the files VisRad (Prism Computational Sciences, VISRAD 21.1.1)
imports and exports — NIF SST pointing, the Laser Beam Summary table, time-dependent
beam powers, element and power-summary tables, workspaces (`.vrw`: target components
with bundled CAD, and beams) — and a benchmark kit that writes
everything needed to rebuild a Briareus case in VisRad together with Briareus's results
for comparison. Formats follow the [VisRad manual](https://prism-cs.com/Manuals/VisRad/);
the provenance and confidence of each detail are in `reference/research-notes.md`
("VisRad interchange formats").

```@autodocs
Modules = [Briareus]
Pages = ["src/visrad.jl", "src/visrad_workspace.jl"]
Private = false
```
