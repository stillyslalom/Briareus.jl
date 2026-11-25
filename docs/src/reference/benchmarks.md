# Published data and benchmark fixtures

The digitized figures and tables of the papers that the models are fitted to or checked
against live in the submodule `Briareus.Benchmarks`, with the fitted CBET models of those
shots and the mock-up pointing of the examples. Its names are not exported from
`Briareus`; `using Briareus.Benchmarks` brings them in. The data are in `data/*.csv`,
each with its sources in the header, and every function cites its paper.

```@autodocs
Modules = [Briareus.Benchmarks]
Pages = ["src/benchmarks.jl"]
Private = false
```
