## How to run the benchmarks

The files in this folder define a benchmark suite with the tools provided by [PkgBenchmark](https://github.com/JuliaCI/PkgBenchmark.jl) and [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl).

To run the benchmarks, execute:

```julia
julia> using PkgBenchmark

julia> results = benchmarkpkg("Femshop")
```

## How to compare benchmarks

To compare current version to another tagged version, commit or branch:

```julia
julia> results = judge("TaylorModels", <tagged-version-or-branch>)
```

## Exporting results

To export the benchmark results to a Markdown file:

```julia
julia> export_markdown("results.md", results)
```

To export the benchmark results to a JSON file:

```julia
julia> writeresults("results.json", results)
```