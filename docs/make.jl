include(joinpath(dirname(@__DIR__), "src", "MultidimensionalTools.jl"))
using Documenter, .MultidimensionalTools

Documenter.makedocs(
    clean = true,
    doctest = true,
    modules = Module[MultidimensionalTools],
    repo = "",
    highlightsig = true,
    sitename = "MultidimensionalTools Documentation",
    expandfirst = [],
    pages = [
        "Index" => "index.md",
    ]
)

deploydocs(;
    repo  =  "github.com/jakewilliami/MultidimensionalTools.jl.git",
)
