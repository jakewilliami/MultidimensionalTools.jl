include(joinpath(dirname(@__DIR__), "src", "MultidimensionalTools.jl"))
using Documenter, .MyCoolPackage

Documenter.makedocs(
    clean = true,
    doctest = true,
    modules = Module[MyCoolPackage],
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
