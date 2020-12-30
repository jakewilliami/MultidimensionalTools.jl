#!/usr/bin/env bash
    #=
    exec julia --project="$(realpath $(dirname $(dirname $0)))" --color=yes --startup-file=no -e "include(popfirst!(ARGS))" \
    "${BASH_SOURCE[0]}" "$@"
    =#

include(joinpath(dirname(dirname(@__FILE__)), "src", "MultidimensionalTools.jl"))

using .ComputabilityTheory
using Test

@time @testset "MultidimensionalTools.jl" begin
    
end # end testset
