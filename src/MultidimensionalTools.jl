module MultidimensionalTools

# type wrappers
export AbstractIndex, AbstractIndices, AbstractIndexOrIndices
# helper functions
export n_adjacencies, extrema_indices, tryindex, get_directions
# main functions
export append_n_times, append_n_times_backwards, promote_to_nD, promote_to_3D,
    expand_as_required, adjacencies, expanded_adjacencies, n_adjacent_to,
    expanded_n_adjacent_to, global_adjacencies_indices, global_adjacencies,
    global_n_adjacent_to

include("utils.jl")
include("generic.jl")
include("adjacent.jl")
# include("hyperdet.jl")

end # end module
