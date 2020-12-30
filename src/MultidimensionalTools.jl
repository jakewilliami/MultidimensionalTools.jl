module MultidimensionalTools

# helper functions
export n_adjacencies, extrema_indices, tryindex, get_directions
# main functions
export append_n_times, append_n_times_backwards, promote_to_nD, promote_to_3D,
    reshape_as_required, adjacencies, expanded_adjacencies, n_adjacent_to,
    expanded_n_adjacent_to, global_adjacencies_indices, global_adjacencies,
    global_n_adjacent_to

"""
    n_adjacencies(dim::Int) -> Int
    n_adjacencies(M::AbstractArray) -> Int

Given a matrix or dimension, returns the number of elements adjacent to any given element in an infinite lattice/matrix (i.e., not including edges).  For edges, see ``).
"""
n_adjacencies(dim::Int) = 3^dim - 1
n_adjacencies(M::AbstractArray) = n_adjacencies(ndims(M))

"""
Returns an array of tuples of directions.

```julia
    get_directions(dim::Int; include_zero::Bool = false)

julia> get_directions(3)
27-element Array{Tuple{Int64,Int64,Int64},1}:
 (-1, -1, -1)
 (0, -1, -1)
 (1, -1, -1)
 (-1, 0, -1)
 (0, 0, -1)
 (1, 0, -1)
 (-1, 1, -1)
 (0, 1, -1)
 (1, 1, -1)
 (-1, -1, 0)
 (0, -1, 0)
 (1, -1, 0)
 (-1, 0, 0)
 (0, 0, 0)
 (1, 0, 0)
 (-1, 1, 0)
 (0, 1, 0)
 (1, 1, 0)
 (-1, -1, 1)
 (0, -1, 1)
 (1, -1, 1)
 (-1, 0, 1)
 (0, 0, 1)
 (1, 0, 1)
 (-1, 1, 1)
 (0, 1, 1)
 (1, 1, 1)
```
"""
function get_directions(dim::Int; include_zero::Bool = false)
    𝟎 = ntuple(_ -> zero(Int), dim)
    D = reshape(collect(Iterators.product([-1:1 for _ in 1:dim]...)), :)
    return include_zero ? D : filter(d -> d ≠ 𝟎, D)
end

get_directions(M::AbstractArray; include_zero::Bool = false) = get_directions(ndims(M), include_zero = include_zero)

"""
    tryindex(M::Matrix{T}, inds::NTuple{N, Int}...)

For each index specified, gets the index or nothing if the index is unavailable.
"""
function tryindex(M::Array{T, N}, inds::NTuple{N, Int}...) where {T, N}
    indices = Vector{Union{T, Nothing}}()
    
    for idx in inds
        try
            push!(indices, getindex(M, idx...))
        catch
            push!(indices, nothing)
        end
    end
    
    return indices
end

"""
Given indices, returns an NTuple.  Each element in the NTuple represents a (min, max) tuple for each dimension

```julia
julia> T = [(3, 1), (4, 5), (6, 2), (1, 1)]
4-element Array{Tuple{Int64,Int64},1}:
 (3, 1)
 (4, 5)
 (6, 2)
 (1, 1)
 
julia> extrema_indices(T)
1×2 Array{Tuple{Int64,Int64},2}:
(1, 6)  (1, 5)
```
"""
function extrema_indices(I::Union{Vector{NTuple{N, Int}}, Vector{CartesianIndex{N}}}) where N
    return vec(CartesianIndex.(extrema(reduce(vcat, permutedims(collect(i)) for i in I), dims = 1)))
end
"""
Repeats a specified value n many times along the specified dimension.

```julia
julia> A = rand(Int8, 2, 2)
2×2 Array{Int8,2}:
  28  23
 -47  54

julia> append_n_times(A, 2, Int8(3), dims = 1) # repeat the value 3 twice along the first dimension
4×2 Array{Int8,2}:
  28  23
 -47  54
   3   3
   3   3
```
"""
function append_n_times(M::Array{T, N}, n::Int, fill_elem::T; dims::Int = 1) where {T, N}
    sz = ntuple(d -> d == dims ? n : size(M, d), max(N, dims))
    return cat(M, fill(fill_elem, sz); dims = dims)
end

"""
```julia
julia> A = rand(Int8, 2, 2)
2×2 Array{Int8,2}:
  28  23
 -47  54

julia> append_n_times_backwards(A, 2, Int8(3), dims = 1) # repeat the value 3 twice along the first dimension
4×2 Array{Int8,2}:
  3   3
  3   3
  28  23
 -47  54
```
"""
function append_n_times_backwards(M::Array{T, N}, n::Int, fill_elem::T; dims::Int = 1) where {T, N}
    sz = ntuple(d -> d == dims ? n : size(M, d), max(N, dims))
    return cat(fill(fill_elem, sz), M; dims = dims)
end

"""
    function promote_to_nD(M::Array{T, N}, n::Int, fill_elem::T)

Assumes the given matrix M is an (n - 1) dimensional slice of an n-dimensional structure, and fills in the array to n dimensions with `fill_elem`.
"""
function promote_to_nD(M::Array{T, N}, n::Int, fill_elem::T) where {T, N}
    ndims(M) == n && return M
    n < ndims(M) && throw(error("Cannot reduce the number of dimensions this array has.  See `resize`."))
    
    for d in (ndims(M) + 1):n
        M = append_n_times(M, 1, fill_elem, dims = d)
        M = append_n_times_backwards(M, 1, fill_elem, dims = d)
    end
    
    return M
end

"""
    promote_to_3D(M::Array{T, N}, fill_elem::T)
    
A simple wrapper for `promote_to_nD`.
"""
promote_to_3D(M::Array{T, N}, fill_elem::T) where {T, N} = promote_to_nD(M, 3, fill_elem)

"""
    reshape_as_required(M::Array{T, N}, expand_by::T, inds::Union{Vector{NTuple{N, Int}}, Vector{CartesianIndex{N}}})
    reshape_as_required(M::Array{T, N}, expand_by::T, inds::Union{NTuple{N, Int}, CartesianIndex{N}}...)

Given indices, `reshape_as_required` will fill in a matrix with `expand_by` where needed if such indices are not currently accessible.

```julia
julia> A = rand(Int8, 2, 2)
2×2 Array{Int8,2}:
 -96   83
 -94  -39

julia> reshape_as_required(A, zero(Int8), [(1, 1), (2, 2), (3, 3)])
3×3 Array{Int8,2}:
 -96   83  0
 -94  -39  0
   0    0  0
```
"""
function reshape_as_required(M::Array{T, N}, expand_by::T, inds::Union{Vector{NTuple{N, Int}}, Vector{CartesianIndex{N}}}) where {T, N}
    indices = Vector{Union{T, Nothing}}()
    ind_extrema = extrema_indices(inds)
    
    for d in 1:ndims(M)
        ith_extrema = Tuple(ind_extrema[d])
        if ! all(map(m -> m ≤ size(M, d), ith_extrema))
            for invalid_idx in filter(i -> i > size(M, i) || i < 1, ith_extrema)
                difference = invalid_idx - size(M, d)
                M = sign(difference) == 1 ? append_n_times(M, abs(difference), expand_by, dims = d) : append_n_times_backwards(M, abs(difference), expand_by, dims = d)
            end
        end
    end
    
    return M
end

function reshape_as_required(M::Array{T, N}, expand_by::T, inds::Union{NTuple{N, Int}, CartesianIndex{N}}...) where {T, N}
    return reshape_as_required(M, expand_by, Int[inds...])
end

function adjacencies(M::Array{T, N}, idx::NTuple{N, Int}) where {T, N}
    𝟎 = ntuple(_ -> zero(Int), N)
    return T[k for k in tryindex(M, NTuple{N, Int}[idx .+ j for j in NTuple{N, Int}[t for t in Base.Iterators.product([-one(Int):one(Int) for i in one(Int):N]...)] if j ≠ 𝟎]...) if ! isnothing(k)]
end

"""
    expanded_adjacencies(M::Array{T, N}, expand_by::T, idx::Union{NTuple{N, Int}, CartesianIndex{N}})

`expanded_adjacencies` will get all adjacencies of an index in the matrix M, given that the matrix is infinitely expanding by a single element beyond the specified M.  See also `reshape_as_required`.
"""
function expanded_adjacencies(M::Array{T, N}, expand_by::T, idx::Union{NTuple{N, Int}, CartesianIndex{N}}) where {T, N}
    𝟎 = ntuple(_ -> zero(Int), N)
    idx = Tuple(idx)
    adjacent_multipliers = Iterators.product([-1:1 for _ in 1:N]...)
    adjacent_indices = NTuple{N, Int}[idx .+ j for j in NTuple{N, Int}[t for t in adjacent_multipliers] if j ≠ 𝟎]
    M = reshape_as_required(M, expand_by, adjacent_indices)
    idx_shift = abs.(size(M) .- idx) # in case we have had to shift the array for non-positive indices
    return T[M[(i .+ idx_shift .- 1)...] for i in adjacent_indices]
end

"""
    n_adjacent_to(M::Array{T, N}, idx::Union{NTuple{N, Int}, CartesianIndex{N}}, adj_elem::T)

Given a matrix M, counts the number of adjacent elements to index `idx` that are exactly the `adj_elem`.
"""
function n_adjacent_to(M::Array{T, N}, idx::Union{NTuple{N, Int}, CartesianIndex{N}}, adj_elem::T) where {N, T}
    n = 0
    
    for i in adjacencies(M, idx)
        if i == adj_elem
            n += 1
        end
    end
    
    return n
end

"""
    expanded_n_adjacent_to(M::Array{T, N}, expand_by::T, idx::Union{NTuple{N, Int}, CartesianIndex{N}}, adj_elem::T)

Given a matrix M, counts the number of adjacent elements to index `idx` that are exactly the `adj_elem`, expanding the matrix if needed (see `expanded_adjacencies`).
"""
function expanded_n_adjacent_to(M::Array{T, N}, expand_by::T, idx::Union{NTuple{N, Int}, CartesianIndex{N}}, adj_elem::T) where {N, T}
    n = 0
    
    for i in expanded_adjacencies(M, expand_by, idx)
        if i == adj_elem
            n += 1
        end
    end
    
    return n
end

"""
    global_adjacencies(M::Matrix{T}, idx::Union{NTuple{N, Int}, CartesianIndex{N}}, ignored_elem::T)

Returns the "globally" adjacent indices in arbitrary positions in the cardinal directions from a specified index in matrix `M`, ignoring certain adjacent elements (i.e., skipping over them).
"""
function global_adjacencies_indices(M::Array{T, N}, idx::Union{NTuple{N, Int}, CartesianIndex{N}}, ignored_elem::T) where {T, N}
    adjacent_indices, 𝟎 = Union{NTuple{N, Int}, CartesianIndex{N}}[], ntuple(_ -> zero(Int), N)
    directional_shifts = NTuple{N, Int}[i for i in NTuple{N, Int}[t for t in Base.Iterators.product([-one(Int):one(Int) for i in one(Int):N]...)] if i ≠ 𝟎]
    n_adjacent, adjacent_count = n_adjacencies(ndims(M)), 0

    while adjacent_count < n_adjacent
        for directional_shift in directional_shifts
            while true
                adj_index = idx .+ directional_shift

                if nothing ∈ tryindex(M, adj_index)
                    n_adjacent -= 1
                    break
                end

                if M[adj_index...] ≠ ignored_elem
                    adjacent_count += 1
                    push!(adjacent_indices, adj_index)
                    break
                else
                    directional_shift = (abs.(directional_shift) .+ 1) .* sign.(directional_shift)
                end
            end
        end
    end

    return adjacent_indices
end

"""
    global_adjacencies(M::Matrix{T}, idx::Union{NTuple{N, Int}, CartesianIndex{N}}, ignored_elem::T)

Using `global_adjacencies_indices`, returns the elements of each globally adjacent index.  This function is mainly used for `global_n_adjacent_to`.
"""
function global_adjacencies(M::Array{T, N}, idx::Union{NTuple{N, Int}, CartesianIndex{N}}, ignored_elem::T) where {N, T}
    return T[M[i...] for i in global_adjacencies_indices(M, idx, ignored_elem)]
end

"""
    global_n_adjacent_to(M::Array{T, N},idx::Union{NTuple{N, Int}, CartesianIndex{N}}, adj_elem::T) where {N, T}

Given a matrix M, counts the number of adjacent elements to index `idx` that are exactly the `adj_elem`, skipping over certain elements if needed (see `expanded_adjacencies`).
"""
function global_n_adjacent_to(M::Array{T, N},idx::Union{NTuple{N, Int}, CartesianIndex{N}}, ignored_elem::T) where {N, T}
    n = 0
    
    for i in global_adjacencies(M, idx, ignored_elem)
        if i == adj_elem
            n += 1
        end
    end
    
    return n
end

end # end module
