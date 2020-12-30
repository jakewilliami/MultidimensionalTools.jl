module MultidimensionalTools

# type wrappers
export AbstractIndex, AbstractIndices, AbstractIndexOrIndices
# helper functions
export n_adjacencies, extrema_indices, tryindex, get_directions
# main functions
export append_n_times, append_n_times_backwards, promote_to_nD, promote_to_3D,
    reshape_as_required, adjacencies, expanded_adjacencies, n_adjacent_to,
    expanded_n_adjacent_to, global_adjacencies_indices, global_adjacencies,
    global_n_adjacent_to

"""
    AbstractIndex{N} = Union{NTuple{N, T}, CartesianIndex{N}}
"""
AbstractIndex{N} = Union{NTuple{N, T}, CartesianIndex{N}} where {T}

"""
    AbstractIndices{N}  = Union{AbstractArray{I, M}, NTuple{M, I}, CartesianIndices{N, NTuple{N, Base.OneTo{T}}}} where I <: AbstractIndex{N}
"""
AbstractIndices{N}  = Union{AbstractArray{I, M}, NTuple{M, I}, CartesianIndices{N, NTuple{N, Base.OneTo{T}}}} where {I <: AbstractIndex{N}, T, M}

"""
    AbstractIndexOrIndices{N} = Union{AbstractIndex{N}, AbstractIndices{N}}
"""
AbstractIndexOrIndices{N} = Union{AbstractIndex{N}, AbstractIndices{N}}

"""
    n_adjacencies(dim::Integer) -> Integer
    n_adjacencies(M::AbstractArray{T, N}) -> Integer
    n_adjacencies(I::AbstractIndexOrIndices{T, N}) -> Integer

Given a matrix or dimension, returns the number of elements adjacent to any given element in an infinite lattice/matrix (i.e., not including edges).  For edges, see `adjacencies`).
"""
n_adjacencies(dim::Integer) = 3^dim - 1
n_adjacencies(M::AbstractArray{T, N}) where {T, N} = n_adjacencies(ndims(M))
n_adjacencies(I::AbstractIndexOrIndices{N}) where {N} = n_adjacencies(length(first(I)))
# n_adjacencies(M::AbstractArray{T, N}) where {T, N} = n_adjacencies(N)

"""
    𝟘(n::Integer) -> NTuple{n, Integer}
    \\bbzero[TAB](n::Integer) -> NTuple{n, Integer}

Simply returns the zero tuple of a given dimension.
"""
𝟘(n::Integer) = ntuple(_ -> zero(typeof(n)), n)

"""
Returns an array of tuples of directions.

```julia
    get_directions(dim::Integer; include_zero::Bool = false) -> Array{NTuple{N, Integer}, 1}

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
function get_directions(dim::T; include_zero::Bool = false) where T <: Integer
    D = reshape(collect(Iterators.product([-one(T):one(T) for _ in 1:dim]...)), :)
    return include_zero ? D : filter(d -> d ≠ 𝟘(dim), D)
end

get_directions(M::AbstractArray{T, N}; include_zero::Bool = false) where {T, N} = get_directions(ndims(M), include_zero = include_zero)
# get_directions(M::AbstractArray{T, N}; include_zero::Bool = false) where {T, N} = get_directions(N, include_zero = include_zero)

function __tryindex_internal(M::AbstractArray{T, N}, inds::AbstractIndexOrIndices{N}...) where {T, N}
    indices = Vector{Union{T, Nothing}}()
    for idx in inds
        try
            push!(indices, getindex(M, idx...))
        catch
            push!(indices, nothing)
        end
    end

    return Tuple(indices)
end

"""
    tryindex(M::AbstractArray{T, N}, inds::AbstractIndex{N}}...)
    tryindex(M::AbstractArray{T, N}, inds::AbstractIndices{N}})

For each index specified, gets the index or missing if the index is unavailable.
"""
tryindex(M::AbstractArray{T, N}, inds::AbstractIndex{N}...) where {T, N} = __tryindex_internal(M, inds...)

tryindex(M::AbstractArray{T, N}, inds::AbstractIndices{N}) where {T, N} = __tryindex_internal(M, inds...)

function __extrema_indices_internal(I::AbstractIndexOrIndices{N}...) where {N}
    return Tuple(extrema(reduce(vcat, permutedims(collect(i)) for i in I...), dims = 1))
end
"""
    extrema_indices(I::AbstractIndices{N}) -> NTuple{N, NTuple{N, Integer}}
    extrema_indices(i::AbstractIndex{N}...) -> NTuple{N, NTuple{N, Integer}}

Given indices, returns an `NTuple`.  Each element in the `NTuple` represents a (min, max) tuple for each dimension.

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

julia> # that is, the (min, max) for dim 2 is (1, 5); etc.
```
"""
extrema_indices(I::AbstractIndices{N}) where {N} = __extrema_indices_internal(I)
extrema_indices(i::AbstractIndex{N}...) where {N} = __extrema_indices_internal(i)

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
function append_n_times(M::AbstractArray{T, N}, n::Integer, fill_elem::T; dims::Integer = 1) where {T, N}
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
function append_n_times_backwards(M::AbstractArray{T, N}, n::Integer, fill_elem::T; dims::Integer = 1) where {T, N}
    sz = ntuple(d -> d == dims ? n : size(M, d), max(N, dims))
    return cat(fill(fill_elem, sz), M; dims = dims)
end

"""
    function promote_to_nD(M::AbstractArray{T, N}, n::Integer, fill_elem::T)

Assumes the given matrix M is an (n - 1) dimensional slice of an n-dimensional structure, and fills in the array to n dimensions with `fill_elem`.
"""
function promote_to_nD(M::AbstractArray{T, N}, n::Integer, fill_elem::T) where {T, N}
    ndims(M) == n && return M
    n < ndims(M) && throw(error("Cannot reduce the number of dimensions this array has.  See `resize`."))
    
    for d in (ndims(M) + 1):n
        M = append_n_times(M, 1, fill_elem, dims = d)
        M = append_n_times_backwards(M, 1, fill_elem, dims = d)
    end
    
    return M
end

"""
    promote_to_3D(M::AbstractArray{T, N}, fill_elem::T)
    
A simple, self-evident wrapper for `promote_to_nD`.
"""
promote_to_3D(M::AbstractArray{T, N}, fill_elem::T) where {T, N} = promote_to_nD(M, 3, fill_elem)

function __reshape_as_required_internal(M::AbstractArray{T, N}, expand_by::T, inds::AbstractIndexOrIndices{N}...) where {T, N}
    indices = Vector{Union{T, Nothing}}()
    ind_extrema = extrema_indices(inds...)
    
    for d in 1:ndims(M)
        ith_extrema = ind_extrema[d]
        if ! all(map(m -> m ≤ size(M, d), ith_extrema))
            for invalid_idx in filter(i -> i > size(M, i) || i < 1, ith_extrema)
                difference = invalid_idx - size(M, d)
                M = sign(difference) == 1 ? append_n_times(M, abs(difference), expand_by, dims = d) : append_n_times_backwards(M, abs(difference), expand_by, dims = d)
            end
        end
    end
    
    return M
end

"""
    reshape_as_required(M::AbstractArray{T, N}, expand_by::T, inds::AbstractIndexOrIndices{N})

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
reshape_as_required(M::AbstractArray{T, N}, expand_by::T, inds::AbstractIndices{N}) where {T, N} = __reshape_as_required_internal(M, expand_by, inds)
reshape_as_required(M::AbstractArray{T, N}, expand_by::T, inds::AbstractIndex{N}...) where {T, N} = __reshape_as_required_internal(M, expand_by, inds)

"""
    adjacencies(M::AbstractArray{T, N}, idx::AbstractIndex{N})

`adjacencies` will get all adjacent element to an index in a given matrix.
"""
function adjacencies(M::AbstractArray{T, N}, idx::AbstractIndex{N}) where {T, N}
    direction_multipliers = get_directions(M, include_zero = false)
    adjacent_indices = AbstractIndex{N}[idx .+ j for j in direction_multipliers]
    return T[k for k in tryindex(M, adjacent_indices) if ! isnothing(k)]
end

"""
    expanded_adjacencies(M::Array{T, N}, expand_by::T, idx::AbstractIndex{N})

`expanded_adjacencies` will get all adjacencies of an index in the matrix M, given that the matrix is infinitely expanding by a single element beyond the specified M.  See also `adjacencies` and `reshape_as_required`.
"""
function expanded_adjacencies(M::Array{T, N}, expand_by::T, idx::AbstractIndex{N}) where {T, N}
    𝟎 = ntuple(_ -> zero(Int), N)
    idx = Tuple(idx)
    direction_multipliers = get_directions(M, include_zero = false)
    adjacent_indices = AbstractIndex{N}[idx .+ j for j in direction_multipliers]
    M = reshape_as_required(M, expand_by, adjacent_indices)
    idx_shift = abs.(size(M) .- idx) # in case we have had to shift the array for non-positive indices
    return T[M[(i .+ idx_shift .- 1)...] for i in adjacent_indices]
end

"""
    n_adjacent_to(M::AbstractArray{T, N}, idx::AbstractIndex{N}, adj_elem::T)

Given a matrix M, counts the number of adjacent elements to index `idx` that are exactly the `adj_elem`.
"""
function n_adjacent_to(M::AbstractArray{T, N}, idx::AbstractIndex{N}, adj_elem::T) where {N, T}
    n = 0
    
    for i in adjacencies(M, idx)
        if i == adj_elem
            n += 1
        end
    end
    
    return n
end

"""
    expanded_n_adjacent_to(M::AbstractArray{T, N}, expand_by::T, idx::AbstractIndex{N}, adj_elem::T)

Given a matrix M, counts the number of adjacent elements to index `idx` that are exactly the `adj_elem`, expanding the matrix if needed (see `expanded_adjacencies`).
"""
function expanded_n_adjacent_to(M::AbstractArray{T, N}, expand_by::T, idx::AbstractIndex{N}, adj_elem::T) where {N, T}
    n = 0
    
    for i in expanded_adjacencies(M, expand_by, idx)
        if i == adj_elem
            n += 1
        end
    end
    
    return n
end

"""
    global_adjacencies(M::AbstractArray{T, N}, idx::AbstractIndex{N}, ignored_elem::T)

Returns the "globally" adjacent indices in arbitrary positions in the cardinal directions from a specified index in matrix `M`, ignoring certain adjacent elements (i.e., skipping over them).
"""
function global_adjacencies_indices(M::AbstractArray{T, N}, idx::AbstractIndex{N}, ignored_elem::T) where {T, N}
    adjacent_indices = AbstractIndex{N}[]
    directional_shifts = get_directions(M, include_zero = false)
    n_adjacent, adjacent_count = n_adjacencies(M), 0

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
    global_adjacencies(M::AbstractArray{T}, idx::AbstractIndex{N}, ignored_elem::T)

Using `global_adjacencies_indices`, returns the elements of each globally adjacent index.  This function is mainly used for `global_n_adjacent_to`.
"""
function global_adjacencies(M::AbstractArray{T, N}, idx::AbstractIndex{N}, ignored_elem::T) where {N, T}
    return T[M[i...] for i in global_adjacencies_indices(M, idx, ignored_elem)]
end

"""
    global_n_adjacent_to(M::AbstractArray{T, N},idx::AbstractIndex{N}, adj_elem::T, ignored_elem::T) where {N, T}

Given a matrix M, counts the number of adjacent elements to index `idx` that are exactly the `adj_elem`, skipping over certain elements if needed (see `expanded_adjacencies`).
"""
function global_n_adjacent_to(M::AbstractArray{T, N},idx::AbstractIndex{N}, adj_elem::T, ignored_elem::T) where {N, T}
    n = 0
    
    for i in global_adjacencies(M, idx, ignored_elem)
        if i == adj_elem
            n += 1
        end
    end
    
    return n
end

end # end module
