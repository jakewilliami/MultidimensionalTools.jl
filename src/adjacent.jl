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
function get_directions(dim::T; include_zero::Bool = false) where {T <: Integer}
    D = reshape(collect(Iterators.product([(-one(T)):one(T) for _ in 1:dim]...)), :)
    return include_zero ? D : filter(d -> d ≠ 𝟘(dim), D)
end
get_directions(M::AbstractArray{T, N}; include_zero::Bool = false) where {T, N} =
    get_directions(ndims(M), include_zero = include_zero)
# get_directions(M::AbstractArray{T, N}; include_zero::Bool = false) where {T, N} = get_directions(N, include_zero = include_zero)

function _expand_as_required(
    M::AbstractArray{T, N}, expand_by::T, inds::AbstractIndexOrIndices{N}...
) where {T, N}
    indices = Vector{Union{T, Nothing}}()
    ind_extrema = extrema_indices(inds...)

    for d in 1:ndims(M)
        ith_extrema = ind_extrema[d]
        if !all(map(m -> m ≤ size(M, d), ith_extrema))
            for invalid_idx in filter(i -> i > size(M, i) || i < 1, ith_extrema)
                difference = invalid_idx - size(M, d)
                M =
                    sign(difference) == 1 ?
                    append_n_times(M, abs(difference), expand_by, dims = d) :
                    append_n_times_backwards(M, abs(difference), expand_by, dims = d)
            end
        end
    end

    return M
end

"""
    expand_as_required(M::AbstractArray{T, N}, expand_by::T, inds::AbstractIndexOrIndices{N})

Given indices, `expand_as_required` will fill in a matrix with `expand_by` where needed if such indices are not currently accessible.  If no `expand_by` is given, and the element type of the matrix is a number, fills with zero.

```julia
julia> A = rand(Int8, 2, 2)
2×2 Array{Int8,2}:
 -96   83
 -94  -39

julia> expand_as_required(A, zero(Int8), [(1, 1), (2, 2), (3, 3)])
3×3 Array{Int8,2}:
 -96   83  0
 -94  -39  0
   0    0  0
```
"""
expand_as_required(
    M::AbstractArray{T, N}, expand_by::T, inds::AbstractIndices{N}
) where {T, N} = _expand_as_required(M, expand_by, inds)
expand_as_required(
    M::AbstractArray{T, N}, expand_by::T, inds::AbstractIndex{N}...
) where {T, N} = _expand_as_required(M, expand_by, inds)
expand_as_required(
    M::AbstractArray{T, N}, inds::AbstractIndices{N}
) where {T <: Number, N} = _expand_as_required(M, zero(T), inds)
expand_as_required(
    M::AbstractArray{T, N}, inds::AbstractIndex{N}...
) where {T <: Number, N} = _expand_as_required(M, zero(T), inds)

"""
    adjacencies(M::AbstractArray{T, N}, idx::AbstractIndex{N})

`adjacencies` will get all adjacent element to an index in a given matrix.
"""
function adjacencies(M::AbstractArray{T, N}, idx::AbstractIndex{N}) where {T, N}
    direction_multipliers = get_directions(M, include_zero = false)
    adjacent_indices = AbstractIndex{N}[idx .+ j for j in direction_multipliers]
    return T[k for k in tryindex(M, adjacent_indices) if !isnothing(k)]
end

"""
    expanded_adjacencies(M::Array{T, N}, expand_by::T, idx::AbstractIndex{N})

`expanded_adjacencies` will get all adjacencies of an index in the matrix M, given that the matrix is infinitely expanding by a single element beyond the specified M.  See also `adjacencies` and `expand_as_required`.  If no `expand_by` is given, and the element type of the matrix is a number, fills with zero.
"""
function expanded_adjacencies(
    M::Array{T, N}, expand_by::T, idx::AbstractIndex{N}
) where {T, N}
    𝟎 = ntuple(_ -> zero(Int), N)
    idx = Tuple(idx)
    direction_multipliers = get_directions(M, include_zero = false)
    adjacent_indices = AbstractIndex{N}[idx .+ j for j in direction_multipliers]
    M = expand_as_required(M, expand_by, adjacent_indices)
    idx_shift = abs.(size(M) .- idx) # in case we have had to shift the array for non-positive indices
    return T[M[(i .+ idx_shift .- 1)...] for i in adjacent_indices]
end
expanded_adjacencies(M::Array{T, N}, idx::AbstractIndex{N}) where {T <: Number, N} =
    expanded_adjacencies(M, zero(T), idx)

"""
    n_adjacent_to(M::AbstractArray{T, N}, idx::AbstractIndex{N}, adj_elem::T)

Given a matrix M, counts the number of adjacent elements to index `idx` that are exactly the `adj_elem`.
"""
function n_adjacent_to(
    M::AbstractArray{T, N}, idx::AbstractIndex{N}, adj_elem::T
) where {N, T}
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

Given a matrix M, counts the number of adjacent elements to index `idx` that are exactly the `adj_elem`, expanding the matrix if needed (see `expanded_adjacencies`).  If no `expand_by` is given, and the element type of the matrix is a number, fills with zero.
"""
function expanded_n_adjacent_to(
    M::AbstractArray{T, N}, expand_by::T, idx::AbstractIndex{N}, adj_elem::T
) where {N, T}
    n = 0

    for i in expanded_adjacencies(M, expand_by, idx)
        if i == adj_elem
            n += 1
        end
    end

    return n
end
expanded_n_adjacent_to(
    M::AbstractArray{T, N}, idx::AbstractIndex{N}, adj_elem::T
) where {N, T <: Number} = expanded_n_adjacent_to(M, zero(T), idx, adj_elem)

"""
    global_adjacencies(M::AbstractArray{T, N}, idx::AbstractIndex{N}, ignored_elem::T)

Returns the "globally" adjacent indices in arbitrary positions in the cardinal directions from a specified index in matrix `M`, ignoring certain adjacent elements (i.e., skipping over them).
"""
function global_adjacencies_indices(
    M::AbstractArray{T, N}, idx::AbstractIndex{N}, ignored_elem::T
) where {T, N}
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
                    directional_shift =
                        (abs.(directional_shift) .+ 1) .* sign.(directional_shift)
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
function global_adjacencies(
    M::AbstractArray{T, N}, idx::AbstractIndex{N}, ignored_elem::T
) where {N, T}
    return T[M[i...] for i in global_adjacencies_indices(M, idx, ignored_elem)]
end

"""
    global_n_adjacent_to(M::AbstractArray{T, N},idx::AbstractIndex{N}, adj_elem::T, ignored_elem::T) where {N, T}

Given a matrix M, counts the number of adjacent elements to index `idx` that are exactly the `adj_elem`, skipping over certain elements if needed (see `expanded_adjacencies`).
"""
function global_n_adjacent_to(
    M::AbstractArray{T, N}, idx::AbstractIndex{N}, adj_elem::T, ignored_elem::T
) where {N, T}
    n = 0

    for i in global_adjacencies(M, idx, ignored_elem)
        if i == adj_elem
            n += 1
        end
    end

    return n
end
