## NOTE:
# n-dim sparse array?
# each code word in CodingTheory is a point in an n-dim array

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
Repeats a specified value n many times along the specified dimension.  If no `fill_elem` is given, and the element type of the matrix is a number, fills with zero.

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
append_n_times(M::AbstractArray{T, N}, n::Integer, dims::Integer = 1) where {T <: Number, N} = append_n_times(M, n, zero(T); dims = dims)

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
append_n_times_backwards(M::AbstractArray{T, N}, n::Integer, dims::Integer = 1) where {T <: Number, N} = append_n_times(M, n, zero(T); dims = dims)

"""
    function promote_to_nD(M::AbstractArray{T, N}, n::Integer, fill_elem::T)

Assumes the given matrix M is an (n - 1) dimensional slice of an n-dimensional structure, and fills in the array to n dimensions with `fill_elem`.  If no `fill_elem` is given, and the element type of the matrix is a number, fills with zero.
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
promote_to_nD(M::AbstractArray{T, N}, n::Integer) where {T <: Number, N} = promote_to_nD(M, n, zero(T))

"""
    promote_to_3D(M::AbstractArray{T, N}, fill_elem::T)
    promote_to_3D(M::AbstractArray{T, N}) where T <: Number
    
A simple, self-evident wrapper for `promote_to_nD`.   No `fill_elem` assumes that the value is zero (if it can be).
"""
promote_to_3D(M::AbstractArray{T, N}, fill_elem::T) where {T, N} = promote_to_nD(M, 3, fill_elem)
promote_to_3D(M::AbstractArray{T, N}) where {T <: Number, N} = promote_to_3D(M, zero(T))

"""
    reduce_dim(A::AbstractArray{T, N}; reduce_by::Int = 1)

Reduce an `N`-dimensional array `A` to dimensions `N - reduce_by`.  Can reduce arbitrarily.  Currently unsafe (i.e., doesn't check that `reduce_by < ndims(A)`.

```julia
julia> A = zeros(Int, 3, 3, 3);

julia> reduce_dims(A)
3×3 Matrix{Vector{Int64}}:
 [0, 0, 0]  [0, 0, 0]  [0, 0, 0]
 [0, 0, 0]  [0, 0, 0]  [0, 0, 0]
 [0, 0, 0]  [0, 0, 0]  [0, 0, 0]
```
"""
function reduce_dims(A::AbstractArray{T, N}; reduce_by::Int = 1) where {T, N}
    sz = size(A)
    out_A = Array{Array{T, reduce_by}, N - reduce_by}(undef, sz[1:(end - reduce_by)])
    for i in CartesianIndices(sz[1:(end - reduce_by)])
        out_A[i] = A[i, (Colon() for _ in reduce_by)...]
    end
    return out_A
end
