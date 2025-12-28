module MultidimensionalTools

using LinearAlgebra

# type wrappers
export AbstractIndex, AbstractIndices, AbstractIndexOrIndices
# helper functions
export n_adjacencies, n_faces, extrema_indices, tryindex, get_directions
# main functions
export append_n_times,
    append_n_times_backwards,
    promote_to_nD,
    promote_to_3D,
    expand_as_required,
    adjacencies,
    expanded_adjacencies,
    n_adjacent_to,
    expanded_n_adjacent_to,
    global_adjacencies_indices,
    global_adjacencies,
    global_n_adjacent_to,
    reduce_dims

_IterableList{T} = Union{AbstractVector{T}, NTuple{N, T}} where {N}
AbstractIndex{N} = Union{NTuple{N, T}, CartesianIndex{N}} where {T <: Integer}
AbstractIndices{N} = Union{
    AbstractVector{I}, NTuple{M, I}, CartesianIndices{N, NTuple{N, Base.OneTo{T}}}
} where {I <: AbstractIndex{N}, T <: Integer, M}
AbstractIndexOrIndices{N} = Union{AbstractIndex{N}, AbstractIndices{N}}

"""
    𝟘(n::Integer) -> NTuple{n, Integer}
    \\bbzero<tab>(n::Integer) -> NTuple{n, Integer}

Simply returns the zero tuple of a given dimension.
"""
𝟘(n::T) where {T <: Integer} = ntuple(_ -> zero(T), n)

_tryindex(M::AbstractArray, inds::CartesianIndex...) =
    Tuple((checkbounds(Bool, M, i) ? M[i] : nothing) for i in inds)
_tryindex(M::AbstractArray, inds::AbstractArray{CartesianIndex}) = _tryindex(M, inds...)
_tryindex(M::AbstractArray, inds::NTuple...) = _tryindex(M, CartesianIndex.(inds)...)
_tryindex(M::AbstractArray, inds::AbstractArray{NTuple}) = _tryindex(M, inds...)

"""
    tryindex(M::AbstractArray{T, N}, inds::AbstractIndex{N}}...)
    tryindex(M::AbstractArray{T, N}, inds::AbstractIndices{N}})

For each index specified, gets the index or missing if the index is unavailable.
"""
tryindex(M::AbstractArray, inds::AbstractIndex...) = _tryindex(M, inds...)
tryindex(M::AbstractArray, inds::AbstractIndices) = _tryindex(M, inds...)

function _extrema_indices(I::CartesianIndex...)
    Tuple(extrema(reduce(vcat, permutedims(collect(i)) for i in I...), dims = 1))
end
function _extrema_indices(I::AbstractIndexOrIndices{N}...) where {N}
    return Tuple(extrema(reduce(vcat, permutedims(collect(i)) for i in I...), dims = 1))
end

"""
    extrema_indices(I::AbstractIndices{N}) -> NTuple{N, NTuple{N, Integer}}
    extrema_indices(i::AbstractIndex{N}...) -> NTuple{N, NTuple{N, Integer}}
    extrema_indices(A::AbstractArray{T, N}...) -> NTuple{N, NTuple{N, Integer}}

Given indices or arrays, returns an `NTuple`.  Each element in the `NTuple` represents a (min, max) tuple for each dimension.

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
extrema_indices(I::AbstractIndices{N}) where {N} = _extrema_indices(I)
extrema_indices(i::AbstractIndex{N}...) where {N} = _extrema_indices(i)
extrema_indices(A::AbstractArray{T, N}) where {T, N} = size(A)
extrema_indices(A::AbstractArray{T, N}...) where {T, N} = _extrema_indices(size.(A))

@inline function allequal(A::_IterableList{T}) where {T}
    length(A) < 2 && return true

    @inbounds for i in 2:length(A)
        first(A) ≠ A[i] && return false
    end

    return true
end
@inline allequal(a::T...) where {T} = allequal(a)

"""
Check whether the given array is a cube (i.e., all its dimensions are the same size).  Despite the function name, this will work for dimension < 3
"""
ishypercube(::Number) = true # this is the case where there is zero dimensions: the point in space
ishypercube(A::AbstractVector{T}) where {T} = length(A) == 1 ? true : false # a vector of length one is a square, hence it is a lower-dimensional cube
ishypercube(A::AbstractArray{T, N}) where {T, N} = allequal(size(A))
function ishypercube_alt(A::AbstractArray{T, N}) where {T, N}
    ax1 = axes(A, 1)
    for d in 2:ndims(A)
        axes(A, d) == ax1 || return false
    end
    return true
end
ishypercube(A::AbstractArray{T, N}...) where {T, N} = Bool[ishypercube(a) for a in A]

"""
Check whether the 2-dimensional aspect of a given array is a square.
"""
is2Dsquare(::Union{Number, AbstractVector{T}}) where {T} = false # there is no 2-dimensional aspect of zero- or one-dimensional objects
is2Dsquare(A::AbstractArray{T, N}) where {T, N} = size(A, 1) == size(A, 2) ? true : false
is2Dsquare(A::AbstractArray{T, N}...) where {T, N} = Bool[is2Dsquare(a) for a in A]

include("generic.jl")
include("adjacent.jl")
# include("multilinalg.jl")

end # end module
