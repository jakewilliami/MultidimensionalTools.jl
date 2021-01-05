"""
    AbstractIndex{N} = Union{NTuple{N, T}, CartesianIndex{N}}
"""
AbstractIndex{N} = Union{NTuple{N, T}, CartesianIndex{N}} where {T <: Integer}

"""
    AbstractIndices{N}  = Union{AbstractArray{I, M}, NTuple{M, I}, CartesianIndices{N, NTuple{N, Base.OneTo{T}}}} where I <: AbstractIndex{N}
"""
AbstractIndices{N}  = Union{AbstractArray{I, M}, NTuple{M, I}, CartesianIndices{N, NTuple{N, Base.OneTo{T}}}} where {I <: AbstractIndex{N}, T <: Integer, M}

"""
    AbstractIndexOrIndices{N} = Union{AbstractIndex{N}, AbstractIndices{N}}
"""
AbstractIndexOrIndices{N} = Union{AbstractIndex{N}, AbstractIndices{N}}

"""
    𝟘(n::Integer) -> NTuple{n, Integer}
    \\bbzero<tab>(n::Integer) -> NTuple{n, Integer}

Simply returns the zero tuple of a given dimension.
"""
𝟘(n::T) where {T <: Integer} = ntuple(_ -> zero(T), n)

_tryindex(M::AbstractArray, inds::CartesianIndex...) where {T, N} = Tuple((checkbounds(Bool, M, i) ? M[i] : nothing) for i in inds)
_tryindex(M::AbstractArray, inds::AbstractArray{CartesianIndex}) where {T, N} = _tryindex(M, inds...)
_tryindex(M::AbstractArray, inds::NTuple...) = _tryindex(M, CartesianIndex.(inds)...)
_tryindex(M::AbstractArray, inds::AbstractArray{NTuple}) where {T, N} = _tryindex(M, inds...)

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
extrema_indices(I::AbstractIndices{N}) where {N} = _extrema_indices(I)
extrema_indices(i::AbstractIndex{N}...) where {N} = _extrema_indices(i)
