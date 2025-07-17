# struct Diagonal{T, A <: AbstractVector{T}} <: AbstractArray{T}
#     diag::A
#
#     function Diagonal{T, A}(diag) where {T, A <: AbstractArray{T}}
#         Base.require_one_based_indexing(diag)
#         return new{T, A}(diag)
#     end
# end
# Diagonal(a::AbstractArray{T, N}) where {T, N} = Diagonal{T, typeof(a)}(a)
# Diagonal{T, N}(a::AbstractArray) where {T, N} = Diagonal(convert(AbstractArray{T, N}, a)::AbstractArray{T, N})

function LinearAlgebra.diagind(A::AbstractArray)
    return CartesianIndex{ndims(A)}[i for i in CartesianIndices(A) if allequal(i.I)]
end

LinearAlgebra.diag(A::AbstractArray) = A[diagind(A)]

using Einsum

#=
If I were to extend LinearAlgebra.diag for an arbitrary dimension, what would it look like?  Would it be all of the indices that are the same (thus returning a vector)?  I.e., (1, 1, …, 1), (2, 2, …, 2), …, (n, n, …, n)?  Or would it be a diagonal hyperplane that is one dimension less than the given matrix?
=#

"""
```julia
diag_alt(A::AbstractArray)
```

Given an array of arbitrary dimesions, select all of the indices that are the same, thus returning a vector.  E.g.:
```julia
julia> A = rand(Int8, 3, 3)
3×3 Array{Int8,2}:
 93  123  -32
 71  -68  -75
 99   -7  -57

julia> diag_alt(A)
3-element Array{Int8,1}:
  93
 -68
 -57

julia> B = rand(Int8, 3, 3, 3)
3×3×3 Array{Int8,3}:
[:, :, 1] =
  16   -23    53
  13    -7  -115
 -96  -105    29

[:, :, 2] =
 -80  17    53
  88  40    15
  37  51  -109

[:, :, 3] =
  -81    88   37
   38  -115  -29
 -112    73  -41

julia> diag_alt(B)
3-element Array{Int8,1}:
  16
  40
 -41
```
"""
@generated function diag_alt(A::AbstractArray)
    I = Symbol[:i for _ in 1:ndims(A)]
    return :(@einsum B[i] := A[$(I...)]) # e.g., @einsum y[i] := x[i,i,i]
end

"""
```julia
diag_hyperplaneplane(A::AbstractArray)
```

Given an array of arbitrary dimesions, return the diagonal hyperplane.
"""
@generated function diag_hyperplaneplane(A::AbstractArray{T, N}) where {T, N}
    N < 3 && return :(diag_alt(A))
    I = vcat(Symbol[:i for _ in 1:(ndims(A) - 1)], :k)
    K = vcat(Symbol[:i for _ in 1:(ndims(A) - 2)], :k)
    return :(@einsum B[$(K...)] := A[$(I...)]) # e.g., @einsum y[i,k] := x[i,i,k]
end

function diagind_alt(A::AbstractArray)
    ishypercube(A) ||
        throw(DimensionMismatch("Error trying to access the diagonal indices from a non-diagonalisable array (i.e., input array is not a hypercube)"))
    return map(i -> CartesianIndex(ntuple(d -> i, ndims(A))), axes(A, 1))
end

function diag_alt_alt(A::AbstractArray)
    return A[diagind_alt(A)]
end

#=
julia> @pretty @matmul M[i,j] := sum(k) A[i, k] * B[j, k]
begin
    local lobster = permutedims(B, (2, 1))
    local boar = A * lobster
    M = boar
end
julia> @pretty @reduce M[i,j] := sum(k) A[i, k] * B[j, k]
begin
    local lobster = orient(PermuteDims(B), (*, :, :))
    M = dropdims(sum(@__dot__(A * lobster), dims = 2), dims = 2)
end
=#
