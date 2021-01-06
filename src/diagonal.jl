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
