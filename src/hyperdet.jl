using LinearAlgebra

# LinearAlgebra's definition of `det`:
#=
function det(A::AbstractMatrix{T}) where T
    if istriu(A) || istril(A) # is triangular upper or is triangular lower
        S = typeof((one(T)*zero(T) + zero(T))/one(T))
        return convert(S, det(UpperTriangular(A)))
    end
    return det(lu(A; check = false)) # LU factorisation: https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/lu.jl
end

det(x::Number) = x

# https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/lu.jl#L442
function det(F::LU{T}) where T
    n = checksquare(F)
    issuccess(F) || return zero(T)
    P = one(T)
    c = 0
    @inbounds for i = 1:n
        P *= F.factors[i,i]
        if F.ipiv[i] != i
            c += 1
        end
    end
    s = (isodd(c) ? -one(T) : one(T))
    return P * s
end
=#

ishypercube(::Number) && return true # this is the case where there is zero dimensions: the point in space
function ishypercube(A::AbstractArray{T, N}) where {T, N}
    
end

function _first_cayley_hyperdet(A::AbstractArray{T, N}) where {T <: Real, N} # usually denoted "det₀"

end

function _second_cayley_hyperdet(A::AbstractArray{T, N}) where {T <: Real, N} # usually denoted "Det."
    
end

function LinearAlgebra.det(A::AbstractArray{T, N}) where {T <: Real, N}
    N < 3 && return det(A)
    iseven(ndims(A)) && return _first_cayley_hyperdet(A)
    # size(A) == ntuple(_ -> 2, ndims(A)) && return _first_cayley_hyperdet(A)
    
    return "hello"
end
