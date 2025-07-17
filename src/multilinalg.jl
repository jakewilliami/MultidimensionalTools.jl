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

# LinearAlgebra's definition of `det`:
#=
# https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/dense.jl
function inv(A::StridedMatrix{T}) where T
    checksquare(A)
    S = typeof((one(T)*zero(T) + one(T)*zero(T))/one(T))
    AA = convert(AbstractArray{S}, A)
    if istriu(AA)
        Ai = triu!(parent(inv(UpperTriangular(AA))))
    elseif istril(AA)
        Ai = tril!(parent(inv(LowerTriangular(AA))))
    else
        Ai = inv!(lu(AA))
        Ai = convert(typeof(parent(Ai)), Ai)
    end
    return Ai
end

=#
# code lowered:
#=
julia> @code_typed inv(A)
CodeInfo(
1 ── %1   = Base.arraysize(A, 1)::Int64
│    %2   = Base.arraysize(A, 2)::Int64
│    %3   = (%1 === %2)::Bool
└───        goto #3 if not %3
2 ──        goto #4
3 ── %6   = Base.arraysize(A, 1)::Int64
│    %7   = Base.arraysize(A, 2)::Int64
│    %8   = Core.tuple(%6, %7)::Tuple{Int64,Int64}
│    %9   = invoke Base.print_to_string("matrix is not square: dimensions are "::String, %8::Vararg{Any,N} where N)::String
│    %10  = %new(Base.DimensionMismatch, %9)::DimensionMismatch
│           LinearAlgebra.throw(%10)::Union{}
└───        $(Expr(:unreachable))::Union{}
4 ┄─ %13  = invoke AbstractArray{Float64,2}(_2::Array{Int8,2})::Array{Float64,2}
│    %14  = LinearAlgebra.istriu::Core.Compiler.Const(LinearAlgebra.istriu, false)
│    %15  = invoke %14(%13::Array{Float64,2}, 0::Int64)::Bool
└───        goto #23 if not %15
5 ──        Base.arraysize(%13, 1)::Int64
│           Base.arraysize(%13, 2)::Int64
│    %19  = Base.or_int(false, false)::Bool
│    %20  = Base.or_int(%19, false)::Bool
│    %21  = Base.not_int(%20)::Bool
└───        goto #7 if not %21
6 ──        goto #8
7 ── %24  = %new(Core.ArgumentError, "offset arrays are not supported but got an array with index other than 1")::ArgumentError
│           Base.throw(%24)::Union{}
└───        $(Expr(:unreachable))::Union{}
8 ┄─ %27  = Base.arraysize(%13, 1)::Int64
│    %28  = Base.arraysize(%13, 2)::Int64
│    %29  = (%27 === %28)::Bool
└───        goto #10 if not %29
9 ──        goto #11
10 ─ %32  = Base.arraysize(%13, 1)::Int64
│    %33  = Base.arraysize(%13, 2)::Int64
│    %34  = Core.tuple(%32, %33)::Tuple{Int64,Int64}
│    %35  = invoke Base.print_to_string("matrix is not square: dimensions are "::String, %34::Vararg{Any,N} where N)::String
│    %36  = %new(Base.DimensionMismatch, %35)::DimensionMismatch
│           LinearAlgebra.throw(%36)::Union{}
└───        $(Expr(:unreachable))::Union{}
11 ┄        goto #12
12 ─        goto #13
13 ─ %41  = Base.arraysize(%13, 1)::Int64
│    %42  = Base.arraysize(%13, 1)::Int64
│    %43  = LinearAlgebra.I::UniformScaling{Bool}
│    %44  = Core.tuple(%41, %42)::Tuple{Int64,Int64}
│    %45  = invoke Array{Float64,2}(%43::UniformScaling{Bool}, %44::Tuple{Int64,Int64})::Array{Float64,2}
│    %46  = LinearAlgebra.LAPACK.trtrs!::typeof(LinearAlgebra.LAPACK.trtrs!)
│    %47  = invoke %46('U'::Char, 'N'::Char, 'N'::Char, %13::Array{Float64,2}, %45::Array{Float64,2})::Array{Float64,2}
│           Base.arraysize(%47, 1)::Int64
│           Base.arraysize(%47, 2)::Int64
│    %50  = Base.or_int(false, false)::Bool
│    %51  = Base.or_int(%50, false)::Bool
│    %52  = Base.not_int(%51)::Bool
└───        goto #15 if not %52
14 ─        goto #16
15 ─ %55  = %new(Core.ArgumentError, "offset arrays are not supported but got an array with index other than 1")::ArgumentError
│           Base.throw(%55)::Union{}
└───        $(Expr(:unreachable))::Union{}
16 ┄ %58  = Base.arraysize(%47, 1)::Int64
│    %59  = Base.arraysize(%47, 2)::Int64
│    %60  = (%58 === %59)::Bool
└───        goto #18 if not %60
17 ─        goto #19
18 ─ %63  = Base.arraysize(%47, 1)::Int64
│    %64  = Base.arraysize(%47, 2)::Int64
│    %65  = Core.tuple(%63, %64)::Tuple{Int64,Int64}
│    %66  = invoke Base.print_to_string("matrix is not square: dimensions are "::String, %65::Vararg{Any,N} where N)::String
│    %67  = %new(Base.DimensionMismatch, %66)::DimensionMismatch
│           LinearAlgebra.throw(%67)::Union{}
└───        $(Expr(:unreachable))::Union{}
19 ┄        goto #20
20 ─        goto #21
21 ─        goto #22
22 ─ %73  = invoke LinearAlgebra.triu!(%47::Array{Float64,2}, 0::Int64)::Array{Float64,2}
└───        goto #56
23 ─ %75  = LinearAlgebra.istril::Core.Compiler.Const(LinearAlgebra.istril, false)
│    %76  = invoke %75(%13::Array{Float64,2}, 0::Int64)::Bool
└───        goto #42 if not %76
24 ─        Base.arraysize(%13, 1)::Int64
│           Base.arraysize(%13, 2)::Int64
│    %80  = Base.or_int(false, false)::Bool
│    %81  = Base.or_int(%80, false)::Bool
│    %82  = Base.not_int(%81)::Bool
└───        goto #26 if not %82
25 ─        goto #27
26 ─ %85  = %new(Core.ArgumentError, "offset arrays are not supported but got an array with index other than 1")::ArgumentError
│           Base.throw(%85)::Union{}
└───        $(Expr(:unreachable))::Union{}
27 ┄ %88  = Base.arraysize(%13, 1)::Int64
│    %89  = Base.arraysize(%13, 2)::Int64
│    %90  = (%88 === %89)::Bool
└───        goto #29 if not %90
28 ─        goto #30
29 ─ %93  = Base.arraysize(%13, 1)::Int64
│    %94  = Base.arraysize(%13, 2)::Int64
│    %95  = Core.tuple(%93, %94)::Tuple{Int64,Int64}
│    %96  = invoke Base.print_to_string("matrix is not square: dimensions are "::String, %95::Vararg{Any,N} where N)::String
│    %97  = %new(Base.DimensionMismatch, %96)::DimensionMismatch
│           LinearAlgebra.throw(%97)::Union{}
└───        $(Expr(:unreachable))::Union{}
30 ┄        goto #31
31 ─        goto #32
32 ─ %102 = Base.arraysize(%13, 1)::Int64
│    %103 = Base.arraysize(%13, 1)::Int64
│    %104 = LinearAlgebra.I::UniformScaling{Bool}
│    %105 = Core.tuple(%102, %103)::Tuple{Int64,Int64}
│    %106 = invoke Array{Float64,2}(%104::UniformScaling{Bool}, %105::Tuple{Int64,Int64})::Array{Float64,2}
│    %107 = LinearAlgebra.LAPACK.trtrs!::typeof(LinearAlgebra.LAPACK.trtrs!)
│    %108 = invoke %107('L'::Char, 'N'::Char, 'N'::Char, %13::Array{Float64,2}, %106::Array{Float64,2})::Array{Float64,2}
│           Base.arraysize(%108, 1)::Int64
│           Base.arraysize(%108, 2)::Int64
│    %111 = Base.or_int(false, false)::Bool
│    %112 = Base.or_int(%111, false)::Bool
│    %113 = Base.not_int(%112)::Bool
└───        goto #34 if not %113
33 ─        goto #35
34 ─ %116 = %new(Core.ArgumentError, "offset arrays are not supported but got an array with index other than 1")::ArgumentError
│           Base.throw(%116)::Union{}
└───        $(Expr(:unreachable))::Union{}
35 ┄ %119 = Base.arraysize(%108, 1)::Int64
│    %120 = Base.arraysize(%108, 2)::Int64
│    %121 = (%119 === %120)::Bool
└───        goto #37 if not %121
36 ─        goto #38
37 ─ %124 = Base.arraysize(%108, 1)::Int64
│    %125 = Base.arraysize(%108, 2)::Int64
│    %126 = Core.tuple(%124, %125)::Tuple{Int64,Int64}
│    %127 = invoke Base.print_to_string("matrix is not square: dimensions are "::String, %126::Vararg{Any,N} where N)::String
│    %128 = %new(Base.DimensionMismatch, %127)::DimensionMismatch
│           LinearAlgebra.throw(%128)::Union{}
└───        $(Expr(:unreachable))::Union{}
38 ┄        goto #39
39 ─        goto #40
40 ─        goto #41
41 ─ %134 = invoke LinearAlgebra.tril!(%108::Array{Float64,2}, 0::Int64)::Array{Float64,2}
└───        goto #56
42 ─ %136 = $(Expr(:foreigncall, :(:jl_array_copy), Ref{Array{Float64,2}}, svec(Any), 0, :(:ccall), :(%13)))::Array{Float64,2}
│    %137 = LinearAlgebra.lu!::typeof(lu!)
│    %138 = Base.sle_int(1, 1)::Bool
└───        goto #44 if not %138
43 ─ %140 = Base.sle_int(1, 0)::Bool
└───        goto #45
44 ─        nothing::Nothing
45 ┄ %143 = φ (#43 => %140, #44 => false)::Bool
└───        goto #47 if not %143
46 ─        invoke Base.getindex(()::Tuple, 1::Int64)::Union{}
└───        $(Expr(:unreachable))::Union{}
47 ┄        goto #48
48 ─        goto #49
49 ─        goto #50
50 ─        goto #51
51 ─ %151 = invoke LinearAlgebra.:(var"#lu!#132")(true::Bool, %137::typeof(lu!), %136::Array{Float64,2}, $(QuoteNode(Val{true}()))::Val{true})::LU{Float64,Array{Float64,2}}
└───        goto #52
52 ─        goto #53
53 ─        goto #54
54 ─        goto #55
55 ─ %156 = LinearAlgebra.LAPACK.getri!::typeof(LinearAlgebra.LAPACK.getri!)
│    %157 = LinearAlgebra.getfield(%151, :factors)::Array{Float64,2}
│           Base.arraysize(%157, 1)::Int64
│           Base.arraysize(%157, 2)::Int64
│    %160 = LinearAlgebra.getfield(%151, :factors)::Array{Float64,2}
│    %161 = LinearAlgebra.getfield(%151, :factors)::Array{Float64,2}
│           Base.arraysize(%161, 1)::Int64
│           Base.arraysize(%161, 2)::Int64
│    %164 = LinearAlgebra.getfield(%151, :ipiv)::Array{Int64,1}
└─── %165 = invoke %156(%160::Array{Float64,2}, %164::Array{Int64,1})::Array{Float64,2}
56 ┄ %166 = φ (#22 => %73, #41 => %134, #55 => %165)::Array{Float64,2}
└───        return %166
) => Array{Float64,2}
=#

#=
# *(A::AbstractArray{T,2} where T, B::AbstractArray{T,2} where T) in LinearAlgebra at /Applications/Julia-1.5.app/Contents/Resources/julia/share/julia/stdlib/v1.5/LinearAlgebra/src/matmul.jl:151
function (*)(A::AbstractMatrix, B::AbstractMatrix)
    TS = promote_op(matprod, eltype(A), eltype(B))
    mul!(similar(B, TS, (size(A,1), size(B,2))), A, B)
end
=#

function get_I(S::NTuple{N, T}) where {T <: Integer, N}
    A = spzeros(S...)
    for i in CartesianIndices(A)
        if allequal(i.I)
            A[i] = 1
        end
    end

    return A
end

function _promote_I(S::NTuple{N, T}) where {T <: Integer, N}
    return T.(reshape(
        repeat(I(first(S)), inner = ntuple(i -> i == 2 ? first(S) : 1, first(S))),
        ntuple(_ -> first(S), N),
    ))
end

function LinearAlgebra.I(S::NTuple{N, T}; dims::Int = 2) where {T <: Integer, N}
    allequal(S) || (s = rand(S);
    throw(DimensionMismatch("""
    The identity array should be square.  Size received: $(size(A)).  Example size: $(ntuple(_ -> s, length(A))).
""")))

    dims == 2 && return _promote_I(S)
    return I(S, dims = (dims - 1))

    # return reshape(repeat(I(first(S)), inner = ntuple(i -> i == 2 ? first(S) : 1, first(S))), ntuple(_ -> first(S), N))

    # reshape(repeat(I(3), inner = (1, 3, 1)), (3, 3, 3))

    # cat(Iterators.repeated(I(3), 3)..., dims = 1)
    #
    # sz = ntuple(d -> d == dims ? n : size(M, d), max(N, dims))
    # return cat(M, fill(fill_elem, sz); dims = dims)
end
LinearAlgebra.I(A::Integer...) = LinearAlgebra.I(A)

LinearAlgebra.I(N::Integer...) =
    function _inv_nD(A::AbstractArray{T, N}) where {T <: Real, N}
        issquare(A) || (s = rand(size(A));
        throw(DimensionMismatch("""
    The inverse is defined as the matrix required to multiply with the input matrix to produce an identity matrix.  The identity matrix is square, so you need to ensure every dimension of the input array has the same size.  Size received: $(size(A)).  Example size: $(ntuple(_ -> s, ndims(A))).
""")))
    end

function LinearAlgebra.inv(A::AbstractArray{T, N}) where {T <: Real, N}
    N < 3 && return inv(A)
    return _inv_nD(A)
end

import Base: *
@generated function *(A::AbstractArray{T, N}, B::AbstractArray{T, N}) where {T, N}
    N < 3 && return :(A * B)
    C = transpose(B)
    return :(@matmul M[i, j] := sum(k) A[i, k] * C[j, k])
end
