struct ILU0{T <: Any, N <: Integer}
    L::LowerTriangular{T}
    U::UpperTriangular{T}
    p::AbstractVector{N}
    q::AbstractVector{N}
end

ILU0(A::Symmetric) = ILU0(A.data)
ILU0(A::Symmetric, p) = ILU0(A.data, p, p)

function tril(F::ILUZero.ILU0Precon)
    LowerTriangular(I + SparseMatrixCSC(F.m, min(F.m, F.n), F.l_colptr, F.l_rowval, F.l_nzval))
end

function triu(F::ILUZero.ILU0Precon)
    UpperTriangular(SparseMatrixCSC(min(F.m, F.n), F.n, F.u_colptr, F.u_rowval, F.u_nzval))
end

Base.iterate(F::ILUZero.ILU0Precon, args...) = iterate((tril(F), triu(F)), args...)

function ILU0(F::ILUZero.ILU0Precon,
    p::AbstractVector{<:Integer} = 1:F.n,
    q::AbstractVector{<:Integer} = 1:F.m)

    L, U = F
    
    ILU0(L, U, p, q)
end

ILU0(A::SparseMatrixCSC,
p::AbstractVector{<:Integer} = 1:size(A, 1),
q::AbstractVector{<:Integer} = 1:size(A, 2)) = ILU0(ILUZero.ilu0(permute(A, p, q)), p, q)

function ldiv!(Y, A::ILU0, B)
    permute!(B, A.p)
    ldiv!(Y, A.L, B)
    invpermute!(B, A.q)
    ldiv!(A.U, Y)
    invpermute!(Y, A.q)
end

function \(A::ILU0, B)
    X = zero(B)
    ldiv!(X, A, B)
    
    X
end

Base.iterate(F::ILU0, args...) = iterate((F.L, F.U, F.p, F.q), args...)

function distance(A::Symmetric{T, SparseMatrixCSC{T, V}},
                    F::ILU0{T, V}) where {T <: Real, V <: Integer}
    
    L, U, p = F
    norm(permute(A, p, p) - L * U)
end

distance(A, F) = distance(F, A)