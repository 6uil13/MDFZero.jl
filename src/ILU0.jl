struct ILU{T <: Any, N <: Integer}
    L::LowerTriangular{T}
    U::UpperTriangular{T}
    p::AbstractVector{N}
    q::AbstractVector{N}
end

ILU(A::Symmetric) = ILU(A.data)
ILU(A::Symmetric, p) = ILU(A.data, p, p)

function tril(F::ILUZero.ILU0Precon)
    LowerTriangular(I + SparseMatrixCSC(F.m, min(F.m, F.n), F.l_colptr, F.l_rowval, F.l_nzval))
end

function triu(F::ILUZero.ILU0Precon)
    UpperTriangular(SparseMatrixCSC(min(F.m, F.n), F.n, F.u_colptr, F.u_rowval, F.u_nzval))
end

Base.iterate(F::ILUZero.ILU0Precon, args...) = iterate((tril(F), triu(F)), args...)

function ILU(F::ILUZero.ILU0Precon,
    p::AbstractVector{<:Integer} = 1:F.n,
    q::AbstractVector{<:Integer} = 1:F.m)

    L, U = F
    
    ILU(L, U, p, q)
end

ILU(A::SparseMatrixCSC,
p::AbstractVector{<:Integer} = 1:size(A, 1),
q::AbstractVector{<:Integer} = 1:size(A, 2)) = ILU(ILUZero.ilu0(permute(A, p, q)), p, q)

function ldiv!(Y, A::ILU, B)
    permute!(B, A.p)
    ldiv!(Y, A.L, B)
    invpermute!(B, A.q)
    ldiv!(A.U, Y)
    invpermute!(Y, A.q)
end

function \(A::ILU, B)
    X = zero(B)
    ldiv!(X, A, B)
    
    X
end

Base.iterate(F::ILU, args...) = iterate((F.L, F.U, F.p, F.q), args...)

function distance(A::Symmetric{T, SparseMatrixCSC{T, V}},
                    F::ILU{T, V}) where {T <: Real, V <: Integer}
    
    L, U, p = F
    norm(permute(A, p, p) - L * U)
end

distance(A, F) = distance(F, A)