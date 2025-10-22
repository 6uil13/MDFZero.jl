module MDFZero

import ILUZero
import LinearAlgebra: Symmetric, LowerTriangular, UpperTriangular, I, ldiv!, \, norm, tril, triu, one
import SparseArrays: SparseMatrixCSC, permute, spdiagm
import Base: permute!

include("ILU0.jl")
include("Laplacian.jl")

export mdf0, mdf0!, ILU0, laplacian, distance, permute, permute!

function update!(A::SparseMatrixCSC{T, <:Integer}, m::Int, a::Vector{T}) where T <: Real
    colptrm = A.colptr[m]:A.colptr[m + 1] - 1
    
    a .= 0.0
    Amm = 0.0
    for r in colptrm
        i = A.rowval[r]
        a[r - colptrm[1] + 1] = A.nzval[r]
        if i == m
            Amm = a[r - colptrm[1] + 1]
        end
    end

    for r in colptrm
        i = A.rowval[r]
        Aim = a[r - colptrm[1] + 1]

        𝓘 = A.colptr[i]:A.colptr[i + 1] - 1
        rows = view(A.rowval, 𝓘)
        for s in colptrm
            j = A.rowval[s]
            # find Aji = Aij != 0
            t = findfirst(isequal(j), rows)

            if !isnothing(t)
                Ajm = a[s - colptrm[1] + 1] # = Amj
                A.nzval[t + 𝓘[1] - 1] -= (Aim * Ajm) / Amm
            end
        end
    end

    nothing
end

function discardedfill(A::SparseMatrixCSC{<:Real, <:Integer}, m::Int)
    f = 0.0
    defficiency = 0
    degree = 0

    Amm = 0.0
    colptrm = A.colptr[m]:A.colptr[m + 1] - 1
    for r in colptrm
        i = A.rowval[r]
        Aim = A.nzval[r]
        if i == m
            Amm = Aim
        end

        𝓘 = A.colptr[i]:A.colptr[i + 1] - 1
        rows = view(A.rowval, 𝓘)
        for s in colptrm
            j = A.rowval[s]
            t = findfirst(isequal(j), rows)
            null_Aij = ! (!isnothing(t) &&
            abs(A.nzval[t + 𝓘[1] - 1]) > eps(1e2))

            if null_Aij
                Ajm = A.nzval[s]
                f += (Aim * Ajm)^2
                defficiency += 1
            else
                degree += 1
            end
        end
    end

    f /= Amm^2

    f, defficiency, degree, m
end

Base.@propagate_inbounds function mdf0!(S::Symmetric{T, SparseMatrixCSC{T, U}}, 
    discard::AbstractVector{T} = T[]) where {T <: Real, U <: Integer}

    A = S.data
    n = size(A, 1)

    colnnz = maximum(diff(A.colptr))
    a = zeros(T, colnnz)

    # 4 measures: discard, defficiency, degree, label
    fillin = Dict{Int, Tuple{T, Int, Int, Int}}()

    # initial discard
    for m in 1:n
        fillin[m] = discardedfill(A, m)
    end

    # main decomposition loop
    σ = zeros(Int, n)
    for k = 1:n-1
        m = argmin(fillin)
        σ[k] = m
        if !isempty(discard)
            discard[m] = √fillin[m][1]
        end

    	update!(A, m, a)
        
        pop!(fillin, m)
        𝓝 = A.rowval[A.colptr[m]:A.colptr[m + 1] - 1]
        for v in 𝓝
            if haskey(fillin, v)
                fillin[v] = discardedfill(A, v)
            end
        end
    end
    σ[n] = argmin(fillin)
    
    σ
end

mdf0(A, discard = eltype(A)[]) = mdf0!(copy(A), discard)

end
