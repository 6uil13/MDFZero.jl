module MDF

using SparseArrays, LinearAlgebra
export mdf, mdf!

function update!(A::SparseMatrixCSC{T, Int}, m::Int) where T <: Real
    colptrm = A.colptr[m]:A.colptr[m + 1] - 1
    a = A[:, m]

    amm = a[m]
    for r in colptrm
        i = A.rowval[r]
        aim = a[i]

        𝓘 = A.colptr[i]:A.colptr[i + 1] - 1
        rows = view(A.rowval, 𝓘)
        for s in colptrm
            j = A.rowval[s]
            # find aji = aij != 0
            t = findfirst(rows .== j)

            if !isnothing(t)
                ajm = a[j] # = amj
                A.nzval[t + 𝓘[1] - 1] -= (aim * ajm) / amm
            end
        end
    end

    nothing
end

function discardedfill(A::SparseMatrixCSC{T, Int}, m::Int) where T <: Real
    amm = A[m, m]
    f = 0.0
    defficiency = 0
    degree = 0
    for r = A.colptr[m]:A.colptr[m + 1] - 1
        i = A.rowval[r]
        aim = A.nzval[r]

        𝓘 = A.colptr[i]:A.colptr[i + 1] - 1
        rows = view(A.rowval, 𝓘)
        for s = A.colptr[m]:A.colptr[m + 1] - 1
            j = A.rowval[s]
            t = findfirst(rows .== j)
            null_aij = ! (!isnothing(t) &&
            abs(A.nzval[t + 𝓘[1] - 1]) > eps(1e2))

            if null_aij
                ajm = A.nzval[s]
                f += (aim * ajm)^2
                defficiency += 1
            else
                degree += 1
            end
        end
    end

    f /= amm^2

    f, defficiency, degree, m
end

Base.@propagate_inbounds function mdf!(S::Symmetric{T, SparseMatrixCSC{T, Int}};
              finaldiscard::Bool = false) where T <: Real

    A = S.data
    n = size(A, 1)

    # 4 measures: discard, defficiency, degree, label
    fillin = Dict{Int, Tuple{T, T, T, Int}}()

    # initial discard
    for m in 1:n
        fillin[m] = discardedfill(A, m)
    end

    # main decomposition loop
    discard = zeros(n)
    σ = zeros(Int, n)
    for k = 1:n-1
        m = argmin(fillin)
        σ[k] = m
        if finaldiscard
            discard[m] = fillin[m][1]
        end

    	update!(A, m)
        
        pop!(fillin, m)
        𝓝 = A.rowval[A.colptr[m]:A.colptr[m + 1] - 1]
        for v in 𝓝
            if !haskey(fillin, v)
                continue
            end
            fillin[v] = discardedfill(A, v)
        end
    end
    σ[n] = argmin(fillin)
    
    if finaldiscard
        return σ, .√discard
    else
        return σ
    end
end

mdf(A; finaldiscard::Bool = false) = mdf!(copy(A); finaldiscard)

end
