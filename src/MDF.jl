module MDF

using SparseArrays, LinearAlgebra
export mdf, mdf!

function tiebreacking(discard::AbstractMatrix, candidates::Vector)
    n, l = size(discard)
    k = length(candidates)
    j = 1
    while j ≤ l && length(candidates) > 1
        discardj = view(discard, candidates, j)
        d = minimum(discardj)
        candidates = candidates[abs.(discardj .- d) .< eps(0.0)]
        j += 1
    end

    candidates[1]
end

function update!(A::SparseMatrixCSC{T, Int64}, m::Int64) where T <: Real
    a = A[:, m]
    amm = a[m]
    for r = A.colptr[m]:A.colptr[m + 1] - 1
        i = A.rowval[r]
        aim = a[i]

        𝓘 = A.colptr[i]:A.colptr[i + 1] - 1
        rows = view(A.rowval, 𝓘)
        for s = A.colptr[m]:A.colptr[m + 1] - 1
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

function discardedfill(A::SparseMatrixCSC{T, Int64}, m::Int64) where T <: Real
    amm = A[m, m]
    f = 0.0
    colnorm = 0.0
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
    colnorm += abs(amm)

    [f, defficiency, degree]
end

Base.@propagate_inbounds function mdf!(S::Symmetric{T, SparseMatrixCSC{T, Int64}};
              finaldiscard::Bool = false) where T <: Real

    A = S.data
    n = size(A, 1)

    # 4 measures: discard, defficiency, degree, label
    discard = zeros(n, 4)
    σ = zeros(Int64, n)
    𝓥 = collect(1:n)
    discard[:, 4] = 𝓥

    # initial discard
    for m in 𝓥
        discard[m, 1:3] = discardedfill(A, m)
    end

    # main decomposition loop
    for k = 1:n-1
        m = tiebreacking(discard, 𝓥) # most expensive
        σ[k] = m
        setdiff!(𝓥, m)

        𝓝 = A.rowval[A.colptr[m]:A.colptr[m + 1] - 1]
        setdiff!(𝓝, σ)

    	update!(A, m)
        
        Base.Threads.@threads for v in 𝓝
            discard[v, 1:3] = discardedfill(A, v)
        end
    end
    σ[n] = 𝓥[1]
    
    if finaldiscard
        return σ, .√discard[:, 1]
    else
        return σ
    end
end

mdf(A; finaldiscard::Bool = false) = mdf!(copy(A); finaldiscard)

end
