using MDF
using Test

include("ILU0.jl")
using Random
import Base.permute!, SparseArrays.permute

permute!(A::Symmetric, p, q) = permute!(A.data, p, q)
permute(A::Symmetric, p, q) = permute(A.data, p, q)
ILU(A::Symmetric) = ILU(A.data)
ILU(A::Symmetric, p, q) = ILU(A.data, p, q)

function symsprand(n, ε = 1.0)
    R = sprand(n, n, ε)
    A = R + R' + 10*I

    Symmetric(A)
end

function errorilu(A, F)
    L, U, p, q = F

    norm(permute(A, p, q) - L * U)
end

struct laplacian{d}
end

⊗ = kron

function laplacian{d}(n::Int64) where d
    if d == 1
        e = ones(n)
        A = spdiagm(0 => 2 * e, 1 => -e[1:n-1], -1 => -e[1:n-1])
    else
        A₁ = laplacian{d - 1}(n)
        𝓘 = one(A₁)
        A = A₁ ⊗ 𝓘 + 𝓘 ⊗ A₁  
    end

    Symmetric(A)
end

function papertest()
    A = laplacian{2}(4)
    σ, δ = MDF.mdf(A; finaldiscard = true)
    p = [1, 4, 13, 16, 2, 3, 5, 9, 8, 12, 14, 15, 6, 7, 10, 11]
    corner = √0.125
    sidefirst = √32 / 15
    interiorfirst = δ[p[13]]
    d = zeros(16)
    d[1:4] .= corner
    d[5:2:12] .= sidefirst
    d[13] = interiorfirst # unknown

    @test σ == p
    @test δ[p] ≈ d

    permute!(A, p, p)
    q = MDF.mdf(A)
    @test q == 1:16
end

@testset "MDF.jl" begin
    papertest()
    
    n = 10
    A = symsprand(n, 0.5)
    p = MDF.mdf(A)
    F = ILU(A)
    Fp = ILU(A, p, p)
    @test errorilu(A, Fp) <= errorilu(A, F)

    n = 100
    A = symsprand(n, 0.1)
    p = MDF.mdf(A)
    F = ILU(A)
    Fp = ILU(A, p, p)
    @test errorilu(A, Fp) <= errorilu(A, F)
end
