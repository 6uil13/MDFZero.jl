using MDFZero, Test, Random, SparseArrays, LinearAlgebra

function symsprand(n, ε = 1.0)
    R = sprand(n, n, ε)
    A = R + R' + 10*I

    Symmetric(A)
end

function laplacian_test()
    A = laplacian{2}(4)
    δ = zeros(size(A, 1))
    σ = mdf0(A, δ)
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
    q = mdf0(A)
    @test q == 1:16
end

@testset "MDFZero.jl" begin
    laplacian_test()
    
    n = 10
    A = symsprand(n, 0.5)
    p = mdf0(A)
    F = ILU0(A)
    Fp = ILU0(A, p)
    @test distance(A, Fp) <= distance(A, F)

    n = 100
    A = symsprand(n, 0.1)
    p = mdf0(A)
    F = ILU0(A)
    Fp = ILU0(A, p)
    @test distance(A, Fp) < distance(A, F)
end
