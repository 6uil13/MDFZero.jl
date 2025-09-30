# MDFZero.jl
The [minimum discard fill](https://doi.org/10.1137/0613057) (MDF) permutation heuristically minimizes the discarded values of an incomplete LU (iLU) factorization.

`MDFZero.jl` is a Julia implementation of the minimum discard fill ordering for the iLU factorization with zero level of fill in.

Installation
-------------
```
julia> ]
julia> add MDFZero
```

Purpose
-------------
- When solving large sparse symmetric systems, iterative methods like PCG (Preconditioned Conjugate Gradients) are sensitive to the ordering of the degrees of freedom.
- This is specially important for systems involved in anisotropic problems, where matrix entries may contain large ratios.
- As a more accurate preconditioner, 
the MDF ordering aims to minimize the discarded fill between incomplete LU factorization and the original matrix. 

How to use
-------------
```
julia> using MDFZero, LinearAlgebra, SparseArrays
julia> A = Symmetric(sprand(1000, 1000, 1e-4) + 10I)
julia> p = mdf0(A)
julia> Fp = ILU0(A, p)
julia> F = ILU0(A)
julia> distance(A, Fp) < distance(A, F)
```
- `mdf0(A)`: MDF permutation based on a symmetric sparse matrix A, creates a copy of A
- `mdf0!(A)`: MDF permutation, updates matrix A
- `ILU0` uses `ILUZero.ilu0` function
- `distance`: measures the Frobenius norm of the difference between two matrices

Performance
-------------
```
julia> using MDFZero, LinearAlgebra, SparseArrays, BenchmarkTools
julia> A = Symmetric(sprand(1000, 1000, 1e-4) + 10I);
julia> p = @btime mdf0(A);
  6.240 ms (13049 allocations: 663.89 KiB)
```

Related
-------------
[`ILUZero.jl`](https://github.com/mcovalt/ILUZero.jl)

References
-------------
[D’Azevedo, E. F. and Forsyth, P. A. and Tang, Wei-Pai - Ordering Methods for Preconditioned Conjugate Gradient Methods Applied to Unstructured Grid Problems](https://doi.org/10.1137/0613057).