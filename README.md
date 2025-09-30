# MDF

`MDF.jl` is a Julia implementation of the minimum discard fill ordering for the incomplete LU factorization with zero level of fill in.

The [minimum discard fill](https://doi.org/10.1137/0613057) (MDF) permutation heuristically minimizes the discarded values of an (ILU) factorization. See tests as an example.

Installation
-------------
```
julia> ]
julia> add MDF
```

How to use
-------------
```
julia> using MDF, ILUZero
```
-`p = mdf(A)`: MDF permutation based on a symmetric sparse matrix A\
-`p = mdf!(A)`\
-`ilu0(A)`

Performance
-------------


Related
-------------
[`ILUZero.jl`](https://github.com/mcovalt/ILUZero.jl)

References
-------------
[D’Azevedo, E. F. and Forsyth, P. A. and Tang, Wei-Pai - Ordering Methods for Preconditioned Conjugate Gradient Methods Applied to Unstructured Grid Problems](https://doi.org/10.1137/0613057).