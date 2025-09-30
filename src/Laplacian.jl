struct laplacian{d}
end

⊗ = kron

function laplacian{d}(n::Int) where d
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