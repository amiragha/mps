"""
    correlationmatrix(hopmatrix, n_occupied)

calculates the two-body correlation matrix lambda ``Λ = ⟨a^†_i a_j⟩``
by full diagonalization of given Hamiltonian matrix.

Here, the correlation matrix Λ is made from its diagonal form Γ, which
has n_occupied 1s--corresponding to the lowest eigenvalues--and the
rest of the diagonal are 0s. We have ``Λ = U^{*} Γ U^{T}``

"""
function correlationmatrix(hopmatrix::Matrix{Float64},
                           n_occupied ::Int64)

    lx = size(hopmatrix)[1]
    @assert 0 < n_occupied && n_occupied <= lx

    fact = eigen(Symmetric(hopmatrix))
    println(fact.values[n_occupied-2:n_occupied+2])
    e_fermi = fact.values[n_occupied]
    n=1
    while fact.values[n_occupied+n] - e_fermi < 1.e-10
        n=n+1
    end
    n > 1 && @warn "Hamiltonian is $n-fold degenerate at the given filling factor $n_occupied/$lx"
    vecs = fact.vectors[:,1:n_occupied]
    #vecs = eigvecs(Symmetric(hopmatrix))[:, 1:n_occupied]
    return vecs * transpose(vecs)
end

function correlationmatrix(hopmatrix::Matrix{ComplexF64},
                           n_occupied ::Int64)

    lx = size(hopmatrix)[1]
    @assert 0 < n_occupied && n_occupied <= lx

    fact = eigen(Hermitian(hopmatrix))
    e_fermi = fact.values[n_occupied]
    n=1
    while fact.values[n_occupied+n] - e_fermi < 1.e-10
        n=n+1
    end
    n > 1 && @warn "Hamiltonian is $n-fold degenerate at the given filling factor $n_occupied/$lx"
    vecs = fact.vectors[:,1:n_occupied]
    #vecs = eigvecs(Hermitian(hopmatrix))[:, 1:n_occupied]
    return conj(vecs) * transpose(vecs)
end
