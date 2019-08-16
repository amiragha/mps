### A bunch of helper functions to apply things to mps
function _applymps0site(Λ::Matrix{T},
                        envL::Array{T, 3},
                        envR::Array{T, 3}) where {T<:Number}
    @tensor A[d,d'] := envL[u,m,d] * Λ[u, u'] * envR[u',m,d']
    A
end

function _dmrgupdateright(renv::Array{T, 3},
                          mat::Array{T, 3},
                          hmpo::Array{T,4}) where {T<:Number}
    @tensor R[u,m,d] :=
        ((renv[u',m',d'] * mat[u,o,u']) * hmpo[m,o',m',o]) * conj(mat)[d,o',d']
    R
end

function _dmrgupdateleft(lenv::Array{T, 3},
                         mat::Array{T, 3},
                         hmpo::Array{T,4}) where {T<:Number}
    @tensor L[u,m,d] :=
        ((lenv[u',m',d'] * mat[u',o,u]) * hmpo[m',o',m,o]) * conj(mat)[d',o',d]
    L
end

function _dmrg1sitematvec(v,
                          envL::Array{T,3},
                          envR::Array{T,3},
                          hmpo::Array{T,4}) where {T<:Number}

    @tensor v[l,o,r] := ((envL[l',ml,l] * v[l',o',r']) *
                         hmpo[ml,o,mr,o']) * envR[r',mr,r]
    v
end
