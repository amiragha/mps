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
    v#reshape(v, prod(size(v)))
end

"""
    dmrg1sitesweep!(mps, mpo, env, maxdim [; verbose])

Performs a single sweep (left to right and back) of 1site DMRG on
`mps` where the hamiltonian is given in `mpo` and environment matrices
for each step are stored in `env`. The singular values are truncated
to `maxdim`.

Returns the tuple of (energies, truncations) and updates, `mps` and `env`
"""
function dmrg1sitesweep!(mps::MatrixProductState{T},
                         mpo::MatrixProductOperator{T},
                         env::Vector{Array{T, 3}};
                         verbose::Bool=false) where {T<:Number}
    lx = mps.lx
    d = mps.d

    if mps.center == 1
        move_center!(mps, 1)
    end

    mat = mps.matrices[1]
    for l = 1:lx-1
        es, vs, info = eigsolve(v->_dmrg1sitematvec(v, env[l], env[l+2], mpo.tensors[l]),
                                mat, 1, :SR, ishermitian=true)
        v = vs[1]
        e = es[1]
        verbose && println("Sweep L2R: site $l -> energy $e")
        #push!(energies, e)

        Q, Λ = qr(reshape(v, size(mat,1)*d, :))

        ##NOTE: the Matrix function here is needed to return the thin Q!
        mps.matrices[l] = reshape(Matrix(Q), size(mat))
        env[l+1] = _dmrgupdateleft(env[l], mps.matrices[l], mpo.tensors[l])

        @tensor mat[-1,-2,-3] := Λ[-1,1] * mps.matrices[l+1][1,-2,-3]
    end
    l = lx
    es, vs, info = eigsolve(v->_dmrg1sitematvec(v, env[l], env[l+2], mpo.tensors[l]),
                            mat, 1, :SR, ishermitian=true)
    e=es[1]
    verbose && println("Sweep L2R: site $l -> energy $e")
    #push!(energies, e)

    for l = lx-1:-1:1
        Q, Λ = qr(transpose(reshape(mat, size(mat,1), :)))

        ##NOTE: the matrix function here is needed to return the thin Q!
        Q = transpose(Matrix(Q))

        mps.matrices[l+1] = reshape(Q, size(mat))
        env[l+2] = _dmrgupdateright(env[l+3], mps.matrices[l+1], mpo.tensors[l+1])
        @tensor mat[-1,-2,-3] := mps.matrices[l][-1,-2, 1] * Λ[1,-3]

        es, vs, info = eigsolve(v->_dmrg1sitematvec(v, env[l], env[l+2], mpo.tensors[l]),
                                mat, 1, :SR, ishermitian=true)

        v = vs[1]
        e = es[1]
        verbose && println("Sweep L2R: site $l -> energy $e")
        #push!(energies, e)
    end
    mps.matrices[1] = mat

    return e
end

"""
    initialenv(mps, mpo)

initialize and environment for a dmrg algorithm that is about to start
from left!

"""
function initialenv(mps::MatrixProductState{T},
                    mpo::MatrixProductOperator{T}) where {T<:Number}
    lx = mps.lx
    env = Vector{Array{T, 3}}(undef, lx+2)
    env[1] = ones(T, 1,1,1)
    env[lx+2] = ones(T, 1,1,1)
    for l = lx:-1:1
        env[l+1] = _dmrgupdateright(env[l+2], mps.matrices[l], mpo.tensors[l])
    end
    return env
end
