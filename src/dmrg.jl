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

function initialenv(mps::MatrixProductState{ComplexF64},
                    mpo::MatrixProductOperator{Float64})
    initialenv(mps, convert(MatrixProductOperator{ComplexF64}, mpo))
end
"""
    dmrg1sitesweep!(mps, mpo, env, maxdim [; verbose])

Performs a single sweep (left to right and back) of 1site DMRG on
`mps` where the hamiltonian is given in `mpo` and environment matrices
for each step are stored in `env`.

Returns the energy and updates, `mps` and `env`
"""
function dmrg1sitesweep!(mps::MatrixProductState{T},
                         mpo::MatrixProductOperator{T},
                         env::Vector{Array{T, 3}};
                         verbose::Bool=false) where {T<:Number}
    lx = mps.lx
    d = mps.d

    if mps.center != 1
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

        @tensor mat[l,o,r] := Λ[l,r] * mps.matrices[l+1][m,o,r]
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
        Λ = transpose(Λ)

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
    dmrg2sitesweep!(mps, mpo, env, maxdim [; verbose])

Performs a single sweep (left to right and back) of 2site DMRG on
`mps` where the hamiltonian is given in `mpo` and environment matrices
for each step are stored in `env`. The singular values are truncated
upto truncation error or up to number of `maxdim`.

Returns the tuple of (energies, truncations) and updates, `mps` and `env`
"""
function dmrg2sitesweep!(mps::MatrixProductState{T},
                         mpo::MatrixProductOperator{T},
                         env::Vector{Array{T, 3}},
                         maxdim::Int=200,
                         tol::Float64=1.e-8;
                         verbose::Bool=false) where {T<:Number}
    lx = mps.lx
    d = mps.d

    if mps.center != 1
        move_center!(mps, 1)
    end

    mat = mps.matrices[1]
    for l = 1:lx-2
        @tensor vmat[-1,-2,-3,-4] := mat[-1,-2,1] * mps.matrices[l+1][1,-3,-4]
        es, vs, info = eigsolve(v->_dmrg2sitematvec(v, env[l], env[l+3],
                                                    mpo.tensors[l], mpo.tensors[l+1]),
                                vmat, 1, :SR, ishermitian=true)
        v = vs[1]
        e = es[1]
        verbose && println("Sweep L2R: site $l -> energy $e")
        #push!(energies, e)

        U, S, Vt = svdtrunc(reshape(v, size(v,1)*d, :), maxdim=maxdim, tol=tol)

        mps.matrices[l] = reshape(U, size(mat,1), d, :)
        mps.dims[l+1] = size(S, 1)

        env[l+1] = _dmrgupdateleft(env[l], mps.matrices[l], mpo.tensors[l])

        mat = reshape(S*Vt, size(S,1), d, size(mps.matrices[l+1],3))
    end

    l = lx-1
    @tensor vmat[-1,-2,-3,-4] := mat[-1,-2,1] * mps.matrices[l+1][1,-3,-4]
    es, vs, info = eigsolve(v->_dmrg2sitematvec(v, env[l], env[l+3],
                                                mpo.tensors[l], mpo.tensors[l+1]),
                            vmat, 1, :SR, ishermitian=true)
    v = vs[1]
    e = es[1]
    verbose && println("Sweep L2R: site $l -> energy $e")
    #push!(energies, e)

    for l = lx-1:-1:2
        U, S, Vt = svdtrunc(reshape(v, size(v,1)*d, :), maxdim=maxdim, tol=tol)

        mps.matrices[l+1] = reshape(Vt, size(Vt, 1), d, size(v,4))
        mps.dims[l+1] = size(S, 1)

        env[l+2] = _dmrgupdateright(env[l+3], mps.matrices[l+1], mpo.tensors[l+1])

        mat = reshape(U*S, size(v, 1), d, size(S, 2))
        @tensor vmat[-1,-2,-3,-4] := mps.matrices[l-1][-1,-2,1] * mat[1,-3,-4]

        es, vs, info = eigsolve(v->_dmrg2sitematvec(v, env[l-1], env[l+2],
                                                    mpo.tensors[l-1], mpo.tensors[l]),
                                vmat, 1, :SR, ishermitian=true)
        v = vs[1]
        e = es[1]
        verbose && println("Sweep L2R: site $l -> energy $e")
        #push!(energies, e)
    end
    l = 1
    U, S, Vt = svdtrunc(reshape(v, size(v,1)*d, :), maxdim=maxdim, tol=tol)

    mps.matrices[l+1] = reshape(Vt, size(Vt, 1), d, size(v,4))
    mps.dims[l+1] = size(S, 1)
    env[l+2] = _dmrgupdateright(env[l+3], mps.matrices[l+1], mpo.tensors[l+1])

    mps.matrices[l] = reshape(U*S, size(v, 1), d, size(S, 2))

    return e
end
