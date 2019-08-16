"""
    tdvp1sitesweep!(dt, mps, mpo, env [; verbose])

Performs a single sweep (left to right and back) of the 1site TDVP on
`mps` where the `mpo` is the Hamiltonian and environment matrices are
stored in `env`. This process symplectically evolves the `mps` in time
for the times step `dt`.

See arXiv:1408.5056 for details of the algorithm and explanations

"""
function tdvp1sitesweep!(dt::Float64,
                         mps::MatrixProductState{T},
                         mpo::MatrixProductOperator{T},
                         env::Vector{Array{T, 3}};
                         verbose::Bool=false) where {T<:Number}

    if !(T <: Complex)
        error("TDVP only accpets complex MPSs. Convert first!")
    end
    lx = mps.lx
    d = mps.d

    if mps.center != 1
        move_center!(mps, 1)
    end

    A = mps.matrices[1]
    for l = 1:lx-1
        # forward evolution of mps at site l
        A, info = exponentiate(v->_dmrg1sitematvec(v, env[l], env[l+2], mpo.tensors[l]),
                               -im*dt, A, ishermitian=true)

        if verbose
            e = dot(A, _dmrg1sitematvec(A, env[l], env[l+2], mpo.tensors[l]))
            println("Sweep L2R: mps site $l -> energy $e")
        end

        Q, Λ = qr(reshape(A, size(A,1)*size(A,2), size(A,3)))
        mps.matrices[l] = reshape(Matrix(Q), size(A))
        env[l+1] = _dmrgupdateleft(env[l], mps.matrices[l], mpo.tensors[l])

        # backward evolution of Λ
        Λ, info = exponentiate(v->_applymps0site(v, env[l+1], env[l+2]),
                               +im*dt, Λ)

        if verbose
            e = dot(Λ, _applymps0site(Λ, env[l+1], env[l+2]))
            println("Sweep L2R: Λ between  site $l and $(l+1) -> energy $e")
        end

        @tensor A[l,o,r] := Λ[l,m] * mps.matrices[l+1][m,o,r]
    end

    l = lx
    A, info = exponentiate(v->_dmrg1sitematvec(v, env[l], env[l+2], mpo.tensors[l]),
                           -im*dt, A, ishermitian=true)

    if verbose
        e = dot(A, _dmrg1sitematvec(A, env[l], env[l+2], mpo.tensors[l]))
        println("Sweep L2R: mps site $l -> energy $e")
    end

    for l = lx-1:-1:1
        Q, Λ = qr(transpose(reshape(A, size(A, 1), size(A,2)*size(A,3))))
        Λ = Matrix(transpose(Λ))

        mps.matrices[l+1] = reshape(transpose(Matrix(Q)), size(A))
        env[l+2] = _dmrgupdateright(env[l+3], mps.matrices[l+1], mpo.tensors[l+1])

        # backward evolution of Λ
        Λ, info = exponentiate(v->_applymps0site(v, env[l+1], env[l+2]),
                               +im*dt, Λ)
        if verbose
            e = dot(Λ, _applymps0site(Λ, env[l+1], env[l+2]))
            println("Sweep L2R: Λ between  site $l and $(l+1) -> energy $e")
        end

        @tensor A[l,o,r] := mps.matrices[l][l,o,m] * Λ[m,r]

        # forward evolution of mps at site l
        A, info = exponentiate(v->_dmrg1sitematvec(v, env[l], env[l+2], mpo.tensors[l]),
                               -im*dt, A; ishermitian = true)

        if verbose
            e = dot(A, _dmrg1sitematvec(A, env[l], env[l+2], mpo.tensors[l]))
            println("Sweep R2L: mps site $l -> energy $e")
        end
    end
    mps.matrices[1] = A
    return nothing
end
