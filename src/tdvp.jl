"""
    tdvp1sitesweep!(dt, mps, mpo, env [; verbose])

Performs a single sweep (left to right and back) of the 1site TDVP on
`mps` where the `mpo` is the Hamiltonian and environment matrices are
stored in `env`. This process symplectically evolves the `mps` in time
for the times step `dt`.

See arXiv:1408.5056 for details of the algorithm and explanations

"""
function tdvp1sitesweep!(dt::Float64,
                         mps::MPS{T},
                         mpo::MPO{T},
                         env::Vector{Array{T, 3}};
                         verbose::Bool=false) where {T<:Number}

    !(T <: Complex) && error("TDVP only accpets complex MPSs. Convert first!")

    lx = length(mps)

    center_at!(mps, 1)

    A = mps.As[1]
    for l = 1:lx-1
        # forward evolution of mps at site l
        A, info = exponentiate(v->_applymps1site(v, env[l], env[l+2], mpo.Ws[l]),
                               -im*dt, A, ishermitian=true)

        if verbose
            e = dot(A, _applymps1site(A, env[l], env[l+2], mpo.Ws[l]))
            println("Sweep L2R: mps site $l -> energy $e")
        end

        Q, Λ = qr(reshape(A, size(A,1)*size(A,2), size(A,3)))
        mps.As[l] = reshape(Matrix(Q), size(A))
        env[l+1] = _mpsupdateleft(env[l], mps.As[l], mpo.Ws[l])

        # backward evolution of Λ
        Λ, info = exponentiate(v->_applymps0site(v, env[l+1], env[l+2]),
                               +im*dt, Λ)

        if verbose
            e = dot(Λ, _applymps0site(Λ, env[l+1], env[l+2]))
            println("Sweep L2R: Λ between  site $l and $(l+1) -> energy $e")
        end

        @tensor A[l,o,r] := Λ[l,m] * mps.As[l+1][m,o,r]
    end

    l = lx
    A, info = exponentiate(v->_applymps1site(v, env[l], env[l+2], mpo.Ws[l]),
                           -im*dt, A, ishermitian=true)

    if verbose
        e = dot(A, _applymps1site(A, env[l], env[l+2], mpo.Ws[l]))
        println("Sweep L2R: mps site $l -> energy $e")
    end

    for l = lx-1:-1:1
        Q, Λ = qr(transpose(reshape(A, size(A, 1), size(A,2)*size(A,3))))
        Λ = Matrix(transpose(Λ))

        mps.As[l+1] = reshape(transpose(Matrix(Q)), size(A))
        env[l+2] = _mpsupdateright(env[l+3], mps.As[l+1], mpo.Ws[l+1])

        # backward evolution of Λ
        Λ, info = exponentiate(v->_applymps0site(v, env[l+1], env[l+2]),
                               +im*dt, Λ)
        if verbose
            e = dot(Λ, _applymps0site(Λ, env[l+1], env[l+2]))
            println("Sweep L2R: Λ between  site $l and $(l+1) -> energy $e")
        end

        @tensor A[l,o,r] := mps.As[l][l,o,m] * Λ[m,r]

        # forward evolution of mps at site l
        A, info = exponentiate(v->_applymps1site(v, env[l], env[l+2], mpo.Ws[l]),
                               -im*dt, A; ishermitian = true)

        if verbose
            e = dot(A, _applymps1site(A, env[l], env[l+2], mpo.Ws[l]))
            println("Sweep R2L: mps site $l -> energy $e")
        end
    end
    mps.As[1] = A
    mps.center = 1
    mps
end

"""
    tdvp2sitesweep!(dt, mps, mpo, env [; maxdim, tol, verbose])

Performs a single sweep (left to right and back) of the 2site TDVP on
`mps` where the `mpo` is the Hamiltonian and environment As are
stored in `env`. This process is not a result of direct formulation of
tangent spaces, because there are now a hierarchy of manifold of MPS
with different system sizes, but it is a direct suggestion from the
one-site version.

See arXiv:1408.5056 for details of the algorithm and explanations/
justification of the two-site proposed algorithm.

"""
function tdvp2sitesweep!(dt::Float64,
                         mps::MPState{Array{T, 3}},
                         mpo::MPOperator{Array{T, 4}},
                         env::Vector{Array{T, 3}};
                         maxdim::Int=200,
                         tol::Float64=1.e-14,
                         verbose::Bool=false) where {T<:Number}

    !(T <: Complex) && error("TDVP only accpets complex MPSs. Convert first!")
    lx = length(mps)

    center_at!(mps, 1)

    A = mps.As[1]
    for l = 1:lx-2
        # forward evolution of mps at site l and site l+1
        @tensor AA[-1,-2,-3,-4] := A[-1,-2,1] * mps.As[l+1][1,-3,-4]
        AA, info = exponentiate(v->_applymps2site(v, env[l], env[l+3],
                                                  mpo.Ws[l], mpo.Ws[l+1]),
                                -im*dt, AA, ishermitian=true)

        if verbose
            e = dot(AA, _applymps2site(AA, env[l], env[l+3],
                                       mpo.Ws[l], mpo.Ws[l+1]))
            println("Sweep L2R: mps site $l, $(l+1) -> energy $e")
        end

        u,s,v = svdtrunc(reshape(AA, prod(size(AA)[1:2]), prod(size(AA)[3:4])),
                         maxdim=maxdim, tol=tol)
        mps.As[l] = reshape(Matrix(u), (size(AA)[1:2]..., size(s,1)))
        env[l+1] = _mpsupdateleft(env[l], mps.As[l], mpo.Ws[l])

        # backward evolution of Λ at site l+1
        Λ = reshape(s*v, (size(s,1), size(AA)[3:4]...))
        Λ, info = exponentiate(v->_applymps1site(v, env[l+1], env[l+3], mpo.Ws[l+1]),
                               +im*dt, Λ)

        if verbose
            e = dot(Λ, _applymps1site(Λ, env[l+1], env[l+3], mpo.Ws[l+1]))
            println("Sweep L2R: Λ on site $(l+1) -> energy $e")
        end

        A = Λ
        #@tensor A[l,o,r] := Λ[l,m] * mps.As[l+1][m,o,r]
    end

    l = lx-1
    @tensor AA[-1,-2,-3,-4] := A[-1,-2,1] * mps.As[l+1][1,-3,-4]
    AA, info = exponentiate(v->_applymps2site(v, env[l], env[l+3],
                                              mpo.Ws[l], mpo.Ws[l+1]),
                            -im*dt, AA, ishermitian=true)

    if verbose
        e = dot(AA, _applymps2site(AA, env[l], env[l+3],
                                   mpo.Ws[l], mpo.Ws[l+1]))
        println("Sweep L2R: mps site $l, $(l+1) -> energy $e")
    end

    for l = lx-1:-1:2

        u, s, v = svdtrunc(reshape(AA, prod(size(AA)[1:2]), prod(size(AA)[3:4])),
                           maxdim=maxdim, tol=tol)
        mps.As[l+1] = reshape(v, (size(s,1),size(AA)[3:4]...))
        env[l+2] = _mpsupdateright(env[l+3], mps.As[l+1], mpo.Ws[l+1])

        # backward evolution of Λ at site l
        Λ = reshape(u*s, (size(AA)[1:2]...,size(s,2)))
        Λ, info = exponentiate(v->_applymps1site(v, env[l], env[l+2], mpo.Ws[l]),
                               +im*dt, Λ)
        if verbose
            e = dot(Λ, _applymps1site(Λ, env[l], env[l+2], mpo.Ws[l]))
            println("Sweep R2L: Λ on site $l -> energy $e")
        end

        A = Λ
        #@tensor A[l,o,r] := mps.As[l-1][l,o,m] * Λ[m,r]
        @tensor AA[-1,-2,-3,-4] := mps.As[l-1][-1,-2,1] * A[1,-3,-4]

        # forward evolution of mps at site l
        AA, info = exponentiate(v->_applymps2site(v, env[l-1], env[l+2],
                                                  mpo.Ws[l-1], mpo.Ws[l]),
                                -im*dt, AA; ishermitian = true)

        if verbose
            e = dot(AA, _applymps2site(AA, env[l-1], env[l+2],
                                       mpo.Ws[l-1], mpo.Ws[l]))
            println("Sweep R2L: mps site $(l-1), $(l) -> energy $e")
        end
    end
    l=1
    u, s, v = svdtrunc(reshape(AA, prod(size(AA)[1:2]), prod(size(AA)[3:4])),
                       maxdim=maxdim, tol=tol)
    mps.As[l+1] = reshape(v, (size(s,1),size(AA)[3:4]...))
    env[l+2] = _mpsupdateright(env[l+3], mps.As[l+1], mpo.Ws[l+1])


    Λ = reshape(u*s, (size(AA)[1:2]..., size(s,2)))
    ### (IMPORTANT) QUESTION? do I need a final backward evolution on site 1?
    # Λ, info = exponentiate(v->_applymps1site(v, env[l], env[l+2], mpo.Ws[l]),
    #                        +im*dt, Λ)

    mps.As[1] = Λ
    mps.center = 1

    mps
end
