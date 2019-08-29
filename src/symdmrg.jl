function initialenv(mps::SymMatrixProductState{Tv},
                    mpo::SymMatrixProductOperator{Tv}) where {Tv<:Number}
    lx = mps.lx
    env = Vector{SymTensor{Tv, 3}}(undef, lx+2)

    lchr = mps.matrices[1].legs[1].chrs[1]
    rchr = mps.matrices[lx].legs[3].chrs[1]

    ##NOTE: assuming mpo has charge 0 at both ends!
    # mpo.matrices[1].legs[1].chrs[1]
    # mpo.matrices[lx].legs[3].chrs[1]

    env[1] = fillSymTensor(one(Tv), 0, (STLeg(-1,[lchr],[1]),
                                        STLeg(-1,[0],[1]),
                                        STLeg(+1,[lchr],[1])))
    env[lx+2] = fillSymTensor(one(Tv), 0, (STLeg(+1,[rchr],[1]),
                                           STLeg(+1,[0],[1]),
                                           STLeg(-1,[rchr],[1])))
    for l = lx:-1:1
        env[l+1] = _mpsupdateright(env[l+2], mps.matrices[l], mpo.tensors[l])
    end
    return env
end

function initialenv(mps::SymMatrixProductState{ComplexF64},
                    mpo::SymMatrixProductOperator{Float64})
    initialenv(mps, convert(MatrixProductOperator{ComplexF64}, mpo))
end


function dmrg2sitesweep!(mps::SymMatrixProductState{Tv},
                         mpo::SymMatrixProductOperator{Tv},
                         env::Vector{SymTensor{Tv, 3}},
                         maxdim::Int=200,
                         tol::Float64=1.e-8;
                         verbose::Bool=false) where {Tv<:Number}
    lx = mps.lx
    d = mps.d

    if mps.center != 1
        move_center!(mps, 1)
    end

    A = mps.matrices[1]
    for l = 1:lx-2
        AA = contract(A, (1,2,-1), mps.matrices[l+1], (-1,3,4))
        es, vs, info = eigsolve(v->_applymps2site(v, env[l], env[l+3],
                                                  mpo.tensors[l], mpo.tensors[l+1]),
                                AA, 1, :SR, ishermitian=true)
        v = vs[1]
        e = es[1]
        verbose && println("Sweep L2R: site $l -> energy $e")
        #push!(energies, e)

        vreleg = fuselegs(fuselegs(v, +1, 1, 2), -1, 2, 2)
        U, S, Vt = svdtrunc(vreleg, maxdim=maxdim, tol=tol)

        mps.matrices[l] = defuse_leg(U, 1, A.legs[1:2])
        mps.dims[l+1] = size(S, 1)

        env[l+1] = _mpsupdateleft(env[l], mps.matrices[l], mpo.tensors[l])

        A = defuse_leg(S*Vt, 2, AA.legs[3:4])
    end

    l = lx-1
    AA = contract(A, (1,2,-1), mps.matrices[l+1], (-1,3,4))
    es, vs, info = eigsolve(v->_applymps2site(v, env[l], env[l+3],
                                              mpo.tensors[l], mpo.tensors[l+1]),
                            AA, 1, :SR, ishermitian=true)
    v = vs[1]
    e = es[1]
    verbose && println("Sweep L2R: site $l -> energy $e")
    #push!(energies, e)

    for l = lx-1:-1:2
        vreleg = fuselegs(fuselegs(v, +1, 1, 2), -1, 2, 2)
        U, S, Vt = svdtrunc(vreleg, maxdim=maxdim, tol=tol)

        mps.matrices[l+1] = defuse_leg(Vt, 2, AA.legs[3:4])
        mps.dims[l+1] = size(S, 2)

        env[l+2] = _mpsupdateright(env[l+3], mps.matrices[l+1], mpo.tensors[l+1])

        A = defuse_leg(U*S, 1, AA.legs[1:2])
        AA = contract(mps.matrices[l-1], (1,2,-1), A, (-1,3,4))

        es, vs, info = eigsolve(v->_applymps2site(v, env[l-1], env[l+2],
                                                  mpo.tensors[l-1], mpo.tensors[l]),
                                AA, 1, :SR, ishermitian=true)
        v = vs[1]
        e = es[1]
        verbose && println("Sweep R2L: site $l -> energy $e")
        #push!(energies, e)
    end
    l = 1
    vreleg = fuselegs(fuselegs(v, +1, 1, 2), -1, 2, 2)
    U, S, Vt = svdtrunc(vreleg, maxdim=maxdim, tol=tol)

    mps.matrices[l+1] = defuse_leg(Vt, 2, AA.legs[3:4])
    mps.dims[l+1] = size(S, 2)

    env[l+2] = _mpsupdateright(env[l+3], mps.matrices[l+1], mpo.tensors[l+1])

    mps.matrices[l] = defuse_leg(U*S, 1, AA.legs[1:2])
    return e
end