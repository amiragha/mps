function initialenv(mps::SymMatrixProductState{Tv},
                    mpo::SymMatrixProductOperator{Tv};
                    isometry::Symbol=:R) where {Tv<:Number}
    lx = mps.lx
    env = Vector{SymTensor{Tv, 3}}(undef, lx+2)

    lchr = mps.matrices[1].legs[1].chrs[1]
    rchr = mps.matrices[lx].legs[3].chrs[1]

    ##NOTE: assuming mpo has charge 0 at both ends!
    # mpo.matrices[1].legs[1].chrs[1]
    # mpo.matrices[lx].legs[3].chrs[1]

    env[1] = fill(one(Tv), 0, (STLeg(-1,[lchr],[1]),
                               STLeg(-1,[0],[1]),
                               STLeg(+1,[lchr],[1])))
    env[lx+2] = fill(one(Tv), 0, (STLeg(+1,[rchr],[1]),
                                  STLeg(+1,[0],[1]),
                                  STLeg(-1,[rchr],[1])))

    if isometry == :R
        for l = lx:-1:1
            mps.center == 1 || error("MPS center has to be 1")
            env[l+1] = _mpsupdateright(env[l+2], mps.matrices[l], mpo.tensors[l])
        end
    elseif isometry == :L
        for l = 1:lx
            mps.center == lx || error("MPS center has to be $lx")
            env[l+1] = _mpsupdateleft(env[l], mps.matrices[l], mpo.tensors[l])
        end
    else
        error("Invalid isometry value $isometry")
    end
    return env
end

function initialenv(mps::SymMatrixProductState{ComplexF64},
                    mpo::SymMatrixProductOperator{Float64};
                    isometry::Symbol=:R)
    initialenv(mps,
               convert(MatrixProductOperator{ComplexF64}, mpo),
               isometry=isometry)
end


function dmrg2sitesweep!(mps::SymMatrixProductState{Tv},
                         mpo::SymMatrixProductOperator{Tv},
                         env::Vector{SymTensor{Tv, 3}};
                         maxdim::Int=200,
                         tol::Float64=1.e-9,
                         lanczostol::Float64=1.e-7,
                         verbose::Bool=false) where {Tv<:Number}
    lx = mps.lx
    d = mps.d

    mps.center == 1 || error("The center of MPS has to be 1 for dmrg L->R->L")

    energies = Vector{Float64}()

    A = mps.matrices[1]
    for l = 1:lx-2
        AA = contract(A, (1,2,-1), mps.matrices[l+1], (-1,3,4))
        es, vs, info = eigsolve(v->_applymps2site(v, env[l], env[l+3],
                                                  mpo.tensors[l], mpo.tensors[l+1]),
                                AA, 1, :SR, ishermitian=true,
                                tol=lanczostol)#, krylovdim=5, maxiter=8)

        v = vs[1]
        e = es[1]
        verbose && println("Sweep L2R: bond $l, $(l+1) -> energy $e")
        verbose && println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)\n")
        push!(energies, e)

        vreleg = fuselegs(fuselegs(v, +1, 1, 2), -1, 2, 2)
        #vreleg = SymMatrix(v, [1,2], [3,4])
        U, S, Vt = svdtrunc(vreleg, maxdim=maxdim, tol=tol)

        mps.matrices[l] = unfuseleg(U, 1, A.legs[1:2])
        mps.dims[l+1] = size(S, 1)

        env[l+1] = _mpsupdateleft(env[l], mps.matrices[l], mpo.tensors[l])

        A = unfuseleg(S*Vt, 2, AA.legs[3:4])
    end

    l = lx-1
    AA = contract(A, (1,2,-1), mps.matrices[l+1], (-1,3,4))
    es, vs, info = eigsolve(v->_applymps2site(v, env[l], env[l+3],
                                              mpo.tensors[l], mpo.tensors[l+1]),
                            AA, 1, :SR, ishermitian=true,
                            tol=lanczostol)#, krylovdim=5, maxiter=8)
    v = vs[1]
    e = es[1]
    verbose && println("Sweep L2R: bond $l, $(l+1) -> energy $e")
    verbose && println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)\n")
    push!(energies, e)

    for l = lx-1:-1:2
        vreleg = fuselegs(fuselegs(v, +1, 1, 2), -1, 2, 2)
        #vreleg = SymMatrix(v, [1,2], [3,4])
        U, S, Vt = svdtrunc(vreleg, maxdim=maxdim, tol=tol)

        mps.matrices[l+1] = unfuseleg(Vt, 2, AA.legs[3:4])
        mps.dims[l+1] = size(S, 2)

        env[l+2] = _mpsupdateright(env[l+3], mps.matrices[l+1], mpo.tensors[l+1])

        A = unfuseleg(U*S, 1, AA.legs[1:2])
        AA = contract(mps.matrices[l-1], (1,2,-1), A, (-1,3,4))

        es, vs, info = eigsolve(v->_applymps2site(v, env[l-1], env[l+2],
                                                  mpo.tensors[l-1], mpo.tensors[l]),
                                AA, 1, :SR, ishermitian=true,
                                tol=lanczostol)#, krylovdim=5, maxiter=8)
        v = vs[1]
        e = es[1]
        verbose && println("Sweep R2L: bond $(l-1), $l -> energy $e")
        verbose && println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)\n")
        push!(energies, e)
    end
    l = 1
    vreleg = fuselegs(fuselegs(v, +1, 1, 2), -1, 2, 2)
    #vreleg = SymMatrix(v, [1,2], [3,4])
    U, S, Vt = svdtrunc(vreleg, maxdim=maxdim, tol=tol)

    mps.matrices[l+1] = unfuseleg(Vt, 2, AA.legs[3:4])
    mps.dims[l+1] = size(S, 2)

    env[l+2] = _mpsupdateright(env[l+3], mps.matrices[l+1], mpo.tensors[l+1])

    mps.matrices[l] = unfuseleg(U*S, 1, AA.legs[1:2])
    return energies
end

# function dmrg2site!(mps::SymMatrixProductState{Tv},
#                     mpo::SymMatrixProductOperator{Tv};
#                     infostring::String="",
#                     finalmeasurements::Vector{Measurement}=[],
#                     sweepmeasurements::Vector{Measurement}=[],
#                     n_sweeps::Int=4,
#                     maxdim::Int=200,
#                     tol::Float64=1.e-8,
#                     verbose::Bool=false)

#     n < 1 && return

#     env = initialenv(mps, mpo)
#     energies = zeros(Float64, 2*mps.lx-3, n_sweeps)
#     for n in 1:n_sweeps
#         verbose && println("Starting sweeps ")
#         energies[:, n] =
#             dmrg2sitesweep!(mps, mpo, env, maxdim=maxdim, tol=tol, verbose=verbose)
#         if length(sweepmeasurements) > 0
#             # do the measurements
#         end
#     end

# end
