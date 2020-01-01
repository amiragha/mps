function initialenv(mps::MPState{S, T},
                    mpo::MPOperator{S, T};
                    isometry::Symbol=:R) where {S, T}
    lx = length(mps)
    env = Vector{SymTensor{S,T,3}}(undef, lx+2)

    lchr = charges(mps.As[1].space[1])[1]
    rchr = charges(mps.As[lx].space[3])[1]

    ##NOTE: assuming mpo has charge 0 at both ends!
    lVa = VectorSpace{S}(lchr => 1)
    rVa = VectorSpace{S}(rchr => 1)
    Vw = VectorSpace{S}(0 => 1)

    env[1] = fill(one(T), zero(S), (dual(lVa), dual(Vw), lVa))
    env[lx+2] = fill(one(T), zero(S), (rVa, Vw, dual(rVa)))

    if isometry == :R
        for l = lx:-1:1
            mps.center == 1 || error("MPS center has to be 1")
            # println(env[l+2].space)
            # println(mps.As[l].space)
            # println(mpo.Ws[l].space)
            env[l+1] = _mpsupdateright(env[l+2], mps.As[l], mpo.Ws[l])
        end
    elseif isometry == :L
        for l = 1:lx
            mps.center == lx || error("MPS center has to be $lx")
            env[l+1] = _mpsupdateleft(env[l], mps.As[l], mpo.Ws[l])
        end
    else
        error("Invalid isometry value $isometry")
    end
    return env
end

function initialenv(mps::MPState{S,T1},
                    mpo::MPOperator{S,T2};
                    isometry::Symbol=:R) where {S,T1,T2}
    T = promote_type(T1, T2)
    initialenv(convert(MPState{S,T}, mps),
               convert(MPOperator{S,T}, mpo),
               isometry=isometry)
end


function dmrg2sitesweep!(mps::MPState{S, T},
                         mpo::MPOperator{S, T},
                         env::Vector{SymTensor{S,T,3}};
                         maxdim::Int=200,
                         tol::Float64=1.e-9,
                         lanczostol::Float64=1.e-7,
                         krylovdim::Int=5,
                         krylovmaxiter::Int=8,
                         verbose::Bool=false) where {S,T<:Number}

    lx = length(mps)
    mps.center == 1 || error("The center of MPS has to be 1 for dmrg L->R->L")

    energies = Vector{Float64}()
    A = mps.As[1]
    for l = 1:lx-2
        AA = contract(A, (1,2,-1), mps.As[l+1], (-1,3,4))
        es, vs, info = eigsolve(v->_applymps2site(v, env[l], env[l+3],
                                                  mpo.Ws[l], mpo.Ws[l+1]),
                                AA, 1, :SR, ishermitian=true,
                                tol=lanczostol,
                                krylovdim=krylovdim,
                                maxiter=krylovmaxiter)
        v = vs[1]
        e = es[1]
        verbose && println("Sweep L2R: bond $l, $(l+1) -> energy $e")
        verbose && println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)\n")
        push!(energies, e)
        u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)

        mps.As[l] = splitleg(u, 1, A.space[1:2])

        env[l+1] = _mpsupdateleft(env[l], mps.As[l], mpo.Ws[l])

        A = splitleg(s*v, 2, AA.space[3:4])
    end

    l = lx-1
    AA = contract(A, (1,2,-1), mps.As[l+1], (-1,3,4))
    es, vs, info = eigsolve(v->_applymps2site(v, env[l], env[l+3],
                                              mpo.Ws[l], mpo.Ws[l+1]),
                            AA, 1, :SR, ishermitian=true,
                            tol=lanczostol,
                            krylovdim=krylovdim,
                            maxiter=krylovmaxiter)

    v = vs[1]
    e = es[1]
    verbose && println("Sweep L2R: bond $l, $(l+1) -> energy $e")
    verbose && println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)\n")
    push!(energies, e)

    for l = lx-1:-1:2
        u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)

        mps.As[l+1] = splitleg(v, 2, AA.space[3:4])
        env[l+2] = _mpsupdateright(env[l+3], mps.As[l+1], mpo.Ws[l+1])

        A = splitleg(u*s, 1, AA.space[1:2])
        AA = contract(mps.As[l-1], (1,2,-1), A, (-1,3,4))

        es, vs, info = eigsolve(v->_applymps2site(v, env[l-1], env[l+2],
                                                  mpo.Ws[l-1], mpo.Ws[l]),
                                AA, 1, :SR, ishermitian=true,
                                tol=lanczostol,
                                krylovdim=krylovdim,
                                maxiter=krylovmaxiter)
        v = vs[1]
        e = es[1]
        verbose && println("Sweep R2L: bond $(l-1), $l -> energy $e")
        verbose && println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)\n")
        push!(energies, e)
    end
    l = 1
    u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)

    mps.As[l+1] = splitleg(v, 2, AA.space[3:4])
    env[l+2] = _mpsupdateright(env[l+3], mps.As[l+1], mpo.Ws[l+1])

    mps.As[l] = splitleg(u*s, 1, AA.space[1:2])
    return energies
end

# function dmrg2site!(mps::MPState{S,Tv},
#                     mpo::MPOperator{Tv};
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
