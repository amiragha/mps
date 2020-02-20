"""
    initialenv(mps, mpo)

construct and initialize the environment for a dmrg algorithm that is
about to start at positon `at` (defaulted to `1`)! The environment is
made from `mps` and `mpo` by the _mpsupdateright function.

"""
function initialenv(mps::MPState{Ys},
                    mpo::MPOperator{Yo};
                    at::Int=1) where {Ys, Yo}
    vtype(Ys) == vtype(Yo) && eltype(Ys) == eltype(Yo) ||
        error("MPS, MPO not the same type!")
    center(mps) == at || error("MPS center has to be $at")

    lx = length(mps)
    env = Vector{Ys}(undef, lx+2)

    lVa = leftspace(mps)
    rVa = rightspace(mps)
    lVw = leftspace(mpo)
    rVw = rightspace(mpo)

    T = eltype(Ys)
    env[1] = fill(one(T), (dual(lVa), dual(lVw), lVa))
    env[lx+2] = fill(one(T), (rVa, rVw, dual(rVa)))

    for l = lx:-1:at
        env[l+1] = _mpsupdateright(env[l+2], mps.As[l], mpo.Ws[l])
    end

    for l = 1:at-1
        env[l+1] = _mpsupdateleft(env[l], mps.As[l], mpo.Ws[l])
    end
    return env
end

# function _initialenv(mps::MPState{Y1},
#                      mpo::MPOperator{Y2};
#                      at::Int=1) where {Y1, Y2}
#     Y = promote_type(Y1, Y2)
#     initialenv(convert(MPState{Y}, mps),
#                convert(MPOperator{Y}, mpo),
#                at=at)
# end


function dmrg2sitesweep!(mps::MPState{Ys},
                         mpo::MPOperator{Yo},
                         env::Vector{Ys};
                         maxdim::Int=200,
                         tol::Float64=1.e-9,
                         lanczostol::Float64=1.e-7,
                         krylovdim::Int=5,
                         krylovmaxiter::Int=8,
                         verbose::Bool=false) where {Ys, Yo}

    vtype(Ys) == vtype(Yo) && eltype(Ys) == eltype(Yo) ||
        error("MPS, MPO not the same type!")
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
        #e = es[1]
        verbose && println("Sweep L2R: bond $l, $(l+1) -> energy $(es[1])")
        verbose && println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)")

        u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)
        verbose && println("truncation error = $(1-norm(s))")
        normalize!(s)
        AA = splitleg(splitleg(u*s*v, 2, space(AA)[3:4]), 1, space(AA)[1:2])
        e = dot(AA, _applymps2site(AA, env[l], env[l+3],
                                   mpo.Ws[l], mpo.Ws[l+1]))
        verbose && println("energy after truncation $e")
        push!(energies, e)
        if verbose
            values = diag(s)
            println("Entanglement S1 = $(entropy(values.^2))")
            println()
        end
        mps.As[l] = splitleg(u, 1, space(A)[1:2])

        env[l+1] = _mpsupdateleft(env[l], mps.As[l], mpo.Ws[l])

        A = splitleg(s*v, 2, space(AA)[3:4])
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
    #e = es[1]
    verbose && println("Sweep L2R: bond $l, $(l+1) -> energy $(es[1])")
    verbose && println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)")
    # verbose && println("truncation error = 1")
    # push!(energies, e)

    for l = lx-1:-1:2
        u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)
        verbose && println("truncation error = $(1-norm(s))")
        normalize!(s)
        AA = splitleg(splitleg(u*s*v, 2, space(AA)[3:4]), 1, space(AA)[1:2])
        e = dot(AA, _applymps2site(AA, env[l], env[l+3],
                                   mpo.Ws[l], mpo.Ws[l+1]))
        verbose && println("energy after truncation $e")
        push!(energies, e)

        if verbose
            values = diag(s)
            println("Entanglement S1 = $(entropy(values.^2))")
            println()
        end

        mps.As[l+1] = splitleg(v, 2, space(AA)[3:4])
        env[l+2] = _mpsupdateright(env[l+3], mps.As[l+1], mpo.Ws[l+1])

        A = splitleg(u*s, 1, space(AA)[1:2])
        AA = contract(mps.As[l-1], (1,2,-1), A, (-1,3,4))

        es, vs, info = eigsolve(v->_applymps2site(v, env[l-1], env[l+2],
                                                  mpo.Ws[l-1], mpo.Ws[l]),
                                AA, 1, :SR, ishermitian=true,
                                tol=lanczostol,
                                krylovdim=krylovdim,
                                maxiter=krylovmaxiter)
        v = vs[1]
        #e = es[1]
        verbose && println("Sweep R2L: bond $(l-1), $l -> energy $(es[1])")
        verbose && println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)")
        push!(energies, e)
    end
    l = 1
    u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)
    verbose && println("truncation error = $(1-norm(s))")
    normalize!(s)
    AA = splitleg(splitleg(u*s*v, 2, space(AA)[3:4]), 1, space(AA)[1:2])
    e = dot(AA, _applymps2site(AA, env[l], env[l+3],
                               mpo.Ws[l], mpo.Ws[l+1]))
    verbose && println("energy after truncation $e")
    push!(energies, e)

    if verbose
        values = diag(s)
        println("Entanglement S1 = $(entropy(values.^2))")
        println()
    end

    mps.As[l+1] = splitleg(v, 2, space(AA)[3:4])
env[l+2] = _mpsupdateright(env[l+3], mps.As[l+1], mpo.Ws[l+1])

mps.As[l] = splitleg(u*s, 1, space(AA)[1:2])
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
