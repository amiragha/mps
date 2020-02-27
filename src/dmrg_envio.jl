"""
Generates the environment and saves them on infividual `.dat` files
with `prefix` and numbered from 1 to lx+2
"""
function _initialenv_tempfiles(mps::MPState{Ys},
                               mpo::MPOperator{Yo};
                               at::Int=1,
                               prefix="$(Date(now())), temporary_env") where {Ys, Yo}
    vtype(Ys) == vtype(Yo) && eltype(Ys) == eltype(Yo) ||
        error("MPS, MPO not the same type!")
    center(mps) == at || error("MPS center has to be $at")

    lx = length(mps)
    env = ["$(prefix)#$n.dat" for n in 1:lx+2]

    lVa = leftspace(mps)
    rVa = rightspace(mps)
    lVw = leftspace(mpo)
    rVw = rightspace(mpo)

    T = eltype(Ys)

    envL = fill(one(T), (dual(lVa), dual(lVw), lVa))
    @save env[1] envL

    envR = fill(one(T), (rVa, rVw, dual(rVa)))
    @save env[lx+2] envR

    for l = lx:-1:at
        envR = _mpsupdateright(envR, mps.As[l], mpo.Ws[l])
        @save env[l+1] envR
    end

    for l = 2:at-1
        envL = _mpsupdateleft(envL, mps.As[l], mpo.Ws[l])
        @save env[l+1] envL
    end
    return env
end

function dmrg2sitesweep_envio!(mps::MPState{Ys},
                               mpo::MPOperator{Yo},
                               env::Vector{String};
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
    envL = load(env[1], "envL")
    env_next = load(env[4], "envR")
    for l = 1:lx-2
        save(env[l], "envL", envL)
        envR = env_next
        env_next = load(env[l+4], "envR")

        AA = contract(A, (1,2,-1), mps.As[l+1], (-1,3,4))
        es, vs, info = eigsolve(v->_applymps2site(v, envL, envR,
                                                  mpo.Ws[l], mpo.Ws[l+1]),
                                AA, 1, :SR, ishermitian=true,
                                tol=lanczostol,
                                krylovdim=krylovdim,
                                maxiter=krylovmaxiter)
        v = vs[1]

        u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)
        truncerror = 1-norm(s)
        normalize!(s)
        AA = splitleg(splitleg(u*s*v, 2, space(AA)[3:4]), 1, space(AA)[1:2])
        e = dot(AA, _applymps2site(AA, envL, envR,
                                   mpo.Ws[l], mpo.Ws[l+1]))
        push!(energies, e)
        if verbose
            println("Sweep L2R: bond $l, $(l+1) -> energy $(es[1])")
            println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)")
            println("truncation error = $truncerror")
            println("energy after truncation $e")
            values = diag(s)
            println("Entanglement S1 = $(entropy(values.^2))")
            println()
        end
        mps.As[l] = splitleg(u, 1, space(A)[1:2])

        envL = _mpsupdateleft(envL, mps.As[l], mpo.Ws[l])

        A = splitleg(s*v, 2, space(AA)[3:4])
    end
    # # saving env[lx-1] which is 1...(lx-2) (no need to save!)
    # save(env[l], "envL", envL)
    envR = env_next
    env_next = envL
    for l = lx-1:-1:1
        save(env[l+3], "envR", envR)
        envL = env_next
        if l > 1
            env_next = load(env[l-1], "envL")
        end
        if l == lx-1
            AA = contract(A, (1,2,-1), mps.As[l+1], (-1,3,4))
        else
            AA = contract(mps.As[l], (1,2,-1), A, (-1,3,4))
        end
        es, vs, info = eigsolve(v->_applymps2site(v, envL, envR,
                                                  mpo.Ws[l], mpo.Ws[l+1]),
                                AA, 1, :SR, ishermitian=true,
                                tol=lanczostol,
                                krylovdim=krylovdim,
                                maxiter=krylovmaxiter)

        v = vs[1]
        u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)
        truncerror = 1-norm(s)
        normalize!(s)
        AA = splitleg(splitleg(u*s*v, 2, space(AA)[3:4]), 1, space(AA)[1:2])

        e = dot(AA, _applymps2site(AA, envL, envR,
                                   mpo.Ws[l], mpo.Ws[l+1]))
        push!(energies, e)

        if verbose
            println("Sweep L2R: bond $l, $(l+1) -> energy $(es[1])")
            println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)")
            println("truncation error = $truncerror")
            println("energy after truncation $e")
            values = diag(s)
            println("Entanglement S1 = $(entropy(values.^2))")
            println()
        end

        mps.As[l+1] = splitleg(v, 2, space(AA)[3:4])

        if l > 1
            envR = _mpsupdateright(envR, mps.As[l+1], mpo.Ws[l+1])
            A = splitleg(u*s, 1, space(AA)[1:2])
        else
            mps.As[1] = splitleg(u*s, 1, space(AA)[1:2])
        end
    end

    return energies
end

function dmrg2sitesweep_envasyncio!(mps::MPState{Ys},
                                    mpo::MPOperator{Yo},
                                    env::Vector{String};
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
    envL = load(env[1], "envL")
    env_next = load(env[4], "envR")
    for l = 1:lx-2
        @sync begin
            @async save(env[l], "envL", envL)
            envR = env_next
            @async env_next = load(env[l+4], "envR")

            AA = contract(A, (1,2,-1), mps.As[l+1], (-1,3,4))
            es, vs, info = eigsolve(v->_applymps2site(v, envL, envR,
                                                      mpo.Ws[l], mpo.Ws[l+1]),
                                    AA, 1, :SR, ishermitian=true,
                                    tol=lanczostol,
                                    krylovdim=krylovdim,
                                    maxiter=krylovmaxiter)
            v = vs[1]

            u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)
            truncerror = 1-norm(s)
            normalize!(s)
            AA = splitleg(splitleg(u*s*v, 2, space(AA)[3:4]), 1, space(AA)[1:2])
            e = dot(AA, _applymps2site(AA, envL, envR,
                                       mpo.Ws[l], mpo.Ws[l+1]))
            push!(energies, e)

            if verbose
                println("Sweep L2R: bond $l, $(l+1) -> energy $(es[1])")
                println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)")
                println("truncation error = $truncerror")
                println("energy after truncation $e")
                values = diag(s)
                println("Entanglement S1 = $(entropy(values.^2))")
                println()
            end
            mps.As[l] = splitleg(u, 1, space(A)[1:2])
            A = splitleg(s*v, 2, space(AA)[3:4])
        end
        envL = _mpsupdateleft(envL, mps.As[l], mpo.Ws[l])
    end
    # # saving env[lx-1] which is 1...(lx-2) (no need to save!)
    # save(env[l], "envL", envL)
    envR = env_next
    env_next = envL
    for l = lx-1:-1:1
        @sync begin
            @async save(env[l+3], "envR", envR)
            envL = env_next
            @async begin
                if l > 1
                    env_next = load(env[l-1], "envL")
                end
            end

            if l == lx-1
                AA = contract(A, (1,2,-1), mps.As[l+1], (-1,3,4))
            else
                AA = contract(mps.As[l], (1,2,-1), A, (-1,3,4))
            end
            es, vs, info = eigsolve(v->_applymps2site(v, envL, envR,
                                                      mpo.Ws[l], mpo.Ws[l+1]),
                                    AA, 1, :SR, ishermitian=true,
                                    tol=lanczostol,
                                    krylovdim=krylovdim,
                                    maxiter=krylovmaxiter)

            v = vs[1]
            u,s,v = svdtrunc(SymMatrix(v, [1,2], [3,4]), maxdim=maxdim, tol=tol)
            truncerror = 1-norm(s)
            normalize!(s)
            AA = splitleg(splitleg(u*s*v, 2, space(AA)[3:4]), 1, space(AA)[1:2])

            e = dot(AA, _applymps2site(AA, envL, envR,
                                       mpo.Ws[l], mpo.Ws[l+1]))
            push!(energies, e)

            if verbose
                println("Sweep L2R: bond $l, $(l+1) -> energy $(es[1])")
                println("normres = $(info.normres[1]), #iterations = $(info.numiter), #applications = $(info.numops)")
                println("truncation error = $truncerror")
                println("energy after truncation $e")
                values = diag(s)
                println("Entanglement S1 = $(entropy(values.^2))")
                println()
            end
            mps.As[l+1] = splitleg(v, 2, space(AA)[3:4])
            A = splitleg(u*s, 1, space(AA)[1:2])
        end
        if l > 1
            envR = _mpsupdateright(envR, mps.As[l+1], mpo.Ws[l+1])
        end
    end
    mps.As[1] = A
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
# """
#     dmrg1sitesweep!(mps, mpo, env, maxdim [; verbose])

# Performs a single sweep (left to right and back) of 1site DMRG on
# `mps` where the hamiltonian is given in `mpo` and environment matrices
# for each step are stored in `env`.

# Returns the energy and updates, `mps` and `env`
# """

# function dmrg1sitesweep!(mps::MatrixProductState{T},
#                          mpo::MatrixProductOperator{T},
#                          env::Vector{Array{T, 3}};
#                          verbose::Bool=false) where {T<:Number}
#     lx = mps.lx
#     d = mps.d

#     mps.center == 1 || error("The center of MPS has to be 1 for dmrg L->R->L")

#     mat = mps.matrices[1]
#     for l = 1:lx-1
#         es, vs, info = eigsolve(v->_applymps1site(v, env[l], env[l+2], mpo.tensors[l]),
#                                 mat, 1, :SR, ishermitian=true)
#         v = vs[1]
#         e = es[1]
#         verbose && println("Sweep L2R: site $l -> energy $e")
#         #push!(energies, e)

#         Q, Λ = qr(reshape(v, size(mat,1)*d, :))

#         ##NOTE: the Matrix function here is needed to return the thin Q!
#         mps.matrices[l] = reshape(Matrix(Q), size(mat))
#         env[l+1] = _mpsupdateleft(env[l], mps.matrices[l], mpo.tensors[l])

#         @tensor mat[l,o,r] := Λ[l,r] * mps.matrices[l+1][m,o,r]
#     end

#     l = lx
#     es, vs, info = eigsolve(v->_applymps1site(v, env[l], env[l+2], mpo.tensors[l]),
#                             mat, 1, :SR, ishermitian=true)
#     e=es[1]
#     verbose && println("Sweep L2R: site $l -> energy $e")
#     #push!(energies, e)

#     for l = lx-1:-1:1
#         Q, Λ = qr(transpose(reshape(mat, size(mat,1), :)))

#         ##NOTE: the matrix function here is needed to return the thin Q!
#         Q = transpose(Matrix(Q))
#         Λ = transpose(Λ)

#         mps.matrices[l+1] = reshape(Q, size(mat))
#         env[l+2] = _mpsupdateright(env[l+3], mps.matrices[l+1], mpo.tensors[l+1])
#         @tensor mat[-1,-2,-3] := mps.matrices[l][-1,-2, 1] * Λ[1,-3]

#         es, vs, info = eigsolve(v->_applymps1site(v, env[l], env[l+2], mpo.tensors[l]),
#                                 mat, 1, :SR, ishermitian=true)

#         v = vs[1]
#         e = es[1]
#         verbose && println("Sweep L2R: site $l -> energy $e")
#         #push!(energies, e)
#     end
#     mps.matrices[1] = mat

#     return e
# end
