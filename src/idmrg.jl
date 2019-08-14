function idmrg2site(mpo::MatrixProductOperator{T},
                    maxdim::Int=100;
                    tol::Float64=tol,
                    rounds::Int=20,
                    envelope::Vector{Float64}=ones(rounds),
                    verbose::Bool=false) where{T<:Number}

    @assert mpo.lx == 4

    energy = Float64[]
    guessfidelity = Float64[]
    convergence = Float64[]
    truncation_errors = Float64[]
    Λs = Vector{Float64}[]
    dims = Int[]

    d = size(mpo.tensors[2], 2)

    @tensor h[l, d1', d2', r, d1, d2] :=
        mpo.tensors[1][l, d1', m, d1] * mpo.tensors[4][m,d2',r,d2]
    eigfact = eigen(reshape(h, d*d, d*d))

    v = eigfact.vectors[:, 1]
    e = eigfact.values[1]
    verbose && println("IDMRG round 1 -> Energy = $e")

    push!(energy, e)
    push!(guessfidelity, 1)
    push!(convergence, 1)

    fact = svd(reshape(v, d,d), full=false)
    A1 = fact.U
    Λ = fact.S
    B1 = fact.Vt

    push!(truncation_errors, 1-norm(Λ))
    push!(Λs, Λ)
    push!(dims, length(Λ))

    lmpo = reshape(mpo.tensors[1], d, size(mpo.tensors[1],3), d)
    @tensor envL[u,m,d] := conj(A1)[o',d] * lmpo[o', m, o] * A1[o,u]

    rmpo = reshape(mpo.tensors[4], size(mpo.tensors[4], 1), d, d)
    @tensor envR[u,m,d] := conj(B1)[d,o'] * rmpo[m, o', o] * B1[u,o]

    ## and now exactly find the eigenvalue of the left-mpo-right for
    ## the four-site lattice : A_0 A_1 Λ_1 B_1 B_0
    @tensor h4[dl,o1',o2',dr, ul,o1,o2,ur] :=
        (envL[ul,ml,dl] * mpo.tensors[2][ml,o1',mm,o1]) *
        (mpo.tensors[3][mm, o2',mr, o2] * envR[ur, mr, dr])

    eigfact = eigen(reshape(h4, 2*d*d*2, 2*d*d*2))

    v = eigfact.vectors[:, 1]
    e = eigfact.values[1]/3
    verbose && println("IDMRG round 2 -> Energy = $e")

    push!(energy, e)
    push!(guessfidelity, 1)
    push!(convergence, 1)

    fact = svd(reshape(v, 4,4), full=false)
    matAn = fact.U
    Λ = fact.S
    matBn = fact.Vt

    push!(truncation_errors, 1-norm(Λ))
    push!(Λs, Λ)
    push!(dims, length(Λ))


    for n = 2:rounds

        ten1 = reshape(matAn, dims[n-1], d, dims[n])
        @tensor envL[u,m,d] := ((envL[u',m',d'] * conj(ten1)[d', o', d]) *
                                mpo.tensors[2][m',o',m,o]) * ten1[u',o,u]

        ten2 = reshape(matBn, dims[n], d, dims[n-1])
        @tensor envR[u,m,d] := ((envR[u',m',d'] * conj(ten2)[d, o', d']) *
                                mpo.tensors[3][m, o', m',o]) * ten2[u,o,u']

        ## 2. rotate the center to the left to get Λ^L_n B_{n+1}
        U, S, Vt = svdtrunc(reshape(matAn * Diagonal(Λs[n]), dims[n-1], d*dims[n]),
                            maxdim=maxdim, tol=tol)
        Bnp1 = reshape(Vt, size(S, 1), d, dims[n])
        Λln = U * S

        ## 3. rotate the center to the right to get A_{n+1} Λ^R_n
        U, S, Vt = svdtrunc(reshape(Diagonal(Λs[n]) * matBn, dims[n]*d, dims[n-1]),
                            maxdim=maxdim, tol=tol)
        Anp1 = reshape(U, dims[n], d, size(S,1))
        Λrn = S * Vt

        ## 4. trial wavefunction for increased two-size is then:
        ## ... A_{n+1} Λ^R_n Λ_{n-1}^-1 Λ^L_n B_{n+1} ...
        core = Λrn * Diagonal(1 ./Λs[n-1]) * Λln
        @tensor guess[l,o1,o2,r] := Anp1[l,o1,mr] * core[mr,mm] * Bnp1[mm,o2,r]

        ## 5. use the trial as the initial guess and eigensolve to get
        ## the wavefunction A_n+1 Λ_n+1 B_n+1
        es, vs, info = eigsolve(v->_dmrg2sitematvec(v, envL, envR,
                                                    mpo.tensors[2], mpo.tensors[3]),
                                guess, 1, :SR, ishermitian=true)

        v = vs[1]
        e = es[1]/(2*n-1)
        verbose && println("IDMRG round $(n+1) -> Energy = $e")

        push!(energy, e)
        push!(guessfidelity, 1)
        push!(convergence, 1)

        ## 6. truncate to desired size
        U, S, Vt = svdtrunc(reshape(v, dims[n]*d, d*dims[n]), maxdim=maxdim, tol=tol)

        matAn = U
        Λ = diag(S)
        matBn = Vt

        push!(truncation_errors, 1-norm(Λ))
        push!(Λs, Λ./norm(Λ))
        push!(dims, length(Λ))

        ## 7. check if the the fixed point has been reached (tol)
        # println("energy ", energy ./ [2*n-1 for n=1:length(energy)])
        # println("guessfiledlity ", guessfidelity)
        # println("convergence ", convergence)
        # println("truncation_errors", truncation_errors)
        # println("Lambdas : ")
        # println(Λs)
    end

    #make and imps
    n = rounds
    M = reshape(matAn * Diagonal(Λ) * matBn, dims[n], d^2, dims[n])
    imps = InfiniteMatrixProductState(M)

    return imps, energy, truncation_errors
end
