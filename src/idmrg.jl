##TODO: make this into a step function!
function idmrg2site(mpo::MatrixProductOperator{T},
                    maxdim::Int=100;
                    tol::Float64=1.e-8,
                    rounds::Int=20,
                    envelope::Vector{Float64}=ones(rounds),
                    verbose::Bool=false) where{T<:Number}

    @assert mpo.lx == 4

    energy = Float64[]
    guessfidelity = Float64[]
    convergence = Float64[]
    truncation_errors = Float64[]
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
    push!(dims, length(Λ))

    lmpo = reshape(mpo.tensors[1], d, size(mpo.tensors[1],3), d)
    @tensor envL[u,m,d] := conj(A1)[o',d] * lmpo[o', m, o] * A1[o,u]

    rmpo = reshape(mpo.tensors[4], size(mpo.tensors[4], 1), d, d)
    @tensor envR[u,m,d] := conj(B1)[d,o'] * rmpo[m, o', o] * B1[u,o]

    #envL = envelope[1]/envelope[2] .* envL
    #envR = envelope[1]/envelope[2] .* envR

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
    Λex = Λ
    Λ = fact.S
    matBn = fact.Vt

    push!(truncation_errors, 1-norm(Λ))
    Λ = Λ./norm(Λ)
    push!(dims, length(Λ))


    for n = 2:rounds-1
        A = reshape(matAn, dims[n-1], d, dims[n])
        @tensor envL[u,m,d] := ((envL[u',m',d'] * conj(A)[d', o', d]) *
                                mpo.tensors[2][m',o',m,o]) * A[u',o,u]

        B = reshape(matBn, dims[n], d, dims[n-1])
        @tensor envR[u,m,d] := ((envR[u',m',d'] * conj(B)[d, o', d']) *
                                mpo.tensors[3][m, o', m',o]) * B[u,o,u']

        # envL = envelope[n]/envelope[n+1] .* envL
        # envR = envelope[n]/envelope[n+1] .* envR

        ## 2. rotate the center to the left to get Λ^L_n B_{n+1}
        # U, S, Vt = svdtrunc(reshape(matAn * Diagonal(Λ), dims[n-1], d*dims[n]),
        #                     maxdim=maxdim, tol=tol)
        # Bnp1 = reshape(Vt, size(S, 1), d, dims[n])
        # Λln = U * S
        Q, R = qr(transpose(reshape(matAn * Diagonal(Λ), dims[n-1], d*dims[n])))
        Bnp1 = reshape(transpose(Matrix(Q)), dims[n-1], d, dims[n])
        Λln = transpose(R)

        ## 3. rotate the center to the right to get A_{n+1} Λ^R_n
        # U, S, Vt = svdtrunc(reshape(Diagonal(Λ) * matBn, dims[n]*d, dims[n-1]),
        #                     maxdim=maxdim, tol=tol)
        # Anp1 = reshape(U, dims[n], d, size(S,1))
        # Λrn = S * Vt
        Q, R = qr(reshape(Diagonal(Λ) * matBn, dims[n]*d, dims[n-1]))
        Anp1 = reshape(Matrix(Q), dims[n], d, dims[n-1])
        Λrn = R

        ## 4. trial wavefunction for increased two-size is then:
        ## ... A_{n+1} Λ^R_n Λ_{n-1}^-1 Λ^L_n B_{n+1} ...
        core = Λrn * Diagonal(1 ./Λex) * Λln
        @tensor guess[l,o1,o2,r] := Anp1[l,o1,mr] * core[mr,mm] * Bnp1[mm,o2,r]

        ## 5. use the trial as the initial guess and eigensolve to get
        ## the wavefunction A_n+1 Λ_n+1 B_n+1
        es, vs, info = eigsolve(v->_applymps2site(v, envL, envR,
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
        Λex = Λ
        Λ = diag(S)
        matBn = Vt

        push!(truncation_errors, 1-norm(Λ))
        Λ = Λ./norm(Λ)
        push!(dims, length(Λ))

        ## 7. check if the the fixed point has been reached (tol)
        # println("energy ", energy ./ [2*n-1 for n=1:length(energy)])
        # println("guessfiledlity ", guessfidelity)
        # println("convergence ", convergence)
        # println("truncation_errors", truncation_errors)
    end
    #println(dims)
    #make and imps
    n = rounds
    #M = reshape(matAn * Diagonal(Λ) * matBn, dims[n], d^2, dims[n])
    A = reshape(matAn, dims[n],d,dims[n])
    B = reshape(matBn, dims[n],d,dims[n])

return return A, Λ, B, Λex
end
