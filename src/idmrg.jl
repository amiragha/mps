function runiDMRG()

    energy = Float64[]
    guessfidelity = Float64[]
    convergence = Float64[]
    truncation_errors = Float64[]
    Λs = Vector{Float64}[]
    dims = Int[]

    ## 1. obtain the wavefunction for a two-site lattice : A_0 Λ_0 B_0

    physd = 2
    mpo = xxz_mpo(Float64, 4, 2)
    @tensor h[l, d1', d2', r, d1, d2] :=
        mpo.tensors[1][l, d1', m, d1] * mpo.tensors[4][m,d2',r,d2]
    eigfact = eigen(reshape(h, 4, 4))
    #println(eigfact.values)
    push!(energy, eigfact.values[1])
    push!(guessfidelity, 1)
    push!(convergence, 1)

    fact = svd(reshape(eigfact.vectors[:,1], 2,2), full=false)
    A1 = fact.U
    Λ = fact.S
    B1 = fact.Vt

    push!(truncation_errors, 1-norm(Λ))
    push!(Λs, Λ)
    push!(dims, length(Λ))

    lmpo = reshape(mpo.tensors[1], 2, 4, 2)
    @tensor left[dl,ml,ul] := conj(A1)[a,dl] * lmpo[a, ml, b] * A1[b,ul]
    rmpo = reshape(mpo.tensors[4], 4, 2, 2)
    @tensor rght[dr,mr,ur] := conj(B1)[dr,a] * rmpo[mr, a, b] * B1[ur,b]

    ## and now exactly find the eigenvalue of the left-mpo-right for
    ## the four-site lattice : A_0 A_1 Λ_1 B_1 B_0
    @tensor h4[dl,d1',d2',dr, ul,d1',d2',ur] :=
        (left[dl,ml,ul] * mpo.tensors[2][ml,d1',mm,d2]) *
        (mpo.tensors[3][mm, d2',mr, d2] * rght[dr, mr, ur])

    h4mat = reshape(h4, 2*physd*physd*2, 2*physd*physd*2)
    eigfact = eigen(h4)

    push!(energy, eigfact.values[1])
    push!(guessfidelity, 1)
    push!(convergence, 1)

    # n = 2
    fact = svd(reshape(eigfact.vectors[:,1], 4,4), full=false)
    matAn = fact.U
    Λ = fact.S
    matBn = fact.Vt

    push!(truncation_errors, 1-norm(Λ))
    push!(Λs, Λ)
    push!(dims, length(Λ))

    ten1 = reshape(matAn, 2, 2, 4)
    @tensor left[dl,ml,ul] = ((left[dl',ml',rl'] * conj(ten1)[dl', d', dl]) *
                              mpo.tensors[2][ml',d',ml,d]) * ten1[ul',d,ul]
    ten2 = reshape(matBn, 4, 2, 2)
    @tensor rght[dr,mr,ur] = ((rght[dr',mr',ur'] * conj(ten2)[dr, d', dr']) *
                              mpo.tensors[3][mr, d', mr',d]) * ten2[ur,d,ur']

    for n=2:10
        ## 2. rotate the center to the left to get Λ^L_n B_{n+1}
        fact = svd(reshape(matAn * Λs[n], dim[n-1], physd*dim[n]), full=false)
        Bnp1 = fact.Vt
        Λln = fact.U * fact.S

        ## 3. rotate the center to the right to get A_{n+1} Λ^R_n
        fact = svd(reshape(Λs[n] * matBn, dim[n]*physd, dim[n-1]), full=false)
        Anp1 = fact.U
        Λrn = fact.S * fact.Vt

        ## 4. trial wavefunction for increased two-size is then:
        ## ... A_{n+1} Λ^R_n Λ_{n-1}^-1 Λ^L_n B_{n+1} ...
        guess = Anp1 * Λrn * Diagonal(1./Λs[n-1]) * Λln * Bnp1
        guessten = reshape(guess, ldim, physdim, physdim, rdim)

        ## 5. use the trial as the initial guess and eigensolve to get
        ## the wavefunction A_n+1 Λ_n+1 B_n+1

        ## 6. truncate to desired size
        ## 7. check if the the fixed point has been reached (tol)
        println("energy ", energy ./ [2*n-1 for n=1:length(energy)])
        println("guessfiledlity ", guessfidelity)
        println("convergence ", convergence)
        println("truncation_errors", truncation_errors)
        println("Lambdas : ")
        println(Λs)
    end
end
