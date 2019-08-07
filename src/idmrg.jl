function runiDMRG()

    energy = Float64[]
    guessfidelity = Float64[]
    convergence = Float64[]
    truncation_errors = Float64[]
    Λs = Vector{Float64}[]

    ## 1. obtain the wavefunction for a two-site lattice : A_0 Λ_0 B_0
    ## and a four-site lattice : A_0 A_1 Λ_1 B_1 B_0
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
    Λ1 = fact.S
    B1 = fact.Vt

    push!(truncation_errors, 1-norm(Λ1))
    push!(Λs, Λ1)

    @tensor h[l, d1', d2', d3', d4', r, d1, d2, d3, d4] :=
        mpo.tensors[1][l, d1', m1, d1] * mpo.tensors[2][m1,d2',m2,d2] *
        mpo.tensors[3][m2,d3',m3,d3] * mpo.tensors[4][m3,d4',r,d4]
    eigfact = eigen(reshape(h,16,16))
    #println(eigfact.values)

    push!(energy, eigfact.values[1])
    push!(guessfidelity, 1)
    push!(convergence, 1)

    fact = svd(reshape(eigfact.vectors[:,1], 4,4), full=false)
    A2 = fact.U
    Λ2 = fact.S
    B2 = fact.Vt

    push!(truncation_errors, 1-norm(Λ2))
    push!(Λs, Λ2)




    ## 2. rotate the center to the left to get Λ^L_n B_{n+1}
    ## 3. rotate the center to the right to get A_{n+1} Λ^R_n
    ## 4. trial wavefunction for increased two-size is then:
    ## ... A_{n+1} Λ^R_n Λ_{n-1}^-1 Λ^L_n B_{n+1} ...
    ## 5. use the trial as the initial guess and eigensolve to get
    ## the wavefunction A_n+1 Λ_n+1 B_n+1
    ## 6. truncate to desired size
    ## 7. check if the the fixed point has been reached (tol)
    println("energy ", energy ./ [2*n-1 for n=1:2])
    println("guessfiledlity ", guessfidelity)
    println("convergence ", convergence)
    println("truncation_errors", truncation_errors)
    println("Lambdas : ")
    println(Λs)
end
