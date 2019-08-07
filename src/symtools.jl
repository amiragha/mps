function svdsym(smat::SymMatrix{Tv}) where {Tv<:Number}
    @assert signs(smat.legs) == (+1, -1)

    sects_U = Tuple{Int, Int}[]
    nzblks_U = Matrix{Tv}[]
    sects_S = Tuple{Int, Int}[]
    nzblks_S = Vector{Float64}[]
    nzblks_Vt = Matrix{Tv}[]
    sects_Vt = Tuple{Int, Int}[]
    middle_dims = Int[]

    for idx in eachindex(smat.sects)
        sect = smat.sects[idx]
        push!(sects_U, (sect[1], sect[1]))
        push!(sects_S, (sect[1], sect[2]))
        push!(sects_Vt, (sect[2], sect[2]))

        nzblk = smat.nzblks[idx]
        fact = svd(nzblk, full=false)
        push!(nzblks_U, fact.U)
        push!(nzblks_S, fact.S)
        push!(nzblks_Vt, fact.Vt)

        push!(middle_dims, length(fact.S))
    end

    Uleg2 = STLeg(-1, smat.legs[1].chrs, middle_dims)
    Sleg1 = STLeg(+1, smat.legs[1].chrs, middle_dims)
    Sleg2 = STLeg(-1, smat.legs[2].chrs, middle_dims)
    Vtleg1 = STLeg(+1, smat.legs[2].chrs, middle_dims)

    U = SymTensor(0, (smat.legs[1], Uleg2), sects_U, nzblks_U)
    S = SymTensor(smat.charge, (Sleg1, Sleg2), sects_S, [Diagonal(blk) for blk in nzblks_S])
    Vt = SymTensor(0, (Vtleg1, smat.legs[2]), sects_Vt, nzblks_Vt)
    return U, S, Vt
end
