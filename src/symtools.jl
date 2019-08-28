function svdsym(smat::SymMatrix{Tv}; debug::Bool=false) where {Tv<:Number}
    signs(smat.legs) != (+1, -1) &&
        error("svdsym only accespts a SymMatrix (+1,-1) but ", signs(smat.legs))

    debug && println(smat)

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

    ls, rs = zip(sects_S...)
    lchrs = [ls...]
    rchrs = [rs...]

    if debug
        println(middle_dims)
        println(lchrs)
        println(rchrs)
    end

    Uleg2 = STLeg(-1, lchrs, middle_dims)
    Sleg1 = STLeg(+1, lchrs, middle_dims)
    Sleg2 = STLeg(-1, rchrs, middle_dims)
    Vtleg1 = STLeg(+1, rchrs, middle_dims)

    U = SymTensor(0, (smat.legs[1], Uleg2), sects_U, nzblks_U)
    S = SymTensor(smat.charge, (Sleg1, Sleg2), sects_S, [Diagonal(blk) for blk in nzblks_S])
    Vt = SymTensor(0, (Vtleg1, smat.legs[2]), sects_Vt, nzblks_Vt)
    return U, S, Vt
end


"""
    fermionswapgate(l1, l2)

Fot the two given legs, `l1` and `l2` make the fermionic swap gate
which is `-1` if the two legs has charge odd and `+1` otherwise.

The output legs are l2_afterX, l1_afterX, l2, l1
"""
function fermionswapgate(l1::STLeg, l2::STLeg)

    (l1.chrs != [0, 1] || l1.dims != [1,1]) && error("fermionswap oops!")

    l1X = STLeg(-l1.sign, l1.chrs, l1.dims)
    l2X = STLeg(-l2.sign, l2.chrs, l2.dims)
    legs = (l2, l1, l2X, l1X)

    sects = NTuple{4, Int}[]
    nzblks = Array{Float64, 4}[]

    for i1 in eachindex(l1.chrs)
        for i2 in eachindex(l2.chrs)
            c1, c2 = l1.chrs[i1], l2.chrs[i2]
            d1, d2 = l1.dims[i1], l2.dims[i2]
            push!(sects, (c2, c1, c2, c1))
            s = isodd(c1) && isodd(c2) ? -1 : +1
            push!(nzblks,s * reshape(I(d2), d2, d1, d2, d1))
        end
    end
    SymTensor(0, legs, sects, nzblks)
end
