# contract A and B using idxA and idxB. Note idxA and idxB are the
# corresponding leg numbers in each tensor that are to be contracted

###TODO: probably give this function a better name!
function symMatrix(sten::SymTensor{Tv, N},
                   rowidxs::Vector{Int}, colidxs::Vector{Int}) where{Tv<:Number, N}
    idxperm = [rowidxs; colidxs]
    @assert sort(idxperm) == collect(1:N)
    psten = permutelegs(sten, idxperm)
    if length(rowidxs) == 0
        return fuse_conseqlegs(psten, -1, 1, N)
    elseif length(colidxs) == 0
        return fuse_conseqlegs(psten, +1, 1, N)
    else
        return fuse_conseqlegs(
            fuse_conseqlegs(psten, -1, length(rowidxs)+1, length(colidxs)),
            +1, 1, length(rowidxs))
    end
end

function *(sten1::SymTensor{Tv, 1}, sten2::SymTensor{Tv, 1}) where{Tv<:Number}
    _are_contractible(sten1.legs[1], sten2.legs[1]) || error("not contractible!")
    @assert length(sten1.nzblks) == length(sten1.nzblks) == 1
    @assert signs(sten1.legs) == (-1,)
    @assert signs(sten2.legs) == (+1,)
    return sum(sten1.nzblks[1] .* sten2.nzblks[1])
end

function *(sten1::SymTensor{Tv, 2}, sten2::SymTensor{Tv, 2}) where{Tv<:Number}
    _are_contractible(sten1.legs[2], sten2.legs[1]) || error("not contractible!")
    @assert length(sten1.nzblks) == length(sten1.nzblks)
    sects = Tuple{Int, Int}[]
    nzblks = Matrix{Tv}[]

    @assert signs(sten1.legs) == (+1, -1)
    @assert signs(sten2.legs) == (+1, -1)
    # perm1 = _sectors_sortperm(sten1.sects, by=x->(x[2],))
    # perm2 = _sectors_sortperm(sten2.sects, by=x->(-x[1],))

    for i = 1:length(sten1.sects)
        push!(sects, (sten1.sects[i][1],sten2.sects[i][2]))
        push!(nzblks, sten1.nzblks[i]*sten2.nzblks[i])
    end

    perm = _sectors_sortperm(sects)
    matC = SymTensor(sten1.charge+sten2.charge, (sten1.legs[1], sten2.legs[2]),
                     sects[perm], nzblks[perm])
    return matC
end

function _are_contractible(l1::STLeg, l2::STLeg)
    if l1.sign == -l2.sign
        return l1.chrs == l2.chrs && l1.dims == l2.dims
    end
    false
end


function _contract_index_perm(indexset)
    rems = Int[]
    cons = Int[]
    cons_order = Int[]
    tofinals = Int[]
    for i in eachindex(indexset)
        idx = indexset[i]
        if idx < 0
            push!(cons, i)
            push!(cons_order, idx)
        else
            push!(rems, i)
            push!(tofinals, idx)
        end
    end
    perm = sortperm(tofinals)
    rems[perm], cons[sortperm(cons_order, by=x->abs(x))], tofinals[perm]
end

"""
    contract(A, idxA, B, idxB)

contract the symmetric tensor A, B. The to-be-contracted
indexes are shown with negative integers while the remaining
indexes for the result are shown with positive integers. For
example if there are l negative numbers in both A and B (should
be numbered from -1 to -l) then the rest of the indexes in A and B
combined should be 1:N+M-2n.
"""
function contract(A::SymTensor{T, N}, idxA::NTuple{N, Int},
                  B::SymTensor{T, M}, idxB::NTuple{M, Int}) where{T<:Number, N, M}
    remsA, consA, tofinalsA = _contract_index_perm(idxA)
    remsB, consB, tofinalsB = _contract_index_perm(idxB)

    # TODO: check compatibility
    #println(remsA, consA, remsB, consB)
    matA = symMatrix(A, remsA, consA)
    matB = symMatrix(B, consB, remsB)

    if numoflegs(matA) == 1
        @assert numoflegs(matB) == 1
        return matA * matB
    end
    permutelegs(defuse_leg(defuse_leg(matA * matB, 2, B.legs[remsB]), 1, A.legs[remsA]),
                inv_perm([tofinalsA;tofinalsB]))
end
