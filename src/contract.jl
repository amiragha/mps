"""
    contract(A, idxA, B, idxB)

contract the symmetric tensor A, B. The to-be-contracted
indexes are shown with negative integers while the remaining
indexes for the result are shown with positive integers. For
example if there are l negative numbers in both A and B (should
be numbered from -1 to -l) then the rest of the indexes in A and B
combined should be 1:N+M-2n.
"""
function contract(A     ::AbstractSymTensor{T1, N},
                  idxA  ::NTuple{N, Int},
                  B     ::AbstractSymTensor{T2, M},
                  idxB  ::NTuple{M, Int};
                  debug ::Bool=false) where{T1<:Number, T2<:Number, N, M}

    remsA, consA, tofinalsA = _contract_index_perm(idxA)
    remsB, consB, tofinalsB = _contract_index_perm(idxB)

    # techniqually here the contraction vector space should be made
    # once! But now is made twice, so fix this!

    # TODO: check
    # compatibility println(remsA, consA, remsB, consB)
    if length(remsA) == 0
        matA = SymVector(A, -1, consA)
        if length(remsB) == 0
            matB = SymVector(B, +1, consB)
            return matA * matB
        else
            matB = SymMatrix(B, remsB, consB)
            return permutelegs(defuse_leg(matA*matB, 2, B.legs[remsB]),
                               invperm(tofinalsB))
        end
    end

    matA = SymMatrix(A, remsA, consA)

    if length(remsB) == 0
        matB = SymVector(B, +1, consB)
        return permutelegs(defuse_leg(matA*matB, 1, A.legs[remsA]),
                           invperm(tofinalsA))

    end

    matB = SymMatrix(B, consB, remsB)
    return permutelegs(
        defuse_leg(defuse_leg(matA * matB, 1, A.legs[remsA]), 2, B.legs[remsB]),
        invperm([tofinalsA;tofinalsB]))
end


function SymMatrix(A::AbstractSymTensor{T, N},
                   rowidxs::Vector{Int},
                   colidxs::Vector{Int}) where {T<:Number, N}

    idxperm = [rowidxs; colidxs]
    @assert sort(idxperm) == collect(1:N)
    pA = permutelegs(A, idxperm)

    length(rowidxs) == 0 && error("Zero row for matrix not allowed!")
    length(colidxs) == 0 && error("Zero col for matrix not allowed!")

    fusedsectors = NTuple{Int, Int}[]
    patranges = NTuple{UnitRange{Int}, UnitRange{Int}}[]

#     csects = sten.sects
#     fchrdict = Dict{Int, FusedCharge{n}}()
#     signs = Tuple(sten.legs[c].sign for c in l:l+n-1)
#     for i in eachindex(csects)
#         sect = [csects[i]...]

#         pat = Tuple(sect[l:l+n-1])
#         spat = sign .* signs .* pat
#         fchr = sign * sum(signs .* [sect[c] for c in l:l+n-1])
#         patrange = UnitRange{Int}(1,1)
#         if haskey(fchrdict, fchr)
#             if haskey(fchrdict[fchr].pats, spat)
#                 patrange = fchrdict[fchr].pats[spat]
#             else
#                 fdim = prod([getdim(sten.legs[c+l-1], pat[c]) for c in 1:n])
#                 dim = fchrdict[fchr].dim
#                 fchrdict[fchr].dim += fdim
#                 patrange = dim+1:dim+fdim
#                 fchrdict[fchr].pats[spat] = patrange
#             end
#         else
#             fdim = prod([getdim(sten.legs[c+l-1], pat[c]) for c in 1:n])
#             patrange = 1:fdim
#             fchrdict[fchr] = FusedCharge(fchr, fdim, Dict(spat => patrange))
#         end

#         fsect = vcat( sect[1:l-1], fchr, sect[l+n:end])
#         push!(fusedsectors, Tuple(fsect))
#         push!(patranges, patrange)
#     end

#     fusedsectperm = _sectors_sortperm(fusedsectors)
#     new_sectors, refs = uniquesorted(fusedsectors[fusedsectperm])
#     refs = refs[invperm(fusedsectperm)]

#     # the infos are: new_sectors, refs, patranges

#     # make the new leg
#     fleg = STLeg(sign, unzip([(k, fchrdict[k].dim) for k in sort(collect(keys(fchrdict)))])...)

#     new_legs = Tuple(vcat([sten.legs[1:l-1]...], fleg, [sten.legs[l+n:end]...]))

#     new_nzblks = Array{T, N-n+1}[]
#     for sect in new_sectors
#         push!(new_nzblks,Array{T, N-n+1}(undef, [getdim(new_legs[c], sect[c]) for c in 1:N-n+1]...))
#     end
#     #indexing = [1:end for i=1:N-n+1]...
#     for i in eachindex(sten.nzblks)
#         s = size(sten.nzblks[i])
#         new_nzblks[refs[i]][[1:s[i] for i=1:l-1]..., patranges[i],[1:s[i] for i=l+n:N]...] =
#             reshape(sten.nzblks[i], s[1:l-1]..., prod(s[l:l+n-1]), s[l+n:end]...)
#     end

#     SymTensor(sten.charge, new_legs, new_sectors, new_nzblks)


    fuselegs(
            fuselegs(psten, -1, length(rowidxs)+1, length(colidxs)),
            +1, 1, length(rowidxs))

end

# if debug
#         println(A)
#         println(remsA, " ", consA)
#         println(matA)
#         println(B)
#         println(remsB, " ", consB)
#         println(matB)
#         println(matA*matB)
#     end

# function contract(A::SymTensor{ComplexF64, N}, idxA::NTuple{N, Int},
#                   B::SymTensor{Float64, M}, idxB::NTuple{M, Int}) where{N, M}
#     contract(A, idxA, convert(SymTensor{ComplexF64, M}, B), idxB)
# end
# contract A and B using idxA and idxB. Note idxA and idxB are the
# corresponding leg numbers in each tensor that are to be contracted

###TODO: probably give this function a better name!
###TODO: replace this function by a version of fuselegs (that fuses both of them at once!)
function symMatrix(A::AbstractSymTensor{Tv, N},
                   rowidxs::Vector{Int},
                   colidxs::Vector{Int}) where{Tv<:Number, N}

    idxperm = [rowidxs; colidxs]
    @assert sort(idxperm) == collect(1:N)
    psten = permutelegs(sten, idxperm)
    if length(rowidxs) == 0
        return fuselegs(psten, -1, 1, N)
    elseif length(colidxs) == 0
        return fuselegs(psten, +1, 1, N)
    else
        return fuselegs(
            fuselegs(psten, -1, length(rowidxs)+1, length(colidxs)),
            +1, 1, length(rowidxs))
    end
end

function *(sten1::SymVector{Tv}, sten2::SymVector{Tv}) where {Tv<:Number}
    #bool, idx1, idx2 = _are_contractible(sten1.legs[1], sten2.legs[1])
    bool = _are_contractible(sten1.legs[1], sten2.legs[1])
    bool ||
        error("not contractible!", sten1.legs[2], " and ", sten2.legs[1])
    # @assert length(idx1) == length(idx2) == 1
    @assert signs(sten1.legs) == (-1,)
    @assert signs(sten2.legs) == (+1,)
    return sum(sten1.nzblks[1] .* sten2.nzblks[1])
end

function *(sten1::SymTensor{Tv, 2}, sten2::SymTensor{Tv, 2}) where{Tv<:Number}
    #bool, idx1, idx2 = _are_contractible(sten1.legs[2], sten2.legs[1])
    bool = _are_contractible(sten1.legs[2], sten2.legs[1])
    bool ||
        error("not contractible! ", sten1.legs[2], " and ", sten2.legs[1])
    @assert length(sten1.nzblks) == length(sten1.nzblks)
    sects = Tuple{Int, Int}[]
    nzblks = Matrix{Tv}[]

    @assert signs(sten1.legs) == (+1, -1)
    @assert signs(sten2.legs) == (+1, -1)
    # perm1 = _sectors_sortperm(sten1.sects, by=x->(x[2],))
    # perm2 = _sectors_sortperm(sten2.sects, by=x->(-x[1],))

    #println(idx1)
    #println(idx2)
    ##TODO: make the below better and explain!
    i=1
    j=1
    while i <= length(sten1.sects) && j <= length(sten2.sects)
        if sten1.sects[i][2] == sten2.sects[j][1]
            push!(sects, (sten1.sects[i][1],sten2.sects[j][2]))
            push!(nzblks, sten1.nzblks[i]*sten2.nzblks[j])
            i+=1
            j+=1
        elseif sten1.sects[i][2] > sten2.sects[i][1]
            j+=1
        else
            i+=1
        end
    end

    perm = _sectors_sortperm(sects)
    matC = SymTensor(sten1.charge+sten2.charge, (sten1.legs[1], sten2.legs[2]),
                     sects[perm], nzblks[perm])
    return matC
end

##TODO: make these conversion more sane!!!
function *(sten1::SymTensor{ComplexF64, 2},
           sten2::SymTensor{Float64, 2})
    *(sten1, convert(SymTensor{ComplexF64, 2}, sten2))
end

function *(sten1::SymTensor{Float64, 2},
           sten2::SymTensor{ComplexF64, 2})
    *(convert(SymTensor{ComplexF64, 2}, sten1), sten2)
end

function _are_contractible(l1::STLeg, l2::STLeg)
    if l1.sign == -l2.sign
        return l1.chrs == l2.chrs && l1.dims == l2.dims
    end
    false
end

# function _are_contractible(l1::STLeg, l2::STLeg)
#     if l1.sign == -l2.sign
#         # find conciding charges
#         chrs, idx1, idx2 = intersect(l1, l2)
#         for c in chrs
#             getdim(l1, c) != getdim(l2, c) && return false, idx1, idx2
#         end
#         return true, idx1, idx2
#     end
#     false, Int[], Int[]
# end

# This function returns the
function _contract_index_perm(indexset; mode::Symbol=:MINUS)
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
