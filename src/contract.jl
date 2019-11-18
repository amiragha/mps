"""
    contract(A, idxA, B, idxB)

contract the symmetric tensor A, B. The to-be-contracted
indexes are shown with negative integers while the remaining
indexes for the result are shown with positive integers. For
example if there are l negative numbers in both A and B (should
be numbered from -1 to -l) then the rest of the indexes in A and B
combined should be 1:N+M-2n.
"""
function contract(A     :: AbstractSymTensor{T1, N},
                  idxA  :: NTuple{N, Int},
                  B     :: AbstractSymTensor{T2, M},
                  idxB  :: NTuple{M, Int};
                  debug :: Bool=false) where{T1<:Number, T2<:Number, N, M}

    remsA, consA, tofinalsA = _contract_index_perm(idxA)
    remsB, consB, tofinalsB = _contract_index_perm(idxB)

    ## TODO: techniqually here the contraction vector space should be made
    ## once! But now is made twice, so find a way fix this!

    # TODO: check
    # compatibility println(remsA, consA, remsB, consB)
    if length(remsA) == 0
        _A = SymVector(A, -1, consA)
        if length(remsB) == 0
            _B = SymVector(B, +1, consB)
            return _A * _B
        else
            _B = SymMatrix(B, remsB, consB)
            return permutelegs(unfuseleg(_A*_B, 2, B.legs[remsB]),
                               invperm(tofinalsB))
        end
    end

    _A = SymMatrix(A, remsA, consA)

    if length(remsB) == 0
        _B = SymVector(B, +1, consB)
        return permutelegs(unfuseleg(_A*_B, 1, A.legs[remsA]),
                           invperm(tofinalsA))

    end

    _B = SymMatrix(B, consB, remsB)
    return permutelegs(
        unfuseleg(unfuseleg(_A * _B, 1, A.legs[remsA]), length(remsA)+1, B.legs[remsB]),
        invperm([tofinalsA;tofinalsB]))
end

# auxillary function to find indexes for contract
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

function SymVector(A::AbstractSymTensor,
                   sign::Int,
                   perm::Vector{Int})
    N = numoflegs(A)
    abs(sign) == 1 || error("Incorrect sign $sign")
    sort(perm) == collect(1:N) || error("Incorrect index set for conversion to SymVector!")

    pA = permutelegs(A, perm)
    SymVector(fuselegs(pA, sign, 1, N))
end

function SymMatrix(A::AbstractSymTensor,
                   rowidxs::Vector{Int},
                   colidxs::Vector{Int})

    N = numoflegs(A)
    idxperm = [rowidxs; colidxs]
    sort(idxperm) == collect(1:N) ||
        error("Incorrect index set for conversion to SymMatrix!")
    length(rowidxs) == 0 && error("Zero row for matrix not allowed!")
    length(colidxs) == 0 && error("Zero col for matrix not allowed!")

    T = eltype(A)

    n_row = length(rowidxs)
    n_col = length(colidxs)
    rowlegs = A.legs[rowidxs]
    collegs = A.legs[colidxs]
    rsigns = [leg.sign for leg in rowlegs]
    csigns = [-leg.sign for leg in collegs]

    perm1 = _sectors_sortperm(A.sects, by=x->x[idxperm])
    csects = A.sects[perm1]
    fsects = Vector{Tuple{Int, Int}}(undef, length(A.sects))
    for i in eachindex(A.sects)
        sect = csects[i]
        c1 = sum([rsigns[i] * sect[rowidxs][i] for i=1:n_row])
        c2 = sum([csigns[i] * sect[colidxs][i] for i=1:n_col])
        fsects[i] = (c1, c2)
    end

    fsectperm = _sectors_sortperm(fsects)

    legs = (fuse(+1, rowlegs), fuse(-1, collegs))

    sects, sizes = _allsectorsandsizes(A.charge, legs)

    nzblks = Vector{Matrix{T}}()
    pointer = 1
    for index in 1:length(sects)
        blk = Matrix{T}(undef, sizes[index])
        #range = [1:m for m in sizes[index]]

        rlim, clim = sizes[index]
        pc = 1
        while pc <= clim
            pr = 1
            fdc = 0
            while pr <= rlim
                Ablk = permutedims(A.nzblks[perm1][fsectperm][pointer], idxperm)
                s = size(Ablk)
                fdr = prod(s[1:n_row])
                fdc = prod(s[n_row+1:N])
                blk[pr:pr+fdr-1, pc:pc+fdc-1] =
                    reshape(Ablk, fdr, fdc)
                pr += fdr
                pointer += 1
            end
            pc += fdc
        end
        push!(nzblks, blk)
    end
    SymMatrix(A.charge, legs, sects, nzblks)
end

# function SymMatrix(A::AbstractSymTensor,
#                    rowidxs::Vector{Int},
#                    colidxs::Vector{Int})

#     N = numoflegs(A)
#     idxperm = [rowidxs; colidxs]
#     sort(idxperm) == collect(1:N) || error("Incorrect index set for conversion to SymMatrix!")
#     length(rowidxs) == 0 && error("Zero row for matrix not allowed!")
#     length(colidxs) == 0 && error("Zero col for matrix not allowed!")

#     pA = permutelegs(A, idxperm)
#     SymMatrix(fuselegs(
#         fuselegs(pA, -1, length(rowidxs)+1, length(colidxs)),
#         +1, 1, length(rowidxs)))

# end

function _arecontractible(l1::STLeg, l2::STLeg)
    if l1.sign == -l2.sign
        return l1.chrs == l2.chrs && l1.dims == l2.dims
    end
    false
end

function *(A::SymVector, B::SymVector) where {T<:Number}

    _arecontractible(A.legs[1], B.legs[1]) ||
        error("not contractible!", A.legs[1], " and ", B.legs[1])
    signs(A.legs) == (-1,) && signs(B.legs) == (+1,) ||
        error("* for SymVector only defined for contraction!")

    return sum(A.nzblks[1] .* B.nzblks[1])
end

function *(A::AbstractSymMatrix,
           B::AbstractSymMatrix)

    _arecontractible(A.legs[2], B.legs[1]) ||
        error("not contractible! ", A.legs[2], " and ", B.legs[1])

    eltype(A) == eltype(B) || error("eltypes don't match!")
    T = eltype(A)
    ## NOTE: assume the set of charges for the vector space to be
    ## contracted is {c1,..,ci,...,cn} then each sector in A is given
    ## by (CA+ci, ci) and each sector in B by (ci, ci-CB) and since
    ## both are sorted based on the last charge in sector, then they
    ## are sorted the same. So we simply need to multiply them here
    ## but take into account that it is possible for the matrices to
    ## have different number of sectors allowed due to the other
    ## vector space they map from/to.

    ## NOTE: After contraction it is possible that some new sectors
    ## (that weren't allowed with the contraction indexes) are now
    ## allowed. So they have to be manually added and initialized to
    ## zero. I am actually not very happy with this; originally I
    ## decided to not include these new sectors at all but then
    ## realized at some instances the existance of them is necessary,
    ## well they are allowed sectors anyways, right! Still don't know
    ## what is the best approach.

    charge = A.charge + B.charge
    legs = (A.legs[1], B.legs[2])

    sects, sizes = _allsectorsandsizes(charge, legs)
    n_sectors = length(sects)
    nzblks = Vector{Matrix{T}}(undef, n_sectors)

    bmax = length(B.sects)
    amax = length(A.sects)
    a, b = 1, 1
    #println(A)
    #println(B)
    for i=1:n_sectors
        c1, c2 = sects[i]
        a += searchsortedfirst(A.sects[a:end], (c1, c1-A.charge)) - 1
        b += searchsortedfirst(B.sects[b:end], (B.charge+c2, c2)) - 1

        if a <=amax && b <=bmax && c1 == A.sects[a][1] && c2 == B.sects[b][2]
            nzblks[i] = A.nzblks[a] * B.nzblks[b]
        else
            nzblks[i] = zeros(T, sizes[i])
        end
    end
    SymMatrix(charge, legs, sects, nzblks)
end

###TODO: This is the above function with a version of fuselegs (that fuses both of them at once!)
function SymMatrix2(A::AbstractSymTensor{T, N},
                    rowidxs::Vector{Int},
                    colidxs::Vector{Int}) where {T<:Number, N}

    idxperm = [rowidxs; colidxs]
    sort(idxperm) == collect(1:N) || error("Can't convert to SymMatrix!")
    length(rowidxs) == 0 && error("Zero row for matrix not allowed!")
    length(colidxs) == 0 && error("Zero col for matrix not allowed!")

    pA = permutelegs(A, idxperm)

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
#         println(_A)
#         println(B)
#         println(remsB, " ", consB)
#         println(_B)
#         println(_A*_B)
#     end

# function contract(A::SymTensor{ComplexF64, N}, idxA::NTuple{N, Int},
#                   B::SymTensor{Float64, M}, idxB::NTuple{M, Int}) where{N, M}
#     contract(A, idxA, convert(SymTensor{ComplexF64, M}, B), idxB)
# end
# contract A and B using idxA and idxB. Note idxA and idxB are the
# corresponding leg numbers in each tensor that are to be contracted

# ##TODO: make these conversion more sane!!!
# function *(sten1::SymTensor{ComplexF64, 2},
#            sten2::SymTensor{Float64, 2})
#     *(sten1, convert(SymTensor{ComplexF64, 2}, sten2))
# end

# function *(sten1::SymTensor{Float64, 2},
#            sten2::SymTensor{ComplexF64, 2})
#     *(convert(SymTensor{ComplexF64, 2}, sten1), sten2)
# end

# function _arecontractible(l1::STLeg, l2::STLeg)
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
