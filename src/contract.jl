"""
    contract(A, idxA, B, idxB)

contract the symmetric tensor A, B. The to-be-contracted
indexes are shown with negative integers while the remaining
indexes for the result are shown with positive integers. For
example if there are l negative numbers in both A and B (should
be numbered from -1 to -l) then the rest of the indexes in A and B
combined should be 1:N+M-2n.
"""
function contract(A     :: AbstractSymTensor{S,T1,N},
                  idxA  :: NTuple{N, Int},
                  B     :: AbstractSymTensor{S,T2,M},
                  idxB  :: NTuple{M, Int}) where{S,T1,T2, N, M}

    remsA, consA, tofinalsA = _contract_index_perm(idxA)
    remsB, consB, tofinalsB = _contract_index_perm(idxB)

    ## TODO: techniqually here the contraction vector space should be made
    ## once! But now is made twice, so can we do better!

    # TODO: check
    # compatibility println(remsA, consA, remsB, consB)
    if length(remsA) == 0
        _A = SymVector(A, consA, false)
        if length(remsB) == 0
            _B = SymVector(B, consB, true)
            return _A * _B
        else
            _B = SymMatrix(B, consB, remsB)
            return permutelegs(unfuseleg(_A*_B, 2, B.legs[remsB]),
                               invperm(tofinalsB))
        end
    end

    _A = SymMatrix(A, remsA, consA)

    if length(remsB) == 0
        _B = SymVector(B, consB, true)
        return permutelegs(unfuseleg(_A*_B, 1, A.legs[remsA]),
                           invperm(tofinalsA))
    end

    _B = SymMatrix(B, consB, remsB)
    # println(_A)
    # println(_B)
    return permutelegs(SymTensor(_A * _B, A.space[remsA], B.space[remsB]),
                       #unfuseleg(unfuseleg(_A * _B, 1, A.legs[remsA]), length(remsA)+1, B.legs[remsB]),
                       invperm([tofinalsA; tofinalsB]))
end

function contract(A     :: AbstractArray{T1,N},
                  idxA  :: NTuple{N, Int},
                  B     :: AbstractArray{T2,M},
                  idxB  :: NTuple{M, Int}) where{S,T1,T2, N, M}
    remsA, consA, tofinalsA = _contract_index_perm(idxA)
    remsB, consB, tofinalsB = _contract_index_perm(idxB)
    permA = [remsA, consA]
    permB = [consB, remsB]
    _A = reshape(permutedims(A, permA), prod(size(A)[remsA]), prod(size(A)[consA]))
    _B = reshape(permutedims(B, permB), prod(size(B)[consB]), prod(size(B)[remsB]))
    permutedims(reshape(_A*_B, size(A)[remsA]..., size(B)[remsB]),
                invperm([rofinalsA; tofinalsB]))
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
                   perm::Vector{Int},
                   _dual::Bool)
    N = rank(A)
    sort(perm) == collect(1:N) ||
        error("Incorrect index set for conversion to SymVector!")

    pA = permutelegs(A, perm)
    fuselegs(pA, 1, N, _dual)
end

function SymMatrix(A::AbstractSymTensor,
                   rowidxs::Vector{Int},
                   colidxs::Vector{Int})

    N = rank(A)
    idxperm = [rowidxs; colidxs]
    sort(idxperm) == collect(1:N) ||
        error("Incorrect index set for conversion to SymMatrix!")
    length(rowidxs) == 0 && error("empty row for matrix not allowed!")
    length(colidxs) == 0 && error("empty col for matrix not allowed!")
    n = length(rowidxs)

    T = eltype(A)
    S = vtype(A)

    rowdualinfo = isdual.(A.space[rowidxs])
    coldualinfo = isdual.(A.space[colidxs])
    csects = [Sector(x[idxperm]) for x in sectors(A)]
    sperm1 = sortperm(csects)
    fsects = Vector{Sector{S, 2}}(undef, length(csects))
    for i in eachindex(csects)
        sector = csects[sperm1][i]
        c1 = sum(Sector(sector[1:n]), rowdualinfo)
        c2 = inv(sum(Sector(sector[n+1:N]), coldualinfo))
        fsects[i] = Sector(c1, c2)
    end
    sperm2 = sortperm(fsects)
    semits = collect(onlysemitokens(A.blocks))[sperm1][sperm2]

    space = (fuse(false, A.space[rowidxs]), fuse(true, A.space[colidxs]))
    #println(space)

    sects, sizes = _allsectorsandsizes(A.charge, space)

    blocks = SortedDict{Sector{S, 2}, Matrix{T}}()
    pointer = 1
    for index in 1:length(sects)
        #println("SymMat for sector $(sects[index])")
        blk = Matrix{T}(undef, sizes[index])
        #range = [1:m for m in sizes[index]]

        rlim, clim = sizes[index]
        pc = 1
        while pc <= clim
            pr = 1
            fdc = 0
            while pr <= rlim
                #println("sector $(sectors(A)[sperm1][sperm2][pointer])")
                Ablk = permutedims(A.blocks[semits[pointer]], idxperm)
                s = size(Ablk)
                fdr = prod(s[1:n])
                fdc = prod(s[n+1:N])
                blk[pr:pr+fdr-1, pc:pc+fdc-1] =
                    reshape(Ablk, fdr, fdc)
                pr += fdr
                pointer += 1
            end
            pc += fdc
        end
        sector = sects[index]
        blocks[sector] = blk
    end
    SymMatrix(A.charge, space, blocks)
end

# function SymMatrix(A::AbstractSymTensor,
#                    rowidxs::Vector{Int},
#                    colidxs::Vector{Int})

#     N = rank(A)
#     idxperm = [rowidxs; colidxs]
#     sort(idxperm) == collect(1:N) || error("Incorrect index set for conversion to SymMatrix!")
#     length(rowidxs) == 0 && error("Zero row for matrix not allowed!")
#     length(colidxs) == 0 && error("Zero col for matrix not allowed!")

#     pA = permutelegs(A, idxperm)
#     SymMatrix(fuselegs(
#         fuselegs(pA, -1, length(rowidxs)+1, length(colidxs)),
#         +1, 1, length(rowidxs)))

# end

function SymTensor(A     :: SymMatrix,
                   rlegs :: NTuple{N, VectorSpace{S}},
                   clegs :: NTuple{M, VectorSpace{S}}) where {S, N, M}
    T = eltype(A)
    #sperm = sortperm(A.sects, by=x->x[2])
    #csects = A.sects[sectperm]
    blocks = SortedDict{Sector{S, N+M}, Array{T, N+M}}()
    csects = sectors(A)
    #    oldcharge = 0
    rpats, rsizes = Vector{Sector{S, M}}(), Vector{NTuple{M, Int}}()
    cpats, csizes = Vector{Sector{S, M}}(), Vector{NTuple{M, Int}}()
    isduals = isdual.(A.space)
    for (s,blk) in A.blocks
        c1, c2 = s.charges
        c1 = isdual(A.space[1]) ? inv(c1) : c1
        c2 = isdual(A.space[2]) ? inv(c2) : c2
        rpats, rsizes = _allsectorsandsizes(c1, rlegs)
        cpats, csizes = _allsectorsandsizes(c2, clegs)

        pc = 0
        for cpatidx in eachindex(cpats)
            csl = prod(csizes[cpatidx])
            pr = 0
            for rpatidx in eachindex(rpats)
                rsl = prod(rsizes[rpatidx])
                blocks[Sector(rpats[rpatidx]..., cpats[cpatidx]...)] =
                    reshape(blk[pr+1:pr+rsl, pc+1:pc+csl],
                            rsizes[rpatidx]..., csizes[cpatidx]...)
                pr += rsl
            end
            pc += csl
        end
    end

    #TODO: now check to see if the new legs can fuse into the original leg
    space = (rlegs...,clegs...)
    SymTensor(A.charge, space, blocks)
end

function *(A::SymVector, B::SymVector)
    isdual(A.space[1], B.space[1]) ||
        error("not contractible!", A.space[1], " and ", B.space[1])

    semits = collect(onlysemitokens(A.blocks))
    return sum(A.blocks[semits[1]] .* B.blocks[semits[1]])
end

function *(A::AbstractSymMatrix,
           B::AbstractSymMatrix)

    isdual(A.space[2], B.space[1]) ||
        error("not contractible! ", A.space[2], " and ", B.space[1])

    S = vtype(A)
    T = promote_type(eltype(A), eltype(B))
    T == Union{} &&
        error("can't promote_type! $(eltype(A)) , $(eltype(B))")
    ## NOTE: assume the set of charges for the vector space to be
    ## contracted is {c1,..,ci,...,cn} then each sector in A is given
    ## by (CA+ci, -ci) and each sector in B by (ci, -ci+CB) and since
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

    ## TODO: Also it would be nice if there is the option of dropping
    ## the charges that only show up in complete zero sectors! This is
    ## particularly useful is Gutzwiller projection thing

    charge = A.charge + B.charge
    space = (A.space[1], B.space[2])

    sects, sizes = _allsectorsandsizes(charge, space)
    n_sectors = length(sects)
    blocks = SortedDict{Sector{S, 2}, Matrix{T}}()

    semitsA = collect(onlysemitokens(A.blocks))
    semitsB = collect(onlysemitokens(B.blocks))
    bmax = length(semitsB)
    amax = length(semitsA)
    a, b = 1, 1
    #println(A)
    #println(B)
    for i=1:n_sectors
        c1, c2 = sects[i]
        _c1 = isdual(A.space[1]) ? inv(c1) : c1
        _c2 = isdual(B.space[2]) ? inv(c2) : c2
        if isdual(A.space[2])
            sectorA = Sector(c1, inv(A.charge)+_c1)
            sectorB = Sector(B.charge-_c2, c2)
        else
            sectorA = Sector(c1, A.charge-_c1)
            sectorB = Sector(inv(B.charge)+_c2, c2)
        end
        if haskey(A.blocks, sectorA) && haskey(B.blocks, sectorB)
            blocks[sects[i]] = A.blocks[sectorA] * B.blocks[sectorB]
        else
            blocks[sects[i]] = zeros(T, sizes[i])
        end
    end
    SymMatrix(charge, space, blocks)
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

    #     fusedsectperm = sortperm(fusedsectors)
    #     new_sectors, refs = uniquesorted(fusedsectors[fusedsectperm])
    #     refs = refs[invperm(fusedsectperm)]

    #     # the infos are: new_sectors, refs, patranges

    #     # make the new leg
    #     fleg = U1Space(sign, unzip([(k, fchrdict[k].dim) for k in sort(collect(keys(fchrdict)))])...)

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

# function _arecontractible(l1::U1Space, l2::U1Space)
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
