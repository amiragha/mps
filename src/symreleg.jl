"""
    uniquesorted(A)

For a sorted vector `A` return the tuple of (unqA, refs) where unqA is
the unique elements in A and refs are the indexes of unqA for each
element of A

"""
function uniquesorted(A::Vector{T}) where{T}
    unqA = T[]
    refs = Int[]
    push!(unqA, A[1])
    idx = 1
    push!(refs, idx)
    for i=2:length(A)
        if A[i] != unqA[end]
            push!(unqA, A[i])
            idx += 1
            push!(refs, idx)
        else
            push!(refs, idx)
        end
    end
    unqA, refs
end

unzip(a) = map(x->getfield.(a,x), fieldnames(eltype(a)))

################################
#### RELEG #####################
################################

"permute legs, It is important to bring the tensor into some
cononical(convention) tensor product fusion tree"
function permutelegs(A::AbstractSymTensor,
                     perm::Vector{Int})
    N = rank(A)
    sort(perm) == collect(1:N) || error("Not a valid permutation!")
    perm == collect(1:N) && return A

    sectperm = sortperm(A.sects,  by=x->x[perm])
    sects = [sect[perm] for sect in A.sects[sectperm]]

    ## NOTE: The line below because for some reason the permutedims of
    ## Diagonal returns a sparse matrix in Julia 1.2.x
    typeof(A) <: SymDiagonal &&
        return typeof(A)(A.charge, A.legs[perm], sects, A.nzblks)

    typeof(A)(A.charge, A.legs[perm], sects,
              [permutedims(nzblk, perm) for nzblk in A.nzblks[sectperm]])
end

"""
        fuselegs(A, sign, l, n)

fuses two or `n` consequative legs of SymTensor `A` starting at leg `l`
into a new leg with sign (direction) `sign`.

The fuse operation works as follows: We have assumed the sectors are
sorted, so they stay sorted even after the fusion is done!

"""
function fuselegs(A::AbstractSymTensor,
                  l::Int,
                  n::Int=2)

    N = rank(A)
    0 < l < N+2-n  ||
        error("fuselegs $l to $(l+n-1) for size $N")

    if n < 2
        return A
    end

    T = eltype(A)

    csects = collect(keys(A.data))
    semits = collect(onlysemitokens(A.data))
    S = vtype(A)
    fsects = Vector{Sector{S, N-n+1}}(undef, length(csects))
    for i in eachindex(csects)
        csect = csects[i]
        fsects[i] = Sector(csect.charges[1:l-1]...,
                           sum(csect.charges[l:l+n-1]),
                           csect.charges[l+n:end]...)
    end

    fsectperm = sortperm(fsects)
    semits = semits[fsectperm]

    fspace = fuse(A.space[l:l+n-1])
    space = Tuple([A.space[1:l-1]..., fspace, A.space[l+n:end]...])

    sects, sizes = _allsectorsandsizes(A.charge, space)

    data = SortedDict{Sector{S, N-n+1}, Array{T, N-n+1}}()
    pointer = 1
    for index in 1:length(sects)
        block = Array{T, N-n+1}(undef, sizes[index])
        range = [1:m for m in sizes[index]]

        sizel = sizes[index][l]
        p = 1
        while p <= sizel
            s = size(A.data[semits[pointer]])
            fd = prod(s[l:l+n-1])
            range[l] = p:p+fd-1
            block[range...] =
                reshape(A.data[semits[pointer]], s[1:l-1]..., fd, s[l+n:end]...)
            p += fd
            pointer += 1
        end
        sector = sects[index]
        data[sector] = block
    end
    SymTensor(A.charge, space, data)
end

function delinsert(tuple::NTuple{N, T}, items::NTuple{M, T}, index::Int) where{T, N, M}
    (tuple[1:index]..., items..., tuple[index+1:N])
end

# defuses a leg into some given legs
splitleg(A::AbstractSymTensor,  l::Int, space::U1Space...) = splitleg(A, l, space)
function splitleg(A    :: AbstractSymTensor,
                  l    :: Int,
                  space :: NTuple{M, VectorSpace{S}}) where {M, S}
    T = eltype(A)
    N = rank(A)
    0 < l <= N || error("integer l not in range $l, $N")
    vtype(A) == S || throw("SpaceMismatch $(vtype(A)), $S")
    if length(space) == 1
        A.space[l] == first(space) ||
            throw("SpaceMismatch $(A.space[l]), $(first(space))")
        return A
    end

    data = SortedDict{Sector{S, N+M-1}, Array{T, N+M-1}}()

    ##TODO: explain what tensor product structure change the below
    ##procedure amounts to!
    ##We first sort the sectors based on the charge the leg to be
    ##unfused. Then for each charge find all sectors and sizes!
    sects = collect(keys(A.data))
    semits = collect(onlysemitokens(A.data))
    sectperm = sortperm(sects, by=x->x.charges[l])
    csects = sects[sectperm]
    semits = semits[sectperm]
    oldcharge = 0
    pats, sizes = Vector{Sector{S, M}}(), Vector{NTuple{M, Int}}()
    for i in eachindex(csects)
        #charge = sign * csects[i][l]
        charge = csects[i][l]
        if i == 1 || charge != oldcharge
            pats, sizes = _allsectorsandsizes(charge, space)
            oldcharge = charge
        end

        old_nzblock =  A.data[semits[i]]
        s = size(old_nzblock)
        pivot = 0
        for patidx in eachindex(pats)
            sl = prod(sizes[patidx])
            sector = Sector(csects[i][1:l-1]..., pats[patidx]..., csects[i][l+1:N]...)
            data[sector] = reshape(old_nzblock[[1:s[n] for n=1:l-1]...,
                                                pivot+1:pivot+sl,[1:s[n] for n=l+1:N]...],
                                    s[1:l-1]...,sizes[patidx]...,s[l+1:N]...)
            pivot += sl
        end
    end

    #TODO: now check to see if the new legs can fuse into the original leg
    new_space = (A.space[1:l-1]...,space...,A.space[l+1:end]...)

    SymTensor(A.charge, new_space, data)
end

# """
#     fuse_set(sect, parts)

# fuse the values in the `sect` using the binary operation `op` by to
# the partitions given by `parts`

# """
# function fuse_set(op,
#                   sect::NTuple{N, Tc},
#                   partsizes::Vector{Int}) where {N, Tc<:Number}
#     ###TODO: I should be able to code this nicer, probably using some
#     ###functional stuff like map, fold, etc
#     result = Tc[]
#     p = 0
#     for n in partsizes
#         push!(result, foldl(op, sect[p+1:p+n]))
#         p += n
#     end
#     Tuple(result)
# end


# # defuses a leg into some given legs
# function unfuseleg(A    :: AbstractSymTensor{T, N},
#                    l    :: Int,
#                    legs :: NTuple{M, U1Space}) where {T<:Number, N, M}
#     0 < l <= N || error("integer l not in range $l, $N")
#     #println(A.legs)
#     #println(legs)
#     if length(legs) == 1
#         if legs[1].sign == A.legs[l].sign
#             A.legs[l] == legs[1] || error("Not the same legs", A.legs[l], legs[1])
#             return A
#         else
#             Aleg = negate(A.legs[l])
#             Aleg == legs[1] || (error("Not the same legs", Aleg, legs[1]))
#             return negateleg(A, l)
#         end
#     end

#     sects = NTuple{N+M-1, Int}[]
#     nzblks = Array{T, N+M-1}[]

#     sign = A.legs[l].sign
#     csects = A.sects
#     fchrdict = Dict{Int, Tuple{Vector{NTuple{M,Int}}, Vector{NTuple{M,Int}}}}()
#     for i in eachindex(csects)
#         charge = sign * csects[i][l]
#         if haskey(fchrdict, charge)
#             pats, patdims = fchrdict[charge]
#         else
#             pats, patdims = _allsectorsandsizes(charge, legs)
#             #println(charge, pats, patdims)
#             ## NOTE: sorting is not longer required because
#             ## allsectorsandsizes return sorted pats
#             # patperms = sortperm(pats)
#             # pats, patdims = pats[patperms], patdims[patperms]
#             #println(pats, patdims)
#             fchrdict[charge] = (pats, patdims)
#         end
#         old_nzblock =  A.nzblks[i]
#         s = size(old_nzblock)
#         pivot = 0
#         for patidx in eachindex(pats)
#             sl = prod(patdims[patidx])
#             push!(sects, (csects[i][1:l-1]..., pats[patidx]..., csects[i][l+1:N]...))
#             push!(nzblks, reshape(old_nzblock[[1:s[n] for n=1:l-1]...,pivot+1:pivot+sl,[1:s[n] for n=l+1:N]...],
#                                   s[1:l-1]...,patdims[patidx]...,s[l+1:N]...))
#             pivot += sl
#         end
#     end

#     # now check to see if the new legs can fuse into the original leg
#     ###TODO: there should be an options to remove this step for
#     ###efficientcy when it is not needed!
#     for charge in keys(fchrdict)
#         cidx = findall(A.legs[l].chrs .== sign * charge)[1]
#         A.legs[l].dims[cidx] == sum([prod(patdims) for patdims in fchrdict[charge][2]]) ||
#             error("For defuse, dimensions don't match for charge " , sign*charge)
#         #, " ", A.legs[l].dims[cidx], " == ", fchrdict[charge]

#     end

#     #println(fchrdict)
#     new_legs = (A.legs[1:l-1]...,legs...,A.legs[l+1:end]...)
#     #println(A.sects)
#     #println(sects)

#     ## NOTE: It is necessary to sort the sectors here, because while
#     ## the patterns were sorted whitin each fused charge, they are not
#     ## sorted when are the charges are concerned
#     sectperm = sortperm(sects)
#     SymTensor(A.charge, new_legs, sects[sectperm], nzblks[sectperm])
# end
