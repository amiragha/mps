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

function permutelegs(A::AbstractSymTensor{T, N},
                     perm::Vector{Int}) where{T, N}

    sort(perm) == collect(1:N) || error("Not a valid permutation!")
    sectperm = _sectors_sortperm(A.sects,  by=x->x[perm])

    sects = [sect[perm] for sect in A.sects[sectperm]]
    nzblks = [permutedims(nzblk, perm) for nzblk in A.nzblks[sectperm]]

    typeof(A)(A.charge, A.legs[perm], sects, nzblks)
end

"""
    fuse_set(sect, parts)

fuse the values in the `sect` using the binary operation `op` by to
the partitions given by `parts`

"""
function fuse_set(op,
                  sect::NTuple{N, Tc},
                  partsizes::Vector{Int}) where {N, Tc<:Number}
    ###TODO: I should be able to code this nicer, probably using some
    ###functional stuff like map, fold, etc
    result = Tc[]
    p = 0
    for n in partsizes
        push!(result, foldl(op, sect[p+1:p+n]))
        p += n
    end
    Tuple(result)
end



# fuse n consequative legs `l` and `l+n-1` to a new leg with sign `sign`
# function fuse_conseqlegs(sten::SymTensor{T, N},
#                          sign::Int, l::Int, n::Int=2; debug=false) where{T<:Number, N}
#     @assert abs(sign) == 1
#     #@assert 0 < l < N+n-2
#     @assert 0 < l < N+2-n

#     fusedsectors = NTuple{N-n+1, Int}[]
#     patranges = UnitRange{Int}[]

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
# end

# fuse n consequative legs `l` and `l+n-1` to a new leg with sign `sign`

"""
    fuselegs

fuses two or `n` consequative legs of a SymTensor.

The fuse operation works as follows: We have assumed the sectors are
sorted, so they stay sorted even after the fusion is done!

"""
function fuselegs(sten::SymTensor{Tv, N},
                  sign::Int,
                  l::Int,
                  n::Int=2;
                  debug=false) where{Tv<:Number, N}
    abs(sign) == 1 || "sign has to be ± 1!"
    0 < l < N+2-n  || "fuselegs $(l+n-1) vs $N"

    fsects = NTuple{N-n+1, Int}[]
    patranges = UnitRange{Int}[]
    signs = Tuple(sten.legs[i].sign for i in l:l+n-1)
    fleg, fcdict = fuse(sign, sten.legs[l:l+n-1])

    csects = sten.sects
    for i in eachindex(csects)
        csect = csects[i]
        pat = Tuple(csect[l:l+n-1])
        spat = sign .* signs .* pat
        fc = sum(spat)
        patrange = fcdict[fc].pats[spat]
        fsect = (csect[1:l-1]..., fc, csect[l+n:end]...)
        push!(fsects, fsect)
        push!(patranges, patrange)
    end

    fsectperm = _sectors_sortperm(fsects)
    newsects, refs = uniquesorted(fsects[fsectperm])
    refs = refs[invperm(fsectperm)]

    newlegs = Tuple(vcat([sten.legs[1:l-1]...], fleg, [sten.legs[l+n:end]...]))

    newnzblks = Array{Tv, N-n+1}[]
    for sect in newsects
        push!(newnzblks,zeros(Tv, [getdim(newlegs[c], sect[c])
                                   for c in 1:N-n+1]...))
    end

    for i in eachindex(sten.nzblks)
        s = size(sten.nzblks[i])
        newnzblks[refs[i]][[1:s[i] for i=1:l-1]..., patranges[i],[1:s[i] for i=l+n:N]...] =
            reshape(sten.nzblks[i], s[1:l-1]..., prod(s[l:l+n-1]), s[l+n:end]...)
    end

    SymTensor(sten.charge, newlegs, newsects, newnzblks)
end

function delinsert(tuple::NTuple{N, T}, items::NTuple{M, T}, index::Int) where{T, N, M}
    (tuple[1:index]..., items..., tuple[index+1:N])
end

# defuses a leg into some given legs
function defuse_leg(A::AbstractSymTensor{T, N},
                    l::Int,
                    legs::NTuple{M, STLeg}) where {T<:Number, N, M}
    0 < l <= N || error("integer l not in range $l, $N")
    #println(sten.legs)
    #println(legs)
    if length(legs) == 1
        if legs[1].sign == A.legs[l].sign
            A.legs[l] == legs[1] || error("Not the same legs", A.legs[l], legs[1])
            return A
        else
            error("change of leg sign is not longer supported in defuse_leg!")
            # Aleg = change_sign(A.legs[l])
            # Aleg == legs[1] || (error("Not the same legs", Aleg, legs[1]))
            # return change_legsign(A, l)
        end
    end

    sects = NTuple{N+M-1, Int}[]
    nzblks = Array{T, N+M-1}[]

    sign = A.legs[l].sign
    csects = A.sects
    fchrdict = Dict{Int, Tuple{Vector{NTuple{M,Int}}, Vector{NTuple{M,Int}}}}()
    for i in eachindex(csects)
        charge = sign * csects[i][l]
        if haskey(fchrdict, charge)
            pats, patdims = fchrdict[charge]
        else
            pats, patdims = _allsectorsandsizes(charge, legs)
            #println(charge, pats, patdims)
            ## NOTE: sorting is not longer required because
            ## allsectorsandsizes return sorted pats
            # patperms = _sectors_sortperm(pats)
            # pats, patdims = pats[patperms], patdims[patperms]
            #println(pats, patdims)
            fchrdict[charge] = (pats, patdims)
        end
        old_nzblock =  A.nzblks[i]
        s = size(old_nzblock)
        pivot = 0
        for patidx in eachindex(pats)
            sl = prod(patdims[patidx])
            push!(sects, (csects[i][1:l-1]..., pats[patidx]..., csects[i][l+1:N]...))
            push!(nzblks, reshape(old_nzblock[[1:s[n] for n=1:l-1]...,pivot+1:pivot+sl,[1:s[n] for n=l+1:N]...],
                                  s[1:l-1]...,patdims[patidx]...,s[l+1:N]...))
            pivot += sl
        end
    end


    # now check to see if the new legs can fuse into the original leg
    ###TODO: there should be an options to remove this step for
    ###efficientcy when it is not needed!
    for charge in keys(fchrdict)
        cidx = findall(A.legs[l].chrs .== sign * charge)[1]
        A.legs[l].dims[cidx] == sum([prod(patdims) for patdims in fchrdict[charge][2]]) ||
            error("For defuse, dimensions don't match for charge " , sign*charge)
        #, " ", A.legs[l].dims[cidx], " == ", fchrdict[charge]

    end

    #println(fchrdict)
    new_legs = (A.legs[1:l-1]...,legs...,A.legs[l+1:end]...)
    #println(A.sects)
    #println(sects)

    ## NOTE: It is necessary to sort the sectors here, because while
    ## the patterns were sorted whitin each fused charge, they are not
    ## sorted when are the charges are concerned
    sectperm = _sectors_sortperm(sects)
    SymTensor(A.charge, new_legs, sects[sectperm], nzblks[sectperm])
end
