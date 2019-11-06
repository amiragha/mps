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

function permutelegs(A::AbstractSymTensor,
                     perm::Vector{Int})
    N = numoflegs(A)
    sort(perm) == collect(1:N) || error("Not a valid permutation!")
    perm == collect(1:N) && return A

    sectperm = _sectors_sortperm(A.sects,  by=x->x[perm])
    sects = [sect[perm] for sect in A.sects[sectperm]]

    ## NOTE: The line below because for some reason the permutedims of
    ## Diagonal returns a sparse matrix in Julia 1.2.x
    typeof(A) <: SymDiagonal &&
        return typeof(A)(A.charge, A.legs[perm], sects, A.nzblks)

    typeof(A)(A.charge, A.legs[perm], sects,
              [permutedims(nzblk, perm) for nzblk in A.nzblks[sectperm]])
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

"""
    fuselegs

fuses two or `n` consequative legs of a SymTensor.

The fuse operation works as follows: We have assumed the sectors are
sorted, so they stay sorted even after the fusion is done!

"""
function fuselegs(A::AbstractSymTensor,
                  sign::Int,
                  l::Int,
                  n::Int=2;
                  debug=false)

    abs(sign) == 1 || "sign has to be ± 1!"
    N = numoflegs(A)
    0 < l < N+2-n  || "fuselegs $(l+n-1) vs $N"

    if n == 1
        if sign == A.legs[l].sign
            return A
        else
            return negateleg(A, l)
        end
    end

    T = eltype(A)

    signs = [A.legs[i].sign for i in l:l+n-1]

    csects = A.sects
    fsects = Vector{NTuple{N-n+1, Int}}(undef, length(A.sects))
    for i in eachindex(csects)
        csect = csects[i]
        #fc = sign * sum(signs .* csect[l:l+n-1]) # this is slower than below!
        fc = sign * sum([signs[i] * csect[l+i-1] for i=1:n])
        fsects[i] = (csect[1:l-1]..., fc, csect[l+n:end]...)
    end

    fsectperm = _sectors_sortperm(fsects)
    newsects, refs = uniquesorted(fsects[fsectperm])
    refs = refs[invperm(fsectperm)]

    fleg = fuse(sign, A.legs[l:l+n-1])
    newlegs = Tuple([A.legs[1:l-1]..., fleg, A.legs[l+n:end]...])

    newnzblks = Vector{Array{T, N-n+1}}(undef, length(newsects))
    for i in eachindex(newsects)
        sect = newsects[i]
        newnzblks[i] = zeros(T, [getdim(newlegs[c], sect[c])
                                 for c in 1:N-n+1]...)
    end

    pointers = ones(Int, length(newnzblks))
    for i in eachindex(A.nzblks)
        s = size(A.nzblks[i])
        p = pointers[refs[i]]
        fd = prod(s[l:l+n-1])
        # This is slower than below!
        #axesA = axes(A.nzblks[i])
        # newnzblks[refs[i]][axesA[1:l-1]..., p:p+fd-1, axesA[l+n:N]...] =
        #     reshape(A.nzblks[i], s[1:l-1]..., fd, s[l+n:end]...)
        newnzblks[refs[i]][[1:s[i] for i=1:l-1]..., p:p+fd-1, [1:s[i] for i=l+n:N]...] =
            reshape(A.nzblks[i], s[1:l-1]..., fd, s[l+n:end]...)
        pointers[refs[i]] += fd
    end

    SymTensor(A.charge, newlegs, newsects, newnzblks)
end

function delinsert(tuple::NTuple{N, T}, items::NTuple{M, T}, index::Int) where{T, N, M}
    (tuple[1:index]..., items..., tuple[index+1:N])
end

# defuses a leg into some given legs
function unfuseleg(A    :: AbstractSymTensor{T, N},
                   l    :: Int,
                   legs :: NTuple{M, STLeg}) where {T<:Number, N, M}
    0 < l <= N || error("integer l not in range $l, $N")
    #println(A.legs)
    #println(legs)
    if length(legs) == 1
        if legs[1].sign == A.legs[l].sign
            A.legs[l] == legs[1] || error("Not the same legs", A.legs[l], legs[1])
            return A
        else
            Aleg = negate(A.legs[l])
            Aleg == legs[1] || (error("Not the same legs", Aleg, legs[1]))
            return negateleg(A, l)
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
