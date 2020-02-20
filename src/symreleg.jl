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

"permute legs"
function permutelegs(A::AbstractSymTensor,
                     perm::Vector{Int})
    N = rank(A)
    sort(perm) == collect(1:N) || error("Not a valid permutation!")
    perm == collect(1:N) && return A

    # sectperm = sortperm(A.sects,  by=x->x[perm])
    # sects = [sect[perm] for sect in A.sects[sectperm]]

    ## NOTE: The line below because for some reason the permutedims of
    ## Diagonal returns a sparse matrix in Julia 1.2.x
    typeof(A) <: SymDiagonal &&
        return typeof(A)(A.charge, A.space[perm], A.blocks)

    typeof(A)(A.charge, A.space[perm],
              SortedDict([Sector(s[perm]) => permutedims(blk, perm)
                         for (s,blk) in A.blocks]))
end
function permutelegs(A::AbstractArray,
                     perm::Vector{Int})
    permutedims(A, perm)
end

"""
        fuselegs(A, l, n)

fuses two or `n` consequative legs of SymTensor `A` starting at leg `l`
into a new leg.

Below is not the current implementation!----
Fusion or tensor product is only defined for two vector spaces so for
larger ones it has to be preformed in a recursive fashion. Here we
follow the convention of tensor producting pairs from left to right,
which is (...((V1 ⊗ V2) ⊗ V3) ⊗ ...) ⊗ Vn)
----

We just pick the convention of sorting the blocks according to the
above tensor product rule and then lay the the memory for each
sector. Note that this is not equal to the recursive fusion protocol
but should be fine, because here the only thing that matters is
consistency! So we should only fuse and unfuse consistency, mix and
match leads to incorrect results obviously (I should think of a way to
prevent this).

"""
function fuselegs(A::AbstractSymTensor,
                  l::Int,
                  n::Int=2,
                  _dual::Bool=false)

    N = rank(A)
    0 < l < N+2-n  ||
        error("fuselegs $l to $(l+n-1) for size $N")

    if n < 2
        if _dual == isdual(A.space[l])
            return A
        else
            Vf = fuse(_dual, A.space[l])
            space = Tuple([A.space[1:l-1]..., Vf, A.space[l+n:end]...])
            SymTensor(A.charge, space,
                      SortedDict([Sector(s[1:l-1]...,inv(s[l]), s[l+1:N]...)=>b for (s,b) in A.blocks]))
        end
    end

    T = eltype(A)

    dualinfo = isdual.(A.space[l:l+n-1])
    csects = sectors(A)
    S = vtype(A)
    fsects = Vector{Sector{S, N-n+1}}(undef, length(csects))
    for i in eachindex(csects)
        sector = csects[i]
        fc = sum(Sector(sector[l:l+n-1]), dualinfo)
        if _dual
            fc = inv(fc)
        end
        fsects[i] = Sector(sector[1:l-1]...,
                           fc,
                           sector[l+n:N]...)
    end

    sperm = sortperm(fsects)
    semits = collect(onlysemitokens(A.blocks))[sperm]

    Vf = fuse(_dual, A.space[l:l+n-1])
    space = Tuple([A.space[1:l-1]..., Vf, A.space[l+n:end]...])

    sects, sizes = _allsectorsandsizes(A.charge, space)

    blocks = SortedDict{Sector{S, N-n+1}, Array{T, N-n+1}}()
    pointer = 1
    for index in 1:length(sects)
        blk = Array{T, N-n+1}(undef, sizes[index])
        range = [1:m for m in sizes[index]]

        sizel = sizes[index][l]
        p = 1
        while p <= sizel
            s = size(A.blocks[semits[pointer]])
            fd = prod(s[l:l+n-1])
            range[l] = p:p+fd-1
            blk[range...] =
                reshape(A.blocks[semits[pointer]], s[1:l-1]..., fd, s[l+n:end]...)
            p += fd
            pointer += 1
        end
        sector = sects[index]
        blocks[sector] = blk
    end
    #SymTensor(A.charge, space, blocks)
    SymTensor{S,T,N-n+1}(A.charge, space, blocks)
end
function fuselegs(A::AbstractArray,
                  l::Int,
                  n::Int=2,
                  _dual::Bool=false)
    sz = size(A)
    reshape(A, sz[1:l-1]..., prod(sz[l:l+n-1]), sz[l+n:end]...)
end

function delinsert(tuple::NTuple{N, T}, items::NTuple{M, T}, index::Int) where{T, N, M}
    (tuple[1:index]..., items..., tuple[index+1:N])
end

# defuses a leg into some given legs
splitleg(A::AbstractSymTensor,  l::Int, space::U1Space...) = splitleg(A, l, space)
function splitleg(A     :: AbstractSymTensor,
                  l     :: Int,
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

    blocks = SortedDict{Sector{S, N+M-1}, Array{T, N+M-1}}()

    ##TODO: explain what tensor product structure change the below
    ##procedure amounts to!
    ##We first sort the sectors based on the charge the leg to be
    ##unfused. Then for each charge find all sectors and sizes!

    sects = collect(keys(A.blocks))
    semits = collect(onlysemitokens(A.blocks))
    sectperm = sortperm(sects, by=x->x.charges[l])
    csects = sects[sectperm]
    semits = semits[sectperm]
    oldcharge = 0
    pats, sizes = Vector{Sector{S, M}}(), Vector{NTuple{M, Int}}()
    for i in eachindex(csects)
        #charge = sign * csects[i][l]
        charge = isdual(A.space[l]) ? inv(csects[i][l]) : csects[i][l]
        if i == 1 || charge != oldcharge
            pats, sizes = _allsectorsandsizes(charge, space)
            oldcharge = charge
        end

        old_nzblock =  A.blocks[semits[i]]
        s = size(old_nzblock)
        pivot = 0
        for patidx in eachindex(pats)
            sl = prod(sizes[patidx])
            sector = Sector(csects[i][1:l-1]..., pats[patidx]..., csects[i][l+1:N]...)
            blocks[sector] = reshape(old_nzblock[[1:s[n] for n=1:l-1]...,
                                                pivot+1:pivot+sl,[1:s[n] for n=l+1:N]...],
                                    s[1:l-1]...,sizes[patidx]...,s[l+1:N]...)
            pivot += sl
        end
    end

    #TODO: now check to see if the new legs can fuse into the original leg
    new_space = (A.space[1:l-1]...,space...,A.space[l+1:end]...)

    #SymTensor(A.charge, new_space, blocks)
    SymTensor{S,T,N+M-1}(A.charge, new_space, blocks)
end

function splitleg(A     :: AbstractArray,
                  l     :: Int,
                  space :: NTuple{M, TrivialVectorSpace}) where {M}
    sz = size(A)
    reshape(A, sz[1:l-1]..., dim.(space)..., sz[l+1:end]...)
end
