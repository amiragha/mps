mutable struct SymTensor{Tv<:Number, N}
    charge :: Int     # total charge of the tensor
    legs :: NTuple{N, STLeg}
    sects :: Vector{NTuple{N, Int}}    # possible nonzero sectors
    nzblks :: Vector{AbstractArray{Tv, N}}     # stored nonzero blocks

    # TODO: add a full sanity check for tests!
    function SymTensor(charge::Int,
                       legs::NTuple{N, STLeg},
                       sects::Vector{NTuple{N, Int}},
                       nzblks::Vector{<:AbstractArray{Tv, N}}) where{Tv<:Number, N}
        issorted(sects, lt=_sector_less_than) ||  error("sectors not sorted!")
        #legs = _trimlegs(legs, sects)
        new{Tv, N}(charge, legs, sects, nzblks)
    end
end

const SymMatrix{Tv} = SymTensor{Tv, 2} where {Tv<:Number}
@inline numoflegs(s::SymTensor{Tv, N}) where{Tv<:Number, N} = N

@inline size(s::SymTensor) = Tuple(fulldims(l) for l in s.legs)
@inline size(s::SymTensor, l::Int) = fulldims(s.legs[l])

#struct SymDiagonal{Tv} =
"""
    SymTensor(signs, sector, nzblock)

Make a SymTensor with the given possible charges and only one nonzero
sector. Find the dimensions from that.

"""
function SymTensor(signs::NTuple{N, Int},
                   sector::NTuple{N, Int},
                   nzblock::Array{Tv, N}) where {Tv<:Number, N}

    chrs = [0, 1]
    charge = sum(signs .* sector)
    _sector_is_allowed(charge, signs, sector) ||
        error("sector is not allowed ", sector)
    legs = STLeg[]
    legs = Tuple([STLeg(signs[n], [sector[n]], [size(nzblock, n)])
                  for n in eachindex(sector)])
    SymTensor(charge, Tuple(legs), [sector], [nzblock])
end

function SymVector(sign::Int, charge::Int, nzblk::Vector{Tv}) where{Tv<:Number}
    leg = STLeg(sign, [charge], [length(nzblk)])
    SymTensor(charge, (leg,), [(charge,)], [nzblk])
end

function _trimlegs(legs::NTuple{N, STLeg}, sects::Vector{NTuple{N, Int}}) where{N}
    legcharges = Tuple(BitSet() for i=1:N)
    for sect in sects
        for i in 1:N
            push!(legcharges[i], sect[i])
        end
    end
    #print(legcharges)
    trimmedlegs = STLeg[]
    for lidx in 1:N
        leg = legs[lidx]
        indexes = findindexes_sorted(leg.chrs, collect(legcharges[lidx]))
        push!(trimmedlegs, STLeg(leg.sign, leg.chrs[indexes], leg.dims[indexes]))
    end
    # println("trimming :")
    # println(legs)
    # println(trimmedlegs)
    Tuple(trimmedlegs)
end

function convert(::Type{SymTensor{ComplexF64, N}},
                 A::SymTensor{Float64, N}) where {N}
    SymTensor(A.charge, A.legs, A.sects, convert.(Array{ComplexF64, N}, A.nzblks))
end

"""
    findindexes_sorted(A, elements)

Find the index of a set of sorted elements in a sorted vector `A`.
"""
function findindexes_sorted(A::Vector{T}, elements::Vector{T}) where{T}
    #println(A, elements)
    p = 1
    indexes = Int[]
    for i in eachindex(A)
        if A[i] == elements[p]
            push!(indexes, i)
            if p==length(elements) break end
            p += 1
        elseif A[i] < elements[p]
            error("findindexes_sorted ", A, elements)
        end
    end
    indexes
end

function change_legsign(sten::SymTensor{Tv, N}, l::Int) where{Tv<:Number, N}
    legs = (sten.legs[1:l-1]..., change_sign(sten.legs[l]), sten.legs[l+1:N]...)
    sects = NTuple{N, Int}[]
    for sect in sten.sects
        push!(sects, (sect[1:l-1]..., -sect[l], sect[l+1:N]...))
    end
    perm = _sectors_sortperm(sects)
    SymTensor(sten.charge, legs, sects[perm], sten.nzblks[perm])
end

"""
    invlegs(sten)

invert the direction (sign) of all legs of the tensor therefore
negative the total charge of the tensor as well!
"""
function invlegs(sten::SymTensor{Tv, N}) where{Tv<:Number, N}
    legs = Tuple([STLeg(-l.sign, l.chrs, l.dims) for l in sten.legs])
    SymTensor(-sten.charge, legs, sten.sects, sten.nzblks)
end

"""
    combpatterns(charge, legs)

find how many different charge patterns for legs give rise to the given charge.
"""
function _possible_fuse_patterns(charge::Int, legs::NTuple{N, STLeg}) where{N}
    pats = NTuple{N, Int}[]
    patdims = NTuple{N, Int}[]
    if length(legs) > 1
        s = legs[N].sign
        for cidx in eachindex(legs[N].chrs)
            c = legs[N].chrs[cidx]
            d = legs[N].dims[cidx]
            headchr, headdim = _possible_fuse_patterns(charge - s * c, legs[1:N-1])
            for i in eachindex(headchr)
                push!(pats, (headchr[i]..., c))
                push!(patdims, (headdim[i]..., d))
            end
        end
    else
        i = searchsortedfirst(legs[1].chrs, legs[1].sign*charge)
        if i <= length(legs[1].chrs) && legs[1].sign*legs[1].chrs[i] == charge
            return (legs[1].chrs[i]) , (legs[1].dims[i])
        end
    end
    pats, patdims
end

function randSymTensor(Tv::Type, charge::Int, legs::NTuple{N, STLeg}) where{N}
    sects = NTuple{N, Int}[]
    nzblks = Array{Tv, N}[]
    pats, patdims = _possible_fuse_patterns(charge, legs)
    perm = _sectors_sortperm(pats)
    pats, patdims = pats[perm], patdims[perm]
    for patidx in eachindex(pats)
        push!(sects, pats[patidx])
        push!(nzblks, rand(Tv, patdims[patidx]...))
    end
    SymTensor(charge, legs, sects, nzblks)
end

function fill(x::Tv,
              charge::Int,
              legs::NTuple{N, STLeg}) where {Tv<:Number, N}
    sects = NTuple{N, Int}[]
    nzblks = Array{Tv, N}[]
    pats, patdims = _possible_fuse_patterns(charge, legs)
    perm = _sectors_sortperm(pats)
    pats, patdims = pats[perm], patdims[perm]
    for patidx in eachindex(pats)
        push!(sects, pats[patidx])
        push!(nzblks, fill(x, patdims[patidx]...))
    end
    SymTensor(charge, legs, sects, nzblks)
end

fillSymTensor(x, charge, legs) = fill(x, charge, legs)

function fill!(A::SymTensor{Tv, N}, x::Tv) where {Tv<:Number, N}
    for i in eachindex(A.nzblks)
        fill!(A.nzblks[i], x)
    end
    A
end

function eye(Tv::Type,
             charge::Int,
             chrs::Vector{Int},
             dims::Vector{Int})
    l1 = STLeg(+1, chrs, dims)
    l2 = STLeg(-1, chrs.-charge, dims)
    sects = Tuple{Int, Int}[]
    nzblks = Matrix{Tv}[]
    pats, patdims = _possible_fuse_patterns(charge, (l1,l2))
    perm = _sectors_sortperm(pats)
    pats, patdims = pats[perm], patdims[perm]
    for patidx in eachindex(pats)
        push!(sects, pats[patidx])
        push!(nzblks, Matrix{Tv}(I, patdims[patidx]...))
    end
    SymTensor(charge, (l1,l2), sects, nzblks)
end

function isequal(sten1::SymTensor{Tv, N}, sten2::SymTensor{Tv, N}) where{Tv<:Number, N}
    sten1.charge == sten2.charge && sten1.legs == sten2.legs &&
        sten1.sects == sten2.sects && sten1.nzblks == sten2.nzblks
end

==(sten1::SymTensor{Tv, N}, sten2::SymTensor{Tv, N}) where{Tv<:Number, N} = isequal(sten1, sten2)

@inline _sector_is_allowed(total::Int, signs::NTuple{N}, sector::NTuple{N, Int}) where {N} =
    sum(signs .* sector) == total

# @inline size(sten::SymTensor) =
#     Tuple(sum(sten.legs[n].dims) for n in eachindex(sten.legs))


# we always sort sectors based on charges (ignore signs)
### NOTE the sorting is column major. That means the leftmost (first)
### leg changes charge faster (first)!
function _sector_less_than(sect1::NTuple{N, Int},
                           sect2::NTuple{N, Int}
                           #,signs::NTuple{N, Int}
                           ) where {N, Tc<:Integer}
    for n in reverse(eachindex(sect1))
        #sign = signs[n]
        #a, b = sign*sect1[n], sign*sect2[n]
        a, b = sect1[n], sect2[n]
        if a < b
            return true
        elseif a > b
            return false
        end
    end
    false
end

@inline _sectors_sortperm(sects::Vector{NTuple{N, Int}};
                          by::Function=identity) where {N} =
                              sortperm(sects, by=by, lt=_sector_less_than)

function conj(sten::SymTensor{Tv}) where {Tv<:Number}
    SymTensor(sten.charge, sten.legs, sten.sects,
              [conj(blk) for blk in sten.nzblks])
end

function change_nzblk!(sten::SymTensor{Tv, N},
                       sect::NTuple{N, Int},
                       nzblk::Array{Tv, N}) where {Tv<:Number, N}

    # should I use a different search fn?
    index, = searchsorted(sten.sects, sect, lt=_sector_less_than)
    length(index) == 0 && error("sector not found!")

    size(nzblk) != size(sten.nzblks[index]) &&
        error("non-zero block size doesn't match",
              size(nzblk), " == ", size(sten.nzblks[index]))

    sten.nzblks[index] = nzblk
    nothing
end

function mapcharges(f::Function, A::SymTensor{Tv, N}) where{Tv<:Number, N}
    SymTensor(f(A.charge),
              mapcharges.(f, A.legs),
              [f.(s) for s in A.sects],
              A.nzblks)
end

function mapcharges(f::NTuple{N, Function}, A::SymTensor{Tv, N}) where{Tv<:Number, N}
    legs = Tuple(mapcharges(f[i], A.legs[i]) for i in 1:N)
    sects = [Tuple(f[i](s[i]) for i in 1:N) for s in A.sects]
    sgns = signs(legs)
    charges = [sum(sgns .* s) for s in sects]
    !all(charges .== charges[1]) && error("mapcharges function is inconsistent!")
    SymTensor(charges[1], legs, sects, A.nzblks)
end

function array_representation(sten::SymTensor{Tv}) where {Tv<:Number}
    arrep = zeros(sum.(alldims(sten.legs)))
    adims = Tuple([0;s] for s in accdims(sten.legs))
    #println(adims)
    chrs = [0, 1]
    for idx in eachindex(sten.sects)
        sect = sten.sects[idx]
        nzblk = sten.nzblks[idx]
        ranges = []
        for n in eachindex(sect)
            i = findfirst(sect[n].== sten.legs[n].chrs)
            push!(ranges, adims[n][i]+1:adims[n][i+1])
        end
        arrep[Tuple(ranges)...] = nzblk
    end
    arrep
end

function show(ten::SymTensor)
    num_ch = length(ten.legs)
    #m, n = size(bsa)
    println(size(ten), typeof(ten), " with ", length(ten.sects), " blocks: ")
    for i in eachindex(ten.sects)
        print("sector: ", ten.sects[i], ", dims: ")
        show(stdout, "text/plain", ten.nzblks[i])
        println()
    end
    nothing
end

function *(A::SymTensor, a::T) where {T<:Number}
    SymTensor(A.charge, A.legs, A.sects, [a .* blk for blk in A.nzblks])
end

function removedummyleg(sten::SymTensor{Tv, N}, l::Int) where {Tv,N}
    !isdummy(sten.legs[l]) && error("leg is not dummy!")

    SymTensor(sten.charge,
              (sten.legs[1:l-1]..., sten.legs[l+1:N]...),
              [(s[1:l-1]...,s[l+1:N]...) for s in sten.sects],
              [reshape(blk, size(blk)[1:l-1]...,size(blk)[l+1:N]...) for blk in sten.nzblks])

end

"the scalar type (i.e. <:Number) of the data in v"
eltype(v::SymTensor{Tv}) where {Tv<:Number} = Tv

"a way to construct additional similar vectors, possibly with a
different scalar type T."
function similar(v::SymTensor{Tv, N}, T::Type=Tv) where {Tv<:Number, N}
    nzblks = [similar(blk) for blk in v.nzblks]
    SymTensor(v.charge, v.legs, v.sects, nzblks)
end

"copy the contents of v to a preallocated vector w"
function copyto!(w::SymTensor, v::SymTensor)
    w.nzblks = v.nzblks
    w
end

"out of place scalar multiplication; multiply vector v with scalar α
and store the result in w"
function mul!(w::SymTensor, v::SymTensor, α)
    w.nzblks = [α .* blk for blk in v.nzblks]
    w
end

"in-place scalar multiplication of v with α; in particular with α =
false, v is initialized with all zeros"
function rmul!(v::SymTensor, α)
    v.nzblks .*= α
    v
end

####TODO: I probably have to force allocation of all possbile sectors!
####(in order to prevent recreation of the symtensor) but this
####requires a bit of thinking and benchmarking to see which options works best.

###NOTE: We assume that w has more sectors than v!  because all
### possible sectors may not be used! (should I force this to happen!)
" store in w the result of α*v + w"
function axpy!(α, v::SymTensor{T1, N}, w::SymTensor{T2, N}) where{T1,T2,N}
    sects = NTuple{N, Int}[]
    nzblks = Array{T2, N}[]

    nv = length(v.sects)
    nw = length(w.sects)
    iv, iw = 1, 1
    while (iv <= nv || iw <=nw)
        if iv <= nv && iw <=nw && w.sects[iw] == v.sects[iv]
            push!(sects, w.sects[iw])
            push!(nzblks, α .* v.nzblks[iv] + w.nzblks[iw])
            iv+=1
            iw+=1
        elseif iv > nw || _sector_less_than(w.sects[iw], v.sects[iv])
            push!(sects, w.sects[iw])
            push!(nzblks, w.nzblks[iw])
            iw+=1
        else
            push!(sects, v.sects[iv])
            push!(nzblks, α .* v.nzblks[iv])
            iv+=1
        end
    end
    w = SymTensor(w.charge, w.legs, sects, nzblks)
    w
end
###NOTE: We assume that w has more sectors than v!  because all
### possible sectors may not be used! (should I force this to happen!)
" store in w the result of α*v + β*w"
function axpby!(α, v::SymTensor, β, w::SymTensor)
    sects = NTuple{N, Int}[]
    nzblks = Array{T2, N}[]

    nv = length(v.sects)
    nw = length(w.sects)
    iv, iw = 1, 1
    while (iv <= nv || iw <=nw)
        if iv <= nv && iw <=nw && w.sects[iw] == v.sects[iv]
            push!(sects, w.sects[iw])
            push!(nzblks, α .* v.nzblks[iv] + β .* w.nzblks[iw])
            iv+=1
            iw+=1
        elseif iv > nw || _sector_less_than(w.sects[iw], v.sects[iv])
            push!(sects, w.sects[iw])
            push!(nzblks, β .* w.nzblks[iw])
            iw+=1
        else
            push!(sects, v.sects[iv])
            push!(nzblks, α .* v.nzblks[iv])
            iv+=1
        end
    end
    w = SymTensor(w.charge, w.legs, sects, nzblks)
    w
end

" compute the inner product of two vectors"
function dot(v::SymTensor{Tv, N}, w::SymTensor{Tv, N}) where{Tv<:Number,N}
    contract(v, .-Tuple(1:N), invlegs(conj(w)), .-Tuple(1:N))
end

" compute the 2-norm of a vector"
function norm(v::SymTensor)
    sqrt(dot(v,v))
end

function makenormalize!(S::SymTensor)
    #s = norm(vcat([diag(blk) for blk in S.nzblks]...))
    s = norm(S)
    for i in eachindex(S.nzblks)
        S.nzblks[i] = S.nzblks[i] ./ s
    end
    S
end
