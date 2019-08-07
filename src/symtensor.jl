mutable struct FusedCharge{N}
    charge :: Int # total charge
    dim :: Int # total dimension
    #signs :: NTuple{N, Int}
    pats :: Dict{NTuple{N, Int}, UnitRange{Int}} # patterns and positions

    function FusedCharge(charge, dim, pats::Dict{NTuple{N, Int}, UnitRange{Int}}) where {N}
        all(sum.(keys(pats)) .== charge) || error("oops", pats, charge)
        #all(abs.(signs) .== 1) || error("oops")
        new{N}(charge, dim, pats)
    end
end

# struct DefusedCharge{N}
#     charge :: Int
#     dim :: Int

# end
struct STLeg
    sign :: Int    # +1 for ket (ingoing) or -1 for bra (outgoing)
    chrs :: Vector{Int}    # set of possible charges
    dims :: Vector{Int}    # set of dimensions for possible charges

    function STLeg(sign, chrs, dims)
        abs(sign) == 1 || error("unacceptable sign : ", sign)
        all(dims .> 0) || error("space dimensions can only be positive!")
        length(chrs) == length(dims) ||
            error("Number of charges and dimensions don't match!")
        issorted(chrs) || error("charges are not sorted")
        new(sign, chrs, dims)
    end
end

function isequal(l1::STLeg, l2::STLeg)
    l1.sign == l2.sign && l1.chrs == l2.chrs && l1.dims == l2.dims
end

==(l1::STLeg, l2::STLeg) = isequal(l1, l2)

# This function assumes that charge definitely exists!
function getdim(leg::STLeg, charge::Int)
    leg.dims[searchsortedfirst(leg.chrs, charge)]
end

@inline charges(legs::NTuple{N, STLeg}) where {N} =
    Tuple(legs[n].chrs for n in eachindex(legs))

@inline signs(legs::NTuple{N, STLeg}) where {N} =
    Tuple(legs[n].sign for n in eachindex(legs))

@inline alldims(legs::NTuple{N, STLeg}) where {N} =
    Tuple(legs[n].dims for n in eachindex(legs))

@inline accdims(legs::NTuple{N, STLeg}) where {N} =
    Tuple(cumsum(legs[n].dims) for n in eachindex(legs))

@inline fulldims(leg::STLeg) = sum(leg.dims)
@inline fulldims(legs::NTuple{N, STLeg}) where {N} =
    prod([fulldims(legs[n]) for n in eachindex(legs)])

function change_sign(leg::STLeg)
    perm = sortperm(leg.chrs, by=x -> -x)
    STLeg(-leg.sign, -1 .*(leg.chrs[perm]), leg.dims[perm])
end

struct SymTensor{Tv<:Number, N}
    charge :: Int     # total charge of the tensor
    legs :: NTuple{N, STLeg}
    sects :: Vector{NTuple{N, Int}}    # possible nonzero sectors
    nzblks :: Vector{AbstractArray{Tv, N}}     # stored nonzero blocks

    function SymTensor(charge::Int,
                       legs::NTuple{N, STLeg},
                       sects::Vector{NTuple{N, Int}},
                       nzblks::Vector{<:AbstractArray{Tv, N}}) where{Tv<:Number, N}
        issorted(sects, lt=_sector_less_than) ||  error("sectors not sorted!")
        legs = _trimlegs(legs, sects)
        new{Tv, N}(charge, legs, sects, nzblks)
    end
end

const SymMatrix{Tv} = SymTensor{Tv, 2} where {Tv<:Number}
@inline numoflegs(s::SymTensor{Tv, N}) where{Tv<:Number, N} = N

#struct SymDiagonal{Tv} =
"""
    SymTensor(signs, sector, nzblock)

Make a SymTensor with the given possible charges and only one nonzero
sector. Find the dimensions from that.

"""
function SymTensor(signs::NTuple{N, Int},
                   sector::NTuple{N, Int},
                   nzblock::Array{Tv, N},
                   charge::Int=0) where {Tv<:Number, N}

    chrs = [0, 1]
    _sector_is_allowed(charge, signs, sector) ||
        error("sector is not allowed ", sector)
    legs = STLeg[]
    legs = Tuple([STLeg(signs[n], sector[n], size(nzblock, n))
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

function fillSymTensor(value::Tv, charge::Int, legs::NTuple{N, STLeg}) where{Tv<:Number, N}
    sects = NTuple{N, Int}[]
    nzblks = Array{Tv, N}[]
    pats, patdims = _possible_fuse_patterns(charge, legs)
    perm = _sectors_sortperm(pats)
    pats, patdims = pats[perm], patdims[perm]
    for patidx in eachindex(pats)
        push!(sects, pats[patidx])
        push!(nzblks, fill(value, patdims[patidx]...))
    end
    SymTensor(charge, legs, sects, nzblks)
end

function isequal(sten1::SymTensor{Tv, N}, sten2::SymTensor{Tv, N}) where{Tv<:Number, N}
    sten1.charge == sten2.charge && sten1.legs == sten2.legs &&
        sten1.sects == sten2.sects && sten1.nzblks == sten2.nzblks
end

==(sten1::SymTensor{Tv, N}, sten2::SymTensor{Tv, N}) where{Tv<:Number, N} = isequal(sten1, sten2)

@inline _sector_is_allowed(total::Int, signs::NTuple{N}, sector::NTuple{N, Int}) where {N} =
    sum(signs .* sector) == total

@inline size(sten::SymTensor) =
    Tuple(sum(sten.legs[n].dims) for n in eachindex(sten.legs))


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


################################
#### RELEG #####################
################################

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

function permutelegs(sten::SymTensor{T, N}, perm) where{T, N}
    @assert sort(perm) == collect(1:N)
    sectperm = _sectors_sortperm(sten.sects,  by=x->x[perm])
    nzblks = Array{T, N}[]
    sects = [sect[perm] for sect in sten.sects[sectperm]]
    for nzblk in sten.nzblks[sectperm]
        #println(size(nzblk), " ", perm)
        push!(nzblks, permutedims(nzblk, perm))
    end
    SymTensor(sten.charge, sten.legs[perm], sects, nzblks)
end
