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
@inline getdim(leg::STLeg, charge::Int) =
    leg.dims[searchsortedfirst(leg.chrs, charge)]

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

@inline reverseleg(leg::STLeg) = STLeg(-leg.sign, leg.chrs, leg.dims)

function negate(leg::STLeg)
    perm = sortperm(leg.chrs, by=x -> -x)
    STLeg(-leg.sign, -1 .*(leg.chrs[perm]), leg.dims[perm])
end

function intersect(l1::STLeg, l2::STLeg)
    chrs = Int[]
    idxs = (Int[], Int[])
    for i in eachindex(l1.chrs)
        c = l1.chrs[i]
        ix = searchsortedfirst(l2.chrs, c)
        if ix <= length(l2.chrs) && l2.chrs[ix] == c
            push!(chrs, c)
            push!(idxs[1], i)
            push!(idxs[2], ix)
        end
    end
    return chrs, idxs...
end

function fuse(sign::Int, legs::NTuple{N, STLeg}) where {N}

    #signs = sign .* [leg.sign for leg in legs]
    lchrs = sign .* [leg.sign .* leg.chrs for leg in legs]
    ldims = [leg.dims for leg in legs]
    lengs = [length(leg.chrs) for leg in legs]

    fdims = reshape(prod.(Iterators.product(ldims...)), prod(lengs))
    fchrs = reshape(sum.(Iterators.product(lchrs...)), prod(lengs))

    perm = sortperm(fchrs)
    sortedfchrs = fchrs[perm]
    sortedfdims = fdims[perm]

    chrs = Int[]
    dims = Int[]
    c = undef
    for i in eachindex(fchrs)
        if c != sortedfchrs[i]
            c = sortedfchrs[i]
            push!(chrs, c)
            push!(dims, sortedfdims[i])
        else
            dims[end] += sortedfdims[i]
        end
    end
    STLeg(sign, chrs, dims)
end

function mapcharges(f::Function, l::STLeg)
    STLeg(l.sign, f.(l.chrs), l.dims)
end

isdummy(l::STLeg) = l.chrs == [0] && l.dims == [1]

"""
    _allsectorsandsizes(charge, legs)

Generate all sectors for a tuple of legs that add up to a given
`charge` and are thus sectors are the corresponding SymTensor. This is
recursive function calling itself with legs[1:N-1], etc.

since the legs have sorted charges, the output sectors are sorted
"""
function _allsectorsandsizes(charge::Int, legs::NTuple{N, STLeg}) where{N}
    sects = NTuple{N, Int}[]
    sizes = NTuple{N, Int}[]
    if length(legs) > 1
        s = legs[N].sign
        for cidx in eachindex(legs[N].chrs)
            c = legs[N].chrs[cidx]
            d = legs[N].dims[cidx]
            sect_head, size_head = _allsectorsandsizes(charge - s * c, legs[1:N-1])
            for i in eachindex(sect_head)
                push!(sects, (sect_head[i]..., c))
                push!(sizes, (size_head[i]..., d))
            end
        end
    else
        i = searchsortedfirst(legs[1].chrs, legs[1].sign*charge)
        if i <= length(legs[1].chrs) && legs[1].sign*legs[1].chrs[i] == charge
            return [(legs[1].chrs[i],)] , [(legs[1].dims[i],)]
        end
    end
    sects, sizes
end

# interestingly the below method is singnificantly slower!
function _allsectorsandsizes2(charge::Int, legs::NTuple{N, STLeg}) where{N}
    sects = NTuple{N, Int}[]
    sizes = NTuple{N, Int}[]
    signs = [leg.sign for leg in legs]
    for indeces in Iterators.product([eachindex(legs[i].chrs) for i=1:N]...)
        sector = Tuple(legs[i].chrs[indeces[i]] for i in 1:N)
        if sum(signs .* sector) == charge
            push!(sects, sector)
            push!(sizes, Tuple(legs[i].dims[indeces[i]] for i in 1:N))
        end
    end
    sects, sizes
end

@inline _sectorisallowed(total::Int, signs::NTuple{N}, sector::NTuple{N, Int}) where {N} =
    sum(signs .* sector) == total

##TODO: this should move to the documentation of SymTensors!
# we always sort sectors based on charges (ignore signs)
### NOTE the sorting is column major. That means the leftmost (first)
### leg changes charge faster (first)!

function _sectorlessthan(s1::NTuple{N, Int},
                         s2::NTuple{N, Int}) where {N}
    for n in reverse(eachindex(s1))
        a, b = s1[n], s2[n]
        if a < b
            return true
        elseif a > b
            return false
        end
    end
    false
end
