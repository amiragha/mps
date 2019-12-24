abstract type RepVectorSpace end
struct U1Space <: RepVectorSpace
    sectors :: SortedDict{Int, Int}

    function U1Space(s)
        allunique([c for (c,d) in s]) ||
            error("duplicate charges! $s")
        sectors = SortedDict{Int, Int}(s)
        for (c, d) in sectors
            d > 0 || error("space dimension can only be positive!")
        end
        new(sectors)
    end
end
U1Space(sectors...) = U1Space(sectors)

@inline isequal(V1::U1Space, V2::U1Space) = isequal(V1.sectors, V2.sectors)
==(V1::U1Space, V2::U1Space) = isequal(V1, V2)

@inline hascharge(V::U1Space, charge::Int) = haskey(V.sectors)

"Find the dimension of the the given `charge`. Returns 0 if charges
doesn't exist. Return total dimension if charge is not specified."
@inline dim(V::U1Space, charge::Int) = get(V.sectors, charge, 0)
@inline dim(V::U1Space) = sum(chargedims(V))

@inline charges(V::U1Space) = (c for (c,d) in V.sectors)
@inline chargedims(V::U1Space) = (d for (c,d) in V.sectors)

@inline Base.iterate(V::U1Space) = iterate(V.sectors)
@inline Base.iterate(V::U1Space, i) = iterate(V.sectors, i)


#@inline charges(legs::NTuple{N, U1Space}) where {N} = charges.(legs)

#@inline signs(legs::NTuple{N, U1Space}) where {N} =
#    Tuple(legs[n].sign for n in eachindex(legs))

#@inline alldims(legs::NTuple{N, U1Space}) where {N} =
#    Tuple(legs[n].dims for n in eachindex(legs))

#@inline accdims(legs::NTuple{N, U1Space}) where {N} =
#    Tuple(cumsum(legs[n].dims) for n in eachindex(legs))

#@inline fulldims(V::U1Space) = sum(d for (c,d) in V.sectors)
#@inline fulldims(legs::NTuple{N, U1Space}) where {N} = prod(fulldims.legs)

"Return the dual of the vector space"
@inline dual(V::U1Space) = U1Space(-c=>d for (c,d) in V.sectors)

"Find sector intersections between two VSpace"
function intersect(V1::U1Space, V2::U1Space)
    sectors = SortedDict{Int, Int}()
    for (c, d) in V1.sectors
        if d == get(V2, c, 0)
            sectors[c] = d
        elseif d > 0
            error("same charge but different dims!")
        end
    end
    return U1Space(sectors)
end

"""
    fuse(Vs)

fuse any number of VSpaces into the tensor product VSpace. useful for
fuselegs fn. Note that fuse is a binary operation. That is the tensor
product is the bifunctor in the monoidal category of vector spaces. So
it should be by definition a recursive function. And the fusion tree
(order/style of recursion) has to specified. Here I assume a right to
left fusion, (V1 ⊗ (V2 ⊗ (V3 ⊗ ...))).

The order matters because the associator is not identity for
representation vector spaces as the catogory of tensor product of
representation spaces is not an strict monoidal cateogry. So one needs
the associator isomorphism to change the fusion tree.

"""

fuse(Vs::U1Space...) = fuse(Vs)

function fuse(Vs::NTuple{N, U1Space}) where {N}
    if N < 2
        return Vs[1]
    elseif N > 2
        # lets for now generally always fuse from right to left!
        return fuse(Vs[1:N-2]..., fuse(Vs[N-1:N]...))
    end
    # N == 2 case
    sectors = SortedDict{Int, Int}()
    for (c2, d2) in Vs[2]
        for (c1, d1) in Vs[1]
            if haskey(sectors, c1+c2)
                sectors[c1+c2] += d1*d2
            else
                sectors[c1+c2] = d1*d2
            end
        end
    end
    U1Space(sectors)
end

"map the charges of the VSpace by some strictly ascending function `f`"
function mapcharges(f::Function, V::U1Space)
    mappedchrs = f.(charges(V))
    issoretd(mappedchrs) ||
        error("function not strictly ascending!")
    allunique(mappedchrs) ||
        error("mapping charges creates duplicates! $(l.chrs) -> $mappedchrs")
    U1Space(f(c)=>d for (c,d) in V.sectors)
end

isdummy(V::U1Space) = charges(V) == (0,) && V[0] == 1
dummyleg() = U1Space(0 => 1)

Base.show(io::IO, ::Type{U1Space}) = print(io, "U1Space")

function Base.show(io::IO, V::RepVectorSpace)
    show(io, typeof(V))
    print(io, "(")
    seperator = ""
    comma = ", "
    io2 = IOContext(io, :typeinfo => typeof(V))
    for c in charges(V)
        print(io2, seperator, c, "=>", dim(V, c))
        seperator = comma
    end
    print(io, ")")
    return nothing
end
