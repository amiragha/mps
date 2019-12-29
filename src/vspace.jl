abstract type AbstractVectorSpace end
struct VectorSpace{S} <: AbstractVectorSpace
    sectors :: SortedDict{S, Int}

    function VectorSpace{S}(s) where {S}
        allunique([c for (c,d) in s]) ||
            error("duplicate charges! $s")
        sectors = SortedDict{S, Int}(s)
        for (c, d) in sectors
            d > 0 || error("space dimension can only be positive!")
        end
        new{S}(sectors)
    end

    #VectorSpace{S}() where{S} = VectorSpace{S}(SortedDict{S, Int}())
end

VectorSpace{S}(sectors...) where {S} = VectorSpace{S}(sectors)

@inline vtype(::VectorSpace{S}) where {S} = S

@inline isequal(V1::VectorSpace, V2::VectorSpace) = isequal(V1.sectors, V2.sectors)
==(V1::VectorSpace, V2::VectorSpace) = isequal(V1, V2)

@inline hascharge(V::VectorSpace{S}, charge::S) where{S}= haskey(V.sectors)
@inline findcharge(V::VectorSpace{S}, charge::S) where{S}= haskey(V.sectors)
@inline charges(V::VectorSpace) = [c for (c,d) in V.sectors]
@inline dims(V::VectorSpace) = [d for (c,d) in V.sectors]

"Find the dimension of the the given `charge`. Returns 0 if charges
doesn't exist. Return total dimension if charge is not specified."
@inline dim(V::VectorSpace{S}, charge::S) where{S}= get(V.sectors, charge, 0)
@inline dim(V::VectorSpace) = sum(dims(V))

@inline Base.length(V::VectorSpace) = length(V.sectors)
@inline Base.iterate(V::VectorSpace) = iterate(V.sectors)
@inline Base.iterate(V::VectorSpace, i) = iterate(V.sectors, i)

"Return the dual of the vector space"
@inline dual(V::VectorSpace) = VectorSpace{vtype(V)}(-c=>d for (c,d) in V.sectors)

@inline isdual(V1::T, V2::T) where {T<:VectorSpace} =
    charges(V1) == -reverse(charges(V2)) && dims(V1) == reverse(dims(V2))

function layout(V::VectorSpace; rev::Bool=false)
    S = vtype(V)
    dict = SortedDict{S, UnitRange{Int}}()
    if rev
        p = dim(V)
        for (c, d) in V
            dict[c] = p-d+1:p
            p -= d
        end
    else
        p = 0
        for (c,d) in V
            dict[c] = p+1:p+d
            p+=d
        end
    end
    dict
end

const U1Space = VectorSpace{Int}

"Find sector intersections between two VSpace"
function intersect(V1::VectorSpace{S}, V2::VectorSpace{S}) where {S}
    sectors = SortedDict{S, Int}()
    for (c, d) in V1.sectors
        if d == get(V2, c, 0)
            sectors[c] = d
        elseif d > 0
            error("same charge but different dims!")
        end
    end
    return VectorSpace(sectors)
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

fuse(Vs::VectorSpace...) = fuse(Vs)

function fuse(Vs::NTuple{N, VectorSpace{S}}) where {N, S}
    if N < 2
        return Vs[1]
    elseif N > 2
        # lets for now generally always fuse from right to left!
        return fuse(Vs[1:N-2]..., fuse(Vs[N-1:N]...))
    end
    # N == 2 case
    sectors = SortedDict{S, Int}()
    for (c2, d2) in Vs[2]
        for (c1, d1) in Vs[1]
            if haskey(sectors, c1+c2)
                sectors[c1+c2] += d1*d2
            else
                sectors[c1+c2] = d1*d2
            end
        end
    end
    VectorSpace{S}(sectors)
end

"map the charges of the VSpace by some strictly ascending function `f`"
function mapcharges(f::Function, V::VectorSpace)
    mappedchrs = f.(charges(V))
    issoretd(mappedchrs) ||
        error("function not strictly ascending!")
    allunique(mappedchrs) ||
        error("mapping charges creates duplicates! $(l.chrs) -> $mappedchrs")
    VectorSpace(f(c)=>d for (c,d) in V.sectors)
end

isdummy(V::VectorSpace) = charges(V) == (0,) && V[0] == 1
dummyleg() = VectorSpace(0 => 1)

Base.show(io::IO, ::Type{U1Space}) = print(io, "U1Space")

function Base.show(io::IO, V::VectorSpace)
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
