abstract type AbstractVectorSpace end
struct TrivialVectorSpace <:AbstractVectorSpace
    n :: Int
end
struct VectorSpace{S} <: AbstractVectorSpace
    sectors :: SortedDict{S, Int}
    isdual  :: Bool
    function VectorSpace{S}(s, id) where {S}
        allunique([c for (c,d) in s]) ||
            error("duplicate charges! $s")
        sectors = SortedDict{S, Int}(s)
        for (c, d) in sectors
            d > 0 || error("space dimension can only be positive!")
        end
        new{S}(sectors, id)
    end
    VectorSpace{S}(s) where{S} = VectorSpace{S}(s, false)
end
VectorSpace{S}(pairs::Pair{S,Int}...) where {S} = VectorSpace{S}(pairs, false)

@inline vtype(::VectorSpace{S}) where {S} = S

@inline isdual(V::VectorSpace) = V.isdual
@inline isdual(V1::T, V2::T) where {T<:VectorSpace} =
    V1.sectors == V2.sectors && V1.isdual != V2.isdual
@inline isequal(V1::VectorSpace, V2::VectorSpace) =
    isequal(V1.sectors, V2.sectors) && V1.isdual == V2.isdual
==(V1::VectorSpace, V2::VectorSpace) = isequal(V1, V2)

@inline hascharge(V::VectorSpace{S}, charge::S) where{S}= haskey(V.sectors)
@inline findcharge(V::VectorSpace{S}, charge::S) where{S}= haskey(V.sectors)
@inline charges(V::VectorSpace) = [c for (c,d) in V.sectors]
@inline dims(V::VectorSpace) = [d for (c,d) in V.sectors]

"Find the dimension of the the given `charge`. Returns 0 if charges
doesn't exist. Return total dimension if charge is not specified."
@inline dim(V::VectorSpace{S}, charge::S) where{S}= get(V.sectors, charge, 0)
@inline dim(V::VectorSpace{S}, charge::Int) where{S}= get(V.sectors, convert(S, charge), 0)
@inline dim(V::VectorSpace) = sum(dims(V))

@inline Base.length(V::VectorSpace) = length(V.sectors)
@inline Base.iterate(V::VectorSpace) = iterate(V.sectors)
@inline Base.iterate(V::VectorSpace, i) = iterate(V.sectors, i)

"Return the dual of the vector space"
@inline dual(V::VectorSpace) = VectorSpace{vtype(V)}(V.sectors, !V.isdual)

function layout(V::VectorSpace)
    S = vtype(V)
    dict = SortedDict{S, UnitRange{Int}}()
    # if V.isdual
    #     p = dim(V)
    #     for (c, d) in V
    #         dict[c] = p-d+1:p
    #         p -= d
    #     end
    # else
        p = 0
        for (c,d) in V
            dict[c] = p+1:p+d
            p+=d
        end
    #end
    dict
end

const U1Space = VectorSpace{U1Charge}
U1Space(pairs::Pair{Int, Int}...) = U1Space(U1(c)=>d for (c,d) in pairs)

"Find sector intersections between two VSpace"
function intersect(V1::VectorSpace{S}, V2::VectorSpace{S}) where {S}
    V1.isdual != V2.isdual || error("Not the same dual-type")
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

fuse(_dual::Bool, Vs::VectorSpace...) = fuse(_dual, Vs)
#fuse(Vs::VectorSpace...) = fuse(Vs)
#fuse(Vs::NTuple{VectorSpace{S}}) where{N,S} = fuse(false, Vs)
function fuse(_dual::Bool, Vs::NTuple{N, VectorSpace{S}}) where {N, S}
    if N < 2
        if _dual == Vs[1].isdual
            return Vs[1]
        else
            return dual(mapcharges(x->inv(x), Vs[1]))
        end
    elseif N > 2
        # lets for now generally always fuse from right to left!
        return fuse(_dual, Vs[1:N-2]..., fuse(_dual, Vs[N-1:N]...))
    end
    # N == 2 case
    sectors = SortedDict{S, Int}()
    for (c2, d2) in Vs[2]
        c2 = Vs[2].isdual ? inv(c2) : c2
        for (c1, d1) in Vs[1]
            c1 = Vs[1].isdual ? inv(c1) : c1
            c = _dual ? inv(c1+c2) : c1+c2
            if haskey(sectors, c)
                sectors[c] += d1*d2
            else
                sectors[c] = d1*d2
            end
        end
    end
    VectorSpace{S}(sectors, _dual)
end

"map the charges of the VSpace by some strictly ascending function `f`"
function mapcharges(f::Function, V::VectorSpace)
    mappedchrs = f.(charges(V))
    # issorted(mappedchrs) ||
    #     error("function not strictly ascending!")
    allunique(mappedchrs) ||
        error("mapping charges creates duplicates! $(l.chrs) -> $mappedchrs")
    VectorSpace{vtype(V)}(SortedDict(f(c)=>d for (c,d) in V.sectors), V.isdual)
end

isdummy(V::VectorSpace) = charges(V) == (zero(S),) && V[zero(S)] == 1
dummyleg() = VectorSpace{S}(0 => 1)

Base.show(io::IO, ::Type{U1Space}) = print(io, "U1Space")

function Base.show(io::IO, V::VectorSpace)
    show(io, typeof(V))
    V.isdual && print(io, "†")
    print(io, "(")
    seperator = ""
    comma = ", "
    io2 = IOContext(io, :typeinfo => typeof(V))
    for c in charges(V)
        if V.isdual
            print(io2, seperator, c, "'", "=>", dim(V, c))
        else
            print(io2, seperator, c, "=>", dim(V, c))
        end
        seperator = comma
    end
    print(io, ")")
    return nothing
end
