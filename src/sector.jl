struct Sector{S<:AbstractCharge, N}
    charges :: NTuple{N, S}

    Sector{S,N}(cs) where {S,N} = new{S,N}(cs)
    Sector{S,N}(cs::S...) where {S<:AbstractCharge,N} = new{S,N}(cs)
    function Sector{S}(cs...) where {S}
        new{S,length(cs)}(Tuple(convert(S,c) for c in cs))
    end
    Sector(cs::NTuple{N, S}) where {S, N} = new{S,N}(cs)
    Sector(cs::S...) where {S<:AbstractCharge} = new{S, length(cs)}(cs)
end

@inline Base.sum(s::Sector) = sum(s)
@inline inv(s::Sector) = Sector(inv(c) for c in s)
Base.getindex(s::Sector, i) = getindex(s.charges, i)
Base.iterate(s::Sector) = iterate(s.charges)
Base.iterate(s::Sector, i) = iterate(s.charges, i)

@inline layout(space::NTuple{N, VectorSpace{S}},
               s::Sector{S, N};
               rev::NTuple{N,Bool}=zeros(Bool,N)) where {S, N} =
                   Tuple(layout(space[i], rev=rev[i])[s[i]] for i=1:N)

"""
we always sort sectors based on charges The sorting is column
major. That means the leftmost (first) VSpace changes charge faster
(first)!
"""
@inline Base.isless(s1::Sector{S, N}, s2::Sector{S, N}) where {N, S} =
    reverse(s1.charges) < reverse(s2.charges)

function _check_sectordatadict(charge, space, data)
    sects, sizes = _allsectorsandsizes(charge, space)
    i = 1
    for (s,d) in data
        s == sects[i] || throw{SectorMismatch()}
        size(d) == sizes[i] || throw(SizeMismatch())
        i += 1
    end
end

"""
    _allsectorsandsizes(charge, Vs)

Generate all sectors for a tuple of VSpace that add up to a given
`charge` and are thus sectors are the corresponding SymTensor. This is
recursive function calling itself with Vs[1:N-1], etc. The output
sector is sorted by construction.
"""
function _allsectorsandsizes(charge::S, Vs::NTuple{N, VectorSpace{S}}) where{S, N}
    if length(Vs) < 2
        d = dim(Vs[1], charge)
        if d > 0
            return [Sector(charge)], [(d,)]
        end
        return Sector{S, 1}[], NTuple{1, Int}[]
    end

    sects = Sector{S, N}[]
    sizes = NTuple{N, Int}[]
    for (c,d) in Vs[N]
        sect_head, size_head = _allsectorsandsizes(charge - c, Vs[1:N-1])
        for i in eachindex(sect_head)
            push!(sects, Sector(sect_head[i].charges..., c))
            push!(sizes, (size_head[i]..., d))
        end
    end
    sects, sizes
end
@inline _allsectorsandsizes(c::Int, Vs::NTuple{N, VectorSpace{S}}) where{S,N} =
    _allsectorsandsizes(convert(S, c), Vs)

function Base.show(io::IO, s::Sector)
    if !get(io, :compact, false)
        show(io, typeof(s))
    end
    print(io, "(")
    separator = ""
    comma = ", "
    io2 = IOContext(io, :typeinfo => typeof(s))
    for c in s.charges
        print(io2, separator, c)
        separator=comma
    end
    print(io2, ")")
    return nothing
end
