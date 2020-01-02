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

@inline Base.length(::Sector{S,N}) where {S,N} = N
@inline function Base.sum(s::Sector{S, N}, dualinfo::NTuple{N, Bool}) where{S,N}
    r = zero(S)
    for i=1:N
        r += dualinfo[i] ? inv(s[i]) : s[i]
    end
    r
end
@inline Base.inv(s::Sector{S}) where{S} = Sector{S}([inv(c) for c in s]...)
Base.getindex(s::Sector, i) = getindex(s.charges, i)
Base.iterate(s::Sector) = iterate(s.charges)
Base.iterate(s::Sector, i) = iterate(s.charges, i)

@inline layout(space::NTuple{N, VectorSpace{S}},
               s::Sector{S, N}) where {S, N} =
                   Tuple(layout(space[i])[s[i]] for i=1:N)

"""
we always sort sectors based on charges The sorting is column
major. That means the leftmost (first) VSpace changes charge faster
(first)!
"""
@inline Base.isless(s1::Sector{S, N}, s2::Sector{S, N}) where {N, S} =
    reverse(s1.charges) < reverse(s2.charges)

# @inline function Base.isless(s1::Sector{S, N},
#                              s2::Sector{S, N}) where {N, S}
#     for n in N:-1:1
#         @inbounds a, b = s1[n], s2[n]
#         if a < b
#             return true
#         elseif a > b
#             return false
#         end
#     end
#     false
# end

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
    if N < 2
        _c = Vs[1].isdual ? inv(charge) : charge
        d = dim(Vs[1], _c)
        if d > 0
            return [Sector(_c)], [(d,)]
        end
        return Sector{S, 1}[], NTuple{1, Int}[]
    end

    sects = Sector{S, N}[]
    sizes = NTuple{N, Int}[]
    for (c,d) in Vs[N]
        _c = Vs[N].isdual ? inv(c) : c
        sect_head, size_head = _allsectorsandsizes(charge - _c, Vs[1:N-1])
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
