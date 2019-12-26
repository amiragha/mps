struct Sector{C, N}
    charges :: NTuple{N, C}
end
Sector(s...) = Sector(s)
#const Sector{C, N} = NTuple{N, C}
const U1Sector{N} = Sector{Int, N}

@inline Base.sum(S::Sector) = sum(S)
@inline dual(S::Sector) = Sector(dual(s) for s in S)

"""
we always sort sectors based on charges The sorting is column
major. That means the leftmost (first) VSpace changes charge faster
(first)!
"""
@inline Base.isless(s1::Sector{C, N}, s2::Sector{C, N}) where {N, C} =
 reverse(s1.charges) < reverse(s2.charges)

const SymSectorData{S, N, T} = SortedDict{Sector{S, N}, T}

function _check_sectordatadict(charge, space, data)
    sects, sizes = _allsectorsandsizes(charge, space)
    i = 1
    for (s,d) in data
        s == sects[i] || throw{SectorMismatch()}
        size(d) == sizes[i] || throw(SizeMismatch())
        i += 1
    end
end


abstract type SectorData{N, T} end
mutable struct SectorArray{N, T} <: SectorData{N, T}
    data :: SortedDict{U1Sector{N}, Array{T, N}}
end
mutable struct SectorDiagonal{T} <: SectorData{2, T}
    data :: SortedDict{U1Sector{2}, Diagonal{T}}
end

SectorArray(data...) = SectorData(data)
SectorDiagonal(data...) = SectorData(data)

@inline sectors(A::SectorData) = [s for (s,d) in A.data]

"""
    _allsectorsandsizes(charge, Vs)

Generate all sectors for a tuple of VSpace that add up to a given
`charge` and are thus sectors are the corresponding SymTensor. This is
recursive function calling itself with Vs[1:N-1], etc. The output
sector is sorted by construction.
"""
function _allsectorsandsizes(charge::C, Vs::NTuple{N, U1Space}) where{C, N}
    if length(Vs) < 2
        d = dim(Vs[1], charge)
        if d > 0
            return [Sector(charge)], [(d,)]
        end
        return Sector{C, 1}[], NTuple{1, Int}[]
    end

    sects = Sector{C, N}[]
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
