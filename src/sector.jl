const Sector{C, N} = NTuple{N, C}
const U1Sector{N} = Sector{Int, N}

@inline sum(S::Sector) = sum(S)
@inline dual(S::Sector) = Sector(dual(s) for s in S)

"""
we always sort sectors based on charges The sorting is column
major. That means the leftmost (first) VSpace changes charge faster
(first)!
"""
@inline isless(s1::Sector{N, C}, s2::Sector{N, C}) where {N, C} =
 reverse(s1) < reverse(s2)

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
        d = getdim(Vs[1], charge)
        if d > 0
            return [(charge,)], [(d,)]
        end
    end

    sects = NTuple{N, C}[]
    sizes = NTuple{N, Int}[]
    for (c,d) in Vs[N]
        sect_head, size_head = _allsectorsandsizes(charge - c, Vs[1:N-1])
        for i in eachindex(sect_head)
            push!(sects, (sect_head[i]..., c))
            push!(sizes, (size_head[i]..., d))
        end
    end
    sects, sizes
end
