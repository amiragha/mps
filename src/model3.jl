struct UnitCell{D}
    n     :: Int
    sites :: Vector{NTuple{D, Float64}}
    as    :: Vector{NTuple{D, Float64}}

    function UnitCell{D}(n     :: Int,
                         sites :: Vector,
                         as    :: Vector) where{D}
        length(sites) == n || error("number of sites don't match!")
        length(as) == D || error("dimension of unitcell vector don't match")
        new{D}(n, sites, as)
    end
end

UnitCell{1}(n::Int, sites::Vector{Float64}, a::Float64) =
    UnitCell{1}(n, [(site,) for site in sites], [(a,)])

struct QLattice{D}
    unitc :: UnitCell{D}
    sizes :: NTuple{D, Int}
    bc    :: Symbol
end

abstract type AbstractQType end
struct SpinType <: AbstractQType
    d :: Int
end

abstract type AbstractQInteraction end
struct QModelInteraction{D, N, T} <: AbstractQInteraction
    amp     :: T
    sites   :: NTuple{N, Int}
    offsets :: NTuple{N, NTuple{D, Int}}
    repeat  :: Union{NTuple{N, Int}, Nothing}
    terms   :: Vector{NTuple{N, Matrix{T}}}
end

# This is the generic linear QTerm
struct QInteraction{T, N} <: AbstractQInteraction
    amp   :: T
    sites :: NTuple{N, Int}
    terms :: Vector{NTuple{N, Matrix{T}}}
end
support(::QInteraction{T,N}) where{T, N} = N

abstract type AbstractQModel end
struct UnitCellQModel{Q<:AbstractQType, D} <: AbstractQModel{Q, D}
    qtype   :: Q
    lattice :: QLattice{D}
    inters  :: Vector{QModelInteraction}
end

dimension(::AbstractQModel{Q, D}) where {Q, D} == D
