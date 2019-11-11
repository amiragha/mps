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

QLattice(uc::UnitCell{D}, lx::Int, bc::Symbol) where{D}=
QLattice{1}(uc, (lx,), bc)

# find the linear index of a site with a particular unitcell at some
# particular position
function sitelinearindex(lattice::QLattice{D},
                         ucidx::Int,
                         x_uc::NTuple{D, Int}) where{D}
    if lattice.bc == :OBC
        if all([1 <= x_uc[i] <= lattice.sizes[i] for i=1:D])
            return ucidx + lattice.unitc.n *
                sum((x_uc .- 1) .* [1, cumprod([lattice.sizes...])[1:end-1]...])
        end
    else
        error("boundary not supported yet!")
    end
    return 0
end

abstract type AbstractQType end
struct SpinType <: AbstractQType
    d :: Int
end

abstract type AbstractQInteraction end
struct QModelInteraction{D, N, T} <: AbstractQInteraction
    amp     :: T
    ucidxs  :: NTuple{N, Int}
    offsets :: NTuple{N, NTuple{D, Int}}
    #repeat  :: Union{NTuple{D, Int}, Nothing}
    terms   :: Vector{NTuple{N, Matrix{T}}}
end
support(::QModelInteraction{D, N, T}) where{D, N, T} = N
eltype(::QModelInteraction{D, N, T}) where{D, N, T} = T

# This is the generic linear QTerm
struct QInteraction{T, N} <: AbstractQInteraction
    amp   :: T
    sites :: NTuple{N, Int}
    terms :: Vector{NTuple{N, Matrix{T}}}
end
support(::QInteraction{T,N}) where{T, N} = N

removehead(A::QInteraction{T, N}) where {T, N} =
    QInteraction{T, N-1}(A.amp, A.sites[2:N], [term[2:N] for term in A.terms])

abstract type AbstractQModel{Q, D} end
struct UnitCellQModel{Q<:AbstractQType, D} <: AbstractQModel{Q, D}
    qtype   :: Q
    lattice :: QLattice{D}
    inters  :: Vector{QModelInteraction}
end

dimension(::AbstractQModel{Q, D}) where {Q, D} = D
