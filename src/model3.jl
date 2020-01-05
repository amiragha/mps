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

@inline coordsum(A::NTuple{N}, B::NTuple{N}) where {N} = A .+ B

struct QLattice{D}
    unitc :: UnitCell{D}
    sizes :: NTuple{D, Int}
    bcs   :: NTuple{D, Symbol}
end

QLattice(uc::UnitCell{1}, lx::Int, bc::Symbol) where{D} =
    QLattice{1}(uc, (lx,), (bc,))

@inline dimension(lattice::QLattice{D}) where{D} = D

"return a changed boundary version of `lattice` in the `d`th direction"
function changeboundary(lattice::QLattice, d::Int, boundary::Symbol)
    D = dimension(lattice)
    1 <= d <= D || error("can't change boundary!")
    QLattice(lattice.unitc,
             lattice.sizes,
             Tuple([lattice.bcs[1:d-1]..., boundary, lattice.bcs[d+1:D]...]))
end

"return a chaged size version of `lattice` in the `d`th direction"
function changesize(lattice::QLattice, d::Int, dim::Int)
    D = dimension(lattice)
    1 <= d <= D || error("can't change size!")
    QLattice(lattice.unitc,
             Tuple([lattice.sizes[1:d-1]..., dim, lattice.sizes[d+1:D]...]),
             lattice.bcs)
end

"""
    sitelinearindex(lattice, ucidx, x_uc)

for a given lattice `lattice`, unite cell site `ucidx`, and a
corrdinate of unitcell `x_uc` which can be any integer number, find
the site linear index of that site if allowed according to the
boundary condition!

"""
function sitelinearindex(lattice::QLattice{D},
                         ucidx::Int,
                         x_uc::NTuple{D, Int}) where {D}

    x_uc_vec = zeros(Int, D)
    crossings = zeros(Int, D)
    for i in 1:D
        if lattice.bcs[i] == :OBC
            if 1 <= x_uc[i] <= lattice.sizes[i]
                x_uc_vec[i] = x_uc[i]
            else
                return nothing, crossings
            end

        elseif lattice.bcs[i] in [:PBC, :APBC]
            x_uc_vec[i] = mod(x_uc[i] - 1, lattice.sizes[i]) + 1
            crossings[i] = fld(x_uc[i] - 1, lattice.sizes[i])

        elseif lattice.bcs[i] == :INF
            if i < D
                error("only the last dimensions can have INF boundary conition!")
            end
            x_uc_vec[i] = x_uc[i]

        else
            error("Exhaustive check!")
        end
    end
    index = ucidx + lattice.unitc.n *
        sum((x_uc_vec .- 1) .* [1, cumprod([lattice.sizes...])[1:D-1]...])

    return index, Tuple(crossings)
end

abstract type AbstractQType end
struct SpinType <: AbstractQType
    d :: Int
end

# because different species of fermions should be different for
# anticommutation relations
struct Fermion{Name} <: AbstractQType
    Fermion{Name}() where {Name} = new{typeassert(Name, Symbol)}()
end

abstract type QuantumOperator end
abstract type FermionOp end

struct FCreate{F} <: FermionOp
    FCreate{F}() where {F} = new{typeassert(F, Fermion)}()
end

struct FAnnihilate{F} <: FermionOp
    FAnnihilate{F}() where {F}= new{typeassert(F, Fermion)}()
end

function fermionoperators(f::Fermion)
    FCreate{f}(), FAnnihilate{f}()
end
abstract type AbstractQInteraction{T, N} end

struct QAmpTerm{O, N, T}
    amp :: T
    ops :: NTuple{N, O}
end

abstract type AbstractQModelInteraction{D, N, T} end
struct QModelInteraction{D, N, T} <: AbstractQModelInteraction{D, N, T}
    amp     :: T
    ucidxs  :: NTuple{N, Int}
    offsets :: NTuple{N, NTuple{D, Int}}
    terms   :: Vector{QAmpTerm{Symbol,N,T}}
end

# struct SymQModelInteraction{D, N, T} <: AbstractQModelInteraction{D, N, T}
#     amp     :: T
#     ucidxs  :: NTuple{N, Int}
#     offsets :: NTuple{N, NTuple{D, Int}}
#     terms   :: Vector{NTuple{N, SymMatrix{T}}}
# end

struct FermionQModelInteraction{D, N, T} <: AbstractQModelInteraction{D, N, T}
    amp     :: T
    ucidxs  :: NTuple{N, Int}
    offsets :: NTuple{N, NTuple{D, Int}}
    terms   :: Vector{NTuple{N, FermionOp}}
end

function supportrange(inter::AbstractQModelInteraction{D, N, T}) where{D,N,T}
    [(minimum([offset[d] for offset in inter.offsets]),
      maximum([offset[d] for offset in inter.offsets])) for d in 1:D]
end

support(::AbstractQModelInteraction{D, N, T}) where{D, N, T} = N
eltype(::AbstractQModelInteraction{D, N, T}) where{D, N, T} = T

# This is the generic linear QTerm
struct QInteraction{T, N} <: AbstractQInteraction{T, N}
    amp   :: T
    sites :: NTuple{N, Int}
    terms :: Vector{NTuple{N, Matrix{T}}}
end

support(::AbstractQInteraction{T,N}) where{T, N} = N
removehead(A::QInteraction{T, N}) where {T, N} =
    QInteraction{T, N-1}(A.amp, A.sites[2:N], [term[2:N] for term in A.terms])

struct QTerm{OP, N}
    sites :: NTuple{N, Int}
    ops   :: NTuple{N, OP}
end

eltype(term::QTerm{OP}) where {OP} = eltype(OP)
support(::QTerm{OP, N}) where {OP, N} = N
removehead(A::QTerm{OP, N}) where {OP, N} =
    QTerm{OP, N-1}(A.sites[2:N], A.ops[2:N])

# function *(a::T, term::QTerm) where {T<:Number}
#     T == eltype(eltype(term)) || error("oops, $T vs $(eltype(eltype(term)))")
#     N = support(term)
#     QTerm{OP, N}(term.sites, [op for op in term.ops[1:N-1]]..., a*term.ops[N])
# end

abstract type AbstractQModel{Q, D} end
struct UnitCellQModel{Q<:AbstractQType, D, J<:AbstractQModelInteraction} <: AbstractQModel{Q, D}
    qtype    :: Q
    lattice  :: QLattice{D}
    symmetry :: Type{<:AbstractCharge}
    inters   :: Vector{J}
end

@inline dimension(::AbstractQModel{Q, D}) where {Q, D} = D
@inline numofsites(m::UnitCellQModel) = prod(m.lattice.sizes)*m.lattice.unitc.n

function changeboundary(model::UnitCellQModel, d::Int, boundary::Symbol)
    UnitCellQModel(model.qtype,
                   changeboundary(model.lattice, d, boundary),
                   model.inters)
end

function changesize(model::UnitCellQModel, d::Int, dim::Int)
    UnitCellQModel(model.qtype,
                   changesize(model.lattice, d, dim),
                   model.inters)
end

function largestxrange(model::UnitCellQModel)
    D = dimension(model)
    ranges = [supportrange(inter)[D] for inter in model.inters]
    maximum([r[2] for r in ranges]) - minimum([r[1] for r in ranges])
end
