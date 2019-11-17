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

abstract type AbstractQModelInteraction{D, N, T} end
struct QModelInteraction{D, N, T} <: AbstractQModelInteraction{D, N, T}
    amp     :: T
    ucidxs  :: NTuple{N, Int}
    offsets :: NTuple{N, NTuple{D, Int}}
    terms   :: Vector{NTuple{N, Matrix{T}}}
end

struct SymQModelInteraction{D, N, T} <: AbstractQModelInteraction{D, N, T}
    amp     :: T
    ucidxs  :: NTuple{N, Int}
    offsets :: NTuple{N, NTuple{D, Int}}
    terms   :: Vector{NTuple{N, SymMatrix{T}}}
end

struct FermionQModelInteraction{D, N, T} <: AbstractQModelInteraction{D, N, T}
    amp     :: T
    ucidxs  :: NTuple{N, Int}
    offsets :: NTuple{N, NTuple{D, Int}}
    terms   :: Vector{NTuple{N, FermionOp}}
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
struct UnitCellQModel{Q<:AbstractQType, D} <: AbstractQModel{Q, D}
    qtype   :: Q
    lattice :: QLattice{D}
    inters  :: Vector{AbstractQModelInteraction}
end

dimension(::AbstractQModel{Q, D}) where {Q, D} = D

function tikzlattice(model::UnitCellQModel,
                     filename::String)
    isfile(filename) && error("file already exists: $filename")
    D = dimension(model)
    D == 2 || error("only 2D plotting is support for now!")
    lattice = model.lattice

    open(filename, "w") do f
        write(f, """\\documentclass[tikz, border=1mm]{standalone}
  \\usepackage{tikz,pgf}
  %\\usetikzlibrary{calc,arrows,arrows.meta}
  %\\usetikzlibrary{decorations.markings}

  \\begin{document}

  \\begin{tikzpicture}
  """)
        write(f, """
  \\tikzset{
    crbond/.style={thick, blue!70!black},
    bond/.style={thin, double, brown!80!yellow},
    ball/.style={inner color=blue, ball color=green!20!black}
  }

  \\def\\l{1cm}
  \\def\\h{\\l * 0.866}

  \\def\\ya{0 * \\h}
  \\def\\yb{1 * \\h}

  \\def\\radius{\\l/10}\n
""")
        for is in Iterators.product([1:l for l in model.lattice.sizes]...)
            for interaction in model.inters
                support(interaction) == 2 || continue
                n1, n2 = interaction.ucidxs
                off1, off2 = interaction.offsets
                if lattice.bc == :OBC
                    x_uc1 = [off1[i] .+ is[i] for i=1:D]
                    x_uc2 = [off2[i] .+ is[i] for i=1:D]
                    if all([1 <= x_uc1[i] <= lattice.sizes[i] for i=1:D]) &
                        all([1 <= x_uc2[i] <= lattice.sizes[i] for i=1:D])
                        one = reduce(coordsum,
                                     #[(x_uc1 .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                     [reverse(x_uc1 .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                     init=(0,0)).+ lattice.unitc.sites[n1]
                        two = reduce(coordsum,
                                     #[(x_uc2 .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                     [reverse(x_uc2 .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                                     init=(0,0)).+ lattice.unitc.sites[n2]
                        write(f, "\\draw [bond] $one -- $two;\n")
                    end
                else
                    error("boundary not supported yet!")
                end
            end
        end

        for ls in Iterators.product([1:l for l in model.lattice.sizes]...)
            x_uc = reduce(coordsum,
                          #[(ls .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                          [reverse(ls .- 1)[i] .* lattice.unitc.as[i] for i=1:D],
                          init=(0,0))
            for n in 1:lattice.unitc.n
                x_site = x_uc .+ lattice.unitc.sites[n]
                write(f, "\\shade [ball] $x_site circle (\\radius);\n")
            end
        end
        write(f, "\\end{tikzpicture}\n")
        write(f, "\\end{document}")
    end
end
