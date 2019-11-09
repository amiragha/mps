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

chainunitcell = UnitCell{1}(1, [0.], +1.)
squareunitcell = UnitCell{2}(1, [(0, 0)], [(0, 1), (1, 0)])
triangularunitcell = UnitCell{2}(1, [(0, 0)], [(1, 0), (0.5, sin(pi/3))])

struct QLattice{D}
    unitc :: UnitCell{D}
    sizes :: NTuple{D, Int}
    bc    :: Symbol
end

struct QSpecies
    d :: Int
end

spinhalf = QSpecies(2)

abstract type AbstractQTerm end
struct QModelTerm{D, N, T} <: AbstractQTerm
    amp     :: T
    site    :: NTuple{N, Int}
    offsets :: NTuple{N, NTuple{D, Int}}
    repeat  :: Union{NTuple{N, Int}, Nothing}
    terms   :: Vector{NTuple{N, Matrix{T}}}
end

# This is the generic linear QTerm
struct QTerm{T, N} <: AbstractQTerm
    amp   :: T
    sites :: NTuple{N, Int}
    terms :: Vector{NTuple{N, Matrix{T}}}
end

struct QModel{D}
    species :: QSpecies
    lattice :: QLattice{D}
    terms   :: Vector{QTerm}
end

heisj1 = QModelTerm{1, 2, Float64}(
    j1,
    (1, 1),
    ((0,), (1,)),
    (1,),
    [(sz, sz), (sp, 0.5*sm), (sm, 0.5*sp)])
heisj2 = QModelTerm{1, 2, Float64}(
    j2,
    (1, 1),
    ((0,), (2,)),
    (1,),
    [(sz, sz), (sp, 0.5*sm), (sm, 0.5*sp)])

j1j2 = QModel(spinhalf,
              QLattice(chainunitcell, lx),
              [heisj1, heisj2])

function generatempo()
    #D = dimension(model)

    # put all sites in order (we should be able to choose this order somehow!)
    # now generate all terms sorted with the starting site
    # below fterms(finishing terms), sterms(starting terms), lterms(localterms), cterms(crossingterms)
    # now generate the MPO
    d = 2
    tensors = Vector{Array{Float64, 4}}(undef, n_sites)
    ldim = 1
    for n in 1:n_sites
        rdim = 2 + n_sterms + n_cterms
        W = zeros(T, ldim, d, rdim, d)
        W[1,:,1,:] = I(d)
        W[ldim, :, rdim, :] = I(d)
        W[dlim, :, 1, :] = sum(lterms)
        for t in 1:n_fterms
            W[t+1, :, 1, :] = nextop(fterms[t])
        end
        for t in 1:n_startingterms
            w[ldim, :, t+1, :] = nextop(sterms[t])
        end
        for t in 1:n_cterms
            W[1+n_fterms+t, :, 1, :] = nextop(cterms[t])

            tensors.push!(W)
        end

        # return MPO
    end
