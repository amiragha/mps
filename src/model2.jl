# seems like an overkill! isn't it!
abstract type QuantumTerm end
abstract type FermionTerm <: QuantumTerm end

struct ChemicalPotential{T<:Number} <: FermionTerm
    op ::
    amp :: T
end

abstract type UnitCell{D} end

struct UnitCellFermion{D} <: UnitCell{D}
    sites :: Vector{Vector{Float64}}
    vects :: NTuple{D, Vector{Float64}}
    hops ::
    function UnitCellFermion(n, sites vects, bonds)
        length(vects) == N || error("uc n_vectors and dim don't match!")
        all(length.(vects) == N) || error("Dim of vectors and uc don't match!")
        new(sites, vects)
    end
end

struct UnitCellSpin{D} <: UnitCell{D}
    sites :: Vector{Vector{Float64}}
    vects :: NTuple{D, Vector{Float64}}
    twospin ::
    function UnitCellSpin(sites vects, bonds)
        length(vects) == N || error("Number of UnitCell vectors and dimensinos don't match!")
        all(length.(vects) == N) || error("Dimension of vectors and UnitCell don't match!")
        new(sites, vects)
    end
end

@inline n_sites(uc::UnitCell) = length(uc.sites)
@inline dimension(uc::UnitCell{D}) where {D} = D

abstract type QuantumModel end

struct FermionModel{D} <: QuantumModel
    uc :: UnitCell{D}
    Ls :: NTuple{D, Int}

    hops ::
end

@inline size
