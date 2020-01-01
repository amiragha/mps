module QuantumModels

using LinearAlgebra
using SparseArrays

#using DataStructures
using SymTensors

import Base: eltype, *

include("spin_definitions.jl")
include("model3.jl")
include("spinmodels.jl")
include("fermionmodels.jl")
include("modeldraw.jl")

export UnitCell
export QLattice
export AbstractQType, SpinType
export QModelInteraction, QInteraction
export SymQModelInteraction, QInteraction
export FermionQModelInteraction
export UnitCellQModel
export QTerm

export dimension
export changeboundary
export changesize
export sitelinearindex
export supportrange
export support
export removehead
export largestxrange

export Fermion
export fermionoperators

export spinhalf
export spinoperators
export ringexchangeoperator
export nbodyopexpansion
export permutespins

export chainunitcell
export triangularunitcell

# spin model examples
export j1j2model
export triangularspinmodel

# fermion model examples
export t1t1model
export triangularhopping

export tikzlattice
end
