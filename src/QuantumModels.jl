module QuantumModels

using LinearAlgebra
using SparseArrays

using SymTensors

import Base: eltype, *

include("spin_definitions.jl")
include("model3.jl")
include("spinmodels.jl")
include("fermionmodels.jl")

export UnitCell
export QLattice
export AbstractQType, SpinType
export QModelInteraction, QInteraction
export SymQModelInteraction, QInteraction
export FermionQModelInteraction
export UnitCellQModel
export QTerm

export dimension
export sitelinearindex
export support
export removehead

export spinhalf
export spinoperators
export ringexchangeoperator
export nbodyopexpansion
export permutespins

# spin model examples
export chainunitcell
export triangularunitcell
export j1j2model
export triangularspinmodel

end
