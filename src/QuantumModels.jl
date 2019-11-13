module QuantumModels

using LinearAlgebra
using SparseArrays

using SymTensors

import Base: eltype, *

include("spin_definitions.jl")
include("model3.jl")
include("model_examples.jl")

export UnitCell
export QLattice
export AbstractQType, SpinType
export QModelInteraction, QInteraction
export SymQModelInteraction, QInteraction
export UnitCellQModel
export QTerm

export dimension
export sitelinearindex
export support
export removehead

export spinhalf
export sz_half, sp_half, sm_half
export sz_half_U1sym
export sp_half_U1sym
export sm_half_U1sym
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
