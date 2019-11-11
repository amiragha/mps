module QuantumModel

using LinearAlgebra
using SparseArrays

using SymTensors

include("spin_definitions.jl")
include("model3.jl")
include("model_examples.jl")

export UnitCell
export QLattice
export AbstractQType, SpinType
export QModelInteraction, QInteraction
export UnitCellQModel

export dimension

export sz_half, sp_half, sm_half
export sz_half_U1sym
export sp_half_U1sym
export sm_half_U1sym
export spinoperators

export j1j2model

end
