module ExactDiagonalizationTools

using LinearAlgebra
using SparseArrays
using SymTensors
RLorCX = Union{Float64, ComplexF64}

include("spin_definitions.jl")
include("xxz.jl")
include("hopping.jl")
include("correlationmatrix.jl")

export sz_half, sp_half, sm_half
export sz_half_U1sym

export xxz_hamiltonian
export hopping_chain
export correlationmatrix

end
