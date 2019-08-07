module ExactDiagonalizationTools

using LinearAlgebra
using SparseArrays
RLorCX = Union{Float64, ComplexF64}

include("xxz.jl")
include("hopping.jl")
include("correlationmatrix.jl")

export xxz_hamiltonian
export hopping_chain
export correlationmatrix

end
