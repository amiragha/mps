module ExactDiagonalizationTools

using QuantumModels

using LinearAlgebra
using SparseArrays
using DataStructures
using SymTensors
using Random
using QuadGK
RLorCX = Union{Float64, ComplexF64}
⊗ = kron

import LinearAlgebra: normalize!

import Base: +, *
include("model.jl")

include("ket.jl")
include("xxz.jl")
include("longxxz.jl")
include("qitf.jl")
include("hopping.jl")
include("correlationmatrix.jl")

export KetState
export opextend
export normalize!
export randket
export apply1site
export measure1point

export xxz_hamiltonian
export xxz_longrange
export j1j2_explicit
export qitf_hamiltonian, qitf_bondtensor
export qitf_energy_exact

export hopping_chain
export nnhoppingchain
export correlationmatrix

# Model stuff to be separated
export Point2D, Bond2D, UnitCell2D
export triangular_unitcell
export enlargeunitcell
export makehamiltonian
export makemodel
export makemodelJW

export makemodelJW_u1sym

export generatebdg

end
