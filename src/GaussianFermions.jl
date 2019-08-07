module GaussianFermions

using MatrixProductStateTools
using LinearAlgebra

RLorCX = Union{Float64, ComplexF64}

include("fishman.jl")
include("fishman2mps.jl")

export FishmanGateSet
export generate_fishmangates
export fishman2mps

end
