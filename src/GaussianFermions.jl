module GaussianFermions

using MatrixProductStateTools
using SymTensors
using LinearAlgebra

include("gmps.jl")
include("gmps2mps.jl")

export GaussianMPS
export corrmat2gmps
export gmps2mps

end
