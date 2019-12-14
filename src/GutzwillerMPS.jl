module GutzwillerMPS

using MatrixProductStateTools
using SymTensors
using LinearAlgebra
using TensorOperations

RLorCX = Union{Float64, ComplexF64}

include("gutzwiller.jl")

export gutzwillerexact, zipandgutzwiller!

end
