module GutzwillerMPS

using MatrixProductStateTools
using LinearAlgebra
using TensorOperations

RLorCX = Union{Float64, ComplexF64}

include("gutzwiller.jl")

export zipandgutzwiller_exact, zipandgutzwiller!

end
