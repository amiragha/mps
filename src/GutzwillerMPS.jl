module GutzwillerMPS

using MatrixProductStateTools
using SymTensors
using LinearAlgebra
using TensorOperations

include("gutzwiller.jl")

export gutzwillerexact, zipandgutzwiller!

end
