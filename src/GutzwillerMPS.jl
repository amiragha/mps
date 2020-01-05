module GutzwillerMPS

using MatrixProductStateTools
using SymTensors
using LinearAlgebra
using TensorOperations

include("gutzwiller.jl")

export gutzwillerexact
export zipandgutzwiller!
export _tensorproductzip!
export _applygutzwillerzip!

end
