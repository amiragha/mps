using Test

using Random
using LinearAlgebra
using TensorOperations
using KrylovKit

using QuantumModels
using MatrixProductStateTools
using ExactDiagonalizationTools

using SymTensors

function Base.isapprox(A::Tuple, B::Tuple)
    length(A) == length(B) || error("Not equal length")
    all(isapprox(A[i],B[i]) for i=1:length(A))
end

include("modeltests.jl")
include("mpstests.jl")
include("mpotests.jl")
include("symtensortests.jl")
