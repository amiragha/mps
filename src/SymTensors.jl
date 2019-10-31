module SymTensors

using LinearAlgebra
using Random

import Base: convert, size, show, isequal, ==, isapprox, *, conj
import Base: eltype, similar, copyto!
import Base: fill, fill!, rand
import LinearAlgebra: mul!, rmul!, axpy!, axpby!, dot, norm
import Base.intersect

include("leg.jl")
include("symtensor.jl")
include("symreleg.jl")
include("contract.jl")
include("symtools.jl")

export +, *
export isequal, ==
export isapprox
export STLeg, SymTensor, SymMatrix, SymVector
export eye
#export randSymTensor, fillSymTensor
export FusedCharge

export symMatrix
export _allsectorsandsizes
export _allsectorsandsizes2
export _sectorlessthan
export _sectorisallowed
#export _possible_fuse_patterns

#export array_representation
export array
export alldims, accdims, getdim, fulldims, signs
export numoflegs
export _sectors_sortperm

export fuse_set
export fuse_conseqlegs
export fuselegs
export defuse_leg
export permutelegs
export removedummyleg
export contract

export change_sign, change_legsign
export mapcharges
export invlegs

export set_sector!

export svdsym

export eltype, similar, copyto!
export mul!, rmul!, axpy!, axpby!, dot, norm
export fill, fill!

export fermionswapgate
end
