module SymTensors

using LinearAlgebra

import Base: size, show, isequal, ==, *, conj

#include("bsmatrix.jl")
include("symtensor.jl")
include("symreleg.jl")
include("contract.jl")
include("symtools.jl")

export STLeg, SymTensor, SymVector
export randSymTensor, fillSymTensor
export FusedCharge

export symMatrix
export _possible_fuse_patterns

export array_representation
export alldims, accdims, getdim, fulldims, signs
export numoflegs
export _sectors_sortperm

export fuse_set
export fuse_conseqlegs
export defuse_leg
export permutelegs
export contract

export change_sign, change_legsign
export invlegs

export inv_perm

export svdsym
end
