module SymTensors

using LinearAlgebra

import Base: convert, size, show, isequal, ==, *, conj
import Base.intersect

include("leg.jl")
include("symtensor.jl")
include("symreleg.jl")
include("contract.jl")
include("symtools.jl")

export STLeg, SymTensor, SymVector
export eye
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
export fuselegs
export defuse_leg
export permutelegs
export contract

export change_sign, change_legsign
export mapcharges
export invlegs

export change_nzblk!
export inv_perm

export svdsym

end
