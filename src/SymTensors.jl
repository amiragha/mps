module SymTensors

using LinearAlgebra
using Random

import Base: convert, size, show, isequal, ==, isapprox, *, conj
import Base: eltype, similar, copyto!
import Base: fill, fill!, rand
import LinearAlgebra: mul!, rmul!, axpy!, axpby!, dot, norm, normalize!
import Base.intersect

include("leg.jl")
include("symtensor.jl")
include("symreleg.jl")
include("contract.jl")
include("symtools.jl")

export +, *
export isequal, ==
export isapprox
export issimilar
export AbstractSymTensor, AbstractSymMatrix
export STLeg, SymTensor, SymMatrix, SymDiagonal, SymVector
export eye
export FusedCharge

export symMatrix
export _allsectorsandsizes
export _allsectorsandsizes2
export _sectorlessthan
export _sectorisallowed

export array
export alldims, accdims, getdim, fulldims, signs
export numoflegs
export _sectors_sortperm

export fuse
export fuselegs
export unfuseleg
export permutelegs
export removedummyleg
export trimleg!
export contract

export reverseleg
export mapcharges
export invlegs

export index_sector
export get_sector!
export set_sector!

export svdsym

export eltype, similar, copyto!
export mul!, rmul!, axpy!, axpby!, dot, norm
export normalize!
export fill, fill!

export fermionswapgate

export isleftisometry, isrightisometry, isunitary

end
