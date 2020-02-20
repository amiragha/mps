module SymTensors

using LinearAlgebra
using Random

import Base: convert, size, show, isless, isequal, ==, isapprox, *, conj, zero, inv
import Base: eltype, similar, copyto!
import Base: fill, fill!, rand, sum
import LinearAlgebra: mul!, rmul!, axpy!, axpby!, dot, norm, normalize!, svd, diag
import Base.intersect

#using DataStructures
include("sorteddict.jl")
export SortedDict, onlysemitokens

include("group.jl")
include("vspace.jl")
include("sector.jl")
#include("productspace.jl")
include("symtensor.jl")
include("trivialtensor.jl")
include("symreleg.jl")
include("contract.jl")
include("symtools.jl")

export +, *, zero, inv
export isequal, ==
export isapprox
export issimilar
export AbstractSymTensor, AbstractSymMatrix
export VectorSpace#, ProductSpace
export AbstractCharge
export Trivial
export U1Space, U1, U1Charge
export Z2Charge
export Sector
export SymTensor, SymMatrix, SymDiagonal, SymVector
export eye

export vtype

export SectorArray, SectorDiagonal
export _allsectorsandsizes
export _sectorisallowed

export array
export dim, dims, charges
export space
#export alldims, accdims, getdim, fulldims, signs
export rank, sectors, blocks
export fuse


#export _sectors_sortperm
export fuselegs
export splitleg
export permutelegs
export dropdummyleg
export trimleg!
export contract

export memoryrepr
export dual, isdual
export mapcharges
#export invlegs

export index_sector
export get_sector!
export set_sector!

export _svd_

export eltype, similar, copyto!
export mul!, rmul!, axpy!, axpby!, dot, norm
export normalize!
export fill, fill!
export fill_linearindex

export fermionswapgate

export isleftisometry, isrightisometry, isunitary

end
