abstract type AbstractSymTensor{S, T, N} end
const AbstractSymMatrix{S, T} = AbstractSymTensor{S, T, 2}

"""
A symmetric tensor is represented as a tensor map over the space of
V_row (codomain) with tensor product structure of (((V1 ⊗ V2) ⊗ V3) ⊗
...) and the space of V_col (domain) with tensor product structure of
(...⊗ (V_{n-3} ⊗ (V_{n-1} ⊗ V_n))). The blocks is also represented as a
matrix
"""

mutable struct SymTensor{S, T, N} <: AbstractSymTensor{S, T, N}
    charge   :: S
    space    :: NTuple{N, VectorSpace{S}}
    blocks   :: SortedDict{Sector{S, N}, Array{T, N}}

    function SymTensor(charge :: S,
                       space  :: NTuple{N, VectorSpace{S}},
                       blocks   :: SortedDict{Sector{S, N}, Array{T, N}}) where{S,T,N}
        sects, sizes = _allsectorsandsizes(charge, space)
        i = 1
        for (s,d) in blocks
            s == sects[i] || throw("$(keys(blocks)) SectorMismatch(), $s vs $(sects[i])")
            size(d) == sizes[i] || throw(SizeMismatch())
            i += 1
        end
        new{S,T,N}(charge, space, blocks)
    end
    SymTensor{S,T,N}(c,s,d) where{S,T,N} = SymTensor(c,s,d)
end

# function SymTensor{S,T1,N}(A::SymTensor{S,T2,N}) where {S,T1, T2, N}
#     SymTensor{S,T1,N}(A.charge, A.space, SortedDict{Sector{S,N}, T1}(A.blocks))
# end

const SymMatrix{S, T} = SymTensor{S, T, 2}
SymMatrix(c::S, s::NTuple{2,T}, d) where{S, T} = SymMatrix{S,T}(c,s,d)

const SymVector{S, T} = SymTensor{S, T, 1}

const U1Tensor{T, N} = SymTensor{Int, T, N}
const U1Matrix{T} = SymMatrix{Int, T}
const U1Vector{T} = SymVector{Int, T}

mutable struct SymDiagonal{S, T<:Number} <: AbstractSymMatrix{S, T}
    charge :: S
    space  :: NTuple{2, VectorSpace{S}}
    blocks   :: SortedDict{Sector{S, 2}, Diagonal{T,Vector{T}}}

    function SymDiagonal(charge :: S,
                         space  :: NTuple{2, VectorSpace{S}},
                         blocks   :: SortedDict{Sector{S,2}, Diagonal{T,Vector{T}}}) where{S,T}
        sects, sizes = _allsectorsandsizes(charge, space)
        i = 1
        for (s,d) in blocks
            s == sects[i] || throw("$blocks SectorMismatch(), $s vs $(sects[i])")
            size(d) == sizes[i] || throw(SizeMismatch())
            i += 1
        end
        # dims(first(space)) == dims(last(space)) ||
        #     throw("SymDiagonal has to be square $space")
        new{S,T}(charge, space, blocks)
    end
end

SymDiagonal{S, T}(c,s,b) where{S,T} = SymDiagonal(c,s,b)
const U1Diagonal{T} = SymDiagonal{Int, T}

function SymVector(charge::S, v::Vector{T}, isdual::Bool=false) where {S<:AbstractCharge, T}
    V = VectorSpace{S}([charge => length(v)], isdual)
    SymVector{S, T}(charge, (V,), SortedDict(Sector(charge) => v))
end

function convert(::Type{SymTensor{S, T1, N}},
                 A::SymTensor{S, T2, N}) where {S, T1, T2, N}
    SymTensor{S, T1, N}(A.charge, A.space, A.blocks)
end

@inline space(A::AbstractSymTensor) = A.space
@inline space(A::AbstractSymTensor, l::Int) = A.space[l]
@inline LinearAlgebra.rank(::AbstractSymTensor{S, T, N}) where {S, T, N} = N
@inline charge(A::AbstractSymTensor) = A.charge
@inline sectors(A::AbstractSymTensor) = collect(keys(A.blocks))
@inline blocks(A::AbstractSymTensor) = A.blocks

@inline size(A::AbstractSymTensor) = Tuple(dim(V) for V in A.space)
@inline size(A::AbstractSymTensor, l::Int) = dim(A.space[l])

@inline issimilar(A::T, B::T) where {T<:AbstractSymTensor} =
    A.charge == B.charge && A.space == B.space

@inline isequal(A::T, B::T) where {T<:AbstractSymTensor} =
    issimilar(A, B) && isequal(A.blocks, B.blocks)

@inline isapprox(A::AbstractSymTensor, B::AbstractSymTensor) =
    issimilar(A, B) && all(isapprox(A.blocks[c], B.blocks[c]) for c in keys(A.blocks))

@inline ==(A::T, B::T) where {T<:AbstractSymTensor} = isequal(A, B)

@inline eltype(::AbstractSymTensor{S, T}) where {S, T} = T
@inline eltype(::Type{<:AbstractSymTensor{S, T}}) where {S, T} = T

@inline vtype(::AbstractSymTensor{S, T}) where {S, T} = S
@inline vtype(::Type{<:AbstractSymTensor{S, T}}) where {S, T} = S

#@inline size(A::SymVector) = size(nzblks[1])

# function index_sector(A::AbstractSymTensor{T, N},
#                       s::NTuple{N, Int}) where{T<:Number, N}
#     index = searchsortedfirst(A.sects, s, lt=_sectorlessthan)
#     (index > length(A.sects) || s != A.sects[index]) &&
#         error("sector not found!")
#     index
# end

@inline Base.getindex(A::AbstractSymTensor{S,T,N},
                      s::Sector{S, N}) where {S,T,N} = getindex(A.blocks, s)

function Base.setindex!(A::AbstractSymTensor{S,T,N},
                        blk::Array{T, N},
                        s::Sector{S, N}) where {S,T,N}
    haskey(A.blocks, s) || throw("Block doesn't exist")
    A.blocks[s] = blk
    A
end

SymTensor(f::Function, charge::S, space::NTuple{N,VectorSpace{S}}) where {S,N} =
    SymTensor(f, Float64, charge, space)

function SymTensor(f::Function, ::Type{T}, charge::S,
                   space::NTuple{N,VectorSpace{S}}) where {S,T,N}
    sects, sizes = _allsectorsandsizes(charge, space)
    blocks = SortedDict([sects[i] => f(T,sizes[i]) for i in eachindex(sects)])
    SymTensor(charge, space, blocks)
end

function rand(::Type{T}, charge::S, space::NTuple{N, VectorSpace{S}};
              seed::Int=1911) where {S,T,N}
    sects, sizes = _allsectorsandsizes(charge, space)
    rng = MersenneTwister(seed)
    blocks = SortedDict([sects[i] => rand(rng, T, sizes[i]) for i in eachindex(sects)])
    SymTensor(charge, space, blocks)
end

rand(::Type{T}, space::NTuple{N, VectorSpace{S}}) where {S,T,N} =
    rand(T, zero(S), space)
rand(charge::Int, space::NTuple{N, VectorSpace{S}}) where {S,N} =
    rand(Float64, charge, space)
rand(space::NTuple{N, VectorSpace{S}}) where {S,N} =
    rand(Float64, zero(0), space)

fill(x::T, space::NTuple{N, VectorSpace{S}}) where {S,T,N} =
    fill(x, zero(S), space)
function fill(x::T, charge::S, space::NTuple{N, VectorSpace{S}}) where {S,T,N}
    sects, sizes = _allsectorsandsizes(charge, space)
    blocks = SortedDict([sects[i] => fill(x, sizes[i]) for i in eachindex(sects)])
    SymTensor(charge, space, blocks)
end

function fill!(A::SymTensor{S,T,N}, x::T) where {S,T,N}
    for b in values(A.blocks)
        fill!(b, x)
    end
    A
end

function fill_linearindex(charge::S, space::NTuple{N, VectorSpace{S}}) where {S,N}
    sects, sizes = _allsectorsandsizes(charge, space)
    blocks = SortedDict{Sector{S, N}, Array{Int, N}}()
    p = 0
    for i in eachindex(sects)
        sz = sizes[i]
        blocks[sects[i]] = reshape(collect(p+1:p+prod(sz)), sz...)
        p += prod(sz)
    end
    SymTensor(charge, space, blocks)
end
fill_linearindex(space::NTuple{N, VectorSpace{S}}) where{S,N} =
    fill_linearindex(zero(S), space)

##TODO: Mix this with the usual UnitScaling
@inline eye(V::VectorSpace{S}) where{S} = eye(Float64, V)
@inline eye(::Type{T}, V::VectorSpace{S}) where {S, T} =
    SymMatrix(zero(S), (V, dual(V)),
              SortedDict([Sector(c, c)=>Matrix{T}(I,d,d) for (c,d) in V]))

"""
    dual(A)

Make the dual of tensor. This operation negates the: tensor charge,
leg signs, sectors as well! One can set whether or not it also
conjugates the blocks with `conjugate` which is true by default.

"""
function dual(A::AbstractSymTensor; conjugate::Bool=true)
    if !conjugate
        return typeof(A)(inv(A.charge),
                         dual.(A.space),
                         A.blocks)
        #SortedDict([inv(s)=>b for (s,b) in A.blocks]))
    end
    typeof(A)(inv(A.charge),
              dual.(A.space),
              A.blocks)
    #SortedDict([inv(s)=>conj(b) for (s,b) in A.blocks]))
end

function mapcharges(f::Function, A::AbstractSymTensor)
    SymTensor(f(A.charge),
              mapcharges.(f, A.space),
              typeof(A.blocks)(mapcharge(f, s),d for (s,d) in A.blocks))
end

function mapcharges(f::NTuple{N, Function},
                    A::AbstractSymTensor{S,T,N}) where{S,T<:Number,N}
    space = Tuple(mapcharges(f[i], A.space[i]) for i in 1:N)

    sects = [Tuple(f[i](s[i]) for i in 1:N) for s in A.sects]
    charges = [sum(s, isdual.(space)) for s in sectors(A.blocks)]
    !all(charges .== charges[1]) &&
        error("mapcharges functions are inconsistent!")
    SymTensor(charges[1], space,
              typeof(A.blocks)(mapcharge(f, s),d for (s,d) in A.blocks))
end

conj(A::AbstractSymTensor) =
    typeof(A)(A.charge, A.space, conj(A.blocks))

function array(A::AbstractSymTensor{S,T,N}) where {S,T,N}
    arrep = zeros(eltype(A), dim.(space(A))...)
    for (s,b) in A.blocks
        arrep[layout(A.space, s)...] = b
    end
    arrep
end

@inline *(A::AbstractSymTensor, a::T) where {T<:Number} =
    typeof(A)(A.charge, A.space,
              SortedDict([s=>a.*b for (s,b) in A.blocks]))

@inline *(a::T, A::AbstractSymTensor) where {T<:Number} = *(A, a)

"construct an additional SymTensor similar to A, possibly with a
different scalar type T."
function similar(A::AbstractSymTensor, T::Type=eltype(A))
    typeof(A)(A.charge, A.space,
              SortedDict([s=>similar(d, T) for (s,d) in A.blocks]))
end

"copy the contents of A to a preallocated SymTensor B"
function copyto!(B::T, A::T) where {T<:AbstractSymTensor}
    issimilar(B, A) && error("Can not copy to a  non-similar SymTensor!")
    for st in onlysemitokens(B.blocks)
        copyto!(B.blocks[st], A.blocks[st])
    end
    B
end

"out of place scalar multiplication; multiply SymTensor A with scalar α
and store the result in B"
function mul!(B::T, A::T, α) where {T<:AbstractSymTensor}
    issimilar(A, B) || error("oops")
    semitsA = collect(onlysemitokens(A.blocks))
    semitsB = collect(onlysemitokens(B.blocks))
    for i in 1:length(semitsA)
        if length(A.blocks[semitsA[i]]) != length(B.blocks[semitsB[i]])
            println()
            println(A.blocks)
            println(B.blocks)
            println(A.blocks[semitsA[i]])
            println(B.blocks[semitsB[i]])
            println()
        end
    end
    for st in onlysemitokens(B.blocks)
        mul!(B.blocks[st], A.blocks[st], α)
    end
    B
end

"in-place scalar multiplication of A with α; in particular with α =
false, A is initialized with all zeros"
function rmul!(A::T, α) where {T<:AbstractSymTensor}
    for st in onlysemitokens(A.blocks)
        rmul!(A.blocks[st],  α)
    end
    A
end

"Stores in B the result of α*A + B"
function axpy!(α,
               A::AbstractSymTensor{S,T1,N},
               B::AbstractSymTensor{S,T2,N}) where {S,T1<:Number, T2<:Number, N}
    issimilar(A, B) || error("axpy! The two matrices are not simliar!")
    for st in onlysemitokens(B.blocks)
        B.blocks[st] = α .* A.blocks[st] + B.blocks[st]
    end
    B
end

"Stores in B the result of α*A + β*B"
function axpby!(α,
                A::AbstractSymTensor{S,T1,N},
                β,
                B::AbstractSymTensor{S,T2,N}) where {S,T1<:Number, T2<:Number, N}
    issimilar(A, B) || error("axpy! The two matrices are not simliar!")
    for st in onlysemitokens(B.blocks)
        B.blocks[st] = α .* A.blocks[st] + β .* B.blocks[st]
    end
    B
end

"""
compute the inner product of two similar SymTensors which conjugates
the second one
"""
function dot(A::AbstractSymTensor{S,T1,N},
             B::AbstractSymTensor{S,T2,N}) where{S,T1<:Number, T2<:Number, N}
    issimilar(A, B) || error("The two SymTensors have to be similar to dot!")
    sum([dot(A.blocks[st], B.blocks[st]) for st in onlysemitokens(A.blocks)])
end

" compute the 2-norm of a  AbstractSymTensor"
function norm(A::AbstractSymTensor)
    sqrt(Float64(dot(A, A)))
end

normalize!(S::SymDiagonal) = rmul!(S, 1/norm(S))

# function dropdummyleg(A::AbstractSymTensor, l::Int)
#     isdummy(A.space[l]) || error("space is not dummy!")

#     N = rank(A)
#     SymTensor(A.charge,
#               (A.space[1:l-1]..., A.space[l+1:N]...),
#               typeof(A.blocks)((s[1:l-1]...,s[l+1:N]...),
#     reshape(d, size(d)[1:l-1]...,size(d)[l+1:N]...) for (s,d) in A.blocks)
#               end
#               end

function Base.show(io::IO, A::AbstractSymTensor)
    show(io, typeof(A))
    num_ch = length(A.space)
    println(io)
    #println("$(size(A)) $(typeof(A)) with $(length(A.sects)) blocks: ")
    for (s,b) in A.blocks
        print(io, "$s with size: ")
        show(io, "text/plain", b)
        println(io)
    end
    nothing
end

function trimspace!(A::AbstractSymTensor, l::Int)
    N = rank(A)
    1 <= l <= N || error("not a leg index!")
    leg = A.legs[l]
    s = leg.sign
    legcharges = BitSet()
    for sect in A.sects
        push!(legcharges, s*sect[l])
    end
    indexes = _findindexes_sorted(leg.chrs, collect(legcharges))
    leg = U1Space(leg.sign, leg.chrs[indexes], leg.dims[indexes])
    A.legs = (A.legs[1:l-1]..., leg, A.legs[l+1:N]...)
    A
end

# function negateleg(A::AbstractSymTensor, l::Int)
#     N = rank(A)
#     0 < l <= N || error("leg $l doesn't exist!")
#     legs = (A.legs[1:l-1]..., negate(A.legs[l]), A.legs[l+1:N]...)
#     SymTensor{eltype(A),N}(A.charge, legs, sects, A.nzblks)
# end

function _findindexes_sorted(A::Vector{T}, elements::Vector{T}) where{T}
    #println(A, elements)
    p = 1
    indexes = Int[]
    for i in eachindex(A)
        if A[i] == elements[p]
            push!(indexes, i)
            if p==length(elements) break end
            p += 1
        elseif A[i] > elements[p]
            error("_findindexes_sorted ", A, elements)
        end
    end
    indexes
end

# """
#     SymTensor(signs, sector, nzblock)

# Make a SymTensor with the given possible charges and only one nonzero
# sector. Find the dimensions from that.

# """
# function SymTensor(signs::NTuple{N, Int},
#                    sector::NTuple{N, Int},
#                    nzblock::Array{T, N}) where {T<:Number, N}

#     chrs = [0, 1]
#     charge = sum(signs .* sector)
#     _sector_is_allowed(charge, signs, sector) ||
#         error("sector is not allowed ", sector)
#     legs = U1Space[]
#     legs = Tuple([U1Space(signs[n], [sector[n]], [size(nzblock, n)])
#                   for n in eachindex(sector)])
#     SymTensor(charge, Tuple(legs), [sector], [nzblock])
# end

# """
#     findindexes_sorted(A, elements

# Find the index of a set of sorted elements in a sorted vector `A`.
# """
# function _trimlegs(legs::NTuple{N, U1Space}, sects::Vector{NTuple{N, Int}}) where{N}
#     legcharges = Tuple(BitSet() for i=1:N)
#     for sect in sects
#         for i in 1:N
#             push!(legcharges[i], sect[i])
#         end
#     end
#     #print(legcharges)
#     trimmedlegs = U1Space[]
#     for lidx in 1:N
#         leg = legs[lidx]
#         indexes = findindexes_sorted(leg.chrs, collect(legcharges[lidx]))
#         push!(trimmedlegs, U1Space(leg.sign, leg.chrs[indexes], leg.dims[indexes]))
#     end
#     # println("trimming :")
#     # println(legs)
#     # println(trimmedlegs)
#     Tuple(trimmedlegs)
# end
