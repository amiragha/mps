abstract type AbstractSymTensor{S, T, N} end
const AbstractSymMatrix{S, T} = AbstractSymTensor{S, T, 2}

"""
A symmetric tensor is represented as a tensor map over the space of
V_row (codomain) with tensor product structure of (((V1 ⊗ V2) ⊗ V3) ⊗
...) and the space of V_col (domain) with tensor product structure of
(...⊗ (V_{n-3} ⊗ (V_{n-1} ⊗ V_n))). The data is also represented as a
matrix
"""

mutable struct SymTensor{S, T, N} <: AbstractSymTensor{S, T, N}
    charge   :: S
    space    :: NTuple{N, VectorSpace{S}}
    data     :: SortedDict{Sector{S, N}, Array{T, N}}

    function SymTensor(charge :: S,
                       space  :: NTuple{N, VectorSpace{S}},
                       data   :: SortedDict{Sector{S, N}, Array{T, N}}) where{S,T,N}
        sects, sizes = _allsectorsandsizes(charge, space)
        i = 1
        for (s,d) in data
            s == sects[i] || throw("SectorMismatch(), $s, $sects[i]")
            size(d) == sizes[i] || throw(SizeMismatch())
            i += 1
        end
        new{S,T,N}(charge, space, data)
    end
    SymTensor{S,T,N}(c,s,d) where{S,T,N}= SymTensor(c,s,d)
end

# function SymTensor{S,T1,N}(A::SymTensor{S,T2,N}) where {S,T1, T2, N}
#     SymTensor{S,T1,N}(A.charge, A.space, SortedDict{Sector{S,N}, T1}(A.data))
# end

const SymMatrix{S, T} = SymTensor{S, T, 2}
const SymVector{S, T} = SymTensor{S, T, 1}

const U1Tensor{T, N} = SymTensor{Int, T, N}
const U1Matrix{T} = SymMatrix{Int, T}
const U1Vector{T} = SymVector{Int, T}

mutable struct SymDiagonal{S, T<:Number} <: AbstractSymMatrix{S, T}
    charge :: S
    space  :: NTuple{2, VectorSpace{S}}
    data   :: SortedDict{Sector{S, 2}, Diagonal{T}}

    function SymDiagonal(charge :: S,
                         space  :: NTuple{2, VectorSpace{S}},
                         data   :: SortedDict{Sector{S,2}, Diagonal{T}}) where{S,T}
        sects, sizes = _allsectorsandsizes(charge, space)
        i = 1
        for (s,d) in data
            s == sects[i] || throw{SectorMismatch()}
            size(d) == sizes[i] || throw(SizeMismatch())
            i += 1
        end
        chargedims(first(space)) == chargedims(last(space)) ||
            throw("SymDiagonal has to be square")
        new{S,T}(charge, space, data)
    end
end

const U1Diagonal{T} = SymDiagonal{Int, T}

function SymVector(charge::S, v::Vector{T}) where{S, T}
    V = VectorSpace{S}(charge => length(v))
    SymVector{S, T}(charge, (V,), SortedDict(Sector(charge) => v))
end

function convert(::Type{SymTensor{S, T1, N}},
                 A::SymTensor{S, T2, N}) where {S, T1, T2, N}
    SymTensor{S, T1, N}(A.charge, A.space, A.data)
end

@inline rank(::AbstractSymTensor{S, T, N}) where {S, T, N} = N
@inline charge(A::AbstractSymTensor) = A.charge

@inline size(A::AbstractSymTensor) = Tuple(dim(V) for V in A.space)
@inline size(A::AbstractSymTensor, l::Int) = dim(A.space[l])

@inline issimilar(A::T, B::T) where {T<:AbstractSymTensor} =
    A.charge == B.charge && A.space == B.space

@inline isequal(A::T, B::T) where {T<:AbstractSymTensor} =
    issimilar(A, B) && isequal(A.data, B.data)

@inline isapprox(A::AbstractSymTensor, B::AbstractSymTensor) =
    issimilar(A, B) && isapprox(A.nzblks, B.nzblks)

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

@inline get_sector(A::AbstractSymTensor{T, N},
                   S::U1Sector{N}) where {T<:Number, N} = get(A.data, S)
function set_sector!(A::AbstractSymTensor{T, N},
                     S::U1Sector{N},
                     blk::Array{T, N}) where {T<:Number, N}
    A.data[S] = blk
    A
end

function rand(::Type{T}, charge::Int, space::NTuple{N, U1Space};
              seed::Int=1911) where {T<:Number, N}
    sects, sizes = _allsectorsandsizes(charge, space)
    rng = MersenneTwister(seed)
    data = SortedDict(sects[i] => rand(rng, T, sizes[i]) for i=1:length(sects))
    SymTensor(charge, space, data)
end

rand(charge::Int, space::NTuple{N, U1Space}) where {N} = rand(Float64, charge, space)

function fill(x::T, charge::Int, space::NTuple{N, VectorSpace{S}}) where {S,T<:Number,N}
    sects, sizes = _allsectorsandsizes(charge, space)
    data = SortedDict(sects[i] => fill(x, sizes[i]) for i=1:length(sects))
    SymTensor(charge, space, data)
end

function fill!(A::SymTensor{S,T,N}, x::T) where {S,T,N}
    for i in eachindex(A.nzblks)
        fill!(A.nzblks[i], x)
    end
    A
end

function fill_linearindex(charge::Int, space::NTuple{N, U1Space}) where {N}
    sects, sizes = _allsectorsandsizes(charge, legs)
    nzblks = Vector{Array{Int, N}}()
    p = 0
    for dims in sizes
        push!(nzblks, reshape(collect(p+1:p+prod(dims)), dims...))
        p += prod(dims)
    end
    SymTensor(charge, legs, SectorArray(sects[i]=>nzblks[i] for i=1:length(sects)))
end

##TODO: Mix this with the usual UnitScaling
function eye(::Type{T}, chrs::Vector{Int}, dims::Vector{Int}) where {T<:Number}
    V1 = U1Space(chrs, dims)
    V2 = U1Space(-1 .* chrs, dims)

    sects, sizes = _allsectorsandsizes(zero(Int), (V1,V2))
    SymMatrix(zero(Int), (l1,l2), sects,
              SectorArray(sects[i] => I(sizes[i]) for i=1:length(sects)))
end

eye(cs::Vector{Int}, ds::Vector{Int}) = eye(Float64, cs, ds)

"""
    tensordual(A)

Change the direction (sign) of all legs of the tensor to make the dual
of tensor. This operation negates the: tensor charge, leg signs,
sectors as well! One can set whether or not it also conjugates the
data with `conjugate` which is true by default.

"""
# This function seems to be not respecting the (+1, -1) convention for
# SymMatrix objects, should I change the convention then?
function dual(A::AbstractSymTensor; conjugate::Bool=true)
    if !conjugate
        return typeof(A)(-A.charge,
                         dual.(A.space),
                         typeof(A.data)(dual(s), d for (s,d) in A.data))
    end
    typeof(A)(-A.charge,
              dual.(A.space),
              typeof(A.data)(dual(s), conj(d) for (s,d) in A.data))
end

function mapcharges(f::Function, A::AbstractSymTensor)
    SymTensor(f(A.charge),
              mapcharges.(f, A.space),
              typeof(A.data)(mapcharge(f, s),d for (s,d) in A.data))
end

function mapcharges(f::NTuple{N, Function},
                    A::AbstractSymTensor{T, N}) where{T<:Number, N}
    space = Tuple(mapcharges(f[i], A.space[i]) for i in 1:N)

    sects = [Tuple(f[i](s[i]) for i in 1:N) for s in A.sects]
    sgns = signs(legs)
    charges = [sum(s) for s in sectors(A.data)]
    !all(charges .== charges[1]) &&
        error("mapcharges functions are inconsistent!")
    SymTensor(charges[1], space,
              typeof(A.data)(mapcharge(f, s),d for (s,d) in A.data))
end

conj(A::AbstractSymTensor) =
    typeof(A)(A.charge, A.space, conj(A.data))

#TODO: to make this work!
function array(A::AbstractSymTensor)
    arrep = zeros(eltype(A), dims.(A.space)...)
    adims = Tuple([0;s] for s in accdims(A.space))
    chrs = [0, 1]
    for idx in eachindex(A.sects)
        sect = A.sects[idx]
        nzblk = A.nzblks[idx]
        ranges = []
        for n in eachindex(sect)
            i = findfirst(isequal(legs[n].sign*sect[n]), A.legs[n].chrs)
            push!(ranges, adims[n][i]+1:adims[n][i+1])
        end
        arrep[Tuple(ranges)...] = nzblk
    end
    arrep
end


@inline *(A::AbstractSymTensor, a::T) where {T<:Number} =
    typeof(A)(A.charge, A.space, *(A.data, a))

@inline *(a::T, A::AbstractSymTensor) where {T<:Number} = *(A, a)

"construct an additional SymTensor similar to A, possibly with a
different scalar type T."
function similar(A::AbstractSymTensor, T::Type=eltype(A))
    typeof(A)(A.charge, A.legs, typeof(A)((s, similar(d, T)) for (s,d) in A.data))
end

"copy the contents of A to a preallocated SymTensor B"
function copyto!(B::T, A::T) where {T<:AbstractSymTensor}
    issimilar(B, A) && error("Can not copy to a  non-similar SymTensor!")
    for i in eachindex(B.data.values)
        copyto!(B.data.values[i], A.data.values[i])
    end
    B
end

"out of place scalar multiplication; multiply SymTensor A with scalar α
and store the result in B"
function mul!(B::T, A::T, α) where {T<:AbstractSymTensor}
    for i in eachindex(B.data.values)
        mul!(B.data.values[i], A.data.values[i], α)
    end
    B
end

"in-place scalar multiplication of A with α; in particular with α =
false, A is initialized with all zeros"
function rmul!(A::T, α) where {T<:AbstractSymTensor}
    for i in eachindex(A.data.values)
        rmul!(A.data.values[i],  α)
    end
    A
end

"Stores in B the result of α*A + B"
function axpy!(α,
               A::AbstractSymTensor{T1, N},
               B::AbstractSymTensor{T2, N}) where {T1<:Number, T2<:Number, N}
    issimilar(A, B) || error("axpy! The two matrices are not simliar!")
    for i in eachindex(A.data.values)
        B.data.values[i] = α .* A.data.values[i] + B.data.values[i]
    end
    B
end

"Stores in B the result of α*A + β*B"
function axpby!(α,
                A::AbstractSymTensor{T1, N},
                β,
                B::AbstractSymTensor{T2, N}) where {T1<:Number, T2<:Number, N}
    issimilar(A, B) || error("axpy! The two matrices are not simliar!")
    for i in eachindex(A.data.values)
        B.data.values[i] = α .* A.data.values[i] + β .* B.data.values[i]
    end
    B
end

"""
compute the inner product of two similar SymTensors which conjugates
the second one
"""
function dot(A::AbstractSymTensor{T1, N},
             B::AbstractSymTensor{T2, N}) where{T1<:Number, T2<:Number, N}
    issimilar(A, B) || error("The two SymTensors have to be similar to dot!")
    sum([dot(A.data.values[i], B.data.values[i]) for i in eachindex(A.data.values)])
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
#               typeof(A.data)((s[1:l-1]...,s[l+1:N]...),
#     reshape(d, size(d)[1:l-1]...,size(d)[l+1:N]...) for (s,d) in A.data)
#               end
#               end

function show(A::AbstractSymTensor)
    num_ch = length(A.space)
    println("$(size(A)) $(typeof(A)) with $(length(A.sects)) blocks: ")
    for i in eachindex(A.data.values)
        print("sector: $(A.data.keys[i]) size: ")
        show(stdout, "text/plain", A.data.values[i])
        println()
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
#     findindexes_sorted(A, elements)

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
