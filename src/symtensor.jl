abstract type AbstractSymTensor{T<:Number, N} end
abstract type AbstractSymMatrix{T<:Number} <: AbstractSymTensor{T, 2} end

@inline numoflegs(::AbstractSymTensor{T, N}) where {T<:Number, N} = N
@inline charge(A::AbstractSymTensor) = A.charge

@inline size(A::AbstractSymTensor) = Tuple(fulldims(l) for l in A.legs)
@inline size(A::AbstractSymTensor, l::Int) = size(A)[l]

##TODO: should the below methods be "@inline" ?
@inline issimilar(A::T, B::T) where {T<:AbstractSymTensor} =
    A.charge == B.charge && A.legs == B.legs

@inline isequal(A::T, B::T) where {T<:AbstractSymTensor} =
    issimilar(A, B) &&  isequal(A.nzblks, B.nzblks)

@inline isapprox(A::T, B::T) where {T<:AbstractSymTensor} =
    issimilar(A, B) &&  isapprox(A.nzblks, B.nzblks)

@inline ==(A::T, B::T) where {T<:AbstractSymTensor} = isequal(A, B)

mutable struct SymTensor{T<:Number, N} <: AbstractSymTensor{T, N}
    charge :: Int
    legs   :: NTuple{N, STLeg}
    sects  :: Vector{NTuple{N, Int}}
    nzblks :: Vector{Array{T, N}}
    function SymTensor{T, N}(c,l,s,n) where {T<:Number, N}
        new{T, N}(c,l,s,n)
    end
end

function SymTensor(charge :: Int,
                   legs   :: NTuple{N, STLeg},
                   sects  :: Vector{NTuple{N, Int}},
                   nzblks :: Vector{<:AbstractArray{T, N}}) where{T<:Number, N}
    #issorted(sects, lt=_sectorlessthan) ||  error("sectors not sorted!")
    _allsectorsandsizes(charge, legs) == (sects, size.(nzblks)) ||
        error("The given sectors and block-sizes do not match legs!")#, charge, legs, sects, size.(nzblks))
    #length(sects) == length(nzblks) || error("sects don't match nzblks!")
    SymTensor{T, N}(charge, legs, sects, nzblks)
end

mutable struct SymMatrix{T<:Number} <: AbstractSymMatrix{T}
    charge :: Int
    legs   :: Tuple{STLeg, STLeg}
    sects  :: Vector{Tuple{Int, Int}}
    nzblks :: Vector{Matrix{T}}
    function SymMatrix{T}(c,l,s,n) where {T<:Number}
        new{T}(c,l,s,n)
    end
end

function SymMatrix(charge :: Int,
                   legs   :: Tuple{STLeg, STLeg},
                   sects  :: Vector{Tuple{Int, Int}},
                   nzblks :: Vector{Matrix{T}}) where {T<:Number}
    signs(legs) == (+1, -1) || error("SymMatrix signs should be (+1,-1) ", signs(legs))
    #issorted(sects, lt=_sectorlessthan) || error("sectors not sorted!")
    _allsectorsandsizes(charge, legs) == (sects, size.(nzblks)) ||
        error("The given sectors and block-sizes do not match legs!")
    #lenth(sects) == length(nzblks) || error("sects don't match nzblks!")
    SymMatrix{T}(charge, legs, sects, nzblks)
end

##TODO: see if this should be a convert function!
SymMatrix(A::AbstractSymMatrix) = A
function SymMatrix(A::SymTensor{T, 2}) where {T<:Number}
    signs(A.legs) == (+1, -1) ||
        error("SymTensor{T, 2} doesn't have correct signs to convert to SymMatrix")
    SymMatrix{T}(A.charge, A.legs, A.sects, A.nzblks)
end

mutable struct SymDiagonal{T<:Number} <: AbstractSymMatrix{T}
    charge :: Int
    legs   :: Tuple{STLeg, STLeg}
    sects  :: Vector{Tuple{Int, Int}}
    nzblks :: Vector{Diagonal{T,Vector{T}}}
    function SymDiagonal{T}(c,l,s,n) where {T<:Number}
        new{T}(c,l,s,n)
    end
end

function SymDiagonal(charge :: Int,
                     legs   :: Tuple{STLeg, STLeg},
                     sects  :: Vector{Tuple{Int, Int}},
                     nzblks :: Vector{Diagonal{T,Vector{T}}}) where {T<:Number}
    signs(legs) == (+1, -1) || error("SyDiagonal signs should be (+1,-1) ", signs(legs))
    #issorted(sects, lt=_sectorlessthan) || error("sectors not sorted!")
    _allsectorsandsizes(charge, legs) == (sects, size.(nzblks)) ||
        error("The given sectors and block-sizes do not match legs!")
    #lenth(sects) == length(nzblks) || error("sects don't match nzblks!")
    legs[1].dims == legs[2].dims || error("SymDiagonal is not square!")
    SymDiagonal{T}(charge, legs, sects, nzblks)
end

struct SymVector{T<:Number} <: AbstractSymTensor{T, 1}
    charge :: Int
    legs   :: Tuple{STLeg}
    sects  :: Vector{Tuple{Int}}
    nzblks :: Vector{Vector{T}}
    function SymVector{T}(c,l,s,n) where {T<:Number}
        new{T}(c,l,s,n)
    end
end

function SymVector(charge :: Int,
                   legs   :: Tuple{STLeg},
                   sects  :: Vector{Tuple{Int}},
                   nzblks :: Vector{Vector{T}}) where {T<:Number}
    length(sects) == length(nzblks) == 1 ||
        error("SymVector accepts only one sector or nzblk!")
    size(nzblks, 1) == legs[1].dims[1] || error("")
    legs[1].chrs = [charge]
    legs[1].dims = [size(nzblks[1], 1)]
    SymVector{T}(charge, legs, sects, nzblks)
end

SymVector(A::SymVector{T}) where{T<:Number} = A
function SymVector(sign::Int, charge::Int, v::Vector{T}) where{T<:Number}
    leg = STLeg(sign, [charge], [length(v)])
    SymVector{T}(charge, (leg,), [(charge,)], [v])
end

##TODO: see if this should be a convert function!
function SymVector(A::SymTensor{T, 1}) where {T<:Number}
    SymVector{T}(A.charge, A.legs, A.sects, A.nzblks)
end

@inline size(s::SymVector) = size(nzblks[1])

# ##TODO: make the below three into one
# function convert(::Type{SymTensor{T1, N}},
#                  A::SymTensor{T2, N}) where {T1<:Number, T2<:Number, N}
#     SymTensor(A.charge, A.legs, A.sects, convert.(Array{T2, N}, A.nzblks))
# end

# function convert(::Type{SymMatrix{T1}},
#                  A::SymMatrix{T2}) where {T1<:Number, T2<:Number}
#     SymMatrix(A.charge, A.legs, A.sects, convert.(Matrix{T2}, A.nzblks))
# end

# function convert(::Type{SymVector{T1}},
#                  A::SymVector{T2}) where {T1<:Number, T2<:Number}
#     SymVector(A.charge, A.legs, A.sects, convert.(Array{T2, N}, A.nzblks))
# end

@inline _sectors_sortperm(sects::Vector{NTuple{N, Int}};
                          by::Function=identity) where {N} =
                              sortperm(sects, by=by, lt=_sectorlessthan)

function index_sector(A::AbstractSymTensor{T, N},
                      s::NTuple{N, Int}) where{T<:Number, N}
    index = searchsortedfirst(A.sects, s, lt=_sectorlessthan)
    index > length(A.sects) && s != A.sects[index] && error("sector not found!")
    index
end

function get_sector(A::AbstractSymTensor{T, N},
                    s::NTuple{N, Int}) where {T<:Number, N}
    A.nzblks[index_sector(A, s)]
end

function set_sector!(A::AbstractSymTensor{T, N},
                     s::NTuple{N, Int},
                     nzblk::Array{T, N}) where {T<:Number, N}
    A.nzblks[index_sector(A, s)] = nzblk
    A
end

function rand(::Type{T}, charge::Int, legs::NTuple{N, STLeg};
              seed::Int=1911) where {T<:Number, N}
    sects, sizes = _allsectorsandsizes(charge, legs)
    rng = MersenneTwister(seed)
    nzblks = [rand(rng, T, dims) for dims in sizes]
    SymTensor(charge, legs, sects, nzblks)
end

rand(charge::Int, legs::NTuple{N, STLeg}) where {N} = rand(Float64, charge, legs)

function fill(x::T, charge::Int, legs::NTuple{N, STLeg}) where {T<:Number, N}
    sects, sizes = _allsectorsandsizes(charge, legs)
    nzblks = [fill(x, dims) for dims in sizes]
    SymTensor(charge, legs, sects, nzblks)
end

function fill!(A::SymTensor{T, N}, x::T) where {T<:Number, N}
    for i in eachindex(A.nzblks)
        fill!(A.nzblks[i], x)
    end
    A
end

##TODO: Mix this with the usual UnitScaling
function eye(::Type{T}, chrs::Vector{Int}, dims::Vector{Int}) where {T<:Number}
    l1 = STLeg(+1, chrs, dims)
    l2 = STLeg(-1, chrs, dims)

    sects, sizes = _allsectorsandsizes(zero(Int), (l1,l2))
    nzblks = [Matrix{T}(I, dims) for dims in sizes]
    SymMatrix(zero(Int), (l1,l2), sects, nzblks)
end

eye(cs::Vector{Int}, ds::Vector{Int}) = eye(Float64, cs, ds)

"""
    invlegs(sten)

invert the direction (sign) of all legs of the tensor therefore
negative the total charge of the tensor as well!
"""
# This function seems to be not respecting the (+1, -1) convention for
# SymMatrix objects, should I change the convention then?
function invlegs(A::AbstractSymTensor)
    legs = Tuple([STLeg(-l.sign, l.chrs, l.dims) for l in A.legs])
    typeof(A)(-A.charge, legs, A.sects, A.nzblks)
end

function mapcharges(f::Function, A::AbstractSymTensor)
    SymTensor(f(A.charge),
              mapcharges.(f, A.legs),
              [f.(s) for s in A.sects],
              A.nzblks)
end

function mapcharges(f::NTuple{N, Function},
                    A::AbstractSymTensor{T, N}) where{T<:Number, N}
    legs = Tuple(mapcharges(f[i], A.legs[i]) for i in 1:N)
    sects = [Tuple(f[i](s[i]) for i in 1:N) for s in A.sects]
    sgns = signs(legs)
    charges = [sum(sgns .* s) for s in sects]
    !all(charges .== charges[1]) && error("mapcharges function is inconsistent!")
    SymTensor(charges[1], legs, sects, A.nzblks)
end

conj(A::AbstractSymTensor) =
    typeof(A)(A.charge, A.legs, A.sects, [conj(blk) for blk in A.nzblks])

function array(A::AbstractSymTensor)
    arrep = zeros(eltype(A), sum.(alldims(A.legs)))
    adims = Tuple([0;s] for s in accdims(A.legs))
    chrs = [0, 1]
    for idx in eachindex(A.sects)
        sect = A.sects[idx]
        nzblk = A.nzblks[idx]
        ranges = []
        for n in eachindex(sect)
            i = findfirst(sect[n].== A.legs[n].chrs)
            push!(ranges, adims[n][i]+1:adims[n][i+1])
        end
        arrep[Tuple(ranges)...] = nzblk
    end
    arrep
end


@inline *(A::AbstractSymTensor, a::T) where {T<:Number} =
    SymTensor(A.charge, A.legs, A.sects, [a .* blk for blk in A.nzblks])
@inline *(a::T, A::AbstractSymTensor) where {T<:Number} = *(A, a)

function removedummyleg(A::AbstractSymTensor, l::Int)
    isdummy(A.legs[l]) || error("leg is not dummy!")

    N = numoflegs(A)
    SymTensor(A.charge,
              (A.legs[1:l-1]..., A.legs[l+1:N]...),
              [(s[1:l-1]...,s[l+1:N]...) for s in A.sects],
              [reshape(blk, size(blk)[1:l-1]...,size(blk)[l+1:N]...) for blk in A.nzblks])
end

##TODO: this stuff should go up
@inline eltype(::AbstractSymTensor{T}) where {T<:Number} = T
@inline eltype(::Type{<:AbstractSymTensor{T}}) where {T<:Number} = T

"construct an additional SymTensor similar to A, possibly with a
different scalar type T."
function similar(A::AbstractSymTensor, T::Type=eltype(A))
    typeof(A)(A.charge, A.legs, A.sects, [similar(blk, T) for blk in A.nzblks])
end

"copy the contents of A to a preallocated SymTensor B"
function copyto!(B::T, A::T) where {T<:AbstractSymTensor}
    issimilar(B, A) && error("Can not copy to a  non-similar SymTensor!")
    for i in eachindex(A.nzblks)
        copyto!(B.nzblks[i], A.nzblks[i])
    end
    B
end

"out of place scalar multiplication; multiply SymTensor A with scalar α
and store the result in B"
function mul!(B::T, A::T, α) where {T<:AbstractSymTensor}
    B.nzblks = [α .* blk for blk in A.nzblks]
    B
end

"in-place scalar multiplication of A with α; in particular with α =
false, A is initialized with all zeros"
function rmul!(A::T, α) where {T<:AbstractSymTensor}
    for i in eachindex(A.nzblks)
        rmul!(A.nzblks[i],  α)
    end
    A
end

"Stores in B the result of α*A + B"
function axpy!(α,
               A::AbstractSymTensor{T1, N},
               B::AbstractSymTensor{T2, N}) where {T1<:Number, T2<:Number, N}
    issimilar(A, B) || error("axpy! The two matrices are not simliar!")
    for i in eachindex(A.nzblks)
        B.nzblks[i] = α .* A.nzblks[i] + B.nzblks[i]
    end
    B
end

"Stores in B the result of α*A + β*B"
function axpby!(α,
                A::AbstractSymTensor{T1, N},
                β,
                B::AbstractSymTensor{T2, N}) where {T1<:Number, T2<:Number, N}
    issimilar(A, B) || error("axpy! The two matrices are not simliar!")
    for i in eachindex(A.nzblks)
        B.nzblks[i] = α .* A.nzblks[i] + B.nzblks[i]
    end
    B
end

"""
compute the inner product of two similar SymTensors which conjugates
the second one
"""
function dot(A::AbstractSymTensor{T1, N},
             B::AbstractSymTensor{T2, N}) where{T1<:Number, T2<:Number, N}
    #contract(A, .-Tuple(1:N), invlegs(conj(B)), .-Tuple(1:N))
    issimilar(A, B) || error("The two SymTensors have to be similar to dot!")
    sum([dot(A.nzblks[i], B.nzblks[i]) for i in eachindex(A.nzblks)])
end

" compute the 2-norm of a  AbstractSymTensor"
function norm(A::AbstractSymTensor)
    sqrt(dot(A, A))
end

function normalize!(S::SymDiagonal)
    s = norm(S)
    for i in eachindex(S.nzblks)
        rmul!(S.nzblks[i], 1/s)
    end
    S
end

function show(A::AbstractSymTensor)
    num_ch = length(A.legs)
    println("$(size(A)) $(typeof(A)) with $(length(A.sects)) blocks: ")
    for i in eachindex(A.sects)
        print("sector: $(A.sects[i]) size: ")
        show(stdout, "text/plain", A.nzblks[i])
        println()
    end
    nothing
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
#     legs = STLeg[]
#     legs = Tuple([STLeg(signs[n], [sector[n]], [size(nzblock, n)])
#                   for n in eachindex(sector)])
#     SymTensor(charge, Tuple(legs), [sector], [nzblock])
# end

# """
#     findindexes_sorted(A, elements)

# Find the index of a set of sorted elements in a sorted vector `A`.
# """
# function findindexes_sorted(A::Vector{T}, elements::Vector{T}) where{T}
#     #println(A, elements)
#     p = 1
#     indexes = Int[]
#     for i in eachindex(A)
#         if A[i] == elements[p]
#             push!(indexes, i)
#             if p==length(elements) break end
#             p += 1
#         elseif A[i] < elements[p]
#             error("findindexes_sorted ", A, elements)
#         end
#     end
#     indexes
# end

# function _trimlegs(legs::NTuple{N, STLeg}, sects::Vector{NTuple{N, Int}}) where{N}
#     legcharges = Tuple(BitSet() for i=1:N)
#     for sect in sects
#         for i in 1:N
#             push!(legcharges[i], sect[i])
#         end
#     end
#     #print(legcharges)
#     trimmedlegs = STLeg[]
#     for lidx in 1:N
#         leg = legs[lidx]
#         indexes = findindexes_sorted(leg.chrs, collect(legcharges[lidx]))
#         push!(trimmedlegs, STLeg(leg.sign, leg.chrs[indexes], leg.dims[indexes]))
#     end
#     # println("trimming :")
#     # println(legs)
#     # println(trimmedlegs)
#     Tuple(trimmedlegs)
# end

function negateleg(A::AbstractSymTensor, l::Int)
    N = numoflegs(A)
    0 < l <= N || error("leg $l doesn't exist!")
    legs = (A.legs[1:l-1]..., negate(A.legs[l]), A.legs[l+1:N]...)

    n_sectors = length(A.sects)
    sects = Vector{NTuple{N, Int}}(undef, n_sectors)
    for i = 1:n_sectors
        sect = A.sects[i]
        sects[i] = (sect[1:l-1]..., -sect[l], sect[l+1:N]...)
    end
    perm = _sectors_sortperm(sects)
    SymTensor{eltype(A),N}(A.charge, legs, sects[perm], A.nzblks[perm])
end
