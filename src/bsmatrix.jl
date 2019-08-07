struct STLeg
    sign :: Int    # +1 for ket (ingoing) or -1 for bra (outgoing)
    chrs :: Vector{Charge}    # set of possible charges
    dims :: Vector{Int}    # set of dimensions for possible charges

    function STleg(sign, chrs, dims)
        abs(sign) == 1 || error("unacceptable sign : ", sign)
        all(dims .=> 0) || error("space dimensions can only be non-negative!")
        length(chrs) == length(dims) ||
            error("Number of charges and dimensions don't match!")
        issorted(chrs) || error("charges are not sorted")
        new(sign, chrs, dims)
    end
end

# struct Z2Leg
#     sign :: Int    # +1 for ket (ingoing) or -1 for bra (outgoing)
#     dims :: Vector{Int}    # set of dimensions for possible charges

#     function Z2Leg(sign, dims)
#         abs(sign) == 1 || error("unacceptable sign : ", sign)
#         length(dims) == 2 || error("Z2 has only two possible charges! ", dims)
#         all(dims .>= 0) || error("Dimensions should be non-negative! ", dims)
#         new(sign, dims)
#     end
# end

@inline charges(legs::NTuple{N, STLeg}) where {N} =
    Tuple(legs[n].chrs for n in eachindex(legs))

@inline signs(legs::NTuple{N, STLeg}) where {N} =
    Tuple(legs[n].sign for n in eachindex(legs))

@inline alldims(legs::NTuple{N, STLeg}) where {N} =
    Tuple(legs[n].dims for n in eachindex(legs))

@inline accdims(legs::NTuple{N, STLeg}) where {N} =
    Tuple(cumsum(legs[n].dims) for n in eachindex(legs))

# @inline signs(legs::NTuple{N, Z2Leg}) where {N} =
#     Tuple(legs[n].sign for n in eachindex(legs))

# @inline alldims(legs::NTuple{N, Z2Leg}) where {N} =
#     Tuple(legs[n].dims for n in eachindex(legs))

# @inline accdims(legs::NTuple{N, Z2Leg}) where {N} =
#     Tuple(cumsum(legs[n].dims) for n in eachindex(legs))


# # simplest possible symmetry matrix, has only charges {+1, -1}
# mutable struct Z2Matrix{Tv}
#     rowdims :: Tuple{Int, Int}
#     coldims :: Tuple{Int, Int}
#     nzblocks :: Vector{Matrix{Tv}}
# end

# struct Z2Tensor{Tv<:Number, N}
#     legs :: NTuple{N, Z2Leg}
#     sects :: Vector{NTuple{N, Int}}    # possible nonzero sectors
#     nzblks :: Vector{Array{Tv, N}}     # stored nonzero blocks

#     function Z2Tensor(legs::NTuple{N, Z2Leg},
#                       sects::Vector{NTuple{N, Int}},
#                       nzblks::Vector{Array{Tv, N}}) where{Tv<:Number, N}
#         issorted(sects, lt=_sector_less_than) ||  error("sectors not sorted!")
#         new{Tv, N}(legs, sects, nzblks)
#     end
# end

struct SymTensor{Tv<:Number, N}
    legs :: NTuple{N, Z2Leg}
    sects :: Vector{NTuple{N, Int}}    # possible nonzero sectors
    nzblks :: Vector{Array{Tv, N}}     # stored nonzero blocks

    function SymTensor(legs::NTuple{N, Z2Leg},
                      sects::Vector{NTuple{N, Int}},
                      nzblks::Vector{Array{Tv, N}}) where{Tv<:Number, N}
        issorted(sects, lt=_sector_less_than) ||  error("sectors not sorted!")
        new{Tv, N}(legs, sects, nzblks)
    end
end



"""
    Z2Tensor(signs, sector, nzblock)

Make Z2Tensor with the given possible charges and only one nonzero
sector. Find the dimensions from that.

"""
function Z2Tensor(signs::NTuple{N, Int},
                  sector::NTuple{N, Int},
                  nzblock::Array{Tv, N}) where {Tv<:Number, N}

    chrs = [0, 1]
    _sector_is_allowed(0, signs, sector, [0, 1]) ||
        error("sector is not allowed ", sector)
    legs = Z2Leg[]
    for n in eachindex(sector)
        dims = zeros(Int, length(chrs))
        dims[findfirst(sector[n] .== chrs)] = size(nzblock, n)
        push!(legs, Z2Leg(signs[n], dims))
    end

    return Z2Tensor{Tv, N}(Tuple(legs), [sector], [nzblock])
end

@inline size(bst::Z2Tensor) =
    Tuple(sum(bst.legs[n].dims) for n in eachindex(bst.legs))


# we always sort sectors based on charges (ignore signs)
### NOTE the sorting is column major. That means the leftmost (first)
### leg changes charge faster (first)!
function _sector_less_than(sect1::NTuple{N, Int},
                           sect2::NTuple{N, Int}
                           #,signs::NTuple{N, Int}
                           ) where {N, Tc<:Integer}
    for n in reverse(eachindex(sect1))
        #sign = signs[n]
        #a, b = sign*sect1[n], sign*sect2[n]
        a, b = sect1[n], sect2[n]
        if a < b
            return true
        elseif a > b
            return false
        end
    end
    false
end

@inline _sectors_sortperm(sects::Vector{NTuple{N, Int}};
                          by::Function=identity) where {N} =
    sortperm(sects, by=by, lt=_sector_less_than)

function _sector_is_allowed(total::Int,
                            signs::NTuple{N},
                            sector::NTuple{N, Int},
                            chrs::Vector{Int}) where {N}
    #@assert length(sector) == length(chrs)
    for n in eachindex(sector)
        if !(sector[n] in chrs)
            return false
        end
    end
    return sum(signs .* sector) % 2 == total
end

"""
    fuse_legs(legs, parts)

fuse a sect of `legs` based on the partitions given by `parts` with new
signs given by `parts_signs`.

"""
function fuse_legs(legs::NTuple{N, Z2Leg},
                   parts::Vector{Int},
                   part_signs::Vector{Int}) where {N, Tc<:Integer}
    @assert sum(parts) == N
    @assert abs.(part_sgns) == ones(Int, length(parts))

    chrs = charges(legs)
    dims = alldims(legs)
    sects, sect_sizes = _possible_charge_sectors_and_sizes(legs)

    # fuse the sectors and find the fused info
    flegs_info = [Tuple{Tc, Vector{Tc}, Int}[] for n in 1:length(parts)]
    sector_map = Vector{Int}[]
    for i in eachindex(sects)
        sect = sects[i]
        sect_size =sect_sizes[i]
        fsect = fuse_set(+, Tuple(sect), parts)
        fsize = fuse_set(*, Tuple(sect_size), parts)

        # the information for the new legs are stored in
        pivot = 0
        map = zeros(Int, length(parts))

        fcharges = [Tc[] for n in length(parts)]
        fchcombos = [Vector{Tc}[] for n in length(parts)]
        fchsizes = [Int[] for n in length(parts)]

        for n in eachindex(fsect)
            fcharge  = fsect[n]
            fchcombo = sect[pivot+1:pivot+parts[n]]
            fchsize  = fsize[n]

            #info = (fcharge, fchcombo, fchsize)

            cidx = findfirst(fcharges, fcharge)

            if didx > 0
                # this charge already exists
                combo_idx = findfirst(fchcombos[cidx], fchcombo)
                if combo_idx > 0
                    # the combination for this charge also exist
                    fchsizes[cdix][combo_idx]
                else
                    # the combination is new
                    # Add combo and size and save the map
                end
            else
                # the charge is new

            end

            idx = searchsortedfirst(flegs_info[n], info)
            if length(flegs_info[n]) > idx && flegs_info[n][idx] == info
                map[n] = index
            else
                push!(flegs_info[n], info)
                map[n] = length(flegs_info[n])
            end
            pivot += parts[n]
        end
        push!(sector_map, map)
    end

    return flegs_info, sector_map
end

function symtensor2array(ten::Z2Tensor{T}) where {T<:Number}
    arrep = zeros(sum.(alldims(ten.legs)))
    adims = Tuple([0;s] for s in accdims(ten.legs))
    println(adims)
    chrs = [0, 1]
    for idx in eachindex(ten.sects)
        sect = ten.sects[idx]
        nzblk = ten.nzblks[idx]
        ranges = []
        for n in eachindex(sect)
            i = findfirst(sect[n].== chrs)
            push!(ranges, adims[n][i]+1:adims[n][i+1])
        end
        arrep[Tuple(ranges)...] = nzblk
    end
    arrep
end

function show(ten::Z2Tensor)
    num_ch = length(ten.legs)
    #m, n = size(bsa)
    println(size(ten), typeof(ten), " with ", length(ten.sects), " blocks: ")
    for i in eachindex(ten.sects)
        print("sector: ", ten.sects[i], ", dims: ")
        show(stdout, "text/plain", ten.nzblks[i])
        println()
    end
    nothing
end
