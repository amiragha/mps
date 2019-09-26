abstract type AbstractCharge end

# mutable struct FusedCharge{N}
#     charge :: Int # total charge
#     dim :: Int # total dimension
#     #signs :: NTuple{N, Int}
#     pats :: Dict{NTuple{N, Int}, UnitRange{Int}} # patterns and positions

#     function FusedCharge(charge, dim, pats::Dict{NTuple{N, Int}, UnitRange{Int}}) where {N}
#         all(sum.(keys(pats)) .== charge) || error("oops", pats, charge)
#         #all(abs.(signs) .== 1) || error("oops")
#         new{N}(charge, dim, pats)
#     end
# end

mutable struct FusedCharge{N}
    charge :: Int # total charge
    dim :: Int # total dimension

    pats :: Vector{NTuple{N, Int}}
    ranges :: Vector{UnitRange{Int}}

    function FusedCharge(charge, dim, pats::Vector{NTuple{N, Int}}, ranges::UnitRange{Int}) where {N}
        #all(sum.(keys(pats)) .== charge) || error("oops", pats, charge)
        #all(abs.(signs) .== 1) || error("oops")
        new{N}(charge, dim, pats)
    end
end

# struct DefusedCharge{N}
#     charge :: Int
#     dim :: Int

# end
