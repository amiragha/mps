### A bunch of helper functions to apply things to mps
function _applymps0site(Λ::Matrix{T},
                        envL::Array{T, 3},
                        envR::Array{T, 3}) where {T<:Number}
    @tensor A[d,d'] := envL[u,m,d] * Λ[u, u'] * envR[u',m,d']
    A
end
