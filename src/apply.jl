### A bunch of helper functions to apply stuff with environments to mps
function _applymps0site(Λ::Matrix{T},
                        envL::Array{T, 3},
                        envR::Array{T, 3}) where {T<:Number}
    @tensor A[d,d'] := envL[u,m,d] * Λ[u, u'] * envR[u',m,d']
    A
end

function _mpsupdateright(renv::Array{T, 3},
                         mat::Array{T, 3},
                         hmpo::Array{T,4}) where {T<:Number}
    @tensor R[u,m,d] :=
        ((renv[u',m',d'] * mat[u,o,u']) * hmpo[m,o',m',o]) * conj(mat)[d,o',d']
    R
end

function _mpsupdateright(renv::SymTensor{S,T,3},
                         mat::SymTensor{S,T,3},
                         hmpo::SymTensor{S,T,4}) where {S,T}

    R = contract(contract(contract(renv, (-1,3,4), mat, (1,2,-1)),
                          (1, -1, -2, 4), hmpo, (3,2,-2,-1)),
                 (1,-1,2,-2), dual(mat), (3,-1,-2))
    R
end

# function _mpsupdateright(renv::SymTensor{Tv, 4},
#                          mat::SymTensor{Tv, 3},
#                          hmpo::SymTensor{Tv, 4}) where {Tv<:Number}

#     # println(renv.legs)
#     # println(mat.legs)
#     # println(hmpo.legs)
#     # renv is W_l x M_d x W_r x M_u
#     R = contract(contract(contract(renv, (1,2,3,-1), mat, (4,5,-1)),
#                           (-1,3,4,5,-2), hmpo, (1,2,-1,-2)),
#                  (1, -1,-2,3,4), invlegs(conj(mat)), (2, -1,-2))
#     R
# end

function _mpsupdateleft(lenv::Array{T, 3},
                        mat::Array{T, 3},
                        hmpo::Array{T,4}) where {T<:Number}
    @tensor L[u,m,d] :=
        ((lenv[u',m',d'] * mat[u',o,u]) * hmpo[m',o',m,o]) * conj(mat)[d',o',d]
    L
end

function _mpsupdateleft(lenv::SymTensor{S,T, 3},
                        mat::SymTensor{S,T, 3},
                        hmpo::SymTensor{S,T, 4}) where {S,T}

    L = contract(contract(contract(lenv, (-1,3,4), mat, (-1,2,1)),
                          (1, -1, -2, 4), hmpo, (-2,2,3,-1)),
                 (1,-1,2,-2), dual(mat), (-2,-1,3))

    L
end

# function _mpsupdateleft(lenv::SymTensor{Tv, 4},
#                         mat::SymTensor{Tv, 3},
#                         hmpo::SymTensor{Tv, 4}) where {Tv<:Number}

#     L = contract(contract(contract(lenv, (1,2,3,-1), mat, (-1,4,5)),
#                           (1,2,-1,-2,5), hmpo, (-1,3,4,-2)),
#                  (1,-1,-2,3,4), invlegs(conj(mat)), (-1,-2,2))
#     L
# end

function _applymps1site(v,
                        envL::Array{T,3},
                        envR::Array{T,3},
                        hmpo::Array{T,4}) where {T}

    @tensor v[l,o,r] := ((envL[l',ml,l] * v[l',o',r']) *
                         hmpo[ml,o,mr,o']) * envR[r',mr,r]
    v
end

function _applymps2site(v,
                        envL::Array{T,3},
                        envR::Array{T,3},
                        hmpoL::Array{T,4},
                        hmpoR::Array{T,4}) where {T}

    @tensor v[l,o1,o2,r] := (((envL[l',ml,l] * v[l',o1',o2',r']) *
                              hmpoL[ml,o1,mm,o1']) *
                             hmpoR[mm,o2,mr,o2']) *
                             envR[r',mr,r]
    v
end

function _applymps2site(v::SymTensor{S,T,4},
                        envL::SymTensor{S,T,3},
                        envR::SymTensor{S,T,3},
                        hmpoL::SymTensor{S,T,4},
                        hmpoR::SymTensor{S,T,4}) where {S,T}

    v = contract(contract(contract(contract(
        envL, (-1,2,1), v, (-1,3,4,5)),
                                   (1,-1,-2,4,5), hmpoL, (-1,2,3,-2)),
                          (1,2,-1,-2,5), hmpoR, (-1,3,4,-2)),
                 (1,2,3,-1,-2), envR, (-2,-1,4))
    v
end

# function _applymps2site(v::SymTensor{Tv,4},
#                         envL::SymTensor{Tv,4},
#                         envR::SymTensor{Tv,4},
#                         hmpoL::SymTensor{Tv,4},
#                         hmpoR::SymTensor{Tv,4}) where {Tv<:Number}

#     v = contract(contract(contract(contract(
#         envL, (1,2,3,-1), v, (-1,4,5,6)),
#                                    (1,2,-1,-2,5,6), hmpoL, (-1,3,4,-2)),
#                           (1,2,3,-1,-2,6), hmpoR, (-1,4,5,-2)),
#                  (-1,1,2,3,-2,-3), envR, (-2,4,-1,-3))
#     v
# end
