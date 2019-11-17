"""
    zipandgutzwiller_exact(mps1, mps2, mode)

given two MPS, with the same number of sites and 2 physical
dimensions perform the exact gutzwiller projection of the two. The
final result is an MPS with 2 physical dimension where the final
matrices are the tensor product of the given ones according to the
`mode`.

This function is just for testing. It is very inefficient because it
generates the tensor product matrices first and then applies svd on
them, the function to use is `zipandgutzwiller`

"""
function zipandgutzwiller_exact(mps1::MatrixProductState{T},
                                mps2::MatrixProductState{T};
                                mode::Symbol=:B14,
                                maxdim::Int64=200) where {T<:RLorCX}
    @assert mps1.d = mps2.d = 2
    lx = mps1.length
    @assert mps2.length = lx
    center = mps1.center
    @assert mps2.center = center

    # what should I do here exactly?!
    # what do I mean by exact?

end

"""
    zipandgutzwiller(mps1, mps2, mode)

given two MPS, with the same number of sites and two physical
dimensions perform gutzwiller projection of the two by making the
matrices using iterative SVD. The final result is an MPS with 2
physical dimension where the final matrices are the tensor product of
the given ones according to the `mode`.

"""

function zipandgutzwiller!(mps1::MatrixProductState{T},
                           mps2::MatrixProductState{T};
                           mode::Symbol=:B14,
                           maxdim::Int64=200) where {T <: RLorCX}
    @assert mps1.d == mps2.d == 2
    lx = mps1.lx
    @assert mps2.lx == lx

    if mps1.center != 1
        move_center!(mps1, 1)
    end
    if mps2.center != 1
        move_center!(mps2, 1)
    end

    dims = ones(Int64, lx+1)
    matrices = Array{T, 3}[]

    # Gutzwiller projector for mode :B14
    gutzp = zeros(T, 2,2,2)
    gutzp[1,1,1] = 1
    gutzp[2,2,2] = 1

    E = ones(T, 1,1,1)
    for site=1:lx-1
        mat1 = mps1.matrices[site]
        mat2 = mps2.matrices[site]

        diml1, dimr1 = size(mat1, 1), size(mat1, 3)
        diml2, dimr2 = size(mat2, 1), size(mat2, 3)
        dimle = size(E, 1)

        @tensor C[le, dg, r2, r1] := ((E[le, l2, l1] * mat1[l1,d1,r1]) *
                                      gutzp[dg,d1,d2]) * mat2[l2,d2, r2]

        fact = svd(reshape(C, dimle*2, dimr2*dimr1), full=false)
        S, n, ratio = MatrixProductStateTools.truncate(fact.S, maxdim=maxdim)
        dims[site+1] = n
        U = fact.U[:, 1:n]
        push!(matrices, reshape(U, dimle, 2, n))
        E = reshape(Diagonal(S) * fact.Vt[1:n, :], n, dimr2, dimr1)
    end
    mat1 = mps1.matrices[lx]
    mat2 = mps2.matrices[lx]
    @tensor C[le, dg, r2, r1] := ((E[le, l2, l1] * mat1[l1,d1,r1]) *
                                  gutzp[dg,d1,d2]) * mat2[l2,d2,r2]
    push!(matrices, reshape(C, size(E, 1), 2, 1))

    return MatrixProductState{T}(lx, 2, dims, matrices, lx)
end

function zipandgutzwiller!(mps1::SymMatrixProductState{Tv},
                           mps2::SymMatrixProductState{Tv};
                           mode::Symbol=:B14,
                           maxdim::Int64=200) where {Tv<:RLorCX}
    if mode==:B14
        return _zipandgutzwiller_B14!(mps1, mps2, maxdim=maxdim)
    elseif mode==:F23
        return _zipandgutzwiller_F23!(mps1, mps2, maxdim=maxdim)
    else
        error("mode not defined : ", mode)
    end
end

#bosonic version
function _zipandgutzwiller_B14!(mps1::SymMatrixProductState{Tv},
                                mps2::SymMatrixProductState{Tv};
                                maxdim::Int64=200) where {Tv<:RLorCX}
    @assert mps1.d == mps2.d == 2
    lx = mps1.lx
    @assert mps2.lx == lx

    if mps1.center != 1
        move_center!(mps1, 1)
    end
    if mps2.center != 1
        move_center!(mps2, 1)
    end

    dims = ones(Int64, lx+1)
    matrices = SymTensor{Tv, 3}[]

    ## NOTE: in order to make the gutzwiller projector respect the U1
    ## symmetry we need to do the follwoing. Assume that ↑↑
    ## corresponds to ↑ or 2 and ↓↓ corresponds or ↓ or 0.
    G = fill(one(Tv), 0,
             (STLeg(+1, [0,2], [1,1]),
              STLeg(-1, [0,1], [1,1]),
              STLeg(-1, [0,1], [1,1])))

    E = fill(one(Tv), 0, (STLeg(+1, [0], [1]),
                          STLeg(-1, [0], [1]),
                          STLeg(-1, [0], [1])))
    for l=1:lx-1
        A = mps1.matrices[l]
        B = mps2.matrices[l]

        ##NOTE:
        # contracting (E2, A1)
        # First tensor has EA = E1 E3 A2 A3
        # contracting (A2, G2)
        # Scond tensor has EAG = E1 E3 G1 A3 G3
        # contracting (E3, B1) and (G3, B2)
        # Final tensor has E1 G1 A3 B3
        C = contract(contract(contract(E, (1,-1, 2), A, (-1, 3, 4)),
                              (1, 2, -1, 4), G, (3, -1, 5)),
                     (1,-1, 2, 3,-2), B, (-1,-2, 4))


        U, S, Vt = svdtrunc(fuselegs(fuselegs(C, +1, 1, 2), -1, 2, 2), maxdim=maxdim)
        dims[l+1] = size(S, 1)
        push!(matrices,
              mapcharges(x->div(x,2),
                         unfuseleg(U, 1, (E.legs[1], G.legs[1]))))
        E = unfuseleg(S * Vt, 2, (A.legs[3], B.legs[3]))
    end
    A = mps1.matrices[lx]
    B = mps2.matrices[lx]
    C = contract(contract(contract(E, (1,-1, 2), A, (-1, 3, 4)),
                          (1, 2, -1, 4), G, (3, -1, 5)),
                 (1,-1, 2, 3,-2), B, (-1,-2, 4))

    push!(matrices, mapcharges(x->div(x,2), fuselegs(C, -1, 3, 2)))

    return SymMatrixProductState{Tv}(lx, 2, dims, matrices, lx)
end

#fermionic version (have to yet implement the minus sign and charge transform)
function _zipandgutzwiller_F23!(mps1::SymMatrixProductState{Tv},
                                mps2::SymMatrixProductState{Tv};
                                maxdim::Int64=200) where {Tv<:RLorCX}
    @assert mps1.d == mps2.d == 2
    lx = mps1.lx
    @assert mps2.lx == lx

    if mps1.center != 1
        move_center!(mps1, 1)
    end
    if mps2.center != 1
        move_center!(mps2, 1)
    end

    dims = ones(Int64, lx+1)
    matrices = SymTensor{Tv, 3}[]

    ## NOTE: in order to make the gutzwiller projector respect the U1
    ## symmetry we need to do the follwoing. Assume the first mps
    ## corresponds to ↑ or 1 and second mps to ↓ or -1.
    G = fill(one(Tv), 0,
             (STLeg(+1, [-1,1], [1,1]),
              STLeg(-1, [0,1], [1,1]),
              STLeg(+1, [0,1], [1,1])))

    E = fill(one(Tv), 0, (STLeg(+1, [0], [1]),
                          STLeg(-1, [0], [1]),
                          STLeg(+1, [0], [1])))
    for l=1:lx-1
        A = mps1.matrices[l]
        B = invlegs(mps2.matrices[l])

        fswap = fermionswapgate(A.legs[2], B.legs[1])
        ##NOTE:
        # contracting (E2, A1)
        # First tensor has EA = E1 E3 A2 A3
        # contracting with swap
        # second tensor has EAX = E1 A2 E3 A3
        # contracting (A2, G2)
        # Third tensor has EAG = E1 G1 G3 E3 A3
        # contracting (E3, B1) and (G3, B2)
        # Final tensor has E1 G1 A3 B3
        C = contract(contract(contract(contract(E, (1,-1, 2), A, (-1, 3, 4)),
                                       (1, -1, -2, 4), fswap, (-1, 2, 3, -2)),
                              (1, -1, 4, 5), G, (2, -1, 3)),
                     (1,2, -2, -1,3), B, (-1,-2, 4))

        U, S, Vt = svdtrunc(fuselegs(fuselegs(C, +1, 1, 2), -1, 2, 2), maxdim=maxdim)
        dims[l+1] = size(S, 1)

        fnl = x->div(x+l-1, 2)
        fnd = x->div(x+1, 2)
        fnr = x->div(x+l, 2)
        push!(matrices,
              mapcharges((fnl,fnd,fnr),
                         unfuseleg(U, 1, (E.legs[1], G.legs[1]))))
        E = unfuseleg(S * Vt, 2, (A.legs[3], B.legs[3]))
    end
    A = mps1.matrices[lx]
    B = invlegs(mps2.matrices[lx])
    fswap = fermionswapgate(A.legs[2], B.legs[1])
    C = contract(contract(contract(contract(E, (1,-1, 2), A, (-1, 3, 4)),
                                   (1, -1, -2, 4), fswap, (-1, 2, 3, -2)),
                          (1, -1, 4, 5), G, (2, -1, 3)),
                 (1,2, -2, -1,3), B, (-1,-2, 4))

    fnl = x->div(x+lx-1, 2)
    fnd = x->div(x+1, 2)
    fnr = x->div(x+lx, 2)
    push!(matrices, mapcharges((fnl,fnd,fnr), fuselegs(C, -1, 3, 2)))

    return SymMatrixProductState{Tv}(lx, 2, dims, matrices, lx)
end
