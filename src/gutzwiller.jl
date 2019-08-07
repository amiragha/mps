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

    if mps1.center == 1
        move_center!(mps1, 1)
    end
    if mps2.center == 1
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
