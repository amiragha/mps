# DetMatrix type to represent determinantal wavefunctions
mutable struct DetMatrix{T}
    matrix      :: Matrix{T}
    determinant :: T
    inverse     :: Matrix{T}
end

@inline Base.eltype(::DetMatrix{T}) where {T} = T

# constructor with just a matrix
function DetMatrix(mat::Matrix)
    fact = lu(mat)
    #return DetMatrix(mat, det(mat), inv(mat))
    return DetMatrix(mat, det(fact), inv(fact))
end

"""
    det_ratio_1col(dmat, col, c)

Making use of the Matrix Determinant Lemma (MDL), finds the ratio of
the new determinant to old one given by the change of `c`th
column to `col` in the `dmat`.

returns determinant ratio.
"""
function det_ratio_1col(dmat :: DetMatrix,
                        col  :: Vector,
                        c    :: Int)

    ### TODO: use better dot product probably `vecdot`
    return (dmat.inverse[[c],:] * col)[1]
end

"""
    det_ratio_2col(dmat, cols, c1, c2)

Making use of the Matrix Determinant Lemma (MDL), finds the matrix 2x2
that its determinant is ratio of the new determinant to old one given
by the change of `c1, c2`th columns to the two columns in `cols`
respectively.

returns determinant ratio matrix (2x2)
"""
function det_ratio_2cols(dmat      :: DetMatrix,
                         cols      :: Matrix,
                         c1        :: Int,
                         c2        :: Int)

    ### TODO: find a vectorize way to do this?
    T = eltype(dmat)
    detratio_mat = zeros(T, 2, 2)
    detratio_mat[1, 1] = (dmat.inverse[[c1],:] * cols[:,1])[1]
    detratio_mat[1, 2] = (dmat.inverse[[c1],:] * cols[:,2])[1]
    detratio_mat[2, 1] = (dmat.inverse[[c2],:] * cols[:,1])[1]
    detratio_mat[2, 2] = (dmat.inverse[[c2],:] * cols[:,2])[1]
    return detratio_mat
end

"""
    update_detmatrix_1col!(dmat, new_col, col_index, detratio)

Making use of the Sherman-Morrison formula, updates the matrix,
determinant and inverse of the new DetMatrix given by the change of
`col_index` column to `new_col`.

"""
function update_detmatrix_1col!(dmat      :: DetMatrix,
                                col   :: Vector,
                                c :: Int,
                                detratio)
    # update matrix
    dmat.matrix[:,c] = col

    # update determinant
    dmat.determinant *= detratio

    # update inverse, note: column is an auxiliary vector
    # TODO: elegance!
    column = dmat.inverse * col
    column[c] -= 1

    dmat.inverse -= (column * dmat.inverse[[c],:])/detratio
    return nothing
end

"""
    update_detmatrix_2cols!(dmat, cols, c1, c1, detratio)

Making use of the generalized Sherman-Morrison formula (Woodbury?!),
updates the matrix, determinant and inverse of the new DetMatrix given
by the change of `c1,c2` columns to the two columns in `cols` respectively.

updates dmat.
"""
function update_detmatrix_2cols!(dmat      :: DetMatrix,
                                 cols      :: Matrix,
                                 c1        :: Int,
                                 c2        :: Int,
                                 detratio_mat  ::  Matrix)
    # update matrix
    dmat.matrix[:,c1] = cols[:, 1]
    dmat.matrix[:,c2] = cols[:, 2]

    # update determinant
    dmat.determinant *= det(detratio_mat)

    # update inverse, note: columns is an auxiliary 2 column matrix
    # TODO: elegance!
    columns = dmat.inverse * cols
    columns[c1, 1] -= 1
    columns[c2, 2] -= 1

    dmat.inverse -= columns * inv(detratio_mat) * dmat.inverse[[c1,c2], :]
    return nothing
end

"""
    check_and_update_detmatrix!(dmat)

check the inverse and determinant in the detmat and if they are off,
replace them with exact values. TODO: What should the tolerance be here?!
"""
function check_and_update_detmatrix!(dmat::DetMatrix)
    fact = lu(dmat.matrix)
    exact_det = det(fact)
    exact_inv = inv(fact)
    if !isapprox(dmat.determinant, exact_det, rtol=1.e-14)
        print("d")
    end
    # if !isapprox(dmat.matrix * dmat.inverse,
    #             eye(ComplexF64, size(dmat.matrix)[1]), rtol=1.e-13)
    if !isapprox(dmat.inverse, exact_inv, rtol=1.e-13)
        print("i")
    end
    dmat.determinant = exact_det
    dmat.inverse = exact_inv
    return nothing
end
