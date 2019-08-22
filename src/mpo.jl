"""
    MatrixProductOperator

`lx` is the length of the MPO or the number of sites.  `d` is the
physical on-site dimension, `dims` is the dimenions of bonds, starting
and ending with 1.  `tensors` is the local tensors, in four
dimenionsional tensor format, where the indeces are ordered starting
from the left one using the counterclockwise convention for the
indeces, so left, down, right, up.

"""
struct MatrixProductOperator{T<:RLorCX}
    lx :: Int64
    d :: Int64
    dims :: Vector{Int64}
    tensors :: Vector{Array{T,4}}
end

### constructors
################

function xxz_mpo(T::DataType, lx::Int64, d::Int64, delta::Float64=1.0) #where {T<:RLorCX}

    Sp = T[0.  1.; 0.  0.]
    Sm = T[0.  0.; 1.  0.]
    Sz = T[.5  0.; 0. -.5]
    I2 = Matrix{T}(I, 2, 2)

    mat = zeros(T, 5, d, 5, d)
    mat[1,:,1,:] = I2
    mat[2,:,1,:] = 0.5 * Sp
    mat[3,:,1,:] = 0.5 * Sm
    mat[4,:,1,:] = delta * Sz

    mat[5,:,2,:] = Sm
    mat[5,:,3,:] = Sp
    mat[5,:,4,:] = Sz
    mat[5,:,5,:] = I2

    tensors = Array{T,4}[]
    dims = zeros(Int64, lx+1)

    dims[1] = 1
    push!(tensors, mat[5:5,:,:,:])
    dims[2] = 5
    for site=2:lx-1
        push!(tensors, mat)
        dims[site+1] = 5
    end
    push!(tensors, mat[:,:,1:1,:])
    dims[lx+1] = 1

    MatrixProductOperator{T}(lx, d, dims, tensors)
end

function xxzlong_mpo(T::DataType, lx::Int64, d::Int64, delta::Float64=1.0, r::Float64=0.5)
    Sp = T[0.  1.; 0.  0.]
    Sm = T[0.  0.; 1.  0.]
    Sz = T[.5  0.; 0. -.5]
    I2 = Matrix{T}(I, 2, 2)

    mat = zeros(T, 5, d, 5, d)
    mat[1,:,1,:] = I2
    mat[2,:,1,:] = 0.5 * Sp
    mat[3,:,1,:] = 0.5 * Sm
    mat[4,:,1,:] = delta * Sz

    mat[5,:,2,:] = Sm
    mat[5,:,3,:] = Sp
    mat[5,:,4,:] = Sz
    mat[5,:,5,:] = I2

    mat[2,:,2,:] = r*I2
    mat[3,:,3,:] = r*I2
    mat[4,:,4,:] = r*I2

    tensors = Array{T,4}[]
    dims = zeros(Int64, lx+1)

    dims[1] = 1
    push!(tensors, mat[5:5,:,:,:])
    dims[2] = 5
    for site=2:lx-1
        push!(tensors, mat)
        dims[site+1] = 5
    end
    push!(tensors, mat[:,:,1:1,:])
    dims[lx+1] = 1

    MatrixProductOperator{T}(lx, d, dims, tensors)

end

function qitf_mpo(T::DataType, lx::Int64, d::Int64,
                  g::Float64=1.0, J::Float64=-1.0) #where {T<:RLorCX}

    Sp = T[0.  1.; 0.  0.]
    Sm = T[0.  0.; 1.  0.]
    Sz = T[1.  0.; 0. -1]
    Sx = Sp + Sm
    I2 = Matrix{T}(I, 2, 2)

    mat = zeros(T, 3, d, 3, d)
    mat[1,:,1,:] = I2
    mat[2,:,1,:] = J * Sz
    mat[3,:,1,:] = g/2. * Sx

    mat[3,:,2,:] = Sz
    mat[3,:,3,:] = I2

    tensors = Array{T,4}[]
    dims = zeros(Int64, lx+1)

    dims[1] = 1
    push!(tensors, mat[3:3,:,:,:])
    dims[2] = 3
    for site=2:lx-1
        push!(tensors, mat)
        dims[site+1] = 3
    end
    push!(tensors, mat[:,:,1:1,:])
    dims[lx+1] = 1

    MatrixProductOperator{T}(lx, d, dims, tensors)
end

function j1j2_mpo(Lx::Int64, j1::T, j2::T, d::Int64=2) where {T<:RLorCX}
    Id = Matrix{T}(I, d, d)
    mat = zeros(T, 8, d, 8, d)

    mat[1,:,1,:] = Id
    mat[2,:,1,:] = sz_half
    mat[3,:,1,:] = sp_half
    mat[4,:,1,:] = sm_half

    mat[5,:,2,:] = Id
    mat[6,:,3,:] = Id
    mat[7,:,4,:] = Id

    mat[8,:,2,:] = j1 * sz_half
    mat[8,:,3,:] = j1 * 0.5 * sm_half
    mat[8,:,4,:] = j1 * 0.5 * sp_half
    mat[8,:,5,:] = j2 * sz_half
    mat[8,:,6,:] = j2 * 0.5 * sm_half
    mat[8,:,7,:] = j2 * 0.5 * sp_half
    mat[8,:,8,:] = Id

    tensors = Array{T,4}[]
    dims = zeros(Int64, lx+1)

    dims[1] = 1
    push!(tensors, mat[8:8,:,:,:])
    dims[2] = 8
    for site=2:lx-1
        push!(tensors, mat)
        dims[site+1] = 8
    end
    push!(tensors, mat[:,:,1:1,:])
    dims[lx+1] = 1

    return MatrixProductOperator{T}(lx, d, dims, tensors)
end

### conversions
###############

function convert(::Type{MatrixProductOperator{ComplexF64}},
                 mpo::MatrixProductOperator{Float64})
    MatrixProductOperator{ComplexF64}(mpo.lx, mpo.d, mpo.dims,
                                      convert(Vector{Array{ComplexF64,4}}, mpo.tensors))
end
