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
    lx :: Int
    d :: Int
    dims :: Vector{Int64}
    tensors :: Vector{Array{T,4}}
end

struct InfMatrixProductOperator{T<:Number, OP<:Union{Array{T, 4}, SymTensor{T, 4}}}
    ly :: Int
    ds :: Vector{Int}
    Ws :: Vector{OP}
end
### constructors
################

function xxz_mpo(T::DataType, lx::Int64, d::Int64, delta::Float64=1.0) #where {T<:RLorCX}

    Sz, Sp, Sm = spinoperators(1/2)
    I2 = I(2)

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
    Sz, Sp, Sm = spinoperators(1/2)
    I2 = I(2)

    mat = zeros(T, 5, d, 5, d)
    mat[1,:,1,:] = I2
    mat[2,:,1,:] = 0.5 * Sp
    mat[3,:,1,:] = 0.5 * Sm
    mat[4,:,1,:] = delta * Sz

    mat[5,:,2,:] = Sm
    mat[5,:,3,:] = Sp
    mat[5,:,4,:] = Sz
    mat[5,:,5,:] = I2

    mat[2,:,2,:] = r.*I2
    mat[3,:,3,:] = r.*I2
    mat[4,:,4,:] = r.*I2

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

function j1j2_mpo(lx::Int64, j1::T, j2::T, d::Int64=2) where {T<:RLorCX}
    mat = zeros(T, 8, d, 8, d)
    sz, sp, sm = spinoperators(1/2)
    Id = I(d)
    mat[1,:,1,:] = Id
    mat[2,:,1,:] = sz
    mat[3,:,1,:] = sp
    mat[4,:,1,:] = sm

    mat[5,:,2,:] = Id
    mat[6,:,3,:] = Id
    mat[7,:,4,:] = Id

    mat[8,:,2,:] = j1 * sz
    mat[8,:,3,:] = j1 * 0.5 * sm
    mat[8,:,4,:] = j1 * 0.5 * sp
    mat[8,:,5,:] = j2 * sz
    mat[8,:,6,:] = j2 * 0.5 * sm
    mat[8,:,7,:] = j2 * 0.5 * sp
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

function mpo2hamiltonian(mpo::MatrixProductOperator)
    lx = mpo.lx
    d = mpo.d
    mpo.d^lx > 10000 && error("model is too large for explicit Hamiltonian!")
    L = mpo.tensors[1]
    for two in mpo.tensors[2:end]
        @tensor L[l,o1,o2,r,o1',o2'] := L[l,o1,m,o1'] * two[m, o2, r, o2']
        dd = size(L,2)*d
        L = reshape(L, 1, dd, size(two,3), dd)
    end
    reshape(L, d^lx, d^lx)
end

function reducempo!(mpo::MatrixProductOperator)
    A = mpo.tensors[1]
    for l=1:mpo.lx-1
        U, S, Vt = svdtrunc(reshape(permutedims(A, [1,2,4,3]),
                                     prod(size(A)[[1,2,4]]), size(A, 3)))
        mpo.dims[l+1] = size(S, 1)
        mpo.tensors[l] = permutedims(reshape(U*S,
                                               size(A)[[1,2,4]]..., size(S, 1)),
                                     [1,2,4,3])
        @tensor A[l,o,r,o'] := Vt[l, m] * mpo.tensors[l+1][m,o,r,o']
    end

    for l=mpo.lx:-1:2
        U, S, Vt = svdtrunc(reshape(A, size(A, 1), prod(size(A)[2:4])))
        mpo.dims[l] = size(S, 1)
        mpo.tensors[l] = reshape(S*Vt, size(S, 1), size(A)[2:4]...)

        @tensor A[l,o,r,o'] := mpo.tensors[l-1][l,o,m,o'] *  U[m, r]
    end
    mpo.tensors[1] = A
    mpo
end
