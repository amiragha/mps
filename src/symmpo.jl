struct MPOperator{S, T}
    Ws :: Vector{SymTensor{S,T,4}}
end

const U1MPO{T} = MPOperator{Int, T}

"function to explicitly make the xxz mpo (for test and example)"
function xxz_u1mpo(T::DataType, lx::Int, d::Int, delta::Float64=1.0)

    V0 = U1Space(0=>1)
    V = U1Space(-1=>1, 0=>3, +1=>1)
    d = U1Space(0=>1, 1=>1)

    space = (V, d, dual(V), dual(d))

    Wbeg = fill(zero(T), 0, (V0, space[2:4]...))
    W    = fill(zero(T), 0, space)
    Wend = fill(zero(T), 0, (space[1:2]...,V0, space[4]))

    # making the generic MPO tensor
    mat0 = [
        1.   0.   0.;
        delta*-.5  0.   0.;
        0.  -0.5  1
    ]
    mat1 = [
        1.   0.   0.;
        delta*.5   0.   0.;
        0.   0.5  1
    ]
    W[Sector{U1}(0,0,0,0)] =  reshape(mat0, 3,1,3,1)
    W[Sector{U1}(0,1,0,1)] =  reshape(mat1, 3,1,3,1)

    W[Sector{U1}(0,1,1,0)] = reshape([0. 0. 1.], 3,1,1,1)
    W[Sector{U1}(0,0,-1,1)] = reshape([0. 0. 1.], 3,1,1,1)
    W[Sector{U1}(1,0,0,1)] = 0.5*reshape([1. 0. 0.], 1,1,3,1)
    W[Sector{U1}(-1,1,0,0)] = 0.5*reshape([1. 0. 0.], 1,1,3,1)

    # making the begging tensor
    Wbeg[Sector{U1}(0,0,0,0)] =  reshape(mat0[3,:], 1,1,3,1)
    Wbeg[Sector{U1}(0,1,0,1)] =  reshape(mat1[3,:], 1,1,3,1)

    Wbeg[Sector{U1}(0,1,1,0)] =  ones(T,1,1,1,1)
    Wbeg[Sector{U1}(0,0,-1,1)] =  ones(T,1,1,1,1)

    # making the end tensor
    Wend[Sector{U1}(0,0,0,0)] =  reshape(mat0[:,1], 3,1,1,1)
    Wend[Sector{U1}(0,1,0,1)] =  reshape(mat1[:,1], 3,1,1,1)

    Wend[Sector{U1}(-1,1,0,0)] =  0.5*ones(T,1,1,1,1)
    Wend[Sector{U1}(1,0,0,1)] =  0.5*ones(T,1,1,1,1)

    mpo = MPOperator{S,T}()
    push!(mpo, Wbeg)
    for site=2:lx-1
        push!(mpo, W)
    end
    push!(mpo, Wend)
    mpo
end

### conversions
###############

MPOperator{T}(mpo::MPOperator{T}) where {T<:Number} = mpo
function MPOperator{S,T}(mpo::MPOperator) where {S,T}
    MPOperator{S,T}([SymTensor{S, T, 4}(W) for W in mpo])
end

# function MatrixProductOperator(mpo::MPOperator{T}) where {T<:Number}
#     MatrixProductOperator{T}(mpo.lx, mpo.d, mpo.dims,
#                              [array(mat) for mat in mpo.tensors])
# end

function mpo2hamiltonian(mpo::MPOperator)
    lx = mpo.lx
    d = mpo.d
    mpo.d^lx > 10000 && error("model is too large for explicit Hamiltonian!")
    L = mpo.tensors[1]
    for two in mpo.tensors[2:lx-1]
        L = fuselegs(fuselegs(contract(L, (1, 2, -1, 5),
                                       two, (-1, 3, 4, 6)),
                              -1, 5, 2),
                     +1, 2, 2)
    end
    L = fuselegs(fuselegs(contract(L, (-2, 1, -1, 3),
                                       mpo.tensors[lx], (-1, 2, -2, 4)),
                              -1, 3, 2),
                 +1, 1, 2)
    L
end

function reducempo!(mpo::MPOperator)
    # this definitely needs to be better
    if size(mpo.tensors[1], 1) != 1
        return _reduceinfinitempo!(mpo)
    end
    A = mpo.tensors[1]
    for l=1:mpo.lx-1
        U, S, Vt = svdtrunc(fuselegs(permutelegs(A,
                                                 [1,2,4,3]),
                                     +1, 1, 3))
        mpo.dims[l+1] = size(S, 1)
        mpo.tensors[l] = permutelegs(unfuseleg(U*S, 1,
                                               A.legs[[1,2,4]]),
                                     [1,2,4,3])
        A = contract(Vt, (1, -1), mpo.tensors[l+1], (-1, 2, 3, 4))
    end

    for l=mpo.lx:-1:2
        U, S, Vt = svdtrunc(fuselegs(A, -1, 2, 3))
        mpo.dims[l] = size(S, 1)
        mpo.tensors[l] = unfuseleg(S*Vt, 2, A.legs[2:4])

        A = contract(mpo.tensors[l-1], (1,2,-1,4), U, (-1, 3))
    end
    mpo.tensors[1] = A
    mpo
end

function _reduceinfinitempo!(mpo::MPOperator)
    A = mpo.tensors[1]
    for l=1:mpo.lx-1
        U, S, Vt = svdtrunc(fuselegs(permutelegs(A,
                                                 [1,2,4,3]),
                                     +1, 1, 3))
        mpo.dims[l+1] = size(S, 1)
        mpo.tensors[l] = permutelegs(unfuseleg(U*S, 1,
                                               A.legs[[1,2,4]]),
                                     [1,2,4,3])
        A = contract(Vt, (1, -1), mpo.tensors[l+1], (-1, 2, 3, 4))
    end
    U, S, Vt = svdtrunc(fuselegs(permutelegs(A,
                                             [1,2,4,3]),
                                 +1, 1, 3))
    mpo.dims[mpo.lx+1] = size(S, 1)
    mpo.dims[1] = size(S, 1)
    mpo.tensors[mpo.lx] = permutelegs(unfuseleg(U*S, 1,
                                            A.legs[[1,2,4]]),
                                  [1,2,4,3])
    A = contract(Vt, (1, -1), mpo.tensors[1], (-1, 2, 3, 4))

    U, S, Vt = svdtrunc(fuselegs(A, -1, 2, 3))
    mpo.dims[mpo.lx+1] = size(S, 1)
    mpo.dims[1] = size(S, 1)
    mpo.tensors[1] = unfuseleg(S*Vt, 2, A.legs[2:4])

    A = contract(mpo.tensors[mpo.lx], (1,2,-1,4), U, (-1, 3))

    for l=mpo.lx:-1:2
        U, S, Vt = svdtrunc(fuselegs(A, -1, 2, 3))
        mpo.dims[l] = size(S, 1)
        mpo.tensors[l] = unfuseleg(S*Vt, 2, A.legs[2:4])

        A = contract(mpo.tensors[l-1], (1,2,-1,4), U, (-1, 3))
    end
    mpo.tensors[1] = A#mpo.tensors[mpo.lx] = A

    mpo
end
