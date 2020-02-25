struct MPOperator{Y<:Tensor{T,4} where T}
    Ws :: Vector{Y}
end

@inline tensortype(::MPOperator{Y}) where {Y} = Y
@inline Base.eltype(::MPOperator{Y}) where {Y} = eltype(Y)
@inline Base.length(mpo::MPOperator) = length(mpo.Ws)
#@inline checkbounds(mpo, l) = 0 < l <= L || throw(BoundsError(mpo.vectors), l)

@inline bondspace(mpo::MPOperator, l::Int) = dual(space(mpo.Ws[l], 3))
@inline bonddim(mpo::MPOperator, l::Int) = size(mpo.Ws[l], 3)
@inline bonddim(mpo::MPOperator) =
    [size(mpo.Ws[1], 1); [bonddim(mpo, l) for l in 1:length(mpo)]]

@inline leftspace(mpo::MPOperator) = space(mpo.Ws[1], 1)
@inline rightspace(mpo::MPOperator) = bondspace(mpo, length(mpo))

@inline sitespace(mpo::MPOperator, l::Int) = space(mpo.Ws[l], 2)

@inline function Base.push!(mpo::MPOperator{Y},
                            W::Y) where{Y}
    lx = length(mpo)

    isdual(space(W, 2), space(W, 4))
    if lx > 0
        isequal(bondspace(mpo, lx), space(W,1)) || throw("SpaceMismatch()")
    else
        size(W,1) == 1 || throw("SpaceMismatch()")
    end
    push!(mpo.Ws, W)
    mpo
end

const U1MPO{T} = MPOperator{SymTensor{U1,T,4}}
const MPO{T} = MPOperator{Array{T,4}}

MPOperator{Y}() where{Y} = MPOperator{Y}(Vector{Y}())

function convert(::Type{MPOperator{Y1}},
                 mpo::MPOperator{Y2}) where {Y1, Y2}
    MPOperator{Y1}(convert(Vector{Y1}, mpo.Ws))
end

"function to explicitly make the xxz mpo (for test and example)"
function xxz_u1mpo(T::DataType, lx::Int, d::Int, delta::Float64=1.0)

    V0 = U1Space(0=>1)
    V = U1Space(-1=>1, 0=>3, +1=>1)
    d = U1Space(0=>1, 1=>1)

    space = (V, d, dual(V), dual(d))

    Wbeg = fill(zero(T), U1(0), (V0, space[2:4]...))
    W    = fill(zero(T), U1(0), space)
    Wend = fill(zero(T), U1(0), (space[1:2]...,V0, space[4]))

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

    mpo = U1MPO{T}()
    push!(mpo, Wbeg)
    for site=2:lx-1
        push!(mpo, W)
    end
    push!(mpo, Wend)
    mpo
end

### conversions
###############

MPOperator{Y}(mpo::MPOperator{Y}) where {Y} = mpo
function MPOperator{Y}(mpo::MPOperator) where {Y}
    MPOperator{Y}([Y(W) for W in mpo])
end

# function MatrixProductOperator(mpo::MPOperator{T}) where {T<:Number}
#     MatrixProductOperator{T}(mpo.lx, mpo.d, mpo.dims,
#                              [array(mat) for mat in mpo.tensors])
# end

function mpo2hamiltonian(mpo::MPOperator)
    lx = length(mpo)
    ds = [dim(sitespace(mpo, l)) for l in 1:lx]
    prod(ds) > 10000 && error("model is too large for explicit Hamiltonian!")
    H = mpo.Ws[1]
    for two in mpo.Ws[2:lx-1]
        H = fuselegs(fuselegs(contract(H, (1, 2, -1, 5),
                                       two, (-1, 3, 4, 6)),
                              5, 2),
                     2, 2)
    end
    H = fuselegs(fuselegs(contract(H, (-2, 1, -1, 3),
                                       mpo.Ws[lx], (-1, 2, -2, 4)),
                              3, 2),
                 1, 2)
    H
end

function reducempo!(mpo::MPOperator)
    if dim(rightspace(mpo)) != 1
        println(rightspace(mpo))
        return _reduceinfinitempo!(mpo)
    end
    lx = length(mpo)
    W = mpo.Ws[1]
    for l=1:lx-1
        u,s,v = svdtrunc(fuselegs(permutelegs(W,
                                              [1,2,4,3]),
                                  1, 3))
        mpo.Ws[l] = permutelegs(splitleg(u*s, 1,
                                         W.space[[1,2,4]]),
                                [1,2,4,3])
        W = contract(v, (1, -1), mpo.Ws[l+1], (-1, 2, 3, 4))
    end

    for l=lx:-1:2
        u,s,v = svdtrunc(fuselegs(W, 2, 3, true))
        mpo.Ws[l] = splitleg(s*v, 2, W.space[2:4])

        W = contract(mpo.Ws[l-1], (1,2,-1,4), u, (-1, 3))
    end
    mpo.Ws[1] = W
    mpo
end

function _reduceinfinitempo!(mpo::MPOperator)
    A = mpo.Ws[1]
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
