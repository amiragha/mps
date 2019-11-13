struct SymMatrixProductOperator{Tv<:RLorCX}
    lx :: Int
    d  :: Int
    dims :: Vector{Int}
    tensors :: Vector{SymTensor{Tv, 4}}
end

function xxz_symmpo(Tv::DataType, lx::Int, d::Int, delta::Float64=1.0)
    legs = (
        STLeg(+1, [-1,0,+1], [1,3,1]),
        STLeg(+1, [0,+1], [1,1]),
        STLeg(-1, [-1,0,+1], [1,3,1]),
        STLeg(-1, [0,+1], [1,1]),
    )

    legbeg = STLeg(+1, [0], [1])
    legend = STLeg(-1, [0], [1])

    Abeg = fill(zero(Tv), 0, (legbeg, legs[2:4]...))
    A = fill(zero(Tv), 0, legs)
    Aend = fill(zero(Tv), 0, (legs[1:2]...,legend, legs[4]))

    ###NOTE: here we need nice functions for symtensor which we don't have!!!
    ###TODO: Make the nice functions (slicing and indexing, etc) for SymTensors!

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
    set_sector!(A, (0,0,0,0), reshape(mat0, 3,1,3,1))
    set_sector!(A, (0,1,0,1), reshape(mat1, 3,1,3,1))

    set_sector!(A, (0,1,1,0), reshape([0. 0. 1.], 3,1,1,1))
    set_sector!(A, (0,0,-1,1), reshape([0. 0. 1.], 3,1,1,1))
    set_sector!(A, (1,0,0,1), 0.5*reshape([1. 0. 0.], 1,1,3,1))
    set_sector!(A, (-1,1,0,0), 0.5*reshape([1. 0. 0.], 1,1,3,1))

    # making the begging tensor
    set_sector!(Abeg, (0,0,0,0), reshape(mat0[3,:], 1,1,3,1))
    set_sector!(Abeg, (0,1,0,1), reshape(mat1[3,:], 1,1,3,1))

    set_sector!(Abeg, (0,1,1,0), ones(Tv,1,1,1,1))
    set_sector!(Abeg, (0,0,-1,1), ones(Tv,1,1,1,1))

    # making the end tensor
    set_sector!(Aend, (0,0,0,0), reshape(mat0[:,1], 3,1,1,1))
    set_sector!(Aend, (0,1,0,1), reshape(mat1[:,1], 3,1,1,1))

    set_sector!(Aend, (-1,1,0,0), 0.5*ones(Tv,1,1,1,1))
    set_sector!(Aend, (1,0,0,1), 0.5*ones(Tv,1,1,1,1))

    tensors = SymTensor{Tv, 4}[]
    dims = zeros(Int, lx+1)

    dims[1] = 1
    push!(tensors, Abeg)
    dims[2] = 5
    for site=2:lx-1
        push!(tensors, A)
        dims[site+1] = 5
    end
    push!(tensors, Aend)
    dims[lx+1] = 1

    SymMatrixProductOperator{Tv}(lx, d, dims, tensors)
end

### conversions
###############

function MatrixProductOperator(mpo::SymMatrixProductOperator{T}) where {T<:Number}
    MatrixProductOperator{T}(mpo.lx, mpo.d, mpo.dims,
                          [array(mat) for mat in mpo.tensors])
end

function mpo2hamiltonian(mpo::SymMatrixProductOperator)
    lx = mpo.lx
    d = mpo.d
    mpo.d^lx > 10000 && error("model is too large for explicit Hamiltonian!")
    L = mpo.tensors[1]
    for two in mpo.tensors[2:end]
        L = fuselegs(fuselegs(contract(L, (1, 2, -1, 5),
                                       two, (-1, 3, 4, 6)),
                              -1, 5, 2),
                     +1, 2, 2)
    end
    removedummyleg(removedummyleg(L, 3), 1)
end
