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

    Abeg = fillSymTensor(zero(Tv), 0, (legbeg, legs[2:4]...))
    A = fillSymTensor(zero(Tv), 0, legs)
    Aend = fillSymTensor(zero(Tv), 0, (legs[1:2]...,legend, legs[4]))

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
    change_nzblk!(A, (0,0,0,0), reshape(mat0, 3,1,3,1))
    change_nzblk!(A, (0,1,0,1), reshape(mat1, 3,1,3,1))

    change_nzblk!(A, (0,1,1,0), reshape([0. 0. 1.], 3,1,1,1))
    change_nzblk!(A, (0,0,-1,1), reshape([0. 0. 1.], 3,1,1,1))
    change_nzblk!(A, (1,0,0,1), 0.5*reshape([1. 0. 0.], 1,1,3,1))
    change_nzblk!(A, (-1,1,0,0), 0.5*reshape([1. 0. 0.], 1,1,3,1))

    # making the begging tensor
    change_nzblk!(Abeg, (0,0,0,0), reshape(mat0[3,:], 1,1,3,1))
    change_nzblk!(Abeg, (0,1,0,1), reshape(mat1[3,:], 1,1,3,1))

    change_nzblk!(Abeg, (0,1,1,0), ones(Tv,1,1,1,1))
    change_nzblk!(Abeg, (0,0,-1,1), ones(Tv,1,1,1,1))

    # making the end tensor
    change_nzblk!(Aend, (0,0,0,0), reshape(mat0[:,1], 3,1,1,1))
    change_nzblk!(Aend, (0,1,0,1), reshape(mat1[:,1], 3,1,1,1))

    change_nzblk!(Aend, (-1,1,0,0), 0.5*ones(Tv,1,1,1,1))
    change_nzblk!(Aend, (1,0,0,1), 0.5*ones(Tv,1,1,1,1))

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