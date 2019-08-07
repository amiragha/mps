mutable struct SymMatrixProductState{Tv<:RLorCX}
    lx       :: Int64
    d        :: Int64
    dims     :: Vector{Int64}
    matrices :: Vector{SymTensor{Tv, 3}}
    center   :: Int64
end

# constructor with intial configuration vector
function SymMatrixProductState{Tv}(lx::Int64,
                                   initconf::Vector{Int64}) where{Tv<:RLorCX}
    @assert all((0 .<= initconf) .& (initconf .< 2))
    matrices = SymTensor{Tv,3}[]
    lchr = 0    # lchr is the "left charge"
    for site=1:lx
        if initconf[site] == 0
            push!(matrices, SymTensor((+1,+1,-1), (lchr,0,lchr), ones(Tv,1,1,1)))
        else
            newlchr = lchr + 1
            push!(matrices, SymTensor((+1,+1,-1), (lchr,1,newlchr), ones(Tv,1,1,1)))
            lchr = newlchr
        end
    end
    dims = ones(Int64, lx+1)
    #canonicalize_at!(matrices, lx)
    SymMatrixProductState{Tv}(lx, 2, dims, matrices, lx)
end

#constructor from a ketstate given in Ising basis
function SymMatrixProductState(lx::Int64, d::Int64,
                               ketstate::Vector{Tv};
                               svtruncation::Bool=false) where {Tv<:RLorCX}

    @assert length(ketstate) == binomial(lx, div(lx,2))
    charge = div(lx,2)

    matrices = SymTensor{Tv,3}[]

    dims = zeros(Int64, lx+1)
    dims[1] = 1

    statematrix = SymVector(+1, charge, ketstate)

    ℓ = 1
    ex_leg_ℓ = STLeg(1, [0], [1])
    physleg = STLeg(1, [0, 1], [1, 1])
    charge_range = min(charge, lx-ℓ):-1:max(0, charge-ℓ)
    smleg_ℓ = STLeg(-1, -1*charge_range, [binomial(lx-ℓ, c) for c in charge_range])
    legs = (STLeg(1, [0, 1], [1, 1]), smleg_ℓ)

    statematrix = defuse_leg(statematrix, 1, legs)
    U, S, Vt = svdsym(statematrix)
    dims[ℓ+1] = fulldims(U.legs[2])

    push!(matrices, defuse_leg(U, 1, (ex_leg_ℓ, physleg)))

    for ℓ = 2:lx-1
        #@show ℓ
        physleg = STLeg(1, [0, 1], [1, 1])
        ex_leg_ℓ = STLeg(+1, U.legs[2].chrs, U.legs[2].dims)

        charge_range = min(charge, lx-ℓ):-1:max(0, charge-ℓ)
        smleg_ℓ = STLeg(-1, -1*charge_range, [binomial(lx-ℓ, c) for c in charge_range])
        legs = (STLeg(1, [0, 1], [1, 1]), smleg_ℓ)
        statematrix = fuse_conseqlegs(defuse_leg(S*Vt, 2, legs), 1, 1, 2)
        U, S, Vt = svdsym(statematrix)
        dims[ℓ+1] = fulldims(U.legs[2])

        push!(matrices, defuse_leg(U, 1, (ex_leg_ℓ, physleg)))
    end

    dims[lx+1] = 1
    legs = (STLeg(1,[0, 1],[1, 1]), STLeg(-1, [0], [1]))
    push!(matrices, defuse_leg(S*Vt,2,legs))

    return SymMatrixProductState{Tv}(lx, d, dims, matrices, lx)
end

### conversions
###############

### TOOLS
#########

function isometrize_push_right!(matrices::Vector{SymTensor{Tv, 3}},
                                site::Int64;
                                svtruncation::Bool=false) where {Tv<:RLorCX}
    lx = length(matrices)
    if site < lx
        a = matrices[site]
        dims = size(a)
        U, S, Vt = svdsym(fuse_conseqlegs(a, +1, 1, 2))

        if svtruncation
            ## TODO
            # S, n, ratio =  truncate(fact.S, threshold=1.e-15)
            # U = fact.U[:,1:n]
            # Vt = fact.Vt[1:n,:]
        else
            S, n, ratio = S, size(S)[1], 1.
        end

        legs = a.legs[1:2]
        matrices[site] = defuse_leg(U, 1, legs)

        matrices[site+1] = contract(S*Vt, (1, -1), matrices[site+1], (-1, 2, 3))
    end
    return nothing
end

function isometrize_push_left!(matrices::Vector{SymTensor{Tv, 3}},
                               site::Int64;
                               svtruncation::Bool=false) where {Tv<:RLorCX}
    if site > 0
        a = matrices[site]
        dims = size(a)
        U, S, Vt = svdsym(fuse_conseqlegs(a, -1, 2, 2))

        if svtruncation
            ### TODO
            # S, n, ratio =  truncate(fact.S)
            # U = fact.U[:,1:n]
            # Vt = fact.Vt[1:n,:]
        else
            S, n, ratio = S, size(S)[1], 1.
        end

        legs = a.legs[2:3]
        matrices[site] = defuse_leg(Vt, 2, legs)

        matrices[site-1] = contract(matrices[site-1], (1, 2, -1), U*S, (-1, 3))
    end
    nothing
end

###
###NOTE: the following funcitons are exactly the same as the generic case
###
function correctdims!(mps::SymMatrixProductState{Tv}) where {Tv<:RLorCX}
    dims = Int64[]
    push!(dims, size(mps.matrices[1])[1])
    for n=1:mps.lx
        push!(dims, size(mps.matrices[n])[3])
    end
    mps.dims = dims
    nothing
end

function move_center!(mps::SymMatrixProductState{Tv},
                      new_center::Int64;
                      svnormalize::Bool=false) where {Tv<:RLorCX}
    lx = mps.lx
    @assert new_center > 0 && new_center <= lx

    center = mps.center
    if new_center > center
        for p=center:new_center-1
            isometrize_push_right!(mps.matrices, p)
        end
    elseif new_center < center
        for p=center:-1:new_center+1
            isometrize_push_left!(mps.matrices, p)
        end
    end
    correctdims!(mps)
    mps.center = new_center
    nothing
end
###
### NOTE: end of exactly the same functions
###

function measure_1point(mps::SymMatrixProductState{Tv}, op::SymTensor{Tv,2},
                        site::Int) where {Tv<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert (site <= lx) && (site > 0)
    @assert signs(op.legs) == (+1, -1)
    # @assert op.legs[1] == STLeg(+1,tuple(collect(1:d)), tuple(ones(Int, d)) )
    # @assert op.legs[2] == STLeg(-1,tuple(collect(1:d)), tuple(ones(Int, d)) )

    move_center!(mps, site)
    mat = mps.matrices[site]

    temp = contract(invlegs(conj(mat)), (1, -1, 3), op, (-1, 2))
    v = contract(mat, (-1, -2, -3), temp, (-1, -2, -3))
    v
end

function measure_1point(mps::SymMatrixProductState{Tv},
                        op::SymTensor{Tv,2}) where {Tv<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert signs(op.legs) == (+1, -1)
    # @assert op.legs[1] == STLeg(+1,tuple(collect(1:d)), tuple(ones(Int, d)) )
    # @assert op.legs[2] == STLeg(-1,tuple(collect(1:d)), tuple(ones(Int, d)) )

    result = Vector{Tv}(undef, lx)

    left = SymTensor(0,
                     (STLeg(+1, [0], [1]), STLeg(-1, [0], [1])),
                     [(0, 0)],
                     [ones(1,1)])
    move_center!(mps, 1)
    for site = 1:lx
        mat = mps.matrices[site]
        v = contract(contract(contract(left, (-1, 1), invlegs(conj(mat)), (-1, 2, 3)),
                              (1, -1, 3), op, (-1, 2)),
                     (-1,-2,-3), mat, (-1,-2,-3))

        result[site] =  v
        left = contract(contract(left, (-1, 1), invlegs(conj(mat)), (-1, 2, 3)),
                        (-1, -2, 1), mat, (-1,-2, 2))
    end
    result
end

function measure_2point(mps::SymMatrixProductState{Tv},
                        op1::SymTensor{Tv,2}, op2::SymTensor{Tv,2},
                        site1::Int64, site2::Int64) where {Tv<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert (site1 <= lx) && (site1 > 0)
    @assert (site2 <= lx) && (site2 > site1)
    @assert signs(op1.legs) == signs(op2.legs) == (+1, -1)

    move_center!(mps, site1)
    mat = mps.matrices[site1]

    left = contract(contract(invlegs(conj(mat)), (1, -1, 3), op1, (-1, 2)),
                    (-1, -2, 1), mat, (-1, -2, 2))
    for site=site1+1:site2-1
        mat = mps.matrices[site]
        left = contract(contract(left, (-1, 1), invlegs(conj(mat)), (-1, 2, 3)),
                        (-1, -2, 1), mat, (-1,-2, 2))
    end
    mat = mps.matrices[site2]
    v = contract(contract(contract(left, (-1, 1), invlegs(conj(mat)), (-1, 2, 3)),
                          (1, -1, 3), op2, (-1, 2)),
                 (-1,-2,-3), mat, (-1,-2,-3))
    v
end

function measure_2point(mps::SymMatrixProductState{Tv},
                        op1::SymTensor{Tv,2}, op2::SymTensor{Tv,2}) where {Tv<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert signs(op1.legs) == signs(op2.legs) == (+1, -1)

    result = Tv[]
    for site1=1:lx-1
        move_center!(mps, site1)
        mat = mps.matrices[site1]

        left = contract(contract(invlegs(conj(mat)), (1, -1, 3), op1, (-1, 2)),
                        (-1, -2, 1), mat, (-1, -2, 2))
        for site=site1+1:lx-1
            mat = mps.matrices[site]
            v = contract(contract(contract(left, (-1, 1), invlegs(conj(mat)), (-1, 2, 3)),
                              (1, -1, 3), op2, (-1, 2)),
                         (-1,-2,-3), mat, (-1,-2,-3))
            push!(result, v)
            left = contract(contract(left, (-1, 1), invlegs(conj(mat)), (-1, 2, 3)),
                            (-1, -2, 1), mat, (-1,-2, 2))
        end
        mat = mps.matrices[lx]
        v = contract(contract(contract(left, (-1, 1), invlegs(conj(mat)), (-1, 2, 3)),
                              (1, -1, 3), op2, (-1, 2)),
                     (-1,-2,-3), mat, (-1,-2,-3))
        push!(result,v)
    end
    result
end
