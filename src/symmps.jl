mutable struct SymMatrixProductState{Tv<:RLorCX}
    lx       :: Int
    d        :: Int
    dims     :: Vector{Int}
    matrices :: Vector{SymTensor{Tv, 3}}
    center   :: Int64
end

# equal probability constructor
function SymMatrixProductState{Tv}(
    lx::Int,
    d::Int,
    m::Int;
    noise::Float64=0.0) where{Tv <:RLorCX}

    matrices = SymTensor{Tv, 3}[]
    dims = zeros(Int, lx+1)

    # Just find all the possible sectors for each tensor at each site
    # and fill them with one(1,1,1). This gives an MPS with a norm of
    # binomial(lx,m). How to make it unit?
    ##TODO: there sure is a better way to the following!
    lchrs = [0]
    rvals = [1]
    dims[1] = 1
    for site in 1:lx
        legl = STLeg(+1, lchrs, ones(Int, length(lchrs)))
        legd = STLeg(+1, collect(0:d-1), ones(Int, d))
        cmin = max(0, m-(d-1)*(lx-site))
        cmax = min(m, (d-1)*site)
        rchrs = collect(cmin:cmax)
        dims[site+1] = length(rchrs)
        legr = STLeg(-1, rchrs, ones(Int, length(rchrs)))
        sects = NTuple{3, Int}[]
        for c in lchrs
            for n in 0:d-1
                if cmax >= c+n >= cmin
                    push!(sects, (c, n, c+n))
                end
            end
        end
        sects = sort(sects, lt=SymTensors._sectorlessthan)
        l = length(sects)
        A = SymTensor(0,
                      (legl,legd,legr),
                      sects,
                      #[1/sqrt(length(rchrs)) .*
                      [(ones(Tv,1,1,1) + noise * randn(Tv,1,1,1)) for n in 1:l])
        push!(matrices, A)
        lchrs = rchrs
    end

    canonicalize_at!(matrices, lx)
    mps = SymMatrixProductState{Tv}(lx,d,dims,matrices,lx)
    mps
end

# constructor with intial configuration vector
function SymMatrixProductState{Tv}(
    lx::Int,
    d::Int,
    initconf::Vector{Int}) where{Tv<:RLorCX}

    @assert all((0 .<= initconf) .& (initconf .< 2))
    matrices = SymTensor{Tv,3}[]
    lchr = 0
    for site=1:lx
        legl = STLeg(+1, [lchr], [1])
        legd = STLeg(+1, [0, 1], [1,1])
        if initconf[site] == 0
            legr = STLeg(-1, [lchr], [1])
            A = fill(one(Tv), 0, (legl,legd,legr))
            push!(matrices, A)
        else
            legr = STLeg(-1, [lchr+1], [1])
            A = fill(one(Tv), 0, (legl,legd,legr))
            push!(matrices, A)
            lchr += 1
        end
    end
    dims = ones(Int64, lx+1)
    #canonicalize_at!(matrices, lx)
    SymMatrixProductState{Tv}(lx, 2, dims, matrices, lx)
end

#constructor from a ketstate given in Ising basis
function SymMatrixProductState(
    lx::Int64, d::Int64,
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

    statematrix = unfuseleg(statematrix, 1, legs)
    U, S, Vt = svdsym(statematrix)
    dims[ℓ+1] = fulldims(U.legs[2])

    push!(matrices, unfuseleg(U, 1, (ex_leg_ℓ, physleg)))

    for ℓ = 2:lx-1
        #@show ℓ
        physleg = STLeg(1, [0, 1], [1, 1])
        ex_leg_ℓ = STLeg(+1, U.legs[2].chrs, U.legs[2].dims)

        charge_range = min(charge, lx-ℓ):-1:max(0, charge-ℓ)
        smleg_ℓ = STLeg(-1, -1*charge_range, [binomial(lx-ℓ, c) for c in charge_range])
        legs = (STLeg(1, [0, 1], [1, 1]), smleg_ℓ)
        statematrix = fuselegs(unfuseleg(S*Vt, 2, legs), 1, 1, 2)
        U, S, Vt = svdsym(statematrix)
        dims[ℓ+1] = fulldims(U.legs[2])

        push!(matrices, unfuseleg(U, 1, (ex_leg_ℓ, physleg)))
    end

    dims[lx+1] = 1
    legs = (STLeg(1,[0, 1],[1, 1]), STLeg(-1, [0], [1]))
    push!(matrices, unfuseleg(S*Vt,2,legs))

    return SymMatrixProductState{Tv}(lx, d, dims, matrices, lx)
end

### conversions
###############

function MatrixProductState(mps::SymMatrixProductState{T}) where {T<:Number}
    MatrixProductState{T}(mps.lx, mps.d, mps.dims,
                          [array(mat) for mat in mps.matrices], mps.center)
end

### TOOLS
#########

function normalize!(mps::SymMatrixProductState{Tv}) where {Tv<:RLorCX}
    A = mps.matrices[mps.center]
    if mps.center > div(mps.lx, 2)
        U,S,Vt = svdsym(fuselegs(A, +1, 1, 2))
        ss = vcat(diag.(S.nzblks)...)
        n = norm(ss)

        S_nzblks = S.nzblks ./ n
        S = SymDiagonal(S.charge, S.legs, S.sects, S_nzblks)

        mps.matrices[mps.center] = unfuseleg(U * S * Vt, 1, A.legs[1:2])
    else
        U,S,Vt = svdsym(fuselegs(A, -1, 2, 2))
        ss = vcat(diag.(S.nzblks)...)
        n = norm(ss)

        S_nzblks = S.nzblks ./ n
        S = SymDiagonal(S.charge, S.legs, S.sects, S_nzblks)

        mps.matrices[mps.center] = unfuseleg(U * S * Vt, 2, A.legs[2:3])
    end
    sort(ss./n, rev=true)
end

function canonicalize_at!(matrices::Vector{SymTensor{Tv, 3}},
                          center::Int) where {Tv<:RLorCX}
    lx = length(matrices)
    @assert center > 0 && center < lx + 1

    for site=1:center-1
        isometrize_push_right!(matrices, site, svtruncation=false)
    end

    for site=lx:-1:center+1
        isometrize_push_left!(matrices, site, svtruncation=false)
    end
end

function isometrize_push_right!(matrices::Vector{SymTensor{Tv, 3}},
                                site::Int64;
                                svtruncation::Bool=false) where {Tv<:RLorCX}
    lx = length(matrices)
    if site < lx
        A = matrices[site]
        dims = size(A)
        U, S, Vt = svdsym(fuselegs(A, +1, 1, 2))

        if svtruncation
            ## TODO
            # S, n, ratio =  truncate(fact.S, threshold=1.e-15)
            # U = fact.U[:,1:n]
            # Vt = fact.Vt[1:n,:]
        else
            S, n, ratio = S, size(S)[1], 1.
        end

        legs = A.legs[1:2]
        matrices[site] = unfuseleg(U, 1, legs)

        matrices[site+1] = contract(S*Vt, (1, -1), matrices[site+1], (-1, 2, 3))
    end
    return nothing
end

function isometrize_push_left!(matrices::Vector{SymTensor{Tv, 3}},
                               site::Int64;
                               svtruncation::Bool=false) where {Tv<:RLorCX}
    if site > 0
        A = matrices[site]
        dims = size(A)
        U, S, Vt = svdsym(fuselegs(A, -1, 2, 2))

        if svtruncation
            ### TODO
            # S, n, ratio =  truncate(fact.S)
            # U = fact.U[:,1:n]
            # Vt = fact.Vt[1:n,:]
        else
            S, n, ratio = S, size(S)[1], 1.
        end

        legs = A.legs[2:3]
        matrices[site] = unfuseleg(Vt, 2, legs)

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

function entanglementspectrum(mps::SymMatrixProductState)
    lx = mps.lx
    result = Vector{Vector{Float64}}(undef, lx-1)
    move_center!(mps, 1)
    A = mps.matrices[1]

    for l = 1:lx-1
        U, S, Vt = svdsym(fuselegs(A, +1, 1, 2))
        mps.matrices[l] = unfuseleg(U, 1, A.legs[1:2])
        spectrum = sort(vcat([diag(blk) for blk in S.nzblks]...), rev=true)
        result[l] = spectrum
        A = contract(S*Vt, (1, -1), mps.matrices[l+1], (-1,2,3))
    end

    mps.matrices[lx] = A
    mps.center = lx
    result
end

function entanglemententropy(mps::SymMatrixProductState;
                             alpha::Int=1)
    [entropy(spectrum.^2) for spectrum in entanglementspectrum(mps)]
end

function measure_1point(mps::SymMatrixProductState{Tv}, op::SymMatrix{Tv},
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
                        op::SymMatrix{Tv}) where {Tv<:RLorCX}
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
                        op1::SymMatrix{Tv}, op2::SymMatrix{Tv},
                        site1::Int64, site2::Int64) where {Tv<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert (site1 <= lx) && (site1 > 0)
    @assert (site2 <= lx) && (site2 > site1)
    @assert signs(op1.legs) == signs(op2.legs) == (+1, -1)
    (op1.charge + op2.charge != 0) &&
        error("operator charges don't add up to zero! $(op1.charge), $(op2.charge)")

    # why move mps center?
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
                        op1::SymMatrix{Tv}, op2::SymMatrix{Tv}) where {Tv<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert signs(op1.legs) == signs(op2.legs) == (+1, -1)
    (op1.charge + op2.charge != 0) &&
        error("operator charges don't add up to zero! $(op1.charge), $(op2.charge)")

    result = Tv[]
    for site1=1:lx-1
        # why do I need to move the center to where my operator measures?
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

function measure_2point(mps::SymMatrixProductState{ComplexF64},
                        op1::SymTensor{Float64,2},
                        op2::SymTensor{Float64,2}) where {Tv<:RLorCX}
    measure_2point(mps,
                   convert(SymTensor{ComplexF64, 2}, op1),
                   convert(SymTensor{ComplexF64, 2}, op2))
end

function measure_mpo(mps::SymMatrixProductState{Tv},
                     mpo::SymMatrixProductOperator{Tv}) where{Tv<:RLorCX}
    @assert mps.d == mpo.d && mps.lx == mpo.lx

    L = fill(one(Tv), 0, (STLeg(-1, [0], [1]),
                          STLeg(-1, [0], [1]),
                          STLeg(+1, [0], [1])))

    for site=1:mps.lx
        A = mps.matrices[site]
        W = mpo.tensors[site]
        L = contract(contract(contract(
            L, (-1,3,4), A, (-1,2,1)),
                              (1, -1, -2, 4), W, (-2, 3, 2,-1)),
                     (1,2,-1,-2), invlegs(conj(A)), (-2,-1,3))
    end
    L.nzblks[1][1,1,1]
end

"""
    apply_2siteoperator!(mps, l, operator, maxdim, pushto)

applies the `operator` which is `d x d x d x d` sym tensor to site `l`
and `l+1` of the `mps`. The order of indeces start from the bottom
left, bottom right, top left, top right. So the counterclockwise
convention is not assumed here!

The `maxdim` operator chooses the max possible size of dimension of
the new mps at bond between `l` and `l+1`, The singular values are
push to either left `:L` or right `:R` (default) matrices using the
argument `pushto`.

"""
function apply_2siteoperator!(mps        :: SymMatrixProductState{Tv},
                              l          :: Int,
                              op         :: SymTensor{Tv, 4};
                              maxdim     :: Int=mps.dims[l+1],
                              pushto     :: Symbol=:R,
                              normalizeS :: Bool=false) where {Tv<:RLorCX}

    #@assert mps_dims_are_consistent(mps)
    @assert 0 < l < mps.lx

    d = mps.d

    if mps.center < l
        move_center!(mps, l)
    elseif mps.center > l+1
        move_center!(mps, l+1)
    end

    one = mps.matrices[l]
    two = mps.matrices[l+1]

    dim_l = mps.dims[l]
    dim_m = mps.dims[l+1]
    dim_r = mps.dims[l+2]

    #display(one)
    #display(two)
    RR = contract(one, (1,2,-1), two, (-1,3,4))
    #display(RR)
    #display(op)
    R = contract(RR,
                 (1,-1,-2,4), op, (2,3,-1,-2))

    #display(fuselegs(fuselegs(R, -1, 3, 2), +1, 1, 2))
    U, S, Vt = svdtrunc(fuselegs(fuselegs(R, -1, 3, 2), +1, 1, 2),
                        maxdim=maxdim, tol=1.e-14)

    normalizeS && normalize!(S)
    mps.dims[l+1] = size(U, 2)

    if (pushto == :R)
        mps.matrices[l] = unfuseleg(U, 1, (R.legs[1], R.legs[2]))
        mps.matrices[l+1] = unfuseleg(S * Vt, 2, (R.legs[3], R.legs[4]))
        mps.center = l+1
    elseif (pushto == :L)
        mps.matrices[l] = unfuseleg(U * S, 1, (R.legs[1], R.legs[2]))
        mps.matrices[l+1] = unfuseleg(Vt, 2, (R.legs[3], R.legs[4]))
        mps.center = l
    else
        error("invalid push_to :", pushto)
    end
    nothing
end

function apply_2siteoperator!(mps        :: SymMatrixProductState{ComplexF64},
                              l          :: Int64,
                              op         :: SymTensor{Float64, 4};
                              maxdim     :: Int64=mps.dims[l+1],
                              pushto     :: Symbol=:R,
                              normalizeS :: Bool=false)

    apply_2siteoperator!(mps, l, convert(SymTensor{ComplexF64, 4}, op),
                         maxdim=maxdim, pushto=pushto, normalizeS=normalizeS)
end

"""
    overlap(mps1, mps2)

calculates the overlap between two matrix product states `mps1` and
`mps2` that is to run the tensor contraction corresponding to
``⟨ψ_2|ψ_1⟩``. Note that it is not divided by the norm of the two
MPSs, so it returned value is the overlap of the two states multiplied
by the norm of each.

"""
function overlap(mps1::SymMatrixProductState{Tv},
                 mps2::SymMatrixProductState{Tv}) where {Tv<:RLorCX}

    lx = mps1.lx
    @assert mps1.d == mps2.d && mps2.lx == lx

    A = mps1.matrices[1]
    B = mps2.matrices[1]
    L = contract(A, (-1, -2, 1), conj(invlegs(B)), (-1, -2, 2))

    for l=2:lx-1
        A = mps1.matrices[l]
        B = mps2.matrices[l]
        L = contract(contract(L, (-1, 1), A, (-1, 2, 3)),
                     (-1,-2,1), conj(invlegs(B)), (-1,-2,2))
    end
    A = mps1.matrices[lx]
    B = mps2.matrices[lx]
    v = contract(contract(L, (-1, 1), A, (-1, 2, 3)),
                 (-1,-2,-3), conj(invlegs(B)), (-1,-2,-3))
    v
end

overlap(mps1::SymMatrixProductState{ComplexF64}, mps2::SymMatrixProductState{Float64}) =
    overlap(mps1, convert(SymMatrixProductState{ComplexF64}, mps2))
overlap(mps1::SymMatrixProductState{Float64}, mps2::SymMatrixProductState{ComplexF64}) =
    overlap(convert(SymMatrixProductState{ComplexF64}, mps1), mps2)

"""
    norm2(mps)

calculates the norm of a matrix product state `mps` that is to
calculate the tensor contraction corresponding to ``⟨ψ|ψ⟩``.

"""
norm2(mps::SymMatrixProductState{Tv}) where {Tv<:RLorCX} =
    real(overlap(mps, mps))


"""
    mps_dims_are_consistent(mps)

check if are dimensions are consistent in an mps. This is made for
testing and double checks; in principle all operations must not break
the consistency of the bond dimensions of MPS.

"""
function mps_dims_are_consistent(mps::SymMatrixProductState{Tv}) where {Tv<:RLorCX}
    for n=1:mps.lx-1
        dim1 = size(mps.matrices[n], 3)
        dim2 = size(mps.matrices[n+1], 1)
        if dim1 != dim2 || dim1 != mps.dims[n+1]
            return false
        end
    end
    return size(mps.matrices[mps.lx], 3) == mps.dims[mps.lx+1]
end
