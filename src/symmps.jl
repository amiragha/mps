abstract type AbstractMPState end
mutable struct MPState{S<:AbstractCharge, T<:Number} <: AbstractMPState
    As     :: Vector{SymTensor{S,T,3}}
    center :: Int
end

@inline Base.eltype(mps::MPState) = eltype(mps.As[1])
@inline Base.length(mps::MPState) = length(mps.As)
@inline checkbounds(mps, l) = 0 < l <= L || throw(BoundsError(mps.vectors), l)

@inline bondspace(mps::MPState, l::Int) = space(mps.As[l], 3)
@inline bonddim(mps::MPState, l::Int) = dim(mps.As[l], 3)
@inline bonddim(mps::MPState) =
    [dim(mps.As[1], 1); [bonddim(mps, l) for l in 1:length(mps)]]

@inline sitespace(mps::MPState, l::Int) = space(mps.As[l], 2)

@inline function Base.push!(mps::MPState{S},
                            A::SymTensor{S}) where{S}
    lx = length(mps)
    Vl = space(A, 1)
    if lx > 0
        bondspace(mps, lx) == dual(Vl) || throw("SpaceMismatch()")
    else
        dim(Vl) == 1 || throw("SpaceMismatch()")
    end
    push!(mps.As, A)
    mps.center = 0
    mps
end

const U1MPState{T} = MPState{U1, T}
####################
### CONSTRUCTORS ###
####################
MPState{S,T}() where{S,T} = MPState{S,T}(Vector{SymTensor{S,T,3}}(), 0)

# equal probability constructor
function U1MPState(::Type{T},
                   lx::Int,
                   d::Int,
                   m::U1Charge;
                   noise::T=0.0) where {T}

    mps = MPState{U1, T}()

    # Just find all the possible sectors for each tensor at each site
    # and fill them with one(1,1,1). This gives an MPS with a norm of
    # binomial(lx,m).
    ##TODO: there sure is a better way to the following!
    Vd = U1Space(c=>1 for c in 0:d-1)
    Vl = U1Space(0=>1)
    for site in 1:lx
        cmin = max(0, m-(d-1)*(lx-site))
        cmax = min(m, (d-1)*site)
        rchrs = collect(cmin:cmax)
        Vr = U1Space(c=>1 for c in rchrs)

        A = SymTensor(ones + noise * randn, zero(U1), (Vl, Vd, dual(Vr)))
        push!(mps, A)
        Vl = Vr
    end
    center_at!(mps, 1)
    mps
end

# constructor with intial configuration vector
function U1MPState(::Type{T},
                   lx::Int,
                   d::Int,
                   initconf::Vector{Int}) where {T}

    @assert all((0 .<= initconf) .& (initconf .< 2))
    mps = MPState{U1, T}()
    Vd = U1Space(c=>1 for c in 0:d-1)
    lchr = 0
    for site=1:lx
        Vl = U1Space(lchr=>1)
        if initconf[site] == 0
            Vr = U1Space(lchr=>1)
            A = fill(one(T), zero(U1), (Vl, Vd, dual(Vr)))
            push!(mps, A)
        else
            Vr = U1Space(lchr+1=>1)
            A = fill(one(T), zero(U1), (Vl, Vd, dual(Vr)))
            push!(mps, A)
            lchr += 1
        end
    end
    center_at!(mps, 1)
    mps
end

#constructor from a ketstate given in Ising basis
function U1MPState(lx  ::Int,
                   d   ::Int,
                   ketstate::Vector{T};
                   svtruncation::Bool=false) where {T}

    @assert length(ketstate) == binomial(lx, div(lx,2))
    charge = div(lx,2)

    mps = MPState{U1, T}()

    A = SymVector(U1(charge), ketstate)
    #println(A)

    ℓ = 1
    ex_leg_ℓ = U1Space(0=>1)
    Vd = U1Space(0=>1, 1=>1)
    charge_range = min(charge, lx-ℓ):-1:max(0, charge-ℓ)
    smleg_ℓ = U1Space(c=>binomial(lx-ℓ, c) for c in charge_range)

    A = splitleg(A, 1, (Vd, smleg_ℓ))
    #println(A)
    u,s,v = svd(A)

    # println(u)
    # println(s)
    # println(v)
    #println(u*s*v ≈ A)
    push!(mps, splitleg(u, 1, (ex_leg_ℓ, Vd)))

    for ℓ = 2:lx-1
        #@show ℓ

        ex_leg_ℓ = space(s, 1)
        charge_range = min(charge, lx-ℓ):-1:max(0, charge-ℓ)
        smleg_ℓ = U1Space(c=>binomial(lx-ℓ, c) for c in charge_range)
        #println(s*v)
        A = fuselegs(splitleg(s*v, 2, (Vd, smleg_ℓ)), 1, 2)
        #println(A)
        u,s,v = svd(A)
        #println(u*s*v ≈ A)

        push!(mps, splitleg(u, 1, (ex_leg_ℓ, Vd)))
    end

    push!(mps, splitleg(s*v, 2, (Vd, U1Space(0=>1))))

    mps
end

### conversions
###############

MPState{T}(mps::MPState{T}) where {T<:Number} = mps
function MPState{T}(mps::MPState) where {T<:Number}
    MPState{T}(mps.lx, mps.d, mps.dims,
               [SymTensor{T, 3}(mat) for mat in mps.matrices],
               mps.center)
end

function MatrixProductState(mps::MPState{T}) where {T<:Number}
    MatrixProductState{T}(mps.lx, mps.d, mps.dims,
                          [array(mat) for mat in mps.matrices], mps.center)
end

### TOOLS
#########

function normalize!(mps::MPState{T}) where {T}
    A = mps.matrices[mps.center]
    if mps.center > div(mps.lx, 2)
        U,S,Vt = svd(fuselegs(A, +1, 1, 2))
        ss = vcat(diag.(S.nzblks)...)
        n = norm(ss)

        S_nzblks = S.nzblks ./ n
        S = SymDiagonal(S.charge, S.legs, S.sects, S_nzblks)

        mps.matrices[mps.center] = unfuseleg(U * S * Vt, 1, A.legs[1:2])
    else
        U,S,Vt = svd(fuselegs(A, -1, 2, 2))
        ss = vcat(diag.(S.nzblks)...)
        n = norm(ss)

        S_nzblks = S.nzblks ./ n
        S = SymDiagonal(S.charge, S.legs, S.sects, S_nzblks)

        mps.matrices[mps.center] = unfuseleg(U * S * Vt, 2, A.legs[2:3])
    end
    sort(ss./n, rev=true)
end

@inline function _pushleft!(mps::MPState, l::Int)
    @boundscheck 1 < l <= length(mps) ||
        throw("can't pushleft at $l of $(length(mps))")
    u,s,v = svd(fuselegs(mps.As[l], 2, 2))
    mps.As[l] = splitleg(v, 2, space(mps.As[l])[2:3])
    mps.As[l-1] = contract(mps.As[l-1], (1, 2, -1), u*s, (-1, 3))
    nothing
end

@inline function _pushright!(mps::MPState, l::Int)
    @boundscheck 0 < l < length(mps) ||
        throw("can't pushright at $l of $(length(mps))")
    u,s,v = svd(fuselegs(mps.As[l], 1, 2))
    mps.As[l] = splitleg(u, 1, space(mps.As[l])[1:2])
    mps.As[l+1] = contract(s*v, (1,-1), mps.As[l+1], (-1,2,3))
    nothing
end

function center_at!(mps::MPState, center::Int)
    x1 = 1
    xL = length(mps)
    if mps.center > 0
        x1 = mps.center
        xL = mps.center
    end
    for p=x1:center-1
        _pushright!(mps, p)
    end
    for p=xL:-1:center+1
        _pushleft!(mps, p)
    end
    mps.center = center
    mps
end

function schmidtvalues(mps::MPState)
    lx = mps.lx
    values = Vector{Vector{Float64}}(undef, lx-1)
    center_at!(mps, 1)
    A = mps.matrices[1]
    for l = 1:lx-1
        u,s,v = svd(fuselegs(A, 1, 2))
        values[l] = sort(vcat([diag(blk) for (c,blk) in s.blocks]...), rev=true)
        A = contract(s*v, (1, -1), mps.As[l+1], (-1,2,3))
    end
    values
end

entanglementspectrum(mps::MPState) = -1 .* log.(schmidtvalues(mps))
entanglemententropy(mps::MPState) = [entropy(x.^2) for x in schmidtvalues]

"""
    measure(mps, ops...[, xs...])

Measures the expectation values of `ops` for the given `mps`. If the
sites `xs` are not indicated then the functions measures
everywhere. The `ops` can be single site operators (matrices) or MPOs as in
```
measure(mps, mpo)
```

It returns a dictionary of position of measurements as keys and
results as values.

"""
function measure end

function measure(mps::MPState{S,T}, op::SymMatrix{S,T}, l::Int) where {S,T}
    lx = length(mps)
    @assert (l <= lx) && (l > 0)
    @assert sitespace(mps, l) == space(op, 2) == dual(space(op, 1))

    center_at!(mps, site)
    A = mps.As[site]
    v = contract(A, (-1, -2, -3),
                 contract(dual(A), (1, -1, 3), op, (-1, 2)), (-1, -2, -3))
    Dict{Int, T}(l=>v)
end

function measure(mps::MPState{S,T}, op::SymMatrix{S,T}) where {S,T}

    Vd = space(op, 1)
    @assert space(op, 2) == Vd

    lx = length(mps)
    values = Dict{Int, T}
    center_at!(mps, 1)
    L = SymTensor(ones, zero(S), (U1Space(0=>1), U1Space(0=>1)))

    for x = 1:lx
        A = mps.As[x]
        v = contract(contract(contract(L, (-1, 1), dual(A), (-1, 2, 3)),
                              (1, -1, 3), op, (-1, 2)),
                     (-1,-2,-3), A, (-1,-2,-3))
        values[x] = v
        L = contract(contract(L, (-1, 1), dual(A), (-1, 2, 3)),
                     (-1, -2, 1), A, (-1,-2, 2))
    end
    values
end

function measure(mps::MPState{S,T},
                 op1::SymMatrix{S,T},
                 op2::SymMatrix{S,T},
                 x1::Int,
                 x2::Int) where {S,T}
    lx = length(mps)
    @assert (0 < x1 < x2) && (x2 <= lx)

    Vd1 = space(op1, 1)
    @assert space(op1, 2) == Vd1
    Vd2 = space(op2, 1)
    @assert space(op2, 2) == Vd2
    (op1.charge + op2.charge != zero(S)) &&
        error("operator charges don't add up to zero! $(op1.charge), $(op2.charge)")

    values = Dict{NTuple{2, Int}, T}
    center_at!(mps, site1)
    A = mps.As[x1]

    L = contract(contract(dual(A), (1, -1, 3), op1, (-1, 2)),
                 (-1, -2, 1), A, (-1, -2, 2))
    for x = x1+1:x2-1
        A = mps.As[x]
        L = contract(contract(L, (-1, 1), dual(A), (-1, 2, 3)),
                     (-1, -2, 1), A, (-1,-2, 2))
    end
    A = mps.As[x2]
    v = contract(contract(contract(L, (-1, 1), dual(A), (-1, 2, 3)),
                          (1, -1, 3), op2, (-1, 2)),
                 (-1,-2,-3), A, (-1,-2,-3))
    values[(x1, x2)] = v
end

function measure(mps::MPState{S,T},
                 op1::SymMatrix{S,T},
                 op2::SymMatrix{S,T}) where {S,T}

    lx = length(mps)
    Vd1 = space(op1, 1)
    @assert space(op1, 2) == dual(Vd1)
    Vd2 = space(op2, 1)
    @assert space(op2, 2) == dual(Vd2)
    (op1.charge + op2.charge != zero(S)) &&
        error("operator charges don't add up to zero! $(op1.charge), $(op2.charge)")

    values = Dict{NTuple{2,Int}, T}()
    for x1=1:lx-1
        center_at!(mps, x1)
        A = mps.As[x1]
        L = contract(contract(dual(A), (1, -1, 3), op1, (-1, 2)),
                     (-1, -2, 1), A, (-1, -2, 2))
        for x=x1+1:lx-1
            A = mps.As[x]
            v = contract(contract(contract(L, (-1, 1), dual(A), (-1, 2, 3)),
                                  (1, -1, 3), op2, (-1, 2)),
                         (-1,-2,-3), A, (-1,-2,-3))
            values[(x1, x)] = v
            L = contract(contract(L, (-1, 1), dual(A), (-1, 2, 3)),
                         (-1, -2, 1), A, (-1,-2, 2))
        end
        A = mps.As[lx]
        v = contract(contract(contract(L, (-1, 1), dual(A), (-1, 2, 3)),
                              (1, -1, 3), op2, (-1, 2)),
                     (-1,-2,-3), A, (-1,-2,-3))
        values[x1, lx] = v
    end
    values
end

function measure(mps::MPState{S,T1},
                 op1::SymMatrix{S,T2},
                 op2::SymMatrix{S,T3}) where {S,T1,T2,T3}
    T = promote_type(T1,T2,T3)
    measure(convert(MPState{S,T}, mps),
            convert(SymMatrix{S,T}, op1),
            convert(SymMatrix{S,T}, op2))
end

function measure(mps::MPState{S,T},
                 mpo::MPOperator{S,T}) where{S,T}
    @assert mps.d == mpo.d && mps.lx == mpo.lx
    values = Dict{Int, T}()
    dummy = VectorSpace{S}(0=>1)
    L = fill(one(T), zero(S), (dummy, dummy, dummy))

    for x=1:mps.lx
        A = mps.As[x]
        W = mpo.Ws[x]
        L = contract(contract(contract(
            L, (-1,3,4), A, (-1,2,1)),
                              (1, -1, -2, 4), W, (-2, 3, 2,-1)),
                     (1,2,-1,-2), dual(A), (-2,-1,3))
    end
    values[1] = contract(L, (-1,-2,-3), fill(one(T), (dummy,dummy,dummy)), (-1,-2,-3))
    values
end

"""
    apply!(mps, op, x[, maxdim, pushto, svnormalize])

applies the `operator` which is `d x d x d x d` sym tensor to site `x`
and `x+1` of the `mps`. The order of indeces start from the bottom
left, bottom right, top left, top right. So the counterclockwise
convention is not assumed here!

The `maxdim` operator chooses the max possible size of dimension of
the new mps at bond between `x` and `x+1`, The singular values are
push to either left `:L` or right `:R` (default) matrices using the
argument `pushto`.

"""
function apply!(mps         :: MPState{S,T},
                op          :: SymTensor{S,T,4},
                x           :: Int;
                maxdim      :: Int=bonddim(mps, l+1),
                pushto      :: Symbol=:R,
                svnormalize :: Bool=false) where {S,T}

    @boundscheck 0 < x < length(mps) || throw()

    if mps.center < x
        center_at!(mps, x)
    elseif mps.center > x+1
        center_at!(mps, x+1)
    end

    A1 = mps.As[l]
    A2 = mps.As[l+1]

    RR = contract(A1, (1,2,-1), A2, (-1,3,4))
    R = contract(RR,
                 (1,-1,-2,4), op, (2,3,-1,-2))

    #display(fuselegs(fuselegs(R, -1, 3, 2), +1, 1, 2))
    u,s,v = svdtrunc(fuselegs(fuselegs(R, -1, 3, 2), +1, 1, 2),
                     maxdim=maxdim, tol=1.e-14)

    svnormalize && normalize!(s)

    space = space(R)
    if (pushto == :R)
        mps.As[x] = splitleg(u, 1, space[1:2])
        mps.As[x+1] = splitleg(s*v, 2, space[3:4])
        mps.center = x+1
    elseif (pushto == :L)
        mps.As[x] = splitleg(u*s, 1, space[1:2])
        mps.As[x+1] = splitleg(v, 2, space[3:4])
        mps.center = x
    else
        error("invalid push_to :", pushto)
    end
    mps
end

function apply!(mps         :: MPState{S,T1},
                op          :: SymTensor{S,T2,4},
                x           :: Int;
                maxdim      :: Int=bonddim(mps, l+1),
                pushto      :: Symbol=:R,
                svnormalize :: Bool=false) where {S,T1,T2}

    T = promote_type(T1, T2)
    apply!(convert(MPState{S,T}, mps),
           convert(SymTensor{S,T,4}, op),
           x, maxdim=maxdim, pushto=pushto, svnormalize=svnormalize)
end

"""
    overlap(mps1, mps2)

calculates the overlap between two matrix product states `mps1` and
`mps2` that is to run the tensor contraction corresponding to
``⟨ψ_2|ψ_1⟩``. Note that it is not divided by the norm of the two
MPSs, so it returned value is the overlap of the two states multiplied
by the norm of each.

"""
function overlap(mps1::MPState{S,T},
                 mps2::MPState{S,T}) where {S,T}

    lx = length(mps1)
    @assert length(mps2) == lx

    A = mps1.As[1]
    B = mps2.As[1]
    L = contract(A, (-1, -2, 1), dual(B), (-1, -2, 2))

    for x=2:lx-1
        A = mps1.As[x]
        B = mps2.As[x]
        L = contract(contract(L, (-1, 1), A, (-1, 2, 3)),
                     (-1,-2,1), dual(B), (-1,-2,2))
    end
    A = mps1.As[lx]
    B = mps2.As[lx]
    v = contract(contract(L, (-1, 1), A, (-1, 2, 3)),
                 (-1,-2,-3), dual(B), (-1,-2,-3))
    v
end

function overlap(mps1::MPState{S,T1}, mps2::MPState{S,T2}) where {S,T1,T2}
    T = promote_type(T1, T2)
    overlap(convert(MPSteat{S,T}, mps1), convert(MPSteat{S,T}, mps2))
end

"""
    norm(mps)

calculates the norm of a matrix product state `mps` that is to
calculate the sqrt tensor contraction corresponding to ``⟨ψ|ψ⟩``.

"""
@inline norm(mps::MPState) = sqrt(real(overlap(mps, mps)))
