mutable struct MPState{Y<:Tensor{T,3} where T}
    As     :: Vector{Y}
    center :: Int
end

@inline center(mps::MPState) = mps.center
@inline tensortype(::MPState{Y}) where {Y} = Y
@inline Base.eltype(::MPState{Y}) where {Y} = eltype(Y)
@inline Base.length(mps::MPState) = length(mps.As)
#@inline checkbounds(mps, l) = 0 < l <= L || throw(BoundsError(mps.vectors), l)

@inline bondspace(mps::MPState, l::Int) = dual(space(mps.As[l], 3))
@inline bonddim(mps::MPState, l::Int) = size(mps.As[l], 3)
@inline bonddim(mps::MPState) =
    [dim(mps.As[1], 1); [bonddim(mps, l) for l in 1:length(mps)]]

@inline leftspace(mps::MPState) = space(mps.As[1], 1)
@inline rightspace(mps::MPState) = bondspace(mps, length(mps))

@inline sitespace(mps::MPState, l::Int) = space(mps.As[l], 2)

@inline function Base.push!(mps::MPState{Y}, A::Y) where{Y}
    lx = length(mps)
    Vl = space(A, 1)
    if lx > 0
        bondspace(mps, lx) == Vl || throw("SpaceMismatch()")
    else
        dim(Vl) == 1 || throw("SpaceMismatch()")
    end
    push!(mps.As, A)
    mps.center = 0
    mps
end

const U1MPS{T} = MPState{SymTensor{U1,T,3}}
const MPS{T} = MPState{Array{T, 3}}

####################
### CONSTRUCTORS ###
####################
MPState{Y}() where{Y} = MPState{Y}(Vector{Y}(), 0)

# equal probability constructor
function U1MPS(::Type{T},
               lx::Int,
               d::Int,
               m::U1Charge;
               noise::T=0.0) where {T}

    mps = U1MPS{T}()
    M = m.charge

    # Just find all the possible sectors for each site and put one(T)
    # This gives an MPS with a norm of binomial(lx,m) if no noise.
    Vd = U1Space(c=>1 for c in 0:d-1)
    Vl = U1Space(0=>1)
    for site in 1:lx
        cmin = max(0, M-(d-1)*(lx-site))
        cmax = min(M, (d-1)*site)
        Vr = U1Space(c=>1 for c in cmin:cmax)
        A = SymTensor((x,y)->ones(x,y) .+ noise .* randn(x,y),
                      zero(U1), (Vl, Vd, dual(Vr)))
        push!(mps, A)
        Vl = Vr
    end
    center_at!(mps, 1)
    mps
end

# constructor with intial configuration vector
function U1MPS(::Type{T},
               lx::Int,
               d::Int,
               initconf::Vector{Int}) where {T}

    @assert all((0 .<= initconf) .& (initconf .< 2))
    mps = U1MPS{T}()
    Vd = U1Space(c=>1 for c in 0:d-1)
    lchr = 0
    for site=1:lx
        Vl = U1Space(lchr=>1)
        if initconf[site] == 0
            Vr = U1Space(lchr=>1)
        else
            Vr = U1Space(lchr+1=>1)
            lchr += 1
        end
        push!(mps, fill(one(T), (Vl, Vd, dual(Vr))))
    end
    center_at!(mps, 1)
    mps
end

#constructor from a ketstate given in Ising basis
function U1MPS(lx  ::Int,
               d   ::Int,
               ketstate::Vector{T};
               svtruncation::Bool=false) where {T}

    @assert length(ketstate) == binomial(lx, div(lx,2))
    charge = div(lx,2)
    A = SymVector(U1(charge), ketstate)

    mps = U1MPS{T}()
    ℓ = 1
    ex_leg_ℓ = U1Space(0=>1)
    Vd = U1Space(0=>1, 1=>1)
    charge_range = min(charge, lx-ℓ):-1:max(0, charge-ℓ)
    smleg_ℓ = U1Space(c=>binomial(lx-ℓ, c) for c in charge_range)

    A = splitleg(A, 1, (Vd, smleg_ℓ))
    u,s,v = _svd_(A)

    push!(mps, splitleg(u, 1, (ex_leg_ℓ, Vd)))

    for ℓ = 2:lx-1
        ex_leg_ℓ = space(s, 1)
        charge_range = min(charge, lx-ℓ):-1:max(0, charge-ℓ)
        smleg_ℓ = U1Space(c=>binomial(lx-ℓ, c) for c in charge_range)
        A = fuselegs(splitleg(s*v, 2, (Vd, smleg_ℓ)), 1, 2)
        u,s,v = _svd_(A)

        push!(mps, splitleg(u, 1, (ex_leg_ℓ, Vd)))
    end
    push!(mps, splitleg(s*v, 2, (Vd, U1Space(0=>1))))
    mps
end

### conversions
###############

MPState{Y}(mps::MPState{Y}) where {Y} = mps
function MPState{Y}(mps::MPState) where {Y}
    MPState{Y}(mps.As, mps.center)
end

function MPState{Y1}(mps::MPState{Y2}) where {Y1,Y2}
    MPState{Y1}([convert(Y1, mps.As[i]) for i in eachindex(mps.As)], mps.center)
end

### TOOLS
#########

function normalize!(mps::MPState)
    if mps.center < 1
        center_at!(mps, 1)
    end
    A = mps.As[mps.center]
    if mps.center > div(length(mps), 2)
        u,s,v = _svd_(fuselegs(A, 1, 2, false))
        normalize!(s)
        mps.As[mps.center] = splitleg(u*s*v, 1, A.space[1:2])
    else
        u,s,v = _svd_(fuselegs(A, 2, 2, true))
        normalize!(s)
        mps.As[mps.center] = splitleg(u*s*v, 2, A.space[2:3])
    end
    mps
end

@inline function _pushleft!(mps::MPState, l::Int)
    @boundscheck 1 < l <= length(mps) ||
        throw("can't pushleft at $l of $(length(mps))")
    u,s,v = _svd_(fuselegs(mps.As[l], 2, 2))
    mps.As[l] = splitleg(v, 2, space(mps.As[l])[2:3])
    mps.As[l-1] = contract(mps.As[l-1], (1, 2, -1), u*s, (-1, 3))
    nothing
end

@inline function _pushright!(mps::MPState, l::Int)
    @boundscheck 0 < l < length(mps) ||
        throw("can't pushright at $l of $(length(mps))")
    u,s,v = _svd_(fuselegs(mps.As[l], 1, 2))
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
    lx = length(mps)
    values = Vector{Vector{Float64}}(undef, lx-1)
    center_at!(mps, 1)
    A = mps.As[1]
    for l = 1:lx-1
        u,s,v = _svd_(fuselegs(A, 1, 2))
        if typeof(s) <: Diagonal
            values[l] = diag(s)
        elseif typeof(s) <: SymDiagonal
            values[l] = sort(vcat([diag(blk) for (c,blk) in s.blocks]...), rev=true)
        else
            throw("unknown type for singular values!")
        end
        A = contract(s*v, (1, -1), mps.As[l+1], (-1,2,3))
    end
    values
end

entanglementspectrum(mps::MPState) = -1 .* log.(schmidtvalues(mps))
entanglemententropy(mps::MPState) = [entropy(x.^2) for x in schmidtvalues(mps)]

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

function measure(mps::MPState{SymTensor{S,T,3}},
                 op::SymMatrix{S,T}, l::Int) where {S,T}
    lx = length(mps)
    @assert (l <= lx) && (l > 0)
    @assert sitespace(mps, l) == space(op, 2) == dual(space(op, 1))

    center_at!(mps, site)
    A = mps.As[site]
    v = contract(A, (-1, -2, -3),
                 contract(dual(A), (1, -1, 3), op, (-1, 2)), (-1, -2, -3))
    Dict{Int, T}(l=>v)
end

function measure(mps::MPState{SymTensor{S,T,3}},
                 op::SymMatrix{S,T}) where {S,T}

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

function measure(mps::MPState{SymTensor{S,T,3}},
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

function measure(mps::MPState{SymTensor{S,T,3}},
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

function measure(mps::MPState{SymTensor{S,T1,3}},
                 op1::SymMatrix{S,T2},
                 op2::SymMatrix{S,T3}) where {S,T1,T2,T3}
    T = promote_type(T1,T2,T3)
    measure(convert(MPState{SymTensor{S,T,3}}, mps),
            convert(SymMatrix{S,T}, op1),
            convert(SymMatrix{S,T}, op2))
end

function measure(mps::MPState{SymTensor{S,T,3}},
                 mpo::MPOperator{SymTensor{S,T,4}}) where{S,T}
    @assert length(mpo) == length(mps)
    values = Dict{Int, T}()
    dummy = VectorSpace{S}(0=>1)
    L = fill(one(T), zero(S), (dummy, dummy, dummy))

    for x=1:length(mps)
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
function apply!(mps         :: MPState{SymTensor{S,T,3}},
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

    A1 = mps.As[x]
    A2 = mps.As[x+1]

    RR = contract(A1, (1,2,-1), A2, (-1,3,4))
    R = contract(RR,
                 (1,-1,-2,4), op, (2,3,-1,-2))

    u,s,v = svdtrunc(SymMatrix(R, [1,2], [3,4]),
                     maxdim=maxdim, tol=1.e-14)

    svnormalize && normalize!(s)

    _space = space(R)
    if (pushto == :R)
        mps.As[x] = splitleg(u, 1, _space[1:2])
        mps.As[x+1] = splitleg(s*v, 2, _space[3:4])
        mps.center = x+1
    elseif (pushto == :L)
        mps.As[x] = splitleg(u*s, 1, _space[1:2])
        mps.As[x+1] = splitleg(v, 2, _space[3:4])
        mps.center = x
    else
        error("invalid push_to :", pushto)
    end
    mps
end

function apply!(mps         :: MPState{SymTensor{S,T1,3}},
                op          :: SymTensor{S,T2,4},
                x           :: Int;
                maxdim      :: Int=bonddim(mps, l+1),
                pushto      :: Symbol=:R,
                svnormalize :: Bool=false) where {S,T1,T2}

    T = promote_type(T1, T2)
    apply!(convert(MPState{SymTensor{S,T,3}}, mps),
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
function overlap(mps1::MPState{Y},
                 mps2::MPState{Y}) where {Y}

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

function overlap(mps1::MPState{Y1}, mps2::MPState{Y2}) where {Y1,Y2}
    Y = promote_type(Y1, Y2)
    overlap(convert(MPState{Y}, mps1), convert(MPSteat{Y}, mps2))
end

"""
    norm(mps)

calculates the norm of a matrix product state `mps` that is to
calculate the sqrt tensor contraction corresponding to ``⟨ψ|ψ⟩``.

"""
@inline norm(mps::MPState) = sqrt(real(overlap(mps, mps)))
