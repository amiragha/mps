# equal probability constructor
function MPS{T}(lx::Int, d::Int=2;
                noise::Float64=0.0) where {T<:Number}

    As = [ sqrt(1/d) * ones(T, 1, d, 1) + noise * randn(T, 1, d, 1)
           for i=1:lx ]
    mps = MPS{T}(As, 0)
    center_at!(mps, 1)
    mps
end

# # constructor with intial configuration vector
# function MatrixProductState{T}(lx::Int64, d::Int64,
#                                initconf::Vector{Int64}) where{T<:RLorCX}
#     @assert all((0 .<= initconf) .& (initconf .< d))
#     matrices = [ zeros(T, 1, d, 1) for i=1:lx ]
#     for site=1:lx
#         matrices[site][1, initconf[site]+1, 1] = 1.
#     end
#     dims = ones(Int64, lx+1)
#     canonicalize_at!(matrices, lx)
#     MatrixProductState{T}(lx, d, dims, matrices, lx)
# end

# # constructor with intial configuration vector and noise
# function MatrixProductState{T}(lx::Int64, d::Int64,
#                                initconf::Vector{Int64},
#                                maxdim::Int,
#                                noise::Float64) where{T<:RLorCX}
#     @assert all((0 .<= initconf) .& (initconf .< d))

#     dims = ones(Int, lx+1)
#     ### this can be wrapped into a function itself
#     for l = 2:lx
#         dims[l] = min(maxdim, dims[l-1]*d)
#     end
#     for l = lx:-1:1
#         dims[l] = min(dims[l], dims[l+1]*d)
#     end
#     ### till here!

#     matrices = [ noise * randn(T, dims[i], d, dims[i+1]) for i=1:lx ]

#     for site=1:lx
#         matrices[site][1, initconf[site]+1, 1] = 1.
#     end

#     canonicalize_at!(matrices, lx)
#     MatrixProductState{T}(lx, d, dims, matrices, lx)
# end

#constructor from a ketstate given in Ising basis
function MPS(lx::Int, d::Int,
             ketstate::Vector{T};
             svtruncation::Bool=false) where {T<:Number}

    @assert length(ketstate) == d^lx
    mps = MPS{T}()
    #As = Array{T,3}[]

    dims = zeros(Int64, lx+1)
    dims[1] = 1
    rdim = d^(lx-1)

    # state_matrix is now : dims[1]*d x rdim
    state_matrix = transpose(reshape(ketstate, rdim, dims[1]*d))

    for link=2:lx-1
        fact = svd(state_matrix, full=false)

        if svtruncation
            S, n, ratio =  truncate(fact.S, threshold=1.e-15)
        else
            S, n, ratio = fact.S, length(fact.S), 1.
        end
        dims[link] = n

        # U is : dims[link-1] * d, dims[link]
        U = fact.U[:,1:n]
        push!(mps, reshape(U, dims[link-1], d, dims[link]))

        # S*Vt is : dims[link] x rdim where rdim is the new rdim x d
        rdim = div(rdim , d)
        state_matrix = reshape(permutedims(
            reshape(Diagonal(S)*fact.Vt[1:n,:], dims[link], rdim, d),
            [1,3,2]), dims[link]*d, rdim)
        # final results is dims[link]*d x rdim
    end

    fact = svd(state_matrix, full=false)

    if svtruncation
        S, n, ratio =  truncate(fact.S)
    else
        S, n, ratio = fact.S, length(fact.S), 1.
    end

    dims[lx] = n

    U = fact.U[:,1:n]
    push!(mps, reshape(U, dims[lx-1], d, dims[lx]))

    dims[lx+1] = 1
    push!(mps, permutedims(reshape(Diagonal(S) * fact.Vt[1:n,:],
                                   dims[lx], dims[lx+1], d), [1,3,2]))

    return mps
end

"""
    mps2ketstate(mps)

make a ketstate from a `mps` by multiplication the matrices
corresponding to each Ising configuration. Note that this function is
very expensive and is only meant for testing.

"""
function mps2ketstate(mps::MPS{T}) where {T<:Number}

    lx = length(mps)
    ds = [dim(sitespace(mps, l)) for l=1:lx]

    ketstate = zeros(T, prod(ds))

    for ketindex::Int = 0:prod(ds)-1
        resolve::Int = ketindex
        amplitude = T[1]
        ## TODO: find a better way!
        for i::Int64 = lx:-1:1
            d = ds[i]
            A = mps.As[i]
            amplitude = reshape(A[:,(resolve % d) + 1,:],
                                size(A, 1), size(A, 3)) * amplitude
            resolve = div(resolve, d)
        end
        ketstate[ketindex+1] = amplitude[1]
    end
    ketstate
end

# function randmps(T::Type{<:RLorCX}, rng::AbstractRNG, lx::Int, d::Int, maxdim::Int)
#     dims = Vector{Int}(undef, lx+1)
#     dims[1] = 1
#     dims[lx+1] = 1

#     ### this can be wrapped into a function itself
#     for l = 2:lx
#         dims[l] = min(maxdim, dims[l-1]*d)
#     end
#     for l = lx:-1:1
#         dims[l] = min(dims[l], dims[l+1]*d)
#     end
#     ### till here!

#     matrices = Vector{Array{T, 3}}(undef, lx)
#     for l = 1:lx
#         Dr = dims[l+1]
#         Dl = dims[l]
#         matrices[l] = reshape(randisometry(T,rng, Dl, d*Dr), Dl, d, Dr)
#     end
#     return MatrixProductState{T}(lx, d, dims, matrices, lx)
# end

# ### conversions
# ###############

# function convert(::Type{MatrixProductState{ComplexF64}},
#                  mps::MatrixProductState{Float64})
#     MatrixProductState{ComplexF64}(mps.lx, mps.d, mps.dims,
#                                    convert(Vector{Array{ComplexF64,3}}, mps.matrices),
#                                    mps.center)
# end

# ### TOOLS
# #########

# function truncate(svector::Vector{Float64};
#                   maxdim::Int64=length(svector),
#                   threshold::Float64=1.e-15)

#     ## NOTE: S is assumed be all non-negative and sorted in
#     ## descending order
#     n = min(maxdim, sum(svector .> svector[1]*threshold))
#     return svector[1:n], n, sum(svector[1:n])/sum(svector)
# end

# function normalize!(mps::MatrixProductState{T}) where{T<:RLorCX}
#     A = mps.matrices[mps.center]
#     if mps.center > div(mps.lx, 2)
#         U,S,Vt = svd(reshape(A, prod(size(A)[1:2]), size(A,3)))
#     else
#         U,S,Vt = svd(reshape(A, size(A,1), prod(size(A)[2:3])))
#     end
#     normalize!(S)
#     mps.matrices[mps.center] = reshape(U * S * Vt, size(A))
#     S
# end

# """
#     canonicalize_at!(matrices, center)

# canonicalize at `center` from the both ends of an MPS. This function
# is mainly used during the initialization so one should not truncate
# during the canonicalization process!

# """
# function canonicalize_at!(matrices::Vector{Array{T, 3}},
#                           center::Int64) where {T<:RLorCX}
#     lx = length(matrices)
#     @assert center > 0 && center < lx + 1

#     for site=1:center-1
#         isometrize_push_right!(matrices, site, svtruncation=false)
#     end

#     for site=lx:-1:center+1
#         isometrize_push_left!(matrices, site, svtruncation=false)
#     end
# end

# """
#    isometrize_push_right!(matrices, site)

# perform a right step of the canonicalization procedure of a vector of
# matrices at `site`. This is to perform an SVD on the matrices at
# `site`, and then multiply the singluar values and right unitary
# matricx ``V^{†}`` to matrices of `site+1`. This ensures the matrices
# of `site` are left isometric, because U^{†}U = I.

# """
# function isometrize_push_right!(matrices::Vector{Array{T, 3}},
#                                 site::Int64;
#                                 svtruncation::Bool=false) where {T<:RLorCX}
#     lx = length(matrices)
#     if site < lx
#         a = matrices[site]
#         dims = size(a)
#         fact = svd(reshape(a, dims[1]*dims[2], dims[3]), full=false)

#         if svtruncation
#             S, n, ratio =  truncate(fact.S, threshold=1.e-15)
#             U = fact.U[:,1:n]
#             Vt = fact.Vt[1:n,:]
#         else
#             S, n, ratio = fact.S, length(fact.S), 1.
#             U = fact.U
#             Vt = fact.Vt
#         end

#         matrices[site] = reshape(U, dims[1], dims[2], n)

#         @tensor matrices[site+1][i,d,j] := (Diagonal(S) * Vt)[i, k] *
#             matrices[site+1][k,d,j]
#     end
#     return nothing
# end

# """
#     isometrize_push_left!(matrices, site)

# perform a left step of the canonicalization procedure of a vector of
# matrices at `site`. This is perform an SVD on the matrices at `site`,
# and then multiply the left unitary matricx ``U`` and singular values
# to matrices of `site-1`. This ensures the matrices of `site` are
# right isometric, because VV^{†} = I.

# """
# function isometrize_push_left!(matrices::Vector{Array{T, 3}},
#                                site::Int64;
#                                svtruncation::Bool=false) where {T<:RLorCX}
#     if site > 0
#         a = matrices[site]
#         dims = size(a)
#         fact = svd(reshape(a, dims[1], dims[2] * dims[3]), full=false)

#         if svtruncation
#             S, n, ratio =  truncate(fact.S)
#             U = fact.U[:,1:n]
#             Vt = fact.Vt[1:n,:]
#         else
#             S, n, ratio = fact.S, length(fact.S), 1.
#             U = fact.U
#             Vt = fact.Vt
#         end

#         #matrices[site] = permutedims(reshape(Vt, n, dims[3], dims[2]), [1,3,2])
#         matrices[site] = reshape(Vt, n, dims[2], dims[3])

#         @tensor matrices[site-1][i,d,j] := matrices[site-1][i,d,k] *
#             (U * Diagonal(S))[k,j]
#     end
#     nothing
# end
# """
#     correctdims!(mps)

# correctly store the bond dimensions of the MPS. This is needed when
# some function manipulate the matrices of the MPS but don't see/fix the
# dims field.

# """
# function correctdims!(mps::MatrixProductState{T}) where {T<:RLorCX}
#     dims = Int64[]
#     push!(dims,size(mps.matrices[1])[1])
#     for n=1:mps.lx
#         push!(dims, size(mps.matrices[n])[3])
#     end
#     mps.dims = dims
#     nothing
# end

# """
#     move_center!(mps, new_center)

# move center of mps to a new location at `new_center`.

# """
# function move_center!(mps::MatrixProductState{T},
#                       new_center::Int64;
#                       svnormalize::Bool=false) where {T<:RLorCX}
#     lx= mps.lx
#     @assert new_center > 0 && new_center <= lx

#     center = mps.center
#     if new_center > center
#         for p=center:new_center-1
#             isometrize_push_right!(mps.matrices, p)
#         end
#     elseif new_center < center
#         for p=center:-1:new_center+1
#             isometrize_push_left!(mps.matrices, p)
#         end
#     end
#     correctdims!(mps)
#     mps.center = new_center
#     nothing
# end

## TODO: combine these with the symmetric cases!
# specific measure functions
function measure(mps::MPS{T},
                 op::Matrix{T}, l::Int) where {T}
    lx = length(mps)
    @assert (l <= lx) && (l > 0)
    @assert sitespace(mps, l) == space(op, 2) == dual(space(op, 1))

    center_at!(mps, l)
    A = mps.As[l]
    @tensor v = scalar((A[l,d',r] * conj(A)[l,d,r]) * op[d,d'])
    Dict{Int, T}(l=>v)
end

function measure(mps::MPS{T},
                 op::Matrix{T}) where {T}
    Vd = space(op, 1)
    @assert space(op, 2) == Vd

    lx = length(mps)
    values = Dict{Int, T}()
    center_at!(mps, 1)
    L = ones(T, 1, 1)

    for x = 1:lx
        A = mps.As[x]
        @tensor v = scalar((L[lu, ld] * A[lu, d', r]) * op[d, d'] * conj(A)[ld, d, r])
        values[x] = v
        @tensor L[ru, rd] := (L[lu, ld] * A[lu, d, ru]) * conj(A)[ld, d, rd]
    end
    values
end
function measure(mps::MPS{T1},
                 op::Matrix{T2}) where {T1,T2}
    measure(mps, convert(Matrix{T1}, op))
end


function measure(mps::MPS{T},
                 op1::Matrix{T},
                 op2::Matrix{T},
                 x1::Int,
                 x2::Int) where {T}
    lx = length(mps)
    @assert (0 < x1 < x2) && (x2 <= lx)

    Vd1 = space(op1, 1)
    @assert space(op1, 2) == Vd1
    Vd2 = space(op2, 1)
    @assert space(op2, 2) == Vd2

    values = Dict{NTuple{2, Int}, T}()
    center_at!(mps, x1)
    A = mps.As[x1]

    @tensor L[ru,rd] := A[l,d',ru] * op1[d,d'] * conj(A)[l,d,rd]
    for x=x1+1:x2-1
        A = mps.As[site]
        @tensor L[ru,rd] := L[lu, ld] * A[lu,d,ru] * conj(A)[ld,d,rd]
    end
    A = mps.As[x2]
    @tensor v = scalar((L[lu, ld] * A[lu,d',r] * conj(A)[ld,d,r]) * op2[d,d'])
    values[(x1, x2)] = v
end

function measure(mps::MPS{T},
                 op1::Matrix{T},
                 op2::Matrix{T}) where {T}
    lx = length(mps)
    Vd1 = space(op1, 1)
    @assert space(op1, 2) == dual(Vd1)
    Vd2 = space(op2, 1)
    @assert space(op2, 2) == dual(Vd2)

    values = Dict{NTuple{2,Int}, T}()
    for x1=1:lx-1
        center_at!(mps, x1)
        A = mps.As[x1]
        @tensor L[ru,rd] := A[l,d',ru] * op1[d,d'] * conj(A)[l,d,rd]
        for x=x1+1:lx-1
            A = mps.As[x]
            @tensor v = scalar((L[lu, ld] * A[lu,d',r] * conj(A)[ld,d,r]) * op2[d,d'])
            values[(x1, x)] = v
            @tensor L[ru,rd] := L[lu, ld] * A[lu,d,ru] * conj(A)[ld,d,rd]
        end
        A = mps.As[lx]
        @tensor v = scalar((L[lu, ld] * A[lu,d',r] * conj(A)[ld,d,r]) * op2[d,d'])
        values[x1, lx] = v
    end
    values
end

function measure(mps::MPS{T1},
                 op1::Matrix{T2},
                 op2::Matrix{T3}) where {T1,T2,T3}
    measure(mps, convert(Matrix{T1}, op1), convert(Matrix{T1}, op2))
end

function measure(mps::MPS{T},
                 mpo::MPO{T}) where {T}

    @assert length(mpo) == length(mps)
    values = Dict{Int, T}()
    L = ones(T, 1, 1, 1)
    for x=1:length(mps)
        A = mps.As[x]
        W = mpo.Ws[x]
        @tensor L[ru,rm,rd] := (L[lu,lm,ld] * A[lu,d',ru]) *
            W[lm, d, rm, d'] * conj(A)[ld,d,rd]
    end
    values[1] = L[1,1,1]
    values
end

function apply!(mps         :: MPS{T},
                op          :: Array{T,4},
                x           :: Int;
                maxdim      :: Int=bonddim(mps, x+1),
                pushto      :: Symbol=:R,
                svnormalize :: Bool=false) where {T}

    @boundscheck 0 < x < length(mps) || throw()

    if mps.center < x
        center_at!(mps, x)
    elseif mps.center > x+1
        center_at!(mps, x+1)
    end

    A1 = mps.As[x]
    A2 = mps.As[x+1]

    @tensor R[a,i,j,b] := op[i,j,k,l] * (A1[a,k,c] * A2[c,l,b])

    u,s,v = svdtrunc(reshape(R, prod(size(R)[1:2]), prod(size(R)[3:4])),
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

function apply!(mps         :: MPS{T1},
                op          :: Array{T2,4},
                x           :: Int;
                maxdim      :: Int=bonddim(mps, x+1),
                pushto      :: Symbol=:R,
                svnormalize :: Bool=false) where {T1,T2}

    T = promote_type(T1, T2)
    apply!(convert(MPS{T}, mps),
           convert(Array{T,4}, op),
           x, maxdim=maxdim, pushto=pushto, svnormalize=svnormalize)
end

function apply!(mps      :: MPS{T},
                op       :: Array{T, 2},
                x        :: Int) where {T}

    @boundscheck 0 < x < length(mps) || throw()

    center_at!(mps, x)

    A = mps.As[x]
    @tensor A[l,o,r] := A[l,o',r] * op[o,o']
    mps.As[x] = A
    mps
end
