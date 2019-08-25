"""

`lx` is the length of the MPS or the number of sites.  `d` is the
physical on-site dimension, `dims` is the dimenions of bonds, starting
and ending with 1.  `mats` is the actual matrices, in three
dimenionsional tensor format, where the 2nd dimenisons is used for the
physical leg using the counterclockwise convention for the indeces.

"""
mutable struct MatrixProductState{T<:RLorCX}
    lx       :: Int64
    d        :: Int64
    dims     :: Vector{Int64}
    matrices :: Vector{Array{T, 3}}
    center   :: Int64
end

### CONSTRUCTORS

# x-direction! constructor
function MatrixProductState{T}(lx::Int64, d::Int64=2;
                               noise::Float64=0.0) where {T<:RLorCX}

    matrices = [ sqrt(1/d) * ones(T, 1, d, 1) + noise * randn(T, 1, d, 1)
                 for i=1:lx ]
    dims = ones(Int64, lx+1)
    canonicalize_at!(matrices, lx)
    MatrixProductState{T}(lx, d, dims, matrices, lx)
end

# constructor with intial configuration vector
function MatrixProductState{T}(lx::Int64, d::Int64,
                               initconf::Vector{Int64}) where{T<:RLorCX}
    @assert all((0 .<= initconf) .& (initconf .< d))
    matrices = [ zeros(T, 1, d, 1) for i=1:lx ]
    for site=1:lx
        matrices[site][1, initconf[site]+1, 1] = 1.
    end
    dims = ones(Int64, lx+1)
    canonicalize_at!(matrices, lx)
    MatrixProductState{T}(lx, d, dims, matrices, lx)
end

# constructor with intial configuration vector and noise
function MatrixProductState{T}(lx::Int64, d::Int64,
                               initconf::Vector{Int64},
                               maxdim::Int,
                               noise::Float64) where{T<:RLorCX}
    @assert all((0 .<= initconf) .& (initconf .< d))

    dims = ones(Int, lx+1)
    ### this can be wrapped into a function itself
    for l = 2:lx
        dims[l] = min(maxdim, dims[l-1]*d)
    end
    for l = lx:-1:1
        dims[l] = min(dims[l], dims[l+1]*d)
    end
    ### till here!

    matrices = [ noise * randn(T, dims[i], d, dims[i+1]) for i=1:lx ]

    for site=1:lx
        matrices[site][1, initconf[site]+1, 1] = 1.
    end

    canonicalize_at!(matrices, lx)
    MatrixProductState{T}(lx, d, dims, matrices, lx)
end

#constructor from a ketstate given in Ising basis
function MatrixProductState(lx::Int64, d::Int64,
                            ketstate::Vector{T};
                            svtruncation::Bool=false) where {T<:RLorCX}

    @assert length(ketstate) == d^lx
    matrices = Array{T,3}[]

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
        push!(matrices, reshape(U, dims[link-1], d, dims[link]))

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
    push!(matrices, reshape(U,dims[lx-1], d, dims[lx]))

    dims[lx+1] = 1
    push!(matrices, permutedims(reshape(Diagonal(S) * fact.Vt[1:n,:],
                                        dims[lx], dims[lx+1], d), [1,3,2]))

    return MatrixProductState{T}(lx, d, dims, matrices, lx)
end

function randmps(T::Type{<:RLorCX}, rng::AbstractRNG, lx::Int, d::Int, maxdim::Int)
    dims = Vector{Int}(undef, lx+1)
    dims[1] = 1
    dims[lx+1] = 1

    ### this can be wrapped into a function itself
    for l = 2:lx
        dims[l] = min(maxdim, dims[l-1]*d)
    end
    for l = lx:-1:1
        dims[l] = min(dims[l], dims[l+1]*d)
    end
    ### till here!

    matrices = Vector{Array{T, 3}}(undef, lx)
    for l = 1:lx
        Dr = dims[l+1]
        Dl = dims[l]
        matrices[l] = reshape(randisometry(T,rng, Dl, d*Dr), Dl, d, Dr)
    end
    return MatrixProductState{T}(lx, d, dims, matrices, lx)
end

### conversions
###############

function convert(::Type{MatrixProductState{ComplexF64}},
                 mps::MatrixProductState{Float64})
    MatrixProductState{ComplexF64}(mps.lx, mps.d, mps.dims,
                                   convert(Vector{Array{ComplexF64,3}}, mps.matrices),
                                   mps.center)
end

### TOOLS
#########

function truncate(svector::Vector{Float64};
                  maxdim::Int64=length(svector),
                  threshold::Float64=1.e-15)

    ## NOTE: S is assumed be all non-negative and sorted in
    ## descending order
    n = min(maxdim, sum(svector .> svector[1]*threshold))
    return svector[1:n], n, sum(svector[1:n])/sum(svector)
end

function normalize!(mps::MatrixProductState{T}) where{T<:RLorCX}
    A = mps.matrices[mps.center]
    if mps.center > div(mps.lx, 2)
        U,S,Vt = svd(reshape(A, prod(size(A)[1:2]), size(A,3)))
    else
        U,S,Vt = svd(reshape(A, size(A,1), prod(size(A)[2:3])))
    end
    normalize!(S)
    mps.matrices[mps.center] = reshape(U * S * Vt, size(A))
    S
end

"""
    canonicalize_at!(matrices, center)

canonicalize at `center` from the both ends of an MPS. This function
is mainly used during the initialization so one should not truncate
during the canonicalization process!

"""
function canonicalize_at!(matrices::Vector{Array{T, 3}},
                          center::Int64) where {T<:RLorCX}
    lx = length(matrices)
    @assert center > 0 && center < lx + 1

    for site=1:center-1
        isometrize_push_right!(matrices, site, svtruncation=false)
    end

    for site=lx:-1:center+1
        isometrize_push_left!(matrices, site, svtruncation=false)
    end
end

"""
   isometrize_push_right!(matrices, site)

perform a right step of the canonicalization procedure of a vector of
matrices at `site`. This is to perform an SVD on the matrices at
`site`, and then multiply the singluar values and right unitary
matricx ``V^{†}`` to matrices of `site+1`. This ensures the matrices
of `site` are left isometric, because U^{†}U = I.

"""
function isometrize_push_right!(matrices::Vector{Array{T, 3}},
                                site::Int64;
                                svtruncation::Bool=false) where {T<:RLorCX}
    lx = length(matrices)
    if site < lx
        a = matrices[site]
        dims = size(a)
        fact = svd(reshape(a, dims[1]*dims[2], dims[3]), full=false)

        if svtruncation
            S, n, ratio =  truncate(fact.S, threshold=1.e-15)
            U = fact.U[:,1:n]
            Vt = fact.Vt[1:n,:]
        else
            S, n, ratio = fact.S, length(fact.S), 1.
            U = fact.U
            Vt = fact.Vt
        end

        matrices[site] = reshape(U, dims[1], dims[2], n)

        @tensor matrices[site+1][i,d,j] := (Diagonal(S) * Vt)[i, k] *
            matrices[site+1][k,d,j]
    end
    return nothing
end

"""
    isometrize_push_left!(matrices, site)

perform a left step of the canonicalization procedure of a vector of
matrices at `site`. This is perform an SVD on the matrices at `site`,
and then multiply the left unitary matricx ``U`` and singular values
to matrices of `site-1`. This ensures the matrices of `site` are
right isometric, because VV^{†} = I.

"""
function isometrize_push_left!(matrices::Vector{Array{T, 3}},
                               site::Int64;
                               svtruncation::Bool=false) where {T<:RLorCX}
    if site > 0
        a = matrices[site]
        dims = size(a)
        fact = svd(reshape(a, dims[1], dims[2] * dims[3]), full=false)

        if svtruncation
            S, n, ratio =  truncate(fact.S)
            U = fact.U[:,1:n]
            Vt = fact.Vt[1:n,:]
        else
            S, n, ratio = fact.S, length(fact.S), 1.
            U = fact.U
            Vt = fact.Vt
        end

        #matrices[site] = permutedims(reshape(Vt, n, dims[3], dims[2]), [1,3,2])
        matrices[site] = reshape(Vt, n, dims[2], dims[3])

        @tensor matrices[site-1][i,d,j] := matrices[site-1][i,d,k] *
            (U * Diagonal(S))[k,j]
    end
    nothing
end
"""
    correctdims!(mps)

correctly store the bond dimensions of the MPS. This is needed when
some function manipulate the matrices of the MPS but don't see/fix the
dims field.

"""
function correctdims!(mps::MatrixProductState{T}) where {T<:RLorCX}
    dims = Int64[]
    push!(dims,size(mps.matrices[1])[1])
    for n=1:mps.lx
        push!(dims, size(mps.matrices[n])[3])
    end
    mps.dims = dims
    nothing
end

"""
    move_center!(mps, new_center)

move center of mps to a new location at `new_center`.

"""
function move_center!(mps::MatrixProductState{T},
                      new_center::Int64;
                      svnormalize::Bool=false) where {T<:RLorCX}
    lx= mps.lx
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

"""
    measure_1point(mps, operator, site)

local single site operator measurement at site `site`.
"""
function measure_1point(mps::MatrixProductState{T}, op::Matrix{T},
                        site::Int64) where {T<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert (site <= lx) && (site > 0)
    @assert size(op) == (d, d)

    move_center!(mps, site)
    mat = mps.matrices[site]

    @tensor v = scalar((mat[l,d',r] * conj(mat)[l,d,r]) * op[d,d'])
    #@tensor v = (mat[l,d',r] * conj(mat)[l,d,r]) * op[d,d']
    v
end

"""
    measure_1point(mps, operator)

local single site operator measurement at every site.
"""
function measure_1point(mps::MatrixProductState{T},
                        op::Matrix{T}) where {T<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert size(op) == (d, d)

    result = Vector{T}(undef, lx)

    left = ones(T, 1, 1)
    move_center!(mps, 1)
    for site = 1:lx
        mat = mps.matrices[site]
        @tensor v = scalar((left[lu, ld] * mat[lu, d', r]) * op[d, d'] * conj(mat)[ld, d, r])
        result[site] =  v
        @tensor left[ru, rd] := (left[lu, ld] * mat[lu, d, ru]) * conj(mat)[ld, d, rd]
    end
    result
end

## TODO: make these conversions more general -- learn more about the type system
function measure_1point(mps::MatrixProductState{ComplexF64},
                        op::Matrix{Float64})

    measure_1point(mps, convert(Matrix{ComplexF64}, op))
end

"""
    measure_2point(mps, op1, op2, site1, site2)

two local single site operators `op1` and `op2` measurement at sites,
`site1` and `site2`.

"""
function measure_2point(mps::MatrixProductState{T},
                        op1::Matrix{T}, op2::Matrix{T},
                        site1::Int64, site2::Int64) where {T<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert (site1 <= lx) && (site1 > 0)
    @assert (site2 <= lx) && (site2 > site1)
    @assert size(op1) == size(op2) == (d, d)

    move_center!(mps, site1)
    mat = mps.matrices[site1]

    @tensor left[ru,rd] := mat[l,d',ru] * op1[d,d'] * conj(mat)[l,d,rd]
    for site=site1+1:site2-1
        mat = mps.matrices[site]
        @tensor left[ru,rd] := left[lu, ld] * mat[lu,d,ru] * conj(mat)[ld,d,rd]
    end
    mat = mps.matrices[site2]
    @tensor v = scalar((left[lu, ld] * mat[lu,d',r] * conj(mat)[ld,d,r]) * op2[d,d'])
    v
end

"""
    measure_2point(mps, op1, op2)

two local single site operators `op1` and `op2` measurement at all
possible choices of two sites where site1 is less than site2. The
output the is the vector of measuremets having the measurement on
site1, site2 at index given by the `half_measurement_index` auxilary
function.

"""
function measure_2point(mps::MatrixProductState{T},
                        op1::Matrix{T}, op2::Matrix{T}) where {T<:RLorCX}
    d = mps.d
    lx = mps.lx
    @assert size(op1) == size(op2) == (d, d)

    result = T[]
    for site1=1:lx-1
        move_center!(mps, site1)
        mat = mps.matrices[site1]

        @tensor left[ru,rd] := mat[l,d',ru] * op1[d,d'] * conj(mat)[l,d,rd]
        for site=site1+1:lx-1
            mat = mps.matrices[site]
            @tensor v = scalar((left[lu, ld] * mat[lu,d',r] * conj(mat)[ld,d,r]) * op2[d,d'])
            push!(result, v)
            @tensor left[ru,rd] := left[lu, ld] * mat[lu,d,ru] * conj(mat)[ld,d,rd]
        end
        mat = mps.matrices[lx]
        @tensor v = scalar((left[lu, ld] * mat[lu,d',r] * conj(mat)[ld,d,r]) * op2[d,d'])
        push!(result, v)
    end
    result
end

## TODO: make these conversions more general -- learn more about the type system
function measure_2point(mps::MatrixProductState{ComplexF64},
                        op1::Matrix{Float64},
                        op2::Matrix{Float64})
    measure_2point(mps, convert(Matrix{ComplexF64}, op1), convert(Matrix{ComplexF64}, op2))
end


"""
    half_measurement_index(lx, m, n)

Returns the index corresponding to the half_measurement of lattice
size `lx`. That means the index in the vector When `n > m` for all `m` and
`n`.

"""
function half_measurement_index(lx::Int64, m::Int64, n::Int64)
    @assert 1 <= m < n <= lx
    # more readable version : div(lx*(lx-1),2) - div((lx-m+1)*(lx-m),2) + (n-m)
    div((2*lx-m)*(m-1), 2) + (n-m)
end

"""
    measure_mpo(mps, mpo)

The output is the result of measurement of a full MPO that is a
number.

"""
function measure_mpo(mps::MatrixProductState{T},
                     mpo::MatrixProductOperator{T}) where{T<:RLorCX}
    @assert mps.d == mpo.d && mps.lx == mpo.lx

    left = ones(T, 1, 1, 1)

    for site=1:mps.lx
        mat = mps.matrices[site]
        ten = mpo.tensors[site]
        @tensor left[ru,rm,rd] := (left[lu,lm,ld] * mat[lu,d',ru]) *
            ten[lm, d, rm, d'] * conj(mat)[ld,d,rd]
    end
    left[1,1,1]
end

"""
    entanglemententropy(A)

Measure the entanglement entropy for an MPS at a given cut `l` or if
ommited at every bond of the MPS.

"""
function entanglemententropy(mps::MatrixProductState{T};
                             alpha::Int=1) where{T}

    lx = mps.lx
    result = Vector{Float64}(undef, lx-1)
    move_center!(mps, 1)
    A = mps.matrices[1]

    for l = 1:lx-1
        U, S, Vt = svd(reshape(A, size(A, 1)*size(A,2), size(A,3)))
        mps.matrices[l] = reshape(U, size(A))
        result[l] = entropy(S, alpha=alpha)
        @tensor A[l,o,r] := (Diagonal(S)*Vt)[l,m] * mps.matrices[l+1][m,o,r]
    end

    mps.matrices[lx] = A
    mps.center = lx

    result
end


### USEFULL FUNCTIONS FOR TESTING
"""
    mps2ketstate(mps)

make a ketstate from a `mps` by multiplication the matrices
corresponding to each Ising configuration. Note that this function is
very expensive and is only meant for testing.

"""
function mps2ketstate(mps::MatrixProductState{T}) where {T<:RLorCX}

    lx = mps.lx
    d = mps.d

    ketstate = zeros(T, d^lx)

    for ketindex::Int64 = 0:d^lx-1
        resolve::Int64 = ketindex
        amplitude = T[1]
        ## TODO: find a better way!
        for i::Int64 = lx:-1:1
            amplitude = reshape(mps.matrices[i][:,(resolve % d) + 1,:],
                                mps.dims[i], mps.dims[i+1]) * amplitude
            resolve = div(resolve, d)
        end
        ketstate[ketindex+1] = amplitude[1]
    end

    return ketstate
end

"""
    apply_2siteoperator!(mps, l, operator, maxdim, pushto)

applies the `operator` which is `d x d x d x d` tensor to site `l` and
`l+1` of the `mps`. The order of indeces start from the bottom left
(resulting on site l) and are counterclockwise, so it is bottom left,
bottom right, top right, top left.

The `maxdim` operator chooses the max possible size of dimension of
the new mps at bond between `l` and `l+1`, The singular values are
push to either left `:L` or right `:R` (default) matrices using the
argument `pushto`.

"""
function apply_2siteoperator!(mps      ::MatrixProductState{T},
                              l        ::Int64,
                              operator ::Array{T, 4};
                              maxdim  ::Int64=mps.dims[l+1],
                              pushto  ::Symbol=:R) where {T<:RLorCX}

    @assert mps_dims_are_consistent(mps)
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

    @tensor R[a,i,j,b] := operator[i,j,l,k] * (one[a,k,c] * two[c,l,b])

    fact = svd(reshape(R, dim_l*d, d*dim_r), full=false)

    #println(fact[:S])
    S, n, ratio = truncate(fact.S, maxdim=maxdim)
    #println(S)

    U = fact.U[:,1:n]
    Vt = fact.Vt[1:n,:]
    ## QQQ? do we need to normalize S here?

    mps.dims[l+1] = n

    if (pushto == :R)
        mps.matrices[l] = reshape(U, dim_l, d, n)
        mps.matrices[l+1] = reshape(Diagonal(S) * Vt, n, d, dim_r)
        mps.center = l+1
    elseif (pushto == :L)
        mps.matrices[l] = reshape(U * Diagonal(S), dim_l, d, n)
        mps.matrices[l+1] = reshape(Diagonal(S) * Vt, n, d, dim_r)
        mps.center = l
    else
        error("invalid push_to :", pushto)
    end
    nothing
end

function apply_2siteoperator!(mps      ::MatrixProductState{ComplexF64},
                              l        ::Int64,
                              op ::Array{Float64, 4};
                              max_dim  ::Int64=mps.dims[l+1],
                              push_to  ::Symbol=:R)

    apply_2siteoperator!(mps, l, convert(Array{ComplexF64, 4}, op),
                         max_dim, push_to)
end

function twosite_tensor(op1::Matrix{T}, op2::Matrix{T}) where {T<:RLorCX}
    a, b = size(op1)
    c, d = size(op2)
    @assert a == b && c == d
    permutedims(reshape(kron(op1, op2),a,c,a,c), [1,2,4,3])
end

"""
    apply_1siteoperator!(mps, l, op)

applies the `op` which is `d x d` tensor to site `l` of the `mps`. The
order of indeces are bottom, top.

"""
function apply_1siteoperator!(mps      ::MatrixProductState{T},
                              l        ::Int64,
                              op       ::Matrix{T}) where {T<:RLorCX}

    @assert 0 < l <= mps.lx

    d = mps.d
    move_center!(mps, l)
    A = mps.matrices[l]
    @tensor A[l,o,r] := A[l,o',r] * op[o,o']
    mps.matrices[l] = A

    nothing
end

### Multi-MPS functions
#######################

"""
    overlap(mps1, mps2)

calculates the overlap between two matrix product states `mps1` and
`mps2` that is to run the tensor contraction corresponding to
``⟨ψ_2|ψ_1⟩``. Note that it is not divided by the norm of the two
MPSs, so it returned value is the overlap of the two states multiplied
by the norm of each.

"""
function overlap(mps1::MatrixProductState{T},
                 mps2::MatrixProductState{T}) where {T<:RLorCX}

    @assert mps1.d == mps2.d && mps1.lx == mps2.lx

    left = ones(T, 1, 1)

    for site=1:mps1.lx
        mat1 = mps1.matrices[site]
        mat2 = mps2.matrices[site]
        @tensor left[ru,rd] := (left[lu,ld] * mat1[lu,d,ru]) * conj(mat2)[ld,d,rd]
    end
    left[1,1]
end

overlap(mps1::MatrixProductState{ComplexF64}, mps2::MatrixProductState{Float64}) =
    overlap(mps1, convert(MatrixProductState{ComplexF64}, mps2))
overlap(mps1::MatrixProductState{Float64}, mps2::MatrixProductState{ComplexF64}) =
    overlap(convert(MatrixProductState{ComplexF64}, mps1), mps2)

"""
    norm2(mps)

calculates the norm of a matrix product state `mps` that is to
calculate the tensor contraction corresponding to ``⟨ψ|ψ⟩``.

"""
norm2(mps::MatrixProductState{T}) where {T<:RLorCX} = real(overlap(mps, mps))

### useful tools
################

"""
    display_matrices(mps, range)

display all the matrices of the `mps` in `range`. This function is
useful for testing or educational purposes.

"""
function display_matrices(mps::MatrixProductState{T},
                          range::UnitRange{Int64}=1:mps.length) where {T<:RLorCX}
    for site=range
        for s=1:mps.d
            display("Matrix $site $(s-1)")
            display(mps.matrices[site][:,s,:])
        end
    end
end

"""
    mps_dims_are_consistent(mps)

check if are dimensions are consistent in an mps. This is made for
testing and double checks; in principle all operations must not break
the consistency of the bond dimensions of MPS.

"""
function mps_dims_are_consistent(mps::MatrixProductState{T}) where {T<:RLorCX}
    for n=1:mps.lx-1
        dim1 = size(mps.matrices[n])[3]
        dim2 = size(mps.matrices[n+1])[1]
        if dim1 != dim2 || dim1 != mps.dims[n+1]
            return false
        end
    end
    return size(mps.matrices[mps.lx])[3] == mps.dims[mps.lx+1]
end
