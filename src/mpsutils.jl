"""
    half_measurement_index(lx, m, n)

Returns the index corresponding to the half_measurement of lattice
size `lx`. That means the index in the vector When `n > m` for all `m` and
`n`.

"""
# more readable version : div(lx*(lx-1),2) - div((lx-m+1)*(lx-m),2) + (n-m)
half_measurement_index(lx::Int, m::Int, n::Int) =
    1 <= m < n <= lx ? div((2*lx-m)*(m-1), 2) + (n-m) : error("1 <= $m < $n <= $lx")

"""
    randisometry(d1, d2)

Returns a uniformly chosen random isometry matrix of size `d1` × `d2`
"""
randisometry(T,rng::AbstractRNG, d1, d2) =
    d1 >= d2 ? Matrix(qr!(randn(rng, T, d1, d2)).Q) : Matrix(lq!(randn(rng, T, d1, d2)).Q)
randisometry(rng, d1, d2) = randisometry(Float64, rng, d1, d2)

"""
    svdtrunc(A [;maxdim, tol])

Performs the singular value decomposition (SVD) and truncated the
singular values (SV) accoding to the tolerance `tol` or `maxdim`
specified, whichever results in a lower number of SVs.

Return value is the tuple of truncated `U`, `S`, `Vt`. Where `U`, and
`Vt` are left and right isometries and `S` is a Diagonal matrix
containing the SVs.

"""
function svdtrunc(A::Matrix{T};
                  maxdim::Int=200,
                  tol::Float64=1.e-14) where{T<:Number}
    fact = svd(A, full=false)
    n = min(maxdim, sum(fact.S .> fact.S[1]*tol))
    n == 0 && error("The largest singular value is smaller than tol! $m < $tol")

    fact.U[:, 1:n], Diagonal(fact.S[1:n]), fact.Vt[1:n, :]
end

"""
    picklargests(Ss, maxdim)

for a vector of descending sorted vectors `Ss`, pick the `maxdim`
largest numbers and discard the rest.

"""
##TODO: make this better and find the best version!
function picklargests(Ss::Vector{Vector{Float64}},
                      maxdim::Int=200, tol::Float64=1.e-14)
    #vs = Vector{Vector{Float64}}(Float64[], length(Ss))
    vs = [Float64[] for i=1:length(Ss)]

    pointers = ones(Int, length(Ss))
    c = [Ss[i][pointers[i]] for i in eachindex(Ss)]

    for n=1:min(maxdim, sum([length(s) for s in Ss]))
        #println(c, pointers)
        m, index = findmax(c)
        if (m < tol)
            n==1 && error("The largest singular value is smaller than tol! $m < $tol")
            break
        end
        push!(vs[index], Ss[index][pointers[index]])
        if pointers[index] == length(Ss[index])
            #println("here")
            c[index] = 0.0
        else
            pointers[index] += 1
            c[index] = Ss[index][pointers[index]]
        end
    end
    return vs
end

######NOTE
###TODO: this needs to be tested!
function svdtrunc(A::SymTensor{Tv, 2};
                  maxdim::Int=200,
                  tol::Float64=1.e-14) where {Tv<:Number}
    @assert signs(A.legs) == (+1, -1)

    sects_U = Tuple{Int, Int}[]
    nzblks_U = Matrix{Tv}[]
    sects_S = Tuple{Int, Int}[]
    nzblks_S = Vector{Float64}[]
    nzblks_Vt = Matrix{Tv}[]
    sects_Vt = Tuple{Int, Int}[]
    middle_dims = Int[]

    for idx in eachindex(A.sects)
        sect = A.sects[idx]
        push!(sects_U, (sect[1], sect[1]))
        push!(sects_S, (sect[1], sect[2]))
        push!(sects_Vt, (sect[2], sect[2]))

        nzblk = A.nzblks[idx]
        fact = svd(nzblk, full=false)
        push!(nzblks_U, fact.U)
        push!(nzblks_S, fact.S)
        push!(nzblks_Vt, fact.Vt)

        push!(middle_dims, length(fact.S))
    end

    nzblks_S = picklargests(nzblks_S, maxdim, tol)
    vblkindexes = Int[]
    for i in eachindex(nzblks_S)
        n = size(nzblks_S[i], 1)
        if n != 0
            push!(vblkindexes, i)
        else
            continue
        end
        middle_dims[i] = n
        #println(n)
        #println(nzblks_U[i])
        nzblks_U[i] = nzblks_U[i][:, 1:n]
        nzblks_Vt[i] = nzblks_Vt[i][1:n, :]
    end
    middle_dims = middle_dims[vblkindexes]
    nzblks_U = nzblks_U[vblkindexes]
    nzblks_Vt = nzblks_Vt[vblkindexes]
    nzblks_S = nzblks_S[vblkindexes]
    sects_U = sects_U[vblkindexes]
    sects_S = sects_S[vblkindexes]
    sects_Vt = sects_Vt[vblkindexes]

    ls, rs = zip(sects_S...)
    lchrs = [ls...]
    rchrs = [rs...]

    Uleg2 = STLeg(-1, lchrs, middle_dims)
    Sleg1 = STLeg(+1, rchrs, middle_dims)
    Sleg2 = STLeg(-1, rchrs, middle_dims)
    Vtleg1 = STLeg(+1, lchrs, middle_dims)

    U = SymTensor(0, (A.legs[1], Uleg2), sects_U, nzblks_U)
    S = SymTensor(A.charge, (Sleg1, Sleg2), sects_S, [Diagonal(blk) for blk in nzblks_S])
    Vt = SymTensor(0, (Vtleg1, A.legs[2]), sects_Vt, nzblks_Vt)
    return U, S, Vt
end

"""
    entorpy(P [; alpha])

Calculate the entropy of a vector of numbers `P`. The numbers are
assumed to form a probability distribution which means they add up to
one.

 If `alpha=1` (the default) it calculates the usual Shannon
(Von-Neumann) entropy and if `alpha > 1` calculates the Renyi entorpy.

"""
function entropy(P::Vector{T};
                 alpha::Int=1) where {T<:Number}

    any(P .< 0.0) && error("Negative value in P vector ", P)
    sum(P) ≈ 1 || error("P vector doesn't add up to one : $(sum(P))")
    if alpha == 1
        return -sum(P .* log.(P))
    else
        return log(sum(P.^alpha))/(1-alpha)
    end
end

function entropy(spectrums::Vector{Vector{T}};
                 alpha::Int=1) where {T<:Number}
    result = Vector{T}(undef, length(spectrums))
    for i in eachindex(spectrums)
        result[i] =  entropy(spectrums[i], alpha=alpha)
    end
    return result
end
