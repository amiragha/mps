randisometry(T, d1, d2; rng::AbstractRNG=GLOBAL_RNG) =
    d1 >= d2 ? Matrix(qr!(randn(rng, T, d1, d2)).Q) : Matrix(lq!(randn(rng, T, d1, d2)).Q)
randisometry(d1, d2; rng=rng) = randisometry(Float64, d1, d2, rng=rng)

"""
    svdtrunc(A [;maxdim, tol])

Performs the singular value decomposition and truncated the singular
values accoding to the tolerance `tol` or `maxdim` specified,
whichever is lower.

Return value is the tuple of truncated U, S, Vt. Here S is a diagonal
matrix instead of a vector.

"""
function svdtrunc(A::Matrix{T};
                  maxdim::Int=200,
                  tol::Float64=1.e-14) where{T<:Number}
    fact = svd(A, full=false)
    n = min(maxdim, sum(fact.S .> fact.S[1]*tol))

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
        (m < tol) && break
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
    entorpy(spectrum [; alpha])

calculate the entropy of a vector of numbers `spectrum`. The numbers
are assumed to be probabilities which means they are positive. The
numbers will be normalized so that they add up to 1. If `alpha=1` (the
default) it calculates the usual Shannon (Von-Neumann) entropy and if
`alpha > 1` calculates the Renyi entorpy.

If the input is a vector of vectors, then the entropy of each
individual vector is calculated.

"""
function entropy(spectrum::Vector{T};
                 alpha::Int64=1) where {T<:Number}

    s = spectrum ./ sum(spectrum)
    if alpha == 1
        return - sum(s .* log2.(s))
    else
        return log2(sum(s.^alpha))/(1-alpha)
    end
end

function entropy(spectrums::Vector{Vector{T}};
                 alpha::Int64=1) where {T<:Number}
    result = Vector{T}(undef, length(spectrums))
    for i in eachindex(spectrums)
        result[n] =  entropy(spectrums[i], alpha=alpha)
    end
    return result
end
