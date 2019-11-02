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

"""
    truncatedsizes(S [; maxdim, tol])

for a vector of float vectors `S` of singular values (or basically
it is just assumed that each sector is sorted descendingly), returns
the vector of new sizes when the singular values are to be truncated
to the `maxdim` largest numbers or when `tol` is reached whichever is
lower.

"""
function truncatedsizes(S::Vector{Vector{Float64}};
                        maxdim::Int=200,
                        tol::Float64=1.e-14)

    pointers = ones(Int, length(S))

    c = [S[i][pointers[i]] for i in eachindex(S)]
    for n=1:min(maxdim, sum([length(s) for s in S]))
        #println(c, pointers)
        m, index = findmax(c)
        if (m < tol)
            n==1 && error("The largest singular value is smaller than tol! $m < $tol")
            break
        end
        if pointers[index] >= size(S[index], 1)
            pointers[index] += 1
            c[index] = 0.0
        else
            pointers[index] += 1
            c[index] = S[index][pointers[index]]
        end
    end
    pointers .- 1
end

######NOTE
###TODO: this needs to be tested!
function svdtrunc(A::AbstractSymTensor;
                  maxdim::Int=200,
                  tol::Float64=1.e-14)
    numoflegs(A) == 2 ||
        error("svd only defined for matrix like objects N = ", numoflegs(A))
    signs(A.legs) == (+1, -1) ||
        error("svdsym only accepts a SymMatrix (+1,-1) but ", signs(A.legs))
    T = eltype(A)

    n_sects = length(A.sects)
    sects_U  = Vector{Tuple{Int, Int}}(undef, n_sects)
    sects_S  = Vector{Tuple{Int, Int}}(undef, n_sects)
    sects_Vt = Vector{Tuple{Int, Int}}(undef, n_sects)

    blks_U  = Vector{Matrix{T}}(undef, n_sects)
    blks_S  = Vector{Vector{Float64}}(undef, n_sects)
    blks_Vt = Vector{Matrix{T}}(undef, n_sects)

    for idx in eachindex(A.sects)
        c1, c2 = A.sects[idx]
        sects_U[idx]  = (c1, c1)
        sects_S[idx]  = (c1, c2)
        sects_Vt[idx] = (c2, c2)

        fact = svd(A.nzblks[idx], full=false)
        blks_U[idx] = fact.U
        blks_S[idx] = fact.S
        blks_Vt[idx] = fact.Vt
    end

    ns = truncatedsizes(blks_S, maxdim=maxdim, tol=tol)
    indices = findall(x->x>0, ns)
    #n_sects_new = sum(ns .> 0)

    sects_U  = sects_U[indices]
    sects_S  = sects_S[indices]
    sects_Vt = sects_Vt[indices]

    blks_U = [blks_U[index][:, 1:ns[index]] for index in indices]
    blks_S = [blks_S[index][1:ns[index]] for index in indices]
    blks_Vt = [blks_Vt[index][1:ns[index], :] for index in indices]

    middledims = ns[indices]

    ls, rs = zip(sects_S...)
    lchrs = [ls...]
    rchrs = [rs...]

    Uleg2  = STLeg(-1, lchrs, middledims)
    Sleg1  = STLeg(+1, lchrs, middledims)
    Sleg2  = STLeg(-1, rchrs, middledims)
    Vtleg1 = STLeg(+1, rchrs, middledims)

    U  = SymMatrix{T}(0, (A.legs[1], Uleg2), sects_U, blks_U)
    S  = SymDiagonal{Float64}(A.charge, (Sleg1, Sleg2), sects_S, [Diagonal(blk) for blk in blks_S])
    Vt = SymMatrix{T}(0, (Vtleg1, A.legs[2]), sects_Vt, blks_Vt)

    U, S, Vt
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
    sum(P) ≈ 1 || println("P vector doesn't add up to one : $(sum(P))")
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
