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
                  maxdim::Int=size(A, 2),
                  tol::Float64=1.e-12) where{T<:Number}
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

        m, index = findmax(c)
        if (m < tol)
            n==1 && error("The largest singular value is smaller than tol! $m < $tol")
            break
        end
        push!(vs[index], Ss[index][pointers[index]])
        if pointers[index] == length(Ss[index])

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
                  maxdim::Int=maximum(size(A)),
                  tol::Float64=1.e-12)
    T = eltype(A)
    S = vtype(A)

    n_sects = length(A.blocks)
    sects_u  = Vector{Sector{S, 2}}(undef, n_sects)
    sects_s  = Vector{Sector{S, 2}}(undef, n_sects)
    sects_v = Vector{Sector{S, 2}}(undef, n_sects)

    blks_u  = Vector{Matrix{T}}(undef, n_sects)
    blks_s  = Vector{Vector{Float64}}(undef, n_sects)
    blks_v = Vector{Matrix{T}}(undef, n_sects)

    semits = collect(onlysemitokens(A.blocks))
    i=1
    fact = svd(zeros(T, 1,1), full=false)
    for (sect, blk) in A.blocks
        c1, c2 = sect
        sects_u[i] = Sector{S}(c1, c1)
        sects_s[i] = Sector{S}(c1, c2)
        sects_v[i] = Sector{S}(c2, c2)

        try
            fact = svd(blk, full=false)
        catch e
            @warn "error $e"
            fact = svd(blk, full=false, alg=LinearAlgebra.QRIteration())
        end
        blks_u[i] = fact.U
        blks_s[i] = fact.S
        blks_v[i] = fact.Vt
        i += 1
    end

    ns = truncatedsizes(blks_s, maxdim=maxdim, tol=tol)
    indices = findall(x->x>0, ns)

    sects_u  = sects_u[indices]
    sects_s  = sects_s[indices]
    sects_v = sects_v[indices]

    blks_u = [blks_u[index][:, 1:ns[index]] for index in indices]
    blks_s = [blks_s[index][1:ns[index]] for index in indices]
    blks_v = [blks_v[index][1:ns[index], :] for index in indices]

    middledims = ns[indices]

    n_sectors =length(middledims)
    Vl = VectorSpace{S}([sects_s[i][1]=>middledims[i] for i=1:n_sectors],
                        isdual(A.space[1]))
    Vr = VectorSpace{S}([sects_s[i][2]=>middledims[i] for i=1:n_sectors],
                        isdual(A.space[2]))

    u  = SymMatrix{S,T}(zero(S), (A.space[1], dual(Vl)),
                        SortedDict(zip(sects_u,blks_u)))
    s  = SymDiagonal{S,Float64}(A.charge, (Vl, Vr),
                                SortedDict([sects_s[i]=>Diagonal(blks_s[i]) for i in 1:n_sectors]))
    v = SymMatrix{S,T}(zero(S), (dual(Vr), A.space[2]),
                       SortedDict(zip(sects_v,blks_v)))

    u, s, v
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

    alpha < 1 && error("Renyi entropy only defined for positive alpha")
    any(P .< 0.0) && error("Negative value in P vector ", P)
    sum(P) ≈ 1 || println("P vector doesn't add up to one : $(sum(P))")
    Pnz = filter(x->x>0.0, P)
    if alpha == 1
        return -sum(Pnz .* log.(Pnz))
    else
        return log(sum(Pnz.^alpha))/(1-alpha)
    end
end

function _realwithcheck(a::Number, tol::Float64=1.e-12)
    imag(a) > tol && error("Not real! Im = $(imag(a))")
    real(a)
end
