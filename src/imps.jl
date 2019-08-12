mutable struct InfiniteMatrixProductState{T<:RLorCX}
    d  :: Int
    Γ  :: Array{T, 3}
    Λ  :: Vector{Float64}
end

function InfiniteMatrixProductState(M::Array{T, 3}) where{T<:RLorCX}
    diml, d, dimr = size(M)
    @assert (diml == dimr)
    fact = svd(reshape(M, diml, d*dimr), full=false)
    Vt = reshape(fact.Vt, diml, d, dimr)

    @tensor Γ[l,d,r] := Vt[l,d,m] * fact.U[m,r]

    InfiniteMatrixProductState(d, Γ, fact.S/norm(fact.S))
end

function ketstate2imps(vstate::Vector{T}, lx::Int, d::Int;
                       maxdim::Int=d^div(lx, 2)) where{T<:RLorCX}
    @assert lx % 2 == 1
    @assert length(vstate) == d^lx

    l = div(lx, 2)
    fact = svd(reshape(vstate, d^(l+1), d^l), full=false)

    S, n, ratio = truncate(fact.S, maxdim=maxdim)
    Vt = fact.Vt[1:n, :]

    ##NOTE: Using inversion symmetry the left and right l sites have
    ##the same vector, so we only need to put an identity in the
    ##middle, conjugate and then contract to get rid of all the left
    ##and right tensors except the one in the middle (which we are
    ##looking for)

    psi_b = reshape(Vt, n, 1, d*ones(Int, l)...)
    psi_a = permutedims(psi_b, [1,2,(l+2:-1:3)...])

    ten_b = reshape(conj(psi_b), n, 1, d^l)
    ten_a = reshape(conj(psi_a), n, 1, d^l)
    I_ext = reshape(Matrix(1.0I, d, d), 1, d, 1, d)
    vstate_aib = reshape(vstate, d^l, 2, d^l)

    @tensor M[a,d,b] := ((vstate_aib[a',d', b'] * ten_a[a, il, a']) *
                         ten_b[b,ir, b']) * I_ext[il, d, ir, d']

    InfiniteMatrixProductState(M)
end

function measure_bond(imps::InfiniteMatrixProductState{T},
                      bond::Array{T, 4}) where{T<:Number}
    d = imps.d
    @assert size(bond) == (d,d,d,d)

    @tensor sbb[l,d1,d2,r] :=
        Diagonal(imps.Λ)[l, ml] * imps.Γ[ml, d1, mr] * imps.Γ[mr, d2, r]
    @tensor v = scalar(sbb[l,d1',d2',r] * bond[d1,d2,d1',d2'] * conj(sbb)[l,d1,d2,r])
    real(v)
end
