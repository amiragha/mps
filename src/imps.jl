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

    InfiniteMatrixProductState{T}(d, (Γ,), fact.S/norm(fact.S))
end

"""
    poshermitsqrt(A)

return the square root of the positvie hermitian matrix `A` using
eigenvalue decomposition.

"""
function poshermitsqrt(A::Matrix{T}) where{T<:Number}
    !isapprox(A, A') && error("Matrix is not hermitian! ", norm(A-A'))
    S, U = eigen(Hermitian(A))
    if !all(S .> 0)
        !all(S .< 0) && error("Matrix is not positive!", S)
        S = -S
    end

    U * Diagonal(sqrt.(S)) * U'
end

function _transfermatvec(v::Matrix{T}, uc::Array{T, 4}, direction::Symbol) where{T<:Number}
    if direction == :LEFT
        @tensor w[u,d] := (v[u',d'] * uc[u',o1,o2,u]) * conj(uc)[d',o1,o2,d]
        w
    elseif direction == :RIGHT
        @tensor w[u,d] := (v[u',d'] * uc[u,o1,o2,u']) * conj(uc)[d,o1,o2,d']
        w
    else
        error("unknown direction!")
    end
end

function InfiniteMatrixProductState(A::Array{T, 3},
                                    Λ::Vector{Float64},
                                    B::Array{T, 3},
                                    Λex::Vector{Float64}) where {T<:RLorCX}
    n = size(A, 1)
    d = size(A, 2)
    @assert size(A) == (n, d, n)
    @assert length(Λ) == n
    @assert size(B) == (n, d, n)
    @assert length(Λex) == n

    # @tensor LB[l,o,r] := Diagonal(Λ)[l,m] * B[m,o,r]
    # Q, Λr = qr(reshape(LB, n*d, n))
    # Q = Matrix(Q)
    @tensor Luc[l,o1,o2,r] := (A[l,o1,ml] * Diagonal(Λ)[ml,mm]) *
        B[mm,o2,mr] * Diagonal(1 ./Λex)[mr,r]

    es, vs, info = eigsolve(v->_transfermatvec(v, Luc, :LEFT),
                            rand(n,n), 1, :LR)
    if !(es[1] ≈ 1.0)
       println("The largest left eigenvalue of transfer matrix is not unit : ", es[1])
    end
    v = vs[1]
    if !(v ≈ conj(v))
        println("The left eigenvector is not real!", norm(v-conj(v)))
        v = real.(v)
    end
    X = poshermitsqrt(real.(v))

    # @tensor AL[l,o,r] := A[l,o,m] * Diagonal(Λ)[m,r]
    # Q, Λl = qr(transpose(reshape(AL, n, d*n)))
    # Q = transpose(Matrix(Q))
    # Λl = transpose(ΛL)
    @tensor Ruc[l,o1,o2,r] := (Diagonal(1 ./Λex)[l,ml] * A[ml,o1,mm]) *
        Diagonal(Λ)[mm, mr] * B[mr,o2,r]

    es, vs, info = eigsolve(v->_transfermatvec(v, Ruc, :RIGHT),
                            rand(n,n), 1, :LR)
    if !(es[1] ≈ 1.0)
        println("The largest right eigenvalue of transfer matrix is not unit : ", es[1])
    end
    v = vs[1]
    if !(v ≈ conj(v))
        println("The right eigenvector is not real!", norm(v-conj(v)))
    end
    Y = poshermitsqrt(real.(v))

    ###NOTE:
    ## Left normalzied canonical unitcell is (XA) Λ (BY) inv(Y Λex X)
    ## Right noralized canonical unitcell is inv(Y Λex X) (XA) Λ (BY)
    ## so naming A <- XA and B <-BY and Λex <- Y Λex X works!

    #println(X - conj(X))
    @tensor XA[l,o,r] := X[l,m] * A[m,o,r]
    Q, R = qr(reshape(XA, n*d, n))
    Γ1 = reshape(Matrix(Q), n,d,n)

    @tensor BY[l,o,r] := (R[l,ml] * Diagonal(Λ)[ml,mm]) * B[mm,o,mr] * Y[mr,r]
    Q, R = qr(reshape(BY, n*d, n))
    Γ2 = reshape(Matrix(Q), n,d,n)

    @show norm(R - Matrix(I, n, n))

    L = Y*Diagonal(Λex)*X
    norm(L - Diagonal(diag(L)))
    InfiniteMatrixProductState{T}(d, Γ1, inv(Y*Diagonal(Λex)*X))
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

    Γ = imps.Γs#[1]
    @tensor sbb[l,d1,d2,r] :=
        Diagonal(imps.Λ)[l, ml] * Γ[ml, d1, mr] * Γ[mr, d2, r]
    @tensor v = scalar(sbb[l,d1',d2',r] * bond[d1,d2,d1',d2'] * conj(sbb)[l,d1,d2,r])
    real(v)
end
