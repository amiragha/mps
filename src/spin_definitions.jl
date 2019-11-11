eye(n) = sparse(1.0I, n, n)

"""
    spinoperators(S)

For spin `S` returns the operators (`Sz`,`Sp`, `Sm`)
"""

function spinoperators(S)
    d = Int(round(2*S + 1))
    ms = [S-n+1 for n in 1:d]
    amps = [sqrt((S+m)*(S-m+1)) for m in ms[1:d-1]]

    Sz = diagm(0 => ms)
    Sp = diagm(1 => amps)
    Sm = diagm(-1 => amps)

    Sz, Sp, Sm
end

# generic spin_half
const sz_half = sparse(Float64[.5 0; 0 -.5])
const sp_half = sparse(Float64[ 0 1; 0  0])
const sm_half = sparse(Float64[ 0 0; 1  0])
const I2 = eye(2)

# U1 symmetric spin_half
const sz_half_U1sym = SymTensor(0,
                                (STLeg(+1, [0, 1], [1, 1]),
                                 STLeg(-1, [0, 1], [1, 1])),
                                [(0, 0), (1, 1)],
                                [-.5*ones(1,1), .5*ones(1,1)])

const sp_half_U1sym = SymTensor(+1,
                                (STLeg(+1, [0, 1], [1, 1]),
                                 STLeg(-1, [0, 1], [1, 1])),
                                [(1, 0)],
                                [ones(1,1)])

const sm_half_U1sym = SymTensor(-1,
                                (STLeg(+1, [0, 1], [1, 1]),
                                 STLeg(-1, [0, 1], [1, 1])),
                                [(0, 1)],
                                [ones(1,1)])


"""
    ringexchangeoperator(n)

defined based on its action on Ising basis, which is a cricular shift
to right plus the hermitian conjugate.

With this definition heiseberg is 1/4 *(R(2)-I(4))
"""
function ringexchangeoperator(n::Int)
    @assert 1 < n < 7
    Pn = spzeros(2^n, 2^n)
    for j=1:2^n
        i = ((j-1) >> 1) | (((j-1) & 1) << (n-1))
        Pn[i+1, j] = 1
    end
    Pn + Pn'
end

# only for spin-half right now
function nbodyopexpansion(n::Int,
                          H::AbstractMatrix;
                          verbose::Bool=false)
    @assert 1 < n < 7
    @assert size(H) == (2^n, 2^n)
    ⊗ = kron
    Z = sz_half
    P = sp_half
    M = sm_half
    E = I(2)
    ops = [Z, P, M, E]
    ols = ['Z', 'P', 'M', 'E']

    T = eltype(H)
    amps = Vector{T}()
    terms = Vector{NTuple{n, Matrix{T}}}()
    for is in Iterators.product([1:4 for i=1:n]...)
        op = reduce(⊗, [ops[i] for i in is], init=I(1))
        opnorm = norm(op)
        amp = dot(H, op)/norm(op)^2
        if amp != 0
            term = [Matrix(ops[i]) for i in is]
            term[n] = term[n] * amp
            push!(terms, Tuple(term))
            verbose && println([ols[i] for i in is]..., "  $amp")
            is == Tuple([4 for i=1:n]) &&
                @warn "Operator has identity with amplitude $amp"
        end
    end
    terms
end
