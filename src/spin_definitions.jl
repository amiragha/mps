"""
    spinoperators(S)

For spin `S` returns the operators (`Sz`,`Sp`, `Sm`). `S` has to be
integer or half-integer
"""
function spinoperators(S; symmetry::Symbol=:NONE)
    d = Int(round(2*S + 1))
    ms = [S-n+1 for n in 1:d]
    amps = [sqrt((S+m)*(S-m+1)) for m in ms[1:d-1]]

    if symmetry == :NONE
        Sz = diagm(0 => ms)
        Sp = diagm(1 => amps)
        Sm = diagm(-1 => amps)
        return Sz, Sp, Sm
    elseif symmetry == :U1
        chrs = collect(0:d-1)
        dims = ones(Int, d)
        legs = (STLeg(+1, chrs, dims), STLeg(-1, chrs, dims))
        Sz = SymMatrix(0, legs, [(c,c) for c in 0:d-1],
                       [m*ones(1,1) for m in reverse(ms)])
        Sp = SymMatrix(+1, legs, [(c+1, c) for c in 0:d-2],
                       [m * ones(1,1) for m in amps])
        Sm = SymMatrix(-1, legs, [(c, c+1) for c in 0:d-2],
                       [m * ones(1,1) for m in amps])
        return Sz, Sp, Sm
    else
        error("Unkown symmetry: $symmetry")
    end
end

# generic spin_half
const sz_half = sparse(Float64[.5 0; 0 -.5])
const sp_half = sparse(Float64[ 0 1; 0  0])
const sm_half = sparse(Float64[ 0 0; 1  0])
const I2 = I(2)

# U1 symmetric spin_half
const sz_half_U1sym = SymMatrix(0,
                                (STLeg(+1, [0, 1], [1, 1]),
                                 STLeg(-1, [0, 1], [1, 1])),
                                [(0, 0), (1, 1)],
                                [-.5*ones(1,1), .5*ones(1,1)])

const sp_half_U1sym = SymMatrix(+1,
                                (STLeg(+1, [0, 1], [1, 1]),
                                 STLeg(-1, [0, 1], [1, 1])),
                                [(1, 0)],
                                [ones(1,1)])

const sm_half_U1sym = SymMatrix(-1,
                                (STLeg(+1, [0, 1], [1, 1]),
                                 STLeg(-1, [0, 1], [1, 1])),
                                [(0, 1)],
                                [ones(1,1)])

######################
### spin-1/2 stuff ###

"""
    permutespins(op, perm)
Given an operator `op` on spins 1...n return the same operator on the
given permutation `perm` of those spins.

"""
function permutespins(op::Matrix, perm::Vector{Int})
    size(op, 1) == size(op, 2) || error("square matrix please!")
    n = Int(log2(size(op,1)))
    sort(perm) == collect(1:n) || error("not a valid permutation!")

    mask = 2^n - 1
    for num=0:mask
        0 # to do using some bindary represenation of integers
    end
    op
end

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
                          symmetry::Symbol=:NONE,
                          verbose::Bool=false)
    1 < n < 7 || error("too large!")
    size(H) == (2^n, 2^n) || error("operator size doesn't match n")

    ops = [spinoperators(1/2)..., I(2)]
    ols = ['Z', 'P', 'M', 'I']

    T = eltype(H)
    symops = [spinoperators(1/2, symmetry=:U1)..., eye(T, [0,1], [1,1])]

    amps = Vector{T}()
    if symmetry == :NONE
        terms = Vector{NTuple{n, Matrix{T}}}()
    elseif symmetry == :U1
        terms = Vector{NTuple{n, SymMatrix{T}}}()
    end
    ⊗ = kron
    for is in Iterators.product([1:4 for i=1:n]...)
        op = reduce(⊗, [ops[i] for i in is], init=I(1))
        opnorm = norm(op)
        amp = dot(H, op)/norm(op)^2
        if amp != 0
            verbose && println([ols[i] for i in is]..., "  $amp")
            is == Tuple([4 for i=1:n]) &&
                @warn "Operator has identity with amplitude $amp"
            if symmetry == :NONE
                term = [Matrix(ops[i]) for i in is]
                term[n] = term[n] * amp
                push!(terms, Tuple(term))
            elseif symmetry == :U1
                term = [symops[i] for i in is]
                term[n] = term[n] * amp
                push!(terms, Tuple(term))
            else
                error("symmetry not yet implemented: $symmetry")
            end
        end
    end
    terms
end
