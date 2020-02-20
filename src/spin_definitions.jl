"""
    spinoperators(s)

For spin `s` returns the operators (`sz`,`sp`, `sm`). `s` has to be
integer or half-integer
"""
function spinoperators(s;
                       symmetry::Type{<:AbstractCharge}=Trivial)
    d = Int(round(2*s + 1))
    ms = [s-n+1 for n in 1:d]
    amps = [sqrt((s+m)*(s-m+1)) for m in ms[1:d-1]]

    if symmetry == Trivial
        sz = diagm(0 => ms)
        sp = diagm(1 => amps)
        sm = diagm(-1 => amps)
        return sz, sp, sm
    elseif symmetry == U1
        V = U1Space(c=>1 for c in 0:d-1)
        sz_blocks = SortedDict{Sector{U1, 2}, Matrix{Float64}}()
        sp_blocks = SortedDict{Sector{U1, 2}, Matrix{Float64}}()
        sm_blocks = SortedDict{Sector{U1, 2}, Matrix{Float64}}()
        for c in 0:d-1
            sz_blocks[Sector{U1}(c, c)] = reverse(ms)[c+1] * ones(1,1)
        end
        for c in 0:d-2
            sp_blocks[Sector{U1}(c+1, c)] = amps[c+1] * ones(1,1)
            sm_blocks[Sector{U1}(c, c+1)] = amps[c+1] * ones(1,1)
        end
        sz = SymMatrix(zero(U1), (V,dual(V)), sz_blocks)
        sp = SymMatrix(U1(+1), (V,dual(V)), sp_blocks)
        sm = SymMatrix(U1(-1), (V,dual(V)), sm_blocks)
        return sz, sp, sm
    else
        throw(ArgumentError("Unkown symmetry : $symmetry"))
    end
end

######################
### spin-1/2 stuff ###

"""
    permutespins(op, perm)
Given an operator `op` on spins 1...n return the same operator on the
given permutation `perm` of those spins.

"""
function permutespins(op::AbstractMatrix, perm::Vector{Int})
    size(op, 1) == size(op, 2) || error("square matrix please!")
    n = Int(log2(size(op,1)))
    sort(perm) == collect(1:n) || error("not a valid permutation!")

    mask = 2^n - 1
    indexes = Vector{Int}(undef, 2^n)
    for num=0:mask
        new = 0
        for i=1:n
            ni = perm[i]
            new += (num & (1 << (ni-1))) >> (ni - i)
        end
        indexes[num+1] = new+1
    end
    op[indexes, indexes]
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
                          mode::Symbol=:SYMBOL,
                          verbose::Bool=false)
    1 < n < 7 || error("too large!")
    size(H) == (2^n, 2^n) || error("operator size doesn't match n")

    ops = [spinoperators(1/2)..., I(2)]
    ols = ['Z', 'P', 'M', 'I']
    symbols = [:Z, :P, :M, :I]
    T = eltype(H)
    symops = [spinoperators(1/2, symmetry=U1)..., eye(T, U1Space(0=>1, 1=>1))]

    amps = Vector{T}()
    if mode == :SYMBOL
        MType = Symbol
    elseif mode == :Trivial
        MType = Matrix{T}
    elseif mode == :U1
        MType = SymTensor{U1,T,2}
    else
        error("unrecognized mode!")
    end
    terms = Vector{QAmpTerm{MType, n, T}}()
    ⊗ = kron
    for is in Iterators.product([1:4 for i=1:n]...)
        op = reduce(⊗, [ops[i] for i in is], init=I(1))
        opnorm = norm(op)
        amp = dot(H, op)/norm(op)^2
        if amp != 0
            verbose && println([ols[i] for i in is]..., "  $amp")
            is == Tuple([4 for i=1:n]) &&
                @warn "Operator has identity with amplitude $amp"
            if mode == :SYMBOL
                term = [symbols[i] for i in is]
            elseif mode == :Trivial
                term = [Matrix{T}(ops[i]) for i in is]
            elseif mode == :U1
                term = [symops[i] for i in is]
            else
                error("unrecognized mode : $mode")
            end
            push!(terms, QAmpTerm(amp, Tuple(term)))
        end
    end
    terms
end
