mutable struct KetState{T<:RLorCX}
    lx :: Int
    d  :: Int
    v  :: Vector{T}
end

function randket(T::Type, lx::Int, d::Int; rng::AbstractRNG=GLOBAL_RNG)
    v = randn(rng, T, d^lx)
    KetState{T}(lx, d, v)
end

randket(lx::Int, d::Int; rng::AbstractRNG=GLOBAL_RNG) = randket(Float64, lx, d; rng=rng)

function normalize!(v::KetState{T}) where {T<:RLorCX}
    v.v = v.v ./ norm(v.v)
    nothing
end

@inline opextend(op, lx::Int, l::Int) = eye(2^(l-1)) ⊗ op ⊗ eye(2^(v.lx-l))

function apply1site(v::KetState{T}, l::Int, op::SparseMatrixCSC{T,Int}) where {T<:RLorCX}
    size(op) != (v.d, v.d) && error("Operator size not correct!")
    (0 < l <= v.lx) && error("position l is not within the lattice")

    KetState(v.lx, v.d, opext(op, v.lx, l) * v.v)
end

function measure1point(v::KetState{T}, op::SparseMatrixCSC{T,Int}) where {T<:RLorCX}
    result = Vector{T}(undef, lx)
    for l in 1:v.lx
        result[l] = dot(v.v, opextend(op, v.lx, l) * v.v)
    end
    result
end
