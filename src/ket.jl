mutable struct KetState{T<:RLorCX}
    lx :: Int
    d  :: Int
    v  :: Vector{T}
end

function randket(T::Type, rng::AbstractRNG, lx::Int, d::Int)# where {T<:RLorCX}
    v = randn(rng, T, d^lx)
    KetState{T}(lx, d, v)
end

randket(rng::AbstractRNG, lx::Int, d::Int) = randket(Float64, rng, lx, d)

function normalize!(v::KetState{T}) where {T<:RLorCX}
    v.v = v.v ./ norm(v.v)
    nothing
end

@inline opextend(op, lx::Int, l::Int) = eye(2^(l-1)) ⊗ op ⊗ eye(2^(lx-l))

function apply1site(v::KetState{T}, l::Int, op::SparseMatrixCSC{T,Int}) where {T<:RLorCX}
    size(op) != (v.d, v.d) && error("Operator size not correct!")
    (0 < l <= v.lx) || error("position l is not within the lattice")

    KetState(v.lx, v.d, opextend(op, v.lx, l) * v.v)
end

apply1site(v::KetState{ComplexF64}, l::Int, op::SparseMatrixCSC{Float64, Int}) =
    apply1site(v,l, convert(SparseMatrixCSC{ComplexF64,Int}, op))

function measure1point(v::KetState{T}, op::SparseMatrixCSC{T,Int}) where {T<:RLorCX}
    result = Vector{T}(undef, v.lx)
    for l in 1:v.lx
        result[l] = dot(v.v, opextend(op, v.lx, l) * v.v)
    end
    result
end

measure1point(v::KetState{ComplexF64}, op::SparseMatrixCSC{Float64, Int}) =
    measure1point(v, convert(SparseMatrixCSC{ComplexF64,Int}, op))
