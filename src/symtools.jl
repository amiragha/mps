function svd(A::AbstractSymMatrix)
    T = eltype(A)
    S = vtype(A)

    u_blocks = SortedDict{Sector{S, 2}, Matrix{T}}()
    s_blocks = SortedDict{Sector{S, 2}, Diagonal{Float64, Vector{Float64}}}()
    v_blocks = SortedDict{Sector{S, 2}, Matrix{T}}()

    midl = SortedDict{S, Int}()
    midr = SortedDict{S, Int}()
    for (sect,blk) in A.blocks
        c1, c2 = sect
        fact = svd(blk, full=false)
        u_blocks[Sector(c1, c1)] = fact.U
        s_blocks[Sector(c1, c2)] = Diagonal(fact.S)
        v_blocks[Sector(c2, c2)] = fact.Vt
        midl[c1] = length(fact.S)
        midr[c2] = length(fact.S)
    end

    Vl = VectorSpace{S}(midl, A.space[1].isdual)
    Vr = VectorSpace{S}(midr, A.space[2].isdual)

    # println(A.space[1])
    # println(dual(Vl))
    # println(Vl)
    # println(Vr)
    # println(dual(Vr))
    # println(A.space[2])
    u  = SymMatrix{S,T}(zero(S), (A.space[1], dual(Vl)), u_blocks)
    s  = SymDiagonal{S,Float64}(A.charge, (Vl, Vr), s_blocks)
    v = SymMatrix{S,T}(zero(S), (dual(Vr), A.space[2]), v_blocks)

    u, s, v
end


"""
    fermionswapgate(l1, l2)

Fot the two given legs, `l1` and `l2` make the fermionic swap gate
which is `-1` if the two legs has charge odd and `+1` otherwise.

The output legs are l2_afterX, l1_afterX, l2, l1
"""
function fermionswapgate(l1::VectorSpace{S}, l2::VectorSpace{S}) where{S}

    space = (l2, l1, dual(l2), dual(l1))

    sects, sizes = _allsectorsandsizes(0, space)
    blocks = SortedDict{Sector{S, 4}, Array{Float64, 4}}()

    for index in eachindex(sects)
        c2, c1, c2X, c1X = sects[index]
        if c1 == c1X && c2 == c2X
            d2, d1 = sizes[index][1:2]
            s = isodd(c1) && isodd(c2) ? -1 : +1
            blocks[sects] =  s * reshape(I(d1*d2), d2, d1, d2, d1)
        else
            blocks[sects] = zeros(sizes[index])
        end
    end
    SymTensor(zero(s), space, blocks)
end

# function isrightisometry(A::SymTensor)
#     N = rank(A)

# end

function isrightisometry(A::SymMatrix)
    eye(eltype(A), A.space[1]) ≈ SymMatrix(contract(A, (1, -1), dual(A), (2, -1)))
end

function isleftisometry(A::SymMatrix)
    eye(eltype(A), A.space[2]) ≈ SymMatrix(contract(dual(A), (-1, 1), A, (-1, 2)))
end

isunitary(A::SymMatrix) = isrightisometry(A) && isleftisometry(A)
