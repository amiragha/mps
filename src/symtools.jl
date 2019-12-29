function svd(A::AbstractSymTensor)
    rank(A) == 2 ||
        error("svdsym only defined for matrix like objects rank is $(rank(A))")
    T = eltype(A)
    S = vtype(A)

    u_blocks = SortedDict{Sector{S, 2}, Matrix{T}}()
    s_blocks = SortedDict{Sector{S, 2}, Diagonal{T}}()
    v_blocks = SortedDict{Sector{S, 2}, Matrix{T}}()

    mid = SortedDict{S, Int}()
    for (sect,blk) in A.blocks
        c1, c2 = sect
        fact = svd(blk, full=false)
        u_blocks[Sector(c1, -c1)] = fact.U
        s_blocks[Sector(c1,  c2)] = Diagonal(fact.S)
        v_blocks[Sector(-c2, c2)] = fact.Vt
        mid[c2] = length(fact.S)
    end

    Vr = VectorSpace{S}(mid)
    Vl = mapcharges(x->x+A.charge, dual(Vr))

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
function fermionswapgate(l1::U1Space, l2::U1Space)
    (l1.chrs != [0, 1] || l1.dims != [1,1]) && error("fermionswap oops!")

    #l1X = U1Space(-l1.sign, l1.chrs, l1.dims)
    l1X = legdual(l1)
    #l2X = U1Space(-l2.sign, l2.chrs, l2.dims)
    l2X = legdual(l2)
    legs = (l2, l1, l2X, l1X)

    sects, sizes = _allsectorsandsizes(0, legs)
    nzblks = Vector{Array{Float64, 4}}(undef, length(sizes))

    for index in eachindex(sects)
        c2, c1, c2X, c1X = sects[index]
        if c1 == c1X && c2 == c2X
            d2, d1 = sizes[index][1:2]
            s = isodd(c1) && isodd(c2) ? -1 : +1
            nzblks[index] =  s * reshape(I(d1*d2), d2, d1, d2, d1)
        else
            nzblks[index] = zeros(sizes[index])
        end
    end
    SymTensor(0, legs, sects, nzblks)
end

# function isrightisometry(A::SymTensor)
#     N = rank(A)

# end

function isrightisometry(A::SymMatrix)
    leg = A.legs[1]
    eye(eltype(A), leg.chrs, leg.dims) ≈ SymMatrix(contract(A, (1, -1), tensordual(A), (2, -1)))
end

function isleftisometry(A::SymMatrix)
    leg = A.legs[2]
    eye(eltype(A), leg.chrs, leg.dims) ≈ SymMatrix(contract(tensordual(A), (-1, 1), A, (-1, 2)))
end

isunitary(A::SymMatrix) = isrightisometry(A) && isleftisometry(A)
