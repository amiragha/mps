struct Point1D
    x :: Float64
end

struct Point2D
    x :: Float64
    y :: Float64
end

struct Bond1D{T<:Number}
    one :: Int
    two :: Int
    offx :: Int
    t   :: T
end

struct Bond2D{T<:Number}
    one :: Int
    two :: Int
    offx :: Int
    offy :: Int
    t   :: T
end

struct UnitCell1D{T<:Number}
    n :: Int
    sites :: Vector{Point1D}
    vs :: Vector{Point1D}
    bonds :: Vector{Bond1D{T}}
    mus :: Vector{T}
end

struct UnitCell2D{T<:Number}
    n :: Int
    sites :: Vector{Point2D}
    vs :: Vector{Point2D}
    bonds :: Vector{Bond2D{T}}
    mus :: Vector{T}
end

function enlargeunitcell(uc::UnitCell2D{T}, lx::Int, ly::Int) where{T}

    sites = Point2D[]
    bonds = Point2D[]
    for (y, x) in Iterators.product(1:ly, 1:lx)
        for s in uc.sites
            push!(sites, (y-1)*uc.vs[2] + (x-1)*uc.vs[1] + s)
        end
        for b in uc.bonds
            index1 = l+b.one
            if 0 < y+b.offy <= ly && 0 < x+b.offx <= lx
                index2 = l+b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                push!(bonds, Bond2D(index1, index2))
            elseif 0 < y+b.offy <= ly
                # x wrap
            elseif  0 < x+b.offx <= lx
                # y wrap
            else
                # double wrap
            end
        end
    end
    UnitCell2D(
        uc.n*ly*lx,
        sites,
        [multiply((lx,ly), uc.vs[1]), multiply((lx,ly), uc.vs[2])],
        bonds
    )
end

function chainunitcell(t::T, mu::T=zero(T)) where {T<:Number}
    uc = UnitCell1D(
        1,
        [Point1D(0,0)],
        [Point1D(+1)],
        [Bond2D(1, 1, +1, t)],
        [mu])
    uc
end

function squareunitcell(tx::T, ty::T, mu::T=zero(T)) where {T<:Number}
    uc = UnitCell2D(
        1,
        [Point2D(0,0)],
        [Point2D(0,1), Point2D(1,0)],
        [Bond2D(1, 1, +1, 0, tx), Bond2D(1, 1, 0, +1, ty)],
        [mu])
    uc
end

function triangular_unitcell(t1::T, t2::T, mu::T=zero(T)) where {T<:Number}
    uc = UnitCell2D(
        1,
        [Point2D(0,0)],
        [Point2D(0.5,sin(pi/3)), Point2D(1,0)],
        [Bond2D(1, 1, +1, 0, t2), Bond2D(1, 1, 0, +1, t1), Bond2D(1, 1, -1, +1, t1)],
        [mu])
    uc
end

abstract type AbstractHamiltonian end
mutable struct BdGHamiltonian{T} <: AbstractHamiltonian
    nsites :: Int
    mat :: Matrix{T}
end
function BdGHamiltonian{T}(nsites::Int) where {T<:Number}
    BdGHamiltonian{T}(nsites, zeros(nsites, nsites))
end

function initializehamilt(::Type{T}, mode::Symbol, nsites::Int) where {T<:Number}
    if mode == :BdG
        BdGHamiltonian{T}(nsites)
    elseif mode == :SPIN
        SpinHamiltonian(nsites)
    elseif mode == :SPINU1
        U1SpinHamiltonian(nsites, div(nsites, 2))
    else
        error("unrecognized mode : $mode")
    end
end

function addterm!(H::BdGHamiltonian{T}, i::Int, mu::T) where {T<:Number}
    H.mat[i, i] += mu
    H
end

function addterm!(H::BdGHamiltonian{T}, i::Int, j::Int, t::T) where {T<:Number}
    H.mat[i, j] += t
    H.mat[j, i] += conj(t)
    H
end

mutable struct SpinHamiltonian <: AbstractHamiltonian
    nsites :: Int
    mat :: SparseMatrixCSC{Float64, Int}
end
function SpinHamiltonian(nsites::Int)
    SpinHamiltonian(nsites, spzeros(2^nsites, 2^nsites))
end

function addterm!(H::SpinHamiltonian, i::Int, mu::Float64)
    H.mat += mu * opextend(sz_half+1/2*I2, H.nsites, i)
    H
end

function addterm!(H::SpinHamiltonian, i::Int, j::Int, J::Float64)
    i == j && error("Add term for the same site!")
    idx1st = j > i ? i : j
    dis = abs(i - j)
    if dis > 1
        between = foldl(⊗, [2*sz_half for i=1:dis-1])
        H.mat += J * opextend(2*sp_half*sz_half ⊗ between ⊗ sm_half, H.nsites-dis, idx1st)
        H.mat += J * opextend(2*sz_half*sm_half ⊗ between ⊗ sp_half, H.nsites-dis, idx1st)
    else
        H.mat += J * opextend(2*sp_half*sz_half ⊗ sm_half, H.nsites-dis, idx1st)
        H.mat += J * opextend(2*sz_half*sm_half ⊗ sp_half, H.nsites-dis, idx1st)
    end
    H
end

mutable struct U1SpinHamiltonian <: AbstractHamiltonian
    nsites  :: Int
    mat     :: SparseMatrixCSC{Float64}
    I       :: Vector{Int32}
    J       :: Vector{Int32}
    V       :: Vector{Float64}
    sector  :: Int
    states  :: Vector{Int}
end

function U1SpinHamiltonian(nsites, sector)
    blocksize = binomial(nsites, sector)
    states = Vector{Int}(undef, blocksize)
    index = 0
    for state=0:2^nsites-1
        if count_ones(state) == sector
            index += 1
            states[index] = state
        end
    end
    U1SpinHamiltonian(nsites,
                      spzeros(blocksize, blocksize),
                      [], [], [],
                      sector, states)
end

function addterm!(H::U1SpinHamiltonian, i::Int, mu::Float64)
    for idx in eachindex(H.states)
        state = H.states[idx]
        if (state >> i) & 1 != 0
            push!(H.I, idx)
            push!(H.J, idx)
            push!(H.V, mu)
        end
    end
    H
end

function addterm!(H::U1SpinHamiltonian, i::Int, j::Int, J::Float64)
    i == j && error("Add term for the same site!")
    idx1st = j > i ? i : j
    dis = abs(i - j)
    for idx in eachindex(H.states)
        state = H.states[idx]
        if xor((state >> (i-1)) & 1, (state >> (j-1)) & 1) != 0
            if dis > 1
                sign = iseven(dis-1 - count_ones(state & (sum([2^n for n=0:dis-2])*2^idx1st))) ? +1 : -1
            else
                sign = +1
            end
            flippedstate = xor(state, 1 << (i-1) | 1 << (j-1))
            idx_flipped = searchsortedfirst(H.states, flippedstate)
            push!(H.I, idx)
            push!(H.J, idx_flipped)
            push!(H.V, -sign*J)
        end
    end
    H
end

function makehamiltonian(
    uc::UnitCell2D{T},
    lx::Int,
    ly::Int;
    boundary::Symbol=:OBC,
    mode::Symbol=:BdG) where {T<:Number}

    nsites = uc.n*ly*lx
    H = initializehamilt(T, mode, nsites)
    for (y,x) in Iterators.product(1:ly, 1:lx)
        l = (x-1) * ly * uc.n + (y-1) * uc.n
        for n in 1:uc.n
            addterm!(H, l+n, uc.mus[n])
        end
        for b in uc.bonds
            index1 = l+b.one
            if 0 < y+b.offy <= ly && 0 < x+b.offx <= lx
                index2 = l+b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                addterm!(H, index1, index2, -b.t)
            elseif boundary == :PBCX && 0 < y+b.offy <= ly
                x+b.offx > lx ? lnew = l - lx * ly * uc.n : lnew = l + lx * ly * uc.n
                index2 = lnew + b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                addterm!(H, index1, index2, -b.t)
            elseif boundary == :APBCX && 0 < y+b.offy <= ly
                x+b.offx > lx ? lnew = l - lx * ly * uc.n : lnew = l + lx * ly * uc.n
                index2 = lnew + b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                addterm!(H, index1, index2, b.t)
            elseif !(boundary in [:OBC, :PBCX, :APBCX])
                error("unrecognized boundary condition : ", boundary)
            end
        end
    end
    H
end

##TODO: combine this stuff into a function that works for every dimension.
function makemodel(
    uc::UnitCell2D{T},
    lx::Int,
    ly::Int;
    boundary::Symbol=:OBC) where{T<:Number}

    nsites = uc.n*ly*lx
    H = zeros(T, nsites, nsites)
    for (y,x) in Iterators.product(1:ly,1:lx)
        l = (x-1) * ly * uc.n + (y-1) * uc.n
        for n in 1:uc.n
            H[l+n,l+n] = uc.mus[n]
        end
        for b in uc.bonds
            index1 = l+b.one
            if 0 < y+b.offy <= ly && 0 < x+b.offx <= lx
                index2 = l+b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                #println("$x, $y => ($index1, $index2)")
                H[index1, index2] += -b.t
                H[index2, index1] += -conj(b.t)
            elseif boundary == :PBCX && 0 < y+b.offy <= ly
                x+b.offx > lx ? lnew = l - lx * ly * uc.n : lnew = l + lx * ly * uc.n
                index2 = lnew + b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                #println("wrap x, $y => ($index1, $index2)")
                H[index1, index2] += -b.t
                H[index2, index1] += -conj(b.t)
            elseif !(boundary in [:OBC, :PBCX])
                error("unrecognized boundary condition : ", boundary)
            end
        end
    end
    H
end

function makemodelJW(
    uc::UnitCell2D{T},
    lx::Int,
    ly::Int;
    boundary::Symbol=:OBC) where{T<:Number}

    @assert uc.n == 1
    nsites = uc.n*ly*lx
    H = spzeros(T, 2^nsites, 2^nsites)
    for (y,x) in Iterators.product(1:ly,1:lx)
        l = (x-1) * ly * uc.n + (y-1) * uc.n
        for n in 1:uc.n
            H += uc.mus[n] .* opextend(sz_half+1/2*I2, lx*ly*uc.n, l+n)
        end
        for b in uc.bonds
            index1 = l+b.one
            if 0 < y+b.offy <= ly && 0 < x+b.offx <= lx
                index2 = l+b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                println("$x, $y => ($index1, $index2)")
                if index2 > index1
                    idx1st = index1
                    dis = index2 - index1
                else
                    dis = index1 - index2
                    idx1st = index2
                end
                if dis > 1
                    between = foldl(⊗, [2*sz_half for i=1:dis-1])
                    H += -b.t * opextend(2*sp_half*sz_half ⊗ between ⊗ sm_half, nsites-dis, idx1st)
                    H += -conj(b.t) * opextend(2*sz_half*sm_half ⊗ between ⊗ sp_half, nsites-dis, idx1st)
                else
                    H += -b.t * opextend(2*sp_half*sz_half ⊗ sm_half, nsites-dis, idx1st)
                    H += -conj(b.t) * opextend(2*sz_half*sm_half ⊗ sp_half, nsites-dis, idx1st)
                end
            elseif boundary == :PBCX && 0 < y+b.offy <= ly
                x+b.offx > lx ? lnew = l - lx * ly * uc.n : lnew = l + lx * ly * uc.n
                index2 = lnew + b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                println("wrap x, $y => ($index1, $index2)")
                if index2 > index1
                    idx1st = index1
                    dis = index2 - index1
                else
                    dis = index1 - index2
                    idx1st = index2
                end
                if dis > 1
                    between = foldl(⊗, [2*sz_half for i=1:dis-1])
                    H += -b.t * opextend(2*sp_half*sz_half ⊗ between ⊗ sm_half, nsites-dis, idx1st)
                    H += -conj(b.t) * opextend(2*sz_half*sm_half ⊗ between ⊗ sp_half, nsites-dis, idx1st)
                else
                    H += -b.t * opextend(2*sp_half*sz_half ⊗ sm_half, nsites-dis, idx1st)
                    H += -conj(b.t) * opextend(2*sz_half*sm_half ⊗ sp_half, nsites-dis, idx1st)
                end
            elseif !(boundary in [:OBC, :PBCX])
                error("unrecognized boundary condition : ", boundary)
            end
        end
    end
    H
end

function makemodelJW_u1sym(
    uc::UnitCell2D{T},
    lx::Int,
    ly::Int;
    boundary::Symbol=:OBC,
    zsector::Int=((lx*ly*uc.n) % 2)) where{T<:Number}

    @assert uc.n == 1
    nsites = uc.n*ly*lx
    #H = spzeros(T, 2^nsites, 2^nsites)

    M = div(nsites+zsector, 2)
    block_size = binomial(nsites, M)

    #println("$M, $block_size")
    I = Int32[]
    J = Int32[]
    V = T[]

    szblock_states = Vector{Int}(undef, block_size)
    index = 0
    for state=0:2^(nsites)-1
        if count_ones(state) == M
            #println("$state, $(string(state, base=2))")
            index += 1
            szblock_states[index] = state
        end
    end

    # for all states
    for state_index=1:block_size
        state = szblock_states[state_index]
        # for all terms in Hamiltonian
        for (y,x) in Iterators.product(1:ly,1:lx)
            l = (x-1) * ly * uc.n + (y-1) * uc.n
            for n in 1:uc.n
                if (state >> (l+n)) & 1 != 0
                    push!(I, state_index)
                    push!(J, state_index)
                    push!(V, uc.mus[n])
                end
            end
            for b in uc.bonds
                index1 = l+b.one
                if 0 < y+b.offy <= ly && 0 < x+b.offx <= lx
                    index2 = l+b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                    #println("$x, $y => ($index1, $index2)")
                    if index2 > index1
                        idx1st = index1
                        dis = index2 - index1
                    else
                        dis = index1 - index2
                        idx1st = index2
                    end
                    if xor((state >> (index1-1)) & 1, (state >> (index2-1)) & 1) != 0
                        if dis > 1
                            #println("dis = $dis, state = $(string(state, base=2)), stuff = $(string(sum([2^i for i=1:dis-1])*2^idx1st, base=2)), and = $(state & (sum([2^i for i=1:dis-1])*2^idx1st))")
                            sign = iseven(dis-1 - count_ones(state & (sum([2^i for i=0:dis-2])*2^idx1st))) ? +1 : -1
                        else
                            sign = +1
                        end

                        flipped_state = xor(state, 1 << (index1-1) | 1 << (index2-1))
                        flipped_state_index = searchsortedfirst(szblock_states, flipped_state)

                        #println("$index1, $index2, $(string(state, base=2)), $(string(flipped_state, base=2)), indeces : $state_index, $flipped_state_index, sign = $sign")
                        if flipped_state_index > block_size
                            error("Fucked, $(string(state, base=2)), $(string(flipped_state, base=2))")
                        end
                        push!(I, state_index)
                        push!(J, flipped_state_index)
                        push!(V, sign*b.t)
                    end

                elseif boundary == :PBCX && 0 < y+b.offy <= ly
                    x+b.offx > lx ? lnew = l - lx * ly * uc.n : lnew = l + lx * ly * uc.n
                    index2 = lnew + b.two + (b.offx * uc.n * ly) + (b.offy * uc.n)
                    #println("wrap x, $y => ($index1, $index2)")
                    if index2 > index1
                        idx1st = index1
                        dis = index2 - index1
                    else
                        dis = index1 - index2
                        idx1st = index2
                    end
                    if xor((state >> (index1-1)) & 1, (state >> (index2-1)) & 1) != 0
                        if dis > 1
                            sign = iseven(dis-1-count_ones(state & (sum([2^i for i=0:dis-2])*2^idx1st))) ? +1 : -1
                        else
                            sign = +1
                        end
                        flipped_state = xor(state, 1 << (index1-1) | 1 << (index2-1))
                        flipped_state_index = searchsortedfirst(szblock_states, flipped_state)
                        push!(I, state_index)
                        push!(J, flipped_state_index)
                        push!(V, sign*b.t)
                    end
                elseif !(boundary in [:OBC, :PBCX])
                    error("unrecognized boundary condition : ", boundary)
                end
            end
        end
    end
    #println(I)
    #println(J)
    return sparse(I, J, V, block_size, block_size, +), (szblock_states .+ 1)
end
