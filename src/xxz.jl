function xxz_hamiltonian(lx::Int64, delta::Float64=1.0;
                         boundary::Symbol=:OBC, mode::Symbol=:FULL)
    @assert lx > 1
    if mode == :FULL
        return xxz_explicit(lx, delta, boundary)
    elseif mode == :U1
        return xxz_u1sym(lx, delta, boundary=boundary)
    else
        error("unrecognized mode :", mode)
    end
end

function xxz_explicit(lx::Int64, delta::Float64, boundary::Symbol)

    heis_term = 0.5 * (kron(sp_half, sm_half) + kron(sm_half, sp_half)) +
        delta * kron(sz_half,sz_half)

    Hmat = heis_term
    if lx > 2
        # recursive generation of open chain
        for i=3:lx
            Hmat = kron(Hmat, I(2)) + kron(I(2^(i-2)), heis_term)
        end
    end

    if boundary == :OBC
        return Hmat
    elseif boundary == :PBC
        Hmat = Hmat +
            0.5 * kron(kron(sp_half, I(2^(lx-2))), sm_half) +
            0.5 * kron(kron(sm_half, I(2^(lx-2))), sp_half) +
            delta * kron(kron(sz_half, I(2^(lx-2))), sz_half)
        return Hmat
    else
        error("unrecognized boundary condition :", boundary)
    end
end

"""
    xxz_u1sym(Lx, delta; boundary, zsector)
"""

# # Generates the XXZ Hamiltonian of size `lx` in the given `zsector`
# # block. This is possible because ``\left[ S^z_{\text{tot}},
# # H\right]=0`` and its eigenvalues can be used to label separate
# # magnetization blocks in the Hamiltonian matrix and each block can be
# # diagonalized separately to significantly reduce time complexity. Each
# # `zsector` block consists of Ising vectors with exactly `M = (zsector
# # + lx)/2` spins up, so the size of the block is ``\binom{lx}{M}``. Using
# # Sterling formula the size of the largest block, ``M=lx/2``, is ``√{π
# # lx/2}`` times smaller than the full Hilbert space, which For typical
# # sizes of ED calculation ,``lx∼ 30``, is roughly about ``∼\!7`` times.

# """
function xxz_u1sym(lx::Int, delta::Float64=1.0;
                   boundary::Symbol=:OBC, zsector::Int=(lx % 2))

    @assert lx > 1
    @assert (lx-zsector) % 2 == 0

    if boundary == :OBC
        bond_range = 0:lx-2
    elseif boundary == :PBC
        bond_range = 0:lx-1
    else
        error("unrecognized boundary condition :", boundary)
    end

    M = div(lx + zsector, 2)
    block_size = binomial(lx, M)

    I = Int32[]
    J = Int32[]
    V = Float64[]

    szblock_states = Vector{Int}(undef, block_size)
    index = 0
    for state=0:2^lx-1
        if count_ones(state) == M
            index += 1
            szblock_states[index] = state
        end
    end

    # For all states (of magnetization M):
    for state_index=1:block_size
        state = szblock_states[state_index]
        # For all bond terms of the Hamiltonian
        for site=bond_range
            site_nn = (site + 1) % lx
            # if spin are parallel:
            if xor((state >> site) & 1, (state >> site_nn) & 1) == 0
                push!(I, state_index)
                push!(J, state_index)
                push!(V, delta * 0.25)
                # if spin are not parallel:
            else
                push!(I, state_index)
                push!(J, state_index)
                push!(V, - delta * 0.25)

                flipped_state = xor(state, 1 << site | 1 << site_nn)
                flipped_state_index = searchsortedfirst(szblock_states, flipped_state)
                push!(I, state_index)
                push!(J, flipped_state_index)
                push!(V, 0.5)
            end
        end
    end

    return sparse(I, J, V, block_size, block_size, +), (szblock_states .+ 1)
 end

"""
    j1j2_explicit(Lx, j1, j2)

explicit construction of j1j2 sparse full Hamiltonian.

"""
function j1j2_explicit(Lx::Int64,
                       j1::Float64,
                       j2::Float64;
                       boundary::Symbol=:open)
    @assert Lx > 2
    Sz = sparse(sz_half)
    Sp = sparse(sp_half)
    Sm = sparse(sm_half)
    I2 = I(2)

    heis_term1 = j1 * 0.5 * (kron(Sp, Sm) + kron(Sm, Sp)) + j1 * kron(Sz,Sz)
    heis_term2 = j2 * 0.5 * (kron(kron(Sp, I2), Sm) + kron(kron(Sm, I2), Sp)) +
        j2 * kron(kron(Sz, I2), Sz)

    Hmat = heis_term1
    for i=3:Lx
        Hmat = kron(Hmat, I2) + kron(I(2^(i-2)), heis_term1) +
            kron(I(2^(i-3)), heis_term2)
    end

    if boundary == :open
        return Hmat
    elseif boundary == :periodic
        Hmat = Hmat +
            j1 * 0.5 * kron(kron(Sp,I(2^(Lx-2))), Sm) +
            j1 * 0.5 * kron(kron(Sm,I(2^(Lx-2))), Sp) +
            j1 * kron(kron(Sz,I(2^(Lx-2))), Sz) +
            j2 * 0.5 * kron(I2, kron(kron(Sp,I(2^(Lx-3))), Sm)) +
            j2 * 0.5 * kron(I2, kron(kron(Sm,I(2^(Lx-3))), Sp)) +
            j2 * kron(I2, kron(kron(Sz,I(2^(Lx-3))), Sz)) +
            j2 * 0.5 * kron(kron(kron(Sp,I(2^(Lx-3))), Sm), I2) +
            j2 * 0.5 * kron(kron(kron(Sm,I(2^(Lx-3))), Sp), I2) +
            j2 * kron(kron(kron(Sz,I(2^(Lx-3))), Sz), I2)
        return Hmat
    else
        error("unrecognized boundary conditions :", boundary)
    end
end
