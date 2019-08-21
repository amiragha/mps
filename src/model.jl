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
    I2 = speye(Float64,2)

    heis_term1 = j1 * 0.5 * (kron(Sp, Sm) + kron(Sm, Sp)) + j1 * kron(Sz,Sz)
    heis_term2 = j2 * 0.5 * (kron(kron(Sp, I2), Sm) + kron(kron(Sm, I2), Sp)) +
        j2 * kron(kron(Sz, I2), Sz)

    Hmat = heis_term1
    for i=3:Lx
        Hmat = kron(Hmat, I2) + kron(eye(2^(i-2)), heis_term1) +
            kron(eye(2^(i-3)), heis_term2)
    end

    if boundary == :open
        return Hmat
    elseif boundary == :periodic
        Hmat = Hmat +
            j1 * 0.5 * kron(kron(Sp,eye(2^(Lx-2))), Sm) +
            j1 * 0.5 * kron(kron(Sm,eye(2^(Lx-2))), Sp) +
            j1 * kron(kron(Sz,eye(2^(Lx-2))), Sz) +
            j2 * 0.5 * kron(I2, kron(kron(Sp,eye(2^(Lx-3))), Sm)) +
            j2 * 0.5 * kron(I2, kron(kron(Sm,eye(2^(Lx-3))), Sp)) +
            j2 * kron(I2, kron(kron(Sz,eye(2^(Lx-3))), Sz)) +
            j2 * 0.5 * kron(kron(kron(Sp,eye(2^(Lx-3))), Sm), I2) +
            j2 * 0.5 * kron(kron(kron(Sm,eye(2^(Lx-3))), Sp), I2) +
            j2 * kron(kron(kron(Sz,eye(2^(Lx-3))), Sz), I2)
        return Hmat
    else
        error("unrecognized boundary conditions :", boundary)
    end
end
