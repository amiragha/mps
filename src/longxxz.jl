function xxz_longrange(lx::Int, delta::Float64=1.0, r::Float64=0.5)
    (lx < 2) && error("The system size is too small : ", lx)

    sz, sp, sm = sz_half, sp_half, sm_half

    heis_term = 0.5 * (sp ⊗ sm + sm ⊗ sp) + delta * sz ⊗ sz

    Hmat = heis_term
    # recursive generation of Hamiltonian
    if lx > 2
        for n=3:lx
            Hmat = (Hmat ⊗ I(2)) + (I(2^(n-2)) ⊗ heis_term)
            for l in n-2:-1:1
                d = n-l-1
                Hmat += (r^d) * I(2^(l-1)) ⊗ (
                    delta * (sz ⊗ I(2^d) ⊗ sz) +
                    0.5 * (sp ⊗ I(2^d) ⊗ sm + sm ⊗I(2^d) ⊗ sp)
                )
            end
        end
    end
    Hmat
end
