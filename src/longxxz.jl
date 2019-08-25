function xxz_longrange(lx::Int, delta::Float64=1.0, r::Float64=0.5)
    (lx < 2) && error("The system size is too small : ", lx)

    sz, sp, sm = sz_half, sp_half, sm_half

    heis_term = 0.5 * (sp ⊗ sm + sm ⊗ sp) + delta * sz ⊗ sz

    Hmat = heis_term
    # recursive generation of Hamiltonian
    if lx > 2
        for n=3:lx
            Hmat = (Hmat ⊗ I2) + (eye(2^(n-2)) ⊗ heis_term)
            for l = n-2:1
                d = n-l-1
                Hmat += (r^d) * eye(2^(l-1)) ⊗ (delta * (sz ⊗ eye(2^d) ⊗ sz) +
                    0.5 * (sp ⊗ eye(2^d) ⊗ sm + sm ⊗eye(2^d) ⊗ sp))
            end
        end
    end
    Hmat
end
