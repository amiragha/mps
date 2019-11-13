@testset "MPO" begin

    @testset "xxz" begin
        smpo = xxz_symmpo(Float64, 4, 2, 1.0)
        mpo = xxz_mpo(Float64, 4, 2, 1.0);

        h2sym = contract(smpo.tensors[1], (1,2,-1,5),
                         smpo.tensors[4], (-1,3,4,6))
        h2sym = fuselegs(fuselegs(h2sym,+1, 2, 2), -1, 4, 2)
        h2sym = removedummyleg(h2sym, 3)
        h2sym = removedummyleg(h2sym, 1)

        h4sym = contract(
            contract(
                contract(smpo.tensors[1], (1,2,-1,5),
                         smpo.tensors[2], (-1,3,4,6)), (1,2,3,-1,6,7),
                smpo.tensors[3], (-1,4,5,8)), (1,2,3,4,-1,7,8,9),
            smpo.tensors[4], (-1,5,6,10))

        h4sym = fuselegs(fuselegs(h4sym, +1, 2, 4), -1, 4, 4)
        h4sym = removedummyleg(h4sym, 3)
        h4sym = removedummyleg(h4sym, 1)

        @tensor h2[:] :=
            mpo.tensors[1][-1,-2,1,-5] * mpo.tensors[4][1,-3,-4,-6]
        h2 = reshape(h2,4,4)

        @tensor h3[:] :=
            mpo.tensors[1][-1,-2,1,-6] * mpo.tensors[2][1,-3,2,-7] *
            mpo.tensors[4][2,-4,-5,-8]
        h3 = reshape(h3,8,8)

        @tensor h4[:] :=
            mpo.tensors[1][-1,-2,1,-7] * mpo.tensors[2][1,-3,2,-8] *
            mpo.tensors[3][2,-4,3,-9] * mpo.tensors[4][3,-5,-6,-10]
        h4 = reshape(h4,16,16)

        H2 = Matrix(xxz_hamiltonian(2, 1.0))
        H3 = Matrix(xxz_hamiltonian(3, 1.0))
        H4 = Matrix(xxz_hamiltonian(4, 1.0))

        @test h2 ≈ H2
        @test h3 ≈ H3
        @test h4 ≈ H4

        nums = collect(0:15)
        index4 = count_ones.(nums) .== 4
        index3 = count_ones.(nums) .== 3
        index2 = count_ones.(nums) .== 2
        index1 = count_ones.(nums) .== 1
        index0 = count_ones.(nums) .== 0
        @test h4[index0,index0] == h4sym.nzblks[1]
        @test h4[index1,index1] == h4sym.nzblks[2]
        @test h4[index2,index2] == h4sym.nzblks[3]
        @test h4[index3,index3] == h4sym.nzblks[4]
        @test h4[index4,index4] == h4sym.nzblks[5]
    end


    @testset "xxzlong" begin
        delta, r = 2.0, 0.6

        mpo = xxzlong_mpo(Float64,4,2,delta,r)
        @tensor h2[:] :=
            mpo.tensors[1][-1,-2,1,-5] * mpo.tensors[4][1,-3,-4,-6]
        h2 = reshape(h2,4,4)

        @tensor h3[:] :=
            mpo.tensors[1][-1,-2,1,-6] * mpo.tensors[2][1,-3,2,-7] *
            mpo.tensors[4][2,-4,-5,-8]
        h3 = reshape(h3,8,8)

        @tensor h4[:] :=
            mpo.tensors[1][-1,-2,1,-7] * mpo.tensors[2][1,-3,2,-8] *
            mpo.tensors[3][2,-4,3,-9] * mpo.tensors[4][3,-5,-6,-10]
        h4 = reshape(h4,16,16)

        H2 = xxz_longrange(2,delta,r)
        H3 = xxz_longrange(3,delta,r)
        H4 = xxz_longrange(4,delta,r)

        @test h2 ≈ Matrix(H2)
        @test h3 ≈ Matrix(H3)
        @test h4 ≈ Matrix(H4)
    end
end

@testset "MPOgen" begin

    @testset "j1j2" begin
        n = 4
        lx = 2*n
        j1, j2 = 1.0, 0.3
        model = j1j2model(lx, j1, j2)
        mpo = generatempo(model)
        legacympo = MatrixProductStateTools.j1j2_mpo(lx, j1, j2, 2)
        @test mpo2hamiltonian(legacympo) ≈ mpo2hamiltonian(mpo)

        tmodel = triangularspinmodel((2, n),
                                     0.3, 1.0, 1.0, 0.0, 0.0, 0.0)
        tmpo = generatempo(tmodel)
        @test mpo2hamiltonian(tmpo) ≈ mpo2hamiltonian(mpo)
    end

    @testset "j1j2sym" begin
        j1, j2 = 1.0, 0.3
        for lx in [3, 4]
            smodel = j1j2model(lx, j1, j2, symmetry=:U1)
            smpo = generatesymmpo(smodel)

            model = j1j2model(lx, j1, j2)
            mpo = generatempo(model)

            sH = mpo2hamiltonian(smpo)
            H = mpo2hamiltonian(mpo)

            nums = collect(0:2^lx-1)
            indexes = [count_ones.(nums) .== m for m in 0:lx]
            for n = 1:lx
                @test H[indexes[n], indexes[n]] ≈ sH.nzblks[n]
            end
        end

        for n in [2, 3]
            smodel = j1j2model(2*n, j1, j2, symmetry=:U1)
            smpo = generatesymmpo(smodel)

            tmodel = triangularspinmodel((2, n),
                                         0.3, 1.0, 1.0, 0.0, 0.0, 0.0,
                                         symmetry=:U1)
            tmpo = generatesymmpo(tmodel)
            @test mpo2hamiltonian(tmpo) ≈ mpo2hamiltonian(smpo)
        end
    end

    @testset "2legtriangle" begin
        ⊗ = kron
        n = 4
        j1, j2, k1, k2 = 1.0, 0.3, 0.4, 0.7
        k = 0.3
        H = j1j2_explicit(2*n, j1, j2)
        for i=1:2*n-3
            if i % 2 == 1
                H += k1 * I(2^(i-1)) ⊗ permutespins(ringop, [1,3,4,2]) ⊗ I(2^(2*n-i-3))
            else
                H += k2 * I(2^(i-1)) ⊗ permutespins(ringop, [1,3,4,2]) ⊗ I(2^(2*n-i-3))
            end
        end

        model = triangularspinmodel((2, n), j2, j1, j1, k1, k2, 0.)
        mpo = generatempo(model)
        Hmpo = mpo2hamiltonian(mpo)

        nk1s = div(2*n-3, 2) + (2*n-3) % 2
        nk2s = div(2*n-3, 2)
        @test Hmpo + 0.25*I(2^(2*n))*k1*nk1s + 0.25*I(2^(2*n))*k2*nk2s ≈ H
    end

    @testset "2legtrianglesym" begin
        j1, j2, j3, k1, k2, k3 = rand(6)
        for n in [2, 3]
            smodel = triangularspinmodel((2,n), j1, j2, j3, k1, k2, k3, symmetry=:U1)
            smpo = generatesymmpo(smodel)

            model = triangularspinmodel((2,n), j1, j2, j3, k1, k2, k3)
            mpo = generatempo(model)

            sH = mpo2hamiltonian(smpo)
            H = mpo2hamiltonian(mpo)

            lx = 2*n
            nums = collect(0:2^lx-1)
            indexes = [count_ones.(nums) .== m for m in 0:lx]
            for n = 1:lx
                @test H[indexes[n], indexes[n]] ≈ sH.nzblks[n]
            end
        end
    end
end
