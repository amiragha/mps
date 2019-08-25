sz_half = [0.5 0.0; 0.0 -.5]
@testset "MPS" begin
    d, lx = 2, 8

    # A uniformly chosen random state as MPS
    ketstate = complex.(randn(d^lx), randn(d^lx))
    randmps = MatrixProductState(lx, d, ketstate)

    # The Bethe chain ground state as MPS
    H = xxz_hamiltonian(lx)
    eheis, vheis = eigsolve(H, 1, :SR)
    mps = MatrixProductState(lx, 2, vheis[1])

    @testset "ketstate to MPS to ketstate" begin
        @test norm(ketstate)^2 ≈ norm2(randmps)
        @test ketstate ≈ mps2ketstate(randmps)
    end

    @testset "overlap of two MPS" begin
        @test overlap(mps, mps) ≈ norm2(mps)
        @test overlap(mps, randmps) ≈ conj.(overlap(randmps, mps))
    end

    @testset "operator measurements" begin
        @test measure_1point(mps, sz_half) ≈ zeros(lx) atol=1.e-12

        mpo = xxz_mpo(Float64, lx, 2)
        @test measure_mpo(mps, mpo) ≈ eheis[1]
    end

    @testset "MPS center manipulations" begin
        ### TODO: write more tests for here!
        move_center!(randmps, lx-1)
        move_center!(randmps, lx-2)
        move_center!(randmps, 1)
        move_center!(randmps, lx)
        @test mps2ketstate(randmps) ≈ ketstate

        canonicalize_at!(randmps.matrices, lx)
        @test mps2ketstate(randmps) ≈ ketstate
    end

end

@testset "apply" begin
    @testset "twosite" begin
        lx, l = 8, 3

        H = xxz_hamiltonian(lx)
        eheis, vheis = eigsolve(H, 1, :SR)
        mps = MatrixProductState(lx, 2, vheis[1])

        mpscopy = deepcopy(mps)
        move_center!(mps, l)
        M = twosite_tensor(sz_half, sz_half)
        apply_2siteoperator!(mps, l, M)

        @test mps_dims_are_consistent(mps)
        @test overlap(mpscopy, mps) ≈ measure_2point(mpscopy, sz_half, sz_half, l, l+1)
    end
end
