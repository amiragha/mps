sz_half = [0.5 0.0; 0.0 -.5]
@testset "MPS" begin
    d, lx = 2, 8

    # A uniformly chosen random state as MPS
    ketstate = complex.(randn(d^lx), randn(d^lx))
    randmps = MPS(lx, d, ketstate)

    # The Bethe chain ground state as MPS
    H = xxz_hamiltonian(lx)
    eheis, vheis = eigsolve(H, 1, :SR)
    mps = MPS(lx, 2, vheis[1])

    @testset "ketstate to MPS to ketstate" begin
        @test norm(ketstate) ≈ norm(randmps)
        @test ketstate ≈ mps2ketstate(randmps)
    end

    @testset "overlap of two MPS" begin
        @test sqrt(overlap(mps, mps)) ≈ norm(mps)
        @test overlap(mps, randmps) ≈ conj.(overlap(randmps, mps))
    end

    @testset "operator measurements" begin
        sz = spinoperators(1/2)[1]
        szdata = measure(mps, sz)
        @test all([isapprox(szdata[i], 0.0, atol=1.e-12) for i=1:lx])

        mpo = xxz_mpo(Float64, lx, 2)
        @test measure(mps, mpo)[1] ≈ eheis[1]
    end

    @testset "MPS center manipulations" begin
        ### TODO: write more tests for here!
        center_at!(randmps, lx-1)
        center_at!(randmps, lx-2)
        center_at!(randmps, 1)
        center_at!(randmps, lx)
        @test mps2ketstate(randmps) ≈ ketstate

        center_at!(randmps, lx)
        @test mps2ketstate(randmps) ≈ ketstate
    end

end

@testset "apply" begin
    @testset "twosite" begin
        lx, l = 8, 3

        H = xxz_hamiltonian(lx)
        eheis, vheis = eigsolve(H, 1, :SR)
        mps = MPS(lx, 2, vheis[1])

        mpscopy = deepcopy(mps)
        center_at!(mps, l)
        sz = spinoperators(1/2)[1]
        M = reshape(kron(sz, sz), 2, 2, 2, 2)
        apply!(mps, M, l)

        #@test mps_dims_are_consistent(mps)
        @test overlap(mpscopy, mps) ≈ measure(mpscopy, sz, sz, l, l+1)
    end
end
