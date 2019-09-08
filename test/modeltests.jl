@testset "Model" begin

    @testset "fermion" begin
        lx = 4
        t1 = 1.0
        t2 = 0.7
        mu = 0.3

        uc = triangular_unitcell(t1, t2, 0.3)

        ly = 1
        H1 = sparse(makemodel(uc, lx, ly, boundary=:PBCX))
        @test H1 ≈ hopping_chain(lx*ly, t2, mu, boundary=:PBC)

        ly = 2
        H2 = sparse(makemodel(uc, lx, ly, boundary=:PBCX))
        @test H2 ≈ nnhoppingchain(lx*ly, t1, t2, mu, boundary=:PBC)
    end

    @testset "JW" begin

    end
end
