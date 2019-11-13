@testset "Model" begin

    # @testset "fermion" begin
    #     lx = 4
    #     t1 = 1.0
    #     t2 = 0.7
    #     mu = 0.3

    #     uc = triangular_unitcell(t1, t2, 0.3)

    #     ly = 1
    #     H1 = sparse(makemodel(uc, lx, ly, boundary=:PBCX))
    #     @test H1 ≈ hopping_chain(lx*ly, t2, mu, boundary=:PBC)

    #     ly = 2
    #     H2 = sparse(makemodel(uc, lx, ly, boundary=:PBCX))
    #     @test H2 ≈ nnhoppingchain(lx*ly, t1, t2, mu, boundary=:PBC)
    # end

    # @testset "JW" begin

    # end

    @testset "opexpansion" begin
        terms = nbodyopexpansion(4, ringop)
        h = zeros(16,16)
        for t in terms
            h += reduce(⊗, t, init=I(1))
        end
        @test h == ringop
    end

    @testset "ringop" begin
        symmetry = :NONE
        k = 0.3
        R4 = nbodyopexpansion(4,
                              ringexchangeoperator(4) - 0.25 * I(16),
                              symmetry=symmetry)

        ring1 = QModelInteraction{1, 4, Float64}(
            0.3,
            (1, 1, 1, 1),
            ((0,), (1,), (2,), (3,)),
            R4)

        model = UnitCellQModel{SpinType, 1}(spinhalf,
                                            QLattice(chainunitcell, lx, :OBC),
                                            [ring1])
        mpo = generatempo(model)
        hmpo = mpo2hamiltonian(mpo)

        @test hmpo + 0.25 *k* I(16) ≈ k*ringop
    end
end
