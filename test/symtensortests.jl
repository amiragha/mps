@testset "releg" begin

    @testset "generic" begin
        legs = STLeg(+1, [0,1], [2,3]), STLeg(+1, [0,1], [1,1]), STLeg(-1, [0,1], [3,4])

        # Note that all the possible sectors are given!
        sects = [(0, 0, 0), (0, 1, 1), (1, 0, 1)]
        nzblks = [rand(1:9, 2,1,3), rand(1:9, 2,1,4), rand(1:9, 3,1,4)]

        perm = _sectors_sortperm(sects)

        test = SymTensor(0, legs, sects[perm], nzblks[perm])
        ftest = fuselegs(test, +1, 1, 2)
        dftest = unfuseleg(ftest, 1, legs[1:2])
        @test dftest == test
    end

    @testset "2siteop s=1/2" begin
        A = eye(Float64, [0,1,2], [1,2,1])
        set_sector!(A, (1,1), [2. 3;4 5])
        rlegs = (STLeg(+1, [0,1], [1,1]), STLeg(+1, [0,1], [1,1]))
        clegs = (STLeg(-1, [0,1], [1,1]), STLeg(-1, [0,1], [1,1]))
        B = unfuseleg(unfuseleg(A, 2, clegs), 1, rlegs)
        @test B.sects == [(0,0,0,0), (1,0,1,0), (0,1,1,0),
                          (1,0,0,1), (0,1,0,1), (1,1,1,1)]
        i = ones(1,1,1,1)
        @test B.nzblks == [i,2*i, 4*i, 3*i, 5*i, 1*i]
    end
end

@testset "contract" begin

    @testset "matrices" begin
        leglist = (
            STLeg(+1, [0, 1, 2], [2,3,4]),
            STLeg(-1, [0,1,2,3], [4,2,2,3]),
            STLeg(+1, [0,1,2,3], [4,2,2,3]),
            STLeg(-1, [0, 1, 2], [3,5,2]),
        )
        A = rand(Float64, 0, leglist[1:2])
        B = fill(1.0, 0, leglist[3:4])

        A_ = array(A)
        B_ = array(B)

        C = contract(A, (1, -1), B, (-1, 2))

        @test A_*B_ == array(C)
    end

    @testset "multi leg tensors" begin
        leglist = (
            STLeg(+1, [0, 1, 2], [2,3,4]),
            STLeg(-1, [0,1,2,3], [4,2,2,3]),
            STLeg(+1, [0, 1, 2], [3,5,2]),
            STLeg(-1, [0, 1], [6,7]),
            STLeg(+1, [0, 1], [6,7]),
            STLeg(+1, [0,1,2,3,4], [2,3,4,5,6]),
            STLeg(+1, [0,1], [2,3]),
            STLeg(-1, [0,1,2], [3,4,5])
        )
        A = rand(Float64, 0, leglist[1:4])
        B = rand(Float64, 0, leglist[5:8])

        @test unfuseleg(fuselegs(A, +1, 1, 3), 1, leglist[1:3]) ≈ A

        A_ = array(A)
        B_ = array(B)

        C = contract(A, (1,2,3,-1), B, (-1,4,5,6))
        @tensor C_[a,b,c,e,f,g] := A_[a,b,c,d] * B_[d,e,f,g]

        @test C_ ≈ array(C);
    end
end

@testset "LA fns" begin
    leglist = (
        STLeg(+1, [0, 1, 2], [2,3,4]),
        STLeg(-1, [0,1,2,3], [4,2,2,3]),
        STLeg(+1, [0, 1, 2], [3,5,2]),
        STLeg(-1, [0, 1], [6,7]),
        STLeg(+1, [0, 1], [6,7]),
        STLeg(+1, [0,1,2,3,4], [2,3,4,5,6]),
        STLeg(+1, [0,1], [2,3]),
        STLeg(-1, [0,1,2], [3,4,5])
    )
    @testset "dot" begin
        A = rand(ComplexF64, 0, leglist[1:4])
        A_ = array(A)
        @test dot(A_,A_) ≈ dot(A, A)
    end
end
