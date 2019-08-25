@testset "releg" begin

    @testset "generic" begin
        legs = STLeg(+1, [0,1], [2,3]), STLeg(+1, [0,1], [1,1]), STLeg(-1, [0,1], [3,4])

        sects = [(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)]
        nzblks = [rand(1:9, 2,1,3), rand(1:9, 2,1,4), rand(1:9, 3,1,4), rand(1:9, 3,1,3)]

        perm = _sectors_sortperm(sects)

        test = SymTensor(0, legs, sects[perm], nzblks[perm])
        ftest = fuselegs(test, +1, 1, 2)
        dftest = defuse_leg(ftest, 1, legs[1:2])
        @test dftest == test
    end

    @testset "2siteop s=1/2" begin
        smat = eye(Float64, 0, [0,1,2], [1,2,1])
        change_nzblk!(smat, (1,1), [2. 3;4 5])
        rlegs = (STLeg(+1, [0,1], [1,1]), STLeg(+1, [0,1], [1,1]))
        clegs = (STLeg(-1, [0,1], [1,1]), STLeg(-1, [0,1], [1,1]))
        smat2 = defuse_leg(defuse_leg(smat, 2, clegs), 1, rlegs)
        @test smat2.sects == [(0,0,0,0), (1,0,1,0), (0,1,1,0),
                              (1,0,0,1), (0,1,0,1), (1,1,1,1)]
        i = ones(1,1,1,1)
        @test smat2.nzblks == [i,2*i, 4*i, 3*i, 5*i, 1*i]
    end
end

@testset "contract" begin

    @testset "simple" begin
        leglist = (
            STLeg(+1, [0, 1, 2], [2,3,4]),
            STLeg(-1, [0,1,2,3], [4,2,2,3]),
            STLeg(+1, [0,1,2,3], [4,2,2,3]),
            STLeg(-1, [0, 1, 2], [3,5,2]),
        )
        sten1 = randSymTensor(Float64, 0, leglist[1:2])
        sten2 = fillSymTensor(1.0, 0, leglist[3:4])

        arr1 = array_representation(sten1)
        arr2 = array_representation(sten2)

        contest = contract(sten1, (1, -1), sten2, (-1, 2))

        @test arr1*arr2 == array_representation(contest)
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
        sten1 = randSymTensor(Float64, 0, leglist[1:4]);
        sten2 = randSymTensor(Float64, 0, leglist[5:8]);

        arr1 = array_representation(sten1)
        arr2 = array_representation(sten2)

        contest = contract(sten1, (1,2,3,-1), sten2, (-1,4,5,6));
        @tensor arr[a,b,c,e,f,g] := arr1[a,b,c,d] * arr2[d,e,f,g];

        @test arr == array_representation(contest);
    end
end
