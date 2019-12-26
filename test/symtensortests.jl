@testset "releg" begin

    @testset "fuse/unfuse" begin
        space = (
            U1Space(0=>2, 1=>3),
            U1Space(0=>1, 1=>1),
            U1Space(0=>3, -1=>4),
        )

        # Note that all the possible sectors are given!
        sects = [Sector(0, 0, 0), Sector(0, 1, -1), Sector(1, 0, -1)]
        nzblks = [reshape(collect(1:6), 2,1,3),
                  reshape(collect(7:14), 2,1,4),
                  reshape(collect(15:26), 3,1,4)]

        A0 = SymTensor(0, space, SortedDict(zip(sects,nzblks)))
        A1 = fuselegs(A0, 1, 2)
        A2 = splitleg(A1, 1, space[1:2])
        @test A0 == A2
    end

    @testset "U1 unitary" begin
        A = eye(Float64, [0,1,2], [1,2,1])
        set_sector!(A, (1,-1), [2. 3;4 5])
        rlegs = (U1Space(+1, [0,1], [1,1]), U1Space(+1, [0,1], [1,1]))
        clegs = (U1Space(-1, [0,1], [1,1]), U1Space(-1, [0,1], [1,1]))
        B = unfuseleg(unfuseleg(A, 2, clegs), 1, rlegs)
        @test B.sects == [(1,1,-1,-1),(1,0,0,-1), (0,1,0,-1),
                          (1,0,-1,0), (0,1,-1,0), (0,0,0,0)]
        i = ones(1,1,1,1)
        @test B.nzblks == [i,2*i, 4*i, 3*i, 5*i, 1*i]
    end

    @testset "Associativity" begin
        legs = (
            U1Space(+1, [-1, 0, 1], [1, 2, 1]),
            U1Space(-1, [-1, 0, 1], [1, 2, 1]),
            U1Space(+1, [0, 1, 2], [1, 2, 1]),
            U1Space(-1, [0, 1, 2], [1, 2, 1]),
            # U1Space(+1, [-2,-1, 0], [1, 2, 3]),
            # U1Space(-1, [-2,-1, 0], [1, 2, 3]),
            # U1Space(+1, [0], [1, 2, 3]),
            # U1Space(-1, [0], [1, 2, 3]),
        )

        A = rand(0,  legs[1:4])
        A1 = fuselegs(A, -1, 2, 3)
        A2 = fuselegs(A, -1, 3, 2)
        A3 = fuselegs(A, +1, 1, 3)
        A4 = fuselegs(A, +1, 1, 2)

        @test A2 ≈ unfuseleg(A1, 2, legs[2], fuse(-1, legs[3], legs[4]))
        @test A4 ≈ unfuseleg(A3, 1, fuse(+1, legs[1], legs[2]), legs[3])
        @test A1 ≈ fuselegs(A2, -1, 2, 2)
        @test A3 ≈ fuselegs(A4, +1, 1, 2)
    end
end

@testset "contract" begin

    @testset "matrices" begin
        leglist = (
            U1Space(+1, [0, 1, 2], [2,3,4]),
            U1Space(-1, [0,1,2,3], [4,2,2,3]),
            U1Space(+1, [0,1,2,3], [4,2,2,3]),
            U1Space(-1, [0, 1, 2], [3,5,2]),
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
            U1Space(+1, [0, 1, 2], [2,3,4]),
            U1Space(-1, [0,1,2,3], [4,2,2,3]),
            U1Space(+1, [0, 1, 2], [3,5,2]),
            U1Space(-1, [0, 1], [6,7]),
            U1Space(+1, [0, 1], [6,7]),
            U1Space(+1, [0,1,2,3,4], [2,3,4,5,6]),
            U1Space(+1, [0,1], [2,3]),
            U1Space(-1, [0,1,2], [3,4,5])
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
        U1Space(+1, [0, 1, 2], [2,3,4]),
        U1Space(-1, [0,1,2,3], [4,2,2,3]),
        U1Space(+1, [0, 1, 2], [3,5,2]),
        U1Space(-1, [0, 1], [6,7]),
        U1Space(+1, [0, 1], [6,7]),
        U1Space(+1, [0,1,2,3,4], [2,3,4,5,6]),
        U1Space(+1, [0,1], [2,3]),
        U1Space(-1, [0,1,2], [3,4,5])
    )
    @testset "dot" begin
        A = rand(ComplexF64, 0, leglist[1:4])
        A_ = array(A)
        @test dot(A_,A_) ≈ dot(A, A)
    end
end
