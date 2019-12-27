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
        A = eye(Float64, U1Space(0=>1, 1=>2, 2=>1))
        A[Sector(1,-1)] =  [2. 3;4 5]
        d = U1Space(0=>1, 1=>1)
        cod = (d, d)
        dom = (dual(d), dual(d))
        B = splitleg(splitleg(A, 2, dom), 1, cod)
        @test sectors(B) == [Sector(1,1,-1,-1), Sector(1,0,0,-1),
                             Sector(0,1,0,-1), Sector(1,0,-1,0),
                             Sector(0,1,-1,0), Sector(0,0,0,0)]
        i = ones(1,1,1,1)
        @test collect(values(B.blocks)) == [i,2*i, 4*i, 3*i,5*i, 1*i]
    end

    @testset "Associativity" begin
        d = U1Space(-1=>1, 0=>2, 1=>1)
        space = (d, dual(d), d, dual(d))

        A = fill_linearindex(0,  space[1:4])
        A1 = fuselegs(A, 2, 3)
        A2 = fuselegs(A, 3, 2)
        A3 = fuselegs(A, 1, 3)
        A4 = fuselegs(A, 1, 2)

        #@test A2 ≈ splitleg(A1, 2, space[2], fuse(space[3], space[4]))
        #@test A4 ≈ splitleg(A3, 1, fuse(space[1], space[2]), space[3])
        #@test A1 ≈ fuselegs(A2, 2, 2)
        #@test A3 ≈ fuselegs(A4, 1, 2)
    end
end

@testset "contract" begin

    @testset "matrices" begin
        vlist = (
            U1Space(0=>2, 1=>3, 2=>4),
            U1Space(0=>4, 1=>2, 2=>2, 3=>3),
            U1Space(0=>4, 1=>2, 2=>2, 3=>3),
            U1Space(0=>3, 1=>5, 2=>2),
        )
        A = rand(Float64, 0, vlist[1:2])
        B = fill(1.0, 0, vlist[3:4])

        A_ = array(A)
        B_ = array(B)

        C = contract(A, (1, -1), B, (-1, 2))

        @test A_*B_ == array(C)
    end

    @testset "multi leg tensors" begin
        vlist = (
            U1Space(+1, [0, 1, 2], [2,3,4]),
            U1Space(-1, [0,1,2,3], [4,2,2,3]),
            U1Space(+1, [0, 1, 2], [3,5,2]),
            U1Space(-1, [0, 1], [6,7]),
            U1Space(+1, [0, 1], [6,7]),
            U1Space(+1, [0,1,2,3,4], [2,3,4,5,6]),
            U1Space(+1, [0,1], [2,3]),
            U1Space(-1, [0,1,2], [3,4,5])
        )
        A = rand(Float64, 0, vlist[1:4])
        B = rand(Float64, 0, vlist[5:8])

        @test unfuseleg(fuselegs(A, +1, 1, 3), 1, vlist[1:3]) ≈ A

        A_ = array(A)
        B_ = array(B)

        C = contract(A, (1,2,3,-1), B, (-1,4,5,6))
        @tensor C_[a,b,c,e,f,g] := A_[a,b,c,d] * B_[d,e,f,g]

        @test C_ ≈ array(C);
    end
end

@testset "LA fns" begin
    vlist = (
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
        A = rand(ComplexF64, 0, vlist[1:4])
        A_ = array(A)
        @test dot(A_,A_) ≈ dot(A, A)
    end
end
