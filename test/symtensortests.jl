@testset "releg" begin

    @testset "fuse/split" begin
        space = (
            U1Space(0=>2, 1=>3),
            U1Space(0=>1, 1=>1),
            dual(U1Space(0=>3, 1=>4)),
        )

        # Note that all the possible sectors are given!
        sects = [Sector{U1}(0, 0, 0),
                 Sector{U1}(0, 1, 1),
                 Sector{U1}(1, 0, 1)]
        nzblks = [reshape(collect(1:6), 2,1,3),
                  reshape(collect(7:14), 2,1,4),
                  reshape(collect(15:26), 3,1,4)]

        A0 = SymTensor(zero(U1), space, SortedDict(zip(sects,nzblks)))
        #A0 = SymTensor(zero(U1), space,
        #SortedDict{Sector{U1,3}, Array{Int,3}}(zip(sects,nzblks)))

        A1 = fuselegs(A0, 1, 2)
        A2 = splitleg(A1, 1, space[1:2])
        @test A0 == A2

        A1 = fuselegs(A0, 2, 2, true)
        A2 = splitleg(A1, 2, space[2:3])
        @test A0 == A2
    end

    @testset "dual" begin
        A = eye(U1Space(0=>1, 1=>2, 2=>1))
        A[Sector{U1}(2, 2)] *= -1
        A[Sector{U1}(1, 1)] *= 1/sqrt(2) .* [-1. 1; 1 1]

        B = splitleg(A, 1, U1Space(0=>1, 1=>1), U1Space(0=>1, 1=>1))
        B = splitleg(A, 1, U1Space(0=>1, 1=>1), U1Space(0=>1, 1=>1))
        fuselegs(dual(B), 1, 2, false) == dual(A)
    end

    @testset "U1 unitary" begin
        A = eye(Float64, U1Space(0=>1, 1=>2, 2=>1))
        A[Sector{U1}(1, 1)] =  [2. 3;4 5]
        d = U1Space(0=>1, 1=>1)
        cod = (d, d)
        dom = (dual(d), dual(d))
        B = splitleg(splitleg(A, 2, dom), 1, cod)
        @test sectors(B) == sort([Sector{U1}(1,1,1,1), Sector{U1}(1,0,0,1),
                                  Sector{U1}(0,1,0,1), Sector{U1}(1,0,1,0),
                                  Sector{U1}(0,1,1,0), Sector{U1}(0,0,0,0)])
        i = ones(1,1,1,1)
        @test collect(values(B.blocks)) == [i,2*i, 4*i, 3*i,5*i, 1*i]
    end

    @testset "Associativity" begin
        d = U1Space(-1=>1, 0=>2, 1=>1)
        space = (d, dual(d), d, dual(d))

        A = fill_linearindex(zero(U1),  space[1:4])
        A1 = fuselegs(A, 2, 3)
        A2 = fuselegs(A, 3, 2)
        A3 = fuselegs(A, 1, 3)
        A4 = fuselegs(A, 1, 2)
        A5 = fuselegs(A, 2 ,2)

        #@test A2 ≈ splitleg(A1, 2, space[2], fuse(space[3], space[4]))
        #@test A4 ≈ splitleg(A3, 1, fuse(space[1], space[2]), space[3])
        #@test A1 ≈ fuselegs(A2, 2, 2)
        #@test A1 ≈ fuselegs(A5, 2, 2)
        #@test A3 ≈ fuselegs(A4, 1, 2)
    end
end

@testset "contract" begin

    @testset "converts" begin
        vlist = (
            U1Space(0=>2, 1=>3, 2=>4),
            U1Space(-2=>4, -1=>3, 0=>2, 1=>1),
            U1Space(-1=>1, 0=>2, 1=>3, 2=>4),
            U1Space(-2=>3, -1=>2, 0=>1),
        )
        A = fill_linearindex(vlist)
        @test SymMatrix(A, [1,2], [3,4]) == fuselegs(fuselegs(A, 3, 2, true), 1, 2)
        @test SymMatrix(A, [1,2], [3,4]) == fuselegs(fuselegs(A, 1, 2), 2, 2, true)
        #@test SymVector(A, [1,2,3,4]) == fuselegs(A,1,4)

        perm = randperm(4)
        @test SymMatrix(A, perm[1:2], perm[3:4]) ==
            fuselegs(fuselegs(permutelegs(A, perm), 3, 2, true), 1, 2)
        #@test SymVector(A, perm) == fuselegs(permutelegs(A, perm), 1,  4)
    end

    @testset "matrices" begin
        vlist = (
            U1Space(0=>2, 1=>3, 2=>4),
            dual(U1Space(0=>4, 1=>2, 2=>2, 3=>3)),
            U1Space(0=>4, 1=>2, 2=>2, 3=>3),
            dual(U1Space(0=>3, 1=>5, 2=>2)),
        )
        A = fill_linearindex(vlist[1:2])
        B = fill(1.0, zero(U1), vlist[3:4])

        A_ = array(A)
        B_ = array(B)

        C = contract(A, (1, -1), B, (-1, 2))

        @test A_*B_ == array(C)
    end

    @testset "multi leg tensors" begin
        vlist = (
            U1Space(0=>2, 1=>2, 2=>2),
            U1Space(0=>2, -1=>2, -2=>2, -3=>2),
            U1Space(0=>2, 1=>2, 2=>2),
            dual(U1Space(0=>2, 1=>2)),
            U1Space(0=>2, 1=>2),
            U1Space(0=>2, 1=>2, 2=>2, 3=>2, 4=>2),
            U1Space(0=>2, 1=>2),
            U1Space(0=>2, -1=>2, -2=>2)
        )
        A = rand(Float64, zero(U1), vlist[1:4])
        B = rand(Float64, zero(U1), vlist[5:8])

        @test splitleg(fuselegs(A, 1, 3), 1, vlist[1:3]) ≈ A
        @test splitleg(fuselegs(A, 2, 2), 2, vlist[2:3]) ≈ A
        @test splitleg(fuselegs(A, 2, 3), 2, vlist[2:4]) ≈ A

        A_ = array(A)
        B_ = array(B)

        C = contract(A, (1,2,3,-1), B, (-1,4,5,6))
        @tensor C_[a,b,c,e,f,g] := A_[a,b,c,d] * B_[d,e,f,g]

        @test C_ ≈ array(C);
    end
end

@testset "LA fns" begin
    vlist = (
        U1Space(0=>2, 1=>3, 2=>4),
        U1Space(0=>4, -1=>2, -2=>2, -3=>3),
        U1Space(0=>3, 1=>5, 2=>2),
        U1Space(0=>6, -1=>7),
        U1Space(0=>6, 1=>7),
        U1Space(0=>2, 1=>3, 2=>4, 3=>5, 4=>6),
        U1Space(0=>2, 1=>3),
        U1Space(0=>3, -1=>4, -2=>5)
    )
    @testset "dot" begin
        A = rand(ComplexF64, zero(U1), vlist[1:4])
        A_ = array(A)
        @test dot(A_,A_) ≈ dot(A, A)
    end
end

@testset "factorization" begin

    V1 = U1Space(0=>2, 1=>3, 2=>4)
    V2 = U1Space(0=>4, -1=>2, -2=>2, -3=>3)
    V3 = U1Space(0=>4, 1=>2, 2=>2, 3=>3)
    V4 = U1Space(0=>3, -1=>5, -2=>2)

    @testset "svd" begin
        A1 = rand(Float64, zero(U1), (V1,V2))
        A2 = rand(Float64, U1(1), (V1,V2))
        A3 = rand(Float64, U1(-1), (V1,V2))
        u,s,v = _svd_(A1)
        @test u*s*v ≈ A1
        u,s,v = _svd_(A2)
        @test u*s*v ≈ A2
        u,s,v = _svd_(A3)
        @test u*s*v ≈ A3

        B1 = rand(ComplexF64, zero(U1), (V1,dual(V3)))
        B2 = rand(ComplexF64, U1(1), (V1,dual(V3)))
        B3 = rand(ComplexF64, U1(-1), (V1,dual(V3)))
        u,s,v = _svd_(B1)
        @test u*s*v ≈ B1
        u,s,v = _svd_(B2)
        @test u*s*v ≈ B2
        u,s,v = _svd_(B3)
        @test u*s*v ≈ B3
    end

    @testset "svdtrunc" begin
        A1 = SymTensor(rand, ComplexF64, zero(U1), (V1,V2))
        @test _svd_(A1) ≈ svdtrunc(A1)
    end
end
