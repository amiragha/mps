@testset "MPO" begin

    smpo = xxz_symmpo(Float64, 4, 2, 1.0)
    mpo = xxz_mpo(Float64, 4, 2, 1.0);

    h2ten = contract(smpo.tensors[1], (1,2,-1,5), smpo.tensors[4], (-1,3,4,6))
    h2ten = fuselegs(fuselegs(h2ten,+1, 2, 2), -1, 4, 2)
    h2ten = removedummyleg(h2ten, 3)
    h2ten = removedummyleg(h2ten, 1)

    h4ten = contract(
        contract(
            contract(smpo.tensors[1], (1,2,-1,5), smpo.tensors[2], (-1,3,4,6)),
            (1,2,3,-1,6,7), smpo.tensors[3], (-1,4,5,8)),
        (1,2,3,4,-1,7,8,9), smpo.tensors[4], (-1,5,6,10))

    h4ten = fuselegs(fuselegs(h4ten, +1, 2, 4), -1, 4, 4)
    h4ten = removedummyleg(h4ten, 3)
    h4ten = removedummyleg(h4ten, 1)

    @tensor h2[:] := mpo.tensors[1][-1,-2,1,-5] * mpo.tensors[4][1,-3,-4,-6]
    h2 = reshape(h2,4,4)

    @tensor h4[:] := mpo.tensors[1][-1,-2,1,-7] * mpo.tensors[2][1,-3,2,-8] *
        mpo.tensors[3][2,-4,3,-9] * mpo.tensors[4][3,-5,-6,-10]
    h4 = reshape(h4,16,16)

    @testset "mpo to Hamiltonian" begin
        @test h4 ≈ Matrix(xxz_hamiltonian(4, 1.0))
    end

    @testset "mpo symmpo" begin
        nums = collect(0:15)
        index4 = count_ones.(nums) .== 4
        index3 = count_ones.(nums) .== 3
        index2 = count_ones.(nums) .== 2
        index1 = count_ones.(nums) .== 1
        index0 = count_ones.(nums) .== 0
        @test h4[index0,index0] == h4ten.nzblks[1]
        @test h4[index1,index1] == h4ten.nzblks[2]
        @test h4[index2,index2] == h4ten.nzblks[3]
        @test h4[index3,index3] == h4ten.nzblks[4]
        @test h4[index4,index4] == h4ten.nzblks[5]
    end
end
