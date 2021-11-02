using Network_ADER
using Test

@testset "Computation of c constants" begin
    # Testing right boundary
    # zero order right boundary
    @test Network_ADER.compute_c(1,0,0,0) == 1
    # 3 rd order right boundary
    @test Network_ADER.compute_c(3,1,1,2) == 1
    @test Network_ADER.compute_c(3,2,0,2) == 11/6
    
    # Testing left boundary
    # left boundary seond order
    @test Network_ADER.compute_c(2,0,0,1,true) == -1/2
    # left boundary 4th order
    @test Network_ADER.compute_c(4,1,1,3,true) == -1/12
    @test Network_ADER.compute_c(4,2,2,1,true) == -1/2
end

@testset "Computation of C matrices" begin
    # Testing right boundary
    # Order 3, r=0 test
    @test isapprox(Network_ADER.compute_C(3,0) , [1/3 5/6 -1/6;-1 1 0;1 -2 1])
  
    # Testing left boundary
    # Order 3, r=2 test
    @test isapprox(Network_ADER.compute_C(3,2,true) , [-1/6 5/6 1/3;0 -1 1;1 -2 1])
end

@testset "Computation of B matrices" begin
    # Testing right boundary
    # Order 2, r=0 test
    @test isapprox(Network_ADER.compute_B(0,[3/2 -1/2;-1 1]) , [1 -1; -1 1])
    # Order 3, r=2 test
    @test isapprox(Network_ADER.compute_B(2,[1/3 -7/6 11/6;1 -3 2;1 -2 1]) , [4/3 -19/6 11/6; -19/6 25/3 -31/6;11/6 -31/6 10/3])
end

@testset "Computation of protoype D coeffiencents" begin
    C= Network_ADER.compute_C_constants(5)
    @test isapprox(Network_ADER.compute_D(2,C) , [2/3; 1/3])
    @test isapprox(Network_ADER.compute_D(3,C) , [3/10; 3/5 ; 1/10])
end

@testset "Computation of a complete protoype D coeffiencent set" begin
    C= Network_ADER.compute_C_constants(5)
    D_expected = Array{Array{Float64,2},1}(undef,2);
    D_expected[1]=zeros(1,1)
    D_expected[1][1]=1
    D_expected[2]=zeros(2,1)
    D_expected[2][:]=[2/3;1/3]
    @test isapprox(Network_ADER.compute_D_constants(2,C) , D_expected)
end

@testset "Computation of a complete B Matrix set" begin
    C= Network_ADER.compute_C_constants(5)
    @test isapprox(Network_ADER.compute_B_constants(3,C) , D_expected)
end

K=Network_ADER.PR_constants(2)