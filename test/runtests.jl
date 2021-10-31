using Network_ADER
using Test


@testset "Network_ADER.jl" begin
    # Write your tests here.
    @test Network_ADER.my_f(2,1) == 5
    @test Network_ADER.my_f(3,3) == 9    
end
