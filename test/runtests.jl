using RealPolynomialRoots
using Test

@testset "RealPolynomialRoots.jl" begin
    ps = [-1, 254, -16129, 0, 0, 0, 0, 1]
    @test length(real_roots_sqfree(ps)) == 3

    # prod( (x-i) for i in 1:10
    ps = [3628800, -10628640, 12753576, -8409500, 3416930, -902055, 157773, -18150, 1320, -55, 1]
    @test length(real_roots_sqfree(ps)) == 10
    
end
