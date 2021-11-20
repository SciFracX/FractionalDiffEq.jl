# This test file test the auxillary functions in FractionalDiffEq.jl

using FractionalDiffEq
using Test
using MittagLeffler


@testset "Test D" begin
    @test isapprox(D(3, 0.5, 0.01), [10.0 0.0 0.0; -5.0 10.0 0.0; -1.25 -5.0 10.0]; atol=0.0001)
    @test isapprox(D(2, 2, 0.01), [10000.0 0.0; -20000.0 10000.0])
end

@testset "Test Riesz Matrix" begin
    @test isapprox(RieszMatrix(0.5, 3, 0.01), [-5.0 4.375 -0.3125; 4.375 -5.0 4.375; -0.3125 4.375 -5.0]; atol=0.00001)
end

