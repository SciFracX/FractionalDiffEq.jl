# This test file test the auxillary functions in FractionalDiffEq.jl

using FractionalDiffEq
using Test

@testset "Test Riesz Matrix" begin
    @test isapprox(RieszMatrix(0.5, 3, 0.01), [-5.0 4.375 -0.3125; 4.375 -5.0 4.375; -0.3125 4.375 -5.0]; atol=0.00001)
end

@testset "Test ω function" begin
    @test isapprox(omega(3, 0.5), [1.0; -0.5; -0.125; -0.0625]; atol=1e-5)
    @test isapprox(omega(5, 0.5), [1.0; -0.5; -0.125; -0.0625; -0.0390625; -0.02734375]; atol=1e-5)
end

@testset "Test Mittag Leffler function" begin
    @test isapprox(mittleff(1, 2, 1), exp(1)-1; atol=1e-5)
    @test isapprox(mittleff(2, 2, 1), sinh(1); atol=1e-5)
end