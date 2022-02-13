# This test file test detailed models in FractionalDiffEq.jl

using FractionalDiffEq
using Test

@testset "Test Bagley Torvik equation" begin
    @test isapprox(bagleytorvik(1, 1, 1, 1, 1, 0.5), [0; 0; 0.12773958089728293]; atol=1e-2)
end