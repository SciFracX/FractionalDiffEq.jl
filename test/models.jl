# This test file test detailed models in FractionalDiffEq.jl

using FractionalDiffEq
using Test

@testset "Test Bagley Torvik equation" begin
    @test isapprox(bagleytorvik(1, 1, 1, 1, 1, 0.5), [0; 0; 0.12773958089728293]; atol=1e-2)
end

@testset "Test Fractional Lorenz system" begin
    a=40; b=3; c=10; d=15
    q=0.96
    h=0.5
    T=1
    x0=[1, 2, 3]
    prob = FractionalLorenz(a, b, c, d, q, x0)
    x, y, z = solve(prob, h, T, LorenzADM())
    @test isapprox(x, [1.0, 25685.227909653906, -8.14836468585539e22])
    @test isapprox(y, [2.0, 2778.844648448374, -3.1218405233973532e25])
    @test isapprox(z, [3.0, 4518.851378723348, -5.266480656610201e25])
end