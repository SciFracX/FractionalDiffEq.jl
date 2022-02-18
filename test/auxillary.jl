# This test file test the auxillary functions in FractionalDiffEq.jl

using FractionalDiffEq
using Test

@testset "Test eliminator" begin
    @test isapprox(eliminator(3, 2), [1.0 0.0 0.0; 0.0 0.0 1.0]; atol=1e-5)
end

@testset "Test Riesz Matrix" begin
    @test isapprox(RieszMatrix(0.5, 3, 0.01), [-5.0 4.375 -0.3125; 4.375 -5.0 4.375; -0.3125 4.375 -5.0]; atol=0.00001)
end

@testset "Test ω function" begin
    @test isapprox(omega(3, 0.5), [1.0; -0.5; -0.125; -0.0625]; atol=1e-5)
    @test isapprox(omega(5, 0.5), [1.0; -0.5; -0.125; -0.0625; -0.0390625; -0.02734375]; atol=1e-5)
end

@testset "Test distributed order utils function" begin
    @test isapprox(DOB(x->6*x*(1-x), [0, 1], 0.1, 1, 0.1), [3.547014602981364]; atol=1e-5)
    @test isapprox(DOF(x->6*x*(1-x), [0, 1], 0.1, 1, 0.1), [-3.547014602981364]; atol=1e-5)
    @test isapprox(DORANORT(x->6*x*(1-x), [0, 1], 0.1, 1, 0.1), [2.120468497843435]; atol=1e-5)
end

#FIXME: Add array type tests
@testset "Test Mittag Leffler function" begin
    @test isapprox(mittleff(1, 2, 1), exp(1)-1; atol=1e-5)
    @test isapprox(mittleff(2, 2, 1), sinh(1); atol=1e-5)
    @test isapprox(mittleff(1, 1, 0), 1; atol=1e-5)
    @test isapprox(mittleff(1, 2, [1, 2, 3]), [1.718281828459045; 3.194528049465325; 6.361845641062556]; atol=1e-2)


    @test isapprox(mittlefferr(1, 1, 1, 1), 2.718281828459045; atol=1e-5)
    @test isapprox(mittlefferr(1, 1, 1), 2.718281828459045; atol=1e-5)

    @test isapprox(mittleffderiv(1, 1, 1), 2.718281828459045; atol=1e-5)
    @test isapprox(mittleffderiv(1, 1), 2.718281828459045; atol=1e-5)

    @test isapprox(mittleffderiv(1, 1, 1/3), 1.39561242508609; atol=1e-5)
    @test isapprox(mittleffderiv(2, 2, 1), 0.18393972058572117; atol=1e-5)

    @test isapprox(mittleff(1, 2, 1, 1), 1.7182818284590438; atol=1e-5)
    @test isapprox(mittleff(2, 2, 1, 1), 1.1752011936438016; atol=1e-3)
    @test isapprox(mittleff(2.1, 2, 1, 1), 1.1527981788744044; atol=1e-4)
end

@testset "Test meshgrid" begin
    a, b = meshgrid(0:1, 0:1)
    isapprox(a, [0 1; 0 1]; atol=1e-2)
    isapprox(b, [0 0; 1 1]; atol=1e-2)
end

@testset "Test isFunction" begin
    @test isFunction(x->x)==true
    @test isFunction("Hello")==false
end

@testset "FLMM Weihgts function" begin
    a, b, c = Weights(0.5, 2)
    @test a ≈ [0.816497  0.544331  0.408248]
    @test b ≈ [-0.816497  -0.396913  -0.333742     
    0.0        0.393173   0.370502     
   -0.408248  -0.228708  -0.210067]
   @test c ≈ 2
end

@testset "Test FastConv" begin
    result = FastConv([1 2 3], [1 2 3])
    @test result ≈ [13 13 10]
end