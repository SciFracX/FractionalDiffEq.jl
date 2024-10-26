# This test file test the auxillary functions in FractionalDiffEq.jl

using FractionalDiffEq
using Test

#FIXME: Add array type tests
@testset "Test Mittag Leffler function" begin
    @test isapprox(mittleff(1, 2, 1), exp(1) - 1; atol = 1e-5)
    @test isapprox(mittleff(2, 2, 1), sinh(1); atol = 1e-5)
    @test isapprox(mittleff(1, 1, 0), 1; atol = 1e-5)
    @test isapprox(mittleff(1, 2, [1, 2, 3]),
        [1.718281828459045; 3.194528049465325; 6.361845641062556]; atol = 1e-2)

    @test isapprox(mittleffderiv(1, 1, 1), 2.718281828459045; atol = 1e-5)
    @test isapprox(mittleffderiv(1, 1), 2.718281828459045; atol = 1e-5)

    @test isapprox(mittleffderiv(1, 1, 1 / 3), 1.39561242508609; atol = 1e-5)
    @test isapprox(mittleffderiv(2, 2, 1), 0.18393972058572117; atol = 1e-5)

    @test isapprox(mittleff(1, 2, 1, 1), 1.7182818284590438; atol = 1e-5)
    @test isapprox(mittleff(2, 2, 1, 1), 1.1752011936438016; atol = 1e-3)
    @test isapprox(mittleff(2.1, 2, 1, 1), 1.1527981788744044; atol = 1e-4)

    # Also tests for Complex z:

    m = mittleff(1, 2, im)
    @test isapprox(real(m), 0.8414709848078965; atol = 1e-5)
    @test isapprox(imag(m), 0.4596976941318611; atol = 1e-5)
end

@testset "Test isFunction" begin
    @test isFunction(x -> x) == true
    @test isFunction("Hello") == false
end
