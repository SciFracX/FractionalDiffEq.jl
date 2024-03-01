# This test file test the auxillary functions in FractionalDiffEq.jl

using FractionalDiffEq
using Test

#FIXME: Add array type tests
@testset "Test Mittag Leffler function" begin
    @test isapprox(mittleff(1, 2, 1), exp(1)-1; atol=1e-5)
    @test isapprox(mittleff(2, 2, 1), sinh(1); atol=1e-5)
    @test isapprox(mittleff(1, 1, 0), 1; atol=1e-5)
    @test isapprox(mittleff(1, 2, [1, 2, 3]), [1.718281828459045; 3.194528049465325; 6.361845641062556]; atol=1e-2)

    @test isapprox(mittleffderiv(1, 1, 1), 2.718281828459045; atol=1e-5)
    @test isapprox(mittleffderiv(1, 1), 2.718281828459045; atol=1e-5)

    @test isapprox(mittleffderiv(1, 1, 1/3), 1.39561242508609; atol=1e-5)
    @test isapprox(mittleffderiv(2, 2, 1), 0.18393972058572117; atol=1e-5)

    @test isapprox(mittleff(1, 2, 1, 1), 1.7182818284590438; atol=1e-5)
    @test isapprox(mittleff(2, 2, 1, 1), 1.1752011936438016; atol=1e-3)
    @test isapprox(mittleff(2.1, 2, 1, 1), 1.1527981788744044; atol=1e-4)

    # Also tests for Complex z:

    m = mittleff(1, 2, im)
    @test isapprox(real(m), 0.8414709848078965; atol=1e-5)
    @test isapprox(imag(m), 0.4596976941318611; atol=1e-5)
end

@testset "Test isFunction" begin
    @test isFunction(x->x)==true
    @test isFunction("Hello")==false
end

@testset "Test DODEProblem show method" begin
    h = 0.5; t = collect(0:h:1)
    dodefun(t)=0
    dodeprob = DODEProblem([1, 0.1], [x->6*x*(1-x), 0], (0, 1), dodefun, 1, t)
    dodesol = solve(dodeprob, h, DOMatrixDiscrete())

    @test_nowarn show(dodeprob)
    @test_nowarn show(dodesol)
end

@testset "Test FractionalDiscreteProblem" begin
    discretefun(x) = 0.5*x+1
    α=0.5;x0=1;
    T=1; h=0.1
    discreteprob = FractionalDiscreteProblem(discretefun, α, x0, T)
    discretesol = solve(discreteprob, h, PECE())

    @test_nowarn show(discreteprob)
    @test_nowarn show(discretesol)
end

@testset "Test FFMODEProblem show method" begin
    α=1;β=1;h=0.1;tfinal=0.5
    u0=[-2, 1, -1]
    function fun(du, u, p, t)
        a=10;b=28;c=8/3
        du[1] = a*(u[2]-u[1])
        du[2] = (b-u[3])*u[1]-u[2]
        du[3] = u[1]*u[2]-c*u[3]
    end
    prob = FFMODEProblem(fun, [α, β], u0, (0, tfinal))

    @test_nowarn show(prob)
end