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

@testset "Test singleterm show method" begin
    # SingleTermFODEProblem
    singlefun(x, y) = 1-y
    singletermprob = SingleTermFODEProblem(singlefun, 1.8, 0, 5)
    @test_nowarn show(singletermprob)
end

@testset "Test Multiterms show method" begin
    # MultiTermsFODEProblem
    rightfun(x, y) = 172/125*cos(4/5*x)
    multitermsprob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0 0 0 0 0 0], 0, 10)
    multitermssol = solve(multitermsprob, 0.5, PIEX())

    @test_nowarn show(multitermsprob)
    @test_nowarn show(multitermssol)
end

@testset "Test FODESystem show method" begin
    # FODESystem
    h=0.5; tf=1
    alpha = [0.99, 0.99, 0.99]
    x0 = [1, 0, 1]
    function testf!(du, u, p, t)
        a, b, c = 10, 28, 8/3
        du[1] = a*(u[2]-u[1])
        du[2] = u[1]*(b-u[3])-u[2]
        du[3] = u[1]*u[2]-c*u[3]
    end
    fodesystemprob = FODESystem(testf!, alpha, x0, tf)
    fodesystemsol = solve(fodesystemprob, h, GL())

    @test_nowarn show(fodesystemprob)
    @test_nowarn show(fodesystemsol)
end

@testset "Test FDDEProblem show method" begin
    # FDDEProblem
    function ϕ(x)
        if x == 0
            return 19.00001
        else
            return 19.0
        end
    end
    function delayf(t, y, ϕ)
        return 3.5*y*(1-ϕ/19)
    end
    fddeprob = FDDEProblem(delayf, ϕ, 0.97, 0.8, 2, 0)
    fddesol = solve(fddeprob, 0.5, DelayPI())

    @test_nowarn show(fddeprob)
    @test_nowarn show(fddesol)
end

@testset "Test FDDESystem show method" begin
    function EnzymeKinetics!(dy, y, ϕ, t)
        dy[1] = 10.5-y[1]/(1+0.0005*ϕ[4]^3)
        dy[2] = y[1]/(1+0.0005*ϕ[4]^3)-y[2]
        dy[3] = y[2]-y[3]
        dy[4] = y[3]-0.5*y[4]
    end
    q = [60, 10, 10, 20]; α = [0.95, 0.95, 0.95, 0.95]
    fddesystemprob = FDDESystem(EnzymeKinetics!, q, α, 4, 6)

    @test_nowarn show(fddesystemprob)
end

@testset "Test DODEProblem show method" begin
    h = 0.5; t = collect(0:h:1)
    dodefun(t)=0
    dodeprob = DODEProblem([1, 0.1], [x->6*x*(1-x), 0], (0, 1), dodefun, 1, t)
    dodesol = solve(dodeprob, h, DOMatrixDiscrete())

    @test_nowarn show(dodeprob)
    @test_nowarn show(dodesol)
end

@testset "Test FractionalDifferenceProblem" begin
    differencefun(x) = 0.5*x+1
    α=0.5;x0=1;
    T=1; h=0.1
    differenceprob = FractionalDifferenceProblem(differencefun, α, x0)
    differencesol = solve(differenceprob, T, h, PECEDifference())

    @test_nowarn show(differenceprob)
    @test_nowarn show(differencesol)
end

@testset "Test FIEProblem show method" begin
    fietspan = LinRange(-1, 1, 5)
    fieprob = FIEProblem([1, 1], [1, 0.5], 1, fietspan)
    fiesol = solve(fieprob, 20, SpectralUltraspherical())

    @test_nowarn show(fieprob)
    @test_nowarn show(fiesol)
end