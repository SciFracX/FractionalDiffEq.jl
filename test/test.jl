using FractionalDiffEq
using Test

@testset "Test Diethelm PECE algorithms" begin
    fun(x, y) = 1-y
    prob = SingleTermFODEProblem(fun, 1.8, 0, 5)
    result = solve(prob, 0.01, PECE())
    tspan = collect(0:0.01:5)

    target = []

    for i in 0:0.01:5
        push!(target, i^1.8*mittleff(1.8, 2.8,-i^1.8))
    end

    @test isapprox(result, target; atol=1)
end

@testset "Test Matrix discrete method" begin
    fun(x, y) = 1-y
    tspan=collect(0.01:0.01:5)
    target = []
    for i in 0.01:0.01:5
        push!(target, i^0.5*mittleff(0.5, 1.5,-i^0.5))
    end

    singleprob = MultiTermsFODEProblem([1, 1], [0.5, 0], 1)
    result = solve(singleprob, 0.01, 5, FODEMatrixDiscrete())

    @test isapprox(result, target; atol=1)

    ########### Yet another test case ##########
    yafun(x, y) = 1-y

    #Analytical solution
    yatarget = []

    #MittagLeffler.jl doesn't support array argument
    for i in 0.01:0.01:20
        push!(yatarget, i^1.8*mittleff(1.8,2.8,-i^1.8))
    end

    highsingleprob = MultiTermsFODEProblem([1, 1], [1.8, 0], 1)
    yaresult = solve(highsingleprob, 0.01, 20, FODEMatrixDiscrete())

    @test isapprox(yaresult, yatarget; atol=1)

end

@testset "Test Closed Form method" begin
    t=collect(0:0.002:10);

    rightfun(x) = sin(x)

    prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90], [1 0.3])

    result=solve(prob, t, ClosedForm())

end

@testset "Test GLWithMemory algorithm" begin
    h=0.5
    alpha = [0.99, 0.99, 0.99]
    x0 = [1, 0, 1]
    tf=1
    function f(t, x, y, z, k)
        a, b, c = 10, 28, 8/3
        if k == 1
            return a*(y-x)
        elseif k == 2
            return x*(b-z)-y
        elseif k == 3
            return x*y-c*z
        end
    end
    prob = FODESystem(f, alpha, x0)
    x, y, z = solve(prob, h, tf, GLWithMemory())
    @test x≈[1.0, -4.044777750283594, 84.80744193501619]
    @test y≈[0.0, 13.593899925765704, -51.12509457411144]
    @test z≈[1.0, -0.3526074000756252, -27.554093040332816]
end

@testset "Test FPDEMatrixDiscrete" begin
    @test isapprox(solve(0.5, 0.5, 3, 2, 2, FPDEMatrixDiscrete()), [0 0; 0 0]; atol=1e-2)
end

@testset "Test Bagley Torvik equation" begin
    @test isapprox(bagleytorvik(1, 1, 1, 1, 1, 0.5), [0; 0; 0.12773958089728293]; atol=1e-2)
end