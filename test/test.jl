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