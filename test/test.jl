using FractionalDiffEq
using Test

@testset "Test Diethelm PECE algorithms" begin
    fun(x, y) = 1-y
    prob=FODEProblem(fun, 1.8, 0, 5, 0.01)
    result=solve(prob, PECE())
    tspan=collect(0:0.01:5)

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

    equation = D(500, 0.5, 0.01) + D(500, 0, 0.01)

    result = solve(equation, 1, 1, 0.01, 5, FODEMatrixDiscrete())

    @test isapprox(result, target; atol=1)

    ########### Yet another test case ##########
    yafun(x, y) = 1-y

    #Analytical solution
    yatarget = []

    #MittagLeffler.jl doesn't support array argument
    for i in 0.01:0.01:20
        push!(yatarget, i^1.8*mittleff(1.8,2.8,-i^1.8))
    end

    yaequation = D(2000, 1.8, 0.01) + D(2000, 0, 0.01)

    yaresult = solve(yaequation, 1, 2, 0.01, 20, FODEMatrixDiscrete())

    @test isapprox(yaresult, yatarget; atol=1)

end