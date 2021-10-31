using FractionalDiffEq
using Test
using MittagLeffler


@testset "Test Diethelm PECE algorithms" begin
    fun(x, y) = 1-y
    prob=FDEProblem(fun, 1.8, 0, 5, 0.01)
    result=solve(prob, PECE())
    tspan=collect(0:0.01:5)

    target = []

    for i in 0:0.01:5
        push!(target, i^1.8*mittleff(1.8, 2.8,-i^1.8))
    end

    @test isapprox(result, target; atol=1e-5)
end