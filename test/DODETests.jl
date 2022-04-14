@testset "Test Distributed Order Matrix method" begin    
    h = 0.5; t = collect(h:h:1);
    fun(t)=-0.1
    prob = DODEProblem([1, 0.1], [x->6*x*(1-x), 0], [0, 1], fun, t)
    result = solve(prob, h, DOMatrixDiscrete())
    @test isapprox(result, [0; -0.06531153193179398]; atol=1e-3)
end