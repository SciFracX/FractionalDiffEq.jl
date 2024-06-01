@testset "Test Distributed Order Relaxation" begin
    h = 0.5
    t = collect(0:h:1)
    fun(t) = 0
    prob = DODEProblem([1, 0.1], [x -> 6 * x * (1 - x), 0], (0, 1), fun, 1, t)
    result = solve(prob, h, DOMatrixDiscrete())
    @test isapprox(
        result.u, [0.934688468068206, 0.9020578767132177, 0.879670766076252]; atol = 1e-3)
end

@testset "Test Distributed Order Oscillator" begin
    h = 0.5
    t = collect(0:h:1)
    dooscillator(t) = 0 ≤ t ≤ 1 ? (return 8) : (return 0)
    prob = DODEProblem(
        [1, 1, 1], [2, x -> 6 * x * (1 - x), 0], (0, 1), dooscillator, [0; 0], t)
    sol = solve(prob, h, DOMatrixDiscrete())
    @test isapprox(sol.u, [0.0, 0.0, 1.243950673139022]; atol = 1e-3)
end
