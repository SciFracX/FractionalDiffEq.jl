@testset "Test Fractional Integral Equations SpectralUltraspherical method" begin
    tspan = LinRange(-1, 1, 5)
    prob = FIEProblem([1, 1], [1, 0.5], 1, tspan)
    sol = solve(prob, 20, SpectralUltraspherical())

    @test isapprox(sol.u, [ 1.0
    0.5231565837302468
    0.42758357615580733
    0.37316567427801584
    0.3362040024463395]; atol=1e-3)
end