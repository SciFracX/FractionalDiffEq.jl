@testset "Test Fractional Differences Equations PECE method" begin
    differencefun(x) = 0.5*x+1
    α=0.5;x0=1;
    T=1; h=0.1
    prob = FractionalDifferenceProblem(differencefun, α, x0)
    sol = solve(prob, T, h, PECEDifference())
    @test isapprox(sol.u, [1.0, 2.5, 4.25, 5.125, 5.5625, 5.78125, 5.890625, 5.9453125, 5.97265625, 5.986328125, 5.9931640625]; atol=1e-3)
end