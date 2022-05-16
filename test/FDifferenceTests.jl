@testset "Test Fractional Differences Equations PECE method" begin
    differencefun(x) = 0.5*x+1
    α=0.5;x0=1;
    T=1; h=0.1
    prob = FractionalDifferenceProblem(differencefun, α, x0)
    sol = solve(prob, T, h, PECEDifference())
    @test isapprox(sol.u, [1.0, 2.5, 4.25, 5.125, 5.5625, 5.78125, 5.890625, 5.9453125, 5.97265625, 5.986328125, 5.9931640625]; atol=1e-3)
end

@testset "Test Grunwald Letnikov method for fractional difference equatoins" begin
    function sys!(du, u, p, t)
        du[1] = -0.05*u[2] - 0.05*u[3] + 0.01*tanh(u[2])
        du[2] = 0.05*u[1] + 0.02*u[2] + 0.01*tanh(u[1])
        du[3] = 0.1 - 0.2*u[3] + 0.05*u[1]*u[3] + 0.01*tanh(u[3])
    end
    
    prob = FractionalDifferenceSystem(sys!, 0.98, [1, -1, 0])
    result = solve(prob, 7, GL())

    @test isapprox(result, [1.0  0.0423841  0.0350316  0.0241793  0.0103323  -0.00594294  -0.0241551
    -1.0  0.0376159  0.0401587  0.0426291  0.0445987   0.0458335    0.0461619
     0.0  0.1        0.179209   0.24285    0.294188    0.335624     0.369018]; atol=1e-5)
end