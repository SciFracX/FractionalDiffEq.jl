@testset "Test FdeSolver wrapper" begin
    using FractionalDiffEq, FdeSolver, Test

    function chua!(du, x, p, t)
        a = 10.725
        b = 10.593
        c = 0.268
        m0 = -1.1726
        m1 = -0.7872
        du[1] = a * (x[2] - x[1] -
                 (m1 * x[1] + 0.5 * (m0 - m1) * (abs(x[1] + 1) - abs(x[1] - 1))))
        du[2] = x[1] - x[2] + x[3]
        du[3] = -b * x[2] - c * x[3]
    end
    α = [0.93, 0.99, 0.92]
    u0 = [0.2; -0.1; 0.1]
    tspan = (0, 0.5)
    prob = FODEProblem(chua!, α, u0, tspan)
    @test_no_warn solve(prob, FdeSolverPECE(), dt = 0.01)
end