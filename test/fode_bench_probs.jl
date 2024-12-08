# from "BENCHMARK PROBLEMS FOR CAPUTO  FRACTIONAL-ORDER ORDINARY  DIFFERENTIAL EQUATIONS"
using SpecialFunctions
@testset "1st problem" begin
    function f!(u, p, t)
        return if 0 ≤ t ≤ 1
            1/gamma(1.3)*t^0.3
        else
            1/gamma(1.3)*t^0.3 - 2/gamma(2.3)*(t-1)^1.3
        end
    end
    analytical_solution(u, p, t) = if 0 ≤ t ≤ 1
        return t
    else
        return t - (t-1)^2
    end
    u0 = 0.0
    tspan = (0.0, 2.0)
    order = 0.7
    fun = ODEFunction(f!, analytic = analytical_solution)
    prob = FODEProblem(fun, order, u0, tspan)
    sol = solve(prob, BDF(), dt=0.05)
end

@testset "5th problem" begin
    function f!(du, u, p, t)
        du[1] = 1/sqrt(pi)*(((u[2]-0.5)*(u[3]-0.3))^(1/6) + sqrt(t))
        du[2] = gamma(2.2)*(u[1]-1)
        du[3] = gamma(2.8)/gamma(2.2)*(u[2]-0.5)
    end
    analytical_solution(u, p, t) = [t+1, t^1.2+0.5, t^1.8+0.3]
    u0 = [1.0, 0.5, 0.3]
    tspan = (0.0, 5.0)
    order = [0.5, 0.2, 0.6]
    fun = ODEFunction(f!, analytic = analytical_solution)
    prob = FODEProblem(fun, order, u0, tspan)
    sol = solve(prob, BDF(), dt=0.05)
end

