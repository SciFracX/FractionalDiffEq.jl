module FractionalDiffEqFdeSolverExt

using DiffEqBase, SciMLBase, FractionalDiffEq
using FdeSolver

# const default_values = (2^-6, 1, nothing, 1e-6, 100)
function SciMLBase.__solve(prob::FODEProblem, alg::FdeSolverPECE; dt = 0.0, abstol = 1e-3, maxiters = 1000, kwargs...)
    (; f, order, tspan, u0, p) = prob
    
    tSpan = [first(tspan), last(tspan)]
    # FdeSolver only supports out-of-place computing
    newf = function (t, y, par)
        du = similar(y)
        f(du, y, par, t)
        return du
    end
    t, y = FDEsolver(newf, tSpan, u0, order, p, JF = prob.f.jac, h = dt, tol = abstol)
    u = eachrow(y)

    return DiffEqBase.build_solution(prob, alg, t, u)
end

export FdeSolverPECE

end
