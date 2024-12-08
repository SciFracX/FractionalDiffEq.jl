module FractionalDiffEqFdeSolverExt

using DiffEqBase, SciMLBase, FractionalDiffEq
using FdeSolver

# const default_values = (2^-6, 1, nothing, 1e-6, 100)
function SciMLBase.__solve(prob::FODEProblem, alg::FdeSolverPECE; dt = 0.0, abstol = 1e-3, maxiters = 1000, kwargs...)
    (; f, order, tspan, u0, p) = prob

    iip = isinplace(prob)
    
    tSpan = [first(tspan), last(tspan)]
    # FdeSolver only supports out-of-place computing
    newf = if iip
        function (t, y, par)
            du = similar(y)
            f(du, y, par, t)
            return du
        end
    else
        function (t, y, par)
            return f.(y, par, t)
        end
    end
    par = p isa SciMLBase.NullParameters ? nothing : p
    length(u0) == 1 && (u0 = first(u0))
    t, y = FDEsolver(newf, tSpan, u0, order, par, JF = prob.f.jac, h = dt, tol = abstol)
    u = collect(Vector{eltype(y)}, eachrow(y))

    return DiffEqBase.build_solution(prob, alg, t, u)
end

export FdeSolverPECE

end
