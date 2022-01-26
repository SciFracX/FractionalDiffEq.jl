using FractionalDiffEq

function ϕ(x)
    if x == 0
        return 19.00001
    else
        return 19.0
    end
end

function f(t, y, ϕ)
    return 3.5*y*(1-ϕ/19)
end

h = 0.05
α = 0.97
τ = 0.8
T = 56
fddeprob = FDDEProblem(f, ϕ, α, τ)
V, y = solve(fddeprob, T, h, DelayPECE())

using Plots
plot(V, y)