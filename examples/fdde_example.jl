function ϕ(x)
    if x == 0
        return 19.00001
    else
        return 19.0
    end
end

function f(t, y, v)
    return 3.5*y*(1-v/19)
end

h=0.05
α=0.97
τ=0.8
T=56
fddeprob = FDDEProblem(f, α, τ)
V, y = solve(fddeprob, T, h, DelayPECE())

using Plots, SpecialFunctions
plot(V, y)