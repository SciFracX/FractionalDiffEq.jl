using FractionalDiffEq

function ϕ(x)
    if x == 0
        return 19.00001
    else
        return 19.0
    end
end

f(t, y, ϕ) = 3.5 * y * (1 - ϕ / 19)

h = 0.05
α = 0.97
τ = 0.8
T = 56
fddeprob = FDDEProblem(f, ϕ, α, τ, T)
V, y = solve(fddeprob, h, DelayPECE())

using Plots
plot(y, V, xlabel = "y(t)", ylabel = "y(t-τ)")
