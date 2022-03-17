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
T = 56; h=0.05
prob = FDDEProblem(f, ϕ, 0.97, 0.8, 0, T)

result = solve(prob, h, DelayPI())
using Plots
tspan = collect(0:h:T)
plot(tspan, result[:])