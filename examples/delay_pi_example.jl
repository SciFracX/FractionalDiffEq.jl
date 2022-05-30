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
t0=0; T=2
tspan = (t0, T); h=0.5
prob = FDDEProblem(f, ϕ, 0.97, 0.8, tspan)

result = solve(prob, h, DelayPI())
using Plots
tspan = collect(0:h:T)
plot(tspan, result[:])