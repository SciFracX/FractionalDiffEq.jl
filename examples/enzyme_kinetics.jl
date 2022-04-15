using FractionalDiffEq
function EnzymeKinetics(t, ϕ, y, k)
    if k == 1
        return 10.5-y[1]/(1+0.0005*ϕ[4]^3)
    elseif k == 2
        return y[1]/(1+0.0005*ϕ[4]^3)-y[2]
    elseif k == 3
        return y[2]-y[3]
    elseif k == 4
        return y[3]-0.5*y[4]
    end
end

q = [60, 10, 10, 20]; α = [0.95, 0.95, 0.95, 0.95]
prob = FDDESystem(f, q, α, 4, 6)
sold, sol = solve(prob, 0.1, DelayABM())
tspan = collect(0:0.1:2)
using Plots
plot(tspan, sol[:, 1])
plot!(tspan, sol[:, 2])
plot!(tspan, sol[:, 3])
plot!(tspan, sol[:, 4])