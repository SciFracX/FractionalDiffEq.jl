using FractionalDiffEq
function EnzymeKinetics!(dy, y, ϕ, t)
    dy[1] = 10.5-y[1]/(1+0.0005*ϕ[4]^3)
    dy[2] = y[1]/(1+0.0005*ϕ[4]^3)-y[2]
    dy[3] = y[2]-y[3]
    dy[4] = y[3]-0.5*y[4]
end
q = [60, 10, 10, 20]; α = [0.95, 0.95, 0.95, 0.95]
prob = FDDESystem(EnzymeKinetics!, q, α, 4, 6)
sold, sol = solve(prob, 0.1, DelayABM())
tspan = collect(0:0.1:2)
using Plots
plot(tspan, sol[:, 1])
plot!(tspan, sol[:, 2])
plot!(tspan, sol[:, 3])
plot!(tspan, sol[:, 4])