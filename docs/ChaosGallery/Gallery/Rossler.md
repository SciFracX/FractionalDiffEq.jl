# Fractional Order Rossler System

```julia
using FractionalDiffEq, Plots

h=0.005
alpha = [0.9, 0.85, 0.95]
x0 = [0.5, 1.5, 0.1]
tspan = (0, 120)
function Rossler!(du, u, p, t)
    a, b, c = 0.5, 0.2, 10
    du[1] = -u[2]-u[3]
    du[2] = u[1]+a*u[2]
    du[3] = b+u[3]*(u[1]-c)
end
prob = FODESystem(Rossler!, alpha, x0, tspan)
sol = solve(prob, h, GL())

plot(sol, vars=(1,2,3), title="Fractional Order Rossler System")

```

![Rossler](./assets/Rossler.png)