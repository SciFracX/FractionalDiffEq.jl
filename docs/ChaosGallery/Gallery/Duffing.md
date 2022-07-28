# Fractional Order Duffing System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.9, 1]
x0 = [0.21, 0.31]
tapan = (0, 100)
function Duffing!(du, u, p, t)
    α, δ, ω = 0.15, 0.3, 1
    du[1] = u[2]
    du[2] = u[1]-α*u[2]-u[1]^3+δ*cos(ω*t)
end
prob = FODESystem(Duffing!, alpha, x0, tspan)
result = solve(prob, h, GL())

using Plots
plot(result, vars=(1,2), title="Fractional Order Duffing System")
```

![Duffing](./assets/Duffing.png)