# Fractional Order Arneodo System

```julia
using FractionalDiffEq, Plots

h=0.005
alpha = [0.97, 0.97, 0.96]
x0 = [-0.2, 0.5, 0.2]
tspan = (0, 100)
function Arneodo!(du, u, p, t)
    b1, b2, b3, b4 = -5.5, 3.5, 0.8, -1.0
    du[1] = u[2]
    du[2] = u[3]
    du[3] = -b1*u[1]-b2*u[2]-b3*u[3]+b4*u[1]^3
end
prob = FODESystem(Arneodo!, alpha, x0, tspan)
sol = solve(prob, h, GL())

plot(sol, vars=(1,2,3), title="Fractional Order Arneodo System")
```

![Arneodo](./assets/Arneodo.png)