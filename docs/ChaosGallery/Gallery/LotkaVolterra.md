# Fractional Order Lotka Volterra System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [1, 1.4, 1]
tspan = (0, 50)
function LotkaVolterra!(du, u, p, t)
    a, b, c, d, e, p, s = 1, 1, 1, 1, 2, 3, 2.7
    du[1] = a*u[1]+e*u[1]^2-b*u[1]*u[2]-s*u[3]*u[1]^2
    du[2] = -c*u[2]+d*u[1]*u[2]
    du[3] = -p*u[3]+s*u[3]*u[1]^2
end
prob = FODESystem(LotkaVolterra!, alpha, x0, tspan)
result = solve(prob, h, GL())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Lotka Volterra System")
```

![LotkaVolterra](./assets/LotkaVolterra.png)