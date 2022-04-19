# Fractional Order Liu System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [0.2, 0, 0.5]
tf=100
function Liu!(du, u, p, t)
    a, b, c, e, n, m = 1, 2.5, 5, 1, 4, 4
    du[1] = -a*u[1]-e*u[2]^2
    du[2] = b*u[2]-n*u[1]*u[3]
    du[3] = -c*u[3]+m*u[1]*u[2]
end
prob = FODESystem(Liu!, alpha, x0, tf)
result = solve(prob, h, GL())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Liu System")
```

![Liu](./assets/Liu.png)