# Fractional Order Genesio-Tesi System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [-0.1, 0.5, 0.2]
tf=100
function GenesioTesi!(du, u, p, t)
    b1, b2, b3, b4 = 1, 1.1, 0.4, 1.0
    du[1] = u[2]
    du[2] = u[3]
    du[3] -b1*u[1]-b2*u[2]-b3*u[3]+b4*u[1]^2
end
prob = FODESystem(GenesioTesi!, alpha, x0, tf)
result = solve(prob, h, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Genesio-Tesi System")
```

![Genesio-Tesi](./assets/Genesio-Tesi.png)