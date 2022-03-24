# Fractional Order Chen System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.9, 0.9, 0.9]
x0 = [-9, -5, 14]
tf=100
function Chen!(du, u, p, t)
    a, b, c, d = 35, 3, 28, -7
    du[1] = a*(u[2]-u[1])
    du[2] = d*u[1]-u[1]*u[3]+c*u[2]
    du[3] = u[1]*u[2]-b*u[3]
end
prob = FODESystem(Chen!, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Chen System")
```

![Chen](./assets/Chen.png)