# Fractional Order Chen System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.9, 0.9, 0.9]
x0 = [-9, -5, 14]
tf=100
function f(t, x, y, z, k)
    a, b, c, d = 35, 3, 28, -7
    if k == 1
        return a*(y-x)
    elseif k == 2
        return d*x-x*z+c*y
    elseif k == 3
        return x*y-b*z
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Chen System")

```

![Chen](./assets/ChenCover.png)