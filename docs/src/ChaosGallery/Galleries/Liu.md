# Fractional Order Liu System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [0.2, 0, 0.5]
tf=100
function f(t, x, y, z, k)
    a, b, c, e, n, m = 1, 2.5, 5, 1, 4, 4
    if k == 1
        return -a*x-e*y^2
    elseif k == 2
        return b*y-n*x*z
    elseif k == 3
        return -c*z+m*x*y
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Liu System")
```

![Liu](./assets/Liu.png)