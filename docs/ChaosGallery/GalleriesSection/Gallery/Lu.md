# Fractional Order Lu System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.985, 0.99, 0.98]
x0 = [0.2, 0.5, 0.3]
tf=60
function f(t, x, y, z, k)
    a, b, c = 36, 3, 20
    if k == 1
        return a*(y-x)
    elseif k == 2
        return -x*z+c*y
    elseif k == 3
        return x*y-b*z
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Lu System")
```

![Lu](./assets/Lu.png)