# Fractional Order Arneodo System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.97, 0.97, 0.96]
x0 = [-0.2, 0.5, 0.2]
tf=100
function f(t, x, y, z, k)
    b1, b2, b3, b4 = -5.5, 3.5, 0.8, -1.0
    if k == 1
        return y
    elseif k == 2
        return z
    elseif k == 3
        return -b1*x-b2*y-b3*z+b4*x^3
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Arneodo System")
```

![Arneodo](./assets/Arneodo.png)