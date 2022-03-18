# Fractional Order Genesio-Tesi System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.95, 0.95, 0.95]
x0 = [-0.1, 0.5, 0.2]
tf=100
function f(t, x, y, z, k)
    b1, b2, b3, b4 = 1, 1.1, 0.4, 1.0
    if k == 1
        return y
    elseif k == 2
        return z
    elseif k == 3
        return -b1*x-b2*y-b3*z+b4*x^2
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Genesio-Tesi System")
```

![Genesio-Tesi](./assets/Genesio-Tesi.png)