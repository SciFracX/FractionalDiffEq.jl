# Fractional Order Qi System

```julia
using FractionalDiffEq, Plots

function qi(t, x, y, z, k)
    a, b, c, d, r = 35, 8/3, 80, -1, 1
    if k == 1
        return -a*x+a*y+r*y*z
    elseif k == 2
        return c*x+d*y-x*z
    elseif k == 3
        return -b*z+x*y
    end
end

alpha = [0.98, 0.98, 0.98]
h=0.001
T=50
x0=[0.1, 0.2, 0.3]
prob=FODESystem(qi, alpha, x0)
result = solve(prob, h, T, GLWithMemory())

plot3d(result[:, 1], result[:, 2], result[:, 3], title="Fractional Order Qi System")
```

![Rossler](./assets/Qi.png)