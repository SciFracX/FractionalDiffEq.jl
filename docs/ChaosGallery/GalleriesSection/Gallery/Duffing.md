# Fractional Order Duffing System

```julia
using FractionalDiffEq

h=0.005
alpha = [0.9, 1]
x0 = [0.21, 0.31]
tf=100
function f(t, x, y, k)
    α, δ, ω = 0.15, 0.3, 1
    if k == 1
        return y
    elseif k == 2
        return x-α*y-x^3+δ*cos(ω*t)
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot(result[:, 1], result[:, 2], title="Fractional Order Duffing System")
```

![Duffing](./assets/Duffing.png)