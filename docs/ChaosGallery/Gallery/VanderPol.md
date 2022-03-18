# Fractional Order Van der Pol Oscillator

```julia
using FractionalDiffEq

h=0.005
alpha = [1.2, 0.8]
x0 = [0.2, -0.2]
tf=60
function f(t, x, y, k)
    ϵ = 1
    if k == 1
        return y
    elseif k == 2
        return -x-ϵ*(x^2-1)*y
    end
end
prob = FODESystem(f, alpha, x0)
result = solve(prob, h, tf, GLWithMemory())

using Plots
plot(result[:, 1], result[:, 2], title="Fractional Order Van der Pol Oscillator")
```

![VanderPol](./assets/VanderPol.png)