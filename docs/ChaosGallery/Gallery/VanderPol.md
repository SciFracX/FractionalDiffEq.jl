# Fractional Order Van der Pol Oscillator

```julia
using FractionalDiffEq

h=0.005
alpha = [1.2, 0.8]
x0 = [0.2, -0.2]
tf=60
function VanderPol!(du, u, p, t)
    ϵ = 1
    du[1] = u[2]
    du[2] = -u[1]-ϵ*(u[1]^2-1)*u[2]
end
prob = FODESystem(VanderPol!, alpha, x0, tf)
result = solve(prob, h, GLWithMemory())

using Plots
plot(result[:, 1], result[:, 2], title="Fractional Order Van der Pol Oscillator")
```

![VanderPol](./assets/VanderPol.png)