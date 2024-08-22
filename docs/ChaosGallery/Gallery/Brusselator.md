# Fractional Order Brusselator System

```math
\begin{cases}
{_{t_0}D_t^\alpha}y_1(t)=a-(\mu+1)y_1(t)+(y_1(t))^2y_2(t)\\
{_{t_0}D_t^\alpha}y_2(t)=\mu y_1(t)-(y_1(t))^2y_2(t)\\
y_1(t_0)=y_{1,0}\\
y_2(t_0)=y_{2,0}
\end{cases}
```

```julia
using FractionalDiffEq, Plots

α = [0.8, 0.8]
u0 = [0.2, 0.03]
h = 0.001
tspan = (0, 100)

function Brusselator!(du, u, p, t)
    a, μ = 1, 4
    du[1] = a-(μ+1)*u[1]+(u[1])^2*u[2]
    du[2] = μ*u[1]-(u[1])^2*u[2]
end

prob = FODEProblem(Brusselator!, α, u0, tspan)
sol = solve(prob, h, GL())

# Phase plane
plot(sol, vars=(1,2))

# Time plane
plot(sol, vars=(0,1,2))
```

![BrusselatorPhase](./assets/Brusselator.png)

![BrusselatorTime](./assets/BrusselatorTime.png)
