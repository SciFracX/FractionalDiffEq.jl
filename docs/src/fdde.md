# Fractional Order Delayed Differential Equations

In real world systems, delay is very often encountered in many practical systems, such as automatic control, biology and hydraulic networks, economics and long transmission lines. The delayed differential equation is used to describe these dynamical systems. Fractional order delayed differential equations as the generalization of the delayed differential equations, provide more freedom when we describing these systems, let's see how we can use FractionalDiffEq.jl to accelerate the simulation of delayed differential equations.

The delayed fractional differential equations has the general form:

```math
D^\alpha_ty(t)=f(t,\ y(t),\ y(t-\tau)),\quad t\geq\xi
```

```math
y(t)=\phi(t),\quad t\in[\xi-\tau,\ \xi]
```

While only given the initial condition is not enough to solve the delayed differential equations, a history function ``\phi(t)`` must be provided to describe the history of the system(``\phi(t)`` should be a continuous function).

All we need to do is to pass the function ``f(t,\ y(t),\ y(t-tau))``, and history function ``\phi(t)`` to the ```FDDEProblem``` and choose an algorithm to solve problem:

```julia
using FractionalDiffEq

function ϕ(x)
    if x == 0
        return 19.00001
    else
        return 19.0
    end
end

function f(t, y, ϕ)
    return 3.5*y*(1-ϕ/19)
end

h = 0.05
α = 0.97
τ = 0.8
T = 56
fddeprob = FDDEProblem(f, ϕ, α, τ)
V, y = solve(fddeprob, T, h, DelayPECE())

using Plots
plot(y, V, xlabel="y(t)", ylabel="y(t-τ)")
```

![Delayed](./assets/fdde_example.png)