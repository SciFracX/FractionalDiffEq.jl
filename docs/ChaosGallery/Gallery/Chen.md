# Fractional Order Chen System

Fractional order[Chen system](https://www.worldscientific.com/doi/abs/10.1142/s0218127499001024):

```math
\begin{cases}
D^{\alpha_1}x=a(y-x)\\
D^{\alpha_2}y=(c-a)x-xz+cy\\
D^{\alpha_3}z=xy-bz
\end{cases}
```

Behave chaotic property when ``a=35``, ``b=3``, ``c=28``.


```julia
using FractionalDiffEq, Plots

h=0.005
alpha = [0.9, 0.9, 0.9]
x0 = [-9, -5, 14]
tspan = (0, 100)
function Chen!(du, u, p, t)
    a, b, c, d = 35, 3, 28, -7
    du[1] = a*(u[2]-u[1])
    du[2] = d*u[1]-u[1]*u[3]+c*u[2]
    du[3] = u[1]*u[2]-b*u[3]
end
prob = FODESystem(Chen!, alpha, x0, tspan)
sol = solve(prob, h, GL())

plot(sol, vars=(1,2,3), title="Fractional Order Chen System")
```

![Chen](./assets/Chen.png)