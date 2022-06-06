# Fractional Order Labyrinth System

```math
\begin{cases}
{_{t_0}^{CF}D_t^\alpha}x(t)=-\sin(y)-bx\\
{_{t_0}^{CF}D_t^\alpha}y(t)=-\sin(z)-by\\
{_{t_0}^{CF}D_t^\alpha}z(t)=-\sin(x)-bz\\
\end{cases}
```

```julia
using FractionalDiffEq

t0=0; tfinal=200; h=0.01;
α = [0.98, 0.98, 0.98]
u0 = [-1, 1, 1]
function fun(du, u, p, t)
    b=0.1;
    du[1] = -sin(u[2])-b*u[1]
    du[2] = -sin(u[3])-b*u[2]
    du[3] = -sin(u[1])-b*u[3]
end
prob = FODESystem(fun, α, u0, (t0, tfinal))
sol = solve(prob, h, NewtonPolynomial())

using Plots
plot3d(sol[1, :], sol[2, :], sol[3, :], title="Fractional Order Labyrinth System")
```

![Labyrinth](./assets/Labyrinth.png)