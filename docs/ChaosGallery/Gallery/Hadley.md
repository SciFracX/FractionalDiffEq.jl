# Fractional Order Hadley System

```math
\begin{cases}
{_{t_0}^{ABC}D_t^\alpha}x(t)=-y^2-z^2-ax+a\zeta\\
{_{t_0}^{ABC}D_t^\alpha}y(t)=xy-bxz-y+\delta\\
{_{t_0}^{ABC}D_t^\alpha}z(t)=bxy+xz-z\\
\end{cases}
```

Fractional differential equations with Atangana-Baleanu fractional operator ``{_{t_0}^{ABC}D_t^\alpha}`` in the sense of Caputo.

```julia
using FractionalDiffEq

t0=0;tfinal=50;h=0.01
α = [0.99, 0.99, 0.99]
u0 = [-0.1; 0.1; -0.1]
function fun(du, u, p, t)
    a=0.2;b=4;ζ=8;δ=1;
    du[1] = -u[2]^2-u[3]^2-a*u[1]+a*ζ
    du[2] = u[1]*u[2]-b*u[1]*u[3]-u[2]+δ
    du[3] = b*u[1]*u[2]+u[1]*u[3]-u[3]
end
prob = FODESystem(fun, α, u0, (t0, tfinal))
sol = solve(prob, h, AtanganaSedaAB())

using Plots
plot3d(sol[1, :], sol[2, :], sol[3, :], title="Fractional Order Hadley System")
```

![Hadley](./assets/Hadley.png)