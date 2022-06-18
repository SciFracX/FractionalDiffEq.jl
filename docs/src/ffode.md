# Fractal-Fractional Order Ordinary Differential Equations

The concept of fractal–fractional differentiation has appeared as a combination of two
mathematical concepts: fractal differentiation and fractional differentiation. These operators are helpful for modeling more complex problems, also
they allow us to understand many physical problems which have fractal properties.

## Fractal-fractional ODE system

To construct fractal-fractional differential problems, we need to use ```FFODESystem``` to define our problem:

```julia-repl
FFODESystem(f, [α, β], u0, tspan)
```

Let's see the Lorenz system in Atangana-Baleanu-Caputo sense:

$$
\begin{cases}
{^{FFM}D^{\alpha,\beta}}x=a(y-x)\\
{^{FFM}D^{\alpha,\beta}}y=(b-z)x-y\\
{^{FFM}D^{\alpha,\beta}}z=xy-cz
\end{cases}
$$

```julia
using FractionalDiffEq, Plots
α=1;β=1;h=0.01;tfinal=50
u0 = [-2, 1, -1]
function fun(du, u, p, t)
    a=10;b=28;c=8/3
    du[1] = a*(u[2]-u[1])
    du[2] = (b-u[3])*u[1]-u[2]
    du[3] = u[1]*u[2]-c*u[3]
end
prob = FFODESystem(fun, [α, β], u0, (0, tfinal))
sol = solve(prob, h, AtanganaSeda())
plot3d(sol[1, :], sol[2, :], sol[3, :], title="Fractal-fractional Order Lorenz System")
```

![ABCLorenz](./assets/LorenzABC.png)