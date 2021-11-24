# ODE Example

It is noteworthy that some differential equation solvers in FractionalDiffEq.jl is also capable of solving **Ordinary Differential Equations**, let's directly see an example here!!

If the IVP is:

```math
\frac{d^2y}{dx^2}+\frac{dy}{dx}=\sin(x)
```
```math
y(0)=0
```

We already know the anlytical solution is

```math
\frac{1}{2}(-e^{-x}-\sin(x)-\cos(x)+2)
```

We can use the **FODEMatrixDiscrete** algorithm to solve this ODE:

```julia
using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ODE\\ Example\$"

T = 30
h=0.05
tspan = collect(0.05:h:T)

f(x)=1/2*(-exp(-x)-sin(x)-cos(x)+2)
target=f.(tspan)

eq = D(600, 2, h)+D(600, 1, h)
rightfun(x) = sin(x)
result = solve(eq, rightfun, 2, h, T, FODEMatrixDiscrete())

plot(tspan, result, title=s, legend=:bottomright, label="ODE Numerical Solution!")

plot!(tspan, target, lw=3,ls=:dash,label="ODE Analytical Solution!")
```

And by plotting the numerical and analytical solution, we can see the matrix discrete algorithm in FractionalDiffEq.jl is quite powerful!

![ODE Example](../assets/ode_example.png)


!!! tip "Better Choice"
    While the solver in FractionalDiffEq.jl can solve ordinary differential equations, we still strongly recommend users to use [SciML/OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) to solve ODEs instead, for various, robust and perfornant algorithms