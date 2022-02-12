# Get Start

We would use the simple example -- [Relaxation Oscillation Process](https://encyclopediaofmath.org/wiki/Relaxation_oscillation) example to show you how to use FractionalDiffEq.jl🙂

The mathematical model of the Relaxation Oscillation can be abstracted as IVP:

```math
D^{1.8}y(t)+y(t)=1,\ (t>0)
```

```math
y^{(k)}(0)=0
```

While we can know the analytical solution of this equation is:

```math
u(t)=t^{1.8}E_{1.8,\ 2.8}(-t^{1.8})
```

We can solve this problem by the following code using FractionalDiffEq.jl:

```julia
using FractionalDiffEq
using Plots, LaTeXStrings

# Analytical solution
analytical(x) = x.^1.8 .*mittleff(1.8, 2.8, -x.^1.8)

s="\$D^{1.8}y(x)=1-y(x),\\ y(0)=0\$"

# Numerical solution
fun(x, y) = 1-y
h=0.01;T=20;u0=0
prob = SingleTermFODEProblem(fun, 1.8, u0, T)
result = solve(prob, h, T, PECE())
tspan = collect(0:0.01:20)
target = analytical(tspan)

plot(tspan, result, title=s, linewidth=5, label="Numerical", legend=:bottomright)
plot!(tspan, target, lw=3, ls=:dash, label="Analytical")
```

By plotting the numerical result, we can get the approximation result:

![Relaxation Oscillation](./assets/example.png)

To provide users a simple way to solve fractional differential equations, we follow the design pattern of [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

## Step 1: Defining a Problem

First, we need to specify the problem we want to solve. Just by passing the parameter —— describing function, order, step size and time span:

```julia
using FractionalDiffEq

fun(x, y) = 1-y
prob = SingleTermFODEProblem(fun, 1.8, 0.01, 20)
```

The ```SingleTermFODEProblem``` is a class of fractional differential equation, describing equations with ``D^{\alpha}u=f(t, u)`` pattern. For other patterns of fractional differential equation, please refer to [Problem types](@ref problems)

## Step 2: Solving a Problem

After defining a problem, we can solve it by calling the ```solve``` function:

```julia
result = solve(prob, h, T, Alg())
```

Note that there are different algorithms for differential fractional differential equations, such as FODE, FPDE, FDDE, we need to choose a properiate algorithm for specific problem. For all the algorithms, please refer to [algorithms documentation](@ref algorithms).

## Step3 : Analyzing the Solution

Simply call plot to visualize the solution:

```julia
using Plots
plot(tspan, result)
```

![Relaxation Oscillation](./assets/example.png)