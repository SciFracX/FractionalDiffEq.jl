# Multi-term FODE

By specifying different orders in the equation, we can handle multi-terms FODE now!

Let's see if we have an initial value problem with multiple terms derivative containing both fractional and integer, we can use the **FODEMatrixDiscrete** algorithm to solve the equation.

All we need to do is passing the parameters and orders of the fractional ordinary differential equation to the API ```solve``` as two arrays.

!!! warning "The parameters and orders array must have the same length"
    When we are solving the multi-terms FODE problem, please note we should keep the parameters array and orders array have the same length.

## Detailed Usage

Let's see if we have an equation like:

```math
2y''(t)+4D^{1.5}y(t)=1
```

To solve the multi-terms fractional order  equation, you can use the code:

```julia
rightside = 1
prob = MultiTermsFODEProblem([2, 4], [2, 1.5], rightside, [0, 0], 30)
sol = solve(prob, 0.01, FODEMatrixDiscrete())
```

Bingo! The result ```sol``` is the numerical solution of this equation!!!!

## Example

We have an initial problem:

```math
y'''(t)+\frac{1}{16} {^C_0D^{2.5}_ty(t)}+\frac{4}{5}y''(t)+\frac{3}{2}y'(t)+\frac{1}{25}{^C_0D^{0.5}_ty(t)}+\frac{6}{5}y(t)=\frac{172}{125}\cos(\frac{4t}{5})
```

```math
y(0)=0,\ y'(0)=0,\ y''(0)=0
```

Model this problem and solve the equation:

```julia
using FractionalDiffEq, Plots
T = 30; h = 0.01
rightfun(x, y) = 172/125*cos(4/5*x)
prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], T)
sol = solve(prob, h, PIEx())
plot(sol, legend=:bottomright)
```

By solving the equation and plotting the result, we can see its solution here:

![Solution](./assets/complicated_example.png)


!!! info "Example in GitHub"
    This example is an official example in the source code, please see the [example folder](https://github.com/SciFracX/FractionalDiffEq.jl/blob/master/examples/complicated_example.jl)

## Solvers for Ordinary differential equations

FractionalDiffEq.jl is also able to solve ordinary differential equations~ Let's see an example here:

```math
y''(x) + y'(x) = \sin(x)
```

```math
y(0) = 0
```

```julia
using FractionalDiffEq, Plots
T = 30; h = 0.05
tspan = collect(h:h:T)
f(x) = 1/2*(-exp(-x)-sin(x)-cos(x)+2)
target =f.(tspan)
rightfun(x) = sin(x)
prob = MultiTermsFODEProblem([1, 1], [2, 1], rightfun, T)
sol = solve(prob, h, FODEMatrixDiscrete())
plot(sol, title=s, legend=:bottomright, label="ODE Numerical Solution!")
plot!(tspan, target, lw=3,ls=:dash,label="ODE Analytical Solution!")
```

![ODE Example](./assets/ode_example.png)