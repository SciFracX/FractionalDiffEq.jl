# Multi term FDE

By specifying different order in the equation, we can handle multi-term differential equations now!

Let's see if we have an initial value problem with multiple terms derivative containing both fractional and integer, we can use the **FODEMatrixDiscrete** algorithm to solve the equation.

All we have to do is passing the parameters and orders of the fractional ordinary differential equation to the API ```solve``` as two arrays.

!!! warning "The parameters and orders array must have the same length"
    When we are using ```FODEMatrixDiscrete``` to solve the problem, please note we should keep the parameters array and orders array must have the same length.

## Detailed Usage

Let's see if we have an equation like:

```math
2y''(t)+4D^{1.5}y(t)=1
```

To solve this equation, you can use the code:

```julia
rightside = 1
solve([2, 4], [2, 1.5], rightside, 30, 0.01)
```

Bingo! the result would represent the numerical solution of this equation!!!!

## Example

We have an initla problem:

```math
y'''(t)+\frac{1}{16} {^C_0D^{2.5}_ty(t)}+\frac{4}{5}y''(t)+\frac{3}{2}y'(t)+\frac{1}{25}{^C_0D^{0.5}_ty(t)}+\frac{6}{5}y(t)=\frac{172}{125}\cos(\frac{4t}{5})
```

```math
y(0)=0,\ y'(0)=0,\ y''(0)=0
```

Model this problem and solve the equation:

```julia
using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ A\\ complicated\\ example \$"

T=30
h=0.05
tspan = collect(0.05:h:T)

rightfun(x)=172/125*cos(4/5*x)
result = solve([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 1], rightfun, h, T, FODEMatrixDiscrete())

plot(tspan, result, title=s, legend=:bottomright)
```


By solving the equation and plotting the result, we can see its solution here:

![Solution](./assets/complicated_example.png)


!!! info "Example in GitHub"
    This example is an official example in the source code, please see the [example folder](https://github.com/SciFracX/FractionalDiffEq.jl/blob/master/examples/complicated_example.jl)