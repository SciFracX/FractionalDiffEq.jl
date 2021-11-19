# Multi term FDE

By specifying different order in the equation, we can handle multi-term differential equations now!

Let's see if we have a initial value problem with multiple terms derivative containing both fractional and integer, we can use the **FODEMatrixDiscrete** algorithm to solve the equation.

All we have to do is use the general derivative representing function ```D(size, order, step)``` to represent different derivative, for example, ```D(30, 2, 0.01)``` represent the second order derivative $y''(t)$ term and ```D(30, 2.5, 0.01)``` represent the 2.5 order derivative $D^{2.5}y(t)$ term.

!!! warning "Keep the parameter unanimous"
When we are use ```D``` to represent different order deriavtives, please note we should keep the first parameter and third parameter unanimous, which represent the size of the discrete matrix and step size.

## Detailed Usage

Let's see if we have a equation like:

```math
2y''(t)+4D^{1.5}y(t)=1
```

To solve this equation, you can use the code:

```julia
equation = 2*D(30, 2, 0.01) + 4*D(30, 1.5, 0.01)
rightside = 1
solve(equation, rightside, 2, 30, 0.01)
```

Bingo! the result would represent the numerical solution of this equation!!!!

## Example

We have an initla problem:

![LaTeX](./assets/complicated_example_latex.png)

Model this problem and solve the equation:

```julia
using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ An\\ complicated\\ example \$"

T=30
h=0.05
tspan = collect(0.05:h:T)

equation = D(600, 3, h)+1/16*D(600, 2.5, h)+4/5*D(600, 2, h)+3/2*D(600, 1, h)+1/25*D(600, 0.5, h)+6/5*D(600, 1, h);
rightfun(x)=172/125*cos(4/5*x)
result=solve(equation, rightfun, 3, h, T)

result=result.-1 .-4/5 .*tspan .+16/25 .*tspan.^2

plot(tspan, result, title=s, legend=:bottomright)
savefig("./complicated_example.png")
```

!!! info "Example in GitHub"
This example is an official example in the source code, please see the [example folder](https://github.com/SciFracX/FractionalDiffEq.jl/blob/master/examples/complicated_example.jl)

By solving the equation and plotting the result, we can see its solution here:

![Solution](./assets/complicated_example.png)