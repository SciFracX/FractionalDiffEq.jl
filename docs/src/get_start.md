# Get Start

## What is Fractional Differential Equations?

While the Ordinary Differential Equations and Partial Differential Equations are widely used in enormous areas and play inportant roles in their theoretical analysis, someone may asks, ODE and PDE are enough for nowadays modeling, has FDE any usage in our life?

Well, fractional differential equation can be seen as the generalization of ODE and PDE. In our daily life, models usually are better described in fractional differential equations.

A special applying case for fractional differential equations is viscoelasty, which research the property of a subject with both [Viscosity](https://en.wikipedia.org/wiki/Viscosity) and [Elasticity](https://en.wikipedia.org/wiki/Elasticity_(physics)). 

And also the CRONE controller and $PI^\lambda D^\mu$ controller deploy fractional derivative to better describe system. 

## A simple example —— Relaxation Oscillation Equation

Let's see a simple model involving fractional differential equations: [Relaxation Oscillation Process]()

```math
D^{1.8}y(t)+y(t)=1,\ (t>0)
```

```math
y^{(k)}(0)=0
```

We can solve the Relaxation Oscillation Equation using FractionalDiffEq.jl:

```julia
using FractionalDiffEq
using Plots, LaTeXStrings

s="\$Bagley\\ Torvik\\ Equation\$"

T=30
h=0.05
tspan = collect(0:h:T)
result = bagleytorvik(1, 1, 1, T, h)

plot(tspan, result, title=s, legend=:bottomright)
```

By ploting the numerical result, we can get the approximation result:

![Relaxation Oscillation](./assets/example.png)

## FDE with specific initial value

While the former examples we only use the zero initial value problem, here, we can look at some problems with non-zero examples.

## Some Algorithms explanation:

As a matter of fact, to solve a fractional differential equations is to solve a volterra integral equation:

