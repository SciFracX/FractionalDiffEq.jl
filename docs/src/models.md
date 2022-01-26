# Detailed Models

## Diffusion equation.

### Diffusion equation

**Diffusion equation** is a classical partial equation widely used in Physics to describe **density fluctuations** in a material undergoing diffusion.

```math
\frac{\partial u}{\partial t}=\nabla\cdot(\kappa \nabla u)
```

Here, $\kappa$ is the [diffusion coefficient](https://en.wikipedia.org/wiki/Mass_diffusivity) which is spatially-dependent. When $\kappa$ is a constant, **Diffusion equation** evolves to **[Heat equation](https://en.wikipedia.org/wiki/Heat_equation)**:

```math
\frac{\partial u}{\partial t}=\kappa\nabla^2u
```

The Diffusion Equation belongs to [Partial Differential Equations](https://en.wikipedia.org/wiki/Partial_differential_equation), and there are many amazing organizations and packages are dedicate to solving partial differential equations using high performance and advancing algorithms:

See [Survey of PDE Packages](https://github.com/JuliaPDE/SurveyofPDEPackages)

### Differential Equations with spatial fractional derivative

When we are modeling our problems, **spatial fractional derivative** maybe more suitable for our model, which changes the diffusion equation to **Diffusion equation with spatial fractional derivative**

```math
\frac{\partial u}{\partial t}-\frac{\partial^\beta u}{\partial |x|^\beta}=f(t, u)
```


### Diffusion Equations with time fractional derivative

When we are modeling our problems, **time fractional derivative** maybe more suitable for our model, which changes the diffusion equation to **Diffusion equation with time fractional derivative**

```math
^C_0D^\alpha_tu=\frac{\partial^2u}{\partial x^2}
```

### General fractional diffusion equation

Well, time fractional derivative and spatial fractional derivative are both need to take into consideration:

```math
^C_0D^\alpha_tu-\frac{\partial^\beta u}{\partial|x|^\beta}=f(t, u)
```

## Bagley Torvik Equation

The Bagley Torvik can be used to describe the moving of a damped object.

![bagleytorvik](./assets/damped.png)

```math
Ay''(t)+BD^{\frac{3}{2}}_ty(t)+Cy(t)=f(t)
```

In FractionalDiffEq.jl, we can specify the parameters and solve the equation:

```julia
using FractionalDiffEq
using Plots, LaTeXStrings

s="\$Bagley\\ Torvik\\ Equation\$"

T=30
h=0.05
tspan = collect(0:h:T)
result = bagleytorvik(1, 1, 1, 1, T, h)

plot(tspan, result, title=s, legend=:bottomright)
```



## Relaxation Oscillation Equation

![Relaxo](./assets/Relaxo.png)

In relaxation oscillation process, the fractional calculus can be introduced to better describe the model.

https://en.wikipedia.org/wiki/Relaxation_oscillator