# Fractional Order Partial Differential Equations

Several models of physical and biological processes are better described using fractional PDEs than the corresponding integer-order PDEs. They serve as a generalization of the integer-order PDEs and give some degree of freedom in varying the rate of change of these physical and biological processes.

FractionalDiffEq.jl has support for solving fractional partial differential equations. Let's see an example here:

Time fractional differential equation:

```math
_{0}^{C}\!D_{t}^{\alpha}y - \frac{\partial^\beta y}{\partial |x|^\beta} = f(x,t)
```

With initial and boundry conditions:

```math
y(0,t) = 0, \quad y(1,t) = 0 \qquad  \quad y(x,0) = 0
```

We can use the ```FPDEMatrixDiscrete``` algorithm to solve this problem:

```julia
using FractionalDiffEq
using Plots

tmp = solve(0.7, 1.8, 1, 21, 148, FPDEMatrixDiscrete())



YS = reshape(tmp, 19, 147)
YS = reverse(YS, dims=2)
U = YS


rows, columns = size(U)

U = [zeros(1, columns); U; zeros(1, columns)]
U=[zeros(1, 21)' U]

XX, YY = meshgrid(0.05^2/6 .*(0:147), 0:0.05:1)

plotly()

plot(XX, YY, U, st=:surface)
```

![Diffusion](./assets/diffusion.png)

## Time fractional diffusion equation

When we are modeling our problems, **time fractional derivative** maybe more suitable for our model, which changes the diffusion equation to **Diffusion equation with time fractional derivative**

```math
^C_0D^\alpha_tu=\kappa\frac{\partial^2u}{\partial x^2}
```
Here, ``\kappa`` is the [diffusion coefficient](https://en.wikipedia.org/wiki/Mass_diffusivity)

```julia
using FractionalDiffEq

K = 1
fdorder = 1.9
dx = pi/20
dt = 0.01
n = 2
xStart = 0
xEnd = pi

using Plots
plotlyjs()
U=solve(fdorder, dx, dt, xStart, xEnd, n, K, FiniteDiffEx())
plot(x, t, U, st=:surface)
```

![DiffusionEx](./assets/finitediffex.png)