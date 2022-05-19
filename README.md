# FractionalDiffEq.jl

<p align="center">
<img width="250px" src="https://raw.githubusercontent.com/SciFracX/FractionalDiffEq.jl/master/docs/src/assets/logo.svg"/>
</p>


<p align="center">
  <a href="https://github.com/SciFracX/FractionalDiffEq.jl/actions?query=workflow%3ACI">
    <img alt="building" src="https://github.com/SciFracX/FractionalDiffEq.jl/workflows/CI/badge.svg">
  </a>
  <a href="https://codecov.io/gh/SciFracX/FractionalDiffEq.jl">
    <img alt="codecov" src="https://codecov.io/gh/SciFracX/FractionalDiffEq.jl/branch/master/graph/badge.svg">
  </a>
  <a href="https://scifracx.github.io/FractionalDiffEq.jl/dev/">
    <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="license">
  </a>
  <a href="https://github.com/SciFracX/FractionalDiffEq.jl/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/SciFracX/FractionalDiffEq.jl?style=flat-square" alt="license">
  </a>
  <a href="https://zenodo.org/badge/latestdoi/420992306">
  	<img src="https://zenodo.org/badge/420992306.svg" alt="DOI">
  </a>
</p>

<p align="center">
  <a href="https://github.com/SciFracX/FractionalDiffEq.jl/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/SciFracX/FractionalDiffEq.jl?style=flat-square">
  </a>
  <a href="#">
    <img alt="GitHub stars" src="https://img.shields.io/github/stars/SciFracX/FractionalDiffEq.jl?style=flat-square">
  </a>
  <a href="https://github.com/SciFracX/FractionalDiffEq.jl/network">
    <img alt="GitHub forks" src="https://img.shields.io/github/forks/SciFracX/FractionalDiffEq.jl?style=flat-square">
  </a>
</p>

FractionalDiffEq.jl provides FDE solvers to [DifferentialEquations.jl](https://diffeq.sciml.ai/dev/) ecosystem. There are many performant solvers available, capable of solving many kinds of fractional differential equations.

# Installation

If you have already installed Julia, you can install FractionalDiffEq.jl in REPL using Julia package manager:

```julia
pkg> add FractionalDiffEq
```

# Quick start

### Fractional ordinary differential equations

Let's see if we have an initial value problem:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?D^{0.5}y(x)=1-y" title="D^{0.5}y(x)=1-y" />

</p>

<p align="center">

<img src="https://latex.codecogs.com/svg.image?y(0)=0" title="y(0)=0" />

</p>


So we can use FractionalDiffEq.jl to solve the problem:

```julia
using FractionalDiffEq, Plots
fun(x, y) = 1-y
u0 = 0; T = 5; h = 0.001;
prob = SingleTermFODEProblem(fun, 0.5, u0, T)
sol = solve(prob, h, PECE())
plot(sol)
```

And if you plot the result, you can see the result of the fractional differential equation:

![Example](/docs/src/assets/simple_example.png)

### A sophisticated example

Let's see if the initial value problem like:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?y'''(t)&plus;\frac{1}{16}{^C_0D^{2.5}_t}y(t)&plus;\frac{4}{5}y''(t)&plus;\frac{3}{2}y'(t)&plus;\frac{1}{25}{^C_0D^{0.5}_t}y(t)&plus;\frac{6}{5}y(t)=\frac{172}{125}\cos(\frac{4t}{5})" title="y'''(t)+\frac{1}{16}{^C_0D^{2.5}_t}y(t)+\frac{4}{5}y''(t)+\frac{3}{2}y'(t)+\frac{1}{25}{^C_0D^{0.5}_t}y(t)+\frac{6}{5}y(t)=\frac{172}{125}\cos(\frac{4t}{5})" />

</p>

<p align="center">

<img src="https://latex.codecogs.com/svg.image?y(0)=0,\&space;y'(0)=0,\&space;y''(0)=0" title="y(0)=0,\ y'(0)=0,\ y''(0)=0" />

</p>

```julia
using FractionalDiffEq, Plots
h=0.01; tspan = (0, 30)
rightfun(x, y) = 172/125*cos(4/5*x)
prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], (0, T))
sol = solve(prob, h, PIEX())
plot(sol, legend=:bottomright)
```

Or use the [example file](https://github.com/SciFracX/FractionalDiffEq.jl/blob/master/examples/complicated_example.jl) to plot the numerical approximation, we can see the FDE solver in FractionalDiffEq.jl is amazingly powerful:

![Example](docs/src/assets/complicated_example.png)

### System of Fractional Differential Equations:

FractionalDiffEq.jl is a powerful tool to solve system of fractional differential equations, if you are familiar with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), it would be just like out of the box.

Let's see if we have a Chua chaos system:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?\begin{cases}D^{\alpha_1}x=10.725[y-1.7802x-[0.1927(|x&plus;1|-|x-1|)]\\D^{\alpha_2}y=x-y&plus;z\\D^{\alpha_3}z=-10.593y-0.268z\end{cases}" title="\begin{cases}D^{\alpha_1}x=10.725[y-1.7802x-[0.1927(|x+1|-|x-1|)]\\D^{\alpha_2}y=x-y+z\\D^{\alpha_3}z=-10.593y-0.268z\end{cases}" />

</p>

By using the ```NonLinearAlg``` algorithms to solve this problem:

```julia
using FractionalDiffEq; Plots
function chua!(du, x, p, t)
    a = 10.725; b = 10.593; c = 0.268; m0 = -1.1726; m1 = -0.7872
    du[1] = a*(x[2]-x[1]-(m1*x[1]+0.5*(m0-m1)*(abs(x[1]+1)-abs(x[1]-1))))
    du[2] = x[1]-x[2]+x[3]
    du[3] = -b*x[2]-c*x[3]
end
α = [0.93, 0.99, 0.92];
x0 = [0.2; -0.1; 0.1];
h = 0.001; tn = 0.5;
prob = FODESystem(chua!, α, x0, tn)
result = solve(prob, h, NonLinearAlg())
plot(result[:, 1], result[:, 2], title="Chua System", legend=:bottomright)
```

And plot the result:

![Chua](docs/src/assets/chua.png)

## Fractional Partial Differential Equations

FractionalDiffEq.jl provides powerful algorithms to solve fractional partial differential equations, let's see a diffusion equation here:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?_{0}^{C}\!D_{t}^{\alpha}y-&space;\frac{\partial^\beta&space;y}{\partial&space;|x|^\beta}&space;=&space;f(x,t)" title="_{0}^{C}\!D_{t}^{\alpha}y- \frac{\partial^\beta y}{\partial |x|^\beta} = f(x,t)" />

</p>

With initial and boundry conditions:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?y(0,t)&space;=&space;0,&space;\quad&space;y(1,t)&space;=&space;0&space;\qquad&space;&space;\quad&space;y(x,0)&space;=&space;0" title="y(0,t) = 0, \quad y(1,t) = 0 \qquad  \quad y(x,0) = 0" />

</p>

Use the FPDE solvers in FractionalDiffEq.jl and plot the numerical approximation:

![diffusion](docs/src/assets/diffusion.png)


## Fractional Delay Differential Equations

There are also many powerful solvers for solving fractional delay differential equations.

<p align="center">

<img src="https://latex.codecogs.com/svg.image?D^\alpha_ty(t)=3.5y(t)(1-\frac{y(t-0.74)}{19}),\&space;y(0)=19.00001" title="https://latex.codecogs.com/svg.image?D^\alpha_ty(t)=3.5y(t)(1-\frac{y(t-0.74)}{19}),\ y(0)=19.00001" />

</p>

With history function:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?y(t)=19,\&space;t<0" title="https://latex.codecogs.com/svg.image?y(t)=19,\ t<0" />

</p>

```julia
using FractionalDiffEq, Plots
ϕ(x) = x == 0 ? (return 19.00001) : (return 19.0)
f(t, y, ϕ) = 3.5*y*(1-ϕ/19)
h = 0.05; α = 0.97; τ = 0.8; T = 56
fddeprob = FDDEProblem(f, ϕ, α, τ, T)
V, y = solve(fddeprob, h, DelayPECE())
plot(y, V, xlabel="y(t)", ylabel="y(t-τ)")
```

![Delayed](docs/src/assets/fdde_example.png)

## Fractional Integral Equations

To solve fractional integral equations with FractionalDiffEq.jl, we only need to follow the previous process:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?u(x)&plus;{_{-1}I_x^{1/2}}u(x)=1" title="https://latex.codecogs.com/svg.image?u(x)+{_{-1}I_x^{1/2}}u(x)=1" />

</p>

```julia
using FractionalDiffEq, Plots, SpecialFunctions
analytical(x)=exp(1+x)*erfc(sqrt(1+x))
tspan = LinRange(-1, 1, 100)
prob = FIEProblem([1, 1], [1, 0.5], 1, tspan)
sol = solve(prob, 20, SpectralUltraspherical())
solanalytical = analytical.(xx)
plot(sol, title="Second kind Abel integral equation", label="Numerical")
plot!(tspan, solanalytical, ls=:dash, label="Analytical")
```

![Second kind Abel IE](docs/src/assets/abelinteqexample.png)


# Available Solvers

For more performant solvers, please refer to the [FractionalDiffEq.jl Solvers](https://scifracx.org/FractionalDiffEq.jl/dev/algorithms/) page.
# Road map

* More performant algorithms
* Better docs
* More interesting ideas~

# Contributing

If you are interested in Fractional Differential Equations and Julia, welcome to raise an issue or file a Pull Request!!