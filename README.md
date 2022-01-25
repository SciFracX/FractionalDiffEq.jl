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

> See our talk on JuliaCN 2021 Winter Conf: [Slide](https://julia-cn-conf2021.vercel.app/1), [YouTube](https://www.youtube.com/watch?v=oVvrW7EgEwg), [BiliBili](https://www.bilibili.com/video/BV1vY411W7Dw?p=18)

# Installation

If you have already installed Julia, you can install FractionalDiffEq.jl in REPL using Julia package manager:

```julia
pkg> add FractionalDiffEq
```

Or if you want to experience the latest version of FractionalDiffEq.jl:

```julia
pkg> add FractionalDiffEq#master
```

# Quick start

### An easy example

Let's see if we have an initial value problem:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?D^{0.5}y(x)=1-y" title="D^{0.5}y(x)=1-y" />

</p>

<p align="center">

<img src="https://latex.codecogs.com/svg.image?y(0)=0" title="y(0)=0" />

</p>


So we can use FractionalDiffEq.jl to solve the problem:

```julia
fun(x, y) = 1-y
prob = SingleTermFODEProblem(fun, 0.5, 0, 5)
result = solve(prob, 0.001, PECE())
tspan = collect(0:0.001:5)
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
using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ A\\ complicated\\ example \$"

T = 30
h = 0.05
tspan = collect(0.05:h:T)

rightfun(x) = 172/125*cos(4/5*x)

prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 1], rightfun)

result = solve(prob, h, T, FODEMatrixDiscrete())
plot(tspan, result, title=s, legend=:bottomright)
```

Or use the [example file](https://github.com/SciFracX/FractionalDiffEq.jl/blob/master/examples/complicated_example.jl) to plot the numerical approximation, we can see the FDE solver in FractionalDiffEq.jl is amazingly powerful:

![Example](docs/src/assets/complicated_example.png)

### System of Fractional Differential Equations:

FractionalDiffEq.jl is a powerful tool to solve system of fractional differential equations:

Let's see if we have a Chua chaos system:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?\begin{cases}D^{\alpha_1}x=10.725[y-1.7802x-[0.1927(|x&plus;1|-|x-1|)]\\D^{\alpha_2}y=x-y&plus;z\\D^{\alpha_3}z=-10.593y-0.268z\end{cases}" title="\begin{cases}D^{\alpha_1}x=10.725[y-1.7802x-[0.1927(|x+1|-|x-1|)]\\D^{\alpha_2}y=x-y+z\\D^{\alpha_3}z=-10.593y-0.268z\end{cases}" />

</p>

By using the ```NonLinearAlg``` algorithms to solve this problem:

```julia
using FractionalDiffEq
using Plots

function chua(t, x, k)
    a=10.725
    b=10.593
    c=0.268
    m0=-1.1726
    m1=-0.7872

    if k==1
        f=m1*x[1]+0.5*(m0-m1)*(abs(x[1]+1)-abs(x[1]-1))
        y=a*(x[2]-x[1]-f)
        return y
    elseif k==2
        y=x[1]-x[2]+x[3]
        return y
    elseif k==3
        y=-b*x[2]-c*x[3]
        return y
    end
end

alpha = [0.93, 0.99, 0.92];
x0 = [0.2; -0.1; 0.1];
h = 0.01;
tn = 200;
result = solve(chua, alpha, x0, h, tn, NonLinearAlg())

gr()
plot(result[:, 1], result[:, 2], title="Chua System", legend=:bottomright)
```

And plot the result:

![Chua](docs/src/assets/chua.png)

## Fractional Partial Differential Equations

Fractional provide powerful algorithms to solve fractional partial differential equations, let's see a diffusion example here:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?_{0}^{C}\!D_{t}^{\alpha}y-&space;\frac{\partial^\beta&space;y}{\partial&space;|x|^\beta}&space;=&space;f(x,t)" title="_{0}^{C}\!D_{t}^{\alpha}y- \frac{\partial^\beta y}{\partial |x|^\beta} = f(x,t)" />

</p>

With initial and boundry conditions:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?y(0,t)&space;=&space;0,&space;\quad&space;y(1,t)&space;=&space;0&space;\qquad&space;&space;\quad&space;y(x,0)&space;=&space;0" title="y(0,t) = 0, \quad y(1,t) = 0 \qquad  \quad y(x,0) = 0" />

</p>

By using the FPDE solvers in FractionalDiffEq.jl and plot the numerical approximation:

![diffusion](docs/src/assets/diffusion.png)

### ODE Example

FractionalDiffEq.jl is also able to solve ordinary differential equations~ Let's see an example here:

<p align="center">

<img src="https://latex.codecogs.com/svg.image?y''(x)&plus;y'(x)=\sin(x)" title="y''(x)+y'(x)=\sin(x)" />

</p>

<p align="center">

<img src="https://latex.codecogs.com/svg.image?y(0)=0" title="y(0)=0" />

</p>


```julia
using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ODE\\ Example\$"

T = 30
h = 0.05
tspan = collect(0.05:h:T)

f(x) = 1/2*(-exp(-x)-sin(x)-cos(x)+2)
target =f .(tspan)

rightfun(x) = sin(x)
prob = MultiTermsFODEProblem([1, 1], [2, 1], rightfun)
result = solve(prob, h, T, FODEMatrixDiscrete())

plot(tspan, result, title=s, legend=:bottomright, label="ODE Numerical Solution!")

plot!(tspan, target, lw=3,ls=:dash,label="ODE Analytical Solution!")
```

![ODE Example](docs/src/assets/ode_example.png)

## Road map

* More performant algorithms
* Better docs
* More interesting ideas~

## Contributing

If you are interested in Fractional Differential Equations and Julia, welcome to raise an issue or file a Pull Request!!