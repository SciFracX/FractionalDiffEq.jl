# FractionalDiffEq

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
  <a href="https://www.erikqqy.xyz/FractionalDiffEq.jl/dev/">
    <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="license">
  </a>
  <a href="https://github.com/SciFracX/FractionalDiffEq.jl/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/SciFracX/FractionalDiffEq.jl?style=flat-square" alt="license">
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


# Installation

If you have already installed Julia, you can install FractionalDiffEq.jl in REPL using Julia package manager:

```julia
Pkg> add FractionalDiffEq
```

Or if you want to experience the latest version of FractionalDiffEq.jl:

```julia
Pkg> add FractionalDiffEq#master
```

# Quick start

### An easy example

Let's see if we have an initial value problem:

<p align="center">
<img width="250px" src="https://raw.githubusercontent.com/SciFracX/FractionalDiffEq.jl/master/docs/src/assets/simple_example_latex.png"/>
</p>

So we can use FractionalDiffEq.jl to solve the problem:

```julia
fun(x, y) = 1-y
prob = FODEProblem(fun, 0.5, 0, 5, 0.001)
result=solve(prob, PECE())
tspan=collect(0:0.001:5)
```

And if you plot the result, you can see the result of the fractional differential equation:

![Example](/docs/src/assets/simple_example.png)

### A sophisticated example

Let's see if the initial value problem like:

<p align="center">
<img width="600px" src="https://raw.githubusercontent.com/SciFracX/FractionalDiffEq.jl/master/docs/src/assets/complicated_example_latex.png"/>
</p>

Use the [example file](https://github.com/SciFracX/FractionalDiffEq.jl/blob/master/examples/complicated_example.jl) to plot the numerical approximation, we can see the FDE solver in FractionalDiffEq.jl is amazingly powerful:

![Example](docs/src/assets/complicated_example.png)

## Road map

* Trapezoid PI rule
* Fractional Linear Multi-step Method

## Contributing

If you are interested in Fractional Differential Equations and Julia, welcome to raise an issue or file a Pull Request!!