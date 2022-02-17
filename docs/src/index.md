```@meta
CurrentModule = FractionalDiffEq
```

# FractionalDiffEq.jl

Hello there👋!

FractionalDiffEq.jl is a Julia package aiming at solving Fractional Differential Equations using high performance numerical methods.

While the Ordinary Differential Equations and Partial Differential Equations are widely used in enormous areas and play important roles in their theoretical analysis, someone may asks, ODE and PDE are enough for nowadays modeling, has FDE any usage in our life?

Well, fractional differential equation can be seen as the generalization of ODE and PDE. In our daily life, models usually are better described in fractional differential equations.

```@contents
Pages = ["index.md"]
```

!!! tip "Star Us"
	If you think **FractionalDiffEq.jl** is useful, please [star us in GitHub](httpd://github.com/SciFracX/FractionalDiffEq.jl), it means a lot to us!

## Installation

If you have already installed Julia, you can install FractionalDiffEq.jl in REPL using Julia package manager:

```julia
pkg> add FractionalDiffEq
```

Or if you want to experience the latest version of FractionalDiffEq.jl:

```julia
pkg> add FractionalDiffEq#master
```

## Features

* While most fractional differential equations solvers are implemented using Matlab, **FractionalDiffEq.jl** is totally driven by [Julia](https://julialang.org/) and licensed with [MIT License](https://en.wikipedia.org/wiki/MIT_License), ensuring its everlasting development and open source.

* Benefit from the advancing features of JuliaLang, FractionalDiffEq.jl has impressive performance and high speed, help the model more efficient and robust.

* Capable of solving both linear and nonlinear fractional differential equations. Including fractional ordinary differential equations, fractional partial differential equations, fractional delayed differential equations, distributed order differential equations and system of fractional differential equations.

* Detailed models supporting, such as Bagley Torvik equations, Relaxation Oscillation equations and Diffusion equations many more.

## Roadmap

* More algorithms are planned to support.

* Improve benchmark.

* Connect with SciML ecosystem.

* More interesting ideas~


## Contributing

Just by using FractionalDiffEq.jl you're already contributing!

The development of FractionalDiffEq.jl is on [GitHub](https://github.com/SciFracX/FractionalDiffEq.jl), please [report bugs](https://github.com/SciFracX/FractionalDiffEq.jl/issues) or [send pull requests](https://github.com/SciFracX/FractionalDiffEq.jl/pulls) to help FractionalDiffEq.jl grow.

## Acknowledge

**FractionalDiffEq.jl** is built upon the hard work of many scientific researchers, I sincerely appreciate what they have done to help the development of science and technology.

!!! info "WIP"
	FractionalDiffEq.jl is under heavy construction, some APIs or docs might change a lot.