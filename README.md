# FractionalDiffEq.jl

  <a href="https://github.com/SciFracX/FractionalDiffEq.jl/actions?query=workflow%3ACI">
    <img alt="building" src="https://github.com/SciFracX/FractionalDiffEq.jl/workflows/CI/badge.svg">
  </a>
  <a href="https://codecov.io/gh/SciFracX/FractionalDiffEq.jl">
    <img alt="codecov" src="https://codecov.io/gh/SciFracX/FractionalDiffEq.jl/branch/master/graph/badge.svg">
  </a>
  <a href="https://scifracx.github.io/FractionalDiffEq.jl/dev/">
    <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="license">
  </a>
  <a href="https://zenodo.org/badge/latestdoi/420992306">
  	<img src="https://zenodo.org/badge/420992306.svg" alt="DOI">
  </a>

  [![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

FractionalDiffEq.jl provides FDE solvers to [DifferentialEquations.jl](https://diffeq.sciml.ai/dev/) ecosystem, including FODE(Fractional Ordinary Differential Equations), FDDE(Fractional Delay Differential Equations) and many more. There are many performant solvers available, capable of solving many kinds of fractional differential equations.

# Installation

If you have already installed Julia, you can install FractionalDiffEq.jl in REPL using the Julia package manager:

```julia
pkg> add FractionalDiffEq
```

# Quick start

### Fractional ordinary differential equations

Let's see if we have an fractional order initial value problem:

$$ D^{1.8}y(x)=1-y $$


$$ y(0)=0 $$

So we can use FractionalDiffEq.jl to solve the problem:

```julia
using FractionalDiffEq, Plots
fun(u, p, t) = 1-u
u0 = [0, 0]; tspan = (0, 20);
prob = FODEProblem(fun, 1.8, u0, tspan)
sol = solve(prob, PECE(), dt=0.001)
plot(sol)
```

And if you plot the result, you can see the result of the above IVP:

![Example](/docs/src/assets/example.png)

### A sophisticated example

Let's see if we have an initial value problem:

$$ y'''(t)+\frac{1}{16}{^C_0D^{2.5}_t}y(t)+\frac{4}{5}y''(t)+\frac{3}{2}y'(t)+\frac{1}{25}{^C_0D^{0.5}_t}y(t)+\frac{6}{5}y(t)=\frac{172}{125}\cos(\frac{4t}{5}) $$

$$ y(0)=0,\ y'(0)=0,\ y''(0)=0 $$

```julia
using FractionalDiffEq, Plots
tspan = (0, 30)
rightfun(x, y) = 172/125*cos(4/5*x)
prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 0], rightfun, [0, 0, 0, 0, 0, 0], tspan)
sol = solve(prob, PIEX(), dt=0.01)
plot(sol, legend=:bottomright)
```

Or use the [example file](https://github.com/SciFracX/FractionalDiffEq.jl/blob/master/examples/complicated_example.jl) to plot the numerical approximation, we can see the FDE solver in FractionalDiffEq.jl is amazingly powerful:

![Example](docs/src/assets/complicated_example.png)

### System of Fractional Differential Equations:

FractionalDiffEq.jl is a powerful tool to solve a system of fractional differential equations, if you are familiar with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), it would be just like out of the box.

Let's see if we have a Chua chaos system:

$$ \begin{cases}D^{\alpha_1}x=10.725[y-1.7802x-[0.1927(|x+1|-|x-1|)]\\
D^{\alpha_2}y=x-y+z\\
D^{\alpha_3}z=-10.593y-0.268z\end{cases} $$

By using the ```NonLinearAlg``` algorithms to solve this problem:

```julia
using FractionalDiffEq, Plots
function chua!(du, x, p, t)
    a, b, c, m0, m1 = p
    du[1] = a*(x[2]-x[1]-(m1*x[1]+0.5*(m0-m1)*(abs(x[1]+1)-abs(x[1]-1))))
    du[2] = x[1]-x[2]+x[3]
    du[3] = -b*x[2]-c*x[3]
end
α = [0.93, 0.99, 0.92];
x0 = [0.2; -0.1; 0.1];
tspan = (0, 100);
p = [10.725, 10.593, 0.268, -1.1726, -0.7872]
prob = FODEProblem(chua!, α, x0, tspan, p)
sol = solve(prob, BDF(), dt=0.01)
plot(sol, vars=(1, 2), title="Chua System", legend=:bottomright)
```

And plot the result:

![Chua](docs/src/assets/chua.png)

## Fractional Delay Differential Equations

There are also many powerful solvers for solving fractional delay differential equations.

$$ D^\alpha_ty(t)=3.5y(t)(1-\frac{y(t-0.74)}{19}) $$

$$ y(0)=19.00001 $$


With history function:

$$ y(t)=19,\ t<0 $$

```julia
using FractionalDiffEq, Plots
ϕ(x) = x == 0 ? (return 19.00001) : (return 19.0)
f(t, y, ϕ) = 3.5*y*(1-ϕ/19)
α = 0.97; τ = 0.8; tspan = (0,56)
fddeprob = FDDEProblem(f, ϕ, α, constant_lags=[τ], tspan)
sol = solve(fddeprob, DelayPECE(), dt=0.05)
plot(sol, xlabel="y(t)", ylabel="y(t-τ)")
```

![Delayed](docs/src/assets/fdde_example.png)

# Available Solvers

For more performant solvers, please refer to the [FractionalDiffEq.jl Solvers](https://scifracx.org/FractionalDiffEq.jl/dev/algorithms/) page.
# Road map

* More performant algorithms
* Better docs
* More interesting ideas~

# Contributing

If you are interested in Fractional Differential Equations and Julia, welcome to raise an issue or file a Pull Request!!