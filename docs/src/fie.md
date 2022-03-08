# Fractional Integral Equations

!!! warning "WIP"
    Supporting for fractional integral equatoins is still under heavy development.

The second-kind Abel integral equation:

```math
u(x)+{_{-1}Q_x^{1/2}}u(x)=1
```

While the solution is

```math
u(x)=e^{1+x}erfc(\sqrt{1+x})
```

```julia
using FractionalDiffEq, Plots

e(x)=1+0*x
f(x)=0*x
xx = LinRange(-1, 1, 100)
result=solve(f, e, 20, SpectralUltraspherical())
sol=myeval(result, xx, 2)
plot(xx, sol)
```

![Second kind Abel IE](./assets/abelinteqexample.png)