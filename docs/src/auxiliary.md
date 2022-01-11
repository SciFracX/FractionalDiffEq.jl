# Auxiliary functions in FractionalDiffEq.jl

There are some build-in auxiliary functions in FractionalDiffEq.jl.

## Mittag Leffler function

> The mittag leffler function is adapted from [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl) implemented by [John Lapeyre](https://github.com/jlapeyre). We build in the mittag leffler function and add a few more functionalities.

The [mittag leffler function](https://en.wikipedia.org/wiki/Mittag-Leffler_function) is defined as:

```math
E_{\alpha,\ \beta}(z)=\sum_{k=0}^{\infty}\frac{z^k}{\Gamma(\alpha k+\beta)}
```

And single parameter version:

```math
E_\alpha(z)=E_{\alpha,\ 1}(z)
```

In FractionalDiffEq.jl, you can compute the mittag leffler function by calling:

```julia-repl
julia> mittleff(α, β, z)
julia> mittleff(α, z)
```

Different order single parameter plot(``0<\alpha<1``):

![MittLeff](./assets/mittlefffun.png)

And also ``1<\alpha<2``:

![MittagLeffler](./assets/mittlefffunhigh.png)