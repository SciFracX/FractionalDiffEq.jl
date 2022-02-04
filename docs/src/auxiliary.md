# Auxiliary functions in FractionalDiffEq.jl

There are some build-in auxiliary functions in FractionalDiffEq.jl.

## Generalized exponential function

### Mittag Leffler function

!!! info
    The mittag leffler function is adapted from [MittagLeffler.jl](https://github.com/jlapeyre/MittagLeffler.jl) implemented by [John Lapeyre](https://github.com/jlapeyre). We built in the mittag leffler function and added a few more functionalities.

The two-parametric [mittag leffler function](https://en.wikipedia.org/wiki/Mittag-Leffler_function) is defined by Gösta Magnus Mittag-Leffler bu a power series as:

```math
E_{\alpha,\ \beta}(z)=\sum_{k=0}^{\infty}\frac{z^k}{\Gamma(\alpha k+\beta)}
```

And single-parametric version:

```math
E_\alpha(z)=E_{\alpha,\ 1}(z)
```

And also the three-parametric Mittag Leffler function:

```math
E_{\alpha,\ \beta}^\gamma(z)=\sum_{k=0}^{\infty}\frac{(\gamma)_k}{k!\Gamma(\alpha k+\beta)}
```

Here ``(\gamma)_k=\frac{\Gamma(\gamma+k)}{\Gamma(\gamma)}`` is the [Pochhammer symbol](https://en.wikipedia.org/wiki/Falling_and_rising_factorials).

> Funny anecdote: The Mittag Leffler function can be seen as the Queen Function of Fractional Calculus.[^1]

In FractionalDiffEq.jl, you can compute the mittag leffler function by calling:

```julia-repl
julia> mittleff(α, β, z)
julia> mittleff(α, z)
```

Different order single parameter plot(``0<\alpha<1``):

![MittLeff](./assets/mittlefffun.png)

And also ``1<\alpha<2``:

![MittagLeffler](./assets/mittlefffunhigh.png)

For more interesting topics and applications about Mittag Leffler function, we recommend you to read [Mittag-Leffler Functions, Related Topics and Applications](https://link.springer.com/book/10.1007/978-3-662-61550-8)


### The generalized ``\alpha`` exponential function

The generalized ``\alpha`` exponential function ``e_\alpha^{\lambda z}`` is defined as:

```math
e_\alpha^{\lambda z}=z^{\alpha-1}E_{\alpha,\ \alpha}(\lambda z^\alpha)=\sum_{k=0}^{\infty}\frac{\lambda^k z^{\alpha(k-1)}}{\Gamma(\alpha(k+1))}\quad (z,\ \lambda\in\mathbb{C},\ \mathfrak{Re}(\alpha)>0)
```


[^1]: Page 2 of [Mittag-Leffler Functions, Related Topics and Applications](https://link.springer.com/book/10.1007/978-3-662-61550-8)