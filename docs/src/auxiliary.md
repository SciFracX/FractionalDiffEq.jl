# Auxiliary functions in FractionalDiffEq.jl

There are some build-in auxiliary functions in FractionalDiffEq.jl.

## Mittag Leffler function

> For more intreresting topics and applications about Mittag Leffler function, we recommend you to read [Mittag-Leffler Functions, Related Topics and Applications](https://link.springer.com/book/10.1007/978-3-662-43930-2)

Being called "The Queen Function of the Fractional Calculus, Mittag Leffler play an important role in fractional order computing and modeling. Here, we classify the Mittag Leffler function into three types:

### Classical Mittag Leffler function(Single-Parametric version)

The [mittag leffler function](https://en.wikipedia.org/wiki/Mittag-Leffler_function) is defined by Gösta Magnus Mittag-Leffler by a power series as:

```math
E_\alpha(z)=\sum_{k=0}^{\infty}\frac{z^k}{\Gamma(\alpha k+1)},\ (\alpha\in\mathbb{C})
```

### Two-Parametric Mittag Leffler function

The two-parametric Mittag Leffler function as the generalization of the classical Mittag Leffler function :

```math
E_{\alpha,\ \beta}(z)=\sum_{k=0}^{\infty}\frac{z^k}{\Gamma(\alpha k+\beta)},\ (Re(\alpha)>0,\ \beta\in\mathbb{C})
```

### Three-Parametric Mittag Leffler function

And also the three-parametric Mittag Leffler function:

```math
E_{\alpha,\ \beta}^\gamma(z)=\sum_{k=0}^{\infty}\frac{(\gamma)_k}{k!\Gamma(\alpha k+\beta)},\ (Re(\alpha)>0,\ Re(\beta)>0,\ \gamma\in\mathbb{C})
```

Here ``(\gamma)_k=\frac{\Gamma(\gamma+k)}{\Gamma(\gamma)}`` is the [Pochhammer symbol](https://en.wikipedia.org/wiki/Falling_and_rising_factorials).

In FractionalDiffEq.jl, you can compute the mittag leffler function by calling:

```julia
julia> mittleff(α, z)
julia> mittleff(α, β, z)
julia> mittleff(α, β, γ, z)
```

Different order classical Mittag Leffler function plot(``0<\alpha<1``):

![MittLeff](./assets/mittlefffun.png)

And also ``1<\alpha<2``:

![MittagLeffler](./assets/mittlefffunhigh.png)

### The generalized ``\alpha`` exponential function

The generalized ``\alpha`` exponential function ``e_\alpha^{\lambda z}`` is defined as:

```math
e_\alpha^{\lambda z}=z^{\alpha-1}E_{\alpha,\ \alpha}(\lambda z^\alpha)=\sum_{k=0}^{\infty}\frac{\lambda^k z^{\alpha(k-1)}}{\Gamma(\alpha(k+1))}\quad (z,\ \lambda\in\mathbb{C},\ \mathfrak{Re}(\alpha)>0)
```