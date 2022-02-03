# Distributed Order Differential Equations

We integrate ``{_0D_t^\alpha f(t)}`` with respect to the order, we can obtain distributed-order differential/integral equations. Since we normally adopted the notation as:

```math
{_0D_t^{\omega(\alpha)} f(t)} := \int_{\gamma_1}^{\gamma_2}\omega(\alpha){_0D_t^\alpha f(t)}d\alpha
```

We can write the general form of distributed order differential equations as:

```math
\int_0^m \mathscr{A}(r,\ D_*^r u(t))dr = f(t)
```

Similar with what we have learned about single-term and multi-term fractional differential equations in linear fractional differential equations, we can also write the single-term distributed order differential equations:

```math
D_*^ru(t)=f(t,\ u(t))
```

And multi-term distributed order differential equations

```math
\sum_{i=1}^k \gamma_i D_*^{r_i}u(t) = f(t,\ u(t))
```

## Distributed Order Relaxation



!!! tip
    Please see [Distributed-Order Dynamic Systems](https://link.springer.com/book/10.1007/978-1-4471-2852-6) for systematic introduction and knowledge.