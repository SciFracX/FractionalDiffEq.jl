"""
    solve(FDDE::FDDEProblem, T, h, DelayPECE())

Using the delayed predictor-corrector method to solve the delayed fractional differential equation problem.

```tex
@article{Wang2013ANM,
  title={A Numerical Method for Delayed Fractional-Order Differential Equations},
  author={Zhen Wang},
  journal={J. Appl. Math.},
  year={2013},
  volume={2013},
  pages={256071:1-256071:7}
}
```
"""
struct DelayPECE <: FractionalDiffEqAlgorithm end

function solve(FDDE::FDDEProblem, T, h, ::DelayPECE)
    f, ϕ, α, τ = FDDE.f, FDDE.ϕ, FDDE.α, FDDE.τ
    t = collect(0:h:T)
    maxn = size(t, 1)
    yp = collect(0:h:T+h)
    y = copy(t)
    y[1] = ϕ(0)

    for n in 1:maxn-1
        yp[n+1] = 0
        for j = 1:n
            yp[n+1] = yp[n+1]+b(j-1, n-1, α, h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp))
        end
        yp[n+1] = yp[n+1]/gamma(α)+ϕ(0)

        y[n+1] = 0

        for j=1:n
            y[n+1] = y[n+1]+a(j-1, n-1, α, h)*f(t[j], y[j], v(ϕ, j, τ, h, y, yp))
        end

        y[n+1] = y[n+1]/gamma(α)+h^α*f(t[n+1], yp[n+1], v(ϕ, n+1, τ, h, y, yp))/gamma(α+2)+ϕ(0)
    end

    V = copy(t)
    for n = 1:maxn-1
        V[n] = v(ϕ, n, τ, h, y, yp)
    end
    return V, y
end

function a(j, n, α, h)
    if j == n+1
        result = 1
    elseif j == 0
        result = n^(α+1)-(n-α)*(n+1)^α
    else
        result = (n-j+2)^(α+1) + (n-j)^(α+1) - 2*(n-j+1)^(α+1)
    end
    return result*h^α / (α * (α + 1))
end

function b(j, n, α, h)
    return h^α/α*((n-j+1)^α - (n-j)^α)
end

function v(ϕ, n, τ, h, y, yp)
    if τ >= (n-1)*h
        return ϕ((n-1)*h-τ)
    else
        m = Int64(floor(τ/h))
        delta = m-τ/h

        if m>1
            return delta*y[n-m+2] + (1-delta)*y[n-m+1]
        else
            return delta*yp[n+1] + (1-delta)*y[n]
        end
    end
end