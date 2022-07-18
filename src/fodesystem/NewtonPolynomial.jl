"""
    solve(prob::FODESystem, h, NewtonPolynomial())

Solve FODE system using Newton Polynomials methods.

!!! tip
    Used for the Caputo Fabrizio fractional differential operators.

```tex
https://doi.org/10.1016/c2020-0-02711-8
```
"""
struct NewtonPolynomial <: AbstractFDEAlgorithm end

function solve(prob::FODESystem, h, ::NewtonPolynomial)
    @unpack f, α, u0, tspan, p = prob
    t0 = tspan[1]; tfinal = tspan[2]
    α = α[1]
    t = collect(Float64, t0:h:tfinal)
    M = 1-α+α/gamma(α)
    N::Int = ceil(Int, (tfinal-t0)/h)
    l = length(u0)
    result = zeros(l, N+1)

    result[:, 1]=u0
    temp = zeros(l)
    f(temp, result[:, 1], p, t[1])
    result[:, 2] = result[:, 1]+h.*temp
    temptemp = zeros(l)
    f(temptemp, result[:, 2], p, t[2])
    result[:, 3] = result[:, 2]+(h/2).*(3 .*temptemp-temp)

    temp1 = zeros(Float64, l)
    temp2 = zeros(Float64, l)
    temp3 = zeros(Float64, l)

    for n=3:N
        f(temp1, result[:, n], p, t[n])
        f(temp2, result[:, n-1], p, t[n-1])
        f(temp3, result[:, n-2], p, t[n-2])
        result[:, n+1] = result[:, n] + (1-α)/M*(temp1-temp2)+α.*M.*h.*(23/12*temp1-4/3*temp2+5/12*temp3)
    end
    return FODESystemSolution(t, result)
end