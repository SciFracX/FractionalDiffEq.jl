"""
    solve(prob::FODESystem, h, AtanganaSedaCF())

Atagana Seda method for Caputo-Fabrizio fractional order differential equations.

!!! tip
    Used for the Caputo Fabrizio fractional differential operators.

```tex
https://doi.org/10.1016/c2020-0-02711-8
```
"""
struct AtanganaSedaCF <: AbstractFDEAlgorithm end
#FIXME: Tests
function solve(prob::FODESystem, h, ::AtanganaSedaCF)
    @unpack f, α, u0, tspan, p = prob
    t0 = tspan[1]; tfinal = tspan[2]
    t=collect(Float64, t0:h:tfinal)
    α=α[1]
    M=1-α+α/gamma(α)
    N=ceil(Int, (tfinal-t0)/h)
    l=length(u0)
    result = zeros(l, N+1)
    temp1 = zeros(l); temp2 = zeros(l)
    f(temp1, u0, p, t[1])
    result[:, 2] = u0+temp1
    f(temp2, result[:, 2], p, t[2])
    result[:, 3] = result[:, 2] + (h/2)*(3*temp2-temp1)
    tempn, tempn1, tempn2 = zeros(l), zeros(l), zeros(l)
    for n=3:N
        f(tempn, result[:, n], p, t[n])
        f(tempn1, result[:, n-1], p, t[n-1])

        f(tempn2, result[:, n-2], p, t[n-2])
        result[:, n+1] = result[:, n] + (1-α)/M*(tempn-tempn1)+α*M*h*(23/12*tempn-4/3*tempn1+5/12*tempn2)
    end
    return result
end