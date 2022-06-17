"""
    solve(prob, h, Euler())

The basic forward Euler method for fractional ordinary differential equations.

```tex
@inproceedings{Li2015NumericalMF,
  title={Numerical Methods for Fractional Calculus},
  author={Changpin Li and Fanhai Zeng},
  year={2015}
}
```
"""
struct Euler <: FractionalDiffEqAlgorithm end

function solve(prob::SingleTermFODEProblem, h, ::Euler)
    @unpack f, α, u0, tspan = prob
    t0 = tspan[1]; tfinal = tspan[2]
    t = collect(t0:h:tfinal)
    N::Int = ceil(Int, (tfinal-t0)/h)
    y = zeros(Float64, N+1)
    y[1] = sum(u0)
    for n=1:N
        temp = zero(Float64)
        for j=1:n
            temp += ((n-j+1)^α-(n-j)^α)*f(t[j], y[j])
        end
        y[n+1] = y[1] + h^α/gamma(α+1)*temp
    end
    return FODESolution(t, y)
end