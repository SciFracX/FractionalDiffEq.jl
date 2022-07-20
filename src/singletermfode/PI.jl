"""
# Usage

    solve(prob::SingleTermFODEProblem, h, PIEX())

### References

```tex
@inproceedings{Garrappa2018NumericalSO,
  title={Numerical Solution of Fractional Differential Equations: A Survey and a Software Tutorial},
  author={Roberto Garrappa},
  year={2018}
}
```
"""
struct PIEX <: SingleTermFODEAlgorithm end

function solve(FODE::SingleTermFODEProblem, h::Float64, ::PIEX)
    @unpack f, α, u0, tspan = FODE
    t0 = tspan[1]; T = tspan[2]
    N::Int64 = round(Int, (T-t0)/h)
    y = zeros(Float64, N+1)

    y[1]=u0
    for j in range(2, N+1, step=1)
        middle=0
        @turbo for i=0:j-1
            middle += bcoefficients(j-i, α)*f(t0+i*h, y[i+1])
        end
        middle = middle/gamma(α+1)
        y[j] = u0 + middle*h^α
    end
    t = collect(Float64, t0:h:T)
    return FODESolution(t, y)
end
#=
function acoefficients(n, α)
    if n == 0
        return 1
    else
        return ((n-1)^(α+1)-2*n^(α+1)+(n+1)^(α+1))
    end
end
=#

bcoefficients(n, α) = ((n+1)^α-n^α)