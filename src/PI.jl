"""
# Usage

    solve(f, α, u0, T, h, PIEX())

### References

```tex
@inproceedings{Garrappa2018NumericalSO,
  title={Numerical Solution of Fractional Differential Equations: A Survey and a Software Tutorial},
  author={Roberto Garrappa},
  year={2018}
}
```
"""
struct PIEX <: FractionalDiffEqAlgorithm end

"""
# Usage

    solve(f, α, u0, T, h, PIEX())

### References

```tex
@inproceedings{Garrappa2018NumericalSO,
  title={Numerical Solution of Fractional Differential Equations: A Survey and a Software Tutorial},
  author={Roberto Garrappa},
  year={2018}
}
```
"""
struct PIIM <: FractionalDiffEqAlgorithm end

function solve(f, α, u0, T, h, ::PIEX)
    N::Int64 = T/h
    y = zeros(N)

    y[1]=u0
    for j in range(2, N, step=1)
        middle=0
        for i=0:j-1
            middle += bcoefficients(j-i, α)*f(i*h, y[i+1])
        end
        y[j] = u0 + middle*h^α
    end
    return y
end

function solve(f, α, u0, T, h, ::PIIM)
    N::Int64 = T/h
    y = zeros(N)

    y[1]=u0
    for j in range(2, N, step=1)
        middle=0
        for i=1:j
            middle += bcoefficients(j-i, α)*f(i*h, y[i])
        end
        y[j] = u0 + middle*h^α
    end
    return y
end

function bcoefficients(n, α)
    return ((n+1)^α-n^α)/gamma(α+1)
end


fun(t, y)=1-y
sol=solve(fun, 0.5, 0, 5, 0.5, PIEX())