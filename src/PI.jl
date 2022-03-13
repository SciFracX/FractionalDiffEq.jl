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
struct PIEx <: FractionalDiffEqAlgorithm end

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
struct PIIm <: FractionalDiffEqAlgorithm end

"""
# Usage

    solve(prob::SingleTermFODEProblem, h, PITrap())

### References

```tex
@inproceedings{Garrappa2018NumericalSO,
  title={Numerical Solution of Fractional Differential Equations: A Survey and a Software Tutorial},
  author={Roberto Garrappa},
  year={2018}
}
```
"""
struct PITrap <: FractionalDiffEqAlgorithm end

function solve(FODE::SingleTermFODEProblem, h, ::PIEx)
    @unpack f, α, u0, T = FODE
    N::Int64 = round(T/h)
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

function solve(FODE::SingleTermFODEProblem, h, ::PIIm)
    @unpack f, α, u0, T = FODE
    N::Int64 = round(T/h)
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

function solve(FODE::SingleTermFODEProblem, h, ::PITrap)
    @unpack f, α, u0, T = FODE
    N::Int64 = round(T/h)
    y = zeros(N)
    y[1] = u0

    for n in range(2, N, step=1)
        middle=0
        for j=1:n
            middle += acoefficients(n-j, α)*f(j*h, y[j])
        end
        y[n] = u0 + h^α*(tacoefficients(n, α)+middle)
    end
    return y
end

function acoefficients(n, α)
    if n == 0
        return 1/gamma(α+2)
    else
        return ((n-1)^(α+1)-2*n^(α+1)+(n+1)^(α+1))/gamma(α+2)
    end
end

function tacoefficients(n, α)
    return ((n-1)^(α+1)-n^α*(n-α-1))/gamma(α+2)
end

function bcoefficients(n, α)
    return ((n+1)^α-n^α)/gamma(α+1)
end