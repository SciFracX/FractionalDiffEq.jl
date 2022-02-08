"""
    SingleTermDODEProblem(ω, t, h, B, rightfun)

Define a single term distributed order differential equation problem.
"""
struct SingleTermDODEProblem <: FDEProblem
    ω
    t
    h
    B
    rightfun
end

"""
# Usage

    solve()

Use distributed order strip matrix algorithm to solve distriubted order problem.

### References

```tex
@inproceedings{Jiao2012DistributedOrderDS,
  title={Distributed-Order Dynamic Systems - Stability, Simulation, Applications and Perspectives},
  author={Zhuang Jiao and Yang Quan Chen and Igor Podlubny},
  booktitle={Springer Briefs in Electrical and Computer Engineering},
  year={2012}
}
```
"""
struct DOMatrixDiscrete <: FractionalDiffEqAlgorithm end

function solve(ω, t, h, fun, ::DOMatrixDiscrete)
    N = length(t)
    M = zeros(N, N)

    M = DOB(ω, [0, 1], 0.01, N, h)
    F = fun.(t)

    M = eliminator(N, 1)*M*eliminator(N, 1)'
    F = eliminator(N, 1)*F

    Y = M\F

    Y0 = [0; Y]

    return Y0.+1
end

#=
h = 0.01; t = collect(0:h:5)
fun(t)=sin(t)
result=testsolve(x->6*x*(1-x), t, h, fun)
using Plots
plot(t, result)
=#