"""
# Usage

    solve(prob, DOMatrixDiscrete())

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


isFunction(x) = isa(x, Function) ? true : false

function solve(prob::DODEProblem, h, ::DOMatrixDiscrete)
    @unpack parameters, orders, interval, rightfun, tspan = prob
    N = length(tspan)
    # find the index of the distributed order
    DOid = findall(isFunction, orders)

    ϕ = orders[DOid][1]
    ω = parameters[DOid][1]
    modifiedorders = deleteat!(orders, DOid)
    modifiedparameters = deleteat!(parameters, DOid)

    highestorder = Int64(findmax(ceil.(modifiedorders))[1])
    rows=collect(1:highestorder)

    equation = zeros(N, N)

    # Construct systems using matrices
    for (i, j) in zip(modifiedparameters, modifiedorders)
        equation += i*D(N, j, h)
    end

    # Don't forget the distributed order term!
    equation += ω.*DOB(ϕ, interval, 0.01, N, h)

    F = eliminator(N, rows)*rightfun.(tspan)
    M = eliminator(N, rows)*equation*eliminator(N, 1)'

    Y = M\F

    Y0 = [0; Y]
    return Y0
end

#=
h = 0.01; t = collect(h:h:5);
fun(t)=-0.1
prob = DODEProblem([1, 0.1], [x->6*x*(1-x), 0], [0, 1], t, h, fun)

result = solve(prob, DOMatrixDiscrete())
using Plots
plot(t, result)
=#