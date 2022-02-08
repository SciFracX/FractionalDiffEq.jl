"""
    SingleTermDODEProblem(ω, t, h, B, rightfun)

Define a single term distributed order differential equation problem.
"""
struct DODEProblem <: FDEProblem
    parameters::AbstractArray
    orders::AbstractArray
    interval
    tspan
    h::Float64
    rightfun::Function
end

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


isfunction(x) = isa(x, Function) ? true : false

function solve(prob::DODEProblem, ::DOMatrixDiscrete)
    parameters, orders, interval, tspan, h, rightfun = prob.parameters, prob.orders, prob.interval, prob.tspan, prob.h, prob.rightfun
    N = length(tspan)  
    DOid = findall(isfunction, orders)
    ϕ = orders[DOid][1]
    ω = parameters[DOid][1]
    modifiedorders = deleteat!(orders, DOid)
    modifiedparameters = deleteat!(parameters, DOid)

    highestorder = Int64(findmax(ceil.(modifiedorders))[1])
    rows=collect(1:highestorder)

    equation = zeros(N, N)

    for (i, j) in zip(modifiedparameters, modifiedorders)
        equation += i*D(N, j, h)
    end

    equation += ω.*DOB(ϕ, interval, 0.01, N, h)

    F = rightfun.(t)
    M = eliminator(N, rows)*equation*eliminator(N, 1)'

    F = eliminator(N, rows)*F

    Y = M\F

    Y0 = [0; Y]
    return Y0.+1
end
#=
h = 0.01; t = collect(h:h:5);
fun(t)=cos(t)
prob = DODEProblem([1, 1, 2, 3], [x->6*x*(1-x), 0, 1, 0.5], [0, 1], t, h, fun)

result = solve(prob, DOMatrixDiscrete())
using Plots
plot(t, result)
=#