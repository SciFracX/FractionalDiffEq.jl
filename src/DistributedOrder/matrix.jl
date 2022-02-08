"""
    SingleTermDODEProblem(ω, t, h, B, rightfun)

Define a single term distributed order differential equation problem.
"""
struct SingleTermDODEProblem <: FDEProblem
    parameters::AbstractArray
    orders::AbstractArray
    ω::Function
    interval
    tspan
    h::Float64
    B
    rightfun::Function
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

isfunction(x) = isa(x, Function) ? true : false

function testsolve(M, t, h, B, rightfun)
    N = length(t)
    F = rightfun.(t)

    M = eliminator(N, 1)*M*eliminator(N, 1)'

    F = eliminator(N, 1)*F

    Y = M\F

    Y0 = [0; Y]

    return Y0.+1
end
#=
h = 0.01; t = collect(h:h:5);
#prob = SingleTermDODEProblem(x->6*x(1-x), [0, 1], t, h, 1, fun)
fun(t)=sin(t)
equation = DOB(x->6*x*(1-x), [0, 1], 0.01, 500, h) + D(500, 0, 0.01);
result=testsolve(equation, t, h, 1, fun);
using Plots
plot(t, result)
=#