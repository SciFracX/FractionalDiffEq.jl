using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ A\\ complicated\\ example \$"

T = 30
h = 0.05
tspan = collect(0.05:h:T)

rightfun(x) = 172/125*cos(4/5*x)

prob = MultiTermsFODEProblem([1, 1/16, 4/5, 3/2, 1/25, 6/5], [3, 2.5, 2, 1, 0.5, 1], rightfun)

result = solve(prob, h, T, FODEMatrixDiscrete())
plot(tspan, result, title=s, legend=:bottomright)