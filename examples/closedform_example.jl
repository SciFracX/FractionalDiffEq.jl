using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ Use\\ closed-form\\ solution\$"

t = collect(0:0.002:10);
rightfun(x)=sin(x^2)
prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90 3], [1 0.3 0], t)
sol = solve(prob, ClosedForm())
plot(sol)