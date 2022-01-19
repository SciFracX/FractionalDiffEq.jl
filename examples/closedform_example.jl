using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ Use\\ closed-form\\ solution\$"

t = collect(0:0.002:10);

u = ones(length(t));

prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], [30 90 3], [1 0.3 0])


#prob = MultiTermsFODEProblem([1 5 9 7 2], [1.2 0.9 0.6 0.3 0], 1, 0)

result=solve(prob, u, t, ClosedForm())

plot(t, result)