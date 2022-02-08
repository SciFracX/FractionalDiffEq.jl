using FractionalDiffEq, Plots, LaTeXStrings

s = "\$Distributed\\ Order\\ Relaxation\$"

h = 0.01; t = collect(h:h:5);
fun(t)=-0.1
prob = DODEProblem([1, 0.1], [x->6*x*(1-x), 0], [0, 1], t, h, fun)
result = solve(prob, DOMatrixDiscrete())
plot(t, result, title=s)