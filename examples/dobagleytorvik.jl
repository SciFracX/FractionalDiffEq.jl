using FractionalDiffEq, Plots, LaTeXStrings

h=0.075
t=collect(h:h:30)
f(t)=1
prob=DODEProblem([1, 1, 1], [2, x->6*x*(1-x), 0], [0.01, 1], t, h, f)
result=solve(prob, DOMatrixDiscrete())
plot(t, result)