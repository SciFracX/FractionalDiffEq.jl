using FractionalDiffEq, Plots
h = 0.05; t = collect(0:h:5); u0 = 1
fun(t)=0
prob = DODEProblem([1, 0.1], [x->6*x*(1-x), 0], (0, 1), fun, u0, t)
sol = solve(prob, h, DOMatrixDiscrete())
plot(sol)