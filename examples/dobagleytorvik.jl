using FractionalDiffEq, Plots
h=0.075; t=collect(0:h:30)
f(t) = 0 ≤ t ≤ 1 ? (return 8) : (return 0)
prob=DODEProblem([1, 1, 1], [2, x->6*x*(1-x), 0], (0, 1), f, [0; 0], t)
sol=solve(prob, h, DOMatrixDiscrete())
plot(sol)