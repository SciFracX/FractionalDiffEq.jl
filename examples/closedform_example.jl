using FractionalDiffEq#, Plots

rightfun(x)=sin(x^2)
prob = MultiTermsFODEProblem([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], rightfun, [30 90 3], [1 0.3 0], [0, 0], (0, 5))
sol = solve(prob, 0.5, ClosedForm())
plot(sol)