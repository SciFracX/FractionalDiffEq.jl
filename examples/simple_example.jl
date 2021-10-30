using FractionalDiffEq, Plots

fun(x, y) = 1-y
prob=FDEProblem(fun, 1.8, 0, 10, 0.01)
result=solve(prob)
tspan=collect(0:0.01:10)

plot(tspan, result, title="D^0.5 y(x)=1-y, y(0)=0", linewidth=2, legend=:bottomright)