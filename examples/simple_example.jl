using FractionalDiffEq, Plots

fun(x, y) = 1-y
prob=FDEProblem(fun, 0.5, 0, 5, 0.001)
result=solve(prob)
tspan=collect(0:0.001:5)

plot(tspan, result, title="D^0.5 y(x)=1-y, y(0)=0", linewidth=2, legend=:bottomright)