using FractionalDiffEq
using Plots
using LaTeXStrings

fun(x, y) = 1-y
prob = FDEProblem(fun, 0.5, 0, 5, 0.001)
result=solve(prob)
tspan=collect(0:0.001:5)

#Want to use a more elegant plotting backend, especially for the LaTeX rendering, failed anyway (T_T)
#pgfplotsx()
gr()

plot(tspan, result, title="D^α y(x)=1-y, y(0)=0", linewidth=2, legend=:bottomright)