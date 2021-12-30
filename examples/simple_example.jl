using FractionalDiffEq, Plots, LaTeXStrings

s="\$D^{0.5}y(x)=1-y,\\ y(0)=0\$"

fun(x, y) = 1-y
prob = SingleTermFODEProblem(fun, 0.5, 0.01)
result=solve(prob, 0, 5, PECE())
tspan=collect(0:0.01:5)

plot(tspan, result, title=s, linewidth=2, legend=:bottomright)