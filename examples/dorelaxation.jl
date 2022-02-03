using FractionalDiffEq, Plots, LaTeXStrings

s = "\$Distributed\\ Order\\ Relaxation\$"

h = 0.01
t = collect(0:h:5)
result = solve(x->6*x*(1-x), t, h, 0.1, DOMatrixDiscrete())

using Plots
plot(t, result, title=s)