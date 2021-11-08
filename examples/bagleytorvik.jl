using FractionalDiffEq
using Plots, LaTeXStrings

s="\$Bagley\\ Torvik\\ Equation\$"

T=30
h=0.05
tspan = collect(0:h:T)
result = bagleytorvik(1, 1, 1, T, h)

plot(tspan, result, title=s, legend=:bottomright)
savefig("./bagleytorvik.png")