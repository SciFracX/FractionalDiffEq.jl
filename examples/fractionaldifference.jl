using FractionalDiffEq, Plots

fun(x) = 0.5*x+1
α=0.5;x0=1;
T=1; h=0.1
t, y=solve(fun, α, x0, T, h, PECEDifference())
plot(t, y, xlims=(0, 10), legend=:bottomright)