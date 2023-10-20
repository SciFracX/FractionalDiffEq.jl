using FractionalDiffEq, Plots

fun(x) = 0.5*x+1
α=0.5;x0=1;
T=1; h=0.1
prob = FractionalDifferenceProblem(fun, α, x0)
sol=solve(prob, T, h, PECEDifference())
plot(sol, seriestype=:scatter, legend=:bottomright)