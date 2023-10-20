using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ODE\\ Example\$"

T = 30
h = 0.05
tspan = collect(0.05:h:T)

f(x)=1/2*(-exp(-x)-sin(x)-cos(x)+2)
target=f.(tspan)

rightfun(x) = sin(x)
prob = MultiTermsFODEProblem([1, 1], [2, 1], rightfun, [0, 0], T)
sol = solve(prob, h, FODEMatrixDiscrete())

plot(sol, title=s, legend=:bottomright, label="ODE Numerical Solution!")

plot!(tspan, target, lw=3,ls=:dash,label="ODE Analytical Solution!")