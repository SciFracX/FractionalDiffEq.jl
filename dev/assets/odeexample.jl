using FractionalDiffEq
using Plots, LaTeXStrings

T = 30
h=0.05
tspan = collect(0.05:h:T)

f(x)=1/2*(-exp(-x)-sin(x)-cos(x)+2)
target=f.(tspan)

eq = D(600, 2, h)+D(600, 1, h)
rightfun(x) = sin(x)
result = solve(eq, rightfun, 2, h, T, FODEMatrixDiscrete())

plot(tspan, result, legend=:bottomright, label="ODE numerical solution!")

plot!(tspan, target, lw=3,ls=:dash,label="ODE analytical Solution!")