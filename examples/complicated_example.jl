using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ A\\ complicated\\ example \$"

T=30
h=0.05
tspan = collect(0.05:h:T)

equation = D(600, 3, h)+1/16*D(600, 2.5, h)+4/5*D(600, 2, h)+3/2*D(600, 1, h)+1/25*D(600, 0.5, h)+6/5*D(600, 1, h);
rightfun(x)=172/125*cos(4/5*x)
result=solve(equation, rightfun, 3, h, T)

plot(tspan, result, title=s, legend=:bottomright)