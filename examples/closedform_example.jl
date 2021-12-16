using FractionalDiffEq
using Plots, LaTeXStrings

s="\$ Use\\ closed-form\\ solution\$"

t=collect(0:0.002:10);

u=ones(length(t));

#result=solve([1 8 26 73 90], [3.5 3.1 2.3 1.2 0.5], [30 90], [1 0.3], u, t, ClosedForm())
result = solve([1 1], [1.8 0], [])
plot(t, result)