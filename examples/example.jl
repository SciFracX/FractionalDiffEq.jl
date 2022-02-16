using FractionalDiffEq
using Plots, LaTeXStrings

# Analytical solution
analytical(x) = x.^1.8 .*mittleff(1.8, 2.8, -x.^1.8)



s="\$D^{1.8}y(x)=1-y(x),\\ y(0)=0\$"

# Numerical solution
fun(x, y) = 1-y
prob = SingleTermFODEProblem(fun, 1.8, 0, 20)
result = solve(prob, 0.01, PECE())
tspan = collect(0:0.01:20)
target = analytical(tspan)

gr()

plot(tspan, result, title=s, linewidth=5, label="Numerical", legend=:bottomright)
plot!(tspan, target, lw=3, ls=:dash, label="Analytical")