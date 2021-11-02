using FractionalDiffEq
using Plots, LaTeXStrings
using MittagLeffler

#Analytical solution
target = []

#MittagLeffler.jl doesn't support array argument
for i in 0:0.01:20
    push!(target, i^1.8*mittleff(1.8,2.8,-i^1.8))
end

s="\$D^{0.5}y(x)=1-y,\\ y(0)=0\$"

#Numerical solution
fun(x, y) = 1-y
prob = FDEProblem(fun, 1.8, 0, 20, 0.01)
result=solve(prob, PECE())
tspan=collect(0:0.01:20)

gr()

plot(tspan, result, title=s, linewidth=5, label="Numerical", legend=:bottomright)
plot!(tspan, target, lw=3, ls=:dash, label="Analytical")